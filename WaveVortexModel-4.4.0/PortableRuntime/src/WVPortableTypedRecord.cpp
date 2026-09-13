#include "WaveVortexRuntime/WVPortableTypedRecord.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <type_traits>

namespace wavevortex::runtime {
namespace {

constexpr std::size_t encodedLengthBytes = sizeof(std::uint64_t);

std::size_t saturatingAdd(std::size_t left, std::size_t right) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left)
    return std::numeric_limits<std::size_t>::max();
  return left + right;
}

std::size_t saturatingMultiply(std::size_t left, std::size_t right) noexcept {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    return std::numeric_limits<std::size_t>::max();
  return left * right;
}

std::size_t expectedCount(const std::vector<std::size_t> &dimensions,
                          bool &overflow) noexcept {
  overflow = false;
  if (dimensions.empty())
    return 1;
  std::size_t result = 1;
  for (const auto dimension : dimensions) {
    if (dimension != 0 &&
        result > std::numeric_limits<std::size_t>::max() / dimension) {
      overflow = true;
      return 0;
    }
    result *= dimension;
  }
  return result;
}

} // namespace

namespace {

template <typename Value>
void appendScalar(std::vector<std::uint8_t> &bytes, Value value) {
  const auto *source = reinterpret_cast<const std::uint8_t *>(&value);
  bytes.insert(bytes.end(), source, source + sizeof(Value));
}

void appendString(std::vector<std::uint8_t> &bytes, const std::string &value) {
  appendScalar<std::uint64_t>(bytes, value.size());
  bytes.insert(bytes.end(), value.begin(), value.end());
}

template <typename Value>
bool readScalar(const std::vector<std::uint8_t> &bytes, std::size_t &offset,
                Value &value) {
  if (offset > bytes.size() || sizeof(Value) > bytes.size() - offset)
    return false;
  std::memcpy(&value, bytes.data() + offset, sizeof(Value));
  offset += sizeof(Value);
  return true;
}

bool readString(const std::vector<std::uint8_t> &bytes, std::size_t &offset,
                std::string &value) {
  std::uint64_t count = 0;
  if (!readScalar(bytes, offset, count) || count > bytes.size() - offset)
    return false;
  value.assign(reinterpret_cast<const char *>(bytes.data() + offset),
               static_cast<std::size_t>(count));
  offset += static_cast<std::size_t>(count);
  return true;
}

} // namespace

WVPortableValueType WVPortableNamedValue::valueType() const noexcept {
  switch (storage.index()) {
  case 0:
    return WVPortableValueType::boolean;
  case 1:
    return WVPortableValueType::integer;
  case 2:
    return WVPortableValueType::real;
  case 3:
    return WVPortableValueType::text;
  }
  return WVPortableValueType::real;
}

std::size_t WVPortableNamedValue::valueCount() const noexcept {
  return std::visit([](const auto &values) { return values.size(); }, storage);
}

std::size_t WVPortableNamedValue::encodedBytes() const noexcept {
  std::size_t result = encodedLengthBytes + name.size() + sizeof(std::uint8_t) +
                       encodedLengthBytes;
  result = saturatingAdd(
      result, saturatingMultiply(dimensions.size(), sizeof(std::uint64_t)));
  result = saturatingAdd(result, encodedLengthBytes);
  std::visit(
      [&](const auto &values) {
        using Value = typename std::decay_t<decltype(values)>::value_type;
        if constexpr (std::is_same_v<Value, std::string>) {
          for (const auto &value : values) {
            result = saturatingAdd(result, encodedLengthBytes);
            result = saturatingAdd(result, value.size());
          }
        } else {
          result = saturatingAdd(
              result, saturatingMultiply(values.size(), sizeof(Value)));
        }
      },
      storage);
  return result;
}

std::size_t WVPortableNamedValue::persistentBytes() const noexcept {
  std::size_t result = sizeof(*this) + name.capacity() +
                       dimensions.capacity() * sizeof(std::size_t);
  std::visit(
      [&](const auto &values) {
        using Value = typename std::decay_t<decltype(values)>::value_type;
        result = saturatingAdd(result, values.capacity() * sizeof(Value));
        if constexpr (std::is_same_v<Value, std::string>)
          for (const auto &value : values)
            result = saturatingAdd(result, value.capacity());
      },
      storage);
  return result;
}

const WVPortableNamedValue *
WVPortableTypedRecord::value(const std::string &name) const noexcept {
  const auto found =
      std::find_if(values.begin(), values.end(), [&](const auto &candidate) {
        return candidate.name == name;
      });
  return found == values.end() ? nullptr : &*found;
}

std::size_t WVPortableTypedRecord::encodedBytes() const noexcept {
  std::size_t result = encodedLengthBytes + schemaIdentifier.size() +
                       sizeof(schemaVersion) + encodedLengthBytes;
  for (const auto &value : values)
    result = saturatingAdd(result, value.encodedBytes());
  return result;
}

std::size_t WVPortableTypedRecord::persistentBytes() const noexcept {
  std::size_t result = sizeof(*this) + schemaIdentifier.capacity() +
                       values.capacity() * sizeof(WVPortableNamedValue);
  for (const auto &value : values)
    result = saturatingAdd(result, value.persistentBytes() - sizeof(value));
  return result;
}

WVKernelStatus
validatePortableTypedRecord(const WVPortableTypedRecord &record,
                            WVPortableTypedRecordValidation validation) {
  if (record.schemaIdentifier.empty() || record.schemaVersion == 0)
    return {WVKernelStatusCode::invalidConfiguration,
            "A portable typed record requires a schema identifier and "
            "positive schema version."};
  std::set<std::string> names;
  for (const auto &value : record.values) {
    if (value.name.empty() || !names.insert(value.name).second)
      return {WVKernelStatusCode::invalidConfiguration,
              "Portable typed-record value names must be nonempty and "
              "unique."};
    bool overflow = false;
    const auto count = expectedCount(value.dimensions, overflow);
    if (overflow)
      return {WVKernelStatusCode::sizeOverflow,
              "Portable typed-record dimensions overflow size_t."};
    if (count != value.valueCount())
      return {WVKernelStatusCode::invalidShape,
              "A portable typed-record value does not match its declared "
              "dimensions."};
    if (value.valueType() == WVPortableValueType::text && !validation.allowText)
      return {WVKernelStatusCode::invalidConfiguration,
              "Text values are not permitted in this portable typed record."};
    if (value.valueType() == WVPortableValueType::boolean) {
      const auto &values = std::get<std::vector<std::uint8_t>>(value.storage);
      if (std::any_of(values.begin(), values.end(),
                      [](const auto item) { return item > 1; }))
        return {WVKernelStatusCode::invalidConfiguration,
                "Portable Boolean values must be encoded as zero or one."};
    }
    if (value.valueType() == WVPortableValueType::real &&
        validation.requireFiniteReals) {
      const auto &values = std::get<std::vector<double>>(value.storage);
      if (std::any_of(values.begin(), values.end(),
                      [](const auto item) { return !std::isfinite(item); }))
        return {WVKernelStatusCode::invalidConfiguration,
                "Portable real values must be finite in this record."};
    }
  }
  const auto bytes = record.encodedBytes();
  if (bytes == std::numeric_limits<std::size_t>::max())
    return {WVKernelStatusCode::sizeOverflow,
            "Portable typed-record encoded size overflowed size_t."};
  if (bytes > validation.maximumEncodedBytes)
    return {WVKernelStatusCode::sizeOverflow,
            "Portable typed record exceeds its encoded-size limit."};
  return WVKernelStatus::ok();
}

WVKernelStatus encodePortableTypedRecord(const WVPortableTypedRecord &record,
                                         std::vector<std::uint8_t> &bytes) {
  auto status = validatePortableTypedRecord(record);
  if (!status)
    return status;
  try {
    std::vector<std::uint8_t> candidate;
    candidate.reserve(record.encodedBytes());
    appendString(candidate, record.schemaIdentifier);
    appendScalar(candidate, record.schemaVersion);
    appendScalar<std::uint64_t>(candidate, record.values.size());
    for (const auto &value : record.values) {
      appendString(candidate, value.name);
      appendScalar(candidate, static_cast<std::uint8_t>(value.valueType()));
      appendScalar<std::uint64_t>(candidate, value.dimensions.size());
      for (const auto dimension : value.dimensions)
        appendScalar<std::uint64_t>(candidate, dimension);
      appendScalar<std::uint64_t>(candidate, value.valueCount());
      std::visit(
          [&](const auto &values) {
            using Element = typename std::decay_t<decltype(values)>::value_type;
            if constexpr (std::is_same_v<Element, std::string>) {
              for (const auto &item : values)
                appendString(candidate, item);
            } else if (!values.empty()) {
              const auto *source =
                  reinterpret_cast<const std::uint8_t *>(values.data());
              candidate.insert(candidate.end(), source,
                               source + values.size() * sizeof(Element));
            }
          },
          value.storage);
    }
    bytes = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Portable typed-record encoding allocation failed."};
  }
}

WVKernelStatus
decodePortableTypedRecord(const std::vector<std::uint8_t> &bytes,
                          WVPortableTypedRecord &record,
                          WVPortableTypedRecordValidation validation) {
  if (bytes.size() > validation.maximumEncodedBytes)
    return {WVKernelStatusCode::sizeOverflow,
            "Portable typed record exceeds its encoded-size limit."};
  try {
    WVPortableTypedRecord candidate;
    std::size_t offset = 0;
    std::uint64_t valueCount = 0;
    if (!readString(bytes, offset, candidate.schemaIdentifier) ||
        !readScalar(bytes, offset, candidate.schemaVersion) ||
        !readScalar(bytes, offset, valueCount) || valueCount > bytes.size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Portable typed-record encoding is truncated or malformed."};
    candidate.values.reserve(static_cast<std::size_t>(valueCount));
    for (std::uint64_t valueIndex = 0; valueIndex < valueCount; ++valueIndex) {
      WVPortableNamedValue value;
      std::uint8_t type = 0;
      std::uint64_t dimensionCount = 0, count = 0;
      if (!readString(bytes, offset, value.name) ||
          !readScalar(bytes, offset, type) || type > 3 ||
          !readScalar(bytes, offset, dimensionCount) ||
          dimensionCount > bytes.size())
        return {WVKernelStatusCode::invalidConfiguration,
                "Portable typed-record value header is malformed."};
      value.dimensions.reserve(static_cast<std::size_t>(dimensionCount));
      for (std::uint64_t dimension = 0; dimension < dimensionCount;
           ++dimension) {
        std::uint64_t extent = 0;
        if (!readScalar(bytes, offset, extent) ||
            extent > std::numeric_limits<std::size_t>::max())
          return {WVKernelStatusCode::sizeOverflow,
                  "Portable typed-record dimension is invalid."};
        value.dimensions.push_back(static_cast<std::size_t>(extent));
      }
      if (!readScalar(bytes, offset, count) || count > bytes.size())
        return {WVKernelStatusCode::invalidConfiguration,
                "Portable typed-record value count is malformed."};
      if (type == static_cast<std::uint8_t>(WVPortableValueType::text)) {
        std::vector<std::string> values;
        values.reserve(static_cast<std::size_t>(count));
        for (std::uint64_t item = 0; item < count; ++item) {
          std::string text;
          if (!readString(bytes, offset, text))
            return {WVKernelStatusCode::invalidConfiguration,
                    "Portable typed-record text value is truncated."};
          values.push_back(std::move(text));
        }
        value.storage = std::move(values);
      } else {
        const std::size_t widths[] = {sizeof(std::uint8_t),
                                      sizeof(std::int64_t), sizeof(double)};
        const auto width = widths[type];
        if (count > (bytes.size() - offset) / width)
          return {WVKernelStatusCode::invalidConfiguration,
                  "Portable typed-record numeric value is truncated."};
        if (type == 0) {
          std::vector<std::uint8_t> values(count);
          if (count != 0)
            std::memcpy(values.data(), bytes.data() + offset, count * width);
          value.storage = std::move(values);
        } else if (type == 1) {
          std::vector<std::int64_t> values(count);
          if (count != 0)
            std::memcpy(values.data(), bytes.data() + offset, count * width);
          value.storage = std::move(values);
        } else {
          std::vector<double> values(count);
          if (count != 0)
            std::memcpy(values.data(), bytes.data() + offset, count * width);
          value.storage = std::move(values);
        }
        offset += static_cast<std::size_t>(count) * width;
      }
      candidate.values.push_back(std::move(value));
    }
    if (offset != bytes.size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Portable typed-record encoding has trailing bytes."};
    auto status = validatePortableTypedRecord(candidate, validation);
    if (!status)
      return status;
    record = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Portable typed-record decoding allocation failed."};
  }
}

bool samePortableTypedRecordValue(const WVPortableTypedRecord &left,
                                  const WVPortableTypedRecord &right) noexcept {
  if (left.schemaIdentifier != right.schemaIdentifier ||
      left.schemaVersion != right.schemaVersion ||
      left.values.size() != right.values.size())
    return false;
  for (std::size_t index = 0; index < left.values.size(); ++index) {
    const auto &first = left.values[index];
    const auto &second = right.values[index];
    if (first.name != second.name || first.dimensions != second.dimensions ||
        first.storage != second.storage)
      return false;
  }
  return true;
}

bool isCanonicalEmptyPortableTypedRecord(
    const WVPortableTypedRecord &record) noexcept {
  return record.schemaIdentifier.empty() && record.schemaVersion == 0 &&
         record.values.empty();
}

} // namespace wavevortex::runtime
