#include "WaveVortexRuntime/WVOutputSchedule.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <utility>

namespace wavevortex::runtime {
namespace {

double tolerance(double first, double second) noexcept {
  return 8.0 * std::numeric_limits<double>::epsilon() *
         std::max({1.0, std::abs(first), std::abs(second)});
}

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

constexpr std::uint64_t fnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

void hashBytes(std::uint64_t &hash, const void *bytes,
               std::size_t count) noexcept {
  const auto *values = static_cast<const std::uint8_t *>(bytes);
  for (std::size_t index = 0; index < count; ++index) {
    hash ^= values[index];
    hash *= fnvPrime;
  }
}

std::size_t scalarBytes(WVOutputSchedulePayloadType type) noexcept {
  switch (type) {
  case WVOutputSchedulePayloadType::real64:
    return sizeof(double);
  case WVOutputSchedulePayloadType::integer64:
    return sizeof(std::int64_t);
  case WVOutputSchedulePayloadType::boolean8:
    return sizeof(std::uint8_t);
  }
  return 0;
}

std::size_t alignedOffset(std::size_t offset, std::size_t alignment) noexcept {
  const auto remainder = offset % alignment;
  return remainder == 0 ? offset : offset + alignment - remainder;
}

WVKernelStatus validatePayloadAccess(
    const WVOutputSchedulePayloadSchema &schema,
    const WVOutputSchedulePayload &payload, std::size_t slot,
    WVOutputSchedulePayloadType type, std::size_t count) noexcept {
  if (payload.schemaFingerprint() != schema.fingerprint() ||
      payload.byteCount() != schema.payloadBytes())
    return invalid("An occurrence payload does not match its resolved schema.");
  if (slot >= schema.slotCount())
    return invalid("An occurrence-payload slot is out of range.");
  const auto &resolved = schema.slots()[slot];
  if (resolved.scalarType != type || resolved.elementCount != count)
    return invalid("An occurrence-payload slot has the wrong type or extent.");
  return WVKernelStatus::ok();
}

const WVOutputSchedulePayloadSchema &makeEmptyPayloadSchema() {
  static const auto schema = [] {
    WVOutputSchedulePayloadSchema result;
    const auto status = WVOutputSchedulePayloadSchema::create(
        "empty-occurrence-payload-v1", 1, {}, result);
    (void)status;
    return result;
  }();
  return schema;
}

class WVEvenlySpacedSchedule final : public WVOutputSchedule {
public:
  explicit WVEvenlySpacedSchedule(const WVOutputScheduleRecord &record)
      : anchor_(record.initialTime), interval_(record.outputInterval),
        final_(record.finalTime) {}

  const char *typeIdentifier() const noexcept override {
    return WVEvenlySpacedOutputScheduleType;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  const WVOutputSchedulePayloadSchema &payloadSchema() const noexcept override {
    return emptyOutputSchedulePayloadSchema();
  }

  WVKernelStatus
  validateCursor(const WVOutputScheduleCursor &cursor) const override {
    if (cursor.committedOrdinal < WVNoCommittedOutputOrdinal)
      return invalid("An output-schedule ordinal cannot be less than -1.");
    if (!cursor.values.schemaIdentifier.empty() ||
        cursor.values.schemaVersion != 0 || !cursor.values.values.empty())
      return invalid("The evenly-spaced schedule cursor contains unexpected "
                     "provider state.");
    if (cursor.committedOrdinal >= 0) {
      const auto time = std::fma(static_cast<double>(cursor.committedOrdinal),
                                 interval_, anchor_);
      if (!std::isfinite(time) || time > final_ + tolerance(time, final_))
        return invalid("Committed output progress lies beyond its schedule.");
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus committedTime(const WVOutputScheduleCursor &cursor,
                               double &time, bool &available) const override {
    auto status = validateCursor(cursor);
    if (!status)
      return status;
    available = cursor.committedOrdinal >= 0;
    if (available)
      time = std::fma(static_cast<double>(cursor.committedOrdinal), interval_,
                      anchor_);
    return WVKernelStatus::ok();
  }

  WVKernelStatus peek(const WVOutputScheduleCursor &cursor, double lowerBound,
                      double upperBound, WVOutputScheduleOccurrence &occurrence,
                      bool &available) const override {
    available = false;
    auto status = validateCursor(cursor);
    if (!status)
      return status;
    if (!std::isfinite(lowerBound) || !std::isfinite(upperBound) ||
        upperBound < lowerBound)
      return invalid(
          "Output-schedule bounds must be finite and nondecreasing.");
    const auto lower = std::max(lowerBound, anchor_);
    const auto upper = std::min(upperBound, final_);
    if (upper + tolerance(lower, upper) < lower)
      return WVKernelStatus::ok();
    const long double ratio =
        (static_cast<long double>(lower) - static_cast<long double>(anchor_)) /
        static_cast<long double>(interval_);
    const long double allowance =
        static_cast<long double>(tolerance(anchor_, lower) / interval_);
    const auto candidate = std::ceil(ratio - allowance);
    if (candidate > static_cast<long double>(
                        std::numeric_limits<WVOutputScheduleOrdinal>::max()))
      return {WVKernelStatusCode::sizeOverflow,
              "Output schedule ordinal exceeds int64 capacity."};
    auto ordinal = candidate <= 0.0L
                       ? WVOutputScheduleOrdinal{0}
                       : static_cast<WVOutputScheduleOrdinal>(candidate);
    if (cursor.committedOrdinal >= ordinal) {
      if (cursor.committedOrdinal ==
          std::numeric_limits<WVOutputScheduleOrdinal>::max())
        return {WVKernelStatusCode::sizeOverflow,
                "Output schedule ordinal overflowed int64 capacity."};
      ordinal = cursor.committedOrdinal + 1;
    }
    const auto time =
        std::fma(static_cast<double>(ordinal), interval_, anchor_);
    if (!std::isfinite(time))
      return {WVKernelStatusCode::sizeOverflow,
              "Output schedule time overflowed finite precision."};
    if (time > upper + tolerance(time, upper))
      return WVKernelStatus::ok();
    occurrence = {};
    occurrence.scheduledTime = time;
    occurrence.ordinal = ordinal;
    occurrence.proposedCursor = {ordinal, {}};
    status = occurrence.payload.reset(payloadSchema());
    if (!status)
      return status;
    occurrence.cursorIdentity = static_cast<std::uint64_t>(ordinal) + 1ULL;
    available = true;
    return WVKernelStatus::ok();
  }

  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }

private:
  double anchor_ = 0.0;
  double interval_ = 0.0;
  double final_ = 0.0;
};

} // namespace

WVKernelStatus WVOutputSchedulePayloadSchema::create(
    std::string identifier, std::uint32_t version,
    std::vector<WVOutputSchedulePayloadField> fields,
    WVOutputSchedulePayloadSchema &schema) {
  if (identifier.empty() || version == 0)
    return invalid("An occurrence-payload schema needs an identity and version.");
  try {
    WVOutputSchedulePayloadSchema candidate;
    candidate.identifier_ = std::move(identifier);
    candidate.version_ = version;
    candidate.slots_.reserve(fields.size());
    std::set<std::string> names;
    std::size_t offset = 0;
    std::uint64_t fingerprint = fnvOffset;
    hashBytes(fingerprint, candidate.identifier_.data(),
              candidate.identifier_.size());
    hashBytes(fingerprint, &candidate.version_, sizeof(candidate.version_));
    for (auto &field : fields) {
      if (field.name.empty() || !names.insert(field.name).second)
        return invalid("Occurrence-payload field names must be nonempty and unique.");
      std::size_t count = 1;
      for (const auto extent : field.dimensions) {
        if (extent == 0 ||
            count > std::numeric_limits<std::size_t>::max() / extent)
          return {WVKernelStatusCode::sizeOverflow,
                  "An occurrence-payload field extent is invalid."};
        count *= extent;
      }
      const auto width = scalarBytes(field.scalarType);
      if (width == 0 || count > WVMaximumOutputSchedulePayloadBytes / width)
        return {WVKernelStatusCode::sizeOverflow,
                "An occurrence-payload field exceeds 4 KiB."};
      offset = alignedOffset(offset, width);
      const auto bytes = count * width;
      if (offset > WVMaximumOutputSchedulePayloadBytes - bytes)
        return {WVKernelStatusCode::sizeOverflow,
                "An encoded occurrence payload exceeds 4 KiB."};
      hashBytes(fingerprint, field.name.data(), field.name.size());
      hashBytes(fingerprint, &field.scalarType, sizeof(field.scalarType));
      for (const auto extent : field.dimensions)
        hashBytes(fingerprint, &extent, sizeof(extent));
      candidate.slots_.push_back(
          {std::move(field.name), field.scalarType,
           std::move(field.dimensions), count, offset, bytes});
      offset += bytes;
    }
    candidate.payloadBytes_ = offset;
    candidate.fingerprint_ = fingerprint;
    schema = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate an occurrence-payload schema."};
  }
}

std::size_t WVOutputSchedulePayloadSchema::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) + identifier_.capacity() +
                      slots_.capacity() *
                          sizeof(WVResolvedOutputSchedulePayloadSlot);
  for (const auto &slot : slots_)
    bytes += slot.name.capacity() +
             slot.dimensions.capacity() * sizeof(std::size_t);
  return bytes;
}

WVKernelStatus WVOutputSchedulePayload::reset(
    const WVOutputSchedulePayloadSchema &schema) noexcept {
  if (schema.payloadBytes() > bytes_.size())
    return {WVKernelStatusCode::sizeOverflow,
            "An encoded occurrence payload exceeds 4 KiB."};
  byteCount_ = schema.payloadBytes();
  schemaFingerprint_ = schema.fingerprint();
  std::fill_n(bytes_.begin(), byteCount_, std::uint8_t{0});
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::setReal(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    const double *values, std::size_t count) noexcept {
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::real64, count);
  if (!status)
    return status;
  if (count != 0 && values == nullptr)
    return invalid("A numeric occurrence-payload value is missing.");
  for (std::size_t index = 0; index < count; ++index)
    if (!std::isfinite(values[index]))
      return invalid("Occurrence-payload numeric values must be finite.");
  std::memcpy(bytes_.data() + schema.slots()[slot].byteOffset, values,
              count * sizeof(double));
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::setInteger(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    const std::int64_t *values, std::size_t count) noexcept {
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::integer64, count);
  if (!status)
    return status;
  if (count != 0 && values == nullptr)
    return invalid("An integer occurrence-payload value is missing.");
  std::memcpy(bytes_.data() + schema.slots()[slot].byteOffset, values,
              count * sizeof(std::int64_t));
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::setBoolean(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    const std::uint8_t *values, std::size_t count) noexcept {
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::boolean8, count);
  if (!status)
    return status;
  if (count != 0 && values == nullptr)
    return invalid("A Boolean occurrence-payload value is missing.");
  for (std::size_t index = 0; index < count; ++index)
    if (values[index] > 1)
      return invalid("Occurrence-payload Boolean values must be zero or one.");
  std::memcpy(bytes_.data() + schema.slots()[slot].byteOffset, values, count);
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::real(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    WVOutputSchedulePayloadRealView &view) const noexcept {
  const auto count = slot < schema.slotCount()
                         ? schema.slots()[slot].elementCount
                         : 0;
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::real64, count);
  if (!status)
    return status;
  view = {reinterpret_cast<const double *>(
              bytes_.data() + schema.slots()[slot].byteOffset),
          count};
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::integer(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    WVOutputSchedulePayloadIntegerView &view) const noexcept {
  const auto count = slot < schema.slotCount()
                         ? schema.slots()[slot].elementCount
                         : 0;
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::integer64, count);
  if (!status)
    return status;
  view = {reinterpret_cast<const std::int64_t *>(
              bytes_.data() + schema.slots()[slot].byteOffset),
          count};
  return WVKernelStatus::ok();
}

WVKernelStatus WVOutputSchedulePayload::boolean(
    const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
    WVOutputSchedulePayloadBooleanView &view) const noexcept {
  const auto count = slot < schema.slotCount()
                         ? schema.slots()[slot].elementCount
                         : 0;
  auto status = validatePayloadAccess(
      schema, *this, slot, WVOutputSchedulePayloadType::boolean8, count);
  if (!status)
    return status;
  view = {bytes_.data() + schema.slots()[slot].byteOffset, count};
  return WVKernelStatus::ok();
}

std::uint64_t WVOutputSchedulePayload::valueFingerprint() const noexcept {
  std::uint64_t result = fnvOffset;
  hashBytes(result, &schemaFingerprint_, sizeof(schemaFingerprint_));
  hashBytes(result, bytes_.data(), byteCount_);
  return result;
}

bool WVOutputSchedulePayload::sameValue(
    const WVOutputSchedulePayload &other) const noexcept {
  return schemaFingerprint_ == other.schemaFingerprint_ &&
         byteCount_ == other.byteCount_ &&
         std::memcmp(bytes_.data(), other.bytes_.data(), byteCount_) == 0;
}

bool sameOutputSchedulePayloadSchema(
    const WVOutputSchedulePayloadSchema &first,
    const WVOutputSchedulePayloadSchema &second) noexcept {
  if (first.identifier() != second.identifier() ||
      first.version() != second.version() ||
      first.payloadBytes() != second.payloadBytes() ||
      first.slotCount() != second.slotCount())
    return false;
  for (std::size_t index = 0; index < first.slotCount(); ++index) {
    const auto &left = first.slots()[index];
    const auto &right = second.slots()[index];
    if (left.name != right.name || left.scalarType != right.scalarType ||
        left.dimensions != right.dimensions ||
        left.elementCount != right.elementCount ||
        left.byteOffset != right.byteOffset || left.byteCount != right.byteCount)
      return false;
  }
  return true;
}

const WVOutputSchedulePayloadSchema &emptyOutputSchedulePayloadSchema() {
  return makeEmptyPayloadSchema();
}

std::shared_ptr<const WVOutputSchedule> makeEvenlySpacedOutputSchedule(
    const WVOutputScheduleRecord &record, WVKernelStatus &status) {
  if (!std::isfinite(record.initialTime) ||
      !std::isfinite(record.outputInterval) || record.outputInterval <= 0.0 ||
      std::isnan(record.finalTime) || record.finalTime < record.initialTime) {
    status = invalid("The evenly-spaced output schedule is invalid.");
    return {};
  }
  if (!record.configuration.schemaIdentifier.empty() ||
      record.configuration.schemaVersion != 0 ||
      !record.configuration.values.empty()) {
    status = invalid("A legacy evenly-spaced schedule cannot carry a typed "
                     "configuration record.");
    return {};
  }
  if (record.finalTime > record.initialTime) {
    const long double ratio = (static_cast<long double>(record.finalTime) -
                               static_cast<long double>(record.initialTime)) /
                              static_cast<long double>(record.outputInterval);
    const auto lastCandidate = std::floor(ratio);
    if (lastCandidate >= 1.0L &&
        lastCandidate <=
            static_cast<long double>(
                std::numeric_limits<WVOutputScheduleOrdinal>::max())) {
      const auto last = static_cast<WVOutputScheduleOrdinal>(lastCandidate);
      const auto current = std::fma(static_cast<double>(last),
                                    record.outputInterval, record.initialTime);
      const auto previous = std::fma(static_cast<double>(last - 1),
                                     record.outputInterval, record.initialTime);
      if (!std::isfinite(current) || !(current > previous)) {
        status = invalid("Output schedule ordinals are not distinguishable in "
                         "finite precision.");
        return {};
      }
    }
  }
  try {
    status = WVKernelStatus::ok();
    return std::make_shared<WVEvenlySpacedSchedule>(record);
  } catch (const std::bad_alloc &) {
    status = {WVKernelStatusCode::allocationFailure,
              "Unable to allocate an output schedule."};
    return {};
  }
}

} // namespace wavevortex::runtime
