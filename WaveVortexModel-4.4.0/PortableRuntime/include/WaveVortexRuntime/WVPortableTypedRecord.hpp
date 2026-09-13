#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <variant>
#include <vector>

namespace wavevortex::runtime {

enum class WVPortableValueType : std::uint8_t { boolean, integer, real, text };

using WVPortableValueStorage =
    std::variant<std::vector<std::uint8_t>, std::vector<std::int64_t>,
                 std::vector<double>, std::vector<std::string>>;

// Construction-only named data. Empty dimensions denote a scalar; otherwise
// their product is the number of stored values. Runtime implementations resolve
// these records to concrete typed state before integration begins.
struct WVPortableNamedValue {
  std::string name;
  std::vector<std::size_t> dimensions;
  WVPortableValueStorage storage = std::vector<double>{};

  WVPortableValueType valueType() const noexcept;
  std::size_t valueCount() const noexcept;
  std::size_t encodedBytes() const noexcept;
  std::size_t persistentBytes() const noexcept;
};

struct WVPortableTypedRecord {
  std::string schemaIdentifier;
  std::uint32_t schemaVersion = 0;
  std::vector<WVPortableNamedValue> values;

  const WVPortableNamedValue *value(const std::string &name) const noexcept;
  std::size_t encodedBytes() const noexcept;
  std::size_t persistentBytes() const noexcept;
};

struct WVPortableTypedRecordValidation {
  std::size_t maximumEncodedBytes = std::numeric_limits<std::size_t>::max();
  bool requireFiniteReals = false;
  bool allowText = true;
};

WVKernelStatus
validatePortableTypedRecord(const WVPortableTypedRecord &record,
                            WVPortableTypedRecordValidation validation = {});

WVKernelStatus encodePortableTypedRecord(const WVPortableTypedRecord &record,
                                         std::vector<std::uint8_t> &bytes);
WVKernelStatus
decodePortableTypedRecord(const std::vector<std::uint8_t> &bytes,
                          WVPortableTypedRecord &record,
                          WVPortableTypedRecordValidation validation = {});

// Exact, allocation-free value comparison for already validated construction
// or cursor-boundary records. Named event code may carry these records but must
// not look values up by string.
bool samePortableTypedRecordValue(const WVPortableTypedRecord &left,
                                  const WVPortableTypedRecord &right) noexcept;
bool isCanonicalEmptyPortableTypedRecord(
    const WVPortableTypedRecord &record) noexcept;

} // namespace wavevortex::runtime
