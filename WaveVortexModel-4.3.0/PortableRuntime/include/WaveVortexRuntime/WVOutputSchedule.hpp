#pragma once

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

using WVOutputScheduleOrdinal = std::int64_t;
inline constexpr WVOutputScheduleOrdinal WVNoCommittedOutputOrdinal = -1;
inline constexpr std::size_t WVMaximumOutputScheduleCursorBytes = 4096;
inline constexpr std::size_t WVMaximumOutputSchedulePayloadBytes = 4096;
inline constexpr const char *WVEvenlySpacedOutputScheduleType =
    "WVEvenlySpacedOutputSchedule";
inline constexpr const char *WVStateTriggeredOutputScheduleType =
    "WVStateTriggeredOutputSchedule";

struct WVOutputScheduleCursor {
  WVOutputScheduleOrdinal committedOrdinal = WVNoCommittedOutputOrdinal;
  WVPortableTypedRecord values;
};

enum class WVOutputSchedulePayloadType : std::uint8_t {
  real64,
  integer64,
  boolean8
};

// Construction-time declaration for one compact event-payload field. Names
// and dimensions are metadata only; schedules and observers use resolved slot
// ordinals after construction.
struct WVOutputSchedulePayloadField {
  std::string name;
  WVOutputSchedulePayloadType scalarType =
      WVOutputSchedulePayloadType::real64;
  std::vector<std::size_t> dimensions;
};

struct WVResolvedOutputSchedulePayloadSlot {
  std::string name;
  WVOutputSchedulePayloadType scalarType =
      WVOutputSchedulePayloadType::real64;
  std::vector<std::size_t> dimensions;
  std::size_t elementCount = 0;
  std::size_t byteOffset = 0;
  std::size_t byteCount = 0;
};

class WVOutputSchedulePayloadSchema final {
public:
  static WVKernelStatus
  create(std::string identifier, std::uint32_t version,
         std::vector<WVOutputSchedulePayloadField> fields,
         WVOutputSchedulePayloadSchema &schema);

  const std::string &identifier() const noexcept { return identifier_; }
  std::uint32_t version() const noexcept { return version_; }
  const std::vector<WVResolvedOutputSchedulePayloadSlot> &slots() const
      noexcept {
    return slots_;
  }
  std::size_t slotCount() const noexcept { return slots_.size(); }
  std::size_t payloadBytes() const noexcept { return payloadBytes_; }
  std::uint64_t fingerprint() const noexcept { return fingerprint_; }
  std::size_t persistentBytes() const noexcept;

private:
  std::string identifier_;
  std::uint32_t version_ = 0;
  std::vector<WVResolvedOutputSchedulePayloadSlot> slots_;
  std::size_t payloadBytes_ = 0;
  std::uint64_t fingerprint_ = 0;
};

struct WVOutputSchedulePayloadRealView {
  const double *data = nullptr;
  std::size_t count = 0;
};

struct WVOutputSchedulePayloadIntegerView {
  const std::int64_t *data = nullptr;
  std::size_t count = 0;
};

struct WVOutputSchedulePayloadBooleanView {
  const std::uint8_t *data = nullptr;
  std::size_t count = 0;
};

// Fixed-capacity hot-path payload. It performs no heap allocation and contains
// no names. The schema resolves every name, type, shape, and aligned byte
// offset once while the schedule is constructed.
class WVOutputSchedulePayload final {
public:
  WVKernelStatus reset(const WVOutputSchedulePayloadSchema &schema) noexcept;
  WVKernelStatus setReal(const WVOutputSchedulePayloadSchema &schema,
                         std::size_t slot, const double *values,
                         std::size_t count) noexcept;
  WVKernelStatus setInteger(const WVOutputSchedulePayloadSchema &schema,
                            std::size_t slot, const std::int64_t *values,
                            std::size_t count) noexcept;
  WVKernelStatus setBoolean(const WVOutputSchedulePayloadSchema &schema,
                            std::size_t slot, const std::uint8_t *values,
                            std::size_t count) noexcept;
  WVKernelStatus real(const WVOutputSchedulePayloadSchema &schema,
                      std::size_t slot,
                      WVOutputSchedulePayloadRealView &view) const noexcept;
  WVKernelStatus
  integer(const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
          WVOutputSchedulePayloadIntegerView &view) const noexcept;
  WVKernelStatus
  boolean(const WVOutputSchedulePayloadSchema &schema, std::size_t slot,
          WVOutputSchedulePayloadBooleanView &view) const noexcept;

  std::size_t byteCount() const noexcept { return byteCount_; }
  std::uint64_t schemaFingerprint() const noexcept {
    return schemaFingerprint_;
  }
  std::uint64_t valueFingerprint() const noexcept;
  bool sameValue(const WVOutputSchedulePayload &other) const noexcept;

private:
  alignas(double)
      std::array<std::uint8_t, WVMaximumOutputSchedulePayloadBytes> bytes_{};
  std::size_t byteCount_ = 0;
  std::uint64_t schemaFingerprint_ = 0;
};

bool sameOutputSchedulePayloadSchema(
    const WVOutputSchedulePayloadSchema &first,
    const WVOutputSchedulePayloadSchema &second) noexcept;
const WVOutputSchedulePayloadSchema &emptyOutputSchedulePayloadSchema();

struct WVOutputScheduleOccurrence {
  double scheduledTime = 0.0;
  WVOutputScheduleOrdinal ordinal = WVNoCommittedOutputOrdinal;
  WVOutputScheduleCursor proposedCursor;
  WVOutputSchedulePayload payload;
  // Provider-resolved identity of the complete proposed cursor state. This is
  // compared only at event granularity and avoids named cursor inspection in
  // the observation hot path.
  std::uint64_t cursorIdentity = 0;
};

// Stable source API v1 boundary for one exact identity/version schedule.
// Implementations are immutable and retained by shared ownership after
// construction. peek() is transactional: the caller owns and commits the
// proposed cursor only after the corresponding output route is delivered.
// Schedule identity/version, configuration, payload-schema identity/version,
// and typed cursor are independent exact data contracts rather than
// consequences of the C++ source API version.
class WVOutputSchedule {
public:
  virtual ~WVOutputSchedule() = default;
  virtual const char *typeIdentifier() const noexcept = 0;
  virtual std::uint32_t contractVersion() const noexcept = 0;
  virtual const WVOutputSchedulePayloadSchema &payloadSchema() const noexcept =
      0;
  virtual WVKernelStatus
  validateCursor(const WVOutputScheduleCursor &cursor) const = 0;
  virtual WVKernelStatus committedTime(const WVOutputScheduleCursor &cursor,
                                       double &time, bool &available) const = 0;
  virtual WVKernelStatus peek(const WVOutputScheduleCursor &cursor,
                              double lowerBound, double upperBound,
                              WVOutputScheduleOccurrence &occurrence,
                              bool &available) const = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
};

using WVOutputScheduleFactory = std::shared_ptr<const WVOutputSchedule> (*)(
    const WVOutputScheduleRecord &, WVKernelStatus &);

std::shared_ptr<const WVOutputSchedule>
makeEvenlySpacedOutputSchedule(const WVOutputScheduleRecord &record,
                              WVKernelStatus &status);

} // namespace wavevortex::runtime
