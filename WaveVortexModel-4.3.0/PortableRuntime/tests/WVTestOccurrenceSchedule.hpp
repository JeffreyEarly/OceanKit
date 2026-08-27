#pragma once

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime::test {

inline constexpr const char *occurrenceScheduleType =
    "WVTestOccurrenceOutputSchedule";
inline constexpr const char *occurrenceScheduleConfigurationIdentifier =
    "wv-test-occurrence-schedule-v1";
inline constexpr const char *occurrencePayloadSchemaIdentifier =
    "wv-test-occurrence-payload-v1";

struct WVTestOccurrenceScheduleCounters {
  std::size_t constructionCount = 0;
  std::size_t peekCount = 0;
};

namespace occurrence_schedule_detail {

inline std::atomic<std::size_t> constructionCount{0};
inline std::atomic<std::size_t> peekCount{0};

inline WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

inline const WVPortableNamedValue *scalar(
    const WVPortableTypedRecord &record, const char *name,
    WVPortableValueType type) noexcept {
  const auto *value = record.value(name);
  return value != nullptr && value->valueType() == type &&
                 value->valueCount() == 1 && value->dimensions.empty()
             ? value
             : nullptr;
}

inline std::uint64_t occurrenceIdentity(
    WVOutputScheduleOrdinal ordinal,
    const WVOutputSchedulePayload &payload) noexcept {
  auto value = payload.valueFingerprint();
  value ^= static_cast<std::uint64_t>(ordinal) + 0x9e3779b97f4a7c15ULL +
           (value << 6U) + (value >> 2U);
  return value;
}

} // namespace occurrence_schedule_detail

inline void resetOccurrenceScheduleCounters() noexcept {
  occurrence_schedule_detail::constructionCount.store(
      0, std::memory_order_relaxed);
  occurrence_schedule_detail::peekCount.store(0, std::memory_order_relaxed);
}

inline WVTestOccurrenceScheduleCounters
occurrenceScheduleCounters() noexcept {
  return {occurrence_schedule_detail::constructionCount.load(
              std::memory_order_relaxed),
          occurrence_schedule_detail::peekCount.load(
              std::memory_order_relaxed)};
}

inline WVKernelStatus
createOccurrencePayloadSchema(WVOutputSchedulePayloadSchema &schema) {
  return WVOutputSchedulePayloadSchema::create(
      occurrencePayloadSchemaIdentifier, 1,
      {{"sourceReal", WVOutputSchedulePayloadType::real64, {}},
       {"sourceInteger", WVOutputSchedulePayloadType::integer64, {}},
       {"sourceBoolean", WVOutputSchedulePayloadType::boolean8, {}}},
      schema);
}

class WVTestOccurrenceSchedule final : public WVOutputSchedule {
public:
  WVTestOccurrenceSchedule(double interval, double initialTime,
                           double finalTime, double sourceReal,
                           std::int64_t sourceInteger,
                           std::uint8_t sourceBoolean,
                           WVOutputSchedulePayloadSchema payloadSchema)
      : interval_(interval), initialTime_(initialTime), finalTime_(finalTime),
        sourceReal_(sourceReal), sourceInteger_(sourceInteger),
        sourceBoolean_(sourceBoolean),
        payloadSchema_(std::move(payloadSchema)) {}

  const char *typeIdentifier() const noexcept override {
    return occurrenceScheduleType;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  const WVOutputSchedulePayloadSchema &payloadSchema() const noexcept override {
    return payloadSchema_;
  }

  WVKernelStatus
  validateCursor(const WVOutputScheduleCursor &cursor) const override {
    if (cursor.committedOrdinal < WVNoCommittedOutputOrdinal)
      return occurrence_schedule_detail::invalid(
          "test occurrence cursor ordinal is less than -1");
    if (!cursor.values.schemaIdentifier.empty() ||
        cursor.values.schemaVersion != 0 || !cursor.values.values.empty())
      return occurrence_schedule_detail::invalid(
          "test occurrence cursor contains unexpected provider state");
    if (cursor.committedOrdinal >= 0) {
      const auto time =
          std::fma(static_cast<double>(cursor.committedOrdinal), interval_,
                   initialTime_);
      if (!std::isfinite(time) || time > finalTime_)
        return occurrence_schedule_detail::invalid(
            "test occurrence cursor lies beyond the source schedule");
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus committedTime(const WVOutputScheduleCursor &cursor,
                               double &time,
                               bool &available) const override {
    const auto status = validateCursor(cursor);
    if (!status)
      return status;
    available = cursor.committedOrdinal >= 0;
    if (available)
      time = std::fma(static_cast<double>(cursor.committedOrdinal), interval_,
                      initialTime_);
    return WVKernelStatus::ok();
  }

  WVKernelStatus peek(const WVOutputScheduleCursor &cursor, double lowerBound,
                      double upperBound, WVOutputScheduleOccurrence &occurrence,
                      bool &available) const override {
    occurrence_schedule_detail::peekCount.fetch_add(
        1, std::memory_order_relaxed);
    available = false;
    auto status = validateCursor(cursor);
    if (!status)
      return status;
    if (!std::isfinite(lowerBound) || !std::isfinite(upperBound) ||
        upperBound < lowerBound)
      return occurrence_schedule_detail::invalid(
          "test occurrence bounds must be finite and nondecreasing");

    if (cursor.committedOrdinal ==
        std::numeric_limits<WVOutputScheduleOrdinal>::max())
      return {WVKernelStatusCode::sizeOverflow,
              "test occurrence ordinal exceeds int64 capacity"};
    auto ordinal = cursor.committedOrdinal + 1;
    if (lowerBound > initialTime_) {
      const auto estimate = std::ceil((lowerBound - initialTime_) / interval_);
      if (estimate > static_cast<double>(
                         std::numeric_limits<WVOutputScheduleOrdinal>::max()))
        return {WVKernelStatusCode::sizeOverflow,
                "test occurrence ordinal exceeds int64 capacity"};
      ordinal = std::max(ordinal,
                         static_cast<WVOutputScheduleOrdinal>(estimate));
    }
    const auto time = std::fma(static_cast<double>(ordinal), interval_,
                               initialTime_);
    if (!std::isfinite(time))
      return {WVKernelStatusCode::sizeOverflow,
              "test occurrence time is not finite"};
    if (time > std::min(upperBound, finalTime_))
      return WVKernelStatus::ok();
    if (sourceInteger_ >
        std::numeric_limits<std::int64_t>::max() - ordinal)
      return {WVKernelStatusCode::sizeOverflow,
              "test occurrence integer payload overflowed"};

    occurrence = {};
    occurrence.scheduledTime = time;
    occurrence.ordinal = ordinal;
    occurrence.proposedCursor = {ordinal, {}};
    status = occurrence.payload.reset(payloadSchema_);
    if (!status)
      return status;
    const auto realValue =
        std::fma(static_cast<double>(ordinal), 0.5, sourceReal_);
    const auto integerValue = sourceInteger_ + ordinal;
    const std::uint8_t booleanValue = static_cast<std::uint8_t>(
        sourceBoolean_ ^ static_cast<std::uint8_t>(ordinal & 1));
    status = occurrence.payload.setReal(payloadSchema_, 0, &realValue, 1);
    if (status)
      status = occurrence.payload.setInteger(payloadSchema_, 1,
                                             &integerValue, 1);
    if (status)
      status = occurrence.payload.setBoolean(payloadSchema_, 2,
                                             &booleanValue, 1);
    if (!status)
      return status;
    occurrence.cursorIdentity =
        occurrence_schedule_detail::occurrenceIdentity(ordinal,
                                                       occurrence.payload);
    available = true;
    return WVKernelStatus::ok();
  }

  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + payloadSchema_.persistentBytes() -
           sizeof(payloadSchema_);
  }

private:
  double interval_ = 0.0;
  double initialTime_ = 0.0;
  double finalTime_ = 0.0;
  double sourceReal_ = 0.0;
  std::int64_t sourceInteger_ = 0;
  std::uint8_t sourceBoolean_ = 0;
  WVOutputSchedulePayloadSchema payloadSchema_;
};

inline std::shared_ptr<const WVOutputSchedule>
makeOccurrenceSchedule(const WVOutputScheduleRecord &record,
                       WVKernelStatus &status) {
  using namespace occurrence_schedule_detail;
  const auto validation =
      validatePortableTypedRecord(record.configuration, {4096, true, false});
  const auto *real = scalar(record.configuration, "sourceReal",
                            WVPortableValueType::real);
  const auto *integer = scalar(record.configuration, "sourceInteger",
                               WVPortableValueType::integer);
  const auto *boolean = scalar(record.configuration, "sourceBoolean",
                               WVPortableValueType::boolean);
  if (!validation ||
      record.configuration.schemaIdentifier !=
          occurrenceScheduleConfigurationIdentifier ||
      record.configuration.schemaVersion != 1 ||
      record.configuration.values.size() != 3 || real == nullptr ||
      integer == nullptr || boolean == nullptr ||
      !std::isfinite(record.outputInterval) || record.outputInterval <= 0.0 ||
      !std::isfinite(record.initialTime) || !std::isfinite(record.finalTime) ||
      record.finalTime < record.initialTime) {
    status = invalid("invalid test occurrence schedule configuration");
    return {};
  }
  const auto sourceReal = std::get<std::vector<double>>(real->storage)[0];
  const auto sourceInteger =
      std::get<std::vector<std::int64_t>>(integer->storage)[0];
  const auto sourceBoolean =
      std::get<std::vector<std::uint8_t>>(boolean->storage)[0];
  if (!std::isfinite(sourceReal) || sourceBoolean > 1) {
    status = invalid("invalid test occurrence schedule source value");
    return {};
  }
  WVOutputSchedulePayloadSchema payloadSchema;
  status = createOccurrencePayloadSchema(payloadSchema);
  if (!status)
    return {};
  try {
    auto result = std::make_shared<WVTestOccurrenceSchedule>(
        record.outputInterval, record.initialTime, record.finalTime, sourceReal,
        sourceInteger, sourceBoolean, std::move(payloadSchema));
    constructionCount.fetch_add(1, std::memory_order_relaxed);
    status = WVKernelStatus::ok();
    return result;
  } catch (const std::bad_alloc &) {
    status = {WVKernelStatusCode::allocationFailure,
              "unable to allocate a test occurrence schedule"};
    return {};
  }
}

inline WVOutputScheduleRecord occurrenceSchedule(
    double sourceReal, std::int64_t sourceInteger,
    std::uint8_t sourceBoolean, double finalTime = 1.0,
    double initialTime = 0.0, double interval = 1.0) {
  WVOutputScheduleRecord result;
  result.outputInterval = interval;
  result.initialTime = initialTime;
  result.finalTime = finalTime;
  result.typeIdentifier = occurrenceScheduleType;
  result.contractVersion = 1;
  result.configuration.schemaIdentifier =
      occurrenceScheduleConfigurationIdentifier;
  result.configuration.schemaVersion = 1;
  result.configuration.values = {
      {"sourceReal", {}, std::vector<double>{sourceReal}},
      {"sourceInteger", {}, std::vector<std::int64_t>{sourceInteger}},
      {"sourceBoolean", {}, std::vector<std::uint8_t>{sourceBoolean}}};
  return result;
}

inline WVKernelStatus registerOccurrenceSchedule(
    WVExtensionCatalogBuilder &builder) {
  return builder.addOutputScheduleFactory(
      {occurrenceScheduleType, 1, &makeOccurrenceSchedule});
}

} // namespace wavevortex::runtime::test
