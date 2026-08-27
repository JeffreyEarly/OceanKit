#pragma once

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace wavevortex::runtime::test {

inline constexpr const char *quadraticScheduleType =
    "WVTestQuadraticOutputSchedule";

inline double scalarReal(const WVPortableTypedRecord &record,
                         const char *name) {
  const auto *value = record.value(name);
  if (value == nullptr || value->valueType() != WVPortableValueType::real ||
      value->valueCount() != 1)
    return std::numeric_limits<double>::quiet_NaN();
  return std::get<std::vector<double>>(value->storage)[0];
}

class WVTestQuadraticSchedule final : public WVOutputSchedule {
public:
  WVTestQuadraticSchedule(double anchor, double scale, double finalTime)
      : anchor_(anchor), scale_(scale), finalTime_(finalTime) {}
  const char *typeIdentifier() const noexcept override {
    return quadraticScheduleType;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  const WVOutputSchedulePayloadSchema &payloadSchema() const noexcept override {
    return emptyOutputSchedulePayloadSchema();
  }
  WVKernelStatus
  validateCursor(const WVOutputScheduleCursor &cursor) const override {
    if (cursor.committedOrdinal < -1)
      return {WVKernelStatusCode::invalidConfiguration,
              "quadratic cursor ordinal"};
    if (cursor.values.schemaIdentifier.empty())
      return cursor.committedOrdinal == WVNoCommittedOutputOrdinal
                 ? WVKernelStatus::ok()
                 : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                                  "quadratic committed cursor payload"};
    const auto status = validatePortableTypedRecord(
        cursor.values, {WVMaximumOutputScheduleCursorBytes, true, false});
    if (!status)
      return status;
    const auto *next = cursor.values.value("nextOrdinal");
    if (cursor.values.schemaIdentifier != "quadratic-cursor-v1" ||
        cursor.values.schemaVersion != 1 || next == nullptr ||
        next->valueType() != WVPortableValueType::integer ||
        next->valueCount() != 1 ||
        std::get<std::vector<std::int64_t>>(next->storage)[0] !=
            cursor.committedOrdinal + 1)
      return {WVKernelStatusCode::invalidConfiguration,
              "quadratic cursor payload"};
    return WVKernelStatus::ok();
  }
  WVKernelStatus committedTime(const WVOutputScheduleCursor &cursor,
                               double &time, bool &available) const override {
    const auto status = validateCursor(cursor);
    if (!status)
      return status;
    available = cursor.committedOrdinal >= 0;
    if (available) {
      const auto n = static_cast<double>(cursor.committedOrdinal);
      time = std::fma(scale_, n * n, anchor_);
    }
    return WVKernelStatus::ok();
  }
  WVKernelStatus peek(const WVOutputScheduleCursor &cursor, double lowerBound,
                      double upperBound, WVOutputScheduleOccurrence &occurrence,
                      bool &available) const override {
    const auto status = validateCursor(cursor);
    if (!status)
      return status;
    available = false;
    auto next = cursor.committedOrdinal + 1;
    if (lowerBound > anchor_) {
      const auto estimate = static_cast<WVOutputScheduleOrdinal>(
          std::ceil(std::sqrt((lowerBound - anchor_) / scale_)));
      next = std::max(next, estimate);
    }
    const auto n = static_cast<double>(next);
    const auto time = std::fma(scale_, n * n, anchor_);
    if (!std::isfinite(time) || time > std::min(upperBound, finalTime_))
      return WVKernelStatus::ok();
    WVPortableTypedRecord proposed;
    proposed.schemaIdentifier = "quadratic-cursor-v1";
    proposed.schemaVersion = 1;
    proposed.values.push_back(
        {"nextOrdinal", {}, std::vector<std::int64_t>{next + 1}});
    occurrence = {};
    occurrence.scheduledTime = time;
    occurrence.ordinal = next;
    occurrence.proposedCursor = {next, std::move(proposed)};
    const auto payloadStatus = occurrence.payload.reset(payloadSchema());
    if (!payloadStatus)
      return payloadStatus;
    occurrence.cursorIdentity = static_cast<std::uint64_t>(next) + 1ULL;
    available = true;
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }

private:
  double anchor_;
  double scale_;
  double finalTime_;
};

inline std::shared_ptr<const WVOutputSchedule>
makeQuadraticSchedule(const WVOutputScheduleRecord &record,
                      WVKernelStatus &status) {
  const auto validation =
      validatePortableTypedRecord(record.configuration, {4096, true, false});
  const auto anchor = scalarReal(record.configuration, "anchor");
  const auto scale = scalarReal(record.configuration, "scale");
  const auto finalTime = scalarReal(record.configuration, "finalTime");
  if (!validation || !std::isfinite(anchor) || !std::isfinite(scale) ||
      scale <= 0.0 || !std::isfinite(finalTime) || finalTime < anchor) {
    status = {WVKernelStatusCode::invalidConfiguration,
              "invalid quadratic schedule configuration"};
    return {};
  }
  status = WVKernelStatus::ok();
  return std::make_shared<WVTestQuadraticSchedule>(anchor, scale, finalTime);
}

inline WVOutputScheduleRecord quadraticSchedule(double finalTime,
                                                double initialTime = 0.0,
                                                double scale = 1.0) {
  WVOutputScheduleRecord result;
  result.typeIdentifier = quadraticScheduleType;
  result.contractVersion = 1;
  result.configuration.schemaIdentifier = "quadratic-schedule-v1";
  result.configuration.schemaVersion = 1;
  result.configuration.values = {
      {"anchor", {}, std::vector<double>{initialTime}},
      {"scale", {}, std::vector<double>{scale}},
      {"finalTime", {}, std::vector<double>{finalTime}}};
  return result;
}

inline WVKernelStatus registerQuadraticSchedule(
    WVExtensionCatalogBuilder &builder) {
  return builder.addOutputScheduleFactory(
      {quadraticScheduleType, 1, &makeQuadraticSchedule});
}

} // namespace wavevortex::runtime::test
