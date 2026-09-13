#pragma once

#include "WVTestOccurrenceSchedule.hpp"

#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime::test {

inline constexpr const char *fixedGeometryObservationType =
    "WVTestFixedGeometryObservation";
inline constexpr const char *eventVariableOneDimensionalGeometryType =
    "WVTestEventVariableOneDimensionalGeometry";
inline constexpr const char *variableSamplesByFixedDepthBinsType =
    "WVTestVariableSamplesByFixedDepthBins";
inline constexpr const char *stateCoupledIrregularGeometryType =
    "WVTestStateCoupledIrregularGeometry";
inline constexpr const char *nestedRaggedGeometryType =
    "WVTestNestedRaggedGeometry";

enum class WVTestObservationOccurrenceProviderKind : std::size_t {
  fixedGeometry,
  eventVariableOneDimensionalGeometry,
  variableSamplesByFixedDepthBins,
  stateCoupledIrregularGeometry,
  nestedRaggedGeometry,
  count
};

struct WVTestObservationOccurrenceProviderCounters {
  std::size_t constructionCount = 0;
  std::size_t outputPlanCount = 0;
  std::size_t occurrencePreparationCount = 0;
  std::size_t batchCount = 0;
  // These proof providers are passive. A nonzero count would demonstrate an
  // accidental per-stage observer call.
  std::size_t integrationStageCount = 0;
};

namespace observation_occurrence_provider_detail {

inline constexpr std::size_t providerCount = static_cast<std::size_t>(
    WVTestObservationOccurrenceProviderKind::count);

struct AtomicCounters {
  std::atomic<std::size_t> constructionCount{0};
  std::atomic<std::size_t> outputPlanCount{0};
  std::atomic<std::size_t> occurrencePreparationCount{0};
  std::atomic<std::size_t> batchCount{0};
  std::atomic<std::size_t> integrationStageCount{0};
};

inline std::array<AtomicCounters, providerCount> counters{};

inline WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

inline WVKernelStatus allocationFailure() {
  return {WVKernelStatusCode::allocationFailure,
          "Unable to allocate a test observation occurrence."};
}

template <WVTestObservationOccurrenceProviderKind Kind>
constexpr std::size_t kindIndex() noexcept {
  return static_cast<std::size_t>(Kind);
}

template <WVTestObservationOccurrenceProviderKind Kind>
const std::string &typeIdentifier() noexcept {
  if constexpr (Kind ==
                WVTestObservationOccurrenceProviderKind::fixedGeometry) {
    static const std::string value = fixedGeometryObservationType;
    return value;
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  eventVariableOneDimensionalGeometry) {
    static const std::string value =
        eventVariableOneDimensionalGeometryType;
    return value;
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  variableSamplesByFixedDepthBins) {
    static const std::string value = variableSamplesByFixedDepthBinsType;
    return value;
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  stateCoupledIrregularGeometry) {
    static const std::string value = stateCoupledIrregularGeometryType;
    return value;
  } else {
    static const std::string value = nestedRaggedGeometryType;
    return value;
  }
}

inline std::uint64_t magnitude(std::int64_t value) noexcept {
  return value < 0
             ? static_cast<std::uint64_t>(-(value + 1)) + std::uint64_t{1}
             : static_cast<std::uint64_t>(value);
}

struct PayloadScalars {
  double real = 0.0;
  std::int64_t integer = 0;
  std::uint8_t boolean = 0;
};

inline WVKernelStatus payloadScalars(
    const WVObserverOccurrencePreparationContext &context,
    PayloadScalars &output) {
  if (context.payloadSchema == nullptr || context.payload == nullptr)
    return invalid("Test occurrence payload is absent.");
  WVOutputSchedulePayloadRealView real;
  WVOutputSchedulePayloadIntegerView integer;
  WVOutputSchedulePayloadBooleanView boolean;
  auto status = context.payload->real(*context.payloadSchema, 0, real);
  if (status)
    status = context.payload->integer(*context.payloadSchema, 1, integer);
  if (status)
    status = context.payload->boolean(*context.payloadSchema, 2, boolean);
  if (!status)
    return status;
  if (real.count != 1 || integer.count != 1 || boolean.count != 1 ||
      real.data == nullptr || integer.data == nullptr ||
      boolean.data == nullptr || !std::isfinite(real.data[0]) ||
      boolean.data[0] > 1)
    return invalid("Test occurrence payload is malformed.");
  output = {real.data[0], integer.data[0], boolean.data[0]};
  return WVKernelStatus::ok();
}

inline WVObservationAxis axis(std::string identifier,
                              WVObservationAxisKind kind,
                              std::size_t extent,
                              WVObservationCoordinateRole role =
                                  WVObservationCoordinateRole::none) {
  WVObservationAxis result;
  result.identifier = identifier;
  result.name = std::move(identifier);
  result.kind = kind;
  result.extent = extent;
  result.coordinateRole = role;
  return result;
}

inline WVObservationVariable variable(
    std::string identifier, WVObservationScalarType scalarType,
    std::vector<std::string> dimensions, WVObservationValueLayout layout,
    WVObservationCoordinateRole coordinateRole =
        WVObservationCoordinateRole::none,
    WVObservationRaggedRole raggedRole = WVObservationRaggedRole::none,
    std::string raggedChild = {}) {
  WVObservationVariable result;
  result.identifier = identifier;
  result.name = std::move(identifier);
  result.scalarType = scalarType;
  result.dimensionIdentifiers = std::move(dimensions);
  result.layout = layout;
  result.coordinateRole = coordinateRole;
  result.raggedRole = raggedRole;
  result.raggedChildAxisIdentifier = std::move(raggedChild);
  return result;
}

inline void addOccurrenceValue(
    WVObserverOutputPlan &plan, std::string identifier,
    WVObservationScalarType scalarType, std::vector<std::string> dimensions,
    WVObservationValueLayout layout,
    WVObservationCoordinateRole coordinateRole =
        WVObservationCoordinateRole::none,
    WVObservationRaggedRole raggedRole = WVObservationRaggedRole::none,
    std::string raggedChild = {}) {
  plan.schema.variables.push_back(
      variable(identifier, scalarType, std::move(dimensions), layout,
               coordinateRole, raggedRole, std::move(raggedChild)));
  plan.occurrenceValues.push_back({identifier, 0});
  WVObserverOutputChannel channel;
  channel.variableIdentifier = identifier;
  channel.source = WVObserverOutputChannelSource::occurrenceValue;
  channel.sourceIdentifier = std::move(identifier);
  plan.channels.push_back(std::move(channel));
}

inline void addOccurrenceField(WVObserverOutputPlan &plan,
                               std::string identifier,
                               std::string fieldName,
                               std::vector<std::string> dimensions,
                               std::size_t positionSetSlot,
                               WVPositionInterpolation interpolation =
                                   WVPositionInterpolation::linear,
                               WVObservationValueLayout layout =
                                   WVObservationValueLayout::flat) {
  plan.schema.variables.push_back(variable(
      identifier, WVObservationScalarType::real64, std::move(dimensions),
      layout));
  WVObserverOutputChannel channel;
  channel.variableIdentifier = std::move(identifier);
  channel.source = WVObserverOutputChannelSource::occurrenceField;
  channel.sourceIdentifier = std::move(fieldName);
  channel.positionSetSlot = positionSetSlot;
  channel.sampling.interpolation = interpolation;
  plan.channels.push_back(std::move(channel));
}

inline void addOccurrencePositionSet(WVObserverOutputPlan &plan,
                                     std::string identifier) {
  plan.occurrencePositionSets.push_back(
      {std::move(identifier), "sample-time", "x", "y", "z"});
}

inline WVKernelStatus addPayloadSchema(WVObserverOutputPlan &plan) {
  return createOccurrencePayloadSchema(plan.occurrencePayloadSchema);
}

template <WVTestObservationOccurrenceProviderKind Kind>
WVKernelStatus buildOutputPlan(const WVObserverRecord &observer,
                               const WVObserverOutputPlanningContext &context,
                               WVObserverOutputPlan &output) {
  if (context.configuration == nullptr)
    return invalid("Test occurrence output planning requires a model configuration.");
  WVObserverOutputPlan plan;
  plan.schema.identifier = observer.identifier + "-occurrence-observation-v1";
  plan.schema.version = 1;
  plan.schema.metadata.attributes = {
      {"provider", typeIdentifier<Kind>()},
      {"contract", "generic-event-geometry-v1"}};

  if constexpr (Kind ==
                WVTestObservationOccurrenceProviderKind::fixedGeometry) {
    if (observer.x.empty() || observer.x.size() != observer.y.size() ||
        observer.x.size() != observer.z.size())
      return invalid("Fixed test occurrence coordinates are invalid.");
    plan.occurrencePayloadSchema = emptyOutputSchedulePayloadSchema();
    plan.schema.axes.push_back(axis(
        "sample", WVObservationAxisKind::fixed, observer.x.size(),
        WVObservationCoordinateRole::identifier));
    for (const auto &declaration :
         std::array<std::pair<const char *, WVObservationCoordinateRole>, 4>{{
             {"sample-time", WVObservationCoordinateRole::sampleTime},
             {"x", WVObservationCoordinateRole::x},
             {"y", WVObservationCoordinateRole::y},
             {"z", WVObservationCoordinateRole::z}}})
      addOccurrenceValue(plan, declaration.first,
                         WVObservationScalarType::real64, {"sample"},
                         WVObservationValueLayout::record,
                         declaration.second);
    addOccurrencePositionSet(plan, "fixed-samples");
    addOccurrenceField(plan, "sampled-u", "u", {"sample"}, 0,
                       WVPositionInterpolation::linear,
                       WVObservationValueLayout::record);
    addOccurrenceField(plan, "sampled-ssh", "ssh", {"sample"}, 0,
                       WVPositionInterpolation::linear,
                       WVObservationValueLayout::record);
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  eventVariableOneDimensionalGeometry) {
    auto status = addPayloadSchema(plan);
    if (!status)
      return status;
    plan.schema.axes.push_back(axis(
        "sample", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::identifier));
    for (const auto &declaration :
         std::array<std::pair<const char *, WVObservationCoordinateRole>, 4>{{
             {"sample-time", WVObservationCoordinateRole::sampleTime},
             {"x", WVObservationCoordinateRole::x},
             {"y", WVObservationCoordinateRole::y},
             {"z", WVObservationCoordinateRole::z}}})
      addOccurrenceValue(plan, declaration.first,
                         WVObservationScalarType::real64, {"sample"},
                         WVObservationValueLayout::flat, declaration.second);
    addOccurrenceValue(plan, "sample-identifier",
                       WVObservationScalarType::integer64, {"sample"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    addOccurrenceValue(plan, "sample-direction",
                       WVObservationScalarType::boolean8, {"sample"},
                       WVObservationValueLayout::flat);
    addOccurrenceValue(plan, "payload-real",
                       WVObservationScalarType::real64, {},
                       WVObservationValueLayout::record);
    addOccurrenceValue(plan, "payload-integer",
                       WVObservationScalarType::integer64, {},
                       WVObservationValueLayout::record);
    addOccurrenceValue(plan, "payload-boolean",
                       WVObservationScalarType::boolean8, {},
                       WVObservationValueLayout::record);
    addOccurrencePositionSet(plan, "event-line");
    addOccurrenceField(plan, "sampled-u", "u", {"sample"}, 0);
    addOccurrenceField(plan, "sampled-ssh", "ssh", {"sample"}, 0);
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  variableSamplesByFixedDepthBins) {
    auto status = addPayloadSchema(plan);
    if (!status)
      return status;
    if (observer.z.size() != 3)
      return invalid("Depth-bin proof requires three configured depths.");
    plan.schema.axes.push_back(axis(
        "depth", WVObservationAxisKind::fixed, observer.z.size(),
        WVObservationCoordinateRole::depth));
    plan.schema.axes.push_back(axis(
        "sample", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::identifier));
    for (const auto &declaration :
         std::array<std::pair<const char *, WVObservationCoordinateRole>, 4>{{
             {"sample-time", WVObservationCoordinateRole::sampleTime},
             {"x", WVObservationCoordinateRole::x},
             {"y", WVObservationCoordinateRole::y},
             {"z", WVObservationCoordinateRole::z}}})
      addOccurrenceValue(plan, declaration.first,
                         WVObservationScalarType::real64,
                         {"depth", "sample"},
                         WVObservationValueLayout::flat, declaration.second);
    addOccurrenceValue(plan, "depth-identifier",
                       WVObservationScalarType::integer64,
                       {"depth", "sample"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    addOccurrencePositionSet(plan, "depth-bin-samples");
    addOccurrenceField(plan, "sampled-rho-e", "rho_e",
                       {"depth", "sample"}, 0);
    addOccurrenceField(plan, "sampled-ssu", "ssu",
                       {"depth", "sample"}, 0,
                       WVPositionInterpolation::spline);
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  stateCoupledIrregularGeometry) {
    auto status = addPayloadSchema(plan);
    if (!status)
      return status;
    if (observer.stateBlockIdentifiers.size() != 3)
      return invalid(
          "State-coupled proof requires x, y, and z coordinate state blocks.");
    for (const auto &identifier : observer.stateBlockIdentifiers)
      plan.occurrenceStateBlocks.push_back({identifier});
    plan.schema.axes.push_back(axis(
        "sample", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::identifier));
    for (const auto &declaration :
         std::array<std::pair<const char *, WVObservationCoordinateRole>, 4>{{
             {"sample-time", WVObservationCoordinateRole::sampleTime},
             {"x", WVObservationCoordinateRole::x},
             {"y", WVObservationCoordinateRole::y},
             {"z", WVObservationCoordinateRole::z}}})
      addOccurrenceValue(plan, declaration.first,
                         WVObservationScalarType::real64, {"sample"},
                         WVObservationValueLayout::flat, declaration.second);
    addOccurrenceValue(plan, "sample-identifier",
                       WVObservationScalarType::integer64, {"sample"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    addOccurrencePositionSet(plan, "state-coupled-samples");
    addOccurrenceField(plan, "sampled-v", "v", {"sample"}, 0);
    addOccurrenceField(plan, "sampled-ssv", "ssv", {"sample"}, 0);
  } else {
    auto status = addPayloadSchema(plan);
    if (!status)
      return status;
    plan.schema.axes.push_back(axis(
        "pass", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::pass));
    plan.schema.axes.push_back(axis(
        "profile", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::profile));
    plan.schema.axes.push_back(axis(
        "sample", WVObservationAxisKind::unlimited, 0,
        WVObservationCoordinateRole::identifier));
    addOccurrenceValue(plan, "profile-count-by-pass",
                       WVObservationScalarType::integer64, {"pass"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::none,
                       WVObservationRaggedRole::rowCount, "profile");
    addOccurrenceValue(plan, "sample-count-by-profile",
                       WVObservationScalarType::integer64, {"profile"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::none,
                       WVObservationRaggedRole::rowCount, "sample");
    addOccurrenceValue(plan, "pass-identifier",
                       WVObservationScalarType::integer64, {"pass"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    addOccurrenceValue(plan, "profile-identifier",
                       WVObservationScalarType::integer64, {"profile"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    addOccurrenceValue(plan, "sample-identifier",
                       WVObservationScalarType::integer64, {"sample"},
                       WVObservationValueLayout::flat,
                       WVObservationCoordinateRole::identifier);
    for (const auto &declaration :
         std::array<std::pair<const char *, WVObservationCoordinateRole>, 4>{{
             {"sample-time", WVObservationCoordinateRole::sampleTime},
             {"x", WVObservationCoordinateRole::x},
             {"y", WVObservationCoordinateRole::y},
             {"z", WVObservationCoordinateRole::z}}})
      addOccurrenceValue(plan, declaration.first,
                         WVObservationScalarType::real64, {"sample"},
                         WVObservationValueLayout::flat, declaration.second);
    addOccurrenceValue(plan, "sample-direction",
                       WVObservationScalarType::boolean8, {"sample"},
                       WVObservationValueLayout::flat);
    addOccurrenceValue(plan, "payload-real",
                       WVObservationScalarType::real64, {},
                       WVObservationValueLayout::record);
    addOccurrenceValue(plan, "payload-integer",
                       WVObservationScalarType::integer64, {},
                       WVObservationValueLayout::record);
    addOccurrenceValue(plan, "payload-boolean",
                       WVObservationScalarType::boolean8, {},
                       WVObservationValueLayout::record);
    addOccurrencePositionSet(plan, "nested-ragged-samples");
    addOccurrenceField(plan, "sampled-u", "u", {"sample"}, 0);
    addOccurrenceField(plan, "sampled-ssh", "ssh", {"sample"}, 0,
                       WVPositionInterpolation::spline);
  }

  const auto status = validateObservationSchema(plan.schema);
  if (!status)
    return status;
  output = std::move(plan);
  return WVKernelStatus::ok();
}

template <class Value>
inline void copyValues(Value *destination, const std::vector<Value> &source) {
  std::copy(source.begin(), source.end(), destination);
}

inline WVKernelStatus storeReal(WVObserverOccurrenceWorkspace &workspace,
                                std::size_t slot,
                                std::vector<std::size_t> extents,
                                const std::vector<double> &values) {
  double *data = nullptr;
  auto status = workspace.resizeReal(slot, std::move(extents), data);
  if (status)
    copyValues(data, values);
  return status;
}

inline WVKernelStatus storeInteger(
    WVObserverOccurrenceWorkspace &workspace, std::size_t slot,
    std::vector<std::size_t> extents,
    const std::vector<std::int64_t> &values) {
  std::int64_t *data = nullptr;
  auto status = workspace.resizeInteger(slot, std::move(extents), data);
  if (status)
    copyValues(data, values);
  return status;
}

inline WVKernelStatus storeBoolean(
    WVObserverOccurrenceWorkspace &workspace, std::size_t slot,
    std::vector<std::size_t> extents,
    const std::vector<std::uint8_t> &values) {
  std::uint8_t *data = nullptr;
  auto status = workspace.resizeBoolean(slot, std::move(extents), data);
  if (status)
    copyValues(data, values);
  return status;
}

inline WVKernelStatus storeScalarPayload(
    WVObserverOccurrenceWorkspace &workspace, std::size_t firstSlot,
    const PayloadScalars &payload) {
  auto status = storeReal(workspace, firstSlot, {}, {payload.real});
  if (status)
    status =
        storeInteger(workspace, firstSlot + 1, {}, {payload.integer});
  if (status)
    status =
        storeBoolean(workspace, firstSlot + 2, {}, {payload.boolean});
  return status;
}

inline WVKernelStatus storeCoordinateValues(
    WVObserverOccurrenceWorkspace &workspace, std::size_t firstSlot,
    const std::vector<std::size_t> &extents,
    const std::vector<double> &sampleTimes, const std::vector<double> &x,
    const std::vector<double> &y, const std::vector<double> &z) {
  auto status = storeReal(workspace, firstSlot, extents, sampleTimes);
  if (status)
    status = storeReal(workspace, firstSlot + 1, extents, x);
  if (status)
    status = storeReal(workspace, firstSlot + 2, extents, y);
  if (status)
    status = storeReal(workspace, firstSlot + 3, extents, z);
  return status;
}

inline void setPositionSet(WVObserverOccurrenceWorkspace &workspace,
                           std::size_t slot,
                           std::vector<std::size_t> extents,
                           std::vector<double> sampleTimes,
                           std::vector<double> x, std::vector<double> y,
                           std::vector<double> z) {
  auto &set = workspace.positionSets[slot];
  set.extents = std::move(extents);
  set.sampleTimes = std::move(sampleTimes);
  set.x = std::move(x);
  set.y = std::move(y);
  set.z = std::move(z);
}

template <WVTestObservationOccurrenceProviderKind Kind>
WVKernelStatus prepareWorkspace(
    const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
    const WVObserverOccurrencePreparationContext &context,
    WVObserverOccurrenceWorkspace &workspace) {
  workspace.prepareFor(plan);

  if constexpr (Kind ==
                WVTestObservationOccurrenceProviderKind::fixedGeometry) {
    const std::vector<std::size_t> extents{observer.x.size()};
    const std::vector<double> sampleTimes(observer.x.size(),
                                          context.scheduledTime);
    auto status = storeCoordinateValues(workspace, 0, extents, sampleTimes,
                                        observer.x, observer.y, observer.z);
    if (!status)
      return status;
    setPositionSet(workspace, 0, extents, sampleTimes, observer.x, observer.y,
                   observer.z);
    return WVKernelStatus::ok();
  }

  PayloadScalars payload;
  auto status = payloadScalars(context, payload);
  if (!status)
    return status;
  const auto direction = payload.boolean == 0 ? -1.0 : 1.0;
  const auto baseX = observer.x.empty() ? 1000.0 : observer.x.front();
  const auto baseY = observer.y.empty() ? 1200.0 : observer.y.front();
  const auto baseZ = observer.z.empty() ? -300.0 : observer.z.front();

  if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  eventVariableOneDimensionalGeometry) {
    const std::size_t count =
        payload.boolean == 0 ? 0 : 1 + magnitude(payload.integer) % 5;
    const std::vector<std::size_t> extents{count};
    std::vector<double> sampleTimes(count), x(count), y(count), z(count);
    std::vector<std::int64_t> identifiers(count);
    std::vector<std::uint8_t> directions(count, payload.boolean);
    for (std::size_t index = 0; index < count; ++index) {
      sampleTimes[index] = context.scheduledTime +
                           0.01 * (static_cast<double>(index) -
                                   0.5 * static_cast<double>(count));
      x[index] = baseX + direction *
                             (40.0 * static_cast<double>(index) +
                              0.25 * payload.real);
      y[index] = baseY + 25.0 * static_cast<double>(index) -
                 0.125 * payload.real;
      z[index] = baseZ - 15.0 * static_cast<double>(index);
      identifiers[index] = static_cast<std::int64_t>(index + 1);
    }
    status = storeCoordinateValues(workspace, 0, extents, sampleTimes, x, y,
                                   z);
    if (status)
      status = storeInteger(workspace, 4, extents, identifiers);
    if (status)
      status = storeBoolean(workspace, 5, extents, directions);
    if (status)
      status = storeScalarPayload(workspace, 6, payload);
    if (!status)
      return status;
    setPositionSet(workspace, 0, extents, sampleTimes, x, y, z);
    return WVKernelStatus::ok();
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  variableSamplesByFixedDepthBins) {
    const std::size_t perDepth =
        payload.boolean == 0 ? 0 : 1 + magnitude(payload.integer) % 3;
    const std::size_t depthCount = observer.z.size();
    const std::size_t count = depthCount * perDepth;
    const std::vector<std::size_t> extents{depthCount, perDepth};
    std::vector<double> sampleTimes(count), x(count), y(count), z(count);
    std::vector<std::int64_t> depthIdentifiers(count);
    for (std::size_t sample = 0; sample < perDepth; ++sample) {
      for (std::size_t depth = 0; depth < depthCount; ++depth) {
        const auto index = depth + depthCount * sample;
        sampleTimes[index] = context.scheduledTime +
                             0.02 * static_cast<double>(sample);
        x[index] = baseX + direction *
                               (60.0 * static_cast<double>(sample) +
                                0.1 * payload.real);
        y[index] = baseY + 35.0 * static_cast<double>(sample) +
                   5.0 * static_cast<double>(depth);
        z[index] = observer.z[depth];
        depthIdentifiers[index] = static_cast<std::int64_t>(depth + 1);
      }
    }
    status = storeCoordinateValues(workspace, 0, extents, sampleTimes, x, y,
                                   z);
    if (status)
      status = storeInteger(workspace, 4, extents, depthIdentifiers);
    if (!status)
      return status;
    setPositionSet(workspace, 0, extents, sampleTimes, x, y, z);
    return WVKernelStatus::ok();
  } else if constexpr (
      Kind == WVTestObservationOccurrenceProviderKind::
                  stateCoupledIrregularGeometry) {
    if (context.observerStateBlocks == nullptr ||
        context.observerStateBlockCount != 3)
      return invalid("State-coupled occurrence has no observer state.");
    const auto &xBlock = context.observerStateBlocks[0];
    const auto &yBlock = context.observerStateBlocks[1];
    const auto &zBlock = context.observerStateBlocks[2];
    if (xBlock.layout == nullptr || yBlock.layout == nullptr ||
        zBlock.layout == nullptr || xBlock.realData == nullptr ||
        yBlock.realData == nullptr || zBlock.realData == nullptr ||
        xBlock.layout->scalarType != WVStateScalarType::real64 ||
        yBlock.layout->scalarType != WVStateScalarType::real64 ||
        zBlock.layout->scalarType != WVStateScalarType::real64 ||
        xBlock.layout->dimensions.size() != 1 ||
        yBlock.layout->dimensions != xBlock.layout->dimensions ||
        zBlock.layout->dimensions != xBlock.layout->dimensions)
      return invalid(
          "State-coupled occurrence coordinate state is malformed.");
    const auto available = xBlock.layout->dimensions[0];
    const std::size_t count =
        payload.boolean == 0 || available == 0
            ? 0
            : 1 + magnitude(payload.integer) % available;
    const std::vector<std::size_t> extents{count};
    std::vector<double> sampleTimes(count), x(count), y(count), z(count);
    std::vector<std::int64_t> identifiers(count);
    for (std::size_t index = 0; index < count; ++index) {
      sampleTimes[index] = context.scheduledTime;
      x[index] = xBlock.realData[index] + 0.05 * payload.real;
      y[index] = yBlock.realData[index] - 0.05 * payload.real;
      z[index] = zBlock.realData[index];
      identifiers[index] = static_cast<std::int64_t>(index + 1);
    }
    status = storeCoordinateValues(workspace, 0, extents, sampleTimes, x, y,
                                   z);
    if (status)
      status = storeInteger(workspace, 4, extents, identifiers);
    if (!status)
      return status;
    setPositionSet(workspace, 0, extents, sampleTimes, x, y, z);
    return WVKernelStatus::ok();
  } else {
    const std::size_t passCount =
        payload.boolean == 0 ? 0 : 1 + magnitude(payload.integer) % 2;
    std::vector<std::int64_t> profileCounts(passCount);
    std::size_t profileCount = 0;
    for (std::size_t pass = 0; pass < passCount; ++pass) {
      profileCounts[pass] = static_cast<std::int64_t>(
          1 + (pass + magnitude(payload.integer)) % 2);
      profileCount += static_cast<std::size_t>(profileCounts[pass]);
    }
    std::vector<std::int64_t> sampleCounts(profileCount);
    std::size_t sampleCount = 0;
    for (std::size_t profile = 0; profile < profileCount; ++profile) {
      sampleCounts[profile] = static_cast<std::int64_t>(
          1 + (profile + magnitude(payload.integer)) % 3);
      sampleCount += static_cast<std::size_t>(sampleCounts[profile]);
    }
    std::vector<std::int64_t> passIdentifiers(passCount),
        profileIdentifiers(profileCount), sampleIdentifiers(sampleCount);
    std::vector<double> sampleTimes(sampleCount), x(sampleCount),
        y(sampleCount), z(sampleCount);
    std::vector<std::uint8_t> directions(sampleCount, payload.boolean);
    for (std::size_t pass = 0; pass < passCount; ++pass)
      passIdentifiers[pass] = static_cast<std::int64_t>(pass + 1);
    std::size_t sample = 0;
    for (std::size_t profile = 0; profile < profileCount; ++profile) {
      profileIdentifiers[profile] = static_cast<std::int64_t>(profile + 1);
      for (std::size_t local = 0;
           local < static_cast<std::size_t>(sampleCounts[profile]); ++local) {
        sampleIdentifiers[sample] = static_cast<std::int64_t>(sample + 1);
        sampleTimes[sample] = context.scheduledTime +
                              0.03 * static_cast<double>(local);
        x[sample] = baseX + direction *
                                (90.0 * static_cast<double>(profile) +
                                 12.0 * static_cast<double>(local) +
                                 0.1 * payload.real);
        y[sample] = baseY + 55.0 * static_cast<double>(profile) +
                    8.0 * static_cast<double>(local);
        z[sample] = baseZ - 80.0 * static_cast<double>(profile) -
                    10.0 * static_cast<double>(local);
        ++sample;
      }
    }
    status = storeInteger(workspace, 0, {passCount}, profileCounts);
    if (status)
      status = storeInteger(workspace, 1, {profileCount}, sampleCounts);
    if (status)
      status = storeInteger(workspace, 2, {passCount}, passIdentifiers);
    if (status)
      status =
          storeInteger(workspace, 3, {profileCount}, profileIdentifiers);
    if (status)
      status =
          storeInteger(workspace, 4, {sampleCount}, sampleIdentifiers);
    if (status)
      status = storeCoordinateValues(workspace, 5, {sampleCount}, sampleTimes,
                                     x, y, z);
    if (status)
      status = storeBoolean(workspace, 9, {sampleCount}, directions);
    if (status)
      status = storeScalarPayload(workspace, 10, payload);
    if (!status)
      return status;
    setPositionSet(workspace, 0, {sampleCount}, sampleTimes, x, y, z);
    return WVKernelStatus::ok();
  }
}

template <WVTestObservationOccurrenceProviderKind Kind>
class ProviderBase : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept final {
    return observation_occurrence_provider_detail::typeIdentifier<Kind>();
  }
  std::uint32_t contractVersion() const noexcept final { return 1; }

  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &blocks,
      std::map<std::string, std::size_t> &ownerCounts) const final {
    if constexpr (
        Kind == WVTestObservationOccurrenceProviderKind::
                    stateCoupledIrregularGeometry) {
      if (observer.stateBlockIdentifiers.size() != 3 || observer.x.empty() ||
          observer.x.size() != observer.y.size() ||
          observer.x.size() != observer.z.size())
        return invalid(
            "State-coupled proof requires three coordinate vectors of equal "
            "nonzero length.");
      for (const auto &identifier : observer.stateBlockIdentifiers) {
        const auto found = blocks.find(identifier);
        if (found == blocks.end() || found->second == nullptr ||
            found->second->scalarType != WVStateScalarType::real64 ||
            found->second->ownership != WVStateOwnership::integratorOwned ||
            found->second->dimensions !=
                std::vector<std::size_t>({observer.x.size()}))
          return invalid("State-coupled proof state blocks must be "
                         "integrator-owned coordinate vectors.");
        ++ownerCounts.at(identifier);
      }
    } else if (!observer.stateBlockIdentifiers.empty()) {
      return invalid("Passive occurrence proof unexpectedly owns state.");
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus bindIntegration(
      const WVObserverRecord &observer,
      const WVObserverIntegrationBinder &binder) const final {
    if constexpr (
        Kind == WVTestObservationOccurrenceProviderKind::
                    stateCoupledIrregularGeometry) {
      if (!binder.advectedPositions)
        return invalid(
            "State-coupled proof cannot bind its integrated coordinates.");
      return binder.advectedPositions(observer);
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const final {
    plan = {};
    plan.persistedName = observer.name;
    return WVKernelStatus::ok();
  }

  WVKernelStatus outputPlan(
      const WVObserverRecord &observer,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const final {
    counters[kindIndex<Kind>()].outputPlanCount.fetch_add(
        1, std::memory_order_relaxed);
    return buildOutputPlan<Kind>(observer, context, plan);
  }

  WVKernelStatus prepareOccurrence(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOccurrencePreparationContext &context,
      WVObserverOccurrenceWorkspace &workspace) const final {
    counters[kindIndex<Kind>()].occurrencePreparationCount.fetch_add(
        1, std::memory_order_relaxed);
    try {
      return prepareWorkspace<Kind>(observer, plan, context, workspace);
    } catch (const std::bad_alloc &) {
      return allocationFailure();
    }
  }

  WVKernelStatus observationBatch(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &batch) const final {
    counters[kindIndex<Kind>()].batchCount.fetch_add(
        1, std::memory_order_relaxed);
    return WVObservingSystem::observationBatch(observer, plan, context, kind,
                                               batch);
  }

  std::size_t persistentBytes() const noexcept final { return sizeof(*this); }
};

template <WVTestObservationOccurrenceProviderKind Kind, class Provider>
WVKernelStatus addProvider(WVExtensionCatalogBuilder &builder) {
  return builder.addObserverFactory(
      {typeIdentifier<Kind>(), 1,
       [](const WVObserverRecord &, const WVPortableTypedRecord &,
          std::shared_ptr<const WVObservingSystem> &result) {
         try {
           result = std::make_shared<Provider>();
           counters[kindIndex<Kind>()].constructionCount.fetch_add(
               1, std::memory_order_relaxed);
           return WVKernelStatus::ok();
         } catch (const std::bad_alloc &) {
           return allocationFailure();
         }
       },
       {},
       [](const WVObserverRecord &observer,
          const WVObserverOutputPlanningContext &context,
          WVObserverOutputPlan &plan) {
         counters[kindIndex<Kind>()].outputPlanCount.fetch_add(
             1, std::memory_order_relaxed);
         return buildOutputPlan<Kind>(observer, context, plan);
       }});
}

inline WVObserverRecord record(std::string identifier, std::string name,
                               const char *type) {
  WVObserverRecord result;
  result.identifier = std::move(identifier);
  result.name = std::move(name);
  result.typeIdentifier = type;
  result.contractVersion = 1;
  result.configuration.schemaIdentifier =
      std::string("wv-test-") + type + "-configuration-v1";
  result.configuration.schemaVersion = 1;
  result.x = {1000.0};
  result.y = {1200.0};
  result.z = {-300.0};
  return result;
}

} // namespace observation_occurrence_provider_detail

inline void resetObservationOccurrenceProviderCounters() noexcept {
  for (auto &value : observation_occurrence_provider_detail::counters) {
    value.constructionCount.store(0, std::memory_order_relaxed);
    value.outputPlanCount.store(0, std::memory_order_relaxed);
    value.occurrencePreparationCount.store(0, std::memory_order_relaxed);
    value.batchCount.store(0, std::memory_order_relaxed);
    value.integrationStageCount.store(0, std::memory_order_relaxed);
  }
}

inline WVTestObservationOccurrenceProviderCounters
observationOccurrenceProviderCounters(
    WVTestObservationOccurrenceProviderKind kind) noexcept {
  const auto &value = observation_occurrence_provider_detail::counters
      [static_cast<std::size_t>(kind)];
  return {value.constructionCount.load(std::memory_order_relaxed),
          value.outputPlanCount.load(std::memory_order_relaxed),
          value.occurrencePreparationCount.load(std::memory_order_relaxed),
          value.batchCount.load(std::memory_order_relaxed),
          value.integrationStageCount.load(std::memory_order_relaxed)};
}

class WVTestFixedGeometryObservation final
    : public observation_occurrence_provider_detail::ProviderBase<
          WVTestObservationOccurrenceProviderKind::fixedGeometry> {};

class WVTestEventVariableOneDimensionalGeometry final
    : public observation_occurrence_provider_detail::ProviderBase<
          WVTestObservationOccurrenceProviderKind::
              eventVariableOneDimensionalGeometry> {};

class WVTestVariableSamplesByFixedDepthBins final
    : public observation_occurrence_provider_detail::ProviderBase<
          WVTestObservationOccurrenceProviderKind::
              variableSamplesByFixedDepthBins> {};

class WVTestStateCoupledIrregularGeometry final
    : public observation_occurrence_provider_detail::ProviderBase<
          WVTestObservationOccurrenceProviderKind::
              stateCoupledIrregularGeometry> {};

class WVTestNestedRaggedGeometry final
    : public observation_occurrence_provider_detail::ProviderBase<
          WVTestObservationOccurrenceProviderKind::nestedRaggedGeometry> {};

inline WVKernelStatus registerObservationOccurrenceProviders(
    WVExtensionCatalogBuilder &builder) {
  using namespace observation_occurrence_provider_detail;
  auto status =
      addProvider<WVTestObservationOccurrenceProviderKind::fixedGeometry,
                  WVTestFixedGeometryObservation>(builder);
  if (status)
    status = addProvider<
        WVTestObservationOccurrenceProviderKind::
            eventVariableOneDimensionalGeometry,
        WVTestEventVariableOneDimensionalGeometry>(builder);
  if (status)
    status = addProvider<
        WVTestObservationOccurrenceProviderKind::
            variableSamplesByFixedDepthBins,
        WVTestVariableSamplesByFixedDepthBins>(builder);
  if (status)
    status = addProvider<
        WVTestObservationOccurrenceProviderKind::stateCoupledIrregularGeometry,
        WVTestStateCoupledIrregularGeometry>(builder);
  if (status)
    status =
        addProvider<WVTestObservationOccurrenceProviderKind::
                        nestedRaggedGeometry,
                    WVTestNestedRaggedGeometry>(builder);
  return status;
}

inline WVObserverRecord fixedGeometryObservationRecord(
    std::string identifier = "fixed-geometry",
    std::string name = "fixed_geometry") {
  auto result = observation_occurrence_provider_detail::record(
      std::move(identifier), std::move(name), fixedGeometryObservationType);
  result.x = {900.0, 2400.0, 5100.0};
  result.y = {800.0, 3200.0, 4700.0};
  result.z = {-150.0, -450.0, -850.0};
  return result;
}

inline WVObserverRecord eventVariableOneDimensionalGeometryRecord(
    std::string identifier = "event-variable-line",
    std::string name = "event_variable_line") {
  return observation_occurrence_provider_detail::record(
      std::move(identifier), std::move(name),
      eventVariableOneDimensionalGeometryType);
}

inline WVObserverRecord variableSamplesByFixedDepthBinsRecord(
    std::string identifier = "variable-depth-bins",
    std::string name = "variable_depth_bins") {
  auto result = observation_occurrence_provider_detail::record(
      std::move(identifier), std::move(name),
      variableSamplesByFixedDepthBinsType);
  result.z = {-120.0, -460.0, -920.0};
  return result;
}

inline WVObserverRecord stateCoupledIrregularGeometryRecord(
    std::string identifier = "state-coupled-irregular",
    std::string name = "state_coupled_irregular") {
  auto result = observation_occurrence_provider_detail::record(
      std::move(identifier), std::move(name),
      stateCoupledIrregularGeometryType);
  result.x = {700.0, 1900.0, 3600.0, 6200.0};
  result.y = {900.0, 2200.0, 4100.0, 5200.0};
  result.z = {-100.0, -350.0, -700.0, -1050.0};
  result.isXYOnly = false;
  result.horizontalAbsoluteTolerance = 1e-8;
  result.verticalAbsoluteTolerance = 1e-8;
  result.stateBlockIdentifiers = {"event-coordinate-x", "event-coordinate-y",
                                  "event-coordinate-z"};
  return result;
}

inline WVObserverRecord nestedRaggedGeometryRecord(
    std::string identifier = "nested-ragged",
    std::string name = "nested_ragged") {
  return observation_occurrence_provider_detail::record(
      std::move(identifier), std::move(name), nestedRaggedGeometryType);
}

} // namespace wavevortex::runtime::test
