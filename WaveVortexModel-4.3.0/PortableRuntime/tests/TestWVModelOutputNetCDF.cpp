#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"
#include "WVLegacyObserverCompatibility.hpp"

#include "WVReferenceFFTEngine.hpp"
#include "WVTestQuadraticSchedule.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::detail;
using namespace wavevortex::runtime::test;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

WVKernelStatus makeTestOccurrenceIdentity(
    const void *preparationOwner, std::uint64_t preparationGeneration,
    double scheduledTime,
    const WVOutputRouteView &route, const WVOutputObserverView &observer,
    WVObservationOccurrenceIdentity &output) {
  if (route.schedulePayload == nullptr)
    return {WVKernelStatusCode::invalidConfiguration,
            "A test observation route has no schedule payload."};
  output = {};
  output.preparationOwner = preparationOwner;
  output.preparationGeneration = preparationGeneration;
  output.preparedOccurrenceSlot = route.semanticScheduleOrdinal;
  output.observerOrdinal = observer.observerOrdinal;
  output.semanticScheduleOrdinal = route.semanticScheduleOrdinal;
  output.scheduleOrdinal = route.scheduleOrdinal;
  output.scheduledTime = scheduledTime;
  output.scheduleCursorIdentity = route.scheduleCursorIdentity;
  output.payloadFingerprint = route.schedulePayload->valueFingerprint();
  return WVKernelStatus::ok();
}

const std::shared_ptr<const WVExtensionCatalog> &modelOutputCatalog();

class WVTestFieldsImplementation final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestFields";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &record,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.fieldListAttribute = "fieldNames";
    plan.outputFields = record.fieldNames;
    for (const auto *family : {"Ap", "Am", "A0"})
      if (std::find(record.fieldNames.begin(), record.fieldNames.end(),
                    family) != record.fieldNames.end())
        plan.coefficientRestartFamilies.emplace_back(family);
    return WVKernelStatus::ok();
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &record,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const override {
    if (context.configuration == nullptr)
      return {WVKernelStatusCode::invalidConfiguration,
              "Test field output requires a transform configuration."};
    WVObserverOutputPlan candidate;
    candidate.schema.identifier =
        "legacy-" + record.identifier + "-observation-v1";
    candidate.schema.preservesLegacyEncoding = true;
    candidate.schema.metadata.attributes = {
        {"AnnotatedClass", typeIdentifier()},
        {"portableIdentifier", record.identifier}};
    candidate.schema.metadata.stringListAttributes = {
        {"fieldNames", record.fieldNames}};
    for (const auto &field : record.fieldNames) {
      const auto *metadata = findPortableVariable(field);
      if (metadata == nullptr)
        return {WVKernelStatusCode::invalidConfiguration,
                "Test field output is unsupported."};
      std::vector<std::string> names;
      for (std::size_t dimension = 0; dimension < metadata->dimensionCount;
           ++dimension)
        names.emplace_back(metadata->dimensions[dimension]);
      std::vector<std::size_t> extents;
      WVObserverOutputChannel channel;
      channel.sourceIdentifier = field;
      if (metadata->kind == WVPortableVariableKind::coefficient) {
        WVTransformConstantStratificationDescriptor transform;
        auto status = WVTransformConstantStratificationDescriptor::create(
            *context.configuration, transform);
        if (!status)
          return status;
        extents = {context.configuration->Nj, transform.Nkl()};
        channel.source = WVObserverOutputChannelSource::coefficient;
        channel.coefficientFamily =
            metadata->identifier == WVPortableVariable::Ap
                ? 0
                : metadata->identifier == WVPortableVariable::Am ? 1 : 2;
      } else {
        channel.source = WVObserverOutputChannelSource::sampledField;
        if (metadata->naturalRank == WVPortableNaturalRank::vertical)
          extents = {context.configuration->Nz};
        else if (metadata->naturalRank == WVPortableNaturalRank::horizontal)
          extents = {context.configuration->Nx, context.configuration->Ny};
        else if (metadata->naturalRank == WVPortableNaturalRank::scalar)
          extents = {};
        else
          extents = {context.configuration->Nx, context.configuration->Ny,
                     context.configuration->Nz};
      }
      for (std::size_t index = 0; index < names.size(); ++index) {
        const auto found = std::find_if(
            candidate.schema.axes.begin(), candidate.schema.axes.end(),
            [&](const auto &axis) { return axis.identifier == names[index]; });
        if (found == candidate.schema.axes.end())
          candidate.schema.axes.push_back(
              {names[index], names[index], WVObservationAxisKind::fixed,
               extents[index], WVObservationCoordinateRole::none});
      }
      WVObservationVariable variable;
      variable.identifier = "derived-" + field;
      variable.name = field;
      variable.scalarType =
          metadata->kind == WVPortableVariableKind::coefficient
              ? WVObservationScalarType::complex64
              : WVObservationScalarType::real64;
      variable.dimensionIdentifiers = names;
      variable.layout =
          context.isDynamicsLinear &&
                  !metadata->isVariableWithLinearTimeStep
              ? WVObservationValueLayout::initialValue
              : WVObservationValueLayout::record;
      variable.units = metadata->units;
      variable.description = metadata->description;
      channel.variableIdentifier = variable.identifier;
      candidate.schema.variables.push_back(std::move(variable));
      candidate.channels.push_back(std::move(channel));
    }
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  }
  WVKernelStatus
  validate(const WVObserverRecord &record,
           const std::map<std::string, const WVStateBlockRecord *> &,
           std::map<std::string, std::size_t> &) const override {
    return record.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                                "Test fields cannot own state blocks."};
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
};

class WVTestPortablePointDiagnosticImplementation final
    : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestPortablePointDiagnostic";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &record,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.fieldListAttribute = "fieldNames";
    plan.persistedName = record.name;
    plan.outputFields = record.fieldNames;
    return WVKernelStatus::ok();
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &record,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const override {
    if (context.configuration == nullptr || record.fieldNames.size() != 1)
      return {WVKernelStatusCode::invalidConfiguration,
              "Point diagnostic output planning is invalid."};
    const auto *metadata = findPortableVariable(record.fieldNames.front());
    if (metadata == nullptr ||
        metadata->kind != WVPortableVariableKind::field)
      return {WVKernelStatusCode::invalidConfiguration,
              "Point diagnostic field is unsupported."};
    WVObserverOutputPlan candidate;
    candidate.schema.identifier =
        "legacy-" + record.identifier + "-observation-v1";
    candidate.schema.preservesLegacyEncoding = true;
    candidate.schema.metadata.attributes = {
        {"AnnotatedClass", typeIdentifier()},
        {"portableIdentifier", record.identifier}, {"name", record.name},
        {"trackedVarInterpolation", "linear"}};
    candidate.schema.metadata.stringListAttributes = {
        {"fieldNames", record.fieldNames}};
    candidate.schema.metadata.variables.push_back(
        {"outputScale",
         WVObservationValue::ownReal("outputScale", {}, {record.outputScale}),
         false});
    candidate.schema.metadata.variables.push_back(
        {"outputOffset",
         WVObservationValue::ownReal("outputOffset", {}, {record.outputOffset}),
         false});
    const std::string idName = record.name + "_id";
    candidate.schema.axes.push_back(
        {idName, idName, WVObservationAxisKind::fixed, record.x.size(),
         WVObservationCoordinateRole::identifier});
    const auto addConstant = [&](std::string identifier, std::string name,
                                 const std::vector<double> &values,
                                 WVObservationCoordinateRole role,
                                 std::string description) {
      candidate.schema.variables.push_back(
          {identifier, std::move(name), WVObservationScalarType::real64,
           {idName}, WVObservationValueLayout::staticValue, "m",
           std::move(description), {}, role,
           WVObservationRaggedRole::none, {}});
      candidate.constantValues.push_back(WVObservationValue::ownReal(
          std::move(identifier), {values.size()}, values));
    };
    std::vector<double> identifiers(record.x.size());
    for (std::size_t index = 0; index < identifiers.size(); ++index)
      identifiers[index] = static_cast<double>(index + 1);
    addConstant("static-" + idName, idName, identifiers,
                WVObservationCoordinateRole::identifier, "");
    addConstant("static-x", record.name + "_x", record.x,
                WVObservationCoordinateRole::x,
                "x position of fixed observation");
    addConstant("static-y", record.name + "_y", record.y,
                WVObservationCoordinateRole::y,
                "y position of fixed observation");
    addConstant("static-z", record.name + "_z", record.z,
                WVObservationCoordinateRole::z,
                "z position of fixed observation");
    WVObservationVariable value;
    value.identifier = "derived-" + record.fieldNames.front();
    value.name = record.name + "_value";
    value.dimensionIdentifiers = {idName};
    value.layout =
        context.isDynamicsLinear && !metadata->isVariableWithLinearTimeStep
            ? WVObservationValueLayout::initialValue
            : WVObservationValueLayout::record;
    value.units = metadata->units;
    value.description = std::string(metadata->description) +
                        ", sampled and affinely transformed by the observing system";
    candidate.schema.variables.push_back(value);
    WVObserverOutputChannel channel;
    channel.variableIdentifier = value.identifier;
    channel.source = WVObserverOutputChannelSource::sampledField;
    channel.sourceIdentifier = record.fieldNames.front();
    channel.sampling.kind = WVFieldSamplingKind::positions;
    channel.sampling.x = record.x;
    channel.sampling.y = record.y;
    channel.sampling.z = record.z;
    channel.scale = record.outputScale;
    channel.offset = record.outputOffset;
    candidate.channels.push_back(std::move(channel));
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  }
  WVKernelStatus
  validate(const WVObserverRecord &record,
           const std::map<std::string, const WVStateBlockRecord *> &,
           std::map<std::string, std::size_t> &) const override {
    if (!record.stateBlockIdentifiers.empty() ||
        record.fieldNames.size() != 1 || record.x.empty() ||
        record.x.size() != record.y.size() ||
        record.x.size() != record.z.size() ||
        !std::isfinite(record.outputScale) ||
        !std::isfinite(record.outputOffset))
      return {WVKernelStatusCode::invalidConfiguration,
              "Point diagnostic configuration is invalid."};
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
};

class WVTestObservationBatchesImplementation final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestObservationBatches";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &record,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.persistedName = record.name;
    return WVKernelStatus::ok();
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &record,
      const WVObserverOutputPlanningContext &,
      WVObserverOutputPlan &output) const override {
    WVObserverOutputPlan plan;
    plan.schema.identifier = "synthetic-variable-observation";
    plan.schema.metadata.attributes = {
        {"AnnotatedClass", record.typeIdentifier},
        {"portableIdentifier", record.identifier}, {"name", record.name}};
    plan.schema.axes = {
        {"depth", "synthetic_depth", WVObservationAxisKind::fixed, 2,
         WVObservationCoordinateRole::depth},
        {"profile", "synthetic_profile", WVObservationAxisKind::unlimited, 0,
         WVObservationCoordinateRole::profile},
        {"sample", "synthetic_sample", WVObservationAxisKind::unlimited, 0,
         WVObservationCoordinateRole::identifier}};
    plan.schema.variables = {
        {"depth", "synthetic_depth", WVObservationScalarType::real64,
         {"depth"}, WVObservationValueLayout::staticValue, "m",
         "fixed depth-bin axis", {}, WVObservationCoordinateRole::depth,
         WVObservationRaggedRole::none, {}},
        {"profile-id", "synthetic_profile_id",
         WVObservationScalarType::integer64, {"profile"},
         WVObservationValueLayout::flat, "1", "profile identifier", {},
         WVObservationCoordinateRole::profile,
         WVObservationRaggedRole::none, {}},
        {"row-count", "synthetic_row_count",
         WVObservationScalarType::integer64, {"profile"},
         WVObservationValueLayout::flat, "1", "samples in each profile", {},
         WVObservationCoordinateRole::none,
         WVObservationRaggedRole::rowCount, "sample"},
        {"time", "synthetic_time", WVObservationScalarType::real64,
         {"sample"}, WVObservationValueLayout::flat, "s", "sample time", {},
         WVObservationCoordinateRole::sampleTime,
         WVObservationRaggedRole::none, {}},
        {"x", "synthetic_x", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving x coordinate", {},
         WVObservationCoordinateRole::x, WVObservationRaggedRole::none, {}},
        {"y", "synthetic_y", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving y coordinate", {},
         WVObservationCoordinateRole::y, WVObservationRaggedRole::none, {}},
        {"z", "synthetic_z", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving z coordinate", {},
         WVObservationCoordinateRole::z, WVObservationRaggedRole::none, {}},
        {"bins", "synthetic_bins", WVObservationScalarType::complex64,
         {"depth", "sample"}, WVObservationValueLayout::flat, "m s-1",
         "fixed-depth-bin samples", {{"provider", "fixed-bin"}},
         WVObservationCoordinateRole::none, WVObservationRaggedRole::none, {}},
        {"valid", "synthetic_valid", WVObservationScalarType::boolean8,
         {"sample"}, WVObservationValueLayout::flat, "1", "sample validity",
         {}, WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}},
        {"label", "synthetic_label", WVObservationScalarType::text,
         {"profile"}, WVObservationValueLayout::flat, "", "profile label",
         {}, WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}},
        {"pass", "synthetic_pass", WVObservationScalarType::integer64, {},
         WVObservationValueLayout::record, "1", "pass identifier", {},
         WVObservationCoordinateRole::pass,
         WVObservationRaggedRole::none, {}}};
    plan.constantValues.push_back(WVObservationValue::ownReal(
        "depth", {2}, std::vector<double>{-10.0, -20.0}));
    output = std::move(plan);
    return WVKernelStatus::ok();
  }
  WVKernelStatus observationBatch(
      const WVObserverRecord &record, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &output) const override {
    if (kind == WVObservationBatchKind::initial)
      return WVObservingSystem::observationBatch(record, plan, context, kind,
                                                  output);
    ++batchCounts_[record.identifier];
    const bool empty = context.scheduleOrdinal() % 2 == 1;
    const std::size_t profileCount = empty ? 0 : 2;
    const std::size_t sampleCount = empty ? 0 : 3;
    WVObservationBatch batch;
    batch.schemaIdentifier = plan.schema.identifier;
    batch.schemaVersion = plan.schema.version;
    batch.values.push_back(WVObservationValue::ownInteger(
        "profile-id", {profileCount},
        empty ? std::vector<std::int64_t>{}
              : std::vector<std::int64_t>{10, 11}));
    batch.values.push_back(WVObservationValue::ownInteger(
        "row-count", {profileCount},
        empty ? std::vector<std::int64_t>{}
              : std::vector<std::int64_t>{1, 2}));
    batch.values.push_back(WVObservationValue::ownReal(
        "time", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{context.scheduledTime(),
                                    context.scheduledTime() + 0.1,
                                    context.scheduledTime() + 0.2}));
    batch.values.push_back(WVObservationValue::ownReal(
        "x", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{1.0, 2.0, 3.0}));
    batch.values.push_back(WVObservationValue::ownReal(
        "y", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{4.0, 5.0, 6.0}));
    batch.values.push_back(WVObservationValue::ownReal(
        "z", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{-1.0, -2.0, -3.0}));
    batch.values.push_back(WVObservationValue::ownComplex(
        "bins", {2, sampleCount},
        empty ? std::vector<WVComplex64>{}
              : std::vector<WVComplex64>{{1.0, -1.0}, {2.0, -2.0},
                                         {3.0, -3.0}, {4.0, -4.0},
                                         {5.0, -5.0}, {6.0, -6.0}}));
    batch.values.push_back(WVObservationValue::ownBoolean(
        "valid", {sampleCount},
        empty ? std::vector<std::uint8_t>{}
              : std::vector<std::uint8_t>{1, 0, 1}));
    batch.values.push_back(WVObservationValue::ownText(
        "label", {profileCount},
        empty ? std::vector<std::string>{}
              : std::vector<std::string>{"", "pass-b"}));
    batch.values.push_back(WVObservationValue::ownInteger(
        "pass", {}, {context.scheduleOrdinal()}));
    for (std::size_t index = 0; index < batch.values.size(); ++index)
      batch.values[index].resolvedVariableIndex = index + 1;
    output = std::move(batch);
    return WVKernelStatus::ok();
  }
  WVKernelStatus
  validate(const WVObserverRecord &record,
           const std::map<std::string, const WVStateBlockRecord *> &,
           std::map<std::string, std::size_t> &) const override {
    return record.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                                "Synthetic batches cannot own state blocks."};
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
  std::size_t batchCount(const std::string &identifier) const {
    const auto found = batchCounts_.find(identifier);
    return found == batchCounts_.end() ? 0 : found->second;
  }
  void resetCounts() const { batchCounts_.clear(); }

private:
  mutable std::map<std::string, std::size_t> batchCounts_;
};

std::shared_ptr<WVTestObservationBatchesImplementation> registeredBatches;

enum class TestTopology {
  fixedBin,
  movingTrack,
  variablePass,
  raggedProfile,
  nestedRagged
};

class WVTestTopologyImplementation final : public WVObservingSystem {
public:
  WVTestTopologyImplementation(std::string identifier, TestTopology topology)
      : identifier_(std::move(identifier)), topology_(topology) {}
  const std::string &typeIdentifier() const noexcept override {
    return identifier_;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    return WVKernelStatus::ok();
  }
  WVKernelStatus validate(
      const WVObserverRecord &record,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    return record.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                                "Test topology providers are sample-only."};
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &record,
      const WVObserverOutputPlanningContext &,
      WVObserverOutputPlan &output) const override {
    WVObserverOutputPlan plan;
    plan.schema.identifier = schemaIdentifier();
    plan.schema.metadata.attributes = {
        {"AnnotatedClass", identifier_},
        {"portableIdentifier", record.identifier}, {"name", record.name}};
    if (topology_ == TestTopology::fixedBin) {
      plan.schema.axes = {
          {"depth", "fixed_bin_depth", WVObservationAxisKind::fixed, 3,
           WVObservationCoordinateRole::depth}};
      plan.schema.variables = {
          {"depth", "fixed_bin_depth", WVObservationScalarType::real64,
           {"depth"}, WVObservationValueLayout::staticValue, "m",
           "fixed depth bins", {}, WVObservationCoordinateRole::depth,
           WVObservationRaggedRole::none, {}},
          {"value", "fixed_bin_value", WVObservationScalarType::complex64,
           {"depth"}, WVObservationValueLayout::record, "m s-1",
           "fixed-bin values", {{"ordered-a", "1"}, {"ordered-b", "2"}},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::none, {}}};
      plan.constantValues.push_back(WVObservationValue::ownReal(
          "depth", {3}, {-5.0, -15.0, -25.0}));
    } else if (topology_ == TestTopology::movingTrack) {
      plan.schema.axes = {
          {"sample", "moving_track_sample",
           WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::identifier}};
      plan.schema.variables = {
          {"time", "moving_track_time", WVObservationScalarType::real64,
           {"sample"}, WVObservationValueLayout::flat, "s", "sample time",
           {}, WVObservationCoordinateRole::sampleTime,
           WVObservationRaggedRole::none, {}},
          {"x", "moving_track_x", WVObservationScalarType::real64,
           {"sample"}, WVObservationValueLayout::flat, "m", "moving x", {},
           WVObservationCoordinateRole::x, WVObservationRaggedRole::none, {}},
          {"y", "moving_track_y", WVObservationScalarType::real64,
           {"sample"}, WVObservationValueLayout::flat, "m", "moving y", {},
           WVObservationCoordinateRole::y, WVObservationRaggedRole::none, {}},
          {"z", "moving_track_z", WVObservationScalarType::real64,
           {"sample"}, WVObservationValueLayout::flat, "m", "moving z", {},
           WVObservationCoordinateRole::z, WVObservationRaggedRole::none, {}}};
    } else if (topology_ == TestTopology::variablePass) {
      plan.schema.axes = {
          {"pass", "variable_pass", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::pass}};
      plan.schema.variables = {
          {"pass", "variable_pass_id", WVObservationScalarType::integer64,
           {"pass"}, WVObservationValueLayout::flat, "1", "pass id", {},
           WVObservationCoordinateRole::pass,
           WVObservationRaggedRole::none, {}},
          {"label", "variable_pass_label", WVObservationScalarType::text,
           {"pass"}, WVObservationValueLayout::flat, "", "pass label", {},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::none, {}}};
    } else if (topology_ == TestTopology::raggedProfile) {
      plan.schema.axes = {
          {"profile", "ragged_profile", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::profile},
          {"sample", "ragged_sample", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::identifier}};
      plan.schema.variables = {
          {"profile", "ragged_profile_id",
           WVObservationScalarType::integer64, {"profile"},
           WVObservationValueLayout::flat, "1", "profile id", {},
           WVObservationCoordinateRole::profile,
           WVObservationRaggedRole::none, {}},
          {"row-count", "ragged_row_count",
           WVObservationScalarType::integer64, {"profile"},
           WVObservationValueLayout::flat, "1", "row count", {},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::rowCount, "sample"},
          {"value", "ragged_value", WVObservationScalarType::real64,
           {"sample"}, WVObservationValueLayout::flat, "1", "sample value",
           {}, WVObservationCoordinateRole::none,
           WVObservationRaggedRole::none, {}}};
    } else {
      plan.schema.axes = {
          {"pass", "nested_pass", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::pass},
          {"profile", "nested_profile", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::profile},
          {"depth", "nested_depth", WVObservationAxisKind::fixed, 2,
           WVObservationCoordinateRole::depth},
          {"sample", "nested_sample", WVObservationAxisKind::unlimited, 0,
           WVObservationCoordinateRole::identifier}};
      plan.schema.variables = {
          {"pass-profile-count", "nested_pass_profile_count",
           WVObservationScalarType::integer64, {"pass"},
           WVObservationValueLayout::flat, "1", "profiles in each pass", {},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::rowCount, "profile"},
          {"profile-sample-offset", "nested_profile_sample_offset",
           WVObservationScalarType::integer64, {"profile"},
           WVObservationValueLayout::flat, "1",
           "first sample in each profile", {},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::rowOffset, "sample"},
          {"sample-value", "nested_sample_value",
           WVObservationScalarType::real64, {"depth", "sample"},
           WVObservationValueLayout::flat, "1",
           "sample value by fixed depth", {},
           WVObservationCoordinateRole::none,
           WVObservationRaggedRole::none, {}}};
    }
    output = std::move(plan);
    return WVKernelStatus::ok();
  }
  WVKernelStatus observationBatch(
      const WVObserverRecord &record, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &output) const override {
    if (kind == WVObservationBatchKind::initial)
      return WVObservingSystem::observationBatch(record, plan, context, kind,
                                                  output);
    ++batchCounts_[record.identifier];
    WVObservationBatch batch;
    batch.schemaIdentifier = plan.schema.identifier;
    batch.schemaVersion = plan.schema.version;
    if (topology_ == TestTopology::fixedBin) {
      batch.values.push_back(WVObservationValue::ownComplex(
          "value", {3}, {{1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0}}));
    } else if (topology_ == TestTopology::movingTrack) {
      const std::size_t count = context.scheduleOrdinal() % 2 == 0 ? 2 : 1;
      batch.values.push_back(WVObservationValue::ownReal(
          "time", {count},
          count == 2 ? std::vector<double>{context.scheduledTime(),
                                           context.scheduledTime() + 0.25}
                     : std::vector<double>{context.scheduledTime()}));
      batch.values.push_back(WVObservationValue::ownReal(
          "x", {count}, count == 2 ? std::vector<double>{1.0, 2.0}
                                     : std::vector<double>{3.0}));
      batch.values.push_back(WVObservationValue::ownReal(
          "y", {count}, count == 2 ? std::vector<double>{4.0, 5.0}
                                     : std::vector<double>{6.0}));
      batch.values.push_back(WVObservationValue::ownReal(
          "z", {count}, count == 2 ? std::vector<double>{-1.0, -2.0}
                                     : std::vector<double>{-3.0}));
    } else if (topology_ == TestTopology::variablePass) {
      const std::size_t count = context.scheduleOrdinal() % 2 == 0 ? 2 : 1;
      batch.values.push_back(WVObservationValue::ownInteger(
          "pass", {count}, count == 2 ? std::vector<std::int64_t>{7, 8}
                                       : std::vector<std::int64_t>{9}));
      batch.values.push_back(WVObservationValue::ownText(
          "label", {count}, count == 2
                                ? std::vector<std::string>{"out", "back"}
                                : std::vector<std::string>{"final"}));
    } else if (topology_ == TestTopology::raggedProfile) {
      const bool malformed = record.name == "malformed";
      batch.values.push_back(WVObservationValue::ownInteger(
          "profile", {2}, {1, 2}));
      batch.values.push_back(WVObservationValue::ownInteger(
          "row-count", {2}, malformed ? std::vector<std::int64_t>{2, 2}
                                        : std::vector<std::int64_t>{1, 2}));
      batch.values.push_back(WVObservationValue::ownReal(
          "value", {3}, {10.0, 20.0, 30.0}));
    } else {
      batch.values.push_back(WVObservationValue::ownInteger(
          "pass-profile-count", {2}, {2, 1}));
      batch.values.push_back(WVObservationValue::ownInteger(
          "profile-sample-offset", {3}, {0, 2, 2}));
      batch.values.push_back(WVObservationValue::ownReal(
          "sample-value", {2, 4},
          {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}));
    }
    for (std::size_t index = 0; index < batch.values.size(); ++index)
      batch.values[index].resolvedVariableIndex =
          index + (topology_ == TestTopology::fixedBin ? 1 : 0);
    output = std::move(batch);
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + identifier_.capacity();
  }
  std::size_t batchCount(const std::string &identifier) const {
    const auto found = batchCounts_.find(identifier);
    return found == batchCounts_.end() ? 0 : found->second;
  }
  void resetCounts() const { batchCounts_.clear(); }

private:
  std::string schemaIdentifier() const {
    if (topology_ == TestTopology::fixedBin)
      return "test-fixed-bin-v1";
    if (topology_ == TestTopology::movingTrack)
      return "test-moving-track-v1";
    if (topology_ == TestTopology::variablePass)
      return "test-variable-pass-v1";
    if (topology_ == TestTopology::raggedProfile)
      return "test-ragged-profile-v1";
    return "test-nested-ragged-v1";
  }
  std::string identifier_;
  TestTopology topology_;
  mutable std::map<std::string, std::size_t> batchCounts_;
};

std::shared_ptr<WVTestTopologyImplementation> fixedBinProvider;
std::shared_ptr<WVTestTopologyImplementation> movingTrackProvider;
std::shared_ptr<WVTestTopologyImplementation> variablePassProvider;
std::shared_ptr<WVTestTopologyImplementation> raggedProfileProvider;
std::shared_ptr<WVTestTopologyImplementation> nestedRaggedProvider;

const std::shared_ptr<const WVExtensionCatalog> &modelOutputCatalog() {
  static const auto catalog = [] {
    WVExtensionCatalogBuilder builder;
    auto status = addBuiltInExtensions(builder);
    const auto add = [&](std::string identifier, WVObserverFactory factory,
                         WVLegacyObserverOperationResolver legacy = {},
                         WVLegacyObserverPersistenceMetadata persistence = {},
                         WVObserverOutputPlanResolver outputPlan = {}) {
      if (!status)
        return;
      WVObserverFactoryRegistration registration(
          std::move(identifier), 1, std::move(factory), {},
          std::move(outputPlan));
      if (legacy)
        status =
            WVObserverFactoryRegistrationAccess::attachLegacyCompatibility(
                registration, std::move(legacy), std::move(persistence));
      if (status)
        status = builder.addObserverFactory(std::move(registration));
    };
    add("WVTestFields",
        [](const WVObserverRecord &, const WVPortableTypedRecord &,
           std::shared_ptr<const WVObservingSystem> &result) {
          result = std::make_shared<WVTestFieldsImplementation>();
          return WVKernelStatus::ok();
        },
        [](const WVObserverRecord &,
           const WVLegacyObserverOperationBinder &binder) {
          return binder.fullField();
        },
        {"fieldNames", {}, {}, false},
        [](const WVObserverRecord &record,
           const WVObserverOutputPlanningContext &context,
           WVObserverOutputPlan &plan) {
          return WVTestFieldsImplementation().outputPlan(record, context,
                                                         plan);
        });
    add("WVTestPortablePointDiagnostic",
        [](const WVObserverRecord &, const WVPortableTypedRecord &,
           std::shared_ptr<const WVObservingSystem> &result) {
          result =
              std::make_shared<WVTestPortablePointDiagnosticImplementation>();
          return WVKernelStatus::ok();
        },
        [](const WVObserverRecord &,
           const WVLegacyObserverOperationBinder &binder) {
          return binder.fixedPositions();
        },
        {"fieldNames", {}, {}, false},
        [](const WVObserverRecord &record,
           const WVObserverOutputPlanningContext &context,
           WVObserverOutputPlan &plan) {
          return WVTestPortablePointDiagnosticImplementation().outputPlan(
              record, context, plan);
        });
    add("WVTestObservationBatches",
        [](const WVObserverRecord &, const WVPortableTypedRecord &,
           std::shared_ptr<const WVObservingSystem> &result) {
          registeredBatches =
              std::make_shared<WVTestObservationBatchesImplementation>();
          result = registeredBatches;
          return WVKernelStatus::ok();
        }, {}, {},
        [](const WVObserverRecord &record,
           const WVObserverOutputPlanningContext &context,
           WVObserverOutputPlan &plan) {
          return WVTestObservationBatchesImplementation().outputPlan(
              record, context, plan);
        });
    const auto addTopology = [&](const char *identifier,
                                 TestTopology topology,
                                 std::shared_ptr<WVTestTopologyImplementation>
                                     *lastCreated) {
      add(identifier,
          [identifier, topology, lastCreated](
              const WVObserverRecord &, const WVPortableTypedRecord &,
              std::shared_ptr<const WVObservingSystem> &result) {
            *lastCreated = std::make_shared<WVTestTopologyImplementation>(
                identifier, topology);
            result = *lastCreated;
            return WVKernelStatus::ok();
          }, {}, {},
          [identifier, topology](
              const WVObserverRecord &record,
              const WVObserverOutputPlanningContext &context,
              WVObserverOutputPlan &plan) {
            return WVTestTopologyImplementation(identifier, topology)
                .outputPlan(record, context, plan);
          });
    };
    addTopology("WVTestFixedBin", TestTopology::fixedBin, &fixedBinProvider);
    addTopology("WVTestMovingTrack", TestTopology::movingTrack,
                &movingTrackProvider);
    addTopology("WVTestVariablePass", TestTopology::variablePass,
                &variablePassProvider);
    addTopology("WVTestRaggedProfile", TestTopology::raggedProfile,
                &raggedProfileProvider);
    addTopology("WVTestNestedRagged", TestTopology::nestedRagged,
                &nestedRaggedProvider);
    if (status)
      status = registerQuadraticSchedule(builder);
    std::shared_ptr<const WVExtensionCatalog> result;
    if (status)
      status = builder.freeze(result);
    if (!status)
      throw std::runtime_error(status.message);
    return result;
  }();
  return catalog;
}

struct InspectionFactoryCounts {
  std::size_t observers = 0;
  std::size_t schedules = 0;
  std::size_t forcings = 0;
};

InspectionFactoryCounts *activeInspectionFactoryCounts = nullptr;

std::shared_ptr<const WVOutputSchedule> inspectionTrapSchedule(
    const WVOutputScheduleRecord &, WVKernelStatus &status) {
  ++activeInspectionFactoryCounts->schedules;
  status = {WVKernelStatusCode::invalidConfiguration,
            "raw inspection constructed an output schedule"};
  return {};
}

std::shared_ptr<const WVExtensionCatalog>
inspectionTrapCatalog(InspectionFactoryCounts &counts) {
  activeInspectionFactoryCounts = &counts;
  WVExtensionCatalogBuilder builder;
  for (auto registration : modelOutputCatalog()->observers().registrations()) {
    registration.factory =
        [&counts](const WVObserverRecord &, const WVPortableTypedRecord &,
                  std::shared_ptr<const WVObservingSystem> &) {
      ++counts.observers;
      return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                            "raw inspection constructed an observer"};
    };
    const auto status = builder.addObserverFactory(std::move(registration));
    if (!status)
      throw std::runtime_error(status.message);
  }
  for (auto registration :
       modelOutputCatalog()->outputSchedules().registrations()) {
    registration.factory = &inspectionTrapSchedule;
    const auto status =
        builder.addOutputScheduleFactory(std::move(registration));
    if (!status)
      throw std::runtime_error(status.message);
  }
  for (auto registration : modelOutputCatalog()->forcings().registrations()) {
    if (registration.isSupported)
      registration.factory =
          [&counts](const WVFrozenForcingEntry &,
                    const WVTransformConstantStratificationDescriptor &, bool,
                    std::unique_ptr<WVForcing> &) {
        ++counts.forcings;
        return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                              "raw inspection constructed a forcing"};
      };
    const auto status = builder.addForcingFactory(std::move(registration));
    if (!status)
      throw std::runtime_error(status.message);
  }
  std::shared_ptr<const WVExtensionCatalog> catalog;
  const auto status = builder.freeze(catalog);
  if (!status)
    throw std::runtime_error(status.message);
  return catalog;
}

class VariableBatchSource final : public WVObserverSampleSource {
public:
  explicit VariableBatchSource(std::uint32_t version = 1,
                               bool failFirstBatch = false)
      : version_(version), failFirstBatch_(failFirstBatch) {}

  WVKernelStatus observationSchema(const WVObserverRecord &observer,
                                   WVObservationSchema &output) override {
    WVObservationSchema schema;
    if (observer.typeIdentifier != "WVTestObservationBatches") {
      schema.identifier = "legacy-coefficient-observation-v1";
      schema.preservesLegacyEncoding = true;
      schema.metadata.attributes = {{"AnnotatedClass", observer.typeIdentifier},
                                    {"portableIdentifier", observer.identifier},
                                    {"name", observer.name}};
      WVObservationMetadataVariable tolerance;
      tolerance.name = "absTolerance";
      tolerance.value =
          WVObservationValue::ownReal("absTolerance", {}, {1e-6});
      schema.metadata.variables.push_back(std::move(tolerance));
      output = std::move(schema);
      return WVKernelStatus::ok();
    }
    schema.identifier = "synthetic-variable-observation";
    schema.version = version_;
    schema.metadata.attributes = {{"AnnotatedClass", observer.typeIdentifier},
                                  {"portableIdentifier", observer.identifier},
                                  {"name", observer.name}};
    schema.axes = {
        {"depth", "synthetic_depth", WVObservationAxisKind::fixed, 2,
         WVObservationCoordinateRole::depth},
        {"profile", "synthetic_profile", WVObservationAxisKind::unlimited, 0,
         WVObservationCoordinateRole::profile},
        {"sample", "synthetic_sample", WVObservationAxisKind::unlimited, 0,
         WVObservationCoordinateRole::identifier}};
    schema.variables = {
        {"depth", "synthetic_depth", WVObservationScalarType::real64,
         {"depth"}, WVObservationValueLayout::staticValue, "m",
         "fixed depth-bin axis", {}, WVObservationCoordinateRole::depth,
         WVObservationRaggedRole::none, {}},
        {"profile-id", "synthetic_profile_id",
         WVObservationScalarType::integer64, {"profile"},
         WVObservationValueLayout::flat, "1", "profile identifier", {},
         WVObservationCoordinateRole::profile,
         WVObservationRaggedRole::none, {}},
        {"row-count", "synthetic_row_count",
         WVObservationScalarType::integer64, {"profile"},
         WVObservationValueLayout::flat, "1", "samples in each profile", {},
         WVObservationCoordinateRole::none,
         WVObservationRaggedRole::rowCount, "sample"},
        {"time", "synthetic_time", WVObservationScalarType::real64,
         {"sample"}, WVObservationValueLayout::flat, "s", "sample time", {},
         WVObservationCoordinateRole::sampleTime,
         WVObservationRaggedRole::none, {}},
        {"x", "synthetic_x", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving x coordinate", {},
         WVObservationCoordinateRole::x, WVObservationRaggedRole::none, {}},
        {"y", "synthetic_y", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving y coordinate", {},
         WVObservationCoordinateRole::y, WVObservationRaggedRole::none, {}},
        {"z", "synthetic_z", WVObservationScalarType::real64, {"sample"},
         WVObservationValueLayout::flat, "m", "moving z coordinate", {},
         WVObservationCoordinateRole::z, WVObservationRaggedRole::none, {}},
        {"bins", "synthetic_bins", WVObservationScalarType::complex64,
         {"depth", "sample"}, WVObservationValueLayout::flat, "m s-1",
         "fixed-depth-bin samples", {}, WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}},
        {"valid", "synthetic_valid", WVObservationScalarType::boolean8,
         {"sample"}, WVObservationValueLayout::flat, "1", "sample validity",
         {}, WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}},
        {"label", "synthetic_label", WVObservationScalarType::text,
         {"profile"}, WVObservationValueLayout::flat, "", "profile label",
         {}, WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}},
        {"pass", "synthetic_pass", WVObservationScalarType::integer64, {},
         WVObservationValueLayout::record, "1", "pass identifier", {},
         WVObservationCoordinateRole::pass,
         WVObservationRaggedRole::none, {}}};
    output = std::move(schema);
    return WVKernelStatus::ok();
  }

  WVKernelStatus initialObservationBatch(
      const WVObserverRecord &observer, WVObservationBatch &output) override {
    WVObservationSchema schema;
    auto status = observationSchema(observer, schema);
    if (!status)
      return status;
    WVObservationBatch batch;
    batch.schemaIdentifier = schema.identifier;
    batch.schemaVersion = schema.version;
    batch.kind = WVObservationBatchKind::initial;
    if (observer.typeIdentifier == "WVTestObservationBatches")
      batch.values.push_back(WVObservationValue::ownReal(
          "depth", {2}, std::vector<double>{-10.0, -20.0}));
    output = std::move(batch);
    return WVKernelStatus::ok();
  }

  WVKernelStatus prepare(const WVOutputEvent &event) override {
    ++prepareCount_;
    ++preparationGeneration_;
    if (preparationGeneration_ == 0)
      ++preparationGeneration_;
    scheduledTime_ = event.scheduledTime;
    return WVKernelStatus::ok();
  }

  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const override {
    return makeTestOccurrenceIdentity(this, preparationGeneration_,
                                      scheduledTime_, route, observer, output);
  }

  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &identity,
      const WVObserverRecord &observer, WVObservationBatch &output) override {
    ++batchCounts_[observer.identifier];
    WVObservationSchema schema;
    auto status = observationSchema(observer, schema);
    if (!status)
      return status;
    WVObservationBatch batch;
    batch.schemaIdentifier = schema.identifier;
    batch.schemaVersion = schema.version;
    if (observer.typeIdentifier != "WVTestObservationBatches") {
      output = std::move(batch);
      return WVKernelStatus::ok();
    }
    const bool empty = identity.scheduleOrdinal % 2 == 1;
    const std::size_t profileCount = empty ? 0 : 2;
    const std::size_t sampleCount = empty ? 0 : 3;
    batch.values.push_back(WVObservationValue::ownInteger(
        "profile-id", {profileCount},
        empty ? std::vector<std::int64_t>{}
              : std::vector<std::int64_t>{10, 11}));
    batch.values.push_back(WVObservationValue::ownInteger(
        "row-count", {profileCount},
        empty ? std::vector<std::int64_t>{}
              : failFirstBatch_ && !failedOnce_
                    ? std::vector<std::int64_t>{2, 2}
                    : std::vector<std::int64_t>{1, 2}));
    batch.values.push_back(WVObservationValue::ownReal(
        "time", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{identity.scheduledTime,
                                    identity.scheduledTime + 0.1,
                                    identity.scheduledTime + 0.2}));
    batch.values.push_back(WVObservationValue::ownReal(
        "x", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{1.0, 2.0, 3.0}));
    batch.values.push_back(WVObservationValue::ownReal(
        "y", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{4.0, 5.0, 6.0}));
    batch.values.push_back(WVObservationValue::ownReal(
        "z", {sampleCount},
        empty ? std::vector<double>{}
              : std::vector<double>{-1.0, -2.0, -3.0}));
    batch.values.push_back(WVObservationValue::ownComplex(
        "bins", {2, sampleCount},
        empty ? std::vector<WVComplex64>{}
              : std::vector<WVComplex64>{{1.0, -1.0}, {2.0, -2.0},
                                         {3.0, -3.0}, {4.0, -4.0},
                                         {5.0, -5.0}, {6.0, -6.0}}));
    batch.values.push_back(WVObservationValue::ownBoolean(
        "valid", {sampleCount},
        empty ? std::vector<std::uint8_t>{}
              : std::vector<std::uint8_t>{1, 0, 1}));
    batch.values.push_back(WVObservationValue::ownText(
        "label", {profileCount},
        empty ? std::vector<std::string>{}
              : std::vector<std::string>{"", "pass-b"}));
    batch.values.push_back(WVObservationValue::ownInteger(
        "pass", {},
        std::vector<std::int64_t>{static_cast<std::int64_t>(
            identity.scheduleOrdinal)}));
    for (std::size_t index = 0; index < batch.values.size(); ++index)
      batch.values[index].resolvedVariableIndex = index + 1;
    if (failFirstBatch_ && !failedOnce_)
      failedOnce_ = true;
    output = std::move(batch);
    return WVKernelStatus::ok();
  }

  std::size_t prepareCount() const noexcept { return prepareCount_; }

  std::size_t batchCount(const std::string &observerIdentifier) const {
    const auto found = batchCounts_.find(observerIdentifier);
    return found == batchCounts_.end() ? 0 : found->second;
  }

private:
  std::uint32_t version_ = 1;
  bool failFirstBatch_ = false;
  bool failedOnce_ = false;
  double scheduledTime_ = 0.0;
  std::size_t prepareCount_ = 0;
  std::uint64_t preparationGeneration_ = 0;
  std::map<std::string, std::size_t> batchCounts_;
};

class ConflictingSchemaSource final : public WVObserverSampleSource {
public:
  explicit ConflictingSchemaSource(bool conflictingAxis)
      : conflictingAxis_(conflictingAxis) {}
  WVKernelStatus observationSchema(const WVObserverRecord &observer,
                                   WVObservationSchema &output) override {
    WVObservationSchema schema;
    if (observer.typeIdentifier != "WVTestObservationBatches") {
      schema.identifier = "legacy-coefficient-observation-v1";
      schema.preservesLegacyEncoding = true;
      output = std::move(schema);
      return WVKernelStatus::ok();
    }
    schema.identifier = "conflict-" + observer.identifier;
    const std::size_t extent =
        conflictingAxis_ && observer.identifier == "conflict-b" ? 3 : 2;
    schema.axes = {{"shared", "shared_axis", WVObservationAxisKind::fixed,
                    extent, WVObservationCoordinateRole::identifier}};
    schema.variables = {
        {"value",
         conflictingAxis_ ? observer.identifier + "_value" : "shared_value",
         WVObservationScalarType::real64, {"shared"},
         WVObservationValueLayout::record, "1", "conflict probe", {},
         WVObservationCoordinateRole::none,
         WVObservationRaggedRole::none, {}}};
    output = std::move(schema);
    return WVKernelStatus::ok();
  }
  WVKernelStatus initialObservationBatch(
      const WVObserverRecord &observer, WVObservationBatch &output) override {
    WVObservationSchema schema;
    auto status = observationSchema(observer, schema);
    if (!status)
      return status;
    WVObservationBatch batch;
    batch.schemaIdentifier = schema.identifier;
    batch.schemaVersion = schema.version;
    batch.kind = WVObservationBatchKind::initial;
    output = std::move(batch);
    return WVKernelStatus::ok();
  }
  WVKernelStatus prepare(const WVOutputEvent &event) override {
    ++preparationGeneration_;
    if (preparationGeneration_ == 0)
      ++preparationGeneration_;
    scheduledTime_ = event.scheduledTime;
    return WVKernelStatus::ok();
  }
  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const override {
    return makeTestOccurrenceIdentity(this, preparationGeneration_,
                                      scheduledTime_, route, observer, output);
  }
  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &,
      const WVObserverRecord &observer, WVObservationBatch &output) override {
    WVObservationSchema schema;
    auto status = observationSchema(observer, schema);
    if (!status)
      return status;
    WVObservationBatch batch;
    batch.schemaIdentifier = schema.identifier;
    batch.schemaVersion = schema.version;
    const auto extent = schema.axes.empty() ? 0 : schema.axes.front().extent;
    if (!schema.variables.empty())
      batch.values.push_back(WVObservationValue::ownReal(
          "value", {extent}, std::vector<double>(extent, 0.0)));
    if (!batch.values.empty())
      batch.values.front().resolvedVariableIndex = 0;
    output = std::move(batch);
    return WVKernelStatus::ok();
  }

private:
  bool conflictingAxis_ = false;
  double scheduledTime_ = 0.0;
  std::uint64_t preparationGeneration_ = 0;
};

std::filesystem::path fixture(const std::string &name) {
  return std::filesystem::path(WV_CHECKPOINT_FIXTURE_DIR) / name;
}

std::string shellQuote(const std::filesystem::path &path) {
  std::string result = "'";
  for (const auto character : path.string())
    result += character == '\'' ? "'\\''" : std::string(1, character);
  return result + "'";
}

std::vector<char> fileBytes(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

std::string fileText(const std::filesystem::path &path) {
  const auto bytes = fileBytes(path);
  return {bytes.begin(), bytes.end()};
}

double jsonNumber(const std::string &document, const std::string &key) {
  const auto position = document.find("\"" + key + "\":");
  require(position != std::string::npos,
          "JSON report field is missing: " + key);
  const char *begin = document.c_str() + position + key.size() + 3;
  char *end = nullptr;
  const auto value = std::strtod(begin, &end);
  require(end != begin, "JSON report field is not numeric: " + key);
  return value;
}

double coefficientDifference(const WVCheckpoint &first,
                             const WVCheckpoint &second) {
  double difference = 0.0;
  const std::vector<WVComplex64> *left[] = {
      &first.state.coefficients.Ap, &first.state.coefficients.Am,
      &first.state.coefficients.A0};
  const std::vector<WVComplex64> *right[] = {
      &second.state.coefficients.Ap, &second.state.coefficients.Am,
      &second.state.coefficients.A0};
  for (std::size_t family = 0; family < 3; ++family) {
    require(left[family]->size() == right[family]->size(),
            "checkpoint coefficient sizes differ");
    for (std::size_t index = 0; index < left[family]->size(); ++index)
      difference = std::max(
          difference,
          std::hypot((*left[family])[index].real -
                         (*right[family])[index].real,
                     (*left[family])[index].imag -
                         (*right[family])[index].imag));
  }
  return difference;
}

class FactoryPreflightRejectingSampleSource final
    : public WVObserverSampleSource {
public:
  WVKernelStatus preflight(const WVOutputPlan &) override {
    ++preflightCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "injected direct-factory source preflight rejection"};
  }
  WVKernelStatus observationSchema(const WVObserverRecord &,
                                   WVObservationSchema &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "schema discovery followed rejected source preflight"};
  }
  WVKernelStatus prepare(const WVOutputEvent &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "event preparation followed rejected source preflight"};
  }
  WVKernelStatus
  preparedOccurrenceIdentity(const WVOutputRouteView &,
                             const WVOutputObserverView &,
                             WVObservationOccurrenceIdentity &) const override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "occurrence lookup followed rejected source preflight"};
  }
  WVKernelStatus observationBatch(const WVObservationOccurrenceIdentity &,
                                  const WVObserverRecord &,
                                  WVObservationBatch &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "batch construction followed rejected source preflight"};
  }

  std::size_t preflightCount() const noexcept { return preflightCount_; }
  std::size_t unexpectedCallCount() const noexcept {
    return unexpectedCallCount_;
  }

private:
  std::size_t preflightCount_ = 0;
  mutable std::size_t unexpectedCallCount_ = 0;
};

class FailOnceNetCDFSink final : public WVOutputSink {
public:
  struct Attempt {
    std::size_t eventOrdinal = 0;
    std::string fileIdentifier;
    std::string groupIdentifier;
  };

  FailOnceNetCDFSink(WVModelOutputNetCDFSink &sink, std::size_t failAttempt)
      : sink_(sink), failAttempt_(failAttempt) {}

  WVKernelStatus preflight(const WVOutputPlan &plan) override {
    return sink_.preflight(plan);
  }

  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    ++attemptCount_;
    attempts_.push_back({event.eventOrdinal,
                         std::string(route.fileIdentifier),
                         std::string(route.groupIdentifier)});
    if (attemptCount_ == failAttempt_) {
      failAttempt_ = std::numeric_limits<std::size_t>::max();
      return {WVKernelStatusCode::invalidConfiguration,
              "Injected one-shot NetCDF route failure."};
    }
    return sink_.deliver(event, route, result);
  }

  const std::vector<Attempt> &attempts() const noexcept { return attempts_; }

private:
  WVModelOutputNetCDFSink &sink_;
  std::size_t failAttempt_ = 0;
  std::size_t attemptCount_ = 0;
  std::vector<Attempt> attempts_;
};

struct TemporaryDirectory {
  std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("wave-vortex-model-output-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  TemporaryDirectory() { std::filesystem::create_directories(path); }
  ~TemporaryDirectory() {
    if (std::getenv("WV_KEEP_MODEL_OUTPUT_FIXTURE") != nullptr) {
      std::cout << "Retained model-output fixture at " << path << '\n';
      return;
    }
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
};

WVCheckpoint checkpointTemplate() {
  WVCheckpoint checkpoint;
  const char *overridePath = std::getenv("WV_RUNTIME_CHECKPOINT_FIXTURE");
  const auto path = overridePath == nullptr
                        ? fixture("forcing-nonlinear.nc").string()
                        : std::string(overridePath);
  const auto result = WVCheckpointReader::read(path, *modelOutputCatalog(), checkpoint);
  require(static_cast<bool>(result), result.message);
  return checkpoint;
}

WVPortableObserverRecord recordFor(const WVCheckpoint &checkpoint,
                                   const std::filesystem::path &path) {
  WVPortableObserverRecord record;
  const std::vector<std::size_t> shape = {
      checkpoint.state.coefficients.shape.rows,
      checkpoint.state.coefficients.shape.columns};
  for (const auto *name : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back({name, WVStateScalarType::complex64, shape,
                                  WVToleranceKind::coefficientEnergyScaled,
                                  1e-6, WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "wave-vortex coefficient flux";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers.push_back(coefficients);
  record.outputFiles = {{"primary",
                         path.string(),
                         {{"restart",
                           "wave-vortex",
                           {1.0, checkpoint.state.t, checkpoint.state.t + 2.0},
                           {"coefficients"},
                           true}}}};
  return record;
}

WVPortableObserverDescriptor
descriptorFor(const WVPortableObserverRecord &record) {
  WVPortableObserverDescriptor descriptor;
  const auto status = WVPortableObserverDescriptor::create(record, modelOutputCatalog(), descriptor);
  require(static_cast<bool>(status), status.message);
  return descriptor;
}

WVIntegrationState eventState(
    const WVCheckpoint &checkpoint, double scheduledTime,
    const WVAdditionalStateBlockConstView *additionalBlocks = nullptr,
    std::size_t additionalBlockCount = 0) {
  auto waveVortex = checkpoint.state.view();
  waveVortex.t = scheduledTime;
  return {waveVortex, additionalBlocks, additionalBlockCount};
}

std::vector<WVOutputScheduleContinuation> continuationsFrom(
    const std::vector<WVOutputDestinationProgress> &progress) {
  std::vector<WVOutputScheduleContinuation> result;
  result.reserve(progress.size());
  for (const auto &entry : progress)
    result.push_back({entry.fileIdentifier, entry.groupIdentifier,
                      entry.committedScheduleCursor});
  return result;
}

struct RestoredInspectionState {
  WVCheckpoint checkpoint;
  WVIntegrationStateLayout layout;
  WVAdditionalStateStorage additionalState;
};

RestoredInspectionState restoreInspection(
    const WVModelOutputNetCDFInspection &inspection,
    const WVPortableObserverDescriptor &descriptor) {
  RestoredInspectionState result;
  auto status = WVIntegrationStateLayout::create(
      inspection.latestRestart.coefficientShape, descriptor, result.layout);
  require(static_cast<bool>(status), status.message);
  const auto restored = WVModelOutputNetCDFSink::restoreState(
      inspection, *modelOutputCatalog(), result.layout, result.checkpoint,
      result.additionalState);
  require(static_cast<bool>(restored), restored.message);
  return result;
}

class ZeroSampleSource final : public WVObserverSampleSource {
public:
  explicit ZeroSampleSource(
      WVTransformConstantStratificationConfiguration configuration)
      : configuration_(configuration) {}

  WVKernelStatus observationSchema(
      const WVObserverRecord &observer, WVObservationSchema &output) override {
    WVObservationSchema schema;
    schema.identifier = "legacy-" + observer.identifier + "-observation-v1";
    schema.preservesLegacyEncoding = true;
    WVObservationVariable variable;
    std::vector<std::size_t> extents;
    if (observer.typeIdentifier == "WVEulerianFields" &&
        std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  "u") != observer.fieldNames.end()) {
      variable.identifier = observer.identifier + "-u";
      variable.name = "u";
      variable.dimensionIdentifiers = {"x", "y", "z"};
      extents = {configuration_.Nx, configuration_.Ny, configuration_.Nz};
      variable.units = "m s-1";
      variable.description = "x-component of the fluid velocity";
    } else if (observer.typeIdentifier == "WVLagrangianParticles" &&
               !observer.fieldNames.empty()) {
      variable.identifier = observer.identifier + "-u";
      variable.name = observer.name + "_u";
      variable.dimensionIdentifiers = {observer.name + "_id"};
      extents = {observer.x.size()};
      variable.units = "m s-1";
      variable.description =
          "x-component of the fluid velocity, recorded along the particle "
          "trajectory";
    } else if (observer.typeIdentifier == "WVMooring" &&
               !observer.fieldNames.empty()) {
      variable.identifier = observer.identifier + "-u";
      variable.name = observer.name + "_u";
      variable.dimensionIdentifiers = {observer.name + "_z",
                                       observer.name + "_id"};
      extents = {observer.z.size(), observer.x.size()};
      variable.units = "m s-1";
      variable.description =
          "x-component of the fluid velocity, recorded at the mooring";
    } else if (observer.typeIdentifier ==
               "WVTestPortablePointDiagnostic") {
      variable.identifier = observer.identifier + "-value";
      variable.name = observer.name + "_value";
      variable.dimensionIdentifiers = {observer.name + "_id"};
      extents = {observer.x.size()};
      variable.units = "m s-1";
      variable.description = "affinely transformed point diagnostic";
    }
    if (!variable.identifier.empty()) {
      variable.scalarType = WVObservationScalarType::real64;
      variable.layout = WVObservationValueLayout::record;
      for (std::size_t index = 0;
           index < variable.dimensionIdentifiers.size(); ++index)
        schema.axes.push_back({variable.dimensionIdentifiers[index],
                               variable.dimensionIdentifiers[index],
                               WVObservationAxisKind::fixed, extents[index],
                               WVObservationCoordinateRole::none});
      std::sort(schema.axes.begin(), schema.axes.end(),
                [](const auto &left, const auto &right) {
                  return left.identifier < right.identifier;
                });
      schema.variables.push_back(std::move(variable));
    }
    const auto status = validateObservationSchema(schema);
    if (!status)
      return status;
    output = std::move(schema);
    return WVKernelStatus::ok();
  }

  WVKernelStatus prepare(const WVOutputEvent &event) override {
    ++preparationGeneration_;
    if (preparationGeneration_ == 0)
      ++preparationGeneration_;
    scheduledTime_ = event.scheduledTime;
    return WVKernelStatus::ok();
  }

  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const override {
    return makeTestOccurrenceIdentity(this, preparationGeneration_,
                                      scheduledTime_, route, observer, output);
  }

  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &,
      const WVObserverRecord &observer, WVObservationBatch &output) override {
    WVObservationSchema schema;
    auto status = observationSchema(observer, schema);
    if (!status)
      return status;
    WVObservationBatch batch;
    batch.schemaIdentifier = schema.identifier;
    batch.schemaVersion = schema.version;
    batch.kind = WVObservationBatchKind::event;
    for (std::size_t variableIndex = 0;
         variableIndex < schema.variables.size(); ++variableIndex) {
      const auto &variable = schema.variables[variableIndex];
      std::vector<std::size_t> variableExtents;
      variableExtents.reserve(variable.dimensionIdentifiers.size());
      for (const auto &dimensionIdentifier : variable.dimensionIdentifiers) {
        const auto axis = std::find_if(
            schema.axes.begin(), schema.axes.end(), [&](const auto &candidate) {
              return candidate.identifier == dimensionIdentifier;
            });
        if (axis == schema.axes.end())
          return {WVKernelStatusCode::invalidConfiguration,
                  "Zero sample source has an unresolved observation axis."};
        variableExtents.push_back(axis->extent);
      }
      std::size_t count = 1;
      for (const auto extent : variableExtents)
        count *= extent;
      auto &storage = values_[variable.identifier];
      storage.assign(count, 0.0);
      batch.values.push_back(WVObservationValue::borrowReal(
          variable.identifier, std::move(variableExtents), storage.data()));
      batch.values.back().resolvedVariableIndex = variableIndex;
    }
    status = validateObservationBatch(schema, batch);
    if (!status)
      return status;
    output = std::move(batch);
    return WVKernelStatus::ok();
  }

private:
  WVTransformConstantStratificationConfiguration configuration_;
  double scheduledTime_ = 0.0;
  std::uint64_t preparationGeneration_ = 0;
  std::map<std::string, std::vector<double>> values_;
};

void deliverPlannedEvent(WVModelOutputNetCDFSink &sink,
                         const WVOutputPlan &plan, std::size_t eventIndex,
                         const WVCheckpoint &checkpoint) {
  const auto planned = plan.event(eventIndex);
  require(planned.routeCount == 1, "test plan route count");
  WVOutputEvent event;
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state = eventState(checkpoint, planned.scheduledTime);
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  WVOutputDeliveryResult delivery;
  const auto status = sink.deliver(event, planned.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  require(delivery.writeCount == 7 || delivery.writeCount == 9,
          "coefficient delivery did not report six components and time");
}

void requireTimeSeries(const std::filesystem::path &path,
                       const std::vector<double> &expected) {
  int root = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &root) == NC_NOERR,
          "open output for time inspection");
  int group = -1;
  require(nc_inq_ncid(root, "wave-vortex", &group) == NC_NOERR,
          "locate output group");
  int variable = -1;
  require(nc_inq_varid(group, "t", &variable) == NC_NOERR,
          "locate time variable");
  std::vector<double> actual(expected.size());
  require(nc_get_var_double(group, variable, actual.data()) == NC_NOERR,
          "read time variable");
  require(actual == expected, "output time lattice changed");
  require(nc_close(root) == NC_NOERR, "close time inspection");
}

void requireLegacyScheduleSchema(const std::filesystem::path &path) {
  int root = -1, group = -1, variable = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &root) == NC_NOERR &&
              nc_inq_ncid(root, "wave-vortex", &group) == NC_NOERR,
          "open legacy schedule schema");
  require(nc_inq_att(group, NC_GLOBAL, "portableScheduleTypeIdentifier",
                     nullptr, nullptr) == NC_ENOTATT &&
              nc_inq_varid(group, "portableScheduleOrdinal", &variable) ==
                  NC_ENOTVAR &&
              nc_inq_varid(group, "portableScheduleCursor", &variable) ==
                  NC_ENOTVAR,
          "legacy evenly-spaced schema gained algorithmic metadata");
  for (const auto *name : {"outputInterval", "initialTime", "finalTime"})
    require(nc_inq_varid(group, name, &variable) == NC_NOERR,
            "legacy schedule scalar was removed");
  require(nc_close(root) == NC_NOERR, "close legacy schedule schema");
}

void testCreateReadAndAppend() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto path = directory.path / "multigroup.nc";
  auto record = recordFor(checkpoint, path);
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan firstPlan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(),
                                checkpoint.state.t,
                                checkpoint.state.t + 1.0, {}, firstPlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFConfiguration mismatchedConfiguration{
      test::extensionCatalog(), checkpoint, false};
  WVModelOutputNetCDFSink mismatchedSink;
  auto mismatchedStatus = WVModelOutputNetCDFSink::createNew(
      mismatchedConfiguration, descriptor, firstPlan, layout, nullptr,
      mismatchedSink);
  require(mismatchedStatus.code == WVCheckpointStatusCode::schemaMismatch,
          "new output sink accepted a descriptor from another catalog");
  mismatchedStatus = WVModelOutputNetCDFSink::replaceExisting(
      mismatchedConfiguration, descriptor, firstPlan, layout, nullptr,
      mismatchedSink);
  require(mismatchedStatus.code == WVCheckpointStatusCode::schemaMismatch,
          "replacement output sink accepted a descriptor from another "
          "catalog");
  mismatchedStatus = WVModelOutputNetCDFSink::openAppend(
      mismatchedConfiguration, descriptor, firstPlan, layout, nullptr, {},
      mismatchedSink);
  require(mismatchedStatus.code == WVCheckpointStatusCode::schemaMismatch &&
              !std::filesystem::exists(path),
          "append output sink accepted a descriptor from another catalog");
  WVModelOutputNetCDFConfiguration configuration{modelOutputCatalog(), checkpoint, false};
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      configuration, descriptor, firstPlan, layout, nullptr, sink);
  require(static_cast<bool>(persistence), persistence.message);
  auto incompatibleRecord = record;
  incompatibleRecord.outputFiles.front().groups.front().identifier =
      "different-group";
  auto incompatibleDescriptor = descriptorFor(incompatibleRecord);
  WVOutputPlan incompatiblePlan;
  status = WVOutputPlan::create(incompatibleDescriptor, modelOutputCatalog(), checkpoint.state.t,
                                checkpoint.state.t + 1.0, {}, incompatiblePlan);
  require(static_cast<bool>(status), status.message);
  status = sink.preflight(incompatiblePlan);
  require(!status, "NetCDF sink accepted an incompatible output plan");
  status = sink.preflight(firstPlan);
  require(static_cast<bool>(status), status.message);
  checkpoint.state.t = firstPlan.event(0).scheduledTime;
  deliverPlannedEvent(sink, firstPlan, 0, checkpoint);
  checkpoint.state.t = firstPlan.event(1).scheduledTime;
  deliverPlannedEvent(sink, firstPlan, 1, checkpoint);
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);
  require(sink.metrics().synchronizationCount == 2 &&
              sink.metrics().writtenBytes > 0 &&
              sink.metrics().payloadWriteSeconds >= 0.0 &&
              sink.metrics().synchronizationSeconds >= 0.0,
          "model-output timing and byte metrics are incomplete");
  requireTimeSeries(path, {firstPlan.event(0).scheduledTime,
                           firstPlan.event(1).scheduledTime});
  requireLegacyScheduleSchema(path);

  WVCheckpoint restored;
  persistence = WVCheckpointReader::read(path.string(), *modelOutputCatalog(), restored);
  require(static_cast<bool>(persistence), persistence.message);
  require(restored.state.t == checkpoint.state.t,
          "latest output state was not restartable");

  WVOutputPlan appendPlan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), checkpoint.state.t,
                                checkpoint.state.t + 1.0,
                                continuationsFrom(
                                    sink.destinationProgress()),
                                appendPlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, appendPlan, layout, nullptr,
      sink.destinationProgress(), append);
  require(static_cast<bool>(persistence), persistence.message);
  require(appendPlan.eventCount() == 1,
          "append plan repeated a committed output time");
  status = append.preflight(appendPlan);
  require(static_cast<bool>(status), status.message);
  checkpoint.state.t = appendPlan.event(0).scheduledTime;
  deliverPlannedEvent(append, appendPlan, 0, checkpoint);
  persistence = append.close();
  require(static_cast<bool>(persistence), persistence.message);
  requireTimeSeries(
      path, {restored.state.t - 1.0, restored.state.t, checkpoint.state.t});

  WVModelOutputNetCDFInspection inspection;
  InspectionFactoryCounts factoryCounts;
  const auto trapCatalog = inspectionTrapCatalog(factoryCounts);
  WVModelOutputNetCDFInspection rawInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {path.string()}, *trapCatalog, rawInspection);
  activeInspectionFactoryCounts = nullptr;
  require(static_cast<bool>(persistence) && factoryCounts.observers == 0 &&
              factoryCounts.schedules == 0 && factoryCounts.forcings == 0,
          "raw inspection constructed semantic extension implementations");
  persistence = WVModelOutputNetCDFSink::inspect({path.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence), persistence.message);
  require(inspection.latestRestart.t == checkpoint.state.t,
          "output inspection selected the wrong restart state");
  require(inspection.observerRecord.outputFiles.size() == 1 &&
              inspection.observerRecord.observers.size() == 1 &&
              inspection.scheduleContinuations.size() == 1 &&
              inspection.destinationProgress.size() == 1,
          "output inspection did not reconstruct the observer graph");

  int file = -1;
  require(nc_open(path.c_str(), NC_WRITE, &file) == NC_NOERR,
          "open output for interrupted-record injection");
  int group = -1;
  require(nc_inq_ncid(file, "wave-vortex", &group) == NC_NOERR,
          "locate interrupted output group");
  int variable = -1;
  require(nc_inq_varid(group, "Ap_real", &variable) == NC_NOERR,
          "locate interrupted output variable");
  const std::size_t start[] = {2, 0, 0};
  const std::size_t count[] = {1, 1, 1};
  const double fill = NC_FILL_DOUBLE;
  require(nc_put_vara_double(group, variable, start, count, &fill) == NC_NOERR,
          "inject interrupted output value");
  require(nc_close(file) == NC_NOERR, "close interrupted-record injection");
  WVModelOutputNetCDFSink rejected;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, appendPlan, layout, nullptr,
      inspection.destinationProgress, rejected);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::incompleteRecord,
          "append accepted an incomplete committed output record");
  WVModelOutputNetCDFInspection rejectedInspection;
  persistence =
      WVModelOutputNetCDFSink::inspect({path.string()}, *modelOutputCatalog(), rejectedInspection);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::incompleteRecord,
          "inspection accepted an incomplete committed output record");
}

void testLinearInitialCoefficientsAndPassiveFields() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  checkpoint.state.coefficients.Ap[1] = {1.0, 2.0};
  checkpoint.state.coefficients.Am[1] = {3.0, -4.0};
  checkpoint.state.coefficients.A0[1] = {5.0, 6.0};
  const auto path = directory.path / "linear-passive.nc";
  auto record = recordFor(checkpoint, path);
  record.observers.clear();
  WVObserverRecord fields;
  fields.identifier = "eulerian-fields-Ap-Am-A0-u";
  fields.name = "WVTestFields";
  fields.typeIdentifier = "WVTestFields";
  fields.fieldNames = {"Ap", "Am", "A0", "u", "psi"};
  record.observers.push_back(fields);
  WVObserverRecord mooring;
  mooring.identifier = "mooring-central";
  mooring.name = "central";
  mooring.typeIdentifier = "WVMooring";
  mooring.fieldNames = {"u"};
  mooring.x = {-1.0, checkpoint.configuration.Lx};
  mooring.y = {checkpoint.configuration.Ly, 0.5 * checkpoint.configuration.Ly};
  record.observers.push_back(mooring);
  WVObserverRecord diagnostic;
  diagnostic.identifier = "point-diagnostic";
  diagnostic.name = "diagnostic";
  diagnostic.typeIdentifier = "WVTestPortablePointDiagnostic";
  diagnostic.fieldNames = {"u"};
  diagnostic.x = {100.0, 200.0};
  diagnostic.y = {300.0, 400.0};
  diagnostic.z = {-100.0, -200.0};
  diagnostic.outputScale = 2.5;
  diagnostic.outputOffset = -1.25;
  record.observers.push_back(diagnostic);
  record.outputFiles.front().groups.front().observerIdentifiers = {
      fields.identifier, mooring.identifier, diagnostic.identifier};
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, true, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFConfiguration configuration{modelOutputCatalog(), checkpoint, true};
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(),
                                checkpoint.state.t, checkpoint.state.t, {},
                                plan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      configuration, descriptor, plan, layout, source.get(), sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  const auto planned = plan.event(0);
  WVOutputEvent event;
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state = eventState(checkpoint, planned.scheduledTime);
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  WVOutputDeliveryResult delivery;
  status = sink.deliver(event, planned.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  require(delivery.writeCount == 4,
          "linear delivery must write Eulerian, mooring, diagnostic, and time");
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  int root = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &root) == NC_NOERR,
          "open linear passive output");
  int group = -1;
  require(nc_inq_ncid(root, "wave-vortex", &group) == NC_NOERR,
          "locate linear passive group");
  int variable = -1;
  require(nc_inq_varid(group, "Ap_real", &variable) == NC_NOERR,
          "locate initial coefficient");
  int rank = 0;
  require(nc_inq_varndims(group, variable, &rank) == NC_NOERR && rank == 2,
          "linear coefficient must be initial-only [kl,j]");
  require(nc_inq_varid(group, "u", &variable) == NC_NOERR, "locate Eulerian u");
  require(nc_inq_varndims(group, variable, &rank) == NC_NOERR && rank == 4,
          "linear u must remain [t,z,y,x]");
  require(nc_inq_varid(group, "psi", &variable) == NC_NOERR &&
              nc_inq_varndims(group, variable, &rank) == NC_NOERR && rank == 3,
          "linear psi must remain initial-only [z,y,x]");
  require(nc_inq_varid(group, "central_u", &variable) == NC_NOERR &&
              nc_inq_varndims(group, variable, &rank) == NC_NOERR && rank == 3,
          "mooring field must use [t,id,z] NetCDF order");
  require(nc_inq_varid(group, "central_x", &variable) == NC_NOERR,
          "mooring x coordinate is absent");
  std::vector<double> x(2);
  require(nc_get_var_double(group, variable, x.data()) == NC_NOERR &&
              x[0] == checkpoint.configuration.Lx - 1.0 && x[1] == 0.0,
          "mooring periodic x coordinates changed");
  require(nc_inq_varid(group, "diagnostic_value", &variable) == NC_NOERR &&
              nc_inq_varndims(group, variable, &rank) == NC_NOERR && rank == 2,
          "point-diagnostic output must use [t,id] NetCDF order");
  require(nc_close(root) == NC_NOERR, "close linear passive output");

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect({path.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence),
          persistence.message + " at " + persistence.location);
  require(inspection.observerRecord.observers.size() == 3 &&
              inspection.observerRecord.observers.front().typeIdentifier ==
                  "WVTestFields" &&
              inspection.observerRecord.observers.front().name ==
                  "WVTestFields",
          "registered field observer did not round-trip through NetCDF");
  const auto restoredDiagnostic = std::find_if(
      inspection.observerRecord.observers.begin(),
      inspection.observerRecord.observers.end(), [](const auto &observer) {
        return observer.typeIdentifier == "WVTestPortablePointDiagnostic";
      });
  require(restoredDiagnostic != inspection.observerRecord.observers.end() &&
              restoredDiagnostic->fieldNames ==
                  std::vector<std::string>({"u"}) &&
              restoredDiagnostic->x == diagnostic.x &&
              restoredDiagnostic->y == diagnostic.y &&
              restoredDiagnostic->z == diagnostic.z &&
              restoredDiagnostic->outputScale == diagnostic.outputScale &&
              restoredDiagnostic->outputOffset == diagnostic.outputOffset,
          "point-diagnostic configuration did not round-trip through NetCDF");

  require(nc_open(path.c_str(), NC_WRITE, &root) == NC_NOERR &&
              nc_inq_ncid(root, "wave-vortex", &group) == NC_NOERR &&
              nc_redef(root) == NC_NOERR,
          "open legacy MATLAB attribute spelling injection");
  for (const auto &[name, description] :
       std::array<std::pair<const char *, const char *>, 2>{{
           {"central_x", "x coordinate position of mooring"},
           {"central_y", "y coordinate position of mooring"}}}) {
    require(nc_inq_varid(group, name, &variable) == NC_NOERR &&
                nc_put_att_text(group, variable, "long_name",
                                std::char_traits<char>::length(description),
                                description) == NC_NOERR,
            "inject legacy MATLAB long_name spelling");
  }
  require(nc_enddef(root) == NC_NOERR && nc_close(root) == NC_NOERR,
          "close legacy MATLAB attribute spelling injection");

  auto appendDescriptor = descriptorFor(inspection.observerRecord);
  auto restoredInspection = restoreInspection(inspection, appendDescriptor);
  require(restoredInspection.checkpoint.state.coefficients.Ap[1].real == 1.0 &&
              restoredInspection.checkpoint.state.coefficients.Ap[1].imag ==
                  2.0,
          "linear initial coefficient did not round-trip");
  std::unique_ptr<WVObserverOutputEvaluationService> appendSource;
  status = WVObserverOutputEvaluationService::create(
      inspection.latestRestart.configuration, true, appendDescriptor,
      std::make_unique<WVReferenceFFTEngine>(), appendSource);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFConfiguration appendConfiguration{
      modelOutputCatalog(), restoredInspection.checkpoint, true};
  WVOutputPlan appendPlan;
  status = WVOutputPlan::create(
      appendDescriptor, modelOutputCatalog(), inspection.latestRestart.t,
      inspection.latestRestart.t + 1.0, inspection.scheduleContinuations,
      appendPlan);
  require(static_cast<bool>(status) && appendPlan.eventCount() == 1,
          "linear append plan mismatch");
  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      appendConfiguration, appendDescriptor, appendPlan,
      restoredInspection.layout, appendSource.get(),
      inspection.destinationProgress, append);
  require(static_cast<bool>(persistence), persistence.message);
  status = append.preflight(appendPlan);
  require(static_cast<bool>(status), status.message);
  const auto next = appendPlan.event(0);
  event.eventOrdinal = next.eventOrdinal;
  event.scheduledTime = next.scheduledTime;
  event.state = eventState(restoredInspection.checkpoint, next.scheduledTime);
  event.routes = next.routes;
  event.routeCount = next.routeCount;
  delivery = {};
  status = append.deliver(event, next.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  persistence = append.close();
  require(static_cast<bool>(persistence), persistence.message);
  WVCheckpoint reread;
  persistence = WVCheckpointReader::read(path.string(), *modelOutputCatalog(), reread);
  require(static_cast<bool>(persistence), persistence.message);
  require(reread.state.t == next.scheduledTime &&
              reread.state.coefficients.Ap[1].imag == 2.0,
          "linear append changed initial coefficients or latest time");
}

void testTransactionalRefusal() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto path = directory.path / "existing.nc";
  {
    int file = -1;
    require(nc_create(path.c_str(), NC_NETCDF4, &file) == NC_NOERR,
            "create collision sentinel");
    require(nc_close(file) == NC_NOERR, "close collision sentinel");
  }
  auto record = recordFor(checkpoint, path);
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(),
                                checkpoint.state.t,
                                checkpoint.state.t + 2.0, {}, plan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink sink;
  const auto result = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, descriptor, plan, layout,
      nullptr, sink);
  require(!result && result.code == WVCheckpointStatusCode::commitFailure,
          "create-new output replaced an existing destination");
}

void testDirectFactoryPreflightBeforeDestinationAccess() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const WVModelOutputNetCDFConfiguration configuration{modelOutputCatalog(),
                                                       checkpoint, false};

  const auto createPath = directory.path / "rejected-create.nc";
  auto createRecord = recordFor(checkpoint, createPath);
  auto createDescriptor = descriptorFor(createRecord);
  WVIntegrationStateLayout createLayout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, createDescriptor, createLayout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan createPlan;
  status = WVOutputPlan::create(createDescriptor, modelOutputCatalog(),
                                checkpoint.state.t, checkpoint.state.t + 2.0,
                                {}, createPlan);
  require(static_cast<bool>(status), status.message);
  auto mismatchedCreateRecord = createRecord;
  mismatchedCreateRecord.outputFiles.front().groups.front().identifier =
      "mismatched-create-group";
  auto mismatchedCreateDescriptor = descriptorFor(mismatchedCreateRecord);
  WVOutputPlan mismatchedCreatePlan;
  status = WVOutputPlan::create(
      mismatchedCreateDescriptor, modelOutputCatalog(), checkpoint.state.t,
      checkpoint.state.t + 2.0, {}, mismatchedCreatePlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink mismatchedCreate;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      configuration, createDescriptor, mismatchedCreatePlan, createLayout,
      nullptr, mismatchedCreate);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::schemaMismatch &&
              !std::filesystem::exists(createPath),
          "direct create factory touched its destination for a mismatched "
          "compiled plan");
  FactoryPreflightRejectingSampleSource createSource;
  WVModelOutputNetCDFSink rejectedCreate;
  persistence = WVModelOutputNetCDFSink::createNew(
      configuration, createDescriptor, createPlan, createLayout, &createSource,
      rejectedCreate);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::unsupportedObserver &&
              createSource.preflightCount() == 1 &&
              createSource.unexpectedCallCount() == 0 &&
              !std::filesystem::exists(createPath),
          "direct create factory touched its destination after rejected "
          "source preflight");

  const auto existingPath = directory.path / "rejected-existing.nc";
  auto existingRecord = recordFor(checkpoint, existingPath);
  auto existingDescriptor = descriptorFor(existingRecord);
  WVIntegrationStateLayout existingLayout;
  status = WVIntegrationStateLayout::create(checkpoint.state.coefficients.shape,
                                            existingDescriptor, existingLayout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan existingPlan;
  status = WVOutputPlan::create(existingDescriptor, modelOutputCatalog(),
                                checkpoint.state.t, checkpoint.state.t + 2.0,
                                {}, existingPlan);
  require(static_cast<bool>(status), status.message);
  auto mismatchedExistingRecord = existingRecord;
  mismatchedExistingRecord.outputFiles.front().groups.front().identifier =
      "mismatched-existing-group";
  auto mismatchedExistingDescriptor = descriptorFor(mismatchedExistingRecord);
  WVOutputPlan mismatchedExistingPlan;
  status = WVOutputPlan::create(
      mismatchedExistingDescriptor, modelOutputCatalog(), checkpoint.state.t,
      checkpoint.state.t + 2.0, {}, mismatchedExistingPlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink seed;
  persistence = WVModelOutputNetCDFSink::createNew(
      configuration, existingDescriptor, existingPlan, existingLayout, nullptr,
      seed);
  require(static_cast<bool>(persistence), persistence.message);
  const auto destinationProgress = seed.destinationProgress();
  persistence = seed.close();
  require(static_cast<bool>(persistence), persistence.message);
  const auto originalBytes = fileBytes(existingPath);

  WVModelOutputNetCDFSink mismatchedReplacement;
  persistence = WVModelOutputNetCDFSink::replaceExisting(
      configuration, existingDescriptor, mismatchedExistingPlan, existingLayout,
      nullptr, mismatchedReplacement);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::schemaMismatch &&
              fileBytes(existingPath) == originalBytes,
          "direct replace factory changed its destination for a mismatched "
          "compiled plan");
  FactoryPreflightRejectingSampleSource replaceSource;
  WVModelOutputNetCDFSink rejectedReplacement;
  persistence = WVModelOutputNetCDFSink::replaceExisting(
      configuration, existingDescriptor, existingPlan, existingLayout,
      &replaceSource, rejectedReplacement);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::unsupportedObserver &&
              replaceSource.preflightCount() == 1 &&
              replaceSource.unexpectedCallCount() == 0 &&
              fileBytes(existingPath) == originalBytes,
          "direct replace factory changed its destination after rejected "
          "source preflight");

  WVOutputPlan appendPlan;
  status =
      WVOutputPlan::create(existingDescriptor, modelOutputCatalog(),
                           checkpoint.state.t, checkpoint.state.t + 2.0,
                           continuationsFrom(destinationProgress), appendPlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink mismatchedAppend;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, existingDescriptor, mismatchedExistingPlan, existingLayout,
      nullptr, destinationProgress, mismatchedAppend);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::schemaMismatch &&
              fileBytes(existingPath) == originalBytes,
          "direct append factory changed its destination for a mismatched "
          "compiled plan");
  FactoryPreflightRejectingSampleSource appendSource;
  WVModelOutputNetCDFSink rejectedAppend;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, existingDescriptor, appendPlan, existingLayout,
      &appendSource, destinationProgress, rejectedAppend);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::unsupportedObserver &&
              appendSource.preflightCount() == 1 &&
              appendSource.unexpectedCallCount() == 0 &&
              fileBytes(existingPath) == originalBytes,
          "direct append factory changed its destination after rejected "
          "source preflight");
}

void testMultipleFilesGroupsAndSharedState() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  WVPortableObserverRecord record;
  const std::vector<std::size_t> coefficientShape = {
      checkpoint.state.coefficients.shape.rows,
      checkpoint.state.coefficients.shape.columns};
  for (const auto *name : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back({name, WVStateScalarType::complex64,
                                  coefficientShape,
                                  WVToleranceKind::coefficientEnergyScaled,
                                  1e-6, WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
  for (const auto *name : {"particles-x", "particles-y"})
    record.stateBlocks.push_back({name,
                                  WVStateScalarType::real64,
                                  {2},
                                  WVToleranceKind::uniformAbsolute,
                                  1e-4,
                                  WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
  record.stateBlocks.push_back(
      {"tracer-state",
       WVStateScalarType::real64,
       {checkpoint.configuration.Nx, checkpoint.configuration.Ny,
        checkpoint.configuration.Nz},
       WVToleranceKind::uniformAbsolute,
       1e-5,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  WVObserverRecord particles;
  particles.identifier = "particles";
  particles.name = "particles";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.stateBlockIdentifiers = {"particles-x", "particles-y"};
  particles.x = {0.1, 0.2};
  particles.y = {0.3, 0.4};
  particles.z = {-100.0, -300.0};
  particles.isXYOnly = true;
  particles.horizontalAbsoluteTolerance = 1e-4;
  WVObserverRecord tracer;
  tracer.identifier = "tracer";
  tracer.name = "tracer";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"tracer-state"};
  tracer.isXYOnly = false;
  tracer.shouldAntialias = true;
  WVObserverRecord fields;
  fields.identifier = "fields";
  fields.name = "WVEulerianFields";
  fields.typeIdentifier = "WVEulerianFields";
  fields.fieldNames = {"u", "v", "rho_e"};
  WVObserverRecord mooring;
  mooring.identifier = "mooring";
  mooring.name = "mooring";
  mooring.typeIdentifier = "WVMooring";
  mooring.fieldNames = {"u", "v"};
  mooring.x = {0.0, 0.5 * checkpoint.configuration.Lx};
  mooring.y = {0.0, 0.5 * checkpoint.configuration.Ly};
  record.observers = {coefficients, fields, mooring, particles, tracer};
  const auto outputInterval = 1e-7;
  const auto end = checkpoint.state.t + 2.0 * outputInterval;
  const char *exportPath = std::getenv("WV_RUNTIME_MODEL_OUTPUT_EXPORT");
  const auto first = exportPath == nullptr ? directory.path / "first.nc"
                                           : std::filesystem::path(exportPath);
  const auto second =
      exportPath == nullptr
          ? directory.path / "second.nc"
          : std::filesystem::path(std::string(exportPath) + ".second.nc");
  record.outputFiles = {
      {"first",
       first.string(),
       {{"restart",
         "wave-vortex",
         {outputInterval, checkpoint.state.t, end},
         {"coefficients", "fields", "mooring"},
         true},
        {"shared",
         "shared",
         {outputInterval, checkpoint.state.t, end},
         {"fields", "mooring", "particles", "tracer"},
         false}}},
      {"second",
       second.string(),
       {{"restart",
         "restart",
         {outputInterval, checkpoint.state.t, end},
         {"coefficients", "fields", "mooring", "particles", "tracer"},
         true}}}};
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  WVAdditionalStateStorage additional;
  status = additional.initialize(layout);
  require(static_cast<bool>(status), status.message);
  for (std::size_t block = 0; block < additional.blockCount(); ++block)
    std::fill_n(additional.mutableBlocks()[block].realData,
                additional.mutableBlocks()[block].layout->elementCount,
                static_cast<double>(block + 1));
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(),
                                checkpoint.state.t, end, {}, plan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, descriptor, plan, layout,
      source.get(), sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  WVOutputEvent event;
  const auto planned = plan.event(0);
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state = eventState(checkpoint, planned.scheduledTime,
                           additional.constBlocks(), additional.blockCount());
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  require(planned.routeCount == 3,
          "coincident multi-file routes were not grouped");
  for (std::size_t route = 0; route < planned.routeCount; ++route) {
    WVOutputDeliveryResult delivery;
    status = sink.deliver(event, planned.routes[route], delivery);
    require(static_cast<bool>(status), status.message);
  }
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);
  {
    int file = -1;
    int group = -1;
    int ids = -1;
    int z = -1;
    int tracerVariable = -1;
    require(nc_open(first.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
                nc_inq_ncid(file, "shared", &group) == NC_NOERR &&
                nc_inq_varid(group, "particles_id", &ids) == NC_NOERR &&
                nc_inq_varid(group, "particles_z", &z) == NC_NOERR &&
                nc_inq_varid(group, "tracer", &tracerVariable) == NC_NOERR,
            "particle MATLAB schema variables missing");
    std::array<double, 2> idValues{};
    require(nc_get_var_double(group, ids, idValues.data()) == NC_NOERR &&
                idValues == std::array<double, 2>{{1.0, 2.0}},
            "particle identifiers are not one-based");
    const std::size_t start[] = {0, 0};
    const std::size_t count[] = {1, 2};
    std::array<double, 2> zValues{};
    require(nc_get_vara_double(group, z, start, count, zValues.data()) ==
                    NC_NOERR &&
                zValues == std::array<double, 2>{{-100.0, -300.0}},
            "fixed particle z values were not written on [id,t]");
    char attribute[64]{};
    require(nc_get_att_text(group, z, "particleName", attribute) == NC_NOERR &&
                std::string(attribute, std::string("particles").size()) ==
                    "particles",
            "particle metadata attributes missing");
    std::fill(std::begin(attribute), std::end(attribute), '\0');
    require(nc_get_att_text(group, tracerVariable, "isTracer", attribute) ==
                    NC_NOERR &&
                attribute[0] == '1',
            "tracer marker attribute missing");
    require(nc_close(file) == NC_NOERR, "close particle schema inspection");
  }
  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {first.string(), second.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence), persistence.message);
  require(inspection.observerRecord.outputFiles.size() == 2 &&
              inspection.scheduleContinuations.size() == 3 &&
              inspection.destinationProgress.size() == 3 &&
              inspection.observerRecord.observers.size() == 5,
          "multi-file graph or shared observer identity was not reconstructed");
  WVCheckpoint restoredCheckpoint;
  WVAdditionalStateStorage restoredAdditionalState;
  persistence = WVModelOutputNetCDFSink::restoreState(
      inspection, *modelOutputCatalog(), layout, restoredCheckpoint,
      restoredAdditionalState);
  require(static_cast<bool>(persistence), persistence.message);
  require(restoredAdditionalState.blockCount() == 3,
          "shared dynamic state was duplicated during reconstruction");
  const auto particleX =
      std::find_if(restoredAdditionalState.constBlocks(),
                   restoredAdditionalState.constBlocks() +
                       restoredAdditionalState.blockCount(),
                   [](const auto &view) {
                     return view.layout->identifier == "particles-x";
                   });
  const auto tracerState =
      std::find_if(restoredAdditionalState.constBlocks(),
                   restoredAdditionalState.constBlocks() +
                       restoredAdditionalState.blockCount(),
                   [](const auto &view) {
                     return view.layout->identifier == "tracer-state";
                   });
  require(particleX != restoredAdditionalState.constBlocks() +
                           restoredAdditionalState.blockCount() &&
              particleX->realData[0] == 1.0 &&
              tracerState != restoredAdditionalState.constBlocks() +
                                 restoredAdditionalState.blockCount() &&
              tracerState->realData[0] == 3.0,
          "dynamic state stored outside the coefficient group was not "
          "restored");
  for (std::size_t block = 0;
       block < restoredAdditionalState.blockCount();
       ++block) {
    const auto &view = restoredAdditionalState.constBlocks()[block];
    require(view.layout != nullptr && view.realData != nullptr,
            "restored dynamic-state view is invalid");
    require(std::all_of(view.realData,
                        view.realData + view.layout->elementCount,
                        [](double value) { return std::isfinite(value); }),
            "restored dynamic state contains nonfinite values");
  }
  const auto restoredParticles = std::find_if(
      inspection.observerRecord.observers.begin(),
      inspection.observerRecord.observers.end(),
      [](const auto &observer) { return observer.identifier == "particles"; });
  require(restoredParticles != inspection.observerRecord.observers.end() &&
              restoredParticles->z == std::vector<double>({-100.0, -300.0}),
          "fixed XY particle z configuration was not reconstructed");
  {
    const auto sourceFirstBytes = fileBytes(first);
    const auto sourceSecondBytes = fileBytes(second);
    const auto requestPath = directory.path / "portable-run.json";
    const auto requestedFirst = directory.path / "requested-first.nc";
    const auto requestedSecond = directory.path / "requested-second.nc";
    const auto requestedReport = directory.path / "portable-run-report.json";
    std::ofstream request(requestPath, std::ios::binary | std::ios::trunc);
    request << std::setprecision(17)
            << "{\"schemaIdentifier\":\"wave-vortex-run-request-v1\","
               "\"schemaVersion\":1,\"modelFiles\":["
            << "\"" << first.string() << "\",\"" << second.string()
            << "\"],\"integration\":{\"method\":\"fixed-rk4\","
               "\"finalTime\":"
            << end << ",\"initialStep\":" << outputInterval
            << "},\"output\":{\"policy\":\"create\",\"destinations\":{"
               "\"first\":\""
            << requestedFirst.string() << "\",\"second\":\""
            << requestedSecond.string()
            << "\"}},\"execution\":{\"fftProvider\":\"reference\","
               "\"threads\":1},\"report\":\""
            << requestedReport.string() << "\"}";
    request.close();
    std::ostringstream command;
    command << shellQuote(WV_RUNTIME_RUNNER) << " --request "
            << shellQuote(requestPath) << " >/dev/null";
    require(std::system(command.str().c_str()) == 0,
            "MATLAB-authored multi-file run request failed");
    WVModelOutputNetCDFInspection requested;
    persistence = WVModelOutputNetCDFSink::inspect(
        {requestedFirst.string(), requestedSecond.string()}, *modelOutputCatalog(), requested);
    require(static_cast<bool>(persistence), persistence.message);
    auto requestedDescriptor = descriptorFor(requested.observerRecord);
    auto requestedState = restoreInspection(requested, requestedDescriptor);
    require(std::abs(requested.latestRestart.t - end) <= 1e-14 &&
                requested.observerRecord.observers.size() == 5 &&
                requestedState.additionalState.blockCount() == 3,
            "run request did not preserve the complete model graph");
    require(fileBytes(first) == sourceFirstBytes &&
                fileBytes(second) == sourceSecondBytes,
            "create request mutated its MATLAB-authored source bundle");
    require(std::filesystem::exists(requestedReport),
            "run request did not write its report");
    const auto v1Report = fileText(requestedReport);
    require(v1Report.find(
                "\"schemaIdentifier\":\"wave-vortex-run-request-v1\"") !=
                std::string::npos &&
                v1Report.find("\"schemaVersion\":1") != std::string::npos &&
                v1Report.find("\"id\":\"fixed-rk4\"") !=
                    std::string::npos &&
                v1Report.find("\"integrationRequest\"") ==
                    std::string::npos &&
                v1Report.find("\"parse\":") == std::string::npos,
            "v1 execution or report semantics changed under v2 support");

    struct V2Run {
      std::filesystem::path request;
      std::filesystem::path first;
      std::filesystem::path second;
      std::filesystem::path report;
      std::string reportText;
    };
    const auto runV2 = [&](const std::string &name,
                           const std::string &integration) {
      V2Run result;
      result.request = directory.path / (name + "-request.json");
      result.first = directory.path / (name + "-first.nc");
      result.second = directory.path / (name + "-second.nc");
      result.report = directory.path / (name + "-report.json");
      std::ofstream document(result.request,
                             std::ios::binary | std::ios::trunc);
      document << "{\"schemaIdentifier\":\"wave-vortex-run-request-v2\","
                  "\"schemaVersion\":2,\"modelFiles\":[\""
               << first.string() << "\",\"" << second.string()
               << "\"],\"integration\":" << integration
               << ",\"output\":{\"policy\":\"create\","
                  "\"destinations\":{\"first\":\""
               << result.first.string() << "\",\"second\":\""
               << result.second.string()
               << "\"}},\"execution\":{\"fftProvider\":\"reference\","
                  "\"threads\":1},\"report\":\""
               << result.report.string() << "\"}";
      document.close();
      const auto immutableRequest = fileBytes(result.request);
      require(std::system((shellQuote(WV_RUNTIME_RUNNER) + " --request " +
                           shellQuote(result.request) + " >/dev/null")
                              .c_str()) == 0,
              "v2 request execution failed: " + name);
      require(fileBytes(result.request) == immutableRequest &&
                  fileBytes(first) == sourceFirstBytes &&
                  fileBytes(second) == sourceSecondBytes,
              "v2 execution mutated its request or source bundle: " + name);
      result.reportText = fileText(result.report);
      require(result.reportText.find(
                  "\"schemaIdentifier\":\"wave-vortex-run-request-v2\"") !=
                      std::string::npos &&
                  result.reportText.find("\"schemaVersion\":2") !=
                      std::string::npos &&
                  result.reportText.find("\"requestedMethod\":") !=
                      std::string::npos &&
                  result.reportText.find("\"activeMethod\":") !=
                      std::string::npos &&
                  result.reportText.find("\"noFallback\":true") !=
                      std::string::npos &&
                  result.reportText.find("\"parse\":") !=
                      std::string::npos &&
                  result.reportText.find("\"preflight\":") !=
                      std::string::npos &&
                  result.reportText.find("\"provider\":") !=
                      std::string::npos &&
                  result.reportText.find("\"startup\":") !=
                      std::string::npos,
              "v2 report omitted contract or lifecycle evidence: " + name);
      return result;
    };

    std::unique_ptr<WVConstantStratificationIntegrationSystem> cflSystem;
    status = WVConstantStratificationIntegrationSystem::create(
        inspection.latestRestart.configuration,
        inspection.latestRestart.forcingSchedule,modelOutputCatalog(),
        std::make_unique<WVReferenceFFTEngine>(),cflSystem);
    require(static_cast<bool>(status),status.message);
    WVFixedTimeStepCandidates unitCandidates;
    status = cflSystem->evaluateFixedTimeStepCandidates(
        {restoredCheckpoint.state.view(),nullptr,0},1.0,unitCandidates);
    require(static_cast<bool>(status),status.message);
    const auto unitSelected =
        std::min(unitCandidates.advective,unitCandidates.oscillatory);
    require(std::isfinite(unitSelected) && unitSelected > 0.0,
            "v2 trajectory fixture has no finite CFL candidate");
    const auto matchedCFL = outputInterval/unitSelected;
    std::ostringstream explicitIntegration;
    explicitIntegration << std::setprecision(17)
                        << "{\"method\":\"fixed-rk4\",\"finalTime\":"
                        << end << ",\"initialStep\":" << outputInterval
                        << '}';
    const auto explicitV2 = runV2("v2-fixed-explicit",
                                  explicitIntegration.str());
    std::ostringstream cflIntegration;
    cflIntegration << std::setprecision(17)
                   << "{\"method\":\"fixed-rk4\",\"finalTime\":" << end
                   << ",\"cfl\":" << matchedCFL
                   << ",\"timeStepConstraint\":\"min\"}";
    const auto cflV2 = runV2("v2-fixed-cfl",cflIntegration.str());
    require(cflV2.reportText.find("\"stepPolicy\":\"cfl\"") !=
                std::string::npos &&
                cflV2.reportText.find("\"timeStepConstraint\":\"min\"") !=
                    std::string::npos &&
                jsonNumber(cflV2.reportText,"transientWorkspaceMaximumLiveBytes") >
                    0.0 &&
                std::abs(jsonNumber(cflV2.reportText,"selectedStep") -
                         outputInterval) <= 1e-12*outputInterval,
            "v2 CFL report omitted matched selection or transient storage");
    require(explicitV2.reportText.find("\"stepPolicy\":\"explicit\"") !=
                std::string::npos &&
                std::abs(jsonNumber(explicitV2.reportText,"selectedStep") -
                         outputInterval) <= 1e-12*outputInterval &&
                std::abs(jsonNumber(explicitV2.reportText,
                                    "requestedInitialStep") -
                         outputInterval) <= 1e-12*outputInterval,
            "v2 explicit-step report omitted its requested selection");
    WVCheckpoint explicitCheckpoint;
    WVCheckpoint cflCheckpoint;
    require(WVCheckpointReader::read(explicitV2.first.string(),
                                     *modelOutputCatalog(),
                                     explicitCheckpoint) &&
                WVCheckpointReader::read(cflV2.first.string(),
                                         *modelOutputCatalog(),cflCheckpoint),
            "matched explicit/CFL outputs are unreadable");
    require(coefficientDifference(explicitCheckpoint,cflCheckpoint) <= 1e-12,
            "matched explicit and CFL RK4 requests changed the trajectory");

    struct AdaptiveFixture {
      const char *method;
      const char *controller;
    };
    for (const auto &adaptive : {
             AdaptiveFixture{"adaptive-rk23","matlab-ode23-v1"},
             AdaptiveFixture{"adaptive-rk45","matlab-ode45-v1"},
             AdaptiveFixture{"adaptive-rk78","matlab-ode78-v1"}}) {
      std::ostringstream integration;
      integration << std::setprecision(17) << "{\"method\":\""
                  << adaptive.method << "\",\"finalTime\":" << end
                  << ",\"initialStep\":" << 2.0*outputInterval
                  << ",\"maximumStep\":" << 2.0*outputInterval
                  << ",\"relativeTolerance\":1e-6,"
                     "\"absoluteToleranceScale\":1e-8}";
      const auto result = runV2(std::string("v2-")+adaptive.method,
                                integration.str());
      require(result.reportText.find(std::string("\"requestedMethod\":\"")+
                                         adaptive.method+"\"") !=
                      std::string::npos &&
                  result.reportText.find(std::string("\"activeMethod\":\"")+
                                         adaptive.method+"\"") !=
                      std::string::npos &&
                  result.reportText.find(std::string("\"controller\":\"")+
                                         adaptive.controller+"\"") !=
                      std::string::npos &&
                  std::abs(jsonNumber(result.reportText,
                                      "requestedInitialStep") -
                           2.0*outputInterval) <= 2e-12*outputInterval &&
                  std::abs(jsonNumber(result.reportText,
                                      "requestedMaximumStep") -
                           2.0*outputInterval) <= 2e-12*outputInterval &&
                  jsonNumber(result.reportText,"relativeTolerance") == 1e-6 &&
                  jsonNumber(result.reportText,"absoluteTolerance") == 1e-8,
              std::string("v2 request did not execute exact method ")+
                  adaptive.method);
      if (std::string(adaptive.method) == "adaptive-rk78")
        require(jsonNumber(result.reportText,
                           "continuousExtensionRightHandSideEvaluationCount") ==
                        4.0 &&
                    jsonNumber(result.reportText,"denseOutputCacheBuildCount") ==
                        1.0,
                "v2 RK78 model output did not use its lazy continuous "
                "extension exactly once");
    }

    const auto writeRejectedV2 = [&](const std::string &name,
                                     const std::string &integration,
                                     const std::string &provider) {
      const auto requestPath = directory.path/(name+".json");
      const auto rejectedFirst = directory.path/(name+"-first.nc");
      const auto rejectedSecond = directory.path/(name+"-second.nc");
      std::ofstream document(requestPath,std::ios::binary|std::ios::trunc);
      document << "{\"schemaIdentifier\":\"wave-vortex-run-request-v2\","
                  "\"schemaVersion\":2,\"modelFiles\":[\""
               << first.string() << "\",\"" << second.string()
               << "\"],\"integration\":" << integration
               << ",\"output\":{\"policy\":\"create\","
                  "\"destinations\":{\"first\":\""
               << rejectedFirst.string() << "\",\"second\":\""
               << rejectedSecond.string()
               << "\"}},\"execution\":{\"fftProvider\":\"" << provider
               << "\",\"threads\":1},\"report\":\""
               << (directory.path/(name+"-report.json")).string() << "\"}";
      document.close();
      require(std::system((shellQuote(WV_RUNTIME_RUNNER)+" --request "+
                           shellQuote(requestPath)+" >/dev/null 2>&1").c_str()) !=
                      0 &&
                  !std::filesystem::exists(rejectedFirst) &&
                  !std::filesystem::exists(rejectedSecond),
              "rejected v2 request mutated a destination: "+name);
    };
    writeRejectedV2(
        "v2-mixed-controls",
        "{\"method\":\"fixed-rk4\",\"finalTime\":"+
            std::to_string(end)+
            ",\"initialStep\":1e-7,\"cfl\":0.25,"
            "\"timeStepConstraint\":\"min\"}",
        "reference");
    writeRejectedV2("v2-unavailable-provider",
                    explicitIntegration.str(),"native-fftw");

    const auto partialPath = directory.path / "partial-run.json";
    std::ofstream partial(partialPath, std::ios::binary | std::ios::trunc);
    partial << std::setprecision(17)
            << "{\"schemaIdentifier\":\"wave-vortex-run-request-v1\","
               "\"schemaVersion\":1,\"modelFiles\":[\""
            << first.string() << "\",\"" << second.string()
            << "\"],\"integration\":{\"method\":\"fixed-rk4\","
               "\"finalTime\":"
            << end << ",\"initialStep\":" << outputInterval
            << "},\"output\":{\"policy\":\"create\",\"destinations\":{"
               "\"first\":\""
            << (directory.path / "partial.nc").string()
            << "\"}},\"execution\":{\"fftProvider\":\"native-fftw\","
               "\"threads\":1},\"report\":\""
            << (directory.path / "partial-report.json").string() << "\"}";
    partial.close();
    require(std::system((shellQuote(WV_RUNTIME_RUNNER) + " --request " +
                         shellQuote(partialPath) + " >/dev/null 2>&1")
                            .c_str()) != 0,
            "incomplete destination map passed graph preflight");
    require(
        !std::filesystem::exists(directory.path / "partial.nc"),
        "failed graph preflight mutated output before provider construction");
    require(
        std::system((shellQuote(WV_RUNTIME_RUNNER) + " --request " +
                     shellQuote(requestPath) + " --threads 1 >/dev/null 2>&1")
                        .c_str()) != 0,
        "request mode accepted legacy semantic overrides");
  }
  {
    int file = -1;
    int group = -1;
    int variable = -1;
    require(nc_open(second.c_str(), NC_WRITE, &file) == NC_NOERR &&
                nc_inq_ncid(file, "restart", &group) == NC_NOERR &&
                nc_inq_varid(group, "particles_x", &variable) == NC_NOERR,
            "open duplicate particle state for ambiguity test");
    const std::size_t start[] = {0, 0};
    const std::size_t count[] = {1, 1};
    double conflicting = 99.0;
    require(nc_put_vara_double(group, variable, start, count, &conflicting) ==
                    NC_NOERR &&
                nc_close(file) == NC_NOERR,
            "write conflicting duplicate particle state");
    WVModelOutputNetCDFInspection rejected;
    const auto rejectedStatus = WVModelOutputNetCDFSink::inspect(
        {first.string(), second.string()}, *modelOutputCatalog(), rejected);
    require(!rejectedStatus &&
                rejectedStatus.code == WVCheckpointStatusCode::ambiguousState,
            "conflicting cross-group restart state was not rejected");
  }
  std::ostringstream fixedCommand;
  fixedCommand << shellQuote(WV_RUNTIME_RUNNER) << ' ' << shellQuote(first)
               << " --restart-mode model --output-policy append --delta-t "
               << outputInterval << " --final-time " << std::setprecision(17)
               << end << " --fft-provider reference >/dev/null";
  require(std::system(fixedCommand.str().c_str()) == 0,
          "fixed-step full-model continuation failed");
  WVModelOutputNetCDFInspection continued;
  persistence = WVModelOutputNetCDFSink::inspect({first.string()}, *modelOutputCatalog(), continued);
  require(static_cast<bool>(persistence), persistence.message);
  WVCheckpoint continuedCheckpoint;
  WVAdditionalStateStorage continuedAdditionalState;
  persistence = WVModelOutputNetCDFSink::restoreState(
      continued, *modelOutputCatalog(), layout, continuedCheckpoint,
      continuedAdditionalState);
  require(static_cast<bool>(persistence), persistence.message);
  require(std::abs(continued.latestRestart.t - end) <= 1e-14 &&
              continued.observerRecord.observers.size() == 5 &&
              continuedAdditionalState.blockCount() == 3,
          "fixed-step continuation did not preserve the complete model graph");

  std::cout << "OUTPUT_METRICS files=" << sink.metrics().fileCount
            << " groups=" << sink.metrics().groupCount
            << " records=" << sink.metrics().committedRecordCount
            << " syncs=" << sink.metrics().synchronizationCount
            << " bytes=" << sink.metrics().writtenBytes
            << " payload_seconds=" << sink.metrics().payloadWriteSeconds
            << " sync_seconds=" << sink.metrics().synchronizationSeconds
            << " sink_bytes=" << sink.metrics().retainedStorageBytes
            << " observer_bytes=" << source->persistentBytes() << '\n';
}

void testOptionalMatlabFixture() {
  const char *path = std::getenv("WV_MATLAB_MODEL_OUTPUT_FIXTURE");
  if (path == nullptr)
    return;
  std::vector<std::string> paths{path};
  const char *secondPath = std::getenv("WV_MATLAB_MODEL_OUTPUT_FIXTURE_SECOND");
  if (secondPath != nullptr)
    paths.emplace_back(secondPath);
  WVModelOutputNetCDFInspection inspection;
  const auto inspectStatus =
      WVModelOutputNetCDFSink::inspect(paths, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(inspectStatus),
          inspectStatus.message + " at " + inspectStatus.location);
  const auto hasType = [&](const std::string &typeIdentifier) {
    return std::any_of(inspection.observerRecord.observers.begin(),
                       inspection.observerRecord.observers.end(),
                       [&](const auto &observer) {
                         return observer.typeIdentifier == typeIdentifier;
                       });
  };
  const bool hasCoefficients = hasType("WVCoefficients");
  const bool hasEulerian = hasType("WVEulerianFields");
  const bool hasMooring = hasType("WVMooring");
  const bool hasParticles = hasType("WVLagrangianParticles");
  const bool hasTracer = hasType("WVTracer");
  require(hasEulerian && hasMooring && hasParticles && hasTracer,
          "MATLAB observer graph kinds: coefficients=" +
              std::to_string(hasCoefficients) +
              " eulerian=" + std::to_string(hasEulerian) +
              " mooring=" + std::to_string(hasMooring) +
              " particles=" + std::to_string(hasParticles) +
              " tracer=" + std::to_string(hasTracer));
  require(inspection.observerRecord.outputFiles.size() == paths.size() &&
              inspection.observerRecord.outputFiles.front().groups.size() == 2,
          "MATLAB multi-group graph was not reconstructed");
  TemporaryDirectory directory;
  const auto appendPath = directory.path / "matlab-append.nc";
  std::filesystem::copy_file(path, appendPath);
  WVModelOutputNetCDFInspection appendInspection;
  auto persistence =
      WVModelOutputNetCDFSink::inspect({appendPath.string()}, *modelOutputCatalog(), appendInspection);
  require(static_cast<bool>(persistence), persistence.message);
  appendInspection.observerRecord.outputFiles.front().destination =
      appendPath.string();
  auto descriptor = descriptorFor(appendInspection.observerRecord);
  auto restoredInspection = restoreInspection(appendInspection, descriptor);
  auto &layout = restoredInspection.layout;
  auto status = WVKernelStatus::ok();
  require(restoredInspection.additionalState.blockCount() == 4,
          "MATLAB dynamic particle/tracer state was not reconstructed");
  ZeroSampleSource samples(appendInspection.latestRestart.configuration);
  WVOutputPlan plan;
  status = WVOutputPlan::create(
      descriptor, modelOutputCatalog(), appendInspection.latestRestart.t,
      appendInspection.latestRestart.t + 1.0,
      appendInspection.scheduleContinuations, plan);
  require(static_cast<bool>(status), status.message);
  require(plan.eventCount() >= 1,
          "MATLAB append plan did not select a future schedule point");
  WVModelOutputNetCDFSink sink;
  persistence = WVModelOutputNetCDFSink::openAppend(
      {modelOutputCatalog(), restoredInspection.checkpoint, false}, descriptor,
      plan, layout, &samples, appendInspection.destinationProgress, sink);
  require(static_cast<bool>(persistence),
          persistence.message + " at " + persistence.location);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  std::vector<WVAdditionalStateBlockConstView> blocks;
  blocks.reserve(restoredInspection.additionalState.blockCount());
  for (std::size_t index = 0;
       index < restoredInspection.additionalState.blockCount(); ++index)
    blocks.push_back(restoredInspection.additionalState.constBlocks()[index]);
  WVOutputEvent event;
  const auto planned = plan.event(0);
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state = eventState(restoredInspection.checkpoint,
                           planned.scheduledTime, blocks.data(), blocks.size());
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  for (std::size_t route = 0; route < planned.routeCount; ++route) {
    WVOutputDeliveryResult delivery;
    status = sink.deliver(event, planned.routes[route], delivery);
    require(static_cast<bool>(status), status.message);
  }
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);
  WVCheckpoint appended;
  persistence = WVCheckpointReader::read(appendPath.string(), *modelOutputCatalog(), appended);
  require(static_cast<bool>(persistence), persistence.message);
  require(appended.state.t == planned.scheduledTime,
          "runtime did not append the next MATLAB output time");
}

void testOptionalMatlabLinearFixture() {
  const char *path = std::getenv("WV_MATLAB_LINEAR_OUTPUT_FIXTURE");
  if (path == nullptr)
    return;
  WVModelOutputNetCDFInspection inspection;
  auto persistence = WVModelOutputNetCDFSink::inspect({path}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence),
          persistence.message + " at " + persistence.location);
  require(inspection.observerRecord.observers.size() == 1 &&
              inspection.observerRecord.observers.front().typeIdentifier ==
                  "WVEulerianFields",
          "MATLAB linear Eulerian observer was not reconstructed");
  require(inspection.latestRestart.coefficientShape.elementCount() > 0,
          "MATLAB linear coefficients were not reconstructed");

  TemporaryDirectory directory;
  const auto appendPath = directory.path / "matlab-linear-append.nc";
  std::filesystem::copy_file(path, appendPath);
  inspection.observerRecord.outputFiles.front().destination =
      appendPath.string();
  auto descriptor = descriptorFor(inspection.observerRecord);
  auto restoredInspection = restoreInspection(inspection, descriptor);
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  auto status = WVObserverOutputEvaluationService::create(
      inspection.latestRestart.configuration, true, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(
      descriptor, modelOutputCatalog(), inspection.latestRestart.t,
      inspection.latestRestart.t + 1.0, inspection.scheduleContinuations,
      plan);
  require(static_cast<bool>(status) && plan.eventCount() == 1,
          "MATLAB linear append plan mismatch");
  WVModelOutputNetCDFSink sink;
  persistence = WVModelOutputNetCDFSink::openAppend(
      {modelOutputCatalog(), restoredInspection.checkpoint, true}, descriptor,
      plan, restoredInspection.layout, source.get(),
      inspection.destinationProgress, sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  const auto planned = plan.event(0);
  WVOutputEvent event;
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state =
      eventState(restoredInspection.checkpoint, planned.scheduledTime);
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  WVOutputDeliveryResult delivery;
  status = sink.deliver(event, planned.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);
  WVCheckpoint appended;
  persistence = WVCheckpointReader::read(appendPath.string(), *modelOutputCatalog(), appended);
  require(static_cast<bool>(persistence) &&
              appended.state.t == planned.scheduledTime,
          "MATLAB linear file did not append");
}

void testOptionalMatlabPassiveFixture() {
  const char *path = std::getenv("WV_MATLAB_PASSIVE_OUTPUT_FIXTURE");
  if (path == nullptr)
    return;
  WVModelOutputNetCDFInspection inspection;
  const auto persistence = WVModelOutputNetCDFSink::inspect({path}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence),
          persistence.message + " at " + persistence.location);
  require(inspection.observerRecord.observers.size() == 2,
          "MATLAB passive observer graph was not reconstructed");
  const auto coefficients = std::count_if(
      inspection.observerRecord.observers.begin(),
      inspection.observerRecord.observers.end(), [](const auto &observer) {
        return observer.typeIdentifier == "WVCoefficients";
      });
  const auto eulerian = std::count_if(
      inspection.observerRecord.observers.begin(),
      inspection.observerRecord.observers.end(), [](const auto &observer) {
        return observer.typeIdentifier == "WVEulerianFields";
      });
  require(coefficients == 1 && eulerian == 1,
          "MATLAB coefficient and Eulerian metadata changed");
}

void testAlgorithmicSchedulePersistence() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  const auto path = directory.path / "algorithmic.nc";
  auto record = recordFor(checkpoint, path);
  record.outputFiles[0].groups[0].schedule =
      quadraticSchedule(initialTime + 4.0, initialTime);
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFConfiguration configuration{modelOutputCatalog(), checkpoint, false};
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime, initialTime + 1.0, {},
                                plan);
  require(static_cast<bool>(status) && plan.eventCount() == 2,
          "algorithmic schedule create plan");
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      configuration, descriptor, plan, layout, nullptr, sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  for (std::size_t index = 0; index < plan.eventCount(); ++index) {
    checkpoint.state.t = plan.event(index).scheduledTime;
    deliverPlannedEvent(sink, plan, index, checkpoint);
  }
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect({path.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence), persistence.message);
  const auto &restoredSchedule =
      inspection.observerRecord.outputFiles[0].groups[0].schedule;
  const auto *next =
      inspection.scheduleContinuations[0].cursor.values.value("nextOrdinal");
  require(restoredSchedule.typeIdentifier == quadraticScheduleType &&
              restoredSchedule.contractVersion == 1 &&
              restoredSchedule.configuration.value("scale") != nullptr &&
              inspection.scheduleContinuations[0].cursor.committedOrdinal ==
                  1 &&
              next != nullptr &&
              std::get<std::vector<std::int64_t>>(next->storage)[0] == 2,
          "algorithmic schedule configuration and cursor round trip");

  WVOutputPlan resumed;
  status = WVOutputPlan::create(
      descriptor, modelOutputCatalog(), initialTime + 1.0, initialTime + 4.0,
      inspection.scheduleContinuations, resumed);
  require(static_cast<bool>(status) && resumed.eventCount() == 1 &&
              resumed.event(0).scheduledTime == initialTime + 4.0,
          "algorithmic schedule append cursor");

  const auto beforeRejectedAppend = fileBytes(path);
  auto mismatchedProgress = inspection.destinationProgress;
  auto &mismatchedCursorValues =
      mismatchedProgress[0].committedScheduleCursor.values.values;
  const auto mismatchedNext = std::find_if(
      mismatchedCursorValues.begin(), mismatchedCursorValues.end(),
      [](const auto &value) { return value.name == "nextOrdinal"; });
  require(mismatchedNext != mismatchedCursorValues.end(),
          "algorithmic destination cursor payload is missing");
  ++std::get<std::vector<std::int64_t>>(mismatchedNext->storage)[0];
  WVModelOutputNetCDFSink rejectedAppend;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, resumed, layout, nullptr, mismatchedProgress,
      rejectedAppend);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::appendConflict &&
              fileBytes(path) == beforeRejectedAppend,
          "append accepted typed destination-cursor drift or mutated output");
  const auto rejectProgress = [&](WVOutputDestinationProgress progress,
                                  const std::string &description) {
    WVModelOutputNetCDFSink rejected;
    const auto rejectedStatus = WVModelOutputNetCDFSink::openAppend(
        configuration, descriptor, resumed, layout, nullptr,
        {std::move(progress)}, rejected);
    require(!rejectedStatus &&
                rejectedStatus.code == WVCheckpointStatusCode::appendConflict &&
                fileBytes(path) == beforeRejectedAppend,
            "append accepted " + description + " or mutated output");
  };
  auto recordCountMismatch = inspection.destinationProgress[0];
  ++recordCountMismatch.recordCount;
  rejectProgress(std::move(recordCountMismatch),
                 "destination record-count drift");
  auto markerMismatch = inspection.destinationProgress[0];
  markerMismatch.hasCommittedTime = false;
  rejectProgress(std::move(markerMismatch),
                 "destination time-marker drift");
  auto lastTimeMismatch = inspection.destinationProgress[0];
  lastTimeMismatch.lastCommittedTime += 1.0;
  rejectProgress(std::move(lastTimeMismatch),
                 "destination time-last drift");

  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, resumed, layout, nullptr,
      inspection.destinationProgress, append);
  require(static_cast<bool>(persistence), persistence.message);
  status = append.preflight(resumed);
  require(static_cast<bool>(status), status.message);
  checkpoint.state.t = resumed.event(0).scheduledTime;
  deliverPlannedEvent(append, resumed, 0, checkpoint);
  persistence = append.close();
  require(static_cast<bool>(persistence), persistence.message);
}

void testAlgorithmicScheduleThroughModelAndRequestRunner() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  constexpr double scale = 1e-7;
  const auto scheduledTime = [&](WVOutputScheduleOrdinal ordinal) {
    const auto n = static_cast<double>(ordinal);
    return std::fma(scale, n * n, initialTime);
  };
  const auto sourcePath = directory.path / "algorithmic-model-source.nc";
  auto record = recordFor(checkpoint, sourcePath);
  record.outputFiles[0].groups[0].schedule =
      quadraticSchedule(scheduledTime(5), initialTime, scale);
  record.outputFiles[0].groups.push_back(
      {"evenly-spaced-diagnostics", "evenly-spaced-diagnostics",
       {scale, initialTime, scheduledTime(5)}, {}, false});
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan seedPlan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime,
                                scheduledTime(1), {}, seedPlan);
  require(static_cast<bool>(status) && seedPlan.eventCount() == 2,
          "algorithmic WVModel seed plan mismatch");
  WVModelOutputNetCDFSink seed;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, descriptor, seedPlan, layout,
      nullptr, seed);
  require(static_cast<bool>(persistence), persistence.message);
  status = seed.preflight(seedPlan);
  require(static_cast<bool>(status), status.message);
  for (std::size_t eventIndex = 0; eventIndex < seedPlan.eventCount();
       ++eventIndex) {
    const auto planned = seedPlan.event(eventIndex);
    checkpoint.state.t = planned.scheduledTime;
    WVOutputEvent event;
    event.eventOrdinal = planned.eventOrdinal;
    event.scheduledTime = planned.scheduledTime;
    event.state = eventState(checkpoint, planned.scheduledTime);
    event.routes = planned.routes;
    event.routeCount = planned.routeCount;
    for (std::size_t routeIndex = 0; routeIndex < planned.routeCount;
         ++routeIndex) {
      WVOutputDeliveryResult delivery;
      status = seed.deliver(event, planned.routes[routeIndex], delivery);
      require(static_cast<bool>(status), status.message);
    }
  }
  persistence = seed.close();
  require(static_cast<bool>(persistence), persistence.message);
  const auto immutableSource = fileBytes(sourcePath);

  WVModelOutputNetCDFInspection sourceInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {sourcePath.string()}, *modelOutputCatalog(), sourceInspection);
  require(static_cast<bool>(persistence), persistence.message);
  InspectionFactoryCounts rawInspectionCounts;
  const auto rawInspectionCatalog =
      inspectionTrapCatalog(rawInspectionCounts);
  WVModelOutputNetCDFInspection trappedAlgorithmicInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {sourcePath.string()}, *rawInspectionCatalog,
      trappedAlgorithmicInspection);
  activeInspectionFactoryCounts = nullptr;
  require(static_cast<bool>(persistence) &&
              rawInspectionCounts.observers == 0 &&
              rawInspectionCounts.schedules == 0 &&
              rawInspectionCounts.forcings == 0,
          "algorithmic raw inspection constructed a runtime provider");
  WVModelOutputRequest rejectedRequest;
  rejectedRequest.policy = WVModelOutputPolicy::create;
  rejectedRequest.destinations = {
      {"primary", (directory.path / "rejected-output.nc").string()}};
  rejectedRequest.finalTime = scheduledTime(2);
  WVModelOutputConfiguration rejectedConfiguration;

  std::shared_ptr<const WVExtensionCatalog> builtInCatalog;
  status = makeBuiltInExtensionCatalog(builtInCatalog);
  require(static_cast<bool>(status), status.message);
  status = WVModel::prepareModelOutput(
      builtInCatalog, sourceInspection, rejectedRequest,
      rejectedConfiguration);
  require(!status &&
              !std::filesystem::exists(directory.path / "rejected-output.nc"),
          "missing algorithmic factory passed preflight or created output");

  InspectionFactoryCounts factoryCounts;
  const auto trapCatalog = inspectionTrapCatalog(factoryCounts);
  status = WVModel::prepareModelOutput(
      trapCatalog, sourceInspection, rejectedRequest, rejectedConfiguration);
  activeInspectionFactoryCounts = nullptr;
  require(!status && factoryCounts.schedules == 1 &&
              factoryCounts.observers == 0 && factoryCounts.forcings == 0,
          "schedule preflight constructed observers or forcing before the "
          "schedule failure");

  auto incompatibleVersion = sourceInspection;
  incompatibleVersion.observerRecord.outputFiles[0]
      .groups[0]
      .schedule.contractVersion = 99;
  status = WVModel::prepareModelOutput(
      modelOutputCatalog(), incompatibleVersion, rejectedRequest,
      rejectedConfiguration);
  require(!status,
          "an incompatible algorithmic schedule version passed preflight");

  auto incompatibleConfiguration = sourceInspection;
  auto &scheduleValues = incompatibleConfiguration.observerRecord.outputFiles[0]
                             .groups[0]
                             .schedule.configuration.values;
  const auto scaleValue = std::find_if(
      scheduleValues.begin(), scheduleValues.end(),
      [](const auto &value) { return value.name == "scale"; });
  require(scaleValue != scheduleValues.end(),
          "algorithmic schedule scale configuration is missing");
  std::get<std::vector<double>>(scaleValue->storage)[0] = -scale;
  status = WVModel::prepareModelOutput(
      modelOutputCatalog(), incompatibleConfiguration, rejectedRequest,
      rejectedConfiguration);
  require(!status,
          "an incompatible algorithmic schedule configuration passed "
          "preflight");

  auto incompatibleCursor = sourceInspection;
  auto &cursorValues =
      incompatibleCursor.scheduleContinuations[0].cursor.values.values;
  const auto nextOrdinal = std::find_if(
      cursorValues.begin(), cursorValues.end(),
      [](const auto &value) { return value.name == "nextOrdinal"; });
  require(nextOrdinal != cursorValues.end(),
          "algorithmic cursor payload is missing");
  ++std::get<std::vector<std::int64_t>>(nextOrdinal->storage)[0];
  status = WVModel::prepareModelOutput(
      modelOutputCatalog(), incompatibleCursor, rejectedRequest,
      rejectedConfiguration);
  require(!status,
          "an incompatible full typed schedule cursor passed preflight");

  const auto fixedPath = directory.path / "algorithmic-model-fixed.nc";
  WVModelOutputRequest fixedRequest;
  fixedRequest.policy = WVModelOutputPolicy::create;
  fixedRequest.destinations = {{"primary", fixedPath.string()}};
  fixedRequest.finalTime = scheduledTime(2);
  WVModel fixedModel;
  WVModelState fixedState;
  status = WVModel::createFromModelOutputFiles(
      modelOutputCatalog(), {sourcePath.string()}, fixedRequest,
      std::make_unique<WVReferenceFFTEngine>(), {}, fixedModel, fixedState);
  require(static_cast<bool>(status), status.message);
  status = fixedModel.prepareStateAfterRestart(fixedState);
  require(static_cast<bool>(status), status.message);
  status = fixedModel.advanceToTime(fixedState, fixedRequest.finalTime, scale);
  require(static_cast<bool>(status), status.message);
  const auto fixedMetrics = fixedModel.metrics(&fixedState);
  persistence = fixedModel.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);
  require(fixedMetrics.outputDriver.acceptedEndpointStateEventCount == 3 &&
              fixedMetrics.outputDriver.interpolatedStateEvaluationCount == 0,
          "fixed RK4 did not deliver the mixed-graph events at exact "
          "accepted endpoints");
  requireTimeSeries(fixedPath, {scheduledTime(2)});

  WVModelOutputNetCDFInspection appendGraphInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {fixedPath.string()}, *modelOutputCatalog(), appendGraphInspection);
  require(static_cast<bool>(persistence), persistence.message);
  const auto fixedBeforeRejectedAppend = fileBytes(fixedPath);
  appendGraphInspection.observerRecord.observers[0].outputScale += 1.0;
  WVModelOutputRequest rejectedAppendRequest;
  rejectedAppendRequest.policy = WVModelOutputPolicy::append;
  rejectedAppendRequest.finalTime = scheduledTime(3);
  status = WVModel::prepareModelOutput(
      modelOutputCatalog(), appendGraphInspection, rejectedAppendRequest,
      rejectedConfiguration);
  require(!status && fileBytes(fixedPath) == fixedBeforeRejectedAppend,
          "complete observer-configuration drift passed append preflight or "
          "mutated the destination");

  const auto adaptivePath = directory.path / "algorithmic-model-adaptive.nc";
  WVModelOutputRequest adaptiveRequest;
  adaptiveRequest.policy = WVModelOutputPolicy::create;
  adaptiveRequest.destinations = {{"primary", adaptivePath.string()}};
  adaptiveRequest.finalTime = scheduledTime(4);
  WVModelIntegratorConfiguration adaptiveConfiguration;
  adaptiveConfiguration.kind = WVModelIntegratorKind::adaptiveRK23;
  adaptiveConfiguration.adaptive.maximumStepSize = 15.0 * scale;
  WVModel adaptiveModel;
  WVModelState adaptiveState;
  status = WVModel::createFromModelOutputFiles(
      modelOutputCatalog(), {sourcePath.string()}, adaptiveRequest,
      std::make_unique<WVReferenceFFTEngine>(), adaptiveConfiguration,
      adaptiveModel, adaptiveState);
  require(static_cast<bool>(status), status.message);
  status = adaptiveModel.prepareStateAfterRestart(adaptiveState);
  require(static_cast<bool>(status), status.message);
  status = adaptiveModel.advanceToTime(adaptiveState,
                                       adaptiveRequest.finalTime,
                                       15.0 * scale);
  require(static_cast<bool>(status), status.message);
  const auto adaptiveMetrics = adaptiveModel.metrics(&adaptiveState);
  persistence = adaptiveModel.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);
  require(adaptiveMetrics.outputDriver.interpolatedStateEvaluationCount ==
                  14 &&
              adaptiveMetrics.outputDriver.acceptedEndpointStateEventCount ==
                  1,
          "adaptive RK3(2) did not deliver dense and exact mixed-graph "
          "events through WVModel");
  requireTimeSeries(adaptivePath,
                    {scheduledTime(2), scheduledTime(3), scheduledTime(4)});

  WVModelOutputRequest appendRequest;
  appendRequest.policy = WVModelOutputPolicy::append;
  appendRequest.finalTime = scheduledTime(5);
  WVModel appendModel;
  WVModelState appendState;
  status = WVModel::createFromModelOutputFiles(
      modelOutputCatalog(), {adaptivePath.string()}, appendRequest,
      std::make_unique<WVReferenceFFTEngine>(), {}, appendModel, appendState);
  require(static_cast<bool>(status), status.message);
  status = appendModel.prepareStateAfterRestart(appendState);
  require(static_cast<bool>(status), status.message);
  status = appendModel.advanceToTime(appendState, appendRequest.finalTime,
                                     scale);
  require(static_cast<bool>(status), status.message);
  persistence = appendModel.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);
  requireTimeSeries(adaptivePath,
                    {scheduledTime(2), scheduledTime(3), scheduledTime(4),
                     scheduledTime(5)});

  const auto replacePath = directory.path / "algorithmic-model-replace.nc";
  std::filesystem::copy_file(sourcePath, replacePath);
  WVModelOutputRequest replaceRequest = fixedRequest;
  replaceRequest.policy = WVModelOutputPolicy::replace;
  replaceRequest.destinations = {{"primary", replacePath.string()}};
  WVModel replaceModel;
  WVModelState replaceState;
  status = WVModel::createFromModelOutputFiles(
      modelOutputCatalog(), {sourcePath.string()}, replaceRequest,
      std::make_unique<WVReferenceFFTEngine>(), {}, replaceModel,
      replaceState);
  require(static_cast<bool>(status), status.message);
  status = replaceModel.prepareStateAfterRestart(replaceState);
  require(static_cast<bool>(status), status.message);
  status = replaceModel.advanceToTime(replaceState, replaceRequest.finalTime,
                                      scale);
  require(static_cast<bool>(status), status.message);
  persistence = replaceModel.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);
  requireTimeSeries(replacePath, {scheduledTime(2)});

  const auto requestPath = directory.path / "algorithmic-run-request.json";
  const auto requestDestination =
      directory.path / "algorithmic-request-output.nc";
  const auto requestReport = directory.path / "algorithmic-request-report.json";
  std::ofstream request(requestPath, std::ios::binary | std::ios::trunc);
  request << std::setprecision(17)
          << "{\"schemaIdentifier\":\"wave-vortex-run-request-v1\","
             "\"schemaVersion\":1,\"modelFiles\":[\""
          << sourcePath.string()
          << "\"],\"integration\":{\"method\":\"fixed-rk4\","
             "\"finalTime\":"
          << scheduledTime(2) << ",\"initialStep\":" << scale
          << "},\"output\":{\"policy\":\"create\",\"destinations\":{"
             "\"primary\":\""
          << requestDestination.string()
          << "\"}},\"execution\":{\"fftProvider\":\"reference\","
             "\"threads\":1},\"report\":\""
          << requestReport.string() << "\"}";
  request.close();
  const auto command = shellQuote(WV_RUNTIME_EXTENDED_RUNNER) + " --request " +
                       shellQuote(requestPath) + " >/dev/null";
  require(std::system(command.c_str()) == 0,
          "explicit-catalog request runner rejected an algorithmic graph");
  requireTimeSeries(requestDestination, {scheduledTime(2)});
  require(std::filesystem::exists(requestReport),
          "algorithmic request runner omitted its report");
  require(fileBytes(sourcePath) == immutableSource,
          "create/replace/request continuations mutated their source model "
          "output");
}

void testVariableObservationBatches() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  const auto path = directory.path / "variable-observations.nc";
  WVObserverRecord synthetic;
  synthetic.identifier = "synthetic-observations";
  synthetic.name = "synthetic";
  synthetic.typeIdentifier = "WVTestObservationBatches";
  const auto configuredRecord = [&](const std::filesystem::path &destination) {
    auto configured = recordFor(checkpoint, destination);
    configured.observers.push_back(synthetic);
    configured.outputFiles[0].groups[0].observerIdentifiers.push_back(
        synthetic.identifier);
    return configured;
  };
  const auto invalidPath = directory.path / "invalid-observations.nc";
  auto invalidRecord = configuredRecord(invalidPath);
  auto invalidDescriptor = descriptorFor(invalidRecord);
  WVIntegrationStateLayout invalidLayout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, invalidDescriptor, invalidLayout);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFConfiguration configuration{modelOutputCatalog(), checkpoint, false};
  VariableBatchSource malformedSource(1, true);
  WVOutputPlan invalidPlan;
  status = WVOutputPlan::create(invalidDescriptor, modelOutputCatalog(), initialTime,
                                initialTime + 1.0, {}, invalidPlan);
  require(static_cast<bool>(status) && invalidPlan.eventCount() == 2,
          "malformed variable-batch output plan");
  WVModelOutputNetCDFSink malformedSink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      configuration, invalidDescriptor, invalidPlan, invalidLayout,
      &malformedSource, malformedSink);
  require(static_cast<bool>(persistence), persistence.message);
  status = malformedSink.preflight(invalidPlan);
  require(static_cast<bool>(status), status.message);

  WVOutputEvent invalidEvent;
  invalidEvent.eventOrdinal = invalidPlan.event(0).eventOrdinal;
  invalidEvent.scheduledTime = invalidPlan.event(0).scheduledTime;
  invalidEvent.state = eventState(checkpoint, invalidEvent.scheduledTime);
  invalidEvent.routes = invalidPlan.event(0).routes;
  invalidEvent.routeCount = invalidPlan.event(0).routeCount;
  WVOutputDeliveryResult delivery;
  status = malformedSink.deliver(invalidEvent, invalidEvent.routes[0], delivery);
  require(!status, "malformed ragged batch was written");
  delivery = {};
  status = malformedSink.deliver(invalidEvent, invalidEvent.routes[0], delivery);
  require(!status, "an exact-event retry accepted its malformed batch");
  require(malformedSource.prepareCount() == 1,
          "an exact-event retry repeated source preparation: " +
              std::to_string(malformedSource.prepareCount()));
  require(malformedSource.batchCount(synthetic.identifier) == 1,
          "an exact-event retry regenerated its malformed observation batch: " +
              std::to_string(
                  malformedSource.batchCount(synthetic.identifier)));
  int file = -1, group = -1, timeDimension = -1;
  std::size_t timeCount = 99;
  require(nc_open(invalidPath.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
              nc_inq_ncid(file, "wave-vortex", &group) == NC_NOERR &&
              nc_inq_dimid(group, "t", &timeDimension) == NC_NOERR &&
              nc_inq_dimlen(group, timeDimension, &timeCount) == NC_NOERR &&
              timeCount == 0 && nc_close(file) == NC_NOERR,
          "failed batch mutated the output before validation");
  persistence = malformedSink.close();
  require(static_cast<bool>(persistence), persistence.message);

  auto record = configuredRecord(path);
  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  registeredBatches->resetCounts();
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime, initialTime + 1.0, {},
                                plan);
  require(static_cast<bool>(status) && plan.eventCount() == 2,
          "variable-batch output plan");
  WVModelOutputNetCDFSink sink;
  persistence = WVModelOutputNetCDFSink::createNew(
      configuration, descriptor, plan, layout, source.get(), sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);

  WVOutputEvent event;
  event.eventOrdinal = plan.event(0).eventOrdinal;
  event.scheduledTime = plan.event(0).scheduledTime;
  event.state = eventState(checkpoint, event.scheduledTime);
  event.routes = plan.event(0).routes;
  event.routeCount = plan.event(0).routeCount;
  delivery = {};
  status = sink.deliver(event, event.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  event.eventOrdinal = plan.event(1).eventOrdinal;
  event.scheduledTime = plan.event(1).scheduledTime;
  event.state = eventState(checkpoint, event.scheduledTime);
  event.routes = plan.event(1).routes;
  event.routeCount = plan.event(1).routeCount;
  delivery = {};
  status = sink.deliver(event, event.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  require(sink.metrics().batchMaximumLiveBytes > 0 &&
              sink.metrics().batchRetainedStorageBytes > 0 &&
              sink.metrics().failureCount == 0 &&
              source->metrics().preparedEventCount == 3 &&
              registeredBatches->batchCount(synthetic.identifier) == 2,
          "variable-batch storage/failure metrics are incomplete");
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  file = -1;
  group = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
              nc_inq_ncid(file, "wave-vortex", &group) == NC_NOERR,
          "open variable observation output");
  int sampleDimension = -1, profileDimension = -1;
  std::size_t sampleCount = 0, profileCount = 0;
  require(nc_inq_dimid(group, "synthetic_sample", &sampleDimension) ==
                  NC_NOERR &&
              nc_inq_dimlen(group, sampleDimension, &sampleCount) == NC_NOERR &&
              sampleCount == 3 &&
              nc_inq_dimid(group, "synthetic_profile", &profileDimension) ==
                  NC_NOERR &&
              nc_inq_dimlen(group, profileDimension, &profileCount) ==
                  NC_NOERR &&
              profileCount == 2,
          "flat variable observation axes have the wrong committed extent");
  int progressVariable = -1;
  std::array<long long, 2> sampleProgress{};
  require(nc_inq_varid(group, "portableCommitted_synthetic_sample",
                       &progressVariable) == NC_NOERR &&
              nc_get_var_longlong(group, progressVariable,
                                  sampleProgress.data()) == NC_NOERR &&
              sampleProgress == std::array<long long, 2>{3, 3},
          "zero-length batch changed committed sample progress");
  int integerVariable = -1, booleanVariable = -1, textVariable = -1;
  nc_type integerType = NC_NAT, booleanType = NC_NAT, textType = NC_NAT;
  require(nc_inq_varid(group, "synthetic_pass", &integerVariable) == NC_NOERR &&
              nc_inq_vartype(group, integerVariable, &integerType) == NC_NOERR &&
              integerType == NC_INT64 &&
              nc_inq_varid(group, "synthetic_valid", &booleanVariable) ==
                  NC_NOERR &&
              nc_inq_vartype(group, booleanVariable, &booleanType) == NC_NOERR &&
              booleanType == NC_UBYTE &&
              nc_inq_varid(group, "synthetic_label", &textVariable) ==
                  NC_NOERR &&
              nc_inq_vartype(group, textVariable, &textType) == NC_NOERR &&
              textType == NC_STRING,
          "integer/Boolean/text observation types were not preserved");
  std::array<char *, 2> textPayload{};
  require(nc_get_var_string(group, textVariable, textPayload.data()) ==
                  NC_NOERR &&
              textPayload[0] != nullptr && *textPayload[0] == 0 &&
              textPayload[1] != nullptr &&
              std::string(textPayload[1]) == "pass-b" &&
              nc_free_string(textPayload.size(), textPayload.data()) ==
                  NC_NOERR,
          "empty text observation values were not preserved exactly");
  int metadataRoot = -1, metadata = -1;
  std::size_t schemaLength = 0;
  require(nc_inq_ncid(group, "observingSystems", &metadataRoot) == NC_NOERR &&
              nc_inq_ncid(metadataRoot, "observingSystems-2", &metadata) ==
                  NC_NOERR &&
              nc_inq_attlen(metadata, NC_GLOBAL,
                            "portableObservationSchemaIdentifier",
                            &schemaLength) == NC_NOERR &&
              schemaLength == std::string("synthetic-variable-observation").size() &&
              nc_close(file) == NC_NOERR,
          "observation schema identity was not persisted");

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect({path.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence) &&
              inspection.observerRecord.observers.size() == 2,
          "generic graph reader reconstructed " +
              std::to_string(inspection.observerRecord.observers.size()) +
              " observers: " + persistence.message);
  require(inspection.observationSchemas.size() == 1 &&
              inspection.observationSchemas[0].observerIdentifier ==
                  synthetic.identifier &&
              inspection.observationSchemas[0].schema.identifier ==
                  "synthetic-variable-observation" &&
              inspection.observationSchemas[0].schema.axes.size() == 3 &&
              inspection.observationSchemas[0].schema.variables.size() == 11,
          "generic graph reader did not reconstruct the portable schema");
  std::size_t prematureObserverConstructions = 0;
  WVExtensionCatalogBuilder driftBuilder;
  for (auto registration : modelOutputCatalog()->observers().registrations()) {
    const auto factory = registration.factory;
    registration.factory =
        [factory, &prematureObserverConstructions](
            const WVObserverRecord &observer,
            const WVPortableTypedRecord &configuration,
            std::shared_ptr<const WVObservingSystem> &result) {
          ++prematureObserverConstructions;
          return factory(observer, configuration, result);
        };
    if (registration.typeIdentifier == "WVTestObservationBatches") {
      const auto resolver = registration.outputPlanResolver;
      registration.outputPlanResolver =
          [resolver](const WVObserverRecord &observer,
                     const WVObserverOutputPlanningContext &context,
                     WVObserverOutputPlan &plan) {
            auto status = resolver(observer, context, plan);
            if (status)
              plan.schema.identifier += "-drift";
            return status;
          };
    }
    status = driftBuilder.addObserverFactory(std::move(registration));
    require(static_cast<bool>(status), status.message);
  }
  for (auto registration :
       modelOutputCatalog()->outputSchedules().registrations()) {
    status = driftBuilder.addOutputScheduleFactory(std::move(registration));
    require(static_cast<bool>(status), status.message);
  }
  for (auto registration : modelOutputCatalog()->forcings().registrations()) {
    status = driftBuilder.addForcingFactory(std::move(registration));
    require(static_cast<bool>(status), status.message);
  }
  std::shared_ptr<const WVExtensionCatalog> driftCatalog;
  status = driftBuilder.freeze(driftCatalog);
  require(static_cast<bool>(status), status.message);
  const auto driftDestination = directory.path / "schema-drift-output.nc";
  WVModelOutputRequest driftRequest;
  driftRequest.policy = WVModelOutputPolicy::create;
  driftRequest.destinations = {{"primary", driftDestination.string()}};
  driftRequest.finalTime = initialTime + 2.0;
  WVModelOutputConfiguration driftConfiguration;
  status = WVModel::prepareModelOutput(
      driftCatalog, inspection, driftRequest, driftConfiguration);
  require(!status && prematureObserverConstructions == 0 &&
              !std::filesystem::exists(driftDestination),
          "provider-schema drift was not rejected before observer "
          "construction or output mutation");
  const auto &inspectedSchema = inspection.observationSchemas[0].schema;
  const auto inspectedBins = std::find_if(
      inspectedSchema.variables.begin(), inspectedSchema.variables.end(),
      [](const auto &variable) { return variable.identifier == "bins"; });
  const auto inspectedRows = std::find_if(
      inspectedSchema.variables.begin(), inspectedSchema.variables.end(),
      [](const auto &variable) { return variable.identifier == "row-count"; });
  require(inspectedBins != inspectedSchema.variables.end() &&
              inspectedBins->scalarType ==
                  WVObservationScalarType::complex64 &&
              inspectedBins->dimensionIdentifiers ==
                  std::vector<std::string>{"depth", "sample"} &&
              inspectedRows != inspectedSchema.variables.end() &&
              inspectedRows->raggedRole ==
                  WVObservationRaggedRole::rowCount &&
              inspectedRows->raggedChildAxisIdentifier == "sample",
          "generic graph reader changed typed or ragged schema declarations");

  for (const auto *attribute : {"portableObserverContractVersion",
                                "portableObserverConfiguration", "name",
                                "portableIdentifier"}) {
    const auto malformedPath =
        directory.path / (std::string("missing-") + attribute + ".nc");
    std::filesystem::copy_file(path, malformedPath);
    int malformedFile = -1;
    int malformedGroup = -1;
    int metadataRoot = -1;
    int metadata = -1;
    require(nc_open(malformedPath.c_str(), NC_WRITE, &malformedFile) ==
                    NC_NOERR &&
                nc_inq_ncid(malformedFile, "wave-vortex",
                            &malformedGroup) == NC_NOERR &&
                nc_inq_ncid(malformedGroup, "observingSystems",
                            &metadataRoot) == NC_NOERR &&
                nc_inq_ncid(metadataRoot, "observingSystems-2", &metadata) ==
                    NC_NOERR &&
                nc_redef(malformedFile) == NC_NOERR &&
                nc_del_att(metadata, NC_GLOBAL, attribute) == NC_NOERR &&
                nc_enddef(malformedFile) == NC_NOERR &&
                nc_close(malformedFile) == NC_NOERR,
            std::string("failed to remove canonical observer attribute ") +
                attribute);
    WVModelOutputNetCDFInspection malformedInspection;
    const auto malformed = WVModelOutputNetCDFSink::inspect(
        {malformedPath.string()}, *modelOutputCatalog(),
        malformedInspection);
    require(!malformed,
            std::string("canonical observer record accepted missing ") +
                attribute);
  }

  const auto requireRejectedFill =
      [&](const std::string &label, const char *variableName, nc_type type) {
        const auto malformedPath = directory.path / (label + ".nc");
        std::filesystem::copy_file(path, malformedPath);
        int malformedFile = -1;
        int malformedGroup = -1;
        int malformedVariable = -1;
        const std::size_t position[] = {0};
        require(nc_open(malformedPath.c_str(), NC_WRITE, &malformedFile) ==
                        NC_NOERR &&
                    nc_inq_ncid(malformedFile, "wave-vortex",
                                &malformedGroup) == NC_NOERR &&
                    nc_inq_varid(malformedGroup, variableName,
                                 &malformedVariable) == NC_NOERR,
                "failed to open " + label + " fixture");
        int writeStatus = NC_EBADTYPE;
        if (type == NC_INT64) {
          const long long value = NC_FILL_INT64;
          writeStatus = nc_put_var1_longlong(
              malformedGroup, malformedVariable, position, &value);
        } else if (type == NC_UBYTE) {
          const unsigned char value = NC_FILL_UBYTE;
          writeStatus = nc_put_var1_uchar(
              malformedGroup, malformedVariable, position, &value);
        } else {
          const char *value = NC_FILL_STRING;
          writeStatus = nc_put_var1_string(
              malformedGroup, malformedVariable, position, &value);
        }
        require(writeStatus == NC_NOERR && nc_close(malformedFile) == NC_NOERR,
                "failed to write " + label + " fixture");
        WVModelOutputNetCDFInspection malformedInspection;
        const auto malformed = WVModelOutputNetCDFSink::inspect(
            {malformedPath.string()}, *modelOutputCatalog(),
            malformedInspection);
        require(!malformed &&
                    malformed.code == WVCheckpointStatusCode::incompleteRecord,
                "raw inspection accepted " + label);
      };
  requireRejectedFill("integer-fill", "synthetic_pass", NC_INT64);
  requireRejectedFill("boolean-fill", "synthetic_valid", NC_UBYTE);

  const auto raggedMismatchPath = directory.path / "ragged-mismatch.nc";
  std::filesystem::copy_file(path, raggedMismatchPath);
  int raggedFile = -1;
  int raggedGroup = -1;
  int raggedProgress = -1;
  const std::size_t raggedRecordIndex = 1;
  const long long inconsistentCommittedOffset = 2;
  require(nc_open(raggedMismatchPath.c_str(), NC_WRITE, &raggedFile) ==
                  NC_NOERR &&
              nc_inq_ncid(raggedFile, "wave-vortex", &raggedGroup) ==
                  NC_NOERR &&
              nc_inq_varid(raggedGroup,
                           "portableCommitted_synthetic_sample",
                           &raggedProgress) == NC_NOERR &&
              nc_put_var1_longlong(raggedGroup, raggedProgress,
                                   &raggedRecordIndex,
                                   &inconsistentCommittedOffset) == NC_NOERR &&
              nc_close(raggedFile) == NC_NOERR,
          "failed to create ragged-progress mismatch fixture");
  WVModelOutputNetCDFInspection raggedMismatchInspection;
  const auto raggedMismatch = WVModelOutputNetCDFSink::inspect(
      {raggedMismatchPath.string()}, *modelOutputCatalog(),
      raggedMismatchInspection);
  require(!raggedMismatch,
          "raw graph inspection accepted committed/physical ragged-offset "
          "drift");

  const auto requireRejectedTextAttribute =
      [&](const std::string &label, const std::string &variableName,
          const char *attributeName, const std::string &value) {
        const auto malformedPath = directory.path / (label + ".nc");
        std::filesystem::copy_file(path, malformedPath);
        int malformedFile = -1;
        int malformedGroup = -1;
        int malformedVariable = -1;
        require(nc_open(malformedPath.c_str(), NC_WRITE, &malformedFile) ==
                    NC_NOERR &&
                    nc_inq_ncid(malformedFile, "wave-vortex",
                                &malformedGroup) == NC_NOERR &&
                    nc_inq_varid(malformedGroup, variableName.c_str(),
                                 &malformedVariable) == NC_NOERR &&
                    nc_redef(malformedFile) == NC_NOERR &&
                    nc_put_att_text(malformedGroup, malformedVariable,
                                    attributeName, value.size(),
                                    value.c_str()) == NC_NOERR &&
                    nc_enddef(malformedFile) == NC_NOERR &&
                    nc_close(malformedFile) == NC_NOERR,
                "failed to create malformed " + label + " fixture");
        WVModelOutputNetCDFInspection malformedInspection;
        const auto malformedStatus = WVModelOutputNetCDFSink::inspect(
            {malformedPath.string()}, *modelOutputCatalog(), malformedInspection);
        require(!malformedStatus,
                "reader accepted malformed " + label + " metadata");
      };
  requireRejectedTextAttribute(
      "unknown-layout", "synthetic_pass",
      "portableObservationValueLayout", "unknown");
  requireRejectedTextAttribute(
      "unknown-coordinate", "synthetic_time", "portableCoordinateRole",
      "unknown");
  requireRejectedTextAttribute(
      "unknown-ragged", "synthetic_row_count", "portableRaggedRole",
      "unknown");
  requireRejectedTextAttribute("custom-attribute-drift",
                               "synthetic_bins_real", "provider", "changed");
  const auto axisRolePath = directory.path / "unknown-axis-role.nc";
  std::filesystem::copy_file(path, axisRolePath);
  int axisFile = -1;
  int axisGroup = -1;
  int axisVariable = -1;
  const char *unknownRole = "unknown";
  require(nc_open(axisRolePath.c_str(), NC_WRITE, &axisFile) == NC_NOERR &&
              nc_inq_ncid(axisFile, "wave-vortex", &axisGroup) == NC_NOERR &&
              nc_inq_varid(axisGroup, "synthetic_time", &axisVariable) ==
                  NC_NOERR &&
              nc_redef(axisFile) == NC_NOERR &&
              nc_put_att_string(axisGroup, axisVariable,
                                "portableObservationAxisCoordinateRoles", 1,
                                &unknownRole) == NC_NOERR &&
              nc_enddef(axisFile) == NC_NOERR &&
              nc_close(axisFile) == NC_NOERR,
          "failed to create malformed axis-role fixture");
  WVModelOutputNetCDFInspection malformedAxisInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {axisRolePath.string()}, *modelOutputCatalog(), malformedAxisInspection);
  require(!persistence, "reader accepted an unknown axis coordinate role");

  WVOutputPlan appendPlan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime + 1.0,
                                initialTime + 2.0,
                                inspection.scheduleContinuations,
                                appendPlan);
  require(static_cast<bool>(status) && appendPlan.eventCount() == 1,
          "variable-batch append plan");
  VariableBatchSource driftedSource(2);
  WVModelOutputNetCDFSink drifted;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, appendPlan, layout, &driftedSource,
      inspection.destinationProgress, drifted);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::appendConflict,
          "append accepted observation schema drift");

  std::unique_ptr<WVObserverOutputEvaluationService> appendSource;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), appendSource);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      configuration, descriptor, appendPlan, layout, appendSource.get(),
      inspection.destinationProgress, append);
  require(static_cast<bool>(persistence), persistence.message);
  status = append.preflight(appendPlan);
  require(static_cast<bool>(status), status.message);
  event.eventOrdinal = appendPlan.event(0).eventOrdinal;
  event.scheduledTime = appendPlan.event(0).scheduledTime;
  event.state = eventState(checkpoint, event.scheduledTime);
  event.routes = appendPlan.event(0).routes;
  event.routeCount = appendPlan.event(0).routeCount;
  delivery = {};
  status = append.deliver(event, event.routes[0], delivery);
  require(static_cast<bool>(status), status.message);
  persistence = append.close();
  require(static_cast<bool>(persistence), persistence.message);
}

void testRegisteredTopologyProviders() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  const auto path = directory.path / "registered-topologies.nc";
  auto record = recordFor(checkpoint, path);
  const auto addObserver = [&](std::string identifier, std::string name,
                               std::string type) {
    WVObserverRecord observer;
    observer.identifier = std::move(identifier);
    observer.name = std::move(name);
    observer.typeIdentifier = std::move(type);
    record.observers.push_back(std::move(observer));
  };
  addObserver("fixed-bin", "fixed-bin", "WVTestFixedBin");
  addObserver("moving-track", "moving-track", "WVTestMovingTrack");
  addObserver("variable-pass", "variable-pass", "WVTestVariablePass");
  addObserver("ragged-a", "ragged-a", "WVTestRaggedProfile");
  addObserver("ragged-b", "ragged-b", "WVTestRaggedProfile");
  addObserver("nested-ragged", "nested-ragged", "WVTestNestedRagged");
  auto &primary = record.outputFiles[0].groups[0];
  primary.observerIdentifiers.insert(primary.observerIdentifiers.end(),
                                     {"fixed-bin", "moving-track",
                                      "variable-pass", "ragged-a",
                                      "nested-ragged"});
  auto secondary = primary;
  secondary.identifier = "ragged-secondary";
  secondary.name = "ragged-secondary";
  secondary.containsCompleteCoefficientRestart = false;
  secondary.observerIdentifiers = {"ragged-b"};
  record.outputFiles[0].groups.push_back(std::move(secondary));

  const auto requireNestedPersistence =
      [&](std::size_t occurrenceCount) {
        int file = -1;
        int group = -1;
        int passDimension = -1;
        int profileDimension = -1;
        int sampleDimension = -1;
        std::size_t passCount = 0;
        std::size_t profileCount = 0;
        std::size_t sampleCount = 0;
        int profileCountVariable = -1;
        int sampleOffsetVariable = -1;
        std::vector<long long> profileCounts(2 * occurrenceCount);
        std::vector<long long> sampleOffsets(3 * occurrenceCount);
        require(
            nc_open(path.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
                nc_inq_ncid(file, "wave-vortex", &group) == NC_NOERR &&
                nc_inq_dimid(group, "nested_pass", &passDimension) ==
                    NC_NOERR &&
                nc_inq_dimlen(group, passDimension, &passCount) == NC_NOERR &&
                nc_inq_dimid(group, "nested_profile", &profileDimension) ==
                    NC_NOERR &&
                nc_inq_dimlen(group, profileDimension, &profileCount) ==
                    NC_NOERR &&
                nc_inq_dimid(group, "nested_sample", &sampleDimension) ==
                    NC_NOERR &&
                nc_inq_dimlen(group, sampleDimension, &sampleCount) ==
                    NC_NOERR &&
                passCount == 2 * occurrenceCount &&
                profileCount == 3 * occurrenceCount &&
                sampleCount == 4 * occurrenceCount &&
                nc_inq_varid(group, "nested_pass_profile_count",
                             &profileCountVariable) == NC_NOERR &&
                nc_get_var_longlong(group, profileCountVariable,
                                    profileCounts.data()) == NC_NOERR &&
                nc_inq_varid(group, "nested_profile_sample_offset",
                             &sampleOffsetVariable) == NC_NOERR &&
                nc_get_var_longlong(group, sampleOffsetVariable,
                                    sampleOffsets.data()) == NC_NOERR &&
                nc_close(file) == NC_NOERR,
            "nested pass/profile/sample persistence is incomplete");
        for (std::size_t occurrence = 0; occurrence < occurrenceCount;
             ++occurrence) {
          require(profileCounts[2 * occurrence] == 2 &&
                      profileCounts[2 * occurrence + 1] == 1 &&
                      sampleOffsets[3 * occurrence] == 0 &&
                      sampleOffsets[3 * occurrence + 1] == 2 &&
                      sampleOffsets[3 * occurrence + 2] == 2,
                  "nested ragged counts or occurrence-local offsets changed "
                  "during persistence");
        }
      };

  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime, initialTime + 1.0, {},
                                plan);
  require(static_cast<bool>(status) && plan.groupCount() == 2,
          "registered topology plan did not retain both groups");
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, descriptor, plan, layout,
      source.get(), sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  for (std::size_t eventIndex = 0; eventIndex < plan.eventCount(); ++eventIndex) {
    const auto planned = plan.event(eventIndex);
    WVOutputEvent event;
    event.eventOrdinal = planned.eventOrdinal;
    event.scheduledTime = planned.scheduledTime;
    event.state = eventState(checkpoint, planned.scheduledTime);
    event.routes = planned.routes;
    event.routeCount = planned.routeCount;
    for (std::size_t route = 0; route < planned.routeCount; ++route) {
      WVOutputDeliveryResult delivery;
      status = sink.deliver(event, planned.routes[route], delivery);
      require(static_cast<bool>(status), status.message);
    }
  }
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);
  requireNestedPersistence(2);

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {path.string()}, *modelOutputCatalog(), inspection);
  require(static_cast<bool>(persistence),
          persistence.message + " at " + persistence.location);
  const auto raggedSchemas = std::count_if(
      inspection.observationSchemas.begin(),
      inspection.observationSchemas.end(), [](const auto &candidate) {
        return candidate.schema.identifier == "test-ragged-profile-v1" &&
               (candidate.observerIdentifier == "ragged-a" ||
                candidate.observerIdentifier == "ragged-b");
      });
  const auto fixedSchema = std::find_if(
      inspection.observationSchemas.begin(),
      inspection.observationSchemas.end(), [](const auto &candidate) {
        return candidate.observerIdentifier == "fixed-bin";
      });
  const auto nestedSchema = std::find_if(
      inspection.observationSchemas.begin(),
      inspection.observationSchemas.end(), [](const auto &candidate) {
        return candidate.observerIdentifier == "nested-ragged" &&
               candidate.schema.identifier == "test-nested-ragged-v1";
      });
  require(raggedSchemas == 2 &&
              inspection.observationSchemas.size() == 6 &&
              fixedSchema != inspection.observationSchemas.end() &&
              nestedSchema != inspection.observationSchemas.end() &&
              nestedSchema->schema.variables.size() == 3 &&
              nestedSchema->schema.variables[0]
                      .raggedChildAxisIdentifier == "profile" &&
              nestedSchema->schema.variables[1]
                      .raggedChildAxisIdentifier == "sample",
          "registered providers or same-provider ownership did not round trip");
  const auto fixedValue = std::find_if(
      fixedSchema->schema.variables.begin(), fixedSchema->schema.variables.end(),
      [](const auto &variable) { return variable.identifier == "value"; });
  require(fixedValue != fixedSchema->schema.variables.end() &&
              fixedValue->attributes.size() == 2 &&
              fixedValue->attributes[0].name == "ordered-a" &&
              fixedValue->attributes[1].name == "ordered-b",
          "ordered schema attributes did not round trip");

  std::unique_ptr<WVObserverOutputEvaluationService> appendSource;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), appendSource);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan appendPlan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime + 1.0,
                                initialTime + 2.0,
                                inspection.scheduleContinuations,
                                appendPlan);
  require(static_cast<bool>(status) && appendPlan.eventCount() == 1,
          "same-provider append plan is incorrect");
  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      {modelOutputCatalog(), checkpoint, false}, descriptor, appendPlan, layout,
      appendSource.get(), inspection.destinationProgress, append);
  require(static_cast<bool>(persistence), persistence.message);
  status = append.preflight(appendPlan);
  require(static_cast<bool>(status), status.message);
  const auto planned = appendPlan.event(0);
  WVOutputEvent appendEvent;
  appendEvent.eventOrdinal = planned.eventOrdinal;
  appendEvent.scheduledTime = planned.scheduledTime;
  appendEvent.state = eventState(checkpoint, planned.scheduledTime);
  appendEvent.routes = planned.routes;
  appendEvent.routeCount = planned.routeCount;
  for (std::size_t route = 0; route < planned.routeCount; ++route) {
    WVOutputDeliveryResult delivery;
    status = append.deliver(appendEvent, planned.routes[route], delivery);
    require(static_cast<bool>(status), status.message);
  }
  persistence = append.close();
  require(static_cast<bool>(persistence), persistence.message);
  requireNestedPersistence(3);

  const auto invalidPath = directory.path / "registered-ragged-retry.nc";
  auto invalidRecord = recordFor(checkpoint, invalidPath);
  WVObserverRecord malformed;
  malformed.identifier = "malformed-ragged";
  malformed.name = "malformed";
  malformed.typeIdentifier = "WVTestRaggedProfile";
  invalidRecord.observers.push_back(malformed);
  invalidRecord.outputFiles[0].groups[0].observerIdentifiers.push_back(
      malformed.identifier);
  auto invalidDescriptor = descriptorFor(invalidRecord);
  WVIntegrationStateLayout invalidLayout;
  status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, invalidDescriptor, invalidLayout);
  require(static_cast<bool>(status), status.message);
  raggedProfileProvider->resetCounts();
  std::unique_ptr<WVObserverOutputEvaluationService> invalidSource;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, invalidDescriptor,
      std::make_unique<WVReferenceFFTEngine>(), invalidSource);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan invalidPlan;
  status = WVOutputPlan::create(invalidDescriptor, modelOutputCatalog(), initialTime, initialTime, {},
                                invalidPlan);
  require(static_cast<bool>(status), status.message);
  WVModelOutputNetCDFSink invalidSink;
  persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, invalidDescriptor,
      invalidPlan, invalidLayout, invalidSource.get(), invalidSink);
  require(static_cast<bool>(persistence), persistence.message);
  status = invalidSink.preflight(invalidPlan);
  require(static_cast<bool>(status), status.message);
  const auto invalidPlanned = invalidPlan.event(0);
  WVOutputEvent invalidEvent;
  invalidEvent.eventOrdinal = invalidPlanned.eventOrdinal;
  invalidEvent.scheduledTime = invalidPlanned.scheduledTime;
  invalidEvent.state = eventState(checkpoint, invalidPlanned.scheduledTime);
  invalidEvent.routes = invalidPlanned.routes;
  invalidEvent.routeCount = invalidPlanned.routeCount;
  WVOutputDeliveryResult delivery;
  status = invalidSink.deliver(invalidEvent, invalidPlanned.routes[0], delivery);
  require(!status, "malformed registered ragged batch was accepted");
  delivery = {};
  status = invalidSink.deliver(invalidEvent, invalidPlanned.routes[0], delivery);
  require(!status &&
              raggedProfileProvider->batchCount(malformed.identifier) == 1,
          "registered exact-event retry regenerated its retained batch");
  persistence = invalidSink.close();
  require(static_cast<bool>(persistence), persistence.message);
}

void testObservationGraphCollisionPreflight() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto requireRejected = [&](bool conflictingAxis,
                                   const std::string &filename) {
    const auto path = directory.path / filename;
    auto record = recordFor(checkpoint, path);
    for (const auto &identifier : {"conflict-a", "conflict-b"}) {
      WVObserverRecord observer;
      observer.identifier = identifier;
      observer.name = identifier;
      observer.typeIdentifier = "WVTestObservationBatches";
      record.observers.push_back(observer);
      record.outputFiles[0].groups[0].observerIdentifiers.push_back(identifier);
    }
    auto descriptor = descriptorFor(record);
    WVIntegrationStateLayout layout;
    auto status = WVIntegrationStateLayout::create(
        checkpoint.state.coefficients.shape, descriptor, layout);
    require(static_cast<bool>(status), status.message);
    WVOutputPlan plan;
    status = WVOutputPlan::create(
        descriptor, modelOutputCatalog(), checkpoint.state.t,
        checkpoint.state.t + 2.0, {}, plan);
    require(static_cast<bool>(status), status.message);
    ConflictingSchemaSource source(conflictingAxis);
    WVModelOutputNetCDFSink sink;
    const auto persistence = WVModelOutputNetCDFSink::createNew(
        {modelOutputCatalog(), checkpoint, false}, descriptor, plan, layout,
        &source, sink);
    require(!persistence && !std::filesystem::exists(path),
            conflictingAxis
                ? "incompatible shared axes survived complete-graph preflight"
                : "duplicate persisted variable names survived complete-graph preflight");
  };
  requireRejected(false, "variable-collision.nc");
  requireRejected(true, "axis-collision.nc");
}

void testCoincidentRoutesKeepDistinctLogicalBatches() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  const auto primaryPath = directory.path / "coincident-primary.nc";
  const auto secondaryPath = directory.path / "coincident-secondary.nc";
  auto record = recordFor(checkpoint, primaryPath);
  WVObserverRecord synthetic;
  synthetic.identifier = "coincident-observations";
  synthetic.name = "coincident";
  synthetic.typeIdentifier = "WVTestObservationBatches";
  record.observers.push_back(synthetic);
  record.outputFiles[0].groups[0].observerIdentifiers.push_back(
      synthetic.identifier);
  auto secondary = record.outputFiles[0];
  secondary.identifier = "secondary";
  secondary.destination = secondaryPath.string();
  secondary.groups[0].identifier = "secondary-restart";
  secondary.groups[0].name = "wave-vortex-secondary";
  record.outputFiles.push_back(std::move(secondary));

  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, descriptor, layout);
  require(static_cast<bool>(status), status.message);
  registeredBatches->resetCounts();
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      checkpoint.configuration, false, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), source);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, modelOutputCatalog(), initialTime, initialTime, {}, plan);
  require(static_cast<bool>(status) && plan.eventCount() == 1 &&
              plan.event(0).routeCount == 2,
          "coincident-route output plan did not preserve two routes");
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, descriptor, plan, layout,
      source.get(), sink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);

  WVOutputEvent event;
  event.eventOrdinal = plan.event(0).eventOrdinal;
  event.scheduledTime = plan.event(0).scheduledTime;
  event.state = eventState(checkpoint, event.scheduledTime);
  event.routes = plan.event(0).routes;
  event.routeCount = plan.event(0).routeCount;
  WVOutputDeliveryResult delivery;
  status = sink.deliver(event, event.routes[0], delivery);
  require(static_cast<bool>(status) &&
              source->metrics().preparedEventCount == 2 &&
              registeredBatches->batchCount(synthetic.identifier) == 2,
          "coincident distinct logical groups shared an observation batch");
  delivery = {};
  status = sink.deliver(event, event.routes[1], delivery);
  require(static_cast<bool>(status) &&
              source->metrics().preparedEventCount == 2 &&
              registeredBatches->batchCount(synthetic.identifier) == 2,
          "second coincident route regenerated an exact-event batch");
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  WVObserverOutputPlan driftedPlan;
  status = registeredBatches->outputPlan(
      synthetic, WVObserverOutputPlanningContext{}, driftedPlan);
  require(static_cast<bool>(status), status.message);
  driftedPlan.schema.metadata.attributes.push_back(
      {"schema-drift-probe", "secondary-only"});
  std::vector<std::uint8_t> driftedManifest;
  status = encodeObservationSchemaManifest(driftedPlan.schema,
                                           driftedManifest);
  require(static_cast<bool>(status), status.message);
  const auto hex = [](const std::vector<std::uint8_t> &bytes) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string result;
    result.reserve(2 * bytes.size());
    for (const auto byte : bytes) {
      result.push_back(digits[byte >> 4]);
      result.push_back(digits[byte & 0x0f]);
    }
    return result;
  };
  const auto driftedManifestText = hex(driftedManifest);
  int driftedFile = -1;
  int driftedGroup = -1;
  int driftedMetadataRoot = -1;
  int driftedMetadata = -1;
  require(nc_open(secondaryPath.c_str(), NC_WRITE, &driftedFile) ==
                  NC_NOERR &&
              nc_inq_ncid(driftedFile, "wave-vortex-secondary",
                          &driftedGroup) == NC_NOERR &&
              nc_inq_ncid(driftedGroup, "observingSystems",
                          &driftedMetadataRoot) == NC_NOERR,
          "failed to open secondary observer metadata");
  int metadataCount = 0;
  require(nc_inq_grps(driftedMetadataRoot, &metadataCount, nullptr) ==
              NC_NOERR,
          "failed to enumerate secondary observer metadata");
  std::vector<int> metadataGroups(static_cast<std::size_t>(metadataCount));
  require(nc_inq_grps(driftedMetadataRoot, &metadataCount,
                      metadataGroups.data()) == NC_NOERR,
          "failed to read secondary observer metadata groups");
  for (const int candidate : metadataGroups) {
    std::size_t identifierLength = 0;
    if (nc_inq_attlen(candidate, NC_GLOBAL, "portableIdentifier",
                      &identifierLength) != NC_NOERR)
      continue;
    std::string identifier(identifierLength, '\0');
    if (nc_get_att_text(candidate, NC_GLOBAL, "portableIdentifier",
                        identifier.data()) == NC_NOERR &&
        identifier == synthetic.identifier) {
      driftedMetadata = candidate;
      break;
    }
  }
  require(driftedMetadata >= 0 && nc_redef(driftedFile) == NC_NOERR &&
          nc_put_att_text(driftedMetadata, NC_GLOBAL,
                          "portableObservationSchemaManifest",
                          driftedManifestText.size(),
                          driftedManifestText.c_str()) == NC_NOERR &&
          nc_enddef(driftedFile) == NC_NOERR &&
          nc_close(driftedFile) == NC_NOERR,
      "failed to create metadata-only shared-schema drift");
  WVModelOutputNetCDFInspection primaryInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {primaryPath.string()}, *modelOutputCatalog(), primaryInspection);
  require(static_cast<bool>(persistence), persistence.message);
  WVModelOutputNetCDFInspection secondaryInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {secondaryPath.string()}, *modelOutputCatalog(), secondaryInspection);
  require(static_cast<bool>(persistence), persistence.message);
  const auto schemaFor = [&](const WVModelOutputNetCDFInspection &candidate)
      -> const WVObservationSchema * {
    const auto found = std::find_if(
        candidate.observationSchemas.begin(),
        candidate.observationSchemas.end(), [&](const auto &schema) {
          return schema.observerIdentifier == synthetic.identifier;
        });
    return found == candidate.observationSchemas.end() ? nullptr
                                                        : &found->schema;
  };
  const auto *primarySchema = schemaFor(primaryInspection);
  const auto *secondarySchema = schemaFor(secondaryInspection);
  std::vector<std::uint8_t> primaryManifest;
  std::vector<std::uint8_t> secondaryManifest;
  require(primarySchema != nullptr && secondarySchema != nullptr &&
              encodeObservationSchemaManifest(*primarySchema,
                                              primaryManifest) &&
              encodeObservationSchemaManifest(*secondarySchema,
                                              secondaryManifest) &&
              primaryManifest != secondaryManifest,
          "metadata-only schema drift was not persisted in the secondary "
          "file");
  WVModelOutputNetCDFInspection driftedInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {primaryPath.string(), secondaryPath.string()}, *modelOutputCatalog(),
      driftedInspection);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::schemaMismatch,
          "raw inspection accepted metadata-only shared-schema drift");
}

void testWVModelRetainsFailedNetCDFRouteForRetry() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto initialTime = checkpoint.state.t;
  const auto sourcePath = directory.path / "retry-algorithmic-source.nc";
  auto sourceRecord = recordFor(checkpoint, sourcePath);
  sourceRecord.outputFiles[0].groups[0].schedule =
      quadraticSchedule(initialTime + 4.0, initialTime);
  auto sourceDescriptor = descriptorFor(sourceRecord);
  WVIntegrationStateLayout sourceLayout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, sourceDescriptor, sourceLayout);
  require(static_cast<bool>(status), status.message);
  WVOutputPlan sourcePlan;
  status = WVOutputPlan::create(sourceDescriptor, modelOutputCatalog(),
                                initialTime, initialTime + 1.0, {}, sourcePlan);
  require(static_cast<bool>(status) && sourcePlan.eventCount() == 2,
          "algorithmic retry source plan is incomplete");
  WVModelOutputNetCDFSink sourceSink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), checkpoint, false}, sourceDescriptor, sourcePlan,
      sourceLayout, nullptr, sourceSink);
  require(static_cast<bool>(persistence), persistence.message);
  status = sourceSink.preflight(sourcePlan);
  require(static_cast<bool>(status), status.message);
  for (std::size_t eventIndex = 0; eventIndex < sourcePlan.eventCount();
       ++eventIndex) {
    checkpoint.state.t = sourcePlan.event(eventIndex).scheduledTime;
    deliverPlannedEvent(sourceSink, sourcePlan, eventIndex, checkpoint);
  }
  persistence = sourceSink.close();
  require(static_cast<bool>(persistence), persistence.message);
  const auto immutableSource = fileBytes(sourcePath);
  WVModelOutputNetCDFInspection sourceInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {sourcePath.string()}, *modelOutputCatalog(), sourceInspection);
  require(static_cast<bool>(persistence), persistence.message);

  const auto finalTime = initialTime + 4.0;
  const auto primaryPath = directory.path / "retry-primary.nc";
  const auto secondaryPath = directory.path / "retry-secondary.nc";
  auto record = sourceInspection.observerRecord;
  record.outputFiles[0].destination = primaryPath.string();
  auto secondary = record.outputFiles[0];
  secondary.identifier = "secondary";
  secondary.destination = secondaryPath.string();
  secondary.groups[0].identifier = "secondary-restart";
  secondary.groups[0].name = "wave-vortex-secondary";
  record.outputFiles.push_back(std::move(secondary));

  auto descriptor = descriptorFor(record);
  WVIntegrationStateLayout layout;
  status = WVIntegrationStateLayout::create(
      sourceInspection.latestRestart.coefficientShape, descriptor, layout);
  require(static_cast<bool>(status), status.message);

  WVModel model;
  status = WVModel::create(
      modelOutputCatalog(), sourceInspection.latestRestart.configuration,
      sourceInspection.latestRestart.forcingSchedule, descriptor,
      std::make_unique<WVReferenceFFTEngine>(), {}, model);
  require(static_cast<bool>(status), status.message);
  WVCheckpoint restoredCheckpoint;
  WVAdditionalStateStorage restoredAdditionalState;
  persistence = WVModelOutputNetCDFSink::restoreState(
      sourceInspection, *modelOutputCatalog(), model.stateLayout(),
      restoredCheckpoint, restoredAdditionalState);
  require(static_cast<bool>(persistence), persistence.message);
  auto continuations = sourceInspection.scheduleContinuations;
  auto secondaryContinuation = continuations[0];
  secondaryContinuation.fileIdentifier = "secondary";
  secondaryContinuation.groupIdentifier = "secondary-restart";
  continuations.push_back(std::move(secondaryContinuation));
  WVOutputPlan plan;
  status = WVOutputPlan::create(
      descriptor, modelOutputCatalog(), sourceInspection.latestRestart.t,
      finalTime, continuations, plan);
  require(static_cast<bool>(status) && plan.eventCount() == 1 &&
              plan.event(0).routeCount == 2,
          "WVModel retry plan did not preserve the algorithmic continuation "
          "and two coincident routes");
  WVModelOutputNetCDFSink netCDFSink;
  persistence = WVModelOutputNetCDFSink::createNew(
      {modelOutputCatalog(), restoredCheckpoint,
       sourceInspection.isDynamicsLinear},
      descriptor, plan, layout, nullptr, netCDFSink);
  require(static_cast<bool>(persistence), persistence.message);
  FailOnceNetCDFSink sink(netCDFSink, 2);
  status = netCDFSink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  WVModelState state;
  status = WVModelState::create(std::move(restoredCheckpoint),
                                model.stateLayout(), state,
                                &restoredAdditionalState);
  require(static_cast<bool>(status), status.message);
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);
  const auto timeCount = [](const std::filesystem::path &path,
                            const char *groupName) {
    int file = -1;
    int group = -1;
    int dimension = -1;
    std::size_t count = 0;
    require(nc_open(path.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
                nc_inq_ncid(file, groupName, &group) == NC_NOERR &&
                nc_inq_dimid(group, "t", &dimension) == NC_NOERR &&
                nc_inq_dimlen(group, dimension, &count) == NC_NOERR &&
                nc_close(file) == NC_NOERR,
            "failed to inspect retry destination progress");
    return count;
  };

  status = model.advanceToTime(state, finalTime, 1.0, plan, sink);
  require(!status && state.checkpoint().state.t == finalTime &&
              timeCount(primaryPath, "wave-vortex") == 1 &&
              timeCount(secondaryPath, "wave-vortex-secondary") == 0 &&
              fileBytes(sourcePath) == immutableSource,
          "WVModel route failure lost the accepted state, changed failed "
          "offsets, or mutated the restart source");

  status = model.advanceToTime(state, finalTime, 1.0, plan, sink);
  require(static_cast<bool>(status) &&
              timeCount(primaryPath, "wave-vortex") == 1 &&
              timeCount(secondaryPath, "wave-vortex-secondary") == 1 &&
              fileBytes(sourcePath) == immutableSource,
          "WVModel route retry did not complete exact NetCDF continuation");
  const auto finalEvent = plan.event(0).eventOrdinal;
  const auto primaryAttempts = std::count_if(
      sink.attempts().begin(), sink.attempts().end(),
      [&](const auto &attempt) {
        return attempt.eventOrdinal == finalEvent &&
               attempt.fileIdentifier == "primary";
      });
  const auto secondaryAttempts = std::count_if(
      sink.attempts().begin(), sink.attempts().end(),
      [&](const auto &attempt) {
        return attempt.eventOrdinal == finalEvent &&
               attempt.fileIdentifier == "secondary";
      });
  require(primaryAttempts == 1 && secondaryAttempts == 2,
          "WVModel retry repeated a successful route or skipped the failed "
          "route replay");
  persistence = netCDFSink.close();
  require(static_cast<bool>(persistence), persistence.message);
  WVModelOutputNetCDFInspection retryInspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {primaryPath.string(), secondaryPath.string()}, *modelOutputCatalog(),
      retryInspection);
  require(static_cast<bool>(persistence) &&
              retryInspection.destinationProgress.size() == 2,
          persistence.message);
  for (const auto &progress : retryInspection.destinationProgress) {
    const auto *next =
        progress.committedScheduleCursor.values.value("nextOrdinal");
    require(progress.recordCount == 1 && progress.hasCommittedTime &&
                progress.lastCommittedTime == finalTime && next != nullptr &&
                std::get<std::vector<std::int64_t>>(next->storage)[0] == 3,
            "WVModel retry did not commit the exact typed algorithmic cursor");
  }
}

} // namespace

int main() {
  try {
    (void)modelOutputCatalog();
    testCreateReadAndAppend();
    testLinearInitialCoefficientsAndPassiveFields();
    testTransactionalRefusal();
    testDirectFactoryPreflightBeforeDestinationAccess();
    testMultipleFilesGroupsAndSharedState();
    testOptionalMatlabFixture();
    testOptionalMatlabLinearFixture();
    testOptionalMatlabPassiveFixture();
    testAlgorithmicSchedulePersistence();
    testAlgorithmicScheduleThroughModelAndRequestRunner();
    testVariableObservationBatches();
    testRegisteredTopologyProviders();
    testObservationGraphCollisionPreflight();
    testCoincidentRoutesKeepDistinctLogicalBatches();
    testWVModelRetainsFailedNetCDFRouteForRetry();
    std::cout << "PASS: MATLAB-compatible model-output persistence\n";
    return 0;
  } catch (const std::exception &exception) {
    std::cerr << "FAIL: " << exception.what() << '\n';
    return 1;
  }
}
