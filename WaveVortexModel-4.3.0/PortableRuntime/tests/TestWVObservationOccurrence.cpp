#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"

#include "WVReferenceFFTEngine.hpp"
#include "WVTestObservationOccurrenceProviders.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::test;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

struct OccurrenceTemporaryDirectory {
  std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("wave-vortex-observation-occurrence-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));

  OccurrenceTemporaryDirectory() {
    std::filesystem::create_directories(path);
  }

  ~OccurrenceTemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
};

WVTransformConstantStratificationConfiguration configuration() {
  WVTransformConstantStratificationConfiguration value;
  value.Nx = 8;
  value.Ny = 6;
  value.Nz = 7;
  value.Nj = 6;
  value.Lx = 8000.0;
  value.Ly = 6000.0;
  value.Lz = 1200.0;
  value.N0 = 5.2e-3;
  value.rho0 = 1025.0;
  value.g = 9.81;
  value.planetaryRadius = 6.371e6;
  value.rotationRate = 7.2921e-5;
  value.latitude = 35.0;
  value.shouldAntialias = true;
  return value;
}

struct OwnedState {
  WVShape2D shape;
  std::vector<WVComplex64> Ap;
  std::vector<WVComplex64> Am;
  std::vector<WVComplex64> A0;

  WVState view(double time) const {
    return {time,
            0.0,
            {{Ap.data(), shape}, {Am.data(), shape}, {A0.data(), shape}}};
  }
};

OwnedState state(const WVTransformConstantStratificationConfiguration &config) {
  WVTransformConstantStratificationDescriptor descriptor;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   descriptor)),
          "transform descriptor construction failed");
  OwnedState result;
  result.shape = descriptor.spectralShape();
  result.Ap.resize(result.shape.elementCount());
  result.Am.resize(result.shape.elementCount());
  result.A0.resize(result.shape.elementCount());
  for (std::size_t index = 0; index < result.Ap.size(); ++index) {
    const auto value = static_cast<double>(index + 1);
    result.Ap[index] = {1.1e-3 * std::sin(0.17 * value),
                        8.0e-4 * std::cos(0.11 * value)};
    result.Am[index] = {-7.0e-4 * std::cos(0.09 * value),
                        9.0e-4 * std::sin(0.13 * value)};
    result.A0[index] = {6.0e-4 * std::sin(0.07 * value),
                        5.0e-4 * std::cos(0.19 * value)};
  }
  return result;
}

WVCheckpoint checkpointTemplate(
    const WVTransformConstantStratificationConfiguration &config) {
  const auto owned = state(config);
  WVCheckpoint checkpoint;
  checkpoint.configuration = config;
  checkpoint.state.t = 0.0;
  checkpoint.state.t0 = 0.0;
  checkpoint.state.coefficients.shape = owned.shape;
  checkpoint.state.coefficients.Ap = owned.Ap;
  checkpoint.state.coefficients.Am = owned.Am;
  checkpoint.state.coefficients.A0 = owned.A0;
  checkpoint.metadata.modelVersion = "4.2.1";
  checkpoint.metadata.transformClass = "WVTransformConstantStratification";
  checkpoint.metadata.stateGroupPath = "/";
  checkpoint.metadata.stateCount = 1;
  return checkpoint;
}

constexpr const char *mismatchedCoordinateObservationType =
    "WVTestMismatchedOccurrenceCoordinates";
constexpr const char *affineOccurrenceFieldObservationType =
    "WVTestAffineOccurrenceField";

WVKernelStatus affineOccurrenceFieldOutputPlan(
    const WVObserverRecord &observer,
    const WVObserverOutputPlanningContext &context,
    WVObserverOutputPlan &plan) {
  auto status = observation_occurrence_provider_detail::buildOutputPlan<
      WVTestObservationOccurrenceProviderKind::
          eventVariableOneDimensionalGeometry>(observer, context, plan);
  if (!status)
    return status;
  const auto field = std::find_if(
      plan.channels.begin(), plan.channels.end(), [](const auto &channel) {
        return channel.source ==
               WVObserverOutputChannelSource::occurrenceField;
      });
  if (field == plan.channels.end())
    return {WVKernelStatusCode::invalidConfiguration,
            "Affine occurrence-field proof has no occurrence field."};
  field->scale = observer.outputScale;
  field->offset = observer.outputOffset;
  return WVKernelStatus::ok();
}

class AffineOccurrenceFieldObservation final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept final {
    static const std::string identifier = affineOccurrenceFieldObservationType;
    return identifier;
  }
  std::uint32_t contractVersion() const noexcept final { return 1; }

  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const final {
    return observer.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : WVKernelStatus{
                     WVKernelStatusCode::invalidConfiguration,
                     "Affine occurrence-field proof unexpectedly owns state."};
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
    return affineOccurrenceFieldOutputPlan(observer, context, plan);
  }

  std::size_t persistentBytes() const noexcept final { return sizeof(*this); }
};

class MismatchedCoordinateObservation final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept final {
    static const std::string identifier =
        mismatchedCoordinateObservationType;
    return identifier;
  }
  std::uint32_t contractVersion() const noexcept final { return 1; }

  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const final {
    return observer.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : WVKernelStatus{
                     WVKernelStatusCode::invalidConfiguration,
                     "Mismatched-coordinate proof unexpectedly owns state."};
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
    return observation_occurrence_provider_detail::buildOutputPlan<
        WVTestObservationOccurrenceProviderKind::
            eventVariableOneDimensionalGeometry>(observer, context, plan);
  }

  WVKernelStatus prepareOccurrence(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOccurrencePreparationContext &context,
      WVObserverOccurrenceWorkspace &workspace) const final {
    auto status = observation_occurrence_provider_detail::prepareWorkspace<
        WVTestObservationOccurrenceProviderKind::
            eventVariableOneDimensionalGeometry>(observer, plan, context,
                                                  workspace);
    if (status && !workspace.positionSets.empty() &&
        !workspace.positionSets.front().x.empty())
      workspace.positionSets.front().x.front() += 1.0;
    return status;
  }

  WVKernelStatus observationBatch(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &batch) const final {
    return WVObservingSystem::observationBatch(observer, plan, context, kind,
                                               batch);
  }

  std::size_t persistentBytes() const noexcept final { return sizeof(*this); }
};

std::shared_ptr<const WVExtensionCatalog> extensionCatalog() {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (status)
    status = registerOccurrenceSchedule(builder);
  if (status)
    status = registerObservationOccurrenceProviders(builder);
  if (status)
    status = builder.addObserverFactory(
        {mismatchedCoordinateObservationType, 1,
         [](const WVObserverRecord &, const WVPortableTypedRecord &,
            std::shared_ptr<const WVObservingSystem> &result) {
           result = std::make_shared<MismatchedCoordinateObservation>();
           return WVKernelStatus::ok();
         },
         {},
         [](const WVObserverRecord &observer,
            const WVObserverOutputPlanningContext &context,
            WVObserverOutputPlan &plan) {
           return observation_occurrence_provider_detail::buildOutputPlan<
               WVTestObservationOccurrenceProviderKind::
                   eventVariableOneDimensionalGeometry>(observer, context,
                                                         plan);
         }});
  if (status)
    status = builder.addObserverFactory(
        {affineOccurrenceFieldObservationType, 1,
         [](const WVObserverRecord &, const WVPortableTypedRecord &,
            std::shared_ptr<const WVObservingSystem> &result) {
           result = std::make_shared<AffineOccurrenceFieldObservation>();
           return WVKernelStatus::ok();
         },
         {}, affineOccurrenceFieldOutputPlan});
  std::shared_ptr<const WVExtensionCatalog> result;
  if (status)
    status = builder.freeze(result);
  require(static_cast<bool>(status),
          "occurrence extension catalog construction failed: " +
              status.message);
  return result;
}

WVPortableObserverDescriptor descriptor(
    const std::shared_ptr<const WVExtensionCatalog> &catalog) {
  WVPortableObserverRecord record;
  const auto addCoordinateBlock = [&](std::string identifier,
                                      std::size_t count) {
    WVStateBlockRecord block;
    block.identifier = std::move(identifier);
    block.scalarType = WVStateScalarType::real64;
    block.dimensions = {count};
    block.toleranceKind = WVToleranceKind::uniformAbsolute;
    block.absoluteTolerance = 1e-8;
    block.ownership = WVStateOwnership::integratorOwned;
    block.restartRequirement = WVRestartRequirement::requiredDynamicState;
    record.stateBlocks.push_back(std::move(block));
  };
  addCoordinateBlock("unrelated-particle-x", 2);
  addCoordinateBlock("unrelated-particle-y", 2);
  addCoordinateBlock("event-coordinate-x", 4);
  addCoordinateBlock("event-coordinate-y", 4);
  addCoordinateBlock("event-coordinate-z", 4);
  record.observers = {fixedGeometryObservationRecord(),
                      eventVariableOneDimensionalGeometryRecord(),
                      variableSamplesByFixedDepthBinsRecord(),
                      stateCoupledIrregularGeometryRecord(),
                      nestedRaggedGeometryRecord()};
  WVObserverRecord unrelatedParticles;
  unrelatedParticles.identifier = "unrelated-particles";
  unrelatedParticles.name = "unrelated_particles";
  unrelatedParticles.typeIdentifier = "WVLagrangianParticles";
  unrelatedParticles.stateBlockIdentifiers = {"unrelated-particle-x",
                                               "unrelated-particle-y"};
  unrelatedParticles.x = {11.0, 22.0};
  unrelatedParticles.y = {33.0, 44.0};
  unrelatedParticles.z = {-55.0, -66.0};
  unrelatedParticles.isXYOnly = true;
  unrelatedParticles.horizontalAbsoluteTolerance = 1e-8;
  record.observers.push_back(std::move(unrelatedParticles));
  WVPortableObserverDescriptor result;
  const auto status =
      WVPortableObserverDescriptor::create(record, catalog, result);
  require(static_cast<bool>(status),
          "occurrence observer descriptor construction failed: " +
              status.message);
  return result;
}

WVPortableObserverDescriptor protocolDescriptor(
    const std::shared_ptr<const WVExtensionCatalog> &catalog,
    double finalTime, bool useOccurrenceSchedule = true) {
  const auto config = configuration();
  WVTransformConstantStratificationDescriptor transform;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   transform)),
          "protocol transform descriptor construction failed");
  const auto shape = transform.spectralShape();
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::complex64,
         {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled,
         0.0,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers = {std::move(coefficients),
                      eventVariableOneDimensionalGeometryRecord()};

  const auto outputFile = [&](std::string identifier,
                              std::string destination) {
    WVOutputGroupRecord restart;
    restart.identifier = "shared-restart";
    restart.name = "Shared restart";
    restart.schedule = {0.5, 0.0, finalTime};
    restart.observerIdentifiers = {"coefficients"};
    restart.containsCompleteCoefficientRestart = true;
    WVOutputGroupRecord observations;
    observations.identifier = "shared-occurrence";
    observations.name = "Shared occurrence";
    observations.schedule =
        useOccurrenceSchedule
            ? occurrenceSchedule(3.5, 3, 1, finalTime, 0.0, 0.25)
            : WVOutputScheduleRecord{0.25, 0.0, finalTime};
    observations.observerIdentifiers = {"event-variable-line"};
    WVOutputFileRecord file;
    file.identifier = std::move(identifier);
    file.destination = std::move(destination);
    file.groups = {std::move(restart), std::move(observations)};
    return file;
  };
  record.outputFiles = {outputFile("primary", "primary.memory"),
                        outputFile("secondary", "secondary.memory")};
  WVPortableObserverDescriptor result;
  const auto status =
      WVPortableObserverDescriptor::create(record, catalog, result);
  require(static_cast<bool>(status),
          "protocol observer descriptor construction failed: " +
              status.message);
  return result;
}

WVPortableObserverDescriptor stateCoupledProtocolDescriptor(
    const std::shared_ptr<const WVExtensionCatalog> &catalog,
    double finalTime) {
  const auto config = configuration();
  WVTransformConstantStratificationDescriptor transform;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   transform)),
          "state-coupled protocol transform descriptor construction failed");
  const auto shape = transform.spectralShape();
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::complex64,
         {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled,
         0.0,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});

  auto stateCoupled = stateCoupledIrregularGeometryRecord();
  for (const auto &identifier : stateCoupled.stateBlockIdentifiers)
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::real64,
         {stateCoupled.x.size()},
         WVToleranceKind::uniformAbsolute,
         1e-8,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});

  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers = {std::move(coefficients), std::move(stateCoupled)};

  WVOutputGroupRecord restart;
  restart.identifier = "restart";
  restart.name = "Restart";
  restart.schedule = {0.5, 0.0, finalTime};
  restart.observerIdentifiers = {"coefficients"};
  restart.containsCompleteCoefficientRestart = true;
  WVOutputGroupRecord observations;
  observations.identifier = "state-coupled-occurrence";
  observations.name = "State-coupled occurrence";
  observations.schedule =
      occurrenceSchedule(4.0, 3, 1, finalTime, 0.0, 0.25);
  observations.observerIdentifiers = {"state-coupled-irregular"};
  WVOutputFileRecord file;
  file.identifier = "state-coupled";
  file.destination = "state-coupled.memory";
  file.groups = {std::move(restart), std::move(observations)};
  record.outputFiles.push_back(std::move(file));

  WVPortableObserverDescriptor result;
  const auto status =
      WVPortableObserverDescriptor::create(record, catalog, result);
  require(static_cast<bool>(status),
          "state-coupled protocol descriptor construction failed: " +
              status.message);
  return result;
}

WVPortableObserverDescriptor nestedRaggedNetCDFDescriptor(
    const std::shared_ptr<const WVExtensionCatalog> &catalog,
    const std::filesystem::path &destination, double finalTime) {
  const auto config = configuration();
  WVTransformConstantStratificationDescriptor transform;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   transform)),
          "nested-ragged transform descriptor construction failed");
  const auto shape = transform.spectralShape();
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::complex64,
         {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled,
         1e-8,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});

  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers = {std::move(coefficients),
                      nestedRaggedGeometryRecord()};

  WVOutputGroupRecord restart;
  restart.identifier = "restart";
  restart.name = "restart";
  restart.schedule = {0.5, 0.0, finalTime};
  restart.observerIdentifiers = {"coefficients"};
  restart.containsCompleteCoefficientRestart = true;
  WVOutputGroupRecord nested;
  nested.identifier = "nested-ragged";
  nested.name = "nested-ragged";
  nested.schedule = occurrenceSchedule(4.0, 1, 1, finalTime, 0.0, 0.25);
  nested.observerIdentifiers = {"nested-ragged"};
  WVOutputFileRecord file;
  file.identifier = "nested-ragged-file";
  file.destination = destination.string();
  file.groups = {std::move(restart), std::move(nested)};
  record.outputFiles.push_back(std::move(file));

  WVPortableObserverDescriptor result;
  const auto status =
      WVPortableObserverDescriptor::create(record, catalog, result);
  require(static_cast<bool>(status),
          "nested-ragged NetCDF descriptor construction failed: " +
              status.message);
  return result;
}

WVPortableObserverDescriptor primitiveReuseDescriptor(
    const std::shared_ptr<const WVExtensionCatalog> &catalog) {
  const auto config = configuration();
  WVTransformConstantStratificationDescriptor transform;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   transform)),
          "primitive-reuse transform descriptor construction failed");
  const auto shape = transform.spectralShape();
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::complex64,
         {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled,
         0.0,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers = {std::move(coefficients),
                      variableSamplesByFixedDepthBinsRecord()};

  WVOutputGroupRecord restart;
  restart.identifier = "restart";
  restart.name = "Restart";
  restart.schedule = {1.0, 0.0, 0.0};
  restart.observerIdentifiers = {"coefficients"};
  restart.containsCompleteCoefficientRestart = true;
  const auto occurrenceGroup = [](std::string identifier,
                                  double sourceReal) {
    WVOutputGroupRecord group;
    group.name = "Coincident variable-depth occurrence " + identifier;
    group.identifier = std::move(identifier);
    group.schedule = occurrenceSchedule(sourceReal, 1, 1, 0.0);
    group.observerIdentifiers = {"variable-depth-bins"};
    return group;
  };
  WVOutputFileRecord file;
  file.identifier = "primitive-reuse";
  file.destination = "primitive-reuse.memory";
  file.groups = {std::move(restart), occurrenceGroup("first", 2.0),
                 occurrenceGroup("second", 7.0)};
  record.outputFiles.push_back(std::move(file));

  WVPortableObserverDescriptor result;
  const auto status =
      WVPortableObserverDescriptor::create(record, catalog, result);
  require(static_cast<bool>(status),
          "primitive-reuse observer descriptor construction failed: " +
              status.message);
  return result;
}

const WVObservationValue &value(const WVObservationBatch &batch,
                                const char *identifier) {
  const auto found =
      std::find_if(batch.values.begin(), batch.values.end(),
                   [&](const auto &candidate) {
                     return candidate.variableIdentifier == identifier;
                   });
  require(found != batch.values.end(),
          std::string("missing occurrence value ") + identifier);
  return *found;
}

std::size_t sumIntegers(const WVObservationValue &input) {
  require(input.scalarType == WVObservationScalarType::integer64,
          "ragged relationship is not integer-valued");
  const auto *data = input.integer64Data();
  std::size_t sum = 0;
  for (std::size_t index = 0; index < input.elementCount(); ++index) {
    require(data[index] >= 0, "ragged relationship is negative");
    sum += static_cast<std::size_t>(data[index]);
  }
  return sum;
}

struct ResolvedOccurrence {
  std::shared_ptr<const WVOutputSchedule> schedule;
  WVOutputScheduleOccurrence value;
};

ResolvedOccurrence occurrence(
    const std::shared_ptr<const WVExtensionCatalog> &catalog,
    const WVOutputScheduleRecord &record) {
  ResolvedOccurrence result;
  auto status = catalog->outputSchedules().resolve(record, result.schedule);
  require(static_cast<bool>(status), "schedule resolution failed");
  bool available = false;
  status = result.schedule->peek({}, record.initialTime, record.initialTime,
                                 result.value, available);
  require(static_cast<bool>(status) && available,
          "schedule did not produce its initial occurrence");
  return result;
}

struct PreparedBatch {
  WVObservationSchema schema;
  WVObservationBatch batch;
  WVObservationOccurrenceIdentity identity;
};

PreparedBatch prepareBatch(
    WVObserverOutputEvaluationService &service,
    const WVPortableObserverDescriptor &observers, std::size_t observerOrdinal,
    const WVOutputSchedulePayloadSchema &payloadSchema,
    const WVOutputScheduleOccurrence &occurrenceValue,
    const WVIntegrationState &stateView, std::size_t eventOrdinal,
    std::size_t semanticScheduleOrdinal) {
  const auto &observer = observers.observers().at(observerOrdinal);
  const auto *resolved = observers.resolvedObserver(observer);
  require(resolved != nullptr, "proof observer is unresolved");
  WVOutputObserverView observerView{observerOrdinal, &observer, resolved};
  WVOutputRouteView route;
  route.scheduleOrdinal = occurrenceValue.ordinal;
  route.observers = &observerView;
  route.observerCount = 1;
  route.semanticScheduleOrdinal = semanticScheduleOrdinal;
  route.proposedScheduleCursor = &occurrenceValue.proposedCursor.values;
  route.schedulePayloadSchema = &payloadSchema;
  route.schedulePayload = &occurrenceValue.payload;
  route.scheduleCursorIdentity = occurrenceValue.cursorIdentity;
  WVOutputEvent event;
  event.eventOrdinal = eventOrdinal;
  event.scheduledTime = occurrenceValue.scheduledTime;
  event.kind = WVOutputEventKind::acceptedEndpoint;
  event.state = stateView;
  event.state.waveVortex.t = occurrenceValue.scheduledTime;
  event.routes = &route;
  event.routeCount = 1;

  auto status = service.prepare(event);
  require(static_cast<bool>(status),
          "occurrence preparation failed: " + status.message);
  const auto preparedCount = service.metrics().occurrencePreparationCount;
  status = service.prepare(event);
  require(static_cast<bool>(status),
          "idempotent occurrence preparation failed");
  require(service.metrics().occurrencePreparationCount == preparedCount,
          "retry rediscovered an already prepared occurrence");

  PreparedBatch result;
  status = service.preparedOccurrenceIdentity(route, observerView,
                                              result.identity);
  require(static_cast<bool>(status),
          "prepared occurrence identity lookup failed");
  status = service.observationSchema(observer, result.schema);
  require(static_cast<bool>(status), "occurrence schema lookup failed");
  status = service.observationBatch(result.identity, observer, result.batch);
  require(static_cast<bool>(status),
          "occurrence batch construction failed: " + status.message);
  status = validateObservationBatch(result.schema, result.batch);
  require(static_cast<bool>(status),
          "occurrence batch validation failed: " + status.message);
  return result;
}

void compareField(const WVObservationBatch &batch,
                  WVFieldEvaluationService &fieldService, const WVState &state,
                  const std::string &variableIdentifier,
                  const std::string &fieldName,
                  WVPositionInterpolation interpolation =
                      WVPositionInterpolation::linear) {
  const auto &observed = value(batch, variableIdentifier.c_str());
  const auto &x = value(batch, "x");
  const auto &y = value(batch, "y");
  const auto &z = value(batch, "z");
  require(observed.extents == x.extents && x.extents == y.extents &&
              y.extents == z.extents,
          "dynamic field and coordinate extents differ");
  if (observed.elementCount() == 0)
    return;
  WVFieldSamplingRequest sampling;
  sampling.kind = WVFieldSamplingKind::positions;
  sampling.x.assign(x.real64Data(), x.real64Data() + x.elementCount());
  sampling.y.assign(y.real64Data(), y.real64Data() + y.elementCount());
  sampling.z.assign(z.real64Data(), z.real64Data() + z.elementCount());
  sampling.interpolation = interpolation;
  WVFieldEvaluationPlan plan;
  auto status = fieldService.createPlan(
      {{"expected", fieldName, std::move(sampling)}}, plan);
  require(static_cast<bool>(status),
          "comparison field plan construction failed for " + fieldName +
              ": " + status.message);
  require(plan.outputCount() == 1 &&
              plan.outputs().front().elementCount == observed.elementCount(),
          "comparison field plan shape differs");
  std::vector<double> expected(observed.elementCount());
  WVFieldOutputView output{expected.data(), expected.size()};
  status = fieldService.evaluate(plan, state, &output, 1);
  require(static_cast<bool>(status),
          "comparison field evaluation failed for " + fieldName + ": " +
              status.message);
  const auto *actual = observed.real64Data();
  for (std::size_t index = 0; index < expected.size(); ++index)
    require(std::abs(actual[index] - expected[index]) <= 1e-12,
            variableIdentifier + " differs from fixed-geometry evaluation");
}

void requireCoordinatesFinite(const WVObservationBatch &batch) {
  for (const char *identifier : {"sample-time", "x", "y", "z"}) {
    const auto &coordinate = value(batch, identifier);
    const auto *data = coordinate.real64Data();
    for (std::size_t index = 0; index < coordinate.elementCount(); ++index)
      require(std::isfinite(data[index]),
              std::string(identifier) + " contains a nonfinite value");
  }
}

class ZeroErrorPolicy final : public WVIntegrationErrorPolicy {
public:
  explicit ZeroErrorPolicy(const WVIntegrationStateLayout &layout)
      : layout_(layout) {}

  std::size_t componentCount() const noexcept override {
    return 3 + layout_.additionalBlocks().size();
  }
  std::size_t elementCount(std::size_t component) const noexcept override {
    return component < 3
               ? layout_.coefficientShape().elementCount()
               : layout_.additionalBlocks()[component - 3].elementCount;
  }
  double absoluteTolerance(std::size_t, std::size_t) const noexcept override {
    return 1e-10;
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }

private:
  const WVIntegrationStateLayout &layout_;
};

class ZeroIntegrationSystem final : public WVIntegrationSystem {
public:
  explicit ZeroIntegrationSystem(WVIntegrationStateLayout layout)
      : layout_(std::move(layout)) {}

  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(
      const WVIntegrationState &, WVIntegrationFlux &rightHandSide) override {
    const auto coefficientCount = layout_.coefficientShape().elementCount();
    for (auto output : {rightHandSide.waveVortex.Fp,
                        rightHandSide.waveVortex.Fm,
                        rightHandSide.waveVortex.F0})
      std::fill_n(output.data, coefficientCount, WVComplex64{});
    for (std::size_t block = 0;
         block < rightHandSide.additionalBlockCount; ++block) {
      const auto &metadata = layout_.additionalBlocks()[block];
      if (metadata.scalarType == WVStateScalarType::real64)
        std::fill_n(rightHandSide.additionalBlocks[block].realData,
                    metadata.elementCount, 0.0);
      else
        std::fill_n(rightHandSide.additionalBlocks[block].complexData,
                    metadata.elementCount, WVComplex64{});
    }
    return WVKernelStatus::ok();
  }
  WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &) override {
    return {WVKernelStatus::ok(), 0, true};
  }
  WVKernelStatus createErrorPolicy(
      double, std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override {
    policy = std::make_unique<ZeroErrorPolicy>(layout_);
    return WVKernelStatus::ok();
  }

private:
  WVIntegrationStateLayout layout_;
};

struct MutableState {
  WVShape2D shape;
  std::vector<WVComplex64> Ap;
  std::vector<WVComplex64> Am;
  std::vector<WVComplex64> A0;
  WVAdditionalStateStorage additional;
  WVMutableIntegrationState view;

  explicit MutableState(const WVIntegrationStateLayout &layout)
      : shape(layout.coefficientShape()), Ap(shape.elementCount()),
        Am(shape.elementCount()), A0(shape.elementCount()) {
    require(static_cast<bool>(additional.initialize(layout)),
            "protocol state storage initialization failed");
    for (std::size_t index = 0; index < shape.elementCount(); ++index) {
      const auto scalar = static_cast<double>(index + 1);
      Ap[index] = {1e-3 * std::sin(0.13 * scalar),
                   1e-3 * std::cos(0.17 * scalar)};
      Am[index] = {-8e-4 * std::cos(0.11 * scalar),
                   7e-4 * std::sin(0.19 * scalar)};
      A0[index] = {6e-4 * std::sin(0.07 * scalar),
                   5e-4 * std::cos(0.23 * scalar)};
    }
    view = {{0.0,
             0.0,
             {{Ap.data(), shape}, {Am.data(), shape}, {A0.data(), shape}}},
            additional.mutableBlocks(), additional.blockCount()};
  }
};

class OccurrenceCachingSink final : public WVOutputSink {
public:
  explicit OccurrenceCachingSink(WVObserverOutputEvaluationService &source)
      : source_(source) {}

  WVKernelStatus preflight(const WVOutputPlan &plan) override {
    return source_.preflight(plan);
  }

  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    auto status = source_.prepare(event);
    if (!status)
      return status;
    if (!hasEvent_ || event.eventOrdinal != eventOrdinal_ ||
        event.scheduledTime != scheduledTime_) {
      cached_.clear();
      hasEvent_ = true;
      eventOrdinal_ = event.eventOrdinal;
      scheduledTime_ = event.scheduledTime;
    }
    for (std::size_t observerIndex = 0;
         observerIndex < route.observerCount; ++observerIndex) {
      const auto &observer = route.observers[observerIndex];
      WVObservationOccurrenceIdentity identity;
      status = source_.preparedOccurrenceIdentity(route, observer, identity);
      if (!status)
        return status;
      const auto found = std::find_if(
          cached_.begin(), cached_.end(), [&](const auto &candidate) {
            return samePreparedObservationOccurrenceIdentity(
                candidate.identity, identity);
          });
      if (found == cached_.end()) {
        CachedOccurrence occurrence;
        occurrence.identity = identity;
        status = source_.observationBatch(identity, *observer.record,
                                          occurrence.batch);
        if (!status)
          return status;
        const auto sampleTimes = std::find_if(
            occurrence.batch.values.begin(), occurrence.batch.values.end(),
            [](const auto &candidate) {
              return candidate.variableIdentifier == "sample-time";
            });
        if (sampleTimes != occurrence.batch.values.end())
          for (std::size_t sample = 0;
               sample < sampleTimes->elementCount(); ++sample)
            sawDistinctSampleTime_ =
                sawDistinctSampleTime_ ||
                sampleTimes->real64Data()[sample] != event.scheduledTime;
        cached_.push_back(std::move(occurrence));
        ++batchBuildCount_;
      } else {
        ++batchReuseCount_;
      }
    }
    ++attemptCount_;
    if (event.eventOrdinal == 0 &&
        route.fileOrdinal < firstEventAttempts_.size() &&
        route.groupOrdinal < firstEventAttempts_[route.fileOrdinal].size())
      ++firstEventAttempts_[route.fileOrdinal][route.groupOrdinal];
    if (attemptCount_ == failureAttempt_)
      return {WVKernelStatusCode::numericalFailure,
              "injected later-route occurrence failure"};
    result.writeCount = route.observerCount;
    ++successfulRouteCount_;
    return WVKernelStatus::ok();
  }

  void disableFailure() noexcept {
    failureAttempt_ = std::numeric_limits<std::size_t>::max();
  }
  void failAt(std::size_t attempt) noexcept { failureAttempt_ = attempt; }
  std::size_t batchBuildCount() const noexcept { return batchBuildCount_; }
  std::size_t batchReuseCount() const noexcept { return batchReuseCount_; }
  std::size_t successfulRouteCount() const noexcept {
    return successfulRouteCount_;
  }
  const std::array<std::array<std::size_t, 2>, 2> &
  firstEventAttempts() const noexcept {
    return firstEventAttempts_;
  }
  bool sawDistinctSampleTime() const noexcept {
    return sawDistinctSampleTime_;
  }

private:
  struct CachedOccurrence {
    WVObservationOccurrenceIdentity identity;
    WVObservationBatch batch;
  };
  WVObserverOutputEvaluationService &source_;
  std::vector<CachedOccurrence> cached_;
  std::size_t eventOrdinal_ = 0;
  double scheduledTime_ = 0.0;
  std::size_t attemptCount_ = 0;
  std::size_t failureAttempt_ = std::numeric_limits<std::size_t>::max();
  std::size_t batchBuildCount_ = 0;
  std::size_t batchReuseCount_ = 0;
  std::size_t successfulRouteCount_ = 0;
  std::array<std::array<std::size_t, 2>, 2> firstEventAttempts_{};
  bool hasEvent_ = false;
  bool sawDistinctSampleTime_ = false;
};

struct StateCoupledOccurrenceSnapshot {
  double scheduledTime = 0.0;
  WVOutputEventKind kind = WVOutputEventKind::acceptedEndpoint;
  std::vector<double> stateX;
  std::vector<double> stateY;
  std::vector<double> stateZ;
  std::vector<double> batchX;
  std::vector<double> batchY;
  std::vector<double> batchZ;
  std::vector<double> sampledV;
  std::vector<double> sampledSsv;
};

class StateCoupledOccurrenceSink final : public WVOutputSink {
public:
  explicit StateCoupledOccurrenceSink(
      WVObserverOutputEvaluationService &source)
      : source_(source) {}

  WVKernelStatus preflight(const WVOutputPlan &plan) override {
    return source_.preflight(plan);
  }

  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    auto status = source_.prepare(event);
    if (!status)
      return status;
    for (std::size_t observerIndex = 0;
         observerIndex < route.observerCount; ++observerIndex) {
      const auto &observerView = route.observers[observerIndex];
      if (observerView.record == nullptr ||
          observerView.record->identifier != "state-coupled-irregular")
        continue;

      WVObservationOccurrenceIdentity identity;
      status = source_.preparedOccurrenceIdentity(route, observerView,
                                                  identity);
      if (!status)
        return status;
      WVObservationBatch batch;
      status = source_.observationBatch(identity, *observerView.record, batch);
      if (!status)
        return status;
      WVObservationSchema schema;
      status = source_.observationSchema(*observerView.record, schema);
      if (status)
        status = validateObservationBatch(schema, batch);
      if (!status)
        return status;

      StateCoupledOccurrenceSnapshot snapshot;
      snapshot.scheduledTime = event.scheduledTime;
      snapshot.kind = event.kind;
      const auto copyBatchValue = [&](const char *identifier,
                                      std::vector<double> &output) {
        const auto found = std::find_if(
            batch.values.begin(), batch.values.end(), [&](const auto &item) {
              return item.variableIdentifier == identifier;
            });
        if (found == batch.values.end() ||
            found->scalarType != WVObservationScalarType::real64)
          return WVKernelStatus{
              WVKernelStatusCode::invalidConfiguration,
              std::string("State-coupled batch has no real ") + identifier +
                  " value."};
        const auto count = found->elementCount();
        const auto *data = found->real64Data();
        if (count != 0 && data == nullptr)
          return WVKernelStatus{
              WVKernelStatusCode::invalidConfiguration,
              std::string("State-coupled batch has no ") + identifier +
                  " storage."};
        output.clear();
        if (count != 0)
          output.assign(data, data + count);
        return WVKernelStatus::ok();
      };
      const auto copyStateBlock = [&](const std::string &identifier,
                                      std::vector<double> &output) {
        for (std::size_t block = 0;
             block < event.state.additionalBlockCount; ++block) {
          const auto &candidate = event.state.additionalBlocks[block];
          if (candidate.layout != nullptr &&
              candidate.layout->identifier == identifier) {
            if (candidate.realData == nullptr)
              return WVKernelStatus{
                  WVKernelStatusCode::invalidConfiguration,
                  "State-coupled event block has no real storage."};
            output.assign(candidate.realData,
                          candidate.realData + candidate.layout->elementCount);
            return WVKernelStatus::ok();
          }
        }
        return WVKernelStatus{
            WVKernelStatusCode::invalidConfiguration,
            "State-coupled event block is unavailable."};
      };

      status = copyStateBlock(observerView.record->stateBlockIdentifiers[0],
                              snapshot.stateX);
      if (status)
        status = copyStateBlock(observerView.record->stateBlockIdentifiers[1],
                                snapshot.stateY);
      if (status)
        status = copyStateBlock(observerView.record->stateBlockIdentifiers[2],
                                snapshot.stateZ);
      if (status)
        status = copyBatchValue("x", snapshot.batchX);
      if (status)
        status = copyBatchValue("y", snapshot.batchY);
      if (status)
        status = copyBatchValue("z", snapshot.batchZ);
      if (status)
        status = copyBatchValue("sampled-v", snapshot.sampledV);
      if (status)
        status = copyBatchValue("sampled-ssv", snapshot.sampledSsv);
      if (!status)
        return status;
      snapshots_.push_back(std::move(snapshot));
    }
    result.writeCount = route.observerCount;
    return WVKernelStatus::ok();
  }

  const std::vector<StateCoupledOccurrenceSnapshot> &snapshots() const
      noexcept {
    return snapshots_;
  }

private:
  WVObserverOutputEvaluationService &source_;
  std::vector<StateCoupledOccurrenceSnapshot> snapshots_;
};

struct NestedRaggedPersistence {
  std::size_t passCount = 0;
  std::size_t profileCount = 0;
  std::size_t sampleCount = 0;
  std::vector<double> times;
  std::vector<long long> scheduleOrdinals;
  std::vector<long long> committedPass;
  std::vector<long long> committedProfile;
  std::vector<long long> committedSample;
  std::vector<long long> profileCounts;
  std::vector<long long> sampleCounts;
  std::vector<long long> passIdentifiers;
  std::vector<long long> profileIdentifiers;
  std::vector<long long> sampleIdentifiers;
  std::vector<double> sampleTimes;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<unsigned char> directions;
  std::vector<double> payloadReal;
  std::vector<long long> payloadInteger;
  std::vector<unsigned char> payloadBoolean;
  std::vector<double> sampledU;
  std::vector<double> sampledSsh;
};

NestedRaggedPersistence
readNestedRaggedPersistence(const std::filesystem::path &path) {
  int file = -1;
  int group = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &file) == NC_NOERR &&
              nc_inq_ncid(file, "nested-ragged", &group) == NC_NOERR,
          "unable to open nested-ragged output group");
  const auto dimensionLength = [&](const char *name) {
    int dimension = -1;
    std::size_t length = 0;
    require(nc_inq_dimid(group, name, &dimension) == NC_NOERR &&
                nc_inq_dimlen(group, dimension, &length) == NC_NOERR,
            std::string("unable to inspect nested-ragged dimension ") + name);
    return length;
  };
  const auto readDouble = [&](const char *name, std::size_t count) {
    int variable = -1;
    std::vector<double> result(count);
    require(nc_inq_varid(group, name, &variable) == NC_NOERR &&
                (count == 0 ||
                 nc_get_var_double(group, variable, result.data()) ==
                     NC_NOERR),
            std::string("unable to read nested-ragged variable ") + name);
    return result;
  };
  const auto readInteger = [&](const char *name, std::size_t count) {
    int variable = -1;
    std::vector<long long> result(count);
    require(nc_inq_varid(group, name, &variable) == NC_NOERR &&
                (count == 0 ||
                 nc_get_var_longlong(group, variable, result.data()) ==
                     NC_NOERR),
            std::string("unable to read nested-ragged variable ") + name);
    return result;
  };
  const auto readBoolean = [&](const char *name, std::size_t count) {
    int variable = -1;
    std::vector<unsigned char> result(count);
    require(nc_inq_varid(group, name, &variable) == NC_NOERR &&
                (count == 0 ||
                 nc_get_var_uchar(group, variable, result.data()) ==
                     NC_NOERR),
            std::string("unable to read nested-ragged variable ") + name);
    return result;
  };

  NestedRaggedPersistence result;
  const auto recordCount = dimensionLength("t");
  result.passCount = dimensionLength("pass");
  result.profileCount = dimensionLength("profile");
  result.sampleCount = dimensionLength("sample");
  result.times = readDouble("t", recordCount);
  result.scheduleOrdinals =
      readInteger("portableScheduleOrdinal", recordCount);
  result.committedPass = readInteger("portableCommitted_pass", recordCount);
  result.committedProfile =
      readInteger("portableCommitted_profile", recordCount);
  result.committedSample =
      readInteger("portableCommitted_sample", recordCount);
  result.profileCounts =
      readInteger("profile-count-by-pass", result.passCount);
  result.sampleCounts =
      readInteger("sample-count-by-profile", result.profileCount);
  result.passIdentifiers = readInteger("pass-identifier", result.passCount);
  result.profileIdentifiers =
      readInteger("profile-identifier", result.profileCount);
  result.sampleIdentifiers =
      readInteger("sample-identifier", result.sampleCount);
  result.sampleTimes = readDouble("sample-time", result.sampleCount);
  result.x = readDouble("x", result.sampleCount);
  result.y = readDouble("y", result.sampleCount);
  result.z = readDouble("z", result.sampleCount);
  result.directions =
      readBoolean("sample-direction", result.sampleCount);
  result.payloadReal = readDouble("payload-real", recordCount);
  result.payloadInteger = readInteger("payload-integer", recordCount);
  result.payloadBoolean = readBoolean("payload-boolean", recordCount);
  result.sampledU = readDouble("sampled-u", result.sampleCount);
  result.sampledSsh = readDouble("sampled-ssh", result.sampleCount);
  require(nc_close(file) == NC_NOERR,
          "unable to close nested-ragged output group");
  return result;
}

template <class Value>
bool isPrefix(const std::vector<Value> &prefix,
              const std::vector<Value> &values) {
  return prefix.size() <= values.size() &&
         std::equal(prefix.begin(), prefix.end(), values.begin());
}

void testProviders() {
  resetObservationOccurrenceProviderCounters();
  resetOccurrenceScheduleCounters();
  const auto config = configuration();
  const auto catalog = extensionCatalog();
  const auto observers = descriptor(catalog);
  auto ownedState = state(config);

  std::unique_ptr<WVObserverOutputEvaluationService> service;
  auto status = WVObserverOutputEvaluationService::create(
      config, false, observers, std::make_unique<WVReferenceFFTEngine>(),
      service);
  require(static_cast<bool>(status),
          "occurrence evaluation service construction failed: " +
              status.message);
  std::unique_ptr<WVFieldEvaluationService> comparisonFields;
  status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), comparisonFields);
  require(static_cast<bool>(status),
          "comparison field service construction failed");

  const auto stateLayout = [](std::string identifier, std::size_t count) {
    WVAdditionalStateBlockLayout result;
    result.identifier = std::move(identifier);
    result.scalarType = WVStateScalarType::real64;
    result.dimensions = {count};
    result.ownership = WVStateOwnership::integratorOwned;
    result.restartRequirement = WVRestartRequirement::requiredDynamicState;
    result.elementCount = count;
    return result;
  };
  const std::array<WVAdditionalStateBlockLayout, 5> stateLayouts{{
      stateLayout("unrelated-particle-x", 2),
      stateLayout("unrelated-particle-y", 2),
      stateLayout("event-coordinate-x", 4),
      stateLayout("event-coordinate-y", 4),
      stateLayout("event-coordinate-z", 4)}};
  const std::array<double, 2> unrelatedX{{-9876.5, -8765.4}};
  const std::array<double, 2> unrelatedY{{-7654.3, -6543.2}};
  const std::array<double, 4> coordinateX{{700.0, 1900.0, 3600.0, 6200.0}};
  const std::array<double, 4> coordinateY{{900.0, 2200.0, 4100.0, 5200.0}};
  const std::array<double, 4> coordinateZ{{-100.0, -350.0, -700.0,
                                           -1050.0}};
  const std::array<WVAdditionalStateBlockConstView, 5> stateBlocks{{
      {&stateLayouts[0], unrelatedX.data(), nullptr},
      {&stateLayouts[1], unrelatedY.data(), nullptr},
      {&stateLayouts[2], coordinateX.data(), nullptr},
      {&stateLayouts[3], coordinateY.data(), nullptr},
      {&stateLayouts[4], coordinateZ.data(), nullptr}}};
  WVIntegrationState integrationState;
  integrationState.waveVortex = ownedState.view(2.0);
  integrationState.additionalBlocks = stateBlocks.data();
  integrationState.additionalBlockCount = stateBlocks.size();

  WVOutputScheduleRecord fixedRecord{1.0, 2.0, 2.0};
  const auto fixedOccurrence = occurrence(catalog, fixedRecord);
  auto fixed = prepareBatch(*service, observers, 0,
                            fixedOccurrence.schedule->payloadSchema(),
                            fixedOccurrence.value,
                            integrationState, 1, 10);
  require(value(fixed.batch, "x").extents ==
              std::vector<std::size_t>({3}),
          "fixed proof geometry changed shape");
  requireCoordinatesFinite(fixed.batch);
  compareField(fixed.batch, *comparisonFields,
               ownedState.view(fixedOccurrence.value.scheduledTime), "sampled-u",
               "u");
  compareField(fixed.batch, *comparisonFields,
               ownedState.view(fixedOccurrence.value.scheduledTime), "sampled-ssh",
               "ssh");

  const auto lineRecord = occurrenceSchedule(3.5, 3, 1, 3.0, 3.0, 1.0);
  const auto lineOccurrence = occurrence(catalog, lineRecord);
  auto line = prepareBatch(*service, observers, 1,
                           lineOccurrence.schedule->payloadSchema(),
                           lineOccurrence.value,
                           integrationState, 2, 11);
  require(value(line.batch, "x").extents ==
              std::vector<std::size_t>({4}),
          "event-variable line did not use the integer payload");
  require(value(line.batch, "payload-real").real64Data()[0] == 3.5 &&
              value(line.batch, "payload-integer").integer64Data()[0] == 3 &&
              value(line.batch, "payload-boolean").boolean8Data()[0] == 1,
          "typed payload values did not reach occurrence storage");
  const auto &lineSampleTimes = value(line.batch, "sample-time");
  require(lineSampleTimes.elementCount() > 1 &&
              lineSampleTimes.real64Data()[0] !=
                  lineSampleTimes.real64Data()[lineSampleTimes.elementCount() -
                                               1],
          "per-sample times were not retained as occurrence metadata");
  requireCoordinatesFinite(line.batch);
  compareField(line.batch, *comparisonFields,
               ownedState.view(lineOccurrence.value.scheduledTime), "sampled-u",
               "u");
  compareField(line.batch, *comparisonFields,
               ownedState.view(lineOccurrence.value.scheduledTime), "sampled-ssh",
               "ssh");

  const auto depthRecord = occurrenceSchedule(2.0, 1, 1, 4.0, 4.0, 1.0);
  const auto depthOccurrence = occurrence(catalog, depthRecord);
  auto depth = prepareBatch(*service, observers, 2,
                            depthOccurrence.schedule->payloadSchema(),
                            depthOccurrence.value,
                            integrationState, 3, 12);
  require(value(depth.batch, "x").extents ==
              std::vector<std::size_t>({3, 2}),
          "fixed-depth proof did not preserve multidimensional extents");
  requireCoordinatesFinite(depth.batch);
  compareField(depth.batch, *comparisonFields,
               ownedState.view(depthOccurrence.value.scheduledTime),
               "sampled-rho-e", "rho_e");
  compareField(depth.batch, *comparisonFields,
               ownedState.view(depthOccurrence.value.scheduledTime), "sampled-ssu",
               "ssu", WVPositionInterpolation::spline);

  const auto stateRecord = occurrenceSchedule(4.0, 2, 1, 5.0, 5.0, 1.0);
  const auto stateOccurrence = occurrence(catalog, stateRecord);
  auto stateBatch = prepareBatch(
      *service, observers, 3, stateOccurrence.schedule->payloadSchema(),
      stateOccurrence.value,
      integrationState, 4, 13);
  require(value(stateBatch.batch, "x").extents ==
              std::vector<std::size_t>({3}),
          "state-coupled proof did not use event state cardinality");
  require(std::abs(value(stateBatch.batch, "x").real64Data()[0] - 700.2) <=
              1e-12 &&
              std::abs(value(stateBatch.batch, "y").real64Data()[0] - 899.8) <=
                  1e-12,
          "state-coupled proof did not read the event state");
  requireCoordinatesFinite(stateBatch.batch);
  compareField(stateBatch.batch, *comparisonFields,
               ownedState.view(stateOccurrence.value.scheduledTime), "sampled-v",
               "v");
  compareField(stateBatch.batch, *comparisonFields,
               ownedState.view(stateOccurrence.value.scheduledTime), "sampled-ssv",
               "ssv");

  const auto raggedRecord = occurrenceSchedule(1.25, 1, 1, 6.0, 6.0, 1.0);
  const auto raggedOccurrence = occurrence(catalog, raggedRecord);
  auto ragged = prepareBatch(*service, observers, 4,
                             raggedOccurrence.schedule->payloadSchema(),
                             raggedOccurrence.value,
                             integrationState, 5, 14);
  const auto passCount = value(ragged.batch, "pass-identifier").elementCount();
  const auto profileCount =
      value(ragged.batch, "profile-identifier").elementCount();
  const auto sampleCount =
      value(ragged.batch, "sample-identifier").elementCount();
  require(passCount == 2 &&
              sumIntegers(value(ragged.batch, "profile-count-by-pass")) ==
                  profileCount &&
              sumIntegers(value(ragged.batch, "sample-count-by-profile")) ==
                  sampleCount,
          "nested ragged relationships do not span their child axes");
  requireCoordinatesFinite(ragged.batch);
  compareField(ragged.batch, *comparisonFields,
               ownedState.view(raggedOccurrence.value.scheduledTime), "sampled-u",
               "u");
  compareField(ragged.batch, *comparisonFields,
               ownedState.view(raggedOccurrence.value.scheduledTime), "sampled-ssh",
               "ssh", WVPositionInterpolation::spline);

  const auto emptyRecord = occurrenceSchedule(7.0, 4, 0, 7.0, 7.0, 1.0);
  const auto emptyOccurrence = occurrence(catalog, emptyRecord);
  auto empty = prepareBatch(*service, observers, 1,
                            emptyOccurrence.schedule->payloadSchema(),
                            emptyOccurrence.value,
                            integrationState, 6, 15);
  require(value(empty.batch, "x").elementCount() == 0 &&
              value(empty.batch, "sampled-u").elementCount() == 0 &&
              value(empty.batch, "sampled-ssh").elementCount() == 0,
          "zero-length occurrence retained phantom samples");
  require(value(empty.batch, "payload-boolean").elementCount() == 1 &&
              value(empty.batch, "payload-boolean").boolean8Data()[0] == 0,
          "zero-length occurrence lost its scalar payload");

  for (std::size_t index = 0;
       index < static_cast<std::size_t>(
                   WVTestObservationOccurrenceProviderKind::count);
       ++index) {
    const auto counters = observationOccurrenceProviderCounters(
        static_cast<WVTestObservationOccurrenceProviderKind>(index));
    require(counters.constructionCount == 1,
            "proof provider was not constructed exactly once");
    require(counters.outputPlanCount == 1,
            "proof provider output plan was not resolved once");
    const auto expectedOccurrences = index == 1 ? 2U : 1U;
    require(counters.occurrencePreparationCount == expectedOccurrences,
            "proof provider occurrence count differs from event count");
    require(counters.batchCount == expectedOccurrences,
            "proof provider batch count differs from event count");
    require(counters.integrationStageCount == 0,
            "passive proof provider was dispatched during an integration stage");
  }
  require(service->metrics().occurrencePreparationCount == 6,
          "evaluation service occurrence count is incorrect");
  require(service->metrics().fieldEvaluationCount == 6,
          "central fields were not evaluated once per occurrence");
  require(service->metrics().occurrenceWorkspaceRetainedBytes > 0 &&
              service->metrics().occurrenceWorkspaceMaximumLiveBytes > 0,
          "occurrence workspace ownership was not measured");
  const auto scheduleCounters = occurrenceScheduleCounters();
  require(scheduleCounters.constructionCount == 5 &&
              scheduleCounters.peekCount == 5,
          "typed occurrence schedules were not resolved at coarse granularity");
}

void testIntegratedExactDenseAndRetry(bool adaptive) {
  resetObservationOccurrenceProviderCounters();
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  WVOutputPlan plan;
  {
    const auto catalog = extensionCatalog();
    const auto observers = protocolDescriptor(catalog, 1.0);
    auto status = WVOutputPlan::create(observers, catalog, 0.0, 1.0, {}, plan);
    require(static_cast<bool>(status),
            "integrated occurrence output planning failed: " +
                status.message);
    status = WVObserverOutputEvaluationService::create(
        configuration(), false, observers,
        std::make_unique<WVReferenceFFTEngine>(), source);
    require(static_cast<bool>(status),
            "integrated occurrence source construction failed: " +
                status.message);
  }
  WVOutputPlan movedPlan = std::move(plan);
  ZeroIntegrationSystem system(movedPlan.stateLayout());
  MutableState state(system.stateLayout());
  std::unique_ptr<WVTimeIntegrator> integrator;
  if (adaptive) {
    WVAdaptiveRK23Options options;
    options.maximumStepSize = 0.5;
    integrator = std::make_unique<WVAdaptiveRK23>(system, options);
  } else {
    integrator = std::make_unique<WVFixedStepRK4>(
        system, WVFixedStepRK4Options{true});
  }
  auto status = integrator->prepareStateAfterRestart(state.view);
  require(static_cast<bool>(status),
          "integrated occurrence restart preparation failed");
  OccurrenceCachingSink sink(*source);
  if (!adaptive)
    sink.failAt(4);
  WVOutputDriver driver(*integrator, movedPlan);
  status = driver.advanceToTime(state.view, 1.0, 0.5, sink);
  if (!adaptive) {
    require(!status && driver.hasPendingDelivery() && state.view.waveVortex.t == 0.0,
            "later-route occurrence failure did not retain the exact event");
    sink.disableFailure();
    status = driver.advanceToTime(state.view, 1.0, 0.5, sink);
  }
  require(static_cast<bool>(status) && !driver.hasPendingDelivery() &&
              state.view.waveVortex.t == 1.0,
          "integrated occurrence run did not reach its final time");

  const auto &driverMetrics = driver.metrics();
  const auto &sourceMetrics = source->metrics();
  require(driverMetrics.initialStateEventCount == 1 &&
              driverMetrics.interpolatedStateEvaluationCount == 2 &&
              driverMetrics.acceptedEndpointStateEventCount == 2 &&
              driverMetrics.generatedSemanticOccurrenceCount == 8 &&
              driverMetrics.committedDeliveryCount == 16,
          "fixed/adaptive occurrence delivery did not span initial, dense, "
          "and exact events");
  require(sourceMetrics.occurrencePreparationCount == 8 &&
              sourceMetrics.occurrenceReuseCount == 8 &&
              sourceMetrics.occurrenceBatchBuildCount == 8 &&
              sourceMetrics.fieldEvaluationCount == 5 &&
              sink.batchBuildCount() == 8 && sink.batchReuseCount() >= 8 &&
              sink.successfulRouteCount() == 16,
          "compatible destination routes did not reuse one exact occurrence "
          "workspace, field result, and batch");
  require(sink.sawDistinctSampleTime(),
          "per-sample time metadata did not remain distinct from the single "
          "trigger-state time");
  if (!adaptive)
    require(sink.firstEventAttempts() ==
                    std::array<std::array<std::size_t, 2>, 2>{{
                        {{1, 1}}, {{1, 2}}}} &&
                driverMetrics.deliveryAttemptCount == 17 &&
                driverMetrics.failureCount == 1,
            "later-route retry repeated a successful route or regenerated "
            "the occurrence");
  const auto counters = observationOccurrenceProviderCounters(
      WVTestObservationOccurrenceProviderKind::
          eventVariableOneDimensionalGeometry);
  require(counters.occurrencePreparationCount == 5 &&
              counters.batchCount == 5 &&
              counters.integrationStageCount == 0,
          "the integrated proof provider was called outside event/batch "
          "granularity");
}

void testIntegratedStateCoupledGeometry(bool adaptive) {
  resetObservationOccurrenceProviderCounters();
  const auto config = configuration();
  const auto catalog = extensionCatalog();
  const auto observers = stateCoupledProtocolDescriptor(catalog, 1.0);
  WVOutputPlan plan;
  auto status =
      WVOutputPlan::create(observers, catalog, 0.0, 1.0, {}, plan);
  require(static_cast<bool>(status),
          "state-coupled integration output planning failed: " +
              status.message);

  WVFrozenForcingSchedule forcing;
  std::unique_ptr<WVConstantStratificationIntegrationSystem> system;
  status = WVConstantStratificationIntegrationSystem::create(
      config, forcing, observers, catalog,
      std::make_unique<WVReferenceFFTEngine>(), system);
  require(static_cast<bool>(status) && system != nullptr &&
              system->particles().size() == 1 &&
              system->particles().front().record().typeIdentifier ==
                  stateCoupledIrregularGeometryType,
          "behavior-named state-coupled observer did not bind as advected "
          "positions: " +
              status.message);

  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      config, false, observers, std::make_unique<WVReferenceFFTEngine>(),
      source);
  require(static_cast<bool>(status),
          "state-coupled occurrence source construction failed: " +
              status.message);

  MutableState state(system->stateLayout());
  status = system->initializeParticleState(state.view);
  require(static_cast<bool>(status),
          "state-coupled position initialization failed: " + status.message);
  std::unique_ptr<WVTimeIntegrator> integrator;
  if (adaptive) {
    WVAdaptiveRK23Options options;
    options.maximumStepSize = 0.5;
    integrator = std::make_unique<WVAdaptiveRK23>(*system, options);
  } else {
    integrator = std::make_unique<WVFixedStepRK4>(
        *system, WVFixedStepRK4Options{true});
  }
  status = integrator->prepareStateAfterRestart(state.view);
  require(static_cast<bool>(status),
          "state-coupled integration restart preparation failed: " +
              status.message);

  StateCoupledOccurrenceSink sink(*source);
  WVOutputDriver driver(*integrator, plan);
  status = driver.advanceToTime(state.view, 1.0, 0.5, sink);
  require(static_cast<bool>(status) && state.view.waveVortex.t == 1.0 &&
              !driver.hasPendingDelivery(),
          "state-coupled fixed/adaptive output run did not finish: " +
              status.message);

  const auto &snapshots = sink.snapshots();
  require(snapshots.size() == 5 &&
              snapshots.front().scheduledTime == 0.0 &&
              snapshots.front().kind == WVOutputEventKind::initial &&
              snapshots.back().scheduledTime == 1.0 &&
              snapshots.back().kind == WVOutputEventKind::acceptedEndpoint,
          "state-coupled output did not span the complete event lattice");
  const bool sawInterpolated = std::any_of(
      snapshots.begin(), snapshots.end(), [](const auto &snapshot) {
        return snapshot.kind == WVOutputEventKind::interpolated;
      });
  require(sawInterpolated,
          "state-coupled output never consumed an integrated dense state");

  for (const auto &snapshot : snapshots) {
    const auto ordinal = static_cast<std::size_t>(
        std::llround(snapshot.scheduledTime / 0.25));
    const bool emitsSamples = (ordinal & 1U) == 0;
    const auto expectedCount =
        emitsSamples ? 1U + (3U + ordinal) % 4U : 0U;
    const auto coordinateOffset = 0.05 * (4.0 + 0.5 * ordinal);
    require(snapshot.stateX.size() == 4 &&
                snapshot.stateY.size() == snapshot.stateX.size() &&
                snapshot.stateZ.size() == snapshot.stateX.size() &&
                snapshot.batchX.size() == expectedCount &&
                snapshot.batchY.size() == expectedCount &&
                snapshot.batchZ.size() == expectedCount &&
                snapshot.sampledV.size() == expectedCount &&
                snapshot.sampledSsv.size() == expectedCount,
            "state-coupled occurrence lost an advected position or field "
            "sample");
    for (std::size_t index = 0; index < expectedCount; ++index) {
      require(std::abs(snapshot.batchX[index] -
                       (snapshot.stateX[index] + coordinateOffset)) <= 1e-12 &&
                  std::abs(snapshot.batchY[index] -
                           (snapshot.stateY[index] - coordinateOffset)) <=
                      1e-12 &&
                  snapshot.batchZ[index] == snapshot.stateZ[index] &&
                  std::isfinite(snapshot.sampledV[index]) &&
                  std::isfinite(snapshot.sampledSsv[index]),
              "source-built occurrence geometry did not come from the "
              "driver event state");
    }
  }

  const auto configured = stateCoupledIrregularGeometryRecord();
  require(snapshots.front().stateX == configured.x &&
              snapshots.front().stateY == configured.y &&
              snapshots.front().stateZ == configured.z,
          "bindIntegration position initialization did not reach the initial "
          "occurrence");
  bool moved = false;
  for (std::size_t index = 0; index < snapshots.front().stateX.size(); ++index)
    moved = moved ||
            std::abs(snapshots.back().stateX[index] -
                     snapshots.front().stateX[index]) > 1e-12 ||
            std::abs(snapshots.back().stateY[index] -
                     snapshots.front().stateY[index]) > 1e-12 ||
            std::abs(snapshots.back().stateZ[index] -
                     snapshots.front().stateZ[index]) > 1e-12;
  require(moved && system->metrics().rightHandSideEvaluationCount > 0 &&
              system->metrics().velocityFieldEvaluationCount > 0 &&
              system->metrics().particleValueWriteCount > 0,
          "bound state-coupled positions were not advanced by the real "
          "integration system");

  const auto finalBlock = [&](const std::string &identifier) {
    for (std::size_t block = 0; block < state.view.additionalBlockCount;
         ++block)
      if (state.view.additionalBlocks[block].layout != nullptr &&
          state.view.additionalBlocks[block].layout->identifier == identifier)
        return state.view.additionalBlocks + block;
    return static_cast<WVAdditionalStateBlockView *>(nullptr);
  };
  const auto *finalX = finalBlock(configured.stateBlockIdentifiers[0]);
  const auto *finalY = finalBlock(configured.stateBlockIdentifiers[1]);
  const auto *finalZ = finalBlock(configured.stateBlockIdentifiers[2]);
  require(finalX != nullptr && finalY != nullptr && finalZ != nullptr,
          "integrated state lost a bound coordinate block");
  for (std::size_t index = 0; index < snapshots.back().stateX.size(); ++index)
    require(snapshots.back().stateX[index] == finalX->realData[index] &&
                snapshots.back().stateY[index] == finalY->realData[index] &&
                snapshots.back().stateZ[index] == finalZ->realData[index],
            "final source batch did not use the accepted integrated position "
            "state");

  const auto counters = observationOccurrenceProviderCounters(
      WVTestObservationOccurrenceProviderKind::
          stateCoupledIrregularGeometry);
  require(counters.occurrencePreparationCount == snapshots.size() &&
              counters.batchCount == snapshots.size(),
          "state-coupled provider did not execute once per routed occurrence");
}

void testNestedRaggedNetCDFPersistenceAndAppend() {
  resetObservationOccurrenceProviderCounters();
  OccurrenceTemporaryDirectory directory;
  const auto path = directory.path / "nested-ragged.nc";
  const auto config = configuration();
  const auto catalog = extensionCatalog();
  const auto observers =
      nestedRaggedNetCDFDescriptor(catalog, path, 1.0);
  const auto checkpoint = checkpointTemplate(config);
  WVIntegrationStateLayout layout;
  auto status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, observers, layout);
  require(static_cast<bool>(status),
          "nested-ragged output state layout failed: " + status.message);
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      config, false, observers, std::make_unique<WVReferenceFFTEngine>(),
      source);
  require(static_cast<bool>(status),
          "nested-ragged output source construction failed: " +
              status.message);

  const auto deliverPlan = [&](WVModelOutputNetCDFSink &sink,
                               const WVOutputPlan &plan) {
    auto deliveryStatus = sink.preflight(plan);
    require(static_cast<bool>(deliveryStatus),
            "nested-ragged sink preflight failed: " +
                deliveryStatus.message);
    for (std::size_t eventIndex = 0; eventIndex < plan.eventCount();
         ++eventIndex) {
      const auto planned = plan.event(eventIndex);
      auto waveVortex = checkpoint.state.view();
      waveVortex.t = planned.scheduledTime;
      WVOutputEvent event;
      event.eventOrdinal = planned.eventOrdinal;
      event.scheduledTime = planned.scheduledTime;
      event.kind = planned.scheduledTime == 0.0
                       ? WVOutputEventKind::initial
                       : WVOutputEventKind::acceptedEndpoint;
      event.state.waveVortex = waveVortex;
      event.routes = planned.routes;
      event.routeCount = planned.routeCount;
      for (std::size_t route = 0; route < planned.routeCount; ++route) {
        WVOutputDeliveryResult delivery;
        deliveryStatus = sink.deliver(event, planned.routes[route], delivery);
        require(static_cast<bool>(deliveryStatus),
                "nested-ragged delivery failed: " +
                    deliveryStatus.message);
      }
    }
  };

  WVOutputPlan firstPlan;
  status = WVOutputPlan::create(observers, catalog, 0.0, 0.5, {}, firstPlan);
  require(static_cast<bool>(status) && firstPlan.eventCount() == 3,
          "nested-ragged create plan has the wrong event lattice");
  WVModelOutputNetCDFConfiguration outputConfiguration{catalog, checkpoint,
                                                        false};
  WVModelOutputNetCDFSink sink;
  auto persistence = WVModelOutputNetCDFSink::createNew(
      outputConfiguration, observers, firstPlan, layout, source.get(), sink);
  require(static_cast<bool>(persistence),
          "nested-ragged output creation failed: " + persistence.message);
  deliverPlan(sink, firstPlan);
  persistence = sink.close();
  require(static_cast<bool>(persistence),
          "nested-ragged output close failed: " + persistence.message);
  auto counters = observationOccurrenceProviderCounters(
      WVTestObservationOccurrenceProviderKind::nestedRaggedGeometry);
  require(counters.occurrencePreparationCount == 3 &&
              counters.batchCount == 4,
          "nested-ragged create phase regenerated an occurrence batch");
  const auto firstPersistence = readNestedRaggedPersistence(path);

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {path.string()}, *catalog, inspection);
  require(static_cast<bool>(persistence),
          "nested-ragged output inspection failed: " + persistence.message);
  WVOutputPlan appendPlan;
  status = WVOutputPlan::create(observers, catalog, 0.5, 1.0,
                                inspection.scheduleContinuations,
                                appendPlan);
  require(static_cast<bool>(status) && appendPlan.eventCount() == 2,
          "nested-ragged append plan repeated committed occurrences");
  WVModelOutputNetCDFSink append;
  persistence = WVModelOutputNetCDFSink::openAppend(
      outputConfiguration, observers, appendPlan, layout, source.get(),
      inspection.destinationProgress, append);
  require(static_cast<bool>(persistence),
          "nested-ragged append open failed: " + persistence.message);
  deliverPlan(append, appendPlan);
  persistence = append.close();
  require(static_cast<bool>(persistence),
          "nested-ragged append close failed: " + persistence.message);
  counters = observationOccurrenceProviderCounters(
      WVTestObservationOccurrenceProviderKind::nestedRaggedGeometry);
  require(counters.occurrencePreparationCount == 5 &&
              counters.batchCount == 6,
          "nested-ragged append regenerated committed occurrence batches");

  const auto persisted = readNestedRaggedPersistence(path);
  const bool preservedPrefix =
      isPrefix(firstPersistence.times, persisted.times) &&
      isPrefix(firstPersistence.scheduleOrdinals,
               persisted.scheduleOrdinals) &&
      isPrefix(firstPersistence.committedPass, persisted.committedPass) &&
      isPrefix(firstPersistence.committedProfile,
               persisted.committedProfile) &&
      isPrefix(firstPersistence.committedSample,
               persisted.committedSample) &&
      isPrefix(firstPersistence.profileCounts, persisted.profileCounts) &&
      isPrefix(firstPersistence.sampleCounts, persisted.sampleCounts) &&
      isPrefix(firstPersistence.passIdentifiers,
               persisted.passIdentifiers) &&
      isPrefix(firstPersistence.profileIdentifiers,
               persisted.profileIdentifiers) &&
      isPrefix(firstPersistence.sampleIdentifiers,
               persisted.sampleIdentifiers) &&
      isPrefix(firstPersistence.sampleTimes, persisted.sampleTimes) &&
      isPrefix(firstPersistence.x, persisted.x) &&
      isPrefix(firstPersistence.y, persisted.y) &&
      isPrefix(firstPersistence.z, persisted.z) &&
      isPrefix(firstPersistence.directions, persisted.directions) &&
      isPrefix(firstPersistence.payloadReal, persisted.payloadReal) &&
      isPrefix(firstPersistence.payloadInteger, persisted.payloadInteger) &&
      isPrefix(firstPersistence.payloadBoolean, persisted.payloadBoolean) &&
      isPrefix(firstPersistence.sampledU, persisted.sampledU) &&
      isPrefix(firstPersistence.sampledSsh, persisted.sampledSsh);
  require(preservedPrefix,
          "nested-ragged append rewrote an already committed value or offset");

  NestedRaggedPersistence expected;
  for (std::size_t ordinal = 0; ordinal < 5; ++ordinal) {
    const auto scheduledTime = 0.25 * static_cast<double>(ordinal);
    const auto payloadReal =
        std::fma(static_cast<double>(ordinal), 0.5, 4.0);
    const auto payloadInteger = 1LL + static_cast<long long>(ordinal);
    const auto payloadBoolean = static_cast<unsigned char>(1U ^ (ordinal & 1U));
    expected.times.push_back(scheduledTime);
    expected.scheduleOrdinals.push_back(static_cast<long long>(ordinal));
    expected.payloadReal.push_back(payloadReal);
    expected.payloadInteger.push_back(payloadInteger);
    expected.payloadBoolean.push_back(payloadBoolean);

    const auto passCount =
        payloadBoolean == 0 ? 0U
                            : 1U + static_cast<std::size_t>(payloadInteger) % 2U;
    std::vector<long long> occurrenceProfileCounts(passCount);
    std::size_t profileCount = 0;
    for (std::size_t pass = 0; pass < passCount; ++pass) {
      occurrenceProfileCounts[pass] = static_cast<long long>(
          1U + (pass + static_cast<std::size_t>(payloadInteger)) % 2U);
      profileCount +=
          static_cast<std::size_t>(occurrenceProfileCounts[pass]);
      expected.passIdentifiers.push_back(static_cast<long long>(pass + 1));
    }
    expected.profileCounts.insert(expected.profileCounts.end(),
                                  occurrenceProfileCounts.begin(),
                                  occurrenceProfileCounts.end());
    std::size_t occurrenceSampleCount = 0;
    for (std::size_t profile = 0; profile < profileCount; ++profile) {
      const auto sampleCount = static_cast<std::size_t>(
          1U + (profile + static_cast<std::size_t>(payloadInteger)) % 3U);
      expected.sampleCounts.push_back(static_cast<long long>(sampleCount));
      expected.profileIdentifiers.push_back(
          static_cast<long long>(profile + 1));
      for (std::size_t local = 0; local < sampleCount; ++local) {
        expected.sampleIdentifiers.push_back(
            static_cast<long long>(occurrenceSampleCount + 1));
        expected.sampleTimes.push_back(scheduledTime + 0.03 * local);
        expected.x.push_back(1000.0 + 90.0 * profile + 12.0 * local +
                             0.1 * payloadReal);
        expected.y.push_back(1200.0 + 55.0 * profile + 8.0 * local);
        expected.z.push_back(-300.0 - 80.0 * profile - 10.0 * local);
        expected.directions.push_back(payloadBoolean);
        ++occurrenceSampleCount;
      }
    }
    expected.committedPass.push_back(static_cast<long long>(
        expected.passIdentifiers.size()));
    expected.committedProfile.push_back(static_cast<long long>(
        expected.profileIdentifiers.size()));
    expected.committedSample.push_back(static_cast<long long>(
        expected.sampleIdentifiers.size()));
  }
  expected.passCount = expected.passIdentifiers.size();
  expected.profileCount = expected.profileIdentifiers.size();
  expected.sampleCount = expected.sampleIdentifiers.size();

  require(persisted.passCount == expected.passCount &&
              persisted.profileCount == expected.profileCount &&
              persisted.sampleCount == expected.sampleCount &&
              persisted.times == expected.times &&
              persisted.scheduleOrdinals == expected.scheduleOrdinals &&
              persisted.committedPass == expected.committedPass &&
              persisted.committedProfile == expected.committedProfile &&
              persisted.committedSample == expected.committedSample &&
              persisted.profileCounts == expected.profileCounts &&
              persisted.sampleCounts == expected.sampleCounts &&
              persisted.passIdentifiers == expected.passIdentifiers &&
              persisted.profileIdentifiers == expected.profileIdentifiers &&
              persisted.sampleIdentifiers == expected.sampleIdentifiers &&
              persisted.sampleTimes == expected.sampleTimes &&
              persisted.x == expected.x && persisted.y == expected.y &&
              persisted.z == expected.z &&
              persisted.directions == expected.directions &&
              persisted.payloadReal == expected.payloadReal &&
              persisted.payloadInteger == expected.payloadInteger &&
              persisted.payloadBoolean == expected.payloadBoolean,
          "persisted nested pass/profile/sample counts, offsets, or values "
          "differ from the provider occurrences");
  require(persisted.committedPass[1] == persisted.committedPass[0] &&
              persisted.committedProfile[1] ==
                  persisted.committedProfile[0] &&
              persisted.committedSample[1] ==
                  persisted.committedSample[0] &&
              persisted.committedPass[3] == persisted.committedPass[2] &&
              persisted.committedProfile[3] ==
                  persisted.committedProfile[2] &&
              persisted.committedSample[3] ==
                  persisted.committedSample[2],
          "zero-length nested occurrences advanced a ragged-axis offset");
  const auto fieldsArePersisted = [](const std::vector<double> &values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
      return std::isfinite(value) && value != NC_FILL_DOUBLE;
    });
  };
  require(persisted.sampledU.size() == expected.sampleCount &&
              persisted.sampledSsh.size() == expected.sampleCount &&
              fieldsArePersisted(persisted.sampledU) &&
              fieldsArePersisted(persisted.sampledSsh),
          "nested-ragged sampled fields were not persisted for every sample");
}

void testOccurrenceCompatibilityIdentity() {
  resetObservationOccurrenceProviderCounters();
  const auto catalog = extensionCatalog();
  const auto observers = protocolDescriptor(catalog, 0.0);
  std::unique_ptr<WVObserverOutputEvaluationService> service;
  auto status = WVObserverOutputEvaluationService::create(
      configuration(), false, observers,
      std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status),
          "occurrence identity source construction failed");
  const auto stateStorage = state(configuration());
  const auto sameSchedule = occurrenceSchedule(3.5, 3, 1, 0.0);
  const auto differentSchedule = occurrenceSchedule(9.5, 3, 1, 0.0);
  const auto same = occurrence(catalog, sameSchedule);
  const auto different = occurrence(catalog, differentSchedule);
  const auto groupRecord = [](std::string identifier,
                              WVOutputScheduleRecord schedule) {
    WVOutputGroupRecord group;
    group.identifier = std::move(identifier);
    group.name = group.identifier;
    group.schedule = std::move(schedule);
    group.observerIdentifiers = {"event-variable-line"};
    return group;
  };
  const auto sameGroup = groupRecord("same-occurrence", sameSchedule);
  const auto differentPayloadGroup =
      groupRecord("same-occurrence", differentSchedule);
  const auto otherLogicalGroup =
      groupRecord("other-occurrence", sameSchedule);
  const auto &record = observers.observers().at(1);
  const auto *resolved = observers.resolvedObserver(record);
  require(resolved != nullptr, "occurrence identity observer is unresolved");
  WVOutputObserverView observer{1, &record, resolved};
  const auto route = [&](const ResolvedOccurrence &item,
                         const WVOutputGroupRecord &group,
                         std::size_t semanticOrdinal,
                         std::size_t fileOrdinal) {
    WVOutputRouteView result;
    result.fileOrdinal = fileOrdinal;
    result.scheduleOrdinal = item.value.ordinal;
    result.observers = &observer;
    result.observerCount = 1;
    result.semanticScheduleOrdinal = semanticOrdinal;
    result.proposedScheduleCursor = &item.value.proposedCursor.values;
    result.schedulePayloadSchema = &item.schedule->payloadSchema();
    result.schedulePayload = &item.value.payload;
    result.scheduleCursorIdentity = item.value.cursorIdentity;
    result.semanticScheduleRecord = &group;
    return result;
  };
  std::array<WVOutputRouteView, 2> routes{{route(same, sameGroup, 20, 0),
                                          route(same, sameGroup, 20, 1)}};
  WVOutputEvent event;
  event.eventOrdinal = 7;
  event.scheduledTime = 0.0;
  event.state.waveVortex = stateStorage.view(0.0);
  event.routes = routes.data();
  event.routeCount = routes.size();
  status = service->prepare(event);
  require(static_cast<bool>(status),
          "compatible occurrence preparation failed");
  WVObservationOccurrenceIdentity first;
  WVObservationOccurrenceIdentity second;
  require(static_cast<bool>(service->preparedOccurrenceIdentity(
              routes[0], observer, first)) &&
              static_cast<bool>(service->preparedOccurrenceIdentity(
                  routes[1], observer, second)) &&
              sameObservationOccurrenceIdentity(first, second) &&
              samePreparedObservationOccurrenceIdentity(first, second) &&
              service->metrics().occurrencePreparationCount == 1 &&
              service->metrics().occurrenceReuseCount == 1,
          "compatible destinations received different semantic occurrence "
          "identities");
  auto diagnosticOnly = second;
  diagnosticOnly.observerOrdinal += 10;
  diagnosticOnly.semanticScheduleOrdinal += 10;
  diagnosticOnly.scheduleCursorIdentity ^= 0x1;
  diagnosticOnly.payloadFingerprint ^= 0x2;
  diagnosticOnly.geometryFingerprint ^= 0x4;
  diagnosticOnly.fieldPlanFingerprint ^= 0x8;
  require(sameObservationOccurrenceIdentity(first, diagnosticOnly),
          "plan-local ordinals or diagnostic fingerprints became semantic "
          "occurrence keys");

  routes[1] = route(different, differentPayloadGroup, 20, 1);
  status = service->prepare(event);
  require(static_cast<bool>(status),
          "payload-isolated occurrence preparation failed");
  require(static_cast<bool>(service->preparedOccurrenceIdentity(
              routes[0], observer, first)) &&
              static_cast<bool>(service->preparedOccurrenceIdentity(
                  routes[1], observer, second)) &&
              !sameObservationOccurrenceIdentity(first, second) &&
              !samePreparedObservationOccurrenceIdentity(first, second) &&
              first.payloadFingerprint != second.payloadFingerprint &&
              first.geometryFingerprint != second.geometryFingerprint &&
              service->metrics().occurrencePreparationCount == 3,
          "different payload/cursor geometry shared a prepared occurrence");

  routes = {{route(same, sameGroup, 30, 0),
             route(same, otherLogicalGroup, 31, 1)}};
  event.eventOrdinal = 8;
  event.routes = routes.data();
  status = service->prepare(event);
  require(static_cast<bool>(status),
          "logical-schedule-isolated occurrence preparation failed");
  require(static_cast<bool>(service->preparedOccurrenceIdentity(
              routes[0], observer, first)) &&
              static_cast<bool>(service->preparedOccurrenceIdentity(
                  routes[1], observer, second)) &&
              !sameObservationOccurrenceIdentity(first, second) &&
              !samePreparedObservationOccurrenceIdentity(first, second) &&
              first.semanticScheduleOrdinal !=
                  second.semanticScheduleOrdinal &&
              service->metrics().occurrencePreparationCount == 5,
          "coincident logical schedule instances shared an occurrence");
}

void testSegmentedOccurrenceIdentity() {
  const auto catalog = extensionCatalog();
  const auto observers = protocolDescriptor(catalog, 1.0);
  WVOutputPlan uninterrupted;
  auto status =
      WVOutputPlan::create(observers, catalog, 0.0, 1.0, {}, uninterrupted);
  require(static_cast<bool>(status),
          "uninterrupted occurrence planning failed");

  auto continuations = uninterrupted.initialContinuations();
  for (std::size_t eventIndex = 0;
       eventIndex < uninterrupted.eventCount(); ++eventIndex) {
    const auto event = uninterrupted.event(eventIndex);
    if (event.scheduledTime > 0.5)
      break;
    for (std::size_t routeIndex = 0; routeIndex < event.routeCount;
         ++routeIndex) {
      const auto &route = event.routes[routeIndex];
      const auto found = std::find_if(
          continuations.begin(), continuations.end(), [&](const auto &item) {
            return item.fileIdentifier == route.fileIdentifier &&
                   item.groupIdentifier == route.groupIdentifier;
          });
      require(found != continuations.end() &&
                  route.proposedScheduleCursor != nullptr,
              "segment continuation route is incomplete");
      found->cursor.committedOrdinal = route.scheduleOrdinal;
      found->cursor.values = *route.proposedScheduleCursor;
    }
  }

  WVPortableObserverDescriptor resumedObservers;
  status = WVPortableObserverDescriptor::create(
      observers.record(), catalog, resumedObservers);
  require(static_cast<bool>(status),
          "resumed observer descriptor recompilation failed");
  WVOutputPlan resumed;
  status = WVOutputPlan::create(resumedObservers, catalog, 0.5, 1.0,
                                continuations, resumed);
  require(static_cast<bool>(status) && resumed.eventCount() > 0,
          "resumed occurrence planning failed");
  const auto resumedEvent = resumed.event(0);
  require(resumedEvent.scheduledTime == 0.75,
          "resumed occurrence duplicated the committed boundary");

  WVOutputPlannedEventView uninterruptedEvent;
  bool foundEvent = false;
  for (std::size_t eventIndex = 0;
       eventIndex < uninterrupted.eventCount(); ++eventIndex) {
    const auto candidate = uninterrupted.event(eventIndex);
    if (candidate.scheduledTime == resumedEvent.scheduledTime) {
      uninterruptedEvent = candidate;
      foundEvent = true;
      break;
    }
  }
  require(foundEvent,
          "uninterrupted plan omitted the resumed occurrence");
  const auto occurrenceRoute = [](const WVOutputPlannedEventView &event) {
    const WVOutputRouteView *result = nullptr;
    for (std::size_t index = 0; index < event.routeCount; ++index)
      if (event.routes[index].fileIdentifier == "primary" &&
          event.routes[index].groupIdentifier == "shared-occurrence") {
        result = event.routes + index;
        break;
      }
    return result;
  };
  const auto *uninterruptedRoute = occurrenceRoute(uninterruptedEvent);
  const auto *resumedRoute = occurrenceRoute(resumedEvent);
  require(uninterruptedRoute != nullptr && resumedRoute != nullptr &&
              uninterruptedRoute->observerCount == 1 &&
              resumedRoute->observerCount == 1 &&
              uninterruptedRoute->semanticScheduleRecord !=
                  resumedRoute->semanticScheduleRecord &&
              uninterruptedRoute->observers[0].record !=
                  resumedRoute->observers[0].record &&
              uninterruptedRoute->semanticScheduleOrdinal ==
                  resumedRoute->semanticScheduleOrdinal &&
              uninterruptedRoute->scheduleOrdinal ==
                  resumedRoute->scheduleOrdinal &&
              uninterruptedRoute->scheduleCursorIdentity ==
                  resumedRoute->scheduleCursorIdentity &&
              uninterruptedRoute->schedulePayload != nullptr &&
              resumedRoute->schedulePayload != nullptr &&
              uninterruptedRoute->schedulePayload->sameValue(
                  *resumedRoute->schedulePayload),
          "segmented occurrence changed its schedule, cursor, or payload "
          "identity");

  std::unique_ptr<WVObserverOutputEvaluationService> uninterruptedSource;
  std::unique_ptr<WVObserverOutputEvaluationService> resumedSource;
  status = WVObserverOutputEvaluationService::create(
      configuration(), false, observers,
      std::make_unique<WVReferenceFFTEngine>(), uninterruptedSource);
  if (status)
    status = WVObserverOutputEvaluationService::create(
        configuration(), false, resumedObservers,
        std::make_unique<WVReferenceFFTEngine>(), resumedSource);
  require(static_cast<bool>(status),
          "segmented occurrence source construction failed");
  status = uninterruptedSource->preflight(uninterrupted);
  if (status)
    status = resumedSource->preflight(resumed);
  require(static_cast<bool>(status),
          "segmented occurrence preflight failed");
  const auto stateStorage = state(configuration());
  const auto prepareIdentity = [&](WVObserverOutputEvaluationService &source,
                                   const WVOutputPlannedEventView &planned,
                                   const WVOutputRouteView &route,
                                   WVObservationOccurrenceIdentity &identity) {
    WVOutputEvent event;
    event.eventOrdinal = planned.eventOrdinal;
    event.scheduledTime = planned.scheduledTime;
    event.state.waveVortex = stateStorage.view(planned.scheduledTime);
    event.routes = planned.routes;
    event.routeCount = planned.routeCount;
    auto prepareStatus = source.prepare(event);
    require(static_cast<bool>(prepareStatus),
            "segmented occurrence preparation failed");
    require(route.observerCount == 1,
            "segmented occurrence route has an unexpected observer count");
    prepareStatus = source.preparedOccurrenceIdentity(
        route, route.observers[0], identity);
    require(static_cast<bool>(prepareStatus),
            "segmented prepared-occurrence identity lookup failed");
  };
  WVObservationOccurrenceIdentity uninterruptedIdentity;
  WVObservationOccurrenceIdentity resumedIdentity;
  prepareIdentity(*uninterruptedSource, uninterruptedEvent, *uninterruptedRoute,
                  uninterruptedIdentity);
  prepareIdentity(*resumedSource, resumedEvent, *resumedRoute, resumedIdentity);
  require(sameObservationOccurrenceIdentity(uninterruptedIdentity,
                                             resumedIdentity),
          "restart/segment reconstruction changed semantic occurrence "
          "identity");

  const auto eventAt = [](const WVOutputPlan &plan, double scheduledTime) {
    WVOutputPlannedEventView result;
    for (std::size_t index = 0; index < plan.eventCount(); ++index) {
      const auto candidate = plan.event(index);
      if (candidate.scheduledTime == scheduledTime) {
        result = candidate;
        break;
      }
    }
    return result;
  };

  // A plan-local semantic ordinal and provider fingerprints are not an exact
  // schedule identity. A different construction schedule can emit the same
  // observer, time, ordinal, cursor, payload, geometry, and field plan.
  const auto differentScheduleObservers = protocolDescriptor(catalog, 2.0);
  WVOutputPlan differentSchedulePlan;
  status = WVOutputPlan::create(differentScheduleObservers, catalog, 0.0, 2.0,
                                {}, differentSchedulePlan);
  require(static_cast<bool>(status),
          "different-schedule occurrence planning failed");
  const auto differentScheduleEvent = eventAt(differentSchedulePlan, 0.75);
  const auto *differentScheduleRoute = occurrenceRoute(differentScheduleEvent);
  require(differentScheduleRoute != nullptr,
          "different-schedule plan omitted the comparison occurrence");
  std::unique_ptr<WVObserverOutputEvaluationService> differentScheduleSource;
  status = WVObserverOutputEvaluationService::create(
      configuration(), false, differentScheduleObservers,
      std::make_unique<WVReferenceFFTEngine>(), differentScheduleSource);
  if (status)
    status = differentScheduleSource->preflight(differentSchedulePlan);
  require(static_cast<bool>(status),
          "different-schedule occurrence source construction failed");
  WVObservationOccurrenceIdentity differentScheduleIdentity;
  prepareIdentity(*differentScheduleSource, differentScheduleEvent,
                  *differentScheduleRoute, differentScheduleIdentity);
  require(uninterruptedIdentity.observerOrdinal ==
                  differentScheduleIdentity.observerOrdinal &&
              uninterruptedIdentity.semanticScheduleOrdinal ==
                  differentScheduleIdentity.semanticScheduleOrdinal &&
              uninterruptedIdentity.scheduleOrdinal ==
                  differentScheduleIdentity.scheduleOrdinal &&
              uninterruptedIdentity.scheduleCursorIdentity ==
                  differentScheduleIdentity.scheduleCursorIdentity &&
              uninterruptedIdentity.payloadFingerprint ==
                  differentScheduleIdentity.payloadFingerprint &&
              uninterruptedIdentity.geometryFingerprint ==
                  differentScheduleIdentity.geometryFingerprint &&
              uninterruptedIdentity.fieldPlanFingerprint ==
                  differentScheduleIdentity.fieldPlanFingerprint &&
              !sameObservationOccurrenceIdentity(
                  uninterruptedIdentity, differentScheduleIdentity),
          "different construction schedules aliased through plan-local "
          "ordinals or fingerprints");

  // Resolved observer semantics are likewise value-based across independent
  // descriptors, not the plan-local observer ordinal.
  auto differentObserverRecord = observers.record();
  const auto observerFound = std::find_if(
      differentObserverRecord.observers.begin(),
      differentObserverRecord.observers.end(), [](const auto &observer) {
        return observer.identifier == "event-variable-line";
      });
  require(observerFound != differentObserverRecord.observers.end(),
          "observer-identity comparison record is incomplete");
  observerFound->outputScale = 2.0;
  WVPortableObserverDescriptor differentObserverDescriptor;
  status = WVPortableObserverDescriptor::create(
      differentObserverRecord, catalog, differentObserverDescriptor);
  require(static_cast<bool>(status),
          "different-observer descriptor construction failed");
  WVOutputPlan differentObserverPlan;
  status = WVOutputPlan::create(differentObserverDescriptor, catalog, 0.0, 1.0,
                                {}, differentObserverPlan);
  require(static_cast<bool>(status),
          "different-observer occurrence planning failed");
  const auto differentObserverEvent = eventAt(differentObserverPlan, 0.75);
  const auto *differentObserverRoute = occurrenceRoute(differentObserverEvent);
  require(differentObserverRoute != nullptr,
          "different-observer plan omitted the comparison occurrence");
  std::unique_ptr<WVObserverOutputEvaluationService> differentObserverSource;
  status = WVObserverOutputEvaluationService::create(
      configuration(), false, differentObserverDescriptor,
      std::make_unique<WVReferenceFFTEngine>(), differentObserverSource);
  if (status)
    status = differentObserverSource->preflight(differentObserverPlan);
  require(static_cast<bool>(status),
          "different-observer occurrence source construction failed");
  WVObservationOccurrenceIdentity differentObserverIdentity;
  prepareIdentity(*differentObserverSource, differentObserverEvent,
                  *differentObserverRoute, differentObserverIdentity);
  require(uninterruptedIdentity.observerOrdinal ==
                  differentObserverIdentity.observerOrdinal &&
              uninterruptedIdentity.semanticScheduleOrdinal ==
                  differentObserverIdentity.semanticScheduleOrdinal &&
              uninterruptedIdentity.scheduleOrdinal ==
                  differentObserverIdentity.scheduleOrdinal &&
              uninterruptedIdentity.scheduleCursorIdentity ==
                  differentObserverIdentity.scheduleCursorIdentity &&
              uninterruptedIdentity.payloadFingerprint ==
                  differentObserverIdentity.payloadFingerprint &&
              uninterruptedIdentity.geometryFingerprint ==
                  differentObserverIdentity.geometryFingerprint &&
              uninterruptedIdentity.fieldPlanFingerprint ==
                  differentObserverIdentity.fieldPlanFingerprint &&
              !sameObservationOccurrenceIdentity(
                  uninterruptedIdentity, differentObserverIdentity),
          "different resolved observer semantics aliased through a plan-local "
          "ordinal");
}

void testCoincidentPrimitiveReconstructionReuse() {
  const auto config = configuration();
  const auto catalog = extensionCatalog();
  const auto observers = primitiveReuseDescriptor(catalog);
  WVOutputPlan plan;
  auto status = WVOutputPlan::create(observers, catalog, 0.0, 0.0, {}, plan);
  require(static_cast<bool>(status) && plan.eventCount() == 1,
          "coincident primitive-reuse output planning failed");
  const auto planned = plan.event(0);
  require(planned.routeCount == 3,
          "coincident primitive-reuse event omitted a route");

  std::unique_ptr<WVFieldEvaluationService> fields;
  status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), fields);
  require(static_cast<bool>(status),
          "coincident primitive-reuse field service construction failed");
  std::unique_ptr<WVObserverOutputEvaluationService> source;
  status = WVObserverOutputEvaluationService::create(
      config, false, observers, {}, source, fields.get());
  require(static_cast<bool>(status),
          "coincident primitive-reuse source construction failed");
  status = source->preflight(plan);
  require(static_cast<bool>(status),
          "coincident primitive-reuse preflight failed");

  const auto stateStorage = state(config);
  WVOutputEvent event;
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state.waveVortex = stateStorage.view(planned.scheduledTime);
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  const auto before = fields->metrics();
  status = source->prepare(event);
  require(static_cast<bool>(status),
          "coincident primitive-reuse preparation failed");
  const auto &after = fields->metrics();
  require(after.eventGeometryPreparationCount -
                  before.eventGeometryPreparationCount ==
              2 &&
              after.eventEvaluationCount - before.eventEvaluationCount == 2 &&
              after.linearInterpolationCount -
                      before.linearInterpolationCount ==
                  12 &&
              after.splineInterpolationCount -
                      before.splineInterpolationCount ==
                  12,
          "distinct coincident occurrences did not retain separate geometry "
          "and interpolation work");
  require(after.transformCount - before.transformCount == 1 &&
              after.primitiveFieldEvaluationCount -
                      before.primitiveFieldEvaluationCount ==
                  2,
          "distinct coincident occurrences repeated primitive "
          "reconstruction");
}

struct WindowStorageMetrics {
  std::size_t persistent = 0;
  std::size_t retained = 0;
  std::size_t maximumLive = 0;
  std::size_t batchRetained = 0;
  std::size_t batchMaximumLive = 0;
};

WindowStorageMetrics runOccurrenceWindow(std::size_t eventCount) {
  const auto config = configuration();
  const auto catalog = extensionCatalog();
  const auto observers = protocolDescriptor(catalog, 0.0);
  std::unique_ptr<WVObserverOutputEvaluationService> service;
  auto status = WVObserverOutputEvaluationService::create(
      config, false, observers, std::make_unique<WVReferenceFFTEngine>(),
      service);
  require(static_cast<bool>(status),
          "bounded-window source construction failed");
  WVKernelStatus scheduleStatus;
  const auto schedule = makeOccurrenceSchedule(
      occurrenceSchedule(3.5, 0, 1,
                         static_cast<double>(eventCount - 1), 0.0, 1.0),
      scheduleStatus);
  require(static_cast<bool>(scheduleStatus) && schedule != nullptr,
          "bounded-window schedule construction failed");
  const auto stateStorage = state(config);
  const auto &record = observers.observers().at(1);
  const auto *resolved = observers.resolvedObserver(record);
  require(resolved != nullptr, "bounded-window observer is unresolved");
  WVOutputObserverView observer{1, &record, resolved};
  WVOutputScheduleCursor cursor;
  for (std::size_t index = 0; index < eventCount; ++index) {
    WVOutputScheduleOccurrence occurrenceValue;
    bool available = false;
    status = schedule->peek(cursor, static_cast<double>(index),
                            static_cast<double>(index), occurrenceValue,
                            available);
    require(static_cast<bool>(status) && available,
            "bounded-window schedule omitted an occurrence");
    WVOutputRouteView route;
    route.scheduleOrdinal = occurrenceValue.ordinal;
    route.observers = &observer;
    route.observerCount = 1;
    route.semanticScheduleOrdinal = 90;
    route.proposedScheduleCursor =
        &occurrenceValue.proposedCursor.values;
    route.schedulePayloadSchema = &schedule->payloadSchema();
    route.schedulePayload = &occurrenceValue.payload;
    route.scheduleCursorIdentity = occurrenceValue.cursorIdentity;
    WVOutputEvent event;
    event.eventOrdinal = index;
    event.scheduledTime = occurrenceValue.scheduledTime;
    event.state.waveVortex = stateStorage.view(occurrenceValue.scheduledTime);
    event.routes = &route;
    event.routeCount = 1;
    status = service->prepare(event);
    require(static_cast<bool>(status),
            "bounded-window occurrence preparation failed");
    WVObservationOccurrenceIdentity identity;
    status = service->preparedOccurrenceIdentity(route, observer, identity);
    require(static_cast<bool>(status),
            "bounded-window occurrence identity lookup failed");
    WVObservationBatch batch;
    status = service->observationBatch(identity, record, batch);
    require(static_cast<bool>(status),
            "bounded-window occurrence batch construction failed");
    cursor = occurrenceValue.proposedCursor;
  }
  return {service->persistentBytes(),
          service->metrics().occurrenceWorkspaceRetainedBytes,
          service->metrics().occurrenceWorkspaceMaximumLiveBytes,
          service->metrics().batchRetainedStorageBytes,
          service->metrics().batchMaximumLiveBytes};
}

void testBoundedStorageAndMoveOwnership() {
  const auto shortWindow = runOccurrenceWindow(20);
  const auto longWindow = runOccurrenceWindow(200);
  require(shortWindow.persistent == longWindow.persistent &&
              shortWindow.retained == longWindow.retained &&
              shortWindow.maximumLive == longWindow.maximumLive &&
              shortWindow.batchRetained == longWindow.batchRetained &&
              shortWindow.batchMaximumLive ==
                  longWindow.batchMaximumLive &&
              shortWindow.retained <= shortWindow.maximumLive &&
              shortWindow.batchRetained > 0 &&
              shortWindow.batchMaximumLive > 0,
          "occurrence storage grew with complete future event count");

  WVObserverOutputPlan plan;
  WVObserverOccurrencePositionSetPlan movedPositionSet;
  movedPositionSet.identifier = "moved-samples";
  plan.occurrencePositionSets.push_back(std::move(movedPositionSet));
  plan.occurrenceValues.push_back({"moved-value", 0});
  WVObserverOccurrenceWorkspace workspace;
  workspace.prepareFor(plan);
  workspace.positionSets[0].extents = {2};
  workspace.positionSets[0].sampleTimes = {1.0, 1.25};
  workspace.positionSets[0].x = {2.0, 3.0};
  workspace.positionSets[0].y = {4.0, 5.0};
  workspace.positionSets[0].z = {-6.0, -7.0};
  double *stored = nullptr;
  require(static_cast<bool>(workspace.resizeReal(0, {2}, stored)) &&
              stored != nullptr,
          "occurrence move workspace allocation failed");
  stored[0] = 8.0;
  stored[1] = 9.0;
  const auto fingerprint = workspace.geometryFingerprint();
  const auto retained = workspace.retainedBytes();
  WVObserverOccurrenceWorkspace moved = std::move(workspace);
  require(moved.geometryFingerprint() == fingerprint &&
              moved.retainedBytes() == retained &&
              moved.positionSets[0].x == std::vector<double>({2.0, 3.0}) &&
              moved.values[0].real64 == std::vector<double>({8.0, 9.0}),
          "moved occurrence workspace lost evaluator-owned geometry or data");
}

void testPayloadSchemaPreflight() {
  const auto catalog = extensionCatalog();
  const auto observers = protocolDescriptor(catalog, 1.0, false);
  WVOutputPlan plan;
  auto status = WVOutputPlan::create(observers, catalog, 0.0, 1.0, {}, plan);
  require(static_cast<bool>(status),
          "payload-mismatch output planning failed unexpectedly");
  std::unique_ptr<WVObserverOutputEvaluationService> service;
  status = WVObserverOutputEvaluationService::create(
      configuration(), false, observers,
      std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status),
          "payload-mismatch source construction failed");
  status = service->preflight(plan);
  require(!status && service->metrics().occurrencePreparationCount == 0 &&
              service->metrics().fieldEvaluationCount == 0,
          "incompatible schedule/observer payload schema reached event "
          "preparation");
}

void testOccurrenceStateTimeAndCoordinateBinding() {
  const auto catalog = extensionCatalog();
  WVPortableObserverRecord source;
  auto observer = eventVariableOneDimensionalGeometryRecord(
      "mismatched-coordinate-proof", "mismatched_coordinate_proof");
  observer.typeIdentifier = mismatchedCoordinateObservationType;
  source.observers.push_back(std::move(observer));
  WVPortableObserverDescriptor observers;
  auto status =
      WVPortableObserverDescriptor::create(source, catalog, observers);
  require(static_cast<bool>(status),
          "mismatched-coordinate descriptor construction failed: " +
              status.message);

  std::unique_ptr<WVObserverOutputEvaluationService> service;
  status = WVObserverOutputEvaluationService::create(
      configuration(), false, observers,
      std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status),
          "mismatched-coordinate source construction failed: " +
              status.message);

  const auto resolved = occurrence(
      catalog, occurrenceSchedule(3.5, 3, 1, 2.0, 2.0, 1.0));
  const auto &record = observers.observers().front();
  const auto *implementation = observers.resolvedObserver(record);
  require(implementation != nullptr,
          "mismatched-coordinate observer is unresolved");
  WVOutputObserverView observerView{0, &record, implementation};
  WVOutputRouteView route;
  route.scheduleOrdinal = resolved.value.ordinal;
  route.observers = &observerView;
  route.observerCount = 1;
  route.semanticScheduleOrdinal = 9;
  route.proposedScheduleCursor = &resolved.value.proposedCursor.values;
  route.schedulePayloadSchema = &resolved.schedule->payloadSchema();
  route.schedulePayload = &resolved.value.payload;
  route.scheduleCursorIdentity = resolved.value.cursorIdentity;
  auto ownedState = state(configuration());
  WVOutputEvent event;
  event.eventOrdinal = 4;
  event.scheduledTime = resolved.value.scheduledTime;
  event.kind = WVOutputEventKind::acceptedEndpoint;
  event.state.waveVortex = ownedState.view(event.scheduledTime + 0.25);
  event.routes = &route;
  event.routeCount = 1;

  status = service->prepare(event);
  require(!status &&
              status.message.find("scheduled trigger time") !=
                  std::string::npos &&
              service->metrics().occurrencePreparationCount == 0 &&
              service->metrics().fieldEvaluationCount == 0,
          "occurrence preparation accepted a state from another model time");

  event.state.waveVortex.t = event.scheduledTime;
  status = service->prepare(event);
  require(!status &&
              status.message.find("differ from interpolation geometry") !=
                  std::string::npos &&
              service->metrics().occurrencePreparationCount == 0 &&
              service->metrics().fieldEvaluationCount == 0,
          "occurrence preparation accepted independently persisted "
          "coordinates that differ from interpolation geometry");
}

void testOccurrenceFieldAffineRejectedAtConstruction() {
  const auto catalog = extensionCatalog();
  for (const auto &[scale, offset] :
       std::array<std::pair<double, double>, 2>{{{2.0, 0.0},
                                                 {1.0, -1.0}}}) {
    WVPortableObserverRecord source;
    auto observer = eventVariableOneDimensionalGeometryRecord(
        "affine-occurrence-field-proof", "affine_occurrence_field_proof");
    observer.typeIdentifier = affineOccurrenceFieldObservationType;
    observer.configuration.schemaIdentifier =
        std::string("wv-test-") + affineOccurrenceFieldObservationType +
        "-configuration-v1";
    observer.outputScale = scale;
    observer.outputOffset = offset;
    source.observers.push_back(std::move(observer));
    WVPortableObserverDescriptor observers;
    auto status =
        WVPortableObserverDescriptor::create(source, catalog, observers);
    require(static_cast<bool>(status),
            "affine occurrence-field descriptor construction failed: " +
                status.message);

    std::unique_ptr<WVObserverOutputEvaluationService> service;
    status = WVObserverOutputEvaluationService::create(
        configuration(), false, observers,
        std::make_unique<WVReferenceFFTEngine>(), service);
    require(!status &&
                status.code == WVKernelStatusCode::invalidConfiguration &&
                status.message ==
                    "Occurrence-field affine transforms are unsupported." &&
                service == nullptr,
            "occurrence-field affine transform reached runtime evaluation");
  }
}

} // namespace

int main() {
  try {
    testProviders();
    testIntegratedExactDenseAndRetry(false);
    testIntegratedExactDenseAndRetry(true);
    testIntegratedStateCoupledGeometry(false);
    testIntegratedStateCoupledGeometry(true);
    testNestedRaggedNetCDFPersistenceAndAppend();
    testOccurrenceCompatibilityIdentity();
    testSegmentedOccurrenceIdentity();
    testBoundedStorageAndMoveOwnership();
    testPayloadSchemaPreflight();
    testOccurrenceStateTimeAndCoordinateBinding();
    testOccurrenceFieldAffineRejectedAtConstruction();
    testCoincidentPrimitiveReconstructionReuse();
    std::cout << "Observation-occurrence provider tests passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Test failure: " << error.what() << '\n';
    return 1;
  }
}
