#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"

#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

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
  WVKernelStatus validate(
      const WVObserverRecord &record,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    if (!record.stateBlockIdentifiers.empty() || record.fieldNames.size() != 1 ||
        record.x.empty() || record.x.size() != record.y.size() ||
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

const std::shared_ptr<const WVExtensionCatalog> &extensionCatalog() {
  static const auto catalog = [] {
    WVExtensionCatalogBuilder builder;
    auto status = addBuiltInExtensions(builder);
    if (status)
      status = builder.addObserverFactory(
          {"WVTestPortablePointDiagnostic", 1,
           [](const WVObserverRecord &, const WVPortableTypedRecord &,
              std::shared_ptr<const WVObservingSystem> &result) {
             result = std::make_shared<
                 WVTestPortablePointDiagnosticImplementation>();
             return WVKernelStatus::ok();
           }});
    std::shared_ptr<const WVExtensionCatalog> result;
    if (status)
      status = builder.freeze(result);
    if (!status)
      throw std::runtime_error(status.message);
    return result;
  }();
  return catalog;
}

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
  WVState view() const {
    return {3.0,
            0.0,
            {{Ap.data(), shape}, {Am.data(), shape}, {A0.data(), shape}}};
  }
};

OwnedState state(const WVTransformConstantStratificationConfiguration &config) {
  WVTransformConstantStratificationDescriptor descriptor;
  require(static_cast<bool>(
              WVTransformConstantStratificationDescriptor::create(config,
                                                                   descriptor)),
          "descriptor construction failed");
  OwnedState result;
  result.shape = descriptor.spectralShape();
  result.Ap.resize(result.shape.elementCount());
  result.Am.resize(result.shape.elementCount());
  result.A0.resize(result.shape.elementCount());
  for (std::size_t index = 0; index < result.Ap.size(); ++index) {
    const double x = static_cast<double>(index + 1);
    result.Ap[index] = {1e-3 * std::sin(x), 8e-4 * std::cos(x)};
    result.Am[index] = {-7e-4 * std::cos(0.3 * x), 9e-4 * std::sin(0.2 * x)};
    result.A0[index] = {6e-4 * std::sin(0.4 * x), 5e-4 * std::cos(0.7 * x)};
  }
  return result;
}

WVPortableObserverDescriptor descriptor() {
  WVPortableObserverRecord record;
  for (const char *name : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {name,
         WVStateScalarType::complex64,
         {6, 1},
         WVToleranceKind::coefficientEnergyScaled,
         1e-6,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  for (const char *name : {"particle-x", "particle-y"})
    record.stateBlocks.push_back(
        {name, WVStateScalarType::real64, {2},
         WVToleranceKind::uniformAbsolute, 1e-5,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "WVCoefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers.push_back(coefficients);
  WVObserverRecord fields;
  fields.identifier = "fields-a";
  fields.name = "WVEulerianFields";
  fields.typeIdentifier = "WVEulerianFields";
  fields.fieldNames = {"Ap", "u", "psi"};
  record.observers.push_back(fields);
  fields.identifier = "fields-b";
  fields.fieldNames = {"u"};
  record.observers.push_back(fields);
  WVObserverRecord mooring;
  mooring.identifier = "mooring";
  mooring.name = "central";
  mooring.typeIdentifier = "WVMooring";
  mooring.fieldNames = {"u", "eta"};
  mooring.x = {-1.0, 8000.0};
  mooring.y = {6000.0, 2999.0};
  record.observers.push_back(mooring);
  WVObserverRecord particles;
  particles.identifier = "particles";
  particles.name = "drifters";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.stateBlockIdentifiers = {"particle-x", "particle-y"};
  particles.fieldNames = {"u", "rho_e"};
  particles.x = {-10.0, 8100.0};
  particles.y = {20.0, -30.0};
  particles.z = {-200.0, -600.0};
  particles.isXYOnly = true;
  particles.horizontalAbsoluteTolerance = 1e-5;
  particles.trackedFieldInterpolation = WVPositionInterpolation::spline;
  record.observers.push_back(particles);
  WVObserverRecord diagnostic;
  diagnostic.identifier = "point-diagnostic";
  diagnostic.name = "diagnostic";
  diagnostic.typeIdentifier = "WVTestPortablePointDiagnostic";
  diagnostic.fieldNames = {"u"};
  diagnostic.x = {1000.0, 3500.0};
  diagnostic.y = {1200.0, 4200.0};
  diagnostic.z = {-250.0, -750.0};
  diagnostic.outputScale = 2.5;
  diagnostic.outputOffset = -1.25;
  record.observers.push_back(diagnostic);
  WVPortableObserverDescriptor result;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(record, extensionCatalog(), result)),
          "observer descriptor construction failed");
  return result;
}

const WVObservationVariable &findVariable(const WVObservationSchema &schema,
                                          const char *name) {
  const auto found = std::find_if(
      schema.variables.begin(), schema.variables.end(),
      [&](const auto &variable) { return variable.name == name; });
  require(found != schema.variables.end(),
          "missing schema variable " + std::string(name));
  return *found;
}

std::vector<std::size_t>
variableExtents(const WVObservationSchema &schema,
                const WVObservationVariable &variable) {
  std::vector<std::size_t> extents;
  extents.reserve(variable.dimensionIdentifiers.size());
  for (const auto &identifier : variable.dimensionIdentifiers) {
    const auto axis = std::find_if(
        schema.axes.begin(), schema.axes.end(), [&](const auto &candidate) {
          return candidate.identifier == identifier;
        });
    require(axis != schema.axes.end(), "missing schema axis " + identifier);
    extents.push_back(axis->extent);
  }
  return extents;
}

const WVObservationValue &findValue(const WVObservationSchema &schema,
                                    const WVObservationBatch &batch,
                                    const WVObservationVariable &variable) {
  const auto variableIndex =
      static_cast<std::size_t>(&variable - schema.variables.data());
  const auto found =
      std::find_if(batch.values.begin(), batch.values.end(),
                   [&](const auto &value) {
                     return value.resolvedVariableIndex == variableIndex;
                   });
  require(found != batch.values.end(),
          "missing batch value " + variable.identifier);
  return *found;
}

void testService(bool linear) {
  auto config = configuration();
  auto observers = descriptor();
  std::unique_ptr<WVObserverOutputEvaluationService> service;
  auto status = WVObserverOutputEvaluationService::create(
      config, linear, observers, std::make_unique<WVReferenceFFTEngine>(),
      service);
  require(static_cast<bool>(status), "service construction failed");
  std::unique_ptr<WVFieldEvaluationService> sharedFields;
  status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), sharedFields);
  require(static_cast<bool>(status), "shared field service creation failed");
  status = service->useFieldEvaluationService(*sharedFields);
  require(static_cast<bool>(status), "shared field service binding failed");
  const auto requireRejectedConfiguration = [&](auto mutation,
                                                const char *message) {
    auto incompatibleConfiguration = config;
    mutation(incompatibleConfiguration);
    std::unique_ptr<WVFieldEvaluationService> incompatibleFields;
    auto incompatibleStatus = WVFieldEvaluationService::create(
        incompatibleConfiguration, std::make_unique<WVReferenceFFTEngine>(),
        incompatibleFields);
    require(static_cast<bool>(incompatibleStatus),
            "incompatible field service construction failed");
    incompatibleStatus =
        service->useFieldEvaluationService(*incompatibleFields);
    require(incompatibleStatus.code ==
                WVKernelStatusCode::invalidConfiguration,
            message);
  };
  requireRejectedConfiguration(
      [](auto &value) { value.planetaryRadius += 1.0; },
      "borrowed field service accepted a different planetary radius");
  requireRejectedConfiguration(
      [](auto &value) { value.rotationRate += 1e-8; },
      "borrowed field service accepted a different rotation rate");
  requireRejectedConfiguration(
      [](auto &value) { value.latitude += 1.0; },
      "borrowed field service accepted a different latitude");
  require(service->metrics().uniqueFieldOutputCount == 5,
          "field requests were not deduplicated");
  require(service->metrics().sharedFieldReuseCount == 1,
          "shared full-grid field was not recorded");

  WVObservationSchema fieldSchema;
  status =
      service->observationSchema(observers.observers()[1], fieldSchema);
  require(static_cast<bool>(status), "field schema query failed");
  const auto &coefficientVariable = findVariable(fieldSchema, "Ap");
  const auto &psiVariable = findVariable(fieldSchema, "psi");
  const auto &fieldVariable = findVariable(fieldSchema, "u");
  require(coefficientVariable.scalarType ==
              WVObservationScalarType::complex64,
          "coefficient type mismatch");
  require(coefficientVariable.layout ==
              (linear ? WVObservationValueLayout::initialValue
                      : WVObservationValueLayout::record),
          "coefficient layout mismatch");
  require(psiVariable.layout ==
              (linear ? WVObservationValueLayout::initialValue
                      : WVObservationValueLayout::record),
          "psi layout mismatch");
  require(fieldVariable.layout == WVObservationValueLayout::record,
          "u must remain a record value");

  WVObservationSchema mooringSchema;
  status =
      service->observationSchema(observers.observers()[3], mooringSchema);
  require(static_cast<bool>(status), "mooring schema query failed");
  const auto &mooringVariable = findVariable(mooringSchema, "central_u");
  require(variableExtents(mooringSchema, mooringVariable) ==
              std::vector<std::size_t>({config.Nz, 2}),
          "mooring schema mismatch");

  WVObservationSchema particleSchema;
  status =
      service->observationSchema(observers.observers()[4], particleSchema);
  require(static_cast<bool>(status), "particle schema query failed");
  const auto &particleVariable = findVariable(particleSchema, "drifters_u");
  require(variableExtents(particleSchema, particleVariable) ==
              std::vector<std::size_t>({2}) &&
              particleVariable.attributes.size() == 3,
          "particle tracked-field schema mismatch");

  WVObservationSchema diagnosticSchema;
  status =
      service->observationSchema(observers.observers()[5], diagnosticSchema);
  require(static_cast<bool>(status),
          "point-diagnostic schema query failed");
  const auto &diagnosticVariable =
      findVariable(diagnosticSchema, "diagnostic_value");
  require(variableExtents(diagnosticSchema, diagnosticVariable) ==
              std::vector<std::size_t>({2}),
          "point-diagnostic schema mismatch");

  auto owned = state(config);
  if (linear) {
    status = service->prepareInitial(owned.view());
    require(static_cast<bool>(status), "initial observer preparation failed");
    WVObservationBatch initialBatch;
    status = service->initialObservationBatch(observers.observers()[1],
                                              initialBatch);
    require(static_cast<bool>(status), "initial observation batch failed");
    const auto &initialPsi =
        findValue(fieldSchema, initialBatch, psiVariable);
    require(initialPsi.real64Data() != nullptr,
            "initial-only field value failed");
    const auto &initialCoefficient =
        findValue(fieldSchema, initialBatch, coefficientVariable);
    require(initialCoefficient.complex64Data() == owned.Ap.data(),
            "initial coefficient must remain a borrowed view");
  }
  WVOutputEvent event;
  event.eventOrdinal = 2;
  event.scheduledTime = 3.0;
  event.state.waveVortex = owned.view();
  WVAdditionalStateBlockLayout xLayout;
  xLayout.identifier = "particle-x";
  xLayout.scalarType = WVStateScalarType::real64;
  xLayout.elementCount = 2;
  WVAdditionalStateBlockLayout yLayout = xLayout;
  yLayout.identifier = "particle-y";
  const std::array<double, 2> particleX{{-10.0, 8100.0}};
  const std::array<double, 2> particleY{{20.0, -30.0}};
  const std::array<WVAdditionalStateBlockConstView, 2> particleBlocks{{
      {&xLayout, particleX.data(), nullptr},
      {&yLayout, particleY.data(), nullptr}}};
  event.state.additionalBlocks = particleBlocks.data();
  event.state.additionalBlockCount = particleBlocks.size();
  WVOutputSchedulePayload occurrencePayload;
  status = occurrencePayload.reset(emptyOutputSchedulePayloadSchema());
  require(static_cast<bool>(status), "empty occurrence payload setup failed");
  WVPortableTypedRecord occurrenceCursor;
  WVOutputGroupRecord activeScheduleRecord;
  activeScheduleRecord.identifier = "active-observer-output";
  activeScheduleRecord.observerIdentifiers = {
      observers.observers()[1].identifier, observers.observers()[3].identifier,
      observers.observers()[4].identifier, observers.observers()[5].identifier};
  const std::array<WVOutputObserverView, 4> activeObservers{{
      {1, &observers.observers()[1],
       observers.resolvedObserver(observers.observers()[1])},
      {3, &observers.observers()[3],
       observers.resolvedObserver(observers.observers()[3])},
      {4, &observers.observers()[4],
       observers.resolvedObserver(observers.observers()[4])},
      {5, &observers.observers()[5],
       observers.resolvedObserver(observers.observers()[5])}}};
  WVOutputRouteView activeRoute;
  activeRoute.observers = activeObservers.data();
  activeRoute.observerCount = activeObservers.size();
  activeRoute.scheduleOrdinal = 2;
  activeRoute.semanticScheduleOrdinal = 0;
  activeRoute.proposedScheduleCursor = &occurrenceCursor;
  activeRoute.schedulePayloadSchema = &emptyOutputSchedulePayloadSchema();
  activeRoute.schedulePayload = &occurrencePayload;
  activeRoute.scheduleCursorIdentity = 3;
  activeRoute.semanticScheduleRecord = &activeScheduleRecord;
  event.routes = &activeRoute;
  event.routeCount = 1;
  status = service->prepare(event);
  require(static_cast<bool>(status), "observer preparation failed");
  require(service->metrics().fieldEvaluationCount == (linear ? 3U : 2U),
          "initial and time-series plans were not evaluated independently");

  const auto eventBatch = [&](std::size_t routedObserverIndex,
                              WVObservationBatch &batch) {
    WVObservationOccurrenceIdentity identity;
    auto batchStatus = service->preparedOccurrenceIdentity(
        activeRoute, activeObservers[routedObserverIndex], identity);
    require(static_cast<bool>(batchStatus),
            "prepared occurrence identity query failed");
    batchStatus = service->observationBatch(
        identity, *activeObservers[routedObserverIndex].record, batch);
    require(static_cast<bool>(batchStatus), "event observation batch failed");
  };

  WVObservationBatch fieldsBatch;
  eventBatch(0, fieldsBatch);
  if (!linear) {
    const auto &coefficientValue =
        findValue(fieldSchema, fieldsBatch, coefficientVariable);
    require(coefficientValue.complex64Data() == owned.Ap.data(),
            "coefficient output must remain a borrowed view");
  }
  const auto &fieldValue = findValue(fieldSchema, fieldsBatch, fieldVariable);
  require(fieldValue.real64Data() != nullptr &&
              fieldValue.elementCount() == config.Nx * config.Ny * config.Nz,
          "Eulerian field shape mismatch");

  WVObservationBatch mooringBatch;
  eventBatch(1, mooringBatch);
  const auto &mooringValue =
      findValue(mooringSchema, mooringBatch, mooringVariable);
  require(mooringValue.elementCount() == config.Nz * 2,
          "mooring value mismatch");

  WVObservationBatch particleBatch;
  eventBatch(2, particleBatch);
  const auto &particleValue =
      findValue(particleSchema, particleBatch, particleVariable);
  const auto *particleData = particleValue.real64Data();
  require(particleValue.elementCount() == 2 && particleData != nullptr &&
              std::isfinite(particleData[0]) &&
              std::isfinite(particleData[1]),
          "particle tracked-field value mismatch");

  WVObservationBatch diagnosticBatch;
  eventBatch(3, diagnosticBatch);
  const auto &diagnosticValue =
      findValue(diagnosticSchema, diagnosticBatch, diagnosticVariable);
  const auto *diagnosticData = diagnosticValue.real64Data();
  require(diagnosticValue.elementCount() == 2 && diagnosticData != nullptr,
          "point-diagnostic value failed");
  WVFieldSamplingRequest sampling;
  sampling.kind = WVFieldSamplingKind::positions;
  sampling.x = observers.observers()[5].x;
  sampling.y = observers.observers()[5].y;
  sampling.z = observers.observers()[5].z;
  WVFieldEvaluationPlan diagnosticPlan;
  status = sharedFields->createPlan({{"diagnostic-u", "u", sampling}},
                                    diagnosticPlan);
  require(static_cast<bool>(status),
          "point-diagnostic reference plan failed");
  std::array<double, 2> raw{};
  WVFieldOutputView rawView{raw.data(), raw.size()};
  status = sharedFields->evaluate(diagnosticPlan, owned.view(), &rawView, 1);
  require(static_cast<bool>(status),
          "point-diagnostic reference evaluation failed");
  for (std::size_t index = 0; index < raw.size(); ++index)
    require(std::abs(diagnosticData[index] - (2.5 * raw[index] - 1.25)) <
                1e-12,
            "point-diagnostic affine output changed");
  WVOutputObserverView routedObserver{
      1, &observers.observers()[1],
      observers.resolvedObserver(observers.observers()[1])};
  WVOutputRouteView passiveRoute;
  WVOutputSchedulePayload emptyPayload;
  WVPortableTypedRecord emptyCursor;
  status = emptyPayload.reset(emptyOutputSchedulePayloadSchema());
  require(static_cast<bool>(status), "empty occurrence payload setup failed");
  passiveRoute.observers = &routedObserver;
  passiveRoute.observerCount = 1;
  passiveRoute.scheduleOrdinal = 3;
  passiveRoute.semanticScheduleOrdinal = 0;
  passiveRoute.proposedScheduleCursor = &emptyCursor;
  passiveRoute.schedulePayloadSchema = &emptyOutputSchedulePayloadSchema();
  passiveRoute.schedulePayload = &emptyPayload;
  passiveRoute.scheduleCursorIdentity = 4;
  event.eventOrdinal = 3;
  event.routes = &passiveRoute;
  event.routeCount = 1;
  status = service->prepare(event);
  require(static_cast<bool>(status) &&
              service->metrics().skippedParticleEvaluationCount == 1,
          "route-aware preparation evaluated unrouted particles");
}

} // namespace

int main() {
  try {
    testService(false);
    testService(true);
    std::cout << "Passive observer output evaluation contracts passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
