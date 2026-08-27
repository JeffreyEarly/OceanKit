#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WVTestExtensionCatalog.hpp"

#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
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

class WVTestPortableTracerImplementation final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestPortableTracer";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.persistedName = observer.name;
    return WVKernelStatus::ok();
  }
  WVKernelStatus bindIntegration(
      const WVObserverRecord &observer,
      const WVObserverIntegrationBinder &binder) const override {
    return binder.advectedScalar
               ? binder.advectedScalar(observer)
               : WVKernelStatus{WVKernelStatusCode::unsupportedOperation,
                                "Test tracer binding is unavailable."};
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &blocks,
      std::map<std::string, std::size_t> &owners) const override {
    if (observer.stateBlockIdentifiers.size() != 1)
      return {WVKernelStatusCode::invalidConfiguration,
              "Test tracer requires one state block."};
    const auto found = blocks.find(observer.stateBlockIdentifiers.front());
    if (found == blocks.end() ||
        found->second->scalarType != WVStateScalarType::real64 ||
        found->second->ownership != WVStateOwnership::integratorOwned)
      return {WVKernelStatusCode::invalidConfiguration,
              "Test tracer state block is incompatible."};
    ++owners.at(observer.stateBlockIdentifiers.front());
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
          {"WVTestPortableTracer", 1,
           [](const WVObserverRecord &, const WVPortableTypedRecord &,
              std::shared_ptr<const WVObservingSystem> &result) {
             result = std::make_shared<WVTestPortableTracerImplementation>();
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

WVTransformConstantStratificationConfiguration configuration(bool hydrostatic) {
  WVTransformConstantStratificationConfiguration value;
  value.Nx = 6;
  value.Ny = 5;
  value.Nz = 7;
  value.Nj = 4;
  value.Lx = 15000.0;
  value.Ly = 12000.0;
  value.Lz = 1300.0;
  value.N0 = 5.2e-3;
  value.rho0 = 1027.0;
  value.g = 9.80665;
  value.planetaryRadius = 6.3712e6;
  value.rotationRate = 7.292115e-5;
  value.latitude = 33.0;
  value.isHydrostatic = hydrostatic;
  value.shouldAntialias = true;
  return value;
}

WVPortableObserverDescriptor descriptorFor(
    const WVTransformConstantStratificationConfiguration &configuration) {
  WVTransformConstantStratificationDescriptor transform;
  require(static_cast<bool>(WVTransformConstantStratificationDescriptor::create(
              configuration, transform)),
          "transform descriptor");
  WVPortableObserverRecord record;
  const std::vector<std::size_t> shape{configuration.Nj, transform.Nkl()};
  for (const auto *name : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {name, WVStateScalarType::complex64, shape,
         WVToleranceKind::coefficientEnergyScaled, 1e-6,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  const auto addReal = [&](const char *name, std::size_t count,
                           double tolerance) {
    record.stateBlocks.push_back(
        {name, WVStateScalarType::real64, {count},
         WVToleranceKind::uniformAbsolute, tolerance,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  };
  addReal("surface-x", 3, 1e-5);
  addReal("surface-y", 3, 1e-5);
  addReal("volume-x", 2, 1e-5);
  addReal("volume-y", 2, 1e-5);
  addReal("volume-z", 2, 1e-6);
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers.push_back(coefficients);
  WVObserverRecord surface;
  surface.identifier = "surface";
  surface.name = "surfaceParticles";
  surface.typeIdentifier = "WVLagrangianParticles";
  surface.stateBlockIdentifiers = {"surface-x", "surface-y"};
  surface.x = {-10.0, 1000.0, configuration.Lx + 3.0};
  surface.y = {20.0, -100.0, 2000.0};
  surface.z = {-100.0, -300.0, -500.0};
  surface.isXYOnly = true;
  surface.advectionInterpolation = WVPositionInterpolation::linear;
  surface.trackedFieldInterpolation = WVPositionInterpolation::spline;
  surface.fieldNames = {"u", "rho_e"};
  surface.horizontalAbsoluteTolerance = 1e-5;
  record.observers.push_back(surface);
  WVObserverRecord volume;
  volume.identifier = "volume";
  volume.name = "volumeParticles";
  volume.typeIdentifier = "WVLagrangianParticles";
  volume.stateBlockIdentifiers = {"volume-x", "volume-y", "volume-z"};
  volume.x = {200.0, 3000.0};
  volume.y = {400.0, 5000.0};
  volume.z = {-200.0, -900.0};
  volume.isXYOnly = false;
  volume.advectionInterpolation = WVPositionInterpolation::spline;
  volume.trackedFieldInterpolation = WVPositionInterpolation::linear;
  volume.fieldNames = {"eta"};
  volume.horizontalAbsoluteTolerance = 1e-5;
  volume.verticalAbsoluteTolerance = 1e-6;
  record.observers.push_back(volume);
  WVPortableObserverDescriptor descriptor;
  const auto status = WVPortableObserverDescriptor::create(record, extensionCatalog(), descriptor);
  require(static_cast<bool>(status), status.message);
  return descriptor;
}

WVPortableObserverDescriptor descriptorWithTracers(
    const WVTransformConstantStratificationConfiguration &configuration) {
  auto record = descriptorFor(configuration).record();
  const std::vector<std::size_t> shape{
      configuration.Nx, configuration.Ny, configuration.Nz};
  for (const auto *identifier : {"dye", "temperature"}) {
    record.stateBlocks.push_back(
        {identifier, WVStateScalarType::real64, shape,
         WVToleranceKind::uniformAbsolute, 1e-8,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
    WVObserverRecord tracer;
    tracer.identifier = identifier;
    tracer.name = identifier;
    tracer.typeIdentifier = std::string(identifier) == "temperature"
                                ? "WVTestPortableTracer"
                                : "WVTracer";
    tracer.stateBlockIdentifiers = {identifier};
    tracer.shouldAntialias = std::string(identifier) == "dye";
    record.observers.push_back(std::move(tracer));
  }
  WVPortableObserverDescriptor descriptor;
  const auto status = WVPortableObserverDescriptor::create(record, extensionCatalog(), descriptor);
  require(static_cast<bool>(status), status.message);
  return descriptor;
}

struct Fixture {
  WVShape2D shape;
  std::vector<WVComplex64> coefficients;
  std::vector<WVComplex64> flux;
  WVAdditionalStateStorage stateStorage;
  WVAdditionalStateStorage fluxStorage;
  WVMutableIntegrationState state;
  WVIntegrationFlux rhs;
  std::vector<WVAdditionalStateBlockConstView> constViews;

  explicit Fixture(const WVIntegrationStateLayout &layout)
      : shape(layout.coefficientShape()),
        coefficients(3 * shape.elementCount()),
        flux(3 * shape.elementCount()) {
    require(static_cast<bool>(stateStorage.initialize(layout)),
            "state storage");
    require(static_cast<bool>(fluxStorage.initialize(layout)), "flux storage");
    const auto M = shape.elementCount();
    for (std::size_t index = 0; index < M; ++index) {
      const double value = static_cast<double>(index + 1);
      coefficients[index] = {1e-4 * std::sin(0.17 * value),
                             1e-4 * std::cos(0.11 * value)};
      coefficients[M + index] = {1e-4 * std::cos(0.13 * value),
                                 -1e-4 * std::sin(0.19 * value)};
      coefficients[2 * M + index] = {1e-4 * std::sin(0.07 * value),
                                     1e-4 * std::cos(0.23 * value)};
    }
    state = {{0.25,
              -0.5,
              {{coefficients.data(), shape},
               {coefficients.data() + M, shape},
               {coefficients.data() + 2 * M, shape}}},
             stateStorage.mutableBlocks(), stateStorage.blockCount()};
    rhs = {{{flux.data(), shape}, {flux.data() + M, shape},
            {flux.data() + 2 * M, shape}},
           fluxStorage.mutableBlocks(), fluxStorage.blockCount()};
  }

  WVIntegrationState constView() {
    return integrationConstView(state, constViews);
  }
};

void initializeParticles(Fixture &fixture,
                         const WVPortableObserverDescriptor &descriptor) {
  for (const auto &observer : descriptor.observers()) {
    if (observer.typeIdentifier != "WVLagrangianParticles")
      continue;
    const std::array<const std::vector<double> *, 3> values{
        {&observer.x, &observer.y, &observer.z}};
    for (std::size_t coordinate = 0;
         coordinate < observer.stateBlockIdentifiers.size(); ++coordinate)
      for (std::size_t block = 0; block < fixture.state.additionalBlockCount;
           ++block)
        if (fixture.state.additionalBlocks[block].layout->identifier ==
            observer.stateBlockIdentifiers[coordinate])
          std::copy(values[coordinate]->begin(), values[coordinate]->end(),
                    fixture.state.additionalBlocks[block].realData);
  }
}

std::size_t blockIndex(const WVMutableIntegrationState &state,
                       const std::string &identifier) {
  for (std::size_t index = 0; index < state.additionalBlockCount; ++index)
    if (state.additionalBlocks[index].layout->identifier == identifier)
      return index;
  throw std::runtime_error("missing state block " + identifier);
}

void initializeTracer(Fixture &fixture, const std::string &identifier,
                      const WVTransformConstantStratificationConfiguration &config,
                      double scale) {
  auto &block = fixture.state.additionalBlocks[blockIndex(fixture.state, identifier)];
  const double pi = std::acos(-1.0);
  for (std::size_t z = 0; z < config.Nz; ++z)
    for (std::size_t y = 0; y < config.Ny; ++y)
      for (std::size_t x = 0; x < config.Nx; ++x) {
        const auto index = x + config.Nx * (y + config.Ny * z);
        block.realData[index] =
            scale * (std::sin(2.0 * pi * static_cast<double>(x) /
                              static_cast<double>(config.Nx)) +
                     0.25 * std::cos(2.0 * pi * static_cast<double>(y) /
                                     static_cast<double>(config.Ny)) +
                     0.1 * std::cos(pi * static_cast<double>(z) /
                                    static_cast<double>(config.Nz - 1)));
      }
}

void testTracers(bool hydrostatic) {
  const auto config = configuration(hydrostatic);
  const auto descriptor = descriptorWithTracers(config);
  std::unique_ptr<WVConstantStratificationIntegrationSystem> system;
  auto status = WVConstantStratificationIntegrationSystem::create(
      config, {}, descriptor, extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(),
      system);
  require(static_cast<bool>(status), status.message);
  require(system->tracers().size() == 2, "tracers were not resolved");
  Fixture fixture(system->stateLayout());
  status = system->initializeParticleState(fixture.state);
  require(static_cast<bool>(status), status.message);
  initializeTracer(fixture, "dye", config, 1.0);
  initializeTracer(fixture, "temperature", config, -0.5);
  const auto scratchCapacityBytes =
      system->kernelMetrics().scratchCapacityBytes;
  const auto planCount = system->kernelMetrics().planCount;
  require(planCount == 18,
          "antialiased tracer did not prepare the scalar inverse plan");
  status = system->evaluateRightHandSide(fixture.constView(), fixture.rhs);
  require(static_cast<bool>(status), status.message);

  for (const auto *identifier : {"dye", "temperature"}) {
    const auto index = blockIndex(fixture.state, identifier);
    const auto &flux = fixture.rhs.additionalBlocks[index];
    bool hasNonzeroValue = false;
    for (std::size_t value = 0; value < flux.layout->elementCount; ++value) {
      require(std::isfinite(flux.realData[value]),
              std::string(identifier) + " tracer flux was not finite");
      hasNonzeroValue = hasNonzeroValue || std::abs(flux.realData[value]) > 0.0;
    }
    require(hasNonzeroValue,
            std::string(identifier) + " tracer flux was not evaluated");
  }
  require(system->metrics().tracerEvaluationCount == 2 &&
              system->metrics().antialiasedTracerEvaluationCount == 1,
          "mixed tracer antialias dispatch changed");
  require(system->kernelMetrics().advectionVelocityReconstructionCount == 1 &&
              system->kernelMetrics().scalarAdvectionCount == 2 &&
              system->kernelMetrics().scalarAntialiasCount == 1,
          "tracers did not share exactly one RHS velocity reconstruction");
  require(system->kernelMetrics().scratchCapacityBytes == scratchCapacityBytes &&
              system->kernelMetrics().planCount == planCount,
          "tracer evaluation added array-sized scratch or lazily created a plan");
  require(system->kernelMetrics().scalarForwardSeconds >= 0.0 &&
              system->kernelMetrics().scalarDerivativeAssemblySeconds >= 0.0 &&
              system->kernelMetrics().scalarVerticalDerivativeSeconds >= 0.0 &&
              system->kernelMetrics().scalarInverseSeconds >= 0.0 &&
              system->kernelMetrics().scalarProductSeconds >= 0.0 &&
              system->kernelMetrics().scalarAntialiasSeconds >= 0.0,
          "tracer stage timing metrics were invalid");
  require(system->fieldEvaluationService()->metrics().movingPrimitiveTransformCount == 0,
          "tracer RHS invoked the particle interpolation transform");

  const auto persistentBytes = system->persistentBytes();
  WVFixedStepRK4 rk4(*system, {true});
  status = rk4.prepareStateAfterRestart(fixture.state);
  require(static_cast<bool>(status), status.message);
  status = rk4.step(fixture.state, 1e-4);
  require(static_cast<bool>(status), status.message);
  require(rk4.lastAcceptedStep() != nullptr &&
              rk4.lastAcceptedStep()->denseOutput != nullptr,
          "tracer RK4 dense output missing");
  require(system->persistentBytes() == persistentBytes,
          "tracer integration grew persistent storage");
}

void testScalarAdvectionOperator(bool shouldAntialias) {
  const auto config = configuration(true);
  std::unique_ptr<WVTransformConstantStratificationKernel> kernel;
  auto status = WVTransformConstantStratificationKernel::create(
      config, std::make_unique<WVReferenceFFTEngine>(), kernel);
  require(static_cast<bool>(status), status.message);
  const WVShape3D shape{config.Nx, config.Ny, config.Nz};
  const auto R = shape.elementCount();
  std::vector<double> scalar(R);
  std::vector<double> fields(3 * R, 0.0);
  std::vector<double> flux(R, 0.0);
  const double pi = std::acos(-1.0);
  const double u = 0.75;
  const double v = -0.4;
  const double w = 0.2;
  const auto verticalMode = config.Nz - 2;
  for (std::size_t z = 0; z < config.Nz; ++z)
    for (std::size_t y = 0; y < config.Ny; ++y)
      for (std::size_t x = 0; x < config.Nx; ++x) {
        const auto index = x + config.Nx * (y + config.Ny * z);
        scalar[index] =
            std::sin(2.0 * pi * static_cast<double>(x) /
                     static_cast<double>(config.Nx)) +
            0.25 * std::cos(2.0 * pi * static_cast<double>(y) /
                            static_cast<double>(config.Ny)) +
            0.1 * std::cos(pi * static_cast<double>(verticalMode * z) /
                           static_cast<double>(config.Nz - 1));
        fields[index] = u;
        fields[R + index] = v;
        fields[2 * R + index] = w;
      }
  const WVRealVolumeConstView scalarView{scalar.data(), shape};
  const WVRealFieldBundleConstView fieldView{
      fields.data(), {config.Nx, config.Ny, config.Nz, 3}};
  WVRealVolumeView fluxView{flux.data(), shape};
  status = kernel->advectFGridScalar(scalarView, fieldView, shouldAntialias,
                                     fluxView);
  require(static_cast<bool>(status), status.message);
  const double k = 2.0 * pi / config.Lx;
  const double l = 2.0 * pi / config.Ly;
  double maximumError = 0.0;
  for (std::size_t z = 0; z < config.Nz; ++z)
    for (std::size_t y = 0; y < config.Ny; ++y)
      for (std::size_t x = 0; x < config.Nx; ++x) {
        const auto index = x + config.Nx * (y + config.Ny * z);
        const double expected =
            -u * k * std::cos(2.0 * pi * static_cast<double>(x) /
                              static_cast<double>(config.Nx)) +
            v * 0.25 * l *
                std::sin(2.0 * pi * static_cast<double>(y) /
                         static_cast<double>(config.Ny)) +
            w * 0.1 * pi * static_cast<double>(verticalMode) / config.Lz *
                std::sin(pi * static_cast<double>(verticalMode * z) /
                         static_cast<double>(config.Nz - 1));
        maximumError = std::max(maximumError,
                                std::abs(flux[index] - expected));
      }
  require(maximumError <= 1e-12,
          "shared scalar differential operator changed its sign or scaling");
}

void testIntegratedObservers(bool hydrostatic) {
  const auto config = configuration(hydrostatic);
  const auto descriptor = descriptorFor(config);
  WVFrozenForcingSchedule emptySchedule;
  std::unique_ptr<WVConstantStratificationIntegrationSystem> system;
  auto status = WVConstantStratificationIntegrationSystem::create(
      config, emptySchedule, descriptor, extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(), system);
  require(static_cast<bool>(status), status.message);
  require(system->particles().size() == 2, "particle systems not resolved");
  Fixture fixture(system->stateLayout());
  status = system->initializeParticleState(fixture.state);
  require(static_cast<bool>(status), status.message);
  const double unwrapped = fixture.state.additionalBlocks[0].realData[0];
  status = system->evaluateRightHandSide(fixture.constView(), fixture.rhs);
  require(static_cast<bool>(status), status.message);
  require(fixture.state.additionalBlocks[0].realData[0] == unwrapped,
          "particle evaluation wrapped accepted state");
  for (std::size_t block = 0; block < fixture.rhs.additionalBlockCount; ++block)
    for (std::size_t index = 0;
         index < fixture.rhs.additionalBlocks[block].layout->elementCount;
         ++index)
      require(std::isfinite(fixture.rhs.additionalBlocks[block].realData[index]),
              "particle RHS was not completely written");
  require(system->metrics().velocityFieldEvaluationCount == 1,
          "particle systems did not share one primitive transform");
  require(system->metrics().sharedRightHandSideContextCount == 1,
          "particle systems did not use one shared RHS context");
  require(system->fieldEvaluationService()->metrics().movingPrimitiveTransformCount == 0 &&
              system->fieldEvaluationService()->metrics().primitiveFieldReuseCount == 1,
          "particle interpolation recomputed prepared advection fields");
  require(system->fieldEvaluationService()->metrics().movingPositionCount == 5,
          "moving position count changed");
  const auto persistentBytes = system->persistentBytes();
  const auto interpolationBytes = system->fieldEvaluationService()
                                      ->metrics()
                                      .movingInterpolationWorkspaceBytes;
  status = system->evaluateRightHandSide(fixture.constView(), fixture.rhs);
  require(static_cast<bool>(status), status.message);
  require(system->persistentBytes() == persistentBytes &&
              system->fieldEvaluationService()
                      ->metrics()
                      .movingInterpolationWorkspaceBytes == interpolationBytes,
          "repeated particle RHS changed bounded persistent storage");

  WVFixedStepRK4 rk4(*system, {true});
  status = rk4.prepareStateAfterRestart(fixture.state);
  require(static_cast<bool>(status), status.message);
  status = rk4.step(fixture.state, 1e-3);
  require(static_cast<bool>(status), status.message);
  require(rk4.lastAcceptedStep() != nullptr &&
              rk4.lastAcceptedStep()->denseOutput != nullptr,
          "particle RK4 dense output missing");

  std::vector<double> acceptedParticleState;
  for (std::size_t block = 0; block < fixture.state.additionalBlockCount;
       ++block)
    acceptedParticleState.insert(
        acceptedParticleState.end(),
        fixture.state.additionalBlocks[block].realData,
        fixture.state.additionalBlocks[block].realData +
            fixture.state.additionalBlocks[block].layout->elementCount);

  Fixture dense(system->stateLayout());
  initializeParticles(dense, descriptor);
  status = rk4.evaluateDenseOutput(fixture.state.waveVortex.t - 5e-4,
                                   dense.state);
  require(static_cast<bool>(status), status.message);
  require(dense.state.waveVortex.t == fixture.state.waveVortex.t - 5e-4,
          "particle dense-output time changed");
  std::vector<double> afterDenseOutput;
  for (std::size_t block = 0; block < fixture.state.additionalBlockCount;
       ++block)
    afterDenseOutput.insert(
        afterDenseOutput.end(), fixture.state.additionalBlocks[block].realData,
        fixture.state.additionalBlocks[block].realData +
            fixture.state.additionalBlocks[block].layout->elementCount);
  require(afterDenseOutput == acceptedParticleState,
          "particle dense output mutated accepted state");

  WVAdaptiveRK23Options strictOptions;
  strictOptions.relativeTolerance = 1e-12;
  strictOptions.absoluteToleranceScale = 1e-12;
  WVAdaptiveRK23 adaptive(*system, strictOptions);
  Fixture adaptiveFixture(system->stateLayout());
  status = system->initializeParticleState(adaptiveFixture.state);
  require(static_cast<bool>(status), status.message);
  status = adaptive.prepareStateAfterRestart(adaptiveFixture.state);
  require(static_cast<bool>(status), status.message);
  status = adaptive.step(adaptiveFixture.state, 100.0);
  require(static_cast<bool>(status), status.message);
  require(adaptive.lastAcceptedStep() != nullptr,
          "particle adaptive step was not accepted");
  require(adaptive.metrics().rejectedStepCount > 0,
          "strict particle adaptive step did not exercise rejection");
  require(adaptive.lastAcceptedStep()->denseOutput != nullptr,
          "accepted particle adaptive step omitted dense output");
}

void testValidation() {
  auto config = configuration(true);
  auto descriptor = descriptorFor(config);
  std::unique_ptr<WVConstantStratificationIntegrationSystem> system;
  auto status = WVConstantStratificationIntegrationSystem::create(
      config, {}, descriptor, extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(),
      system);
  require(static_cast<bool>(status) && system,
          "valid integration system construction failed");
  std::unique_ptr<WVIntegrationErrorPolicy> errorPolicy;
  status = system->createErrorPolicy(1e-12, errorPolicy);
  require(static_cast<bool>(status) && errorPolicy,
          "particle error policy construction failed");
  require(std::abs(errorPolicy->absoluteTolerance(3, 0) - 1e-5) < 1e-20,
          "particle absolute tolerance was rescaled by the coefficient tolerance");
  auto record = descriptor.record();
  auto particle = std::find_if(record.observers.begin(), record.observers.end(),
                               [](const auto &observer) {
                                 return observer.identifier == "surface";
                               });
  particle->z.clear();
  WVPortableObserverDescriptor genericDescriptor;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(
              record, extensionCatalog(), genericDescriptor)),
          "generic descriptor should retain future 2-D XY allowance");
  status = WVConstantStratificationIntegrationSystem::create(
      config, {}, genericDescriptor, extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(),
      system);
  require(status.code == WVKernelStatusCode::invalidConfiguration && !system,
          "constant-stratification XY particles accepted missing fixed z");

  auto tracerRecord = descriptor.record();
  tracerRecord.stateBlocks.push_back(
      {"tracer", WVStateScalarType::real64, {config.Nx, config.Ny},
       WVToleranceKind::uniformAbsolute, 1e-6,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState});
  WVObserverRecord tracer;
  tracer.identifier = "tracer";
  tracer.name = "tracer";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"tracer"};
  tracer.isXYOnly = true;
  tracerRecord.observers.push_back(tracer);
  WVPortableObserverDescriptor tracerDescriptor;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(
              tracerRecord, extensionCatalog(), tracerDescriptor)),
          "generic tracer descriptor");
  const auto tracerStatus = WVConstantStratificationIntegrationSystem::create(
      config, {}, tracerDescriptor, extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(),
      system);
  require(tracerStatus.code == WVKernelStatusCode::unsupportedOperation,
          "constant-stratification system did not defer two-dimensional tracers");
}

} // namespace

int main() {
  try {
    testIntegratedObservers(true);
    testIntegratedObservers(false);
    testTracers(true);
    testTracers(false);
    testScalarAdvectionOperator(false);
    testScalarAdvectionOperator(true);
    testValidation();
    std::cout << "PASS TestWVLagrangianParticles\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
