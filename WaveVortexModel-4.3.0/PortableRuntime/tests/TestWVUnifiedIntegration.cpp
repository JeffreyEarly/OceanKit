#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WVTestExtensionCatalog.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

WVPortableObserverRecord record() {
  WVPortableObserverRecord result;
  const std::vector<std::size_t> coefficientShape{2, 3};
  for (const auto *identifier : {"Ap", "Am", "A0"})
    result.stateBlocks.push_back({identifier, WVStateScalarType::complex64,
                                  coefficientShape,
                                  WVToleranceKind::coefficientEnergyScaled, 0.0,
                                  WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
  result.stateBlocks.push_back({"particleX",
                                WVStateScalarType::real64,
                                {3},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
  result.stateBlocks.push_back({"particleY",
                                WVStateScalarType::real64,
                                {3},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
  result.stateBlocks.push_back({"tracerAmplitude",
                                WVStateScalarType::real64,
                                {1, 2, 2},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
  result.stateBlocks.push_back({"complexAuxiliary",
                                WVStateScalarType::complex64,
                                {2, 2},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::observerDerived,
                                WVRestartRequirement::derivedState});
  result.stateBlocks.push_back({"sampledVelocity",
                                WVStateScalarType::real64,
                                {3, 2},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::observerDerived,
                                WVRestartRequirement::derivedState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  result.observers.push_back(coefficients);
  WVObserverRecord particles;
  particles.identifier = "particles";
  particles.name = "Particles";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.stateBlockIdentifiers = {"particleX", "particleY"};
  particles.x = {0, 1, 2};
  particles.y = {3, 4, 5};
  particles.isXYOnly = true;
  particles.horizontalAbsoluteTolerance = 1e-8;
  result.observers.push_back(particles);
  WVObserverRecord tracer;
  tracer.identifier = "tracer";
  tracer.name = "Tracer";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"tracerAmplitude"};
  result.observers.push_back(tracer);
  result.outputFiles.push_back({"history",
                                "history.nc",
                                {{"state",
                                  "State",
                                  {1.0, 0.0, 10.0},
                                  {"coefficients", "particles", "tracer"},
                                  true}}});
  return result;
}

class LinearIntegrationSystem final : public WVIntegrationSystem {
public:
  class ErrorPolicy final : public WVIntegrationErrorPolicy {
  public:
    explicit ErrorPolicy(const WVIntegrationStateLayout &layout)
        : layout_(layout) {}
    std::size_t componentCount() const noexcept override {
      return 3 + layout_.additionalBlocks().size();
    }
    std::size_t elementCount(std::size_t component) const noexcept override {
      return component < 3 ? layout_.coefficientShape().elementCount()
                           : layout_.additionalBlocks()[component - 3].elementCount;
    }
    double absoluteTolerance(std::size_t, std::size_t) const noexcept override {
      return 1e-10;
    }
    std::size_t persistentBytes() const noexcept override { return 0; }
  private:
    const WVIntegrationStateLayout &layout_;
  };
  explicit LinearIntegrationSystem(WVIntegrationStateLayout layout,
                                   bool zeroDerivative = false,
                                   bool constraintsFSALCompatible = true,
                                   bool nonfiniteDerivative = false)
      : layout_(std::move(layout)), zeroDerivative_(zeroDerivative),
        constraintsFSALCompatible_(constraintsFSALCompatible),
        nonfiniteDerivative_(nonfiniteDerivative) {}
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(const WVIntegrationState &state,
                                       WVIntegrationFlux &rhs) override {
    if (rhs.additionalBlockCount != state.additionalBlockCount)
      return {WVKernelStatusCode::invalidShape, "RHS layout mismatch."};
    const WVComplexConstView source[] = {state.waveVortex.coefficients.Ap,
                                         state.waveVortex.coefficients.Am,
                                         state.waveVortex.coefficients.A0};
    WVComplexView destination[] = {rhs.waveVortex.Fp, rhs.waveVortex.Fm,
                                   rhs.waveVortex.F0};
    if (nonfiniteDerivative_) {
      const auto value = std::numeric_limits<double>::quiet_NaN();
      for (auto &component : destination)
        std::fill_n(component.data, component.shape.elementCount(),
                    WVComplex64{value, value});
      for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
        const auto &metadata = *state.additionalBlocks[block].layout;
        if (metadata.scalarType == WVStateScalarType::real64)
          std::fill_n(rhs.additionalBlocks[block].realData,
                      metadata.elementCount, value);
        else
          std::fill_n(rhs.additionalBlocks[block].complexData,
                      metadata.elementCount, WVComplex64{value, value});
      }
      ++evaluations;
      return WVKernelStatus::ok();
    }
    if (zeroDerivative_) {
      for (auto &component : destination)
        std::fill_n(component.data, component.shape.elementCount(),
                    WVComplex64{});
      for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
        const auto &metadata = *state.additionalBlocks[block].layout;
        if (metadata.scalarType == WVStateScalarType::real64)
          std::fill_n(rhs.additionalBlocks[block].realData,
                      metadata.elementCount, 0.0);
        else
          std::fill_n(rhs.additionalBlocks[block].complexData,
                      metadata.elementCount, WVComplex64{});
      }
      ++evaluations;
      return WVKernelStatus::ok();
    }
    for (std::size_t component = 0; component < 3; ++component)
      for (std::size_t i = 0; i < source[component].shape.elementCount(); ++i)
        destination[component].data[i] = {-source[component].data[i].real,
                                          -source[component].data[i].imag};
    for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
      const auto &metadata = *state.additionalBlocks[block].layout;
      if (metadata.scalarType == WVStateScalarType::real64)
        for (std::size_t i = 0; i < metadata.elementCount; ++i)
          rhs.additionalBlocks[block].realData[i] =
              -2.0 * state.additionalBlocks[block].realData[i];
      else
        for (std::size_t i = 0; i < metadata.elementCount; ++i) {
          const auto value = state.additionalBlocks[block].complexData[i];
          rhs.additionalBlocks[block].complexData[i] = {-value.imag,
                                                        value.real};
        }
    }
    ++evaluations;
    return WVKernelStatus::ok();
  }
  WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &state) override {
    std::size_t modified = 0;
    for (std::size_t block = 0; block < state.additionalBlockCount; ++block)
      if (state.additionalBlocks[block].layout->identifier == "particleX" ||
          state.additionalBlocks[block].layout->identifier == "particleY")
        for (std::size_t i = 0;
             i < state.additionalBlocks[block].layout->elementCount; ++i)
          if (state.additionalBlocks[block].realData[i] < 0) {
            state.additionalBlocks[block].realData[i] = 0;
            ++modified;
          }
    return {WVKernelStatus::ok(), modified,
            modified == 0 && constraintsFSALCompatible_};
  }
  WVKernelStatus createErrorPolicy(
      double, std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override {
    policy = std::make_unique<ErrorPolicy>(layout_);
    return WVKernelStatus::ok();
  }
  std::size_t evaluations = 0;

private:
  WVIntegrationStateLayout layout_;
  bool zeroDerivative_ = false;
  bool constraintsFSALCompatible_ = true;
  bool nonfiniteDerivative_ = false;
};

class A0OnlyIntegrationSystem final : public WVIntegrationSystem {
public:
  class ErrorPolicy final : public WVIntegrationErrorPolicy {
  public:
    explicit ErrorPolicy(std::size_t count) : count_(count) {}
    std::size_t componentCount() const noexcept override { return 1; }
    std::size_t elementCount(std::size_t component) const noexcept override {
      return component == 0 ? count_ : 0;
    }
    double absoluteTolerance(std::size_t,
                             std::size_t) const noexcept override {
      return 1e-10;
    }
    std::size_t persistentBytes() const noexcept override {
      return sizeof(*this);
    }

  private:
    std::size_t count_ = 0;
  };

  A0OnlyIntegrationSystem() {
    WVTransformStateDescription description{
        "WVTransformBarotropicQG", {8, 6},
        {{"A0", {24}, WVToleranceKind::coefficientEnergyScaled}}};
    const auto status = WVIntegrationStateLayout::createCoefficientOnly(
        std::move(description), layout_);
    require(static_cast<bool>(status), "A0-only transform preflight");
  }

  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(const WVIntegrationState &state,
                                       WVIntegrationFlux &rhs) override {
    auto status = validateIntegrationState(layout_, state);
    if (!status)
      return status;
    if (rhs.coefficientFamilyCount != 1 ||
        rhs.coefficientFamilies == nullptr ||
        rhs.coefficientFamilies[0].layout !=
            &layout_.coefficientFamilies()[0])
      return {WVKernelStatusCode::invalidShape,
              "A0-only RHS family layout mismatch."};
    const auto source = coefficientFamilyView(layout_, state, 0);
    auto destination = coefficientFamilyView(layout_, rhs, 0);
    for (std::size_t index = 0;
         index < layout_.coefficientFamilies()[0].elementCount; ++index)
      destination.data[index] = {-source.data[index].real,
                                 -source.data[index].imag};
    ++rightHandSideCount;
    return WVKernelStatus::ok();
  }
  WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &) override {
    return {WVKernelStatus::ok(), 0, true};
  }
  WVKernelStatus createErrorPolicy(
      double, std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override {
    policy = std::make_unique<ErrorPolicy>(layout_.coefficientElementCount());
    return WVKernelStatus::ok();
  }

  std::size_t rightHandSideCount = 0;

private:
  WVIntegrationStateLayout layout_;
};

struct StateFixture {
  WVShape2D shape{2, 3};
  std::vector<WVComplex64> coefficients =
      std::vector<WVComplex64>(18, {1.0, 0.5});
  WVAdditionalStateStorage extra;
  WVMutableIntegrationState state;
  std::vector<WVComplex64> outputCoefficients = std::vector<WVComplex64>(18);
  WVAdditionalStateStorage outputExtra;
  WVMutableIntegrationState output;
  explicit StateFixture(const WVIntegrationStateLayout &layout) {
    require(static_cast<bool>(extra.initialize(layout)),
            "initialize additional state");
    require(static_cast<bool>(outputExtra.initialize(layout)),
            "initialize output state");
    state = {{0.0,
              0.0,
              {{coefficients.data(), shape},
               {coefficients.data() + 6, shape},
               {coefficients.data() + 12, shape}}},
             extra.mutableBlocks(),
             extra.blockCount()};
    output = {{0.0,
               0.0,
               {{outputCoefficients.data(), shape},
                {outputCoefficients.data() + 6, shape},
                {outputCoefficients.data() + 12, shape}}},
              outputExtra.mutableBlocks(),
              outputExtra.blockCount()};
    for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
      const auto &metadata = *state.additionalBlocks[block].layout;
      if (metadata.scalarType == WVStateScalarType::real64)
        std::fill_n(state.additionalBlocks[block].realData,
                    metadata.elementCount, 1.0);
      else
        std::fill_n(state.additionalBlocks[block].complexData,
                    metadata.elementCount, WVComplex64{1.0, 0.0});
    }
  }
};

void testContracts(WVPortableObserverDescriptor &descriptor,
                   WVIntegrationStateLayout &layout) {
  auto source = record();
  require(static_cast<bool>(
              WVPortableObserverDescriptor::create(source, test::extensionCatalog(), descriptor)),
          "valid observer descriptor");
  const auto roundTrip = descriptor.record();
  require(roundTrip.schemaIdentifier == source.schemaIdentifier &&
              roundTrip.stateBlocks.size() == source.stateBlocks.size() &&
              roundTrip.outputFiles[0].groups[0].observerIdentifiers ==
                  source.outputFiles[0].groups[0].observerIdentifiers,
          "deterministic descriptor record");
  require(test::extensionCatalog()->observers().registration(
              "WVLagrangianParticles", WVPortablePairContractVersion) != nullptr,
          "factory identity");
  require(test::extensionCatalog()->observers().registration(
              "WVUnknownObserver", WVPortablePairContractVersion) == nullptr,
          "unknown identity rejected");
  auto duplicate = source;
  duplicate.stateBlocks.push_back(duplicate.stateBlocks.front());
  WVPortableObserverDescriptor ignored;
  require(!WVPortableObserverDescriptor::create(duplicate, test::extensionCatalog(), ignored),
          "duplicate block rejected");
  auto badReference = source;
  badReference.observers.back().stateBlockIdentifiers = {"missing"};
  require(!WVPortableObserverDescriptor::create(badReference, test::extensionCatalog(), ignored),
          "unknown state reference rejected");
  auto orphan = source;
  orphan.observers.erase(orphan.observers.begin() + 2);
  require(!WVPortableObserverDescriptor::create(orphan, test::extensionCatalog(), ignored),
          "orphan integrator-owned block rejected");
  auto sharedTracer = source;
  auto secondTracer = sharedTracer.observers[2];
  secondTracer.identifier = "secondTracer";
  secondTracer.name = "Second tracer";
  sharedTracer.observers.push_back(secondTracer);
  require(!WVPortableObserverDescriptor::create(sharedTracer, test::extensionCatalog(), ignored),
          "state block shared by two tracers rejected");
  auto sharedParticles = source;
  auto secondParticles = sharedParticles.observers[1];
  secondParticles.identifier = "secondParticles";
  secondParticles.name = "Second particles";
  sharedParticles.observers.push_back(secondParticles);
  require(!WVPortableObserverDescriptor::create(sharedParticles, test::extensionCatalog(), ignored),
          "state blocks shared by two particle systems rejected");
  auto mixedOwners = source;
  mixedOwners.observers[2].stateBlockIdentifiers = {"particleX"};
  require(!WVPortableObserverDescriptor::create(mixedOwners, test::extensionCatalog(), ignored),
          "state block shared by particle and tracer observers rejected");
  require(static_cast<bool>(
              WVIntegrationStateLayout::create({2, 3}, descriptor, layout)),
          "integration layout");
  require(layout.additionalBlocks().size() == 3 &&
              layout.additionalBlocks()[0].identifier == "particleX" &&
              layout.additionalBlocks()[1].identifier == "particleY" &&
              layout.additionalBlocks()[2].identifier == "tracerAmplitude",
          "derived block excluded and order frozen");
  require(layout.realElementCount() == 10 && layout.complexElementCount() == 0,
          "integration-state counts");
  std::size_t expectedLayoutStorage =
      layout.transformIdentifier().capacity() +
      layout.spatialDimensions().capacity() * sizeof(std::size_t) +
      layout.coefficientFamilies().capacity() *
          sizeof(WVCoefficientFamilyLayout) +
      layout.additionalBlocks().capacity() *
          sizeof(WVAdditionalStateBlockLayout) +
      layout.stateBlockRecords().capacity() * sizeof(WVStateBlockRecord) +
      layout.observerRecords().capacity() * sizeof(WVObserverRecord);
  for (const auto &family : layout.coefficientFamilies())
    expectedLayoutStorage += family.identifier.capacity() +
                             family.spectralDimensions.capacity() *
                                 sizeof(std::size_t);
  for (const auto &block : layout.additionalBlocks())
    expectedLayoutStorage += block.identifier.capacity() +
                             block.dimensions.capacity() * sizeof(std::size_t);
  for (const auto &block : layout.stateBlockRecords())
    expectedLayoutStorage += block.identifier.capacity() +
                             block.dimensions.capacity() * sizeof(std::size_t);
  for (const auto &observer : layout.observerRecords()) {
    expectedLayoutStorage +=
        observer.identifier.capacity() + observer.name.capacity() +
        observer.typeIdentifier.capacity() +
        observer.configuration.persistentBytes() -
            sizeof(WVPortableTypedRecord) +
        observer.stateBlockIdentifiers.capacity() * sizeof(std::string) +
        observer.fieldNames.capacity() * sizeof(std::string) +
        (observer.x.capacity() + observer.y.capacity() +
         observer.z.capacity()) *
            sizeof(double);
    for (const auto &identifier : observer.stateBlockIdentifiers)
      expectedLayoutStorage += identifier.capacity();
    for (const auto &field : observer.fieldNames)
      expectedLayoutStorage += field.capacity();
  }
  require(layout.persistentBytes() == expectedLayoutStorage,
          "integration layout exact retained-storage ledger");
  WVIntegrationStateLayout badLayout;
  require(!WVIntegrationStateLayout::create({3, 2}, descriptor, badLayout),
          "coefficient shape mismatch rejected");
}

void testA0OnlyTransformContract() {
  A0OnlyIntegrationSystem system;
  const auto &layout = system.stateLayout();
  require(layout.transformIdentifier() == "WVTransformBarotropicQG" &&
              layout.spatialDimensions() ==
                  std::vector<std::size_t>({8, 6}) &&
              layout.coefficientFamilyCount() == 1 &&
              layout.coefficientFamilies()[0].identifier == "A0" &&
              layout.coefficientFamilies()[0].spectralDimensions ==
                  std::vector<std::size_t>({24}) &&
              !layout.hasLegacyCoefficientTriple(),
          "A0-only transform identity and spatial/spectral ranks");

  WVCoefficientStateStorage storage;
  require(static_cast<bool>(storage.initialize(layout)),
          "A0-only state allocation");
  require(storage.familyCount() == 1 &&
              layout.coefficientElementCount() == 24 &&
              layout.integratedScalarCount() == 48 &&
              storage.capacityBytes() ==
                  24 * sizeof(WVComplex64) +
                      sizeof(WVCoefficientFamilyView) +
                      sizeof(WVCoefficientFamilyConstView),
          "A0-only allocation has no dummy wave families");
  for (std::size_t index = 0; index < 24; ++index)
    storage.mutableFamilies()[0].data[index] =
        {1.0 + static_cast<double>(index), -0.5};
  WVMutableIntegrationState state;
  state.waveVortex.t = 0.0;
  state.waveVortex.t0 = -2.0;
  state.coefficientFamilies = storage.mutableFamilies();
  state.coefficientFamilyCount = storage.familyCount();
  require(state.waveVortex.coefficients.Ap.data == nullptr &&
              state.waveVortex.coefficients.Am.data == nullptr &&
              state.waveVortex.coefficients.A0.data == nullptr,
          "A0-only state retained no legacy state-sized compatibility views");

  WVCoefficientStateStorage fluxStorage;
  require(static_cast<bool>(fluxStorage.initialize(layout)),
          "A0-only RHS allocation");
  WVIntegrationFlux flux;
  flux.coefficientFamilies = fluxStorage.mutableFamilies();
  flux.coefficientFamilyCount = fluxStorage.familyCount();
  std::vector<WVCoefficientFamilyConstView> coefficientViews;
  std::vector<WVAdditionalStateBlockConstView> blockViews;
  const auto constState =
      integrationConstView(state, coefficientViews, blockViews);
  require(static_cast<bool>(system.evaluateRightHandSide(constState, flux)) &&
              fluxStorage.mutableFamilies()[0].data[0].real == -1.0 &&
              fluxStorage.mutableFamilies()[0].data[0].imag == 0.5,
          "A0-only transform RHS execution");

  WVFixedStepRK4 integrator(system, {true});
  auto status = integrator.prepareStateAfterRestart(state);
  require(static_cast<bool>(status),
          "A0-only RK4 restart preparation: " + status.message);
  status = integrator.step(state, 0.1);
  require(static_cast<bool>(status),
          "A0-only RK4 execution: " + status.message);
  require(integrator.metrics().workspaceCapacityBytes ==
              4 * layout.coefficientElementCount() * sizeof(WVComplex64),
          "A0-only RK4 allocated exactly one family per retained workspace");

  WVCoefficientStateStorage denseStorage;
  require(static_cast<bool>(denseStorage.initialize(layout)),
          "A0-only dense-output allocation");
  WVMutableIntegrationState dense;
  dense.coefficientFamilies = denseStorage.mutableFamilies();
  dense.coefficientFamilyCount = denseStorage.familyCount();
  status = integrator.evaluateDenseOutput(0.05, dense);
  require(static_cast<bool>(status),
          "A0-only dense-output execution: " + status.message);
  require(std::abs(denseStorage.mutableFamilies()[0].data[0].real -
                   std::exp(-0.05)) < 3e-6,
          "A0-only dense-output value: " +
              std::to_string(
                  denseStorage.mutableFamilies()[0].data[0].real));

  WVAdaptiveRK23Options rk23Options;
  rk23Options.relativeTolerance = 1e-6;
  rk23Options.maximumStepSize = 0.01;
  WVAdaptiveRK23 rk23(system, rk23Options);
  status = rk23.prepareStateAfterRestart(state);
  require(static_cast<bool>(status),
          "A0-only RK23 restart preparation: " + status.message);
  status = rk23.step(state, 0.01);
  require(static_cast<bool>(status),
          "A0-only RK23 execution: " + status.message);
  require(rk23.metrics().workspaceStateEquivalentCount == 5 &&
              rk23.metrics().acceptedStepCount == 1 &&
              rk23.metrics().rejectedStepCount == 0 &&
              rk23.metrics().rightHandSideEvaluationCount == 4,
          "A0-only RK23 preserved adaptive work and workspace contracts");

  WVAdaptiveRK45Options rk45Options;
  rk45Options.relativeTolerance = 1e-6;
  rk45Options.maximumStepSize = 0.01;
  WVAdaptiveRK45 rk45(system, rk45Options);
  status = rk45.prepareStateAfterRestart(state);
  require(static_cast<bool>(status),
          "A0-only RK45 restart preparation: " + status.message);
  status = rk45.step(state, 0.01);
  require(static_cast<bool>(status),
          "A0-only RK45 execution: " + status.message);
  require(rk45.metrics().workspaceStateEquivalentCount == 7 &&
              rk45.metrics().acceptedStepCount == 1 &&
              rk45.metrics().rejectedStepCount == 0 &&
              rk45.metrics().rightHandSideEvaluationCount == 7 &&
              state.coefficientFamilyCount == 1 &&
              state.waveVortex.coefficients.Ap.data == nullptr &&
              state.waveVortex.coefficients.Am.data == nullptr &&
              state.waveVortex.coefficients.A0.data == nullptr,
          "A0-only RK45 preserved adaptive work and transform-neutral state");

  WVAdaptiveRK78Options rk78Options;
  rk78Options.relativeTolerance = 1e-6;
  rk78Options.maximumStepSize = 0.01;
  WVAdaptiveRK78 rk78(system, rk78Options);
  status = rk78.prepareStateAfterRestart(state);
  require(static_cast<bool>(status),
          "A0-only RK78 restart preparation: " + status.message);
  status = rk78.step(state, 0.01);
  require(static_cast<bool>(status),
          "A0-only RK78 execution: " + status.message);
  require(rk78.metrics().workspaceStateEquivalentCount == 11 &&
              rk78.metrics().acceptedStepCount == 1 &&
              rk78.metrics().rejectedStepCount == 0 &&
              rk78.metrics().rightHandSideEvaluationCount == 13 &&
              state.coefficientFamilyCount == 1 &&
              state.waveVortex.coefficients.Ap.data == nullptr &&
              state.waveVortex.coefficients.Am.data == nullptr &&
              state.waveVortex.coefficients.A0.data == nullptr,
          "A0-only RK78 preserved adaptive work and transform-neutral state");

  coefficientViews.clear();
  blockViews.clear();
  const auto checkpointSource =
      integrationConstView(state, coefficientViews, blockViews);
  WVTransformStateCheckpoint checkpoint;
  require(static_cast<bool>(captureTransformStateCheckpoint(
              layout, checkpointSource, checkpoint)) &&
              checkpoint.coefficientFamilies.size() == 1 &&
              checkpoint.coefficientFamilies[0].identifier == "A0" &&
              checkpoint.coefficientFamilies[0].values.size() == 24,
          "A0-only checkpoint capture");
  WVCoefficientStateStorage restoredStorage;
  WVMutableIntegrationState restored;
  require(static_cast<bool>(restoreTransformStateCheckpoint(
              checkpoint, layout, restoredStorage, restored)) &&
              restored.coefficientFamilyCount == 1 &&
              restored.waveVortex.coefficients.Ap.data == nullptr &&
              restored.waveVortex.coefficients.Am.data == nullptr &&
              restoredStorage.mutableFamilies()[0].data[0].real ==
                  storage.mutableFamilies()[0].data[0].real,
          "A0-only checkpoint round trip without dummy storage");
}

void testRK4(LinearIntegrationSystem &system) {
  StateFixture leanFixture(system.stateLayout());
  WVFixedStepRK4 leanRK4(system, {false});
  leanFixture.state.additionalBlocks[0].realData[0] = -1.0;
  require(
      static_cast<bool>(leanRK4.prepareStateAfterRestart(leanFixture.state)),
      "lean RK4 restart preparation");
  require(leanFixture.state.additionalBlocks[0].realData[0] == 0.0,
          "restart reconstruction applies integration-state constraints");
  StateFixture fixture(system.stateLayout());
  WVFixedStepRK4 rk4(system, {true});
  require(static_cast<bool>(rk4.prepareStateAfterRestart(fixture.state)),
          "RK4 restart preparation");
  require(static_cast<bool>(rk4.step(fixture.state, 0.01)), "RK4 step");
  require(std::abs(fixture.coefficients[0].real - std::exp(-0.01)) < 1e-10,
          "RK4 coefficient result");
  require(std::abs(fixture.state.additionalBlocks[0].realData[0] -
                   std::exp(-0.02)) < 1e-9,
          "RK4 real block result");
  require(static_cast<bool>(rk4.evaluateDenseOutput(0.005, fixture.output)),
          "RK4 integration-state dense output");
  require(std::abs(fixture.outputCoefficients[0].real - std::exp(-0.005)) <
              1e-7,
          "RK4 dense coefficient result");
  require(rk4.metrics().workspaceCapacityBytes > 0 &&
              rk4.metrics().workspaceMaximumLiveBytes ==
                  rk4.metrics().workspaceCapacityBytes &&
              rk4.metrics().workspaceStateEquivalentCount == 4 &&
              rk4.metrics().workspaceMaximumLiveStateEquivalentCount == 4 &&
              rk4.metrics().denseHistoryStateEquivalentCount == 1 &&
              leanRK4.metrics().workspaceStateEquivalentCount == 3 &&
              leanRK4.metrics().workspaceMaximumLiveStateEquivalentCount == 3 &&
              leanRK4.metrics().denseHistoryStateEquivalentCount == 0,
          "RK4 storage accounting");
  require(WVFixedStepRK4::stageBufferLastUseRecordCount() == 4 &&
              std::string(WVFixedStepRK4::stageBufferLastUseRecords()[0]
                              .producer) == "stage-state construction" &&
              std::string(WVFixedStepRK4::stageBufferLastUseRecords()[3]
                              .lastUse) == "dense-output interpolation",
          "RK4 explicit buffer producer and last-consumer schedule");
  require(rk4.persistentBytes() >
              sizeof(rk4) + rk4.metrics().workspaceCapacityBytes,
          "RK4 retained ledger omitted its workspace object or accepted views");
  require(leanRK4.metrics().workspaceCapacityBytes <
              rk4.metrics().workspaceCapacityBytes,
          "RK4 dense history allocated only when requested");
}

void testRK23(LinearIntegrationSystem &system) {
  StateFixture fixture(system.stateLayout());
  WVAdaptiveRK23Options options;
  options.relativeTolerance = 1e-8;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepFactor = 2.0;
  WVAdaptiveRK23 rk23(system, options);
  require(static_cast<bool>(rk23.prepareStateAfterRestart(fixture.state)),
          "RK23 restart preparation");
  require(static_cast<bool>(rk23.step(fixture.state, 0.5)),
          "RK23 adaptive step");
  require(rk23.metrics().rejectedStepCount > 0, "RK23 rejection");
  const auto accepted = rk23.lastAcceptedStep();
  require(accepted && accepted->finalTime > 0, "RK23 accepted step");
  require(accepted->methodStatistics.rejectedStepCount > 0 &&
              std::abs(accepted->methodStatistics.nextStepSize -
                       accepted->methodStatistics.stepSize) < 1e-15,
          "RK23 MATLAB controller does not grow immediately after rejection");
  const auto midpoint = 0.5 * (accepted->initialTime + accepted->finalTime);
  require(static_cast<bool>(rk23.evaluateDenseOutput(midpoint, fixture.output)),
          "RK23 integration-state dense output");
  require(std::abs(fixture.outputCoefficients[0].real - std::exp(-midpoint)) <
              1e-6,
          "RK23 dense coefficient result");
  require(rk23.metrics().workspaceCapacityBytes > 0 &&
              rk23.metrics().workspaceStateEquivalentCount == 5 &&
              rk23.metrics().workspaceMaximumLiveStateEquivalentCount == 5 &&
              WVAdaptiveRK23::stageBufferLastUseRecordCount() == 5 &&
              std::string(WVAdaptiveRK23::stageBufferLastUseRecords()[0]
                              .producer) == "stage-state construction",
          "RK23 exact workspace and buffer-liveness accounting");
  const auto adaptiveArrayAndPolicyBytes =
      sizeof(rk23) + rk23.metrics().workspaceCapacityBytes +
      rk23.metrics().errorPolicyBytes +
      rk23.stepDiagnostics().capacity() *
          sizeof(WVAdaptiveRK23StepDiagnostic) +
      rk23.toleranceComponentHashes().capacity() * sizeof(std::uint64_t);
  require(rk23.persistentBytes() > adaptiveArrayAndPolicyBytes,
          "RK23 retained ledger omitted its workspace object or accepted views");
}

void testRK23MatlabControllerWork() {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {2, 3}, layout)),
          "coefficient-only parity layout");
  LinearIntegrationSystem system(std::move(layout), true);
  StateFixture fixture(system.stateLayout());
  WVAdaptiveRK23Options options;
  options.maximumStepSize = 0.1;
  WVAdaptiveRK23 rk23(system, options);
  require(static_cast<bool>(rk23.prepareStateAfterRestart(fixture.state)),
          "RK23 MATLAB-controller restart preparation");
  require(static_cast<bool>(rk23.advanceToTime(fixture.state, 1.0, 0.5)),
          "RK23 MATLAB-controller interval");
  require(rk23.metrics().acceptedStepCount == 10,
          "RK23 MATLAB-controller maximum-step count");
  require(rk23.metrics().rejectedStepCount == 0,
          "RK23 MATLAB-controller rejection count");
  require(rk23.metrics().rightHandSideEvaluationCount == 31,
          "RK23 MATLAB-controller FSAL work count");
  require(rk23.stepDiagnostics().size() == 10,
          "RK23 MATLAB-controller step diagnostics");
  for (const auto &diagnostic : rk23.stepDiagnostics()) {
    require(std::abs(diagnostic.acceptedStepSize - 0.1) < 1e-14,
            "RK23 MATLAB-controller accepted step size");
    require(diagnostic.rejectedAttemptCount == 0 &&
                diagnostic.normalizedError == 0.0,
            "RK23 MATLAB-controller accepted-step metadata");
  }
  require(!rk23.stepDiagnostics().front().reusedFSALDerivative &&
              rk23.stepDiagnostics().back().reusedFSALDerivative,
          "RK23 MATLAB-controller FSAL metadata");
  require(std::string(WVAdaptiveRK23::controllerIdentifier()) ==
              "matlab-ode23-v1",
          "RK23 MATLAB controller identity");
}

void testRK23MatlabOde23ParityFixture() {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "MATLAB ode23 parity layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK23Options options;
  options.relativeTolerance = 1e-8;
  options.absoluteToleranceScale = 1e-10;
  options.maximumStepSize = 0.5;
  WVAdaptiveRK23 rk23(system, options);
  require(static_cast<bool>(rk23.prepareStateAfterRestart(state)),
          "MATLAB ode23 parity restart preparation");
  require(static_cast<bool>(rk23.advanceToTime(state, 1.0, 0.5)),
          "MATLAB ode23 parity interval");
  require(rk23.metrics().acceptedStepCount == 159 &&
              rk23.metrics().rejectedStepCount == 6 &&
              rk23.metrics().rightHandSideEvaluationCount == 496,
          "MATLAB ode23 parity work counts");
  require(rk23.stepDiagnosticsComplete() &&
              rk23.stepDiagnostics().size() == 159,
          "MATLAB ode23 parity diagnostics");
  require(rk23.toleranceComponentHashes().size() == 3 &&
              rk23.toleranceHash() != 0,
          "MATLAB ode23 parity tolerance audit");
  require(std::abs(rk23.stepDiagnostics().front().acceptedStepSize -
                   0.0078125) < 1e-15,
          "MATLAB ode23 parity first accepted step");
  constexpr double matlabReal = 0.36787943731476008;
  constexpr double matlabImaginary = 0.18393971865738004;
  require(std::abs(coefficients[0].real - matlabReal) < 1e-14 &&
              std::abs(coefficients[0].imag - matlabImaginary) < 1e-14,
          "MATLAB ode23 parity endpoint");
}

void testRK45(LinearIntegrationSystem &system) {
  StateFixture fixture(system.stateLayout());
  WVAdaptiveRK45Options options;
  options.relativeTolerance = 1e-8;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepFactor = 2.0;
  WVAdaptiveRK45 rk45(system, options);
  require(static_cast<bool>(rk45.prepareStateAfterRestart(fixture.state)),
          "RK45 restart preparation");
  require(static_cast<bool>(rk45.step(fixture.state, 0.5)),
          "RK45 adaptive step");
  require(rk45.metrics().rejectedStepCount > 0, "RK45 rejection");
  const auto accepted = rk45.lastAcceptedStep();
  require(accepted && accepted->finalTime > 0.0, "RK45 accepted step");
  require(accepted->methodStatistics.rejectedStepCount > 0 &&
              std::abs(accepted->methodStatistics.nextStepSize -
                       accepted->methodStatistics.stepSize) < 1e-15,
          "RK45 MATLAB controller does not grow immediately after rejection");
  const auto midpoint = 0.5 * (accepted->initialTime + accepted->finalTime);
  require(static_cast<bool>(rk45.evaluateDenseOutput(midpoint, fixture.output)),
          "RK45 integration-state dense output");
  require(std::abs(fixture.outputCoefficients[0].real - std::exp(-midpoint)) <
              1e-9,
          "RK45 dense coefficient result");
  require(rk45.metrics().workspaceStateEquivalentCount == 7 &&
              rk45.metrics().denseHistoryStateEquivalentCount == 6 &&
              rk45.metrics().stateCapacityBytes > 0 &&
              rk45.metrics().workspaceMaximumLiveBytes ==
                  rk45.metrics().workspaceCapacityBytes,
          "RK45 exact workspace and dense-history ledger");
  require(rk45.metrics().diagnosticCapacityBytes ==
              rk45.stepDiagnostics().capacity() *
                  sizeof(WVAdaptiveRK45StepDiagnostic),
          "RK45 exact diagnostic ledger");
  require(WVAdaptiveRK45::stageBufferLastUseRecordCount() == 7 &&
              std::string(WVAdaptiveRK45::stageBufferLastUseRecords()[2]
                              .bufferIdentifier) == "k2/k7" &&
              std::string(WVAdaptiveRK45::stageBufferLastUseRecords()[2]
                              .producer) ==
                  "second-stage then endpoint right-hand side",
          "RK45 explicit stage-buffer liveness schedule");
  require(rk45.persistentBytes() >
              sizeof(rk45) + rk45.metrics().workspaceCapacityBytes +
                  rk45.metrics().errorPolicyBytes +
                  rk45.metrics().diagnosticCapacityBytes,
          "RK45 retained ledger omitted internal records");
}

void testRK45MatlabOde45ParityFixture() {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "MATLAB ode45 parity layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
  std::vector<WVComplex64> denseCoefficients(3);
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVMutableIntegrationState denseState{
      {0.0,
       0.0,
       {{denseCoefficients.data(), {1, 1}},
        {denseCoefficients.data() + 1, {1, 1}},
        {denseCoefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK45Options options;
  options.relativeTolerance = 1e-8;
  options.absoluteToleranceScale = 1e-10;
  options.maximumStepSize = 0.5;
  WVAdaptiveRK45 rk45(system, options);
  require(static_cast<bool>(rk45.prepareStateAfterRestart(state)),
          "MATLAB ode45 parity restart preparation");
  double stepSize = 0.5;
  bool evaluatedDenseOutput = false;
  while (state.waveVortex.t < 1.0) {
    const auto remaining = 1.0 - state.waveVortex.t;
    const auto use = 1.1 * stepSize >= remaining
                         ? remaining
                         : std::min(stepSize, remaining);
    require(static_cast<bool>(rk45.step(state, use)),
            "MATLAB ode45 parity accepted step");
    const auto accepted = rk45.lastAcceptedStep();
    if (!evaluatedDenseOutput && accepted->initialTime <= 0.5 &&
        accepted->finalTime >= 0.5) {
      require(static_cast<bool>(rk45.evaluateDenseOutput(0.5, denseState)),
              "MATLAB ode45 parity dense output");
      evaluatedDenseOutput = true;
    }
    stepSize = rk45.nextStepSize();
  }
  require(rk45.metrics().acceptedStepCount == 13 &&
              rk45.metrics().rejectedStepCount == 1 &&
              rk45.metrics().rightHandSideEvaluationCount == 85,
          "MATLAB ode45 parity work counts");
  require(rk45.stepDiagnosticsComplete() &&
              rk45.stepDiagnostics().size() == 13,
          "MATLAB ode45 parity diagnostics");
  require(std::abs(rk45.stepDiagnostics().front().acceptedStepSize -
                   0.080303422113484) < 2e-14,
          "MATLAB ode45 parity first accepted step");
  constexpr double matlabReal = 0.367879441616479;
  constexpr double matlabImaginary = 0.183939720808240;
  require(std::abs(coefficients[0].real - matlabReal) < 2e-14 &&
              std::abs(coefficients[0].imag - matlabImaginary) < 2e-14,
          "MATLAB ode45 parity endpoint");
  constexpr double matlabDenseReal = 0.606530659859205;
  constexpr double matlabDenseImaginary = 0.303265329929602;
  require(evaluatedDenseOutput &&
              std::abs(denseCoefficients[0].real - matlabDenseReal) < 2e-13 &&
              std::abs(denseCoefficients[0].imag - matlabDenseImaginary) <
                  2e-13,
          "MATLAB ode45 parity continuous extension");
  require(std::string(WVAdaptiveRK45::controllerIdentifier()) ==
              "matlab-ode45-v1",
          "RK45 MATLAB controller identity");
}

double rk45SingleStepError(double stepSize) {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "RK45 order-test layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.0});
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK45Options options;
  options.relativeTolerance = 1.0;
  options.maximumStepSize = stepSize;
  WVAdaptiveRK45 rk45(system, options);
  require(static_cast<bool>(rk45.prepareStateAfterRestart(state)) &&
              static_cast<bool>(rk45.step(state, stepSize)) &&
              rk45.metrics().rejectedStepCount == 0,
          "RK45 order-test step");
  return std::abs(coefficients[0].real - std::exp(-stepSize));
}

void testRK45OrderConstraintsAndSegmentation() {
  const auto coarseError = rk45SingleStepError(0.2);
  const auto fineError = rk45SingleStepError(0.1);
  require(coarseError / fineError > 25.0,
          "RK45 fifth-order accepted solution convergence");

  WVIntegrationStateLayout constrainedLayout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, constrainedLayout)),
          "RK45 constraint layout");
  LinearIntegrationSystem constrainedSystem(std::move(constrainedLayout),
                                            false, false);
  std::vector<WVComplex64> constrainedCoefficients(3, {1.0, 0.0});
  WVMutableIntegrationState constrainedState{
      {0.0,
       0.0,
       {{constrainedCoefficients.data(), {1, 1}},
        {constrainedCoefficients.data() + 1, {1, 1}},
        {constrainedCoefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK45Options constrainedOptions;
  constrainedOptions.relativeTolerance = 1.0;
  WVAdaptiveRK45 constrainedRK45(constrainedSystem, constrainedOptions);
  require(static_cast<bool>(
              constrainedRK45.prepareStateAfterRestart(constrainedState)) &&
              static_cast<bool>(constrainedRK45.step(constrainedState, 0.01)) &&
              static_cast<bool>(constrainedRK45.step(constrainedState, 0.01)),
          "RK45 constrained steps");
  require(constrainedRK45.metrics().fsalReuseCount == 0 &&
              constrainedRK45.metrics().fsalInvalidationCount == 2 &&
              constrainedRK45.metrics().rightHandSideEvaluationCount == 14,
          "RK45 constraint invalidation disables FSAL reuse");

  const auto integrate = [](bool segmented) {
    WVIntegrationStateLayout layout;
    require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
                {1, 1}, layout)),
            "RK45 segmentation layout");
    LinearIntegrationSystem system(std::move(layout));
    std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
    WVMutableIntegrationState state{
        {0.0,
         0.0,
         {{coefficients.data(), {1, 1}},
          {coefficients.data() + 1, {1, 1}},
          {coefficients.data() + 2, {1, 1}}}},
        nullptr,
        0};
    WVAdaptiveRK45Options options;
    options.relativeTolerance = 1e-6;
    options.absoluteToleranceScale = 1e-8;
    options.maximumStepSize = 0.1;
    WVAdaptiveRK45 rk45(system, options);
    require(static_cast<bool>(rk45.prepareStateAfterRestart(state)),
            "RK45 segmentation restart preparation");
    if (segmented)
      require(static_cast<bool>(rk45.advanceToTime(state, 0.5, 0.1)) &&
                  static_cast<bool>(rk45.advanceToTime(state, 1.0, 0.1)),
              "RK45 segmented integration");
    else
      require(static_cast<bool>(rk45.advanceToTime(state, 1.0, 0.1)),
              "RK45 continuous integration");
    return coefficients[0];
  };
  const auto continuous = integrate(false);
  const auto segmented = integrate(true);
  // The segmented run closes its first interval with 0.5 - t, whereas the
  // continuous run can use the stored 0.1 maximum step directly.  Those
  // mathematically equivalent step sizes can differ by one representable
  // value after accumulated-time rounding, so require roundoff-level rather
  // than bitwise agreement across compilers.
  const auto segmentationScale =
      std::max({1.0, std::abs(continuous.real), std::abs(continuous.imag),
                std::abs(segmented.real), std::abs(segmented.imag)});
  const auto segmentationTolerance =
      64.0 * std::numeric_limits<double>::epsilon() * segmentationScale;
  require(std::abs(continuous.real - segmented.real) <=
                  segmentationTolerance &&
              std::abs(continuous.imag - segmented.imag) <=
                  segmentationTolerance,
          "RK45 segmentation preserves the accepted trajectory");
}

void testRK78(LinearIntegrationSystem &system) {
  StateFixture fixture(system.stateLayout());
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1e-10;
  options.absoluteToleranceScale = 1e-10;
  options.maximumStepFactor = 2.0;
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(fixture.state)),
          "RK78 restart preparation");
  require(static_cast<bool>(rk78.step(fixture.state, 0.5)),
          "RK78 adaptive step");
  const auto accepted = rk78.lastAcceptedStep();
  require(accepted != nullptr && accepted->finalTime > 0.0 &&
              accepted->denseOutput == nullptr,
          "RK78 endpoint-only accepted step");
  require(std::abs(fixture.coefficients[0].real -
                   std::exp(-accepted->finalTime)) < 1e-12,
          "RK78 composite-state coefficient result");
  require(std::abs(fixture.state.additionalBlocks[0].realData[0] -
                   std::exp(-2.0 * accepted->finalTime)) < 1e-12,
          "RK78 composite-state real result");
  require(rk78.metrics().workspaceStateEquivalentCount == 11 &&
              rk78.metrics().denseHistoryStateEquivalentCount == 0 &&
              rk78.metrics().denseHistoryCapacityBytes == 0 &&
              rk78.metrics().retainedBaseStageStateEquivalentCount == 0 &&
              rk78.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0 &&
              rk78.metrics()
                      .continuousExtensionWorkspaceStateEquivalentCount == 0 &&
              rk78.metrics().workspaceMaximumLiveBytes ==
                  rk78.metrics().workspaceCapacityBytes &&
              rk78.metrics().workspaceMaximumLiveStateEquivalentCount == 11 &&
              rk78.metrics().errorPolicyBytes == 0,
          "RK78 exact endpoint-only workspace ledger");
  require(rk78.metrics().diagnosticCapacityBytes ==
              rk78.stepDiagnostics().capacity() *
                  sizeof(WVAdaptiveRK78StepDiagnostic),
          "RK78 exact diagnostic ledger");
  require(WVAdaptiveRK78::stageBufferLastUseRecordCount() == 15 &&
              std::string(WVAdaptiveRK78::stageBufferLastUseRecords()[2]
                              .bufferIdentifier) == "k2/k3/k5" &&
              std::string(WVAdaptiveRK78::stageBufferLastUseRecords()[3]
                              .bufferIdentifier) == "k4/k13" &&
              std::string(WVAdaptiveRK78::stageBufferLastUseRecords()[14]
                              .bufferIdentifier) == "k17" &&
              std::string(WVAdaptiveRK78::stageBufferLastUseRecords()[14]
                              .producer) ==
                  "lazy continuous-extension right-hand side",
          "RK78 explicit stage-buffer liveness schedule");
  require(std::string(WVAdaptiveRK78::controllerIdentifier()) ==
                  "matlab-ode78-v1" &&
              std::string(WVAdaptiveRK78::methodIdentifier()) ==
                  "adaptive-rk78",
          "RK78 stable method and controller identities");
  require(rk78.persistentBytes() >
              sizeof(rk78) + rk78.metrics().workspaceCapacityBytes +
                  rk78.metrics().errorPolicyBytes +
                  rk78.metrics().diagnosticCapacityBytes,
          "RK78 retained ledger omitted internal records");
}

double rk78DenseSingleStepError(double stepSize, double theta) {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "RK78 dense order-test layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
  std::vector<WVComplex64> outputCoefficients(3);
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVMutableIntegrationState output{
      {0.0,
       0.0,
       {{outputCoefficients.data(), {1, 1}},
        {outputCoefficients.data() + 1, {1, 1}},
        {outputCoefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1.0;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepSize = stepSize;
  options.retainDenseOutput = true;
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(state)) &&
              static_cast<bool>(rk78.step(state, stepSize)) &&
              static_cast<bool>(
                  rk78.evaluateDenseOutput(theta * stepSize, output)),
          "RK78 dense order-test step");
  return std::abs(outputCoefficients[0].real -
                  std::exp(-theta * stepSize));
}

void testRK78ContinuousExtension(LinearIntegrationSystem &system) {
  StateFixture fixture(system.stateLayout());
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1.0;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepSize = 0.1;
  options.retainDenseOutput = true;
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(fixture.state)) &&
              static_cast<bool>(rk78.step(fixture.state, 0.1)),
          "RK78 dense accepted step");
  const auto *accepted = rk78.lastAcceptedStep();
  require(accepted != nullptr && accepted->denseOutput == &rk78 &&
              rk78.metrics().retainedBaseStageStateEquivalentCount == 8 &&
              rk78.metrics().denseHistoryStateEquivalentCount == 8 &&
              rk78.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0 &&
              rk78.metrics()
                      .continuousExtensionWorkspaceStateEquivalentCount == 0,
          "RK78 dense base-stage retention remains extension-lazy");
  require(static_cast<bool>(rk78.evaluateDenseOutput(0.0, fixture.output)) &&
              fixture.outputCoefficients[0].real == 1.0 &&
              fixture.outputCoefficients[0].imag == 0.5 &&
              rk78.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0,
          "RK78 exact initial endpoint recovery is extension-free");
  require(static_cast<bool>(rk78.evaluateDenseOutput(
              accepted->finalTime, fixture.output)) &&
              fixture.outputCoefficients[0].real ==
                  fixture.coefficients[0].real &&
              fixture.outputCoefficients[0].imag ==
                  fixture.coefficients[0].imag &&
              rk78.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0,
          "RK78 exact final endpoint recovery is extension-free");
  const auto acceptedCoefficient = fixture.coefficients[0];
  const auto aliasStatus =
      rk78.evaluateDenseOutput(0.5 * accepted->finalTime, fixture.state);
  require(!aliasStatus &&
              aliasStatus.code == WVKernelStatusCode::invalidConfiguration &&
              fixture.coefficients[0].real == acceptedCoefficient.real &&
              fixture.coefficients[0].imag == acceptedCoefficient.imag &&
              rk78.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0,
          "RK78 interpolated output cannot alias accepted integration state");
  const auto firstInterior = 0.31 * accepted->finalTime;
  require(static_cast<bool>(
              rk78.evaluateDenseOutput(firstInterior, fixture.output)) &&
              std::abs(fixture.outputCoefficients[0].real -
                       std::exp(-firstInterior)) < 1e-12 &&
              std::abs(fixture.outputCoefficients[0].imag -
                       0.5 * std::exp(-firstInterior)) < 1e-12 &&
              std::abs(fixture.output.additionalBlocks[0].realData[0] -
                       std::exp(-2.0 * firstInterior)) < 1e-12,
          "RK78 real, complex, and composite interior state: coefficient=" +
              std::to_string(fixture.outputCoefficients[0].real) + "+" +
              std::to_string(fixture.outputCoefficients[0].imag) +
              "i realBlock=" +
              std::to_string(
                  fixture.output.additionalBlocks[0].realData[0]));
  const auto metricsAfterBuild = rk78.metrics();
  require(metricsAfterBuild.continuousExtensionRightHandSideEvaluationCount ==
                  4 &&
              metricsAfterBuild.denseOutputCacheBuildCount == 1 &&
              metricsAfterBuild.denseOutputCacheReuseCount == 0 &&
              metricsAfterBuild.continuousExtensionWorkspaceStateEquivalentCount ==
                  4 &&
              metricsAfterBuild.workspaceMaximumLiveStateEquivalentCount == 15 &&
              metricsAfterBuild.workspaceMaximumLiveBytes >
                  metricsAfterBuild.workspaceCapacityBytes,
          "RK78 lazy extension build and maximum-live ledger");
  const auto secondInterior = 0.73 * accepted->finalTime;
  require(static_cast<bool>(
              rk78.evaluateDenseOutput(secondInterior, fixture.output)) &&
              rk78.metrics().continuousExtensionRightHandSideEvaluationCount ==
                  4 &&
              rk78.metrics().denseOutputCacheBuildCount == 1 &&
              rk78.metrics().denseOutputCacheReuseCount == 1,
          "RK78 multiple interior samples reuse one extension cache");
  require(static_cast<bool>(rk78.step(fixture.state, 0.1)) &&
              rk78.metrics()
                      .continuousExtensionWorkspaceStateEquivalentCount == 0 &&
              rk78.metrics().workspaceLiveBytes ==
                  rk78.metrics().workspaceCapacityBytes,
          "RK78 releases extension-only state before advancing");

  const auto coarse = rk78DenseSingleStepError(0.8, 0.37);
  const auto fine = rk78DenseSingleStepError(0.4, 0.37);
  require(coarse / fine > 150.0,
          "RK78 seventh-order continuous-extension convergence: " +
              std::to_string(coarse / fine));
}

void testRK78MatlabOde78ParityFixture() {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "MATLAB ode78 parity layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1e-12;
  options.absoluteToleranceScale = 1e-10;
  options.maximumStepSize = 0.5;
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(state)),
          "MATLAB ode78 parity restart preparation");
  require(static_cast<bool>(rk78.advanceToTime(state, 1.0, 0.5)),
          "MATLAB ode78 parity interval");
  require(rk78.metrics().acceptedStepCount == 4 &&
              rk78.metrics().rejectedStepCount == 1 &&
              rk78.metrics().rightHandSideEvaluationCount == 64,
          "MATLAB ode78 parity work counts");
  require(rk78.stepDiagnosticsComplete() &&
              rk78.stepDiagnostics().size() == 4,
          "MATLAB ode78 parity diagnostics");
  require(std::abs(rk78.stepDiagnostics().front().acceptedStepSize -
                   0.292534500487804) < 1e-5,
          "MATLAB ode78 parity first accepted step: " +
              std::to_string(
                  rk78.stepDiagnostics().front().acceptedStepSize));
  constexpr double matlabReal = 0.367879441171209;
  constexpr double matlabImaginary = 0.183939720585604;
  require(std::abs(coefficients[0].real - matlabReal) < 1e-12 &&
              std::abs(coefficients[0].imag - matlabImaginary) < 1e-12,
          "MATLAB ode78 parity endpoint");
  require(!rk78.stepDiagnostics().front().reusedFSALDerivative &&
              !rk78.stepDiagnostics().back().reusedFSALDerivative &&
              rk78.metrics().fsalReuseCount == 0,
          "RK78 correctly avoids unsafe endpoint derivative reuse");
}

struct RK78RequestedTimeEvidence {
  WVComplex64 endpoint;
  std::vector<WVComplex64> requestedValues;
  std::vector<WVAdaptiveRK78StepDiagnostic> diagnostics;
  WVIntegratorMetrics metrics;
};

RK78RequestedTimeEvidence
runRK78RequestedTimes(const std::vector<double> &requestedTimes) {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "RK78 requested-time layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.0});
  std::vector<WVComplex64> outputCoefficients(3);
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVMutableIntegrationState output{
      {0.0,
       0.0,
       {{outputCoefficients.data(), {1, 1}},
        {outputCoefficients.data() + 1, {1, 1}},
        {outputCoefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1e-12;
  options.absoluteToleranceScale = 1e-10;
  options.maximumStepSize = 0.5;
  options.retainDenseOutput = !requestedTimes.empty();
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(state)),
          "RK78 requested-time preparation");
  std::vector<WVComplex64> values;
  std::size_t nextOutput = 0;
  if (!requestedTimes.empty() && requestedTimes.front() == 0.0) {
    values.push_back(coefficients[0]);
    ++nextOutput;
  }
  double stepSize = 0.5;
  while (state.waveVortex.t < 1.0) {
    const auto remaining = 1.0 - state.waveVortex.t;
    const auto use = 1.1 * stepSize >= remaining
                         ? remaining
                         : std::min(stepSize, remaining);
    require(static_cast<bool>(rk78.step(state, use)),
            "RK78 requested-time accepted step");
    const auto *accepted = rk78.lastAcceptedStep();
    require(accepted != nullptr, "RK78 requested-time accepted history");
    while (nextOutput < requestedTimes.size() &&
           requestedTimes[nextOutput] <=
               accepted->finalTime +
                   8.0 * std::numeric_limits<double>::epsilon()) {
      const auto requested = requestedTimes[nextOutput];
      if (std::abs(requested - accepted->finalTime) <=
          8.0 * std::numeric_limits<double>::epsilon())
        values.push_back(coefficients[0]);
      else {
        require(static_cast<bool>(
                    rk78.evaluateDenseOutput(requested, output)),
                "RK78 requested-time interpolation");
        values.push_back(outputCoefficients[0]);
      }
      ++nextOutput;
    }
    stepSize = rk78.nextStepSize();
  }
  require(nextOutput == requestedTimes.size(),
          "RK78 requested-time fixture emitted every request");
  return {coefficients[0], values, rk78.stepDiagnostics(), rk78.metrics()};
}

void testRK78MatlabContinuousParityAndTrajectoryIndependence() {
  const std::vector<double> matlabTimes{0.0, 0.1, 0.2, 0.4, 0.7, 1.0};
  const double matlabValues[] = {1.0,
                                 0.904837418056636,
                                 0.818730753079412,
                                 0.670320046055729,
                                 0.496585303810493,
                                 0.367879441171218};
  const auto matlabFixture = runRK78RequestedTimes(matlabTimes);
  require(matlabFixture.requestedValues.size() == matlabTimes.size(),
          "RK78 MATLAB requested-time fixture size");
  for (std::size_t index = 0; index < matlabTimes.size(); ++index)
    require(std::abs(matlabFixture.requestedValues[index].real -
                     matlabValues[index]) <= 1e-12,
            "RK78 MATLAB ode78 requested-time parity at index " +
                std::to_string(index));
  require(matlabFixture.metrics.acceptedStepCount == 4 &&
              matlabFixture.metrics.rejectedStepCount == 1 &&
              matlabFixture.metrics.baseRightHandSideEvaluationCount == 64 &&
              matlabFixture.metrics
                      .continuousExtensionRightHandSideEvaluationCount == 12 &&
              matlabFixture.metrics.rightHandSideEvaluationCount == 76 &&
              matlabFixture.metrics.denseOutputCacheBuildCount == 3 &&
              matlabFixture.metrics.denseOutputCacheReuseCount == 1,
          "RK78 MATLAB ode78 controller, RHS, and cache parity");

  const auto endpointOnly = runRK78RequestedTimes({});
  const auto regrouped = runRK78RequestedTimes(
      {0.05, 0.1, 0.2, 0.21, 0.4, 0.7, 0.9, 1.0});
  const auto sameTrajectory = [](const RK78RequestedTimeEvidence &first,
                                 const RK78RequestedTimeEvidence &second) {
    if (first.endpoint.real != second.endpoint.real ||
        first.endpoint.imag != second.endpoint.imag ||
        first.diagnostics.size() != second.diagnostics.size())
      return false;
    for (std::size_t index = 0; index < first.diagnostics.size(); ++index) {
      const auto &a = first.diagnostics[index];
      const auto &b = second.diagnostics[index];
      if (a.initialTime != b.initialTime ||
          a.acceptedStepSize != b.acceptedStepSize ||
          a.normalizedError != b.normalizedError ||
          a.nextStepSize != b.nextStepSize ||
          a.rejectedAttemptCount != b.rejectedAttemptCount ||
          a.rightHandSideEvaluationCount !=
              b.rightHandSideEvaluationCount)
        return false;
    }
    return true;
  };
  require(sameTrajectory(endpointOnly, matlabFixture) &&
              sameTrajectory(matlabFixture, regrouped) &&
              endpointOnly.metrics
                      .continuousExtensionRightHandSideEvaluationCount == 0 &&
              endpointOnly.metrics
                      .continuousExtensionWorkspaceMaximumLiveBytes == 0,
          "RK78 adding, removing, or regrouping output times preserves the "
          "accepted trajectory");
}

struct RK78ConvergenceEvidence {
  double acceptedError = 0.0;
  double embeddedErrorEstimate = 0.0;
};

RK78ConvergenceEvidence rk78SingleStepEvidence(double stepSize) {
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, layout)),
          "RK78 order-test layout");
  LinearIntegrationSystem system(std::move(layout));
  std::vector<WVComplex64> coefficients(3, {1.0, 0.0});
  WVMutableIntegrationState state{
      {0.0,
       0.0,
       {{coefficients.data(), {1, 1}},
        {coefficients.data() + 1, {1, 1}},
        {coefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1.0;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepSize = stepSize;
  WVAdaptiveRK78 rk78(system, options);
  require(static_cast<bool>(rk78.prepareStateAfterRestart(state)) &&
              static_cast<bool>(rk78.step(state, stepSize)) &&
              rk78.metrics().rejectedStepCount == 0,
          "RK78 order-test step");
  return {std::abs(coefficients[0].real - std::exp(-stepSize)),
          rk78.metrics().normalizedError};
}

void testRK78OrdersFailuresSegmentationAndRestart() {
  const auto coarse = rk78SingleStepEvidence(0.8);
  const auto fine = rk78SingleStepEvidence(0.4);
  require(coarse.acceptedError / fine.acceptedError > 350.0,
          "RK78 eighth-order accepted-solution convergence");
  require(coarse.embeddedErrorEstimate / fine.embeddedErrorEstimate > 180.0,
          "RK78 seventh-order embedded-estimate convergence");

  WVIntegrationStateLayout failureLayout;
  require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
              {1, 1}, failureLayout)),
          "RK78 failure layout");
  LinearIntegrationSystem failureSystem(std::move(failureLayout), false, true,
                                        true);
  std::vector<WVComplex64> failureCoefficients(3, {1.0, 0.0});
  WVMutableIntegrationState failureState{
      {1.0,
       0.0,
       {{failureCoefficients.data(), {1, 1}},
        {failureCoefficients.data() + 1, {1, 1}},
        {failureCoefficients.data() + 2, {1, 1}}}},
      nullptr,
      0};
  WVAdaptiveRK78 failureRK78(failureSystem);
  require(static_cast<bool>(failureRK78.prepareStateAfterRestart(failureState)),
          "RK78 failure restart preparation");
  const auto minimumStep =
      16.0 * (std::nextafter(1.0, std::numeric_limits<double>::infinity()) -
              1.0);
  const auto failureStatus = failureRK78.step(failureState, minimumStep);
  require(!failureStatus &&
              failureStatus.code == WVKernelStatusCode::numericalFailure &&
              failureRK78.metrics().normalizedError ==
                  std::numeric_limits<double>::infinity() &&
              failureRK78.metrics().rejectedStepCount == 1 &&
              failureRK78.metrics().rightHandSideEvaluationCount == 13 &&
              failureState.waveVortex.t == 1.0,
          "RK78 minimum-step and nonfinite-error failure contract");

  const auto integrate = [](bool segmented, bool restartAtMidpoint) {
    WVIntegrationStateLayout layout;
    require(static_cast<bool>(WVIntegrationStateLayout::createCoefficientOnly(
                {1, 1}, layout)),
            "RK78 segmentation layout");
    LinearIntegrationSystem system(std::move(layout));
    std::vector<WVComplex64> coefficients(3, {1.0, 0.5});
    WVMutableIntegrationState state{
        {0.0,
         0.0,
         {{coefficients.data(), {1, 1}},
          {coefficients.data() + 1, {1, 1}},
          {coefficients.data() + 2, {1, 1}}}},
        nullptr,
        0};
    WVAdaptiveRK78Options options;
    options.relativeTolerance = 1e-8;
    options.absoluteToleranceScale = 1e-10;
    options.maximumStepSize = 0.1;
    WVAdaptiveRK78 first(system, options);
    require(static_cast<bool>(first.prepareStateAfterRestart(state)),
            "RK78 segmentation restart preparation");
    if (!segmented) {
      require(static_cast<bool>(first.advanceToTime(state, 1.0, 0.1)),
              "RK78 continuous integration");
    } else {
      require(static_cast<bool>(first.advanceToTime(state, 0.5, 0.1)),
              "RK78 first segment");
      if (restartAtMidpoint) {
        WVAdaptiveRK78 reconstructed(system, options);
        require(static_cast<bool>(reconstructed.prepareStateAfterRestart(state)) &&
                    reconstructed.lastAcceptedStep() == nullptr &&
                    reconstructed.stepDiagnostics().empty() &&
                    static_cast<bool>(
                        reconstructed.advanceToTime(state, 1.0, 0.1)),
                "RK78 restart reconstruction");
      } else {
        require(static_cast<bool>(first.advanceToTime(state, 1.0, 0.1)),
                "RK78 second segment");
      }
    }
    return coefficients[0];
  };
  const auto continuous = integrate(false, false);
  const auto segmented = integrate(true, false);
  const auto restarted = integrate(true, true);
  const auto tolerance =
      128.0 * std::numeric_limits<double>::epsilon();
  require(std::abs(continuous.real - segmented.real) <= tolerance &&
              std::abs(continuous.imag - segmented.imag) <= tolerance &&
              std::abs(segmented.real - restarted.real) <= tolerance &&
              std::abs(segmented.imag - restarted.imag) <= tolerance,
          "RK78 segmentation and restart reconstruction preserve trajectory");
}

} // namespace

int main() {
  WVPortableObserverDescriptor descriptor;
  WVIntegrationStateLayout layout;
  testContracts(descriptor, layout);
  LinearIntegrationSystem system(std::move(layout));
  testRK4(system);
  testRK23(system);
  testRK45(system);
  testRK23MatlabControllerWork();
  testRK23MatlabOde23ParityFixture();
  testRK45MatlabOde45ParityFixture();
  testRK45OrderConstraintsAndSegmentation();
  testRK78(system);
  testRK78ContinuousExtension(system);
  testRK78MatlabOde78ParityFixture();
  testRK78MatlabContinuousParityAndTrajectoryIndependence();
  testRK78OrdersFailuresSegmentationAndRestart();
  testA0OnlyTransformContract();
  std::cout <<
      "PASS: portable observer contracts and unified RK4/RK23/RK45/RK78 "
      "integration\n";
  return 0;
}
