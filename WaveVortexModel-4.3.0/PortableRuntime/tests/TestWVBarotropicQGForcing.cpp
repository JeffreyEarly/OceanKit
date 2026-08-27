#include "WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

WVTransformBarotropicQGConfiguration configuration(
    std::size_t Nx = 8, std::size_t Ny = 6, std::uint32_t j = 1,
    bool shouldAntialias = true) {
  WVTransformBarotropicQGConfiguration value;
  value.Nx = Nx;
  value.Ny = Ny;
  value.Lx = 17000.0;
  value.Ly = 11000.0;
  value.h = 0.8;
  value.j = j;
  value.g = 9.81;
  value.planetaryRadius = 6.371e6;
  value.rotationRate = 7.2921e-5;
  value.latitude = 33.0;
  value.shouldAntialias = shouldAntialias;
  return value;
}

WVPortableTypedRecord emptyConfiguration() {
  return {"wave-vortex-forcing-configuration-v1", 1, {}};
}

WVPortableNamedValue realValue(std::string name, std::vector<double> values) {
  const auto count = values.size();
  return {std::move(name),
          count == 1 ? std::vector<std::size_t>{}
                     : std::vector<std::size_t>{count},
          std::move(values)};
}

WVPortableNamedValue integerValue(std::string name,
                                  std::vector<std::int64_t> values) {
  const auto count = values.size();
  return {std::move(name), {count}, std::move(values)};
}

WVPortableTypedRecord scalarConfiguration(const char *name, double value) {
  auto result = emptyConfiguration();
  result.values.push_back(realValue(name, {value}));
  return result;
}

WVPortableTypedRecord fixedConfiguration(
    std::vector<std::int64_t> indices,
    std::vector<WVComplex64> values) {
  auto result = emptyConfiguration();
  std::vector<double> real(values.size());
  std::vector<double> imag(values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    real[index] = values[index].real;
    imag[index] = values[index].imag;
  }
  result.values.push_back(integerValue("A0Indices", std::move(indices)));
  result.values.push_back(realValue("A0ValuesReal", std::move(real)));
  result.values.push_back(realValue("A0ValuesImag", std::move(imag)));
  return result;
}

WVFrozenForcingEntry entry(const char *identifier, const char *name,
                           WVForcingStage stage, std::uint8_t priority,
                           std::size_t ordinal,
                           WVPortableTypedRecord values =
                               emptyConfiguration()) {
  return {identifier, WVPortablePairContractVersion, name, stage, priority,
          ordinal, "test", std::move(values)};
}

WVFrozenForcingSchedule schedule(
    std::vector<WVFrozenForcingEntry> entries) {
  WVFrozenForcingSchedule result;
  result.entries = std::move(entries);
  return result;
}

std::shared_ptr<const WVExtensionCatalog> catalog() {
  std::shared_ptr<const WVExtensionCatalog> result;
  const auto status = makeBuiltInExtensionCatalog(result);
  require(static_cast<bool>(status), status.message);
  return result;
}

std::vector<WVComplex64> deterministicState(
    const WVTransformBarotropicQGDescriptor &descriptor) {
  std::vector<WVComplex64> values(descriptor.Nkl());
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index] = {
        2e-5 * std::sin(0.31 * static_cast<double>(index + 1)),
        1e-5 * std::cos(0.17 * static_cast<double>(index + 3))};
    const auto &mode = descriptor.fourierModes()[index];
    if (mode.Kh == 0.0)
      values[index] = {};
    if (mode.dftPrimaryIndex == mode.dftConjugateIndex)
      values[index].imag = 0.0;
  }
  return values;
}

std::unique_ptr<WVBarotropicQGForcingEngine> createEngine(
    const WVTransformBarotropicQGConfiguration &value,
    const WVFrozenForcingSchedule &forcingSchedule,
    std::unique_ptr<WVFFTEngine> fft =
        std::make_unique<WVReferenceFFTEngine>()) {
  std::unique_ptr<WVBarotropicQGForcingEngine> result;
  const auto status = WVBarotropicQGForcingEngine::create(
      value, forcingSchedule, catalog(), std::move(fft), result);
  require(static_cast<bool>(status), status.message);
  return result;
}

std::vector<WVComplex64> evaluate(
    WVBarotropicQGForcingEngine &engine,
    const std::vector<WVComplex64> &A0) {
  std::vector<WVComplex64> result(A0.size(),
                                  {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()});
  WVComplexConstView input{A0.data(), engine.kernel().descriptor().spectralShape()};
  WVComplexView output{result.data(), engine.kernel().descriptor().spectralShape()};
  const auto status = engine.evaluateRightHandSide(input, output);
  require(static_cast<bool>(status), status.message);
  for (const auto value : result)
    require(std::isfinite(value.real) && std::isfinite(value.imag),
            "RHS did not overwrite every compact A0 tendency");
  return result;
}

double maximumRelativeError(const std::vector<WVComplex64> &actual,
                            const std::vector<WVComplex64> &expected) {
  require(actual.size() == expected.size(), "coefficient length mismatch");
  double error = 0.0;
  double scale = 0.0;
  for (std::size_t index = 0; index < actual.size(); ++index) {
    error = std::max(error, std::hypot(actual[index].real - expected[index].real,
                                       actual[index].imag - expected[index].imag));
    scale = std::max(scale,
                     std::hypot(expected[index].real, expected[index].imag));
  }
  return scale == 0.0 ? error : error / scale;
}

std::vector<WVComplex64> add(const std::vector<WVComplex64> &left,
                             const std::vector<WVComplex64> &right) {
  require(left.size() == right.size(), "sum length mismatch");
  auto result = left;
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index].real += right[index].real;
    result[index].imag += right[index].imag;
  }
  return result;
}

std::vector<double> adaptiveDamping(
    const WVTransformBarotropicQGDescriptor &descriptor) {
  double maximumComponent = 0.0;
  for (const auto &mode : descriptor.fourierModes())
    maximumComponent = std::max(
        maximumComponent, std::max(std::abs(mode.k), std::abs(mode.l)));
  const auto &value = descriptor.configuration();
  const double effectiveResolution = pi / maximumComponent;
  const double dklMinimum =
      std::min(2.0 * pi / value.Lx, 2.0 * pi / value.Ly);
  const double cutoff =
      dklMinimum * std::pow(maximumComponent / dklMinimum, 0.75);
  std::vector<double> result(descriptor.Nkl());
  for (std::size_t index = 0; index < result.size(); ++index) {
    const auto &mode = descriptor.fourierModes()[index];
    double filter = 0.0;
    if (mode.Kh > maximumComponent)
      filter = 1.0;
    else if (mode.Kh >= cutoff) {
      const double ratio =
          (mode.Kh - maximumComponent) / (mode.Kh - cutoff);
      filter = std::exp(-(ratio * ratio));
    }
    result[index] = -effectiveResolution / (pi * pi) * filter *
                    (mode.k * mode.k + mode.l * mode.l);
  }
  return result;
}

WVFrozenForcingEntry nonlinear(std::size_t ordinal = 0) {
  return entry("WVNonlinearAdvection", "nonlinear", WVForcingStage::spatial,
               127, ordinal);
}

WVFrozenForcingEntry damping(std::size_t ordinal = 0) {
  return entry("WVAdaptiveDamping", "damping", WVForcingStage::spectral, 255,
               ordinal);
}

WVFrozenForcingEntry linear(double r = 2.5e-7, std::size_t ordinal = 0) {
  return entry("WVBottomFrictionLinear", "linear", WVForcingStage::spatial,
               220, ordinal, scalarConfiguration("r", r));
}

WVFrozenForcingEntry quadratic(double Cd = 1.5e-3,
                               std::size_t ordinal = 0) {
  return entry("WVBottomFrictionQuadratic", "quadratic",
               WVForcingStage::spatial, 210, ordinal,
               scalarConfiguration("Cd", Cd));
}

WVFrozenForcingEntry beta(std::size_t ordinal = 0) {
  return entry("WVBetaPlanePVAdvection", "beta", WVForcingStage::spatial,
               230, ordinal);
}

void testNumericalForcingMatrix() {
  const std::vector<std::pair<std::size_t, std::size_t>> grids = {
      {8, 6}, {9, 7}, {8, 7}, {9, 6}};
  for (const auto &grid : grids) {
    for (std::uint32_t j = 0; j <= 1; ++j) {
      for (const bool shouldAntialias : {false, true}) {
        const auto value =
            configuration(grid.first, grid.second, j, shouldAntialias);
        auto nonlinearEngine = createEngine(value, schedule({nonlinear()}));
        const auto A0 = deterministicState(nonlinearEngine->kernel().descriptor());
        const auto nonlinearActual = evaluate(*nonlinearEngine, A0);
        std::vector<WVComplex64> nonlinearExpected(A0.size());
        WVComplexConstView input{A0.data(),
                                 nonlinearEngine->kernel().descriptor().spectralShape()};
        WVComplexView nonlinearOutput{
            nonlinearExpected.data(),
            nonlinearEngine->kernel().descriptor().spectralShape()};
        require(static_cast<bool>(nonlinearEngine->kernel().nonlinearFlux(
                    input, nonlinearOutput)),
                "direct nonlinear kernel failed");
        require(maximumRelativeError(nonlinearActual, nonlinearExpected) <= 1e-13,
                "resolved nonlinear forcing changed #280 output");

        auto dampingEngine = createEngine(value, schedule({damping()}));
        const auto dampingActual = evaluate(*dampingEngine, A0);
        double uvMax = 0.0;
        require(static_cast<bool>(dampingEngine->kernel().uvMax(input, uvMax)),
                "uvMax failed");
        const auto dampingOperator =
            adaptiveDamping(dampingEngine->kernel().descriptor());
        std::vector<WVComplex64> dampingExpected(A0.size());
        for (std::size_t index = 0; index < A0.size(); ++index) {
          const double factor = uvMax * dampingOperator[index];
          dampingExpected[index] =
              {factor * A0[index].real, factor * A0[index].imag};
        }
        require(maximumRelativeError(dampingActual, dampingExpected) <= 2e-15,
                "adaptive damping operator/tendency mismatch");

        constexpr double r = 2.5e-7;
        auto linearEngine = createEngine(value, schedule({linear(r)}));
        const auto linearActual = evaluate(*linearEngine, A0);
        std::vector<WVComplex64> linearExpected(A0.size());
        const auto &zeta = linearEngine->kernel().descriptor().modes().zetaZFactor;
        for (std::size_t index = 0; index < A0.size(); ++index)
          linearExpected[index] = {-r * zeta[index] * A0[index].real,
                                   -r * zeta[index] * A0[index].imag};
        require(maximumRelativeError(linearActual, linearExpected) <= 1e-12,
                "linear bottom-friction tendency mismatch");

        auto betaEngine = createEngine(value, schedule({beta()}));
        const auto betaActual = evaluate(*betaEngine, A0);
        const double betaValue =
            2.0 * value.rotationRate *
            std::cos(value.latitude * pi / 180.0) / value.planetaryRadius;
        const auto &vFactor = betaEngine->kernel().descriptor().modes().vFactor;
        std::vector<WVComplex64> betaExpected(A0.size());
        for (std::size_t index = 0; index < A0.size(); ++index) {
          const WVComplex64 product{
              vFactor[index].real * A0[index].real -
                  vFactor[index].imag * A0[index].imag,
              vFactor[index].real * A0[index].imag +
                  vFactor[index].imag * A0[index].real};
          betaExpected[index] =
              {-betaValue * product.real, -betaValue * product.imag};
        }
        require(maximumRelativeError(betaActual, betaExpected) <= 1e-12,
                "beta-plane sign/units/tendency mismatch");

        auto quadraticEngine = createEngine(value, schedule({quadratic()}));
        const auto quadraticActual = evaluate(*quadraticEngine, A0);
        std::vector<WVComplex64> quadraticExpected(A0.size());
        WVComplexView quadraticOutput{
            quadraticExpected.data(),
            quadraticEngine->kernel().descriptor().spectralShape()};
        WVBarotropicQGOperationWorkspace workspace;
        require(static_cast<bool>(
                    quadraticEngine->kernel().addQuadraticBottomFriction(
                        input, 1.5e-3 / 4000.0, quadraticOutput, false,
                        workspace)),
                "direct quadratic bottom-friction operation failed");
        require(maximumRelativeError(quadraticActual, quadraticExpected) <=
                    1e-13,
                "quadratic bottom-friction factory scaling mismatch");

        auto composed = createEngine(
            value,
            schedule({damping(7), beta(6), linear(r, 5), nonlinear(9),
                      quadratic(1.5e-3, 4)}));
        const auto composedActual = evaluate(*composed, A0);
        auto composedExpected = add(nonlinearActual, dampingActual);
        composedExpected = add(composedExpected, linearActual);
        composedExpected = add(composedExpected, quadraticActual);
        composedExpected = add(composedExpected, betaActual);
        const auto compositionError =
            maximumRelativeError(composedActual, composedExpected);
        require(compositionError <= 2e-12,
                "ordered forcing composition is not additive: " +
                    std::to_string(compositionError));
        const auto &metrics = composed->metrics();
        require(metrics.evaluationCount == 1 && metrics.forcingCallCount == 5 &&
                    metrics.physicalFieldReconstructionCount == 1 &&
                    metrics.physicalFieldReuseCount == 4 &&
                    metrics.spatialTendencyProjectionCount == 4 &&
                    metrics.workspaceCapacityBytes == 0,
                "RHS-scoped field reuse/projection metrics mismatch");
        require(composed->scheduleIdentifier() ==
                    "wave-vortex-forcing-v1:WVNonlinearAdvection,"
                    "WVBottomFrictionQuadratic,WVBottomFrictionLinear,"
                    "WVBetaPlanePVAdvection,WVAdaptiveDamping",
                "stage/priority/original-ordinal order mismatch");
        const auto &kernelMetrics = composed->kernel().metrics();
        const auto halfRows =
            composed->kernel().descriptor().halfSpectrumMappings().NxHalf *
            value.Ny;
        const auto exactScratch =
            4 * halfRows * sizeof(WVComplex64) +
            5 * value.Nx * value.Ny * sizeof(double);
        require(composed->kernel().scratchBytes() == exactScratch &&
                    kernelMetrics.planCount == 3 &&
                    kernelMetrics.persistentFullHermitianBytes == 0 &&
                    composed->kernel().engineIdentifier() ==
                        "reference-direct",
                "compact provider/storage invariants changed");
      }
    }
  }
}

void testFixedAmplitudeAndNarrowBand() {
  const auto value = configuration(9, 6, 1, true);
  auto base = createEngine(value, schedule({nonlinear()}));
  const auto A0 = deterministicState(base->kernel().descriptor());
  const std::vector<std::int64_t> indices = {1, 3};
  const std::vector<WVComplex64> fixedValues = {{3e-6, -2e-6}, {-4e-6, 0.0}};
  auto fixed = entry("WVFixedAmplitudeForcing", "fixed",
                     WVForcingStage::spectralAmplitude, 255, 8,
                     fixedConfiguration(indices, fixedValues));
  auto engine = createEngine(value, schedule({fixed, nonlinear(3)}));
  auto tendency = evaluate(*engine, A0);
  require(tendency[1].real == 0.0 && tendency[1].imag == 0.0 &&
              tendency[3].real == 0.0 && tendency[3].imag == 0.0,
          "fixed-amplitude stage did not zero selected tendencies");
  auto constrained = A0;
  WVComplexView view{constrained.data(), engine->kernel().descriptor().spectralShape()};
  const auto result = engine->restoreForcingAmplitudes(view);
  require(static_cast<bool>(result) && result.modifiedCoefficientCount == 2 &&
              !result.fsalCompatible &&
              constrained[1].real == fixedValues[0].real &&
              constrained[1].imag == fixedValues[0].imag &&
              constrained[3].real == fixedValues[1].real &&
              constrained[3].imag == fixedValues[1].imag &&
              engine->metrics().constraintOperationCount == 1 &&
              engine->metrics().restoredCoefficientCount == 2 &&
              engine->metrics().stateConstraintElementWrites == 2,
          "fixed-amplitude restoration/metrics mismatch");

  auto narrow = entry("WVNarrowBandGeostrophicForcing", "narrow",
                      WVForcingStage::spectralAmplitude, 255, 2,
                      fixedConfiguration(indices, fixedValues));
  auto narrowEngine = createEngine(value, schedule({nonlinear(), narrow}));
  const auto narrowTendency = evaluate(*narrowEngine, A0);
  require(narrowTendency[1].real == 0.0 && narrowTendency[1].imag == 0.0 &&
              narrowTendency[3].real == 0.0 &&
              narrowTendency[3].imag == 0.0,
          "narrow-band record did not resolve through fixed-amplitude contract");
}

struct ProviderCounts {
  std::size_t plansCreated = 0;
  std::size_t plansDestroyed = 0;
  std::size_t executions = 0;
  std::size_t failAtExecution = 0;
};

class CountingPlan final : public WVFFTPlan {
public:
  CountingPlan(std::unique_ptr<WVFFTPlan> inner,
               std::shared_ptr<ProviderCounts> counts)
      : inner_(std::move(inner)), counts_(std::move(counts)) {}
  ~CountingPlan() override { ++counts_->plansDestroyed; }
  WVKernelStatus execute(const void *input, void *output) override {
    ++counts_->executions;
    if (counts_->failAtExecution != 0 &&
        counts_->executions == counts_->failAtExecution)
      return {WVKernelStatusCode::fftExecutionFailure,
              "injected Barotropic QG forcing FFT failure"};
    return inner_->execute(input, output);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + inner_->persistentBytes();
  }

private:
  std::unique_ptr<WVFFTPlan> inner_;
  std::shared_ptr<ProviderCounts> counts_;
};

class CountingEngine final : public WVFFTEngine {
public:
  explicit CountingEngine(std::shared_ptr<ProviderCounts> counts)
      : counts_(std::move(counts)) {}
  std::string identifier() const override { return "reference-counting"; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + reference_.persistentBytes();
  }
  WVKernelStatus createPlan(const WVFFTPlanSpecification &specification,
                            std::unique_ptr<WVFFTPlan> &plan) override {
    std::unique_ptr<WVFFTPlan> inner;
    auto status = reference_.createPlan(specification, inner);
    if (status) {
      ++counts_->plansCreated;
      plan = std::make_unique<CountingPlan>(std::move(inner), counts_);
    }
    return status;
  }

private:
  std::shared_ptr<ProviderCounts> counts_;
  WVReferenceFFTEngine reference_;
};

WVKernelStatus tryCreate(const WVFrozenForcingSchedule &forcingSchedule,
                         std::shared_ptr<ProviderCounts> counts) {
  std::unique_ptr<WVBarotropicQGForcingEngine> engine;
  return WVBarotropicQGForcingEngine::create(
      configuration(), forcingSchedule, catalog(),
      std::make_unique<CountingEngine>(std::move(counts)), engine);
}

void testPreflightAndLifecycle() {
  auto expectPreflightFailure = [](const char *label,
                                   WVFrozenForcingSchedule value) {
    auto counts = std::make_shared<ProviderCounts>();
    const auto status = tryCreate(value, counts);
    require(!status && counts->plansCreated == 0 &&
                counts->plansDestroyed == 0 && counts->executions == 0,
            std::string("invalid forcing was not rejected before numerical allocation: ") +
                label);
  };

  auto duplicate = schedule({nonlinear(0), beta(1)});
  duplicate.entries[1].name = duplicate.entries[0].name;
  expectPreflightFailure("duplicate name", std::move(duplicate));
  auto version = schedule({nonlinear()});
  ++version.entries[0].contractVersion;
  expectPreflightFailure("version", std::move(version));
  auto incompatible = schedule({linear()});
  incompatible.entries[0].stage = WVForcingStage::spectral;
  expectPreflightFailure("stage", std::move(incompatible));
  auto malformed = schedule({linear()});
  malformed.entries[0].configuration.values[0].storage =
      std::vector<double>{2e-7, 3e-7};
  malformed.entries[0].configuration.values[0].dimensions = {2};
  expectPreflightFailure("malformed scalar", std::move(malformed));
  auto duplicateIndex = schedule({entry(
      "WVFixedAmplitudeForcing", "fixed", WVForcingStage::spectralAmplitude,
      255, 0, fixedConfiguration({1, 1}, {{1.0, 0.0}, {2.0, 0.0}}))});
  expectPreflightFailure("duplicate fixed index", std::move(duplicateIndex));
  auto waveFamily = schedule({entry(
      "WVFixedAmplitudeForcing", "fixed", WVForcingStage::spectralAmplitude,
      255, 0, fixedConfiguration({1}, {{1.0, 0.0}}))});
  waveFamily.entries[0].configuration.values.push_back(
      integerValue("ApIndices", {1}));
  expectPreflightFailure("wave family", std::move(waveFamily));
  auto pseudo = schedule({entry("WVPseudoTopographicWaveGeneration", "pseudo",
                                WVForcingStage::spectral, 255, 0)});
  expectPreflightFailure("pseudo", std::move(pseudo));
  auto later = schedule({entry("WVVerticalDiffusivity", "later",
                               WVForcingStage::spectral, 255, 0)});
  expectPreflightFailure("later closure", std::move(later));

  auto counts = std::make_shared<ProviderCounts>();
  {
    auto engine = createEngine(configuration(), schedule({nonlinear()}),
                               std::make_unique<CountingEngine>(counts));
    require(counts->plansCreated == 3 && counts->plansDestroyed == 0,
            "provider did not retain exactly three balanced plans");
    const auto A0 = deterministicState(engine->kernel().descriptor());
    counts->failAtExecution = 1;
    std::vector<WVComplex64> output(A0.size(), {19.0, -7.0});
    WVComplexConstView input{A0.data(), engine->kernel().descriptor().spectralShape()};
    WVComplexView destination{output.data(), engine->kernel().descriptor().spectralShape()};
    auto status = engine->evaluateRightHandSide(input, destination);
    require(status.code == WVKernelStatusCode::fftExecutionFailure,
            "provider execution failure did not propagate");
    counts->failAtExecution = 0;
    status = engine->evaluateRightHandSide(input, destination);
    require(static_cast<bool>(status),
            "forcing engine did not clean up execution guard after failure");
  }
  require(counts->plansCreated == 3 && counts->plansDestroyed == 3,
          "provider plans were not balanced on destruction");
}

struct StateStorage {
  WVCoefficientStateStorage coefficients;
  WVMutableIntegrationState state;
  std::vector<WVCoefficientFamilyConstView> coefficientViews;
  std::vector<WVAdditionalStateBlockConstView> blockViews;

  explicit StateStorage(const WVIntegrationStateLayout &layout) {
    require(static_cast<bool>(coefficients.initialize(layout)),
            "compact state allocation failed");
    state.coefficientFamilies = coefficients.mutableFamilies();
    state.coefficientFamilyCount = coefficients.familyCount();
  }
  WVIntegrationState constView() {
    coefficientViews.clear();
    blockViews.clear();
    return integrationConstView(state, coefficientViews, blockViews);
  }
};

std::unique_ptr<WVBarotropicQGIntegrationSystem> createSystem(
    const WVFrozenForcingSchedule &forcingSchedule) {
  std::unique_ptr<WVBarotropicQGIntegrationSystem> result;
  const auto status = WVBarotropicQGIntegrationSystem::create(
      configuration(6, 5, 1, true), forcingSchedule, catalog(),
      std::make_unique<WVReferenceFFTEngine>(), result);
  require(static_cast<bool>(status), status.message);
  return result;
}

void initializeState(StateStorage &storage,
                     const WVTransformBarotropicQGDescriptor &descriptor) {
  const auto values = deterministicState(descriptor);
  std::copy(values.begin(), values.end(),
            storage.coefficients.mutableFamilies()[0].data);
  storage.state.waveVortex.t = 0.0;
  storage.state.waveVortex.t0 = 0.0;
}

void requireFixed(const StateStorage &state, std::size_t index,
                  WVComplex64 expected, const std::string &message) {
  const auto actual = state.coefficients.constFamilies()[0].data[index];
  require(actual.real == expected.real && actual.imag == expected.imag, message);
}

void testGenericIntegratorsConstraintsAndRestart() {
  const WVComplex64 fixedValue{3e-6, -2e-6};
  const auto forcingSchedule = schedule(
      {nonlinear(), linear(), damping(),
       entry("WVFixedAmplitudeForcing", "fixed",
             WVForcingStage::spectralAmplitude, 255, 3,
             fixedConfiguration({1}, {fixedValue}))});

  {
    auto system = createSystem(forcingSchedule);
    StateStorage state(system->stateLayout());
    initializeState(state, system->kernel().descriptor());
    WVFixedStepRK4 integrator(*system, {false});
    require(static_cast<bool>(integrator.prepareStateAfterRestart(state.state)),
            "RK4 restart preparation failed");
    require(static_cast<bool>(integrator.step(state.state, 1e-3)),
            "RK4 QG forcing step failed");
    requireFixed(state, 1, fixedValue, "RK4 fixed amplitude not restored");
  }
  {
    auto system = createSystem(forcingSchedule);
    StateStorage state(system->stateLayout());
    initializeState(state, system->kernel().descriptor());
    WVAdaptiveRK23Options options;
    options.relativeTolerance = 1e-8;
    options.absoluteToleranceScale = 1e-10;
    options.maximumStepSize = 1e-3;
    options.retainDenseOutput = false;
    WVAdaptiveRK23 integrator(*system, options);
    require(static_cast<bool>(integrator.prepareStateAfterRestart(state.state)),
            "RK23 restart preparation failed");
    require(static_cast<bool>(integrator.step(state.state, 1e-3)),
            "RK23 QG forcing step failed");
    requireFixed(state, 1, fixedValue, "RK23 fixed amplitude not restored");
  }
  {
    auto system = createSystem(forcingSchedule);
    StateStorage state(system->stateLayout());
    initializeState(state, system->kernel().descriptor());
    for (std::size_t index = 0;
         index < system->stateLayout().coefficientElementCount(); ++index) {
      state.coefficients.mutableFamilies()[0].data[index].real *= 1e5;
      state.coefficients.mutableFamilies()[0].data[index].imag *= 1e5;
    }
    state.coefficients.mutableFamilies()[0].data[1] = fixedValue;
    WVAdaptiveRK23Options options;
    options.relativeTolerance = 1e-12;
    options.absoluteToleranceScale = 1e-14;
    options.retainDenseOutput = false;
    WVAdaptiveRK23 integrator(*system, options);
    require(static_cast<bool>(integrator.prepareStateAfterRestart(state.state)),
            "rejection-path RK23 restart preparation failed");
    require(static_cast<bool>(integrator.step(state.state, 10.0)),
            "rejection-path RK23 QG forcing step failed");
    require(integrator.metrics().rejectedStepCount > 0,
            "real QG RK23 fixture did not exercise a rejected trial");
    requireFixed(state, 1, fixedValue,
                 "rejected RK23 trial changed fixed-amplitude state");
  }
  {
    auto system = createSystem(forcingSchedule);
    StateStorage state(system->stateLayout());
    initializeState(state, system->kernel().descriptor());
    WVAdaptiveRK45Options options;
    options.relativeTolerance = 1e-8;
    options.absoluteToleranceScale = 1e-10;
    options.maximumStepSize = 1e-3;
    options.retainDenseOutput = false;
    WVAdaptiveRK45 integrator(*system, options);
    require(static_cast<bool>(integrator.prepareStateAfterRestart(state.state)),
            "RK45 restart preparation failed");
    require(static_cast<bool>(integrator.step(state.state, 1e-3)),
            "RK45 QG forcing step failed");
    requireFixed(state, 1, fixedValue, "RK45 fixed amplitude not restored");
  }
  {
    auto system = createSystem(forcingSchedule);
    StateStorage state(system->stateLayout());
    initializeState(state, system->kernel().descriptor());
    WVAdaptiveRK78Options options;
    options.relativeTolerance = 1e-8;
    options.absoluteToleranceScale = 1e-10;
    options.maximumStepSize = 1e-3;
    options.retainDenseOutput = false;
    WVAdaptiveRK78 integrator(*system, options);
    require(static_cast<bool>(integrator.prepareStateAfterRestart(state.state)),
            "RK78 restart preparation failed");
    require(static_cast<bool>(integrator.step(state.state, 1e-3)),
            "endpoint-only RK78 QG forcing step failed");
    requireFixed(state, 1, fixedValue, "RK78 fixed amplitude not restored");
    require(integrator.lastAcceptedStep() != nullptr &&
                integrator.lastAcceptedStep()->denseOutput == nullptr &&
                integrator.metrics().denseHistoryCapacityBytes == 0 &&
                integrator.metrics().continuousExtensionRightHandSideEvaluationCount ==
                    0,
            "endpoint-only RK78 allocated or evaluated dense output");
  }

  auto uninterruptedSystem = createSystem(forcingSchedule);
  auto restartedSystem = createSystem(forcingSchedule);
  StateStorage uninterrupted(uninterruptedSystem->stateLayout());
  StateStorage restarted(restartedSystem->stateLayout());
  initializeState(uninterrupted, uninterruptedSystem->kernel().descriptor());
  initializeState(restarted, restartedSystem->kernel().descriptor());
  WVFixedStepRK4 uninterruptedIntegrator(*uninterruptedSystem, {false});
  WVFixedStepRK4 firstSegment(*restartedSystem, {false});
  require(static_cast<bool>(
              uninterruptedIntegrator.prepareStateAfterRestart(uninterrupted.state)) &&
              static_cast<bool>(
                  firstSegment.prepareStateAfterRestart(restarted.state)),
          "segmented restart preparation failed");
  require(static_cast<bool>(uninterruptedIntegrator.step(uninterrupted.state,
                                                        1e-3)) &&
              static_cast<bool>(uninterruptedIntegrator.step(
                  uninterrupted.state, 1e-3)) &&
              static_cast<bool>(firstSegment.step(restarted.state, 1e-3)),
          "first segmented continuation failed");
  WVTransformStateCheckpoint checkpoint;
  require(static_cast<bool>(captureTransformStateCheckpoint(
              restartedSystem->stateLayout(), restarted.constView(), checkpoint)),
          "forcing-system checkpoint capture failed");
  StateStorage restored(restartedSystem->stateLayout());
  require(static_cast<bool>(restoreTransformStateCheckpoint(
              checkpoint, restartedSystem->stateLayout(), restored.coefficients,
              restored.state)),
          "forcing-system checkpoint restore failed");
  WVFixedStepRK4 secondSegment(*restartedSystem, {false});
  require(static_cast<bool>(secondSegment.prepareStateAfterRestart(restored.state)) &&
              static_cast<bool>(secondSegment.step(restored.state, 1e-3)),
          "second segmented continuation failed");
  std::vector<WVComplex64> uninterruptedValues(
      uninterrupted.coefficients.mutableFamilies()[0].data,
      uninterrupted.coefficients.mutableFamilies()[0].data +
          uninterruptedSystem->stateLayout().coefficientElementCount());
  std::vector<WVComplex64> restoredValues(
      restored.coefficients.mutableFamilies()[0].data,
      restored.coefficients.mutableFamilies()[0].data +
          restartedSystem->stateLayout().coefficientElementCount());
  require(maximumRelativeError(restoredValues, uninterruptedValues) <= 1e-14,
          "checkpoint/restart segmented continuation diverged");
}

} // namespace

int main() {
  try {
    testNumericalForcingMatrix();
    testFixedAmplitudeAndNarrowBand();
    testPreflightAndLifecycle();
    testGenericIntegratorsConstraintsAndRestart();
    std::cout << "Barotropic QG forcing tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
