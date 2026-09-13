#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

double argument(char **values, int index) {
  return std::stod(values[index]);
}

void require(const WVKernelStatus &status) {
  if (!status)
    throw std::runtime_error(status.message);
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

WVPortableTypedRecord fixedConfiguration() {
  auto result = emptyConfiguration();
  result.values.push_back(integerValue("A0Indices", {1, 3}));
  result.values.push_back(realValue("A0ValuesReal", {3e-6, -4e-6}));
  result.values.push_back(realValue("A0ValuesImag", {-2e-6, 0.0}));
  return result;
}

WVFrozenForcingEntry entry(const char *identifier, const char *name,
                           WVForcingStage stage, std::uint8_t priority,
                           std::size_t ordinal,
                           WVPortableTypedRecord values =
                               emptyConfiguration()) {
  return {identifier, WVPortablePairContractVersion, name, stage, priority,
          ordinal, "matlab-parity", std::move(values)};
}

WVFrozenForcingEntry nonlinear(std::size_t ordinal = 0) {
  return entry("WVNonlinearAdvection", "nonlinear", WVForcingStage::spatial,
               127, ordinal);
}

WVFrozenForcingEntry damping(std::size_t ordinal = 0) {
  return entry("WVAdaptiveDamping", "damping", WVForcingStage::spectral, 255,
               ordinal);
}

WVFrozenForcingEntry linear(std::size_t ordinal = 0) {
  return entry("WVBottomFrictionLinear", "linear", WVForcingStage::spatial,
               220, ordinal, scalarConfiguration("r", 2.5e-7));
}

WVFrozenForcingEntry quadratic(std::size_t ordinal = 0) {
  return entry("WVBottomFrictionQuadratic", "quadratic",
               WVForcingStage::spatial, 210, ordinal,
               scalarConfiguration("Cd", 1.5e-3));
}

WVFrozenForcingEntry beta(std::size_t ordinal = 0) {
  return entry("WVBetaPlanePVAdvection", "beta", WVForcingStage::spatial,
               230, ordinal);
}

WVFrozenForcingEntry fixed(const char *identifier =
                               "WVFixedAmplitudeForcing",
                           std::size_t ordinal = 0) {
  return entry(identifier, "fixed", WVForcingStage::spectralAmplitude, 255,
               ordinal, fixedConfiguration());
}

WVFrozenForcingSchedule schedule(
    std::vector<WVFrozenForcingEntry> entries) {
  WVFrozenForcingSchedule result;
  result.entries = std::move(entries);
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

struct Evaluation {
  std::vector<WVComplex64> tendency;
  WVBarotropicQGForcingEngineMetrics metrics;
  std::string scheduleIdentifier;
};

Evaluation evaluate(const WVTransformBarotropicQGConfiguration &configuration,
                    const WVFrozenForcingSchedule &forcingSchedule,
                    const std::vector<WVComplex64> &A0,
                    const std::shared_ptr<const WVExtensionCatalog> &catalog) {
  std::unique_ptr<WVBarotropicQGForcingEngine> engine;
  require(WVBarotropicQGForcingEngine::create(
      configuration, forcingSchedule, catalog,
      std::make_unique<WVReferenceFFTEngine>(), engine));
  Evaluation result;
  result.tendency.resize(A0.size());
  WVComplexConstView input{A0.data(), engine->kernel().descriptor().spectralShape()};
  WVComplexView output{result.tendency.data(),
                       engine->kernel().descriptor().spectralShape()};
  require(engine->evaluateRightHandSide(input, output));
  result.metrics = engine->metrics();
  result.scheduleIdentifier = engine->scheduleIdentifier();
  return result;
}

void complexComponent(const std::vector<WVComplex64> &values,
                      bool imaginary) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0)
      std::cout << ',';
    std::cout << (imaginary ? values[index].imag : values[index].real);
  }
  std::cout << ']';
}

void namedComplex(const char *name, const Evaluation &value) {
  std::cout << ",\"" << name << "Real\":";
  complexComponent(value.tendency, false);
  std::cout << ",\"" << name << "Imag\":";
  complexComponent(value.tendency, true);
}

struct StateStorage {
  WVCoefficientStateStorage coefficients;
  WVMutableIntegrationState state;

  explicit StateStorage(const WVIntegrationStateLayout &layout) {
    require(coefficients.initialize(layout));
    state.coefficientFamilies = coefficients.mutableFamilies();
    state.coefficientFamilyCount = coefficients.familyCount();
  }
};

template <typename Integrator>
std::vector<WVComplex64> endpoint(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &forcingSchedule,
    const std::shared_ptr<const WVExtensionCatalog> &catalog,
    Integrator &&makeIntegrator) {
  std::unique_ptr<WVBarotropicQGIntegrationSystem> system;
  require(WVBarotropicQGIntegrationSystem::create(
      configuration, forcingSchedule, catalog,
      std::make_unique<WVReferenceFFTEngine>(), system));
  StateStorage storage(system->stateLayout());
  auto initial = deterministicState(system->kernel().descriptor());
  initial[1] = {3e-6, -2e-6};
  initial[3] = {-4e-6, 0.0};
  std::copy(initial.begin(), initial.end(),
            storage.coefficients.mutableFamilies()[0].data);
  storage.state.waveVortex.t = 0.0;
  storage.state.waveVortex.t0 = 0.0;
  auto integrator = makeIntegrator(*system);
  require(integrator->prepareStateAfterRestart(storage.state));
  require(integrator->advanceToTime(storage.state, 0.01, 0.005));
  const auto *family = storage.coefficients.constFamilies();
  return {family[0].data,
          family[0].data + system->stateLayout().coefficientElementCount()};
}

void namedComplex(const char *name,
                  const std::vector<WVComplex64> &values) {
  std::cout << ",\"" << name << "Real\":";
  complexComponent(values, false);
  std::cout << ",\"" << name << "Imag\":";
  complexComponent(values, true);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 12) {
    std::cerr << "usage: forcing-dump Nx Ny Lx Ly h j g rotationRate "
                 "latitude shouldAntialias planetaryRadius\n";
    return 2;
  }
  try {
    WVTransformBarotropicQGConfiguration configuration;
    configuration.Nx = static_cast<std::size_t>(argument(argv, 1));
    configuration.Ny = static_cast<std::size_t>(argument(argv, 2));
    configuration.Lx = argument(argv, 3);
    configuration.Ly = argument(argv, 4);
    configuration.h = argument(argv, 5);
    configuration.j = static_cast<std::uint32_t>(argument(argv, 6));
    configuration.g = argument(argv, 7);
    configuration.rotationRate = argument(argv, 8);
    configuration.latitude = argument(argv, 9);
    configuration.shouldAntialias = argument(argv, 10) != 0.0;
    configuration.planetaryRadius = argument(argv, 11);

    std::shared_ptr<const WVExtensionCatalog> extensions;
    require(makeBuiltInExtensionCatalog(extensions));
    std::unique_ptr<WVBarotropicQGForcingEngine> descriptorOwner;
    require(WVBarotropicQGForcingEngine::create(
        configuration, schedule({}), extensions,
        std::make_unique<WVReferenceFFTEngine>(), descriptorOwner));
    const auto A0 = deterministicState(descriptorOwner->kernel().descriptor());
    const auto nonlinearValue =
        evaluate(configuration, schedule({nonlinear()}), A0, extensions);
    const auto dampingValue =
        evaluate(configuration, schedule({damping()}), A0, extensions);
    const auto linearValue =
        evaluate(configuration, schedule({linear()}), A0, extensions);
    const auto quadraticValue =
        evaluate(configuration, schedule({quadratic()}), A0, extensions);
    const auto betaValue =
        evaluate(configuration, schedule({beta()}), A0, extensions);
    const auto nonlinearDamping = evaluate(
        configuration, schedule({damping(2), nonlinear(1)}), A0, extensions);
    const auto nonlinearLinear = evaluate(
        configuration, schedule({linear(2), nonlinear(1)}), A0, extensions);
    const auto nonlinearQuadratic = evaluate(
        configuration, schedule({quadratic(2), nonlinear(1)}), A0,
        extensions);
    const auto betaDamping = evaluate(
        configuration, schedule({damping(2), beta(1)}), A0, extensions);
    const auto all = evaluate(
        configuration,
        schedule({damping(5), beta(4), linear(3), nonlinear(1), quadratic(2)}),
        A0, extensions);
    const auto fixedValue = evaluate(
        configuration, schedule({fixed("WVFixedAmplitudeForcing", 2),
                                 nonlinear(1)}),
        A0, extensions);
    const auto narrowValue = evaluate(
        configuration,
        schedule({fixed("WVNarrowBandGeostrophicForcing", 2), nonlinear(1)}),
        A0, extensions);
    const auto integrationSchedule = schedule(
        {damping(5), beta(4), linear(3), nonlinear(1), quadratic(2),
         fixed("WVFixedAmplitudeForcing", 6)});
    const auto rk4Endpoint = endpoint(
        configuration, integrationSchedule, extensions,
        [](WVIntegrationSystem &system) {
          return std::make_unique<WVFixedStepRK4>(
              system, WVFixedStepRK4Options{false});
        });
    const auto rk23Endpoint = endpoint(
        configuration, integrationSchedule, extensions,
        [](WVIntegrationSystem &system) {
          WVAdaptiveRK23Options options;
          options.relativeTolerance = 1e-8;
          options.absoluteToleranceScale = 1e-10;
          options.retainDenseOutput = false;
          return std::make_unique<WVAdaptiveRK23>(system, options);
        });
    const auto rk45Endpoint = endpoint(
        configuration, integrationSchedule, extensions,
        [](WVIntegrationSystem &system) {
          WVAdaptiveRK45Options options;
          options.relativeTolerance = 1e-8;
          options.absoluteToleranceScale = 1e-10;
          options.retainDenseOutput = false;
          return std::make_unique<WVAdaptiveRK45>(system, options);
        });
    const auto rk78Endpoint = endpoint(
        configuration, integrationSchedule, extensions,
        [](WVIntegrationSystem &system) {
          WVAdaptiveRK78Options options;
          options.relativeTolerance = 1e-8;
          options.absoluteToleranceScale = 1e-10;
          options.retainDenseOutput = false;
          return std::make_unique<WVAdaptiveRK78>(system, options);
        });

    std::cout << std::setprecision(17) << "{\"Nkl\":" << A0.size();
    namedComplex("nonlinear", nonlinearValue);
    namedComplex("damping", dampingValue);
    namedComplex("linear", linearValue);
    namedComplex("quadratic", quadraticValue);
    namedComplex("beta", betaValue);
    namedComplex("nonlinearDamping", nonlinearDamping);
    namedComplex("nonlinearLinear", nonlinearLinear);
    namedComplex("nonlinearQuadratic", nonlinearQuadratic);
    namedComplex("betaDamping", betaDamping);
    namedComplex("all", all);
    namedComplex("fixed", fixedValue);
    namedComplex("narrow", narrowValue);
    namedComplex("rk4Endpoint", rk4Endpoint);
    namedComplex("rk23Endpoint", rk23Endpoint);
    namedComplex("rk45Endpoint", rk45Endpoint);
    namedComplex("rk78Endpoint", rk78Endpoint);
    const auto &metrics = all.metrics;
    std::cout << ",\"fieldReconstructionCount\":"
              << metrics.physicalFieldReconstructionCount
              << ",\"fieldReuseCount\":" << metrics.physicalFieldReuseCount
              << ",\"projectionCount\":"
              << metrics.spatialTendencyProjectionCount
              << ",\"forcingCallCount\":" << metrics.forcingCallCount
              << ",\"workspaceCapacityBytes\":"
              << metrics.workspaceCapacityBytes << "}\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 3;
  }
}
