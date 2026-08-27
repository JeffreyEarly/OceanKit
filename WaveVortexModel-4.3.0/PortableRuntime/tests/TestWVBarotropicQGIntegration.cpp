#include "WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
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

void requireRelative(double actual, double expected, double tolerance,
                     const std::string &message) {
  const auto scale = std::max(std::abs(expected), 1.0);
  require(std::abs(actual - expected) <= tolerance * scale, message);
}

WVBarotropicQGPersistedNumericalRecord persistedRecord() {
  WVBarotropicQGPersistedNumericalRecord record;
  record.Lx = 17000.0;
  record.Ly = 11000.0;
  record.x.resize(9);
  record.y.resize(6);
  for (std::size_t index = 0; index < record.x.size(); ++index)
    record.x[index] = record.Lx * static_cast<double>(index) /
                      static_cast<double>(record.x.size());
  for (std::size_t index = 0; index < record.y.size(); ++index)
    record.y[index] = record.Ly * static_cast<double>(index) /
                      static_cast<double>(record.y.size());
  record.h = 0.8;
  record.j = 1;
  record.g = 9.81;
  record.planetaryRadius = 6.371e6;
  record.rotationRate = 7.2921e-5;
  record.latitude = 33.0;
  record.shouldAntialias = true;
  record.t = 17.0;
  record.t0 = -3.0;
  return record;
}

struct StateStorage {
  WVCoefficientStateStorage coefficients;
  WVMutableIntegrationState state;
  std::vector<WVCoefficientFamilyConstView> coefficientViews;
  std::vector<WVAdditionalStateBlockConstView> blockViews;

  explicit StateStorage(const WVIntegrationStateLayout &layout) {
    require(static_cast<bool>(coefficients.initialize(layout)),
            "compact coefficient allocation");
    state.waveVortex.t = 0.0;
    state.waveVortex.t0 = 0.0;
    state.coefficientFamilies = coefficients.mutableFamilies();
    state.coefficientFamilyCount = coefficients.familyCount();
  }

  WVIntegrationState constView() {
    coefficientViews.clear();
    blockViews.clear();
    return integrationConstView(state, coefficientViews, blockViews);
  }
};

void initializeState(StateStorage &storage) {
  auto &family = storage.coefficients.mutableFamilies()[0];
  for (std::size_t index = 0; index < family.layout->elementCount; ++index)
    family.data[index] = {
        2e-5 * std::sin(0.31 * static_cast<double>(index + 1)),
        1e-5 * std::cos(0.17 * static_cast<double>(index + 3))};
  family.data[0].imag = 0.0;
}

void testConfigurationDecode() {
  auto record = persistedRecord();
  WVBarotropicQGNumericalConfiguration decoded;
  auto status = decodeBarotropicQGNumericalConfiguration(record, decoded);
  require(static_cast<bool>(status), "persisted numerical decode");
  require(decoded.transform.Nx == record.x.size() &&
              decoded.transform.Ny == record.y.size() &&
              decoded.transform.h == record.h &&
              decoded.transform.j == record.j &&
              decoded.transform.g == record.g &&
              decoded.transform.planetaryRadius == record.planetaryRadius &&
              decoded.transform.rotationRate == record.rotationRate &&
              decoded.transform.latitude == record.latitude &&
              decoded.transform.shouldAntialias == record.shouldAntialias &&
              decoded.t == record.t && decoded.t0 == record.t0,
          "complete persisted physical/model parameter decode");
  require(decoded.stateDescription.transformIdentifier ==
                  "WVTransformBarotropicQG" &&
              decoded.stateDescription.spatialDimensions ==
                  std::vector<std::size_t>({9, 6}) &&
              decoded.stateDescription.coefficientFamilies.size() == 1 &&
              decoded.stateDescription.coefficientFamilies[0].identifier ==
                  "A0" &&
              decoded.stateDescription.coefficientFamilies[0]
                      .spectralDimensions.size() == 1,
          "allocation-light compact state description");
  require(decoded.persistentBytes() >= sizeof(decoded),
          "decoded storage accounting");

  auto invalidCoordinate = record;
  invalidCoordinate.x[3] += 1.0;
  status =
      decodeBarotropicQGNumericalConfiguration(invalidCoordinate, decoded);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "invalid persisted coordinate rejected");
  auto invalidJ = record;
  invalidJ.j = 2;
  status = decodeBarotropicQGNumericalConfiguration(invalidJ, decoded);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "invalid persisted j rejected");
}

void testSystemAndIntegrators() {
  WVBarotropicQGNumericalConfiguration decoded;
  require(static_cast<bool>(decodeBarotropicQGNumericalConfiguration(
              persistedRecord(), decoded)),
          "decode before system creation");
  std::unique_ptr<WVBarotropicQGIntegrationSystem> system;
  auto status = WVBarotropicQGIntegrationSystem::create(
      decoded.transform, std::make_unique<WVReferenceFFTEngine>(), system);
  require(static_cast<bool>(status) && system,
          "Barotropic QG integration-system creation");
  const auto &layout = system->stateLayout();
  require(layout.transformIdentifier() == "WVTransformBarotropicQG" &&
              layout.spatialDimensions() ==
                  std::vector<std::size_t>({9, 6}) &&
              layout.coefficientFamilyCount() == 1 &&
              layout.coefficientFamilies()[0].identifier == "A0" &&
              layout.coefficientFamilies()[0].spectralDimensions ==
                  std::vector<std::size_t>(
                      {system->kernel().descriptor().Nkl()}) &&
              !layout.hasLegacyCoefficientTriple(),
          "real system declares compact A0-only layout");

  StateStorage state(layout);
  StateStorage flux(layout);
  const auto expectedStateCapacity =
      layout.coefficientElementCount() * sizeof(WVComplex64) +
      sizeof(WVCoefficientFamilyView) +
      sizeof(WVCoefficientFamilyConstView);
  require(state.coefficients.capacityBytes() == expectedStateCapacity &&
              layout.coefficientElementCount() ==
                  system->kernel().descriptor().Nkl(),
          "exact compact A0-only state byte accounting");
  initializeState(state);
  WVIntegrationFlux rightHandSide;
  rightHandSide.coefficientFamilies = flux.coefficients.mutableFamilies();
  rightHandSide.coefficientFamilyCount = flux.coefficients.familyCount();
  status = system->evaluateRightHandSide(state.constView(), rightHandSide);
  require(static_cast<bool>(status), "real compact nonlinear RHS");
  require(state.state.waveVortex.coefficients.Ap.data == nullptr &&
              state.state.waveVortex.coefficients.Am.data == nullptr &&
              state.state.waveVortex.coefficients.A0.data == nullptr &&
              state.coefficients.familyCount() == 1 &&
              flux.coefficients.familyCount() == 1,
          "no legacy or dummy wave coefficient storage");

  WVTransformStateCheckpoint checkpoint;
  status = captureTransformStateCheckpoint(layout, state.constView(),
                                            checkpoint);
  require(static_cast<bool>(status) &&
              checkpoint.transformIdentifier ==
                  "WVTransformBarotropicQG" &&
              checkpoint.coefficientFamilies.size() == 1 &&
              checkpoint.coefficientFamilies[0].identifier == "A0" &&
              checkpoint.coefficientFamilies[0].spectralDimensions ==
                  std::vector<std::size_t>(
                      {system->kernel().descriptor().Nkl()}),
          "transform-neutral checkpoint captures only A0");
  StateStorage restored(layout);
  status = restoreTransformStateCheckpoint(checkpoint, layout,
                                            restored.coefficients,
                                            restored.state);
  require(static_cast<bool>(status) &&
              restored.coefficients.familyCount() == 1 &&
              restored.state.waveVortex.coefficients.Ap.data == nullptr &&
              restored.state.waveVortex.coefficients.Am.data == nullptr,
          "transform-neutral checkpoint restores only A0");
  for (std::size_t index = 0;
       index < layout.coefficientElementCount(); ++index)
    require(restored.coefficients.mutableFamilies()[0].data[index].real ==
                    state.coefficients.mutableFamilies()[0].data[index].real &&
                restored.coefficients.mutableFamilies()[0].data[index].imag ==
                    state.coefficients.mutableFamilies()[0].data[index].imag,
            "A0 checkpoint value round trip");

  std::unique_ptr<WVIntegrationErrorPolicy> policy;
  status = system->createErrorPolicy(1e-8, policy);
  require(static_cast<bool>(status) && policy &&
              policy->componentCount() == 1 &&
              policy->elementCount(0) == layout.coefficientElementCount() &&
              policy->absoluteTolerance(0, 0) == 1.0 &&
              std::isfinite(policy->absoluteTolerance(0, 1)),
          "MATLAB-compatible one-family adaptive error policy");

  state.coefficients.mutableFamilies()[0].data[0].imag = 1.0;
  const auto constrained = system->enforceStateConstraints(state.state);
  require(static_cast<bool>(constrained) &&
              constrained.modifiedCoefficientCount == 1 &&
              !constrained.fsalCompatible &&
              state.coefficients.mutableFamilies()[0].data[0].real == 0.0 &&
              state.coefficients.mutableFamilies()[0].data[0].imag == 0.0,
          "masked zero mode and self-conjugate reality constraint");

  WVFixedStepRK4 rk4(*system, {true});
  status = rk4.prepareStateAfterRestart(state.state);
  require(static_cast<bool>(status), "generic RK4 restart preparation");
  status = rk4.step(state.state, 0.01);
  require(static_cast<bool>(status) &&
              rk4.metrics().workspaceCapacityBytes ==
                  4 * layout.coefficientElementCount() *
                      sizeof(WVComplex64),
          "generic RK4 advances one compact family");
  StateStorage dense(layout);
  status = rk4.evaluateDenseOutput(0.005, dense.state);
  require(static_cast<bool>(status) &&
              dense.state.waveVortex.coefficients.Ap.data == nullptr &&
              dense.state.coefficientFamilyCount == 1,
          "generic dense output preserves A0-only state");

  WVAdaptiveRK23Options rk23Options;
  rk23Options.relativeTolerance = 1e-6;
  rk23Options.maximumStepSize = 0.001;
  WVAdaptiveRK23 rk23(*system, rk23Options);
  status = rk23.prepareStateAfterRestart(state.state);
  require(static_cast<bool>(status), "generic RK23 restart preparation");
  status = rk23.step(state.state, 0.001);
  require(static_cast<bool>(status) &&
              rk23.metrics().workspaceStateEquivalentCount == 5,
          "generic RK23 advances one compact family");
  status = rk23.evaluateDenseOutput(state.state.waveVortex.t - 0.0005,
                                    dense.state);
  require(static_cast<bool>(status) &&
              dense.state.coefficientFamilyCount == 1,
          "generic RK23 dense output preserves A0-only state");

  WVAdaptiveRK45Options rk45Options;
  rk45Options.relativeTolerance = 1e-6;
  rk45Options.maximumStepSize = 0.001;
  WVAdaptiveRK45 rk45(*system, rk45Options);
  status = rk45.prepareStateAfterRestart(state.state);
  require(static_cast<bool>(status), "generic RK45 restart preparation");
  status = rk45.step(state.state, 0.001);
  require(static_cast<bool>(status) &&
              rk45.metrics().workspaceStateEquivalentCount == 7 &&
              state.state.waveVortex.coefficients.Ap.data == nullptr &&
              state.state.waveVortex.coefficients.Am.data == nullptr &&
              state.state.waveVortex.coefficients.A0.data == nullptr,
          "generic RK45 has no Ap/Am or transform dispatch assumption");
  status = rk45.evaluateDenseOutput(state.state.waveVortex.t - 0.0005,
                                    dense.state);
  require(static_cast<bool>(status) &&
              dense.state.coefficientFamilyCount == 1,
          "generic RK45 dense output preserves A0-only state");

  WVAdaptiveRK78Options rk78Options;
  rk78Options.relativeTolerance = 1e-6;
  rk78Options.maximumStepSize = 0.001;
  WVAdaptiveRK78 rk78(*system, rk78Options);
  status = rk78.prepareStateAfterRestart(state.state);
  require(static_cast<bool>(status), "real-system RK78 restart preparation");
  status = rk78.step(state.state, 0.001);
  const auto *accepted = rk78.lastAcceptedStep();
  require(static_cast<bool>(status) && accepted != nullptr &&
              accepted->denseOutput == nullptr &&
              rk78.metrics().workspaceStateEquivalentCount == 11 &&
              rk78.metrics().denseHistoryStateEquivalentCount == 0 &&
              rk78.metrics().denseHistoryCapacityBytes == 0 &&
              rk78.metrics().acceptedStepCount == 1 &&
              state.state.coefficientFamilyCount == 1 &&
              state.state.waveVortex.coefficients.Ap.data == nullptr &&
              state.state.waveVortex.coefficients.Am.data == nullptr &&
              state.state.waveVortex.coefficients.A0.data == nullptr &&
              std::string(WVAdaptiveRK78::methodIdentifier()) ==
                  "adaptive-rk78" &&
              std::string(WVAdaptiveRK78::controllerIdentifier()) ==
                  "matlab-ode78-v1",
          "endpoint-only RK78 advances the real compact QG system");
  require(system->persistentBytes() >= system->kernel().persistentBytes() &&
              system->kernel().metrics().persistentFullHermitianBytes == 0,
          "system retained-storage and compact-spectrum evidence");
}

void testMatlabCFLFixture() {
  WVTransformBarotropicQGConfiguration configuration;
  configuration.Nx = 6;
  configuration.Ny = 5;
  configuration.Lx = 15000.0;
  configuration.Ly = 12000.0;
  configuration.h = 0.8;
  configuration.j = 1;
  configuration.g = 9.80665;
  configuration.planetaryRadius = 6.3712e6;
  configuration.rotationRate = 7.292115e-5;
  configuration.latitude = 33.0;
  configuration.shouldAntialias = true;
  std::unique_ptr<WVBarotropicQGIntegrationSystem> system;
  auto status = WVBarotropicQGIntegrationSystem::create(
      configuration, std::make_unique<WVReferenceFFTEngine>(), system);
  require(static_cast<bool>(status), "QG CFL system creation");
  StateStorage state(system->stateLayout());
  auto &family = state.coefficients.mutableFamilies()[0];
  for (std::size_t index = 0; index < family.layout->elementCount; ++index) {
    const double p = static_cast<double>(index + 1);
    family.data[index] = {1e-5 * std::sin(0.19 * (p + 4.0)),
                          1e-5 * std::cos(0.05 * (p + 5.0))};
  }
  const auto persistentBefore = system->persistentBytes();
  WVFixedTimeStepCandidates candidates;
  status = system->evaluateFixedTimeStepCandidates(state.constView(), 0.4,
                                                    candidates);
  require(static_cast<bool>(status) &&
              system->supportsFixedTimeStepSelection(),
          "QG CFL evaluation");
  requireRelative(candidates.effectiveHorizontalGridResolution,
                  6000.0000000000009, 1e-13,
                  "QG MATLAB effective horizontal resolution fixture");
  requireRelative(candidates.maximumHorizontalSpeed,
                  0.11956917129957827, 1e-12,
                  "QG MATLAB uvMax fixture");
  requireRelative(candidates.horizontalAdvective,
                  20072.063508635067, 1e-12,
                  "QG MATLAB advective CFL fixture");
  require(candidates.advective == candidates.horizontalAdvective &&
              std::isinf(candidates.verticalAdvective) &&
              std::isinf(candidates.oscillatory) &&
              candidates.highestActiveWaveFrequency == 0.0 &&
              candidates.transientWorkspaceMaximumLiveBytes == 0 &&
              system->persistentBytes() == persistentBefore,
          "QG has no vertical or wave limit and retains no CFL workspace");

  std::fill_n(family.data, family.layout->elementCount, WVComplex64{});
  status = system->evaluateFixedTimeStepCandidates(state.constView(), 0.4,
                                                    candidates);
  require(static_cast<bool>(status) &&
              candidates.maximumHorizontalSpeed == 0.0 &&
              std::isinf(candidates.advective) &&
              std::isinf(candidates.oscillatory),
          "zero-velocity QG candidates are infinite");
}

} // namespace

int main() {
  testConfigurationDecode();
  testSystemAndIntegrators();
  testMatlabCFLFixture();
  std::cout << "Barotropic QG integration tests passed\n";
  return 0;
}
