#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

double difference(const WVCheckpoint &left, const WVCheckpoint &right) {
  double maximum = 0.0;
  const std::vector<WVComplex64> *a[] = {
      &left.state.coefficients.Ap, &left.state.coefficients.Am,
      &left.state.coefficients.A0};
  const std::vector<WVComplex64> *b[] = {
      &right.state.coefficients.Ap, &right.state.coefficients.Am,
      &right.state.coefficients.A0};
  for (std::size_t component = 0; component < 3; ++component)
    for (std::size_t index = 0; index < a[component]->size(); ++index)
      maximum = std::max(maximum,
                         std::hypot((*a[component])[index].real -
                                        (*b[component])[index].real,
                                    (*a[component])[index].imag -
                                        (*b[component])[index].imag));
  return maximum;
}

WVCheckpoint readFixture() {
  WVCheckpoint checkpoint;
  const auto status = WVCheckpointReader::read(
      std::string(WV_RUNTIME_FIXTURE_DIR) + "/forcing-nonlinear.nc",
      *test::extensionCatalog(), checkpoint);
  require(static_cast<bool>(status), status.message);
  return checkpoint;
}

void fixedFacadeMatchesDirectIntegrator() {
  auto directCheckpoint = readFixture();
  auto facadeCheckpoint = directCheckpoint;

  std::unique_ptr<WVConstantStratificationIntegrationSystem> directSystem;
  auto status = WVConstantStratificationIntegrationSystem::create(
      directCheckpoint.configuration, directCheckpoint.forcingSchedule,
      test::extensionCatalog(),
      std::make_unique<WVReferenceFFTEngine>(), directSystem);
  require(static_cast<bool>(status), status.message);
  WVAdditionalStateStorage directAdditional;
  status = directAdditional.initialize(directSystem->stateLayout());
  require(static_cast<bool>(status), status.message);
  const auto shape = directCheckpoint.state.coefficients.shape;
  WVMutableIntegrationState directState{
      {directCheckpoint.state.t,
       directCheckpoint.state.t0,
       {{directCheckpoint.state.coefficients.Ap.data(), shape},
        {directCheckpoint.state.coefficients.Am.data(), shape},
        {directCheckpoint.state.coefficients.A0.data(), shape}}},
      directAdditional.mutableBlocks(), directAdditional.blockCount()};
  WVFixedStepRK4 direct(*directSystem);
  status = direct.prepareStateAfterRestart(directState);
  require(static_cast<bool>(status), status.message);

  WVModel model;
  status = WVModel::create(
      test::extensionCatalog(), facadeCheckpoint.configuration,
      facadeCheckpoint.forcingSchedule,
      std::make_unique<WVReferenceFFTEngine>(), {}, model);
  require(static_cast<bool>(status), status.message);
  WVModelState state;
  status = WVModelState::create(std::move(facadeCheckpoint),
                                model.stateLayout(), state);
  require(static_cast<bool>(status), status.message);
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);

  constexpr double step = 1e-5;
  status = direct.step(directState, step);
  require(static_cast<bool>(status), status.message);
  directCheckpoint.state.t = directState.waveVortex.t;
  status = model.step(state, step);
  require(static_cast<bool>(status), status.message);
  require(difference(directCheckpoint, state.checkpoint()) == 0.0,
          "WVModel changed fixed-RK4 results");
  require(model.metrics(&state).statePersistentBytes ==
              state.persistentBytes(),
          "WVModel state accounting mismatch");
}

void adaptiveFacadeAdvances() {
  auto checkpoint = readFixture();
  WVModelIntegratorConfiguration options;
  options.kind = WVModelIntegratorKind::adaptiveRK23;
  options.adaptive.maximumStepSize = 1e-5;
  WVModel model;
  auto status = WVModel::create(
      test::extensionCatalog(), checkpoint.configuration,
      checkpoint.forcingSchedule,
      std::make_unique<WVReferenceFFTEngine>(), options, model);
  require(static_cast<bool>(status), status.message);
  WVModelState state;
  status = WVModelState::create(std::move(checkpoint), model.stateLayout(),
                                state);
  require(static_cast<bool>(status), status.message);
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);
  const auto target = state.checkpoint().state.t + 2e-5;
  status = model.advanceToTime(state, target, 1e-5);
  require(static_cast<bool>(status), status.message);
  require(std::abs(state.checkpoint().state.t - target) <= 1e-14,
          "WVModel adaptive integration stopped at the wrong time");
  require(model.integratorKind() == WVModelIntegratorKind::adaptiveRK23,
          "WVModel lost the active integrator identity");
}

void adaptiveRK78FacadeAdvances() {
  auto checkpoint = readFixture();
  WVModelIntegratorConfiguration options;
  options.kind = WVModelIntegratorKind::adaptiveRK78;
  options.adaptiveRK78.maximumStepSize = 1e-5;
  WVModel model;
  auto status = WVModel::create(
      test::extensionCatalog(), checkpoint.configuration,
      checkpoint.forcingSchedule,
      std::make_unique<WVReferenceFFTEngine>(), options, model);
  require(static_cast<bool>(status), status.message);
  WVModelState state;
  status = WVModelState::create(std::move(checkpoint), model.stateLayout(),
                                state);
  require(static_cast<bool>(status), status.message);
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);
  const auto target = state.checkpoint().state.t + 2e-5;
  status = model.advanceToTime(state, target, 1e-5);
  require(static_cast<bool>(status), status.message);
  const auto metrics = model.metrics(&state);
  require(std::abs(state.checkpoint().state.t - target) <= 1e-14 &&
              model.integratorKind() == WVModelIntegratorKind::adaptiveRK78 &&
              metrics.integrator.workspaceStateEquivalentCount == 11 &&
              metrics.integrator.denseHistoryCapacityBytes == 0,
          "WVModel lost the endpoint-only adaptive-rk78 contract");
}

} // namespace

int main() {
  try {
    fixedFacadeMatchesDirectIntegrator();
    adaptiveFacadeAdvances();
    adaptiveRK78FacadeAdvances();
    std::cout << "WVModel façade tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
