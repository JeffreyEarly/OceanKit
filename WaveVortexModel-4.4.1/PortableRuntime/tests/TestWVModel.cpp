#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVAllocationProbe.hpp"

#include <algorithm>
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

void controlledCFLSelectionPreservesFixedSteps() {
  WVModel baseline, stopped;
  WVModelState baselineState, stoppedState;
  auto status = WVModel::createFromCheckpoint(
      test::extensionCatalog(), readFixture(),
      std::make_unique<WVReferenceFFTEngine>(), {}, baseline, baselineState);
  require(bool(status), status.message);
  status = WVModel::createFromCheckpoint(
      test::extensionCatalog(), readFixture(),
      std::make_unique<WVReferenceFFTEngine>(), {}, stopped, stoppedState);
  require(bool(status), status.message);
  require(bool(baseline.prepareStateAfterRestart(baselineState)) &&
              bool(stopped.prepareStateAfterRestart(stoppedState)),
          "CFL stop prepare");
  WVFixedTimeStepCandidates candidates;
  status =
      stopped.evaluateFixedTimeStepCandidates(stoppedState, 1e-7, candidates);
  require(bool(status), status.message);
  const double h = std::min(candidates.advective, candidates.oscillatory);
  require(std::isfinite(h) && h > 0.0,
          "CFL must select a finite positive step");
  const double initialTime = stoppedState.checkpoint().state.t;
  const double target = initialTime + 2.0 * h;
  status = baseline.advanceToTime(baselineState, target, h);
  require(bool(status), status.message);
  const auto result = stopped.advanceToTime(
      stoppedState, target, h, {[](const auto &progress) {
        return progress.boundary == WVIntegrationBoundary::acceptedStep;
      }});
  require(bool(result) && result.termination.stopped() &&
              result.metrics.integrator.acceptedStepCount == 1 &&
              result.metrics.integrator.rightHandSideEvaluationCount == 4 &&
              result.termination.finalAcceptedTime ==
                  stoppedState.checkpoint().state.t &&
              stoppedState.checkpoint().state.t == initialTime + h,
          "CFL stop must return the last accepted state and complete metrics");
  status = stopped.advanceToTime(stoppedState, target, h);
  require(bool(status) &&
              difference(baselineState.checkpoint(),
                         stoppedState.checkpoint()) == 0.0 &&
              baseline.metrics(&baselineState)
                      .integrator.rightHandSideEvaluationCount ==
                  stopped.metrics(&stoppedState)
                      .integrator.rightHandSideEvaluationCount,
          "CFL controlled stop and resume changed the fixed trajectory");
}

void variablePolicyChangeIsTransactional() {
  auto checkpoint=readFixture();
  const auto nonlinear=std::find_if(checkpoint.forcingSchedule.entries.begin(),
      checkpoint.forcingSchedule.entries.end(),[](const auto& entry) {
        return entry.typeIdentifier=="WVNonlinearAdvection";
      });
  require(nonlinear!=checkpoint.forcingSchedule.entries.end(),
      "Transactional policy fixture lacks nonlinear forcing");
  auto repeated=*nonlinear;
  repeated.name="second nonlinear advection";
  repeated.ordinal=0;
  for(const auto& entry:checkpoint.forcingSchedule.entries)
    repeated.ordinal=std::max(repeated.ordinal,entry.ordinal+1);
  checkpoint.forcingSchedule.entries.push_back(std::move(repeated));

  WVPortableObserverRecord record;
  const auto shape=checkpoint.state.coefficients.shape;
  for(const char* identifier:{"Ap","Am","A0"})
    record.stateBlocks.push_back({identifier,WVStateScalarType::complex64,
        {shape.rows,shape.columns},WVToleranceKind::coefficientEnergyScaled,
        1e-10,WVStateOwnership::integratorOwned,
        WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier="coefficients";
  coefficients.name="Wave-vortex coefficients";
  coefficients.typeIdentifier="WVCoefficients";
  coefficients.stateBlockIdentifiers={"Ap","Am","A0"};
  record.observers.push_back(std::move(coefficients));
  WVPortableObserverDescriptor descriptor;
  auto status=WVPortableObserverDescriptor::create(
      record,test::extensionCatalog(),descriptor);
  require(bool(status),status.message);

  WVModel model;
  status=WVModel::create(test::extensionCatalog(),checkpoint.configuration,
      checkpoint.forcingSchedule,descriptor,
      std::make_unique<WVReferenceFFTEngine>(),{},model);
  require(bool(status),status.message);
  require(bool(model.setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy::lowMemory)),
      "Transactional model low-memory setup failed");
  const auto lowMetrics=model.metrics();

  // The forcing cache is the first allocation. Fail the following field-arena
  // preparation and require the model-level rollback to restore both policies.
  allocationProbe::calls=0;
  allocationProbe::counting=true;
  allocationProbe::failAfter=1;
  status=model.setVariableEvaluationPolicy(WVVariableEvaluationPolicy::reuse);
  allocationProbe::failAfter=-1;
  allocationProbe::counting=false;
  const auto failedAllocationCalls=allocationProbe::calls.load();
  const auto failedMetrics=model.metrics();
  require(status.code==WVKernelStatusCode::allocationFailure &&
      failedAllocationCalls>=2 &&
      failedMetrics.forcing.workspaceCapacityBytes==
          lowMetrics.forcing.workspaceCapacityBytes &&
      failedMetrics.integrationSystemPersistentBytes==
          lowMetrics.integrationSystemPersistentBytes,
      "Failed field preparation left model policy participants mismatched");
  require(bool(model.setVariableEvaluationPolicy(
              WVVariableEvaluationPolicy::reuse)) &&
          bool(model.setVariableEvaluationPolicy(
              WVVariableEvaluationPolicy::lowMemory)) &&
          model.metrics().integrationSystemPersistentBytes==
              lowMetrics.integrationSystemPersistentBytes,
      "Transactional model policy retry changed low-memory storage");
}

} // namespace

int main() {
  try {
    variablePolicyChangeIsTransactional();
    controlledCFLSelectionPreservesFixedSteps();
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
