#include "WVDensityEventEvaluation.hpp"
#include "WVAllocationProbe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::detail;

namespace {
void require(bool value, const char *message) {
  if (!value) throw std::runtime_error(message);
}
void close(double value, double expected, double tolerance, const char *message) {
  require(std::isfinite(value) && std::abs(value - expected) <= tolerance, message);
}
void untouched(const std::vector<double> &values) {
  require(std::all_of(values.begin(), values.end(), [](double x) { return x == -719; }),
          "failed preparation changed caller output");
}
struct Fixture {
  std::vector<double> z{-2, -1, 0};
  std::vector<double> weights{.5, 1, .5};
  std::vector<double> initial{1027, 1026, 1025};
  std::vector<double> density{1027, 1027, 1026.5, 1026.5, 1025, 1025};
  double gravity = 10;
  WVDensityEventGeometry geometry() const { return {&z, &weights, &initial, 2, gravity, 1000}; }
  WVRealVolumeConstView volume() const { return {density.data(), {2, 1, 3}}; }
  void bind(WVDensityEventEvaluation &event, WVNoMotionReference reference,
            const WVNoMotionRecoveryOptions &options = {}) const {
    WVDensityDiagnosticContract contract;
    contract.reference = reference;
    require(bool(WVDensityEventEvaluation::create(volume(), geometry(), contract, event, options)),
            "valid event binding failed");
  }
};
std::array<WVDensityEventOutput, 3> outputs(std::vector<double> &rho,
                                           std::vector<double> &eta,
                                           std::vector<double> &ape) {
  return {{{WVDensityEventField::rhoNm, rho.data(), rho.size()},
           {WVDensityEventField::etaTrue, eta.data(), eta.size()},
           {WVDensityEventField::ape, ape.data(), ape.size()}}};
}

void mixedReferencesAndDemandOrder() {
  for (const auto reference : {WVNoMotionReference::actual, WVNoMotionReference::initial}) {
    Fixture fixture;
    WVDensityEventEvaluation event;
    fixture.bind(event, reference);
    require(event.metrics().liveBytes == 0 && event.metrics().recoveryAttemptCount == 0,
            "binding eagerly allocated or recovered density");
    std::vector<double> rho(3, -719), eta(6, -719), ape(6, -719);
    auto requests = outputs(rho, eta, ape);
    std::reverse(requests.begin(), requests.end());
    require(bool(event.evaluate(requests.data(), requests.size())), "mixed density request failed");
    require(rho == std::vector<double>({1027, 1026.5, 1025}),
            "rho_nm incorrectly followed the selected initial reference");
    for (std::size_t i = 0; i < eta.size(); ++i) {
      const bool middle = i == 2 || i == 3;
      const double expectedEta = reference == WVNoMotionReference::initial && middle ? .5 : 0;
      close(eta[i], expectedEta, 2e-15, "mixed-reference displacement is incorrect");
      close(ape[i], .5 * .01 * expectedEta * expectedEta, 2e-17,
            "mixed-reference APE is incorrect");
    }
    const auto &m = event.metrics();
    require(m.recoveryCount == 1 && m.profileConstructionCount == 1 &&
                m.inversePassCount == 1 && m.inverseSampleCount == 6 &&
                m.apePassCount == 1 && m.apeSampleCount == 6,
            "combined diagnostics repeated expensive work");
    require(m.liveBytes > 0 && m.highWaterBytes >= m.liveBytes &&
                m.profileConstructionWorkspaceHighWaterBytes > 0,
            "density storage was not accounted");
  }

  Fixture fixture;
  WVDensityEventEvaluation event;
  fixture.bind(event, WVNoMotionReference::initial);
  std::vector<double> rho(3), eta(6), ape(6);
  auto requests = outputs(rho, eta, ape);
  require(bool(event.evaluate(&requests[1], 1)), "initial-only eta failed");
  require(event.metrics().recoveryCount == 0 && event.view(WVDensityEventField::rhoNm).data == nullptr,
          "initial-only eta unnecessarily recovered actual density");
  const auto *material = event.materialHeights().data;
  const auto *etaStorage = event.view(WVDensityEventField::etaTrue).data;
  require(bool(event.evaluate(&requests[2], 1)), "sequential APE demand failed");
  require(event.metrics().inversePassCount == 1 && event.metrics().recoveryCount == 0,
          "APE repeated inversion or recovered an unused actual profile");
  require(bool(event.evaluate(&requests[0], 1)), "sequential actual profile demand failed");
  require(event.metrics().profileConstructionCount == 1 && event.metrics().recoveryCount == 1 &&
              event.materialHeights().data == material &&
              event.view(WVDensityEventField::etaTrue).data == etaStorage,
          "demand extension invalidated prepared views or rebuilt the selected profile");
  const auto bytes = event.metrics().liveBytes;
  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  bool passed = true;
  for (int i = 0; i < 1000; ++i)
    passed = bool(event.evaluate(requests.data(), requests.size())) && passed;
  allocationProbe::counting = false;
  require(passed && allocationProbe::calls == 0 && event.metrics().liveBytes == bytes &&
              event.metrics().recoveryCount == 1 && event.metrics().inversePassCount == 1 &&
              event.metrics().apePassCount == 1,
          "prepared replay allocated or repeated numerical work");

  WVDensityEventEvaluation profileOnly;
  fixture.bind(profileOnly, WVNoMotionReference::actual);
  require(bool(profileOnly.evaluate(&requests[0], 1)), "rho_nm-only request failed");
  require(profileOnly.metrics().recoveryCount == 1 &&
              profileOnly.metrics().profileConstructionCount == 0 &&
              profileOnly.metrics().inversePassCount == 0 && !profileOnly.materialHeights().data,
          "rho_nm-only request allocated parcel calculus");
}

void preflightAndTransactionalFailures() {
  Fixture fixture;
  WVDensityEventEvaluation event;
  fixture.bind(event, WVNoMotionReference::initial);
  std::vector<double> rho(3, -719), eta(6, -719), ape(6, -719);
  auto requests = outputs(rho, eta, ape);
  requests.back().elementCount = 5;
  auto status = event.evaluate(requests.data(), requests.size());
  require(status.code == WVKernelStatusCode::invalidShape && event.metrics().recoveryCount == 0,
          "late invalid output was not rejected before recovery");
  untouched(rho); untouched(eta); untouched(ape);
  requests = outputs(rho, eta, ape);
  requests.back().data = eta.data();
  require(event.evaluate(requests.data(), requests.size()).code == WVKernelStatusCode::overlappingArrays,
          "overlapping diagnostic outputs accepted");
  untouched(rho); untouched(eta); untouched(ape);
  WVDensityEventOutput inputAlias{WVDensityEventField::rhoNm, fixture.initial.data(), 3};
  require(event.evaluate(&inputAlias, 1).code == WVKernelStatusCode::overlappingArrays,
          "output alias of borrowed profile accepted");
  require(fixture.initial == std::vector<double>({1027, 1026, 1025}), "input profile was modified");
  require(bool(event.prepare(WVDensityEventEvaluation::rhoNmDemand)), "profile preparation failed");
  auto cached = event.view(WVDensityEventField::rhoNm);
  WVDensityEventOutput cacheAlias{WVDensityEventField::rhoNm, const_cast<double *>(cached.data), cached.elementCount};
  require(event.evaluate(&cacheAlias, 1).code == WVKernelStatusCode::overlappingArrays,
          "output alias of prepared profile accepted");
  require(!event.prepare(128), "unknown demand accepted");

  // Recovery succeeds, then initial-reference inversion rejects changed extrema.
  Fixture outside;
  outside.density[0] = outside.density[1] = 1028;
  WVDensityEventEvaluation failed;
  outside.bind(failed, WVNoMotionReference::initial);
  requests = outputs(rho, eta, ape);
  require(!failed.evaluate(requests.data(), requests.size()), "out-of-range inversion succeeded");
  untouched(rho); untouched(eta); untouched(ape);
  require(failed.metrics().recoveryCount == 1 && failed.metrics().inversePassCount == 0,
          "late inverse failure lost recovery evidence");
  require(bool(failed.evaluate(&requests[0], 1)) && rho.front() == 1028 &&
              failed.metrics().recoveryCount == 1,
          "successful prior stage was lost after a later failure");

  // The last parcel fails, so no earlier parcel can have reached caller output.
  Fixture nonfinite;
  nonfinite.density.back() = std::numeric_limits<double>::quiet_NaN();
  WVDensityEventEvaluation late;
  nonfinite.bind(late, WVNoMotionReference::initial);
  require(!late.evaluate(&requests[1], 1), "nonfinite last parcel accepted");
  untouched(eta);
  require(late.metrics().inverseSampleCount == 6 && late.metrics().inversePassCount == 0,
          "late inverse attempt count was lost");

  Fixture overflow;
  overflow.gravity = 1e308;
  WVDensityEventEvaluation energyFailure;
  overflow.bind(energyFailure, WVNoMotionReference::initial);
  require(!energyFailure.evaluate(&requests[1], 2), "overflowing APE succeeded");
  untouched(eta); untouched(ape);
  require(energyFailure.metrics().inversePassCount == 1 && energyFailure.metrics().apePassCount == 0,
          "APE failure discarded or published the wrong stage");
  require(bool(energyFailure.evaluate(&requests[1], 1)) && eta[2] == .5 &&
              energyFailure.metrics().inversePassCount == 1,
          "valid eta could not reuse an inverse after APE failure");
}

void recoveryFitAndBudgetFailure() {
  Fixture fixture;
  fixture.weights = {2.0 / 3, 2.0 / 3, 2.0 / 3};
  fixture.density = {1027, 1025, 1026.5, 1026.5, 1025, 1027};
  WVDensityEventEvaluation event;
  fixture.bind(event, WVNoMotionReference::actual);
  std::vector<double> rho(3, -719);
  WVDensityEventOutput request{WVDensityEventField::rhoNm, rho.data(), rho.size()};
  require(bool(event.evaluate(&request, 1)), "nontrivial moment recovery failed");
  close(rho[1], 1026.5, 2e-7, "helper fitted the wrong rearranged density distribution");
  require(event.recoveryReport().qualified && event.recoveryReport().evaluations > 1 &&
              event.metrics().recoveryCount == 1 && event.metrics().inversePassCount == 0 &&
              event.metrics().recoveryWorkspaceHighWaterBytes == event.recoveryReport().workspaceBytes,
          "helper recovery diagnostics or demand isolation changed");
  WVNoMotionRecoveryOptions exhausted;
  exhausted.maximumIterations = 0;
  WVDensityEventEvaluation failure;
  fixture.bind(failure, WVNoMotionReference::actual, exhausted);
  std::fill(rho.begin(), rho.end(), -719);
  require(!failure.evaluate(&request, 1), "exhausted recovery published a profile");
  untouched(rho);
  require(failure.recoveryReport().exitFlag == 0 && !failure.recoveryReport().qualified &&
              failure.metrics().recoveryAttemptCount == 1 && failure.metrics().recoveryCount == 0 &&
              failure.metrics().highWaterBytes >= failure.recoveryReport().workspaceBytes,
          "unqualified fit lost failure diagnostics or peak workspace");
}

void profileConstructionMetricsAndAllocationFailure() {
  Fixture fixture;
  WVNoMotionProfile profile;
  std::size_t peak = 0;
  require(bool(WVNoMotionProfile::create(fixture.z, fixture.initial, profile, &peak)),
          "profile with construction accounting failed");
  const auto retainedVectors = profile.retainedBytes() - sizeof(profile);
  require(peak == retainedVectors + sizeof(double) * (3 * fixture.z.size() - 2),
          "profile construction omitted temporary vector capacities");
  double rho = 0;
  require(bool(profile.density(-.5, rho)) && rho == 1025.5,
          "profile metrics extension changed density arithmetic");
  allocationProbe::failAfter = 3;
  const auto failed = WVNoMotionProfile::create(fixture.z, fixture.initial, profile, &peak);
  allocationProbe::failAfter = -1;
  require(failed.code == WVKernelStatusCode::allocationFailure && peak == 3 * 3 * sizeof(double),
          "partial profile allocation was not reported");
  require(bool(profile.density(-.5, rho)) && rho == 1025.5,
          "failed construction changed the previous profile");

  WVDensityEventEvaluation event;
  fixture.bind(event, WVNoMotionReference::initial);
  std::vector<double> eta(6, -719);
  WVDensityEventOutput request{WVDensityEventField::etaTrue, eta.data(), eta.size()};
  allocationProbe::failAfter = 3;
  auto status = event.evaluate(&request, 1);
  allocationProbe::failAfter = -1;
  require(status.code == WVKernelStatusCode::allocationFailure && event.metrics().liveBytes == 0 &&
              event.metrics().highWaterBytes == 3 * 3 * sizeof(double),
          "helper omitted failed nested construction peak");
  untouched(eta);
  require(bool(event.evaluate(&request, 1)) && eta[2] == .5 &&
              event.metrics().profileConstructionAttemptCount == 2 &&
              event.metrics().profileConstructionCount == 1,
          "profile allocation failure prevented safe retry");

  WVDensityEventEvaluation materialFailure;
  fixture.bind(materialFailure, WVNoMotionReference::initial);
  std::fill(eta.begin(), eta.end(), -719);
  allocationProbe::failAfter = 7;
  status = materialFailure.evaluate(&request, 1);
  allocationProbe::failAfter = -1;
  require(status.code == WVKernelStatusCode::allocationFailure &&
              materialFailure.metrics().profileConstructionCount == 1 &&
              materialFailure.metrics().inversePassCount == 0,
          "material-height allocation failure corrupted completed profile stage");
  untouched(eta);
  require(bool(materialFailure.evaluate(&request, 1)) &&
              materialFailure.metrics().profileConstructionCount == 1 &&
              materialFailure.metrics().inverseAttemptCount == 2,
          "retry rebuilt an already completed selected profile");
}

void bindingAndLifecycle() {
  Fixture fixture;
  WVDensityEventEvaluation event;
  require(bool(event.prepare(0)) && bool(event.evaluate(nullptr, 0)) &&
              event.metrics().liveBytes == 0,
          "empty density demand performed work");
  std::vector<double> rho(3, -719);
  WVDensityEventOutput request{WVDensityEventField::rhoNm, rho.data(), rho.size()};
  require(!event.evaluate(&request, 1), "uninitialized event evaluated density");
  fixture.bind(event, WVNoMotionReference::actual);
  auto invalidGeometry = fixture.geometry();
  invalidGeometry.depth = 0;
  require(!WVDensityEventEvaluation::create(fixture.volume(), invalidGeometry, {}, event) &&
              event.initialized(), "failed binding destroyed previous valid event");
  require(bool(event.evaluate(&request, 1)) && rho[1] == 1026.5, "first event failed");
  const auto peak = event.metrics().highWaterBytes;
  const auto recoveries = event.metrics().recoveryCount;
  const auto *samePointer = fixture.density.data();
  event.release();
  require(!event.initialized() && event.metrics().liveBytes == 0 &&
              event.metrics().highWaterBytes == peak && !event.view(WVDensityEventField::rhoNm).data,
          "event release retained derived state or lost metrics");
  fixture.density[2] = fixture.density[3] = 1026.75;
  require(fixture.density.data() == samePointer, "lifecycle fixture changed its data address");
  fixture.bind(event, WVNoMotionReference::actual);
  require(bool(event.evaluate(&request, 1)) && rho[1] == 1026.75 &&
              event.metrics().recoveryCount == recoveries + 1,
          "fresh event reused stale same-pointer density");
  for (int i = 0; i < 100; ++i) {
    event.release();
    fixture.bind(event, WVNoMotionReference::actual);
    require(bool(event.prepare(7)), "bounded repeated event preparation failed");
    const auto metrics = event.metrics();
    require(metrics.liveBytes <= 128 * sizeof(double) && metrics.highWaterBytes <= 128 * sizeof(double),
            "event history grew retained workspace");
  }
  event.release();
  require(event.metrics().liveBytes == 0, "last event retained workspace");
}

void lowMemoryPreparedScratch() {
  Fixture fixture;
  WVDensityEventEvaluation event;
  constexpr auto derived=static_cast<std::uint8_t>(
      WVDensityEventEvaluation::etaTrueDemand|
      WVDensityEventEvaluation::apeDemand);
  require(bool(event.reserveLowMemoryStorage(
      fixture.density.size(),fixture.z.size(),derived,
      WVNoMotionReference::initial)),
      "low-memory density scratch preparation failed");
  fixture.bind(event,WVNoMotionReference::initial);
  std::vector<double> eta(6),ape(6),etaReplay(6),apeReplay(6);
  WVDensityEventOutput first[]={{WVDensityEventField::etaTrue,
      eta.data(),eta.size()},{WVDensityEventField::ape,ape.data(),ape.size()}};
  require(bool(event.evaluate(first,2)),
      "prepared low-memory density evaluation failed");
  const auto* etaStorage=event.view(WVDensityEventField::etaTrue).data;
  const auto* apeStorage=event.view(WVDensityEventField::ape).data;
  require(etaStorage && apeStorage && etaStorage!=apeStorage &&
          !event.materialHeights().data,
      "low-memory density results did not consume the inverse scratch in place");
  const auto peak=event.metrics().highWaterBytes;
  event.discardDerived();
  require(!event.view(WVDensityEventField::etaTrue).data &&
          !event.view(WVDensityEventField::ape).data &&
          event.metrics().liveBytes==2*fixture.density.size()*sizeof(double),
      "low-memory density discard released capacity or retained a ready value");
  WVDensityEventOutput second[]={{WVDensityEventField::etaTrue,
      etaReplay.data(),etaReplay.size()},{WVDensityEventField::ape,
      apeReplay.data(),apeReplay.size()}};
  require(bool(event.evaluate(second,2)) && etaReplay==eta && apeReplay==ape &&
          event.view(WVDensityEventField::etaTrue).data==etaStorage &&
          event.view(WVDensityEventField::ape).data==apeStorage &&
          event.metrics().inversePassCount==2 &&
          event.metrics().apePassCount==2 &&
          event.metrics().highWaterBytes==peak,
      "low-memory density replay allocated, changed values or reused evicted work");
  require(!event.prepare(WVDensityEventEvaluation::rhoNmDemand),
      "low-memory density accepted an unprepared demand");
  event.release();
  require(event.metrics().liveBytes==0,
      "low-memory density release retained prepared capacity");

  WVDensityEventEvaluation rhoOnly;
  require(bool(rhoOnly.reserveLowMemoryStorage(
      fixture.density.size(),fixture.z.size(),
      WVDensityEventEvaluation::rhoNmDemand,WVNoMotionReference::actual)),
      "rho_nm-only low-memory preparation failed");
  fixture.bind(rhoOnly,WVNoMotionReference::actual);
  require(bool(rhoOnly.prepare(WVDensityEventEvaluation::rhoNmDemand)) &&
          rhoOnly.metrics().inversePassCount==0,
      "rho_nm-only low-memory preparation allocated parcel calculus");
  rhoOnly.discardDerived();
  require(rhoOnly.metrics().liveBytes==fixture.z.size()*sizeof(double),
      "rho_nm-only low-memory scratch retained a volume");
}
} // namespace

int main() {
  try {
    lowMemoryPreparedScratch();
    mixedReferencesAndDemandOrder();
    preflightAndTransactionalFailures();
    recoveryFitAndBudgetFailure();
    profileConstructionMetricsAndAllocationFailure();
    bindingAndLifecycle();
    std::cout << "Density event helper passed: mixed reference, demand extension, "
                 "transactionality, qualified recovery, allocation failures, exact workspace and lifecycle.\n";
    return 0;
  } catch (const std::exception &error) {
    allocationProbe::counting = false;
    allocationProbe::failAfter = -1;
    std::cerr << error.what() << '\n';
    return 1;
  }
}
