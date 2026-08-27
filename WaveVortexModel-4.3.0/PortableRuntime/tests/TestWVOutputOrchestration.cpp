#include "WVTestQuadraticSchedule.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVOutputOrchestration.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::test;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

WVPortableObserverRecord outputRecord(double finalTime = 1.0) {
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back({identifier,
                                  WVStateScalarType::complex64,
                                  {1, 1},
                                  WVToleranceKind::coefficientEnergyScaled,
                                  0.0,
                                  WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
  record.stateBlocks.push_back({"tracerAmplitude",
                                WVStateScalarType::real64,
                                {1, 1, 2},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
  record.stateBlocks.push_back({"diagnostic",
                                WVStateScalarType::real64,
                                {1},
                                WVToleranceKind::uniformAbsolute,
                                1e-10,
                                WVStateOwnership::observerDerived,
                                WVRestartRequirement::derivedState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers.push_back(coefficients);
  WVObserverRecord tracer;
  tracer.identifier = "sharedTracer";
  tracer.name = "Shared tracer";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"tracerAmplitude"};
  record.observers.push_back(tracer);
  record.outputFiles = {{"alpha",
                         "alpha.record",
                         {{"fast",
                           "Fast state",
                           {0.25, 0.0, finalTime},
                           {"coefficients", "sharedTracer"},
                           true},
                          {"slow",
                           "Slow tracer",
                           {0.5, 0.0, finalTime},
                           {"sharedTracer"},
                           false}}},
                        {"beta",
                         "beta.record",
                         {{"offset",
                           "Offset state",
                           {0.5, 0.125, finalTime},
                           {"sharedTracer", "coefficients"},
                           true}}}};
  return record;
}

WVPortableObserverDescriptor
descriptorFrom(const WVPortableObserverRecord &record) {
  WVPortableObserverDescriptor descriptor;
  const auto status = WVPortableObserverDescriptor::create(record, test::extensionCatalog(), descriptor);
  require(static_cast<bool>(status),
          "observer descriptor creation: " + status.message);
  return descriptor;
}

WVPortableObserverRecord
singleScheduleRecord(double interval, double initialTime, double finalTime) {
  auto record = outputRecord();
  record.outputFiles = {{"single",
                         "single.record",
                         {{"group",
                           "Single group",
                           {interval, initialTime, finalTime},
                           {"coefficients"},
                           true}}}};
  return record;
}

inline constexpr const char *emptyCursorScheduleType =
    "WVTestEmptyCursorOutputSchedule";

class EmptyCursorSchedule final : public WVOutputSchedule {
public:
  const char *typeIdentifier() const noexcept override {
    return emptyCursorScheduleType;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  const WVOutputSchedulePayloadSchema &payloadSchema() const noexcept override {
    return emptyOutputSchedulePayloadSchema();
  }
  WVKernelStatus
  validateCursor(const WVOutputScheduleCursor &cursor) const override {
    if (cursor.committedOrdinal == WVNoCommittedOutputOrdinal)
      return WVKernelStatus::ok();
    return {WVKernelStatusCode::invalidConfiguration,
            "test algorithmic schedule requires a typed cursor"};
  }
  WVKernelStatus committedTime(const WVOutputScheduleCursor &, double &,
                               bool &available) const override {
    available = false;
    return WVKernelStatus::ok();
  }
  WVKernelStatus peek(const WVOutputScheduleCursor &cursor, double lowerBound,
                      double, WVOutputScheduleOccurrence &occurrence,
                      bool &available) const override {
    if (cursor.committedOrdinal >= 0) {
      available = false;
      return WVKernelStatus::ok();
    }
    occurrence = {};
    occurrence.scheduledTime = lowerBound;
    occurrence.ordinal = 0;
    occurrence.proposedCursor = {0, {}};
    const auto status = occurrence.payload.reset(payloadSchema());
    if (!status)
      return status;
    occurrence.cursorIdentity = 1;
    available = true;
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
};

std::shared_ptr<const WVOutputSchedule>
makeEmptyCursorSchedule(const WVOutputScheduleRecord &,
                        WVKernelStatus &status) {
  status = WVKernelStatus::ok();
  return std::make_shared<EmptyCursorSchedule>();
}

class LinearSystem final : public WVIntegrationSystem {
public:
  class ErrorPolicy final : public WVIntegrationErrorPolicy {
  public:
    explicit ErrorPolicy(const WVIntegrationStateLayout &layout)
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
    std::size_t persistentBytes() const noexcept override { return 0; }

  private:
    const WVIntegrationStateLayout &layout_;
  };
  explicit LinearSystem(WVIntegrationStateLayout layout)
      : layout_(std::move(layout)) {}
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(const WVIntegrationState &state,
                                       WVIntegrationFlux &rhs) override {
    const WVComplexConstView source[] = {state.waveVortex.coefficients.Ap,
                                         state.waveVortex.coefficients.Am,
                                         state.waveVortex.coefficients.A0};
    WVComplexView destination[] = {rhs.waveVortex.Fp, rhs.waveVortex.Fm,
                                   rhs.waveVortex.F0};
    for (std::size_t component = 0; component < 3; ++component)
      for (std::size_t index = 0;
           index < source[component].shape.elementCount(); ++index)
        destination[component].data[index] = {
            -source[component].data[index].real,
            -source[component].data[index].imag};
    for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
      const auto count = state.additionalBlocks[block].layout->elementCount;
      for (std::size_t index = 0; index < count; ++index)
        rhs.additionalBlocks[block].realData[index] =
            -2.0 * state.additionalBlocks[block].realData[index];
    }
    return WVKernelStatus::ok();
  }
  WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &) override {
    return {WVKernelStatus::ok(), 0, true};
  }
  WVKernelStatus createErrorPolicy(
      double,
      std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override {
    policy = std::make_unique<ErrorPolicy>(layout_);
    return WVKernelStatus::ok();
  }

private:
  WVIntegrationStateLayout layout_;
};

struct StateFixture {
  WVShape2D shape{1, 1};
  std::vector<WVComplex64> coefficients{{1.0, 0.5}, {0.8, -0.2}, {0.3, 0.1}};
  WVAdditionalStateStorage additional;
  WVMutableIntegrationState state;

  StateFixture(const WVIntegrationStateLayout &layout, double time = 0.0) {
    require(static_cast<bool>(additional.initialize(layout)),
            "state storage initialization");
    state = {{time,
              0.0,
              {{coefficients.data(), shape},
               {coefficients.data() + 1, shape},
               {coefficients.data() + 2, shape}}},
             additional.mutableBlocks(),
             additional.blockCount()};
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

  std::vector<double> values() const {
    std::vector<double> result;
    for (const auto value : coefficients) {
      result.push_back(value.real);
      result.push_back(value.imag);
    }
    for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
      const auto &metadata = *state.additionalBlocks[block].layout;
      for (std::size_t index = 0; index < metadata.elementCount; ++index) {
        if (metadata.scalarType == WVStateScalarType::real64)
          result.push_back(state.additionalBlocks[block].realData[index]);
        else {
          result.push_back(
              state.additionalBlocks[block].complexData[index].real);
          result.push_back(
              state.additionalBlocks[block].complexData[index].imag);
        }
      }
    }
    result.push_back(state.waveVortex.t);
    return result;
  }
};

struct DeliveredRoute {
  std::size_t eventOrdinal = 0;
  std::size_t fileOrdinal = 0;
  std::size_t groupOrdinal = 0;
  WVOutputScheduleOrdinal scheduleOrdinal = 0;
  double time = 0.0;
  WVOutputEventKind kind = WVOutputEventKind::acceptedEndpoint;
  const WVObserverRecord *firstObserver = nullptr;
  double firstCoefficientReal = 0.0;
  WVOutputScheduleCursor proposedCursor;
};

class RecordingSink final : public WVOutputSink {
public:
  WVKernelStatus preflight(const WVOutputPlan &plan) override {
    ++preflightAttempts;
    preflightEventCount = plan.eventCount();
    if (preflightFailure)
      return {WVKernelStatusCode::allocationFailure,
              "simulated staging allocation failure"};
    return WVKernelStatus::ok();
  }
  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    ++attempts;
    delivered.push_back(
        {event.eventOrdinal, route.fileOrdinal, route.groupOrdinal,
         route.scheduleOrdinal, event.scheduledTime, event.kind,
         route.observerCount ? route.observers[0].record : nullptr,
         event.state.waveVortex.coefficients.Ap.data[0].real,
         {route.scheduleOrdinal,
          route.proposedScheduleCursor == nullptr
              ? WVPortableTypedRecord{}
              : *route.proposedScheduleCursor}});
    if (failAtAttempt == attempts)
      return {failureCode, failureMessage};
    result.writeCount = route.observerCount;
    result.writtenBytes = route.observerCount * 16;
    if (terminateAtAttempt == attempts)
      result.action = WVOutputDeliveryResult::Action::terminate;
    return WVKernelStatus::ok();
  }

  bool preflightFailure = false;
  std::size_t failAtAttempt = std::numeric_limits<std::size_t>::max();
  std::size_t terminateAtAttempt = std::numeric_limits<std::size_t>::max();
  WVKernelStatusCode failureCode = WVKernelStatusCode::numericalFailure;
  std::string failureMessage = "simulated interruption";
  std::size_t preflightEventCount = 0;
  std::size_t preflightAttempts = 0;
  std::size_t attempts = 0;
  std::vector<DeliveredRoute> delivered;
};

// Bounded test sink used to prove that driver storage is independent of the
// number of delivered occurrences. It retains scalar diagnostics only.
class CountingSink final : public WVOutputSink {
public:
  WVKernelStatus preflight(const WVOutputPlan &) override {
    ++preflightAttempts;
    return WVKernelStatus::ok();
  }
  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    ++attempts;
    lastEventOrdinal = event.eventOrdinal;
    lastScheduleOrdinal = route.scheduleOrdinal;
    if (driver != nullptr)
      maximumDriverRetainedBytes =
          std::max(maximumDriverRetainedBytes, driver->persistentBytes());
    result.writeCount = route.observerCount;
    result.writtenBytes = route.observerCount * sizeof(double);
    return WVKernelStatus::ok();
  }

  WVOutputDriver *driver = nullptr;
  std::size_t preflightAttempts = 0;
  std::size_t attempts = 0;
  std::size_t lastEventOrdinal = 0;
  WVOutputScheduleOrdinal lastScheduleOrdinal = WVNoCommittedOutputOrdinal;
  std::size_t maximumDriverRetainedBytes = 0;
};

class ReentrantPreflightSink final : public WVOutputSink {
public:
  WVKernelStatus preflight(const WVOutputPlan &) override {
    ++preflightAttempts;
    nestedStatus = driver->advanceToTime(*state, finalTime, stepSize, *this);
    return WVKernelStatus::ok();
  }
  WVKernelStatus deliver(const WVOutputEvent &, const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override {
    ++deliveries;
    result.writeCount = route.observerCount;
    return WVKernelStatus::ok();
  }

  WVOutputDriver *driver = nullptr;
  WVMutableIntegrationState *state = nullptr;
  double finalTime = 0.0;
  double stepSize = 0.0;
  WVKernelStatus nestedStatus;
  std::size_t preflightAttempts = 0;
  std::size_t deliveries = 0;
};

// Test-only method proving that output orchestration depends only on the
// method-neutral integrator contract. The method intentionally leaves every
// state value unchanged while advancing accepted endpoints.
class EndpointOnlyIntegrator final : public WVTimeIntegrator {
public:
  explicit EndpointOnlyIntegrator(const WVIntegrationStateLayout &layout)
      : layout_(layout) {}
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) override {
    accepted_ = {};
    return validateMutableIntegrationState(layout_, state);
  }
  WVKernelStatus step(WVMutableIntegrationState &state,
                      double stepSize) override {
    auto status = validateMutableIntegrationState(layout_, state);
    if (!status)
      return status;
    if (!(stepSize > 0.0) || !std::isfinite(stepSize))
      return {WVKernelStatusCode::invalidConfiguration,
              "Endpoint-only test step must be finite and positive."};
    const double initial = state.waveVortex.t;
    state.waveVortex.t += stepSize;
    const auto endpoint = integrationConstView(state, views_);
    accepted_ = {initial,
                 state.waveVortex.t,
                 endpoint,
                 {++acceptedCount_, 0, 0, stepSize, stepSize, stepSize, 0.0},
                 nullptr};
    nextStepSize_ = stepSize;
    return WVKernelStatus::ok();
  }
  WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                               double finalTime, double stepSize) override {
    while (state.waveVortex.t < finalTime) {
      const auto status =
          step(state, std::min(stepSize, finalTime - state.waveVortex.t));
      if (!status)
        return status;
    }
    return WVKernelStatus::ok();
  }
  const WVAcceptedStep *lastAcceptedStep() const noexcept override {
    return acceptedCount_ == 0 ? nullptr : &accepted_;
  }
  double nextStepSize() const noexcept override { return nextStepSize_; }
  std::size_t persistentBytes() const noexcept override { return 0; }

private:
  const WVIntegrationStateLayout &layout_;
  std::vector<WVAdditionalStateBlockConstView> views_;
  WVAcceptedStep accepted_;
  std::size_t acceptedCount_ = 0;
  double nextStepSize_ = 0.0;
};

struct Context {
  WVPortableObserverDescriptor descriptor;
  WVIntegrationStateLayout layout;
  LinearSystem system;

  Context()
      : descriptor(descriptorFrom(outputRecord())), layout([&] {
          WVIntegrationStateLayout value;
          require(static_cast<bool>(WVIntegrationStateLayout::create(
                      {1, 1}, descriptor, value)),
                  "integration layout creation");
          return value;
        }()),
        system(std::move(layout)) {}
};

void testPlanningOrderingIdentityAndMetrics(Context &context) {
  WVOutputPlan plan;
  auto status = WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan);
  require(static_cast<bool>(status), "complete output plan");
  require(plan.metrics().fileCount == 2 && plan.metrics().groupCount == 3 &&
              plan.metrics().distinctObserverCount == 2,
          "plan inventory metrics");
  require(plan.eventCount() == 7,
          "diagnostic schedule view preserves event count");
  require(plan.event(0).scheduledTime == 0.0 &&
              plan.event(6).scheduledTime == 1.0,
          "inclusive schedule endpoints");
  const auto initial = plan.event(0);
  require(initial.routeCount == 2 && initial.routes[0].fileOrdinal == 0 &&
              initial.routes[0].groupOrdinal == 0 &&
              initial.routes[1].groupOrdinal == 1,
          "deterministic file/group ordering");
  const auto sharedFromFast = initial.routes[0].observers[1].record;
  const auto sharedFromSlow = initial.routes[1].observers[0].record;
  const auto sharedFromOtherFile = plan.event(1).routes[0].observers[0].record;
  require(sharedFromFast == sharedFromSlow &&
              sharedFromSlow == sharedFromOtherFile,
          "stable shared observer identity");
  require(plan.persistentBytes() == plan.metrics().retainedStorageBytes,
          "exact plan retained storage");

  auto emptyRecord = outputRecord(3.0);
  emptyRecord.outputFiles[0].groups[0].schedule.initialTime = 2.0;
  emptyRecord.outputFiles[0].groups[1].schedule.initialTime = 2.0;
  emptyRecord.outputFiles[1].groups[0].schedule.initialTime = 2.125;
  auto emptyDescriptor = descriptorFrom(emptyRecord);
  WVOutputPlan empty;
  status = WVOutputPlan::create(emptyDescriptor, test::extensionCatalog(), 0.0, 1.0, {}, empty);
  require(static_cast<bool>(status) && empty.eventCount() == 0,
          "empty bounded schedules");
  StateFixture emptyFixture(context.system.stateLayout());
  WVFixedStepRK4 emptyIntegrator(context.system, {true});
  require(static_cast<bool>(
              emptyIntegrator.prepareStateAfterRestart(emptyFixture.state)),
          "empty-schedule preparation");
  RecordingSink emptySink;
  WVOutputDriver emptyDriver(emptyIntegrator, empty);
  status = emptyDriver.advanceToTime(emptyFixture.state, 1.0, 0.2, emptySink);
  require(static_cast<bool>(status) && emptySink.attempts == 0 &&
              emptyFixture.state.waveVortex.t == 1.0,
          "empty schedules advance without deliveries");

  const double nextAfterOne =
      std::nextafter(1.0, std::numeric_limits<double>::infinity());
  auto distinctRecord = outputRecord();
  distinctRecord.outputFiles = {
      {"distinct",
       "distinct.record",
       {{"first", "First", {1.0, 1.0, 1.0}, {"coefficients"}, true},
        {"second",
         "Second",
         {1.0, nextAfterOne, nextAfterOne},
         {"sharedTracer"},
         false}}}};
  auto distinctDescriptor = descriptorFrom(distinctRecord);
  WVOutputPlan distinct;
  status =
      WVOutputPlan::create(distinctDescriptor, test::extensionCatalog(), 1.0, nextAfterOne, {}, distinct);
  require(static_cast<bool>(status) && distinct.eventCount() == 2 &&
              distinct.event(0).routeCount == 1 &&
              distinct.event(1).routeCount == 1,
          "distinct representable times are not tolerance-aggregated");

  const double largeAnchor = 1.0e16;
  WVPortableObserverDescriptor largeAnchorDescriptor;
  status = WVPortableObserverDescriptor::create(
      singleScheduleRecord(1.0, largeAnchor, largeAnchor + 4.0),
      test::extensionCatalog(), largeAnchorDescriptor);
  WVOutputPlan indistinguishable;
  if (status)
    status = WVOutputPlan::create(largeAnchorDescriptor, test::extensionCatalog(), largeAnchor,
                                  largeAnchor + 4.0, {}, indistinguishable);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "large-anchor indistinguishable ordinals are rejected");

  const double tinyInterval = std::numeric_limits<double>::epsilon() / 4.0;
  WVPortableObserverDescriptor tinyDescriptor;
  status = WVPortableObserverDescriptor::create(
      singleScheduleRecord(tinyInterval, 1.0, nextAfterOne),
      test::extensionCatalog(), tinyDescriptor);
  if (status)
    status = WVOutputPlan::create(tinyDescriptor, test::extensionCatalog(), 1.0, nextAfterOne, {},
                                  indistinguishable);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "tiny-interval indistinguishable ordinals are rejected");
}

void testFixedDeliveryAndExactMetrics(Context &context) {
  StateFixture fixture(context.system.stateLayout());
  WVFixedStepRK4 rk4(context.system, {true});
  require(static_cast<bool>(rk4.prepareStateAfterRestart(fixture.state)),
          "fixed restart preparation");
  WVOutputPlan plan;
  require(static_cast<bool>(
              WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan)),
          "fixed plan");
  RecordingSink sink;
  WVOutputDriver driver(rk4, plan);
  require(driver.metrics().retainedStorageBytes == driver.persistentBytes(),
          "new output-driver retained storage is exact");
  const auto status = driver.advanceToTime(fixture.state, 1.0, 0.2, sink);
  require(static_cast<bool>(status), "fixed scheduled delivery");
  const auto &metrics = driver.metrics();
  require(metrics.acceptedStepCount == 5 &&
              metrics.outputStateEvaluationCount == 7 &&
              metrics.initialStateEventCount == 1 &&
              metrics.interpolatedStateEvaluationCount == 5 &&
              metrics.acceptedEndpointStateEventCount == 1,
          "fixed state-evaluation aggregation metrics");
  require(metrics.deliveryAttemptCount == 10 &&
              metrics.committedDeliveryCount == 10 &&
              metrics.writeCount == 17 && metrics.writtenBytes == 272 &&
              metrics.failureCount == 0,
          "exact delivery/write metrics");
  require(metrics.files[0].committedDeliveryCount == 8 &&
              metrics.files[0].groups[0].committedDeliveryCount == 5 &&
              metrics.files[0].groups[1].committedDeliveryCount == 3 &&
              metrics.files[1].groups[0].committedDeliveryCount == 2,
          "per-file/group delivery metrics");
  require(driver.records().size() == 2 &&
              std::all_of(driver.records().begin(), driver.records().end(),
                          [](const auto &record) {
                            return record.attempted && record.committed &&
                                   record.failure.empty();
                          }) &&
              driver.records()[0].eventOrdinal == 6 &&
              driver.records()[0].routeOrdinal == 8 &&
              driver.records()[1].routeOrdinal == 9,
          "bounded latest-event route records");
  require(metrics.interpolationBufferCapacityBytes > 0 &&
              metrics.interpolationBufferMaximumLiveBytes ==
                  metrics.interpolationBufferCapacityBytes &&
              metrics.retainedStorageBytes == driver.persistentBytes(),
          "bounded staging and exact retained storage");
  std::cout << "METRICS schedule_events=" << metrics.generatedEventCount
            << " schedule_routes=" << metrics.generatedRouteCount
            << " maximum_coincident_routes="
            << plan.metrics().maximumCoincidentRouteCount
            << " plan_retained_bytes=" << plan.metrics().retainedStorageBytes
            << " output_state_evaluations="
            << metrics.outputStateEvaluationCount
            << " interpolations=" << metrics.interpolatedStateEvaluationCount
            << " deliveries=" << metrics.committedDeliveryCount
            << " writes=" << metrics.writeCount
            << " written_bytes=" << metrics.writtenBytes
            << " failures=" << metrics.failureCount
            << " interpolation_capacity_bytes="
            << metrics.interpolationBufferCapacityBytes
            << " driver_retained_bytes=" << metrics.retainedStorageBytes
            << " alpha_deliveries=" << metrics.files[0].committedDeliveryCount
            << " alpha_fast_deliveries="
            << metrics.files[0].groups[0].committedDeliveryCount
            << " alpha_slow_deliveries="
            << metrics.files[0].groups[1].committedDeliveryCount
            << " beta_deliveries=" << metrics.files[1].committedDeliveryCount
            << '\n';
}

void testSegmentedContinuation(Context &context) {
  StateFixture fixture(context.system.stateLayout());
  WVFixedStepRK4 firstIntegrator(context.system, {true});
  require(static_cast<bool>(
              firstIntegrator.prepareStateAfterRestart(fixture.state)),
          "first segment preparation");
  WVOutputPlan firstPlan;
  require(static_cast<bool>(WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 0.5,
                                                 {}, firstPlan)),
          "first segment plan");
  RecordingSink firstSink;
  WVOutputDriver firstDriver(firstIntegrator, firstPlan);
  require(static_cast<bool>(
              firstDriver.advanceToTime(fixture.state, 0.5, 0.2, firstSink)),
          "first segment delivery");

  WVOutputPlan secondPlan;
  require(static_cast<bool>(WVOutputPlan::create(
              context.descriptor, test::extensionCatalog(), 0.5, 1.0,
              firstDriver.committedContinuations(),
              secondPlan)),
          "second segment plan");
  require(secondPlan.eventCount() > 0 &&
              secondPlan.event(0).scheduledTime > 0.5,
          "committed boundary is not duplicated");
  bool foundAnchoredFastRoute = false;
  for (std::size_t eventIndex = 0; eventIndex < secondPlan.eventCount();
       ++eventIndex) {
    const auto event = secondPlan.event(eventIndex);
    for (std::size_t routeIndex = 0; routeIndex < event.routeCount;
         ++routeIndex)
      if (event.routes[routeIndex].fileIdentifier == "alpha" &&
          event.routes[routeIndex].groupIdentifier == "fast") {
        foundAnchoredFastRoute =
            event.routes[routeIndex].scheduleOrdinal == 3 &&
            event.scheduledTime == 0.75;
        break;
      }
    if (foundAnchoredFastRoute)
      break;
  }
  require(foundAnchoredFastRoute,
          "segment remains anchored to original lattice ordinal");
  WVFixedStepRK4 secondIntegrator(context.system, {true});
  require(static_cast<bool>(
              secondIntegrator.prepareStateAfterRestart(fixture.state)),
          "second segment preparation");
  RecordingSink secondSink;
  WVOutputDriver secondDriver(secondIntegrator, secondPlan);
  require(static_cast<bool>(
              secondDriver.advanceToTime(fixture.state, 1.0, 0.2, secondSink)),
          "second segment delivery");
  require(secondDriver.committedContinuations()[0]
                      .cursor.committedOrdinal == 4 &&
              secondDriver.committedContinuations()[1]
                      .cursor.committedOrdinal == 2 &&
              secondDriver.committedContinuations()[2]
                      .cursor.committedOrdinal == 1,
          "segmented committed progress");
}

void testPreflightAndMalformedProgress(Context &context) {
  WVOutputPlan plan;
  require(static_cast<bool>(
              WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan)),
          "preflight plan");
  StateFixture fixture(context.system.stateLayout());
  const auto before = fixture.values();
  WVFixedStepRK4 rk4(context.system, {true});
  require(static_cast<bool>(rk4.prepareStateAfterRestart(fixture.state)),
          "allocation-failure preparation");
  RecordingSink sink;
  sink.preflightFailure = true;
  WVOutputDriver driver(rk4, plan);
  auto status = driver.advanceToTime(fixture.state, 1.0, 0.2, sink);
  require(status.code == WVKernelStatusCode::allocationFailure &&
              fixture.values() == before &&
              rk4.metrics().acceptedStepCount == 0,
          "allocation failure occurs before state mutation");
  sink.preflightFailure = false;
  status = driver.advanceToTime(fixture.state, 1.0, 0.2, sink);
  require(static_cast<bool>(status) && fixture.state.waveVortex.t == 1.0 &&
              sink.preflightAttempts == 2,
          "ordinary preflight failure permits a later retry");

  auto requireDescriptorMismatch = [&](const WVPortableObserverRecord &record,
                                       const std::string &label) {
    auto mismatchDescriptor = descriptorFrom(record);
    WVIntegrationStateLayout mismatchLayout;
    require(static_cast<bool>(WVIntegrationStateLayout::create(
                {1, 1}, mismatchDescriptor, mismatchLayout)),
            label + " layout creation");
    LinearSystem mismatchSystem(std::move(mismatchLayout));
    StateFixture mismatchFixture(mismatchSystem.stateLayout());
    const auto mismatchBefore = mismatchFixture.values();
    WVFixedStepRK4 mismatchIntegrator(mismatchSystem, {true});
    require(static_cast<bool>(mismatchIntegrator.prepareStateAfterRestart(
                mismatchFixture.state)),
            label + " integrator preparation");
    RecordingSink mismatchSink;
    WVOutputDriver mismatchDriver(mismatchIntegrator, plan);
    const auto mismatchStatus = mismatchDriver.advanceToTime(
        mismatchFixture.state, 1.0, 0.2, mismatchSink);
    require(mismatchStatus.code == WVKernelStatusCode::invalidConfiguration &&
                mismatchSink.preflightAttempts == 0 &&
                mismatchIntegrator.metrics().acceptedStepCount == 0 &&
                mismatchFixture.values() == mismatchBefore,
            label + " fails before callbacks or state mutation");
  };

  auto identifierMismatch = outputRecord();
  identifierMismatch.stateBlocks[3].identifier = "otherTracerAmplitude";
  identifierMismatch.observers[1].stateBlockIdentifiers = {
      "otherTracerAmplitude"};
  requireDescriptorMismatch(identifierMismatch,
                            "state-block identifier mismatch");

  auto orderMismatch = outputRecord();
  std::swap(orderMismatch.stateBlocks[0], orderMismatch.stateBlocks[1]);
  requireDescriptorMismatch(orderMismatch, "state-block order mismatch");

  auto typeMismatch = outputRecord();
  typeMismatch.stateBlocks[4].scalarType = WVStateScalarType::complex64;
  requireDescriptorMismatch(typeMismatch, "state-block scalar-type mismatch");

  auto dimensionMismatch = outputRecord();
  dimensionMismatch.stateBlocks[3].dimensions = {2, 1, 1};
  requireDescriptorMismatch(dimensionMismatch,
                            "state-block dimension mismatch");

  auto observerMismatch = outputRecord();
  observerMismatch.observers[1].name = "Different observer identity";
  requireDescriptorMismatch(observerMismatch, "observer descriptor mismatch");

  StateFixture reentrantFixture(context.system.stateLayout());
  WVFixedStepRK4 reentrantIntegrator(context.system, {true});
  require(static_cast<bool>(reentrantIntegrator.prepareStateAfterRestart(
              reentrantFixture.state)),
          "reentrant-preflight preparation");
  WVOutputDriver reentrantDriver(reentrantIntegrator, plan);
  ReentrantPreflightSink reentrantSink;
  reentrantSink.driver = &reentrantDriver;
  reentrantSink.state = &reentrantFixture.state;
  reentrantSink.finalTime = 1.0;
  reentrantSink.stepSize = 0.2;
  status = reentrantDriver.advanceToTime(reentrantFixture.state, 1.0, 0.2,
                                         reentrantSink);
  require(static_cast<bool>(status) &&
              reentrantSink.nestedStatus.code ==
                  WVKernelStatusCode::reentrantExecution &&
              reentrantSink.preflightAttempts == 1 &&
              reentrantSink.deliveries ==
                  reentrantDriver.metrics().generatedRouteCount,
          "reentrant sink preflight is rejected while outer run succeeds");

  auto malformed = plan.initialContinuations();
  malformed[0].fileIdentifier = "wrong";
  WVOutputPlan ignored;
  status =
      WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, malformed, ignored);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "misnamed progress rejected");
  malformed = plan.initialContinuations();
  malformed[0].cursor.committedOrdinal = 3;
  status =
      WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, malformed, ignored);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "future progress rejected");
  malformed.pop_back();
  status =
      WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, malformed, ignored);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "incomplete progress rejected");
}

void testTerminationInterruptionAndLaterRouteFailure(Context &context) {
  WVOutputPlan plan;
  require(static_cast<bool>(
              WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan)),
          "failure plan");

  StateFixture terminated(context.system.stateLayout());
  WVFixedStepRK4 terminatingIntegrator(context.system, {true});
  require(static_cast<bool>(
              terminatingIntegrator.prepareStateAfterRestart(terminated.state)),
          "termination preparation");
  RecordingSink terminatingSink;
  terminatingSink.terminateAtAttempt = 1;
  WVOutputDriver terminatingDriver(terminatingIntegrator, plan);
  require(static_cast<bool>(terminatingDriver.advanceToTime(
              terminated.state, 1.0, 0.2, terminatingSink)) &&
              terminated.state.waveVortex.t == 0.0 &&
              terminatingSink.attempts == 2,
          "termination completes coincident routes then stops cleanly");

  StateFixture laterFailure(context.system.stateLayout());
  WVFixedStepRK4 laterIntegrator(context.system, {true});
  require(static_cast<bool>(
              laterIntegrator.prepareStateAfterRestart(laterFailure.state)),
          "later-route failure preparation");
  RecordingSink laterSink;
  laterSink.failAtAttempt = 2;
  laterSink.failureMessage = "later route failed";
  WVOutputDriver laterDriver(laterIntegrator, plan);
  auto status =
      laterDriver.advanceToTime(laterFailure.state, 1.0, 0.2, laterSink);
  require(!status && laterFailure.state.waveVortex.t == 0.0 &&
              laterDriver.records().size() == 2 &&
              laterDriver.records()[0].committed &&
              laterDriver.records()[1].attempted &&
              !laterDriver.records()[1].committed &&
              laterDriver.records()[1].failure == "later route failed",
          "later coincident route failure records partial delivery");
  const auto failedRouteOrdinal = laterDriver.records()[1].routeOrdinal;
  laterSink.failAtAttempt = std::numeric_limits<std::size_t>::max();
  laterSink.terminateAtAttempt = 3;
  status =
      laterDriver.advanceToTime(laterFailure.state, 1.0, 0.2, laterSink);
  require(static_cast<bool>(status) &&
              laterDriver.records().size() == 2 &&
              laterDriver.records()[0].attemptCount == 1 &&
              laterDriver.records()[1].routeOrdinal == failedRouteOrdinal &&
              laterDriver.records()[1].committed &&
              laterDriver.records()[1].attemptCount == 2 &&
              laterDriver.records()[1].failureCount == 1 &&
              laterDriver.records()[1].failure == "later route failed",
          "bounded latest-event records preserve exact retry diagnostics");

  StateFixture interrupted(context.system.stateLayout());
  WVFixedStepRK4 interruptedIntegrator(context.system, {true});
  require(static_cast<bool>(interruptedIntegrator.prepareStateAfterRestart(
              interrupted.state)),
          "interruption preparation");
  RecordingSink interruptedSink;
  interruptedSink.failAtAttempt = 4;
  interruptedSink.failureMessage = "interrupted";
  WVOutputDriver interruptedDriver(interruptedIntegrator, plan);
  status = interruptedDriver.advanceToTime(interrupted.state, 1.0, 0.2,
                                           interruptedSink);
  require(
      !status && interrupted.state.waveVortex.t == 0.4 &&
          std::abs(interrupted.coefficients[0].real - std::exp(-0.4)) < 1e-5 &&
          interruptedDriver.metrics().committedDeliveryCount == 3 &&
          interruptedDriver.metrics().failureCount == 1,
      "interruption preserves the latest accepted state and partial records");

  StateFixture resumed(context.system.stateLayout());
  WVFixedStepRK4 resumedIntegrator(context.system, {true});
  require(static_cast<bool>(
              resumedIntegrator.prepareStateAfterRestart(resumed.state)),
          "interior resume preparation");
  RecordingSink resumedSink;
  resumedSink.failAtAttempt = 6;
  resumedSink.failureMessage = "interior sibling failed";
  WVOutputDriver resumedDriver(resumedIntegrator, plan);
  status = resumedDriver.advanceToTime(resumed.state, 1.0, 0.4, resumedSink);
  require(!status && resumed.state.waveVortex.t == 0.8 &&
              resumedDriver.hasPendingDelivery() &&
              resumedDriver.records().size() == 2 &&
              resumedDriver.records()[0].committed &&
              resumedDriver.records()[0].attemptCount == 1 &&
              !resumedDriver.records()[1].committed &&
              resumedDriver.records()[1].attemptCount == 1 &&
              resumedDriver.records()[1].failureCount == 1,
          "interior coincident failure retains event state and route cursor");
  resumedSink.failAtAttempt = std::numeric_limits<std::size_t>::max();
  status = resumedDriver.advanceToTime(resumed.state, 1.0, 0.4, resumedSink);
  require(
      static_cast<bool>(status) && !resumedDriver.hasPendingDelivery() &&
          resumed.state.waveVortex.t == 1.0 &&
          resumedDriver.records().size() == 2 &&
          resumedDriver.records()[0].committed &&
          resumedDriver.records()[1].committed &&
          resumedDriver.metrics().deliveryAttemptCount == 11 &&
          resumedDriver.metrics().committedDeliveryCount == 10 &&
          resumedDriver.metrics().failureCount == 1 &&
          resumedDriver.metrics().outputStateEvaluationCount == 7 &&
          resumedDriver.metrics().interpolatedStateEvaluationCount == 5,
      "retry replays only failed route and continues without reinterpolation");
  std::size_t committedSiblingDeliveries = 0;
  std::vector<double> failedRouteValues;
  for (const auto &delivery : resumedSink.delivered) {
    if (delivery.fileOrdinal == 0 && delivery.groupOrdinal == 0 &&
        delivery.scheduleOrdinal == 2)
      ++committedSiblingDeliveries;
    if (delivery.fileOrdinal == 0 && delivery.groupOrdinal == 1 &&
        delivery.scheduleOrdinal == 1)
      failedRouteValues.push_back(delivery.firstCoefficientReal);
  }
  require(committedSiblingDeliveries == 1 && failedRouteValues.size() == 2 &&
              failedRouteValues[0] == failedRouteValues[1],
          "resume preserves exact event state without duplicating sibling");

  StateFixture resumeControl(context.system.stateLayout());
  WVFixedStepRK4 resumeControlIntegrator(context.system, {true});
  require(static_cast<bool>(resumeControlIntegrator.prepareStateAfterRestart(
              resumeControl.state)) &&
              static_cast<bool>(resumeControlIntegrator.advanceToTime(
                  resumeControl.state, 1.0, 0.4)) &&
              resumeControl.values() == resumed.values() &&
              resumeControlIntegrator.metrics().acceptedStepCount ==
                  resumedIntegrator.metrics().acceptedStepCount,
          "failed-route resume preserves accepted steps and final state");
}

struct RunResult {
  std::vector<double> values;
  std::size_t accepted = 0;
  std::size_t rejected = 0;
};

RunResult fixedRun(Context &context, bool withOutput) {
  StateFixture fixture(context.system.stateLayout());
  WVFixedStepRK4 integrator(context.system, {true});
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "fixed invariance preparation");
  if (withOutput) {
    WVOutputPlan plan;
    require(static_cast<bool>(
                WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan)),
            "fixed invariance plan");
    RecordingSink sink;
    WVOutputDriver driver(integrator, plan);
    require(
        static_cast<bool>(driver.advanceToTime(fixture.state, 1.0, 0.2, sink)),
        "fixed output invariance run");
  } else {
    require(
        static_cast<bool>(integrator.advanceToTime(fixture.state, 1.0, 0.2)),
        "fixed control run");
  }
  return {fixture.values(), integrator.metrics().acceptedStepCount,
          integrator.metrics().rejectedStepCount};
}

RunResult adaptiveRun(Context &context, bool withOutput) {
  StateFixture fixture(context.system.stateLayout());
  WVAdaptiveRK23Options options;
  options.relativeTolerance = 1e-7;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepFactor = 2.0;
  WVAdaptiveRK23 integrator(context.system, options);
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "adaptive invariance preparation");
  if (withOutput) {
    WVOutputPlan plan;
    require(static_cast<bool>(
                WVOutputPlan::create(context.descriptor, test::extensionCatalog(), 0.0, 1.0, {}, plan)),
            "adaptive invariance plan");
    RecordingSink sink;
    WVOutputDriver driver(integrator, plan);
    require(
        static_cast<bool>(driver.advanceToTime(fixture.state, 1.0, 0.5, sink)),
        "adaptive output invariance run");
  } else {
    require(
        static_cast<bool>(integrator.advanceToTime(fixture.state, 1.0, 0.5)),
        "adaptive control run");
  }
  return {fixture.values(), integrator.metrics().acceptedStepCount,
          integrator.metrics().rejectedStepCount};
}

void testSolverInvariance(Context &context) {
  const auto fixedControl = fixedRun(context, false);
  const auto fixedOutput = fixedRun(context, true);
  require(fixedControl.values == fixedOutput.values &&
              fixedControl.accepted == fixedOutput.accepted &&
              fixedControl.rejected == fixedOutput.rejected,
          "fixed RK4 output does not change steps or final state");
  const auto adaptiveControl = adaptiveRun(context, false);
  const auto adaptiveOutput = adaptiveRun(context, true);
  require(adaptiveControl.values == adaptiveOutput.values &&
              adaptiveControl.accepted == adaptiveOutput.accepted &&
              adaptiveControl.rejected == adaptiveOutput.rejected &&
              adaptiveOutput.rejected > 0,
          "adaptive RK3(2) output does not change acceptance or final state");
  std::cout << "SOLVER_METRICS fixed_accepted=" << fixedOutput.accepted
            << " fixed_rejected=" << fixedOutput.rejected
            << " adaptive_accepted=" << adaptiveOutput.accepted
            << " adaptive_rejected=" << adaptiveOutput.rejected << '\n';
}

void testRK78LazyDenseOutputAndTransactionalRetry(Context &context) {
  StateFixture fixture(context.system.stateLayout());
  WVAdaptiveRK78Options options;
  options.relativeTolerance = 1.0;
  options.absoluteToleranceScale = 1.0;
  options.maximumStepSize = 0.3;
  options.retainDenseOutput = true;
  WVAdaptiveRK78 integrator(context.system, options);
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "RK78 output-driver preparation");
  WVOutputPlan plan;
  require(static_cast<bool>(WVOutputPlan::createExplicit(
              context.system.stateLayout(), test::extensionCatalog(), 0.0,
              0.3,
              {{0.1, "rk78-first"}, {0.2, "rk78-second"},
               {0.3, "rk78-endpoint"}},
              plan)),
          "RK78 explicit output plan");
  RecordingSink sink;
  sink.failAtAttempt = 1;
  sink.failureMessage = "RK78 transactional interruption";
  WVOutputDriver driver(integrator, plan);
  auto status = driver.advanceToTime(fixture.state, 0.3, 0.3, sink);
  require(!status && fixture.state.waveVortex.t == 0.3 &&
              driver.hasPendingDelivery() &&
              integrator.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 4 &&
              integrator.metrics().denseOutputCacheBuildCount == 1 &&
              integrator.metrics().denseOutputCacheReuseCount == 0,
          "RK78 failed delivery retains one computed interpolant transaction");
  const auto failedValue = sink.delivered.back().firstCoefficientReal;
  sink.failAtAttempt = std::numeric_limits<std::size_t>::max();
  status = driver.advanceToTime(fixture.state, 0.3, 0.3, sink);
  require(static_cast<bool>(status) && !driver.hasPendingDelivery() &&
              sink.delivered.size() == 4 &&
              sink.delivered[1].firstCoefficientReal == failedValue &&
              integrator.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 4 &&
              integrator.metrics().denseOutputCacheBuildCount == 1 &&
              integrator.metrics().denseOutputCacheReuseCount == 1 &&
              integrator.metrics().denseOutputEvaluationCount == 2 &&
              driver.metrics().interpolatedStateEvaluationCount == 2,
          "RK78 retry reuses staged output and cached extension stages");

  StateFixture control(context.system.stateLayout());
  WVAdaptiveRK78Options endpointOptions = options;
  endpointOptions.retainDenseOutput = false;
  WVAdaptiveRK78 endpointIntegrator(context.system, endpointOptions);
  require(static_cast<bool>(
              endpointIntegrator.prepareStateAfterRestart(control.state)) &&
              static_cast<bool>(
                  endpointIntegrator.advanceToTime(control.state, 0.3, 0.3)) &&
              control.values() == fixture.values() &&
              endpointIntegrator.metrics().acceptedStepCount ==
                  integrator.metrics().acceptedStepCount &&
              endpointIntegrator.metrics().rejectedStepCount ==
                  integrator.metrics().rejectedStepCount &&
              endpointIntegrator.metrics()
                      .continuousExtensionRightHandSideEvaluationCount == 0 &&
              endpointIntegrator.metrics()
                      .continuousExtensionWorkspaceMaximumLiveBytes == 0,
          "RK78 requested output does not alter its accepted trajectory");
}

void testIntegratorExtensionBoundary(Context &context) {
  StateFixture fixture(context.system.stateLayout());
  const auto initialValues = fixture.values();
  EndpointOnlyIntegrator integrator(context.system.stateLayout());
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "test-only integrator preparation");
  WVOutputPlan plan;
  WVOutputPlan nullCatalogPlan;
  require(WVOutputPlan::createExplicit(
              context.system.stateLayout(), {}, 0.0, 0.4,
              {{0.2, "unreachable"}}, nullCatalogPlan)
                  .code == WVKernelStatusCode::invalidConfiguration,
          "explicit output planning dereferenced a null catalog");
  require(static_cast<bool>(WVOutputPlan::createExplicit(
              context.system.stateLayout(), test::extensionCatalog(), 0.0, 0.4,
              {{0.2, "test-method-first"}, {0.4, "test-method-second"}}, plan)),
          "test-only integrator output plan");
  RecordingSink sink;
  WVOutputDriver driver(integrator, plan);
  require(
      static_cast<bool>(driver.advanceToTime(fixture.state, 0.4, 0.2, sink)) &&
          sink.delivered.size() == 2 &&
          driver.metrics().acceptedStepCount == 2 &&
          fixture.state.additionalBlockCount != 0,
      "test-only integrator works with the unchanged output driver and "
      "additional state blocks");
  auto finalValues = fixture.values();
  require(initialValues.size() == finalValues.size() &&
              std::equal(initialValues.begin(), initialValues.end() - 1,
                         finalValues.begin()),
          "test-only integrator preserves all coefficient and observer state");
}

struct DriverStorageRun {
  std::size_t generatedEventCount = 0;
  std::size_t generatedRouteCount = 0;
  std::size_t retainedBytes = 0;
  std::size_t maximumObservedRetainedBytes = 0;
  std::size_t interpolationCapacityBytes = 0;
  std::size_t interpolationMaximumLiveBytes = 0;
  std::size_t routeStagingCapacityBytes = 0;
  std::size_t routeStagingMaximumLiveBytes = 0;
};

DriverStorageRun runBoundedDriverWindow(double finalTime) {
  const auto descriptor =
      descriptorFrom(singleScheduleRecord(0.25, 0.0, finalTime));
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::create(
              {1, 1}, descriptor, layout)),
          "bounded-window state layout");
  StateFixture fixture(layout);
  EndpointOnlyIntegrator integrator(layout);
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "bounded-window integrator preparation");
  WVOutputPlan plan;
  require(static_cast<bool>(WVOutputPlan::create(
              descriptor, test::extensionCatalog(), 0.0, finalTime, {}, plan)),
          "bounded-window output plan");
  CountingSink sink;
  WVOutputDriver driver(integrator, plan);
  sink.driver = &driver;
  require(static_cast<bool>(driver.advanceToTime(
              fixture.state, finalTime, 0.25, sink)),
          "bounded-window driver execution");
  const auto &metrics = driver.metrics();
  require(sink.preflightAttempts == 1 &&
              sink.attempts == metrics.generatedRouteCount &&
              driver.records().size() == 1 &&
              driver.records()[0].attempted &&
              driver.records()[0].committed &&
              driver.records()[0].eventOrdinal == sink.lastEventOrdinal &&
              driver.records()[0].scheduleOrdinal ==
                  sink.lastScheduleOrdinal &&
              metrics.retainedStorageBytes == driver.persistentBytes(),
          "bounded-window latest record and exact retained bytes");
  return {metrics.generatedEventCount,
          metrics.generatedRouteCount,
          metrics.retainedStorageBytes,
          sink.maximumDriverRetainedBytes,
          metrics.interpolationBufferCapacityBytes,
          metrics.interpolationBufferMaximumLiveBytes,
          metrics.routeStagingCapacityBytes,
          metrics.routeStagingMaximumLiveBytes};
}

void testDriverStorageBoundedByConfiguration() {
  const auto shortWindow = runBoundedDriverWindow(1.0);
  const auto longWindow = runBoundedDriverWindow(100.0);
  require(shortWindow.generatedEventCount == 5 &&
              longWindow.generatedEventCount == 401 &&
              shortWindow.generatedRouteCount == 5 &&
              longWindow.generatedRouteCount == 401,
          "short and long windows delivered their complete schedules");
  require(shortWindow.retainedBytes == longWindow.retainedBytes &&
              shortWindow.maximumObservedRetainedBytes ==
                  longWindow.maximumObservedRetainedBytes &&
              shortWindow.interpolationCapacityBytes ==
                  longWindow.interpolationCapacityBytes &&
              shortWindow.interpolationMaximumLiveBytes ==
                  longWindow.interpolationMaximumLiveBytes &&
              shortWindow.routeStagingCapacityBytes ==
                  longWindow.routeStagingCapacityBytes &&
              shortWindow.routeStagingMaximumLiveBytes ==
                  longWindow.routeStagingMaximumLiveBytes,
          "driver retained and maximum-live storage is independent of "
          "delivered event count");
}

void testLazyQuadraticSchedule() {
  auto record = outputRecord();
  record.outputFiles.resize(1);
  record.outputFiles[0].groups.resize(1);
  record.outputFiles[0].groups[0].schedule = quadraticSchedule(1.0e12);
  auto descriptor = descriptorFrom(record);
  WVOutputPlan plan;
  auto status = WVOutputPlan::create(descriptor, test::extensionCatalog(), 0.0, 9.0, {}, plan);
  require(static_cast<bool>(status) && plan.eventCount() == 4 &&
              plan.event(0).scheduledTime == 0.0 &&
              plan.event(3).scheduledTime == 9.0,
          "quadratic schedule generates exact nonuniform occurrences");
  const auto retained = plan.persistentBytes();
  record.outputFiles[0].groups[0].schedule = quadraticSchedule(1.0e15);
  descriptor = descriptorFrom(record);
  WVOutputPlan longer;
  status = WVOutputPlan::create(descriptor, test::extensionCatalog(), 0.0, 9.0, {}, longer);
  require(static_cast<bool>(status) && longer.eventCount() == 4 &&
              longer.event(0).scheduledTime == 0.0 &&
              longer.event(3).scheduledTime == 9.0 &&
              longer.persistentBytes() == retained,
          "lazy schedule storage is independent of future event count");

  WVIntegrationStateLayout retryLayout;
  require(static_cast<bool>(WVIntegrationStateLayout::create(
              {1, 1}, descriptor, retryLayout)),
          "quadratic retry layout");
  LinearSystem retrySystem(std::move(retryLayout));
  StateFixture retryState(retrySystem.stateLayout());
  WVFixedStepRK4 retryIntegrator(retrySystem, {true});
  require(static_cast<bool>(
              retryIntegrator.prepareStateAfterRestart(retryState.state)),
          "quadratic retry integrator preparation");
  RecordingSink retrySink;
  retrySink.failAtAttempt = 2;
  WVOutputDriver retryDriver(retryIntegrator, longer);
  status = retryDriver.advanceToTime(retryState.state, 9.0, 2.0, retrySink);
  require(!status && retryDriver.hasPendingDelivery() &&
              retryDriver.committedContinuations()[0]
                      .cursor.committedOrdinal == 0 &&
              retryDriver.committedContinuations()[0]
                      .cursor.values.value("nextOrdinal") != nullptr,
          "algorithmic route failure changed its committed destination "
          "cursor");
  const auto failedCursor = retrySink.delivered.back().proposedCursor;
  retrySink.failAtAttempt = std::numeric_limits<std::size_t>::max();
  status = retryDriver.advanceToTime(retryState.state, 9.0, 2.0, retrySink);
  require(static_cast<bool>(status) && !retryDriver.hasPendingDelivery(),
          "algorithmic route retry did not complete");
  std::size_t firstSuccessfulAttempts = 0;
  std::vector<WVOutputScheduleCursor> retriedCursors;
  for (const auto &delivery : retrySink.delivered) {
    if (delivery.scheduleOrdinal == 0)
      ++firstSuccessfulAttempts;
    if (delivery.scheduleOrdinal == 1)
      retriedCursors.push_back(delivery.proposedCursor);
  }
  const auto sameTypedCursor = [](const WVOutputScheduleCursor &left,
                                  const WVOutputScheduleCursor &right) {
    if (left.committedOrdinal != right.committedOrdinal ||
        left.values.schemaIdentifier != right.values.schemaIdentifier ||
        left.values.schemaVersion != right.values.schemaVersion ||
        left.values.values.size() != right.values.values.size())
      return false;
    for (std::size_t index = 0; index < left.values.values.size(); ++index) {
      const auto &a = left.values.values[index];
      const auto &b = right.values.values[index];
      if (a.name != b.name || a.dimensions != b.dimensions ||
          a.storage != b.storage)
        return false;
    }
    return true;
  };
  require(firstSuccessfulAttempts == 1 && retriedCursors.size() == 2 &&
              sameTypedCursor(failedCursor, retriedCursors[0]) &&
              sameTypedCursor(retriedCursors[0], retriedCursors[1]),
          "algorithmic retry repeated a successful route or changed the full "
          "typed cursor replay");

  WVOutputScheduleContinuation progress{
      record.outputFiles[0].identifier,
      record.outputFiles[0].groups[0].identifier, {2, {}}};
  progress.cursor.values.schemaIdentifier = "quadratic-cursor-v1";
  progress.cursor.values.schemaVersion = 1;
  progress.cursor.values.values.push_back(
      {"nextOrdinal", {}, std::vector<std::int64_t>{3}});
  WVOutputPlan resumed;
  status = WVOutputPlan::create(descriptor, test::extensionCatalog(), 4.0, 9.0, {progress}, resumed);
  require(static_cast<bool>(status) && resumed.eventCount() == 1 &&
              resumed.event(0).scheduledTime == 9.0,
          "quadratic cursor resumes from a large committed ordinal");

  progress.cursor.values.values[0].storage =
      std::vector<std::int64_t>(1024, 3);
  status = WVOutputPlan::create(descriptor, test::extensionCatalog(), 4.0, 9.0, {progress}, resumed);
  require(!status, "oversized or malformed cursors fail preflight");

  record.outputFiles[0].groups[0].schedule = {};
  record.outputFiles[0].groups[0].schedule.typeIdentifier =
      WVStateTriggeredOutputScheduleType;
  record.outputFiles[0].groups[0].schedule.contractVersion = 1;
  WVPortableObserverDescriptor unsupported;
  status = WVPortableObserverDescriptor::create(record, test::extensionCatalog(), unsupported);
  require(status.code == WVKernelStatusCode::unsupportedOperation,
          "state-triggered schedules fail before allocation");
}

void testScheduleOwnsProposedCursorContract() {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (status)
    status = builder.addOutputScheduleFactory(
        {emptyCursorScheduleType, 1, &makeEmptyCursorSchedule});
  std::shared_ptr<const WVExtensionCatalog> catalog;
  if (status)
    status = builder.freeze(catalog);
  require(static_cast<bool>(status), "test schedule catalog construction");

  auto record = outputRecord();
  record.outputFiles.resize(1);
  record.outputFiles[0].groups.resize(1);
  record.outputFiles[0].groups[0].schedule = {};
  record.outputFiles[0].groups[0].schedule.typeIdentifier =
      emptyCursorScheduleType;
  record.outputFiles[0].groups[0].schedule.contractVersion = 1;
  WVPortableObserverDescriptor descriptor;
  require(static_cast<bool>(
              WVPortableObserverDescriptor::create(record, catalog,
                                                   descriptor)),
          "test schedule descriptor construction");
  WVOutputPlan plan;
  require(static_cast<bool>(
              WVOutputPlan::create(descriptor, catalog, 0.0, 1.0, {}, plan)),
          "test schedule output planning");
  WVIntegrationStateLayout layout;
  require(static_cast<bool>(WVIntegrationStateLayout::create(
              {1, 1}, descriptor, layout)),
          "test schedule state layout");
  LinearSystem system(std::move(layout));
  StateFixture fixture(system.stateLayout());
  const auto before = fixture.values();
  WVFixedStepRK4 integrator(system, {true});
  require(static_cast<bool>(integrator.prepareStateAfterRestart(fixture.state)),
          "test schedule integrator preparation");
  RecordingSink sink;
  WVOutputDriver driver(integrator, plan);
  status = driver.advanceToTime(fixture.state, 1.0, 0.2, sink);
  require(status.code == WVKernelStatusCode::invalidConfiguration &&
              fixture.values() == before &&
              integrator.metrics().acceptedStepCount == 0,
          "resolved schedule validates its proposed cursor without identity "
          "dispatch");
}

} // namespace

int main() {
  Context context;
  testPlanningOrderingIdentityAndMetrics(context);
  testFixedDeliveryAndExactMetrics(context);
  testSegmentedContinuation(context);
  testPreflightAndMalformedProgress(context);
  testTerminationInterruptionAndLaterRouteFailure(context);
  testSolverInvariance(context);
  testRK78LazyDenseOutputAndTransactionalRetry(context);
  testIntegratorExtensionBoundary(context);
  testDriverStorageBoundedByConfiguration();
  testLazyQuadraticSchedule();
  testScheduleOwnsProposedCursorContract();
  std::cout << "PASS: unified multi-file/multi-group output orchestration\n";
  return 0;
}
