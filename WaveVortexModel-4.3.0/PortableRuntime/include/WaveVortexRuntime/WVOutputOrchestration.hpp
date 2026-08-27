#pragma once

#include "WaveVortexRuntime/WVCheckpointWriter.hpp"
#include "WaveVortexRuntime/WVOutputSchedule.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace wavevortex::runtime {

// Complete schedule state from which one group resumes. This is independent
// of records already committed in any particular output destination.
struct WVOutputScheduleContinuation {
  WVOutputScheduleContinuation() = default;
  WVOutputScheduleContinuation(std::string file, std::string group,
                               WVOutputScheduleCursor value)
      : fileIdentifier(std::move(file)), groupIdentifier(std::move(group)),
        cursor(std::move(value)) {}
  std::string fileIdentifier;
  std::string groupIdentifier;
  WVOutputScheduleCursor cursor;
};

// One sink-owned unlimited-axis offset. physicalCount is retained separately
// so a partially committed ragged destination cannot be mistaken for valid
// schedule continuation state.
struct WVOutputDestinationAxisProgress {
  std::string axisIdentifier;
  std::size_t committedCount = 0;
  std::size_t physicalCount = 0;
};

// Complete commit state of one destination group. The schedule cursor stored
// with the last committed record is evidence about that destination; it is
// never used as the execution continuation unless the compiler first proves
// that the two complete typed cursors agree.
struct WVOutputDestinationProgress {
  std::string fileIdentifier;
  std::string groupIdentifier;
  std::size_t recordCount = 0;
  bool hasCommittedTime = false;
  double lastCommittedTime = 0.0;
  WVOutputScheduleCursor committedScheduleCursor;
  std::vector<WVOutputDestinationAxisProgress> unlimitedAxes;
};

// A shared observer identity resolved once by the output plan. The same
// record pointer is used whenever an observer is routed to several groups or
// destinations.
struct WVOutputObserverView {
  std::size_t observerOrdinal = 0;
  const WVObserverRecord *record = nullptr;
  const WVResolvedObserver *resolved = nullptr;
};

// Immutable route occurrence for one group on one original schedule ordinal.
// All strings and observer arrays are owned by the plan and remain valid for
// the plan's lifetime.
struct WVOutputRouteView {
  std::size_t fileOrdinal = 0;
  std::size_t groupOrdinal = 0;
  WVOutputScheduleOrdinal scheduleOrdinal = WVNoCommittedOutputOrdinal;
  std::string_view fileIdentifier;
  std::string_view destination;
  std::string_view groupIdentifier;
  std::string_view groupName;
  const WVOutputObserverView *observers = nullptr;
  std::size_t observerCount = 0;
  const WVPortableTypedRecord *proposedScheduleCursor = nullptr;
  // Construction-canonical logical group/schedule identity. Compatible
  // groups in different destinations share this ordinal; destination route
  // identity remains fileOrdinal/groupOrdinal.
  std::size_t semanticScheduleOrdinal = 0;
  const WVOutputSchedulePayloadSchema *schedulePayloadSchema = nullptr;
  const WVOutputSchedulePayload *schedulePayload = nullptr;
  std::uint64_t scheduleCursorIdentity = 0;
  // Borrowed immutable construction record for the logical group/schedule.
  // Equivalent groups in different destinations compare by value, never by
  // this pointer or by the plan-local semanticScheduleOrdinal alone.
  const WVOutputGroupRecord *semanticScheduleRecord = nullptr;
};

// Exact, allocation-free semantic comparisons for construction-resolved
// observer and logical group/schedule records. These compare typed records by
// value and perform no named lookup.
bool sameOutputObserverSemanticIdentity(const WVObserverRecord &left,
                                        const WVObserverRecord &right) noexcept;
bool sameLogicalOutputScheduleIdentity(
    const WVOutputGroupRecord &left,
    const WVOutputGroupRecord &right) noexcept;

struct WVOutputPlannedEventView {
  std::size_t eventOrdinal = 0;
  double scheduledTime = 0.0;
  const WVOutputRouteView *routes = nullptr;
  std::size_t routeCount = 0;
};

struct WVOutputPlanMetrics {
  std::size_t fileCount = 0;
  std::size_t groupCount = 0;
  std::size_t distinctObserverCount = 0;
  // Generated counts are driver metrics. The immutable plan retains no
  // complete-window occurrence list.
  std::size_t maximumCoincidentRouteCount = 0;
  std::size_t retainedStorageBytes = 0;
};

// One arbitrary output occurrence. This is the coefficient-checkpoint bridge
// used by the command-line runtime; regular observing systems use descriptor
// output groups and their evenly spaced schedule lattices.
struct WVExplicitOutputTarget {
  double requestedTime = 0.0;
  std::string destination;
};

// Fully resolved, immutable multi-file/multi-group schedule. create()
// validates every route and resolves each schedule before an integrator or
// output sink can mutate state. Future occurrences are generated lazily.
class WVOutputPlan final {
public:
  WVOutputPlan();
  ~WVOutputPlan();
  WVOutputPlan(WVOutputPlan &&) noexcept;
  WVOutputPlan &operator=(WVOutputPlan &&) noexcept;
  WVOutputPlan(const WVOutputPlan &) = delete;
  WVOutputPlan &operator=(const WVOutputPlan &) = delete;

  static WVKernelStatus
  create(const WVPortableObserverDescriptor &descriptor,
         std::shared_ptr<const WVExtensionCatalog> catalog, double initialTime,
         double finalTime,
         const std::vector<WVOutputScheduleContinuation> &continuations,
         WVOutputPlan &plan);
  static WVKernelStatus
  create(const WVIntegrationStateLayout &layout,
         const WVPortableObserverDescriptor &descriptor,
         std::shared_ptr<const WVExtensionCatalog> catalog, double initialTime,
         double finalTime,
         const std::vector<WVOutputScheduleContinuation> &continuations,
         WVOutputPlan &plan);
  static WVKernelStatus
  createExplicit(const WVIntegrationStateLayout &layout,
                 std::shared_ptr<const WVExtensionCatalog> catalog,
                 double initialTime,
                 double finalTime,
                 const std::vector<WVExplicitOutputTarget> &targets,
                 WVOutputPlan &plan);

  double initialTime() const noexcept;
  double finalTime() const noexcept;
  std::size_t groupCount() const noexcept;
  WVOutputRouteView groupRoute(std::size_t index) const noexcept;
  // Authoring/test compatibility views generated on demand. Production
  // orchestration never calls these methods or retains the full window.
  std::size_t eventCount() const noexcept;
  WVOutputPlannedEventView event(std::size_t index) const noexcept;
  const std::vector<WVOutputScheduleContinuation> &
  initialContinuations() const noexcept;
  const WVOutputPlanMetrics &metrics() const noexcept;
  const WVIntegrationStateLayout &stateLayout() const noexcept;
  std::size_t persistentBytes() const noexcept;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
  friend class WVOutputDriver;
};

enum class WVOutputEventKind : std::uint8_t {
  initial,
  interpolated,
  acceptedEndpoint
};

struct WVOutputEvent {
  std::size_t eventOrdinal = 0;
  double scheduledTime = 0.0;
  WVOutputEventKind kind = WVOutputEventKind::acceptedEndpoint;
  WVIntegrationState state;
  const WVOutputRouteView *routes = nullptr;
  std::size_t routeCount = 0;
};

struct WVOutputDeliveryResult {
  enum class Action : std::uint8_t { continueIntegration, terminate };
  Action action = Action::continueIntegration;
  std::size_t writeCount = 0;
  std::size_t writtenBytes = 0;
};

// Abstract orchestration sink. preflight() must validate all destinations and
// resources without writing output. deliver() receives immutable state and
// route views. A failed delivery must leave its route uncommitted and safe to
// retry with the same event state; a successfully returned delivery is never
// repeated by the driver. NetCDF persistence is intentionally outside this
// issue.
class WVOutputSink {
public:
  virtual ~WVOutputSink() = default;
  virtual WVKernelStatus preflight(const WVOutputPlan &plan) = 0;
  virtual WVKernelStatus deliver(const WVOutputEvent &event,
                                 const WVOutputRouteView &route,
                                 WVOutputDeliveryResult &result) = 0;
};

struct WVOutputDeliveryRecord {
  std::size_t eventOrdinal = 0;
  std::size_t routeOrdinal = 0;
  std::size_t fileOrdinal = 0;
  std::size_t groupOrdinal = 0;
  WVOutputScheduleOrdinal scheduleOrdinal = WVNoCommittedOutputOrdinal;
  double scheduledTime = 0.0;
  WVOutputEventKind eventKind = WVOutputEventKind::acceptedEndpoint;
  std::string fileIdentifier;
  std::string destination;
  std::string groupIdentifier;
  std::size_t observerCount = 0;
  std::size_t writeCount = 0;
  std::size_t writtenBytes = 0;
  std::size_t attemptCount = 0;
  std::size_t failureCount = 0;
  bool attempted = false;
  bool committed = false;
  WVKernelStatusCode failureCode = WVKernelStatusCode::success;
  std::string failure;
};

struct WVOutputGroupMetrics {
  std::string fileIdentifier;
  std::string groupIdentifier;
  std::size_t scheduledDeliveryCount = 0;
  std::size_t attemptedDeliveryCount = 0;
  std::size_t committedDeliveryCount = 0;
  std::size_t writeCount = 0;
  std::size_t writtenBytes = 0;
  std::size_t failureCount = 0;
};

struct WVOutputFileMetrics {
  std::string fileIdentifier;
  std::string destination;
  std::size_t scheduledDeliveryCount = 0;
  std::size_t attemptedDeliveryCount = 0;
  std::size_t committedDeliveryCount = 0;
  std::size_t writeCount = 0;
  std::size_t writtenBytes = 0;
  std::size_t failureCount = 0;
  std::vector<WVOutputGroupMetrics> groups;
};

struct WVOutputDriverMetrics {
  std::size_t generatedEventCount = 0;
  std::size_t generatedRouteCount = 0;
  std::size_t generatedSemanticOccurrenceCount = 0;
  std::size_t maximumCoincidentRouteCount = 0;
  std::size_t acceptedStepCount = 0;
  std::size_t outputStateEvaluationCount = 0;
  std::size_t initialStateEventCount = 0;
  std::size_t interpolatedStateEvaluationCount = 0;
  std::size_t acceptedEndpointStateEventCount = 0;
  std::size_t deliveryAttemptCount = 0;
  std::size_t committedDeliveryCount = 0;
  std::size_t writeCount = 0;
  std::size_t writtenBytes = 0;
  std::size_t failureCount = 0;
  std::size_t interpolationBufferCapacityBytes = 0;
  std::size_t interpolationBufferMaximumLiveBytes = 0;
  std::size_t routeStagingCapacityBytes = 0;
  std::size_t routeStagingMaximumLiveBytes = 0;
  double interpolationSeconds = 0.0;
  std::size_t retainedStorageBytes = 0;
  std::vector<WVOutputFileMetrics> files;
};

// Method-neutral scheduled-output driver. Solver steps are selected without
// consulting output times. One reusable integration state is staged for every
// event, and coincident routes share that single state evaluation. If a route
// fails, the driver retains that immutable state and route cursor; calling
// advanceToTime() again replays only the failed route before integration
// continues from the later accepted state. Delivery records describe only the
// latest selected event (including all attempts of a pending retry); cumulative
// delivery history is reported by metrics instead of retained per occurrence.
class WVOutputDriver final {
public:
  WVOutputDriver(WVTimeIntegrator &integrator, const WVOutputPlan &plan);
  ~WVOutputDriver();
  WVOutputDriver(const WVOutputDriver &) = delete;
  WVOutputDriver &operator=(const WVOutputDriver &) = delete;

  WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                               double finalTime, double initialStepSize,
                               WVOutputSink &sink);

  const std::vector<WVOutputScheduleContinuation> &
  committedContinuations() const noexcept;
  // Bounded by the immutable group count. The returned records belong to the
  // latest selected event and remain stable while that event is pending.
  const std::vector<WVOutputDeliveryRecord> &records() const noexcept;
  const WVOutputDriverMetrics &metrics() const noexcept;
  bool hasPendingDelivery() const noexcept;
  std::size_t persistentBytes() const noexcept;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

struct WVCheckpointOutputRecord {
  std::size_t ordinal = 0;
  double requestedTime = 0.0;
  double emittedTime = 0.0;
  WVOutputEventKind eventKind = WVOutputEventKind::acceptedEndpoint;
  std::string destination;
  double writeSeconds = 0.0;
  bool committed = false;
  std::string failure;
};

struct WVCheckpointOutputSinkMetrics {
  std::size_t receivedEventCount = 0;
  std::size_t checkpointWriteCount = 0;
  std::size_t copiedCoefficientBytes = 0;
  double checkpointWriteSeconds = 0.0;
};

// Transactional coefficient-checkpoint sink for explicit output plans. The
// route owns the destination and scheduled time; the sink owns one reusable
// checkpoint-sized staging buffer. Additional integrated blocks are not yet
// part of the legacy coefficient checkpoint format.
class WVCheckpointOutputSink final : public WVOutputSink {
public:
  WVCheckpointOutputSink(std::shared_ptr<const WVExtensionCatalog> catalog,
                         WVCheckpoint checkpointTemplate);

  WVKernelStatus preflight(const WVOutputPlan &plan) override;
  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override;

  const std::vector<WVCheckpointOutputRecord> &records() const noexcept {
    return records_;
  }
  const WVCheckpointOutputSinkMetrics &metrics() const noexcept {
    return metrics_;
  }
  std::size_t persistentBytes() const noexcept;

private:
  std::shared_ptr<const WVExtensionCatalog> catalog_;
  WVCheckpoint checkpoint_;
  std::vector<WVCheckpointOutputRecord> records_;
  WVCheckpointOutputSinkMetrics metrics_;
  const WVIntegrationStateLayout *layout_ = nullptr;
  bool preflighted_ = false;
};

} // namespace wavevortex::runtime
