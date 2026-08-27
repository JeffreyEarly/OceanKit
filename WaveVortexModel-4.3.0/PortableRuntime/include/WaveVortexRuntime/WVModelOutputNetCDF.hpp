#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVObservation.hpp"
#include "WaveVortexRuntime/WVOutputOrchestration.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

// Destination-independent identity of one prepared observer occurrence.
// Global event ordinals, preparation/discovery order, and file/group route
// ordinals deliberately do not participate.
struct WVObservationOccurrenceIdentity {
  // Event-scoped, collision-free cache token minted by the sample source
  // after it has compared the complete semantic occurrence, prepared geometry,
  // and resolved field plan. The owner/generation/slot triple is authoritative
  // for in-flight reuse; the fields below remain the destination-independent
  // semantic identity used by diagnostics and segmented-run comparisons.
  const void *preparationOwner = nullptr;
  std::uint64_t preparationGeneration = 0;
  std::size_t preparedOccurrenceSlot = 0;
  // Exact semantic views borrowed from the immutable compiled plan and the
  // currently prepared event. They remain valid only while that prepared
  // event is in flight. Callers that need longer-lived diagnostics may retain
  // the scalar fingerprints below, but must not use fingerprints as exact
  // cache keys.
  const WVObserverRecord *resolvedObserverRecord = nullptr;
  const WVOutputGroupRecord *logicalScheduleRecord = nullptr;
  const WVOutputSchedulePayloadSchema *schedulePayloadSchema = nullptr;
  const WVPortableTypedRecord *proposedScheduleCursor = nullptr;
  const WVOutputSchedulePayload *resolvedSchedulePayload = nullptr;
  // Plan-local ordinals and fingerprints are retained diagnostics only. The
  // semantic comparator uses the exact borrowed views above together with the
  // schedule occurrence ordinal and time below.
  std::size_t observerOrdinal = 0;
  std::size_t semanticScheduleOrdinal = 0;
  WVOutputScheduleOrdinal scheduleOrdinal = WVNoCommittedOutputOrdinal;
  double scheduledTime = 0.0;
  std::uint64_t scheduleCursorIdentity = 0;
  std::uint64_t payloadFingerprint = 0;
  std::uint64_t geometryFingerprint = 0;
  std::uint64_t fieldPlanFingerprint = 0;
};

bool sameObservationOccurrenceIdentity(
    const WVObservationOccurrenceIdentity &left,
    const WVObservationOccurrenceIdentity &right) noexcept;

bool samePreparedObservationOccurrenceIdentity(
    const WVObservationOccurrenceIdentity &left,
    const WVObservationOccurrenceIdentity &right) noexcept;

// Stable source-linked schema/batch boundary. Observer evaluation remains
// independent of NetCDF and may evaluate all coincident routes once in
// prepare().
class WVObserverSampleSource {
public:
  virtual ~WVObserverSampleSource() = default;
  virtual WVKernelStatus observationSchema(
      const WVObserverRecord &observer, WVObservationSchema &output) = 0;
  virtual WVKernelStatus initialObservationBatch(
      const WVObserverRecord &observer, WVObservationBatch &output);
  virtual WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const = 0;
  virtual WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &identity,
      const WVObserverRecord &observer, WVObservationBatch &output) = 0;

  virtual WVKernelStatus preflight(const WVOutputPlan &) {
    return WVKernelStatus::ok();
  }
  virtual WVKernelStatus prepareInitial(const WVState &) {
    return WVKernelStatus::ok();
  }
  virtual WVKernelStatus prepareInitial(const WVIntegrationState &state) {
    return prepareInitial(state.waveVortex);
  }
  virtual WVKernelStatus prepare(const WVOutputEvent &event) = 0;
  // Called exactly once after every destination route for the prepared event
  // has committed. Sources release event-scoped geometry and evaluated data
  // here; failed events deliberately retain them for exact retry.
  virtual void complete(const WVOutputEvent &) noexcept {}
  virtual std::size_t occurrenceWorkspaceRetainedBytes() const noexcept {
    return 0;
  }
  virtual std::size_t occurrenceWorkspaceLiveBytes() const noexcept {
    return 0;
  }
};

struct WVModelOutputNetCDFConfiguration {
  std::shared_ptr<const WVExtensionCatalog> catalog;
  WVCheckpoint checkpointTemplate;
  bool isDynamicsLinear = false;
};

struct WVModelOutputNetCDFMetrics {
  std::size_t fileCount = 0;
  std::size_t groupCount = 0;
  std::size_t initializedFileCount = 0;
  std::size_t committedRecordCount = 0;
  std::size_t synchronizationCount = 0;
  std::size_t writtenBytes = 0;
  double payloadWriteSeconds = 0.0;
  double synchronizationSeconds = 0.0;
  std::size_t failureCount = 0;
  std::size_t retainedStorageBytes = 0;
  std::size_t batchRetainedStorageBytes = 0;
  std::size_t batchMaximumLiveBytes = 0;
  std::size_t occurrenceWorkspaceRetainedBytes = 0;
  std::size_t occurrenceWorkspaceMaximumLiveBytes = 0;
};

struct WVInspectedObservationSchema {
  std::string observerIdentifier;
  WVObservationSchema schema;
};

struct WVModelOutputNetCDFInspection {
  // Allocation-light latest complete coefficient restart among paths. Raw
  // inspection never loads coefficient arrays or constructs implementations.
  WVCheckpointInspection latestRestart;
  std::string latestRestartPath;
  bool isDynamicsLinear = false;
  // Reconstructed canonical observer/output records and declared schemas.
  WVPortableObserverRecord observerRecord;
  std::vector<WVInspectedObservationSchema> observationSchemas;
  // Schedule state at the selected restart is independent of the tail and
  // offsets committed in each destination.
  std::vector<WVOutputScheduleContinuation> scheduleContinuations;
  std::vector<WVOutputDestinationProgress> destinationProgress;
  std::vector<std::string> paths;
};

// MATLAB-compatible multi-file/multi-group persistence.
//
// createNew() fully defines, writes, synchronizes, and closes every sibling
// staging file before the destination set becomes visible. openAppend()
// validates the complete supplied graph, schedules, shapes, record counts,
// time-last markers, cursor state, ragged offsets, and committed payloads
// read-only before reopening the accepted set for mutation. inspect()
// reconstructs the observer graph and allocation-light restart metadata;
// restoreState() later loads the selected coefficient and observer state.
//
// deliver() writes all payloads at the group's next record index, writes time
// last as the commit marker, then synchronizes the file. A failed call leaves
// the route uncommitted and safe to retry with the same immutable event.
class WVModelOutputNetCDFSink final : public WVOutputSink {
public:
  WVModelOutputNetCDFSink();
  ~WVModelOutputNetCDFSink() override;
  WVModelOutputNetCDFSink(WVModelOutputNetCDFSink &&) noexcept;
  WVModelOutputNetCDFSink &operator=(WVModelOutputNetCDFSink &&) noexcept;
  WVModelOutputNetCDFSink(const WVModelOutputNetCDFSink &) = delete;
  WVModelOutputNetCDFSink &operator=(const WVModelOutputNetCDFSink &) = delete;

  // Every public factory capability-preflights sampleSource against the
  // supplied compiled plan before discovering observer schemas or touching a
  // destination.
  static WVCheckpointStatus createNew(
      const WVModelOutputNetCDFConfiguration &configuration,
      const WVPortableObserverDescriptor &descriptor, const WVOutputPlan &plan,
      const WVIntegrationStateLayout &stateLayout,
      WVObserverSampleSource *sampleSource, WVModelOutputNetCDFSink &sink);

  // Stages the complete file set before replacing any destination. Failure
  // restores every original destination byte-for-byte.
  static WVCheckpointStatus replaceExisting(
      const WVModelOutputNetCDFConfiguration &configuration,
      const WVPortableObserverDescriptor &descriptor, const WVOutputPlan &plan,
      const WVIntegrationStateLayout &stateLayout,
      WVObserverSampleSource *sampleSource, WVModelOutputNetCDFSink &sink);

  static WVCheckpointStatus
  openAppend(const WVModelOutputNetCDFConfiguration &configuration,
             const WVPortableObserverDescriptor &descriptor,
             const WVOutputPlan &plan,
             const WVIntegrationStateLayout &stateLayout,
             WVObserverSampleSource *sampleSource,
             const std::vector<WVOutputDestinationProgress>
                 &expectedDestinationProgress,
             WVModelOutputNetCDFSink &sink);

  static WVCheckpointStatus inspect(const std::vector<std::string> &paths,
                                    const WVExtensionCatalog &catalog,
                                    WVModelOutputNetCDFInspection &inspection);

  // Load the selected coefficient and observer-owned restart state only after
  // the canonical output graph has completed capability preflight.
  static WVCheckpointStatus restoreState(
      const WVModelOutputNetCDFInspection &inspection,
      const WVExtensionCatalog &catalog,
      const WVIntegrationStateLayout &stateLayout, WVCheckpoint &checkpoint,
      WVAdditionalStateStorage &additionalState);

  WVKernelStatus preflight(const WVOutputPlan &plan) override;
  WVKernelStatus deliver(const WVOutputEvent &event,
                         const WVOutputRouteView &route,
                         WVOutputDeliveryResult &result) override;

  const std::vector<WVOutputDestinationProgress> &
  destinationProgress() const noexcept;
  const WVModelOutputNetCDFMetrics &metrics() const noexcept;
  WVCheckpointStatus close() noexcept;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace wavevortex::runtime
