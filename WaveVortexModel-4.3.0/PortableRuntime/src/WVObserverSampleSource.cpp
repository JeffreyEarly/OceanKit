#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"

#include <utility>

namespace wavevortex::runtime {

bool sameObservationOccurrenceIdentity(
    const WVObservationOccurrenceIdentity &left,
    const WVObservationOccurrenceIdentity &right) noexcept {
  if (left.resolvedObserverRecord == nullptr ||
      right.resolvedObserverRecord == nullptr ||
      left.logicalScheduleRecord == nullptr ||
      right.logicalScheduleRecord == nullptr ||
      left.schedulePayloadSchema == nullptr ||
      right.schedulePayloadSchema == nullptr ||
      left.proposedScheduleCursor == nullptr ||
      right.proposedScheduleCursor == nullptr ||
      left.resolvedSchedulePayload == nullptr ||
      right.resolvedSchedulePayload == nullptr)
    return false;
  return sameOutputObserverSemanticIdentity(*left.resolvedObserverRecord,
                                            *right.resolvedObserverRecord) &&
         sameLogicalOutputScheduleIdentity(*left.logicalScheduleRecord,
                                           *right.logicalScheduleRecord) &&
         sameOutputSchedulePayloadSchema(*left.schedulePayloadSchema,
                                         *right.schedulePayloadSchema) &&
         samePortableTypedRecordValue(*left.proposedScheduleCursor,
                                      *right.proposedScheduleCursor) &&
         left.resolvedSchedulePayload->sameValue(
             *right.resolvedSchedulePayload) &&
         left.scheduleOrdinal == right.scheduleOrdinal &&
         left.scheduledTime == right.scheduledTime;
}

bool samePreparedObservationOccurrenceIdentity(
    const WVObservationOccurrenceIdentity &left,
    const WVObservationOccurrenceIdentity &right) noexcept {
  return left.preparationOwner != nullptr &&
         left.preparationOwner == right.preparationOwner &&
         left.preparationGeneration != 0 &&
         left.preparationGeneration == right.preparationGeneration &&
         left.preparedOccurrenceSlot == right.preparedOccurrenceSlot &&
         left.observerOrdinal == right.observerOrdinal &&
         left.semanticScheduleOrdinal == right.semanticScheduleOrdinal &&
         left.scheduleOrdinal == right.scheduleOrdinal &&
         left.scheduledTime == right.scheduledTime;
}

WVKernelStatus WVObserverSampleSource::initialObservationBatch(
    const WVObserverRecord &observer, WVObservationBatch &output) {
  WVObservationSchema schema;
  auto status = observationSchema(observer, schema);
  if (!status)
    return status;
  WVObservationBatch candidate;
  candidate.schemaIdentifier = schema.identifier;
  candidate.schemaVersion = schema.version;
  candidate.kind = WVObservationBatchKind::initial;
  status = validateObservationBatch(schema, candidate);
  if (!status)
    return status;
  output = std::move(candidate);
  return WVKernelStatus::ok();
}

} // namespace wavevortex::runtime
