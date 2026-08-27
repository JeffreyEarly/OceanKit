#pragma once

#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexRuntime/WVObservation.hpp"
#include "WaveVortexRuntime/WVOutputSchedule.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {

// Borrowed construction context shared by every resolved observer. Observer
// implementations declare their schema and coarse evaluation channels here;
// they never receive a persistence adapter or issue NetCDF calls.
struct WVObserverOutputPlanningContext {
  const WVTransformConstantStratificationConfiguration *configuration =
      nullptr;
  const WVStateBlockRecord *stateBlocks = nullptr;
  std::size_t stateBlockCount = 0;
  bool isDynamicsLinear = false;
  // Transform-neutral coefficient and rank contract. The stabilized
  // constant-stratification configuration remains available above for v1
  // observer implementations that require its physical parameters.
  const WVIntegrationStateLayout *stateLayout = nullptr;

  const WVStateBlockRecord *
  stateBlock(const std::string &identifier) const noexcept;
};

enum class WVObserverOutputChannelSource : std::uint8_t {
  coefficient,
  sampledField,
  movingField,
  additionalState,
  occurrenceValue,
  occurrenceField
};

// One coarse value channel requested by a source-linked observer. The central
// evaluator resolves and shares the expensive field work, then makes the
// resulting immutable value available by variableIdentifier once per batch.
struct WVObserverOutputChannel {
  std::string variableIdentifier;
  WVObserverOutputChannelSource source =
      WVObserverOutputChannelSource::sampledField;
  std::string sourceIdentifier;
  WVFieldSamplingRequest sampling;
  std::size_t coefficientFamily = 0;
  double scale = 1.0;
  double offset = 0.0;
  // Resolved once by WVObserverOutputEvaluationService::create().
  std::size_t resolvedVariableIndex = 0;
  std::size_t resolvedValueSlot = 0;
  std::size_t occurrenceValueSlot = 0;
  std::size_t positionSetSlot = 0;
};

struct WVObserverMovingPositionSource {
  std::vector<std::string> stateBlockIdentifiers;
  std::vector<double> fixedZ;
  std::size_t positionCount = 0;
  bool isXYOnly = false;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
};

struct WVObserverOccurrencePositionSetPlan {
  std::string identifier;
  // Construction metadata identifying the observation values that persist
  // this interpolation geometry. The evaluator resolves them to numeric
  // occurrence-value slots once and compares only those slots at events.
  std::string sampleTimeVariableIdentifier;
  std::string xVariableIdentifier;
  std::string yVariableIdentifier;
  std::string zVariableIdentifier;
  std::size_t resolvedSampleTimeValueSlot = WVNoResolvedObservationVariable;
  std::size_t resolvedXValueSlot = WVNoResolvedObservationVariable;
  std::size_t resolvedYValueSlot = WVNoResolvedObservationVariable;
  std::size_t resolvedZValueSlot = WVNoResolvedObservationVariable;
};

struct WVObserverOccurrenceStateBlockPlan {
  // Construction metadata only. Event preparation receives the corresponding
  // observer-scoped views in this declared order.
  std::string identifier;
  std::size_t resolvedAdditionalStateBlockIndex =
      WVNoResolvedObservationVariable;
};

struct WVObserverOccurrenceValuePlan {
  // Construction-time schema identity. It is resolved to the corresponding
  // schema-variable ordinal before integration.
  std::string variableIdentifier;
  std::size_t resolvedVariableIndex = 0;
};

// Immutable per-record declaration retained by the evaluation service. Values
// in constantValues own their bounded configuration storage; batches borrow it.
struct WVObserverOutputPlan {
  WVObservationSchema schema;
  std::vector<WVObservationValue> constantValues;
  std::vector<WVObserverOutputChannel> channels;
  WVObserverMovingPositionSource movingPositions;
  WVOutputSchedulePayloadSchema occurrencePayloadSchema =
      emptyOutputSchedulePayloadSchema();
  std::vector<WVObserverOccurrenceStateBlockPlan> occurrenceStateBlocks;
  std::vector<WVObserverOccurrencePositionSetPlan> occurrencePositionSets;
  std::vector<WVObserverOccurrenceValuePlan> occurrenceValues;
};

struct WVObserverOccurrencePositionSet {
  std::vector<std::size_t> extents;
  std::vector<double> sampleTimes;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  void clearForReuse() noexcept;
  std::size_t elementCount() const noexcept;
  std::size_t retainedBytes() const noexcept;
  std::size_t liveBytes() const noexcept;
};

struct WVObserverOccurrenceValueStorage {
  WVObservationScalarType scalarType = WVObservationScalarType::real64;
  std::vector<std::size_t> extents;
  std::vector<double> real64;
  std::vector<WVComplex64> complex64;
  std::vector<std::int64_t> integer64;
  std::vector<std::uint8_t> boolean8;
  std::vector<std::string> text;

  void clearForReuse() noexcept;
  std::size_t elementCount() const noexcept;
  std::size_t retainedBytes() const noexcept;
  std::size_t liveBytes() const noexcept;
};

// Evaluator-owned, event-scoped storage. Observer implementations fill this
// object by resolved numeric slots. Storage belongs only to the current
// in-flight occurrence set and is released after every consuming route commits.
struct WVObserverOccurrenceWorkspace {
  std::vector<WVObserverOccurrencePositionSet> positionSets;
  std::vector<WVObserverOccurrenceValueStorage> values;

  void prepareFor(const WVObserverOutputPlan &plan);
  WVKernelStatus resizeReal(std::size_t slot,
                            std::vector<std::size_t> extents,
                            double *&data);
  WVKernelStatus resizeComplex(std::size_t slot,
                               std::vector<std::size_t> extents,
                               WVComplex64 *&data);
  WVKernelStatus resizeInteger(std::size_t slot,
                               std::vector<std::size_t> extents,
                               std::int64_t *&data);
  WVKernelStatus resizeBoolean(std::size_t slot,
                               std::vector<std::size_t> extents,
                               std::uint8_t *&data);
  std::uint64_t geometryFingerprint() const noexcept;
  bool sameGeometry(const WVObserverOccurrenceWorkspace &other) const noexcept;
  std::size_t retainedBytes() const noexcept;
  std::size_t liveBytes() const noexcept;
};

struct WVObserverOccurrencePreparationContext {
  double scheduledTime = 0.0;
  WVOutputScheduleOrdinal scheduleOrdinal = WVNoCommittedOutputOrdinal;
  const WVOutputSchedulePayloadSchema *payloadSchema = nullptr;
  const WVOutputSchedulePayload *payload = nullptr;
  const WVAdditionalStateBlockConstView *observerStateBlocks = nullptr;
  std::size_t observerStateBlockCount = 0;
};

struct WVObserverBorrowedValueView {
  WVObservationScalarType scalarType = WVObservationScalarType::real64;
  const double *real64 = nullptr;
  const WVComplex64 *complex64 = nullptr;
  const std::int64_t *integer64 = nullptr;
  const std::uint8_t *boolean8 = nullptr;
  const std::string *text = nullptr;
  const std::size_t *extents = nullptr;
  std::size_t extentCount = 0;
  std::size_t elementCount = 0;
};

// Event-scoped borrowed context. Calls are coarse (one per declared variable),
// never per grid point, particle, profile sample, or other element.
class WVObserverOutputEvaluationContext {
public:
  virtual ~WVObserverOutputEvaluationContext() = default;
  virtual WVOutputScheduleOrdinal scheduleOrdinal() const noexcept = 0;
  virtual double scheduledTime() const noexcept = 0;
  virtual WVKernelStatus
  value(std::size_t resolvedValueSlot,
        WVObserverBorrowedValueView &output) const = 0;
};

std::size_t
observerOutputPlanRetainedBytes(const WVObserverOutputPlan &plan) noexcept;

} // namespace wavevortex::runtime
