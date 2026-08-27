#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <sstream>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

bool sameTime(double first, double second) noexcept {
  return std::isfinite(first) && std::isfinite(second) &&
         std::abs(first - second) <=
             8.0 * std::numeric_limits<double>::epsilon() *
                 std::max({1.0, std::abs(first), std::abs(second)});
}

std::string samplingKey(const std::string &field,
                        const WVFieldSamplingRequest &sampling) {
  std::ostringstream stream;
  stream << field << ':' << static_cast<int>(sampling.kind);
  for (const auto value : sampling.xIndices)
    stream << ":x" << value;
  for (const auto value : sampling.yIndices)
    stream << ":y" << value;
  for (const auto value : sampling.x)
    stream << ":X" << value;
  for (const auto value : sampling.y)
    stream << ":Y" << value;
  for (const auto value : sampling.z)
    stream << ":Z" << value;
  stream << ":i" << static_cast<int>(sampling.interpolation);
  return stream.str();
}

const WVObservationVariable *
findVariable(const WVObservationSchema &schema,
             const std::string &identifier) noexcept {
  const auto found =
      std::find_if(schema.variables.begin(), schema.variables.end(),
                   [&](const auto &candidate) {
                     return candidate.identifier == identifier;
                   });
  return found == schema.variables.end() ? nullptr : &*found;
}

WVKernelStatus fixedExtents(const WVObservationSchema &schema,
                            const WVObservationVariable &variable,
                            std::vector<std::size_t> &extents) {
  extents.clear();
  extents.reserve(variable.dimensionIdentifiers.size());
  for (const auto &identifier : variable.dimensionIdentifiers) {
    const auto axis = std::find_if(
        schema.axes.begin(), schema.axes.end(), [&](const auto &candidate) {
          return candidate.identifier == identifier;
        });
    if (axis == schema.axes.end() ||
        axis->kind != WVObservationAxisKind::fixed)
      return invalid("Shared observer evaluation channels require fixed axes.");
    extents.push_back(axis->extent);
  }
  return WVKernelStatus::ok();
}

std::size_t elementCount(const std::vector<std::size_t> &extents) noexcept {
  std::size_t count = 1;
  for (const auto extent : extents)
    count *= extent;
  return count;
}

WVKernelStatus checkedElementCount(const std::vector<std::size_t> &extents,
                                   std::size_t &count) {
  count = 1;
  for (const auto extent : extents) {
    if (extent != 0 &&
        count > std::numeric_limits<std::size_t>::max() / extent)
      return {WVKernelStatusCode::sizeOverflow,
              "Observation occurrence extents overflow size_t."};
    count *= extent;
  }
  return WVKernelStatus::ok();
}

WVKernelStatus validateOccurrenceWorkspace(
    const WVObserverOutputPlan &plan,
    const WVObserverOccurrenceWorkspace &workspace) {
  if (workspace.positionSets.size() != plan.occurrencePositionSets.size() ||
      workspace.values.size() != plan.occurrenceValues.size())
    return invalid("An observer occurrence returned the wrong workspace shape.");
  std::size_t structuralBytes = 0;
  for (const auto &set : workspace.positionSets) {
    if (set.extents.empty())
      return invalid("Occurrence position extents must be explicit.");
    if (set.extents.size() >
        WVMaximumOutputSchedulePayloadBytes / sizeof(std::size_t) ||
        structuralBytes > WVMaximumOutputSchedulePayloadBytes -
                              set.extents.size() * sizeof(std::size_t))
      return invalid("Observation occurrence metadata exceeds 4 KiB.");
    structuralBytes += set.extents.size() * sizeof(std::size_t);
    std::size_t count = 0;
    auto status = checkedElementCount(set.extents, count);
    if (!status)
      return status;
    if (set.x.size() != count || set.y.size() != count ||
        (!set.z.empty() && set.z.size() != count) ||
        (!set.sampleTimes.empty() && set.sampleTimes.size() != count))
      return invalid("Occurrence coordinates do not match their extents.");
    const auto finite = [](const std::vector<double> &values) {
      return std::all_of(values.begin(), values.end(),
                         [](double value) { return std::isfinite(value); });
    };
    if (!finite(set.x) || !finite(set.y) || !finite(set.z) ||
        !finite(set.sampleTimes))
      return invalid("Occurrence coordinates and sample times must be finite.");
  }
  for (std::size_t slot = 0; slot < workspace.values.size(); ++slot) {
    const auto &value = workspace.values[slot];
    const auto variableIndex = plan.occurrenceValues[slot].resolvedVariableIndex;
    if (variableIndex >= plan.schema.variables.size() ||
        plan.schema.variables[variableIndex].scalarType != value.scalarType)
      return invalid("An occurrence value has the wrong resolved scalar type.");
    if (value.extents.size() >
            WVMaximumOutputSchedulePayloadBytes / sizeof(std::size_t) ||
        structuralBytes > WVMaximumOutputSchedulePayloadBytes -
                              value.extents.size() * sizeof(std::size_t))
      return invalid("Observation occurrence metadata exceeds 4 KiB.");
    structuralBytes += value.extents.size() * sizeof(std::size_t);
    std::size_t count = 0;
    auto status = checkedElementCount(value.extents, count);
    if (!status)
      return status;
    const std::size_t storedCount =
        value.scalarType == WVObservationScalarType::real64
            ? value.real64.size()
            : value.scalarType == WVObservationScalarType::complex64
                  ? value.complex64.size()
                  : value.scalarType == WVObservationScalarType::integer64
                        ? value.integer64.size()
                        : value.scalarType == WVObservationScalarType::boolean8
                              ? value.boolean8.size()
                              : value.text.size();
    if (storedCount != count)
      return invalid("An occurrence value does not match its extents.");
    if (value.scalarType == WVObservationScalarType::boolean8 &&
        std::any_of(value.boolean8.begin(), value.boolean8.end(),
                    [](std::uint8_t item) { return item > 1; }))
      return invalid("Occurrence Boolean values must be zero or one.");
    for (const auto &text : value.text) {
      if (text.size() > WVMaximumOutputSchedulePayloadBytes ||
          structuralBytes > WVMaximumOutputSchedulePayloadBytes - text.size())
        return invalid("Observation occurrence metadata exceeds 4 KiB.");
      structuralBytes += text.size();
    }
  }
  const auto sameCoordinate = [&](std::size_t valueSlot,
                                  const std::vector<std::size_t> &extents,
                                  const std::vector<double> &coordinates,
                                  const char *name) -> WVKernelStatus {
    if (valueSlot == WVNoResolvedObservationVariable)
      return WVKernelStatus::ok();
    if (valueSlot >= workspace.values.size())
      return invalid(std::string("Occurrence ") + name +
                     " coordinate has no resolved value slot.");
    const auto &value = workspace.values[valueSlot];
    if (value.scalarType != WVObservationScalarType::real64 ||
        value.extents != extents || value.real64 != coordinates)
      return invalid(std::string("Persisted occurrence ") + name +
                     " coordinates differ from interpolation geometry.");
    return WVKernelStatus::ok();
  };
  for (std::size_t slot = 0; slot < plan.occurrencePositionSets.size();
       ++slot) {
    const auto &binding = plan.occurrencePositionSets[slot];
    const auto &set = workspace.positionSets[slot];
    auto status = sameCoordinate(binding.resolvedSampleTimeValueSlot,
                                 set.extents, set.sampleTimes, "sample-time");
    if (status)
      status = sameCoordinate(binding.resolvedXValueSlot, set.extents, set.x,
                              "x");
    if (status)
      status = sameCoordinate(binding.resolvedYValueSlot, set.extents, set.y,
                              "y");
    if (status)
      status = sameCoordinate(binding.resolvedZValueSlot, set.extents, set.z,
                              "z");
    if (!status)
      return status;
  }
  return WVKernelStatus::ok();
}

} // namespace

class WVObserverOutputEvaluationService::Impl {
public:
  struct Output {
    std::vector<std::size_t> extents;
    WVObservationScalarType scalarType = WVObservationScalarType::real64;
    WVObserverOutputChannelSource source =
        WVObserverOutputChannelSource::sampledField;
    std::size_t coefficientFamily = 0;
    std::size_t fieldOutput = 0;
    std::size_t additionalStateBlockIndex =
        std::numeric_limits<std::size_t>::max();
    bool initialField = false;
    std::vector<double> affineStorage;
    double scale = 1.0;
    double offset = 0.0;
  };

  struct PreparedOccurrence {
    WVObservationOccurrenceIdentity identity;
    const WVPortableTypedRecord *proposedCursor = nullptr;
    WVOutputSchedulePayload payload;
    const WVObserverRecord *observer = nullptr;
    const WVObserverOutputPlan *plan = nullptr;
    const WVEventFieldEvaluationPlan *eventFieldPlan = nullptr;
    std::vector<Output> *outputs = nullptr;
    std::shared_ptr<const WVObservingSystem> implementation;
    WVObserverOccurrenceWorkspace workspace;
    std::vector<WVAdditionalStateBlockConstView> observerStateViews;
    std::vector<WVEventPositionSetView> positionViews;
    WVPreparedFieldGeometry fieldGeometry;
    std::vector<std::vector<double>> fieldStorage;
    std::vector<WVFieldOutputView> fieldViews;

    std::size_t retainedBytes() const noexcept {
      std::size_t bytes = sizeof(*this) + workspace.retainedBytes() -
                          sizeof(workspace) + fieldGeometry.retainedBytes() -
                          sizeof(fieldGeometry) +
                          observerStateViews.capacity() *
                              sizeof(WVAdditionalStateBlockConstView) +
                          positionViews.capacity() *
                              sizeof(WVEventPositionSetView) +
                          fieldStorage.capacity() * sizeof(std::vector<double>) +
                          fieldViews.capacity() * sizeof(WVFieldOutputView);
      for (const auto &storage : fieldStorage)
        bytes += storage.capacity() * sizeof(double);
      return bytes;
    }

    std::size_t liveBytes() const noexcept {
      std::size_t bytes = sizeof(*this) + workspace.liveBytes() -
                          sizeof(workspace) + fieldGeometry.liveBytes() -
                          sizeof(fieldGeometry) -
                          fieldGeometry.borrowedCoordinateBytes() +
                          observerStateViews.size() *
                              sizeof(WVAdditionalStateBlockConstView) +
                          positionViews.size() *
                              sizeof(WVEventPositionSetView) +
                          fieldStorage.size() * sizeof(std::vector<double>) +
                          fieldViews.size() * sizeof(WVFieldOutputView);
      for (const auto &storage : fieldStorage)
        bytes += storage.size() * sizeof(double);
      return bytes;
    }
  };

  struct ObserverBinding {
    const WVObserverRecord *record = nullptr;
    const WVObserverOutputPlan *plan = nullptr;
    std::vector<Output> *outputs = nullptr;
    std::shared_ptr<const WVObservingSystem> implementation;
    const WVEventFieldEvaluationPlan *eventFieldPlan = nullptr;
    bool hasMovingChannels = false;
  };

  struct MovingCoordinates {
    std::size_t xBlockIndex = 0;
    std::size_t yBlockIndex = 0;
    std::size_t zBlockIndex = 0;
    std::vector<double> fixedZ;
    std::size_t offset = 0;
    std::size_t count = 0;
    bool isXYOnly = false;
  };

  WVIntegrationStateLayout stateLayout;
  bool isDynamicsLinear = false;
  WVPortableObserverDescriptor descriptor;
  std::unique_ptr<WVFieldEvaluationService> ownedFields;
  WVFieldEvaluationService *fields = nullptr;
  WVFieldEvaluationPlan initialFieldPlan;
  WVFieldEvaluationPlan timeSeriesFieldPlan;
  std::vector<std::vector<double>> initialFieldStorage;
  std::vector<std::vector<double>> timeSeriesFieldStorage;
  std::vector<WVFieldOutputView> initialFieldViews;
  std::vector<WVFieldOutputView> timeSeriesFieldViews;
  WVMovingFieldEvaluationPlan movingFieldPlan;
  std::vector<std::vector<double>> movingFieldStorage;
  std::vector<WVFieldOutputView> movingFieldViews;
  std::vector<MovingCoordinates> movingCoordinates;
  std::vector<double> movingX;
  std::vector<double> movingY;
  std::vector<double> movingZ;
  std::vector<WVObserverOutputPlan> observerPlans;
  std::vector<WVEventFieldEvaluationPlan> eventFieldPlans;
  std::vector<std::vector<Output>> observerOutputs;
  WVIntegrationState preparedIntegrationState;
  std::vector<WVCoefficientFamilyConstView> preparedLegacyCoefficientViews;
  std::size_t preparedEventOrdinal = 0;
  double preparedScheduledTime = 0.0;
  std::uint64_t preparationGeneration = 0;
  std::vector<PreparedOccurrence> preparedOccurrences;
  std::vector<WVEventFieldEvaluationBatchEntry> eventFieldBatchEntries;
  std::vector<ObserverBinding> observerBindings;
  bool prepared = false;
  bool preparedOutputEvent = false;
  bool running = false;

  ObserverBinding *binding(const WVObserverRecord &observer) noexcept {
    const auto found = std::find_if(
        observerBindings.begin(), observerBindings.end(),
        [&](const auto &candidate) {
          return candidate.record == &observer ||
                 (candidate.record != nullptr &&
                  candidate.record->identifier == observer.identifier);
        });
    return found == observerBindings.end() ? nullptr : &*found;
  }

  const PreparedOccurrence *occurrence(
      const WVOutputRouteView &route,
      const WVOutputObserverView &observer) const noexcept {
    if (route.schedulePayload == nullptr ||
        route.proposedScheduleCursor == nullptr)
      return nullptr;
    for (const auto &candidate : preparedOccurrences)
      if (candidate.identity.observerOrdinal == observer.observerOrdinal &&
          candidate.identity.semanticScheduleOrdinal ==
              route.semanticScheduleOrdinal &&
          candidate.identity.scheduleOrdinal == route.scheduleOrdinal &&
          candidate.identity.scheduledTime == preparedScheduledTime &&
          candidate.proposedCursor != nullptr &&
          samePortableTypedRecordValue(*candidate.proposedCursor,
                                       *route.proposedScheduleCursor) &&
          candidate.payload.sameValue(*route.schedulePayload))
        return &candidate;
    return nullptr;
  }

  PreparedOccurrence *occurrence(
      const WVObservationOccurrenceIdentity &identity) noexcept {
    const auto found = std::find_if(
        preparedOccurrences.begin(), preparedOccurrences.end(),
        [&](const auto &candidate) {
          return samePreparedObservationOccurrenceIdentity(candidate.identity,
                                                            identity);
        });
    return found == preparedOccurrences.end() ? nullptr : &*found;
  }

  void updateOccurrenceMetrics(
      WVObserverOutputEvaluationMetrics &metrics) const noexcept {
    std::size_t retained =
        preparedOccurrences.capacity() * sizeof(PreparedOccurrence) +
        eventFieldBatchEntries.capacity() *
            sizeof(WVEventFieldEvaluationBatchEntry);
    std::size_t live =
        preparedOccurrences.size() * sizeof(PreparedOccurrence) +
        eventFieldBatchEntries.size() *
            sizeof(WVEventFieldEvaluationBatchEntry);
    for (const auto &occurrence : preparedOccurrences) {
      retained += occurrence.retainedBytes() - sizeof(PreparedOccurrence);
      live += occurrence.liveBytes() - sizeof(PreparedOccurrence);
    }
    metrics.occurrenceWorkspaceRetainedBytes = retained;
    metrics.occurrenceWorkspaceLiveBytes = live;
    metrics.occurrenceWorkspaceMaximumLiveBytes =
        std::max(metrics.occurrenceWorkspaceMaximumLiveBytes,
                 std::max(retained, live));
  }

  WVKernelStatus evaluate(const WVIntegrationState &state, bool initial,
                          bool evaluateMoving,
                          WVObserverOutputEvaluationMetrics &metrics) {
    if (running)
      return invalid("Observer evaluation is not reentrant.");
    running = true;
    const auto finish = [this]() { running = false; };
    auto &views = initial ? initialFieldViews : timeSeriesFieldViews;
    const auto &plan = initial ? initialFieldPlan : timeSeriesFieldPlan;
    if (!views.empty()) {
      const auto status = fields->evaluate(plan, state, views.data(), views.size());
      if (!status) {
        finish();
        return status;
      }
      ++metrics.fieldEvaluationCount;
    }
    if (!initial && !movingFieldViews.empty() && evaluateMoving) {
      for (const auto &coordinates : movingCoordinates) {
        if (coordinates.xBlockIndex >=
                state.additionalBlockCount ||
            coordinates.yBlockIndex >=
                state.additionalBlockCount ||
            (!coordinates.isXYOnly &&
             coordinates.zBlockIndex >=
                 state.additionalBlockCount)) {
          finish();
          return invalid("Moving observer coordinate state is unavailable.");
        }
        const auto *x = state.additionalBlocks +
                        coordinates.xBlockIndex;
        const auto *y = state.additionalBlocks +
                        coordinates.yBlockIndex;
        const auto *z =
            coordinates.isXYOnly
                ? nullptr
                : state.additionalBlocks +
                      coordinates.zBlockIndex;
        std::copy_n(x->realData, coordinates.count,
                    movingX.data() + coordinates.offset);
        std::copy_n(y->realData, coordinates.count,
                    movingY.data() + coordinates.offset);
        if (coordinates.isXYOnly)
          std::copy(coordinates.fixedZ.begin(), coordinates.fixedZ.end(),
                    movingZ.data() + coordinates.offset);
        else
          std::copy_n(z->realData, coordinates.count,
                      movingZ.data() + coordinates.offset);
      }
      const auto status = fields->evaluateMoving(
          movingFieldPlan, state,
          {movingX.data(), movingY.data(), movingZ.data(), movingX.size()},
          movingFieldViews.data(), movingFieldViews.size());
      if (!status) {
        finish();
        return status;
      }
      ++metrics.fieldEvaluationCount;
      ++metrics.routeAwareParticleEvaluationCount;
    } else if (!initial && !movingFieldViews.empty()) {
      ++metrics.skippedParticleEvaluationCount;
    }
    preparedIntegrationState = state;
    preparedLegacyCoefficientViews.clear();
    preparedLegacyCoefficientViews.reserve(stateLayout.coefficientFamilyCount());
    for (std::size_t family = 0;
         family < stateLayout.coefficientFamilyCount(); ++family)
      preparedLegacyCoefficientViews.push_back(
          coefficientFamilyView(stateLayout, state, family));
    preparedIntegrationState.coefficientFamilies =
        preparedLegacyCoefficientViews.data();
    preparedIntegrationState.coefficientFamilyCount =
        preparedLegacyCoefficientViews.size();
    prepared = true;
    ++metrics.preparedEventCount;
    finish();
    return WVKernelStatus::ok();
  }

  WVKernelStatus borrowedValue(
      std::vector<Output> &outputs, std::size_t resolvedValueSlot,
      const PreparedOccurrence *preparedOccurrence,
      WVObserverBorrowedValueView &value,
      WVObserverOutputEvaluationMetrics &metrics) {
    if (resolvedValueSlot >= outputs.size())
      return invalid("Observer output value slot is not part of this service.");
    auto &entry = outputs[resolvedValueSlot];
    if (entry.source == WVObserverOutputChannelSource::occurrenceValue) {
      if (preparedOccurrence == nullptr ||
          entry.fieldOutput >= preparedOccurrence->workspace.values.size())
        return invalid("Occurrence-value storage is unavailable.");
      const auto &storage =
          preparedOccurrence->workspace.values[entry.fieldOutput];
      value = {};
      value.scalarType = storage.scalarType;
      value.extents = storage.extents.data();
      value.extentCount = storage.extents.size();
      value.elementCount = storage.elementCount();
      value.real64 = storage.real64.data();
      value.complex64 = storage.complex64.data();
      value.integer64 = storage.integer64.data();
      value.boolean8 = storage.boolean8.data();
      value.text = storage.text.data();
      return WVKernelStatus::ok();
    }
    if (entry.source == WVObserverOutputChannelSource::occurrenceField) {
      if (preparedOccurrence == nullptr ||
          entry.fieldOutput >= preparedOccurrence->fieldStorage.size() ||
          entry.fieldOutput >= preparedOccurrence->fieldGeometry.outputs().size())
        return invalid("Occurrence-field storage is unavailable.");
      const auto &storage = preparedOccurrence->fieldStorage[entry.fieldOutput];
      value = {};
      value.scalarType = WVObservationScalarType::real64;
      const auto &dimensions = preparedOccurrence->fieldGeometry.outputs()
                                   [entry.fieldOutput]
                                       .dimensions;
      value.extents = dimensions.data();
      value.extentCount = dimensions.size();
      value.elementCount = storage.size();
      value.real64 = storage.data();
      return WVKernelStatus::ok();
    }
    if (!prepared)
      return invalid("Observer values were requested before prepare().");
    value = {};
    value.scalarType = entry.scalarType;
    value.extents = entry.extents.data();
    value.extentCount = entry.extents.size();
    value.elementCount = elementCount(entry.extents);
    if (entry.source == WVObserverOutputChannelSource::coefficient) {
      const auto family = coefficientFamilyView(
          stateLayout, preparedIntegrationState, entry.coefficientFamily);
      if (family.layout == nullptr || family.data == nullptr ||
          family.layout->elementCount != value.elementCount)
        return invalid("Observer coefficient output does not match the "
                       "resolved transform-family layout.");
      value.complex64 = family.data;
      ++metrics.borrowedCoefficientViewCount;
      return WVKernelStatus::ok();
    }
    if (entry.source == WVObserverOutputChannelSource::additionalState) {
      if (entry.additionalStateBlockIndex >=
          preparedIntegrationState.additionalBlockCount)
        return invalid("Observer state output is unavailable.");
      const auto &block = preparedIntegrationState.additionalBlocks
          [entry.additionalStateBlockIndex];
      if (block.realData == nullptr)
        return invalid("Observer state output is not real-valued.");
      value.real64 = block.realData;
      return WVKernelStatus::ok();
    }
    const auto &storage =
        entry.source == WVObserverOutputChannelSource::movingField
            ? movingFieldStorage
            : entry.initialField ? initialFieldStorage
                                 : timeSeriesFieldStorage;
    if (entry.fieldOutput >= storage.size())
      return invalid("Observer field output binding is invalid.");
    if (entry.scale != 1.0 || entry.offset != 0.0) {
      entry.affineStorage.resize(storage[entry.fieldOutput].size());
      std::transform(storage[entry.fieldOutput].begin(),
                     storage[entry.fieldOutput].end(),
                     entry.affineStorage.begin(), [&](double input) {
                       return entry.scale * input + entry.offset;
                     });
      value.real64 = entry.affineStorage.data();
      value.elementCount = entry.affineStorage.size();
    } else {
      value.real64 = storage[entry.fieldOutput].data();
      value.elementCount = storage[entry.fieldOutput].size();
    }
    return WVKernelStatus::ok();
  }
};

WVObserverOutputEvaluationService::~WVObserverOutputEvaluationService() =
    default;

WVKernelStatus WVObserverOutputEvaluationService::create(
    bool isDynamicsLinear, const WVPortableObserverDescriptor &descriptor,
    WVFieldEvaluationService &fieldEvaluationService,
    std::unique_ptr<WVObserverOutputEvaluationService> &service) {
  return create({}, isDynamicsLinear, descriptor, nullptr, service,
                &fieldEvaluationService);
}

WVKernelStatus WVObserverOutputEvaluationService::create(
    const WVTransformConstantStratificationConfiguration &configuration,
    bool isDynamicsLinear, const WVPortableObserverDescriptor &descriptor,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVObserverOutputEvaluationService> &service,
    WVFieldEvaluationService *borrowedFieldEvaluationService) {
  try {
    auto candidate = std::unique_ptr<WVObserverOutputEvaluationService>(
        new WVObserverOutputEvaluationService());
    candidate->impl_ = std::make_unique<Impl>();
    auto &impl = *candidate->impl_;
    impl.isDynamicsLinear = isDynamicsLinear;
    impl.descriptor = descriptor;
    const auto &descriptorRecord = impl.descriptor.record();
    WVKernelStatus status;
    if (borrowedFieldEvaluationService == nullptr) {
      status = WVFieldEvaluationService::create(configuration, std::move(engine),
                                                impl.ownedFields);
      if (!status)
        return status;
      impl.fields = impl.ownedFields.get();
    } else {
      impl.fields = borrowedFieldEvaluationService;
    }
    status = impl.fields->createStateLayout(descriptor, impl.stateLayout);
    if (!status)
      return status;

    WVObserverOutputPlanningContext planningContext;
    planningContext.configuration = impl.fields->hasLegacyConfiguration()
                                        ? &impl.fields->configuration()
                                        : nullptr;
    planningContext.stateBlocks = descriptorRecord.stateBlocks.data();
    planningContext.stateBlockCount = descriptorRecord.stateBlocks.size();
    planningContext.isDynamicsLinear = isDynamicsLinear;
    planningContext.stateLayout = &impl.stateLayout;
    std::vector<WVFieldRequest> initialRequests;
    std::vector<WVFieldRequest> timeSeriesRequests;
    std::map<std::string, std::size_t> initialRequestIndex;
    std::map<std::string, std::size_t> timeSeriesRequestIndex;
    std::vector<WVMovingFieldRequest> movingRequests;
    impl.observerBindings.reserve(descriptorRecord.observers.size());
    impl.observerPlans.reserve(descriptorRecord.observers.size());
    impl.eventFieldPlans.reserve(descriptorRecord.observers.size());
    impl.observerOutputs.reserve(descriptorRecord.observers.size());
    const auto additionalBlockIndex =
        [&](const std::string &identifier, std::size_t &resolvedIndex) {
          resolvedIndex = 0;
          for (const auto &block : descriptorRecord.stateBlocks) {
            const bool coefficientFamily = std::any_of(
                impl.stateLayout.coefficientFamilies().begin(),
                impl.stateLayout.coefficientFamilies().end(),
                [&](const auto &family) {
                  return family.identifier == block.identifier;
                });
            if (coefficientFamily ||
                block.ownership != WVStateOwnership::integratorOwned)
              continue;
            if (block.identifier == identifier)
              return true;
            ++resolvedIndex;
          }
          resolvedIndex = std::numeric_limits<std::size_t>::max();
          return false;
        };

    for (const auto &observer : descriptorRecord.observers) {
      const auto *resolved = impl.descriptor.resolvedObserver(observer);
      if (resolved == nullptr)
        return {WVKernelStatusCode::unsupportedOperation,
                "Observer output evaluation received an unresolved observer."};
      WVObserverOutputPlan plan;
      status = resolved->outputPlan(observer, planningContext, plan);
      if (!status)
        return status;
      status = validateObservationSchema(plan.schema);
      if (!status)
        return status;
      impl.observerPlans.push_back(std::move(plan));
      auto &storedPlan = impl.observerPlans.back();
      impl.observerOutputs.emplace_back();
      auto &outputs = impl.observerOutputs.back();
      for (auto &constant : storedPlan.constantValues) {
        const auto *variable =
            findVariable(storedPlan.schema, constant.variableIdentifier);
        if (variable == nullptr)
          return invalid("Observer constant value references an unknown variable.");
        constant.resolvedVariableIndex = static_cast<std::size_t>(
            variable - storedPlan.schema.variables.data());
      }
      for (auto &occurrenceValue : storedPlan.occurrenceValues) {
        const auto *variable = findVariable(
            storedPlan.schema, occurrenceValue.variableIdentifier);
        if (variable == nullptr)
          return invalid("Occurrence-value plan references an unknown variable.");
        occurrenceValue.resolvedVariableIndex = static_cast<std::size_t>(
            variable - storedPlan.schema.variables.data());
      }
      for (auto &stateBlock : storedPlan.occurrenceStateBlocks) {
        if (stateBlock.identifier.empty() ||
            std::find(observer.stateBlockIdentifiers.begin(),
                      observer.stateBlockIdentifiers.end(),
                      stateBlock.identifier) ==
                observer.stateBlockIdentifiers.end() ||
            !additionalBlockIndex(
                stateBlock.identifier,
                stateBlock.resolvedAdditionalStateBlockIndex))
          return invalid("Occurrence state must resolve to an integrator-owned "
                         "block declared by that observer.");
      }
      std::set<std::string> positionSetIdentifiers;
      for (auto &positionSet : storedPlan.occurrencePositionSets) {
        if (positionSet.identifier.empty() ||
            !positionSetIdentifiers.insert(positionSet.identifier).second)
          return invalid(
              "Occurrence position-set identifiers must be nonempty and unique.");
        const auto resolveCoordinate =
            [&](const std::string &identifier,
                WVObservationCoordinateRole expectedRole, bool required,
                std::size_t &slot) -> WVKernelStatus {
          slot = WVNoResolvedObservationVariable;
          if (identifier.empty())
            return required
                       ? invalid("Occurrence interpolation geometry is missing "
                                 "a required persisted coordinate binding.")
                       : WVKernelStatus::ok();
          const auto found = std::find_if(
              storedPlan.occurrenceValues.begin(),
              storedPlan.occurrenceValues.end(), [&](const auto &candidate) {
                return candidate.variableIdentifier == identifier;
              });
          if (found == storedPlan.occurrenceValues.end())
            return invalid("Occurrence interpolation geometry references an "
                           "unknown persisted coordinate.");
          const auto variableIndex = found->resolvedVariableIndex;
          if (variableIndex >= storedPlan.schema.variables.size())
            return invalid("Occurrence coordinate has no resolved schema slot.");
          const auto &variable = storedPlan.schema.variables[variableIndex];
          if (variable.scalarType != WVObservationScalarType::real64 ||
              variable.coordinateRole != expectedRole)
            return invalid("Occurrence coordinate binding has the wrong type "
                           "or coordinate role.");
          slot = static_cast<std::size_t>(
              found - storedPlan.occurrenceValues.begin());
          return WVKernelStatus::ok();
        };
        status = resolveCoordinate(
            positionSet.sampleTimeVariableIdentifier,
            WVObservationCoordinateRole::sampleTime, false,
            positionSet.resolvedSampleTimeValueSlot);
        if (status)
          status = resolveCoordinate(positionSet.xVariableIdentifier,
                                     WVObservationCoordinateRole::x, true,
                                     positionSet.resolvedXValueSlot);
        if (status)
          status = resolveCoordinate(positionSet.yVariableIdentifier,
                                     WVObservationCoordinateRole::y, true,
                                     positionSet.resolvedYValueSlot);
        if (status)
          status = resolveCoordinate(positionSet.zVariableIdentifier,
                                     WVObservationCoordinateRole::z, false,
                                     positionSet.resolvedZValueSlot);
        if (!status)
          return status;
        const auto &xVariable = storedPlan.schema.variables
            [storedPlan.occurrenceValues[positionSet.resolvedXValueSlot]
                 .resolvedVariableIndex];
        const auto sameDimensions = [&](std::size_t slot) {
          return slot == WVNoResolvedObservationVariable ||
                 storedPlan.schema.variables
                         [storedPlan.occurrenceValues[slot]
                              .resolvedVariableIndex]
                             .dimensionIdentifiers ==
                     xVariable.dimensionIdentifiers;
        };
        if (!sameDimensions(positionSet.resolvedSampleTimeValueSlot) ||
            !sameDimensions(positionSet.resolvedYValueSlot) ||
            !sameDimensions(positionSet.resolvedZValueSlot))
          return invalid("Persisted occurrence coordinates must share one "
                         "logical shape per position set.");
      }
      std::vector<WVEventFieldRequest> eventRequests;

      std::size_t movingOffset = 0;
      const bool hasMovingChannels = std::any_of(
          storedPlan.channels.begin(), storedPlan.channels.end(),
          [](const auto &channel) {
            return channel.source ==
                   WVObserverOutputChannelSource::movingField;
          });
      if (hasMovingChannels) {
        const auto &positions = storedPlan.movingPositions;
        if (positions.positionCount == 0 ||
            positions.stateBlockIdentifiers.size() < 2 ||
            (!positions.isXYOnly &&
             positions.stateBlockIdentifiers.size() < 3) ||
            (!positions.fixedZ.empty() &&
             positions.fixedZ.size() != positions.positionCount))
          return invalid("Moving observer output plan has invalid coordinates.");
        movingOffset = impl.movingX.size();
        impl.movingX.resize(movingOffset + positions.positionCount);
        impl.movingY.resize(movingOffset + positions.positionCount);
        impl.movingZ.resize(movingOffset + positions.positionCount);
        std::size_t xBlockIndex = 0;
        std::size_t yBlockIndex = 0;
        std::size_t zBlockIndex = 0;
        if (!additionalBlockIndex(positions.stateBlockIdentifiers[0],
                                  xBlockIndex) ||
            !additionalBlockIndex(positions.stateBlockIdentifiers[1],
                                  yBlockIndex) ||
            (!positions.isXYOnly &&
             !additionalBlockIndex(positions.stateBlockIdentifiers[2],
                                   zBlockIndex)))
          return invalid("Moving observer coordinate state was not resolved.");
        impl.movingCoordinates.push_back(
            {xBlockIndex, yBlockIndex, zBlockIndex, positions.fixedZ,
             movingOffset, positions.positionCount, positions.isXYOnly});
      }

      for (auto &channel : storedPlan.channels) {
        const auto *variable =
            findVariable(storedPlan.schema, channel.variableIdentifier);
        if (variable == nullptr)
          return invalid("Observer output channel references an unknown variable.");
        const bool occurrenceChannel =
            channel.source == WVObserverOutputChannelSource::occurrenceValue ||
            channel.source == WVObserverOutputChannelSource::occurrenceField;
        if (!occurrenceChannel &&
            variable->scalarType != WVObservationScalarType::real64 &&
            variable->scalarType != WVObservationScalarType::complex64)
          return invalid("Shared evaluation channels support real or complex values.");
        std::vector<std::size_t> extents;
        if (!occurrenceChannel) {
          status = fixedExtents(storedPlan.schema, *variable, extents);
          if (!status)
            return status;
        }
        channel.resolvedVariableIndex = static_cast<std::size_t>(
            variable - storedPlan.schema.variables.data());
        channel.resolvedValueSlot = outputs.size();
        Impl::Output output;
        output.scalarType = variable->scalarType;
        if (!occurrenceChannel)
          output.extents = std::move(extents);
        output.source = channel.source;
        output.coefficientFamily = channel.coefficientFamily;
        if (channel.source ==
                WVObserverOutputChannelSource::additionalState &&
            !additionalBlockIndex(channel.sourceIdentifier,
                                  output.additionalStateBlockIndex))
          return invalid("Observer state output has no resolved block slot.");
        output.scale = channel.scale;
        output.offset = channel.offset;
        if (channel.source == WVObserverOutputChannelSource::sampledField) {
          const bool initial =
              variable->layout == WVObservationValueLayout::initialValue ||
              variable->layout == WVObservationValueLayout::staticValue;
          auto &requests = initial ? initialRequests : timeSeriesRequests;
          auto &requestIndex =
              initial ? initialRequestIndex : timeSeriesRequestIndex;
          const auto key =
              samplingKey(channel.sourceIdentifier, channel.sampling);
          const auto found = requestIndex.find(key);
          if (found == requestIndex.end()) {
            output.fieldOutput = requests.size();
            requestIndex.emplace(key, output.fieldOutput);
            requests.push_back({"field-" + std::to_string(output.fieldOutput),
                                channel.sourceIdentifier, channel.sampling});
          } else {
            output.fieldOutput = found->second;
            ++candidate->metrics_.sharedFieldReuseCount;
          }
          output.initialField = initial;
        } else if (channel.source ==
                   WVObserverOutputChannelSource::movingField) {
          output.fieldOutput = impl.movingFieldStorage.size();
          impl.movingFieldStorage.emplace_back(elementCount(output.extents));
          impl.movingFieldViews.push_back(
              {impl.movingFieldStorage.back().data(),
               impl.movingFieldStorage.back().size()});
          movingRequests.push_back(
              {observer.identifier + '-' + channel.sourceIdentifier,
               channel.sourceIdentifier, movingOffset,
               storedPlan.movingPositions.positionCount,
               storedPlan.movingPositions.interpolation});
        } else if (channel.source ==
                   WVObserverOutputChannelSource::occurrenceValue) {
          const auto sourceIdentifier = channel.sourceIdentifier.empty()
                                            ? channel.variableIdentifier
                                            : channel.sourceIdentifier;
          const auto found = std::find_if(
              storedPlan.occurrenceValues.begin(),
              storedPlan.occurrenceValues.end(), [&](const auto &candidate) {
                return candidate.variableIdentifier == sourceIdentifier;
              });
          if (found == storedPlan.occurrenceValues.end())
            return invalid("Occurrence-value channel has no resolved storage slot.");
          channel.occurrenceValueSlot = static_cast<std::size_t>(
              found - storedPlan.occurrenceValues.begin());
          output.fieldOutput = channel.occurrenceValueSlot;
          if (output.scale != 1.0 || output.offset != 0.0)
            return invalid("Occurrence-value affine transforms are unsupported.");
        } else if (channel.source ==
                   WVObserverOutputChannelSource::occurrenceField) {
          if (channel.positionSetSlot >=
              storedPlan.occurrencePositionSets.size())
            return invalid("Occurrence-field position-set slot is out of range.");
          if (variable->scalarType != WVObservationScalarType::real64)
            return invalid("Occurrence fields must have real-valued outputs.");
          if (output.scale != 1.0 || output.offset != 0.0)
            return invalid("Occurrence-field affine transforms are unsupported.");
          output.fieldOutput = eventRequests.size();
          eventRequests.push_back(
              {channel.variableIdentifier, channel.sourceIdentifier,
               channel.positionSetSlot, channel.sampling.interpolation});
        }
        outputs.push_back(std::move(output));
      }
      WVEventFieldEvaluationPlan eventPlan;
      if (!eventRequests.empty()) {
        status = impl.fields->createEventPlan(eventRequests, eventPlan);
        if (!status)
          return status;
      }
      impl.eventFieldPlans.push_back(std::move(eventPlan));
      impl.observerBindings.push_back(
          {&observer, &storedPlan, &outputs, resolved->implementationHandle(),
           &impl.eventFieldPlans.back(), hasMovingChannels});
    }

    const auto buildPlan = [&](const std::vector<WVFieldRequest> &requests,
                               WVFieldEvaluationPlan &plan,
                               std::vector<std::vector<double>> &storage,
                               std::vector<WVFieldOutputView> &views) {
      if (requests.empty())
        return WVKernelStatus::ok();
      auto planStatus = impl.fields->createPlan(requests, plan);
      if (!planStatus)
        return planStatus;
      storage.resize(plan.outputCount());
      views.resize(plan.outputCount());
      for (std::size_t index = 0; index < plan.outputCount(); ++index) {
        const auto count = plan.outputs()[index].elementCount;
        storage[index].resize(count);
        views[index] = {storage[index].data(), count};
        candidate->metrics_.outputCapacityBytes +=
            storage[index].capacity() * sizeof(double);
      }
      return WVKernelStatus::ok();
    };
    status = buildPlan(initialRequests, impl.initialFieldPlan,
                       impl.initialFieldStorage, impl.initialFieldViews);
    if (!status)
      return status;
    status = buildPlan(timeSeriesRequests, impl.timeSeriesFieldPlan,
                       impl.timeSeriesFieldStorage,
                       impl.timeSeriesFieldViews);
    if (!status)
      return status;
    if (!movingRequests.empty()) {
      status = impl.fields->createMovingPlan(movingRequests,
                                             impl.movingFieldPlan);
      if (!status)
        return status;
      for (const auto &storage : impl.movingFieldStorage)
        candidate->metrics_.outputCapacityBytes +=
            storage.capacity() * sizeof(double);
    }
    candidate->metrics_.uniqueFieldOutputCount =
        initialRequests.size() + timeSeriesRequests.size();
    candidate->metrics_.retainedStorageBytes = candidate->persistentBytes();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate observer-output evaluation storage."};
  }
}

WVKernelStatus WVObserverOutputEvaluationService::observationSchema(
    const WVObserverRecord &observer, WVObservationSchema &output) {
  const auto *binding = impl_->binding(observer);
  if (binding == nullptr || binding->plan == nullptr)
    return invalid("Observer is not part of this evaluation service.");
  output = binding->plan->schema;
  return WVKernelStatus::ok();
}

WVKernelStatus WVObserverOutputEvaluationService::observationBatchForKind(
    const WVObservationOccurrenceIdentity *identity,
    const WVObserverRecord &observer,
    WVObservationBatchKind kind,
    WVObservationBatch &output) {
  if (!impl_->prepared)
    return invalid("Observation batch was requested before prepare().");
  const WVObserverOutputPlan *plan = nullptr;
  std::shared_ptr<const WVObservingSystem> implementation;
  std::vector<Impl::Output> *outputs = nullptr;
  Impl::PreparedOccurrence *preparedOccurrence = nullptr;
  if (identity != nullptr) {
    preparedOccurrence = impl_->occurrence(*identity);
    if (preparedOccurrence == nullptr ||
        preparedOccurrence->observer != &observer ||
        preparedOccurrence->plan == nullptr ||
        preparedOccurrence->outputs == nullptr)
      return invalid("Observation occurrence identity is not prepared.");
    plan = preparedOccurrence->plan;
    implementation = preparedOccurrence->implementation;
    outputs = preparedOccurrence->outputs;
  } else {
    auto *binding = impl_->binding(observer);
    if (binding == nullptr || binding->plan == nullptr ||
        binding->implementation == nullptr || binding->outputs == nullptr)
      return invalid("Observer is not part of this evaluation service.");
    plan = binding->plan;
    implementation = binding->implementation;
    outputs = binding->outputs;
  }
  class Context final : public WVObserverOutputEvaluationContext {
  public:
    Context(Impl &impl, WVObserverOutputEvaluationMetrics &metrics,
            std::vector<Impl::Output> &outputs,
            const Impl::PreparedOccurrence *occurrence)
        : impl_(impl), metrics_(metrics), outputs_(outputs),
          occurrence_(occurrence) {}
    WVOutputScheduleOrdinal scheduleOrdinal() const noexcept override {
      return occurrence_ == nullptr ? WVNoCommittedOutputOrdinal
                                    : occurrence_->identity.scheduleOrdinal;
    }
    double scheduledTime() const noexcept override {
      return impl_.preparedScheduledTime;
    }
    WVKernelStatus value(std::size_t resolvedValueSlot,
                         WVObserverBorrowedValueView &value) const override {
      return impl_.borrowedValue(outputs_, resolvedValueSlot, occurrence_, value,
                                 metrics_);
    }

  private:
    Impl &impl_;
    WVObserverOutputEvaluationMetrics &metrics_;
    std::vector<Impl::Output> &outputs_;
    const Impl::PreparedOccurrence *occurrence_;
  } context(*impl_, metrics_, *outputs, preparedOccurrence);
  WVObservationBatch batch;
  auto status = implementation->observationBatch(observer, *plan, context,
                                                  kind, batch);
  if (!status)
    return status;
  // Initial/static delivery retains the legacy named-schema contract. Event
  // delivery is validated by the sink's construction-compiled numeric schema
  // bindings before any file mutation; invoking the named validator here would
  // restore per-occurrence identifier and axis lookup to the hot path.
  if (identity == nullptr) {
    status = validateObservationBatch(plan->schema, batch);
    if (!status)
      return status;
  }
  const auto batchMetrics = batch.metrics();
  metrics_.batchRetainedStorageBytes =
      std::max(metrics_.batchRetainedStorageBytes,
               batchMetrics.retainedStorageBytes);
  metrics_.batchMaximumLiveBytes =
      std::max(metrics_.batchMaximumLiveBytes, batchMetrics.liveBytes);
  if (identity != nullptr)
    ++metrics_.occurrenceBatchBuildCount;
  output = std::move(batch);
  metrics_.retainedStorageBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVObserverOutputEvaluationService::initialObservationBatch(
    const WVObserverRecord &observer, WVObservationBatch &output) {
  return observationBatchForKind(nullptr, observer,
                                 WVObservationBatchKind::initial, output);
}

WVKernelStatus WVObserverOutputEvaluationService::preparedOccurrenceIdentity(
    const WVOutputRouteView &route, const WVOutputObserverView &observer,
    WVObservationOccurrenceIdentity &output) const {
  if (!impl_ || !impl_->preparedOutputEvent)
    return invalid("No output occurrence is prepared.");
  const auto *prepared = impl_->occurrence(route, observer);
  if (prepared == nullptr)
    return invalid("The requested observation occurrence was not prepared.");
  output = prepared->identity;
  return WVKernelStatus::ok();
}

WVKernelStatus WVObserverOutputEvaluationService::observationBatch(
    const WVObservationOccurrenceIdentity &identity,
    const WVObserverRecord &observer, WVObservationBatch &output) {
  return observationBatchForKind(&identity, observer,
                                 WVObservationBatchKind::event, output);
}

WVKernelStatus WVObserverOutputEvaluationService::preflight(
    const WVOutputPlan &plan) {
  for (std::size_t groupIndex = 0; groupIndex < plan.groupCount(); ++groupIndex) {
    const auto route = plan.groupRoute(groupIndex);
    for (std::size_t observerIndex = 0;
         observerIndex < route.observerCount; ++observerIndex) {
      const auto &resolved = route.observers[observerIndex];
      const auto ordinal = resolved.observerOrdinal;
      if (resolved.record == nullptr || resolved.resolved == nullptr ||
          ordinal >= impl_->observerBindings.size() ||
          impl_->observerBindings[ordinal].record != resolved.record ||
          impl_->observerBindings[ordinal].plan == nullptr)
        return invalid("Output plan references an observer outside this "
                       "evaluation service.");
      const auto *observerPlan = impl_->observerBindings[ordinal].plan;
      if (route.schedulePayloadSchema == nullptr ||
          !sameOutputSchedulePayloadSchema(
              observerPlan->occurrencePayloadSchema,
              *route.schedulePayloadSchema))
        return invalid("Output schedule and observer occurrence payload schemas "
                       "are incompatible.");
    }
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVObserverOutputEvaluationService::useFieldEvaluationService(
    WVFieldEvaluationService &fieldEvaluationService) {
  const bool compatible =
      impl_->fields != nullptr &&
      impl_->fields->isCompatibleWith(fieldEvaluationService) &&
      fieldEvaluationService.isCompatibleWith(impl_->stateLayout);
  if (!compatible)
    return invalid("Borrowed field-evaluation service uses an incompatible "
                   "resolved transform layout.");
  impl_->ownedFields.reset();
  impl_->fields = &fieldEvaluationService;
  metrics_.retainedStorageBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus
WVObserverOutputEvaluationService::prepareInitial(const WVState &state) {
  return prepareInitial(WVIntegrationState{state, nullptr, 0, nullptr, 0});
}

WVKernelStatus WVObserverOutputEvaluationService::prepareInitial(
    const WVIntegrationState &state) {
  const auto started = std::chrono::steady_clock::now();
  impl_->preparedOccurrences.clear();
  impl_->eventFieldBatchEntries.clear();
  impl_->preparedOutputEvent = false;
  ++impl_->preparationGeneration;
  if (impl_->preparationGeneration == 0)
    ++impl_->preparationGeneration;
  impl_->updateOccurrenceMetrics(metrics_);
  impl_->preparedEventOrdinal = 0;
  impl_->preparedScheduledTime = state.waveVortex.t;
  const auto status = impl_->evaluate(state, true, false, metrics_);
  metrics_.evaluationSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  metrics_.retainedStorageBytes = persistentBytes();
  return status;
}

WVKernelStatus WVObserverOutputEvaluationService::prepare(
    const WVOutputEvent &event) {
  const auto started = std::chrono::steady_clock::now();
  if (!sameTime(event.state.waveVortex.t, event.scheduledTime))
    return invalid("Observation occurrence state is not evaluated at its "
                   "scheduled trigger time.");
  if (impl_->preparedOutputEvent &&
      impl_->preparedEventOrdinal == event.eventOrdinal &&
      impl_->preparedScheduledTime == event.scheduledTime) {
    bool complete = true;
    for (std::size_t routeIndex = 0;
         routeIndex < event.routeCount && complete; ++routeIndex)
      for (std::size_t observerIndex = 0;
           observerIndex < event.routes[routeIndex].observerCount;
           ++observerIndex)
        if (impl_->occurrence(
                event.routes[routeIndex],
                event.routes[routeIndex].observers[observerIndex]) == nullptr) {
          complete = false;
          break;
        }
    if (complete)
      return WVKernelStatus::ok();
  }

  impl_->preparedOccurrences.clear();
  impl_->eventFieldBatchEntries.clear();
  impl_->preparedOutputEvent = false;
  ++impl_->preparationGeneration;
  if (impl_->preparationGeneration == 0)
    ++impl_->preparationGeneration;
  impl_->updateOccurrenceMetrics(metrics_);
  bool needsMoving = event.routes == nullptr;
  for (std::size_t route = 0; route < event.routeCount && !needsMoving; ++route)
    for (std::size_t observer = 0; observer < event.routes[route].observerCount;
         ++observer) {
      const auto ordinal =
          event.routes[route].observers[observer].observerOrdinal;
      if (ordinal < impl_->observerBindings.size() &&
          impl_->observerBindings[ordinal].hasMovingChannels) {
        needsMoving = true;
        break;
      }
    }
  impl_->preparedEventOrdinal = event.eventOrdinal;
  impl_->preparedScheduledTime = event.scheduledTime;
  auto status = impl_->evaluate(event.state, false, needsMoving, metrics_);
  if (status) {
    std::size_t requestedOccurrenceCount = 0;
    for (std::size_t routeIndex = 0; routeIndex < event.routeCount;
         ++routeIndex)
      requestedOccurrenceCount += event.routes[routeIndex].observerCount;
    impl_->preparedOccurrences.reserve(requestedOccurrenceCount);

    for (std::size_t routeIndex = 0; routeIndex < event.routeCount && status;
         ++routeIndex) {
      const auto &route = event.routes[routeIndex];
      if (route.schedulePayloadSchema == nullptr ||
          route.schedulePayload == nullptr ||
          route.proposedScheduleCursor == nullptr) {
        status = invalid("An output route has no resolved occurrence payload.");
        break;
      }
      for (std::size_t observerIndex = 0; observerIndex < route.observerCount;
           ++observerIndex) {
        const auto &observerView = route.observers[observerIndex];
        if (impl_->occurrence(route, observerView) != nullptr) {
          ++metrics_.occurrenceReuseCount;
          continue;
        }
        if (observerView.observerOrdinal >= impl_->observerBindings.size()) {
          status = invalid("An output route has an invalid observer ordinal.");
          break;
        }
        const auto &binding =
            impl_->observerBindings[observerView.observerOrdinal];
        if (binding.record != observerView.record || binding.plan == nullptr ||
            binding.implementation == nullptr ||
            binding.eventFieldPlan == nullptr) {
          status =
              invalid("An output route has an incompatible observer binding.");
          break;
        }
        if (binding.plan->occurrencePayloadSchema.fingerprint() !=
                route.schedulePayloadSchema->fingerprint() ||
            route.schedulePayload->schemaFingerprint() !=
                route.schedulePayloadSchema->fingerprint()) {
          status = invalid("An output occurrence payload is incompatible with "
                           "its observer plan.");
          break;
        }

        Impl::PreparedOccurrence occurrence;
        occurrence.proposedCursor = route.proposedScheduleCursor;
        occurrence.payload = *route.schedulePayload;
        occurrence.observer = binding.record;
        occurrence.plan = binding.plan;
        occurrence.eventFieldPlan = binding.eventFieldPlan;
        occurrence.outputs = binding.outputs;
        occurrence.implementation = binding.implementation;
        occurrence.observerStateViews.reserve(
            binding.plan->occurrenceStateBlocks.size());
        for (const auto &stateBlock : binding.plan->occurrenceStateBlocks) {
          if (stateBlock.resolvedAdditionalStateBlockIndex >=
                  event.state.additionalBlockCount ||
              event.state.additionalBlocks == nullptr) {
            status = invalid(
                "Resolved occurrence observer state is unavailable.");
            break;
          }
          occurrence.observerStateViews.push_back(
              event.state.additionalBlocks
                  [stateBlock.resolvedAdditionalStateBlockIndex]);
        }
        if (!status)
          break;
        WVObserverOccurrencePreparationContext context;
        context.scheduledTime = event.scheduledTime;
        context.scheduleOrdinal = route.scheduleOrdinal;
        context.payloadSchema = route.schedulePayloadSchema;
        context.payload = route.schedulePayload;
        context.observerStateBlocks = occurrence.observerStateViews.data();
        context.observerStateBlockCount =
            occurrence.observerStateViews.size();
        status = occurrence.implementation->prepareOccurrence(
            *binding.record, *binding.plan, context, occurrence.workspace);
        if (!status)
          break;
        status =
            validateOccurrenceWorkspace(*binding.plan, occurrence.workspace);
        if (!status)
          break;

        const auto &eventPlan = *binding.eventFieldPlan;
        if (eventPlan.outputCount() != 0) {
          if (eventPlan.positionSetCount() >
              occurrence.workspace.positionSets.size()) {
            status =
                invalid("An occurrence field plan has unavailable position "
                        "sets.");
            break;
          }
          occurrence.positionViews.reserve(eventPlan.positionSetCount());
          for (std::size_t slot = 0; slot < eventPlan.positionSetCount();
               ++slot) {
            const auto &set = occurrence.workspace.positionSets[slot];
            occurrence.positionViews.push_back(
                {set.x.data(), set.y.data(),
                 set.z.empty() ? nullptr : set.z.data(), set.x.size(),
                 set.extents.data(), set.extents.size()});
          }
          status = impl_->fields->prepareEventGeometry(
              eventPlan, occurrence.positionViews.data(),
              occurrence.positionViews.size(), occurrence.fieldGeometry);
          if (!status)
            break;
          occurrence.fieldStorage.resize(eventPlan.outputCount());
          occurrence.fieldViews.resize(eventPlan.outputCount());
          for (std::size_t outputIndex = 0;
               outputIndex < eventPlan.outputCount(); ++outputIndex) {
            const auto count =
                occurrence.fieldGeometry.outputs()[outputIndex].elementCount;
            occurrence.fieldStorage[outputIndex].resize(count);
            occurrence.fieldViews[outputIndex] = {
                occurrence.fieldStorage[outputIndex].data(), count};
          }
        }

        occurrence.identity.observerOrdinal = observerView.observerOrdinal;
        occurrence.identity.preparationOwner = impl_.get();
        occurrence.identity.preparationGeneration =
            impl_->preparationGeneration;
        occurrence.identity.preparedOccurrenceSlot =
            impl_->preparedOccurrences.size();
        occurrence.identity.resolvedObserverRecord = binding.record;
        occurrence.identity.logicalScheduleRecord =
            route.semanticScheduleRecord;
        occurrence.identity.schedulePayloadSchema =
            route.schedulePayloadSchema;
        occurrence.identity.proposedScheduleCursor =
            route.proposedScheduleCursor;
        occurrence.identity.resolvedSchedulePayload =
            route.schedulePayload;
        occurrence.identity.semanticScheduleOrdinal =
            route.semanticScheduleOrdinal;
        occurrence.identity.scheduleOrdinal = route.scheduleOrdinal;
        occurrence.identity.scheduledTime = event.scheduledTime;
        occurrence.identity.scheduleCursorIdentity =
            route.scheduleCursorIdentity;
        occurrence.identity.payloadFingerprint =
            route.schedulePayload->valueFingerprint();
        occurrence.identity.geometryFingerprint =
            occurrence.workspace.geometryFingerprint();
        occurrence.identity.fieldPlanFingerprint =
            eventPlan.fieldPlanFingerprint();
        impl_->preparedOccurrences.push_back(std::move(occurrence));
        ++metrics_.occurrencePreparationCount;
        impl_->updateOccurrenceMetrics(metrics_);
      }
    }
  }
  if (status) {
    impl_->eventFieldBatchEntries.clear();
    impl_->eventFieldBatchEntries.reserve(impl_->preparedOccurrences.size());
    for (auto &occurrence : impl_->preparedOccurrences)
      if (occurrence.eventFieldPlan != nullptr &&
          occurrence.eventFieldPlan->outputCount() != 0)
        impl_->eventFieldBatchEntries.push_back(
            {occurrence.eventFieldPlan, &occurrence.fieldGeometry,
             occurrence.fieldViews.data(), occurrence.fieldViews.size()});
    impl_->updateOccurrenceMetrics(metrics_);
    if (!impl_->eventFieldBatchEntries.empty()) {
      status = impl_->fields->evaluateEventBatch(
          event.state, impl_->eventFieldBatchEntries.data(),
          impl_->eventFieldBatchEntries.size());
      if (status)
        ++metrics_.fieldEvaluationCount;
    }
    impl_->eventFieldBatchEntries.clear();
    impl_->updateOccurrenceMetrics(metrics_);
  }
  if (status)
    impl_->preparedOutputEvent = true;
  else
    impl_->prepared = false;
  metrics_.evaluationSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  metrics_.retainedStorageBytes = persistentBytes();
  return status;
}

void WVObserverOutputEvaluationService::complete(
    const WVOutputEvent &event) noexcept {
  if (!impl_->preparedOutputEvent ||
      impl_->preparedEventOrdinal != event.eventOrdinal ||
      impl_->preparedScheduledTime != event.scheduledTime)
    return;
  impl_->preparedOccurrences.clear();
  impl_->eventFieldBatchEntries.clear();
  impl_->preparedOutputEvent = false;
  impl_->updateOccurrenceMetrics(metrics_);
  metrics_.retainedStorageBytes = persistentBytes();
}

std::size_t WVObserverOutputEvaluationService::occurrenceWorkspaceRetainedBytes()
    const noexcept {
  return metrics_.occurrenceWorkspaceRetainedBytes;
}

std::size_t WVObserverOutputEvaluationService::occurrenceWorkspaceLiveBytes()
    const noexcept {
  return metrics_.occurrenceWorkspaceLiveBytes;
}

std::size_t WVObserverOutputEvaluationService::persistentBytes() const noexcept {
  if (!impl_)
    return sizeof(*this);
  std::size_t bytes =
      sizeof(*this) + sizeof(Impl) +
      impl_->stateLayout.persistentBytes() +
      (impl_->ownedFields ? impl_->ownedFields->persistentBytes() : 0) +
      impl_->initialFieldPlan.persistentBytes() -
          sizeof(impl_->initialFieldPlan) +
      impl_->timeSeriesFieldPlan.persistentBytes() -
          sizeof(impl_->timeSeriesFieldPlan) +
      impl_->movingFieldPlan.persistentBytes() -
          sizeof(impl_->movingFieldPlan) +
      impl_->initialFieldStorage.capacity() *
          sizeof(std::vector<double>) +
      impl_->timeSeriesFieldStorage.capacity() *
          sizeof(std::vector<double>) +
      impl_->movingFieldStorage.capacity() * sizeof(std::vector<double>) +
      impl_->initialFieldViews.capacity() * sizeof(WVFieldOutputView) +
      impl_->timeSeriesFieldViews.capacity() * sizeof(WVFieldOutputView) +
      impl_->movingFieldViews.capacity() * sizeof(WVFieldOutputView) +
      impl_->movingCoordinates.capacity() * sizeof(Impl::MovingCoordinates) +
      impl_->observerPlans.capacity() * sizeof(WVObserverOutputPlan) +
      impl_->eventFieldPlans.capacity() *
          sizeof(WVEventFieldEvaluationPlan) +
      impl_->observerOutputs.capacity() *
          sizeof(std::vector<Impl::Output>) +
      impl_->observerBindings.capacity() * sizeof(Impl::ObserverBinding) +
      impl_->preparedOccurrences.capacity() *
          sizeof(Impl::PreparedOccurrence) +
      impl_->eventFieldBatchEntries.capacity() *
          sizeof(WVEventFieldEvaluationBatchEntry);
  for (const auto &storage : impl_->initialFieldStorage)
    bytes += storage.capacity() * sizeof(double);
  for (const auto &storage : impl_->timeSeriesFieldStorage)
    bytes += storage.capacity() * sizeof(double);
  for (const auto &storage : impl_->movingFieldStorage)
    bytes += storage.capacity() * sizeof(double);
  bytes += (impl_->movingX.capacity() + impl_->movingY.capacity() +
            impl_->movingZ.capacity()) *
           sizeof(double);
  for (const auto &coordinates : impl_->movingCoordinates)
    bytes += coordinates.fixedZ.capacity() * sizeof(double);
  for (const auto &plan : impl_->observerPlans)
    bytes += observerOutputPlanRetainedBytes(plan);
  for (const auto &plan : impl_->eventFieldPlans)
    bytes += plan.persistentBytes() - sizeof(plan);
  for (const auto &outputs : impl_->observerOutputs) {
    bytes += outputs.capacity() * sizeof(Impl::Output);
    for (const auto &output : outputs)
      bytes += output.extents.capacity() * sizeof(std::size_t) +
               output.affineStorage.capacity() * sizeof(double);
  }
  for (const auto &occurrence : impl_->preparedOccurrences)
    bytes += occurrence.retainedBytes() - sizeof(Impl::PreparedOccurrence);
  return bytes;
}

} // namespace wavevortex::runtime
