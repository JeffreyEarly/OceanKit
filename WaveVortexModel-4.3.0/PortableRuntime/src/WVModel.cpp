#include "WaveVortexRuntime/WVModel.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVModelInternalAccess.hpp"
#include "WVModelTransformAdapters.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <new>
#include <set>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

std::size_t outputDriverMetricsDynamicBytes(
    const WVOutputDriverMetrics &metrics) noexcept {
  std::size_t bytes =
      metrics.files.capacity() * sizeof(WVOutputFileMetrics);
  for (const auto &file : metrics.files) {
    bytes += file.fileIdentifier.capacity() + file.destination.capacity() +
             file.groups.capacity() * sizeof(WVOutputGroupMetrics);
    for (const auto &group : file.groups)
      bytes += group.fileIdentifier.capacity() +
               group.groupIdentifier.capacity();
  }
  return bytes;
}

#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
WVKernelStatus fromCheckpoint(const WVCheckpointStatus &status) {
  if (status)
    return WVKernelStatus::ok();
  return {status.code == WVCheckpointStatusCode::unsupportedObserver
              ? WVKernelStatusCode::unsupportedOperation
              : WVKernelStatusCode::invalidConfiguration,
          status.message +
              (status.location.empty() ? "" : " [" + status.location + "]")};
}

const std::string *remappedDestination(
    const std::vector<WVModelOutputDestination> &destinations,
    const std::string &identifier) {
  const auto found =
      std::find_if(destinations.begin(), destinations.end(),
                   [&](const auto &candidate) {
                     return candidate.fileIdentifier == identifier;
                   });
  return found == destinations.end() ? nullptr : &found->path;
}

bool sameObservationSchema(const WVObservationSchema &left,
                           const WVObservationSchema &right) {
  std::vector<std::uint8_t> leftBytes;
  std::vector<std::uint8_t> rightBytes;
  return encodeObservationSchemaManifest(left, leftBytes) &&
         encodeObservationSchemaManifest(right, rightBytes) &&
         leftBytes == rightBytes;
}

WVKernelStatus compileOutputConfiguration(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVModelOutputNetCDFInspection &inspection,
    const WVModelOutputRequest &request,
    WVModelOutputConfiguration &configuration) {
  if (!std::isfinite(request.finalTime) ||
      request.finalTime < inspection.latestRestart.t)
    return invalid("The output request final time must be finite and cannot "
                   "precede the selected restart state.");
  std::set<std::string> requestedIdentifiers;
  for (const auto &destination : request.destinations) {
    if (destination.fileIdentifier.empty() || destination.path.empty() ||
        !requestedIdentifiers.insert(destination.fileIdentifier).second)
      return invalid("Output destination remapping requires unique nonempty "
                     "file identifiers and paths.");
  }
  if ((request.policy == WVModelOutputPolicy::create ||
       request.policy == WVModelOutputPolicy::replace ||
       !request.destinations.empty()) &&
      request.destinations.size() != inspection.observerRecord.outputFiles.size())
    return invalid("Output destination remapping requires one destination for "
                   "every output-file identifier.");

  auto observerRecord = inspection.observerRecord;
  std::set<std::string> normalizedDestinations;
  for (auto &record : observerRecord.outputFiles) {
    const auto *mapped = remappedDestination(request.destinations,
                                             record.identifier);
    if (mapped == nullptr && !request.destinations.empty())
      return invalid("The output destination map is incomplete.");
    const std::string &path = mapped == nullptr ? record.destination : *mapped;
    try {
      record.destination =
          std::filesystem::absolute(std::filesystem::path(path))
              .lexically_normal()
              .string();
    } catch (const std::exception &error) {
      return invalid(std::string("The output destination cannot be resolved: ") +
                     error.what());
    }
    if (!normalizedDestinations.insert(record.destination).second)
      return invalid("Output destinations must not alias each other.");
    if (request.policy != WVModelOutputPolicy::append) {
      for (const auto &source : inspection.paths) {
        const auto normalizedSource =
            std::filesystem::absolute(std::filesystem::path(source))
                .lexically_normal()
                .string();
        if (record.destination == normalizedSource)
          return invalid("Create and replace destinations must not alias a "
                         "source model-output file.");
      }
    }
  }
  if (requestedIdentifiers.size() != request.destinations.size())
    return invalid("The output destination map is invalid.");
  for (const auto &identifier : requestedIdentifiers) {
    const auto found = std::find_if(
        observerRecord.outputFiles.begin(), observerRecord.outputFiles.end(),
        [&](const auto &file) {
          return file.identifier == identifier;
        });
    if (found == observerRecord.outputFiles.end())
      return invalid("The output destination map contains an unknown file "
                     "identifier.");
  }
  return WVModelOutputConfiguration::compile(
      std::move(observerRecord), inspection.observationSchemas,
      inspection.scheduleContinuations, request.policy, std::move(catalog),
      inspection.latestRestart.t, request.finalTime, configuration,
      detail::legacyModelOutputPlanningConfiguration(
          inspection.latestRestart),
      inspection.isDynamicsLinear,
      &inspection.latestRestart.stateDescription);
}
#endif

WVKernelStatus copyAdditionalState(const WVAdditionalStateStorage &source,
                                   WVAdditionalStateStorage &destination) {
  const auto *sourceBlocks = source.constBlocks();
  for (std::size_t destinationIndex = 0;
       destinationIndex < destination.blockCount(); ++destinationIndex) {
    auto &target = destination.mutableBlocks()[destinationIndex];
    const WVAdditionalStateBlockConstView *matching = nullptr;
    for (std::size_t sourceIndex = 0; sourceIndex < source.blockCount();
         ++sourceIndex) {
      if (sourceBlocks[sourceIndex].layout->identifier ==
          target.layout->identifier) {
        matching = sourceBlocks + sourceIndex;
        break;
      }
    }
    if (matching == nullptr ||
        matching->layout->scalarType != target.layout->scalarType ||
        matching->layout->elementCount != target.layout->elementCount)
      return invalid("Restored observer state does not match the model state "
                     "layout.");
    if (target.realData != nullptr)
      std::copy_n(matching->realData, target.layout->elementCount,
                  target.realData);
    else
      std::copy_n(matching->complexData, target.layout->elementCount,
                  target.complexData);
  }
  return WVKernelStatus::ok();
}

} // namespace

WVKernelStatus WVModelState::create(
    WVCheckpoint checkpoint, const WVIntegrationStateLayout &layout,
    WVModelState &state,
    const WVAdditionalStateStorage *restoredAdditionalState) {
  auto checkpointStatus =
      detail::validateModelCheckpointState(checkpoint, layout);
  if (!checkpointStatus)
    return checkpointStatus;
  try {
    WVModelState candidate;
    auto status = candidate.additionalState_.initialize(layout);
    if (!status)
      return status;
    if (restoredAdditionalState != nullptr) {
      status = copyAdditionalState(*restoredAdditionalState,
                                   candidate.additionalState_);
      if (!status)
        return status;
    }
    candidate.checkpoint_ = std::move(checkpoint);
    candidate.rebuildCoefficientViews(layout);
    state = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel state allocation failed."};
  }
}

void WVModelState::rebuildCoefficientViews(
    const WVIntegrationStateLayout &layout) {
  layout_ = &layout;
  mutableCoefficientViews_.clear();
  constCoefficientViews_.clear();
  mutableCoefficientViews_.reserve(layout.coefficientFamilyCount());
  constCoefficientViews_.reserve(layout.coefficientFamilyCount());
  for (std::size_t family = 0; family < layout.coefficientFamilyCount();
       ++family) {
    auto *values = detail::modelCheckpointCoefficientData(
        checkpoint_, layout, family);
    mutableCoefficientViews_.push_back(
        {&layout.coefficientFamilies()[family], values});
    constCoefficientViews_.push_back(
        {&layout.coefficientFamilies()[family], values});
  }
}

void WVModelState::setTimes(double t, double t0) noexcept {
  detail::setModelCheckpointTimes(checkpoint_, t, t0);
}

WVMutableIntegrationState WVModelState::mutableView() noexcept {
  auto legacy = detail::modelCheckpointLegacyView(checkpoint_);
  return {legacy,
          additionalState_.mutableBlocks(), additionalState_.blockCount(),
          mutableCoefficientViews_.data(), mutableCoefficientViews_.size()};
}

WVIntegrationState WVModelState::constView() {
  auto mutableState = mutableView();
  return integrationConstView(mutableState, constCoefficientViews_,
                              constViews_);
}

std::size_t WVModelState::persistentBytes() const noexcept {
  return sizeof(*this) + checkpointRetainedBytes(checkpoint_) -
         sizeof(checkpoint_) +
         additionalState_.capacityBytes() +
         mutableCoefficientViews_.capacity() *
             sizeof(WVCoefficientFamilyView) +
         constCoefficientViews_.capacity() *
             sizeof(WVCoefficientFamilyConstView) +
         constViews_.capacity() * sizeof(WVAdditionalStateBlockConstView);
}

class WVModel::Impl final {
public:
  WVKernelStatus configureIntegrator(
      const WVModelIntegratorConfiguration &configuration);

  std::unique_ptr<detail::WVResolvedModelSystem> resolvedSystem;
  std::shared_ptr<const WVExtensionCatalog> catalog;
  std::unique_ptr<WVTimeIntegrator> integrator;
  WVModelIntegratorKind integratorKind = WVModelIntegratorKind::fixedRK4;
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  std::unique_ptr<WVModelOutputConfiguration> outputConfiguration;
  std::unique_ptr<WVObserverOutputEvaluationService> outputEvaluation;
  WVModelOutputNetCDFSink outputSink;
  WVModelOutputNetCDFMetrics outputMetrics;
  std::unique_ptr<WVOutputDriver> outputDriver;
  const WVOutputPlan *outputDriverPlan = nullptr;
  WVOutputSink *outputDriverSink = nullptr;
  double outputFinalTime = 0.0;
  bool outputOpen = false;
#endif
  WVOutputDriverMetrics outputDriverMetrics;
};

WVModel::WVModel() : impl_(std::make_unique<Impl>()) {}
WVModel::~WVModel() = default;
WVModel::WVModel(WVModel &&) noexcept = default;
WVModel &WVModel::operator=(WVModel &&) noexcept = default;

WVKernelStatus WVModel::Impl::configureIntegrator(
    const WVModelIntegratorConfiguration &configuration) {
  try {
    if (configuration.kind == WVModelIntegratorKind::adaptiveRK23)
      integrator =
          std::make_unique<WVAdaptiveRK23>(resolvedSystem->integrationSystem(),
                                           configuration.adaptive);
    else if (configuration.kind == WVModelIntegratorKind::adaptiveRK45)
      integrator = std::make_unique<WVAdaptiveRK45>(
          resolvedSystem->integrationSystem(), configuration.adaptiveRK45);
    else if (configuration.kind == WVModelIntegratorKind::adaptiveRK78)
      integrator = std::make_unique<WVAdaptiveRK78>(
          resolvedSystem->integrationSystem(), configuration.adaptiveRK78);
    else
      integrator =
          std::make_unique<WVFixedStepRK4>(resolvedSystem->integrationSystem(),
                                           configuration.fixed);
    integratorKind = configuration.kind;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel integrator allocation failed."};
  }
}

WVKernelStatus WVModel::create(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &forcingSchedule,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model) {
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  try {
    auto candidate = std::make_unique<Impl>();
    auto status = detail::createConstantStratificationModelSystem(
        configuration, forcingSchedule, nullptr, catalog, std::move(engine),
        candidate->resolvedSystem);
    if (!status)
      return status;
    status = candidate->configureIntegrator(integratorConfiguration);
    if (!status)
      return status;
    candidate->catalog = std::move(catalog);
    model.impl_ = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel allocation failed."};
  }
}

WVKernelStatus WVModel::createFromCheckpoint(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    WVCheckpoint checkpoint, std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model, WVModelState &state) {
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  try {
    auto candidate = std::make_unique<Impl>();
    const auto inspection = detail::modelCheckpointInspection(checkpoint);
    auto status = detail::createPersistedModelSystem(
        inspection, checkpoint.forcingSchedule, nullptr, catalog,
        std::move(engine), candidate->resolvedSystem);
    if (!status)
      return status;
    status = candidate->configureIntegrator(integratorConfiguration);
    if (!status)
      return status;
    WVModelState candidateState;
    status = WVModelState::create(std::move(checkpoint),
                                  candidate->resolvedSystem
                                      ->integrationSystem()
                                      .stateLayout(),
                                  candidateState);
    if (!status)
      return status;
    candidate->catalog = std::move(catalog);
    model.impl_ = std::move(candidate);
    state = std::move(candidateState);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel allocation failed."};
  }
}

WVKernelStatus WVModel::validateCheckpointForcingSchedule(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVExtensionCatalog &catalog) {
  return detail::validatePersistedModelForcingSchedule(inspection, schedule,
                                                        catalog);
}

WVKernelStatus WVModel::create(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &forcingSchedule,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model) {
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  try {
    auto candidate = std::make_unique<Impl>();
    auto status = detail::createBarotropicQGModelSystem(
        configuration, forcingSchedule, nullptr, catalog, std::move(engine),
        candidate->resolvedSystem);
    if (!status)
      return status;
    status = candidate->configureIntegrator(integratorConfiguration);
    if (!status)
      return status;
    candidate->catalog = std::move(catalog);
    model.impl_ = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel allocation failed."};
  }
}

WVKernelStatus WVModel::createFromModelOutputFiles(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const std::vector<std::string> &paths,
    const WVModelOutputRequest &outputRequest,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model, WVModelState &state) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  if (paths.empty())
    return invalid("At least one model-output NetCDF path is required.");
  WVModelOutputNetCDFInspection inspection;
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  auto checkpointStatus = WVModelOutputNetCDFSink::inspect(paths, *catalog,
                                                           inspection);
  if (!checkpointStatus)
    return fromCheckpoint(checkpointStatus);

  return createFromModelOutputInspection(
      std::move(catalog), std::move(inspection), outputRequest, std::move(engine),
      integratorConfiguration, model, state);
#else
  (void)paths;
  (void)outputRequest;
  (void)engine;
  (void)integratorConfiguration;
  (void)model;
  (void)state;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include model-output persistence."};
#endif
}

WVKernelStatus WVModel::createFromModelOutputInspection(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    WVModelOutputNetCDFInspection inspection,
    const WVModelOutputRequest &outputRequest,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model, WVModelState &state) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT

  WVModelOutputConfiguration outputConfiguration;
  auto status = prepareModelOutput(catalog, inspection, outputRequest,
                                   outputConfiguration);
  if (!status)
    return status;
  return createFromModelOutputInspection(
      std::move(catalog), std::move(inspection), std::move(outputConfiguration),
      std::move(engine), integratorConfiguration, model, state);
#else
  (void)inspection;
  (void)outputRequest;
  (void)engine;
  (void)integratorConfiguration;
  (void)model;
  (void)state;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include model-output persistence."};
#endif
}

WVKernelStatus WVModel::prepareModelOutput(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVModelOutputNetCDFInspection &inspection,
    const WVModelOutputRequest &outputRequest,
    WVModelOutputConfiguration &configuration) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  return compileOutputConfiguration(std::move(catalog), inspection, outputRequest,
                                    configuration);
#else
  (void)inspection;
  (void)outputRequest;
  (void)configuration;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include model-output persistence."};
#endif
}

WVKernelStatus WVModel::createFromModelOutputInspection(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    WVModelOutputNetCDFInspection inspection,
    WVModelOutputConfiguration outputConfiguration,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model, WVModelState &state) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT

  if (!catalog || outputConfiguration.catalog() != catalog)
    return invalid("Prepared output and model must share one extension catalog.");
  auto forcingSchedule = inspection.latestRestart.forcingSchedule;
  if (inspection.isDynamicsLinear)
    forcingSchedule.entries.clear();
  const auto &descriptor = outputConfiguration.descriptor();

  WVModel candidate;
  auto candidateImpl = std::make_unique<Impl>();
  auto status = detail::createPersistedModelSystem(
      inspection.latestRestart, forcingSchedule, &descriptor, catalog,
      std::move(engine), candidateImpl->resolvedSystem);
  if (status)
    status = candidateImpl->configureIntegrator(integratorConfiguration);
  if (status) {
    candidateImpl->catalog = catalog;
    candidate.impl_ = std::move(candidateImpl);
  }
  if (!status)
    return status;

  WVCheckpoint restoredCheckpoint;
  WVAdditionalStateStorage restoredAdditionalState;
  auto checkpointStatus = WVModelOutputNetCDFSink::restoreState(
      inspection, *catalog, candidate.stateLayout(), restoredCheckpoint,
      restoredAdditionalState);
  if (!checkpointStatus)
    return fromCheckpoint(checkpointStatus);

  WVModelState candidateState;
  status = WVModelState::create(std::move(restoredCheckpoint),
                                candidate.stateLayout(), candidateState,
                                &restoredAdditionalState);
  if (!status)
    return status;

  status = candidate.openOutput(candidateState, std::move(outputConfiguration),
                                inspection.isDynamicsLinear);
  if (!status)
    return status;
  state = std::move(candidateState);
  model = std::move(candidate);
  return WVKernelStatus::ok();
#else
  (void)inspection;
  (void)outputConfiguration;
  (void)engine;
  (void)integratorConfiguration;
  (void)model;
  (void)state;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include model-output persistence."};
#endif
}

WVKernelStatus WVModel::create(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &forcingSchedule,
    const WVPortableObserverDescriptor &observerDescriptor,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model) {
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  if (observerDescriptor.catalog() != catalog)
    return invalid("WVModel and its observer descriptor require the same "
                   "extension catalog.");
  try {
    auto candidate = std::make_unique<Impl>();
    auto status = detail::createConstantStratificationModelSystem(
        configuration, forcingSchedule, &observerDescriptor, catalog,
        std::move(engine), candidate->resolvedSystem);
    if (!status)
      return status;
    status = candidate->configureIntegrator(integratorConfiguration);
    if (!status)
      return status;
    candidate->catalog = std::move(catalog);
    model.impl_ = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel allocation failed."};
  }
}

WVKernelStatus WVModel::create(
    std::shared_ptr<const WVExtensionCatalog> catalog,
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &forcingSchedule,
    const WVPortableObserverDescriptor &observerDescriptor,
    std::unique_ptr<WVFFTEngine> engine,
    const WVModelIntegratorConfiguration &integratorConfiguration,
    WVModel &model) {
  if (!catalog)
    return invalid("WVModel requires an extension catalog.");
  if (observerDescriptor.catalog() != catalog)
    return invalid("WVModel and its observer descriptor require the same "
                   "extension catalog.");
  try {
    auto candidate = std::make_unique<Impl>();
    auto status = detail::createBarotropicQGModelSystem(
        configuration, forcingSchedule, &observerDescriptor, catalog,
        std::move(engine), candidate->resolvedSystem);
    if (!status)
      return status;
    status = candidate->configureIntegrator(integratorConfiguration);
    if (!status)
      return status;
    candidate->catalog = std::move(catalog);
    model.impl_ = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "WVModel allocation failed."};
  }
}

WVKernelStatus WVModel::prepareStateAfterRestart(WVModelState &state) {
  auto view = state.mutableView();
  return impl_->integrator->prepareStateAfterRestart(view);
}

WVKernelStatus WVModel::initializeObserverState(WVModelState &state) {
  auto view = state.mutableView();
  return impl_->resolvedSystem->initializeObserverState(view);
}

WVKernelStatus WVModel::openOutput(
    WVModelState &state, WVModelOutputConfiguration outputConfiguration,
    bool isDynamicsLinear) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  if (impl_->outputOpen)
    return invalid("WVModel already has an open output graph.");
  if (outputConfiguration.catalog() != impl_->catalog)
    return invalid("Prepared output and model must share one extension catalog.");
  auto modelState = state.constView();
  auto status = validateIntegrationState(stateLayout(), modelState);
  if (!status)
    return status;
  auto *fieldEvaluationService =
      impl_->resolvedSystem->integrationSystem().fieldEvaluationService();
  if (fieldEvaluationService == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "The selected transform does not provide field evaluation for "
            "model output."};
  std::unique_ptr<WVObserverOutputEvaluationService> outputEvaluation;
  status = WVObserverOutputEvaluationService::create(
      isDynamicsLinear, outputConfiguration.descriptor(),
      *fieldEvaluationService, outputEvaluation);
  if (!status)
    return status;
  for (const auto &declared : outputConfiguration.observationSchemas()) {
    const auto &observers = outputConfiguration.descriptor().observers();
    const auto observer = std::find_if(
        observers.begin(), observers.end(), [&](const auto &candidate) {
          return candidate.identifier == declared.observerIdentifier;
        });
    if (observer == observers.end())
      return invalid("A declared observation schema references an unknown "
                     "compiled observer.");
    WVObservationSchema resolvedSchema;
    status = outputEvaluation->observationSchema(*observer, resolvedSchema);
    if (!status)
      return status;
    if (!sameObservationSchema(declared.schema, resolvedSchema))
      return invalid("The resolved observer schema differs from the canonical "
                     "schema restored from NetCDF.");
  }
  WVModelOutputNetCDFSink outputSink;
  WVModelOutputNetCDFConfiguration sinkConfiguration{
      impl_->catalog, state.checkpoint(), isDynamicsLinear};
  const auto checkpointStatus = outputConfiguration.openNetCDFSink(
      sinkConfiguration, stateLayout(), outputEvaluation.get(), outputSink);
  if (!checkpointStatus)
    return fromCheckpoint(checkpointStatus);
  impl_->outputFinalTime = outputConfiguration.plan().finalTime();
  impl_->outputSink = std::move(outputSink);
  impl_->outputEvaluation = std::move(outputEvaluation);
  impl_->outputConfiguration =
      std::make_unique<WVModelOutputConfiguration>(
          std::move(outputConfiguration));
  impl_->outputOpen = true;
  return WVKernelStatus::ok();
#else
  (void)state;
  (void)outputConfiguration;
  (void)isDynamicsLinear;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include model-output persistence."};
#endif
}

WVKernelStatus WVModel::evaluateRightHandSide(
    const WVIntegrationState &state, WVIntegrationFlux &rightHandSide) {
  return impl_->resolvedSystem->integrationSystem().evaluateRightHandSide(
      state, rightHandSide);
}

bool WVModel::supportsFixedTimeStepSelection() const noexcept {
  return impl_->resolvedSystem->integrationSystem()
      .supportsFixedTimeStepSelection();
}

WVKernelStatus WVModel::evaluateFixedTimeStepCandidates(
    WVModelState &state, double cfl,
    WVFixedTimeStepCandidates &candidates) {
  return impl_->resolvedSystem->integrationSystem()
      .evaluateFixedTimeStepCandidates(state.constView(), cfl, candidates);
}

WVKernelStatus WVModel::step(WVModelState &state, double proposedStepSize) {
  auto view = state.mutableView();
  const auto status = impl_->integrator->step(view, proposedStepSize);
  if (status)
    state.setTimes(view.waveVortex.t, view.waveVortex.t0);
  return status;
}

WVKernelStatus WVModel::advanceToTime(WVModelState &state, double finalTime,
                                      double initialStepSize) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  if (impl_->outputOpen) {
    if (finalTime != impl_->outputFinalTime)
      return invalid("The requested final time differs from the compiled "
                     "output graph.");
    return advanceToTime(state, finalTime, initialStepSize,
                         impl_->outputConfiguration->plan(),
                         impl_->outputSink);
  }
#endif
  auto view = state.mutableView();
  const auto status =
      impl_->integrator->advanceToTime(view, finalTime, initialStepSize);
  if (status)
    state.setTimes(view.waveVortex.t, view.waveVortex.t0);
  return status;
}

WVCheckpointStatus WVModel::closeOutput() noexcept {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  if (!impl_->outputOpen)
    return WVCheckpointStatus::ok();
  const auto status = impl_->outputSink.close();
  impl_->outputMetrics = impl_->outputSink.metrics();
  if (status) {
    impl_->outputDriver.reset();
    impl_->outputDriverPlan = nullptr;
    impl_->outputDriverSink = nullptr;
    impl_->outputOpen = false;
  }
  return status;
#else
  return WVCheckpointStatus::ok();
#endif
}

bool WVModel::hasOutput() const noexcept {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  return impl_->outputOpen;
#else
  return false;
#endif
}

WVKernelStatus WVModel::advanceToTime(WVModelState &state, double finalTime,
                                      double initialStepSize,
                                      const WVOutputPlan &plan,
                                      WVOutputSink &sink) {
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  auto view = state.mutableView();
  if (impl_->outputDriver == nullptr) {
    try {
      impl_->outputDriver =
          std::make_unique<WVOutputDriver>(*impl_->integrator, plan);
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "WVModel output-driver allocation failed."};
    }
    impl_->outputDriverPlan = &plan;
    impl_->outputDriverSink = &sink;
  } else if (impl_->outputDriverPlan != &plan ||
             impl_->outputDriverSink != &sink) {
    return invalid("A failed output delivery must be retried with the same "
                   "output plan and sink.");
  }
  const auto status = impl_->outputDriver->advanceToTime(
      view, finalTime, initialStepSize, sink);
  impl_->outputDriverMetrics = impl_->outputDriver->metrics();
  state.setTimes(view.waveVortex.t, view.waveVortex.t0);
  if (status) {
    impl_->outputDriver.reset();
    impl_->outputDriverPlan = nullptr;
    impl_->outputDriverSink = nullptr;
  }
  return status;
#else
  (void)state;
  (void)finalTime;
  (void)initialStepSize;
  (void)plan;
  (void)sink;
  return {WVKernelStatusCode::unsupportedOperation,
          "This adapter does not include output orchestration."};
#endif
}

const WVIntegrationStateLayout &WVModel::stateLayout() const noexcept {
  return impl_->resolvedSystem->integrationSystem().stateLayout();
}

double WVModel::nextStepSize() const noexcept {
  return impl_->integrator->nextStepSize();
}

WVModelIntegratorKind WVModel::integratorKind() const noexcept {
  return impl_->integratorKind;
}

const std::string &WVModel::forcingScheduleIdentifier() const noexcept {
  return impl_->resolvedSystem->forcingScheduleIdentifier();
}

const std::string &WVModel::kernelProviderIdentifier() const noexcept {
  return impl_->resolvedSystem->kernelProviderIdentifier();
}

const std::string &WVModel::kernelProviderLibraryIdentity() const noexcept {
  return impl_->resolvedSystem->kernelProviderLibraryIdentity();
}

WVModelMetrics WVModel::metrics(const WVModelState *state) const noexcept {
  WVModelMetrics result;
  result.modelPersistentBytes =
      sizeof(*this) + sizeof(Impl) +
      outputDriverMetricsDynamicBytes(impl_->outputDriverMetrics);
  result.catalogPersistentBytes =
      impl_->catalog == nullptr ? 0 : impl_->catalog->persistentBytes();
  result.statePersistentBytes = state == nullptr ? 0 : state->persistentBytes();
  result.integrationSystemPersistentBytes =
      impl_->resolvedSystem->integrationSystem().persistentBytes();
  result.integratorPersistentBytes = impl_->integrator->persistentBytes();
  impl_->resolvedSystem->populateMetrics(result);
  if (impl_->integratorKind == WVModelIntegratorKind::adaptiveRK23)
    result.integrator =
        static_cast<const WVAdaptiveRK23 &>(*impl_->integrator).metrics();
  else if (impl_->integratorKind == WVModelIntegratorKind::adaptiveRK45)
    result.integrator =
        static_cast<const WVAdaptiveRK45 &>(*impl_->integrator).metrics();
  else if (impl_->integratorKind == WVModelIntegratorKind::adaptiveRK78)
    result.integrator =
        static_cast<const WVAdaptiveRK78 &>(*impl_->integrator).metrics();
  else
    result.integrator =
        static_cast<const WVFixedStepRK4 &>(*impl_->integrator).metrics();
  result.outputDriver = impl_->outputDriverMetrics;
#if !defined(WV_MODEL_ENABLE_OUTPUT) || WV_MODEL_ENABLE_OUTPUT
  if (impl_->outputEvaluation != nullptr)
    result.outputEvaluation = impl_->outputEvaluation->metrics();
  result.output = impl_->outputOpen ? impl_->outputSink.metrics()
                                    : impl_->outputMetrics;
  if (impl_->outputConfiguration != nullptr)
    result.outputConfigurationPersistentBytes =
        impl_->outputConfiguration->persistentBytes();
  result.outputEvaluationPersistentBytes =
      impl_->outputEvaluation == nullptr
          ? 0
          : impl_->outputEvaluation->persistentBytes();
  result.outputSinkPersistentBytes =
      impl_->outputSink.metrics().retainedStorageBytes;
  result.outputPersistentBytes =
      result.outputConfigurationPersistentBytes +
      result.outputEvaluationPersistentBytes +
      result.outputSinkPersistentBytes;
#endif
  return result;
}

WVIntegrationSystem &
detail::WVModelInternalAccess::integrationSystem(WVModel &model) noexcept {
  return model.impl_->resolvedSystem->integrationSystem();
}

WVTransformConstantStratificationKernel &
detail::WVModelInternalAccess::constantStratificationKernel(
    WVModel &model) noexcept {
  return *model.impl_->resolvedSystem->constantStratificationKernel();
}

WVTimeIntegrator &
detail::WVModelInternalAccess::integrator(WVModel &model) noexcept {
  return *model.impl_->integrator;
}

} // namespace wavevortex::runtime
