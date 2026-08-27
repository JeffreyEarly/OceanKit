#include "WVModelTransformAdapters.hpp"

#include <new>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

class ConstantStratificationModelSystem final
    : public WVResolvedModelSystem {
public:
  explicit ConstantStratificationModelSystem(
      std::unique_ptr<WVConstantStratificationIntegrationSystem> system)
      : system_(std::move(system)) {}

  WVIntegrationSystem &integrationSystem() noexcept override {
    return *system_;
  }
  const WVIntegrationSystem &integrationSystem() const noexcept override {
    return *system_;
  }
  const std::string &forcingScheduleIdentifier() const noexcept override {
    return system_->scheduleIdentifier();
  }
  const std::string &kernelProviderIdentifier() const noexcept override {
    return system_->kernel().engineIdentifier();
  }
  const std::string &kernelProviderLibraryIdentity() const noexcept override {
    return system_->kernel().engineLibraryIdentity();
  }
  WVKernelStatus initializeObserverState(
      WVMutableIntegrationState &state) const override {
    return system_->initializeParticleState(state);
  }
  void populateMetrics(WVModelMetrics &metrics) const noexcept override {
    metrics.kernel = system_->kernelMetrics();
    metrics.forcing = system_->forcingMetrics();
    metrics.integratedObservers = system_->metrics();
  }
  WVTransformConstantStratificationKernel *
  constantStratificationKernel() noexcept override {
    return &system_->kernel();
  }

private:
  std::unique_ptr<WVConstantStratificationIntegrationSystem> system_;
};

class BarotropicQGModelSystem final : public WVResolvedModelSystem {
public:
  explicit BarotropicQGModelSystem(
      std::unique_ptr<WVBarotropicQGIntegrationSystem> system)
      : system_(std::move(system)) {}

  WVIntegrationSystem &integrationSystem() noexcept override {
    return *system_;
  }
  const WVIntegrationSystem &integrationSystem() const noexcept override {
    return *system_;
  }
  const std::string &forcingScheduleIdentifier() const noexcept override {
    return system_->forcingScheduleIdentifier();
  }
  const std::string &kernelProviderIdentifier() const noexcept override {
    return system_->kernel().engineIdentifier();
  }
  const std::string &kernelProviderLibraryIdentity() const noexcept override {
    return system_->kernel().engineLibraryIdentity();
  }
  WVKernelStatus initializeObserverState(
      WVMutableIntegrationState &state) const override {
    return system_->initializeParticleState(state);
  }
  void populateMetrics(WVModelMetrics &metrics) const noexcept override {
    metrics.barotropicQGKernel = system_->kernel().metrics();
    metrics.barotropicQGForcing = system_->forcingMetrics();
    metrics.integratedObservers = system_->metrics();
    const auto &kernel = metrics.barotropicQGKernel;
    metrics.kernel.descriptorBytes = kernel.descriptorBytes;
    metrics.kernel.planBytes = kernel.planBytes;
    metrics.kernel.engineBytes = kernel.engineBytes;
    metrics.kernel.kernelManagementBytes = kernel.kernelManagementBytes;
    metrics.kernel.scratchCapacityBytes = kernel.scratchCapacityBytes;
    metrics.kernel.scratchHighWaterBytes = kernel.scratchHighWaterBytes;
    metrics.kernel.halfSpectrumScratchCapacityBytes =
        kernel.halfSpectrumScratchCapacityBytes;
    metrics.kernel.realScratchCapacityBytes = kernel.realScratchCapacityBytes;
    metrics.kernel.planCount = kernel.planCount;
    metrics.kernel.executionCount = kernel.executionCount;
    metrics.kernel.horizontalExecutionCount =
        kernel.forwardExecutionCount + kernel.inverseExecutionCount;
    metrics.kernel.nonlinearFluxCallCount = kernel.nonlinearFluxCallCount;
    metrics.kernel.bytesCopied = kernel.bytesCopied;
    const auto &forcing = metrics.barotropicQGForcing;
    metrics.forcing.scheduleBytes = forcing.scheduleBytes;
    metrics.forcing.derivedOperatorBytes = forcing.derivedOperatorBytes;
    metrics.forcing.workspaceCapacityBytes = forcing.workspaceCapacityBytes;
    metrics.forcing.evaluationCount = forcing.evaluationCount;
    metrics.forcing.restoredCoefficientCount =
        forcing.restoredCoefficientCount;
    metrics.forcing.resolvedSpatialCount = forcing.resolvedSpatialCount;
    metrics.forcing.resolvedSpectralCount = forcing.resolvedSpectralCount;
    metrics.forcing.resolvedAmplitudeCount = forcing.resolvedAmplitudeCount;
    metrics.forcing.physicalFieldReconstructionCount =
        forcing.physicalFieldReconstructionCount;
    metrics.forcing.physicalFieldReuseCount =
        forcing.physicalFieldReuseCount;
    metrics.forcing.spatialTendencyProjectionCount =
        forcing.spatialTendencyProjectionCount;
    metrics.forcing.stateConstraintElementWrites =
        forcing.stateConstraintElementWrites;
  }

private:
  std::unique_ptr<WVBarotropicQGIntegrationSystem> system_;
};

} // namespace

WVKernelStatus createConstantStratificationModelSystem(
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system) {
  std::unique_ptr<WVConstantStratificationIntegrationSystem> numerical;
  auto status = descriptor == nullptr
                    ? WVConstantStratificationIntegrationSystem::create(
                          configuration, schedule, std::move(catalog),
                          std::move(engine), numerical)
                    : WVConstantStratificationIntegrationSystem::create(
                          configuration, schedule, *descriptor,
                          std::move(catalog), std::move(engine), numerical);
  if (!status)
    return status;
  try {
    system = std::make_unique<ConstantStratificationModelSystem>(
        std::move(numerical));
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the constant-stratification model adapter."};
  }
}

WVKernelStatus createBarotropicQGModelSystem(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system) {
  std::unique_ptr<WVBarotropicQGIntegrationSystem> numerical;
  auto status = descriptor == nullptr
                    ? WVBarotropicQGIntegrationSystem::create(
                          configuration, schedule, std::move(catalog),
                          std::move(engine), numerical)
                    : WVBarotropicQGIntegrationSystem::create(
                          configuration, schedule, *descriptor,
                          std::move(catalog), std::move(engine), numerical);
  if (!status)
    return status;
  try {
    system =
        std::make_unique<BarotropicQGModelSystem>(std::move(numerical));
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Barotropic QG model adapter."};
  }
}

WVKernelStatus createPersistedModelSystem(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system) {
  if (inspection.transformKind == WVPersistedTransformKind::barotropicQG)
    return createBarotropicQGModelSystem(
        inspection.barotropicQGConfiguration, schedule, descriptor,
        std::move(catalog), std::move(engine), system);
  return createConstantStratificationModelSystem(
      inspection.configuration, schedule, descriptor, std::move(catalog),
      std::move(engine), system);
}

WVKernelStatus validatePersistedModelForcingSchedule(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVExtensionCatalog &catalog) {
  if (inspection.transformKind == WVPersistedTransformKind::barotropicQG)
    return WVBarotropicQGForcingEngine::validateSchedule(
        inspection.barotropicQGConfiguration, schedule,
        inspection.coefficientShape.elementCount(), catalog);
  return WVConstantStratificationForcingEngine::validateSchedule(
      inspection.configuration, schedule, inspection.coefficientShape,
      catalog);
}

WVCheckpointInspection modelCheckpointInspection(
    const WVCheckpoint &checkpoint) {
  WVCheckpointInspection inspection;
  inspection.transformKind = checkpoint.transformKind;
  inspection.configuration = checkpoint.configuration;
  inspection.barotropicQGConfiguration =
      checkpoint.barotropicQGConfiguration;
  inspection.stateDescription = checkpoint.stateDescription;
  inspection.coefficientShape =
      checkpoint.transformKind == WVPersistedTransformKind::barotropicQG
          ? WVShape2D{1,
                      checkpoint.transformState.coefficientFamilies.empty()
                          ? 0
                          : checkpoint.transformState.coefficientFamilies[0]
                                .values.size()}
          : checkpoint.state.coefficients.shape;
  inspection.t = checkpoint.state.t;
  inspection.t0 = checkpoint.state.t0;
  inspection.metadata = checkpoint.metadata;
  inspection.forcingSchedule = checkpoint.forcingSchedule;
  return inspection;
}

const WVTransformConstantStratificationConfiguration *
legacyModelOutputPlanningConfiguration(
    const WVCheckpointInspection &inspection) noexcept {
  return inspection.transformKind ==
                 WVPersistedTransformKind::constantStratification
             ? &inspection.configuration
             : nullptr;
}

WVKernelStatus validateModelCheckpointState(
    WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout) {
  if (checkpoint.transformKind ==
      WVPersistedTransformKind::constantStratification) {
    if (!layout.hasLegacyCoefficientTriple())
      return invalid("A legacy checkpoint requires the constant-stratification "
                     "coefficient layout.");
    const auto expectedShape = layout.coefficientShape();
    if (checkpoint.state.coefficients.shape.rows != expectedShape.rows ||
        checkpoint.state.coefficients.shape.columns != expectedShape.columns)
      return invalid("Checkpoint coefficients do not match the model state "
                     "layout.");
    const auto count = expectedShape.elementCount();
    if (checkpoint.state.coefficients.Ap.size() != count ||
        checkpoint.state.coefficients.Am.size() != count ||
        checkpoint.state.coefficients.A0.size() != count)
      return invalid("Checkpoint coefficient storage is incomplete.");
    return WVKernelStatus::ok();
  }
  if (layout.hasLegacyCoefficientTriple() ||
      checkpoint.transformState.transformIdentifier !=
          layout.transformIdentifier() ||
      checkpoint.transformState.coefficientFamilies.size() !=
          layout.coefficientFamilyCount())
    return invalid("Transform checkpoint identity does not match the model "
                   "state layout.");
  for (std::size_t family = 0; family < layout.coefficientFamilyCount();
       ++family) {
    const auto &expected = layout.coefficientFamilies()[family];
    const auto &actual = checkpoint.transformState.coefficientFamilies[family];
    if (actual.identifier != expected.identifier ||
        actual.spectralDimensions != expected.spectralDimensions ||
        actual.values.size() != expected.elementCount)
      return invalid("Transform checkpoint coefficient families do not match "
                     "the model state layout.");
  }
  checkpoint.state.t = checkpoint.transformState.t;
  checkpoint.state.t0 = checkpoint.transformState.t0;
  return WVKernelStatus::ok();
}

WVComplex64 *modelCheckpointCoefficientData(
    WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout,
    std::size_t family) noexcept {
  if (family >= layout.coefficientFamilyCount())
    return nullptr;
  if (!layout.hasLegacyCoefficientTriple())
    return family < checkpoint.transformState.coefficientFamilies.size()
               ? checkpoint.transformState.coefficientFamilies[family]
                     .values.data()
               : nullptr;
  WVComplex64 *values[] = {checkpoint.state.coefficients.Ap.data(),
                           checkpoint.state.coefficients.Am.data(),
                           checkpoint.state.coefficients.A0.data()};
  return family < 3 ? values[family] : nullptr;
}

const WVComplex64 *modelCheckpointCoefficientData(
    const WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout,
    std::size_t family) noexcept {
  return modelCheckpointCoefficientData(
      const_cast<WVCheckpoint &>(checkpoint), layout, family);
}

WVMutableState modelCheckpointLegacyView(WVCheckpoint &checkpoint) noexcept {
  WVMutableState state;
  state.t = checkpoint.state.t;
  state.t0 = checkpoint.state.t0;
  if (checkpoint.transformKind !=
      WVPersistedTransformKind::constantStratification)
    return state;
  const auto shape = checkpoint.state.coefficients.shape;
  state.coefficients = {{checkpoint.state.coefficients.Ap.data(), shape},
                        {checkpoint.state.coefficients.Am.data(), shape},
                        {checkpoint.state.coefficients.A0.data(), shape}};
  return state;
}

void setModelCheckpointTimes(WVCheckpoint &checkpoint, double t,
                             double t0) noexcept {
  checkpoint.state.t = t;
  checkpoint.state.t0 = t0;
  if (checkpoint.transformKind == WVPersistedTransformKind::barotropicQG) {
    checkpoint.transformState.t = t;
    checkpoint.transformState.t0 = t0;
  }
}

} // namespace wavevortex::runtime::detail
