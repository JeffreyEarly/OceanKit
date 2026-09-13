#include "WVModelTransformAdapters.hpp"
#include "WaveVortexRuntime/WVStratifiedQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVHydrostaticIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVBoussinesqIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"

#include <new>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

template<class System>
WVKernelStatus setVariableEvaluationPolicyTransaction(
    System& system,WVVariableEvaluationPolicy policy) {
  auto* fields=system.fieldEvaluationService();
  auto status=system.validateVariableEvaluationPolicyChange(policy);
  if(!status) return status;
  if(fields) {
    status=fields->validateVariableEvaluationPolicyChange(policy);
    if(!status) return status;
  }
  if(policy==WVVariableEvaluationPolicy::lowMemory) {
    if(fields) {
      status=fields->setVariableEvaluationPolicy(policy);
      if(!status) return status;
    }
    // Every forcing family validates above and commits low memory without
    // allocation, so this second commit cannot strand mismatched policies.
    return system.setVariableEvaluationPolicy(policy);
  }
  status=system.setVariableEvaluationPolicy(policy);
  if(!status) return status;
  if(!fields) return status;
  status=fields->setVariableEvaluationPolicy(policy);
  if(!status) {
    // Reuse preparation is transactional. Its rollback releases caches and
    // therefore cannot allocate or fail after the successful preflight.
    (void)system.setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::lowMemory);
  }
  return status;
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
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) override {
    return setVariableEvaluationPolicyTransaction(*system_,policy);
  }
  void setLinearDynamics(bool linear) noexcept override {
    system_->setLinearDynamics(linear);
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
    metrics.variableEvaluation=system_->variableEvaluationMetrics();
    if (const auto* fields=system_->fieldEvaluationService()) {
      const auto outputProducers=fields->producerMetrics();
      metrics.variableProducers.horizontalSpeedReductions=outputProducers.horizontalSpeedReductions;
      metrics.variableProducers.verticalSpeedReductions=outputProducers.verticalSpeedReductions;
      metrics.variableProducers.energyReductions=outputProducers.energyReductions;
      metrics.variableProducers.reconstructions=outputProducers.reconstructions;
    }
    metrics.variableProducers.horizontalSpeedReductions += system_->forcingMetrics().horizontalSpeedReductionCount;
    metrics.variableProducers.verticalSpeedReductions += system_->forcingMetrics().verticalSpeedReductionCount;
    const auto& producers=system_->kernelMetrics();
    metrics.variableProducers.stateValidations=producers.stateValidationCount;
    metrics.variableProducers.phasePreparations=producers.phasePreparationCount;
    metrics.variableProducers.derivedValidations=producers.derivedValidationCount;
    metrics.variableProducers.tendencyReconstructions=producers.tendencyReconstructionCount;
    if(!system_->fieldEvaluationService())
      metrics.variableProducers.reconstructions=producers.reconstructionCount;
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
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) override {
    return setVariableEvaluationPolicyTransaction(*system_,policy);
  }
  void setLinearDynamics(bool linear) noexcept override {
    system_->setLinearDynamics(linear);
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
    metrics.variableEvaluation=system_->variableEvaluationMetrics();
    if (const auto* fields=system_->fieldEvaluationService()) {
      const auto outputProducers=fields->producerMetrics();
      metrics.variableProducers.horizontalSpeedReductions=outputProducers.horizontalSpeedReductions;
      metrics.variableProducers.verticalSpeedReductions=outputProducers.verticalSpeedReductions;
      metrics.variableProducers.energyReductions=outputProducers.energyReductions;
      metrics.variableProducers.reconstructions=outputProducers.reconstructions;
    }
    const auto& producers=system_->kernel().metrics();
    if (!system_->fieldEvaluationService())
      metrics.variableProducers.horizontalSpeedReductions=producers.horizontalSpeedMaximumReductionCount;
    metrics.variableProducers.stateValidations=producers.stateValidationCount;
    if(!system_->fieldEvaluationService()) {
      for(std::size_t field=0;field<producers.componentReconstructionCount.size();++field)
        for(std::size_t derivative=0;derivative<producers.componentReconstructionCount[field].size();++derivative)
          metrics.variableProducers.reconstructions[field][derivative]=
              producers.componentReconstructionCount[field][derivative];
    }
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
    metrics.forcing.horizontalSpeedReductionCount=forcing.horizontalSpeedReductionCount;
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

class StratifiedQGModelSystem final : public WVResolvedModelSystem {
public:
  explicit StratifiedQGModelSystem(
      std::unique_ptr<WVStratifiedQGIntegrationSystem> system)
      : system_(std::move(system)) {}

  WVIntegrationSystem &integrationSystem() noexcept override {
    return *system_;
  }
  const WVIntegrationSystem &integrationSystem() const noexcept override {
    return *system_;
  }
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) override {
    return setVariableEvaluationPolicyTransaction(*system_,policy);
  }
  void setLinearDynamics(bool linear) noexcept override {
    system_->setLinearDynamics(linear);
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
    metrics.variableEvaluation=system_->variableEvaluationMetrics();
    if (const auto* fields=system_->fieldEvaluationService()) {
      const auto outputProducers=fields->producerMetrics();
      metrics.variableProducers.horizontalSpeedReductions=outputProducers.horizontalSpeedReductions;
      metrics.variableProducers.verticalSpeedReductions=outputProducers.verticalSpeedReductions;
      metrics.variableProducers.energyReductions=outputProducers.energyReductions;
      metrics.variableProducers.reconstructions=outputProducers.reconstructions;
    }
    const auto& producers=system_->kernel().metrics();
    if (!system_->fieldEvaluationService())
      metrics.variableProducers.horizontalSpeedReductions=producers.horizontalSpeedMaximumReductionCount;
    metrics.variableProducers.horizontalSpeedReductions += system_->forcingMetrics().horizontalSpeedMaximumReductionCount;
    metrics.variableProducers.stateValidations=producers.stateValidationCount;
    if(!system_->fieldEvaluationService()) {
      for(std::size_t field=0;field<producers.componentReconstructionCount.size();++field)
        for(std::size_t derivative=0;derivative<producers.componentReconstructionCount[field].size();++derivative)
          metrics.variableProducers.reconstructions[field][derivative]=
              producers.componentReconstructionCount[field][derivative];
    }
    metrics.integratedObservers = system_->metrics();
    const auto& storage = system_->kernel().storage();
    metrics.kernel.descriptorBytes = storage.sharedScientificBytes + storage.factorBytes;
    metrics.kernel.planCount = 2;
    metrics.kernel.planBytes = storage.planBytesLowerBound;
    metrics.kernel.engineBytes = storage.providerBytesLowerBound;
    metrics.kernel.kernelManagementBytes = storage.preparedBytes;
    metrics.kernel.scratchCapacityBytes = storage.workspaceBytes + storage.spectralScratchBytes + storage.realScratchBytes;
    metrics.kernel.scratchHighWaterBytes = metrics.kernel.scratchCapacityBytes;
    metrics.kernel.realScratchCapacityBytes = storage.realScratchBytes;
    const auto &forcing = system_->forcingMetrics();
    metrics.forcing.horizontalSpeedReductionCount = forcing.horizontalSpeedMaximumReductionCount;
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
  std::unique_ptr<WVStratifiedQGIntegrationSystem> system_;
};
class HydrostaticModelSystem final : public WVResolvedModelSystem {
public:
  explicit HydrostaticModelSystem(
      std::unique_ptr<WVHydrostaticIntegrationSystem> system)
      : system_(std::move(system)) {}

  WVIntegrationSystem &integrationSystem() noexcept override {
    return *system_;
  }
  const WVIntegrationSystem &integrationSystem() const noexcept override {
    return *system_;
  }
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) override {
    return setVariableEvaluationPolicyTransaction(*system_,policy);
  }
  void setLinearDynamics(bool linear) noexcept override {
    system_->setLinearDynamics(linear);
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
    metrics.variableEvaluation=system_->variableEvaluationMetrics();
    if (const auto* fields=system_->fieldEvaluationService()) {
      const auto outputProducers=fields->producerMetrics();
      metrics.variableProducers.horizontalSpeedReductions=outputProducers.horizontalSpeedReductions;
      metrics.variableProducers.verticalSpeedReductions=outputProducers.verticalSpeedReductions;
      metrics.variableProducers.energyReductions=outputProducers.energyReductions;
      metrics.variableProducers.reconstructions=outputProducers.reconstructions;
    }
    metrics.variableProducers.horizontalSpeedReductions += system_->forcingMetrics().horizontalSpeedReductionCount;
    metrics.variableProducers.verticalSpeedReductions += system_->forcingMetrics().verticalSpeedReductionCount;
    const auto& producers=system_->kernel().metrics();
    metrics.variableProducers.stateValidations=producers.stateValidationCount;
    metrics.variableProducers.phasePreparations=producers.phasePreparationCount;
    metrics.variableProducers.derivedValidations=producers.derivedValidationCount;
    metrics.variableProducers.coefficientAssemblies=producers.coefficientAssemblyCount;
    metrics.variableProducers.verticalOperatorExecutions=producers.verticalOperatorExecutionCount;
    metrics.variableProducers.verticalPreparations=producers.verticalPreparationCount;
    metrics.variableProducers.horizontalSpectrumReuses=producers.horizontalSpectrumReuseCount;
    metrics.variableProducers.derivativeAdvectionConsumers=producers.derivativeAdvectionConsumerCount;
    metrics.variableProducers.tiledNonlinearExecutions=producers.tiledNonlinearCount;
    metrics.variableProducers.tiledColumnInverses=producers.tiledColumnInverseCount;
    metrics.variableProducers.tiledRowInverses=producers.tiledRowInverseCount;
    metrics.variableProducers.tiledReusedColumns=producers.tiledReusedColumnCount;

    metrics.variableProducers.preparedVerticalDerivatives=producers.preparedVerticalDerivativeCount;
    metrics.variableProducers.tendencyReconstructions=producers.tendencyReconstructionCount;
    if(!system_->fieldEvaluationService())
      metrics.variableProducers.reconstructions=producers.reconstructionCount;
    metrics.integratedObservers = system_->metrics();
    const auto& storage = system_->kernel().storage();
    metrics.kernel.descriptorBytes = storage.sharedScientificBytes + storage.factorBytes;
    metrics.kernel.planCount = 2;
    metrics.kernel.planBytes = storage.planBytesLowerBound;
    metrics.kernel.engineBytes = storage.providerBytesLowerBound;
    metrics.kernel.kernelManagementBytes = storage.preparedBytes;
    metrics.kernel.scratchCapacityBytes = storage.workspaceBytes + storage.spectralScratchBytes + storage.realScratchBytes;
    metrics.kernel.scratchHighWaterBytes = metrics.kernel.scratchCapacityBytes;
    metrics.kernel.realScratchCapacityBytes = storage.realScratchBytes;
    const auto &forcing = system_->forcingMetrics();
    metrics.forcing.nonlinearProducerCount = forcing.nonlinearProducerCount;
    metrics.forcing.horizontalSpeedReductionCount = forcing.horizontalSpeedReductionCount;
    metrics.forcing.verticalSpeedReductionCount = forcing.verticalSpeedReductionCount;
    metrics.forcing.gridCalculusProducerCount = forcing.gridCalculusProducerCount;
    metrics.forcing.constantLaplacianProducerCount = forcing.constantLaplacianProducerCount;
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
  std::unique_ptr<WVHydrostaticIntegrationSystem> system_;
};
class BoussinesqModelSystem final : public WVResolvedModelSystem {
public:
  explicit BoussinesqModelSystem(
      std::unique_ptr<WVBoussinesqIntegrationSystem> system)
      : system_(std::move(system)) {}

  WVIntegrationSystem &integrationSystem() noexcept override {
    return *system_;
  }
  const WVIntegrationSystem &integrationSystem() const noexcept override {
    return *system_;
  }
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) override {
    return setVariableEvaluationPolicyTransaction(*system_,policy);
  }
  void setLinearDynamics(bool linear) noexcept override {
    system_->setLinearDynamics(linear);
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
    metrics.variableEvaluation=system_->variableEvaluationMetrics();
    if (const auto* fields=system_->fieldEvaluationService()) {
      const auto outputProducers=fields->producerMetrics();
      metrics.variableProducers.horizontalSpeedReductions=outputProducers.horizontalSpeedReductions;
      metrics.variableProducers.verticalSpeedReductions=outputProducers.verticalSpeedReductions;
      metrics.variableProducers.energyReductions=outputProducers.energyReductions;
      metrics.variableProducers.reconstructions=outputProducers.reconstructions;
    }
    metrics.variableProducers.horizontalSpeedReductions += system_->forcingMetrics().horizontalSpeedReductionCount;
    metrics.variableProducers.verticalSpeedReductions += system_->forcingMetrics().verticalSpeedReductionCount;
    const auto& producers=system_->kernel().metrics();
    metrics.variableProducers.stateValidations=producers.stateValidationCount;
    metrics.variableProducers.phasePreparations=producers.phasePreparationCount;
    metrics.variableProducers.derivedValidations=producers.derivedValidationCount;
    metrics.variableProducers.coefficientAssemblies=producers.coefficientAssemblyCount;
    metrics.variableProducers.verticalOperatorExecutions=producers.verticalOperatorExecutionCount;
    metrics.variableProducers.verticalMatrixGroupExecutions=producers.verticalMatrixGroupExecutionCount;
    metrics.variableProducers.verticalPreparations=producers.verticalPreparationCount;
    metrics.variableProducers.horizontalSpectrumReuses=producers.horizontalSpectrumReuseCount;
    metrics.variableProducers.derivativeAdvectionConsumers=producers.derivativeAdvectionConsumerCount;
    metrics.variableProducers.tiledNonlinearExecutions=producers.tiledNonlinearCount;
    metrics.variableProducers.tiledColumnInverses=producers.tiledColumnInverseCount;
    metrics.variableProducers.tiledRowInverses=producers.tiledRowInverseCount;
    metrics.variableProducers.tiledReusedColumns=producers.tiledReusedColumnCount;

    metrics.variableProducers.preparedVerticalDerivatives=producers.preparedVerticalDerivativeCount;
    metrics.variableProducers.tendencyReconstructions=producers.tendencyReconstructionCount;
    if(!system_->fieldEvaluationService())
      metrics.variableProducers.reconstructions=producers.reconstructionCount;
    metrics.integratedObservers = system_->metrics();
    const auto& storage = system_->kernel().storage();
    metrics.kernel.descriptorBytes = storage.sharedScientificBytes + storage.factorBytes;
    metrics.kernel.planCount = 2;
    metrics.kernel.planBytes = storage.planBytesLowerBound;
    metrics.kernel.engineBytes = storage.providerBytesLowerBound;
    metrics.kernel.kernelManagementBytes = storage.preparedBytes;
    metrics.kernel.scratchCapacityBytes = storage.workspaceBytes + storage.spectralScratchBytes + storage.realScratchBytes;
    metrics.kernel.scratchHighWaterBytes = metrics.kernel.scratchCapacityBytes;
    metrics.kernel.realScratchCapacityBytes = storage.realScratchBytes;
    const auto &forcing = system_->forcingMetrics();
    metrics.forcing.nonlinearProducerCount = forcing.nonlinearProducerCount;
    metrics.forcing.horizontalSpeedReductionCount = forcing.horizontalSpeedReductionCount;
    metrics.forcing.verticalSpeedReductionCount = forcing.verticalSpeedReductionCount;
    metrics.forcing.gridCalculusProducerCount = forcing.gridCalculusProducerCount;
    metrics.forcing.constantLaplacianProducerCount = forcing.constantLaplacianProducerCount;
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
  std::unique_ptr<WVBoussinesqIntegrationSystem> system_;
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

WVKernelStatus createStratifiedQGModelSystem(
    std::shared_ptr<const WVStratifiedModalSource> source,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system,
    const WVVariableKernelServices &services) {
  std::unique_ptr<WVStratifiedQGIntegrationSystem> numerical;
  auto status = descriptor == nullptr
                    ? WVStratifiedQGIntegrationSystem::create(
                          std::move(source), schedule, std::move(catalog),
                          std::move(engine), numerical, services)
                    : WVStratifiedQGIntegrationSystem::create(
                          std::move(source), schedule, *descriptor,
                          std::move(catalog), std::move(engine), numerical, services);
  if (!status)
    return status;
  try {
    system =
        std::make_unique<StratifiedQGModelSystem>(std::move(numerical));
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Stratified QG model adapter."};
  }
}

WVKernelStatus createHydrostaticModelSystem(
    std::shared_ptr<const WVStratifiedModalSource> source,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system,
    const WVVariableKernelServices &services) {
  std::unique_ptr<WVHydrostaticIntegrationSystem> numerical;
  auto status = descriptor == nullptr
                    ? WVHydrostaticIntegrationSystem::create(
                          std::move(source), schedule, std::move(catalog),
                          std::move(engine), numerical, services)
                    : WVHydrostaticIntegrationSystem::create(
                          std::move(source), schedule, *descriptor,
                          std::move(catalog), std::move(engine), numerical, services);
  if (!status)
    return status;
  try {
    system =
        std::make_unique<HydrostaticModelSystem>(std::move(numerical));
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Stratified QG model adapter."};
  }
}
WVKernelStatus createBoussinesqModelSystem(
    std::shared_ptr<const WVStratifiedModalSource> source,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system,
    const WVVariableKernelServices &services) {
  std::unique_ptr<WVBoussinesqIntegrationSystem> numerical;
  auto status = descriptor == nullptr
                    ? WVBoussinesqIntegrationSystem::create(
                          std::move(source), schedule, std::move(catalog),
                          std::move(engine), numerical, services)
                    : WVBoussinesqIntegrationSystem::create(
                          std::move(source), schedule, *descriptor,
                          std::move(catalog), std::move(engine), numerical, services);
  if (!status)
    return status;
  try {
    system =
        std::make_unique<BoussinesqModelSystem>(std::move(numerical));
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Boussinesq model adapter."};
  }
}

WVKernelStatus createPersistedModelSystem(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system,
    const WVVariableKernelServices &services) {
  if (inspection.transformKind==WVPersistedTransformKind::hydrostatic) {
    if (!inspection.stratifiedModalSource || inspection.stratifiedModalSource->N2FunctionPayload().empty()) return invalid("MATLAB-compatible Hydrostatic output requires its N2Function payload.");
    return createHydrostaticModelSystem(inspection.stratifiedModalSource,schedule,descriptor,std::move(catalog),std::move(engine),system,services);
  }
  if (inspection.transformKind==WVPersistedTransformKind::boussinesq) {
    if (!inspection.stratifiedModalSource || inspection.stratifiedModalSource->N2FunctionPayload().empty()) return invalid("MATLAB-compatible Boussinesq output requires its N2Function payload.");
    return createBoussinesqModelSystem(inspection.stratifiedModalSource,schedule,descriptor,std::move(catalog),std::move(engine),system,services);
  }
  if (inspection.transformKind == WVPersistedTransformKind::stratifiedQG) {
    if (!inspection.stratifiedModalSource || inspection.stratifiedModalSource->N2FunctionPayload().empty())
      return invalid("MATLAB-compatible SQG model output requires its opaque N2Function persistence payload.");
    return createStratifiedQGModelSystem(inspection.stratifiedModalSource,schedule,descriptor,std::move(catalog),std::move(engine),system,services);
  }
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
  if (inspection.transformKind==WVPersistedTransformKind::hydrostatic) {
    if (!inspection.stratifiedModalSource) return invalid("Hydrostatic scientific source is absent.");
    return WVHydrostaticForcingEngine::validateSchedule(inspection.stratifiedModalSource->geometry(),schedule,inspection.coefficientShape,catalog);
  }
  if (inspection.transformKind==WVPersistedTransformKind::boussinesq) {
    if (!inspection.stratifiedModalSource) return invalid("Boussinesq scientific source is absent.");
    return WVBoussinesqForcingEngine::validateSchedule(inspection.stratifiedModalSource->geometry(),schedule,inspection.coefficientShape,catalog);
  }
  if (inspection.transformKind == WVPersistedTransformKind::stratifiedQG) {
    if (!inspection.stratifiedModalSource) return invalid("SQG scientific source is absent.");
    return WVStratifiedQGForcingEngine::validateSchedule(inspection.stratifiedModalSource->geometry(),schedule,inspection.coefficientShape.elementCount(),catalog);
  }
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
  inspection.stratifiedModalSource = checkpoint.stratifiedModalSource;
  inspection.stateDescription = checkpoint.stateDescription;
  inspection.coefficientShape =
      checkpoint.transformKind == WVPersistedTransformKind::barotropicQG
          ? WVShape2D{1,
                      checkpoint.transformState.coefficientFamilies.empty()
                          ? 0
                          : checkpoint.transformState.coefficientFamilies[0]
                                .values.size()}
          : checkpoint.state.coefficients.shape;
  if ((checkpoint.transformKind == WVPersistedTransformKind::stratifiedQG || (checkpoint.transformKind == WVPersistedTransformKind::hydrostatic || checkpoint.transformKind == WVPersistedTransformKind::boussinesq)) && checkpoint.stratifiedModalSource) {
    const auto& g = checkpoint.stratifiedModalSource->geometry();
    inspection.coefficientShape = {g.Nj,g.Nkl};
  }
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
  if (checkpoint.transformState.transformIdentifier !=
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
  if (checkpoint.transformKind!=WVPersistedTransformKind::constantStratification)
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
  if ((checkpoint.transformKind==WVPersistedTransformKind::hydrostatic || checkpoint.transformKind==WVPersistedTransformKind::boussinesq) && checkpoint.transformState.coefficientFamilies.size()==3) {
    const auto shape=WVShape2D{checkpoint.stratifiedModalSource->geometry().Nj,checkpoint.stratifiedModalSource->geometry().Nkl};
    auto& families=checkpoint.transformState.coefficientFamilies;
    state.coefficients={{families[0].values.data(),shape},{families[1].values.data(),shape},{families[2].values.data(),shape}};
    return state;
  }
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
  if (checkpoint.transformKind != WVPersistedTransformKind::constantStratification) {
    checkpoint.transformState.t = t;
    checkpoint.transformState.t0 = t0;
  }
}

} // namespace wavevortex::runtime::detail
