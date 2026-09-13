#pragma once

#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVLagrangianParticles.hpp"
#include "WaveVortexRuntime/WVTracer.hpp"
#include "WaveVortexRuntime/WVVariableKernelServices.hpp"

#include <memory>
#include <vector>

namespace wavevortex::runtime {

// Hydrostatic integration numerical system. The canonical
// WaveVortex coefficient RHS is delegated to the frozen forcing engine while
// integrated observers consume one shared per-RHS evaluation context.
class WVHydrostaticIntegrationSystem final
    : public WVIntegrationSystem {
public:
  static WVKernelStatus create(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVHydrostaticIntegrationSystem> &system,
      const WVVariableKernelServices &services = {});
  static WVKernelStatus create(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor &descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVHydrostaticIntegrationSystem> &system,
      const WVVariableKernelServices &services = {});

  ~WVHydrostaticIntegrationSystem() override;
  WVHydrostaticIntegrationSystem(
      const WVHydrostaticIntegrationSystem &) = delete;
  WVHydrostaticIntegrationSystem &operator=(
      const WVHydrostaticIntegrationSystem &) = delete;

  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(const WVIntegrationState &state,
                                       WVIntegrationFlux &rightHandSide) override;
  WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &state) override;
  WVKernelStatus
  initializeParticleState(WVMutableIntegrationState &state) const;
  WVKernelStatus createErrorPolicy(
      double absoluteToleranceScale,
      std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override;
  bool supportsFixedTimeStepSelection() const noexcept override { return true; }
  WVKernelStatus evaluateFixedTimeStepCandidates(
      const WVIntegrationState &state, double cfl,
      WVFixedTimeStepCandidates &candidates) override;

  const std::vector<WVLagrangianParticles> &particles() const noexcept {
    return particles_;
  }
  const std::vector<WVTracer> &tracers() const noexcept { return tracers_; }
  WVFieldEvaluationService *fieldEvaluationService() noexcept override {
    return fields_.get();
  }
  const WVIntegratedObserverMetrics &metrics() const noexcept {
    return metrics_;
  }
  const WVTransformHydrostaticKernel &kernel() const noexcept {
    return forcing_->kernel();
  }
  WVTransformHydrostaticKernel &kernel() noexcept {
    return forcing_->kernel();
  }
  const WVForcingEngineMetrics &forcingMetrics() const noexcept {
    return forcing_->metrics();
  }
  const std::string &scheduleIdentifier() const noexcept {
    return forcing_->scheduleIdentifier();
  }
  std::size_t persistentBytes() const noexcept override;

  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) {
    return forcing_->setVariableEvaluationPolicy(policy);
  }
  WVKernelStatus validateVariableEvaluationPolicyChange(
      WVVariableEvaluationPolicy policy) const noexcept {
    return forcing_->validateVariableEvaluationPolicyChange(policy);
  }
  const WVVariableEvaluationMetrics& variableEvaluationMetrics() const noexcept {
    return forcing_->variableEvaluationMetrics();
  }
  void setLinearDynamics(bool linear) noexcept { forcing_->setLinearDynamics(linear); }

private:
  static WVKernelStatus createImpl(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor *descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVHydrostaticIntegrationSystem> &system,
      const WVVariableKernelServices &services);
  WVHydrostaticIntegrationSystem() = default;
  WVIntegrationStateLayout layout_;
  std::unique_ptr<WVHydrostaticForcingEngine> forcing_;
  std::unique_ptr<WVFieldEvaluationService> fields_;
  WVMovingFieldEvaluationPlan velocityPlan_;
  std::vector<WVLagrangianParticles> particles_;
  std::vector<WVTracer> tracers_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> z_;
  std::vector<std::vector<double>> velocityStorage_;
  std::vector<WVFieldOutputView> velocityViews_;
  WVIntegratedObserverMetrics metrics_;
  bool executing_ = false;
};

} // namespace wavevortex::runtime
