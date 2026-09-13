#pragma once

#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationContracts.hpp"
#include "WaveVortexRuntime/WVVariableKernelServices.hpp"

#include <memory>

namespace wavevortex::runtime {

// Transform-specific numerical system over the #279 coefficient-family
// boundary. It owns only the compact A0 family declared by stateLayout().
class WVStratifiedQGIntegrationSystem final : public WVIntegrationSystem {
public:
  static WVKernelStatus create(
      std::shared_ptr<const WVStratifiedModalSource> source,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVStratifiedQGIntegrationSystem> &system,
      const WVVariableKernelServices &services = {});
  static WVKernelStatus create(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor &descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVStratifiedQGIntegrationSystem> &system,
      const WVVariableKernelServices &services = {});
  static WVKernelStatus create(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVStratifiedQGIntegrationSystem> &system,
      const WVVariableKernelServices &services = {});

  ~WVStratifiedQGIntegrationSystem() override = default;
  WVStratifiedQGIntegrationSystem(const WVStratifiedQGIntegrationSystem &) =
      delete;
  WVStratifiedQGIntegrationSystem &operator=(
      const WVStratifiedQGIntegrationSystem &) = delete;

  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return layout_;
  }
  WVKernelStatus evaluateRightHandSide(
      const WVIntegrationState &state,
      WVIntegrationFlux &rightHandSide) override;
  WVStateConstraintResult enforceStateConstraints(
      WVMutableIntegrationState &state) override;
  WVKernelStatus createErrorPolicy(
      double absoluteToleranceScale,
      std::unique_ptr<WVIntegrationErrorPolicy> &policy) const override;
  bool supportsFixedTimeStepSelection() const noexcept override { return true; }
  WVKernelStatus evaluateFixedTimeStepCandidates(
      const WVIntegrationState &state, double cfl,
      WVFixedTimeStepCandidates &candidates) override;
  std::size_t persistentBytes() const noexcept override;
  WVKernelStatus initializeParticleState(
      WVMutableIntegrationState &state) const;
  WVFieldEvaluationService *fieldEvaluationService() noexcept override {
    return fields_.get();
  }
  const WVIntegratedObserverMetrics &metrics() const noexcept {
    return observerMetrics_;
  }

  const WVTransformStratifiedQGKernel &kernel() const noexcept {
    return forcingEngine_->kernel();
  }
  WVTransformStratifiedQGKernel &kernel() noexcept {
    return forcingEngine_->kernel();
  }
  const WVStratifiedQGForcingEngineMetrics &forcingMetrics() const noexcept {
    return forcingEngine_->metrics();
  }
  const std::string &forcingScheduleIdentifier() const noexcept {
    return forcingEngine_->scheduleIdentifier();
  }

  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) {
    return forcingEngine_->setVariableEvaluationPolicy(policy);
  }
  WVKernelStatus validateVariableEvaluationPolicyChange(
      WVVariableEvaluationPolicy policy) const noexcept {
    return forcingEngine_->validateVariableEvaluationPolicyChange(policy);
  }
  const WVVariableEvaluationMetrics& variableEvaluationMetrics() const noexcept {
    return forcingEngine_->variableEvaluationMetrics();
  }
  void setLinearDynamics(bool linear) noexcept { forcingEngine_->setLinearDynamics(linear); }

private:
  static WVKernelStatus createImpl(
      std::shared_ptr<const WVStratifiedModalSource> source,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor *descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVStratifiedQGIntegrationSystem> &system,
      const WVVariableKernelServices &services);
  struct Particle {
    WVObserverRecord record;
    std::size_t xBlock = 0;
    std::size_t yBlock = 0;
    std::size_t positionOffset = 0;
    std::size_t particleCount = 0;
    std::size_t uOutput = 0;
    std::size_t vOutput = 0;
  };
  struct Tracer {
    WVObserverRecord record;
    std::size_t stateBlock = 0;
  };
  WVStratifiedQGIntegrationSystem() = default;
  WVIntegrationStateLayout layout_;
  std::unique_ptr<WVStratifiedQGForcingEngine> forcingEngine_;
  std::unique_ptr<WVFieldEvaluationService> fields_;
  WVMovingFieldEvaluationPlan velocityPlan_;
  std::vector<Particle> particles_;
  std::vector<Tracer> tracers_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> z_;
  std::vector<double> advectionStorage_;
  std::vector<std::vector<double>> velocityStorage_;
  std::vector<WVFieldOutputView> velocityViews_;
  WVIntegratedObserverMetrics observerMetrics_;
  bool executing_ = false;
};

} // namespace wavevortex::runtime
