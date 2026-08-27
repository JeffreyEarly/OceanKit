#pragma once

#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVLagrangianParticles.hpp"
#include "WaveVortexRuntime/WVTracer.hpp"

#include <memory>
#include <vector>

namespace wavevortex::runtime {

struct WVIntegratedObserverMetrics {
  std::size_t rightHandSideEvaluationCount = 0;
  std::size_t velocityFieldEvaluationCount = 0;
  std::size_t particleValueWriteCount = 0;
  std::size_t tracerEvaluationCount = 0;
  std::size_t tracerValueWriteCount = 0;
  std::size_t antialiasedTracerEvaluationCount = 0;
  std::size_t positionCapacityBytes = 0;
  std::size_t velocityCapacityBytes = 0;
  std::size_t persistentBytes = 0;
  std::size_t sharedRightHandSideContextCount = 0;
  double rightHandSideSeconds = 0.0;
  double waveVortexFluxSeconds = 0.0;
  double additionalStateClearSeconds = 0.0;
  double tracerAdvectionSeconds = 0.0;
  double particleAdvectionSeconds = 0.0;
};

using WVLagrangianParticleMetrics = WVIntegratedObserverMetrics;

// Constant-stratification integration numerical system. The canonical
// WaveVortex coefficient RHS is delegated to the frozen forcing engine while
// integrated observers consume one shared per-RHS evaluation context.
class WVConstantStratificationIntegrationSystem final
    : public WVIntegrationSystem {
public:
  static WVKernelStatus create(
      const WVTransformConstantStratificationConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVConstantStratificationIntegrationSystem> &system);
  static WVKernelStatus create(
      const WVTransformConstantStratificationConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor &descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVConstantStratificationIntegrationSystem> &system);

  ~WVConstantStratificationIntegrationSystem() override;
  WVConstantStratificationIntegrationSystem(
      const WVConstantStratificationIntegrationSystem &) = delete;
  WVConstantStratificationIntegrationSystem &operator=(
      const WVConstantStratificationIntegrationSystem &) = delete;

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
  const WVKernelMetrics &kernelMetrics() const noexcept {
    return forcing_->kernel().metrics();
  }
  const WVTransformConstantStratificationKernel &kernel() const noexcept {
    return forcing_->kernel();
  }
  WVTransformConstantStratificationKernel &kernel() noexcept {
    return forcing_->kernel();
  }
  const WVForcingEngineMetrics &forcingMetrics() const noexcept {
    return forcing_->metrics();
  }
  const std::string &scheduleIdentifier() const noexcept {
    return forcing_->scheduleIdentifier();
  }
  std::size_t persistentBytes() const noexcept override;

private:
  static WVKernelStatus createImpl(
      const WVTransformConstantStratificationConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor *descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVConstantStratificationIntegrationSystem> &system);
  WVConstantStratificationIntegrationSystem() = default;
  WVIntegrationStateLayout layout_;
  std::unique_ptr<WVConstantStratificationForcingEngine> forcing_;
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
