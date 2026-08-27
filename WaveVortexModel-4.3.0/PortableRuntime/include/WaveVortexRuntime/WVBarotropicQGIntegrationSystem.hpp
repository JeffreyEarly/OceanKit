#pragma once

#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationContracts.hpp"

#include <memory>

namespace wavevortex::runtime {

// Allocation-light numerical record decoded by a transform-specific
// persistence adapter before compact coefficient storage is loaded. A later
// composition layer may connect it to NetCDF without changing this API.
struct WVBarotropicQGPersistedNumericalRecord {
  std::vector<double> x;
  std::vector<double> y;
  double Lx = 0.0;
  double Ly = 0.0;
  double h = 0.0;
  std::uint32_t j = 1;
  double g = 0.0;
  double planetaryRadius = 0.0;
  double rotationRate = 0.0;
  double latitude = 0.0;
  bool shouldAntialias = true;
  double t = 0.0;
  double t0 = 0.0;
};

struct WVBarotropicQGNumericalConfiguration {
  WVTransformBarotropicQGConfiguration transform;
  WVTransformStateDescription stateDescription;
  double t = 0.0;
  double t0 = 0.0;
  std::size_t persistentBytes() const noexcept;
};

WVKernelStatus decodeBarotropicQGNumericalConfiguration(
    const WVBarotropicQGPersistedNumericalRecord &record,
    WVBarotropicQGNumericalConfiguration &configuration);

// Transform-specific numerical system over the #279 coefficient-family
// boundary. It owns only the compact A0 family declared by stateLayout().
class WVBarotropicQGIntegrationSystem final : public WVIntegrationSystem {
public:
  static WVKernelStatus create(
      const WVTransformBarotropicQGConfiguration &configuration,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVBarotropicQGIntegrationSystem> &system);
  static WVKernelStatus create(
      const WVTransformBarotropicQGConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor &descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVBarotropicQGIntegrationSystem> &system);
  static WVKernelStatus create(
      const WVTransformBarotropicQGConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVBarotropicQGIntegrationSystem> &system);

  ~WVBarotropicQGIntegrationSystem() override = default;
  WVBarotropicQGIntegrationSystem(const WVBarotropicQGIntegrationSystem &) =
      delete;
  WVBarotropicQGIntegrationSystem &operator=(
      const WVBarotropicQGIntegrationSystem &) = delete;

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

  const WVTransformBarotropicQGKernel &kernel() const noexcept {
    return forcingEngine_->kernel();
  }
  WVTransformBarotropicQGKernel &kernel() noexcept {
    return forcingEngine_->kernel();
  }
  const WVBarotropicQGForcingEngineMetrics &forcingMetrics() const noexcept {
    return forcingEngine_->metrics();
  }
  const std::string &forcingScheduleIdentifier() const noexcept {
    return forcingEngine_->scheduleIdentifier();
  }

private:
  static WVKernelStatus createImpl(
      const WVTransformBarotropicQGConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      const WVPortableObserverDescriptor *descriptor,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVBarotropicQGIntegrationSystem> &system);
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
  WVBarotropicQGIntegrationSystem() = default;
  WVIntegrationStateLayout layout_;
  std::unique_ptr<WVBarotropicQGForcingEngine> forcingEngine_;
  std::unique_ptr<WVFieldEvaluationService> fields_;
  WVMovingFieldEvaluationPlan velocityPlan_;
  std::vector<Particle> particles_;
  std::vector<Tracer> tracers_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<std::vector<double>> velocityStorage_;
  std::vector<WVFieldOutputView> velocityViews_;
  WVIntegratedObserverMetrics observerMetrics_;
  bool executing_ = false;
};

} // namespace wavevortex::runtime
