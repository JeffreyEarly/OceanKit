#pragma once

#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"

#include <memory>

namespace wavevortex::runtime {

struct WVObserverOutputEvaluationMetrics {
  std::size_t preparedEventCount = 0;
  std::size_t fieldEvaluationCount = 0;
  std::size_t uniqueFieldOutputCount = 0;
  std::size_t sharedFieldReuseCount = 0;
  std::size_t borrowedCoefficientViewCount = 0;
  std::size_t outputCapacityBytes = 0;
  std::size_t retainedStorageBytes = 0;
  std::size_t routeAwareParticleEvaluationCount = 0;
  std::size_t skippedParticleEvaluationCount = 0;
  std::size_t batchRetainedStorageBytes = 0;
  std::size_t batchMaximumLiveBytes = 0;
  std::size_t occurrencePreparationCount = 0;
  std::size_t occurrenceReuseCount = 0;
  std::size_t occurrenceBatchBuildCount = 0;
  std::size_t occurrenceWorkspaceRetainedBytes = 0;
  std::size_t occurrenceWorkspaceLiveBytes = 0;
  std::size_t occurrenceWorkspaceMaximumLiveBytes = 0;
  double evaluationSeconds = 0.0;
};

// Evaluates passive MATLAB-compatible observing systems independently of
// NetCDF. Coefficients remain borrowed views of the immutable event state;
// coincident Eulerian and mooring requests share one field-evaluation plan.
class WVObserverOutputEvaluationService final : public WVObserverSampleSource {
public:
  ~WVObserverOutputEvaluationService() override;
  static WVKernelStatus
  create(const WVTransformConstantStratificationConfiguration &configuration,
         bool isDynamicsLinear,
         const WVPortableObserverDescriptor &descriptor,
         std::unique_ptr<WVFFTEngine> engine,
         std::unique_ptr<WVObserverOutputEvaluationService> &service,
         WVFieldEvaluationService *borrowedFieldEvaluationService = nullptr);
  static WVKernelStatus
  create(bool isDynamicsLinear,
         const WVPortableObserverDescriptor &descriptor,
         WVFieldEvaluationService &fieldEvaluationService,
         std::unique_ptr<WVObserverOutputEvaluationService> &service);

  WVKernelStatus observationSchema(
      const WVObserverRecord &observer,
      WVObservationSchema &output) override;
  WVKernelStatus initialObservationBatch(
      const WVObserverRecord &observer,
      WVObservationBatch &output) override;
  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const override;
  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &identity,
      const WVObserverRecord &observer,
      WVObservationBatch &output) override;
  WVKernelStatus preflight(const WVOutputPlan &plan) override;
  WVKernelStatus useFieldEvaluationService(
      WVFieldEvaluationService &fieldEvaluationService);
  WVKernelStatus prepareInitial(const WVState &state) override;
  WVKernelStatus prepareInitial(const WVIntegrationState &state) override;
  WVKernelStatus prepare(const WVOutputEvent &event) override;
  void complete(const WVOutputEvent &event) noexcept override;
  std::size_t occurrenceWorkspaceRetainedBytes() const noexcept override;
  std::size_t occurrenceWorkspaceLiveBytes() const noexcept override;

  const WVObserverOutputEvaluationMetrics &metrics() const noexcept {
    return metrics_;
  }
  std::size_t persistentBytes() const noexcept;

private:
  WVObserverOutputEvaluationService() = default;
  WVKernelStatus observationBatchForKind(
      const WVObservationOccurrenceIdentity *identity,
      const WVObserverRecord &observer, WVObservationBatchKind kind,
      WVObservationBatch &output);
  class Impl;
  std::unique_ptr<Impl> impl_;
  WVObserverOutputEvaluationMetrics metrics_;
};

} // namespace wavevortex::runtime
