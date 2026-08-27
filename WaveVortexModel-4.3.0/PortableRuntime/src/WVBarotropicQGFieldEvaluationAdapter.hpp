#pragma once

#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"

namespace wavevortex::runtime::detail {

// Named transform adapter behind the resolved WVFieldEvaluationService
// boundary. Plans retain only immutable interpolation metadata; the adapter
// owns one horizontal real scratch field and never a full Hermitian spectrum.
class WVBarotropicQGFieldEvaluationAdapter final {
public:
  ~WVBarotropicQGFieldEvaluationAdapter();

  static WVKernelStatus create(
      const WVTransformBarotropicQGConfiguration &configuration,
      std::unique_ptr<WVFFTEngine> engine,
      std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter> &adapter);
  static WVKernelStatus createBorrowing(
      WVTransformBarotropicQGKernel &kernel,
      std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter> &adapter);

  WVKernelStatus createPlan(const std::vector<WVFieldRequest> &requests,
                            WVFieldEvaluationPlan &plan) const;
  WVKernelStatus evaluate(const WVFieldEvaluationPlan &plan,
                          const WVIntegrationState &state,
                          WVFieldOutputView *outputs,
                          std::size_t outputCount);
  WVKernelStatus createMovingPlan(
      const std::vector<WVMovingFieldRequest> &requests,
      WVMovingFieldEvaluationPlan &plan) const;
  WVKernelStatus evaluateMoving(const WVMovingFieldEvaluationPlan &plan,
                                const WVIntegrationState &state,
                                WVMovingPositionView positions,
                                WVFieldOutputView *outputs,
                                std::size_t outputCount);
  WVKernelStatus evaluateMovingFromAdvectionFields(
      const WVMovingFieldEvaluationPlan &plan,
      const WVIntegrationState &state,
      const WVRealFieldBundleConstView &advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount);
  WVKernelStatus createEventPlan(
      const std::vector<WVEventFieldRequest> &requests,
      WVEventFieldEvaluationPlan &plan);
  WVKernelStatus prepareEventGeometry(
      const WVEventFieldEvaluationPlan &plan,
      const WVEventPositionSetView *positionSets,
      std::size_t positionSetCount,
      WVPreparedFieldGeometry &geometry);
  WVKernelStatus evaluateEventBatch(
      const WVIntegrationState &state,
      const WVEventFieldEvaluationBatchEntry *entries,
      std::size_t entryCount);

  bool isCompatibleWith(const WVIntegrationStateLayout &layout) const noexcept;
  const WVTransformBarotropicQGConfiguration &configuration() const noexcept;
  const WVFieldEvaluationMetrics &metrics() const noexcept { return metrics_; }
  std::size_t persistentBytes() const noexcept;

private:
  struct MovingInterpolationWorkspace;
  WVKernelStatus evaluateMovingImpl(
      const WVMovingFieldEvaluationPlan &plan,
      const WVIntegrationState &state,
      const WVRealFieldBundleConstView *advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount);
  WVBarotropicQGFieldEvaluationAdapter() = default;
  std::unique_ptr<WVTransformBarotropicQGKernel> ownedKernel_;
  WVTransformBarotropicQGKernel *kernel_ = nullptr;
  std::vector<double> fieldScratch_;
  std::unique_ptr<MovingInterpolationWorkspace> movingInterpolation_;
  WVFieldEvaluationMetrics metrics_;
  bool executing_ = false;
};

} // namespace wavevortex::runtime::detail
