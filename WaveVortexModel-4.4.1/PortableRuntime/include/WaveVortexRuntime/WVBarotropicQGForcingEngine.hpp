#pragma once
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"

#include "WaveVortexRuntime/WVForcingSchedule.hpp"
#include "WaveVortexRuntime/WVForcingTendency.hpp"
#include "WaveVortexRuntime/WVIntegrationContracts.hpp"
#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVExtensionCatalog;
struct WVForcingEvaluationDependencies;
class WVBarotropicQGForcingEngine;

struct WVBarotropicQGFixedAmplitudeConfiguration {
  std::vector<std::size_t> A0Indices;
  std::vector<WVComplex64> A0Values;
};

struct WVBarotropicQGForcingEngineMetrics {
  std::size_t scheduleBytes = 0;
  std::size_t derivedOperatorBytes = 0;
  std::size_t workspaceCapacityBytes = 0;
  std::size_t evaluationCount = 0;
  std::size_t forcingCallCount = 0;
  std::size_t constraintOperationCount = 0;
  std::size_t restoredCoefficientCount = 0;
  std::size_t resolvedSpatialCount = 0;
  std::size_t resolvedSpectralCount = 0;
  std::size_t resolvedAmplitudeCount = 0;
  std::size_t physicalFieldReconstructionCount = 0;
  std::size_t physicalFieldReuseCount = 0;
  std::size_t horizontalSpeedReductionCount = 0;
  std::size_t spatialTendencyProjectionCount = 0;
  std::size_t stateConstraintElementWrites = 0;
};

class WVBarotropicQGForcingExecutionContext final {
public:
  WVKernelStatus nonlinearAdvection();
  void filterTendency(const std::vector<std::size_t> &indices);
  WVKernelStatus adaptiveDamping(const std::vector<double> &dampingOperator);
  WVKernelStatus linearBottomFriction(double rate);
  WVKernelStatus quadraticBottomFriction(double drag);
  WVKernelStatus betaPlanePVAdvection(double beta);
  void zeroSelectedTendencies(
      const WVBarotropicQGFixedAmplitudeConfiguration &configuration);

private:
  WVBarotropicQGForcingEngine *engine_ = nullptr;
  WVComplexConstView A0_;
  WVComplexView F0_;
  bool outputInitialized_ = false;
  WVBarotropicQGOperationWorkspace workspace_;
  detail::WVForcingDiagnosticWorkspace* diagnosticWorkspace_ = nullptr;
  friend class WVBarotropicQGForcingEngine;
};

// One exact identity/version implementation resolved before integration.
// Calls cross this interface once per forcing stage operation; all modal and
// grid loops execute below the coarse QG operation service.
class WVBarotropicQGForcing {
public:
  virtual ~WVBarotropicQGForcing() = default;
  virtual const std::string &typeIdentifier() const noexcept = 0;
  virtual std::uint32_t contractVersion() const noexcept = 0;
  virtual const std::string &name() const noexcept = 0;
  virtual WVForcingStage stage() const noexcept = 0;
  virtual std::uint8_t priority() const noexcept = 0;
  virtual std::size_t ordinal() const noexcept = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
  virtual bool supportsTendencyDiagnostics() const noexcept { return false; }
  virtual bool requiresDiagnosticPhysicalFields() const noexcept { return false; }
  virtual std::size_t constraintWriteCount() const noexcept { return 0; }
  virtual WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const = 0;
  virtual WVStateConstraintResult applyConstraint(WVComplexView &) const {
    return {WVKernelStatus::ok(), 0, true};
  }
};

class WVBarotropicQGForcingEngine final {
public:
  static WVKernelStatus validateSchedule(
      const WVTransformBarotropicQGConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule, std::size_t coefficientCount,
      const WVExtensionCatalog &catalog);

  static WVKernelStatus create(
      const WVTransformBarotropicQGConfiguration &configuration,
      const WVFrozenForcingSchedule &schedule,
      std::shared_ptr<const WVExtensionCatalog> catalog,
      std::unique_ptr<WVFFTEngine> fftEngine,
      std::unique_ptr<WVBarotropicQGForcingEngine> &forcingEngine);

  ~WVBarotropicQGForcingEngine();
  WVBarotropicQGForcingEngine(const WVBarotropicQGForcingEngine &) = delete;
  WVBarotropicQGForcingEngine &
  operator=(const WVBarotropicQGForcingEngine &) = delete;

  WVKernelStatus evaluateRightHandSide(const WVComplexConstView &A0,
                                       WVComplexView &F0,
                                       WVRealFieldBundleConstView *
                                           advectionFields = nullptr);
  // The borrowed A0 array must remain immutable until endStateEvaluation().
  WVKernelStatus beginStateEvaluation(const WVComplexConstView &A0);
  WVKernelStatus endStateEvaluation();
  bool stateEvaluationActive() const noexcept { return evaluation_.active(); }
  WVKernelStatus
  validateStateEvaluation(const WVComplexConstView &A0) const noexcept;
  WVKernelStatus horizontalSpeedMaximum(const WVComplexConstView &A0,
                                        double &maximum);
  WVStateConstraintResult restoreForcingAmplitudes(WVComplexView &A0);

  const WVTransformBarotropicQGKernel &kernel() const noexcept {
    return *kernel_;
  }
  WVTransformBarotropicQGKernel &kernel() noexcept { return *kernel_; }
  const WVBarotropicQGForcingEngineMetrics &metrics() const noexcept {
    return metrics_;
  }
  const std::string &scheduleIdentifier() const noexcept {
    return scheduleIdentifier_;
  }
  std::size_t persistentBytes() const noexcept;
  std::size_t forcingCount() const noexcept { return forcing_.size(); }
  const WVBarotropicQGForcing* forcingInstance(std::size_t index) const noexcept {
    return index<forcing_.size() ? forcing_[index].get() : nullptr;
  }
  const WVForcingEvaluationDependencies*
  forcingEvaluationDependencies(std::size_t index) const noexcept;
  // Optional u/v fields must describe this exact state and time.
  // They are borrowed for this invocation and must not alias state or outputs.
  WVKernelStatus evaluateForcingTendencies(const WVComplexConstView&,
      const WVForcingTendencyOutput*,std::size_t,
      const WVRealFieldBundleConstView* preparedPhysical = nullptr,
      detail::WVForcingDiagnosticWorkspace* session = nullptr);
  const WVForcingTendencyMetrics& tendencyMetrics() const noexcept { return tendencyMetrics_; }

  // Linear evolution retains instances for diagnostics and amplitude constraints.
  // Only their ordinary coefficient RHS contributions are disabled.
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy);
  WVKernelStatus validateVariableEvaluationPolicyChange(
      WVVariableEvaluationPolicy policy) const noexcept;
  const WVVariableEvaluationMetrics& variableEvaluationMetrics() const noexcept { return evaluation_.metrics(); }
  void setLinearDynamics(bool linear) noexcept { linearDynamics_ = linear; }

private:
  bool linearDynamics_ = false;
  WVVariableEvaluationContext evaluation_;
  WVVariableEvaluationPolicy evaluationPolicy_=WVVariableEvaluationPolicy::reuse;
  WVComplexConstView evaluationState_{};
  bool evaluationOwnsKernelScope_=false;
  double horizontalSpeedMaximum_ = 0.0;
  std::vector<WVComplex64> nonlinearScratch_;
  WVBarotropicQGForcingEngine() = default;
  WVKernelStatus initialize(const WVFrozenForcingSchedule &schedule);
  void initializeOutputWithZeros(WVComplexView &F0);

  std::unique_ptr<WVTransformBarotropicQGKernel> kernel_;
  std::shared_ptr<const WVExtensionCatalog> catalog_;
  std::vector<std::unique_ptr<WVBarotropicQGForcing>> forcing_;
  WVBarotropicQGForcingEngineMetrics metrics_;
  std::string scheduleIdentifier_;
  WVForcingTendencyMetrics tendencyMetrics_;
  bool executing_ = false;
  friend class WVBarotropicQGForcingExecutionContext;
};

} // namespace wavevortex::runtime
