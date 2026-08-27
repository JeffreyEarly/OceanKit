#pragma once

#include "WaveVortexRuntime/WVIntegrationContracts.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace wavevortex::runtime {

struct WVIntegratorMetrics {
  std::size_t workspaceCapacityBytes = 0;
  std::size_t workspaceMaximumLiveBytes = 0;
  std::size_t workspaceMaximumLiveStateEquivalentCount = 0;
  std::size_t rightHandSideEvaluationCount = 0;
  std::size_t baseRightHandSideEvaluationCount = 0;
  std::size_t continuousExtensionRightHandSideEvaluationCount = 0;
  std::size_t acceptedStepCount = 0;
  std::size_t rejectedStepCount = 0;
  std::size_t denseOutputEvaluationCount = 0;
  std::size_t denseOutputCacheBuildCount = 0;
  std::size_t denseOutputCacheReuseCount = 0;
  std::size_t errorPolicyBytes = 0;
  std::size_t diagnosticCapacityBytes = 0;
  std::size_t denseHistoryCapacityBytes = 0;
  std::size_t stateCapacityBytes = 0;
  std::size_t workspaceStateEquivalentCount = 0;
  std::size_t denseHistoryStateEquivalentCount = 0;
  std::size_t retainedBaseStageCapacityBytes = 0;
  std::size_t retainedBaseStageStateEquivalentCount = 0;
  std::size_t continuousExtensionWorkspaceCapacityBytes = 0;
  std::size_t continuousExtensionWorkspaceMaximumLiveBytes = 0;
  std::size_t continuousExtensionWorkspaceStateEquivalentCount = 0;
  std::size_t continuousExtensionWorkspaceMaximumLiveStateEquivalentCount = 0;
  double lastStepSize = 0.0;
  double lastProposedStepSize = 0.0;
  double lastAcceptedStepSize = 0.0;
  double normalizedError = 0.0;
  double nextStepSize = 0.0;
  std::size_t stepCount = 0;
  std::size_t fsalReuseCount = 0;
  std::size_t fsalInvalidationCount = 0;
  std::size_t rejectedInitialDerivativeReuseCount = 0;
  std::size_t constraintModifiedCoefficientCount = 0;
  std::size_t stageStateConstructionElementReads = 0;
  std::size_t stageStateConstructionElementWrites = 0;
  std::size_t stageFluxClearElementWrites = 0;
  std::size_t weightedFluxClearElementWrites = 0;
  std::size_t weightedFluxInitializationElementReads = 0;
  std::size_t weightedFluxInitializationElementWrites = 0;
  std::size_t weightedAccumulationElementReads = 0;
  std::size_t weightedAccumulationElementWrites = 0;
  std::size_t finalStateUpdateElementReads = 0;
  std::size_t finalStateUpdateElementWrites = 0;
  std::size_t acceptedStateCommitElementReads = 0;
  std::size_t acceptedStateCommitElementWrites = 0;
  std::size_t workspaceLiveBytes = 0;
  std::size_t denseOutputElementReads = 0;
  std::size_t denseOutputElementWrites = 0;
  double denseOutputSeconds = 0.0;
  double continuousExtensionSeconds = 0.0;
};

using WVFixedStepRK4Metrics = WVIntegratorMetrics;
using WVAdaptiveRK23Metrics = WVIntegratorMetrics;
using WVAdaptiveRK45Metrics = WVIntegratorMetrics;
using WVAdaptiveRK78Metrics = WVIntegratorMetrics;

struct WVFixedStepRK4Options {
  bool retainDenseOutput = false;
};

struct WVAdaptiveRKStageBufferLastUse;

class WVFixedStepRK4 final : public WVTimeIntegrator,
                                      public WVDenseOutput {
public:
  explicit WVFixedStepRK4(WVIntegrationSystem &system,
                          WVFixedStepRK4Options options = {});
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return system_.stateLayout();
  }
  WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) override;
  WVKernelStatus step(WVMutableIntegrationState &state,
                      double stepSize) override;
  WVKernelStatus advanceToTime(WVMutableIntegrationState &state, double finalTime,
                               double stepSize) override;
  WVKernelStatus evaluateDenseOutput(double time,
                                     WVMutableIntegrationState &output) const;
  WVKernelStatus
  evaluateState(double time, WVMutableIntegrationState &output) const override {
    return evaluateDenseOutput(time, output);
  }
  double initialTime() const noexcept override {
    return hasAcceptedStep_ ? acceptedStep_.initialTime : 0.0;
  }
  double finalTime() const noexcept override {
    return hasAcceptedStep_ ? acceptedStep_.finalTime : 0.0;
  }
  const WVAcceptedStep *lastAcceptedStep() const noexcept override {
    return hasAcceptedStep_ ? &acceptedStep_ : nullptr;
  }
  double nextStepSize() const noexcept override { return nextStepSize_; }
  std::size_t persistentBytes() const noexcept override;
  const WVIntegratorMetrics &metrics() const noexcept {
    return metrics_;
  }
  static const WVAdaptiveRKStageBufferLastUse *stageBufferLastUseRecords()
      noexcept;
  static std::size_t stageBufferLastUseRecordCount() noexcept;

private:
  class Workspace;
  WVKernelStatus ensureWorkspace(const WVMutableIntegrationState &state);
  WVIntegrationSystem &system_;
  WVFixedStepRK4Options options_;
  Workspace *workspace_ = nullptr;
  WVAcceptedStep acceptedStep_;
  mutable WVIntegratorMetrics metrics_;
  double nextStepSize_ = 0.0;
  bool hasAcceptedStep_ = false;
  bool acceptedStateConstrained_ = false;
  bool stepping_ = false;
  mutable bool evaluatingDenseOutput_ = false;

public:
  ~WVFixedStepRK4();
  WVFixedStepRK4(const WVFixedStepRK4 &) = delete;
  WVFixedStepRK4 &operator=(const WVFixedStepRK4 &) = delete;
};

struct WVAdaptiveRK23Options {
  double relativeTolerance = 1e-3;
  double absoluteToleranceScale = 1e-6;
  double safetyFactor = 0.8;
  double rejectionFloorFactor = 0.5;
  double maximumStepFactor = 5.0;
  double maximumStepSize = std::numeric_limits<double>::infinity();
  std::size_t maximumRecordedStepDiagnostics = 512;
  bool retainDenseOutput = true;
};

struct WVAdaptiveRKStepDiagnostic {
  double initialTime = 0.0;
  double acceptedStepSize = 0.0;
  double normalizedError = 0.0;
  double nextStepSize = 0.0;
  std::size_t rejectedAttemptCount = 0;
  std::size_t rightHandSideEvaluationCount = 0;
  bool reusedFSALDerivative = false;
};

using WVAdaptiveRK23StepDiagnostic = WVAdaptiveRKStepDiagnostic;

struct WVAdaptiveRKStageBufferLastUse {
  const char *bufferIdentifier = nullptr;
  const char *producer = nullptr;
  const char *lastUse = nullptr;
  std::size_t stage = 0;
};

class WVAdaptiveRK23 final : public WVTimeIntegrator,
                             public WVDenseOutput {
public:
  explicit WVAdaptiveRK23(WVIntegrationSystem &system,
                          WVAdaptiveRK23Options options = {});
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return system_.stateLayout();
  }
  WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) override;
  WVKernelStatus step(WVMutableIntegrationState &state,
                      double proposedStepSize) override;
  WVKernelStatus advanceToTime(WVMutableIntegrationState &state, double finalTime,
                               double initialStepSize) override;
  WVKernelStatus evaluateDenseOutput(double time,
                                     WVMutableIntegrationState &output) const;
  WVKernelStatus
  evaluateState(double time, WVMutableIntegrationState &output) const override {
    return evaluateDenseOutput(time, output);
  }
  double initialTime() const noexcept override {
    return hasAcceptedStep_ ? acceptedStep_.initialTime : 0.0;
  }
  double finalTime() const noexcept override {
    return hasAcceptedStep_ ? acceptedStep_.finalTime : 0.0;
  }
  const WVAcceptedStep *lastAcceptedStep() const noexcept override {
    return hasAcceptedStep_ ? &acceptedStep_ : nullptr;
  }
  const WVIntegratorMetrics &metrics() const noexcept { return metrics_; }
  const std::vector<WVAdaptiveRK23StepDiagnostic> &stepDiagnostics() const
      noexcept {
    return stepDiagnostics_;
  }
  std::uint64_t toleranceHash() const noexcept { return toleranceHash_; }
  const std::vector<std::uint64_t> &toleranceComponentHashes() const noexcept {
    return toleranceComponentHashes_;
  }
  bool stepDiagnosticsComplete() const noexcept {
    return stepDiagnostics_.size() == metrics_.acceptedStepCount;
  }
  static constexpr const char *controllerIdentifier() noexcept {
    return "matlab-ode23-v1";
  }
  static constexpr const char *methodIdentifier() noexcept {
    return "adaptive-rk23";
  }
  static const WVAdaptiveRKStageBufferLastUse *stageBufferLastUseRecords()
      noexcept;
  static std::size_t stageBufferLastUseRecordCount() noexcept;
  double nextStepSize() const noexcept override { return nextStepSize_; }
  std::size_t persistentBytes() const noexcept override;

private:
  class Workspace;
  WVKernelStatus ensureWorkspace(const WVMutableIntegrationState &state);
  WVIntegrationSystem &system_;
  WVAdaptiveRK23Options options_;
  std::unique_ptr<WVIntegrationErrorPolicy> errorPolicy_;
  Workspace *workspace_ = nullptr;
  WVAcceptedStep acceptedStep_;
  mutable WVIntegratorMetrics metrics_;
  std::vector<WVAdaptiveRK23StepDiagnostic> stepDiagnostics_;
  std::vector<std::uint64_t> toleranceComponentHashes_;
  std::uint64_t toleranceHash_ = 0;
  double nextStepSize_ = 0.0;
  bool hasAcceptedStep_ = false;
  bool fsalAvailable_ = false;
  bool stepping_ = false;
  mutable bool evaluatingDenseOutput_ = false;

public:
  ~WVAdaptiveRK23();
  WVAdaptiveRK23(const WVAdaptiveRK23 &) = delete;
  WVAdaptiveRK23 &operator=(const WVAdaptiveRK23 &) = delete;
};

struct WVAdaptiveRK45Options {
  double relativeTolerance = 1e-3;
  double absoluteToleranceScale = 1e-6;
  double safetyFactor = 0.8;
  double rejectionFloorFactor = 0.1;
  double repeatedRejectionFactor = 0.5;
  double maximumStepFactor = 5.0;
  double maximumStepSize = std::numeric_limits<double>::infinity();
  std::size_t maximumRecordedStepDiagnostics = 512;
  bool retainDenseOutput = true;
};

using WVAdaptiveRK45StepDiagnostic = WVAdaptiveRKStepDiagnostic;

class WVAdaptiveRK45 final : public WVTimeIntegrator,
                             public WVDenseOutput {
public:
  explicit WVAdaptiveRK45(WVIntegrationSystem &system,
                          WVAdaptiveRK45Options options = {});
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return system_.stateLayout();
  }
  WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) override;
  WVKernelStatus step(WVMutableIntegrationState &state,
                      double proposedStepSize) override;
  WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                               double finalTime,
                               double initialStepSize) override;
  WVKernelStatus evaluateDenseOutput(double time,
                                     WVMutableIntegrationState &output) const;
  WVKernelStatus
  evaluateState(double time, WVMutableIntegrationState &output) const override {
    return evaluateDenseOutput(time, output);
  }
  double initialTime() const noexcept override;
  double finalTime() const noexcept override;
  const WVAcceptedStep *lastAcceptedStep() const noexcept override;
  const WVIntegratorMetrics &metrics() const noexcept;
  const std::vector<WVAdaptiveRK45StepDiagnostic> &stepDiagnostics() const
      noexcept;
  std::uint64_t toleranceHash() const noexcept;
  const std::vector<std::uint64_t> &toleranceComponentHashes() const noexcept;
  bool stepDiagnosticsComplete() const noexcept;
  static constexpr const char *controllerIdentifier() noexcept {
    return "matlab-ode45-v1";
  }
  static constexpr const char *methodIdentifier() noexcept {
    return "adaptive-rk45";
  }
  static const WVAdaptiveRKStageBufferLastUse *stageBufferLastUseRecords()
      noexcept;
  static std::size_t stageBufferLastUseRecordCount() noexcept;
  double nextStepSize() const noexcept override;
  std::size_t persistentBytes() const noexcept override;

private:
  class Workspace;
  WVKernelStatus ensureWorkspace(const WVMutableIntegrationState &state);
  WVKernelStatus stepImplementation(WVMutableIntegrationState &state,
                                    double proposedStepSize,
                                    bool allowFinalStepStretch);
  WVIntegrationSystem &system_;
  WVAdaptiveRK45Options options_;
  std::unique_ptr<WVIntegrationErrorPolicy> errorPolicy_;
  Workspace *workspace_ = nullptr;
  WVAcceptedStep acceptedStep_;
  mutable WVIntegratorMetrics metrics_;
  std::vector<WVAdaptiveRK45StepDiagnostic> stepDiagnostics_;
  std::vector<std::uint64_t> toleranceComponentHashes_;
  std::uint64_t toleranceHash_ = 0;
  double nextStepSize_ = 0.0;
  bool hasAcceptedStep_ = false;
  bool fsalAvailable_ = false;
  bool stepping_ = false;
  mutable bool evaluatingDenseOutput_ = false;

public:
  ~WVAdaptiveRK45();
  WVAdaptiveRK45(const WVAdaptiveRK45 &) = delete;
  WVAdaptiveRK45 &operator=(const WVAdaptiveRK45 &) = delete;
};

struct WVAdaptiveRK78Options {
  double relativeTolerance = 1e-3;
  double absoluteToleranceScale = 1e-6;
  double safetyFactor = 0.8;
  double rejectionFloorFactor = 0.1;
  double repeatedRejectionFactor = 0.5;
  double maximumStepFactor = 5.0;
  double maximumStepSize = std::numeric_limits<double>::infinity();
  std::size_t maximumRecordedStepDiagnostics = 512;
  bool retainDenseOutput = false;
};

using WVAdaptiveRK78StepDiagnostic = WVAdaptiveRKStepDiagnostic;

// Verner's most-efficient Runge--Kutta 8(7) pair used by MATLAB ode78.
// Endpoint-only execution shares the base-step workspace but does not retain
// an initial-state view or allocate/evaluate the four extension-only stages.
class WVAdaptiveRK78 final : public WVTimeIntegrator,
                             public WVDenseOutput {
public:
  explicit WVAdaptiveRK78(WVIntegrationSystem &system,
                          WVAdaptiveRK78Options options = {});
  const WVIntegrationStateLayout &stateLayout() const noexcept override {
    return system_.stateLayout();
  }
  WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) override;
  WVKernelStatus step(WVMutableIntegrationState &state,
                      double proposedStepSize) override;
  WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                               double finalTime,
                               double initialStepSize) override;
  WVKernelStatus evaluateDenseOutput(
      double time, WVMutableIntegrationState &output) const;
  WVKernelStatus evaluateState(
      double time, WVMutableIntegrationState &output) const override {
    return evaluateDenseOutput(time, output);
  }
  double initialTime() const noexcept override;
  double finalTime() const noexcept override;
  const WVAcceptedStep *lastAcceptedStep() const noexcept override;
  const WVIntegratorMetrics &metrics() const noexcept;
  const std::vector<WVAdaptiveRK78StepDiagnostic> &stepDiagnostics() const
      noexcept;
  std::uint64_t toleranceHash() const noexcept;
  const std::vector<std::uint64_t> &toleranceComponentHashes() const noexcept;
  bool stepDiagnosticsComplete() const noexcept;
  static const char *controllerIdentifier() noexcept;
  static const char *methodIdentifier() noexcept;
  static const WVAdaptiveRKStageBufferLastUse *stageBufferLastUseRecords()
      noexcept;
  static std::size_t stageBufferLastUseRecordCount() noexcept;
  double nextStepSize() const noexcept override;
  std::size_t persistentBytes() const noexcept override;

private:
  class Workspace;
  class ContinuousExtensionWorkspace;
  WVKernelStatus ensureWorkspace(const WVMutableIntegrationState &state);
  WVKernelStatus ensureContinuousExtension() const;
  void releaseContinuousExtension() const noexcept;
  WVKernelStatus stepImplementation(WVMutableIntegrationState &state,
                                    double proposedStepSize,
                                    bool allowFinalStepStretch);
  WVIntegrationSystem &system_;
  WVAdaptiveRK78Options options_;
  std::unique_ptr<WVIntegrationErrorPolicy> errorPolicy_;
  Workspace *workspace_ = nullptr;
  mutable ContinuousExtensionWorkspace *continuousExtension_ = nullptr;
  WVAcceptedStep acceptedStep_;
  mutable WVIntegratorMetrics metrics_;
  std::vector<WVAdaptiveRK78StepDiagnostic> stepDiagnostics_;
  std::vector<std::uint64_t> toleranceComponentHashes_;
  std::uint64_t toleranceHash_ = 0;
  double nextStepSize_ = 0.0;
  bool hasAcceptedStep_ = false;
  bool derivativeReuseAvailable_ = false;
  bool stepping_ = false;
  mutable bool evaluatingDenseOutput_ = false;
  mutable bool continuousExtensionReady_ = false;

public:
  ~WVAdaptiveRK78();
  WVAdaptiveRK78(const WVAdaptiveRK78 &) = delete;
  WVAdaptiveRK78 &operator=(const WVAdaptiveRK78 &) = delete;
};

} // namespace wavevortex::runtime
