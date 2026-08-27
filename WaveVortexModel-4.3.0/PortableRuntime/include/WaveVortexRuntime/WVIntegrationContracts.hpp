#pragma once

#include "WaveVortexRuntime/WVIntegrationState.hpp"

#include <cstddef>
#include <limits>
#include <memory>

namespace wavevortex::runtime {

class WVFieldEvaluationService;

// Method-neutral absolute-error policy. Components are addressed by the
// integration layout's stable coefficient-family order, then each additional
// state block. Numerical methods never inspect transform or observer identity.
class WVIntegrationErrorPolicy {
public:
  virtual ~WVIntegrationErrorPolicy() = default;
  virtual std::size_t componentCount() const noexcept = 0;
  virtual std::size_t elementCount(std::size_t component) const noexcept = 0;
  virtual double absoluteTolerance(std::size_t component,
                                   std::size_t index) const noexcept = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
};

struct WVStateConstraintResult {
  WVKernelStatus status;
  std::size_t modifiedCoefficientCount = 0;
  bool fsalCompatible = true;

  explicit operator bool() const noexcept { return static_cast<bool>(status); }
};

// MATLAB WVModel.timeStepForCFL candidates produced by a resolved numerical
// system. Selection policy remains caller-owned and transform-neutral.
struct WVFixedTimeStepCandidates {
  double effectiveHorizontalGridResolution =
      std::numeric_limits<double>::infinity();
  double maximumHorizontalSpeed = 0.0;
  double horizontalAdvective = std::numeric_limits<double>::infinity();
  double verticalAdvective = std::numeric_limits<double>::infinity();
  double advective = std::numeric_limits<double>::infinity();
  double highestActiveWaveFrequency = 0.0;
  double oscillatory = std::numeric_limits<double>::infinity();
  std::size_t transientWorkspaceMaximumLiveBytes = 0;
  double evaluationSeconds = 0.0;
};

// Model boundary used by every numerical method. A successful RHS evaluation
// completely overwrites every coefficient and additional-state tendency.
class WVIntegrationSystem {
public:
  virtual ~WVIntegrationSystem() = default;
  virtual const WVIntegrationStateLayout &stateLayout() const noexcept = 0;
  virtual WVKernelStatus
  evaluateRightHandSide(const WVIntegrationState &state,
                        WVIntegrationFlux &rightHandSide) = 0;
  virtual WVStateConstraintResult
  enforceStateConstraints(WVMutableIntegrationState &state) = 0;
  virtual WVKernelStatus
  createErrorPolicy(double absoluteToleranceScale,
                    std::unique_ptr<WVIntegrationErrorPolicy> &policy) const = 0;
  virtual bool supportsFixedTimeStepSelection() const noexcept {
    return false;
  }
  virtual WVKernelStatus evaluateFixedTimeStepCandidates(
      const WVIntegrationState &, double,
      WVFixedTimeStepCandidates &) {
    return {WVKernelStatusCode::unsupportedOperation,
            "The numerical system does not provide fixed-step CFL "
            "candidates."};
  }
  // Optional transform-selected service consumed through the neutral model
  // boundary. Integrators do not depend on field evaluation.
  virtual WVFieldEvaluationService *fieldEvaluationService() noexcept {
    return nullptr;
  }
  virtual std::size_t persistentBytes() const noexcept { return 0; }
};

// Method-owned continuous extension over the most recent accepted interval.
class WVDenseOutput {
public:
  virtual ~WVDenseOutput() = default;
  virtual double initialTime() const noexcept = 0;
  virtual double finalTime() const noexcept = 0;
  virtual const WVIntegrationStateLayout &stateLayout() const noexcept = 0;
  virtual WVKernelStatus
  evaluateState(double time, WVMutableIntegrationState &output) const = 0;
};

struct WVAcceptedStep {
  struct MethodStatistics {
    std::size_t acceptedStepCount = 0;
    std::size_t rejectedStepCount = 0;
    std::size_t rightHandSideEvaluationCount = 0;
    double stepSize = 0.0;
    double proposedStepSize = 0.0;
    double nextStepSize = 0.0;
    double normalizedError = 0.0;
  };

  double initialTime = 0.0;
  double finalTime = 0.0;
  WVIntegrationState endpoint;
  MethodStatistics methodStatistics;
  // Valid only until the owning integrator advances or is prepared again.
  const WVDenseOutput *denseOutput = nullptr;
};

class WVTimeIntegrator {
public:
  virtual ~WVTimeIntegrator() = default;
  virtual const WVIntegrationStateLayout &stateLayout() const noexcept = 0;
  virtual WVKernelStatus
  prepareStateAfterRestart(WVMutableIntegrationState &state) = 0;
  virtual WVKernelStatus step(WVMutableIntegrationState &state,
                              double stepSize) = 0;
  virtual WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                                       double finalTime,
                                       double stepSize) = 0;
  virtual const WVAcceptedStep *lastAcceptedStep() const noexcept = 0;
  virtual double nextStepSize() const noexcept = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
};

} // namespace wavevortex::runtime
