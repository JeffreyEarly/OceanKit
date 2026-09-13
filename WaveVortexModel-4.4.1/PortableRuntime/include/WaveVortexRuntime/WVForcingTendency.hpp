#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"

#include <cstddef>

namespace wavevortex::runtime {

// Construction-resolved forcing index and caller-owned full-grid output.
// One entry requests all applicable spatial channels for one forcing instance:
// [Fu,Fv,Feta], [Fu,Fv,Fw,Feta], or [Fqgpv]. Numerical evaluation never resolves
// names. Sampling and individual channel selection belong to the field plan.
struct WVForcingTendencyOutput {
  std::size_t executionIndex = 0;
  WVRealFieldBundleView fields;
};

struct WVForcingTendencyMetrics {
  WVVariableEvaluationMetrics standaloneVariableEvaluation;
  std::size_t evaluationCount = 0;
  std::size_t forcingEvaluationCount = 0;
  std::size_t spatialProjectionCount = 0;
  std::size_t spectralReconstructionCount = 0;
  std::size_t workspaceLiveBytes = 0;
  std::size_t workspaceLastPeakBytes = 0;
  std::size_t workspaceHighWaterBytes = 0;
};

namespace detail {
class WVForcingDiagnosticWorkspace;
}

} // namespace wavevortex::runtime
