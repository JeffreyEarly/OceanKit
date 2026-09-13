#pragma once

#include "WaveVortexKernel/WVVariableExecutionOptions.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"

#include <functional>
#include <memory>

namespace wavevortex::runtime {

// Author/runtime construction services. The default keeps the portable scalar
// matrix backend and established variable execution schedule.
struct WVVariableKernelServices {
  std::function<WVKernelStatus(std::unique_ptr<WVVerticalMatrixBackend>&)>
      matrixBackendFactory = WVCreateScalarMatrixBackend;
  WVVariableExecutionOptions execution;
  WVVariableEvaluationPolicy variableEvaluationPolicy=WVVariableEvaluationPolicy::reuse;
};

} // namespace wavevortex::runtime
