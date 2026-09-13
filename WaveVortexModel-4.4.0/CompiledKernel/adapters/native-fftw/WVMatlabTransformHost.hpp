#pragma once
#include "mex.h"
#include <cstddef>
#include <string>

// Dispatches the transform-neutral, call-scoped MATLAB bridge. The legacy
// constant preview gateway remains available while family adapters migrate.
bool WVDispatchMatlabTransform(const std::string& command, int nlhs,
    mxArray* plhs[], int nrhs, const mxArray* prhs[]);
std::size_t WVMatlabTransformCount() noexcept;
void WVCleanupMatlabTransforms() noexcept;
