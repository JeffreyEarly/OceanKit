#pragma once
#include "WaveVortexKernel/WVSpectralOperators.hpp"
namespace wavevortex {
// Adapter factory; returns unsupportedOperation on non-Apple builds. Does not
// change process-wide BLAS settings or create application worker pools.
WVKernelStatus WVCreateAccelerateMatrixBackend(std::unique_ptr<WVVerticalMatrixBackend>&);
}
