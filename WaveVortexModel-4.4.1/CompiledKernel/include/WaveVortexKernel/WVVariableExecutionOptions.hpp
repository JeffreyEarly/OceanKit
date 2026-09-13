#pragma once
#include "WVSpectralOperators.hpp"

namespace wavevortex {
// Optional event-owned bridge between a nonlinear transform and typed output
// derivative storage. An empty lookup result is a cache miss. Capture may be a
// no-op when the output plan proves the derivative cannot be consumed again.
struct WVStateDerivativeAccess {
    void* context = nullptr;
    WVKernelStatus (*lookup)(void*,std::size_t field,std::size_t derivative,
        WVRealVolumeConstView&) = nullptr;
    WVKernelStatus (*capture)(void*,std::size_t field,std::size_t derivative,
        WVRealVolumeConstView) = nullptr;
};

enum class WVVariableSpectralSchedule {
    establishedInterleaved,
    compactSplitFusedViews
};

// Execution policy only: canonical coefficients, modal records and scientific
// contracts are independent of execution scheduling. Kernel constructors retain
// the established path; eligible native runners may supply a qualified schedule.
struct WVVariableExecutionOptions {
    WVRetainedHorizontalSchedule horizontalSchedule = WVRetainedHorizontalSchedule::fullFFT;
    std::size_t horizontalWorkers = 1;
    bool streamedNonlinear = false;
    // The candidate compact variable path owns split-complex operator storage
    // and passes direct family views between vertical and horizontal services.
    // This remains explicit for direct kernel callers.
    WVVariableSpectralSchedule spectralSchedule = WVVariableSpectralSchedule::establishedInterleaved;
    // Prepared pointwise workers partition independent physical-grid cells.
    // The default preserves established serial arithmetic for direct kernel callers.
    std::size_t pointwiseWorkers = 1;
    // Reuse modal and horizontal-spectrum preparation within one evaluation.
    // Disabling this retains the independent reconstruction path for parity
    // qualification; it is not a separate scientific or memory policy.
    bool sharedFieldGradients = true;
    // Consume completed inverse planes in nonlinear advection, retaining the
    // ordinary derivative output and success-only evaluator capture.
    bool fusedDerivativeAdvection = false;
    // Fuse first physical-field production and nonlinear advection in short
    // native horizontal tiles. Direct callers opt in; other paths remain valid.
    bool tiledNonlinear = false;
    // Boussinesq projection needs the grouped wave-F result only at the exact
    // zero horizontal mode. Disable for a full-coverage qualification oracle.
    bool inertialOnlyProjection = true;
    // Persistent outer workers for independent prepared vertical matrix groups.
    // One preserves serial backend execution for direct kernel callers.
    std::size_t verticalGroupWorkers = 1;
    bool usesCompactSplitViews() const noexcept {
        return spectralSchedule==WVVariableSpectralSchedule::compactSplitFusedViews;
    }
};
} // namespace wavevortex
