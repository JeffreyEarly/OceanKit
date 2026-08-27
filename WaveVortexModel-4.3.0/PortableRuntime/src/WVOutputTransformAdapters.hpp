#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"

namespace wavevortex::runtime {

class WVPortableObserverDescriptor;

namespace detail {

// Compatibility and persistence details for output orchestration live behind
// this named boundary. Generic scheduling and delivery code consumes only the
// resolved integration layout and ordered coefficient-family views.
WVKernelStatus createLegacyOutputStateLayout(
    const WVPortableObserverDescriptor &descriptor,
    WVIntegrationStateLayout &layout);

void bindLegacyOutputStateView(
    const WVIntegrationStateLayout &layout,
    WVCoefficientStateStorage &storage,
    WVMutableIntegrationState &state) noexcept;

WVKernelStatus validateOutputCheckpointTemplate(
    const WVIntegrationStateLayout &layout, const WVCheckpoint &checkpoint);

WVKernelStatus stageOutputCheckpointState(
    const WVIntegrationStateLayout &layout, WVCheckpoint &checkpoint,
    const WVIntegrationState &state);

} // namespace detail
} // namespace wavevortex::runtime
