#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexKernel/WVVariableExecutionOptions.hpp"

#include <cstddef>
#include <string_view>

namespace wavevortex::runtime {

enum class WVNativeMatrixBackend { scalar, accelerate };

struct WVNativeHostTopology {
    std::size_t logicalWorkers = 1;
    std::size_t performanceWorkers = 1;
};

struct WVNativeVariablePolicy {
    bool buildEnabled = false;
    bool variableTransform = false;
    bool compact = false;
    WVPersistedTransformKind transformKind = WVPersistedTransformKind::constantStratification;
    WVNativeMatrixBackend matrixBackend = WVNativeMatrixBackend::scalar;
    WVVariableExecutionOptions execution;
    std::size_t effectiveFFTThreads = 1;
    std::string_view selection = "build-disabled";
};

WVNativeHostTopology nativeHostTopology() noexcept;

WVNativeVariablePolicy selectNativeVariablePolicy(
    bool buildEnabled,
    WVPersistedTransformKind transformKind,
    std::string_view provider,
    std::size_t fftThreads,
    bool hasRequestedThreads,
    WVNativeHostTopology topology) noexcept;

const char* nativeMatrixBackendIdentifier(WVNativeMatrixBackend) noexcept;
const char* nativeTransformKindIdentifier(WVPersistedTransformKind) noexcept;
const char* nativeSpectralScheduleIdentifier(WVVariableSpectralSchedule) noexcept;
const char* nativeHorizontalScheduleIdentifier(WVRetainedHorizontalSchedule) noexcept;

} // namespace wavevortex::runtime
