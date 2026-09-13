#pragma once

#include "WaveVortexRuntime/WVNativeVariablePolicy.hpp"

#include <cstddef>
#include <string_view>

namespace wavevortex::runtime::cli {

using WVRunnerMatrixBackend = WVNativeMatrixBackend;
using WVRunnerHostTopology = WVNativeHostTopology;
using WVRunnerVariablePolicy = WVNativeVariablePolicy;

WVRunnerHostTopology runnerHostTopology() noexcept;

WVRunnerVariablePolicy selectRunnerVariablePolicy(
    bool buildEnabled,
    WVPersistedTransformKind transformKind,
    std::string_view provider,
    std::size_t fftThreads,
    bool hasRequestedThreads,
    WVRunnerHostTopology topology) noexcept;

const char* runnerMatrixBackendIdentifier(WVRunnerMatrixBackend) noexcept;
const char* runnerTransformKindIdentifier(WVPersistedTransformKind) noexcept;
const char* runnerSpectralScheduleIdentifier(WVVariableSpectralSchedule) noexcept;
const char* runnerHorizontalScheduleIdentifier(WVRetainedHorizontalSchedule) noexcept;

} // namespace wavevortex::runtime::cli
