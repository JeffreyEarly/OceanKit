#include "WVRunnerVariablePolicy.hpp"

#include <utility>

namespace wavevortex::runtime::cli {

WVRunnerHostTopology runnerHostTopology() noexcept {
    return nativeHostTopology();
}

WVRunnerVariablePolicy selectRunnerVariablePolicy(
    bool buildEnabled,
    WVPersistedTransformKind transformKind,
    std::string_view provider,
    std::size_t fftThreads,
    bool hasRequestedThreads,
    WVRunnerHostTopology topology) noexcept {
    return selectNativeVariablePolicy(buildEnabled,transformKind,provider,fftThreads,hasRequestedThreads,topology);
}

const char* runnerMatrixBackendIdentifier(WVRunnerMatrixBackend backend) noexcept {
    return nativeMatrixBackendIdentifier(backend);
}

const char* runnerTransformKindIdentifier(WVPersistedTransformKind kind) noexcept {
    return nativeTransformKindIdentifier(kind);
}

const char* runnerSpectralScheduleIdentifier(WVVariableSpectralSchedule schedule) noexcept {
    return nativeSpectralScheduleIdentifier(schedule);
}

const char* runnerHorizontalScheduleIdentifier(WVRetainedHorizontalSchedule schedule) noexcept {
    return nativeHorizontalScheduleIdentifier(schedule);
}

} // namespace wavevortex::runtime::cli
