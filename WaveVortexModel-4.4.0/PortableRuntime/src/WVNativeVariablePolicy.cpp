#include "WaveVortexRuntime/WVNativeVariablePolicy.hpp"

#include <algorithm>
#include <cstdint>
#include <thread>

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace wavevortex::runtime {
namespace {

constexpr std::size_t maximumHorizontalWorkers=12;
constexpr std::size_t maximumPointwiseWorkers=8;

bool variableTransform(WVPersistedTransformKind kind) noexcept {
    return kind==WVPersistedTransformKind::stratifiedQG ||
        kind==WVPersistedTransformKind::hydrostatic ||
        kind==WVPersistedTransformKind::boussinesq;
}

#if defined(__APPLE__)
std::size_t sysctlCount(const char* name) noexcept {
    std::uint32_t value=0;
    std::size_t size=sizeof(value);
    return ::sysctlbyname(name,&value,&size,nullptr,0)==0 && size==sizeof(value) && value ? value : 0;
}
#endif

} // namespace

WVNativeHostTopology nativeHostTopology() noexcept {
    const auto logical=std::max<std::size_t>(1,std::thread::hardware_concurrency());
    std::size_t performance=logical;
#if defined(__APPLE__)
    performance=sysctlCount("hw.perflevel0.physicalcpu");
    if (!performance) performance=sysctlCount("hw.physicalcpu");
    if (!performance) performance=logical;
#endif
    return {logical,std::min(performance,logical)};
}

WVNativeVariablePolicy selectNativeVariablePolicy(
    bool buildEnabled,
    WVPersistedTransformKind transformKind,
    std::string_view provider,
    std::size_t fftThreads,
    bool hasRequestedThreads,
    WVNativeHostTopology topology) noexcept {
    WVNativeVariablePolicy policy;
    policy.buildEnabled=buildEnabled;
    policy.transformKind=transformKind;
    policy.variableTransform=variableTransform(transformKind);
    policy.effectiveFFTThreads=std::max<std::size_t>(1,fftThreads);
    if (!buildEnabled) return policy;
    if (!policy.variableTransform) {
        policy.selection="non-variable-transform";
        return policy;
    }
    if (provider=="reference") {
        policy.selection="reference-provider";
        return policy;
    }
    if (provider!="native-fftw") {
        policy.selection="unsupported-provider";
        return policy;
    }

    policy.matrixBackend=WVNativeMatrixBackend::accelerate;
    // --threads continues to mean FFTW-internal threads. An explicit parallel
    // FFT request keeps the established variable schedule and is reported as
    // such; the compact topology owns parallelism outside FFTW.
    if (hasRequestedThreads && fftThreads>1) {
        policy.selection="established-explicit-fft-threads";
        return policy;
    }

    const auto performanceWorkers=std::max<std::size_t>(1,
        std::min({topology.performanceWorkers,
            std::max<std::size_t>(1,topology.logicalWorkers),maximumHorizontalWorkers}));
    policy.compact=true;
    policy.selection="compact-native-accelerate";
    policy.effectiveFFTThreads=1;
    policy.execution.horizontalSchedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
    policy.execution.horizontalWorkers=performanceWorkers;
    policy.execution.streamedNonlinear=true;
    policy.execution.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
    policy.execution.pointwiseWorkers=std::min(maximumPointwiseWorkers,performanceWorkers);
    if (transformKind==WVPersistedTransformKind::boussinesq)
        policy.execution.verticalGroupWorkers=std::min<std::size_t>(8,performanceWorkers);
    policy.execution.fusedDerivativeAdvection=
        transformKind==WVPersistedTransformKind::hydrostatic ||
        transformKind==WVPersistedTransformKind::boussinesq;
    policy.execution.tiledNonlinear=policy.execution.fusedDerivativeAdvection;
    return policy;
}

const char* nativeMatrixBackendIdentifier(WVNativeMatrixBackend backend) noexcept {
    return backend==WVNativeMatrixBackend::accelerate ? "accelerate" : "scalar";
}

const char* nativeTransformKindIdentifier(WVPersistedTransformKind kind) noexcept {
    switch (kind) {
        case WVPersistedTransformKind::constantStratification: return "constant-stratification";
        case WVPersistedTransformKind::barotropicQG: return "barotropic-qg";
        case WVPersistedTransformKind::stratifiedQG: return "stratified-qg";
        case WVPersistedTransformKind::hydrostatic: return "hydrostatic";
        case WVPersistedTransformKind::boussinesq: return "boussinesq";
    }
    return "unknown";
}

const char* nativeSpectralScheduleIdentifier(WVVariableSpectralSchedule schedule) noexcept {
    return schedule==WVVariableSpectralSchedule::compactSplitFusedViews ?
        "compact-split-fused-views" : "established-interleaved";
}

const char* nativeHorizontalScheduleIdentifier(WVRetainedHorizontalSchedule schedule) noexcept {
    return schedule==WVRetainedHorizontalSchedule::streamingPrunedTile16 ?
        "streaming-pruned-tile16" : "full-fft";
}

} // namespace wavevortex::runtime
