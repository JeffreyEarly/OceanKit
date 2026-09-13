#include "WVRunnerVariablePolicy.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::cli;

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void established(const WVRunnerVariablePolicy& policy,const char* message) {
    require(!policy.compact,message);
    require(policy.execution.horizontalSchedule==WVRetainedHorizontalSchedule::fullFFT,message);
    require(policy.execution.horizontalWorkers==1,message);
    require(!policy.execution.streamedNonlinear,message);
    require(policy.execution.spectralSchedule==WVVariableSpectralSchedule::establishedInterleaved,message);
    require(policy.execution.pointwiseWorkers==1,message);
    require(policy.execution.verticalGroupWorkers==1,message);
    require(!policy.execution.fusedDerivativeAdvection,message);
    require(!policy.execution.tiledNonlinear,message);
}

} // namespace

int main() {
    try {
        const WVRunnerHostTopology host{18,12};
        auto policy=selectRunnerVariablePolicy(false,WVPersistedTransformKind::boussinesq,
            "native-fftw",18,false,host);
        established(policy,"Disabled policy changed variable execution");
        require(policy.selection=="build-disabled" && policy.effectiveFFTThreads==18 &&
            policy.matrixBackend==WVRunnerMatrixBackend::scalar,"Disabled policy report differs");

        for (const auto kind:{WVPersistedTransformKind::constantStratification,
                WVPersistedTransformKind::barotropicQG}) {
            policy=selectRunnerVariablePolicy(true,kind,"native-fftw",18,false,host);
            established(policy,"Variable policy changed a non-variable transform");
            require(policy.selection=="non-variable-transform" && policy.effectiveFFTThreads==18 &&
                policy.matrixBackend==WVRunnerMatrixBackend::scalar,"Non-variable selection report differs");
        }

        for (const auto kind:{WVPersistedTransformKind::stratifiedQG,
                WVPersistedTransformKind::hydrostatic,WVPersistedTransformKind::boussinesq}) {
            policy=selectRunnerVariablePolicy(true,kind,"native-fftw",18,false,host);
            require(policy.compact && policy.variableTransform &&
                policy.selection=="compact-native-accelerate" &&
                policy.matrixBackend==WVRunnerMatrixBackend::accelerate &&
                policy.effectiveFFTThreads==1,"Native compact selection differs");
            require(policy.execution.horizontalSchedule==WVRetainedHorizontalSchedule::streamingPrunedTile16 &&
                policy.execution.horizontalWorkers==12 && policy.execution.streamedNonlinear &&
                policy.execution.spectralSchedule==WVVariableSpectralSchedule::compactSplitFusedViews &&
                policy.execution.pointwiseWorkers==8 &&
                policy.execution.verticalGroupWorkers==
                    (kind==WVPersistedTransformKind::boussinesq ? 8 : 1) &&
                policy.execution.tiledNonlinear==policy.execution.fusedDerivativeAdvection &&
                policy.execution.fusedDerivativeAdvection==
                    (kind==WVPersistedTransformKind::hydrostatic ||
                     kind==WVPersistedTransformKind::boussinesq),
                "Native compact topology differs");
        }

        for (const auto kind:{WVPersistedTransformKind::stratifiedQG,
                WVPersistedTransformKind::hydrostatic,WVPersistedTransformKind::boussinesq}) {
            policy=selectRunnerVariablePolicy(true,kind,"native-fftw",4,true,host);
            established(policy,"Explicit FFT threads did not retain established execution");
            require(policy.selection=="established-explicit-fft-threads" &&
                policy.matrixBackend==WVRunnerMatrixBackend::accelerate &&
                policy.effectiveFFTThreads==4,"Explicit FFT selection was silently changed");
        }

        policy=selectRunnerVariablePolicy(true,WVPersistedTransformKind::stratifiedQG,
            "reference",1,false,host);
        established(policy,"Reference provider did not retain portable execution");
        require(policy.selection=="reference-provider" &&
            policy.matrixBackend==WVRunnerMatrixBackend::scalar &&
            policy.effectiveFFTThreads==1,"Reference provider selection report differs");

        policy=selectRunnerVariablePolicy(true,WVPersistedTransformKind::boussinesq,
            "native-fftw",8,false,{4,2});
        require(policy.execution.horizontalWorkers==2 && policy.execution.pointwiseWorkers==2 &&
            policy.execution.verticalGroupWorkers==2 && policy.execution.fusedDerivativeAdvection,
            "Compact topology exceeded host performance workers");
        policy=selectRunnerVariablePolicy(true,WVPersistedTransformKind::boussinesq,
            "native-fftw",24,false,{24,16});
        require(policy.execution.horizontalWorkers==12 && policy.execution.pointwiseWorkers==8 &&
            policy.execution.verticalGroupWorkers==8 && policy.execution.fusedDerivativeAdvection,
            "Compact topology exceeded its calibrated worker bounds");
        require(std::string(runnerMatrixBackendIdentifier(WVRunnerMatrixBackend::accelerate))=="accelerate" &&
            std::string(runnerTransformKindIdentifier(WVPersistedTransformKind::boussinesq))=="boussinesq" &&
            std::string(runnerSpectralScheduleIdentifier(WVVariableSpectralSchedule::compactSplitFusedViews))=="compact-split-fused-views" &&
            std::string(runnerHorizontalScheduleIdentifier(WVRetainedHorizontalSchedule::streamingPrunedTile16))=="streaming-pruned-tile16",
            "Policy report identifiers differ");
        std::cout<<"Runner variable policy tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
