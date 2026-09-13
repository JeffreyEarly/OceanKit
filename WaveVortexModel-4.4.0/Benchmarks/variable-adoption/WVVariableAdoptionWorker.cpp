#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#include "nlohmann/json.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {
template<class Status> void require(Status status) {
    if (!status) throw std::runtime_error(status.message);
}
std::size_t count(const char* value) {
    const std::string text(value);
    if (text.empty() || text.find_first_not_of("0123456789")!=std::string::npos)
        throw std::runtime_error("Expected a nonnegative integer.");
    return std::stoull(text);
}
template<class Kernel, class Execute>
json measure(Kernel& kernel,Execute execute,std::size_t warmups,std::size_t samples,std::chrono::steady_clock::time_point preparationStart) {
    const double preparationSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-preparationStart).count();
    for (std::size_t i=0;i<warmups;++i) require(execute());
    kernel.resetMetrics();
    const auto bytes=kernel.persistentBytes();
    std::vector<double> times(samples);
    for (auto& time:times) {
        const auto start=std::chrono::steady_clock::now();
        const auto status=execute();
        time=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        // Journal completed timing before validation, outside the measured
        // interval, so a later failure cannot erase earlier samples.
        std::cerr<<json{{"sampleSeconds",time},{"success",bool(status)}}.dump()<<std::endl;
        require(status);
    }
    if (kernel.persistentBytes()!=bytes) throw std::runtime_error("Persistent storage changed during execution.");
    const auto& s=kernel.storage();
    return {{"preparationSeconds",preparationSeconds},{"matrixBackendIdentifier",kernel.matrixBackendIdentifier()},{"compactSplitViews",kernel.executionOptions().usesCompactSplitViews()},{"streamedNonlinear",kernel.executionOptions().streamedNonlinear},{"samplesSeconds",times},{"horizontalSchedule",kernel.horizontalScheduleIdentifier()},
        {"persistentBytes",bytes},{"realScratchBytes",s.realScratchBytes},
        {"spectralScratchBytes",s.spectralScratchBytes},{"preparedBytes",s.preparedBytes},
        {"workspaceBytes",s.workspaceBytes},{"providerBytesLowerBound",s.providerBytesLowerBound},
        {"planBytesLowerBound",s.planBytesLowerBound}};
}
}
int main(int argc,char** argv) {
    try {
        if (argc<7 || argc>9) throw std::runtime_error("Usage: worker INPUT frozen|pruned|streamed|pruned-streamed|compact WORKERS WARMUPS SAMPLES OUTPUT [POINTWISE_WORKERS] [VERTICAL_GROUP_WORKERS]");
        const std::string selection=argv[2];
        if (selection!="frozen" && selection!="pruned" && selection!="streamed" && selection!="pruned-streamed" && selection!="compact")
            throw std::runtime_error("Unknown execution selection.");
        const auto workers=count(argv[3]),warmups=count(argv[4]),samples=count(argv[5]);
        if (!workers || !samples) throw std::runtime_error("Workers and samples must be positive.");
        if (std::filesystem::exists(argv[6])) throw std::runtime_error("Refusing to overwrite an existing payload.");
        WVVariableExecutionOptions options;
        if (argc>=8) options.pointwiseWorkers=count(argv[7]);
        if (argc==9) options.verticalGroupWorkers=count(argv[8]);
        if (!options.pointwiseWorkers || !options.verticalGroupWorkers)
            throw std::runtime_error("Pointwise and vertical group workers must be positive.");
        if (selection=="frozen" && (options.pointwiseWorkers!=1 || options.verticalGroupWorkers!=1))
            throw std::runtime_error("The frozen selection requires one pointwise and vertical group worker.");
        if (selection=="frozen") options.inertialOnlyProjection=false;
        if (selection=="pruned" || selection=="pruned-streamed" || selection=="compact") {
            options.horizontalSchedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
            options.horizontalWorkers=workers;
        }
        options.streamedNonlinear=selection=="streamed" || selection=="pruned-streamed" || selection=="compact";
        if (selection=="compact") options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
        std::shared_ptr<const WVExtensionCatalog> catalog; require(makeBuiltInExtensionCatalog(catalog));
        WVCheckpoint checkpoint; require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
        if (!checkpoint.stratifiedModalSource) throw std::runtime_error("A variable-stratification source is required.");
        if (checkpoint.forcingSchedule.entries.size()!=1 || checkpoint.forcingSchedule.entries.front().typeIdentifier!="WVNonlinearAdvection")
            throw std::runtime_error("Fixture must contain only nonlinear advection.");
        std::unique_ptr<WVFFTEngine> fft;
        // Keep the frozen FFT provider policy fixed. Pruned workers belong to
        // its prepared outer executor; no process-wide BLAS settings change.
        require(WVFFTWEngine::create(1,fft));
        const auto& g=checkpoint.stratifiedModalSource->geometry();
        const WVShape2D shape{g.Nj,g.Nkl}; const auto S=g.Nj*g.Nkl;
        const bool qg=g.transformClass=="WVTransformStratifiedQG";
        std::vector<WVComplex64> fp(qg ? 0 : S),fm(qg ? 0 : S),f0(S);
        const auto coefficient=[&](const std::string& name) -> WVComplexConstView {
            for (const auto& family:checkpoint.transformState.coefficientFamilies) if (family.identifier==name) {
                if (family.values.size()!=S || family.spectralDimensions!=std::vector<std::size_t>{g.Nj,g.Nkl})
                    throw std::runtime_error("Fixture coefficient shape differs from the modal record.");
                return {family.values.data(),shape};
            }
            throw std::runtime_error("Fixture is missing coefficient family "+name);
        };
        const WVState state{checkpoint.transformState.t,checkpoint.transformState.t0,
            {qg ? WVComplexConstView{} : coefficient("Ap"),qg ? WVComplexConstView{} : coefficient("Am"),coefficient("A0")}};
        WVFlux flux{{fp.data(),shape},{fm.data(),shape},{f0.data(),shape}};
        json result;
        const auto preparationStart=std::chrono::steady_clock::now();
        if (qg) {
            std::unique_ptr<WVTransformStratifiedQGKernel> kernel;
            require(WVTransformStratifiedQGKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&] { return kernel->nonlinearFlux(state.coefficients.A0,flux.F0); },warmups,samples,preparationStart);
        } else if (g.transformClass=="WVTransformHydrostatic") {
            std::unique_ptr<WVTransformHydrostaticKernel> kernel;
            require(WVTransformHydrostaticKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&] { return kernel->nonlinearFlux(state,flux); },warmups,samples,preparationStart);
        } else if (g.transformClass=="WVTransformBoussinesq") {
            std::unique_ptr<WVTransformBoussinesqKernel> kernel;
            require(WVTransformBoussinesqKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&] { return kernel->nonlinearFlux(state,flux); },warmups,samples,preparationStart);
            if (kernel->metrics().verticalOperatorExecutionCount%samples ||
                kernel->metrics().verticalMatrixGroupExecutionCount%samples)
                throw std::runtime_error("Boussinesq vertical execution counts differ between timed RHS samples.");
            result["verticalOperatorExecutionCount"]=kernel->metrics().verticalOperatorExecutionCount;
            result["verticalMatrixGroupExecutionCount"]=kernel->metrics().verticalMatrixGroupExecutionCount;
            result["verticalOperatorsPerRHS"]=kernel->metrics().verticalOperatorExecutionCount/samples;
            result["verticalMatrixGroupsPerRHS"]=kernel->metrics().verticalMatrixGroupExecutionCount/samples;
        } else throw std::runtime_error("Unsupported transform family.");
        rusage usage{}; if (getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("Cannot measure peak RSS.");
        result["peakRSSBytes"]=static_cast<std::size_t>(usage.ru_maxrss);
        result["schema"]="wvm-variable-screening-v1"; result["family"]=g.transformClass;
        result["pointwiseWorkers"]=options.pointwiseWorkers;
        result["verticalGroupWorkers"]=options.verticalGroupWorkers;
        result["selection"]=selection; result["horizontalWorkers"]=options.horizontalWorkers;
        result["grid"]={g.Nx,g.Ny,g.Nz}; result["warmups"]=warmups;
        result["matrixBackend"]="accelerate";
        const auto identity=WVFFTWEngine::linkedLibraries();
        result["provider"]={{"version",identity.version},{"baseLibrary",identity.baseLibrary},{"threadLibrary",identity.threadLibrary},{"fftThreads",1}};
        std::ofstream output(argv[6],std::ios::binary);
        for (const auto* a:{&fp,&fm,&f0}) if (!a->empty())
            output.write(reinterpret_cast<const char*>(a->data()),static_cast<std::streamsize>(a->size()*sizeof(WVComplex64)));
        output.close(); if (!output) throw std::runtime_error("Cannot write flux payload.");
        std::cout<<result.dump()<<'\n';
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
