#include "WaveVortexRuntime/WVModel.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {
template<class T> void require(T status) { if (!status) throw std::runtime_error(status.message); }
void check(bool condition,const char* message) { if(!condition) throw std::runtime_error(message); }
std::unique_ptr<WVFFTEngine> provider(const std::string& name) {
  if(name=="reference") return std::make_unique<WVReferenceFFTEngine>();
#if WV_TEST_NATIVE_FFTW
  if(name=="native-fftw") { std::unique_ptr<WVFFTEngine> engine; require(WVFFTWEngine::create(1,engine)); return engine; }
#endif
  throw std::runtime_error("Requested lifecycle provider unavailable.");
}
}
int main(int argc,char** argv) {
  try {
    check(argc==4,"Usage: WVStratifiedLifecycleProbe source.nc report.json reference|native-fftw");
    std::shared_ptr<const WVExtensionCatalog> catalog;require(makeBuiltInExtensionCatalog(catalog));
    std::weak_ptr<const WVStratifiedModalRecord> sourceOwner;
    std::size_t maximumGrowth=0,maximumAllocations=0,completed=0;
    json lifecycles=json::array();
    json grid;
    WVPersistedTransformKind transformKind=WVPersistedTransformKind::stratifiedQG;
    std::size_t modes=0;
    const std::filesystem::path path=std::string(argv[2])+".nc";
    struct Cleanup { std::filesystem::path path; ~Cleanup(){std::error_code error;std::filesystem::remove(path,error);} } cleanup{path};
    for(std::size_t cycle=0;cycle<6;++cycle) {
      std::filesystem::copy_file(argv[1],path,std::filesystem::copy_options::overwrite_existing);
      {
        const auto setupStart=std::chrono::steady_clock::now();
        WVModel model;WVModelState state;
        WVModelOutputRequest output;output.finalTime=617;
        require(WVModel::createFromModelOutputFiles(catalog,{path.string()},output,provider(argv[3]),{},model,state));
        const double setupSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-setupStart).count();
        transformKind=state.checkpoint().transformKind;
        sourceOwner=state.checkpoint().stratifiedModalSource;
        const auto& geometry=state.checkpoint().stratifiedModalSource->geometry();
        grid={geometry.Nx,geometry.Ny,geometry.Nz};
        modes=geometry.Nj;
        require(model.step(state,.125));
        const auto before=model.metrics(&state);
        allocationProbe::calls=0;allocationProbe::counting=true;
        const auto start=std::chrono::steady_clock::now();
        WVKernelStatus status=WVKernelStatus::ok();
        for(int step=0;step<16 && status;++step) status=model.step(state,.125);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        allocationProbe::counting=false;require(status);
        const auto allocations=allocationProbe::calls.load();
        const auto after=model.metrics(&state);
        const auto view=state.constView();
        const auto finiteComplex=[](WVComplex64 value){return std::isfinite(value.real) && std::isfinite(value.imag);};
        for(std::size_t family=0;family<view.coefficientFamilyCount;++family) {
          const auto& values=view.coefficientFamilies[family];
          check(std::all_of(values.data,values.data+values.layout->elementCount,finiteComplex),"Nonfinite lifecycle coefficients.");
        }
        for(std::size_t block=0;block<view.additionalBlockCount;++block) {
          const auto& values=view.additionalBlocks[block];
          if(values.realData) check(std::all_of(values.realData,values.realData+values.layout->elementCount,[](double value){return std::isfinite(value);}),"Nonfinite lifecycle observer state.");
          if(values.complexData) check(std::all_of(values.complexData,values.complexData+values.layout->elementCount,finiteComplex),"Nonfinite lifecycle observer state.");
        }
        const auto retained=[](const WVModelMetrics& m){return m.integrationSystemPersistentBytes+m.integratorPersistentBytes+m.statePersistentBytes+m.outputPersistentBytes;};
        check(retained(after)>=retained(before),"Prepared retained accounting unexpectedly shrank.");
        const auto growth=retained(after)-retained(before);
        maximumGrowth=std::max(maximumGrowth,growth);maximumAllocations=std::max(maximumAllocations,allocations);
        lifecycles.push_back({{"setupSeconds",setupSeconds},{"retainedBytesBefore",retained(before)},{"retainedBytesAfter",retained(after)},{"preparedStepAllocations",allocations},{"sampleSteps",16},{"sampleSeconds",seconds},{"rhsEvaluations",after.forcing.evaluationCount-before.forcing.evaluationCount},{"fieldReconstructions",after.forcing.physicalFieldReconstructionCount-before.forcing.physicalFieldReconstructionCount},{"spatialProjections",after.forcing.spatialTendencyProjectionCount-before.forcing.spatialTendencyProjectionCount},{"scientificBytes",state.checkpoint().stratifiedModalSource->persistentBytes()}});
        require(model.closeOutput());
      }
      check(sourceOwner.expired(),"Scientific modal owner survived model/state destruction.");
      ++completed;
    }
    json report={{"grid",grid},{"Nj",modes},{"finalStateFinite",true},{"schemaIdentifier",transformKind==WVPersistedTransformKind::boussinesq ? "wave-vortex-boussinesq-lifecycle-v1" : (transformKind==WVPersistedTransformKind::hydrostatic ? "wave-vortex-hydrostatic-lifecycle-v1" : "wave-vortex-sqg-lifecycle-v1")},{"provider",argv[3]},{"completedLifecycles",completed},{"scientificOwnersReleased",sourceOwner.expired()},{"retainedGrowthBytes",maximumGrowth},{"preparedStepAllocations",maximumAllocations},{"measurements",lifecycles},{"accountingScope","reported C++ capacities; excludes allocator metadata and opaque provider/NetCDF internals"}};
    std::ofstream out(argv[2]);out<<report.dump(2)<<'\n';check(static_cast<bool>(out),"Unable to write lifecycle evidence.");
  } catch(const std::exception& e) {allocationProbe::counting=false;std::cerr<<e.what()<<'\n';return 1;}
}
