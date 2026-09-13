#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
#include "WVFieldEvaluationEventWorkspace.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <fstream>
#include <iostream>
#include <cstring>
using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {template<class T> void require(T status) {if(!status) throw std::runtime_error(status.message);}}
int main(int argc,char** argv) {
  try {
    if(argc!=5) throw std::runtime_error("Usage: WVDensityEventDump model.nc request.json output.json reference|native");
    WVExtensionCatalogBuilder builder;require(addBuiltInExtensions(builder));
    std::shared_ptr<const WVExtensionCatalog> catalog;require(builder.freeze(catalog));
    WVCheckpoint checkpoint;require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
    std::ifstream requestFile(argv[2]);json request;requestFile>>request;
    std::unique_ptr<WVFFTEngine> fft;
    if(std::string(argv[4])=="reference") fft=std::make_unique<WVReferenceFFTEngine>();
    else if(std::string(argv[4])=="native") {
#if WV_TEST_NATIVE_FFTW
      require(WVFFTWEngine::create(1,fft));
#else
      throw std::runtime_error("Native provider not built.");
#endif
    } else throw std::runtime_error("Unknown FFT provider.");
    // Explicit scalar setup precedes the retained-storage baseline. Keep the
    // borrowed constant kernel alive until after the field service is destroyed.
    std::unique_ptr<WVTransformConstantStratificationKernel> constantKernel;
    std::unique_ptr<WVFieldEvaluationService> fields;
    if(checkpoint.stratifiedModalSource) require(WVFieldEvaluationService::create(checkpoint.stratifiedModalSource,std::move(fft),fields));
    else if(checkpoint.transformKind==WVPersistedTransformKind::barotropicQG)
      require(WVFieldEvaluationService::create(checkpoint.barotropicQGConfiguration,std::move(fft),fields));
    else {
      require(WVTransformConstantStratificationKernel::create(checkpoint.configuration,std::move(fft),constantKernel));
      require(constantKernel->prepareScalarAdvection());
      require(WVFieldEvaluationService::createBorrowing(*constantKernel,fields));
    }
    WVIntegrationStateLayout layout;require(WVIntegrationStateLayout::createCoefficientOnly(checkpoint.stateDescription,layout));
    std::vector<WVCoefficientFamilyConstView> families;
    std::vector<std::vector<WVComplex64>> before;
    for(std::size_t i=0;i<layout.coefficientFamilyCount();++i) {
      const auto& values=checkpoint.transformKind==WVPersistedTransformKind::constantStratification ?
          (i==0 ? checkpoint.state.coefficients.Ap : i==1 ? checkpoint.state.coefficients.Am : checkpoint.state.coefficients.A0) :
          checkpoint.transformState.coefficientFamilies[i].values;
      families.push_back({&layout.coefficientFamilies()[i],values.data()});
      before.push_back(values);
    }
    WVIntegrationState state;state.coefficientFamilies=families.data();state.coefficientFamilyCount=families.size();
    state.waveVortex=checkpoint.state.view();
    const auto names=request.at("fields").get<std::vector<std::string>>();
    std::vector<WVFieldRequest> requests;
    for(const auto& name:names) requests.push_back({name,name,{}});
    WVDensityDiagnosticContract contract;
    const auto selection=request.value("reference",std::string("actual"));
    if(selection=="initial") contract.reference=WVNoMotionReference::initial;
    else if(selection!="actual") throw std::runtime_error("Invalid density reference.");
    WVFieldEvaluationPlan plan;
    require(fields->createPlan(requests,plan,contract));
    std::vector<WVFieldEvaluationPlan> singlePlans(names.size());
    for(std::size_t i=0;i<names.size();++i)
      require(fields->createPlan({requests[i]},singlePlans[i],contract));
    std::vector<std::vector<double>> values(names.size());
    std::vector<WVFieldOutputView> views;
    for(std::size_t i=0;i<names.size();++i) {
      values[i].resize(plan.outputs()[i].elementCount);
      views.push_back({values[i].data(),values[i].size()});
    }
    const auto retainedBytes=[&] {
      auto bytes=fields->persistentBytes()+plan.persistentBytes()+(constantKernel ? constantKernel->persistentBytes() : 0);
      for(const auto& single:singlePlans) bytes+=single.persistentBytes();
      return bytes;
    };
    {
      detail::WVFieldEvaluationEventScope event(*fields,state,true,false);
      require(event.status());
      // Separate coincident consumers prepare each dependency at most once.
      for(std::size_t i=0;i<names.size();++i) {
        require(fields->evaluate(singlePlans[i],state,&views[i],1));
      }
      const auto prior=values;
      require(fields->evaluate(plan,state,views.data(),views.size()));
      if(prior!=values) throw std::runtime_error("Coincident density batching changed values.");
    }
    const auto metrics=fields->metrics();
    if(metrics.densityWorkspaceLiveBytes || metrics.eventFieldWorkspaceLiveBytes)
      throw std::runtime_error("Evaluation left event scratch live.");
    // The first event establishes bounded lazy prepared-field cache capacity.
    // Replays must reuse that capacity while event scratch remains transient.
    const auto retained=retainedBytes();
    json result;
    for(std::size_t i=0;i<names.size();++i) {
      result["fields"][names[i]]["dimensions"]=plan.outputs()[i].dimensions;
      result["fields"][names[i]]["real"]=values[i];
    }
    result["metrics"]={{"recoveryCount",metrics.densityRecoveryCount},
      {"profileConstructionCount",metrics.densityProfileConstructionCount},
      {"inversePassCount",metrics.densityInversePassCount},{"apePassCount",metrics.densityAPEPassCount},
      {"apvPassCount",metrics.densityAPVPassCount},{"apvReuseCount",metrics.densityAPVReuseCount},
      {"reuseCount",metrics.densityReuseCount},{"liveBytes",metrics.densityWorkspaceLiveBytes},
      {"highWaterBytes",metrics.densityWorkspaceHighWaterBytes},
      {"eventLiveBytes",metrics.eventFieldWorkspaceLiveBytes}};
    result["reference"]=selection;
    result["retainedBytes"]=retained;
    for(std::size_t i=0;i<before.size();++i)
      if(std::memcmp(before[i].data(),families[i].data,before[i].size()*sizeof(WVComplex64)))
        throw std::runtime_error("Density evaluation modified coefficient state.");
    // A later event recomputes and reproduces the same values.
    const auto prior=values;
    require(fields->evaluate(plan,state,views.data(),views.size()));
    if(prior!=values) throw std::runtime_error("Density replay changed values.");
    if(retained!=retainedBytes()) throw std::runtime_error("Density replay grew retained storage.");
    if(fields->metrics().densityWorkspaceLiveBytes || fields->metrics().eventFieldWorkspaceLiveBytes)
      throw std::runtime_error("Density replay left event scratch live.");
    result["replayRecoveryCount"]=fields->metrics().densityRecoveryCount;
    std::ofstream out(argv[3]);out<<result.dump()<<'\n';
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
