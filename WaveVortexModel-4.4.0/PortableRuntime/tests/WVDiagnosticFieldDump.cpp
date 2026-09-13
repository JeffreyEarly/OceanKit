#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
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
    if(argc!=5) throw std::runtime_error("Usage: WVDiagnosticFieldDump model.nc requests.json output.json reference|native");
    WVExtensionCatalogBuilder builder;require(addBuiltInExtensions(builder));
    std::shared_ptr<const WVExtensionCatalog> catalog;require(builder.freeze(catalog));
    WVCheckpoint checkpoint;require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
    std::ifstream requestFile(argv[2]);json request;requestFile>>request;
    std::unique_ptr<WVFFTEngine> fft;
    if(std::string(argv[4])=="reference") fft=std::make_unique<WVReferenceFFTEngine>();
    else {
#if WV_TEST_NATIVE_FFTW
      require(WVFFTWEngine::create(1,fft));
#else
      throw std::runtime_error("Native provider not built.");
#endif
    }
    std::unique_ptr<WVFieldEvaluationService> fields;
    if(checkpoint.stratifiedModalSource) require(WVFieldEvaluationService::create(checkpoint.stratifiedModalSource,std::move(fft),fields));
    else if(checkpoint.transformKind==WVPersistedTransformKind::barotropicQG)
      require(WVFieldEvaluationService::create(checkpoint.barotropicQGConfiguration,std::move(fft),fields));
    else require(WVFieldEvaluationService::create(checkpoint.configuration,std::move(fft),fields));
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
    WVFieldEvaluationPlan plan;require(fields->createPlan(requests,plan));
    std::vector<std::vector<double>> real(names.size());
    std::vector<std::vector<WVComplex64>> complex(names.size());
    std::vector<WVFieldOutputView> views;
    for(std::size_t i=0;i<names.size();++i) {
      const auto& spec=plan.outputs()[i];
      if(spec.isComplex) {complex[i].resize(spec.elementCount);views.push_back({nullptr,spec.elementCount,complex[i].data()});}
      else {real[i].resize(spec.elementCount);views.push_back({real[i].data(),spec.elementCount});}
    }
    require(fields->evaluate(plan,state,views.data(),views.size()));
    const auto metrics=fields->metrics();
    if(metrics.diagnosticWorkspaceLiveBytes) throw std::runtime_error("Evaluation left diagnostic scratch live.");
    json result;
    for(std::size_t i=0;i<names.size();++i) {
      auto& value=result["fields"][names[i]];
      value["dimensions"]=plan.outputs()[i].dimensions;
      if(plan.outputs()[i].isComplex) {
        std::vector<double> re,im;
        for(auto z:complex[i]) {re.push_back(z.real);im.push_back(z.imag);}
        value["real"]=re;value["imag"]=im;
      } else value["real"]=real[i];
    }
    result["diagnosticEvaluations"]=metrics.diagnosticEvaluationCount;
    result["primitiveOutputs"]=metrics.diagnosticPrimitiveOutputCount;
    result["scratchHighWaterBytes"]=metrics.diagnosticWorkspaceHighWaterBytes;
    result["scratchLiveBytes"]=metrics.diagnosticWorkspaceLiveBytes;
    for(std::size_t i=0;i<before.size();++i)
      if(std::memcmp(before[i].data(),families[i].data,before[i].size()*sizeof(WVComplex64)))
        throw std::runtime_error("Diagnostic evaluation modified coefficient state.");
    const auto evaluateSingles=[&] {
      for(std::size_t i=0;i<names.size();++i) {
        WVFieldEvaluationPlan single;require(fields->createPlan({requests[i]},single));
        auto re=real[i];auto im=complex[i];
        WVFieldOutputView view{re.data(),views[i].elementCount,im.data()};
        require(fields->evaluate(single,state,&view,1));
        if(fields->metrics().diagnosticWorkspaceLiveBytes)
          throw std::runtime_error("Independent diagnostic request left scratch live.");
        if(re!=real[i] || (!im.empty() && std::memcmp(im.data(),complex[i].data(),im.size()*sizeof(WVComplex64))))
          throw std::runtime_error("Diagnostic batching changed values: "+names[i]);
      }
    };
    // Combined and independently planned requests can establish different
    // bounded prepared-field peaks. Warm both before fixing retained capacity.
    evaluateSingles();
    const auto retained=fields->persistentBytes()+plan.persistentBytes();
    result["retainedBytes"]=retained;
    require(fields->evaluate(plan,state,views.data(),views.size()));
    if(retained!=fields->persistentBytes()+plan.persistentBytes())
      throw std::runtime_error("Diagnostic replay grew retained storage.");
    if(fields->metrics().diagnosticWorkspaceLiveBytes)
      throw std::runtime_error("Diagnostic replay left scratch live.");
    evaluateSingles();
    if(retained!=fields->persistentBytes()+plan.persistentBytes())
      throw std::runtime_error("Independent diagnostic replay grew retained storage.");
    auto invalidFamilies=families;
    auto malformed=layout.coefficientFamilies()[0];malformed.elementCount=1;
    invalidFamilies[0].layout=&malformed;
    auto invalidState=state;invalidState.coefficientFamilies=invalidFamilies.data();
    if(fields->evaluate(plan,invalidState,views.data(),views.size())) throw std::runtime_error("Malformed coefficient layout accepted.");
    for(const auto& name:request.value("rejected",std::vector<std::string>{})) {
      WVFieldEvaluationPlan rejected;const auto status=fields->createPlan({{name,name,{}}},rejected);
      if(status) throw std::runtime_error("Unsupported diagnostic accepted: "+name);
      result["rejected"][name]=status.message;
    }
    std::ofstream out(argv[3]);out<<result.dump()<<'\n';
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
