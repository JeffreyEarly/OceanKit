#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <fstream>
#include <iostream>
#include <set>
using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {template<class T> void require(T status) {if(!status) throw std::runtime_error(status.message);}}
int main(int argc,char** argv) {
  try {
    if(argc!=5) throw std::runtime_error("Usage: WVStratifiedQGFieldDump model.nc positions.json output.json reference|native");
    WVExtensionCatalogBuilder builder;require(addBuiltInExtensions(builder));std::shared_ptr<const WVExtensionCatalog> catalog;require(builder.freeze(catalog));
    WVCheckpoint checkpoint;require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
    std::ifstream requestFile(argv[2]);json request;requestFile>>request;
    auto x=request.at("x").get<std::vector<double>>(),y=request.at("y").get<std::vector<double>>(),z=request.at("z").get<std::vector<double>>();
    const auto names=request.at("fields").get<std::vector<std::string>>();
    std::unique_ptr<WVFFTEngine> fft;
    if(std::string(argv[4])=="reference") fft=std::make_unique<WVReferenceFFTEngine>();
    else if(std::string(argv[4])=="native") {
#if WV_TEST_NATIVE_FFTW
      require(WVFFTWEngine::create(1,fft));
#else
      throw std::runtime_error("Native provider not built.");
#endif
    } else throw std::runtime_error("Unknown FFT provider.");
    std::unique_ptr<WVConstantStratificationForcingEngine> constantEngine;
    std::unique_ptr<WVBarotropicQGForcingEngine> barotropicEngine;
    std::unique_ptr<WVStratifiedQGForcingEngine> sqgEngine;
    std::unique_ptr<WVHydrostaticForcingEngine> hydrostaticEngine;
    std::unique_ptr<WVBoussinesqForcingEngine> boussinesqEngine;
    std::unique_ptr<WVFieldEvaluationService> fields;
    const bool legacyStratifiedProbe=static_cast<bool>(checkpoint.stratifiedModalSource);
    if(checkpoint.transformKind==WVPersistedTransformKind::stratifiedQG) {
      require(WVStratifiedQGForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),sqgEngine));
      require(WVFieldEvaluationService::createBorrowing(*sqgEngine,fields));
    } else if(checkpoint.transformKind==WVPersistedTransformKind::hydrostatic) {
      require(WVHydrostaticForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),hydrostaticEngine));
      require(WVFieldEvaluationService::createBorrowing(*hydrostaticEngine,fields));
    } else if(checkpoint.transformKind==WVPersistedTransformKind::boussinesq) {
      require(WVBoussinesqForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),boussinesqEngine));
      require(WVFieldEvaluationService::createBorrowing(*boussinesqEngine,fields));
    } else if(checkpoint.transformKind==WVPersistedTransformKind::barotropicQG) {
      require(WVBarotropicQGForcingEngine::create(checkpoint.barotropicQGConfiguration,checkpoint.forcingSchedule,catalog,std::move(fft),barotropicEngine));
      require(WVFieldEvaluationService::createBorrowing(*barotropicEngine,fields));
    } else {
      require(WVConstantStratificationForcingEngine::create(checkpoint.configuration,checkpoint.forcingSchedule,catalog,std::move(fft),constantEngine));
      require(constantEngine->kernel().prepareScalarAdvection());
      require(WVFieldEvaluationService::createBorrowing(*constantEngine,fields));
    }
    WVIntegrationStateLayout layout;require(WVIntegrationStateLayout::createCoefficientOnly(checkpoint.stateDescription,layout));
    std::vector<WVCoefficientFamilyConstView> families;
    for (std::size_t i=0;i<layout.coefficientFamilyCount();++i) {
      const auto& values=checkpoint.transformKind==WVPersistedTransformKind::constantStratification ?
          (i==0 ? checkpoint.state.coefficients.Ap : i==1 ? checkpoint.state.coefficients.Am : checkpoint.state.coefficients.A0) :
          checkpoint.transformState.coefficientFamilies[i].values;
      families.push_back({&layout.coefficientFamilies()[i],values.data()});
    }
    WVIntegrationState state;state.coefficientFamilies=families.data();state.coefficientFamilyCount=families.size();
    state.waveVortex=checkpoint.state.view();
    json result;
    const auto profileRequest=request.value("profiles",json::object());
    const auto profileNames=profileRequest.value("fields",std::vector<std::string>{});
    const auto profileX=profileRequest.value("xIndices",std::vector<std::size_t>{});
    const auto profileY=profileRequest.value("yIndices",std::vector<std::size_t>{});
    if(!profileNames.empty() && profileX.size()!=profileNames.size()) throw std::runtime_error("Profile x-index count differs from profile fields.");
    if(!profileNames.empty() && profileY.size()!=profileNames.size()) throw std::runtime_error("Profile y-index count differs from profile fields.");
    WVDensityDiagnosticContract densityContract;
    const auto reference=request.value("reference",std::string("actual"));
    if(reference=="initial") densityContract.reference=WVNoMotionReference::initial;
    else if(reference!="actual") throw std::runtime_error("Unknown density reference.");
    // Supply the complete numerical fields independently of the sampling
    // routes, so MATLAB can test interpolation separately from reconstruction.
    std::set<std::string> fullNames(names.begin(),names.end());
    fullNames.insert(profileNames.begin(),profileNames.end());
    std::vector<WVFieldRequest> fullRequests;
    for(const auto& name:fullNames) fullRequests.push_back({name,name,{}});
    WVFieldEvaluationPlan fullPlan;require(fields->createPlan(fullRequests,fullPlan,densityContract));
    std::vector<std::vector<double>> fullStorage(fullRequests.size());
    std::vector<WVFieldOutputView> fullOutputs;
    for(std::size_t i=0;i<fullRequests.size();++i) {
      fullStorage[i].resize(fullPlan.outputs()[i].elementCount);
      fullOutputs.push_back({fullStorage[i].data(),fullStorage[i].size()});
    }
    require(fields->evaluate(fullPlan,state,fullOutputs.data(),fullOutputs.size()));
    for(std::size_t i=0;i<fullRequests.size();++i) result["full"][fullRequests[i].fieldName]=fullStorage[i];
    for(auto interpolation:{WVPositionInterpolation::linear,WVPositionInterpolation::spline}) {
      const std::string method=interpolation==WVPositionInterpolation::linear?"linear":"spline";
      std::vector<WVFieldRequest> requests;
      std::vector<WVMovingFieldRequest> movingRequests;
      std::vector<std::string> movingNames;
      std::vector<WVEventFieldRequest> eventRequests;
      for(const auto& name:names) {
        WVFieldSamplingRequest sampling;sampling.kind=WVFieldSamplingKind::positions;sampling.interpolation=interpolation;sampling.x=x;sampling.y=y;sampling.z=z;
        requests.push_back({name,name,sampling});
        movingNames.push_back(name);
        movingRequests.push_back({name,name,0,x.size(),interpolation});
        eventRequests.push_back({name,name,0,interpolation});
      }
      WVFieldEvaluationPlan plan;require(fields->createPlan(requests,plan,densityContract));
      std::vector<std::vector<double>> storage(names.size(),std::vector<double>(x.size()));std::vector<WVFieldOutputView> output;
      for(auto& values:storage) output.push_back({values.data(),values.size()});
      require(fields->evaluate(plan,state,output.data(),output.size()));
      for(std::size_t i=0;i<names.size();++i) result[method]["fixed"][names[i]]=storage[i];
      if(!profileNames.empty()) {
        std::vector<WVFieldRequest> profileRequests;
        for(std::size_t i=0;i<profileNames.size();++i) {
          WVFieldSamplingRequest sampling; sampling.kind=WVFieldSamplingKind::fixedVerticalProfiles;
          sampling.xIndices={profileX[i]}; sampling.yIndices={profileY[i]};
          profileRequests.push_back({profileNames[i],profileNames[i],sampling});
        }
        WVFieldEvaluationPlan profilePlan;require(fields->createPlan(profileRequests,profilePlan,densityContract));
        std::vector<std::vector<double>> profileStorage(profileNames.size());std::vector<WVFieldOutputView> profileOutput;
        for(std::size_t i=0;i<profileNames.size();++i) {profileStorage[i].resize(profilePlan.outputs()[i].elementCount);profileOutput.push_back({profileStorage[i].data(),profileStorage[i].size()});}
        require(fields->evaluate(profilePlan,state,profileOutput.data(),profileOutput.size()));
        for(std::size_t i=0;i<profileNames.size();++i) result[method]["profiles"][profileNames[i]]=profileStorage[i];
      }
      WVMovingFieldEvaluationPlan moving;require(fields->createMovingPlan(movingRequests,moving,densityContract));
      std::vector<std::vector<double>> movingStorage(movingNames.size(),std::vector<double>(x.size()));
      std::vector<WVFieldOutputView> movingOutput;
      for(auto& values:movingStorage) movingOutput.push_back({values.data(),values.size()});
      allocationProbe::calls=0;allocationProbe::counting=true;
      auto status=fields->evaluateMoving(moving,state,{x.data(),y.data(),z.data(),x.size()},movingOutput.data(),movingOutput.size());
      allocationProbe::counting=false;require(status);
      result[method]["allocations"]=allocationProbe::calls.load();
      for(std::size_t i=0;i<movingNames.size();++i) result[method]["moving"][movingNames[i]]=movingStorage[i];
      WVEventFieldEvaluationPlan event;require(fields->createEventPlan(eventRequests,event,densityContract));
      WVEventPositionSetView positions{x.data(),y.data(),z.data(),x.size(),nullptr,0};WVPreparedFieldGeometry geometry;
      require(fields->prepareEventGeometry(event,&positions,1,geometry));
      WVEventFieldEvaluationBatchEntry occurrence{&event,&geometry,output.data(),output.size()};
      require(fields->evaluateEventBatch(state,&occurrence,1));
      for(std::size_t i=0;i<names.size();++i) result[method]["event"][names[i]]=storage[i];
      if(legacyStratifiedProbe) {
        WVEventFieldEvaluationPlan otherEvent;require(fields->createEventPlan(eventRequests,otherEvent,densityContract));
        occurrence.plan=&otherEvent;
        if(fields->evaluateEventBatch(state,&occurrence,1)) throw std::runtime_error("Mismatched event geometry was accepted.");
        auto invalidPositions=positions;invalidPositions.extentCount=1;invalidPositions.extents=nullptr;
        if(fields->prepareEventGeometry(event,&invalidPositions,1,geometry)) throw std::runtime_error("Null event extents were accepted.");
      }
    }
    std::ofstream out(argv[3]);out<<result.dump()<<'\n';
  } catch(const std::exception& e) {allocationProbe::counting=false;std::cerr<<e.what()<<'\n';return 1;}
}
