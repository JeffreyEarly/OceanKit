#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {
template<class T> void require(T status) { if (!status) throw std::runtime_error(status.message); }
json values(const std::vector<WVComplex64>& a) {
  std::vector<double> real,imag;
  for (const auto x:a) { real.push_back(x.real); imag.push_back(x.imag); }
  return {{"real",real},{"imag",imag}};
}
template<class Engine,class State>
void tendencies(Engine& engine,const State& state,WVShape4D shape,json& result) {
  std::vector<std::vector<double>> storage(engine.forcingCount());
  std::vector<WVForcingTendencyOutput> outputs;
  for (std::size_t index=0;index<storage.size();++index) {
    storage[index].resize(shape.elementCount());
    outputs.push_back({index,{storage[index].data(),shape}});
  }
  const auto bytes=engine.persistentBytes();
  require(engine.evaluateForcingTendencies(state,outputs.data(),outputs.size()));
  if (engine.persistentBytes()!=bytes) throw std::runtime_error("Diagnostic evaluation retained workspace.");
  result["tendencies"]=json::array();
  const auto R=shape.first*shape.second*shape.third;
  const char* hydrostatic[]={"Fu","Fv","Feta"};
  const char* nonhydrostatic[]={"Fu","Fv","Fw","Feta"};
  const char* qg[]={"Fqgpv"};
  const auto* names=shape.fourth==1 ? qg : (shape.fourth==3 ? hydrostatic : nonhydrostatic);
  for (std::size_t index=0;index<storage.size();++index) {
    const auto* forcing=engine.forcingInstance(index);
    json entry={{"name",forcing->name()},{"type",forcing->typeIdentifier()},
      {"ordinal",forcing->ordinal()},{"stage",forcingStageName(forcing->stage())}};
    for (std::size_t channel=0;channel<shape.fourth;++channel)
      entry["fields"][names[channel]]=std::vector<double>(storage[index].begin()+channel*R,storage[index].begin()+(channel+1)*R);
    result["tendencies"].push_back(std::move(entry));
  }
  result["diagnosticWorkspaceHighWaterBytes"]=engine.tendencyMetrics().workspaceHighWaterBytes;
  result["diagnosticWorkspaceLiveBytes"]=engine.tendencyMetrics().workspaceLiveBytes;
  result["diagnosticForcingEvaluationCount"]=engine.tendencyMetrics().forcingEvaluationCount;
}
}
int main(int argc,char** argv) {
  try {
    if (argc != 4 && argc != 5) throw std::runtime_error("Usage: WVStableForcingDump model.nc output.json reference|native [tendencies]");
    const bool diagnostics=argc==5;
    if (diagnostics && std::string(argv[4])!="tendencies") throw std::runtime_error("Unknown diagnostic request.");
    WVExtensionCatalogBuilder builder; require(addBuiltInExtensions(builder));
    std::shared_ptr<const WVExtensionCatalog> catalog; require(builder.freeze(catalog));
    WVCheckpoint checkpoint; require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
    std::unique_ptr<WVFFTEngine> fft;
    const std::string provider(argv[3]);
    if (provider == "reference") fft=std::make_unique<WVReferenceFFTEngine>();
    else if (provider == "native") {
#if WV_TEST_NATIVE_FFTW
      require(WVFFTWEngine::create(1,fft));
#else
      throw std::runtime_error("Native adapter unavailable.");
#endif
    } else throw std::runtime_error("Unknown provider.");
    json result;
    if (checkpoint.transformKind == WVPersistedTransformKind::barotropicQG) {
      std::unique_ptr<WVBarotropicQGForcingEngine> engine;
      require(WVBarotropicQGForcingEngine::create(checkpoint.barotropicQGConfiguration,checkpoint.forcingSchedule,catalog,std::move(fft),engine));
      auto& a=checkpoint.transformState.coefficientFamilies.at(0).values;
      std::vector<WVComplex64> f(a.size());
      WVComplexConstView input{a.data(),engine->kernel().descriptor().spectralShape()};
      WVComplexView output{f.data(),input.shape};
      require(engine->evaluateRightHandSide(input,output));
      allocationProbe::calls=0; allocationProbe::counting=true;
      const auto status=engine->evaluateRightHandSide(input,output);
      allocationProbe::counting=false; require(status);
      result["F0"]=values(f);
      if (diagnostics) {
        const auto plane=engine->kernel().descriptor().spatialShape();
        tendencies(*engine,input,{plane.rows,plane.columns,1,1},result);
      }
    } else if (checkpoint.transformKind==WVPersistedTransformKind::stratifiedQG) {
      std::unique_ptr<WVStratifiedQGForcingEngine> engine;
      require(WVStratifiedQGForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),engine));
      auto& a=checkpoint.transformState.coefficientFamilies.at(0).values;
      std::vector<WVComplex64> f(a.size());
      WVComplexConstView input{a.data(),engine->kernel().spectralShape()}; WVComplexView output{f.data(),input.shape};
      require(engine->evaluateRightHandSide(input,output));
      allocationProbe::calls=0; allocationProbe::counting=true;
      const auto status=engine->evaluateRightHandSide(input,output);
      allocationProbe::counting=false; require(status);
      result["F0"]=values(f);
      if (diagnostics) {
        const auto volume=engine->kernel().spatialShape();
        tendencies(*engine,input,{volume.first,volume.second,volume.third,1},result);
      }
    } else if (checkpoint.transformKind==WVPersistedTransformKind::hydrostatic) {
      std::unique_ptr<WVHydrostaticForcingEngine> engine;
      require(WVHydrostaticForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),engine));
      const auto shape=engine->kernel().spectralShape();
      auto& families=checkpoint.transformState.coefficientFamilies;
      WVState state{checkpoint.state.t,checkpoint.state.t0,{{families[0].values.data(),shape},{families[1].values.data(),shape},{families[2].values.data(),shape}}};
      std::vector<WVComplex64> fp(shape.elementCount()),fm(fp.size()),f0(fp.size());
      WVFlux flux{{fp.data(),shape},{fm.data(),shape},{f0.data(),shape}};
      require(engine->nonlinearFlux(state,flux));
      allocationProbe::calls=0; allocationProbe::counting=true;
      const auto status=engine->nonlinearFlux(state,flux);
      allocationProbe::counting=false; require(status);
      result["Fp"]=values(fp); result["Fm"]=values(fm); result["F0"]=values(f0);
      if (diagnostics) {
        const auto shape=engine->kernel().spatialShape();
        tendencies(*engine,state,{shape.first,shape.second,shape.third,3},result);
      }
    } else if (checkpoint.transformKind==WVPersistedTransformKind::boussinesq) {
      std::unique_ptr<WVBoussinesqForcingEngine> engine;
      require(WVBoussinesqForcingEngine::create(checkpoint.stratifiedModalSource,checkpoint.forcingSchedule,catalog,std::move(fft),engine));
      const auto shape=engine->kernel().spectralShape();
      auto& families=checkpoint.transformState.coefficientFamilies;
      WVState state{checkpoint.state.t,checkpoint.state.t0,{{families[0].values.data(),shape},{families[1].values.data(),shape},{families[2].values.data(),shape}}};
      std::vector<WVComplex64> fp(shape.elementCount()),fm(fp.size()),f0(fp.size());
      WVFlux flux{{fp.data(),shape},{fm.data(),shape},{f0.data(),shape}};
      require(engine->nonlinearFlux(state,flux));
      allocationProbe::calls=0; allocationProbe::counting=true;
      const auto status=engine->nonlinearFlux(state,flux);
      allocationProbe::counting=false; require(status);
      result["Fp"]=values(fp); result["Fm"]=values(fm); result["F0"]=values(f0);
      if (diagnostics) {
        const auto shape=engine->kernel().spatialShape();
        tendencies(*engine,state,{shape.first,shape.second,shape.third,4},result);
      }
    } else {
      std::unique_ptr<WVConstantStratificationForcingEngine> engine;
      require(WVConstantStratificationForcingEngine::create(checkpoint.configuration,checkpoint.forcingSchedule,catalog,std::move(fft),engine));
      const auto shape=checkpoint.state.coefficients.shape;
      std::vector<WVComplex64> fp(shape.elementCount()),fm(fp.size()),f0(fp.size());
      WVFlux flux{{fp.data(),shape},{fm.data(),shape},{f0.data(),shape}};
      require(engine->nonlinearFlux(checkpoint.state.view(),flux));
      allocationProbe::calls=0; allocationProbe::counting=true;
      const auto status=engine->nonlinearFlux(checkpoint.state.view(),flux);
      allocationProbe::counting=false; require(status);
      result["Fp"]=values(fp); result["Fm"]=values(fm); result["F0"]=values(f0);
      if (diagnostics) {
        const auto& c=checkpoint.configuration;
        tendencies(*engine,checkpoint.state.view(),{c.Nx,c.Ny,c.Nz,c.isHydrostatic ? 3U : 4U},result);
      }
    }
    result["applicationAllocations"]=allocationProbe::calls.load();
    std::ofstream out(argv[2]); out<<result.dump()<<'\n';
    if (!out) throw std::runtime_error("Output write failed.");
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
