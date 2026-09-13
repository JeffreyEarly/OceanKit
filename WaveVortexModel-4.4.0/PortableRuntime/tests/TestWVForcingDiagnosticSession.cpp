#include "WVForcingDiagnosticWorkspace.hpp"

#include <array>
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::detail;

namespace {
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
struct Forcing {
  std::size_t index;
  WVForcingStage stage() const {return index<2 ? WVForcingStage::spatial : WVForcingStage::spectral;}
};
void exercise(WVVariableEvaluationPolicy policy,bool injectFailure) {
  std::vector<std::unique_ptr<Forcing>> forcing;
  std::vector<WVForcingStage> stages;
  for(std::size_t index=0;index<4;++index) {
    forcing.push_back(std::make_unique<Forcing>(Forcing{index}));
    stages.push_back(forcing.back()->stage());
  }
  WVVariableEvaluationContext evaluation;
  require(bool(evaluation.prepare(WVForcingDiagnosticWorkspace::dependencyKeys(4))),"prefix key preparation failed");
  int owner=0,foreign=0;
  require(bool(evaluation.begin(&owner,policy)),"prefix scope failed");
  WVForcingDiagnosticWorkspace work({1,1},{1,1,1,1},1,0);
  require(bool(work.beginScopedEvaluation(evaluation,stages)),"prefix storage preparation failed");
  WVComplex64 coefficient{1,0};
  const WVState state{0,0,{{},{},{&coefficient,{1,1}}}};
  require(bool(work.bind(&owner,state)),"prefix binding failed");
  require(!work.bind(&foreign,state),"foreign prefix owner accepted");
  work.markInitialized();
  WVForcingTendencyMetrics metrics;
  std::array<std::size_t,4> successes{};
  std::size_t projections=0,reconstructions=0;
  bool shouldFail=injectFailure;
  const auto execute=[&](const Forcing& item,WVFlux& flux) {
    if(item.index<2) {
      const double value=static_cast<double>(item.index+1);
      const auto status=work.addSpatial({&value,{1,1,1,1}});
      if(!status) return status;
    } else if(item.index==2) flux.F0.data[0].real+=5;
    else {
      flux.F0.data[0].real*=0.5;
      if(shouldFail) {shouldFail=false;throw std::runtime_error("injected producer failure");}
    }
    ++successes[item.index];
    return WVKernelStatus::ok();
  };
  const auto project=[&](WVRealFieldBundleConstView value,WVFlux& flux) {
    ++projections;flux.F0.data[0]={10*value.data[0],0};return WVKernelStatus::ok();
  };
  const auto reconstruct=[&](std::vector<WVComplex64>& value,WVRealFieldBundleView output) {
    ++reconstructions;output.data[0]=0.1*value[0].real;return WVKernelStatus::ok();
  };
  double output=-99;
  WVForcingTendencyOutput request{0,{&output,{1,1,1,1}}};
  auto run=[&] {return evaluateForcingTendencySequence(forcing,work,&request,1,metrics,execute,project,reconstruct);};
  require(bool(run()) && output==1 && successes==std::array<std::size_t,4>{1,0,0,0},"first query executed an unneeded suffix");
  request.executionIndex=3;output=-99;
  if(injectFailure) {
    bool caught=false;
    try {run();} catch(const std::runtime_error&) {caught=true;}
    require(caught && output==-99,"failed batch published output");
    if(policy==WVVariableEvaluationPolicy::reuse)
      require(work.nextIndex==3 && work.flux[0].real==35,"failed suffix did not preserve the completed prefix");
  }
  require(bool(run()) && output==-1.75,"ordered suffix result changed");
  request.executionIndex=2;
  require(bool(run()) && output==0.5,"late request lost an earlier spectral difference");
  const auto before=successes;
  require(bool(run()) && output==0.5,"repeated forcing result changed");
  if(policy==WVVariableEvaluationPolicy::reuse) {
    require(successes==std::array<std::size_t,4>{1,1,1,1} && before==successes,"reuse repeated an already completed forcing");
    require(projections==1 && reconstructions==2,"reuse repeated projection or physicalization");
    require(evaluation.metrics().recomputations==0 && evaluation.metrics().duplicateExecutions==0,"reuse ledger reports duplicate work");
  } else {
    require(successes[0]>1 && evaluation.metrics().recomputations>0 && evaluation.metrics().evictions>0,"low-memory replay was not explicit");
    require(work.prefix.empty(),"low-memory retained per-forcing differences");
  }
  evaluation.end();
  require(!work.bind(&owner,state),"expired prefix accepted");
  const auto retained=work.bytes();
  work.endScopedEvaluation();
  require(work.bytes()==retained && work.flux.empty() && work.physical.empty(),
      "completed scope lost capacity or kept live fields");
  require(bool(evaluation.begin(&owner,policy)),"second context failed");
  require(bool(work.beginScopedEvaluation(evaluation,stages)) && work.bytes()==retained,
      "prepared forcing storage grew on the next evaluation");
  require(bool(work.bind(&owner,state)) && !work.initialized() && work.nextIndex==0,
      "next scope retained the old prefix identity");
  work.endScopedEvaluation();evaluation.end();
}
void sharedProducers(WVVariableEvaluationPolicy policy) {
  WVVariableEvaluationContext evaluation;
  int owner=0;
  require(bool(evaluation.prepare(WVForcingDiagnosticWorkspace::dependencyKeys(2))),"shared keys failed");
  require(bool(evaluation.begin(&owner,policy)),"shared context failed");
  WVForcingDiagnosticWorkspace work({1,1},{1,1,1,1},1,0);
  work.nonlinearUseCount=2;
  require(bool(work.beginScopedEvaluation(evaluation,{WVForcingStage::spatial,WVForcingStage::spatial})),"shared workspace failed");
  std::size_t reductions=0,nonlinear=0;
  const auto reduce=[&](double& value) {++reductions;value=7;return WVKernelStatus::ok();};
  double maximum=0;
  require(bool(work.evaluateHorizontalMaximum(maximum,reduce)) && maximum==7,"first reduction failed");
  require(bool(work.evaluateHorizontalMaximum(maximum,reduce)) && reductions==1,"shared reduction repeated");
  const auto produce=[&] {++nonlinear;work.raw[0]=3;return WVKernelStatus::ok();};
  require(bool(work.evaluateNonlinearRaw(produce)),"first nonlinear producer failed");
  work.raw[0]=-99;
  require(bool(work.evaluateNonlinearRaw(produce)) && work.raw[0]==3,"shared nonlinear result changed");
  require(nonlinear==(policy==WVVariableEvaluationPolicy::reuse ? 1u : 2u),"wrong nonlinear production count");
  require(evaluation.metrics().recomputations==(policy==WVVariableEvaluationPolicy::reuse ? 0u : 1u),"nonlinear replay not counted");
  work.endScopedEvaluation();evaluation.end();
}

}
int main() {
  try {
    for(auto policy:{WVVariableEvaluationPolicy::reuse,WVVariableEvaluationPolicy::lowMemory}) {
      for(bool failure:{false,true}) exercise(policy,failure);
      sharedProducers(policy);
    }
  } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
  std::cout<<"Forcing diagnostic prefix lifecycle passed\n";
}
