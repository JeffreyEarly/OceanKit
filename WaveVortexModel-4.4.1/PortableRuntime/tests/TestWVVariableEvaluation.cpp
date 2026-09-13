#include "WaveVortexRuntime/WVVariableEvaluation.hpp"
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
void lifecycle(WVVariableEvaluationPolicy policy) {
  WVVariableEvaluationContext context;
  const WVVariableEvaluationKey u{WVVariableEvaluationNode::physicalField,0};
  const WVVariableEvaluationKey v{WVVariableEvaluationNode::physicalField,1};
  const WVVariableEvaluationKey speed{WVVariableEvaluationNode::reduction,0};
  require(bool(context.prepare({u,v,speed})),"prepare failed");
  int state=3,owner=0,otherOwner=0,executions=0,computed=0;
  auto producer=[&] {++executions; computed=state; return WVKernelStatus::ok();};
  require(!context.evaluate(u,sizeof(double),producer),"evaluation outside scope succeeded");
  require(bool(context.begin(&owner,policy)),"begin failed");
  auto generation=context.generation();
  require(!context.begin(&owner,policy),"nested begin succeeded");
  require(!context.prepare({u}),"active plan mutation succeeded");
  require(!context.validate(&otherOwner,generation),"foreign owner accepted");
  require(bool(context.evaluate(u,sizeof(double),producer)),"first producer failed");
  require(bool(context.evaluate(u,sizeof(double),producer)),"cache hit failed");
  require(executions==1 && computed==3,"duplicate production");
  require(!context.evaluate(u,2*sizeof(double),producer),"cached extent mismatch accepted");
  require(executions==1,"extent mismatch invoked producer");
  require(context.pin(u),"pin failed");
  require(!context.evict(u),"pinned field evicted");
  context.unpin(u);
  const bool evicted=context.evict(u);
  require(evicted==(policy==WVVariableEvaluationPolicy::lowMemory),"policy eviction mismatch");
  require(bool(context.evaluate(u,sizeof(double),producer)),"post-eviction access failed");
  require(executions==(evicted ? 2 : 1),"unreported recomputation");
  const auto failure=WVKernelStatus{WVKernelStatusCode::invalidConfiguration,"injected"};
  require(!context.evaluate(v,8,[&]{return failure;}),"failed producer accepted");
  require(!context.ready(v),"failed value published");
  require(!context.evaluate(v,8,[&]{return context.evaluate(v,8,producer);}),"cycle accepted");
  try {context.evaluate(v,8,[]()->WVKernelStatus{throw std::runtime_error("injected");});}
  catch(const std::runtime_error&) {}
  require(bool(context.evaluate(v,8,producer)),"failed node did not recover");
  require(bool(context.evaluate(speed,8,[&]{return context.evaluate(u,8,producer);})),"dependency evaluation failed");
  require(context.metrics().duplicateExecutions==0,"default duplicate executions reported");
  require(context.metrics().recomputations==(evicted ? 1u : 0u),"recomputation accounting wrong");
  context.end();
  require(context.metrics().liveBytes==0 && !context.ready(u),"ended scope retained validity");
  require(!context.validate(&owner,generation),"ended view accepted");
  state=7; // Same owner and storage, changed coefficients between scopes.
  require(bool(context.begin(&owner,policy)),"new state begin failed");
  require(!context.validate(&owner,generation),"old generation accepted");
  require(bool(context.evaluate(u,8,producer)) && computed==7,"new state reused stale value");
  context.end();
}
void fusedLifecycle() {
  WVVariableEvaluationContext context;
  const WVVariableEvaluationKey a{WVVariableEvaluationNode::forcingTendency,0};
  const WVVariableEvaluationKey b{WVVariableEvaluationNode::forcingTendency,1};
  const WVVariableEvaluationKey c{WVVariableEvaluationNode::forcingTendency,2};
  const std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> group{{a,8},{b,16}};
  int owner=0,calls=0;
  require(bool(context.prepare({a,b,c})),"fused prepare failed");
  require(bool(context.begin(&owner,WVVariableEvaluationPolicy::lowMemory)),"fused begin failed");
  auto produce=[&]{++calls;return WVKernelStatus::ok();};
  require(!context.evaluateGroup({{a,8},{a,8}},produce),"duplicate group key accepted");
  require(!context.evaluateGroup(group,[&]{return context.evaluate(a,8,produce);}),"group dependency cycle accepted");
  require(!context.ready(a) && !context.ready(b),"failed group published partial result");
  try {context.evaluateGroup(group,[]()->WVKernelStatus{throw std::runtime_error("injected");});}
  catch(const std::runtime_error&) {}
  require(!context.ready(a) && !context.ready(b),"throwing group published result");
  require(bool(context.evaluateGroup(group,[&]{
    auto status=context.evaluate(c,4,produce); if(!status) return status;
    return produce();
  })),"group with independent dependency failed");
  require(calls==2 && context.metrics().producerExecutions==2,"fused siblings counted as separate producers");
  require(bool(context.evaluateGroup(group,produce)) && calls==2,"fused cache hit recomputed");
  require(!context.evaluateGroup({{a,16},{b,16}},produce),"fused extent mismatch accepted");
  require(context.pin(a) && !context.evict(a),"pinned group output evicted");
  context.unpin(a);
  require(context.evict(a),"group member eviction failed");
  require(!context.evaluateGroup(group,produce) && calls==2,"partial ready group recomputed silently");
  require(context.evict(b),"remaining group eviction failed");
  require(bool(context.evaluateGroup(group,produce)),"explicit group recomputation failed");
  require(calls==3 && context.metrics().recomputations==1 && context.metrics().evictions==2,"fused recomputation metrics wrong");
  require(context.metrics().liveBytes==28,"fused live bytes wrong");
  context.end();
}

}
int main() {
  try {lifecycle(WVVariableEvaluationPolicy::reuse); lifecycle(WVVariableEvaluationPolicy::lowMemory); fusedLifecycle();}
  catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return 1;}
  std::cout<<"Variable evaluation lifecycle passed\n";
}
