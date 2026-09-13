#include "WVForcingDiagnosticBinding.hpp"
#include "WVFieldEvaluationEventWorkspace.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WVDiagnosticFieldPlan.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WVBarotropicQGFieldEvaluationAdapter.hpp"
#include "WVStratifiedFieldEvaluationAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime {
namespace detail {
WVFieldEvaluationEventScope::WVFieldEvaluationEventScope(WVFieldEvaluationService& service,const WVIntegrationState& state,bool enabled,bool retainPrimitiveFields)
    : workspace_(state,service.variableEvaluationContext_,*service.eventArena_,
          service.variableEvaluationPolicy_,enabled) {
  if(!enabled) return;
  if(!service.variableEvaluationPreparationStatus_) {
    workspace_.endEvaluation();
    status_=service.variableEvaluationPreparationStatus_;
    return;
  }
  if(!workspace_.status()) {status_=workspace_.status(); return;}
  if(service.eventWorkspace_) {
    status_={WVKernelStatusCode::reentrantExecution,"An output field event is already active."};
    return;
  }
  service_=&service;
  workspace_.retainPrimitiveFields_=retainPrimitiveFields;
  service.eventWorkspace_=&workspace_;
  if(service.stratified_) {
    service.stratified_->eventWorkspace_=&workspace_;
    workspace_.metrics_=&service.stratified_->metrics_;
  } else if(service.barotropicQG_) {
    service.barotropicQG_->eventWorkspace_=&workspace_;
    workspace_.metrics_=&service.barotropicQG_->metrics_;
  } else workspace_.metrics_=&service.metrics_;
  status_=service.beginStateEvaluation(state,&workspace_);
  if(!status_) {
    release();
    return;
  }
  service.stateEvaluationActive_=true;
  kernelEvaluationActive_=true;
}
void WVFieldEvaluationEventScope::release() noexcept {
  if(!service_) return;
  auto variableMetrics=workspace_.evaluationMetrics();
  workspace_.endForcingEvaluation();
  if(kernelEvaluationActive_) {
    service_->endStateEvaluation();
    service_->stateEvaluationActive_=false;
    kernelEvaluationActive_=false;
  }
  workspace_.endEvaluation();
  service_->eventWorkspace_=nullptr;
  if(service_->stratified_) service_->stratified_->eventWorkspace_=nullptr;
  if(service_->barotropicQG_) service_->barotropicQG_->eventWorkspace_=nullptr;
  for(auto& field:workspace_.realFields_) {
    field->values.clear();
    field->assigned=false;
  }
  for(auto& field:workspace_.complexFields_) {
    field->values.clear();
    field->assigned=false;
  }
  for(auto& density:workspace_.density_) density.resetRetainingCapacity();
  for(auto& values:workspace_.apv_) values.clear();
  workspace_.apvReady_={};
  workspace_.densitySource_.clear();
  workspace_.densityHeights_.clear();
  workspace_.densityWeights_.clear();
  workspace_.densityInitial_.clear();
  workspace_.metrics_->densityWorkspaceLiveBytes=0;
  workspace_.metrics_->eventFieldWorkspaceLiveBytes=0;
  auto& aggregate=workspace_.metrics_->variableEvaluation;
  aggregate.contexts+=variableMetrics.contexts;
  aggregate.producerExecutions+=variableMetrics.producerExecutions;
  aggregate.cacheHits+=variableMetrics.cacheHits;
  aggregate.evictions+=variableMetrics.evictions;
  aggregate.recomputations+=variableMetrics.recomputations;
  aggregate.duplicateExecutions+=variableMetrics.duplicateExecutions;
  aggregate.liveBytes=0;
  aggregate.highWaterBytes=std::max(aggregate.highWaterBytes,
                                    variableMetrics.highWaterBytes);
  service_=nullptr;
}
} // namespace detail

WVFieldEvaluationService::WVFieldEvaluationService() {
  eventArena_=std::make_unique<detail::WVFieldEvaluationArena>();
  variableEvaluationPreparationStatus_=
      detail::WVFieldEvaluationEventWorkspace::prepareEvaluationContext(
          variableEvaluationContext_);
}

WVKernelStatus WVFieldEvaluationService::prepareEventArena(
    std::size_t requestCount) const {
  auto status=eventArena_->prepareGroupNodes(requestCount);
  if(status) status=eventArena_->preparePolicy(variableEvaluationPolicy_);
  if(status && forcing_ && eventArenaForcingPrepared_)
    status=eventArena_->prepareForcing(*forcing_,variableEvaluationPolicy_);
  auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
  metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
  metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
  metrics.servicePersistentBytes=persistentBytes();
  return status;
}

WVKernelStatus WVFieldEvaluationService::prepareEventField(
    const WVVariableEvaluationKey& key,std::size_t elements,bool complex) const {
  const auto status=complex ? eventArena_->prepareComplex(key,elements) :
      eventArena_->prepareReal(key,elements);
  if(status) {
    if(!complex && key.node==WVVariableEvaluationNode::derivative) {
      WVShape3D shape;
      if(stratified_) {
        const auto& geometry=stratified_->configuration();
        shape={geometry.Nx,geometry.Ny,geometry.Nz};
      } else if(barotropicQG_) {
        const auto& configuration=barotropicQG_->configuration();
        shape={configuration.Nx,configuration.Ny,1};
      } else shape=transform_->descriptor().spatialShape();
      eventArena_->setPreparedVolumeShape(key,shape);
    }
    auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
    metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
    metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
    metrics.servicePersistentBytes=persistentBytes();
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::prepareDensityEventArena(
    std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
    WVNoMotionReference reference,bool apvNeeded) const {
  const auto status=eventArena_->prepareDensity(
      sampleCount,profileCount,demands,reference,apvNeeded);
  if(status) {
    auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
    metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
    metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
    metrics.servicePersistentBytes=persistentBytes();
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::prepareEventArena(
    const WVFieldEvaluationPlan& plan,std::uint32_t componentIdentity) const {
  return prepareEventArenaImpl(plan,componentIdentity,false);
}

WVKernelStatus WVFieldEvaluationService::prepareEventArenaForActiveEvaluation(
    const WVFieldEvaluationPlan& plan,std::uint32_t componentIdentity) const {
  return prepareEventArenaImpl(plan,componentIdentity,true);
}

WVKernelStatus WVFieldEvaluationService::prepareEventArenaImpl(
    const WVFieldEvaluationPlan& plan,std::uint32_t componentIdentity,
    bool activePreparation) const {
  constexpr std::uint64_t primitiveValueMask=1ULL<<0;
  constexpr std::uint64_t pressureHeightMask=1ULL<<1;
  constexpr std::uint64_t streamfunctionMask=1ULL<<2;
  constexpr std::uint64_t potentialVorticityMask=1ULL<<3;
  constexpr std::uint64_t uDerivativeMask=1ULL<<4;
  constexpr std::uint64_t vDerivativeMask=1ULL<<5;
  constexpr std::uint64_t wDerivativeMask=1ULL<<6;
  if(plan.diagnosticPlan_) {
    if(componentIdentity)
      return {WVKernelStatusCode::invalidConfiguration,
          "A diagnostic plan cannot be nested as a component dependency."};
    return activePreparation ?
        plan.diagnosticPlan_->prepareEventArenaForActiveEvaluation(*this) :
        plan.diagnosticPlan_->prepareEventArena(*this);
  }
  std::size_t R=0;
  if(stratified_) {
    const auto& geometry=stratified_->configuration();
    R=geometry.Nx*geometry.Ny*geometry.Nz;
  } else if(barotropicQG_) {
    const auto& configuration=barotropicQG_->configuration();
    R=configuration.Nx*configuration.Ny;
  } else R=transform_->descriptor().spatialShape().elementCount();
  const auto prepareReal=[&](WVVariableEvaluationKey key,std::size_t elements) {
    return activePreparation ?
        prepareEventFieldForActiveEvaluation(key,elements,false) :
        prepareEventField(key,elements,false);
  };
  std::size_t fusedGroupSize=0;
  const auto finishPreparation=[&]() {
    if(!activePreparation) return prepareEventArena(fusedGroupSize);
    auto status=eventArena_->prepareGroupNodesForActive(fusedGroupSize);
    if(status && forcing_ && eventArenaForcingPrepared_)
      status=eventArena_->prepareForcingForActive(
          *forcing_,variableEvaluationPolicy_);
    if(status) {
      auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
      metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
      metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
      metrics.servicePersistentBytes=persistentBytes();
    }
    return status;
  };
  if(transform_) {
    if(plan.dependencyMask_&primitiveValueMask) {
      auto status=prepareReal({WVVariableEvaluationNode::physicalField,
          static_cast<std::uint32_t>(WVPortableVariable::u),componentIdentity,
          0,0,0,1},4*R);
      if(!status) return status;
    }
    const struct {std::uint64_t mask; WVPortableVariable field;} derivatives[]={{
        uDerivativeMask,WVPortableVariable::u},{vDerivativeMask,WVPortableVariable::v},
        {wDerivativeMask,WVPortableVariable::w}};
    for(const auto& dependency:derivatives) if(plan.dependencyMask_&dependency.mask) {
      fusedGroupSize=3;
      for(std::uint32_t axis=1;axis<=3;++axis) {
        const auto status=prepareReal({WVVariableEvaluationNode::derivative,
            static_cast<std::uint32_t>(dependency.field),componentIdentity,axis},R);
        if(!status) return status;
      }
    }
    const struct {std::uint64_t mask; WVPortableVariable field;} fused[]={{
        pressureHeightMask,WVPortableVariable::pi},{streamfunctionMask,WVPortableVariable::psi},
        {potentialVorticityMask,WVPortableVariable::qgpv}};
    for(const auto& dependency:fused) if(plan.dependencyMask_&dependency.mask) {
      const auto status=prepareReal({WVVariableEvaluationNode::physicalField,
          static_cast<std::uint32_t>(dependency.field),componentIdentity},4*R);
      if(!status) return status;
    }
    for(const auto& output:plan.outputs_) {
      const auto* metadata=findPortableVariable(output.fieldName);
      if(!metadata) continue;
      const auto field=metadata->identifier;
      if(field==WVPortableVariable::p || field==WVPortableVariable::rhoE ||
          field==WVPortableVariable::rhoTotal ||
          field==WVPortableVariable::zetaX ||
          field==WVPortableVariable::zetaY ||
          field==WVPortableVariable::zetaZ) {
        const auto status=prepareReal({WVVariableEvaluationNode::physicalField,
            static_cast<std::uint32_t>(field),componentIdentity},R);
        if(!status) return status;
      } else if(field==WVPortableVariable::energy) {
        const auto status=prepareReal({WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(field),componentIdentity},1);
        if(!status) return status;
      }
    }
    return finishPreparation();
  }
  for(const auto& output:plan.outputs_) {
    auto* metadata=findPortableVariable(output.fieldName);
    if(!metadata) continue;
    auto field=metadata->identifier;
    if(field==WVPortableVariable::ssu) field=WVPortableVariable::u;
    else if(field==WVPortableVariable::ssv) field=WVPortableVariable::v;
    else if(field==WVPortableVariable::ssh) field=WVPortableVariable::pi;
    if(field==WVPortableVariable::totalEnergySpatiallyIntegrated ||
        field==WVPortableVariable::energy) {
      const auto status=prepareReal({WVVariableEvaluationNode::reduction,
          static_cast<std::uint32_t>(field),
          field==WVPortableVariable::energy ? componentIdentity : 0u},1);
      if(!status) return status;
      continue;
    }
    if(field==WVPortableVariable::uvMax) {
      auto status=prepareReal({WVVariableEvaluationNode::reduction,
          static_cast<std::uint32_t>(field)},1);
      if(!status) return status;
      for(const auto velocity:{WVPortableVariable::u,WVPortableVariable::v}) {
        status=prepareReal({WVVariableEvaluationNode::physicalField,
            static_cast<std::uint32_t>(velocity),componentIdentity},R);
        if(!status) return status;
      }
      continue;
    }
    if(field==WVPortableVariable::wMax) {
      auto status=prepareReal({WVVariableEvaluationNode::reduction,
          static_cast<std::uint32_t>(field)},1);
      if(!status) return status;
      field=WVPortableVariable::w;
    }
    auto status=prepareReal({WVVariableEvaluationNode::physicalField,
        static_cast<std::uint32_t>(field),componentIdentity},R);
    if(!status) return status;
    if((field==WVPortableVariable::zetaX || field==WVPortableVariable::zetaY) &&
        stratified_ && (stratified_->configuration().transformClass==
          "WVTransformHydrostatic" || stratified_->configuration().transformClass==
          "WVTransformBoussinesq")) {
      const auto first=field==WVPortableVariable::zetaX ?
          WVPortableVariable::w : WVPortableVariable::u;
      const auto firstAxis=field==WVPortableVariable::zetaX ? 2u : 3u;
      const auto second=field==WVPortableVariable::zetaX ?
          WVPortableVariable::v : WVPortableVariable::w;
      const auto secondAxis=field==WVPortableVariable::zetaX ? 3u : 1u;
      status=prepareReal({WVVariableEvaluationNode::derivative,
          static_cast<std::uint32_t>(first),componentIdentity,firstAxis},R);
      if(!status) return status;
      status=prepareReal({WVVariableEvaluationNode::derivative,
          static_cast<std::uint32_t>(second),componentIdentity,secondAxis},R);
      if(!status) return status;
    }
  }
  return finishPreparation();
}

WVKernelStatus WVFieldEvaluationService::prepareEventArena(
    const WVEventFieldEvaluationPlan& plan) const {
  WVFieldEvaluationPlan fields;
  fields.dependencyMask_=plan.dependencyMask_;
  fields.requestedFieldMask_=plan.requestedFieldMask_;
  try {
    fields.outputs_.reserve(plan.outputs_.size());
    for(const auto& output:plan.outputs_)
      fields.outputs_.push_back({output.identifier,output.fieldName,
          WVFieldSamplingKind::fullGrid,{},0,false});
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,
        "Unable to prepare sampled output-event dependencies."};
  }
  return prepareEventArena(fields);
}

WVKernelStatus WVFieldEvaluationService::prepareEventFieldForActiveEvaluation(
    const WVVariableEvaluationKey& key,std::size_t elements,bool complex) const {
  const auto status=complex ? eventArena_->prepareComplexForActive(key,elements) :
      eventArena_->prepareRealForActive(key,elements);
  if(status) {
    if(!complex && key.node==WVVariableEvaluationNode::derivative) {
      WVShape3D shape;
      if(stratified_) {
        const auto& geometry=stratified_->configuration();
        shape={geometry.Nx,geometry.Ny,geometry.Nz};
      } else if(barotropicQG_) {
        const auto& configuration=barotropicQG_->configuration();
        shape={configuration.Nx,configuration.Ny,1};
      } else shape=transform_->descriptor().spatialShape();
      eventArena_->setPreparedVolumeShape(key,shape);
    }
    auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
    metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
    metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
    metrics.servicePersistentBytes=persistentBytes();
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::prepareDensityEventArenaForActiveEvaluation(
    std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
    WVNoMotionReference reference,bool apvNeeded) const {
  const auto status=eventArena_->prepareDensityForActive(sampleCount,profileCount,
      demands,reference,apvNeeded);
  if(status) {
    auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
    metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
    metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
    metrics.servicePersistentBytes=persistentBytes();
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::beginStateEvaluation(
    const WVIntegrationState& state,const void* owner) {
  if(stratified_) return stratified_->beginStateEvaluation(state,owner);
  if(barotropicQG_) return barotropicQG_->beginStateEvaluation(state,owner);
  return transform_->beginStateEvaluation(state.waveVortex,owner);
}

WVKernelStatus WVFieldEvaluationService::addStateEvaluationView(
    const WVIntegrationState& state,std::size_t componentIdentity) {
  if(!stateEvaluationActive_) return WVKernelStatus::ok();
  if(stratified_) return stratified_->addStateEvaluationView(
      state,eventWorkspace_,componentIdentity);
  if(barotropicQG_) return barotropicQG_->addStateEvaluationView(
      state,eventWorkspace_,componentIdentity);
  return transform_->addStateEvaluationView(
      state.waveVortex,eventWorkspace_,componentIdentity);
}

WVKernelStatus WVFieldEvaluationService::removeStateEvaluationView(
    const WVIntegrationState& state,std::size_t componentIdentity) {
  if(!stateEvaluationActive_ || !eventWorkspace_)
    return {WVKernelStatusCode::invalidConfiguration,
        "No field state evaluation is active."};
  if(stratified_) return stratified_->removeStateEvaluationView(
      state,eventWorkspace_,componentIdentity);
  if(barotropicQG_) return barotropicQG_->removeStateEvaluationView(
      state,eventWorkspace_,componentIdentity);
  return transform_->removeStateEvaluationView(
      state.waveVortex,eventWorkspace_,componentIdentity);
}

void WVFieldEvaluationService::endStateEvaluation() noexcept {
  if(stratified_) stratified_->endStateEvaluation();
  else if(barotropicQG_) barotropicQG_->endStateEvaluation();
  else (void)transform_->endStateEvaluation();
}

WVKernelStatus WVFieldEvaluationService::prepareForcingEvaluationContext() {
  std::vector<WVVariableEvaluationKey> keys;
  try {
    keys.reserve(forcingVariableBindings().size()*2);
    for(std::uint32_t index=0;index<forcingVariableBindings().size();++index) {
      keys.push_back({WVVariableEvaluationNode::forcingTendency,
          static_cast<std::uint32_t>(WVPortableVariable::Fu_portable_catalog_forcing),
          0,0,0,index,1});
      keys.push_back({WVVariableEvaluationNode::forcingTendency,
          static_cast<std::uint32_t>(WVPortableVariable::Fqgpv_portable_catalog_forcing),
          0,0,0,index,1});
    }
    const auto dependencies=forcing_->dependencyKeys();
    keys.insert(keys.end(),dependencies.begin(),dependencies.end());
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,
        "Unable to prepare forcing variable evaluation keys."};
  }
  return variableEvaluationContext_.prepare(keys);
}

class WVFieldEvaluationSession::Impl final {
public:
  Impl(WVFieldEvaluationService& service,const WVIntegrationState& state)
      : scope(service,state,true,true) {}
  detail::WVFieldEvaluationEventScope scope;
};

WVFieldEvaluationSession::WVFieldEvaluationSession()=default;
WVFieldEvaluationSession::~WVFieldEvaluationSession()=default;
WVFieldEvaluationSession::WVFieldEvaluationSession(WVFieldEvaluationSession&&) noexcept=default;
WVFieldEvaluationSession& WVFieldEvaluationSession::operator=(WVFieldEvaluationSession&&) noexcept=default;
bool WVFieldEvaluationSession::active() const noexcept {
  return impl_ && static_cast<bool>(impl_->scope.status());
}
const WVKernelStatus& WVFieldEvaluationSession::status() const noexcept {
  static const WVKernelStatus inactive{WVKernelStatusCode::invalidConfiguration,
      "Field evaluation session is inactive."};
  return impl_ ? impl_->scope.status() : inactive;
}

WVKernelStatus WVFieldEvaluationService::validateVariableEvaluationPolicyChange(
    WVVariableEvaluationPolicy policy) const noexcept {
  if(eventWorkspace_)
    return {WVKernelStatusCode::invalidConfiguration,
            "Cannot change variable evaluation policy during an active session."};
  if(policy!=WVVariableEvaluationPolicy::reuse &&
      policy!=WVVariableEvaluationPolicy::lowMemory)
    return {WVKernelStatusCode::invalidConfiguration,
            "Unknown variable evaluation policy."};
  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::setVariableEvaluationPolicy(
    WVVariableEvaluationPolicy policy) {
  auto validation=validateVariableEvaluationPolicyChange(policy);
  if(!validation) return validation;
  if(policy==variableEvaluationPolicy_) return WVKernelStatus::ok();
  if(policy==WVVariableEvaluationPolicy::lowMemory) {
    if(forcing_ && eventArenaForcingPrepared_) {
      const auto status=eventArena_->prepareForcing(*forcing_,policy);
      if(!status) return status;
    }
    eventArena_->clearReal();
    variableEvaluationPolicy_=policy;
  } else {
    variableEvaluationPolicy_=policy;
    const auto status=prepareEventArena(0);
    if(!status) {
      variableEvaluationPolicy_=WVVariableEvaluationPolicy::lowMemory;
      eventArena_->clearReal();
      if(forcing_ && eventArenaForcingPrepared_) (void)eventArena_->prepareForcing(
          *forcing_,WVVariableEvaluationPolicy::lowMemory);
      return status;
    }
  }
  auto& metrics=mutableMetrics();
  metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
  metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
  metrics.servicePersistentBytes=persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::beginEvaluationSession(
    const WVIntegrationState& state,WVFieldEvaluationSession& session) {
  if(!variableEvaluationPreparationStatus_)
    return variableEvaluationPreparationStatus_;
  if(session.impl_)
    return {WVKernelStatusCode::invalidConfiguration,
            "Field evaluation session is already active."};
  if(eventWorkspace_)
    return {WVKernelStatusCode::reentrantExecution,
            "A field evaluation session is already active."};
  try {
    auto candidate=std::make_unique<WVFieldEvaluationSession::Impl>(*this,state);
    const auto status=candidate->scope.status();
    if(!status) return status;
    session.impl_=std::move(candidate);
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to create a field evaluation session."};
  }
}

WVKernelStatus WVFieldEvaluationService::prepareBuiltinNonlinearCoefficientEvaluation() {
  if(builtinNonlinearPrepared_) return WVKernelStatus::ok();
  if(eventWorkspace_)
    return {WVKernelStatusCode::reentrantExecution,
        "Built-in nonlinear coefficient preparation requires an idle service."};
  if(!forcing_ || !forcing_->supportsExactBuiltinNonlinear())
    return {WVKernelStatusCode::unsupportedOperation,
        "The field service does not own the exact built-in nonlinear forcing schedule."};
  WVFieldEvaluationPlan candidate;
  const std::vector<WVFieldRequest> requests{{"builtin-u","u",{}},
      {"builtin-v","v",{}},{"builtin-w","w",{}},{"builtin-eta","eta",{}}};
  auto status=createPlan(requests,candidate); if(!status) return status;
  std::size_t R=0;
  std::size_t S=0;
  if(stratified_) {const auto& g=stratified_->configuration(); R=g.Nx*g.Ny*g.Nz; S=g.Nj*g.Nkl;}
  else {const auto& g=transform_->descriptor().configuration(); R=g.Nx*g.Ny*g.Nz;
    S=transform_->descriptor().spectralShape().elementCount();}
  try {builtinNonlinearPhysical_.resize(4*R); builtinNonlinearCoefficients_.resize(3*S);}
  catch(const std::bad_alloc&) {return {WVKernelStatusCode::allocationFailure,
      "Unable to prepare built-in nonlinear physical storage."};}
  const bool previous=eventArenaForcingPrepared_;
  eventArenaForcingPrepared_=true;
  status=prepareEventArena(3);
  if(!status) {eventArenaForcingPrepared_=previous; return status;}
  builtinNonlinearFieldsPlan_=std::move(candidate);
  builtinNonlinearPrepared_=true;
  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::evaluateBuiltinNonlinearCoefficients(
    const WVIntegrationState& state,WVFlux& flux) {
  if(!builtinNonlinearPrepared_ || !forcing_ ||
      !forcing_->supportsExactBuiltinNonlinear())
    return {WVKernelStatusCode::unsupportedOperation,
        "Built-in nonlinear coefficient evaluation was not prepared."};
  if(!eventWorkspace_ || !stateEvaluationActive_ ||
      variableEvaluationPolicy_!=WVVariableEvaluationPolicy::reuse)
    return {WVKernelStatusCode::invalidConfiguration,
        "Built-in nonlinear coefficients require an active reuse-policy evaluation."};
  if(executing_ || !eventWorkspace_->preparationIdle())
    return {WVKernelStatusCode::reentrantExecution,
        "Built-in nonlinear coefficients require an idle field producer."};
  auto status=eventWorkspace_->validateState(state); if(!status) return status;
  const WVShape2D expected=stratified_ ?
      WVShape2D{stratified_->configuration().Nj,stratified_->configuration().Nkl} :
      transform_->descriptor().spectralShape();
  const WVComplexView coefficientOutputs[]{flux.Fp,flux.Fm,flux.F0};
  const WVComplexConstView coefficientInputs[]{state.waveVortex.coefficients.Ap,
      state.waveVortex.coefficients.Am,state.waveVortex.coefficients.A0};
  const auto coefficientBytes=expected.elementCount()*sizeof(WVComplex64);
  const auto addressFits=[&](const void* pointer) {
    const auto address=reinterpret_cast<std::uintptr_t>(pointer);
    return address && address%alignof(WVComplex64)==0 &&
        coefficientBytes<=std::numeric_limits<std::uintptr_t>::max()-address;
  };
  const auto overlaps=[](const void* first,std::size_t firstBytes,
      const void* second,std::size_t secondBytes) {
    const auto a=reinterpret_cast<std::uintptr_t>(first);
    const auto b=reinterpret_cast<std::uintptr_t>(second);
    return a<=b ? b-a<firstBytes : a-b<secondBytes;
  };
  for(const auto input:coefficientInputs)
    if(input.shape.rows!=expected.rows || input.shape.columns!=expected.columns)
      return {WVKernelStatusCode::invalidShape,
          "Built-in nonlinear coefficient state has the wrong shape."};
  for(std::size_t output=0;output<3;++output) {
    const auto view=coefficientOutputs[output];
    if(view.shape.rows!=expected.rows ||
        view.shape.columns!=expected.columns)
      return {WVKernelStatusCode::invalidShape,
          "Built-in nonlinear coefficient output has the wrong shape."};
    if(!addressFits(view.data))
      return {WVKernelStatusCode::invalidPointer,
          "Built-in nonlinear coefficient output has an invalid address."};
    for(const auto input:coefficientInputs)
      if(overlaps(view.data,coefficientBytes,input.data,coefficientBytes))
        return {WVKernelStatusCode::overlappingArrays,
            "Built-in nonlinear coefficient output overlaps the immutable state."};
    for(std::size_t previous=0;previous<output;++previous)
      if(overlaps(view.data,coefficientBytes,
          coefficientOutputs[previous].data,coefficientBytes))
        return {WVKernelStatusCode::overlappingArrays,
            "Built-in nonlinear coefficient outputs overlap."};
  }
  std::size_t R=0;
  if(stratified_) {const auto& g=stratified_->configuration(); R=g.Nx*g.Ny*g.Nz;}
  else {const auto& g=transform_->descriptor().configuration(); R=g.Nx*g.Ny*g.Nz;}
  if(builtinNonlinearPhysical_.size()<4*R)
    return {WVKernelStatusCode::invalidConfiguration,
        "Built-in nonlinear physical workspace was not prepared."};
  WVFieldOutputView outputs[]={{builtinNonlinearPhysical_.data(),R},
      {builtinNonlinearPhysical_.data()+R,R},{builtinNonlinearPhysical_.data()+2*R,R},
      {builtinNonlinearPhysical_.data()+3*R,R}};
  detail::WVForcingDiagnosticWorkspace* workspace=nullptr;
  status=eventWorkspace_->forcingWorkspace(*forcing_,workspace); if(!status) return status;
  const WVVariableEvaluationKey primitiveKey{WVVariableEvaluationNode::physicalField,
      static_cast<std::uint32_t>(WVPortableVariable::u),eventWorkspace_->component(),0,0,0,1};
  const bool physicalReady=eventWorkspace_->ready(primitiveKey);
  bool individualPhysicalReady=false;
  for(const auto variable:{WVPortableVariable::u,WVPortableVariable::v,
      WVPortableVariable::w,WVPortableVariable::eta}) {
    const WVVariableEvaluationKey key{WVVariableEvaluationNode::physicalField,
        static_cast<std::uint32_t>(variable),eventWorkspace_->component()};
    individualPhysicalReady|=eventWorkspace_->ready(key);
  }
  if(transform_ || physicalReady || individualPhysicalReady) {
    status=evaluate(builtinNonlinearFieldsPlan_,state,outputs,4); if(!status) return status;
  }
  bool preparedPhysical=transform_ || physicalReady;
  if(!transform_ && !physicalReady && individualPhysicalReady) {
    bool reused=false;
    status=eventWorkspace_->evaluate(primitiveKey,builtinNonlinearPhysical_.data(),4*R,
        [] {return WVKernelStatus::ok();},reused);
    if(!status) return status;
    preparedPhysical=true;
  }
  executing_=true;
  struct ResetExecution {bool& executing; ~ResetExecution(){executing=false;}}
      resetExecution{executing_};
  const std::size_t channels=transform_ ? 3u : 4u;
  const auto Nx=stratified_ ? stratified_->configuration().Nx :
      transform_->descriptor().configuration().Nx;
  const auto Ny=stratified_ ? stratified_->configuration().Ny :
      transform_->descriptor().configuration().Ny;
  const auto Nz=stratified_ ? stratified_->configuration().Nz :
      transform_->descriptor().configuration().Nz;
  const WVRealFieldBundleConstView prepared{builtinNonlinearPhysical_.data(),
      {Nx,Ny,Nz,channels}};
  const auto S=expected.elementCount();
  WVFlux staged{{builtinNonlinearCoefficients_.data(),expected},
      {builtinNonlinearCoefficients_.data()+S,expected},
      {builtinNonlinearCoefficients_.data()+2*S,expected}};
  status=forcing_->evaluateBuiltinNonlinear(state.waveVortex,
      preparedPhysical ? &prepared : nullptr,workspace,staged);
  if(!status) return status;
  if(!transform_ && !preparedPhysical) {
    std::copy_n(workspace->physical.data(),4*R,builtinNonlinearPhysical_.data());
    bool reused=false;
    status=eventWorkspace_->evaluate(primitiveKey,builtinNonlinearPhysical_.data(),4*R,
        [] {return WVKernelStatus::ok();},reused);
  }
  if(!status) return status;
  for(std::size_t family=0;family<3;++family)
    std::copy_n(builtinNonlinearCoefficients_.data()+family*S,S,
        coefficientOutputs[family].data);
  return WVKernelStatus::ok();
}

bool WVFieldEvaluationPlan::hasDensityDiagnostics() const noexcept {
  return diagnosticPlan_ && diagnosticPlan_->hasDensityDiagnostics();
}
namespace {

enum Dependency : std::uint64_t {
  primitiveValues = 1ULL << 0,
  pressureHeight = 1ULL << 1,
  streamfunction = 1ULL << 2,
  potentialVorticity = 1ULL << 3,
  uDerivatives = 1ULL << 4,
  vDerivatives = 1ULL << 5,
  wDerivatives = 1ULL << 6,
  spectralEnergy = 1ULL << 7
};

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

constexpr std::uint64_t fingerprintOffset = 1469598103934665603ULL;
constexpr std::uint64_t fingerprintPrime = 1099511628211ULL;

void appendFingerprint(std::uint64_t &fingerprint, const void *data,
                       std::size_t byteCount) noexcept {
  const auto *bytes = static_cast<const std::uint8_t *>(data);
  for (std::size_t index = 0; index < byteCount; ++index) {
    fingerprint ^= bytes[index];
    fingerprint *= fingerprintPrime;
  }
}

template <typename Value>
void appendFingerprint(std::uint64_t &fingerprint,
                       const Value &value) noexcept {
  appendFingerprint(fingerprint, &value, sizeof(value));
}

void appendConfigurationFingerprint(
    std::uint64_t &fingerprint,
    const WVTransformConstantStratificationConfiguration &configuration)
    noexcept {
  appendFingerprint(fingerprint, configuration.contractVersion);
  appendFingerprint(fingerprint, configuration.Nx);
  appendFingerprint(fingerprint, configuration.Ny);
  appendFingerprint(fingerprint, configuration.Nz);
  appendFingerprint(fingerprint, configuration.Nj);
  appendFingerprint(fingerprint, configuration.Lx);
  appendFingerprint(fingerprint, configuration.Ly);
  appendFingerprint(fingerprint, configuration.Lz);
  appendFingerprint(fingerprint, configuration.N0);
  appendFingerprint(fingerprint, configuration.rho0);
  appendFingerprint(fingerprint, configuration.g);
  appendFingerprint(fingerprint, configuration.planetaryRadius);
  appendFingerprint(fingerprint, configuration.rotationRate);
  appendFingerprint(fingerprint, configuration.latitude);
  appendFingerprint(fingerprint, configuration.isHydrostatic);
  appendFingerprint(fingerprint, configuration.shouldAntialias);
}

std::size_t checkedProduct(std::size_t first, std::size_t second) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first)
    throw std::overflow_error("field-evaluation size overflow");
  return first * second;
}

bool memoryOverlaps(const void *first, std::size_t firstBytes,
                    const void *second, std::size_t secondBytes) {
  if (!first || !second || firstBytes == 0 || secondBytes == 0)
    return false;
  const auto firstAddress = reinterpret_cast<std::uintptr_t>(first);
  const auto secondAddress = reinterpret_cast<std::uintptr_t>(second);
  return firstAddress < secondAddress + secondBytes &&
         secondAddress < firstAddress + firstBytes;
}

WVComplex64 subtract(WVComplex64 first, WVComplex64 second) noexcept {
  return {first.real - second.real, first.imag - second.imag};
}

WVComplex64 multiply(WVComplex64 first, WVComplex64 second) noexcept {
  return {first.real * second.real - first.imag * second.imag,
          first.real * second.imag + first.imag * second.real};
}

WVComplex64 multiply(WVComplex64 value, double scale) noexcept {
  return {value.real * scale, value.imag * scale};
}

WVComplex64 conjugate(WVComplex64 value) noexcept {
  return {value.real, -value.imag};
}

double squaredMagnitude(WVComplex64 value) noexcept {
  return value.real * value.real + value.imag * value.imag;
}

double wrapped(double coordinate, double length) noexcept {
  double value = std::fmod(coordinate, length);
  if (value < 0.0)
    value += length;
  if (value >= length)
    value = 0.0;
  return value;
}

class SplineSystem final {
public:
  explicit SplineSystem(std::size_t count) : count_(count) {
    if (count < 2)
      throw std::invalid_argument(
          "Spline interpolation requires at least two points per axis.");
    if (count < 4)
      return;
    lu_.assign(checkedProduct(count, count), 0.0);
    pivots_.resize(count);
    // Factor the transpose of the uniform-grid not-a-knot system. The
    // common grid-spacing factor cancels between the system and its RHS.
    auto setSystem = [&](std::size_t row, std::size_t column, double value) {
      lu_[column + count * row] = value;
    };
    setSystem(0, 0, -1.0);
    setSystem(0, 1, 2.0);
    setSystem(0, 2, -1.0);
    for (std::size_t row = 1; row + 1 < count; ++row) {
      setSystem(row, row - 1, 1.0);
      setSystem(row, row, 4.0);
      setSystem(row, row + 1, 1.0);
    }
    setSystem(count - 1, count - 3, -1.0);
    setSystem(count - 1, count - 2, 2.0);
    setSystem(count - 1, count - 1, -1.0);
    std::vector<double> transpose(checkedProduct(count, count));
    for (std::size_t row = 0; row < count; ++row)
      for (std::size_t column = 0; column < count; ++column)
        transpose[column + count * row] = lu_[row + count * column];
    lu_.swap(transpose);
    factor();
  }

  std::vector<double> weights(double origin, double spacing, double query,
                              std::size_t circularShift = 0) const {
    std::vector<double> result;
    weightsInto(origin, spacing, query, result, circularShift);
    return result;
  }

  void weightsInto(double origin, double spacing, double query,
                   std::vector<double> &result,
                   std::size_t circularShift = 0,
                   std::vector<double> *rightHandSideWorkspace = nullptr,
                   std::vector<double> *shiftedWorkspace = nullptr) const {
    const double normalized = (query - origin) / spacing;
    auto restoreOriginalOrdering = [&](const std::vector<double> &shifted) {
      result.assign(count_, 0.0);
      for (std::size_t shiftedIndex = 0; shiftedIndex < count_;
           ++shiftedIndex) {
        const auto originalIndex =
            (shiftedIndex + count_ - circularShift % count_) % count_;
        result[originalIndex] += shifted[shiftedIndex];
      }
    };
    if (count_ == 2) {
      const std::vector<double> linear{1.0 - normalized, normalized};
      restoreOriginalOrdering(linear);
      return;
    }
    if (count_ == 3) {
      std::vector<double> quadratic = {
          (normalized - 1.0) * (normalized - 2.0) / 2.0,
          -normalized * (normalized - 2.0),
          normalized * (normalized - 1.0) / 2.0};
      restoreOriginalOrdering(quadratic);
      return;
    }
    std::size_t interval = normalized <= 0.0
                               ? 0
                               : static_cast<std::size_t>(std::floor(normalized));
    interval = std::min(interval, count_ - 2);
    const double fraction = std::clamp(normalized - static_cast<double>(interval),
                                       0.0, 1.0);
    const double first = 1.0 - fraction;
    const double second = fraction;
    std::vector<double> localRightHandSide;
    auto &rhs = rightHandSideWorkspace == nullptr ? localRightHandSide
                                                   : *rightHandSideWorkspace;
    rhs.assign(count_, 0.0);
    rhs[interval] = (first * first * first - first) * spacing * spacing / 6.0;
    rhs[interval + 1] =
        (second * second * second - second) * spacing * spacing / 6.0;
    solve(rhs);
    std::vector<double> localShifted;
    auto &shifted = shiftedWorkspace == nullptr ? localShifted
                                                : *shiftedWorkspace;
    shifted.assign(count_, 0.0);
    shifted[interval] += first;
    shifted[interval + 1] += second;
    const double rhsScale = 6.0 / (spacing * spacing);
    for (std::size_t row = 1; row + 1 < count_; ++row) {
      shifted[row - 1] += rhsScale * rhs[row];
      shifted[row] -= 2.0 * rhsScale * rhs[row];
      shifted[row + 1] += rhsScale * rhs[row];
    }
    restoreOriginalOrdering(shifted);
  }

  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + lu_.capacity() * sizeof(double) +
           pivots_.capacity() * sizeof(std::size_t);
  }

private:
  void factor() {
    for (std::size_t column = 0; column < count_; ++column) {
      std::size_t pivot = column;
      double maximum = std::abs(lu_[column + count_ * column]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        const double candidate = std::abs(lu_[column + count_ * row]);
        if (candidate > maximum) {
          maximum = candidate;
          pivot = row;
        }
      }
      if (maximum == 0.0)
        throw std::invalid_argument("Spline interpolation system is singular.");
      pivots_[column] = pivot;
      if (pivot != column)
        for (std::size_t entry = 0; entry < count_; ++entry)
          std::swap(lu_[entry + count_ * column],
                    lu_[entry + count_ * pivot]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        lu_[column + count_ * row] /= lu_[column + count_ * column];
        const double multiplier = lu_[column + count_ * row];
        for (std::size_t entry = column + 1; entry < count_; ++entry)
          lu_[entry + count_ * row] -=
              multiplier * lu_[entry + count_ * column];
      }
    }
  }

  void solve(std::vector<double> &rightHandSide) const {
    for (std::size_t column = 0; column < count_; ++column) {
      if (pivots_[column] != column)
        std::swap(rightHandSide[column], rightHandSide[pivots_[column]]);
      for (std::size_t row = column + 1; row < count_; ++row)
        rightHandSide[row] -=
            lu_[column + count_ * row] * rightHandSide[column];
    }
    for (std::size_t reverse = 0; reverse < count_; ++reverse) {
      const std::size_t row = count_ - 1 - reverse;
      for (std::size_t column = row + 1; column < count_; ++column)
        rightHandSide[row] -=
            lu_[column + count_ * row] * rightHandSide[column];
      rightHandSide[row] /= lu_[row + count_ * row];
    }
  }

  std::size_t count_ = 0;
  std::vector<double> lu_;
  std::vector<std::size_t> pivots_;
};

class ExecutionGuard final {
public:
  explicit ExecutionGuard(bool &executing)
      : executing_(executing), entered_(!executing) {
    if (entered_)
      executing_ = true;
  }
  ~ExecutionGuard() {
    if (entered_)
      executing_ = false;
  }
  bool entered() const noexcept { return entered_; }

private:
  bool &executing_;
  bool entered_ = false;
};

} // namespace

class WVFieldEvaluationService::MovingWorkspace final {
public:
  explicit MovingWorkspace(
      const WVTransformConstantStratificationConfiguration &configuration)
      : xSpline(configuration.Nx), ySpline(configuration.Ny),
        zSpline(configuration.Nz) {
    xWeights.reserve(configuration.Nx);
    yWeights.reserve(configuration.Ny);
    zWeights.reserve(configuration.Nz);
    xRightHandSide.reserve(configuration.Nx);
    yRightHandSide.reserve(configuration.Ny);
    zRightHandSide.reserve(configuration.Nz);
    xShifted.reserve(configuration.Nx);
    yShifted.reserve(configuration.Ny);
    zShifted.reserve(configuration.Nz);
  }
  SplineSystem xSpline;
  SplineSystem ySpline;
  SplineSystem zSpline;
  std::vector<double> xWeights;
  std::vector<double> yWeights;
  std::vector<double> zWeights;
  std::vector<double> xRightHandSide;
  std::vector<double> yRightHandSide;
  std::vector<double> zRightHandSide;
  std::vector<double> xShifted;
  std::vector<double> yShifted;
  std::vector<double> zShifted;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) +
           xSpline.persistentBytes() - sizeof(SplineSystem) +
           ySpline.persistentBytes() - sizeof(SplineSystem) +
           zSpline.persistentBytes() - sizeof(SplineSystem) +
           (xWeights.capacity() + yWeights.capacity() + zWeights.capacity() +
            xRightHandSide.capacity() + yRightHandSide.capacity() +
            zRightHandSide.capacity() + xShifted.capacity() +
            yShifted.capacity() + zShifted.capacity()) *
               sizeof(double);
  }
};

namespace detail {
class WVSampledMovingFieldPlan final {
public:
  struct Request {
    std::string identifier;
    std::string fieldName;
    WVPortableNaturalRank naturalRank = WVPortableNaturalRank::volume;
    std::size_t positionOffset = 0;
    std::size_t positionCount = 0;
    WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
    std::size_t outputIndex = 0;
  };
  std::string configurationIdentifier;
  const WVFieldEvaluationService *owner = nullptr;
  WVDensityDiagnosticContract densityContract;
  WVFieldEvaluationPlan fullGridPlan;
  std::vector<Request> requests;
  std::size_t persistentBytes() const noexcept {
    std::size_t bytes = sizeof(*this) + configurationIdentifier.capacity() +
                        requests.capacity() * sizeof(Request) +
                        fullGridPlan.persistentBytes() - sizeof(fullGridPlan);
    for (const auto &request : requests)
      bytes += request.identifier.capacity() + request.fieldName.capacity();
    return bytes;
  }
};
} // namespace detail

class WVFieldEvaluationService::SampledMovingWorkspace final {
public:
  WVKernelStatus prepare(const WVFieldEvaluationPlan& fullGridPlan,
      const std::vector<detail::WVSampledMovingFieldPlan::Request>& requests) {
    try {
      const auto count=fullGridPlan.outputCount();
      fullStorage.resize(count);
      fullViews.resize(count);
      selection.resize(count);
      sampledStorage.resize(count);
      for(const auto& request:requests) {
        fullStorage[request.outputIndex].reserve(
            fullGridPlan.outputs()[request.outputIndex].elementCount);
        sampledStorage[request.outputIndex].reserve(request.positionCount);
      }
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare sampled moving-field storage."};
    }
  }
  std::vector<std::vector<double>> fullStorage;
  std::vector<WVFieldOutputView> fullViews;
  std::vector<std::uint8_t> selection;
  std::vector<std::vector<double>> sampledStorage;
  std::size_t persistentBytes() const noexcept {
    std::size_t bytes=sizeof(*this)+
        fullStorage.capacity()*sizeof(std::vector<double>)+
        fullViews.capacity()*sizeof(WVFieldOutputView)+
        selection.capacity()*sizeof(std::uint8_t)+
        sampledStorage.capacity()*sizeof(std::vector<double>);
    for(const auto& values:fullStorage) bytes+=values.capacity()*sizeof(double);
    for(const auto& values:sampledStorage) bytes+=values.capacity()*sizeof(double);
    return bytes;
  }
};

std::size_t WVMovingFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      requests_.capacity() * sizeof(ResolvedRequest) +
                      outputs_.capacity() * sizeof(WVFieldOutputSpecification);
  for (const auto &output : outputs_)
    bytes += output.identifier.capacity() + output.fieldName.capacity() +
             output.dimensions.capacity() * sizeof(std::size_t);
  return bytes + transformPlanBytes_ +
         (sampledPlan_ ? sampledPlan_->persistentBytes() : 0);
}

std::size_t WVEventFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t bytes =
      sizeof(*this) + requests_.capacity() * sizeof(ResolvedRequest) +
      outputs_.capacity() * sizeof(WVEventFieldOutputSpecification) +
      requiresZByPositionSet_.capacity() * sizeof(std::uint8_t);
  for (const auto &output : outputs_)
    bytes += output.identifier.capacity() + output.fieldName.capacity();
  return bytes + transformPlanBytes_ + configurationIdentifier_.capacity() +
         (planIdentity_ ? sizeof(std::uint8_t) : 0);
}

std::size_t WVFieldEvaluationPlan::PositionWeights::persistentBytes() const
    noexcept {
  return xSplineWeights.capacity() * sizeof(double) +
         ySplineWeights.capacity() * sizeof(double) +
         zSplineWeights.capacity() * sizeof(double);
}

std::size_t WVFieldEvaluationPlan::ResolvedRequest::persistentBytes() const
    noexcept {
  std::size_t value = profileXIndices.capacity() * sizeof(std::size_t) +
                      profileYIndices.capacity() * sizeof(std::size_t) +
                      positionWeights.capacity() * sizeof(PositionWeights);
  for (const auto &weights : positionWeights)
    value += weights.persistentBytes();
  return value;
}

std::size_t WVFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t value = sizeof(*this) +
                      requests_.capacity() * sizeof(ResolvedRequest) +
                      outputs_.capacity() * sizeof(WVFieldOutputSpecification);
  for (const auto &request : requests_)
    value += request.persistentBytes();
  for (const auto &output : outputs_) {
    value += output.identifier.capacity() + output.fieldName.capacity() +
             output.dimensions.capacity() * sizeof(std::size_t);
  }
  return value + transformPlanBytes_ + (diagnosticPlan_ ? diagnosticPlan_->persistentBytes() : 0);
}

WVEventPositionSetView
WVPreparedFieldGeometry::positionSet(std::size_t slot) const noexcept {
  if (slot >= positionSets_.size())
    return {};
  const auto &set = positionSets_[slot];
  return {set.x, set.y, set.z, set.positionCount, set.extents.data(),
          set.extents.size()};
}

bool WVPreparedFieldGeometry::sameGeometry(
    const WVPreparedFieldGeometry &other) const noexcept {
  if (fieldPlanFingerprint_ != other.fieldPlanFingerprint_ ||
      geometryFingerprint_ != other.geometryFingerprint_ ||
      positionSets_.size() != other.positionSets_.size())
    return false;
  for (std::size_t slot = 0; slot < positionSets_.size(); ++slot) {
    const auto &first = positionSets_[slot];
    const auto &second = other.positionSets_[slot];
    if (first.positionCount != second.positionCount ||
        first.extents != second.extents ||
        (first.z == nullptr) != (second.z == nullptr))
      return false;
    const auto bytes = first.positionCount * sizeof(double);
    if (bytes != 0 &&
        (std::memcmp(first.x, second.x, bytes) != 0 ||
         std::memcmp(first.y, second.y, bytes) != 0 ||
         (first.z != nullptr &&
          std::memcmp(first.z, second.z, bytes) != 0)))
      return false;
  }
  return true;
}

std::size_t WVPreparedFieldGeometry::retainedBytes() const noexcept {
  std::size_t bytes =
      sizeof(*this) + positionSets_.capacity() * sizeof(PositionSet) +
      outputs_.capacity() * sizeof(WVPreparedFieldOutputSpecification);
  for (const auto &set : positionSets_)
    bytes += set.extents.capacity() * sizeof(std::size_t);
  for (const auto &output : outputs_)
    bytes += output.dimensions.capacity() * sizeof(std::size_t);
  const auto evaluationBytes = evaluationPlan_.persistentBytes();
  if (evaluationBytes >= sizeof(evaluationPlan_))
    bytes += evaluationBytes - sizeof(evaluationPlan_);
  return bytes + transformGeometryBytes_;
}

std::size_t WVPreparedFieldGeometry::liveBytes() const noexcept {
  return retainedBytes() + borrowedCoordinateBytes_;
}

WVPreparedFieldGeometryMetrics
WVPreparedFieldGeometry::metrics() const noexcept {
  return {positionSets_.size(), positionCount_, retainedBytes(), liveBytes()};
}

WVKernelStatus WVFieldEvaluationService::create(
    const WVTransformConstantStratificationConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status = WVTransformConstantStratificationKernel::create(
        configuration, std::move(engine), candidate->ownedTransform_);
    if (!status)
      return status;
    candidate->transform_ = candidate->ownedTransform_.get();
    status = candidate->initializeScratch();
    if (!status)
      return status;
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate bounded field-evaluation storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformConstantStratificationKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    candidate->transform_ = &transform;
    const auto status = candidate->initializeScratch();
    if (!status)
      return status;
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed field-evaluation storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  }
}

WVKernelStatus WVFieldEvaluationService::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status = detail::WVBarotropicQGFieldEvaluationAdapter::create(
        configuration, std::move(engine), candidate->barotropicQG_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->barotropicQG_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Barotropic QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformBarotropicQGKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status =
        detail::WVBarotropicQGFieldEvaluationAdapter::createBorrowing(
            transform, candidate->barotropicQG_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->barotropicQG_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the borrowed Barotropic QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::initializeScratch() {
  const auto fieldElements = transform_->descriptor().spatialShape().elementCount();
  const auto coefficientElements =
      transform_->descriptor().spectralShape().elementCount();
  realScratch_.resize(checkedProduct(6, fieldElements));
  complexScratch_.resize(checkedProduct(2, coefficientElements));
  movingWorkspace_ = std::make_unique<MovingWorkspace>(
      transform_->descriptor().configuration());
  metrics_.transformPersistentBytes = transform_->persistentBytes();
  metrics_.scratchCapacityBytes =
      realScratch_.capacity() * sizeof(double) +
      complexScratch_.capacity() * sizeof(WVComplex64);
  metrics_.movingInterpolationWorkspaceBytes =
      movingWorkspace_->persistentBytes();
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVFieldEvaluationService::~WVFieldEvaluationService() = default;
WVKernelStatus WVFieldEvaluationService::create(
    std::shared_ptr<const WVStratifiedModalSource> source,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status = detail::WVStratifiedFieldEvaluationAdapter::create(
        source, std::move(engine), candidate->stratified_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->stratified_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Stratified QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformStratifiedQGKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status =
        detail::WVStratifiedFieldEvaluationAdapter::createBorrowing(
            transform, candidate->stratified_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->stratified_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the borrowed Stratified QG field boundary."};
  }
}
WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformHydrostaticKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status =
        detail::WVStratifiedFieldEvaluationAdapter::createBorrowing(
            transform, candidate->stratified_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->stratified_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the borrowed Stratified QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformBoussinesqKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status =
        detail::WVStratifiedFieldEvaluationAdapter::createBorrowing(
            transform, candidate->stratified_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->stratified_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the borrowed Boussinesq field boundary."};
  }
}

std::vector<std::string> WVFieldEvaluationService::supportedFieldNames() {
  std::vector<std::string> result;
  result.reserve(WVPortableVariableCatalog.size());
  for (const auto &variable : WVPortableVariableCatalog)
    if (variable.kind != WVPortableVariableKind::coefficient && findExecutablePortableVariable(variable.name))
      result.emplace_back(variable.name);
  return result;
}

const std::vector<WVPortableForcingVariableBinding>& WVFieldEvaluationService::forcingVariableBindings() const noexcept {
  static const std::vector<WVPortableForcingVariableBinding> empty;
  return forcing_ ? forcing_->bindings() : empty;
}
std::string WVFieldEvaluationService::portableVariableConfiguration() const {
  return detail::WVDiagnosticFieldPlan::configurationIdentifier(*this);
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(WVConstantStratificationForcingEngine& engine,std::unique_ptr<WVFieldEvaluationService>& result) {
  std::unique_ptr<WVFieldEvaluationService> candidate;
  auto status=createBorrowing(engine.kernel(),candidate); if(!status) return status;
  status=detail::WVForcingDiagnosticBinding::create<WVConstantStratificationForcingEngine,false>(engine,candidate->forcing_); if(!status) return status;
  status=candidate->prepareForcingEvaluationContext(); if(!status) return status;
  result=std::move(candidate); return WVKernelStatus::ok();
}
WVKernelStatus WVFieldEvaluationService::createBorrowing(WVBarotropicQGForcingEngine& engine,std::unique_ptr<WVFieldEvaluationService>& result) {
  std::unique_ptr<WVFieldEvaluationService> candidate;
  auto status=createBorrowing(engine.kernel(),candidate); if(!status) return status;
  status=detail::WVForcingDiagnosticBinding::create<WVBarotropicQGForcingEngine,true>(engine,candidate->forcing_); if(!status) return status;
  status=candidate->prepareForcingEvaluationContext(); if(!status) return status;
  result=std::move(candidate); return WVKernelStatus::ok();
}
WVKernelStatus WVFieldEvaluationService::createBorrowing(WVStratifiedQGForcingEngine& engine,std::unique_ptr<WVFieldEvaluationService>& result) {
  std::unique_ptr<WVFieldEvaluationService> candidate;
  auto status=createBorrowing(engine.kernel(),candidate); if(!status) return status;
  status=detail::WVForcingDiagnosticBinding::create<WVStratifiedQGForcingEngine,true>(engine,candidate->forcing_); if(!status) return status;
  status=candidate->prepareForcingEvaluationContext(); if(!status) return status;
  result=std::move(candidate); return WVKernelStatus::ok();
}
WVKernelStatus WVFieldEvaluationService::createBorrowing(WVHydrostaticForcingEngine& engine,std::unique_ptr<WVFieldEvaluationService>& result) {
  std::unique_ptr<WVFieldEvaluationService> candidate;
  auto status=createBorrowing(engine.kernel(),candidate); if(!status) return status;
  status=detail::WVForcingDiagnosticBinding::create<WVHydrostaticForcingEngine,false>(engine,candidate->forcing_); if(!status) return status;
  status=candidate->prepareForcingEvaluationContext(); if(!status) return status;
  result=std::move(candidate); return WVKernelStatus::ok();
}
WVKernelStatus WVFieldEvaluationService::createBorrowing(WVBoussinesqForcingEngine& engine,std::unique_ptr<WVFieldEvaluationService>& result) {
  std::unique_ptr<WVFieldEvaluationService> candidate;
  auto status=createBorrowing(engine.kernel(),candidate); if(!status) return status;
  status=detail::WVForcingDiagnosticBinding::create<WVBoussinesqForcingEngine,false>(engine,candidate->forcing_); if(!status) return status;
  status=candidate->prepareForcingEvaluationContext(); if(!status) return status;
  result=std::move(candidate); return WVKernelStatus::ok();
}
WVKernelStatus WVFieldEvaluationService::createPlan(
    const std::vector<WVFieldRequest> &requests,
    WVFieldEvaluationPlan &plan, WVDensityDiagnosticContract densityContract) const {
  return createPlanImpl(requests,plan,densityContract,true);
}

WVKernelStatus WVFieldEvaluationService::createPlanForActiveEvaluation(
    const std::vector<WVFieldRequest>& requests,WVFieldEvaluationPlan& plan,
    WVDensityDiagnosticContract densityContract) const {
  if(!eventWorkspace_ || !stateEvaluationActive_)
    return {WVKernelStatusCode::invalidConfiguration,
        "Active field-plan creation requires an evaluation session."};
  if(variableEvaluationPolicy_!=WVVariableEvaluationPolicy::reuse)
    return {WVKernelStatusCode::invalidConfiguration,
        "Active field-plan creation requires the reuse policy."};
  if(executing_ || !eventWorkspace_->preparationIdle())
    return {WVKernelStatusCode::reentrantExecution,
        "Active field-plan creation is allowed only between producer calls."};
  WVFieldEvaluationPlan candidate;
  auto status=createPlanImpl(requests,candidate,densityContract,false);
  if(!status) return status;
  const bool forcingPreparedBefore=eventArenaForcingPrepared_;
  status=prepareEventArenaForActiveEvaluation(candidate);
  eventArena_->refreshAccounting();
  auto& metrics=const_cast<WVFieldEvaluationService*>(this)->mutableMetrics();
  metrics.eventFieldArenaPlannedBytes=eventArena_->plannedBytes;
  metrics.eventFieldArenaPeakBytes=eventArena_->peakBytes;
  metrics.servicePersistentBytes=persistentBytes();
  if(!status) {
    eventArenaForcingPrepared_=forcingPreparedBefore;
    return status;
  }
  plan=std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::createPlanImpl(
    const std::vector<WVFieldRequest> &requests,
    WVFieldEvaluationPlan &plan,WVDensityDiagnosticContract densityContract,
    bool prepareScientificDependencies) const {
  if (densityContract.reference != WVNoMotionReference::actual &&
      densityContract.reference != WVNoMotionReference::initial)
    return {WVKernelStatusCode::invalidConfiguration,"Invalid density diagnostic reference."};
  if (detail::WVDiagnosticFieldPlan::required(requests,stratified_ != nullptr))
  {
    WVFieldEvaluationPlan candidate;
    const auto status=detail::WVDiagnosticFieldPlan::create(
        *this,requests,candidate,densityContract,prepareScientificDependencies);
    if(!status) return status;
    const bool prepareForcing=candidate.diagnosticPlan_ &&
        candidate.diagnosticPlan_->hasForcingDiagnostics();
    const bool forcingPreparedBefore=eventArenaForcingPrepared_;
    if(prepareScientificDependencies) eventArenaForcingPrepared_|=prepareForcing;
    if(prepareScientificDependencies) {
      const auto arenaStatus=prepareEventArena(candidate);
      if(!arenaStatus) {
        eventArenaForcingPrepared_=forcingPreparedBefore;
        return arenaStatus;
      }
    }
    plan=std::move(candidate);
    return WVKernelStatus::ok();
  }
  if (barotropicQG_) {
    WVFieldEvaluationPlan candidate;
    const auto status=barotropicQG_->createPlan(requests, candidate);
    if(!status) return status;
    if(prepareScientificDependencies) {
      const auto arenaStatus=prepareEventArena(candidate);
      if(!arenaStatus) return arenaStatus;
    }
    plan=std::move(candidate);
    return WVKernelStatus::ok();
  }
  if (stratified_) {
    WVFieldEvaluationPlan candidate;
    const auto status=stratified_->createPlan(requests, candidate);
    if(!status) return status;
    if(prepareScientificDependencies) {
      const auto arenaStatus=prepareEventArena(candidate);
      if(!arenaStatus) return arenaStatus;
    }
    plan=std::move(candidate);
    return WVKernelStatus::ok();
  }
  try {
    WVFieldEvaluationPlan candidate;
    const auto &configuration = transform_->descriptor().configuration();
    candidate.configuration_ = configuration;
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    bool needsSpline = false;
    for (const auto &request : requests)
      needsSpline = needsSpline ||
                    (request.sampling.kind == WVFieldSamplingKind::positions &&
                     request.sampling.interpolation ==
                         WVPositionInterpolation::spline);
    std::unique_ptr<SplineSystem> xSpline;
    std::unique_ptr<SplineSystem> ySpline;
    std::unique_ptr<SplineSystem> zSpline;
    if (needsSpline) {
      xSpline = std::make_unique<SplineSystem>(configuration.Nx);
      ySpline = std::make_unique<SplineSystem>(configuration.Ny);
    }
    const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
    const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
    const double dz = configuration.Lz /
                      static_cast<double>(configuration.Nz - 1);
    for (std::size_t outputIndex = 0; outputIndex < requests.size();
         ++outputIndex) {
      const auto &request = requests[outputIndex];
      if (request.identifier.empty())
        return invalid("Every field request must have a nonempty identifier.");
      if (!identifiers.insert(request.identifier).second)
        return invalid("Field request identifiers must be unique: " +
                       request.identifier + ".");
      const auto *metadata = findExecutablePortableVariable(request.fieldName);
      if (metadata == nullptr || metadata->kind != WVPortableVariableKind::field)
        return invalid("Unknown or unsupported field: " + request.fieldName + ".");
      const auto fieldIndex = static_cast<std::size_t>(metadata->ordinal);
      const auto field = metadata->identifier;
      const auto rank = metadata->naturalRank;

      WVFieldEvaluationPlan::ResolvedRequest resolved;
      if (request.sampling.kind != WVFieldSamplingKind::fullGrid &&
          request.sampling.kind !=
              WVFieldSamplingKind::fixedVerticalProfiles &&
          request.sampling.kind != WVFieldSamplingKind::positions)
        return invalid("Unknown field-sampling kind.");
      const std::uint8_t requestedSampling =
          request.sampling.kind == WVFieldSamplingKind::fullGrid
              ? portableFullGridSampling
              : request.sampling.kind ==
                        WVFieldSamplingKind::fixedVerticalProfiles
                    ? portableFixedVerticalProfileSampling
                    : portablePositionSampling;
      if ((metadata->samplingMask & requestedSampling) == 0)
        return invalid("Sampling mode is unsupported for field " +
                       request.fieldName + ".");
      if (request.sampling.kind == WVFieldSamplingKind::positions &&
          request.sampling.interpolation != WVPositionInterpolation::linear &&
          request.sampling.interpolation != WVPositionInterpolation::spline)
        return invalid("Unknown position-interpolation method.");
      resolved.field = field;
      resolved.nativeRank = rank;
      resolved.samplingKind = request.sampling.kind;
      resolved.interpolation = request.sampling.interpolation;
      resolved.outputIndex = outputIndex;
      WVFieldOutputSpecification output;
      output.identifier = request.identifier;
      output.fieldName = request.fieldName;
      output.samplingKind = request.sampling.kind;
      if (request.sampling.kind == WVFieldSamplingKind::fullGrid) {
        if (rank == WVFieldEvaluationPlan::NativeRank::volume)
          output.dimensions = {configuration.Nx, configuration.Ny,
                               configuration.Nz};
        else if (rank == WVFieldEvaluationPlan::NativeRank::horizontal)
          output.dimensions = {configuration.Nx, configuration.Ny};
        else if (rank == WVFieldEvaluationPlan::NativeRank::vertical)
          output.dimensions = {configuration.Nz};
      } else if (request.sampling.kind ==
                 WVFieldSamplingKind::fixedVerticalProfiles) {
        if (rank != WVFieldEvaluationPlan::NativeRank::volume)
          return invalid("Fixed vertical profiles require a three-dimensional field: " +
                         request.fieldName + ".");
        if (request.sampling.xIndices.empty() ||
            request.sampling.xIndices.size() !=
                request.sampling.yIndices.size())
          return invalid("Fixed-profile xIndices and yIndices must be nonempty and have equal length.");
        resolved.profileXIndices.reserve(request.sampling.xIndices.size());
        resolved.profileYIndices.reserve(request.sampling.yIndices.size());
        for (std::size_t index = 0; index < request.sampling.xIndices.size();
             ++index) {
          const auto xIndex = request.sampling.xIndices[index];
          const auto yIndex = request.sampling.yIndices[index];
          if (xIndex == 0 || xIndex > configuration.Nx || yIndex == 0 ||
              yIndex > configuration.Ny)
            return invalid("Fixed-profile indices must use MATLAB one-based values within the model grid.");
          resolved.profileXIndices.push_back(xIndex - 1);
          resolved.profileYIndices.push_back(yIndex - 1);
        }
        output.dimensions = {configuration.Nz,
                             request.sampling.xIndices.size()};
      } else {
        if (rank == WVFieldEvaluationPlan::NativeRank::scalar ||
            rank == WVFieldEvaluationPlan::NativeRank::vertical)
          return invalid("Position sampling is unsupported for field " +
                         request.fieldName + ".");
        const auto positionCount = request.sampling.x.size();
        if (positionCount == 0 || request.sampling.y.size() != positionCount)
          return invalid("Position x and y arrays must be nonempty and have equal length.");
        if (rank == WVFieldEvaluationPlan::NativeRank::volume &&
            request.sampling.z.size() != positionCount)
          return invalid("Three-dimensional position sampling requires equal-length x, y, and z arrays.");
        if (rank == WVFieldEvaluationPlan::NativeRank::horizontal &&
            !request.sampling.z.empty() &&
            request.sampling.z.size() != positionCount)
          return invalid("Optional z positions must be empty or match x and y.");
        resolved.positionWeights.reserve(positionCount);
        for (std::size_t position = 0; position < positionCount; ++position) {
          const double x = request.sampling.x[position];
          const double y = request.sampling.y[position];
          const double z = rank == WVFieldEvaluationPlan::NativeRank::volume
                               ? request.sampling.z[position]
                               : 0.0;
          if (!std::isfinite(x) || !std::isfinite(y) ||
              (rank == WVFieldEvaluationPlan::NativeRank::volume &&
               !std::isfinite(z)))
            return invalid("Position coordinates must be finite.");
          WVFieldEvaluationPlan::PositionWeights weights;
          const double xWrapped = wrapped(x, configuration.Lx);
          const double yWrapped = wrapped(y, configuration.Ly);
          const auto xLower = std::min(
              static_cast<std::size_t>(std::floor(xWrapped / dx)),
              configuration.Nx - 1);
          const auto yLower = std::min(
              static_cast<std::size_t>(std::floor(yWrapped / dy)),
              configuration.Ny - 1);
          if (request.sampling.interpolation ==
              WVPositionInterpolation::linear) {
            weights.xLinearIndices = {xLower, (xLower + 1) % configuration.Nx};
            weights.yLinearIndices = {yLower, (yLower + 1) % configuration.Ny};
            const double xFraction =
                (xWrapped - static_cast<double>(xLower) * dx) / dx;
            const double yFraction =
                (yWrapped - static_cast<double>(yLower) * dy) / dy;
            weights.xLinearWeights = {1.0 - xFraction, xFraction};
            weights.yLinearWeights = {1.0 - yFraction, yFraction};
            if (rank == WVFieldEvaluationPlan::NativeRank::volume) {
              weights.outsideInterpolationDomain =
                  z < -configuration.Lz || z > 0.0;
              if (!weights.outsideInterpolationDomain) {
                const double normalizedZ = (z + configuration.Lz) / dz;
                const auto zLower = std::min(
                    static_cast<std::size_t>(std::max(0.0, std::floor(normalizedZ))),
                    configuration.Nz - 2);
                const double zFraction =
                    std::clamp(normalizedZ - static_cast<double>(zLower), 0.0,
                               1.0);
                weights.zLinearIndices = {zLower, zLower + 1};
                weights.zLinearWeights = {1.0 - zFraction, zFraction};
              }
            }
          } else {
            const bool xBoundary = xLower < 3 || xLower > configuration.Nx - 4;
            const bool yBoundary = yLower < 3 || yLower > configuration.Ny - 4;
            const std::size_t xShift = xBoundary ? 4 : 0;
            const std::size_t yShift = yBoundary ? 4 : 0;
            const double xQuery = xBoundary
                                      ? wrapped(x + 4.0 * dx, configuration.Lx)
                                      : xWrapped;
            const double yQuery = yBoundary
                                      ? wrapped(y + 4.0 * dy, configuration.Ly)
                                      : yWrapped;
            // MATLAB's interpn(...,"spline",0) applies the zero fill value
            // when a boundary-shifted query lies beyond the last stored
            // periodic grid point (the duplicated x=L/y=L points are absent).
            weights.outsideInterpolationDomain =
                xQuery > static_cast<double>(configuration.Nx - 1) * dx ||
                yQuery > static_cast<double>(configuration.Ny - 1) * dy;
            if (!weights.outsideInterpolationDomain) {
              weights.xSplineWeights =
                  xSpline->weights(0.0, dx, xQuery, xShift);
              weights.ySplineWeights =
                  ySpline->weights(0.0, dy, yQuery, yShift);
            }
            if (rank == WVFieldEvaluationPlan::NativeRank::volume) {
              weights.outsideInterpolationDomain =
                  weights.outsideInterpolationDomain ||
                  z < -configuration.Lz || z > 0.0;
              if (!weights.outsideInterpolationDomain) {
                if (!zSpline)
                  zSpline = std::make_unique<SplineSystem>(configuration.Nz);
                weights.zSplineWeights =
                    zSpline->weights(-configuration.Lz, dz, z);
              }
            }
          }
          resolved.positionWeights.push_back(std::move(weights));
        }
        output.dimensions = {positionCount};
      }
      output.elementCount = 1;
      for (const auto dimension : output.dimensions)
        output.elementCount = checkedProduct(output.elementCount, dimension);
      candidate.requestedFieldMask_ |= 1ULL << fieldIndex;
      if (field == WVFieldEvaluationPlan::Field::psi &&
          transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Streamfunction evaluation is undefined when the Coriolis "
                "frequency is zero."};
      resolved.dependencyMask = metadata->primitiveDependencyMask;
      candidate.dependencyMask_ |= metadata->primitiveDependencyMask;
      candidate.requests_.push_back(std::move(resolved));
      candidate.outputs_.push_back(std::move(output));
    }
    if(prepareScientificDependencies) {
      const auto arenaStatus=prepareEventArena(candidate);
      if(!arenaStatus) return arenaStatus;
    }
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate immutable field-evaluation plan storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::invalid_argument &error) {
    return invalid(error.what());
  }
}

WVKernelStatus WVFieldEvaluationService::evaluate(
    const WVFieldEvaluationPlan &plan, const WVState &state,
    WVFieldOutputView *outputs, std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const bool anyActive=outputCount && (!activeOutputs ||
      std::any_of(activeOutputs,activeOutputs+outputCount,
          [](std::uint8_t selected) {return selected!=0;}));
  if(!eventWorkspace_ && anyActive) {
    if (!transform_)
      return {WVKernelStatusCode::unsupportedOperation,
          "This transform requires coefficient-family state views."};
    detail::WVFieldEvaluationEventScope scope(*this,{state});
    if(!scope.status()) return scope.status();
    return evaluate(plan,state,outputs,outputCount,activeOutputs);
  }
  if(eventWorkspace_) {
    const auto status=eventWorkspace_->validateState({state});
    if(!status) return status;
  }
  if (plan.diagnosticPlan_) {
    auto& metrics=stratified_ ? stratified_->mutableMetrics() :
        barotropicQG_ ? barotropicQG_->mutableMetrics() : metrics_;
    const auto evaluationsBefore=metrics.evaluationCount;
    const auto batchesBefore=metrics.coincidentBatchCount;
    const auto status=plan.diagnosticPlan_->evaluate(
        *this,{state},outputs,outputCount,activeOutputs);
    if(status) {
      std::size_t activeCount=0;
      for(std::size_t i=0;i<outputCount;++i)
        activeCount+=!activeOutputs || activeOutputs[i];
      metrics.evaluationCount=evaluationsBefore+1;
      metrics.coincidentBatchCount=batchesBefore+(activeCount>1 ? 1 : 0);
      metrics.lastPlanBytes=plan.persistentBytes();
      metrics.maximumPlanBytes=std::max(
          metrics.maximumPlanBytes,metrics.lastPlanBytes);
    }
    return status;
  }
  if (!transform_) return {WVKernelStatusCode::unsupportedOperation,"This transform requires coefficient-family state views."};
  const PlanInvocation invocation{&plan, outputs, outputCount, activeOutputs};
  return evaluatePlanBatch(&invocation, 1, state);
}

WVKernelStatus WVFieldEvaluationService::evaluate(
    const WVFieldEvaluationPlan &plan, const WVIntegrationState &state,
    WVFieldOutputView *outputs, std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const bool anyActive=outputCount && (!activeOutputs ||
      std::any_of(activeOutputs,activeOutputs+outputCount,
          [](std::uint8_t selected) {return selected!=0;}));
  if(!eventWorkspace_ && anyActive) {
    detail::WVFieldEvaluationEventScope scope(*this,state);
    if(!scope.status()) return scope.status();
    return evaluate(plan,state,outputs,outputCount,activeOutputs);
  }
  if(eventWorkspace_) {
    const auto status=eventWorkspace_->validateState({state});
    if(!status) return status;
  }
  if (plan.diagnosticPlan_) {
    auto& metrics=stratified_ ? stratified_->mutableMetrics() :
        barotropicQG_ ? barotropicQG_->mutableMetrics() : metrics_;
    const auto evaluationsBefore=metrics.evaluationCount;
    const auto batchesBefore=metrics.coincidentBatchCount;
    const auto status=plan.diagnosticPlan_->evaluate(
        *this,state,outputs,outputCount,activeOutputs);
    if(status) {
      std::size_t activeCount=0;
      for(std::size_t i=0;i<outputCount;++i)
        activeCount+=!activeOutputs || activeOutputs[i];
      metrics.evaluationCount=evaluationsBefore+1;
      metrics.coincidentBatchCount=batchesBefore+(activeCount>1 ? 1 : 0);
      metrics.lastPlanBytes=plan.persistentBytes();
      metrics.maximumPlanBytes=std::max(
          metrics.maximumPlanBytes,metrics.lastPlanBytes);
    }
    return status;
  }
  if (barotropicQG_)
    return barotropicQG_->evaluate(plan, state, outputs, outputCount, activeOutputs);
  if (stratified_)
    return stratified_->evaluate(plan, state, outputs, outputCount, activeOutputs);
  return evaluate(plan, state.waveVortex, outputs, outputCount, activeOutputs);
}

WVKernelStatus
WVFieldEvaluationService::evaluatePlanBatch(const PlanInvocation *invocations,
                                            std::size_t invocationCount,
                                            const WVState &state) {
  if (invocationCount == 0)
    return WVKernelStatus::ok();
  if (invocations == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Field-evaluation plan invocations have a null pointer."};

  std::size_t totalOutputCount = 0;
  std::size_t planBytes = 0;
  std::uint64_t requestedFieldMask = 0;
  std::uint64_t dependencyMask = 0;
  bool allOutputsEmpty = true;
  for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
       ++invocationIndex) {
    const auto &invocation = invocations[invocationIndex];
    if (invocation.plan == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation invocation has a null plan."};
    const auto &plan = *invocation.plan;
    if (!sameTransformConfiguration(plan.configuration_,
                                    transform_->descriptor().configuration()))
      return invalid("The evaluation plan was created for a different "
                     "transform configuration.");
    if (invocation.outputCount != plan.outputs_.size())
      return {WVKernelStatusCode::invalidShape,
              "The output-view count must match the evaluation plan."};
    if (invocation.outputCount != 0 && invocation.outputs == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation outputs have a null view pointer."};
    if (totalOutputCount >
        std::numeric_limits<std::size_t>::max() - invocation.outputCount)
      return {WVKernelStatusCode::sizeOverflow,
              "Field-evaluation batch output count overflows size_t."};
    totalOutputCount += invocation.outputCount;
    const auto invocationPlanBytes = plan.persistentBytes();
    planBytes = planBytes > std::numeric_limits<std::size_t>::max() -
                                invocationPlanBytes
                    ? std::numeric_limits<std::size_t>::max()
                    : planBytes + invocationPlanBytes;
    if (!invocation.activeOutputs) {
      requestedFieldMask |= plan.requestedFieldMask_;
      dependencyMask |= plan.dependencyMask_;
    } else for (const auto &request : plan.requests_) {
      if (!invocation.activeOutputs[request.outputIndex]) continue;
      requestedFieldMask |= 1ULL << static_cast<std::size_t>(request.field);
      dependencyMask |= request.dependencyMask;
    }
    for (std::size_t output = 0; output < plan.outputs_.size(); ++output)
      if (!invocation.activeOutputs || invocation.activeOutputs[output])
        allOutputsEmpty = allOutputsEmpty && plan.outputs_[output].elementCount == 0;
  }
  if (!std::isfinite(state.t) || !std::isfinite(state.t0))
    return invalid("Field-evaluation state times must be finite.");
  const auto spectral = transform_->descriptor().spectralShape();
  const auto coefficientBytes = spectral.elementCount() * sizeof(WVComplex64);
  const void *coefficientInputs[] = {state.coefficients.Ap.data,
                                     state.coefficients.Am.data,
                                     state.coefficients.A0.data};
  const WVComplexConstView coefficientViews[] = {
      state.coefficients.Ap, state.coefficients.Am, state.coefficients.A0};
  for (const auto &view : coefficientViews) {
    if (view.shape.rows != spectral.rows ||
        view.shape.columns != spectral.columns)
      return {WVKernelStatusCode::invalidShape,
              "Field-evaluation coefficients must have shape [Nj,Nkl]."};
    if (view.data == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation coefficients have a null pointer."};
  }
  for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
       ++invocationIndex) {
    const auto &invocation = invocations[invocationIndex];
    const auto &plan = *invocation.plan;
    for (std::size_t outputIndex = 0; outputIndex < invocation.outputCount;
         ++outputIndex) {
      if (invocation.activeOutputs && !invocation.activeOutputs[outputIndex]) continue;
      const auto &output = invocation.outputs[outputIndex];
      if (output.elementCount != plan.outputs_[outputIndex].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Caller-owned output has the wrong element count for request " +
                    plan.outputs_[outputIndex].identifier + "."};
      if (output.elementCount != 0 && output.data == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Caller-owned output has a null pointer for request " +
                    plan.outputs_[outputIndex].identifier + "."};
      const auto bytes = output.elementCount * sizeof(double);
      for (const auto *input : coefficientInputs)
        if (memoryOverlaps(output.data, bytes, input, coefficientBytes))
          return {WVKernelStatusCode::overlappingArrays,
                  "Field outputs must not overlap coefficient inputs."};
      for (std::size_t otherInvocation = invocationIndex;
           otherInvocation < invocationCount; ++otherInvocation) {
        const auto &other = invocations[otherInvocation];
        const std::size_t firstOther =
            otherInvocation == invocationIndex ? outputIndex + 1 : 0;
        for (std::size_t otherOutput = firstOther;
             otherOutput < other.outputCount; ++otherOutput)
          if ((!other.activeOutputs || other.activeOutputs[otherOutput]) && memoryOverlaps(
                  output.data, bytes, other.outputs[otherOutput].data,
                  other.outputs[otherOutput].elementCount * sizeof(double)))
            return {WVKernelStatusCode::overlappingArrays,
                    "Caller-owned field outputs must not overlap each other."};
      }
    }
  }
  ExecutionGuard guard(executing_);
  if (!guard.entered())
    return {WVKernelStatusCode::reentrantExecution,
            "Field evaluation is not reentrant."};

  ++metrics_.evaluationCount;
  if (totalOutputCount > 1)
    ++metrics_.coincidentBatchCount;
  metrics_.lastPlanBytes = planBytes;
  metrics_.maximumPlanBytes =
      std::max(metrics_.maximumPlanBytes, metrics_.lastPlanBytes);

  if (allOutputsEmpty)
    return WVKernelStatus::ok();

  const auto &configuration = transform_->descriptor().configuration();
  const auto spatial = transform_->descriptor().spatialShape();
  const auto fieldElements = spatial.elementCount();
  const auto horizontalElements = configuration.Nx * configuration.Ny;
  const auto coefficientElements = spectral.elementCount();
  auto updateScratchHighWater = [&](std::size_t realElements,
                                    std::size_t complexElements) {
    const auto bytes =
        realElements * sizeof(double) + complexElements * sizeof(WVComplex64);
    metrics_.scratchHighWaterBytes =
        std::max(metrics_.scratchHighWaterBytes, bytes);
  };
  auto invokeTransform = [&](auto &&operation) {
    const auto before = transform_->metrics().executionCount;
    ++metrics_.transformCount;
    const auto status = operation();
    metrics_.fftExecutionCount += transform_->metrics().executionCount - before;
    return status;
  };

  auto fieldRequested = [&](WVFieldEvaluationPlan::Field field) {
    return (requestedFieldMask & (1ULL << static_cast<std::size_t>(field))) !=
           0;
  };

  auto writeField = [&](WVFieldEvaluationPlan::Field field,
                        WVFieldEvaluationPlan::NativeRank rank,
                        const double *source) -> WVKernelStatus {
    std::size_t sourceElements = 1;
    if (rank == WVFieldEvaluationPlan::NativeRank::volume)
      sourceElements = fieldElements;
    else if (rank == WVFieldEvaluationPlan::NativeRank::horizontal)
      sourceElements = horizontalElements;
    else if (rank == WVFieldEvaluationPlan::NativeRank::vertical)
      sourceElements = configuration.Nz;
    const double *samplingSource = source;
    const WVFieldEvaluationPlan::ResolvedRequest *firstFullRequest = nullptr;
    std::size_t firstFullInvocation = invocationCount;
    std::size_t consumer = 0;
    for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
         ++invocationIndex) {
      const auto &invocation = invocations[invocationIndex];
      for (const auto &request : invocation.plan->requests_) {
        if (invocation.activeOutputs && !invocation.activeOutputs[request.outputIndex]) continue;
        if (request.field != field)
          continue;
        if (request.samplingKind == WVFieldSamplingKind::fullGrid &&
            firstFullRequest == nullptr) {
          std::copy(source, source + sourceElements,
                    invocation.outputs[request.outputIndex].data);
          samplingSource = invocation.outputs[request.outputIndex].data;
          firstFullRequest = &request;
          firstFullInvocation = invocationIndex;
          ++metrics_.fullGridWriteCount;
          metrics_.outputElementWriteCount += sourceElements;
        }
        ++consumer;
      }
    }
    if (consumer > 1)
      metrics_.primitiveFieldReuseCount += consumer - 1;
    for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
         ++invocationIndex) {
      const auto &invocation = invocations[invocationIndex];
      for (const auto &request : invocation.plan->requests_) {
        if (invocation.activeOutputs && !invocation.activeOutputs[request.outputIndex]) continue;
        if (request.field != field)
          continue;
        auto &output = invocation.outputs[request.outputIndex];
        if (request.samplingKind == WVFieldSamplingKind::fullGrid) {
          if (invocationIndex != firstFullInvocation ||
              &request != firstFullRequest) {
            std::copy(samplingSource, samplingSource + sourceElements,
                      output.data);
            ++metrics_.fullGridWriteCount;
            metrics_.outputElementWriteCount += sourceElements;
          }
          continue;
        }
        if (request.samplingKind ==
            WVFieldSamplingKind::fixedVerticalProfiles) {
          for (std::size_t profile = 0;
               profile < request.profileXIndices.size(); ++profile) {
            const auto horizontalIndex =
                request.profileXIndices[profile] +
                configuration.Nx * request.profileYIndices[profile];
            for (std::size_t z = 0; z < configuration.Nz; ++z)
              output.data[z + configuration.Nz * profile] =
                  samplingSource[horizontalIndex + horizontalElements * z];
          }
          ++metrics_.profileWriteCount;
          metrics_.outputElementWriteCount += output.elementCount;
          continue;
        }
        for (std::size_t position = 0;
             position < request.positionWeights.size(); ++position) {
          const auto &weights = request.positionWeights[position];
          double value = 0.0;
          if (!weights.outsideInterpolationDomain) {
            if (request.interpolation == WVPositionInterpolation::linear) {
              const std::size_t zCount =
                  rank == WVFieldEvaluationPlan::NativeRank::volume ? 2 : 1;
              for (std::size_t iz = 0; iz < zCount; ++iz)
                for (std::size_t iy = 0; iy < 2; ++iy)
                  for (std::size_t ix = 0; ix < 2; ++ix) {
                    const auto zIndex =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zLinearIndices[iz]
                            : 0;
                    const auto sourceIndex =
                        weights.xLinearIndices[ix] +
                        configuration.Nx * weights.yLinearIndices[iy] +
                        horizontalElements * zIndex;
                    const double zWeight =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zLinearWeights[iz]
                            : 1.0;
                    value += samplingSource[sourceIndex] *
                             weights.xLinearWeights[ix] *
                             weights.yLinearWeights[iy] * zWeight;
                  }
              ++metrics_.linearInterpolationCount;
            } else {
              const std::size_t zCount =
                  rank == WVFieldEvaluationPlan::NativeRank::volume
                      ? configuration.Nz
                      : 1;
              for (std::size_t iz = 0; iz < zCount; ++iz)
                for (std::size_t iy = 0; iy < configuration.Ny; ++iy)
                  for (std::size_t ix = 0; ix < configuration.Nx; ++ix) {
                    const auto sourceIndex =
                        ix + configuration.Nx * iy + horizontalElements * iz;
                    const double zWeight =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zSplineWeights[iz]
                            : 1.0;
                    value += samplingSource[sourceIndex] *
                             weights.xSplineWeights[ix] *
                             weights.ySplineWeights[iy] * zWeight;
                  }
              ++metrics_.splineInterpolationCount;
            }
          } else if (request.interpolation == WVPositionInterpolation::linear) {
            ++metrics_.linearInterpolationCount;
          } else {
            ++metrics_.splineInterpolationCount;
          }
          output.data[position] = value;
        }
        metrics_.outputElementWriteCount += output.elementCount;
      }
    }
    return WVKernelStatus::ok();
  };

  if ((dependencyMask & primitiveValues) != 0) {
    updateScratchHighWater(4 * fieldElements, 0);
    WVRealFieldBundleView primitiveBundle{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    bool reused=false;
    const auto operation=[&]() {return invokeTransform([&]() { return transform_->transformWaveVortexToUVWEta(state, primitiveBundle); });};
    const WVVariableEvaluationKey primitiveKey{
        WVVariableEvaluationNode::physicalField,
        static_cast<std::uint32_t>(WVPortableVariable::u),
        eventWorkspace_ ? eventWorkspace_->component() : 0u,0,0,0,1};
    auto status = eventWorkspace_ ? eventWorkspace_->evaluate(primitiveKey,primitiveBundle.data,4*fieldElements,operation,reused) : operation();
    if(reused) ++metrics_.primitiveFieldReuseCount;
    if (!status)
      return status;
    const bool primitiveNeeded[] = {
        fieldRequested(WVFieldEvaluationPlan::Field::u) ||
            fieldRequested(WVFieldEvaluationPlan::Field::ssu) ||
            fieldRequested(WVFieldEvaluationPlan::Field::uvMax),
        fieldRequested(WVFieldEvaluationPlan::Field::v) ||
            fieldRequested(WVFieldEvaluationPlan::Field::ssv) ||
            fieldRequested(WVFieldEvaluationPlan::Field::uvMax),
        fieldRequested(WVFieldEvaluationPlan::Field::w) ||
            fieldRequested(WVFieldEvaluationPlan::Field::wMax),
        fieldRequested(WVFieldEvaluationPlan::Field::eta) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoE) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoBar)};
    if(!reused) metrics_.primitiveFieldEvaluationCount +=
        static_cast<std::size_t>(primitiveNeeded[0]) +
        static_cast<std::size_t>(primitiveNeeded[1]) +
        static_cast<std::size_t>(primitiveNeeded[2]) +
        static_cast<std::size_t>(primitiveNeeded[3]);
    const double *u = realScratch_.data();
    const double *v = u + fieldElements;
    const double *w = v + fieldElements;
    const double *eta = w + fieldElements;
    if (fieldRequested(WVFieldEvaluationPlan::Field::u))
      writeField(WVFieldEvaluationPlan::Field::u,
                 WVFieldEvaluationPlan::NativeRank::volume, u);
    if (fieldRequested(WVFieldEvaluationPlan::Field::v))
      writeField(WVFieldEvaluationPlan::Field::v,
                 WVFieldEvaluationPlan::NativeRank::volume, v);
    if (fieldRequested(WVFieldEvaluationPlan::Field::w))
      writeField(WVFieldEvaluationPlan::Field::w,
                 WVFieldEvaluationPlan::NativeRank::volume, w);
    if (fieldRequested(WVFieldEvaluationPlan::Field::eta))
      writeField(WVFieldEvaluationPlan::Field::eta,
                 WVFieldEvaluationPlan::NativeRank::volume, eta);
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssu))
      writeField(WVFieldEvaluationPlan::Field::ssu,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 u + fieldElements - horizontalElements);
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssv))
      writeField(WVFieldEvaluationPlan::Field::ssv,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 v + fieldElements - horizontalElements);
    if (fieldRequested(WVFieldEvaluationPlan::Field::uvMax)) {
      double maximum = 0.0;
      for (std::size_t index = 0; index < fieldElements; ++index)
        maximum = std::max(maximum, std::sqrt(u[index] * u[index] +
                                              v[index] * v[index]));
      writeField(WVFieldEvaluationPlan::Field::uvMax,
                 WVFieldEvaluationPlan::NativeRank::scalar, &maximum);
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::wMax)) {
      double maximum = 0.0;
      for (std::size_t index = 0; index < fieldElements; ++index)
        maximum = std::max(maximum, std::abs(w[index]));
      writeField(WVFieldEvaluationPlan::Field::wMax,
                 WVFieldEvaluationPlan::NativeRank::scalar, &maximum);
    }
    double *derived = realScratch_.data() + 4 * fieldElements;
    const double densityScale =
        configuration.rho0 * configuration.N0 * configuration.N0 /
        configuration.g;
    if (fieldRequested(WVFieldEvaluationPlan::Field::rhoE) ||
        fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal)) {
      updateScratchHighWater(5 * fieldElements, 0);
      if (fieldRequested(WVFieldEvaluationPlan::Field::rhoE)) {
        const auto produce=[&]() {
          for (std::size_t index = 0; index < fieldElements; ++index)
            derived[index] = densityScale * eta[index];
          ++outputProducerMetrics_.reconstructions[
              static_cast<std::size_t>(WVHydrostaticField::rhoE)][0]
              [eventWorkspace_ ? eventWorkspace_->component() : 0u];
          return WVKernelStatus::ok();
        };
        bool reusedDensity=false;
        const auto densityStatus=eventWorkspace_ ? eventWorkspace_->evaluate(
            {WVVariableEvaluationNode::physicalField,
              static_cast<std::uint32_t>(WVPortableVariable::rhoE),
              eventWorkspace_->component()},derived,fieldElements,produce,
            reusedDensity) : produce();
        if(!densityStatus) return densityStatus;
        writeField(WVFieldEvaluationPlan::Field::rhoE,
                   WVFieldEvaluationPlan::NativeRank::volume, derived);
      }
      if (fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal)) {
        const auto produce=[&]() {
          for (std::size_t z = 0; z < configuration.Nz; ++z) {
            const double rhoNoMotion=
                configuration.rho0-densityScale*
                    transform_->descriptor().verticalModes().z[z];
            for (std::size_t horizontal=0;
                horizontal<horizontalElements;++horizontal) {
              const auto index=horizontal+horizontalElements*z;
              derived[index]=rhoNoMotion+densityScale*eta[index];
            }
          }
          ++outputProducerMetrics_.reconstructions[
              static_cast<std::size_t>(WVHydrostaticField::rhoTotal)][0]
              [eventWorkspace_ ? eventWorkspace_->component() : 0u];
          return WVKernelStatus::ok();
        };
        bool reusedDensity=false;
        const auto densityStatus=eventWorkspace_ ? eventWorkspace_->evaluate(
            {WVVariableEvaluationNode::physicalField,
              static_cast<std::uint32_t>(WVPortableVariable::rhoTotal),
              eventWorkspace_->component()},derived,fieldElements,produce,
            reusedDensity) : produce();
        if(!densityStatus) return densityStatus;
        writeField(WVFieldEvaluationPlan::Field::rhoTotal,
                   WVFieldEvaluationPlan::NativeRank::volume, derived);
      }
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::rhoBar)) {
      updateScratchHighWater(4 * fieldElements + configuration.Nz, 0);
      for (std::size_t z = 0; z < configuration.Nz; ++z) {
        double meanEta = 0.0;
        for (std::size_t horizontal = 0; horizontal < horizontalElements;
             ++horizontal)
          meanEta += eta[horizontal + horizontalElements * z];
        meanEta /= static_cast<double>(horizontalElements);
        derived[z] =
            configuration.rho0 -
            densityScale * transform_->descriptor().verticalModes().z[z] +
            densityScale * meanEta;
      }
      writeField(WVFieldEvaluationPlan::Field::rhoBar,
                 WVFieldEvaluationPlan::NativeRank::vertical, derived);
    }
  }

  auto fillFFieldCoefficients = [&](WVFieldEvaluationPlan::Field field) {
    const auto &modes = transform_->descriptor().verticalModes();
    const auto &horizontalModes = transform_->descriptor().fourierModes();
    auto *wave = complexScratch_.data();
    auto *zeroFrequency = wave + coefficientElements;
    WVComplexConstView phases;
    if(field==WVFieldEvaluationPlan::Field::pi) {
      const auto phaseStatus=transform_->preparedPhase(state,phases);
      if(!phaseStatus) return phaseStatus;
    }
    const double f = modes.coriolisFrequency;
    for (std::size_t mode = 0; mode < horizontalModes.size(); ++mode) {
      const double Kh2 = horizontalModes[mode].Kh * horizontalModes[mode].Kh;
      for (std::size_t j = 0; j < configuration.Nj; ++j) {
        const auto index = j + configuration.Nj * mode;
        wave[index] = {};
        zeroFrequency[index] = {};
        if (field == WVFieldEvaluationPlan::Field::pi) {
          const auto phase=phases.data[index];
          const auto positive = multiply(state.coefficients.Ap.data[index], phase);
          const auto negative =
              multiply(state.coefficients.Am.data[index], conjugate(phase));
          const double rawNAp =
              modes.NApField[index] / modes.gWaveScale[j];
          wave[index] = multiply(subtract(positive, negative), rawNAp);
        }
        if (horizontalModes[mode].Kh > 0.0) {
          const double inverseRossbySquared =
              j == 0 ? 0.0 : f * f / (configuration.g * modes.h0[j]);
          const double denominator = Kh2 + inverseRossbySquared;
          if (field == WVFieldEvaluationPlan::Field::pi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index],
                         -(f / configuration.g) / denominator);
          else if (field == WVFieldEvaluationPlan::Field::psi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], -1.0 / denominator);
          else if (field == WVFieldEvaluationPlan::Field::qgpv)
            zeroFrequency[index] = state.coefficients.A0.data[index];
        } else if (j > 0) {
          if (field == WVFieldEvaluationPlan::Field::pi)
            zeroFrequency[index] = state.coefficients.A0.data[index];
          else if (field == WVFieldEvaluationPlan::Field::psi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], configuration.g / f);
          else if (field == WVFieldEvaluationPlan::Field::qgpv)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], -f / modes.h0[j]);
        }
      }
    }
    return WVKernelStatus::ok();
  };

  bool fFieldReused=false;
  auto evaluateFField = [&](WVFieldEvaluationPlan::Field field) {
    fFieldReused=false;
    updateScratchHighWater(4 * fieldElements, 2 * coefficientElements);
    WVRealFieldBundleView fieldAndDerivatives{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    const WVComplexConstView wave{complexScratch_.data(), spectral};
    const WVComplexConstView zeroFrequency{
        complexScratch_.data() + coefficientElements, spectral};
    const auto identifiedField=field==WVFieldEvaluationPlan::Field::pi ?
        WVConstantFField::pi : field==WVFieldEvaluationPlan::Field::psi ?
        WVConstantFField::psi : WVConstantFField::qgpv;
    const auto operation=[&]() {
      const auto coefficientStatus=fillFFieldCoefficients(field);
      if(!coefficientStatus) return coefficientStatus;
      return invokeTransform([&]() {
        return transform_->transformToSpatialDomainWithFAllDerivatives(
            identifiedField,wave,zeroFrequency,fieldAndDerivatives,
            eventWorkspace_ ? eventWorkspace_->component() : 0);
      });
    };
    const auto status=eventWorkspace_ ? eventWorkspace_->evaluate(
        static_cast<std::size_t>(field),state,realScratch_.data(),
        4*fieldElements,operation,fFieldReused) : operation();
    if(fFieldReused) ++metrics_.primitiveFieldReuseCount;
    return status;
  };

  if ((dependencyMask & pressureHeight) != 0) {
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::pi);
    if (!status)
      return status;
    if(!fFieldReused) ++metrics_.primitiveFieldEvaluationCount;
    const double *pressureHeightField = realScratch_.data();
    if (fieldRequested(WVFieldEvaluationPlan::Field::pi))
      writeField(WVFieldEvaluationPlan::Field::pi,
                 WVFieldEvaluationPlan::NativeRank::volume,
                 pressureHeightField);
    if (fieldRequested(WVFieldEvaluationPlan::Field::p)) {
      double *pressure = realScratch_.data() + 4 * fieldElements;
      updateScratchHighWater(5 * fieldElements, 2 * coefficientElements);
      const double pressureScale = configuration.rho0 * configuration.g;
      const auto producePressure=[&]() {
        for (std::size_t index = 0; index < fieldElements; ++index)
          pressure[index] = pressureScale * pressureHeightField[index];
        ++outputProducerMetrics_.reconstructions[
            static_cast<std::size_t>(WVHydrostaticField::p)][0]
            [eventWorkspace_ ? eventWorkspace_->component() : 0u];
        return WVKernelStatus::ok();
      };
      bool pressureReused=false;
      const auto pressureStatus=eventWorkspace_ ? eventWorkspace_->evaluate(
          {WVVariableEvaluationNode::physicalField,
            static_cast<std::uint32_t>(WVPortableVariable::p),
            eventWorkspace_->component()},pressure,fieldElements,
          producePressure,pressureReused) : producePressure();
      if(!pressureStatus) return pressureStatus;
      writeField(WVFieldEvaluationPlan::Field::p,
                 WVFieldEvaluationPlan::NativeRank::volume, pressure);
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssh))
      writeField(WVFieldEvaluationPlan::Field::ssh,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 pressureHeightField + fieldElements - horizontalElements);
  }

  if ((dependencyMask & streamfunction) != 0) {
    if (transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
      return {WVKernelStatusCode::unsupportedOperation,
              "Streamfunction evaluation is undefined when the Coriolis frequency is zero."};
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::psi);
    if (!status)
      return status;
    if(!fFieldReused) ++metrics_.primitiveFieldEvaluationCount;
    writeField(WVFieldEvaluationPlan::Field::psi,
               WVFieldEvaluationPlan::NativeRank::volume,
               realScratch_.data());
  }

  if ((dependencyMask & potentialVorticity) != 0) {
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::qgpv);
    if (!status)
      return status;
    if(!fFieldReused) ++metrics_.primitiveFieldEvaluationCount;
    writeField(WVFieldEvaluationPlan::Field::qgpv,
               WVFieldEvaluationPlan::NativeRank::volume,
               realScratch_.data());
  }

  if ((dependencyMask & spectralEnergy) != 0) {
    const auto &modes = transform_->descriptor().verticalModes();
    const auto &horizontalModes = transform_->descriptor().fourierModes();
    const double f = modes.coriolisFrequency;
    double energy = 0.0;
    const auto produceEnergy=[&]() {
      energy=0.0;
      for (std::size_t mode = 0; mode < horizontalModes.size(); ++mode) {
      const double Kh = horizontalModes[mode].Kh;
      const double Kh2 = Kh * Kh;
      for (std::size_t j = 0; j < configuration.Nj; ++j) {
        const auto index = j + configuration.Nj * mode;
        const double M = modes.verticalWavenumber[j];
        double hWave = 1.0;
        if (j > 0)
          hWave = configuration.isHydrostatic
                      ? configuration.N0 * configuration.N0 /
                            (configuration.g * M * M)
                      : (configuration.N0 * configuration.N0 - f * f) /
                            (configuration.g * (M * M + Kh2));
        double waveFactor = 0.0;
        if (Kh > 0.0 && j > 0)
          waveFactor = 2.0 * hWave;
        else if (Kh == 0.0)
          waveFactor = j == 0 ? configuration.Lz : hWave;
        double zeroFactor = 0.0;
        if (Kh > 0.0) {
          if (j == 0)
            zeroFactor = configuration.Lz / Kh2;
          else {
            const double denominator =
                Kh2 + f * f / (configuration.g * modes.h0[j]);
            zeroFactor = modes.h0[j] / denominator;
          }
        } else if (j > 0) {
          zeroFactor = configuration.g / 2.0;
        }
        energy += waveFactor *
                      (squaredMagnitude(state.coefficients.Ap.data[index]) +
                       squaredMagnitude(state.coefficients.Am.data[index])) +
                  zeroFactor *
                      squaredMagnitude(state.coefficients.A0.data[index]);
      }
      }
      ++outputProducerMetrics_.energyReductions;
      return WVKernelStatus::ok();
    };
    bool energyReused=false;
    const auto energyStatus=eventWorkspace_ ? eventWorkspace_->evaluate(
        {WVVariableEvaluationNode::reduction,
          static_cast<std::uint32_t>(WVPortableVariable::energy),
          eventWorkspace_->component()},
        &energy,1,produceEnergy,energyReused) : produceEnergy();
    if(!energyStatus) return energyStatus;
    writeField(WVFieldEvaluationPlan::Field::energy,
               WVFieldEvaluationPlan::NativeRank::scalar, &energy);
  }

  const auto derivativeDependencies =
      uDerivatives | vDerivatives | wDerivatives;
  if ((dependencyMask & derivativeDependencies) != 0) {
    updateScratchHighWater(6 * fieldElements, 0);
    double *zetaX = realScratch_.data();
    double *zetaY = zetaX + fieldElements;
    double *zetaZ = zetaY + fieldElements;
    double *derivatives = zetaZ + fieldElements;
    const bool requestedZeta[]={
        fieldRequested(WVFieldEvaluationPlan::Field::zetaX),
        fieldRequested(WVFieldEvaluationPlan::Field::zetaY),
        fieldRequested(WVFieldEvaluationPlan::Field::zetaZ)};
    const WVPortableVariable zetaVariables[]={WVPortableVariable::zetaX,
        WVPortableVariable::zetaY,WVPortableVariable::zetaZ};
    double* zetaFields[]={zetaX,zetaY,zetaZ};
    bool cachedZeta[3]{};
    for(std::size_t index=0;index<3;++index) {
      if(!requestedZeta[index]) continue;
      const WVVariableEvaluationKey key{WVVariableEvaluationNode::physicalField,
          static_cast<std::uint32_t>(zetaVariables[index]),
          eventWorkspace_ ? eventWorkspace_->component() : 0u};
      if(eventWorkspace_ && eventWorkspace_->ready(key)) {
        bool reused=false;
        const auto status=eventWorkspace_->evaluate(key,zetaFields[index],
            fieldElements,[]() {
              return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                  "A ready vorticity field unexpectedly requested production."};
            },reused);
        if(!status) return status;
        cachedZeta[index]=true;
      } else std::fill_n(zetaFields[index],fieldElements,0.0);
    }
    WVRealFieldBundleView derivativeBundle{
        derivatives,
        {configuration.Nx, configuration.Ny, configuration.Nz, 3}};
    const auto evaluateDerivatives=[&](WVDynamicalField field) {
      bool reused=false;
      const auto portableField=field==WVDynamicalField::u ? WVPortableVariable::u :
          field==WVDynamicalField::v ? WVPortableVariable::v :
          field==WVDynamicalField::w ? WVPortableVariable::w :
          WVPortableVariable::eta;
      const auto operation=[&]() {return invokeTransform([&]() {
        return transform_->transformStateFieldDerivatives(state,field,derivativeBundle);
      });};
      std::array<WVVariableEvaluationKey,3> keys;
      std::array<WVFieldOutputView,3> derivativeViews;
      for(std::uint32_t axis=1;axis<=3;++axis) {
        keys[axis-1]={WVVariableEvaluationNode::derivative,
            static_cast<std::uint32_t>(portableField),
            eventWorkspace_ ? eventWorkspace_->component() : 0u,axis};
        derivativeViews[axis-1]={derivatives+(axis-1)*fieldElements,
            fieldElements};
      }
      const auto status=eventWorkspace_ ? eventWorkspace_->evaluateGroup(
          keys,derivativeViews,operation,reused) : operation();
      if(status) {
        if(reused) ++metrics_.primitiveFieldReuseCount;
        else ++metrics_.primitiveFieldEvaluationCount;
      }
      return status;
    };
    if ((dependencyMask & uDerivatives) != 0) {
      auto status = evaluateDerivatives(WVDynamicalField::u);
      if (!status)
        return status;
      const double *uy = derivatives + fieldElements;
      const double *uz = uy + fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        if(requestedZeta[1] && !cachedZeta[1]) zetaY[index] += uz[index];
        if(requestedZeta[2] && !cachedZeta[2]) zetaZ[index] -= uy[index];
      }
    }
    if ((dependencyMask & vDerivatives) != 0) {
      auto status = evaluateDerivatives(WVDynamicalField::v);
      if (!status)
        return status;
      const double *vx = derivatives;
      const double *vz = derivatives + 2 * fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        if(requestedZeta[0] && !cachedZeta[0]) zetaX[index] -= vz[index];
        if(requestedZeta[2] && !cachedZeta[2]) zetaZ[index] += vx[index];
      }
    }
    if ((dependencyMask & wDerivatives) != 0) {
      auto status = evaluateDerivatives(WVDynamicalField::w);
      if (!status)
        return status;
      const double *wx = derivatives;
      const double *wy = derivatives + fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        if(requestedZeta[0] && !cachedZeta[0]) zetaX[index] += wy[index];
        if(requestedZeta[1] && !cachedZeta[1]) zetaY[index] -= wx[index];
      }
    }
    for(std::size_t index=0;index<3;++index) {
      if(!requestedZeta[index] || cachedZeta[index]) continue;
      const auto produced=[&]() {
        ++outputProducerMetrics_.reconstructions[10+index][0]
            [eventWorkspace_ ? eventWorkspace_->component() : 0u];
        return WVKernelStatus::ok();
      };
      bool reused=false;
      const auto status=eventWorkspace_ ? eventWorkspace_->evaluate(
          {WVVariableEvaluationNode::physicalField,
            static_cast<std::uint32_t>(zetaVariables[index]),
            eventWorkspace_->component()},zetaFields[index],fieldElements,
          produced,reused) : produced();
      if(!status) return status;
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaX))
      writeField(WVFieldEvaluationPlan::Field::zetaX,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaX);
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaY))
      writeField(WVFieldEvaluationPlan::Field::zetaY,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaY);
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaZ))
      writeField(WVFieldEvaluationPlan::Field::zetaZ,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaZ);
  }

  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::samplePreparedField(
    const WVFieldEvaluationPlan &plan, const double *source,
    WVFieldOutputView output) {
  if (barotropicQG_)
    return barotropicQG_->samplePreparedField(plan, source, output);
  if (stratified_)
    return stratified_->samplePreparedField(plan, source, output);
  if (plan.requests_.size() != 1 || plan.outputs_.size() != 1 || !source ||
      !output.data || output.elementCount != plan.outputs_.front().elementCount)
    return {WVKernelStatusCode::invalidShape,
            "Prepared field sampler storage has the wrong shape."};
  if (!sameTransformConfiguration(plan.configuration_,
                                  transform_->descriptor().configuration()))
    return invalid("Prepared field sampler belongs to another transform.");
  const auto &request = plan.requests_.front();
  const auto &configuration = transform_->descriptor().configuration();
  const auto plane = configuration.Nx * configuration.Ny;
  if (request.samplingKind == WVFieldSamplingKind::fixedVerticalProfiles) {
    for (std::size_t profile = 0; profile < request.profileXIndices.size();
         ++profile) {
      const auto horizontal = request.profileXIndices[profile] +
                              configuration.Nx * request.profileYIndices[profile];
      for (std::size_t z = 0; z < configuration.Nz; ++z)
        output.data[z + configuration.Nz * profile] =
            source[horizontal + plane * z];
    }
    ++metrics_.profileWriteCount;
  } else if (request.samplingKind == WVFieldSamplingKind::positions) {
    const auto zCount = request.nativeRank == WVPortableNaturalRank::volume
                            ? configuration.Nz
                            : 1;
    for (std::size_t position = 0;
         position < request.positionWeights.size(); ++position) {
      const auto &weights = request.positionWeights[position];
      double value = 0.0;
      if (!weights.outsideInterpolationDomain) {
        if (request.interpolation == WVPositionInterpolation::linear) {
          const auto verticalCount = zCount == 1 ? 1u : 2u;
          for (std::size_t z = 0; z < verticalCount; ++z)
            for (std::size_t y = 0; y < 2; ++y)
              for (std::size_t x = 0; x < 2; ++x)
                value += source[weights.xLinearIndices[x] + configuration.Nx *
                                (weights.yLinearIndices[y] + configuration.Ny *
                                 weights.zLinearIndices[z])] *
                         weights.xLinearWeights[x] *
                         weights.yLinearWeights[y] *
                         (zCount == 1 ? 1.0 : weights.zLinearWeights[z]);
        } else {
          for (std::size_t z = 0; z < zCount; ++z)
            for (std::size_t y = 0; y < configuration.Ny; ++y)
              for (std::size_t x = 0; x < configuration.Nx; ++x)
                value += source[x + configuration.Nx *
                                (y + configuration.Ny * z)] *
                         weights.xSplineWeights[x] *
                         weights.ySplineWeights[y] *
                         (zCount == 1 ? 1.0 : weights.zSplineWeights[z]);
        }
      }
      output.data[position] = value;
    }
    if (request.interpolation == WVPositionInterpolation::linear)
      metrics_.linearInterpolationCount += request.positionWeights.size();
    else
      metrics_.splineInterpolationCount += request.positionWeights.size();
  } else {
    std::copy_n(source, output.elementCount, output.data);
    ++metrics_.fullGridWriteCount;
  }
  metrics_.outputElementWriteCount += output.elementCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::createEventPlan(
    const std::vector<WVEventFieldRequest> &requests,
    WVEventFieldEvaluationPlan &plan,
    WVDensityDiagnosticContract densityContract) {
  if (densityContract.reference != WVNoMotionReference::actual &&
      densityContract.reference != WVNoMotionReference::initial)
    return invalid("Invalid event-field density diagnostic reference.");
  std::vector<WVFieldRequest> diagnosticRequests;
  try {
    diagnosticRequests.reserve(requests.size());
    for (const auto &request : requests)
      diagnosticRequests.push_back({request.identifier, request.fieldName, {}});
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to inspect event-field requests."};
  }
  if (detail::WVDiagnosticFieldPlan::required(diagnosticRequests,
                                              stratified_ != nullptr)) {
    try {
      WVEventFieldEvaluationPlan candidate;
      candidate.genericSampling_ = true;
      candidate.owner_ = this;
      candidate.planIdentity_ = std::make_shared<const std::uint8_t>(0);
      candidate.configurationIdentifier_ = portableVariableConfiguration();
      candidate.densityContract_ = densityContract;
      candidate.requests_.reserve(requests.size());
      candidate.outputs_.reserve(requests.size());
      std::set<std::string> identifiers;
      for (std::size_t index = 0; index < requests.size(); ++index) {
        const auto &request = requests[index];
        if (request.identifier.empty() ||
            !identifiers.insert(request.identifier).second)
          return invalid("Event-field request identifiers must be nonempty and unique.");
        if (request.positionSetSlot == std::numeric_limits<std::size_t>::max())
          return {WVKernelStatusCode::sizeOverflow,
                  "An event-field position-set slot overflows its plan."};
        if (request.interpolation != WVPositionInterpolation::linear &&
            request.interpolation != WVPositionInterpolation::spline)
          return invalid("Event-field interpolation method is invalid.");
        const WVPortableVariableMetadata *metadata =
            findPortableVariable(request.fieldName);
        const WVPortableVariableContract *contract = nullptr;
        if (metadata) {
          contract = portableVariableContract(metadata->identifier,
                                               candidate.configurationIdentifier_);
        } else if (forcing_) {
          detail::WVForcingDiagnosticBinding::Output bound;
          const auto status = forcing_->resolve(
              request.fieldName, candidate.configurationIdentifier_,
              portablePositionSampling, bound);
          if (!status)
            return status;
          contract = bound.contract;
          metadata = contract ? &contract->metadata : nullptr;
        }
        if (!metadata ||
            (metadata->naturalRank != WVPortableNaturalRank::volume &&
             metadata->naturalRank != WVPortableNaturalRank::horizontal) ||
            (metadata->samplingMask & portablePositionSampling) == 0)
          return {WVKernelStatusCode::unsupportedOperation,
                  "Event-position sampling does not support field " +
                      request.fieldName + "."};
        if (metadata->ordinal >= 23 && !contract)
          return {WVKernelStatusCode::unsupportedOperation,
                  "Diagnostic is unavailable on this transform: " +
                      request.fieldName};
        candidate.requests_.push_back(
            {metadata->identifier, metadata->naturalRank,
             metadata->primitiveDependencyMask, request.positionSetSlot,
             request.interpolation, index});
        candidate.outputs_.push_back(
            {request.identifier, request.fieldName, metadata->identifier,
             metadata->naturalRank, metadata->primitiveDependencyMask,
             request.positionSetSlot, request.interpolation});
        candidate.positionSetCount_ = std::max(
            candidate.positionSetCount_, request.positionSetSlot + 1);
      }
      candidate.requiresZByPositionSet_.assign(candidate.positionSetCount_, 0);
      for (const auto &request : candidate.requests_)
        if (request.nativeRank == WVPortableNaturalRank::volume)
          candidate.requiresZByPositionSet_[request.positionSetSlot] = 1;
      std::uint64_t fingerprint = fingerprintOffset;
      appendFingerprint(fingerprint, candidate.configurationIdentifier_.data(),
                        candidate.configurationIdentifier_.size());
      appendFingerprint(fingerprint, candidate.positionSetCount_);
      for (std::size_t index = 0; index < requests.size(); ++index) {
        appendFingerprint(fingerprint, requests[index].fieldName.data(),
                          requests[index].fieldName.size());
        appendFingerprint(fingerprint, requests[index].positionSetSlot);
        appendFingerprint(fingerprint, requests[index].interpolation);
      }
      candidate.fingerprint_ = fingerprint;
      const auto planBytes = candidate.persistentBytes();
      const auto arenaStatus=prepareEventArena(candidate);
      if(!arenaStatus) return arenaStatus;
      plan = std::move(candidate);
      auto &eventMetrics = mutableMetrics();
      ++eventMetrics.eventPlanCreationCount;
      eventMetrics.eventPlanFieldResolutionCount += requests.size();
      eventMetrics.lastEventPlanBytes = planBytes;
      eventMetrics.maximumEventPlanBytes = std::max(
          eventMetrics.maximumEventPlanBytes, planBytes);
      return WVKernelStatus::ok();
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate a sampled event-field plan."};
    }
  }
  if (barotropicQG_) {
    const auto status = barotropicQG_->createEventPlan(requests, plan);
    if (status) {
      const auto arenaStatus=prepareEventArena(plan);
      if(!arenaStatus) return arenaStatus;
      plan.owner_ = this;
    }
    return status;
  }
  if (stratified_) {
    const auto status = stratified_->createEventPlan(requests, plan);
    if (status) {
      const auto arenaStatus=prepareEventArena(plan);
      if(!arenaStatus) return arenaStatus;
      plan.owner_ = this;
    }
    return status;
  }
  try {
    WVEventFieldEvaluationPlan candidate;
    candidate.owner_ = this;
    const auto &configuration = transform_->descriptor().configuration();
    candidate.configuration_ = configuration;
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    std::size_t positionSetCount = 0;
    for (std::size_t outputIndex = 0; outputIndex < requests.size();
         ++outputIndex) {
      const auto &request = requests[outputIndex];
      if (request.identifier.empty() ||
          !identifiers.insert(request.identifier).second)
        return invalid(
            "Event-field request identifiers must be nonempty and unique.");
      if (request.positionSetSlot ==
          std::numeric_limits<std::size_t>::max())
        return {WVKernelStatusCode::sizeOverflow,
                "An event-field position-set slot overflows its plan."};
      if (request.interpolation != WVPositionInterpolation::linear &&
          request.interpolation != WVPositionInterpolation::spline)
        return invalid("Event-field interpolation method is invalid.");
      const auto *metadata = findExecutablePortableVariable(request.fieldName);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field ||
          (metadata->samplingMask & portablePositionSampling) == 0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Event-position sampling does not support field " +
                    request.fieldName + "."};
      if (metadata->naturalRank != WVPortableNaturalRank::volume &&
          metadata->naturalRank != WVPortableNaturalRank::horizontal)
        return {WVKernelStatusCode::unsupportedOperation,
                "Event-position sampling requires a volume or horizontal "
                "field: " +
                    request.fieldName + "."};
      if (metadata->identifier == WVPortableVariable::psi &&
          transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Streamfunction evaluation is undefined when the Coriolis "
                "frequency is zero."};

      candidate.requests_.push_back(
          {metadata->identifier, metadata->naturalRank,
           metadata->primitiveDependencyMask, request.positionSetSlot,
           request.interpolation, outputIndex});
      candidate.outputs_.push_back(
          {request.identifier, request.fieldName, metadata->identifier,
           metadata->naturalRank, metadata->primitiveDependencyMask,
           request.positionSetSlot, request.interpolation});
      candidate.requestedFieldMask_ |=
          1ULL << static_cast<std::size_t>(metadata->ordinal);
      candidate.dependencyMask_ |= metadata->primitiveDependencyMask;
      positionSetCount =
          std::max(positionSetCount, request.positionSetSlot + 1);
    }
    candidate.positionSetCount_ = positionSetCount;
    candidate.requiresZByPositionSet_.assign(positionSetCount, 0);
    for (const auto &request : candidate.requests_)
      if (request.nativeRank == WVPortableNaturalRank::volume)
        candidate.requiresZByPositionSet_[request.positionSetSlot] = 1;

    std::uint64_t fingerprint = fingerprintOffset;
    appendConfigurationFingerprint(fingerprint, configuration);
    appendFingerprint(fingerprint, candidate.positionSetCount_);
    appendFingerprint(fingerprint, candidate.requestedFieldMask_);
    appendFingerprint(fingerprint, candidate.dependencyMask_);
    const auto requestCount = candidate.requests_.size();
    appendFingerprint(fingerprint, requestCount);
    for (std::size_t index = 0; index < candidate.requests_.size(); ++index) {
      const auto &resolved = candidate.requests_[index];
      const auto &output = candidate.outputs_[index];
      appendFingerprint(fingerprint, resolved.field);
      appendFingerprint(fingerprint, resolved.nativeRank);
      appendFingerprint(fingerprint, resolved.dependencyMask);
      appendFingerprint(fingerprint, resolved.positionSetSlot);
      appendFingerprint(fingerprint, resolved.interpolation);
      appendFingerprint(fingerprint, output.identifier.data(),
                        output.identifier.size());
      appendFingerprint(fingerprint, output.fieldName.data(),
                        output.fieldName.size());
    }
    candidate.fingerprint_ = fingerprint;
    const auto planBytes = candidate.persistentBytes();
    const auto arenaStatus=prepareEventArena(candidate);
    if(!arenaStatus) return arenaStatus;
    plan = std::move(candidate);
    ++metrics_.eventPlanCreationCount;
    metrics_.eventPlanFieldResolutionCount += requests.size();
    metrics_.lastEventPlanBytes = planBytes;
    metrics_.maximumEventPlanBytes =
        std::max(metrics_.maximumEventPlanBytes, planBytes);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate an event-field evaluation plan."};
  }
}

WVKernelStatus WVFieldEvaluationService::prepareEventGeometry(
    const WVEventFieldEvaluationPlan &plan,
    const WVEventPositionSetView *positionSets,
    std::size_t positionSetCount, WVPreparedFieldGeometry &geometry) {
  if (plan.owner_ != this)
    return invalid("The event-field plan belongs to another field service.");
  if (plan.genericSampling_) {
    if (plan.owner_ != this ||
        plan.configurationIdentifier_ != portableVariableConfiguration())
      return invalid("The sampled event-field plan belongs to another transform.");
    if (positionSetCount != plan.positionSetCount_)
      return {WVKernelStatusCode::invalidShape,
              "Event position-set count must match the sampled field plan."};
    if (positionSetCount && !positionSets)
      return {WVKernelStatusCode::invalidPointer,
              "Event position sets have a null view pointer."};
    try {
      WVPreparedFieldGeometry candidate;
      candidate.fieldPlanFingerprint_ = plan.fingerprint_;
      candidate.planIdentity_ = plan.planIdentity_;
      candidate.positionSets_.reserve(positionSetCount);
      std::uint64_t geometryFingerprint = fingerprintOffset;
      appendFingerprint(geometryFingerprint, plan.fingerprint_);
      for (std::size_t slot = 0; slot < positionSetCount; ++slot) {
        const auto &view = positionSets[slot];
        if (view.extentCount && !view.extents)
          return {WVKernelStatusCode::invalidPointer,
                  "Event position-set extents have a null pointer."};
        if (view.positionCount &&
            (!view.x || !view.y ||
             (plan.requiresZByPositionSet_[slot] && !view.z)))
          return {WVKernelStatusCode::invalidPointer,
                  "Event coordinates required by a sampled field are missing."};
        WVPreparedFieldGeometry::PositionSet stored;
        stored.x = view.x;
        stored.y = view.y;
        stored.z = view.z;
        stored.positionCount = view.positionCount;
        if (view.extentCount)
          stored.extents.assign(view.extents, view.extents + view.extentCount);
        else
          stored.extents = {view.positionCount};
        std::size_t extentProduct = 1;
        for (const auto extent : stored.extents)
          extentProduct = checkedProduct(extentProduct, extent);
        if (extentProduct != view.positionCount)
          return {WVKernelStatusCode::invalidShape,
                  "Event position-set extents do not match its sample count."};
        for (std::size_t position = 0; position < view.positionCount; ++position)
          if (!std::isfinite(view.x[position]) ||
              !std::isfinite(view.y[position]) ||
              (view.z && !std::isfinite(view.z[position])))
            return invalid("Event position coordinates must be finite.");
        if (candidate.positionCount_ >
            std::numeric_limits<std::size_t>::max() - view.positionCount)
          return {WVKernelStatusCode::sizeOverflow,
                  "Event position count overflows its prepared geometry."};
        candidate.positionCount_ += view.positionCount;
        const auto coordinateCount = view.z ? 3u : 2u;
        candidate.borrowedCoordinateBytes_ +=
            checkedProduct(checkedProduct(view.positionCount, coordinateCount),
                           sizeof(double));
        appendFingerprint(geometryFingerprint, slot);
        appendFingerprint(geometryFingerprint, view.positionCount);
        for (const auto extent : stored.extents)
          appendFingerprint(geometryFingerprint, extent);
        const auto bytes = view.positionCount * sizeof(double);
        appendFingerprint(geometryFingerprint, view.x, bytes);
        appendFingerprint(geometryFingerprint, view.y, bytes);
        if (view.z)
          appendFingerprint(geometryFingerprint, view.z, bytes);
        candidate.positionSets_.push_back(std::move(stored));
      }
      std::vector<WVFieldRequest> requests;
      requests.reserve(plan.requests_.size());
      candidate.outputs_.reserve(plan.requests_.size());
      for (const auto &request : plan.requests_) {
        const auto &set = candidate.positionSets_[request.positionSetSlot];
        WVFieldSamplingRequest sampling;
        sampling.kind = WVFieldSamplingKind::positions;
        sampling.interpolation = request.interpolation;
        if (set.positionCount) {
          sampling.x.assign(set.x, set.x + set.positionCount);
          sampling.y.assign(set.y, set.y + set.positionCount);
          if (request.nativeRank == WVPortableNaturalRank::volume)
            sampling.z.assign(set.z, set.z + set.positionCount);
        }
        requests.push_back(
            {plan.outputs_[request.outputIndex].identifier,
             plan.outputs_[request.outputIndex].fieldName,
             std::move(sampling)});
        candidate.outputs_.push_back(
            {request.outputIndex, request.positionSetSlot, set.extents,
             set.positionCount});
      }
      auto status = createPlan(requests, candidate.evaluationPlan_,
                               plan.densityContract_);
      if (!status)
        return status;
      candidate.geometryFingerprint_ = geometryFingerprint;
      const auto retainedBytes = candidate.retainedBytes();
      const auto liveBytes = candidate.liveBytes();
      geometry = std::move(candidate);
      auto &eventMetrics = mutableMetrics();
      ++eventMetrics.eventGeometryPreparationCount;
      eventMetrics.eventPositionSetCount += positionSetCount;
      eventMetrics.eventPositionCount += geometry.positionCount_;
      eventMetrics.lastPreparedGeometryRetainedBytes = retainedBytes;
      eventMetrics.maximumPreparedGeometryRetainedBytes = std::max(
          eventMetrics.maximumPreparedGeometryRetainedBytes, retainedBytes);
      eventMetrics.lastPreparedGeometryLiveBytes = liveBytes;
      eventMetrics.maximumPreparedGeometryLiveBytes = std::max(
          eventMetrics.maximumPreparedGeometryLiveBytes, liveBytes);
      return WVKernelStatus::ok();
    } catch (const std::overflow_error &) {
      return {WVKernelStatusCode::sizeOverflow,
              "Event geometry extents or storage overflow size_t."};
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate sampled event geometry."};
    }
  }
  if (barotropicQG_)
    return barotropicQG_->prepareEventGeometry(
        plan, positionSets, positionSetCount, geometry);
  if (stratified_)
    return stratified_->prepareEventGeometry(
        plan, positionSets, positionSetCount, geometry);
  if (!sameTransformConfiguration(
          plan.configuration_, transform_->descriptor().configuration()))
    return invalid(
        "The event-field plan belongs to a different transform configuration.");
  if (positionSetCount != plan.positionSetCount_)
    return {WVKernelStatusCode::invalidShape,
            "Event position-set count must match the resolved field plan."};
  if (positionSetCount != 0 && positionSets == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Event position sets have a null view pointer."};

  try {
    WVPreparedFieldGeometry candidate;
    candidate.fieldPlanFingerprint_ = plan.fingerprint_;
    candidate.positionSets_.reserve(positionSetCount);
    std::uint64_t geometryFingerprint = fingerprintOffset;
    appendFingerprint(geometryFingerprint, plan.fingerprint_);
    appendFingerprint(geometryFingerprint, positionSetCount);
    for (std::size_t slot = 0; slot < positionSetCount; ++slot) {
      const auto &view = positionSets[slot];
      if (view.extentCount != 0 && view.extents == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Event position-set extents have a null pointer."};
      if (view.positionCount != 0 &&
          (view.x == nullptr || view.y == nullptr ||
           (plan.requiresZByPositionSet_[slot] != 0 && view.z == nullptr)))
        return {WVKernelStatusCode::invalidPointer,
                "Event coordinates required by a resolved field are missing."};

      WVPreparedFieldGeometry::PositionSet prepared;
      prepared.x = view.x;
      prepared.y = view.y;
      prepared.z = view.z;
      prepared.positionCount = view.positionCount;
      if (view.extentCount == 0)
        prepared.extents.push_back(view.positionCount);
      else
        prepared.extents.assign(view.extents,
                                view.extents + view.extentCount);
      std::size_t extentProduct = 1;
      for (const auto extent : prepared.extents)
        extentProduct = checkedProduct(extentProduct, extent);
      if (extentProduct != view.positionCount)
        return {WVKernelStatusCode::invalidShape,
                "Event position-set extents do not match its sample count."};
      for (std::size_t position = 0; position < view.positionCount;
           ++position) {
        if (!std::isfinite(view.x[position]) ||
            !std::isfinite(view.y[position]) ||
            (view.z != nullptr && !std::isfinite(view.z[position])))
          return invalid("Event position coordinates must be finite.");
      }
      if (candidate.positionCount_ >
          std::numeric_limits<std::size_t>::max() - view.positionCount)
        return {WVKernelStatusCode::sizeOverflow,
                "Event position count overflows its prepared geometry."};
      candidate.positionCount_ += view.positionCount;
      const std::size_t coordinateArrayCount = view.z == nullptr ? 2 : 3;
      const auto coordinateBytes = checkedProduct(
          checkedProduct(view.positionCount, coordinateArrayCount),
          sizeof(double));
      if (candidate.borrowedCoordinateBytes_ >
          std::numeric_limits<std::size_t>::max() - coordinateBytes)
        return {WVKernelStatusCode::sizeOverflow,
                "Event coordinate storage overflows its metrics."};
      candidate.borrowedCoordinateBytes_ += coordinateBytes;

      appendFingerprint(geometryFingerprint, slot);
      appendFingerprint(geometryFingerprint, view.positionCount);
      const auto extentCount = prepared.extents.size();
      appendFingerprint(geometryFingerprint, extentCount);
      for (const auto extent : prepared.extents)
        appendFingerprint(geometryFingerprint, extent);
      const bool hasZ = view.z != nullptr;
      appendFingerprint(geometryFingerprint, hasZ);
      const auto coordinateByteCount = view.positionCount * sizeof(double);
      appendFingerprint(geometryFingerprint, view.x, coordinateByteCount);
      appendFingerprint(geometryFingerprint, view.y, coordinateByteCount);
      if (hasZ)
        appendFingerprint(geometryFingerprint, view.z, coordinateByteCount);
      candidate.positionSets_.push_back(std::move(prepared));
    }
    candidate.geometryFingerprint_ = geometryFingerprint;

    const auto &configuration = transform_->descriptor().configuration();
    candidate.evaluationPlan_.configuration_ = configuration;
    candidate.evaluationPlan_.requestedFieldMask_ = plan.requestedFieldMask_;
    candidate.evaluationPlan_.dependencyMask_ = plan.dependencyMask_;
    candidate.evaluationPlan_.requests_.reserve(plan.requests_.size());
    candidate.evaluationPlan_.outputs_.reserve(plan.outputs_.size());
    candidate.outputs_.reserve(plan.outputs_.size());
    const double dx =
        configuration.Lx / static_cast<double>(configuration.Nx);
    const double dy =
        configuration.Ly / static_cast<double>(configuration.Ny);
    const double dz =
        configuration.Lz / static_cast<double>(configuration.Nz - 1);

    for (const auto &eventRequest : plan.requests_) {
      const auto &set =
          candidate.positionSets_[eventRequest.positionSetSlot];
      WVFieldEvaluationPlan::ResolvedRequest request;
      request.field = eventRequest.field;
      request.dependencyMask = eventRequest.dependencyMask;
      request.nativeRank = eventRequest.nativeRank;
      request.samplingKind = WVFieldSamplingKind::positions;
      request.interpolation = eventRequest.interpolation;
      request.outputIndex = eventRequest.outputIndex;
      request.positionWeights.reserve(set.positionCount);
      for (std::size_t position = 0; position < set.positionCount;
           ++position) {
        const double x = set.x[position];
        const double y = set.y[position];
        const double z = eventRequest.nativeRank ==
                                 WVPortableNaturalRank::volume
                             ? set.z[position]
                             : 0.0;
        const double xWrapped = wrapped(x, configuration.Lx);
        const double yWrapped = wrapped(y, configuration.Ly);
        const auto xLower = std::min(
            static_cast<std::size_t>(std::floor(xWrapped / dx)),
            configuration.Nx - 1);
        const auto yLower = std::min(
            static_cast<std::size_t>(std::floor(yWrapped / dy)),
            configuration.Ny - 1);
        WVFieldEvaluationPlan::PositionWeights weights;
        if (eventRequest.interpolation == WVPositionInterpolation::linear) {
          weights.xLinearIndices = {xLower,
                                    (xLower + 1) % configuration.Nx};
          weights.yLinearIndices = {yLower,
                                    (yLower + 1) % configuration.Ny};
          const double xFraction =
              (xWrapped - static_cast<double>(xLower) * dx) / dx;
          const double yFraction =
              (yWrapped - static_cast<double>(yLower) * dy) / dy;
          weights.xLinearWeights = {1.0 - xFraction, xFraction};
          weights.yLinearWeights = {1.0 - yFraction, yFraction};
          if (eventRequest.nativeRank == WVPortableNaturalRank::volume) {
            weights.outsideInterpolationDomain =
                z < -configuration.Lz || z > 0.0;
            if (!weights.outsideInterpolationDomain) {
              const double normalizedZ = (z + configuration.Lz) / dz;
              const auto zLower = std::min(
                  static_cast<std::size_t>(
                      std::max(0.0, std::floor(normalizedZ))),
                  configuration.Nz - 2);
              const double zFraction = std::clamp(
                  normalizedZ - static_cast<double>(zLower), 0.0, 1.0);
              weights.zLinearIndices = {zLower, zLower + 1};
              weights.zLinearWeights = {1.0 - zFraction, zFraction};
            }
          }
        } else {
          const bool xBoundary =
              xLower < 3 || xLower > configuration.Nx - 4;
          const bool yBoundary =
              yLower < 3 || yLower > configuration.Ny - 4;
          const std::size_t xShift = xBoundary ? 4 : 0;
          const std::size_t yShift = yBoundary ? 4 : 0;
          const double xQuery =
              xBoundary ? wrapped(x + 4.0 * dx, configuration.Lx) : xWrapped;
          const double yQuery =
              yBoundary ? wrapped(y + 4.0 * dy, configuration.Ly) : yWrapped;
          weights.outsideInterpolationDomain =
              xQuery > static_cast<double>(configuration.Nx - 1) * dx ||
              yQuery > static_cast<double>(configuration.Ny - 1) * dy;
          if (!weights.outsideInterpolationDomain) {
            movingWorkspace_->xSpline.weightsInto(
                0.0, dx, xQuery, weights.xSplineWeights, xShift,
                &movingWorkspace_->xRightHandSide,
                &movingWorkspace_->xShifted);
            movingWorkspace_->ySpline.weightsInto(
                0.0, dy, yQuery, weights.ySplineWeights, yShift,
                &movingWorkspace_->yRightHandSide,
                &movingWorkspace_->yShifted);
          }
          if (eventRequest.nativeRank == WVPortableNaturalRank::volume) {
            weights.outsideInterpolationDomain =
                weights.outsideInterpolationDomain ||
                z < -configuration.Lz || z > 0.0;
            if (!weights.outsideInterpolationDomain)
              movingWorkspace_->zSpline.weightsInto(
                  -configuration.Lz, dz, z, weights.zSplineWeights, 0,
                  &movingWorkspace_->zRightHandSide,
                  &movingWorkspace_->zShifted);
          }
        }
        request.positionWeights.push_back(std::move(weights));
      }
      candidate.evaluationPlan_.requests_.push_back(std::move(request));
      candidate.evaluationPlan_.outputs_.push_back(
          {std::string{}, std::string{}, WVFieldSamplingKind::positions, {},
           set.positionCount});
      candidate.outputs_.push_back(
          {eventRequest.outputIndex, eventRequest.positionSetSlot,
           set.extents, set.positionCount});
    }

    const auto retainedBytes = candidate.retainedBytes();
    const auto liveBytes = candidate.liveBytes();
    geometry = std::move(candidate);
    ++metrics_.eventGeometryPreparationCount;
    metrics_.eventPositionSetCount += positionSetCount;
    metrics_.eventPositionCount += geometry.positionCount_;
    metrics_.lastPreparedGeometryRetainedBytes = retainedBytes;
    metrics_.maximumPreparedGeometryRetainedBytes =
        std::max(metrics_.maximumPreparedGeometryRetainedBytes,
                 retainedBytes);
    metrics_.lastPreparedGeometryLiveBytes = liveBytes;
    metrics_.maximumPreparedGeometryLiveBytes =
        std::max(metrics_.maximumPreparedGeometryLiveBytes, liveBytes);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &) {
    return {WVKernelStatusCode::sizeOverflow,
            "Event geometry extents or storage overflow size_t."};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate prepared event geometry."};
  }
}

WVKernelStatus WVFieldEvaluationService::evaluateEvent(
    const WVEventFieldEvaluationPlan &plan,
    const WVPreparedFieldGeometry &geometry, const WVState &state,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  const WVEventFieldEvaluationBatchEntry entry{&plan, &geometry, outputs,
                                               outputCount};
  return evaluateEventBatch(state, &entry, 1);
}

WVKernelStatus WVFieldEvaluationService::evaluateEvent(
    const WVEventFieldEvaluationPlan &plan,
    const WVPreparedFieldGeometry &geometry,
    const WVIntegrationState &state, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  const WVEventFieldEvaluationBatchEntry entry{&plan, &geometry, outputs,
                                               outputCount};
  return evaluateEventBatch(state, &entry, 1);
}

WVKernelStatus WVFieldEvaluationService::evaluateEventBatch(
    const WVState &state, const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entries && std::any_of(entries, entries + entryCount, [](const auto &entry) {
        return entry.plan && entry.plan->genericSampling_;
      }))
    return evaluateSampledEventBatch({state}, entries, entryCount);
  if (!transform_) return {WVKernelStatusCode::unsupportedOperation,"This transform requires coefficient-family state views."};
  if (entryCount != 0 && entries == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Event field batch has a null entry pointer."};
  std::size_t batchOutputCount = 0;
  for (std::size_t entryIndex = 0; entryIndex < entryCount; ++entryIndex) {
    const auto &entry = entries[entryIndex];
    if (entry.plan == nullptr || entry.geometry == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Event field batch entry has a null plan or geometry."};
    const auto &plan = *entry.plan;
    const auto &geometry = *entry.geometry;
    if (!sameTransformConfiguration(plan.configuration_,
                                    transform_->descriptor().configuration()))
      return invalid("The event-field plan belongs to a different transform "
                     "configuration.");
    if (geometry.fieldPlanFingerprint_ != plan.fingerprint_ ||
        geometry.outputCount() != plan.outputCount())
      return invalid(
          "Prepared event geometry does not match the resolved field plan.");
    if (entry.outputCount != geometry.outputs_.size())
      return {WVKernelStatusCode::invalidShape,
              "Event output-view count must match the prepared geometry."};
    if (entry.outputCount != 0 && entry.outputs == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Event outputs have a null view pointer."};
    for (std::size_t outputIndex = 0; outputIndex < entry.outputCount;
         ++outputIndex)
      if (entry.outputs[outputIndex].elementCount !=
          geometry.outputs_[outputIndex].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Event output has the wrong element count for request " +
                    plan.outputs_[outputIndex].identifier + "."};
    batchOutputCount =
        batchOutputCount >
                std::numeric_limits<std::size_t>::max() - entry.outputCount
            ? std::numeric_limits<std::size_t>::max()
            : batchOutputCount + entry.outputCount;
  }

  try {
    eventBatchInvocations_.clear();
    eventBatchInvocations_.reserve(entryCount);
    for (std::size_t entryIndex = 0; entryIndex < entryCount; ++entryIndex) {
      const auto &entry = entries[entryIndex];
      eventBatchInvocations_.push_back(
          {&entry.geometry->evaluationPlan_, entry.outputs, entry.outputCount});
    }
  } catch (const std::bad_alloc &) {
    eventBatchInvocations_.clear();
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate event field batch invocation storage."};
  }
  metrics_.eventBatchInvocationWorkspaceBytes =
      eventBatchInvocations_.capacity() * sizeof(PlanInvocation);
  metrics_.servicePersistentBytes = persistentBytes();
  const auto status = evaluatePlanBatch(eventBatchInvocations_.data(),
                                        eventBatchInvocations_.size(), state);
  eventBatchInvocations_.clear();
  if (status) {
    metrics_.eventEvaluationCount += entryCount;
    ++metrics_.eventBatchEvaluationCount;
    metrics_.eventBatchOccurrenceCount += entryCount;
    metrics_.eventBatchOutputCount += batchOutputCount;
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::evaluateEventBatch(
    const WVIntegrationState &state,
    const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entries && std::any_of(entries, entries + entryCount, [](const auto &entry) {
        return entry.plan && entry.plan->genericSampling_;
      }))
    return evaluateSampledEventBatch(state, entries, entryCount);
  if (barotropicQG_)
    return barotropicQG_->evaluateEventBatch(state, entries, entryCount);
  if (stratified_)
    return stratified_->evaluateEventBatch(state, entries, entryCount);
  return evaluateEventBatch(state.waveVortex, entries, entryCount);
}

WVKernelStatus WVFieldEvaluationService::evaluateSampledEventBatch(
    const WVIntegrationState &state,
    const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entryCount && !entries)
    return {WVKernelStatusCode::invalidPointer,
            "Sampled event field batch has a null entry pointer."};
  if ((state.coefficientFamilyCount && !state.coefficientFamilies) ||
      (state.additionalBlockCount && !state.additionalBlocks))
    return {WVKernelStatusCode::invalidPointer,
            "Sampled event state views have a null array pointer."};
  std::size_t outputCount = 0;
  for (std::size_t index = 0; index < entryCount; ++index) {
    const auto &entry = entries[index];
    if (!entry.plan || !entry.geometry)
      return {WVKernelStatusCode::invalidPointer,
              "Sampled event field batch entry has a null plan or geometry."};
    if (entry.plan->owner_ != this)
      return invalid("Sampled event field batch contains a foreign plan.");
    if (entry.plan->genericSampling_ &&
        (
         entry.plan->configurationIdentifier_ != portableVariableConfiguration()))
      return invalid("Sampled event field batch contains an incompatible plan.");
    if ((entry.plan->genericSampling_ &&
         entry.geometry->planIdentity_ != entry.plan->planIdentity_) ||
        entry.geometry->fieldPlanFingerprint_ != entry.plan->fingerprint_ ||
        entry.geometry->outputCount() != entry.plan->outputCount())
      return invalid("Prepared sampled event geometry does not match its plan.");
    if (entry.outputCount != entry.geometry->outputs_.size() ||
        (entry.outputCount && !entry.outputs))
      return {WVKernelStatusCode::invalidShape,
              "Sampled event outputs must match prepared geometry."};
    for (std::size_t output = 0; output < entry.outputCount; ++output) {
      const auto &view = entry.outputs[output];
      if (view.elementCount !=
          entry.geometry->outputs_[output].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "A sampled event output has the wrong shape."};
      if (view.elementCount && !view.data)
        return {WVKernelStatusCode::invalidPointer,
                "A sampled event output has a null pointer."};
      const auto bytes = view.elementCount * sizeof(double);
      for (const auto coefficient :
           {state.waveVortex.coefficients.Ap, state.waveVortex.coefficients.Am,
            state.waveVortex.coefficients.A0})
        if (memoryOverlaps(view.data, bytes, coefficient.data,
                           coefficient.shape.elementCount() *
                               sizeof(WVComplex64)))
          return {WVKernelStatusCode::overlappingArrays,
                  "Sampled event outputs must not overlap coefficient state."};
      for (std::size_t family = 0; family < state.coefficientFamilyCount;
           ++family) {
        const auto &coefficient = state.coefficientFamilies[family];
        if (coefficient.layout &&
            memoryOverlaps(view.data, bytes, coefficient.data,
                           coefficient.layout->elementCount *
                               sizeof(WVComplex64)))
          return {WVKernelStatusCode::overlappingArrays,
                  "Sampled event outputs must not overlap coefficient state."};
      }
      for (std::size_t block = 0; block < state.additionalBlockCount; ++block) {
        const auto &additional = state.additionalBlocks[block];
        if (additional.layout &&
            (memoryOverlaps(view.data, bytes, additional.realData,
                            additional.layout->elementCount * sizeof(double)) ||
             memoryOverlaps(view.data, bytes, additional.complexData,
                            additional.layout->elementCount *
                                sizeof(WVComplex64))))
          return {WVKernelStatusCode::overlappingArrays,
                  "Sampled event outputs must not overlap additional state."};
      }
      for (std::size_t priorEntry = 0; priorEntry <= index; ++priorEntry) {
        const auto lastOutput =
            priorEntry == index ? output : entries[priorEntry].outputCount;
        for (std::size_t priorOutput = 0; priorOutput < lastOutput;
             ++priorOutput) {
          const auto &prior = entries[priorEntry].outputs[priorOutput];
          if (memoryOverlaps(view.data, bytes, prior.data,
                             prior.elementCount * sizeof(double)))
            return {WVKernelStatusCode::overlappingArrays,
                    "Sampled event outputs must not overlap each other."};
        }
      }
    }
    if (outputCount > std::numeric_limits<std::size_t>::max() -
                          entry.outputCount)
      return {WVKernelStatusCode::sizeOverflow,
              "Sampled event output count overflows size_t."};
    outputCount += entry.outputCount;
  }
  if (entryCount == 0)
    return WVKernelStatus::ok();
  const bool ownsEventScope = entryCount != 0 && eventWorkspace_ == nullptr;
  detail::WVFieldEvaluationEventScope scope(*this, state, ownsEventScope,
                                            true);
  if (!scope.status())
    return scope.status();
  if (entryCount != 0) {
    const auto stateStatus = eventWorkspace_->validateState(state);
    if (!stateStatus)
      return stateStatus;
  }
  try {
    std::vector<std::vector<std::vector<double>>> staged(entryCount);
    std::vector<std::vector<WVFieldOutputView>> stagedViews(entryCount);
    std::size_t stagingBytes =
        staged.capacity() * sizeof(std::vector<std::vector<double>>) +
        stagedViews.capacity() * sizeof(std::vector<WVFieldOutputView>);
    for (std::size_t index = 0; index < entryCount; ++index) {
      const auto &entry = entries[index];
      staged[index].resize(entry.outputCount);
      stagedViews[index].resize(entry.outputCount);
      stagingBytes +=
          staged[index].capacity() * sizeof(std::vector<double>) +
          stagedViews[index].capacity() * sizeof(WVFieldOutputView);
      for (std::size_t output = 0; output < entry.outputCount; ++output) {
        staged[index][output].resize(entry.outputs[output].elementCount);
        stagingBytes += staged[index][output].capacity() * sizeof(double);
        stagedViews[index][output] = {staged[index][output].data(),
                                      staged[index][output].size()};
      }
    }
    auto &workspaceMetrics = mutableMetrics();
    const auto previousExternalBytes =
        eventWorkspace_->externalWorkspaceBytes();
    if (previousExternalBytes > std::numeric_limits<std::size_t>::max() -
                                    stagingBytes)
      return {WVKernelStatusCode::sizeOverflow,
              "Sampled event staging storage overflows size_t."};
    struct RestoreExternalWorkspace {
      detail::WVFieldEvaluationEventWorkspace *workspace;
      std::size_t bytes;
      ~RestoreExternalWorkspace() {
        workspace->setExternalWorkspaceBytes(bytes);
      }
    } restoreExternal{eventWorkspace_, previousExternalBytes};
    eventWorkspace_->setExternalWorkspaceBytes(previousExternalBytes +
                                               stagingBytes);
    for (std::size_t index = 0; index < entryCount; ++index) {
      const auto &entry = entries[index];
      const auto status = evaluate(entry.geometry->evaluationPlan_, state,
                                   stagedViews[index].data(),
                                   stagedViews[index].size());
      if (!status)
        return status;
    }
    for (std::size_t index = 0; index < entryCount; ++index)
      for (std::size_t output = 0; output < entries[index].outputCount;
           ++output)
        std::copy(staged[index][output].begin(), staged[index][output].end(),
                  entries[index].outputs[output].data);
    workspaceMetrics.eventEvaluationCount += entryCount;
    ++workspaceMetrics.eventBatchEvaluationCount;
    workspaceMetrics.eventBatchOccurrenceCount += entryCount;
    workspaceMetrics.eventBatchOutputCount += outputCount;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate sampled event batch staging."};
  }
}

WVKernelStatus WVFieldEvaluationService::createMovingPlan(
    const std::vector<WVMovingFieldRequest> &requests,
    WVMovingFieldEvaluationPlan &plan,
    WVDensityDiagnosticContract densityContract) const {
  if (densityContract.reference != WVNoMotionReference::actual &&
      densityContract.reference != WVNoMotionReference::initial)
    return invalid("Invalid moving-field density diagnostic reference.");
  std::vector<WVFieldRequest> diagnosticRequests;
  try {
    diagnosticRequests.reserve(requests.size());
    for (const auto &request : requests)
      diagnosticRequests.push_back({request.identifier, request.fieldName, {}});
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to inspect moving-field requests."};
  }
  const bool sampledPlanRequired =
      detail::WVDiagnosticFieldPlan::required(diagnosticRequests,
                                              stratified_ != nullptr) ||
      std::any_of(requests.begin(), requests.end(), [](const auto &request) {
        const auto *metadata = findPortableVariable(request.fieldName);
        return metadata && metadata->movingPrimitiveChannel < 0;
      });
  if (sampledPlanRequired) {
    try {
      auto implementation =
          std::make_shared<detail::WVSampledMovingFieldPlan>();
      implementation->configurationIdentifier = portableVariableConfiguration();
      implementation->owner = this;
      implementation->densityContract = densityContract;
      implementation->requests.reserve(requests.size());
      std::set<std::string> identifiers;
      std::size_t positionCount = 0;
      for (std::size_t index = 0; index < requests.size(); ++index) {
        const auto &request = requests[index];
        if (request.identifier.empty() ||
            !identifiers.insert(request.identifier).second ||
            request.positionCount == 0 ||
            request.positionOffset > std::numeric_limits<std::size_t>::max() -
                                         request.positionCount)
          return invalid("Sampled moving-field request is invalid.");
        if (request.interpolation != WVPositionInterpolation::linear &&
            request.interpolation != WVPositionInterpolation::spline)
          return invalid("Moving-field interpolation method is invalid.");
        const WVPortableVariableMetadata *metadata =
            findPortableVariable(request.fieldName);
        const WVPortableVariableContract *contract = nullptr;
        if (metadata) {
          contract = portableVariableContract(metadata->identifier,
                                               implementation->configurationIdentifier);
        } else if (forcing_) {
          detail::WVForcingDiagnosticBinding::Output bound;
          const auto status = forcing_->resolve(
              request.fieldName, implementation->configurationIdentifier,
              portablePositionSampling, bound);
          if (!status)
            return status;
          contract = bound.contract;
          metadata = contract ? &contract->metadata : nullptr;
        }
        if (!metadata || metadata->kind == WVPortableVariableKind::coefficient ||
            (metadata->naturalRank != WVPortableNaturalRank::volume &&
             metadata->naturalRank != WVPortableNaturalRank::horizontal) ||
            (metadata->samplingMask & portablePositionSampling) == 0)
          return {WVKernelStatusCode::unsupportedOperation,
                  "Moving-position sampling does not support field " +
                      request.fieldName + "."};
        if (metadata->ordinal >= 23 && !contract)
          return {WVKernelStatusCode::unsupportedOperation,
                  "Diagnostic is unavailable on this transform: " +
                      request.fieldName};
        implementation->requests.push_back(
            {request.identifier, request.fieldName, metadata->naturalRank,
             request.positionOffset, request.positionCount,
             request.interpolation, index});
        positionCount = std::max(positionCount,
                                 request.positionOffset + request.positionCount);
      }
      auto status = createPlan(diagnosticRequests,
                               implementation->fullGridPlan, densityContract);
      if (!status)
        return status;
      if(!sampledMovingWorkspace_)
        sampledMovingWorkspace_=std::make_unique<SampledMovingWorkspace>();
      status=sampledMovingWorkspace_->prepare(
          implementation->fullGridPlan,implementation->requests);
      if(!status) return status;
      WVMovingFieldEvaluationPlan candidate;
      candidate.positionCount_ = positionCount;
      candidate.sampledPlan_ = implementation;
      candidate.outputs_.reserve(requests.size());
      for (const auto &request : requests)
        candidate.outputs_.push_back(
            {request.identifier, request.fieldName,
             WVFieldSamplingKind::positions, {request.positionCount},
             request.positionCount});
      plan = std::move(candidate);
      return WVKernelStatus::ok();
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate a sampled moving-field plan."};
    }
  }
  if (barotropicQG_ || stratified_) {
    const auto status=barotropicQG_ ?
        barotropicQG_->createMovingPlan(requests,plan) :
        stratified_->createMovingPlan(requests,plan);
    if(!status) return status;
    WVFieldEvaluationPlan fields;
    try {
      fields.outputs_.reserve(requests.size());
      for(const auto& request:requests)
        fields.outputs_.push_back({request.identifier,request.fieldName,
            WVFieldSamplingKind::fullGrid,{},0,false});
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare moving output-event dependencies."};
    }
    return prepareEventArena(fields);
  }
  try {
    WVMovingFieldEvaluationPlan candidate;
    candidate.configuration_ = transform_->descriptor().configuration();
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &request = requests[index];
      if (request.identifier.empty() ||
          !identifiers.insert(request.identifier).second)
        return invalid(
            "Moving-field request identifiers must be nonempty and unique.");
      if (request.positionCount == 0 ||
          request.positionOffset >
              std::numeric_limits<std::size_t>::max() - request.positionCount)
        return invalid("Moving-field request position range is invalid.");
      if (request.interpolation != WVPositionInterpolation::linear &&
          request.interpolation != WVPositionInterpolation::spline)
        return invalid("Moving-field interpolation method is invalid.");
      const auto *metadata = findExecutablePortableVariable(request.fieldName);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field ||
          metadata->movingPrimitiveChannel < 0 ||
          (metadata->samplingMask & portablePositionSampling) == 0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Moving-position sampling does not support field " +
                    request.fieldName + "."};
      const auto channel =
          static_cast<std::size_t>(metadata->movingPrimitiveChannel);
      candidate.positionCount_ =
          std::max(candidate.positionCount_,
                   request.positionOffset + request.positionCount);
      candidate.requests_.push_back({channel, request.positionOffset,
                                     request.positionCount,
                                     request.interpolation, index});
      candidate.outputs_.push_back(
          {request.identifier, request.fieldName,
           WVFieldSamplingKind::positions, {request.positionCount},
           request.positionCount});
    }
    WVFieldEvaluationPlan fields;
    fields.dependencyMask_=primitiveValues;
    const auto arenaStatus=prepareEventArena(fields);
    if(!arenaStatus) return arenaStatus;
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a moving-field evaluation plan."};
  }
}

WVKernelStatus WVFieldEvaluationService::evaluateMoving(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const bool anyActive=outputCount && (!activeOutputs ||
      std::any_of(activeOutputs,activeOutputs+outputCount,
          [](std::uint8_t selected) {return selected!=0;}));
  if(!eventWorkspace_ && anyActive) {
    if (!transform_)
      return {WVKernelStatusCode::unsupportedOperation,
          "This transform requires coefficient-family state views."};
    detail::WVFieldEvaluationEventScope scope(*this,{state});
    if(!scope.status()) return scope.status();
    return evaluateMoving(plan,state,positions,outputs,outputCount,
        activeOutputs);
  }
  if (plan.sampledPlan_)
    return evaluateSampledMovingImpl(plan, {state}, positions, outputs,
                                     outputCount, activeOutputs);
  return evaluateMovingImpl(plan,state,nullptr,positions,outputs,outputCount,activeOutputs);
}

WVKernelStatus WVFieldEvaluationService::evaluateMoving(
    const WVMovingFieldEvaluationPlan &plan,
    const WVIntegrationState &state, WVMovingPositionView positions,
    WVFieldOutputView *outputs, std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const bool anyActive=outputCount && (!activeOutputs ||
      std::any_of(activeOutputs,activeOutputs+outputCount,
          [](std::uint8_t selected) {return selected!=0;}));
  if(!eventWorkspace_ && anyActive) {
    detail::WVFieldEvaluationEventScope scope(*this,state);
    if(!scope.status()) return scope.status();
    return evaluateMoving(plan,state,positions,outputs,outputCount,
        activeOutputs);
  }
  if(eventWorkspace_) {
    const auto status=eventWorkspace_->validateState(state);
    if(!status) return status;
  }
  if (plan.sampledPlan_)
    return evaluateSampledMovingImpl(plan, state, positions, outputs,
                                     outputCount, activeOutputs);
  if (barotropicQG_)
    return barotropicQG_->evaluateMoving(plan, state, positions, outputs,
                                         outputCount, activeOutputs);
  if (stratified_)
    return stratified_->evaluateMoving(plan, state, positions, outputs,
                                         outputCount, activeOutputs);
  return evaluateMoving(plan, state.waveVortex, positions, outputs,
                        outputCount, activeOutputs);
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  if(eventWorkspace_) {
    const auto status=eventWorkspace_->validateState({state});
    if(!status) return status;
  }
  if (plan.sampledPlan_)
    return {WVKernelStatusCode::unsupportedOperation,
            "Derived moving fields cannot be evaluated from prepared "
            "advection fields."};
  return evaluateMovingImpl(plan,state,&advectionFields,positions,outputs,outputCount,activeOutputs);
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &plan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  if (plan.sampledPlan_)
    return {WVKernelStatusCode::unsupportedOperation,
            "Derived moving fields cannot be evaluated from prepared "
            "advection fields."};
  if (barotropicQG_)
    return barotropicQG_->evaluateMovingFromAdvectionFields(
        plan, state, advectionFields, positions, outputs, outputCount, activeOutputs);
  if (stratified_)
    return stratified_->evaluateMovingFromAdvectionFields(
        plan, state, advectionFields, positions, outputs, outputCount, activeOutputs);
  return evaluateMovingFromAdvectionFields(
      plan, state.waveVortex, advectionFields, positions, outputs,
      outputCount, activeOutputs);
}

WVKernelStatus WVFieldEvaluationService::evaluateSampledMovingImpl(
    const WVMovingFieldEvaluationPlan &plan,
    const WVIntegrationState &state, WVMovingPositionView positions,
    WVFieldOutputView *outputs, std::size_t outputCount,
    const std::uint8_t *activeOutputs) {
  const auto implementation = plan.sampledPlan_;
  if (!implementation || implementation->owner != this ||
      implementation->configurationIdentifier != portableVariableConfiguration())
    return invalid("The sampled moving-field plan belongs to another transform.");
  if (positions.positionCount != plan.positionCount_ ||
      (positions.positionCount &&
       (!positions.x || !positions.y)))
    return {WVKernelStatusCode::invalidShape,
            "Moving coordinates must match the sampled field plan."};
  if (outputCount != plan.outputs_.size() || (outputCount && !outputs))
    return {WVKernelStatusCode::invalidShape,
            "Moving-field outputs must match the sampled field plan."};
  bool anyActive = false;
  for (const auto &request : implementation->requests) {
    if (activeOutputs && !activeOutputs[request.outputIndex])
      continue;
    anyActive = true;
    const auto &output = outputs[request.outputIndex];
    if (!output.data || output.elementCount != request.positionCount)
      return {WVKernelStatusCode::invalidShape,
              "A sampled moving-field output has the wrong shape."};
    for (std::size_t local = 0; local < request.positionCount; ++local) {
      const auto position = request.positionOffset + local;
      if (!std::isfinite(positions.x[position]) ||
          !std::isfinite(positions.y[position]) ||
          (request.naturalRank == WVPortableNaturalRank::volume &&
           (!positions.z || !std::isfinite(positions.z[position]))))
        return invalid("Moving coordinates must be finite.");
    }
  }
  if (!anyActive)
    return WVKernelStatus::ok();
  auto &workspaceMetrics = mutableMetrics();
  const auto outerWorkspaceBytes =
      workspaceMetrics.diagnosticWorkspaceLiveBytes;
  struct ResetSampledMovingWorkspace {
    WVFieldEvaluationMetrics &metrics;
    std::size_t outerBytes;
    ~ResetSampledMovingWorkspace() {
      metrics.diagnosticWorkspaceLiveBytes = outerBytes;
    }
  } resetWorkspace{workspaceMetrics, outerWorkspaceBytes};
  try {
    if(!sampledMovingWorkspace_ ||
        sampledMovingWorkspace_->fullStorage.size()<outputCount)
      return invalid("Sampled moving-field storage was not prepared.");
    auto& fullStorage=sampledMovingWorkspace_->fullStorage;
    auto& fullViews=sampledMovingWorkspace_->fullViews;
    auto& selection=sampledMovingWorkspace_->selection;
    auto& sampledStorage=sampledMovingWorkspace_->sampledStorage;
    for(std::size_t index=0;index<outputCount;++index) {
      fullStorage[index].clear();
      sampledStorage[index].clear();
      fullViews[index]={};
      selection[index]=0;
    }
    std::size_t transientSamplerBytes = 0;
    const auto accountWorkspace = [&]() {
      const std::size_t bytes=transientSamplerBytes;
      workspaceMetrics.diagnosticWorkspaceLiveBytes =
          outerWorkspaceBytes + bytes;
      workspaceMetrics.diagnosticWorkspaceHighWaterBytes = std::max(
          workspaceMetrics.diagnosticWorkspaceHighWaterBytes,
          workspaceMetrics.diagnosticWorkspaceLiveBytes);
      workspaceMetrics.additionalTransientHighWaterBytes=std::max(
          workspaceMetrics.additionalTransientHighWaterBytes,
          transientSamplerBytes);
    };
    for (const auto &request : implementation->requests) {
      if (activeOutputs && !activeOutputs[request.outputIndex])
        continue;
      selection[request.outputIndex] = 1;
      const auto elements = implementation->fullGridPlan.outputs()
                                [request.outputIndex].elementCount;
      fullStorage[request.outputIndex].resize(elements);
      fullViews[request.outputIndex] = {fullStorage[request.outputIndex].data(),
                                        elements};
    }
    accountWorkspace();
    auto status = evaluate(implementation->fullGridPlan, state,
                           fullViews.data(), fullViews.size(), selection.data());
    if (!status)
      return status;
    for (const auto &request : implementation->requests) {
      if (activeOutputs && !activeOutputs[request.outputIndex])
        continue;
      WVFieldSamplingRequest sampling;
      sampling.kind = WVFieldSamplingKind::positions;
      sampling.interpolation = request.interpolation;
      sampling.x.assign(positions.x + request.positionOffset,
                        positions.x + request.positionOffset +
                            request.positionCount);
      sampling.y.assign(positions.y + request.positionOffset,
                        positions.y + request.positionOffset +
                            request.positionCount);
      if (request.naturalRank == WVPortableNaturalRank::volume)
        sampling.z.assign(positions.z + request.positionOffset,
                          positions.z + request.positionOffset +
                              request.positionCount);
      const char *proxy = request.naturalRank == WVPortableNaturalRank::volume
                              ? "u"
                              : (barotropicQG_ ? "qgpv" : "ssu");
      WVFieldEvaluationPlan sampler;
      status = createPlanImpl(
          {{"moving-sampler", proxy, std::move(sampling)}},sampler,{},false);
      if (!status)
        return status;
      transientSamplerBytes = sampler.persistentBytes();
      auto &sampled = sampledStorage[request.outputIndex];
      sampled.resize(request.positionCount);
      accountWorkspace();
      status = samplePreparedField(
          sampler, fullStorage[request.outputIndex].data(),
          {sampled.data(), sampled.size()});
      if (!status)
        return status;
      transientSamplerBytes = 0;
    }
    for (const auto &request : implementation->requests) {
      if (activeOutputs && !activeOutputs[request.outputIndex])
        continue;
      std::copy(sampledStorage[request.outputIndex].begin(),
                sampledStorage[request.outputIndex].end(),
                outputs[request.outputIndex].data);
    }
    if (barotropicQG_)
      barotropicQG_->recordSampledMoving(positions.positionCount);
    else if (stratified_)
      stratified_->recordSampledMoving(positions.positionCount);
    else {
      ++metrics_.movingEvaluationCount;
      metrics_.movingPositionCount += positions.positionCount;
    }
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate sampled moving-field workspace."};
  }
}

WVRealFieldBundleView WVFieldEvaluationService::advectionFieldStorage() noexcept {
  if (barotropicQG_ || stratified_)
    return {};
  const auto spatial = transform_->descriptor().spatialShape();
  return {realScratch_.data(),{spatial.first,spatial.second,spatial.third,3}};
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingImpl(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    const WVRealFieldBundleConstView *preparedAdvectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  if(eventWorkspace_) {
    const auto status=eventWorkspace_->validateState({state});
    if(!status) return status;
  }
  if (!transform_) return {WVKernelStatusCode::unsupportedOperation,"This transform requires coefficient-family state views."};
  if (!sameTransformConfiguration(
          plan.configuration_, transform_->descriptor().configuration()))
    return invalid(
        "The moving-field plan belongs to a different transform configuration.");
  if (positions.positionCount != plan.positionCount_ ||
      (positions.positionCount != 0 &&
       (positions.x == nullptr || positions.y == nullptr ||
        positions.z == nullptr)))
    return {WVKernelStatusCode::invalidShape,
            "Moving coordinates must match the plan's shared position count."};
  if (outputCount != plan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Moving-field outputs must match the plan."};
  const auto spectral = transform_->descriptor().spectralShape();
  for (const auto view : {state.coefficients.Ap, state.coefficients.Am,
                          state.coefficients.A0})
    if (view.data == nullptr || view.shape.rows != spectral.rows ||
        view.shape.columns != spectral.columns)
      return {WVKernelStatusCode::invalidShape,
              "Moving-field coefficients must have shape [Nj,Nkl]."};
  for (std::size_t index = 0; index < outputCount; ++index)
    if ((!activeOutputs || activeOutputs[index]) && (outputs[index].data == nullptr ||
        outputs[index].elementCount != plan.outputs_[index].elementCount))
      return {WVKernelStatusCode::invalidShape,
              "Moving-field output shape does not match its request."};
  const auto finitePosition=[&](std::size_t index) {return std::isfinite(positions.x[index]) && std::isfinite(positions.y[index]) && std::isfinite(positions.z[index]);};
  if(!activeOutputs) {
    for(std::size_t index=0;index<positions.positionCount;++index)
      if(!finitePosition(index)) return invalid("Moving coordinates must be finite.");
  } else {
    bool anyActive=false;
    for(const auto& request:plan.requests_) if(activeOutputs[request.outputIndex]) {
      anyActive=true;
      for(std::size_t index=request.positionOffset;index<request.positionOffset+request.positionCount;++index)
        if(!finitePosition(index)) return invalid("Moving coordinates must be finite.");
    }
    if(!anyActive) return WVKernelStatus::ok();
  }
  ExecutionGuard guard(executing_);
  if (!guard.entered())
    return {WVKernelStatusCode::reentrantExecution,
            "Field evaluation is not reentrant."};

  const auto &configuration = transform_->descriptor().configuration();
  const auto spatial = transform_->descriptor().spatialShape();
  const auto R = spatial.elementCount();
  const auto horizontalCount = configuration.Nx * configuration.Ny;
  const double *primitiveFields = realScratch_.data();
  bool primitiveReused=false;
  if (preparedAdvectionFields == nullptr) {
    WVRealFieldBundleView fields{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    const auto before = transform_->metrics().executionCount;
    const auto operation=[&](){return transform_->transformWaveVortexToUVWEta(state, fields);};
    const WVVariableEvaluationKey primitiveKey{
        WVVariableEvaluationNode::physicalField,
        static_cast<std::uint32_t>(WVPortableVariable::u),
        eventWorkspace_ ? eventWorkspace_->component() : 0u,0,0,0,1};
    const auto status = eventWorkspace_ ? eventWorkspace_->evaluate(primitiveKey,fields.data,4*R,operation,primitiveReused) : operation();
    if (!status)
      return status;
    metrics_.fftExecutionCount += transform_->metrics().executionCount - before;
    if(primitiveReused) ++metrics_.primitiveFieldReuseCount;
    else ++metrics_.movingPrimitiveTransformCount;
  } else {
    if (preparedAdvectionFields->data == nullptr ||
        preparedAdvectionFields->shape.first != configuration.Nx ||
        preparedAdvectionFields->shape.second != configuration.Ny ||
        preparedAdvectionFields->shape.third != configuration.Nz ||
        preparedAdvectionFields->shape.fourth != 3)
      return {WVKernelStatusCode::invalidShape,
              "Prepared advection fields must have shape [Nx,Ny,Nz,3]."};
    if (std::any_of(plan.requests_.begin(),plan.requests_.end(),[&](const auto &request) {
          return (!activeOutputs || activeOutputs[request.outputIndex]) && request.primitiveChannel > 2;
        }))
      return {WVKernelStatusCode::unsupportedOperation,
              "Prepared advection fields support only u, v, and w requests."};
    primitiveFields = preparedAdvectionFields->data;
    ++metrics_.primitiveFieldReuseCount;
  }
  ++metrics_.evaluationCount;
  ++metrics_.movingEvaluationCount;
  metrics_.movingPositionCount += positions.positionCount;
  if(!primitiveReused) {
    ++metrics_.transformCount;
    metrics_.primitiveFieldEvaluationCount += preparedAdvectionFields == nullptr ? 4 : 3;
  }
  metrics_.scratchHighWaterBytes =
      std::max(metrics_.scratchHighWaterBytes, 4 * R * sizeof(double));

  const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
  const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
  const double dz = configuration.Lz / static_cast<double>(configuration.Nz - 1);
  const double densityScale =
      configuration.rho0 * configuration.N0 * configuration.N0 /
      configuration.g;
  const bool needsTotalDensity = std::any_of(
      plan.requests_.begin(), plan.requests_.end(), [&](const auto &request) {
        return (!activeOutputs || activeOutputs[request.outputIndex]) && request.primitiveChannel == 5;
      });
  double *totalDensity = nullptr;
  if (needsTotalDensity) {
    totalDensity = realScratch_.data() + 4 * R;
    const auto &vertical = transform_->descriptor().verticalModes().z;
    const double *eta = primitiveFields + 3 * R;
    for (std::size_t iz = 0; iz < configuration.Nz; ++iz) {
      const double densityNoMotion =
          configuration.rho0 - densityScale * vertical[iz];
      for (std::size_t horizontal = 0; horizontal < horizontalCount;
           ++horizontal) {
        const auto index = horizontal + horizontalCount * iz;
        totalDensity[index] = densityNoMotion + densityScale * eta[index];
      }
    }
    metrics_.scratchHighWaterBytes =
        std::max(metrics_.scratchHighWaterBytes, 5 * R * sizeof(double));
  }
  for (const auto &request : plan.requests_) {
    if(activeOutputs && !activeOutputs[request.outputIndex]) continue;
    auto &output = outputs[request.outputIndex];
    for (std::size_t local = 0; local < request.positionCount; ++local) {
      const auto position = request.positionOffset + local;
      const double x = wrapped(positions.x[position], configuration.Lx);
      const double y = wrapped(positions.y[position], configuration.Ly);
      const double z = positions.z[position];
      double value = 0.0;
      if (z >= -configuration.Lz && z <= 0.0) {
        const auto channel = std::min<std::size_t>(request.primitiveChannel, 3);
        const double *source = request.primitiveChannel == 5
                                   ? totalDensity
                                   : primitiveFields + channel * R;
        if (request.interpolation == WVPositionInterpolation::linear) {
          const auto x0 = std::min(
              static_cast<std::size_t>(std::floor(x / dx)),
              configuration.Nx - 1);
          const auto y0 = std::min(
              static_cast<std::size_t>(std::floor(y / dy)),
              configuration.Ny - 1);
          const double normalizedZ = (z + configuration.Lz) / dz;
          const auto z0 = std::min(
              static_cast<std::size_t>(
                  std::max(0.0, std::floor(normalizedZ))),
              configuration.Nz - 2);
          const double wx1 = (x - static_cast<double>(x0) * dx) / dx;
          const double wy1 = (y - static_cast<double>(y0) * dy) / dy;
          const double wz1 = std::clamp(normalizedZ - static_cast<double>(z0),
                                        0.0, 1.0);
          const std::array<std::size_t, 2> xi{{x0,
                                               (x0 + 1) % configuration.Nx}};
          const std::array<std::size_t, 2> yi{{y0,
                                               (y0 + 1) % configuration.Ny}};
          const std::array<double, 2> wx{{1.0 - wx1, wx1}};
          const std::array<double, 2> wy{{1.0 - wy1, wy1}};
          const std::array<double, 2> wz{{1.0 - wz1, wz1}};
          for (std::size_t iz = 0; iz < 2; ++iz)
            for (std::size_t iy = 0; iy < 2; ++iy)
              for (std::size_t ix = 0; ix < 2; ++ix)
                value += source[xi[ix] + configuration.Nx * yi[iy] +
                                horizontalCount * (z0 + iz)] *
                         wx[ix] * wy[iy] * wz[iz];
          ++metrics_.linearInterpolationCount;
        } else {
          const auto xLower = std::min(
              static_cast<std::size_t>(std::floor(x / dx)),
              configuration.Nx - 1);
          const auto yLower = std::min(
              static_cast<std::size_t>(std::floor(y / dy)),
              configuration.Ny - 1);
          const bool xBoundary =
              xLower < 3 || xLower > configuration.Nx - 4;
          const bool yBoundary =
              yLower < 3 || yLower > configuration.Ny - 4;
          const std::size_t xShift = xBoundary ? 4 : 0;
          const std::size_t yShift = yBoundary ? 4 : 0;
          const double xQuery =
              xBoundary ? wrapped(x + 4.0 * dx, configuration.Lx) : x;
          const double yQuery =
              yBoundary ? wrapped(y + 4.0 * dy, configuration.Ly) : y;
          if (xQuery <= static_cast<double>(configuration.Nx - 1) * dx &&
              yQuery <= static_cast<double>(configuration.Ny - 1) * dy) {
            movingWorkspace_->xSpline.weightsInto(
                0.0, dx, xQuery, movingWorkspace_->xWeights, xShift,
                &movingWorkspace_->xRightHandSide,
                &movingWorkspace_->xShifted);
            movingWorkspace_->ySpline.weightsInto(
                0.0, dy, yQuery, movingWorkspace_->yWeights, yShift,
                &movingWorkspace_->yRightHandSide,
                &movingWorkspace_->yShifted);
            movingWorkspace_->zSpline.weightsInto(
                -configuration.Lz, dz, z, movingWorkspace_->zWeights, 0,
                &movingWorkspace_->zRightHandSide,
                &movingWorkspace_->zShifted);
            for (std::size_t iz = 0; iz < configuration.Nz; ++iz)
              for (std::size_t iy = 0; iy < configuration.Ny; ++iy)
                for (std::size_t ix = 0; ix < configuration.Nx; ++ix)
                  value += source[ix + configuration.Nx * iy +
                                  horizontalCount * iz] *
                           movingWorkspace_->xWeights[ix] *
                           movingWorkspace_->yWeights[iy] *
                           movingWorkspace_->zWeights[iz];
          }
          ++metrics_.splineInterpolationCount;
        }
        if (request.primitiveChannel == 4)
          value *= densityScale;
      } else if (request.interpolation == WVPositionInterpolation::linear) {
        ++metrics_.linearInterpolationCount;
      } else {
        ++metrics_.splineInterpolationCount;
      }
      output.data[local] = value;
    }
    metrics_.outputElementWriteCount += output.elementCount;
  }
  return WVKernelStatus::ok();
}

const WVStratifiedModalGeometry* WVFieldEvaluationService::stratifiedGeometry() const noexcept {
  return stratified_ ? &stratified_->configuration() : nullptr;
}

const WVTransformConstantStratificationConfiguration &
WVFieldEvaluationService::configuration() const noexcept {
  return transform_->descriptor().configuration();
}

WVKernelStatus WVFieldEvaluationService::createStateLayout(
    const WVPortableObserverDescriptor &descriptor,
    WVIntegrationStateLayout &layout) const {
  (void)descriptor;
  if (stratified_) {
    const auto& g=stratified_->configuration();
    WVTransformStateDescription description{g.transformClass,{g.Nx,g.Ny,g.Nz},{},true};
    if ((g.transformClass=="WVTransformHydrostatic" || g.transformClass=="WVTransformBoussinesq")) for (const char* name:{"Ap","Am"}) description.coefficientFamilies.push_back({name,{g.Nj,g.Nkl},WVToleranceKind::coefficientEnergyScaled});
    description.coefficientFamilies.push_back({"A0",{g.Nj,g.Nkl},WVToleranceKind::coefficientEnergyScaled});
    return WVIntegrationStateLayout::createCoefficientOnly(std::move(description),layout);
  }
  if (barotropicQG_) {
    const auto &configuration = barotropicQG_->configuration();
    WVTransformBarotropicQGDescriptor transform;
    auto status = WVTransformBarotropicQGDescriptor::create(configuration,
                                                             transform);
    if (!status)
      return status;
    return WVIntegrationStateLayout::createCoefficientOnly(
        {"WVTransformBarotropicQG",
         {configuration.Nx, configuration.Ny},
         {{"A0", {transform.Nkl()},
           WVToleranceKind::coefficientEnergyScaled}},
         true},
        layout);
  }
  const auto &configuration = transform_->descriptor().configuration();
  const auto shape = transform_->descriptor().spectralShape();
  return WVIntegrationStateLayout::createCoefficientOnly(
      {"WVTransformConstantStratification",
       {configuration.Nx, configuration.Ny, configuration.Nz},
       {{"Ap", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled},
        {"Am", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled},
        {"A0", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled}},
       true},
      layout);
}

bool WVFieldEvaluationService::isCompatibleWith(
    const WVIntegrationStateLayout &layout) const noexcept {
  if (barotropicQG_)
    return barotropicQG_->isCompatibleWith(layout);
  if (stratified_)
    return stratified_->isCompatibleWith(layout);
  if (!transform_ || !layout.hasLegacyCoefficientTriple())
    return false;
  const auto shape = transform_->descriptor().spectralShape();
  const auto &configuration = transform_->descriptor().configuration();
  return layout.coefficientShape().rows == shape.rows &&
         layout.coefficientShape().columns == shape.columns &&
         layout.spatialDimensions() ==
             std::vector<std::size_t>{configuration.Nx, configuration.Ny,
                                      configuration.Nz};
}

bool WVFieldEvaluationService::isCompatibleWith(
    const WVFieldEvaluationService &other) const noexcept {
  if(static_cast<bool>(stratified_)!=static_cast<bool>(other.stratified_)) return false;
  if(stratified_) {
    const auto& a=stratified_->configuration(); const auto& b=other.stratified_->configuration();
    return a.Nx==b.Nx && a.Ny==b.Ny && a.Nj==b.Nj && a.Nkl==b.Nkl && a.Lx==b.Lx && a.Ly==b.Ly && a.z==b.z && a.j==b.j && a.k==b.k && a.l==b.l;
  }

  if (static_cast<bool>(barotropicQG_) !=
      static_cast<bool>(other.barotropicQG_))
    return false;
  if (barotropicQG_)
    return sameTransformConfiguration(barotropicQG_->configuration(),
                                      other.barotropicQG_->configuration());
  return transform_ != nullptr && other.transform_ != nullptr &&
         sameTransformConfiguration(transform_->descriptor().configuration(),
                                    other.transform_->descriptor()
                                        .configuration());
}

const WVFieldEvaluationMetrics &
WVFieldEvaluationService::metrics() const noexcept {
  if (stratified_) metrics_ = stratified_->metrics();
  else if (barotropicQG_) metrics_ = barotropicQG_->metrics();
  metrics_.servicePersistentBytes = persistentBytes();
  return metrics_;
}

WVVariableEvaluationMetrics
WVFieldEvaluationService::activeVariableEvaluationMetrics() const noexcept {
  return eventWorkspace_ ? eventWorkspace_->evaluationMetrics() :
      WVVariableEvaluationMetrics{};
}

bool WVFieldEvaluationService::activeDensityRecoveryReport(
    WVNoMotionRecoveryReport& report) const noexcept {
  return eventWorkspace_ && eventWorkspace_->densityRecoveryReport(report);
}

WVVariableProducerMetrics
WVFieldEvaluationService::producerMetrics() const noexcept {
  if(stratified_) return stratified_->producerMetrics();
  if(barotropicQG_) return barotropicQG_->producerMetrics();
  WVVariableProducerMetrics result;
  const auto& metrics=transform_->metrics();
  result.stateValidations=metrics.stateValidationCount;
  result.phasePreparations=metrics.phasePreparationCount;
  result.derivedValidations=metrics.derivedValidationCount;
  result.tendencyReconstructions=metrics.tendencyReconstructionCount;
  result.reconstructions=metrics.reconstructionCount;
  result.horizontalSpeedReductions=
      outputProducerMetrics_.horizontalSpeedReductions;
  result.verticalSpeedReductions=outputProducerMetrics_.verticalSpeedReductions;
  result.energyReductions=outputProducerMetrics_.energyReductions;
  for(std::size_t field=0;field<result.reconstructions.size();++field)
    for(std::size_t derivative=0;
        derivative<result.reconstructions[field].size();++derivative)
      for(std::size_t component=0;
          component<result.reconstructions[field][derivative].size();
          ++component)
        result.reconstructions[field][derivative][component]+=
            outputProducerMetrics_.reconstructions[field][derivative][component];
  return result;
}

WVFieldEvaluationMetrics &WVFieldEvaluationService::mutableMetrics() noexcept {
  if (stratified_)
    return stratified_->mutableMetrics();
  if (barotropicQG_)
    return barotropicQG_->mutableMetrics();
  return metrics_;
}

std::size_t WVFieldEvaluationService::persistentBytes() const noexcept {
  const auto forcingBytes=forcing_ ? forcing_->persistentBytes() : 0;
  const auto evaluationBytes=variableEvaluationContext_.persistentBytes();
  const auto arenaBytes=eventArena_ ? eventArena_->persistentBytes() : 0;
  const auto sampledMovingBytes=sampledMovingWorkspace_ ?
      sampledMovingWorkspace_->persistentBytes() : 0;
  const auto builtinBytes=builtinNonlinearFieldsPlan_.persistentBytes()+
      builtinNonlinearPhysical_.capacity()*sizeof(double)+
      builtinNonlinearCoefficients_.capacity()*sizeof(WVComplex64);
  if(stratified_) return sizeof(*this)+stratified_->persistentBytes()+forcingBytes+evaluationBytes+arenaBytes+sampledMovingBytes+builtinBytes;
  if (barotropicQG_)
    return sizeof(*this) + barotropicQG_->persistentBytes()+forcingBytes+evaluationBytes+arenaBytes+sampledMovingBytes+builtinBytes;
  return sizeof(*this) + forcingBytes + evaluationBytes + arenaBytes +
         sampledMovingBytes+builtinBytes+
         (ownedTransform_ ? transform_->persistentBytes() : 0) +
         realScratch_.capacity() * sizeof(double) +
         complexScratch_.capacity() * sizeof(WVComplex64) +
         eventBatchInvocations_.capacity() * sizeof(PlanInvocation) +
         (movingWorkspace_ ? movingWorkspace_->persistentBytes() : 0);
}

} // namespace wavevortex::runtime
