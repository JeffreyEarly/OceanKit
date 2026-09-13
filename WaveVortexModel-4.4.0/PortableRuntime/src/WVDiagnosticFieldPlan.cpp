#include "WVDiagnosticFieldPlan.hpp"
#include "WVForcingDiagnosticBinding.hpp"
#include "WVBarotropicQGFieldEvaluationAdapter.hpp"
#include "WVStratifiedFieldEvaluationAdapter.hpp"
#include "WVFieldEvaluationEventWorkspace.hpp"
#include "WVDensityEventEvaluation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace wavevortex::runtime::detail {
namespace {
using Variable = WVPortableVariable;
WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration,std::move(message)};
}
WVKernelStatus unsupported(std::string message) {
  return {WVKernelStatusCode::unsupportedOperation,std::move(message)};
}
bool sameSampling(const WVFieldSamplingRequest& a,const WVFieldSamplingRequest& b) {
  return a.kind==b.kind && a.interpolation==b.interpolation && a.xIndices==b.xIndices &&
      a.yIndices==b.yIndices && a.x==b.x && a.y==b.y && a.z==b.z;
}
std::size_t samplingBit(WVFieldSamplingKind kind) {
  switch(kind) {
    case WVFieldSamplingKind::fullGrid: return portableFullGridSampling;
    case WVFieldSamplingKind::fixedVerticalProfiles: return portableFixedVerticalProfileSampling;
    case WVFieldSamplingKind::positions: return portablePositionSampling;
  }
  return 0;
}
bool overlap(const void* a,std::size_t n,const void* b,std::size_t m) {
  if (!a || !b || !n || !m) return false;
  const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
  return x<=y ? y-x<n : x-y<m;
}
template<class Visitor>
void visitPortableExecution(const WVPortableVariablePlan& plan,Visitor&& visitor) {
  for(std::size_t index=0;index<plan.count;++index)
    visitor(plan.order[index],index+1==plan.count);
}
WVVariableEvaluationKey scratchKey(std::uint64_t signature,std::size_t group,
    std::size_t slot,std::uint32_t stage) {
  return {WVVariableEvaluationNode::registeredVariable,
      static_cast<std::uint32_t>(slot),static_cast<std::uint32_t>(group),
      0,0,0,stage,signature | (1ULL<<63)};
}
}

bool WVDiagnosticFieldPlan::required(const std::vector<WVFieldRequest>& requests,bool stratified) noexcept {
  for(const auto& request:requests) {
    if(isPortableForcingVariableName(request.fieldName)) return true;
    const auto* metadata=findPortableVariable(request.fieldName);
    if(metadata && (metadata->ordinal>=23 ||
        metadata->identifier==Variable::uvMax ||
        metadata->identifier==Variable::wMax ||
        (stratified && metadata->identifier==Variable::rhoBar))) return true;
  }
  return false;
}

WVKernelStatus WVDiagnosticFieldPlan::configure(const WVFieldEvaluationService& service) {
  owner_=&service;
  constant_=service.transform_;
  bool antialias=false;
  if(constant_) {
    const auto& c=constant_->descriptor().configuration();
    spatial_=constant_->descriptor().spatialShape(); spectral_=constant_->descriptor().spectralShape();
    isHydrostatic_=c.isHydrostatic; antialias=c.shouldAntialias; Lz_=c.Lz; constantN2_=c.N0*c.N0;
    configuration_=isHydrostatic_ ? "constant-hydrostatic" : "constant-nonhydrostatic";
  } else if(service.stratified_) {
    modal_=&service.stratified_->configuration();
    hydrostatic_=service.stratified_->hydrostaticKernel_;
    boussinesq_=service.stratified_->boussinesqKernel_;
    isQG_=!hydrostatic_ && !boussinesq_; isHydrostatic_=!boussinesq_;
    spatial_={modal_->Nx,modal_->Ny,modal_->Nz}; spectral_={modal_->Nj,modal_->Nkl};
    Lz_=modal_->Lz; antialias=modal_->shouldAntialias;
    configuration_=isQG_ ? "stratified-qg" : hydrostatic_ ? "hydrostatic" : "boussinesq";
  } else if(service.barotropicQG_) {
    const auto& c=service.barotropicQG_->configuration();
    spatial_={c.Nx,c.Ny,1};
    WVIntegrationStateLayout layout;
    auto status=service.createStateLayout({},layout); if(!status) return status;
    spectral_={1,layout.coefficientElementCount()};
    isBarotropic_=isQG_=true; Lz_=c.h; barotropicG_=c.g; antialias=c.shouldAntialias;
    configuration_="barotropic";
  } else return invalid("Diagnostic evaluation has no resolved transform.");
  configuration_+=antialias ? "-aa1" : "-aa0";
  return WVKernelStatus::ok();
}

std::string WVDiagnosticFieldPlan::configurationIdentifier(const WVFieldEvaluationService& service) {
  WVDiagnosticFieldPlan plan;
  return plan.configure(service) ? plan.configuration_ : std::string{};
}

WVKernelStatus WVDiagnosticFieldPlan::create(const WVFieldEvaluationService& service,
    const std::vector<WVFieldRequest>& requests,WVFieldEvaluationPlan& result,
    WVDensityDiagnosticContract contract,bool prepareScientificDependencies) {
  return createImpl(service,requests,result,false,contract,
      prepareScientificDependencies);
}

WVKernelStatus WVDiagnosticFieldPlan::createDensityQualification(
    const WVFieldEvaluationService& service,const std::vector<WVFieldRequest>& requests,
    WVDensityDiagnosticContract contract,WVFieldEvaluationPlan& result) {
  if(requests.empty()) return invalid("Density qualification requires at least one output.");
  for(const auto& request:requests)
    if(request.fieldName!="rho_nm" && request.fieldName!="eta_true" && request.fieldName!="ape" && request.fieldName!="apv")
      return unsupported("Density qualification accepts only rho_nm, eta_true, ape and apv.");
  if(contract.reference!=WVNoMotionReference::actual && contract.reference!=WVNoMotionReference::initial)
    return invalid("Density reference selection is invalid.");
  return createImpl(service,requests,result,true,contract,true);
}

WVKernelStatus WVDiagnosticFieldPlan::createDensityQualificationForActiveEvaluation(
    const WVFieldEvaluationService& service,
    const std::vector<WVFieldRequest>& requests,
    WVDensityDiagnosticContract contract,WVFieldEvaluationPlan& result) {
  if(!service.eventWorkspace_ || !service.stateEvaluationActive_)
    return {WVKernelStatusCode::invalidConfiguration,
        "Active density qualification requires an evaluation session."};
  if(service.variableEvaluationPolicy_!=WVVariableEvaluationPolicy::reuse)
    return {WVKernelStatusCode::invalidConfiguration,
        "Active density qualification requires the reuse policy."};
  if(service.executing_ || !service.eventWorkspace_->preparationIdle())
    return {WVKernelStatusCode::reentrantExecution,
        "Active density qualification is allowed only between producer calls."};
  WVFieldEvaluationPlan candidate;
  auto status=createImpl(service,requests,candidate,true,contract,false);
  if(!status) return status;
  status=service.prepareEventArenaForActiveEvaluation(candidate);
  service.eventArena_->refreshAccounting();
  if(!status) return status;
  result=std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus WVDiagnosticFieldPlan::createImpl(const WVFieldEvaluationService& service,
    const std::vector<WVFieldRequest>& requests,WVFieldEvaluationPlan& result,
    bool densityQualification,WVDensityDiagnosticContract densityContract,
    bool prepareScientificDependencies) {
  if(densityContract.reference!=WVNoMotionReference::actual && densityContract.reference!=WVNoMotionReference::initial)
    return invalid("Density reference selection is invalid.");
  try {
    auto plan=std::shared_ptr<WVDiagnosticFieldPlan>(new WVDiagnosticFieldPlan);
    for(const auto& request:requests) {
      for(const auto character:request.fieldName) {
        plan->scratchSignature_^=static_cast<unsigned char>(character);
        plan->scratchSignature_*=1099511628211ULL;
      }
      plan->scratchSignature_^=static_cast<std::uint8_t>(request.sampling.kind);
      plan->scratchSignature_*=1099511628211ULL;
    }
    auto status=plan->configure(service); if(!status) return status;
    plan->densityQualification_=densityQualification;
    plan->densityContract_=densityContract;
    plan->hasDensity_=std::any_of(requests.begin(),requests.end(),[](const auto& request) {
      return request.fieldName=="rho_nm" || request.fieldName=="eta_true" || request.fieldName=="ape" || request.fieldName=="apv";
    });
    if(plan->hasDensity_) {
      if(plan->isQG_) return unsupported("Density evaluation requires a wave-bearing transform.");
      if(plan->modal_) {
        plan->densityHeights_=plan->modal_->z;
        plan->densityWeights_=plan->modal_->z_int;
        plan->densityInitial_=plan->modal_->rho_nm0;
        plan->densityGravity_=plan->modal_->g;
        plan->densityReference_=plan->modal_->rho0;
      } else {
        const auto& c=plan->constant_->descriptor().configuration();
        plan->densityHeights_=plan->constant_->descriptor().verticalModes().z;
        plan->densityWeights_.resize(c.Nz);
        plan->densityInitial_.resize(c.Nz);
        const double scale=c.rho0*c.N0*c.N0/c.g;
        for(std::size_t z=0;z<c.Nz;++z) {
          plan->densityWeights_[z]=plan->weight(z);
          plan->densityInitial_[z]=c.rho0-scale*plan->densityHeights_[z];
        }
        plan->densityGravity_=c.g; plan->densityReference_=c.rho0;
      }
    }
    std::set<std::string> identifiers;
    const auto dependency=[&](std::size_t group,std::string name,const WVFieldSamplingRequest& sampling) {
      auto& list=plan->groups_[group].requests;
      for(std::size_t index=0;index<list.size();++index)
        if(list[index].fieldName==name && sameSampling(list[index].sampling,sampling)) return index;
      const auto index=list.size();
      list.push_back({"diagnostic-"+std::to_string(group)+"-"+std::to_string(index),std::move(name),sampling});
      return index;
    };
    const auto prepareSampler = [&](Output &output,
                                    WVPortableNaturalRank rank,
                                    const WVFieldSamplingRequest &sampling) {
      if (sampling.kind == WVFieldSamplingKind::fullGrid)
        return WVKernelStatus::ok();
      if (rank != WVPortableNaturalRank::volume &&
          rank != WVPortableNaturalRank::horizontal)
        return invalid("Only gridded diagnostics support sampled output: " +
                       output.specification.fieldName);
      const char *proxy = rank == WVPortableNaturalRank::volume
                              ? "u"
                              : (plan->isBarotropic_ ? "qgpv" : "ssu");
      WVFieldEvaluationPlan sampler;
      auto samplerStatus = service.createPlanImpl(
          {{"diagnostic-sampler", proxy, sampling}}, sampler,{},
          prepareScientificDependencies);
      if (!samplerStatus)
        return samplerStatus;
      output.sampled = true;
      output.sampling = sampling;
      output.sampler = std::move(sampler);
      output.specification.dimensions = output.sampler.outputs().front().dimensions;
      output.specification.elementCount =
          output.sampler.outputs().front().elementCount;
      return WVKernelStatus::ok();
    };
    for(const auto& request:requests) {
      if(request.identifier.empty() || !identifiers.insert(request.identifier).second)
        return invalid("Diagnostic output identifiers must be nonempty and unique.");
      const auto* m=findPortableVariable(request.fieldName);
      const auto* contract=m ? portableVariableContract(m->identifier,plan->configuration_) : nullptr;
      if(!m || (contract && std::string_view(contract->authority)=="forcing-instance-template")) {
        if(!service.forcing_) return unsupported("Field service has no resolved forcing schedule: "+request.fieldName);
        WVForcingDiagnosticBinding::Output bound;
        auto binding=service.forcing_->resolve(request.fieldName,plan->configuration_,static_cast<std::uint8_t>(samplingBit(request.sampling.kind)),bound);
        if(!binding) return binding;
        plan->forcingPhysicalChannels_=std::max(plan->forcingPhysicalChannels_,
            bound.physicalChannels);
        const auto found=std::find(plan->forcingIndices_.begin(),plan->forcingIndices_.end(),bound.executionIndex);
        const auto slot=static_cast<std::size_t>(found-plan->forcingIndices_.begin());
        if(found==plan->forcingIndices_.end())
          plan->forcingIndices_.push_back(bound.executionIndex);
        Output output; output.variable=bound.contract->metadata.identifier;
        output.execution=bound.plan; output.forcing=true;
        output.forcingSlot=slot; output.forcingChannel=bound.channel; output.forcingPhysicalChannels=bound.physicalChannels;
        output.specification.identifier=request.identifier; output.specification.fieldName=request.fieldName;
        output.specification.samplingKind=request.sampling.kind;
        output.specification.dimensions=plan->isBarotropic_ ? std::vector<std::size_t>{plan->spatial_.first,plan->spatial_.second} :
            std::vector<std::size_t>{plan->spatial_.first,plan->spatial_.second,plan->spatial_.third};
        output.specification.elementCount=plan->spatial_.elementCount();
        auto samplerStatus=prepareSampler(output,bound.contract->metadata.naturalRank,request.sampling);
        if(!samplerStatus) return samplerStatus;
        plan->outputs_.push_back(std::move(output));
        continue;
      }
      // Constant-stratification rho_bar is a legacy natural-grid field.  It has
      // no generated diagnostic contract, but it remains a valid dependency
      // when a request batch also contains a canonical reduction.
      const bool legacyConstantRhoBar=!contract && plan->constant_ &&
          m->identifier==Variable::rhoBar;
      if(!contract && !legacyConstantRhoBar)
        return unsupported("Diagnostic is unavailable on this transform: "+request.fieldName);
      if(m->ordinal>=23 && (contract->metadata.samplingMask & (m->naturalRank==WVPortableNaturalRank::coefficient && request.sampling.kind==WVFieldSamplingKind::fullGrid ? static_cast<std::size_t>(portableCoefficientSampling) : samplingBit(request.sampling.kind)))==0)
        return unsupported("Diagnostic sampling contract is unsupported: "+request.fieldName);
      WVPortableVariablePlan resolved;
      WVPortableVariableOptions options; options.source=WVPortableOperationSource::builtIn;
      options.requireEvaluator=m->ordinal>=23 && !densityQualification;
      options.noMotionSolver=WVPortableNoMotionSolver::dampedLeastSquares;
      options.shouldUseTrueNoMotionProfile=densityContract.reference==WVNoMotionReference::actual;
      const auto sampling=m->naturalRank==WVPortableNaturalRank::coefficient ? static_cast<std::size_t>(portableCoefficientSampling) : samplingBit(request.sampling.kind);
      if(!legacyConstantRhoBar) {
        const auto resolution=resolvePortableVariablePlan(request.fieldName,plan->configuration_,
            static_cast<std::uint8_t>(sampling),options,resolved);
        if(resolution!=WVPortableVariableStatus::supported)
          return unsupported(std::string(contract->configurationRestriction)+": "+request.fieldName);
      }
      Output output; output.variable=m->identifier; output.execution=resolved;
      output.specification.identifier=request.identifier; output.specification.fieldName=request.fieldName;
      output.specification.samplingKind=request.sampling.kind;
      switch(m->identifier) {
        case Variable::rho_nm: case Variable::eta_true: case Variable::ape: case Variable::apv:
          output.density=true;
          output.dependency=dependency(0,"rho_total",{});
          plan->densityDependency_=output.dependency;
          if(m->identifier==Variable::apv) {
            output.auxiliaries[0]=dependency(0,"zeta_x",{});
            output.auxiliaries[1]=dependency(0,"zeta_y",{});
            output.auxiliaries[2]=dependency(0,"zeta_z",{});
          }
          output.specification.dimensions=m->identifier==Variable::rho_nm ?
              std::vector<std::size_t>{plan->spatial_.third} :
              std::vector<std::size_t>{plan->spatial_.first,plan->spatial_.second,plan->spatial_.third};
          output.specification.elementCount=m->identifier==Variable::rho_nm ? plan->spatial_.third : plan->spatial_.elementCount();
          break;
        case Variable::A0t: case Variable::Apt: case Variable::Amt:
        case Variable::phase: case Variable::conjPhase:
          output.specification.isComplex=true;
          output.specification.dimensions=plan->isBarotropic_ ? std::vector<std::size_t>{plan->spectral_.columns} :
              std::vector<std::size_t>{plan->spectral_.rows,plan->spectral_.columns};
          output.specification.elementCount=plan->spectral_.elementCount();
          break;
        case Variable::totalEnergySpatiallyIntegrated:
          output.auxiliaries[0]=dependency(0,"u",{});
          output.auxiliaries[1]=dependency(0,"v",{});
          output.auxiliaries[2]=dependency(0,"eta",{});
          if(!plan->isHydrostatic_) output.auxiliaries[3]=dependency(0,"w",{});
          output.specification.elementCount=1;
          break;
        default: {
          if(m->identifier==Variable::rhoBar && plan->modal_) {
            if(request.sampling.kind!=WVFieldSamplingKind::fullGrid) return unsupported("rho_bar requires its full vertical grid.");
            output.verticalMean=true; output.dependency=dependency(0,"rho_total",{});
          plan->densityDependency_=output.dependency; break;
          }
          std::string base=request.fieldName;
          if(m->identifier==Variable::totalEnergy) base="energy";
          else if(m->identifier==Variable::geostrophicEnergy) {
            base="energy"; output.group=plan->isQG_ ? 0 : 1;
          } else if(m->ordinal>=23) {
            const auto component=contract->component;
            if(component==WVPortableFlowComponent::total)
              return unsupported("No numerical diagnostic implementation: "+request.fieldName);
            output.group=static_cast<std::size_t>(component);
            base=base.substr(0,base.rfind('_'));
          }
          if((m->ordinal>=23 || request.sampling.kind==WVFieldSamplingKind::fullGrid) &&
              (base=="ssu" || base=="ssv" || base=="ssh")) {
            output.surface=true;
            base=base=="ssu" ? "u" : base=="ssv" ? "v" : "pi";
          } else if(request.sampling.kind==WVFieldSamplingKind::fullGrid &&
                    (base=="uvMax" || base=="wMax")) {
              output.extrema=true;
              output.auxiliaries[0]=dependency(0,base=="uvMax" ? "u" : "w",{});
              if(base=="uvMax") output.auxiliaries[1]=dependency(0,"v",{});
              output.specification.elementCount=1;
              break;
          }
          output.dependency=dependency(output.group,base,
              m->ordinal>=23 && request.sampling.kind!=WVFieldSamplingKind::fullGrid ?
                  WVFieldSamplingRequest{} : request.sampling);
          break;
        }
      }
      if(m->ordinal>=23) {
        auto samplerStatus=prepareSampler(output,m->naturalRank,request.sampling);
        if(!samplerStatus) return samplerStatus;
      }
      plan->outputs_.push_back(std::move(output));
    }
    const char* physicalNames[]={"u","v","w","eta"};
    for(std::size_t channel=0;channel<plan->forcingPhysicalChannels_;++channel)
      plan->forcingPhysicalDependencies_[channel]=dependency(0,physicalNames[channel],{});
    for(auto& group:plan->groups_) {
      if(group.requests.empty()) continue;
      status=service.createPlanImpl(group.requests,group.fields,{},
          prepareScientificDependencies); if(!status) return status;
    }
    WVFieldEvaluationPlan candidate;
    for(auto& output:plan->outputs_) {
      if(!output.sampled && !output.density && !output.forcing && !output.specification.isComplex && !output.extrema && output.variable!=Variable::totalEnergySpatiallyIntegrated) {
        const auto& primitive=plan->groups_[output.group].fields.outputs()[output.dependency];
        output.specification.dimensions=primitive.dimensions;
        output.specification.elementCount=primitive.elementCount;
        if(output.verticalMean) {
          output.specification.dimensions={plan->spatial_.third};
          output.specification.elementCount=plan->spatial_.third;
        }
        if(output.surface) {
          output.specification.dimensions={plan->spatial_.first,plan->spatial_.second};
          output.specification.elementCount=plan->spatial_.first*plan->spatial_.second;
        }
      }
      candidate.outputs_.push_back(output.specification);
    }
    if(prepareScientificDependencies) {
      status=plan->prepareEventArena(service);
      if(!status) return status;
    }
    candidate.diagnosticPlan_=std::move(plan);
    result=std::move(candidate);
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,"Unable to allocate the diagnostic dependency plan."};
  }
}

WVKernelStatus WVDiagnosticFieldPlan::prepareEventArena(
    const WVFieldEvaluationService& service) const {
  return prepareEventArenaImpl(service,false);
}

WVKernelStatus WVDiagnosticFieldPlan::prepareEventArenaForActiveEvaluation(
    const WVFieldEvaluationService& service) const {
  return prepareEventArenaImpl(service,true);
}

WVKernelStatus WVDiagnosticFieldPlan::prepareEventArenaImpl(
    const WVFieldEvaluationService& service,bool active) const {
  const auto R=spatial_.elementCount();
  const auto S=spectral_.elementCount();
  const auto preparePlan=[&](const WVFieldEvaluationPlan& plan,
      std::uint32_t component=0) {
    return active ? service.prepareEventArenaForActiveEvaluation(plan,component) :
        service.prepareEventArena(plan,component);
  };
  const auto prepareField=[&](const WVVariableEvaluationKey& key,
      std::size_t elements,bool complex) {
    return active ? service.prepareEventFieldForActiveEvaluation(
        key,elements,complex) : service.prepareEventField(key,elements,complex);
  };
  for(std::size_t group=0;group<groups_.size();++group) {
    if(groups_[group].fields.outputCount()) {
      const auto status=preparePlan(groups_[group].fields,
          static_cast<std::uint32_t>(group));
      if(!status) return status;
    }
    for(std::size_t output=0;output<groups_[group].fields.outputCount();++output) {
      const auto status=prepareField(scratchKey(
          scratchSignature_,group,output,16),
          groups_[group].fields.outputs()[output].elementCount,false);
      if(!status) return status;
    }
    if(group && !groups_[group].requests.empty())
      for(std::uint32_t family=0;family<3;++family) {
        const auto status=prepareField(
            {WVVariableEvaluationNode::componentCoefficients,family,
              static_cast<std::uint32_t>(group)},S,true);
        if(!status) return status;
      }
  }
  bool phase=false;
  std::uint8_t densityDemands=0;
  bool apvNeeded=false;
  for(const auto& output:outputs_) {
    if(active && output.sampled) {
      const auto status=preparePlan(output.sampler);
      if(!status) return status;
    }
    if(output.density) {
      densityDemands|=output.variable==Variable::rho_nm ?
          WVDensityEventEvaluation::rhoNmDemand :
          (output.variable==Variable::eta_true || output.variable==Variable::apv) ?
              WVDensityEventEvaluation::etaTrueDemand :
              WVDensityEventEvaluation::apeDemand;
      apvNeeded|=output.variable==Variable::apv;
    }
    visitPortableExecution(output.execution,[&](Variable node,bool) {
      phase|=node==Variable::Apt || node==Variable::Amt ||
          node==Variable::phase || node==Variable::conjPhase;
    });
    if(output.specification.isComplex) {
      const auto finalNode=output.execution.count ?
          output.execution.order[output.execution.count-1] : output.variable;
      const auto status=prepareField(
          {WVVariableEvaluationNode::registeredVariable,
            static_cast<std::uint32_t>(finalNode),
            static_cast<std::uint32_t>(output.group)},S,true);
      if(!status) return status;
    }
    if(output.variable==Variable::totalEnergySpatiallyIntegrated ||
        output.verticalMean || output.extrema) {
      const auto status=prepareField(
          {WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(output.variable)},
          output.specification.elementCount,false);
      if(!status) return status;
    }
  }
  if(phase) {
    const auto status=prepareField(
        {WVVariableEvaluationNode::phaseFactors},S,true);
    if(!status) return status;
  }
  if(densityDemands) {
    const auto status=active ? service.prepareDensityEventArenaForActiveEvaluation(
        R,densityHeights_.size(),densityDemands,densityContract_.reference,
        apvNeeded) : service.prepareDensityEventArena(R,densityHeights_.size(),
        densityDemands,densityContract_.reference,apvNeeded);
    if(!status) return status;
  }
  if(!forcingIndices_.empty()) {
    service.eventArenaForcingPrepared_=true;
    if(service.forcing_->horizontalMaximumNeeded()) {
      const auto status=prepareField(
          {WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(WVPortableVariable::uvMax)},1,false);
      if(!status) return status;
    }
    const auto channels=isQG_ ? 1u : isHydrostatic_ ? 3u : 4u;
    for(const auto index:forcingIndices_) {
      const auto status=prepareField(
          {WVVariableEvaluationNode::forcingTendency,
            static_cast<std::uint32_t>(isQG_ ?
              WVPortableVariable::Fqgpv_portable_catalog_forcing :
              WVPortableVariable::Fu_portable_catalog_forcing),0,0,0,
            static_cast<std::uint32_t>(index),1},channels*R,false);
      if(!status) return status;
    }
    for(std::size_t slot=0;slot<forcingIndices_.size();++slot) {
      const auto status=prepareField(scratchKey(
          scratchSignature_,0,slot,17),
          channels*R,false);
      if(!status) return status;
    }
    if(forcingPhysicalChannels_) {
      const auto status=prepareField(scratchKey(
          scratchSignature_,0,0,18),
          forcingPhysicalChannels_*R,false);
      if(!status) return status;
    }
  }
  if(std::any_of(outputs_.begin(),outputs_.end(),[](const auto& output) {
      return output.variable==Variable::apv;
    })) {
    const auto status=prepareField(scratchKey(
        scratchSignature_,0,0,19),3*R,false);
    if(!status) return status;
  }
  if(active) {
    auto status=service.eventArena_->prepareGroupNodesForActive(
        forcingIndices_.size());
    if(status && service.forcing_ && service.eventArenaForcingPrepared_)
      status=service.eventArena_->prepareForcingForActive(
          *service.forcing_,service.variableEvaluationPolicy_);
    return status;
  }
  return service.prepareEventArena(forcingIndices_.size());
}

WVKernelStatus WVDiagnosticFieldPlan::rebind(const WVFieldEvaluationService& service,WVFieldEvaluationPlan& result) const {
  try {
    std::vector<WVFieldRequest> requests;
    requests.reserve(outputs_.size());
    for(const auto& output:outputs_) {
      if(output.forcing) {
        if(!owner_->forcing_ || !service.forcing_) return unsupported("Rebinding forcing diagnostics requires a resolved forcing schedule.");
        if(configuration_!=service.portableVariableConfiguration())
          return invalid("Rebinding would change forcing diagnostic transform applicability.");
        WVForcingDiagnosticBinding::Output original,replacement;
        auto status=owner_->forcing_->resolve(output.specification.fieldName,configuration_,portableFullGridSampling,original);
        if(!status) return status;
        status=service.forcing_->resolve(output.specification.fieldName,configuration_,portableFullGridSampling,replacement);
        if(!status) return status;
        if(!owner_->forcing_->hasSamePrefix(*service.forcing_,original.executionIndex,replacement.executionIndex))
          return invalid("Rebinding would change the resolved forcing identity, stage, or ordered prefix.");
      }
      WVFieldSamplingRequest sampling;
      if(output.sampled)
        sampling=output.sampling;
      else if(!output.forcing && !output.specification.isComplex && !output.extrema && output.variable!=Variable::totalEnergySpatiallyIntegrated)
        sampling=groups_[output.group].requests[output.dependency].sampling;
      requests.push_back({output.specification.identifier,output.specification.fieldName,std::move(sampling)});
    }
    return densityQualification_ ? createDensityQualification(service,requests,densityContract_,result) : create(service,requests,result,densityContract_);
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,"Unable to rebind the diagnostic plan."};
  }
}

double WVDiagnosticFieldPlan::omega(std::size_t index) const noexcept {
  if(constant_) return constant_->descriptor().verticalModes().omega[index];
  if(hydrostatic_) return hydrostatic_->factors()[index].omega;
  if(boussinesq_) return boussinesq_->factors()[index].omega;
  return 0;
}
bool WVDiagnosticFieldPlan::keep(std::size_t group,std::size_t family,std::size_t index) const noexcept {
  const auto horizontal=index/spectral_.rows,vertical=index%spectral_.rows;
  const bool mean=constant_ ? constant_->descriptor().fourierModes()[horizontal].Kh==0 :
      modal_->modes[horizontal].k==0 && modal_->modes[horizontal].l==0;
  const bool zeroVertical=constant_ ? constant_->descriptor().verticalModes().j[vertical]==0 : modal_->j[vertical]==0;
  switch(group) {
    case 1: return family==2 && !mean;
    case 2: return family!=2 && !mean && !zeroVertical;
    case 3: return family!=2 && mean;
    case 4: return family==2 && mean && !zeroVertical;
    default: return true;
  }
}
double WVDiagnosticFieldPlan::weight(std::size_t z) const noexcept {
  if(isBarotropic_) return Lz_;
  if(modal_) return modal_->z_int[z];
  const double dz=Lz_/static_cast<double>(spatial_.third-1);
  return (z==0 || z+1==spatial_.third) ? dz/2 : dz;
}
double WVDiagnosticFieldPlan::stratification(std::size_t z) const noexcept {
  if(isBarotropic_) return barotropicG_/Lz_;
  return modal_ ? modal_->N2[z] : constantN2_;
}

WVKernelStatus WVDiagnosticFieldPlan::evaluate(WVFieldEvaluationService& service,const WVIntegrationState& input,
    WVFieldOutputView* outputs,std::size_t count,const std::uint8_t* activeOutputs) const {
  if(owner_!=&service) return invalid("The diagnostic plan belongs to a different field service.");
  const bool preparedState=service.eventWorkspace_ && service.stateEvaluationActive_;
  if(service.eventWorkspace_) {
    const auto stateStatus=service.eventWorkspace_->validateState(input);
    if(!stateStatus) return stateStatus;
  }
  if(count!=outputs_.size() || (count && !outputs))
    return {WVKernelStatusCode::invalidShape,
        "Diagnostic output-view count is invalid."};
  const auto active=[&](std::size_t index) {return !activeOutputs || activeOutputs[index];};
  WVState amplitudes=input.waveVortex;
  if(input.coefficientFamilyCount) {
    if(input.coefficientFamilyCount!=(isQG_ ? 1u : 3u) || !input.coefficientFamilies)
      return invalid("Diagnostic coefficient families do not match the transform.");
    for(std::size_t family=0;family<input.coefficientFamilyCount;++family) {
      const auto* layout=input.coefficientFamilies[family].layout;
      const char* expected=isQG_ ? "A0" : family==0 ? "Ap" : family==1 ? "Am" : "A0";
      if(!layout || layout->identifier!=expected || layout->elementCount!=spectral_.elementCount() ||
          (isBarotropic_ ? layout->spectralDimensions.size()!=1 || layout->spectralDimensions[0]!=spectral_.columns :
            layout->spectralDimensions.size()!=2 || layout->spectralDimensions[0]!=spectral_.rows || layout->spectralDimensions[1]!=spectral_.columns))
        return invalid("Diagnostic coefficient family layout does not match the transform.");
    }
    if(isQG_) amplitudes.coefficients.A0={input.coefficientFamilies[0].data,spectral_};
    else amplitudes.coefficients={{input.coefficientFamilies[0].data,spectral_},
        {input.coefficientFamilies[1].data,spectral_},{input.coefficientFamilies[2].data,spectral_}};
  }
  const WVComplexConstView coefficientViews[]={amplitudes.coefficients.Ap,amplitudes.coefficients.Am,amplitudes.coefficients.A0};
  for(std::size_t family=isQG_ ? 2 : 0;family<3;++family)
    if(coefficientViews[family].shape.rows!=spectral_.rows || coefficientViews[family].shape.columns!=spectral_.columns)
      return invalid("Diagnostic coefficient view has the wrong shape.");
  const WVComplex64* coefficients[]={amplitudes.coefficients.Ap.data,amplitudes.coefficients.Am.data,amplitudes.coefficients.A0.data};
  const auto n=spectral_.elementCount();
  if(!preparedState && (!std::isfinite(amplitudes.t) || !std::isfinite(amplitudes.t0) ||
      !std::isfinite(amplitudes.t-amplitudes.t0)))
    return invalid("Diagnostic state times and elapsed time must be finite.");
  if(input.additionalBlockCount && !input.additionalBlocks) return invalid("Diagnostic additional state views are missing.");
  for(std::size_t block=0;block<input.additionalBlockCount;++block)
    if(!input.additionalBlocks[block].layout) return invalid("Diagnostic additional state layout is missing.");
  if(!preparedState && !isQG_) for(std::size_t index=0;index<n;++index)
    if(!std::isfinite(omega(index)*(amplitudes.t-amplitudes.t0))) return invalid("Diagnostic phase would overflow.");
  for(std::size_t family=isQG_ ? 2 : 0;family<3;++family) {
    if(!coefficients[family]) return invalid("Diagnostic coefficients are missing.");
    if(!preparedState) for(std::size_t index=0;index<n;++index)
      if(!std::isfinite(coefficients[family][index].real) || !std::isfinite(coefficients[family][index].imag))
        return invalid("Diagnostic coefficients must be finite.");
  }
  for(std::size_t index=0;index<count;++index) {
    if(!active(index)) continue;
    const auto& spec=outputs_[index].specification;
    const void* pointer=spec.isComplex ? static_cast<void*>(outputs[index].complexData) : outputs[index].data;
    const auto bytes=spec.elementCount*(spec.isComplex ? sizeof(WVComplex64) : sizeof(double));
    if(outputs[index].elementCount!=spec.elementCount || (bytes && !pointer)) return invalid("Diagnostic output storage has the wrong shape or type.");
    for(const auto* coefficient:coefficients)
      if(overlap(pointer,bytes,coefficient,n*sizeof(WVComplex64)))
        return {WVKernelStatusCode::overlappingArrays,"Diagnostic output overlaps coefficient state."};
    for(std::size_t block=0;block<input.additionalBlockCount;++block) {
      const auto& view=input.additionalBlocks[block];
      if(overlap(pointer,bytes,view.realData,view.layout->elementCount*sizeof(double)) ||
          overlap(pointer,bytes,view.complexData,view.layout->elementCount*sizeof(WVComplex64)))
        return {WVKernelStatusCode::overlappingArrays,"Diagnostic output overlaps additional state."};
    }
    for(std::size_t previous=0;previous<index;++previous) {
      const auto& other=outputs_[previous].specification;
      if(overlap(pointer,bytes,other.isComplex ? static_cast<void*>(outputs[previous].complexData) : outputs[previous].data,
          other.elementCount*(other.isComplex ? sizeof(WVComplex64) : sizeof(double))))
        return {WVKernelStatusCode::overlappingArrays,"Diagnostic outputs overlap."};
    }
  }
  auto& metrics=service.stratified_ ? service.stratified_->metrics_ : service.barotropicQG_ ? service.barotropicQG_->metrics_ : service.metrics_;
  auto& producerMetrics=service.stratified_ ?
      service.stratified_->outputProducerMetrics_ : service.barotropicQG_ ?
      service.barotropicQG_->outputProducerMetrics_ :
      service.outputProducerMetrics_;
  const auto outputWritesBefore=metrics.outputElementWriteCount;
  std::size_t activeOutputElements=0;
  for(std::size_t index=0;index<count;++index)
    if(active(index)) activeOutputElements+=outputs_[index].specification.elementCount;
  const auto outerWorkspaceBytes=metrics.diagnosticWorkspaceLiveBytes;
  struct ResetLive {
    WVFieldEvaluationMetrics& metrics;
    std::size_t outerBytes;
    ~ResetLive() {metrics.diagnosticWorkspaceLiveBytes=outerBytes;}
  } reset{metrics,outerWorkspaceBytes};
  std::uint8_t densityDemands=0;
  bool needsAPV=false;
  for(std::size_t index=0;index<count;++index) if(active(index) && outputs_[index].density) {
    const auto variable=outputs_[index].variable;
    needsAPV|=variable==Variable::apv;
    densityDemands|=variable==Variable::rho_nm ? WVDensityEventEvaluation::rhoNmDemand :
        (variable==Variable::eta_true || variable==Variable::apv) ? WVDensityEventEvaluation::etaTrueDemand : WVDensityEventEvaluation::apeDemand;
  }
  WVFieldEvaluationEventScope densityScope(service,input,densityDemands && !service.eventWorkspace_,false);
  if(!densityScope.status()) return densityScope.status();
  if(densityDemands) {
    const auto status=service.eventWorkspace_->validateDensityBinding(amplitudes,densityContract_);
    if(!status) return status;
  }
  struct DensityUseGuard {
    WVFieldEvaluationEventWorkspace* workspace;
    WVDensityDiagnosticContract contract;
    ~DensityUseGuard() {
      if(workspace) workspace->finishDensityUse(contract);
    }
  } densityUseGuard{densityDemands ? service.eventWorkspace_ : nullptr,
      densityContract_};
  try {
    std::array<std::vector<std::uint8_t>,5> activeDependencies;
    for(std::size_t group=0;group<groups_.size();++group)
      activeDependencies[group].resize(groups_[group].fields.outputCount());
    std::vector<std::uint8_t> activeForcing(forcingIndices_.size());
    std::size_t physicalChannels=0;
    for(std::size_t index=0;index<count;++index) {
      if(!active(index)) continue;
      const auto& output=outputs_[index];
      if(output.density) {
        if(!service.eventWorkspace_->hasDensitySource()) activeDependencies[0][output.dependency]=1;
        if(output.variable==Variable::apv &&
            !service.eventWorkspace_->hasAPV(densityContract_.reference))
          for(std::size_t axis=0;axis<3;++axis) activeDependencies[0][output.auxiliaries[axis]]=1;
      } else if(output.forcing) {
        activeForcing[output.forcingSlot]=1;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::forcingTendency,
            static_cast<std::uint32_t>(isQG_ ?
                WVPortableVariable::Fqgpv_portable_catalog_forcing :
                WVPortableVariable::Fu_portable_catalog_forcing),0,0,0,
            static_cast<std::uint32_t>(forcingIndices_[output.forcingSlot]),1};
        if(!service.eventWorkspace_ || !service.eventWorkspace_->ready(key))
          physicalChannels=std::max(physicalChannels,output.forcingPhysicalChannels);
      } else if(!output.specification.isComplex) {
        if(output.variable==Variable::totalEnergySpatiallyIntegrated)
          for(std::size_t channel=0;channel<(isHydrostatic_ ? 3u : 4u);++channel)
            activeDependencies[0][output.auxiliaries[channel]]=1;
        else if(output.extrema) {
          activeDependencies[0][output.auxiliaries[0]]=1;
          if(output.variable==Variable::uvMax) activeDependencies[0][output.auxiliaries[1]]=1;
        } else activeDependencies[output.group][output.dependency]=1;
      }
    }
    for(std::size_t channel=0;channel<physicalChannels;++channel)
      activeDependencies[0][forcingPhysicalDependencies_[channel]]=1;
    std::array<std::vector<std::vector<double>>,5> fields;
    std::array<std::vector<WVComplex64>,3> masked;
    std::vector<WVComplex64> phases;
    std::vector<std::vector<double>> forcingFields;
    std::vector<WVForcingTendencyOutput> forcingViews;
    std::vector<double> forcingPhysical;
    std::vector<double> densityDerivatives;
    struct ScratchUse {WVVariableEvaluationKey key; std::vector<double>* storage;};
    std::vector<ScratchUse> scratchUses;
    std::size_t scratchUseCount=forcingIndices_.size()+
        (forcingPhysicalChannels_ ? 1u : 0u)+(needsAPV ? 1u : 0u);
    for(const auto& group:groups_) scratchUseCount+=group.fields.outputCount();
    scratchUses.reserve(scratchUseCount);
    struct ScratchGuard {
      WVFieldEvaluationEventWorkspace* workspace;
      std::vector<ScratchUse>& uses;
      ~ScratchGuard() {
        if(workspace) for(auto use=uses.rbegin();use!=uses.rend();++use)
          workspace->returnScratch(use->key,*use->storage);
      }
    } scratchGuard{service.eventWorkspace_,scratchUses};
    const auto checkout=[&](const WVVariableEvaluationKey& key,
        std::size_t elements,std::vector<double>& storage) {
      WVKernelStatus status=WVKernelStatus::ok();
      if(service.eventWorkspace_)
        status=service.eventWorkspace_->checkoutScratch(key,elements,storage);
      else try {storage.resize(elements);}
      catch(const std::bad_alloc&) {
        status={WVKernelStatusCode::allocationFailure,
            "Unable to allocate standalone diagnostic scratch storage."};
      }
      if(status) scratchUses.push_back({key,&storage});
      return status;
    };
    std::size_t primitiveCount=0;
    const auto account=[&](std::size_t viewBytes) {
      std::size_t bytes=viewBytes+activeForcing.capacity()*sizeof(std::uint8_t);
      std::size_t additional=bytes;
      bytes+=scratchUses.capacity()*sizeof(ScratchUse);
      additional+=scratchUses.capacity()*sizeof(ScratchUse);
      for(const auto& selection:activeDependencies) bytes+=selection.capacity()*sizeof(std::uint8_t);
      for(const auto& selection:activeDependencies)
        additional+=selection.capacity()*sizeof(std::uint8_t);
      for(const auto& group:fields) {
        bytes+=group.capacity()*sizeof(std::vector<double>);
        additional+=group.capacity()*sizeof(std::vector<double>);
        for(const auto& buffer:group) bytes+=buffer.capacity()*sizeof(double);
      }
      for(const auto& buffer:masked) bytes+=buffer.capacity()*sizeof(WVComplex64);
      bytes+=phases.capacity()*sizeof(WVComplex64);
      bytes+=(forcingPhysical.capacity()+densityDerivatives.capacity())*sizeof(double);
      bytes+=forcingFields.capacity()*sizeof(std::vector<double>)+
          forcingViews.capacity()*sizeof(WVForcingTendencyOutput);
      additional+=forcingFields.capacity()*sizeof(std::vector<double>)+
          forcingViews.capacity()*sizeof(WVForcingTendencyOutput);
      for(const auto& buffer:forcingFields) bytes+=buffer.capacity()*sizeof(double);
      if(!service.eventWorkspace_) {
        for(const auto& group:fields)
          for(const auto& buffer:group)
            additional+=buffer.capacity()*sizeof(double);
        for(const auto& buffer:masked)
          additional+=buffer.capacity()*sizeof(WVComplex64);
        additional+=phases.capacity()*sizeof(WVComplex64)+
            (forcingPhysical.capacity()+densityDerivatives.capacity())*
                sizeof(double);
        for(const auto& buffer:forcingFields)
          additional+=buffer.capacity()*sizeof(double);
      } else additional+=service.eventWorkspace_->externalWorkspaceBytes();
      metrics.diagnosticWorkspaceLiveBytes=outerWorkspaceBytes+bytes;
      metrics.diagnosticWorkspaceHighWaterBytes=std::max(
          metrics.diagnosticWorkspaceHighWaterBytes,
          metrics.diagnosticWorkspaceLiveBytes+
              metrics.densityWorkspaceLiveBytes);
      metrics.additionalTransientHighWaterBytes=std::max(
          metrics.additionalTransientHighWaterBytes,additional);
    };
    for(std::size_t group=0;group<groups_.size();++group) {
      const auto& plan=groups_[group].fields;
      if(std::none_of(activeDependencies[group].begin(),activeDependencies[group].end(),[](auto value){return value!=0;})) continue;
      std::vector<WVFieldOutputView> views(plan.outputCount());
      fields[group].resize(plan.outputCount());
      for(std::size_t output=0;output<plan.outputCount();++output) {
        if(!activeDependencies[group][output]) continue;
        auto& buffer=fields[group][output];
        const auto scratchStatus=checkout(scratchKey(
            scratchSignature_,group,output,16),
            plan.outputs()[output].elementCount,buffer);
        if(!scratchStatus) return scratchStatus;
        views[output]={buffer.data(),buffer.size()};
        ++primitiveCount;
      }
      WVIntegrationState selected=input;
      selected.waveVortex=amplitudes;
      std::array<WVCoefficientFamilyConstView,3> selectedFamilies{};
      struct ComponentUseGuard {
        WVFieldEvaluationService& service;
        WVFieldEvaluationEventWorkspace* workspace;
        WVIntegrationState& state;
        std::uint32_t component;
        std::size_t acquired=0;
        bool selected=false,registered=false,armed=true;
        void releasePins() noexcept {
          while(acquired) {
            --acquired;
            workspace->releaseComponentCoefficients(component,
                static_cast<std::uint32_t>(acquired));
          }
        }
        void fallback() noexcept {
          if(!armed || !workspace) return;
          if(selected) workspace->setComponent(0);
          // A failed unregister must keep the coefficient storage pinned. This
          // preserves the kernel's borrowed immutable view until event teardown.
          if(registered && workspace->policy()==WVVariableEvaluationPolicy::lowMemory &&
              !service.removeStateEvaluationView(state,component)) return;
          releasePins();
          armed=false;
        }
        WVKernelStatus release() {
          if(!armed || !workspace) return WVKernelStatus::ok();
          if(selected) workspace->setComponent(0);
          if(registered && workspace->policy()==WVVariableEvaluationPolicy::lowMemory) {
            const auto status=service.removeStateEvaluationView(state,component);
            if(!status) return status;
          }
          releasePins();
          armed=false;
          return WVKernelStatus::ok();
        }
        ~ComponentUseGuard() {fallback();}
      } componentUse{service,service.eventWorkspace_,selected,
          static_cast<std::uint32_t>(group)};
      if(group) {
        for(std::size_t family=0;family<3;++family) {
          if(service.eventWorkspace_) {
            const WVComplex64* prepared=nullptr;
            const auto componentStatus=service.eventWorkspace_->componentCoefficients(
                static_cast<std::uint32_t>(group),static_cast<std::uint32_t>(family),n,
                [&](WVComplex64* destination) {
                  for(std::size_t index=0;index<n;++index)
                    destination[index]=keep(group,family,index) ? coefficients[family][index] : WVComplex64{};
                },prepared);
            if(!componentStatus) return componentStatus;
            ++componentUse.acquired;
            masked[family].clear();
            selected.waveVortex.coefficients.Ap = family==0 ? WVComplexConstView{prepared,spectral_} : selected.waveVortex.coefficients.Ap;
            selected.waveVortex.coefficients.Am = family==1 ? WVComplexConstView{prepared,spectral_} : selected.waveVortex.coefficients.Am;
            selected.waveVortex.coefficients.A0 = family==2 ? WVComplexConstView{prepared,spectral_} : selected.waveVortex.coefficients.A0;
          } else {
            masked[family].resize(n);
            for(std::size_t index=0;index<n;++index)
              masked[family][index]=keep(group,family,index) ? coefficients[family][index] : WVComplex64{};
          }
        }
        if(!service.eventWorkspace_)
          selected.waveVortex={amplitudes.t,amplitudes.t0,{{masked[0].data(),spectral_},{masked[1].data(),spectral_},{masked[2].data(),spectral_}}};
        if(input.coefficientFamilyCount) {
          for(std::size_t family=0;family<3;++family) selectedFamilies[family]={
              input.coefficientFamilies[family].layout,
              service.eventWorkspace_ ? (family==0 ? selected.waveVortex.coefficients.Ap.data :
                  family==1 ? selected.waveVortex.coefficients.Am.data : selected.waveVortex.coefficients.A0.data) :
                  masked[family].data()};
          selected.coefficientFamilies=selectedFamilies.data();
        }
      }
      if(service.eventWorkspace_ && service.stateEvaluationActive_ && group) {
        const auto registrationStatus=service.addStateEvaluationView(
            selected,static_cast<std::size_t>(group));
        if(!registrationStatus) return registrationStatus;
        componentUse.registered=true;
      }
      if(service.eventWorkspace_) {
        service.eventWorkspace_->setComponent(static_cast<std::uint32_t>(group),&selected);
        componentUse.selected=true;
      }
      account(views.capacity()*sizeof(WVFieldOutputView));
      const auto status=service.evaluate(plan,selected,views.data(),views.size(),activeDependencies[group].data());
      if(!status) return status;
      const auto componentReleaseStatus=componentUse.release();
      if(!componentReleaseStatus) return componentReleaseStatus;
    }
    if(densityDemands) {
      auto* workspace=service.eventWorkspace_;
      const WVDensityEventGeometry geometry{&densityHeights_,&densityWeights_,&densityInitial_,Lz_,densityGravity_,densityReference_};
      bool preserveSource=false;
      for(std::size_t index=0;index<count;++index) {
        const auto& output=outputs_[index];
        preserveSource|=active(index) && !output.density && !output.forcing && !output.specification.isComplex &&
            output.group==0 && output.dependency==densityDependency_;
      }
      std::vector<double>* densitySource=nullptr;
      if(!workspace->hasDensitySource()) {
        if(densityDependency_>=fields[0].size())
          return invalid("Density source storage was not prepared for this event.");
        densitySource=&fields[0][densityDependency_];
      }
      const auto densityStatus=workspace->bindDensity(
          densitySource,spatial_,geometry,densityContract_,preserveSource);
      if(!densityStatus) return densityStatus;
      auto status=workspace->prepareDensity(densityDemands,densityContract_);
      if(!status) return status;
      if(needsAPV) {
        const auto R=spatial_.elementCount();
        status=workspace->prepareAPV(R,densityContract_.reference,[&](WVDensityEventView eta,double* apv) {
          const auto scratchStatus=checkout(scratchKey(
              scratchSignature_,0,0,19),3*R,
              densityDerivatives);
          if(!scratchStatus) return scratchStatus;
          account(0);
          WVRealVolumeConstView displacement{eta.data,spatial_};
          if(constant_) {
            WVRealFieldBundleView derivatives{densityDerivatives.data(),{spatial_.first,spatial_.second,spatial_.third,3}};
            const auto derivativeStatus=service.transform_->transformGGridScalarDerivatives(displacement,derivatives);
            if(!derivativeStatus) return derivativeStatus;
          } else {
            for(std::size_t axis=0;axis<3;++axis) {
              WVRealVolumeView derivative{densityDerivatives.data()+axis*R,spatial_};
              const auto derivativeStatus=hydrostatic_ ?
                (axis<2 ? service.stratified_->hydrostaticKernel_->differentiateHorizontal(displacement,axis==0,derivative) :
                  service.stratified_->hydrostaticKernel_->differentiateVertical(displacement,WVHydrostaticFamily::G,1,derivative)) :
                (axis<2 ? service.stratified_->boussinesqKernel_->differentiateHorizontal(displacement,axis==0,derivative) :
                  service.stratified_->boussinesqKernel_->differentiateVertical(displacement,WVBoussinesqFamily::G,1,derivative));
              if(!derivativeStatus) return derivativeStatus;
            }
          }
          const auto found=std::find_if(outputs_.begin(),outputs_.end(),[](const auto& output){return output.variable==Variable::apv;});
          const auto& zx=fields[0][found->auxiliaries[0]];
          const auto& zy=fields[0][found->auxiliaries[1]];
          const auto& zz=fields[0][found->auxiliaries[2]];
          const double f=constant_ ? constant_->descriptor().verticalModes().coriolisFrequency : 2*modal_->rotationRate*std::sin(modal_->latitude*std::acos(-1.0)/180.0);
          for(std::size_t i=0;i<R;++i) {
            apv[i]=zz[i]-zx[i]*densityDerivatives[i]-zy[i]*densityDerivatives[R+i]-(zz[i]+f)*densityDerivatives[2*R+i];
            if(!std::isfinite(apv[i])) return WVKernelStatus{WVKernelStatusCode::numericalFailure,"Available potential vorticity is not finite."};
          }
          return WVKernelStatus::ok();
        });
        if(!status) return status;
      }
    }
    if(std::any_of(activeForcing.begin(),activeForcing.end(),[](auto value){return value!=0;})) {
      const std::size_t channels=isQG_ ? 1 : isHydrostatic_ ? 3 : 4;
      forcingFields.resize(forcingIndices_.size());
      WVRealFieldBundleConstView prepared;
      if(physicalChannels) {
        const auto R=spatial_.elementCount();
        const auto scratchStatus=checkout(scratchKey(
            scratchSignature_,0,0,18),
            physicalChannels*R,forcingPhysical);
        if(!scratchStatus) return scratchStatus;
        for(std::size_t channel=0;channel<physicalChannels;++channel)
          std::copy_n(fields[0][forcingPhysicalDependencies_[channel]].data(),R,forcingPhysical.data()+channel*R);
        prepared={forcingPhysical.data(),{spatial_.first,spatial_.second,spatial_.third,physicalChannels}};
      }
      account(0);
      std::vector<WVVariableEvaluationKey> forcingKeys;
      std::vector<WVFieldOutputView> forcingCacheViews;
      forcingKeys.reserve(forcingIndices_.size());
      forcingCacheViews.reserve(forcingIndices_.size());
      forcingViews.reserve(forcingIndices_.size());
      for(std::size_t slot=0;slot<forcingIndices_.size();++slot) {
        if(!activeForcing[slot]) continue;
        auto& buffer=forcingFields[slot];
        const auto scratchStatus=checkout(scratchKey(
            scratchSignature_,0,slot,17),
            channels*spatial_.elementCount(),buffer);
        if(!scratchStatus) return scratchStatus;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::forcingTendency,
            static_cast<std::uint32_t>(isQG_ ?
                WVPortableVariable::Fqgpv_portable_catalog_forcing :
                WVPortableVariable::Fu_portable_catalog_forcing),0,0,0,
            static_cast<std::uint32_t>(forcingIndices_[slot]),1};
        if(service.eventWorkspace_ && service.eventWorkspace_->ready(key)) {
          bool reused=false;
          const auto copied=service.eventWorkspace_->evaluate(
              key,buffer.data(),buffer.size(),[]() {
                return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                    "A ready forcing tendency unexpectedly requested production."};
              },reused);
          if(!copied) return copied;
          continue;
        }
        forcingKeys.push_back(key);
        forcingCacheViews.push_back({buffer.data(),buffer.size()});
        forcingViews.push_back({forcingIndices_[slot],{buffer.data(),
            {spatial_.first,spatial_.second,spatial_.third,channels}}});
      }
      const auto operation=[&]() {
        WVForcingDiagnosticWorkspace* forcingWorkspace=nullptr;
        if(service.eventWorkspace_) {
          const auto workspaceStatus=service.eventWorkspace_->forcingWorkspace(
              *service.forcing_,forcingWorkspace);
          if(!workspaceStatus) return workspaceStatus;
        }
        return service.forcing_->evaluate(amplitudes,forcingViews.data(),forcingViews.size(),
            physicalChannels ? &prepared : nullptr,forcingWorkspace);
      };
      account(forcingKeys.capacity()*sizeof(WVVariableEvaluationKey)+
          forcingCacheViews.capacity()*sizeof(WVFieldOutputView));
      bool forcingReused=false;
      const auto status=forcingKeys.empty() ? WVKernelStatus::ok() :
          service.eventWorkspace_ ? service.eventWorkspace_->evaluateGroup(
              forcingKeys,forcingCacheViews,operation,forcingReused) : operation();
      metrics.diagnosticWorkspaceHighWaterBytes=std::max(metrics.diagnosticWorkspaceHighWaterBytes,
          metrics.diagnosticWorkspaceLiveBytes+service.forcing_->metrics().workspaceLastPeakBytes);
      if(!status) return status;
    }
    bool needsPhase=false;
    for(std::size_t index=0;index<count;++index) if(active(index))
      visitPortableExecution(outputs_[index].execution,[&](Variable node,bool) {
        needsPhase|=node==Variable::Apt || node==Variable::Amt ||
            node==Variable::phase || node==Variable::conjPhase;
      });
    const std::vector<WVComplex64>* phaseValues=nullptr;
    const WVComplex64* phaseData=nullptr;
    const WVVariableEvaluationKey phaseKey{WVVariableEvaluationNode::phaseFactors};
    struct ComplexUseGuard {
      WVFieldEvaluationEventWorkspace* workspace;
      WVVariableEvaluationKey key;
      bool armed=false;
      ~ComplexUseGuard() {if(armed) workspace->releaseComplex(key);}
    } phaseUse{service.eventWorkspace_,phaseKey};
    if(needsPhase) {
      const auto preparePhases=[&](WVComplex64* destination) {
        WVComplexConstView prepared;
        WVKernelStatus status;
        if(constant_) status=service.transform_->preparedPhase(amplitudes,prepared);
        else if(hydrostatic_)
          status=service.stratified_->hydrostaticKernel_->preparedPhase(
              amplitudes,prepared);
        else if(boussinesq_)
          status=service.stratified_->boussinesqKernel_->preparedPhase(
              amplitudes,prepared);
        else return WVKernelStatus{WVKernelStatusCode::unsupportedOperation,
            "This diagnostic transform has no phase factors."};
        if(!status) return status;
        std::copy_n(prepared.data,n,destination);
        return WVKernelStatus::ok();
      };
      if(service.eventWorkspace_) {
        const auto phaseStatus=service.eventWorkspace_->complexView(
            phaseKey,n,preparePhases,phaseValues);
        if(!phaseStatus) return phaseStatus;
        phaseUse.armed=true;
        phaseData=phaseValues->data();
      } else {
        phases.resize(n);
        const auto phaseStatus=preparePhases(phases.data());
        if(!phaseStatus) return phaseStatus;
        phaseValues=&phases;
        phaseData=phases.data();
      }
    }
    account(0);
    for(std::size_t index=0;index<count;++index) {
      if(!active(index)) continue;
      const auto& output=outputs_[index];
      if(output.sampled) {
        const double* source=nullptr;
        if(output.density) {
          const auto field=output.variable==Variable::rho_nm ? WVDensityEventField::rhoNm :
              output.variable==Variable::eta_true ? WVDensityEventField::etaTrue : WVDensityEventField::ape;
          const auto values=output.variable==Variable::apv ? service.eventWorkspace_->apvView(densityContract_.reference) :
              service.eventWorkspace_->densityView(field,densityContract_.reference);
          source=values.data;
        } else if(output.forcing) {
          source=forcingFields[output.forcingSlot].data()+
              output.forcingChannel*spatial_.elementCount();
        } else {
          const auto& field=fields[output.group][output.dependency];
          const auto offset=output.surface ?
              (spatial_.third-1)*spatial_.first*spatial_.second : 0;
          source=field.data()+offset;
        }
        const auto status=service.samplePreparedField(output.sampler,source,outputs[index]);
        if(!status) return status;
      } else if(output.density) {
        const auto field=output.variable==Variable::rho_nm ? WVDensityEventField::rhoNm :
            output.variable==Variable::eta_true ? WVDensityEventField::etaTrue : WVDensityEventField::ape;
        const auto values=output.variable==Variable::apv ? service.eventWorkspace_->apvView(densityContract_.reference) :
            service.eventWorkspace_->densityView(field,densityContract_.reference);
        std::copy_n(values.data,values.elementCount,outputs[index].data);
      } else if(output.forcing) {
        const auto R=spatial_.elementCount();
        std::copy_n(forcingFields[output.forcingSlot].data()+output.forcingChannel*R,R,outputs[index].data);
      } else if(output.specification.isComplex) {
        const auto produce=[&](WVComplex64* destination) {
          if(output.variable==Variable::phase || output.variable==Variable::conjPhase) {
            for(std::size_t coefficient=0;coefficient<n;++coefficient) {
              auto phase=phaseData[coefficient];
              if(output.variable==Variable::conjPhase) phase.imag=-phase.imag;
              destination[coefficient]=phase;
            }
          } else {
            const auto family=output.variable==Variable::Apt ? 0 : output.variable==Variable::Amt ? 1 : 2;
            for(std::size_t coefficient=0;coefficient<n;++coefficient) {
              const auto a=coefficients[family][coefficient];
              if(family==2) destination[coefficient]=a;
              else {
                auto phase=phaseData[coefficient]; if(family==1) phase.imag=-phase.imag;
                destination[coefficient]={a.real*phase.real-a.imag*phase.imag,a.real*phase.imag+a.imag*phase.real};
              }
            }
            }
          return WVKernelStatus::ok();
        };
        if(service.eventWorkspace_) {
          const auto finalNode=output.execution.count ?
              output.execution.order[output.execution.count-1] : output.variable;
          const WVVariableEvaluationKey key{WVVariableEvaluationNode::registeredVariable,
              static_cast<std::uint32_t>(finalNode),static_cast<std::uint32_t>(output.group)};
          const std::vector<WVComplex64>* values=nullptr;
          const auto status=service.eventWorkspace_->complexView(key,n,produce,values);
          if(!status) return status;
          ComplexUseGuard outputUse{service.eventWorkspace_,key,true};
          std::copy(values->begin(),values->end(),outputs[index].complexData);
        } else {
          const auto status=produce(outputs[index].complexData);
          if(!status) return status;
        }
      } else if(output.variable==Variable::totalEnergySpatiallyIntegrated) {
        const auto produce=[&]() {
          const auto& u=fields[0][output.auxiliaries[0]];
          const auto& v=fields[0][output.auxiliaries[1]];
          const auto& eta=fields[0][output.auxiliaries[2]];
          const auto horizontal=spatial_.first*spatial_.second;
          double integral=0;
          for(std::size_t z=0;z<spatial_.third;++z) {
            double horizontalSum=0;
            for(std::size_t point=0;point<horizontal;++point) {
              const auto i=point+horizontal*z;
              double value=u[i]*u[i]+v[i]*v[i]+stratification(z)*eta[i]*eta[i];
              if(!isHydrostatic_) { const double w=fields[0][output.auxiliaries[3]][i]; value+=w*w; }
              horizontalSum+=value;
            }
            integral+=weight(z)*horizontalSum/static_cast<double>(horizontal);
          }
          outputs[index].data[0]=integral/2;
          ++producerMetrics.energyReductions;
          return WVKernelStatus::ok();
        };
        bool reused=false;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(output.variable)};
        const auto status=service.eventWorkspace_ ? service.eventWorkspace_->evaluate(
            key,outputs[index].data,1,produce,reused) : produce();
        if(!status) return status;
      } else if(output.verticalMean) {
        const auto produce=[&]() {
          const auto& source=fields[0][output.dependency];
          const auto plane=spatial_.first*spatial_.second;
          for(std::size_t z=0;z<spatial_.third;++z) {
            double sum=0;
            for(std::size_t i=0;i<plane;++i) sum+=source[z*plane+i];
            outputs[index].data[z]=sum/static_cast<double>(plane);
          }
          return WVKernelStatus::ok();
        };
        bool reused=false;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(output.variable)};
        const auto status=service.eventWorkspace_ ? service.eventWorkspace_->evaluate(
            key,outputs[index].data,outputs[index].elementCount,produce,reused) : produce();
        if(!status) return status;
      } else if(output.extrema) {
        const auto produce=[&]() {
          const auto& first=fields[0][output.auxiliaries[0]];
          double maximum=0;
          for(std::size_t point=0;point<first.size();++point) {
            const double value=output.variable==Variable::uvMax ?
                (constant_ || isBarotropic_ ? std::sqrt(first[point]*first[point]+fields[0][output.auxiliaries[1]][point]*fields[0][output.auxiliaries[1]][point]) : std::hypot(first[point],fields[0][output.auxiliaries[1]][point])) : std::abs(first[point]);
            maximum=std::max(maximum,value);
          }
          outputs[index].data[0]=maximum;
          if(output.variable==Variable::uvMax) {
            ++producerMetrics.horizontalSpeedReductions;
          } else ++producerMetrics.verticalSpeedReductions;
          return WVKernelStatus::ok();
        };
        bool reused=false;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(output.variable)};
        const auto status=service.eventWorkspace_ ? service.eventWorkspace_->evaluate(
            key,outputs[index].data,1,produce,reused) : produce();
        if(!status) return status;
      } else {
        const auto& source=fields[output.group][output.dependency];
        const auto offset=output.surface ? (spatial_.third-1)*spatial_.first*spatial_.second : 0;
        std::copy_n(source.data()+offset,output.specification.elementCount,outputs[index].data);
      }
    }
    ++metrics.diagnosticEvaluationCount;
    metrics.diagnosticPrimitiveOutputCount+=primitiveCount;
    std::size_t primitiveReferences=physicalChannels,phaseReferences=0;
    for(std::size_t index=0;index<count;++index) {
      if(!active(index)) continue;
      const auto& output=outputs_[index];
      if(output.forcing) continue;
      if(output.specification.isComplex) {
        if(output.variable!=Variable::A0t) ++phaseReferences;
      } else if(output.variable==Variable::totalEnergySpatiallyIntegrated) primitiveReferences+=isHydrostatic_ ? 3 : 4;
      else if(output.extrema) primitiveReferences+=output.variable==Variable::uvMax ? 2 : 1;
      else ++primitiveReferences;
    }
    metrics.diagnosticIntermediateReuseCount+=primitiveReferences-primitiveCount+(phaseReferences ? phaseReferences-1 : 0);
    metrics.diagnosticWorkspaceLiveBytes=outerWorkspaceBytes;
    // Nested primitive plans write into event scratch.  Public write accounting
    // describes caller-owned outputs, independent of how many cached
    // dependencies were materialized to produce them.
    metrics.outputElementWriteCount=outputWritesBefore+activeOutputElements;
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,"Unable to allocate event-scoped diagnostic scratch."};
  }
}

std::size_t WVDiagnosticFieldPlan::persistentBytes() const noexcept {
  std::size_t bytes=sizeof(*this)+configuration_.capacity()+outputs_.capacity()*sizeof(Output)+forcingIndices_.capacity()*sizeof(std::size_t);
  bytes+=(densityHeights_.capacity()+densityWeights_.capacity()+densityInitial_.capacity())*sizeof(double);
  for(const auto& output:outputs_)
    bytes+=output.specification.identifier.capacity()+output.specification.fieldName.capacity()+
        output.specification.dimensions.capacity()*sizeof(std::size_t)+
        (output.sampled ? output.sampler.persistentBytes()-sizeof(output.sampler) : 0)+
        (output.sampling.x.capacity()+output.sampling.y.capacity()+output.sampling.z.capacity())*sizeof(double)+
        (output.sampling.xIndices.capacity()+output.sampling.yIndices.capacity())*sizeof(std::size_t);
  for(const auto& group:groups_) {
    bytes+=group.fields.persistentBytes()-sizeof(group.fields)+group.requests.capacity()*sizeof(WVFieldRequest);
    for(const auto& request:group.requests)
      bytes+=request.identifier.capacity()+request.fieldName.capacity()+
          (request.sampling.x.capacity()+request.sampling.y.capacity()+request.sampling.z.capacity())*sizeof(double)+
          (request.sampling.xIndices.capacity()+request.sampling.yIndices.capacity())*sizeof(std::size_t);
  }
  return bytes;
}
} // namespace wavevortex::runtime::detail
