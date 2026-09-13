#pragma once
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WaveVortexRuntime/WVForcingTendency.hpp"
#include "WaveVortexRuntime/WVForcing.hpp"
#include "WaveVortexRuntime/WVPortableVariablePlan.hpp"
#include "WVForcingDiagnosticWorkspace.hpp"
#include <array>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace wavevortex::runtime {
class WVConstantStratificationForcingEngine;
class WVBarotropicQGForcingEngine;
class WVStratifiedQGForcingEngine;
class WVHydrostaticForcingEngine;
class WVBoussinesqForcingEngine;
}

namespace wavevortex::runtime::detail {

// Borrows an existing forcing engine and its immutable, already resolved
// instances. This is a coarse service binding, never a second forcing registry.
class WVForcingDiagnosticBinding final {
public:
  struct Identity {
    std::string_view type,name;
    std::uint32_t version;
    WVForcingStage stage;
    std::uint8_t priority;
    std::size_t ordinal;
  };
  struct Output {
    const WVPortableVariableContract* contract=nullptr;
    WVPortableVariablePlan plan;
    std::size_t executionIndex=0,channel=0,physicalChannels=0;
    std::string_view instanceName;
    std::uint32_t instanceOrdinal=0;
  };
  template<class Engine,bool QG>
  static WVKernelStatus create(Engine& engine,std::unique_ptr<WVForcingDiagnosticBinding>& result) {
    try {
      auto candidate=std::unique_ptr<WVForcingDiagnosticBinding>(new WVForcingDiagnosticBinding);
      candidate->engine_=&engine;
      bool prefixQualified=true,physicalNeeded=false;
      for(std::size_t index=0;index<engine.forcingCount();++index) {
        const auto* instance=engine.forcingInstance(index);
        const auto* dependencies=engine.forcingEvaluationDependencies(index);
        if(!dependencies)
          return {WVKernelStatusCode::invalidConfiguration,
              "Forcing evaluation dependencies are unavailable."};
        candidate->nonlinearUseCount_+=dependencies->nonlinearUseCount;
        candidate->horizontalMaximumNeeded_|=dependencies->horizontalMaximumNeeded;
        if constexpr(std::is_same_v<Engine,
            WVConstantStratificationForcingEngine>) {
          for(std::size_t i=0;i<candidate->constantLaplacianUseCount_.size();++i)
            candidate->constantLaplacianUseCount_[i]+=
                dependencies->constantLaplacianUseCount[i];
        } else if constexpr(std::is_same_v<Engine,
            WVHydrostaticForcingEngine>) {
          for(std::size_t i=0;i<candidate->gridCalculusUseCount_.size();++i)
            candidate->gridCalculusUseCount_[i]+=
                dependencies->hydrostaticGridCalculusUseCount[i];
        } else if constexpr(std::is_same_v<Engine,
            WVBoussinesqForcingEngine>) {
          for(std::size_t i=0;i<candidate->gridCalculusUseCount_.size();++i)
            candidate->gridCalculusUseCount_[i]+=
                dependencies->boussinesqGridCalculusUseCount[i];
        }
        if(instance->ordinal()>std::numeric_limits<std::uint32_t>::max())
          return {WVKernelStatusCode::invalidConfiguration,"Forcing instance ordinal exceeds the graph contract."};
        for(const auto& previous:candidate->bindings_) if(previous.instanceOrdinal==instance->ordinal())
          return {WVKernelStatusCode::invalidConfiguration,"Forcing graph ordinals must be unique for diagnostic binding."};
        prefixQualified&=instance->supportsTendencyDiagnostics();
        physicalNeeded|=instance->requiresDiagnosticPhysicalFields();
        candidate->physicalPrefix_.push_back(physicalNeeded);
        candidate->bindings_.push_back({instance->name(),static_cast<std::uint32_t>(instance->ordinal()),prefixQualified});
        candidate->stages_.push_back(instance->stage());
      }
      candidate->physicalChannels_=QG ? 2 : std::is_same_v<Engine,WVConstantStratificationForcingEngine> ? 3 : 4;
      candidate->evaluate_=[](void* pointer,const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,const WVRealFieldBundleConstView* prepared,WVForcingDiagnosticWorkspace* session) {
        auto& resolved=*static_cast<Engine*>(pointer);
        if constexpr(QG) return resolved.evaluateForcingTendencies(state.coefficients.A0,outputs,count,prepared,session);
        else return resolved.evaluateForcingTendencies(state,outputs,count,prepared,session);
      };
      if constexpr(!QG) candidate->evaluateBuiltin_=[](void* pointer,const WVState& state,
          const WVRealFieldBundleConstView* prepared,WVForcingDiagnosticWorkspace* session,
          WVFlux& flux) {
        return static_cast<Engine*>(pointer)->evaluateForcingTendenciesImpl(
            state,nullptr,0,prepared,session,&flux);
      };
      candidate->metrics_=[](const void* pointer)->const WVForcingTendencyMetrics& {return static_cast<const Engine*>(pointer)->tendencyMetrics();};
      candidate->instance_=[](const void* pointer,std::size_t index)->Identity {
        const auto* instance=static_cast<const Engine*>(pointer)->forcingInstance(index);
        return {instance->typeIdentifier(),instance->name(),instance->contractVersion(),
            instance->stage(),instance->priority(),instance->ordinal()};
      };
      // The engine owns an immutable resolved schedule, so this setup-time
      // identity remains valid for the lifetime of the borrowed binding.
      if(candidate->evaluateBuiltin_ && candidate->bindings_.size()==1 &&
          candidate->stages_.size()==1) {
        const auto defaultNonlinear=defaultNonlinearAdvectionSchedule();
        if(defaultNonlinear.entries.size()==1) {
          const auto actual=candidate->instance_(candidate->engine_,0);
          const auto& expected=defaultNonlinear.entries.front();
          candidate->exactBuiltinNonlinear_=
              actual.type==expected.typeIdentifier &&
              actual.version==expected.contractVersion &&
              actual.name==expected.name && actual.stage==expected.stage &&
              candidate->stages_[0]==expected.stage &&
              actual.priority==expected.priority &&
              actual.ordinal==expected.ordinal;
        }
      }
      if constexpr(std::is_same_v<Engine,
          WVConstantStratificationForcingEngine>) {
        const auto& descriptor=engine.kernel().descriptor();
        const auto shape=descriptor.spatialShape();
        candidate->spectral_=descriptor.spectralShape();
        const auto& configuration=descriptor.configuration();
        candidate->requiresFourChannelTendencySelection_=configuration.isHydrostatic;
        candidate->spatial_={shape.first,shape.second,shape.third,
            configuration.isHydrostatic ? 3u : 4u};
      } else if constexpr(std::is_same_v<Engine,
          WVBarotropicQGForcingEngine>) {
        const auto& descriptor=engine.kernel().descriptor();
        const auto shape=descriptor.spatialShape();
        candidate->spectral_=descriptor.spectralShape();
        candidate->spatial_={shape.rows,shape.columns,1,1};
      } else if constexpr(std::is_same_v<Engine,
          WVStratifiedQGForcingEngine>) {
        const auto shape=engine.kernel().spatialShape();
        candidate->spectral_=engine.kernel().spectralShape();
        candidate->spatial_={shape.first,shape.second,shape.third,1};
      } else if constexpr(std::is_same_v<Engine,
          WVHydrostaticForcingEngine>) {
        const auto shape=engine.kernel().spatialShape();
        candidate->spectral_=engine.kernel().spectralShape();
        candidate->spatial_={shape.first,shape.second,shape.third,3};
      } else {
        const auto shape=engine.kernel().spatialShape();
        candidate->spectral_=engine.kernel().spectralShape();
        candidate->spatial_={shape.first,shape.second,shape.third,4};
      }
      candidate->coefficientFamilies_=QG ? 1 : 3;
      result=std::move(candidate);
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,"Unable to bind resolved forcing diagnostics."};
    }
  }
  WVKernelStatus resolve(std::string_view name,std::string_view configuration,std::uint8_t sampling,Output& output) const {
    WVPortableVariableOptions options; options.source=WVPortableOperationSource::builtIn;
    // The source-linked engine implementation qualifies numerical execution.
    // Catalog delivery status records the qualified output implementation.
    WVPortableVariablePlan plan;
    const auto status=resolvePortableForcingVariablePlan(name,configuration,sampling,bindings_.data(),bindings_.size(),options,plan);
    if(status!=WVPortableVariableStatus::supported)
      return {WVKernelStatusCode::unsupportedOperation,"Forcing diagnostic identity, sampling, or instance binding is unsupported or ambiguous: "+std::string(name)};
    for(std::size_t index=0;index<bindings_.size();++index) if(bindings_[index].instanceOrdinal==plan.forcingInstanceOrdinal) {
      output.contract=plan.output; output.plan=plan; output.executionIndex=index;
      output.instanceName=bindings_[index].instanceName; output.instanceOrdinal=bindings_[index].instanceOrdinal;
      output.physicalChannels=physicalPrefix_[index] ? physicalChannels_ : 0;
      switch(plan.output->metadata.identifier) {
        case WVPortableVariable::Fu_portable_catalog_forcing: case WVPortableVariable::Fqgpv_portable_catalog_forcing: output.channel=0; break;
        case WVPortableVariable::Fv_portable_catalog_forcing: output.channel=1; break;
        case WVPortableVariable::Fw_portable_catalog_forcing: output.channel=2; break;
        case WVPortableVariable::Feta_portable_catalog_forcing:
          output.channel=configuration.rfind("constant-nonhydrostatic-",0)==0 || configuration.rfind("boussinesq-",0)==0 ? 3 : 2; break;
        default: return {WVKernelStatusCode::unsupportedOperation,"Unknown forcing diagnostic channel."};
      }
      return WVKernelStatus::ok();
    }
    return {WVKernelStatusCode::invalidConfiguration,"Resolved forcing instance is absent from its execution schedule."};
  }
  WVKernelStatus evaluate(const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,const WVRealFieldBundleConstView* prepared,WVForcingDiagnosticWorkspace* session) const {
    return evaluate_(engine_,state,outputs,count,prepared,session);
  }
  bool supportsExactBuiltinNonlinear() const noexcept {
    return exactBuiltinNonlinear_;
  }
  WVKernelStatus evaluateBuiltinNonlinear(const WVState& state,
      const WVRealFieldBundleConstView* prepared,
      WVForcingDiagnosticWorkspace* session,WVFlux& flux) const {
    if(!supportsExactBuiltinNonlinear())
      return {WVKernelStatusCode::unsupportedOperation,
          "The resolved forcing schedule is not the exact built-in nonlinear advection."};
    return evaluateBuiltin_(engine_,state,prepared,session,flux);
  }
  std::unique_ptr<WVForcingDiagnosticWorkspace> createWorkspace() const {
    auto workspace=std::make_unique<WVForcingDiagnosticWorkspace>(
        spectral_,spatial_,coefficientFamilies_,coefficientFamilies_==1 ? 2 : 4);
    workspace->nonlinearUseCount=nonlinearUseCount_;
    workspace->requiresFourChannelTendencySelection=requiresFourChannelTendencySelection_;
    workspace->gridCalculusUseCount=gridCalculusUseCount_;
    workspace->constantLaplacianUseCount=constantLaplacianUseCount_;
    return workspace;
  }
  const std::vector<WVForcingStage>& stages() const noexcept {return stages_;}
  bool horizontalMaximumNeeded() const noexcept {
    return horizontalMaximumNeeded_;
  }
  std::vector<WVVariableEvaluationKey> dependencyKeys() const {
    return WVForcingDiagnosticWorkspace::dependencyKeys(bindings_.size(),
        &gridCalculusUseCount_,&constantLaplacianUseCount_);
  }
  bool hasSamePrefix(const WVForcingDiagnosticBinding& other,std::size_t index,std::size_t otherIndex) const noexcept {
    if(index!=otherIndex || index>=bindings_.size() || otherIndex>=other.bindings_.size()) return false;
    for(std::size_t position=0;position<=index;++position) {
      const auto first=instance_(engine_,position);
      const auto second=other.instance_(other.engine_,position);
      if(first.type!=second.type || first.version!=second.version || first.stage!=second.stage ||
          first.priority!=second.priority || first.ordinal!=second.ordinal || first.name!=second.name) return false;
    }
    return true;
  }
  const WVForcingTendencyMetrics& metrics() const {return metrics_(engine_);}
  const std::vector<WVPortableForcingVariableBinding>& bindings() const noexcept {return bindings_;}
  std::size_t fullSchedulePhysicalChannels() const noexcept {
    return !physicalPrefix_.empty() && physicalPrefix_.back() ? physicalChannels_ : 0;
  }
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this)+bindings_.capacity()*sizeof(WVPortableForcingVariableBinding)+
        physicalPrefix_.capacity()*sizeof(std::uint8_t)+
        stages_.capacity()*sizeof(WVForcingStage);
  }
private:
  void* engine_=nullptr;
  WVKernelStatus (*evaluate_)(void*,const WVState&,const WVForcingTendencyOutput*,std::size_t,const WVRealFieldBundleConstView*,WVForcingDiagnosticWorkspace*)=nullptr;
  WVKernelStatus (*evaluateBuiltin_)(void*,const WVState&,
      const WVRealFieldBundleConstView*,WVForcingDiagnosticWorkspace*,WVFlux&)=nullptr;
  const WVForcingTendencyMetrics& (*metrics_)(const void*)=nullptr;
  Identity (*instance_)(const void*,std::size_t)=nullptr;
  std::vector<WVPortableForcingVariableBinding> bindings_;
  std::vector<std::uint8_t> physicalPrefix_;
  std::vector<WVForcingStage> stages_;
  WVShape2D spectral_{};
  WVShape4D spatial_{};
  std::size_t coefficientFamilies_=3;
  std::size_t physicalChannels_=0;
  std::size_t nonlinearUseCount_=0;
  std::array<std::size_t,16> gridCalculusUseCount_{};
  std::array<std::size_t,2> constantLaplacianUseCount_{};
  bool requiresFourChannelTendencySelection_=false;
  bool horizontalMaximumNeeded_=false;
  bool exactBuiltinNonlinear_=false;
};
}
