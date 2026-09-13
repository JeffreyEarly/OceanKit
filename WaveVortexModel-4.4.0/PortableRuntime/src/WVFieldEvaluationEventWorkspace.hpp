#pragma once

#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"
#include "WVDensityEventEvaluation.hpp"
#include "WVForcingDiagnosticBinding.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wavevortex::runtime::detail {

struct WVFieldEvaluationRealEntry {
  WVVariableEvaluationKey key;
  WVVariableEvaluationKey preparedKey;
  std::vector<double> values;
  WVShape3D volumeShape{};
  std::size_t preparedElements=0;
  bool assigned=false,prepared=false;
};
struct WVFieldEvaluationComplexEntry {
  WVVariableEvaluationKey key;
  WVVariableEvaluationKey preparedKey;
  std::vector<WVComplex64> values;
  std::size_t preparedElements=0;
  bool assigned=false,prepared=false;
};
class WVFieldEvaluationArena final {
public:
  WVKernelStatus prepareGroupNodesForActive(std::size_t count) {
    if(!groupNodes.empty())
      return {WVKernelStatusCode::reentrantExecution,
          "Cannot extend field preparation during a fused producer."};
    return prepareGroupNodes(count);
  }
  WVKernelStatus prepareRealForActive(const WVVariableEvaluationKey& key,
      std::size_t elements) {
    try {
      for(auto& entry:realFields) if(entry->prepared && entry->preparedKey==key) {
        if(entry->assigned && entry->values.size()!=elements)
          return {WVKernelStatusCode::invalidConfiguration,
              "Active field preparation would resize an assigned real field."};
        if(!entry->assigned) entry->values.reserve(elements);
        entry->preparedElements=std::max(entry->preparedElements,elements);
        notePrepared(); return WVKernelStatus::ok();
      }
      for(auto& entry:realFields) if(entry->assigned && entry->key==key) {
        if(entry->values.size()!=elements)
          return {WVKernelStatusCode::invalidConfiguration,
              "Active field preparation would resize an assigned real field."};
        entry->prepared=true; entry->preparedKey=key;
        entry->preparedElements=elements;
        notePrepared(); return WVKernelStatus::ok();
      }
      auto entry=std::make_unique<WVFieldEvaluationRealEntry>();
      entry->prepared=true; entry->preparedKey=key;
      entry->preparedElements=elements; entry->values.reserve(elements);
      realFields.push_back(std::move(entry)); notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to extend active real field preparation."};
    } catch(const std::length_error&) {
      return {WVKernelStatusCode::sizeOverflow,
          "Active real field preparation exceeds vector capacity."};
    }
  }
  WVKernelStatus prepareComplexForActive(const WVVariableEvaluationKey& key,
      std::size_t elements) {
    try {
      for(auto& entry:complexFields) if(entry->prepared && entry->preparedKey==key) {
        if(entry->assigned && entry->values.size()!=elements)
          return {WVKernelStatusCode::invalidConfiguration,
              "Active field preparation would resize an assigned complex field."};
        if(!entry->assigned) entry->values.reserve(elements);
        entry->preparedElements=std::max(entry->preparedElements,elements);
        notePrepared(); return WVKernelStatus::ok();
      }
      for(auto& entry:complexFields) if(entry->assigned && entry->key==key) {
        if(entry->values.size()!=elements)
          return {WVKernelStatusCode::invalidConfiguration,
              "Active field preparation would resize an assigned complex field."};
        entry->prepared=true; entry->preparedKey=key;
        entry->preparedElements=elements;
        notePrepared(); return WVKernelStatus::ok();
      }
      auto entry=std::make_unique<WVFieldEvaluationComplexEntry>();
      entry->prepared=true; entry->preparedKey=key;
      entry->preparedElements=elements; entry->values.reserve(elements);
      complexFields.push_back(std::move(entry)); notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to extend active complex field preparation."};
    } catch(const std::length_error&) {
      return {WVKernelStatusCode::sizeOverflow,
          "Active complex field preparation exceeds vector capacity."};
    }
  }
  WVKernelStatus prepareForcingForActive(
      const WVForcingDiagnosticBinding& binding,
      WVVariableEvaluationPolicy policy) {
    try {
      if(!forcingWorkspace) forcingWorkspace=binding.createWorkspace();
      const auto status=forcingWorkspace->prepareScopedStorageForActive(
          policy,binding.stages());
      if(status) notePrepared();
      return status;
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to extend active forcing diagnostic preparation."};
    }
  }
  WVKernelStatus prepareDensityForActive(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference,bool apvNeeded) {
    try {
      const bool bound=density[0].initialized() || density[1].initialized();
      if(bound) {
        if(densitySource.size()!=sampleCount ||
            densityHeights.size()!=profileCount ||
            densityWeights.size()!=profileCount ||
            densityInitial.size()!=profileCount)
          return {WVKernelStatusCode::invalidConfiguration,
              "Active density preparation has incompatible event geometry."};
      } else {
        densitySource.reserve(sampleCount);
        densityHeights.reserve(profileCount);
        densityWeights.reserve(profileCount);
        densityInitial.reserve(profileCount);
      }
      const auto selected=reference==WVNoMotionReference::initial ? 1u : 0u;
      const auto derivedDemands=static_cast<std::uint8_t>(demands&
          (WVDensityEventEvaluation::etaTrueDemand|
           WVDensityEventEvaluation::apeDemand));
      if(demands&WVDensityEventEvaluation::rhoNmDemand) {
        const auto status=density[0].reserveAdditionalStorage(sampleCount,
            profileCount,WVDensityEventEvaluation::rhoNmDemand,
            WVNoMotionReference::actual);
        if(!status) return status;
      }
      if(derivedDemands) {
        const auto status=density[selected].reserveAdditionalStorage(
            sampleCount,profileCount,derivedDemands,reference);
        if(!status) return status;
      }
      if(apvNeeded) apv[selected].reserve(sampleCount);
      notePrepared(); return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to extend active density preparation."};
    } catch(const std::length_error&) {
      return {WVKernelStatusCode::sizeOverflow,
          "Active density preparation exceeds vector capacity."};
    }
  }
  WVKernelStatus prepareGroupNodes(std::size_t count) {
    try {
      groupNodes.reserve(count);
      notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare fused output-event metadata."};
    }
  }
  WVKernelStatus prepareReal(const WVVariableEvaluationKey& key,
      std::size_t elements) {
    try {
      for(auto& entry:realFields) if(entry->prepared && entry->preparedKey==key) {
        entry->preparedElements=std::max(entry->preparedElements,elements);
        if(policy==WVVariableEvaluationPolicy::reuse ||
            (key.geometry&(1ULL<<63)))
          entry->values.reserve(entry->preparedElements);
        notePrepared(); return WVKernelStatus::ok();
      }
      auto entry=std::make_unique<WVFieldEvaluationRealEntry>();
      entry->prepared=true; entry->preparedKey=key;
      entry->preparedElements=elements;
      if(policy==WVVariableEvaluationPolicy::reuse ||
          (key.geometry&(1ULL<<63)))
        entry->values.reserve(elements);
      realFields.push_back(std::move(entry)); notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare a real output-event field."};
    }
  }
  WVKernelStatus prepareComplex(const WVVariableEvaluationKey& key,
      std::size_t elements) {
    try {
      for(auto& entry:complexFields) if(entry->prepared && entry->preparedKey==key) {
        entry->preparedElements=std::max(entry->preparedElements,elements);
        entry->values.reserve(entry->preparedElements);
        notePrepared(); return WVKernelStatus::ok();
      }
      auto entry=std::make_unique<WVFieldEvaluationComplexEntry>();
      entry->prepared=true; entry->preparedKey=key;
      entry->preparedElements=elements; entry->values.reserve(elements);
      complexFields.push_back(std::move(entry)); notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare a complex output-event field."};
    }
  }
  WVKernelStatus prepareForcing(const WVForcingDiagnosticBinding& binding,
      WVVariableEvaluationPolicy policy) {
    try {
      if(!forcingWorkspace) forcingWorkspace=binding.createWorkspace();
      const auto status=forcingWorkspace->prepareScopedStorage(
          policy,binding.stages());
      if(!status) return status;
      plannedBytes=persistentBytes();
      peakBytes=std::max(peakBytes,plannedBytes);
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare the output-event forcing arena."};
    }
  }
  WVKernelStatus prepareDensity(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference,bool apvNeeded) {
    try {
      densitySource.reserve(sampleCount);
      densityHeights.reserve(profileCount);
      densityWeights.reserve(profileCount);
      densityInitial.reserve(profileCount);
      const auto selected=reference==WVNoMotionReference::initial ? 1u : 0u;
      const auto derivedDemands=static_cast<std::uint8_t>(demands&
          (WVDensityEventEvaluation::etaTrueDemand|
           WVDensityEventEvaluation::apeDemand));
      if(policy==WVVariableEvaluationPolicy::reuse) {
        if(demands&WVDensityEventEvaluation::rhoNmDemand) {
          const auto status=density[0].reserveStorage(sampleCount,profileCount,
              WVDensityEventEvaluation::rhoNmDemand,
              WVNoMotionReference::actual);
          if(!status) return status;
        }
        if(derivedDemands) {
          const auto status=density[selected].reserveStorage(
              sampleCount,profileCount,derivedDemands,reference);
          if(!status) return status;
        }
        if(apvNeeded) apv[selected].reserve(sampleCount);
      } else {
        const auto actualDemands=static_cast<std::uint8_t>(
            (demands&WVDensityEventEvaluation::rhoNmDemand)|
            (reference==WVNoMotionReference::actual ? derivedDemands : 0));
        const auto initialDemands=static_cast<std::uint8_t>(
            reference==WVNoMotionReference::initial ? derivedDemands : 0);
        auto status=density[0].reserveLowMemoryStorage(sampleCount,profileCount,
            actualDemands,WVNoMotionReference::actual);
        if(!status) return status;
        status=density[1].reserveLowMemoryStorage(sampleCount,profileCount,
            initialDemands,WVNoMotionReference::initial);
        if(!status) return status;
        for(std::size_t index=0;index<apv.size();++index) {
          if(apvNeeded && index==selected) apv[index].reserve(sampleCount);
          else std::vector<double>{}.swap(apv[index]);
        }
      }
      notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare the density output-event arena."};
    }
  }
  bool hasPreparedReal(const WVVariableEvaluationKey& key) const noexcept {
    return std::any_of(realFields.begin(),realFields.end(),[&](const auto& entry) {
      return entry->prepared && entry->preparedKey==key;
    });
  }
  void setPreparedVolumeShape(const WVVariableEvaluationKey& key,
      WVShape3D shape) noexcept {
    for(auto& entry:realFields) if(entry->prepared && entry->preparedKey==key) {
      entry->volumeShape=shape;
      return;
    }
  }
  void clearReal() noexcept {
    policy=WVVariableEvaluationPolicy::lowMemory;
    for(auto& entry:realFields)
      if(!(entry->preparedKey.geometry&(1ULL<<63)))
        std::vector<double>{}.swap(entry->values);
    for(auto& value:density) value.release();
    for(auto& value:apv) std::vector<double>{}.swap(value);
    plannedBytes=persistentBytes();
    peakBytes=plannedBytes;
  }
  WVKernelStatus preparePolicy(WVVariableEvaluationPolicy selected) {
    policy=selected;
    if(selected==WVVariableEvaluationPolicy::lowMemory) {
      clearReal();
      return WVKernelStatus::ok();
    }
    try {
      for(auto& entry:realFields) if(entry->prepared)
        entry->values.reserve(entry->preparedElements);
      notePrepared();
      return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to restore the prepared output-event arena."};
    }
  }
  std::size_t persistentBytes() const noexcept {
    std::size_t bytes=sizeof(*this)+
        realFields.capacity()*sizeof(std::unique_ptr<WVFieldEvaluationRealEntry>)+
        complexFields.capacity()*sizeof(std::unique_ptr<WVFieldEvaluationComplexEntry>)+
        groupNodes.capacity()*sizeof(decltype(groupNodes)::value_type);
    for(const auto& entry:realFields)
      bytes+=sizeof(*entry)+entry->values.capacity()*sizeof(double);
    for(const auto& entry:complexFields)
      bytes+=sizeof(*entry)+entry->values.capacity()*sizeof(WVComplex64);
    bytes+=(densitySource.capacity()+densityHeights.capacity()+
        densityWeights.capacity()+densityInitial.capacity()+
        apv[0].capacity()+apv[1].capacity())*sizeof(double);
    bytes+=density[0].metrics().liveBytes+density[1].metrics().liveBytes;
    if(forcingWorkspace) bytes+=forcingWorkspace->bytes();
    return bytes;
  }
  void notePeak() noexcept {peakBytes=std::max(peakBytes,persistentBytes());}
  void refreshAccounting() noexcept {notePrepared();}
  std::vector<std::unique_ptr<WVFieldEvaluationRealEntry>> realFields;
  std::vector<std::unique_ptr<WVFieldEvaluationComplexEntry>> complexFields;
  std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> groupNodes;
  std::unique_ptr<WVForcingDiagnosticWorkspace> forcingWorkspace;
  std::vector<double> densitySource;
  std::vector<double> densityHeights,densityWeights,densityInitial;
  std::array<WVDensityEventEvaluation,2> density;
  std::array<std::vector<double>,2> apv;
  std::size_t plannedBytes=0,peakBytes=0;
  WVVariableEvaluationPolicy policy=WVVariableEvaluationPolicy::reuse;
private:
  void notePrepared() noexcept {
    plannedBytes=persistentBytes();
    peakBytes=std::max(peakBytes,plannedBytes);
  }
};

// Expensive canonical DAG nodes shared only within one output preparation.
// The owning service fixes the transform; the borrowed state fixes time and
// coefficients. Component identity is explicit, so masked reconstructions can
// share within the event without making pointer identity part of the cache key.
// Nothing survives the scope, so an in-place mutation needs a new session.
class WVFieldEvaluationEventWorkspace final {
public:
  static constexpr std::size_t dynamicalDerivativeKeyBase=32;
  static WVKernelStatus prepareEvaluationContext(WVVariableEvaluationContext& context) {
    std::vector<WVVariableEvaluationKey> keys;
    try {
      keys.reserve(WVPortableVariableCatalog.size()*18+32);
      for(const auto& variable:WVPortableVariableCatalog) {
        for(std::uint32_t component=0;component<5;++component) {
          keys.push_back({WVVariableEvaluationNode::registeredVariable,variable.ordinal,component});
          keys.push_back({WVVariableEvaluationNode::physicalField,variable.ordinal,component});
          keys.push_back({WVVariableEvaluationNode::physicalField,variable.ordinal,component,0,0,0,1});
          if(component && variable.ordinal==static_cast<std::uint32_t>(WVPortableVariable::energy))
            keys.push_back({WVVariableEvaluationNode::reduction,variable.ordinal,component});
        }
        keys.push_back({WVVariableEvaluationNode::reduction,variable.ordinal});
        for(std::uint32_t reference=0;reference<2;++reference)
          keys.push_back({WVVariableEvaluationNode::registeredVariable,variable.ordinal,0,0,reference});
      }
      for(std::uint32_t component=0;component<5;++component)
        for(std::uint32_t family=0;family<3;++family)
          keys.push_back({WVVariableEvaluationNode::componentCoefficients,family,component});
      for(std::uint32_t field=0;field<16;++field)
        for(std::uint32_t component=0;component<5;++component)
          for(std::uint32_t derivative=1;derivative<4;++derivative)
            keys.push_back({WVVariableEvaluationNode::derivative,field,
                component,derivative});
      keys.push_back({WVVariableEvaluationNode::phaseFactors});
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to prepare output variable evaluation keys."};
    }
    return context.prepare(keys);
  }
  explicit WVFieldEvaluationEventWorkspace(const WVIntegrationState& state,
      WVVariableEvaluationContext& evaluation,WVFieldEvaluationArena& arena,
      WVVariableEvaluationPolicy policy=WVVariableEvaluationPolicy::reuse,
      bool enabled=true)
      : t_(state.waveVortex.t), t0_(state.waveVortex.t0),
        realFields_(arena.realFields),complexFields_(arena.complexFields),
        groupNodes_(arena.groupNodes),
        evaluation_(evaluation),arena_(arena),
        evaluationMetricsBefore_(evaluation.metrics()),
        densitySource_(arena.densitySource),
        densityHeights_(arena.densityHeights),
        densityWeights_(arena.densityWeights),
        densityInitial_(arena.densityInitial),density_(arena.density),
        apv_(arena.apv) {
    coefficients_={state.waveVortex.coefficients.Ap.data,state.waveVortex.coefficients.Am.data,state.waveVortex.coefficients.A0.data};
    if(state.coefficientFamilies) {
      if(state.coefficientFamilyCount==1) coefficients_={nullptr,nullptr,state.coefficientFamilies[0].data};
      else if(state.coefficientFamilyCount==3)
        for(std::size_t family=0;family<3;++family) coefficients_[family]=state.coefficientFamilies[family].data;
    }
    for(std::size_t reference=0;reference<density_.size();++reference)
      densityMetrics_[reference]=density_[reference].metrics();
    if(!enabled) return;
    prepareStatus_=evaluation_.begin(this,policy);
  }

  const WVKernelStatus& status() const noexcept {return prepareStatus_;}
  WVVariableEvaluationMetrics evaluationMetrics() const noexcept {
    const auto& now=evaluation_.metrics();
    return {now.contexts-evaluationMetricsBefore_.contexts,
        now.producerExecutions-evaluationMetricsBefore_.producerExecutions,
        now.cacheHits-evaluationMetricsBefore_.cacheHits,
        now.evictions-evaluationMetricsBefore_.evictions,
        now.recomputations-evaluationMetricsBefore_.recomputations,
        now.duplicateExecutions-evaluationMetricsBefore_.duplicateExecutions,
        now.liveBytes,variableHighWaterBytes_};
  }
  bool densityRecoveryReport(WVNoMotionRecoveryReport& report) const noexcept {
    if(density_[0].metrics().recoveryAttemptCount==0) return false;
    report=density_[0].recoveryReport();
    return true;
  }
  WVVariableEvaluationPolicy policy() const noexcept {return evaluation_.policy();}
  std::uint64_t generation() const noexcept {return evaluation_.generation();}
  bool preparationIdle() const noexcept {return groupNodes_.empty();}
  bool ready(const WVVariableEvaluationKey& key) const noexcept {
    return evaluation_.ready(key);
  }
  WVKernelStatus forcingWorkspace(const WVForcingDiagnosticBinding& binding,
      WVForcingDiagnosticWorkspace*& result) {
    if(!forcingWorkspace_) {
      forcingWorkspace_=arena_.forcingWorkspace.get();
      if(!forcingWorkspace_)
        return {WVKernelStatusCode::invalidConfiguration,
            "The output-event forcing arena was not prepared."};
      const auto status=forcingWorkspace_->beginScopedEvaluation(
          evaluation_,binding.stages());
      if(!status) {forcingWorkspace_=nullptr; return status;}
      forcingWorkspace_->scalarEvaluationOwner=this;
      forcingWorkspace_->horizontalMaximumEvaluator=[](void* context,
          double& result,void* producerOwner,
          WVForcingDiagnosticWorkspace::ScalarProducer producer) {
        auto& workspace=*static_cast<WVFieldEvaluationEventWorkspace*>(context);
        bool reused=false;
        const WVVariableEvaluationKey key{WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(WVPortableVariable::uvMax)};
        return workspace.evaluate(key,&result,1,
            [&]() {return producer(producerOwner,result);},reused);
      };
      forcingWorkspace_->derivativeAccess.context=this;
      forcingWorkspace_->derivativeAccess.lookup=[](void* context,
          std::size_t field,std::size_t derivative,
          WVRealVolumeConstView& result) {
        return static_cast<WVFieldEvaluationEventWorkspace*>(context)->
            lookupDerivative(field,derivative,result);
      };
      forcingWorkspace_->derivativeAccess.capture=[](void* context,
          std::size_t field,std::size_t derivative,
          WVRealVolumeConstView value) {
        return static_cast<WVFieldEvaluationEventWorkspace*>(context)->
            captureDerivative(field,derivative,value);
      };
      account();
    }
    result=forcingWorkspace_;
    return WVKernelStatus::ok();
  }
  void endForcingEvaluation() noexcept {
    if(forcingWorkspace_) forcingWorkspace_->endScopedEvaluation();
    forcingWorkspace_=nullptr;
  }
  void endEvaluation() noexcept {evaluation_.end();}
  std::uint32_t component() const noexcept {return component_;}
  void setComponent(std::uint32_t component,
      const WVIntegrationState* state=nullptr) noexcept {
    component_=component;
    componentCoefficients_={};
    if(component && state)
      componentCoefficients_={state->waveVortex.coefficients.Ap.data,
          state->waveVortex.coefficients.Am.data,
          state->waveVortex.coefficients.A0.data};
  }

  WVKernelStatus validateState(const WVIntegrationState& state) const {
    std::array<const WVComplex64*,3> coefficients{
        state.waveVortex.coefficients.Ap.data,
        state.waveVortex.coefficients.Am.data,
        state.waveVortex.coefficients.A0.data};
    if(state.coefficientFamilies) {
      if(state.coefficientFamilyCount==1)
        coefficients={nullptr,nullptr,state.coefficientFamilies[0].data};
      else if(state.coefficientFamilyCount==3)
        for(std::size_t family=0;family<3;++family)
          coefficients[family]=state.coefficientFamilies[family].data;
    }
    const auto& expected=component_ ? componentCoefficients_ : coefficients_;
    if(state.waveVortex.t!=t_ || state.waveVortex.t0!=t0_ || coefficients!=expected)
      return {WVKernelStatusCode::invalidConfiguration,
              "Field evaluation belongs to a different active event state."};
    return WVKernelStatus::ok();
  }

  std::size_t externalWorkspaceBytes() const noexcept {
    return externalWorkspaceBytes_;
  }
  void setExternalWorkspaceBytes(std::size_t bytes) noexcept {
    externalWorkspaceBytes_=bytes;
    if(metrics_) account();
  }

  template<class Operation>
  WVKernelStatus evaluate(const WVVariableEvaluationKey& key,double* output,
      std::size_t count,Operation&& operation,bool& reused) {
    reused=false;
    if(!retainPrimitiveFields_)
      return operation();
    if(policy()==WVVariableEvaluationPolicy::lowMemory) {
      const auto status=evaluation_.evaluate(key,0,std::forward<Operation>(operation));
      if(!status) return status;
      noteVariableBytes();
      (void)evaluation_.evict(key);
      account();
      return WVKernelStatus::ok();
    }
    const auto firstPhysical=static_cast<std::uint32_t>(WVPortableVariable::u);
    const auto lastPhysical=static_cast<std::uint32_t>(WVPortableVariable::eta);
    if(key.node==WVVariableEvaluationNode::physicalField &&
        key.variable>=firstPhysical && key.variable<=lastPhysical &&
        key.derivative==0 && key.reference==0 && key.forcingOrdinal==0 &&
        key.stage==0 && key.geometry==0) {
      const WVVariableEvaluationKey canonical{
          WVVariableEvaluationNode::physicalField,firstPhysical,key.component,
          0,0,0,1};
      if(evaluation_.ready(canonical)) {
        auto* values=findReal(canonical);
        if(!values || values->size()%4!=0 || values->size()/4!=count)
          return {WVKernelStatusCode::invalidConfiguration,
              "A cached canonical physical bundle has a different extent."};
        const auto status=evaluation_.evaluate(canonical,
            values->capacity()*sizeof(double),
            [](){return WVKernelStatus::ok();});
        if(!status) return status;
        const auto channel=key.variable-firstPhysical;
        std::copy_n(values->data()+channel*count,count,output);
        reused=true;
        ++metrics_->eventFieldReuseCount;
        noteVariableBytes();
        return WVKernelStatus::ok();
      }
    }
    auto* cached=findReal(key);
    if(evaluation_.ready(key)) {
      if(!cached || cached->size()!=count)
        return {WVKernelStatusCode::invalidConfiguration,
            "A cached output variable was requested with a different extent."};
      const auto status=evaluation_.evaluate(key,cached->capacity()*sizeof(double),[](){return WVKernelStatus::ok();});
      if(!status) return status;
      const auto& values=*cached;
      std::copy_n(values.data(),count,output);
      reused=true;
      ++metrics_->eventFieldReuseCount;
      noteVariableBytes();
      return WVKernelStatus::ok();
    }
    auto& prepared=real(key,count);
    try {prepared.reserve(count);} catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to reserve a reconstructed output field."};
    }
    const auto bytes=prepared.capacity()*sizeof(double);
    const auto status=evaluation_.evaluate(key,bytes,[&]() {
      const auto produced=operation();
      if(!produced) return produced;
      auto* values=findReal(key);
      if(!values)
        return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
            "Output variable storage changed while its producer was running."};
      try {values->assign(output,output+count);}
      catch(const std::bad_alloc&) {
        return WVKernelStatus{WVKernelStatusCode::allocationFailure,
          "Unable to retain a reconstructed field for this output event."};
      }
      return WVKernelStatus::ok();
    });
    if(!status) return status;
    noteVariableBytes();
    account();
    if(evaluation_.evict(key)) if(auto* values=findReal(key)) std::vector<double>{}.swap(*values);
    return WVKernelStatus::ok();
  }

  template<class Operation>
  WVKernelStatus evaluate(std::size_t field,const WVState&,double* output,
      std::size_t count,Operation&& operation,bool& reused) {
    return evaluate({field>=dynamicalDerivativeKeyBase ? WVVariableEvaluationNode::derivative : WVVariableEvaluationNode::physicalField,
        static_cast<std::uint32_t>(field>=dynamicalDerivativeKeyBase ? field-dynamicalDerivativeKeyBase : field),
        component_,field>=dynamicalDerivativeKeyBase ? 1u : 0u},output,count,
        std::forward<Operation>(operation),reused);
  }

  template<class Keys,class Outputs,class Operation>
  WVKernelStatus evaluateGroup(const Keys& keys,const Outputs& outputs,Operation&& operation,
      bool& reused) {
    reused=false;
    if(keys.size()!=outputs.size() || keys.empty())
      return {WVKernelStatusCode::invalidConfiguration,
          "A fused output producer requires one nonempty output per key."};
    if(!groupNodes_.empty() || keys.size()>groupNodes_.capacity())
      return {WVKernelStatusCode::invalidConfiguration,
          "Fused output-event metadata was not prepared for this workload."};
    struct ClearNodes {
      std::vector<std::pair<WVVariableEvaluationKey,std::size_t>>& nodes;
      ~ClearNodes() {nodes.clear();}
    } clearNodes{groupNodes_};
    if(policy()==WVVariableEvaluationPolicy::lowMemory) {
      for(const auto& key:keys) groupNodes_.push_back({key,0});
      const auto status=evaluation_.evaluateGroup(groupNodes_,
          std::forward<Operation>(operation));
      if(!status) return status;
      noteVariableBytes();
      for(const auto& key:keys) (void)evaluation_.evict(key);
      account();
      return WVKernelStatus::ok();
    }
    try {
      for(std::size_t index=0;index<keys.size();++index) {
        auto& values=real(keys[index],outputs[index].elementCount);
        values.reserve(outputs[index].elementCount);
        groupNodes_.push_back({keys[index],values.capacity()*sizeof(double)});
      }
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,
          "Unable to reserve fused output variable storage."};
    }
    bool produced=false;
    const auto status=evaluation_.evaluateGroup(groupNodes_,[&]() {
      const auto result=operation();
      if(!result) return result;
      for(std::size_t index=0;index<keys.size();++index) {
        auto* values=findReal(keys[index]);
        if(!values) return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
            "Fused output storage changed while its producer was running."};
        try {values->assign(outputs[index].data,
            outputs[index].data+outputs[index].elementCount);}
        catch(const std::bad_alloc&) {
          return WVKernelStatus{WVKernelStatusCode::allocationFailure,
              "Unable to retain a fused output variable."};
        }
      }
      produced=true;
      return WVKernelStatus::ok();
    });
    if(!status) return status;
    noteVariableBytes();
    if(!produced) {
      for(std::size_t index=0;index<keys.size();++index) {
        const auto* values=findReal(keys[index]);
        if(!values || values->size()!=outputs[index].elementCount)
          return {WVKernelStatusCode::invalidConfiguration,
              "A cached fused output variable has a different extent."};
        std::copy(values->begin(),values->end(),outputs[index].data);
      }
      reused=true;
      ++metrics_->eventFieldReuseCount;
    }
    account();
    for(const auto& key:keys)
      if(evaluation_.evict(key))
        if(auto* values=findReal(key)) std::vector<double>{}.swap(*values);
    return WVKernelStatus::ok();
  }

  template<class Operation>
  WVKernelStatus complexView(const WVVariableEvaluationKey& key,std::size_t count,
      Operation&& operation,const std::vector<WVComplex64>*& view) {
    auto* cached=findComplex(key);
    if(evaluation_.ready(key)) {
      if(!cached || cached->size()!=count)
        return {WVKernelStatusCode::invalidConfiguration,
            "A cached complex output variable was requested with a different extent."};
      const auto status=evaluation_.evaluate(key,cached->capacity()*sizeof(WVComplex64),
          [](){return WVKernelStatus::ok();});
      if(!status) return status;
      view=cached;
      evaluation_.pin(key);
      noteVariableBytes();
      return WVKernelStatus::ok();
    }
    auto& prepared=complex(key);
    try {prepared.resize(count);} catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,"Unable to retain an event complex field."};
    }
    const auto bytes=prepared.capacity()*sizeof(WVComplex64);
    auto* destination=prepared.data();
    const auto status=evaluation_.evaluate(key,bytes,[&]() {
      return operation(destination);
    });
    if(!status) return status;
    view=findComplex(key);
    if(!view) return {WVKernelStatusCode::invalidConfiguration,
        "Complex output variable storage changed while its producer was running."};
    evaluation_.pin(key); noteVariableBytes(); account(); return WVKernelStatus::ok();
  }
  void releaseComplex(const WVVariableEvaluationKey& key) noexcept {
    evaluation_.unpin(key);
    if(evaluation_.evict(key)) for(auto& entry:complexFields_)
      if(entry->assigned && entry->key==key) {
        if(entry->prepared) entry->values.clear();
        else std::vector<WVComplex64>{}.swap(entry->values);
        break;
      }
    account();
  }
  WVKernelStatus checkoutScratch(const WVVariableEvaluationKey& key,
      std::size_t count,std::vector<double>& storage) {
    if(!storage.empty())
      return {WVKernelStatusCode::invalidConfiguration,
          "Output-event scratch destination is already in use."};
    RealEntry* selected=nullptr;
    for(auto& entry:realFields_) if(!entry->assigned && entry->prepared &&
        entry->preparedKey==key) {selected=entry.get(); break;}
    if(!selected) {
      const auto matching=std::count_if(realFields_.begin(),realFields_.end(),
          [&](const auto& entry) {return entry->prepared && entry->preparedKey==key;});
      return {WVKernelStatusCode::invalidConfiguration,
          "Output-event scratch storage was not prepared or remained in use: stage="+
              std::to_string(key.stage)+" group="+std::to_string(key.component)+
              " slot="+std::to_string(key.variable)+" matching="+
              std::to_string(matching)+"."};
    }
    selected->key=key; selected->assigned=true;
    storage.swap(selected->values);
    try {storage.resize(count);} catch(const std::bad_alloc&) {
      storage.swap(selected->values); selected->assigned=false;
      return {WVKernelStatusCode::allocationFailure,
          "Unable to size prepared output-event scratch storage."};
    }
    account();
    return WVKernelStatus::ok();
  }
  void returnScratch(const WVVariableEvaluationKey& key,
      std::vector<double>& storage) noexcept {
    for(auto& entry:realFields_) if(entry->assigned && entry->key==key) {
      storage.clear(); storage.swap(entry->values); entry->assigned=false;
      account(); return;
    }
  }

  WVKernelStatus lookupDerivative(std::size_t field,std::size_t derivative,
      WVRealVolumeConstView& result) {
    result={};
    static constexpr WVPortableVariable fields[]={WVPortableVariable::u,
        WVPortableVariable::v,WVPortableVariable::w,WVPortableVariable::eta};
    if(field>=std::size(fields))
      return {WVKernelStatusCode::invalidConfiguration,
          "A state derivative requested an unknown field."};
    const WVVariableEvaluationKey key{WVVariableEvaluationNode::derivative,
        static_cast<std::uint32_t>(fields[field]),component_,
        static_cast<std::uint32_t>(derivative)};
    if(!evaluation_.ready(key)) return WVKernelStatus::ok();
    auto* values=findReal(key);
    auto* entry=findRealEntry(key);
    if(!values || !entry || values->size()!=entry->volumeShape.elementCount())
      return {WVKernelStatusCode::invalidConfiguration,
          "A cached state derivative has incompatible storage."};
    const auto status=evaluation_.evaluate(key,values->capacity()*sizeof(double),
        [](){return WVKernelStatus::ok();});
    if(status) result={values->data(),entry->volumeShape};
    return status;
  }
  WVKernelStatus captureDerivative(std::size_t field,std::size_t derivative,
      WVRealVolumeConstView value) {
    static constexpr WVPortableVariable fields[]={WVPortableVariable::u,
        WVPortableVariable::v,WVPortableVariable::w,WVPortableVariable::eta};
    if(field>=std::size(fields))
      return {WVKernelStatusCode::invalidConfiguration,
          "A state derivative captured an unknown field."};
    const WVVariableEvaluationKey key{WVVariableEvaluationNode::derivative,
        static_cast<std::uint32_t>(fields[field]),component_,
        static_cast<std::uint32_t>(derivative)};
    if(!arena_.hasPreparedReal(key)) return WVKernelStatus::ok();
    if(!value.data)
      return {WVKernelStatusCode::invalidPointer,
          "A captured state derivative has null storage."};
    if(policy()==WVVariableEvaluationPolicy::lowMemory) {
      const auto status=evaluation_.evaluate(key,0,[](){return WVKernelStatus::ok();});
      if(status) (void)evaluation_.evict(key);
      return status;
    }
    auto& values=real(key,value.shape.elementCount());
    auto* entry=findRealEntry(key);
    const auto status=evaluation_.evaluate(key,
        values.capacity()*sizeof(double),[&]() {
          values.assign(value.data,value.data+value.shape.elementCount());
          entry->volumeShape=value.shape;
          return WVKernelStatus::ok();
        });
    if(status) {noteVariableBytes(); account();}
    return status;
  }

  template<class Operation>
  WVKernelStatus componentCoefficients(std::uint32_t component,std::uint32_t family,
      std::size_t count,Operation&& operation,const WVComplex64*& view) {
    const WVVariableEvaluationKey key{WVVariableEvaluationNode::componentCoefficients,
      family,component};
    const std::vector<WVComplex64>* values=nullptr;
    const auto status=complexView(key,count,[&](WVComplex64* output) {
      operation(output); return WVKernelStatus::ok();
    },values);
    if(status) view=values->data();
    return status;
  }
  void releaseComponentCoefficients(std::uint32_t component,
      std::uint32_t family) noexcept {
    releaseComplex({WVVariableEvaluationNode::componentCoefficients,family,component});
  }

  WVKernelStatus validateDensityBinding(const WVState& state,
      WVDensityDiagnosticContract) const {
    const std::array<const WVComplex64*,3> coefficients{state.coefficients.Ap.data,state.coefficients.Am.data,state.coefficients.A0.data};
    if(state.t!=t_ || state.t0!=t0_ || coefficients!=coefficients_)
      return {WVKernelStatusCode::invalidConfiguration,"Density evaluation belongs to a different active event state."};
    return WVKernelStatus::ok();
  }
  bool hasDensitySource() const noexcept {return !densitySource_.empty();}
  WVKernelStatus bindDensity(std::vector<double>* source,WVShape3D shape,
      WVDensityEventGeometry geometry,WVDensityDiagnosticContract contract,bool preserveSource=false) {
    if(densitySource_.empty()) {
      if(!source)
        return {WVKernelStatusCode::invalidConfiguration,
            "Density source storage was not prepared for this event."};
      if(!geometry.heights || !geometry.integrationWeights ||
          !geometry.initialProfile)
        return {WVKernelStatusCode::invalidPointer,
            "Density event geometry is incomplete."};
      bool movedSource=false;
      try {
        // Complete every potentially allocating geometry copy before moving
        // the produced source. Both source buffers belong to the prepared
        // arena, so swapping also preserves their capacities for later events.
        densityHeights_=*geometry.heights;
        densityWeights_=*geometry.integrationWeights;
        densityInitial_=*geometry.initialProfile;
        if(preserveSource) densitySource_=*source;
        else {densitySource_.swap(*source); movedSource=true;}
        densityShape_=shape;
        densityGeometry_=geometry;
        densityGeometry_.heights=&densityHeights_;
        densityGeometry_.integrationWeights=&densityWeights_;
        densityGeometry_.initialProfile=&densityInitial_;
      } catch(...) {
        if(movedSource) densitySource_.swap(*source);
        else densitySource_.clear();
        densityHeights_.clear();
        densityWeights_.clear();
        densityInitial_.clear();
        densityShape_={};
        densityGeometry_={};
        account();
        throw;
      }
    }
    WVKernelStatus status=WVKernelStatus::ok();
    for(const auto reference:{WVNoMotionReference::actual,WVNoMotionReference::initial}) {
      const auto selected=densityIndex(reference);
      if(density_[selected].initialized()) continue;
      auto selectedContract=contract;
      selectedContract.reference=reference;
      status=WVDensityEventEvaluation::create(
          {densitySource_.data(),densityShape_},densityGeometry_,selectedContract,density_[selected]);
      if(!status) break;
    }
    account();
    return status;
  }
  WVKernelStatus prepareDensity(std::uint8_t demands,
      WVDensityDiagnosticContract contract) {
    WVKernelStatus status=WVKernelStatus::ok();
    const struct {std::uint8_t demand; WVPortableVariable variable;} requests[]={{
        WVDensityEventEvaluation::rhoNmDemand,WVPortableVariable::rho_nm},{
        WVDensityEventEvaluation::etaTrueDemand,WVPortableVariable::eta_true},{
        WVDensityEventEvaluation::apeDemand,WVPortableVariable::ape}};
    for(const auto& request:requests) if(demands&request.demand) {
      const auto count=request.variable==WVPortableVariable::rho_nm ?
          densityHeights_.size() : densitySource_.size();
      const auto selected=request.variable==WVPortableVariable::rho_nm ? 0u :
          densityIndex(contract.reference);
      const WVVariableEvaluationKey key{WVVariableEvaluationNode::registeredVariable,
          static_cast<std::uint32_t>(request.variable),0,0,
          request.variable==WVPortableVariable::rho_nm ? 0u :
              static_cast<std::uint32_t>(contract.reference)};
      const bool reused=evaluation_.ready(key);
      status=evaluation_.evaluate(key,count*sizeof(double),[&]() {
        return density_[selected].prepare(request.demand);
      });
      if(!status) break;
      if(reused) ++metrics_->densityReuseCount;
      noteVariableBytes();
    }
    for(std::size_t reference=0;reference<density_.size();++reference) {
      const auto& now=density_[reference].metrics();
      metrics_->densityRecoveryCount+=now.recoveryCount-densityMetrics_[reference].recoveryCount;
      metrics_->densityProfileConstructionCount+=now.profileConstructionCount-densityMetrics_[reference].profileConstructionCount;
      metrics_->densityInversePassCount+=now.inversePassCount-densityMetrics_[reference].inversePassCount;
      metrics_->densityAPEPassCount+=now.apePassCount-densityMetrics_[reference].apePassCount;
      metrics_->densityReuseCount+=now.reuseCount-densityMetrics_[reference].reuseCount;
      densityMetrics_[reference]=now;
    }
    account();
    return status;
  }
  void finishDensityUse(WVDensityDiagnosticContract contract) noexcept {
    if(policy()!=WVVariableEvaluationPolicy::lowMemory) return;
    for(const auto variable:{WVPortableVariable::rho_nm,WVPortableVariable::eta_true,
                             WVPortableVariable::ape})
      evaluation_.evict({WVVariableEvaluationNode::registeredVariable,
          static_cast<std::uint32_t>(variable),0,0,
          variable==WVPortableVariable::rho_nm ? 0u :
              static_cast<std::uint32_t>(contract.reference)});
    for(auto& density:density_) density.discardDerived();
    releaseAPV();
    account();
  }
  WVDensityEventView densityView(WVDensityEventField field,
      WVNoMotionReference reference) const noexcept {
    const auto selected=field==WVDensityEventField::rhoNm ? 0u : densityIndex(reference);
    return density_[selected].view(field);
  }
  bool hasAPV(WVNoMotionReference reference) const noexcept {return apvReady_[densityIndex(reference)];}
  WVDensityEventView apvView(WVNoMotionReference reference) const noexcept {
    const auto selected=densityIndex(reference);
    return apvReady_[selected] ? WVDensityEventView{apv_[selected].data(),apv_[selected].size()} : WVDensityEventView{};
  }
  template<class Operation>
  WVKernelStatus prepareAPV(std::size_t count,WVNoMotionReference reference,
      Operation&& operation) {
    const auto selected=densityIndex(reference);
    if(apvReady_[selected]) {++metrics_->densityAPVReuseCount; return WVKernelStatus::ok();}
    apv_[selected].resize(count); account();
    const auto status=operation(density_[selected].view(WVDensityEventField::etaTrue),apv_[selected].data());
    if(!status) return status;
    apvReady_[selected]=true;
    ++metrics_->densityAPVPassCount;
    return WVKernelStatus::ok();
  }
  void releaseAPV(WVNoMotionReference reference) noexcept {
    const auto selected=densityIndex(reference);
    apv_[selected].clear(); apvReady_[selected]=false;
  }
  void releaseAPV() noexcept {
    releaseAPV(WVNoMotionReference::actual);
    releaseAPV(WVNoMotionReference::initial);
  }


private:
  friend class WVFieldEvaluationEventScope;
  double t_=0,t0_=0;
  std::array<const WVComplex64*,3> coefficients_{};
  std::array<const WVComplex64*,3> componentCoefficients_{};
  using RealEntry=WVFieldEvaluationRealEntry;
  using ComplexEntry=WVFieldEvaluationComplexEntry;
  std::vector<std::unique_ptr<RealEntry>>& realFields_;
  std::vector<std::unique_ptr<ComplexEntry>>& complexFields_;
  std::vector<std::pair<WVVariableEvaluationKey,std::size_t>>& groupNodes_;
  WVVariableEvaluationContext& evaluation_;
  WVFieldEvaluationArena& arena_;
  WVVariableEvaluationMetrics evaluationMetricsBefore_;
  std::size_t variableHighWaterBytes_=0;
  WVKernelStatus prepareStatus_=WVKernelStatus::ok();
  std::uint32_t component_=0;
  WVFieldEvaluationMetrics* metrics_=nullptr;
  bool retainPrimitiveFields_=true;
  std::vector<double>& densitySource_;
  std::vector<double>& densityHeights_;
  std::vector<double>& densityWeights_;
  std::vector<double>& densityInitial_;
  WVShape3D densityShape_{};
  WVDensityEventGeometry densityGeometry_{};
  std::array<WVDensityEventEvaluation,2>& density_;
  std::array<WVDensityEventMetrics,2> densityMetrics_;
  std::array<std::vector<double>,2>& apv_;
  std::array<bool,2> apvReady_{};
  std::size_t externalWorkspaceBytes_=0;
  WVForcingDiagnosticWorkspace* forcingWorkspace_=nullptr;
  static constexpr std::size_t densityIndex(WVNoMotionReference reference) noexcept {
    return reference==WVNoMotionReference::initial ? 1u : 0u;
  }
  void noteVariableBytes() noexcept {
    variableHighWaterBytes_=std::max(variableHighWaterBytes_,
        evaluation_.metrics().liveBytes);
  }
  void account() noexcept {
    const auto sourceBytes=(densitySource_.capacity()+densityHeights_.capacity()+
        densityWeights_.capacity()+densityInitial_.capacity()+apv_[0].capacity()+
        apv_[1].capacity())*sizeof(double);
    const auto live=sourceBytes+density_[0].metrics().liveBytes+
        density_[1].metrics().liveBytes;
    const auto peak=sourceBytes+density_[0].metrics().highWaterBytes+
        density_[1].metrics().highWaterBytes;
    metrics_->densityWorkspaceLiveBytes=live;
    metrics_->densityWorkspaceHighWaterBytes=std::max(metrics_->densityWorkspaceHighWaterBytes,peak);
    arena_.notePeak();
    metrics_->eventFieldArenaPlannedBytes=arena_.plannedBytes;
    metrics_->eventFieldArenaPeakBytes=std::max(
        metrics_->eventFieldArenaPeakBytes,arena_.peakBytes);
    const auto logicalVariableBytes=evaluation_.metrics().liveBytes;
    metrics_->eventFieldWorkspaceLiveBytes=externalWorkspaceBytes_+live+
        logicalVariableBytes;
    metrics_->eventFieldWorkspaceHighWaterBytes=std::max(
        metrics_->eventFieldWorkspaceHighWaterBytes,
        externalWorkspaceBytes_+peak+variableHighWaterBytes_);
    metrics_->additionalTransientHighWaterBytes=std::max(
        metrics_->additionalTransientHighWaterBytes,externalWorkspaceBytes_);
  }
  std::vector<double>* findReal(const WVVariableEvaluationKey& key) noexcept {
    for(auto& entry:realFields_) if(entry->assigned && entry->key==key) return &entry->values;
    return nullptr;
  }
  RealEntry* findRealEntry(const WVVariableEvaluationKey& key) noexcept {
    for(auto& entry:realFields_) if(entry->assigned && entry->key==key)
      return entry.get();
    return nullptr;
  }
  std::vector<WVComplex64>* findComplex(const WVVariableEvaluationKey& key) noexcept {
    for(auto& entry:complexFields_) if(entry->assigned && entry->key==key) return &entry->values;
    return nullptr;
  }
  std::vector<double>& real(const WVVariableEvaluationKey& key,
      std::size_t count) {
    if(auto* values=findReal(key)) return *values;
    for(auto& entry:realFields_) if(!entry->assigned && entry->prepared &&
        entry->preparedKey==key) {
      entry->key=key; entry->assigned=true; return entry->values;
    }
    RealEntry* best=nullptr;
    for(auto& entry:realFields_) if(!entry->assigned && !entry->prepared &&
        entry->values.capacity()>=count &&
        (!best || entry->values.capacity()<best->values.capacity())) {
      best=entry.get();
    }
    if(best) {
      best->key=key;
      best->assigned=true;
      return best->values;
    }
    for(auto& entry:realFields_) if(!entry->assigned && !entry->prepared) {
      entry->key=key; entry->assigned=true; return entry->values;
    }
    auto entry=std::make_unique<RealEntry>();
    entry->key=key; entry->assigned=true;
    realFields_.push_back(std::move(entry));
    return realFields_.back()->values;
  }
  std::vector<WVComplex64>& complex(const WVVariableEvaluationKey& key) {
    if(auto* values=findComplex(key)) return *values;
    for(auto& entry:complexFields_) if(!entry->assigned && entry->prepared &&
        entry->preparedKey==key) {
      entry->key=key; entry->assigned=true; return entry->values;
    }
    for(auto& entry:complexFields_) if(!entry->assigned && !entry->prepared) {
      entry->key=key; entry->assigned=true; return entry->values;
    }
    auto entry=std::make_unique<ComplexEntry>();
    entry->key=key; entry->assigned=true;
    complexFields_.push_back(std::move(entry));
    return complexFields_.back()->values;
  }
};

class WVFieldEvaluationEventScope final {
public:
  WVFieldEvaluationEventScope(WVFieldEvaluationService&,const WVIntegrationState&,bool enabled=true,bool retainPrimitiveFields=true);
  ~WVFieldEvaluationEventScope() {release();}
  WVFieldEvaluationEventScope(const WVFieldEvaluationEventScope&)=delete;
  WVFieldEvaluationEventScope& operator=(const WVFieldEvaluationEventScope&)=delete;
  const WVKernelStatus& status() const noexcept {return status_;}
  void setExternalWorkspaceBytes(std::size_t bytes) noexcept {
    workspace_.setExternalWorkspaceBytes(bytes);
  }
  void release() noexcept;
private:
  WVFieldEvaluationService* service_=nullptr;
  WVFieldEvaluationEventWorkspace workspace_;
  WVKernelStatus status_=WVKernelStatus::ok();
  bool kernelEvaluationActive_=false;
};

} // namespace wavevortex::runtime::detail
