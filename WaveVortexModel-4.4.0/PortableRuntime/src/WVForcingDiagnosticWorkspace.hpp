#pragma once

#include "WaveVortexRuntime/WVForcing.hpp"
#include "WaveVortexKernel/WVVariableExecutionOptions.hpp"
#include "WaveVortexRuntime/WVForcingTendency.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace wavevortex::runtime::detail {

enum class WVGridCalculusKind : std::size_t {
  horizontalSecondX,
  horizontalSecondY,
  verticalSecondF,
  verticalSecondG
};

struct WVStateCalculusAccess {
  void* context=nullptr;
  WVKernelStatus (*lookup)(void*,std::size_t,WVGridCalculusKind,
      WVRealVolumeConstView&)=nullptr;
  WVKernelStatus (*capture)(void*,std::size_t,WVGridCalculusKind,
      WVRealVolumeConstView)=nullptr;
};

// Standalone diagnostic calls own a fresh ledger; explicit event sessions use
// their enclosing ledger instead. Retain only cumulative metrics after return.
class WVForcingDiagnosticLedger final {
public:
  explicit WVForcingDiagnosticLedger(WVForcingTendencyMetrics& metrics):metrics_(metrics) {}
  ~WVForcingDiagnosticLedger() {
    if(!context.active()) return;
    const auto value=context.metrics();
    auto& sum=metrics_.standaloneVariableEvaluation;
    sum.contexts+=value.contexts; sum.producerExecutions+=value.producerExecutions;
    sum.cacheHits+=value.cacheHits; sum.evictions+=value.evictions;
    sum.recomputations+=value.recomputations; sum.duplicateExecutions+=value.duplicateExecutions;
    sum.liveBytes=0; sum.highWaterBytes=std::max(sum.highWaterBytes,value.highWaterBytes);
    context.end();
  }
  WVVariableEvaluationContext context;
private:
  WVForcingTendencyMetrics& metrics_;
};

// Owned by one diagnostic event (or one standalone invocation).
// The ordinary RHS has no dependency on these state-sized arrays.
class WVForcingDiagnosticWorkspace final {
public:
  WVForcingDiagnosticWorkspace(WVShape2D spectral, WVShape4D spatial,
      std::size_t coefficientFamilies=3,std::size_t physicalChannels=4)
      : spectral(spectral), spatial(spatial), flux(coefficientFamilies*spectral.elementCount()),
        previous(flux.size()), temporary(flux.size()),
        cumulative(spatial.elementCount()), raw(spatial.elementCount()),
        physical(physicalChannels*spatial.first*spatial.second*spatial.third),
        coefficientElements_(flux.size()),physicalElements_(physical.size()) {}

  static WVVariableEvaluationKey prefixKey(std::size_t index) {
    return {WVVariableEvaluationNode::forcingTendency,0,0,0,0,
        static_cast<std::uint32_t>(index),2};
  }
  static WVVariableEvaluationKey outputKey(std::size_t index) {
    return {WVVariableEvaluationNode::forcingTendency,0,0,0,0,
        static_cast<std::uint32_t>(index),3};
  }
  static WVVariableEvaluationKey projectionKey() {
    return {WVVariableEvaluationNode::forcingTendency,0,0,0,0,0,4};
  }
  static WVVariableEvaluationKey nonlinearKey() {
    return {WVVariableEvaluationNode::forcingTendency,0,0,0,0,0,6};
  }
  static WVVariableEvaluationKey gridCalculusKey(
      std::size_t field,WVGridCalculusKind kind) {
    return {WVVariableEvaluationNode::gridCalculus,
        static_cast<std::uint32_t>(field),0,
        static_cast<std::uint32_t>(kind)};
  }
  static WVVariableEvaluationKey constantLaplacianKey(std::size_t direction) {
    return {WVVariableEvaluationNode::gridCalculus,16,0,
        static_cast<std::uint32_t>(direction)};
  }
  std::size_t nonlinearUseCount=1;
  std::array<std::size_t,16> gridCalculusUseCount{};
  std::array<std::size_t,2> constantLaplacianUseCount{};
  WVStateDerivativeAccess derivativeAccess;
  WVStateDerivativeAccess* stateDerivativeAccess() noexcept {
    return derivativeAccess.lookup || derivativeAccess.capture ? &derivativeAccess : nullptr;
  }
  WVStateCalculusAccess* stateCalculusAccess() noexcept {
    return std::any_of(gridCalculusUseCount.begin(),gridCalculusUseCount.end(),
        [](std::size_t count){return count>1;}) ? &calculusAccess_ : nullptr;
  }
  template<class Producer>
  WVKernelStatus evaluateNonlinearRaw(Producer&& producer) {
    if(!evaluation_ || (reusesPrefix() && nonlinearUseCount<2)) return producer();
    const bool retain=reusesPrefix();
    const auto status=evaluation_->evaluate(nonlinearKey(),
        retain ? nonlinearRaw_.size()*sizeof(double) : 0,[&] {
      auto result=producer();
      if(result && retain) std::copy(raw.begin(),raw.end(),nonlinearRaw_.begin());
      return result;
    });
    if(!status) return status;
    if(retain) std::copy(nonlinearRaw_.begin(),nonlinearRaw_.end(),raw.begin());
    else evaluation_->evict(nonlinearKey());
    return status;
  }
  static WVVariableEvaluationKey horizontalMaximumKey() {
    return {WVVariableEvaluationNode::reduction,0,0,0,0,0,5};
  }
  using ScalarProducer=WVKernelStatus(*)(void*,double&);
  using ScalarEvaluator=WVKernelStatus(*)(void*,double&,void*,ScalarProducer);
  void* scalarEvaluationOwner=nullptr;
  ScalarEvaluator horizontalMaximumEvaluator=nullptr;
  template<class Producer>
  WVKernelStatus evaluateHorizontalMaximum(double& result,Producer&& producer) {
    const auto invoke=[](void* context,double& value) {
      return (*static_cast<std::remove_reference_t<Producer>*>(context))(value);
    };
    if(horizontalMaximumEvaluator)
      return horizontalMaximumEvaluator(scalarEvaluationOwner,result,
          const_cast<void*>(static_cast<const void*>(&producer)),invoke);
    if(!evaluation_) return producer(result);
    const auto status=evaluation_->evaluate(horizontalMaximumKey(),sizeof(double),[&] {
      return producer(horizontalMaximum_);
    });
    if(status) result=horizontalMaximum_;
    return status;
  }
  static std::vector<WVVariableEvaluationKey> dependencyKeys(std::size_t count,
      const std::array<std::size_t,16>* gridUses=nullptr,
      const std::array<std::size_t,2>* constantUses=nullptr) {
    std::vector<WVVariableEvaluationKey> keys;
    keys.reserve(2*count+21);
    for(std::size_t index=0;index<count;++index) {
      keys.push_back(prefixKey(index)); keys.push_back(outputKey(index));
    }
    keys.push_back(projectionKey()); keys.push_back(horizontalMaximumKey()); keys.push_back(nonlinearKey());
    for(std::size_t field=0;gridUses && field<4;++field)
      for(std::size_t kind=0;kind<4;++kind) if((*gridUses)[4*field+kind]>1)
        keys.push_back(gridCalculusKey(field,
            static_cast<WVGridCalculusKind>(kind)));
    for(std::size_t direction=0;constantUses && direction<2;++direction)
      if((*constantUses)[direction]>1)
      keys.push_back(constantLaplacianKey(direction));
    return keys;
  }
  // Preparation is allowed only between evaluations. Slots retain capacity,
  // but no value from a completed evaluation remains logically live.
  WVKernelStatus prepareScopedStorage(WVVariableEvaluationPolicy policy,
      const std::vector<WVForcingStage>& stages) {
    if(evaluation_ || initialized_ || owner_)
      return {WVKernelStatusCode::invalidConfiguration,"Cannot prepare active forcing diagnostic storage."};
    try {
      const bool stagesChanged=!storagePrepared_ || preparedStages_!=stages;
      std::vector<WVForcingStage> preparedStages;
      if(stagesChanged) preparedStages=stages;
      // Low memory needs this one prefix scratch allocation. Reserve it before
      // releasing any reuse-policy storage so allocation failure leaves the
      // current policy's prepared workspace intact.
      if(policy==WVVariableEvaluationPolicy::lowMemory)
        staged.reserve(stages.size()*spatial.elementCount());
      flux.resize(coefficientElements_); previous.resize(coefficientElements_);
      temporary.resize(coefficientElements_);
      cumulative.resize(spatial.elementCount()); raw.resize(spatial.elementCount());
      physical.resize(physicalElements_);
      const bool laplacianNeeded=constantLaplacianUseCount[0] || constantLaplacianUseCount[1];
      if(laplacianNeeded) laplacianCoefficients.resize(coefficientElements_);
      else std::vector<WVComplex64>().swap(laplacianCoefficients);
      // Constant hydrostatic diagnostics select u/v/eta from a four-channel
      // coefficient transform even when no Laplacian forcing is present.
      const bool selectTendency=requiresFourChannelTendencySelection &&
          std::any_of(stages.begin(),stages.end(),[](WVForcingStage stage) {
            return stage!=WVForcingStage::spatial;
          });
      if(laplacianNeeded || selectTendency)
        laplacianFields.resize(4*spatial.first*spatial.second*spatial.third);
      else std::vector<double>().swap(laplacianFields);
      if(policy==WVVariableEvaluationPolicy::reuse && nonlinearUseCount>1)
        nonlinearRaw_.resize(spatial.elementCount());
      else std::vector<double>().swap(nonlinearRaw_);
      gridCalculusSlots_.fill(-1);
      constantLaplacianSlots_.fill(-1);
      std::size_t calculusSlots=0,laplacianSlots=0;
      if(policy==WVVariableEvaluationPolicy::reuse) {
        for(std::size_t index=0;index<gridCalculusUseCount.size();++index)
          if(gridCalculusUseCount[index]>1)
            gridCalculusSlots_[index]=static_cast<int>(calculusSlots++);
        for(std::size_t index=0;index<constantLaplacianUseCount.size();++index)
          if(constantLaplacianUseCount[index]>1)
            constantLaplacianSlots_[index]=static_cast<int>(laplacianSlots++);
        const auto R=spatial.first*spatial.second*spatial.third;
        gridCalculusValues_.resize(calculusSlots*R);
        constantLaplacianValues_.resize(laplacianSlots*4*R);
      } else {
        std::vector<double>().swap(gridCalculusValues_);
        std::vector<double>().swap(constantLaplacianValues_);
      }
      if(policy==WVVariableEvaluationPolicy::reuse) {
        builtinNodes.reserve(2);
        prefix.resize(stages.size());
        for(std::size_t index=0;index<stages.size();++index) {
          prefix[index].spatial=stages[index]==WVForcingStage::spatial;
          prefix[index].fields.resize(raw.size());
          if(!prefix[index].spatial) prefix[index].coefficients.resize(flux.size());
          else prefix[index].coefficients.clear();
        }
      } else {
        std::vector<Prefix>().swap(prefix);
      }
      preparedPolicy_=policy;
      if(stagesChanged) preparedStages_.swap(preparedStages);
      storagePrepared_=true;
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,"Unable to prepare forcing prefix storage."};
    } catch(const std::length_error&) {
      return {WVKernelStatusCode::sizeOverflow,"Forcing prefix storage exceeds vector capacity."};
    }
    return WVKernelStatus::ok();
  }
  WVKernelStatus prepareScopedStorageForActive(
      WVVariableEvaluationPolicy policy,
      const std::vector<WVForcingStage>& stages) {
    if(!evaluation_ && !initialized_ && !owner_)
      return prepareScopedStorage(policy,stages);
    if(policy!=WVVariableEvaluationPolicy::reuse || !storagePrepared_ ||
        preparedPolicy_!=policy || preparedStages_!=stages)
      return {WVKernelStatusCode::invalidConfiguration,
          "Active forcing diagnostic storage has an incompatible preparation signature."};
    return WVKernelStatus::ok();
  }
  WVKernelStatus beginScopedEvaluation(WVVariableEvaluationContext& context,
      const std::vector<WVForcingStage>& stages) {
    if(!context.active())
      return {WVKernelStatusCode::invalidConfiguration,"Forcing prefix requires an active evaluation."};
    auto status=prepareScopedStorage(context.policy(),stages);
    if(!status) return status;
    evaluation_=&context; generation_=context.generation();
    calculusAccess_={this,lookupGridCalculus,captureGridCalculus};
    return WVKernelStatus::ok();
  }
  void endScopedEvaluation() noexcept {
    owner_=nullptr; state_={}; evaluation_=nullptr; generation_=0;
    initialized_=false; nextIndex=0; projected=false;
    horizontalMaximum_=0; scalarEvaluationOwner=nullptr; horizontalMaximumEvaluator=nullptr;
    derivativeAccess={}; calculusAccess_={};
    physicalPrepared=false; spatialCaptured=false; captureBuiltinProjection=false;
    flux.clear(); previous.clear(); temporary.clear(); laplacianCoefficients.clear();
    nonlinearRaw_.clear(); gridCalculusValues_.clear();
    constantLaplacianValues_.clear();
    cumulative.clear(); raw.clear(); physical.clear(); laplacianFields.clear(); staged.clear();
    // Keep the outer slots alive: destroying them would free their arenas.
    for(auto& value:prefix) {value.coefficients.clear(); value.fields.clear();}
  }
  WVKernelStatus bind(const void* engine,const WVState& state) {
    if(evaluation_ && (!evaluation_->active() || evaluation_->generation()!=generation_))
      return {WVKernelStatusCode::invalidConfiguration,"Forcing prefix belongs to an expired evaluation."};
    if(!owner_) {owner_=engine; state_=state; return WVKernelStatus::ok();}
    if(!evaluation_)
      return {WVKernelStatusCode::invalidConfiguration,"Reusing forcing diagnostic storage requires an explicit evaluation scope."};
    if(owner_!=engine || state.t!=state_.t || state.t0!=state_.t0)
      return {WVKernelStatusCode::invalidConfiguration,"Forcing prefix belongs to a different owner or state."};
    const WVComplexConstView expected[]={state_.coefficients.Ap,state_.coefficients.Am,state_.coefficients.A0};
    const WVComplexConstView actual[]={state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0};
    for(std::size_t index=0;index<3;++index)
      if(expected[index].data!=actual[index].data ||
          expected[index].shape.rows!=actual[index].shape.rows ||
          expected[index].shape.columns!=actual[index].shape.columns)
        return {WVKernelStatusCode::invalidConfiguration,"Forcing prefix coefficients do not match its immutable state."};
    return WVKernelStatus::ok();
  }
  bool initialized() const noexcept {return initialized_;}
  void markInitialized() noexcept {initialized_=true;}
  bool reusesPrefix() const noexcept {
    return evaluation_ && evaluation_->policy()==WVVariableEvaluationPolicy::reuse;
  }
  WVVariableEvaluationContext* evaluation() const noexcept {return evaluation_;}

  template<class Producer>
  WVKernelStatus evaluateConstantLaplacian(std::size_t direction,
      WVRealFieldBundleConstView& result,Producer&& producer) {
    result={};
    if(direction>=constantLaplacianUseCount.size())
      return {WVKernelStatusCode::invalidConfiguration,
          "Unknown constant-stratification Laplacian direction."};
    const auto R=spatial.first*spatial.second*spatial.third;
    const auto produceScratch=[&]() {
      if(laplacianFields.size()!=4*R)
        return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
            "Constant Laplacian scratch was not prepared."};
      const auto status=producer(WVRealFieldBundleView{laplacianFields.data(),
          {spatial.first,spatial.second,spatial.third,4}});
      if(status) result={laplacianFields.data(),
          {spatial.first,spatial.second,spatial.third,4}};
      return status;
    };
    if(!evaluation_ || constantLaplacianUseCount[direction]<2)
      return produceScratch();
    const auto key=constantLaplacianKey(direction);
    if(!reusesPrefix()) {
      const auto status=evaluation_->evaluate(key,0,produceScratch);
      if(status) (void)evaluation_->evict(key);
      return status;
    }
    const auto slot=constantLaplacianSlots_[direction];
    if(slot<0) return {WVKernelStatusCode::invalidConfiguration,
        "Constant Laplacian cache was not prepared."};
    auto* destination=constantLaplacianValues_.data()+
        static_cast<std::size_t>(slot)*4*R;
    const auto status=evaluation_->evaluate(key,4*R*sizeof(double),[&]() {
      return producer(WVRealFieldBundleView{destination,
          {spatial.first,spatial.second,spatial.third,4}});
    });
    if(status) result={destination,
        {spatial.first,spatial.second,spatial.third,4}};
    return status;
  }

  WVFlux fluxView() { return views(flux); }
  WVFlux temporaryView() { return views(temporary); }
  WVRealFieldBundleView rawView() { return {raw.data(),spatial}; }
  WVRealFieldBundleConstView cumulativeView() const { return {cumulative.data(),spatial}; }
  WVKernelStatus addSpatial(WVRealFieldBundleConstView value) {
    if (value.shape.first!=spatial.first || value.shape.second!=spatial.second ||
        value.shape.third!=spatial.third || value.shape.fourth!=spatial.fourth)
      return {WVKernelStatusCode::invalidShape,"A forcing supplied invalid spatial diagnostic channels."};
    if (!value.data)
      return {WVKernelStatusCode::invalidPointer,"A forcing supplied null spatial diagnostic storage."};
    spatialCaptured=true;
    for (std::size_t i=0;i<raw.size();++i) raw[i]+=value.data[i];
    return WVKernelStatus::ok();
  }
  std::size_t bytes() const noexcept {
    std::size_t bytes=sizeof(*this)+(flux.capacity()+previous.capacity()+temporary.capacity()+
        laplacianCoefficients.capacity())*sizeof(WVComplex64)+
        (cumulative.capacity()+raw.capacity()+physical.capacity()+laplacianFields.capacity()+
         staged.capacity()+nonlinearRaw_.capacity()+gridCalculusValues_.capacity()+
         constantLaplacianValues_.capacity())*sizeof(double);
    bytes+=prefix.capacity()*sizeof(Prefix)+
        builtinNodes.capacity()*sizeof(decltype(builtinNodes)::value_type);
    bytes+=preparedStages_.capacity()*sizeof(WVForcingStage);
    for(const auto& value:prefix)
      bytes+=value.fields.capacity()*sizeof(double)+value.coefficients.capacity()*sizeof(WVComplex64);
    return bytes;
  }

  WVShape2D spectral;
  WVShape4D spatial;
  std::vector<WVComplex64> flux,previous,temporary,laplacianCoefficients;
  std::vector<double> cumulative,raw,physical,laplacianFields,staged;
  std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> builtinNodes;
  bool physicalPrepared=false,spatialCaptured=false;
  bool captureBuiltinProjection=false;
  bool requiresFourChannelTendencySelection=false;
  struct Prefix {
    std::vector<WVComplex64> coefficients;
    std::vector<double> fields;
    bool spatial=false;
  };
  std::vector<Prefix> prefix;
  std::size_t nextIndex=0;
  bool projected=false;

private:
  std::vector<WVForcingStage> preparedStages_;
  WVVariableEvaluationPolicy preparedPolicy_=WVVariableEvaluationPolicy::reuse;
  bool storagePrepared_=false;
  std::vector<double> nonlinearRaw_;
  std::vector<double> gridCalculusValues_,constantLaplacianValues_;
  std::array<int,16> gridCalculusSlots_{};
  std::array<int,2> constantLaplacianSlots_{};
  WVStateCalculusAccess calculusAccess_{};
  double horizontalMaximum_=0;
  std::size_t coefficientElements_=0,physicalElements_=0;
  const void* owner_=nullptr;
  WVState state_;
  WVVariableEvaluationContext* evaluation_=nullptr;
  std::uint64_t generation_=0;
  bool initialized_=false;
  static WVKernelStatus lookupGridCalculus(void* context,std::size_t field,
      WVGridCalculusKind kind,WVRealVolumeConstView& result) {
    auto& workspace=*static_cast<WVForcingDiagnosticWorkspace*>(context);
    result={};
    const auto kindIndex=static_cast<std::size_t>(kind);
    if(field>=4 || kindIndex>=4)
      return {WVKernelStatusCode::invalidConfiguration,
          "Unknown grid-calculus derivative identity."};
    const auto index=4*field+kindIndex;
    if(workspace.gridCalculusUseCount[index]<2 || !workspace.evaluation_ ||
        !workspace.evaluation_->ready(gridCalculusKey(field,kind)))
      return WVKernelStatus::ok();
    const auto slot=workspace.gridCalculusSlots_[index];
    const auto R=workspace.spatial.first*workspace.spatial.second*workspace.spatial.third;
    if(slot<0 || workspace.gridCalculusValues_.size()<
        (static_cast<std::size_t>(slot)+1)*R)
      return {WVKernelStatusCode::invalidConfiguration,
          "Grid-calculus cache has incompatible storage."};
    const auto status=workspace.evaluation_->evaluate(gridCalculusKey(field,kind),
        R*sizeof(double),[](){return WVKernelStatus::ok();});
    if(status) result={workspace.gridCalculusValues_.data()+
        static_cast<std::size_t>(slot)*R,
        {workspace.spatial.first,workspace.spatial.second,workspace.spatial.third}};
    return status;
  }
  static WVKernelStatus captureGridCalculus(void* context,std::size_t field,
      WVGridCalculusKind kind,WVRealVolumeConstView value) {
    auto& workspace=*static_cast<WVForcingDiagnosticWorkspace*>(context);
    const auto kindIndex=static_cast<std::size_t>(kind);
    if(field>=4 || kindIndex>=4)
      return {WVKernelStatusCode::invalidConfiguration,
          "Unknown grid-calculus derivative identity."};
    const auto index=4*field+kindIndex;
    if(workspace.gridCalculusUseCount[index]<2 || !workspace.evaluation_)
      return WVKernelStatus::ok();
    const WVShape3D expected{workspace.spatial.first,workspace.spatial.second,
        workspace.spatial.third};
    if(!value.data) return {WVKernelStatusCode::invalidPointer,
        "Captured grid-calculus storage is null."};
    if(value.shape.first!=expected.first || value.shape.second!=expected.second ||
        value.shape.third!=expected.third)
      return {WVKernelStatusCode::invalidShape,
          "Captured grid-calculus storage has an incompatible shape."};
    const auto key=gridCalculusKey(field,kind);
    if(!workspace.reusesPrefix()) {
      const auto status=workspace.evaluation_->evaluate(key,0,
          [](){return WVKernelStatus::ok();});
      if(status) (void)workspace.evaluation_->evict(key);
      return status;
    }
    const auto slot=workspace.gridCalculusSlots_[index];
    const auto R=expected.elementCount();
    if(slot<0 || workspace.gridCalculusValues_.size()<
        (static_cast<std::size_t>(slot)+1)*R)
      return {WVKernelStatusCode::invalidConfiguration,
          "Grid-calculus cache was not prepared."};
    auto* destination=workspace.gridCalculusValues_.data()+
        static_cast<std::size_t>(slot)*R;
    return workspace.evaluation_->evaluate(key,R*sizeof(double),[&]() {
      std::copy_n(value.data,R,destination); return WVKernelStatus::ok();
    });
  }
  WVFlux views(std::vector<WVComplex64>& data) {
    const auto S=spectral.elementCount();
    if (data.size()==S) return {{},{},{data.data(),spectral}};
    return {{data.data(),spectral},{data.data()+S,spectral},{data.data()+2*S,spectral}};
  }
};

// MATLAB's inverse Fourier transform takes the real part at self-conjugate
// modes. The modal record excludes Nyquist modes, so only Kh=0 needs this
// projection. Averaging the inertial pair preserves its real u/v contribution;
// copying Ap into Am (the model-state constraint) would double an Ap-only delta.
// Only a temporary per-instance difference is changed, never the accumulator.
template<class Geometry>
void projectRealMeanTendency(std::vector<WVComplex64>& difference,const Geometry& g) {
  const auto S=g.Nj*g.Nkl;
  for (std::size_t mode=0;mode<g.Nkl;++mode) if (g.k[mode]==0 && g.l[mode]==0)
    for (std::size_t j=0;j<g.Nj;++j) {
      const auto i=j+g.Nj*mode;
      const auto ap=difference[i],am=difference[S+i];
      difference[i]={.5*ap.real+.5*am.real,.5*ap.imag-.5*am.imag};
      difference[S+i]={difference[i].real,-difference[i].imag};
      difference[2*S+i].imag=0;
    }
}

inline bool forcingArraysOverlap(const void* a,std::size_t n,const void* b,std::size_t m) {
  const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
  return a && b && n && m && (x<=y ? y-x<n : x-y<m);
}

inline WVKernelStatus validatePreparedDiagnosticFields(
    const WVRealFieldBundleConstView* prepared,WVShape4D spatial,std::size_t channels,
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count) {
  if (!prepared) return WVKernelStatus::ok();
  const auto shape=prepared->shape;
  if(shape.first!=spatial.first || shape.second!=spatial.second || shape.third!=spatial.third || shape.fourth!=channels)
    return {WVKernelStatusCode::invalidShape,"Prepared forcing diagnostic fields have the wrong channels."};
  const auto elements=shape.elementCount();
  const auto bytes=elements*sizeof(double),address=reinterpret_cast<std::uintptr_t>(prepared->data);
  if(!address || address%alignof(double) || bytes>UINTPTR_MAX-address)
    return {WVKernelStatusCode::invalidPointer,"Invalid prepared forcing diagnostic field storage."};
  for(const auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0})
    if(forcingArraysOverlap(prepared->data,bytes,input.data,input.shape.elementCount()*sizeof(WVComplex64)))
      return {WVKernelStatusCode::overlappingArrays,"Prepared forcing fields overlap coefficient state."};
  for(std::size_t index=0;index<count;++index)
    if(forcingArraysOverlap(prepared->data,bytes,outputs[index].fields.data,spatial.elementCount()*sizeof(double)))
      return {WVKernelStatusCode::overlappingArrays,"Prepared forcing fields overlap diagnostic output."};
  for(std::size_t index=0;index<elements;++index)
    if(!std::isfinite(prepared->data[index]))
      return {WVKernelStatusCode::invalidConfiguration,"Prepared forcing fields must be finite."};
  return WVKernelStatus::ok();
}

template<class Forcing>
WVKernelStatus validateForcingTendencyOutputs(
    const std::vector<std::unique_ptr<Forcing>>& forcing,WVShape2D spectral,WVShape4D spatial,
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count) {
  if (count && !outputs)
    return {WVKernelStatusCode::invalidPointer,"Missing forcing diagnostic output descriptors."};
  const auto bytes=spatial.elementCount()*sizeof(double);
  for (std::size_t i=0;i<count;++i) {
    const auto& output=outputs[i];
    if (output.executionIndex>=forcing.size())
      return {WVKernelStatusCode::invalidConfiguration,"Forcing diagnostic index is not bound to this schedule."};
    const auto shape=output.fields.shape;
    if (shape.first!=spatial.first || shape.second!=spatial.second ||
        shape.third!=spatial.third || shape.fourth!=spatial.fourth)
      return {WVKernelStatusCode::invalidShape,"Forcing diagnostic output has the wrong spatial channels."};
    const auto address=reinterpret_cast<std::uintptr_t>(output.fields.data);
    if (!address || address%alignof(double) || bytes>UINTPTR_MAX-address)
      return {WVKernelStatusCode::invalidPointer,"Invalid forcing diagnostic output storage."};
    for (const auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0})
      if (forcingArraysOverlap(output.fields.data,bytes,input.data,spectral.elementCount()*sizeof(WVComplex64)))
        return {WVKernelStatusCode::overlappingArrays,"Forcing diagnostic output overlaps model coefficients."};
    for (std::size_t j=0;j<i;++j) {
      if (output.executionIndex==outputs[j].executionIndex)
        return {WVKernelStatusCode::invalidConfiguration,"Repeated forcing diagnostic index must share one output."};
      if (forcingArraysOverlap(output.fields.data,bytes,outputs[j].fields.data,bytes))
        return {WVKernelStatusCode::overlappingArrays,"Forcing diagnostic outputs overlap."};
    }
  }
  if (count) {
    std::size_t last=0;
    for (std::size_t i=0;i<count;++i) last=std::max(last,outputs[i].executionIndex);
    for (std::size_t i=0;i<=last;++i)
      if (!forcing[i]->supportsTendencyDiagnostics())
        return {WVKernelStatusCode::unsupportedOperation,"The resolved forcing has no qualified diagnostic implementation."};
  }
  return WVKernelStatus::ok();
}

// The existing resolved instances remain the only forcing dispatch. Spatial
// operations are accumulated before the stage-boundary projection, matching
// SpatialForcingOperation. Later operations observe the preceding accumulator.
template<class Forcing,class Execute,class Project,class Reconstruct>
WVKernelStatus evaluateForcingTendencySequence(
    const std::vector<std::unique_ptr<Forcing>>& forcing,
    WVForcingDiagnosticWorkspace& work,const WVForcingTendencyOutput* outputs,
    std::size_t count,WVForcingTendencyMetrics& metrics,
    Execute execute,Project project,Reconstruct reconstruct,
    WVFlux* builtinNonlinearProjection=nullptr) {
  if (!count && !builtinNonlinearProjection) return WVKernelStatus::ok();
  if(work.reusesPrefix()) {
    std::size_t last=0;
    for(std::size_t output=0;output<count;++output)
      last=std::max(last,outputs[output].executionIndex);
    if(work.prefix.size()!=forcing.size())
      return {WVKernelStatusCode::invalidConfiguration,"Prepared forcing prefix does not match the resolved schedule."};
    auto& evaluation=*work.evaluation();
    auto flux=work.fluxView();
    struct FluxRollback {
      WVForcingDiagnosticWorkspace& work;
      bool committed=false;
      ~FluxRollback() {
        if(!committed) std::copy(work.previous.begin(),work.previous.end(),work.flux.begin());
      }
    };
    if(builtinNonlinearProjection && work.nextIndex && !work.projected) {
      std::copy(work.flux.begin(),work.flux.end(),work.previous.begin());
      FluxRollback rollback{work};
      const auto status=evaluation.evaluate(work.projectionKey(),0,[&] {
        return project(work.cumulativeView(),flux);
      });
      if(!status) return status;
      rollback.committed=true; work.projected=true;
      ++metrics.spatialProjectionCount;
    }
    while(work.nextIndex<=last) {
      const auto index=work.nextIndex;
      auto& cached=work.prefix[index];
      const bool spatial=forcing[index]->stage()==WVForcingStage::spatial;
      if(cached.spatial!=spatial)
        return {WVKernelStatusCode::invalidConfiguration,"Prepared forcing stage differs from its resolved schedule."};
      if(!spatial && !work.projected) {
        std::copy(work.flux.begin(),work.flux.end(),work.previous.begin());
        FluxRollback rollback{work};
        const auto status=evaluation.evaluate(work.projectionKey(),0,[&] {
          return project(work.cumulativeView(),flux);
        });
        if(!status) return status;
        rollback.committed=true;
        ++metrics.spatialProjectionCount;
        work.projected=true;
      }
      std::copy(work.flux.begin(),work.flux.end(),work.previous.begin());
      FluxRollback rollback{work};
      const auto bytes=spatial ? cached.fields.capacity()*sizeof(double) :
          cached.coefficients.capacity()*sizeof(WVComplex64);
      const bool combined=builtinNonlinearProjection && index==0 && spatial;
      if(combined) {
        work.builtinNodes.clear();
        work.builtinNodes.push_back({work.prefixKey(index),bytes});
        work.builtinNodes.push_back({work.projectionKey(),0});
        work.captureBuiltinProjection=true;
      }
      const auto operation=[&]() -> WVKernelStatus {
        std::fill(work.raw.begin(),work.raw.end(),0.0);
        work.spatialCaptured=false;
        auto produced=execute(*forcing[index],flux);
        if(!produced) return produced;
        if(spatial && !work.spatialCaptured)
          return {WVKernelStatusCode::unsupportedOperation,"The resolved spatial forcing does not expose a diagnostic contribution."};
        ++metrics.forcingEvaluationCount;
        if(spatial) {
          for(std::size_t element=0;element<work.raw.size();++element) {
            const auto before=work.cumulative[element];
            work.cumulative[element]+=work.raw[element];
            cached.fields[element]=work.cumulative[element]-before;
          }
        } else {
          for(std::size_t element=0;element<work.flux.size();++element)
            cached.coefficients[element]={work.flux[element].real-work.previous[element].real,
                                         work.flux[element].imag-work.previous[element].imag};
        }
        return WVKernelStatus::ok();
      };
      const auto status=combined ? evaluation.evaluateGroup(work.builtinNodes,operation) :
          evaluation.evaluate(work.prefixKey(index),bytes,operation);
      work.captureBuiltinProjection=false;
      if(!status) return status;
      rollback.committed=true;
      if(combined) {work.projected=true; ++metrics.spatialProjectionCount;}
      ++work.nextIndex;
    }
    // Materialize only requested physical outputs. Earlier unrequested spectral
    // differences remain compact until a later query actually needs them.
    for(std::size_t output=0;output<count;++output) {
      const auto index=outputs[output].executionIndex;
      auto& cached=work.prefix[index];
      if(cached.spatial) continue;
      const auto status=evaluation.evaluate(work.outputKey(index),
          cached.fields.capacity()*sizeof(double),[&]() {
        // A reconstruction may apply a real-mean projection to its input.
        std::copy(cached.coefficients.begin(),cached.coefficients.end(),work.previous.begin());
        const auto produced=reconstruct(work.previous,{cached.fields.data(),work.spatial});
        if(produced) ++metrics.spectralReconstructionCount;
        return produced;
      });
      if(!status) return status;
    }
    // Publish the requested batch only after every producer succeeds.
    for(std::size_t output=0;output<count;++output) {
      const auto& cached=work.prefix[outputs[output].executionIndex];
      std::copy(cached.fields.begin(),cached.fields.end(),outputs[output].fields.data);
    }
    ++metrics.evaluationCount;
    if(builtinNonlinearProjection) {
      const auto S=work.spectral.elementCount();
      for(std::size_t i=0;i<S;++i) {
        builtinNonlinearProjection->Fp.data[i]=work.flux[i];
        builtinNonlinearProjection->Fm.data[i]=work.flux[S+i];
        builtinNonlinearProjection->F0.data[i]=work.flux[2*S+i];
      }
    }
    return WVKernelStatus::ok();
  }
  // Low-memory calls intentionally replay the needed prefix. Every repeated
  // producer is ledger-visible; no unrequested later forcing is executed.
  if(work.evaluation()) {
    std::fill(work.flux.begin(),work.flux.end(),WVComplex64{});
    std::fill(work.cumulative.begin(),work.cumulative.end(),0.0);
  }
  if (count>work.raw.max_size()/work.raw.size())
    return {WVKernelStatusCode::sizeOverflow,"Forcing diagnostic staging size overflow."};
  work.staged.resize(count*work.raw.size());
  std::size_t last=0;
  for (std::size_t i=0;i<count;++i) last=std::max(last,outputs[i].executionIndex);
  bool projected=false;
  auto flux=work.fluxView();
  for (std::size_t index=0;index<=last;++index) {
    const auto spatial=forcing[index]->stage()==WVForcingStage::spatial;
    if (!spatial && !projected) {
      const auto operation=[&] {return project(work.cumulativeView(),flux);};
      auto status=work.evaluation() ? work.evaluation()->evaluate(work.projectionKey(),0,operation) : operation();
      if (!status) return status;
      if(work.evaluation()) work.evaluation()->evict(work.projectionKey());
      ++metrics.spatialProjectionCount;
      projected=true;
    }
    if (!spatial) std::copy(work.flux.begin(),work.flux.end(),work.previous.begin());
    std::fill(work.raw.begin(),work.raw.end(),0.0);
    work.spatialCaptured=false;
    const auto operation=[&] {return execute(*forcing[index],flux);};
    auto status=work.evaluation() ? work.evaluation()->evaluate(work.prefixKey(index),0,operation) : operation();
    if (!status) return status;
    if(work.evaluation()) work.evaluation()->evict(work.prefixKey(index));
    if (spatial && !work.spatialCaptured)
      return {WVKernelStatusCode::unsupportedOperation,"The resolved spatial forcing does not expose a diagnostic contribution."};
    ++metrics.forcingEvaluationCount;
    WVRealFieldBundleView* destination=nullptr;
    WVRealFieldBundleView destinationView;
    for (std::size_t output=0;output<count;++output)
      if (outputs[output].executionIndex==index) {
        destinationView=outputs[output].fields;
        destinationView.data=work.staged.data()+output*work.raw.size();
        destination=&destinationView;
        break;
      }
    if (spatial) {
      for (std::size_t i=0;i<work.raw.size();++i) {
        const auto previous=work.cumulative[i];
        work.cumulative[i]+=work.raw[i];
        if (destination) destination->data[i]=work.cumulative[i]-previous;
      }
    } else if (destination) {
      for (std::size_t i=0;i<work.flux.size();++i)
        work.previous[i]={work.flux[i].real-work.previous[i].real,
                          work.flux[i].imag-work.previous[i].imag};
      const auto reconstruction=[&] {return reconstruct(work.previous,*destination);};
      status=work.evaluation() ? work.evaluation()->evaluate(work.outputKey(index),0,reconstruction) : reconstruction();
      if (!status) return status;
      if(work.evaluation()) work.evaluation()->evict(work.outputKey(index));
      ++metrics.spectralReconstructionCount;
    }
  }
  for (std::size_t output=0;output<count;++output)
    std::copy_n(work.staged.data()+output*work.raw.size(),work.raw.size(),outputs[output].fields.data);
  ++metrics.evaluationCount;
  return WVKernelStatus::ok();
}

} // namespace wavevortex::runtime::detail
