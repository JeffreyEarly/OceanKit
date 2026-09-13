#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingDiagnosticWorkspace.hpp"
#include "WVScopedStateEvaluation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
namespace wavevortex::runtime {
namespace {
constexpr double pi=3.14159265358979323846;
WVComplex64 add(WVComplex64 a,WVComplex64 b) { return {a.real+b.real,a.imag+b.imag}; }
WVComplex64 multiply(WVComplex64 a,WVComplex64 b) { return {a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real}; }
WVComplex64 scale(WVComplex64 a,double b) { return {a.real*b,a.imag*b}; }
WVComplex64 conjugate(WVComplex64 a) { return {a.real,-a.imag}; }
WVKernelStatus invalid(const char* message) { return {WVKernelStatusCode::invalidConfiguration,message}; }
class BoussinesqErrorPolicy final : public WVIntegrationErrorPolicy {
public:
    BoussinesqErrorPolicy(std::vector<double> wave,std::vector<double> vortex):wave_(std::move(wave)),vortex_(std::move(vortex)) {}
    std::size_t componentCount() const noexcept override { return 3; }
    std::size_t elementCount(std::size_t component) const noexcept override { return component<3 ? wave_.size() : 0; }
    double absoluteTolerance(std::size_t component,std::size_t index) const noexcept override { return component<3 && index<wave_.size() ? (component<2 ? wave_[index] : vortex_[index]) : std::numeric_limits<double>::quiet_NaN(); }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+(wave_.capacity()+vortex_.capacity())*sizeof(double); }
private:
    std::vector<double> wave_,vortex_;
};
}
WVKernelStatus WVBoussinesqForcingEngine::validateSchedule(const WVStratifiedModalGeometry& g,const WVFrozenForcingSchedule& schedule,WVShape2D shape,const WVExtensionCatalog& catalog) {
    if (g.transformClass!="WVTransformBoussinesq" || shape.rows!=g.Nj || shape.columns!=g.Nkl) return invalid("Boussinesq forcing geometry and coefficients differ.");
    WVTransformConstantStratificationConfiguration compatibility;
    compatibility.Nj=g.Nj; compatibility.shouldAntialias=g.shouldAntialias; compatibility.isHydrostatic=false;
    auto s=WVConstantStratificationForcingEngine::validateSchedule(compatibility,schedule,shape,catalog); if (!s) return s;
    for (const auto& entry:schedule.entries) {
        const auto* registration=catalog.forcings().registration(entry.typeIdentifier,entry.contractVersion);
        if (!registration || !registration->boussinesqFactory) return {WVKernelStatusCode::unsupportedOperation,"No paired Boussinesq forcing implementation."};
        const char* type=entry.stage==WVForcingStage::spatial ? "NonhydrostaticSpatial" : entry.stage==WVForcingStage::spectral ? "Spectral" : "SpectralAmplitude";
        if (std::find(registration->forcingTypes.begin(),registration->forcingTypes.end(),type)==registration->forcingTypes.end()) return {WVKernelStatusCode::unsupportedOperation,"Forcing is incompatible with Boussinesq dynamics."};
        for (const char* name:{"ApIndices","AmIndices","A0Indices"}) {
            const auto* value=entry.configuration.value(name); if (!value) continue;
            const auto* indices=std::get_if<std::vector<std::int64_t>>(&value->storage);
            if (!indices) return invalid("Forcing indices must be integers.");
            for (auto index:*indices) if (index<0 || static_cast<std::size_t>(index)>=shape.elementCount()) return invalid("Forcing index outside Boussinesq coefficients.");
        }
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::create(std::shared_ptr<const WVStratifiedModalSource> source,const WVFrozenForcingSchedule& schedule,std::shared_ptr<const WVExtensionCatalog> catalog,std::unique_ptr<WVFFTEngine> fft,std::unique_ptr<WVBoussinesqForcingEngine>& result,const WVVariableKernelServices& services) {
    if (!source || !catalog || !fft) return invalid("Boussinesq forcing requires scientific source, catalog and FFT engine.");
    try {
        const auto& g=source->geometry();
        auto s=validateSchedule(g,schedule,{g.Nj,g.Nkl},*catalog); if (!s) return s;
        auto candidate=std::unique_ptr<WVBoussinesqForcingEngine>(new WVBoussinesqForcingEngine);
        candidate->catalog_=std::move(catalog);
        candidate->evaluationPolicy_=services.variableEvaluationPolicy;
        auto execution=services.execution;
        bool hasNonlinear=false;
        for (const auto& entry:schedule.entries)
            hasNonlinear|=candidate->catalog_->forcings().registration(entry.typeIdentifier,entry.contractVersion)->evaluationDependencies.nonlinearUseCount!=0;
        // Prepare tiled scratch only for a workload that declares its consumer.
        execution.tiledNonlinear=execution.tiledNonlinear && hasNonlinear &&
            services.variableEvaluationPolicy==WVVariableEvaluationPolicy::reuse;
        s=WVTransformBoussinesqKernel::create(std::move(source),std::move(fft),candidate->kernel_,services.matrixBackendFactory,execution); if (!s) return s;
        s=candidate->initialize(schedule); if (!s) return s;
        result=std::move(candidate); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Boussinesq forcing allocation failed."}; }
}
WVKernelStatus WVBoussinesqForcingEngine::initialize(const WVFrozenForcingSchedule& schedule) {
    const auto& g=kernel().geometry();
    preparation_.maximumVerticalMode=static_cast<std::size_t>(*std::max_element(g.j.begin(),g.j.end()));
    for (std::size_t m=0;m<g.Nkl;++m) preparation_.maximumHorizontalComponent=std::max({preparation_.maximumHorizontalComponent,std::abs(g.k[m]),std::abs(g.l[m])});
    std::vector<const WVFrozenForcingEntry*> entries;
    for (const auto& e:schedule.entries) { entries.push_back(&e); preparation_.hasAdaptiveDamping|=catalog_->forcings().registration(e.typeIdentifier,e.contractVersion)->providesAdaptiveDamping; }
    for(const auto& e:schedule.entries) {
        const auto* registration=catalog_->forcings().registration(
            e.typeIdentifier,e.contractVersion);
        if(!registration) return invalid("Forcing evaluation dependencies are unavailable.");
        for(std::size_t i=0;i<gridCalculusUseCount_.size();++i)
            gridCalculusUseCount_[i]+=
                registration->evaluationDependencies.boussinesqGridCalculusUseCount[i];
    }
    for(const auto* entry:entries) {
        const auto& prepare=catalog_->forcings().registration(entry->typeIdentifier,entry->contractVersion)->prepareBoussinesqResolution;
        if(prepare) { const auto status=prepare(*entry,kernel(),preparation_); if(!status) return status; }
    }
    std::stable_sort(entries.begin(),entries.end(),[](const auto* a,const auto* b) { if (a->stage!=b->stage) return a->stage<b->stage; if (a->priority!=b->priority) return a->priority<b->priority; return a->ordinal<b->ordinal; });
    std::ostringstream identifier; identifier<<WVForcingScheduleProfileIdentifier<<':';
    for (const auto* e:entries) {
        std::unique_ptr<WVForcing> forcing;
        auto s=catalog_->forcings().registration(e->typeIdentifier,e->contractVersion)->boussinesqFactory(*e,kernel(),preparation_,forcing); if (!s) return s;
        if (!forcing || forcing->stage()!=e->stage) return invalid("Boussinesq forcing factory returned an incompatible stage.");
        identifier<<e->typeIdentifier<<',';
        if (e->stage==WVForcingStage::spatial) ++metrics_.resolvedSpatialCount;
        else if(e->stage==WVForcingStage::spectral) ++metrics_.resolvedSpectralCount;
        else ++metrics_.resolvedAmplitudeCount;
        metrics_.derivedOperatorBytes+=forcing->persistentBytes(); forcing_.push_back(std::move(forcing));
    }
    scheduleIdentifier_=identifier.str();
    const auto R=kernel().spatialShape().elementCount(),S=kernel().spectralShape().elementCount();
    physical_.resize(4*R); spatial_.resize(4*R); derivative_.resize(2*R); temporary_.resize(3*S);
    for (std::uint32_t channel=0;channel<4;++channel)
        physicalGroup_.push_back({{WVVariableEvaluationNode::physicalField,channel},R*sizeof(double)});
    gridCalculusSlots_.fill(-1);
    std::size_t calculusSlots=0;
    for(std::size_t index=0;index<gridCalculusUseCount_.size();++index)
        if(gridCalculusUseCount_[index]>1)
            gridCalculusSlots_[index]=static_cast<int>(calculusSlots++);
    if(evaluationPolicy_==WVVariableEvaluationPolicy::reuse)
        gridCalculus_.resize(calculusSlots*R);
    std::vector<WVVariableEvaluationKey> evaluationKeys{
        {WVVariableEvaluationNode::physicalField,0}, {WVVariableEvaluationNode::physicalField,1},
        {WVVariableEvaluationNode::physicalField,2}, {WVVariableEvaluationNode::physicalField,3},
        {WVVariableEvaluationNode::reduction,0}, {WVVariableEvaluationNode::reduction,1},
        {WVVariableEvaluationNode::forcingTendency,0}};
    for(std::size_t index=0;index<gridCalculusUseCount_.size();++index)
        if(gridCalculusUseCount_[index]>1)
            evaluationKeys.push_back({WVVariableEvaluationNode::gridCalculus,
                static_cast<std::uint32_t>(index/4),0,
                static_cast<std::uint32_t>(index%4)});
    auto prepared=evaluation_.prepare(evaluationKeys);
    if(!prepared) return prepared;
    metrics_.scheduleBytes=scheduleIdentifier_.capacity()+forcing_.capacity()*sizeof(std::unique_ptr<WVForcing>);
    metrics_.workspaceCapacityBytes=(physical_.capacity()+spatial_.capacity()+derivative_.capacity()+gridCalculus_.capacity())*sizeof(double)+(temporary_.capacity()+nonlinearCache_.capacity())*sizeof(WVComplex64);
    metrics_.workspaceHighWaterBytes=metrics_.workspaceCapacityBytes;
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::validateVariableEvaluationPolicyChange(
    WVVariableEvaluationPolicy policy) const noexcept {
    if(evaluation_.active() || executing_) return {WVKernelStatusCode::reentrantExecution,"Cannot change an active evaluation policy."};
    if(policy!=WVVariableEvaluationPolicy::reuse && policy!=WVVariableEvaluationPolicy::lowMemory)
        return {WVKernelStatusCode::invalidConfiguration,"Unknown variable evaluation policy."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) {
    auto status=validateVariableEvaluationPolicyChange(policy); if(!status) return status;
    const auto wholeFluxCount=std::count_if(forcing_.begin(),forcing_.end(),[](const auto& forcing){return forcing->producesCompleteFlux();});
    std::vector<WVComplex64> preparedNonlinear;
    std::vector<double> preparedCalculus;
    try {
        if(policy==WVVariableEvaluationPolicy::reuse && wholeFluxCount>1 && nonlinearCache_.empty())
            preparedNonlinear.resize(3*kernel().spectralShape().elementCount());
        if(policy==WVVariableEvaluationPolicy::reuse && gridCalculus_.empty()) {
            const auto slots=std::count_if(gridCalculusUseCount_.begin(),gridCalculusUseCount_.end(),
                [](std::size_t count){return count>1;});
            preparedCalculus.resize(slots*kernel().spatialShape().elementCount());
        }
    } catch(const std::bad_alloc&) {return {WVKernelStatusCode::allocationFailure,"Unable to allocate shared nonlinear result."};}
    if(policy==WVVariableEvaluationPolicy::reuse) {
        if(!preparedNonlinear.empty()) nonlinearCache_.swap(preparedNonlinear);
        if(!preparedCalculus.empty()) gridCalculus_.swap(preparedCalculus);
    } else {
        std::vector<WVComplex64>{}.swap(nonlinearCache_);
        std::vector<double>{}.swap(gridCalculus_);
    }
    evaluationPolicy_=policy;
    metrics_.workspaceCapacityBytes=(physical_.capacity()+spatial_.capacity()+derivative_.capacity()+gridCalculus_.capacity())*sizeof(double)+(temporary_.capacity()+nonlinearCache_.capacity())*sizeof(WVComplex64);
    metrics_.workspaceHighWaterBytes=std::max(metrics_.workspaceHighWaterBytes,metrics_.workspaceCapacityBytes);
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::beginStateEvaluation(const WVState& state) {
    if (evaluation_.active() || executing_) return {WVKernelStatusCode::reentrantExecution,"State evaluation is already active."};
    auto status=kernel().beginStateEvaluation(state); if(!status) return status;
    status=evaluation_.begin(this,evaluationPolicy_);
    if(!status) {kernel().endStateEvaluation(); return status;}
    evaluationState_=state;
    return WVKernelStatus::ok();
}
void WVBoussinesqForcingEngine::endStateEvaluation() noexcept {
    evaluation_.end(); kernel().endStateEvaluation(); evaluationState_={};
}
WVKernelStatus WVBoussinesqForcingEngine::validateStateEvaluation(const WVState& state) const {
    if(!evaluation_.active() || state.t!=evaluationState_.t || state.t0!=evaluationState_.t0)
        return invalid("State does not belong to the active evaluation.");
    const WVComplexConstView a[]={state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0};
    const WVComplexConstView b[]={evaluationState_.coefficients.Ap,evaluationState_.coefficients.Am,evaluationState_.coefficients.A0};
    for(std::size_t i=0;i<3;++i) if(a[i].data!=b[i].data || a[i].shape.rows!=b[i].shape.rows || a[i].shape.columns!=b[i].shape.columns)
        return invalid("Coefficients do not belong to the active evaluation.");
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::ensurePhysicalField(const WVState& state,std::size_t channel) {
    const auto shape=kernel().spatialShape(); const auto R=shape.elementCount();
    const WVBoussinesqField names[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    return evaluation_.evaluate({WVVariableEvaluationNode::physicalField,static_cast<std::uint32_t>(channel)},R*sizeof(double),[&] {
        return kernel().transformStateField(state,names[channel],{physical_.data()+channel*R,shape});
    });
}
WVKernelStatus WVBoussinesqForcingEngine::nonlinearFlux(const WVState& state,WVFlux& flux) {
    if (executing_) return {WVKernelStatusCode::reentrantExecution,"Boussinesq forcing workspace is active."};
    detail::WVScopedStateEvaluation<WVBoussinesqForcingEngine> evaluation(*this,state);
    if(!evaluation.status()) return evaluation.status();
    const auto bytes=kernel().spectralShape().elementCount()*sizeof(WVComplex64);
    for (auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0})
        for (auto output:{flux.Fp,flux.Fm,flux.F0}) {
            const auto a=reinterpret_cast<std::uintptr_t>(input.data),b=reinterpret_cast<std::uintptr_t>(output.data);
            if (a && b && (a<=b ? b-a<bytes : a-b<bytes)) return {WVKernelStatusCode::overlappingArrays,"Boussinesq RHS output overlaps its input state."};
        }
    // Validate output storage before clearing it; state was validated at scope entry.
    auto s=kernel().validateFluxOutput(state,flux); if (!s) return s;
    executing_=true; struct Guard { bool& x; ~Guard(){x=false;} } guard{executing_};
    for (auto x:{flux.Fp,flux.Fm,flux.F0}) std::fill_n(x.data,x.shape.elementCount(),WVComplex64{});
    WVForcingExecutionContext context; bool initialized=true;
    context.boussinesq_=this; context.state_=&state; context.flux_=&flux; context.outputInitialized_=&initialized;
    if (!linearDynamics_) for (const auto& forcing:forcing_) { s=forcing->addRightHandSide(context); if (!s) return s; }
    ++metrics_.evaluationCount; return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::physicalFields(const WVState& state,WVRealFieldBundleConstView& fields) {
    if (diagnosticWorkspace_) {
        auto& work=*diagnosticWorkspace_;
        const auto shape=kernel().spatialShape(); const auto R=shape.elementCount();
        if (!work.physicalPrepared) {
            const WVBoussinesqField names[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
            for (std::size_t channel=0;channel<4;++channel) {
                auto status=kernel().transformStateField(state,names[channel],{work.physical.data()+channel*R,shape});
                if (!status) return status;
            }
            work.physicalPrepared=true; ++metrics_.physicalFieldReconstructionCount;
        } else ++metrics_.physicalFieldReuseCount;
        fields={work.physical.data(),{shape.first,shape.second,shape.third,4}};
        return WVKernelStatus::ok();
    }
    detail::WVScopedStateEvaluation<WVBoussinesqForcingEngine> evaluation(*this,state);
    if(!evaluation.status()) return evaluation.status();
    const auto shape=kernel().spatialShape();
    bool reconstructed=false;
    for(std::size_t channel=0;channel<4;++channel) {
        reconstructed|=!evaluation_.ready({WVVariableEvaluationNode::physicalField,static_cast<std::uint32_t>(channel)});
        auto status=ensurePhysicalField(state,channel); if(!status) return status;
    }
    if(reconstructed) ++metrics_.physicalFieldReconstructionCount;
    else ++metrics_.physicalFieldReuseCount;
    fields={physical_.data(),{shape.first,shape.second,shape.third,4}};
    return WVKernelStatus::ok();
}

WVRealFieldBundleView WVBoussinesqForcingEngine::clearedSpatialTendency() {
    std::fill(spatial_.begin(),spatial_.end(),0); const auto g=kernel().spatialShape(); return {spatial_.data(),{g.first,g.second,g.third,4}};
}
WVKernelStatus WVBoussinesqForcingEngine::addProjectedSpatialTendency(const WVState& state,WVRealFieldBundleConstView fields,WVFlux& flux) {
    if (diagnosticWorkspace_) {
        return diagnosticWorkspace_->addSpatial(fields);
    }
    const auto R=kernel().spatialShape().elementCount(),S=kernel().spectralShape().elementCount(); const auto shape=kernel().spectralShape();
    if (fields.shape.first!=kernel().geometry().Nx || fields.shape.second!=kernel().geometry().Ny || fields.shape.third!=kernel().geometry().Nz || fields.shape.fourth!=4 || !fields.data) return invalid("Boussinesq tendency requires [Nx,Ny,Nz,4].");
    WVMutableCoefficients out{{temporary_.data(),shape},{temporary_.data()+S,shape},{temporary_.data()+2*S,shape}};
    auto s=kernel().transformUVWEtaToWaveVortex({fields.data,kernel().spatialShape()},{fields.data+R,kernel().spatialShape()},{fields.data+2*R,kernel().spatialShape()},{fields.data+3*R,kernel().spatialShape()},state.t,state.t0,out); if (!s) return s;
    const WVComplexView target[]={flux.Fp,flux.Fm,flux.F0};
    for(int j=0;j<3;++j) for(std::size_t i=0;i<S;++i) target[j].data[i]=add(target[j].data[i],temporary_[j*S+i]);
    ++metrics_.spatialTendencyProjectionCount; return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::addNonlinearFlux(const WVState& state,WVFlux& flux) {
    if (diagnosticWorkspace_) {
        if(diagnosticWorkspace_->captureBuiltinProjection &&
            !diagnosticWorkspace_->physicalPrepared && kernel().supportsTiledNonlinear()) {
            const auto shape=kernel().spatialShape();
            WVRealFieldBundleView fields{diagnosticWorkspace_->physical.data(),
                {shape.first,shape.second,shape.third,4}};
            auto raw=diagnosticWorkspace_->rawView();
            const auto status=kernel().nonlinearFluxAndFields(state,flux,fields,&raw);
            if(status) {
                diagnosticWorkspace_->physicalPrepared=true;
                diagnosticWorkspace_->spatialCaptured=true;
                ++metrics_.nonlinearProducerCount;
                ++metrics_.physicalFieldReconstructionCount;
            }
            return status;
        }
        WVRealFieldBundleConstView fields;
        auto status=physicalFields(state,fields); if (!status) return status;
        diagnosticWorkspace_->spatialCaptured=true;
        auto raw=diagnosticWorkspace_->rawView();
        auto temporary=diagnosticWorkspace_->captureBuiltinProjection ? flux :
            diagnosticWorkspace_->temporaryView();
        return diagnosticWorkspace_->evaluateNonlinearRaw([&] {
            ++metrics_.nonlinearProducerCount;
            return kernel().nonlinearFlux(state,temporary,&raw,&fields,
                diagnosticWorkspace_->captureBuiltinProjection,
                diagnosticWorkspace_->stateDerivativeAccess());
        });
    }
    const auto S=kernel().spectralShape().elementCount(); const auto shape=kernel().spectralShape();
    auto& storage=nonlinearCache_.empty() ? temporary_ : nonlinearCache_;
    WVFlux tmp{{storage.data(),shape},{storage.data()+S,shape},{storage.data()+2*S,shape}};
    const auto produce=[&] {
        bool allEmpty=true;
        for (const auto& node:physicalGroup_) allEmpty&=!evaluation_.ready(node.first);
        if (evaluationPolicy_==WVVariableEvaluationPolicy::reuse && allEmpty && kernel().supportsTiledNonlinear()) {
            const auto grid=kernel().spatialShape();
            const auto result=evaluation_.evaluateGroup(physicalGroup_,[&] {
                ++metrics_.nonlinearProducerCount;
                return kernel().nonlinearFluxAndFields(state,tmp,
                    {physical_.data(),{grid.first,grid.second,grid.third,4}});
            });
            if (result) ++metrics_.physicalFieldReconstructionCount;
            return result;
        }
        WVRealFieldBundleConstView fields;
        auto status=physicalFields(state,fields); if(!status) return status;
        ++metrics_.nonlinearProducerCount;
        return kernel().nonlinearFlux(state,tmp,nullptr,&fields);
    };
    const auto status=(nonlinearCache_.empty() && evaluationPolicy_==WVVariableEvaluationPolicy::reuse) ? produce() : evaluation_.evaluate({WVVariableEvaluationNode::forcingTendency,0},3*S*sizeof(WVComplex64),produce);
    if(!status) return status;
    const WVComplexView target[]={flux.Fp,flux.Fm,flux.F0};
    for(int j=0;j<3;++j) for(std::size_t i=0;i<S;++i) target[j].data[i]=add(target[j].data[i],storage[j*S+i]);
    evaluation_.evict({WVVariableEvaluationNode::forcingTendency,0});
    return WVKernelStatus::ok();
}

WVKernelStatus WVBoussinesqForcingEngine::gridSecondDerivative(
    WVRealVolumeConstView input,std::size_t field,std::size_t kind,
    WVBoussinesqFamily family,const double*& result) {
    result=nullptr;
    if(field>=4 || kind>=4) return invalid("Unknown Boussinesq grid-calculus identity.");
    const auto shape=kernel().spatialShape(); const auto R=shape.elementCount();
    auto* access=diagnosticWorkspace_ ? diagnosticWorkspace_->stateCalculusAccess() : nullptr;
    const auto calculusKind=static_cast<detail::WVGridCalculusKind>(kind);
    WVRealVolumeConstView cached{};
    if(access && access->lookup) {
        auto status=access->lookup(access->context,field,calculusKind,cached);
        if(!status) return status;
    }
    if(cached.data) {
        if(cached.shape.first!=shape.first || cached.shape.second!=shape.second ||
            cached.shape.third!=shape.third)
            return invalid("Cached Boussinesq grid derivative has an incompatible shape.");
        result=cached.data; return WVKernelStatus::ok();
    }
    const auto produce=[&](double* destination) {
        WVKernelStatus status;
        if(kind<2) {
            const bool x=kind==0;
            status=kernel().differentiateHorizontal(input,x,{derivative_.data(),shape});
            if(status) status=kernel().differentiateHorizontal(
                {derivative_.data(),shape},x,{destination,shape});
        } else status=kernel().differentiateVertical(input,family,2,{destination,shape});
        if(status) ++metrics_.gridCalculusProducerCount[field][kind];
        return status;
    };
    double* scratch=derivative_.data()+(kind<2 ? R : 0);
    const auto index=4*field+kind;
    if(access) {
        auto status=produce(scratch); if(!status) return status;
        if(access->capture) {
            status=access->capture(access->context,field,calculusKind,{scratch,shape});
            if(!status) return status;
        }
        result=scratch; return WVKernelStatus::ok();
    }
    if(gridCalculusUseCount_[index]>1 && evaluation_.active()) {
        const auto key=detail::WVForcingDiagnosticWorkspace::gridCalculusKey(field,calculusKind);
        if(evaluationPolicy_==WVVariableEvaluationPolicy::reuse) {
            const auto slot=gridCalculusSlots_[index];
            if(slot<0 || gridCalculus_.size()<(static_cast<std::size_t>(slot)+1)*R)
                return invalid("Boussinesq grid-calculus cache was not prepared.");
            auto* destination=gridCalculus_.data()+static_cast<std::size_t>(slot)*R;
            const auto status=evaluation_.evaluate(key,R*sizeof(double),[&] {
                return produce(destination);
            });
            if(status) result=destination;
            return status;
        }
        const auto status=evaluation_.evaluate(key,0,[&] {return produce(scratch);});
        if(status) evaluation_.evict(key);
        if(status) result=scratch;
        return status;
    }
    const auto status=produce(scratch); if(status) result=scratch; return status;
}

WVKernelStatus WVBoussinesqForcingEngine::addLaplacianDamping(const WVState& state,double nu,double kappa,WVLaplacianDirection direction,WVFlux& flux) {
    WVRealFieldBundleConstView fields; auto s=physicalFields(state,fields); if (!s) return s;
    const auto shape=kernel().spatialShape(); const auto R=shape.elementCount(); auto tendency=clearedSpatialTendency();
    for (std::size_t j=0;j<4;++j) {
        const double coefficient=j==3 ? kappa : nu; const WVRealVolumeConstView a{fields.data+j*R,shape};
        if (direction==WVLaplacianDirection::vertical) {
            const double* derivative=nullptr;
            s=gridSecondDerivative(a,j,j>=2 ? 3 : 2,
                j>=2 ? WVBoussinesqFamily::G : WVBoussinesqFamily::F,derivative); if(!s) return s;
            for(std::size_t i=0;i<R;++i) tendency.data[j*R+i]=coefficient*(derivative[i]-(j==3 ? kernel().geometry().dLnN2[i/(shape.first*shape.second)] : 0));
        } else {
            for(std::size_t kind=0;kind<2;++kind) {
                const double* derivative=nullptr;
                s=gridSecondDerivative(a,j,kind,WVBoussinesqFamily::F,derivative); if(!s) return s;
                for(std::size_t i=0;i<R;++i) tendency.data[j*R+i]+=coefficient*derivative[i];
            }
        }
    }
    return addProjectedSpatialTendency(state,{tendency.data,tendency.shape},flux);
}
WVKernelStatus WVBoussinesqForcingEngine::addVerticalDiffusivity(const WVState& state,double kappa,bool forceMean,WVFlux& flux) {
    WVRealFieldBundleConstView fields; auto s=physicalFields(state,fields); if (!s) return s;
    const auto shape=kernel().spatialShape(); const auto R=shape.elementCount(); auto tendency=clearedSpatialTendency();
    const double* derivative=nullptr;
    s=gridSecondDerivative({fields.data+3*R,shape},3,3,WVBoussinesqFamily::G,derivative); if(!s) return s;
    for (std::size_t i=0;i<R;++i) tendency.data[3*R+i]=kappa*(derivative[i]-(forceMean ? kernel().geometry().dLnN2[i/(shape.first*shape.second)] : 0));
    return addProjectedSpatialTendency(state,{tendency.data,tendency.shape},flux);
}
WVKernelStatus WVBoussinesqForcingEngine::horizontalSpeedMaximum(const WVState& state,double& uv) {
    if(diagnosticWorkspace_) {
        return diagnosticWorkspace_->evaluateHorizontalMaximum(uv,[&](double& maximum) {
            WVRealFieldBundleConstView fields; auto status=physicalFields(state,fields); if(!status) return status;
            const auto R=kernel().spatialShape().elementCount(); maximum=0;
            for(std::size_t i=0;i<R;++i) maximum=std::max(maximum,std::hypot(fields.data[i],fields.data[R+i]));
            ++metrics_.horizontalSpeedReductionCount;
            return WVKernelStatus::ok();
        });
    }
    detail::WVScopedStateEvaluation<WVBoussinesqForcingEngine> scope(*this,state);
    if(!scope.status()) return scope.status();
    auto status=evaluation_.evaluate({WVVariableEvaluationNode::reduction,0},sizeof(double),[&] {
        for(std::size_t c=0;c<2;++c) {auto s=ensurePhysicalField(state,c); if(!s) return s;}
        const auto R=kernel().spatialShape().elementCount(); horizontalMaximum_=0;
        for(std::size_t i=0;i<R;++i) horizontalMaximum_=std::max(horizontalMaximum_,std::hypot(physical_[i],physical_[R+i]));
        ++metrics_.horizontalSpeedReductionCount;
        return WVKernelStatus::ok();
    });
    if(status) uv=horizontalMaximum_;
    return status;
}
WVKernelStatus WVBoussinesqForcingEngine::speedMaxima(const WVState& state,double& uv,double& w) {
    detail::WVScopedStateEvaluation<WVBoussinesqForcingEngine> scope(*this,state);
    if(!scope.status()) return scope.status();
    auto status=horizontalSpeedMaximum(state,uv); if(!status) return status;
    status=evaluation_.evaluate({WVVariableEvaluationNode::reduction,1},sizeof(double),[&] {
        auto s=ensurePhysicalField(state,2); if(!s) return s;
        const auto R=kernel().spatialShape().elementCount(); verticalMaximum_=0;
        for(std::size_t i=0;i<R;++i) verticalMaximum_=std::max(verticalMaximum_,std::abs(physical_[2*R+i]));
        ++metrics_.verticalSpeedReductionCount;
        return WVKernelStatus::ok();
    });
    if(status) w=verticalMaximum_;
    return status;
}
WVKernelStatus WVBoussinesqForcingEngine::addAdaptiveDamping(const WVState& state,const std::vector<double>& damping,WVFlux& flux) {
    double uv; auto s=horizontalSpeedMaximum(state,uv); if (!s) return s;
    const WVComplexConstView a[]={state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0}; const WVComplexView b[]={flux.Fp,flux.Fm,flux.F0};
    for (int j=0;j<3;++j) for(std::size_t i=0;i<damping.size();++i) b[j].data[i]=add(b[j].data[i],scale(a[j].data[i],uv*damping[i]));
    return WVKernelStatus::ok();
}
WVKernelStatus WVBoussinesqForcingEngine::addPseudoTopographicGeneration(const WVState& state,const WVPseudoTopographicOperators& o,WVFlux& flux) {
    const auto& c=o.configuration; const double elapsed=state.t-c.startTime;
    if(elapsed<0) return WVKernelStatus::ok();
    const double ramp=c.rampDuration==0 || elapsed>=c.rampDuration ? 1 : .5*(1-std::cos(pi*elapsed/c.rampDuration));
    const WVComplex64 oscillation{std::cos(c.frequency*elapsed),-std::sin(c.frequency*elapsed)};
    const double u=ramp*multiply(c.barotropicVelocityAmplitude[0],oscillation).real,v=ramp*multiply(c.barotropicVelocityAmplitude[1],oscillation).real;
    WVComplexConstView phases;
    const auto phaseStatus=kernel().preparedPhase(state,phases); if(!phaseStatus) return phaseStatus;
    for(std::size_t i=0;i<o.responsePlusX.size();++i) {
        const auto phase=phases.data[i];
        flux.Fp.data[i]=add(flux.Fp.data[i],multiply(add(scale(o.responsePlusX[i],u),scale(o.responsePlusY[i],v)),conjugate(phase)));
        flux.Fm.data[i]=add(flux.Fm.data[i],multiply(add(scale(o.responseMinusX[i],u),scale(o.responseMinusY[i],v)),phase));
    }
    return WVKernelStatus::ok();
}
WVStateConstraintResult WVBoussinesqForcingEngine::restoreForcingAmplitudes(WVMutableCoefficients& a) {
    if(evaluation_.active()) return {{WVKernelStatusCode::invalidConfiguration,"Cannot constrain coefficients during an immutable evaluation."},0,false};
    // Match MATLAB's forcing-owned accepted-state restoration. The transform's
    // modal validity is established by its input/projection contracts.
    for(auto view:{a.Ap,a.Am,a.A0}) if (!view.data || view.shape.rows!=kernel().geometry().Nj || view.shape.columns!=kernel().geometry().Nkl) return {invalid("Invalid Boussinesq accepted-state shape."),0,false};
    std::size_t modified=0; bool fsal=true;
    for (const auto& forcing:forcing_) { auto result=forcing->applyConstraint(a); if (!result) return result; modified+=result.modifiedCoefficientCount; fsal&=result.fsalCompatible; metrics_.restoredCoefficientCount+=forcing->constraintWriteCount(); }
    return {WVKernelStatus::ok(),modified,fsal};
}
WVKernelStatus WVBoussinesqForcingEngine::createErrorPolicy(double tolerance,std::unique_ptr<WVIntegrationErrorPolicy>& result) const {
    if (!(tolerance>0) || !std::isfinite(tolerance)) return invalid("Adaptive tolerance must be finite and positive.");
    try {
        const auto& g=kernel().geometry(); std::vector<double> kh(g.Nkl),unique;
        for(std::size_t i=0;i<g.Nkl;++i) kh[i]=std::hypot(g.k[i],g.l[i]);
        unique=kh; std::sort(unique.begin(),unique.end()); unique.erase(std::unique(unique.begin(),unique.end()),unique.end());
        if(unique.size()<2) return invalid("Adaptive tolerances require resolved radial modes.");
        double dk=0; for(std::size_t i=1;i<unique.size();++i) dk=std::max(dk,unique[i]-unique[i-1]);
        std::vector<double> centers; for(double k=0;k<=unique.back()+.5*dk;k+=dk) centers.push_back(k);
        std::vector<std::size_t> bins(g.Nkl,centers.size()),counts(centers.size());
        for(std::size_t m=0;m<g.Nkl;++m) for(std::size_t b=0;b<centers.size();++b) if(centers[b]-.5*dk<kh[m] && kh[m]<=centers[b]+.5*dk) { bins[m]=b; ++counts[b]; break; }
        std::vector<double> wave(g.Nj*g.Nkl,1),vortex(wave.size(),1);
        for(std::size_t i=0;i<wave.size();++i) {
            const auto b=bins[i/g.Nj]; if (b==centers.size()) return invalid("Unresolved adaptive radial bin.");
            const double energy=(centers[b]+.5*dk-std::max(centers[b]-.5*dk,0.0))/counts[b]; const auto& f=kernel().factors()[i];
            if (f.waveEnergy>0) wave[i]=tolerance*std::sqrt(energy/f.waveEnergy);
            if (f.balancedEnergy>0) vortex[i]=tolerance*std::sqrt(energy/f.balancedEnergy);
        }
        result=std::make_unique<BoussinesqErrorPolicy>(std::move(wave),std::move(vortex)); return WVKernelStatus::ok();
    } catch(const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Boussinesq tolerance allocation failed."}; }
}
std::size_t WVBoussinesqForcingEngine::persistentBytes() const noexcept {
    return sizeof(*this)+physicalGroup_.capacity()*sizeof(decltype(physicalGroup_)::value_type)+evaluation_.persistentBytes()+(kernel_ ? kernel_->persistentBytes() : 0)+metrics_.scheduleBytes+metrics_.derivedOperatorBytes+metrics_.workspaceCapacityBytes;
}

const WVForcingEvaluationDependencies*
WVBoussinesqForcingEngine::forcingEvaluationDependencies(
    std::size_t index) const noexcept {
    const auto* forcing=forcingInstance(index);
    if(!forcing || !catalog_) return nullptr;
    const auto* registration=catalog_->forcings().registration(
        forcing->typeIdentifier(),forcing->contractVersion());
    return registration ? &registration->evaluationDependencies : nullptr;
}
WVKernelStatus WVBoussinesqForcingEngine::evaluateForcingTendencies(
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,
    const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session) {
    return evaluateForcingTendenciesImpl(state,outputs,count,preparedPhysical,
        session,nullptr);
}
WVKernelStatus WVBoussinesqForcingEngine::evaluateForcingTendenciesImpl(
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,
    const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session,
    WVFlux* builtinNonlinearProjection) {
    if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing diagnostics require an idle engine."};
    tendencyMetrics_.workspaceLastPeakBytes=0;
    const auto shape=kernel().spatialShape();
    const auto R=shape.elementCount();
    const WVShape4D spatial{shape.first,shape.second,shape.third,4};
    auto status=detail::validateForcingTendencyOutputs(forcing_,kernel().spectralShape(),spatial,state,outputs,count);
    if (!status || (!count && !builtinNonlinearProjection)) return status;
    status=detail::validatePreparedDiagnosticFields(preparedPhysical,spatial,4,state,outputs,count);
    if (!status) return status;
    detail::WVScopedStateEvaluation<WVTransformBoussinesqKernel> kernelScope(kernel(),state);
    if(!kernelScope.status()) return kernelScope.status();
    try {
        detail::WVForcingDiagnosticLedger localLedger(tendencyMetrics_);
        std::unique_ptr<detail::WVForcingDiagnosticWorkspace> local;
        if (!session) {
            local=std::make_unique<detail::WVForcingDiagnosticWorkspace>(kernel().spectralShape(),spatial);
            local->gridCalculusUseCount=gridCalculusUseCount_;
            std::vector<WVForcingStage> stages;
            local->nonlinearUseCount=0;
            for(std::size_t index=0;index<forcing_.size();++index) {
                const auto& forcing=forcing_[index];
                stages.push_back(forcing->stage());
                const auto* dependencies=forcingEvaluationDependencies(index);
                if(!dependencies) return invalid("Forcing evaluation dependencies are unavailable.");
                local->nonlinearUseCount+=dependencies->nonlinearUseCount;
            }
            status=localLedger.context.prepare(detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                forcing_.size(),&gridCalculusUseCount_));
            if(!status) return status;
            status=localLedger.context.begin(this,evaluationPolicy_); if(!status) return status;
            status=local->beginScopedEvaluation(localLedger.context,stages); if(!status) return status;
            session=local.get();
        }
        auto& work=*session;
        status=work.bind(this,state); if (!status) return status;
        const auto spectral=kernel().spectralShape(); const auto S=spectral.elementCount();
        if (work.spectral.rows!=spectral.rows || work.spectral.columns!=spectral.columns ||
            work.spatial.first!=spatial.first || work.spatial.second!=spatial.second ||
            work.spatial.third!=spatial.third || work.spatial.fourth!=spatial.fourth ||
            work.flux.size()!=3*S || work.previous.size()!=3*S || work.temporary.size()!=3*S ||
            work.cumulative.size()!=spatial.elementCount() || work.raw.size()!=spatial.elementCount() ||
            work.physical.size()!=4*R)
            return {WVKernelStatusCode::invalidShape,"Forcing diagnostic session has incompatible Boussinesq storage."};
        auto flux=work.fluxView();
        if (!work.initialized()) {
            if (preparedPhysical) {
                std::copy_n(preparedPhysical->data,preparedPhysical->shape.elementCount(),work.physical.data());
                work.physicalPrepared=true;
            }
            status=kernel().evolveCoefficients(state,{flux.Fp,flux.Fm,flux.F0}); if (!status) return status;
            std::fill(work.flux.begin(),work.flux.end(),WVComplex64{});
            work.markInitialized();
        }
        executing_=true; diagnosticWorkspace_=&work;
        struct Guard {
            WVBoussinesqForcingEngine& engine;
            ~Guard() {
                engine.tendencyMetrics_.workspaceLastPeakBytes=engine.diagnosticWorkspace_->bytes();
                engine.tendencyMetrics_.workspaceHighWaterBytes=std::max(engine.tendencyMetrics_.workspaceHighWaterBytes,engine.diagnosticWorkspace_->bytes());
                engine.tendencyMetrics_.workspaceLiveBytes=0;
                engine.diagnosticWorkspace_=nullptr; engine.executing_=false;
            }
        } guard{*this};
        tendencyMetrics_.workspaceLiveBytes=work.bytes();
        bool initialized=true;
        WVForcingExecutionContext context; context.boussinesq_=this;
        context.state_=&state; context.outputInitialized_=&initialized;
        return detail::evaluateForcingTendencySequence(forcing_,work,outputs,count,tendencyMetrics_,
            [&](const WVForcing& forcing,WVFlux& destination) {
                context.flux_=&destination;
                return forcing.addRightHandSide(context);
            },
            [&](WVRealFieldBundleConstView fields,WVFlux& destination) {
                return kernel().transformUVWEtaToWaveVortex({fields.data+0*R,shape},{fields.data+1*R,shape},{fields.data+2*R,shape},{fields.data+3*R,shape},state.t,state.t0,
                    {destination.Fp,destination.Fm,destination.F0});
            },
            [&](std::vector<WVComplex64>& difference,WVRealFieldBundleView destination) {
                detail::projectRealMeanTendency(difference,kernel().geometry());
                const auto spectral=kernel().spectralShape(); const auto S=spectral.elementCount();
                const WVState delta{state.t,state.t0,{{difference.data(),spectral},
                    {difference.data()+S,spectral},{difference.data()+2*S,spectral}}};
                return kernel().transformCoefficientTendencyToUVWEta(delta,destination);
            },builtinNonlinearProjection);
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Unable to allocate event forcing diagnostic workspace."};
    }
}

}
