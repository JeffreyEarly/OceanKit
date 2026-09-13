#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WVReferenceFFTEngine.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#include "WVStratifiedModalTestFixture.hpp"
#include <iostream>
#include <limits>
using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::test_fixture;
static_assert(!std::is_copy_assignable<WVStratifiedModalRecord>::value,"Published scientific records must not be reassigned.");
static_assert(!std::is_move_assignable<WVStratifiedModalRecord>::value,"Published scientific records must not be moved over.");
namespace {
struct Counters { int plans=0,engines=0,created=0,failAt=-1,executed=0,failExecuteAt=-1; };
struct DerivativeAccessProbe {
    std::size_t cachedField=0,cachedDerivative=0,captures=0;
    std::size_t captureField=std::numeric_limits<std::size_t>::max();
    std::size_t captureDerivative=std::numeric_limits<std::size_t>::max();
    WVRealVolumeConstView cached{};
    std::vector<double>* captured=nullptr;
    bool failLookup=false,failCapture=false;
    static WVKernelStatus lookup(void* context,std::size_t field,std::size_t derivative,
        WVRealVolumeConstView& result) {
        auto& probe=*static_cast<DerivativeAccessProbe*>(context);
        if(probe.failLookup)
            return {WVKernelStatusCode::invalidConfiguration,"Injected derivative lookup failure."};
        result=field==probe.cachedField && derivative==probe.cachedDerivative ?
            probe.cached : WVRealVolumeConstView{};
        return WVKernelStatus::ok();
    }
    static WVKernelStatus capture(void* context,std::size_t field,std::size_t derivative,
        WVRealVolumeConstView values) {
        auto& probe=*static_cast<DerivativeAccessProbe*>(context);
        ++probe.captures;
        if(probe.captured && field==probe.captureField && derivative==probe.captureDerivative)
            probe.captured->assign(values.data,values.data+values.shape.first*values.shape.second*values.shape.third);
        return probe.failCapture ?
            WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                "Injected derivative capture failure."} : WVKernelStatus::ok();
    }
    WVStateDerivativeAccess access() {return {this,lookup,capture};}
};
class Plan final : public WVFFTPlan {
    std::unique_ptr<WVFFTPlan> plan_; Counters& counters_;
public:
    Plan(std::unique_ptr<WVFFTPlan> p,Counters& c):plan_(std::move(p)),counters_(c) { ++counters_.plans; }
    ~Plan() override { --counters_.plans; }
    WVKernelStatus execute(const void* a,void* b) override {
        if(counters_.executed++==counters_.failExecuteAt)
            return {WVKernelStatusCode::fftExecutionFailure,"Injected FFT execution failure."};
        return plan_->execute(a,b);
    }
    std::size_t persistentBytes() const noexcept override { return plan_->persistentBytes(); }
};
class Engine final : public WVFFTEngine {
    WVReferenceFFTEngine engine_; Counters& counters_;
public:
    explicit Engine(Counters& c):counters_(c) { ++counters_.engines; }
    ~Engine() override { --counters_.engines; }
    std::string identifier() const override { return "failure-probe"; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    WVKernelStatus createPlan(const WVFFTPlanSpecification& s,std::unique_ptr<WVFFTPlan>& p) override {
        if (counters_.created++==counters_.failAt) return {WVKernelStatusCode::fftPlanFailure,"Injected plan failure."};
        std::unique_ptr<WVFFTPlan> q; auto status=engine_.createPlan(s,q); if (!status) return status;
        p=std::make_unique<Plan>(std::move(q),counters_); return WVKernelStatus::ok();
    }
};
void contracts(const std::shared_ptr<const WVStratifiedModalRecord>& source) {
    Counters counters;
    std::unique_ptr<WVTransformHydrostaticKernel> kernel;
    require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<Engine>(counters),kernel)),"Create failed");
    const auto& g=source->geometry(); const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    const auto shape=kernel->spectralShape(); const auto volume=kernel->spatialShape();
    std::array<std::vector<WVComplex64>,3> a,b;
    for (auto& x:a) x.resize(S,{.001,-.002});
    for (auto& x:b) x.resize(S,{17,19});
    WVMutableCoefficients amplitudes{{a[0].data(),shape},{a[1].data(),shape},{a[2].data(),shape}};
    require(bool(kernel->constrainCoefficients(amplitudes)),"Constraints failed");
    WVState state{83,17,{{a[0].data(),shape},{a[1].data(),shape},{a[2].data(),shape}}};
    WVMutableCoefficients out{{b[0].data(),shape},{b[1].data(),shape},{b[2].data(),shape}};
    WVFlux flux{out.Ap,out.Am,out.A0};
    std::vector<double> spatial(R,29),scratch(R,31); WVRealVolumeView field{spatial.data(),volume};
    kernel->resetMetrics();
    int evaluationOwner=0,foreignOwner=0;
    require(bool(kernel->beginStateEvaluation(state,&evaluationOwner)),"Begin scoped state evaluation failed");
    require(kernel->metrics().stateValidationCount==1 && kernel->metrics().phasePreparationCount==1,
        "Scoped setup did not validate and prepare phase exactly once");
    require(!kernel->beginStateEvaluation(state),"Nested state evaluation was accepted");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::u,field)),"Scoped field reconstruction failed");
    WVComplexConstView scopedPhase;
    require(bool(kernel->preparedPhase(state,scopedPhase)) && scopedPhase.data!=nullptr && scopedPhase.shape.rows==shape.rows && scopedPhase.shape.columns==shape.columns,
        "Scoped prepared phase was unavailable");
    require(bool(kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},state.t,state.t0,out)),
        "Matching scoped projection was rejected");
    require(kernel->evolveCoefficients(state,amplitudes).code==WVKernelStatusCode::overlappingArrays,
        "Scoped evolution mutated the active immutable state");
    require(kernel->constrainCoefficients(amplitudes).code==WVKernelStatusCode::overlappingArrays,
        "Scoped constraints mutated the active immutable state");
    require(kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},state.t,state.t0,amplitudes).code==WVKernelStatusCode::overlappingArrays,
        "Scoped projection mutated the active immutable state");
    require(kernel->metrics().stateValidationCount==1 && kernel->metrics().phasePreparationCount==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::u)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1,
        "Scoped calls repeated preparation or lost reconstruction metrics");
    require(!kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},state.t+1,state.t0,out),
        "Mismatched scoped projection time was accepted");
    auto registeredCoefficients=a;
    auto foreignState=state;
    foreignState.coefficients={{registeredCoefficients[0].data(),shape},{registeredCoefficients[1].data(),shape},{registeredCoefficients[2].data(),shape}};
    require(!kernel->transformStateField(foreignState,WVHydrostaticField::u,field),"Foreign scoped state was accepted");
    std::vector<double> tendencyFields(3*R);
    WVRealFieldBundleView tendencyBundle{tendencyFields.data(),{g.Nx,g.Ny,g.Nz,3}};
    require(bool(kernel->transformCoefficientTendencyToUVEta(foreignState,tendencyBundle)),
        "Derived coefficient tendency was rejected by the active primary-state scope");
    require(kernel->metrics().stateValidationCount==1 &&
        kernel->metrics().derivedValidationCount==1 &&
        kernel->metrics().tendencyReconstructionCount==std::array<std::size_t,4>{1,1,1,0} &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::u)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1,
        "Derived tendency validation or production was attributed to the primary state");
    require(!kernel->addStateEvaluationView(foreignState,&foreignOwner),"Foreign evaluation owner registered a state view");
    require(bool(kernel->addStateEvaluationView(foreignState,&evaluationOwner,2)),"Additional immutable state view was rejected");
    const auto registeredStatus=kernel->transformStateField(foreignState,WVHydrostaticField::u,field);
    require(bool(registeredStatus),registeredStatus.message.c_str());
    require(kernel->metrics().stateValidationCount==2 && kernel->metrics().phasePreparationCount==1,
        "Additional state view skipped validation or repeated phase preparation");
    require(kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::u)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)][2]==1,
        "Registered Hydrostatic component production lost its component identity");
    require(!kernel->removeStateEvaluationView(state,&evaluationOwner,0),
        "Primary Hydrostatic state view was removed");
    require(!kernel->removeStateEvaluationView(foreignState,&foreignOwner,2),
        "Foreign owner removed a Hydrostatic state view");
    require(bool(kernel->removeStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Hydrostatic state view removal failed");
    require(!kernel->validateStateEvaluation(foreignState) &&
                !kernel->removeStateEvaluationView(foreignState,&evaluationOwner,2),
        "Removed Hydrostatic state view remained registered");
    require(bool(kernel->addStateEvaluationView(foreignState,&evaluationOwner,3)) &&
                bool(kernel->validateStateEvaluation(foreignState)),
        "Hydrostatic state view storage could not be re-registered");
    require(bool(kernel->endStateEvaluation()),"End scoped state evaluation failed");
    const auto savedScopedCoefficient=a[0][0]; a[0][0].real=std::numeric_limits<double>::infinity();
    require(!kernel->beginStateEvaluation(state),"New scope skipped validation for mutated same-pointer state");
    a[0][0]=savedScopedCoefficient;
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVHydrostaticField::u,field)),"Standalone field reconstruction failed");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::v,field)),"Second standalone field reconstruction failed");
    require(kernel->metrics().stateValidationCount==2 && kernel->metrics().phasePreparationCount==2,
        "Standalone calls did not perform fresh state preparation");
    std::vector<double> compoundReference(R),compoundFirst(R),compoundSecond(R),compoundResult(R);
    require(bool(kernel->beginStateEvaluation(state)),"Compound producer scope begin failed");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::zetaX,
        {compoundReference.data(),volume})),"Reference horizontal vorticity failed");
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVHydrostaticField::w,
        {compoundFirst.data(),volume},WVHydrostaticDerivative::y)),"Prepared w_y failed");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::v,
        {compoundSecond.data(),volume},WVHydrostaticDerivative::z)),"Prepared v_z failed");
    require(bool(kernel->combinePreparedHorizontalVorticity(WVHydrostaticField::zetaX,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume})),"Prepared horizontal vorticity combination failed");
    require(compoundResult==compoundReference &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::w)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::y)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::v)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::z)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::zetaX)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1,
        "Prepared horizontal vorticity changed values or repeated a producer");
    require(bool(kernel->combinePreparedHorizontalVorticity(WVHydrostaticField::zetaX,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundFirst.data(),volume},WVHydrostaticComponent::geostrophic)) &&
        compoundFirst==compoundReference &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::zetaX)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
            [static_cast<std::size_t>(WVHydrostaticComponent::geostrophic)]==1,
        "Exact prepared vorticity alias or component attribution failed");
    std::vector<double> partialCompound(R+1);
    std::copy_n(compoundReference.data(),R,partialCompound.data());
    require(kernel->combinePreparedHorizontalVorticity(WVHydrostaticField::zetaX,
        {partialCompound.data(),volume},{compoundSecond.data(),volume},
        {partialCompound.data()+1,volume}).code==WVKernelStatusCode::overlappingArrays,
        "Partial prepared vorticity alias was accepted");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::rhoTotal,
        {compoundReference.data(),volume},WVHydrostaticDerivative::z)),"Reference density derivative failed");
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVHydrostaticField::eta,
        {compoundFirst.data(),volume},WVHydrostaticDerivative::z)),"Prepared eta_z failed");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::eta,
        {compoundSecond.data(),volume})),"Prepared eta failed");
    require(bool(kernel->combinePreparedDensityZDerivative(WVHydrostaticField::rhoTotal,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume})),"Prepared density derivative combination failed");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundResult[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Prepared density derivative changed values");
    require(kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::eta)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::z)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::eta)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVHydrostaticField::rhoTotal)]
            [static_cast<std::size_t>(WVHydrostaticDerivative::z)]
            [static_cast<std::size_t>(WVHydrostaticComponent::all)]==1,
        "Prepared density derivative repeated a producer");
    require(bool(kernel->combinePreparedDensityZDerivative(WVHydrostaticField::rhoTotal,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundFirst.data(),volume})),
        "Exact prepared density alias was rejected");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundFirst[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Exact prepared density alias changed values");
    std::copy_n(compoundReference.data(),R,partialCompound.data());
    require(kernel->combinePreparedDensityZDerivative(WVHydrostaticField::rhoTotal,
        {partialCompound.data(),volume},{compoundSecond.data(),volume},
        {partialCompound.data()+1,volume}).code==WVKernelStatusCode::overlappingArrays,
        "Partial prepared density alias was accepted");
    require(!kernel->combinePreparedDensityZDerivative(WVHydrostaticField::u,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume}),"Invalid prepared density target was accepted");
    require(bool(kernel->endStateEvaluation()),"Compound producer scope end failed");
    std::fill(spatial.begin(),spatial.end(),29);
    allocationProbe::calls=0; allocationProbe::counting=true;
    require(bool(kernel->beginStateEvaluation(state)),"Allocation probe scope begin failed");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::u,field)),"Allocation probe scoped field failed");
    require(bool(kernel->endStateEvaluation()),"Allocation probe scope end failed");
    allocationProbe::counting=false;
    require(allocationProbe::calls==0,"Scoped state preparation or reuse allocated");
    std::fill(spatial.begin(),spatial.end(),29);
    auto status=kernel->nonlinearFlux(state,flux); require(bool(status),status.message.c_str());
    const auto referenceFlux=b,referenceState=a;
    const auto retainedBytes=kernel->persistentBytes();
    std::vector<double> raw(3*R),physical(4*R),gradient(R),expected(3*R,0);
    WVRealFieldBundleView rawView{raw.data(),{g.Nx,g.Ny,g.Nz,3}};
    const WVHydrostaticField dynamical[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::w,WVHydrostaticField::eta};
    for(std::size_t channel=0;channel<4;++channel)
        require(bool(kernel->transformStateField(state,dynamical[channel],{physical.data()+channel*R,volume})),"Prepare shared physical fields");
    const WVRealFieldBundleConstView prepared{physical.data(),{g.Nx,g.Ny,g.Nz,4}};
    const auto physicalBefore=physical;
    const auto fullStart=counters.executed;
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared)),"Observe raw nonlinear tendency");
    const auto fullExecutions=counters.executed-fullStart;
    for(std::size_t channel=0;channel<3;++channel) {
        const auto target=channel==2 ? 3 : channel;
        for(std::size_t axis=0;axis<3;++axis) {
            require(bool(kernel->transformStateField(state,dynamical[target],{gradient.data(),volume},static_cast<WVHydrostaticDerivative>(axis+1))),"Independent tendency derivative");
            for(std::size_t i=0;i<R;++i) {
                const double correction=target==3 && axis==2 ? physical[3*R+i]*g.dLnN2[i/(g.Nx*g.Ny)] : 0;
                expected[channel*R+i]-=physical[axis*R+i]*(gradient[i]+correction);
            }
        }
    }
    require(raw==expected && physical==physicalBefore,"Raw tendencies precede projection and preserve borrowed fields");
    require(bool(kernel->transformStateField(state,WVHydrostaticField::u,
        {gradient.data(),volume},WVHydrostaticDerivative::z)),
        "Prepare cached Hydrostatic derivative");
    kernel->resetMetrics();
    DerivativeAccessProbe derivativeProbe;
    derivativeProbe.cachedField=static_cast<std::size_t>(WVHydrostaticField::u);
    derivativeProbe.cachedDerivative=static_cast<std::size_t>(WVHydrostaticDerivative::z);
    derivativeProbe.cached={gradient.data(),volume};
    auto derivativeAccess=derivativeProbe.access();
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess)) && raw==expected,
        "Hydrostatic nonlinear derivative reuse changed the spatial tendency");
    require(kernel->metrics().reconstructionCount
                [static_cast<std::size_t>(WVHydrostaticField::u)]
                [static_cast<std::size_t>(WVHydrostaticDerivative::z)]
                [static_cast<std::size_t>(WVHydrostaticComponent::all)]==0 &&
            derivativeProbe.captures==8,
        "Hydrostatic nonlinear path reproduced a cached derivative or missed capture");
    derivativeProbe.failLookup=true;
    require(!kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess),
        "Hydrostatic nonlinear derivative lookup failure was ignored");
    derivativeProbe.failLookup=false;
    derivativeProbe.failCapture=true;
    derivativeProbe.cachedField=99;
    require(!kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess),
        "Hydrostatic nonlinear derivative capture failure was ignored");
    for(std::size_t channel=0;channel<3;++channel) for(std::size_t i=0;i<S;++i)
        require(b[channel][i].real==referenceFlux[channel][i].real && b[channel][i].imag==referenceFlux[channel][i].imag &&
                a[channel][i].real==referenceState[channel][i].real && a[channel][i].imag==referenceState[channel][i].imag,"Observation changed flux or input state");
    for (auto& values:b) std::fill(values.begin(),values.end(),WVComplex64{17,19});
    const auto rawStart=counters.executed;
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared,false)),"Spatial-only nonlinear tendency");
    require(raw==expected && physical==physicalBefore && counters.executed-rawStart<fullExecutions,
        "Spatial-only evaluation changed the tendency or retained redundant projection work");
    for(const auto& values:b) for(const auto value:values)
        require(value.real==17 && value.imag==19,"Spatial-only evaluation wrote spectral flux");
    const auto rejectedStart=counters.executed;
    require(kernel->nonlinearFlux(state,flux,nullptr,&prepared,false).code==WVKernelStatusCode::invalidConfiguration &&
        counters.executed==rejectedStart,"Spatial-only evaluation accepted missing output or rejected after FFT work");
    auto badRaw=rawView; badRaw.data=physical.data();
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::overlappingArrays,"Observed tendency aliases shared fields");
    badRaw=rawView; badRaw.shape.fourth=4;
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::invalidShape,"Wrong tendency channel count accepted");
    badRaw=rawView; badRaw.data=nullptr;
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::invalidPointer,"Null tendency output accepted");
    require(raw==expected && physical==physicalBefore,"Invalid observation mutated output");
    allocationProbe::calls=0; allocationProbe::counting=true;
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared)),"Prepared observation failed");
    allocationProbe::counting=false;
    require(allocationProbe::calls==0 && kernel->persistentBytes()==retainedBytes,"Observation allocated persistent storage");
    for (auto& x:b) std::fill(x.begin(),x.end(),WVComplex64{17,19});
    auto badFlux=flux; badFlux.Fp=amplitudes.Ap;
    status=kernel->nonlinearFlux(state,badFlux); require(status.code==WVKernelStatusCode::overlappingArrays,"Aliased flux accepted");
    badFlux=flux; badFlux.Fm=flux.Fp;
    status=kernel->nonlinearFlux(state,badFlux); require(status.code==WVKernelStatusCode::overlappingArrays,"Overlapping outputs accepted");
    badFlux=flux; badFlux.F0.shape={S,1};
    status=kernel->nonlinearFlux(state,badFlux); require(status.code==WVKernelStatusCode::invalidShape,"Bad shape accepted");
    badFlux=flux; badFlux.Fm.data=nullptr;
    status=kernel->nonlinearFlux(state,badFlux); require(status.code==WVKernelStatusCode::invalidPointer,"Null output accepted");
    auto invalid=state; invalid.t=std::numeric_limits<double>::quiet_NaN();
    require(!kernel->nonlinearFlux(invalid,flux),"NaN time accepted");
    invalid=state; invalid.t=std::numeric_limits<double>::max(); invalid.t0=-invalid.t;
    require(!kernel->evolveCoefficients(invalid,out),"Elapsed-time overflow accepted");
    require(!kernel->transformStateField(state,static_cast<WVHydrostaticField>(99),field),"Invalid field accepted");
    require(!kernel->transformStateField(state,WVHydrostaticField::u,field,static_cast<WVHydrostaticDerivative>(99)),"Invalid derivative accepted");
    require(!kernel->transformStateField(state,WVHydrostaticField::u,field,WVHydrostaticDerivative::value,static_cast<WVHydrostaticComponent>(99)),"Invalid component accepted");
    require(!kernel->transformStateField(state,WVHydrostaticField::zetaX,field,WVHydrostaticDerivative::x),"Unsupported vorticity derivative accepted");
    require(!kernel->transformStateField(state,WVHydrostaticField::rhoTotal,field,WVHydrostaticDerivative::value,WVHydrostaticComponent::wave),"Masked background accepted");
    require(!kernel->transformStateField(state,WVHydrostaticField::ssh,field),"Wrong surface shape accepted");
    require(!kernel->differentiateVertical({scratch.data(),volume},WVHydrostaticFamily::F,0,field),"Order zero accepted");
    require(!kernel->differentiateVertical({scratch.data(),volume},WVHydrostaticFamily::F,5,field),"Order five accepted");
    require(!kernel->integrateVertical({scratch.data(),volume},static_cast<WVHydrostaticFamily>(99),field),"Invalid family accepted");
    require(!kernel->differentiateVertical({spatial.data()+1,volume},WVHydrostaticFamily::G,1,field),"Partial calculus alias accepted");
    auto aliased=out; aliased.Ap=amplitudes.Am;
    require(!kernel->evolveCoefficients(state,aliased),"Cross-family evolution alias accepted");
    const auto saved=a;
    a[2][S-1].real=std::numeric_limits<double>::infinity();
    require(!kernel->transformStateField(state,WVHydrostaticField::u,field),"Infinite coefficient accepted");
    double energy=71;
    require(!kernel->totalEnergy(state.coefficients,energy) && energy==71,"Invalid energy mutated scalar");
    a[2][S-1]=saved[2][S-1];
    for (auto& x:b) for (auto value:x) require(value.real==17 && value.imag==19,"Invalid call mutated coefficients");
    for (auto value:spatial) require(value==29,"Invalid call mutated field");
    require(bool(kernel->evolveCoefficients(state,out)),"Evolution failed");
    require(bool(kernel->evolveCoefficients(state,amplitudes)),"In-place evolution failed");
    for (int j=0;j<3;++j) for (std::size_t i=0;i<S;++i) require(a[j][i].real==b[j][i].real && a[j][i].imag==b[j][i].imag,"In-place evolution differs");
    const auto bytes=kernel->persistentBytes();
    allocationProbe::calls=0; allocationProbe::counting=true;
    for (int i=0;i<5;++i) {
        require(bool(kernel->nonlinearFlux(state,flux)),"Prepared flux failed");
        require(bool(kernel->transformStateField(state,WVHydrostaticField::rhoTotal,field,WVHydrostaticDerivative::z)),"Prepared density failed");
        require(bool(kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},83,17,out)),"Prepared projection failed");
        require(bool(kernel->differentiateVertical({scratch.data(),volume},WVHydrostaticFamily::F,4,field)),"Prepared derivative failed");
        require(bool(kernel->integrateVertical({scratch.data(),volume},WVHydrostaticFamily::G,field)),"Prepared integral failed");
    }
    allocationProbe::counting=false; require(allocationProbe::calls==0 && kernel->persistentBytes()==bytes,"Prepared execution allocates or changes storage");
    const auto* old=kernel.get();
    status=WVTransformHydrostaticKernel::create({},std::make_unique<WVReferenceFFTEngine>(),kernel);
    require(!status && kernel.get()==old,"Failed setup replaced kernel");
    for (int fail=0;fail<4;++fail) {
        int count=0; Counters counters;
        status=WVTransformHydrostaticKernel::create(source,std::make_unique<Engine>(counters),kernel,[&](std::unique_ptr<WVVerticalMatrixBackend>& backend) {
            if (count++==fail) return WVKernelStatus{WVKernelStatusCode::allocationFailure,"Injected matrix factory failure."};
            return WVCreateScalarMatrixBackend(backend);
        });
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"Matrix setup failure leaked or replaced kernel");
    }
    for (int fail=0;fail<2;++fail) {
        Counters counters; counters.failAt=fail;
        status=WVTransformHydrostaticKernel::create(source,std::make_unique<Engine>(counters),kernel);
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"FFT setup failure leaked or replaced kernel");
    }
    bool succeeded=false;
    for (long fail=0;fail<2048;++fail) {
        Counters counters; auto engine=std::make_unique<Engine>(counters);
        allocationProbe::failAfter=fail;
        status=WVTransformHydrostaticKernel::create(source,std::move(engine),kernel);
        allocationProbe::failAfter=-1;
        if (status) { kernel.reset(); require(!counters.engines && !counters.plans,"Successful destruction leaked"); succeeded=true; break; }
        require(kernel.get()==old && !counters.engines && !counters.plans,"Allocation failure leaked or replaced kernel");
    }
    require(succeeded,"Allocation sweep never succeeded");
}

void variableScheduleParity(const std::shared_ptr<const WVStratifiedModalRecord>& source,bool compact = false) {
    Counters candidateCounters; // Outlives the instrumented plans owned by candidate.
    std::unique_ptr<WVTransformHydrostaticKernel> frozen, candidate;
    require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),frozen)),"Frozen schedule setup failed");
    WVVariableExecutionOptions options{WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true};
    if (compact) options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
    options.pointwiseWorkers=2;
    options.fusedDerivativeAdvection=true;
    require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<Engine>(candidateCounters),candidate, WVCreateScalarMatrixBackend, options)),"Candidate schedule setup failed");
    require(std::string(candidate->horizontalScheduleIdentifier())=="full-fft-gather","Reference provider fallback was not reported");
    const auto& g=source->geometry(); const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    std::array<std::vector<WVComplex64>,3> input,frozenOut,candidateOut;
    for (std::size_t j=0;j<3;++j) { input[j].resize(S); frozenOut[j].resize(S); candidateOut[j].resize(S); for (std::size_t i=0;i<S;++i) input[j][i]={.001*std::sin(.17*i+j),-.002*std::cos(.11*i+j)}; }
    WVMutableCoefficients validInput{{input[0].data(),{g.Nj,g.Nkl}},{input[1].data(),{g.Nj,g.Nkl}},{input[2].data(),{g.Nj,g.Nkl}}};
    require(bool(frozen->constrainCoefficients(validInput)),"Parity input constraints failed");
    const auto inputBefore=input;
    const WVState state{83,17,{{input[0].data(),{g.Nj,g.Nkl}},{input[1].data(),{g.Nj,g.Nkl}},{input[2].data(),{g.Nj,g.Nkl}}}};
    WVFlux frozenFlux{{frozenOut[0].data(),{g.Nj,g.Nkl}},{frozenOut[1].data(),{g.Nj,g.Nkl}},{frozenOut[2].data(),{g.Nj,g.Nkl}}};
    WVFlux candidateFlux{{candidateOut[0].data(),{g.Nj,g.Nkl}},{candidateOut[1].data(),{g.Nj,g.Nkl}},{candidateOut[2].data(),{g.Nj,g.Nkl}}};
    require(bool(frozen->nonlinearFlux(state,frozenFlux)),"Frozen nonlinear flux failed");
    require(bool(candidate->nonlinearFlux(state,candidateFlux)),"Candidate nonlinear flux failed");
    require(!frozen->executionOptions().fusedDerivativeAdvection &&
        candidate->executionOptions().fusedDerivativeAdvection &&
        candidate->metrics().derivativeAdvectionConsumerCount==9,
        "Hydrostatic fused derivative consumers were not explicit or complete");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i) require(std::abs(frozenOut[j][i].real-candidateOut[j][i].real)<1e-12 && std::abs(frozenOut[j][i].imag-candidateOut[j][i].imag)<1e-12,"Candidate nonlinear flux differs");
    auto serialOptions=options; serialOptions.pointwiseWorkers=1;
    std::unique_ptr<WVTransformHydrostaticKernel> serial; std::array<std::vector<WVComplex64>,3> serialOut;
    for (auto& values:serialOut) values.resize(S);
    require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),serial,
        WVCreateScalarMatrixBackend,serialOptions)),"Single-worker candidate setup failed");
    WVFlux serialFlux{{serialOut[0].data(),{g.Nj,g.Nkl}},{serialOut[1].data(),{g.Nj,g.Nkl}},{serialOut[2].data(),{g.Nj,g.Nkl}}};
    require(bool(serial->nonlinearFlux(state,serialFlux)),"Single-worker candidate flux failed");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i)
        require(serialOut[j][i].real==candidateOut[j][i].real && serialOut[j][i].imag==candidateOut[j][i].imag,
            "Pointwise worker partition changed hydrostatic arithmetic");
    // Phases are recomputed from the current time pair, including repeated,
    // reversed and extreme finite times, without recurrence or new storage.
    for (double time:{17.0,83.0,-83.0,1e12,1e300,83.0,17.0}) {
        auto timedState=state; timedState.t=time;
        for (auto* target:{serial.get(),candidate.get()}) {
            require(bool(target->beginStateEvaluation(timedState)),"Phase partition scope failed");
            WVComplexConstView phases;
            require(bool(target->preparedPhase(timedState,phases)),"Prepared phases were unavailable");
            for (std::size_t i=0;i<S;++i) {
                const double angle=target->factors()[i].omega*(time-timedState.t0);
                require(phases.data[i].real==std::cos(angle) && phases.data[i].imag==std::sin(angle),
                    "Phase workers changed scalar trigonometry or reused an old time");
            }
            require(bool(target->endStateEvaluation()),"Phase partition scope release failed");
        }
    }
    // Compare reduction with the serial hypot/max oracle, including magnitudes
    // for which an unscaled square/sum/sqrt would overflow or underflow.
    std::vector<double> speedU(R),speedV(R);
    const double scales[]={1.0,1e308,std::numeric_limits<double>::denorm_min()};
    for (double magnitude:scales) {
        for (std::size_t i=0;i<R;++i) {
            speedU[i]=magnitude*(i%2 ? 1.0 : -.5);
            speedV[i]=magnitude*(i%3 ? -.75 : .25);
        }
        double expected=0;
        for (std::size_t i=0;i<R;++i) expected=std::max(expected,std::hypot(speedU[i],speedV[i]));
        for (auto* target:{serial.get(),candidate.get()}) {
            double actual=-1;
            allocationProbe::calls=0; allocationProbe::counting=true;
            auto status=target->reduceHorizontalSpeedMaximum(
                {speedU.data(),{g.Nx,g.Ny,g.Nz}},{speedV.data(),{g.Nx,g.Ny,g.Nz}},actual);
            allocationProbe::counting=false;
            require(bool(status) && actual==expected,"Speed reduction changed robust hypot arithmetic");
            require(allocationProbe::calls==0,"Prepared speed reduction allocated");
        }
    }
    std::fill(speedU.begin(),speedU.end(),std::numeric_limits<double>::quiet_NaN());
    std::fill(speedV.begin(),speedV.end(),0.0);
    double maximum=-1;
    require(bool(candidate->reduceHorizontalSpeedMaximum({speedU.data(),{g.Nx,g.Ny,g.Nz}},
        {speedV.data(),{g.Nx,g.Ny,g.Nz}},maximum)) && maximum==0,
        "NaN speed contributions changed the existing max policy");
    speedU.back()=std::numeric_limits<double>::infinity();
    require(bool(candidate->reduceHorizontalSpeedMaximum({speedU.data(),{g.Nx,g.Ny,g.Nz}},
        {speedV.data(),{g.Nx,g.Ny,g.Nz}},maximum)) && std::isinf(maximum),"Infinite speed was lost");
    maximum=17;
    require(!candidate->reduceHorizontalSpeedMaximum({nullptr,{g.Nx,g.Ny,g.Nz}},
        {speedV.data(),{g.Nx,g.Ny,g.Nz}},maximum) && maximum==17,"Invalid speed input published a result");
    // Exercise every field dispatch and component with the same spectral layout
    // on one and two workers, including surface aliases and compound fields.
    std::vector<double> serialField(R),parallelField(R);
    for (int field=static_cast<int>(WVHydrostaticField::u);field<=static_cast<int>(WVHydrostaticField::ssv);++field) {
        const auto name=static_cast<WVHydrostaticField>(field);
        const bool compound=name==WVHydrostaticField::zetaX || name==WVHydrostaticField::zetaY;
        for (int derivative=0;derivative<(compound ? 1 : 4);++derivative) for (int component=0;component<5;++component) {
            if (name==WVHydrostaticField::rhoTotal && component!=0) continue;
            const auto d=static_cast<WVHydrostaticDerivative>(derivative);
            const auto c=static_cast<WVHydrostaticComponent>(component);
            const WVShape3D shape{g.Nx,g.Ny,name>=WVHydrostaticField::ssh ? 1 : g.Nz};
            require(bool(serial->transformStateField(state,name,{serialField.data(),shape},d,c)),"Serial field partition reference failed");
            require(bool(candidate->transformStateField(state,name,{parallelField.data(),shape},d,c)),"Parallel field partition failed");
            require(serialField==parallelField,"Coefficient worker partition changed field arithmetic");
        }
    }
    std::vector<double> ff(4*R),cf(4*R),fr(3*R),cr(3*R);
    const WVHydrostaticField fields[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::w,WVHydrostaticField::eta};
    for (std::size_t j=0;j<4;++j) { require(bool(frozen->transformStateField(state,fields[j],{ff.data()+j*R,{g.Nx,g.Ny,g.Nz}})),"Frozen borrowed fields failed"); require(bool(candidate->transformStateField(state,fields[j],{cf.data()+j*R,{g.Nx,g.Ny,g.Nz}})),"Candidate borrowed fields failed"); }
    const auto ffBefore=ff,cfBefore=cf;
    WVRealFieldBundleConstView ffields{ff.data(),{g.Nx,g.Ny,g.Nz,4}},cfields{cf.data(),{g.Nx,g.Ny,g.Nz,4}}; WVRealFieldBundleView frv{fr.data(),{g.Nx,g.Ny,g.Nz,3}},crv{cr.data(),{g.Nx,g.Ny,g.Nz,3}};
    require(bool(frozen->nonlinearFlux(state,frozenFlux,&frv,&ffields)),"Frozen raw flux failed"); allocationProbe::calls=0; allocationProbe::counting=true; require(bool(candidate->nonlinearFlux(state,candidateFlux,&crv,&cfields)),"Candidate raw flux failed"); allocationProbe::counting=false; require(allocationProbe::calls==0,"Candidate prepared nonlinear flux allocated");
    require(ff==ffBefore && cf==cfBefore,"Candidate mutated borrowed fields");
    for (std::size_t i=0;i<fr.size();++i) require(std::abs(fr[i]-cr[i])<1e-12,"Candidate raw tendency differs");
    for (std::size_t i=0;i<ff.size();++i) require(std::abs(ff[i]-cf[i])<1e-12,"Candidate reconstructed field differs");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i) require(std::abs(frozenOut[j][i].real-candidateOut[j][i].real)<1e-12 && std::abs(frozenOut[j][i].imag-candidateOut[j][i].imag)<1e-12,"Candidate prepared nonlinear flux differs");
    for (auto* outputs:{&frozenOut,&candidateOut}) for (auto& values:*outputs) std::fill(values.begin(),values.end(),WVComplex64{17,19});
    require(bool(frozen->nonlinearFlux(state,frozenFlux,&frv,&ffields,false)),"Frozen spatial-only flux failed");
    const auto uncachedStart=candidateCounters.executed;
    require(bool(candidate->nonlinearFlux(state,candidateFlux,&crv,&cfields,false)),"Candidate spatial-only flux failed");
    const auto uncachedExecutions=candidateCounters.executed-uncachedStart;
    require(ff==ffBefore && cf==cfBefore,"Candidate spatial-only evaluation mutated borrowed fields");
    for (std::size_t i=0;i<fr.size();++i) require(std::abs(fr[i]-cr[i])<1e-12,"Candidate spatial-only tendency differs");
    for (const auto* outputs:{&frozenOut,&candidateOut}) for (const auto& values:*outputs) for (const auto value:values)
        require(value.real==17 && value.imag==19,"Spatial-only evaluation wrote spectral flux");
    std::vector<double> cachedDerivative(R),captureReference(R),capturedDerivative;
    require(bool(candidate->transformStateField(state,WVHydrostaticField::u,
        {cachedDerivative.data(),{g.Nx,g.Ny,g.Nz}},WVHydrostaticDerivative::z)) &&
        bool(candidate->transformStateField(state,WVHydrostaticField::u,
        {captureReference.data(),{g.Nx,g.Ny,g.Nz}},WVHydrostaticDerivative::x)),
        "Hydrostatic fused derivative cache fixtures failed");
    DerivativeAccessProbe derivativeProbe;
    derivativeProbe.cachedField=static_cast<std::size_t>(WVHydrostaticField::u);
    derivativeProbe.cachedDerivative=static_cast<std::size_t>(WVHydrostaticDerivative::z);
    derivativeProbe.cached={cachedDerivative.data(),{g.Nx,g.Ny,g.Nz}};
    derivativeProbe.captureField=static_cast<std::size_t>(WVHydrostaticField::u);
    derivativeProbe.captureDerivative=static_cast<std::size_t>(WVHydrostaticDerivative::x);
    derivativeProbe.captured=&capturedDerivative;
    auto derivativeAccess=derivativeProbe.access();
    candidate->resetMetrics();
    const auto cachedStart=candidateCounters.executed;
    require(bool(candidate->nonlinearFlux(state,candidateFlux,&crv,&cfields,false,&derivativeAccess)),
        "Hydrostatic fused cached-derivative evaluation failed");
    require(candidateCounters.executed-cachedStart+1==uncachedExecutions &&
        candidate->metrics().derivativeAdvectionConsumerCount==8 &&
        derivativeProbe.captures==8 && capturedDerivative==captureReference &&
        cfields.data==cf.data() && cf==cfBefore,
        "Hydrostatic cached derivative repeated an inverse/consumer or lost capture");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i)
        require(input[j][i].real==inputBefore[j][i].real && input[j][i].imag==inputBefore[j][i].imag,"Parity evaluation mutated input coefficients");
    require(frozen->storage().realScratchBytes==10*R*sizeof(double) && candidate->storage().realScratchBytes==6*R*sizeof(double),"Hydrostatic streamed scratch accounting differs");
    require(frozen->storage().spectralScratchBytes==candidate->storage().spectralScratchBytes,"Hydrostatic streamed spectral scratch changed");
    require(candidate->executionOptions().horizontalWorkers==2 && candidate->executionOptions().streamedNonlinear &&
        candidate->executionOptions().pointwiseWorkers==2 && candidate->executionOptions().fusedDerivativeAdvection,
        "Candidate options were not retained");
    auto* retained=candidate.get(); WVVariableExecutionOptions invalidOptions; invalidOptions.pointwiseWorkers=0;
    const auto status=WVTransformHydrostaticKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),candidate,
        WVCreateScalarMatrixBackend,invalidOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && candidate.get()==retained,
        "Zero pointwise workers were accepted or replaced the retained kernel");
}

void sharedFieldGradientParity(const std::shared_ptr<const WVStratifiedModalRecord>& source,
    bool compact) {
    WVVariableExecutionOptions sharedOptions{WVRetainedHorizontalSchedule::fullFFT,2,true};
    sharedOptions.pointwiseWorkers=2;
    sharedOptions.fusedDerivativeAdvection=true;
    if (compact) {
        sharedOptions={WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true};
        sharedOptions.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
        sharedOptions.pointwiseWorkers=2;
        sharedOptions.fusedDerivativeAdvection=true;
    }
    auto independentOptions=sharedOptions;
    independentOptions.sharedFieldGradients=false;
    std::unique_ptr<WVTransformHydrostaticKernel> shared,independent;
    require(bool(WVTransformHydrostaticKernel::create(source,
        std::make_unique<WVReferenceFFTEngine>(),shared,WVCreateScalarMatrixBackend,
        sharedOptions)),"Shared field-gradient kernel setup failed");
    require(bool(WVTransformHydrostaticKernel::create(source,
        std::make_unique<WVReferenceFFTEngine>(),independent,WVCreateScalarMatrixBackend,
        independentOptions)),"Independent field-gradient kernel setup failed");

    const auto& g=source->geometry();
    const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    const WVShape2D spectralShape{g.Nj,g.Nkl};
    const WVShape3D spatialShape{g.Nx,g.Ny,g.Nz};
    std::array<std::vector<WVComplex64>,3> coefficients;
    for(std::size_t family=0;family<coefficients.size();++family) {
        coefficients[family].resize(S);
        for(std::size_t i=0;i<S;++i)
            coefficients[family][i]={.003*std::sin(.13*i+family),
                -.002*std::cos(.19*i+2*family)};
    }
    WVMutableCoefficients mutableCoefficients{{coefficients[0].data(),spectralShape},
        {coefficients[1].data(),spectralShape},{coefficients[2].data(),spectralShape}};
    require(bool(shared->constrainCoefficients(mutableCoefficients)),
        "Shared field-gradient input constraints failed");
    const auto coefficientsBefore=coefficients;
    const WVState state{83,17,{{coefficients[0].data(),spectralShape},
        {coefficients[1].data(),spectralShape},{coefficients[2].data(),spectralShape}}};
    std::vector<double> sharedField(R),independentField(R);
    require(bool(shared->beginStateEvaluation(state)),
        "Shared field-gradient scope begin failed");
    require(bool(independent->beginStateEvaluation(state)),
        "Independent field-gradient scope begin failed");
    shared->resetMetrics();
    independent->resetMetrics();
    struct FieldComponent { WVHydrostaticField field; WVHydrostaticComponent component; };
    const FieldComponent keys[]={
        {WVHydrostaticField::u,WVHydrostaticComponent::all},
        {WVHydrostaticField::u,WVHydrostaticComponent::wave},
        {WVHydrostaticField::v,WVHydrostaticComponent::all},
        {WVHydrostaticField::v,WVHydrostaticComponent::inertial},
        {WVHydrostaticField::w,WVHydrostaticComponent::all},
        {WVHydrostaticField::w,WVHydrostaticComponent::geostrophic},
        {WVHydrostaticField::eta,WVHydrostaticComponent::all},
        {WVHydrostaticField::eta,WVHydrostaticComponent::meanDensityAnomaly}};
    for(const auto key:keys) for(std::size_t derivative=0;derivative<4;++derivative) {
            const auto d=static_cast<WVHydrostaticDerivative>(derivative);
            require(bool(shared->transformStateField(state,key.field,
                {sharedField.data(),spatialShape},d,key.component)),
                "Shared field-gradient reconstruction failed");
            require(bool(independent->transformStateField(state,key.field,
                {independentField.data(),spatialShape},d,key.component)),
                "Independent field-gradient reconstruction failed");
            for(std::size_t i=0;i<R;++i)
                require(std::abs(sharedField[i]-independentField[i])<=1e-12*
                    std::max(1.0,std::abs(independentField[i])),
                    "Shared field-gradient reconstruction differs from the independent path");
        }
    const auto sharedFieldMetrics=shared->metrics();
    const auto independentFieldMetrics=independent->metrics();
    require(sharedFieldMetrics.coefficientAssemblyCount==8 &&
            sharedFieldMetrics.verticalPreparationCount==16 &&
            sharedFieldMetrics.verticalOperatorExecutionCount==16 &&
            sharedFieldMetrics.horizontalSpectrumReuseCount==16 &&
            sharedFieldMetrics.preparedVerticalDerivativeCount==8,
        "Shared field-gradient producers did not execute once per field/component key");
    require(independentFieldMetrics.coefficientAssemblyCount==32 &&
            independentFieldMetrics.verticalPreparationCount==32 &&
            independentFieldMetrics.verticalOperatorExecutionCount==32 &&
            independentFieldMetrics.horizontalSpectrumReuseCount==0 &&
            independentFieldMetrics.preparedVerticalDerivativeCount==0,
        "Independent field-gradient qualification path unexpectedly reused preparation");
    require(bool(shared->endStateEvaluation()),"Shared field-gradient scope end failed");
    require(bool(independent->endStateEvaluation()),
        "Independent field-gradient scope end failed");

    // Removing and re-adding the exact coefficient storage must assign a new
    // registered-view generation. Reuse within either generation remains valid.
    auto foreignCoefficients=coefficients;
    WVState foreignState=state;
    foreignState.coefficients={{foreignCoefficients[0].data(),spectralShape},
        {foreignCoefficients[1].data(),spectralShape},
        {foreignCoefficients[2].data(),spectralShape}};
    int evaluationOwner=0;
    require(bool(shared->beginStateEvaluation(state,&evaluationOwner)),
        "Registered-view field-gradient scope begin failed");
    require(bool(shared->addStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Registered-view field-gradient setup failed");
    shared->resetMetrics();
    require(bool(shared->transformStateField(foreignState,WVHydrostaticField::u,
        {sharedField.data(),spatialShape})),"First registered-view preparation failed");
    require(bool(shared->transformStateField(foreignState,WVHydrostaticField::u,
        {sharedField.data(),spatialShape},WVHydrostaticDerivative::x)),
        "Registered-view horizontal reuse failed");
    require(bool(shared->removeStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Prepared registered-view removal failed");
    require(bool(shared->addStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Same-storage registered-view re-add failed");
    require(bool(shared->transformStateField(foreignState,WVHydrostaticField::u,
        {sharedField.data(),spatialShape})),"Re-added registered-view preparation failed");
    require(shared->metrics().coefficientAssemblyCount==2 &&
            shared->metrics().verticalPreparationCount==2 &&
            shared->metrics().horizontalSpectrumReuseCount==1,
        "Re-added same-storage view reused an invalid generation or lost in-view reuse");
    require(bool(shared->endStateEvaluation()),
        "Registered-view field-gradient scope end failed");

    shared->resetMetrics();
    require(bool(shared->transformStateField(state,WVHydrostaticField::u,
        {sharedField.data(),spatialShape})),"First standalone shared preparation failed");
    require(bool(shared->transformStateField(state,WVHydrostaticField::u,
        {sharedField.data(),spatialShape})),"Second standalone shared preparation failed");
    require(shared->metrics().coefficientAssemblyCount==2 &&
            shared->metrics().verticalPreparationCount==2,
        "Standalone operations published shared preparation across call boundaries");

    std::array<std::vector<WVComplex64>,3> sharedFluxStorage,independentFluxStorage;
    for(auto& values:sharedFluxStorage) values.resize(S);
    for(auto& values:independentFluxStorage) values.resize(S);
    WVFlux sharedFlux{{sharedFluxStorage[0].data(),spectralShape},
        {sharedFluxStorage[1].data(),spectralShape},{sharedFluxStorage[2].data(),spectralShape}};
    WVFlux independentFlux{{independentFluxStorage[0].data(),spectralShape},
        {independentFluxStorage[1].data(),spectralShape},
        {independentFluxStorage[2].data(),spectralShape}};
    std::vector<double> sharedRaw(3*R),independentRaw(3*R);
    WVRealFieldBundleView sharedRawView{sharedRaw.data(),{g.Nx,g.Ny,g.Nz,3}};
    WVRealFieldBundleView independentRawView{independentRaw.data(),{g.Nx,g.Ny,g.Nz,3}};
    require(bool(shared->beginStateEvaluation(state)),
        "Shared nonlinear preparation scope begin failed");
    require(bool(independent->beginStateEvaluation(state)),
        "Independent nonlinear preparation scope begin failed");
    shared->resetMetrics();
    independent->resetMetrics();
    require(bool(shared->nonlinearFlux(state,sharedFlux,&sharedRawView)),
        "Shared nonlinear preparation failed");
    require(bool(independent->nonlinearFlux(state,independentFlux,&independentRawView)),
        "Independent nonlinear preparation failed");
    for(std::size_t i=0;i<sharedRaw.size();++i)
        require(std::abs(sharedRaw[i]-independentRaw[i])<=1e-12*
            std::max(1.0,std::abs(independentRaw[i])),
            "Fused Hydrostatic raw tendency differs from the independent fallback");
    for(std::size_t family=0;family<coefficients.size();++family)
        for(std::size_t i=0;i<S;++i)
            require(coefficients[family][i].real==coefficientsBefore[family][i].real &&
                coefficients[family][i].imag==coefficientsBefore[family][i].imag,
                "Fused Hydrostatic nonlinear evaluation mutated input coefficients");
    for(std::size_t family=0;family<3;++family) for(std::size_t i=0;i<S;++i) {
        require(std::abs(sharedFluxStorage[family][i].real-
                    independentFluxStorage[family][i].real)<=1e-12*
                    std::max(1.0,std::abs(independentFluxStorage[family][i].real)) &&
                std::abs(sharedFluxStorage[family][i].imag-
                    independentFluxStorage[family][i].imag)<=1e-12*
                    std::max(1.0,std::abs(independentFluxStorage[family][i].imag)),
            "Shared nonlinear flux differs from the independent path");
    }
    require(shared->metrics().coefficientAssemblyCount==4 &&
            shared->metrics().verticalPreparationCount==7 &&
            shared->metrics().verticalOperatorExecutionCount==10 &&
            shared->metrics().horizontalSpectrumReuseCount==6 &&
            shared->metrics().preparedVerticalDerivativeCount==3 &&
            shared->metrics().derivativeAdvectionConsumerCount==9,
        "Scoped nonlinear evaluation did not share four assemblies and seven vertical products");
    require(independent->metrics().coefficientAssemblyCount==13 &&
            independent->metrics().verticalPreparationCount==13 &&
            independent->metrics().derivativeAdvectionConsumerCount==9,
        "Independent nonlinear qualification path did not execute all thirteen producers");
    require(bool(shared->endStateEvaluation()),
        "Shared nonlinear preparation scope end failed");
    require(bool(independent->endStateEvaluation()),
        "Independent nonlinear preparation scope end failed");

    shared->resetMetrics();
    require(bool(shared->nonlinearFlux(state,sharedFlux)),
        "Standalone shared nonlinear preparation failed");
    require(shared->metrics().coefficientAssemblyCount==4 &&
            shared->metrics().verticalPreparationCount==7 &&
            shared->metrics().verticalOperatorExecutionCount==10 &&
            shared->metrics().horizontalSpectrumReuseCount==6 &&
            shared->metrics().preparedVerticalDerivativeCount==3 &&
            shared->metrics().derivativeAdvectionConsumerCount==9,
        "Standalone nonlinear evaluation did not retain operation-local sharing");

    // The horizontal inverse consumes an already completed vertical product.
    // A consumer failure must retain that successfully published grid so a
    // retry does not repeat coefficient assembly or the vertical multiply.
    Counters failureCounters;
    std::unique_ptr<WVTransformHydrostaticKernel> failureKernel;
    require(bool(WVTransformHydrostaticKernel::create(source,
        std::make_unique<Engine>(failureCounters),failureKernel,
        WVCreateScalarMatrixBackend,sharedOptions)),
        "Hydrostatic consumer-failure kernel setup failed");
    require(bool(failureKernel->beginStateEvaluation(state)),
        "Hydrostatic consumer-failure scope begin failed");
    failureKernel->resetMetrics();
    failureCounters.failExecuteAt=failureCounters.executed;
    require(failureKernel->transformStateField(state,WVHydrostaticField::u,
                {sharedField.data(),spatialShape}).code==
            WVKernelStatusCode::fftExecutionFailure,
        "Injected Hydrostatic horizontal inverse failure was not preserved");
    failureCounters.failExecuteAt=-1;
    require(bool(failureKernel->transformStateField(state,WVHydrostaticField::u,
        {sharedField.data(),spatialShape})),
        "Hydrostatic retry after horizontal inverse failure failed");
    require(failureKernel->metrics().coefficientAssemblyCount==1 &&
            failureKernel->metrics().verticalPreparationCount==1 &&
            failureKernel->metrics().verticalOperatorExecutionCount==1 &&
            failureKernel->metrics().horizontalSpectrumReuseCount==1,
        "Hydrostatic retry repeated a successfully published vertical product");
    require(bool(failureKernel->endStateEvaluation()),
        "Hydrostatic consumer-failure scope end failed");
    for(auto& values:sharedFluxStorage) std::fill(values.begin(),values.end(),WVComplex64{17,19});
    failureKernel->resetMetrics();
    failureCounters.failExecuteAt=failureCounters.executed+4;
    require(failureKernel->nonlinearFlux(state,sharedFlux).code==WVKernelStatusCode::fftExecutionFailure &&
        failureKernel->metrics().derivativeAdvectionConsumerCount==0,
        "Failed Hydrostatic derivative inverse published a fused consumer");
    for(const auto& values:sharedFluxStorage) for(const auto value:values)
        require(value.real==17 && value.imag==19,
            "Failed Hydrostatic derivative inverse published spectral output");
}
#if WV_TEST_NATIVE_FFTW
void tiledNonlinearParity(const std::shared_ptr<const WVStratifiedModalRecord>& source) {
    const auto& g=source->geometry();const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    for (bool compact:{false,true}) for (std::size_t workers:{std::size_t{1},std::size_t{3},std::size_t{12}}) {
        WVVariableExecutionOptions options{WVRetainedHorizontalSchedule::streamingPrunedTile16,workers,true};
        options.pointwiseWorkers=2;options.fusedDerivativeAdvection=true;
        if (compact) options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
        const auto create=[&](bool tiled,std::unique_ptr<WVTransformHydrostaticKernel>& kernel) {
            std::unique_ptr<WVFFTEngine> engine;require(bool(WVFFTWEngine::create(1,engine)),"Native tiled engine setup failed");
            auto selected=options;selected.tiledNonlinear=tiled;
            require(bool(WVTransformHydrostaticKernel::create(source,std::move(engine),kernel,WVCreateScalarMatrixBackend,selected)),"Tiled family setup failed");
        };
        std::unique_ptr<WVTransformHydrostaticKernel> baseline,candidate;create(false,baseline);create(true,candidate);
        require(candidate->supportsTiledNonlinear() && !baseline->supportsTiledNonlinear(),"Tiled opt-in capability mismatch");
        auto referenceOptions=options;referenceOptions.tiledNonlinear=true;
        std::unique_ptr<WVTransformHydrostaticKernel> fallback;
        require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),fallback,WVCreateScalarMatrixBackend,referenceOptions)) && !fallback->supportsTiledNonlinear(),"Reference provider failed tiled fallback");
        std::array<std::vector<WVComplex64>,3> input,first,second;
        for (std::size_t j=0;j<3;++j) {
            input[j].resize(S);first[j].resize(S);second[j].resize(S);
            for (std::size_t i=0;i<S;++i) input[j][i]={.001*std::sin(.17*i+j),-.002*std::cos(.11*i+j)};
        }
        const auto shape=candidate->spectralShape();const auto volume=candidate->spatialShape();
        WVMutableCoefficients valid{{input[0].data(),shape},{input[1].data(),shape},{input[2].data(),shape}};
        require(bool(candidate->constrainCoefficients(valid)),"Tiled input constraints failed");
        const WVState state{83,17,{{input[0].data(),shape},{input[1].data(),shape},{input[2].data(),shape}}};
        WVFlux a{{first[0].data(),shape},{first[1].data(),shape},{first[2].data(),shape}};
        WVFlux b{{second[0].data(),shape},{second[1].data(),shape},{second[2].data(),shape}};
        std::vector<double> oldFields(4*R),newFields(4*R),oldRaw(3*R),newRaw(3*R),scratch(R);
        const WVRealFieldBundleView output{newFields.data(),{g.Nx,g.Ny,g.Nz,4}};
        WVRealFieldBundleView oldRawOutput{oldRaw.data(),{g.Nx,g.Ny,g.Nz,3}};
        WVRealFieldBundleView newRawOutput{newRaw.data(),{g.Nx,g.Ny,g.Nz,3}};
        const WVRealFieldBundleConstView borrowed{oldFields.data(),{g.Nx,g.Ny,g.Nz,4}};
        const WVHydrostaticField names[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::w,WVHydrostaticField::eta};
        for (std::size_t pass=0;pass<3;++pass) {
            // Same time and addresses, new evaluation with modified coefficients.
            if (pass) {input[0][g.Nj].real*=1.1;require(bool(candidate->constrainCoefficients(valid)),"Tiled changed-state constraints failed");}
            const auto snapshot=input;
            baseline->resetMetrics();candidate->resetMetrics();
            require(bool(baseline->beginStateEvaluation(state)) && bool(candidate->beginStateEvaluation(state)),"Tiled scope begin failed");
            for (std::size_t f=0;f<4;++f) require(bool(baseline->transformStateField(state,names[f],{oldFields.data()+f*R,volume})),"Legacy fields failed");
            require(bool(baseline->nonlinearFlux(state,a,&oldRawOutput,&borrowed)),"Legacy nonlinear failed");
            const auto bytes=candidate->persistentBytes();
            allocationProbe::calls=0;allocationProbe::counting=pass>0;
            const auto status=candidate->nonlinearFluxAndFields(state,b,output,&newRawOutput);
            allocationProbe::counting=false;
            require(bool(status),status.message.c_str());
            if (pass) require(allocationProbe::calls==0 && candidate->persistentBytes()==bytes,"Warmed tiled execution allocated");
            const auto& cm=candidate->metrics();const auto& bm=baseline->metrics();
            require(cm.tiledNonlinearCount==1 && cm.tiledColumnInverseCount==g.Nz*(4+2*3) &&
                cm.tiledRowInverseCount==g.Nz*(4+3*3) && cm.tiledReusedColumnCount==g.Nz*3,"Tiled producer counts mismatch");
            require(cm.coefficientAssemblyCount==bm.coefficientAssemblyCount && cm.verticalOperatorExecutionCount==bm.verticalOperatorExecutionCount && cm.reconstructionCount==bm.reconstructionCount,"Tiled producer or vertical work changed");
            for (std::size_t i=0;i<4*R;++i) require(oldFields[i]==newFields[i],"Tiled physical field changed");
            for (std::size_t i=0;i<3*R;++i)
                require(std::abs(oldRaw[i]-newRaw[i])<=1e-12*(1+std::abs(oldRaw[i])),
                    "Tiled raw physical tendency differs");
            for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i) {
                require(std::abs(first[j][i].real-second[j][i].real)<=1e-12*(1+std::abs(first[j][i].real)) &&
                    std::abs(first[j][i].imag-second[j][i].imag)<=1e-12*(1+std::abs(first[j][i].imag)),"Tiled projected flux differs");
                require(input[j][i].real==snapshot[j][i].real && input[j][i].imag==snapshot[j][i].imag,"Tiled input mutated");
            }
            require(bool(baseline->endStateEvaluation()) && bool(candidate->endStateEvaluation()),"Tiled scope end failed");
        }
        // No public flux API may overwrite another registered immutable view.
        int owner=0;
        require(bool(candidate->beginStateEvaluation(state,&owner)),"Tiled alias scope failed");
        const WVState alternate{state.t,state.t0,{{first[0].data(),shape},{first[1].data(),shape},{first[2].data(),shape}}};
        require(bool(candidate->addStateEvaluationView(alternate,&owner,1)),"Tiled alternate view registration failed");
        WVFlux forbidden{{first[0].data(),shape},b.Fm,b.F0};
        candidate->resetMetrics();
        require(candidate->nonlinearFluxAndFields(state,forbidden,output).code==WVKernelStatusCode::overlappingArrays &&
            candidate->nonlinearFlux(state,forbidden).code==WVKernelStatusCode::overlappingArrays &&
            candidate->metrics().tiledNonlinearCount==0,"Flux overwrote registered immutable view");
        require(bool(candidate->endStateEvaluation()),"Tiled alias scope end failed");
        // Invalid public output is rejected before FFT execution, then retry works.
        candidate->resetMetrics();auto bad=output;bad.shape.fourth=3;
        require(!candidate->nonlinearFluxAndFields(state,b,bad) && candidate->metrics().tiledNonlinearCount==0,"Invalid tiled output executed");
        auto badRaw=newRawOutput;badRaw.shape.fourth=4;
        require(!candidate->nonlinearFluxAndFields(state,b,output,&badRaw) && candidate->metrics().tiledNonlinearCount==0,
            "Invalid tiled tendency output executed");
        auto aliasedRaw=newRawOutput;aliasedRaw.data=newFields.data();
        require(candidate->nonlinearFluxAndFields(state,b,output,&aliasedRaw).code==WVKernelStatusCode::overlappingArrays &&
            candidate->metrics().tiledNonlinearCount==0,"Aliased tiled tendency output executed");
        require(bool(candidate->nonlinearFluxAndFields(state,b,output)),"Tiled failure recovery failed");
        // Calling the established API remains valid even when tiling is prepared.
        candidate->resetMetrics();require(bool(candidate->nonlinearFlux(state,b,nullptr,&borrowed)) && candidate->metrics().tiledNonlinearCount==0,"Borrowed fields unexpectedly tiled");
        require(baseline->nonlinearFluxAndFields(state,a,output).code==WVKernelStatusCode::unsupportedOperation,"Disabled tiled API did not reject");
    }
}
#endif

}
int main() {
    try {
        Temporary file; fixture(file.path);
        { File f(file.path); for (const auto* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformHydrostatic"),"WVTransformHydrostatic")); }
        std::shared_ptr<const WVStratifiedModalRecord> source; auto status=WVStratifiedModalReader::read(file.path.string(),source); require(bool(status),status.message.c_str());
        contracts(source);
#if WV_TEST_NATIVE_FFTW
        { Temporary tiledFile; fixture(tiledFile.path);
          {File f(tiledFile.path);for(const auto* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformHydrostatic"),"WVTransformHydrostatic"));}
          for(std::size_t z=0;z<7;++z) change(tiledFile.path,"dLnN2",.001*(z+1),z);
          std::shared_ptr<const WVStratifiedModalRecord> tiledSource;
          auto readStatus=WVStratifiedModalReader::read(tiledFile.path.string(),tiledSource);require(bool(readStatus),readStatus.message.c_str());
          tiledNonlinearParity(tiledSource); }
#endif
        variableScheduleParity(source);
        variableScheduleParity(source,true);
        sharedFieldGradientParity(source,false);
        sharedFieldGradientParity(source,true);
        std::unique_ptr<WVTransformHydrostaticKernel> kernel; require(bool(WVTransformHydrostaticKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),kernel)),"Lifetime setup failed");
        std::weak_ptr<const WVStratifiedModalRecord> weak=source; source.reset(); require(!weak.expired(),"Kernel lost scientific owner"); kernel.reset(); require(weak.expired(),"Scientific owner leaked");
        std::cout<<"Hydrostatic kernel contracts passed\n"; return 0;
    } catch(const std::exception& e) { allocationProbe::failAfter=-1; allocationProbe::counting=false; std::cerr<<e.what()<<'\n'; return 1; }
}
