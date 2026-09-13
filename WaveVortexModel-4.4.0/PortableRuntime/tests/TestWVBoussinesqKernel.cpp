#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WVReferenceFFTEngine.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#include "WVBoussinesqModalTestFixture.hpp"
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
    std::unique_ptr<WVTransformBoussinesqKernel> kernel;
    require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<Engine>(counters),kernel)),"Create failed");
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
    require(bool(kernel->transformStateField(state,WVBoussinesqField::u,field)),"Scoped field reconstruction failed");
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
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::u)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1,
        "Scoped calls repeated preparation or lost reconstruction metrics");
    require(!kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},state.t+1,state.t0,out),
        "Mismatched scoped projection time was accepted");
    auto registeredCoefficients=a;
    auto foreignState=state;
    foreignState.coefficients={{registeredCoefficients[0].data(),shape},{registeredCoefficients[1].data(),shape},{registeredCoefficients[2].data(),shape}};
    require(!kernel->transformStateField(foreignState,WVBoussinesqField::u,field),"Foreign scoped state was accepted");
    std::vector<double> tendencyFields(4*R);
    WVRealFieldBundleView tendencyBundle{tendencyFields.data(),{g.Nx,g.Ny,g.Nz,4}};
    require(bool(kernel->transformCoefficientTendencyToUVWEta(foreignState,tendencyBundle)),
        "Derived coefficient tendency was rejected by the active primary-state scope");
    require(kernel->metrics().stateValidationCount==1 &&
        kernel->metrics().derivedValidationCount==1 &&
        kernel->metrics().tendencyReconstructionCount==std::array<std::size_t,4>{1,1,1,1} &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::u)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1,
        "Derived tendency validation or production was attributed to the primary state");
    require(!kernel->addStateEvaluationView(foreignState,&foreignOwner),"Foreign evaluation owner registered a state view");
    require(bool(kernel->addStateEvaluationView(foreignState,&evaluationOwner,2)),"Additional immutable state view was rejected");
    const auto registeredStatus=kernel->transformStateField(foreignState,WVBoussinesqField::u,field);
    require(bool(registeredStatus),registeredStatus.message.c_str());
    require(kernel->metrics().stateValidationCount==2 && kernel->metrics().phasePreparationCount==1,
        "Additional state view skipped validation or repeated phase preparation");
    require(kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::u)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)][2]==1,
        "Registered Boussinesq component production lost its component identity");
    require(!kernel->removeStateEvaluationView(state,&evaluationOwner,0),
        "Primary Boussinesq state view was removed");
    require(!kernel->removeStateEvaluationView(foreignState,&foreignOwner,2),
        "Foreign owner removed a Boussinesq state view");
    require(bool(kernel->removeStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Boussinesq state view removal failed");
    require(!kernel->validateStateEvaluation(foreignState) &&
                !kernel->removeStateEvaluationView(foreignState,&evaluationOwner,2),
        "Removed Boussinesq state view remained registered");
    require(bool(kernel->addStateEvaluationView(foreignState,&evaluationOwner,3)) &&
                bool(kernel->validateStateEvaluation(foreignState)),
        "Boussinesq state view storage could not be re-registered");
    require(bool(kernel->endStateEvaluation()),"End scoped state evaluation failed");
    const auto savedScopedCoefficient=a[0][0]; a[0][0].real=std::numeric_limits<double>::infinity();
    require(!kernel->beginStateEvaluation(state),"New scope skipped validation for mutated same-pointer state");
    a[0][0]=savedScopedCoefficient;
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVBoussinesqField::u,field)),"Standalone field reconstruction failed");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::v,field)),"Second standalone field reconstruction failed");
    require(kernel->metrics().stateValidationCount==2 && kernel->metrics().phasePreparationCount==2,
        "Standalone calls did not perform fresh state preparation");
    std::vector<double> compoundReference(R),compoundFirst(R),compoundSecond(R),compoundResult(R);
    require(bool(kernel->beginStateEvaluation(state)),"Compound producer scope begin failed");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::zetaX,
        {compoundReference.data(),volume})),"Reference horizontal vorticity failed");
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVBoussinesqField::w,
        {compoundFirst.data(),volume},WVBoussinesqDerivative::y)),"Prepared w_y failed");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::v,
        {compoundSecond.data(),volume},WVBoussinesqDerivative::z)),"Prepared v_z failed");
    require(bool(kernel->combinePreparedHorizontalVorticity(WVBoussinesqField::zetaX,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume})),"Prepared horizontal vorticity combination failed");
    require(compoundResult==compoundReference &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::w)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::y)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::v)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::z)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::zetaX)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1,
        "Prepared horizontal vorticity changed values or repeated a producer");
    require(bool(kernel->combinePreparedHorizontalVorticity(WVBoussinesqField::zetaX,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundFirst.data(),volume},WVBoussinesqComponent::geostrophic)) &&
        compoundFirst==compoundReference &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::zetaX)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
            [static_cast<std::size_t>(WVBoussinesqComponent::geostrophic)]==1,
        "Exact prepared vorticity alias or component attribution failed");
    std::vector<double> partialCompound(R+1);
    std::copy_n(compoundReference.data(),R,partialCompound.data());
    require(kernel->combinePreparedHorizontalVorticity(WVBoussinesqField::zetaX,
        {partialCompound.data(),volume},{compoundSecond.data(),volume},
        {partialCompound.data()+1,volume}).code==WVKernelStatusCode::overlappingArrays,
        "Partial prepared vorticity alias was accepted");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::rhoTotal,
        {compoundReference.data(),volume},WVBoussinesqDerivative::z)),"Reference density derivative failed");
    kernel->resetMetrics();
    require(bool(kernel->transformStateField(state,WVBoussinesqField::eta,
        {compoundFirst.data(),volume},WVBoussinesqDerivative::z)),"Prepared eta_z failed");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::eta,
        {compoundSecond.data(),volume})),"Prepared eta failed");
    require(bool(kernel->combinePreparedDensityZDerivative(WVBoussinesqField::rhoTotal,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume})),"Prepared density derivative combination failed");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundResult[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Prepared density derivative changed values");
    require(kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::eta)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::z)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::eta)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVBoussinesqField::rhoTotal)]
            [static_cast<std::size_t>(WVBoussinesqDerivative::z)]
            [static_cast<std::size_t>(WVBoussinesqComponent::all)]==1,
        "Prepared density derivative repeated a producer");
    require(bool(kernel->combinePreparedDensityZDerivative(WVBoussinesqField::rhoTotal,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundFirst.data(),volume})),
        "Exact prepared density alias was rejected");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundFirst[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Exact prepared density alias changed values");
    std::copy_n(compoundReference.data(),R,partialCompound.data());
    require(kernel->combinePreparedDensityZDerivative(WVBoussinesqField::rhoTotal,
        {partialCompound.data(),volume},{compoundSecond.data(),volume},
        {partialCompound.data()+1,volume}).code==WVKernelStatusCode::overlappingArrays,
        "Partial prepared density alias was accepted");
    require(!kernel->combinePreparedDensityZDerivative(WVBoussinesqField::u,
        {compoundFirst.data(),volume},{compoundSecond.data(),volume},
        {compoundResult.data(),volume}),"Invalid prepared density target was accepted");
    require(bool(kernel->endStateEvaluation()),"Compound producer scope end failed");
    std::fill(spatial.begin(),spatial.end(),29);
    allocationProbe::calls=0; allocationProbe::counting=true;
    require(bool(kernel->beginStateEvaluation(state)),"Allocation probe scope begin failed");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::u,field)),"Allocation probe scoped field failed");
    require(bool(kernel->endStateEvaluation()),"Allocation probe scope end failed");
    allocationProbe::counting=false;
    require(allocationProbe::calls==0,"Scoped state preparation or reuse allocated");
    std::fill(spatial.begin(),spatial.end(),29);
    auto status=kernel->nonlinearFlux(state,flux); require(bool(status),status.message.c_str());
    const auto referenceFlux=b,referenceState=a;
    const auto retainedBytes=kernel->persistentBytes();
    std::vector<double> raw(4*R),physical(4*R),gradient(R),expected(4*R,0),
        absoluteContributions(4*R,0);
    WVRealFieldBundleView rawView{raw.data(),{g.Nx,g.Ny,g.Nz,4}};
    const WVBoussinesqField dynamical[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    for(std::size_t channel=0;channel<4;++channel)
        require(bool(kernel->transformStateField(state,dynamical[channel],{physical.data()+channel*R,volume})),"Prepare shared physical fields");
    const WVRealFieldBundleConstView prepared{physical.data(),{g.Nx,g.Ny,g.Nz,4}};
    const auto physicalBefore=physical;
    const auto fullStart=counters.executed;
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared)),"Observe raw nonlinear tendency");
    const auto fullExecutions=counters.executed-fullStart;
    for(std::size_t channel=0;channel<4;++channel) {
        const auto target=channel;
        for(std::size_t axis=0;axis<3;++axis) {
            require(bool(kernel->transformStateField(state,dynamical[target],{gradient.data(),volume},static_cast<WVBoussinesqDerivative>(axis+1))),"Independent tendency derivative");
            for(std::size_t i=0;i<R;++i) {
                const double correction=target==3 && axis==2 ? physical[3*R+i]*g.dLnN2[i/(g.Nx*g.Ny)] : 0;
                const double contribution=physical[axis*R+i]*(gradient[i]+correction);
                expected[channel*R+i]-=contribution;
                absoluteContributions[channel*R+i]+=std::abs(contribution);
            }
        }
    }
    const double roundoffFactor=8*std::numeric_limits<double>::epsilon();
    const auto rawMatchesExpected=[&]() {
        for(std::size_t i=0;i<raw.size();++i)
            if(!(std::abs(raw[i]-expected[i])<=
                    roundoffFactor*absoluteContributions[i])) return false;
        return true;
    };
    const bool rawWithinRoundoff=rawMatchesExpected();
    if(!rawWithinRoundoff || physical!=physicalBefore) {
        std::size_t firstRaw=raw.size(),firstPhysical=physical.size();
        double firstRawBound=0,maximumRawError=0,maximumPhysicalError=0;
        for(std::size_t i=0;i<raw.size();++i) {
            const auto error=std::abs(raw[i]-expected[i]);
            const auto bound=roundoffFactor*absoluteContributions[i];
            if(!(error<=bound) && firstRaw==raw.size()) {
                firstRaw=i;
                firstRawBound=bound;
            }
            if(!std::isfinite(error) || error>maximumRawError) maximumRawError=error;
        }
        for(std::size_t i=0;i<physical.size();++i) {
            const auto error=std::abs(physical[i]-physicalBefore[i]);
            if(physical[i]!=physicalBefore[i] && firstPhysical==physical.size()) firstPhysical=i;
            if(!std::isfinite(error) || error>maximumPhysicalError) maximumPhysicalError=error;
        }
        const auto precision=std::cerr.precision(17);
        std::cerr<<"Boussinesq raw mismatch: firstRaw="<<firstRaw;
        if(firstRaw<raw.size())
            std::cerr<<" actual="<<raw[firstRaw]<<" expected="<<expected[firstRaw]
                <<" bound="<<firstRawBound;
        std::cerr<<" maxRawAbs="<<maximumRawError
            <<" borrowedEqual="<<(physical==physicalBefore)
            <<" firstBorrowed="<<firstPhysical;
        if(firstPhysical<physical.size())
            std::cerr<<" actualBorrowed="<<physical[firstPhysical]
                <<" expectedBorrowed="<<physicalBefore[firstPhysical];
        std::cerr<<" maxBorrowedAbs="<<maximumPhysicalError<<'\n';
        std::cerr.precision(precision);
    }
    require(rawWithinRoundoff && physical==physicalBefore,
        "Raw tendencies precede projection and preserve borrowed fields");
    require(bool(kernel->transformStateField(state,WVBoussinesqField::u,
        {gradient.data(),volume},WVBoussinesqDerivative::z)),
        "Prepare cached Boussinesq derivative");
    kernel->resetMetrics();
    DerivativeAccessProbe derivativeProbe;
    derivativeProbe.cachedField=static_cast<std::size_t>(WVBoussinesqField::u);
    derivativeProbe.cachedDerivative=static_cast<std::size_t>(WVBoussinesqDerivative::z);
    derivativeProbe.cached={gradient.data(),volume};
    auto derivativeAccess=derivativeProbe.access();
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess)) && rawMatchesExpected(),
        "Boussinesq nonlinear derivative reuse changed the spatial tendency");
    require(kernel->metrics().reconstructionCount
                [static_cast<std::size_t>(WVBoussinesqField::u)]
                [static_cast<std::size_t>(WVBoussinesqDerivative::z)]
                [static_cast<std::size_t>(WVBoussinesqComponent::all)]==0 &&
            derivativeProbe.captures==11,
        "Boussinesq nonlinear path reproduced a cached derivative or missed capture");
    derivativeProbe.failLookup=true;
    require(!kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess),
        "Boussinesq nonlinear derivative lookup failure was ignored");
    derivativeProbe.failLookup=false;
    derivativeProbe.failCapture=true;
    derivativeProbe.cachedField=99;
    require(!kernel->nonlinearFlux(state,flux,&rawView,&prepared,false,
                &derivativeAccess),
        "Boussinesq nonlinear derivative capture failure was ignored");
    for(std::size_t channel=0;channel<3;++channel) for(std::size_t i=0;i<S;++i)
        require(b[channel][i].real==referenceFlux[channel][i].real && b[channel][i].imag==referenceFlux[channel][i].imag &&
                a[channel][i].real==referenceState[channel][i].real && a[channel][i].imag==referenceState[channel][i].imag,"Observation changed flux or input state");
    for (auto& values:b) std::fill(values.begin(),values.end(),WVComplex64{17,19});
    const auto rawStart=counters.executed;
    require(bool(kernel->nonlinearFlux(state,flux,&rawView,&prepared,false)),"Spatial-only nonlinear tendency");
    require(rawMatchesExpected() && physical==physicalBefore && counters.executed-rawStart<fullExecutions,
        "Spatial-only evaluation changed the tendency or retained redundant projection work");
    const auto rawBeforeInvalid=raw;
    for(const auto& values:b) for(const auto value:values)
        require(value.real==17 && value.imag==19,"Spatial-only evaluation wrote spectral flux");
    const auto rejectedStart=counters.executed;
    require(kernel->nonlinearFlux(state,flux,nullptr,&prepared,false).code==WVKernelStatusCode::invalidConfiguration &&
        counters.executed==rejectedStart,"Spatial-only evaluation accepted missing output or rejected after FFT work");
    auto badRaw=rawView; badRaw.data=physical.data();
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::overlappingArrays,"Observed tendency aliases shared fields");
    badRaw=rawView; badRaw.shape.fourth=5;
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::invalidShape,"Wrong tendency channel count accepted");
    badRaw=rawView; badRaw.data=nullptr;
    require(kernel->nonlinearFlux(state,flux,&badRaw,&prepared).code==WVKernelStatusCode::invalidPointer,"Null tendency output accepted");
    require(raw==rawBeforeInvalid && physical==physicalBefore,"Invalid observation mutated output");
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
    require(!kernel->transformStateField(state,static_cast<WVBoussinesqField>(99),field),"Invalid field accepted");
    require(!kernel->transformStateField(state,WVBoussinesqField::u,field,static_cast<WVBoussinesqDerivative>(99)),"Invalid derivative accepted");
    require(!kernel->transformStateField(state,WVBoussinesqField::u,field,WVBoussinesqDerivative::value,static_cast<WVBoussinesqComponent>(99)),"Invalid component accepted");
    require(!kernel->transformStateField(state,WVBoussinesqField::zetaX,field,WVBoussinesqDerivative::x),"Unsupported vorticity derivative accepted");
    require(!kernel->transformStateField(state,WVBoussinesqField::rhoTotal,field,WVBoussinesqDerivative::value,WVBoussinesqComponent::wave),"Masked background accepted");
    require(!kernel->transformStateField(state,WVBoussinesqField::ssh,field),"Wrong surface shape accepted");
    require(!kernel->differentiateVertical({scratch.data(),volume},WVBoussinesqFamily::F,0,field),"Order zero accepted");
    require(!kernel->differentiateVertical({scratch.data(),volume},WVBoussinesqFamily::F,5,field),"Order five accepted");
    require(!kernel->integrateVertical({scratch.data(),volume},static_cast<WVBoussinesqFamily>(99),field),"Invalid family accepted");
    require(!kernel->differentiateVertical({spatial.data()+1,volume},WVBoussinesqFamily::G,1,field),"Partial calculus alias accepted");
    auto aliased=out; aliased.Ap=amplitudes.Am;
    require(!kernel->evolveCoefficients(state,aliased),"Cross-family evolution alias accepted");
    auto invalidVolume=volume; invalidVolume.third=1;
    require(!kernel->transformUVWEtaToWaveVortex({scratch.data(),volume},{scratch.data(),volume},{scratch.data(),invalidVolume},{scratch.data(),volume},83,17,out),"Invalid vertical-velocity shape accepted");
    require(!kernel->transformUVWEtaToWaveVortex({scratch.data(),volume},{scratch.data(),volume},{reinterpret_cast<const double*>(out.Ap.data),volume},{scratch.data(),volume},83,17,out),"Aliased vertical velocity accepted");
    const auto saved=a;
    a[2][S-1].real=std::numeric_limits<double>::infinity();
    require(!kernel->transformStateField(state,WVBoussinesqField::u,field),"Infinite coefficient accepted");
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
        require(bool(kernel->transformStateField(state,WVBoussinesqField::rhoTotal,field,WVBoussinesqDerivative::z)),"Prepared density failed");
        require(bool(kernel->transformUVEtaToWaveVortex({spatial.data(),volume},{spatial.data(),volume},{spatial.data(),volume},83,17,out)),"Prepared projection failed");
        require(bool(kernel->transformUVWEtaToWaveVortex({scratch.data(),volume},{scratch.data(),volume},{scratch.data(),volume},{scratch.data(),volume},83,17,out)),"Prepared four-field projection failed");
        require(bool(kernel->differentiateVertical({scratch.data(),volume},WVBoussinesqFamily::F,4,field)),"Prepared derivative failed");
        require(bool(kernel->integrateVertical({scratch.data(),volume},WVBoussinesqFamily::G,field)),"Prepared integral failed");
    }
    allocationProbe::counting=false; require(allocationProbe::calls==0 && kernel->persistentBytes()==bytes,"Prepared execution allocates or changes storage");
    const auto* old=kernel.get();
    status=WVTransformBoussinesqKernel::create({},std::make_unique<WVReferenceFFTEngine>(),kernel);
    require(!status && kernel.get()==old,"Failed setup replaced kernel");
    for (int fail=0;fail<11;++fail) {
        int count=0; Counters counters;
        status=WVTransformBoussinesqKernel::create(source,std::make_unique<Engine>(counters),kernel,[&](std::unique_ptr<WVVerticalMatrixBackend>& backend) {
            if (count++==fail) return WVKernelStatus{WVKernelStatusCode::allocationFailure,"Injected matrix factory failure."};
            return WVCreateScalarMatrixBackend(backend);
        });
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"Matrix setup failure leaked or replaced kernel");
    }
    for (int fail=0;fail<2;++fail) {
        Counters counters; counters.failAt=fail;
        status=WVTransformBoussinesqKernel::create(source,std::make_unique<Engine>(counters),kernel);
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"FFT setup failure leaked or replaced kernel");
    }
    bool succeeded=false;
    for (long fail=0;fail<2048;++fail) {
        Counters counters; auto engine=std::make_unique<Engine>(counters);
        allocationProbe::failAfter=fail;
        status=WVTransformBoussinesqKernel::create(source,std::move(engine),kernel);
        allocationProbe::failAfter=-1;
        if (status) { kernel.reset(); require(!counters.engines && !counters.plans,"Successful destruction leaked"); succeeded=true; break; }
        require(kernel.get()==old && !counters.engines && !counters.plans,"Allocation failure leaked or replaced kernel");
    }
    require(succeeded,"Allocation sweep never succeeded");
}
void variableScheduleParity(const std::shared_ptr<const WVStratifiedModalRecord>& source,bool compact = false) {
    Counters candidateCounters; // Outlives the instrumented plans owned by candidate.
    std::unique_ptr<WVTransformBoussinesqKernel> frozen, candidate;
    WVVariableExecutionOptions frozenOptions;
    frozenOptions.inertialOnlyProjection=false;
    require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),frozen,
        WVCreateScalarMatrixBackend,frozenOptions)),"Frozen schedule setup failed");
    WVVariableExecutionOptions options{WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true};
    if (compact) options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
    options.pointwiseWorkers=2;
    options.fusedDerivativeAdvection=true;
    options.verticalGroupWorkers=2;
    require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<Engine>(candidateCounters),candidate, WVCreateScalarMatrixBackend, options)),"Candidate schedule setup failed");
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
    frozen->resetMetrics(); candidate->resetMetrics();
    require(bool(frozen->nonlinearFlux(state,frozenFlux)),"Frozen nonlinear flux failed");
    require(bool(candidate->nonlinearFlux(state,candidateFlux)),"Candidate nonlinear flux failed");
    require(!frozen->executionOptions().fusedDerivativeAdvection &&
        candidate->executionOptions().fusedDerivativeAdvection &&
        candidate->metrics().derivativeAdvectionConsumerCount==12,
        "Boussinesq fused derivative consumers were not explicit or complete");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i) require(std::abs(frozenOut[j][i].real-candidateOut[j][i].real)<1e-12 && std::abs(frozenOut[j][i].imag-candidateOut[j][i].imag)<1e-12,"Candidate nonlinear flux differs");
    require(frozen->metrics().verticalOperatorExecutionCount==candidate->metrics().verticalOperatorExecutionCount &&
        frozen->metrics().verticalMatrixGroupExecutionCount==
            candidate->metrics().verticalMatrixGroupExecutionCount+2*(source->groups().size()-1),
        "Inertial-only projection did not replace two full grouped wave projections with two exact columns");
    auto serialOptions=options; serialOptions.pointwiseWorkers=1;
    std::unique_ptr<WVTransformBoussinesqKernel> serial; std::array<std::vector<WVComplex64>,3> serialOut;
    for (auto& values:serialOut) values.resize(S);
    require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),serial,
        WVCreateScalarMatrixBackend,serialOptions)),"Single-worker candidate setup failed");
    WVFlux serialFlux{{serialOut[0].data(),{g.Nj,g.Nkl}},{serialOut[1].data(),{g.Nj,g.Nkl}},{serialOut[2].data(),{g.Nj,g.Nkl}}};
    require(bool(serial->nonlinearFlux(state,serialFlux)),"Single-worker candidate flux failed");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i)
        require(serialOut[j][i].real==candidateOut[j][i].real && serialOut[j][i].imag==candidateOut[j][i].imag,
            "Pointwise worker partition changed boussinesq arithmetic");
    std::vector<double> ff(4*R),cf(4*R),fr(4*R),cr(4*R);
    const WVBoussinesqField fields[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    for (std::size_t j=0;j<4;++j) { require(bool(frozen->transformStateField(state,fields[j],{ff.data()+j*R,{g.Nx,g.Ny,g.Nz}})),"Frozen borrowed fields failed"); require(bool(candidate->transformStateField(state,fields[j],{cf.data()+j*R,{g.Nx,g.Ny,g.Nz}})),"Candidate borrowed fields failed"); }
    const auto ffBefore=ff,cfBefore=cf;
    bool nonzeroDensityCorrection=false;
    for(std::size_t i=0;i<R;++i)
        nonzeroDensityCorrection|=cf[2*R+i]*cf[3*R+i]*g.dLnN2[i/(g.Nx*g.Ny)]!=0;
    require(nonzeroDensityCorrection,"Boussinesq fused parity did not exercise the density-gradient correction");
    WVRealFieldBundleConstView ffields{ff.data(),{g.Nx,g.Ny,g.Nz,4}},cfields{cf.data(),{g.Nx,g.Ny,g.Nz,4}}; WVRealFieldBundleView frv{fr.data(),{g.Nx,g.Ny,g.Nz,4}},crv{cr.data(),{g.Nx,g.Ny,g.Nz,4}};
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
    require(bool(candidate->transformStateField(state,WVBoussinesqField::u,
        {cachedDerivative.data(),{g.Nx,g.Ny,g.Nz}},WVBoussinesqDerivative::z)) &&
        bool(candidate->transformStateField(state,WVBoussinesqField::u,
        {captureReference.data(),{g.Nx,g.Ny,g.Nz}},WVBoussinesqDerivative::x)),
        "Boussinesq fused derivative cache fixtures failed");
    DerivativeAccessProbe derivativeProbe;
    derivativeProbe.cachedField=static_cast<std::size_t>(WVBoussinesqField::u);
    derivativeProbe.cachedDerivative=static_cast<std::size_t>(WVBoussinesqDerivative::z);
    derivativeProbe.cached={cachedDerivative.data(),{g.Nx,g.Ny,g.Nz}};
    derivativeProbe.captureField=static_cast<std::size_t>(WVBoussinesqField::u);
    derivativeProbe.captureDerivative=static_cast<std::size_t>(WVBoussinesqDerivative::x);
    derivativeProbe.captured=&capturedDerivative;
    auto derivativeAccess=derivativeProbe.access();
    candidate->resetMetrics();
    const auto cachedStart=candidateCounters.executed;
    require(bool(candidate->nonlinearFlux(state,candidateFlux,&crv,&cfields,false,&derivativeAccess)),
        "Boussinesq fused cached-derivative evaluation failed");
    require(candidateCounters.executed-cachedStart+1==uncachedExecutions &&
        candidate->metrics().derivativeAdvectionConsumerCount==11 &&
        derivativeProbe.captures==11 && capturedDerivative==captureReference &&
        cf==cfBefore,
        "Boussinesq cached derivative repeated an inverse/consumer or lost capture");
    for (std::size_t j=0;j<3;++j) for (std::size_t i=0;i<S;++i)
        require(input[j][i].real==inputBefore[j][i].real && input[j][i].imag==inputBefore[j][i].imag,"Parity evaluation mutated input coefficients");
    require(frozen->storage().realScratchBytes==11*R*sizeof(double) && candidate->storage().realScratchBytes==6*R*sizeof(double),
        "Boussinesq streamed scratch accounting differs");
    require(candidate->executionOptions().horizontalWorkers==2 && candidate->executionOptions().streamedNonlinear &&
        candidate->executionOptions().usesCompactSplitViews()==compact && candidate->executionOptions().pointwiseWorkers==2 &&
        candidate->executionOptions().fusedDerivativeAdvection &&
        candidate->executionOptions().inertialOnlyProjection && !frozen->executionOptions().inertialOnlyProjection &&
        candidate->executionOptions().verticalGroupWorkers==2 && frozen->executionOptions().verticalGroupWorkers==1,
        "Candidate options were not retained");
    require(candidate->persistentBytes()>=candidate->storage().workspaceBytes,"Candidate storage ledger under-reports workspace");

    auto* retained=candidate.get(); WVVariableExecutionOptions invalidWorkerOptions; invalidWorkerOptions.pointwiseWorkers=0;
    auto status=WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),candidate,
        WVCreateScalarMatrixBackend,invalidWorkerOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && candidate.get()==retained,
        "Zero pointwise workers were accepted or replaced the retained kernel");
    WVVariableExecutionOptions invalidVerticalWorkerOptions;
    invalidVerticalWorkerOptions.verticalGroupWorkers=0;
    status=WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),candidate,
        WVCreateScalarMatrixBackend,invalidVerticalWorkerOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && candidate.get()==retained,
        "Zero vertical group workers were accepted or replaced the retained kernel");
    std::unique_ptr<WVTransformBoussinesqKernel> invalid;
    WVVariableExecutionOptions invalidOptions;
    invalidOptions.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
    status=WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),invalid,
        WVCreateScalarMatrixBackend,invalidOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && !invalid,
        "Compact split views accepted a non-streaming horizontal schedule");
    invalidOptions.spectralSchedule=static_cast<WVVariableSpectralSchedule>(99);
    status=WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),invalid,
        WVCreateScalarMatrixBackend,invalidOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && !invalid,
        "Unknown variable spectral schedule was accepted");
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
    std::unique_ptr<WVTransformBoussinesqKernel> shared,independent;
    require(bool(WVTransformBoussinesqKernel::create(source,
        std::make_unique<WVReferenceFFTEngine>(),shared,WVCreateScalarMatrixBackend,
        sharedOptions)),"Shared Boussinesq field-gradient kernel setup failed");
    require(bool(WVTransformBoussinesqKernel::create(source,
        std::make_unique<WVReferenceFFTEngine>(),independent,WVCreateScalarMatrixBackend,
        independentOptions)),"Independent Boussinesq field-gradient kernel setup failed");

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
        "Shared Boussinesq field-gradient input constraints failed");
    const auto coefficientsBefore=coefficients;
    const WVState state{83,17,{{coefficients[0].data(),spectralShape},
        {coefficients[1].data(),spectralShape},{coefficients[2].data(),spectralShape}}};
    std::vector<double> sharedField(R),independentField(R);
    require(bool(shared->beginStateEvaluation(state)),
        "Shared Boussinesq field-gradient scope begin failed");
    require(bool(independent->beginStateEvaluation(state)),
        "Independent Boussinesq field-gradient scope begin failed");
    shared->resetMetrics();
    independent->resetMetrics();
    struct FieldComponent { WVBoussinesqField field; WVBoussinesqComponent component; };
    const FieldComponent keys[]={
        {WVBoussinesqField::u,WVBoussinesqComponent::all},
        {WVBoussinesqField::u,WVBoussinesqComponent::wave},
        {WVBoussinesqField::v,WVBoussinesqComponent::all},
        {WVBoussinesqField::v,WVBoussinesqComponent::inertial},
        {WVBoussinesqField::w,WVBoussinesqComponent::all},
        {WVBoussinesqField::w,WVBoussinesqComponent::geostrophic},
        {WVBoussinesqField::eta,WVBoussinesqComponent::all},
        {WVBoussinesqField::eta,WVBoussinesqComponent::meanDensityAnomaly}};
    for(const auto key:keys) for(std::size_t derivative=0;derivative<4;++derivative) {
        const auto d=static_cast<WVBoussinesqDerivative>(derivative);
        require(bool(shared->transformStateField(state,key.field,
            {sharedField.data(),spatialShape},d,key.component)),
            "Shared Boussinesq field-gradient reconstruction failed");
        require(bool(independent->transformStateField(state,key.field,
            {independentField.data(),spatialShape},d,key.component)),
            "Independent Boussinesq field-gradient reconstruction failed");
        for(std::size_t i=0;i<R;++i)
            require(std::abs(sharedField[i]-independentField[i])<=1e-12*
                std::max(1.0,std::abs(independentField[i])),
                "Shared Boussinesq reconstruction differs from the independent path");
    }
    const auto sharedFieldMetrics=shared->metrics();
    const auto independentFieldMetrics=independent->metrics();
    require(sharedFieldMetrics.coefficientAssemblyCount==8 &&
            sharedFieldMetrics.verticalPreparationCount==16 &&
            sharedFieldMetrics.verticalOperatorExecutionCount==32 &&
            sharedFieldMetrics.horizontalSpectrumReuseCount==24 &&
            sharedFieldMetrics.preparedVerticalDerivativeCount==8,
        "Shared Boussinesq producers did not execute once per field/component key");
    require(independentFieldMetrics.coefficientAssemblyCount==32 &&
            independentFieldMetrics.verticalPreparationCount==64 &&
            independentFieldMetrics.horizontalSpectrumReuseCount==0 &&
            independentFieldMetrics.preparedVerticalDerivativeCount==0,
        "Independent Boussinesq qualification path unexpectedly reused preparation");
    require(bool(shared->endStateEvaluation()),
        "Shared Boussinesq field-gradient scope end failed");
    require(bool(independent->endStateEvaluation()),
        "Independent Boussinesq field-gradient scope end failed");

    auto foreignCoefficients=coefficients;
    WVState foreignState=state;
    foreignState.coefficients={{foreignCoefficients[0].data(),spectralShape},
        {foreignCoefficients[1].data(),spectralShape},
        {foreignCoefficients[2].data(),spectralShape}};
    int evaluationOwner=0;
    require(bool(shared->beginStateEvaluation(state,&evaluationOwner)),
        "Registered Boussinesq view scope begin failed");
    require(bool(shared->addStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Registered Boussinesq view setup failed");
    shared->resetMetrics();
    require(bool(shared->transformStateField(foreignState,WVBoussinesqField::u,
        {sharedField.data(),spatialShape})),"First registered Boussinesq view failed");
    require(bool(shared->transformStateField(foreignState,WVBoussinesqField::u,
        {sharedField.data(),spatialShape},WVBoussinesqDerivative::x)),
        "Registered Boussinesq view horizontal reuse failed");
    require(bool(shared->removeStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Prepared Boussinesq view removal failed");
    require(bool(shared->addStateEvaluationView(foreignState,&evaluationOwner,2)),
        "Same-storage Boussinesq view re-add failed");
    require(bool(shared->transformStateField(foreignState,WVBoussinesqField::u,
        {sharedField.data(),spatialShape})),"Re-added Boussinesq view failed");
    require(shared->metrics().coefficientAssemblyCount==2 &&
            shared->metrics().verticalPreparationCount==4 &&
            shared->metrics().horizontalSpectrumReuseCount==1,
        "Re-added Boussinesq view reused an invalid generation or lost in-view reuse");
    require(bool(shared->endStateEvaluation()),"Registered Boussinesq view scope end failed");

    shared->resetMetrics();
    require(bool(shared->transformStateField(state,WVBoussinesqField::u,
        {sharedField.data(),spatialShape})),"First standalone Boussinesq preparation failed");
    require(bool(shared->transformStateField(state,WVBoussinesqField::u,
        {sharedField.data(),spatialShape})),"Second standalone Boussinesq preparation failed");
    require(shared->metrics().coefficientAssemblyCount==2 &&
            shared->metrics().verticalPreparationCount==4,
        "Standalone Boussinesq operations published preparation across calls");

    std::array<std::vector<WVComplex64>,3> sharedFluxStorage,independentFluxStorage;
    for(auto& values:sharedFluxStorage) values.resize(S);
    for(auto& values:independentFluxStorage) values.resize(S);
    WVFlux sharedFlux{{sharedFluxStorage[0].data(),spectralShape},
        {sharedFluxStorage[1].data(),spectralShape},{sharedFluxStorage[2].data(),spectralShape}};
    WVFlux independentFlux{{independentFluxStorage[0].data(),spectralShape},
        {independentFluxStorage[1].data(),spectralShape},
        {independentFluxStorage[2].data(),spectralShape}};
    std::vector<double> sharedRaw(4*R),independentRaw(4*R);
    WVRealFieldBundleView sharedRawView{sharedRaw.data(),{g.Nx,g.Ny,g.Nz,4}};
    WVRealFieldBundleView independentRawView{independentRaw.data(),{g.Nx,g.Ny,g.Nz,4}};
    require(bool(shared->beginStateEvaluation(state)),
        "Shared Boussinesq nonlinear scope begin failed");
    require(bool(independent->beginStateEvaluation(state)),
        "Independent Boussinesq nonlinear scope begin failed");
    shared->resetMetrics();
    independent->resetMetrics();
    require(bool(shared->nonlinearFlux(state,sharedFlux,&sharedRawView)),
        "Shared Boussinesq nonlinear preparation failed");
    require(bool(independent->nonlinearFlux(state,independentFlux,&independentRawView)),
        "Independent Boussinesq nonlinear preparation failed");
    for(std::size_t i=0;i<sharedRaw.size();++i)
        require(std::abs(sharedRaw[i]-independentRaw[i])<=1e-12*
            std::max(1.0,std::abs(independentRaw[i])),
            "Fused Boussinesq raw tendency differs from the independent dz fallback");
    for(std::size_t family=0;family<coefficients.size();++family)
        for(std::size_t i=0;i<S;++i)
            require(coefficients[family][i].real==coefficientsBefore[family][i].real &&
                coefficients[family][i].imag==coefficientsBefore[family][i].imag,
                "Fused Boussinesq nonlinear evaluation mutated input coefficients");
    for(std::size_t family=0;family<3;++family) for(std::size_t i=0;i<S;++i) {
        require(std::abs(sharedFluxStorage[family][i].real-
                    independentFluxStorage[family][i].real)<=1e-12*
                    std::max(1.0,std::abs(independentFluxStorage[family][i].real)) &&
                std::abs(sharedFluxStorage[family][i].imag-
                    independentFluxStorage[family][i].imag)<=1e-12*
                    std::max(1.0,std::abs(independentFluxStorage[family][i].imag)),
            "Shared Boussinesq nonlinear flux differs from the independent path");
    }
    require(shared->metrics().coefficientAssemblyCount==4 &&
            shared->metrics().verticalPreparationCount==8 &&
            shared->metrics().verticalOperatorExecutionCount==24 &&
            shared->metrics().horizontalSpectrumReuseCount==12 &&
            shared->metrics().preparedVerticalDerivativeCount==4 &&
            shared->metrics().derivativeAdvectionConsumerCount==12,
        "Scoped Boussinesq nonlinear evaluation did not share the expected producers");
    require(independent->metrics().coefficientAssemblyCount==16 &&
            independent->metrics().verticalPreparationCount==32 &&
            independent->metrics().derivativeAdvectionConsumerCount==12,
        "Independent Boussinesq nonlinear path did not execute all producers");
    require(bool(shared->endStateEvaluation()),"Shared Boussinesq nonlinear scope end failed");
    require(bool(independent->endStateEvaluation()),
        "Independent Boussinesq nonlinear scope end failed");

    shared->resetMetrics();
    require(bool(shared->nonlinearFlux(state,sharedFlux)),
        "Standalone shared Boussinesq nonlinear preparation failed");
    require(shared->metrics().coefficientAssemblyCount==4 &&
            shared->metrics().verticalPreparationCount==8 &&
            shared->metrics().verticalOperatorExecutionCount==24 &&
            shared->metrics().horizontalSpectrumReuseCount==12 &&
            shared->metrics().preparedVerticalDerivativeCount==4 &&
            shared->metrics().derivativeAdvectionConsumerCount==12,
        "Standalone Boussinesq nonlinear evaluation lost operation-local sharing");

    Counters failureCounters;
    std::unique_ptr<WVTransformBoussinesqKernel> failureKernel;
    require(bool(WVTransformBoussinesqKernel::create(source,
        std::make_unique<Engine>(failureCounters),failureKernel,
        WVCreateScalarMatrixBackend,sharedOptions)),
        "Boussinesq consumer-failure kernel setup failed");
    require(bool(failureKernel->beginStateEvaluation(state)),
        "Boussinesq consumer-failure scope begin failed");
    failureKernel->resetMetrics();
    failureCounters.failExecuteAt=failureCounters.executed;
    require(failureKernel->transformStateField(state,WVBoussinesqField::u,
                {sharedField.data(),spatialShape}).code==
            WVKernelStatusCode::fftExecutionFailure,
        "Injected Boussinesq horizontal inverse failure was not preserved");
    failureCounters.failExecuteAt=-1;
    require(bool(failureKernel->transformStateField(state,WVBoussinesqField::u,
        {sharedField.data(),spatialShape})),
        "Boussinesq retry after horizontal inverse failure failed");
    require(failureKernel->metrics().coefficientAssemblyCount==1 &&
            failureKernel->metrics().verticalPreparationCount==2 &&
            failureKernel->metrics().verticalOperatorExecutionCount==2 &&
            failureKernel->metrics().horizontalSpectrumReuseCount==1,
        "Boussinesq retry repeated a successfully published vertical product");
    require(bool(failureKernel->endStateEvaluation()),
        "Boussinesq consumer-failure scope end failed");
    for(auto& values:sharedFluxStorage) std::fill(values.begin(),values.end(),WVComplex64{17,19});
    failureKernel->resetMetrics();
    failureCounters.failExecuteAt=failureCounters.executed+4;
    require(failureKernel->nonlinearFlux(state,sharedFlux).code==WVKernelStatusCode::fftExecutionFailure &&
        failureKernel->metrics().derivativeAdvectionConsumerCount==0,
        "Failed Boussinesq derivative inverse published a fused consumer");
    for(const auto& values:sharedFluxStorage) for(const auto value:values)
        require(value.real==17 && value.imag==19,
            "Failed Boussinesq derivative inverse published spectral output");
}

void modalContracts(const std::filesystem::path& path) {
    std::shared_ptr<const WVStratifiedModalRecord> original; auto status=WVStratifiedModalReader::read(path.string(),original); require(bool(status),status.message.c_str());
    require(original->groups().size()==3 && original->groups()[2].columns==std::vector<std::size_t>{2},"Persisted wave group was lost");
    WVStratifiedModalInspection inspection; require(bool(WVStratifiedModalReader::inspect(path.string(),inspection)),"Wave inspection failed");
    require(inspection.scientificMatrixBytes==(4*7*3+3*(4*7*3+3*3))*sizeof(double),"Wave matrix accounting differs");
    const auto reject=[&]() {
        auto retained=original; WVStratifiedModalInspection observed; observed.scientificMatrixBytes=17;
        require(!WVStratifiedModalReader::inspect(path.string(),observed) && observed.scientificMatrixBytes==17,"Bad wave metadata accepted or partial inspection published");
        require(!WVStratifiedModalReader::read(path.string(),retained) && retained==original,"Bad wave record replaced immutable source");
    };
    const double f=2*7.2921e-5*std::sin(33*std::acos(-1.0)/180);
    struct Bad {const char* name; double value; std::size_t first,second;};
    for (const auto& bad:{Bad{"iK2unique",0,1,0},{"iK2unique",1.5,1,0},{"iK2unique",4,1,0},{"iK2unique",3,1,0},
        {"K2unique",0,1,0},{"Ppm",0,1,1},{"Qpm",NC_FILL_DOUBLE,2,2},{"h_pm",-1,1,1},{"N2",f*f,2,0}}) {
        boussinesqFixture(path); change(path,bad.name,bad.value,bad.first,bad.second); reject();
    }
    for (const auto* matrix:{"PFpmInv","QGpmInv","PFpm","QGpm","QGwg"}) {
        boussinesqFixture(path); changeWave(path,matrix,std::numeric_limits<double>::quiet_NaN(),2,1,1); reject();
        boussinesqFixture(path); { File f(path); int v; nc(nc_inq_varid(f.id,matrix,&v)); nc(nc_rename_var(f.id,v,"missing_wave_matrix")); } reject();
    }
    boussinesqFixture(path); changeWave(path,"QGpmInv",1,2,1,0); reject();
    boussinesqFixture(path); changeWave(path,"PFpmInv",2,2,0,0); reject();
    boussinesqFixture(path); changeWave(path,"QGwg",1,2,0,1); reject();
    boussinesqFixture(path);
    { File f(path); int v,z,j,k; nc(nc_inq_varid(f.id,"PFpmInv",&v)); nc(nc_rename_var(f.id,v,"old_wave")); nc(nc_inq_dimid(f.id,"z",&z)); nc(nc_inq_dimid(f.id,"j",&j)); nc(nc_inq_dimid(f.id,"K2unique",&k)); int dims[]={k,z,j}; nc(nc_def_var(f.id,"PFpmInv",NC_DOUBLE,3,dims,&v)); } reject();
    boussinesqFixture(path);
    // Allocations during the grouped scientific read must remain transactional.
    const auto name=path.string(); bool success=false;
    for (long i=0;i<512;++i) {
        auto retained=original; allocationProbe::failAfter=i; status=WVStratifiedModalReader::read(name,retained); allocationProbe::failAfter=-1;
        if (status) { success=true; break; }
        require(status.code==WVCheckpointStatusCode::allocationFailure && retained==original,"Wave read allocation failure escaped or published partial state");
    }
    require(success,"Wave record allocation sweep never succeeded");
}

#if WV_TEST_NATIVE_FFTW
void tiledNonlinearParity(const std::shared_ptr<const WVStratifiedModalRecord>& source) {
    const auto& g=source->geometry();const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    for (bool compact:{false,true}) for (std::size_t workers:{std::size_t{1},std::size_t{3},std::size_t{12}}) {
        WVVariableExecutionOptions options{WVRetainedHorizontalSchedule::streamingPrunedTile16,workers,true};
        options.pointwiseWorkers=2;options.fusedDerivativeAdvection=true;
        if (compact) options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
        const auto create=[&](bool tiled,std::unique_ptr<WVTransformBoussinesqKernel>& kernel) {
            std::unique_ptr<WVFFTEngine> engine;require(bool(WVFFTWEngine::create(1,engine)),"Native tiled engine setup failed");
            auto selected=options;selected.tiledNonlinear=tiled;
            require(bool(WVTransformBoussinesqKernel::create(source,std::move(engine),kernel,WVCreateScalarMatrixBackend,selected)),"Tiled family setup failed");
        };
        std::unique_ptr<WVTransformBoussinesqKernel> baseline,candidate;create(false,baseline);create(true,candidate);
        require(candidate->supportsTiledNonlinear() && !baseline->supportsTiledNonlinear(),"Tiled opt-in capability mismatch");
        auto referenceOptions=options;referenceOptions.tiledNonlinear=true;
        std::unique_ptr<WVTransformBoussinesqKernel> fallback;
        require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),fallback,WVCreateScalarMatrixBackend,referenceOptions)) && !fallback->supportsTiledNonlinear(),"Reference provider failed tiled fallback");
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
        std::vector<double> oldFields(4*R),newFields(4*R),oldRaw(4*R),newRaw(4*R),scratch(R);
        const WVRealFieldBundleView output{newFields.data(),{g.Nx,g.Ny,g.Nz,4}};
        WVRealFieldBundleView oldRawOutput{oldRaw.data(),{g.Nx,g.Ny,g.Nz,4}};
        WVRealFieldBundleView newRawOutput{newRaw.data(),{g.Nx,g.Ny,g.Nz,4}};
        const WVRealFieldBundleConstView borrowed{oldFields.data(),{g.Nx,g.Ny,g.Nz,4}};
        const WVBoussinesqField names[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
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
            require(cm.tiledNonlinearCount==1 && cm.tiledColumnInverseCount==g.Nz*(4+2*4) &&
                cm.tiledRowInverseCount==g.Nz*(4+3*4) && cm.tiledReusedColumnCount==g.Nz*4,"Tiled producer counts mismatch");
            require(cm.coefficientAssemblyCount==bm.coefficientAssemblyCount && cm.verticalOperatorExecutionCount==bm.verticalOperatorExecutionCount && cm.reconstructionCount==bm.reconstructionCount,"Tiled producer or vertical work changed");
            for (std::size_t i=0;i<4*R;++i) require(oldFields[i]==newFields[i],"Tiled physical field changed");
            for (std::size_t i=0;i<4*R;++i)
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
        auto badRaw=newRawOutput;badRaw.shape.fourth=3;
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
        Temporary file; boussinesqFixture(file.path);
        for (std::size_t z=0;z<7;++z) change(file.path,"dLnN2",.001*(z+1),z);
        std::shared_ptr<const WVStratifiedModalRecord> source; auto status=WVStratifiedModalReader::read(file.path.string(),source); require(bool(status),status.message.c_str());
        modalContracts(file.path);
        contracts(source);
#if WV_TEST_NATIVE_FFTW
        tiledNonlinearParity(source);
#endif
        variableScheduleParity(source);
        variableScheduleParity(source,true);
        sharedFieldGradientParity(source,false);
        sharedFieldGradientParity(source,true);
        std::unique_ptr<WVTransformBoussinesqKernel> kernel; require(bool(WVTransformBoussinesqKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),kernel)),"Lifetime setup failed");
        std::weak_ptr<const WVStratifiedModalRecord> weak=source; source.reset(); require(!weak.expired(),"Kernel lost scientific owner"); kernel.reset(); require(weak.expired(),"Scientific owner leaked");
        std::cout<<"Boussinesq kernel contracts passed\n"; return 0;
    } catch(const std::exception& e) { allocationProbe::failAfter=-1; allocationProbe::counting=false; std::cerr<<e.what()<<'\n'; return 1; }
}
