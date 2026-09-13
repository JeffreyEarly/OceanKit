#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WVReferenceFFTEngine.hpp"
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
struct Counters { int plans=0,engines=0,created=0,failAt=-1; };
class Plan final : public WVFFTPlan {
    std::unique_ptr<WVFFTPlan> plan_; Counters& counters_;
public:
    Plan(std::unique_ptr<WVFFTPlan> p,Counters& c):plan_(std::move(p)),counters_(c) { ++counters_.plans; }
    ~Plan() override { --counters_.plans; }
    WVKernelStatus execute(const void* a,void* b) override { return plan_->execute(a,b); }
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
    std::unique_ptr<WVTransformStratifiedQGKernel> kernel;
    require(bool(WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),kernel)),"Create failed");
    const auto& g=source->geometry(); const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
    std::vector<WVComplex64> a(S,{0.001,-0.002}),b(S,{17,19}); const auto saved=a;
    std::vector<double> spatial(R,29); WVComplexConstView input{a.data(),{g.Nj,g.Nkl}}; WVComplexView output{b.data(),{g.Nj,g.Nkl}};
    auto status=kernel->nonlinearFlux(input,{a.data(),input.shape}); require(status.code==WVKernelStatusCode::overlappingArrays,"Aliased flux accepted");
    status=kernel->nonlinearFlux(input,{b.data(),{S,1}}); require(status.code==WVKernelStatusCode::invalidShape,"Bad shape accepted");
    status=kernel->nonlinearFlux(input,output,std::numeric_limits<double>::quiet_NaN()); require(status.code==WVKernelStatusCode::invalidConfiguration,"NaN beta accepted");
    require(b[0].real==17 && spatial[0]==29,"Invalid call mutated output");
    status=kernel->transformA0ToField(input,static_cast<WVStratifiedQGField>(99),{spatial.data(),kernel->spatialShape()}); require(status.code==WVKernelStatusCode::unsupportedOperation,"Unknown field accepted");
    status=kernel->transformA0ToField(input,WVStratifiedQGField::u,{nullptr,kernel->spatialShape()}); require(status.code==WVKernelStatusCode::invalidPointer,"Null field accepted");
    using Closure=std::function<WVKernelStatus(WVComplexConstView,double,WVComplexView)>;
    for (const auto& operation : std::array<Closure,3>{
        [&](WVComplexConstView a,double value,WVComplexView b){return kernel->verticalDiffusivityFlux(a,value,b);},
        [&](WVComplexConstView a,double value,WVComplexView b){return kernel->linearBottomFrictionFlux(a,value,b);},
        [&](WVComplexConstView a,double value,WVComplexView b){return kernel->quadraticBottomFrictionFlux(a,value,b);}}) {
        status=operation(input,-1,output); require(status.code==WVKernelStatusCode::invalidConfiguration,"Negative closure coefficient accepted");
        status=operation(input,std::numeric_limits<double>::infinity(),output); require(status.code==WVKernelStatusCode::invalidConfiguration,"Infinite closure coefficient accepted");
        status=operation(input,1,{a.data(),input.shape}); require(status.code==WVKernelStatusCode::overlappingArrays,"Aliased closure accepted");
        require(b[0].real==17,"Invalid closure mutated output");
        require(bool(operation(input,0,output)),"Zero closure failed");
        for (const auto x:b) require(x.real==0 && x.imag==0,"Zero closure has a tendency");
        std::fill(b.begin(),b.end(),WVComplex64{17,19});
    }
    std::vector<double> tracer(R), tracerFlux(R), velocity(3*R,0), expectedTracer(R);
    const double pi=std::acos(-1.0);
    for (std::size_t z=0;z<g.Nz;++z) for(std::size_t y=0;y<g.Ny;++y) for(std::size_t x=0;x<g.Nx;++x) {
        const auto i=x+g.Nx*(y+g.Ny*z); const double scale=1+.2*z;
        tracer[i]=scale*(std::sin(6*pi*x/g.Nx)+std::cos(4*pi*y/g.Ny)+std::cos(pi*x)+std::cos(pi*y));
        velocity[i]=1; velocity[R+i]=.5;
        expectedTracer[i]=scale*(-6*pi/g.Lx*std::cos(6*pi*x/g.Nx)+2*pi/g.Ly*std::sin(4*pi*y/g.Ny));
    }
    WVRealVolumeConstView tracerInput{tracer.data(),kernel->spatialShape()};
    WVRealVolumeView tracerOutput{tracerFlux.data(),kernel->spatialShape()};
    WVRealFieldBundleConstView advection{velocity.data(),{g.Nx,g.Ny,g.Nz,3}};
    status=kernel->advectScalarWithAdvectionFields(tracerInput,advection,false,tracerOutput); require(bool(status),"Full-grid tracer derivative failed");
    for (std::size_t i=0;i<R;++i) require(std::abs(tracerFlux[i]-expectedTracer[i])<1e-14,"Full-grid tracer derivative or Nyquist handling differs");
    status=kernel->advectScalarWithAdvectionFields(tracerInput,advection,true,tracerOutput); require(bool(status),"Tracer projection failed");
    for(double value:tracerFlux) require(std::abs(value)<1e-14,"Tracer flux failed retained-mode filtering");
    status=kernel->advectScalarWithAdvectionFields(tracerInput,advection,false,{tracer.data()+1,kernel->spatialShape()});
    require(status.code==WVKernelStatusCode::overlappingArrays,"Partial tracer alias accepted");
    require(bool(kernel->nonlinearFlux(input,output)),"Flux failed");
    std::vector<double> compoundReference(R),compoundEtaZ(R),compoundEta(R),compoundResult(R);
    require(bool(kernel->beginStateEvaluation(input)),"Compound producer scope begin failed");
    require(bool(kernel->transformA0ToField(input,WVStratifiedQGField::rhoTotal,
        {compoundReference.data(),kernel->spatialShape()},WVStratifiedQGDerivative::z)),
        "Reference QG density derivative failed");
    kernel->resetMetrics();
    require(bool(kernel->transformA0ToField(input,WVStratifiedQGField::eta,
        {compoundEtaZ.data(),kernel->spatialShape()},WVStratifiedQGDerivative::z)),
        "Prepared QG eta_z failed");
    require(bool(kernel->transformA0ToField(input,WVStratifiedQGField::eta,
        {compoundEta.data(),kernel->spatialShape()})),"Prepared QG eta failed");
    require(bool(kernel->combinePreparedDensityZDerivative(WVStratifiedQGField::rhoTotal,
        {compoundEtaZ.data(),kernel->spatialShape()},{compoundEta.data(),kernel->spatialShape()},
        {compoundResult.data(),kernel->spatialShape()})),"Prepared QG density combination failed");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundResult[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Prepared QG density derivative changed values");
    require(kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::eta)]
            [static_cast<std::size_t>(WVStratifiedQGDerivative::z)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::eta)]
            [static_cast<std::size_t>(WVStratifiedQGDerivative::value)]==1 &&
        kernel->metrics().reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::rhoTotal)]
            [static_cast<std::size_t>(WVStratifiedQGDerivative::z)]==1 &&
        kernel->metrics().componentReconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::rhoTotal)]
            [static_cast<std::size_t>(WVStratifiedQGDerivative::z)][0]==1,
        "Prepared QG density derivative repeated a producer");
    require(bool(kernel->combinePreparedDensityZDerivative(WVStratifiedQGField::rhoTotal,
        {compoundEtaZ.data(),kernel->spatialShape()},{compoundEta.data(),kernel->spatialShape()},
        {compoundEtaZ.data(),kernel->spatialShape()},2)),
        "Exact prepared QG density alias was rejected");
    for(std::size_t i=0;i<R;++i)
        require(std::abs(compoundEtaZ[i]-compoundReference[i])<=1e-12*
            std::max(1.0,std::abs(compoundReference[i])),
            "Exact prepared QG density alias changed values");
    require(kernel->metrics().componentReconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::rhoTotal)]
            [static_cast<std::size_t>(WVStratifiedQGDerivative::z)][2]==1,
        "Prepared QG density component attribution failed");
    std::vector<double> partialCompound(R+1);
    std::copy_n(compoundReference.data(),R,partialCompound.data());
    require(kernel->combinePreparedDensityZDerivative(WVStratifiedQGField::rhoTotal,
        {partialCompound.data(),kernel->spatialShape()},{compoundEta.data(),kernel->spatialShape()},
        {partialCompound.data()+1,kernel->spatialShape()}).code==WVKernelStatusCode::overlappingArrays,
        "Partial prepared QG density alias was accepted");
    require(!kernel->combinePreparedDensityZDerivative(WVStratifiedQGField::u,
        {compoundEtaZ.data(),kernel->spatialShape()},{compoundEta.data(),kernel->spatialShape()},
        {compoundResult.data(),kernel->spatialShape()}),"Invalid prepared QG density target was accepted");
    require(bool(kernel->endStateEvaluation()),"Compound producer scope end failed");
    allocationProbe::calls=0; allocationProbe::counting=true;
    for (int i=0;i<5;++i) {
        require(bool(kernel->advectScalarWithAdvectionFields(tracerInput,advection,true,tracerOutput)),"Prepared tracer allocation test failed");
        require(bool(kernel->nonlinearFlux(input,output)),"Prepared flux failed");
        require(bool(kernel->verticalDiffusivityFlux(input,.002,output)),"Prepared diffusivity failed");
        require(bool(kernel->linearBottomFrictionFlux(input,1e-5,output)),"Prepared linear friction failed");
        require(bool(kernel->quadraticBottomFrictionFlux(input,.003,output)),"Prepared quadratic friction failed");
        require(bool(kernel->transformA0ToField(input,WVStratifiedQGField::rhoTotal,{spatial.data(),kernel->spatialShape()},WVStratifiedQGDerivative::z)),"Prepared density derivative failed");
        require(bool(kernel->transformQGPVToA0({spatial.data(),kernel->spatialShape()},output)),"Prepared projection failed");
    }
    allocationProbe::counting=false; require(allocationProbe::calls==0,"Prepared allocation");
    for (std::size_t i=0;i<S;++i) require(a[i].real==saved[i].real && a[i].imag==saved[i].imag,"Input mutated");
    kernel->resetMetrics();
    int evaluationOwner=0;
    require(bool(kernel->beginStateEvaluation(input,&evaluationOwner)),"Begin Stratified QG state evaluation failed");
    require(kernel->evolveA0(input,19,{a.data(),input.shape}).code==WVKernelStatusCode::overlappingArrays,
        "Scoped Stratified QG evolution mutated the active immutable state");
    require(kernel->transformQGPVToA0({spatial.data(),kernel->spatialShape()},{a.data(),input.shape}).code==WVKernelStatusCode::overlappingArrays,
        "Scoped Stratified QG projection mutated the active immutable state");
    auto registered=a;
    WVComplexConstView registeredInput{registered.data(),input.shape};
    require(bool(kernel->addStateEvaluationView(registeredInput,&evaluationOwner,2)),"Register Stratified QG state view failed");
    require(bool(kernel->transformA0ToField(registeredInput,WVStratifiedQGField::u,{spatial.data(),kernel->spatialShape()})),
        "Registered Stratified QG state view was rejected");
    require(kernel->metrics().stateValidationCount==2,"Stratified QG registered state view skipped validation");
    int foreignOwner=0;
    require(!kernel->removeStateEvaluationView(input,&evaluationOwner,0),
        "Primary Stratified QG state view was removed");
    require(!kernel->removeStateEvaluationView(registeredInput,&foreignOwner,2),
        "Foreign owner removed a Stratified QG state view");
    require(bool(kernel->removeStateEvaluationView(registeredInput,&evaluationOwner,2)),
        "Stratified QG state view removal failed");
    require(!kernel->validateStateEvaluation(registeredInput) &&
                !kernel->removeStateEvaluationView(registeredInput,&evaluationOwner,2),
        "Removed Stratified QG state view remained registered");
    require(bool(kernel->addStateEvaluationView(registeredInput,&evaluationOwner,3)) &&
                bool(kernel->validateStateEvaluation(registeredInput)),
        "Stratified QG state view storage could not be re-registered");
    require(bool(kernel->endStateEvaluation()),"End Stratified QG state evaluation failed");
    require(bool(kernel->evolveA0(input,19,{a.data(),input.shape})),"In-place stationary evolution failed");
    const auto* old=kernel.get();
    status=WVTransformStratifiedQGKernel::create({},std::make_unique<WVReferenceFFTEngine>(),kernel);
    require(!status && kernel.get()==old,"Failed setup replaced kernel");
    for (int fail=0;fail<4;++fail) {
        int count=0; Counters counters;
        status=WVTransformStratifiedQGKernel::create(source,std::make_unique<Engine>(counters),kernel,[&](std::unique_ptr<WVVerticalMatrixBackend>& backend) {
            if (count++==fail) return WVKernelStatus{WVKernelStatusCode::allocationFailure,"Injected matrix factory failure."};
            return WVCreateScalarMatrixBackend(backend);
        });
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"Matrix setup failure leaked or replaced kernel");
    }
    for (int fail=0;fail<2;++fail) {
        Counters counters; counters.failAt=fail;
        status=WVTransformStratifiedQGKernel::create(source,std::make_unique<Engine>(counters),kernel);
        require(!status && kernel.get()==old && !counters.engines && !counters.plans,"FFT setup failure leaked or replaced kernel");
    }
    bool succeeded=false;
    for (long fail=0;fail<2048;++fail) {
        Counters counters; auto engine=std::make_unique<Engine>(counters);
        allocationProbe::failAfter=fail;
        status=WVTransformStratifiedQGKernel::create(source,std::move(engine),kernel);
        allocationProbe::failAfter=-1;
        if (status) { kernel.reset(); require(!counters.engines && !counters.plans,"Successful destruction leaked"); succeeded=true; break; }
        require(kernel.get()==old && !counters.engines && !counters.plans,"Allocation failure leaked or replaced kernel");
    }
    require(succeeded,"Allocation sweep never succeeded");
}

void variableScheduleParity(const std::shared_ptr<const WVStratifiedModalRecord>& source) {
    std::unique_ptr<WVTransformStratifiedQGKernel> frozen, candidate;
    require(bool(WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),frozen)),"Frozen schedule setup failed");
    WVVariableExecutionOptions options{WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true,
        WVVariableSpectralSchedule::compactSplitFusedViews};
    options.pointwiseWorkers=2;
    require(bool(WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),candidate, WVCreateScalarMatrixBackend, options)),"Candidate schedule setup failed");
    require(std::string(candidate->horizontalScheduleIdentifier())=="full-fft-gather","Reference provider fallback was not reported");
    const auto& g=source->geometry(); const auto S=g.Nj*g.Nkl;
    std::vector<WVComplex64> input(S), frozenOut(S), candidateOut(S);
    for (std::size_t i=0;i<S;++i) input[i]={0.001*std::sin(.17*i),-0.002*std::cos(.11*i)};
    const auto inputBefore=input;
    const WVComplexConstView in{input.data(),{g.Nj,g.Nkl}};
    require(bool(frozen->nonlinearFlux(in,{frozenOut.data(),in.shape})),"Frozen nonlinear flux failed");
    require(bool(candidate->nonlinearFlux(in,{candidateOut.data(),in.shape})),"Candidate nonlinear flux failed");
    for (std::size_t i=0;i<S;++i) require(std::abs(frozenOut[i].real-candidateOut[i].real)<1e-12 && std::abs(frozenOut[i].imag-candidateOut[i].imag)<1e-12,"Candidate nonlinear flux differs");
    auto serialOptions=options; serialOptions.pointwiseWorkers=1;
    std::unique_ptr<WVTransformStratifiedQGKernel> serial; std::vector<WVComplex64> serialOut(S);
    require(bool(WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),serial,
        WVCreateScalarMatrixBackend,serialOptions)),"Single-worker candidate setup failed");
    require(bool(serial->nonlinearFlux(in,{serialOut.data(),in.shape})),"Single-worker candidate flux failed");
    for (std::size_t i=0;i<S;++i) require(serialOut[i].real==candidateOut[i].real && serialOut[i].imag==candidateOut[i].imag,
        "Pointwise worker partition changed QG arithmetic");
    const auto R=g.Nx*g.Ny*g.Nz; std::vector<double> frozenUV(2*R),candidateUV(2*R),frozenRaw(R),candidateRaw(R);
    require(bool(frozen->transformA0ToField(in,WVStratifiedQGField::u,{frozenUV.data(),{g.Nx,g.Ny,g.Nz}})) && bool(frozen->transformA0ToField(in,WVStratifiedQGField::v,{frozenUV.data()+R,{g.Nx,g.Ny,g.Nz}})),"Frozen borrowed fields failed");
    require(bool(candidate->transformA0ToField(in,WVStratifiedQGField::u,{candidateUV.data(),{g.Nx,g.Ny,g.Nz}})) && bool(candidate->transformA0ToField(in,WVStratifiedQGField::v,{candidateUV.data()+R,{g.Nx,g.Ny,g.Nz}})),"Candidate borrowed fields failed");
    const auto frozenUVBefore=frozenUV,candidateUVBefore=candidateUV;
    WVRealFieldBundleConstView frozenFields{frozenUV.data(),{g.Nx,g.Ny,g.Nz,2}},candidateFields{candidateUV.data(),{g.Nx,g.Ny,g.Nz,2}};
    WVRealVolumeView frozenRawView{frozenRaw.data(),{g.Nx,g.Ny,g.Nz}},candidateRawView{candidateRaw.data(),{g.Nx,g.Ny,g.Nz}};
    WVComplexView frozenFluxView{frozenOut.data(),in.shape},candidateFluxView{candidateOut.data(),in.shape};
    std::fill(frozenOut.begin(),frozenOut.end(),WVComplex64{17,19}); std::fill(candidateOut.begin(),candidateOut.end(),WVComplex64{17,19});
    require(bool(frozen->nonlinearFlux(in,frozenFluxView,0,&frozenRawView,&frozenFields)),"Frozen raw flux failed");
    allocationProbe::calls=0; allocationProbe::counting=true;
    require(bool(candidate->nonlinearFlux(in,candidateFluxView,0,&candidateRawView,&candidateFields)),"Candidate raw flux failed");
    allocationProbe::counting=false; require(allocationProbe::calls==0,"Candidate prepared nonlinear flux allocated");
    for (std::size_t i=0;i<R;++i) {
        require(std::abs(frozenRaw[i]-candidateRaw[i])<1e-12,
            "Candidate raw borrowed-field tendency differs");
        require(std::abs(frozenUV[i]-candidateUV[i])<1e-12 &&
            std::abs(frozenUV[R+i]-candidateUV[R+i])<1e-12,
            "Candidate borrowed velocity differs");
    }
    require(frozenUV==frozenUVBefore && candidateUV==candidateUVBefore,
        "Raw evaluation mutated borrowed fields");
    for (const auto* output:{&frozenOut,&candidateOut}) for (const auto value:*output)
        require(value.real==17 && value.imag==19,"Raw evaluation wrote spectral flux");
    for (std::size_t i=0;i<S;++i) require(input[i].real==inputBefore[i].real && input[i].imag==inputBefore[i].imag,"Parity evaluation mutated input coefficients");
    require(candidate->executionOptions().horizontalWorkers==2 && candidate->executionOptions().streamedNonlinear &&
        candidate->executionOptions().usesCompactSplitViews() && candidate->executionOptions().pointwiseWorkers==2,
        "Candidate options were not retained");
    require(candidate->persistentBytes()>=candidate->storage().workspaceBytes,"Candidate storage ledger under-reports workspace");

    auto* retained=candidate.get(); WVVariableExecutionOptions invalidWorkerOptions; invalidWorkerOptions.pointwiseWorkers=0;
    auto status=WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),candidate,
        WVCreateScalarMatrixBackend,invalidWorkerOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && candidate.get()==retained,
        "Zero pointwise workers were accepted or replaced the retained kernel");
    std::unique_ptr<WVTransformStratifiedQGKernel> invalid;
    WVVariableExecutionOptions invalidOptions;
    invalidOptions.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews;
    status=WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),invalid,
        WVCreateScalarMatrixBackend,invalidOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && !invalid,
        "Compact split views accepted a non-streaming horizontal schedule");
    invalidOptions.spectralSchedule=static_cast<WVVariableSpectralSchedule>(99);
    status=WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),invalid,
        WVCreateScalarMatrixBackend,invalidOptions);
    require(status.code==WVKernelStatusCode::invalidConfiguration && !invalid,
        "Unknown variable spectral schedule was accepted");
}
}
int main() {
    try {
        Temporary file; fixture(file.path);
        std::shared_ptr<const WVStratifiedModalRecord> source; auto status=WVStratifiedModalReader::read(file.path.string(),source); require(bool(status),status.message.c_str());
        contracts(source);
        variableScheduleParity(source);
        std::unique_ptr<WVTransformStratifiedQGKernel> kernel; require(bool(WVTransformStratifiedQGKernel::create(source,std::make_unique<WVReferenceFFTEngine>(),kernel)),"Lifetime setup failed");
        std::weak_ptr<const WVStratifiedModalRecord> weak=source; source.reset(); require(!weak.expired(),"Kernel lost scientific owner"); kernel.reset(); require(weak.expired(),"Scientific owner leaked");
        std::cout<<"Stratified QG kernel contracts passed\n"; return 0;
    } catch(const std::exception& e) { allocationProbe::failAfter=-1; allocationProbe::counting=false; std::cerr<<e.what()<<'\n'; return 1; }
}
