#include "../src/WVFieldEvaluationEventWorkspace.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVBoussinesqModalTestFixture.hpp"

#include <array>
#include <iostream>
#include <numeric>
#include <limits>
#include <sstream>
#include <type_traits>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::test_fixture;

namespace {
int injectedQGFactoryCalls = 0;
WVKernelStatus countingQGScalarBackend(std::unique_ptr<WVVerticalMatrixBackend> &backend) {
    ++injectedQGFactoryCalls;
    return WVCreateScalarMatrixBackend(backend);
}

void testInjectedQG(std::shared_ptr<const WVStratifiedModalRecord> source,
                    const WVFrozenForcingSchedule &scheduleValue) {
    auto catalog = wavevortex::runtime::test::extensionCatalog();
    WVVariableKernelServices services;
    services.matrixBackendFactory = countingQGScalarBackend;
    services.execution = {WVRetainedHorizontalSchedule::streamingPrunedTile16, 2,
                          true, WVVariableSpectralSchedule::compactSplitFusedViews};
    std::unique_ptr<WVStratifiedQGForcingEngine> baseline, injected;
    require(bool(WVStratifiedQGForcingEngine::create(
                    source, scheduleValue, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), baseline)),
            "Baseline QG injection fixture failed");
    require(bool(WVStratifiedQGForcingEngine::create(
                    source, scheduleValue, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), injected,
                    services)),
            "Injected QG fixture failed");
    require(injectedQGFactoryCalls == 4 &&
                injected->kernel().executionOptions().usesCompactSplitViews(),
            "Injected QG services were not retained");
    const auto shape = baseline->kernel().spectralShape();
    const auto count = shape.elementCount();
    std::vector<WVComplex64> coefficients(count, {.001, -.002});
    std::vector<WVComplex64> referenceValues(count), injectedValues(count);
    WVComplexConstView state{coefficients.data(), shape};
    WVComplexView reference{referenceValues.data(), shape};
    WVComplexView actual{injectedValues.data(), shape};
    require(bool(baseline->evaluateRightHandSide(state, reference)), "Baseline QG RHS failed");
    require(bool(injected->evaluateRightHandSide(state, actual)), "Injected QG RHS failed");
    for (std::size_t i = 0; i < count; ++i) {
        require(std::abs(referenceValues[i].real - injectedValues[i].real) < 1e-12,
                "Injected QG real RHS differs");
        require(std::abs(referenceValues[i].imag - injectedValues[i].imag) < 1e-12,
                "Injected QG imaginary RHS differs");
    }
    WVVariableKernelServices rejected;
    rejected.matrixBackendFactory = {};
    auto *old = injected.get();
    require(!WVStratifiedQGForcingEngine::create(
                source, scheduleValue, catalog,
                std::make_unique<WVReferenceFFTEngine>(), injected, rejected) &&
                injected.get() == old,
            "Empty QG backend factory replaced the existing engine");
}

struct FailureCounter { std::size_t calls=0,failAt=0; };
class FailingPlan final : public WVFFTPlan {
public:
    FailingPlan(std::unique_ptr<WVFFTPlan> plan,std::shared_ptr<FailureCounter> count)
        :plan_(std::move(plan)),count_(std::move(count)) {}
    WVKernelStatus execute(const void* input,void* output) override {
        if (++count_->calls==count_->failAt)
            return {WVKernelStatusCode::fftExecutionFailure,"Injected diagnostic FFT failure."};
        return plan_->execute(input,output);
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+plan_->persistentBytes(); }
private:
    std::unique_ptr<WVFFTPlan> plan_;
    std::shared_ptr<FailureCounter> count_;
};
class FailingEngine final : public WVFFTEngine {
public:
    explicit FailingEngine(std::shared_ptr<FailureCounter> count):count_(std::move(count)) {}
    std::string identifier() const override { return "forcing-diagnostic-reference"; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    WVKernelStatus createPlan(const WVFFTPlanSpecification& specification,std::unique_ptr<WVFFTPlan>& result) override {
        std::unique_ptr<WVFFTPlan> plan;
        auto status=reference_.createPlan(specification,plan); if (!status) return status;
        result=std::make_unique<FailingPlan>(std::move(plan),count_); return WVKernelStatus::ok();
    }
private:
    WVReferenceFFTEngine reference_;
    std::shared_ptr<FailureCounter> count_;
};
WVFrozenForcingSchedule schedule(std::size_t S,bool unpairedMean=false) {
    WVPortableTypedRecord fixed{"wave-vortex-forcing-configuration-v1",1,{}};
    std::vector<std::int64_t> indices(S); std::iota(indices.begin(),indices.end(),0);
    for (const char* name:{"Ap","Am","A0"}) {
        fixed.values.push_back({std::string(name)+"Indices",{S},indices});
        fixed.values.push_back({std::string(name)+"ValuesReal",{S},std::vector<double>(S,0)});
        fixed.values.push_back({std::string(name)+"ValuesImag",{S},std::vector<double>(S,0)});
    }
    if (unpairedMean) {
        fixed.values.clear();
        fixed.values={{"ApIndices",{1},std::vector<std::int64_t>{1}},
            {"ApValuesReal",{1},std::vector<double>{0}},
            {"ApValuesImag",{1},std::vector<double>{0}}};
    }
    WVPortableTypedRecord friction{"wave-vortex-forcing-configuration-v1",1,
        {{"r",{},std::vector<double>{2.5e-7}}}};
    WVPortableTypedRecord empty{"wave-vortex-forcing-configuration-v1",1,{}};
    WVPortableTypedRecord filter{"wave-vortex-forcing-configuration-v1",1,
        {{"Nj",{},std::vector<double>{unpairedMean?2.0:1.0}}}};
    WVFrozenForcingSchedule result;
    // Deliberately not in execution order, with actual user instance names.
    result.entries={{"WVFixedAmplitudeForcing",1,"held amplitudes",WVForcingStage::spectralAmplitude,255,0,"",fixed},
        {"WVAntialiasing",1,"mode filter",WVForcingStage::spectral,127,1,"",filter},
        {"WVBottomFrictionLinear",1,"bottom-drag",WVForcingStage::spatial,255,2,"",friction},
        {"WVNonlinearAdvection",1,"nonlinear advection",WVForcingStage::spatial,127,3,"",empty}};
    return result;
}

bool equal(const std::vector<WVComplex64>& a,const std::vector<WVComplex64>& b) {
    if (a.size()!=b.size()) return false;
    for (std::size_t i=0;i<a.size();++i)
        if (a[i].real!=b[i].real || a[i].imag!=b[i].imag) return false;
    return true;
}
void relative(const std::vector<double>& a,const std::vector<double>& b,const char* label="Forcing stage difference") {
    double scale=0,error=0;
    for (std::size_t i=0;i<a.size();++i) {
        scale=std::max(scale,std::abs(b[i])); error=std::max(error,std::abs(a[i]-b[i]));
    }
    if (error>1e-12*std::max(scale,1e-30)) std::cerr<<label<<": error="<<error<<" scale="<<scale<<" relative="<<error/std::max(scale,1e-30)<<"\n";
    require(error<=1e-12*std::max(scale,1e-30),"Forcing stage difference does not match its accumulated tendency");
}

template<class Engine>
void observerService(Engine& engine,WVFieldEvaluationService& fields,const WVIntegrationStateLayout& layout,
    const WVIntegrationState& state,const std::vector<WVFieldRequest>& requests,std::size_t forcingCount,
    const std::vector<std::vector<double>>& expected,const std::shared_ptr<FailureCounter>& counter) {
    auto catalog=wavevortex::runtime::test::extensionCatalog();
    WVPortableObserverRecord record;
    for(const auto& family:layout.coefficientFamilies())
        record.stateBlocks.push_back({family.identifier,WVStateScalarType::complex64,family.spectralDimensions,
            WVToleranceKind::coefficientEnergyScaled,1e-6,WVStateOwnership::integratorOwned,WVRestartRequirement::requiredDynamicState});
    WVObserverRecord first; first.identifier="forcing-fields"; first.name="forcing fields"; first.typeIdentifier="WVEulerianFields";
    for(std::size_t index=0;index<forcingCount;++index) first.fieldNames.push_back(requests[index].fieldName);
    first.fieldNames.push_back("u"); record.observers.push_back(first);
    auto second=first; second.identifier="shared-fields"; second.name="shared fields";
    second.fieldNames={requests[0].fieldName,"v"}; record.observers.push_back(second);
    auto ordinary=first; ordinary.identifier="ordinary-only"; ordinary.name="ordinary only";
    ordinary.fieldNames={"u"}; record.observers.push_back(ordinary);
    WVPortableObserverDescriptor descriptor;
    auto status=WVPortableObserverDescriptor::create(record,catalog,descriptor);
    if(!status) throw std::runtime_error("Forcing observer descriptor: "+status.message);
    std::unique_ptr<WVFieldEvaluationService> rebound;
    std::unique_ptr<WVObserverOutputEvaluationService> service;
    status=WVObserverOutputEvaluationService::create(false,descriptor,fields,service);
    if(!status) throw std::runtime_error("Forcing observer service: "+status.message);
    WVObservationSchema schema;
    require(bool(service->observationSchema(descriptor.observers()[0],schema)),"Forcing observer schema");
    require(service->metrics().sharedFieldReuseCount==2,"Forcing observer outputs were not shared");
    WVFrozenForcingSchedule scheduleValue;
    for(std::size_t index=0;index<engine.forcingCount();++index) {
        const auto* instance=engine.forcingInstance(index);
        scheduleValue.entries.push_back({instance->typeIdentifier(),instance->contractVersion(),instance->name(),
            instance->stage(),instance->priority(),instance->ordinal(),"",{}});
    }
    std::vector<WVPortableForcingVariableBinding> bindings;
    const auto preflightCalls=counter->calls;
    require(bool(catalog->forcings().diagnosticBindings(scheduleValue,fields.portableVariableConfiguration(),bindings)),
        "Data-only forcing metadata binding");
    WVObserverOutputPlanningContext context;
    context.configuration=fields.hasLegacyConfiguration() ? &fields.configuration() : nullptr;
    context.stratifiedGeometry=fields.stratifiedGeometry(); context.stateLayout=&layout;
    context.stateBlocks=record.stateBlocks.data(); context.stateBlockCount=record.stateBlocks.size();
    context.forcingConfiguration=fields.portableVariableConfiguration(); context.forcingBindings=bindings.data(); context.forcingBindingCount=bindings.size();
    WVObserverOutputPlan preflight;
    require(bool(catalog->observers().resolveOutputPlan(descriptor.observers()[0],context,preflight)),"Data-only forcing observer preflight");
    std::vector<std::uint8_t> runtimeManifest,preflightManifest;
    require(bool(encodeObservationSchemaManifest(schema,runtimeManifest)) && bool(encodeObservationSchemaManifest(preflight.schema,preflightManifest)) &&
        runtimeManifest==preflightManifest && counter->calls==preflightCalls,"Forcing runtime/preflight schemas differ or preflight executed FFTs");
    for(int invalidContract=0;invalidContract<5;++invalidContract) {
        auto invalidSchedule=scheduleValue;
        auto configuration=fields.portableVariableConfiguration();
        if(invalidContract==0) invalidSchedule.entries.front().typeIdentifier="UnknownCustomForcing";
        else if(invalidContract==1) ++invalidSchedule.entries.front().contractVersion;
        else if(invalidContract==2) invalidSchedule.entries.front().stage=WVForcingStage::spectralAmplitude;
        else if(invalidContract==3) configuration="unknown-transform";
        else invalidSchedule.entries.back().ordinal=invalidSchedule.entries.front().ordinal;
        std::vector<WVPortableForcingVariableBinding> invalidBindings;
        const auto bindingStatus=catalog->forcings().diagnosticBindings(invalidSchedule,configuration,invalidBindings);
        if(bindingStatus) {
            auto invalidContext=context;
            invalidContext.forcingBindings=invalidBindings.data(); invalidContext.forcingBindingCount=invalidBindings.size();
            WVObserverOutputPlan invalidPlan;
            require(!catalog->observers().resolveOutputPlan(descriptor.observers()[0],invalidContext,invalidPlan),
                "Forcing preflight accepted an unqualified identity, version, stage or transform");
        }
        require(counter->calls==preflightCalls,"Rejected forcing metadata preflight executed FFTs");
    }
    for(std::size_t index=0;index<forcingCount;++index) {
        const auto& variable=schema.variables[index];
        const auto& instance=fields.forcingVariableBindings()[index/(forcingCount/engine.forcingCount())];
        require(variable.name==requests[index].fieldName && variable.description.find(instance.instanceName)!=std::string::npos &&
            variable.description.find("portable_catalog_forcing")==std::string::npos && variable.layout==WVObservationValueLayout::record,
            "Observer forcing metadata did not preserve actual instance identity");
    }
    WVOutputSchedulePayload payload; require(bool(payload.reset(emptyOutputSchedulePayloadSchema())),"Forcing observer payload");
    WVPortableTypedRecord cursor;
    WVOutputGroupRecord group; group.identifier="forcing-occurrence"; group.observerIdentifiers={first.identifier,second.identifier};
    std::array<WVOutputObserverView,2> observers{{{0,&descriptor.observers()[0],descriptor.resolvedObserver(descriptor.observers()[0])},
        {1,&descriptor.observers()[1],descriptor.resolvedObserver(descriptor.observers()[1])}}};
    WVOutputRouteView route; route.observers=observers.data(); route.observerCount=observers.size(); route.scheduleOrdinal=1;
    route.proposedScheduleCursor=&cursor; route.schedulePayloadSchema=&emptyOutputSchedulePayloadSchema(); route.schedulePayload=&payload;
    route.scheduleCursorIdentity=1; route.semanticScheduleRecord=&group;
    WVOutputEvent event; event.eventOrdinal=1; event.scheduledTime=state.waveVortex.t; event.state=state; event.routes=&route; event.routeCount=1;
    const auto calls=engine.tendencyMetrics().forcingEvaluationCount,fftStart=counter->calls;
    require(service->metrics().outputCapacityBytes==0,"Forcing observer retained state-sized output arrays before an event");
    require(bool(service->prepare(event)),"Forcing observer event preparation");
    const auto fftCalls=counter->calls-fftStart;
    require(service->occurrenceWorkspaceLiveBytes()>=forcingCount*expected[0].size()*sizeof(double),"Live observer output workspace was not accounted");
    require(!service->useFieldEvaluationService(fields),"Observer rebind accepted an unfinished event");
    require(engine.tendencyMetrics().forcingEvaluationCount==calls+engine.forcingCount(),"Observers evaluated shared forcing contributions twice");
    WVObservationOccurrenceIdentity identity;
    require(bool(service->preparedOccurrenceIdentity(route,observers[0],identity)),"Forcing observer occurrence identity");
    WVObservationBatch batch;
    require(bool(service->observationBatch(identity,descriptor.observers()[0],batch)),"Forcing observer batch");
    for(std::size_t index=0;index<forcingCount;++index) {
        const auto found=std::find_if(batch.values.begin(),batch.values.end(),[&](const auto& value){return value.resolvedVariableIndex==index;});
        require(found!=batch.values.end() && found->real64Data(),"Forcing observer omitted a contribution");
        relative(std::vector<double>(found->real64Data(),found->real64Data()+found->elementCount()),expected[index],"Observer forcing field");
    }
    service->complete(event);
    require(service->occurrenceWorkspaceLiveBytes()==0 && fields.metrics().diagnosticWorkspaceLiveBytes==0,"Observer completion retained live forcing scratch");
    require(service->metrics().outputCapacityBytes==0,"Completed forcing observer retained output arrays");
    // Only the second observer is due: the first nonlinear contribution and v.
    route.observers=observers.data()+1; route.observerCount=1;
    const auto selectedCalls=engine.tendencyMetrics().forcingEvaluationCount;
    require(bool(service->prepare(event)),"Selected forcing observer preparation");
    require(engine.tendencyMetrics().forcingEvaluationCount==selectedCalls+1 &&
        service->metrics().outputCapacityBytes==2*expected[0].size()*sizeof(double),
        "Inactive forcing observers allocated output or evaluated contributions");
    require(!service->preparedOccurrenceIdentity(route,observers[0],identity),"Inactive observer exposed an occurrence");
    require(bool(service->preparedOccurrenceIdentity(route,observers[1],identity)) &&
        bool(service->observationBatch(identity,descriptor.observers()[1],batch)),"Selected forcing observer batch");
    const auto selectedValue=std::find_if(batch.values.begin(),batch.values.end(),[](const auto& value){return value.resolvedVariableIndex==0;});
    require(selectedValue!=batch.values.end() && selectedValue->real64Data(),"Selected forcing field missing");
    relative(std::vector<double>(selectedValue->real64Data(),selectedValue->real64Data()+selectedValue->elementCount()),expected[0],"Selected observer forcing");
    service->complete(event);
    WVOutputObserverView ordinaryView{2,&descriptor.observers()[2],descriptor.resolvedObserver(descriptor.observers()[2])};
    route.observers=&ordinaryView;
    const auto ordinaryCalls=engine.tendencyMetrics().forcingEvaluationCount;
    const auto ordinaryPrimitives=fields.metrics().primitiveFieldEvaluationCount;
    require(bool(service->prepare(event)),"Ordinary-only event preparation");
    require(engine.tendencyMetrics().forcingEvaluationCount==ordinaryCalls &&
        fields.metrics().primitiveFieldEvaluationCount==ordinaryPrimitives+1 &&
        service->metrics().outputCapacityBytes==expected[0].size()*sizeof(double),
        "Ordinary event evaluated inactive forcing or unrelated fields");
    service->complete(event);
    route.observerCount=0;
    const auto emptyCalls=counter->calls;
    require(bool(service->prepare(event)) && counter->calls==emptyCalls && service->metrics().outputCapacityBytes==0,
        "An event with no observers performed field work");
    service->complete(event);
    route.observers=observers.data(); route.observerCount=observers.size();
    counter->failAt=counter->calls+fftCalls;
    require(service->prepare(event).code==WVKernelStatusCode::fftExecutionFailure,"Observer expected late forcing FFT failure");
    require(service->occurrenceWorkspaceLiveBytes()==0 && service->metrics().outputCapacityBytes==0,
        "Failed observer preparation retained forcing output arrays");
    require(!service->preparedOccurrenceIdentity(route,observers[0],identity),"Failed observer exposed a stale occurrence");
    counter->failAt=0;
    require(bool(WVFieldEvaluationService::createBorrowing(engine,rebound)) && bool(service->useFieldEvaluationService(*rebound)),
        "Equivalent forcing field service could not rebind after failed preparation");
    const auto retryStatus=service->prepare(event);
    if(!retryStatus) throw std::runtime_error(
        "Observer forcing retry after rebind failed: "+retryStatus.message);
    service->complete(event);
    require(service->metrics().outputCapacityBytes==0,"Observer retry retained forcing output arrays");
}

template<class Engine,class Fields>
void fieldService(Engine& engine,const WVState& state,WVShape4D spatial,
    const Fields& reference,const std::shared_ptr<FailureCounter>& counter) {
    const bool qg=spatial.fourth==1;
    const auto R=spatial.first*spatial.second*spatial.third;
    const std::array<const char*,4> prefixes=qg ? std::array<const char*,4>{"Fqgpv_","","",""} :
        spatial.fourth==3 ? std::array<const char*,4>{"Fu_","Fv_","Feta_",""} :
        std::array<const char*,4>{"Fu_","Fv_","Fw_","Feta_"};
    std::unique_ptr<WVFieldEvaluationService> service,unbound;
    require(bool(WVFieldEvaluationService::createBorrowing(engine,service)),"Bind forcing field service");
    require(bool(WVFieldEvaluationService::createBorrowing(engine.kernel(),unbound)),"Unbound kernel service");
    WVIntegrationStateLayout layout;
    require(bool(service->createStateLayout({},layout)),"Forcing field state layout");
    std::vector<WVCoefficientFamilyConstView> families;
    const WVComplex64* statePointers[]={state.coefficients.Ap.data,state.coefficients.Am.data,state.coefficients.A0.data};
    for(std::size_t index=0;index<layout.coefficientFamilyCount();++index)
        families.push_back({&layout.coefficientFamilies()[index],statePointers[qg ? 2 : index]});
    WVIntegrationState integrationState{state,nullptr,0,families.data(),families.size()};
    std::vector<WVFieldRequest> requests;
    for(std::size_t index=0;index<engine.forcingCount();++index) {
        std::string name(engine.forcingInstance(index)->name());
        for(auto& c:name) if(c==' ' || c=='-') c='_';
        for(std::size_t channel=0;channel<spatial.fourth;++channel)
            requests.push_back({"output-"+std::to_string(requests.size()),prefixes[channel]+name,{}});
    }
    const auto forcingOutputCount=requests.size();
    requests.push_back({"repeated",requests.front().fieldName,{}});
    requests.push_back({"ordinary-u","u",{}});
    requests.push_back({"ordinary-v","v",{}});
    WVFieldSamplingRequest profile; profile.kind=WVFieldSamplingKind::fixedVerticalProfiles;
    profile.xIndices={1}; profile.yIndices={1};
    if constexpr(!std::is_same_v<Engine,WVBarotropicQGForcingEngine>)
        requests.push_back({"ordinary-profile","u",profile});
    WVFieldSamplingRequest position; position.kind=WVFieldSamplingKind::positions;
    position.x={0}; position.y={0}; position.z={0};
    requests.push_back({"ordinary-position","u",position});
    requests.push_back({"ordinary-scalar","uvMax",{}});
    WVFieldEvaluationPlan plan,rejected;
    auto status=service->createPlan(requests,plan);
    if(!status) throw std::runtime_error("Forcing field plan: "+status.message);
    require(!unbound->createPlan({requests.front()},rejected),"Kernel-only service accepted an unbound forcing");
    auto bad=requests.front(); bad.fieldName+="_missing";
    require(!service->createPlan({bad},rejected),"Unknown forcing instance accepted");
    require(plan.outputCount()==requests.size(),"Forcing field output count");
    std::vector<std::vector<double>> data(plan.outputCount());
    std::vector<WVFieldOutputView> views;
    for(std::size_t index=0;index<plan.outputCount();++index) {
        const auto& output=plan.outputs()[index];
        require(output.fieldName==requests[index].fieldName && output.samplingKind==requests[index].sampling.kind,
            "Actual forcing output identity/sampling lost");
        if(index<=forcingOutputCount) require(output.elementCount==R,"Forcing output natural shape lost");
        data[index].resize(output.elementCount,99);
        views.push_back({data[index].data(),data[index].size()});
    }
    const auto persistent=service->persistentBytes(),planBytes=plan.persistentBytes();
    const auto count=engine.tendencyMetrics().forcingEvaluationCount,start=counter->calls;
    const auto reconstructions=engine.metrics().physicalFieldReconstructionCount;
    const auto primitiveEvaluations=service->metrics().primitiveFieldEvaluationCount;
    const auto diagnosticReuse=service->metrics().diagnosticIntermediateReuseCount;
    status=service->evaluate(plan,integrationState,views.data(),views.size());
    if(!status) throw std::runtime_error("Forcing field evaluation: "+status.message);
    const auto fftCalls=counter->calls-start;
    constexpr bool stratifiedQG=std::is_same_v<Engine,WVStratifiedQGForcingEngine>;
    require(engine.metrics().physicalFieldReconstructionCount==reconstructions+(stratifiedQG ? 3 : 0),
        "Forcing diagnostics reconstructed velocity already prepared for field outputs");
    const auto expectedPhysical=qg ? 2 : std::is_same_v<Engine,WVConstantStratificationForcingEngine> ? 3 : 4;
    require(service->metrics().primitiveFieldEvaluationCount==primitiveEvaluations+expectedPhysical &&
        service->metrics().diagnosticIntermediateReuseCount==diagnosticReuse+4,
        "Ordinary and forcing outputs failed to share their physical dependencies");
    for(std::size_t index=0;index<forcingOutputCount;++index) {
        const auto instance=index/spatial.fourth,channel=index%spatial.fourth;
        relative(data[index],std::vector<double>(reference[instance].begin()+channel*R,
            reference[instance].begin()+(channel+1)*R),"Bound forcing field");
    }
    require(data[forcingOutputCount]==data.front(),"Repeated forcing channel differs");
    require(engine.tendencyMetrics().forcingEvaluationCount==count+engine.forcingCount(),
        "Field service repeated a contribution across channels/outputs");
    if(!(service->persistentBytes()==persistent &&
        plan.persistentBytes()==planBytes &&
        service->metrics().servicePersistentBytes==persistent &&
        service->metrics().diagnosticWorkspaceLiveBytes==0))
      throw std::runtime_error(
          std::string("Bound ")+
          (std::is_same_v<Engine,WVStratifiedQGForcingEngine> ? "stratified-qg" :
           std::is_same_v<Engine,WVBarotropicQGForcingEngine> ? "barotropic-qg" :
           std::is_same_v<Engine,WVConstantStratificationForcingEngine> ? "constant" :
           std::is_same_v<Engine,WVHydrostaticForcingEngine> ? "hydrostatic" : "boussinesq")+
          " field service retained diagnostic workspace or miscounted persistent bytes before="+
          std::to_string(persistent)+" after="+
          std::to_string(service->persistentBytes())+" metric="+
          std::to_string(service->metrics().servicePersistentBytes)+" live="+
          std::to_string(service->metrics().diagnosticWorkspaceLiveBytes));
    if(engine.forcingCount()>1) {
      WVFieldEvaluationPlan firstForcing,firstTwoForcings;
      require(bool(service->createPlan({requests[0]},firstForcing)) &&
          bool(service->createPlan({requests[0],requests[spatial.fourth]},
              firstTwoForcings)),
          "Incremental forcing plans failed");
      std::vector<double> firstValue(R),reusedValue(R),secondValue(R);
      WVFieldOutputView firstView{firstValue.data(),R};
      WVFieldOutputView firstTwoViews[]={{reusedValue.data(),R},
          {secondValue.data(),R}};
      const auto incrementalBefore=
          engine.tendencyMetrics().forcingEvaluationCount;
      const auto validationBefore=service->producerMetrics().stateValidations;
      WVFieldEvaluationSession session;
      auto incrementalStatus=service->beginEvaluationSession(
          integrationState,session);
      if(incrementalStatus)
        incrementalStatus=service->evaluate(firstForcing,integrationState,
            &firstView,1);
      if(incrementalStatus)
        incrementalStatus=service->evaluate(firstTwoForcings,integrationState,
            firstTwoViews,2);
      if(!incrementalStatus)
        throw std::runtime_error("Incremental forcing session failed: "+
            incrementalStatus.message);
      if(!(engine.tendencyMetrics().forcingEvaluationCount==
              incrementalBefore+2 && firstValue==data[0] &&
          reusedValue==data[0] && secondValue==data[spatial.fourth] &&
          (!qg || service->producerMetrics().stateValidations==
              validationBefore+1)))
        throw std::runtime_error(
            "Incremental forcing session replayed a completed prefix before="+
            std::to_string(incrementalBefore)+" after="+
            std::to_string(engine.tendencyMetrics().forcingEvaluationCount)+
            " validations-before="+std::to_string(validationBefore)+
            " validations-after="+
            std::to_string(service->producerMetrics().stateValidations));
    }
    WVFieldSamplingRequest forcingPosition;
    forcingPosition.kind=WVFieldSamplingKind::positions;
    forcingPosition.x={0};forcingPosition.y={0};forcingPosition.z={0};
    WVFieldEvaluationPlan sampledForcing;
    require(bool(service->createPlan({{"sampled-forcing",requests.front().fieldName,forcingPosition}},sampledForcing)),
        "Forcing diagnostic position plan failed");
    double sampledForcingValue=99;
    WVFieldOutputView sampledForcingView{&sampledForcingValue,1};
    require(bool(service->evaluate(sampledForcing,integrationState,&sampledForcingView,1)) &&
        sampledForcingValue==data.front()[(spatial.third-1)*spatial.first*spatial.second],
        "Forcing diagnostic position differs from its full-grid field");
    WVMovingFieldEvaluationPlan movingForcing;
    require(bool(service->createMovingPlan({{"moving-forcing",requests.front().fieldName,0,1,
        WVPositionInterpolation::linear}},movingForcing)),"Forcing diagnostic moving plan failed");
    double movingForcingValue=99;
    WVFieldOutputView movingForcingView{&movingForcingValue,1};
    const double forcingX=0,forcingY=0,forcingZ=0;
    const auto movingForcingStatus=service->evaluateMoving(movingForcing,
        integrationState,{&forcingX,&forcingY,&forcingZ,1},
        &movingForcingView,1);
    if(!movingForcingStatus)
      throw std::runtime_error("Forcing diagnostic moving evaluation failed: "+
          movingForcingStatus.message);
    require(
        movingForcingValue==sampledForcingValue,
        "Forcing diagnostic moving sample differs from fixed position");
    std::unique_ptr<WVFieldEvaluationService> metricService;
    require(bool(WVFieldEvaluationService::createBorrowing(engine,metricService)),
        "Forcing metric field service creation failed");
    WVMovingFieldEvaluationPlan metricMoving;
    require(bool(metricService->createMovingPlan(
        {{"metric-moving-forcing",requests.front().fieldName,0,1,
          WVPositionInterpolation::linear}},metricMoving)),
        "Forcing metric moving plan failed");
    double metricMovingValue=99;
    WVFieldOutputView metricMovingView{&metricMovingValue,1};
    require(bool(metricService->evaluateMoving(
        metricMoving,integrationState,{&forcingX,&forcingY,&forcingZ,1},
        &metricMovingView,1)) && metricMovingValue==sampledForcingValue,
        "Forcing metric moving evaluation failed");
    require(metricService->metrics().diagnosticWorkspaceLiveBytes==0 &&
        metricService->metrics().diagnosticWorkspaceHighWaterBytes>=
            R*sizeof(double)+engine.tendencyMetrics().workspaceLastPeakBytes,
        "Moving forcing metrics omit outer or nested forcing workspace");
    WVEventFieldEvaluationPlan forcingEvent;
    require(bool(service->createEventPlan({{"event-forcing",requests.front().fieldName,0,
        WVPositionInterpolation::linear}},forcingEvent)),"Forcing diagnostic event plan failed");
    WVPreparedFieldGeometry forcingGeometry;
    const WVEventPositionSetView forcingPositionSet{&forcingX,&forcingY,&forcingZ,1};
    require(bool(service->prepareEventGeometry(forcingEvent,&forcingPositionSet,1,forcingGeometry)),
        "Forcing diagnostic event geometry failed");
    double eventForcingValue=99;
    WVFieldOutputView eventForcingView{&eventForcingValue,1};
    const auto forcingEventStatus=service->evaluateEvent(forcingEvent,forcingGeometry,integrationState,&eventForcingView,1);
    require(bool(forcingEventStatus) &&
        eventForcingValue==sampledForcingValue,
        "Forcing diagnostic event sample differs from fixed position");
    WVFieldSamplingRequest forcingProfile;
    forcingProfile.kind=WVFieldSamplingKind::fixedVerticalProfiles;
    forcingProfile.xIndices={1};forcingProfile.yIndices={1};
    if constexpr(std::is_same_v<Engine,WVBarotropicQGForcingEngine>) {
        require(!service->createPlan({{"profile-forcing",requests.front().fieldName,forcingProfile}},sampledForcing),
            "Barotropic forcing accepted a vertical profile");
    } else {
        require(bool(service->createPlan({{"profile-forcing",requests.front().fieldName,forcingProfile}},sampledForcing)),
            "Forcing diagnostic profile plan failed");
        std::vector<double> forcingProfileValues(spatial.third);
        WVFieldOutputView forcingProfileView{forcingProfileValues.data(),forcingProfileValues.size()};
        require(bool(service->evaluate(sampledForcing,integrationState,&forcingProfileView,1)),
            "Forcing diagnostic profile evaluation failed");
        for(std::size_t level=0;level<spatial.third;++level)
            require(forcingProfileValues[level]==data.front()[level*spatial.first*spatial.second],
                "Forcing diagnostic profile differs from its full-grid field");
    }
    const auto successful=data;
    std::vector<std::uint8_t> active(plan.outputCount());
    std::vector<WVFieldOutputView> selectedViews(plan.outputCount());
    active[0]=1; selectedViews[0]=views[0];
    for(auto& field:data) std::fill(field.begin(),field.end(),99);
    const auto selectedCalls=engine.tendencyMetrics().forcingEvaluationCount;
    require(bool(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data())),
        "Selected diagnostic field evaluation");
    require(engine.tendencyMetrics().forcingEvaluationCount==selectedCalls+1 &&
        data[0]==successful[0],
        "Selected field lost its contribution or evaluated inactive forcings");
    for(std::size_t index=1;index<data.size();++index)
        for(auto value:data[index]) require(value==99,"Inactive diagnostic output was written");
    active[0]=0;
    const auto emptyCalls=counter->calls;
    require(bool(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data())) &&
        counter->calls==emptyCalls,"Empty diagnostic selection performed FFT work");
    // A scalar selection shares just its u/v dependencies and no forcing work.
    active.back()=1; selectedViews.back()=views.back();
    const auto scalarCalls=engine.tendencyMetrics().forcingEvaluationCount;
    const auto scalarPrimitives=service->metrics().primitiveFieldEvaluationCount;
    require(bool(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data())) &&
        data.back()==successful.back() && engine.tendencyMetrics().forcingEvaluationCount==scalarCalls &&
        service->metrics().primitiveFieldEvaluationCount==scalarPrimitives+2,
        "Selected scalar evaluated unrelated dependencies or forcing");
    active.back()=0; active[forcingOutputCount-1]=1;
    selectedViews[forcingOutputCount-1]=views[forcingOutputCount-1];
    const auto lastCalls=engine.tendencyMetrics().forcingEvaluationCount;
    const auto selectedFftStart=counter->calls;
    require(bool(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data())) &&
        data[forcingOutputCount-1]==successful[forcingOutputCount-1] &&
        engine.tendencyMetrics().forcingEvaluationCount==lastCalls+engine.forcingCount(),
        "Selecting a late forcing lost the necessary preceding stages");
    const auto selectedFftCount=counter->calls-selectedFftStart;
    for(auto& field:data) std::fill(field.begin(),field.end(),99);
    counter->failAt=counter->calls+selectedFftCount;
    require(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data()).code==WVKernelStatusCode::fftExecutionFailure,
        "Selected forcing expected late FFT failure");
    for(const auto& field:data) for(auto value:field) require(value==99,"Selected diagnostic exposed partial output on failure");
    counter->failAt=0;
    require(bool(service->evaluate(plan,integrationState,selectedViews.data(),selectedViews.size(),active.data())) &&
        data[forcingOutputCount-1]==successful[forcingOutputCount-1] && service->metrics().diagnosticWorkspaceLiveBytes==0,
        "Selected diagnostic retry differed or retained live scratch");
    WVFieldEvaluationPlan primitivePlan;
    require(bool(service->createPlan({{"selected-u","u",{}},{"inactive-v","v",{}}},primitivePlan)),"Primitive selection plan");
    std::vector<double> primitiveU(R,99);
    WVFieldOutputView primitiveViews[]={{primitiveU.data(),R},{}};
    const std::uint8_t primitiveSelection[]={1,0};
    const auto primitiveBefore=service->metrics().primitiveFieldEvaluationCount;
    require(bool(service->evaluate(primitivePlan,integrationState,primitiveViews,2,primitiveSelection)) &&
        service->metrics().primitiveFieldEvaluationCount==primitiveBefore+1 && primitiveU==successful[forcingOutputCount+1],
        "Primitive plan selection performed unrelated reconstruction or produced incorrect values");
    for(auto& field:data) std::fill(field.begin(),field.end(),99);
    const auto failurePersistent=service->persistentBytes();
    counter->failAt=counter->calls+fftCalls;
    require(service->evaluate(plan,integrationState,views.data(),views.size()).code==WVKernelStatusCode::fftExecutionFailure,
        "Bound field service expected late FFT failure");
    for(const auto& field:data) for(auto value:field) require(value==99,"Partial mixed field output escaped failure");
    require(service->metrics().diagnosticWorkspaceLiveBytes==0 &&
        service->persistentBytes()==failurePersistent,
        "Bound field failure retained workspace");
    counter->failAt=0;
    require(bool(service->evaluate(plan,integrationState,views.data(),views.size())) && data==successful,
        "Bound field retry differs from successful evaluation");
    const auto preparedChannels=static_cast<std::size_t>(expectedPhysical);
    std::vector<double> preparedValues(preparedChannels*R),result(spatial.elementCount(),99);
    WVRealFieldBundleConstView prepared{preparedValues.data(),{spatial.first,spatial.second,spatial.third,preparedChannels}};
    WVForcingTendencyOutput output{0,{result.data(),spatial}};
    const auto evaluatePrepared=[&]() {
        if constexpr(std::is_same_v<Engine,WVBarotropicQGForcingEngine> || std::is_same_v<Engine,WVStratifiedQGForcingEngine>)
            return engine.evaluateForcingTendencies(state.coefficients.A0,&output,1,&prepared);
        else return engine.evaluateForcingTendencies(state,&output,1,&prepared);
    };
    const auto rejectCalls=counter->calls;
    prepared.shape.fourth++;
    require(evaluatePrepared().code==WVKernelStatusCode::invalidShape,"Wrong prepared channel shape accepted");
    prepared.shape.fourth--;
    preparedValues[0]=std::numeric_limits<double>::quiet_NaN();
    require(evaluatePrepared().code==WVKernelStatusCode::invalidConfiguration,"Nonfinite prepared fields accepted");
    preparedValues[0]=0;
    prepared.data=reinterpret_cast<const double*>(state.coefficients.A0.data);
    require(evaluatePrepared().code==WVKernelStatusCode::overlappingArrays,"Prepared fields may alias scientific state");
    prepared.data=result.data();
    require(evaluatePrepared().code==WVKernelStatusCode::overlappingArrays,"Prepared fields may alias diagnostic output");
    require(counter->calls==rejectCalls,"Prepared input rejected after numerical execution");
    for(auto value:result) require(value==99,"Prepared input preflight changed output");
    std::vector<WVFieldRequest> mixedRequests{{"phase",qg ? "A0t" : "Apt",{}},{"velocity","u",{}},requests[0]};
    if(!qg) mixedRequests.push_back({"masked-velocity","u_w",{}});
    WVFieldEvaluationPlan mixed;
    require(bool(service->createPlan(mixedRequests,mixed)),"Mixed selected diagnostic plan");
    std::vector<WVComplex64> phase(mixed.outputs()[0].elementCount);
    std::vector<std::vector<double>> mixedReal(mixed.outputCount(),std::vector<double>(R));
    std::vector<WVFieldOutputView> mixedViews(mixed.outputCount());
    mixedViews[0]={nullptr,phase.size(),phase.data()};
    for(std::size_t index=1;index<mixedViews.size();++index) mixedViews[index]={mixedReal[index].data(),R};
    status=service->evaluate(mixed,integrationState,mixedViews.data(),mixedViews.size());
    if(!status) throw std::runtime_error("Mixed reference evaluation: "+status.message);
    const auto expectedPhase=phase; const auto expectedMixed=mixedReal;
    for(std::size_t selected=0;selected<mixed.outputCount();++selected) {
        std::vector<std::uint8_t> selection(mixed.outputCount()); selection[selected]=1;
        std::vector<WVFieldOutputView> destinations(mixed.outputCount()); destinations[selected]=mixedViews[selected];
        std::fill(phase.begin(),phase.end(),WVComplex64{99,99});
        for(auto& values:mixedReal) std::fill(values.begin(),values.end(),99);
        const auto callsBefore=engine.tendencyMetrics().forcingEvaluationCount,fftBefore=counter->calls;
        require(bool(service->evaluate(mixed,integrationState,destinations.data(),destinations.size(),selection.data())),"Mixed selected evaluation");
        require(engine.tendencyMetrics().forcingEvaluationCount==callsBefore+(selected==2 ? 1 : 0),
            "Mixed selection evaluated an inactive forcing");
        if(selected==0) require(equal(phase,expectedPhase) && counter->calls==fftBefore,"Complex-only selection reconstructed physical fields");
        else require(mixedReal[selected]==expectedMixed[selected],"Selected primitive or masked field differs");
        for(std::size_t index=1;index<mixed.outputCount();++index) if(index!=selected)
            for(auto value:mixedReal[index]) require(value==99,"Mixed selection wrote an inactive output");
    }
    WVMovingFieldEvaluationPlan moving;
    require(bool(service->createMovingPlan({{"moving-u","u",0,1,WVPositionInterpolation::linear},
        {"moving-v","v",0,1,WVPositionInterpolation::linear}},moving)),"Shared moving plan");
    double x=0,y=0,z=0; WVMovingPositionView positions{&x,&y,&z,1};
    double movingU=0,movingV=0;
    WVFieldOutputView movingViews[]={{&movingU,1},{&movingV,1}};
    require(bool(service->evaluateMoving(moving,integrationState,positions,movingViews,2)),"Independent moving reference");
    const double expectedMovingU=movingU,expectedMovingV=movingV;
    WVMovingFieldEvaluationPlan selectedMoving;
    require(bool(service->createMovingPlan({{"active-u","u",0,1,WVPositionInterpolation::linear},
        {"inactive-v","v",1,1,WVPositionInterpolation::linear},{"shared-u","u",0,1,WVPositionInterpolation::linear},
        {"inactive-eta","eta",1,1,WVPositionInterpolation::linear}},selectedMoving)),"Selected moving plan");
    const double nan=std::numeric_limits<double>::quiet_NaN();
    double selectedX[]={0,nan},selectedY[]={0,nan},selectedZ[]={0,nan};
    WVMovingPositionView selectedPositions{selectedX,selectedY,selectedZ,2};
    double firstU=99,secondU=99;
    WVFieldOutputView selectedMovingViews[]={{&firstU,1},{},{&secondU,1},{}};
    std::uint8_t movingSelection[]={1,0,1,0};
    const auto writesBefore=service->metrics().outputElementWriteCount;
    require(bool(service->evaluateMoving(selectedMoving,integrationState,selectedPositions,selectedMovingViews,4,movingSelection)) &&
        firstU==expectedMovingU && secondU==expectedMovingU && service->metrics().outputElementWriteCount==writesBefore+2,
        "Moving selection evaluated inactive coordinates or failed to share selected output");
    movingSelection[0]=movingSelection[2]=0;
    const auto emptyMoving=counter->calls;
    require(bool(service->evaluateMoving(selectedMoving,integrationState,selectedPositions,selectedMovingViews,4,movingSelection)) &&
        counter->calls==emptyMoving,"Empty moving selection reconstructed fields");
    movingSelection[1]=1; double selectedV=99; selectedMovingViews[1]={&selectedV,1};
    require(!service->evaluateMoving(selectedMoving,integrationState,selectedPositions,selectedMovingViews,4,movingSelection) &&
        counter->calls==emptyMoving && selectedV==99,"Invalid active moving coordinates escaped preflight");
    selectedX[1]=selectedY[1]=0;
    if constexpr(!std::is_same_v<Engine,WVBarotropicQGForcingEngine>)
        require(!service->evaluateMoving(selectedMoving,integrationState,selectedPositions,selectedMovingViews,4,movingSelection) &&
            counter->calls==emptyMoving && selectedV==99,"Invalid active moving z escaped preflight");
    selectedZ[1]=0;
    require(bool(service->evaluateMoving(selectedMoving,integrationState,selectedPositions,selectedMovingViews,4,movingSelection)) &&
        selectedV==expectedMovingV,"Selected moving retry differs");
    movingSelection[0]=movingSelection[2]=1; movingSelection[1]=0;
    const std::size_t velocityChannels=std::is_same_v<Engine,WVBarotropicQGForcingEngine> ? 2 : 3;
    std::vector<double> preparedVelocity(velocityChannels*R);
    std::copy(successful[forcingOutputCount+1].begin(),successful[forcingOutputCount+1].end(),preparedVelocity.begin());
    std::copy(successful[forcingOutputCount+2].begin(),successful[forcingOutputCount+2].end(),preparedVelocity.begin()+R);
    if(!qg) {
        WVFieldEvaluationPlan verticalVelocity;
        require(bool(service->createPlan({{"vertical-velocity","w",{}}},verticalVelocity)),"Prepared moving w plan");
        WVFieldOutputView vertical{preparedVelocity.data()+2*R,R};
        require(bool(service->evaluate(verticalVelocity,integrationState,&vertical,1)),"Prepared moving w field");
    }
    const WVRealFieldBundleConstView preparedMoving{preparedVelocity.data(),{spatial.first,spatial.second,spatial.third,velocityChannels}};
    const auto preparedBefore=counter->calls;
    require(bool(service->evaluateMovingFromAdvectionFields(selectedMoving,integrationState,preparedMoving,selectedPositions,selectedMovingViews,4,movingSelection)) &&
        firstU==expectedMovingU && secondU==expectedMovingU && counter->calls==preparedBefore,
        "Prepared moving selection reconstructed fields or rejected an inactive nonvelocity channel");
    WVEventFieldEvaluationPlan eventPlan;
    require(bool(service->createEventPlan({{"event-u","u",0,WVPositionInterpolation::linear},
        {"event-v","v",0,WVPositionInterpolation::linear},{"event-qgpv","qgpv",0,WVPositionInterpolation::linear}},eventPlan)),"Shared occurrence plan");
    WVPreparedFieldGeometry geometry;
    const WVEventPositionSetView positionSet{&x,&y,&z,1};
    require(bool(service->prepareEventGeometry(eventPlan,&positionSet,1,geometry)),"Shared occurrence geometry");
    std::array<std::array<double,3>,2> eventValues{};
    std::array<std::array<WVFieldOutputView,3>,2> eventViews;
    std::array<WVEventFieldEvaluationBatchEntry,2> eventEntries;
    for(std::size_t entry=0;entry<2;++entry) {
        for(std::size_t field=0;field<3;++field) eventViews[entry][field]={&eventValues[entry][field],1};
        eventEntries[entry]={&eventPlan,&geometry,eventViews[entry].data(),3};
    }
    require(bool(service->evaluateEventBatch(integrationState,eventEntries.data(),eventEntries.size())),"Independent occurrence reference");
    const auto expectedEvents=eventValues;
    WVFieldEvaluationPlan surfaceBases,surfaceAliases;
    std::vector<WVFieldRequest> aliasRequests{{"surface-height","ssh",{}}};
    if constexpr(!std::is_same_v<Engine,WVBarotropicQGForcingEngine>) {
        aliasRequests.push_back({"surface-u","ssu",{}});
        aliasRequests.push_back({"surface-v","ssv",{}});
    }
    require(bool(service->createPlan({{"base-u","u",{}},{"base-v","v",{}},{"base-pi","pi",{}}},surfaceBases)) &&
        bool(service->createPlan(aliasRequests,surfaceAliases)),"Surface sharing plans");
    std::array<std::vector<double>,3> baseValues;
    std::array<WVFieldOutputView,3> baseViews;
    for(std::size_t index=0;index<3;++index) {baseValues[index].resize(R); baseViews[index]={baseValues[index].data(),R};}
    std::vector<std::vector<double>> aliasValues(aliasRequests.size());
    std::vector<WVFieldOutputView> aliasViews(aliasRequests.size());
    for(std::size_t index=0;index<aliasRequests.size();++index) {
        aliasValues[index].resize(spatial.first*spatial.second);
        aliasViews[index]={aliasValues[index].data(),aliasValues[index].size()};
    }
    require(bool(service->evaluate(surfaceAliases,integrationState,aliasViews.data(),aliasViews.size())),"Independent surface reference");
    const auto expectedAliases=aliasValues;
    const auto retainedBefore=service->persistentBytes();
    {
        detail::WVFieldEvaluationEventScope scope(*service,integrationState);
        require(bool(scope.status()),"Shared field event scope");
        detail::WVFieldEvaluationEventScope nested(*service,integrationState);
        require(nested.status().code==WVKernelStatusCode::reentrantExecution,"Nested field event scope accepted");
        require(bool(service->evaluate(plan,integrationState,views.data(),views.size())) && data==successful,"Shared diagnostic evaluation differs");
        const auto beforeMoving=counter->calls;
        const auto reuseBefore=service->metrics().eventFieldReuseCount;
        require(bool(service->evaluateMoving(moving,integrationState,positions,movingViews,2)) &&
            movingU==expectedMovingU && movingV==expectedMovingV && counter->calls==beforeMoving &&
            service->metrics().eventFieldReuseCount>reuseBefore,"Moving fields reconstructed shared diagnostic velocities");
        const auto beforeEvents=service->metrics().primitiveFieldEvaluationCount;
        require(bool(service->evaluateEventBatch(integrationState,eventEntries.data(),eventEntries.size())) && eventValues==expectedEvents &&
            service->metrics().primitiveFieldEvaluationCount==beforeEvents+1,"Occurrence batch repeated shared fields or qgpv");
        if(!(service->metrics().eventFieldWorkspaceLiveBytes>=(qg ? 3u : 5u)*R*sizeof(double) &&
            service->persistentBytes()==retainedBefore))
          throw std::runtime_error("Shared event fields are missing from live metrics or became persistent state live="+
              std::to_string(service->metrics().eventFieldWorkspaceLiveBytes)+" retained="+
              std::to_string(retainedBefore)+" now="+std::to_string(service->persistentBytes()));
        require(bool(service->evaluate(surfaceBases,integrationState,baseViews.data(),baseViews.size())),"Shared surface base fields");
        if(service->persistentBytes()!=retainedBefore)
          throw std::runtime_error("Surface bases grew the prepared event arena retained="+
              std::to_string(retainedBefore)+" now="+
              std::to_string(service->persistentBytes())+" planned="+
              std::to_string(service->metrics().eventFieldArenaPlannedBytes)+" peak="+
              std::to_string(service->metrics().eventFieldArenaPeakBytes));
        const auto beforeAliases=counter->calls;
        const auto bytesBeforeAliases=service->metrics().eventFieldWorkspaceLiveBytes;
        require(bool(service->evaluate(surfaceAliases,integrationState,aliasViews.data(),aliasViews.size())) && aliasValues==expectedAliases &&
            counter->calls==beforeAliases && service->metrics().eventFieldWorkspaceLiveBytes==bytesBeforeAliases,
            "Surface aliases reconstructed or retained duplicate volume fields");
        require(service->persistentBytes()==retainedBefore,
            "Surface aliases grew the prepared event arena");
        // Masked coefficients must not reuse full-state fields in this event.
        require(bool(service->evaluate(mixed,integrationState,mixedViews.data(),mixedViews.size())) &&
            equal(phase,expectedPhase),"Shared complex coefficients differ");
        require(service->persistentBytes()==retainedBefore,
            "Mixed components grew the prepared event arena");
        for(std::size_t index=1;index<mixedReal.size();++index)
            require(mixedReal[index]==expectedMixed[index],"Component masks reused full-state event fields");
    }
    if(!(service->metrics().eventFieldWorkspaceLiveBytes==0 &&
        service->persistentBytes()==retainedBefore))
      throw std::runtime_error("Completed event retained shared fields retained="+
          std::to_string(retainedBefore)+" now="+
          std::to_string(service->persistentBytes()));
    // Reuse the same coefficient allocation after changing its contents between events.
    std::vector<std::vector<WVComplex64>> changedCoefficients(families.size());
    auto changedFamilies=families;
    for(std::size_t family=0;family<families.size();++family) {
        changedCoefficients[family].assign(families[family].data,families[family].data+families[family].layout->elementCount);
        changedFamilies[family].data=changedCoefficients[family].data();
    }
    auto changedState=integrationState; changedState.coefficientFamilies=changedFamilies.data();
    const auto changedShape=state.coefficients.A0.shape;
    if(qg) changedState.waveVortex.coefficients={{},{},{changedCoefficients[0].data(),changedShape}};
    else changedState.waveVortex.coefficients={{changedCoefficients[0].data(),changedShape},
        {changedCoefficients[1].data(),changedShape},{changedCoefficients[2].data(),changedShape}};
    std::vector<double> changedU(R);
    WVFieldOutputView changedViews[]={{changedU.data(),R},{}};
    const std::uint8_t onlyU[]={1,0};
    for(double scale:{1.0,0.5}) {
        if(scale==0.5) for(auto& family:changedCoefficients) for(auto& value:family) {value.real*=0.5; value.imag*=0.5;}
        detail::WVFieldEvaluationEventScope scope(*service,changedState);
        const auto before=counter->calls;
        require(bool(service->evaluate(primitivePlan,changedState,changedViews,2,onlyU)) && counter->calls>before,
            "A new event reused fields from an earlier state");
        auto expected=successful[forcingOutputCount+1]; for(auto& value:expected) value*=scale;
        relative(changedU,expected,"Changed event coefficients");
        auto shifted=changedState; shifted.waveVortex.t+=1;
        const auto shiftedBefore=counter->calls;
        require(!service->evaluate(primitivePlan,shifted,changedViews,2,onlyU) &&
            counter->calls==shiftedBefore,
            "Active immutable event accepted a changed time");
        const auto originalBefore=counter->calls;
        require(bool(service->evaluate(primitivePlan,changedState,changedViews,2,onlyU)) && counter->calls==originalBefore,
            "A different phase replaced the original event fields");
        relative(changedU,expected,"Original event fields after a phase change");
    }
    for(auto& field:data) std::fill(field.begin(),field.end(),99);
    {
        detail::WVFieldEvaluationEventScope scope(*service,integrationState);
        counter->failAt=counter->calls+fftCalls;
        require(service->evaluate(plan,integrationState,views.data(),views.size()).code==WVKernelStatusCode::fftExecutionFailure,
            "Shared-field event expected late FFT failure");
        for(const auto& field:data) for(auto value:field) require(value==99,"Shared-field failure exposed partial diagnostics");
    }
    counter->failAt=0;
    require(service->metrics().eventFieldWorkspaceLiveBytes==0,"Failed event retained shared fields");
    {
        detail::WVFieldEvaluationEventScope scope(*service,integrationState);
        require(bool(service->evaluate(plan,integrationState,views.data(),views.size())) && data==successful,"Shared event retry differs");
    }
    if constexpr(std::is_same_v<Engine,WVConstantStratificationForcingEngine>) {
        WVFieldEvaluationPlan vorticity;
        require(bool(service->createPlan({{"zx","zeta_x",{}},{"zy","zeta_y",{}},{"zz","zeta_z",{}}},vorticity)),"Vorticity sharing plan");
        std::array<std::vector<double>,3> values;
        std::array<WVFieldOutputView,3> destinations;
        for(std::size_t channel=0;channel<3;++channel) {
            values[channel].resize(R); destinations[channel]={values[channel].data(),R};
        }
        require(bool(service->evaluate(vorticity,integrationState,destinations.data(),3)),"Independent vorticity reference");
        const auto reference=values;
        const auto retained=service->persistentBytes();
        {
            detail::WVFieldEvaluationEventScope scope(*service,integrationState);
            const std::uint8_t verticalOnly[]={0,0,1},horizontalOnly[]={1,1,0};
            require(bool(service->evaluate(vorticity,integrationState,destinations.data(),3,verticalOnly)),"Shared vertical vorticity");
            require(service->metrics().eventFieldWorkspaceLiveBytes>=6*R*sizeof(double),"Vertical vorticity derivative storage is missing");
            const auto reuse=service->metrics().eventFieldReuseCount;
            require(bool(service->evaluate(vorticity,integrationState,destinations.data(),3,horizontalOnly)) &&
                service->metrics().eventFieldReuseCount==reuse+2 && values==reference,
                "Horizontal vorticity did not reuse the same-event u/v derivatives");
            const auto calls=counter->calls;
            require(bool(service->evaluate(vorticity,integrationState,destinations.data(),3)) && counter->calls==calls && values==reference &&
                service->metrics().eventFieldWorkspaceLiveBytes>=9*R*sizeof(double),
                "Repeated vorticity reconstructed derivatives or retained unexpected storage");
        }
        require(service->metrics().eventFieldWorkspaceLiveBytes==0 && service->persistentBytes()==retained,
            "Vorticity derivatives survived their output event");
    }
    observerService(engine,*service,layout,integrationState,requests,forcingOutputCount,successful,counter);
}

template<class Engine>
void nonlinearVorticityEvaluationSession(Engine& engine,const WVState& state,
    WVShape4D spatial) {
    const auto R=spatial.first*spatial.second*spatial.third;
    std::unique_ptr<WVFieldEvaluationService> service;
    require(bool(WVFieldEvaluationService::createBorrowing(engine,service)),
        "Bind nonlinear/vorticity field service");
    WVIntegrationStateLayout layout;
    require(bool(service->createStateLayout({},layout)),
        "Nonlinear/vorticity state layout");
    const WVComplex64* statePointers[]={state.coefficients.Ap.data,
        state.coefficients.Am.data,state.coefficients.A0.data};
    std::vector<WVCoefficientFamilyConstView> families;
    for(std::size_t index=0;index<layout.coefficientFamilyCount();++index)
        families.push_back({&layout.coefficientFamilies()[index],statePointers[index]});
    const WVIntegrationState integrationState{state,nullptr,0,
        families.data(),families.size()};

    WVFieldEvaluationPlan nonlinear,vorticity;
    require(bool(service->createPlan({{"nonlinear","Fu_nonlinear_advection",{}}},
            nonlinear)) &&
        bool(service->createPlan({{"zeta-x","zeta_x",{}},
            {"zeta-y","zeta_y",{}}},vorticity)),
        "Nonlinear/vorticity plans");
    struct Result {
        std::vector<double> nonlinear,zetaX,zetaY;
    } reference;
    bool hasReference=false;
    for(auto policy:{WVVariableEvaluationPolicy::reuse,
            WVVariableEvaluationPolicy::lowMemory}) {
        require(bool(service->setVariableEvaluationPolicy(policy)),
            "Nonlinear/vorticity policy");
        for(bool vorticityFirst:{false,true}) {
            Result actual{{},std::vector<double>(R),std::vector<double>(R)};
            actual.nonlinear.resize(nonlinear.outputs()[0].elementCount);
            WVFieldOutputView nonlinearView{actual.nonlinear.data(),
                actual.nonlinear.size()};
            WVFieldOutputView vorticityViews[]={{actual.zetaX.data(),R},
                {actual.zetaY.data(),R}};
            const auto before=service->producerMetrics();
            const auto ledgerBefore=service->metrics().variableEvaluation;
            {
                WVFieldEvaluationSession session;
                require(bool(service->beginEvaluationSession(integrationState,session)),
                    "Nonlinear/vorticity evaluation session");
                if(vorticityFirst) {
                    auto status=service->evaluate(vorticity,integrationState,
                        vorticityViews,2);
                    if(!status) throw std::runtime_error("Vorticity first: "+status.message);
                    status=service->evaluate(nonlinear,integrationState,
                        &nonlinearView,1);
                    if(!status) throw std::runtime_error("Nonlinear second: "+status.message);
                } else {
                    auto status=service->evaluate(nonlinear,integrationState,
                        &nonlinearView,1);
                    if(!status) throw std::runtime_error("Nonlinear first: "+status.message);
                    status=service->evaluate(vorticity,integrationState,
                        vorticityViews,2);
                    if(!status) throw std::runtime_error("Vorticity second: "+status.message);
                }
            }
            const auto after=service->producerMetrics();
            const auto ledgerAfter=service->metrics().variableEvaluation;
            const auto delta=[&](std::size_t field,std::size_t axis,
                    std::size_t component=0) {
                return after.reconstructions[field][axis][component]-
                    before.reconstructions[field][axis][component];
            };
            const bool low=policy==WVVariableEvaluationPolicy::lowMemory;
            const auto sharedExpected=low ? 2u : 1u;
            require(delta(0,3)==sharedExpected && delta(1,3)==sharedExpected &&
                    delta(2,2)==(spatial.fourth==3 ? 1u : sharedExpected) &&
                    delta(2,1)==(spatial.fourth==3 ? 1u : sharedExpected),
                "Nonlinear/vorticity session repeated or lost a total-component derivative producer");
            for(std::size_t component=1;component<5;++component)
                require(delta(0,3,component)==0 && delta(1,3,component)==0 &&
                        delta(2,2,component)==0 && delta(2,1,component)==0,
                    "Total nonlinear/vorticity outputs contaminated a distinct component producer");
            const auto sharedCount=spatial.fourth==3 ? 2u : 4u;
            if(low) {
                const auto recomputations=ledgerAfter.recomputations-
                    ledgerBefore.recomputations;
                const auto evictions=ledgerAfter.evictions-ledgerBefore.evictions;
                // Variable-N kernels expose exact derivative keys, so only
                // u_z/v_z (and Boussinesq w_x/w_y) are recomputed. The constant
                // kernel retains its established fused derivative group counts.
                const auto expectedRecomputations=
                    std::is_same_v<Engine,
                        WVConstantStratificationForcingEngine> ?
                    (vorticityFirst ? (spatial.fourth==3 ? 6u : 9u) :
                        (spatial.fourth==3 ? 2u : 3u)) : sharedCount;
                if(recomputations!=expectedRecomputations ||
                    evictions<recomputations)
                    throw std::runtime_error("Low-memory nonlinear/vorticity groups: channels="+
                        std::to_string(spatial.fourth)+" recomputations="+
                        std::to_string(recomputations)+" evictions="+
                        std::to_string(evictions)+" reverse="+
                        std::to_string(vorticityFirst));
            } else {
                require(ledgerAfter.recomputations==ledgerBefore.recomputations &&
                        ledgerAfter.duplicateExecutions==ledgerBefore.duplicateExecutions &&
                        ledgerAfter.cacheHits-ledgerBefore.cacheHits>=sharedCount,
                    "Reuse nonlinear/vorticity session did not cache every shared derivative");
            }
            if(!hasReference) {reference=actual; hasReference=true;}
            else require(actual.nonlinear==reference.nonlinear &&
                    actual.zetaX==reference.zetaX && actual.zetaY==reference.zetaY,
                "Nonlinear/vorticity query order or policy changed values");
        }
    }
    require(bool(service->setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::reuse)),
        "Nonlinear/vorticity reuse restore");

    WVFieldEvaluationPlan laterChannels;
    std::vector<WVFieldRequest> laterRequests{{"later-v",
        "Fv_nonlinear_advection",{}},{"later-eta",
        "Feta_nonlinear_advection",{}}};
    if(spatial.fourth==4)
        laterRequests.push_back({"later-w","Fw_nonlinear_advection",{}});
    std::vector<double> firstValues(nonlinear.outputs()[0].elementCount);
    WVFieldOutputView firstView{firstValues.data(),firstValues.size()};
    const auto nonlinearBefore=engine.metrics().nonlinearProducerCount;
    {
        WVFieldEvaluationSession session;
        require(bool(service->beginEvaluationSession(integrationState,session)),
            "Active forcing demand session");
        require(bool(service->evaluate(nonlinear,integrationState,&firstView,1)),
            "Active forcing demand Fu query");
        const auto status=service->createPlanForActiveEvaluation(
            laterRequests,laterChannels);
        if(!status) throw std::runtime_error(
            "Active forcing channel preparation: "+status.message);
        std::vector<std::vector<double>> values(laterChannels.outputCount());
        std::vector<WVFieldOutputView> views(laterChannels.outputCount());
        for(std::size_t index=0;index<views.size();++index) {
            values[index].resize(laterChannels.outputs()[index].elementCount);
            views[index]={values[index].data(),values[index].size()};
        }
        require(bool(service->evaluate(laterChannels,integrationState,
                    views.data(),views.size())),
            "Active forcing demand later-channel query");
        for(const auto& field:values) for(const auto value:field)
            require(std::isfinite(value),
                "Active forcing demand produced a nonfinite channel");
    }
    require(engine.metrics().nonlinearProducerCount==nonlinearBefore+1,
        "Active forcing channels repeated the nonlinear producer");
}

template<class Engine>
void qgAdaptiveMaximumEvaluationSession(Engine& engine,
    WVComplexConstView A0) {
    std::size_t adaptiveIndex=engine.forcingCount();
    for(std::size_t index=0;index<engine.forcingCount();++index)
        if(engine.forcingInstance(index)->typeIdentifier()=="WVAdaptiveDamping")
            adaptiveIndex=index;
    require(adaptiveIndex<engine.forcingCount(),
        "QG adaptive damping fixture is missing");
    std::string adaptiveName(engine.forcingInstance(adaptiveIndex)->name());
    for(auto& c:adaptiveName) if(c==' ' || c=='-') c='_';

    std::unique_ptr<WVFieldEvaluationService> service;
    require(bool(WVFieldEvaluationService::createBorrowing(engine,service)),
        "Bind QG adaptive/maximum field service");
    WVIntegrationStateLayout layout;
    require(bool(service->createStateLayout({},layout)) &&
            layout.coefficientFamilyCount()==1,
        "QG adaptive/maximum state layout");
    const WVCoefficientFamilyConstView family{
        &layout.coefficientFamilies()[0],A0.data};
    const WVState state{0,0,{{},{},A0}};
    const WVIntegrationState integrationState{state,nullptr,0,&family,1};
    WVFieldEvaluationPlan adaptive,maximum;
    require(bool(service->createPlan({{"adaptive",
                "Fqgpv_"+adaptiveName,{}}},adaptive)) &&
            bool(service->createPlan({{"maximum","uvMax",{}}},maximum)),
        "QG adaptive/maximum plans");

    std::vector<double> adaptiveReference,maximumReference;
    bool hasReference=false;
    const auto reductionCount=[&]() {
        auto count=service->producerMetrics().horizontalSpeedReductions;
        if constexpr(std::is_same_v<Engine,WVStratifiedQGForcingEngine>)
            count+=engine.metrics().horizontalSpeedMaximumReductionCount;
        return count;
    };
    for(auto policy:{WVVariableEvaluationPolicy::reuse,
            WVVariableEvaluationPolicy::lowMemory}) {
        require(bool(service->setVariableEvaluationPolicy(policy)),
            "QG adaptive/maximum policy");
        for(bool maximumFirst:{false,true}) {
            std::vector<double> adaptiveValues(
                adaptive.outputs()[0].elementCount),maximumValues(1);
            WVFieldOutputView adaptiveView{adaptiveValues.data(),
                adaptiveValues.size()};
            WVFieldOutputView maximumView{maximumValues.data(),1};
            const auto reductionsBefore=reductionCount();
            const auto ledgerBefore=service->metrics().variableEvaluation;
            {
                WVFieldEvaluationSession session;
                require(bool(service->beginEvaluationSession(
                            integrationState,session)),
                    "QG adaptive/maximum evaluation session");
                auto first=maximumFirst ?
                    service->evaluate(maximum,integrationState,&maximumView,1) :
                    service->evaluate(adaptive,integrationState,&adaptiveView,1);
                if(!first) throw std::runtime_error(
                    "QG adaptive/maximum first query: "+first.message);
                auto second=maximumFirst ?
                    service->evaluate(adaptive,integrationState,&adaptiveView,1) :
                    service->evaluate(maximum,integrationState,&maximumView,1);
                if(!second) throw std::runtime_error(
                    "QG adaptive/maximum second query: "+second.message);
            }
            const auto ledgerAfter=service->metrics().variableEvaluation;
            const auto reductionsAfter=reductionCount();
            const bool streamedScalar=
                policy==WVVariableEvaluationPolicy::lowMemory;
            const auto expectedReductions=streamedScalar ? 2u : 1u;
            const bool recomputationRecorded=streamedScalar ?
                ledgerAfter.recomputations>ledgerBefore.recomputations :
                ledgerAfter.recomputations==ledgerBefore.recomputations;
            if(reductionsAfter!=reductionsBefore+expectedReductions ||
                !recomputationRecorded)
                throw std::runtime_error(std::string("QG adaptive/maximum family=")+
                    (std::is_same_v<Engine,WVBarotropicQGForcingEngine> ? "barotropic" : "stratified")+
                    " policy="+variableEvaluationPolicyIdentifier(policy)+
                    " reverse="+std::to_string(maximumFirst)+
                    " reduction="+std::to_string(reductionsAfter-reductionsBefore)+
                    " recomputation="+std::to_string(ledgerAfter.recomputations-ledgerBefore.recomputations));
            if(!hasReference) {
                adaptiveReference=adaptiveValues;
                maximumReference=maximumValues;
                hasReference=true;
            } else require(adaptiveValues==adaptiveReference &&
                    maximumValues==maximumReference,
                "QG adaptive/maximum query order or policy changed values");
        }
    }
    require(bool(service->setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::reuse)),
        "QG adaptive/maximum reuse restore");
}

template<class Engine,class Project,class Reconstruct>
void exercise(Engine& engine,WVShape2D spectral,WVShape4D spatial,
    const std::shared_ptr<FailureCounter>& counter,Project project,Reconstruct reconstruct,bool unpairedMean=false) {
    const auto S=spectral.elementCount(),R=spatial.elementCount();
    std::vector<WVComplex64> coefficients(3*S),rhs(3*S),projected(3*S);
    for (std::size_t i=0;i<coefficients.size();++i)
        coefficients[i]={1e-4*std::sin(.17*(i+1)),1e-4*std::cos(.23*(i+1))};
    WVState state{83,17,{{coefficients.data(),spectral},{coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}}};
    if constexpr (!std::is_same_v<Engine,WVConstantStratificationForcingEngine>)
        require(bool(engine.kernel().constrainCoefficients({{coefficients.data(),spectral},
            {coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}})),"Fixture coefficient constraints");
    WVFlux flux{{rhs.data(),spectral},{rhs.data()+S,spectral},{rhs.data()+2*S,spectral}};
    const auto baselineStatus=engine.nonlinearFlux(state,flux);
    if (!baselineStatus) throw std::runtime_error("Baseline forcing RHS: "+baselineStatus.message);
    const auto stateBefore=coefficients,rhsBefore=rhs;
    const auto bytesBefore=engine.persistentBytes();
    const auto evaluationsBefore=engine.metrics().evaluationCount;
    const auto reconstructionsBefore=engine.metrics().physicalFieldReconstructionCount;
    std::array<std::vector<double>,4> fields;
    for (auto& field:fields) field.resize(R,99);
    WVForcingTendencyOutput outputs[]={ {3,{fields[3].data(),spatial}}, {2,{fields[2].data(),spatial}},
        {0,{fields[0].data(),spatial}},{1,{fields[1].data(),spatial}} };
    require(engine.forcingCount()==4 && engine.forcingInstance(0)->name()=="nonlinear advection" &&
        engine.forcingInstance(1)->name()=="bottom-drag" && engine.forcingInstance(2)->name()=="mode filter" && engine.forcingInstance(3)->name()=="held amplitudes" &&
        !engine.forcingInstance(4),"Diagnostic binding does not use the resolved stage order and instance names");
    const auto callStart=counter->calls;
    auto status=engine.evaluateForcingTendencies(state,outputs,4);
    if (!status) throw std::runtime_error(status.message);
    const auto diagnosticFFTCalls=counter->calls-callStart;
    require(equal(coefficients,stateBefore) && engine.metrics().evaluationCount==evaluationsBefore,
        "Diagnostics changed model coefficients or RHS evaluation count");
    require(engine.metrics().physicalFieldReconstructionCount==reconstructionsBefore+1,
        "Same-event forcings did not share physical fields");
    require(engine.tendencyMetrics().forcingEvaluationCount==4 &&
        engine.tendencyMetrics().spatialProjectionCount==1 &&
        engine.tendencyMetrics().spectralReconstructionCount==2 &&
        engine.tendencyMetrics().workspaceLiveBytes==0 &&
        engine.tendencyMetrics().workspaceHighWaterBytes>0 && engine.persistentBytes()==bytesBefore,
        "Forcing diagnostics repeated contributions or retained event workspace");
    std::vector<double> sum(R),expected(R);
    for (std::size_t i=0;i<R;++i) sum[i]=fields[0][i]+fields[1][i];
    require(bool(project(state,sum,projected)),"Independent accumulated tendency projection failed");
    auto removed=projected;
    for (std::size_t i=0;i<projected.size();++i) {
        const auto mode=(i%S)/spectral.rows;
        bool filtered=i%spectral.rows>(unpairedMean?1U:0U);
        if constexpr (std::is_same_v<Engine,WVConstantStratificationForcingEngine>) {
            const auto& descriptor=engine.kernel().descriptor();
            const auto& c=descriptor.configuration();
            filtered|=descriptor.fourierModes()[mode].Kh>(2.0/3.0)*2*std::acos(-1.0)*(c.Nx/2)/c.Lx;
        } else {
            const auto& g=engine.kernel().geometry();
            filtered|=std::hypot(g.k[mode],g.l[mode])>(2.0/3.0)*2*std::acos(-1.0)*(g.Nx/2)/g.Lx;
        }
        removed[i]=filtered ? WVComplex64{-projected[i].real,-projected[i].imag} : WVComplex64{};
        projected[i]=filtered ? WVComplex64{} : WVComplex64{-projected[i].real,-projected[i].imag};
    }
    require(bool(reconstruct(state,removed,expected)),"Independent filter diagnostic reconstruction failed");
    relative(fields[2],expected);
    if (unpairedMean) {
        // A forcing may independently hold just Ap at the horizontal mean.
        // Its real u/v are half the corresponding complete inertial pair.
        const auto ap=projected[1];
        std::fill(projected.begin(),projected.end(),WVComplex64{});
        projected[1]={.5*ap.real,.5*ap.imag};
        projected[S+1]={.5*ap.real,-.5*ap.imag};
    }
    require(bool(reconstruct(state,projected,expected)),"Independent amplitude diagnostic reconstruction failed");
    relative(fields[3],expected);
    require(bool(engine.nonlinearFlux(state,flux)) && equal(rhs,rhsBefore) && equal(coefficients,stateBefore),
        "RHS or scientific state changed after requesting diagnostics");

    auto invalid=outputs[0]; invalid.executionIndex=4;
    require(engine.evaluateForcingTendencies(state,&invalid,1).code==WVKernelStatusCode::invalidConfiguration,"Unbound forcing index accepted");
    invalid=outputs[0]; invalid.fields.data=reinterpret_cast<double*>(coefficients.data());
    require(engine.evaluateForcingTendencies(state,&invalid,1).code==WVKernelStatusCode::overlappingArrays,"Diagnostic output may not alias state");
    invalid=outputs[0]; invalid.fields.shape.fourth++;
    require(engine.evaluateForcingTendencies(state,&invalid,1).code==WVKernelStatusCode::invalidShape,"Wrong diagnostic channel count accepted");
    auto duplicate=std::array<WVForcingTendencyOutput,2>{outputs[0],outputs[0]};
    require(!engine.evaluateForcingTendencies(state,duplicate.data(),duplicate.size()),"Duplicate diagnostic contributions accepted");

    const auto reference=fields;
    for (auto& field:fields) std::fill(field.begin(),field.end(),99);
    // Fail after the spatial contributions have been evaluated, during the
    // final amplitude reconstruction. No partially prepared output may escape.
    counter->failAt=counter->calls+diagnosticFFTCalls;
    status=engine.evaluateForcingTendencies(state,outputs,4);
    require(status.code==WVKernelStatusCode::fftExecutionFailure,"Expected late diagnostic FFT failure");
    for (const auto& field:fields) for (auto value:field) require(value==99,"Failed diagnostic evaluation changed user output");
    require(equal(coefficients,stateBefore) && engine.tendencyMetrics().workspaceLiveBytes==0 && engine.persistentBytes()==bytesBefore,
        "Failed diagnostic evaluation changed state or retained workspace");
    counter->failAt=0;
    require(bool(engine.evaluateForcingTendencies(state,outputs,4)) && fields==reference,
        "Diagnostic retry did not reproduce the successful result");
    std::vector<double> selected(R);
    const WVForcingTendencyOutput single{3,{selected.data(),spatial}};
    require(bool(engine.evaluateForcingTendencies(state,&single,1)) && selected==reference[3],
        "Requesting only the amplitude contribution lost preceding stages");
    nonlinearVorticityEvaluationSession(engine,state,spatial);
    engine.setLinearDynamics(true);
    require(bool(engine.nonlinearFlux(state,flux)) &&
        equal(rhs,std::vector<WVComplex64>(rhs.size())) && equal(coefficients,stateBefore),
        "Linear coefficient RHS must vanish without discarding diagnostic instances");
    auto constrained=coefficients;
    WVMutableCoefficients constraintView{{constrained.data(),spectral},
        {constrained.data()+S,spectral},{constrained.data()+2*S,spectral}};
    const auto constraint=engine.restoreForcingAmplitudes(constraintView);
    require(bool(constraint) && constraint.modifiedCoefficientCount>0 && !equal(constrained,coefficients),
        "Linear evolution discarded fixed-amplitude constraints");
    if constexpr (std::is_same_v<Engine,WVConstantStratificationForcingEngine>) {
        const auto volume=engine.kernel().descriptor().spatialShape();
        std::vector<double> velocity(3*volume.elementCount());
        WVRealFieldBundleView view{velocity.data(),{volume.first,volume.second,volume.third,3}};
        WVConstantStratificationRightHandSideContext context;
        require(bool(engine.beginStateEvaluation(state)),"Linear advection scope failed");
        require(bool(engine.evaluateRightHandSideWithContext(state,flux,view,context)) &&
            context.advectionFields().data==velocity.data() &&
            equal(rhs,std::vector<WVComplex64>(rhs.size())) &&
            std::any_of(velocity.begin(),velocity.end(),[](double value){ return value!=0; }),
            "Linear evolution lost the shared advection context");
        engine.endStateEvaluation();
    }
    fieldService(engine,state,spatial,reference,counter);
    engine.setLinearDynamics(false);
    require(bool(engine.nonlinearFlux(state,flux)) && equal(rhs,rhsBefore) && equal(coefficients,stateBefore),
        "Field-service diagnostics changed subsequent RHS or scientific state");
}

void constant(bool hydrostatic) {
    WVTransformConstantStratificationConfiguration c;
    c.Nx=6; c.Ny=5; c.Nz=7; c.Nj=4; c.Lx=15000; c.Ly=12000; c.Lz=1300;
    c.N0=5.2e-3; c.rho0=1027; c.g=9.80665; c.planetaryRadius=6.3712e6;
    c.rotationRate=7.292115e-5; c.latitude=33; c.isHydrostatic=hydrostatic; c.shouldAntialias=false;
    std::unique_ptr<WVConstantStratificationForcingEngine> initial;
    auto catalog=wavevortex::runtime::test::extensionCatalog();
    require(bool(WVConstantStratificationForcingEngine::create(c,{},catalog,std::make_unique<WVReferenceFFTEngine>(),initial)),"Constant setup");
    const auto spectral=initial->stateShape(); const auto S=spectral.elementCount(),R=c.Nx*c.Ny*c.Nz;
    auto counter=std::make_shared<FailureCounter>();
    std::unique_ptr<WVConstantStratificationForcingEngine> engine;
    require(bool(WVConstantStratificationForcingEngine::create(c,schedule(S),catalog,std::make_unique<FailingEngine>(counter),engine)),"Constant forcing setup");
    const WVShape4D spatial{c.Nx,c.Ny,c.Nz,hydrostatic?3U:4U};
    exercise(*engine,spectral,spatial,counter,
        [&](const WVState& state,const std::vector<double>& sum,std::vector<WVComplex64>& out) {
            WVMutableCoefficients output{{out.data(),spectral},{out.data()+S,spectral},{out.data()+2*S,spectral}};
            return hydrostatic ? engine->kernel().transformUVEtaToWaveVortex({sum.data(),spatial},state.t,state.t0,output) :
                engine->kernel().transformUVWEtaToWaveVortex({sum.data(),spatial},state.t,state.t0,output);
        },
        [&](const WVState& state,const std::vector<WVComplex64>& delta,std::vector<double>& out) {
            WVState d{state.t,state.t0,{{delta.data(),spectral},{delta.data()+S,spectral},{delta.data()+2*S,spectral}}};
            std::vector<double> values(4*R); WVRealFieldBundleView result{values.data(),{c.Nx,c.Ny,c.Nz,4}};
            auto status=engine->kernel().transformWaveVortexToUVWEta(d,result);
            if (status) { std::copy_n(values.data(),2*R,out.data());
                if (hydrostatic) std::copy_n(values.data()+3*R,R,out.data()+2*R);
                else std::copy_n(values.data()+2*R,2*R,out.data()+2*R); }
            return status;
        });
    WVPortableTypedRecord customConfiguration{"wave-vortex-forcing-configuration-v1",1,
        {{"rate",{},std::vector<double>{-.01}}}};
    WVFrozenForcingSchedule custom;
    custom.entries={{wavevortex::runtime::test::LinearCoefficientForcingIdentifier,1,
        "unqualified",WVForcingStage::spectral,90,0,"",customConfiguration}};
    std::unique_ptr<WVConstantStratificationForcingEngine> unsupported;
    require(bool(WVConstantStratificationForcingEngine::create(c,custom,catalog,
        std::make_unique<FailingEngine>(counter),unsupported)),"Custom diagnostic setup");
    std::vector<WVComplex64> coefficients(3*S);
    const WVState state{0,0,{{coefficients.data(),spectral},{coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}}};
    std::vector<double> output(spatial.elementCount(),99);
    const WVForcingTendencyOutput request{0,{output.data(),spatial}};
    const auto calls=counter->calls;
    require(unsupported->evaluateForcingTendencies(state,&request,1).code==WVKernelStatusCode::unsupportedOperation &&
        counter->calls==calls && unsupported->tendencyMetrics().workspaceHighWaterBytes==0,
        "Unqualified extension was not rejected before evaluation and allocation");
    for (auto value:output) require(value==99,"Preflight rejection changed output");
    std::unique_ptr<WVFieldEvaluationService> customService;
    require(bool(WVFieldEvaluationService::createBorrowing(*unsupported,customService)),"Bind custom forcing service");
    WVFieldEvaluationPlan rejected;
    require(customService->createPlan({{"unsupported","Fu_unqualified",{}}},rejected).code==WVKernelStatusCode::unsupportedOperation,
        "Field service accepted an unqualified custom forcing");
    auto collisionSchedule=schedule(S);
    collisionSchedule.entries[2].name="same-name";
    collisionSchedule.entries[3].name="same name";
    std::unique_ptr<WVConstantStratificationForcingEngine> collision;
    require(bool(WVConstantStratificationForcingEngine::create(c,collisionSchedule,catalog,
        std::make_unique<FailingEngine>(counter),collision)),"Colliding forcing names fixture");
    std::unique_ptr<WVFieldEvaluationService> collisionService;
    require(bool(WVFieldEvaluationService::createBorrowing(*collision,collisionService)),"Bind colliding forcing names");
    require(collisionService->createPlan({{"ambiguous","Fu_same_name",{}}},rejected).code==WVKernelStatusCode::unsupportedOperation,
        "Sanitized forcing-name collision resolved arbitrarily");
    require(bool(collisionService->createPlan({{"ordinary","u",{}}},rejected)),"Forcing ambiguity disabled ordinary fields");

    std::unique_ptr<WVFieldEvaluationService> originalFields;
    require(bool(WVFieldEvaluationService::createBorrowing(*engine,originalFields)),"Original rebind fields");
    WVPortableObserverRecord bindingRecord;
    WVObserverRecord bindingObserver;
    bindingObserver.identifier="filter-diagnostic"; bindingObserver.name="filter diagnostic";
    bindingObserver.typeIdentifier="WVEulerianFields"; bindingObserver.fieldNames={"Fu_mode_filter"};
    bindingRecord.observers.push_back(bindingObserver);
    WVPortableObserverDescriptor bindingDescriptor;
    require(bool(WVPortableObserverDescriptor::create(bindingRecord,catalog,bindingDescriptor)),"Rebind descriptor");
    std::unique_ptr<WVObserverOutputEvaluationService> boundObserver;
    require(bool(WVObserverOutputEvaluationService::create(false,bindingDescriptor,*originalFields,boundObserver)),"Rebind observer");
    for(int change=0;change<3;++change) {
        auto replacement=schedule(S);
        if(change==0) {
            // Preserve the requested filter and every graph name/ordinal,
            // but replace an earlier contribution with a different class.
            replacement.entries.back().typeIdentifier="WVBottomFrictionLinear";
            replacement.entries.back().configuration=replacement.entries[2].configuration;
        } else if(change==1) {
            replacement.entries.back().priority=126;
        } else {
            replacement.entries.erase(replacement.entries.begin()+2);
        }
        std::unique_ptr<WVConstantStratificationForcingEngine> replacementEngine;
        require(bool(WVConstantStratificationForcingEngine::create(c,replacement,catalog,
            std::make_unique<FailingEngine>(counter),replacementEngine)),"Replacement forcing schedule");
        std::unique_ptr<WVFieldEvaluationService> replacementFields;
        require(bool(WVFieldEvaluationService::createBorrowing(*replacementEngine,replacementFields)),"Replacement fields");
        const auto callsBefore=counter->calls;
        require(!boundObserver->useFieldEvaluationService(*replacementFields) && counter->calls==callsBefore,
            "Observer rebind accepted a changed forcing class, priority, or ordered prefix");
        require(bool(boundObserver->useFieldEvaluationService(*originalFields)),"Rejected rebind changed the original observer plan");
    }

    WVFrozenForcingSchedule filterOnly; filterOnly.entries={schedule(S).entries[1]};
    std::unique_ptr<WVConstantStratificationForcingEngine> filterEngine;
    require(bool(WVConstantStratificationForcingEngine::create(c,filterOnly,catalog,
        std::make_unique<FailingEngine>(counter),filterEngine)),"Spectral-only forcing fixture");
    std::unique_ptr<WVFieldEvaluationService> filterService;
    require(bool(WVFieldEvaluationService::createBorrowing(*filterEngine,filterService)) &&
        bool(filterService->createPlan({{"filter","Fu_mode_filter",{}}},rejected)),"Spectral-only forcing field plan");
    WVFieldOutputView filtered{output.data(),R};
    require(bool(filterService->evaluate(rejected,state,&filtered,1)),"Spectral-only forcing field evaluation");
    require(filterEngine->metrics().physicalFieldReconstructionCount==0 &&
        filterService->metrics().primitiveFieldEvaluationCount==0 && filterService->metrics().diagnosticIntermediateReuseCount==0,
        "Spectral-only diagnostic unnecessarily reconstructed physical state");
    for(std::size_t index=0;index<R;++index) require(output[index]==0,"Filtering an empty accumulator produced forcing");

}

void builtinNonlinearCoefficientBoundary() {
    WVTransformConstantStratificationConfiguration c;
    c.Nx=5; c.Ny=4; c.Nz=6; c.Nj=4; c.Lx=9000; c.Ly=8000; c.Lz=900;
    c.N0=5e-3; c.rho0=1027; c.g=9.80665; c.planetaryRadius=6.3712e6;
    c.rotationRate=7.292115e-5; c.latitude=31; c.isHydrostatic=true;
    auto counter=std::make_shared<FailureCounter>();
    std::unique_ptr<WVConstantStratificationForcingEngine> engine;
    require(bool(WVConstantStratificationForcingEngine::create(c,
        defaultNonlinearAdvectionSchedule(),wavevortex::runtime::test::extensionCatalog(),
        std::make_unique<FailingEngine>(counter),engine)),
        "Built-in coefficient engine setup");
    std::unique_ptr<WVFieldEvaluationService> service;
    require(bool(WVFieldEvaluationService::createBorrowing(*engine,service)) &&
        bool(service->prepareBuiltinNonlinearCoefficientEvaluation()),
        "Built-in coefficient service preparation");
    WVIntegrationStateLayout layout;
    require(bool(service->createStateLayout({},layout)),"Built-in coefficient layout");
    const auto spectral=engine->stateShape(); const auto S=spectral.elementCount();
    std::vector<WVComplex64> coefficients(3*S);
    for(std::size_t i=0;i<coefficients.size();++i)
        coefficients[i]={2e-4*std::sin(.19*(i+1)),2e-4*std::cos(.11*(i+1))};
    WVState state{7,2,{{coefficients.data(),spectral},{coefficients.data()+S,spectral},
        {coefficients.data()+2*S,spectral}}};
    std::array<WVCoefficientFamilyConstView,3> families;
    for(std::size_t i=0;i<3;++i)
        families[i]={&layout.coefficientFamilies()[i],coefficients.data()+i*S};
    WVIntegrationState integration{state,nullptr,0,families.data(),families.size()};
    std::vector<WVComplex64> values(3*S,{91,92});
    WVFlux flux{{values.data(),spectral},{values.data()+S,spectral},
        {values.data()+2*S,spectral}};
    require(service->evaluateBuiltinNonlinearCoefficients(integration,flux).code==
        WVKernelStatusCode::invalidConfiguration,
        "Built-in coefficient evaluation accepted no active session");
    require(bool(service->setVariableEvaluationPolicy(WVVariableEvaluationPolicy::lowMemory)),
        "Built-in coefficient low-memory policy");
    {
        WVFieldEvaluationSession session;
        require(bool(service->beginEvaluationSession(integration,session)) &&
            service->evaluateBuiltinNonlinearCoefficients(integration,flux).code==
                WVKernelStatusCode::invalidConfiguration,
            "Built-in coefficient evaluation accepted low-memory policy");
    }
    require(bool(service->setVariableEvaluationPolicy(WVVariableEvaluationPolicy::reuse)) &&
        bool(service->prepareBuiltinNonlinearCoefficientEvaluation()),
        "Built-in coefficient reuse restore");
    const auto sentinel=values;
    std::vector<WVComplex64> coefficientReference;
    counter->failAt=counter->calls+1;
    {
        WVFieldEvaluationSession session;
        require(bool(service->beginEvaluationSession(integration,session)),
            "Built-in coefficient failure session");
        require(service->evaluateBuiltinNonlinearCoefficients(integration,flux).code==
            WVKernelStatusCode::fftExecutionFailure && equal(values,sentinel),
            "Built-in coefficient failure exposed partial output");
        counter->failAt=0;
        require(bool(service->evaluateBuiltinNonlinearCoefficients(integration,flux)),
            "Built-in coefficient retry");
        coefficientReference=values;
        require(bool(service->evaluateBuiltinNonlinearCoefficients(integration,flux)) &&
            equal(values,coefficientReference),"Built-in coefficient repeat changed cached result");
        auto overlapping=flux; overlapping.Fp.data=coefficients.data();
        require(service->evaluateBuiltinNonlinearCoefficients(integration,overlapping).code==
            WVKernelStatusCode::overlappingArrays,
            "Built-in coefficient output accepted immutable-state overlap");
        std::vector<std::uint8_t> misalignedBytes(S*sizeof(WVComplex64)+1);
        auto misaligned=flux;
        misaligned.Fp.data=reinterpret_cast<WVComplex64*>(misalignedBytes.data()+1);
        require(service->evaluateBuiltinNonlinearCoefficients(integration,misaligned).code==
            WVKernelStatusCode::invalidPointer,
            "Built-in coefficient output accepted a misaligned address");
        auto forged=integration;
        forged.waveVortex.coefficients.Ap.shape.rows++;
        require(service->evaluateBuiltinNonlinearCoefficients(forged,flux).code==
            WVKernelStatusCode::invalidShape,
            "Built-in coefficient evaluation trusted a forged caller shape");
    }
    WVFieldEvaluationPlan rawPlan;
    require(bool(service->createPlan({{"builtin-raw","Fu_nonlinear_advection",{}}},rawPlan)),
        "Built-in coefficient raw plan");
    std::vector<double> raw(rawPlan.outputs()[0].elementCount);
    WVFieldOutputView rawView{raw.data(),raw.size()};
    const auto producersBefore=engine->metrics().nonlinearProducerCount;
    {
        WVFieldEvaluationSession session;
        require(bool(service->beginEvaluationSession(integration,session)) &&
            bool(service->evaluate(rawPlan,integration,&rawView,1)) &&
            bool(service->evaluateBuiltinNonlinearCoefficients(integration,flux)),
            "Raw-first coefficient evaluation changed values");
        double scale=0,error=0;
        for(std::size_t i=0;i<values.size();++i) {
            scale=std::max({scale,std::abs(coefficientReference[i].real),
                std::abs(coefficientReference[i].imag)});
            error=std::max({error,std::abs(values[i].real-coefficientReference[i].real),
                std::abs(values[i].imag-coefficientReference[i].imag)});
        }
        require(error<=1e-12*std::max(scale,1e-30),
            "Raw-first coefficient projection differs from fused projection");
        if(engine->metrics().nonlinearProducerCount!=producersBefore+1)
            throw std::runtime_error("Raw-first nonlinear producers: before="+
                std::to_string(producersBefore)+" after="+
                std::to_string(engine->metrics().nonlinearProducerCount));
    }
    auto foreignCoefficients=coefficients;
    auto foreign=integration;
    foreign.waveVortex.coefficients.Ap.data=foreignCoefficients.data();
    auto foreignFamilies=families;
    for(std::size_t i=0;i<3;++i)
        foreignFamilies[i].data=foreignCoefficients.data()+i*S;
    foreign.coefficientFamilies=foreignFamilies.data();
    {
        WVFieldEvaluationSession session;
        require(bool(service->beginEvaluationSession(integration,session)) &&
            !service->evaluateBuiltinNonlinearCoefficients(foreign,flux),
            "Built-in coefficient evaluation accepted a foreign state");
    }
}

void barotropic() {
    WVTransformBarotropicQGConfiguration c;
    c.Nx=8; c.Ny=6; c.Lx=17000; c.Ly=11000; c.h=.8; c.j=1;
    c.g=9.81; c.planetaryRadius=6.371e6; c.rotationRate=7.2921e-5; c.latitude=33; c.shouldAntialias=false;
    auto catalog=wavevortex::runtime::test::extensionCatalog();
    auto counter=std::make_shared<FailureCounter>();
    std::unique_ptr<WVBarotropicQGForcingEngine> engine;
    require(bool(WVBarotropicQGForcingEngine::create(c,{},catalog,std::make_unique<FailingEngine>(counter),engine)),"QG setup");
    const auto spectral=engine->kernel().descriptor().spectralShape(),plane=engine->kernel().descriptor().spatialShape();
    const auto S=spectral.elementCount(),R=plane.elementCount();
    const WVShape4D spatial{plane.rows,plane.columns,1,1};
    auto scheduleValue=schedule(S);
    auto& entries=scheduleValue.entries;
    entries[0].configuration.values.erase(entries[0].configuration.values.begin(),entries[0].configuration.values.begin()+6);
    entries.push_back({"WVBottomFrictionQuadratic",1,"quadratic drag",WVForcingStage::spatial,255,4,"",
        {"wave-vortex-forcing-configuration-v1",1,{{"Cd",{},std::vector<double>{.002}}}}});
    const WVPortableTypedRecord empty{"wave-vortex-forcing-configuration-v1",1,{}};
    entries.push_back({"WVBetaPlanePVAdvection",1,"beta advection",WVForcingStage::spatial,255,5,"",empty});
    entries.push_back({"WVAdaptiveDamping",1,"adaptive damping",WVForcingStage::spectral,255,6,"",empty});
    require(bool(WVBarotropicQGForcingEngine::create(c,scheduleValue,catalog,std::make_unique<FailingEngine>(counter),engine)),"QG diagnostic schedule");
    std::vector<WVComplex64> a(S),f(S);
    for (std::size_t i=0;i<S;++i) a[i]={1e-6*std::sin(.17*i),1e-6*std::cos(.23*i)};
    WVComplexView mutableA{a.data(),spectral}; engine->kernel().enforceReality(mutableA);
    const WVComplexConstView state{a.data(),spectral}; WVComplexView flux{f.data(),spectral};
    require(bool(engine->evaluateRightHandSide(state,flux)),"QG baseline RHS");
    const auto before=a,rhs=f; const auto bytes=engine->persistentBytes();
    const auto calls=engine->metrics().evaluationCount,reconstructions=engine->metrics().physicalFieldReconstructionCount;
    std::vector<std::vector<double>> values(7,std::vector<double>(R,99));
    std::vector<WVForcingTendencyOutput> outputs;
    for (std::size_t i=0;i<values.size();++i) outputs.push_back({i,{values[i].data(),spatial}});
    const auto start=counter->calls;
    auto status=engine->evaluateForcingTendencies(state,outputs.data(),outputs.size());
    if (!status) throw std::runtime_error(status.message);
    const auto fftCalls=counter->calls-start;
    require(engine->forcingCount()==7 && engine->forcingInstance(0)->name()=="nonlinear advection" &&
        engine->forcingInstance(6)->name()=="held amplitudes" && !engine->forcingInstance(7),"QG instance binding/order");
    require(engine->tendencyMetrics().forcingEvaluationCount==7 && engine->tendencyMetrics().spatialProjectionCount==1 &&
        engine->tendencyMetrics().spectralReconstructionCount==3 && engine->tendencyMetrics().workspaceLiveBytes==0 &&
        engine->persistentBytes()==bytes && engine->metrics().evaluationCount==calls && equal(a,before),"QG diagnostic state/metrics");
    require(engine->metrics().physicalFieldReconstructionCount==reconstructions+1,"QG same-event fields were reconstructed twice");
    // Independent nonlinear product, without projecting out unretained modes.
    std::vector<double> u(R),v(R),q(3*R),expected(R);
    WVRealView uv{u.data(),plane},vv{v.data(),plane};
    require(bool(engine->kernel().transformA0ToField(state,WVBarotropicQGField::u,uv)) &&
        bool(engine->kernel().transformA0ToField(state,WVBarotropicQGField::v,vv)),"QG independent velocity");
    WVRealFieldBundleView derivatives{q.data(),{plane.rows,plane.columns,1,3}};
    require(bool(engine->kernel().transformA0ToFieldWithDerivatives(state,WVBarotropicQGField::qgpv,derivatives)),"QG independent derivatives");
    for (std::size_t i=0;i<R;++i) expected[i]=-u[i]*q[R+i]-v[i]*q[2*R+i];
    relative(values[0],expected);
    std::vector<double> sum(R);
    for (std::size_t n=0;n<6;++n) for (std::size_t i=0;i<R;++i) sum[i]+=values[n][i];
    // Spatial raw sums need projection before comparison with fixed amplitudes.
    require(bool(engine->kernel().transformQGPVToA0({sum.data(),plane},flux)),"QG cumulative projection");
    WVRealView expectedView{expected.data(),plane};
    require(bool(engine->kernel().transformSpectralTendencyToSpatial({f.data(),spectral},expectedView)),"QG cumulative inverse");
    for (auto& x:expected) x=-x;
    relative(values[6],expected);
    require(bool(engine->evaluateRightHandSide(state,flux)) && equal(f,rhs) && equal(a,before),"QG diagnostics changed later RHS");
    const auto successful=values;
    qgAdaptiveMaximumEvaluationSession(*engine,state);
    {
        WVVariableEvaluationContext evaluation;
        require(bool(evaluation.prepare(
                    detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                        engine->forcingCount()))) &&
                bool(evaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::reuse)),
            "QG persistent diagnostic evaluation setup");
        detail::WVForcingDiagnosticWorkspace session(spectral,spatial,1,2);
        std::vector<WVForcingStage> stages;
        for(std::size_t index=0;index<engine->forcingCount();++index)
            stages.push_back(engine->forcingInstance(index)->stage());
        require(bool(session.beginScopedEvaluation(evaluation,stages)),
            "QG persistent diagnostic workspace setup");
        std::vector<double> first(R),reusedFirst(R),last(R);
        const WVForcingTendencyOutput firstRequest{0,{first.data(),spatial}};
        const WVForcingTendencyOutput lastRequest{6,{last.data(),spatial}};
        const WVForcingTendencyOutput orderedRequests[] = {
            {0,{reusedFirst.data(),spatial}},lastRequest};
        const auto forcingBefore=engine->tendencyMetrics().forcingEvaluationCount;
        require(bool(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&session)) &&
                first==successful[0] && session.initialized() &&
                session.physicalPrepared,
            "QG persistent diagnostic first prefix");
        require(bool(engine->evaluateForcingTendencies(
                    state,orderedRequests,2,nullptr,&session)) &&
                reusedFirst==successful[0] && last==successful[6] &&
                engine->tendencyMetrics().forcingEvaluationCount-forcingBefore==7,
            "QG persistent diagnostic extended prefix replayed prior forcing");
        auto foreign=a;
        const WVComplexConstView foreignState{foreign.data(),spectral};
        std::fill(last.begin(),last.end(),99);
        require(engine->evaluateForcingTendencies(
                    foreignState,&lastRequest,1,nullptr,&session).code==
                    WVKernelStatusCode::invalidConfiguration &&
                std::all_of(last.begin(),last.end(),[](double value) {
                    return value==99;
                }),
            "QG persistent diagnostic accepted a foreign immutable state");
        evaluation.end();
        detail::WVForcingDiagnosticWorkspace incompatible(
            spectral,spatial,1,1);
        std::fill(first.begin(),first.end(),99);
        require(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&incompatible).code==
                    WVKernelStatusCode::invalidShape &&
                std::all_of(first.begin(),first.end(),[](double value) {
                    return value==99;
                }),
            "QG persistent diagnostic accepted incompatible channel storage");
        WVVariableEvaluationContext lowEvaluation;
        require(bool(lowEvaluation.prepare(
                    detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                        engine->forcingCount()))) &&
                bool(lowEvaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::lowMemory)),
            "QG low-memory diagnostic evaluation setup");
        detail::WVForcingDiagnosticWorkspace lowSession(
            spectral,spatial,1,2);
        require(bool(lowSession.beginScopedEvaluation(lowEvaluation,stages)),
            "QG low-memory diagnostic workspace setup");
        std::fill(first.begin(),first.end(),0);
        std::fill(reusedFirst.begin(),reusedFirst.end(),0);
        std::fill(last.begin(),last.end(),0);
        const auto lowForcingBefore=
            engine->tendencyMetrics().forcingEvaluationCount;
        const auto lowReductionBefore=
            engine->kernel().metrics().horizontalSpeedMaximumReductionCount;
        const auto engineBytes=engine->persistentBytes();
        require(bool(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&lowSession)) &&
                lowSession.physicalPrepared,
            "QG low-memory diagnostic did not prepare its persistent velocity");
        const auto lowPhysical=lowSession.physical;
        require(bool(engine->evaluateForcingTendencies(
                    state,orderedRequests,2,nullptr,&lowSession)) &&
                first==successful[0] && reusedFirst==successful[0] &&
                last==successful[6] && lowSession.physical==lowPhysical,
            "QG low-memory incremental diagnostics changed ordered output");
        const auto lowMetrics=lowEvaluation.metrics();
        const auto lowBytes=lowSession.bytes();
        require(engine->tendencyMetrics().forcingEvaluationCount-
                    lowForcingBefore==8 &&
                lowMetrics.producerExecutions==11 &&
                lowMetrics.recomputations==1 &&
                lowMetrics.evictions==10 && lowMetrics.cacheHits==0 &&
                engine->kernel().metrics().horizontalSpeedMaximumReductionCount==
                    lowReductionBefore+1 &&
                lowMetrics.liveBytes==sizeof(double) && lowSession.prefix.empty() &&
                lowSession.nextIndex==0 && !lowSession.projected &&
                lowSession.physicalPrepared &&
                engine->persistentBytes()==engineBytes,
            "QG low-memory prefix replay was retained or not fully counted");
        lowEvaluation.end();
        require(lowEvaluation.metrics().liveBytes==0,
            "Closed QG low-memory context retained its scalar reduction");
        std::fill(last.begin(),last.end(),99);
        require(engine->evaluateForcingTendencies(
                    state,&lastRequest,1,nullptr,&lowSession).code==
                    WVKernelStatusCode::invalidConfiguration &&
                std::all_of(last.begin(),last.end(),[](double value) {
                    return value==99;
                }),
            "QG closed diagnostic context retained an active view");
        require(bool(lowEvaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::lowMemory)),
            "QG low-memory diagnostic context reopen");
        detail::WVForcingDiagnosticWorkspace reopened(
            spectral,spatial,1,2);
        require(bool(reopened.beginScopedEvaluation(lowEvaluation,stages)) &&
                bool(engine->evaluateForcingTendencies(
                    state,&lastRequest,1,nullptr,&reopened)) &&
                last==successful[6] && reopened.bytes()==lowBytes &&
                reopened.prefix.empty() &&
                lowEvaluation.metrics().liveBytes==sizeof(double) &&
                engine->persistentBytes()==engineBytes,
            "QG reopened low-memory context grew or retained a stale prefix");
        lowEvaluation.end();
        require(lowEvaluation.metrics().liveBytes==0,
            "Reclosed QG low-memory context retained its scalar reduction");
    }
    for (auto& field:values) std::fill(field.begin(),field.end(),99);
    counter->failAt=counter->calls+fftCalls;
    require(engine->evaluateForcingTendencies(state,outputs.data(),outputs.size()).code==WVKernelStatusCode::fftExecutionFailure,"QG late FFT injection");
    for (const auto& field:values) for (auto x:field) require(x==99,"QG partial diagnostic escaped failed evaluation");
    require(equal(a,before) && engine->persistentBytes()==bytes && engine->tendencyMetrics().workspaceLiveBytes==0,"QG failed diagnostic changed state/storage");
    counter->failAt=0;
    require(bool(engine->evaluateForcingTendencies(state,outputs.data(),outputs.size())) && values==successful,"QG retry differs");
    std::vector<double> selected(R); const WVForcingTendencyOutput request{6,{selected.data(),spatial}};
    require(bool(engine->evaluateForcingTendencies(state,&request,1)) && selected==successful[6],"QG selected-only prefix lost previous stages");
    const auto beforeInvalid=counter->calls;
    const auto saved=a[0]; a[0].real=std::numeric_limits<double>::quiet_NaN();
    std::fill(selected.begin(),selected.end(),99);
    require(engine->evaluateForcingTendencies(state,&request,1).code==WVKernelStatusCode::invalidConfiguration &&
        counter->calls==beforeInvalid,"Nonfinite QG diagnostic state accepted");
    for (auto x:selected) require(x==99,"Invalid QG state changed output");
    a[0]=saved;
    engine->setLinearDynamics(true);
    const auto forcingCalls=engine->metrics().forcingCallCount;
    require(bool(engine->evaluateRightHandSide(state,flux)) &&
        equal(f,std::vector<WVComplex64>(f.size())) && equal(a,before) &&
        engine->metrics().forcingCallCount==forcingCalls,
        "Linear QG RHS executed forcing or changed coefficients");
    auto constrained=a;
    WVComplexView constraintView{constrained.data(),spectral};
    const auto constraint=engine->restoreForcingAmplitudes(constraintView);
    require(bool(constraint) && constraint.modifiedCoefficientCount>0 && !equal(constrained,a),
        "Linear QG evolution discarded fixed-amplitude constraints");
    fieldService(*engine,WVState{0,0,{{},{},state}},spatial,successful,counter);
    engine->setLinearDynamics(false);
    require(bool(engine->evaluateRightHandSide(state,flux)) && equal(f,rhs) && equal(a,before),
        "QG field-service diagnostics changed subsequent RHS or scientific state");
    // Raw inverse must preserve the mean, while ordinary qgpv retains its mask.
    std::fill(f.begin(),f.end(),WVComplex64{});
    for (std::size_t i=0;i<S;++i) if (engine->kernel().descriptor().fourierModes()[i].Kh==0) f[i]={.25,7};
    require(bool(engine->kernel().transformSpectralTendencyToSpatial({f.data(),spectral},expectedView)),"QG raw mean inverse");
    for (auto x:expected) require(std::abs(x-.25)<1e-14,"QG diagnostic mean was masked or imaginary part escaped");
    require(bool(engine->kernel().transformA0ToQGPV({f.data(),spectral},expectedView)),"QG ordinary mean inverse");
    for (auto x:expected) require(x==0,"QG ordinary mean mask changed");
    WVBarotropicQGOperationWorkspace workspace;
    WVRealView raw{reinterpret_cast<double*>(a.data()),plane}; workspace.spatialTendency=&raw;
    const auto beforeReject=counter->calls;
    require(engine->kernel().addPotentialVorticityAdvection(state,flux,false,workspace).code==WVKernelStatusCode::overlappingArrays &&
        counter->calls==beforeReject,"QG raw overlap accepted or rejected after FFT");
}

void stratifiedQG() {
    Temporary file; fixture(file.path); change(file.path,"shouldAntialias",0);
    std::shared_ptr<const WVStratifiedModalRecord> source;
    require(bool(WVStratifiedModalReader::read(file.path.string(),source)),"QG diagnostic modal source");
    const auto& g=source->geometry();
    const WVShape2D spectral{g.Nj,g.Nkl}; const WVShape3D volume{g.Nx,g.Ny,g.Nz};
    const WVShape4D spatial{g.Nx,g.Ny,g.Nz,1}; const auto S=spectral.elementCount(),R=volume.elementCount();
    auto scheduleValue=schedule(S);
    auto& entries=scheduleValue.entries;
    entries[0].configuration.values.erase(entries[0].configuration.values.begin(),entries[0].configuration.values.begin()+6);
    entries.push_back({"WVBottomFrictionQuadratic",1,"quadratic drag",WVForcingStage::spatial,255,4,"",
        {"wave-vortex-forcing-configuration-v1",1,{{"Cd",{},std::vector<double>{.002}}}}});
    entries.push_back({"WVVerticalDiffusivity",1,"vertical mixing",WVForcingStage::spatial,255,5,"",
        {"wave-vortex-forcing-configuration-v1",1,{{"kappa_z",{},std::vector<double>{.002}},
            {"shouldForceMeanDensityAnomaly",{},std::vector<std::uint8_t>{0}}}}});
    const WVPortableTypedRecord empty{"wave-vortex-forcing-configuration-v1",1,{}};
    entries.push_back({"WVBetaPlanePVAdvection",1,"beta advection",WVForcingStage::spatial,255,6,"",empty});
    entries.push_back({"WVAdaptiveDamping",1,"adaptive damping",WVForcingStage::spectral,255,7,"",empty});
    auto counter=std::make_shared<FailureCounter>();
    std::unique_ptr<WVStratifiedQGForcingEngine> engine;
    auto status=WVStratifiedQGForcingEngine::create(source,scheduleValue,wavevortex::runtime::test::extensionCatalog(),
        std::make_unique<FailingEngine>(counter),engine);
    if (!status) throw std::runtime_error("Stratified QG schedule: "+status.message);
    std::vector<WVComplex64> a(S),f(S);
    for (std::size_t i=0;i<S;++i) if (g.k[i/g.Nj]!=0 || g.l[i/g.Nj]!=0)
        a[i]={1e-6*std::sin(.17*i),1e-6*std::cos(.23*i)};
    const WVComplexConstView state{a.data(),spectral}; WVComplexView flux{f.data(),spectral};
    require(bool(engine->evaluateRightHandSide(state,flux)),"Stratified QG baseline RHS");
    const auto& kernelMetrics=engine->kernel().metrics();
    require(kernelMetrics.stateValidationCount==1 &&
        kernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::u)][0]==1 &&
        kernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::v)][0]==1 &&
        kernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::qgpv)][1]==1 &&
        kernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::qgpv)][2]==1,
        "Stratified QG repeated validation or shared field producers");
    require(engine->variableEvaluationMetrics().contexts==1 &&
        engine->variableEvaluationMetrics().producerExecutions==4 &&
        engine->variableEvaluationMetrics().cacheHits==4,
        "Stratified QG default evaluation did not reuse velocity and speed");
    {
      auto lowCounter=std::make_shared<FailureCounter>();
      std::unique_ptr<WVStratifiedQGForcingEngine> lowMemory;
      status=WVStratifiedQGForcingEngine::create(source,scheduleValue,
          wavevortex::runtime::test::extensionCatalog(),
          std::make_unique<FailingEngine>(lowCounter),lowMemory);
      require(bool(status) && bool(lowMemory->setVariableEvaluationPolicy(
          WVVariableEvaluationPolicy::lowMemory)),"Stratified QG low-memory setup");
      std::vector<WVComplex64> lowFlux(S);
      WVComplexView lowFluxView{lowFlux.data(),spectral};
      require(bool(lowMemory->evaluateRightHandSide(state,lowFluxView)) &&
          equal(lowFlux,f),"Stratified QG policies changed RHS arithmetic");
      const auto& lowKernelMetrics=lowMemory->kernel().metrics();
      require(lowKernelMetrics.stateValidationCount==1 &&
          lowKernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::u)][0]==3 &&
          lowKernelMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::v)][0]==3 &&
          lowMemory->variableEvaluationMetrics().producerExecutions==7 &&
          lowMemory->variableEvaluationMetrics().recomputations==4 &&
          lowMemory->variableEvaluationMetrics().evictions==7 &&
          lowMemory->metrics().workspaceCapacityBytes+S*sizeof(WVComplex64)==
              engine->metrics().workspaceCapacityBytes,
          "Stratified QG low-memory policy did not release and recompute velocity");
    }
    {
      auto singleSchedule=scheduleValue;
      const auto nonlinearEntry=std::find_if(
          singleSchedule.entries.begin(),singleSchedule.entries.end(),
          [](const auto& entry) {
            return entry.typeIdentifier=="WVNonlinearAdvection";
          });
      require(nonlinearEntry!=singleSchedule.entries.end(),
          "Stratified QG nonlinear fixture missing");
      singleSchedule.entries={*nonlinearEntry};
      std::vector<WVComplex64> singleFlux(S);
      {
        auto singleCounter=std::make_shared<FailureCounter>();
        std::unique_ptr<WVStratifiedQGForcingEngine> single;
        status=WVStratifiedQGForcingEngine::create(source,singleSchedule,
            wavevortex::runtime::test::extensionCatalog(),
            std::make_unique<FailingEngine>(singleCounter),single);
        WVComplexView singleView{singleFlux.data(),spectral};
        require(bool(status) && bool(single->evaluateRightHandSide(state,singleView)),
            "Single Stratified QG nonlinear reference");
      }
      auto repeatedSchedule=singleSchedule;
      auto repeatedEntry=repeatedSchedule.entries.front();
      repeatedEntry.name="second nonlinear advection";
      repeatedEntry.ordinal=1;
      repeatedSchedule.entries.push_back(std::move(repeatedEntry));
      auto repeatedCounter=std::make_shared<FailureCounter>();
      std::unique_ptr<WVStratifiedQGForcingEngine> repeated;
      status=WVStratifiedQGForcingEngine::create(source,repeatedSchedule,
          wavevortex::runtime::test::extensionCatalog(),
          std::make_unique<FailingEngine>(repeatedCounter),repeated);
      std::vector<WVComplex64> repeatedFlux(S);
      WVComplexView repeatedView{repeatedFlux.data(),spectral};
      require(bool(status) && bool(repeated->evaluateRightHandSide(state,repeatedView)),
          "Repeated Stratified QG nonlinear evaluation");
      for (std::size_t i=0;i<S;++i) {
        require(repeatedFlux[i].real==2*singleFlux[i].real &&
            repeatedFlux[i].imag==2*singleFlux[i].imag,
            "Repeated Stratified QG nonlinear arithmetic changed");
      }
      const auto& repeatedMetrics=repeated->kernel().metrics();
      require(repeatedMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::qgpv)][1]==1 &&
          repeatedMetrics.reconstructionCount[static_cast<std::size_t>(WVStratifiedQGField::qgpv)][2]==1,
          "Repeated Stratified QG nonlinear forcing rebuilt PV derivatives");
      require(repeated->variableEvaluationMetrics().producerExecutions==3,
          "Repeated Stratified QG nonlinear forcing ran an unexpected producer count");
      require(repeated->variableEvaluationMetrics().cacheHits==1,
          "Repeated Stratified QG nonlinear forcing missed its cached tendency");
      auto foreign=a;
      WVComplexConstView foreignState{foreign.data(),spectral};
      WVComplexView mutableState{a.data(),spectral};
      require(bool(repeated->beginStateEvaluation(state)),
          "Stratified QG state evaluation begin failed");
      require(repeated->setVariableEvaluationPolicy(
              WVVariableEvaluationPolicy::lowMemory).code==
              WVKernelStatusCode::reentrantExecution &&
          repeated->restoreForcingAmplitudes(mutableState).status.code==
              WVKernelStatusCode::reentrantExecution &&
          repeated->evaluateRightHandSide(foreignState,repeatedView).code==
              WVKernelStatusCode::invalidConfiguration,
          "Stratified QG active-state ownership or mutation guard failed");
      require(bool(repeated->endStateEvaluation()),
          "Stratified QG state evaluation end failed");
    }
    const auto before=a,rhs=f; const auto bytes=engine->persistentBytes(),calls=engine->metrics().evaluationCount;
    const auto reconstructions=engine->metrics().physicalFieldReconstructionCount;
    const auto reuse=engine->metrics().physicalFieldReuseCount;
    std::vector<std::vector<double>> values(8,std::vector<double>(R,99));
    std::vector<WVForcingTendencyOutput> outputs;
    for (std::size_t i=0;i<values.size();++i) outputs.push_back({i,{values[i].data(),spatial}});
    const auto start=counter->calls;
    status=engine->evaluateForcingTendencies(state,outputs.data(),outputs.size());
    if (!status) throw std::runtime_error("Stratified QG diagnostics: "+status.message);
    const auto fftCalls=counter->calls-start;
    require(engine->forcingInstance(0)->name()=="nonlinear advection" && engine->forcingInstance(7)->name()=="held amplitudes" &&
        !engine->forcingInstance(8),"Stratified QG instance binding");
    require(engine->tendencyMetrics().forcingEvaluationCount==8 && engine->tendencyMetrics().spatialProjectionCount==1 &&
        engine->tendencyMetrics().spectralReconstructionCount==3 && engine->tendencyMetrics().workspaceLiveBytes==0 &&
        engine->persistentBytes()==bytes && engine->metrics().evaluationCount==calls && equal(a,before),"Stratified QG diagnostic state/metrics");
    require(engine->metrics().physicalFieldReconstructionCount==reconstructions+5 &&
        engine->metrics().physicalFieldReuseCount==reuse+3,"Stratified QG did not share velocity fields across stages");
    std::vector<double> u(R),v(R),qx(R),qy(R),expected(R),sum(R);
    require(bool(engine->kernel().transformA0ToField(state,WVStratifiedQGField::u,{u.data(),volume})) &&
        bool(engine->kernel().transformA0ToField(state,WVStratifiedQGField::v,{v.data(),volume})) &&
        bool(engine->kernel().transformA0ToField(state,WVStratifiedQGField::qgpv,{qx.data(),volume},WVStratifiedQGDerivative::x)) &&
        bool(engine->kernel().transformA0ToField(state,WVStratifiedQGField::qgpv,{qy.data(),volume},WVStratifiedQGDerivative::y)),"Independent QG derivative fields");
    for (std::size_t i=0;i<R;++i) expected[i]=-(u[i]*qx[i]+v[i]*qy[i]);
    relative(values[0],expected,"Stratified QG nonlinear");
    for (std::size_t n=0;n<5;++n) for (std::size_t i=0;i<R;++i) sum[i]+=values[n][i];
    require(bool(engine->kernel().transformQGPVToA0({sum.data(),volume},flux)),"QG cumulative modal projection");
    const double cutoff=(2.0/3.0)*2*std::acos(-1.0)*(g.Nx/2)/g.Lx;
    for (std::size_t i=0;i<S;++i)
        if (i%g.Nj!=0 || std::hypot(g.k[i/g.Nj],g.l[i/g.Nj])>cutoff) f[i]={};
    // The modal fixture intentionally has non-inverse projection matrices.
    // Build the expected accumulator in coefficient space; inverse-transforming
    // a spectral delta and projecting it again is not an identity here.
    WVFrozenForcingSchedule spectralSchedule;
    spectralSchedule.entries={entries[1],entries.back()};
    std::unique_ptr<WVStratifiedQGForcingEngine> spectralEngine;
    require(bool(WVStratifiedQGForcingEngine::create(source,spectralSchedule,wavevortex::runtime::test::extensionCatalog(),
        std::make_unique<WVReferenceFFTEngine>(),spectralEngine)),"Independent spectral schedule");
    std::vector<WVComplex64> spectralContribution(S); WVComplexView spectralFlux{spectralContribution.data(),spectral};
    require(bool(spectralEngine->evaluateRightHandSide(state,spectralFlux)),"Independent adaptive contribution");
    for (std::size_t i=0;i<S;++i) { f[i].real+=spectralContribution[i].real; f[i].imag+=spectralContribution[i].imag; }
    require(bool(engine->kernel().transformSpectralTendencyToSpatial({f.data(),spectral},{expected.data(),volume})),"QG cumulative modal inverse");
    for (auto& x:expected) x=-x;
    relative(values[7],expected,"Stratified QG cumulative fixed");
    // Raw bottom drag stays at the bottom cell; inverse-projecting its modal
    // flux would instead spread the profile through the truncated basis.
    for (std::size_t i=g.Nx*g.Ny;i<R;++i) require(values[1][i]==0 && values[2][i]==0,"QG raw bottom drag leaked vertically");
    require(bool(engine->evaluateRightHandSide(state,flux)) && equal(f,rhs) && equal(a,before),"Stratified QG later RHS changed");
    const auto successful=values;
    qgAdaptiveMaximumEvaluationSession(*engine,state);
    {
        WVVariableEvaluationContext evaluation;
        require(bool(evaluation.prepare(
                    detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                        engine->forcingCount()))) &&
                bool(evaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::reuse)),
            "Stratified QG persistent diagnostic evaluation setup");
        detail::WVForcingDiagnosticWorkspace session(spectral,spatial,1,2);
        std::vector<WVForcingStage> stages;
        for(std::size_t index=0;index<engine->forcingCount();++index)
            stages.push_back(engine->forcingInstance(index)->stage());
        require(bool(session.beginScopedEvaluation(evaluation,stages)),
            "Stratified QG persistent diagnostic workspace setup");
        std::vector<double> first(R),reusedFirst(R),last(R);
        const WVForcingTendencyOutput firstRequest{0,{first.data(),spatial}};
        const WVForcingTendencyOutput lastRequest{7,{last.data(),spatial}};
        const WVForcingTendencyOutput orderedRequests[] = {
            {0,{reusedFirst.data(),spatial}},lastRequest};
        const auto forcingBefore=engine->tendencyMetrics().forcingEvaluationCount;
        require(bool(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&session)) &&
                first==successful[0] && session.initialized() &&
                session.physicalPrepared,
            "Stratified QG persistent diagnostic first prefix");
        require(bool(engine->evaluateForcingTendencies(
                    state,orderedRequests,2,nullptr,&session)) &&
                reusedFirst==successful[0] && last==successful[7] &&
                engine->tendencyMetrics().forcingEvaluationCount-forcingBefore==8,
            "Stratified QG persistent diagnostic extended prefix replayed prior forcing");
        auto foreign=a;
        const WVComplexConstView foreignState{foreign.data(),spectral};
        std::fill(last.begin(),last.end(),99);
        require(engine->evaluateForcingTendencies(
                    foreignState,&lastRequest,1,nullptr,&session).code==
                    WVKernelStatusCode::invalidConfiguration &&
                std::all_of(last.begin(),last.end(),[](double value) {
                    return value==99;
                }),
            "Stratified QG persistent diagnostic accepted a foreign immutable state");
        evaluation.end();
        detail::WVForcingDiagnosticWorkspace incompatible(
            spectral,spatial,1,1);
        std::fill(first.begin(),first.end(),99);
        require(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&incompatible).code==
                    WVKernelStatusCode::invalidShape &&
                std::all_of(first.begin(),first.end(),[](double value) {
                    return value==99;
                }),
            "Stratified QG persistent diagnostic accepted incompatible channel storage");
        WVVariableEvaluationContext lowEvaluation;
        require(bool(lowEvaluation.prepare(
                    detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                        engine->forcingCount()))) &&
                bool(lowEvaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::lowMemory)),
            "Stratified QG low-memory diagnostic evaluation setup");
        detail::WVForcingDiagnosticWorkspace lowSession(
            spectral,spatial,1,2);
        require(bool(lowSession.beginScopedEvaluation(lowEvaluation,stages)),
            "Stratified QG low-memory diagnostic workspace setup");
        std::fill(first.begin(),first.end(),0);
        std::fill(reusedFirst.begin(),reusedFirst.end(),0);
        std::fill(last.begin(),last.end(),0);
        const auto lowForcingBefore=
            engine->tendencyMetrics().forcingEvaluationCount;
        const auto lowReductionBefore=
            engine->metrics().horizontalSpeedMaximumReductionCount;
        const auto engineBytes=engine->persistentBytes();
        require(bool(engine->evaluateForcingTendencies(
                    state,&firstRequest,1,nullptr,&lowSession)) &&
                lowSession.physicalPrepared,
            "Stratified QG low-memory diagnostic did not prepare its persistent velocity");
        const auto lowPhysical=lowSession.physical;
        require(bool(engine->evaluateForcingTendencies(
                    state,orderedRequests,2,nullptr,&lowSession)) &&
                first==successful[0] && reusedFirst==successful[0] &&
                last==successful[7] && lowSession.physical==lowPhysical,
            "Stratified QG low-memory incremental diagnostics changed ordered output");
        const auto lowMetrics=lowEvaluation.metrics();
        const auto lowBytes=lowSession.bytes();
        require(engine->tendencyMetrics().forcingEvaluationCount-
                    lowForcingBefore==9 &&
                lowMetrics.producerExecutions==12 &&
                lowMetrics.recomputations==1 &&
                lowMetrics.evictions==11 && lowMetrics.cacheHits==0 &&
                engine->metrics().horizontalSpeedMaximumReductionCount==
                    lowReductionBefore+1 &&
                lowMetrics.liveBytes==sizeof(double) && lowSession.prefix.empty() &&
                lowSession.nextIndex==0 && !lowSession.projected &&
                lowSession.physicalPrepared &&
                engine->persistentBytes()==engineBytes,
            "Stratified QG low-memory prefix replay was retained or not fully counted");
        lowEvaluation.end();
        require(lowEvaluation.metrics().liveBytes==0,
            "Closed Stratified QG low-memory context retained its scalar reduction");
        std::fill(last.begin(),last.end(),99);
        require(engine->evaluateForcingTendencies(
                    state,&lastRequest,1,nullptr,&lowSession).code==
                    WVKernelStatusCode::invalidConfiguration &&
                std::all_of(last.begin(),last.end(),[](double value) {
                    return value==99;
                }),
            "Stratified QG closed diagnostic context retained an active view");
        require(bool(lowEvaluation.begin(engine.get(),
                    WVVariableEvaluationPolicy::lowMemory)),
            "Stratified QG low-memory diagnostic context reopen");
        detail::WVForcingDiagnosticWorkspace reopened(
            spectral,spatial,1,2);
        require(bool(reopened.beginScopedEvaluation(lowEvaluation,stages)) &&
                bool(engine->evaluateForcingTendencies(
                    state,&lastRequest,1,nullptr,&reopened)) &&
                last==successful[7] && reopened.bytes()==lowBytes &&
                reopened.prefix.empty() &&
                lowEvaluation.metrics().liveBytes==sizeof(double) &&
                engine->persistentBytes()==engineBytes,
            "Stratified QG reopened low-memory context grew or retained a stale prefix");
        lowEvaluation.end();
        require(lowEvaluation.metrics().liveBytes==0,
            "Reclosed Stratified QG low-memory context retained its scalar reduction");
    }
    for (auto& field:values) std::fill(field.begin(),field.end(),99);
    counter->failAt=counter->calls+fftCalls;
    require(engine->evaluateForcingTendencies(state,outputs.data(),outputs.size()).code==WVKernelStatusCode::fftExecutionFailure,"Stratified QG late failure");
    for (const auto& field:values) for (auto x:field) require(x==99,"Stratified QG partial output escaped failure");
    require(equal(a,before) && engine->persistentBytes()==bytes && engine->tendencyMetrics().workspaceLiveBytes==0,"Stratified QG failure changed state/storage");
    counter->failAt=0;
    require(bool(engine->evaluateForcingTendencies(state,outputs.data(),outputs.size())) && values==successful,"Stratified QG retry mismatch");
    std::vector<double> selected(R); const WVForcingTendencyOutput request{7,{selected.data(),spatial}};
    require(bool(engine->evaluateForcingTendencies(state,&request,1)) && selected==successful[7],"Stratified QG selected-only prefix mismatch");
    engine->setLinearDynamics(true);
    const auto forcingCalls=engine->metrics().forcingCallCount;
    require(bool(engine->evaluateRightHandSide(state,flux)) &&
        equal(f,std::vector<WVComplex64>(f.size())) && equal(a,before) &&
        engine->metrics().forcingCallCount==forcingCalls,
        "Linear QG RHS executed forcing or changed coefficients");
    auto constrained=a;
    WVComplexView constraintView{constrained.data(),spectral};
    const auto constraint=engine->restoreForcingAmplitudes(constraintView);
    require(bool(constraint) && constraint.modifiedCoefficientCount>0 && !equal(constrained,a),
        "Linear QG evolution discarded fixed-amplitude constraints");
    fieldService(*engine,WVState{0,0,{{},{},state}},spatial,successful,counter);
    engine->setLinearDynamics(false);
    require(bool(engine->evaluateRightHandSide(state,flux)) && equal(f,rhs) && equal(a,before),
        "QG field-service diagnostics changed subsequent RHS or scientific state");
    std::fill(f.begin(),f.end(),WVComplex64{});
    for (std::size_t i=0;i<S;++i) if (g.k[i/g.Nj]==0 && g.l[i/g.Nj]==0) f[i]={.25,7};
    const auto meanBefore=f;
    require(bool(engine->kernel().transformSpectralTendencyToSpatial({f.data(),spectral},{expected.data(),volume})),"Stratified QG mean inverse");
    require(std::any_of(expected.begin(),expected.end(),[](double x){return std::abs(x)>1e-3;}) && equal(f,meanBefore),"Stratified QG mean masked or input mutated");
    require(bool(engine->kernel().transformA0ToField({f.data(),spectral},WVStratifiedQGField::qgpv,{expected.data(),volume})),"Stratified QG ordinary mean inverse");
    for (auto x:expected) require(x==0,"Stratified QG ordinary mean mask changed");
    std::vector<double> prepared(2*R),captured(R);
    std::copy(u.begin(),u.end(),prepared.begin()); std::copy(v.begin(),v.end(),prepared.begin()+R);
    const auto preparedBefore=prepared; const auto fluxBefore=f;
    WVRealFieldBundleConstView uv{prepared.data(),{g.Nx,g.Ny,g.Nz,2}};
    WVRealVolumeView capture{captured.data(),volume};
    require(bool(engine->kernel().nonlinearFlux(state,flux,0,&capture,&uv)),"Stratified raw prepared nonlinear capture");
    for (std::size_t i=0;i<R;++i) expected[i]=-(u[i]*qx[i]+v[i]*qy[i]);
    relative(captured,expected,"Stratified prepared nonlinear capture");
    require(prepared==preparedBefore && equal(f,fluxBefore),"Raw capture changed prepared velocity or spectral output");
    uv.shape.fourth=3;
    const auto badShapeCalls=counter->calls;
    require(engine->kernel().nonlinearFlux(state,flux,0,&capture,&uv).code==WVKernelStatusCode::invalidShape && counter->calls==badShapeCalls,
        "Invalid prepared QG velocity accepted or rejected after FFT");
    WVRealVolumeView raw{reinterpret_cast<double*>(a.data()),volume};
    const auto rejectCalls=counter->calls;
    require(engine->kernel().nonlinearFlux(state,flux,0,&raw).code==WVKernelStatusCode::overlappingArrays && counter->calls==rejectCalls,
        "Stratified QG raw alias accepted or rejected after FFT");
    for (const auto evaluationPolicy : {WVVariableEvaluationPolicy::reuse,
                                        WVVariableEvaluationPolicy::lowMemory}) {
        std::unique_ptr<WVStratifiedQGIntegrationSystem> system;
        require(bool(WVStratifiedQGIntegrationSystem::create(
                    source,scheduleValue,
                    wavevortex::runtime::test::extensionCatalog(),
                    std::make_unique<WVReferenceFFTEngine>(),system)) &&
                    bool(system->setVariableEvaluationPolicy(evaluationPolicy)),
                "Stratified QG integration lifecycle setup failed");
        WVCoefficientStateStorage storage,denseStorage;
        require(bool(storage.initialize(system->stateLayout())) &&
                    bool(denseStorage.initialize(system->stateLayout())),
                "Stratified QG integration lifecycle storage failed");
        std::copy_n(a.data(),S,storage.mutableFamilies()[0].data);
        WVMutableIntegrationState integrationState;
        integrationState.waveVortex.t=83;
        integrationState.waveVortex.t0=17;
        integrationState.coefficientFamilies=storage.mutableFamilies();
        integrationState.coefficientFamilyCount=storage.familyCount();
        WVMutableIntegrationState denseState;
        denseState.coefficientFamilies=denseStorage.mutableFamilies();
        denseState.coefficientFamilyCount=denseStorage.familyCount();
        WVFixedStepRK4 integrator(*system,{true});
        require(bool(integrator.prepareStateAfterRestart(integrationState)) &&
                    system->variableEvaluationMetrics().liveBytes==0,
                "Stratified QG restart preparation retained an evaluation");
        const auto validationsBefore=system->kernel().metrics().stateValidationCount;
        const auto contextsBefore=system->variableEvaluationMetrics().contexts;
        require(bool(integrator.step(integrationState,1e-4)),
                "Stratified QG integration lifecycle step failed");
        const auto validationsAfter=system->kernel().metrics().stateValidationCount;
        const auto contextsAfter=system->variableEvaluationMetrics().contexts;
        require(validationsAfter>validationsBefore &&
                    validationsAfter-validationsBefore==contextsAfter-contextsBefore &&
                    system->variableEvaluationMetrics().liveBytes==0,
                "Stratified QG RHS validation and evaluation lifecycles diverged");
        require(bool(integrator.evaluateDenseOutput(
                    integrationState.waveVortex.t-5e-5,denseState)) &&
                    system->kernel().metrics().stateValidationCount-validationsAfter==
                        system->variableEvaluationMetrics().contexts-contextsAfter &&
                    system->variableEvaluationMetrics().liveBytes==0,
                "Stratified QG dense output retained or reopened an RHS evaluation");
    }
    testInjectedQG(source, scheduleValue);
}

template<bool Hydrostatic>
void builtinNonlinearPhysicalOrder(
    const std::shared_ptr<const WVStratifiedModalSource>& source) {
    using Engine=std::conditional_t<Hydrostatic,WVHydrostaticForcingEngine,
        WVBoussinesqForcingEngine>;
    using Field=std::conditional_t<Hydrostatic,WVHydrostaticField,
        WVBoussinesqField>;
    std::unique_ptr<Engine> engine;
    require(bool(Engine::create(source,defaultNonlinearAdvectionSchedule(),
        wavevortex::runtime::test::extensionCatalog(),
        std::make_unique<WVReferenceFFTEngine>(),engine)),
        "Built-in physical-order engine setup");
    std::unique_ptr<WVFieldEvaluationService> service;
    require(bool(WVFieldEvaluationService::createBorrowing(*engine,service)) &&
        bool(service->prepareBuiltinNonlinearCoefficientEvaluation()),
        "Built-in physical-order service setup");
    WVFieldEvaluationPlan maximum,fields,waveField;
    std::vector<WVFieldRequest> maximumRequests(1),fieldRequests(4),waveRequests(1);
    maximumRequests[0].identifier="maximum";
    maximumRequests[0].fieldName="uvMax";
    const std::array<const char*,4> physicalFieldNames{{"u","v","w","eta"}};
    for(std::size_t field=0;field<physicalFieldNames.size();++field) {
      fieldRequests[field].identifier=physicalFieldNames[field];
      fieldRequests[field].fieldName=physicalFieldNames[field];
    }
    waveRequests[0].identifier="wave-u";
    waveRequests[0].fieldName="u_w";
    require(bool(service->createPlan(maximumRequests,maximum)) &&
        bool(service->createPlan(fieldRequests,fields)) &&
        bool(service->createPlan(waveRequests,waveField)),
        "Built-in physical-order field plans");
    WVIntegrationStateLayout layout;
    require(bool(service->createStateLayout({},layout)),
        "Built-in physical-order state layout");
    const auto spectral=engine->kernel().spectralShape();
    const auto spatial=engine->kernel().spatialShape();
    const auto S=spectral.elementCount(),R=spatial.elementCount();
    std::vector<WVComplex64> coefficients(3*S);
    for(std::size_t index=0;index<coefficients.size();++index)
      coefficients[index]={1e-5*std::sin(.17*(index+1)),
          1e-5*std::cos(.23*(index+1))};
    WVMutableCoefficients mutableCoefficients{{coefficients.data(),spectral},
        {coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}};
    require(bool(engine->kernel().constrainCoefficients(mutableCoefficients)),
        "Built-in physical-order coefficient constraints");
    const WVState state{3,-1,{{coefficients.data(),spectral},
        {coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}}};
    std::array<WVCoefficientFamilyConstView,3> families;
    for(std::size_t family=0;family<3;++family)
      families[family]={&layout.coefficientFamilies()[family],
          coefficients.data()+family*S};
    const WVIntegrationState integration{state,nullptr,0,families.data(),3};
    std::array<std::vector<double>,4> reference;
    constexpr Field fieldNames[]={Field::u,Field::v,Field::w,Field::eta};
    for(std::size_t field=0;field<4;++field) {
      reference[field].resize(R);
      require(bool(engine->kernel().transformStateField(state,fieldNames[field],
          {reference[field].data(),spatial})),
          "Built-in physical-order field reference");
    }
    std::vector<double> waveReference(R),waveValue(R);
    if constexpr(Hydrostatic) {
      require(bool(engine->kernel().transformStateField(state,Field::u,
          {waveReference.data(),spatial},WVHydrostaticDerivative::value,
          WVHydrostaticComponent::wave)),
          "Built-in physical-order wave reference");
    } else {
      require(bool(engine->kernel().transformStateField(state,Field::u,
          {waveReference.data(),spatial},WVBoussinesqDerivative::value,
          WVBoussinesqComponent::wave)),
          "Built-in physical-order wave reference");
    }
    auto reconstructionCount=[&] {
      const auto& counts=engine->kernel().metrics().reconstructionCount;
      return counts[0][0][0]+counts[1][0][0]+counts[2][0][0]+counts[3][0][0];
    };
    std::vector<WVComplex64> coefficientValues(3*S);
    WVFlux flux{{coefficientValues.data(),spectral},
        {coefficientValues.data()+S,spectral},
        {coefficientValues.data()+2*S,spectral}};
    std::array<std::vector<double>,4> actual;
    std::array<WVFieldOutputView,4> fieldViews;
    for(std::size_t field=0;field<4;++field) {
      actual[field].resize(R);
      fieldViews[field]={actual[field].data(),R};
    }
    double maximumValue=0;
    WVFieldOutputView maximumView{&maximumValue,1};
    {
      WVFieldEvaluationSession session;
      require(bool(service->beginEvaluationSession(integration,session)),
          "Built-in coefficient-first session");
      const auto before=reconstructionCount();
      require(bool(service->evaluateBuiltinNonlinearCoefficients(integration,flux)),
          "Built-in coefficient-first evaluation");
      const auto afterCoefficient=reconstructionCount();
      if(afterCoefficient!=before+4)
        throw std::runtime_error("Built-in coefficient-first reconstructions before="+
            std::to_string(before)+" after="+std::to_string(afterCoefficient));
      require(bool(service->evaluate(fields,integration,fieldViews.data(),4)) &&
          bool(service->evaluate(maximum,integration,&maximumView,1)) &&
          reconstructionCount()==afterCoefficient,
          "Canonical physical slices were reconstructed by later field consumers");
      for(std::size_t field=0;field<4;++field)
        relative(actual[field],reference[field],
            "Canonical physical bundle channel mapping");
      WVFieldOutputView waveView{waveValue.data(),R};
      const auto waveStatus=service->evaluate(waveField,integration,&waveView,1);
      require(bool(waveStatus),"Canonical component field evaluation");
      relative(waveValue,waveReference,"Canonical physical component identity");
      double componentDifference=0;
      for(std::size_t index=0;index<R;++index)
        componentDifference=std::max(componentDifference,
            std::abs(waveValue[index]-reference[0][index]));
      require(componentDifference>0,
          "Canonical component fixture did not distinguish total and wave fields");
    }
    {
      WVFieldEvaluationSession session;
      require(bool(service->beginEvaluationSession(integration,session)),
          "Built-in maximum-first session");
      const auto before=reconstructionCount();
      require(bool(service->evaluate(maximum,integration,&maximumView,1)) &&
          reconstructionCount()==before+2,
          "Maximum-first path did not retain its two velocity fields");
      require(bool(service->evaluateBuiltinNonlinearCoefficients(integration,flux)) &&
          reconstructionCount()==before+4,
          "Built-in coefficients repeated fields already retained by maximum speed");
      require(bool(service->evaluate(fields,integration,fieldViews.data(),4)) &&
          reconstructionCount()==before+4,
          "Gathered canonical fields were not published for later consumers");
      for(std::size_t field=0;field<4;++field)
        relative(actual[field],reference[field],
            "Gathered canonical physical bundle channel mapping");
    }
}

template<bool Hydrostatic> void stratified(bool unpairedMean=false) {
    Temporary file;
    if constexpr (Hydrostatic) {
        fixture(file.path); File f(file.path);
        for (const char* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformHydrostatic"),"WVTransformHydrostatic"));
    } else boussinesqFixture(file.path);
    change(file.path,"shouldAntialias",0);
    std::shared_ptr<const WVStratifiedModalRecord> source;
    require(bool(WVStratifiedModalReader::read(file.path.string(),source)),"Stratified fixture read");
    using Engine=std::conditional_t<Hydrostatic,WVHydrostaticForcingEngine,WVBoussinesqForcingEngine>;
    const auto& g=source->geometry(); const WVShape2D spectral{g.Nj,g.Nkl};
    const auto S=spectral.elementCount(),R=g.Nx*g.Ny*g.Nz;
    const WVShape4D spatial{g.Nx,g.Ny,g.Nz,Hydrostatic?3U:4U}; const WVShape3D volume{g.Nx,g.Ny,g.Nz};
    auto counter=std::make_shared<FailureCounter>(); std::unique_ptr<Engine> engine;
    builtinNonlinearPhysicalOrder<Hydrostatic>(source);
    require(bool(Engine::create(source,schedule(S,unpairedMean),wavevortex::runtime::test::extensionCatalog(),std::make_unique<FailingEngine>(counter),engine)),"Stratified forcing setup");
    exercise(*engine,spectral,spatial,counter,
        [&](const WVState& state,const std::vector<double>& sum,std::vector<WVComplex64>& out) {
            WVMutableCoefficients output{{out.data(),spectral},{out.data()+S,spectral},{out.data()+2*S,spectral}};
            if constexpr (Hydrostatic) return engine->kernel().transformUVEtaToWaveVortex({sum.data(),volume},{sum.data()+R,volume},{sum.data()+2*R,volume},state.t,state.t0,output);
            else return engine->kernel().transformUVWEtaToWaveVortex({sum.data(),volume},{sum.data()+R,volume},{sum.data()+2*R,volume},{sum.data()+3*R,volume},state.t,state.t0,output);
        },
        [&](const WVState& state,const std::vector<WVComplex64>& delta,std::vector<double>& out) {
            const WVState d{state.t,state.t0,{{delta.data(),spectral},{delta.data()+S,spectral},{delta.data()+2*S,spectral}}};
            for (std::size_t channel=0;channel<spatial.fourth;++channel) {
                WVKernelStatus status;
                if constexpr (Hydrostatic) {
                    const WVHydrostaticField fields[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::eta};
                    status=engine->kernel().transformStateField(d,fields[channel],{out.data()+channel*R,volume});
                } else {
                    const WVBoussinesqField fields[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
                    status=engine->kernel().transformStateField(d,fields[channel],{out.data()+channel*R,volume});
                }
                if (!status) return status;
            }
            return WVKernelStatus::ok();
        },unpairedMean);
}

WVFrozenForcingSchedule repeatedGridCalculusSchedule() {
    WVPortableTypedRecord damping{"wave-vortex-forcing-configuration-v1",1,
        {{"nu",{},std::vector<double>{0.125}},
         {"kappa",{},std::vector<double>{0.25}}}};
    WVPortableTypedRecord diffusivity{"wave-vortex-forcing-configuration-v1",1,
        {{"kappa_z",{},std::vector<double>{0.375}},
         {"shouldForceMeanDensityAnomaly",{},std::vector<std::uint8_t>{1}}}};
    WVFrozenForcingSchedule result;
    result.entries={
        {"WVHorizontalDamping",1,"horizontal one",WVForcingStage::spatial,120,0,"",damping},
        {"WVHorizontalDamping",1,"horizontal two",WVForcingStage::spatial,121,1,"",damping},
        {"WVVerticalDamping",1,"vertical one",WVForcingStage::spatial,122,2,"",damping},
        {"WVVerticalDamping",1,"vertical two",WVForcingStage::spatial,123,3,"",damping},
        {"WVVerticalDiffusivity",1,"vertical diffusivity",WVForcingStage::spatial,124,4,"",diffusivity}};
    return result;
}

template<bool Hydrostatic>
void repeatedGridCalculus() {
    Temporary file;
    if constexpr(Hydrostatic) {
        fixture(file.path); File handle(file.path);
        for(const char* name:{"WVTransform","AnnotatedClass"})
            nc(nc_put_att_text(handle.id,NC_GLOBAL,name,
                std::char_traits<char>::length("WVTransformHydrostatic"),
                "WVTransformHydrostatic"));
    } else boussinesqFixture(file.path);
    change(file.path,"shouldAntialias",0);
    std::shared_ptr<const WVStratifiedModalRecord> source;
    require(bool(WVStratifiedModalReader::read(file.path.string(),source)),
        "Repeated grid-calculus fixture read failed");
    using Engine=std::conditional_t<Hydrostatic,WVHydrostaticForcingEngine,
        WVBoussinesqForcingEngine>;
    std::unique_ptr<Engine> engine;
    const auto scheduleValue=repeatedGridCalculusSchedule();
    require(bool(Engine::create(source,scheduleValue,
        wavevortex::runtime::test::extensionCatalog(),
        std::make_unique<WVReferenceFFTEngine>(),engine)),
        "Repeated grid-calculus engine setup failed");
    const auto spectral=engine->kernel().spectralShape();
    const auto spatial=engine->kernel().spatialShape();
    const auto S=spectral.elementCount(),R=spatial.elementCount();
    std::vector<WVComplex64> coefficients(3*S);
    for(std::size_t index=0;index<coefficients.size();++index)
        coefficients[index]={1e-5*std::sin(.17*(index+1)),
            1e-5*std::cos(.23*(index+1))};
    require(bool(engine->kernel().constrainCoefficients(
        {{coefficients.data(),spectral},{coefficients.data()+S,spectral},
         {coefficients.data()+2*S,spectral}})),
        "Repeated grid-calculus coefficient constraint failed");
    const WVState state{.5,-.25,{{coefficients.data(),spectral},
        {coefficients.data()+S,spectral},{coefficients.data()+2*S,spectral}}};
    const auto channels=Hydrostatic ? 3U : 4U;
    const WVShape4D outputShape{spatial.first,spatial.second,spatial.third,channels};

    auto expected=[&](WVVariableEvaluationPolicy policy) {
        std::array<std::array<std::size_t,4>,4> value{};
        const auto repeated=policy==WVVariableEvaluationPolicy::reuse ? 1U : 2U;
        for(std::size_t field=0;field<4;++field)
            if(!Hydrostatic || field!=2)
                for(std::size_t kind=0;kind<2;++kind) value[field][kind]=repeated;
        value[0][2]=repeated; value[1][2]=repeated;
        if constexpr(!Hydrostatic) value[2][3]=repeated;
        value[3][3]=policy==WVVariableEvaluationPolicy::reuse ? 1U : 3U;
        return value;
    };
    auto checkDelta=[&](const auto& before,WVVariableEvaluationPolicy policy,
        const char* label) {
        const auto desired=expected(policy);
        const auto& after=engine->metrics().gridCalculusProducerCount;
        for(std::size_t field=0;field<4;++field)
            for(std::size_t kind=0;kind<4;++kind)
                if(after[field][kind]-before[field][kind]!=desired[field][kind]) {
                    std::ostringstream message;
                    message<<label<<" grid-calculus producer count mismatch field="
                        <<field<<" kind="<<kind<<" actual="
                        <<after[field][kind]-before[field][kind]<<" expected="
                        <<desired[field][kind];
                    throw std::runtime_error(message.str());
                }
    };
    auto evaluateRHS=[&](WVVariableEvaluationPolicy policy,
        std::vector<WVComplex64>& values) {
        values.assign(3*S,{});
        WVFlux flux{{values.data(),spectral},{values.data()+S,spectral},
            {values.data()+2*S,spectral}};
        const auto before=engine->metrics().gridCalculusProducerCount;
        const auto status=engine->nonlinearFlux(state,flux);
        if(!status) throw std::runtime_error(
            "Repeated grid-calculus RHS failed: "+status.message);
        checkDelta(before,policy,"RHS");
    };
    auto evaluateDiagnostics=[&](WVVariableEvaluationPolicy policy,
        std::vector<double>& values) {
        values.assign(scheduleValue.entries.size()*channels*R,0);
        std::vector<WVForcingTendencyOutput> outputs;
        outputs.reserve(scheduleValue.entries.size());
        for(std::size_t index=0;index<scheduleValue.entries.size();++index)
            outputs.push_back({index,{values.data()+index*channels*R,outputShape}});
        const auto before=engine->metrics().gridCalculusProducerCount;
        require(bool(engine->evaluateForcingTendencies(
            state,outputs.data(),outputs.size())),
            "Repeated grid-calculus diagnostics failed");
        checkDelta(before,policy,"diagnostic");
    };

    std::vector<WVComplex64> reusedRHS,lowRHS;
    std::vector<double> reusedDiagnostics,lowDiagnostics;
    evaluateRHS(WVVariableEvaluationPolicy::reuse,reusedRHS);
    evaluateDiagnostics(WVVariableEvaluationPolicy::reuse,reusedDiagnostics);
    require(bool(engine->setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::lowMemory)),
        "Repeated grid-calculus low-memory setup failed");
    evaluateRHS(WVVariableEvaluationPolicy::lowMemory,lowRHS);
    evaluateDiagnostics(WVVariableEvaluationPolicy::lowMemory,lowDiagnostics);
    require(equal(reusedRHS,lowRHS) && reusedDiagnostics==lowDiagnostics,
        "Grid-calculus reuse policy changed forcing values");
}
}
int main() {
    try {
        builtinNonlinearCoefficientBoundary();
        repeatedGridCalculus<true>(); repeatedGridCalculus<false>();
        stratifiedQG(); barotropic(); constant(false); constant(true); stratified<true>(); stratified<false>(); stratified<true>(true); stratified<false>(true);
        std::cout<<"PASS: ordered forcing diagnostics, state isolation, and transactional retry\n";
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
