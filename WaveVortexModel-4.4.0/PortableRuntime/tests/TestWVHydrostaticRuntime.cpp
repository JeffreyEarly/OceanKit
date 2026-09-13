#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVHydrostaticIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVStratifiedModalTestFixture.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <stdexcept>
#include <array>
#include <cmath>
using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::test_fixture;
namespace {
int injectedFactoryCalls = 0;
WVKernelStatus countingScalarBackend(std::unique_ptr<WVVerticalMatrixBackend> &backend) {
    ++injectedFactoryCalls;
    return WVCreateScalarMatrixBackend(backend);
}

#if WV_TEST_NATIVE_FFTW
void nativeTiledLifecycle(std::shared_ptr<const WVStratifiedModalRecord> source,
                          std::shared_ptr<const WVExtensionCatalog> catalog,
                          const WVFrozenForcingSchedule &schedule) {
    const auto createEngine = [&](bool tiled,bool nonlinear=true) {
        WVVariableKernelServices services;
        services.execution.horizontalSchedule =
            WVRetainedHorizontalSchedule::streamingPrunedTile16;
        services.execution.horizontalWorkers = 2;
        services.execution.streamedNonlinear = true;
        services.execution.spectralSchedule =
            WVVariableSpectralSchedule::compactSplitFusedViews;
        services.execution.pointwiseWorkers = 2;
        services.execution.sharedFieldGradients = true;
        services.execution.fusedDerivativeAdvection = true;
        services.execution.tiledNonlinear = tiled;
        std::unique_ptr<WVFFTEngine> fft;
        require(bool(WVFFTWEngine::create(1,fft)),
                "Native Hydrostatic tiled FFT fixture failed");
        std::unique_ptr<WVHydrostaticForcingEngine> result;
        auto selectedSchedule=schedule;
        if (!nonlinear) selectedSchedule.entries.clear();
        const auto status = WVHydrostaticForcingEngine::create(
            source,selectedSchedule,catalog,std::move(fft),result,services);
        require(bool(status),status.message.c_str());
        return result;
    };
    auto noNonlinear=createEngine(true,false);
    require(!noNonlinear->kernel().supportsTiledNonlinear(),"Unused tiled workspace was prepared");
    auto baseline = createEngine(false);
    auto candidate = createEngine(true);
    require(!baseline->kernel().supportsTiledNonlinear() &&
                candidate->kernel().supportsTiledNonlinear(),
            "Hydrostatic tiled option did not control native capability");

    const auto spectral = candidate->kernel().spectralShape();
    const auto spatial = candidate->kernel().spatialShape();
    const auto S = spectral.elementCount(), R = spatial.elementCount();
    std::array<std::vector<WVComplex64>,3> coefficients;
    for (std::size_t family=0;family<coefficients.size();++family) {
        coefficients[family].resize(S);
        for (std::size_t index=0;index<S;++index)
            coefficients[family][index]={
                1e-3*std::sin(.17*static_cast<double>(index+3*family+1)),
                1e-3*std::cos(.11*static_cast<double>(2*index+family+1))};
    }
    WVMutableCoefficients mutableState{{coefficients[0].data(),spectral},
        {coefficients[1].data(),spectral},{coefficients[2].data(),spectral}};
    require(bool(baseline->kernel().constrainCoefficients(mutableState)),
            "Hydrostatic tiled coefficient constraint failed");
    const auto stateView = [&] {
        return WVState{.37,.11,{{coefficients[0].data(),spectral},
            {coefficients[1].data(),spectral},
            {coefficients[2].data(),spectral}}};
    };
    const auto fluxView = [&](std::array<std::vector<WVComplex64>,3>& values) {
        for (auto& family:values) family.resize(S);
        return WVFlux{{values[0].data(),spectral},{values[1].data(),spectral},
            {values[2].data(),spectral}};
    };
    const auto requireClose = [](const auto& actual,const auto& expected,
                                 const char* message) {
        require(actual.size()==expected.size(),message);
        for (std::size_t index=0;index<actual.size();++index)
            require(std::abs(actual[index]-expected[index])<1e-12,message);
    };
    const auto requireComplexClose = [](const auto& actual,const auto& expected,
                                        const char* message) {
        for (std::size_t family=0;family<actual.size();++family) {
            require(actual[family].size()==expected[family].size(),message);
            for (std::size_t index=0;index<actual[family].size();++index) {
                require(std::abs(actual[family][index].real-
                    expected[family][index].real)<1e-12,message);
                require(std::abs(actual[family][index].imag-
                    expected[family][index].imag)<1e-12,message);
            }
        }
    };

    std::array<std::vector<WVComplex64>,3> baselineFlux,candidateFlux;
    auto baselineOutput=fluxView(baselineFlux);
    auto candidateOutput=fluxView(candidateFlux);
    auto state=stateView();
    require(bool(baseline->beginStateEvaluation(state)) &&
                bool(baseline->nonlinearFlux(state,baselineOutput)),
            "Hydrostatic option-off scoped RHS failed");
    WVRealFieldBundleConstView baselineFields;
    require(bool(baseline->physicalFields(state,baselineFields)),
            "Hydrostatic option-off physical fields failed");
    std::vector<double> expectedFields(baselineFields.data,
        baselineFields.data+4*R);
    baseline->endStateEvaluation();
    require(baseline->kernel().metrics().tiledNonlinearCount==0 &&
                baseline->variableEvaluationMetrics().producerExecutions==4,
            "Hydrostatic option-off path executed tiled nonlinear work");

    candidate->kernel().resetMetrics();
    require(bool(candidate->beginStateEvaluation(state)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)),
            "Hydrostatic grouped tiled RHS failed");
    const auto groupedProducers=candidate->variableEvaluationMetrics().producerExecutions;
    const auto groupedMetrics=candidate->kernel().metrics();
    require(groupedProducers==1 && groupedMetrics.tiledNonlinearCount==1 &&
                groupedMetrics.tiledColumnInverseCount==10*spatial.third &&
                groupedMetrics.tiledRowInverseCount==13*spatial.third &&
                groupedMetrics.tiledReusedColumnCount==3*spatial.third,
            "Hydrostatic grouped tiled producer counts changed");
    WVRealFieldBundleConstView groupedFields;
    require(bool(candidate->physicalFields(state,groupedFields)) &&
                candidate->variableEvaluationMetrics().producerExecutions==
                    groupedProducers,
            "Hydrostatic grouped physical fields were not cached");
    std::vector<double> actualFields(groupedFields.data,groupedFields.data+4*R);
    requireClose(actualFields,expectedFields,
                 "Hydrostatic tiled physical fields differ");
    requireComplexClose(candidateFlux,baselineFlux,
                        "Hydrostatic tiled RHS differs");
    double uv=0,w=0;
    require(bool(candidate->speedMaxima(state,uv,w)),
            "Hydrostatic grouped speed maxima failed");
    const auto reducedProducers=candidate->variableEvaluationMetrics().producerExecutions;
    require(bool(candidate->speedMaxima(state,uv,w)) &&
                candidate->variableEvaluationMetrics().producerExecutions==
                    reducedProducers &&
                candidate->kernel().metrics().tiledNonlinearCount==1,
            "Hydrostatic grouped fields were not reused by reductions");
    candidate->endStateEvaluation();
    require(candidate->variableEvaluationMetrics().liveBytes==0,
            "Hydrostatic grouped fields survived their scope");

    candidate->kernel().resetMetrics();
    const auto readyProducers=
        candidate->variableEvaluationMetrics().producerExecutions;
    require(bool(candidate->beginStateEvaluation(state)),
            "Hydrostatic ready-first scope failed");
    WVRealFieldBundleConstView readyFields;
    require(bool(candidate->physicalFields(state,readyFields)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)) &&
                candidate->kernel().metrics().tiledNonlinearCount==0 &&
                candidate->variableEvaluationMetrics().producerExecutions-
                    readyProducers==4,
            "Hydrostatic ready-first request did not use the fallback");
    requireComplexClose(candidateFlux,baselineFlux,
                        "Hydrostatic ready-first fallback RHS differs");
    candidate->endStateEvaluation();

    candidate->kernel().resetMetrics();
    const auto partialProducers=
        candidate->variableEvaluationMetrics().producerExecutions;
    require(bool(candidate->beginStateEvaluation(state)) &&
                bool(candidate->speedMaxima(state,uv,w)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)) &&
                candidate->kernel().metrics().tiledNonlinearCount==0 &&
                candidate->variableEvaluationMetrics().producerExecutions-
                    partialProducers==6,
            "Hydrostatic partial-ready request did not use the fallback");
    WVRealFieldBundleConstView completedFields;
    require(bool(candidate->physicalFields(state,completedFields)),
            "Hydrostatic partial-ready fallback did not complete fields");
    requireComplexClose(candidateFlux,baselineFlux,
                        "Hydrostatic partial-ready fallback RHS differs");
    requireClose(std::vector<double>(completedFields.data,
                     completedFields.data+4*R),expectedFields,
                 "Hydrostatic partial-ready fallback fields differ");
    candidate->endStateEvaluation();

    require(bool(candidate->setVariableEvaluationPolicy(
                WVVariableEvaluationPolicy::lowMemory)),
            "Hydrostatic tiled low-memory policy failed");
    candidate->kernel().resetMetrics();
    const auto lowMemoryProducers=
        candidate->variableEvaluationMetrics().producerExecutions;
    require(bool(candidate->beginStateEvaluation(state)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)) &&
                candidate->kernel().metrics().tiledNonlinearCount==0 &&
                candidate->variableEvaluationMetrics().producerExecutions-
                    lowMemoryProducers==5,
            "Hydrostatic low-memory request executed tiled nonlinear work");
    requireComplexClose(candidateFlux,baselineFlux,
                        "Hydrostatic low-memory fallback RHS differs");
    candidate->endStateEvaluation();
    require(bool(candidate->setVariableEvaluationPolicy(
                WVVariableEvaluationPolicy::reuse)),
            "Hydrostatic tiled reuse policy restoration failed");

    for (auto& family:coefficients)
        for (auto& value:family) {value.real*=.5;value.imag*=.5;}
    state=stateView();
    candidate->kernel().resetMetrics();
    require(bool(candidate->beginStateEvaluation(state)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)),
            "Hydrostatic changed in-place state tiled RHS failed");
    WVRealFieldBundleConstView changedFields;
    require(bool(candidate->physicalFields(state,changedFields)) &&
                candidate->kernel().metrics().tiledNonlinearCount==1,
            "Hydrostatic changed in-place state reused the prior scope");
    bool changed=false;
    for (std::size_t index=0;index<4*R;++index)
        changed|=changedFields.data[index]!=actualFields[index];
    require(changed,"Hydrostatic changed in-place state retained old fields");
    candidate->endStateEvaluation();

    std::size_t meanDensity=S;
    for (std::size_t index=0;index<S;++index)
        if (candidate->kernel().factors()[index].meanDensityAnomaly) {
            meanDensity=index;
            break;
        }
    require(meanDensity<S,
            "Hydrostatic tiled fixture has no mean-density zero mode");
    const auto validMeanImag=coefficients[2][meanDensity].imag;
    coefficients[2][meanDensity].imag=1e-4;
    state=stateView();
    for (auto& family:candidateFlux)
        std::fill(family.begin(),family.end(),WVComplex64{17,19});
    const auto failedProducers=
        candidate->variableEvaluationMetrics().producerExecutions;
    candidate->kernel().resetMetrics();
    require(bool(candidate->beginStateEvaluation(state)),
            "Hydrostatic malformed tiled scope failed to begin");
    const auto malformed=candidate->nonlinearFlux(state,candidateOutput);
    require(malformed.code==WVKernelStatusCode::invalidConfiguration &&
                candidate->variableEvaluationMetrics().producerExecutions==
                    failedProducers &&
                candidate->variableEvaluationMetrics().liveBytes==0 &&
                candidate->kernel().metrics().tiledNonlinearCount==0,
            "Hydrostatic malformed tiled fields were published");
    for (const auto& family:candidateFlux)
        for (const auto& value:family)
            require(value.real==0 && value.imag==0,
                    "Hydrostatic malformed tiled flux was partially accumulated");
    candidate->endStateEvaluation();

    coefficients[2][meanDensity].imag=validMeanImag;
    state=stateView();
    require(bool(baseline->beginStateEvaluation(state)) &&
                bool(baseline->nonlinearFlux(state,baselineOutput)),
            "Hydrostatic restored option-off oracle failed");
    baseline->endStateEvaluation();
    candidate->kernel().resetMetrics();
    require(bool(candidate->beginStateEvaluation(state)) &&
                bool(candidate->nonlinearFlux(state,candidateOutput)) &&
                candidate->kernel().metrics().tiledNonlinearCount==1,
            "Hydrostatic restored tiled state did not recover");
    requireComplexClose(candidateFlux,baselineFlux,
                        "Hydrostatic restored tiled RHS differs");
    candidate->endStateEvaluation();
}
#endif

void injectedServices(std::shared_ptr<const WVStratifiedModalRecord> source,
                      std::shared_ptr<const WVExtensionCatalog> catalog,
                      const WVFrozenForcingSchedule &schedule) {
    WVVariableKernelServices services;
    services.matrixBackendFactory = countingScalarBackend;
    services.execution = {WVRetainedHorizontalSchedule::streamingPrunedTile16, 2,
                          true, WVVariableSpectralSchedule::compactSplitFusedViews};

    std::unique_ptr<WVHydrostaticForcingEngine> baseline, injected;
    require(bool(WVHydrostaticForcingEngine::create(
                    source, schedule, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), baseline)),
            "Baseline Hydrostatic injection fixture failed");
    require(bool(WVHydrostaticForcingEngine::create(
                    source, schedule, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), injected,
                    services)),
            "Injected Hydrostatic fixture failed");
    require(injectedFactoryCalls == 4, "Injected backend factory call count differs");
    require(injected->kernel().executionOptions().usesCompactSplitViews(),
            "Injected Hydrostatic execution options were not retained");
    const auto shape = baseline->kernel().spectralShape();
    const auto count = shape.elementCount();
    std::array<std::vector<WVComplex64>, 3> coefficients, referenceFlux, injectedFlux;
    for (auto &values : coefficients) values.resize(count, {.001, -.002});
    WVMutableCoefficients mutableCoefficients{{coefficients[0].data(), shape},
                                               {coefficients[1].data(), shape},
                                               {coefficients[2].data(), shape}};
    require(bool(baseline->kernel().constrainCoefficients(mutableCoefficients)),
            "Hydrostatic injection coefficient constraint failed");
    for (auto &values : referenceFlux) values.resize(count);
    for (auto &values : injectedFlux) values.resize(count);
    WVState state{.37, .11, {{coefficients[0].data(), shape},
                             {coefficients[1].data(), shape},
                             {coefficients[2].data(), shape}}};
    WVFlux reference{{referenceFlux[0].data(), shape}, {referenceFlux[1].data(), shape},
                     {referenceFlux[2].data(), shape}};
    WVFlux actual{{injectedFlux[0].data(), shape}, {injectedFlux[1].data(), shape},
                  {injectedFlux[2].data(), shape}};
    require(bool(baseline->nonlinearFlux(state, reference)), "Baseline Hydrostatic RHS failed");
    require(bool(injected->nonlinearFlux(state, actual)), "Injected Hydrostatic RHS failed");
    for (std::size_t family = 0; family < 3; ++family)
        for (std::size_t i = 0; i < count; ++i) {
            require(std::abs(referenceFlux[family][i].real - injectedFlux[family][i].real) < 1e-12,
                    "Injected Hydrostatic real RHS differs");
            require(std::abs(referenceFlux[family][i].imag - injectedFlux[family][i].imag) < 1e-12,
                    "Injected Hydrostatic imaginary RHS differs");
        }

    std::unique_ptr<WVHydrostaticIntegrationSystem> system;
    injectedFactoryCalls = 0;
    require(bool(WVHydrostaticIntegrationSystem::create(
                    source, schedule, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), system, services)),
            "Injected Hydrostatic integration fixture failed");
    require(injectedFactoryCalls == 4 &&
                system->kernel().executionOptions().usesCompactSplitViews(),
            "Injected Hydrostatic integration services were not forwarded");
    WVPortableObserverRecord record;
    for (const auto &family : system->stateLayout().coefficientFamilies())
        record.stateBlocks.push_back({family.identifier, WVStateScalarType::complex64,
                                      family.spectralDimensions,
                                      WVToleranceKind::coefficientEnergyScaled, 1e-6,
                                      WVStateOwnership::integratorOwned,
                                      WVRestartRequirement::requiredDynamicState});
    WVPortableObserverDescriptor descriptor;
    require(bool(WVPortableObserverDescriptor::create(record, catalog, descriptor)),
            "Hydrostatic injection descriptor creation failed");
    injectedFactoryCalls = 0;
    require(bool(WVHydrostaticIntegrationSystem::create(
                    source, schedule, descriptor, catalog,
                    std::make_unique<WVReferenceFFTEngine>(), system, services)),
            "Injected Hydrostatic descriptor fixture failed");
    require(injectedFactoryCalls == 4 &&
                system->kernel().executionOptions().usesCompactSplitViews(),
            "Injected Hydrostatic descriptor services were not forwarded");
    WVVariableKernelServices rejected;
    rejected.matrixBackendFactory = {};
    auto *old = injected.get();
    require(!WVHydrostaticForcingEngine::create(
                source, schedule, catalog,
                std::make_unique<WVReferenceFFTEngine>(), injected, rejected) &&
                injected.get() == old,
            "Empty Hydrostatic backend factory replaced the existing engine");
}

void contracts(std::shared_ptr<const WVStratifiedModalRecord> source) {
    WVExtensionCatalogBuilder builder; require(bool(addBuiltInExtensions(builder)),"Built-ins failed");
    std::shared_ptr<const WVExtensionCatalog> catalog; require(bool(builder.freeze(catalog)),"Catalog failed");
    WVFrozenForcingSchedule schedule;
    const auto* registration=catalog->forcings().registration("WVNonlinearAdvection",1);
    WVFrozenForcingEntry entry; entry.typeIdentifier=registration->matlabClassName; entry.contractVersion=1; entry.name=registration->defaultName; entry.stage=registration->stage; entry.priority=registration->priority;
    entry.configuration={"wave-vortex-forcing-configuration-v1",1,{}}; schedule.entries.push_back(entry);
    std::unique_ptr<WVHydrostaticForcingEngine> engine;
    auto status=WVHydrostaticForcingEngine::create(source,schedule,catalog,std::make_unique<WVReferenceFFTEngine>(),engine); require(bool(status),status.message.c_str());
    const auto shape=engine->kernel().spectralShape(); const auto volume=engine->kernel().spatialShape();
    const auto S=shape.elementCount(),R=volume.elementCount();
    std::array<std::vector<WVComplex64>,3> a,b;
    for(auto& x:a) x.resize(S,{.001,.002});
    for(auto& x:b) x.resize(S,{17,19});
    WVMutableCoefficients amplitudes{{a[0].data(),shape},{a[1].data(),shape},{a[2].data(),shape}};
    require(bool(engine->kernel().constrainCoefficients(amplitudes)),"Constraint failed");
    WVState state{83,17,{{a[0].data(),shape},{a[1].data(),shape},{a[2].data(),shape}}}; WVFlux flux{{b[0].data(),shape},{b[1].data(),shape},{b[2].data(),shape}};
    auto invalidState=state; invalidState.t=std::numeric_limits<double>::quiet_NaN();
    require(!engine->nonlinearFlux(invalidState,flux),"Invalid time accepted");
    auto invalidFlux=flux; invalidFlux.F0.shape={S,1}; require(!engine->nonlinearFlux(state,invalidFlux),"Invalid shape accepted");
    invalidFlux=flux; invalidFlux.Fm=flux.Fp; require(!engine->nonlinearFlux(state,invalidFlux),"Aliased outputs accepted");
    for(const auto& x:b) for(auto value:x) require(value.real==17 && value.imag==19,"Invalid RHS mutated output");
    // Exact same-family in-place evolution is a core feature, but an RHS must
    // never overwrite the caller's state before applying forcing.
    invalidFlux=flux; invalidFlux.Fp=amplitudes.Ap;
    require(!engine->nonlinearFlux(state,invalidFlux),"RHS aliased caller state");
    std::vector<double> scalar(R,.3),out(R,29),velocities(3*R);
    for(std::size_t i=0;i<R;++i) scalar[i]+=.1*std::sin(.31*i);
    for(std::size_t i=0;i<3*R;++i) velocities[i]=.01*std::cos(.13*i);
    WVRealFieldBundleConstView fields{velocities.data(),{volume.first,volume.second,volume.third,3}};
    auto badFields=fields; badFields.shape.fourth=2;
    require(!engine->kernel().advectScalarWithAdvectionFields({scalar.data(),volume},badFields,false,{out.data(),volume}),"Invalid tracer field shape accepted");
    require(!engine->kernel().advectScalarWithAdvectionFields({scalar.data(),volume},fields,false,{scalar.data(),volume}),"Aliased tracer output accepted");
    for(auto value:out) require(value==29,"Invalid tracer operation mutated output");
    // One immutable scope spans forcing and subsequent passive consumers.
    for(auto evaluationPolicy:{WVVariableEvaluationPolicy::reuse,WVVariableEvaluationPolicy::lowMemory}) {
        require(bool(engine->setVariableEvaluationPolicy(evaluationPolicy)),"Policy selection failed");
        engine->kernel().resetMetrics();
        require(bool(engine->beginStateEvaluation(state)),"Evaluation scope failed");
        require(!engine->beginStateEvaluation(state),"Nested evaluation accepted");
        require(!engine->setVariableEvaluationPolicy(evaluationPolicy),"Active policy change accepted");
        require(bool(engine->nonlinearFlux(state,flux)),"Scoped RHS failed");
        const auto produced=engine->kernel().metrics().fieldReconstructionCount;
        WVRealFieldBundleConstView reused;
        require(bool(engine->physicalFields(state,reused)),"Passive consumer failed");
        double uv=0,w=0;
        require(bool(engine->speedMaxima(state,uv,w)),"Speed reduction failed");
        const auto afterReduction=engine->variableEvaluationMetrics().producerExecutions;
        require(bool(engine->speedMaxima(state,uv,w)),"Speed reduction reuse failed");
        require(engine->variableEvaluationMetrics().producerExecutions==afterReduction,"Repeated speed reduction executed");
        require(engine->kernel().metrics().fieldReconstructionCount==produced,"Passive consumer reconstructed fields");
        require(engine->kernel().metrics().stateValidationCount==1,"State validated more than once");
        require(engine->kernel().metrics().phasePreparationCount==1,"Phase prepared more than once");
        auto foreign=state; foreign.t+=1;
        require(!engine->physicalFields(foreign,reused),"Foreign state used active cache");
        engine->endStateEvaluation();
        require(engine->variableEvaluationMetrics().liveBytes==0,"Cache validity survived scope");
        a[0][0].real+=.001;
        require(bool(engine->kernel().constrainCoefficients(amplitudes)),"Mutated state constraint failed");
        const auto newStateStatus=engine->physicalFields(state,reused);
        require(bool(newStateStatus),newStateStatus.message.c_str());
        require(engine->kernel().metrics().stateValidationCount==2,"New state reused validation");
    }
    require(bool(engine->setVariableEvaluationPolicy(WVVariableEvaluationPolicy::reuse)),"Default restoration failed");
    auto transactionSchedule=schedule;
    auto secondNonlinear=entry;
    secondNonlinear.name="second nonlinear advection";
    secondNonlinear.ordinal=1;
    transactionSchedule.entries.push_back(std::move(secondNonlinear));
    std::unique_ptr<WVHydrostaticForcingEngine> transactionEngine;
    require(bool(WVHydrostaticForcingEngine::create(source,transactionSchedule,
        catalog,std::make_unique<WVReferenceFFTEngine>(),transactionEngine)) &&
        bool(transactionEngine->setVariableEvaluationPolicy(
            WVVariableEvaluationPolicy::lowMemory)),
        "Transactional policy fixture failed");
    const auto lowPolicyBytes=transactionEngine->persistentBytes();
    allocationProbe::failAfter=0;
    status=transactionEngine->setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::reuse);
    allocationProbe::failAfter=-1;
    require(status.code==WVKernelStatusCode::allocationFailure &&
        transactionEngine->persistentBytes()==lowPolicyBytes &&
        bool(transactionEngine->nonlinearFlux(state,flux)),
        "Failed reuse preparation published caches or damaged low-memory evaluation");
    require(bool(transactionEngine->setVariableEvaluationPolicy(
        WVVariableEvaluationPolicy::reuse)),
        "Transactional policy retry failed");
    // Adaptive damping consumes only the horizontal velocity pair.
    WVFrozenForcingSchedule dampingSchedule;
    dampingSchedule.entries.push_back({"WVAdaptiveDamping",1,"adaptive damping",WVForcingStage::spectral,255,0,"",{"wave-vortex-forcing-configuration-v1",1,{}}});
    std::unique_ptr<WVHydrostaticForcingEngine> dampingEngine;
    require(bool(WVHydrostaticForcingEngine::create(source,dampingSchedule,catalog,std::make_unique<WVReferenceFFTEngine>(),dampingEngine)),"Damping fixture failed");
    require(bool(dampingEngine->nonlinearFlux(state,flux)),"Damping RHS failed");
    require(dampingEngine->kernel().metrics().fieldReconstructionCount[2]==0 &&
            dampingEngine->kernel().metrics().fieldReconstructionCount[3]==0,"Horizontal damping reconstructed w or eta");
    dampingSchedule.entries.insert(dampingSchedule.entries.begin(),entry);
    require(bool(WVHydrostaticForcingEngine::create(source,dampingSchedule,catalog,std::make_unique<WVReferenceFFTEngine>(),dampingEngine)),"Combined forcing fixture failed");
    require(bool(dampingEngine->nonlinearFlux(state,flux)),"Combined forcing RHS failed");
    require(dampingEngine->kernel().metrics().stateValidationCount==1 && dampingEngine->kernel().metrics().phasePreparationCount==1,"Combined forcing duplicated state preparation");
    const auto bytes=engine->persistentBytes();
    allocationProbe::calls=0; allocationProbe::counting=true;
    for(int i=0;i<4;++i) {
        require(bool(engine->nonlinearFlux(state,flux)),"Prepared RHS failed");
        require(bool(engine->kernel().advectScalarWithAdvectionFields({scalar.data(),volume},fields,i%2,{out.data(),volume})),"Prepared tracer failed");
        WVRealFieldBundleConstView physical; require(bool(engine->physicalFields(state,physical)),"Physical fields failed");
    }
    allocationProbe::counting=false;
    require(allocationProbe::calls==0 && bytes==engine->persistentBytes(),"Prepared runtime allocated or changed storage");
    std::unique_ptr<WVIntegrationErrorPolicy> policy; require(bool(engine->createErrorPolicy(1e-10,policy)),"Error policy failed");
    require(policy->componentCount()==3 && policy->elementCount(0)==S,"Wrong adaptive families");
    for(std::size_t i=0;i<S;++i) for(std::size_t j=0;j<3;++j) require(policy->absoluteTolerance(j,i)>0 && std::isfinite(policy->absoluteTolerance(j,i)),"Invalid adaptive tolerance");
    const auto* old=engine.get();
    auto bad=schedule; bad.entries[0].typeIdentifier="unknown";
    require(!WVHydrostaticForcingEngine::create(source,bad,catalog,std::make_unique<WVReferenceFFTEngine>(),engine) && engine.get()==old,"Failed setup replaced engine");
    bool succeeded=false;
    for(long fail=0;fail<2048;++fail) {
        auto fft=std::make_unique<WVReferenceFFTEngine>(); allocationProbe::failAfter=fail;
        status=WVHydrostaticForcingEngine::create(source,schedule,catalog,std::move(fft),engine); allocationProbe::failAfter=-1;
        if(status) { succeeded=true; break; }
        require(engine.get()==old,"Failed allocation replaced engine");
    }
    require(succeeded,"Forcing allocation sweep never succeeded");
    std::unique_ptr<WVHydrostaticIntegrationSystem> system;
    status=WVHydrostaticIntegrationSystem::create(source,schedule,catalog,std::make_unique<WVReferenceFFTEngine>(),system); require(bool(status),status.message.c_str());
    require(system->stateLayout().coefficientFamilyCount()==3 && system->stateLayout().transformIdentifier()=="WVTransformHydrostatic","Wrong integration layout");
    for (const auto evaluationPolicy : {WVVariableEvaluationPolicy::reuse,
                                        WVVariableEvaluationPolicy::lowMemory}) {
        require(bool(system->setVariableEvaluationPolicy(evaluationPolicy)),
                "Hydrostatic integration policy setup failed");
        WVCoefficientStateStorage storage, denseStorage;
        require(bool(storage.initialize(system->stateLayout())) &&
                    bool(denseStorage.initialize(system->stateLayout())),
                "Hydrostatic integration lifecycle storage failed");
        for (std::size_t family=0;family<3;++family)
            std::copy_n(a[family].data(),S,
                        storage.mutableFamilies()[family].data);
        WVMutableIntegrationState integrationState;
        integrationState.waveVortex.t=83;
        integrationState.waveVortex.t0=17;
        integrationState.coefficientFamilies=storage.mutableFamilies();
        integrationState.coefficientFamilyCount=storage.familyCount();
        const auto coefficientShape=system->stateLayout().coefficientShape();
        integrationState.waveVortex.coefficients={
            {storage.mutableFamilies()[0].data,coefficientShape},
            {storage.mutableFamilies()[1].data,coefficientShape},
            {storage.mutableFamilies()[2].data,coefficientShape}};
        WVMutableIntegrationState denseState;
        denseState.coefficientFamilies=denseStorage.mutableFamilies();
        denseState.coefficientFamilyCount=denseStorage.familyCount();
        denseState.waveVortex.coefficients={
            {denseStorage.mutableFamilies()[0].data,coefficientShape},
            {denseStorage.mutableFamilies()[1].data,coefficientShape},
            {denseStorage.mutableFamilies()[2].data,coefficientShape}};
        WVFixedStepRK4 integrator(*system,{true});
        const auto restartStatus=integrator.prepareStateAfterRestart(integrationState);
        if (!restartStatus || system->variableEvaluationMetrics().liveBytes!=0)
            throw std::runtime_error(
                "Hydrostatic restart preparation failed (" +
                restartStatus.message + ") or retained " +
                std::to_string(system->variableEvaluationMetrics().liveBytes) +
                " evaluation bytes");
        const auto validationsBefore=system->kernel().metrics().stateValidationCount;
        const auto contextsBefore=system->variableEvaluationMetrics().contexts;
        require(bool(integrator.step(integrationState,1e-4)),
                "Hydrostatic integration lifecycle step failed");
        const auto validationsAfter=system->kernel().metrics().stateValidationCount;
        const auto contextsAfter=system->variableEvaluationMetrics().contexts;
        require(validationsAfter>validationsBefore &&
                    validationsAfter-validationsBefore==contextsAfter-contextsBefore &&
                    system->variableEvaluationMetrics().liveBytes==0,
                "Hydrostatic RHS validation and evaluation lifecycles diverged");
        require(bool(integrator.evaluateDenseOutput(
                    integrationState.waveVortex.t-5e-5,denseState)) &&
                    system->kernel().metrics().stateValidationCount-validationsAfter==
                        system->variableEvaluationMetrics().contexts-contextsAfter &&
                    system->variableEvaluationMetrics().liveBytes==0,
                "Hydrostatic dense output retained or reopened an RHS evaluation");
    }
    injectedServices(source, catalog, schedule);
#if WV_TEST_NATIVE_FFTW
    nativeTiledLifecycle(source,catalog,schedule);
#endif
}
}
int main() {
    try {
        Temporary file; fixture(file.path);
        { File f(file.path); for(const auto* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformHydrostatic"),"WVTransformHydrostatic")); }
        std::shared_ptr<const WVStratifiedModalRecord> source; auto status=WVStratifiedModalReader::read(file.path.string(),source); require(bool(status),status.message.c_str());
        contracts(source); std::cout<<"Hydrostatic runtime contracts passed\n";
    } catch(const std::exception& e) { allocationProbe::failAfter=-1; allocationProbe::counting=false; std::cerr<<e.what()<<'\n'; return 1; }
}
