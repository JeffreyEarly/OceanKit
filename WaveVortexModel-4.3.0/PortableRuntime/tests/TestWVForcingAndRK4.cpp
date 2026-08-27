#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVTestLinearCoefficientForcing.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

struct FailureCounter {
    std::size_t executionCount = 0;
    std::size_t failAt = 0;
};

class FailingPlan final : public WVFFTPlan {
public:
    FailingPlan(std::unique_ptr<WVFFTPlan> plan, std::shared_ptr<FailureCounter> counter) : plan_(std::move(plan)), counter_(std::move(counter)) {}
    WVKernelStatus execute(const void* input, void* output) override {
        ++counter_->executionCount;
        if (counter_->executionCount == counter_->failAt) return {WVKernelStatusCode::fftExecutionFailure,"injected portable RK4 execution failure"};
        return plan_->execute(input,output);
    }
    std::size_t persistentBytes() const noexcept override { return plan_->persistentBytes()+sizeof(*this); }
private:
    std::unique_ptr<WVFFTPlan> plan_;
    std::shared_ptr<FailureCounter> counter_;
};

class FailingEngine final : public WVFFTEngine {
public:
    explicit FailingEngine(std::shared_ptr<FailureCounter> counter) : counter_(std::move(counter)) {}
    std::string identifier() const override { return "reference-direct-injected-failure"; }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this) + reference_.persistentBytes() - sizeof(reference_) +
               sizeof(FailureCounter);
    }
    WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) override {
        std::unique_ptr<WVFFTPlan> inner;
        auto status = reference_.createPlan(specification,inner);
        if (status) plan = std::make_unique<FailingPlan>(std::move(inner),counter_);
        return status;
    }
private:
    WVReferenceFFTEngine reference_;
    std::shared_ptr<FailureCounter> counter_;
};

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void requireRelative(double actual, double expected, double tolerance,
                     const std::string &message) {
    const auto scale = std::max(std::abs(expected),1.0);
    require(std::abs(actual-expected) <= tolerance*scale,
            message+": actual="+std::to_string(actual)+
                ", expected="+std::to_string(expected));
}

WVTransformConstantStratificationConfiguration configuration(bool hydrostatic) {
    WVTransformConstantStratificationConfiguration value;
    value.Nx = 6; value.Ny = 5; value.Nz = 7; value.Nj = 4;
    value.Lx = 15000.0; value.Ly = 12000.0; value.Lz = 1300.0;
    value.N0 = 5.2e-3; value.rho0 = 1027.0; value.g = 9.80665;
    value.planetaryRadius = 6.3712e6; value.rotationRate = 7.292115e-5; value.latitude = 33.0;
    value.isHydrostatic = hydrostatic; value.shouldAntialias = true;
    return value;
}

WVPortableTypedRecord forcingConfiguration() {
    return {"wave-vortex-forcing-configuration-v1",1,{}};
}

WVPortableNamedValue realValue(std::string name, std::vector<double> values) {
    const auto size = values.size();
    return {std::move(name),size == 1 ? std::vector<std::size_t>{} : std::vector<std::size_t>{size},std::move(values)};
}

WVPortableNamedValue integerValue(std::string name, const std::vector<std::size_t>& values) {
    std::vector<std::int64_t> converted(values.begin(),values.end());
    return {std::move(name),{values.size()},std::move(converted)};
}

WVPortableTypedRecord fixedConfiguration(const WVFixedAmplitudeConfiguration& fixed) {
    auto result = forcingConfiguration();
    auto append = [&](const char* family, const auto& indices, const auto& values) {
        if (indices.empty()) return;
        std::vector<double> real(values.size()), imag(values.size());
        for (std::size_t index = 0; index < values.size(); ++index) { real[index] = values[index].real; imag[index] = values[index].imag; }
        result.values.push_back(integerValue(std::string(family)+"Indices",indices));
        result.values.push_back(realValue(std::string(family)+"ValuesReal",std::move(real)));
        result.values.push_back(realValue(std::string(family)+"ValuesImag",std::move(imag)));
    };
    append("Ap",fixed.ApIndices,fixed.ApValues); append("Am",fixed.AmIndices,fixed.AmValues); append("A0",fixed.A0Indices,fixed.A0Values);
    return result;
}

WVPortableTypedRecord pseudoConfiguration(const WVPseudoTopographicConfiguration& pseudo) {
    auto result = forcingConfiguration();
    result.values.push_back({"topographicHeight",{pseudo.topographicShape.columns,pseudo.topographicShape.rows},pseudo.topographicHeight});
    result.values.push_back(realValue("barotropicVelocityAmplitudeReal",{pseudo.barotropicVelocityAmplitude[0].real,pseudo.barotropicVelocityAmplitude[1].real}));
    result.values.push_back(realValue("barotropicVelocityAmplitudeImag",{pseudo.barotropicVelocityAmplitude[0].imag,pseudo.barotropicVelocityAmplitude[1].imag}));
    result.values.push_back(realValue("frequency",{pseudo.frequency}));
    result.values.push_back({"darwinSymbol",{},std::vector<std::string>{pseudo.darwinSymbol}});
    result.values.push_back(realValue("rampDuration",{pseudo.rampDuration}));
    result.values.push_back(realValue("startTime",{pseudo.startTime}));
    result.values.push_back({"shouldAvoidAdaptiveDamping",{},std::vector<std::uint8_t>{static_cast<std::uint8_t>(pseudo.shouldAvoidAdaptiveDamping)}});
    result.values.push_back(realValue("maximumForcedHorizontalWavenumber",{pseudo.maximumForcedHorizontalWavenumber}));
    result.values.push_back(realValue("maximumForcedVerticalMode",{pseudo.maximumForcedVerticalMode}));
    return result;
}

WVFrozenForcingEntry entry(const char* identifier, const char* name, WVForcingStage stage, std::uint8_t priority, WVPortableTypedRecord values = forcingConfiguration()) {
    return {identifier,WVPortablePairContractVersion,name,stage,priority,0,"test",std::move(values)};
}

WVFrozenForcingSchedule nonlinearSchedule() {
    WVFrozenForcingSchedule schedule;
    schedule.entries.push_back(entry("WVNonlinearAdvection","nonlinear",WVForcingStage::spatial,127));
    return schedule;
}

std::unique_ptr<WVConstantStratificationForcingEngine> createEngine(
    const WVTransformConstantStratificationConfiguration& value,
    const WVFrozenForcingSchedule& schedule) {
    std::unique_ptr<WVConstantStratificationForcingEngine> engine;
    const auto status = WVConstantStratificationForcingEngine::create(value,schedule,test::extensionCatalog(),std::make_unique<WVReferenceFFTEngine>(),engine);
    require(static_cast<bool>(status),status.message);
    return engine;
}

std::unique_ptr<WVConstantStratificationForcingEngine> createEngine(bool hydrostatic, const WVFrozenForcingSchedule& schedule) {
    return createEngine(configuration(hydrostatic),schedule);
}

std::unique_ptr<WVConstantStratificationIntegrationSystem> createSystem(
    bool hydrostatic, const WVFrozenForcingSchedule& schedule,
    std::unique_ptr<WVFFTEngine> fft =
        std::make_unique<WVReferenceFFTEngine>()) {
    std::unique_ptr<WVConstantStratificationIntegrationSystem> system;
    const auto status = WVConstantStratificationIntegrationSystem::create(
        configuration(hydrostatic), schedule, test::extensionCatalog(),
        std::move(fft), system);
    require(static_cast<bool>(status),status.message);
    return system;
}

struct OwnedState {
    WVShape2D shape;
    std::vector<WVComplex64> values;
    double t = 0.5;
    double t0 = -0.25;

    explicit OwnedState(WVShape2D valueShape) : shape(valueShape), values(3*shape.elementCount()) {
        const auto count = shape.elementCount();
        for (std::size_t index = 0; index < count; ++index) {
            values[index] = {1e-5*std::sin(0.17*static_cast<double>(index+1)),1e-5*std::cos(0.13*static_cast<double>(index+2))};
            values[count+index] = {1e-5*std::cos(0.11*static_cast<double>(index+3)),1e-5*std::sin(0.07*static_cast<double>(index+4))};
            values[2*count+index] = {1e-5*std::sin(0.19*static_cast<double>(index+5)),1e-5*std::cos(0.05*static_cast<double>(index+6))};
        }
    }

    WVMutableState mutableView() {
        const auto count = shape.elementCount();
        return {t,t0,{{values.data(),shape},{values.data()+count,shape},{values.data()+2*count,shape}}};
    }
    WVState view() const {
        const auto count = shape.elementCount();
        return {t,t0,{{values.data(),shape},{values.data()+count,shape},{values.data()+2*count,shape}}};
    }
    WVMutableIntegrationState integrationView() {
        return {mutableView(),nullptr,0};
    }
};

WVFlux fluxView(std::vector<WVComplex64>& values, WVShape2D shape) {
    const auto count = shape.elementCount();
    return {{values.data(),shape},{values.data()+count,shape},{values.data()+2*count,shape}};
}

void requireFinite(const std::vector<WVComplex64>& values, const std::string& message) {
    for (const auto value : values) require(std::isfinite(value.real) && std::isfinite(value.imag),message);
}

bool exactlyEqual(const std::vector<WVComplex64>& first, const std::vector<WVComplex64>& second) {
    if (first.size() != second.size()) return false;
    for (std::size_t index = 0; index < first.size(); ++index) {
        if (first[index].real != second[index].real || first[index].imag != second[index].imag) return false;
    }
    return true;
}

double maximumDifference(const std::vector<WVComplex64>& first, const std::vector<WVComplex64>& second) {
    require(first.size() == second.size(),"coefficient vectors have different sizes");
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result = std::max(result,std::abs(first[index].real-second[index].real));
        result = std::max(result,std::abs(first[index].imag-second[index].imag));
    }
    return result;
}

void testNonlinearCompatibility(bool hydrostatic) {
    auto scheduled = createEngine(hydrostatic,nonlinearSchedule());
    std::unique_ptr<WVTransformConstantStratificationKernel> direct;
    auto status = WVTransformConstantStratificationKernel::create(configuration(hydrostatic),std::make_unique<WVReferenceFFTEngine>(),direct);
    require(static_cast<bool>(status),"direct kernel construction failed");
    OwnedState state(direct->descriptor().spectralShape());
    const auto count = state.shape.elementCount();
    const WVComplex64 canary{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    std::vector<WVComplex64> scheduledValues(3*count,canary),directValues(3*count,canary);
    auto scheduledFlux = fluxView(scheduledValues,state.shape);
    auto directFlux = fluxView(directValues,state.shape);
    status = scheduled->nonlinearFlux(state.view(),scheduledFlux);
    require(static_cast<bool>(status),"scheduled nonlinear flux failed");
    requireFinite(scheduledValues,"integration-system RHS did not completely overwrite canary output storage");
    status = direct->nonlinearFlux(state.view(),directFlux);
    require(static_cast<bool>(status),"direct nonlinear flux failed");
    double error = 0.0;
    for (std::size_t index = 0; index < scheduledValues.size(); ++index) {
        error = std::max(error,std::abs(scheduledValues[index].real-directValues[index].real));
        error = std::max(error,std::abs(scheduledValues[index].imag-directValues[index].imag));
    }
    require(error <= 1e-13,"nonlinear-only schedule changed the shared kernel result");
}

void testForcingTrafficAccounting() {
    auto system = createSystem(true,nonlinearSchedule());
    require(system->stateLayout().additionalBlocks().empty(),"coefficient-only integration gained additional state blocks");
    require(system->fieldEvaluationService() == nullptr,"coefficient-only integration retained observing-system field scratch");
    OwnedState owned(system->kernel().descriptor().spectralShape());
    auto state = owned.integrationView();
    WVFixedStepRK4 integrator(*system);
    auto status = integrator.prepareStateAfterRestart(state);
    require(static_cast<bool>(status),"traffic-accounting restart preparation failed");
    status = integrator.step(state,0.1);
    require(static_cast<bool>(status),"traffic-accounting step failed");
    const auto& metrics = system->forcingMetrics();
    require(metrics.accumulatorClearElementWrites == 0,"nonlinear-only forcing still clears an accumulator");
    require(metrics.temporaryFluxClearElementWrites == 0 && metrics.kernelOutputInitializationElementWrites == 0,"nonlinear-only forcing still initializes an intermediate output");
    require(metrics.temporaryAccumulationElementReads == 0 && metrics.temporaryAccumulationElementWrites == 0,"nonlinear-only forcing still accumulates through temporary storage");
    require(metrics.outputCopyElementReads == 0 && metrics.outputCopyElementWrites == 0,"nonlinear-only forcing still copies its completed output");
    require(metrics.workspaceCapacityBytes == 0,"nonlinear-only forcing retained array-sized workspace");
    require(metrics.workspaceLiveBytes == 0 && metrics.workspaceMaximumLiveBytes == 0,"nonlinear-only forcing-workspace liveness is not zero");
}

void testMatlabCFLFixtures() {
    constexpr double cfl = 0.4;
    struct Fixture {
        bool hydrostatic;
        double maximumHorizontalSpeed;
        double horizontalAdvective;
        double verticalAdvective;
        double highestActiveWaveFrequency;
        double oscillatory;
    };
    const Fixture fixtures[] = {
        {false,0.1020286099012897,23522.813868795667,
         3.2273103327544224e8,0.0013924156036910925,1804.9741156372481},
        {true,0.10202859972462885,23522.816215036815,
         3.117452951358881e8,0.0014450221120718154,1739.2634354005847}};
    for (const auto &fixture : fixtures) {
        auto system = createSystem(fixture.hydrostatic,nonlinearSchedule());
        OwnedState owned(system->kernel().descriptor().spectralShape());
        const auto state = owned.integrationView();
        const auto persistentBefore = system->persistentBytes();
        WVFixedTimeStepCandidates candidates;
        const auto status = system->evaluateFixedTimeStepCandidates(
            {state.waveVortex.view(),nullptr,0},cfl,candidates);
        require(static_cast<bool>(status),
                "constant-stratification CFL evaluation: "+status.message);
        require(system->supportsFixedTimeStepSelection(),
                "constant-stratification CFL capability");
        requireRelative(candidates.effectiveHorizontalGridResolution,
                        6000.0000000000009,1e-13,
                        "MATLAB effective horizontal resolution");
        requireRelative(candidates.maximumHorizontalSpeed,
                        fixture.maximumHorizontalSpeed,1e-12,
                        "MATLAB uvMax fixture");
        requireRelative(candidates.horizontalAdvective,
                        fixture.horizontalAdvective,1e-12,
                        "MATLAB horizontal advective CFL fixture");
        requireRelative(candidates.highestActiveWaveFrequency,
                        fixture.highestActiveWaveFrequency,1e-12,
                        "MATLAB highest wave frequency fixture");
        requireRelative(candidates.oscillatory,fixture.oscillatory,1e-12,
                        "MATLAB oscillatory CFL fixture");
        requireRelative(candidates.verticalAdvective,
                        fixture.verticalAdvective,1e-12,
                        "MATLAB vertical CFL fixture");
        require(candidates.advective == candidates.horizontalAdvective,
                "advective minimum fixture");
        const auto expectedWorkspace =
            3*configuration(fixture.hydrostatic).Nx*
            configuration(fixture.hydrostatic).Ny*
            configuration(fixture.hydrostatic).Nz*sizeof(double);
        require(candidates.transientWorkspaceMaximumLiveBytes ==
                    expectedWorkspace &&
                    system->persistentBytes() == persistentBefore,
                "CFL workspace is transient and exactly accounted");
    }

    auto system = createSystem(false,nonlinearSchedule());
    OwnedState zero(system->kernel().descriptor().spectralShape());
    std::fill(zero.values.begin(),zero.values.end(),WVComplex64{});
    const auto state = zero.integrationView();
    WVFixedTimeStepCandidates candidates;
    const auto status = system->evaluateFixedTimeStepCandidates(
        {state.waveVortex.view(),nullptr,0},cfl,candidates);
    require(static_cast<bool>(status) &&
                candidates.maximumHorizontalSpeed == 0.0 &&
                std::isinf(candidates.horizontalAdvective) &&
                std::isinf(candidates.verticalAdvective) &&
                std::isfinite(candidates.oscillatory),
            "zero-velocity CFL limits");
}

void testFixedAmplitudeAndRK4() {
    WVFixedAmplitudeConfiguration fixed;
    fixed.ApIndices = {0}; fixed.ApValues = {{0.125,-0.25}};
    fixed.AmIndices = {1}; fixed.AmValues = {{-0.75,0.5}};
    fixed.A0Indices = {2}; fixed.A0Values = {{0.375,0.625}};
    WVFrozenForcingSchedule schedule;
    schedule.entries.push_back(entry("WVFixedAmplitudeForcing","fixed",WVForcingStage::spectralAmplitude,255,fixedConfiguration(fixed)));
    auto system = createSystem(true,schedule);
    OwnedState owned(system->kernel().descriptor().spectralShape());
    auto state = owned.integrationView();
    WVFixedStepRK4 integrator(*system);
    auto status = integrator.prepareStateAfterRestart(state);
    require(static_cast<bool>(status),"restart amplitude restoration failed");
    require(state.waveVortex.coefficients.Ap.data[0].real == 0.125 && state.waveVortex.coefficients.Ap.data[0].imag == -0.25,"Ap fixed amplitude was not restored");
    require(state.waveVortex.coefficients.Am.data[1].real == -0.75 && state.waveVortex.coefficients.A0.data[2].imag == 0.625,"fixed amplitudes were not restored");
    status = integrator.advanceToTime(state,0.725,0.1);
    require(static_cast<bool>(status),"fixed-amplitude RK4 integration failed");
    require(state.waveVortex.t == 0.725 && integrator.metrics().stepCount == 3 && integrator.metrics().lastStepSize > 0.024999999999 && integrator.metrics().lastStepSize < 0.025000000001,"partial final RK4 step was not deterministic");
    require(state.waveVortex.coefficients.Ap.data[0].real == 0.125 && state.waveVortex.coefficients.Am.data[1].real == -0.75 && state.waveVortex.coefficients.A0.data[2].imag == 0.625,"RK4 did not preserve fixed amplitudes");
    const auto expectedBytes = 9*state.waveVortex.coefficients.Ap.shape.elementCount()*sizeof(WVComplex64);
    require(integrator.metrics().workspaceCapacityBytes == expectedBytes,"RK4 workspace is not the bounded 9M complex schedule");
}

void testRK4DeterminismRestartAndFailure() {
    auto firstSystem = createSystem(false,nonlinearSchedule());
    auto secondSystem = createSystem(false,nonlinearSchedule());
    OwnedState first(firstSystem->kernel().descriptor().spectralShape());
    OwnedState second(secondSystem->kernel().descriptor().spectralShape());
    auto firstState = first.integrationView();
    auto secondState = second.integrationView();
    WVFixedStepRK4 firstIntegrator(*firstSystem);
    WVFixedStepRK4 secondIntegrator(*secondSystem);
    auto status = firstIntegrator.advanceToTime(firstState,0.6,0.1);
    require(static_cast<bool>(status),"first deterministic prefix integration failed");
    status = firstIntegrator.advanceToTime(firstState,0.7,0.1);
    require(static_cast<bool>(status),"first deterministic continuation failed");
    status = secondIntegrator.advanceToTime(secondState,0.6,0.1);
    require(static_cast<bool>(status),"restart prefix integration failed");
    WVFixedStepRK4 restartedIntegrator(*secondSystem);
    status = restartedIntegrator.prepareStateAfterRestart(secondState);
    require(static_cast<bool>(status),"restart preparation failed");
    status = restartedIntegrator.advanceToTime(secondState,0.7,0.1);
    require(static_cast<bool>(status),"restart continuation failed");
    if (firstState.waveVortex.t != secondState.waveVortex.t ||
        !exactlyEqual(first.values, second.values)) {
        std::ostringstream message;
        message.precision(17);
        message << "restart continuation is not bitwise equivalent to "
                   "uninterrupted fixed-step RK4; first time="
                << firstState.waveVortex.t
                << ", restarted time=" << secondState.waveVortex.t
                << ", maximum coefficient difference="
                << maximumDifference(first.values, second.values);
        throw std::runtime_error(message.str());
    }

    auto counter = std::make_shared<FailureCounter>();
    counter->failAt = 1;
    auto failingSystem = createSystem(true,nonlinearSchedule(),std::make_unique<FailingEngine>(counter));
    OwnedState failed(failingSystem->kernel().descriptor().spectralShape());
    const auto originalValues = failed.values;
    const auto originalTime = failed.t;
    auto failedState = failed.integrationView();
    WVFixedStepRK4 failingIntegrator(*failingSystem);
    status = failingIntegrator.step(failedState,0.1);
    require(status.code == WVKernelStatusCode::fftExecutionFailure,"RK4 did not preserve the injected RHS failure");
    require(exactlyEqual(failed.values,originalValues) && failedState.waveVortex.t == originalTime,"failed RK4 step modified the caller-owned state");
}

void testSpectralForcing() {
    WVFrozenForcingSchedule schedule;
    schedule.entries.push_back(entry("WVAdaptiveDamping","adaptive",WVForcingStage::spectral,100));
    schedule.entries.push_back(entry("WVBetaPlanePVAdvection","beta",WVForcingStage::spectral,110));
    auto engine = createEngine(false,schedule);
    OwnedState state(engine->kernel().descriptor().spectralShape());
    const WVComplex64 canary{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    std::vector<WVComplex64> values(3*state.shape.elementCount(),canary);
    auto flux = fluxView(values,state.shape);
    const auto status = engine->nonlinearFlux(state.view(),flux);
    require(static_cast<bool>(status),"adaptive/beta forcing failed");
    requireFinite(values,"adaptive/beta forcing produced non-finite output");
    require(engine->metrics().resolvedSpectralCount == 2,"spectral forcing dispatch was not resolved at construction");
    const auto expectedWorkspace = 4*configuration(false).Nx*configuration(false).Ny*configuration(false).Nz*sizeof(double);
    require(engine->metrics().workspaceCapacityBytes == expectedWorkspace,"spectral forcing did not allocate only its required physical-field workspace");
}

void testCoefficientErrorPolicyStorage() {
    auto engine = createEngine(true,{});
    std::unique_ptr<WVIntegrationErrorPolicy> policy;
    const auto status = engine->createErrorPolicy(1e-12,policy);
    require(static_cast<bool>(status) && policy,
            "coefficient error-policy construction failed");
    const auto coefficientCount =
        engine->kernel().descriptor().spectralShape().elementCount();
    const auto minimumPhysicalBytes =
        2*coefficientCount*sizeof(double) +
        sizeof(WVIntegrationErrorPolicy) + sizeof(WVShape2D) +
        2*sizeof(std::vector<double>);
    require(policy->persistentBytes() >= minimumPhysicalBytes,
            "coefficient error-policy ledger omitted its allocated object");
}

void testSourceLinkedLinearCoefficientExtension() {
    WVFrozenForcingSchedule schedule;
    auto values = forcingConfiguration();
    values.values.push_back(realValue("rate",{0.375}));
    schedule.entries.push_back(entry(test::LinearCoefficientForcingIdentifier,"linear",WVForcingStage::spectral,90,std::move(values)));
    auto engine = createEngine(false,schedule);
    OwnedState state(engine->kernel().descriptor().spectralShape());
    std::vector<WVComplex64> valuesOut(3*state.shape.elementCount());
    auto flux = fluxView(valuesOut,state.shape);
    const auto status = engine->nonlinearFlux(state.view(),flux);
    require(static_cast<bool>(status),"test linear forcing failed");
    const auto count = state.shape.elementCount();
    for (std::size_t index = 0; index < count; ++index) {
        require(valuesOut[index].real == 0.375*state.values[index].real && valuesOut[index].imag == 0.375*state.values[index].imag,
                "test linear forcing changed the Ap tendency formula");
        require(valuesOut[count+index].real == 0.375*state.values[count+index].real && valuesOut[count+index].imag == 0.375*state.values[count+index].imag,
                "test linear forcing changed the Am tendency formula");
        require(valuesOut[2*count+index].real == 0.375*state.values[2*count+index].real && valuesOut[2*count+index].imag == 0.375*state.values[2*count+index].imag,
                "test linear forcing changed the A0 tendency formula");
    }
}

void testLinearBottomFrictionFormula(bool hydrostatic, std::size_t Nz, double rate) {
    auto value = configuration(hydrostatic);
    value.Nz = Nz;
    value.Nj = Nz-1;
    auto drag = forcingConfiguration();
    drag.values.push_back(realValue("r",{rate}));
    WVFrozenForcingSchedule schedule;
    schedule.entries.push_back(entry("WVBottomFrictionLinear","linear bottom friction",WVForcingStage::spatial,255,std::move(drag)));
    auto engine = createEngine(value,schedule);
    OwnedState state(engine->kernel().descriptor().spectralShape());
    const auto spatial = engine->kernel().descriptor().spatialShape();
    const auto R = spatial.elementCount();
    const auto q = hydrostatic ? 3U : 4U;
    std::vector<WVComplex64> actual(3*state.shape.elementCount()),expected(actual.size());
    auto actualFlux = fluxView(actual,state.shape);
    auto status = engine->nonlinearFlux(state.view(),actualFlux);
    require(static_cast<bool>(status),"linear bottom friction failed: "+status.message);

    std::vector<double> fields(4*R),tendency(q*R);
    WVRealFieldBundleView fieldView{fields.data(),{spatial.first,spatial.second,spatial.third,4}};
    status = engine->kernel().transformWaveVortexToUVWEta(state.view(),fieldView);
    require(static_cast<bool>(status),"linear bottom-friction control reconstruction failed");
    const double scaledRate = 2.0*static_cast<double>(Nz-1)*rate;
    const auto horizontalCount = value.Nx*value.Ny;
    for (std::size_t index = 0; index < horizontalCount; ++index) {
        tendency[index] = -scaledRate*fields[index];
        tendency[R+index] = -scaledRate*fields[R+index];
    }
    auto expectedFlux = fluxView(expected,state.shape);
    WVMutableCoefficients expectedCoefficients{expectedFlux.Fp,expectedFlux.Fm,expectedFlux.F0};
    const WVRealFieldBundleConstView tendencyView{tendency.data(),{spatial.first,spatial.second,spatial.third,q}};
    status = hydrostatic ? engine->kernel().transformUVEtaToWaveVortex(tendencyView,state.t,state.t0,expectedCoefficients) : engine->kernel().transformUVWEtaToWaveVortex(tendencyView,state.t,state.t0,expectedCoefficients);
    require(static_cast<bool>(status),"linear bottom-friction control projection failed");
    require(maximumDifference(actual,expected) <= 1e-12,"linear bottom-friction tendency differs from the MATLAB formula");
    if (rate == 0.0) {
        for (const auto coefficient : actual) require(coefficient.real == 0.0 && coefficient.imag == 0.0,"zero linear drag produced a nonzero tendency");
    }
    const double expectedWeight = value.Lz/(2.0*static_cast<double>(Nz-1));
    require(engine->kernel().descriptor().bottomQuadratureWeight() == expectedWeight,"descriptor bottom quadrature weight changed");
    require(engine->kernel().descriptor().verticalModes().bottomQuadratureWeight == expectedWeight,"descriptor did not retain its bottom quadrature weight");
    require(engine->metrics().workspaceCapacityBytes == (4+q)*R*sizeof(double),"linear friction added storage beyond the shared physical and tendency fields");
    require(engine->metrics().physicalFieldReconstructionCount == 1 && engine->metrics().spatialTendencyProjectionCount == 1,"linear friction did not use one reconstruction and one generic projection");
}

void testSharedBottomFrictionOperations(bool hydrostatic) {
    WVFrozenForcingSchedule schedule;
    auto linear = forcingConfiguration(); linear.values.push_back(realValue("r",{2.5e-7}));
    auto quadratic = forcingConfiguration(); quadratic.values.push_back(realValue("Cd",{1.7e-3}));
    schedule.entries.push_back(entry("WVBottomFrictionLinear","linear",WVForcingStage::spatial,129,std::move(linear)));
    schedule.entries.push_back(entry("WVBottomFrictionQuadratic","quadratic",WVForcingStage::spatial,128,std::move(quadratic)));
    schedule.entries.push_back(entry("WVNonlinearAdvection","nonlinear",WVForcingStage::spatial,127));
    auto engine = createEngine(hydrostatic,schedule);
    OwnedState state(engine->kernel().descriptor().spectralShape());
    std::vector<WVComplex64> values(3*state.shape.elementCount());
    auto flux = fluxView(values,state.shape);
    const auto spatial = engine->kernel().descriptor().spatialShape();
    const auto R = spatial.elementCount();
    std::vector<double> fields(3*R);
    WVRealFieldBundleView fieldView{fields.data(),{spatial.first,spatial.second,spatial.third,3}};
    WVConstantStratificationRightHandSideContext context;
    const auto status = engine->evaluateRightHandSideWithContext(state.view(),flux,fieldView,context);
    require(static_cast<bool>(status),"ordered shared bottom-friction evaluation failed: "+status.message);
    requireFinite(values,"shared bottom-friction evaluation produced non-finite output");
    require(engine->scheduleIdentifier() == "wave-vortex-forcing-v1:WVNonlinearAdvection,WVBottomFrictionQuadratic,WVBottomFrictionLinear","stage/priority ordering did not select the exact active pairs");
    require(engine->metrics().physicalFieldReconstructionCount == 1,"nonlinear and bottom-friction forcings redundantly reconstructed physical fields");
    require(engine->metrics().physicalFieldReuseCount >= 2,"bottom-friction forcings did not reuse nonlinear-advection physical fields");
    require(engine->metrics().spatialTendencyProjectionCount == 2,"bottom-friction forcings did not share generic projection");
    const auto q = hydrostatic ? 3U : 4U;
    require(engine->metrics().spatialTendencyClearElementWrites == 2*q*R,"bottom-friction forcings did not clear the shared tendency exactly once each");
}

void testLinearBottomFrictionIntegrationAndFailures(bool hydrostatic) {
    auto drag = forcingConfiguration(); drag.values.push_back(realValue("r",{2.5e-7}));
    auto schedule = nonlinearSchedule();
    schedule.entries.push_back(entry("WVBottomFrictionLinear","linear",WVForcingStage::spatial,255,std::move(drag)));
    auto fixedSystem = createSystem(hydrostatic,schedule);
    OwnedState fixedState(fixedSystem->kernel().descriptor().spectralShape());
    auto fixedView = fixedState.integrationView();
    WVFixedStepRK4 fixed(*fixedSystem);
    auto status = fixed.prepareStateAfterRestart(fixedView);
    require(static_cast<bool>(status),"linear-friction RK4 restart preparation failed");
    status = fixed.advanceToTime(fixedView,0.52,0.01);
    require(static_cast<bool>(status),"linear-friction RK4 continuation failed");
    requireFinite(fixedState.values,"linear-friction RK4 produced non-finite state");

    auto adaptiveSystem = createSystem(hydrostatic,schedule);
    OwnedState adaptiveState(adaptiveSystem->kernel().descriptor().spectralShape());
    auto adaptiveView = adaptiveState.integrationView();
    WVAdaptiveRK23Options options;
    options.relativeTolerance = 1e-7;
    options.absoluteToleranceScale = 1e-10;
    options.maximumStepSize = 0.01;
    WVAdaptiveRK23 adaptive(*adaptiveSystem,options);
    status = adaptive.prepareStateAfterRestart(adaptiveView);
    require(static_cast<bool>(status),"linear-friction RK23 restart preparation failed");
    status = adaptive.advanceToTime(adaptiveView,0.52,0.01);
    require(static_cast<bool>(status),"linear-friction RK23 continuation failed");
    requireFinite(adaptiveState.values,"linear-friction RK23 produced non-finite state");

    auto invalidSchedule = schedule;
    invalidSchedule.entries.back().configuration.values.front() = realValue("r",{-1.0});
    std::unique_ptr<WVConstantStratificationForcingEngine> invalidEngine;
    status = WVConstantStratificationForcingEngine::create(configuration(hydrostatic),invalidSchedule,test::extensionCatalog(),std::make_unique<WVReferenceFFTEngine>(),invalidEngine);
    require(!status && !invalidEngine,"negative linear drag was accepted");
    invalidSchedule = schedule;
    invalidSchedule.entries.back().contractVersion = WVPortablePairContractVersion+1;
    status = WVConstantStratificationForcingEngine::create(configuration(hydrostatic),invalidSchedule,test::extensionCatalog(),std::make_unique<WVReferenceFFTEngine>(),invalidEngine);
    require(!status && !invalidEngine,"linear drag contract-version mismatch was accepted");

    WVFrozenForcingSchedule failureSchedule;
    auto failureDrag = forcingConfiguration(); failureDrag.values.push_back(realValue("r",{2.5e-7}));
    failureSchedule.entries.push_back(entry("WVBottomFrictionLinear","linear",WVForcingStage::spatial,255,std::move(failureDrag)));
    auto counter = std::make_shared<FailureCounter>();
    counter->failAt = 1;
    auto failingSystem = createSystem(hydrostatic,failureSchedule,std::make_unique<FailingEngine>(counter));
    OwnedState failed(failingSystem->kernel().descriptor().spectralShape());
    const auto originalValues = failed.values;
    const auto originalTime = failed.t;
    auto failedView = failed.integrationView();
    WVFixedStepRK4 failingIntegrator(*failingSystem);
    status = failingIntegrator.step(failedView,0.01);
    require(status.code == WVKernelStatusCode::fftExecutionFailure,"linear-friction reconstruction failure was not preserved");
    require(exactlyEqual(failed.values,originalValues) && failedView.waveVortex.t == originalTime,"failed linear-friction step modified accepted state");
}

void testQuadraticAndPseudo(bool hydrostatic) {
    WVPseudoTopographicConfiguration pseudo;
    pseudo.topographicShape = {6,5};
    pseudo.topographicHeight.resize(30);
    for (std::size_t index = 0; index < pseudo.topographicHeight.size(); ++index) pseudo.topographicHeight[index] = 10.0*std::sin(0.2*static_cast<double>(index+1));
    pseudo.barotropicVelocityAmplitude = {WVComplex64{0.12,0.01},WVComplex64{-0.04,0.02}};
    pseudo.frequency = 1.405e-4; pseudo.rampDuration = 900.0; pseudo.startTime = -50.0;
    pseudo.maximumForcedHorizontalWavenumber = std::numeric_limits<double>::infinity();
    pseudo.maximumForcedVerticalMode = std::numeric_limits<double>::infinity();
    WVFrozenForcingSchedule schedule;
    auto drag = forcingConfiguration(); drag.values.push_back(realValue("Cd",{1.7e-3}));
    schedule.entries.push_back(entry("WVBottomFrictionQuadratic","drag",WVForcingStage::spatial,128,std::move(drag)));
    schedule.entries.push_back(entry("WVPseudoTopographicWaveGeneration","topography",WVForcingStage::spectral,127,pseudoConfiguration(pseudo)));
    auto engine = createEngine(hydrostatic,schedule);
    OwnedState state(engine->kernel().descriptor().spectralShape());
    const WVComplex64 canary{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    std::vector<WVComplex64> values(3*state.shape.elementCount(),canary);
    auto flux = fluxView(values,state.shape);
    const auto status = engine->nonlinearFlux(state.view(),flux);
    require(static_cast<bool>(status),"quadratic/pseudo forcing failed: "+status.message);
    requireFinite(values,"quadratic/pseudo forcing produced non-finite output");
    const auto& value = configuration(hydrostatic);
    const auto R = value.Nx*value.Ny*value.Nz;
    const auto q = hydrostatic ? 3U : 4U;
    require(engine->metrics().workspaceCapacityBytes == (4+q)*R*sizeof(double),"quadratic forcing did not allocate only its required real-field workspace");
}

void testMultipleWholeFluxProducers() {
    auto schedule = nonlinearSchedule();
    schedule.entries.push_back(entry("WVNonlinearAdvection","nonlinear-second",WVForcingStage::spatial,128));
    auto engine = createEngine(true,schedule);
    auto direct = createEngine(true,nonlinearSchedule());
    OwnedState state(engine->kernel().descriptor().spectralShape());
    const WVComplex64 canary{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    std::vector<WVComplex64> values(3*state.shape.elementCount(),canary),single(values.size(),canary);
    auto flux = fluxView(values,state.shape);
    auto singleFlux = fluxView(single,state.shape);
    auto status = engine->nonlinearFlux(state.view(),flux);
    require(static_cast<bool>(status),"multiple whole-flux producers failed");
    status = direct->nonlinearFlux(state.view(),singleFlux);
    require(static_cast<bool>(status),"single whole-flux control failed");
    for (std::size_t index = 0; index < values.size(); ++index) {
        require(std::abs(values[index].real-2.0*single[index].real) <= 1e-13 && std::abs(values[index].imag-2.0*single[index].imag) <= 1e-13,"temporary whole-flux accumulation changed the result");
    }
    require(engine->metrics().workspaceCapacityBytes == values.size()*sizeof(WVComplex64),"multiple whole-flux producers did not allocate exactly one temporary tendency");
}

void testRightHandSideContextIdentity() {
    auto first = createEngine(true,nonlinearSchedule());
    auto second = createEngine(true,nonlinearSchedule());
    OwnedState state(first->kernel().descriptor().spectralShape());
    std::vector<WVComplex64> values(3*state.shape.elementCount());
    auto flux = fluxView(values,state.shape);
    const auto spatial = first->kernel().descriptor().spatialShape();
    const auto R = spatial.elementCount();
    std::vector<double> fields(3*R),scalar(R,1.0),scalarFlux(R);
    WVRealFieldBundleView fieldView{fields.data(),{spatial.first,spatial.second,spatial.third,3}};
    WVConstantStratificationRightHandSideContext context;
    auto status = first->evaluateRightHandSideWithContext(state.view(),flux,fieldView,context);
    require(static_cast<bool>(status) && context.hasAdvectionFields(),"RHS context was not prepared");
    const WVRealVolumeConstView scalarView{scalar.data(),spatial};
    WVRealVolumeView scalarFluxView{scalarFlux.data(),spatial};
    status = second->advectFGridScalar(context,scalarView,false,scalarFluxView);
    require(status.code == WVKernelStatusCode::invalidConfiguration,"a foreign RHS context was accepted");
    WVConstantStratificationRightHandSideContext replacement;
    status = first->evaluateRightHandSideWithContext(state.view(),flux,fieldView,replacement);
    require(static_cast<bool>(status),"replacement RHS context failed");
    status = first->advectFGridScalar(context,scalarView,false,scalarFluxView);
    require(status.code == WVKernelStatusCode::invalidConfiguration,"a stale RHS context was accepted");
    status = first->advectFGridScalar(replacement,scalarView,false,scalarFluxView);
    require(static_cast<bool>(status),"current RHS context was rejected");
}

void testValidation() {
    auto schedule = nonlinearSchedule();
    schedule.profileVersion = 99;
    std::unique_ptr<WVConstantStratificationForcingEngine> engine;
    auto status = WVConstantStratificationForcingEngine::create(configuration(true),schedule,test::extensionCatalog(),std::make_unique<WVReferenceFFTEngine>(),engine);
    require(status.code == WVKernelStatusCode::unsupportedOperation && !engine,"unsupported schedule profile was accepted");
    schedule = nonlinearSchedule();
    schedule.entries.front().stage = WVForcingStage::spectral;
    status = WVConstantStratificationForcingEngine::create(configuration(true),schedule,test::extensionCatalog(),std::make_unique<WVReferenceFFTEngine>(),engine);
    require(status.code == WVKernelStatusCode::invalidConfiguration && !engine,"wrong forcing stage was accepted");
}

} // namespace

int main() {
    try {
        testNonlinearCompatibility(true);
        testNonlinearCompatibility(false);
        testFixedAmplitudeAndRK4();
        testForcingTrafficAccounting();
        testMatlabCFLFixtures();
        testRK4DeterminismRestartAndFailure();
        testSpectralForcing();
        testCoefficientErrorPolicyStorage();
        testSourceLinkedLinearCoefficientExtension();
        testLinearBottomFrictionFormula(true,5,0.0);
        testLinearBottomFrictionFormula(false,5,2.5e-7);
        testLinearBottomFrictionFormula(true,9,2.5e-7);
        testLinearBottomFrictionFormula(false,9,2.5e-7);
        testSharedBottomFrictionOperations(true);
        testSharedBottomFrictionOperations(false);
        testLinearBottomFrictionIntegrationAndFailures(true);
        testLinearBottomFrictionIntegrationAndFailures(false);
        testQuadraticAndPseudo(true);
        testQuadraticAndPseudo(false);
        testMultipleWholeFluxProducers();
        testRightHandSideContextIdentity();
        testValidation();
        std::cout << "Portable forcing and RK4 tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
