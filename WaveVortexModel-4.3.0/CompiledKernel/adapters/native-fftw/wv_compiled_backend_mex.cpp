#include "mex.h"

#include "WVNativeFFTWEngine.hpp"
#include "WVModelInternalAccess.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"

#include <cstdint>
#include <chrono>
#include <memory>
#include <string>
#include <unordered_map>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

std::unordered_map<std::uint64_t,std::unique_ptr<WVModel>> models;
std::uint64_t nextHandle = 1;

void fail(const char* identifier, const std::string& message) { mexErrMsgIdAndTxt(identifier,"%s",message.c_str()); }

double scalarField(const mxArray* value, const char* name) {
    const mxArray* field = mxGetField(value,0,name);
    if (field == nullptr || (!mxIsNumeric(field) && !mxIsLogical(field)) || mxIsComplex(field) || mxGetNumberOfElements(field) != 1) fail("WaveVortexModel:CompiledKernelConfiguration",std::string("Missing scalar configuration field ") + name + '.');
    return mxGetScalar(field);
}
WVTransformConstantStratificationConfiguration configuration(const mxArray* value) {
    if (!mxIsStruct(value) || mxGetNumberOfElements(value) != 1) fail("WaveVortexModel:CompiledKernelConfiguration","Configuration must be a scalar struct.");
    WVTransformConstantStratificationConfiguration result;
    result.Nx = static_cast<std::size_t>(scalarField(value,"Nx")); result.Ny = static_cast<std::size_t>(scalarField(value,"Ny"));
    result.Nz = static_cast<std::size_t>(scalarField(value,"Nz")); result.Nj = static_cast<std::size_t>(scalarField(value,"Nj"));
    result.Lx = scalarField(value,"Lx"); result.Ly = scalarField(value,"Ly"); result.Lz = scalarField(value,"Lz");
    result.N0 = scalarField(value,"N0"); result.rho0 = scalarField(value,"rho0"); result.g = scalarField(value,"g");
    result.planetaryRadius = scalarField(value,"planetaryRadius"); result.rotationRate = scalarField(value,"rotationRate"); result.latitude = scalarField(value,"latitude");
    result.isHydrostatic = scalarField(value,"isHydrostatic") != 0.0; result.shouldAntialias = scalarField(value,"shouldAntialias") != 0.0;
    return result;
}

std::uint64_t scalarHandle(const mxArray* value) {
    if (!mxIsUint64(value) || mxGetNumberOfElements(value) != 1) fail("WaveVortexModel:CompiledKernelHandle","Kernel handle must be a uint64 scalar.");
    return *mxGetUint64s(value);
}

WVTransformConstantStratificationKernel& kernel(const mxArray* value) {
    const auto handle = scalarHandle(value);
    const auto iterator = models.find(handle);
    if (iterator == models.end()) fail("WaveVortexModel:CompiledKernelHandle","Kernel handle is invalid or deleted.");
    return detail::WVModelInternalAccess::constantStratificationKernel(*iterator->second);
}

WVModel& model(const mxArray* value) {
    const auto handle = scalarHandle(value);
    const auto iterator = models.find(handle);
    if (iterator == models.end()) fail("WaveVortexModel:CompiledKernelHandle","Kernel handle is invalid or deleted.");
    return *iterator->second;
}

const char* statusIdentifier(WVKernelStatusCode code) {
    switch (code) {
        case WVKernelStatusCode::invalidConfiguration: return "WaveVortexModel:CompiledKernelInvalidConfiguration";
        case WVKernelStatusCode::invalidShape: return "WaveVortexModel:CompiledKernelShape";
        case WVKernelStatusCode::invalidPointer:
        case WVKernelStatusCode::overlappingArrays: return "WaveVortexModel:CompiledKernelOwnership";
        case WVKernelStatusCode::allocationFailure: return "WaveVortexModel:CompiledKernelAllocation";
        case WVKernelStatusCode::fftPlanFailure: return "WaveVortexModel:CompiledKernelPlan";
        case WVKernelStatusCode::fftExecutionFailure: return "WaveVortexModel:CompiledKernelExecution";
        case WVKernelStatusCode::numericalFailure: return "WaveVortexModel:CompiledKernelNumericalFailure";
        case WVKernelStatusCode::unsupportedOperation: return "WaveVortexModel:CompiledKernelUnsupportedOperation";
        case WVKernelStatusCode::reentrantExecution: return "WaveVortexModel:CompiledKernelReentrantExecution";
        case WVKernelStatusCode::sizeOverflow: return "WaveVortexModel:CompiledKernelSizeOverflow";
        case WVKernelStatusCode::success: break;
    }
    return "WaveVortexModel:CompiledKernelFailure";
}

void requireStatus(const WVKernelStatus& status) { if (!status) fail(statusIdentifier(status.code),status.message); }

std::shared_ptr<const WVExtensionCatalog> makeModelCatalog() {
    std::shared_ptr<const WVExtensionCatalog> catalog;
    requireStatus(makeBuiltInExtensionCatalog(catalog));
    return catalog;
}

std::string stringInput(const mxArray* value, const char* name) {
    if (!mxIsChar(value)) fail("WaveVortexModel:CompiledKernelCommand",std::string(name) + " must be a character vector.");
    char buffer[64];
    if (mxGetString(value,buffer,sizeof(buffer)) != 0) fail("WaveVortexModel:CompiledKernelCommand",std::string(name) + " is too long.");
    return buffer;
}

std::string arbitraryLengthStringInput(const mxArray* value, const char* name) {
    if (!mxIsChar(value)) fail("WaveVortexModel:CompiledKernelCommand",std::string(name) + " must be a character vector.");
    char* buffer = mxArrayToString(value);
    if (buffer == nullptr) fail("WaveVortexModel:CompiledKernelCommand",std::string("Unable to read ") + name + '.');
    std::string result(buffer);
    mxFree(buffer);
    return result;
}

class InjectedFailureEngine final : public WVFFTEngine {
public:
    explicit InjectedFailureEngine(WVKernelStatusCode code) : code_(code) {}
    std::string identifier() const override { return "injected-failure"; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    WVKernelStatus createPlan(const WVFFTPlanSpecification&, std::unique_ptr<WVFFTPlan>&) override {
        return {code_,code_ == WVKernelStatusCode::allocationFailure ? "Injected allocation failure." : "Injected plan-creation failure."};
    }
private:
    WVKernelStatusCode code_;
};

class InjectedExecutionPlan final : public WVFFTPlan {
public:
    InjectedExecutionPlan(std::unique_ptr<WVFFTPlan> plan, std::shared_ptr<bool> shouldFail) : plan_(std::move(plan)), shouldFail_(std::move(shouldFail)) {}
    WVKernelStatus execute(const void* input, void* output) override {
        if (*shouldFail_) { *shouldFail_ = false; return {WVKernelStatusCode::fftExecutionFailure,"Injected FFT execution failure."}; }
        return plan_->execute(input,output);
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this) + plan_->persistentBytes(); }
private:
    std::unique_ptr<WVFFTPlan> plan_;
    std::shared_ptr<bool> shouldFail_;
};

class InjectedExecutionEngine final : public WVFFTEngine {
public:
    explicit InjectedExecutionEngine(std::unique_ptr<WVFFTEngine> engine) : engine_(std::move(engine)), shouldFail_(std::make_shared<bool>(true)) {}
    std::string identifier() const override { return "fftw-injected-execution"; }
    std::string libraryIdentity() const override { return engine_->libraryIdentity(); }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this) + engine_->persistentBytes() + sizeof(bool);
    }
    WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) override {
        std::unique_ptr<WVFFTPlan> inner;
        auto status = engine_->createPlan(specification,inner);
        if (!status) return status;
        try { plan = std::make_unique<InjectedExecutionPlan>(std::move(inner),shouldFail_); }
        catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Unable to allocate injected execution plan."}; }
        return WVKernelStatus::ok();
    }
private:
    std::unique_ptr<WVFFTEngine> engine_;
    std::shared_ptr<bool> shouldFail_;
};

std::unique_ptr<WVFFTEngine> injectedEngine(const std::string& mode, std::size_t threads) {
    if (mode == "plan") return std::make_unique<InjectedFailureEngine>(WVKernelStatusCode::fftPlanFailure);
    if (mode == "allocation") return std::make_unique<InjectedFailureEngine>(WVKernelStatusCode::allocationFailure);
    if (mode == "execution") {
        std::unique_ptr<WVFFTEngine> engine;
        requireStatus(WVFFTWEngine::create(threads,engine));
        return std::make_unique<InjectedExecutionEngine>(std::move(engine));
    }
    fail("WaveVortexModel:CompiledKernelCommand","Unknown injected failure mode.");
    return {};
}

WVComplexConstView complexInput(const mxArray* value, WVShape2D shape, const char* name) {
    if (!mxIsDouble(value) || !mxIsComplex(value) || mxGetM(value) != shape.rows || mxGetN(value) != shape.columns) fail("WaveVortexModel:CompiledKernelShape",std::string(name) + " must be a complex [Nj,Nkl] array.");
    return {reinterpret_cast<const WVComplex64*>(mxGetComplexDoubles(value)),shape};
}

mxArray* complexOutput(WVShape2D shape, WVComplexView& view) {
    mxArray* value = mxCreateDoubleMatrix(shape.rows,shape.columns,mxCOMPLEX);
    view = {reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(value)),shape};
    return value;
}

WVRealFieldBundleConstView realBundleInput(const mxArray* value, WVShape3D spatial, std::size_t channels) {
    const auto dimensions = mxGetDimensions(value);
    const auto count = mxGetNumberOfDimensions(value);
    if (!mxIsDouble(value) || mxIsComplex(value) || count != 4 || dimensions[0] != spatial.first || dimensions[1] != spatial.second || dimensions[2] != spatial.third || dimensions[3] != channels) fail("WaveVortexModel:CompiledKernelShape","Spatial fields must be a real [Nx,Ny,Nz,Nfield] array.");
    return {mxGetDoubles(value),{spatial.first,spatial.second,spatial.third,channels}};
}

mxArray* realBundleOutput(WVShape3D spatial, std::size_t channels, WVRealFieldBundleView& view) {
    const mwSize dimensions[] = {spatial.first,spatial.second,spatial.third,channels};
    mxArray* value = mxCreateNumericArray(4,dimensions,mxDOUBLE_CLASS,mxREAL);
    view = {mxGetDoubles(value),{spatial.first,spatial.second,spatial.third,channels}};
    return value;
}

void cleanup() { models.clear(); }

mxArray* scalarString(const std::string& value) { return mxCreateString(value.c_str()); }

mxArray* moduleInfo(const std::string& expectedOpenMPRuntime) {
    const char* names[] = {"engine","version","baseLibrary","threadLibrary","openMPRuntimeLibrary"};
    mxArray* result = mxCreateStructMatrix(1,1,5,names);
    const auto identity = WVFFTWEngine::linkedLibraries(expectedOpenMPRuntime);
    mxSetField(result,0,"engine",mxCreateString("fftw"));
    mxSetField(result,0,"version",scalarString(identity.version));
    mxSetField(result,0,"baseLibrary",scalarString(identity.baseLibrary));
    mxSetField(result,0,"threadLibrary",scalarString(identity.threadLibrary));
    mxSetField(result,0,"openMPRuntimeLibrary",scalarString(identity.openMPRuntimeLibrary));
    return result;
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1 || !mxIsChar(prhs[0])) fail("WaveVortexModel:CompiledKernelCommand","The first input must be a command string.");
    char commandBuffer[64];
    if (mxGetString(prhs[0],commandBuffer,sizeof(commandBuffer)) != 0) fail("WaveVortexModel:CompiledKernelCommand","Command is too long.");
    const std::string command(commandBuffer);
    if (command == "moduleInfo") {
        if ((nrhs != 1 && nrhs != 2) || nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","moduleInfo accepts an optional expected OpenMP runtime path.");
        const auto expectedRuntime = nrhs == 2 ? arbitraryLengthStringInput(prhs[1],"Expected OpenMP runtime") : std::string{};
        plhs[0] = moduleInfo(expectedRuntime);
        return;
    }
    if (command == "moduleMetrics") {
        if (nrhs != 1 || nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","moduleMetrics takes no additional inputs.");
        const char* names[] = {"kernelCount","moduleLocked","activePlans","totalPlansCreated","totalPlansDestroyed","outstandingPlanningBytes","totalPlanningSeconds"};
        plhs[0] = mxCreateStructMatrix(1,1,7,names);
        const auto lifetime = WVFFTWEngine::lifetimeMetrics();
        const double values[] = {static_cast<double>(models.size()),models.empty() ? 0.0 : 1.0,static_cast<double>(lifetime.activePlans),static_cast<double>(lifetime.totalPlansCreated),static_cast<double>(lifetime.totalPlansDestroyed),static_cast<double>(lifetime.outstandingPlanningBytes),lifetime.totalPlanningSeconds};
        for (std::size_t i = 0; i < 7; ++i) mxSetField(plhs[0],0,names[i],mxCreateDoubleScalar(values[i]));
        return;
    }
    if (command == "estimate") {
        if (nrhs != 2 || nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","estimate requires one configuration input.");
        WVTransformConstantStratificationDescriptor descriptor;
        requireStatus(WVTransformConstantStratificationDescriptor::create(configuration(prhs[1]),descriptor));
        const auto& c = descriptor.configuration();
        const auto halfElements = descriptor.halfSpectrumMappings().NxHalf * c.Ny * c.Nz * 4;
        const auto halfSpectrumScratchCapacityBytes = 2 * halfElements * sizeof(double);
        const auto realScratchCapacityBytes = descriptor.spatialShape().elementCount() * 6 * sizeof(double);
        const auto scratchCapacityBytes = halfSpectrumScratchCapacityBytes + realScratchCapacityBytes;
        const auto spectralBytes = descriptor.spectralShape().elementCount() * sizeof(WVComplex64);
        const auto descriptorBytes = descriptor.persistentBytes() +
                                     descriptor.halfSpectrumMappings().NxHalf *
                                         c.Ny * sizeof(std::uint8_t);
        const auto persistentBytesLowerBound = descriptorBytes + scratchCapacityBytes;
        const char* names[] = {"contractVersion","planCount","planMemoryAccounting","descriptorBytes","halfSpectrumScratchCapacityBytes","realScratchCapacityBytes","scratchCapacityBytes","persistentBytesLowerBound","stateInputBytes","fluxOutputBytes","knownMaximumLiveOwnedBytesLowerBound","Nx","Ny","Nz","Nj","Nkl"};
        plhs[0] = mxCreateStructMatrix(1,1,16,names);
        mxSetField(plhs[0],0,"contractVersion",mxCreateDoubleScalar(static_cast<double>(WVKernelContractVersion)));
        mxSetField(plhs[0],0,"planCount",mxCreateDoubleScalar(17.0));
        mxSetField(plhs[0],0,"planMemoryAccounting",mxCreateString("FFTW plan storage is opaque before construction"));
        const double values[] = {static_cast<double>(descriptorBytes),static_cast<double>(halfSpectrumScratchCapacityBytes),static_cast<double>(realScratchCapacityBytes),static_cast<double>(scratchCapacityBytes),static_cast<double>(persistentBytesLowerBound),static_cast<double>(3*spectralBytes),static_cast<double>(3*spectralBytes),static_cast<double>(persistentBytesLowerBound+3*spectralBytes),static_cast<double>(c.Nx),static_cast<double>(c.Ny),static_cast<double>(c.Nz),static_cast<double>(c.Nj),static_cast<double>(descriptor.Nkl())};
        for (std::size_t i = 0; i < 13; ++i) mxSetField(plhs[0],0,names[i+3],mxCreateDoubleScalar(values[i]));
        return;
    }
    if (command == "create") {
        if (nrhs != 3 || nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","create requires configuration and thread count.");
        std::unique_ptr<WVFFTEngine> engine;
        requireStatus(WVFFTWEngine::create(static_cast<std::size_t>(mxGetScalar(prhs[2])),engine));
        auto value = std::make_unique<WVModel>();
        requireStatus(WVModel::create(makeModelCatalog(),configuration(prhs[1]),
                                      defaultNonlinearAdvectionSchedule(),
                                      std::move(engine),{},*value));
        const auto handle = nextHandle++;
        models.emplace(handle,std::move(value));
        if (models.size() == 1) { mexLock(); mexAtExit(cleanup); }
        plhs[0] = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL); *mxGetUint64s(plhs[0]) = handle;
        return;
    }
    if (command == "createInjectedFailure") {
        if (nrhs != 4 || nlhs > 1) fail("WaveVortexModel:CompiledKernelCommand","createInjectedFailure requires configuration, thread count, and failure mode.");
        const auto mode = stringInput(prhs[3],"Failure mode");
        if (mode == "execution" && nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","The execution failure mode returns a handle.");
        auto engine = injectedEngine(mode,static_cast<std::size_t>(mxGetScalar(prhs[2])));
        auto value = std::make_unique<WVModel>();
        requireStatus(WVModel::create(makeModelCatalog(),configuration(prhs[1]),
                                      defaultNonlinearAdvectionSchedule(),
                                      std::move(engine),{},*value));
        const auto handle = nextHandle++;
        models.emplace(handle,std::move(value));
        if (models.size() == 1) { mexLock(); mexAtExit(cleanup); }
        plhs[0] = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL); *mxGetUint64s(plhs[0]) = handle;
        return;
    }
    if (command == "delete") {
        if (nrhs != 2 || nlhs != 0) fail("WaveVortexModel:CompiledKernelCommand","delete requires one handle.");
        if (models.erase(scalarHandle(prhs[1])) == 0) fail("WaveVortexModel:CompiledKernelHandle","Kernel handle is invalid or deleted.");
        if (models.empty()) mexUnlock();
        return;
    }
    auto& value = kernel(nrhs > 1 ? prhs[1] : nullptr);
    const auto spectral = value.descriptor().spectralShape();
    const auto spatial = value.descriptor().spatialShape();
    if (command == "forward" || command == "forwardTimed") {
        const bool timed = command == "forwardTimed";
        if (nrhs != 5 || nlhs != (timed ? 4 : 3)) fail("WaveVortexModel:CompiledKernelCommand","forward requires handle, fields, t, and t0 and returns Ap, Am, A0.");
        const auto channels = value.descriptor().configuration().isHydrostatic ? 3U : 4U;
        const auto fields = realBundleInput(prhs[2],spatial,channels);
        WVMutableCoefficients coefficients;
        plhs[0] = complexOutput(spectral,coefficients.Ap); plhs[1] = complexOutput(spectral,coefficients.Am); plhs[2] = complexOutput(spectral,coefficients.A0);
        const auto start = std::chrono::steady_clock::now();
        requireStatus(value.descriptor().configuration().isHydrostatic ? value.transformUVEtaToWaveVortex(fields,mxGetScalar(prhs[3]),mxGetScalar(prhs[4]),coefficients) : value.transformUVWEtaToWaveVortex(fields,mxGetScalar(prhs[3]),mxGetScalar(prhs[4]),coefficients));
        if (timed) plhs[3] = mxCreateDoubleScalar(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
        return;
    }
    if (command == "inverse" || command == "inverseTimed") {
        const bool timed = command == "inverseTimed";
        if (nrhs != 7 || nlhs != (timed ? 2 : 1)) fail("WaveVortexModel:CompiledKernelCommand","inverse requires handle, Ap, Am, A0, t, and t0.");
        WVState state{mxGetScalar(prhs[5]),mxGetScalar(prhs[6]),{complexInput(prhs[2],spectral,"Ap"),complexInput(prhs[3],spectral,"Am"),complexInput(prhs[4],spectral,"A0")}};
        WVRealFieldBundleView fields; plhs[0] = realBundleOutput(spatial,4,fields);
        const auto start = std::chrono::steady_clock::now();
        requireStatus(value.transformWaveVortexToUVWEta(state,fields));
        if (timed) plhs[1] = mxCreateDoubleScalar(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
        return;
    }
    if (command == "fAll" || command == "gAll" || command == "fAllTimed" || command == "gAllTimed") {
        const bool timed = command == "fAllTimed" || command == "gAllTimed";
        const bool isF = command == "fAll" || command == "fAllTimed";
        if (nrhs != 4 || nlhs != (timed ? 2 : 1)) fail("WaveVortexModel:CompiledKernelCommand","fAll/gAll require handle, Apm, and A0.");
        WVRealFieldBundleView fields; plhs[0] = realBundleOutput(spatial,4,fields);
        const auto Apm = complexInput(prhs[2],spectral,"Apm"); const auto A0 = complexInput(prhs[3],spectral,"A0");
        const auto start = std::chrono::steady_clock::now();
        requireStatus(isF ? value.transformToSpatialDomainWithFAllDerivatives(Apm,A0,fields) : value.transformToSpatialDomainWithGAllDerivatives(Apm,A0,fields));
        if (timed) plhs[1] = mxCreateDoubleScalar(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
        return;
    }
    if (command == "nonlinearFlux" || command == "nonlinearFluxTimed") {
        const bool timed = command == "nonlinearFluxTimed";
        if (nrhs != 7 || nlhs != (timed ? 4 : 3)) fail("WaveVortexModel:CompiledKernelCommand","nonlinearFlux requires handle, Ap, Am, A0, t, and t0 and returns Fp, Fm, F0.");
        WVState state{mxGetScalar(prhs[5]),mxGetScalar(prhs[6]),{complexInput(prhs[2],spectral,"Ap"),complexInput(prhs[3],spectral,"Am"),complexInput(prhs[4],spectral,"A0")}};
        WVFlux flux;
        plhs[0] = complexOutput(spectral,flux.Fp); plhs[1] = complexOutput(spectral,flux.Fm); plhs[2] = complexOutput(spectral,flux.F0);
        const auto start = std::chrono::steady_clock::now();
        WVIntegrationState integrationState{state,nullptr,0};
        WVIntegrationFlux integrationFlux{flux,nullptr,0};
        requireStatus(model(prhs[1]).evaluateRightHandSide(
            integrationState,integrationFlux));
        if (timed) plhs[3] = mxCreateDoubleScalar(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
        return;
    }
    if (command == "setStageInstrumentation") {
        if (nrhs != 3 || nlhs != 0 || !mxIsLogicalScalar(prhs[2])) fail("WaveVortexModel:CompiledKernelCommand","setStageInstrumentation requires handle and a logical scalar.");
        value.setStageInstrumentation(mxIsLogicalScalarTrue(prhs[2]));
        return;
    }
    if (command == "metrics") {
        if (nrhs != 2 || nlhs != 1) fail("WaveVortexModel:CompiledKernelCommand","metrics requires one handle.");
        const char* names[] = {"engine","loadedLibrary","nonlinearFluxSchedule","phaseImplementation","coefficientStorageMode","coefficientArithmeticMode","inverseNormalizationPlacement","optimizationImplementation","coefficientWorkerCount","planMemoryAccounting","contractVersion","planCount","planBytes","engineBytes","kernelManagementBytes","descriptorBytes","scratchCapacityBytes","scratchHighWaterBytes","halfSpectrumScratchCapacityBytes","realScratchCapacityBytes","executionCount","horizontalExecutionCount","verticalExecutionCount","nonlinearFluxCallCount","nonlinearFluxPhaseEvaluationCount","phaseWorkspaceBytes","phaseReservationBytes","persistentBytes","stateInputBytes","fluxOutputBytes","knownMaximumLiveOwnedBytes","persistentFullHermitianBytes","gradientMaskBytes","Nx","Ny","Nz","Nj","Nkl","phaseSeconds","reconstructionSeconds","derivativeReconstructionSeconds","productSeconds","projectionSeconds","coefficientAssemblySeconds","derivativeCoefficientAssemblySeconds","coefficientProjectionSeconds","activeForcingSchedule"};
        plhs[0] = mxCreateStructMatrix(1,1,47,names);
        const auto& metrics = value.metrics();
        const auto& configuration = value.descriptor().configuration();
        const auto spectralBytes = value.descriptor().spectralShape().elementCount() * sizeof(WVComplex64);
        mxSetField(plhs[0],0,"engine",mxCreateString(value.engineIdentifier().c_str()));
        mxSetField(plhs[0],0,"loadedLibrary",mxCreateString(value.engineLibraryIdentity().c_str()));
        mxSetField(plhs[0],0,"nonlinearFluxSchedule",mxCreateString(value.nonlinearFluxScheduleIdentifier()));
        mxSetField(plhs[0],0,"phaseImplementation",mxCreateString(value.phaseImplementationIdentifier()));
        mxSetField(plhs[0],0,"coefficientStorageMode",mxCreateString("natural-dimensional"));
        mxSetField(plhs[0],0,"coefficientArithmeticMode",mxCreateString(value.coefficientArithmeticModeIdentifier()));
        mxSetField(plhs[0],0,"inverseNormalizationPlacement",mxCreateString(value.inverseNormalizationPlacementIdentifier()));
        mxSetField(plhs[0],0,"optimizationImplementation",mxCreateString(value.optimizationImplementationIdentifier()));
        mxSetField(plhs[0],0,"coefficientWorkerCount",mxCreateDoubleScalar(static_cast<double>(value.coefficientWorkerCount())));
        mxSetField(plhs[0],0,"planMemoryAccounting",mxCreateString("wrapper-lower-bound; FFTW-owned plan memory is opaque"));
        mxSetField(plhs[0],0,"activeForcingSchedule",mxCreateString(model(prhs[1]).forcingScheduleIdentifier().c_str()));
        const double numbers[] = {static_cast<double>(WVKernelContractVersion),static_cast<double>(metrics.planCount),static_cast<double>(metrics.planBytes),static_cast<double>(metrics.engineBytes),static_cast<double>(metrics.kernelManagementBytes),static_cast<double>(metrics.descriptorBytes),static_cast<double>(metrics.scratchCapacityBytes),static_cast<double>(metrics.scratchHighWaterBytes),static_cast<double>(metrics.halfSpectrumScratchCapacityBytes),static_cast<double>(metrics.realScratchCapacityBytes),static_cast<double>(metrics.executionCount),static_cast<double>(metrics.horizontalExecutionCount),static_cast<double>(metrics.verticalExecutionCount),static_cast<double>(metrics.nonlinearFluxCallCount),static_cast<double>(metrics.nonlinearFluxPhaseEvaluationCount),static_cast<double>(value.phaseReservationBytes()),static_cast<double>(value.phaseReservationBytes()),static_cast<double>(value.persistentBytes()),static_cast<double>(3*spectralBytes),static_cast<double>(3*spectralBytes),static_cast<double>(value.persistentBytes()+3*spectralBytes),0.0,0.0,static_cast<double>(configuration.Nx),static_cast<double>(configuration.Ny),static_cast<double>(configuration.Nz),static_cast<double>(configuration.Nj),static_cast<double>(value.descriptor().Nkl())};
        for (std::size_t i = 0; i < 28; ++i) mxSetField(plhs[0],0,names[i+10],mxCreateDoubleScalar(numbers[i]));
        const double stageSeconds[] = {metrics.phaseSeconds,metrics.reconstructionSeconds,metrics.derivativeReconstructionSeconds,metrics.productSeconds,metrics.projectionSeconds,metrics.coefficientAssemblySeconds,metrics.derivativeCoefficientAssemblySeconds,metrics.coefficientProjectionSeconds};
        for (std::size_t i = 0; i < 8; ++i) mxSetField(plhs[0],0,names[i+38],mxCreateDoubleScalar(stageSeconds[i]));
        return;
    }
    fail("WaveVortexModel:CompiledKernelCommand","Unknown compiled-kernel command.");
}
