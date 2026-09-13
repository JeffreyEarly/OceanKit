#include "WVMatlabTransformHost.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#include "WaveVortexKernel/WVOwnedStratifiedModalSource.hpp"
#include "WaveVortexRuntime/WVNativeVariablePolicy.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexRuntime/WVPortableVariablePlan.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVScopedStateEvaluation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <variant>
#include <type_traits>

namespace {
using namespace wavevortex;
using namespace wavevortex::runtime;
struct BridgeError : std::runtime_error {
    std::string identifier;
    BridgeError(const char* id, const std::string& message)
        : std::runtime_error(message), identifier(id) {}
};
[[noreturn]] void invalid(const std::string& message) {
    throw BridgeError("WaveVortexModel:CompiledTransformInput",message);
}
void require(const WVKernelStatus& status) {
    if (!status) throw BridgeError("WaveVortexModel:CompiledTransformExecution",status.message);
}
struct ArrayDeleter { void operator()(mxArray* a) const { if(a) mxDestroyArray(a); } };
using Array = std::unique_ptr<mxArray,ArrayDeleter>;
std::string text(const mxArray* a) {
    if (!a || !mxIsChar(a)) invalid("Expected a character vector.");
    char* value=mxArrayToString(a);
    if (!value) throw std::bad_alloc();
    std::string result;
    try { result=value; } catch (...) { mxFree(value); throw; }
    mxFree(value); return result;
}
const mxArray* field(const mxArray* a,const char* name) {
    if (!a || !mxIsStruct(a) || mxGetNumberOfElements(a)!=1)
        invalid("Modal configuration must be a scalar struct.");
    const auto* value=mxGetField(a,0,name);
    if (!value) invalid(std::string("Missing modal field ")+name+'.');
    return value;
}
double scalar(const mxArray* a) {
    if (!a || (!mxIsDouble(a) && !mxIsLogical(a)) || mxIsSparse(a) ||
        mxIsComplex(a) || mxGetNumberOfElements(a)!=1 || !std::isfinite(mxGetScalar(a)))
        invalid("Expected a finite real double or logical scalar.");
    return mxGetScalar(a);
}
std::size_t extent(const mxArray* a) {
    const double n=scalar(a);
    if(n<1 || n!=std::floor(n) || n>static_cast<double>(std::numeric_limits<std::size_t>::max()/sizeof(double)))
        invalid("Modal extent must be a positive addressable integer.");
    return static_cast<std::size_t>(n);
}
std::vector<double> array(const mxArray* a,std::initializer_list<std::size_t> shape) {
    if (!a || !mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a))
        invalid("Modal arrays must be full real doubles.");
    const auto* dims=mxGetDimensions(a); const auto rank=mxGetNumberOfDimensions(a);
    std::size_t i=0;
    for (auto n:shape) {
        if ((i<rank ? dims[i] : 1)!=n) invalid("Modal array dimensions do not match the declared geometry.");
        ++i;
    }
    for (;i<rank;++i) if(dims[i]!=1) invalid("Unexpected modal array dimension.");
    const auto* values=mxGetDoubles(a);
    return {values,values+mxGetNumberOfElements(a)};
}
WVStratifiedModalArrays modalArrays(const mxArray* a) {
    if(text(field(a,"schemaVersion"))!="wv-matlab-stratified-modal-source-v1")
        invalid("Unsupported MATLAB modal-source schema.");
    WVStratifiedModalArrays data; auto& g=data.geometry;
    g.transformClass=text(field(a,"transformClass"));
    if(g.transformClass!="WVTransformHydrostatic" && g.transformClass!="WVTransformBoussinesq" &&
       g.transformClass!="WVTransformStratifiedQG" && g.transformClass!="WVTransformBarotropicQG" && g.transformClass!="WVTransformConstantStratification")
        invalid("This modal bridge requires a Hydrostatic, Boussinesq or Stratified QG transform.");
    g.modelVersion=text(field(a,"modelVersion"));
    g.Nx=extent(field(a,"Nx")); g.Ny=extent(field(a,"Ny"));
    g.Nz=extent(field(a,"Nz")); g.Nj=extent(field(a,"Nj")); g.Nkl=extent(field(a,"Nkl"));
    g.Lx=scalar(field(a,"Lx")); g.Ly=scalar(field(a,"Ly")); g.Lz=scalar(field(a,"Lz"));
    g.g=scalar(field(a,"g"));
    if(g.transformClass!="WVTransformBarotropicQG") g.rho0=scalar(field(a,"rho0"));
    g.latitude=scalar(field(a,"latitude")); g.rotationRate=scalar(field(a,"rotationRate"));
    g.planetaryRadius=scalar(field(a,"planetaryRadius"));
    const auto* antialias=field(a,"shouldAntialias");
    if(!mxIsLogicalScalar(antialias)) invalid("shouldAntialias must be a scalar logical.");
    g.shouldAntialias=mxIsLogicalScalarTrue(antialias);
    if(g.transformClass=="WVTransformConstantStratification") return data;
    if(g.transformClass=="WVTransformBarotropicQG") {
        if(g.Nz!=1 || g.Nj!=1) invalid("Barotropic geometry requires Nz=Nj=1.");
        const auto j=scalar(field(a,"j")); if(j!=0 && j!=1) invalid("Barotropic mode index must be zero or one.");
        g.j={j}; return data;
    }
    g.x=array(field(a,"x"),{g.Nx,1}); g.y=array(field(a,"y"),{g.Ny,1});
    g.z=array(field(a,"z"),{g.Nz,1}); g.j=array(field(a,"j"),{g.Nj,1});
    g.k=array(field(a,"k"),{g.Nkl,1}); g.l=array(field(a,"l"),{g.Nkl,1});
    g.N2=array(field(a,"N2"),{g.Nz,1}); g.rho_nm0=array(field(a,"rho_nm0"),{g.Nz,1});
    g.dLnN2=array(field(a,"dLnN2"),{g.Nz,1}); g.z_int=array(field(a,"z_int"),{g.Nz,1});
    g.P0=array(field(a,"P0"),{g.Nj,1}); g.Q0=array(field(a,"Q0"),{g.Nj,1});
    g.h_0=array(field(a,"h_0"),{g.Nj,1});
    const auto k=array(field(a,"kMode"),{g.Nkl,1}), l=array(field(a,"lMode"),{g.Nkl,1});
    for(std::size_t i=0;i<g.Nkl;++i) {
        if(!std::isfinite(k[i]) || !std::isfinite(l[i]) || k[i]!=std::floor(k[i]) || l[i]!=std::floor(l[i]) ||
           std::abs(k[i])>static_cast<double>(g.Nx) || std::abs(l[i])>static_cast<double>(g.Ny))
            invalid("Horizontal mode keys must be exact integers within the grid.");
        g.modes.push_back({static_cast<std::int64_t>(k[i]),static_cast<std::int64_t>(l[i])});
    }
    data.PF0inv=array(field(a,"PF0inv"),{g.Nz,g.Nj}); data.QG0inv=array(field(a,"QG0inv"),{g.Nz,g.Nj});
    data.PF0=array(field(a,"PF0"),{g.Nj,g.Nz}); data.QG0=array(field(a,"QG0"),{g.Nj,g.Nz});
    if(g.transformClass=="WVTransformBoussinesq") {
        const auto* groups=field(a,"K2unique");
        const auto count=mxGetNumberOfElements(groups);
        if(count==0) invalid("Boussinesq requires wave groups.");
        g.K2unique=array(groups,{count,1});
        g.h_pm=array(field(a,"h_pm"),{g.Nj,g.Nkl});
        g.Ppm=array(field(a,"Ppm"),{g.Nj,count}); g.Qpm=array(field(a,"Qpm"),{g.Nj,count});
        const auto membership=array(field(a,"iK2unique"),{g.Nkl,1});
        for(double group:membership) {
            if(!std::isfinite(group) || group!=std::floor(group) || group<1 || group>count)
                invalid("Boussinesq group indices must be exact one-based integers.");
            g.waveGroup.push_back(static_cast<std::size_t>(group)-1);
        }
        data.PFpmInv=array(field(a,"PFpmInv"),{g.Nz,g.Nj,count});
        data.QGpmInv=array(field(a,"QGpmInv"),{g.Nz,g.Nj,count});
        data.PFpm=array(field(a,"PFpm"),{g.Nj,g.Nz,count});
        data.QGpm=array(field(a,"QGpm"),{g.Nj,g.Nz,count});
        data.QGwg=array(field(a,"QGwg"),{g.Nj,g.Nj,count});
    }
    return data;
}
std::vector<std::string> names(const mxArray* a) {
    if(!a || !mxIsCell(a)) invalid("Variable names must be a cell array of character vectors.");
    std::vector<std::string> result;
    for(std::size_t i=0;i<mxGetNumberOfElements(a);++i) result.push_back(text(mxGetCell(a,i)));
    return result;
}
struct Host {
    std::shared_ptr<const WVOwnedStratifiedModalSource> source;
    using Engines=std::variant<std::unique_ptr<WVHydrostaticForcingEngine>,
        std::unique_ptr<WVBoussinesqForcingEngine>,std::unique_ptr<WVStratifiedQGForcingEngine>,std::unique_ptr<WVBarotropicQGForcingEngine>,std::unique_ptr<WVConstantStratificationForcingEngine>>;
    Engines engine;
    WVStratifiedModalGeometry simpleGeometry;
    const WVStratifiedModalGeometry& geometry() const { return source?source->geometry():simpleGeometry; }
    bool isBarotropic() const { return geometry().transformClass=="WVTransformBarotropicQG"; }
    bool isConstant() const { return geometry().transformClass=="WVTransformConstantStratification"; }
    bool isNonhydrostatic() const { return geometry().transformClass=="WVTransformBoussinesq" || (isConstant() && !constantConfiguration.isHydrostatic); }
    WVTransformConstantStratificationConfiguration constantConfiguration;
    bool isQG() const { return isBarotropic() || geometry().transformClass=="WVTransformStratifiedQG"; }
    template<class F> decltype(auto) visit(F&& f) { return std::visit([&](auto& e)->decltype(auto){return f(*e);},engine); }
    template<class F> decltype(auto) visit(F&& f) const { return std::visit([&](const auto& e)->decltype(auto){return f(*e);},engine); }
    std::unique_ptr<WVFieldEvaluationService> fields;
    WVNativeVariablePolicy policy;
    // One explicitly prepared request. Replacing it bounds retained plan memory.
    std::vector<std::string> requested;
    WVFieldEvaluationPlan plan;
    WVIntegrationStateLayout stateLayout;
    bool planPrepared=false;
    bool planActualDensity=true;
    std::size_t evaluations=0, inputBytes=0, outputBytes=0, planPreparations=0, primitiveExecutions=0;
    std::size_t stateCopyBytes=0, scopedEvaluations=0, coefficientOnlyEvaluations=0;
    std::array<std::vector<WVComplex64>,3> stateStorage;
    std::array<WVCoefficientFamilyConstView,3> stateViews;
    WVIntegrationState ownedState;
    std::uint64_t evaluationToken=0;
    WVNoMotionRecoveryReport densityReport;
    bool densityReportAvailable=false;
    // Destroy the session before its state storage and borrowed field service.
    std::unique_ptr<WVFieldEvaluationSession> session;

    explicit Host(WVStratifiedModalArrays data,const mxArray* configuration,std::size_t fftThreads) {
        if(data.geometry.transformClass=="WVTransformBarotropicQG" || data.geometry.transformClass=="WVTransformConstantStratification") simpleGeometry=std::move(data.geometry);
        else require(WVOwnedStratifiedModalSource::create(std::move(data),source));
        const auto kind=isConstant()?WVPersistedTransformKind::constantStratification:isBarotropic()?WVPersistedTransformKind::barotropicQG:isQG()?WVPersistedTransformKind::stratifiedQG:
            geometry().transformClass=="WVTransformBoussinesq"?WVPersistedTransformKind::boussinesq:WVPersistedTransformKind::hydrostatic;
        policy=selectNativeVariablePolicy(true,kind,
            "native-fftw",fftThreads,false,nativeHostTopology());
        if(!isBarotropic() && !isConstant() && (policy.matrixBackend!=WVNativeMatrixBackend::accelerate || !policy.compact))
            invalid("The MATLAB variable backend requires the qualified native Accelerate policy.");
        WVVariableKernelServices services; services.execution=policy.execution;
        services.matrixBackendFactory=WVCreateAccelerateMatrixBackend;
        std::shared_ptr<const WVExtensionCatalog> catalog;
        require(makeBuiltInExtensionCatalog(catalog));
        const auto* registration=catalog->forcings().registration("WVNonlinearAdvection",1);
        if(!registration) invalid("Built-in nonlinear advection is not registered.");
        WVFrozenForcingSchedule schedule; WVFrozenForcingEntry entry;
        entry.typeIdentifier=registration->matlabClassName; entry.contractVersion=1;
        entry.name=registration->defaultName; entry.stage=registration->stage; entry.priority=registration->priority;
        entry.ordinal=1;
        entry.configuration={"wave-vortex-forcing-configuration-v1",1,{}}; schedule.entries.push_back(entry);
        std::unique_ptr<WVFFTEngine> fft; require(WVFFTWEngine::create(policy.effectiveFFTThreads,fft));
        if(kind==WVPersistedTransformKind::hydrostatic) {
            std::unique_ptr<WVHydrostaticForcingEngine> e;
            require(WVHydrostaticForcingEngine::create(source,schedule,catalog,std::move(fft),e,services)); engine=std::move(e);
        } else if(kind==WVPersistedTransformKind::boussinesq) {
            std::unique_ptr<WVBoussinesqForcingEngine> e;
            require(WVBoussinesqForcingEngine::create(source,schedule,catalog,std::move(fft),e,services)); engine=std::move(e);
        } else if(kind==WVPersistedTransformKind::stratifiedQG) {
            std::unique_ptr<WVStratifiedQGForcingEngine> e;
            require(WVStratifiedQGForcingEngine::create(source,schedule,catalog,std::move(fft),e,services));
            require(e->kernel().prepareMatlabPrimitives()); engine=std::move(e);
        }
        if(isBarotropic()) {
            const auto& g=geometry(); WVTransformBarotropicQGConfiguration config;
            config.Nx=g.Nx; config.Ny=g.Ny; config.Lx=g.Lx; config.Ly=g.Ly; config.h=g.Lz; config.j=static_cast<std::uint32_t>(g.j[0]);
            config.g=g.g; config.latitude=g.latitude; config.rotationRate=g.rotationRate; config.planetaryRadius=g.planetaryRadius; config.shouldAntialias=g.shouldAntialias;
            std::unique_ptr<WVBarotropicQGForcingEngine> e;
            require(WVBarotropicQGForcingEngine::create(config,schedule,catalog,std::move(fft),e));
            if(e->kernel().descriptor().Nkl()!=g.Nkl) invalid("Barotropic retained shape differs from MATLAB.");
            engine=std::move(e);
        }
        if(isConstant()) {
            const auto& g=geometry(); auto& c=constantConfiguration;
            c.Nx=g.Nx; c.Ny=g.Ny; c.Nz=g.Nz; c.Nj=g.Nj; c.Lx=g.Lx; c.Ly=g.Ly; c.Lz=g.Lz;
            c.g=g.g; c.rho0=g.rho0; c.latitude=g.latitude; c.rotationRate=g.rotationRate; c.planetaryRadius=g.planetaryRadius;
            c.shouldAntialias=g.shouldAntialias; c.N0=scalar(field(configuration,"N0"));
            const auto* hydrostatic=field(configuration,"isHydrostatic");
            if(!mxIsLogicalScalar(hydrostatic)) invalid("isHydrostatic must be a scalar logical.");
            c.isHydrostatic=mxIsLogicalScalarTrue(hydrostatic);
            std::unique_ptr<WVConstantStratificationForcingEngine> e;
            require(WVConstantStratificationForcingEngine::create(c,schedule,catalog,std::move(fft),e));
            if(e->kernel().descriptor().spectralShape().columns!=g.Nkl) invalid("Constant retained shape differs from MATLAB.");
            require(e->kernel().prepareMatlabPrimitives()); engine=std::move(e);
        }
        visit([&](auto& e){require(WVFieldEvaluationService::createBorrowing(e,fields));});
        require(fields->createStateLayout({},stateLayout));
        for(std::size_t i=0;i<stateLayout.coefficientFamilyCount();++i)
            stateStorage[i].resize(stateLayout.coefficientFamilies()[i].elementCount);
    }
    void prepare(std::vector<std::string> next,bool actualDensity=true) {
        if(planPrepared && next==requested && actualDensity==planActualDensity) return;
        std::vector<WVFieldRequest> requests;
        for(std::size_t i=0;i<next.size();++i) {
            requests.push_back({std::to_string(i),next[i],{}});
        }
        WVFieldEvaluationPlan candidate;
        WVDensityDiagnosticContract density;
        density.reference=actualDensity?WVNoMotionReference::actual:WVNoMotionReference::initial;
        require(session?fields->createPlanForActiveEvaluation(requests,candidate,density):fields->createPlan(requests,candidate,density));
        plan=std::move(candidate); requested=std::move(next); planPrepared=true; planActualDensity=actualDensity; ++planPreparations;
    }
    void endEvaluation() noexcept {
        refreshDensityReport();
        session.reset(); evaluationToken=0; ownedState={};
        stateViews={};
    }
    void refreshDensityReport() noexcept {
        WVNoMotionRecoveryReport report;
        if(fields->activeDensityRecoveryReport(report)) {
            densityReport=report; densityReportAvailable=true;
        }
    }
    std::size_t stateStorageBytes() const noexcept {
        std::size_t bytes=0;
        for(const auto& values:stateStorage) bytes+=values.capacity()*sizeof(WVComplex64);
        return bytes;
    }
};
std::map<std::uint64_t,std::unique_ptr<Host>> hosts;
std::uint64_t nextHandle=std::uint64_t{1}<<63;
std::uint64_t nextEvaluationToken=1;
Host& host(const mxArray* a) {
    if(!a || !mxIsUint64(a) || mxGetNumberOfElements(a)!=1) invalid("Expected a uint64 transform handle.");
    const auto i=hosts.find(*mxGetUint64s(a));
    if(i==hosts.end()) invalid("The transform handle is foreign, deleted or invalid.");
    return *i->second;
}
mxArray* metadata(const Host& h) {
    const char* keys[]={"transformClass","scope","matrixBackend","policy","effectiveFFTThreads","horizontalWorkers","pointwiseWorkers",
        "sourceBytes","kernelBytes","engineBytes","fieldServiceBytes","preparedPlanBytes","planPreparations","evaluations",
        "stateInputBytes","stateInputCopyBytes","outputBytes","stateValidations","phasePreparations",
        "producerExecutions","cacheHits","duplicateExecutions","liveEvaluationBytes","peakEvaluationBytes","tiledNonlinearExecutions","primitiveExecutions","nonlinearProducerExecutions","physicalReconstructionExecutions",
        "stateStorageBytes","scopedEvaluations","coefficientOnlyEvaluations","evaluationActive","evaluationToken",
        "densityRecoveries","densityProfileConstructions","densityInversePasses","densityAPEPasses","plannedEvaluationBytes","densityRecovery"};
    Array result(mxCreateStructMatrix(1,1,sizeof(keys)/sizeof(keys[0]),keys));
    auto put=[&](const char* key,double v){mxSetField(result.get(),0,key,mxCreateDoubleScalar(v));};
    mxSetField(result.get(),0,"transformClass",mxCreateString(h.geometry().transformClass.c_str()));
    mxSetField(result.get(),0,"scope",mxCreateString("call-scoped stratified fields, nonlinear flux and raw primitives"));
    mxSetField(result.get(),0,"matrixBackend",mxCreateString(nativeMatrixBackendIdentifier(h.policy.matrixBackend)));
    mxSetField(result.get(),0,"policy",mxCreateString(std::string(h.policy.selection).c_str()));
    put("effectiveFFTThreads",h.policy.effectiveFFTThreads);
    if(h.isConstant()) {
        const WVConstantKernelExecutionOptions options;
        put("horizontalWorkers",options.horizontalOuterWorkers); put("pointwiseWorkers",options.pointwiseWorkers);
    } else { put("horizontalWorkers",h.policy.execution.horizontalWorkers); put("pointwiseWorkers",h.policy.execution.pointwiseWorkers); }
    put("sourceBytes",h.source?h.source->persistentBytes():0);
    h.visit([&](const auto& e){put("kernelBytes",e.kernel().persistentBytes()); put("engineBytes",e.persistentBytes());});
    put("fieldServiceBytes",h.fields->persistentBytes()); put("preparedPlanBytes",h.plan.persistentBytes());
    put("planPreparations",h.planPreparations); put("evaluations",h.evaluations);
    put("stateInputBytes",h.inputBytes); put("stateInputCopyBytes",h.stateCopyBytes); put("outputBytes",h.outputBytes);
    put("stateStorageBytes",h.stateStorageBytes()); put("scopedEvaluations",h.scopedEvaluations);
    put("coefficientOnlyEvaluations",h.coefficientOnlyEvaluations);
    put("evaluationActive",h.session?1:0);
    auto* token=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
    *mxGetUint64s(token)=h.evaluationToken; mxSetField(result.get(),0,"evaluationToken",token);
    h.visit([&](const auto& engine) {
        const auto& k=engine.kernel().metrics(); const auto& e=engine.variableEvaluationMetrics();
        const auto& f=h.fields->metrics().variableEvaluation;
        const auto active=h.fields->activeVariableEvaluationMetrics();
        if constexpr(std::is_same_v<std::decay_t<decltype(engine)>,WVStratifiedQGForcingEngine> ||
                     std::is_same_v<std::decay_t<decltype(engine)>,WVBarotropicQGForcingEngine>) {
            put("nonlinearProducerExecutions",std::numeric_limits<double>::quiet_NaN());
            put("physicalReconstructionExecutions",std::numeric_limits<double>::quiet_NaN());
        } else {
            put("nonlinearProducerExecutions",engine.metrics().nonlinearProducerCount);
            std::size_t reconstructions=0;
            for(const auto& field:k.reconstructionCount)
                for(const auto count:field[0]) reconstructions+=count;
            put("physicalReconstructionExecutions",reconstructions);
        }
        put("stateValidations",k.stateValidationCount);
        if constexpr((std::is_same_v<std::decay_t<decltype(engine)>,WVStratifiedQGForcingEngine> || std::is_same_v<std::decay_t<decltype(engine)>,WVBarotropicQGForcingEngine>)) {
            put("phasePreparations",0); put("tiledNonlinearExecutions",0);
        } else {
            put("phasePreparations",k.phasePreparationCount);
            if constexpr(std::is_same_v<std::decay_t<decltype(engine)>,WVConstantStratificationForcingEngine>) put("tiledNonlinearExecutions",0);
            else put("tiledNonlinearExecutions",k.tiledNonlinearCount);
        }
        put("producerExecutions",e.producerExecutions+f.producerExecutions+active.producerExecutions);
        put("cacheHits",e.cacheHits+f.cacheHits+active.cacheHits); put("duplicateExecutions",e.duplicateExecutions+f.duplicateExecutions+active.duplicateExecutions);
        put("liveEvaluationBytes",e.liveBytes+f.liveBytes+active.liveBytes);
        put("peakEvaluationBytes",std::max({e.highWaterBytes,f.highWaterBytes,active.highWaterBytes}));
    });
    put("primitiveExecutions",h.primitiveExecutions);
    const auto& f=h.fields->metrics();
    put("densityRecoveries",f.densityRecoveryCount); put("densityProfileConstructions",f.densityProfileConstructionCount);
    put("densityInversePasses",f.densityInversePassCount); put("densityAPEPasses",f.densityAPEPassCount);
    put("plannedEvaluationBytes",f.eventFieldArenaPlannedBytes);
    const char* reportNames[]={"solver","exitflag","iterations","funcCount","maximumResidual","algorithm","message","qualified"};
    auto* report=mxCreateStructMatrix(h.densityReportAvailable?1:0,1,8,reportNames);
    mxSetField(result.get(),0,"densityRecovery",report);
    if(h.densityReportAvailable) {
        const auto& r=h.densityReport;
        mxSetField(report,0,"solver",mxCreateString("dampedLeastSquares"));
        mxSetField(report,0,"exitflag",mxCreateDoubleScalar(r.exitFlag));
        mxSetField(report,0,"iterations",mxCreateDoubleScalar(r.iterations));
        mxSetField(report,0,"funcCount",mxCreateDoubleScalar(r.evaluations));
        mxSetField(report,0,"maximumResidual",mxCreateDoubleScalar(r.maximumResidual));
        mxSetField(report,0,"algorithm",mxCreateString(r.algorithm));
        mxSetField(report,0,"message",mxCreateString(r.reason));
        mxSetField(report,0,"qualified",mxCreateLogicalScalar(r.qualified));
    }
    return result.release();
}
mxArray* supportedVariables(const Host& h) {
    std::vector<std::string> names;
    for(const auto& row:WVPortableVariableContracts)
        if(row.configuration==h.fields->portableVariableConfiguration() &&
           row.runtime==WVPortableDiagnosticRuntime::implemented &&
           std::string_view(row.authority)!="forcing-instance-template")
            names.emplace_back(row.metadata.name);
    Array result(mxCreateCellMatrix(names.size(),1));
    for(std::size_t i=0;i<names.size();++i) mxSetCell(result.get(),i,mxCreateString(names[i].c_str()));
    return result.release();
}
WVComplexConstView coefficients(const mxArray* a,WVShape2D shape) {
    if(!a || !mxIsDouble(a) || mxIsSparse(a) || !mxIsComplex(a) || mxGetM(a)!=shape.rows || mxGetN(a)!=shape.columns)
        invalid("Coefficients must be full complex doubles in the model's [Nj,Nkl] shape.");
    static_assert(sizeof(mxComplexDouble)==sizeof(WVComplex64),"Complex coefficient layout mismatch");
    return {reinterpret_cast<const WVComplex64*>(mxGetComplexDoubles(a)),shape};
}
void requireEvaluation(const Host& h,const mxArray* token) {
    if(!token || !mxIsUint64(token) || mxGetNumberOfElements(token)!=1 || !h.session ||
        *mxGetUint64s(token)!=h.evaluationToken)
        throw BridgeError("WaveVortexModel:CompiledTransformEvaluation","The evaluation token is foreign, expired or invalid.");
}
mxArray* beginEvaluation(Host& h,const mxArray* const inputs[]) {
    if(h.session) invalid("End the active evaluation before supplying a new state.");
    if(nextEvaluationToken==std::numeric_limits<std::uint64_t>::max()) invalid("Evaluation token space exhausted.");
    const auto& g=h.geometry(); const WVShape2D shape{g.Nj,g.Nkl};
    WVState state{scalar(inputs[3]),scalar(inputs[4]),{{},{},coefficients(inputs[2],shape)}};
    if(h.isQG()) {
        if(mxGetNumberOfElements(inputs[0]) || mxGetNumberOfElements(inputs[1])) invalid("QG has only A0 coefficients.");
    } else { state.coefficients.Ap=coefficients(inputs[0],shape); state.coefficients.Am=coefficients(inputs[1],shape); }
    if(!h.isQG()) require(h.fields->prepareBuiltinNonlinearCoefficientEvaluation());
    const auto& families=h.stateLayout.coefficientFamilies();
    const WVComplexConstView borrowed[]={state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0};
    for(std::size_t i=0;i<families.size();++i) {
        const auto source=borrowed[h.isQG()?2:i];
        std::copy_n(source.data,families[i].elementCount,h.stateStorage[i].data());
        h.stateViews[i]={&families[i],h.stateStorage[i].data()};
    }
    h.inputBytes+=families.size()*shape.elementCount()*sizeof(WVComplex64);
    h.stateCopyBytes+=families.size()*shape.elementCount()*sizeof(WVComplex64);
    h.ownedState={}; h.ownedState.waveVortex=state;
    h.densityReportAvailable=false;
    if(h.isQG()) h.ownedState.waveVortex.coefficients.A0={h.stateStorage[0].data(),shape};
    else h.ownedState.waveVortex.coefficients={{h.stateStorage[0].data(),shape},{h.stateStorage[1].data(),shape},{h.stateStorage[2].data(),shape}};
    h.ownedState.coefficientFamilies=h.stateViews.data(); h.ownedState.coefficientFamilyCount=families.size();
    auto candidate=std::make_unique<WVFieldEvaluationSession>();
    Array result(mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL));
    require(h.fields->beginEvaluationSession(h.ownedState,*candidate));
    h.session=std::move(candidate); h.evaluationToken=nextEvaluationToken++; ++h.scopedEvaluations;
    *mxGetUint64s(result.get())=h.evaluationToken; return result.release();
}
mxArray* queryEvaluation(Host& h,const mxArray* token,const mxArray* requests,const mxArray* actualDensity) {
    requireEvaluation(h,token);
    if(actualDensity && !mxIsLogicalScalar(actualDensity)) invalid("The actual-density selection must be a scalar logical.");
    h.prepare(names(requests),!actualDensity || mxIsLogicalScalarTrue(actualDensity));
    const char* keys[]={"values","metrics"}; Array result(mxCreateStructMatrix(1,1,2,keys));
    auto* values=mxCreateCellMatrix(h.plan.outputCount(),1); mxSetField(result.get(),0,"values",values);
    std::vector<WVFieldOutputView> outputs; outputs.reserve(h.plan.outputCount()); std::size_t bytes=0;
    for(std::size_t i=0;i<h.plan.outputCount();++i) {
        const auto& spec=h.plan.outputs()[i]; std::vector<mwSize> dims(spec.dimensions.begin(),spec.dimensions.end());
        while(dims.size()<2) dims.push_back(1);
        // MATLAB stores the one-dimensional barotropic coefficient family as a row.
        if(h.isBarotropic() && h.requested[i]=="A0t") dims={1,h.geometry().Nkl};
        auto* a=mxCreateNumericArray(dims.size(),dims.data(),mxDOUBLE_CLASS,spec.isComplex?mxCOMPLEX:mxREAL);
        mxSetCell(values,i,a); WVFieldOutputView output; output.elementCount=spec.elementCount;
        if(spec.isComplex) output.complexData=reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(a)); else output.data=mxGetDoubles(a);
        outputs.push_back(output); bytes+=spec.elementCount*(spec.isComplex?sizeof(WVComplex64):sizeof(double));
    }
    const auto status=h.fields->evaluate(h.plan,h.ownedState,outputs.data(),outputs.size());
    h.refreshDensityReport(); require(status);
    ++h.evaluations; h.outputBytes+=bytes;
    mxSetField(result.get(),0,"metrics",metadata(h)); return result.release();
}
mxArray* nonlinearCoefficients(Host& h,const mxArray* token) {
    requireEvaluation(h,token);
    if(h.isQG()) invalid("The nonlinear coefficient boundary requires a wave transform.");
    const auto& g=h.geometry(); const WVShape2D shape{g.Nj,g.Nkl};
    const char* keys[]={"values","metrics"};
    Array result(mxCreateStructMatrix(1,1,2,keys));
    auto* values=mxCreateCellMatrix(3,1); mxSetField(result.get(),0,"values",values);
    WVFlux flux; WVComplexView* channels[]={&flux.Fp,&flux.Fm,&flux.F0};
    for(std::size_t channel=0;channel<3;++channel) {
        auto* value=mxCreateDoubleMatrix(shape.rows,shape.columns,mxCOMPLEX);
        mxSetCell(values,channel,value);
        *channels[channel]={reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(value)),shape};
    }
    require(h.fields->evaluateBuiltinNonlinearCoefficients(h.ownedState,flux));
    ++h.evaluations; h.outputBytes+=3*shape.elementCount()*sizeof(WVComplex64);
    mxSetField(result.get(),0,"metrics",metadata(h));
    return result.release();
}
mxArray* coefficientOnlyRightHandSide(Host& h,const mxArray* const inputs[]) {
    if(h.session) invalid("A coefficient-only RHS requires no active field evaluation.");
    const auto& g=h.geometry(); const WVShape2D shape{g.Nj,g.Nkl};
    WVState state{scalar(inputs[3]),scalar(inputs[4]),{{},{},coefficients(inputs[2],shape)}};
    if(h.isQG()) {
        if(mxGetNumberOfElements(inputs[0]) || mxGetNumberOfElements(inputs[1])) invalid("QG has only A0 coefficients.");
    } else {
        state.coefficients.Ap=coefficients(inputs[0],shape);
        state.coefficients.Am=coefficients(inputs[1],shape);
    }
    const std::size_t families=h.isQG()?1:3;
    const char* keys[]={"values","metrics"};
    Array result(mxCreateStructMatrix(1,1,2,keys));
    auto* values=mxCreateCellMatrix(families,1); mxSetField(result.get(),0,"values",values);
    WVFlux flux; WVComplexView* channels[]={&flux.Fp,&flux.Fm,&flux.F0};
    for(std::size_t family=0;family<families;++family) {
        auto* value=mxCreateDoubleMatrix(shape.rows,shape.columns,mxCOMPLEX);
        mxSetCell(values,family,value);
        *channels[h.isQG()?2:family]={reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(value)),shape};
    }
    // This event has a closed demand set: only coefficient tendencies escape.
    // The same resolved native RHS owns its temporary dependencies and context.
    h.visit([&](auto& engine) {
        using Engine=std::decay_t<decltype(engine)>;
        if constexpr(std::is_same_v<Engine,WVBarotropicQGForcingEngine> ||
                     std::is_same_v<Engine,WVStratifiedQGForcingEngine>)
            require(engine.evaluateRightHandSide(state.coefficients.A0,flux.F0,nullptr));
        else require(engine.nonlinearFlux(state,flux));
    });
    ++h.evaluations; ++h.coefficientOnlyEvaluations;
    const auto bytes=families*shape.elementCount()*sizeof(WVComplex64);
    h.inputBytes+=bytes; h.outputBytes+=bytes;
    mxSetField(result.get(),0,"metrics",metadata(h));
    return result.release();
}
mxArray* evaluate(Host& h,const mxArray* const inputs[]) {
    if(h.session) invalid("Use the active evaluation token to request variables in this scope.");
    const auto& g=h.geometry(); const WVShape2D shape{g.Nj,g.Nkl};
    WVState state{scalar(inputs[3]),scalar(inputs[4]),{{},{},coefficients(inputs[2],shape)}};
    if(h.isQG()) {
        if(mxGetNumberOfElements(inputs[0]) || mxGetNumberOfElements(inputs[1])) invalid("QG has only A0 coefficients.");
    } else { state.coefficients.Ap=coefficients(inputs[0],shape); state.coefficients.Am=coefficients(inputs[1],shape); }
    const auto requested=names(inputs[5]);
    if(!h.planPrepared || requested!=h.requested) invalid("Prepare these variables before evaluating the transform.");
    if(!mxIsLogicalScalar(inputs[6])) invalid("shouldEvaluateFlux must be a scalar logical.");
    const bool fluxRequested=mxIsLogicalScalarTrue(inputs[6]);
    const std::size_t familyCount=h.isQG()?1:3;
    const char* keys[]={"values","flux","metrics"};
    Array result(mxCreateStructMatrix(1,1,3,keys));
    mxSetField(result.get(),0,"values",mxCreateCellMatrix(requested.size(),1));
    mxSetField(result.get(),0,"flux",mxCreateCellMatrix(fluxRequested?familyCount:0,1));
    auto* values=mxGetField(result.get(),0,"values"); auto* fluxValues=mxGetField(result.get(),0,"flux");
    std::vector<WVFieldOutputView> outputs; std::size_t outputBytes=0;
    for(std::size_t i=0;i<h.plan.outputCount();++i) {
        const auto& spec=h.plan.outputs()[i]; std::vector<mwSize> dims(spec.dimensions.begin(),spec.dimensions.end());
        while(dims.size()<2) dims.push_back(1);
        // MATLAB stores the one-dimensional barotropic coefficient family as a row.
        if(h.isBarotropic() && h.requested[i]=="A0t") dims={1,h.geometry().Nkl};
        auto* a=mxCreateNumericArray(dims.size(),dims.data(),mxDOUBLE_CLASS,spec.isComplex?mxCOMPLEX:mxREAL);
        mxSetCell(values,i,a); WVFieldOutputView out; out.elementCount=spec.elementCount;
        if(spec.isComplex) out.complexData=reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(a)); else out.data=mxGetDoubles(a);
        outputs.push_back(out); outputBytes+=spec.elementCount*(spec.isComplex?sizeof(WVComplex64):sizeof(double));
    }
    WVFlux flux;
    if(fluxRequested) {
        WVComplexView* channels[]={&flux.Fp,&flux.Fm,&flux.F0};
        for(std::size_t i=0;i<familyCount;++i) {
            auto* a=mxCreateDoubleMatrix(shape.rows,shape.columns,mxCOMPLEX); mxSetCell(fluxValues,i,a);
            *channels[h.isQG()?2:i]={reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(a)),shape};
        }
        outputBytes+=familyCount*shape.elementCount()*sizeof(WVComplex64);
    }
    // All fields and raw nonlinear tendencies use the same diagnostic event.
    // Projection occurs once after raw spatial contributions are available.
    WVIntegrationState integrationState; integrationState.waveVortex=state;
    const auto& families=h.stateLayout.coefficientFamilies();
    std::vector<WVCoefficientFamilyConstView> views;
    if(h.isQG()) views.push_back({&families[0],state.coefficients.A0.data});
    else { views.push_back({&families[0],state.coefficients.Ap.data}); views.push_back({&families[1],state.coefficients.Am.data}); views.push_back({&families[2],state.coefficients.A0.data}); }
    integrationState.coefficientFamilies=views.data(); integrationState.coefficientFamilyCount=familyCount;
    if(fluxRequested) {
        std::vector<std::string> rawNames=h.isQG()?std::vector<std::string>{"Fqgpv_nonlinear_advection"}:
            h.isNonhydrostatic()?std::vector<std::string>{"Fu_nonlinear_advection","Fv_nonlinear_advection","Fw_nonlinear_advection","Feta_nonlinear_advection"}:
            std::vector<std::string>{"Fu_nonlinear_advection","Fv_nonlinear_advection","Feta_nonlinear_advection"};
        std::vector<WVFieldRequest> requests;
        for(const auto& name:requested) requests.push_back({std::to_string(requests.size()),name,{}});
        for(const auto& name:rawNames) requests.push_back({std::to_string(requests.size()),name,{}});
        WVFieldEvaluationPlan combined; require(h.fields->createPlan(requests,combined));
        const auto count=g.Nx*g.Ny*g.Nz;
        std::vector<double> raw(rawNames.size()*count);
        for(std::size_t i=0;i<rawNames.size();++i) {
            WVFieldOutputView out; out.elementCount=count; out.data=raw.data()+i*count; outputs.push_back(out);
        }
        WVFieldEvaluationSession session; require(h.fields->beginEvaluationSession(integrationState,session));
        require(h.fields->evaluate(combined,integrationState,outputs.data(),outputs.size()));
        h.visit([&](auto& engine) {
            using Engine=std::decay_t<decltype(engine)>;
            if constexpr(std::is_same_v<Engine,WVBarotropicQGForcingEngine>)
                require(engine.kernel().transformQGPVToA0({raw.data(),{g.Nx,g.Ny}},flux.F0));
            else if constexpr(std::is_same_v<Engine,WVStratifiedQGForcingEngine>)
                require(engine.kernel().transformQGPVToA0({raw.data(),{g.Nx,g.Ny,g.Nz}},flux.F0));
            else {
                const WVShape3D volume{g.Nx,g.Ny,g.Nz};
                WVRealVolumeConstView u{raw.data(),volume},v{raw.data()+count,volume},eta{raw.data()+(rawNames.size()-1)*count,volume};
                WVMutableCoefficients out{flux.Fp,flux.Fm,flux.F0};
                if constexpr(std::is_same_v<Engine,WVBoussinesqForcingEngine>)
                    require(engine.kernel().transformUVWEtaToWaveVortex(u,v,{raw.data()+2*count,volume},eta,state.t,state.t0,out));
                else if constexpr(std::is_same_v<Engine,WVConstantStratificationForcingEngine>) {
                    if(h.isNonhydrostatic()) require(engine.kernel().transformUVWEtaToWaveVortex(u,v,{raw.data()+2*count,volume},eta,state.t,state.t0,out));
                    else require(engine.kernel().transformUVEtaToWaveVortex(u,v,eta,state.t,state.t0,out));
                } else require(engine.kernel().transformUVEtaToWaveVortex(u,v,eta,state.t,state.t0,out));
            }
        });
    } else {
        WVFieldEvaluationSession session; require(h.fields->beginEvaluationSession(integrationState,session));
        if(!outputs.empty()) require(h.fields->evaluate(h.plan,integrationState,outputs.data(),outputs.size()));
    }
    ++h.evaluations; h.inputBytes+=familyCount*shape.elementCount()*sizeof(WVComplex64); h.outputBytes+=outputBytes;
    mxSetField(result.get(),0,"metrics",metadata(h));
    return result.release();
}
void dimensions(const mxArray* a,std::initializer_list<std::size_t> shape) {
    const auto* dims=mxGetDimensions(a); const auto rank=mxGetNumberOfDimensions(a);
    std::size_t i=0;
    for(auto n:shape) { if((i<rank?dims[i]:1)!=n) invalid("Primitive input dimensions do not match the operation."); ++i; }
    for(;i<rank;++i) if(dims[i]!=1) invalid("Unexpected primitive input dimension.");
}
WVRealVolumeConstView spatial(const mxArray* a,WVShape3D shape) {
    if(!a || !mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a)) invalid("Spatial inputs must be full real doubles.");
    dimensions(a,{shape.first,shape.second,shape.third}); return {mxGetDoubles(a),shape};
}
WVStratifiedModalOperator verticalOperation(const std::string& name) {
    using Op=WVStratifiedModalOperator;
    const std::pair<const char*,Op> operations[]={
        {"reconstructF",Op::reconstructF},{"projectF",Op::projectF},
        {"reconstructG",Op::reconstructG},{"projectG",Op::projectG},
        {"reconstructFw",Op::reconstructFw},{"projectFw",Op::projectFw},
        {"reconstructGw",Op::reconstructGw},{"projectGw",Op::projectGw},
        {"balancedGToWaveG",Op::balancedGToWaveG},
        {"projectWaveDivergence",Op::projectWaveDivergence},
        {"projectWaveVerticalVelocity",Op::projectWaveVerticalVelocity}};
    for(const auto& item:operations) if(name==item.first) return item.second;
    invalid("Unknown vertical operator.");
}
mxArray* operation(Host& h,const mxArray* nameArray,const mxArray* inputs,const mxArray* options) {
    const auto name=text(nameArray);
    if(!inputs || !mxIsCell(inputs)) invalid("Primitive inputs must be a cell array.");
    if(!options || !mxIsStruct(options) || mxGetNumberOfElements(options)!=1) invalid("Primitive options must be a scalar struct.");
    const auto& g=h.geometry(); const WVShape2D modal{g.Nj,g.Nkl}, grid{g.Nz,g.Nkl};
    const WVShape3D volume{g.Nx,g.Ny,g.Nz};
    const bool projection=name=="toWaveVortex";
    const std::size_t count=projection&&!h.isQG()?3:1;
    const char* keys[]={"values","metrics"}; Array result(mxCreateStructMatrix(1,1,2,keys));
    mxSetField(result.get(),0,"values",mxCreateCellMatrix(count,1)); auto* values=mxGetField(result.get(),0,"values");
    std::size_t bytes=0;
    auto expect=[&](std::size_t n){if(mxGetNumberOfElements(inputs)!=n) invalid("Wrong number of primitive input arrays.");};
    auto input=[&](std::size_t index){return mxGetCell(inputs,index);};
    auto complexOutput=[&](std::size_t index,WVShape2D shape) {
        auto* a=mxCreateDoubleMatrix(shape.rows,shape.columns,mxCOMPLEX); mxSetCell(values,index,a);
        bytes+=shape.elementCount()*sizeof(WVComplex64);
        return WVComplexView{reinterpret_cast<WVComplex64*>(mxGetComplexDoubles(a)),shape};
    };
    auto realOutput=[&]() {
        const mwSize dims[]={g.Nx,g.Ny,g.Nz}; auto* a=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL); mxSetCell(values,0,a);
        bytes+=g.Nx*g.Ny*g.Nz*sizeof(double); return WVRealVolumeView{mxGetDoubles(a),volume};
    };
    if(name=="vertical" || name=="verticalColumn") {
        expect(1); const auto op=verticalOperation(text(field(options,"operator")));
        using Op=WVStratifiedModalOperator;
        const bool reconstruction=op==Op::reconstructF || op==Op::reconstructG || op==Op::reconstructFw || op==Op::reconstructGw;
        const bool modalInput=reconstruction || op==Op::balancedGToWaveG;
        const bool column=name=="verticalColumn";
        const WVShape2D inShape{modalInput?g.Nj:g.Nz,column?1:g.Nkl}, outShape{reconstruction?g.Nz:g.Nj,column?1:g.Nkl};
        const auto in=coefficients(input(0),inShape); auto out=complexOutput(0,outShape);
        if(column) {
            const auto index=extent(field(options,"column"));
            if(index>g.Nkl) invalid("Retained column exceeds the model spectrum.");
            h.visit([&](auto& e){if constexpr(std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>) invalid("Barotropic QG has no vertical matrix primitive."); else require(e.kernel().applyVerticalColumn(op,index-1,in,out));});
        } else h.visit([&](auto& e){if constexpr(std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>) invalid("Barotropic QG has no vertical matrix primitive."); else require(e.kernel().applyVertical(op,in,out));});
    } else if(name=="verticalCalculus") {
        expect(1);
        if(h.isBarotropic()) invalid("Barotropic QG has no vertical calculus.");
        const auto* a=input(0);
        if(!a || !mxIsDouble(a) || mxIsSparse(a) || mxGetNumberOfDimensions(a)!=2 || mxGetM(a)!=g.Nz)
            invalid("Vertical calculus requires a full double [Nz,Ncolumns] matrix.");
        const auto* family=field(options,"inputIsF"); const auto* integral=field(options,"integral");
        if(!mxIsLogicalScalar(family) || !mxIsLogicalScalar(integral)) invalid("Vertical calculus flags must be scalar logicals.");
        const auto order=extent(field(options,"order")); const bool isIntegral=mxIsLogicalScalarTrue(integral);
        if(order>4 || (isIntegral && order!=1)) invalid("Vertical calculus supports derivative orders 1..4 and first antiderivatives.");
        const WVShape2D shape{g.Nz,mxGetN(a)};
        auto* b=mxCreateDoubleMatrix(shape.rows,shape.columns,mxIsComplex(a)?mxCOMPLEX:mxREAL); mxSetCell(values,0,b);
        bytes+=shape.elementCount()*(mxIsComplex(a)?sizeof(WVComplex64):sizeof(double));
        const auto apply=[&](const double* in,double* out) {
            h.visit([&](auto& e) {
                if constexpr(!std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>)
                    require(e.kernel().applyVerticalCalculus({in,shape},mxIsLogicalScalarTrue(family),static_cast<unsigned>(order),isIntegral,{out,shape}));
            });
        };
        if(shape.columns && mxIsComplex(a)) {
            // The scientific operator is real and linear. Preserve MATLAB's
            // accepted complex input without retaining cross-call array views.
            const auto* in=mxGetComplexDoubles(a); auto* out=mxGetComplexDoubles(b);
            std::vector<double> part(shape.elementCount()),transformed(part.size());
            for(std::size_t i=0;i<part.size();++i) part[i]=in[i].real;
            apply(part.data(),transformed.data());
            for(std::size_t i=0;i<part.size();++i) { out[i].real=transformed[i]; part[i]=in[i].imag; }
            apply(part.data(),transformed.data());
            for(std::size_t i=0;i<part.size();++i) out[i].imag=transformed[i];
        } else if(shape.columns) apply(mxGetDoubles(a),mxGetDoubles(b));
    } else if(name=="horizontalForward") {
        expect(1); const auto in=spatial(input(0),volume); auto out=complexOutput(0,grid);
        h.visit([&](auto& e){if constexpr(std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>) require(e.kernel().horizontalForward({in.data,{g.Nx,g.Ny}},out)); else require(e.kernel().horizontalForward(in,out));});
    } else if(name=="horizontalInverse") {
        expect(1); const auto in=coefficients(input(0),grid); auto out=realOutput();
        h.visit([&](auto& e){if constexpr(std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>) { WVRealView target{out.data,{g.Nx,g.Ny}}; require(e.kernel().horizontalInverse(in,target)); } else require(e.kernel().horizontalInverse(in,out));});
    } else if(name=="differentiateHorizontal") {
        expect(1); const auto direction=text(field(options,"direction"));
        if(direction!="x" && direction!="y") invalid("Derivative direction must be x or y.");
        const auto order=extent(field(options,"order"));
        if(order>std::numeric_limits<unsigned>::max()) invalid("Derivative order is not representable.");
        const auto in=spatial(input(0),volume); auto out=realOutput();
        h.visit([&](auto& e){if constexpr(std::is_same_v<std::decay_t<decltype(e)>,WVBarotropicQGForcingEngine>) { WVRealView target{out.data,{g.Nx,g.Ny}}; require(e.kernel().differentiateHorizontal({in.data,{g.Nx,g.Ny}},direction=="x",static_cast<unsigned>(order),target)); } else require(e.kernel().differentiateHorizontal(in,direction=="x",static_cast<unsigned>(order),out));});
    } else if(name=="toWaveVortex") {
        const bool bouss=h.isNonhydrostatic(); expect(bouss?4:3);
        const auto u=spatial(input(0),volume), v=spatial(input(1),volume), eta=spatial(input(bouss?3:2),volume);
        h.visit([&](auto& e) {
            using Engine=std::decay_t<decltype(e)>;
            if constexpr(std::is_same_v<Engine,WVBarotropicQGForcingEngine>) {
                invalid("Barotropic QG supports qgPVToA0 projection.");
            } else if constexpr(std::is_same_v<Engine,WVStratifiedQGForcingEngine>) {
                require(e.kernel().transformUVEtaToA0(u,v,eta,complexOutput(0,modal)));
            } else {
                const double t=scalar(field(options,"t")), t0=scalar(field(options,"t0"));
                WVMutableCoefficients output{complexOutput(0,modal),complexOutput(1,modal),complexOutput(2,modal)};
                if constexpr(std::is_same_v<Engine,WVBoussinesqForcingEngine>)
                    require(e.kernel().transformUVWEtaToWaveVortex(u,v,spatial(input(2),volume),eta,t,t0,output));
                else if constexpr(std::is_same_v<Engine,WVConstantStratificationForcingEngine>) {
                    if(bouss) require(e.kernel().transformUVWEtaToWaveVortex(u,v,spatial(input(2),volume),eta,t,t0,output));
                    else require(e.kernel().transformUVEtaToWaveVortex(u,v,eta,t,t0,output));
                }
                else require(e.kernel().transformUVEtaToWaveVortex(u,v,eta,t,t0,output));
            }
        });
    } else if(name=="qgPVToA0") {
        expect(1); if(!h.isQG()) invalid("qgPVToA0 requires a QG transform.");
        const auto in=spatial(input(0),volume); auto out=complexOutput(0,modal);
        if(h.isBarotropic()) {
            auto& engine=*std::get<std::unique_ptr<WVBarotropicQGForcingEngine>>(h.engine);
            require(engine.kernel().transformQGPVToA0({in.data,{g.Nx,g.Ny}},out));
        } else {
            auto& engine=*std::get<std::unique_ptr<WVStratifiedQGForcingEngine>>(h.engine);
            require(engine.kernel().transformQGPVToA0(in,out));
        }
    } else invalid("Unknown transform primitive.");
    ++h.primitiveExecutions; h.outputBytes+=bytes;
    mxSetField(result.get(),0,"metrics",metadata(h)); return result.release();
}

}

std::size_t WVMatlabTransformCount() noexcept { return hosts.size(); }
void WVCleanupMatlabTransforms() noexcept { hosts.clear(); }

bool WVDispatchMatlabTransform(const std::string& command,int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if(command.rfind("transform",0)!=0) return false;
    // MATLAB error dispatch happens after stack unwinding, never inside a live
    // C++ scope owning scientific resources or partially constructed outputs.
    std::string errorId,errorMessage;
    try {
        if(command=="transformCreate") {
            if(nrhs!=3 || nlhs!=1) invalid("transformCreate requires a modal struct, validated FFT thread count and one output.");
            auto candidate=std::make_unique<Host>(modalArrays(prhs[1]),prhs[1],extent(prhs[2]));
            Array output(mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL));
            if(nextHandle==std::numeric_limits<std::uint64_t>::max()) invalid("Transform handle space exhausted.");
            const auto id=nextHandle++; *mxGetUint64s(output.get())=id;
            hosts.emplace(id,std::move(candidate));
            if(hosts.size()==1) mexLock();
            plhs[0]=output.release(); return true;
        }
        if(nrhs<2) invalid("A transform command requires its handle.");
        auto& h=host(prhs[1]);
        if(command=="transformDelete") {
            if(nrhs!=2 || nlhs!=0) invalid("transformDelete takes one handle and no output.");
            hosts.erase(*mxGetUint64s(prhs[1])); if(hosts.empty()) mexUnlock();
        } else if(command=="transformMetadata") {
            if(nrhs!=2 || nlhs!=1) invalid("transformMetadata takes one handle and returns one struct.");
            plhs[0]=metadata(h);
        } else if(command=="transformVariables") {
            if(nrhs!=2 || nlhs!=1) invalid("transformVariables takes a handle and returns implemented variable names.");
            plhs[0]=supportedVariables(h);
        } else if(command=="transformPrepare") {
            if(nrhs!=3 || nlhs!=0) invalid("transformPrepare takes handle and variable names with no output.");
            h.prepare(names(prhs[2]));
        } else if(command=="transformBeginEvaluation") {
            if(nrhs!=7 || nlhs!=1) invalid("transformBeginEvaluation takes handle, Ap, Am, A0, t and t0.");
            plhs[0]=beginEvaluation(h,prhs+2);
        } else if(command=="transformQueryEvaluation") {
            if((nrhs!=4 && nrhs!=5) || nlhs!=1) invalid("transformQueryEvaluation takes handle, token, variable names and optional density reference.");
            plhs[0]=queryEvaluation(h,prhs[2],prhs[3],nrhs==5?prhs[4]:nullptr);
        } else if(command=="transformNonlinearCoefficients") {
            if(nrhs!=3 || nlhs!=1) invalid("transformNonlinearCoefficients takes handle and token with one output.");
            plhs[0]=nonlinearCoefficients(h,prhs[2]);
        } else if(command=="transformCoefficientOnlyRightHandSide") {
            if(nrhs!=7 || nlhs!=1) invalid("transformCoefficientOnlyRightHandSide takes handle, Ap, Am, A0, t and t0 with one output.");
            plhs[0]=coefficientOnlyRightHandSide(h,prhs+2);
        } else if(command=="transformEndEvaluation") {
            if(nrhs!=3 || nlhs!=0) invalid("transformEndEvaluation takes handle and token with no output.");
            requireEvaluation(h,prhs[2]); h.endEvaluation();
        } else if(command=="transformOperation") {
            if(nrhs!=5 || nlhs!=1) invalid("transformOperation takes handle, operation name, inputs and options.");
            plhs[0]=operation(h,prhs[2],prhs[3],prhs[4]);
        } else if(command=="transformEvaluate") {
            if(nrhs!=9 || nlhs!=1) invalid("transformEvaluate takes handle, Ap, Am, A0, t, t0, names and flux flag.");
            plhs[0]=evaluate(h,prhs+2);
        } else invalid("Unknown transform command.");
        return true;
    } catch(const BridgeError& e) { errorId=e.identifier; errorMessage=e.what(); }
      catch(const std::bad_alloc&) { errorId="WaveVortexModel:CompiledTransformAllocation"; errorMessage="Unable to allocate compiled transform storage."; }
      catch(const std::exception& e) { errorId="WaveVortexModel:CompiledTransformFailure"; errorMessage=e.what(); }
    mexErrMsgIdAndTxt(errorId.c_str(),"%s",errorMessage.c_str());
    return true;
}
