#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WVSpectralValidation.hpp"
#include "WVPreparedModeExecutor.hpp"
#include "WVPreparedFieldCache.hpp"
#include "WVAdvectionConsumer.hpp"
#include "WVVariableComplexBuffer.hpp"
#include "WVStratifiedVerticalCalculus.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <system_error>

namespace wavevortex {
namespace {
using namespace spectral_detail;
constexpr double pi = 3.1415926535897932384626433832795;
WVComplex64 add(WVComplex64 a,WVComplex64 b) { return {a.real+b.real,a.imag+b.imag}; }
WVComplex64 subtract(WVComplex64 a,WVComplex64 b) { return {a.real-b.real,a.imag-b.imag}; }
WVComplex64 scale(WVComplex64 a,double b) { return {a.real*b,a.imag*b}; }
WVComplex64 conjugate(WVComplex64 a) { return {a.real,-a.imag}; }
WVComplex64 multiply(WVComplex64 a,WVComplex64 b) { return {a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real}; }
WVComplexInput input(const WVComplex64* p,std::size_t n) { return {p,nullptr,nullptr,n*sizeof(WVComplex64)}; }
WVComplexOutput output(WVComplex64* p,std::size_t n) { return {p,nullptr,nullptr,n*sizeof(WVComplex64)}; }
void copy(WVComplexInput source,WVComplex64* destination,std::size_t count) {
    for (std::size_t i=0;i<count;++i) destination[i]=read(source,i);
}
void copy(const WVComplex64* source,WVComplexOutput destination,std::size_t count) {
    for (std::size_t i=0;i<count;++i) write(destination,i,source[i]);
}
bool selfConjugate(std::int64_t mode,std::size_t count) {
    return mode==0 || (count%2==0 &&
        (mode==static_cast<std::int64_t>(count/2) || mode==-static_cast<std::int64_t>(count/2)));
}
bool verticalOperation(WVStratifiedModalOperator operation,std::size_t& index,
    bool& inputModal,bool& outputModal) {
    switch(operation) {
        case WVStratifiedModalOperator::reconstructF: index=0; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectF: index=1; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::reconstructG: index=2; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectG: index=3; inputModal=false; outputModal=true; return true;
        default: return false;
    }
}
WVCoefficients view(WVMutableCoefficients a) { return {{a.Ap.data,a.Ap.shape},{a.Am.data,a.Am.shape},{a.A0.data,a.A0.shape}}; }
bool valid(WVHydrostaticFamily f) { return f==WVHydrostaticFamily::F || f==WVHydrostaticFamily::G; }
bool valid(WVHydrostaticComponent c) { return c>=WVHydrostaticComponent::all && c<=WVHydrostaticComponent::meanDensityAnomaly; }
bool surface(WVHydrostaticField f) { return f>=WVHydrostaticField::ssh && f<=WVHydrostaticField::ssv; }
bool selected(const WVHydrostaticModeFactors& f,WVHydrostaticComponent c,bool wave) {
    if (c==WVHydrostaticComponent::all) return wave ? f.wave || f.inertial : f.geostrophic || f.meanDensityAnomaly;
    if (wave) return c==WVHydrostaticComponent::wave ? f.wave : c==WVHydrostaticComponent::inertial && f.inertial;
    return c==WVHydrostaticComponent::geostrophic ? f.geostrophic : c==WVHydrostaticComponent::meanDensityAnomaly && f.meanDensityAnomaly;
}
WVKernelStatus unsupported() { return {WVKernelStatusCode::unsupportedOperation,"Unsupported hydrostatic field, derivative, family or component."}; }
WVKernelStatus reentrant() { return {WVKernelStatusCode::reentrantExecution,"Hydrostatic workspace is already active."}; }
}

WVTransformHydrostaticKernel::~WVTransformHydrostaticKernel() = default;

WVKernelStatus WVTransformHydrostaticKernel::create(std::shared_ptr<const WVStratifiedModalSource> source,
    std::unique_ptr<WVFFTEngine> engine,std::unique_ptr<WVTransformHydrostaticKernel>& result,MatrixBackendFactory factory,WVVariableExecutionOptions options) {
    try {
        if (!source || !engine || !factory) return {WVKernelStatusCode::invalidConfiguration,"Scientific source, FFT engine and matrix backend factory are required."};
        if (options.spectralSchedule!=WVVariableSpectralSchedule::establishedInterleaved && options.spectralSchedule!=WVVariableSpectralSchedule::compactSplitFusedViews)
            return {WVKernelStatusCode::invalidConfiguration,"Unknown variable spectral schedule."};
        if (!options.pointwiseWorkers)
            return {WVKernelStatusCode::invalidConfiguration,"Pointwise worker count must be positive."};
        if (options.usesCompactSplitViews() && (options.horizontalSchedule!=WVRetainedHorizontalSchedule::streamingPrunedTile16 || !options.streamedNonlinear))
            return {WVKernelStatusCode::invalidConfiguration,"Compact views require the pruned horizontal and streamed nonlinear schedules."};
        const auto representation=options.usesCompactSplitViews() ? WVComplexRepresentation::split : WVComplexRepresentation::interleaved;
        const auto& g=source->geometry();
        if (g.transformClass!="WVTransformHydrostatic" || g.Nx<2 || g.Ny<2 || g.Nz<3 || !g.Nj || g.Nj>=g.Nz || !g.Nkl ||
            g.j.size()!=g.Nj || g.h_0.size()!=g.Nj || g.k.size()!=g.Nkl || g.l.size()!=g.Nkl || g.modes.size()!=g.Nkl ||
            g.z.size()!=g.Nz || g.N2.size()!=g.Nz || g.rho_nm0.size()!=g.Nz || g.dLnN2.size()!=g.Nz || g.z_int.size()!=g.Nz ||
            source->sourceIdentity().empty() || source->modeSetIdentity().empty())
            return {WVKernelStatusCode::invalidConfiguration,"Invalid hydrostatic scientific source."};
        auto candidate=std::unique_ptr<WVTransformHydrostaticKernel>(new WVTransformHydrostaticKernel);
        auto& c=*candidate; c.source_=std::move(source);
        c.S_=product(g.Nj,g.Nkl); c.H_=product(g.Nz,g.Nkl); c.R_=product(product(g.Nx,g.Ny),g.Nz);
        product(product(3,c.S_),sizeof(WVComplex64)); product(product(2,c.H_),sizeof(WVComplex64));
        if (3*c.S_>std::numeric_limits<std::size_t>::max()-2*c.H_) throw std::overflow_error("Hydrostatic spectral scratch extent overflow.");
        auto spectralElements=3*c.S_+2*c.H_;
        product(spectralElements,sizeof(WVComplex64));
        product(product(options.streamedNonlinear ? 6 : 10,c.R_),sizeof(double)); product(c.S_,sizeof(WVHydrostaticModeFactors));
        c.engineIdentifier_=engine->identifier(); c.engineLibraryIdentity_=engine->libraryIdentity();
        c.executionOptions_=options;
        WVRetainedHorizontalSpecification horizontal;
        auto status=c.source_->horizontalSpecification(g.Nz,representation,"hydrostatic-grid",horizontal); if (!status) return status;
        horizontal.schedule=options.horizontalSchedule; horizontal.outerWorkers=options.horizontalWorkers;
        status=WVRetainedHorizontalOperator::create(horizontal,std::move(engine),c.horizontal_); if (!status) return status;
        status=c.horizontal_->createWorkspace(c.horizontalWorkspace_); if (!status) return status;
        if (options.tiledNonlinear && options.sharedFieldGradients && options.streamedNonlinear) {
            status=c.horizontal_->prepareAdvection(*c.horizontalWorkspace_,3);
            if (!status && status.code!=WVKernelStatusCode::unsupportedOperation) return status;
            c.tiledNonlinearPrepared_=static_cast<bool>(status);
        }
        const WVStratifiedModalOperator operations[]={WVStratifiedModalOperator::reconstructF,WVStratifiedModalOperator::projectF,WVStratifiedModalOperator::reconstructG,WVStratifiedModalOperator::projectG};
        for (std::size_t i=0;i<4;++i) {
            const bool reconstruction=i%2==0; const std::string family=i<2 ? "F" : "G";
            const auto rows=reconstruction ? g.Nj : g.Nz, outRows=reconstruction ? g.Nz : g.Nj;
            WVComplexLayout in{rows,g.Nkl,1,rows,representation,family+(reconstruction ? "-modal" : "-grid"),c.source_->modeSetIdentity()};
            WVComplexLayout out{outRows,g.Nkl,1,outRows,representation,family+(reconstruction ? "-grid" : "-modal"),c.source_->modeSetIdentity()};
            std::unique_ptr<WVVerticalMatrixBackend> backend; status=factory(backend); if (!status) return status;
            status=c.source_->prepareVertical(operations[i],in,out,std::move(backend),c.vertical_[i]); if (!status) return status;
            status=c.vertical_[i]->createWorkspace(c.verticalWorkspace_[i]); if (!status) return status;
        }
        const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
        if (!std::isfinite(f) || f==0 || !std::isfinite(g.g) || !(g.g>0) || !std::isfinite(g.Lz) || !(g.Lz>0) || !std::isfinite(g.rho0) || !(g.rho0>0))
            return {WVKernelStatusCode::invalidConfiguration,"Invalid hydrostatic physical constants."};
        c.factors_.resize(c.S_);
        for (std::size_t i=0;i<c.S_;++i) {
            const auto j=i%g.Nj,mode=i/g.Nj; auto& a=c.factors_[i];
            const double k=g.k[mode],l=g.l[mode],K2=k*k+l*l,h=g.h_0[j];
            if (!std::isfinite(K2) || !std::isfinite(h) || !(h>0) || !std::isfinite(g.j[j]) || g.j[j]<0)
                return {WVKernelStatusCode::invalidConfiguration,"Invalid hydrostatic mode or equivalent depth."};
            const bool barotropic=g.j[j]==0;
            a.omega=std::sqrt(g.g*h*K2+f*f);
            if (!std::isfinite(a.omega) || !(a.omega>0)) return {WVKernelStatusCode::numericalFailure,"Hydrostatic frequency overflow."};
            a.wave=K2>0 && !barotropic; a.inertial=K2==0;
            a.geostrophic=K2>0; a.meanDensityAnomaly=K2==0 && !barotropic;
            if (a.wave) {
                const double K=std::sqrt(K2),alpha=std::atan2(l,k);
                a.UAp={std::cos(alpha),-f/a.omega*std::sin(alpha)};
                a.VAp={std::sin(alpha),f/a.omega*std::cos(alpha)};
                a.WAp={0,-K*h}; a.NAp=-K*h/a.omega;
                a.ApmD={0,-1/(2*K*h)}; a.ApmN=-a.omega/(2*K*h); a.waveEnergy=2*h;
            }
            if (a.inertial) { a.UAp={1,0}; a.VAp={0,1}; a.waveEnergy=barotropic ? g.Lz : h; }
            if (a.geostrophic) {
                const double denominator=K2+(barotropic ? 0 : f*f/(g.g*h));
                a.UA0={0,l/denominator}; a.VA0={0,-k/denominator};
                a.PA0=-(f/g.g)/denominator; a.NA0=barotropic ? 0 : a.PA0;
                a.A0Z=1; a.A0N=barotropic ? 0 : -f/h;
                a.psi=-1/denominator; a.qgpv=1; a.enstrophy=barotropic ? g.Lz : h;
                a.balancedEnergy=((barotropic ? g.Lz : h)*K2+(barotropic ? 0 : f*f/g.g))/(denominator*denominator);
            }
            if (a.meanDensityAnomaly) {
                a.NA0=1; a.PA0=1; a.A0N=1; a.A0Z=f*f/(2*h);
                a.psi=g.g/f; a.qgpv=-f/h; a.enstrophy=f*f/(2*h); a.balancedEnergy=g.g/2;
            }
            for (double x : {a.UAp.real,a.UAp.imag,a.VAp.real,a.VAp.imag,a.WAp.imag,a.NAp,a.UA0.imag,a.VA0.imag,
                a.NA0,a.PA0,a.ApmD.imag,a.ApmN,a.A0Z,a.A0N,a.waveEnergy,a.balancedEnergy,a.psi,a.qgpv,a.enstrophy})
                if (!std::isfinite(x)) return {WVKernelStatusCode::numericalFailure,"Hydrostatic coefficient factor overflow."};
        }
        if (c.tiledNonlinearPrepared_) {
            if (spectralElements>std::numeric_limits<std::size_t>::max()-c.H_)
                throw std::overflow_error("Tiled Hydrostatic spectral scratch extent overflow.");
            spectralElements+=c.H_;
            product(spectralElements,sizeof(WVComplex64));
        }
        c.spectralStorage_=std::make_unique<WVVariableComplexBuffer>(spectralElements,representation); c.phase_.resize(c.S_); c.real_.resize((options.streamedNonlinear ? 6 : 10)*c.R_);
        if (options.sharedFieldGradients) c.fieldCache_=std::make_unique<kernel_detail::WVPreparedFieldCache>(c.S_,c.H_,representation);
        c.pointwise_=std::make_unique<kernel_detail::WVPreparedModeExecutor>(std::min(options.pointwiseWorkers,c.R_));
        auto& s=c.storage_; s.sharedScientificBytes=c.source_->persistentBytes(); s.preparedBytes=c.horizontal_->persistentBytes();
        s.workspaceBytes=c.horizontalWorkspace_->persistentBytes()+sizeof(WVVariableComplexBuffer)+c.pointwise_->persistentBytes(); s.providerBytesLowerBound=c.horizontal_->providerBytesLowerBound(); s.planBytesLowerBound=c.horizontalWorkspace_->planBytesLowerBound();
        for (std::size_t i=0;i<4;++i) { s.preparedBytes+=c.vertical_[i]->persistentBytes(); s.workspaceBytes+=c.verticalWorkspace_[i]->persistentBytes(); }
        s.spectralScratchBytes=c.spectralStorage_->capacityBytes()+c.phase_.capacity()*sizeof(WVComplex64);
        c.baseSpectralScratchBytes_=s.spectralScratchBytes;
        s.realScratchBytes=c.real_.capacity()*sizeof(double); s.factorBytes=c.factors_.capacity()*sizeof(WVHydrostaticModeFactors);
        result=std::move(candidate); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Hydrostatic setup allocation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
      catch (const std::system_error& e) { return {WVKernelStatusCode::allocationFailure,e.what()}; }
}
const WVHydrostaticStorage& WVTransformHydrostaticKernel::storage() const noexcept {
    storage_.spectralScratchBytes=baseSpectralScratchBytes_+(fieldCache_ ? fieldCache_->capacityBytes() : 0);
    return storage_;
}
std::size_t WVTransformHydrostaticKernel::persistentBytes() const noexcept {
    const auto& s=storage(); return sizeof(*this)+s.sharedScientificBytes+s.preparedBytes+s.workspaceBytes+s.spectralScratchBytes+s.realScratchBytes+s.factorBytes;
}
WVKernelStatus WVTransformHydrostaticKernel::spectral(WVComplexConstView a) const {
    if (a.shape.rows!=geometry().Nj || a.shape.columns!=geometry().Nkl) return {WVKernelStatusCode::invalidShape,"Expected canonical [Nj,Nkl] coefficients."};
    if (!addressFits(a.data,S_*sizeof(WVComplex64),alignof(WVComplex64))) return {WVKernelStatusCode::invalidPointer,"Invalid coefficient storage."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::volume(WVRealVolumeConstView a,bool isSurface) const {
    const auto& g=geometry(); const auto nz=isSurface ? 1 : g.Nz;
    if (a.shape.first!=g.Nx || a.shape.second!=g.Ny || a.shape.third!=nz) return {WVKernelStatusCode::invalidShape,"Unexpected hydrostatic field shape."};
    if (!addressFits(a.data,(R_/g.Nz)*nz*sizeof(double),alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid field storage."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::disjoint(const void* a,std::size_t an,const void* b,std::size_t bn) const {
    if (overlap(a,an,b,bn)) return {WVKernelStatusCode::overlappingArrays,"Hydrostatic arrays must not overlap."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::outputs(WVMutableCoefficients a) const {
    for (auto x : {a.Ap,a.Am,a.A0}) { auto s=spectral({x.data,x.shape}); if (!s) return s; }
    auto s=disjoint(a.Ap.data,S_*sizeof(WVComplex64),a.Am.data,S_*sizeof(WVComplex64)); if (!s) return s;
    s=disjoint(a.Ap.data,S_*sizeof(WVComplex64),a.A0.data,S_*sizeof(WVComplex64)); if (!s) return s;
    return disjoint(a.Am.data,S_*sizeof(WVComplex64),a.A0.data,S_*sizeof(WVComplex64));
}
WVKernelStatus WVTransformHydrostaticKernel::mutableOutputOutsidePreparedState(WVMutableCoefficients a) const {
    for (const auto output:{a.Ap,a.Am,a.A0}) {
        auto status=mutableOutputOutsidePreparedState(output.data,S_*sizeof(WVComplex64));
        if (!status) return status;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::mutableOutputOutsidePreparedState(
    const void* output,std::size_t bytes) const {
    if (!stateEvaluationActive_) return WVKernelStatus::ok();
    for (std::size_t stateIndex=0;stateIndex<preparedStateViewCount_;++stateIndex)
        for (const auto input:{preparedStateViews_[stateIndex].coefficients.Ap,preparedStateViews_[stateIndex].coefficients.Am,preparedStateViews_[stateIndex].coefficients.A0})
            if (overlap(output,bytes,input.data,S_*sizeof(WVComplex64)))
                return {WVKernelStatusCode::overlappingArrays,"Mutable output overlaps an active immutable hydrostatic state view."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::coefficients(const WVCoefficients& a) const {
    for (auto x : {a.Ap,a.Am,a.A0}) {
        auto s=spectral(x); if (!s) return s;
        for (std::size_t i=0;i<S_;++i) if (!std::isfinite(x.data[i].real) || !std::isfinite(x.data[i].imag))
            return {WVKernelStatusCode::numericalFailure,"Nonfinite hydrostatic coefficient."};
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::stateContents(const WVState& a) const {
    auto s=coefficients(a.coefficients); if (!s) return s;
    if (!std::isfinite(a.t) || !std::isfinite(a.t0) || !std::isfinite(a.t-a.t0)) return {WVKernelStatusCode::invalidConfiguration,"Nonfinite hydrostatic time or elapsed time."};
    for (const auto& f:factors_) if (!std::isfinite(f.omega*(a.t-a.t0))) return {WVKernelStatusCode::numericalFailure,"Hydrostatic phase overflow."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::state(const WVState& a) {
    ++metrics_.stateValidationCount;
    return stateContents(a);
}
WVKernelStatus WVTransformHydrostaticKernel::preparePhase(
    double t,double t0,const WVState* validatedState) {
    const bool exactValidatedState=validatedState!=nullptr &&
        stateEvaluationActive_ && matchesStateEvaluation(*validatedState) &&
        validatedState->t==t && validatedState->t0==t0;
    if (!exactValidatedState) {
        if (!std::isfinite(t) || !std::isfinite(t0) || !std::isfinite(t-t0)) return {WVKernelStatusCode::invalidConfiguration,"Nonfinite hydrostatic time or elapsed time."};
        for (const auto& f:factors_) if (!std::isfinite(f.omega*(t-t0))) return {WVKernelStatusCode::numericalFailure,"Hydrostatic phase overflow."};
    }
    ++metrics_.phasePreparationCount;
    pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
        for (std::size_t i=begin;i<end;++i) {
            const double a=factors_[i].omega*(t-t0);
            phase_[i]={std::cos(a),std::sin(a)};
        }
    });
    return WVKernelStatus::ok();
}
bool WVTransformHydrostaticKernel::matchesStateEvaluation(const WVState& a) const noexcept {
    const auto sameView=[](WVComplexConstView x,WVComplexConstView y) {
        return x.data==y.data && x.shape.rows==y.shape.rows && x.shape.columns==y.shape.columns;
    };
    if (!stateEvaluationActive_ || a.t!=preparedState_.t || a.t0!=preparedState_.t0) return false;
    for (std::size_t i=0;i<preparedStateViewCount_;++i) if (
        sameView(a.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(a.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(a.coefficients.A0,preparedStateViews_[i].coefficients.A0)) return true;
    return false;
}
WVKernelStatus WVTransformHydrostaticKernel::validateStateEvaluation(
    const WVState& a) const noexcept {
    if (!matchesStateEvaluation(a))
        return {WVKernelStatusCode::invalidConfiguration,
            "State does not belong to the active Hydrostatic evaluation."};
    return WVKernelStatus::ok();
}
std::size_t WVTransformHydrostaticKernel::stateEvaluationComponent(
    const WVState& a) const noexcept {
    const auto sameView=[](WVComplexConstView x,WVComplexConstView y) {
        return x.data==y.data && x.shape.rows==y.shape.rows && x.shape.columns==y.shape.columns;
    };
    if (!stateEvaluationActive_) return 0;
    for (std::size_t i=0;i<preparedStateViewCount_;++i) if (
        sameView(a.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(a.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(a.coefficients.A0,preparedStateViews_[i].coefficients.A0))
        return preparedStateComponents_[i];
    return 0;
}
std::size_t WVTransformHydrostaticKernel::fieldPreparationView(const WVCoefficients& a) const noexcept {
    for (std::size_t i=0;i<preparedStateViewCount_;++i) {
        const auto& b=preparedStateViews_[i].coefficients;
        if (a.Ap.data==b.Ap.data && a.Am.data==b.Am.data && a.A0.data==b.A0.data)
            return preparedStateViewIds_[i];
    }
    return 0; // Standalone operation; its cache is cleared on entry and exit.
}
WVKernelStatus WVTransformHydrostaticKernel::validateStateForCall(const WVState& a) {
    if (!stateEvaluationActive_) return state(a);
    if (!matchesStateEvaluation(a))
        return {WVKernelStatusCode::invalidConfiguration,"State does not match the active hydrostatic evaluation."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::preparePhaseForCall(const WVState& a) {
    if (stateEvaluationActive_) {
        if (!matchesStateEvaluation(a))
            return {WVKernelStatusCode::invalidConfiguration,"State does not match the active hydrostatic evaluation."};
        return WVKernelStatus::ok();
    }
    if (fieldCache_) fieldCache_->clear();
    return preparePhase(a.t,a.t0);
}
WVKernelStatus WVTransformHydrostaticKernel::prepareProjectionPhaseForCall(double t,double t0) {
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0)
            return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active hydrostatic evaluation."};
        return WVKernelStatus::ok();
    }
    return preparePhase(t,t0);
}
WVKernelStatus WVTransformHydrostaticKernel::beginStateEvaluation(const WVState& a) {
    return beginStateEvaluation(a,nullptr);
}
WVKernelStatus WVTransformHydrostaticKernel::beginStateEvaluation(const WVState& a,const void* evaluationOwner) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (stateEvaluationActive_)
        return {WVKernelStatusCode::reentrantExecution,"Hydrostatic state evaluation is already active."};
    auto s=state(a); if (!s) return s;
    if (fieldCache_) fieldCache_->clear();
    preparedStateViewIds_={}; nextStateViewId_=1;
    preparedState_=a;
    preparedStateViews_[0]=a;
    preparedStateComponents_[0]=0;
    preparedStateViewCount_=1;
    preparedStateOwner_=evaluationOwner;
    stateEvaluationActive_=true;
    // state(a) validated this exact registered view, including every phase
    // product, so phase preparation must not repeat the same finite scan.
    s=preparePhase(a.t,a.t0,&a);
    if (!s) {
        preparedState_={};
        preparedStateViews_[0]={};
        preparedStateComponents_[0]=0;
        preparedStateViewCount_=0;
        preparedStateOwner_=nullptr;
        stateEvaluationActive_=false;
        return s;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::addStateEvaluationView(
    const WVState& a,const void* evaluationOwner,std::size_t componentIdentity) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No hydrostatic state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Hydrostatic state view owner does not match the active evaluation."};
    if (a.t!=preparedState_.t || a.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Additional hydrostatic state view must use the active evaluation times."};
    if (componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"Hydrostatic component identity is out of range."};
    if (matchesStateEvaluation(a))
        return stateEvaluationComponent(a)==componentIdentity ? WVKernelStatus::ok() :
            WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                "Hydrostatic state view is already registered with another component identity."};
    if (preparedStateViewCount_==preparedStateViews_.size())
        return {WVKernelStatusCode::invalidConfiguration,"Hydrostatic state evaluation view capacity exceeded."};
    auto s=state(a); if (!s) return s;
    if (nextStateViewId_==std::numeric_limits<std::size_t>::max())
        return {WVKernelStatusCode::sizeOverflow,"Hydrostatic state view generation overflow."};
    preparedStateViewIds_[preparedStateViewCount_]=nextStateViewId_++;
    preparedStateViews_[preparedStateViewCount_]=a;
    preparedStateComponents_[preparedStateViewCount_++]=componentIdentity;
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::removeStateEvaluationView(
    const WVState& a,const void* evaluationOwner,std::size_t componentIdentity) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No hydrostatic state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Hydrostatic state view owner does not match the active evaluation."};
    if (componentIdentity==0 || componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"The primary hydrostatic state view cannot be removed."};
    if (a.t!=preparedState_.t || a.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Removed hydrostatic state view must use the active evaluation times."};
    const auto sameView=[](WVComplexConstView x,WVComplexConstView y) {
        return x.data==y.data && x.shape.rows==y.shape.rows && x.shape.columns==y.shape.columns;
    };
    for(std::size_t i=1;i<preparedStateViewCount_;++i) if (
        preparedStateComponents_[i]==componentIdentity &&
        sameView(a.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(a.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(a.coefficients.A0,preparedStateViews_[i].coefficients.A0)) {
        if (fieldCache_) fieldCache_->invalidateView(preparedStateViewIds_[i]);
        for(std::size_t j=i+1;j<preparedStateViewCount_;++j) {
            preparedStateViewIds_[j-1]=preparedStateViewIds_[j];
            preparedStateViews_[j-1]=preparedStateViews_[j];
            preparedStateComponents_[j-1]=preparedStateComponents_[j];
        }
        --preparedStateViewCount_;
        preparedStateViews_[preparedStateViewCount_]={};
        preparedStateComponents_[preparedStateViewCount_]=0;
        return WVKernelStatus::ok();
    }
    return {WVKernelStatusCode::invalidConfiguration,
        "Hydrostatic state view is not registered for this component."};
}
WVKernelStatus WVTransformHydrostaticKernel::endStateEvaluation() {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No hydrostatic state evaluation is active."};
    if (fieldCache_) fieldCache_->clear();
    preparedStateViewIds_={};
    preparedState_={};
    for (auto& stateView:preparedStateViews_) stateView={};
    preparedStateComponents_={};
    preparedStateViewCount_=0;
    preparedStateOwner_=nullptr;
    stateEvaluationActive_=false;
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::validateFluxOutput(const WVState& a,const WVFlux& b) {
    auto s=validateStateForCall(a); if (!s) return s;
    WVMutableCoefficients target{b.Fp,b.Fm,b.F0}; s=outputs(target); if (!s) return s;
    for (auto inputView : {a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0})
        for (auto outputView : {b.Fp,b.Fm,b.F0}) {
            s=disjoint(inputView.data,S_*sizeof(WVComplex64),outputView.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
    return mutableOutputOutsidePreparedState(target);
}
WVKernelStatus WVTransformHydrostaticKernel::preparedPhase(const WVState& a,WVComplexConstView& result) {
    result={};
    auto s=validateStateForCall(a); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=preparePhaseForCall(a); if (!s) return s;
    result={phase_.data(),spectralShape()};
    return WVKernelStatus::ok();
}
WVComplexOutput WVTransformHydrostaticKernel::modalView(std::size_t slot) { return spectralStorage_->output(slot*S_,S_); }
WVComplexOutput WVTransformHydrostaticKernel::gridView(std::size_t slot) { return spectralStorage_->output(3*S_+slot*H_,H_); }
WVKernelStatus WVTransformHydrostaticKernel::vertical(std::size_t operation,WVComplexInput a,WVComplexOutput b) {
    ++metrics_.verticalOperatorExecutionCount;
    return vertical_[operation]->execute(*verticalWorkspace_[operation],a,b);
}
WVKernelStatus WVTransformHydrostaticKernel::verticalColumn(std::size_t operation,WVComplexInput a,
    WVComplexOutput b,std::size_t retainedColumn) {
    ++metrics_.verticalOperatorExecutionCount;
    return vertical_[operation]->executeColumn(*verticalWorkspace_[operation],a,b,retainedColumn);
}
WVKernelStatus WVTransformHydrostaticKernel::project(const double* a,WVComplexOutput b,WVHydrostaticFamily family) {
    auto s=horizontal_->forward(*horizontalWorkspace_,{a,R_*sizeof(double)},gridView()); if (!s) return s;
    return vertical(family==WVHydrostaticFamily::F ? 1 : 3,gridView().input(),b);
}
WVKernelStatus WVTransformHydrostaticKernel::transformToSpatial(WVComplexConstView a,WVHydrostaticFamily family,WVRealVolumeView b) {
    if (!valid(family)) return unsupported();
    auto s=spectral(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,S_*sizeof(WVComplex64),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    auto values=input(a.data,S_);
    if (executionOptions_.usesCompactSplitViews()) {
        const auto target=modalView(); for (std::size_t i=0;i<S_;++i) write(target,i,a.data[i]);
        values=target.input();
    }
    s=vertical(family==WVHydrostaticFamily::F ? 0 : 2,values,gridView()); if (!s) return s;
    return horizontal_->inverse(*horizontalWorkspace_,gridView().input(),{b.data,R_*sizeof(double)});
}
WVKernelStatus WVTransformHydrostaticKernel::transformFromSpatial(WVRealVolumeConstView a,WVHydrostaticFamily family,WVComplexView b) {
    if (!valid(family)) return unsupported();
    auto s=volume(a); if (!s) return s; s=spectral({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,S_*sizeof(WVComplex64)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews()) return project(a.data,output(b.data,S_),family);
    s=project(a.data,modalView(),family); if (!s) return s;
    for (std::size_t i=0;i<S_;++i) b.data[i]=read(modalView().input(),i);
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::applyVertical(WVStratifiedModalOperator operation,
    WVComplexConstView a,WVComplexView b) {
    std::size_t index=0; bool inputModal=false,outputModal=false;
    if (!verticalOperation(operation,index,inputModal,outputModal)) return unsupported();
    const auto& g=geometry(); const auto inputRows=inputModal ? g.Nj : g.Nz;
    const auto outputRows=outputModal ? g.Nj : g.Nz;
    if (a.shape.rows!=inputRows || a.shape.columns!=g.Nkl ||
        b.shape.rows!=outputRows || b.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Unexpected hydrostatic vertical-transform shape."};
    const auto inputCount=inputRows*g.Nkl,outputCount=outputRows*g.Nkl;
    if (!addressFits(a.data,inputCount*sizeof(WVComplex64),alignof(WVComplex64)) ||
        !addressFits(b.data,outputCount*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid hydrostatic vertical-transform storage."};
    auto status=mutableOutputOutsidePreparedState(b.data,outputCount*sizeof(WVComplex64)); if (!status) return status;
    status=disjoint(a.data,inputCount*sizeof(WVComplex64),b.data,outputCount*sizeof(WVComplex64)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews()) return vertical(index,input(a.data,inputCount),output(b.data,outputCount));
    const auto in=inputModal ? modalView() : gridView();
    const auto out=outputModal ? modalView(inputModal ? 1 : 0) : gridView(inputModal ? 0 : 1);
    copy(a.data,in,inputCount); status=vertical(index,in.input(),out); if (!status) return status;
    copy(out.input(),b.data,outputCount); return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::applyVerticalColumn(WVStratifiedModalOperator operation,
    std::size_t retainedColumn,WVComplexConstView a,WVComplexView b) {
    std::size_t index=0; bool inputModal=false,outputModal=false;
    if (!verticalOperation(operation,index,inputModal,outputModal)) return unsupported();
    const auto& g=geometry(); const auto inputRows=inputModal ? g.Nj : g.Nz;
    const auto outputRows=outputModal ? g.Nj : g.Nz;
    if (retainedColumn>=g.Nkl || a.shape.rows!=inputRows || a.shape.columns!=1 ||
        b.shape.rows!=outputRows || b.shape.columns!=1)
        return {WVKernelStatusCode::invalidShape,"Unexpected hydrostatic vertical-column shape or retained column."};
    if (!addressFits(a.data,inputRows*sizeof(WVComplex64),alignof(WVComplex64)) ||
        !addressFits(b.data,outputRows*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid hydrostatic vertical-column storage."};
    auto status=mutableOutputOutsidePreparedState(b.data,outputRows*sizeof(WVComplex64)); if (!status) return status;
    status=disjoint(a.data,inputRows*sizeof(WVComplex64),b.data,outputRows*sizeof(WVComplex64)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    const auto in=inputModal ? modalView() : gridView();
    const auto out=outputModal ? modalView(inputModal ? 1 : 0) : gridView(inputModal ? 0 : 1);
    const auto inputOffset=retainedColumn*inputRows,outputOffset=retainedColumn*outputRows;
    for (std::size_t row=0;row<inputRows;++row) write(in,inputOffset+row,a.data[row]);
    status=verticalColumn(index,in.input(),out,retainedColumn); if (!status) return status;
    for (std::size_t row=0;row<outputRows;++row) b.data[row]=read(out.input(),outputOffset+row);
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::horizontalForward(WVRealVolumeConstView a,WVComplexView b) {
    const auto& g=geometry(); auto status=volume(a); if (!status) return status;
    if (b.shape.rows!=g.Nz || b.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Expected hydrostatic Fourier [Nz,Nkl] output."};
    if (!addressFits(b.data,H_*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid hydrostatic Fourier output storage."};
    status=mutableOutputOutsidePreparedState(b.data,H_*sizeof(WVComplex64)); if (!status) return status;
    status=disjoint(a.data,R_*sizeof(double),b.data,H_*sizeof(WVComplex64)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews())
        return horizontal_->forward(*horizontalWorkspace_,{a.data,R_*sizeof(double)},output(b.data,H_));
    auto grid=gridView(); status=horizontal_->forward(*horizontalWorkspace_,{a.data,R_*sizeof(double)},grid); if (!status) return status;
    copy(grid.input(),b.data,H_); return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::horizontalInverse(WVComplexConstView a,WVRealVolumeView b) {
    const auto& g=geometry(); auto status=volume({b.data,b.shape}); if (!status) return status;
    if (a.shape.rows!=g.Nz || a.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Expected hydrostatic Fourier [Nz,Nkl] input."};
    if (!addressFits(a.data,H_*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid hydrostatic Fourier input storage."};
    status=mutableOutputOutsidePreparedState(b.data,R_*sizeof(double)); if (!status) return status;
    status=disjoint(a.data,H_*sizeof(WVComplex64),b.data,R_*sizeof(double)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    auto grid=gridView(); copy(a.data,grid,H_);
    for (std::size_t mode=0;mode<g.Nkl;++mode)
        if (selfConjugate(g.modes[mode].k,g.Nx) && selfConjugate(g.modes[mode].l,g.Ny))
            for (std::size_t z=0;z<g.Nz;++z) {
                const auto i=z+g.Nz*mode;
                const auto value=read(grid.input(),i);
                write(grid,i,{value.real,0});
            }
    return horizontal_->inverse(*horizontalWorkspace_,grid.input(),{b.data,R_*sizeof(double)});
}
WVKernelStatus WVTransformHydrostaticKernel::projectFields(const double* u,const double* v,const double* eta,WVMutableCoefficients b) {
    const auto U=modalView(),V=modalView(1),N=modalView(2);
    auto s=project(u,U,WVHydrostaticFamily::F); if (!s) return s;
    s=project(v,V,WVHydrostaticFamily::F); if (!s) return s;
    s=project(eta,N,WVHydrostaticFamily::G); if (!s) return s;
    return projectedFieldsToCoefficients(U.input(),V.input(),N.input(),b);
}
WVKernelStatus WVTransformHydrostaticKernel::projectedFieldsToCoefficients(
    WVComplexInput U,WVComplexInput V,WVComplexInput N,WVMutableCoefficients b) {
    pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
      for (std::size_t i=begin;i<end;++i) {
        const auto& f=factors_[i]; const auto mode=i/geometry().Nj; const double k=geometry().k[mode],l=geometry().l[mode];
        const auto zeta=subtract(multiply(read(V,i),{0,k}),multiply(read(U,i),{0,l}));
        const auto delta=scale(add(multiply(read(U,i),{0,k}),multiply(read(V,i),{0,l})),geometry().h_0[i%geometry().Nj]);
        const auto a0=add(scale(zeta,f.A0Z),scale(read(N,i),f.A0N));
        const auto density=scale(subtract(read(N,i),scale(a0,f.NA0)),f.ApmN),divergence=multiply(delta,f.ApmD);
        auto ap=add(divergence,density),am=subtract(divergence,density);
        if (f.inertial) { ap=scale(subtract(read(U,i),multiply(read(V,i),{0,1})),.5); am=conjugate(ap); }
        b.Ap.data[i]=multiply(ap,conjugate(phase_[i])); b.Am.data[i]=multiply(am,phase_[i]); b.A0.data[i]=a0;
    }
    });
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::transformUVEtaToWaveVortex(WVRealVolumeConstView u,WVRealVolumeConstView v,WVRealVolumeConstView eta,double t,double t0,WVMutableCoefficients b) {
    auto s=outputs(b); if (!s) return s; s=mutableOutputOutsidePreparedState(b); if (!s) return s;
    for (auto a : {u,v,eta}) { s=volume(a); if (!s) return s; for (auto out : {b.Ap,b.Am,b.A0}) { s=disjoint(a.data,R_*sizeof(double),out.data,S_*sizeof(WVComplex64)); if (!s) return s; } }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=prepareProjectionPhaseForCall(t,t0); if (!s) return s;
    return projectFields(u.data,v.data,eta.data,b);
}

WVKernelStatus WVTransformHydrostaticKernel::reconstruct(const WVCoefficients& a,WVHydrostaticField field,
    WVHydrostaticDerivative derivative,WVHydrostaticComponent component,double* b,bool countPrimary,
    std::size_t metricComponent,const WVRealOutputConsumer* consumer,
    WVComplexInput* preparedSpectrum,const WVComplexOutput* destination) {
    if (countPrimary && !preparedSpectrum) {
        if (metricComponent>=5) metricComponent=static_cast<std::size_t>(component);
        ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(field)];
        ++metrics_.reconstructionCount[static_cast<std::size_t>(field)]
            [static_cast<std::size_t>(derivative)][metricComponent];
    }
    const auto& g=geometry();
    const auto reconstructionGrid=destination ? *destination : gridView();
    if (field==WVHydrostaticField::zetaX || field==WVHydrostaticField::zetaY) {
        const bool x=field==WVHydrostaticField::zetaX;
        auto s=reconstruct(a,x ? WVHydrostaticField::w : WVHydrostaticField::u,x ? WVHydrostaticDerivative::y : WVHydrostaticDerivative::z,component,b,countPrimary,metricComponent); if (!s) return s;
        auto* auxiliary=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 8)*R_;
        s=reconstruct(a,x ? WVHydrostaticField::v : WVHydrostaticField::w,x ? WVHydrostaticDerivative::z : WVHydrostaticDerivative::x,component,auxiliary,countPrimary,metricComponent); if (!s) return s;
        for (std::size_t i=0;i<R_;++i) b[i]-=auxiliary[i];
        return WVKernelStatus::ok();
    }
    const bool density=field==WVHydrostaticField::rhoE || field==WVHydrostaticField::rhoTotal;
    const bool totalDensity=field==WVHydrostaticField::rhoTotal;
    if (density) field=WVHydrostaticField::eta;
    if (field==WVHydrostaticField::ssh) field=WVHydrostaticField::pi;
    if (field==WVHydrostaticField::ssu) field=WVHydrostaticField::u;
    if (field==WVHydrostaticField::ssv) field=WVHydrostaticField::v;
    const bool G=field==WVHydrostaticField::w || field==WVHydrostaticField::eta;
    const bool dz=derivative==WVHydrostaticDerivative::z;
    if (field<WVHydrostaticField::u ||
        (field>WVHydrostaticField::qgpv && field!=WVHydrostaticField::zetaZ)) return unsupported();
    // Independent coefficients preserve their arithmetic order. The prepared
    // executor joins before the vertical operator reuses the modal storage.
    kernel_detail::WVPreparedFieldCache::Entry* prepared=nullptr;
    if (fieldCache_ && countPrimary) {
        auto status=fieldCache_->acquire({fieldPreparationView(a),static_cast<std::size_t>(field),
            static_cast<std::size_t>(component)},prepared);
        if (!status) return status;
    }
    const auto modal=prepared ? prepared->modal(0,S_) : modalView();
    const auto assemble=[&](std::size_t begin,std::size_t end) {
        for (std::size_t i=begin;i<end;++i) {
            const auto& f=factors_[i]; const auto mode=i/g.Nj; WVComplex64 p{},m{},z{};
            switch(field) {
                case WVHydrostaticField::u: p=f.UAp; m=conjugate(p); z=f.UA0; break;
                case WVHydrostaticField::v: p=f.VAp; m=conjugate(p); z=f.VA0; break;
                case WVHydrostaticField::w: p=f.WAp; m=p; break;
                case WVHydrostaticField::eta: p={f.NAp,0}; m={-f.NAp,0}; z={f.NA0,0}; break;
                case WVHydrostaticField::pi: p={f.NAp,0}; m={-f.NAp,0}; z={f.PA0,0}; break;
                case WVHydrostaticField::p: p={g.rho0*g.g*f.NAp,0}; m=scale(p,-1); z={g.rho0*g.g*f.PA0,0}; break;
                case WVHydrostaticField::psi: z={f.psi,0}; break;
                case WVHydrostaticField::qgpv: z={f.qgpv,0}; break;
                case WVHydrostaticField::zetaZ:
                    p=subtract(multiply(f.VAp,{0,g.k[mode]}),multiply(f.UAp,{0,g.l[mode]}));
                    m=subtract(multiply(conjugate(f.VAp),{0,g.k[mode]}),multiply(conjugate(f.UAp),{0,g.l[mode]}));
                    z=subtract(multiply(f.VA0,{0,g.k[mode]}),multiply(f.UA0,{0,g.l[mode]})); break;
                default: break; // Field was validated before dispatch.
            }
            WVComplex64 value{};
            if (selected(f,component,true)) value=add(multiply(multiply(p,a.Ap.data[i]),phase_[i]),multiply(multiply(m,a.Am.data[i]),conjugate(phase_[i])));
            if (selected(f,component,false)) value=add(value,multiply(z,a.A0.data[i]));
            if (!prepared && derivative==WVHydrostaticDerivative::x) value=multiply(value,{0,g.k[mode]});
            if (!prepared && derivative==WVHydrostaticDerivative::y) value=multiply(value,{0,g.l[mode]});
            if (!prepared && dz && G) value=scale(value,1/g.h_0[i%g.Nj]);
            write(modal,i,value);
        }
    };
    if (!prepared || !prepared->modalReady) {
        ++metrics_.coefficientAssemblyCount;
        pointwise_->execute(S_,assemble);
        if (prepared) prepared->modalReady=true;
    }
    WVKernelStatus s;
    auto inverseInput=reconstructionGrid.input();
    if (prepared && !dz) {
        if (!prepared->gridReady) {
            ++metrics_.verticalPreparationCount;
            s=vertical(G ? 2 : 0,modal.input(),prepared->grid()); if (!s) return s;
            prepared->gridReady=true;
        } else ++metrics_.horizontalSpectrumReuseCount;
        inverseInput=prepared->grid().input();
        if (derivative!=WVHydrostaticDerivative::value) {
            pointwise_->execute(H_,[&](std::size_t begin,std::size_t end) {
                for (std::size_t i=begin;i<end;++i) {
                    const auto mode=i/g.Nz;
                    write(reconstructionGrid,i,multiply(read(prepared->grid().input(),i),
                        {0,derivative==WVHydrostaticDerivative::x ? g.k[mode] : g.l[mode]}));
                }
            });
            inverseInput=reconstructionGrid.input();
        }
    } else {
        if (prepared) {
            pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
                for (std::size_t i=begin;i<end;++i)
                    write(modalView(),i,G ? scale(read(modal.input(),i),1/g.h_0[i%g.Nj]) : read(modal.input(),i));
            });
            ++metrics_.preparedVerticalDerivativeCount;
        }
        ++metrics_.verticalPreparationCount;
        s=vertical(G!=dz ? 2 : 0,modalView().input(),reconstructionGrid); if (!s) return s;
    }
    if (dz && !G) for (std::size_t mode=0;mode<g.Nkl;++mode) for (std::size_t z=0;z<g.Nz;++z)
        write(reconstructionGrid,z+g.Nz*mode,scale(read(reconstructionGrid.input(),z+g.Nz*mode),-g.N2[z]/g.g));
    if (preparedSpectrum) { *preparedSpectrum=inverseInput; return WVKernelStatus::ok(); }
    s=consumer ? horizontal_->inverseAndConsume(*horizontalWorkspace_,inverseInput,{b,R_*sizeof(double)},*consumer) :
        horizontal_->inverse(*horizontalWorkspace_,inverseInput,{b,R_*sizeof(double)}); if (!s) return s;
    if (density) {
        double* eta=nullptr;
        if (dz) {
            eta=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 8)*R_;
            s=reconstruct(a,WVHydrostaticField::eta,WVHydrostaticDerivative::value,component,eta,countPrimary,metricComponent); if (!s) return s;
        }
        const auto plane=R_/g.Nz;
        for (std::size_t z=0;z<g.Nz;++z) for (std::size_t xy=0;xy<plane;++xy) {
            const auto i=xy+plane*z;
            b[i]=(g.rho0/g.g)*g.N2[z]*(b[i]+(dz ? g.dLnN2[z]*eta[i] : 0));
            if (totalDensity && derivative==WVHydrostaticDerivative::value) b[i]+=g.rho_nm0[z];
            if (totalDensity && dz) b[i]-=(g.rho0/g.g)*g.N2[z];
        }
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::transformStateField(const WVState& a,WVHydrostaticField field,WVRealVolumeView b,WVHydrostaticDerivative derivative,WVHydrostaticComponent component) {
    if (field<WVHydrostaticField::u || field>WVHydrostaticField::ssv || derivative<WVHydrostaticDerivative::value || derivative>WVHydrostaticDerivative::z || !valid(component) ||
        ((field==WVHydrostaticField::zetaX || field==WVHydrostaticField::zetaY) && derivative!=WVHydrostaticDerivative::value) ||
        (field==WVHydrostaticField::rhoTotal && component!=WVHydrostaticComponent::all)) return unsupported();
    auto s=validateStateForCall(a); if (!s) return s; s=volume({b.data,b.shape},surface(field)); if (!s) return s;
    const auto bytes=(surface(field) ? R_/geometry().Nz : R_)*sizeof(double);
    for (auto x : {a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) { s=disjoint(x.data,S_*sizeof(WVComplex64),b.data,bytes); if (!s) return s; }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    auto* fullField=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 9)*R_;
    const auto metricComponent=component==WVHydrostaticComponent::all ?
        stateEvaluationComponent(a) : static_cast<std::size_t>(component);
    s=reconstruct(a.coefficients,field,derivative,component,
        surface(field) ? fullField : b.data,true,metricComponent); if (!s) return s;
    if (surface(field)) { const auto plane=R_/geometry().Nz; std::copy_n(fullField+R_-plane,plane,b.data); }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::transformCoefficientTendencyToUVEta(
    const WVState& a,WVRealFieldBundleView& fields) {
    const auto& g=geometry();
    if (fields.shape.first!=g.Nx || fields.shape.second!=g.Ny ||
        fields.shape.third!=g.Nz || fields.shape.fourth!=3)
        return {WVKernelStatusCode::invalidShape,"Coefficient tendency requires [Nx,Ny,Nz,3] fields."};
    ++metrics_.derivedValidationCount;
    auto status=stateContents(a); if (!status) return status;
    const auto bytes=3*R_*sizeof(double);
    if (!addressFits(fields.data,bytes,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid coefficient tendency field storage."};
    for (const auto input:{a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) {
        status=disjoint(input.data,S_*sizeof(WVComplex64),fields.data,bytes);
        if (!status) return status;
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    status=prepareProjectionPhaseForCall(a.t,a.t0); if (!status) return status;
    const WVHydrostaticField names[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::eta};
    for (std::size_t channel=0;channel<3;++channel) {
        status=reconstruct(a.coefficients,names[channel],WVHydrostaticDerivative::value,
            WVHydrostaticComponent::all,fields.data+channel*R_,false);
        if (!status) return status;
        ++metrics_.tendencyReconstructionCount[channel];
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::combinePreparedHorizontalVorticity(
    WVHydrostaticField field,WVRealVolumeConstView first,WVRealVolumeConstView second,
    WVRealVolumeView output,WVHydrostaticComponent component) {
    if ((field!=WVHydrostaticField::zetaX && field!=WVHydrostaticField::zetaY) ||
        !valid(component))
        return unsupported();
    auto status=volume(first); if (!status) return status;
    status=volume(second); if (!status) return status;
    status=volume({output.data,output.shape}); if (!status) return status;
    if(first.data!=output.data) {
        status=disjoint(first.data,R_*sizeof(double),output.data,R_*sizeof(double)); if (!status) return status;
    }
    if(second.data!=output.data) {
        status=disjoint(second.data,R_*sizeof(double),output.data,R_*sizeof(double)); if (!status) return status;
    }
    for (std::size_t i=0;i<R_;++i) output.data[i]=first.data[i]-second.data[i];
    ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(field)];
    ++metrics_.reconstructionCount[static_cast<std::size_t>(field)]
        [static_cast<std::size_t>(WVHydrostaticDerivative::value)]
        [static_cast<std::size_t>(component)];
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::combinePreparedDensityZDerivative(
    WVHydrostaticField field,WVRealVolumeConstView etaZ,WVRealVolumeConstView eta,
    WVRealVolumeView output,WVHydrostaticComponent component) {
    if ((field!=WVHydrostaticField::rhoE && field!=WVHydrostaticField::rhoTotal) ||
        !valid(component) ||
        (field==WVHydrostaticField::rhoTotal && component!=WVHydrostaticComponent::all))
        return unsupported();
    auto status=volume(etaZ); if (!status) return status;
    status=volume(eta); if (!status) return status;
    status=volume({output.data,output.shape}); if (!status) return status;
    if(etaZ.data!=output.data) {
        status=disjoint(etaZ.data,R_*sizeof(double),output.data,R_*sizeof(double)); if (!status) return status;
    }
    if(eta.data!=output.data) {
        status=disjoint(eta.data,R_*sizeof(double),output.data,R_*sizeof(double)); if (!status) return status;
    }
    const auto& g=geometry(); const auto plane=R_/g.Nz;
    for (std::size_t z=0;z<g.Nz;++z) for (std::size_t xy=0;xy<plane;++xy) {
        const auto i=xy+plane*z; const double scale=(g.rho0/g.g)*g.N2[z];
        output.data[i]=scale*(etaZ.data[i]+g.dLnN2[z]*eta.data[i]);
        if (field==WVHydrostaticField::rhoTotal) output.data[i]-=scale;
    }
    ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(field)];
    ++metrics_.reconstructionCount[static_cast<std::size_t>(field)]
        [static_cast<std::size_t>(WVHydrostaticDerivative::z)]
        [static_cast<std::size_t>(component)];
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::evolveCoefficients(const WVState& a,WVMutableCoefficients b) {
    auto s=validateStateForCall(a); if (!s) return s; s=outputs(b); if (!s) return s; s=mutableOutputOutsidePreparedState(b); if (!s) return s;
    const WVComplexConstView inputs[]={a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}; const WVComplexView targets[]={b.Ap,b.Am,b.A0};
    for (std::size_t i=0;i<3;++i) for (std::size_t j=0;j<3;++j) if (i!=j || inputs[i].data!=targets[j].data) {
        s=disjoint(inputs[i].data,S_*sizeof(WVComplex64),targets[j].data,S_*sizeof(WVComplex64)); if (!s) return s;
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=preparePhaseForCall(a); if (!s) return s;
    for (std::size_t i=0;i<S_;++i) {
        b.Ap.data[i]=multiply(a.coefficients.Ap.data[i],phase_[i]); b.Am.data[i]=multiply(a.coefficients.Am.data[i],conjugate(phase_[i])); b.A0.data[i]=a.coefficients.A0.data[i];
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::constrainCoefficients(WVMutableCoefficients a) const {
    auto s=outputs(a); if (!s) return s; s=coefficients(view(a)); if (!s) return s; s=mutableOutputOutsidePreparedState(a); if (!s) return s;
    for (std::size_t i=0;i<S_;++i) {
        const auto& f=factors_[i];
        if (!f.wave && !f.inertial) { a.Ap.data[i]={}; a.Am.data[i]={}; }
        if (f.inertial) a.Am.data[i]=conjugate(a.Ap.data[i]);
        if (!f.geostrophic && !f.meanDensityAnomaly) a.A0.data[i]={};
        if (f.meanDensityAnomaly) a.A0.data[i].imag=0;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::nonlinearFluxAndFields(
    const WVState& state,WVFlux& flux,WVRealFieldBundleView fields,
    WVRealFieldBundleView* spatialTendency) {
    if (!tiledNonlinearPrepared_) return unsupported();
    auto status=validateFluxOutput(state,flux); if (!status) return status;
    const auto shape=spatialShape();
    if (fields.shape.first!=shape.first || fields.shape.second!=shape.second ||
        fields.shape.third!=shape.third || fields.shape.fourth!=4)
        return {WVKernelStatusCode::invalidShape,"Tiled nonlinear fields require four full volumes."};
    if (!addressFits(fields.data,4*R_*sizeof(double),alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid tiled physical field storage."};
    for (const auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0}) {
        status=disjoint(fields.data,4*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
    }
    // Registered immutable views are protected even when they differ from the
    // primary state passed to this call.
    for (std::size_t v=0;v<preparedStateViewCount_;++v)
        for (const auto input:{preparedStateViews_[v].coefficients.Ap,preparedStateViews_[v].coefficients.Am,preparedStateViews_[v].coefficients.A0}) {
            status=disjoint(fields.data,4*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
        }
    for (const auto output:{flux.Fp,flux.Fm,flux.F0}) {
        status=disjoint(fields.data,4*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!status) return status;
    }
    if (spatialTendency) {
        if (spatialTendency->shape.first!=shape.first || spatialTendency->shape.second!=shape.second ||
            spatialTendency->shape.third!=shape.third || spatialTendency->shape.fourth!=3)
            return {WVKernelStatusCode::invalidShape,"Tiled nonlinear tendencies require three full volumes."};
        if (!addressFits(spatialTendency->data,3*R_*sizeof(double),alignof(double)))
            return {WVKernelStatusCode::invalidPointer,"Invalid tiled physical tendency storage."};
        status=disjoint(spatialTendency->data,3*R_*sizeof(double),fields.data,4*R_*sizeof(double)); if (!status) return status;
        for (const auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0}) {
            status=disjoint(spatialTendency->data,3*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
        }
        for (std::size_t v=0;v<preparedStateViewCount_;++v)
            for (const auto input:{preparedStateViews_[v].coefficients.Ap,preparedStateViews_[v].coefficients.Am,preparedStateViews_[v].coefficients.A0}) {
                status=disjoint(spatialTendency->data,3*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
            }
        for (const auto output:{flux.Fp,flux.Fm,flux.F0}) {
            status=disjoint(spatialTendency->data,3*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!status) return status;
        }
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    status=preparePhaseForCall(state); if (!status) return status;
    WVRetainedAdvectionWork work;
    work.targets=3;work.fields={fields.data,4*R_*sizeof(double)};
    if (spatialTendency) work.tendencies={spatialTendency->data,3*R_*sizeof(double)};
    work.densityCorrection={geometry().dLnN2.data(),geometry().Nz*sizeof(double)};
    const WVHydrostaticField names[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::w,WVHydrostaticField::eta};
    for (std::size_t f=0;f<4;++f) {
        status=reconstruct(state.coefficients,names[f],WVHydrostaticDerivative::value,WVHydrostaticComponent::all,
            nullptr,true,5,nullptr,&work.base[f]); if (!status) return status;
    }
    const std::size_t channels[]={0,1,3,2},slots[]={1,2,0};
    for (std::size_t t=0;t<work.targets;++t) {
        work.targetSpectra[t]=gridView(slots[t]);
        WVComplexInput prepared;
        status=reconstruct(state.coefficients,names[channels[t]],WVHydrostaticDerivative::z,WVHydrostaticComponent::all,
            nullptr,true,5,nullptr,&prepared,&work.targetSpectra[t]); if (!status) return status;
    }
    WVRetainedAdvectionCounts counts;
    status=horizontal_->advection(*horizontalWorkspace_,work,counts); if (!status) return status;
    ++metrics_.tiledNonlinearCount;
    metrics_.tiledColumnInverseCount+=counts.columnInverses;
    metrics_.tiledRowInverseCount+=counts.rowInverses;
    metrics_.tiledReusedColumnCount+=counts.reusedColumns;
    // Count actual completed physical producers, including streamed derivatives.
    for (std::size_t f=0;f<4;++f) {
        ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(names[f])];
        ++metrics_.reconstructionCount[static_cast<std::size_t>(names[f])][0][0];
    }
    for (std::size_t t=0;t<work.targets;++t) for (std::size_t d=1;d<4;++d) {
        ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(names[channels[t]])];
        ++metrics_.reconstructionCount[static_cast<std::size_t>(names[channels[t]])][d][0];
        ++metrics_.derivativeAdvectionConsumerCount;
    }
    // x and y requests consumed the same prepared base spectrum.
    metrics_.horizontalSpectrumReuseCount+=2*work.targets;
    WVMutableCoefficients target{flux.Fp,flux.Fm,flux.F0};
    status=vertical(1,work.targetSpectra[0].input(),modalView(1)); if (!status) return status;
    status=vertical(1,work.targetSpectra[1].input(),modalView(2)); if (!status) return status;
    status=vertical(3,work.targetSpectra[2].input(),modalView()); if (!status) return status;
    return projectedFieldsToCoefficients(modalView(1).input(),modalView(2).input(),modalView().input(),target);
}
WVKernelStatus WVTransformHydrostaticKernel::nonlinearFlux(const WVState& a,WVFlux& b,
    WVRealFieldBundleView* spatialTendency,const WVRealFieldBundleConstView* preparedFields,
    bool projectFlux,WVStateDerivativeAccess* derivativeAccess) {
    if (!projectFlux && !spatialTendency)
        return {WVKernelStatusCode::invalidConfiguration,"Spatial-only nonlinear evaluation requires output storage."};
    auto s=validateFluxOutput(a,b); if (!s) return s; WVMutableCoefficients target{b.Fp,b.Fm,b.F0};
    const auto validateBundle = [&](const double* data,WVShape4D shape,std::size_t count) {
        if (shape.first!=geometry().Nx || shape.second!=geometry().Ny ||
            shape.third!=geometry().Nz || shape.fourth!=count)
            return WVKernelStatus{WVKernelStatusCode::invalidShape,"Forcing field bundle shape mismatch."};
        for (std::size_t channel=0;channel<count;++channel) {
            auto status=volume({data ? data+channel*R_ : nullptr,spatialShape()});
            if (!status) return status;
        }
        return WVKernelStatus::ok();
    };
    if (preparedFields) {
        s=validateBundle(preparedFields->data,preparedFields->shape,4); if (!s) return s;
        for (auto output:{b.Fp,b.Fm,b.F0}) {
            s=disjoint(preparedFields->data,4*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
    }
    if (spatialTendency) {
        s=validateBundle(spatialTendency->data,spatialTendency->shape,3); if (!s) return s;
        for (auto input:{a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) {
            s=disjoint(spatialTendency->data,3*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
        for (auto output:{b.Fp,b.Fm,b.F0}) {
            s=disjoint(spatialTendency->data,3*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
        if (preparedFields) {
            s=disjoint(spatialTendency->data,3*R_*sizeof(double),preparedFields->data,4*R_*sizeof(double)); if (!s) return s;
        }
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    const WVHydrostaticField fields[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::w,WVHydrostaticField::eta};
    const auto derivativeFor=[&](WVHydrostaticField field,WVHydrostaticDerivative derivative,
        double* scratch,const double*& values,const double* productOutput,
        const WVRealOutputConsumer* consumer = nullptr,bool* consumed = nullptr) {
        WVRealVolumeConstView cached{};
        if (derivativeAccess && derivativeAccess->lookup) {
            auto status=derivativeAccess->lookup(derivativeAccess->context,
                static_cast<std::size_t>(field),static_cast<std::size_t>(derivative),cached);
            if (!status) return status;
        }
        if (cached.data) {
            auto status=volume(cached); if (!status) return status;
            status=disjoint(cached.data,R_*sizeof(double),productOutput,R_*sizeof(double));
            if (!status) return status;
            values=cached.data;
            return WVKernelStatus::ok();
        }
        auto status=reconstruct(a.coefficients,field,derivative,
            WVHydrostaticComponent::all,scratch,true,5,consumer); if (!status) return status;
        if (consumer) {
            ++metrics_.derivativeAdvectionConsumerCount;
            if (consumed) *consumed=true;
        }
        if (derivativeAccess && derivativeAccess->capture) {
            status=derivativeAccess->capture(derivativeAccess->context,
                static_cast<std::size_t>(field),static_cast<std::size_t>(derivative),
                {scratch,spatialShape()});
            if (!status) return status;
        }
        values=scratch;
        return WVKernelStatus::ok();
    };
    if (executionOptions_.streamedNonlinear) {
        const double* advectionFields=preparedFields ? preparedFields->data : real_.data();
        if (!preparedFields) for (std::size_t i=0;i<4;++i) {
            s=reconstruct(a.coefficients,fields[i],WVHydrostaticDerivative::value,WVHydrostaticComponent::all,real_.data()+i*R_); if (!s) return s;
        }
        auto* flux=real_.data()+4*R_; auto* derivative=real_.data()+5*R_;
        for (std::size_t targetIndex=0;targetIndex<3;++targetIndex) {
            const auto field=fields[targetIndex==2 ? 3 : targetIndex];
            std::fill_n(flux,R_,0);
            // These direct modal first derivatives use separate dependency keys
            // from grid Laplacians, which preserve sequential horizontal
            // calculus or the order-two vertical operator from retained values.
            for (std::size_t axis=0;axis<3;++axis) {
                kernel_detail::WVAdvectionConsumer advection{flux,advectionFields+axis*R_,
                    advectionFields+3*R_,geometry().dLnN2.data(),R_/geometry().Nz,targetIndex==2 && axis==2};
                const WVRealOutputConsumer consumer{&advection,kernel_detail::WVAdvectionConsumer::consume};
                bool consumed=false;
                const double* derivativeValues=nullptr;
                s=derivativeFor(field,static_cast<WVHydrostaticDerivative>(axis+1),
                    derivative,derivativeValues,flux,
                    executionOptions_.fusedDerivativeAdvection ? &consumer : nullptr,&consumed); if (!s) return s;
                if (!consumed) pointwise_->execute(R_,[&](std::size_t begin,std::size_t end) {
                    kernel_detail::WVAdvectionConsumer::consume(&advection,begin,end,derivativeValues);
                });
            }
            if (spatialTendency) std::copy_n(flux,R_,spatialTendency->data+targetIndex*R_);
            if (projectFlux) {
                const std::size_t projectedIndex=targetIndex==2 ? 0 : targetIndex+1;
                s=project(flux,modalView(projectedIndex),targetIndex==2 ? WVHydrostaticFamily::G : WVHydrostaticFamily::F); if (!s) return s;
            }
        }
        return projectFlux ? projectedFieldsToCoefficients(modalView(1).input(),modalView(2).input(),modalView().input(),target) : WVKernelStatus::ok();
    }
    if (preparedFields) std::copy_n(preparedFields->data,4*R_,real_.data());
    else for (std::size_t i=0;i<4;++i) { s=reconstruct(a.coefficients,fields[i],WVHydrostaticDerivative::value,WVHydrostaticComponent::all,real_.data()+i*R_); if (!s) return s; }
    for (std::size_t targetIndex=0;targetIndex<3;++targetIndex) {
        const auto field=fields[targetIndex==2 ? 3 : targetIndex]; auto* flux=real_.data()+(4+targetIndex)*R_;
        std::fill_n(flux,R_,0);
        for (std::size_t axis=0;axis<3;++axis) {
            const double* derivativeValues=nullptr;
            s=derivativeFor(field,static_cast<WVHydrostaticDerivative>(axis+1),
                real_.data()+7*R_,derivativeValues,flux); if (!s) return s;
            pointwise_->execute(R_,[&](std::size_t begin,std::size_t end) {
                for (std::size_t i=begin;i<end;++i) {
                    const double correction=targetIndex==2 && axis==2 ? real_[3*R_+i]*geometry().dLnN2[i/(R_/geometry().Nz)] : 0;
                    flux[i]-=real_[axis*R_+i]*(derivativeValues[i]+correction);
                }
            });
        }
    }
    if (spatialTendency) std::copy_n(real_.data()+4*R_,3*R_,spatialTendency->data);
    if (!projectFlux) return WVKernelStatus::ok();
    return projectFields(real_.data()+4*R_,real_.data()+5*R_,real_.data()+6*R_,target);
}
WVKernelStatus WVTransformHydrostaticKernel::totalEnergy(const WVCoefficients& a,double& value,WVHydrostaticComponent component) const {
    if (!valid(component)) return unsupported();
    auto s=coefficients(a); if (!s) return s;
    const auto norm=[](WVComplex64 x) { return x.real*x.real+x.imag*x.imag; }; double sum=0;
    for (std::size_t i=0;i<S_;++i) {
        if (selected(factors_[i],component,true)) sum+=factors_[i].waveEnergy*(norm(a.Ap.data[i])+norm(a.Am.data[i]));
        if (selected(factors_[i],component,false)) sum+=factors_[i].balancedEnergy*norm(a.A0.data[i]);
    }
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Hydrostatic energy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::reduceHorizontalSpeedMaximum(
    WVRealVolumeConstView u,WVRealVolumeConstView v,double& value) {
    auto status=volume(u); if (!status) return status;
    status=volume(v); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    double maximum=0;
    std::mutex merge;
    pointwise_->execute(R_,[&](std::size_t begin,std::size_t end) {
        double local=0;
        for (std::size_t i=begin;i<end;++i)
            local=std::max(local,std::hypot(u.data[i],v.data[i]));
        // Each chunk retains the serial hypot/max semantics, including robust
        // scaling and ignored NaN contributions. Only one merge per chunk.
        std::lock_guard<std::mutex> lock(merge);
        maximum=std::max(maximum,local);
    });
    value=maximum;
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::totalEnstrophy(const WVCoefficients& a,double& value) const {
    auto s=coefficients(a); if (!s) return s; double sum=0;
    for (std::size_t i=0;i<S_;++i) sum+=factors_[i].enstrophy*(a.A0.data[i].real*a.A0.data[i].real+a.A0.data[i].imag*a.A0.data[i].imag);
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Hydrostatic enstrophy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::totalEnergySpatiallyIntegrated(const WVState& a,double& value,WVHydrostaticComponent component) {
    if (!valid(component)) return unsupported();
    auto s=validateStateForCall(a); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    const WVHydrostaticField fields[]={WVHydrostaticField::u,WVHydrostaticField::v,WVHydrostaticField::eta};
    for (std::size_t i=0;i<3;++i) { s=reconstruct(a.coefficients,fields[i],WVHydrostaticDerivative::value,component,real_.data()+i*R_); if (!s) return s; }
    double sum=0; const auto plane=R_/geometry().Nz;
    for (std::size_t i=0;i<R_;++i) sum+=geometry().z_int[i/plane]*(real_[i]*real_[i]+real_[R_+i]*real_[R_+i]+geometry().N2[i/plane]*real_[2*R_+i]*real_[2*R_+i])/(2*plane);
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Hydrostatic spatial energy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformHydrostaticKernel::verticalCalculus(const double* values,WVHydrostaticFamily family,
    unsigned order,bool integral,double* result,std::size_t columns,bool verticalFirst) {
    const auto& g=geometry();
    const auto index=[&](std::size_t column,std::size_t z) {
        return verticalFirst ? z+g.Nz*column : column+columns*z;
    };
    return kernel_detail::applyStratifiedVerticalCalculus(g,columns,family==WVHydrostaticFamily::F,
        order,integral,[&](std::size_t column,std::size_t z) { return values[index(column,z)]; },
        [&](std::size_t column,std::size_t z,double value) { result[index(column,z)]=value; },
        gridView(),gridView(1),modalView(),
        [&](std::size_t operation,WVComplexInput input,WVComplexOutput output) { return vertical(operation,input,output); });
}
WVKernelStatus WVTransformHydrostaticKernel::applyVerticalCalculus(WVRealConstView a,bool inputIsF,
    unsigned order,bool integral,WVRealView b) {
    const auto& g=geometry();
    if (order<1 || order>4 || (integral && order!=1)) return unsupported();
    if (a.shape.rows!=g.Nz || !a.shape.columns || b.shape.rows!=g.Nz || b.shape.columns!=a.shape.columns)
        return {WVKernelStatusCode::invalidShape,"Vertical calculus requires matching nonempty [Nz,Ncolumns] arrays."};
    if (a.shape.columns>std::numeric_limits<std::size_t>::max()/g.Nz ||
        a.shape.columns*g.Nz>std::numeric_limits<std::size_t>::max()/sizeof(double))
        return {WVKernelStatusCode::sizeOverflow,"Vertical calculus extent overflow."};
    const auto bytes=a.shape.columns*g.Nz*sizeof(double);
    if (!addressFits(a.data,bytes,alignof(double)) || !addressFits(b.data,bytes,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid vertical calculus storage."};
    auto status=mutableOutputOutsidePreparedState(b.data,bytes); if (!status) return status;
    status=disjoint(a.data,bytes,b.data,bytes); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return verticalCalculus(a.data,inputIsF ? WVHydrostaticFamily::F : WVHydrostaticFamily::G,
        order,integral,b.data,a.shape.columns,true);
}
WVKernelStatus WVTransformHydrostaticKernel::differentiateVertical(WVRealVolumeConstView a,WVHydrostaticFamily family,unsigned order,WVRealVolumeView b) {
    if (!valid(family) || order<1 || order>4) return unsupported();
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return verticalCalculus(a.data,family,order,false,b.data,R_/geometry().Nz,false);
}
WVKernelStatus WVTransformHydrostaticKernel::integrateVertical(WVRealVolumeConstView a,WVHydrostaticFamily family,WVRealVolumeView b) {
    if (!valid(family)) return unsupported();
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return verticalCalculus(a.data,family,1,true,b.data,R_/geometry().Nz,false);
}
WVKernelStatus WVTransformHydrostaticKernel::differentiateHorizontal(WVRealVolumeConstView a,bool xDerivative,WVRealVolumeView b) {
    return differentiateHorizontal(a,xDerivative,1,b);
}
WVKernelStatus WVTransformHydrostaticKernel::differentiateHorizontal(WVRealVolumeConstView a,
    bool xDerivative,unsigned order,WVRealVolumeView b) {
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=mutableOutputOutsidePreparedState(b.data,R_*sizeof(double)); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{b.data,R_*sizeof(double)},xDerivative,order);
}
WVKernelStatus WVTransformHydrostaticKernel::advectScalarWithAdvectionFields(WVRealVolumeConstView a,WVRealFieldBundleConstView fields,bool antialias,WVRealVolumeView b,bool xyOnly) {
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    const auto& g=geometry();
    if (fields.shape.first!=g.Nx || fields.shape.second!=g.Ny || fields.shape.third!=g.Nz || fields.shape.fourth!=3)
        return {WVKernelStatusCode::invalidShape,"Hydrostatic advection requires [Nx,Ny,Nz,3] velocity fields."};
    if (!addressFits(fields.data,3*R_*sizeof(double),alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid advection field storage."};
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    s=disjoint(fields.data,3*R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{real_.data(),R_*sizeof(double)},true); if (!s) return s;
    s=horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{real_.data()+R_,R_*sizeof(double)},false); if (!s) return s;
    if (!xyOnly) { s=verticalCalculus(a.data,WVHydrostaticFamily::F,1,false,real_.data()+2*R_,R_/geometry().Nz,false); if (!s) return s; }
    for (std::size_t i=0;i<R_;++i) b.data[i]=-fields.data[i]*real_[i]-fields.data[R_+i]*real_[R_+i]-(xyOnly ? 0.0 : fields.data[2*R_+i]*real_[2*R_+i]);
    if (antialias) {
        s=horizontal_->forward(*horizontalWorkspace_,{b.data,R_*sizeof(double)},gridView()); if (!s) return s;
        s=horizontal_->inverse(*horizontalWorkspace_,gridView().input(),{b.data,R_*sizeof(double)});
    }
    return s;
}
} // namespace wavevortex
