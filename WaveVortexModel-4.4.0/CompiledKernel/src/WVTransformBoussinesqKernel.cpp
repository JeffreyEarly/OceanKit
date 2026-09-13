#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WVPreparedFieldCache.hpp"
#include "WVAdvectionConsumer.hpp"
#include "WVSpectralValidation.hpp"
#include "WVPreparedModeExecutor.hpp"
#include "WVVariableComplexBuffer.hpp"
#include "WVStratifiedVerticalCalculus.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
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
WVCoefficients view(WVMutableCoefficients a) { return {{a.Ap.data,a.Ap.shape},{a.Am.data,a.Am.shape},{a.A0.data,a.A0.shape}}; }
bool valid(WVBoussinesqFamily f) { return f>=WVBoussinesqFamily::F && f<=WVBoussinesqFamily::Gw; }
bool verticalOperation(WVStratifiedModalOperator operation,std::size_t& index,
    bool& inputModal,bool& outputModal) {
    switch(operation) {
        case WVStratifiedModalOperator::reconstructF: index=0; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectF: index=1; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::reconstructG: index=2; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectG: index=3; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::reconstructFw: index=4; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectFw: index=5; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::reconstructGw: index=6; inputModal=true; outputModal=false; return true;
        case WVStratifiedModalOperator::projectGw: index=7; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::balancedGToWaveG: index=8; inputModal=true; outputModal=true; return true;
        case WVStratifiedModalOperator::projectWaveDivergence: index=9; inputModal=false; outputModal=true; return true;
        case WVStratifiedModalOperator::projectWaveVerticalVelocity: index=10; inputModal=false; outputModal=true; return true;
        default: return false;
    }
}
bool valid(WVBoussinesqComponent c) { return c>=WVBoussinesqComponent::all && c<=WVBoussinesqComponent::meanDensityAnomaly; }
bool surface(WVBoussinesqField f) { return f>=WVBoussinesqField::ssh && f<=WVBoussinesqField::ssv; }
bool selected(const WVBoussinesqModeFactors& f,WVBoussinesqComponent c,bool wave) {
    if (c==WVBoussinesqComponent::all) return wave ? f.wave || f.inertial : f.geostrophic || f.meanDensityAnomaly;
    if (wave) return c==WVBoussinesqComponent::wave ? f.wave : c==WVBoussinesqComponent::inertial && f.inertial;
    return c==WVBoussinesqComponent::geostrophic ? f.geostrophic : c==WVBoussinesqComponent::meanDensityAnomaly && f.meanDensityAnomaly;
}
WVKernelStatus unsupported() { return {WVKernelStatusCode::unsupportedOperation,"Unsupported boussinesq field, derivative, family or component."}; }
WVKernelStatus reentrant() { return {WVKernelStatusCode::reentrantExecution,"Boussinesq workspace is already active."}; }
}

WVTransformBoussinesqKernel::~WVTransformBoussinesqKernel() = default;

WVKernelStatus WVTransformBoussinesqKernel::create(std::shared_ptr<const WVStratifiedModalSource> source,
    std::unique_ptr<WVFFTEngine> engine,std::unique_ptr<WVTransformBoussinesqKernel>& result,MatrixBackendFactory factory,WVVariableExecutionOptions options) {
    try {
        if (!source || !engine || !factory) return {WVKernelStatusCode::invalidConfiguration,"Scientific source, FFT engine and matrix backend factory are required."};
        if (options.spectralSchedule!=WVVariableSpectralSchedule::establishedInterleaved &&
            options.spectralSchedule!=WVVariableSpectralSchedule::compactSplitFusedViews)
            return {WVKernelStatusCode::invalidConfiguration,"Unknown variable spectral schedule."};
        if (!options.pointwiseWorkers || !options.verticalGroupWorkers)
            return {WVKernelStatusCode::invalidConfiguration,"Pointwise and vertical group worker counts must be positive."};
        if (options.usesCompactSplitViews() &&
            (options.horizontalSchedule!=WVRetainedHorizontalSchedule::streamingPrunedTile16 || !options.streamedNonlinear))
            return {WVKernelStatusCode::invalidConfiguration,"Compact split views require the streaming pruned nonlinear schedule."};
        const auto representation=options.usesCompactSplitViews() ? WVComplexRepresentation::split : WVComplexRepresentation::interleaved;
        const auto& g=source->geometry();
        if (g.transformClass!="WVTransformBoussinesq" || g.Nx<2 || g.Ny<2 || g.Nz<3 || !g.Nj || g.Nj>=g.Nz || !g.Nkl ||
            g.j.size()!=g.Nj || g.h_0.size()!=g.Nj || g.k.size()!=g.Nkl || g.l.size()!=g.Nkl || g.modes.size()!=g.Nkl ||
            g.h_pm.size()!=g.Nj*g.Nkl || g.waveGroup.size()!=g.Nkl || g.K2unique.empty() || g.z.size()!=g.Nz || g.N2.size()!=g.Nz || g.rho_nm0.size()!=g.Nz || g.dLnN2.size()!=g.Nz || g.z_int.size()!=g.Nz ||
            source->sourceIdentity().empty() || source->modeSetIdentity().empty())
            return {WVKernelStatusCode::invalidConfiguration,"Invalid boussinesq scientific source."};
        auto candidate=std::unique_ptr<WVTransformBoussinesqKernel>(new WVTransformBoussinesqKernel);
        auto& c=*candidate; c.source_=std::move(source);
        c.S_=product(g.Nj,g.Nkl); c.H_=product(g.Nz,g.Nkl); c.R_=product(product(g.Nx,g.Ny),g.Nz);
        product(product(6,c.S_),sizeof(WVComplex64)); product(product(5,c.H_),sizeof(WVComplex64));
        product(product(options.streamedNonlinear ? 6 : 11,c.R_),sizeof(double)); product(c.S_,sizeof(WVBoussinesqModeFactors));
        c.engineIdentifier_=engine->identifier(); c.engineLibraryIdentity_=engine->libraryIdentity();
        c.executionOptions_=options;
        WVRetainedHorizontalSpecification horizontal;
        auto status=c.source_->horizontalSpecification(g.Nz,representation,"boussinesq-grid",horizontal); if (!status) return status;
        horizontal.schedule=options.horizontalSchedule; horizontal.outerWorkers=options.horizontalWorkers;
        status=WVRetainedHorizontalOperator::create(horizontal,std::move(engine),c.horizontal_); if (!status) return status;
        status=c.horizontal_->createWorkspace(c.horizontalWorkspace_); if (!status) return status;
        if (options.tiledNonlinear && options.sharedFieldGradients && options.streamedNonlinear) {
            status=c.horizontal_->prepareAdvection(*c.horizontalWorkspace_,4);
            if (!status && status.code!=WVKernelStatusCode::unsupportedOperation) return status;
            c.tiledNonlinearPrepared_=static_cast<bool>(status);
        }
        const WVStratifiedModalOperator operations[]={WVStratifiedModalOperator::reconstructF,WVStratifiedModalOperator::projectF,WVStratifiedModalOperator::reconstructG,WVStratifiedModalOperator::projectG,
            WVStratifiedModalOperator::reconstructFw,WVStratifiedModalOperator::projectFw,WVStratifiedModalOperator::reconstructGw,WVStratifiedModalOperator::projectGw,
            WVStratifiedModalOperator::balancedGToWaveG,WVStratifiedModalOperator::projectWaveDivergence,WVStratifiedModalOperator::projectWaveVerticalVelocity};
        const char* inputs[]={"F-modal","F-grid","G-modal","G-grid","Fw-modal","F-grid","Gw-modal","G-grid","G-modal","F-grid","G-grid"};
        const char* outputs[]={"F-grid","F-modal","G-grid","G-modal","F-grid","Fw-modal","G-grid","Gw-modal","Gw-modal","Gw-modal","Gw-modal"};
        for (std::size_t i=0;i<c.vertical_.size();++i) {
            const auto rows=(i<8 && i%2==0) || i==8 ? g.Nj : g.Nz;
            const auto outRows=i<8 && i%2==0 ? g.Nz : g.Nj;
            WVComplexLayout in{rows,g.Nkl,1,rows,representation,inputs[i],c.source_->modeSetIdentity()};
            WVComplexLayout out{outRows,g.Nkl,1,outRows,representation,outputs[i],c.source_->modeSetIdentity()};
            std::unique_ptr<WVVerticalMatrixBackend> backend; status=factory(backend); if (!status) return status;
            status=c.source_->prepareVertical(operations[i],in,out,std::move(backend),c.vertical_[i]); if (!status) return status;
            if (options.verticalGroupWorkers>1 && !c.vertical_[i]->supportsConcurrentCalls())
                return {WVKernelStatusCode::unsupportedOperation,"Vertical group workers require a concurrent matrix backend."};
            status=c.vertical_[i]->createWorkspace(options.verticalGroupWorkers,c.verticalWorkspace_[i]); if (!status) return status;
        }
        if (options.verticalGroupWorkers>1) {
            status=WVVerticalGroupExecutor::create(options.verticalGroupWorkers,c.verticalGroups_); if (!status) return status;
        }
        const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
        if (!std::isfinite(f) || f==0 || !std::isfinite(g.g) || !(g.g>0) || !std::isfinite(g.Lz) || !(g.Lz>0) || !std::isfinite(g.rho0) || !(g.rho0>0))
            return {WVKernelStatusCode::invalidConfiguration,"Invalid boussinesq physical constants."};
        c.factors_.resize(c.S_);
        for (std::size_t i=0;i<c.S_;++i) {
            const auto j=i%g.Nj,mode=i/g.Nj; auto& a=c.factors_[i];
            const double k=g.k[mode],l=g.l[mode],K2=k*k+l*l,h=g.h_0[j],hw=g.h_pm[i];
            if (!std::isfinite(K2) || !std::isfinite(h) || !(h>0) || !std::isfinite(hw) || !(hw>0) || !std::isfinite(g.j[j]) || g.j[j]<0)
                return {WVKernelStatusCode::invalidConfiguration,"Invalid boussinesq mode or equivalent depth."};
            const bool barotropic=g.j[j]==0;
            a.omega=std::sqrt(g.g*hw*K2+f*f);
            if (!std::isfinite(a.omega) || !(a.omega>0)) return {WVKernelStatusCode::numericalFailure,"Boussinesq frequency overflow."};
            a.wave=K2>0 && !barotropic; a.inertial=K2==0;
            a.geostrophic=K2>0; a.meanDensityAnomaly=K2==0 && !barotropic;
            if (a.wave) {
                const double K=std::sqrt(K2),alpha=std::atan2(l,k);
                a.UAp={std::cos(alpha),-f/a.omega*std::sin(alpha)};
                a.VAp={std::sin(alpha),f/a.omega*std::cos(alpha)};
                a.WAp={0,-K*hw}; a.NAp=-K*hw/a.omega;
                a.ApmD={0,-1/(2*K*hw)}; a.ApmN=-a.omega/(2*K*hw); a.waveEnergy=2*hw;
            }
            if (a.inertial) { a.UAp={1,0}; a.VAp={0,1}; a.waveEnergy=barotropic ? g.Lz : hw; }
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
                if (!std::isfinite(x)) return {WVKernelStatusCode::numericalFailure,"Boussinesq coefficient factor overflow."};
        }
        c.inertialMode_=g.Nkl;
        bool exactInertialIdentity=true;
        for (std::size_t mode=0;mode<g.Nkl;++mode) {
            const bool exactZero=g.modes[mode].k==0 && g.modes[mode].l==0;
            const bool inertial=c.factors_[g.Nj*mode].inertial;
            if (exactZero!=inertial) exactInertialIdentity=false;
            if (exactZero) {
                if (c.inertialMode_!=g.Nkl) exactInertialIdentity=false;
                c.inertialMode_=mode;
            }
        }
        if (!exactInertialIdentity) c.inertialMode_=g.Nkl;
        const auto modalElements=product(6,c.S_),gridElements=product(5,c.H_);
        const auto maximumElements=static_cast<std::size_t>(PTRDIFF_MAX);
        if (modalElements>maximumElements || gridElements>maximumElements-modalElements) throw std::overflow_error("Boussinesq spectral scratch overflow.");
        c.spectralStorage_=std::make_unique<WVVariableComplexBuffer>(modalElements+gridElements,representation);
        if (options.sharedFieldGradients)
            c.fieldCache_=std::make_unique<kernel_detail::WVPreparedFieldCache>(2*c.S_,c.H_,representation);
        c.pointwise_=std::make_unique<kernel_detail::WVPreparedModeExecutor>(std::min(options.pointwiseWorkers,c.R_));
        c.phase_.resize(c.S_); c.real_.resize((options.streamedNonlinear ? 6 : 11)*c.R_);
        auto& s=c.storage_; s.sharedScientificBytes=c.source_->persistentBytes(); s.preparedBytes=c.horizontal_->persistentBytes();
        s.workspaceBytes=c.horizontalWorkspace_->persistentBytes()+sizeof(WVVariableComplexBuffer)+c.pointwise_->persistentBytes()+
            (c.verticalGroups_ ? c.verticalGroups_->persistentBytes() : 0); s.providerBytesLowerBound=c.horizontal_->providerBytesLowerBound(); s.planBytesLowerBound=c.horizontalWorkspace_->planBytesLowerBound();
        for (std::size_t i=0;i<c.vertical_.size();++i) { s.preparedBytes+=c.vertical_[i]->persistentBytes(); s.workspaceBytes+=c.verticalWorkspace_[i]->persistentBytes(); }
        c.baseSpectralScratchBytes_=c.spectralStorage_->capacityBytes()+c.phase_.capacity()*sizeof(WVComplex64);
        s.spectralScratchBytes=c.baseSpectralScratchBytes_;
        s.realScratchBytes=c.real_.capacity()*sizeof(double); s.factorBytes=c.factors_.capacity()*sizeof(WVBoussinesqModeFactors);
        result=std::move(candidate); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Boussinesq setup allocation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
      catch (const std::system_error& e) { return {WVKernelStatusCode::allocationFailure,e.what()}; }
}
const WVBoussinesqStorage& WVTransformBoussinesqKernel::storage() const noexcept {
    storage_.spectralScratchBytes=baseSpectralScratchBytes_+
        (fieldCache_ ? fieldCache_->capacityBytes() : 0);
    return storage_;
}
std::size_t WVTransformBoussinesqKernel::persistentBytes() const noexcept {
    const auto& s=storage(); return sizeof(*this)+s.sharedScientificBytes+s.preparedBytes+s.workspaceBytes+s.spectralScratchBytes+s.realScratchBytes+s.factorBytes;
}
WVKernelStatus WVTransformBoussinesqKernel::spectral(WVComplexConstView a) const {
    if (a.shape.rows!=geometry().Nj || a.shape.columns!=geometry().Nkl) return {WVKernelStatusCode::invalidShape,"Expected canonical [Nj,Nkl] coefficients."};
    if (!addressFits(a.data,S_*sizeof(WVComplex64),alignof(WVComplex64))) return {WVKernelStatusCode::invalidPointer,"Invalid coefficient storage."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::volume(WVRealVolumeConstView a,bool isSurface) const {
    const auto& g=geometry(); const auto nz=isSurface ? 1 : g.Nz;
    if (a.shape.first!=g.Nx || a.shape.second!=g.Ny || a.shape.third!=nz) return {WVKernelStatusCode::invalidShape,"Unexpected boussinesq field shape."};
    if (!addressFits(a.data,(R_/g.Nz)*nz*sizeof(double),alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid field storage."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::disjoint(const void* a,std::size_t an,const void* b,std::size_t bn) const {
    if (overlap(a,an,b,bn)) return {WVKernelStatusCode::overlappingArrays,"Boussinesq arrays must not overlap."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::outputs(WVMutableCoefficients a) const {
    for (auto x : {a.Ap,a.Am,a.A0}) { auto s=spectral({x.data,x.shape}); if (!s) return s; }
    auto s=disjoint(a.Ap.data,S_*sizeof(WVComplex64),a.Am.data,S_*sizeof(WVComplex64)); if (!s) return s;
    s=disjoint(a.Ap.data,S_*sizeof(WVComplex64),a.A0.data,S_*sizeof(WVComplex64)); if (!s) return s;
    return disjoint(a.Am.data,S_*sizeof(WVComplex64),a.A0.data,S_*sizeof(WVComplex64));
}
WVKernelStatus WVTransformBoussinesqKernel::mutableOutputOutsidePreparedState(WVMutableCoefficients a) const {
    for (const auto output:{a.Ap,a.Am,a.A0}) {
        auto status=mutableOutputOutsidePreparedState(output.data,S_*sizeof(WVComplex64));
        if (!status) return status;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::mutableOutputOutsidePreparedState(
    const void* output,std::size_t bytes) const {
    if (!stateEvaluationActive_) return WVKernelStatus::ok();
    for (std::size_t stateIndex=0;stateIndex<preparedStateViewCount_;++stateIndex)
        for (const auto input:{preparedStateViews_[stateIndex].coefficients.Ap,preparedStateViews_[stateIndex].coefficients.Am,preparedStateViews_[stateIndex].coefficients.A0})
            if (overlap(output,bytes,input.data,S_*sizeof(WVComplex64)))
                return {WVKernelStatusCode::overlappingArrays,"Mutable output overlaps an active immutable boussinesq state view."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::coefficients(const WVCoefficients& a) const {
    for (auto x : {a.Ap,a.Am,a.A0}) {
        auto s=spectral(x); if (!s) return s;
        for (std::size_t i=0;i<S_;++i) if (!std::isfinite(x.data[i].real) || !std::isfinite(x.data[i].imag))
            return {WVKernelStatusCode::numericalFailure,"Nonfinite boussinesq coefficient."};
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::stateContents(const WVState& a) const {
    auto s=coefficients(a.coefficients); if (!s) return s;
    if (!std::isfinite(a.t) || !std::isfinite(a.t0) || !std::isfinite(a.t-a.t0)) return {WVKernelStatusCode::invalidConfiguration,"Nonfinite boussinesq time or elapsed time."};
    for (const auto& f:factors_) if (!std::isfinite(f.omega*(a.t-a.t0))) return {WVKernelStatusCode::numericalFailure,"Boussinesq phase overflow."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::state(const WVState& a) {
    ++metrics_.stateValidationCount;
    return stateContents(a);
}
WVKernelStatus WVTransformBoussinesqKernel::preparePhase(
    double t,double t0,const WVState* validatedState) {
    const bool exactValidatedState=validatedState!=nullptr &&
        stateEvaluationActive_ && matchesStateEvaluation(*validatedState) &&
        validatedState->t==t && validatedState->t0==t0;
    if (!exactValidatedState) {
        if (!std::isfinite(t) || !std::isfinite(t0) || !std::isfinite(t-t0)) return {WVKernelStatusCode::invalidConfiguration,"Nonfinite boussinesq time or elapsed time."};
        for (const auto& f:factors_) if (!std::isfinite(f.omega*(t-t0))) return {WVKernelStatusCode::numericalFailure,"Boussinesq phase overflow."};
    }
    ++metrics_.phasePreparationCount;
    for (std::size_t i=0;i<S_;++i) { const double a=factors_[i].omega*(t-t0); phase_[i]={std::cos(a),std::sin(a)}; }
    return WVKernelStatus::ok();
}
bool WVTransformBoussinesqKernel::matchesStateEvaluation(const WVState& a) const noexcept {
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
WVKernelStatus WVTransformBoussinesqKernel::validateStateEvaluation(
    const WVState& a) const noexcept {
    if (!matchesStateEvaluation(a))
        return {WVKernelStatusCode::invalidConfiguration,
            "State does not belong to the active Boussinesq evaluation."};
    return WVKernelStatus::ok();
}
std::size_t WVTransformBoussinesqKernel::stateEvaluationComponent(
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
std::size_t WVTransformBoussinesqKernel::fieldPreparationView(
    const WVCoefficients& a) const noexcept {
    const auto sameView=[](WVComplexConstView x,WVComplexConstView y) {
        return x.data==y.data && x.shape.rows==y.shape.rows && x.shape.columns==y.shape.columns;
    };
    if (!stateEvaluationActive_) return 0;
    for (std::size_t i=0;i<preparedStateViewCount_;++i) if (
        sameView(a.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(a.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(a.A0,preparedStateViews_[i].coefficients.A0))
        return preparedStateViewIds_[i];
    return 0;
}
WVKernelStatus WVTransformBoussinesqKernel::validateStateForCall(const WVState& a) {
    if (!stateEvaluationActive_) return state(a);
    if (!matchesStateEvaluation(a))
        return {WVKernelStatusCode::invalidConfiguration,"State does not match the active Boussinesq evaluation."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::preparePhaseForCall(const WVState& a) {
    if (stateEvaluationActive_) {
        if (!matchesStateEvaluation(a))
            return {WVKernelStatusCode::invalidConfiguration,"State does not match the active Boussinesq evaluation."};
        return WVKernelStatus::ok();
    }
    if (fieldCache_) fieldCache_->clear();
    return preparePhase(a.t,a.t0);
}
WVKernelStatus WVTransformBoussinesqKernel::prepareProjectionPhaseForCall(double t,double t0) {
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0)
            return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active Boussinesq evaluation."};
        return WVKernelStatus::ok();
    }
    return preparePhase(t,t0);
}
WVKernelStatus WVTransformBoussinesqKernel::beginStateEvaluation(const WVState& a) {
    return beginStateEvaluation(a,nullptr);
}
WVKernelStatus WVTransformBoussinesqKernel::beginStateEvaluation(const WVState& a,const void* evaluationOwner) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (stateEvaluationActive_)
        return {WVKernelStatusCode::reentrantExecution,"Boussinesq state evaluation is already active."};
    auto s=state(a); if (!s) return s;
    preparedState_=a;
    preparedStateViews_[0]=a;
    preparedStateComponents_[0]=0;
    preparedStateViewIds_={};
    nextStateViewId_=1;
    preparedStateViewCount_=1;
    preparedStateOwner_=evaluationOwner;
    stateEvaluationActive_=true;
    if (fieldCache_) fieldCache_->clear();
    // state(a) validated this exact registered view, including every phase
    // product, so phase preparation must not repeat the same finite scan.
    s=preparePhase(a.t,a.t0,&a);
    if (!s) {
        preparedState_={};
        preparedStateViews_[0]={};
        preparedStateComponents_[0]=0;
        preparedStateViewIds_[0]=0;
        preparedStateViewCount_=0;
        preparedStateOwner_=nullptr;
        stateEvaluationActive_=false;
        return s;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::addStateEvaluationView(
    const WVState& a,const void* evaluationOwner,std::size_t componentIdentity) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No boussinesq state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Boussinesq state view owner does not match the active evaluation."};
    if (a.t!=preparedState_.t || a.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Additional boussinesq state view must use the active evaluation times."};
    if (componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"Boussinesq component identity is out of range."};
    if (matchesStateEvaluation(a))
        return stateEvaluationComponent(a)==componentIdentity ? WVKernelStatus::ok() :
            WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                "Boussinesq state view is already registered with another component identity."};
    if (preparedStateViewCount_==preparedStateViews_.size())
        return {WVKernelStatusCode::invalidConfiguration,"Boussinesq state evaluation view capacity exceeded."};
    auto s=state(a); if (!s) return s;
    if (nextStateViewId_==std::numeric_limits<std::size_t>::max())
        return {WVKernelStatusCode::sizeOverflow,"Boussinesq state view generation overflow."};
    preparedStateViews_[preparedStateViewCount_]=a;
    preparedStateComponents_[preparedStateViewCount_]=componentIdentity;
    preparedStateViewIds_[preparedStateViewCount_]=nextStateViewId_++;
    ++preparedStateViewCount_;
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::removeStateEvaluationView(
    const WVState& a,const void* evaluationOwner,std::size_t componentIdentity) {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No Boussinesq state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Boussinesq state view owner does not match the active evaluation."};
    if (componentIdentity==0 || componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"The primary Boussinesq state view cannot be removed."};
    if (a.t!=preparedState_.t || a.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Removed Boussinesq state view must use the active evaluation times."};
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
            preparedStateViews_[j-1]=preparedStateViews_[j];
            preparedStateComponents_[j-1]=preparedStateComponents_[j];
            preparedStateViewIds_[j-1]=preparedStateViewIds_[j];
        }
        --preparedStateViewCount_;
        preparedStateViews_[preparedStateViewCount_]={};
        preparedStateComponents_[preparedStateViewCount_]=0;
        preparedStateViewIds_[preparedStateViewCount_]=0;
        return WVKernelStatus::ok();
    }
    return {WVKernelStatusCode::invalidConfiguration,
        "Boussinesq state view is not registered for this component."};
}
WVKernelStatus WVTransformBoussinesqKernel::endStateEvaluation() {
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No Boussinesq state evaluation is active."};
    if (fieldCache_) fieldCache_->clear();
    preparedState_={};
    for (auto& stateView:preparedStateViews_) stateView={};
    preparedStateComponents_={};
    preparedStateViewIds_={};
    preparedStateViewCount_=0;
    preparedStateOwner_=nullptr;
    stateEvaluationActive_=false;
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::validateFluxOutput(const WVState& a,const WVFlux& b) {
    auto s=validateStateForCall(a); if (!s) return s;
    WVMutableCoefficients target{b.Fp,b.Fm,b.F0}; s=outputs(target); if (!s) return s;
    for (auto inputView : {a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0})
        for (auto outputView : {b.Fp,b.Fm,b.F0}) {
            s=disjoint(inputView.data,S_*sizeof(WVComplex64),outputView.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
    return mutableOutputOutsidePreparedState(target);
}
WVKernelStatus WVTransformBoussinesqKernel::preparedPhase(const WVState& a,WVComplexConstView& result) {
    result={};
    auto s=validateStateForCall(a); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=preparePhaseForCall(a); if (!s) return s;
    result={phase_.data(),spectralShape()};
    return WVKernelStatus::ok();
}
WVComplexOutput WVTransformBoussinesqKernel::modalView(std::size_t slot) { return spectralStorage_->output(slot*S_,S_); }
WVComplexOutput WVTransformBoussinesqKernel::gridView(std::size_t slot) { return spectralStorage_->output(6*S_+slot*H_,H_); }
WVKernelStatus WVTransformBoussinesqKernel::vertical(std::size_t operation,WVComplexInput a,WVComplexOutput b) {
    ++metrics_.verticalOperatorExecutionCount;
    auto status=verticalGroups_ ?
        vertical_[operation]->execute(*verticalWorkspace_[operation],*verticalGroups_,a,b) :
        vertical_[operation]->execute(*verticalWorkspace_[operation],a,b);
    if (status) metrics_.verticalMatrixGroupExecutionCount+=vertical_[operation]->preparedGroupCount();
    return status;
}
WVKernelStatus WVTransformBoussinesqKernel::verticalColumn(std::size_t operation,WVComplexInput a,
    WVComplexOutput b,std::size_t retainedColumn) {
    ++metrics_.verticalOperatorExecutionCount;
    auto status=vertical_[operation]->executeColumn(*verticalWorkspace_[operation],a,b,retainedColumn);
    if (status) ++metrics_.verticalMatrixGroupExecutionCount;
    return status;
}
WVKernelStatus WVTransformBoussinesqKernel::project(const double* a,WVComplexOutput b,WVBoussinesqFamily family) {
    auto s=horizontal_->forward(*horizontalWorkspace_,{a,R_*sizeof(double)},gridView()); if (!s) return s;
    return vertical(2*static_cast<std::size_t>(family)+1,gridView().input(),b);
}
WVKernelStatus WVTransformBoussinesqKernel::transformToSpatial(WVComplexConstView a,WVBoussinesqFamily family,WVRealVolumeView b) {
    if (!valid(family)) return unsupported();
    auto s=spectral(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,S_*sizeof(WVComplex64),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    auto values=input(a.data,S_);
    if (executionOptions_.usesCompactSplitViews()) {
        const auto target=modalView(); copy(a.data,target,S_); values=target.input();
    }
    s=vertical(2*static_cast<std::size_t>(family),values,gridView()); if (!s) return s;
    return horizontal_->inverse(*horizontalWorkspace_,gridView().input(),{b.data,R_*sizeof(double)});
}
WVKernelStatus WVTransformBoussinesqKernel::transformFromSpatial(WVRealVolumeConstView a,WVBoussinesqFamily family,WVComplexView b) {
    if (!valid(family)) return unsupported();
    auto s=volume(a); if (!s) return s; s=spectral({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,S_*sizeof(WVComplex64)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews()) return project(a.data,output(b.data,S_),family);
    s=project(a.data,modalView(),family); if (!s) return s;
    copy(modalView().input(),b.data,S_); return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::applyVertical(WVStratifiedModalOperator operation,
    WVComplexConstView a,WVComplexView b) {
    std::size_t index=0; bool inputModal=false,outputModal=false;
    if (!verticalOperation(operation,index,inputModal,outputModal)) return unsupported();
    const auto& g=geometry(); const auto inputRows=inputModal ? g.Nj : g.Nz;
    const auto outputRows=outputModal ? g.Nj : g.Nz;
    if (a.shape.rows!=inputRows || a.shape.columns!=g.Nkl ||
        b.shape.rows!=outputRows || b.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Unexpected Boussinesq vertical-transform shape."};
    const auto inputCount=inputRows*g.Nkl,outputCount=outputRows*g.Nkl;
    if (!addressFits(a.data,inputCount*sizeof(WVComplex64),alignof(WVComplex64)) ||
        !addressFits(b.data,outputCount*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid Boussinesq vertical-transform storage."};
    auto status=mutableOutputOutsidePreparedState(b.data,outputCount*sizeof(WVComplex64)); if (!status) return status;
    status=disjoint(a.data,inputCount*sizeof(WVComplex64),b.data,outputCount*sizeof(WVComplex64)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews()) return vertical(index,input(a.data,inputCount),output(b.data,outputCount));
    const auto in=inputModal ? modalView() : gridView();
    const auto out=outputModal ? modalView(inputModal ? 1 : 0) : gridView(inputModal ? 0 : 1);
    copy(a.data,in,inputCount); status=vertical(index,in.input(),out); if (!status) return status;
    copy(out.input(),b.data,outputCount); return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::applyVerticalColumn(WVStratifiedModalOperator operation,
    std::size_t retainedColumn,WVComplexConstView a,WVComplexView b) {
    std::size_t index=0; bool inputModal=false,outputModal=false;
    if (!verticalOperation(operation,index,inputModal,outputModal)) return unsupported();
    const auto& g=geometry(); const auto inputRows=inputModal ? g.Nj : g.Nz;
    const auto outputRows=outputModal ? g.Nj : g.Nz;
    if (retainedColumn>=g.Nkl || a.shape.rows!=inputRows || a.shape.columns!=1 ||
        b.shape.rows!=outputRows || b.shape.columns!=1)
        return {WVKernelStatusCode::invalidShape,"Unexpected Boussinesq vertical-column shape or retained column."};
    if (!addressFits(a.data,inputRows*sizeof(WVComplex64),alignof(WVComplex64)) ||
        !addressFits(b.data,outputRows*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid Boussinesq vertical-column storage."};
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
WVKernelStatus WVTransformBoussinesqKernel::horizontalForward(WVRealVolumeConstView a,WVComplexView b) {
    const auto& g=geometry(); auto status=volume(a); if (!status) return status;
    if (b.shape.rows!=g.Nz || b.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Expected Boussinesq Fourier [Nz,Nkl] output."};
    if (!addressFits(b.data,H_*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid Boussinesq Fourier output storage."};
    status=mutableOutputOutsidePreparedState(b.data,H_*sizeof(WVComplex64)); if (!status) return status;
    status=disjoint(a.data,R_*sizeof(double),b.data,H_*sizeof(WVComplex64)); if (!status) return status;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    if (!executionOptions_.usesCompactSplitViews())
        return horizontal_->forward(*horizontalWorkspace_,{a.data,R_*sizeof(double)},output(b.data,H_));
    auto grid=gridView(); status=horizontal_->forward(*horizontalWorkspace_,{a.data,R_*sizeof(double)},grid); if (!status) return status;
    copy(grid.input(),b.data,H_); return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::horizontalInverse(WVComplexConstView a,WVRealVolumeView b) {
    const auto& g=geometry(); auto status=volume({b.data,b.shape}); if (!status) return status;
    if (a.shape.rows!=g.Nz || a.shape.columns!=g.Nkl)
        return {WVKernelStatusCode::invalidShape,"Expected Boussinesq Fourier [Nz,Nkl] input."};
    if (!addressFits(a.data,H_*sizeof(WVComplex64),alignof(WVComplex64)))
        return {WVKernelStatusCode::invalidPointer,"Invalid Boussinesq Fourier input storage."};
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
WVKernelStatus WVTransformBoussinesqKernel::projectFields(const double* u,const double* v,const double* w,const double* eta,WVMutableCoefficients b) {
    const auto uh=gridView(),vh=gridView(1),nh=gridView(2),wh=gridView(3),work=gridView(4);
    auto s=horizontal_->forward(*horizontalWorkspace_,{u,R_*sizeof(double)},uh); if (!s) return s;
    s=horizontal_->forward(*horizontalWorkspace_,{v,R_*sizeof(double)},vh); if (!s) return s;
    s=horizontal_->forward(*horizontalWorkspace_,{eta,R_*sizeof(double)},nh); if (!s) return s;
    if (w) { s=horizontal_->forward(*horizontalWorkspace_,{w,R_*sizeof(double)},wh); if (!s) return s; }
    return projectSpectralFields(uh.input(),vh.input(),wh.input(),nh.input(),work,w!=nullptr,b);
}
WVKernelStatus WVTransformBoussinesqKernel::projectSpectralFields(WVComplexInput uh,WVComplexInput vh,
    WVComplexInput wh,WVComplexInput nh,WVComplexOutput work,bool hasW,WVMutableCoefficients b) {
    const auto& g=geometry();
    const auto U=modalView(),V=modalView(1),N=modalView(2),density=modalView(3),divergence=modalView(4),temp=modalView(5);
    auto s=vertical(1,uh,U); if (!s) return s; s=vertical(1,vh,V); if (!s) return s; s=vertical(3,nh,N); if (!s) return s;
    pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
      for (std::size_t i=begin;i<end;++i) {
        const auto& f=factors_[i]; const auto mode=i/g.Nj;
        const auto uValue=read(U.input(),i),vValue=read(V.input(),i),nValue=read(N.input(),i);
        const auto zeta=subtract(multiply(vValue,{0,g.k[mode]}),multiply(uValue,{0,g.l[mode]}));
        b.A0.data[i]=add(scale(zeta,f.A0Z),scale(nValue,f.A0N));
        write(N,i,subtract(nValue,scale(b.A0.data[i],f.NA0)));
        // The has-W projection obtains divergence from the physical spectra below.
        if (!hasW) write(U,i,scale(add(multiply(uValue,{0,g.k[mode]}),multiply(vValue,{0,g.l[mode]})),g.h_0[i%g.Nj]));
    }
    });
    s=vertical(8,N.input(),density); if (!s) return s;
    if (!hasW) {
        s=vertical(8,U.input(),divergence); if (!s) return s;
        for (std::size_t i=0;i<S_;++i) write(divergence,i,multiply(read(divergence.input(),i),factors_[i].ApmD));
    } else {
        pointwise_->execute(g.Nkl,[&](std::size_t begin,std::size_t end) {
          for (std::size_t mode=begin;mode<end;++mode) {
            const double K=std::hypot(g.k[mode],g.l[mode]);
            for (std::size_t z=0;z<g.Nz;++z) { const auto i=z+g.Nz*mode; write(work,i,K==0 ? WVComplex64{} : scale(add(scale(read(uh,i),g.k[mode]),scale(read(vh,i),g.l[mode])),1/(2*K))); }
        }
        });
        s=vertical(9,work.input(),divergence); if (!s) return s;
        s=vertical(10,wh,temp); if (!s) return s;
    }
    // Fio is the zero-wavenumber wave F basis, not the balanced F basis.
    // Only the exact inertial column consumes these two projections.
    if (executionOptions_.inertialOnlyProjection && inertialMode_<g.Nkl) {
        s=verticalColumn(5,uh,U,inertialMode_); if (!s) return s;
        s=verticalColumn(5,vh,V,inertialMode_); if (!s) return s;
    } else {
        s=vertical(5,uh,U); if (!s) return s;
        s=vertical(5,vh,V); if (!s) return s;
    }
    pointwise_->execute(g.Nkl,[&](std::size_t begin,std::size_t end) {
      for (std::size_t mode=begin;mode<end;++mode) {
        const double halfK=hasW ? std::hypot(g.k[mode],g.l[mode])/2 : 0;
        for (std::size_t j=0;j<g.Nj;++j) {
            const auto i=j+g.Nj*mode;
            const auto& f=factors_[i]; const auto n=scale(read(density.input(),i),f.ApmN);
            // Preserve the divergence addition before adding/subtracting density,
            // consuming the vertical-velocity projection without a separate pass.
            const auto d=hasW ? add(read(divergence.input(),i),multiply(read(temp.input(),i),{0,halfK})) : read(divergence.input(),i);
            auto ap=add(d,n),am=subtract(d,n);
            if (f.inertial) { ap=scale(subtract(read(U.input(),i),multiply(read(V.input(),i),{0,1})),.5); am=conjugate(ap); }
            b.Ap.data[i]=multiply(ap,conjugate(phase_[i])); b.Am.data[i]=multiply(am,phase_[i]);
        }
    }
    });
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::transformUVEtaToWaveVortex(WVRealVolumeConstView u,WVRealVolumeConstView v,WVRealVolumeConstView eta,double t,double t0,WVMutableCoefficients b) {
    auto s=outputs(b); if (!s) return s; s=mutableOutputOutsidePreparedState(b); if (!s) return s;
    for (auto a : {u,v,eta}) { s=volume(a); if (!s) return s; for (auto out : {b.Ap,b.Am,b.A0}) { s=disjoint(a.data,R_*sizeof(double),out.data,S_*sizeof(WVComplex64)); if (!s) return s; } }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=prepareProjectionPhaseForCall(t,t0); if (!s) return s;
    return projectFields(u.data,v.data,nullptr,eta.data,b);
}

WVKernelStatus WVTransformBoussinesqKernel::transformUVWEtaToWaveVortex(WVRealVolumeConstView u,WVRealVolumeConstView v,WVRealVolumeConstView w,WVRealVolumeConstView eta,double t,double t0,WVMutableCoefficients b) {
    auto s=outputs(b); if (!s) return s; s=mutableOutputOutsidePreparedState(b); if (!s) return s;
    for (auto a:{u,v,w,eta}) { s=volume(a); if (!s) return s; for (auto out:{b.Ap,b.Am,b.A0}) { s=disjoint(a.data,R_*sizeof(double),out.data,S_*sizeof(WVComplex64)); if (!s) return s; } }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=prepareProjectionPhaseForCall(t,t0); if (!s) return s;
    return projectFields(u.data,v.data,w.data,eta.data,b);
}

WVKernelStatus WVTransformBoussinesqKernel::reconstruct(const WVCoefficients& a,WVBoussinesqField field,
    WVBoussinesqDerivative derivative,WVBoussinesqComponent component,double* b,bool countPrimary,
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
    if (field==WVBoussinesqField::zetaX || field==WVBoussinesqField::zetaY) {
        const bool x=field==WVBoussinesqField::zetaX;
        auto s=reconstruct(a,x ? WVBoussinesqField::w : WVBoussinesqField::u,x ? WVBoussinesqDerivative::y : WVBoussinesqDerivative::z,component,b,countPrimary,metricComponent); if (!s) return s;
        auto* auxiliary=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 8)*R_;
        s=reconstruct(a,x ? WVBoussinesqField::v : WVBoussinesqField::w,x ? WVBoussinesqDerivative::z : WVBoussinesqDerivative::x,component,auxiliary,countPrimary,metricComponent); if (!s) return s;
        for (std::size_t i=0;i<R_;++i) b[i]-=auxiliary[i];
        return WVKernelStatus::ok();
    }
    const bool density=field==WVBoussinesqField::rhoE || field==WVBoussinesqField::rhoTotal;
    const bool totalDensity=field==WVBoussinesqField::rhoTotal;
    if (density) field=WVBoussinesqField::eta;
    if (field==WVBoussinesqField::ssh) field=WVBoussinesqField::pi;
    if (field==WVBoussinesqField::ssu) field=WVBoussinesqField::u;
    if (field==WVBoussinesqField::ssv) field=WVBoussinesqField::v;
    const bool G=field==WVBoussinesqField::w || field==WVBoussinesqField::eta;
    const bool dz=derivative==WVBoussinesqDerivative::z;
    kernel_detail::WVPreparedFieldCache::Entry* prepared=nullptr;
    if (countPrimary && fieldCache_) {
        auto s=fieldCache_->acquire({fieldPreparationView(a),static_cast<std::size_t>(field),
            static_cast<std::size_t>(component)},prepared);
        if (!s) return s;
    }
    const auto wave=prepared ? prepared->modal(0,S_) : modalView();
    const auto balancedView=prepared ? prepared->modal(S_,S_) : modalView(1);
    if (!prepared || !prepared->modalReady) {
        ++metrics_.coefficientAssemblyCount;
        // Validate once before dispatch; workers publish disjoint coefficients.
        if (field<WVBoussinesqField::u ||
            (field>WVBoussinesqField::qgpv && field!=WVBoussinesqField::zetaZ))
            return unsupported();
        pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
          for (std::size_t i=begin;i<end;++i) {
            const auto& f=factors_[i]; const auto mode=i/g.Nj; WVComplex64 p{},m{},z{};
            switch(field) {
                case WVBoussinesqField::u: p=f.UAp; m=conjugate(p); z=f.UA0; break;
                case WVBoussinesqField::v: p=f.VAp; m=conjugate(p); z=f.VA0; break;
                case WVBoussinesqField::w: p=f.WAp; m=p; break;
                case WVBoussinesqField::eta: p={f.NAp,0}; m={-f.NAp,0}; z={f.NA0,0}; break;
                case WVBoussinesqField::pi: p={f.NAp,0}; m={-f.NAp,0}; z={f.PA0,0}; break;
                case WVBoussinesqField::p: p={g.rho0*g.g*f.NAp,0}; m=scale(p,-1); z={g.rho0*g.g*f.PA0,0}; break;
                case WVBoussinesqField::psi: z={f.psi,0}; break;
                case WVBoussinesqField::qgpv: z={f.qgpv,0}; break;
                case WVBoussinesqField::zetaZ:
                    p=subtract(multiply(f.VAp,{0,g.k[mode]}),multiply(f.UAp,{0,g.l[mode]}));
                    m=subtract(multiply(conjugate(f.VAp),{0,g.k[mode]}),multiply(conjugate(f.UAp),{0,g.l[mode]}));
                    z=subtract(multiply(f.VA0,{0,g.k[mode]}),multiply(f.UA0,{0,g.l[mode]})); break;
                default: break; // Rejected before dispatch.
            }
            WVComplex64 value{},balanced{};
            if (selected(f,component,true)) value=add(multiply(multiply(p,a.Ap.data[i]),phase_[i]),multiply(multiply(m,a.Am.data[i]),conjugate(phase_[i])));
            if (selected(f,component,false)) balanced=multiply(z,a.A0.data[i]);
            if (!prepared && derivative==WVBoussinesqDerivative::x) { value=multiply(value,{0,g.k[mode]}); balanced=multiply(balanced,{0,g.k[mode]}); }
            if (!prepared && derivative==WVBoussinesqDerivative::y) { value=multiply(value,{0,g.l[mode]}); balanced=multiply(balanced,{0,g.l[mode]}); }
            write(wave,i,value); write(balancedView,i,balanced);
        }
        });
        if (prepared) prepared->modalReady=true;
    }
    const auto combined=prepared ? prepared->grid() : reconstructionGrid;
    if (!prepared || !prepared->gridReady) {
        const auto balancedGrid=gridView(1);
        ++metrics_.verticalPreparationCount;
        auto s=vertical(G ? 6 : 4,wave.input(),combined); if (!s) return s;
        ++metrics_.verticalPreparationCount;
        s=vertical(G ? 2 : 0,balancedView.input(),balancedGrid); if (!s) return s;
        pointwise_->execute(H_,[&](std::size_t begin,std::size_t end) {
            for (std::size_t i=begin;i<end;++i) write(combined,i,add(read(combined.input(),i),read(balancedGrid.input(),i)));
        });
        if (prepared) prepared->gridReady=true;
    } else ++metrics_.horizontalSpectrumReuseCount;
    auto horizontalInput=combined.input();
    if (prepared && (derivative==WVBoussinesqDerivative::x || derivative==WVBoussinesqDerivative::y)) {
        const auto derivativeGrid=reconstructionGrid;
        pointwise_->execute(H_,[&](std::size_t begin,std::size_t end) {
            for (std::size_t i=begin;i<end;++i) {
                const auto mode=i/g.Nz;
                const WVComplex64 multiplier={0,derivative==WVBoussinesqDerivative::x ? g.k[mode] : g.l[mode]};
                write(derivativeGrid,i,multiply(read(combined.input(),i),multiplier));
            }
        });
        horizontalInput=derivativeGrid.input();
    } else if (prepared && dz) {
        auto s=vertical(G ? 3 : 1,combined.input(),modalView()); if (!s) return s;
        if (G) pointwise_->execute(S_,[&](std::size_t begin,std::size_t end) {
            for (std::size_t i=begin;i<end;++i)
                write(modalView(),i,scale(read(modalView().input(),i),1/g.h_0[i%g.Nj]));
        });
        s=vertical(G ? 0 : 2,modalView().input(),reconstructionGrid); if (!s) return s;
        if (!G) pointwise_->execute(H_,[&](std::size_t begin,std::size_t end) {
            for (std::size_t i=begin;i<end;++i)
                write(reconstructionGrid,i,scale(read(reconstructionGrid.input(),i),-g.N2[i%g.Nz]/g.g));
        });
        ++metrics_.preparedVerticalDerivativeCount;
        horizontalInput=reconstructionGrid.input();
    }
    if (preparedSpectrum) { *preparedSpectrum=horizontalInput; return WVKernelStatus::ok(); }
    const bool consumeInInverse=consumer && (!dz || prepared);
    auto s=consumeInInverse ? horizontal_->inverseAndConsume(*horizontalWorkspace_,horizontalInput,{b,R_*sizeof(double)},*consumer) :
        horizontal_->inverse(*horizontalWorkspace_,horizontalInput,{b,R_*sizeof(double)}); if (!s) return s;
    // v4 defines vertical derivatives through the shared F/G calculus, even
    // for wave fields. Preserve that finite-resolution MATLAB operation.
    if (dz && !prepared) { s=verticalCalculus(b,G ? WVBoussinesqFamily::G : WVBoussinesqFamily::F,
        1,false,b,R_/geometry().Nz,false); if (!s) return s; }
    if (consumer && !consumeInInverse) consumer->consume(consumer->context,0,R_,b);
    if (density) {
        const double* eta=nullptr;
        if (dz) { auto* auxiliary=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 8)*R_; s=reconstruct(a,WVBoussinesqField::eta,WVBoussinesqDerivative::value,component,auxiliary,countPrimary,metricComponent); if (!s) return s; eta=auxiliary; }
        const auto plane=R_/g.Nz;
        for (std::size_t z=0;z<g.Nz;++z) for (std::size_t xy=0;xy<plane;++xy) {
            const auto i=xy+plane*z;
            b[i]=(g.rho0/g.g)*g.N2[z]*(b[i]+(dz ? g.dLnN2[z]*eta[i] : 0));
            if (totalDensity && derivative==WVBoussinesqDerivative::value) b[i]+=g.rho_nm0[z];
            if (totalDensity && dz) b[i]-=(g.rho0/g.g)*g.N2[z];
        }
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::transformStateField(const WVState& a,WVBoussinesqField field,WVRealVolumeView b,WVBoussinesqDerivative derivative,WVBoussinesqComponent component) {
    if (field<WVBoussinesqField::u || field>WVBoussinesqField::ssv || derivative<WVBoussinesqDerivative::value || derivative>WVBoussinesqDerivative::z || !valid(component) ||
        ((field==WVBoussinesqField::zetaX || field==WVBoussinesqField::zetaY) && derivative!=WVBoussinesqDerivative::value) ||
        (field==WVBoussinesqField::rhoTotal && component!=WVBoussinesqComponent::all)) return unsupported();
    auto s=validateStateForCall(a); if (!s) return s; s=volume({b.data,b.shape},surface(field)); if (!s) return s;
    const auto bytes=(surface(field) ? R_/geometry().Nz : R_)*sizeof(double);
    for (auto x : {a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) { s=disjoint(x.data,S_*sizeof(WVComplex64),b.data,bytes); if (!s) return s; }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    auto* fullField=real_.data()+(executionOptions_.streamedNonlinear ? 5 : 9)*R_;
    const auto metricComponent=component==WVBoussinesqComponent::all ?
        stateEvaluationComponent(a) : static_cast<std::size_t>(component);
    s=reconstruct(a.coefficients,field,derivative,component,
        surface(field) ? fullField : b.data,true,metricComponent); if (!s) return s;
    if (surface(field)) { const auto plane=R_/geometry().Nz; std::copy_n(fullField+R_-plane,plane,b.data); }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::transformCoefficientTendencyToUVWEta(
    const WVState& a,WVRealFieldBundleView& fields) {
    const auto& g=geometry();
    if (fields.shape.first!=g.Nx || fields.shape.second!=g.Ny ||
        fields.shape.third!=g.Nz || fields.shape.fourth!=4)
        return {WVKernelStatusCode::invalidShape,"Coefficient tendency requires [Nx,Ny,Nz,4] fields."};
    ++metrics_.derivedValidationCount;
    auto status=stateContents(a); if (!status) return status;
    const auto bytes=4*R_*sizeof(double);
    if (!addressFits(fields.data,bytes,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid coefficient tendency field storage."};
    for (const auto input:{a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) {
        status=disjoint(input.data,S_*sizeof(WVComplex64),fields.data,bytes);
        if (!status) return status;
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    status=prepareProjectionPhaseForCall(a.t,a.t0); if (!status) return status;
    const WVBoussinesqField names[]={WVBoussinesqField::u,WVBoussinesqField::v,
        WVBoussinesqField::w,WVBoussinesqField::eta};
    for (std::size_t channel=0;channel<4;++channel) {
        status=reconstruct(a.coefficients,names[channel],WVBoussinesqDerivative::value,
            WVBoussinesqComponent::all,fields.data+channel*R_,false);
        if (!status) return status;
        ++metrics_.tendencyReconstructionCount[channel];
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::combinePreparedHorizontalVorticity(
    WVBoussinesqField field,WVRealVolumeConstView first,WVRealVolumeConstView second,
    WVRealVolumeView output,WVBoussinesqComponent component) {
    if ((field!=WVBoussinesqField::zetaX && field!=WVBoussinesqField::zetaY) ||
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
        [static_cast<std::size_t>(WVBoussinesqDerivative::value)]
        [static_cast<std::size_t>(component)];
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::combinePreparedDensityZDerivative(
    WVBoussinesqField field,WVRealVolumeConstView etaZ,WVRealVolumeConstView eta,
    WVRealVolumeView output,WVBoussinesqComponent component) {
    if ((field!=WVBoussinesqField::rhoE && field!=WVBoussinesqField::rhoTotal) ||
        !valid(component) ||
        (field==WVBoussinesqField::rhoTotal && component!=WVBoussinesqComponent::all))
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
        if (field==WVBoussinesqField::rhoTotal) output.data[i]-=scale;
    }
    ++metrics_.fieldReconstructionCount[static_cast<std::size_t>(field)];
    ++metrics_.reconstructionCount[static_cast<std::size_t>(field)]
        [static_cast<std::size_t>(WVBoussinesqDerivative::z)]
        [static_cast<std::size_t>(component)];
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::evolveCoefficients(const WVState& a,WVMutableCoefficients b) {
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
WVKernelStatus WVTransformBoussinesqKernel::constrainCoefficients(WVMutableCoefficients a) const {
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
WVKernelStatus WVTransformBoussinesqKernel::nonlinearFluxAndFields(
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
            spatialTendency->shape.third!=shape.third || spatialTendency->shape.fourth!=4)
            return {WVKernelStatusCode::invalidShape,"Tiled nonlinear tendencies require four full volumes."};
        if (!addressFits(spatialTendency->data,4*R_*sizeof(double),alignof(double)))
            return {WVKernelStatusCode::invalidPointer,"Invalid tiled physical tendency storage."};
        status=disjoint(spatialTendency->data,4*R_*sizeof(double),fields.data,4*R_*sizeof(double)); if (!status) return status;
        for (const auto input:{state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0}) {
            status=disjoint(spatialTendency->data,4*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
        }
        for (std::size_t v=0;v<preparedStateViewCount_;++v)
            for (const auto input:{preparedStateViews_[v].coefficients.Ap,preparedStateViews_[v].coefficients.Am,preparedStateViews_[v].coefficients.A0}) {
                status=disjoint(spatialTendency->data,4*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!status) return status;
            }
        for (const auto output:{flux.Fp,flux.Fm,flux.F0}) {
            status=disjoint(spatialTendency->data,4*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!status) return status;
        }
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    status=preparePhaseForCall(state); if (!status) return status;
    WVRetainedAdvectionWork work;
    work.targets=4;work.fields={fields.data,4*R_*sizeof(double)};
    if (spatialTendency) work.tendencies={spatialTendency->data,4*R_*sizeof(double)};
    work.densityCorrection={geometry().dLnN2.data(),geometry().Nz*sizeof(double)};
    const WVBoussinesqField names[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    for (std::size_t f=0;f<4;++f) {
        status=reconstruct(state.coefficients,names[f],WVBoussinesqDerivative::value,WVBoussinesqComponent::all,
            nullptr,true,5,nullptr,&work.base[f]); if (!status) return status;
    }
    const std::size_t channels[]={0,1,3,2},slots[]={2,3,4,1};
    for (std::size_t t=0;t<work.targets;++t) {
        work.targetSpectra[t]=gridView(slots[t]);
        WVComplexInput prepared;
        status=reconstruct(state.coefficients,names[channels[t]],WVBoussinesqDerivative::z,WVBoussinesqComponent::all,
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
    return projectSpectralFields(work.targetSpectra[0].input(),work.targetSpectra[1].input(),
        work.targetSpectra[3].input(),work.targetSpectra[2].input(),gridView(),true,target);
}
WVKernelStatus WVTransformBoussinesqKernel::nonlinearFlux(const WVState& a,WVFlux& b,
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
        s=validateBundle(spatialTendency->data,spatialTendency->shape,4); if (!s) return s;
        for (auto input:{a.coefficients.Ap,a.coefficients.Am,a.coefficients.A0}) {
            s=disjoint(spatialTendency->data,4*R_*sizeof(double),input.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
        for (auto output:{b.Fp,b.Fm,b.F0}) {
            s=disjoint(spatialTendency->data,4*R_*sizeof(double),output.data,S_*sizeof(WVComplex64)); if (!s) return s;
        }
        if (preparedFields) {
            s=disjoint(spatialTendency->data,4*R_*sizeof(double),preparedFields->data,4*R_*sizeof(double)); if (!s) return s;
        }
    }
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    const WVBoussinesqField fields[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    const bool borrowed=executionOptions_.streamedNonlinear && preparedFields;
    const double* advectionFields=borrowed ? preparedFields->data : real_.data();
    const auto derivativeFor=[&](WVBoussinesqField field,WVBoussinesqDerivative derivative,
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
            WVBoussinesqComponent::all,scratch,true,5,consumer); if (!status) return status;
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
        if (!preparedFields) for (std::size_t i=0;i<4;++i) {
            s=reconstruct(a.coefficients,fields[i],WVBoussinesqDerivative::value,WVBoussinesqComponent::all,real_.data()+i*R_); if (!s) return s;
        }
        auto* flux=real_.data()+4*R_; auto* derivative=real_.data()+5*R_;
        const WVBoussinesqField targetFields[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::eta,WVBoussinesqField::w};
        const std::size_t outputChannels[]={0,1,3,2},spectralSlots[]={2,3,4,1};
        for (std::size_t targetIndex=0;targetIndex<4;++targetIndex) {
            const auto field=targetFields[targetIndex]; const auto outputChannel=outputChannels[targetIndex];
            std::fill_n(flux,R_,0);
            // These direct modal first derivatives use separate dependency keys
            // from grid Laplacians, which preserve sequential horizontal
            // calculus or the order-two vertical operator from retained values.
            for (std::size_t axis=0;axis<3;++axis) {
                kernel_detail::WVAdvectionConsumer advection{flux,advectionFields+axis*R_,
                    advectionFields+3*R_,geometry().dLnN2.data(),R_/geometry().Nz,field==WVBoussinesqField::eta && axis==2};
                const WVRealOutputConsumer consumer{&advection,kernel_detail::WVAdvectionConsumer::consume};
                bool consumed=false;
                const double* derivativeValues=nullptr;
                s=derivativeFor(field,static_cast<WVBoussinesqDerivative>(axis+1),
                    derivative,derivativeValues,flux,
                    executionOptions_.fusedDerivativeAdvection ? &consumer : nullptr,&consumed); if (!s) return s;
                if (!consumed) pointwise_->execute(R_,[&](std::size_t begin,std::size_t end) {
                    kernel_detail::WVAdvectionConsumer::consume(&advection,begin,end,derivativeValues);
                });
            }
            if (spatialTendency) std::copy_n(flux,R_,spatialTendency->data+outputChannel*R_);
            if (projectFlux) { s=horizontal_->forward(*horizontalWorkspace_,{flux,R_*sizeof(double)},gridView(spectralSlots[targetIndex])); if (!s) return s; }
        }
        return projectFlux ? projectSpectralFields(gridView(2).input(),gridView(3).input(),gridView(1).input(),
            gridView(4).input(),gridView(),true,target) : WVKernelStatus::ok();
    }
    if (preparedFields) { if (!borrowed) std::copy_n(preparedFields->data,4*R_,real_.data()); }
    else for (std::size_t i=0;i<4;++i) { s=reconstruct(a.coefficients,fields[i],WVBoussinesqDerivative::value,WVBoussinesqComponent::all,real_.data()+i*R_); if (!s) return s; }
    for (std::size_t targetIndex=0;targetIndex<4;++targetIndex) {
        const auto field=fields[targetIndex]; auto* flux=real_.data()+(4+targetIndex)*R_;
        std::fill_n(flux,R_,0);
        for (std::size_t axis=0;axis<3;++axis) {
            const double* derivativeValues=nullptr;
            s=derivativeFor(field,static_cast<WVBoussinesqDerivative>(axis+1),
                real_.data()+10*R_,derivativeValues,flux); if (!s) return s;
            pointwise_->execute(R_,[&](std::size_t begin,std::size_t end) {
                for (std::size_t i=begin;i<end;++i) {
                    const double correction=targetIndex==3 && axis==2 ? advectionFields[3*R_+i]*geometry().dLnN2[i/(R_/geometry().Nz)] : 0;
                    flux[i]-=advectionFields[axis*R_+i]*(derivativeValues[i]+correction);
                }
            });
        }
    }
    if (spatialTendency) std::copy_n(real_.data()+4*R_,4*R_,spatialTendency->data);
    if (!projectFlux) return WVKernelStatus::ok();
    return projectFields(real_.data()+4*R_,real_.data()+5*R_,real_.data()+6*R_,real_.data()+7*R_,target);
}
WVKernelStatus WVTransformBoussinesqKernel::totalEnergy(const WVCoefficients& a,double& value,WVBoussinesqComponent component) const {
    if (!valid(component)) return unsupported();
    auto s=coefficients(a); if (!s) return s;
    const auto norm=[](WVComplex64 x) { return x.real*x.real+x.imag*x.imag; }; double sum=0;
    for (std::size_t i=0;i<S_;++i) {
        if (selected(factors_[i],component,true)) sum+=factors_[i].waveEnergy*(norm(a.Ap.data[i])+norm(a.Am.data[i]));
        if (selected(factors_[i],component,false)) sum+=factors_[i].balancedEnergy*norm(a.A0.data[i]);
    }
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Boussinesq energy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::totalEnstrophy(const WVCoefficients& a,double& value) const {
    auto s=coefficients(a); if (!s) return s; double sum=0;
    for (std::size_t i=0;i<S_;++i) sum+=factors_[i].enstrophy*(a.A0.data[i].real*a.A0.data[i].real+a.A0.data[i].imag*a.A0.data[i].imag);
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Boussinesq enstrophy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::totalEnergySpatiallyIntegrated(const WVState& a,double& value,WVBoussinesqComponent component) {
    if (!valid(component)) return unsupported();
    auto s=validateStateForCall(a); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    kernel_detail::WVPreparedFieldCache::StandaloneScope cacheScope{fieldCache_.get(),!stateEvaluationActive_};
    s=preparePhaseForCall(a); if (!s) return s;
    const WVBoussinesqField fields[]={WVBoussinesqField::u,WVBoussinesqField::v,WVBoussinesqField::w,WVBoussinesqField::eta};
    for (std::size_t i=0;i<4;++i) { s=reconstruct(a.coefficients,fields[i],WVBoussinesqDerivative::value,component,real_.data()+i*R_); if (!s) return s; }
    double sum=0; const auto plane=R_/geometry().Nz;
    for (std::size_t i=0;i<R_;++i) sum+=geometry().z_int[i/plane]*(real_[i]*real_[i]+real_[R_+i]*real_[R_+i]+real_[2*R_+i]*real_[2*R_+i]+geometry().N2[i/plane]*real_[3*R_+i]*real_[3*R_+i])/(2*plane);
    if (!std::isfinite(sum)) return {WVKernelStatusCode::numericalFailure,"Boussinesq spatial energy overflow."};
    value=sum; return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::verticalCalculus(const double* values,WVBoussinesqFamily family,
    unsigned order,bool integral,double* result,std::size_t columns,bool verticalFirst) {
    const auto& g=geometry();
    const auto index=[&](std::size_t column,std::size_t z) {
        return verticalFirst ? z+g.Nz*column : column+columns*z;
    };
    return kernel_detail::applyStratifiedVerticalCalculus(g,columns,family==WVBoussinesqFamily::F,
        order,integral,[&](std::size_t column,std::size_t z) { return values[index(column,z)]; },
        [&](std::size_t column,std::size_t z,double value) { result[index(column,z)]=value; },
        gridView(),gridView(1),modalView(),
        [&](std::size_t operation,WVComplexInput input,WVComplexOutput output) { return vertical(operation,input,output); });
}
WVKernelStatus WVTransformBoussinesqKernel::applyVerticalCalculus(WVRealConstView a,bool inputIsF,
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
    return verticalCalculus(a.data,inputIsF ? WVBoussinesqFamily::F : WVBoussinesqFamily::G,
        order,integral,b.data,a.shape.columns,true);
}
WVKernelStatus WVTransformBoussinesqKernel::differentiateVertical(WVRealVolumeConstView a,WVBoussinesqFamily family,unsigned order,WVRealVolumeView b) {
    if ((family!=WVBoussinesqFamily::F && family!=WVBoussinesqFamily::G) || order<1 || order>4) return unsupported();
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return verticalCalculus(a.data,family,order,false,b.data,R_/geometry().Nz,false);
}
WVKernelStatus WVTransformBoussinesqKernel::integrateVertical(WVRealVolumeConstView a,WVBoussinesqFamily family,WVRealVolumeView b) {
    if (family!=WVBoussinesqFamily::F && family!=WVBoussinesqFamily::G) return unsupported();
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return verticalCalculus(a.data,family,1,true,b.data,R_/geometry().Nz,false);
}
WVKernelStatus WVTransformBoussinesqKernel::differentiateHorizontal(WVRealVolumeConstView a,bool xDerivative,WVRealVolumeView b) {
    return differentiateHorizontal(a,xDerivative,1,b);
}
WVKernelStatus WVTransformBoussinesqKernel::differentiateHorizontal(WVRealVolumeConstView a,
    bool xDerivative,unsigned order,WVRealVolumeView b) {
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    s=mutableOutputOutsidePreparedState(b.data,R_*sizeof(double)); if (!s) return s;
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    return horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{b.data,R_*sizeof(double)},xDerivative,order);
}
WVKernelStatus WVTransformBoussinesqKernel::waveModeVerticalStructureAtIndex(std::size_t z,WVRealView target) {
    if (z>=geometry().Nz || target.shape.rows!=geometry().Nj || target.shape.columns!=geometry().Nkl)
        return {WVKernelStatusCode::invalidShape,"Wave structure requires a valid z index and [Nj,Nkl] output."};
    if (!addressFits(target.data,S_*sizeof(double),alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid wave structure output."};
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    const auto& g=geometry();
    for (std::size_t j=0;j<g.Nj;++j) {
        for (std::size_t i=0;i<S_;++i) write(modalView(),i,{});
        for (std::size_t m=0;m<g.Nkl;++m) write(modalView(),j+g.Nj*m,{1,0});
        auto status=vertical(4,modalView().input(),gridView()); if (!status) return status;
        for (std::size_t m=0;m<g.Nkl;++m) target.data[j+g.Nj*m]=read(gridView().input(),z+g.Nz*m).real;
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVTransformBoussinesqKernel::advectScalarWithAdvectionFields(WVRealVolumeConstView a,WVRealFieldBundleConstView fields,bool antialias,WVRealVolumeView b,bool xyOnly) {
    auto s=volume(a); if (!s) return s; s=volume({b.data,b.shape}); if (!s) return s;
    const auto& g=geometry();
    if (fields.shape.first!=g.Nx || fields.shape.second!=g.Ny || fields.shape.third!=g.Nz || fields.shape.fourth!=3)
        return {WVKernelStatusCode::invalidShape,"Boussinesq advection requires [Nx,Ny,Nz,3] velocity fields."};
    if (!addressFits(fields.data,3*R_*sizeof(double),alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid advection field storage."};
    s=disjoint(a.data,R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    s=disjoint(fields.data,3*R_*sizeof(double),b.data,R_*sizeof(double)); if (!s) return s;
    ActiveCall guard(active_); if (!guard.entered) return reentrant();
    s=horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{real_.data(),R_*sizeof(double)},true); if (!s) return s;
    s=horizontal_->spatialDerivative(*horizontalWorkspace_,{a.data,R_*sizeof(double)},{real_.data()+R_,R_*sizeof(double)},false); if (!s) return s;
    if (!xyOnly) { s=verticalCalculus(a.data,WVBoussinesqFamily::F,1,false,real_.data()+2*R_,R_/geometry().Nz,false); if (!s) return s; }
    for (std::size_t i=0;i<R_;++i) b.data[i]=-fields.data[i]*real_[i]-fields.data[R_+i]*real_[R_+i]-(xyOnly ? 0.0 : fields.data[2*R_+i]*real_[2*R_+i]);
    if (antialias) {
        s=horizontal_->forward(*horizontalWorkspace_,{b.data,R_*sizeof(double)},gridView()); if (!s) return s;
        s=horizontal_->inverse(*horizontalWorkspace_,gridView().input(),{b.data,R_*sizeof(double)});
    }
    return s;
}
} // namespace wavevortex
