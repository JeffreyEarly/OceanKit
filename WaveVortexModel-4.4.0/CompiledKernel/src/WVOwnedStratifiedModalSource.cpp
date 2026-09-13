#include "WaveVortexKernel/WVOwnedStratifiedModalSource.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>

namespace wavevortex {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

WVKernelStatus invalid(const std::string& message) {
    return {WVKernelStatusCode::invalidConfiguration,message};
}

std::size_t product(std::size_t a,std::size_t b) {
    if (a && b>static_cast<std::size_t>(PTRDIFF_MAX)/a)
        throw std::overflow_error("Modal extent exceeds the addressable range.");
    return a*b;
}

bool close(double a,double b,double scale=1) {
    return std::abs(a-b)<=128*std::numeric_limits<double>::epsilon()*
        std::max({scale,std::abs(a),std::abs(b)});
}

WVKernelStatus validateSharedMatrix(const std::vector<double>& values,
    const WVStratifiedModalGeometry& g,bool inverse,bool isG) {
    const auto rows=inverse ? g.Nz : g.Nj;
    const auto columns=inverse ? g.Nj : g.Nz;
    if (values.size()!=product(rows,columns))
        return invalid("Shared modal matrix shape disagrees with the geometry.");
    for (std::size_t c=0;c<columns;++c) {
        double maximum=0;
        for (std::size_t r=0;r<rows;++r) {
            const double value=values[r+rows*c];
            if (!std::isfinite(value)) return invalid("Modal matrix contains a nonfinite value.");
            const auto z=inverse ? r : c,j=inverse ? c : r;
            if (isG && (g.j[j]==0 || z==0 || z+1==g.Nz) && value!=0)
                return invalid("Rigid-lid G matrices require zero boundaries and a zero barotropic mode.");
            if (inverse && !isG && g.j[j]==0 && !close(value,1))
                return invalid("The preconditioned barotropic F column must equal one.");
            maximum=std::max(maximum,std::abs(value));
        }
        if (inverse && !(isG && g.j[c]==0) && !close(maximum,1))
            return invalid("Reconstruction columns must use the persisted unit-maximum preconditioning.");
    }
    return WVKernelStatus::ok();
}

WVKernelStatus validateWaveMatrix(const std::vector<double>& values,
    const WVStratifiedModalGeometry& g,bool inverse,bool isG,bool cross) {
    const auto rows=inverse ? g.Nz : g.Nj;
    const auto columns=cross ? g.Nj : inverse ? g.Nj : g.Nz;
    const auto groups=g.K2unique.size();
    if (values.size()!=product(product(rows,columns),groups))
        return invalid("Boussinesq modal matrix shape disagrees with the geometry.");
    const auto plane=rows*columns;
    for (std::size_t group=0;group<groups;++group)
        for (std::size_t c=0;c<columns;++c) {
            double maximum=0;
            for (std::size_t r=0;r<rows;++r) {
                const double value=values[r+rows*c+plane*group];
                if (!std::isfinite(value)) return invalid("Wave matrix contains a nonfinite value.");
                if (cross) {
                    if ((g.j[r]==0 || g.j[c]==0) && value!=0)
                        return invalid("Wave G cross-projection must exclude the barotropic mode.");
                } else {
                    const auto z=inverse ? r : c,j=inverse ? c : r;
                    if (isG && (g.j[j]==0 || z==0 || z+1==g.Nz) && value!=0)
                        return invalid("Wave G matrices require rigid-lid boundaries and a zero barotropic mode.");
                    if (inverse && !isG && g.j[j]==0 && !close(value,1))
                        return invalid("Wave F barotropic columns must equal one.");
                }
                maximum=std::max(maximum,std::abs(value));
            }
            if (inverse && !(isG && g.j[c]==0) && !close(maximum,1))
                return invalid("Wave reconstruction must retain unit-maximum preconditioning.");
        }
    return WVKernelStatus::ok();
}

WVKernelStatus validate(WVStratifiedModalArrays& arrays) {
    auto& g=arrays.geometry;
    const bool boussinesq=g.transformClass=="WVTransformBoussinesq";
    if (g.modelVersion.rfind("4.",0)!=0)
        return invalid("Owned modal sources require v4 scientific data.");
    if (g.transformClass!="WVTransformStratifiedQG" &&
        g.transformClass!="WVTransformHydrostatic" && !boussinesq)
        return invalid("Unsupported owned stratified transform class.");
    if (g.Nx<2 || g.Ny<2 || g.Nz<3 || !g.Nj || g.Nj>=g.Nz || !g.Nkl ||
        g.Nx>static_cast<std::size_t>(INT64_MAX) ||
        g.Ny>static_cast<std::size_t>(INT64_MAX))
        return invalid("Invalid horizontal or rigid-lid modal dimensions.");
    product(product(product(g.Nx,g.Ny),g.Nz),sizeof(double));
    product(product(g.Nj,g.Nkl),sizeof(WVComplex64));
    if (!std::isfinite(g.Lx) || !std::isfinite(g.Ly) || !std::isfinite(g.Lz) ||
        !std::isfinite(g.g) || !std::isfinite(g.rho0) ||
        !std::isfinite(g.latitude) || !std::isfinite(g.rotationRate) ||
        !std::isfinite(g.planetaryRadius) || g.Lx<=0 || g.Ly<=0 || g.Lz<=0 ||
        g.g<=0 || g.rho0<=0 || g.rotationRate<=0 || g.planetaryRadius<=0)
        return invalid("Invalid owned modal physical scalar.");
    if (std::abs(g.latitude)<5 || std::abs(g.latitude)>85)
        return invalid("Absolute latitude must lie in [5,85].");
    if (g.x.size()!=g.Nx || g.y.size()!=g.Ny || g.z.size()!=g.Nz ||
        g.j.size()!=g.Nj || g.k.size()!=g.Nkl || g.l.size()!=g.Nkl ||
        g.N2.size()!=g.Nz || g.rho_nm0.size()!=g.Nz ||
        g.dLnN2.size()!=g.Nz || g.P0.size()!=g.Nj || g.Q0.size()!=g.Nj ||
        g.h_0.size()!=g.Nj || g.z_int.size()!=g.Nz || g.modes.size()!=g.Nkl)
        return invalid("Owned modal vector shape disagrees with its declared dimensions.");
    for (std::size_t i=0;i<g.Nx;++i)
        if (!std::isfinite(g.x[i]) || !close(g.x[i],g.Lx*static_cast<double>(i)/g.Nx,g.Lx))
            return invalid("x must use the periodic origin and spacing.");
    for (std::size_t i=0;i<g.Ny;++i)
        if (!std::isfinite(g.y[i]) || !close(g.y[i],g.Ly*static_cast<double>(i)/g.Ny,g.Ly))
            return invalid("y must use the periodic origin and spacing.");
    if (std::abs(g.z.front()+g.Lz)>1e-10*g.Lz || std::abs(g.z.back())>1e-10*g.Lz)
        return invalid("Vertical grid must include bottom -Lz and surface zero.");
    long double integral=0;
    for (std::size_t i=0;i<g.Nz;++i) {
        if (!std::isfinite(g.z[i]) || !std::isfinite(g.N2[i]) ||
            !std::isfinite(g.rho_nm0[i]) || !std::isfinite(g.dLnN2[i]) ||
            !std::isfinite(g.z_int[i]) || (i && g.z[i]<=g.z[i-1]) ||
            g.N2[i]<=0 || g.rho_nm0[i]<=0 || (i && g.rho_nm0[i]>g.rho_nm0[i-1]))
            return invalid("Vertical coordinates and stable stratification are inconsistent.");
        integral+=g.z_int[i];
    }
    if (std::abs(integral-g.Lz)>1e-10*g.Lz)
        return invalid("Vertical quadrature must integrate a constant to Lz.");
    for (std::size_t i=0;i<g.Nj;++i) {
        if (!std::isfinite(g.j[i]) || g.j[i]!=std::floor(g.j[i]) || g.j[i]<0 ||
            g.j[i]>=static_cast<double>(g.Nz-1) || (i && g.j[i]<=g.j[i-1]))
            return invalid("Retained vertical keys must be unique increasing rigid-lid mode indices.");
        if (!std::isfinite(g.P0[i]) || !std::isfinite(g.Q0[i]) ||
            !std::isfinite(g.h_0[i]) || g.P0[i]<=0 || g.Q0[i]<=0 || g.h_0[i]<=0)
            return invalid("Preconditioners and equivalent depths must be positive.");
    }
    if (g.j.front()!=0 || !close(g.P0.front(),1) || !close(g.Q0.front(),1) ||
        !close(g.h_0.front(),1))
        return invalid("The barotropic mode requires its documented unit sentinel values.");
    std::set<std::pair<std::int64_t,std::int64_t>> seen;
    const double cutoff=(2.0/3.0)*(2*pi*static_cast<double>(g.Nx/2)/g.Lx);
    for (std::size_t i=0;i<g.Nkl;++i) {
        if (!std::isfinite(g.k[i]) || !std::isfinite(g.l[i]))
            return invalid("Horizontal wavenumbers must be finite.");
        const auto key=g.modes[i];
        const double k=g.k[i]*g.Lx/(2*pi),l=g.l[i]*g.Ly/(2*pi);
        const auto kHalf=static_cast<std::int64_t>((g.Nx-1)/2);
        const auto lHalf=static_cast<std::int64_t>((g.Ny-1)/2);
        if (key.k < -kHalf || key.k > kHalf || key.l < -lHalf || key.l > lHalf ||
            !close(k,static_cast<double>(key.k)) || !close(l,static_cast<double>(key.l)))
            return invalid("Horizontal keys are outside the Nyquist-excluding integer grid.");
        if (key.l<0 || (key.l==0 && key.k<0) || !seen.insert({key.k,key.l}).second)
            return invalid("Horizontal keys must be unique primary Hermitian representatives.");
        const auto radius=std::hypot(g.k[i],g.l[i]);
        if (g.shouldAntialias && radius>cutoff && !close(radius,cutoff,cutoff))
            return invalid("Retained key violates the declared antialias mask.");
    }
    auto status=validateSharedMatrix(arrays.PF0inv,g,true,false); if (!status) return status;
    status=validateSharedMatrix(arrays.QG0inv,g,true,true); if (!status) return status;
    status=validateSharedMatrix(arrays.PF0,g,false,false); if (!status) return status;
    status=validateSharedMatrix(arrays.QG0,g,false,true); if (!status) return status;
    if (!boussinesq) {
        if (!g.K2unique.empty() || !g.waveGroup.empty() || !g.h_pm.empty() ||
            !g.Ppm.empty() || !g.Qpm.empty() || !arrays.PFpmInv.empty() ||
            !arrays.QGpmInv.empty() || !arrays.PFpm.empty() ||
            !arrays.QGpm.empty() || !arrays.QGwg.empty())
            return invalid("Wave modal data require a Boussinesq transform.");
        return WVKernelStatus::ok();
    }
    const auto groupCount=g.K2unique.size();
    if (!groupCount || groupCount>g.Nkl || g.waveGroup.size()!=g.Nkl ||
        g.h_pm.size()!=product(g.Nj,g.Nkl) ||
        g.Ppm.size()!=product(g.Nj,groupCount) ||
        g.Qpm.size()!=product(g.Nj,groupCount))
        return invalid("Boussinesq group arrays disagree with the geometry.");
    for (std::size_t i=0;i<groupCount;++i)
        if (!std::isfinite(g.K2unique[i]) || g.K2unique[i]<0 ||
            (i && g.K2unique[i]<=g.K2unique[i-1]))
            return invalid("Wave group radii must be nonnegative and strictly increasing.");
    if (g.K2unique.front()!=0) return invalid("Wave groups must include the inertial zero wavenumber.");
    std::vector<std::size_t> first(groupCount,g.Nkl);
    for (std::size_t i=0;i<g.Nkl;++i) {
        const auto group=g.waveGroup[i];
        if (group>=groupCount) return invalid("Wave membership index is outside K2unique.");
        first[group]=std::min(first[group],i);
        const double radius=g.k[i]*g.k[i]+g.l[i]*g.l[i];
        if (!close(radius,g.K2unique[group],std::max(radius,g.K2unique[group])))
            return invalid("Persisted wave membership disagrees with its horizontal wavenumber.");
    }
    if (g.waveGroup.front()!=0 || std::find(first.begin(),first.end(),g.Nkl)!=first.end())
        return invalid("Wave groups require complete membership and a first inertial column.");
    for (std::size_t group=0;group<groupCount;++group)
        for (std::size_t j=0;j<g.Nj;++j) {
            const auto scale=j+g.Nj*group;
            if (!std::isfinite(g.Ppm[scale]) || !std::isfinite(g.Qpm[scale]) ||
                g.Ppm[scale]<=0 || g.Qpm[scale]<=0 ||
                (g.j[j]==0 && (g.Ppm[scale]!=1 || g.Qpm[scale]!=1)))
                return invalid("Wave preconditioners require positive values and unit barotropic sentinels.");
        }
    for (std::size_t mode=0;mode<g.Nkl;++mode)
        for (std::size_t j=0;j<g.Nj;++j) {
            const double depth=g.h_pm[j+g.Nj*mode];
            if (!std::isfinite(depth) || depth<=0 || (g.j[j]==0 && depth!=1))
                return invalid("Wave depths require positive values and unit barotropic sentinels.");
            if (depth!=g.h_pm[j+g.Nj*first[g.waveGroup[mode]]])
                return invalid("Equivalent depths differ within one persisted wave group.");
        }
    const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
    for (double N2:g.N2) if (N2==f*f)
        return invalid("Boussinesq wave projection is singular where N2 equals f squared.");
    status=validateWaveMatrix(arrays.PFpmInv,g,true,false,false); if (!status) return status;
    status=validateWaveMatrix(arrays.QGpmInv,g,true,true,false); if (!status) return status;
    status=validateWaveMatrix(arrays.PFpm,g,false,false,false); if (!status) return status;
    status=validateWaveMatrix(arrays.QGpm,g,false,true,false); if (!status) return status;
    return validateWaveMatrix(arrays.QGwg,g,false,true,true);
}

std::uint64_t nextIdentity() {
    static std::atomic<std::uint64_t> sequence{1};
    auto next=sequence.load();
    do {
        if (next==UINT64_MAX) throw std::overflow_error("Modal source identities exhausted.");
    } while (!sequence.compare_exchange_weak(next,next+1));
    return next;
}
} // namespace

WVKernelStatus WVOwnedStratifiedModalSource::create(WVStratifiedModalArrays arrays,
    std::shared_ptr<const WVOwnedStratifiedModalSource>& result) {
    try {
        auto status=validate(arrays); if (!status) return status;
        auto candidate=std::shared_ptr<WVOwnedStratifiedModalSource>(new WVOwnedStratifiedModalSource);
        candidate->arrays_=std::move(arrays);
        const auto& g=candidate->arrays_.geometry;
        candidate->groups_.resize(std::max(std::size_t{1},g.K2unique.size()));
        for (std::size_t i=0;i<candidate->groups_.size();++i) candidate->groups_[i].identity=i;
        for (std::size_t i=0;i<g.Nkl;++i)
            candidate->groups_[g.waveGroup.empty() ? 0 : g.waveGroup[i]].columns.push_back(i);
        candidate->sourceIdentity_="wave-vortex-owned-stratified-modal-source-v1/source/"+
            std::to_string(nextIdentity());
        candidate->modeSetIdentity_=candidate->sourceIdentity_+"/horizontal-modes";
        result=std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Owned modal source allocation failed."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    }
}

std::size_t WVOwnedStratifiedModalSource::persistentBytes() const noexcept {
    const auto& g=arrays_.geometry;
    std::size_t bytes=sizeof(*this)+sourceIdentity_.capacity()+modeSetIdentity_.capacity()+
        g.transformClass.capacity()+g.modelVersion.capacity()+g.modes.capacity()*sizeof(WVRetainedModeKey);
    for (const auto* values:{&g.x,&g.y,&g.z,&g.j,&g.k,&g.l,&g.N2,&g.rho_nm0,
        &g.dLnN2,&g.P0,&g.Q0,&g.h_0,&g.z_int,&arrays_.PF0inv,&arrays_.QG0inv,
        &arrays_.PF0,&arrays_.QG0}) bytes+=values->capacity()*sizeof(double);
    for (const auto* values:{&g.K2unique,&g.h_pm,&g.Ppm,&g.Qpm,&arrays_.PFpmInv,
        &arrays_.QGpmInv,&arrays_.PFpm,&arrays_.QGpm,&arrays_.QGwg})
        bytes+=values->capacity()*sizeof(double);
    bytes+=g.waveGroup.capacity()*sizeof(std::size_t)+groups_.capacity()*sizeof(WVScientificModalGroup);
    for (const auto& group:groups_) bytes+=group.columns.capacity()*sizeof(std::size_t);
    return bytes;
}

WVKernelStatus WVOwnedStratifiedModalSource::matrixView(
    WVStratifiedScientificMatrix matrix,WVVerticalMatrixRecord& result) const {
    try {
        const char* name; const char* modal; const char* grid;
        const std::vector<double>* values; bool inverse;
        switch (matrix) {
            case WVStratifiedScientificMatrix::PF0inv: name="PF0inv"; modal="PF0-modal"; grid="F-grid"; values=&arrays_.PF0inv; inverse=true; break;
            case WVStratifiedScientificMatrix::QG0inv: name="QG0inv"; modal="QG0-modal"; grid="G-grid"; values=&arrays_.QG0inv; inverse=true; break;
            case WVStratifiedScientificMatrix::PF0: name="PF0"; modal="PF0-modal"; grid="F-grid"; values=&arrays_.PF0; inverse=false; break;
            case WVStratifiedScientificMatrix::QG0: name="QG0"; modal="QG0-modal"; grid="G-grid"; values=&arrays_.QG0; inverse=false; break;
            default: return invalid("Unknown persisted modal matrix.");
        }
        const auto& g=geometry();
        const auto rows=inverse ? g.Nz : g.Nj,columns=inverse ? g.Nj : g.Nz;
        WVVerticalMatrixRecord candidate{{sourceIdentity_,name},
            {values->data(),rows,columns,1,rows,values->size()*sizeof(double)},
            inverse ? WVMatrixAction::reconstruction : WVMatrixAction::projection,
            inverse ? modal : grid,inverse ? grid : modal};
        result=std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Scientific matrix view allocation failed."};
    }
}

WVKernelStatus WVOwnedStratifiedModalSource::prepareVertical(
    WVStratifiedModalOperator operation,WVComplexLayout input,WVComplexLayout output,
    std::unique_ptr<WVVerticalMatrixBackend> backend,
    std::unique_ptr<WVPreparedVerticalOperator>& result,WVAccumulation accumulation) const {
    if (operation>=WVStratifiedModalOperator::reconstructFw &&
        operation<=WVStratifiedModalOperator::projectWaveVerticalVelocity)
        return prepareWaveVertical(operation,std::move(input),std::move(output),
            std::move(backend),result,accumulation);
    try {
        const auto& g=geometry(); const char* name; const char* inFamily;
        const char* outFamily; WVMatrixAction action;
        switch (operation) {
            case WVStratifiedModalOperator::reconstructF: name="Finv"; inFamily="F-modal"; outFamily="F-grid"; action=WVMatrixAction::reconstruction; break;
            case WVStratifiedModalOperator::projectF: name="F"; inFamily="F-grid"; outFamily="F-modal"; action=WVMatrixAction::projection; break;
            case WVStratifiedModalOperator::reconstructG: name="Ginv"; inFamily="G-modal"; outFamily="G-grid"; action=WVMatrixAction::reconstruction; break;
            case WVStratifiedModalOperator::projectG: name="G"; inFamily="G-grid"; outFamily="G-modal"; action=WVMatrixAction::projection; break;
            case WVStratifiedModalOperator::GToF: name="F*Ginv"; inFamily="G-modal"; outFamily="F-modal"; action=WVMatrixAction::crossFamily; break;
            case WVStratifiedModalOperator::FToG: name="G*Finv"; inFamily="F-modal"; outFamily="G-modal"; action=WVMatrixAction::crossFamily; break;
            default: return invalid("Unknown modal operator.");
        }
        if (input.family!=inFamily || output.family!=outFamily || input.columns!=g.Nkl ||
            output.columns!=g.Nkl || input.modeSet!=modeSetIdentity_ || output.modeSet!=modeSetIdentity_)
            return invalid("Modal field family or ordered mode identity disagrees with the scientific source.");
        const auto rows=action==WVMatrixAction::reconstruction ? g.Nz : g.Nj;
        const auto columns=action==WVMatrixAction::projection ? g.Nz : g.Nj;
        std::vector<double> values(product(rows,columns));
        for (std::size_t c=0;c<columns;++c) for (std::size_t r=0;r<rows;++r) {
            double value=0;
            switch (operation) {
                case WVStratifiedModalOperator::reconstructF: value=arrays_.PF0inv[r+g.Nz*c]*g.P0[c]; break;
                case WVStratifiedModalOperator::reconstructG: value=arrays_.QG0inv[r+g.Nz*c]*g.Q0[c]; break;
                case WVStratifiedModalOperator::projectF: value=arrays_.PF0[r+g.Nj*c]/g.P0[r]; break;
                case WVStratifiedModalOperator::projectG: value=arrays_.QG0[r+g.Nj*c]/g.Q0[r]; break;
                case WVStratifiedModalOperator::GToF:
                    for (std::size_t z=0;z<g.Nz;++z)
                        value+=(arrays_.PF0[r+g.Nj*z]/g.P0[r])*(arrays_.QG0inv[z+g.Nz*c]*g.Q0[c]);
                    break;
                case WVStratifiedModalOperator::FToG:
                    for (std::size_t z=0;z<g.Nz;++z)
                        value+=(arrays_.QG0[r+g.Nj*z]/g.Q0[r])*(arrays_.PF0inv[z+g.Nz*c]*g.P0[c]);
                    break;
                default: break;
            }
            values[r+rows*c]=value;
        }
        WVVerticalSpecification spec; spec.sourceIdentity=sourceIdentity_; spec.action=action;
        spec.Nz=g.Nz; spec.Nj=g.Nj; spec.input=std::move(input); spec.output=std::move(output);
        spec.accumulation=accumulation;
        spec.matrices={{{sourceIdentity_,name},
            {values.data(),rows,columns,1,rows,values.size()*sizeof(double)},action,inFamily,outFamily}};
        WVVerticalGroup group; group.identity=0; group.matrix=0; group.modes.resize(g.Nkl);
        for (std::size_t i=0;i<g.Nkl;++i) group.modes[i]=i;
        spec.groups.push_back(std::move(group));
        return WVPreparedVerticalOperator::create(spec,std::move(backend),result);
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Modal operator preparation allocation failed."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    }
}

WVKernelStatus WVOwnedStratifiedModalSource::prepareWaveVertical(
    WVStratifiedModalOperator operation,WVComplexLayout input,WVComplexLayout output,
    std::unique_ptr<WVVerticalMatrixBackend> backend,
    std::unique_ptr<WVPreparedVerticalOperator>& result,WVAccumulation accumulation) const {
    try {
        const auto& g=geometry();
        if (g.transformClass!="WVTransformBoussinesq" || groups_.empty())
            return {WVKernelStatusCode::unsupportedOperation,"Wave matrices require a Boussinesq source."};
        const char* name; const char* inFamily; const char* outFamily; WVMatrixAction action;
        switch (operation) {
            case WVStratifiedModalOperator::reconstructFw: name="FwInv"; inFamily="Fw-modal"; outFamily="F-grid"; action=WVMatrixAction::reconstruction; break;
            case WVStratifiedModalOperator::projectFw: name="Fw"; inFamily="F-grid"; outFamily="Fw-modal"; action=WVMatrixAction::projection; break;
            case WVStratifiedModalOperator::reconstructGw: name="GwInv"; inFamily="Gw-modal"; outFamily="G-grid"; action=WVMatrixAction::reconstruction; break;
            case WVStratifiedModalOperator::projectGw: name="Gw"; inFamily="G-grid"; outFamily="Gw-modal"; action=WVMatrixAction::projection; break;
            case WVStratifiedModalOperator::balancedGToWaveG: name="Gwg"; inFamily="G-modal"; outFamily="Gw-modal"; action=WVMatrixAction::crossFamily; break;
            case WVStratifiedModalOperator::projectWaveDivergence: name="GwDdelta"; inFamily="F-grid"; outFamily="Gw-modal"; action=WVMatrixAction::projection; break;
            case WVStratifiedModalOperator::projectWaveVerticalVelocity: name="GwVerticalVelocity"; inFamily="G-grid"; outFamily="Gw-modal"; action=WVMatrixAction::projection; break;
            default: return invalid("Unknown Boussinesq operator.");
        }
        if (input.family!=inFamily || output.family!=outFamily || input.columns!=g.Nkl ||
            output.columns!=g.Nkl || input.modeSet!=modeSetIdentity_ || output.modeSet!=modeSetIdentity_)
            return invalid("Wave field family or exact mode identity disagrees with the source.");
        const auto rows=action==WVMatrixAction::reconstruction ? g.Nz : g.Nj;
        const auto columns=action==WVMatrixAction::projection ? g.Nz : g.Nj;
        const auto plane=product(rows,columns),wavePlane=product(g.Nz,g.Nj);
        std::vector<double> values(product(plane,groups_.size()));
        std::vector<double> delta;
        const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
        if (operation==WVStratifiedModalOperator::projectWaveDivergence) {
            delta.resize(product(g.Nz,g.Nz));
            for (std::size_t c=0;c<g.Nz;++c) for (std::size_t r=0;r<g.Nz;++r) {
                double value=0;
                for (std::size_t j=0;j<g.Nj;++j)
                    value+=arrays_.QG0inv[r+g.Nz*j]*(g.Q0[j]/g.P0[j])*arrays_.PF0[j+g.Nj*c];
                delta[r+g.Nz*c]=g.N2[r]/(g.N2[r]-f*f)*value;
            }
        }
        WVVerticalSpecification spec; spec.sourceIdentity=sourceIdentity_; spec.action=action;
        spec.Nz=g.Nz; spec.Nj=g.Nj; spec.input=std::move(input); spec.output=std::move(output);
        spec.accumulation=accumulation;
        for (std::size_t group=0;group<groups_.size();++group) {
            const auto offset=wavePlane*group,scaleOffset=g.Nj*group;
            for (std::size_t c=0;c<columns;++c) for (std::size_t r=0;r<rows;++r) {
                double value=0;
                switch (operation) {
                    case WVStratifiedModalOperator::reconstructFw: value=arrays_.PFpmInv[r+g.Nz*c+offset]*g.Ppm[c+scaleOffset]; break;
                    case WVStratifiedModalOperator::projectFw: value=arrays_.PFpm[r+g.Nj*c+offset]/g.Ppm[r+scaleOffset]; break;
                    case WVStratifiedModalOperator::reconstructGw: value=arrays_.QGpmInv[r+g.Nz*c+offset]*g.Qpm[c+scaleOffset]; break;
                    case WVStratifiedModalOperator::projectGw: value=arrays_.QGpm[r+g.Nj*c+offset]/g.Qpm[r+scaleOffset]; break;
                    case WVStratifiedModalOperator::balancedGToWaveG: value=arrays_.QGwg[r+g.Nj*c+g.Nj*g.Nj*group]*g.Q0[c]/g.Qpm[r+scaleOffset]; break;
                    case WVStratifiedModalOperator::projectWaveDivergence:
                        for (std::size_t z=0;z<g.Nz;++z)
                            value+=(arrays_.QGpm[r+g.Nj*z+offset]/g.Qpm[r+scaleOffset])*delta[z+g.Nz*c];
                        break;
                    case WVStratifiedModalOperator::projectWaveVerticalVelocity:
                        value=(arrays_.QGpm[r+g.Nj*c+offset]/g.Qpm[r+scaleOffset])*g.g/(g.N2[c]-f*f);
                        break;
                    default: break;
                }
                values[r+rows*c+plane*group]=value;
            }
            spec.matrices.push_back({{sourceIdentity_,std::string(name)+"/"+
                std::to_string(groups_[group].identity)},
                {values.data()+plane*group,rows,columns,1,rows,plane*sizeof(double)},
                action,inFamily,outFamily});
            spec.groups.push_back({groups_[group].identity,group,groups_[group].columns});
        }
        return WVPreparedVerticalOperator::create(spec,std::move(backend),result);
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Wave operator preparation allocation failed."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    }
}

WVKernelStatus WVOwnedStratifiedModalSource::horizontalSpecification(
    std::size_t planes,WVComplexRepresentation representation,const std::string& family,
    WVRetainedHorizontalSpecification& result) const {
    try {
        const auto& g=geometry();
        if (!planes || planes>g.Nz || family.empty() ||
            (representation!=WVComplexRepresentation::split &&
             representation!=WVComplexRepresentation::interleaved))
            return invalid("Invalid modal horizontal field layout.");
        WVRetainedHorizontalSpecification candidate;
        candidate.grid={g.Nx,g.Ny,planes,1,g.Nx,product(g.Nx,g.Ny),family};
        candidate.retained={planes,g.Nkl,1,planes,representation,family,modeSetIdentity_};
        candidate.Lx=g.Lx; candidate.Ly=g.Ly; candidate.modes=g.modes;
        result=std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Modal horizontal specification allocation failed."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    }
}

} // namespace wavevortex
