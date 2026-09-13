#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"
#include "WVForcingDiagnosticWorkspace.hpp"
#include "WVScopedStateEvaluation.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <new>
#include <set>
#include <sstream>
#include <type_traits>

namespace wavevortex::runtime {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

WVComplex64 add(WVComplex64 a, WVComplex64 b) noexcept { return {a.real+b.real,a.imag+b.imag}; }
WVComplex64 multiply(WVComplex64 a, WVComplex64 b) noexcept { return {a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real}; }
WVComplex64 multiply(WVComplex64 a, double b) noexcept { return {a.real*b,a.imag*b}; }
WVComplex64 conjugate(WVComplex64 a) noexcept { return {a.real,-a.imag}; }

std::size_t stageRank(WVForcingStage stage) noexcept { return static_cast<std::size_t>(stage); }

template <typename T>
std::size_t vectorBytes(const std::vector<T>& values) noexcept { return values.capacity()*sizeof(T); }

double vanishingFilter(double value, double cutoff, double maximum) noexcept {
    value = std::abs(value);
    // A collapsed zero-mode interval must preserve the uniform mode.
    if (maximum == 0.0) return value > 0.0 ? 1.0 : 0.0;
    if (value < cutoff) return 0.0;
    if (value > maximum) return 1.0;
    if (maximum == cutoff) return value >= maximum ? 1.0 : 0.0;
    const double ratio = (value-maximum)/(value-cutoff);
    return std::exp(-(ratio*ratio));
}

WVComplex64 normalizedTerrainCoefficient(
    const WVPseudoTopographicConfiguration& record,
    std::int64_t kMode,
    std::int64_t lMode) {
    const auto Nx = record.topographicShape.rows;
    const auto Ny = record.topographicShape.columns;
    std::complex<double> sum(0.0,0.0);
    for (std::size_t iY = 0; iY < Ny; ++iY) {
        for (std::size_t iX = 0; iX < Nx; ++iX) {
            const double angle = -2.0*pi*(static_cast<double>(kMode)*static_cast<double>(iX)/static_cast<double>(Nx) +
                static_cast<double>(lMode)*static_cast<double>(iY)/static_cast<double>(Ny));
            sum += record.topographicHeight[iX+Nx*iY]*std::complex<double>(std::cos(angle),std::sin(angle));
        }
    }
    sum /= static_cast<double>(Nx*Ny);
    return {sum.real(),sum.imag()};
}

WVFlux fluxViews(std::vector<WVComplex64>& storage, WVShape2D shape) {
    const auto count = shape.elementCount();
    return {{storage.data(),shape},{storage.data()+count,shape},{storage.data()+2*count,shape}};
}

void addFlux(const WVFlux& source, WVFlux& destination) {
    const auto count = source.Fp.shape.elementCount();
    for (std::size_t index = 0; index < count; ++index) {
        destination.Fp.data[index] = add(destination.Fp.data[index],source.Fp.data[index]);
        destination.Fm.data[index] = add(destination.Fm.data[index],source.Fm.data[index]);
        destination.F0.data[index] = add(destination.F0.data[index],source.F0.data[index]);
    }
}

WVKernelStatus validateMutableCoefficients(WVShape2D expected, const WVMutableCoefficients& coefficients) {
    for (const auto& value : {coefficients.Ap,coefficients.Am,coefficients.A0}) {
        if (value.data == nullptr) return {WVKernelStatusCode::invalidPointer,"Mutable coefficient storage is null."};
        if (value.shape.rows != expected.rows || value.shape.columns != expected.columns) return {WVKernelStatusCode::invalidShape,"Mutable coefficient shape does not match the forcing engine."};
    }
    return WVKernelStatus::ok();
}

class WVWaveVortexCoefficientErrorPolicy final : public WVIntegrationErrorPolicy {
public:
    static WVKernelStatus create(
        const WVTransformConstantStratificationDescriptor& descriptor,
        double absoluteToleranceScale,
        std::unique_ptr<WVIntegrationErrorPolicy>& result) {
        if (!std::isfinite(absoluteToleranceScale) || absoluteToleranceScale <= 0.0) {
            return {WVKernelStatusCode::invalidConfiguration,"Adaptive absolute-tolerance scale must be finite and positive."};
        }
        try {
            auto policy = std::unique_ptr<WVWaveVortexCoefficientErrorPolicy>(new WVWaveVortexCoefficientErrorPolicy);
            const auto& configuration = descriptor.configuration();
            const auto& horizontalModes = descriptor.fourierModes();
            const auto& verticalModes = descriptor.verticalModes();
            const auto shape = descriptor.spectralShape();
            policy->shape_ = shape;
            policy->waveTolerance_.assign(shape.elementCount(),1.0);
            policy->vortexTolerance_.assign(shape.elementCount(),1.0);

            std::vector<double> uniqueKh;
            uniqueKh.reserve(horizontalModes.size());
            for (const auto& mode : horizontalModes) uniqueKh.push_back(std::abs(mode.Kh));
            std::sort(uniqueKh.begin(),uniqueKh.end());
            uniqueKh.erase(std::unique(uniqueKh.begin(),uniqueKh.end()),uniqueKh.end());
            if (uniqueKh.size() < 2) return {WVKernelStatusCode::invalidConfiguration,"Adaptive tolerances require at least two distinct horizontal radial wavenumbers."};
            double deltaK = 0.0;
            for (std::size_t index = 1; index < uniqueKh.size(); ++index) deltaK = std::max(deltaK,uniqueKh[index]-uniqueKh[index-1]);
            if (!(deltaK > 0.0) || !std::isfinite(deltaK)) return {WVKernelStatusCode::invalidConfiguration,"Adaptive radial tolerance spacing is invalid."};
            const double maximumKh = uniqueKh.back();
            std::vector<double> radialCenters;
            for (double center = 0.0; center <= maximumKh+0.5*deltaK; center += deltaK) radialCenters.push_back(center);
            std::vector<std::size_t> radialBin(horizontalModes.size(),radialCenters.size());
            std::vector<std::size_t> radialCounts(radialCenters.size(),0);
            for (std::size_t modeIndex = 0; modeIndex < horizontalModes.size(); ++modeIndex) {
                const double Kh = horizontalModes[modeIndex].Kh;
                for (std::size_t bin = 0; bin < radialCenters.size(); ++bin) {
                    if (radialCenters[bin]-0.5*deltaK < Kh && Kh <= radialCenters[bin]+0.5*deltaK) {
                        radialBin[modeIndex] = bin;
                        ++radialCounts[bin];
                        break;
                    }
                }
                if (radialBin[modeIndex] == radialCenters.size()) return {WVKernelStatusCode::invalidConfiguration,"A horizontal mode was not assigned to an adaptive radial tolerance bin."};
            }

            const double f = verticalModes.coriolisFrequency;
            const double f2 = f*f;
            const double N02 = configuration.N0*configuration.N0;
            for (std::size_t modeIndex = 0; modeIndex < horizontalModes.size(); ++modeIndex) {
                const double Kh = horizontalModes[modeIndex].Kh;
                const double Kh2 = Kh*Kh;
                const auto bin = radialBin[modeIndex];
                const double center = radialCenters[bin];
                const double radialEnergy = center+0.5*deltaK-std::max(center-0.5*deltaK,0.0);
                const double energyPerCoefficient = radialEnergy/static_cast<double>(radialCounts[bin]);
                for (std::size_t jIndex = 0; jIndex < configuration.Nj; ++jIndex) {
                    const auto coefficientIndex = jIndex+configuration.Nj*modeIndex;
                    const double verticalWavenumber = verticalModes.verticalWavenumber[jIndex];
                    double waveEnergyFactor = 0.0;
                    if (Kh == 0.0) {
                        waveEnergyFactor = jIndex == 0 ? configuration.Lz :
                            (configuration.isHydrostatic ? N02/(configuration.g*verticalWavenumber*verticalWavenumber) :
                                (N02-f2)/(configuration.g*verticalWavenumber*verticalWavenumber));
                    } else if (jIndex > 0) {
                        const double hpm = configuration.isHydrostatic ? N02/(configuration.g*verticalWavenumber*verticalWavenumber) :
                            (N02-f2)/(configuration.g*(verticalWavenumber*verticalWavenumber+Kh2));
                        waveEnergyFactor = 2.0*hpm;
                    }
                    double vortexEnergyFactor = 0.0;
                    if (Kh > 0.0) {
                        const double h0 = verticalModes.h0[jIndex];
                        const double deformationInverse = jIndex == 0 ? 0.0 : f2/(configuration.g*h0);
                        vortexEnergyFactor = h0/(Kh2+deformationInverse);
                    } else if (jIndex > 0) {
                        vortexEnergyFactor = configuration.g/2.0;
                    }
                    if (waveEnergyFactor > 0.0 && std::isfinite(waveEnergyFactor)) {
                        policy->waveTolerance_[coefficientIndex] = absoluteToleranceScale*std::sqrt(energyPerCoefficient/waveEnergyFactor);
                    }
                    if (vortexEnergyFactor > 0.0 && std::isfinite(vortexEnergyFactor)) {
                        policy->vortexTolerance_[coefficientIndex] = absoluteToleranceScale*std::sqrt(energyPerCoefficient/vortexEnergyFactor);
                    }
                }
            }
            result = std::move(policy);
            return WVKernelStatus::ok();
        } catch (const std::bad_alloc&) {
            return {WVKernelStatusCode::allocationFailure,"Adaptive tolerance allocation failed."};
        }
    }

    std::size_t componentCount() const noexcept override { return 3; }
    std::size_t elementCount(std::size_t component) const noexcept override { return component < 3 ? shape_.elementCount() : 0; }
    double absoluteTolerance(std::size_t component, std::size_t index) const noexcept override {
        if (component > 2 || index >= shape_.elementCount()) return std::numeric_limits<double>::quiet_NaN();
        return component < 2 ? waveTolerance_[index] : vortexTolerance_[index];
    }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this) + vectorBytes(waveTolerance_) +
               vectorBytes(vortexTolerance_);
    }

private:
    WVShape2D shape_;
    std::vector<double> waveTolerance_;
    std::vector<double> vortexTolerance_;
};

} // namespace

namespace {

const std::vector<double>* realValues(const WVPortableTypedRecord& record, const char* name) {
    const auto* value = record.value(name);
    return value == nullptr ? nullptr : std::get_if<std::vector<double>>(&value->storage);
}

const std::vector<std::int64_t>* integerValues(const WVPortableTypedRecord& record, const char* name) {
    const auto* value = record.value(name);
    return value == nullptr ? nullptr : std::get_if<std::vector<std::int64_t>>(&value->storage);
}

const std::vector<std::uint8_t>* booleanValues(const WVPortableTypedRecord& record, const char* name) {
    const auto* value = record.value(name);
    return value == nullptr ? nullptr : std::get_if<std::vector<std::uint8_t>>(&value->storage);
}

const std::vector<std::string>* textValues(const WVPortableTypedRecord& record, const char* name) {
    const auto* value = record.value(name);
    return value == nullptr ? nullptr : std::get_if<std::vector<std::string>>(&value->storage);
}

class ResolvedForcing : public WVForcing {
public:
    explicit ResolvedForcing(const WVFrozenForcingEntry& entry)
        : typeIdentifier_(entry.typeIdentifier), name_(entry.name),
          contractVersion_(entry.contractVersion), stage_(entry.stage),
          priority_(entry.priority), ordinal_(entry.ordinal) {}
    const std::string& typeIdentifier() const noexcept override { return typeIdentifier_; }
    std::uint32_t contractVersion() const noexcept override { return contractVersion_; }
    const std::string& name() const noexcept override { return name_; }
    WVForcingStage stage() const noexcept override { return stage_; }
    std::uint8_t priority() const noexcept override { return priority_; }
    std::size_t ordinal() const noexcept override { return ordinal_; }
    bool supportsTendencyDiagnostics() const noexcept override { return true; }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this)+metadataDynamicBytes();
    }
protected:
    std::size_t metadataDynamicBytes() const noexcept {
        return typeIdentifier_.capacity()+name_.capacity();
    }
    std::string typeIdentifier_;
    std::string name_;
    std::uint32_t contractVersion_ = 0;
    WVForcingStage stage_ = WVForcingStage::spatial;
    std::uint8_t priority_ = 255;
    std::size_t ordinal_ = 0;
};

class LaplacianDampingForcing final : public ResolvedForcing {
public:
    LaplacianDampingForcing(const WVFrozenForcingEntry& entry, double nu, double kappa, WVLaplacianDirection direction)
        : ResolvedForcing(entry), nu_(nu), kappa_(kappa), direction_(direction) {}
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override {
        return context.laplacianDamping(nu_,kappa_,direction_);
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes(); }
private:
    double nu_,kappa_;
    WVLaplacianDirection direction_;
};

class ExplicitAntialiasingForcing final : public ResolvedForcing {
public:
    ExplicitAntialiasingForcing(const WVFrozenForcingEntry& entry, std::vector<std::size_t> indices)
        : ResolvedForcing(entry), indices_(std::move(indices)) {}
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override {
        context.filterTendency(indices_);
        return WVKernelStatus::ok();
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes()+vectorBytes(indices_); }
private:
    std::vector<std::size_t> indices_;
};

class NonlinearAdvectionForcing final : public ResolvedForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
    using ResolvedForcing::ResolvedForcing;
    bool producesCompleteFlux() const noexcept override { return true; }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { return context.nonlinearAdvection(); }
};

class QuadraticBottomFrictionForcing final : public ResolvedForcing {
public:
    QuadraticBottomFrictionForcing(WVFrozenForcingEntry entry, double value) : ResolvedForcing(std::move(entry)), drag_(value) {}
    bool requiresPhysicalFields() const noexcept override { return true; }
    bool requiresForcingFields() const noexcept override { return true; }
    bool producesCompleteFlux() const noexcept override { return true; }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override {
        WVRealFieldBundleConstView fields;
        auto status = context.physicalFields(fields);
        if (!status) return status;
        auto tendency = context.clearedSpatialTendency();
        const auto R = fields.shape.first*fields.shape.second*fields.shape.third;
        const auto horizontalCount = fields.shape.first*fields.shape.second;
        for (std::size_t index = 0; index < horizontalCount; ++index) {
            const double u = fields.data[index];
            const double v = fields.data[R+index];
            const double speed = std::sqrt(u*u+v*v);
            tendency.data[index] = -drag_*u*speed;
            tendency.data[R+index] = -drag_*v*speed;
        }
        return context.projectSpatialTendency({tendency.data,tendency.shape});
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes(); }
private:
    double drag_ = 0.0;
};

class LinearBottomFrictionForcing final : public ResolvedForcing {
public:
    LinearBottomFrictionForcing(WVFrozenForcingEntry entry, double value) : ResolvedForcing(std::move(entry)), drag_(value) {}
    bool requiresPhysicalFields() const noexcept override { return true; }
    bool requiresForcingFields() const noexcept override { return true; }
    bool producesCompleteFlux() const noexcept override { return true; }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override {
        WVRealFieldBundleConstView fields;
        auto status = context.physicalFields(fields);
        if (!status) return status;
        auto tendency = context.clearedSpatialTendency();
        const auto R = fields.shape.first*fields.shape.second*fields.shape.third;
        const auto horizontalCount = fields.shape.first*fields.shape.second;
        for (std::size_t index = 0; index < horizontalCount; ++index) {
            tendency.data[index] = -drag_*fields.data[index];
            tendency.data[R+index] = -drag_*fields.data[R+index];
        }
        return context.projectSpatialTendency({tendency.data,tendency.shape});
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes(); }
private:
    double drag_ = 0.0;
};

class AdaptiveDampingForcing final : public ResolvedForcing {
public:
    AdaptiveDampingForcing(WVFrozenForcingEntry entry, std::vector<double> values) : ResolvedForcing(std::move(entry)), damping_(std::move(values)) {}
    bool requiresPhysicalFields() const noexcept override { return true; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes()+vectorBytes(damping_); }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { return context.adaptiveDamping(damping_); }
private:
    std::vector<double> damping_;
};

class BetaPlaneForcing final : public ResolvedForcing {
public:
    BetaPlaneForcing(WVFrozenForcingEntry entry, std::vector<WVComplex64> values) : ResolvedForcing(std::move(entry)), betaA0_(std::move(values)) {}
    std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes()+vectorBytes(betaA0_); }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { return context.betaPlaneAdvection(betaA0_); }
private:
    std::vector<WVComplex64> betaA0_;
};

class PseudoTopographicForcing final : public ResolvedForcing {
public:
    PseudoTopographicForcing(WVFrozenForcingEntry entry, WVPseudoTopographicOperators values) : ResolvedForcing(std::move(entry)), operators_(std::move(values)) {}
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this)+metadataDynamicBytes()+vectorBytes(operators_.configuration.topographicHeight)+operators_.configuration.darwinSymbol.capacity()+vectorBytes(operators_.responsePlusX)+vectorBytes(operators_.responsePlusY)+vectorBytes(operators_.responseMinusX)+vectorBytes(operators_.responseMinusY);
    }
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { return context.pseudoTopographicGeneration(operators_); }
private:
    WVPseudoTopographicOperators operators_;
};

class FixedAmplitudeForcing final : public ResolvedForcing {
public:
    FixedAmplitudeForcing(WVFrozenForcingEntry entry, WVFixedAmplitudeConfiguration values) : ResolvedForcing(std::move(entry)), values_(std::move(values)) {}
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { context.zeroSelectedTendencies(values_); return WVKernelStatus::ok(); }
    WVStateConstraintResult applyConstraint(WVMutableCoefficients& coefficients) const override {
        std::size_t modified = 0;
        const auto restore = [&modified](WVComplexView destination, const auto& indices, const auto& values) {
            for (std::size_t index = 0; index < indices.size(); ++index) {
                const auto previous = destination.data[indices[index]];
                if (previous.real != values[index].real || previous.imag != values[index].imag) ++modified;
                destination.data[indices[index]] = values[index];
            }
        };
        restore(coefficients.Ap,values_.ApIndices,values_.ApValues);
        restore(coefficients.Am,values_.AmIndices,values_.AmValues);
        restore(coefficients.A0,values_.A0Indices,values_.A0Values);
        return {WVKernelStatus::ok(),modified,false};
    }
    std::size_t constraintWriteCount() const noexcept override { return values_.ApIndices.size()+values_.AmIndices.size()+values_.A0Indices.size(); }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this)+metadataDynamicBytes()+vectorBytes(values_.ApIndices)+vectorBytes(values_.ApValues)+vectorBytes(values_.AmIndices)+vectorBytes(values_.AmValues)+vectorBytes(values_.A0Indices)+vectorBytes(values_.A0Values);
    }
private:
    WVFixedAmplitudeConfiguration values_;
};

std::vector<double> adaptiveDampingOperator(const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation& preparation, WVKernelStatus& status) {
    const auto& configuration = descriptor.configuration();
    const auto& modes = descriptor.verticalModes();
    const double f2 = modes.coriolisFrequency*modes.coriolisFrequency;
    if (!(f2 > 0.0)) { status = {WVKernelStatusCode::unsupportedOperation,"Adaptive damping requires nonzero Coriolis frequency."}; return {}; }
    const double maximumComponent = preparation.maximumHorizontalComponent;
    if (!(maximumComponent > 0.0)) { status = {WVKernelStatusCode::invalidConfiguration,"Adaptive damping requires resolved horizontal wavenumbers."}; return {}; }
    const double effectiveResolution = pi/maximumComponent;
    const double dklMinimum = std::min(2.0*pi/configuration.Lx,2.0*pi/configuration.Ly);
    const double klCutoff = dklMinimum*std::pow(maximumComponent/dklMinimum,0.75);
    const double jMaximum = static_cast<double>(preparation.maximumVerticalMode);
    const double jCutoff = configuration.Nj > 2 ? std::pow(jMaximum,0.75) : 0.0;
    const double prefactorXY = effectiveResolution/(pi*pi);
    const double referenceLr2 = configuration.g*modes.h0[preparation.maximumVerticalMode]/f2;
    const double prefactorZ = (pi*pi*referenceLr2/(effectiveResolution*effectiveResolution))*prefactorXY;
    std::vector<double> damping(descriptor.spectralShape().elementCount());
    for (std::size_t iMode = 0; iMode < descriptor.Nkl(); ++iMode) {
        const auto& horizontal = descriptor.fourierModes()[iMode];
        const double Qkl = vanishingFilter(horizontal.Kh,klCutoff,maximumComponent);
        for (std::size_t iJ = 0; iJ < configuration.Nj; ++iJ) {
            const auto index = iJ+configuration.Nj*iMode;
            const double Qj = configuration.Nj > 2 ? vanishingFilter(modes.j[iJ],jCutoff,jMaximum) : 1.0;
            const double Lr2Inverse = f2/(configuration.g*modes.h0[iJ]);
            damping[index] = -prefactorXY*Qkl*horizontal.Kh*horizontal.Kh-prefactorZ*Qj*Lr2Inverse;
        }
    }
    status = WVKernelStatus::ok();
    return damping;
}


std::vector<double> adaptiveDampingOperator(const WVTransformHydrostaticKernel& kernel,const WVForcingPreparation& preparation,WVKernelStatus& status) {
    const auto& g=kernel.geometry();
    const double f=2*g.rotationRate*std::sin(g.latitude*pi/180),maximum=preparation.maximumHorizontalComponent;
    if (!(maximum>0) || f==0) { status={WVKernelStatusCode::invalidConfiguration,"Hydrostatic damping requires nonzero rotation and resolved horizontal modes."}; return {}; }
    const double dx=pi/maximum,dk=std::min(2*pi/g.Lx,2*pi/g.Ly),cutoff=dk*std::pow(maximum/dk,.75);
    const double jm=static_cast<double>(preparation.maximumVerticalMode),dj=g.Nj>1 ? g.j[1]-g.j[0] : 1;
    const double jc=g.Nj>2 ? dj*std::pow(jm/dj,.75) : 0;
    const auto found=std::find(g.j.begin(),g.j.end(),jm);
    if(found==g.j.end()) { status={WVKernelStatusCode::invalidConfiguration,"Effective damping mode is not retained."}; return {}; }
    const auto jindex=static_cast<std::size_t>(found-g.j.begin());
    const double xy=dx/(pi*pi),zz=(pi*pi*g.g*g.h_0[jindex]/(f*f*dx*dx))*xy;
    std::vector<double> result(g.Nj*g.Nkl);
    for (std::size_t i=0;i<result.size();++i) {
        const auto j=i%g.Nj,m=i/g.Nj; const double kh=std::hypot(g.k[m],g.l[m]);
        result[i]=-xy*vanishingFilter(kh,cutoff,maximum)*kh*kh-zz*(g.Nj>2 ? vanishingFilter(g.j[j],jc,jm) : 1)*f*f/(g.g*g.h_0[j]);
    }
    status=WVKernelStatus::ok(); return result;
}
std::vector<double> adaptiveDampingOperator(const WVTransformBoussinesqKernel& kernel,const WVForcingPreparation& preparation,WVKernelStatus& status) {
    const auto& g=kernel.geometry();
    const double f=2*g.rotationRate*std::sin(g.latitude*pi/180),maximum=preparation.maximumHorizontalComponent;
    if (!(maximum>0) || f==0) { status={WVKernelStatusCode::invalidConfiguration,"Boussinesq damping requires nonzero rotation and resolved horizontal modes."}; return {}; }
    const double dx=pi/maximum,dk=std::min(2*pi/g.Lx,2*pi/g.Ly),cutoff=dk*std::pow(maximum/dk,.75);
    const double jm=static_cast<double>(preparation.maximumVerticalMode),dj=g.Nj>1 ? g.j[1]-g.j[0] : 1;
    const double jc=g.Nj>2 ? dj*std::pow(jm/dj,.75) : 0;
    const auto found=std::find(g.j.begin(),g.j.end(),jm);
    if(found==g.j.end()) { status={WVKernelStatusCode::invalidConfiguration,"Effective damping mode is not retained."}; return {}; }
    const auto jindex=static_cast<std::size_t>(found-g.j.begin());
    const double xy=dx/(pi*pi),zz=(pi*pi*g.g*g.h_0[jindex]/(f*f*dx*dx))*xy;
    std::vector<double> result(g.Nj*g.Nkl);
    for (std::size_t i=0;i<result.size();++i) {
        const auto j=i%g.Nj,m=i/g.Nj; const double kh=std::hypot(g.k[m],g.l[m]);
        result[i]=-xy*vanishingFilter(kh,cutoff,maximum)*kh*kh-zz*(g.Nj>2 ? vanishingFilter(g.j[j],jc,jm) : 1)*f*f/(g.g*g.h_0[j]);
    }
    status=WVKernelStatus::ok(); return result;
}
class HydrostaticVerticalDiffusivityForcing final : public ResolvedForcing {
public:
    HydrostaticVerticalDiffusivityForcing(const WVFrozenForcingEntry& e,double k,bool mean):ResolvedForcing(e),kappa_(k),mean_(mean) {}
    WVKernelStatus addRightHandSide(WVForcingExecutionContext& context) const override { return context.verticalDiffusivity(kappa_,mean_); }
private:
    double kappa_; bool mean_;
};

bool emptyConfiguration(const WVFrozenForcingEntry& entry) { return entry.configuration.values.empty(); }

} // namespace

namespace detail {

WVKernelStatus preflightLaplacianDamping(const WVFrozenForcingEntry& entry, bool) {
    for (const auto* name : {"nu","kappa"}) {
        const auto* value = realValues(entry.configuration,name);
        if (value == nullptr || value->size() != 1 || !std::isfinite(value->front()) || value->front() < 0.0)
            return {WVKernelStatusCode::invalidConfiguration,"Damping nu and kappa must be finite nonnegative scalars."};
    }
    return WVKernelStatus::ok();
}
WVKernelStatus preflightVerticalDiffusivity(const WVFrozenForcingEntry& entry, bool) {
    const auto* value = realValues(entry.configuration,"kappa_z");
    const auto* flag = booleanValues(entry.configuration,"shouldForceMeanDensityAnomaly");
    if (value == nullptr || value->size() != 1 || !std::isfinite(value->front()) || value->front() < 0.0 || flag == nullptr || flag->size() != 1 || flag->front() > 1)
        return {WVKernelStatusCode::invalidConfiguration,"Vertical diffusivity requires scalar nonnegative kappa_z and a scalar logical mean-density flag."};
    return WVKernelStatus::ok();
}
WVKernelStatus createHorizontalDamping(const WVFrozenForcingEntry& entry, const WVTransformConstantStratificationDescriptor&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::horizontal);
    return WVKernelStatus::ok();
}
WVKernelStatus createVerticalDamping(const WVFrozenForcingEntry& entry, const WVTransformConstantStratificationDescriptor&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::vertical);
    return WVKernelStatus::ok();
}
WVKernelStatus createVerticalDiffusivity(const WVFrozenForcingEntry& entry, const WVTransformConstantStratificationDescriptor&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightVerticalDiffusivity(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,0.0,realValues(entry.configuration,"kappa_z")->front(),WVLaplacianDirection::vertical);
    return WVKernelStatus::ok();
}

WVKernelStatus preflightExplicitAntialiasing(const WVFrozenForcingEntry& entry, bool shouldAntialias) {
    if (shouldAntialias) return {WVKernelStatusCode::invalidConfiguration,"Explicit antialiasing requires shouldAntialias=false."};
    const auto* values = realValues(entry.configuration,"Nj");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0 || std::floor(values->front()) != values->front())
        return {WVKernelStatusCode::invalidConfiguration,"Explicit antialiasing Nj must be one finite nonnegative integer."};
    return WVKernelStatus::ok();
}

WVKernelStatus prepareExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, WVForcingPreparation& preparation) {
    // MATLAB's transform resolution accessor recognizes this registered name.
    if (entry.name != "antialias filter") return WVKernelStatus::ok();
    const auto& c = descriptor.configuration();
    const double retained = realValues(entry.configuration,"Nj")->front();
    const double cutoff = (2.0/3.0)*2.0*pi*static_cast<double>(c.Nx/2)/c.Lx;
    preparation.maximumHorizontalComponent = 0.0;
    preparation.maximumVerticalMode = 0;
    for (const auto& mode : descriptor.fourierModes()) {
        if (mode.Kh > cutoff) continue;
        for (std::size_t j = 0; j < c.Nj; ++j) {
            if (static_cast<double>(j) > retained-1.0) continue;
            preparation.maximumHorizontalComponent = std::max({preparation.maximumHorizontalComponent,std::abs(mode.k),std::abs(mode.l)});
            preparation.maximumVerticalMode = std::max(preparation.maximumVerticalMode,j);
        }
    }
    return WVKernelStatus::ok();
}

WVKernelStatus createExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto& c = descriptor.configuration();
    auto status = preflightExplicitAntialiasing(entry,c.shouldAntialias);
    if (!status) return status;
    const double Nj = realValues(entry.configuration,"Nj")->front();
    const double maximumK = 2.0*pi*static_cast<double>(c.Nx/2)/c.Lx;
    std::vector<std::size_t> indices;
    for (std::size_t kl = 0; kl < descriptor.Nkl(); ++kl)
        for (std::size_t j = 0; j < c.Nj; ++j)
            if (descriptor.fourierModes()[kl].Kh > 2.0*maximumK/3.0 || descriptor.verticalModes().j[j] > Nj-1.0)
                indices.push_back(j+c.Nj*kl);
    forcing = std::make_unique<ExplicitAntialiasingForcing>(entry,std::move(indices));
    return WVKernelStatus::ok();
}

WVKernelStatus createNonlinearAdvectionForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Nonlinear advection does not accept configuration values."};
    forcing = std::make_unique<NonlinearAdvectionForcing>(entry);
    return WVKernelStatus::ok();
}

WVKernelStatus createAdaptiveDampingForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &preparation,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Adaptive damping does not accept configuration values."};
    WVKernelStatus status;
    auto damping = adaptiveDampingOperator(descriptor,preparation,status);
    if (!status) return status;
    forcing = std::make_unique<AdaptiveDampingForcing>(entry,std::move(damping));
    return WVKernelStatus::ok();
}

WVKernelStatus createQuadraticBottomFriction(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"Cd");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Quadratic drag coefficient must be one finite nonnegative scalar."};
    forcing = std::make_unique<QuadraticBottomFrictionForcing>(entry,values->front()/descriptor.bottomQuadratureWeight());
    return WVKernelStatus::ok();
}

WVKernelStatus createLinearBottomFriction(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"r");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Linear drag rate must be one finite nonnegative scalar."};
    const auto Nz = descriptor.configuration().Nz;
    const double scaledRate = 2.0*static_cast<double>(Nz-1)*values->front();
    forcing = std::make_unique<LinearBottomFrictionForcing>(entry,scaledRate);
    return WVKernelStatus::ok();
}

WVKernelStatus createFixedAmplitudeForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    WVFixedAmplitudeConfiguration configuration;
    const auto decode = [&](const char* indexName, const char* realName, const char* imagName, auto& indices, auto& values) -> bool {
        const auto* sourceIndices = integerValues(entry.configuration,indexName);
        const auto* real = realValues(entry.configuration,realName);
        const auto* imag = realValues(entry.configuration,imagName);
        if (sourceIndices == nullptr && real == nullptr && imag == nullptr) return true;
        if (sourceIndices == nullptr || real == nullptr || imag == nullptr || sourceIndices->size() != real->size() || real->size() != imag->size()) return false;
        std::set<std::size_t> unique;
        for (std::size_t index = 0; index < sourceIndices->size(); ++index) {
            if ((*sourceIndices)[index] < 0 || static_cast<std::size_t>((*sourceIndices)[index]) >= descriptor.spectralShape().elementCount() || !std::isfinite((*real)[index]) || !std::isfinite((*imag)[index])) return false;
            const auto converted = static_cast<std::size_t>((*sourceIndices)[index]);
            if (!unique.insert(converted).second) return false;
            indices.push_back(converted);
            values.push_back({(*real)[index],(*imag)[index]});
        }
        return true;
    };
    if (!decode("ApIndices","ApValuesReal","ApValuesImag",configuration.ApIndices,configuration.ApValues) ||
        !decode("AmIndices","AmValuesReal","AmValuesImag",configuration.AmIndices,configuration.AmValues) ||
        !decode("A0Indices","A0ValuesReal","A0ValuesImag",configuration.A0Indices,configuration.A0Values)) return {WVKernelStatusCode::invalidConfiguration,"Fixed-amplitude configuration is malformed or outside the coefficient shape."};
    forcing = std::make_unique<FixedAmplitudeForcing>(entry,std::move(configuration));
    return WVKernelStatus::ok();
}

WVKernelStatus createBetaPlaneForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Beta-plane advection does not accept configuration values."};
    const auto& configuration = descriptor.configuration();
    const auto& modes = descriptor.verticalModes();
    const double beta = 2.0*configuration.rotationRate*std::cos(configuration.latitude*pi/180.0)/configuration.planetaryRadius;
    std::vector<WVComplex64> betaA0(descriptor.spectralShape().elementCount());
    for (std::size_t iMode = 0; iMode < descriptor.Nkl(); ++iMode)
        for (std::size_t iJ = 0; iJ < configuration.Nj; ++iJ) {
            const auto index = iJ+configuration.Nj*iMode;
            betaA0[index] = iMode == 0 ? WVComplex64{} : multiply(modes.VA0Field[index],-beta/modes.Fg[iJ]);
        }
    forcing = std::make_unique<BetaPlaneForcing>(entry,std::move(betaA0));
    return WVKernelStatus::ok();
}

WVKernelStatus createPseudoTopographicForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor,
    const WVForcingPreparation &preparation, std::unique_ptr<WVForcing>& forcing) {
    WVPseudoTopographicOperators operators;
    auto& record = operators.configuration;
    const auto& configuration = descriptor.configuration();
    const auto& modes = descriptor.verticalModes();
    const auto* height = realValues(entry.configuration,"topographicHeight");
    const auto* velocityReal = realValues(entry.configuration,"barotropicVelocityAmplitudeReal");
    const auto* velocityImag = realValues(entry.configuration,"barotropicVelocityAmplitudeImag");
    const auto* frequency = realValues(entry.configuration,"frequency");
    const auto* ramp = realValues(entry.configuration,"rampDuration");
    const auto* start = realValues(entry.configuration,"startTime");
    const auto* avoid = booleanValues(entry.configuration,"shouldAvoidAdaptiveDamping");
    const auto* maximumK = realValues(entry.configuration,"maximumForcedHorizontalWavenumber");
    const auto* maximumJ = realValues(entry.configuration,"maximumForcedVerticalMode");
    const auto* symbol = textValues(entry.configuration,"darwinSymbol");
    if (height == nullptr || height->size() != configuration.Nx*configuration.Ny || velocityReal == nullptr || velocityImag == nullptr || velocityReal->size() != 2 || velocityImag->size() != 2 || frequency == nullptr || frequency->size() != 1 || ramp == nullptr || ramp->size() != 1 || start == nullptr || start->size() != 1 || avoid == nullptr || avoid->size() != 1 || maximumK == nullptr || maximumK->size() != 1 || maximumJ == nullptr || maximumJ->size() != 1) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing configuration is incomplete."};
    record.topographicShape = {configuration.Nx,configuration.Ny};
    record.topographicHeight = *height;
    record.barotropicVelocityAmplitude = {{{(*velocityReal)[0],(*velocityImag)[0]},{(*velocityReal)[1],(*velocityImag)[1]}}};
    record.frequency = frequency->front(); record.rampDuration = ramp->front(); record.startTime = start->front();
    record.shouldAvoidAdaptiveDamping = avoid->front() != 0;
    record.maximumForcedHorizontalWavenumber = maximumK->front(); record.maximumForcedVerticalMode = maximumJ->front();
    record.darwinSymbol = symbol == nullptr || symbol->empty() ? std::string{} : symbol->front();
    static const std::set<std::string> validSymbols = {"","M2","S2","N2","K1","O1"};
    if (!std::isfinite(record.frequency) || record.frequency <= 0.0 || !std::isfinite(record.rampDuration) || record.rampDuration < 0.0 || !std::isfinite(record.startTime) || std::isnan(record.maximumForcedHorizontalWavenumber) || record.maximumForcedHorizontalWavenumber < 0.0 || std::isnan(record.maximumForcedVerticalMode) || record.maximumForcedVerticalMode < 0.0 || validSymbols.count(record.darwinSymbol) == 0) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing values are invalid."};
    WVKernelStatus dampingStatus = WVKernelStatus::ok();
    const auto damping = preparation.hasAdaptiveDamping && record.shouldAvoidAdaptiveDamping ? adaptiveDampingOperator(descriptor,preparation,dampingStatus) : std::vector<double>{};
    if (!dampingStatus) return dampingStatus;
    const auto count = descriptor.spectralShape().elementCount();
    operators.responsePlusX.assign(count,{}); operators.responsePlusY.assign(count,{}); operators.responseMinusX.assign(count,{}); operators.responseMinusY.assign(count,{});
    for (std::size_t iMode = 0; iMode < descriptor.Nkl(); ++iMode) {
        const auto& horizontal = descriptor.fourierModes()[iMode];
        const auto terrain = normalizedTerrainCoefficient(record,horizontal.kMode,horizontal.lMode);
        const auto dHdx = multiply(WVComplex64{-terrain.imag,terrain.real},horizontal.k);
        const auto dHdy = multiply(WVComplex64{-terrain.imag,terrain.real},horizontal.l);
        for (std::size_t iJ = 1; iJ < configuration.Nj; ++iJ) {
            const auto index = iJ+configuration.Nj*iMode;
            if (!(horizontal.Kh > 0.0) || horizontal.Kh > record.maximumForcedHorizontalWavenumber || modes.j[iJ] > record.maximumForcedVerticalMode || (!damping.empty() && damping[index] != 0.0)) continue;
            const double NAp = modes.NApField[index]/modes.gWaveScale[iJ];
            const double hpm = -NAp*modes.omega[index]/horizontal.Kh;
            if (!(2.0*hpm > 0.0)) continue;
            const WVComplex64 piPlus{configuration.g*modes.fWaveScale[index]*NAp,0.0};
            const WVComplex64 piMinus{-piPlus.real,0.0};
            operators.responsePlusX[index] = multiply(multiply(conjugate(piPlus),dHdx),1.0/(2.0*hpm));
            operators.responsePlusY[index] = multiply(multiply(conjugate(piPlus),dHdy),1.0/(2.0*hpm));
            operators.responseMinusX[index] = multiply(multiply(conjugate(piMinus),dHdx),1.0/(2.0*hpm));
            operators.responseMinusY[index] = multiply(multiply(conjugate(piMinus),dHdy),1.0/(2.0*hpm));
        }
    }
    forcing = std::make_unique<PseudoTopographicForcing>(entry,std::move(operators));
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticHorizontalDamping(const WVFrozenForcingEntry& entry, WVTransformHydrostaticKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::horizontal);
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticVerticalDamping(const WVFrozenForcingEntry& entry, WVTransformHydrostaticKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::vertical);
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticVerticalDiffusivity(const WVFrozenForcingEntry& entry, WVTransformHydrostaticKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightVerticalDiffusivity(entry,false); if (!status) return status;
    forcing = std::make_unique<HydrostaticVerticalDiffusivityForcing>(entry,realValues(entry.configuration,"kappa_z")->front(),booleanValues(entry.configuration,"shouldForceMeanDensityAnomaly")->front()!=0);
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticNonlinearAdvectionForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Nonlinear advection does not accept configuration values."};
    forcing = std::make_unique<NonlinearAdvectionForcing>(entry);
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticAdaptiveDampingForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &preparation,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Adaptive damping does not accept configuration values."};
    WVKernelStatus status;
    auto damping = adaptiveDampingOperator(descriptor,preparation,status);
    if (!status) return status;
    forcing = std::make_unique<AdaptiveDampingForcing>(entry,std::move(damping));
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticQuadraticBottomFriction(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"Cd");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Quadratic drag coefficient must be one finite nonnegative scalar."};
    forcing = std::make_unique<QuadraticBottomFrictionForcing>(entry,values->front()/descriptor.geometry().z_int.front());
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticLinearBottomFriction(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"r");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Linear drag rate must be one finite nonnegative scalar."};
    const double scaledRate = descriptor.geometry().Lz*values->front()/descriptor.geometry().z_int.front();
    forcing = std::make_unique<LinearBottomFrictionForcing>(entry,scaledRate);
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticFixedAmplitudeForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    WVFixedAmplitudeConfiguration configuration;
    const auto decode = [&](const char* indexName, const char* realName, const char* imagName, auto& indices, auto& values) -> bool {
        const auto* sourceIndices = integerValues(entry.configuration,indexName);
        const auto* real = realValues(entry.configuration,realName);
        const auto* imag = realValues(entry.configuration,imagName);
        if (sourceIndices == nullptr && real == nullptr && imag == nullptr) return true;
        if (sourceIndices == nullptr || real == nullptr || imag == nullptr || sourceIndices->size() != real->size() || real->size() != imag->size()) return false;
        std::set<std::size_t> unique;
        for (std::size_t index = 0; index < sourceIndices->size(); ++index) {
            if ((*sourceIndices)[index] < 0 || static_cast<std::size_t>((*sourceIndices)[index]) >= descriptor.spectralShape().elementCount() || !std::isfinite((*real)[index]) || !std::isfinite((*imag)[index])) return false;
            const auto converted = static_cast<std::size_t>((*sourceIndices)[index]);
            if (!unique.insert(converted).second) return false;
            indices.push_back(converted);
            values.push_back({(*real)[index],(*imag)[index]});
        }
        return true;
    };
    if (!decode("ApIndices","ApValuesReal","ApValuesImag",configuration.ApIndices,configuration.ApValues) ||
        !decode("AmIndices","AmValuesReal","AmValuesImag",configuration.AmIndices,configuration.AmValues) ||
        !decode("A0Indices","A0ValuesReal","A0ValuesImag",configuration.A0Indices,configuration.A0Values)) return {WVKernelStatusCode::invalidConfiguration,"Fixed-amplitude configuration is malformed or outside the coefficient shape."};
    forcing = std::make_unique<FixedAmplitudeForcing>(entry,std::move(configuration));
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticBetaPlaneForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Beta-plane advection does not accept configuration values."};
    const auto& configuration = descriptor.geometry();
    const double beta=2*configuration.rotationRate*std::cos(configuration.latitude*pi/180)/configuration.planetaryRadius;
    std::vector<WVComplex64> betaA0(descriptor.spectralShape().elementCount());
    for (std::size_t i=0;i<betaA0.size();++i) if (descriptor.factors()[i].geostrophic) betaA0[i]=multiply(descriptor.factors()[i].VA0,-beta);
    forcing = std::make_unique<BetaPlaneForcing>(entry,std::move(betaA0));
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticPseudoTopographicForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor,
    const WVForcingPreparation &preparation, std::unique_ptr<WVForcing>& forcing) {
    WVPseudoTopographicOperators operators;
    auto& record = operators.configuration;
    const auto& configuration = descriptor.geometry();
    const auto* height = realValues(entry.configuration,"topographicHeight");
    const auto* velocityReal = realValues(entry.configuration,"barotropicVelocityAmplitudeReal");
    const auto* velocityImag = realValues(entry.configuration,"barotropicVelocityAmplitudeImag");
    const auto* frequency = realValues(entry.configuration,"frequency");
    const auto* ramp = realValues(entry.configuration,"rampDuration");
    const auto* start = realValues(entry.configuration,"startTime");
    const auto* avoid = booleanValues(entry.configuration,"shouldAvoidAdaptiveDamping");
    const auto* maximumK = realValues(entry.configuration,"maximumForcedHorizontalWavenumber");
    const auto* maximumJ = realValues(entry.configuration,"maximumForcedVerticalMode");
    const auto* symbol = textValues(entry.configuration,"darwinSymbol");
    if (height == nullptr || height->size() != configuration.Nx*configuration.Ny || velocityReal == nullptr || velocityImag == nullptr || velocityReal->size() != 2 || velocityImag->size() != 2 || frequency == nullptr || frequency->size() != 1 || ramp == nullptr || ramp->size() != 1 || start == nullptr || start->size() != 1 || avoid == nullptr || avoid->size() != 1 || maximumK == nullptr || maximumK->size() != 1 || maximumJ == nullptr || maximumJ->size() != 1) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing configuration is incomplete."};
    record.topographicShape = {configuration.Nx,configuration.Ny};
    record.topographicHeight = *height;
    record.barotropicVelocityAmplitude = {{{(*velocityReal)[0],(*velocityImag)[0]},{(*velocityReal)[1],(*velocityImag)[1]}}};
    record.frequency = frequency->front(); record.rampDuration = ramp->front(); record.startTime = start->front();
    record.shouldAvoidAdaptiveDamping = avoid->front() != 0;
    record.maximumForcedHorizontalWavenumber = maximumK->front(); record.maximumForcedVerticalMode = maximumJ->front();
    record.darwinSymbol = symbol == nullptr || symbol->empty() ? std::string{} : symbol->front();
    static const std::set<std::string> validSymbols = {"","M2","S2","N2","K1","O1"};
    if (!std::isfinite(record.frequency) || record.frequency <= 0.0 || !std::isfinite(record.rampDuration) || record.rampDuration < 0.0 || !std::isfinite(record.startTime) || std::isnan(record.maximumForcedHorizontalWavenumber) || record.maximumForcedHorizontalWavenumber < 0.0 || std::isnan(record.maximumForcedVerticalMode) || record.maximumForcedVerticalMode < 0.0 || validSymbols.count(record.darwinSymbol) == 0) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing values are invalid."};
    WVKernelStatus dampingStatus = WVKernelStatus::ok();
    const auto damping = preparation.hasAdaptiveDamping && record.shouldAvoidAdaptiveDamping ? adaptiveDampingOperator(descriptor,preparation,dampingStatus) : std::vector<double>{};
    if (!dampingStatus) return dampingStatus;
    const auto count = descriptor.spectralShape().elementCount();
    operators.responsePlusX.assign(count,{}); operators.responsePlusY.assign(count,{}); operators.responseMinusX.assign(count,{}); operators.responseMinusY.assign(count,{});
    std::vector<double> bottom(configuration.Nj),physical(configuration.Nx*configuration.Ny*configuration.Nz);
    std::vector<WVComplex64> modal(count);
    for (std::size_t j=0;j<configuration.Nj;++j) {
        std::fill(modal.begin(),modal.end(),WVComplex64{}); modal[j]={1,0};
        const auto status=descriptor.transformToSpatial({modal.data(),descriptor.spectralShape()},WVHydrostaticFamily::F,{physical.data(),descriptor.spatialShape()}); if (!status) return status;
        bottom[j]=physical[0];
    }
    for (std::size_t m=0;m<configuration.Nkl;++m) {
        const auto terrain=normalizedTerrainCoefficient(record,configuration.modes[m].k,configuration.modes[m].l);
        const auto dx=multiply(WVComplex64{-terrain.imag,terrain.real},configuration.k[m]),dy=multiply(WVComplex64{-terrain.imag,terrain.real},configuration.l[m]);
        for (std::size_t j=0;j<configuration.Nj;++j) {
            const auto i=j+configuration.Nj*m; const auto& factor=descriptor.factors()[i];
            if (!factor.wave || std::hypot(configuration.k[m],configuration.l[m])>record.maximumForcedHorizontalWavenumber || configuration.j[j]>record.maximumForcedVerticalMode || (!damping.empty() && damping[i]!=0)) continue;
            const double response=configuration.g*bottom[j]*factor.NAp/factor.waveEnergy;
            operators.responsePlusX[i]=multiply(dx,response); operators.responsePlusY[i]=multiply(dy,response);
            operators.responseMinusX[i]=multiply(dx,-response); operators.responseMinusY[i]=multiply(dy,-response);
        }
    }
    forcing = std::make_unique<PseudoTopographicForcing>(entry,std::move(operators));
    return WVKernelStatus::ok();
}

WVKernelStatus prepareHydrostaticExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    const WVTransformHydrostaticKernel& kernel, WVForcingPreparation& preparation) {
    if (entry.name!="antialias filter") return WVKernelStatus::ok();
    const auto& g=kernel.geometry();
    const double retained=realValues(entry.configuration,"Nj")->front();
    const double cutoff=(2.0/3.0)*2*pi*static_cast<double>(g.Nx/2)/g.Lx;
    preparation.maximumHorizontalComponent=0; preparation.maximumVerticalMode=0;
    for(std::size_t m=0;m<g.Nkl;++m) for(double j:g.j) {
        if(std::hypot(g.k[m],g.l[m])>cutoff || j>retained-1) continue;
        preparation.maximumHorizontalComponent=std::max({preparation.maximumHorizontalComponent,std::abs(g.k[m]),std::abs(g.l[m])});
        preparation.maximumVerticalMode=std::max(preparation.maximumVerticalMode,static_cast<std::size_t>(j));
    }
    return WVKernelStatus::ok();
}

WVKernelStatus createHydrostaticExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    WVTransformHydrostaticKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto& c=descriptor.geometry();
    const double retained=realValues(entry.configuration,"Nj")->front();
    const double cutoff=(2.0/3.0)*2*pi*static_cast<double>(c.Nx/2)/c.Lx;
    std::vector<std::size_t> indices;
    for (std::size_t i=0;i<c.Nj*c.Nkl;++i) if (std::hypot(c.k[i/c.Nj],c.l[i/c.Nj])>cutoff || c.j[i%c.Nj]>retained-1) indices.push_back(i);
    forcing = std::make_unique<ExplicitAntialiasingForcing>(entry,std::move(indices));
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqHorizontalDamping(const WVFrozenForcingEntry& entry, WVTransformBoussinesqKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::horizontal);
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqVerticalDamping(const WVFrozenForcingEntry& entry, WVTransformBoussinesqKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightLaplacianDamping(entry,false); if (!status) return status;
    forcing = std::make_unique<LaplacianDampingForcing>(entry,realValues(entry.configuration,"nu")->front(),realValues(entry.configuration,"kappa")->front(),WVLaplacianDirection::vertical);
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqVerticalDiffusivity(const WVFrozenForcingEntry& entry, WVTransformBoussinesqKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto status = preflightVerticalDiffusivity(entry,false); if (!status) return status;
    forcing = std::make_unique<HydrostaticVerticalDiffusivityForcing>(entry,realValues(entry.configuration,"kappa_z")->front(),booleanValues(entry.configuration,"shouldForceMeanDensityAnomaly")->front()!=0);
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqNonlinearAdvectionForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel&, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Nonlinear advection does not accept configuration values."};
    forcing = std::make_unique<NonlinearAdvectionForcing>(entry);
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqAdaptiveDampingForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &preparation,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Adaptive damping does not accept configuration values."};
    WVKernelStatus status;
    auto damping = adaptiveDampingOperator(descriptor,preparation,status);
    if (!status) return status;
    forcing = std::make_unique<AdaptiveDampingForcing>(entry,std::move(damping));
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqQuadraticBottomFriction(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"Cd");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Quadratic drag coefficient must be one finite nonnegative scalar."};
    forcing = std::make_unique<QuadraticBottomFrictionForcing>(entry,values->front()/descriptor.geometry().z_int.front());
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqLinearBottomFriction(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"r");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Linear drag rate must be one finite nonnegative scalar."};
    const double scaledRate = descriptor.geometry().Lz*values->front()/descriptor.geometry().z_int.front();
    forcing = std::make_unique<LinearBottomFrictionForcing>(entry,scaledRate);
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqFixedAmplitudeForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    WVFixedAmplitudeConfiguration configuration;
    const auto decode = [&](const char* indexName, const char* realName, const char* imagName, auto& indices, auto& values) -> bool {
        const auto* sourceIndices = integerValues(entry.configuration,indexName);
        const auto* real = realValues(entry.configuration,realName);
        const auto* imag = realValues(entry.configuration,imagName);
        if (sourceIndices == nullptr && real == nullptr && imag == nullptr) return true;
        if (sourceIndices == nullptr || real == nullptr || imag == nullptr || sourceIndices->size() != real->size() || real->size() != imag->size()) return false;
        std::set<std::size_t> unique;
        for (std::size_t index = 0; index < sourceIndices->size(); ++index) {
            if ((*sourceIndices)[index] < 0 || static_cast<std::size_t>((*sourceIndices)[index]) >= descriptor.spectralShape().elementCount() || !std::isfinite((*real)[index]) || !std::isfinite((*imag)[index])) return false;
            const auto converted = static_cast<std::size_t>((*sourceIndices)[index]);
            if (!unique.insert(converted).second) return false;
            indices.push_back(converted);
            values.push_back({(*real)[index],(*imag)[index]});
        }
        return true;
    };
    if (!decode("ApIndices","ApValuesReal","ApValuesImag",configuration.ApIndices,configuration.ApValues) ||
        !decode("AmIndices","AmValuesReal","AmValuesImag",configuration.AmIndices,configuration.AmValues) ||
        !decode("A0Indices","A0ValuesReal","A0ValuesImag",configuration.A0Indices,configuration.A0Values)) return {WVKernelStatusCode::invalidConfiguration,"Fixed-amplitude configuration is malformed or outside the coefficient shape."};
    forcing = std::make_unique<FixedAmplitudeForcing>(entry,std::move(configuration));
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqBetaPlaneForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Beta-plane advection does not accept configuration values."};
    const auto& configuration = descriptor.geometry();
    const double beta=2*configuration.rotationRate*std::cos(configuration.latitude*pi/180)/configuration.planetaryRadius;
    std::vector<WVComplex64> betaA0(descriptor.spectralShape().elementCount());
    for (std::size_t i=0;i<betaA0.size();++i) if (descriptor.factors()[i].geostrophic) betaA0[i]=multiply(descriptor.factors()[i].VA0,-beta);
    forcing = std::make_unique<BetaPlaneForcing>(entry,std::move(betaA0));
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqPseudoTopographicForcing(
    const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor,
    const WVForcingPreparation &preparation, std::unique_ptr<WVForcing>& forcing) {
    WVPseudoTopographicOperators operators;
    auto& record = operators.configuration;
    const auto& configuration = descriptor.geometry();
    const auto* height = realValues(entry.configuration,"topographicHeight");
    const auto* velocityReal = realValues(entry.configuration,"barotropicVelocityAmplitudeReal");
    const auto* velocityImag = realValues(entry.configuration,"barotropicVelocityAmplitudeImag");
    const auto* frequency = realValues(entry.configuration,"frequency");
    const auto* ramp = realValues(entry.configuration,"rampDuration");
    const auto* start = realValues(entry.configuration,"startTime");
    const auto* avoid = booleanValues(entry.configuration,"shouldAvoidAdaptiveDamping");
    const auto* maximumK = realValues(entry.configuration,"maximumForcedHorizontalWavenumber");
    const auto* maximumJ = realValues(entry.configuration,"maximumForcedVerticalMode");
    const auto* symbol = textValues(entry.configuration,"darwinSymbol");
    if (height == nullptr || height->size() != configuration.Nx*configuration.Ny || velocityReal == nullptr || velocityImag == nullptr || velocityReal->size() != 2 || velocityImag->size() != 2 || frequency == nullptr || frequency->size() != 1 || ramp == nullptr || ramp->size() != 1 || start == nullptr || start->size() != 1 || avoid == nullptr || avoid->size() != 1 || maximumK == nullptr || maximumK->size() != 1 || maximumJ == nullptr || maximumJ->size() != 1) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing configuration is incomplete."};
    record.topographicShape = {configuration.Nx,configuration.Ny};
    record.topographicHeight = *height;
    record.barotropicVelocityAmplitude = {{{(*velocityReal)[0],(*velocityImag)[0]},{(*velocityReal)[1],(*velocityImag)[1]}}};
    record.frequency = frequency->front(); record.rampDuration = ramp->front(); record.startTime = start->front();
    record.shouldAvoidAdaptiveDamping = avoid->front() != 0;
    record.maximumForcedHorizontalWavenumber = maximumK->front(); record.maximumForcedVerticalMode = maximumJ->front();
    record.darwinSymbol = symbol == nullptr || symbol->empty() ? std::string{} : symbol->front();
    static const std::set<std::string> validSymbols = {"","M2","S2","N2","K1","O1"};
    if (!std::isfinite(record.frequency) || record.frequency <= 0.0 || !std::isfinite(record.rampDuration) || record.rampDuration < 0.0 || !std::isfinite(record.startTime) || std::isnan(record.maximumForcedHorizontalWavenumber) || record.maximumForcedHorizontalWavenumber < 0.0 || std::isnan(record.maximumForcedVerticalMode) || record.maximumForcedVerticalMode < 0.0 || validSymbols.count(record.darwinSymbol) == 0) return {WVKernelStatusCode::invalidConfiguration,"Pseudo-topographic forcing values are invalid."};
    WVKernelStatus dampingStatus = WVKernelStatus::ok();
    const auto damping = preparation.hasAdaptiveDamping && record.shouldAvoidAdaptiveDamping ? adaptiveDampingOperator(descriptor,preparation,dampingStatus) : std::vector<double>{};
    if (!dampingStatus) return dampingStatus;
    const auto count = descriptor.spectralShape().elementCount();
    operators.responsePlusX.assign(count,{}); operators.responsePlusY.assign(count,{}); operators.responseMinusX.assign(count,{}); operators.responseMinusY.assign(count,{});
    std::vector<double> bottom(count);
    auto structure=descriptor.waveModeVerticalStructureAtIndex(0,{bottom.data(),descriptor.spectralShape()}); if (!structure) return structure;
    for (std::size_t m=0;m<configuration.Nkl;++m) {
        const auto terrain=normalizedTerrainCoefficient(record,configuration.modes[m].k,configuration.modes[m].l);
        const auto dx=multiply(WVComplex64{-terrain.imag,terrain.real},configuration.k[m]),dy=multiply(WVComplex64{-terrain.imag,terrain.real},configuration.l[m]);
        for (std::size_t j=0;j<configuration.Nj;++j) {
            const auto i=j+configuration.Nj*m; const auto& factor=descriptor.factors()[i];
            if (!factor.wave || std::hypot(configuration.k[m],configuration.l[m])>record.maximumForcedHorizontalWavenumber || configuration.j[j]>record.maximumForcedVerticalMode || (!damping.empty() && damping[i]!=0)) continue;
            const double response=configuration.g*bottom[i]*factor.NAp/factor.waveEnergy;
            operators.responsePlusX[i]=multiply(dx,response); operators.responsePlusY[i]=multiply(dy,response);
            operators.responseMinusX[i]=multiply(dx,-response); operators.responseMinusY[i]=multiply(dy,-response);
        }
    }
    forcing = std::make_unique<PseudoTopographicForcing>(entry,std::move(operators));
    return WVKernelStatus::ok();
}

WVKernelStatus prepareBoussinesqExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    const WVTransformBoussinesqKernel& kernel, WVForcingPreparation& preparation) {
    if (entry.name!="antialias filter") return WVKernelStatus::ok();
    const auto& g=kernel.geometry();
    const double retained=realValues(entry.configuration,"Nj")->front();
    const double cutoff=(2.0/3.0)*2*pi*static_cast<double>(g.Nx/2)/g.Lx;
    preparation.maximumHorizontalComponent=0; preparation.maximumVerticalMode=0;
    for(std::size_t m=0;m<g.Nkl;++m) for(double j:g.j) {
        if(std::hypot(g.k[m],g.l[m])>cutoff || j>retained-1) continue;
        preparation.maximumHorizontalComponent=std::max({preparation.maximumHorizontalComponent,std::abs(g.k[m]),std::abs(g.l[m])});
        preparation.maximumVerticalMode=std::max(preparation.maximumVerticalMode,static_cast<std::size_t>(j));
    }
    return WVKernelStatus::ok();
}

WVKernelStatus createBoussinesqExplicitAntialiasing(const WVFrozenForcingEntry& entry,
    WVTransformBoussinesqKernel& descriptor, const WVForcingPreparation &,
    std::unique_ptr<WVForcing>& forcing) {
    const auto& c=descriptor.geometry();
    const double retained=realValues(entry.configuration,"Nj")->front();
    const double cutoff=(2.0/3.0)*2*pi*static_cast<double>(c.Nx/2)/c.Lx;
    std::vector<std::size_t> indices;
    for (std::size_t i=0;i<c.Nj*c.Nkl;++i) if (std::hypot(c.k[i/c.Nj],c.l[i/c.Nj])>cutoff || c.j[i%c.Nj]>retained-1) indices.push_back(i);
    forcing = std::make_unique<ExplicitAntialiasingForcing>(entry,std::move(indices));
    return WVKernelStatus::ok();
}

} // namespace detail

WVConstantStratificationForcingEngine::~WVConstantStratificationForcingEngine() = default;

const WVForcingEvaluationDependencies*
WVConstantStratificationForcingEngine::forcingEvaluationDependencies(
    std::size_t index) const noexcept {
    const auto* forcing=forcingInstance(index);
    if(!forcing || !catalog_) return nullptr;
    const auto* registration=catalog_->forcings().registration(
        forcing->typeIdentifier(),forcing->contractVersion());
    return registration ? &registration->evaluationDependencies : nullptr;
}

WVKernelStatus WVConstantStratificationForcingEngine::validateSchedule(
    const WVTransformConstantStratificationConfiguration& configuration,
    const WVFrozenForcingSchedule& schedule,
    WVShape2D coefficientShape, const WVExtensionCatalog& catalog) {
    if (schedule.profileIdentifier != WVForcingScheduleProfileIdentifier || schedule.profileVersion != WVForcingScheduleProfileVersion) {
        return {WVKernelStatusCode::unsupportedOperation,"Unsupported frozen forcing schedule profile."};
    }
    if (coefficientShape.rows != configuration.Nj || coefficientShape.rows == 0 || coefficientShape.columns == 0) {
        return {WVKernelStatusCode::invalidShape,"Frozen forcing coefficient shape is incompatible with the model configuration."};
    }
    std::set<std::string> names;
    for (const auto& entry : schedule.entries) {
        const auto* registration = catalog.forcings().registration(entry.typeIdentifier, entry.contractVersion);
        if (registration == nullptr || !registration->isSupported || !registration->factory || registration->contractVersion != entry.contractVersion) return {WVKernelStatusCode::unsupportedOperation,"The frozen schedule has no matching paired C++ forcing implementation."};
        if (entry.stage != registration->stage) return {WVKernelStatusCode::invalidConfiguration,"A forcing record is assigned to the wrong execution stage."};
        if (entry.name.empty() || !names.insert(entry.name).second) return {WVKernelStatusCode::invalidConfiguration,"Forcing names must be nonempty and unique."};
        const auto recordStatus = catalog.forcings().validateConfiguration(entry);
        if (!recordStatus) return recordStatus;
        if (registration->modelPreflight) {
            const auto status = registration->modelPreflight(entry,configuration.shouldAntialias);
            if (!status) return status;
        }
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::create(
    const WVTransformConstantStratificationConfiguration& configuration,
    const WVFrozenForcingSchedule& schedule,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> fftEngine,
    std::unique_ptr<WVConstantStratificationForcingEngine>& forcingEngine) {
    if (!catalog) return {WVKernelStatusCode::invalidPointer,"Forcing-engine construction requires an extension catalog."};
    if (!fftEngine) return {WVKernelStatusCode::invalidPointer,"Forcing-engine construction requires an FFT engine."};
    try {
        auto candidate = std::unique_ptr<WVConstantStratificationForcingEngine>(new WVConstantStratificationForcingEngine());
        candidate->catalog_ = std::move(catalog);
        auto status = WVTransformConstantStratificationKernel::create(configuration,std::move(fftEngine),candidate->kernel_);
        if (!status) return status;
        status = validateSchedule(configuration,schedule,candidate->kernel_->descriptor().spectralShape(),*candidate->catalog_);
        if (!status) return status;
        status = candidate->initialize(schedule);
        if (!status) return status;
        forcingEngine = std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Forcing-engine allocation failed."};
    } catch (const std::exception& exception) {
        return {WVKernelStatusCode::invalidConfiguration,exception.what()};
    }
}

WVKernelStatus WVConstantStratificationForcingEngine::initialize(const WVFrozenForcingSchedule& schedule) {
    if (schedule.profileIdentifier != WVForcingScheduleProfileIdentifier || schedule.profileVersion != WVForcingScheduleProfileVersion) {
        return {WVKernelStatusCode::unsupportedOperation,"Unsupported frozen forcing schedule profile."};
    }
    std::vector<const WVFrozenForcingEntry*> entries;
    entries.reserve(schedule.entries.size());
    for (const auto& entry : schedule.entries) entries.push_back(&entry);
    for(const auto& entry:schedule.entries) {
        const auto* registration=catalog_->forcings().registration(
            entry.typeIdentifier,entry.contractVersion);
        if(!registration) return {WVKernelStatusCode::invalidConfiguration,
            "Forcing evaluation dependencies are unavailable."};
        for(std::size_t i=0;i<constantLaplacianUseCount_.size();++i)
            constantLaplacianUseCount_[i]+=
                registration->evaluationDependencies.constantLaplacianUseCount[i];
    }
    std::stable_sort(entries.begin(),entries.end(),[](const auto* left,const auto* right) {
        if (stageRank(left->stage) != stageRank(right->stage)) return stageRank(left->stage) < stageRank(right->stage);
        if (left->priority != right->priority) return left->priority < right->priority;
        return left->ordinal < right->ordinal;
    });

    std::set<std::string> names;
    const auto& descriptor = kernel_->descriptor();
    const auto count = descriptor.spectralShape().elementCount();
    const auto& configuration = descriptor.configuration();
    preparation_.maximumVerticalMode = configuration.Nj-1;
    for (const auto& mode : descriptor.fourierModes())
        preparation_.maximumHorizontalComponent = std::max({preparation_.maximumHorizontalComponent,std::abs(mode.k),std::abs(mode.l)});
    preparation_.hasAdaptiveDamping = std::any_of(
        entries.begin(), entries.end(), [this](const auto* entry) {
            const auto* registration = catalog_->forcings().registration(entry->typeIdentifier, entry->contractVersion);
            return registration != nullptr && registration->providesAdaptiveDamping;
        });

    for (const auto* entry : entries) {
        const auto* registration = catalog_->forcings().registration(entry->typeIdentifier,entry->contractVersion);
        if (registration->prepareConstantResolution) {
            const auto status = registration->prepareConstantResolution(*entry,descriptor,preparation_);
            if (!status) return status;
        }
    }

    for (const auto* entryPointer : entries) {
        const auto& entry = *entryPointer;
        if (entry.name.empty() || !names.insert(entry.name).second) return {WVKernelStatusCode::invalidConfiguration,"Forcing names must be nonempty and unique."};
        std::unique_ptr<WVForcing> resolved;
        auto status = catalog_->forcings().create(entry,descriptor,preparation_,resolved);
        if (!status) return status;
        if (!resolved || resolved->stage() != entry.stage) return {WVKernelStatusCode::invalidConfiguration,"A forcing factory returned an incompatible implementation."};
        if (entry.stage == WVForcingStage::spatial) ++metrics_.resolvedSpatialCount;
        else if (entry.stage == WVForcingStage::spectral) ++metrics_.resolvedSpectralCount;
        else ++metrics_.resolvedAmplitudeCount;
        metrics_.derivedOperatorBytes += resolved->persistentBytes();
        forcing_.push_back(std::move(resolved));
    }

    std::ostringstream identifier;
    identifier << WVForcingScheduleProfileIdentifier << ':';
    for (std::size_t index = 0; index < entries.size(); ++index) {
        if (index != 0) identifier << ',';
        identifier << entries[index]->typeIdentifier;
    }
    scheduleIdentifier_ = identifier.str();
    metrics_.scheduleBytes = scheduleIdentifier_.capacity()+
        forcing_.capacity()*sizeof(std::unique_ptr<WVForcing>);
    const auto R = descriptor.spatialShape().elementCount();
    const auto q = configuration.isHydrostatic ? 3U : 4U;
    const auto requiresPhysicalFields = std::any_of(forcing_.begin(),forcing_.end(),[](const auto& forcing) { return forcing->requiresPhysicalFields(); });
    const auto requiresForcingFields = std::any_of(forcing_.begin(),forcing_.end(),[](const auto& forcing) { return forcing->requiresForcingFields(); });
    const auto wholeFluxProducerCount = std::count_if(forcing_.begin(),forcing_.end(),[](const auto& forcing) { return forcing->producesCompleteFlux(); });
    if(requiresPhysicalFields && wholeFluxProducerCount==0) {
        physicalFieldCount_=2;
        auto status=kernel_->prepareHorizontalVelocityTransform(); if(!status) return status;
    }
    if (requiresPhysicalFields || wholeFluxProducerCount>1) physicalFields_.resize(physicalFieldCount_*R);
    if (requiresForcingFields) forcingFields_.resize(q*R);
    if (wholeFluxProducerCount > 1) {temporaryFlux_.resize(3*count); if(evaluationPolicy_==WVVariableEvaluationPolicy::reuse) nonlinearCache_.resize(3*count);}
    auto evaluationPrepared=evaluation_.prepare({{WVVariableEvaluationNode::physicalField,0},{WVVariableEvaluationNode::reduction,0},{WVVariableEvaluationNode::forcingTendency,0}});
    if(!evaluationPrepared) return evaluationPrepared;
    metrics_.workspaceCapacityBytes = vectorBytes(physicalFields_)+vectorBytes(forcingFields_)+vectorBytes(temporaryFlux_)+vectorBytes(nonlinearCache_);
    metrics_.workspaceHighWaterBytes = metrics_.workspaceCapacityBytes;
    metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
    metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::ensurePhysicalFields(
    const WVState& state, WVRealFieldBundleConstView& fields,
    WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared) {
    const auto& c = kernel_->descriptor().configuration();
    WVRealFieldBundleView local{physicalFields_.data(),{c.Nx,c.Ny,c.Nz,physicalFieldCount_}};
    auto& storage=externalFields ? *externalFields : local;
    if(!evaluation_.active()) { // Forcing diagnostics own their separate immutable event workspace.
        if(!externalFieldsPrepared) {
            auto status=storage.shape.fourth==2 ? kernel_->transformWaveVortexToUV(state,storage) : kernel_->transformWaveVortexToUVW(state,storage); if(!status) return status;
            externalFieldsPrepared=true; ++metrics_.physicalFieldReconstructionCount;
        } else ++metrics_.physicalFieldReuseCount;
    } else {
        const bool reused=evaluation_.ready({WVVariableEvaluationNode::physicalField,0});
        auto status=evaluation_.evaluate({WVVariableEvaluationNode::physicalField,0},storage.shape.elementCount()*sizeof(double),[&] {
            return storage.shape.fourth==2 ? kernel_->transformWaveVortexToUV(state,storage) : kernel_->transformWaveVortexToUVW(state,storage);
        });
        if(!status) return status;
        if(reused) ++metrics_.physicalFieldReuseCount; else ++metrics_.physicalFieldReconstructionCount;
        externalFieldsPrepared=true;
    }
    fields={storage.data,storage.shape}; return WVKernelStatus::ok();
}

WVRealFieldBundleView WVConstantStratificationForcingEngine::clearedSpatialTendency() {
    std::fill(forcingFields_.begin(),forcingFields_.end(),0.0);
    const auto& configuration = kernel_->descriptor().configuration();
    metrics_.spatialTendencyClearElementWrites += forcingFields_.size();
    const auto q = configuration.isHydrostatic ? 3U : 4U;
    return {forcingFields_.data(),{configuration.Nx,configuration.Ny,configuration.Nz,q}};
}

WVKernelStatus WVConstantStratificationForcingEngine::projectSpatialTendency(
    const WVState& state, const WVRealFieldBundleConstView& tendency,
    WVFlux& flux) {
    WVMutableCoefficients coefficients{{flux.Fp.data,flux.Fp.shape},{flux.Fm.data,flux.Fm.shape},{flux.F0.data,flux.F0.shape}};
    ++metrics_.spatialTendencyProjectionCount;
    return kernel_->descriptor().configuration().isHydrostatic ? kernel_->transformUVEtaToWaveVortex(tendency,state.t,state.t0,coefficients) : kernel_->transformUVWEtaToWaveVortex(tendency,state.t,state.t0,coefficients);
}

WVKernelStatus WVConstantStratificationForcingEngine::addProjectedSpatialTendency(
    const WVState& state, const WVRealFieldBundleConstView& tendency,
    WVFlux& flux, bool& outputInitialized) {
    if (diagnosticWorkspace_) {
        return diagnosticWorkspace_->addSpatial(tendency);
    }
    if (!outputInitialized) {
        const auto status = projectSpatialTendency(state,tendency,flux);
        if (status) outputInitialized = true;
        return status;
    }
    auto temporary = fluxViews(temporaryFlux_,kernel_->descriptor().spectralShape());
    const auto status = projectSpatialTendency(state,tendency,temporary);
    if (!status) return status;
    addFlux(temporary,flux);
    metrics_.temporaryAccumulationElementReads += 2*temporaryFlux_.size();
    metrics_.temporaryAccumulationElementWrites += temporaryFlux_.size();
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::addAdaptiveDamping(
    const WVState& state, const std::vector<double>& damping, WVFlux& flux,
    WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared) {
    WVRealFieldBundleConstView fields;
    const auto status = ensurePhysicalFields(state,fields,externalFields,externalFieldsPrepared);
    if (!status) return status;
    const auto R = kernel_->descriptor().spatialShape().elementCount();
    const auto calculateMaximum=[&](double& maximum) {
        maximum=0;
        for(std::size_t index=0;index<R;++index) maximum=std::max(maximum,std::sqrt(fields.data[index]*fields.data[index]+fields.data[R+index]*fields.data[R+index]));
        ++metrics_.horizontalSpeedReductionCount;
        return WVKernelStatus::ok();
    };
    const auto maximumStatus=diagnosticWorkspace_ ?
        diagnosticWorkspace_->evaluateHorizontalMaximum(horizontalMaximum_,calculateMaximum) :
        evaluation_.active() ? evaluation_.evaluate({WVVariableEvaluationNode::reduction,0},sizeof(double),[&] {return calculateMaximum(horizontalMaximum_);}) : calculateMaximum(horizontalMaximum_);
    if(!maximumStatus) return maximumStatus;
    const double maximumSpeed=horizontalMaximum_;
    const auto count = damping.size();
    for (std::size_t index = 0; index < count; ++index) {
        flux.Fp.data[index] = add(flux.Fp.data[index],multiply(state.coefficients.Ap.data[index],maximumSpeed*damping[index]));
        flux.Fm.data[index] = add(flux.Fm.data[index],multiply(state.coefficients.Am.data[index],maximumSpeed*damping[index]));
        flux.F0.data[index] = add(flux.F0.data[index],multiply(state.coefficients.A0.data[index],maximumSpeed*damping[index]));
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::addPseudoTopographicGeneration(const WVState& state, const WVPseudoTopographicOperators& operators, WVFlux& flux) {
    const auto& record = operators.configuration;
    const double elapsed = state.t-record.startTime;
    if (elapsed < 0.0) return WVKernelStatus::ok();
    const double ramp = record.rampDuration == 0.0 || elapsed >= record.rampDuration ? 1.0 : 0.5*(1.0-std::cos(pi*elapsed/record.rampDuration));
    const WVComplex64 oscillation{std::cos(record.frequency*elapsed),-std::sin(record.frequency*elapsed)};
    const double velocityX = ramp*multiply(record.barotropicVelocityAmplitude[0],oscillation).real;
    const double velocityY = ramp*multiply(record.barotropicVelocityAmplitude[1],oscillation).real;
    WVComplexConstView phases;
    const auto phaseStatus=kernel_->preparedPhase(state,phases); if(!phaseStatus) return phaseStatus;
    const auto count = operators.responsePlusX.size();
    for (std::size_t index = 0; index < count; ++index) {
        const auto phase=phases.data[index];
        const auto plus = add(multiply(operators.responsePlusX[index],velocityX),multiply(operators.responsePlusY[index],velocityY));
        const auto minus = add(multiply(operators.responseMinusX[index],velocityX),multiply(operators.responseMinusY[index],velocityY));
        flux.Fp.data[index] = add(flux.Fp.data[index],multiply(plus,conjugate(phase)));
        flux.Fm.data[index] = add(flux.Fm.data[index],multiply(minus,phase));
    }
    return WVKernelStatus::ok();
}

void WVConstantStratificationForcingEngine::addBetaPlaneAdvection(const WVState& state, const std::vector<WVComplex64>& betaA0, WVFlux& flux) const {
    for (std::size_t index = 0; index < betaA0.size(); ++index) flux.F0.data[index] = add(flux.F0.data[index],multiply(betaA0[index],state.coefficients.A0.data[index]));
}

WVKernelStatus WVConstantStratificationForcingEngine::addLinearCoefficientTendency(const WVState& state, double rate, WVFlux& flux) const {
    if (!std::isfinite(rate)) return {WVKernelStatusCode::invalidConfiguration,"Linear coefficient forcing rate must be finite."};
    const auto count = kernel_->descriptor().spectralShape().elementCount();
    for (std::size_t index = 0; index < count; ++index) {
        flux.Fp.data[index] = add(flux.Fp.data[index],multiply(state.coefficients.Ap.data[index],rate));
        flux.Fm.data[index] = add(flux.Fm.data[index],multiply(state.coefficients.Am.data[index],rate));
        flux.F0.data[index] = add(flux.F0.data[index],multiply(state.coefficients.A0.data[index],rate));
    }
    return WVKernelStatus::ok();
}

void WVConstantStratificationForcingEngine::initializeOutputWithZeros(WVFlux& flux, bool& outputInitialized) {
    const auto count = kernel_->descriptor().spectralShape().elementCount();
    std::fill_n(flux.Fp.data,count,WVComplex64{});
    std::fill_n(flux.Fm.data,count,WVComplex64{});
    std::fill_n(flux.F0.data,count,WVComplex64{});
    metrics_.accumulatorClearElementWrites += 3*count;
    outputInitialized = true;
}

WVKernelStatus WVConstantStratificationForcingEngine::addNonlinearFlux(
    const WVState& state, WVFlux& flux, bool& outputInitialized,
    WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared) {
    if (diagnosticWorkspace_) {
        WVRealFieldBundleConstView fields;
        auto status=ensurePhysicalFields(state,fields,externalFields,externalFieldsPrepared);
        if (!status) return status;
        diagnosticWorkspace_->spatialCaptured=true;
        auto raw=diagnosticWorkspace_->rawView();
        auto temporary=diagnosticWorkspace_->captureBuiltinProjection ? flux :
            diagnosticWorkspace_->temporaryView();
        return diagnosticWorkspace_->evaluateNonlinearRaw([&] {
            ++metrics_.nonlinearProducerCount;
            return kernel_->nonlinearFluxUsingAdvectionFields(state,temporary,fields,&raw,
                diagnosticWorkspace_->captureBuiltinProjection,
                diagnosticWorkspace_->stateDerivativeAccess());
        });
    }
    const auto evaluate = [&](WVFlux& destination) {
        ++metrics_.nonlinearProducerCount;
        // The resolved single nonlinear consumer can stream its private velocity.
        if(!externalFields && physicalFields_.empty()) return kernel_->nonlinearFlux(state,destination);
        const auto& c=kernel_->descriptor().configuration();
        WVRealFieldBundleView local{physicalFields_.data(),{c.Nx,c.Ny,c.Nz,physicalFieldCount_}};
        auto& storage=externalFields ? *externalFields : local;
        const bool reused=evaluation_.ready({WVVariableEvaluationNode::physicalField,0});
        auto status=reused ? kernel_->nonlinearFluxUsingAdvectionFields(state,destination,{storage.data,storage.shape}) :
            evaluation_.evaluate({WVVariableEvaluationNode::physicalField,0},storage.shape.elementCount()*sizeof(double),[&] {
                return kernel_->nonlinearFluxWithAdvectionFields(state,destination,storage);
            });
        if(status) {
            externalFieldsPrepared=true;
            if(reused) ++metrics_.physicalFieldReuseCount; else ++metrics_.physicalFieldReconstructionCount;
        }
        return status;
    };
    if(!temporaryFlux_.empty()) {
        auto& retained=nonlinearCache_.empty() ? temporaryFlux_ : nonlinearCache_;
        auto cached=fluxViews(retained,kernel_->descriptor().spectralShape());
        const auto status=evaluation_.evaluate({WVVariableEvaluationNode::forcingTendency,0},vectorBytes(retained),[&]{return evaluate(cached);});
        if(!status) return status;
        if(!outputInitialized) {
            const WVComplexConstView input[]={{cached.Fp.data,cached.Fp.shape},{cached.Fm.data,cached.Fm.shape},{cached.F0.data,cached.F0.shape}};
            const WVComplexView output[]={flux.Fp,flux.Fm,flux.F0};
            for(std::size_t family=0;family<3;++family) std::copy_n(input[family].data,input[family].shape.elementCount(),output[family].data);
            outputInitialized=true;
        } else addFlux(cached,flux);
        evaluation_.evict({WVVariableEvaluationNode::forcingTendency,0});
        return WVKernelStatus::ok();
    }
    if (!outputInitialized) {
        const auto status = evaluate(flux);
        if (status) outputInitialized = true;
        return status;
    }
    auto temporary = fluxViews(temporaryFlux_,kernel_->descriptor().spectralShape());
    const auto status = evaluate(temporary);
    if (!status) return status;
    addFlux(temporary,flux);
    metrics_.temporaryAccumulationElementReads += 2*temporaryFlux_.size();
    metrics_.temporaryAccumulationElementWrites += temporaryFlux_.size();
    return WVKernelStatus::ok();
}

void WVForcingExecutionContext::filterTendency(const std::vector<std::size_t>& indices) {
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    for (const auto index : indices) {
        flux_->Fp.data[index] = {};
        flux_->Fm.data[index] = {};
        flux_->F0.data[index] = {};
    }
}

WVKernelStatus WVForcingExecutionContext::laplacianDamping(double nu, double kappa, WVLaplacianDirection direction) {
    if (hydrostatic_) return hydrostatic_->addLaplacianDamping(*state_,nu,kappa,direction,*flux_);
    if (boussinesq_) return boussinesq_->addLaplacianDamping(*state_,nu,kappa,direction,*flux_);
    if (engine_->diagnosticWorkspace_)
        return engine_->diagnosticLaplacian(*state_,nu,kappa,direction,*flux_,*outputInitialized_);
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->kernel().addLaplacianDamping(*state_,nu,kappa,direction,*flux_);
}

WVKernelStatus WVForcingExecutionContext::verticalDiffusivity(double kappa,bool forceMean) {
    if (hydrostatic_) return hydrostatic_->addVerticalDiffusivity(*state_,kappa,forceMean,*flux_);
    if (boussinesq_) return boussinesq_->addVerticalDiffusivity(*state_,kappa,forceMean,*flux_);
    return laplacianDamping(0,kappa,WVLaplacianDirection::vertical);
}

WVKernelStatus WVForcingExecutionContext::nonlinearAdvection() {
    if (hydrostatic_) return hydrostatic_->addNonlinearFlux(*state_,*flux_);
    if (boussinesq_) return boussinesq_->addNonlinearFlux(*state_,*flux_);
    return engine_->addNonlinearFlux(*state_,*flux_,*outputInitialized_,externalFields_,*externalFieldsPrepared_);
}

WVKernelStatus WVForcingExecutionContext::physicalFields(
    WVRealFieldBundleConstView& fields) {
    if (hydrostatic_) return hydrostatic_->physicalFields(*state_,fields);
    if (boussinesq_) return boussinesq_->physicalFields(*state_,fields);
    return engine_->ensurePhysicalFields(*state_,fields,externalFields_,*externalFieldsPrepared_);
}

WVRealFieldBundleView WVForcingExecutionContext::clearedSpatialTendency() {
    if (hydrostatic_) return hydrostatic_->clearedSpatialTendency();
    if (boussinesq_) return boussinesq_->clearedSpatialTendency();
    return engine_->clearedSpatialTendency();
}

WVKernelStatus WVForcingExecutionContext::projectSpatialTendency(
    const WVRealFieldBundleConstView& tendency) {
    if (hydrostatic_) return hydrostatic_->addProjectedSpatialTendency(*state_,tendency,*flux_);
    if (boussinesq_) return boussinesq_->addProjectedSpatialTendency(*state_,tendency,*flux_);
    return engine_->addProjectedSpatialTendency(*state_,tendency,*flux_,*outputInitialized_);
}

WVKernelStatus WVForcingExecutionContext::adaptiveDamping(const std::vector<double>& values) {
    if (hydrostatic_) return hydrostatic_->addAdaptiveDamping(*state_,values,*flux_);
    if (boussinesq_) return boussinesq_->addAdaptiveDamping(*state_,values,*flux_);
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addAdaptiveDamping(*state_,values,*flux_,externalFields_,*externalFieldsPrepared_);
}

WVKernelStatus WVForcingExecutionContext::pseudoTopographicGeneration(const WVPseudoTopographicOperators& operators) {
    if (hydrostatic_) return hydrostatic_->addPseudoTopographicGeneration(*state_,operators,*flux_);
    if (boussinesq_) return boussinesq_->addPseudoTopographicGeneration(*state_,operators,*flux_);
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addPseudoTopographicGeneration(*state_,operators,*flux_);
}

WVKernelStatus WVForcingExecutionContext::betaPlaneAdvection(const std::vector<WVComplex64>& values) {
    if (hydrostatic_ || boussinesq_) {
        for (std::size_t i=0;i<values.size();++i) flux_->F0.data[i]=add(flux_->F0.data[i],multiply(values[i],state_->coefficients.A0.data[i]));
        return WVKernelStatus::ok();
    }
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    engine_->addBetaPlaneAdvection(*state_,values,*flux_);
    return WVKernelStatus::ok();
}

void WVForcingExecutionContext::zeroSelectedTendencies(const WVFixedAmplitudeConfiguration& values) {
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    for (const auto index : values.ApIndices) flux_->Fp.data[index] = {};
    for (const auto index : values.AmIndices) flux_->Fm.data[index] = {};
    for (const auto index : values.A0Indices) flux_->F0.data[index] = {};
}

WVKernelStatus WVForcingExecutionContext::linearCoefficientTendency(double rate) {
    if (hydrostatic_ || boussinesq_) {
        if (!std::isfinite(rate)) return {WVKernelStatusCode::invalidConfiguration,"Nonfinite coefficient rate."};
        const WVComplexConstView a[]={state_->coefficients.Ap,state_->coefficients.Am,state_->coefficients.A0};
        const WVComplexView b[]={flux_->Fp,flux_->Fm,flux_->F0};
        for (int j=0;j<3;++j) for (std::size_t i=0;i<a[j].shape.elementCount();++i) b[j].data[i]=add(b[j].data[i],multiply(a[j].data[i],rate));
        return WVKernelStatus::ok();
    }
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addLinearCoefficientTendency(*state_,rate,*flux_);
}

WVKernelStatus WVConstantStratificationForcingEngine::validateVariableEvaluationPolicyChange(
    WVVariableEvaluationPolicy policy) const noexcept {
    if(evaluation_.active() || executing_) return {WVKernelStatusCode::reentrantExecution,"Cannot change an active evaluation policy."};
    if(policy!=WVVariableEvaluationPolicy::reuse && policy!=WVVariableEvaluationPolicy::lowMemory)
        return {WVKernelStatusCode::invalidConfiguration,"Unknown variable evaluation policy."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVConstantStratificationForcingEngine::setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy) {
    auto status=validateVariableEvaluationPolicyChange(policy); if(!status) return status;
    const auto wholeFluxCount=std::count_if(forcing_.begin(),forcing_.end(),[](const auto& forcing){return forcing->producesCompleteFlux();});
    std::vector<WVComplex64> preparedNonlinear;
    try {
        if(policy==WVVariableEvaluationPolicy::reuse && wholeFluxCount>1 && nonlinearCache_.empty())
            preparedNonlinear.resize(3*kernel_->descriptor().spectralShape().elementCount());
    } catch(const std::bad_alloc&) {return {WVKernelStatusCode::allocationFailure,"Unable to allocate shared nonlinear result."};}
    if(policy==WVVariableEvaluationPolicy::reuse) {
        if(!preparedNonlinear.empty()) nonlinearCache_.swap(preparedNonlinear);
    } else std::vector<WVComplex64>{}.swap(nonlinearCache_);
    evaluationPolicy_=policy;
    metrics_.workspaceCapacityBytes=vectorBytes(physicalFields_)+vectorBytes(forcingFields_)+vectorBytes(temporaryFlux_)+vectorBytes(nonlinearCache_);
    metrics_.workspaceHighWaterBytes=std::max(metrics_.workspaceHighWaterBytes,metrics_.workspaceCapacityBytes);
    metrics_.workspaceLiveBytes=metrics_.workspaceCapacityBytes;
    metrics_.workspaceMaximumLiveBytes=std::max(metrics_.workspaceMaximumLiveBytes,metrics_.workspaceLiveBytes);
    return WVKernelStatus::ok();
}
WVKernelStatus WVConstantStratificationForcingEngine::beginStateEvaluation(const WVState& state) {
    if(evaluation_.active() || executing_) return {WVKernelStatusCode::reentrantExecution,"Constant state evaluation is already active."};
    auto status=kernel_->beginStateEvaluation(state); if(!status) return status;
    status=evaluation_.begin(this,evaluationPolicy_);
    if(!status) {kernel_->endStateEvaluation(); return status;}
    evaluationState_=state; evaluationFields_={}; return WVKernelStatus::ok();
}
void WVConstantStratificationForcingEngine::endStateEvaluation() noexcept {
    evaluation_.end(); kernel_->endStateEvaluation(); evaluationState_={}; evaluationFields_={};
}
WVKernelStatus WVConstantStratificationForcingEngine::validateStateEvaluation(const WVState& state) const {
    if(!evaluation_.active() || state.t!=evaluationState_.t || state.t0!=evaluationState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"State does not belong to the active constant evaluation."};
    const WVComplexConstView a[]={state.coefficients.Ap,state.coefficients.Am,state.coefficients.A0};
    const WVComplexConstView b[]={evaluationState_.coefficients.Ap,evaluationState_.coefficients.Am,evaluationState_.coefficients.A0};
    for(std::size_t i=0;i<3;++i) if(a[i].data!=b[i].data || a[i].shape.rows!=b[i].shape.rows || a[i].shape.columns!=b[i].shape.columns)
        return {WVKernelStatusCode::invalidConfiguration,"Coefficients do not belong to the active constant evaluation."};
    return WVKernelStatus::ok();
}
WVKernelStatus WVConstantStratificationForcingEngine::nonlinearFlux(const WVState& state, WVFlux& flux) {
    return nonlinearFluxImpl(state,flux,nullptr,nullptr);
}

WVKernelStatus WVConstantStratificationForcingEngine::diagnosticLaplacian(
    const WVState& state,double nu,double kappa,WVLaplacianDirection direction,
    WVFlux& flux,bool& initialized) {
    auto& work=*diagnosticWorkspace_;
    const auto S=work.spectral.elementCount();
    const auto R=work.spatial.first*work.spatial.second*work.spatial.third;
    if(work.laplacianCoefficients.size()!=3*S || work.laplacianFields.size()!=4*R)
        return {WVKernelStatusCode::invalidConfiguration,
            "Constant Laplacian diagnostic storage was not prepared."};
    const auto& descriptor=kernel_->descriptor();
    const auto Nj=descriptor.configuration().Nj;
    const auto directionIndex=static_cast<std::size_t>(direction);
    WVRealFieldBundleConstView unscaled;
    auto status=work.evaluateConstantLaplacian(directionIndex,unscaled,
        [&](WVRealFieldBundleView destination) {
            const WVComplexConstView inputs[]={state.coefficients.Ap,
                state.coefficients.Am,state.coefficients.A0};
            for(std::size_t i=0;i<S;++i) {
                const double waveNumber=direction==WVLaplacianDirection::horizontal ?
                    descriptor.fourierModes()[i/Nj].Kh :
                    descriptor.verticalModes().verticalWavenumber[i%Nj];
                for(std::size_t family=0;family<3;++family)
                    work.laplacianCoefficients[family*S+i]=multiply(
                        inputs[family].data[i],-waveNumber*waveNumber);
            }
            const WVState derivative{state.t,state.t0,
                {{work.laplacianCoefficients.data(),work.spectral},
                 {work.laplacianCoefficients.data()+S,work.spectral},
                 {work.laplacianCoefficients.data()+2*S,work.spectral}}};
            const auto result=kernel_->transformCoefficientTendencyToUVWEta(
                derivative,destination);
            if(result) ++metrics_.constantLaplacianProducerCount[directionIndex];
            return result;
        });
    if (!status) return status;
    WVRealFieldBundleView fields{work.laplacianFields.data(),
        {work.spatial.first,work.spatial.second,work.spatial.third,4}};
    for (std::size_t channel=0;channel<work.spatial.fourth;++channel) {
        const bool density=channel==work.spatial.fourth-1;
        const auto source=density ? 3 : channel;
        for (std::size_t i=0;i<R;++i)
            fields.data[channel*R+i]=(density ? kappa : nu)*unscaled.data[source*R+i];
    }
    fields.shape.fourth=work.spatial.fourth;
    return addProjectedSpatialTendency(state,{fields.data,fields.shape},flux,initialized);
}

WVKernelStatus WVConstantStratificationForcingEngine::evaluateForcingTendencies(
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,
    const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session) {
    return evaluateForcingTendenciesImpl(state,outputs,count,preparedPhysical,
        session,nullptr);
}
WVKernelStatus WVConstantStratificationForcingEngine::evaluateForcingTendenciesImpl(
    const WVState& state,const WVForcingTendencyOutput* outputs,std::size_t count,
    const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session,
    WVFlux* builtinNonlinearProjection) {
    if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing diagnostics require an idle engine."};
    tendencyMetrics_.workspaceLastPeakBytes=0;
    const auto& c=kernel_->descriptor().configuration();
    const WVShape4D spatial{c.Nx,c.Ny,c.Nz,c.isHydrostatic ? 3U : 4U};
    auto status=detail::validateForcingTendencyOutputs(forcing_,stateShape(),spatial,state,outputs,count);
    if (!status || (!count && !builtinNonlinearProjection)) return status;
    status=detail::validatePreparedDiagnosticFields(preparedPhysical,spatial,3,state,outputs,count);
    if (!status) return status;
    detail::WVScopedStateEvaluation<WVTransformConstantStratificationKernel> kernelScope(kernel(),state);
    if(!kernelScope.status()) return kernelScope.status();
    try {
        detail::WVForcingDiagnosticLedger localLedger(tendencyMetrics_);
        std::unique_ptr<detail::WVForcingDiagnosticWorkspace> local;
        if (!session) {
            local=std::make_unique<detail::WVForcingDiagnosticWorkspace>(stateShape(),spatial);
            local->constantLaplacianUseCount=constantLaplacianUseCount_;
            local->requiresFourChannelTendencySelection=c.isHydrostatic;
            std::vector<WVForcingStage> stages;
            local->nonlinearUseCount=0;
            for(std::size_t index=0;index<forcing_.size();++index) {
                const auto& forcing=forcing_[index];
                stages.push_back(forcing->stage());
                const auto* dependencies=forcingEvaluationDependencies(index);
                if(!dependencies) return {WVKernelStatusCode::invalidConfiguration,
                    "Forcing evaluation dependencies are unavailable."};
                local->nonlinearUseCount+=dependencies->nonlinearUseCount;
            }
            status=localLedger.context.prepare(detail::WVForcingDiagnosticWorkspace::dependencyKeys(
                forcing_.size(),nullptr,&constantLaplacianUseCount_));
            if(!status) return status;
            status=localLedger.context.begin(this,evaluationPolicy_); if(!status) return status;
            status=local->beginScopedEvaluation(localLedger.context,stages); if(!status) return status;
            session=local.get();
        }
        auto& work=*session;
        status=work.bind(this,state); if (!status) return status;
        const auto spectral=stateShape(); const auto S=spectral.elementCount();
        const auto R=static_cast<std::size_t>(c.Nx)*c.Ny*c.Nz;
        if (work.spectral.rows!=spectral.rows || work.spectral.columns!=spectral.columns ||
            work.spatial.first!=spatial.first || work.spatial.second!=spatial.second ||
            work.spatial.third!=spatial.third || work.spatial.fourth!=spatial.fourth ||
            work.flux.size()!=3*S || work.previous.size()!=3*S || work.temporary.size()!=3*S ||
            work.cumulative.size()!=spatial.elementCount() || work.raw.size()!=spatial.elementCount() ||
            work.physical.size()!=4*R)
            return {WVKernelStatusCode::invalidShape,"Forcing diagnostic session has incompatible constant-stratification storage."};
        auto flux=work.fluxView();
        if (!work.initialized()) {
            if (preparedPhysical) {
                std::copy_n(preparedPhysical->data,preparedPhysical->shape.elementCount(),work.physical.data());
                work.physicalPrepared=true;
            }
            status=validateStateAndFlux(kernel_->descriptor(),state,flux); if (!status) return status;
            work.markInitialized();
        }
        executing_=true; diagnosticWorkspace_=&work;
        struct Guard {
            WVConstantStratificationForcingEngine& engine;
            ~Guard() {
                engine.tendencyMetrics_.workspaceLastPeakBytes=engine.diagnosticWorkspace_->bytes();
                engine.tendencyMetrics_.workspaceHighWaterBytes=std::max(engine.tendencyMetrics_.workspaceHighWaterBytes,engine.diagnosticWorkspace_->bytes());
                engine.tendencyMetrics_.workspaceLiveBytes=0;
                engine.diagnosticWorkspace_=nullptr; engine.executing_=false;
            }
        } guard{*this};
        tendencyMetrics_.workspaceLiveBytes=work.bytes();
        WVRealFieldBundleView physical{work.physical.data(),{c.Nx,c.Ny,c.Nz,3}};
        bool initialized=true;
        WVForcingExecutionContext context;
        context.engine_=this; context.state_=&state; context.outputInitialized_=&initialized;
        context.externalFields_=&physical; context.externalFieldsPrepared_=&work.physicalPrepared;
        return detail::evaluateForcingTendencySequence(forcing_,work,outputs,count,tendencyMetrics_,
            [&](const WVForcing& forcing,WVFlux& destination) {
                context.flux_=&destination;
                return forcing.addRightHandSide(context);
            },
            [&](WVRealFieldBundleConstView fields,WVFlux& destination) {
                return projectSpatialTendency(state,fields,destination);
            },
            [&](const std::vector<WVComplex64>& difference,WVRealFieldBundleView destination) {
                const auto S=stateShape().elementCount();
                const WVState delta{state.t,state.t0,{{difference.data(),stateShape()},
                    {difference.data()+S,stateShape()},{difference.data()+2*S,stateShape()}}};
                if (!c.isHydrostatic) return kernel_->transformCoefficientTendencyToUVWEta(delta,destination);
                const auto R=kernel_->descriptor().spatialShape().elementCount();
                if(work.laplacianFields.size()!=4*R)
                    return WVKernelStatus{WVKernelStatusCode::invalidShape,
                        "Hydrostatic coefficient-tendency selection storage was not prepared."};
                WVRealFieldBundleView fields{work.laplacianFields.data(),{c.Nx,c.Ny,c.Nz,4}};
                auto result=kernel_->transformCoefficientTendencyToUVWEta(delta,fields);
                if (result) {
                    std::copy_n(fields.data,2*R,destination.data);
                    std::copy_n(fields.data+3*R,R,destination.data+2*R);
                }
                return result;
            },builtinNonlinearProjection);
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Unable to allocate event forcing diagnostic workspace."};
    }
}

WVKernelStatus WVConstantStratificationForcingEngine::evaluateRightHandSideWithContext(
    const WVState& state, WVFlux& flux,
    WVRealFieldBundleView& advectionFieldStorage,
    WVConstantStratificationRightHandSideContext& context) {
    if(!evaluation_.active()) return {WVKernelStatusCode::invalidConfiguration,
        "Borrowed RHS fields require an explicit immutable evaluation scope."};
    return nonlinearFluxImpl(state,flux,&advectionFieldStorage,&context);
}

WVKernelStatus WVConstantStratificationForcingEngine::nonlinearFluxImpl(
    const WVState& state, WVFlux& flux, WVRealFieldBundleView* externalFields,
    WVConstantStratificationRightHandSideContext* context) {
    if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing-engine execution is not reentrant."};
    detail::WVScopedStateEvaluation<WVConstantStratificationForcingEngine> evaluation(*this,state);
    if(!evaluation.status()) return evaluation.status();
    const auto validation = validateStateAndFlux(kernel_->descriptor(),state,flux);
    if (!validation) return validation;
    const bool requestedAdvection=externalFields!=nullptr;
    const auto& c=kernel_->descriptor().configuration();
    if(externalFields || !physicalFields_.empty()) {
        if(evaluationFields_.data && externalFields && externalFields->data!=evaluationFields_.data)
            return {WVKernelStatusCode::invalidConfiguration,"Advection storage cannot change inside an immutable evaluation."};
        if(!evaluationFields_.data) evaluationFields_=externalFields ? *externalFields : WVRealFieldBundleView{physicalFields_.data(),{c.Nx,c.Ny,c.Nz,physicalFieldCount_}};
        externalFields=&evaluationFields_;
    }
    if (context != nullptr) *context = {};
    executing_ = true;
    struct Guard { bool& value; ~Guard() { value = false; } } guard{executing_};
    bool outputInitialized = false;
    bool externalFieldsPrepared = evaluation_.ready({WVVariableEvaluationNode::physicalField,0});
    WVForcingExecutionContext forcingContext;
    forcingContext.engine_ = this;
    forcingContext.state_ = &state;
    forcingContext.flux_ = &flux;
    forcingContext.outputInitialized_ = &outputInitialized;
    forcingContext.externalFields_ = externalFields;
    forcingContext.externalFieldsPrepared_ = &externalFieldsPrepared;
    if (!linearDynamics_) for (const auto& forcing : forcing_) {
        const auto status = forcing->addRightHandSide(forcingContext);
        if (!status) return status;
    }
    if (!outputInitialized) initializeOutputWithZeros(flux,outputInitialized);
    if (requestedAdvection && !externalFieldsPrepared) {
        WVRealFieldBundleConstView ignored;
        auto status = ensurePhysicalFields(state,ignored,externalFields,externalFieldsPrepared);
        if (!status) return status;
    }
    if (context != nullptr) {
        context->owner_ = this;
        context->advectionFields_ = {externalFields->data,externalFields->shape};
        context->evaluation_ = &evaluation_;
        context->generation_ = evaluation_.generation();
    }
    ++metrics_.evaluationCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::advectFGridScalar(
    const WVConstantStratificationRightHandSideContext& context,
    const WVRealVolumeConstView& scalar, bool shouldAntialias,
    WVRealVolumeView& rightHandSide) {
    if (context.owner_ != this || context.evaluation_ != &evaluation_ || !context.hasAdvectionFields())
        return {WVKernelStatusCode::invalidConfiguration,"The RHS evaluation context is stale or belongs to another forcing engine."};
    return kernel_->advectFGridScalar(scalar,context.advectionFields_,shouldAntialias,rightHandSide);
}

WVStateConstraintResult WVConstantStratificationForcingEngine::restoreForcingAmplitudes(WVMutableCoefficients& coefficients) {
    if(evaluation_.active()) return {{WVKernelStatusCode::invalidConfiguration,"Cannot constrain coefficients during an immutable evaluation."},0,false};
    const auto status = validateMutableCoefficients(kernel_->descriptor().spectralShape(),coefficients);
    if (!status) return {status,0,false};
    std::size_t modifiedCoefficientCount = 0;
    bool fsalCompatible = true;
    for (const auto& forcing : forcing_) {
        const auto result = forcing->applyConstraint(coefficients);
        if (!result.status) return result;
        modifiedCoefficientCount += result.modifiedCoefficientCount;
        fsalCompatible = fsalCompatible && result.fsalCompatible;
        metrics_.restoredCoefficientCount += forcing->constraintWriteCount();
        metrics_.stateConstraintElementWrites += forcing->constraintWriteCount();
    }
    return {WVKernelStatus::ok(),modifiedCoefficientCount,fsalCompatible};
}

WVKernelStatus WVConstantStratificationForcingEngine::createErrorPolicy(
    double absoluteToleranceScale,
    std::unique_ptr<WVIntegrationErrorPolicy>& policy) const {
    return WVWaveVortexCoefficientErrorPolicy::create(kernel_->descriptor(),absoluteToleranceScale,policy);
}

std::size_t WVConstantStratificationForcingEngine::persistentBytes() const noexcept {
    return sizeof(*this)+evaluation_.persistentBytes()+kernel_->persistentBytes()+metrics_.scheduleBytes+metrics_.derivedOperatorBytes+metrics_.workspaceCapacityBytes;
}

} // namespace wavevortex::runtime
