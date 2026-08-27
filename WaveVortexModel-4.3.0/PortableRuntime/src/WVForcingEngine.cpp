#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"

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

class NonlinearAdvectionForcing final : public ResolvedForcing {
public:
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

std::vector<double> adaptiveDampingOperator(const WVTransformConstantStratificationDescriptor& descriptor, WVKernelStatus& status) {
    const auto& configuration = descriptor.configuration();
    const auto& modes = descriptor.verticalModes();
    const double f2 = modes.coriolisFrequency*modes.coriolisFrequency;
    if (!(f2 > 0.0)) { status = {WVKernelStatusCode::unsupportedOperation,"Adaptive damping requires nonzero Coriolis frequency."}; return {}; }
    double maximumComponent = 0.0;
    for (const auto& horizontal : descriptor.fourierModes()) maximumComponent = std::max(maximumComponent,std::max(std::abs(horizontal.k),std::abs(horizontal.l)));
    if (!(maximumComponent > 0.0)) { status = {WVKernelStatusCode::invalidConfiguration,"Adaptive damping requires resolved horizontal wavenumbers."}; return {}; }
    const double effectiveResolution = pi/maximumComponent;
    const double dklMinimum = std::min(2.0*pi/configuration.Lx,2.0*pi/configuration.Ly);
    const double klCutoff = dklMinimum*std::pow(maximumComponent/dklMinimum,0.75);
    const double jMaximum = static_cast<double>(configuration.Nj-1);
    const double jCutoff = configuration.Nj > 2 ? std::pow(jMaximum,0.75) : 0.0;
    const double prefactorXY = effectiveResolution/(pi*pi);
    const double referenceLr2 = configuration.g*modes.h0.back()/f2;
    const double prefactorZ = (pi*pi*referenceLr2/(effectiveResolution*effectiveResolution))*prefactorXY;
    std::vector<double> damping(descriptor.spectralShape().elementCount());
    for (std::size_t iMode = 0; iMode < descriptor.Nkl(); ++iMode) {
        const auto& horizontal = descriptor.fourierModes()[iMode];
        const double Qkl = vanishingFilter(horizontal.Kh,klCutoff,maximumComponent);
        for (std::size_t iJ = 0; iJ < configuration.Nj; ++iJ) {
            const auto index = iJ+configuration.Nj*iMode;
            const double Qj = configuration.Nj > 2 ? vanishingFilter(modes.j[iJ],jCutoff,jMaximum) : 1.0;
            const double Lr2Inverse = iJ == 0 ? 0.0 : f2/(configuration.g*modes.h0[iJ]);
            damping[index] = -prefactorXY*Qkl*horizontal.Kh*horizontal.Kh-prefactorZ*Qj*Lr2Inverse;
        }
    }
    status = WVKernelStatus::ok();
    return damping;
}

bool emptyConfiguration(const WVFrozenForcingEntry& entry) { return entry.configuration.values.empty(); }

} // namespace

namespace detail {

WVKernelStatus createNonlinearAdvectionForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor&, bool,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Nonlinear advection does not accept configuration values."};
    forcing = std::make_unique<NonlinearAdvectionForcing>(entry);
    return WVKernelStatus::ok();
}

WVKernelStatus createAdaptiveDampingForcing(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, bool,
    std::unique_ptr<WVForcing>& forcing) {
    if (!emptyConfiguration(entry)) return {WVKernelStatusCode::invalidConfiguration,"Adaptive damping does not accept configuration values."};
    WVKernelStatus status;
    auto damping = adaptiveDampingOperator(descriptor,status);
    if (!status) return status;
    forcing = std::make_unique<AdaptiveDampingForcing>(entry,std::move(damping));
    return WVKernelStatus::ok();
}

WVKernelStatus createQuadraticBottomFriction(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, bool,
    std::unique_ptr<WVForcing>& forcing) {
    const auto* values = realValues(entry.configuration,"Cd");
    if (values == nullptr || values->size() != 1 || !std::isfinite(values->front()) || values->front() < 0.0) return {WVKernelStatusCode::invalidConfiguration,"Quadratic drag coefficient must be one finite nonnegative scalar."};
    forcing = std::make_unique<QuadraticBottomFrictionForcing>(entry,values->front()/descriptor.bottomQuadratureWeight());
    return WVKernelStatus::ok();
}

WVKernelStatus createLinearBottomFriction(
    const WVFrozenForcingEntry& entry,
    const WVTransformConstantStratificationDescriptor& descriptor, bool,
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
    const WVTransformConstantStratificationDescriptor& descriptor, bool,
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
    const WVTransformConstantStratificationDescriptor& descriptor, bool,
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
    bool hasAdaptiveDamping, std::unique_ptr<WVForcing>& forcing) {
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
    const auto damping = hasAdaptiveDamping && record.shouldAvoidAdaptiveDamping ? adaptiveDampingOperator(descriptor,dampingStatus) : std::vector<double>{};
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

} // namespace detail

WVConstantStratificationForcingEngine::~WVConstantStratificationForcingEngine() = default;

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
    std::stable_sort(entries.begin(),entries.end(),[](const auto* left,const auto* right) {
        if (stageRank(left->stage) != stageRank(right->stage)) return stageRank(left->stage) < stageRank(right->stage);
        if (left->priority != right->priority) return left->priority < right->priority;
        return left->ordinal < right->ordinal;
    });

    std::set<std::string> names;
    const auto& descriptor = kernel_->descriptor();
    const auto count = descriptor.spectralShape().elementCount();
    const auto& configuration = descriptor.configuration();
    const bool hasAdaptiveDamping = std::any_of(
        entries.begin(), entries.end(), [this](const auto* entry) {
            const auto* registration = catalog_->forcings().registration(entry->typeIdentifier, entry->contractVersion);
            return registration != nullptr && registration->providesAdaptiveDamping;
        });

    for (const auto* entryPointer : entries) {
        const auto& entry = *entryPointer;
        if (entry.name.empty() || !names.insert(entry.name).second) return {WVKernelStatusCode::invalidConfiguration,"Forcing names must be nonempty and unique."};
        std::unique_ptr<WVForcing> resolved;
        auto status = catalog_->forcings().create(entry,descriptor,hasAdaptiveDamping,resolved);
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
    if (requiresPhysicalFields) physicalFields_.resize(4*R);
    if (requiresForcingFields) forcingFields_.resize(q*R);
    if (wholeFluxProducerCount > 1) temporaryFlux_.resize(3*count);
    metrics_.workspaceCapacityBytes = vectorBytes(physicalFields_)+vectorBytes(forcingFields_)+vectorBytes(temporaryFlux_);
    metrics_.workspaceHighWaterBytes = metrics_.workspaceCapacityBytes;
    metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
    metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::ensurePhysicalFields(
    const WVState& state, WVRealFieldBundleConstView& fields,
    WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared) {
    const auto& configuration = kernel_->descriptor().configuration();
    if (externalFields != nullptr) {
        if (!externalFieldsPrepared) {
            auto status = kernel_->transformWaveVortexToUVW(state,*externalFields);
            if (!status) return status;
            externalFieldsPrepared = true;
            ++metrics_.physicalFieldReconstructionCount;
        } else {
            ++metrics_.physicalFieldReuseCount;
        }
        fields = {externalFields->data,externalFields->shape};
        return WVKernelStatus::ok();
    }
    const WVShape4D shape{configuration.Nx,configuration.Ny,configuration.Nz,4};
    if (!physicalFieldsValid_) {
        WVRealFieldBundleView mutableFields{physicalFields_.data(),shape};
        const auto status = kernel_->transformWaveVortexToUVWEta(state,mutableFields);
        if (!status) return status;
        physicalFieldsValid_ = true;
        ++metrics_.physicalFieldReconstructionCount;
    } else {
        ++metrics_.physicalFieldReuseCount;
    }
    fields = {physicalFields_.data(),shape};
    return WVKernelStatus::ok();
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
    double maximumSpeed = 0.0;
    for (std::size_t index = 0; index < R; ++index) maximumSpeed = std::max(maximumSpeed,std::sqrt(fields.data[index]*fields.data[index]+fields.data[R+index]*fields.data[R+index]));
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
    const auto& modes = kernel_->descriptor().verticalModes();
    const auto count = operators.responsePlusX.size();
    for (std::size_t index = 0; index < count; ++index) {
        const double angle = modes.omega[index]*(state.t-state.t0);
        const WVComplex64 phase{std::cos(angle),std::sin(angle)};
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
    const auto evaluate = [&](WVFlux& destination) {
        if (externalFields == nullptr) return kernel_->nonlinearFlux(state,destination);
        const bool reconstructsFields = !externalFieldsPrepared;
        auto status = externalFieldsPrepared ? kernel_->nonlinearFluxUsingAdvectionFields(state,destination,{externalFields->data,externalFields->shape}) : kernel_->nonlinearFluxWithAdvectionFields(state,destination,*externalFields);
        if (status) {
            externalFieldsPrepared = true;
            if (reconstructsFields) ++metrics_.physicalFieldReconstructionCount;
            else ++metrics_.physicalFieldReuseCount;
        }
        return status;
    };
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

WVKernelStatus WVForcingExecutionContext::nonlinearAdvection() {
    return engine_->addNonlinearFlux(*state_,*flux_,*outputInitialized_,externalFields_,*externalFieldsPrepared_);
}

WVKernelStatus WVForcingExecutionContext::physicalFields(
    WVRealFieldBundleConstView& fields) {
    return engine_->ensurePhysicalFields(*state_,fields,externalFields_,*externalFieldsPrepared_);
}

WVRealFieldBundleView WVForcingExecutionContext::clearedSpatialTendency() {
    return engine_->clearedSpatialTendency();
}

WVKernelStatus WVForcingExecutionContext::projectSpatialTendency(
    const WVRealFieldBundleConstView& tendency) {
    return engine_->addProjectedSpatialTendency(*state_,tendency,*flux_,*outputInitialized_);
}

WVKernelStatus WVForcingExecutionContext::adaptiveDamping(const std::vector<double>& values) {
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addAdaptiveDamping(*state_,values,*flux_,externalFields_,*externalFieldsPrepared_);
}

WVKernelStatus WVForcingExecutionContext::pseudoTopographicGeneration(const WVPseudoTopographicOperators& operators) {
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addPseudoTopographicGeneration(*state_,operators,*flux_);
}

WVKernelStatus WVForcingExecutionContext::betaPlaneAdvection(const std::vector<WVComplex64>& values) {
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
    if (!*outputInitialized_) engine_->initializeOutputWithZeros(*flux_,*outputInitialized_);
    return engine_->addLinearCoefficientTendency(*state_,rate,*flux_);
}

WVKernelStatus WVConstantStratificationForcingEngine::nonlinearFlux(const WVState& state, WVFlux& flux) {
    return nonlinearFluxImpl(state,flux,nullptr,nullptr);
}

WVKernelStatus WVConstantStratificationForcingEngine::evaluateRightHandSideWithContext(
    const WVState& state, WVFlux& flux,
    WVRealFieldBundleView& advectionFieldStorage,
    WVConstantStratificationRightHandSideContext& context) {
    return nonlinearFluxImpl(state,flux,&advectionFieldStorage,&context);
}

WVKernelStatus WVConstantStratificationForcingEngine::nonlinearFluxImpl(
    const WVState& state, WVFlux& flux, WVRealFieldBundleView* externalFields,
    WVConstantStratificationRightHandSideContext* context) {
    if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing-engine execution is not reentrant."};
    const auto validation = validateStateAndFlux(kernel_->descriptor(),state,flux);
    if (!validation) return validation;
    ++evaluationGeneration_;
    if (context != nullptr) *context = {};
    executing_ = true;
    struct Guard { bool& value; ~Guard() { value = false; } } guard{executing_};
    clearEvaluationWorkspace();
    bool outputInitialized = false;
    bool externalFieldsPrepared = false;
    WVForcingExecutionContext forcingContext;
    forcingContext.engine_ = this;
    forcingContext.state_ = &state;
    forcingContext.flux_ = &flux;
    forcingContext.outputInitialized_ = &outputInitialized;
    forcingContext.externalFields_ = externalFields;
    forcingContext.externalFieldsPrepared_ = &externalFieldsPrepared;
    for (const auto& forcing : forcing_) {
        const auto status = forcing->addRightHandSide(forcingContext);
        if (!status) return status;
    }
    if (!outputInitialized) initializeOutputWithZeros(flux,outputInitialized);
    if (externalFields != nullptr && !externalFieldsPrepared) {
        WVRealFieldBundleConstView ignored;
        auto status = ensurePhysicalFields(state,ignored,externalFields,externalFieldsPrepared);
        if (!status) return status;
    }
    if (context != nullptr) {
        context->owner_ = this;
        context->advectionFields_ = {externalFields->data,externalFields->shape};
        context->generation_ = evaluationGeneration_;
    }
    ++metrics_.evaluationCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationForcingEngine::advectFGridScalar(
    const WVConstantStratificationRightHandSideContext& context,
    const WVRealVolumeConstView& scalar, bool shouldAntialias,
    WVRealVolumeView& rightHandSide) {
    if (context.owner_ != this || context.generation_ != evaluationGeneration_ ||
        context.advectionFields_.data == nullptr)
        return {WVKernelStatusCode::invalidConfiguration,"The RHS evaluation context is stale or belongs to another forcing engine."};
    return kernel_->advectFGridScalar(scalar,context.advectionFields_,shouldAntialias,rightHandSide);
}

WVStateConstraintResult WVConstantStratificationForcingEngine::restoreForcingAmplitudes(WVMutableCoefficients& coefficients) {
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

void WVConstantStratificationForcingEngine::clearEvaluationWorkspace() noexcept {
    physicalFieldsValid_ = false;
}

std::size_t WVConstantStratificationForcingEngine::persistentBytes() const noexcept {
    return sizeof(*this)+kernel_->persistentBytes()+metrics_.scheduleBytes+metrics_.derivedOperatorBytes+metrics_.workspaceCapacityBytes;
}

} // namespace wavevortex::runtime
