#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"
#include "WVCoefficientFormulas.hpp"
#include "WVPreparedModeExecutor.hpp"
#include "WVHorizontalDerivativeMultiplier.hpp"
#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include <array>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <chrono>
#include <limits>
#include <new>
#include <system_error>
#include <stdexcept>

namespace wavevortex {
namespace kernel_detail {
struct WVCompactHorizontalGroup {
    std::unique_ptr<WVRetainedHorizontalOperator> operation;
    std::unique_ptr<WVRetainedHorizontalWorkspace> workspace;
    std::size_t realBytes = 0, complexBytes = 0;
    std::size_t persistentBytes() const noexcept {
        return sizeof(*this) + operation->persistentBytes() + workspace->persistentBytes();
    }
};
struct WVCompactConstantSchedule {
    WVHalfSpectrumMappings mapping;
    WVConstantKernelExecutionOptions options;
    std::array<std::shared_ptr<WVCompactHorizontalGroup>,4> horizontal;
    std::vector<std::unique_ptr<WVFFTPlan>> scalarPlans;
    std::unique_ptr<WVPreparedModeExecutor> pointwise;
    std::size_t persistentBytes() const noexcept {
        std::size_t bytes = sizeof(*this) - sizeof(mapping) + mapping.persistentBytes() +
            scalarPlans.capacity()*sizeof(std::unique_ptr<WVFFTPlan>);
        if (pointwise) bytes += pointwise->persistentBytes();
        for (std::size_t i=0;i<horizontal.size();++i) if (horizontal[i]) {
            bytes += horizontal[i]->persistentBytes();
            const auto identity=horizontal[i]->workspace->sharedResourceIdentity();
            if (identity) for (std::size_t previous=0;previous<i;++previous)
                if (horizontal[previous] && horizontal[previous]->workspace->sharedResourceIdentity()==identity) {
                    // A provider may share its guarded execution resources
                    // across the different channel layouts of this kernel.
                    bytes -= horizontal[i]->workspace->sharedResourceBytes();
                    break;
                }
        }
        return bytes;
    }
};
}
namespace {
class WVCompactHorizontalPlan final : public WVFFTPlan {
public:
    WVCompactHorizontalPlan(std::shared_ptr<kernel_detail::WVCompactHorizontalGroup> group, bool forward)
        : group_(std::move(group)), forward_(forward) {}
    WVKernelStatus execute(const void* input, void* output) override {
        if (forward_) return group_->operation->forward(*group_->workspace,
            {static_cast<const double*>(input),group_->realBytes},
            {static_cast<WVComplex64*>(output),nullptr,nullptr,group_->complexBytes});
        return group_->operation->inverse(*group_->workspace,
            {static_cast<const WVComplex64*>(input),nullptr,nullptr,group_->complexBytes},
            {static_cast<double*>(output),group_->realBytes});
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
private:
    std::shared_ptr<kernel_detail::WVCompactHorizontalGroup> group_;
    bool forward_;
};

constexpr double pi = 3.141592653589793238462643383279502884;

#ifndef WV_KERNEL_COEFFICIENT_WORKERS
#define WV_KERNEL_COEFFICIENT_WORKERS 2
#endif

#if defined(WV_KERNEL_NATIVE_OPTIMIZATION) && WV_KERNEL_NATIVE_OPTIMIZATION && (defined(__clang__) || defined(__GNUC__))
#define WV_KERNEL_RESTRICT __restrict__
#else
#define WV_KERNEL_RESTRICT
#endif

enum PlanIndex : std::size_t {
    horizontalForward3, horizontalForward4, horizontalInverse2, horizontalInverse3, horizontalInverse4,
    verticalDCT2Storage2,
    verticalDCT2Storage3, verticalDST1Storage3, verticalDST2Storage3, verticalDCT1Storage3,
    verticalDCT2Storage4, verticalDST2Storage4,
    verticalDCT3Storage4, verticalDST3Storage4,
    verticalDCT1Storage4, verticalDST1Storage4,
    horizontalForward1, verticalDCT1Storage1, verticalDST1Storage1,
    planCount
};
enum MatlabPlanIndex : std::size_t {
    matlabRetainedDCT, matlabRetainedDST, matlabColumnDCT, matlabColumnDST,
    matlabPlanCount
};

using detail::EvolvedWaveVortexCoefficients;
using detail::add;
using detail::atReferenceTime;
using detail::buoyancyProjection;
using detail::coefficientValueForField;
using detail::conjugate;
using detail::evolveWaveVortexCoefficients;
using detail::geostrophicCoefficient;
using detail::inertialCoefficientPair;
using detail::multiply;
using detail::phase;
using detail::subtract;
using detail::waveCoefficientPair;

std::size_t checkedProduct(std::size_t first, std::size_t second) {
    if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) throw std::overflow_error("kernel scratch size overflow");
    return first * second;
}

class ExecutionGuard {
public:
    explicit ExecutionGuard(bool& flag) : flag_(flag), entered_(!flag) { if (entered_) flag_ = true; }
    ~ExecutionGuard() { if (entered_) flag_ = false; }
    bool entered() const noexcept { return entered_; }
private:
    bool& flag_;
    bool entered_;
};

WVKernelStatus validateBundle(const WVRealFieldBundleConstView& view, WVShape3D spatial, std::size_t channels, const char* name) {
    if (view.shape.first != spatial.first || view.shape.second != spatial.second || view.shape.third != spatial.third || view.shape.fourth != channels) {
        return {WVKernelStatusCode::invalidShape, std::string(name) + " must have shape [Nx,Ny,Nz,Nfield]."};
    }
    if (view.data == nullptr) return {WVKernelStatusCode::invalidPointer, std::string(name) + " has a null data pointer."};
    return WVKernelStatus::ok();
}

WVKernelStatus validateBundle(const WVRealFieldBundleView& view, WVShape3D spatial, std::size_t channels, const char* name) {
    return validateBundle(WVRealFieldBundleConstView{view.data, view.shape}, spatial, channels, name);
}

WVKernelStatus validateVolume(const WVRealVolumeConstView& view, WVShape3D shape,
    const char* name) {
    if (view.shape.first!=shape.first || view.shape.second!=shape.second ||
        view.shape.third!=shape.third)
        return {WVKernelStatusCode::invalidShape,std::string(name)+" must have shape [Nx,Ny,Nz]."};
    if (!view.data || reinterpret_cast<std::uintptr_t>(view.data)%alignof(double))
        return {WVKernelStatusCode::invalidPointer,std::string(name)+" has invalid storage."};
    return WVKernelStatus::ok();
}
WVKernelStatus validateVolume(const WVRealVolumeView& view, WVShape3D shape,
    const char* name) {
    return validateVolume(WVRealVolumeConstView{view.data,view.shape},shape,name);
}
template <typename T>
WVKernelStatus validateComplexShape(const WVMatrixView<T>& view,WVShape2D shape,
    const char* name) {
    if (view.shape.rows!=shape.rows || view.shape.columns!=shape.columns)
        return {WVKernelStatusCode::invalidShape,std::string(name)+" has the wrong matrix shape."};
    if (!view.data || reinterpret_cast<std::uintptr_t>(view.data)%alignof(WVComplex64))
        return {WVKernelStatusCode::invalidPointer,std::string(name)+" has invalid storage."};
    return WVKernelStatus::ok();
}

template <typename T>
WVKernelStatus validateSpectral(const WVMatrixView<T>& view, WVShape2D spectral, const char* name) {
    if (view.shape.rows != spectral.rows || view.shape.columns != spectral.columns) return {WVKernelStatusCode::invalidShape, std::string(name) + " must have shape [Nj,Nkl]."};
    if (view.data == nullptr) return {WVKernelStatusCode::invalidPointer, std::string(name) + " has a null data pointer."};
    return WVKernelStatus::ok();
}

bool memoryOverlaps(const void* first, std::size_t firstBytes, const void* second, std::size_t secondBytes) {
    if (first == nullptr || second == nullptr || firstBytes == 0 || secondBytes == 0) return false;
    const auto firstAddress = reinterpret_cast<std::uintptr_t>(first);
    const auto secondAddress = reinterpret_cast<std::uintptr_t>(second);
    return firstAddress < secondAddress + secondBytes && secondAddress < firstAddress + firstBytes;
}

WVKernelStatus validateForwardOwnership(const WVRealFieldBundleConstView& fields, const WVMutableCoefficients& coefficients, WVShape2D spectral) {
    const auto coefficientBytes = spectral.elementCount() * sizeof(WVComplex64);
    const auto fieldBytes = fields.shape.elementCount() * sizeof(double);
    const void* outputs[] = {coefficients.Ap.data,coefficients.Am.data,coefficients.A0.data};
    for (std::size_t first = 0; first < 3; ++first) {
        if (memoryOverlaps(outputs[first],coefficientBytes,fields.data,fieldBytes)) return {WVKernelStatusCode::overlappingArrays,"Forward coefficient outputs must not overlap the spatial input."};
        for (std::size_t second = first + 1; second < 3; ++second) if (memoryOverlaps(outputs[first],coefficientBytes,outputs[second],coefficientBytes)) return {WVKernelStatusCode::overlappingArrays,"Forward coefficient outputs must not overlap each other."};
    }
    return WVKernelStatus::ok();
}

WVKernelStatus validateInverseOwnership(const WVState& state, const WVRealFieldBundleView& fields, WVShape2D spectral) {
    const auto coefficientBytes = spectral.elementCount() * sizeof(WVComplex64);
    const auto fieldBytes = fields.shape.elementCount() * sizeof(double);
    for (const void* input : {static_cast<const void*>(state.coefficients.Ap.data),static_cast<const void*>(state.coefficients.Am.data),static_cast<const void*>(state.coefficients.A0.data)}) if (memoryOverlaps(fields.data,fieldBytes,input,coefficientBytes)) return {WVKernelStatusCode::overlappingArrays,"Inverse spatial output must not overlap coefficient inputs."};
    return WVKernelStatus::ok();
}

WVKernelStatus validateDerivativeOwnership(const WVComplexConstView& Apm, const WVComplexConstView& A0, const WVRealFieldBundleView& fields, WVShape2D spectral) {
    const auto coefficientBytes = spectral.elementCount() * sizeof(WVComplex64);
    const auto fieldBytes = fields.shape.elementCount() * sizeof(double);
    if (memoryOverlaps(fields.data,fieldBytes,Apm.data,coefficientBytes) || memoryOverlaps(fields.data,fieldBytes,A0.data,coefficientBytes)) return {WVKernelStatusCode::overlappingArrays,"Derivative output must not overlap coefficient inputs."};
    return WVKernelStatus::ok();
}

WVFFTPlanSpecification horizontalSpecification(const WVTransformConstantStratificationConfiguration& c, std::size_t fields, bool forward) {
    const std::size_t NxHalf = c.Nx / 2 + 1;
    const std::size_t realPlane = c.Nx * c.Ny;
    const std::size_t realField = realPlane * c.Nz;
    const std::size_t halfRows = NxHalf * c.Ny;
    const std::size_t halfField = halfRows * c.Nz;
    WVFFTPlanSpecification specification;
    specification.kind = forward ? WVFFTPlanKind::horizontalRealToComplex2D : WVFFTPlanKind::horizontalComplexToReal2D;
    specification.transformDimensions = {
        {c.Ny, static_cast<std::ptrdiff_t>(c.Nx), static_cast<std::ptrdiff_t>(c.Nz * fields * NxHalf)},
        {c.Nx, 1, static_cast<std::ptrdiff_t>(c.Nz * fields)}};
    specification.batchDimensions = {
        {c.Nz, static_cast<std::ptrdiff_t>(realPlane), 1},
        {fields, static_cast<std::ptrdiff_t>(realField), static_cast<std::ptrdiff_t>(c.Nz)}};
    specification.inputBytes = forward ? fields * realField * sizeof(double) : fields * halfField * sizeof(WVComplex64);
    specification.outputBytes = forward ? fields * halfField * sizeof(WVComplex64) : fields * realField * sizeof(double);
    specification.destroysInput = !forward;
    if (!forward) {
        for (auto& dimension : specification.transformDimensions) std::swap(dimension.inputStride, dimension.outputStride);
        for (auto& dimension : specification.batchDimensions) std::swap(dimension.inputStride, dimension.outputStride);
    }
    return specification;
}

WVFFTPlanSpecification verticalSpecification(const WVTransformConstantStratificationConfiguration& c, std::size_t halfRows, std::size_t storageChannels, std::size_t transformedChannels, bool sine) {
    WVFFTPlanSpecification specification;
    specification.kind = sine ? WVFFTPlanKind::verticalDSTI : WVFFTPlanKind::verticalDCTI;
    const std::size_t transformLength = sine ? c.Nz - 2 : c.Nz;
    specification.transformDimensions = {{transformLength, 2, 2}};
    specification.batchDimensions = {
        {2, 1, 1},
        {halfRows, static_cast<std::ptrdiff_t>(2 * c.Nz * storageChannels), static_cast<std::ptrdiff_t>(2 * c.Nz * storageChannels)},
        {transformedChannels, static_cast<std::ptrdiff_t>(2 * c.Nz), static_cast<std::ptrdiff_t>(2 * c.Nz)}};
    specification.inputBytes = 2 * c.Nz * halfRows * storageChannels * sizeof(double);
    specification.outputBytes = specification.inputBytes;
    specification.inPlace = true;
    return specification;
}

void normalizeForwardDCT(WVComplex64* data, std::size_t Nz, std::size_t rows, std::size_t storageChannels, std::size_t firstChannel, std::size_t channels) {
    const double scale = 1.0 / static_cast<double>(Nz - 1);
    for (std::size_t row = 0; row < rows; ++row) for (std::size_t channel = firstChannel; channel < firstChannel + channels; ++channel) {
        const auto base = Nz * channel + Nz * storageChannels * row;
        for (std::size_t z = 0; z < Nz; ++z) data[base + z] = multiply(data[base + z],scale);
        data[base + Nz - 1] = multiply(data[base + Nz - 1],0.5);
    }
}

void normalizeForwardDST(WVComplex64* data, std::size_t Nz, std::size_t rows, std::size_t storageChannels, std::size_t firstChannel, std::size_t channels) {
    const double scale = 1.0 / static_cast<double>(Nz - 1);
    for (std::size_t row = 0; row < rows; ++row) for (std::size_t channel = firstChannel; channel < firstChannel + channels; ++channel) {
            const auto base = Nz * channel + Nz * storageChannels * row;
            data[base] = {};
            for (std::size_t z = 1; z + 1 < Nz; ++z) data[base + z] = multiply(data[base + z], scale);
            data[base + Nz - 1] = {};
    }
}

void normalizeInverseDCT(WVComplex64* data, std::size_t Nz, std::size_t rows, std::size_t storageChannels, std::size_t firstChannel, std::size_t channels) {
    for (std::size_t row = 0; row < rows; ++row) for (std::size_t channel = firstChannel; channel < firstChannel + channels; ++channel) {
        const auto base = Nz * channel + Nz * storageChannels * row;
        for (std::size_t z = 0; z < Nz; ++z) data[base + z] = multiply(data[base + z],0.5);
    }
}

void normalizeInverseDST(WVComplex64* data, std::size_t Nz, std::size_t rows, std::size_t storageChannels, std::size_t firstChannel, std::size_t channels) {
    for (std::size_t row = 0; row < rows; ++row) for (std::size_t channel = firstChannel; channel < firstChannel + channels; ++channel) {
            const auto base = Nz * channel + Nz * storageChannels * row;
            data[base] = {};
            for (std::size_t z = 1; z + 1 < Nz; ++z) data[base + z] = multiply(data[base + z], 0.5);
            data[base + Nz - 1] = {};
    }
}

template <bool Conjugated>
WVComplex64 storedValue(const WVComplex64* half, const WVHalfSpectrumMappings& mapping, std::size_t Nz, std::size_t channels, std::size_t field, std::size_t mode, std::size_t z) {
    auto value = half[z + Nz * field + Nz * channels * mapping.storageRowsByWVIndex[mode]];
    if constexpr (Conjugated) value = conjugate(value);
    return value;
}

template <bool Conjugated, typename CoefficientSource>
void assembleFieldSpectraForMode(
    WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    CoefficientSource&& source) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto coefficients = source(index);
        WVComplex64 values[] = {
            coefficientValueForField<0>(modes,index,coefficients),
            coefficientValueForField<1>(modes,index,coefficients),
            coefficientValueForField<2>(modes,index,coefficients),
            coefficientValueForField<3>(modes,index,coefficients)};
        for (std::size_t target = 0; target < 4; ++target) {
            // FFTW's DCT-I/DST-I inverse is unnormalized. By linearity, applying
            // the required factor here avoids rewriting the complete dense
            // half spectrum after the vertical transforms.
            values[target] = target >= 2 && (j == 0 || j + 1 == Nz) ? WVComplex64{} : multiply(values[target],0.5);
            half[j + Nz * target + Nz * 4 * storageRow] = Conjugated ? conjugate(values[target]) : values[target];
        }
    }
}

template <std::size_t Target, bool Conjugated, typename CoefficientSource>
void assembleDerivativeSpectraForMode(
    WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    const std::vector<WVFourierMode>& horizontalModes,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    const std::array<bool,3>* derivativeMask,
    CoefficientSource&& source) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    const auto& horizontal = horizontalModes[mode];
    constexpr bool cosine = Target < 2;
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto value = coefficientValueForField<Target>(modes,index,source(index));
        if (!derivativeMask) {
            const auto derivatives = detail::normalizedDerivativeSpectrum<cosine,Conjugated>(value,
                horizontal.k,horizontal.l,modes.verticalWavenumber[j],j == 0 || j + 1 == Nz);
            half[j + Nz * 0 + Nz * 3 * storageRow] = derivatives.x;
            half[j + Nz * 1 + Nz * 3 * storageRow] = derivatives.y;
            half[j + Nz * 2 + Nz * 3 * storageRow] = derivatives.z;
            continue;
        }
        const auto oriented=Conjugated ? conjugate(value) : value;
        if ((*derivativeMask)[0]) {
            auto derivative=multiply(oriented,WVComplex64{0.0,Conjugated ? -horizontal.k : horizontal.k});
            if constexpr (cosine) derivative=multiply(derivative,0.5);
            else if (j == 0 || j + 1 == Nz) derivative={};
            else derivative=multiply(derivative,0.5);
            half[j + Nz * 0 + Nz * 3 * storageRow]=derivative;
        }
        if ((*derivativeMask)[1]) {
            auto derivative=multiply(oriented,WVComplex64{0.0,Conjugated ? -horizontal.l : horizontal.l});
            if constexpr (cosine) derivative=multiply(derivative,0.5);
            else if (j == 0 || j + 1 == Nz) derivative={};
            else derivative=multiply(derivative,0.5);
            half[j + Nz * 1 + Nz * 3 * storageRow]=derivative;
        }
        if ((*derivativeMask)[2]) {
            auto derivative=multiply(oriented,cosine ? -modes.verticalWavenumber[j] : modes.verticalWavenumber[j]);
            if constexpr (cosine) derivative=(j == 0 || j + 1 == Nz) ? WVComplex64{} : multiply(derivative,0.5);
            else derivative=multiply(derivative,0.5);
            half[j + Nz * 2 + Nz * 3 * storageRow]=derivative;
        }
    }
}

template <bool Conjugated, typename CoefficientSource>
void assembleVelocitySpectraForMode(
    WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    CoefficientSource&& source) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto coefficients = source(index);
        WVComplex64 values[] = {
            coefficientValueForField<0>(modes,index,coefficients),
            coefficientValueForField<1>(modes,index,coefficients),
            coefficientValueForField<2>(modes,index,coefficients)};
        values[0] = multiply(values[0],0.5);
        values[1] = multiply(values[1],0.5);
        values[2] = j == 0 || j + 1 == Nz ? WVComplex64{} : multiply(values[2],0.5);
        for (std::size_t target = 0; target < 3; ++target) half[j + Nz * target + Nz * 3 * storageRow] = Conjugated ? conjugate(values[target]) : values[target];
    }
}

template <bool Conjugated, typename CoefficientSource>
void assembleHorizontalVelocitySpectraForMode(
    WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    CoefficientSource&& source) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto coefficients = source(index);
        WVComplex64 values[] = {
            multiply(coefficientValueForField<0>(modes,index,coefficients),0.5),
            multiply(coefficientValueForField<1>(modes,index,coefficients),0.5)};
        for (std::size_t target = 0; target < 2; ++target)
            half[j + Nz * target + Nz * 2 * storageRow] = Conjugated ? conjugate(values[target]) : values[target];
    }
}

template <bool PhaseProvided>
WVComplex64 phaseForIndex(const WVConstantStratificationModes& modes, WVComplexConstView phaseValues, std::size_t index, double elapsed) {
    if constexpr (PhaseProvided) return phaseValues.data[index];
    return phase(modes.omega[index] * elapsed);
}

template <bool Conjugated, bool PhaseProvided, bool InertialMode>
void projectHydrostaticMode(
    const WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    const WVFourierMode& horizontal,
    WVComplexConstView phaseValues,
    double elapsed,
    double horizontalScale,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    WVComplex64* outputAp,
    WVComplex64* outputAm,
    WVComplex64* outputA0) {
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto U = multiply(storedValue<Conjugated>(half,mapping,Nz,3,0,mode,j),horizontalScale);
        const auto V = multiply(storedValue<Conjugated>(half,mapping,Nz,3,1,mode,j),horizontalScale);
        const auto N = multiply(storedValue<Conjugated>(half,mapping,Nz,3,2,mode,j),horizontalScale);
        const auto A0 = geostrophicCoefficient(U,V,N,horizontal.k,horizontal.l,modes.A0FromVorticity[index],modes.A0FromBuoyancy[index]);
        detail::WaveCoefficientPair coefficients;
        if constexpr (InertialMode) {
            coefficients = inertialCoefficientPair(U,V,modes.inertialScale[j]);
        } else {
            const auto divergence = add(multiply(U,WVComplex64{0.0,horizontal.k}),multiply(V,WVComplex64{0.0,horizontal.l}));
            const auto waveContribution = multiply(modes.ApmDProjection[index],divergence);
            coefficients = waveCoefficientPair(waveContribution,buoyancyProjection(N,A0,modes.NA0Field[index],modes.ApmNProjection[index]));
        }
        const auto referenceCoefficients = atReferenceTime(coefficients,phaseForIndex<PhaseProvided>(modes,phaseValues,index,elapsed));
        outputAp[index] = referenceCoefficients.Ap;
        outputAm[index] = referenceCoefficients.Am;
        outputA0[index] = A0;
    }
}

template <bool Conjugated, bool PhaseProvided, bool InertialMode>
void projectNonhydrostaticMode(
    const WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    const WVFourierMode& horizontal,
    WVComplexConstView phaseValues,
    double elapsed,
    double horizontalScale,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    WVComplex64* outputAp,
    WVComplex64* outputAm,
    WVComplex64* outputA0) {
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto U = multiply(storedValue<Conjugated>(half,mapping,Nz,4,0,mode,j),horizontalScale);
        const auto V = multiply(storedValue<Conjugated>(half,mapping,Nz,4,1,mode,j),horizontalScale);
        const auto W = multiply(storedValue<Conjugated>(half,mapping,Nz,4,2,mode,j),horizontalScale);
        const auto N = multiply(storedValue<Conjugated>(half,mapping,Nz,4,3,mode,j),horizontalScale);
        const auto A0 = geostrophicCoefficient(U,V,N,horizontal.k,horizontal.l,modes.A0FromVorticity[index],modes.A0FromBuoyancy[index]);
        detail::WaveCoefficientPair coefficients;
        if constexpr (InertialMode) {
            coefficients = inertialCoefficientPair(U,V,modes.inertialScale[j]);
        } else {
            const auto delta = multiply(add(multiply(U,horizontal.cosAlpha),multiply(V,horizontal.sinAlpha)),modes.ApmDScaled[index]);
            const auto wBar = multiply(WVComplex64{0.0,(horizontal.Kh / 2.0) * modes.apmWProjectionPrefactor[j]},W);
            coefficients = waveCoefficientPair(add(delta,wBar),buoyancyProjection(N,A0,modes.NA0Field[index],modes.ApmNProjection[index]));
        }
        const auto referenceCoefficients = atReferenceTime(coefficients,phaseForIndex<PhaseProvided>(modes,phaseValues,index,elapsed));
        outputAp[index] = referenceCoefficients.Ap;
        outputAm[index] = referenceCoefficients.Am;
        outputA0[index] = A0;
    }
}

template <bool FFamily, bool Conjugated>
void assembleFieldFamilyDerivativesForMode(
    WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    const WVFourierMode& horizontal,
    const WVComplexConstView& Apm,
    const WVComplexConstView& A0,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const double waveScale = FFamily ? modes.fWaveScale[index] : modes.gWaveScale[j];
        const double geostrophicScale = FFamily ? modes.Fg[j] : modes.Gg[j];
        const auto values = FFamily ? detail::fieldFamilyDerivativeSpectrum<true,Conjugated>(Apm.data[index],A0.data[index],waveScale,geostrophicScale,horizontal.k,horizontal.l,modes.verticalWavenumber[j]) : detail::fieldFamilyDerivativeSpectrum<false,Conjugated>(Apm.data[index],A0.data[index],waveScale,geostrophicScale,horizontal.k,horizontal.l,modes.verticalWavenumber[j]);
        half[j + Nz * 0 + Nz * 4 * storageRow] = values.value;
        half[j + Nz * 1 + Nz * 4 * storageRow] = values.x;
        half[j + Nz * 2 + Nz * 4 * storageRow] = values.y;
        half[j + Nz * 3 + Nz * 4 * storageRow] = values.z;
    }
}

template <std::size_t Target, bool Hydrostatic, bool Conjugated, bool InertialMode>
void projectFluxTargetMode(
    const WVComplex64* half,
    const WVHalfSpectrumMappings& mapping,
    const WVConstantStratificationModes& modes,
    const WVFourierMode& horizontal,
    WVComplexConstView phaseValues,
    double horizontalScale,
    std::size_t Nz,
    std::size_t Nj,
    std::size_t mode,
    WVFlux& flux) {
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto value = multiply(storedValue<Conjugated>(half,mapping,Nz,1,0,mode,j),horizontalScale);
        WVComplex64 A0Contribution{};
        if constexpr (Target == 0) A0Contribution = multiply(multiply(value,WVComplex64{0.0,-horizontal.l}),modes.A0FromVorticity[index]);
        if constexpr (Target == 1) A0Contribution = multiply(multiply(value,WVComplex64{0.0,horizontal.k}),modes.A0FromVorticity[index]);
        if constexpr (Target == 3) A0Contribution = multiply(value,modes.A0FromBuoyancy[index]);
        const auto buoyancyValue = subtract(Target == 3 ? value : WVComplex64{},multiply(A0Contribution,modes.NA0Field[index]));
        const auto buoyancyContribution = multiply(buoyancyValue,modes.ApmNProjection[index]);
        WVComplex64 waveContribution{};
        if constexpr (Hydrostatic) {
            if constexpr (Target == 0) waveContribution = multiply(multiply(value,WVComplex64{0.0,horizontal.k}),modes.ApmDProjection[index]);
            if constexpr (Target == 1) waveContribution = multiply(multiply(value,WVComplex64{0.0,horizontal.l}),modes.ApmDProjection[index]);
        } else {
            if constexpr (Target == 0) waveContribution = multiply(value,horizontal.cosAlpha * modes.ApmDScaled[index]);
            if constexpr (Target == 1) waveContribution = multiply(value,horizontal.sinAlpha * modes.ApmDScaled[index]);
            if constexpr (Target == 2) waveContribution = multiply(value,WVComplex64{0.0,(horizontal.Kh / 2.0) * modes.apmWProjectionPrefactor[j]});
        }
        auto contribution = waveCoefficientPair(waveContribution,buoyancyContribution);
        if constexpr (InertialMode) contribution = inertialCoefficientPair<Target>(value,modes.inertialScale[j]);
        detail::accumulateAtReferenceTime(flux.Fp.data[index],flux.Fm.data[index],contribution,phaseValues.data[index]);
        flux.F0.data[index] = add(flux.F0.data[index],A0Contribution);
    }
}

void completeHermitianBoundaries(WVComplex64* half, const WVHalfSpectrumMappings& mapping, std::size_t Nz, std::size_t channels) {
    for (std::size_t field = 0; field < channels; ++field) {
        for (std::size_t i = 0; i < mapping.hermitianCompletionRows.size(); ++i) for (std::size_t z = 0; z < Nz; ++z) half[z + Nz * field + Nz * channels * mapping.hermitianCompletionRows[i]] = conjugate(half[z + Nz * field + Nz * channels * mapping.hermitianSourceRows[i]]);
        for (const auto row : mapping.selfConjugateRows) for (std::size_t z = 0; z < Nz; ++z) half[z + Nz * field + Nz * channels * row].imag = 0.0;
    }
}

} // namespace

WVTransformConstantStratificationKernel::WVTransformConstantStratificationKernel() = default;
WVTransformConstantStratificationKernel::~WVTransformConstantStratificationKernel() = default;

WVKernelStatus WVTransformConstantStratificationKernel::create(
    const WVTransformConstantStratificationConfiguration& configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVTransformConstantStratificationKernel>& kernel, WVConstantKernelExecutionOptions options) {
    if (!engine) return {WVKernelStatusCode::invalidPointer, "An FFT engine is required."};
    if (options.schedule != WVConstantNonlinearFluxSchedule::frozenStreamed && options.schedule != WVConstantNonlinearFluxSchedule::compactCandidate)
        return {WVKernelStatusCode::invalidConfiguration,"Unknown constant nonlinear-flux schedule."};
    if (options.schedule == WVConstantNonlinearFluxSchedule::compactCandidate && (!options.horizontalOuterWorkers || !options.pointwiseWorkers))
        return {WVKernelStatusCode::invalidConfiguration,"Compact worker counts must be positive."};
    try {
        auto candidate = std::unique_ptr<WVTransformConstantStratificationKernel>(new WVTransformConstantStratificationKernel());
        auto status = WVTransformConstantStratificationDescriptor::create(configuration, candidate->descriptor_);
        if (!status) return status;
        candidate->engineIdentifier_ = engine->identifier();
        candidate->engineLibraryIdentity_ = engine->libraryIdentity();
        candidate->engine_ = std::shared_ptr<WVFFTEngine>(std::move(engine));
        if (options.schedule == WVConstantNonlinearFluxSchedule::compactCandidate) {
            candidate->compact_ = std::make_unique<kernel_detail::WVCompactConstantSchedule>();
            auto& compact = *candidate->compact_;
            compact.options = options;
            const auto& original = candidate->descriptor_.halfSpectrumMappings();
            // Keep one original stored orientation per Hermitian orbit. The
            // retained horizontal inverse completes its boundary partner after
            // the real vertical transform; that completion commutes exactly.
            compact.mapping.NxHalf = original.NxHalf;
            compact.mapping.conjugatesStoredValueByWVIndex = original.conjugatesStoredValueByWVIndex;
            for (std::size_t mode=0; mode<candidate->descriptor_.Nkl(); ++mode) {
                compact.mapping.storageRowsByWVIndex.push_back(mode);
                if (std::find(original.selfConjugateRows.begin(),original.selfConjugateRows.end(),original.storageRowsByWVIndex[mode]) != original.selfConjugateRows.end())
                    compact.mapping.selfConjugateRows.push_back(mode);
            }
            compact.pointwise = std::make_unique<kernel_detail::WVPreparedModeExecutor>(options.pointwiseWorkers);
        }
        constexpr std::size_t halfChannels = 4;
        const auto halfElements = checkedProduct(candidate->executionRowCount(), checkedProduct(configuration.Nz, halfChannels));
        const auto halfFieldElements = checkedProduct(checkedProduct(candidate->descriptor_.halfSpectrumMappings().NxHalf,configuration.Ny),configuration.Nz);
        if (candidate->descriptor_.spectralShape().elementCount() > halfFieldElements) return {WVKernelStatusCode::invalidShape,"The streamed phase reservation requires M <= H."};
        constexpr std::size_t realChannels = 6;
        const auto realElements = checkedProduct(candidate->descriptor_.spatialShape().elementCount(),realChannels);
        candidate->halfSpectrumScratch_.resize(2 * halfElements);
        candidate->realScratch_.resize(realElements);
        candidate->preparedPhase_.resize(candidate->descriptor_.spectralShape().elementCount());
        candidate->metrics_.halfSpectrumScratchCapacityBytes = candidate->halfSpectrumScratch_.capacity() * sizeof(double);
        candidate->metrics_.realScratchCapacityBytes = candidate->realScratch_.capacity() * sizeof(double);
        candidate->metrics_.scratchCapacityBytes = candidate->metrics_.halfSpectrumScratchCapacityBytes + candidate->metrics_.realScratchCapacityBytes + candidate->preparedPhase_.capacity()*sizeof(WVComplex64);
        candidate->metrics_.scratchHighWaterBytes = candidate->scratchBytes();
        const auto halfRows = candidate->descriptor_.halfSpectrumMappings().NxHalf * configuration.Ny;
        candidate->scalarAntialiasRows_.assign(halfRows,0);
        const auto& mappings = candidate->descriptor_.halfSpectrumMappings();
        for (const auto row : mappings.storageRowsByWVIndex) candidate->scalarAntialiasRows_[row] = 1;
        for (const auto row : mappings.hermitianCompletionRows) candidate->scalarAntialiasRows_[row] = 1;
        for (const auto row : mappings.hermitianSourceRows) candidate->scalarAntialiasRows_[row] = 1;
        candidate->metrics_.descriptorBytes =
            candidate->descriptor_.persistentBytes() +
            candidate->scalarAntialiasRows_.capacity() * sizeof(std::uint8_t);
        status = candidate->preparePlans();
        if (!status) return status;
        candidate->coefficientExecutor_ = std::make_unique<kernel_detail::WVPreparedModeExecutor>(
            std::max<std::size_t>(1, std::min<std::size_t>(WV_KERNEL_COEFFICIENT_WORKERS, candidate->descriptor_.Nkl())));
        candidate->metrics_.engineBytes = candidate->engine_->persistentBytes();
        candidate->metrics_.kernelManagementBytes =
            sizeof(*candidate) - sizeof(candidate->descriptor_) +
            candidate->engineIdentifier_.capacity() +
            candidate->engineLibraryIdentity_.capacity() +
            candidate->plans_.capacity() *
                sizeof(std::unique_ptr<WVFFTPlan>) +
            candidate->coefficientExecutor_->persistentBytes();
        kernel = std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure, "Unable to allocate the bounded kernel scratch arena."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow, error.what()};
    } catch (const std::system_error& error) {
        return {WVKernelStatusCode::allocationFailure, error.what()};
    }
}

WVKernelStatus WVTransformConstantStratificationKernel::preparePlans() {
    plans_.resize(planCount);
    const auto& c = descriptor_.configuration();
    const auto halfRows = executionRowCount();
    const WVFFTPlanSpecification specifications[] = {
        horizontalSpecification(c, 3, true), horizontalSpecification(c, 4, true), horizontalSpecification(c, 2, false), horizontalSpecification(c, 3, false), horizontalSpecification(c, 4, false),
        verticalSpecification(c,halfRows,2,2,false),
        verticalSpecification(c,halfRows,3,2,false), verticalSpecification(c,halfRows,3,1,true), verticalSpecification(c,halfRows,3,2,true), verticalSpecification(c,halfRows,3,1,false),
        verticalSpecification(c,halfRows,4,2,false), verticalSpecification(c,halfRows,4,2,true),
        verticalSpecification(c,halfRows,4,3,false), verticalSpecification(c,halfRows,4,3,true),
        verticalSpecification(c,halfRows,4,1,false), verticalSpecification(c,halfRows,4,1,true)
        ,horizontalSpecification(c,1,true), verticalSpecification(c,halfRows,1,1,false), verticalSpecification(c,halfRows,1,1,true)
    };
    for (std::size_t i = 0; i < planCount; ++i) {
        if (i==horizontalInverse2 || i==verticalDCT2Storage2) continue;
        WVKernelStatus status;
        const bool horizontal = i == horizontalForward3 || i == horizontalForward4 ||
            i == horizontalInverse2 || i == horizontalInverse3 || i == horizontalInverse4 || i == horizontalForward1;
        if (compact_ && horizontal) {
            const auto channels = i == horizontalForward1 ? 1U :
                (i == horizontalInverse2 ? 2U : (i == horizontalForward4 || i == horizontalInverse4 ? 4U : 3U));
            const auto groupIndex = channels == 1 ? 0U : (channels == 3 ? 1U : (channels == 4 ? 2U : 3U));
            auto& group = compact_->horizontal[groupIndex];
            if (!group) {
                group = std::make_shared<kernel_detail::WVCompactHorizontalGroup>();
                WVRetainedHorizontalSpecification spec;
                const auto planes = c.Nz*channels;
                spec.grid = {c.Nx,c.Ny,planes,1,c.Nx,c.Nx*c.Ny,"constant-fields"};
                spec.retained = {planes,descriptor_.Nkl(),1,planes,WVComplexRepresentation::interleaved,"constant-fields","constant-stored-orbits"};
                spec.Lx=c.Lx; spec.Ly=c.Ly;
                spec.schedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
                spec.outerWorkers=compact_->options.horizontalOuterWorkers;
                const auto& original=descriptor_.halfSpectrumMappings();
                for (const auto row : original.storageRowsByWVIndex) {
                    const auto y=row/original.NxHalf;
                    const auto l=y<(c.Ny+1)/2 ? static_cast<std::int64_t>(y) : static_cast<std::int64_t>(y)-static_cast<std::int64_t>(c.Ny);
                    spec.modes.push_back({static_cast<std::int64_t>(row%original.NxHalf),l});
                }
                status=WVRetainedHorizontalOperator::createShared(spec,engine_,group->operation);
                if (!status) return status;
                status=group->operation->createWorkspace(group->workspace,false);
                if (!status) return status;
                group->realBytes=descriptor_.spatialShape().elementCount()*channels*sizeof(double);
                group->complexBytes=c.Nz*channels*descriptor_.Nkl()*sizeof(WVComplex64);
            }
            plans_[i]=std::make_unique<WVCompactHorizontalPlan>(group,i==horizontalForward1 || i==horizontalForward3 || i==horizontalForward4);
            status=WVKernelStatus::ok();
        } else status = engine_->createPlan(specifications[i], plans_[i]);
        if (!status || !plans_[i]) return status ? WVKernelStatus{WVKernelStatusCode::fftPlanFailure, "FFT engine returned an empty plan."} : status;
        metrics_.planBytes += plans_[i]->persistentBytes();
    }
    metrics_.planCount = plans_.size();
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::prepareHorizontalVelocityTransform() {
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"Prepare the horizontal-velocity transform before beginning a state evaluation."};
    if (plans_[horizontalInverse2] && plans_[verticalDCT2Storage2]) return WVKernelStatus::ok();
    try {
        const auto& c=descriptor_.configuration();
        std::unique_ptr<WVFFTPlan> horizontal;
        std::unique_ptr<WVFFTPlan> vertical;
        std::shared_ptr<kernel_detail::WVCompactHorizontalGroup> compactGroup;
        WVKernelStatus status;
        if (compact_) {
            compactGroup=std::make_shared<kernel_detail::WVCompactHorizontalGroup>();
            WVRetainedHorizontalSpecification specification;
            const auto planes=c.Nz*2;
            specification.grid={c.Nx,c.Ny,planes,1,c.Nx,c.Nx*c.Ny,"constant-horizontal-velocity"};
            specification.retained={planes,descriptor_.Nkl(),1,planes,WVComplexRepresentation::interleaved,
                "constant-horizontal-velocity","constant-stored-orbits"};
            specification.Lx=c.Lx;
            specification.Ly=c.Ly;
            specification.schedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
            specification.outerWorkers=compact_->options.horizontalOuterWorkers;
            const auto& original=descriptor_.halfSpectrumMappings();
            for (const auto row:original.storageRowsByWVIndex) {
                const auto y=row/original.NxHalf;
                const auto l=y<(c.Ny+1)/2 ? static_cast<std::int64_t>(y) :
                    static_cast<std::int64_t>(y)-static_cast<std::int64_t>(c.Ny);
                specification.modes.push_back({static_cast<std::int64_t>(row%original.NxHalf),l});
            }
            status=WVRetainedHorizontalOperator::createShared(specification,engine_,compactGroup->operation);
            if (!status) return status;
            status=compactGroup->operation->createWorkspace(compactGroup->workspace,false);
            if (!status) return status;
            compactGroup->realBytes=descriptor_.spatialShape().elementCount()*2*sizeof(double);
            compactGroup->complexBytes=c.Nz*2*descriptor_.Nkl()*sizeof(WVComplex64);
            horizontal=std::make_unique<WVCompactHorizontalPlan>(compactGroup,false);
        } else {
            status=engine_->createPlan(horizontalSpecification(c,2,false),horizontal);
            if (!status) return status;
        }
        status=engine_->createPlan(verticalSpecification(c,executionRowCount(),2,2,false),vertical);
        if (!status) return status;
        if (!horizontal || !vertical)
            return {WVKernelStatusCode::fftPlanFailure,"FFT engine returned an empty horizontal-velocity plan."};
        plans_[horizontalInverse2]=std::move(horizontal);
        plans_[verticalDCT2Storage2]=std::move(vertical);
        if (compact_) compact_->horizontal[3]=std::move(compactGroup);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Unable to prepare the horizontal-velocity transform."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    } catch (const std::system_error& error) {
        return {WVKernelStatusCode::allocationFailure,error.what()};
    }
}

WVKernelStatus WVTransformConstantStratificationKernel::prepareMatlabPrimitives() {
    ExecutionGuard guard(executing_);
    if (!guard.entered())
        return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,
            "Prepare MATLAB primitives before beginning a state evaluation."};
    if (matlabPrimitivesPrepared_) return WVKernelStatus::ok();
    try {
        auto status=ensureScalarPlans(); if (!status) return status;
        status=ensureScalarInversePlan(); if (!status) return status;
        const auto& c=descriptor_.configuration();
        std::array<std::unique_ptr<WVFFTPlan>,matlabPlanCount> prepared;
        if (!compact_) {
            status=engine_->createPlan(verticalSpecification(c,descriptor_.Nkl(),1,1,false),
                prepared[matlabRetainedDCT]); if (!status) return status;
            status=engine_->createPlan(verticalSpecification(c,descriptor_.Nkl(),1,1,true),
                prepared[matlabRetainedDST]); if (!status) return status;
        }
        status=engine_->createPlan(verticalSpecification(c,1,1,1,false),
            prepared[matlabColumnDCT]); if (!status) return status;
        status=engine_->createPlan(verticalSpecification(c,1,1,1,true),
            prepared[matlabColumnDST]); if (!status) return status;
        const bool missingColumn=!prepared[matlabColumnDCT] || !prepared[matlabColumnDST];
        const bool missingSchedulePlan=!compact_ &&
            (!prepared[matlabRetainedDCT] || !prepared[matlabRetainedDST]);
        if (missingColumn || missingSchedulePlan)
            return {WVKernelStatusCode::fftPlanFailure,
                "FFT engine returned an empty MATLAB primitive plan."};
        matlabPlans_=std::move(prepared);
        matlabPrimitivesPrepared_=true;
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,
            "Unable to prepare MATLAB primitive plans."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    } catch (const std::system_error& error) {
        return {WVKernelStatusCode::allocationFailure,error.what()};
    }
}

std::size_t WVTransformConstantStratificationKernel::persistentBytes() const noexcept {
    const auto& current = metrics();
    return current.kernelManagementBytes + current.engineBytes +
           current.descriptorBytes + current.planBytes +
           current.scratchCapacityBytes;
}

const WVKernelMetrics&
WVTransformConstantStratificationKernel::metrics() const noexcept {
    metrics_.engineBytes = engine_ == nullptr ? 0 : engine_->persistentBytes();
    metrics_.planBytes = 0;
    metrics_.planCount = 0;
    for (const auto& plan : plans_) {
        if (plan != nullptr) {
            metrics_.planBytes += plan->persistentBytes();
            ++metrics_.planCount;
        }
    }
    if (scalarInversePlan_ != nullptr) {
        metrics_.planBytes += scalarInversePlan_->persistentBytes();
        ++metrics_.planCount;
    }
    for (const auto& plan:matlabPlans_) if (plan) {
        metrics_.planBytes+=plan->persistentBytes();
        ++metrics_.planCount;
    }
    metrics_.kernelManagementBytes = sizeof(*this) - sizeof(descriptor_) +
        engineIdentifier_.capacity() + engineLibraryIdentity_.capacity() +
        plans_.capacity()*sizeof(std::unique_ptr<WVFFTPlan>) +
        (coefficientExecutor_ ? coefficientExecutor_->persistentBytes() : 0);
    if (compact_) {
        metrics_.kernelManagementBytes += compact_->persistentBytes();
        for (const auto& group : compact_->horizontal)
            if (group) metrics_.planBytes += group->workspace->planBytesLowerBound();
        for (const auto& plan : compact_->scalarPlans) if (plan) {
            metrics_.planBytes += plan->persistentBytes(); ++metrics_.planCount;
        }
    }
    metrics_.halfSpectrumScratchCapacityBytes=halfSpectrumScratch_.capacity()*sizeof(double);
    metrics_.scratchCapacityBytes=metrics_.halfSpectrumScratchCapacityBytes+metrics_.realScratchCapacityBytes+
        preparedPhase_.capacity()*sizeof(WVComplex64);
    metrics_.scratchHighWaterBytes=std::max(metrics_.scratchHighWaterBytes,scratchBytes());
    return metrics_;
}

const WVHalfSpectrumMappings& WVTransformConstantStratificationKernel::executionMapping() const noexcept {
    return compact_ ? compact_->mapping : descriptor_.halfSpectrumMappings();
}
std::size_t WVTransformConstantStratificationKernel::executionRowCount() const noexcept {
    return compact_ ? descriptor_.Nkl() : descriptor_.halfSpectrumMappings().NxHalf*descriptor_.configuration().Ny;
}
std::size_t WVTransformConstantStratificationKernel::verticalExecutionRowCount() const noexcept { return executionRowCount(); }
std::size_t WVTransformConstantStratificationKernel::pointwiseWorkerCount() const noexcept { return compact_ ? compact_->pointwise->workerCount() : 1; }
std::size_t WVTransformConstantStratificationKernel::retainedHorizontalWorkerCount() const noexcept {
    return compact_ ? compact_->horizontal[1]->workspace->workerCount() : 0;
}
const char* WVTransformConstantStratificationKernel::retainedHorizontalScheduleIdentifier() const noexcept {
    return compact_ ? compact_->horizontal[1]->workspace->scheduleIdentifier() : "full-half-spectrum";
}

const char* WVTransformConstantStratificationKernel::nonlinearFluxScheduleIdentifier() const noexcept {
    return compact_ ? "retained-compact-streamed-target-three-channel-v1" : "streamed-target-three-channel";
}

const char* WVTransformConstantStratificationKernel::phaseImplementationIdentifier() const noexcept {
    return compact_ ? "fused-scalar-sincos" : "scalar-sincos";
}

const char* WVTransformConstantStratificationKernel::coefficientArithmeticModeIdentifier() const noexcept {
    return "natural-dimensional-prescaled";
}

const char* WVTransformConstantStratificationKernel::inverseNormalizationPlacementIdentifier() const noexcept {
    return "coefficient-production";
}

const char* WVTransformConstantStratificationKernel::optimizationImplementationIdentifier() const noexcept {
#if defined(WV_KERNEL_NATIVE_OPTIMIZATION) && WV_KERNEL_NATIVE_OPTIMIZATION
    return "O3-native";
#else
    return "portable";
#endif
}

std::size_t WVTransformConstantStratificationKernel::coefficientWorkerCount() const noexcept { return coefficientExecutor_ ? coefficientExecutor_->workerCount() : 0; }

std::size_t WVTransformConstantStratificationKernel::phaseReservationBytes() const noexcept {
    const auto& c = descriptor_.configuration();
    return executionRowCount() * c.Nz * sizeof(WVComplex64);
}

void WVTransformConstantStratificationKernel::setStageInstrumentation(bool enabled) noexcept {
    stageInstrumentationEnabled_ = enabled;
    if (!enabled) return;
    metrics_.phaseSeconds = 0.0;
    metrics_.reconstructionSeconds = 0.0;
    metrics_.derivativeReconstructionSeconds = 0.0;
    metrics_.productSeconds = 0.0;
    metrics_.projectionSeconds = 0.0;
    metrics_.coefficientAssemblySeconds = 0.0;
    metrics_.derivativeCoefficientAssemblySeconds = 0.0;
    metrics_.coefficientProjectionSeconds = 0.0;
}

WVKernelStatus WVTransformConstantStratificationKernel::validateStateContents(
    const WVState& state,bool allowPreparedExecutor) {
    if (!std::isfinite(state.t) || !std::isfinite(state.t0) || !std::isfinite(state.t-state.t0))
        return {WVKernelStatusCode::invalidConfiguration,"State times and elapsed time must be finite."};
    const auto spectral=descriptor_.spectralShape();
    const WVKernelStatus statuses[]={validateSpectral(state.coefficients.Ap,spectral,"Ap"),
        validateSpectral(state.coefficients.Am,spectral,"Am"),validateSpectral(state.coefficients.A0,spectral,"A0")};
    for (const auto& status:statuses) if (!status) return status;
    const WVComplexConstView views[]={state.coefficients.Ap,state.coefficients.Am,
        state.coefficients.A0};
    const auto count=spectral.elementCount();
    std::atomic<bool> coefficientsFinite{true};
    const auto scanCoefficients=[&](std::size_t begin,std::size_t end) {
        while(begin<end) {
            const auto viewIndex=begin/count;
            const auto viewBegin=begin-viewIndex*count;
            const auto viewEnd=std::min(count,end-viewIndex*count);
            const auto* values=views[viewIndex].data;
            for(std::size_t i=viewBegin;i<viewEnd;++i)
                if(!std::isfinite(values[i].real) || !std::isfinite(values[i].imag)) {
                    coefficientsFinite.store(false,std::memory_order_relaxed);
                    return;
                }
            begin+=viewEnd-viewBegin;
        }
    };
    if(allowPreparedExecutor) coefficientExecutor_->execute(3*count,scanCoefficients);
    else scanCoefficients(0,3*count);
    if(!coefficientsFinite.load(std::memory_order_relaxed))
        return {WVKernelStatusCode::numericalFailure,"Nonfinite constant-stratification coefficient."};
    const double elapsed=state.t-state.t0;
    for (const auto frequency:descriptor_.verticalModes().omega)
        if (!std::isfinite(frequency*elapsed))
            return {WVKernelStatusCode::numericalFailure,"Constant-stratification phase overflow."};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::validateState(
    const WVState& state,bool allowPreparedExecutor) {
    ++metrics_.stateValidationCount;
    return validateStateContents(state,allowPreparedExecutor);
}

bool WVTransformConstantStratificationKernel::matchesStateEvaluation(const WVState& state) const noexcept {
    const auto sameView=[](WVComplexConstView a,WVComplexConstView b) {
        return a.data==b.data && a.shape.rows==b.shape.rows && a.shape.columns==b.shape.columns;
    };
    if (!stateEvaluationActive_ || state.t!=preparedState_.t || state.t0!=preparedState_.t0) return false;
    for (std::size_t i=0;i<preparedStateViewCount_;++i) if (
        sameView(state.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(state.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(state.coefficients.A0,preparedStateViews_[i].coefficients.A0)) return true;
    return false;
}

WVKernelStatus WVTransformConstantStratificationKernel::validateStateEvaluation(
    const WVState& state) const noexcept {
    if (!matchesStateEvaluation(state))
        return {WVKernelStatusCode::invalidConfiguration,
            "State does not belong to the active constant-stratification evaluation."};
    return WVKernelStatus::ok();
}

std::size_t WVTransformConstantStratificationKernel::stateEvaluationComponent(
    const WVState& state) const noexcept {
    const auto sameView=[](WVComplexConstView a,WVComplexConstView b) {
        return a.data==b.data && a.shape.rows==b.shape.rows && a.shape.columns==b.shape.columns;
    };
    if (!stateEvaluationActive_) return 0;
    for (std::size_t i=0;i<preparedStateViewCount_;++i) if (
        sameView(state.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(state.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(state.coefficients.A0,preparedStateViews_[i].coefficients.A0))
        return preparedStateComponents_[i];
    return 0;
}

WVKernelStatus WVTransformConstantStratificationKernel::validateStateForCall(const WVState& state) {
    if (!stateEvaluationActive_) return validateState(state);
    if (!matchesStateEvaluation(state))
        return {WVKernelStatusCode::invalidConfiguration,"State does not match the active constant-stratification evaluation."};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::validateStateAndFluxForCall(const WVState& state,const WVFlux& flux) {
    auto status=validateStateForCall(state); if (!status) return status;
    const auto spectral=descriptor_.spectralShape();
    const WVKernelStatus outputStatuses[]={validateSpectral(flux.Fp,spectral,"Fp"),
        validateSpectral(flux.Fm,spectral,"Fm"),validateSpectral(flux.F0,spectral,"F0")};
    for (const auto& value:outputStatuses) if (!value) return value;
    const auto bytes=spectral.elementCount()*sizeof(WVComplex64);
    const void* inputs[]={state.coefficients.Ap.data,state.coefficients.Am.data,state.coefficients.A0.data};
    const void* outputs[]={flux.Fp.data,flux.Fm.data,flux.F0.data};
    for (std::size_t i=0;i<3;++i) {
        for (std::size_t j=i+1;j<3;++j) if (memoryOverlaps(outputs[i],bytes,outputs[j],bytes))
            return {WVKernelStatusCode::overlappingArrays,"Flux outputs must not overlap."};
        for (const auto* input:inputs) if (memoryOverlaps(outputs[i],bytes,input,bytes))
            return {WVKernelStatusCode::overlappingArrays,"Flux outputs must not overlap state inputs."};
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::validateMutableOutputOutsidePreparedState(
    const WVMutableCoefficients& output) const {
    if (!stateEvaluationActive_) return WVKernelStatus::ok();
    const auto bytes=descriptor_.spectralShape().elementCount()*sizeof(WVComplex64);
    for (const auto target:{output.Ap,output.Am,output.A0}) for (std::size_t stateIndex=0;stateIndex<preparedStateViewCount_;++stateIndex)
        for (const auto input:{preparedStateViews_[stateIndex].coefficients.Ap,preparedStateViews_[stateIndex].coefficients.Am,preparedStateViews_[stateIndex].coefficients.A0})
            if (memoryOverlaps(target.data,bytes,input.data,bytes))
                return {WVKernelStatusCode::overlappingArrays,"Mutable coefficient output overlaps an active immutable constant-stratification state view."};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::validateMutableOutputOutsidePreparedState(
    const void* output,std::size_t outputBytes) const {
    if (!stateEvaluationActive_) return WVKernelStatus::ok();
    const auto bytes=descriptor_.spectralShape().elementCount()*sizeof(WVComplex64);
    for (std::size_t stateIndex=0;stateIndex<preparedStateViewCount_;++stateIndex)
        for (const auto input:{preparedStateViews_[stateIndex].coefficients.Ap,
            preparedStateViews_[stateIndex].coefficients.Am,
            preparedStateViews_[stateIndex].coefficients.A0})
            if (memoryOverlaps(output,outputBytes,input.data,bytes))
                return {WVKernelStatusCode::overlappingArrays,
                    "Mutable output overlaps an active immutable constant-stratification state view."};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::prepareStatePhase(const WVState& state) {
    const double elapsed=state.t-state.t0;
    if (!std::isfinite(state.t) || !std::isfinite(state.t0) || !std::isfinite(elapsed))
        return {WVKernelStatusCode::invalidConfiguration,"State times and elapsed time must be finite."};
    const auto& omega=descriptor_.verticalModes().omega;
    // Every registered active view was validated with these exact times and
    // immutable coefficients. Standalone and foreign-state phase requests keep
    // the full overflow scan.
    if(!stateEvaluationActive_ || !matchesStateEvaluation(state))
        for (std::size_t i=0;i<preparedPhase_.size();++i)
            if (!std::isfinite(omega[i]*elapsed))
                return {WVKernelStatusCode::numericalFailure,"Constant-stratification phase overflow."};
    coefficientExecutor_->execute(preparedPhase_.size(),[&](std::size_t begin,std::size_t end) {
        for (std::size_t i=begin;i<end;++i) preparedPhase_[i]=phase(omega[i]*elapsed);
    });
    ++metrics_.phasePreparationCount;
    if (stateEvaluationActive_) preparedPhaseReady_=true;
    return WVKernelStatus::ok();
}

WVComplexConstView WVTransformConstantStratificationKernel::phaseForPreparedState() const noexcept {
    return {preparedPhase_.data(),descriptor_.spectralShape()};
}

WVKernelStatus WVTransformConstantStratificationKernel::beginStateEvaluation(const WVState& state) {
    return beginStateEvaluation(state,nullptr);
}

WVKernelStatus WVTransformConstantStratificationKernel::beginStateEvaluation(
    const WVState& state,const void* evaluationOwner) {
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (stateEvaluationActive_)
        return {WVKernelStatusCode::reentrantExecution,"Constant-stratification state evaluation is already active."};
    auto status=validateState(state,true); if (!status) return status;
    preparedState_=state;
    preparedStateViews_[0]=state;
    preparedStateComponents_[0]=0;
    preparedStateViewCount_=1;
    preparedStateOwner_=evaluationOwner;
    preparedPhaseReady_=false;
    stateEvaluationActive_=true;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::addStateEvaluationView(
    const WVState& state,const void* evaluationOwner,std::size_t componentIdentity) {
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No constant-stratification state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Constant-stratification state view owner does not match the active evaluation."};
    if (state.t!=preparedState_.t || state.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Additional constant-stratification state view must use the active evaluation times."};
    if (componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"Constant-stratification component identity is out of range."};
    if (matchesStateEvaluation(state))
        return stateEvaluationComponent(state)==componentIdentity ? WVKernelStatus::ok() :
            WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                "Constant-stratification state view is already registered with another component identity."};
    if (preparedStateViewCount_==preparedStateViews_.size())
        return {WVKernelStatusCode::invalidConfiguration,"Constant-stratification state evaluation view capacity exceeded."};
    auto status=validateState(state,true); if (!status) return status;
    preparedStateViews_[preparedStateViewCount_]=state;
    preparedStateComponents_[preparedStateViewCount_++]=componentIdentity;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::removeStateEvaluationView(
    const WVState& state,const void* evaluationOwner,std::size_t componentIdentity) {
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No constant-stratification state evaluation is active."};
    if (evaluationOwner==nullptr || evaluationOwner!=preparedStateOwner_)
        return {WVKernelStatusCode::invalidConfiguration,"Constant-stratification state view owner does not match the active evaluation."};
    if (componentIdentity==0 || componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"The primary constant-stratification state view cannot be removed."};
    if (state.t!=preparedState_.t || state.t0!=preparedState_.t0)
        return {WVKernelStatusCode::invalidConfiguration,"Removed constant-stratification state view must use the active evaluation times."};
    const auto sameView=[](WVComplexConstView a,WVComplexConstView b) {
        return a.data==b.data && a.shape.rows==b.shape.rows && a.shape.columns==b.shape.columns;
    };
    for(std::size_t i=1;i<preparedStateViewCount_;++i) if (
        preparedStateComponents_[i]==componentIdentity &&
        sameView(state.coefficients.Ap,preparedStateViews_[i].coefficients.Ap) &&
        sameView(state.coefficients.Am,preparedStateViews_[i].coefficients.Am) &&
        sameView(state.coefficients.A0,preparedStateViews_[i].coefficients.A0)) {
        for(std::size_t j=i+1;j<preparedStateViewCount_;++j) {
            preparedStateViews_[j-1]=preparedStateViews_[j];
            preparedStateComponents_[j-1]=preparedStateComponents_[j];
        }
        --preparedStateViewCount_;
        preparedStateViews_[preparedStateViewCount_]={};
        preparedStateComponents_[preparedStateViewCount_]=0;
        return WVKernelStatus::ok();
    }
    return {WVKernelStatusCode::invalidConfiguration,
        "Constant-stratification state view is not registered for this component."};
}

WVKernelStatus WVTransformConstantStratificationKernel::endStateEvaluation() {
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!stateEvaluationActive_)
        return {WVKernelStatusCode::invalidConfiguration,"No constant-stratification state evaluation is active."};
    preparedState_={};
    for (auto& stateView:preparedStateViews_) stateView={};
    preparedStateComponents_={};
    preparedStateViewCount_=0;
    preparedStateOwner_=nullptr;
    stateEvaluationActive_=false;
    preparedPhaseReady_=false;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::preparedPhase(
    const WVState& state,WVComplexConstView& result) {
    result={};
    auto status=validateStateForCall(state); if (!status) return status;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!stateEvaluationActive_ || !preparedPhaseReady_) {
        status=prepareStatePhase(state); if (!status) return status;
    }
    result=phaseForPreparedState();
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::applyVertical(
    WVStratifiedModalOperator operation,WVComplexConstView input,WVComplexView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    using Op=WVStratifiedModalOperator;
    const bool projection=operation==Op::projectF || operation==Op::projectG ||
        operation==Op::projectFw || operation==Op::projectGw;
    const bool reconstruction=operation==Op::reconstructF || operation==Op::reconstructG ||
        operation==Op::reconstructFw || operation==Op::reconstructGw;
    const bool balanced=operation==Op::balancedGToWaveG;
    if (!projection && !reconstruction && !balanced)
        return {WVKernelStatusCode::unsupportedOperation,"Unsupported constant vertical operation."};
    const auto& c=descriptor_.configuration();
    const WVShape2D inputShape{projection ? c.Nz : c.Nj,descriptor_.Nkl()};
    const WVShape2D outputShape{reconstruction ? c.Nz : c.Nj,descriptor_.Nkl()};
    auto status=validateComplexShape(input,inputShape,"Constant vertical input"); if (!status) return status;
    status=validateComplexShape(output,outputShape,"Constant vertical output"); if (!status) return status;
    const auto inputBytes=inputShape.elementCount()*sizeof(WVComplex64);
    const auto outputBytes=outputShape.elementCount()*sizeof(WVComplex64);
    status=validateMutableOutputOutsidePreparedState(output.data,outputBytes); if (!status) return status;
    if (memoryOverlaps(input.data,inputBytes,output.data,outputBytes))
        return {WVKernelStatusCode::overlappingArrays,"Constant vertical input and output overlap."};
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& modes=descriptor_.verticalModes();
    if (balanced) {
        for (std::size_t column=0;column<descriptor_.Nkl();++column)
            for (std::size_t j=0;j<c.Nj;++j)
                output.data[j+c.Nj*column]=multiply(input.data[j+c.Nj*column],
                    modes.Gg[j]/modes.gWaveScale[j]);
        return WVKernelStatus::ok();
    }
    const bool sine=operation==Op::projectG || operation==Op::reconstructG ||
        operation==Op::projectGw || operation==Op::reconstructGw;
    const bool wave=operation==Op::projectFw || operation==Op::reconstructFw ||
        operation==Op::projectGw || operation==Op::reconstructGw;
    auto* working=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    const auto workingCount=c.Nz*descriptor_.Nkl();
    std::fill(working,working+workingCount,WVComplex64{});
    if (projection) std::copy(input.data,input.data+workingCount,working);
    else for (std::size_t column=0;column<descriptor_.Nkl();++column)
        for (std::size_t j=0;j<c.Nj;++j) {
            const double factor=sine ? (wave ? modes.gWaveScale[j] : modes.Gg[j]) :
                (wave ? modes.fWaveScale[j+c.Nj*column] : modes.Fg[j]);
            working[j+c.Nz*column]=multiply(input.data[j+c.Nj*column],factor);
        }
    WVFFTPlan* plan=compact_ ? plans_[sine ? verticalDST1Storage1 : verticalDCT1Storage1].get() :
        matlabPlans_[sine ? matlabRetainedDST : matlabRetainedDCT].get();
    if (projection) {
        status=plan->execute(working+(sine ? 1 : 0),working+(sine ? 1 : 0)); if (!status) return status;
        if (sine) normalizeForwardDST(working,c.Nz,descriptor_.Nkl(),1,0,1);
        else normalizeForwardDCT(working,c.Nz,descriptor_.Nkl(),1,0,1);
        for (std::size_t column=0;column<descriptor_.Nkl();++column)
            for (std::size_t j=0;j<c.Nj;++j) {
                const double factor=sine ? (wave ? modes.gWaveScale[j] : modes.Gg[j]) :
                    (wave ? modes.fWaveScale[j+c.Nj*column] : modes.Fg[j]);
                output.data[j+c.Nj*column]=multiply(working[j+c.Nz*column],1.0/factor);
            }
    } else {
        if (sine) normalizeInverseDST(working,c.Nz,descriptor_.Nkl(),1,0,1);
        else normalizeInverseDCT(working,c.Nz,descriptor_.Nkl(),1,0,1);
        status=plan->execute(working+(sine ? 1 : 0),working+(sine ? 1 : 0)); if (!status) return status;
        std::copy(working,working+workingCount,output.data);
    }
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::applyVerticalColumn(
    WVStratifiedModalOperator operation,std::size_t retainedColumn,
    WVComplexConstView input,WVComplexView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    using Op=WVStratifiedModalOperator;
    const bool projection=operation==Op::projectF || operation==Op::projectG ||
        operation==Op::projectFw || operation==Op::projectGw;
    const bool reconstruction=operation==Op::reconstructF || operation==Op::reconstructG ||
        operation==Op::reconstructFw || operation==Op::reconstructGw;
    const bool balanced=operation==Op::balancedGToWaveG;
    if (!projection && !reconstruction && !balanced)
        return {WVKernelStatusCode::unsupportedOperation,"Unsupported constant vertical-column operation."};
    if (retainedColumn>=descriptor_.Nkl())
        return {WVKernelStatusCode::invalidShape,"Constant retained column is out of range."};
    const auto& c=descriptor_.configuration();
    const WVShape2D inputShape{projection ? c.Nz : c.Nj,1};
    const WVShape2D outputShape{reconstruction ? c.Nz : c.Nj,1};
    auto status=validateComplexShape(input,inputShape,"Constant vertical-column input"); if (!status) return status;
    status=validateComplexShape(output,outputShape,"Constant vertical-column output"); if (!status) return status;
    const auto inputBytes=inputShape.elementCount()*sizeof(WVComplex64);
    const auto outputBytes=outputShape.elementCount()*sizeof(WVComplex64);
    status=validateMutableOutputOutsidePreparedState(output.data,outputBytes); if (!status) return status;
    if (memoryOverlaps(input.data,inputBytes,output.data,outputBytes))
        return {WVKernelStatusCode::overlappingArrays,"Constant vertical-column input and output overlap."};
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& modes=descriptor_.verticalModes();
    if (balanced) {
        for (std::size_t j=0;j<c.Nj;++j)
            output.data[j]=multiply(input.data[j],modes.Gg[j]/modes.gWaveScale[j]);
        return WVKernelStatus::ok();
    }
    const bool sine=operation==Op::projectG || operation==Op::reconstructG ||
        operation==Op::projectGw || operation==Op::reconstructGw;
    const bool wave=operation==Op::projectFw || operation==Op::reconstructFw ||
        operation==Op::projectGw || operation==Op::reconstructGw;
    auto* working=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(working,working+c.Nz,WVComplex64{});
    if (projection) std::copy(input.data,input.data+c.Nz,working);
    else for (std::size_t j=0;j<c.Nj;++j) {
        const double factor=sine ? (wave ? modes.gWaveScale[j] : modes.Gg[j]) :
            (wave ? modes.fWaveScale[j+c.Nj*retainedColumn] : modes.Fg[j]);
        working[j]=multiply(input.data[j],factor);
    }
    auto* plan=matlabPlans_[sine ? matlabColumnDST : matlabColumnDCT].get();
    if (projection) {
        status=plan->execute(working+(sine ? 1 : 0),working+(sine ? 1 : 0)); if (!status) return status;
        if (sine) normalizeForwardDST(working,c.Nz,1,1,0,1);
        else normalizeForwardDCT(working,c.Nz,1,1,0,1);
        for (std::size_t j=0;j<c.Nj;++j) {
            const double factor=sine ? (wave ? modes.gWaveScale[j] : modes.Gg[j]) :
                (wave ? modes.fWaveScale[j+c.Nj*retainedColumn] : modes.Fg[j]);
            output.data[j]=multiply(working[j],1.0/factor);
        }
    } else {
        if (sine) normalizeInverseDST(working,c.Nz,1,1,0,1);
        else normalizeInverseDCT(working,c.Nz,1,1,0,1);
        status=plan->execute(working+(sine ? 1 : 0),working+(sine ? 1 : 0)); if (!status) return status;
        std::copy(working,working+c.Nz,output.data);
    }
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::applyVerticalCalculus(
    WVRealConstView input,bool inputIsF,unsigned order,bool integral,WVRealView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    if ((integral && order!=1) || (!integral && (order<1 || order>4)))
        return {WVKernelStatusCode::invalidConfiguration,
            integral ? "Constant vertical integration supports only order one." :
                "Constant vertical differentiation order must be from one through four."};
    const auto& c=descriptor_.configuration();
    if (input.shape.rows!=c.Nz || output.shape.rows!=c.Nz ||
        input.shape.columns==0 || output.shape.columns!=input.shape.columns)
        return {WVKernelStatusCode::invalidShape,
            "Constant vertical calculus requires matching [Nz,Ncolumn] arrays."};
    if (!input.data || reinterpret_cast<std::uintptr_t>(input.data)%alignof(double) ||
        !output.data || reinterpret_cast<std::uintptr_t>(output.data)%alignof(double))
        return {WVKernelStatusCode::invalidPointer,
            "Constant vertical calculus requires aligned real storage."};
    std::size_t count=0,bytes=0;
    try {
        count=checkedProduct(c.Nz,input.shape.columns);
        bytes=checkedProduct(count,sizeof(double));
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow,error.what()};
    }
    auto status=validateMutableOutputOutsidePreparedState(output.data,bytes); if (!status) return status;
    if (memoryOverlaps(input.data,bytes,output.data,bytes))
        return {WVKernelStatusCode::overlappingArrays,
            "Constant vertical-calculus input and output overlap."};
    ExecutionGuard guard(executing_);
    if (!guard.entered())
        return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};

    auto* working=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    const bool outputIsF=integral ? !inputIsF : ((order%2)==0 ? inputIsF : !inputIsF);
    const bool inputIsSine=!inputIsF,outputIsSine=!outputIsF;
    auto* forward=matlabPlans_[inputIsSine ? matlabColumnDST : matlabColumnDCT].get();
    auto* inverse=matlabPlans_[outputIsSine ? matlabColumnDST : matlabColumnDCT].get();
    const double derivativeSignsF[]={-1.0,-1.0,1.0,1.0};
    const double derivativeSignsG[]={1.0,-1.0,-1.0,1.0};
    const double derivativeSign=(inputIsF ? derivativeSignsF : derivativeSignsG)[order-1];
    for (std::size_t column=0;column<input.shape.columns;++column) {
        for (std::size_t z=0;z<c.Nz;++z)
            working[z]={input.data[z+c.Nz*column],0.0};
        status=forward->execute(working+(inputIsSine?1:0),working+(inputIsSine?1:0));
        if (!status) return status;
        if (inputIsSine) normalizeForwardDST(working,c.Nz,1,1,0,1);
        else normalizeForwardDCT(working,c.Nz,1,1,0,1);
        for (std::size_t j=0;j<c.Nz;++j) {
            double multiplier=0.0;
            if (j<c.Nj && j!=0) {
                const double m=pi*static_cast<double>(j)/c.Lz;
                multiplier=integral ? (inputIsF ? 1.0 : -1.0)/m :
                    derivativeSign*std::pow(m,static_cast<int>(order));
            }
            working[j]=multiply(working[j],multiplier);
        }
        if (outputIsSine) normalizeInverseDST(working,c.Nz,1,1,0,1);
        else normalizeInverseDCT(working,c.Nz,1,1,0,1);
        status=inverse->execute(working+(outputIsSine?1:0),working+(outputIsSine?1:0));
        if (!status) return status;
        const double bottom=integral && inputIsSine ? working[0].real : 0.0;
        for (std::size_t z=0;z<c.Nz;++z)
            output.data[z+c.Nz*column]=working[z].real-bottom;
        metrics_.executionCount+=2; metrics_.verticalExecutionCount+=2;
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::horizontalForward(
    WVRealVolumeConstView input,WVComplexView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    const auto spatial=descriptor_.spatialShape(); const WVShape2D shape{descriptor_.configuration().Nz,descriptor_.Nkl()};
    auto status=validateVolume(input,spatial,"Constant horizontal input"); if (!status) return status;
    status=validateComplexShape(output,shape,"Constant horizontal output"); if (!status) return status;
    const auto inputBytes=spatial.elementCount()*sizeof(double),outputBytes=shape.elementCount()*sizeof(WVComplex64);
    status=validateMutableOutputOutsidePreparedState(output.data,outputBytes); if (!status) return status;
    if (memoryOverlaps(input.data,inputBytes,output.data,outputBytes))
        return {WVKernelStatusCode::overlappingArrays,"Constant horizontal input and output overlap."};
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    auto* half=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    status=scalarPlan(horizontalForward1)->execute(input.data,half); if (!status) return status;
    const auto& mapping=descriptor_.halfSpectrumMappings(); const auto Nz=descriptor_.configuration().Nz;
    const double scale=1.0/static_cast<double>(descriptor_.configuration().Nx*descriptor_.configuration().Ny);
    for (std::size_t mode=0;mode<descriptor_.Nkl();++mode) {
        const auto row=mapping.storageRowsByWVIndex[mode];
        const bool self=std::find(mapping.selfConjugateRows.begin(),mapping.selfConjugateRows.end(),row)!=mapping.selfConjugateRows.end();
        for (std::size_t z=0;z<Nz;++z) {
            auto value=half[z+Nz*row];
            if (mapping.conjugatesStoredValueByWVIndex[mode]) value=conjugate(value);
            if (self) value.imag=0;
            output.data[z+Nz*mode]=multiply(value,scale);
        }
    }
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::horizontalInverse(
    WVComplexConstView input,WVRealVolumeView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    const auto spatial=descriptor_.spatialShape(); const WVShape2D shape{descriptor_.configuration().Nz,descriptor_.Nkl()};
    auto status=validateComplexShape(input,shape,"Constant horizontal input"); if (!status) return status;
    status=validateVolume(output,spatial,"Constant horizontal output"); if (!status) return status;
    const auto inputBytes=shape.elementCount()*sizeof(WVComplex64),outputBytes=spatial.elementCount()*sizeof(double);
    status=validateMutableOutputOutsidePreparedState(output.data,outputBytes); if (!status) return status;
    if (memoryOverlaps(input.data,inputBytes,output.data,outputBytes))
        return {WVKernelStatusCode::overlappingArrays,"Constant horizontal input and output overlap."};
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& c=descriptor_.configuration(); const auto& mapping=descriptor_.halfSpectrumMappings();
    const auto rows=mapping.NxHalf*c.Ny; auto* half=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half+c.Nz*rows,WVComplex64{});
    for (std::size_t mode=0;mode<descriptor_.Nkl();++mode) {
        const auto row=mapping.storageRowsByWVIndex[mode];
        const bool self=std::find(mapping.selfConjugateRows.begin(),mapping.selfConjugateRows.end(),row)!=mapping.selfConjugateRows.end();
        for (std::size_t z=0;z<c.Nz;++z) {
            auto value=input.data[z+c.Nz*mode];
            if (self) value.imag=0;
            if (mapping.conjugatesStoredValueByWVIndex[mode]) value=conjugate(value);
            half[z+c.Nz*row]=value;
        }
    }
    completeHermitianBoundaries(half,mapping,c.Nz,1);
    status=scalarInversePlan_->execute(half,output.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::differentiateHorizontal(
    WVRealVolumeConstView input,bool xDerivative,WVRealVolumeView output) {
    return differentiateHorizontal(input,xDerivative,1,output);
}
WVKernelStatus WVTransformConstantStratificationKernel::differentiateHorizontal(
    WVRealVolumeConstView input,bool xDerivative,unsigned order,WVRealVolumeView output) {
    if (!matlabPrimitivesPrepared_)
        return {WVKernelStatusCode::unsupportedOperation,"MATLAB primitives were not prepared."};
    if (!order) return {WVKernelStatusCode::invalidConfiguration,"Horizontal derivative order must be positive."};
    const auto spatial=descriptor_.spatialShape(); auto status=validateVolume(input,spatial,"Constant derivative input"); if (!status) return status;
    status=validateVolume(output,spatial,"Constant derivative output"); if (!status) return status;
    const auto bytes=spatial.elementCount()*sizeof(double); status=validateMutableOutputOutsidePreparedState(output.data,bytes); if (!status) return status;
    if (memoryOverlaps(input.data,bytes,output.data,bytes)) return {WVKernelStatusCode::overlappingArrays,"Constant derivative input and output overlap."};
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& c=descriptor_.configuration();
    auto* values=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    status=scalarPlan(horizontalForward1)->execute(input.data,values); if (!status) return status;
    status=kernel_detail::applyHorizontalDerivativeMultiplier(values,c.Nx,c.Ny,c.Nz,
        1,c.Nz,c.Lx,c.Ly,xDerivative,order); if (!status) return status;
    status=scalarInversePlan_->execute(values,output.data); if (!status) return status;
    metrics_.executionCount+=2; metrics_.horizontalExecutionCount+=2; return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVEtaToWaveVortex(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients) {
    if (!descriptor_.configuration().isHydrostatic) return {WVKernelStatusCode::invalidConfiguration, "transformUVEtaToWaveVortex requires a hydrostatic kernel."};
    const auto status = validateBundle(fields, descriptor_.spatialShape(), 3, "Hydrostatic fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus outputStatuses[] = {validateSpectral(coefficients.Ap, spectral, "Ap"), validateSpectral(coefficients.Am, spectral, "Am"), validateSpectral(coefficients.A0, spectral, "A0")};
    for (const auto& value : outputStatuses) if (!value) return value;
    if (const auto ownership = validateForwardOwnership(fields,coefficients,spectral); !ownership) return ownership;
    if (const auto ownership=validateMutableOutputOutsidePreparedState(coefficients); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    WVComplexConstView phaseValues;
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0)
            return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active constant-stratification evaluation."};
        if (!preparedPhaseReady_) { auto phaseStatus=prepareStatePhase(preparedState_); if (!phaseStatus) return phaseStatus; }
        phaseValues=phaseForPreparedState();
    } else { auto phaseStatus=prepareStatePhase({t,t0,{}}); if (!phaseStatus) return phaseStatus; phaseValues=phaseForPreparedState(); }
    return transformUVEtaToWaveVortexImpl(fields,t,t0,coefficients,phaseValues);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVEtaToWaveVortex(
    WVRealVolumeConstView u,WVRealVolumeConstView v,WVRealVolumeConstView eta,
    double t,double t0,WVMutableCoefficients& coefficients) {
    if (!descriptor_.configuration().isHydrostatic)
        return {WVKernelStatusCode::invalidConfiguration,"transformUVEtaToWaveVortex requires a hydrostatic kernel."};
    const auto spatial=descriptor_.spatialShape();
    for (const auto& item:std::array<std::pair<WVRealVolumeConstView,const char*>,3>{{{u,"u"},{v,"v"},{eta,"eta"}}}) {
        auto status=validateVolume(item.first,spatial,item.second); if (!status) return status;
    }
    const auto spectral=descriptor_.spectralShape();
    for (const auto& item:std::array<std::pair<WVComplexView,const char*>,3>{{{coefficients.Ap,"Ap"},{coefficients.Am,"Am"},{coefficients.A0,"A0"}}}) {
        auto status=validateComplexShape(item.first,spectral,item.second); if (!status) return status;
    }
    const auto realBytes=spatial.elementCount()*sizeof(double),complexBytes=spectral.elementCount()*sizeof(WVComplex64);
    const void* inputs[]={u.data,v.data,eta.data}; const void* outputs[]={coefficients.Ap.data,coefficients.Am.data,coefficients.A0.data};
    for (std::size_t i=0;i<3;++i) {
        for (const auto* input:inputs) if (memoryOverlaps(outputs[i],complexBytes,input,realBytes))
            return {WVKernelStatusCode::overlappingArrays,"Projection outputs overlap physical inputs."};
        for (std::size_t j=i+1;j<3;++j) if (memoryOverlaps(outputs[i],complexBytes,outputs[j],complexBytes))
            return {WVKernelStatusCode::overlappingArrays,"Projection outputs overlap each other."};
    }
    auto status=validateMutableOutputOutsidePreparedState(coefficients); if (!status) return status;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto R=spatial.elementCount(); const WVRealVolumeConstView views[]={u,v,eta};
    for (std::size_t channel=0;channel<3;++channel) std::copy(views[channel].data,views[channel].data+R,realScratch_.data()+channel*R);
    metrics_.bytesCopied+=3*realBytes;
    WVComplexConstView phaseValues;
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0) return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active constant-stratification evaluation."};
        if (!preparedPhaseReady_) { status=prepareStatePhase(preparedState_); if (!status) return status; }
        phaseValues=phaseForPreparedState();
    } else { status=prepareStatePhase({t,t0,{}}); if (!status) return status; phaseValues=phaseForPreparedState(); }
    const WVRealFieldBundleConstView packed{realScratch_.data(),{spatial.first,spatial.second,spatial.third,3}};
    return transformUVEtaToWaveVortexImpl(packed,t,t0,coefficients,phaseValues);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVEtaToWaveVortexImpl(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = executionMapping();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward3]->execute(fields.data, half);
    if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    execute = plans_[verticalDCT2Storage3]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST1Storage3]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeForwardDCT(half,c.Nz,halfRows,3,0,2);
    normalizeForwardDST(half,c.Nz,halfRows,3,2,1);
    const double horizontalScale = compact_ ? 1.0 : 1.0 / static_cast<double>(c.Nx * c.Ny);
    auto* WV_KERNEL_RESTRICT outputAp=coefficients.Ap.data; auto* WV_KERNEL_RESTRICT outputAm=coefficients.Am.data; auto* outputA0=coefficients.A0.data;

    const auto coefficientProjectionStart = std::chrono::steady_clock::now();
    auto projectModes = [&](auto phaseProvidedConstant) {
        constexpr bool phaseProvided = decltype(phaseProvidedConstant)::value;
        coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t iMode = begin; iMode < end; ++iMode) {
                const auto& horizontal = descriptor_.fourierModes()[iMode];
                if (mapping.conjugatesStoredValueByWVIndex[iMode]) {
                    if (iMode == 0) projectHydrostaticMode<true,phaseProvided,true>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                    else projectHydrostaticMode<true,phaseProvided,false>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                } else {
                    if (iMode == 0) projectHydrostaticMode<false,phaseProvided,true>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                    else projectHydrostaticMode<false,phaseProvided,false>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                }
            }
        });
    };
    if (phaseValues.data == nullptr) projectModes(std::false_type{});
    else projectModes(std::true_type{});
    if (stageInstrumentationEnabled_) metrics_.coefficientProjectionSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientProjectionStart).count();
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::addLaplacianDamping(
    const WVState& state, double nu, double kappa, WVLaplacianDirection direction, WVFlux& flux) {
    if (!std::isfinite(nu) || !std::isfinite(kappa) || nu < 0.0 || kappa < 0.0 ||
        (direction != WVLaplacianDirection::horizontal && direction != WVLaplacianDirection::vertical))
        return {WVKernelStatusCode::invalidConfiguration,"Laplacian damping requires finite nonnegative coefficients and a valid direction."};
    auto status = validateStateAndFluxForCall(state,flux);
    if (!status) return status;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!stateEvaluationActive_ || !preparedPhaseReady_) { status=prepareStatePhase(state); if (!status) return status; }
    const auto& c = descriptor_.configuration();
    const auto& modes = descriptor_.verticalModes();
    for (std::size_t kl = 0; kl < descriptor_.Nkl(); ++kl) {
        const auto& horizontal = descriptor_.fourierModes()[kl];
        for (std::size_t j = 0; j < c.Nj; ++j) {
            const auto index = j+c.Nj*kl;
            const double wavenumber = direction == WVLaplacianDirection::horizontal ? horizontal.Kh : modes.verticalWavenumber[j];
            const double velocityRate = -nu*wavenumber*wavenumber;
            const double displacementRate = -kappa*wavenumber*wavenumber;
            const auto phaseValue = phaseForPreparedState().data[index];
            const auto evolved = evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValue);
            const auto U = multiply(coefficientValueForField<0>(modes,index,evolved),velocityRate);
            const auto V = multiply(coefficientValueForField<1>(modes,index,evolved),velocityRate);
            const auto W = multiply(coefficientValueForField<2>(modes,index,evolved),velocityRate);
            const auto N = multiply(coefficientValueForField<3>(modes,index,evolved),displacementRate);
            const auto A0 = geostrophicCoefficient(U,V,N,horizontal.k,horizontal.l,modes.A0FromVorticity[index],modes.A0FromBuoyancy[index]);
            detail::WaveCoefficientPair coefficients;
            if (kl == 0) coefficients = inertialCoefficientPair(U,V,modes.inertialScale[j]);
            else {
                WVComplex64 waveContribution;
                if (c.isHydrostatic) {
                    const auto divergence = add(multiply(U,WVComplex64{0.0,horizontal.k}),multiply(V,WVComplex64{0.0,horizontal.l}));
                    waveContribution = multiply(modes.ApmDProjection[index],divergence);
                } else {
                    const auto delta = multiply(add(multiply(U,horizontal.cosAlpha),multiply(V,horizontal.sinAlpha)),modes.ApmDScaled[index]);
                    const auto wBar = multiply(WVComplex64{0.0,(horizontal.Kh/2.0)*modes.apmWProjectionPrefactor[j]},W);
                    waveContribution = add(delta,wBar);
                }
                coefficients = waveCoefficientPair(waveContribution,buoyancyProjection(N,A0,modes.NA0Field[index],modes.ApmNProjection[index]));
            }
            detail::accumulateAtReferenceTime(flux.Fp.data[index],flux.Fm.data[index],coefficients,phaseValue);
            flux.F0.data[index] = add(flux.F0.data[index],A0);
        }
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVWEtaToWaveVortex(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients) {
    if (descriptor_.configuration().isHydrostatic) return {WVKernelStatusCode::invalidConfiguration, "transformUVWEtaToWaveVortex requires a nonhydrostatic kernel."};
    auto status = validateBundle(fields, descriptor_.spatialShape(), 4, "Nonhydrostatic fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus outputStatuses[] = {validateSpectral(coefficients.Ap, spectral, "Ap"), validateSpectral(coefficients.Am, spectral, "Am"), validateSpectral(coefficients.A0, spectral, "A0")};
    for (const auto& value : outputStatuses) if (!value) return value;
    if (const auto ownership = validateForwardOwnership(fields,coefficients,spectral); !ownership) return ownership;
    if (const auto ownership=validateMutableOutputOutsidePreparedState(coefficients); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    WVComplexConstView phaseValues;
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0)
            return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active constant-stratification evaluation."};
        if (!preparedPhaseReady_) { auto phaseStatus=prepareStatePhase(preparedState_); if (!phaseStatus) return phaseStatus; }
        phaseValues=phaseForPreparedState();
    } else { auto phaseStatus=prepareStatePhase({t,t0,{}}); if (!phaseStatus) return phaseStatus; phaseValues=phaseForPreparedState(); }
    return transformUVWEtaToWaveVortexImpl(fields,t,t0,coefficients,phaseValues);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVWEtaToWaveVortex(
    WVRealVolumeConstView u,WVRealVolumeConstView v,WVRealVolumeConstView w,
    WVRealVolumeConstView eta,double t,double t0,WVMutableCoefficients& coefficients) {
    if (descriptor_.configuration().isHydrostatic)
        return {WVKernelStatusCode::invalidConfiguration,"transformUVWEtaToWaveVortex requires a nonhydrostatic kernel."};
    const auto spatial=descriptor_.spatialShape();
    for (const auto& item:std::array<std::pair<WVRealVolumeConstView,const char*>,4>{{{u,"u"},{v,"v"},{w,"w"},{eta,"eta"}}}) {
        auto status=validateVolume(item.first,spatial,item.second); if (!status) return status;
    }
    const auto spectral=descriptor_.spectralShape();
    for (const auto& item:std::array<std::pair<WVComplexView,const char*>,3>{{{coefficients.Ap,"Ap"},{coefficients.Am,"Am"},{coefficients.A0,"A0"}}}) {
        auto status=validateComplexShape(item.first,spectral,item.second); if (!status) return status;
    }
    const auto realBytes=spatial.elementCount()*sizeof(double),complexBytes=spectral.elementCount()*sizeof(WVComplex64);
    const void* inputs[]={u.data,v.data,w.data,eta.data}; const void* outputs[]={coefficients.Ap.data,coefficients.Am.data,coefficients.A0.data};
    for (std::size_t i=0;i<3;++i) {
        for (const auto* input:inputs) if (memoryOverlaps(outputs[i],complexBytes,input,realBytes))
            return {WVKernelStatusCode::overlappingArrays,"Projection outputs overlap physical inputs."};
        for (std::size_t j=i+1;j<3;++j) if (memoryOverlaps(outputs[i],complexBytes,outputs[j],complexBytes))
            return {WVKernelStatusCode::overlappingArrays,"Projection outputs overlap each other."};
    }
    auto status=validateMutableOutputOutsidePreparedState(coefficients); if (!status) return status;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto R=spatial.elementCount(); const WVRealVolumeConstView views[]={u,v,w,eta};
    for (std::size_t channel=0;channel<4;++channel) std::copy(views[channel].data,views[channel].data+R,realScratch_.data()+channel*R);
    metrics_.bytesCopied+=4*realBytes;
    WVComplexConstView phaseValues;
    if (stateEvaluationActive_) {
        if (t!=preparedState_.t || t0!=preparedState_.t0) return {WVKernelStatusCode::invalidConfiguration,"Projection time does not match the active constant-stratification evaluation."};
        if (!preparedPhaseReady_) { status=prepareStatePhase(preparedState_); if (!status) return status; }
        phaseValues=phaseForPreparedState();
    } else { status=prepareStatePhase({t,t0,{}}); if (!status) return status; phaseValues=phaseForPreparedState(); }
    const WVRealFieldBundleConstView packed{realScratch_.data(),{spatial.first,spatial.second,spatial.third,4}};
    return transformUVWEtaToWaveVortexImpl(packed,t,t0,coefficients,phaseValues);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVWEtaToWaveVortexImpl(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues) {
    const auto& c = descriptor_.configuration(); const auto& mapping = executionMapping(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward4]->execute(fields.data, half); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    execute = plans_[verticalDCT2Storage4]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST2Storage4]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeForwardDCT(half,c.Nz,halfRows,4,0,2); normalizeForwardDST(half,c.Nz,halfRows,4,2,2);
    const double horizontalScale = compact_ ? 1.0 : 1.0 / static_cast<double>(c.Nx * c.Ny);
    auto* WV_KERNEL_RESTRICT outputAp=coefficients.Ap.data; auto* WV_KERNEL_RESTRICT outputAm=coefficients.Am.data; auto* outputA0=coefficients.A0.data;
    const auto coefficientProjectionStart = std::chrono::steady_clock::now();
    auto projectModes = [&](auto phaseProvidedConstant) {
        constexpr bool phaseProvided = decltype(phaseProvidedConstant)::value;
        coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t iMode = begin; iMode < end; ++iMode) {
                const auto& horizontal = descriptor_.fourierModes()[iMode];
                if (mapping.conjugatesStoredValueByWVIndex[iMode]) {
                    if (iMode == 0) projectNonhydrostaticMode<true,phaseProvided,true>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                    else projectNonhydrostaticMode<true,phaseProvided,false>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                } else {
                    if (iMode == 0) projectNonhydrostaticMode<false,phaseProvided,true>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                    else projectNonhydrostaticMode<false,phaseProvided,false>(half,mapping,modes,horizontal,phaseValues,t-t0,horizontalScale,c.Nz,c.Nj,iMode,outputAp,outputAm,outputA0);
                }
            }
        });
    };
    if (phaseValues.data == nullptr) projectModes(std::false_type{});
    else projectModes(std::true_type{});
    if (stageInstrumentationEnabled_) metrics_.coefficientProjectionSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientProjectionStart).count();
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVWEta(const WVState& state, WVRealFieldBundleView& fields) {
    auto status = validateBundle(fields, descriptor_.spatialShape(), 4, "Reconstructed fields"); if (!status) return status;
    status=validateStateForCall(state); if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    if (const auto ownership = validateInverseOwnership(state,fields,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    const bool generatePhase=!stateEvaluationActive_ || !preparedPhaseReady_;
    status=transformWaveVortexToUVWEtaImpl(state,fields,nullptr,
        generatePhase ? WVComplexConstView{} : phaseForPreparedState(),
        generatePhase ? preparedPhase_.data() : nullptr);
    if (status && generatePhase) { ++metrics_.phasePreparationCount; if (stateEvaluationActive_) preparedPhaseReady_=true; }
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformCoefficientTendencyToUVWEta(
    const WVState& state,WVRealFieldBundleView& fields) {
    auto status=validateBundle(fields,descriptor_.spatialShape(),4,"Coefficient tendency fields");
    if (!status) return status;
    ++metrics_.derivedValidationCount;
    status=validateStateContents(state); if (!status) return status;
    if (stateEvaluationActive_ &&
        (state.t!=preparedState_.t || state.t0!=preparedState_.t0))
        return {WVKernelStatusCode::invalidConfiguration,
            "Coefficient tendency times do not match the active evaluation."};
    if (const auto ownership=validateInverseOwnership(state,fields,descriptor_.spectralShape()); !ownership)
        return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered())
        return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const bool generatePhase=!stateEvaluationActive_ || !preparedPhaseReady_;
    status=transformWaveVortexToUVWEtaImpl(state,fields,nullptr,
        generatePhase ? WVComplexConstView{} : phaseForPreparedState(),
        generatePhase ? preparedPhase_.data() : nullptr,false);
    if (status) for (auto& count:metrics_.tendencyReconstructionCount) ++count;
    if (status && generatePhase) {
        ++metrics_.phasePreparationCount;
        if (stateEvaluationActive_) preparedPhaseReady_=true;
    }
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVW(const WVState& state, WVRealFieldBundleView& fields) {
    auto status = validateBundle(fields,descriptor_.spatialShape(),3,"Advection fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    status=validateStateForCall(state); if (!status) return status;
    if (const auto ownership = validateInverseOwnership(state,fields,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const bool generatePhase=!stateEvaluationActive_ || !preparedPhaseReady_;
    ++metrics_.advectionVelocityReconstructionCount;
    status=transformWaveVortexToUVWImpl(state,fields,nullptr,
        generatePhase ? WVComplexConstView{} : phaseForPreparedState(),
        generatePhase ? preparedPhase_.data() : nullptr);
    if (status && generatePhase) { ++metrics_.phasePreparationCount; if (stateEvaluationActive_) preparedPhaseReady_=true; }
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUV(
    const WVState& state, WVRealFieldBundleView& fields) {
    auto status=validateBundle(fields,descriptor_.spatialShape(),2,"Horizontal velocity fields");
    if (!status) return status;
    status=validateStateForCall(state); if (!status) return status;
    if (const auto ownership=validateInverseOwnership(state,fields,descriptor_.spectralShape()); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (!plans_[horizontalInverse2] || !plans_[verticalDCT2Storage2])
        return {WVKernelStatusCode::invalidConfiguration,"Call prepareHorizontalVelocityTransform during model setup before requesting horizontal velocity."};
    const bool generatePhase=!stateEvaluationActive_ || !preparedPhaseReady_;
    status=transformWaveVortexToUVImpl(state,fields,
        generatePhase ? WVComplexConstView{} : phaseForPreparedState(),
        generatePhase ? preparedPhase_.data() : nullptr);
    if (status && generatePhase) { ++metrics_.phasePreparationCount; if (stateEvaluationActive_) preparedPhaseReady_=true; }
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformStateFieldDerivatives(
    const WVState& state, WVDynamicalField field, WVRealFieldBundleView& derivatives) {
    const auto target = static_cast<std::size_t>(field);
    if (target > static_cast<std::size_t>(WVDynamicalField::eta))
        return {WVKernelStatusCode::invalidConfiguration,"Unknown dynamical field derivative target."};
    auto status = validateBundle(derivatives,descriptor_.spatialShape(),3,"Field derivatives");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    status=validateStateForCall(state); if (!status) return status;
    if (const auto ownership = validateInverseOwnership(state,derivatives,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const bool generatePhase=!stateEvaluationActive_ || !preparedPhaseReady_;
    status=transformToSpatialDomainWithDerivativesFromStateImpl(state,
        generatePhase ? WVComplexConstView{} : phaseForPreparedState(),target,derivatives,
        generatePhase ? preparedPhase_.data() : nullptr);
    if (status && generatePhase) { ++metrics_.phasePreparationCount; if (stateEvaluationActive_) preparedPhaseReady_=true; }
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVWEtaImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients, WVComplexConstView phaseValues, WVComplex64* generatedPhase, bool countPrimaryReconstruction) {
    const auto& c = descriptor_.configuration(); const auto& mapping = executionMapping(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount(); auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t iMode = begin; iMode < end; ++iMode) {
        if (evolvedCoefficients == nullptr) {
            auto source = [&](std::size_t index) {
                const auto p = generatedPhase ?
                    (generatedPhase[index]=phase(modes.omega[index]*(state.t-state.t0))) : phaseValues.data[index];
                return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],p);
            };
            if (mapping.conjugatesStoredValueByWVIndex[iMode]) assembleFieldSpectraForMode<true>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
            else assembleFieldSpectraForMode<false>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
        } else {
            auto source = [&](std::size_t index) {
                return EvolvedWaveVortexCoefficients{evolvedCoefficients->Ap.data[index],evolvedCoefficients->Am.data[index],evolvedCoefficients->A0.data[index]};
            };
            if (mapping.conjugatesStoredValueByWVIndex[iMode]) assembleFieldSpectraForMode<true>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
            else assembleFieldSpectraForMode<false>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
        }
        }
    });
    if (stageInstrumentationEnabled_) metrics_.coefficientAssemblySeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,4);
    auto execute = plans_[verticalDCT2Storage4]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST2Storage4]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    execute = plans_[horizontalInverse4]->execute(half, fields.data); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    if (countPrimaryReconstruction) {
        const auto component=stateEvaluationComponent(state);
        for (std::size_t field=0;field<4;++field)
            ++metrics_.reconstructionCount[field][0][component];
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVWImpl(
    const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients, WVComplexConstView phaseValues, WVComplex64* generatedPhase) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = executionMapping();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t iMode = begin; iMode < end; ++iMode) {
            auto source = [&](std::size_t index) {
                if (evolvedCoefficients != nullptr) return EvolvedWaveVortexCoefficients{evolvedCoefficients->Ap.data[index],evolvedCoefficients->Am.data[index],evolvedCoefficients->A0.data[index]};
                const auto phaseValue = generatedPhase ?
                    (generatedPhase[index]=phase(modes.omega[index]*(state.t-state.t0))) : phaseValues.data[index];
                return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValue);
            };
            if (mapping.conjugatesStoredValueByWVIndex[iMode]) assembleVelocitySpectraForMode<true>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
            else assembleVelocitySpectraForMode<false>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
        }
    });
    if (stageInstrumentationEnabled_) metrics_.coefficientAssemblySeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,3);
    auto execute = plans_[verticalDCT2Storage3]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST1Storage3]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    execute = plans_[horizontalInverse3]->execute(half,fields.data); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    const auto component=stateEvaluationComponent(state);
    for (std::size_t field=0;field<3;++field)
        ++metrics_.reconstructionCount[field][0][component];
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVImpl(
    const WVState& state, WVRealFieldBundleView& fields, WVComplexConstView phaseValues,
    WVComplex64* generatedPhase) {
    const auto& c=descriptor_.configuration();
    const auto& mapping=executionMapping();
    const auto& modes=descriptor_.verticalModes();
    const std::size_t halfRows=executionRowCount();
    auto* half=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half+c.Nz*2*halfRows,WVComplex64{});
    const auto coefficientAssemblyStart=std::chrono::steady_clock::now();
    coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin,std::size_t end) {
        for (std::size_t iMode=begin;iMode<end;++iMode) {
            auto source=[&](std::size_t index) {
                const auto phaseValue=generatedPhase ?
                    (generatedPhase[index]=phase(modes.omega[index]*(state.t-state.t0))) : phaseValues.data[index];
                return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValue);
            };
            if (mapping.conjugatesStoredValueByWVIndex[iMode])
                assembleHorizontalVelocitySpectraForMode<true>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
            else
                assembleHorizontalVelocitySpectraForMode<false>(half,mapping,modes,c.Nz,c.Nj,iMode,source);
        }
    });
    if (stageInstrumentationEnabled_)
        metrics_.coefficientAssemblySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,2);
    auto execute=plans_[verticalDCT2Storage2]->execute(half,half); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    execute=plans_[horizontalInverse2]->execute(half,fields.data); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    const auto component=stateEvaluationComponent(state);
    ++metrics_.reconstructionCount[0][0][component];
    ++metrics_.reconstructionCount[1][0][component];
    return WVKernelStatus::ok();
}

namespace {
WVKernelStatus validateAllDerivativeInputs(const WVTransformConstantStratificationKernel& kernel, const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    const auto spectral = kernel.descriptor().spectralShape();
    auto status = validateSpectral(Apm, spectral, "Apm"); if (!status) return status;
    status = validateSpectral(A0, spectral, "A0"); if (!status) return status;
    status = validateBundle(fields, kernel.descriptor().spatialShape(), 4, "Field derivatives"); if (!status) return status;
    status = validateDerivativeOwnership(Apm,A0,fields,spectral); if (!status) return status;
    return WVKernelStatus::ok();
}
} // namespace

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithFAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    auto initial = validateAllDerivativeInputs(*this, Apm, A0, fields); if (!initial) return initial;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration(); const auto& mapping = executionMapping(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount(); auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t mode = begin; mode < end; ++mode) {
            const auto& horizontal = descriptor_.fourierModes()[mode];
            if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleFieldFamilyDerivativesForMode<true,true>(half,mapping,modes,horizontal,Apm,A0,c.Nz,c.Nj,mode);
            else assembleFieldFamilyDerivativesForMode<true,false>(half,mapping,modes,horizontal,Apm,A0,c.Nz,c.Nj,mode);
        }
    });
    if (stageInstrumentationEnabled_) metrics_.derivativeCoefficientAssemblySeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,4);
    auto execute = plans_[verticalDCT3Storage4]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST1Storage4]->execute(half + 3 * c.Nz + 1,half + 3 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeInverseDCT(half,c.Nz,halfRows,4,0,3); normalizeInverseDST(half,c.Nz,halfRows,4,3,1);
    execute = plans_[horizontalInverse4]->execute(half, fields.data); if (!execute) return execute; ++metrics_.executionCount; ++metrics_.horizontalExecutionCount; return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithFAllDerivatives(
    WVConstantFField field,const WVComplexConstView& Apm,const WVComplexConstView& A0,
    WVRealFieldBundleView& fields,std::size_t componentIdentity) {
    if (componentIdentity>=5)
        return {WVKernelStatusCode::invalidConfiguration,"Constant-stratification component identity is out of range."};
    std::size_t target=0;
    switch(field) {
        case WVConstantFField::pi: target=4; break;
        case WVConstantFField::psi: target=6; break;
        case WVConstantFField::qgpv: target=7; break;
        default: return {WVKernelStatusCode::unsupportedOperation,"Unknown constant-stratification F field."};
    }
    const auto status=transformToSpatialDomainWithFAllDerivatives(Apm,A0,fields);
    if (status) for (std::size_t derivative=0;derivative<4;++derivative)
        ++metrics_.reconstructionCount[target][derivative][componentIdentity];
    return status;
}

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithGAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    auto initial = validateAllDerivativeInputs(*this, Apm, A0, fields); if (!initial) return initial;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration(); const auto& mapping = executionMapping(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount(); auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t mode = begin; mode < end; ++mode) {
            const auto& horizontal = descriptor_.fourierModes()[mode];
            if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleFieldFamilyDerivativesForMode<false,true>(half,mapping,modes,horizontal,Apm,A0,c.Nz,c.Nj,mode);
            else assembleFieldFamilyDerivativesForMode<false,false>(half,mapping,modes,horizontal,Apm,A0,c.Nz,c.Nj,mode);
        }
    });
    if (stageInstrumentationEnabled_) metrics_.derivativeCoefficientAssemblySeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,4);
    auto execute = plans_[verticalDST3Storage4]->execute(half + 1,half + 1); if (!execute) return execute;
    execute = plans_[verticalDCT1Storage4]->execute(half + 3 * c.Nz,half + 3 * c.Nz); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeInverseDST(half,c.Nz,halfRows,4,0,3); normalizeInverseDCT(half,c.Nz,halfRows,4,3,1);
    execute = plans_[horizontalInverse4]->execute(half, fields.data); if (!execute) return execute; ++metrics_.executionCount; ++metrics_.horizontalExecutionCount; return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithDerivativesImpl(
    const WVCoefficients& evolvedCoefficients, std::size_t target, WVRealFieldBundleView& derivatives) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = executionMapping();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const bool cosine = target < 2;
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    auto assembleTarget = [&](auto targetConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t mode = begin; mode < end; ++mode) {
                auto source = [&](std::size_t index) { return EvolvedWaveVortexCoefficients{evolvedCoefficients.Ap.data[index],evolvedCoefficients.Am.data[index],evolvedCoefficients.A0.data[index]}; };
                if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleDerivativeSpectraForMode<resolvedTarget,true>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,nullptr,source);
                else assembleDerivativeSpectraForMode<resolvedTarget,false>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,nullptr,source);
            }
        });
    };
    if (target == 0) assembleTarget(std::integral_constant<std::size_t,0>{});
    else if (target == 1) assembleTarget(std::integral_constant<std::size_t,1>{});
    else if (target == 2) assembleTarget(std::integral_constant<std::size_t,2>{});
    else assembleTarget(std::integral_constant<std::size_t,3>{});
    if (stageInstrumentationEnabled_) metrics_.derivativeCoefficientAssemblySeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,3);
    WVKernelStatus execute;
    if (cosine) {
        execute = plans_[verticalDCT2Storage3]->execute(half,half); if (!execute) return execute;
        execute = plans_[verticalDST1Storage3]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    } else {
        execute = plans_[verticalDST2Storage3]->execute(half + 1,half + 1); if (!execute) return execute;
        execute = plans_[verticalDCT1Storage3]->execute(half + 2 * c.Nz,half + 2 * c.Nz); if (!execute) return execute;
    }
    metrics_.executionCount += 2;
    metrics_.verticalExecutionCount += 2;
    execute = plans_[horizontalInverse3]->execute(half,derivatives.data); if (!execute) return execute;
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithDerivativesFromStateImpl(
    const WVState& state, WVComplexConstView phaseValues, std::size_t target,
    WVRealFieldBundleView& derivatives, WVComplex64* generatedPhase,
    const std::array<bool,3>* derivativeMask) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = executionMapping();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const bool cosine = target < 2;
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    auto assembleTarget = [&](auto targetConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t mode = begin; mode < end; ++mode) {
                auto source = [&](std::size_t index) { const auto phaseValue = generatedPhase ?
                    (generatedPhase[index]=phase(modes.omega[index]*(state.t-state.t0))) : phaseValues.data[index];
                return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValue); };
                if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleDerivativeSpectraForMode<resolvedTarget,true>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,derivativeMask,source);
                else assembleDerivativeSpectraForMode<resolvedTarget,false>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,derivativeMask,source);
            }
        });
    };
    if (target == 0) assembleTarget(std::integral_constant<std::size_t,0>{});
    else if (target == 1) assembleTarget(std::integral_constant<std::size_t,1>{});
    else if (target == 2) assembleTarget(std::integral_constant<std::size_t,2>{});
    else assembleTarget(std::integral_constant<std::size_t,3>{});
    if (stageInstrumentationEnabled_) metrics_.derivativeCoefficientAssemblySeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - coefficientAssemblyStart).count();
    completeHermitianBoundaries(half,mapping,c.Nz,3);
    WVKernelStatus execute;
    if (cosine) {
        execute = plans_[verticalDCT2Storage3]->execute(half,half); if (!execute) return execute;
        execute = plans_[verticalDST1Storage3]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    } else {
        execute = plans_[verticalDST2Storage3]->execute(half + 1,half + 1); if (!execute) return execute;
        execute = plans_[verticalDCT1Storage3]->execute(half + 2 * c.Nz,half + 2 * c.Nz); if (!execute) return execute;
    }
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    execute = plans_[horizontalInverse3]->execute(half,derivatives.data); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    const auto component=stateEvaluationComponent(state);
    for (std::size_t derivative=1;derivative<4;++derivative)
        if (!derivativeMask || (*derivativeMask)[derivative-1])
            ++metrics_.reconstructionCount[target][derivative][component];
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::projectSingleFluxTargetImpl(
    const WVRealFieldBundleConstView& field, std::size_t target, WVComplexConstView phaseValues, WVFlux& flux) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = executionMapping();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = executionRowCount();
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward1]->execute(field.data,half); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    if (target < 2) {
        execute = plans_[verticalDCT1Storage1]->execute(half,half); if (!execute) return execute;
        if (!compact_) normalizeForwardDCT(half,c.Nz,halfRows,1,0,1);
    } else {
        execute = plans_[verticalDST1Storage1]->execute(half + 1,half + 1); if (!execute) return execute;
        if (!compact_) normalizeForwardDST(half,c.Nz,halfRows,1,0,1);
        else for (std::size_t row=0;row<halfRows;++row) half[c.Nz*row]={};
    }
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    // Candidate horizontal values already use forward-unit normalization.
    // Apply the Type-I scale only to retained Nj values in the projection.
    // Nj <= Nz-1 excludes the DCT terminal endpoint; the G zero mode is zeroed
    // above because its untransformed boundary sample is not a sine coefficient.
    const double horizontalScale = compact_ ? 1.0/static_cast<double>(c.Nz-1) :
        1.0/static_cast<double>(c.Nx*c.Ny);
    auto projectTarget = [&](auto targetConstant, auto hydrostaticConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        constexpr bool hydrostatic = decltype(hydrostaticConstant)::value;
        coefficientExecutor_->execute(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t iMode = begin; iMode < end; ++iMode) {
                const auto& horizontal = descriptor_.fourierModes()[iMode];
                if (mapping.conjugatesStoredValueByWVIndex[iMode]) {
                    if (iMode == 0) projectFluxTargetMode<resolvedTarget,hydrostatic,true,true>(half,mapping,modes,horizontal,phaseValues,horizontalScale,c.Nz,c.Nj,iMode,flux);
                    else projectFluxTargetMode<resolvedTarget,hydrostatic,true,false>(half,mapping,modes,horizontal,phaseValues,horizontalScale,c.Nz,c.Nj,iMode,flux);
                } else {
                    if (iMode == 0) projectFluxTargetMode<resolvedTarget,hydrostatic,false,true>(half,mapping,modes,horizontal,phaseValues,horizontalScale,c.Nz,c.Nj,iMode,flux);
                    else projectFluxTargetMode<resolvedTarget,hydrostatic,false,false>(half,mapping,modes,horizontal,phaseValues,horizontalScale,c.Nz,c.Nj,iMode,flux);
                }
            }
        });
    };
    auto dispatchTarget = [&](auto hydrostaticConstant) {
        if (target == 0) projectTarget(std::integral_constant<std::size_t,0>{},hydrostaticConstant);
        else if (target == 1) projectTarget(std::integral_constant<std::size_t,1>{},hydrostaticConstant);
        else if (target == 2) projectTarget(std::integral_constant<std::size_t,2>{},hydrostaticConstant);
        else projectTarget(std::integral_constant<std::size_t,3>{},hydrostaticConstant);
    };
    if (c.isHydrostatic) dispatchTarget(std::true_type{});
    else dispatchTarget(std::false_type{});
    return WVKernelStatus::ok();
}

WVFFTPlan* WVTransformConstantStratificationKernel::scalarPlan(std::size_t index) const noexcept {
    return compact_ ? compact_->scalarPlans[index].get() : plans_[index].get();
}
WVKernelStatus WVTransformConstantStratificationKernel::ensureScalarPlans() {
    if (!compact_ || !compact_->scalarPlans.empty()) return WVKernelStatus::ok();
    try {
        const auto& c=descriptor_.configuration();
        const auto rows=descriptor_.halfSpectrumMappings().NxHalf*c.Ny;
        std::vector<std::unique_ptr<WVFFTPlan>> prepared(planCount);
        const std::array<std::pair<std::size_t,WVFFTPlanSpecification>,3> specs{{
            {horizontalForward1,horizontalSpecification(c,1,true)},
            {verticalDCT1Storage1,verticalSpecification(c,rows,1,1,false)},
            {verticalDST1Storage1,verticalSpecification(c,rows,1,1,true)}}};
        for (const auto& item : specs) {
            auto status=engine_->createPlan(item.second,prepared[item.first]);
            if (!status) return status;
            if (!prepared[item.first]) return {WVKernelStatusCode::fftPlanFailure,"Empty full-grid scalar plan."};
        }
        // Scalar calculus retains arbitrary full-grid horizontal modes. It
        // shares the compact arena at max(4*C,H). The immutable scalar Fourier
        // source borrows unused real scratch during scalar-only operations.
        const auto doubles=checkedProduct(checkedProduct(rows,c.Nz),2);
        if (halfSpectrumScratch_.size()<doubles) {
            std::vector<double> expanded(doubles);
            halfSpectrumScratch_.swap(expanded);
        }
        compact_->scalarPlans=std::move(prepared);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Unable to prepare full-grid scalar storage."}; }
      catch (const std::overflow_error& error) { return {WVKernelStatusCode::sizeOverflow,error.what()}; }
}

WVKernelStatus WVTransformConstantStratificationKernel::ensureScalarInversePlan() {
    if (scalarInversePlan_) return WVKernelStatus::ok();
    std::unique_ptr<WVFFTPlan> plan;
    auto status = engine_->createPlan(horizontalSpecification(descriptor_.configuration(),1,false),plan);
    if (!status) return status;
    if (!plan) return {WVKernelStatusCode::fftPlanFailure,"FFT engine returned an empty scalar inverse plan."};
    metrics_.planBytes += plan->persistentBytes();
    ++metrics_.planCount;
    scalarInversePlan_ = std::move(plan);
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::prepareScalarAdvection() {
    auto status=ensureScalarPlans();
    if (!status) return status;
    return ensureScalarInversePlan();
}

WVKernelStatus WVTransformConstantStratificationKernel::antialiasScalarInPlace(WVRealVolumeView& scalar) {
    auto status = ensureScalarInversePlan();
    if (!status) return status;
    const auto& c = descriptor_.configuration();
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    status = scalarPlan(horizontalForward1)->execute(scalar.data,half);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    const double scale = 1.0/static_cast<double>(c.Nx*c.Ny);
    for (std::size_t row = 0; row < halfRows; ++row) {
        const double rowScale = scalarAntialiasRows_[row] == 0 ? 0.0 : scale;
        for (std::size_t z = 0; z < c.Nz; ++z)
            half[z+c.Nz*row] = multiply(half[z+c.Nz*row],rowScale);
    }
    status = scalarInversePlan_->execute(half,scalar.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    ++metrics_.scalarAntialiasCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::prepareCompactScalarDerivatives(
    const WVRealVolumeConstView& scalar, bool sine) {
    auto status=ensureScalarPlans(); if (!status) return status;
    status=ensureScalarInversePlan(); if (!status) return status;
    const auto& c=descriptor_.configuration();
    const auto R=descriptor_.spatialShape().elementCount();
    const auto NxHalf=descriptor_.halfSpectrumMappings().NxHalf;
    const auto rows=NxHalf*c.Ny;
    // All public scalar inputs/velocities are caller-owned. The first3R is
    // otherwise unused by scalar calculus; 2H doubles fit there for Nx>=2.
    // The final3R retains complete derivatives until the public output commits.
    auto* source=reinterpret_cast<WVComplex64*>(realScratch_.data());
    auto* working=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto start=std::chrono::steady_clock::now();
    status=scalarPlan(horizontalForward1)->execute(scalar.data,source);
    if (!status) return status;
    metrics_.scalarForwardSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    const double scale=1.0/static_cast<double>(c.Nx*c.Ny);
    for (std::size_t target=0;target<3;++target) {
        start=std::chrono::steady_clock::now();
        for (std::size_t row=0;row<rows;++row) {
            const auto kIndex=row%NxHalf,lIndex=row/NxHalf;
            const auto lMode=lIndex<(c.Ny+1)/2 ? static_cast<std::int64_t>(lIndex) :
                static_cast<std::int64_t>(lIndex)-static_cast<std::int64_t>(c.Ny);
            const double k=c.Nx%2==0 && kIndex==c.Nx/2 ? 0.0 : 2*pi*static_cast<double>(kIndex)/c.Lx;
            // Preserve the existing G and F full-grid Nyquist conventions.
            const double l=sine && c.Ny%2==0 && lIndex==c.Ny/2 ? 0.0 : 2*pi*static_cast<double>(lMode)/c.Ly;
            for (std::size_t z=0;z<c.Nz;++z) {
                const auto value=multiply(source[z+c.Nz*row],scale);
                working[z+c.Nz*row]=target==2 ? value : multiply(value,WVComplex64{0,target==0?k:l});
            }
        }
        metrics_.scalarDerivativeAssemblySeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if (target==2) {
            start=std::chrono::steady_clock::now();
            status=scalarPlan(sine?verticalDST1Storage1:verticalDCT1Storage1)->execute(working+(sine?1:0),working+(sine?1:0));
            if (!status) return status;
            ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
            if (sine) normalizeForwardDST(working,c.Nz,rows,1,0,1);
            else normalizeForwardDCT(working,c.Nz,rows,1,0,1);
            for (std::size_t row=0;row<rows;++row) {
                const auto base=c.Nz*row;
                working[base]={};
                for (std::size_t j=1;j<c.Nj;++j)
                    working[base+j]=multiply(working[base+j],(sine?1.0:-1.0)*pi*static_cast<double>(j)/c.Lz);
                for (std::size_t j=c.Nj;j<c.Nz;++j) working[base+j]={};
            }
            if (sine) normalizeInverseDCT(working,c.Nz,rows,1,0,1);
            else normalizeInverseDST(working,c.Nz,rows,1,0,1);
            status=scalarPlan(sine?verticalDCT1Storage1:verticalDST1Storage1)->execute(working+(sine?0:1),working+(sine?0:1));
            if (!status) return status;
            ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
            metrics_.scalarVerticalDerivativeSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        }
        for (const auto row:descriptor_.halfSpectrumMappings().selfConjugateRows)
            for (std::size_t z=0;z<c.Nz;++z) working[z+c.Nz*row].imag=0;
        start=std::chrono::steady_clock::now();
        status=scalarInversePlan_->execute(working,realScratch_.data()+(3+target)*R);
        if (!status) return status;
        metrics_.scalarInverseSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformGGridScalarDerivatives(
    const WVRealVolumeConstView& scalar, WVRealFieldBundleView& derivatives) {
    const auto spatial=descriptor_.spatialShape();
    if(scalar.shape.first!=spatial.first || scalar.shape.second!=spatial.second || scalar.shape.third!=spatial.third)
        return {WVKernelStatusCode::invalidShape,"G-grid scalar must have shape [Nx,Ny,Nz]."};
    if(!scalar.data) return {WVKernelStatusCode::invalidPointer,"G-grid scalar storage is missing."};
    auto status=validateBundle(derivatives,spatial,3,"G-grid derivatives");
    if(!status) return status;
    const auto R=spatial.elementCount();
    if(memoryOverlaps(scalar.data,R*sizeof(double),derivatives.data,3*R*sizeof(double)))
        return {WVKernelStatusCode::overlappingArrays,"G-grid scalar and derivative arrays must not overlap."};
    for(std::size_t i=0;i<R;++i) if(!std::isfinite(scalar.data[i]))
        return {WVKernelStatusCode::invalidConfiguration,"G-grid scalar must be finite."};
    ExecutionGuard guard(executing_);
    if(!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    if (compact_) {
        status=prepareCompactScalarDerivatives(scalar,true); if (!status) return status;
        const auto* result=realScratch_.data()+3*R;
        for (std::size_t i=0;i<3*R;++i) if (!std::isfinite(result[i]))
            return {WVKernelStatusCode::numericalFailure,"G-grid scalar derivative overflowed."};
        std::copy_n(result,3*R,derivatives.data);
        return WVKernelStatus::ok();
    }
    status=ensureScalarPlans(); if (!status) return status;
    const auto& c=descriptor_.configuration();
    const auto NxHalf=descriptor_.halfSpectrumMappings().NxHalf;
    const auto halfRows=NxHalf*c.Ny;
    auto* half=reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    status=scalarPlan(horizontalForward1)->execute(scalar.data,half);
    if(!status) return status;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    const double scale=1.0/static_cast<double>(c.Nx*c.Ny);
    // Expand backwards because the single-channel input shares this scratch.
    for(std::size_t row=halfRows;row-- > 0;) {
        const auto kIndex=row%NxHalf,lIndex=row/NxHalf;
        const double k=c.Nx%2==0 && kIndex==c.Nx/2 ? 0.0 : 2*pi*static_cast<double>(kIndex)/c.Lx;
        const auto lMode=lIndex<(c.Ny+1)/2 ? static_cast<std::int64_t>(lIndex) : static_cast<std::int64_t>(lIndex)-static_cast<std::int64_t>(c.Ny);
        const double l=c.Ny%2==0 && lIndex==c.Ny/2 ? 0.0 : 2*pi*static_cast<double>(lMode)/c.Ly;
        for(std::size_t z=c.Nz;z-- > 0;) {
            const auto value=multiply(half[z+c.Nz*row],scale);
            half[z+c.Nz*3*row]=multiply(value,WVComplex64{0,k});
            half[z+c.Nz+3*c.Nz*row]=multiply(value,WVComplex64{0,l});
            half[z+2*c.Nz+3*c.Nz*row]=value;
        }
    }
    // G has sine basis; the first derivative is its retained cosine series.
    status=scalarPlan(verticalDST1Storage3)->execute(half+2*c.Nz+1,half+2*c.Nz+1);
    if(!status) return status;
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    normalizeForwardDST(half,c.Nz,halfRows,3,2,1);
    for(std::size_t row=0;row<halfRows;++row) {
        const auto base=2*c.Nz+3*c.Nz*row;
        for(std::size_t j=1;j<c.Nj;++j) half[base+j]=multiply(half[base+j],pi*static_cast<double>(j)/c.Lz);
        for(std::size_t j=c.Nj;j<c.Nz;++j) half[base+j]={};
    }
    normalizeInverseDCT(half,c.Nz,halfRows,3,2,1);
    status=scalarPlan(verticalDCT1Storage3)->execute(half+2*c.Nz,half+2*c.Nz);
    if(!status) return status;
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    for(const auto row:descriptor_.halfSpectrumMappings().selfConjugateRows)
        for(std::size_t z=0;z<c.Nz;++z) for(std::size_t channel=0;channel<3;++channel)
            half[z+c.Nz*channel+3*c.Nz*row].imag=0;
    auto* result=realScratch_.data()+3*R;
    status=scalarPlan(horizontalInverse3)->execute(half,result);
    if(!status) return status;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    for(std::size_t i=0;i<3*R;++i) if(!std::isfinite(result[i]))
        return {WVKernelStatusCode::numericalFailure,"G-grid scalar derivative overflowed."};
    std::copy_n(result,3*R,derivatives.data);
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::advectFGridScalar(
    const WVRealVolumeConstView& scalar,
    const WVRealFieldBundleConstView& advectionFields,
    bool shouldAntialias,
    WVRealVolumeView& rightHandSide) {
    const auto spatial = descriptor_.spatialShape();
    if (scalar.shape.first != spatial.first || scalar.shape.second != spatial.second || scalar.shape.third != spatial.third ||
        rightHandSide.shape.first != spatial.first || rightHandSide.shape.second != spatial.second || rightHandSide.shape.third != spatial.third)
        return {WVKernelStatusCode::invalidShape,"F-grid scalar input and output must have shape [Nx,Ny,Nz]."};
    if (scalar.data == nullptr || rightHandSide.data == nullptr)
        return {WVKernelStatusCode::invalidPointer,"F-grid scalar input and output pointers must be nonnull."};
    auto status = validateBundle(advectionFields,spatial,3,"Advection fields");
    if (!status) return status;
    const auto R = spatial.elementCount();
    const auto bytes = R*sizeof(double);
    if (memoryOverlaps(scalar.data,bytes,rightHandSide.data,bytes) ||
        memoryOverlaps(advectionFields.data,3*bytes,scalar.data,bytes) ||
        memoryOverlaps(advectionFields.data,3*bytes,rightHandSide.data,bytes))
        return {WVKernelStatusCode::overlappingArrays,"Tracer scalar, velocity, and RHS arrays must not overlap."};
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};

    auto stageStarted=std::chrono::steady_clock::now();
    if (compact_) {
        status=prepareCompactScalarDerivatives(scalar,false); if (!status) return status;
    } else {
    status=ensureScalarPlans(); if (!status) return status;
    const auto& c = descriptor_.configuration();
    const std::size_t NxHalf = descriptor_.halfSpectrumMappings().NxHalf;
    const std::size_t halfRows = NxHalf*c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    stageStarted = std::chrono::steady_clock::now();
    status = scalarPlan(horizontalForward1)->execute(scalar.data,half);
    if (!status) return status;
    metrics_.scalarForwardSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    const double horizontalScale = 1.0/static_cast<double>(c.Nx*c.Ny);
    stageStarted = std::chrono::steady_clock::now();
    for (std::size_t reverse = 0; reverse < halfRows; ++reverse) {
        const std::size_t row = halfRows-1-reverse;
        const std::size_t kIndex = row%NxHalf;
        const std::size_t lIndex = row/NxHalf;
        const double k = c.Nx%2 == 0 && kIndex == c.Nx/2 ? 0.0 : 2.0*pi*static_cast<double>(kIndex)/c.Lx;
        const auto lMode = lIndex < (c.Ny+1)/2 ? static_cast<std::int64_t>(lIndex) : static_cast<std::int64_t>(lIndex)-static_cast<std::int64_t>(c.Ny);
        const double l = 2.0*pi*static_cast<double>(lMode)/c.Ly;
        for (std::size_t z = c.Nz; z-- > 0;) {
            const auto value = multiply(half[z+c.Nz*row],horizontalScale);
            half[z+c.Nz*0+c.Nz*3*row] = multiply(value,WVComplex64{0.0,k});
            half[z+c.Nz*1+c.Nz*3*row] = multiply(value,WVComplex64{0.0,l});
            half[z+c.Nz*2+c.Nz*3*row] = value;
        }
    }
    metrics_.scalarDerivativeAssemblySeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    stageStarted = std::chrono::steady_clock::now();
    status = scalarPlan(verticalDCT1Storage3)->execute(half+2*c.Nz,half+2*c.Nz);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.verticalExecutionCount;
    normalizeForwardDCT(half,c.Nz,halfRows,3,2,1);
    for (std::size_t row = 0; row < halfRows; ++row) {
        const auto base = 2*c.Nz+3*c.Nz*row;
        half[base] = {};
        // MATLAB diffZF uses DCT/iDST restricted to the retained vertical
        // modes. Horizontal derivatives above still retain every grid depth.
        for (std::size_t j = 1; j < c.Nj; ++j) {
            const double verticalWavenumber = pi*static_cast<double>(j)/c.Lz;
            half[base+j] = multiply(half[base+j],-verticalWavenumber);
        }
        for (std::size_t j = c.Nj; j < c.Nz; ++j) half[base+j] = {};
    }
    normalizeInverseDST(half,c.Nz,halfRows,3,2,1);
    status = scalarPlan(verticalDST1Storage3)->execute(half+2*c.Nz+1,half+2*c.Nz+1);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.verticalExecutionCount;
    for (const auto row : descriptor_.halfSpectrumMappings().selfConjugateRows)
        for (std::size_t z = 0; z < c.Nz; ++z) {
            half[z+c.Nz*0+c.Nz*3*row].imag = 0.0;
            half[z+c.Nz*1+c.Nz*3*row].imag = 0.0;
            half[z+c.Nz*2+c.Nz*3*row].imag = 0.0;
        }
    metrics_.scalarVerticalDerivativeSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    stageStarted = std::chrono::steady_clock::now();
    status = scalarPlan(horizontalInverse3)->execute(half,realScratch_.data()+3*R);
    if (!status) return status;
    metrics_.scalarInverseSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
    }
    const double* U = advectionFields.data;
    const double* V = U+R;
    const double* W = V+R;
    const double* dx = realScratch_.data()+3*R;
    const double* dy = dx+R;
    const double* dz = dy+R;
    stageStarted = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < R; ++index)
        rightHandSide.data[index] = -(U[index]*dx[index]+V[index]*dy[index]+W[index]*dz[index]);
    metrics_.scalarProductSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    ++metrics_.scalarAdvectionCount;
    if (shouldAntialias) {
        stageStarted = std::chrono::steady_clock::now();
        status = antialiasScalarInPlace(rightHandSide);
        metrics_.scalarAntialiasSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now()-stageStarted).count();
        return status;
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFlux(
    const WVState& state,WVFlux& flux,WVStateDerivativeAccess* derivativeAccess) {
    auto status = validateStateAndFluxForCall(state,flux);
    if (!status) return status;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,nullptr,false,nullptr,true,
        stateEvaluationActive_ && preparedPhaseReady_ ? phaseForPreparedState() : WVComplexConstView{},
        derivativeAccess);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxWithAdvectionFields(
    const WVState& state,WVFlux& flux,WVRealFieldBundleView& advectionFields,
    WVStateDerivativeAccess* derivativeAccess) {
    auto status = validateStateAndFluxForCall(state,flux);
    if (!status) return status;
    status = validateBundle(advectionFields,descriptor_.spatialShape(),3,"Advection fields");
    if (!status) return status;
    if (const auto ownership = validateInverseOwnership(state,advectionFields,descriptor_.spectralShape()); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,&advectionFields,false,nullptr,true,
        stateEvaluationActive_ && preparedPhaseReady_ ? phaseForPreparedState() : WVComplexConstView{},
        derivativeAccess);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxUsingAdvectionFields(
    const WVState& state, WVFlux& flux, const WVRealFieldBundleConstView& advectionFields,
    WVRealFieldBundleView* spatialTendency,bool projectFlux,
    WVStateDerivativeAccess* derivativeAccess) {
    if (!projectFlux && !spatialTendency)
        return {WVKernelStatusCode::invalidConfiguration,"Spatial-only nonlinear evaluation requires output storage."};
    auto status = validateStateAndFluxForCall(state,flux);
    if (!status) return status;
    status = validateBundle(advectionFields,descriptor_.spatialShape(),3,"Prepared advection fields");
    if (!status) return status;
    WVRealFieldBundleView mutableView{const_cast<double*>(advectionFields.data),advectionFields.shape};
    if (const auto ownership = validateInverseOwnership(state,mutableView,descriptor_.spectralShape()); !ownership) return ownership;
    if (spatialTendency) {
        status = validateBundle(*spatialTendency,descriptor_.spatialShape(),
            descriptor_.configuration().isHydrostatic ? 3 : 4,"Observed spatial tendency");
        if (!status) return status;
        status = validateInverseOwnership(state,*spatialTendency,descriptor_.spectralShape());
        if (!status) return status;
        const auto bytes = spatialTendency->shape.elementCount()*sizeof(double);
        if (memoryOverlaps(spatialTendency->data,bytes,advectionFields.data,
                advectionFields.shape.elementCount()*sizeof(double)))
            return {WVKernelStatusCode::overlappingArrays,"Observed tendency overlaps advection fields."};
        for (const auto output : {flux.Fp,flux.Fm,flux.F0})
            if (memoryOverlaps(spatialTendency->data,bytes,output.data,
                    descriptor_.spectralShape().elementCount()*sizeof(WVComplex64)))
                return {WVKernelStatusCode::overlappingArrays,"Observed tendency overlaps spectral output."};
    }
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,&mutableView,true,spatialTendency,projectFlux,
        stateEvaluationActive_ && preparedPhaseReady_ ? phaseForPreparedState() : WVComplexConstView{},
        derivativeAccess);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxImpl(
    const WVState& state, WVFlux& flux, WVRealFieldBundleView* retainedAdvectionFields, bool advectionFieldsPrepared,
    WVRealFieldBundleView* spatialTendency,bool projectFlux,WVComplexConstView preparedPhaseValues,
    WVStateDerivativeAccess* derivativeAccess) {
    WVKernelStatus status = WVKernelStatus::ok();

    const auto spectral = descriptor_.spectralShape();
    const auto phaseEvaluationCount = spectral.elementCount();
    using Clock = std::chrono::steady_clock;
    auto stageStart = Clock::now();
    for (std::size_t index = 0; index < phaseEvaluationCount; ++index) {
        if (projectFlux) {
            flux.Fp.data[index] = {};
            flux.Fm.data[index] = {};
            flux.F0.data[index] = {};
        }
    }
    if (stageInstrumentationEnabled_) metrics_.phaseSeconds = compact_ ? 0.0 : std::chrono::duration<double>(Clock::now() - stageStart).count();
    ++metrics_.nonlinearFluxCallCount;
    metrics_.nonlinearFluxPhaseEvaluationCount += phaseEvaluationCount;
    WVComplexConstView phaseValues=preparedPhaseValues;

    const auto spatial = descriptor_.spatialShape();
    const auto fieldElements = spatial.elementCount();
    WVRealFieldBundleView internalAdvectionFields{realScratch_.data(),{spatial.first,spatial.second,spatial.third,3}};
    auto& advectingFields = retainedAdvectionFields == nullptr ? internalAdvectionFields : *retainedAdvectionFields;
    stageStart = Clock::now();
    if (!advectionFieldsPrepared) {
        const bool generatePhase=phaseValues.data==nullptr;
        status = transformWaveVortexToUVWImpl(state,advectingFields,nullptr,phaseValues,
            generatePhase ? preparedPhase_.data() : nullptr);
        if (!status) return status;
        if (generatePhase) {
            ++metrics_.phasePreparationCount;
            if (stateEvaluationActive_) preparedPhaseReady_=true;
            phaseValues=phaseForPreparedState();
        }
        ++metrics_.advectionVelocityReconstructionCount;
    } else if (phaseValues.data==nullptr) {
        status=prepareStatePhase(state);
        if (!status) return status;
        phaseValues=phaseForPreparedState();
    }
    if (stageInstrumentationEnabled_) metrics_.reconstructionSeconds = std::chrono::duration<double>(Clock::now() - stageStart).count();

    const std::size_t targetCount = descriptor_.configuration().isHydrostatic ? 3 : 4;
    const std::size_t hydrostaticTargets[] = {0,1,3};
    const std::size_t nonhydrostaticTargets[] = {0,1,2,3};
    const auto* targets = descriptor_.configuration().isHydrostatic ? hydrostaticTargets : nonhydrostaticTargets;
    double derivativeSeconds = 0.0;
    double productSeconds = 0.0;
    double projectionSeconds = 0.0;
    for (std::size_t iTarget = 0; iTarget < targetCount; ++iTarget) {
        double* derivativeData = realScratch_.data() + 3 * fieldElements;
        WVRealFieldBundleView derivatives{derivativeData,{spatial.first,spatial.second,spatial.third,3}};
        std::array<bool,3> missing{true,true,true};
        std::array<const double*,3> derivativeValues{};
        if (derivativeAccess && derivativeAccess->lookup) for (std::size_t axis=0;axis<3;++axis) {
            WVRealVolumeConstView cached{};
            status=derivativeAccess->lookup(derivativeAccess->context,targets[iTarget],axis+1,cached);
            if (!status) return status;
            if (!cached.data) continue;
            if (cached.shape.first!=spatial.first || cached.shape.second!=spatial.second ||
                cached.shape.third!=spatial.third)
                return {WVKernelStatusCode::invalidShape,"Cached constant derivative has an incompatible shape."};
            WVRealFieldBundleView cachedBundle{const_cast<double*>(cached.data),
                {spatial.first,spatial.second,spatial.third,1}};
            status=validateBundle(cachedBundle,spatial,1,"Cached constant derivative");
            if (!status) return status;
            if (memoryOverlaps(cached.data,fieldElements*sizeof(double),derivativeData,
                    fieldElements*sizeof(double)))
                return {WVKernelStatusCode::overlappingArrays,"Cached constant derivative overlaps nonlinear product storage."};
            missing[axis]=false;
            derivativeValues[axis]=cached.data;
        }
        stageStart = Clock::now();
        if (missing[0] || missing[1] || missing[2]) {
            status=transformToSpatialDomainWithDerivativesFromStateImpl(
                state,phaseValues,targets[iTarget],derivatives,nullptr,&missing);
            if (!status) return status;
            for (std::size_t axis=0;axis<3;++axis) if (missing[axis]) {
                derivativeValues[axis]=derivativeData+axis*fieldElements;
                if (derivativeAccess && derivativeAccess->capture) {
                    status=derivativeAccess->capture(derivativeAccess->context,
                        targets[iTarget],axis+1,
                        {derivativeValues[axis],spatial});
                    if (!status) return status;
                }
            }
        }
        if (stageInstrumentationEnabled_) derivativeSeconds += std::chrono::duration<double>(Clock::now() - stageStart).count();
        const double* WV_KERNEL_RESTRICT U = advectingFields.data;
        const double* WV_KERNEL_RESTRICT V = U + fieldElements;
        const double* WV_KERNEL_RESTRICT W = V + fieldElements;
        double* WV_KERNEL_RESTRICT dx = derivativeData;
        const double* WV_KERNEL_RESTRICT sourceDx=derivativeValues[0];
        const double* WV_KERNEL_RESTRICT sourceDy=derivativeValues[1];
        const double* WV_KERNEL_RESTRICT sourceDz=derivativeValues[2];
        stageStart = Clock::now();
        auto product = [&](std::size_t begin, std::size_t end) {
            for (std::size_t i=begin;i<end;++i)
                dx[i]=-(U[i]*sourceDx[i]+V[i]*sourceDy[i]+W[i]*sourceDz[i]);
        };
        if (compact_ && fieldElements >= 4096) compact_->pointwise->execute(fieldElements,product);
        else product(0,fieldElements);
        if (spatialTendency)
            std::copy_n(dx,fieldElements,spatialTendency->data+iTarget*fieldElements);
        if (stageInstrumentationEnabled_) productSeconds += std::chrono::duration<double>(Clock::now() - stageStart).count();
        const WVRealFieldBundleConstView field{dx,{spatial.first,spatial.second,spatial.third,1}};
        stageStart = Clock::now();
        if (projectFlux) {
            status = projectSingleFluxTargetImpl(field,targets[iTarget],phaseValues,flux);
            if (!status) return status;
        }
        if (stageInstrumentationEnabled_) projectionSeconds += std::chrono::duration<double>(Clock::now() - stageStart).count();
    }
    if (stageInstrumentationEnabled_) {
        metrics_.derivativeReconstructionSeconds = derivativeSeconds;
        metrics_.productSeconds = productSeconds;
    }
    if (stageInstrumentationEnabled_) metrics_.projectionSeconds = projectionSeconds;
    return WVKernelStatus::ok();
}

} // namespace wavevortex
