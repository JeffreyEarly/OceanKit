#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"
#include "WVCoefficientFormulas.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <new>
#include <thread>
#include <stdexcept>

namespace wavevortex {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

#ifndef WV_KERNEL_COEFFICIENT_WORKERS
#define WV_KERNEL_COEFFICIENT_WORKERS 2
#endif

#if defined(WV_KERNEL_NATIVE_OPTIMIZATION) && WV_KERNEL_NATIVE_OPTIMIZATION && (defined(__clang__) || defined(__GNUC__))
#define WV_KERNEL_RESTRICT __restrict__
#else
#define WV_KERNEL_RESTRICT
#endif

template <typename Operation>
void forEachModeBlock(std::size_t modeCount, Operation&& operation) {
    constexpr std::size_t requestedWorkers = WV_KERNEL_COEFFICIENT_WORKERS;
    const std::size_t workerCount = std::min(requestedWorkers,modeCount);
    if (workerCount <= 1) {
        operation(0,modeCount);
        return;
    }
    std::vector<std::thread> workers;
    workers.reserve(workerCount - 1);
    const std::size_t blockSize = (modeCount + workerCount - 1) / workerCount;
    for (std::size_t worker = 1; worker < workerCount; ++worker) {
        const std::size_t begin = std::min(worker * blockSize,modeCount);
        const std::size_t end = std::min(begin + blockSize,modeCount);
        workers.emplace_back([begin,end,&operation]() { operation(begin,end); });
    }
    operation(0,std::min(blockSize,modeCount));
    for (auto& worker : workers) worker.join();
}

enum PlanIndex : std::size_t {
    horizontalForward3, horizontalForward4, horizontalInverse3, horizontalInverse4,
    verticalDCT2Storage3, verticalDST1Storage3, verticalDST2Storage3, verticalDCT1Storage3,
    verticalDCT2Storage4, verticalDST2Storage4,
    verticalDCT3Storage4, verticalDST3Storage4,
    verticalDCT1Storage4, verticalDST1Storage4,
    horizontalForward1, verticalDCT1Storage1, verticalDST1Storage1,
    planCount
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
    CoefficientSource&& source) {
    const auto storageRow = mapping.storageRowsByWVIndex[mode];
    const auto& horizontal = horizontalModes[mode];
    constexpr bool cosine = Target < 2;
    for (std::size_t j = 0; j < Nj; ++j) {
        const auto index = j + Nj * mode;
        const auto value = coefficientValueForField<Target>(modes,index,source(index));
        const auto derivatives = detail::normalizedDerivativeSpectrum<cosine,Conjugated>(value,horizontal.k,horizontal.l,modes.verticalWavenumber[j],j == 0 || j + 1 == Nz);
        half[j + Nz * 0 + Nz * 3 * storageRow] = derivatives.x;
        half[j + Nz * 1 + Nz * 3 * storageRow] = derivatives.y;
        half[j + Nz * 2 + Nz * 3 * storageRow] = derivatives.z;
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

WVKernelStatus WVTransformConstantStratificationKernel::create(
    const WVTransformConstantStratificationConfiguration& configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVTransformConstantStratificationKernel>& kernel) {
    if (!engine) return {WVKernelStatusCode::invalidPointer, "An FFT engine is required."};
    try {
        auto candidate = std::unique_ptr<WVTransformConstantStratificationKernel>(new WVTransformConstantStratificationKernel());
        auto status = WVTransformConstantStratificationDescriptor::create(configuration, candidate->descriptor_);
        if (!status) return status;
        candidate->engineIdentifier_ = engine->identifier();
        candidate->engineLibraryIdentity_ = engine->libraryIdentity();
        candidate->engine_ = std::move(engine);
        constexpr std::size_t halfChannels = 4;
        const auto halfElements = checkedProduct(checkedProduct(candidate->descriptor_.halfSpectrumMappings().NxHalf, configuration.Ny), checkedProduct(configuration.Nz, halfChannels));
        const auto halfFieldElements = checkedProduct(checkedProduct(candidate->descriptor_.halfSpectrumMappings().NxHalf,configuration.Ny),configuration.Nz);
        if (candidate->descriptor_.spectralShape().elementCount() > halfFieldElements) return {WVKernelStatusCode::invalidShape,"The streamed phase reservation requires M <= H."};
        constexpr std::size_t realChannels = 6;
        const auto realElements = checkedProduct(candidate->descriptor_.spatialShape().elementCount(),realChannels);
        candidate->halfSpectrumScratch_.resize(2 * halfElements);
        candidate->realScratch_.resize(realElements);
        candidate->metrics_.halfSpectrumScratchCapacityBytes = candidate->halfSpectrumScratch_.capacity() * sizeof(double);
        candidate->metrics_.realScratchCapacityBytes = candidate->realScratch_.capacity() * sizeof(double);
        candidate->metrics_.scratchCapacityBytes = candidate->metrics_.halfSpectrumScratchCapacityBytes + candidate->metrics_.realScratchCapacityBytes;
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
        candidate->metrics_.engineBytes = candidate->engine_->persistentBytes();
        candidate->metrics_.kernelManagementBytes =
            sizeof(*candidate) - sizeof(candidate->descriptor_) +
            candidate->engineIdentifier_.capacity() +
            candidate->engineLibraryIdentity_.capacity() +
            candidate->plans_.capacity() *
                sizeof(std::unique_ptr<WVFFTPlan>);
        kernel = std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure, "Unable to allocate the bounded kernel scratch arena."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow, error.what()};
    }
}

WVKernelStatus WVTransformConstantStratificationKernel::preparePlans() {
    plans_.resize(planCount);
    const auto& c = descriptor_.configuration();
    const auto halfRows = (c.Nx / 2 + 1) * c.Ny;
    const WVFFTPlanSpecification specifications[] = {
        horizontalSpecification(c, 3, true), horizontalSpecification(c, 4, true), horizontalSpecification(c, 3, false), horizontalSpecification(c, 4, false),
        verticalSpecification(c,halfRows,3,2,false), verticalSpecification(c,halfRows,3,1,true), verticalSpecification(c,halfRows,3,2,true), verticalSpecification(c,halfRows,3,1,false),
        verticalSpecification(c,halfRows,4,2,false), verticalSpecification(c,halfRows,4,2,true),
        verticalSpecification(c,halfRows,4,3,false), verticalSpecification(c,halfRows,4,3,true),
        verticalSpecification(c,halfRows,4,1,false), verticalSpecification(c,halfRows,4,1,true)
        ,horizontalSpecification(c,1,true), verticalSpecification(c,halfRows,1,1,false), verticalSpecification(c,halfRows,1,1,true)
    };
    for (std::size_t i = 0; i < planCount; ++i) {
        auto status = engine_->createPlan(specifications[i], plans_[i]);
        if (!status || !plans_[i]) return status ? WVKernelStatus{WVKernelStatusCode::fftPlanFailure, "FFT engine returned an empty plan."} : status;
        metrics_.planBytes += plans_[i]->persistentBytes();
    }
    metrics_.planCount = plans_.size();
    return WVKernelStatus::ok();
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
    return metrics_;
}

const char* WVTransformConstantStratificationKernel::nonlinearFluxScheduleIdentifier() const noexcept {
    return "streamed-target-three-channel";
}

const char* WVTransformConstantStratificationKernel::phaseImplementationIdentifier() const noexcept {
    return "scalar-sincos";
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

std::size_t WVTransformConstantStratificationKernel::coefficientWorkerCount() const noexcept { return WV_KERNEL_COEFFICIENT_WORKERS; }

std::size_t WVTransformConstantStratificationKernel::phaseReservationBytes() const noexcept {
    const auto& c = descriptor_.configuration();
    return descriptor_.halfSpectrumMappings().NxHalf * c.Ny * c.Nz * sizeof(WVComplex64);
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

WVKernelStatus WVTransformConstantStratificationKernel::transformUVEtaToWaveVortex(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients) {
    if (!descriptor_.configuration().isHydrostatic) return {WVKernelStatusCode::invalidConfiguration, "transformUVEtaToWaveVortex requires a hydrostatic kernel."};
    const auto status = validateBundle(fields, descriptor_.spatialShape(), 3, "Hydrostatic fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus outputStatuses[] = {validateSpectral(coefficients.Ap, spectral, "Ap"), validateSpectral(coefficients.Am, spectral, "Am"), validateSpectral(coefficients.A0, spectral, "A0")};
    for (const auto& value : outputStatuses) if (!value) return value;
    if (const auto ownership = validateForwardOwnership(fields,coefficients,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    return transformUVEtaToWaveVortexImpl(fields,t,t0,coefficients);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVEtaToWaveVortexImpl(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = descriptor_.halfSpectrumMappings();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward3]->execute(fields.data, half);
    if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    execute = plans_[verticalDCT2Storage3]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST1Storage3]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeForwardDCT(half,c.Nz,halfRows,3,0,2);
    normalizeForwardDST(half,c.Nz,halfRows,3,2,1);
    const double horizontalScale = 1.0 / static_cast<double>(c.Nx * c.Ny);
    auto* WV_KERNEL_RESTRICT outputAp=coefficients.Ap.data; auto* WV_KERNEL_RESTRICT outputAm=coefficients.Am.data; auto* outputA0=coefficients.A0.data;

    const auto coefficientProjectionStart = std::chrono::steady_clock::now();
    auto projectModes = [&](auto phaseProvidedConstant) {
        constexpr bool phaseProvided = decltype(phaseProvidedConstant)::value;
        forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
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

WVKernelStatus WVTransformConstantStratificationKernel::transformUVWEtaToWaveVortex(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients) {
    if (descriptor_.configuration().isHydrostatic) return {WVKernelStatusCode::invalidConfiguration, "transformUVWEtaToWaveVortex requires a nonhydrostatic kernel."};
    auto status = validateBundle(fields, descriptor_.spatialShape(), 4, "Nonhydrostatic fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus outputStatuses[] = {validateSpectral(coefficients.Ap, spectral, "Ap"), validateSpectral(coefficients.Am, spectral, "Am"), validateSpectral(coefficients.A0, spectral, "A0")};
    for (const auto& value : outputStatuses) if (!value) return value;
    if (const auto ownership = validateForwardOwnership(fields,coefficients,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    return transformUVWEtaToWaveVortexImpl(fields,t,t0,coefficients);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformUVWEtaToWaveVortexImpl(
    const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues) {
    const auto& c = descriptor_.configuration(); const auto& mapping = descriptor_.halfSpectrumMappings(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward4]->execute(fields.data, half); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    execute = plans_[verticalDCT2Storage4]->execute(half,half); if (!execute) return execute;
    execute = plans_[verticalDST2Storage4]->execute(half + 2 * c.Nz + 1,half + 2 * c.Nz + 1); if (!execute) return execute;
    metrics_.executionCount += 2; metrics_.verticalExecutionCount += 2;
    normalizeForwardDCT(half,c.Nz,halfRows,4,0,2); normalizeForwardDST(half,c.Nz,halfRows,4,2,2);
    const double horizontalScale = 1.0 / static_cast<double>(c.Nx * c.Ny);
    auto* WV_KERNEL_RESTRICT outputAp=coefficients.Ap.data; auto* WV_KERNEL_RESTRICT outputAm=coefficients.Am.data; auto* outputA0=coefficients.A0.data;
    const auto coefficientProjectionStart = std::chrono::steady_clock::now();
    auto projectModes = [&](auto phaseProvidedConstant) {
        constexpr bool phaseProvided = decltype(phaseProvidedConstant)::value;
        forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
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
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus inputStatuses[] = {validateSpectral(state.coefficients.Ap, spectral, "Ap"), validateSpectral(state.coefficients.Am, spectral, "Am"), validateSpectral(state.coefficients.A0, spectral, "A0")};
    for (const auto& value : inputStatuses) if (!value) return value;
    if (const auto ownership = validateInverseOwnership(state,fields,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    return transformWaveVortexToUVWEtaImpl(state,fields);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVW(const WVState& state, WVRealFieldBundleView& fields) {
    auto status = validateBundle(fields,descriptor_.spatialShape(),3,"Advection fields");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus inputStatuses[] = {
        validateSpectral(state.coefficients.Ap,spectral,"Ap"),
        validateSpectral(state.coefficients.Am,spectral,"Am"),
        validateSpectral(state.coefficients.A0,spectral,"A0")};
    for (const auto& value : inputStatuses) if (!value) return value;
    if (const auto ownership = validateInverseOwnership(state,fields,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration();
    const auto halfFieldElements = descriptor_.halfSpectrumMappings().NxHalf * c.Ny * c.Nz;
    auto* phaseStorage = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data()) + 3 * halfFieldElements;
    const auto count = spectral.elementCount();
    const auto& omega = descriptor_.verticalModes().omega;
    const double elapsed = state.t-state.t0;
    for (std::size_t index = 0; index < count; ++index) phaseStorage[index] = phase(omega[index]*elapsed);
    ++metrics_.advectionVelocityReconstructionCount;
    return transformWaveVortexToUVWImpl(state,fields,nullptr,{phaseStorage,spectral});
}

WVKernelStatus WVTransformConstantStratificationKernel::transformStateFieldDerivatives(
    const WVState& state, WVDynamicalField field, WVRealFieldBundleView& derivatives) {
    const auto target = static_cast<std::size_t>(field);
    if (target > static_cast<std::size_t>(WVDynamicalField::eta))
        return {WVKernelStatusCode::invalidConfiguration,"Unknown dynamical field derivative target."};
    auto status = validateBundle(derivatives,descriptor_.spatialShape(),3,"Field derivatives");
    if (!status) return status;
    const auto spectral = descriptor_.spectralShape();
    const WVKernelStatus inputStatuses[] = {
        validateSpectral(state.coefficients.Ap,spectral,"Ap"),
        validateSpectral(state.coefficients.Am,spectral,"Am"),
        validateSpectral(state.coefficients.A0,spectral,"A0")};
    for (const auto& value : inputStatuses) if (!value) return value;
    if (const auto ownership = validateInverseOwnership(state,derivatives,spectral); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration();
    const auto halfFieldElements = descriptor_.halfSpectrumMappings().NxHalf * c.Ny * c.Nz;
    auto* phaseStorage = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data()) + 3 * halfFieldElements;
    const auto coefficientCount = spectral.elementCount();
    const auto& omega = descriptor_.verticalModes().omega;
    const double elapsed = state.t-state.t0;
    for (std::size_t index = 0; index < coefficientCount; ++index) phaseStorage[index] = phase(omega[index]*elapsed);
    return transformToSpatialDomainWithDerivativesFromStateImpl(state,{phaseStorage,spectral},target,derivatives);
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVWEtaImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients) {
    const auto& c = descriptor_.configuration(); const auto& mapping = descriptor_.halfSpectrumMappings(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny; auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t iMode = begin; iMode < end; ++iMode) {
        if (evolvedCoefficients == nullptr) {
            auto source = [&](std::size_t index) {
                const auto p = phase(modes.omega[index] * (state.t - state.t0));
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
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::transformWaveVortexToUVWImpl(
    const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients, WVComplexConstView phaseValues) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = descriptor_.halfSpectrumMappings();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
        for (std::size_t iMode = begin; iMode < end; ++iMode) {
            auto source = [&](std::size_t index) {
                if (evolvedCoefficients != nullptr) return EvolvedWaveVortexCoefficients{evolvedCoefficients->Ap.data[index],evolvedCoefficients->Am.data[index],evolvedCoefficients->A0.data[index]};
                return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValues.data[index]);
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
    return WVKernelStatus::ok();
}

namespace {
WVKernelStatus transformAllDerivatives(WVTransformConstantStratificationKernel& kernel, bool cosine, const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    const auto spectral = kernel.descriptor().spectralShape();
    auto status = validateSpectral(Apm, spectral, "Apm"); if (!status) return status;
    status = validateSpectral(A0, spectral, "A0"); if (!status) return status;
    status = validateBundle(fields, kernel.descriptor().spatialShape(), 4, "Field derivatives"); if (!status) return status;
    status = validateDerivativeOwnership(Apm,A0,fields,spectral); if (!status) return status;
    // Implemented as a friendless helper through the public transform path is not possible; the member wrappers below contain the operation.
    return {WVKernelStatusCode::unsupportedOperation, cosine ? "internal F derivative dispatch" : "internal G derivative dispatch"};
}
} // namespace

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithFAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    auto initial = transformAllDerivatives(*this, true, Apm, A0, fields); if (initial.code != WVKernelStatusCode::unsupportedOperation) return initial;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration(); const auto& mapping = descriptor_.halfSpectrumMappings(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny; auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
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

WVKernelStatus WVTransformConstantStratificationKernel::transformToSpatialDomainWithGAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields) {
    auto initial = transformAllDerivatives(*this, false, Apm, A0, fields); if (initial.code != WVKernelStatusCode::unsupportedOperation) return initial;
    ExecutionGuard guard(executing_); if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution, "Kernel operations are not reentrant."};
    const auto& c = descriptor_.configuration(); const auto& mapping = descriptor_.halfSpectrumMappings(); const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny; auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 4 * halfRows,WVComplex64{});
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
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
    const auto& mapping = descriptor_.halfSpectrumMappings();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const bool cosine = target < 2;
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    auto assembleTarget = [&](auto targetConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t mode = begin; mode < end; ++mode) {
                auto source = [&](std::size_t index) { return EvolvedWaveVortexCoefficients{evolvedCoefficients.Ap.data[index],evolvedCoefficients.Am.data[index],evolvedCoefficients.A0.data[index]}; };
                if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleDerivativeSpectraForMode<resolvedTarget,true>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,source);
                else assembleDerivativeSpectraForMode<resolvedTarget,false>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,source);
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
    const WVState& state, WVComplexConstView phaseValues, std::size_t target, WVRealFieldBundleView& derivatives) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = descriptor_.halfSpectrumMappings();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    std::fill(half,half + c.Nz * 3 * halfRows,WVComplex64{});
    const bool cosine = target < 2;
    const auto coefficientAssemblyStart = std::chrono::steady_clock::now();
    auto assembleTarget = [&](auto targetConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
            for (std::size_t mode = begin; mode < end; ++mode) {
                auto source = [&](std::size_t index) { return evolveWaveVortexCoefficients(state.coefficients.Ap.data[index],state.coefficients.Am.data[index],state.coefficients.A0.data[index],phaseValues.data[index]); };
                if (mapping.conjugatesStoredValueByWVIndex[mode]) assembleDerivativeSpectraForMode<resolvedTarget,true>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,source);
                else assembleDerivativeSpectraForMode<resolvedTarget,false>(half,mapping,modes,descriptor_.fourierModes(),c.Nz,c.Nj,mode,source);
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
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformConstantStratificationKernel::projectSingleFluxTargetImpl(
    const WVRealFieldBundleConstView& field, std::size_t target, WVComplexConstView phaseValues, WVFlux& flux) {
    const auto& c = descriptor_.configuration();
    const auto& mapping = descriptor_.halfSpectrumMappings();
    const auto& modes = descriptor_.verticalModes();
    const std::size_t halfRows = mapping.NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto execute = plans_[horizontalForward1]->execute(field.data,half); if (!execute) return execute;
    ++metrics_.executionCount; ++metrics_.horizontalExecutionCount;
    if (target < 2) {
        execute = plans_[verticalDCT1Storage1]->execute(half,half); if (!execute) return execute;
        normalizeForwardDCT(half,c.Nz,halfRows,1,0,1);
    } else {
        execute = plans_[verticalDST1Storage1]->execute(half + 1,half + 1); if (!execute) return execute;
        normalizeForwardDST(half,c.Nz,halfRows,1,0,1);
    }
    ++metrics_.executionCount; ++metrics_.verticalExecutionCount;
    const double horizontalScale = 1.0 / static_cast<double>(c.Nx * c.Ny);
    auto projectTarget = [&](auto targetConstant, auto hydrostaticConstant) {
        constexpr std::size_t resolvedTarget = decltype(targetConstant)::value;
        constexpr bool hydrostatic = decltype(hydrostaticConstant)::value;
        forEachModeBlock(descriptor_.Nkl(),[&](std::size_t begin, std::size_t end) {
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
    return ensureScalarInversePlan();
}

WVKernelStatus WVTransformConstantStratificationKernel::antialiasScalarInPlace(WVRealVolumeView& scalar) {
    auto status = ensureScalarInversePlan();
    if (!status) return status;
    const auto& c = descriptor_.configuration();
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf * c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    status = plans_[horizontalForward1]->execute(scalar.data,half);
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

    const auto& c = descriptor_.configuration();
    const std::size_t NxHalf = descriptor_.halfSpectrumMappings().NxHalf;
    const std::size_t halfRows = NxHalf*c.Ny;
    auto* half = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data());
    auto stageStarted = std::chrono::steady_clock::now();
    status = plans_[horizontalForward1]->execute(scalar.data,half);
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
    status = plans_[verticalDCT1Storage3]->execute(half+2*c.Nz,half+2*c.Nz);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.verticalExecutionCount;
    normalizeForwardDCT(half,c.Nz,halfRows,3,2,1);
    for (std::size_t row = 0; row < halfRows; ++row) {
        const auto base = 2*c.Nz+3*c.Nz*row;
        half[base] = {};
        for (std::size_t j = 1; j+1 < c.Nz; ++j) {
            const double verticalWavenumber = pi*static_cast<double>(j)/c.Lz;
            half[base+j] = multiply(half[base+j],-verticalWavenumber);
        }
        half[base+c.Nz-1] = {};
    }
    normalizeInverseDST(half,c.Nz,halfRows,3,2,1);
    status = plans_[verticalDST1Storage3]->execute(half+2*c.Nz+1,half+2*c.Nz+1);
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
    status = plans_[horizontalInverse3]->execute(half,realScratch_.data()+3*R);
    if (!status) return status;
    metrics_.scalarInverseSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now()-stageStarted).count();
    ++metrics_.executionCount;
    ++metrics_.horizontalExecutionCount;
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

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFlux(const WVState& state, WVFlux& flux) {
    auto status = validateStateAndFlux(descriptor_,state,flux);
    if (!status) return status;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,nullptr,false);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxWithAdvectionFields(
    const WVState& state, WVFlux& flux, WVRealFieldBundleView& advectionFields) {
    auto status = validateStateAndFlux(descriptor_,state,flux);
    if (!status) return status;
    status = validateBundle(advectionFields,descriptor_.spatialShape(),3,"Advection fields");
    if (!status) return status;
    if (const auto ownership = validateInverseOwnership(state,advectionFields,descriptor_.spectralShape()); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,&advectionFields,false);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxUsingAdvectionFields(
    const WVState& state, WVFlux& flux, const WVRealFieldBundleConstView& advectionFields) {
    auto status = validateStateAndFlux(descriptor_,state,flux);
    if (!status) return status;
    status = validateBundle(advectionFields,descriptor_.spatialShape(),3,"Prepared advection fields");
    if (!status) return status;
    WVRealFieldBundleView mutableView{const_cast<double*>(advectionFields.data),advectionFields.shape};
    if (const auto ownership = validateInverseOwnership(state,mutableView,descriptor_.spectralShape()); !ownership) return ownership;
    ExecutionGuard guard(executing_);
    if (!guard.entered()) return {WVKernelStatusCode::reentrantExecution,"Kernel operations are not reentrant."};
    return nonlinearFluxImpl(state,flux,&mutableView,true);
}

WVKernelStatus WVTransformConstantStratificationKernel::nonlinearFluxImpl(
    const WVState& state, WVFlux& flux, WVRealFieldBundleView* retainedAdvectionFields, bool advectionFieldsPrepared) {
    WVKernelStatus status = WVKernelStatus::ok();

    const auto spectral = descriptor_.spectralShape();
    const auto phaseEvaluationCount = spectral.elementCount();
    const auto& modes = descriptor_.verticalModes();
    using Clock = std::chrono::steady_clock;
    auto stageStart = Clock::now();
    const auto& c = descriptor_.configuration();
    const auto halfFieldElements = descriptor_.halfSpectrumMappings().NxHalf * c.Ny * c.Nz;
    auto* phaseStorage = reinterpret_cast<WVComplex64*>(halfSpectrumScratch_.data()) + 3 * halfFieldElements;
    const double elapsed = state.t - state.t0;
    for (std::size_t index = 0; index < phaseEvaluationCount; ++index) {
        phaseStorage[index] = phase(modes.omega[index] * elapsed);
        flux.Fp.data[index] = {};
        flux.Fm.data[index] = {};
        flux.F0.data[index] = {};
    }
    if (stageInstrumentationEnabled_) metrics_.phaseSeconds = std::chrono::duration<double>(Clock::now() - stageStart).count();
    ++metrics_.nonlinearFluxCallCount;
    metrics_.nonlinearFluxPhaseEvaluationCount += phaseEvaluationCount;
    const WVComplexConstView phaseValues{phaseStorage,spectral};

    const auto spatial = descriptor_.spatialShape();
    const auto fieldElements = spatial.elementCount();
    WVRealFieldBundleView internalAdvectionFields{realScratch_.data(),{spatial.first,spatial.second,spatial.third,3}};
    auto& advectingFields = retainedAdvectionFields == nullptr ? internalAdvectionFields : *retainedAdvectionFields;
    stageStart = Clock::now();
    if (!advectionFieldsPrepared) {
        status = transformWaveVortexToUVWImpl(state,advectingFields,nullptr,phaseValues);
        if (!status) return status;
        ++metrics_.advectionVelocityReconstructionCount;
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
        stageStart = Clock::now();
        status = transformToSpatialDomainWithDerivativesFromStateImpl(state,phaseValues,targets[iTarget],derivatives);
        if (!status) return status;
        if (stageInstrumentationEnabled_) derivativeSeconds += std::chrono::duration<double>(Clock::now() - stageStart).count();
        const double* WV_KERNEL_RESTRICT U = advectingFields.data;
        const double* WV_KERNEL_RESTRICT V = U + fieldElements;
        const double* WV_KERNEL_RESTRICT W = V + fieldElements;
        double* WV_KERNEL_RESTRICT dx = derivativeData;
        const double* WV_KERNEL_RESTRICT dy = dx + fieldElements;
        const double* WV_KERNEL_RESTRICT dz = dy + fieldElements;
        stageStart = Clock::now();
        for (std::size_t i = 0; i < fieldElements; ++i) dx[i] = -(U[i] * dx[i] + V[i] * dy[i] + W[i] * dz[i]);
        if (stageInstrumentationEnabled_) productSeconds += std::chrono::duration<double>(Clock::now() - stageStart).count();
        const WVRealFieldBundleConstView field{dx,{spatial.first,spatial.second,spatial.third,1}};
        stageStart = Clock::now();
        status = projectSingleFluxTargetImpl(field,targets[iTarget],phaseValues,flux);
        if (!status) return status;
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
