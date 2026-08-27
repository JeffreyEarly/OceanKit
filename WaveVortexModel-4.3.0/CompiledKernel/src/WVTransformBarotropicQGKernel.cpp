#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wavevortex {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::size_t forwardPlan = 0;
constexpr std::size_t inversePlan = 1;
constexpr std::size_t nonlinearInversePlan = 2;
constexpr std::size_t planCount = 3;

bool finitePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

std::size_t checkedProduct(std::size_t first, std::size_t second) {
    if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::overflow_error("shape element count overflow");
    }
    return first * second;
}

double radialMagnitude(double first, double second) {
    volatile double firstSquared = first * first;
    volatile double secondSquared = second * second;
    return std::sqrt(firstSquared + secondSquared);
}

std::vector<std::int64_t> dftModes(std::size_t count) {
    std::vector<std::int64_t> modes;
    modes.reserve(count);
    const auto positiveCount = (count + 1) / 2;
    for (std::size_t value = 0; value < positiveCount; ++value)
        modes.push_back(static_cast<std::int64_t>(value));
    for (std::size_t value = count / 2; value > 0; --value)
        modes.push_back(-static_cast<std::int64_t>(value));
    return modes;
}

std::size_t dftIndexForMode(std::int64_t mode, std::size_t count) {
    const auto wrapped = mode >= 0 ? mode : static_cast<std::int64_t>(count) + mode;
    return static_cast<std::size_t>(wrapped);
}

bool isSelfConjugate(std::int64_t mode, std::size_t count) {
    return mode == 0 ||
           (count % 2 == 0 &&
            mode == -static_cast<std::int64_t>(count / 2));
}

bool isPrimary(std::int64_t kMode, std::int64_t lMode,
               std::size_t Nx, std::size_t Ny) {
    const bool kSelf = isSelfConjugate(kMode, Nx);
    const bool lSelf = isSelfConjugate(lMode, Ny);
    return lMode > 0 || (lSelf && (kMode > 0 || kSelf));
}

bool isNyquist(std::int64_t kMode, std::int64_t lMode,
                std::size_t Nx, std::size_t Ny) {
    return (Nx % 2 == 0 &&
            kMode == -static_cast<std::int64_t>(Nx / 2)) ||
           (Ny % 2 == 0 &&
            lMode == -static_cast<std::int64_t>(Ny / 2));
}

std::size_t halfRow(std::int64_t kMode, std::int64_t lMode,
                    std::size_t NxHalf, std::size_t Ny) {
    return static_cast<std::size_t>(kMode) +
           NxHalf * dftIndexForMode(lMode, Ny);
}

template <typename T>
std::size_t bytes(const std::vector<T>& values) {
    return values.capacity() * sizeof(T);
}

WVComplex64 conjugate(WVComplex64 value) {
    return {value.real, -value.imag};
}

WVComplex64 multiply(WVComplex64 value, double scale) {
    return {value.real * scale, value.imag * scale};
}

WVComplex64 multiply(WVComplex64 first, WVComplex64 second) {
    return {first.real * second.real - first.imag * second.imag,
            first.real * second.imag + first.imag * second.real};
}

WVKernelStatus validateSpectral(const WVComplexConstView& view,
                                WVShape2D expected,
                                const char* name) {
    if (view.shape.rows != expected.rows ||
        view.shape.columns != expected.columns)
        return {WVKernelStatusCode::invalidShape,
                std::string(name) + " must have shape [1,Nkl]."};
    if (view.data == nullptr && expected.elementCount() != 0)
        return {WVKernelStatusCode::invalidPointer,
                std::string(name) + " has a null data pointer."};
    return WVKernelStatus::ok();
}

WVKernelStatus validateSpectral(const WVComplexView& view,
                                WVShape2D expected,
                                const char* name) {
    const WVComplexConstView constView{view.data, view.shape};
    return validateSpectral(constView, expected, name);
}

WVKernelStatus validateSpatial(const WVRealConstView& view,
                               WVShape2D expected,
                               const char* name) {
    if (view.shape.rows != expected.rows ||
        view.shape.columns != expected.columns)
        return {WVKernelStatusCode::invalidShape,
                std::string(name) + " must have shape [Nx,Ny]."};
    if (view.data == nullptr && expected.elementCount() != 0)
        return {WVKernelStatusCode::invalidPointer,
                std::string(name) + " has a null data pointer."};
    return WVKernelStatus::ok();
}

WVKernelStatus validateSpatial(const WVRealView& view,
                               WVShape2D expected,
                               const char* name) {
    const WVRealConstView constView{view.data, view.shape};
    return validateSpatial(constView, expected, name);
}

bool overlaps(const void* first, std::size_t firstBytes,
              const void* second, std::size_t secondBytes) {
    if (first == nullptr || second == nullptr ||
        firstBytes == 0 || secondBytes == 0)
        return false;
    const auto firstBegin = reinterpret_cast<std::uintptr_t>(first);
    const auto secondBegin = reinterpret_cast<std::uintptr_t>(second);
    return firstBegin < secondBegin + secondBytes &&
           secondBegin < firstBegin + firstBytes;
}

WVFFTPlanSpecification horizontalSpecification(
    const WVTransformBarotropicQGConfiguration& configuration,
    std::size_t fields, bool forward) {
    const std::size_t NxHalf = configuration.Nx / 2 + 1;
    const std::size_t realPlane = configuration.Nx * configuration.Ny;
    const std::size_t halfRows = NxHalf * configuration.Ny;
    WVFFTPlanSpecification specification;
    specification.kind = forward
        ? WVFFTPlanKind::horizontalRealToComplex2D
        : WVFFTPlanKind::horizontalComplexToReal2D;
    specification.transformDimensions = {
        {configuration.Ny, static_cast<std::ptrdiff_t>(configuration.Nx),
         static_cast<std::ptrdiff_t>(fields * NxHalf)},
        {configuration.Nx, 1, static_cast<std::ptrdiff_t>(fields)}};
    specification.batchDimensions = {
        {fields, static_cast<std::ptrdiff_t>(realPlane), 1}};
    specification.inputBytes = forward
        ? fields * realPlane * sizeof(double)
        : fields * halfRows * sizeof(WVComplex64);
    specification.outputBytes = forward
        ? fields * halfRows * sizeof(WVComplex64)
        : fields * realPlane * sizeof(double);
    specification.destroysInput = !forward;
    if (!forward) {
        for (auto& dimension : specification.transformDimensions)
            std::swap(dimension.inputStride, dimension.outputStride);
        for (auto& dimension : specification.batchDimensions)
            std::swap(dimension.inputStride, dimension.outputStride);
    }
    return specification;
}

void completeHermitianBoundaries(WVComplex64* half,
                                 const WVHalfSpectrumMappings& mappings,
                                 std::size_t fields) {
    for (std::size_t field = 0; field < fields; ++field) {
        for (std::size_t index = 0;
             index < mappings.hermitianCompletionRows.size(); ++index) {
            half[field + fields * mappings.hermitianCompletionRows[index]] =
                conjugate(half[field +
                               fields * mappings.hermitianSourceRows[index]]);
        }
        for (const auto row : mappings.selfConjugateRows)
            half[field + fields * row].imag = 0.0;
    }
}

WVComplex64 fieldFactor(const WVBarotropicQGModes& modes,
                        WVBarotropicQGField field,
                        std::size_t index) {
    switch (field) {
        case WVBarotropicQGField::u: return modes.uFactor[index];
        case WVBarotropicQGField::v: return modes.vFactor[index];
        case WVBarotropicQGField::eta:
            return {modes.etaFactor[index], 0.0};
        case WVBarotropicQGField::pi:
        case WVBarotropicQGField::ssh:
            return {modes.piFactor[index], 0.0};
        case WVBarotropicQGField::psi:
            return {modes.psiFactor[index], 0.0};
        case WVBarotropicQGField::qgpv:
            return {modes.qgpvFactor[index], 0.0};
        case WVBarotropicQGField::zetaZ:
            return {modes.zetaZFactor[index], 0.0};
    }
    return {};
}

class ExecutionGuard final {
public:
    explicit ExecutionGuard(bool& value) : value_(value) { value_ = true; }
    ~ExecutionGuard() { value_ = false; }
private:
    bool& value_;
};

} // namespace

bool sameTransformConfiguration(
    const WVTransformBarotropicQGConfiguration& first,
    const WVTransformBarotropicQGConfiguration& second) noexcept {
    return first.contractVersion == second.contractVersion &&
           first.Nx == second.Nx && first.Ny == second.Ny &&
           first.Lx == second.Lx && first.Ly == second.Ly &&
           first.h == second.h && first.j == second.j &&
           first.g == second.g &&
           first.planetaryRadius == second.planetaryRadius &&
           first.rotationRate == second.rotationRate &&
           first.latitude == second.latitude &&
           first.shouldAntialias == second.shouldAntialias;
}

std::size_t WVBarotropicQGModes::persistentBytes() const noexcept {
    return bytes(uFactor) + bytes(vFactor) + bytes(etaFactor) +
           bytes(piFactor) + bytes(psiFactor) + bytes(qgpvFactor) +
           bytes(zetaZFactor) + bytes(energyFactor) +
           bytes(enstrophyFactor);
}

WVKernelStatus WVTransformBarotropicQGDescriptor::create(
    const WVTransformBarotropicQGConfiguration& configuration,
    WVTransformBarotropicQGDescriptor& descriptor) {
    if (configuration.contractVersion != WVKernelContractVersion)
        return {WVKernelStatusCode::invalidConfiguration,
                "Unsupported Barotropic QG kernel contract version."};
    if (configuration.Nx < 2 || configuration.Ny < 2 ||
        configuration.j > 1)
        return {WVKernelStatusCode::invalidConfiguration,
                "Nx and Ny must be at least 2 and j must be 0 or 1."};
    if (!finitePositive(configuration.Lx) ||
        !finitePositive(configuration.Ly) ||
        !finitePositive(configuration.h) ||
        !finitePositive(configuration.g) ||
        !finitePositive(configuration.planetaryRadius) ||
        !finitePositive(configuration.rotationRate) ||
        !std::isfinite(configuration.latitude) ||
        std::abs(configuration.latitude) < 5.0 ||
        std::abs(configuration.latitude) > 85.0)
        return {WVKernelStatusCode::invalidConfiguration,
                "Barotropic QG lengths and physical constants must be finite "
                "and positive, and absolute latitude must lie in [5,85]."};

    try {
        checkedProduct(configuration.Nx, configuration.Ny);
        WVTransformBarotropicQGDescriptor candidate;
        candidate.configuration_ = configuration;
        const auto kModes = dftModes(configuration.Nx);
        const auto lModes = dftModes(configuration.Ny);
        const double maximumK = 2.0 * pi *
            (static_cast<double>(configuration.Nx / 2) / configuration.Lx);
        for (const auto lMode : lModes) {
            for (const auto kMode : kModes) {
                if (!isPrimary(kMode, lMode, configuration.Nx,
                               configuration.Ny) ||
                    isNyquist(kMode, lMode, configuration.Nx,
                               configuration.Ny))
                    continue;
                const double k = 2.0 * pi *
                    (static_cast<double>(kMode) / configuration.Lx);
                const double l = 2.0 * pi *
                    (static_cast<double>(lMode) / configuration.Ly);
                const double Kh = radialMagnitude(k, l);
                if (configuration.shouldAntialias &&
                    Kh > 2.0 * maximumK / 3.0)
                    continue;
                const auto iK = dftIndexForMode(kMode, configuration.Nx);
                const auto iL = dftIndexForMode(lMode, configuration.Ny);
                const auto conjugateK =
                    dftIndexForMode(-kMode, configuration.Nx);
                const auto conjugateL =
                    dftIndexForMode(-lMode, configuration.Ny);
                candidate.fourierModes_.push_back(
                    {kMode, lMode, k, l, Kh,
                     Kh == 0.0 ? 0.0 : k / Kh,
                     Kh == 0.0 ? 0.0 : l / Kh,
                     iK + configuration.Nx * iL,
                     conjugateK + configuration.Nx * conjugateL});
            }
        }
        std::stable_sort(
            candidate.fourierModes_.begin(), candidate.fourierModes_.end(),
            [](const WVFourierMode& first, const WVFourierMode& second) {
                if (first.Kh != second.Kh) return first.Kh < second.Kh;
                if (first.kMode != second.kMode)
                    return first.kMode < second.kMode;
                return first.lMode < second.lMode;
            });
        if (candidate.fourierModes_.empty())
            return {WVKernelStatusCode::invalidConfiguration,
                    "Barotropic QG configuration retained no Fourier modes."};

        auto& mappings = candidate.halfSpectrumMappings_;
        mappings.NxHalf = configuration.Nx / 2 + 1;
        const auto halfRows = checkedProduct(mappings.NxHalf, configuration.Ny);
        mappings.storageRowsByWVIndex.resize(candidate.fourierModes_.size());
        mappings.conjugatesStoredValueByWVIndex.resize(
            candidate.fourierModes_.size());
        std::vector<std::ptrdiff_t> represented(halfRows, -1);
        for (std::size_t iWV = 0; iWV < candidate.fourierModes_.size();
             ++iWV) {
            const auto& mode = candidate.fourierModes_[iWV];
            if (mode.kMode >= 0) {
                const auto row = halfRow(mode.kMode, mode.lMode,
                                         mappings.NxHalf, configuration.Ny);
                mappings.directRows.push_back(row);
                mappings.directWVIndices.push_back(iWV);
                mappings.storageRowsByWVIndex[iWV] = row;
                represented[row] = static_cast<std::ptrdiff_t>(iWV);
            } else {
                const auto row = halfRow(-mode.kMode, -mode.lMode,
                                         mappings.NxHalf, configuration.Ny);
                mappings.conjugatedRows.push_back(row);
                mappings.conjugatedWVIndices.push_back(iWV);
                mappings.storageRowsByWVIndex[iWV] = row;
                mappings.conjugatesStoredValueByWVIndex[iWV] = 1;
                represented[row] = static_cast<std::ptrdiff_t>(iWV);
            }
        }
        for (std::size_t iL = 0; iL < configuration.Ny; ++iL) {
            const std::int64_t lMode = iL <= configuration.Ny / 2
                ? static_cast<std::int64_t>(iL)
                : static_cast<std::int64_t>(iL) -
                      static_cast<std::int64_t>(configuration.Ny);
            const auto destination =
                halfRow(0, lMode, mappings.NxHalf, configuration.Ny);
            const auto source =
                halfRow(0, -lMode, mappings.NxHalf, configuration.Ny);
            if (represented[destination] < 0 && represented[source] >= 0) {
                mappings.hermitianCompletionRows.push_back(destination);
                mappings.hermitianSourceRows.push_back(source);
            }
            if ((lMode == 0 ||
                 (configuration.Ny % 2 == 0 &&
                  std::abs(lMode) ==
                      static_cast<std::int64_t>(configuration.Ny / 2))) &&
                represented[destination] >= 0)
                mappings.selfConjugateRows.push_back(destination);
        }
        std::reverse(mappings.hermitianCompletionRows.begin(),
                     mappings.hermitianCompletionRows.end());
        std::reverse(mappings.hermitianSourceRows.begin(),
                     mappings.hermitianSourceRows.end());

        auto& modes = candidate.modes_;
        modes.coriolisFrequency = 2.0 * configuration.rotationRate *
            std::sin(configuration.latitude * pi / 180.0);
        modes.deformationWavenumberSquared = configuration.j == 0
            ? 0.0
            : modes.coriolisFrequency * modes.coriolisFrequency /
                  (configuration.g * configuration.h);
        const auto count = candidate.fourierModes_.size();
        modes.uFactor.resize(count);
        modes.vFactor.resize(count);
        modes.etaFactor.resize(count);
        modes.piFactor.resize(count);
        modes.psiFactor.resize(count);
        modes.qgpvFactor.resize(count);
        modes.zetaZFactor.resize(count);
        modes.energyFactor.resize(count);
        modes.enstrophyFactor.resize(count);
        for (std::size_t index = 0; index < count; ++index) {
            const auto& mode = candidate.fourierModes_[index];
            const double K2 = mode.k * mode.k + mode.l * mode.l;
            if (K2 == 0.0) continue;
            const double denominator =
                K2 + modes.deformationWavenumberSquared;
            modes.uFactor[index] = {0.0, mode.l / denominator};
            modes.vFactor[index] = {0.0, -mode.k / denominator};
            modes.etaFactor[index] = configuration.j == 0
                ? 0.0 : -(modes.coriolisFrequency / configuration.g) /
                            denominator;
            modes.piFactor[index] =
                -(modes.coriolisFrequency / configuration.g) / denominator;
            modes.psiFactor[index] = -1.0 / denominator;
            modes.qgpvFactor[index] = 1.0;
            modes.zetaZFactor[index] = K2 / denominator;
            modes.energyFactor[index] = configuration.h / denominator;
            modes.enstrophyFactor[index] = configuration.h;
        }
        descriptor = std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,
                "Unable to allocate the Barotropic QG descriptor."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow, error.what()};
    }
}

std::size_t WVTransformBarotropicQGDescriptor::persistentBytes() const noexcept {
    return sizeof(*this) + bytes(fourierModes_) +
           halfSpectrumMappings_.persistentBytes() + modes_.persistentBytes();
}

WVKernelStatus WVTransformBarotropicQGKernel::create(
    const WVTransformBarotropicQGConfiguration& configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVTransformBarotropicQGKernel>& kernel) {
    kernel.reset();
    if (!engine)
        return {WVKernelStatusCode::invalidPointer,
                "A Barotropic QG FFT engine is required."};
    try {
        auto candidate = std::unique_ptr<WVTransformBarotropicQGKernel>(
            new WVTransformBarotropicQGKernel());
        auto status = WVTransformBarotropicQGDescriptor::create(
            configuration, candidate->descriptor_);
        if (!status) return status;
        candidate->engineIdentifier_ = engine->identifier();
        candidate->engineLibraryIdentity_ = engine->libraryIdentity();
        candidate->engine_ = std::move(engine);
        const auto halfRows = checkedProduct(
            candidate->descriptor_.halfSpectrumMappings().NxHalf,
            configuration.Ny);
        const auto realElements = checkedProduct(configuration.Nx,
                                                 configuration.Ny);
        candidate->halfSpectrumScratch_.resize(
            checkedProduct(halfRows, static_cast<std::size_t>(4)));
        candidate->realScratch_.resize(
            checkedProduct(realElements, static_cast<std::size_t>(5)));
        candidate->metrics_.descriptorBytes =
            candidate->descriptor_.persistentBytes();
        candidate->metrics_.halfSpectrumScratchCapacityBytes =
            bytes(candidate->halfSpectrumScratch_);
        candidate->metrics_.realScratchCapacityBytes =
            bytes(candidate->realScratch_);
        candidate->metrics_.scratchCapacityBytes =
            candidate->metrics_.halfSpectrumScratchCapacityBytes +
            candidate->metrics_.realScratchCapacityBytes;
        candidate->metrics_.scratchHighWaterBytes =
            candidate->metrics_.scratchCapacityBytes;
        status = candidate->preparePlans();
        if (!status) return status;
        candidate->metrics_.engineBytes =
            candidate->engine_->persistentBytes();
        candidate->metrics_.kernelManagementBytes =
            sizeof(*candidate) - sizeof(candidate->descriptor_) +
            candidate->engineIdentifier_.capacity() +
            candidate->engineLibraryIdentity_.capacity() +
            candidate->plans_.capacity() *
                sizeof(std::unique_ptr<WVFFTPlan>);
        kernel = std::move(candidate);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,
                "Unable to allocate bounded Barotropic QG scratch."};
    } catch (const std::overflow_error& error) {
        return {WVKernelStatusCode::sizeOverflow, error.what()};
    }
}

WVKernelStatus WVTransformBarotropicQGKernel::preparePlans() {
    plans_.resize(planCount);
    const auto& configuration = descriptor_.configuration();
    const WVFFTPlanSpecification specifications[] = {
        horizontalSpecification(configuration, 1, true),
        horizontalSpecification(configuration, 1, false),
        horizontalSpecification(configuration, 4, false)};
    for (std::size_t index = 0; index < planCount; ++index) {
        auto status = engine_->createPlan(specifications[index], plans_[index]);
        if (!status || !plans_[index])
            return status ? WVKernelStatus{
                                WVKernelStatusCode::fftPlanFailure,
                                "FFT engine returned an empty Barotropic QG plan."}
                          : status;
        metrics_.planBytes += plans_[index]->persistentBytes();
    }
    metrics_.planCount = plans_.size();
    return WVKernelStatus::ok();
}

const WVBarotropicQGKernelMetrics&
WVTransformBarotropicQGKernel::metrics() const noexcept {
    metrics_.engineBytes = engine_ == nullptr ? 0 : engine_->persistentBytes();
    metrics_.planBytes = 0;
    metrics_.planCount = 0;
    for (const auto& plan : plans_) {
        if (plan) {
            metrics_.planBytes += plan->persistentBytes();
            ++metrics_.planCount;
        }
    }
    return metrics_;
}

std::size_t WVTransformBarotropicQGKernel::persistentBytes() const noexcept {
    const auto& current = metrics();
    return current.descriptorBytes + current.planBytes +
           current.engineBytes + current.kernelManagementBytes +
           current.scratchCapacityBytes;
}

std::size_t WVTransformBarotropicQGKernel::scratchBytes() const noexcept {
    return bytes(halfSpectrumScratch_) + bytes(realScratch_);
}

const WVComplex64* WVTransformBarotropicQGKernel::complexFactors(
    WVBarotropicQGField field) const noexcept {
    if (field == WVBarotropicQGField::u)
        return descriptor_.modes().uFactor.data();
    if (field == WVBarotropicQGField::v)
        return descriptor_.modes().vFactor.data();
    return nullptr;
}

const double* WVTransformBarotropicQGKernel::realFactors(
    WVBarotropicQGField field) const noexcept {
    const auto& modes = descriptor_.modes();
    switch (field) {
        case WVBarotropicQGField::eta: return modes.etaFactor.data();
        case WVBarotropicQGField::pi:
        case WVBarotropicQGField::ssh: return modes.piFactor.data();
        case WVBarotropicQGField::psi: return modes.psiFactor.data();
        case WVBarotropicQGField::qgpv: return modes.qgpvFactor.data();
        case WVBarotropicQGField::zetaZ: return modes.zetaZFactor.data();
        case WVBarotropicQGField::u:
        case WVBarotropicQGField::v: return nullptr;
    }
    return nullptr;
}

WVKernelStatus WVTransformBarotropicQGKernel::fillHalfSpectrum(
    const WVComplexConstView& input, const WVComplex64* factors,
    std::size_t field, std::size_t fields) {
    if (factors == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Barotropic QG complex field factors are null."};
    const auto& mappings = descriptor_.halfSpectrumMappings();
    for (std::size_t index = 0; index < mappings.directRows.size(); ++index) {
        const auto iWV = mappings.directWVIndices[index];
        halfSpectrumScratch_[field + fields * mappings.directRows[index]] =
            multiply(input.data[iWV], factors[iWV]);
    }
    for (std::size_t index = 0; index < mappings.conjugatedRows.size();
         ++index) {
        const auto iWV = mappings.conjugatedWVIndices[index];
        halfSpectrumScratch_[field +
                             fields * mappings.conjugatedRows[index]] =
            conjugate(multiply(input.data[iWV], factors[iWV]));
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::fillHalfSpectrum(
    const WVComplexConstView& input, const double* factors,
    std::size_t field, std::size_t fields) {
    if (factors == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Barotropic QG real field factors are null."};
    const auto& mappings = descriptor_.halfSpectrumMappings();
    for (std::size_t index = 0; index < mappings.directRows.size(); ++index) {
        const auto iWV = mappings.directWVIndices[index];
        halfSpectrumScratch_[field + fields * mappings.directRows[index]] =
            multiply(input.data[iWV], factors[iWV]);
    }
    for (std::size_t index = 0; index < mappings.conjugatedRows.size();
         ++index) {
        const auto iWV = mappings.conjugatedWVIndices[index];
        halfSpectrumScratch_[field +
                             fields * mappings.conjugatedRows[index]] =
            conjugate(multiply(input.data[iWV], factors[iWV]));
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::forward(
    const WVRealConstView& input, WVComplexView& output, bool accumulate) {
    std::fill(halfSpectrumScratch_.begin(), halfSpectrumScratch_.end(),
              WVComplex64{});
    auto status = plans_[forwardPlan]->execute(
        input.data, halfSpectrumScratch_.data());
    if (!status) return status;
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const double normalization = 1.0 /
        static_cast<double>(descriptor_.spatialShape().elementCount());
    for (std::size_t index = 0; index < mappings.directRows.size(); ++index) {
        const auto iWV = mappings.directWVIndices[index];
        const auto value = multiply(
            halfSpectrumScratch_[mappings.directRows[index]], normalization);
        output.data[iWV] = accumulate
            ? WVComplex64{output.data[iWV].real + value.real,
                          output.data[iWV].imag + value.imag}
            : value;
    }
    for (std::size_t index = 0; index < mappings.conjugatedRows.size();
         ++index) {
        const auto iWV = mappings.conjugatedWVIndices[index];
        const auto value = multiply(
            conjugate(halfSpectrumScratch_[mappings.conjugatedRows[index]]),
            normalization);
        output.data[iWV] = accumulate
            ? WVComplex64{output.data[iWV].real + value.real,
                          output.data[iWV].imag + value.imag}
            : value;
    }
    ++metrics_.executionCount;
    ++metrics_.forwardExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::inverse(
    const WVComplexConstView& input, const WVComplex64* factors,
    WVRealView& output) {
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf *
                          descriptor_.configuration().Ny;
    std::fill_n(halfSpectrumScratch_.data(), halfRows, WVComplex64{});
    auto status = fillHalfSpectrum(input, factors, 0, 1);
    if (!status) return status;
    completeHermitianBoundaries(halfSpectrumScratch_.data(),
                                descriptor_.halfSpectrumMappings(), 1);
    status = plans_[inversePlan]->execute(halfSpectrumScratch_.data(),
                                          output.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::inverse(
    const WVComplexConstView& input, const double* factors,
    WVRealView& output) {
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf *
                          descriptor_.configuration().Ny;
    std::fill_n(halfSpectrumScratch_.data(), halfRows, WVComplex64{});
    auto status = fillHalfSpectrum(input, factors, 0, 1);
    if (!status) return status;
    completeHermitianBoundaries(halfSpectrumScratch_.data(),
                                descriptor_.halfSpectrumMappings(), 1);
    status = plans_[inversePlan]->execute(halfSpectrumScratch_.data(),
                                          output.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::transformQGPVToA0(
    const WVRealConstView& qgpv, WVComplexView& A0) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpatial(qgpv, descriptor_.spatialShape(), "QGPV");
    if (!status) return status;
    status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    const auto spatialBytes = descriptor_.spatialShape().elementCount() *
                              sizeof(double);
    const auto spectralBytes = descriptor_.spectralShape().elementCount() *
                               sizeof(WVComplex64);
    if (overlaps(qgpv.data, spatialBytes, A0.data, spectralBytes))
        return {WVKernelStatusCode::overlappingArrays,
                "QGPV input must not overlap compact A0 output."};
    ExecutionGuard guard(executing_);
    return forward(qgpv, A0);
}

WVKernelStatus WVTransformBarotropicQGKernel::transformA0ToQGPV(
    const WVComplexConstView& A0, WVRealView& qgpv) {
    return transformA0ToField(A0, WVBarotropicQGField::qgpv, qgpv);
}

WVKernelStatus WVTransformBarotropicQGKernel::transformA0ToField(
    const WVComplexConstView& A0, WVBarotropicQGField field,
    WVRealView& output) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    status = validateSpatial(output, descriptor_.spatialShape(), "Field");
    if (!status) return status;
    const auto spectralBytes = descriptor_.spectralShape().elementCount() *
                               sizeof(WVComplex64);
    const auto spatialBytes = descriptor_.spatialShape().elementCount() *
                              sizeof(double);
    if (overlaps(A0.data, spectralBytes, output.data, spatialBytes))
        return {WVKernelStatusCode::overlappingArrays,
                "A0 input must not overlap spatial field output."};
    ExecutionGuard guard(executing_);
    const auto* complex = complexFactors(field);
    status = complex != nullptr
        ? inverse(A0, complex, output)
        : inverse(A0, realFactors(field), output);
    if (status) ++metrics_.fieldEvaluationCount;
    return status;
}

WVKernelStatus
WVTransformBarotropicQGKernel::transformA0ToFieldWithDerivatives(
    const WVComplexConstView& A0, WVBarotropicQGField field,
    WVRealFieldBundleView& output) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    const auto spatial = descriptor_.spatialShape();
    if (output.shape.first != spatial.rows ||
        output.shape.second != spatial.columns ||
        output.shape.third != 1 || output.shape.fourth != 3)
        return {WVKernelStatusCode::invalidShape,
                "Barotropic field derivatives must have shape [Nx,Ny,1,3]."};
    if (output.data == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Barotropic field derivative output is null."};
    const auto spectralBytes = descriptor_.spectralShape().elementCount() *
                               sizeof(WVComplex64);
    const auto outputBytes = output.shape.elementCount() * sizeof(double);
    if (overlaps(A0.data, spectralBytes, output.data, outputBytes))
        return {WVKernelStatusCode::overlappingArrays,
                "A0 input must not overlap derivative output."};
    ExecutionGuard guard(executing_);
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf *
                          descriptor_.configuration().Ny;
    const auto R = spatial.elementCount();
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const auto& horizontal = descriptor_.fourierModes();
    const auto& modes = descriptor_.modes();
    for (std::size_t derivative = 0; derivative < 3; ++derivative) {
        std::fill_n(halfSpectrumScratch_.data(), halfRows, WVComplex64{});
        const auto valueAt = [&](std::size_t iWV) {
            auto factor = fieldFactor(modes, field, iWV);
            if (derivative == 1)
                factor = multiply(factor, {0.0, horizontal[iWV].k});
            else if (derivative == 2)
                factor = multiply(factor, {0.0, horizontal[iWV].l});
            return multiply(A0.data[iWV], factor);
        };
        for (std::size_t index = 0; index < mappings.directRows.size();
             ++index) {
            const auto iWV = mappings.directWVIndices[index];
            halfSpectrumScratch_[mappings.directRows[index]] = valueAt(iWV);
        }
        for (std::size_t index = 0; index < mappings.conjugatedRows.size();
             ++index) {
            const auto iWV = mappings.conjugatedWVIndices[index];
            halfSpectrumScratch_[mappings.conjugatedRows[index]] =
                conjugate(valueAt(iWV));
        }
        completeHermitianBoundaries(halfSpectrumScratch_.data(), mappings, 1);
        status = plans_[inversePlan]->execute(halfSpectrumScratch_.data(),
                                              output.data + derivative * R);
        if (!status) return status;
        ++metrics_.executionCount;
        ++metrics_.inverseExecutionCount;
    }
    ++metrics_.fieldEvaluationCount;
    ++metrics_.derivativeEvaluationCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::evolveA0(
    const WVComplexConstView& A0, double elapsedTime,
    WVComplexView& evolvedA0) const {
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    status = validateSpectral(evolvedA0, descriptor_.spectralShape(),
                              "Evolved A0");
    if (!status) return status;
    if (!std::isfinite(elapsedTime))
        return {WVKernelStatusCode::invalidConfiguration,
                "Barotropic linear-evolution time must be finite."};
    if (A0.data != evolvedA0.data) {
        std::copy_n(A0.data, descriptor_.spectralShape().elementCount(),
                    evolvedA0.data);
        metrics_.bytesCopied += descriptor_.spectralShape().elementCount() *
                                sizeof(WVComplex64);
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::inverseNonlinearFields(
    const WVComplexConstView& A0) {
    constexpr std::size_t fields = 4;
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf *
                          descriptor_.configuration().Ny;
    std::fill_n(halfSpectrumScratch_.data(), fields * halfRows,
                WVComplex64{});
    auto status = fillHalfSpectrum(A0, descriptor_.modes().uFactor.data(),
                                   0, fields);
    if (!status) return status;
    status = fillHalfSpectrum(A0, descriptor_.modes().vFactor.data(),
                              1, fields);
    if (!status) return status;
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const auto& horizontal = descriptor_.fourierModes();
    const auto& qgpv = descriptor_.modes().qgpvFactor;
    const auto derivativeValue = [&](std::size_t iWV, bool x) {
        return multiply(A0.data[iWV],
                        {0.0, qgpv[iWV] *
                                  (x ? horizontal[iWV].k
                                     : horizontal[iWV].l)});
    };
    for (std::size_t index = 0; index < mappings.directRows.size(); ++index) {
        const auto iWV = mappings.directWVIndices[index];
        const auto row = mappings.directRows[index];
        halfSpectrumScratch_[2 + fields * row] = derivativeValue(iWV, true);
        halfSpectrumScratch_[3 + fields * row] = derivativeValue(iWV, false);
    }
    for (std::size_t index = 0; index < mappings.conjugatedRows.size();
         ++index) {
        const auto iWV = mappings.conjugatedWVIndices[index];
        const auto row = mappings.conjugatedRows[index];
        halfSpectrumScratch_[2 + fields * row] =
            conjugate(derivativeValue(iWV, true));
        halfSpectrumScratch_[3 + fields * row] =
            conjugate(derivativeValue(iWV, false));
    }
    completeHermitianBoundaries(halfSpectrumScratch_.data(), mappings,
                                fields);
    status = plans_[nonlinearInversePlan]->execute(
        halfSpectrumScratch_.data(), realScratch_.data());
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::validateForcingOperation(
    const WVComplexConstView& A0, const WVComplexView& F0) const {
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    status = validateSpectral(F0, descriptor_.spectralShape(), "F0");
    if (!status) return status;
    const auto bytes = descriptor_.spectralShape().elementCount() *
                       sizeof(WVComplex64);
    if (overlaps(A0.data, bytes, F0.data, bytes))
        return {WVKernelStatusCode::overlappingArrays,
                "A0 state and F0 tendency must not overlap."};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::inverseQGPVDerivative(
    const WVComplexConstView& A0, bool xDerivative, WVRealView& output) {
    const auto halfRows = descriptor_.halfSpectrumMappings().NxHalf *
                          descriptor_.configuration().Ny;
    std::fill_n(halfSpectrumScratch_.data(), halfRows, WVComplex64{});
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const auto& horizontal = descriptor_.fourierModes();
    const auto& qgpv = descriptor_.modes().qgpvFactor;
    const auto valueAt = [&](std::size_t iWV) {
        return multiply(A0.data[iWV],
                        {0.0, qgpv[iWV] *
                                  (xDerivative ? horizontal[iWV].k
                                               : horizontal[iWV].l)});
    };
    for (std::size_t index = 0; index < mappings.directRows.size(); ++index) {
        const auto iWV = mappings.directWVIndices[index];
        halfSpectrumScratch_[mappings.directRows[index]] = valueAt(iWV);
    }
    for (std::size_t index = 0; index < mappings.conjugatedRows.size();
         ++index) {
        const auto iWV = mappings.conjugatedWVIndices[index];
        halfSpectrumScratch_[mappings.conjugatedRows[index]] =
            conjugate(valueAt(iWV));
    }
    completeHermitianBoundaries(halfSpectrumScratch_.data(), mappings, 1);
    const auto status = plans_[inversePlan]->execute(
        halfSpectrumScratch_.data(), output.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::ensureForcingFields(
    const WVComplexConstView& A0, bool requireQGPVDerivatives,
    WVBarotropicQGOperationWorkspace& workspace) {
    if (!workspace.physicalFieldsPrepared) {
        auto status = inverseNonlinearFields(A0);
        if (!status) return status;
        workspace.physicalFieldsPrepared = true;
        workspace.qgpvDerivativesPrepared = true;
        ++workspace.physicalFieldReconstructionCount;
        ++metrics_.forcingFieldReconstructionCount;
    } else {
        ++workspace.physicalFieldReuseCount;
        ++metrics_.forcingFieldReuseCount;
    }
    if (requireQGPVDerivatives && !workspace.qgpvDerivativesPrepared) {
        const auto spatial = descriptor_.spatialShape();
        const auto R = spatial.elementCount();
        WVRealView qx{realScratch_.data() + 2 * R, spatial};
        WVRealView qy{realScratch_.data() + 3 * R, spatial};
        auto status = inverseQGPVDerivative(A0, true, qx);
        if (!status) return status;
        status = inverseQGPVDerivative(A0, false, qy);
        if (!status) return status;
        workspace.qgpvDerivativesPrepared = true;
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::spatialDerivative(
    const WVRealConstView& input, bool xDerivative, WVRealView& output) {
    const auto& configuration = descriptor_.configuration();
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const auto halfRows = mappings.NxHalf * configuration.Ny;
    std::fill_n(halfSpectrumScratch_.data(), halfRows, WVComplex64{});
    auto status = plans_[forwardPlan]->execute(
        input.data, halfSpectrumScratch_.data());
    if (!status) return status;
    const double normalization = 1.0 /
        static_cast<double>(descriptor_.spatialShape().elementCount());
    for (std::size_t row = 0; row < halfRows; ++row) {
        const auto iK = row % mappings.NxHalf;
        const auto iL = row / mappings.NxHalf;
        const auto kMode = configuration.Nx % 2 == 0 &&
                                   iK == configuration.Nx / 2
                               ? -static_cast<std::int64_t>(
                                     configuration.Nx / 2)
                               : static_cast<std::int64_t>(iK);
        const auto lMode = iL < (configuration.Ny + 1) / 2
                               ? static_cast<std::int64_t>(iL)
                               : static_cast<std::int64_t>(iL) -
                                     static_cast<std::int64_t>(
                                         configuration.Ny);
        const bool derivativeNyquist =
            xDerivative
                ? configuration.Nx % 2 == 0 && iK == configuration.Nx / 2
                : configuration.Ny % 2 == 0 && iL == configuration.Ny / 2;
        const double wavenumber =
            derivativeNyquist
                ? 0.0
                : 2.0 * pi *
                      static_cast<double>(xDerivative ? kMode : lMode) /
                      (xDerivative ? configuration.Lx : configuration.Ly);
        halfSpectrumScratch_[row] = multiply(
            halfSpectrumScratch_[row], {0.0, normalization * wavenumber});
    }
    completeHermitianBoundaries(halfSpectrumScratch_.data(), mappings, 1);
    status = plans_[inversePlan]->execute(halfSpectrumScratch_.data(),
                                          output.data);
    if (!status) return status;
    ++metrics_.executionCount;
    ++metrics_.forwardExecutionCount;
    ++metrics_.executionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::antialiasScalarInPlace(
    WVRealView& scalar) {
    const auto& configuration = descriptor_.configuration();
    const auto& mappings = descriptor_.halfSpectrumMappings();
    const auto halfRows = mappings.NxHalf * configuration.Ny;
    auto status = plans_[forwardPlan]->execute(
        scalar.data, halfSpectrumScratch_.data());
    if (!status) return status;
    const double normalization = 1.0 /
        static_cast<double>(descriptor_.spatialShape().elementCount());
    const double maximumK = 2.0 * pi *
        static_cast<double>(configuration.Nx / 2) / configuration.Lx;
    const double cutoff = 2.0 * maximumK / 3.0;
    for (std::size_t row = 0; row < halfRows; ++row) {
        const auto iK = row % mappings.NxHalf;
        const auto iL = row / mappings.NxHalf;
        const auto lMode = iL < (configuration.Ny + 1) / 2
            ? static_cast<std::int64_t>(iL)
            : static_cast<std::int64_t>(iL) -
                  static_cast<std::int64_t>(configuration.Ny);
        const double k = 2.0 * pi * static_cast<double>(iK) /
                         configuration.Lx;
        const double l = 2.0 * pi * static_cast<double>(lMode) /
                         configuration.Ly;
        const bool excludedNyquist =
            (configuration.Nx % 2 == 0 && iK == configuration.Nx / 2) ||
            (configuration.Ny % 2 == 0 && iL == configuration.Ny / 2);
        const double scale =
            excludedNyquist ||
                    (configuration.shouldAntialias &&
                     std::hypot(k, l) > cutoff)
                ? 0.0
                : normalization;
        halfSpectrumScratch_[row] = multiply(halfSpectrumScratch_[row], scale);
    }
    completeHermitianBoundaries(halfSpectrumScratch_.data(), mappings, 1);
    status = plans_[inversePlan]->execute(halfSpectrumScratch_.data(),
                                          scalar.data);
    if (!status) return status;
    metrics_.executionCount += 2;
    ++metrics_.forwardExecutionCount;
    ++metrics_.inverseExecutionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::addPotentialVorticityAdvection(
    const WVComplexConstView& A0, WVComplexView& F0, bool accumulate,
    WVBarotropicQGOperationWorkspace& workspace) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateForcingOperation(A0, F0);
    if (!status) return status;
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, true, workspace);
    if (!status) return status;
    const auto R = descriptor_.spatialShape().elementCount();
    const double* u = realScratch_.data();
    const double* v = u + R;
    const double* qx = v + R;
    const double* qy = qx + R;
    double* tendency = realScratch_.data() + 4 * R;
    for (std::size_t index = 0; index < R; ++index)
        tendency[index] = -(u[index] * qx[index] +
                            v[index] * qy[index]);
    status = forward({tendency, descriptor_.spatialShape()}, F0, accumulate);
    if (!status) return status;
    ++workspace.spatialTendencyProjectionCount;
    ++metrics_.forcingSpatialProjectionCount;
    ++metrics_.nonlinearFluxCallCount;
    if (descriptor_.configuration().shouldAntialias)
        ++metrics_.antialiasedNonlinearFluxCallCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::addAdaptiveDamping(
    const WVComplexConstView& A0,
    const std::vector<double>& dampingOperator, WVComplexView& F0,
    bool accumulate, WVBarotropicQGOperationWorkspace& workspace) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateForcingOperation(A0, F0);
    if (!status) return status;
    if (dampingOperator.size() != descriptor_.Nkl())
        return {WVKernelStatusCode::invalidShape,
                "The Barotropic QG damping operator has the wrong length."};
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, false, workspace);
    if (!status) return status;
    const auto R = descriptor_.spatialShape().elementCount();
    const double* u = realScratch_.data();
    const double* v = u + R;
    double maximumSpeed = 0.0;
    for (std::size_t index = 0; index < R; ++index)
        maximumSpeed = std::max(
            maximumSpeed,
            std::sqrt(u[index] * u[index] + v[index] * v[index]));
    for (std::size_t index = 0; index < descriptor_.Nkl(); ++index) {
        const auto value = multiply(
            A0.data[index], maximumSpeed * dampingOperator[index]);
        F0.data[index] = accumulate
            ? WVComplex64{F0.data[index].real + value.real,
                          F0.data[index].imag + value.imag}
            : value;
    }
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::addLinearBottomFriction(
    const WVComplexConstView& A0, double rate, WVComplexView& F0,
    bool accumulate, WVBarotropicQGOperationWorkspace& workspace) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateForcingOperation(A0, F0);
    if (!status) return status;
    if (!std::isfinite(rate) || rate < 0.0)
        return {WVKernelStatusCode::invalidConfiguration,
                "Barotropic linear bottom friction must be finite and nonnegative."};
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, false, workspace);
    if (!status) return status;
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    WVRealView zeta{realScratch_.data() + 2 * R, spatial};
    status = inverse(A0, descriptor_.modes().zetaZFactor.data(), zeta);
    if (!status) return status;
    workspace.qgpvDerivativesPrepared = false;
    double* tendency = realScratch_.data() + 4 * R;
    for (std::size_t index = 0; index < R; ++index)
        tendency[index] = -rate * zeta.data[index];
    status = forward({tendency, spatial}, F0, accumulate);
    if (!status) return status;
    ++workspace.spatialTendencyProjectionCount;
    ++metrics_.forcingSpatialProjectionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::addQuadraticBottomFriction(
    const WVComplexConstView& A0, double drag, WVComplexView& F0,
    bool accumulate, WVBarotropicQGOperationWorkspace& workspace) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateForcingOperation(A0, F0);
    if (!status) return status;
    if (!std::isfinite(drag) || drag < 0.0)
        return {WVKernelStatusCode::invalidConfiguration,
                "Barotropic quadratic bottom friction must be finite and nonnegative."};
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, false, workspace);
    if (!status) return status;
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    const double* u = realScratch_.data();
    const double* v = u + R;
    double* speedV = realScratch_.data() + 2 * R;
    double* speedU = realScratch_.data() + 3 * R;
    for (std::size_t index = 0; index < R; ++index) {
        const double speed =
            std::sqrt(u[index] * u[index] + v[index] * v[index]);
        speedV[index] = speed * v[index];
        speedU[index] = speed * u[index];
    }
    WVRealView dxSpeedV{speedV, spatial};
    WVRealView dySpeedU{speedU, spatial};
    status = spatialDerivative({speedV, spatial}, true, dxSpeedV);
    if (!status) return status;
    status = spatialDerivative({speedU, spatial}, false, dySpeedU);
    if (!status) return status;
    workspace.qgpvDerivativesPrepared = false;
    double* tendency = realScratch_.data() + 4 * R;
    for (std::size_t index = 0; index < R; ++index)
        tendency[index] = -drag * (speedV[index] - speedU[index]);
    status = forward({tendency, spatial}, F0, accumulate);
    if (!status) return status;
    ++workspace.spatialTendencyProjectionCount;
    ++metrics_.forcingSpatialProjectionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::addBetaPlanePVAdvection(
    const WVComplexConstView& A0, double beta, WVComplexView& F0,
    bool accumulate, WVBarotropicQGOperationWorkspace& workspace) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateForcingOperation(A0, F0);
    if (!status) return status;
    if (!std::isfinite(beta))
        return {WVKernelStatusCode::invalidConfiguration,
                "The beta-plane gradient must be finite."};
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, false, workspace);
    if (!status) return status;
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    const double* v = realScratch_.data() + R;
    double* tendency = realScratch_.data() + 4 * R;
    for (std::size_t index = 0; index < R; ++index)
        tendency[index] = -beta * v[index];
    status = forward({tendency, spatial}, F0, accumulate);
    if (!status) return status;
    ++workspace.spatialTendencyProjectionCount;
    ++metrics_.forcingSpatialProjectionCount;
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::nonlinearFlux(
    const WVComplexConstView& A0, WVComplexView& F0) {
    WVBarotropicQGOperationWorkspace workspace;
    return addPotentialVorticityAdvection(A0, F0, false, workspace);
}

WVKernelStatus WVTransformBarotropicQGKernel::totalEnergy(
    const WVComplexConstView& A0, double& energy) const {
    const auto status = validateSpectral(A0, descriptor_.spectralShape(),
                                         "A0");
    if (!status) return status;
    energy = 0.0;
    const auto& factors = descriptor_.modes().energyFactor;
    for (std::size_t index = 0; index < factors.size(); ++index)
        energy += factors[index] *
                  (A0.data[index].real * A0.data[index].real +
                   A0.data[index].imag * A0.data[index].imag);
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::totalEnstrophy(
    const WVComplexConstView& A0, double& enstrophy) const {
    const auto status = validateSpectral(A0, descriptor_.spectralShape(),
                                         "A0");
    if (!status) return status;
    enstrophy = 0.0;
    const auto& factors = descriptor_.modes().enstrophyFactor;
    for (std::size_t index = 0; index < factors.size(); ++index)
        enstrophy += factors[index] *
                     (A0.data[index].real * A0.data[index].real +
                      A0.data[index].imag * A0.data[index].imag);
    return WVKernelStatus::ok();
}

WVKernelStatus
WVTransformBarotropicQGKernel::totalEnergySpatiallyIntegrated(
    const WVComplexConstView& A0, double& energy) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    ExecutionGuard guard(executing_);
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    WVRealView u{realScratch_.data(), spatial};
    WVRealView v{realScratch_.data() + R, spatial};
    WVRealView eta{realScratch_.data() + 2 * R, spatial};
    status = inverse(A0, descriptor_.modes().uFactor.data(), u);
    if (!status) return status;
    status = inverse(A0, descriptor_.modes().vFactor.data(), v);
    if (!status) return status;
    status = inverse(A0, descriptor_.modes().etaFactor.data(), eta);
    if (!status) return status;
    double sum = 0.0;
    for (std::size_t index = 0; index < R; ++index)
        sum += descriptor_.configuration().h *
                   (u.data[index] * u.data[index] +
                    v.data[index] * v.data[index]) /
                   2.0 +
               descriptor_.configuration().g *
                   eta.data[index] * eta.data[index] / 2.0;
    energy = sum / static_cast<double>(R);
    return WVKernelStatus::ok();
}

WVKernelStatus
WVTransformBarotropicQGKernel::totalEnstrophySpatiallyIntegrated(
    const WVComplexConstView& A0, double& enstrophy) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    ExecutionGuard guard(executing_);
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    WVRealView qgpv{realScratch_.data(), spatial};
    status = inverse(A0, descriptor_.modes().qgpvFactor.data(), qgpv);
    if (!status) return status;
    double sum = 0.0;
    for (std::size_t index = 0; index < R; ++index)
        sum += qgpv.data[index] * qgpv.data[index];
    enstrophy = descriptor_.configuration().h * sum /
                (2.0 * static_cast<double>(R));
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::uvMax(
    const WVComplexConstView& A0, double& maximumSpeed) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    ExecutionGuard guard(executing_);
    const auto spatial = descriptor_.spatialShape();
    const auto R = spatial.elementCount();
    WVRealView u{realScratch_.data(), spatial};
    WVRealView v{realScratch_.data() + R, spatial};
    status = inverse(A0, descriptor_.modes().uFactor.data(), u);
    if (!status) return status;
    status = inverse(A0, descriptor_.modes().vFactor.data(), v);
    if (!status) return status;
    maximumSpeed = 0.0;
    for (std::size_t index = 0; index < R; ++index)
        maximumSpeed = std::max(
            maximumSpeed,
            std::sqrt(u.data[index] * u.data[index] +
                      v.data[index] * v.data[index]));
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::advectScalar(
    const WVComplexConstView& A0, const WVRealConstView& scalar,
    bool shouldAntialias, WVRealView& rightHandSide) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    status = validateSpatial(scalar, descriptor_.spatialShape(), "scalar");
    if (!status) return status;
    status = validateSpatial(rightHandSide, descriptor_.spatialShape(),
                             "scalar RHS");
    if (!status) return status;
    const auto R = descriptor_.spatialShape().elementCount();
    if (overlaps(scalar.data, R * sizeof(double), rightHandSide.data,
                 R * sizeof(double)))
        return {WVKernelStatusCode::overlappingArrays,
                "Barotropic QG tracer state and RHS must not overlap."};
    ExecutionGuard guard(executing_);
    WVRealView u{realScratch_.data(), descriptor_.spatialShape()};
    WVRealView v{realScratch_.data() + R, descriptor_.spatialShape()};
    WVRealView dx{realScratch_.data() + 2 * R,
                  descriptor_.spatialShape()};
    WVRealView dy{realScratch_.data() + 3 * R,
                  descriptor_.spatialShape()};
    status = inverse(A0, descriptor_.modes().uFactor.data(), u);
    if (!status) return status;
    status = inverse(A0, descriptor_.modes().vFactor.data(), v);
    if (!status) return status;
    status = spatialDerivative(scalar, true, dx);
    if (!status) return status;
    status = spatialDerivative(scalar, false, dy);
    if (!status) return status;
    for (std::size_t index = 0; index < R; ++index)
        rightHandSide.data[index] =
            -(u.data[index] * dx.data[index] +
              v.data[index] * dy.data[index]);
    if (shouldAntialias)
        return antialiasScalarInPlace(rightHandSide);
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::prepareAdvectionFields(
    const WVComplexConstView& A0,
    WVBarotropicQGOperationWorkspace& workspace,
    WVRealFieldBundleConstView& fields) {
    fields = {};
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpectral(A0, descriptor_.spectralShape(), "A0");
    if (!status) return status;
    ExecutionGuard guard(executing_);
    status = ensureForcingFields(A0, false, workspace);
    if (!status) return status;
    const auto& configuration = descriptor_.configuration();
    fields = {realScratch_.data(),
              {configuration.Nx, configuration.Ny, 1, 2}};
    return WVKernelStatus::ok();
}

WVKernelStatus WVTransformBarotropicQGKernel::advectScalarWithAdvectionFields(
    const WVRealConstView& scalar,
    const WVRealFieldBundleConstView& advectionFields,
    bool shouldAntialias, WVRealView& rightHandSide) {
    if (executing_)
        return {WVKernelStatusCode::reentrantExecution,
                "The Barotropic QG kernel is not reentrant."};
    auto status = validateSpatial(scalar, descriptor_.spatialShape(), "scalar");
    if (!status) return status;
    status = validateSpatial(rightHandSide, descriptor_.spatialShape(),
                             "scalar RHS");
    if (!status) return status;
    const auto& configuration = descriptor_.configuration();
    if (advectionFields.data == nullptr ||
        advectionFields.shape.first != configuration.Nx ||
        advectionFields.shape.second != configuration.Ny ||
        advectionFields.shape.third != 1 ||
        advectionFields.shape.fourth != 2)
        return {WVKernelStatusCode::invalidShape,
                "Barotropic QG advection fields must have shape [Nx,Ny,1,2]."};
    const auto R = descriptor_.spatialShape().elementCount();
    if (overlaps(scalar.data, R * sizeof(double), rightHandSide.data,
                 R * sizeof(double)))
        return {WVKernelStatusCode::overlappingArrays,
                "Barotropic QG tracer state and RHS must not overlap."};
    ExecutionGuard guard(executing_);
    WVRealView dx{realScratch_.data() + 2 * R,
                  descriptor_.spatialShape()};
    WVRealView dy{realScratch_.data() + 3 * R,
                  descriptor_.spatialShape()};
    status = spatialDerivative(scalar, true, dx);
    if (!status) return status;
    status = spatialDerivative(scalar, false, dy);
    if (!status) return status;
    const double* u = advectionFields.data;
    const double* v = u + R;
    for (std::size_t index = 0; index < R; ++index)
        rightHandSide.data[index] =
            -(u[index] * dx.data[index] + v[index] * dy.data[index]);
    if (shouldAntialias)
        return antialiasScalarInPlace(rightHandSide);
    return WVKernelStatus::ok();
}

std::size_t WVTransformBarotropicQGKernel::enforceReality(
    WVComplexView& A0) const noexcept {
    if (A0.data == nullptr ||
        A0.shape.rows != descriptor_.spectralShape().rows ||
        A0.shape.columns != descriptor_.spectralShape().columns)
        return 0;
    std::size_t modified = 0;
    const auto& modes = descriptor_.fourierModes();
    for (std::size_t index = 0; index < modes.size(); ++index) {
        if (descriptor_.modes().qgpvFactor[index] == 0.0 &&
            (A0.data[index].real != 0.0 || A0.data[index].imag != 0.0)) {
            A0.data[index] = {};
            ++modified;
            continue;
        }
        if (modes[index].dftPrimaryIndex == modes[index].dftConjugateIndex &&
            A0.data[index].imag != 0.0) {
            A0.data[index].imag = 0.0;
            ++modified;
        }
    }
    return modified;
}

} // namespace wavevortex
