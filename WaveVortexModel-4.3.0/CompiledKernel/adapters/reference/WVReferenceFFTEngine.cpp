#include "WVReferenceFFTEngine.hpp"

#include <cmath>
#include <complex>
#include <functional>
#include <vector>

namespace wavevortex {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void forEachBatch(const WVFFTPlanSpecification& specification, const std::function<void(std::ptrdiff_t,std::ptrdiff_t)>& action) {
    std::function<void(std::size_t,std::ptrdiff_t,std::ptrdiff_t)> visit;
    visit = [&](std::size_t dimension, std::ptrdiff_t inputOffset, std::ptrdiff_t outputOffset) {
        if (dimension == specification.batchDimensions.size()) { action(inputOffset, outputOffset); return; }
        const auto& batch = specification.batchDimensions[dimension];
        for (std::size_t index = 0; index < batch.count; ++index) visit(dimension + 1, inputOffset + static_cast<std::ptrdiff_t>(index) * batch.inputStride, outputOffset + static_cast<std::ptrdiff_t>(index) * batch.outputStride);
    };
    visit(0, 0, 0);
}

class ReferencePlan final : public WVFFTPlan {
public:
    explicit ReferencePlan(WVFFTPlanSpecification specification) : specification_(std::move(specification)) {}

    WVKernelStatus execute(const void* input, void* output) override {
        if (!input || !output) return {WVKernelStatusCode::invalidPointer, "Reference FFT received a null pointer."};
        switch (specification_.kind) {
            case WVFFTPlanKind::horizontalRealToComplex2D: return r2c(static_cast<const double*>(input), static_cast<WVComplex64*>(output));
            case WVFFTPlanKind::horizontalComplexToReal2D: return c2r(static_cast<const WVComplex64*>(input), static_cast<double*>(output));
            case WVFFTPlanKind::verticalDCTI: return r2r(static_cast<const double*>(input), static_cast<double*>(output), false);
            case WVFFTPlanKind::verticalDSTI: return r2r(static_cast<const double*>(input), static_cast<double*>(output), true);
        }
        return {WVKernelStatusCode::unsupportedOperation, "Unknown reference FFT plan."};
    }

    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this) +
               specification_.transformDimensions.capacity() *
                   sizeof(WVFFTDimension) +
               specification_.batchDimensions.capacity() *
                   sizeof(WVFFTDimension);
    }

private:
    WVKernelStatus r2c(const double* input, WVComplex64* output) {
        if (specification_.transformDimensions.size() != 2) return {WVKernelStatusCode::invalidConfiguration, "Reference r2c requires rank two."};
        const auto& yDimension = specification_.transformDimensions[0];
        const auto& xDimension = specification_.transformDimensions[1];
        forEachBatch(specification_, [&](std::ptrdiff_t inputBase, std::ptrdiff_t outputBase) {
            for (std::size_t ky = 0; ky < yDimension.count; ++ky) for (std::size_t kx = 0; kx <= xDimension.count / 2; ++kx) {
                std::complex<double> sum;
                for (std::size_t y = 0; y < yDimension.count; ++y) for (std::size_t x = 0; x < xDimension.count; ++x) {
                    const double angle = -2.0*pi*(static_cast<double>(kx*x)/static_cast<double>(xDimension.count)+static_cast<double>(ky*y)/static_cast<double>(yDimension.count));
                    sum += input[inputBase+static_cast<std::ptrdiff_t>(y)*yDimension.inputStride+static_cast<std::ptrdiff_t>(x)*xDimension.inputStride]*std::complex<double>(std::cos(angle),std::sin(angle));
                }
                auto& value = output[outputBase+static_cast<std::ptrdiff_t>(ky)*yDimension.outputStride+static_cast<std::ptrdiff_t>(kx)*xDimension.outputStride];
                value = {sum.real(),sum.imag()};
            }
        });
        return WVKernelStatus::ok();
    }

    WVKernelStatus c2r(const WVComplex64* input, double* output) {
        if (specification_.transformDimensions.size() != 2) return {WVKernelStatusCode::invalidConfiguration, "Reference c2r requires rank two."};
        const auto& yDimension = specification_.transformDimensions[0];
        const auto& xDimension = specification_.transformDimensions[1];
        forEachBatch(specification_, [&](std::ptrdiff_t inputBase, std::ptrdiff_t outputBase) {
            for (std::size_t y = 0; y < yDimension.count; ++y) for (std::size_t x = 0; x < xDimension.count; ++x) {
                std::complex<double> sum;
                for (std::size_t ky = 0; ky < yDimension.count; ++ky) for (std::size_t kx = 0; kx < xDimension.count; ++kx) {
                    std::size_t storedKx = kx;
                    std::size_t storedKy = ky;
                    bool shouldConjugate = false;
                    if (kx > xDimension.count/2) { storedKx = xDimension.count-kx; storedKy = (yDimension.count-ky)%yDimension.count; shouldConjugate = true; }
                    const auto raw = input[inputBase+static_cast<std::ptrdiff_t>(storedKy)*yDimension.inputStride+static_cast<std::ptrdiff_t>(storedKx)*xDimension.inputStride];
                    const std::complex<double> coefficient(raw.real,shouldConjugate ? -raw.imag : raw.imag);
                    const double angle = 2.0*pi*(static_cast<double>(kx*x)/static_cast<double>(xDimension.count)+static_cast<double>(ky*y)/static_cast<double>(yDimension.count));
                    sum += coefficient*std::complex<double>(std::cos(angle),std::sin(angle));
                }
                output[outputBase+static_cast<std::ptrdiff_t>(y)*yDimension.outputStride+static_cast<std::ptrdiff_t>(x)*xDimension.outputStride] = sum.real();
            }
        });
        return WVKernelStatus::ok();
    }

    WVKernelStatus r2r(const double* input, double* output, bool sine) {
        if (specification_.transformDimensions.size() != 1) return {WVKernelStatusCode::invalidConfiguration, "Reference r2r requires rank one."};
        const auto& dimension = specification_.transformDimensions[0];
        forEachBatch(specification_, [&](std::ptrdiff_t inputBase, std::ptrdiff_t outputBase) {
            std::vector<double> source(dimension.count);
            for (std::size_t n = 0; n < dimension.count; ++n) source[n] = input[inputBase+static_cast<std::ptrdiff_t>(n)*dimension.inputStride];
            for (std::size_t k = 0; k < dimension.count; ++k) {
                double sum = 0.0;
                if (sine) {
                    for (std::size_t n = 0; n < dimension.count; ++n) sum += 2.0*source[n]*std::sin(pi*static_cast<double>((n+1)*(k+1))/static_cast<double>(dimension.count+1));
                } else {
                    sum = source.front()+(k%2 == 0 ? source.back() : -source.back());
                    for (std::size_t n = 1; n+1 < dimension.count; ++n) sum += 2.0*source[n]*std::cos(pi*static_cast<double>(n*k)/static_cast<double>(dimension.count-1));
                }
                output[outputBase+static_cast<std::ptrdiff_t>(k)*dimension.outputStride] = sum;
            }
        });
        return WVKernelStatus::ok();
    }

    WVFFTPlanSpecification specification_;
};

} // namespace

WVKernelStatus WVReferenceFFTEngine::createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) {
    plan = std::make_unique<ReferencePlan>(specification);
    return WVKernelStatus::ok();
}

} // namespace wavevortex
