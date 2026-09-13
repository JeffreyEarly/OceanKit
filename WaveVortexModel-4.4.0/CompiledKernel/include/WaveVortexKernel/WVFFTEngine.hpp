#pragma once

#include "WVKernelTypes.hpp"

#include <memory>
#include <string>
#include <vector>

namespace wavevortex {

enum class WVFFTPlanKind : std::uint32_t {
    horizontalRealToComplex2D,
    horizontalComplexToReal2D,
    verticalDCTI,
    verticalDSTI
};

struct WVFFTDimension {
    std::size_t count = 0;
    std::ptrdiff_t inputStride = 0;
    std::ptrdiff_t outputStride = 0;
};

struct WVFFTPlanSpecification {
    WVFFTPlanKind kind = WVFFTPlanKind::horizontalRealToComplex2D;
    std::vector<WVFFTDimension> transformDimensions;
    std::vector<WVFFTDimension> batchDimensions;
    std::size_t inputBytes = 0;
    std::size_t outputBytes = 0;
    // Out-of-place execution must preserve input unless explicitly allowed here.
    // In-place execution inherently overwrites its shared input/output storage.
    // Providers must reject unsupported preservation/placement at plan creation.
    bool destroysInput = false;
    bool inPlace = false;
};

class WVFFTPlan {
public:
    virtual ~WVFFTPlan() = default;
    // Execution must use the placement selected at setup: identical addresses
    // for in-place plans, disjoint addressed spans for out-of-place plans.
    virtual WVKernelStatus execute(const void* input, void* output) = 0;
    virtual std::size_t persistentBytes() const noexcept = 0;
};

struct WVRetainedHorizontalSpecification;
struct WVRealInput;
struct WVRealOutput;
struct WVComplexInput;
struct WVComplexOutput;
struct WVRetainedAdvectionWork;
struct WVRetainedAdvectionCounts;

// Synchronous consumer of completed contiguous physical-grid ranges. Different
// ranges may be delivered concurrently. The callback must not throw, reenter
// the transform, mutate its inputs, or publish an evaluation result. Its context
// lives until inverseAndConsume returns; publish only after that call succeeds.
struct WVRealOutputConsumer {
    void* context = nullptr;
    void (*consume)(void*,std::size_t begin,std::size_t end,const double*) noexcept = nullptr;
};

// Optional prepared retained transform. The enclosing operator validates all
// buffers, Hermitian constraints and exclusive workspace use before execution.
class WVRetainedHorizontalPlan {
public:
    virtual ~WVRetainedHorizontalPlan() = default;
    virtual WVKernelStatus forward(WVRealInput, WVComplexOutput) = 0;
    virtual WVKernelStatus inverse(WVComplexInput, WVRealOutput) = 0;
    virtual bool supportsInverseConsumer() const noexcept { return false; }
    virtual WVKernelStatus inverseAndConsume(WVComplexInput, WVRealOutput,
        const WVRealOutputConsumer&);
    virtual bool supportsAdvection(std::size_t) const noexcept { return false; }
    virtual WVKernelStatus prepareAdvection(std::size_t);
    virtual WVKernelStatus advection(const WVRetainedAdvectionWork&,
        WVRetainedAdvectionCounts&);
    virtual std::size_t persistentBytes() const noexcept = 0;
    virtual std::size_t planBytesLowerBound() const noexcept = 0;
    // persistentBytes includes shared dependencies. Owners of multiple plans
    // may subtract repeated sharedResourceBytes with the same nonnull identity.
    virtual const void* sharedResourceIdentity() const noexcept { return nullptr; }
    virtual std::size_t sharedResourceBytes() const noexcept { return 0; }
    virtual std::size_t workerCount() const noexcept = 0;
    virtual const char* identifier() const noexcept = 0;
};

class WVFFTEngine {
public:
    virtual ~WVFFTEngine() = default;
    virtual std::string identifier() const = 0;
    virtual std::string libraryIdentity() const { return {}; }
    // Provider-owned C++ storage retained by the engine itself. Individual
    // plans report their storage separately through WVFFTPlan.
    virtual std::size_t persistentBytes() const noexcept = 0;
    virtual WVKernelStatus createRetainedHorizontalPlan(const WVRetainedHorizontalSpecification&,
        std::unique_ptr<WVRetainedHorizontalPlan>&) {
        return {WVKernelStatusCode::unsupportedOperation,"Provider has no retained horizontal schedule."};
    }
    virtual WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) = 0;
};

} // namespace wavevortex
