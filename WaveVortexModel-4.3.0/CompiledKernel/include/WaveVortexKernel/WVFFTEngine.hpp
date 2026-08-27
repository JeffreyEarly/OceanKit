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
    bool destroysInput = false;
    bool inPlace = false;
};

class WVFFTPlan {
public:
    virtual ~WVFFTPlan() = default;
    virtual WVKernelStatus execute(const void* input, void* output) = 0;
    virtual std::size_t persistentBytes() const noexcept = 0;
};

class WVFFTEngine {
public:
    virtual ~WVFFTEngine() = default;
    virtual std::string identifier() const = 0;
    virtual std::string libraryIdentity() const { return {}; }
    // Provider-owned C++ storage retained by the engine itself. Individual
    // plans report their storage separately through WVFFTPlan.
    virtual std::size_t persistentBytes() const noexcept = 0;
    virtual WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) = 0;
};

} // namespace wavevortex
