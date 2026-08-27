#pragma once

#include "WaveVortexKernel/WVFFTEngine.hpp"

namespace wavevortex {

// Deliberately simple direct transforms for portable correctness testing.
// This provider is not intended for production-size performance work.
class WVReferenceFFTEngine final : public WVFFTEngine {
public:
    std::string identifier() const override { return "reference-direct"; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) override;
};

} // namespace wavevortex
