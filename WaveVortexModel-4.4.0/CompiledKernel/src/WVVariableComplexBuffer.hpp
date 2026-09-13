#pragma once

#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include <vector>

namespace wavevortex::spectral_detail {

// One prepared complex representation with stable storage for the lifetime of
// a variable-stratification kernel.  Callers keep offsets, rather than derived
// pointers, so setup-time allocation cannot leave stale field views.
class WVVariableComplexBuffer final {
public:
    WVVariableComplexBuffer(std::size_t elements,WVComplexRepresentation representation)
        : representation_(representation) {
        if (representation_==WVComplexRepresentation::interleaved) interleaved_.resize(elements);
        else { real_.resize(elements); imag_.resize(elements); }
    }

    std::size_t size() const noexcept {
        return representation_==WVComplexRepresentation::interleaved ? interleaved_.size() : real_.size();
    }
    WVComplexRepresentation representation() const noexcept { return representation_; }
    WVComplexInput input(std::size_t offset,std::size_t count) const noexcept {
        if (representation_==WVComplexRepresentation::interleaved)
            return {interleaved_.data()+offset,nullptr,nullptr,count*sizeof(WVComplex64)};
        return {nullptr,real_.data()+offset,imag_.data()+offset,count*sizeof(double)};
    }
    WVComplexOutput output(std::size_t offset,std::size_t count) noexcept {
        if (representation_==WVComplexRepresentation::interleaved)
            return {interleaved_.data()+offset,nullptr,nullptr,count*sizeof(WVComplex64)};
        return {nullptr,real_.data()+offset,imag_.data()+offset,count*sizeof(double)};
    }
    std::size_t capacityBytes() const noexcept {
        return interleaved_.capacity()*sizeof(WVComplex64)+
            (real_.capacity()+imag_.capacity())*sizeof(double);
    }

private:
    WVComplexRepresentation representation_;
    std::vector<WVComplex64> interleaved_;
    std::vector<double> real_,imag_;
};

} // namespace wavevortex::spectral_detail
