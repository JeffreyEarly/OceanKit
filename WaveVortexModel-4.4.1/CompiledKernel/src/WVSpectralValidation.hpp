#pragma once
#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace wavevortex::spectral_detail {
inline std::size_t product(std::size_t a, std::size_t b) {
    constexpr auto maximum = static_cast<std::size_t>(PTRDIFF_MAX);
    if (a && b > maximum/a) throw std::overflow_error("Spectral extent overflow.");
    return a*b;
}
inline std::size_t span(std::size_t rows, std::size_t columns, std::size_t rs, std::size_t cs, std::size_t size) {
    if (!rows || !columns || !rs || !cs) throw std::invalid_argument("Spectral extents and strides must be positive.");
    const auto a = product(rows-1,rs), b = product(columns-1,cs);
    if (b >= static_cast<std::size_t>(PTRDIFF_MAX)-a) throw std::overflow_error("Spectral span overflow.");
    const auto divisor = std::gcd(rs,cs);
    if (rows > cs/divisor && columns > rs/divisor) throw std::invalid_argument("Spectral layout aliases its own elements.");
    return product(a+b+1,size);
}
inline std::size_t validate(const WVComplexLayout& l) {
    if (l.family.empty() || l.modeSet.empty()) throw std::invalid_argument("Field family and ordered mode-set identity are required.");
    if (l.representation != WVComplexRepresentation::split && l.representation != WVComplexRepresentation::interleaved)
        throw std::invalid_argument("Unknown complex representation.");
    return span(l.rows,l.columns,l.rowStride,l.columnStride,l.representation == WVComplexRepresentation::split ? sizeof(double) : sizeof(WVComplex64));
}
inline bool overlap(const void* a, std::size_t as, const void* b, std::size_t bs) noexcept {
    const auto x = reinterpret_cast<std::uintptr_t>(a), y = reinterpret_cast<std::uintptr_t>(b);
    return x <= y ? y-x < as : x-y < bs;
}
inline bool addressFits(const void* p, std::size_t bytes, std::size_t alignment) noexcept {
    const auto address = reinterpret_cast<std::uintptr_t>(p);
    return p && address % alignment == 0 && bytes <= UINTPTR_MAX-address;
}
inline WVKernelStatus validateStorage(WVComplexRepresentation r, std::size_t bytes, WVComplexInput input) {
    if (input.bytes < bytes) return {WVKernelStatusCode::invalidShape,"Complex buffer capacity is too small."};
    if (r == WVComplexRepresentation::interleaved) {
        if (input.real || input.imag || !addressFits(input.interleaved,bytes,alignof(WVComplex64)))
            return {WVKernelStatusCode::invalidPointer,"Invalid interleaved storage."};
    } else {
        if (input.interleaved || !addressFits(input.real,bytes,alignof(double)) || !addressFits(input.imag,bytes,alignof(double)))
            return {WVKernelStatusCode::invalidPointer,"Invalid split storage."};
        if (overlap(input.real,bytes,input.imag,bytes)) return {WVKernelStatusCode::overlappingArrays,"Split components overlap."};
    }
    return WVKernelStatus::ok();
}
inline bool storageOverlap(WVComplexInput a, std::size_t as, WVComplexInput b, std::size_t bs) noexcept {
    const void* ap[] = {a.interleaved,a.real,a.imag}; const void* bp[] = {b.interleaved,b.real,b.imag};
    for (auto x : ap) for (auto y : bp) if (x && y && overlap(x,as,y,bs)) return true;
    return false;
}
inline bool realOverlap(const void* p, std::size_t bytes, WVComplexInput c, std::size_t cs) noexcept {
    const void* pointers[] = {c.interleaved,c.real,c.imag};
    for (const void* q : pointers) if (q && overlap(p,bytes,q,cs)) return true;
    return false;
}
inline WVComplex64 read(WVComplexInput b, std::size_t index) noexcept {
    return b.interleaved ? b.interleaved[index] : WVComplex64{b.real[index],b.imag[index]};
}
inline void write(WVComplexOutput b, std::size_t index, WVComplex64 value) noexcept {
    if (b.interleaved) b.interleaved[index] = value;
    else { b.real[index] = value.real; b.imag[index] = value.imag; }
}
inline std::size_t identityBytes(const WVComplexLayout& l) noexcept { return l.family.capacity()+l.modeSet.capacity(); }
struct ActiveCall {
    std::atomic<bool>& active;
    bool entered;
    explicit ActiveCall(std::atomic<bool>& flag) : active(flag), entered(!flag.exchange(true)) {}
    ~ActiveCall() { if (entered) active.store(false); }
};
} // namespace wavevortex::spectral_detail
