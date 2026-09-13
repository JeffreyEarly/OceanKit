#pragma once

#include <cstddef>

#if defined(_MSC_VER)
#define WV_NATIVE_RESTRICT __restrict
#elif defined(__clang__) || defined(__GNUC__)
#define WV_NATIVE_RESTRICT __restrict__
#else
#define WV_NATIVE_RESTRICT
#endif

namespace wavevortex::native_detail {

// Private provider contract: every pointer span is disjoint. Flux and derivative
// are separate, unexposed scratch halves; velocity and eta are separate field
// slices, and public preflight keeps fields disjoint from spectral/density input.
inline void accumulateAdvection(
    std::size_t count,double* WV_NATIVE_RESTRICT flux,
    const double* WV_NATIVE_RESTRICT velocity,
    const double* WV_NATIVE_RESTRICT derivative) noexcept {
    for (std::size_t i=0;i<count;++i)
        flux[i]-=velocity[i]*(derivative[i]+0.0);
}

inline void accumulateDensityAdvection(
    std::size_t count,double* WV_NATIVE_RESTRICT flux,
    const double* WV_NATIVE_RESTRICT velocity,
    const double* WV_NATIVE_RESTRICT derivative,
    const double* WV_NATIVE_RESTRICT eta,double densityCorrection) noexcept {
    for (std::size_t i=0;i<count;++i) {
        const double correction=eta[i]*densityCorrection;
        flux[i]-=velocity[i]*(derivative[i]+correction);
    }
}

} // namespace wavevortex::native_detail

#undef WV_NATIVE_RESTRICT
