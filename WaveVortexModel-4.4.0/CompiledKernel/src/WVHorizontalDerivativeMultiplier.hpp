#pragma once
#include "WaveVortexKernel/WVKernelTypes.hpp"
#include <cmath>
#include <cstdint>
namespace wavevortex::kernel_detail {
// Strides allow one arithmetic path for plane-major and z-fast spectra.
inline WVKernelStatus applyHorizontalDerivativeMultiplier(
    WVComplex64* values,std::size_t Nx,std::size_t Ny,std::size_t planes,
    std::size_t planeStride,std::size_t rowStride,double Lx,double Ly,
    bool xDerivative,unsigned order) {
    if (!order) return {WVKernelStatusCode::invalidConfiguration,"Horizontal derivative order must be positive."};
    const auto n=xDerivative ? Nx : Ny; const double length=xDerivative ? Lx : Ly;
    if (order>1) {
        const double maximum=2*std::acos(-1.0)*static_cast<double>(n/2)/length;
        if (!std::isfinite(std::pow(maximum,order))) return {WVKernelStatusCode::numericalFailure,"Horizontal derivative multiplier overflow."};
    }
    const auto half=Nx/2+1; const double scale=1.0/static_cast<double>(Nx*Ny);
    for (std::size_t y=0;y<Ny;++y) for (std::size_t x=0;x<half;++x) {
        const auto i=xDerivative ? x : y;
        const auto mode=i<=n/2 ? static_cast<std::int64_t>(i) : static_cast<std::int64_t>(i)-static_cast<std::int64_t>(n);
        for (std::size_t plane=0;plane<planes;++plane) {
            auto& value=values[plane*planeStride+(x+half*y)*rowStride];
            if (order==1) {
                const double k=n%2==0 && i==n/2 ? 0.0 : 2*std::acos(-1.0)*static_cast<double>(mode)/length*scale;
                value={-k*value.imag,k*value.real}; continue;
            }
            if (n%2==0 && i==n/2 && order%2) { value={}; continue; }
            const double wave=2*std::acos(-1.0)*static_cast<double>(mode)/length;
            const double magnitude=std::pow(wave,order)*scale; WVComplex64 factor;
            switch(order%4) { case 0: factor={magnitude,0}; break; case 1: factor={0,magnitude}; break; case 2: factor={-magnitude,0}; break; default: factor={0,-magnitude}; }
            value={value.real*factor.real-value.imag*factor.imag,value.real*factor.imag+value.imag*factor.real};
        }
    }
    return WVKernelStatus::ok();
}
} // namespace wavevortex::kernel_detail
