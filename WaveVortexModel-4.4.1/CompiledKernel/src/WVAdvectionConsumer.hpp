#pragma once
#include <cstddef>

namespace wavevortex::kernel_detail {
// Each callback owns disjoint output cells. Directional callbacks complete in
// x/y/z order, preserving the original accumulation and density correction.
struct WVAdvectionConsumer {
    double* flux;
    const double* velocity;
    const double* eta;
    const double* dLnN2;
    std::size_t planeSize;
    bool densityCorrection;
    static void consume(void* context,std::size_t begin,std::size_t end,
        const double* derivative) noexcept {
        const auto& a=*static_cast<const WVAdvectionConsumer*>(context);
        for (std::size_t i=begin;i<end;++i) {
            const double correction=a.densityCorrection ? a.eta[i]*a.dLnN2[i/a.planeSize] : 0;
            a.flux[i]-=a.velocity[i]*(derivative[i]+correction);
        }
    }
};
} // namespace wavevortex::kernel_detail
