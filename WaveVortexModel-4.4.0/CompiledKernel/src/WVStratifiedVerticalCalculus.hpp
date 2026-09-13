#pragma once

#include "WVSpectralValidation.hpp"
#include "WaveVortexKernel/WVStratifiedModalSource.hpp"
#include <algorithm>

namespace wavevortex::kernel_detail {

template<class InputValue,class OutputValue,class VerticalOperation>
WVKernelStatus applyStratifiedVerticalCalculus(const WVStratifiedModalGeometry& geometry,
    std::size_t columns,bool inputIsF,unsigned order,bool integral,
    InputValue&& inputValue,OutputValue&& outputValue,WVComplexOutput firstGrid,
    WVComplexOutput secondGrid,WVComplexOutput modal,VerticalOperation&& vertical) {
    using namespace spectral_detail;
    const auto H=geometry.Nz*geometry.Nkl,S=geometry.Nj*geometry.Nkl;
    for (std::size_t begin=0;begin<columns;begin+=geometry.Nkl) {
        auto a=firstGrid,b=secondGrid; const auto count=std::min(geometry.Nkl,columns-begin);
        for (std::size_t i=0;i<H;++i) write(a,i,{});
        for (std::size_t column=0;column<count;++column) for (std::size_t z=0;z<geometry.Nz;++z) {
            double value=inputValue(begin+column,z);
            if (integral && !inputIsF) value/=geometry.N2[z];
            write(a,z+geometry.Nz*column,{value,0});
        }
        // 0: DzF, 1: DzG, 2: DzzG, 3: IntF, 4: IntG after N2 division.
        const auto apply=[&](unsigned operation) -> WVKernelStatus {
            auto status=vertical(operation==0 || operation==3 ? 1 : 3,a.input(),modal); if (!status) return status;
            for (std::size_t i=0;i<S;++i) {
                auto value=read(modal.input(),i);
                if (operation==1 || operation==2) {
                    const double factor=1/geometry.h_0[i%geometry.Nj];
                    write(modal,i,{value.real*factor,value.imag*factor});
                }
                if (operation==3) {
                    const double factor=geometry.h_0[i%geometry.Nj];
                    write(modal,i,{value.real*factor,value.imag*factor});
                }
            }
            status=vertical(operation==1 || operation==4 ? 0 : 2,modal.input(),b); if (!status) return status;
            for (std::size_t column=0;column<geometry.Nkl;++column) {
                const auto bottom=read(b.input(),geometry.Nz*column);
                for (std::size_t z=0;z<geometry.Nz;++z) {
                    auto value=read(b.input(),z+geometry.Nz*column);
                    if (operation==0 || operation==2) {
                        const double factor=-geometry.N2[z]/geometry.g;
                        value={value.real*factor,value.imag*factor};
                    }
                    if (operation==4) {
                        const double factor=-geometry.g;
                        value={(value.real-bottom.real)*factor,(value.imag-bottom.imag)*factor};
                    }
                    write(b,z+geometry.Nz*column,value);
                }
            }
            std::swap(a,b); return WVKernelStatus::ok();
        };
        WVKernelStatus status;
        if (integral) status=apply(inputIsF ? 3 : 4);
        else if (inputIsF) {
            status=apply(0); if (!status) return status;
            if (order>=3) { status=apply(2); if (!status) return status; }
            if (order==2 || order==4) status=apply(1);
        } else {
            if (order==1) status=apply(1);
            else { status=apply(2); if (!status) return status; if (order==3) status=apply(1); if (order==4) status=apply(2); }
        }
        if (!status) return status;
        for (std::size_t column=0;column<count;++column) for (std::size_t z=0;z<geometry.Nz;++z)
            outputValue(begin+column,z,read(a.input(),z+geometry.Nz*column).real);
    }
    return WVKernelStatus::ok();
}

} // namespace wavevortex::kernel_detail
