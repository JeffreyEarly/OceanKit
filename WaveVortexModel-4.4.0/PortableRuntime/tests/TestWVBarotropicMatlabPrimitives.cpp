#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace wavevortex;

namespace {
constexpr double pi=3.141592653589793238462643383279502884;

void require(bool condition,const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

WVTransformBarotropicQGConfiguration configuration() {
    WVTransformBarotropicQGConfiguration value;
    value.Nx=8; value.Ny=6; value.Lx=2*pi; value.Ly=2*pi;
    value.h=.8; value.j=1; value.g=9.81;
    value.planetaryRadius=6.371e6; value.rotationRate=7.2921e-5;
    value.latitude=33; value.shouldAntialias=true;
    return value;
}

double coordinate(std::size_t index,std::size_t count) {
    return 2*pi*static_cast<double>(index)/static_cast<double>(count);
}

std::complex<double> dft(const std::vector<double>& field,
    const WVTransformBarotropicQGConfiguration& c,std::int64_t k,std::int64_t l) {
    std::complex<double> result{};
    for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x) {
        const double phase=-(static_cast<double>(k)*coordinate(x,c.Nx)+
            static_cast<double>(l)*coordinate(y,c.Ny));
        result+=field[x+c.Nx*y]*std::exp(std::complex<double>(0,phase));
    }
    return result/static_cast<double>(c.Nx*c.Ny);
}

double reconstruction(const std::vector<WVComplex64>& spectrum,
    const WVTransformBarotropicQGDescriptor& descriptor,std::size_t x,std::size_t y) {
    const auto& c=descriptor.configuration(); double result=0;
    for (std::size_t i=0;i<spectrum.size();++i) {
        const auto& mode=descriptor.fourierModes()[i];
        const std::complex<double> value{spectrum[i].real,spectrum[i].imag};
        const double phase=static_cast<double>(mode.kMode)*coordinate(x,c.Nx)+
            static_cast<double>(mode.lMode)*coordinate(y,c.Ny);
        const auto contribution=value*std::exp(std::complex<double>(0,phase));
        result+=(mode.dftPrimaryIndex==mode.dftConjugateIndex ? 1.0 : 2.0)*contribution.real();
    }
    return result;
}

double cosineDerivative(double amplitude,unsigned order,double mode,double phase) {
    return amplitude*std::real(std::pow(std::complex<double>(0,mode),
        static_cast<int>(order))*std::exp(std::complex<double>(0,mode*phase)));
}

double sineDerivative(double amplitude,unsigned order,double mode,double phase) {
    return amplitude*std::imag(std::pow(std::complex<double>(0,mode),
        static_cast<int>(order))*std::exp(std::complex<double>(0,mode*phase)));
}

void testRawTransformsAndDerivatives() {
    const auto c=configuration();
    std::unique_ptr<WVTransformBarotropicQGKernel> kernel;
    auto status=WVTransformBarotropicQGKernel::create(c,
        std::make_unique<WVReferenceFFTEngine>(),kernel);
    require(bool(status) && kernel,"Barotropic reference kernel creation");
    const auto spatial=kernel->descriptor().spatialShape();
    const auto spectral=kernel->descriptor().spectralShape();
    const auto R=spatial.elementCount(),S=spectral.elementCount();
    const auto persistent=kernel->persistentBytes(),scratch=kernel->scratchBytes();
    const auto plans=kernel->metrics().planCount,planBytes=kernel->metrics().planBytes;

    std::vector<double> field(R);
    for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x) {
        const double xx=coordinate(x,c.Nx),yy=coordinate(y,c.Ny);
        field[x+c.Nx*y]=3+.7*std::cos(xx)-.4*std::sin(2*yy)+
            .23*std::cos(3*xx)+.11*std::cos(4*xx)+.09*std::cos(3*yy);
    }
    std::vector<WVComplex64> spectrum(S),spectrumBefore;
    WVComplexView spectrumView{spectrum.data(),spectral};
    status=kernel->horizontalForward({field.data(),spatial},spectrumView);
    require(bool(status),"Raw horizontal forward");
    bool sawMean=false,sawRetained=false;
    for (std::size_t i=0;i<S;++i) {
        const auto& mode=kernel->descriptor().fourierModes()[i];
        const auto expected=dft(field,c,mode.kMode,mode.lMode);
        require(std::abs(spectrum[i].real-expected.real())<2e-13 &&
            std::abs(spectrum[i].imag-expected.imag())<2e-13,
            "Raw forward differs from direct retained DFT");
        if (mode.kMode==0 && mode.lMode==0) {
            sawMean=true;
            require(std::abs(spectrum[i].real-3)<2e-13 && std::abs(spectrum[i].imag)<2e-13,
                "Raw forward masked the horizontal mean");
        }
        if (mode.kMode==1 && mode.lMode==0) sawRetained=true;
    }
    require(sawMean && sawRetained,"Expected mean and retained mode are absent");
    spectrumBefore=spectrum;
    std::vector<double> restored(R);
    WVRealView restoredView{restored.data(),spatial};
    status=kernel->horizontalInverse({spectrum.data(),spectral},restoredView);
    require(bool(status),"Raw horizontal inverse");
    require(std::memcmp(spectrum.data(),spectrumBefore.data(),S*sizeof(WVComplex64))==0,
        "Raw inverse modified compact input");
    double omittedDifference=0,restoredMean=0;
    for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x) {
        const auto index=x+c.Nx*y;
        const double expected=reconstruction(spectrum,kernel->descriptor(),x,y);
        require(std::abs(restored[index]-expected)<2e-12,
            "Raw inverse differs from retained reconstruction");
        restoredMean+=restored[index]/static_cast<double>(R);
        omittedDifference=std::max(omittedDifference,std::abs(restored[index]-field[index]));
    }
    require(std::abs(restoredMean-3)<2e-13,"Raw inverse masked the horizontal mean");
    require(omittedDifference>.1,"Raw inverse unexpectedly restored omitted modes");

    std::vector<double> derivative(R),first(R);
    WVRealView derivativeView{derivative.data(),spatial};
    WVRealView firstView{first.data(),spatial};
    for (unsigned order:{1U,2U,4U}) for (bool xDerivative:{true,false}) {
        status=kernel->differentiateHorizontal({field.data(),spatial},xDerivative,order,
            derivativeView);
        require(bool(status),"Raw full-grid derivative");
        for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x) {
            const auto index=x+c.Nx*y;
            const double xx=coordinate(x,c.Nx),yy=coordinate(y,c.Ny);
            double expected=0;
            if (xDerivative) {
                expected+=cosineDerivative(.7,order,1,xx);
                expected+=cosineDerivative(.23,order,3,xx);
                if (order%2==0) expected+=cosineDerivative(.11,order,4,xx);
            } else {
                expected+=sineDerivative(-.4,order,2,yy);
                if (order%2==0) expected+=cosineDerivative(.09,order,3,yy);
            }
            require(std::abs(derivative[index]-expected)<3e-11,
                "Raw derivative disagrees with full-grid DFT semantics");
        }
        if (order==1) {
            status=kernel->differentiateHorizontal({field.data(),spatial},xDerivative,
                firstView);
            require(bool(status) && std::memcmp(first.data(),derivative.data(),R*sizeof(double))==0,
                "Order-one convenience path changed arithmetic");
        }
    }
    status=kernel->differentiateHorizontal({field.data(),spatial},true,0,
        derivativeView);
    require(status.code==WVKernelStatusCode::invalidConfiguration,
        "Zero derivative order was accepted");

    const auto sharedCount=std::max(S,(R*sizeof(double)+sizeof(WVComplex64)-1)/sizeof(WVComplex64));
    std::vector<WVComplex64> shared(sharedCount);
    auto* sharedReal=reinterpret_cast<double*>(shared.data());
    WVComplexView sharedSpectrum{shared.data(),spectral};
    WVRealView sharedSpatial{sharedReal,spatial};
    require(kernel->horizontalForward({sharedReal,spatial},sharedSpectrum).code==
        WVKernelStatusCode::overlappingArrays,"Raw forward accepted aliased storage");
    require(kernel->horizontalInverse({shared.data(),spectral},sharedSpatial).code==
        WVKernelStatusCode::overlappingArrays,"Raw inverse accepted aliased storage");
    require(kernel->differentiateHorizontal({sharedReal,spatial},true,2,
        sharedSpatial).code==WVKernelStatusCode::overlappingArrays,
        "Raw derivative accepted aliased storage");

    require(kernel->persistentBytes()==persistent && kernel->scratchBytes()==scratch &&
        kernel->metrics().planCount==plans && kernel->metrics().planBytes==planBytes,
        "Raw primitives changed persistent storage or FFT plan ownership");
}

void testActiveStateOutputProtection() {
    const auto c=configuration();
    std::unique_ptr<WVTransformBarotropicQGKernel> kernel;
    require(bool(WVTransformBarotropicQGKernel::create(c,
        std::make_unique<WVReferenceFFTEngine>(),kernel)),"State-protection kernel creation");
    const auto spatial=kernel->descriptor().spatialShape();
    const auto spectral=kernel->descriptor().spectralShape();
    const auto stateElements=std::max(spectral.elementCount(),
        (spatial.elementCount()*sizeof(double)+sizeof(WVComplex64)-1)/sizeof(WVComplex64));
    std::vector<WVComplex64> state(stateElements);
    std::vector<double> field(spatial.elementCount());
    require(bool(kernel->beginStateEvaluation({state.data(),spectral})),"Begin active state");
    WVComplexView stateOutput{state.data(),spectral};
    WVRealView overlappingOutput{reinterpret_cast<double*>(state.data()),spatial};
    require(kernel->horizontalForward({field.data(),spatial},stateOutput).code==
        WVKernelStatusCode::overlappingArrays,"Raw forward overwrote active state");
    require(kernel->horizontalInverse({state.data(),spectral},overlappingOutput).code==
        WVKernelStatusCode::overlappingArrays,"Raw inverse overwrote active state");
    require(bool(kernel->endStateEvaluation()),"End active state");
}
}

int main() {
    testRawTransformsAndDerivatives();
    testActiveStateOutputProtection();
    std::cout << "Barotropic MATLAB primitive contracts passed\n";
    return 0;
}
