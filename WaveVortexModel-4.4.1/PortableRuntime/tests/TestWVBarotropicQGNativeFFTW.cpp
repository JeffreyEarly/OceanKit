#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace wavevortex;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

WVTransformBarotropicQGConfiguration configuration(
    std::size_t Nx, std::size_t Ny, std::uint32_t j,
    bool shouldAntialias) {
    WVTransformBarotropicQGConfiguration value;
    value.Nx = Nx;
    value.Ny = Ny;
    value.Lx = 17000.0;
    value.Ly = 11000.0;
    value.h = 0.8;
    value.j = j;
    value.g = 9.81;
    value.planetaryRadius = 6.371e6;
    value.rotationRate = 7.2921e-5;
    value.latitude = 33.0;
    value.shouldAntialias = shouldAntialias;
    return value;
}

std::vector<WVComplex64> coefficients(
    const WVTransformBarotropicQGDescriptor& descriptor) {
    std::vector<WVComplex64> values(descriptor.Nkl());
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] = {
            2e-5 * std::sin(0.31 * static_cast<double>(index + 1)),
            1e-5 * std::cos(0.17 * static_cast<double>(index + 3))};
        const auto& mode = descriptor.fourierModes()[index];
        if (mode.Kh == 0.0) values[index] = {};
        if (mode.dftPrimaryIndex == mode.dftConjugateIndex)
            values[index].imag = 0.0;
    }
    return values;
}

double relativeError(const std::vector<WVComplex64>& actual,
                     const std::vector<WVComplex64>& expected) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        numerator = std::max(
            numerator,
            std::hypot(actual[index].real - expected[index].real,
                       actual[index].imag - expected[index].imag));
        denominator = std::max(
            denominator,
            std::hypot(expected[index].real, expected[index].imag));
    }
    return numerator / std::max(denominator,
                                std::numeric_limits<double>::min());
}

struct ConsumerProbe {
    explicit ConsumerProbe(std::size_t elements,std::size_t planeElements)
        : captured(elements),visited(std::make_unique<std::atomic<unsigned>[]>(elements)),
          elementCount(elements),planeSize(planeElements) {
        for(std::size_t i=0;i<elementCount;++i) visited[i].store(0);
    }
    std::vector<double> captured;
    std::unique_ptr<std::atomic<unsigned>[]> visited;
    const std::size_t elementCount,planeSize;
    std::atomic<std::size_t> callbacks{0};
    std::atomic<bool> valid{true};
    static void consume(void* context,std::size_t begin,std::size_t end,
        const double* output) noexcept {
        auto& probe=*static_cast<ConsumerProbe*>(context);
        if(begin>=end || end>probe.elementCount || begin%probe.planeSize ||
            end-begin!=probe.planeSize) {
            probe.valid.store(false); return;
        }
        probe.callbacks.fetch_add(1);
        for(std::size_t i=begin;i<end;++i) {
            if(probe.visited[i].fetch_add(1)!=0) probe.valid.store(false);
            else probe.captured[i]=output[i];
        }
    }
    WVRealOutputConsumer consumer() { return {this,consume}; }
    bool complete(const std::vector<double>& output,std::size_t planes) const {
        if(!valid.load() || callbacks.load()!=planes || captured!=output) return false;
        for(std::size_t i=0;i<elementCount;++i) if(visited[i].load()!=1) return false;
        return true;
    }
};

WVRetainedHorizontalSpecification retainedSpecification(std::size_t workers) {
    constexpr std::size_t Nx=8,Ny=6,planes=5;
    WVRetainedHorizontalSpecification spec;
    spec.grid={Nx,Ny,planes,1,Nx,Nx*Ny,"consumer-grid"};
    spec.modes={{0,0},{1,0},{0,1},{1,1},{-1,1},{4,0},{0,3}};
    spec.retained={planes,spec.modes.size(),1,planes,
        WVComplexRepresentation::interleaved,"consumer-grid","consumer-modes"};
    spec.Lx=17000; spec.Ly=11000;
    spec.schedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
    spec.outerWorkers=workers;
    return spec;
}

void testInverseConsumers() {
    const auto referenceSpec=retainedSpecification(1);
    const auto plane=referenceSpec.grid.Nx*referenceSpec.grid.Ny;
    const auto elements=plane*referenceSpec.grid.planes;
    const auto spectralElements=referenceSpec.retained.rows*referenceSpec.retained.columns;
    std::vector<double> input(elements);
    for(std::size_t i=0;i<elements;++i)
        input[i]=.03*std::sin(.11*i)+.02*std::cos(.07*i);

    std::unique_ptr<WVRetainedHorizontalOperator> fallback;
    require(bool(WVRetainedHorizontalOperator::create(referenceSpec,
        std::make_unique<WVReferenceFFTEngine>(),fallback)),
        "reference inverse-consumer operator creation");
    std::unique_ptr<WVRetainedHorizontalWorkspace> fallbackWorkspace;
    require(bool(fallback->createWorkspace(fallbackWorkspace,false)) &&
        std::string(fallbackWorkspace->scheduleIdentifier())=="plane-streamed-full-fft-gather",
        "reference inverse-consumer fallback workspace");
    std::vector<WVComplex64> fallbackSpectral(spectralElements);
    require(bool(fallback->forward(*fallbackWorkspace,{input.data(),input.size()*sizeof(double)},
        {fallbackSpectral.data(),nullptr,nullptr,fallbackSpectral.size()*sizeof(WVComplex64)})),
        "reference inverse-consumer forward");
    const auto fallbackSpectralBefore=fallbackSpectral;
    std::vector<double> ordinary(elements),consumed(elements);
    require(bool(fallback->inverse(*fallbackWorkspace,
        {fallbackSpectral.data(),nullptr,nullptr,fallbackSpectral.size()*sizeof(WVComplex64)},
        {ordinary.data(),ordinary.size()*sizeof(double)})),
        "reference ordinary inverse");
    ConsumerProbe fallbackProbe(elements,plane);
    const auto fallbackConsumer=fallbackProbe.consumer();
    require(bool(fallback->inverseAndConsume(*fallbackWorkspace,
        {fallbackSpectral.data(),nullptr,nullptr,fallbackSpectral.size()*sizeof(WVComplex64)},
        {consumed.data(),consumed.size()*sizeof(double)},fallbackConsumer)) &&
        fallbackProbe.complete(consumed,referenceSpec.grid.planes) && consumed==ordinary &&
        std::memcmp(fallbackSpectral.data(),fallbackSpectralBefore.data(),
            fallbackSpectral.size()*sizeof(WVComplex64))==0,
        "reference fallback changed normalization, input, or callback coverage");

    for(const std::size_t requestedWorkers:{1u,2u,8u}) {
        const auto spec=retainedSpecification(requestedWorkers);
        std::unique_ptr<WVFFTEngine> engine;
        require(bool(WVFFTWEngine::create(1,engine)),"native inverse-consumer engine creation");
        std::unique_ptr<WVRetainedHorizontalOperator> native;
        require(bool(WVRetainedHorizontalOperator::create(spec,std::move(engine),native)),
            "native inverse-consumer operator creation");
        std::unique_ptr<WVRetainedHorizontalWorkspace> workspace;
        require(bool(native->createWorkspace(workspace,false)) &&
            workspace->workerCount()==std::min(requestedWorkers,spec.grid.planes),
            "native inverse-consumer compact worker count");
        std::vector<WVComplex64> spectral(spectralElements);
        require(bool(native->forward(*workspace,{input.data(),input.size()*sizeof(double)},
            {spectral.data(),nullptr,nullptr,spectral.size()*sizeof(WVComplex64)})),
            "native inverse-consumer forward");
        const auto spectralBefore=spectral;
        std::vector<double> output(elements);
        ConsumerProbe probe(elements,plane); const auto consumer=probe.consumer();
        require(bool(native->inverseAndConsume(*workspace,
            {spectral.data(),nullptr,nullptr,spectral.size()*sizeof(WVComplex64)},
            {output.data(),output.size()*sizeof(double)},consumer)) &&
            probe.complete(output,spec.grid.planes) &&
            std::memcmp(spectral.data(),spectralBefore.data(),spectral.size()*sizeof(WVComplex64))==0,
            "native inverse consumer changed input or delivered invalid concurrent ranges");
        double scale=0,error=0;
        for(std::size_t i=0;i<elements;++i) {
            scale=std::max(scale,std::abs(ordinary[i]));
            error=std::max(error,std::abs(output[i]-ordinary[i]));
        }
        require(error<=1e-12*std::max(1.0,scale),
            "native inverse consumer differs from normalized fallback output");
    }

    ConsumerProbe rejected(elements,plane); const auto rejectedConsumer=rejected.consumer();
    std::vector<double> untouched(elements,719);
    require(fallback->inverseAndConsume(*fallbackWorkspace,
        {fallbackSpectral.data(),nullptr,nullptr,fallbackSpectral.size()*sizeof(WVComplex64)},
        {untouched.data(),untouched.size()*sizeof(double)-sizeof(double)},rejectedConsumer).code==
            WVKernelStatusCode::invalidShape && rejected.callbacks.load()==0 &&
        std::all_of(untouched.begin(),untouched.end(),[](double value) { return value==719; }),
        "malformed inverse output reached the consumer or mutated output");
}

void testConfiguration(std::size_t Nx, std::size_t Ny, std::uint32_t j,
                       bool shouldAntialias) {
    const auto before = WVFFTWEngine::lifetimeMetrics();
    std::unique_ptr<WVFFTEngine> nativeEngine;
    auto status = WVFFTWEngine::create(1, nativeEngine);
    require(static_cast<bool>(status), "native engine creation");
    std::unique_ptr<WVTransformBarotropicQGKernel> native;
    status = WVTransformBarotropicQGKernel::create(
        configuration(Nx, Ny, j, shouldAntialias),
        std::move(nativeEngine), native);
    require(static_cast<bool>(status) && native,
            "native Barotropic QG kernel creation");
    const auto during = WVFFTWEngine::lifetimeMetrics();
    require(native->engineIdentifier() == "fftw" &&
                native->metrics().planCount == 3 &&
                during.activePlans == before.activePlans + 3 &&
                during.outstandingPlanningBytes == 0,
            "native provider identity and plan ownership");

    std::unique_ptr<WVTransformBarotropicQGKernel> reference;
    status = WVTransformBarotropicQGKernel::create(
        configuration(Nx, Ny, j, shouldAntialias),
        std::make_unique<WVReferenceFFTEngine>(),
        reference);
    require(static_cast<bool>(status) && reference,
            "reference Barotropic QG kernel creation");
    const auto A0 = coefficients(native->descriptor());
    const WVComplexConstView input{
        A0.data(), native->descriptor().spectralShape()};
    std::vector<WVComplex64> nativeFlux(A0.size());
    std::vector<WVComplex64> referenceFlux(A0.size());
    WVComplexView nativeOutput{
        nativeFlux.data(), native->descriptor().spectralShape()};
    WVComplexView referenceOutput{
        referenceFlux.data(), reference->descriptor().spectralShape()};
    status = native->nonlinearFlux(input, nativeOutput);
    require(static_cast<bool>(status), "native nonlinear PV tendency");
    status = reference->nonlinearFlux(input, referenceOutput);
    require(static_cast<bool>(status), "reference nonlinear PV tendency");
    require(relativeError(nativeFlux, referenceFlux) <= 1e-12,
            "native/reference nonlinear parity");
    require(native->metrics().persistentFullHermitianBytes == 0 &&
                native->scratchBytes() == reference->scratchBytes(),
            "native compact scratch and no retained full spectrum");

    native.reset();
    const auto after = WVFFTWEngine::lifetimeMetrics();
    require(after.activePlans == before.activePlans &&
                after.totalPlansCreated - before.totalPlansCreated == 3 &&
                after.totalPlansDestroyed - before.totalPlansDestroyed == 3 &&
                after.outstandingPlanningBytes == 0,
            "native FFTW plan and planning-buffer cleanup");
}

} // namespace

int main() {
    const auto identity = WVFFTWEngine::linkedLibraries();
    const auto expected = std::filesystem::weakly_canonical(
        std::filesystem::path(WV_RUNTIME_EXPECTED_FFTW_ROOT));
    const auto base = std::filesystem::weakly_canonical(identity.baseLibrary);
    const auto threads =
        std::filesystem::weakly_canonical(identity.threadLibrary);
    const auto underExpected = [&expected](const std::filesystem::path& path) {
        return path.string().rfind(
            expected.string() +
                std::string(1, std::filesystem::path::preferred_separator),
            0) == 0;
    };
    require(identity.version.find("3.3.11") != std::string::npos &&
                underExpected(base) && underExpected(threads),
            "pinned FFTW 3.3.11 library identity");
    testConfiguration(8, 6, 0, false);
    testConfiguration(9, 7, 1, false);
    testConfiguration(8, 7, 0, true);
    testConfiguration(9, 6, 1, true);
    testInverseConsumers();
    std::cout << "Native FFTW Barotropic QG tests passed\n";
    return 0;
}
