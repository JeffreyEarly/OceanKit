#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {
using namespace wavevortex;
using Clock = std::chrono::steady_clock;

struct Result {
    std::string provider;
    std::string library;
    double medianSeconds = 0.0;
    std::size_t persistentBytes = 0;
    std::size_t scratchBytes = 0;
    std::size_t stateBytes = 0;
    std::size_t fullSpectrumBytes = 0;
    std::size_t planCount = 0;
    std::vector<WVComplex64> flux;
};

bool parseCount(const char* text, std::size_t& value) {
    if (text == nullptr || *text == '\0' || *text == '-') return false;
    char* end = nullptr;
    const auto parsed = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0' ||
        parsed > std::numeric_limits<std::size_t>::max())
        return false;
    value = static_cast<std::size_t>(parsed);
    return true;
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

WVKernelStatus benchmark(const WVTransformBarotropicQGConfiguration& config,
                         std::unique_ptr<WVFFTEngine> engine,
                         std::size_t warmups, std::size_t samples,
                         Result& result) {
    std::unique_ptr<WVTransformBarotropicQGKernel> kernel;
    auto status = WVTransformBarotropicQGKernel::create(
        config, std::move(engine), kernel);
    if (!status) return status;
    auto A0 = coefficients(kernel->descriptor());
    result.flux.resize(A0.size());
    const WVComplexConstView input{
        A0.data(), kernel->descriptor().spectralShape()};
    WVComplexView output{
        result.flux.data(), kernel->descriptor().spectralShape()};
    for (std::size_t index = 0; index < warmups; ++index) {
        status = kernel->nonlinearFlux(input, output);
        if (!status) return status;
    }
    std::vector<double> timings(samples);
    for (auto& timing : timings) {
        const auto start = Clock::now();
        status = kernel->nonlinearFlux(input, output);
        timing = std::chrono::duration<double>(Clock::now() - start).count();
        if (!status) return status;
    }
    std::sort(timings.begin(), timings.end());
    result.provider = kernel->engineIdentifier();
    result.library = kernel->engineLibraryIdentity();
    result.medianSeconds = timings[timings.size() / 2];
    result.persistentBytes = kernel->persistentBytes();
    result.scratchBytes = kernel->scratchBytes();
    result.stateBytes = A0.size() * sizeof(WVComplex64);
    result.fullSpectrumBytes =
        kernel->metrics().persistentFullHermitianBytes;
    result.planCount = kernel->metrics().planCount;
    return WVKernelStatus::ok();
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

void emitResult(const Result& result) {
    std::cout << "{\"provider\":\"" << result.provider
              << "\",\"library\":\"" << result.library
              << "\",\"medianSeconds\":" << result.medianSeconds
              << ",\"persistentBytes\":" << result.persistentBytes
              << ",\"scratchBytes\":" << result.scratchBytes
              << ",\"compactStateBytes\":" << result.stateBytes
              << ",\"persistentFullHermitianBytes\":"
              << result.fullSpectrumBytes
              << ",\"planCount\":" << result.planCount << '}';
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: benchmark Nx Ny warmups samples\n";
        return 2;
    }
    std::size_t Nx = 0;
    std::size_t Ny = 0;
    std::size_t warmups = 0;
    std::size_t samples = 0;
    if (!parseCount(argv[1], Nx) || !parseCount(argv[2], Ny) ||
        !parseCount(argv[3], warmups) || !parseCount(argv[4], samples) ||
        Nx < 2 || Ny < 2 || samples == 0) {
        std::cerr << "invalid benchmark arguments\n";
        return 2;
    }
    WVTransformBarotropicQGConfiguration configuration;
    configuration.Nx = Nx;
    configuration.Ny = Ny;
    configuration.Lx = 17000.0;
    configuration.Ly = 11000.0;
    configuration.h = 0.8;
    configuration.j = 1;
    configuration.g = 9.81;
    configuration.planetaryRadius = 6.371e6;
    configuration.rotationRate = 7.2921e-5;
    configuration.latitude = 33.0;
    configuration.shouldAntialias = true;

    std::unique_ptr<WVFFTEngine> fftw;
    auto status = WVFFTWEngine::create(1, fftw);
    if (!status) {
        std::cerr << status.message << '\n';
        return 3;
    }
    const auto before = WVFFTWEngine::lifetimeMetrics();
    Result native;
    status = benchmark(configuration, std::move(fftw), warmups, samples,
                       native);
    if (!status) {
        std::cerr << status.message << '\n';
        return 3;
    }
    const auto after = WVFFTWEngine::lifetimeMetrics();
    Result reference;
    status = benchmark(
        configuration,
        std::make_unique<WVReferenceFFTEngine>(),
        warmups, samples, reference);
    if (!status) {
        std::cerr << status.message << '\n';
        return 3;
    }
    const double error = relativeError(native.flux, reference.flux);
    if (error > 1e-12 || after.activePlans != before.activePlans ||
        after.outstandingPlanningBytes != 0) {
        std::cerr << "provider parity or cleanup failed\n";
        return 4;
    }

    std::cout << std::setprecision(17)
              << "{\"schemaVersion\":\"barotropic-qg-same-host-v1\","
              << "\"host\":\"donut\",\"Nx\":" << Nx
              << ",\"Ny\":" << Ny << ",\"j\":1,"
              << "\"shouldAntialias\":true,\"native\":";
    emitResult(native);
    std::cout << ",\"reference\":";
    emitResult(reference);
    std::cout << ",\"nativeToReferenceSpeedup\":"
              << reference.medianSeconds / native.medianSeconds
              << ",\"relativeError\":" << error
              << ",\"nativePlanCleanup\":true}\n";
    return 0;
}
