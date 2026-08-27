#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
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
    std::cout << "Native FFTW Barotropic QG tests passed\n";
    return 0;
}
