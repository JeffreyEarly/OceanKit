#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVCheckpointWriter.hpp"
#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVRungeKutta.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVCheckpointWriterTestHooks.hpp"
#include "WVForcingImplementations.hpp"
#include "WVTestLinearCoefficientForcing.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path fixture(const std::string& name) {
    return std::filesystem::path(WV_CHECKPOINT_FIXTURE_DIR) / name;
}

std::filesystem::path temporaryDirectory() {
    const auto path = std::filesystem::temp_directory_path() / ("wave-vortex-checkpoint-writer-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(path);
    return path;
}

struct DirectoryCleanup {
    std::filesystem::path path;
    ~DirectoryCleanup() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
};

std::vector<char> bytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

bool sameComplex(const std::vector<WVComplex64>& left, const std::vector<WVComplex64>& right) {
    if (left.size() != right.size()) return false;
    for (std::size_t index = 0; index < left.size(); ++index) {
        const bool sameReal = left[index].real == right[index].real || (std::isnan(left[index].real) && std::isnan(right[index].real));
        const bool sameImag = left[index].imag == right[index].imag || (std::isnan(left[index].imag) && std::isnan(right[index].imag));
        if (!sameReal || !sameImag) return false;
    }
    return true;
}

bool sameTypedRecord(const WVPortableTypedRecord& left, const WVPortableTypedRecord& right) {
    if (left.schemaIdentifier != right.schemaIdentifier || left.schemaVersion != right.schemaVersion || left.values.size() != right.values.size()) return false;
    for (std::size_t index = 0; index < left.values.size(); ++index) {
        const auto& a = left.values[index];
        const auto& b = right.values[index];
        if (a.name != b.name || a.dimensions != b.dimensions || a.storage != b.storage) return false;
    }
    return true;
}

void requireSameCheckpoint(const WVCheckpoint& expected, const WVCheckpoint& actual) {
    const auto& a = expected.configuration;
    const auto& b = actual.configuration;
    require(a.Nx == b.Nx && a.Ny == b.Ny && a.Nz == b.Nz && a.Nj == b.Nj && a.Lx == b.Lx && a.Ly == b.Ly && a.Lz == b.Lz && a.N0 == b.N0 && a.rho0 == b.rho0 && a.g == b.g && a.planetaryRadius == b.planetaryRadius && a.rotationRate == b.rotationRate && a.latitude == b.latitude && a.isHydrostatic == b.isHydrostatic && a.shouldAntialias == b.shouldAntialias, "configuration changed during checkpoint writing");
    require(expected.state.t == actual.state.t && expected.state.t0 == actual.state.t0, "checkpoint time changed during writing");
    require(expected.state.coefficients.shape.rows == actual.state.coefficients.shape.rows && expected.state.coefficients.shape.columns == actual.state.coefficients.shape.columns, "coefficient shape changed during writing");
    require(sameComplex(expected.state.coefficients.Ap, actual.state.coefficients.Ap) && sameComplex(expected.state.coefficients.Am, actual.state.coefficients.Am) && sameComplex(expected.state.coefficients.A0, actual.state.coefficients.A0), "coefficient values changed during writing");
    require(expected.forcingSchedule.entries.size() == actual.forcingSchedule.entries.size(), "forcing count changed during writing");
    for (std::size_t index = 0; index < expected.forcingSchedule.entries.size(); ++index) {
        const auto& left = expected.forcingSchedule.entries[index];
        const auto& right = actual.forcingSchedule.entries[index];
        require(left.typeIdentifier == right.typeIdentifier && left.contractVersion == right.contractVersion && left.name == right.name && left.stage == right.stage && left.priority == right.priority && left.ordinal == right.ordinal && sameTypedRecord(left.configuration,right.configuration), "forcing metadata changed during writing");
    }
}

WVCheckpoint read(const std::filesystem::path& path) {
    WVCheckpoint checkpoint;
    const auto result = WVCheckpointReader::read(path.string(), *test::extensionCatalog(), checkpoint);
    require(static_cast<bool>(result), result.message);
    return checkpoint;
}

void testRoundTrips() {
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    const std::array<const char*, 12> names = {
        "root-hydrostatic.nc", "root-nonhydrostatic.nc", "time-series-hydrostatic.nc", "time-series-nonhydrostatic.nc",
        "forcing-nonlinear.nc", "forcing-adaptive-damping.nc", "forcing-fixed-amplitude.nc", "forcing-quadratic-bottom-friction.nc",
        "forcing-pseudo-topographic.nc", "forcing-beta-plane.nc", "forcing-mixed-hydrostatic.nc", "forcing-mixed-nonhydrostatic.nc"};
    for (const char* name : names) {
        const auto expected = read(fixture(name));
        const auto output = directory / name;
        const auto result = WVCheckpointWriter::write(output.string(), *test::extensionCatalog(), expected);
        require(static_cast<bool>(result), result.message);
        const auto actual = read(output);
        requireSameCheckpoint(expected, actual);
        require(actual.metadata.stateGroupPath == "/" && actual.metadata.stateCount == 1 && actual.metadata.selectedStateIndex == 0, "writer did not produce scalar root state");
    }
}

std::unique_ptr<WVConstantStratificationIntegrationSystem> system(const WVCheckpoint& checkpoint) {
    std::unique_ptr<WVConstantStratificationIntegrationSystem> result;
    const auto creation = WVConstantStratificationIntegrationSystem::create(checkpoint.configuration, checkpoint.forcingSchedule, test::extensionCatalog(), std::make_unique<WVReferenceFFTEngine>(), result);
    require(static_cast<bool>(creation), creation.message);
    return result;
}

WVMutableIntegrationState mutableState(WVCheckpoint& checkpoint) {
    const auto shape = checkpoint.state.coefficients.shape;
    return {{checkpoint.state.t, checkpoint.state.t0, {{checkpoint.state.coefficients.Ap.data(), shape}, {checkpoint.state.coefficients.Am.data(), shape}, {checkpoint.state.coefficients.A0.data(), shape}}},nullptr,0};
}

void synchronizeTime(WVCheckpoint& checkpoint, const WVMutableIntegrationState& state) {
    checkpoint.state.t = state.waveVortex.t;
    checkpoint.state.t0 = state.waveVortex.t0;
}

void testRestartContinuation() {
    for (const char* name : {"forcing-mixed-hydrostatic.nc", "forcing-mixed-nonhydrostatic.nc"}) {
        auto uninterrupted = read(fixture(name));
        auto prefix = uninterrupted;
        auto uninterruptedSystem = system(uninterrupted);
        auto prefixSystem = system(prefix);
        WVFixedStepRK4 uninterruptedIntegrator(*uninterruptedSystem);
        WVFixedStepRK4 prefixIntegrator(*prefixSystem);
        auto uninterruptedState = mutableState(uninterrupted);
        auto prefixState = mutableState(prefix);
        auto result = uninterruptedIntegrator.prepareStateAfterRestart(uninterruptedState);
        require(static_cast<bool>(result), result.message);
        result = prefixIntegrator.prepareStateAfterRestart(prefixState);
        require(static_cast<bool>(result), result.message);
        result = uninterruptedIntegrator.advanceToTime(uninterruptedState, uninterruptedState.waveVortex.t + 0.275, 0.1);
        require(static_cast<bool>(result), result.message);
        result = prefixIntegrator.advanceToTime(prefixState, prefixState.waveVortex.t + 0.2, 0.1);
        require(static_cast<bool>(result), result.message);
        synchronizeTime(prefix, prefixState);

        const auto directory = temporaryDirectory();
        DirectoryCleanup cleanup{directory};
        const auto path = directory / "restart.nc";
        const auto writeResult = WVCheckpointWriter::write(path.string(), *test::extensionCatalog(), prefix);
        require(static_cast<bool>(writeResult), writeResult.message);
        auto restarted = read(path);
        auto restartedSystem = system(restarted);
        WVFixedStepRK4 restartedIntegrator(*restartedSystem);
        auto restartedState = mutableState(restarted);
        result = restartedIntegrator.prepareStateAfterRestart(restartedState);
        require(static_cast<bool>(result), result.message);
        result = restartedIntegrator.advanceToTime(restartedState, uninterruptedState.waveVortex.t, 0.1);
        require(static_cast<bool>(result), result.message);
        require(restartedState.waveVortex.t == uninterruptedState.waveVortex.t && sameComplex(restarted.state.coefficients.Ap, uninterrupted.state.coefficients.Ap) && sameComplex(restarted.state.coefficients.Am, uninterrupted.state.coefficients.Am) && sameComplex(restarted.state.coefficients.A0, uninterrupted.state.coefficients.A0), "checkpoint restart continuation differs from uninterrupted RK4");
    }
}

void requireNoTemporaryFiles(const std::filesystem::path& directory) {
    for (const auto& entry : std::filesystem::directory_iterator(directory)) require(entry.path().filename().string().find(".tmp-") == std::string::npos, "checkpoint temporary file was retained after failure");
}

void testTransactionalFailures() {
    auto checkpoint = read(fixture("forcing-mixed-nonhydrostatic.nc"));
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    const auto destination = directory / "checkpoint.nc";
    const std::vector<char> sentinel = {'p','r','e','v','i','o','u','s'};
    {
        std::ofstream output(destination, std::ios::binary);
        output.write(sentinel.data(), static_cast<std::streamsize>(sentinel.size()));
    }
    for (const auto point : {detail::WVCheckpointWriterFailurePoint::afterDefinition, detail::WVCheckpointWriterFailurePoint::afterWrite, detail::WVCheckpointWriterFailurePoint::beforeCommit}) {
        detail::setCheckpointWriterFailurePoint(point);
        const auto result = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(), checkpoint);
        detail::setCheckpointWriterFailurePoint(detail::WVCheckpointWriterFailurePoint::none);
        require(!result, "injected checkpoint writer failure succeeded");
        require(bytes(destination) == sentinel, "injected checkpoint failure changed the prior destination");
        requireNoTemporaryFiles(directory);
    }
    auto malformed = checkpoint;
    malformed.state.coefficients.Ap.pop_back();
    const auto invalid = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(), malformed);
    require(invalid.code == WVCheckpointStatusCode::shapeMismatch, "invalid coefficient shape was not rejected before writing");
    require(bytes(destination) == sentinel, "validation failure changed the prior destination");
    requireNoTemporaryFiles(directory);

    const auto success = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(), checkpoint);
    require(static_cast<bool>(success), success.message);
    requireSameCheckpoint(checkpoint, read(destination));
}

void testCreateNewCommitPolicy() {
    const auto checkpoint = read(fixture("forcing-mixed-nonhydrostatic.nc"));
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    const auto destination = directory / "scheduled.nc";
    auto result = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(),checkpoint,WVCheckpointCommitPolicy::createNew);
    require(static_cast<bool>(result),result.message);
    requireSameCheckpoint(checkpoint,read(destination));
    const auto original = bytes(destination);
    result = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(),checkpoint,WVCheckpointCommitPolicy::createNew);
    require(result.code == WVCheckpointStatusCode::commitFailure,"create-new policy replaced an existing checkpoint");
    require(bytes(destination) == original,"create-new collision changed the existing checkpoint");
    requireNoTemporaryFiles(directory);

    const auto failedDestination = directory / "failed.nc";
    detail::setCheckpointWriterFailurePoint(detail::WVCheckpointWriterFailurePoint::beforeCommit);
    result = WVCheckpointWriter::write(failedDestination.string(), *test::extensionCatalog(),checkpoint,WVCheckpointCommitPolicy::createNew);
    detail::setCheckpointWriterFailurePoint(detail::WVCheckpointWriterFailurePoint::none);
    require(!result && !std::filesystem::exists(failedDestination),"failed create-new commit left a destination");
    requireNoTemporaryFiles(directory);
}

void testRegisteredPairRoundTrip() {
    auto checkpoint = read(fixture("forcing-fixed-amplitude.nc"));
    checkpoint.forcingSchedule.entries.front().typeIdentifier =
        "WVTestPortableFixedAmplitudeForcing";
    checkpoint.metadata.forcingHeaders.front().annotatedClass =
        "WVTestPortableFixedAmplitudeForcing";
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    const auto destination = directory / "paired-fixed-amplitude.nc";
    const auto result = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(), checkpoint);
    require(static_cast<bool>(result), result.message);
    requireSameCheckpoint(checkpoint, read(destination));
}

void testValidation() {
    auto checkpoint = read(fixture("forcing-mixed-nonhydrostatic.nc"));
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    auto malformed = checkpoint;
    malformed.metadata.modelVersion = "5.0.0";
    require(WVCheckpointWriter::write((directory / "version.nc").string(), *test::extensionCatalog(), malformed).code == WVCheckpointStatusCode::unsupportedModelVersion, "unsupported model version was accepted");
    malformed = checkpoint;
    malformed.forcingSchedule.entries.front().ordinal = 9;
    require(WVCheckpointWriter::write((directory / "ordinal.nc").string(), *test::extensionCatalog(), malformed).code == WVCheckpointStatusCode::invalidValue, "noncanonical forcing ordinal was accepted");
    malformed = checkpoint;
    malformed.forcingSchedule.entries.front().stage = WVForcingStage::spectral;
    require(WVCheckpointWriter::write((directory / "stage.nc").string(), *test::extensionCatalog(), malformed).code == WVCheckpointStatusCode::malformedForcing, "noncanonical forcing stage was accepted");
    require(WVCheckpointWriter::write("", *test::extensionCatalog(), checkpoint).code == WVCheckpointStatusCode::writeFailure, "empty destination was accepted");
}

void testRegisteredLinearCoefficientRoundTrip() {
    auto checkpoint = read(fixture("forcing-nonlinear.nc"));
    auto& forcing = checkpoint.forcingSchedule.entries.front();
    forcing.typeIdentifier = test::LinearCoefficientForcingIdentifier;
    forcing.name = "linear coefficient fixture";
    forcing.stage = WVForcingStage::spectral;
    forcing.priority = 90;
    forcing.configuration = {"wave-vortex-forcing-configuration-v1",1,
                             {{"rate",{},std::vector<double>{0.375}}}};
    const auto directory = temporaryDirectory();
    DirectoryCleanup cleanup{directory};
    const auto destination = directory / "linear-coefficient.nc";
    const auto result = WVCheckpointWriter::write(destination.string(), *test::extensionCatalog(),checkpoint);
    require(static_cast<bool>(result),result.message);
    const auto restored = read(destination);
    requireSameCheckpoint(checkpoint,restored);
    const auto* rate = restored.forcingSchedule.entries.front().configuration.value("rate");
    require(rate != nullptr && std::get<std::vector<double>>(rate->storage) == std::vector<double>{0.375},"registered linear coefficient forcing did not round-trip generically");
}

} // namespace

int main() {
    try {
        testRoundTrips();
        testRestartContinuation();
        testTransactionalFailures();
        testCreateNewCommitPolicy();
        testRegisteredPairRoundTrip();
        testRegisteredLinearCoefficientRoundTrip();
        testValidation();
        std::cout << "Portable checkpoint writer tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        detail::setCheckpointWriterFailurePoint(detail::WVCheckpointWriterFailurePoint::none);
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
