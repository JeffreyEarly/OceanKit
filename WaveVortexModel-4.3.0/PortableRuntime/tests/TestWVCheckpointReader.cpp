#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"
#include "WVTestLinearCoefficientForcing.hpp"
#include "WVTestExtensionCatalog.hpp"

#include <netcdf.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#ifndef WV_CHECKPOINT_FIXTURE_DIR
#error "WV_CHECKPOINT_FIXTURE_DIR must identify the committed fixture directory"
#endif

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void requireNetCDF(int status, const std::string& operation) {
    if (status != NC_NOERR) throw std::runtime_error(operation + ": " + nc_strerror(status));
}

void overwriteTextAttribute(int groupId, const char* name, const std::string& value);

std::filesystem::path fixture(const std::string& name) {
    return std::filesystem::path(WV_CHECKPOINT_FIXTURE_DIR) / name;
}

std::filesystem::path temporaryCopy(const std::string& name) {
    static std::size_t sequence = 0;
    const auto path = std::filesystem::temp_directory_path() / ("wv-checkpoint-test-" + std::to_string(++sequence) + ".nc");
    std::filesystem::copy_file(fixture(name), path, std::filesystem::copy_options::overwrite_existing);
    return path;
}

struct TemporaryFile {
    explicit TemporaryFile(std::filesystem::path value) : path(std::move(value)) {}
    ~TemporaryFile() { std::error_code ignored; std::filesystem::remove(path, ignored); }
    std::filesystem::path path;
};

WVCheckpoint read(const std::string& name, WVCheckpointStateSelection selection = WVCheckpointStateSelection::latest()) {
    WVCheckpoint checkpoint;
    const auto result = WVCheckpointReader::read(fixture(name).string(), *test::extensionCatalog(), checkpoint, selection);
    require(static_cast<bool>(result), result.message);
    return checkpoint;
}

void verifyCoefficient(const WVComplex64& value, const std::string& name, double offset, double index) {
    double expectedReal = 0.0;
    double expectedImag = 0.0;
    if (name == "Ap") {
        expectedReal = offset + index / 1000.0;
        expectedImag = -offset - index / 2000.0;
    } else if (name == "Am") {
        expectedReal = -2.0 * offset + index / 1500.0;
        expectedImag = 3.0 * offset - index / 2500.0;
    } else {
        expectedReal = 4.0 * offset - index / 3000.0;
        expectedImag = -5.0 * offset + index / 3500.0;
    }
    require(std::abs(value.real - expectedReal) < 1e-13 && std::abs(value.imag - expectedImag) < 1e-13,
        name + " coefficient did not preserve MATLAB column-major values");
}

template <typename Storage>
const Storage& storedValue(const WVFrozenForcingEntry& entry, const std::string& name) {
    const auto* value = entry.configuration.value(name);
    require(value != nullptr && std::holds_alternative<Storage>(value->storage), "forcing configuration value " + name + " is missing or has the wrong type");
    return std::get<Storage>(value->storage);
}

void verifyCheckpoint(const WVCheckpoint& checkpoint, bool hydrostatic, const std::string& groupPath, std::size_t stateCount, std::size_t selectedIndex, double time, double offset) {
    const auto& configuration = checkpoint.configuration;
    require(checkpoint.metadata.profileIdentifier == "wave-vortex-4x-v1", "profile identifier mismatch");
    require(checkpoint.metadata.modelVersion == "4.2.1", "model version mismatch");
    require(checkpoint.metadata.transformClass == "WVTransformConstantStratification", "transform class mismatch");
    require(checkpoint.metadata.stateGroupPath == groupPath, "state group mismatch");
    require(checkpoint.metadata.stateCount == stateCount && checkpoint.metadata.selectedStateIndex == selectedIndex, "state selection mismatch");
    require(configuration.Nx == 8 && configuration.Ny == 6 && configuration.Nz == 7 && configuration.Nj == 4, "configuration shape mismatch");
    require(configuration.Lx == 15000.0 && configuration.Ly == 12000.0 && configuration.Lz == 1300.0, "domain length mismatch");
    require(configuration.N0 == 5.2e-3 && configuration.rho0 == 1027.0 && configuration.g == 9.80665, "physical configuration mismatch");
    require(configuration.planetaryRadius == 6.3712e6 && configuration.rotationRate == 7.292115e-5 && configuration.latitude == 33.0, "rotating-planet configuration mismatch");
    require(configuration.isHydrostatic == hydrostatic && configuration.shouldAntialias, "logical configuration mismatch");
    require(checkpoint.state.t == time && checkpoint.state.t0 == -3.25, "checkpoint time mismatch");
    require(checkpoint.state.coefficients.shape.rows == 4 && checkpoint.state.coefficients.shape.columns == 9, "coefficient shape mismatch");
    require(checkpoint.state.coefficients.Ap.size() == 36 && checkpoint.state.coefficients.Am.size() == 36 && checkpoint.state.coefficients.A0.size() == 36, "coefficient storage mismatch");
    verifyCoefficient(checkpoint.state.coefficients.Ap.front(), "Ap", offset, 1.0);
    verifyCoefficient(checkpoint.state.coefficients.Ap.back(), "Ap", offset, 36.0);
    verifyCoefficient(checkpoint.state.coefficients.Am.front(), "Am", offset, 1.0);
    verifyCoefficient(checkpoint.state.coefficients.A0.back(), "A0", offset, 36.0);
    require(checkpoint.metadata.forcingHeaders.size() == 1, "default forcing metadata was not recovered");
    require(checkpoint.metadata.forcingHeaders.front().ordinal == 1 && checkpoint.metadata.forcingHeaders.front().groupPath == "/forcing" && checkpoint.metadata.forcingHeaders.front().annotatedClass == "WVNonlinearAdvection", "forcing header mismatch");
    require(checkpoint.forcingSchedule.profileIdentifier == "wave-vortex-forcing-v1" && checkpoint.forcingSchedule.entries.size() == 1, "forcing schedule was not decoded");
    const auto& forcing = checkpoint.forcingSchedule.entries.front();
    require(forcing.typeIdentifier == "WVNonlinearAdvection" && forcing.name == "nonlinear advection" && forcing.stage == WVForcingStage::spatial && forcing.priority == 127, "nonlinear forcing contract mismatch");

    WVTransformConstantStratificationDescriptor descriptor;
    const auto descriptorStatus = WVTransformConstantStratificationDescriptor::create(configuration, descriptor);
    require(static_cast<bool>(descriptorStatus), "checkpoint configuration did not rebuild a descriptor");
    require(descriptor.spectralShape().rows == 4 && descriptor.Nkl() == 9, "rebuilt descriptor shape mismatch");
    const WVState view = checkpoint.state.view();
    require(view.coefficients.Ap.data == checkpoint.state.coefficients.Ap.data() && view.coefficients.Ap.shape.rows == 4, "checkpoint state view mismatch");
}

void testPositiveFixtures() {
    verifyCheckpoint(read("root-nonhydrostatic.nc"), false, "/", 1, 0, 4.5, 10.0);
    verifyCheckpoint(read("root-hydrostatic.nc"), true, "/", 1, 0, 5.5, 30.0);
    verifyCheckpoint(read("time-series-nonhydrostatic.nc"), false, "/wave-vortex", 3, 2, 3.25, 300.0);
    verifyCheckpoint(read("time-series-hydrostatic.nc", WVCheckpointStateSelection::atIndex(0)), true, "/wave-vortex", 3, 0, 0.5, 120.0);
    verifyCheckpoint(read("time-series-hydrostatic.nc", WVCheckpointStateSelection::atIndex(1)), true, "/wave-vortex", 3, 1, 1.75, 220.0);
}

void testAllocationLightInspection() {
    WVCheckpointInspection inspection;
    const auto result = WVCheckpointReader::inspect(fixture("forcing-mixed-nonhydrostatic.nc").string(), *test::extensionCatalog(), inspection);
    require(static_cast<bool>(result), result.message);
    require(inspection.coefficientShape.rows == 4 && inspection.coefficientShape.columns == 9, "inspection coefficient shape mismatch");
    require(inspection.stateDescription.transformIdentifier ==
                    "WVTransformConstantStratification" &&
                inspection.stateDescription.spatialDimensions ==
                    std::vector<std::size_t>({8, 6, 7}) &&
                inspection.stateDescription.coefficientFamilies.size() == 3 &&
                inspection.stateDescription.coefficientFamilies[0].identifier ==
                    "Ap" &&
                inspection.stateDescription.coefficientFamilies[0]
                        .spectralDimensions ==
                    std::vector<std::size_t>({4, 9}),
            "inspection transform-state description mismatch");
    require(inspection.configuration.Nj == 4 && !inspection.configuration.isHydrostatic, "inspection configuration mismatch");
    require(inspection.forcingSchedule.entries.size() == 6, "inspection did not decode the frozen forcing schedule");
    require(inspection.t > inspection.t0, "inspection time metadata mismatch");
}

void testForcingCapabilities() {
    const auto& forcings = test::extensionCatalog()->forcings();
    const auto& capabilities = forcings.registrations();
    require(capabilities.size() == 15, "forcing capability matrix does not cover supplied and test classes");
    std::size_t supported = 0;
    std::set<std::string> identifiers;
    for (const auto& capability : capabilities) {
        require(identifiers.insert(capability.matlabClassName).second, "forcing capability identifier repeated");
        require(!capability.forcingTypes.empty(), "forcing capability omitted WVForcingType names");
        if (capability.isSupported) {
            ++supported;
            require(capability.unavailabilityReason.empty(), "supported forcing has an unavailability reason");
        } else {
            require(!capability.unavailabilityReason.empty(), "unsupported forcing omitted its reason");
        }
    }
    require(supported == 10, "forcing capability matrix must expose eight production pairs and two test pairs");
    require(forcings.capability("WVTestPortableFixedAmplitudeForcing").isSupported(), "registered test forcing pair is unavailable");
    require(forcings.capability("WVTestPortableFixedAmplitudeForcing", 2).status == WVPortableCapabilityStatus::versionMismatch, "forcing pair version mismatch was accepted");
    require(forcings.capability(test::LinearCoefficientForcingIdentifier).isSupported(), "registered linear coefficient forcing pair is unavailable");
    require(forcings.capability("WVUserForcing").status == WVPortableCapabilityStatus::unavailable, "missing forcing pair did not report unavailability");
}

void testRegisteredFixedAmplitudePair() {
    TemporaryFile file(temporaryCopy("forcing-fixed-amplitude.nc"));
    int id = -1;
    requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open paired fixed-amplitude fixture");
    int forcingId = -1;
    requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find paired fixed-amplitude forcing");
    overwriteTextAttribute(forcingId, "AnnotatedClass", "WVTestPortableFixedAmplitudeForcing");
    requireNetCDF(nc_close(id), "close paired fixed-amplitude fixture");
    WVCheckpoint checkpoint;
    const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
    require(static_cast<bool>(result), result.message);
    require(checkpoint.forcingSchedule.entries.size() == 1 &&
                checkpoint.forcingSchedule.entries.front().typeIdentifier == "WVTestPortableFixedAmplitudeForcing",
            "registered fixed-amplitude pair did not reuse the typed payload contract");
}

void testSupportedForcingFixtures() {
    const std::array<std::pair<const char*, const char*>, 6> cases = {{
        {"forcing-nonlinear.nc", "WVNonlinearAdvection"},
        {"forcing-adaptive-damping.nc", "WVAdaptiveDamping"},
        {"forcing-fixed-amplitude.nc", "WVFixedAmplitudeForcing"},
        {"forcing-quadratic-bottom-friction.nc", "WVBottomFrictionQuadratic"},
        {"forcing-pseudo-topographic.nc", "WVPseudoTopographicWaveGeneration"},
        {"forcing-beta-plane.nc", "WVBetaPlanePVAdvection"}
    }};
    for (const auto& testCase : cases) {
        const auto checkpoint = read(testCase.first);
        require(checkpoint.forcingSchedule.entries.size() == 1, std::string(testCase.first) + " did not decode one forcing");
        require(checkpoint.forcingSchedule.entries.front().typeIdentifier == testCase.second, std::string(testCase.first) + " decoded the wrong forcing identity");
    }

    const auto fixed = read("forcing-fixed-amplitude.nc").forcingSchedule.entries.front();
    require(fixed.name == "fixed-amplitude fixture" && fixed.stage == WVForcingStage::spectralAmplitude && fixed.priority == 255, "fixed-amplitude header mismatch");
    require(storedValue<std::vector<std::int64_t>>(fixed,"ApIndices") == std::vector<std::int64_t>({0, 5}) && storedValue<std::vector<std::int64_t>>(fixed,"AmIndices") == std::vector<std::int64_t>({1, 8}) && storedValue<std::vector<std::int64_t>>(fixed,"A0Indices") == std::vector<std::int64_t>({2, 11}), "fixed-amplitude indices were not converted to zero-based offsets");
    require(storedValue<std::vector<double>>(fixed,"ApValuesReal") == std::vector<double>({1.25,-0.5}) && storedValue<std::vector<double>>(fixed,"ApValuesImag")[1] == 0.75, "fixed-amplitude values changed during decoding");

    const auto quadratic = read("forcing-quadratic-bottom-friction.nc").forcingSchedule.entries.front();
    require(storedValue<std::vector<double>>(quadratic,"Cd").front() == 1.7e-3, "quadratic drag coefficient mismatch");

    TemporaryFile linearFile(temporaryCopy("forcing-quadratic-bottom-friction.nc"));
    int linearId = -1;
    requireNetCDF(nc_open(linearFile.path.string().c_str(), NC_WRITE, &linearId), "open linear forcing fixture");
    int linearForcingId = -1;
    requireNetCDF(nc_inq_ncid(linearId, "forcing", &linearForcingId), "find linear forcing group");
    overwriteTextAttribute(linearForcingId, "AnnotatedClass", "WVBottomFrictionLinear");
    overwriteTextAttribute(linearForcingId, "name", "linear bottom friction");
    int rateId = -1;
    requireNetCDF(nc_inq_varid(linearForcingId, "Cd", &rateId), "find source drag variable");
    requireNetCDF(nc_rename_var(linearForcingId, rateId, "r"), "rename linear drag variable");
    const double rate = 2.5e-7;
    requireNetCDF(nc_put_var_double(linearForcingId, rateId, &rate), "write linear drag rate");
    requireNetCDF(nc_close(linearId), "close linear forcing fixture");
    WVCheckpoint linearCheckpoint;
    const auto linearResult = WVCheckpointReader::read(linearFile.path.string(), *test::extensionCatalog(),linearCheckpoint);
    require(static_cast<bool>(linearResult),linearResult.message);
    const auto& linear = linearCheckpoint.forcingSchedule.entries.front();
    require(linear.typeIdentifier == "WVBottomFrictionLinear" && storedValue<std::vector<double>>(linear,"r").front() == rate,"generic forcing persistence did not decode linear drag");

    const auto pseudo = read("forcing-pseudo-topographic.nc").forcingSchedule.entries.front();
    const auto* topography = pseudo.configuration.value("topographicHeight");
    require(pseudo.name == "pseudo-topographic fixture" && topography != nullptr && topography->dimensions == std::vector<std::size_t>({6,8}) && storedValue<std::vector<double>>(pseudo,"topographicHeight").size() == 48, "pseudo-topographic shape mismatch");
    require(storedValue<std::vector<double>>(pseudo,"barotropicVelocityAmplitudeReal")[0] == 0.12 && storedValue<std::vector<double>>(pseudo,"barotropicVelocityAmplitudeImag")[1] == 0.02, "barotropic velocity amplitude mismatch");
    require(storedValue<std::vector<std::string>>(pseudo,"darwinSymbol").front() == "M2" && storedValue<std::vector<double>>(pseudo,"rampDuration").front() == 900.0 && storedValue<std::vector<double>>(pseudo,"startTime").front() == -50.0 && storedValue<std::vector<std::uint8_t>>(pseudo,"shouldAvoidAdaptiveDamping").front(), "pseudo-topographic scalar mismatch");
    require(storedValue<std::vector<double>>(pseudo,"maximumForcedVerticalMode").front() == 2.0, "pseudo-topographic vertical bound mismatch");
}

void testMixedForcingSchedules() {
    for (const std::string file : {"forcing-mixed-hydrostatic.nc", "forcing-mixed-nonhydrostatic.nc"}) {
        const auto checkpoint = read(file);
        const auto& entries = checkpoint.forcingSchedule.entries;
        require(entries.size() == 6, "mixed schedule did not recover all forcing records");
        const std::array<const char*, 6> expected = {"WVNonlinearAdvection", "WVBottomFrictionQuadratic", "WVPseudoTopographicWaveGeneration", "WVAdaptiveDamping", "WVBetaPlanePVAdvection", "WVFixedAmplitudeForcing"};
        for (std::size_t index = 0; index < expected.size(); ++index) {
            require(entries[index].typeIdentifier == expected[index] && entries[index].ordinal == index + 1, "mixed schedule stage/priority/stable ordering mismatch");
        }
        require(entries[2].stage == WVForcingStage::spectral && entries[3].stage == WVForcingStage::spectral && entries[4].stage == WVForcingStage::spectral, "mixed spectral stage mismatch");
    }
}

void overwriteTextAttribute(int groupId, const char* name, const std::string& value) {
    requireNetCDF(nc_redef(groupId), "enter define mode");
    requireNetCDF(nc_put_att_text(groupId, NC_GLOBAL, name, value.size(), value.c_str()), "write text attribute");
    requireNetCDF(nc_enddef(groupId), "leave define mode");
}

void verifyWritableAfterFailure(const std::filesystem::path& path) {
    int id = -1;
    requireNetCDF(nc_open(path.string().c_str(), NC_WRITE, &id), "reopen failed checkpoint");
    requireNetCDF(nc_close(id), "close failed-checkpoint probe");
}

void testUnsupportedVersionAndTransform() {
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open version fixture");
        overwriteTextAttribute(id, "model_version", "5.0.0");
        requireNetCDF(nc_close(id), "close version fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::unsupportedModelVersion, "non-4.x model version was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open transform fixture");
        overwriteTextAttribute(id, "WVTransform", "WVTransformStratifiedQG");
        overwriteTextAttribute(id, "AnnotatedClass", "WVTransformStratifiedQG");
        requireNetCDF(nc_close(id), "close transform fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::unsupportedTransform, "unsupported transform was accepted");
        verifyWritableAfterFailure(file.path);
    }
}

void testMissingPartnerAndWrongType() {
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open partner fixture");
        int variableId = -1;
        requireNetCDF(nc_inq_varid(id, "Ap_imag", &variableId), "find imaginary component");
        requireNetCDF(nc_rename_var(id, variableId, "Ap_imag_missing"), "rename imaginary component");
        requireNetCDF(nc_close(id), "close partner fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::missingComplexPartner, "missing complex partner was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open type fixture");
        requireNetCDF(nc_redef(id), "enter type-fixture define mode");
        int originalId = -1;
        requireNetCDF(nc_inq_varid(id, "Ap_real", &originalId), "find real component");
        requireNetCDF(nc_rename_var(id, originalId, "Ap_real_double"), "rename real component");
        int kl = -1;
        int j = -1;
        requireNetCDF(nc_inq_dimid(id, "kl", &kl), "find kl dimension");
        requireNetCDF(nc_inq_dimid(id, "j", &j), "find j dimension");
        const int dimensions[] = {kl, j};
        int replacement = -1;
        requireNetCDF(nc_def_var(id, "Ap_real", NC_FLOAT, 2, dimensions, &replacement), "define wrong-type component");
        unsigned char one = 1;
        unsigned char zero = 0;
        requireNetCDF(nc_put_att_uchar(id, replacement, "isComplex", NC_UBYTE, 1, &one), "mark wrong-type component complex");
        requireNetCDF(nc_put_att_uchar(id, replacement, "isRealPart", NC_UBYTE, 1, &one), "mark wrong-type component real");
        requireNetCDF(nc_put_att_uchar(id, replacement, "isImaginaryPart", NC_UBYTE, 1, &zero), "mark wrong-type component nonimaginary");
        requireNetCDF(nc_enddef(id), "leave type-fixture define mode");
        requireNetCDF(nc_close(id), "close type fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::typeMismatch, "wrong coefficient type was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open dimension-order fixture");
        requireNetCDF(nc_redef(id), "enter dimension-order define mode");
        int originalId = -1;
        requireNetCDF(nc_inq_varid(id, "Am_real", &originalId), "find ordered real component");
        requireNetCDF(nc_rename_var(id, originalId, "Am_real_ordered"), "rename ordered real component");
        int kl = -1;
        int j = -1;
        requireNetCDF(nc_inq_dimid(id, "kl", &kl), "find ordered kl dimension");
        requireNetCDF(nc_inq_dimid(id, "j", &j), "find ordered j dimension");
        const int reversedDimensions[] = {j, kl};
        int replacement = -1;
        requireNetCDF(nc_def_var(id, "Am_real", NC_DOUBLE, 2, reversedDimensions, &replacement), "define reversed component");
        unsigned char one = 1;
        unsigned char zero = 0;
        requireNetCDF(nc_put_att_uchar(id, replacement, "isComplex", NC_UBYTE, 1, &one), "mark reversed component complex");
        requireNetCDF(nc_put_att_uchar(id, replacement, "isRealPart", NC_UBYTE, 1, &one), "mark reversed component real");
        requireNetCDF(nc_put_att_uchar(id, replacement, "isImaginaryPart", NC_UBYTE, 1, &zero), "mark reversed component nonimaginary");
        requireNetCDF(nc_enddef(id), "leave dimension-order define mode");
        requireNetCDF(nc_close(id), "close dimension-order fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::shapeMismatch, "wrong coefficient dimension order was accepted");
        verifyWritableAfterFailure(file.path);
    }
}

void addComplexVariableDefinition(int groupId, const std::string& baseName, const int* dimensions) {
    for (const bool imaginary : {false, true}) {
        const std::string name = baseName + (imaginary ? "_imag" : "_real");
        int variableId = -1;
        requireNetCDF(nc_def_var(groupId, name.c_str(), NC_DOUBLE, 2, dimensions, &variableId), "define duplicate coefficient");
        unsigned char one = 1;
        unsigned char zero = 0;
        requireNetCDF(nc_put_att_uchar(groupId, variableId, "isComplex", NC_UBYTE, 1, &one), "mark duplicate complex");
        requireNetCDF(nc_put_att_uchar(groupId, variableId, "isRealPart", NC_UBYTE, 1, imaginary ? &zero : &one), "mark duplicate real part");
        requireNetCDF(nc_put_att_uchar(groupId, variableId, "isImaginaryPart", NC_UBYTE, 1, imaginary ? &one : &zero), "mark duplicate imaginary part");
    }
}

void testAmbiguousStateAndIndex() {
    {
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(fixture("time-series-nonhydrostatic.nc").string(), *test::extensionCatalog(), checkpoint, WVCheckpointStateSelection::atIndex(3));
        require(result.code == WVCheckpointStatusCode::stateIndexOutOfRange, "out-of-range state index was accepted");
    }
    {
        TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open duplicate-state fixture");
        requireNetCDF(nc_redef(id), "enter duplicate-state define mode");
        int groupId = -1;
        requireNetCDF(nc_def_grp(id, "duplicate-state", &groupId), "define duplicate state group");
        int kl = -1;
        int j = -1;
        requireNetCDF(nc_inq_dimid(groupId, "kl", &kl), "find inherited kl dimension");
        requireNetCDF(nc_inq_dimid(groupId, "j", &j), "find inherited j dimension");
        int timeId = -1;
        requireNetCDF(nc_def_var(groupId, "t", NC_DOUBLE, 0, nullptr, &timeId), "define duplicate time");
        const int dimensions[] = {kl, j};
        for (const std::string name : {"Ap", "Am", "A0"}) addComplexVariableDefinition(groupId, name, dimensions);
        requireNetCDF(nc_enddef(id), "leave duplicate-state define mode");
        requireNetCDF(nc_close(id), "close duplicate-state fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::ambiguousState, "duplicate state group was accepted");
        verifyWritableAfterFailure(file.path);
    }
}

void testInvalidConfiguration() {
    TemporaryFile file(temporaryCopy("root-hydrostatic.nc"));
    int id = -1;
    requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open invalid-configuration fixture");
    int variableId = -1;
    requireNetCDF(nc_inq_varid(id, "N0", &variableId), "find N0");
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    requireNetCDF(nc_put_var_double(id, variableId, &invalid), "write invalid N0");
    requireNetCDF(nc_close(id), "close invalid-configuration fixture");
    WVCheckpoint checkpoint;
    const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
    require(result.code == WVCheckpointStatusCode::invalidValue, "non-finite configuration was accepted");
    verifyWritableAfterFailure(file.path);
}

void testOrderedForcingHeaders() {
    TemporaryFile file(temporaryCopy("root-nonhydrostatic.nc"));
    int id = -1;
    requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open forcing fixture");
    requireNetCDF(nc_redef(id), "enter forcing define mode");
    int forcingId = -1;
    requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find forcing group");
    requireNetCDF(nc_del_att(forcingId, NC_GLOBAL, "AnnotatedClass"), "remove singleton forcing tag");
    for (const auto& record : std::array<std::pair<const char*, const char*>, 2>{
            std::pair<const char*, const char*>{"forcing-2", "WVAdaptiveDamping"},
            {"forcing-1", "WVNonlinearAdvection"}}) {
        int childId = -1;
        requireNetCDF(nc_def_grp(forcingId, record.first, &childId), "define forcing record");
        requireNetCDF(nc_put_att_text(childId, NC_GLOBAL, "AnnotatedClass", std::char_traits<char>::length(record.second), record.second), "write forcing class");
    }
    requireNetCDF(nc_enddef(id), "leave forcing define mode");
    requireNetCDF(nc_close(id), "close forcing fixture");

    WVCheckpoint checkpoint;
    const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
    require(static_cast<bool>(result), result.message);
    require(checkpoint.metadata.forcingHeaders.size() == 2, "forcing array was not recovered");
    require(checkpoint.metadata.forcingHeaders[0].ordinal == 1 && checkpoint.metadata.forcingHeaders[0].groupPath == "/forcing/forcing-1" && checkpoint.metadata.forcingHeaders[0].annotatedClass == "WVNonlinearAdvection", "first forcing record was not ordered");
    require(checkpoint.metadata.forcingHeaders[1].ordinal == 2 && checkpoint.metadata.forcingHeaders[1].groupPath == "/forcing/forcing-2" && checkpoint.metadata.forcingHeaders[1].annotatedClass == "WVAdaptiveDamping", "second forcing record was not ordered");
    require(checkpoint.forcingSchedule.entries.size() == 2 && checkpoint.forcingSchedule.entries[0].typeIdentifier == "WVNonlinearAdvection" && checkpoint.forcingSchedule.entries[1].typeIdentifier == "WVAdaptiveDamping", "forcing schedule was not ordered by stage and priority");
}

void testUnsupportedForcingClasses() {
    const std::array<const char*, 6> unsupported = {"WVAntialiasing", "WVHorizontalDamping", "WVVerticalDamping", "WVThermalDamping", "WVVerticalDiffusivity", "WVUserForcing"};
    for (const char* typeIdentifier : unsupported) {
        TemporaryFile file(temporaryCopy("forcing-nonlinear.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open unsupported forcing fixture");
        int forcingId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find singleton forcing group");
        overwriteTextAttribute(forcingId, "AnnotatedClass", typeIdentifier);
        requireNetCDF(nc_close(id), "close unsupported forcing fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::unsupportedForcing, std::string(typeIdentifier) + " did not report unsupportedForcing");
        verifyWritableAfterFailure(file.path);
    }
}

void testMalformedForcingRecords() {
    {
        TemporaryFile file(temporaryCopy("forcing-quadratic-bottom-friction.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open invalid r fixture");
        int forcingId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find linear forcing group");
        overwriteTextAttribute(forcingId, "AnnotatedClass", "WVBottomFrictionLinear");
        int variableId = -1;
        requireNetCDF(nc_inq_varid(forcingId, "Cd", &variableId), "find source r variable");
        requireNetCDF(nc_rename_var(forcingId, variableId, "r"), "rename invalid r variable");
        const double invalid = -1.0;
        requireNetCDF(nc_put_var_double(forcingId, variableId, &invalid), "write invalid r");
        requireNetCDF(nc_close(id), "close invalid r fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::malformedForcing, "negative r was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("forcing-quadratic-bottom-friction.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open invalid Cd fixture");
        int forcingId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find quadratic forcing group");
        int variableId = -1;
        requireNetCDF(nc_inq_varid(forcingId, "Cd", &variableId), "find Cd");
        const double invalid = -1.0;
        requireNetCDF(nc_put_var_double(forcingId, variableId, &invalid), "write invalid Cd");
        requireNetCDF(nc_close(id), "close invalid Cd fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::malformedForcing, "negative Cd was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("forcing-fixed-amplitude.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open invalid fixed-index fixture");
        int forcingId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find fixed forcing group");
        int variableId = -1;
        requireNetCDF(nc_inq_varid(forcingId, "Ap_indices", &variableId), "find fixed indices");
        const unsigned long long invalid[] = {0, 6};
        requireNetCDF(nc_put_var_ulonglong(forcingId, variableId, invalid), "write invalid fixed indices");
        requireNetCDF(nc_close(id), "close invalid fixed-index fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::incompatibleForcing, "zero fixed-amplitude index was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("forcing-fixed-amplitude.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open duplicate fixed-index fixture");
        int forcingId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find fixed forcing group");
        int variableId = -1;
        requireNetCDF(nc_inq_varid(forcingId, "Ap_indices", &variableId), "find fixed indices");
        const unsigned long long duplicate[] = {1, 1};
        requireNetCDF(nc_put_var_ulonglong(forcingId, variableId, duplicate), "write duplicate fixed indices");
        requireNetCDF(nc_close(id), "close duplicate fixed-index fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::duplicateForcing, "duplicate fixed-amplitude index was accepted");
        verifyWritableAfterFailure(file.path);
    }
    {
        TemporaryFile file(temporaryCopy("forcing-mixed-nonhydrostatic.nc"));
        int id = -1;
        requireNetCDF(nc_open(file.path.string().c_str(), NC_WRITE, &id), "open duplicate-name fixture");
        int forcingId = -1;
        int pseudoId = -1;
        requireNetCDF(nc_inq_ncid(id, "forcing", &forcingId), "find forcing container");
        requireNetCDF(nc_inq_ncid(forcingId, "forcing-3", &pseudoId), "find pseudo forcing");
        overwriteTextAttribute(pseudoId, "name", "fixed-amplitude fixture");
        requireNetCDF(nc_close(id), "close duplicate-name fixture");
        WVCheckpoint checkpoint;
        const auto result = WVCheckpointReader::read(file.path.string(), *test::extensionCatalog(), checkpoint);
        require(result.code == WVCheckpointStatusCode::duplicateForcing, "duplicate forcing name was accepted");
        verifyWritableAfterFailure(file.path);
    }
}

} // namespace

int main() {
    try {
        testPositiveFixtures();
        testAllocationLightInspection();
        testForcingCapabilities();
        testRegisteredFixedAmplitudePair();
        testSupportedForcingFixtures();
        testMixedForcingSchedules();
        testUnsupportedVersionAndTransform();
        testMissingPartnerAndWrongType();
        testAmbiguousStateAndIndex();
        testInvalidConfiguration();
        testOrderedForcingHeaders();
        testUnsupportedForcingClasses();
        testMalformedForcingRecords();
        std::cout << "WaveVortex checkpoint reader tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
