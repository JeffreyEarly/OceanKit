#pragma once

#include "WaveVortexRuntime/WVForcingSchedule.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVExtensionCatalog;

// Version of the reader contract. This identifier is not written into 4.x files.
inline constexpr std::uint32_t WVCheckpointProfileVersion = 1;
inline constexpr const char* WVCheckpointProfileIdentifier = "wave-vortex-4x-v1";

enum class WVCheckpointStatusCode : std::uint32_t {
    success = 0,
    openFailure,
    netcdfFailure,
    missingAttribute,
    missingDimension,
    missingVariable,
    missingComplexPartner,
    unsupportedModelVersion,
    unsupportedTransform,
    typeMismatch,
    shapeMismatch,
    invalidValue,
    ambiguousState,
    stateIndexOutOfRange,
    descriptorFailure,
    unsupportedForcing,
    malformedForcing,
    duplicateForcing,
    incompatibleForcing,
    writeFailure,
    commitFailure,
    schemaMismatch,
    incompleteRecord,
    appendConflict,
    unsupportedObserver
};

struct WVCheckpointStatus {
    WVCheckpointStatusCode code = WVCheckpointStatusCode::success;
    std::string message;
    std::string location;

    explicit operator bool() const noexcept { return code == WVCheckpointStatusCode::success; }
    static WVCheckpointStatus ok() { return {}; }
};

enum class WVCheckpointStateSelectionKind : std::uint8_t {
    latest,
    index
};

struct WVCheckpointStateSelection {
    WVCheckpointStateSelectionKind kind = WVCheckpointStateSelectionKind::latest;
    std::size_t index = 0;

    static WVCheckpointStateSelection latest() noexcept { return {}; }
    static WVCheckpointStateSelection atIndex(std::size_t value) noexcept {
        return {WVCheckpointStateSelectionKind::index, value};
    }
};

// Owned canonical coefficient storage. The first dimension, Nj, is adjacent.
struct WVCheckpointCoefficients {
    WVShape2D shape;
    std::vector<WVComplex64> Ap;
    std::vector<WVComplex64> Am;
    std::vector<WVComplex64> A0;

    WVCoefficients view() const noexcept;
};

// One selected checkpoint state and its coefficient reference time.
struct WVCheckpointState {
    double t = 0.0;
    double t0 = 0.0;
    WVCheckpointCoefficients coefficients;

    WVState view() const noexcept;
};

// Identity of one existing annotated forcing record. Parameter decoding is separate.
struct WVCheckpointForcingHeader {
    std::size_t ordinal = 0;
    std::string groupPath;
    std::string annotatedClass;
};

// Source and selection details retained for diagnostics and compatibility decisions.
struct WVCheckpointMetadata {
    std::string profileIdentifier = WVCheckpointProfileIdentifier;
    std::uint32_t profileVersion = WVCheckpointProfileVersion;
    std::string modelVersion;
    std::string transformClass;
    std::string stateGroupPath;
    std::size_t selectedStateIndex = 0;
    std::size_t stateCount = 0;
    std::vector<WVCheckpointForcingHeader> forcingHeaders;
};

enum class WVPersistedTransformKind : std::uint8_t {
    constantStratification,
    barotropicQG
};

// Complete owning result needed to rebuild the portable constant-stratification core.
struct WVCheckpoint {
    WVPersistedTransformKind transformKind =
        WVPersistedTransformKind::constantStratification;
    WVTransformConstantStratificationConfiguration configuration;
    WVTransformBarotropicQGConfiguration barotropicQGConfiguration;
    WVTransformStateDescription stateDescription;
    WVCheckpointState state;
    // Populated only for transforms that do not use the stabilized legacy
    // Ap/Am/A0 checkpoint representation. Barotropic QG therefore owns one
    // compact A0 family here while the three legacy vectors above stay empty.
    WVTransformStateCheckpoint transformState;
    WVCheckpointMetadata metadata;
    WVFrozenForcingSchedule forcingSchedule;
};

// Exact capacity owned by the transform's persisted coefficient families.
// Transform-specific storage selection is implemented by the named
// checkpoint adapter rather than by model or CLI orchestration.
std::size_t checkpointCoefficientStorageBytes(
    const WVCheckpoint& checkpoint) noexcept;

// Capacity-based retained storage owned by one complete checkpoint, including
// its object, coefficient arrays, metadata, and frozen forcing records.
inline std::size_t
checkpointRetainedBytes(const WVCheckpoint& checkpoint) noexcept {
    std::size_t bytes = sizeof(checkpoint) +
        (checkpoint.state.coefficients.Ap.capacity() +
         checkpoint.state.coefficients.Am.capacity() +
         checkpoint.state.coefficients.A0.capacity()) * sizeof(WVComplex64) +
        checkpoint.stateDescription.transformIdentifier.capacity() +
        checkpoint.stateDescription.spatialDimensions.capacity() *
            sizeof(std::size_t) +
        checkpoint.stateDescription.coefficientFamilies.capacity() *
            sizeof(WVCoefficientFamilyDescription) +
        checkpoint.transformState.persistentBytes() -
            sizeof(WVTransformStateCheckpoint) +
        checkpoint.metadata.profileIdentifier.capacity() +
        checkpoint.metadata.modelVersion.capacity() +
        checkpoint.metadata.transformClass.capacity() +
        checkpoint.metadata.stateGroupPath.capacity() +
        checkpoint.metadata.forcingHeaders.capacity() *
            sizeof(WVCheckpointForcingHeader) +
        checkpoint.forcingSchedule.profileIdentifier.capacity() +
        checkpoint.forcingSchedule.entries.capacity() *
            sizeof(WVFrozenForcingEntry);
    for (const auto& header : checkpoint.metadata.forcingHeaders) {
        bytes += header.groupPath.capacity() + header.annotatedClass.capacity();
    }
    for (const auto& family : checkpoint.stateDescription.coefficientFamilies) {
        bytes += family.identifier.capacity() +
            family.spectralDimensions.capacity() * sizeof(std::size_t);
    }
    for (const auto& entry : checkpoint.forcingSchedule.entries) {
        bytes += entry.typeIdentifier.capacity() + entry.name.capacity() +
            entry.sourceGroupPath.capacity() +
            entry.configuration.persistentBytes() -
                sizeof(WVPortableTypedRecord);
    }
    return bytes;
}

// Allocation-light checkpoint information used to reject incompatible input
// before loading the three state-sized coefficient arrays.
struct WVCheckpointInspection {
    WVPersistedTransformKind transformKind =
        WVPersistedTransformKind::constantStratification;
    WVTransformConstantStratificationConfiguration configuration;
    WVTransformBarotropicQGConfiguration barotropicQGConfiguration;
    // Resolved before read() allocates or loads any coefficient array. The
    // persistence adapter remains responsible for its transform's encoding.
    WVTransformStateDescription stateDescription;
    WVShape2D coefficientShape;
    double t = 0.0;
    double t0 = 0.0;
    WVCheckpointMetadata metadata;
    WVFrozenForcingSchedule forcingSchedule;
};

// Read and structurally validate existing WaveVortexModel 4.x NetCDF checkpoints.
class WVCheckpointReader final {
public:
    // Inspect structure, configuration, the selected state shape, and forcing
    // without allocating or reading Ap, Am, or A0.
    static WVCheckpointStatus inspect(
        const std::string& path,
        const WVExtensionCatalog& catalog,
        WVCheckpointInspection& inspection,
        WVCheckpointStateSelection selection = WVCheckpointStateSelection::latest());

    // On failure, checkpoint is unchanged and all NetCDF handles are closed.
    static WVCheckpointStatus read(
        const std::string& path,
        const WVExtensionCatalog& catalog,
        WVCheckpoint& checkpoint,
        WVCheckpointStateSelection selection = WVCheckpointStateSelection::latest());
};

} // namespace wavevortex::runtime
