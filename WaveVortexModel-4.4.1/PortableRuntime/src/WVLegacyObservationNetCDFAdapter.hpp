#pragma once

#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace wavevortex::runtime::detail {

// Compatibility boundary for the five legacy MATLAB observing-system
// encodings. The graph reader remains schema-driven; this adapter contains the
// data-only mapping between legacy class metadata and portable records.
WVCheckpointStatus parsePersistedObserver(
    int outputGroup, int metadataGroup, const std::string &outputPath,
    const WVExtensionCatalog &catalog,
    WVPortableObserverRecord &portable, std::string &identifier);

bool persistedObserverCarriesCoefficientState(
    int metadataGroup, const WVExtensionCatalog &catalog) noexcept;

// Historical MATLAB releases used a small number of presentation-only
// attribute spellings that differ from the portable writer's byte-compatible
// legacy encoding. Accept those spellings when appending without changing the
// attributes emitted for newly created files.
bool legacyObservationAttributeMatches(std::string_view name,
                                       std::string_view expected,
                                       std::string_view observed) noexcept;

// Validate availability and cross-group equality of dynamic restart state
// with bounded scratch storage. Only fixed configuration coordinates already
// retained by the observer record are materialized.
WVCheckpointStatus inspectPersistedObserverRestartState(
    const std::vector<WVOutputFileRecord> &files,
    std::vector<WVObserverRecord> &observers, double selectedTime,
    const WVExtensionCatalog &catalog);

WVCheckpointStatus resolvePersistedObserverRestartState(
    const std::vector<WVOutputFileRecord> &files,
    std::vector<WVObserverRecord> &observers, double selectedTime,
    const WVExtensionCatalog &catalog,
    std::map<std::string, std::vector<std::vector<double>>> &resolvedState);

} // namespace wavevortex::runtime::detail
