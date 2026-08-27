#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <string>

namespace wavevortex::runtime {

enum class WVCheckpointCommitPolicy {
    replaceExisting,
    createNew
};

// Write one compatible root-level WaveVortexModel 4.x checkpoint.
//
// The writer validates and writes a same-directory temporary file, reads that
// file through WVCheckpointReader, and commits it only after every operation
// succeeds. The default policy atomically replaces a destination; createNew
// atomically refuses replacement. On failure, an existing destination is
// unchanged and the temporary file is removed.
class WVCheckpointWriter final {
public:
    static WVCheckpointStatus write(
        const std::string& path,
        const WVExtensionCatalog& catalog,
        const WVCheckpoint& checkpoint,
        WVCheckpointCommitPolicy commitPolicy = WVCheckpointCommitPolicy::replaceExisting);
};

} // namespace wavevortex::runtime
