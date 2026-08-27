#pragma once

#include <cstdint>

namespace wavevortex::runtime::detail {

enum class WVCheckpointWriterFailurePoint : std::uint8_t {
    none,
    afterDefinition,
    afterWrite,
    beforeCommit
};

void setCheckpointWriterFailurePoint(WVCheckpointWriterFailurePoint point) noexcept;

} // namespace wavevortex::runtime::detail
