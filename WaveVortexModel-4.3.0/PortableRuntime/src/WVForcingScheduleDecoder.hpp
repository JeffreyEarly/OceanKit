#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace wavevortex::runtime::detail {

struct WVForcingGroupSource {
    int groupId = -1;
    std::size_t ordinal = 0;
    std::string groupPath;
    std::string annotatedClass;
};

WVCheckpointStatus decodeForcingSchedule(
    const std::vector<WVForcingGroupSource>& sources,
    const WVTransformConstantStratificationConfiguration& configuration,
    std::size_t coefficientCount,
    const WVExtensionCatalog& catalog,
    WVFrozenForcingSchedule& schedule);
WVCheckpointStatus decodeForcingSchedule(
    const std::vector<WVForcingGroupSource>& sources,
    const WVTransformBarotropicQGConfiguration& configuration,
    std::size_t coefficientCount,
    const WVExtensionCatalog& catalog,
    WVFrozenForcingSchedule& schedule);

} // namespace wavevortex::runtime::detail
