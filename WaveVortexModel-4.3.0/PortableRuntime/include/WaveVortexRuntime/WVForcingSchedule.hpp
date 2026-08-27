#pragma once

#include "WaveVortexRuntime/WVPortableImplementationContract.hpp"
#include "WaveVortexRuntime/WVPortableTypedRecord.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {

inline constexpr std::uint32_t WVForcingScheduleProfileVersion = 1;
inline constexpr const char *WVForcingScheduleProfileIdentifier =
    "wave-vortex-forcing-v1";

enum class WVForcingStage : std::uint8_t {
  spatial,
  spectral,
  spectralAmplitude
};

inline const char *forcingStageName(WVForcingStage stage) noexcept {
  switch (stage) {
  case WVForcingStage::spatial:
    return "spatial";
  case WVForcingStage::spectral:
    return "spectral";
  case WVForcingStage::spectralAmplitude:
    return "spectral-amplitude";
  }
  return "unknown";
}

// Construction/persistence record. The named configuration is decoded once by
// the source-linked forcing factory; integration retains the resolved instance.
struct WVFrozenForcingEntry {
  std::string typeIdentifier;
  std::uint32_t contractVersion = WVPortablePairContractVersion;
  std::string name;
  WVForcingStage stage = WVForcingStage::spatial;
  std::uint8_t priority = 255;
  std::size_t ordinal = 0;
  std::string sourceGroupPath;
  WVPortableTypedRecord configuration;
};

struct WVFrozenForcingSchedule {
  std::string profileIdentifier = WVForcingScheduleProfileIdentifier;
  std::uint32_t profileVersion = WVForcingScheduleProfileVersion;
  std::vector<WVFrozenForcingEntry> entries;
};

} // namespace wavevortex::runtime
