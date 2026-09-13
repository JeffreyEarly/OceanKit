#pragma once

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVExtensionCatalogBuilder;

namespace detail {

enum class WVMovingFieldChannel : std::uint8_t {
  x,
  y,
  z,
  tracerValue
};

WVKernelStatus addBuiltInObserverFactories(
    wavevortex::runtime::WVExtensionCatalogBuilder &builder);
WVKernelStatus canonicalCoefficientObserver(
    std::string identifier, const WVExtensionCatalog &catalog,
    WVObserverRecord &observer,
    std::vector<std::string> coefficientFamilies = {"Ap", "Am", "A0"});

const char *movingFieldChannelName(WVMovingFieldChannel channel) noexcept;
std::vector<WVMovingFieldChannel>
particlePositionChannels(bool isXYOnly);
std::string movingFieldVariableName(const WVObserverRecord &observer,
                                    WVMovingFieldChannel channel);

} // namespace detail
} // namespace wavevortex::runtime
