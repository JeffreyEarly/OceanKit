#pragma once

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <functional>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime::detail {

// Data-only compatibility actions for the five historical MATLAB observer
// encodings. This is runtime-owned migration machinery, not part of the stable
// source-linked extension API.
struct WVLegacyObserverOperationBinder {
  std::function<WVKernelStatus()> fullField;
  std::function<WVKernelStatus()> fixedVerticalProfiles;
  std::function<WVKernelStatus()> fixedPositions;
  std::function<WVKernelStatus()> movingPositions;
  std::function<WVKernelStatus()> integratedState;
};

using WVLegacyObserverOperationResolver = std::function<WVKernelStatus(
    const WVObserverRecord &, const WVLegacyObserverOperationBinder &)>;

struct WVLegacyObserverPersistenceMetadata {
  std::string fieldListAttribute;
  std::vector<std::string> coefficientRestartFamilies;
  std::string defaultIdentifier;
  bool appendFieldsToDefaultIdentifier = false;
};

struct WVLegacyObserverCompatibility {
  WVLegacyObserverOperationResolver operationResolver;
  WVLegacyObserverPersistenceMetadata persistence;
};

class WVObserverFactoryRegistrationAccess final {
public:
  static WVKernelStatus attachLegacyCompatibility(
      WVObserverFactoryRegistration &registration,
      WVLegacyObserverOperationResolver resolver,
      WVLegacyObserverPersistenceMetadata persistence) {
    try {
      registration.internalCompatibility_ =
          std::make_shared<const WVLegacyObserverCompatibility>(
              WVLegacyObserverCompatibility{std::move(resolver),
                                            std::move(persistence)});
      return WVKernelStatus::ok();
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Legacy observer compatibility allocation failed."};
    }
  }

  static const WVLegacyObserverCompatibility *
  legacyCompatibility(
      const WVObserverFactoryRegistration &registration) noexcept {
    return static_cast<const WVLegacyObserverCompatibility *>(
        registration.internalCompatibility_.get());
  }

  static std::size_t persistentBytes(
      const WVObserverFactoryRegistration &registration) noexcept {
    const auto *compatibility = legacyCompatibility(registration);
    if (compatibility == nullptr)
      return 0;
    std::size_t bytes = sizeof(*compatibility) +
                        compatibility->persistence.fieldListAttribute.capacity() +
                        compatibility->persistence.defaultIdentifier.capacity() +
                        compatibility->persistence.coefficientRestartFamilies.capacity() *
                            sizeof(std::string);
    for (const auto &family :
         compatibility->persistence.coefficientRestartFamilies)
      bytes += family.capacity();
    return bytes;
  }
};

} // namespace wavevortex::runtime::detail
