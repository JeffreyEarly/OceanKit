#pragma once

#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WaveVortexRuntime/WVObserverContracts.hpp"
#include "WaveVortexRuntime/WVOutputSchedule.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

namespace detail {
class WVObserverFactoryRegistrationAccess;
}

using WVObserverFactory = std::function<WVKernelStatus(
    const WVObserverRecord &, const WVPortableTypedRecord &,
    std::shared_ptr<const WVObservingSystem> &)>;
using WVObserverConfigurationResolver = std::function<WVKernelStatus(
    const WVObserverRecord &, WVPortableTypedRecord &)>;
using WVObserverOutputPlanResolver = std::function<WVKernelStatus(
    const WVObserverRecord &, const WVObserverOutputPlanningContext &,
    WVObserverOutputPlan &)>;

struct WVObserverFactoryRegistration {
  // Stable source API v1 has exactly these five public constructor inputs:
  // identity, version, factory, optional data-only configuration resolver, and
  // optional data-only output-plan resolver. Legacy compatibility callbacks
  // and persistence metadata are intentionally not extension inputs.
  WVObserverFactoryRegistration(
      std::string identity, std::uint32_t version, WVObserverFactory make,
      WVObserverConfigurationResolver resolve = {},
      WVObserverOutputPlanResolver resolveOutputPlan = {})
      : typeIdentifier(std::move(identity)), contractVersion(version),
        factory(std::move(make)), configurationResolver(std::move(resolve)),
        outputPlanResolver(std::move(resolveOutputPlan)) {}
  std::string typeIdentifier;
  std::uint32_t contractVersion = WVPortablePairContractVersion;
  WVObserverFactory factory;
  // Data-only compatibility normalization. This callback may reconstruct a
  // typed configuration from legacy record fields, but must not instantiate
  // runtime observer behavior.
  WVObserverConfigurationResolver configurationResolver;
  // Data-only semantic preflight. This must derive the exact output plan
  // without constructing an observing-system implementation.
  WVObserverOutputPlanResolver outputPlanResolver;

private:
  friend class detail::WVObserverFactoryRegistrationAccess;
  // Private ownership for built-in legacy-encoding compatibility. Extensions
  // cannot set or depend on this implementation detail.
  std::shared_ptr<const void> internalCompatibility_;
};

struct WVOutputScheduleFactoryRegistration {
  std::string typeIdentifier;
  std::uint32_t contractVersion = 1;
  WVOutputScheduleFactory factory = nullptr;
};

class WVObserverCatalog final {
public:
  const std::vector<WVObserverFactoryRegistration> &registrations() const
      noexcept {
    return registrations_;
  }
  const WVObserverFactoryRegistration *
  registration(const std::string &typeIdentifier,
               std::uint32_t contractVersion) const noexcept;
  WVPortableCapability capability(
      std::string typeIdentifier,
      std::uint32_t contractVersion = WVPortablePairContractVersion) const;
  WVKernelStatus resolveConfiguration(
      const WVObserverRecord &record,
      WVPortableTypedRecord &configuration) const;
  WVKernelStatus resolveOutputPlan(
      const WVObserverRecord &record,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const;
  WVKernelStatus create(const WVObserverRecord &record,
                        const WVPortableTypedRecord &configuration,
                        std::shared_ptr<const WVObservingSystem> &result) const;
  std::size_t persistentBytes() const noexcept;

private:
  friend class WVExtensionCatalogBuilder;
  std::vector<WVObserverFactoryRegistration> registrations_;
};

class WVOutputScheduleCatalog final {
public:
  const std::vector<WVOutputScheduleFactoryRegistration> &registrations() const
      noexcept {
    return registrations_;
  }
  const WVOutputScheduleFactoryRegistration *
  registration(const std::string &typeIdentifier,
               std::uint32_t contractVersion) const noexcept;
  const WVOutputScheduleFactoryRegistration *
  registration(const WVOutputScheduleRecord &record) const noexcept;
  WVKernelStatus resolve(const WVOutputScheduleRecord &record,
                         std::shared_ptr<const WVOutputSchedule> &result) const;
  std::size_t persistentBytes() const noexcept;

private:
  friend class WVExtensionCatalogBuilder;
  std::vector<WVOutputScheduleFactoryRegistration> registrations_;
};

class WVForcingCatalog final {
public:
  using Registration = WVForcingFactoryRegistration;

  const std::vector<Registration> &registrations() const noexcept {
    return registrations_;
  }
  const Registration *
  registration(const std::string &typeIdentifier,
               std::uint32_t contractVersion) const noexcept;
  WVPortableCapability capability(
      std::string typeIdentifier,
      std::uint32_t contractVersion = WVPortablePairContractVersion) const;
  WVKernelStatus create(
      const WVFrozenForcingEntry &entry,
      const WVTransformConstantStratificationDescriptor &descriptor,
      bool hasAdaptiveDamping, std::unique_ptr<WVForcing> &forcing) const;
  WVKernelStatus createBarotropicQG(
      const WVFrozenForcingEntry &entry,
      const WVTransformBarotropicQGDescriptor &descriptor,
      bool hasAdaptiveDamping,
      std::unique_ptr<WVBarotropicQGForcing> &forcing) const;
  WVKernelStatus
  validateConfiguration(const WVFrozenForcingEntry &entry) const;
  std::size_t persistentBytes() const noexcept;

private:
  friend class WVExtensionCatalogBuilder;
  std::vector<Registration> registrations_;
};

class WVExtensionCatalog final {
public:
  const WVObserverCatalog &observers() const noexcept { return observers_; }
  const WVOutputScheduleCatalog &outputSchedules() const noexcept {
    return outputSchedules_;
  }
  const WVForcingCatalog &forcings() const noexcept { return forcings_; }
  std::size_t persistentBytes() const noexcept;

private:
  friend class WVExtensionCatalogBuilder;
  WVObserverCatalog observers_;
  WVOutputScheduleCatalog outputSchedules_;
  WVForcingCatalog forcings_;
};

// The only mutable source API v1 registration boundary. Add all registrations
// before freeze(); any rejected registration invalidates the builder. A
// builder freezes exactly once into an immutable shared catalog. Runtime
// owners retain that shared catalog, so the builder and caller's original
// handle need not outlive a resolved model, output configuration, or runner.
class WVExtensionCatalogBuilder final {
public:
  WVKernelStatus
  addObserverFactory(WVObserverFactoryRegistration registration);
  WVKernelStatus addOutputScheduleFactory(
      WVOutputScheduleFactoryRegistration registration);
  WVKernelStatus addForcingFactory(WVForcingCatalog::Registration registration);
  WVKernelStatus freeze(
      std::shared_ptr<const WVExtensionCatalog> &catalog);
  bool isFrozen() const noexcept { return frozen_; }

private:
  bool frozen_ = false;
  bool valid_ = true;
  std::string validationError_;
  std::vector<WVObserverFactoryRegistration> observers_;
  std::vector<WVOutputScheduleFactoryRegistration> outputSchedules_;
  std::vector<WVForcingCatalog::Registration> forcings_;
};

WVKernelStatus addBuiltInExtensions(WVExtensionCatalogBuilder &builder);
WVKernelStatus makeBuiltInExtensionCatalog(
    std::shared_ptr<const WVExtensionCatalog> &catalog);

} // namespace wavevortex::runtime
