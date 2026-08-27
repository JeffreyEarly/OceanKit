#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"
#include "WVLegacyObserverCompatibility.hpp"
#include "WVObserverAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <optional>
#include <set>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

template <typename Registration, typename Identity>
bool containsIdentity(const std::vector<Registration> &registrations,
                      const Identity &identity) {
  return std::any_of(registrations.begin(), registrations.end(),
                     [&](const auto &registration) {
                       return identity(registration);
                     });
}

} // namespace

const WVObserverFactoryRegistration *WVObserverCatalog::registration(
    const std::string &typeIdentifier,
    std::uint32_t contractVersion) const noexcept {
  const auto found = std::find_if(
      registrations_.begin(), registrations_.end(), [&](const auto &value) {
        return value.typeIdentifier == typeIdentifier &&
               value.contractVersion == contractVersion;
      });
  return found == registrations_.end() ? nullptr : &*found;
}

WVPortableCapability WVObserverCatalog::capability(
    std::string typeIdentifier, std::uint32_t contractVersion) const {
  std::optional<WVPortableImplementationIdentity> available;
  for (const auto &value : registrations_)
    if (value.typeIdentifier == typeIdentifier &&
        value.contractVersion == contractVersion) {
      available = WVPortableImplementationIdentity{value.typeIdentifier,
                                                    value.contractVersion};
      break;
    }
  if (!available)
    for (const auto &value : registrations_)
      if (value.typeIdentifier == typeIdentifier) {
        available = WVPortableImplementationIdentity{value.typeIdentifier,
                                                      value.contractVersion};
        break;
      }
  return evaluatePortableCapability(
      {std::move(typeIdentifier), contractVersion}, std::move(available));
}

WVKernelStatus WVObserverCatalog::resolveConfiguration(
    const WVObserverRecord &record,
    WVPortableTypedRecord &configuration) const {
  const auto *value = registration(record.typeIdentifier,
                                   record.contractVersion);
  if (value == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported observing-system identity or contract version."};
  if (!value->configurationResolver) {
    WVPortableTypedRecord candidate = record.configuration;
    if (candidate.schemaIdentifier.empty()) {
      candidate.schemaIdentifier =
          "source-linked-" + record.typeIdentifier + "-configuration-v1";
      candidate.schemaVersion = 1;
      const auto addText = [&](std::string name,
                               const std::vector<std::string> &values) {
        if (!values.empty())
          candidate.values.push_back(
              {std::move(name), {values.size()}, values});
      };
      const auto addReal = [&](std::string name,
                               const std::vector<double> &values) {
        if (!values.empty())
          candidate.values.push_back(
              {std::move(name), {values.size()}, values});
      };
      addText("fieldNames", record.fieldNames);
      addReal("x", record.x);
      addReal("y", record.y);
      addReal("z", record.z);
      candidate.values.push_back(
          {"isXYOnly", {},
           std::vector<std::uint8_t>{
               static_cast<std::uint8_t>(record.isXYOnly ? 1 : 0)}});
      candidate.values.push_back(
          {"shouldAntialias", {},
           std::vector<std::uint8_t>{
               static_cast<std::uint8_t>(record.shouldAntialias ? 1 : 0)}});
      candidate.values.push_back(
          {"outputScale", {}, std::vector<double>{record.outputScale}});
      candidate.values.push_back(
          {"outputOffset", {}, std::vector<double>{record.outputOffset}});
    }
    const auto status = validatePortableTypedRecord(
        candidate, {1024 * 1024, true, true});
    if (!status)
      return status;
    configuration = std::move(candidate);
    return WVKernelStatus::ok();
  }
  try {
    return value->configurationResolver(record, configuration);
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to resolve observer configuration."};
  } catch (const std::exception &error) {
    return invalid("Observer configuration resolver failed: " +
                   std::string(error.what()));
  } catch (...) {
    return invalid(
        "Observer configuration resolver failed with an unknown exception.");
  }
}

WVKernelStatus WVObserverCatalog::create(
    const WVObserverRecord &record,
    const WVPortableTypedRecord &configuration,
    std::shared_ptr<const WVObservingSystem> &result) const {
  result.reset();
  const auto *value = registration(record.typeIdentifier,
                                   record.contractVersion);
  if (value == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported observing-system identity or contract version."};
  std::shared_ptr<const WVObservingSystem> candidate;
  WVKernelStatus status;
  try {
    status = value->factory(record, configuration, candidate);
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to construct an observer implementation."};
  } catch (const std::exception &error) {
    return invalid("Observer factory failed: " + std::string(error.what()));
  } catch (...) {
    return invalid("Observer factory failed with an unknown exception.");
  }
  if (!status)
    return status;
  if (!candidate || candidate->typeIdentifier() != record.typeIdentifier ||
      candidate->contractVersion() != record.contractVersion)
    return invalid("An observer factory returned an incompatible implementation.");
  result = std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus WVObserverCatalog::resolveOutputPlan(
    const WVObserverRecord &record,
    const WVObserverOutputPlanningContext &context,
    WVObserverOutputPlan &plan) const {
  const auto *value = registration(record.typeIdentifier,
                                   record.contractVersion);
  if (value == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported observing-system identity or contract version."};
  if (!value->outputPlanResolver)
    return {WVKernelStatusCode::unsupportedOperation,
            "The observing-system registration does not declare a data-only "
            "output-plan resolver."};
  try {
    WVObserverOutputPlan candidate;
    auto status = value->outputPlanResolver(record, context, candidate);
    if (!status)
      return status;
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to resolve an observer output plan."};
  } catch (const std::exception &error) {
    return invalid("Observer output-plan resolver failed: " +
                   std::string(error.what()));
  } catch (...) {
    return invalid(
        "Observer output-plan resolver failed with an unknown exception.");
  }
}

std::size_t WVObserverCatalog::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      registrations_.capacity() *
                          sizeof(WVObserverFactoryRegistration);
  for (const auto &value : registrations_) {
    bytes += value.typeIdentifier.capacity() +
             detail::WVObserverFactoryRegistrationAccess::persistentBytes(
                 value);
  }
  return bytes;
}

const WVOutputScheduleFactoryRegistration *
WVOutputScheduleCatalog::registration(
    const std::string &typeIdentifier,
    std::uint32_t contractVersion) const noexcept {
  const auto found = std::find_if(
      registrations_.begin(), registrations_.end(), [&](const auto &value) {
        return value.typeIdentifier == typeIdentifier &&
               value.contractVersion == contractVersion;
      });
  return found == registrations_.end() ? nullptr : &*found;
}

const WVOutputScheduleFactoryRegistration *
WVOutputScheduleCatalog::registration(
    const WVOutputScheduleRecord &record) const noexcept {
  return record.typeIdentifier.empty()
             ? registration(WVEvenlySpacedOutputScheduleType, 1)
             : registration(record.typeIdentifier, record.contractVersion);
}

WVKernelStatus WVOutputScheduleCatalog::resolve(
    const WVOutputScheduleRecord &record,
    std::shared_ptr<const WVOutputSchedule> &result) const {
  result.reset();
  WVOutputScheduleRecord normalized = record;
  if (normalized.typeIdentifier.empty()) {
    normalized.typeIdentifier = WVEvenlySpacedOutputScheduleType;
    normalized.contractVersion = 1;
  }
  const auto *value =
      registration(normalized.typeIdentifier, normalized.contractVersion);
  if (value == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "No output-schedule implementation is available for " +
                normalized.typeIdentifier + " version " +
                std::to_string(normalized.contractVersion) + "."};
  WVKernelStatus status;
  std::shared_ptr<const WVOutputSchedule> candidate;
  try {
    candidate = value->factory(normalized, status);
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to construct an output schedule."};
  } catch (const std::exception &error) {
    return invalid("Output-schedule factory failed: " +
                   std::string(error.what()));
  } catch (...) {
    return invalid("Output-schedule factory failed with an unknown exception.");
  }
  if (!status)
    return status;
  if (!candidate || candidate->typeIdentifier() != normalized.typeIdentifier ||
      candidate->contractVersion() != normalized.contractVersion)
    return invalid("An output-schedule factory returned an incompatible implementation.");
  result = std::move(candidate);
  return WVKernelStatus::ok();
}

std::size_t WVOutputScheduleCatalog::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      registrations_.capacity() *
                          sizeof(WVOutputScheduleFactoryRegistration);
  for (const auto &value : registrations_)
    bytes += value.typeIdentifier.capacity();
  return bytes;
}

const WVForcingCatalog::Registration *WVForcingCatalog::registration(
    const std::string &typeIdentifier,
    std::uint32_t contractVersion) const noexcept {
  const auto found = std::find_if(
      registrations_.begin(), registrations_.end(), [&](const auto &value) {
        return value.matlabClassName == typeIdentifier &&
               value.contractVersion == contractVersion;
      });
  return found == registrations_.end() ? nullptr : &*found;
}

WVPortableCapability WVForcingCatalog::capability(
    std::string typeIdentifier, std::uint32_t contractVersion) const {
  auto found = std::find_if(
      registrations_.begin(), registrations_.end(), [&](const auto &value) {
        return value.matlabClassName == typeIdentifier &&
               value.contractVersion == contractVersion;
      });
  if (found == registrations_.end())
    found = std::find_if(
        registrations_.begin(), registrations_.end(), [&](const auto &value) {
          return value.matlabClassName == typeIdentifier;
        });
  const auto *value = found == registrations_.end() ? nullptr : &*found;
  std::optional<WVPortableImplementationIdentity> available;
  if (value != nullptr && value->isSupported)
    available = WVPortableImplementationIdentity{value->matlabClassName,
                                                  value->contractVersion};
  auto result = evaluatePortableCapability(
      {std::move(typeIdentifier), contractVersion}, std::move(available));
  if (value != nullptr && !value->isSupported &&
      result.status == WVPortableCapabilityStatus::unavailable)
    result.reason = value->unavailabilityReason;
  return result;
}

WVKernelStatus WVForcingCatalog::create(
    const WVFrozenForcingEntry &entry,
    const WVTransformConstantStratificationDescriptor &descriptor,
    bool hasAdaptiveDamping, std::unique_ptr<WVForcing> &forcing) const {
  forcing.reset();
  const auto *value = registration(entry.typeIdentifier,
                                   entry.contractVersion);
  if (value == nullptr || !value->isSupported || !value->factory)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported forcing identity."};
  try {
    std::unique_ptr<WVForcing> candidate;
    const auto status = value->factory(entry, descriptor, hasAdaptiveDamping,
                                       candidate);
    if (!status)
      return status;
    if (!candidate)
      return invalid("A forcing factory returned no implementation.");
    forcing = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to construct a forcing implementation."};
  } catch (const std::exception &error) {
    return invalid("Forcing factory failed: " + std::string(error.what()));
  } catch (...) {
    return invalid("Forcing factory failed with an unknown exception.");
  }
}

WVKernelStatus WVForcingCatalog::createBarotropicQG(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &descriptor,
    bool hasAdaptiveDamping,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) const {
  forcing.reset();
  const auto *value = registration(entry.typeIdentifier,
                                   entry.contractVersion);
  if (value == nullptr || !value->isSupported ||
      !value->barotropicQGFactory)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported Barotropic QG forcing identity."};
  try {
    std::unique_ptr<WVBarotropicQGForcing> candidate;
    const auto status = value->barotropicQGFactory(
        entry, descriptor, hasAdaptiveDamping, candidate);
    if (!status)
      return status;
    if (!candidate)
      return invalid("A Barotropic QG forcing factory returned no implementation.");
    if (candidate->typeIdentifier() != entry.typeIdentifier ||
        candidate->contractVersion() != entry.contractVersion)
      return invalid(
          "A Barotropic QG forcing factory returned an incompatible identity.");
    forcing = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to construct a Barotropic QG forcing implementation."};
  } catch (const std::exception &error) {
    return invalid("Barotropic QG forcing factory failed: " +
                   std::string(error.what()));
  } catch (...) {
    return invalid(
        "Barotropic QG forcing factory failed with an unknown exception.");
  }
}

WVKernelStatus WVForcingCatalog::validateConfiguration(
    const WVFrozenForcingEntry &entry) const {
  const auto *value = registration(entry.typeIdentifier,
                                   entry.contractVersion);
  if (value == nullptr || !value->isSupported ||
      value->contractVersion != entry.contractVersion)
    return invalid("Forcing configuration has no matching implementation.");
  const auto recordStatus = validatePortableTypedRecord(
      entry.configuration,
      {std::numeric_limits<std::size_t>::max(), false, true});
  if (!recordStatus)
    return recordStatus;
  if (entry.configuration.schemaIdentifier !=
          "wave-vortex-forcing-configuration-v1" ||
      entry.configuration.schemaVersion != 1)
    return invalid("Forcing configuration uses an unsupported schema.");
  std::set<std::string> allowed;
  for (const auto &field : value->persistence.fields) {
    allowed.insert(field.recordName);
    if (!field.imaginaryRecordName.empty())
      allowed.insert(field.imaginaryRecordName);
    const auto *stored = entry.configuration.value(field.recordName);
    if (stored == nullptr) {
      if (field.optional)
        continue;
      return invalid("Required forcing configuration value is missing.");
    }
    if (field.encoding == WVForcingPersistenceEncoding::realVariable) {
      const auto *reals = std::get_if<std::vector<double>>(&stored->storage);
      if (reals == nullptr)
        return invalid("Forcing real configuration has the wrong type.");
      for (const double scalar : *reals)
        if (std::isnan(scalar) ||
            (!field.allowInfinity && !std::isfinite(scalar)) ||
            (field.positive && !(scalar > 0.0)) ||
            (field.nonnegative && scalar < 0.0))
          return invalid("Forcing real configuration violates its bounds.");
    }
  }
  for (const auto &stored : entry.configuration.values)
    if (allowed.count(stored.name) == 0)
      return invalid("Forcing configuration contains an undeclared value.");
  return WVKernelStatus::ok();
}

std::size_t WVForcingCatalog::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      registrations_.capacity() * sizeof(Registration);
  for (const auto &value : registrations_) {
    bytes += value.matlabClassName.capacity() + value.defaultName.capacity() +
             value.unavailabilityReason.capacity() +
             value.forcingTypes.capacity() * sizeof(std::string) +
             value.persistence.fields.capacity() *
                 sizeof(WVForcingPersistenceField);
    for (const auto &type : value.forcingTypes)
      bytes += type.capacity();
    for (const auto &field : value.persistence.fields)
      bytes += field.recordName.capacity() +
               field.imaginaryRecordName.capacity() +
               field.netcdfName.capacity() +
               field.dimensionReference.capacity();
  }
  return bytes;
}

std::size_t WVExtensionCatalog::persistentBytes() const noexcept {
  return sizeof(*this) + observers_.persistentBytes() - sizeof(observers_) +
         outputSchedules_.persistentBytes() - sizeof(outputSchedules_) +
         forcings_.persistentBytes() - sizeof(forcings_);
}

WVKernelStatus WVExtensionCatalogBuilder::addObserverFactory(
    WVObserverFactoryRegistration registrationValue) {
  const auto reject = [this](std::string message) {
    valid_ = false;
    validationError_ = message;
    return invalid(std::move(message));
  };
  if (frozen_)
    return reject("A frozen extension-catalog builder cannot be mutated.");
  if (registrationValue.typeIdentifier.empty() ||
      registrationValue.contractVersion == 0 || !registrationValue.factory)
    return reject("An observer factory requires an identity, positive version, and function.");
  if (containsIdentity(observers_, [&](const auto &value) {
        return value.typeIdentifier == registrationValue.typeIdentifier &&
               value.contractVersion == registrationValue.contractVersion;
      }))
    return reject("An observer factory identity/version pair is duplicated.");
  try {
    observers_.push_back(std::move(registrationValue));
  } catch (const std::bad_alloc &) {
    valid_ = false;
    validationError_ = "Observer registration allocation failed.";
    return {WVKernelStatusCode::allocationFailure, validationError_};
  } catch (const std::exception &error) {
    return reject("Observer registration failed: " +
                  std::string(error.what()));
  } catch (...) {
    return reject("Observer registration failed with an unknown exception.");
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVExtensionCatalogBuilder::addOutputScheduleFactory(
    WVOutputScheduleFactoryRegistration registrationValue) {
  const auto reject = [this](std::string message) {
    valid_ = false;
    validationError_ = message;
    return invalid(std::move(message));
  };
  if (frozen_)
    return reject("A frozen extension-catalog builder cannot be mutated.");
  if (registrationValue.typeIdentifier.empty() ||
      registrationValue.contractVersion == 0 ||
      registrationValue.factory == nullptr)
    return reject("An output-schedule factory requires an identity, positive version, and function.");
  if (containsIdentity(outputSchedules_, [&](const auto &value) {
        return value.typeIdentifier == registrationValue.typeIdentifier &&
               value.contractVersion == registrationValue.contractVersion;
      }))
    return reject("An output-schedule factory identity/version pair is duplicated.");
  try {
    outputSchedules_.push_back(std::move(registrationValue));
  } catch (const std::bad_alloc &) {
    valid_ = false;
    validationError_ = "Output-schedule registration allocation failed.";
    return {WVKernelStatusCode::allocationFailure, validationError_};
  } catch (const std::exception &error) {
    return reject("Output-schedule registration failed: " +
                  std::string(error.what()));
  } catch (...) {
    return reject(
        "Output-schedule registration failed with an unknown exception.");
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVExtensionCatalogBuilder::addForcingFactory(
    WVForcingCatalog::Registration registrationValue) {
  const auto reject = [this](std::string message) {
    valid_ = false;
    validationError_ = message;
    return invalid(std::move(message));
  };
  if (frozen_)
    return reject("A frozen extension-catalog builder cannot be mutated.");
  if (registrationValue.matlabClassName.empty() ||
      registrationValue.contractVersion == 0)
    return reject("A forcing registration requires an identity and positive version.");
  if ((registrationValue.isSupported && !registrationValue.factory &&
       !registrationValue.barotropicQGFactory) ||
      (registrationValue.isSupported &&
       !registrationValue.unavailabilityReason.empty()) ||
      (!registrationValue.isSupported &&
       (registrationValue.factory ||
        registrationValue.barotropicQGFactory)) ||
      (!registrationValue.isSupported &&
       registrationValue.unavailabilityReason.empty()))
    return reject("A forcing registration is incomplete or conflicting.");
  if (containsIdentity(forcings_, [&](const auto &value) {
        return value.matlabClassName == registrationValue.matlabClassName &&
               value.contractVersion == registrationValue.contractVersion;
      }))
    return reject("A forcing registration identity is duplicated.");
  try {
    forcings_.push_back(std::move(registrationValue));
  } catch (const std::bad_alloc &) {
    valid_ = false;
    validationError_ = "Forcing registration allocation failed.";
    return {WVKernelStatusCode::allocationFailure, validationError_};
  } catch (const std::exception &error) {
    return reject("Forcing registration failed: " +
                  std::string(error.what()));
  } catch (...) {
    return reject("Forcing registration failed with an unknown exception.");
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVExtensionCatalogBuilder::freeze(
    std::shared_ptr<const WVExtensionCatalog> &catalog) {
  catalog.reset();
  if (frozen_)
    return invalid("An extension-catalog builder can be frozen only once.");
  if (!valid_)
    return invalid("The extension-catalog builder is invalid: " +
                   validationError_);
  try {
    auto candidate = std::make_shared<WVExtensionCatalog>();
    candidate->observers_.registrations_ = observers_;
    candidate->outputSchedules_.registrations_ = outputSchedules_;
    candidate->forcings_.registrations_ = forcings_;
    catalog = std::move(candidate);
    frozen_ = true;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to freeze the extension catalog."};
  } catch (const std::exception &error) {
    return invalid("Extension-catalog freeze failed: " +
                   std::string(error.what()));
  } catch (...) {
    return invalid("Extension-catalog freeze failed with an unknown exception.");
  }
}

WVKernelStatus addBuiltInExtensions(WVExtensionCatalogBuilder &builder) {
  auto status = detail::addBuiltInObserverFactories(builder);
  if (!status)
    return status;
  status = builder.addOutputScheduleFactory(
      {WVEvenlySpacedOutputScheduleType, 1,
       &makeEvenlySpacedOutputSchedule});
  if (!status)
    return status;
  for (auto registration : builtInForcingFactories()) {
    status = builder.addForcingFactory(std::move(registration));
    if (!status)
      return status;
  }
  return WVKernelStatus::ok();
}

WVKernelStatus makeBuiltInExtensionCatalog(
    std::shared_ptr<const WVExtensionCatalog> &catalog) {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (!status)
    return status;
  return builder.freeze(catalog);
}

} // namespace wavevortex::runtime
