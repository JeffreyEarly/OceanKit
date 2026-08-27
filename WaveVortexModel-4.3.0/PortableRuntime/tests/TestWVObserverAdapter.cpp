#include "WVObserverAdapter.hpp"
#include "WVTestExtensionCatalog.hpp"

#include <cstdlib>
#include <iostream>
#include <set>
#include <string>

using namespace wavevortex::runtime;
using wavevortex::WVKernelStatus;
using wavevortex::WVKernelStatusCode;

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

class WVTestPortablePointDiagnostic final : public WVObservingSystem {
public:
  explicit WVTestPortablePointDiagnostic(WVPortableTypedRecord configuration)
      : configuration_(std::move(configuration)) {}
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestPortablePointDiagnostic";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 1; }
  WVKernelStatus executionPlan(const WVObserverRecord &record,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.fieldListAttribute = "fieldNames";
    plan.persistedName = record.name;
    const auto *field = record.configuration.value("field");
    plan.outputFields =
        field == nullptr ? record.fieldNames
                         : std::get<std::vector<std::string>>(field->storage);
    return WVKernelStatus::ok();
  }
  WVKernelStatus validate(
      const WVObserverRecord &record,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    if (!record.stateBlockIdentifiers.empty() || record.fieldNames.size() != 1 ||
        record.x.empty() || record.x.size() != record.y.size() ||
        record.x.size() != record.z.size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Point diagnostic requires one field and equal fixed x/y/z points."};
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + configuration_.persistentBytes() -
           sizeof(configuration_);
  }

private:
  WVPortableTypedRecord configuration_;
};

class WVTestPortablePointDiagnosticV2 final : public WVObservingSystem {
public:
  const std::string &typeIdentifier() const noexcept override {
    static const std::string value = "WVTestPortablePointDiagnostic";
    return value;
  }
  std::uint32_t contractVersion() const noexcept override { return 2; }
  WVKernelStatus executionPlan(const WVObserverRecord &record,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    plan.fieldListAttribute = "fieldNames";
    plan.persistedName = record.name;
    plan.outputFields = record.fieldNames;
    return WVKernelStatus::ok();
  }
  WVKernelStatus validate(
      const WVObserverRecord &,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
};

} // namespace

int main() {
  WVExtensionCatalogBuilder builder;
  require(static_cast<bool>(addBuiltInExtensions(builder)),
          "built-in extension registration failed");
  const auto pointFactory =
      [](const WVObserverRecord &, const WVPortableTypedRecord &configuration,
         std::shared_ptr<const WVObservingSystem> &result) {
        result =
            std::make_shared<WVTestPortablePointDiagnostic>(configuration);
        return WVKernelStatus::ok();
      };
  require(static_cast<bool>(builder.addObserverFactory(
              {"WVTestPortablePointDiagnostic", 1, pointFactory})),
          "test observer factory registration failed");
  require(!builder.addObserverFactory(
              {"WVTestPortablePointDiagnostic", 1, pointFactory}),
          "duplicate observer factory registration succeeded");

  WVExtensionCatalogBuilder catalogBuilder;
  require(static_cast<bool>(addBuiltInExtensions(catalogBuilder)),
          "built-in extension registration failed");
  require(static_cast<bool>(catalogBuilder.addObserverFactory(
              {"WVTestPortablePointDiagnostic", 1, pointFactory})),
          "test observer factory registration failed");
  require(static_cast<bool>(catalogBuilder.addObserverFactory(
              {"WVTestPortablePointDiagnostic", 2,
               [](const WVObserverRecord &, const WVPortableTypedRecord &,
                  std::shared_ptr<const WVObservingSystem> &result) {
                 result =
                     std::make_shared<WVTestPortablePointDiagnosticV2>();
                 return WVKernelStatus::ok();
               }})),
          "a second observer contract version could not be registered");
  std::shared_ptr<const WVExtensionCatalog> catalog;
  require(static_cast<bool>(catalogBuilder.freeze(catalog)),
          "extension catalog freeze failed");
  require(catalog->observers().capability(
              "WVTestPortablePointDiagnostic", 1)
              .isSupported(),
          "registered observer pair is unavailable");
  require(catalog->observers().capability(
              "WVTestPortablePointDiagnostic", 2)
              .isSupported(),
          "the second observer contract version is unavailable");
  require(catalog->observers().capability(
              "WVTestPortablePointDiagnostic", 3)
              .status == WVPortableCapabilityStatus::versionMismatch,
          "an unavailable observer contract version was accepted");
  require(catalog->observers().capability("WVCustomObserver", 1).status ==
              WVPortableCapabilityStatus::unavailable,
          "missing observer pair did not report unavailability");

  WVPortableObserverRecord record;
  WVObserverRecord observer;
  observer.identifier = "test-point-diagnostic";
  observer.name = "test point diagnostic";
  observer.typeIdentifier = "WVTestPortablePointDiagnostic";
  observer.fieldNames = {"u"};
  observer.x = {0.1, 0.2};
  observer.y = {0.3, 0.4};
  observer.z = {-0.5, -0.6};
  observer.configuration.schemaIdentifier = "test-point-configuration-v1";
  observer.configuration.schemaVersion = 1;
  observer.configuration.values = {
      {"field", {}, std::vector<std::string>{"u"}},
      {"x", {2}, std::vector<double>{0.1, 0.2}},
      {"y", {2}, std::vector<double>{0.3, 0.4}},
      {"z", {2}, std::vector<double>{-0.5, -0.6}}};
  record.observers.push_back(observer);
  observer.identifier = "test-point-diagnostic-2";
  observer.name = "second test point diagnostic";
  record.observers.push_back(observer);
  WVPortableObserverDescriptor descriptor;
  const auto descriptorStatus =
      WVPortableObserverDescriptor::create(record, catalog, descriptor);
  require(static_cast<bool>(descriptorStatus),
          descriptorStatus.message.c_str());
  require(descriptor.implementation(descriptor.observers().front()) != nullptr,
          "descriptor did not retain a resolved implementation");
  const auto *firstResolved =
      descriptor.resolvedObserver(descriptor.observers()[0]);
  const auto *secondResolved =
      descriptor.resolvedObserver(descriptor.observers()[1]);
  require(firstResolved != nullptr && secondResolved != nullptr &&
              firstResolved != secondResolved &&
              &firstResolved->implementation() !=
                  &secondResolved->implementation() &&
              firstResolved->configuration().schemaIdentifier ==
                  "test-point-configuration-v1" &&
              firstResolved->configuration().value("x") != nullptr &&
              firstResolved->executionPlan().outputFields ==
                  std::vector<std::string>{"u"},
          "descriptor did not create immutable per-record resolved observers");
  require(!catalogBuilder.addObserverFactory(
              {"WVLateObserver", 1, pointFactory}),
          "late observer factory registration succeeded");

  WVObserverRecord particles;
  particles.name = "floats";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.isXYOnly = false;
  auto channels = detail::particlePositionChannels(particles.isXYOnly);
  require(channels ==
              std::vector<detail::WVMovingFieldChannel>{
                  detail::WVMovingFieldChannel::x,
                  detail::WVMovingFieldChannel::y,
                  detail::WVMovingFieldChannel::z},
          "three-dimensional particle channels are not x/y/z");
  particles.isXYOnly = true;
  require(detail::particlePositionChannels(particles.isXYOnly).size() == 2,
          "two-dimensional particles must expose x/y channels only");

  WVObserverRecord tracer;
  tracer.name = "dye";
  tracer.typeIdentifier = "WVTracer";
  channels = {detail::WVMovingFieldChannel::tracerValue};
  require(channels == std::vector<detail::WVMovingFieldChannel>{
                          detail::WVMovingFieldChannel::tracerValue},
          "tracer must expose one named value channel");
  require(detail::movingFieldVariableName(tracer, channels.front()) == "dye",
          "tracer variable name changed");
  return 0;
}
