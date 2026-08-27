#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WVTestLinearCoefficientForcing.hpp"
#include "WVTestQuadraticSchedule.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

class TestObserver final : public WVObservingSystem {
public:
  TestObserver(std::string identity, std::uint32_t version,
               WVPortableTypedRecord configuration)
      : identity_(std::move(identity)), version_(version),
        configuration_(std::move(configuration)) {}
  const std::string &typeIdentifier() const noexcept override {
    return identity_;
  }
  std::uint32_t contractVersion() const noexcept override { return version_; }
  WVKernelStatus executionPlan(const WVObserverRecord &,
                               WVObserverExecutionPlan &plan) const override {
    plan = {};
    return WVKernelStatus::ok();
  }
  WVKernelStatus validate(
      const WVObserverRecord &,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + identity_.capacity() +
           configuration_.persistentBytes() - sizeof(configuration_);
  }

private:
  std::string identity_;
  std::uint32_t version_;
  WVPortableTypedRecord configuration_;
};

WVPortableTypedRecord configuration(double value) {
  WVPortableTypedRecord result;
  result.schemaIdentifier = "test-observer-configuration-v1";
  result.schemaVersion = 1;
  result.values.push_back({"value", {}, std::vector<double>{value}});
  return result;
}

WVObserverFactory observerFactory(std::string identity, std::uint32_t version,
                                  int *constructionCount = nullptr) {
  return [identity = std::move(identity), version, constructionCount](
             const WVObserverRecord &, const WVPortableTypedRecord &record,
             std::shared_ptr<const WVObservingSystem> &result) {
    if (constructionCount != nullptr)
      ++*constructionCount;
    result = std::make_shared<TestObserver>(identity, version, record);
    return WVKernelStatus::ok();
  };
}

std::shared_ptr<const WVOutputSchedule>
throwBadAllocationSchedule(const WVOutputScheduleRecord &, WVKernelStatus &) {
  throw std::bad_alloc();
}

std::shared_ptr<const WVOutputSchedule>
throwRuntimeSchedule(const WVOutputScheduleRecord &, WVKernelStatus &) {
  throw std::runtime_error("schedule construction failure");
}

std::shared_ptr<const WVOutputSchedule>
failSchedule(const WVOutputScheduleRecord &, WVKernelStatus &status) {
  status = {WVKernelStatusCode::invalidConfiguration,
            "schedule status failure"};
  return {};
}

WVTransformConstantStratificationConfiguration modelConfiguration() {
  WVTransformConstantStratificationConfiguration configuration;
  configuration.Nx = 4;
  configuration.Ny = 4;
  configuration.Nz = 5;
  configuration.Nj = 3;
  configuration.Lx = 1000.0;
  configuration.Ly = 1000.0;
  configuration.Lz = 100.0;
  configuration.N0 = 5e-3;
  configuration.rho0 = 1025.0;
  configuration.g = 9.81;
  configuration.planetaryRadius = 6.371e6;
  configuration.rotationRate = 7.292115e-5;
  configuration.latitude = 30.0;
  return configuration;
}

WVTransformConstantStratificationDescriptor transformDescriptor() {
  const auto modelConfig = modelConfiguration();
  WVTransformConstantStratificationDescriptor result;
  const auto status =
      WVTransformConstantStratificationDescriptor::create(modelConfig, result);
  if (!status)
    throw std::runtime_error(status.message);
  return result;
}

void addCanonicalBlocks(
    WVPortableObserverRecord &record,
    const WVTransformConstantStratificationConfiguration &configuration) {
  WVTransformConstantStratificationDescriptor transform;
  const auto status =
      WVTransformConstantStratificationDescriptor::create(configuration,
                                                          transform);
  if (!status)
    throw std::runtime_error(status.message);
  for (const char *identifier : {"Ap", "Am", "A0"}) {
    WVStateBlockRecord block;
    block.identifier = identifier;
    block.scalarType = WVStateScalarType::complex64;
    block.dimensions = {configuration.Nj, transform.Nkl()};
    block.toleranceKind = WVToleranceKind::coefficientEnergyScaled;
    block.absoluteTolerance = 1e-6;
    record.stateBlocks.push_back(std::move(block));
  }
}

WVForcingFactoryRegistration throwingForcing(std::string identity,
                                             bool badAllocation) {
  auto registration =
      wavevortex::runtime::test::linearCoefficientRegistration();
  registration.matlabClassName = std::move(identity);
  registration.factory =
      [badAllocation](const WVFrozenForcingEntry &,
                      const WVTransformConstantStratificationDescriptor &, bool,
                      std::unique_ptr<WVForcing> &) -> WVKernelStatus {
    if (badAllocation)
      throw std::bad_alloc();
    throw std::runtime_error("forcing construction failure");
  };
  return registration;
}

} // namespace

int main() {
  std::shared_ptr<const WVExtensionCatalog> builtIns;
  require(static_cast<bool>(makeBuiltInExtensionCatalog(builtIns)),
          "built-in catalog construction failed");

  // A failed add poisons the builder and a failed freeze clears any prior
  // usable output catalog rather than publishing a valid subset.
  WVExtensionCatalogBuilder duplicateBuilder;
  require(static_cast<bool>(duplicateBuilder.addObserverFactory(
              {"TestObserver", 1, observerFactory("TestObserver", 1)})),
          "initial observer registration failed");
  require(!duplicateBuilder.addObserverFactory(
              {"TestObserver", 1, observerFactory("TestObserver", 1)}),
          "duplicate observer registration succeeded");
  auto failedCatalog = builtIns;
  require(!duplicateBuilder.freeze(failedCatalog) && !failedCatalog,
          "failed freeze retained a usable catalog");

  std::shared_ptr<const WVExtensionCatalog> firstCatalog;
  {
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(addBuiltInExtensions(builder)) &&
                static_cast<bool>(builder.addObserverFactory(
                {"TestObserver", 1, observerFactory("TestObserver", 1)})),
            "first catalog observer registration failed");
    require(static_cast<bool>(builder.freeze(firstCatalog)),
            "first catalog freeze failed");
    require(!builder.addObserverFactory(
                {"LateObserver", 1, observerFactory("LateObserver", 1)}),
            "late registration succeeded");
  }
  require(firstCatalog->observers().registration("TestObserver", 1) != nullptr,
          "builder destruction invalidated the frozen catalog");

  WVExtensionCatalogBuilder secondBuilder;
  require(static_cast<bool>(addBuiltInExtensions(secondBuilder)) &&
              static_cast<bool>(secondBuilder.addObserverFactory(
              {"OtherObserver", 1, observerFactory("OtherObserver", 1)})),
          "second catalog observer registration failed");
  std::shared_ptr<const WVExtensionCatalog> secondCatalog;
  require(static_cast<bool>(secondBuilder.freeze(secondCatalog)),
          "second catalog freeze failed");
  require(firstCatalog->observers().registration("OtherObserver", 1) == nullptr &&
              secondCatalog->observers().registration("TestObserver", 1) ==
                  nullptr,
          "independent catalogs leaked capabilities");

  WVPortableObserverRecord firstRuntimeRecord;
  WVObserverRecord firstRuntimeObserver;
  firstRuntimeObserver.identifier = "first-runtime";
  firstRuntimeObserver.name = "first-runtime";
  firstRuntimeObserver.typeIdentifier = "TestObserver";
  firstRuntimeObserver.configuration = configuration(10.0);
  firstRuntimeRecord.observers.push_back(std::move(firstRuntimeObserver));
  const auto modelConfig = modelConfiguration();
  addCanonicalBlocks(firstRuntimeRecord, modelConfig);
  WVPortableObserverDescriptor firstRuntime;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(
              firstRuntimeRecord, firstCatalog, firstRuntime)),
          "first independent resolved runtime failed");
  WVPortableObserverRecord secondRuntimeRecord;
  WVObserverRecord secondRuntimeObserver;
  secondRuntimeObserver.identifier = "second-runtime";
  secondRuntimeObserver.name = "second-runtime";
  secondRuntimeObserver.typeIdentifier = "OtherObserver";
  secondRuntimeObserver.configuration = configuration(20.0);
  secondRuntimeRecord.observers.push_back(std::move(secondRuntimeObserver));
  addCanonicalBlocks(secondRuntimeRecord, modelConfig);
  WVPortableObserverDescriptor secondRuntime;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(
              secondRuntimeRecord, secondCatalog, secondRuntime)),
          "second independent resolved runtime failed");
  WVModelIntegratorConfiguration integratorConfiguration;
  WVModel firstModel;
  WVModel secondModel;
  const auto schedule = defaultNonlinearAdvectionSchedule();
  std::unique_ptr<WVConstantStratificationIntegrationSystem> mismatchedSystem;
  auto modelStatus = WVConstantStratificationIntegrationSystem::create(
      modelConfig, schedule, firstRuntime, secondCatalog,
      std::make_unique<WVReferenceFFTEngine>(), mismatchedSystem);
  require(modelStatus.code == WVKernelStatusCode::invalidConfiguration &&
              !mismatchedSystem,
          "integration system accepted a descriptor from another catalog");
  WVModel mismatchedModel;
  modelStatus = WVModel::create(
      secondCatalog, modelConfig, schedule, firstRuntime,
      std::make_unique<WVReferenceFFTEngine>(), integratorConfiguration,
      mismatchedModel);
  require(modelStatus.code == WVKernelStatusCode::invalidConfiguration,
          "WVModel accepted a descriptor from another catalog");
  modelStatus = WVModel::create(
      firstCatalog, modelConfig, schedule, firstRuntime,
      std::make_unique<WVReferenceFFTEngine>(), integratorConfiguration,
      firstModel);
  if (!modelStatus)
    throw std::runtime_error("first independent WVModel: " +
                             modelStatus.message);
  modelStatus = WVModel::create(
      secondCatalog, modelConfig, schedule, secondRuntime,
      std::make_unique<WVReferenceFFTEngine>(), integratorConfiguration,
      secondModel);
  if (!modelStatus)
    throw std::runtime_error("second independent WVModel: " +
                             modelStatus.message);
  WVModel movedFirstModel = std::move(firstModel);
  firstCatalog.reset();
  secondCatalog.reset();
  const auto secondMetrics = secondModel.metrics();
  const auto firstMetrics = movedFirstModel.metrics();
  require(firstRuntime.catalog() != nullptr &&
              secondRuntime.catalog() != nullptr &&
              firstMetrics.catalogPersistentBytes > 0 &&
              secondMetrics.catalogPersistentBytes > 0 &&
              movedFirstModel.forcingScheduleIdentifier() ==
                  secondModel.forcingScheduleIdentifier(),
          "interleaved runtime/catalog destruction was not independent");

  // Malformed registrations and duplicate identities poison their builder and
  // cannot be frozen into a usable partial catalog.
  const auto requireRejectedBuilder = [](WVExtensionCatalogBuilder &builder) {
    std::shared_ptr<const WVExtensionCatalog> output;
    require(!builder.freeze(output) && !output,
            "invalid builder froze a partial catalog");
  };
  WVExtensionCatalogBuilder zeroVersionBuilder;
  require(!zeroVersionBuilder.addObserverFactory(
              {"ZeroVersion", 0, observerFactory("ZeroVersion", 1)}),
          "zero-version observer registration succeeded");
  requireRejectedBuilder(zeroVersionBuilder);
  WVExtensionCatalogBuilder malformedObserverBuilder;
  require(!malformedObserverBuilder.addObserverFactory(
              {"", 1, observerFactory("Malformed", 1)}),
          "identity-free observer registration succeeded");
  requireRejectedBuilder(malformedObserverBuilder);
  WVExtensionCatalogBuilder duplicateScheduleBuilder;
  require(static_cast<bool>(duplicateScheduleBuilder.addOutputScheduleFactory(
              {wavevortex::runtime::test::quadraticScheduleType, 1,
               &wavevortex::runtime::test::makeQuadraticSchedule})) &&
              !duplicateScheduleBuilder.addOutputScheduleFactory(
                  {wavevortex::runtime::test::quadraticScheduleType, 1,
                   &wavevortex::runtime::test::makeQuadraticSchedule}),
          "duplicate schedule registration was not rejected");
  requireRejectedBuilder(duplicateScheduleBuilder);
  WVExtensionCatalogBuilder duplicateForcingBuilder;
  require(static_cast<bool>(duplicateForcingBuilder.addForcingFactory(
              wavevortex::runtime::test::linearCoefficientRegistration())) &&
              !duplicateForcingBuilder.addForcingFactory(
                  wavevortex::runtime::test::linearCoefficientRegistration()),
          "duplicate forcing registration was not rejected");
  requireRejectedBuilder(duplicateForcingBuilder);
  auto incompleteForcing =
      wavevortex::runtime::test::linearCoefficientRegistration();
  incompleteForcing.factory = {};
  WVExtensionCatalogBuilder incompleteForcingBuilder;
  require(!incompleteForcingBuilder.addForcingFactory(
              std::move(incompleteForcing)),
          "incomplete supported forcing registration succeeded");
  requireRejectedBuilder(incompleteForcingBuilder);
  auto unavailableWithFactory =
      wavevortex::runtime::test::linearCoefficientRegistration();
  unavailableWithFactory.isSupported = false;
  unavailableWithFactory.unavailabilityReason = "intentionally unavailable";
  WVExtensionCatalogBuilder unavailableWithFactoryBuilder;
  require(!unavailableWithFactoryBuilder.addForcingFactory(
              std::move(unavailableWithFactory)),
          "unavailable forcing registration retained a factory");
  requireRejectedBuilder(unavailableWithFactoryBuilder);
  auto supportedWithReason =
      wavevortex::runtime::test::linearCoefficientRegistration();
  supportedWithReason.unavailabilityReason = "conflicts with supported state";
  WVExtensionCatalogBuilder supportedWithReasonBuilder;
  require(!supportedWithReasonBuilder.addForcingFactory(
              std::move(supportedWithReason)),
          "supported forcing registration retained an unavailability reason");
  requireRejectedBuilder(supportedWithReasonBuilder);

  // Exact version pairs coexist and resolve independently.
  WVExtensionCatalogBuilder versionBuilder;
  require(static_cast<bool>(versionBuilder.addObserverFactory(
              {"VersionedObserver", 1,
               observerFactory("VersionedObserver", 1)})) &&
              static_cast<bool>(versionBuilder.addObserverFactory(
                  {"VersionedObserver", 2,
                   observerFactory("VersionedObserver", 2)})),
          "versioned observer registration failed");
  std::shared_ptr<const WVExtensionCatalog> versionCatalog;
  require(static_cast<bool>(versionBuilder.freeze(versionCatalog)) &&
              versionCatalog->observers()
                      .capability("VersionedObserver", 2)
                      .status == WVPortableCapabilityStatus::supported &&
              versionCatalog->observers()
                      .capability("VersionedObserver", 3)
                      .status == WVPortableCapabilityStatus::versionMismatch &&
              versionCatalog->observers()
                      .capability("AbsentObserver", 1)
                      .status == WVPortableCapabilityStatus::unavailable,
          "observer exact-version capability semantics failed");

  WVExtensionCatalogBuilder forcingVersionBuilder;
  auto forcingV1 =
      wavevortex::runtime::test::linearCoefficientRegistration();
  auto forcingV2 = forcingV1;
  forcingV2.contractVersion = 2;
  require(static_cast<bool>(forcingVersionBuilder.addForcingFactory(
              std::move(forcingV1))) &&
              static_cast<bool>(forcingVersionBuilder.addForcingFactory(
                  std::move(forcingV2))),
          "forcing exact-version registration failed");
  auto unavailableForcing =
      wavevortex::runtime::test::linearCoefficientRegistration();
  unavailableForcing.matlabClassName = "UnavailableForcing";
  unavailableForcing.factory = {};
  unavailableForcing.isSupported = false;
  unavailableForcing.unavailabilityReason = "not linked in this catalog";
  require(static_cast<bool>(forcingVersionBuilder.addForcingFactory(
              std::move(unavailableForcing))),
          "unavailable forcing registration failed");
  std::shared_ptr<const WVExtensionCatalog> forcingVersionCatalog;
  require(static_cast<bool>(forcingVersionBuilder.freeze(
              forcingVersionCatalog)) &&
              forcingVersionCatalog->forcings()
                      .capability(
                          wavevortex::runtime::test::
                              LinearCoefficientForcingIdentifier,
                          2)
                      .status == WVPortableCapabilityStatus::supported &&
              forcingVersionCatalog->forcings()
                      .capability(
                          wavevortex::runtime::test::
                              LinearCoefficientForcingIdentifier,
                          3)
                      .status == WVPortableCapabilityStatus::versionMismatch &&
              forcingVersionCatalog->forcings()
                      .capability("UnavailableForcing", 1)
                      .status == WVPortableCapabilityStatus::unavailable,
          "forcing version/unavailable capability semantics failed");

  // Descriptor resolution creates exactly one distinct immutable
  // implementation for each record and retains the catalog lifetime.
  int constructions = 0;
  WVExtensionCatalogBuilder recordBuilder;
  require(static_cast<bool>(recordBuilder.addObserverFactory(
              {"PerRecordObserver", 1,
               observerFactory("PerRecordObserver", 1, &constructions)})),
          "per-record observer registration failed");
  std::shared_ptr<const WVExtensionCatalog> recordCatalog;
  require(static_cast<bool>(recordBuilder.freeze(recordCatalog)),
          "per-record catalog freeze failed");
  WVPortableObserverRecord record;
  WVObserverRecord firstRecord;
  firstRecord.identifier = "first";
  firstRecord.name = "first";
  firstRecord.typeIdentifier = "PerRecordObserver";
  firstRecord.contractVersion = 1;
  firstRecord.configuration = configuration(1.0);
  WVObserverRecord secondRecord = firstRecord;
  secondRecord.identifier = "second";
  secondRecord.name = "second";
  secondRecord.configuration = configuration(2.0);
  record.observers = {std::move(firstRecord), std::move(secondRecord)};
  WVPortableObserverDescriptor descriptor;
  require(static_cast<bool>(WVPortableObserverDescriptor::create(
              record, recordCatalog, descriptor)),
          "per-record descriptor construction failed");
  const auto *first = descriptor.resolvedObserver(descriptor.observers()[0]);
  const auto *second = descriptor.resolvedObserver(descriptor.observers()[1]);
  require(constructions == 2 && first != nullptr && second != nullptr &&
              first->implementationHandle() != second->implementationHandle(),
          "observer records shared one implementation");
  record = {};
  recordCatalog.reset();
  const auto *firstValue =
      descriptor.resolvedObserver(descriptor.observers()[0])
          ->configuration()
          .value("value");
  const auto *secondValue =
      descriptor.resolvedObserver(descriptor.observers()[1])
          ->configuration()
          .value("value");
  require(descriptor.catalog() != nullptr && firstValue != nullptr &&
              secondValue != nullptr &&
              std::get<std::vector<double>>(firstValue->storage).front() == 1.0 &&
              std::get<std::vector<double>>(secondValue->storage).front() == 2.0,
          "descriptor did not retain immutable catalog ownership");

  // Factory failures are status-normalized and output handles are
  // transactional in all three typed subcatalogs.
  for (const bool badAllocation : {true, false}) {
    WVExtensionCatalogBuilder builder;
    WVObserverFactory factory =
        [badAllocation](const WVObserverRecord &,
                        const WVPortableTypedRecord &,
                        std::shared_ptr<const WVObservingSystem> &result)
        -> WVKernelStatus {
      result = std::make_shared<TestObserver>("ThrowingObserver", 1,
                                              configuration(0.0));
      if (badAllocation)
        throw std::bad_alloc();
      throw std::runtime_error("observer construction failure");
    };
    require(static_cast<bool>(builder.addObserverFactory(
                {"ThrowingObserver", 1, std::move(factory)})),
            "throwing observer registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "throwing observer catalog freeze failed");
    WVObserverRecord observer;
    observer.typeIdentifier = "ThrowingObserver";
    observer.contractVersion = 1;
    std::shared_ptr<const WVObservingSystem> output =
        std::make_shared<TestObserver>("prior", 1, configuration(0.0));
    const auto status = catalog->observers().create(
        observer, configuration(0.0), output);
    require(!status && !output &&
                status.code ==
                    (badAllocation ? WVKernelStatusCode::allocationFailure
                                   : WVKernelStatusCode::invalidConfiguration),
            "observer factory failure was not transactional/normalized");
  }

  {
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(builder.addObserverFactory(
                {"StatusObserver", 1,
                 [](const WVObserverRecord &, const WVPortableTypedRecord &,
                    std::shared_ptr<const WVObservingSystem> &result) {
                   result = std::make_shared<TestObserver>(
                       "StatusObserver", 1, configuration(0.0));
                   return WVKernelStatus{
                       WVKernelStatusCode::invalidConfiguration,
                       "observer status failure"};
                 }})),
            "status-failing observer registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "status-failing observer catalog freeze failed");
    WVObserverRecord observer;
    observer.typeIdentifier = "StatusObserver";
    std::shared_ptr<const WVObservingSystem> output =
        std::make_shared<TestObserver>("prior", 1, configuration(0.0));
    const auto status =
        catalog->observers().create(observer, configuration(0.0), output);
    require(!status && !output,
            "observer status failure published a partial implementation");
  }

  for (const bool badAllocation : {true, false}) {
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(builder.addOutputScheduleFactory(
                {badAllocation ? "BadAllocSchedule" : "RuntimeSchedule", 1,
                 badAllocation ? &throwBadAllocationSchedule
                               : &throwRuntimeSchedule})),
            "throwing schedule registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "throwing schedule catalog freeze failed");
    WVOutputScheduleRecord schedule;
    schedule.typeIdentifier =
        badAllocation ? "BadAllocSchedule" : "RuntimeSchedule";
    schedule.contractVersion = 1;
    std::shared_ptr<const WVOutputSchedule> output;
    WVKernelStatus initialStatus;
    output = makeEvenlySpacedOutputSchedule({1.0, 0.0, 2.0}, initialStatus);
    const auto status = catalog->outputSchedules().resolve(schedule, output);
    require(!status && !output &&
                status.code ==
                    (badAllocation ? WVKernelStatusCode::allocationFailure
                                   : WVKernelStatusCode::invalidConfiguration),
            "schedule factory failure was not transactional/normalized");
  }

  {
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(builder.addOutputScheduleFactory(
                {"StatusSchedule", 1, &failSchedule})),
            "status-failing schedule registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "status-failing schedule catalog freeze failed");
    WVOutputScheduleRecord schedule;
    schedule.typeIdentifier = "StatusSchedule";
    std::shared_ptr<const WVOutputSchedule> output;
    WVKernelStatus initialStatus;
    output = makeEvenlySpacedOutputSchedule({1.0, 0.0, 2.0}, initialStatus);
    const auto status = catalog->outputSchedules().resolve(schedule, output);
    require(!status && !output,
            "schedule status failure published a partial implementation");
  }

  const auto transform = transformDescriptor();
  for (const bool badAllocation : {true, false}) {
    const std::string identity =
        badAllocation ? "BadAllocForcing" : "RuntimeForcing";
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(
                builder.addForcingFactory(throwingForcing(identity,
                                                          badAllocation))),
            "throwing forcing registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "throwing forcing catalog freeze failed");
    WVFrozenForcingEntry entry;
    entry.typeIdentifier = identity;
    entry.contractVersion = WVPortablePairContractVersion;
    std::unique_ptr<WVForcing> output =
        std::make_unique<wavevortex::runtime::test::
                             WVTestPortableLinearCoefficientForcing>(entry, 1.0);
    const auto status =
        catalog->forcings().create(entry, transform, false, output);
    require(!status && !output &&
                status.code ==
                    (badAllocation ? WVKernelStatusCode::allocationFailure
                                   : WVKernelStatusCode::invalidConfiguration),
            "forcing factory failure was not transactional/normalized");
  }

  {
    auto registration =
        wavevortex::runtime::test::linearCoefficientRegistration();
    registration.matlabClassName = "StatusForcing";
    registration.factory =
        [](const WVFrozenForcingEntry &entry,
           const WVTransformConstantStratificationDescriptor &, bool,
           std::unique_ptr<WVForcing> &output) {
      output = std::make_unique<wavevortex::runtime::test::
                                    WVTestPortableLinearCoefficientForcing>(
          entry, 1.0);
      return WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                            "forcing status failure"};
    };
    WVExtensionCatalogBuilder builder;
    require(static_cast<bool>(builder.addForcingFactory(
                std::move(registration))),
            "status-failing forcing registration failed");
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(static_cast<bool>(builder.freeze(catalog)),
            "status-failing forcing catalog freeze failed");
    WVFrozenForcingEntry entry;
    entry.typeIdentifier = "StatusForcing";
    std::unique_ptr<WVForcing> output =
        std::make_unique<wavevortex::runtime::test::
                             WVTestPortableLinearCoefficientForcing>(entry, 1.0);
    const auto status =
        catalog->forcings().create(entry, transform, false, output);
    require(!status && !output,
            "forcing status failure published a partial implementation");
  }

  return 0;
}
