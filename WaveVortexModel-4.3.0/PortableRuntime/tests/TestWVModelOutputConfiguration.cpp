#include "WaveVortexRuntime/WVModelOutputConfiguration.hpp"
#include "WVTestExtensionCatalog.hpp"
#include "WVTestQuadraticSchedule.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

struct TemporaryDirectory {
  std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("wv-output-configuration-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  TemporaryDirectory() { std::filesystem::create_directories(path); }
  ~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
};

WVCheckpoint checkpointTemplate() {
  WVCheckpoint checkpoint;
  const auto path = std::filesystem::path(WV_RUNTIME_FIXTURE_DIR) /
                    "forcing-nonlinear.nc";
  const auto status = WVCheckpointReader::read(path.string(), *test::extensionCatalog(), checkpoint);
  require(static_cast<bool>(status), status.message);
  return checkpoint;
}

std::vector<char> fileBytes(const std::filesystem::path &path) {
  std::ifstream stream(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(stream),
          std::istreambuf_iterator<char>()};
}

class FailingSampleSource final : public WVObserverSampleSource {
public:
  WVKernelStatus observationSchema(
      const WVObserverRecord &, WVObservationSchema &) override {
    return {WVKernelStatusCode::allocationFailure,
            "injected output-schema failure"};
  }
  WVKernelStatus prepare(const WVOutputEvent &event) override {
    scheduledTime_ = event.scheduledTime;
    return WVKernelStatus::ok();
  }
  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &route, const WVOutputObserverView &observer,
      WVObservationOccurrenceIdentity &output) const override {
    if (route.schedulePayload == nullptr)
      return {WVKernelStatusCode::invalidConfiguration,
              "A test observation route has no schedule payload."};
    output = {};
    output.observerOrdinal = observer.observerOrdinal;
    output.semanticScheduleOrdinal = route.semanticScheduleOrdinal;
    output.scheduleOrdinal = route.scheduleOrdinal;
    output.scheduledTime = scheduledTime_;
    output.scheduleCursorIdentity = route.scheduleCursorIdentity;
    output.payloadFingerprint = route.schedulePayload->valueFingerprint();
    return WVKernelStatus::ok();
  }
  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &, const WVObserverRecord &,
      WVObservationBatch &) override {
    return {WVKernelStatusCode::allocationFailure,
            "injected observation-batch failure"};
  }
private:
  double scheduledTime_ = 0.0;
};

class PreflightRejectingSampleSource final : public WVObserverSampleSource {
public:
  WVKernelStatus preflight(const WVOutputPlan &) override {
    ++preflightCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "incompatible schedule/observer occurrence payload schemas"};
  }
  WVKernelStatus observationSchema(
      const WVObserverRecord &, WVObservationSchema &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "schema discovery ran after rejected source preflight"};
  }
  WVKernelStatus prepare(const WVOutputEvent &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "event preparation ran after rejected source preflight"};
  }
  WVKernelStatus preparedOccurrenceIdentity(
      const WVOutputRouteView &, const WVOutputObserverView &,
      WVObservationOccurrenceIdentity &) const override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "occurrence lookup ran after rejected source preflight"};
  }
  WVKernelStatus observationBatch(
      const WVObservationOccurrenceIdentity &, const WVObserverRecord &,
      WVObservationBatch &) override {
    ++unexpectedCallCount_;
    return {WVKernelStatusCode::invalidConfiguration,
            "batch construction ran after rejected source preflight"};
  }

  std::size_t preflightCount() const noexcept { return preflightCount_; }
  std::size_t unexpectedCallCount() const noexcept {
    return unexpectedCallCount_;
  }

private:
  std::size_t preflightCount_ = 0;
  mutable std::size_t unexpectedCallCount_ = 0;
};

WVPortableObserverRecord allObserverRecord(const WVCheckpoint &checkpoint) {
  WVPortableObserverRecord record;
  const std::vector<std::size_t> coefficientShape = {
      checkpoint.state.coefficients.shape.rows,
      checkpoint.state.coefficients.shape.columns};
  for (const auto *name : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {name, WVStateScalarType::complex64, coefficientShape,
         WVToleranceKind::coefficientEnergyScaled, 1e-6,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  for (const auto *name : {"particles-x", "particles-y"})
    record.stateBlocks.push_back(
        {name, WVStateScalarType::real64, {2},
         WVToleranceKind::uniformAbsolute, 1e-4,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  record.stateBlocks.push_back(
      {"tracer-state", WVStateScalarType::real64,
       {checkpoint.configuration.Nx, checkpoint.configuration.Ny,
        checkpoint.configuration.Nz},
       WVToleranceKind::uniformAbsolute, 1e-5,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState});

  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  WVObserverRecord fields;
  fields.identifier = "fields";
  fields.name = "WVEulerianFields";
  fields.typeIdentifier = "WVEulerianFields";
  fields.fieldNames = {"u", "v", "rho_e"};
  WVObserverRecord mooring;
  mooring.identifier = "mooring";
  mooring.name = "mooring";
  mooring.typeIdentifier = "WVMooring";
  mooring.fieldNames = {"u"};
  mooring.x = {0.0};
  mooring.y = {0.0};
  WVObserverRecord particles;
  particles.identifier = "particles";
  particles.name = "particles";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.stateBlockIdentifiers = {"particles-x", "particles-y"};
  particles.x = {0.1, 0.2};
  particles.y = {0.3, 0.4};
  particles.z = {-100.0, -200.0};
  particles.isXYOnly = true;
  particles.horizontalAbsoluteTolerance = 1e-4;
  WVObserverRecord tracer;
  tracer.identifier = "tracer";
  tracer.name = "tracer";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"tracer-state"};
  record.observers = {coefficients, fields, mooring, particles, tracer};
  return record;
}

WVPortableObserverRecord coefficientRecord(const WVCheckpoint &checkpoint) {
  auto record = allObserverRecord(checkpoint);
  record.stateBlocks.resize(3);
  record.observers.resize(1);
  return record;
}

bool sameOutputGraph(const std::vector<WVOutputFileRecord> &left,
                     const std::vector<WVOutputFileRecord> &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t file = 0; file < left.size(); ++file) {
    if (left[file].identifier != right[file].identifier ||
        left[file].destination != right[file].destination ||
        left[file].groups.size() != right[file].groups.size())
      return false;
    for (std::size_t group = 0; group < left[file].groups.size(); ++group) {
      const auto &a = left[file].groups[group];
      const auto &b = right[file].groups[group];
      if (a.identifier != b.identifier || a.name != b.name ||
          a.schedule.typeIdentifier != b.schedule.typeIdentifier ||
          a.schedule.contractVersion != b.schedule.contractVersion ||
          a.schedule.configuration.schemaIdentifier !=
              b.schedule.configuration.schemaIdentifier ||
          a.schedule.configuration.schemaVersion !=
              b.schedule.configuration.schemaVersion ||
          a.schedule.configuration.values.size() !=
              b.schedule.configuration.values.size() ||
          a.schedule.outputInterval != b.schedule.outputInterval ||
          a.schedule.initialTime != b.schedule.initialTime ||
          a.schedule.finalTime != b.schedule.finalTime ||
          a.observerIdentifiers != b.observerIdentifiers ||
          a.containsCompleteCoefficientRestart !=
              b.containsCompleteCoefficientRestart)
        return false;
      for (std::size_t value = 0;
           value < a.schedule.configuration.values.size(); ++value) {
        const auto &av = a.schedule.configuration.values[value];
        const auto &bv = b.schedule.configuration.values[value];
        if (av.name != bv.name || av.dimensions != bv.dimensions ||
            av.storage != bv.storage)
          return false;
      }
    }
  }
  return true;
}

WVModelOutputFile fileBuilder(const std::filesystem::path &path,
                              const std::string &identifier,
                              double initialTime, double finalTime,
                              bool allObservers) {
  WVModelOutputFile file;
  auto status = WVModelOutputFile::create(path.string(), file, identifier);
  require(static_cast<bool>(status), status.message);
  WVModelOutputGroup *restart = nullptr;
  status = file.addNewEvenlySpacedOutputGroup(
      "wave-vortex", 1.0, initialTime, finalTime, restart, "restart");
  require(static_cast<bool>(status) && restart != nullptr,
          "failed to add restart group");
  for (const auto *observer : {"coefficients", "fields", "mooring",
                               "particles", "tracer"}) {
    if (allObservers || std::string(observer) == "coefficients") {
      status = restart->addObservingSystem(observer);
      require(static_cast<bool>(status), status.message);
    }
  }
  status = restart->containsCompleteCoefficientRestart(true);
  require(static_cast<bool>(status), status.message);
  return file;
}

void deliverFirstEvent(WVModelOutputNetCDFSink &sink,
                       const WVOutputPlan &plan,
                       const WVCheckpoint &checkpoint) {
  auto status = sink.preflight(plan);
  require(static_cast<bool>(status), status.message);
  require(plan.eventCount() != 0, "test output plan has no event");
  const auto planned = plan.event(0);
  WVOutputEvent event;
  event.eventOrdinal = planned.eventOrdinal;
  event.scheduledTime = planned.scheduledTime;
  event.state = {checkpoint.state.view(), nullptr, 0};
  event.routes = planned.routes;
  event.routeCount = planned.routeCount;
  for (std::size_t route = 0; route < planned.routeCount; ++route) {
    WVOutputDeliveryResult delivery;
    status = sink.deliver(event, planned.routes[route], delivery);
    require(static_cast<bool>(status), status.message);
  }
}

void testStructuralCompilation() {
  TemporaryDirectory directory;
  const auto checkpoint = checkpointTemplate();
  auto record = allObserverRecord(checkpoint);
  std::vector<WVModelOutputFile> files;
  files.push_back(fileBuilder(directory.path / "first.nc", "first",
                              checkpoint.state.t, checkpoint.state.t + 2.0,
                              true));
  files.push_back(fileBuilder(directory.path / "second.nc", "second",
                              checkpoint.state.t, checkpoint.state.t + 2.0,
                              true));
  auto *stable = files.front().outputGroupWithName("wave-vortex");
  require(stable != nullptr, "group lookup failed");
  WVModelOutputGroup extra;
  WVOutputGroupRecord diagnosticRecord{
      "diagnostics", "diagnostics", {},
      {"fields", "mooring", "particles", "tracer"}, false};
  diagnosticRecord.schedule = test::quadraticSchedule(
      checkpoint.state.t + 2.0, checkpoint.state.t, 0.5);
  auto status =
      WVModelOutputGroup::fromRecord(std::move(diagnosticRecord), extra);
  require(static_cast<bool>(status), status.message);
  status = files.front().addOutputGroup(std::move(extra));
  require(static_cast<bool>(status), status.message);
  require(stable == files.front().outputGroupWithName("wave-vortex"),
          "adding a group invalidated a returned group reference");
  require(stable == files.front().outputGroupWithIdentifier("restart"),
          "group identifier lookup failed");
  require(!files.front().addObservingSystem("fields"),
          "file-level addition accepted multiple groups");

  const auto firstPath = std::filesystem::absolute(directory.path / "first.nc")
                             .lexically_normal()
                             .string();
  const auto secondPath =
      std::filesystem::absolute(directory.path / "second.nc")
          .lexically_normal()
          .string();
  const std::vector<std::string> allObservers{
      "coefficients", "fields", "mooring", "particles", "tracer"};
  const std::vector<WVOutputFileRecord> expectedFiles{
      {"first",
       firstPath,
       {{"restart",
         "wave-vortex",
         {1.0, checkpoint.state.t, checkpoint.state.t + 2.0},
         allObservers,
         true},
        {"diagnostics",
         "diagnostics",
         test::quadraticSchedule(checkpoint.state.t + 2.0,
                                 checkpoint.state.t, 0.5),
         {"fields", "mooring", "particles", "tracer"},
         false}}},
      {"second",
       secondPath,
       {{"restart",
         "wave-vortex",
         {1.0, checkpoint.state.t, checkpoint.state.t + 2.0},
         allObservers,
         true}}}};

  WVModelOutputConfiguration configuration;
  status = WVModelOutputConfiguration::build(
      std::move(record), std::move(files), WVModelOutputPolicy::create,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, configuration);
  require(static_cast<bool>(status), status.message);
  require(configuration.descriptor().outputFiles().size() == 2,
          "compiled descriptor lost output files");
  require(sameOutputGraph(configuration.descriptor().outputFiles(),
                          expectedFiles),
          "facade-built and direct low-level output records differ");
  require(configuration.plan().metrics().fileCount == 2 &&
              configuration.plan().metrics().groupCount == 3 &&
              configuration.plan().metrics().distinctObserverCount == 5,
          "compiled plan does not match the builder graph");
  require(configuration.plan().eventCount() == 4,
          "coincident bounded schedule event count changed");
  const std::vector<double> expectedTimes{
      checkpoint.state.t, checkpoint.state.t + 0.5,
      checkpoint.state.t + 1.0, checkpoint.state.t + 2.0};
  const std::vector<std::size_t> expectedRouteCounts{3, 1, 2, 3};
  for (std::size_t event = 0; event < expectedTimes.size(); ++event)
    require(configuration.plan().event(event).scheduledTime ==
                    expectedTimes[event] &&
                configuration.plan().event(event).routeCount ==
                    expectedRouteCounts[event],
            "facade-built plan changed endpoint or coincident routing");
  require(configuration.policy() == WVModelOutputPolicy::create &&
              configuration.persistentBytes() >=
                  configuration.plan().persistentBytes(),
          "compiled configuration accounting is incomplete");
}

void testValidationAndDeterministicIdentifiers() {
  TemporaryDirectory directory;
  const auto checkpoint = checkpointTemplate();
  WVModelOutputGroup unbounded;
  auto unboundedStatus = WVModelOutputGroup::evenlySpaced(
      "unbounded", 1.0, checkpoint.state.t,
      std::numeric_limits<double>::infinity(), unbounded);
  require(static_cast<bool>(unboundedStatus),
          "MATLAB's unbounded evenly spaced schedule was rejected");
  WVModelOutputFile first;
  WVModelOutputFile second;
  auto status = WVModelOutputFile::create(
      (directory.path / "same.nc").string(), first);
  require(static_cast<bool>(status), status.message);
  status = WVModelOutputFile::create(
      (directory.path / "." / "same.nc").string(), second);
  require(static_cast<bool>(status) && first.identifier() == second.identifier(),
          "default file identifiers are not deterministic");
  WVModelOutputGroup *group = nullptr;
  status = first.addNewEvenlySpacedOutputGroup(
      "group", 1.0, checkpoint.state.t, checkpoint.state.t + 1.0, group);
  require(static_cast<bool>(status), status.message);
  require(static_cast<bool>(group->addObservingSystem("coefficients")),
          "observer addition failed");
  require(!group->addObservingSystem("coefficients"),
          "duplicate observer membership was accepted");
  require(static_cast<bool>(group->containsCompleteCoefficientRestart(true)),
          "restart declaration failed");
  const auto generatedGroupIdentifier = group->identifier();
  require(!generatedGroupIdentifier.empty(), "group identifier was not generated");

  auto record = coefficientRecord(checkpoint);
  record.observers.front().identifier = "different";
  std::vector<WVModelOutputFile> files;
  files.push_back(std::move(first));
  WVModelOutputConfiguration rejected;
  status = WVModelOutputConfiguration::build(
      std::move(record), std::move(files), WVModelOutputPolicy::create,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 1.0, rejected);
  require(!status, "an unresolved observing system was accepted");
}

void testCreateReplaceAndAppendPolicies() {
  TemporaryDirectory directory;
  auto checkpoint = checkpointTemplate();
  const auto rejectedCreatePath = directory.path / "preflight-create.nc";
  auto record = coefficientRecord(checkpoint);
  std::vector<WVModelOutputFile> files;
  files.push_back(fileBuilder(rejectedCreatePath, "preflight-create",
                              checkpoint.state.t,
                              checkpoint.state.t + 2.0, false));
  WVModelOutputConfiguration rejectedCreate;
  auto status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::create,
      test::extensionCatalog(), checkpoint.state.t,
      checkpoint.state.t + 2.0, rejectedCreate);
  require(static_cast<bool>(status), status.message);
  WVIntegrationStateLayout rejectedCreateLayout;
  status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, rejectedCreate.descriptor(),
      rejectedCreateLayout);
  require(static_cast<bool>(status), status.message);
  PreflightRejectingSampleSource rejectedCreateSource;
  WVModelOutputNetCDFSink sink;
  auto persistence = rejectedCreate.openNetCDFSink(
      {test::extensionCatalog(), checkpoint, false}, rejectedCreateLayout,
      &rejectedCreateSource, sink);
  require(!persistence &&
              persistence.code == WVCheckpointStatusCode::unsupportedObserver &&
              rejectedCreateSource.preflightCount() == 1 &&
              rejectedCreateSource.unexpectedCallCount() == 0 &&
              !std::filesystem::exists(rejectedCreatePath),
          "payload-schema preflight rejection created an output destination");

  const auto path = directory.path / "policy.nc";
  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 2.0, false));
  WVModelOutputConfiguration create;
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::create,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, create);
  require(static_cast<bool>(status), status.message);
  WVIntegrationStateLayout layout;
  status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, create.descriptor(), layout);
  require(static_cast<bool>(status), status.message);
  persistence = create.openNetCDFSink({test::extensionCatalog(), checkpoint, false}, layout,
                                      nullptr, sink);
  require(static_cast<bool>(persistence), persistence.message);
  deliverFirstEvent(sink, create.plan(), checkpoint);
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 3.0, false));
  WVModelOutputConfiguration rejectedAppend;
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::append,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 3.0, rejectedAppend);
  require(!status, "append accepted a changed schedule graph");

  files.clear();
  auto wrongProgress = fileBuilder(path, "policy", checkpoint.state.t,
                                   checkpoint.state.t + 2.0, false);
  require(static_cast<bool>(
              wrongProgress.outputGroupWithName("wave-vortex")
                  ->scheduleContinuation({1, {}})),
          "explicit append progress setup failed");
  files.push_back(std::move(wrongProgress));
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::append,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, rejectedAppend);
  require(!status, "append accepted mismatched committed progress");

  const auto original = fileBytes(path);
  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 2.0, true));
  WVModelOutputConfiguration failedReplacement;
  status = WVModelOutputConfiguration::build(
      allObserverRecord(checkpoint), std::move(files),
      WVModelOutputPolicy::replace, test::extensionCatalog(), checkpoint.state.t,
      checkpoint.state.t + 2.0, failedReplacement);
  require(static_cast<bool>(status), status.message);
  WVIntegrationStateLayout allLayout;
  status = WVIntegrationStateLayout::create(
      checkpoint.state.coefficients.shape, failedReplacement.descriptor(),
      allLayout);
  require(static_cast<bool>(status), status.message);
  PreflightRejectingSampleSource rejectedReplacementSource;
  persistence = failedReplacement.openNetCDFSink(
      {test::extensionCatalog(), checkpoint, false}, allLayout,
      &rejectedReplacementSource, sink);
  require(!persistence && rejectedReplacementSource.preflightCount() == 1 &&
              rejectedReplacementSource.unexpectedCallCount() == 0 &&
              fileBytes(path) == original,
          "payload-schema preflight rejection changed a replacement destination");
  FailingSampleSource failure;
  persistence = failedReplacement.openNetCDFSink(
      {test::extensionCatalog(), checkpoint, false}, allLayout, &failure, sink);
  require(!persistence && fileBytes(path) == original,
          "failed replacement changed the original destination");

  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 2.0, false));
  WVModelOutputConfiguration replace;
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::replace,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, replace);
  require(static_cast<bool>(status), status.message);
  persistence = replace.openNetCDFSink({test::extensionCatalog(), checkpoint, false}, layout, nullptr,
                                       sink);
  require(static_cast<bool>(persistence), persistence.message);
  deliverFirstEvent(sink, replace.plan(), checkpoint);
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 2.0, false));
  WVModelOutputConfiguration append;
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::append,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, append);
  require(static_cast<bool>(status), status.message);
  require(append.scheduleContinuations().size() == 1 &&
              append.scheduleContinuations()
                      .front()
                      .cursor.committedOrdinal == 0 &&
              append.destinationProgress().size() == 1 &&
              append.destinationProgress()
                      .front()
                      .committedScheduleCursor.committedOrdinal == 0,
          "append progress was not recovered");
  const auto beforeRejectedAppend = fileBytes(path);
  PreflightRejectingSampleSource rejectedAppendSource;
  persistence = append.openNetCDFSink(
      {test::extensionCatalog(), checkpoint, false}, layout,
      &rejectedAppendSource, sink);
  require(!persistence && rejectedAppendSource.preflightCount() == 1 &&
              rejectedAppendSource.unexpectedCallCount() == 0 &&
              fileBytes(path) == beforeRejectedAppend,
          "payload-schema preflight rejection changed an append destination");
  persistence = append.openNetCDFSink({test::extensionCatalog(), checkpoint, false}, layout, nullptr,
                                      sink);
  require(static_cast<bool>(persistence), persistence.message);
  persistence = sink.close();
  require(static_cast<bool>(persistence), persistence.message);

  files.clear();
  files.push_back(fileBuilder(path, "policy", checkpoint.state.t,
                              checkpoint.state.t + 2.0, false));
  WVModelOutputConfiguration collision;
  status = WVModelOutputConfiguration::build(
      record, std::move(files), WVModelOutputPolicy::create,
      test::extensionCatalog(), checkpoint.state.t, checkpoint.state.t + 2.0, collision);
  require(!status, "create policy accepted an existing destination");
}

} // namespace

int main() {
  try {
    testStructuralCompilation();
    testValidationAndDeterministicIdentifiers();
    testCreateReplaceAndAppendPolicies();
    std::cout << "WVModelOutputConfiguration tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
