#include "WVTestOccurrenceSchedule.hpp"
#include "WaveVortexRuntime/WVOutputOrchestration.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::runtime::test;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

std::shared_ptr<const WVExtensionCatalog> occurrenceCatalog() {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (status)
    status = registerOccurrenceSchedule(builder);
  std::shared_ptr<const WVExtensionCatalog> catalog;
  if (status)
    status = builder.freeze(catalog);
  require(static_cast<bool>(status),
          "occurrence schedule catalog: " + status.message);
  return catalog;
}

WVPortableObserverRecord occurrenceRecord() {
  WVPortableObserverRecord record;
  for (const auto *identifier : {"Ap", "Am", "A0"})
    record.stateBlocks.push_back(
        {identifier,
         WVStateScalarType::complex64,
         {1, 1},
         WVToleranceKind::coefficientEnergyScaled,
         0.0,
         WVStateOwnership::integratorOwned,
         WVRestartRequirement::requiredDynamicState});
  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "Wave-vortex coefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"Ap", "Am", "A0"};
  record.observers.push_back(std::move(coefficients));

  const auto outputFile = [](std::string fileIdentifier,
                             std::string destination,
                             std::string groupIdentifier,
                             WVOutputScheduleRecord schedule) {
    WVOutputGroupRecord group;
    group.identifier = std::move(groupIdentifier);
    group.name = "Occurrence group";
    group.schedule = std::move(schedule);
    group.observerIdentifiers = {"coefficients"};
    group.containsCompleteCoefficientRestart = true;
    WVOutputFileRecord file;
    file.identifier = std::move(fileIdentifier);
    file.destination = std::move(destination);
    file.groups.push_back(std::move(group));
    return file;
  };

  record.outputFiles.push_back(outputFile(
      "first", "first.nc", "shared", occurrenceSchedule(1.25, 10, 1, 0.0)));
  record.outputFiles.push_back(outputFile(
      "second", "second.nc", "shared", occurrenceSchedule(1.25, 10, 1, 0.0)));
  record.outputFiles.push_back(outputFile(
      "differentConfig", "different-config.nc", "shared",
      occurrenceSchedule(9.5, 10, 1, 0.0)));
  record.outputFiles.push_back(outputFile(
      "differentGroup", "different-group.nc", "independent",
      occurrenceSchedule(1.25, 10, 1, 0.0)));
  return record;
}

void testResolvedPayloadSchemaAndValidation() {
  WVOutputSchedulePayloadSchema schema;
  auto status = createOccurrencePayloadSchema(schema);
  require(static_cast<bool>(status) &&
              schema.identifier() == occurrencePayloadSchemaIdentifier &&
              schema.version() == 1 && schema.slotCount() == 3 &&
              schema.payloadBytes() == 17,
          "mixed occurrence-payload schema resolves once");
  require(schema.slots()[0].name == "sourceReal" &&
              schema.slots()[0].scalarType ==
                  WVOutputSchedulePayloadType::real64 &&
              schema.slots()[0].elementCount == 1 &&
              schema.slots()[0].byteOffset == 0 &&
              schema.slots()[1].name == "sourceInteger" &&
              schema.slots()[1].scalarType ==
                  WVOutputSchedulePayloadType::integer64 &&
              schema.slots()[1].byteOffset == sizeof(double) &&
              schema.slots()[2].name == "sourceBoolean" &&
              schema.slots()[2].scalarType ==
                  WVOutputSchedulePayloadType::boolean8 &&
              schema.slots()[2].byteOffset == 2 * sizeof(double),
          "payload names, types, extents, and aligned slots are resolved");

  WVOutputSchedulePayload payload;
  require(static_cast<bool>(payload.reset(schema)),
          "mixed occurrence payload reset");
  const double realValue = 4.25;
  const std::int64_t integerValue = -17;
  const std::uint8_t booleanValue = 1;
  require(static_cast<bool>(payload.setReal(schema, 0, &realValue, 1)) &&
              static_cast<bool>(
                  payload.setInteger(schema, 1, &integerValue, 1)) &&
              static_cast<bool>(
                  payload.setBoolean(schema, 2, &booleanValue, 1)),
          "finite numeric, integer, and Boolean payload values are accepted");
  WVOutputSchedulePayloadRealView realView;
  WVOutputSchedulePayloadIntegerView integerView;
  WVOutputSchedulePayloadBooleanView booleanView;
  require(static_cast<bool>(payload.real(schema, 0, realView)) &&
              static_cast<bool>(payload.integer(schema, 1, integerView)) &&
              static_cast<bool>(payload.boolean(schema, 2, booleanView)) &&
              realView.count == 1 && realView.data[0] == realValue &&
              integerView.count == 1 && integerView.data[0] == integerValue &&
              booleanView.count == 1 &&
              booleanView.data[0] == booleanValue,
          "resolved payload slots recover exact typed values");

  WVOutputSchedulePayload same;
  require(static_cast<bool>(same.reset(schema)) &&
              static_cast<bool>(same.setReal(schema, 0, &realValue, 1)) &&
              static_cast<bool>(
                  same.setInteger(schema, 1, &integerValue, 1)) &&
              static_cast<bool>(
                  same.setBoolean(schema, 2, &booleanValue, 1)) &&
              payload.sameValue(same) &&
              payload.valueFingerprint() == same.valueFingerprint(),
          "equal typed payloads have exact value identity");
  const double otherReal = realValue + 1.0;
  require(static_cast<bool>(same.setReal(schema, 0, &otherReal, 1)) &&
              !payload.sameValue(same) &&
              payload.valueFingerprint() != same.valueFingerprint(),
          "different payload values remain isolated");

  const auto nan = std::numeric_limits<double>::quiet_NaN();
  require(payload.setReal(schema, 0, &nan, 1).code ==
                  WVKernelStatusCode::invalidConfiguration &&
              payload.setInteger(schema, 0, &integerValue, 1).code ==
                  WVKernelStatusCode::invalidConfiguration,
          "nonfinite and wrong-type payload writes are rejected");
  const std::uint8_t invalidBoolean = 2;
  require(payload.setBoolean(schema, 2, &invalidBoolean, 1).code ==
                  WVKernelStatusCode::invalidConfiguration &&
              payload.boolean(schema, 0, booleanView).code ==
                  WVKernelStatusCode::invalidConfiguration,
          "invalid Boolean values and wrong-type payload reads are rejected");

  WVOutputSchedulePayloadSchema exactSchema;
  status = WVOutputSchedulePayloadSchema::create(
      "exact-4k-payload-v1", 1,
      {{"values", WVOutputSchedulePayloadType::real64, {512}}}, exactSchema);
  WVOutputSchedulePayload exactPayload;
  std::vector<double> exactValues(512, 1.0);
  require(static_cast<bool>(status) &&
              exactSchema.payloadBytes() ==
                  WVMaximumOutputSchedulePayloadBytes &&
              static_cast<bool>(exactPayload.reset(exactSchema)) &&
              static_cast<bool>(exactPayload.setReal(
                  exactSchema, 0, exactValues.data(), exactValues.size())),
          "an exactly 4 KiB encoded occurrence payload is accepted");

  WVOutputSchedulePayloadSchema oversizedSchema;
  status = WVOutputSchedulePayloadSchema::create(
      "oversized-payload-v1", 1,
      {{"values", WVOutputSchedulePayloadType::real64, {513}}},
      oversizedSchema);
  require(status.code == WVKernelStatusCode::sizeOverflow,
          "an occurrence payload over 4 KiB is rejected");
}

void testSourceLinkedSchedule() {
  resetOccurrenceScheduleCounters();
  auto record = occurrenceSchedule(2.5, 40, 1, 2.0, 0.0, 1.0);
  WVKernelStatus status;
  const auto schedule = makeOccurrenceSchedule(record, status);
  require(static_cast<bool>(status) && schedule != nullptr &&
              occurrenceScheduleCounters().constructionCount == 1 &&
              occurrenceScheduleCounters().peekCount == 0,
          "source-linked occurrence schedule is resolved at construction");

  WVOutputScheduleCursor cursor;
  WVOutputScheduleOccurrence first;
  bool available = false;
  status = schedule->peek(cursor, 0.0, 2.0, first, available);
  WVOutputSchedulePayloadRealView realView;
  WVOutputSchedulePayloadIntegerView integerView;
  WVOutputSchedulePayloadBooleanView booleanView;
  require(static_cast<bool>(status) && available &&
              first.scheduledTime == 0.0 && first.ordinal == 0 &&
              first.proposedCursor.committedOrdinal == 0 &&
              static_cast<bool>(first.payload.real(
                  schedule->payloadSchema(), 0, realView)) &&
              static_cast<bool>(first.payload.integer(
                  schedule->payloadSchema(), 1, integerView)) &&
              static_cast<bool>(first.payload.boolean(
                  schedule->payloadSchema(), 2, booleanView)) &&
              realView.data[0] == 2.5 && integerView.data[0] == 40 &&
              booleanView.data[0] == 1,
          "source values are encoded in deterministic resolved slots");

  WVOutputScheduleOccurrence repeated;
  status = schedule->peek(cursor, 0.0, 2.0, repeated, available);
  require(static_cast<bool>(status) && available &&
              repeated.cursorIdentity == first.cursorIdentity &&
              repeated.payload.sameValue(first.payload) &&
              occurrenceScheduleCounters().peekCount == 2,
          "repeated peeks return exact payload and cursor identity");

  WVOutputScheduleOccurrence second;
  status = schedule->peek(first.proposedCursor, 0.0, 2.0, second, available);
  require(static_cast<bool>(status) && available && second.ordinal == 1 &&
              second.scheduledTime == 1.0 &&
              second.cursorIdentity != first.cursorIdentity &&
              !second.payload.sameValue(first.payload),
          "successive occurrences carry distinct deterministic source data");

  auto invalidRecord = record;
  invalidRecord.configuration.values[0].storage =
      std::vector<double>{std::numeric_limits<double>::infinity()};
  require(makeOccurrenceSchedule(invalidRecord, status) == nullptr &&
              !status &&
              occurrenceScheduleCounters().constructionCount == 1,
          "nonfinite source configuration is rejected before construction");
  invalidRecord = record;
  invalidRecord.configuration.values[1].storage =
      std::vector<double>{40.0};
  require(makeOccurrenceSchedule(invalidRecord, status) == nullptr && !status,
          "source configuration with the wrong scalar type is rejected");
  invalidRecord = record;
  invalidRecord.configuration.values[2].storage =
      std::vector<std::uint8_t>{2};
  require(makeOccurrenceSchedule(invalidRecord, status) == nullptr && !status,
          "source configuration with an invalid Boolean is rejected");
}

void testLogicalScheduleAndPayloadIsolation() {
  resetOccurrenceScheduleCounters();
  const auto catalog = occurrenceCatalog();
  const auto record = occurrenceRecord();
  WVPortableObserverDescriptor descriptor;
  auto status =
      WVPortableObserverDescriptor::create(record, catalog, descriptor);
  require(static_cast<bool>(status) &&
              occurrenceScheduleCounters().constructionCount == 0,
          "descriptor validation does not construct schedule providers");
  WVOutputPlan plan;
  status = WVOutputPlan::create(descriptor, catalog, 0.0, 0.0, {}, plan);
  require(static_cast<bool>(status) && plan.groupCount() == 4 &&
              occurrenceScheduleCounters().constructionCount == 4 &&
              occurrenceScheduleCounters().peekCount == 0,
          "output planning constructs each source-linked route once");

  const auto firstGroup = plan.groupRoute(0);
  const auto secondGroup = plan.groupRoute(1);
  const auto configGroup = plan.groupRoute(2);
  const auto independentGroup = plan.groupRoute(3);
  require(firstGroup.fileOrdinal != secondGroup.fileOrdinal &&
              firstGroup.semanticScheduleOrdinal ==
                  secondGroup.semanticScheduleOrdinal &&
              firstGroup.schedulePayloadSchema != nullptr &&
              secondGroup.schedulePayloadSchema != nullptr &&
              sameOutputSchedulePayloadSchema(
                  *firstGroup.schedulePayloadSchema,
                  *secondGroup.schedulePayloadSchema),
          "compatible logical groups share identity across destinations");
  require(configGroup.semanticScheduleOrdinal !=
                  firstGroup.semanticScheduleOrdinal &&
              independentGroup.semanticScheduleOrdinal !=
                  firstGroup.semanticScheduleOrdinal &&
              configGroup.semanticScheduleOrdinal !=
                  independentGroup.semanticScheduleOrdinal,
          "different configuration or group identity remains isolated");

  const auto event = plan.event(0);
  require(event.routeCount == 4 && event.scheduledTime == 0.0 &&
              occurrenceScheduleCounters().peekCount == 4,
          "one construction-resolved peek per coincident route generates the event");
  const auto &first = event.routes[0];
  const auto &second = event.routes[1];
  const auto &differentConfig = event.routes[2];
  const auto &differentGroup = event.routes[3];
  require(first.schedulePayload != nullptr && second.schedulePayload != nullptr &&
              differentConfig.schedulePayload != nullptr &&
              differentGroup.schedulePayload != nullptr &&
              first.schedulePayload->sameValue(*second.schedulePayload) &&
              first.scheduleCursorIdentity == second.scheduleCursorIdentity &&
              first.semanticScheduleOrdinal ==
                  second.semanticScheduleOrdinal,
          "compatible coincident destinations share complete semantic occurrence data");
  require(!first.schedulePayload->sameValue(
                  *differentConfig.schedulePayload) &&
              first.scheduleCursorIdentity !=
                  differentConfig.scheduleCursorIdentity &&
              first.semanticScheduleOrdinal !=
                  differentConfig.semanticScheduleOrdinal,
          "different configuration and payload cannot share an occurrence");
  require(first.schedulePayload->sameValue(*differentGroup.schedulePayload) &&
              first.scheduleCursorIdentity ==
                  differentGroup.scheduleCursorIdentity &&
              first.semanticScheduleOrdinal !=
                  differentGroup.semanticScheduleOrdinal,
          "equal payloads from different logical groups remain isolated");
}

} // namespace

int main() {
  testResolvedPayloadSchemaAndValidation();
  testSourceLinkedSchedule();
  testLogicalScheduleAndPayloadIsolation();
  std::cout << "PASS: occurrence payload and semantic schedule identity\n";
  return 0;
}
