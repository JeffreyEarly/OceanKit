#include "WaveVortexRuntime/WVObservation.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace wavevortex::runtime;

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

WVObservationSchema schema() {
  WVObservationSchema result;
  result.identifier = "test-observation-v1";
  result.axes = {{"profile", "profile", WVObservationAxisKind::unlimited, 0,
                  WVObservationCoordinateRole::profile},
                 {"sample", "sample", WVObservationAxisKind::unlimited, 0,
                  WVObservationCoordinateRole::identifier},
                 {"bin", "depth_bin", WVObservationAxisKind::fixed, 2,
                  WVObservationCoordinateRole::depth}};
  result.variables = {
      {"profile-count", "profile_sample_count",
       WVObservationScalarType::integer64, {"profile"},
       WVObservationValueLayout::flat, "1", "samples in each profile", {},
       WVObservationCoordinateRole::none,
       WVObservationRaggedRole::rowCount, "sample"},
      {"sample-time", "sample_time", WVObservationScalarType::real64,
       {"sample"}, WVObservationValueLayout::flat, "s", "sample time", {},
       WVObservationCoordinateRole::sampleTime,
       WVObservationRaggedRole::none, {}},
      {"sample-x", "sample_x", WVObservationScalarType::real64, {"sample"},
       WVObservationValueLayout::flat, "m", "sample x", {},
       WVObservationCoordinateRole::x, WVObservationRaggedRole::none, {}},
      {"sample-bin", "sample_bin", WVObservationScalarType::complex64,
       {"bin", "sample"}, WVObservationValueLayout::flat, "m s-1",
       "complex binned sample", {}, WVObservationCoordinateRole::none,
       WVObservationRaggedRole::none, {}},
      {"quality", "quality", WVObservationScalarType::boolean8, {"sample"},
       WVObservationValueLayout::flat, "1", "quality flag", {},
       WVObservationCoordinateRole::none, WVObservationRaggedRole::none, {}},
      {"label", "profile_label", WVObservationScalarType::text, {"profile"},
       WVObservationValueLayout::flat, "", "profile label", {},
       WVObservationCoordinateRole::none, WVObservationRaggedRole::none, {}},
      {"event-id", "event_id", WVObservationScalarType::integer64, {},
       WVObservationValueLayout::record, "1", "event identifier", {},
       WVObservationCoordinateRole::none, WVObservationRaggedRole::none, {}}};
  return result;
}

WVObservationBatch batch() {
  static const double times[] = {1.0, 2.0, 3.0};
  static const wavevortex::WVComplex64 bins[] = {
      {1.0, -1.0}, {2.0, -2.0}, {3.0, -3.0},
      {4.0, -4.0}, {5.0, -5.0}, {6.0, -6.0}};
  WVObservationBatch result;
  result.schemaIdentifier = "test-observation-v1";
  result.values.push_back(WVObservationValue::ownInteger(
      "profile-count", {2}, std::vector<std::int64_t>{1, 2}));
  result.values.push_back(
      WVObservationValue::borrowReal("sample-time", {3}, times));
  result.values.push_back(WVObservationValue::ownReal(
      "sample-x", {3}, std::vector<double>{4.0, 5.0, 6.0}));
  result.values.push_back(
      WVObservationValue::borrowComplex("sample-bin", {2, 3}, bins));
  result.values.push_back(WVObservationValue::ownBoolean(
      "quality", {3}, std::vector<std::uint8_t>{1, 0, 1}));
  result.values.push_back(WVObservationValue::ownText(
      "label", {2}, std::vector<std::string>{"pass-a", "pass-b"}));
  result.values.push_back(WVObservationValue::ownInteger(
      "event-id", {}, std::vector<std::int64_t>{7}));
  return result;
}

WVObservationSchema nestedRaggedSchema() {
  WVObservationSchema result;
  result.identifier = "nested-ragged-observation-v1";
  result.axes = {
      {"pass", "pass", WVObservationAxisKind::unlimited, 0,
       WVObservationCoordinateRole::pass},
      {"profile", "profile", WVObservationAxisKind::unlimited, 0,
       WVObservationCoordinateRole::profile},
      {"sample", "sample", WVObservationAxisKind::unlimited, 0,
       WVObservationCoordinateRole::identifier},
      {"depth", "depth", WVObservationAxisKind::fixed, 2,
       WVObservationCoordinateRole::depth}};
  result.variables = {
      {"pass-profile-count", "pass_profile_count",
       WVObservationScalarType::integer64, {"pass"},
       WVObservationValueLayout::flat, "1", "profiles in each pass", {},
       WVObservationCoordinateRole::none,
       WVObservationRaggedRole::rowCount, "profile"},
      {"profile-sample-offset", "profile_sample_offset",
       WVObservationScalarType::integer64, {"profile"},
       WVObservationValueLayout::flat, "1", "first sample in each profile",
       {}, WVObservationCoordinateRole::none,
       WVObservationRaggedRole::rowOffset, "sample"},
      {"sample-value", "sample_value", WVObservationScalarType::real64,
       {"depth", "sample"}, WVObservationValueLayout::flat, "1",
       "sample value by fixed depth", {}, WVObservationCoordinateRole::none,
       WVObservationRaggedRole::none, {}}};
  return result;
}

WVObservationBatch nestedRaggedBatch() {
  WVObservationBatch result;
  result.schemaIdentifier = "nested-ragged-observation-v1";
  result.values.push_back(WVObservationValue::ownInteger(
      "pass-profile-count", {2}, std::vector<std::int64_t>{2, 1}));
  result.values.push_back(WVObservationValue::ownInteger(
      "profile-sample-offset", {3}, std::vector<std::int64_t>{0, 2, 2}));
  result.values.push_back(WVObservationValue::ownReal(
      "sample-value", {2, 4},
      std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}));
  return result;
}

} // namespace

int main() {
  const auto declaration = schema();
  require(static_cast<bool>(validateObservationSchema(declaration)),
          "valid observation schema was rejected");
  auto values = batch();
  require(static_cast<bool>(validateObservationBatch(declaration, values)),
          "valid mixed typed/ragged observation batch was rejected");
  require(values.values[1].ownership ==
              WVObservationBufferOwnership::borrowed &&
              values.values[2].ownership == WVObservationBufferOwnership::owned,
          "observation buffer ownership was not explicit");
  require(values.metrics().liveBytes > 0 &&
              values.metrics().retainedStorageBytes > 0,
          "observation batch storage metrics were not reported");
  std::vector<std::uint8_t> manifest;
  require(static_cast<bool>(
              encodeObservationSchemaManifest(declaration, manifest)) &&
              !manifest.empty(),
          "observation schema manifest was not encoded");
  WVObservationSchema decoded;
  require(static_cast<bool>(
              decodeObservationSchemaManifest(manifest, decoded)) &&
              decoded.identifier == declaration.identifier &&
              decoded.axes.size() == declaration.axes.size() &&
              decoded.variables.size() == declaration.variables.size() &&
              decoded.variables[0].raggedRole ==
                  WVObservationRaggedRole::rowCount &&
              decoded.variables[0].raggedChildAxisIdentifier == "sample" &&
              decoded.variables[3].scalarType ==
                  WVObservationScalarType::complex64,
          "observation schema manifest did not preserve its contract");

  auto zero = batch();
  zero.values[0] = WVObservationValue::ownInteger("profile-count", {0}, {});
  zero.values[1] = WVObservationValue::borrowReal("sample-time", {0}, nullptr);
  zero.values[2] = WVObservationValue::ownReal("sample-x", {0}, {});
  zero.values[3] =
      WVObservationValue::borrowComplex("sample-bin", {2, 0}, nullptr);
  zero.values[4] = WVObservationValue::ownBoolean("quality", {0}, {});
  zero.values[5] = WVObservationValue::ownText("label", {0}, {});
  require(static_cast<bool>(validateObservationBatch(declaration, zero)),
          "zero-length observation occurrence was rejected");

  auto inconsistent = batch();
  inconsistent.values[2].extents = {2};
  require(!validateObservationBatch(declaration, inconsistent),
          "inconsistent unlimited extents were accepted");
  auto undeclared = batch();
  undeclared.values.back().variableIdentifier = "undeclared";
  require(!validateObservationBatch(declaration, undeclared),
          "undeclared observation variable was accepted");
  auto malformed = batch();
  malformed.values[0].ownedInteger64 = {2, 2};
  require(!validateObservationBatch(declaration, malformed),
          "malformed ragged row counts were accepted");
  auto drifted = batch();
  drifted.schemaVersion = 2;
  require(!validateObservationBatch(declaration, drifted),
          "observation schema drift was accepted");
  auto badBoolean = batch();
  badBoolean.values[4].ownedBoolean8[1] = 2;
  require(!validateObservationBatch(declaration, badBoolean),
          "invalid Boolean observation value was accepted");

  const auto nestedDeclaration = nestedRaggedSchema();
  require(static_cast<bool>(validateObservationSchema(nestedDeclaration)),
          "valid nested ragged relationship DAG was rejected");
  auto nestedValues = nestedRaggedBatch();
  require(static_cast<bool>(
              validateObservationBatch(nestedDeclaration, nestedValues)),
          "valid nested ragged observation batch was rejected");

  std::vector<std::uint8_t> nestedManifest;
  WVObservationSchema decodedNested;
  require(static_cast<bool>(encodeObservationSchemaManifest(
              nestedDeclaration, nestedManifest)) &&
              static_cast<bool>(decodeObservationSchemaManifest(
                  nestedManifest, decodedNested)) &&
              static_cast<bool>(validateObservationSchema(decodedNested)) &&
              decodedNested.variables.size() == 3 &&
              decodedNested.variables[0].raggedChildAxisIdentifier ==
                  "profile" &&
              decodedNested.variables[1].raggedChildAxisIdentifier ==
                  "sample",
          "nested ragged relationship DAG did not survive manifest round trip");

  auto selfEdge = nestedRaggedSchema();
  selfEdge.variables[0].raggedChildAxisIdentifier = "pass";
  require(!validateObservationSchema(selfEdge),
          "ragged self edge was accepted");

  auto cyclic = nestedRaggedSchema();
  cyclic.variables.push_back(
      {"sample-pass-count", "sample_pass_count",
       WVObservationScalarType::integer64, {"sample"},
       WVObservationValueLayout::flat, "1", "passes in each sample", {},
       WVObservationCoordinateRole::none,
       WVObservationRaggedRole::rowCount, "pass"});
  require(!validateObservationSchema(cyclic),
          "cyclic ragged relationship graph was accepted");

  auto multipleParents = nestedRaggedSchema();
  multipleParents.axes.push_back(
      {"transect", "transect", WVObservationAxisKind::unlimited, 0,
       WVObservationCoordinateRole::none});
  multipleParents.variables.push_back(
      {"transect-profile-count", "transect_profile_count",
       WVObservationScalarType::integer64, {"transect"},
       WVObservationValueLayout::flat, "1", "profiles in each transect", {},
       WVObservationCoordinateRole::none,
       WVObservationRaggedRole::rowCount, "profile"});
  require(!validateObservationSchema(multipleParents),
          "ragged child with multiple parents was accepted");

  auto wrongRelationshipRank = nestedRaggedSchema();
  wrongRelationshipRank.variables[0].dimensionIdentifiers = {"pass",
                                                               "profile"};
  require(!validateObservationSchema(wrongRelationshipRank),
          "rank-two ragged relationship variable was accepted");

  auto wrongRelationshipType = nestedRaggedSchema();
  wrongRelationshipType.variables[0].scalarType =
      WVObservationScalarType::real64;
  require(!validateObservationSchema(wrongRelationshipType),
          "non-integer ragged relationship variable was accepted");

  auto fixedParent = nestedRaggedSchema();
  fixedParent.variables[0].dimensionIdentifiers = {"depth"};
  require(!validateObservationSchema(fixedParent),
          "fixed ragged relationship parent axis was accepted");

  auto fixedChild = nestedRaggedSchema();
  fixedChild.variables[0].raggedChildAxisIdentifier = "depth";
  require(!validateObservationSchema(fixedChild),
          "fixed ragged relationship child axis was accepted");

  auto wrongNestedCounts = nestedRaggedBatch();
  wrongNestedCounts.values[0].ownedInteger64 = {1, 1};
  require(!validateObservationBatch(nestedDeclaration, wrongNestedCounts),
          "nested ragged row counts that do not span their child were accepted");

  auto negativeNestedCount = nestedRaggedBatch();
  negativeNestedCount.values[0].ownedInteger64 = {-1, 4};
  require(!validateObservationBatch(nestedDeclaration, negativeNestedCount),
          "negative nested ragged row count was accepted");

  auto nonzeroFirstOffset = nestedRaggedBatch();
  nonzeroFirstOffset.values[1].ownedInteger64 = {1, 2, 2};
  require(!validateObservationBatch(nestedDeclaration, nonzeroFirstOffset),
          "nested ragged row offsets with a nonzero origin were accepted");

  auto decreasingOffsets = nestedRaggedBatch();
  decreasingOffsets.values[1].ownedInteger64 = {0, 3, 2};
  require(!validateObservationBatch(nestedDeclaration, decreasingOffsets),
          "decreasing nested ragged row offsets were accepted");

  auto outOfBoundsOffsets = nestedRaggedBatch();
  outOfBoundsOffsets.values[1].ownedInteger64 = {0, 2, 5};
  require(!validateObservationBatch(nestedDeclaration, outOfBoundsOffsets),
          "out-of-bounds nested ragged row offsets were accepted");

  auto emptyParentWithSamples = nestedRaggedBatch();
  emptyParentWithSamples.values[0] = WVObservationValue::ownInteger(
      "pass-profile-count", {1}, std::vector<std::int64_t>{0});
  emptyParentWithSamples.values[1] = WVObservationValue::ownInteger(
      "profile-sample-offset", {0}, {});
  emptyParentWithSamples.values[2] = WVObservationValue::ownReal(
      "sample-value", {2, 1}, std::vector<double>{1.0, 2.0});
  require(!validateObservationBatch(nestedDeclaration,
                                    emptyParentWithSamples),
          "empty ragged parent with a nonempty child axis was accepted");

  WVObservationBatch emptyNested;
  emptyNested.schemaIdentifier = nestedDeclaration.identifier;
  emptyNested.values.push_back(
      WVObservationValue::ownInteger("pass-profile-count", {0}, {}));
  emptyNested.values.push_back(
      WVObservationValue::ownInteger("profile-sample-offset", {0}, {}));
  emptyNested.values.push_back(
      WVObservationValue::ownReal("sample-value", {2, 0}, {}));
  require(static_cast<bool>(
              validateObservationBatch(nestedDeclaration, emptyNested)),
          "zero-length nested ragged observation batch was rejected");

  std::cout << "observation schema/batch tests passed\n";
  return 0;
}
