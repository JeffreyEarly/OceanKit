#include "WaveVortexRuntime/WVPortableTypedRecord.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace wavevortex::runtime;
using wavevortex::WVKernelStatusCode;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

WVPortableTypedRecord record() {
  WVPortableTypedRecord result;
  result.schemaIdentifier = "test-record-v1";
  result.schemaVersion = 1;
  result.values.push_back({"enabled", {}, std::vector<std::uint8_t>{1}});
  result.values.push_back({"ordinal", {}, std::vector<std::int64_t>{17}});
  result.values.push_back({"position", {2}, std::vector<double>{1.5, 2.5}});
  result.values.push_back({"label", {}, std::vector<std::string>{"sample"}});
  return result;
}

void verifyValidRecord() {
  const auto value = record();
  require(static_cast<bool>(validatePortableTypedRecord(value)),
          "valid typed record was rejected");
  require(value.value("ordinal") != nullptr &&
              value.value("ordinal")->valueType() ==
                  WVPortableValueType::integer,
          "typed record lookup or type changed");
  require(value.value("missing") == nullptr,
          "missing typed-record value resolved");
  require(value.encodedBytes() > 0 && value.persistentBytes() >= sizeof(value),
          "typed-record byte accounting is empty");
}

void verifyMalformedRecords() {
  auto value = record();
  value.schemaIdentifier.clear();
  require(!validatePortableTypedRecord(value),
          "empty schema identifier was accepted");
  value = record();
  value.values.back().name = value.values.front().name;
  require(!validatePortableTypedRecord(value),
          "duplicate value name was accepted");
  value = record();
  value.values[2].dimensions = {3};
  require(validatePortableTypedRecord(value).code ==
              WVKernelStatusCode::invalidShape,
          "dimension/value mismatch was accepted");
  value = record();
  value.values[0].storage = std::vector<std::uint8_t>{2};
  require(!validatePortableTypedRecord(value),
          "invalid Boolean encoding was accepted");
}

void verifyContextRules() {
  auto value = record();
  WVPortableTypedRecordValidation validation;
  validation.allowText = false;
  require(!validatePortableTypedRecord(value, validation),
          "text cursor value was accepted");
  value.values.pop_back();
  value.values[2].storage =
      std::vector<double>{1.5, std::numeric_limits<double>::infinity()};
  validation.requireFiniteReals = true;
  require(!validatePortableTypedRecord(value, validation),
          "nonfinite cursor value was accepted");
  value = record();
  validation = {};
  validation.maximumEncodedBytes = value.encodedBytes() - 1;
  require(validatePortableTypedRecord(value, validation).code ==
              WVKernelStatusCode::sizeOverflow,
          "encoded-size limit was not enforced");
  validation.maximumEncodedBytes = value.encodedBytes();
  require(static_cast<bool>(validatePortableTypedRecord(value, validation)),
          "exact encoded-size boundary was rejected");
}

void verifyCodec() {
  const auto value = record();
  std::vector<std::uint8_t> bytes;
  require(static_cast<bool>(encodePortableTypedRecord(value, bytes)) &&
              !bytes.empty(),
          "typed record did not encode");
  WVPortableTypedRecord decoded;
  require(static_cast<bool>(decodePortableTypedRecord(bytes, decoded)),
          "typed record did not decode");
  std::vector<std::uint8_t> roundTrip;
  require(static_cast<bool>(encodePortableTypedRecord(decoded, roundTrip)) &&
              roundTrip == bytes,
          "typed-record codec is not deterministic");
  const auto preserved = decoded;
  bytes.pop_back();
  require(!decodePortableTypedRecord(bytes, decoded),
          "truncated typed-record encoding was accepted");
  require(samePortableTypedRecordValue(preserved, decoded),
          "failed typed-record decoding changed the output value");
}

void verifyEmptyNumericCodec() {
  WVPortableTypedRecord value;
  value.schemaIdentifier = "empty-numeric-record-v1";
  value.schemaVersion = 1;
  value.values.push_back(
      {"booleans", {0}, std::vector<std::uint8_t>{}});
  value.values.push_back(
      {"integers", {0}, std::vector<std::int64_t>{}});
  value.values.push_back({"reals", {0}, std::vector<double>{}});

  std::vector<std::uint8_t> bytes;
  require(static_cast<bool>(encodePortableTypedRecord(value, bytes)) &&
              !bytes.empty(),
          "empty numeric typed record did not encode");
  WVPortableTypedRecord decoded;
  require(static_cast<bool>(decodePortableTypedRecord(bytes, decoded)) &&
              samePortableTypedRecordValue(value, decoded),
          "empty numeric typed record did not round trip");
  std::vector<std::uint8_t> roundTrip;
  require(static_cast<bool>(encodePortableTypedRecord(decoded, roundTrip)) &&
              roundTrip == bytes,
          "empty numeric typed-record codec is not deterministic");
}

} // namespace

int main() {
  try {
    verifyValidRecord();
    verifyMalformedRecords();
    verifyContextRules();
    verifyCodec();
    verifyEmptyNumericCodec();
    std::cout << "Portable typed-record tests passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Portable typed-record tests failed: " << error.what() << '\n';
    return 1;
  }
}
