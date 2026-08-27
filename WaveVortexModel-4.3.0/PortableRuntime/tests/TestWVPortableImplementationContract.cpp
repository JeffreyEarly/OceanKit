#include "WaveVortexRuntime/WVPortableImplementationContract.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

WVPortableImplementationIdentity identity(std::string type,
                                          std::uint32_t version = 1) {
  return {std::move(type), version};
}

void verifyContractIdentity() {
  static_assert(WVPortableRuntimeSourceAPIMajorVersion == 1);
  static_assert(WVPortableRuntimeSourceAPIMinorVersion == 0);
  static_assert(WVPortableRuntimeSourceAPIIdentifier ==
                "wave-vortex-portable-source-api-v1");
  static_assert(WVPortablePairContractVersion == 1);
  static_assert(WVPortablePairContractIdentifier ==
                "wave-vortex-portable-pair-v1");
  static_assert(!std::is_polymorphic_v<WVPortableImplementationIdentity>);
  static_assert(!std::is_polymorphic_v<WVPortableCapability>);
  static_assert(noexcept(std::declval<const WVPortableCapability &>()
                             .isSupported()));
  static_assert(noexcept(std::declval<const WVPortableCapability &>()
                             .statusIdentifier()));
}

void verifySupported() {
  const auto capability = evaluatePortableCapability(
      identity("WVEulerianFields"), identity("WVEulerianFields"));
  require(capability.isSupported(), "matching identities were not supported");
  require(static_cast<bool>(capability),
          "supported capability bool conversion failed");
  require(capability.statusIdentifier() == "supported",
          "supported status identifier changed");
  require(capability.reason.empty(),
          "supported capability returned a failure reason");
  require(capability.available.has_value(),
          "supported capability lost the available identity");
}

void verifyUnavailable() {
  auto capability = evaluatePortableCapability(identity("WVTracer"),
                                                std::nullopt);
  require(capability.status == WVPortableCapabilityStatus::unavailable,
          "missing implementation did not report unavailable");
  require(!capability.reason.empty(),
          "missing implementation did not provide an actionable reason");
  require(!capability.available.has_value(),
          "missing implementation reported an available identity");

  capability = evaluatePortableCapability(identity("WVTracer"),
                                           identity("WVMooring"));
  require(capability.status == WVPortableCapabilityStatus::unavailable,
          "different type identifier did not report unavailable");
  require(capability.available->typeIdentifier == "WVMooring",
          "unavailable result lost the discovered implementation identity");
}

void verifyVersionMismatch() {
  const auto capability = evaluatePortableCapability(
      identity("WVTracer", 2), identity("WVTracer", 1));
  require(capability.status == WVPortableCapabilityStatus::versionMismatch,
          "different contract versions did not report versionMismatch");
  require(capability.statusIdentifier() == "versionMismatch",
          "version mismatch identifier changed");
  require(!capability.reason.empty(),
          "version mismatch did not provide an actionable reason");
}

void verifyInvalidContracts() {
  auto capability = evaluatePortableCapability(identity("", 1),
                                                identity("WVTracer", 1));
  require(capability.status == WVPortableCapabilityStatus::invalidContract,
          "empty requested type identifier was accepted");
  capability = evaluatePortableCapability(identity("WVTracer", 0),
                                           identity("WVTracer", 1));
  require(capability.status == WVPortableCapabilityStatus::invalidContract,
          "zero requested contract version was accepted");
  capability = evaluatePortableCapability(identity("WVTracer", 1),
                                           identity("", 1));
  require(capability.status == WVPortableCapabilityStatus::invalidContract,
          "empty available type identifier was accepted");
  capability = evaluatePortableCapability(identity("WVTracer", 1),
                                           identity("WVTracer", 0));
  require(capability.status == WVPortableCapabilityStatus::invalidContract,
          "zero available contract version was accepted");
  require(capability.statusIdentifier() == "invalidContract",
          "invalid contract identifier changed");
}

} // namespace

int main() {
  try {
    verifyContractIdentity();
    verifySupported();
    verifyUnavailable();
    verifyVersionMismatch();
    verifyInvalidContracts();
    std::cout << "Portable implementation contract tests passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Portable implementation contract tests failed: "
              << error.what() << '\n';
    return 1;
  }
}
