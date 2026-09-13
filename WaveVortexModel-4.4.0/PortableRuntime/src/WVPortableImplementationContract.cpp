#include "WaveVortexRuntime/WVPortableImplementationContract.hpp"

#include <utility>

namespace wavevortex::runtime {
namespace {

bool valid(const WVPortableImplementationIdentity &identity) noexcept {
  return !identity.typeIdentifier.empty() && identity.contractVersion > 0;
}

WVPortableCapability result(
    WVPortableCapabilityStatus status,
    WVPortableImplementationIdentity requested,
    std::optional<WVPortableImplementationIdentity> available,
    std::string reason) {
  return {status, std::move(requested), std::move(available),
          std::move(reason)};
}

} // namespace

std::string_view WVPortableCapability::statusIdentifier() const noexcept {
  switch (status) {
  case WVPortableCapabilityStatus::supported:
    return "supported";
  case WVPortableCapabilityStatus::unavailable:
    return "unavailable";
  case WVPortableCapabilityStatus::versionMismatch:
    return "versionMismatch";
  case WVPortableCapabilityStatus::invalidContract:
    return "invalidContract";
  }
  return "invalidContract";
}

WVPortableCapability evaluatePortableCapability(
    WVPortableImplementationIdentity requested,
    std::optional<WVPortableImplementationIdentity> available) {
  if (!valid(requested))
    return result(WVPortableCapabilityStatus::invalidContract,
                  std::move(requested), std::move(available),
                  "The requested portable type identifier must be nonempty "
                  "and its contract version must be positive.");
  if (!available)
    return result(WVPortableCapabilityStatus::unavailable,
                  std::move(requested), std::nullopt,
                  "No matching C++ implementation is registered.");
  if (!valid(*available))
    return result(WVPortableCapabilityStatus::invalidContract,
                  std::move(requested), std::move(available),
                  "The registered C++ implementation has an invalid identity "
                  "or contract version.");
  if (requested.typeIdentifier != available->typeIdentifier)
    return result(WVPortableCapabilityStatus::unavailable,
                  std::move(requested), std::move(available),
                  "The registered C++ implementation has a different portable "
                  "type identifier.");
  if (requested.contractVersion != available->contractVersion)
    return result(WVPortableCapabilityStatus::versionMismatch,
                  std::move(requested), std::move(available),
                  "The MATLAB and C++ portable contract versions do not match.");
  return result(WVPortableCapabilityStatus::supported, std::move(requested),
                std::move(available), "");
}

} // namespace wavevortex::runtime
