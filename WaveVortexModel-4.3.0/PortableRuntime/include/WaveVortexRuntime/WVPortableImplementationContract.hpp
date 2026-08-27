#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace wavevortex::runtime {

inline constexpr std::uint32_t WVPortableRuntimeSourceAPIMajorVersion = 1;
inline constexpr std::uint32_t WVPortableRuntimeSourceAPIMinorVersion = 0;
inline constexpr std::string_view WVPortableRuntimeSourceAPIIdentifier =
    "wave-vortex-portable-source-api-v1";

inline constexpr std::uint32_t WVPortablePairContractVersion = 1;
inline constexpr std::string_view WVPortablePairContractIdentifier =
    "wave-vortex-portable-pair-v1";

enum class WVPortableCapabilityStatus : std::uint8_t {
  supported,
  unavailable,
  versionMismatch,
  invalidContract
};

struct WVPortableImplementationIdentity {
  std::string typeIdentifier;
  std::uint32_t contractVersion = 0;
};

// Value returned by construction-time portable implementation discovery.
// Accessors are allocation-free; the result owns its diagnostic strings so it
// can safely outlive the descriptor or registry used to create it.
struct WVPortableCapability {
  WVPortableCapabilityStatus status =
      WVPortableCapabilityStatus::invalidContract;
  WVPortableImplementationIdentity requested;
  std::optional<WVPortableImplementationIdentity> available;
  std::string reason;

  bool isSupported() const noexcept {
    return status == WVPortableCapabilityStatus::supported;
  }
  explicit operator bool() const noexcept { return isSupported(); }
  std::string_view statusIdentifier() const noexcept;
};

// Require an exact source-contract match. Call this during descriptor
// construction or preflight, never from an integration element loop.
WVPortableCapability evaluatePortableCapability(
    WVPortableImplementationIdentity requested,
    std::optional<WVPortableImplementationIdentity> available);

} // namespace wavevortex::runtime
