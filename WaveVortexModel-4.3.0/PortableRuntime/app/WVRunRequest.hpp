#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime::cli {

struct WVRunRequestDestination {
  std::string fileIdentifier;
  std::string path;
};

enum class WVRunRequestIntegrationMethod : std::uint8_t {
  fixedRK4,
  adaptiveRK23,
  adaptiveRK45,
  adaptiveRK78
};

enum class WVRunRequestStepPolicy : std::uint8_t {
  explicitStep,
  cflSelected,
  adaptive
};

enum class WVRunRequestTimeStepConstraint : std::uint8_t {
  notApplicable,
  advective,
  oscillatory,
  minimum
};

const char *serializedIdentifier(WVRunRequestIntegrationMethod method) noexcept;
const char *serializedIdentifier(WVRunRequestStepPolicy policy) noexcept;
const char *serializedIdentifier(
    WVRunRequestTimeStepConstraint constraint) noexcept;

struct WVRunRequestIntegration {
  WVRunRequestIntegrationMethod method =
      WVRunRequestIntegrationMethod::adaptiveRK78;
  WVRunRequestStepPolicy stepPolicy = WVRunRequestStepPolicy::adaptive;
  WVRunRequestTimeStepConstraint timeStepConstraint =
      WVRunRequestTimeStepConstraint::notApplicable;
  double finalTime = 0.0;
  double initialStep = 0.0;
  double cfl = 0.0;
  double maximumStep = 0.0;
  double relativeTolerance = 0.0;
  double absoluteToleranceScale = 0.0;
  bool hasMethod = false;
  bool hasInitialStep = false;
  bool hasMaximumStep = false;
  bool hasRelativeTolerance = false;
  bool hasAbsoluteToleranceScale = false;
};

// A compact execution request. Scientific configuration, state, forcing,
// observers, schedules, and restart progress remain authoritative in the
// referenced NetCDF model bundle.
struct WVRunRequest {
  static constexpr const char *schemaV1Identifier =
      "wave-vortex-run-request-v1";
  static constexpr int schemaV1Version = 1;
  static constexpr const char *schemaV2Identifier =
      "wave-vortex-run-request-v2";
  static constexpr int schemaV2Version = 2;

  std::string requestPath;
  std::string schemaIdentifier;
  int schemaVersion = 0;
  std::vector<std::string> modelFiles;
  WVRunRequestIntegration integration;
  std::string outputPolicy;
  std::vector<WVRunRequestDestination> destinations;
  std::string fftProvider;
  std::size_t threads = 0;
  bool hasFFTProvider = false;
  bool hasThreads = false;
  std::string report;
};

struct WVRunRequestStatus {
  bool success = false;
  std::string message;
  explicit operator bool() const noexcept { return success; }

  static WVRunRequestStatus ok() { return {true, {}}; }
  static WVRunRequestStatus failure(std::string value) {
    return {false, std::move(value)};
  }
};

// Decode, validate, and normalize one request. Every returned path is
// absolute. Existing model files are canonicalized; prospective destinations
// canonicalize their nearest existing parent. Unknown fields are rejected.
WVRunRequestStatus decodeRunRequest(const std::string &path,
                                    WVRunRequest &request);

} // namespace wavevortex::runtime::cli
