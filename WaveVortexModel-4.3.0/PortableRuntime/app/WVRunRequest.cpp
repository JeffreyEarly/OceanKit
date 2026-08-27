#include "WVRunRequest.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace wavevortex::runtime::cli {
namespace {

using json = nlohmann::json;
constexpr std::uintmax_t maximumRequestBytes = UINTMAX_C(1048576);

WVRunRequestStatus invalid(std::string message) {
  return WVRunRequestStatus::failure(std::move(message));
}

WVRunRequestStatus requireObject(
    const json &value, const std::set<std::string> &required,
    const std::set<std::string> &allowed, const std::string &location) {
  if (!value.is_object())
    return invalid(location + " must be a JSON object.");
  for (const auto &[name, unused] : value.items()) {
    (void)unused;
    if (allowed.find(name) == allowed.end())
      return invalid(location + " contains unknown field '" + name + "'.");
  }
  for (const auto &name : required)
    if (!value.contains(name))
      return invalid(location + " is missing required field '" + name + "'.");
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus stringValue(const json &object, const char *name,
                               const std::string &location,
                               std::string &output) {
  const auto &value = object.at(name);
  if (!value.is_string())
    return invalid(location + "." + name + " must be a string.");
  output = value.get<std::string>();
  if (output.empty())
    return invalid(location + "." + name + " must not be empty.");
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus positiveFiniteValue(const json &object, const char *name,
                                       const std::string &location,
                                       double &output) {
  const auto &value = object.at(name);
  if (!value.is_number())
    return invalid(location + "." + name + " must be a number.");
  output = value.get<double>();
  if (!std::isfinite(output) || output <= 0.0)
    return invalid(location + "." + name + " must be finite and positive.");
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus finiteValue(const json &object, const char *name,
                               const std::string &location, double &output) {
  const auto &value = object.at(name);
  if (!value.is_number())
    return invalid(location + "." + name + " must be a number.");
  output = value.get<double>();
  if (!std::isfinite(output))
    return invalid(location + "." + name + " must be finite.");
  return WVRunRequestStatus::ok();
}

std::filesystem::path resolvedPath(const std::filesystem::path &base,
                                   const std::string &raw,
                                   std::error_code &error) {
  auto path = std::filesystem::path(raw);
  if (path.is_relative())
    path = base / path;
  path = std::filesystem::absolute(path, error).lexically_normal();
  if (error)
    return {};
  const auto resolved = std::filesystem::weakly_canonical(path, error);
  return error ? std::filesystem::path{} : resolved;
}

WVRunRequestStatus integrationMethod(const std::string &identifier,
                                     WVRunRequestIntegrationMethod &method) {
  if (identifier == "fixed-rk4")
    method = WVRunRequestIntegrationMethod::fixedRK4;
  else if (identifier == "adaptive-rk23")
    method = WVRunRequestIntegrationMethod::adaptiveRK23;
  else if (identifier == "adaptive-rk45")
    method = WVRunRequestIntegrationMethod::adaptiveRK45;
  else if (identifier == "adaptive-rk78")
    method = WVRunRequestIntegrationMethod::adaptiveRK78;
  else
    return invalid("integration.method names an unsupported integration "
                   "method.");
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus parseV1Integration(const json &value,
                                      WVRunRequestIntegration &integration) {
  auto status = requireObject(
      value, {"method", "finalTime", "initialStep"},
      {"method", "finalTime", "initialStep", "maximumStep",
       "relativeTolerance", "absoluteToleranceScale"},
      "integration");
  if (!status)
    return status;
  std::string method;
  status = stringValue(value, "method", "integration", method);
  if (!status)
    return status;
  if (method != "fixed-rk4" && method != "adaptive-rk23")
    return invalid("integration.method must be fixed-rk4 or adaptive-rk23.");
  status = integrationMethod(method, integration.method);
  if (!status)
    return status;
  integration.hasMethod = true;
  status = finiteValue(value, "finalTime", "integration",
                       integration.finalTime);
  if (!status)
    return status;
  status = positiveFiniteValue(value, "initialStep", "integration",
                               integration.initialStep);
  if (!status)
    return status;
  integration.hasInitialStep = true;
  if (integration.method == WVRunRequestIntegrationMethod::fixedRK4) {
    for (const char *adaptive : {"maximumStep", "relativeTolerance",
                                 "absoluteToleranceScale"})
      if (value.contains(adaptive))
        return invalid(std::string("integration.") + adaptive +
                       " is valid only for adaptive-rk23.");
    integration.maximumStep = integration.initialStep;
    integration.stepPolicy = WVRunRequestStepPolicy::explicitStep;
    return WVRunRequestStatus::ok();
  }
  for (const char *required : {"maximumStep", "relativeTolerance",
                               "absoluteToleranceScale"})
    if (!value.contains(required))
      return invalid(std::string("integration is missing required field '") +
                     required + "'.");
  status = positiveFiniteValue(value, "maximumStep", "integration",
                               integration.maximumStep);
  if (!status)
    return status;
  integration.hasMaximumStep = true;
  status = positiveFiniteValue(value, "relativeTolerance", "integration",
                               integration.relativeTolerance);
  if (!status)
    return status;
  integration.hasRelativeTolerance = true;
  status = positiveFiniteValue(value, "absoluteToleranceScale", "integration",
                               integration.absoluteToleranceScale);
  if (status) {
    integration.hasAbsoluteToleranceScale = true;
    integration.stepPolicy = WVRunRequestStepPolicy::adaptive;
  }
  return status;
}

WVRunRequestStatus parseV2Integration(const json &value,
                                      WVRunRequestIntegration &integration) {
  auto status = requireObject(
      value, {"finalTime"},
      {"method", "finalTime", "initialStep", "cfl", "timeStepConstraint",
       "maximumStep", "relativeTolerance", "absoluteToleranceScale"},
      "integration");
  if (!status)
    return status;
  if (value.contains("method")) {
    std::string method;
    status = stringValue(value, "method", "integration", method);
    if (!status)
      return status;
    status = integrationMethod(method, integration.method);
    if (!status)
      return status;
    integration.hasMethod = true;
  }
  status = finiteValue(value, "finalTime", "integration",
                       integration.finalTime);
  if (!status)
    return status;

  const bool hasInitialStep = value.contains("initialStep");
  const bool hasCFL = value.contains("cfl");
  const bool hasConstraint = value.contains("timeStepConstraint");
  const bool hasMaximumStep = value.contains("maximumStep");
  const bool hasRelativeTolerance = value.contains("relativeTolerance");
  const bool hasAbsoluteTolerance = value.contains("absoluteToleranceScale");

  if (integration.method == WVRunRequestIntegrationMethod::fixedRK4) {
    if (hasMaximumStep || hasRelativeTolerance || hasAbsoluteTolerance)
      return invalid("Adaptive integration controls are not valid for "
                     "fixed-rk4.");
    if (hasInitialStep == hasCFL)
      return invalid("fixed-rk4 requires exactly one of initialStep or cfl.");
    if (hasInitialStep) {
      if (hasConstraint)
        return invalid("integration.timeStepConstraint requires cfl.");
      status = positiveFiniteValue(value, "initialStep", "integration",
                                   integration.initialStep);
      if (!status)
        return status;
      integration.hasInitialStep = true;
      integration.maximumStep = integration.initialStep;
      integration.stepPolicy = WVRunRequestStepPolicy::explicitStep;
      return WVRunRequestStatus::ok();
    }
    if (!hasConstraint)
      return invalid("CFL-selected fixed-rk4 requires "
                     "integration.timeStepConstraint.");
    status = positiveFiniteValue(value, "cfl", "integration",
                                 integration.cfl);
    if (!status)
      return status;
    std::string constraint;
    status = stringValue(value, "timeStepConstraint", "integration",
                         constraint);
    if (!status)
      return status;
    if (constraint == "advective")
      integration.timeStepConstraint =
          WVRunRequestTimeStepConstraint::advective;
    else if (constraint == "oscillatory")
      integration.timeStepConstraint =
          WVRunRequestTimeStepConstraint::oscillatory;
    else if (constraint == "min")
      integration.timeStepConstraint = WVRunRequestTimeStepConstraint::minimum;
    else
      return invalid("integration.timeStepConstraint must be advective, "
                     "oscillatory, or min.");
    integration.stepPolicy = WVRunRequestStepPolicy::cflSelected;
    return WVRunRequestStatus::ok();
  }

  if (hasCFL || hasConstraint)
    return invalid("CFL integration controls are valid only for fixed-rk4.");
  integration.relativeTolerance = 1e-3;
  integration.absoluteToleranceScale = 1e-6;
  if (hasInitialStep) {
    status = positiveFiniteValue(value, "initialStep", "integration",
                                 integration.initialStep);
    if (!status)
      return status;
    integration.hasInitialStep = true;
  }
  if (hasMaximumStep) {
    status = positiveFiniteValue(value, "maximumStep", "integration",
                                 integration.maximumStep);
    if (!status)
      return status;
    integration.hasMaximumStep = true;
  }
  if (hasRelativeTolerance) {
    status = positiveFiniteValue(value, "relativeTolerance", "integration",
                                 integration.relativeTolerance);
    if (!status)
      return status;
    integration.hasRelativeTolerance = true;
  }
  if (hasAbsoluteTolerance) {
    status = positiveFiniteValue(value, "absoluteToleranceScale",
                                 "integration",
                                 integration.absoluteToleranceScale);
    if (!status)
      return status;
    integration.hasAbsoluteToleranceScale = true;
  }
  integration.stepPolicy = WVRunRequestStepPolicy::adaptive;
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus parseOutput(const json &value,
                               const std::filesystem::path &base,
                               WVRunRequest &request) {
  auto status = requireObject(value, {"policy", "destinations"},
                              {"policy", "destinations"}, "output");
  if (!status)
    return status;
  status = stringValue(value, "policy", "output", request.outputPolicy);
  if (!status)
    return status;
  if (request.outputPolicy != "create" && request.outputPolicy != "replace" &&
      request.outputPolicy != "append")
    return invalid("output.policy must be create, replace, or append.");
  const auto &destinations = value.at("destinations");
  if (!destinations.is_object())
    return invalid("output.destinations must be an object keyed by stable "
                   "file identifier.");
  if ((request.outputPolicy == "create" ||
       request.outputPolicy == "replace") &&
      destinations.empty())
    return invalid("Create and replace require a nonempty complete "
                   "destination map.");
  std::set<std::string> identifiers;
  std::set<std::string> paths;
  for (const auto &[identifier, rawValue] : destinations.items()) {
    if (identifier.empty() || !identifiers.insert(identifier).second)
      return invalid("Output destination identifiers must be unique and "
                     "nonempty.");
    if (!rawValue.is_string() || rawValue.get<std::string>().empty())
      return invalid("Every output destination must be a nonempty path "
                     "string.");
    std::error_code error;
    const auto resolved =
        resolvedPath(base, rawValue.get<std::string>(), error);
    if (error || resolved.empty())
      return invalid("An output destination could not be resolved.");
    if (!paths.insert(resolved.string()).second)
      return invalid("Output destinations must not alias each other.");
    request.destinations.push_back({identifier, resolved.string()});
  }
  return WVRunRequestStatus::ok();
}

WVRunRequestStatus parseExecution(const json &value, bool isV1,
                                  WVRunRequest &request) {
  auto status = requireObject(value,
                              isV1 ? std::set<std::string>{"fftProvider", "threads"}
                                   : std::set<std::string>{},
                              {"fftProvider", "threads"}, "execution");
  if (!status)
    return status;
  request.fftProvider = isV1 ? "reference" : "native-fftw";
  if (value.contains("fftProvider")) {
    status = stringValue(value, "fftProvider", "execution",
                         request.fftProvider);
    if (!status)
      return status;
    request.hasFFTProvider = true;
  }
  if (request.fftProvider != "native-fftw" &&
      request.fftProvider != "reference")
    return invalid("execution.fftProvider must be native-fftw or reference.");
  if (value.contains("threads")) {
    const auto &threads = value.at("threads");
    if (!threads.is_number_unsigned() && !threads.is_number_integer())
      return invalid("execution.threads must be a positive integer.");
    const auto raw = threads.get<std::int64_t>();
    if (raw <= 0 || static_cast<std::uint64_t>(raw) >
                        std::numeric_limits<std::size_t>::max())
      return invalid("execution.threads must be a positive integer.");
    request.threads = static_cast<std::size_t>(raw);
    request.hasThreads = true;
  } else if (request.fftProvider == "reference") {
    request.threads = 1;
  }
  if (request.fftProvider == "reference" && request.threads != 1)
    return invalid("The reference FFT provider requires one thread.");
  return WVRunRequestStatus::ok();
}

} // namespace

const char *serializedIdentifier(WVRunRequestIntegrationMethod method) noexcept {
  switch (method) {
  case WVRunRequestIntegrationMethod::fixedRK4:
    return "fixed-rk4";
  case WVRunRequestIntegrationMethod::adaptiveRK23:
    return "adaptive-rk23";
  case WVRunRequestIntegrationMethod::adaptiveRK45:
    return "adaptive-rk45";
  case WVRunRequestIntegrationMethod::adaptiveRK78:
    return "adaptive-rk78";
  }
  return "unknown";
}

const char *serializedIdentifier(WVRunRequestStepPolicy policy) noexcept {
  switch (policy) {
  case WVRunRequestStepPolicy::explicitStep:
    return "explicit";
  case WVRunRequestStepPolicy::cflSelected:
    return "cfl";
  case WVRunRequestStepPolicy::adaptive:
    return "adaptive";
  }
  return "unknown";
}

const char *serializedIdentifier(
    WVRunRequestTimeStepConstraint constraint) noexcept {
  switch (constraint) {
  case WVRunRequestTimeStepConstraint::notApplicable:
    return "not-applicable";
  case WVRunRequestTimeStepConstraint::advective:
    return "advective";
  case WVRunRequestTimeStepConstraint::oscillatory:
    return "oscillatory";
  case WVRunRequestTimeStepConstraint::minimum:
    return "min";
  }
  return "unknown";
}

WVRunRequestStatus decodeRunRequest(const std::string &path,
                                    WVRunRequest &request) {
  try {
    std::error_code error;
    const auto requestPath = resolvedPath(std::filesystem::current_path(), path,
                                          error);
    if (error || requestPath.empty() ||
        !std::filesystem::is_regular_file(requestPath, error) || error)
      return invalid("The run-request JSON file does not exist or is not a "
                     "regular file.");
    const auto requestSize = std::filesystem::file_size(requestPath, error);
    if (error || requestSize > maximumRequestBytes)
      return invalid("The run-request JSON file exceeds the 1 MiB limit.");
    std::ifstream input(requestPath, std::ios::binary);
    if (!input)
      return invalid("The run-request JSON file cannot be opened.");
    const auto document = json::parse(input, nullptr, true, true);
    if (!document.is_object())
      return invalid("run request must be a JSON object.");
    if (!document.contains("schemaIdentifier"))
      return invalid("run request is missing required field "
                     "'schemaIdentifier'.");
    if (!document.contains("schemaVersion"))
      return invalid("run request is missing required field 'schemaVersion'.");
    std::string identifier;
    auto status = stringValue(document, "schemaIdentifier", "run request",
                              identifier);
    if (!status)
      return status;
    if (identifier != WVRunRequest::schemaV1Identifier &&
        identifier != WVRunRequest::schemaV2Identifier)
      return invalid("Unsupported run-request schema identifier '" +
                     identifier + "'.");
    const auto &version = document.at("schemaVersion");
    if (!version.is_number_unsigned() && !version.is_number_integer())
      return invalid("Unsupported run-request schema version.");
    const auto schemaVersion = version.get<std::int64_t>();
    const bool isV1 = identifier == WVRunRequest::schemaV1Identifier &&
                      schemaVersion == WVRunRequest::schemaV1Version;
    const bool isV2 = identifier == WVRunRequest::schemaV2Identifier &&
                      schemaVersion == WVRunRequest::schemaV2Version;
    if (!isV1 && !isV2)
      return invalid("Unsupported run-request schema version.");
    status = requireObject(
        document,
        {"schemaIdentifier", "schemaVersion", "modelFiles", "integration",
         "output", "execution", "report"},
        {"schemaIdentifier", "schemaVersion", "modelFiles", "integration",
         "output", "execution", "report"},
        "run request");
    if (!status)
      return status;

    WVRunRequest candidate;
    candidate.requestPath = requestPath.string();
    candidate.schemaIdentifier = identifier;
    candidate.schemaVersion = static_cast<int>(schemaVersion);
    const auto base = requestPath.parent_path();
    const auto &files = document.at("modelFiles");
    if (!files.is_array() || files.empty())
      return invalid("modelFiles must be a nonempty array of paths.");
    std::set<std::string> modelPaths;
    for (const auto &value : files) {
      if (!value.is_string() || value.get<std::string>().empty())
        return invalid("Every modelFiles entry must be a nonempty string.");
      error.clear();
      const auto resolved =
          resolvedPath(base, value.get<std::string>(), error);
      if (error || resolved.empty() ||
          !std::filesystem::is_regular_file(resolved, error) || error)
        return invalid("A referenced modelFiles path does not exist or is not "
                       "a regular file.");
      if (!modelPaths.insert(resolved.string()).second)
        return invalid("modelFiles entries must not alias each other.");
      candidate.modelFiles.push_back(resolved.string());
    }

    status = isV1
                 ? parseV1Integration(document.at("integration"),
                                      candidate.integration)
                 : parseV2Integration(document.at("integration"),
                                      candidate.integration);
    if (!status)
      return status;
    status = parseOutput(document.at("output"), base, candidate);
    if (!status)
      return status;
    status = parseExecution(document.at("execution"), isV1, candidate);
    if (!status)
      return status;
    std::string rawReport;
    status = stringValue(document, "report", "run request", rawReport);
    if (!status)
      return status;
    error.clear();
    const auto report = resolvedPath(base, rawReport, error);
    if (error || report.empty())
      return invalid("The report path could not be resolved.");
    candidate.report = report.string();

    const auto aliasesReservedPath = [&](const std::string &candidatePath) {
      if (candidatePath == requestPath.string() ||
          modelPaths.find(candidatePath) != modelPaths.end())
        return true;
      return std::any_of(candidate.destinations.begin(),
                         candidate.destinations.end(),
                         [&](const auto &destination) {
                           return destination.path == candidatePath;
                         });
    };
    if (aliasesReservedPath(candidate.report))
      return invalid("The report path must not alias the request, a model "
                     "file, or an output destination.");
    if (candidate.outputPolicy != "append")
      for (const auto &destination : candidate.destinations)
        if (modelPaths.find(destination.path) != modelPaths.end())
          return invalid("Create and replace destinations must not alias a "
                         "source model file.");

    request = std::move(candidate);
    return WVRunRequestStatus::ok();
  } catch (const json::exception &error) {
    return invalid(std::string("Invalid run-request JSON: ") + error.what());
  } catch (const std::exception &error) {
    return invalid(std::string("Run-request decoding failed: ") +
                   error.what());
  }
}

} // namespace wavevortex::runtime::cli
