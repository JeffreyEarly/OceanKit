#include "WVRunRequest.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using wavevortex::runtime::cli::WVRunRequest;
using wavevortex::runtime::cli::WVRunRequestIntegrationMethod;
using wavevortex::runtime::cli::WVRunRequestStepPolicy;
using wavevortex::runtime::cli::WVRunRequestTimeStepConstraint;
using wavevortex::runtime::cli::decodeRunRequest;

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

void write(const std::filesystem::path &path, const std::string &value) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  require(static_cast<bool>(output), "unable to write " + path.string());
  output << value;
}

std::string fixedRequest(const std::string &extra = {}) {
  return R"({
    "schemaIdentifier":"wave-vortex-run-request-v1",
    "schemaVersion":1,
    "modelFiles":["model.nc"],
    "integration":{"method":"fixed-rk4","finalTime":12,"initialStep":0.25},
    "output":{"policy":"create","destinations":{"primary":"output.nc"}},
    "execution":{"fftProvider":"reference","threads":1},
    "report":"report.json")" + extra + "}";
}

std::string v2Request(const std::string &integration,
                      const std::string &extra = {}) {
  return std::string(R"({
    "schemaIdentifier":"wave-vortex-run-request-v2",
    "schemaVersion":2,
    "modelFiles":["model.nc"],
    "integration":)") + integration + R"(,
    "output":{"policy":"append","destinations":{}},
    "execution":{"fftProvider":"reference","threads":1},
    "report":"v2-report.json")" + extra + "}";
}

} // namespace

int main() {
  try {
    const auto root = std::filesystem::temp_directory_path() /
                      "wave-vortex-run-request-contract";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root / "subdirectory");
    write(root / "model.nc", "fixture");
    write(root / "request.json", fixedRequest());

    WVRunRequest request;
    auto status = decodeRunRequest((root / "request.json").string(), request);
    require(static_cast<bool>(status), status.message);
    require(request.requestPath ==
                std::filesystem::weakly_canonical(root / "request.json"),
            "request path was not canonicalized");
    require(request.modelFiles.size() == 1 &&
                request.modelFiles.front() ==
                    std::filesystem::weakly_canonical(root / "model.nc"),
            "model path was not resolved relative to the request");
    require(request.integration.method == WVRunRequestIntegrationMethod::fixedRK4 &&
                request.integration.initialStep == 0.25 &&
                request.integration.maximumStep == 0.25,
            "fixed integration contract was decoded incorrectly");
    require(request.outputPolicy == "create" &&
                request.destinations.size() == 1 &&
                request.destinations.front().fileIdentifier == "primary",
            "output destination contract was decoded incorrectly");
    const auto originalWorkingDirectory = std::filesystem::current_path();
    std::filesystem::current_path(root / "subdirectory");
    status = decodeRunRequest("../request.json", request);
    std::filesystem::current_path(originalWorkingDirectory);
    require(static_cast<bool>(status),
            "request failed outside its authoring directory: " +
                status.message);

    write(root / "adaptive.json", R"({
      "schemaIdentifier":"wave-vortex-run-request-v1",
      "schemaVersion":1,
      "modelFiles":["model.nc"],
      "integration":{"method":"adaptive-rk23","finalTime":12,"initialStep":0.25,"maximumStep":1,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001},
      "output":{"policy":"append","destinations":{}},
      "execution":{"fftProvider":"reference","threads":1},
      "report":"adaptive-report.json"
    })");
    status = decodeRunRequest((root / "adaptive.json").string(), request);
    require(static_cast<bool>(status) &&
                request.integration.method == WVRunRequestIntegrationMethod::adaptiveRK23 &&
                request.destinations.empty(),
            "adaptive append request was not accepted");

    const auto fixtureRoot = std::filesystem::path(WV_RUNTIME_FIXTURE_DIR) /
                             "run-requests";
    struct AcceptedFixture {
      const char *name;
      WVRunRequestIntegrationMethod method;
      WVRunRequestStepPolicy policy;
      WVRunRequestTimeStepConstraint constraint;
      int version;
    };
    const AcceptedFixture accepted[] = {
        {"v1-fixed.json",WVRunRequestIntegrationMethod::fixedRK4,
         WVRunRequestStepPolicy::explicitStep,
         WVRunRequestTimeStepConstraint::notApplicable,1},
        {"v1-adaptive-rk23.json",WVRunRequestIntegrationMethod::adaptiveRK23,
         WVRunRequestStepPolicy::adaptive,
         WVRunRequestTimeStepConstraint::notApplicable,1},
        {"v2-fixed-explicit.json",WVRunRequestIntegrationMethod::fixedRK4,
         WVRunRequestStepPolicy::explicitStep,
         WVRunRequestTimeStepConstraint::notApplicable,2},
        {"v2-fixed-cfl.json",WVRunRequestIntegrationMethod::fixedRK4,
         WVRunRequestStepPolicy::cflSelected,
         WVRunRequestTimeStepConstraint::minimum,2},
        {"v2-adaptive-rk23.json",WVRunRequestIntegrationMethod::adaptiveRK23,
         WVRunRequestStepPolicy::adaptive,
         WVRunRequestTimeStepConstraint::notApplicable,2},
        {"v2-adaptive-rk45.json",WVRunRequestIntegrationMethod::adaptiveRK45,
         WVRunRequestStepPolicy::adaptive,
         WVRunRequestTimeStepConstraint::notApplicable,2},
        {"v2-adaptive-rk78.json",WVRunRequestIntegrationMethod::adaptiveRK78,
         WVRunRequestStepPolicy::adaptive,
         WVRunRequestTimeStepConstraint::notApplicable,2}};
    for (const auto &fixture : accepted) {
      WVRunRequest decoded;
      status = decodeRunRequest((fixtureRoot/fixture.name).string(),decoded);
      require(static_cast<bool>(status),
              std::string("committed request fixture failed: ")+
                  fixture.name+": "+status.message);
      require(decoded.schemaVersion == fixture.version &&
                  decoded.integration.method == fixture.method &&
                  decoded.integration.stepPolicy == fixture.policy &&
                  decoded.integration.timeStepConstraint == fixture.constraint,
              std::string("committed request fixture resolved incorrectly: ")+
                  fixture.name);
    }
    for (const auto &constraint : {"advective","oscillatory","min"}) {
      write(root/"v2-cfl-constraint.json",
            v2Request(std::string("{\"method\":\"fixed-rk4\",\"finalTime\":12,\"cfl\":0.25,\"timeStepConstraint\":\"")+
                      constraint+"\"}"));
      status = decodeRunRequest((root/"v2-cfl-constraint.json").string(),
                                request);
      require(static_cast<bool>(status),
              std::string("v2 CFL constraint was rejected: ")+constraint);
    }
    write(root/"v2-defaults.json",
          v2Request(R"({"finalTime":12})"));
    status = decodeRunRequest((root/"v2-defaults.json").string(),request);
    require(static_cast<bool>(status),
            "v2 omitted defaults were rejected: "+status.message);
    require(request.integration.method ==
                WVRunRequestIntegrationMethod::adaptiveRK78 &&
                request.integration.stepPolicy ==
                    WVRunRequestStepPolicy::adaptive &&
                request.integration.relativeTolerance == 1e-3 &&
                request.integration.absoluteToleranceScale == 1e-6 &&
                !request.integration.hasMethod &&
                !request.integration.hasInitialStep &&
                !request.integration.hasMaximumStep &&
                !request.integration.hasRelativeTolerance &&
                !request.integration.hasAbsoluteToleranceScale,
            "v2 omitted integration defaults resolved incorrectly");
    write(root/"v2-default-execution.json",R"({
      "schemaIdentifier":"wave-vortex-run-request-v2",
      "schemaVersion":2,
      "modelFiles":["model.nc"],
      "integration":{"finalTime":12},
      "output":{"policy":"append","destinations":{}},
      "execution":{},
      "report":"default-report.json"
    })");
    status = decodeRunRequest((root/"v2-default-execution.json").string(),
                              request);
    require(static_cast<bool>(status) &&
                request.fftProvider == "native-fftw" &&
                request.threads == 0 && !request.hasFFTProvider &&
                !request.hasThreads,
            "v2 omitted execution defaults resolved incorrectly");
    write(root/"v2-partial-overrides.json",
          v2Request(R"({"finalTime":12,"relativeTolerance":1e-7,"maximumStep":3})"));
    status = decodeRunRequest((root/"v2-partial-overrides.json").string(),
                              request);
    require(static_cast<bool>(status) &&
                request.integration.method ==
                    WVRunRequestIntegrationMethod::adaptiveRK78 &&
                request.integration.relativeTolerance == 1e-7 &&
                request.integration.maximumStep == 3 &&
                request.integration.hasRelativeTolerance &&
                request.integration.hasMaximumStep &&
                !request.integration.hasInitialStep,
            "v2 partial default overrides resolved incorrectly");

    const auto expectFailure = [&](const std::string &name,
                                   const std::string &contents,
                                   const std::string &message,
                                   const std::string &expected = {}) {
      write(root / name, contents);
      WVRunRequest unchanged;
      unchanged.report = "sentinel";
      const auto failure = decodeRunRequest((root / name).string(), unchanged);
      require(!failure, message);
      if (!expected.empty())
        require(failure.message == expected,
                message+" returned a changed error: "+failure.message);
      require(unchanged.report == "sentinel",
              "failed decoding changed its output object");
    };
    expectFailure("unknown.json", fixedRequest(",\"unexpected\":true"),
                  "unknown field was accepted");
    expectFailure("malformed.json", "{", "malformed JSON was accepted");
    expectFailure(
        "wrong-version.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":2,"modelFiles":["model.nc"],"integration":{},"output":{},"execution":{},"report":"report.json"})",
        "wrong schema version was accepted",
        "Unsupported run-request schema version.");
    expectFailure(
        "missing-model.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["missing.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1},"output":{"policy":"create","destinations":{"primary":"output.nc"}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "missing model file was accepted");
    expectFailure(
        "duplicate-model.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc","./model.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1},"output":{"policy":"create","destinations":{"primary":"output.nc"}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "aliased model files were accepted");
    expectFailure(
        "source-alias.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1},"output":{"policy":"replace","destinations":{"primary":"model.nc"}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "source/output alias was accepted");
    expectFailure(
        "fixed-adaptive-fields.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1,"maximumStep":1},"output":{"policy":"create","destinations":{"primary":"output.nc"}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "fixed request accepted adaptive fields",
        "integration.maximumStep is valid only for adaptive-rk23.");
    expectFailure(
        "rk45-v1.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"adaptive-rk45","finalTime":12,"initialStep":0.25,"maximumStep":1,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001},"output":{"policy":"append","destinations":{}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "run-request v1 accepted unreleased adaptive-rk45 identity",
        "integration.method must be fixed-rk4 or adaptive-rk23.");
    expectFailure(
        "rk78-v1.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"adaptive-rk78","finalTime":12,"initialStep":0.25,"maximumStep":1,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001},"output":{"policy":"append","destinations":{}},"execution":{"fftProvider":"reference","threads":1},"report":"report.json"})",
        "run-request v1 accepted adaptive-rk78 identity",
        "integration.method must be fixed-rk4 or adaptive-rk23.");

    const std::vector<std::pair<std::string,std::string>> invalidV2 = {
        {"fixed-both-steps",R"({"method":"fixed-rk4","finalTime":12,"initialStep":1,"cfl":0.25,"timeStepConstraint":"min"})"},
        {"fixed-neither-step",R"({"method":"fixed-rk4","finalTime":12})"},
        {"explicit-constraint",R"({"method":"fixed-rk4","finalTime":12,"initialStep":1,"timeStepConstraint":"min"})"},
        {"cfl-no-constraint",R"({"method":"fixed-rk4","finalTime":12,"cfl":0.25})"},
        {"cfl-maximum-step",R"({"method":"fixed-rk4","finalTime":12,"cfl":0.25,"timeStepConstraint":"min","maximumStep":1})"},
        {"cfl-relative-tolerance",R"({"method":"fixed-rk4","finalTime":12,"cfl":0.25,"timeStepConstraint":"min","relativeTolerance":0.001})"},
        {"cfl-absolute-tolerance",R"({"method":"fixed-rk4","finalTime":12,"cfl":0.25,"timeStepConstraint":"min","absoluteToleranceScale":0.000001})"},
        {"adaptive-cfl",R"({"method":"adaptive-rk23","finalTime":12,"initialStep":1,"maximumStep":2,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001,"cfl":0.25})"},
        {"adaptive-constraint",R"({"method":"adaptive-rk45","finalTime":12,"initialStep":1,"maximumStep":2,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001,"timeStepConstraint":"min"})"},
        {"unknown-method",R"({"method":"ode45","finalTime":12,"initialStep":1,"maximumStep":2,"relativeTolerance":0.001,"absoluteToleranceScale":0.000001})"},
        {"unknown-integration-field",R"({"method":"fixed-rk4","finalTime":12,"initialStep":1,"deltaT":1})"},
        {"zero-cfl",R"({"method":"fixed-rk4","finalTime":12,"cfl":0,"timeStepConstraint":"min"})"},
        {"negative-initial",R"({"method":"fixed-rk4","finalTime":12,"initialStep":-1})"},
        {"nonfinite-final-time",R"({"method":"fixed-rk4","finalTime":1e999,"initialStep":1})"},
        {"nonfinite-cfl",R"({"method":"fixed-rk4","finalTime":12,"cfl":1e999,"timeStepConstraint":"min"})"}};
    for (const auto &[name,integration] : invalidV2)
      expectFailure(name+".json",v2Request(integration),
                    "invalid v2 integration form was accepted: "+name);
    expectFailure("v2-unknown-root.json",
                  v2Request(R"({"method":"fixed-rk4","finalTime":12,"initialStep":1})",
                            ",\"unexpected\":true"),
                  "v2 unknown root field was accepted");
    expectFailure(
        "report-alias.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1},"output":{"policy":"create","destinations":{"primary":"output.nc"}},"execution":{"fftProvider":"reference","threads":1},"report":"output.nc"})",
        "report/output alias was accepted");
    expectFailure(
        "reference-threads.json",
        R"({"schemaIdentifier":"wave-vortex-run-request-v1","schemaVersion":1,"modelFiles":["model.nc"],"integration":{"method":"fixed-rk4","finalTime":12,"initialStep":1},"output":{"policy":"create","destinations":{"primary":"output.nc"}},"execution":{"fftProvider":"reference","threads":2},"report":"report.json"})",
        "reference provider accepted multiple threads");

    std::filesystem::remove_all(root);
    std::cout << "Run-request decoder tests passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
