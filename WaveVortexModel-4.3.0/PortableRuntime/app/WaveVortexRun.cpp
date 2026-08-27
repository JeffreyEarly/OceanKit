#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVCheckpointWriter.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/WVObserverOutputEvaluationService.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"
#include "WVModelInternalAccess.hpp"
#include "WVRunRequest.hpp"
#ifndef WV_RUNTIME_HAS_DENSE_OUTPUT
#define WV_RUNTIME_HAS_DENSE_OUTPUT 1
#endif
#if WV_RUNTIME_HAS_DENSE_OUTPUT
#include "WaveVortexRuntime/WVOutputOrchestration.hpp"
#endif
#include "WVReferenceFFTEngine.hpp"

#if WV_RUNTIME_HAS_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(__APPLE__)
#include <mach/mach.h>
#endif

namespace {

using namespace wavevortex;
using namespace wavevortex::runtime;
using Clock = std::chrono::steady_clock;

enum class ExitCode : int { success = 0, usage = 2, checkpoint = 3, provider = 4, integration = 5, output = 6 };

struct Options {
    std::string input;
    std::vector<std::string> modelFiles;
    std::string output;
    std::string provider;
    std::string report;
    std::string phaseFile;
    std::string outputDirectory;
    std::string outputPattern = "checkpoint-{index}.nc";
    cli::WVRunRequestIntegrationMethod integrator =
        cli::WVRunRequestIntegrationMethod::fixedRK4;
    cli::WVRunRequestStepPolicy stepPolicy =
        cli::WVRunRequestStepPolicy::explicitStep;
    cli::WVRunRequestTimeStepConstraint timeStepConstraint =
        cli::WVRunRequestTimeStepConstraint::notApplicable;
    std::string requestSchemaIdentifier;
    int requestSchemaVersion = 0;
    std::string restartMode = "model";
    std::string outputPolicy;
    double deltaT = 0.0;
    double initialStep = 0.0;
    double maximumStep = 0.0;
    double relativeTolerance = 1e-3;
    double absoluteTolerance = 1e-6;
    double cfl = 0.0;
    double finalTime = 0.0;
    std::size_t steps = 0;
    std::size_t threads = 0;
    std::size_t benchmarkDenseOutputsPerStep = 0;
    std::size_t benchmarkOutputCount = 0;
    std::size_t benchmarkWarmupSteps = 0;
    bool benchmarkModelDenseOutput = false;
    std::vector<double> outputTimes;
    std::vector<WVModelOutputDestination> outputDestinations;
    bool hasFinalTime = false;
    bool hasSteps = false;
    bool hasRelativeTolerance = false;
    bool hasAbsoluteTolerance = false;
    bool hasInitialStep = false;
    bool hasMaximumStep = false;
    bool hasOutputPattern = false;
    bool requestMode = false;
    bool hasRequestedIntegrator = false;
    bool hasRequestedProvider = false;
    bool hasRequestedThreads = false;
    bool usedDefaultInitialStep = false;
    std::string requestPath;

    bool scheduledOutput() const noexcept { return !outputTimes.empty() || !outputDirectory.empty() || hasOutputPattern; }
};

struct ScheduledOutputPlan {
    std::filesystem::path directory;
    std::vector<WVExplicitOutputTarget> targets;
};

struct Timings {
    double parse = 0.0;
    double preflight = 0.0;
    double provider = 0.0;
    double inspect = 0.0;
    double read = 0.0;
    double construct = 0.0;
    double prepare = 0.0;
    double cfl = 0.0;
    double startup = 0.0;
    double integrate = 0.0;
    double write = 0.0;
    double total = 0.0;
};

double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }

std::string jsonEscape(const std::string& value) {
    std::ostringstream output;
    for (const unsigned char character : value) {
        switch (character) {
            case '\\': output << "\\\\"; break;
            case '"': output << "\\\""; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20) output << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned>(character) << std::dec;
                else output << character;
        }
    }
    return output.str();
}

std::string quoted(const std::string& value) { return "\""+jsonEscape(value)+"\""; }

std::string quoted(cli::WVRunRequestIntegrationMethod method) {
    return quoted(cli::serializedIdentifier(method));
}

std::string jsonNumberOrInfinity(double value) {
    if (std::isinf(value)) return quoted(value > 0.0 ? "infinity" : "-infinity");
    if (!std::isfinite(value)) return "null";
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
}

std::string requestedNumber(bool present, double value) {
    return present ? jsonNumberOrInfinity(value) : "null";
}

std::string timingJSON(const Options &options, const Timings &timings) {
    std::ostringstream output;
    output << std::setprecision(17) << '{';
    if (options.requestSchemaVersion == cli::WVRunRequest::schemaV2Version)
        output << "\"parse\":" << timings.parse
               << ",\"preflight\":" << timings.preflight
               << ",\"provider\":" << timings.provider
               << ",\"startup\":" << timings.startup
               << ",\"cfl\":" << timings.cfl << ',';
    output << "\"inspect\":" << timings.inspect
           << ",\"read\":" << timings.read
           << ",\"construct\":" << timings.construct
           << ",\"prepare\":" << timings.prepare
           << ",\"integrate\":" << timings.integrate
           << ",\"write\":" << timings.write
           << ",\"total\":" << timings.total << '}';
    return output.str();
}

std::string integrationRequestV2JSON(
    const Options &options, const WVFixedTimeStepCandidates &candidates,
    bool evaluatedCandidates, const char *controller,
    double selectedFixedStep, double effectiveIntegrationStep,
    double effectiveMaximumStep) {
    std::ostringstream output;
    const bool isAdaptive =
        options.integrator != cli::WVRunRequestIntegrationMethod::fixedRK4;
    output << std::setprecision(17)
           << "{\"requestedMethod\":"
           << (options.hasRequestedIntegrator ? quoted(options.integrator)
                                              : "null")
           << ",\"activeMethod\":" << quoted(options.integrator)
           << ",\"controller\":" << quoted(controller)
           << ",\"requestedInitialStep\":"
           << requestedNumber(options.hasInitialStep,options.initialStep)
           << ",\"requestedMaximumStep\":"
           << requestedNumber(options.hasMaximumStep,options.maximumStep)
           << ",\"requestedRelativeTolerance\":"
           << requestedNumber(options.hasRelativeTolerance,
                              options.relativeTolerance)
           << ",\"requestedAbsoluteToleranceScale\":"
           << requestedNumber(options.hasAbsoluteTolerance,
                              options.absoluteTolerance)
           << ",\"activeRelativeTolerance\":"
           << (isAdaptive ? jsonNumberOrInfinity(options.relativeTolerance)
                          : "null")
           << ",\"activeAbsoluteToleranceScale\":"
           << (isAdaptive ? jsonNumberOrInfinity(options.absoluteTolerance)
                          : "null")
           << ",\"effectiveMaximumStep\":"
           << (isAdaptive ? jsonNumberOrInfinity(effectiveMaximumStep)
                          : "null")
           << ",\"stepPolicy\":"
           << quoted(cli::serializedIdentifier(options.stepPolicy))
           << ",\"timeStepConstraint\":"
           << quoted(cli::serializedIdentifier(options.timeStepConstraint))
           << ",\"cfl\":"
           << (options.stepPolicy == cli::WVRunRequestStepPolicy::cflSelected
                       || options.usedDefaultInitialStep
                   ? jsonNumberOrInfinity(options.usedDefaultInitialStep
                                              ? 0.5 : options.cfl)
                   : "null")
           << ",\"initialStepPolicy\":"
           << quoted(options.hasInitialStep
                         ? "explicit"
                         : options.stepPolicy ==
                                   cli::WVRunRequestStepPolicy::cflSelected
                               ? "cfl-requested"
                               : "cfl-0.5-default")
           << ",\"candidates\":{\"horizontalAdvective\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(candidates.horizontalAdvective)
                   : "null")
           << ",\"verticalAdvective\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(candidates.verticalAdvective)
                   : "null")
           << ",\"advective\":"
           << (evaluatedCandidates ? jsonNumberOrInfinity(candidates.advective)
                                   : "null")
           << ",\"oscillatory\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(candidates.oscillatory)
                   : "null")
           << ",\"effectiveHorizontalGridResolution\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(
                         candidates.effectiveHorizontalGridResolution)
                   : "null")
           << ",\"maximumHorizontalSpeed\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(candidates.maximumHorizontalSpeed)
                   : "null")
           << ",\"highestActiveWaveFrequency\":"
           << (evaluatedCandidates
                   ? jsonNumberOrInfinity(
                         candidates.highestActiveWaveFrequency)
                   : "null")
           << "},\"selectedStep\":"
           << (evaluatedCandidates ||
                       options.stepPolicy ==
                           cli::WVRunRequestStepPolicy::explicitStep
                   ? jsonNumberOrInfinity(selectedFixedStep)
                   : "null")
           << ",\"effectiveIntegrationStep\":"
           << jsonNumberOrInfinity(effectiveIntegrationStep)
           << ",\"transientWorkspaceMaximumLiveBytes\":"
           << (evaluatedCandidates
                   ? candidates.transientWorkspaceMaximumLiveBytes
                   : 0)
           << ",\"candidateEvaluationSeconds\":"
           << (evaluatedCandidates ? candidates.evaluationSeconds : 0.0)
           << ",\"noFallback\":true}";
    return output.str();
}

bool parseDouble(const std::string& text, double& value) {
    char* end = nullptr;
    value = std::strtod(text.c_str(),&end);
    return end != text.c_str() && *end == '\0' && std::isfinite(value);
}

bool parseSize(const std::string& text, std::size_t& value) {
    if (text.empty() || text.front() == '-') return false;
    char* end = nullptr;
    const auto raw = std::strtoull(text.c_str(),&end,10);
    if (end == text.c_str() || *end != '\0' || raw > std::numeric_limits<std::size_t>::max()) return false;
    value = static_cast<std::size_t>(raw);
    return true;
}

std::size_t automaticThreadCount() noexcept {
    return std::min<std::size_t>(
        18,std::max(1U,std::thread::hardware_concurrency()));
}

bool parseIntegrationMethod(
    const std::string &text,
    cli::WVRunRequestIntegrationMethod &method) noexcept {
    if (text == "fixed-rk4")
        method = cli::WVRunRequestIntegrationMethod::fixedRK4;
    else if (text == "adaptive-rk23")
        method = cli::WVRunRequestIntegrationMethod::adaptiveRK23;
    else if (text == "adaptive-rk45")
        method = cli::WVRunRequestIntegrationMethod::adaptiveRK45;
    else if (text == "adaptive-rk78")
        method = cli::WVRunRequestIntegrationMethod::adaptiveRK78;
    else
        return false;
    return true;
}

bool parseLegacyOptions(int argc, char** argv, Options& options, std::string& error) {
    if (argc < 2) { error = "INPUT is required."; return false; }
    options.input = argv[1];
    int firstOption = 2;
    if (firstOption < argc && std::string(argv[firstOption]).rfind("--",0) != 0) options.output = argv[firstOption++];
    for (int index = firstOption; index < argc; ++index) {
        const std::string name = argv[index];
        if (index+1 >= argc) { error = "Missing value for "+name+"."; return false; }
        const std::string value = argv[++index];
        if (name == "--delta-t") {
            if (!parseDouble(value,options.deltaT) || options.deltaT <= 0.0) { error = "--delta-t must be finite and positive."; return false; }
        } else if (name == "--initial-step") {
            if (!parseDouble(value,options.initialStep) || options.initialStep <= 0.0) { error = "--initial-step must be finite and positive."; return false; }
            options.hasInitialStep = true;
        } else if (name == "--maximum-step") {
            if (!parseDouble(value,options.maximumStep) || options.maximumStep <= 0.0) { error = "--maximum-step must be finite and positive."; return false; }
            options.hasMaximumStep = true;
        } else if (name == "--steps") {
            if (!parseSize(value,options.steps) || options.steps == 0) { error = "--steps must be a positive integer."; return false; }
            options.hasSteps = true;
        } else if (name == "--final-time") {
            if (!parseDouble(value,options.finalTime)) { error = "--final-time must be finite."; return false; }
            options.hasFinalTime = true;
        } else if (name == "--fft-provider") {
            options.provider = value;
        } else if (name == "--integrator") {
            if (!parseIntegrationMethod(value,options.integrator)) { error = "--integrator must be fixed-rk4, adaptive-rk23, adaptive-rk45, or adaptive-rk78."; return false; }
        } else if (name == "--restart-mode") {
            options.restartMode = value;
        } else if (name == "--output-policy") {
            options.outputPolicy = value;
        } else if (name == "--relative-tolerance") {
            if (!parseDouble(value,options.relativeTolerance) || options.relativeTolerance <= 0.0) { error = "--relative-tolerance must be finite and positive."; return false; }
            options.hasRelativeTolerance = true;
        } else if (name == "--absolute-tolerance") {
            if (!parseDouble(value,options.absoluteTolerance) || options.absoluteTolerance <= 0.0) { error = "--absolute-tolerance must be finite and positive."; return false; }
            options.hasAbsoluteTolerance = true;
        } else if (name == "--threads") {
            if (!parseSize(value,options.threads) || options.threads == 0) { error = "--threads must be a positive integer."; return false; }
        } else if (name == "--report") {
            options.report = value;
        } else if (name == "--phase-file") {
            options.phaseFile = value;
        } else if (name == "--output-time") {
            double outputTime = 0.0;
            if (!parseDouble(value,outputTime)) { error = "--output-time must be finite."; return false; }
            options.outputTimes.push_back(outputTime);
        } else if (name == "--output-directory") {
            options.outputDirectory = value;
        } else if (name == "--output-pattern") {
            options.outputPattern = value;
            options.hasOutputPattern = true;
        } else if (name == "--benchmark-dense-outputs-per-step") {
            if (!parseSize(value,options.benchmarkDenseOutputsPerStep) || (options.benchmarkDenseOutputsPerStep != 1 && options.benchmarkDenseOutputsPerStep != 4)) { error = "--benchmark-dense-outputs-per-step must be 1 or 4."; return false; }
        } else if (name == "--benchmark-output-count") {
            if (!parseSize(value,options.benchmarkOutputCount) || options.benchmarkOutputCount == 0) { error = "--benchmark-output-count must be a positive integer."; return false; }
        } else if (name == "--benchmark-warmup-steps") {
            if (!parseSize(value,options.benchmarkWarmupSteps)) { error = "--benchmark-warmup-steps must be a nonnegative integer."; return false; }
        } else if (name == "--benchmark-model-dense-output") {
            std::size_t enabled = 0;
            if (!parseSize(value,enabled) || enabled > 1) { error = "--benchmark-model-dense-output must be 0 or 1."; return false; }
            options.benchmarkModelDenseOutput = enabled == 1;
        } else {
            error = "Unknown option "+name+".";
            return false;
        }
    }
    if (!(options.deltaT > 0.0)) { error = "--delta-t is required."; return false; }
    if (options.hasSteps == options.hasFinalTime) { error = "Exactly one of --steps or --final-time is required."; return false; }
    if (options.restartMode != "model" && options.restartMode != "coefficients") { error = "--restart-mode must be model or coefficients."; return false; }
    if (options.outputPolicy != "create" && options.outputPolicy != "replace" && options.outputPolicy != "append") { error = "--output-policy must be create, replace, or append."; return false; }
    if (options.restartMode == "model" && options.outputPolicy != "append") { error = "Model-graph continuation requires --output-policy append."; return false; }
    if (options.restartMode == "coefficients" && options.outputPolicy == "append") { error = "Coefficient-only continuation supports create or replace, not append."; return false; }
    if (options.restartMode == "model" && !options.output.empty()) { error = "Model-graph continuation appends its restored output and does not accept positional OUTPUT."; return false; }
    if (options.restartMode == "model" && options.scheduledOutput()) { error = "Model-graph continuation uses its restored output schedules."; return false; }
    if (options.restartMode == "model" && options.hasSteps) { error = "Model-graph continuation requires --final-time so restored schedules have a bounded continuation interval."; return false; }
    if (options.scheduledOutput()) {
        if (!options.output.empty()) { error = "Scheduled output cannot be combined with positional OUTPUT."; return false; }
        if (options.outputTimes.empty()) { error = "Scheduled output requires at least one --output-time."; return false; }
        if (options.outputDirectory.empty()) { error = "Scheduled output requires --output-directory."; return false; }
        if (!options.hasFinalTime || options.hasSteps) { error = "Scheduled output requires --final-time and cannot be combined with --steps."; return false; }
    } else if (options.restartMode == "coefficients" && options.output.empty()) {
        error = "OUTPUT is required unless scheduled output is configured.";
        return false;
    }
    if (options.provider != "native-fftw" && options.provider != "reference") { error = "--fft-provider must be native-fftw or reference."; return false; }
    if (options.integrator == cli::WVRunRequestIntegrationMethod::fixedRK4 && (options.hasRelativeTolerance || options.hasAbsoluteTolerance || options.hasInitialStep || options.hasMaximumStep)) { error = "Adaptive tolerance and step-control options require an adaptive integrator."; return false; }
    if (options.provider == "reference" && options.threads > 1) { error = "The reference provider supports only one thread."; return false; }
    if (options.threads == 0)
        options.threads = options.provider == "reference"
            ? 1 : automaticThreadCount();
    if ((options.benchmarkDenseOutputsPerStep != 0 || options.benchmarkWarmupSteps != 0) && !options.hasSteps) { error = "Author-only benchmark controls require --steps."; return false; }
    if (options.benchmarkOutputCount != 0 && (!options.hasFinalTime || options.integrator == cli::WVRunRequestIntegrationMethod::fixedRK4)) { error = "--benchmark-output-count requires an adaptive integrator with --final-time."; return false; }
    if (options.scheduledOutput() && (options.benchmarkDenseOutputsPerStep != 0 || options.benchmarkOutputCount != 0 || options.benchmarkWarmupSteps != 0)) { error = "Scheduled output cannot be combined with author-only output benchmark controls."; return false; }
    if (options.benchmarkModelDenseOutput && options.restartMode != "model") { error = "--benchmark-model-dense-output requires model-graph continuation."; return false; }
    options.modelFiles = {options.input};
    return true;
}

bool parseOptions(int argc, char** argv, Options& options, std::string& error) {
    bool namesRequest = false;
    for (int index = 1; index < argc; ++index)
        namesRequest |= std::string(argv[index]) == "--request";
    if (!namesRequest) return parseLegacyOptions(argc,argv,options,error);
    if (argc != 3 || std::string(argv[1]) != "--request") {
        error = "--request must be the only command-line option and requires one JSON path.";
        return false;
    }
    cli::WVRunRequest request;
    const auto status = cli::decodeRunRequest(argv[2],request);
    if (!status) {
        error = status.message;
        return false;
    }
    options.requestMode = true;
    options.requestPath = request.requestPath;
    options.modelFiles = std::move(request.modelFiles);
    options.input = options.modelFiles.front();
    options.integrator = request.integration.method;
    options.stepPolicy = request.integration.stepPolicy;
    options.timeStepConstraint = request.integration.timeStepConstraint;
    options.requestSchemaIdentifier = request.schemaIdentifier;
    options.requestSchemaVersion = request.schemaVersion;
    options.restartMode = "model";
    options.outputPolicy = request.outputPolicy;
    options.deltaT = request.integration.initialStep;
    options.initialStep = request.integration.initialStep;
    options.cfl = request.integration.cfl;
    options.maximumStep = request.integration.maximumStep;
    options.relativeTolerance = request.integration.relativeTolerance;
    options.absoluteTolerance = request.integration.absoluteToleranceScale;
    options.finalTime = request.integration.finalTime;
    options.threads = request.threads;
    options.provider = request.fftProvider;
    options.report = request.report;
    options.outputDestinations.reserve(request.destinations.size());
    for (auto &destination : request.destinations)
        options.outputDestinations.push_back(
            {std::move(destination.fileIdentifier),std::move(destination.path)});
    options.hasFinalTime = true;
    options.hasRequestedIntegrator = request.integration.hasMethod;
    options.hasRequestedProvider = request.hasFFTProvider;
    options.hasRequestedThreads = request.hasThreads;
    if (options.requestSchemaVersion == cli::WVRunRequest::schemaV2Version) {
        options.hasInitialStep = request.integration.hasInitialStep;
        options.hasMaximumStep = request.integration.hasMaximumStep;
        options.hasRelativeTolerance = request.integration.hasRelativeTolerance;
        options.hasAbsoluteTolerance =
            request.integration.hasAbsoluteToleranceScale;
    } else if (options.stepPolicy == cli::WVRunRequestStepPolicy::adaptive) {
        options.hasInitialStep = true;
        options.hasMaximumStep = true;
        options.hasRelativeTolerance = true;
        options.hasAbsoluteTolerance = true;
    } else if (options.requestSchemaVersion ==
                   cli::WVRunRequest::schemaV2Version &&
               options.stepPolicy ==
                   cli::WVRunRequestStepPolicy::explicitStep) {
        options.hasInitialStep = true;
    }
    if (options.threads == 0)
        options.threads = options.provider == "reference"
            ? 1 : automaticThreadCount();
    return true;
}

#if WV_RUNTIME_HAS_DENSE_OUTPUT
class BenchmarkSink final : public WVOutputSink {
public:
    WVKernelStatus preflight(const WVOutputPlan&) override {
        return WVKernelStatus::ok();
    }
    WVKernelStatus deliver(const WVOutputEvent& event,
                           const WVOutputRouteView&,
                           WVOutputDeliveryResult&) override {
        if (event.kind == WVOutputEventKind::interpolated) ++interpolatedCount;
        return WVKernelStatus::ok();
    }
    std::size_t interpolatedCount = 0;
};

std::vector<double> interiorOutputTimes(double initialTime, double stepSize, std::size_t stepCount, std::size_t outputsPerStep) {
    std::vector<double> times;
    times.reserve(stepCount*outputsPerStep);
    for (std::size_t step = 0; step < stepCount; ++step) {
        for (std::size_t output = 1; output <= outputsPerStep; ++output) times.push_back(initialTime+(static_cast<double>(step)+static_cast<double>(output)/static_cast<double>(outputsPerStep+1))*stepSize);
    }
    return times;
}

std::vector<double> uniformInteriorOutputTimes(double initialTime, double finalTime, std::size_t outputCount) {
    std::vector<double> times;
    times.reserve(outputCount);
    for (std::size_t output = 1; output <= outputCount; ++output) times.push_back(initialTime+(finalTime-initialTime)*static_cast<double>(output)/static_cast<double>(outputCount+1));
    return times;
}
#endif

std::string replaceAll(std::string value, const std::string& token, const std::string& replacement) {
    std::size_t position = 0;
    while ((position = value.find(token,position)) != std::string::npos) {
        value.replace(position,token.size(),replacement);
        position += replacement.size();
    }
    return value;
}

std::string outputTimeText(double value) {
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return output.str();
}

std::filesystem::path normalizedPath(const std::filesystem::path& value) {
    std::error_code error;
    const auto absolute = std::filesystem::absolute(value,error).lexically_normal();
    if (error) return value.lexically_normal();
    const auto canonical = std::filesystem::weakly_canonical(absolute,error);
    return error ? absolute : canonical;
}

bool prepareScheduledOutput(
    const Options& options,
    double initialTime,
    ScheduledOutputPlan& plan,
    ExitCode& failureCode,
    std::string& error) {
    for (std::size_t index = 0; index < options.outputTimes.size(); ++index) {
        if (!std::isfinite(options.outputTimes[index]) ||
            options.outputTimes[index] < initialTime ||
            options.outputTimes[index] > options.finalTime ||
            (index > 0 && !(options.outputTimes[index] > options.outputTimes[index-1]))) {
            failureCode = ExitCode::usage;
            error = "Output times must be finite, strictly increasing, and inside the integration interval.";
            return false;
        }
    }
    const std::filesystem::path patternPath(options.outputPattern);
    const bool hasIndex = options.outputPattern.find("{index}") != std::string::npos;
    const bool hasTime = options.outputPattern.find("{time}") != std::string::npos;
    std::string unknownTokens = replaceAll(replaceAll(options.outputPattern,"{index}",""),"{time}","");
    if (options.outputPattern.empty() || patternPath.has_parent_path() || patternPath.extension() != ".nc" || (!hasIndex && !hasTime) || unknownTokens.find_first_of("{}") != std::string::npos) {
        failureCode = ExitCode::usage;
        error = "--output-pattern must be a .nc filename containing {index} or {time}, without directories or unknown tokens.";
        return false;
    }
    std::error_code filesystemError;
    plan.directory = normalizedPath(options.outputDirectory);
    std::filesystem::create_directories(plan.directory,filesystemError);
    if (filesystemError || !std::filesystem::is_directory(plan.directory,filesystemError)) {
        failureCode = ExitCode::output;
        error = "Unable to create or use scheduled output directory: "+filesystemError.message();
        return false;
    }
    std::set<std::filesystem::path> reserved;
    reserved.insert(normalizedPath(options.input));
    if (!options.report.empty()) reserved.insert(normalizedPath(options.report));
    if (!options.phaseFile.empty()) reserved.insert(normalizedPath(options.phaseFile));
    std::set<std::filesystem::path> destinations;
    plan.targets.clear();
    plan.targets.reserve(options.outputTimes.size());
    for (std::size_t index = 0; index < options.outputTimes.size(); ++index) {
        std::ostringstream ordinal;
        ordinal << std::setw(6) << std::setfill('0') << index+1;
        const std::string filename = replaceAll(replaceAll(options.outputPattern,"{index}",ordinal.str()),"{time}",outputTimeText(options.outputTimes[index]));
        const auto destination = normalizedPath(plan.directory/filename);
        if (!destinations.insert(destination).second) {
            failureCode = ExitCode::usage;
            error = "Scheduled output pattern expands to duplicate destinations.";
            return false;
        }
        if (reserved.find(destination) != reserved.end()) {
            failureCode = ExitCode::usage;
            error = "Scheduled output destination aliases an input, report, or phase file.";
            return false;
        }
        filesystemError.clear();
        const auto fileStatus = std::filesystem::symlink_status(destination,filesystemError);
        if (filesystemError == std::errc::no_such_file_or_directory) filesystemError.clear();
        if (filesystemError || fileStatus.type() != std::filesystem::file_type::not_found) {
            failureCode = ExitCode::output;
            error = filesystemError ? "Unable to inspect scheduled output destination: "+filesystemError.message() : "Scheduled output destination already exists: "+destination.string();
            return false;
        }
        plan.targets.push_back({options.outputTimes[index],destination.string()});
    }
    return true;
}

void phase(const Options& options, const std::string& value) {
    if (options.phaseFile.empty()) return;
    const auto temporary = options.phaseFile+".tmp";
    {
        std::ofstream output(temporary,std::ios::binary|std::ios::trunc);
        if (!output) return;
        output << value << '\n';
    }
    std::error_code ignored;
    std::filesystem::rename(temporary,options.phaseFile,ignored);
}

void phasePlateau(const Options& options, const std::string& value) {
    phase(options,value);
    if (!options.phaseFile.empty()) std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

bool waitForPhaseSample(const std::string& expected) {
    const auto* path = std::getenv("WV_RSS_PHASE_ACK");
    if (path == nullptr || *path == '\0') return true;
    const auto deadline = Clock::now()+std::chrono::seconds(5);
    while (Clock::now() < deadline) {
        std::ifstream input(path);
        std::string sampled;
        if (input && std::getline(input,sampled) && sampled == expected)
            return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return false;
}

class PhaseReportingCheckpointSink final : public WVOutputSink {
public:
    PhaseReportingCheckpointSink(
        std::shared_ptr<const WVExtensionCatalog> catalog,
        WVCheckpoint checkpointTemplate, const Options& options)
        : sink_(std::move(catalog), std::move(checkpointTemplate)),
          options_(options) {}

    WVKernelStatus preflight(const WVOutputPlan& plan) override {
        return sink_.preflight(plan);
    }
    WVKernelStatus deliver(const WVOutputEvent& event,
                           const WVOutputRouteView& route,
                           WVOutputDeliveryResult& result) override {
        const auto previousCount = sink_.metrics().checkpointWriteCount;
        const auto status = sink_.deliver(event,route,result);
        if (status && sink_.metrics().checkpointWriteCount != previousCount)
            phase(options_,"output-committed:"+std::to_string(sink_.metrics().checkpointWriteCount));
        return status;
    }

    WVCheckpointOutputSink& sink() noexcept { return sink_; }
    const WVCheckpointOutputSink& sink() const noexcept { return sink_; }

private:
    WVCheckpointOutputSink sink_;
    const Options& options_;
};

std::size_t currentRSSBytes() {
#if defined(__APPLE__)
    mach_task_basic_info information{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&information),&count) == KERN_SUCCESS) return static_cast<std::size_t>(information.resident_size);
#elif defined(__linux__)
    std::ifstream input("/proc/self/statm");
    std::size_t total = 0;
    std::size_t resident = 0;
    if (input >> total >> resident) return resident*static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
#endif
    return 0;
}

std::size_t peakRSSBytes() {
#if defined(__APPLE__) || defined(__linux__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF,&usage) == 0) {
#if defined(__APPLE__)
        return static_cast<std::size_t>(usage.ru_maxrss);
#else
        return static_cast<std::size_t>(usage.ru_maxrss)*1024U;
#endif
    }
#endif
    return 0;
}

std::string forcingJSON(const WVFrozenForcingSchedule& schedule) {
    std::ostringstream output;
    output << '[';
    for (std::size_t index = 0; index < schedule.entries.size(); ++index) {
        if (index != 0) output << ',';
        const auto& entry = schedule.entries[index];
        output << "{\"ordinal\":" << entry.ordinal << ",\"type\":" << quoted(entry.typeIdentifier) << ",\"name\":" << quoted(entry.name) << ",\"stage\":" << static_cast<unsigned>(entry.stage) << ",\"priority\":" << static_cast<unsigned>(entry.priority) << '}';
    }
    output << ']';
    return output.str();
}

std::string adaptiveStepDiagnosticsJSON(
    const std::vector<WVAdaptiveRK23StepDiagnostic>& diagnostics) {
    std::ostringstream output;
    output << std::setprecision(17) << '[';
    for (std::size_t index = 0; index < diagnostics.size(); ++index) {
        if (index != 0) output << ',';
        const auto& value = diagnostics[index];
        output << "{\"initialTime\":" << value.initialTime
               << ",\"acceptedStepSize\":" << value.acceptedStepSize
               << ",\"normalizedError\":" << value.normalizedError
               << ",\"nextStepSize\":" << value.nextStepSize
               << ",\"rejectedAttemptCount\":" << value.rejectedAttemptCount
               << ",\"rightHandSideEvaluationCount\":" << value.rightHandSideEvaluationCount
               << ",\"reusedFSALDerivative\":" << (value.reusedFSALDerivative ? "true" : "false") << '}';
    }
    output << ']';
    return output.str();
}

std::string adaptiveStageBufferLastUseJSON(
    const WVAdaptiveRKStageBufferLastUse* records, std::size_t count) {
    std::ostringstream output;
    output << '[';
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0) output << ',';
        output << "{\"buffer\":" << quoted(records[index].bufferIdentifier)
               << ",\"producer\":" << quoted(records[index].producer)
               << ",\"stage\":" << records[index].stage
               << ",\"lastUse\":" << quoted(records[index].lastUse) << '}';
    }
    output << ']';
    return output.str();
}

std::string unsignedIntegerArrayJSON(const std::vector<std::uint64_t>& values) {
    std::ostringstream output;
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) output << ',';
        output << quoted(std::to_string(values[index]));
    }
    output << ']';
    return output.str();
}

std::string stringArrayJSON(const std::vector<std::string>& values) {
    std::ostringstream output;
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) output << ',';
        output << quoted(values[index]);
    }
    output << ']';
    return output.str();
}

std::string sizeArrayJSON(const std::vector<std::size_t>& values) {
    std::ostringstream output;
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) output << ',';
        output << values[index];
    }
    output << ']';
    return output.str();
}

std::string destinationMapJSON(
    const std::vector<WVModelOutputDestination>& destinations) {
    std::ostringstream output;
    output << '{';
    for (std::size_t index = 0; index < destinations.size(); ++index) {
        if (index != 0) output << ',';
        output << quoted(destinations[index].fileIdentifier) << ':'
               << quoted(destinations[index].path);
    }
    output << '}';
    return output.str();
}

const char* eventKindName(WVOutputEventKind kind) noexcept {
    switch (kind) {
        case WVOutputEventKind::initial: return "initial";
        case WVOutputEventKind::interpolated: return "interpolated";
        case WVOutputEventKind::acceptedEndpoint: return "accepted-endpoint";
    }
    return "unknown";
}

std::string scheduledOutputJSON(const Options& options, const ScheduledOutputPlan& plan, const WVCheckpointOutputSink* sink) {
    const auto* records = sink == nullptr ? nullptr : &sink->records();
    const auto committedCount = sink == nullptr ? 0 : sink->metrics().checkpointWriteCount;
    const auto writeSeconds = sink == nullptr ? 0.0 : sink->metrics().checkpointWriteSeconds;
    std::ostringstream output;
    output << std::setprecision(17) << "{\"mode\":\"scheduled\",\"directory\":" << quoted(plan.directory.string()) << ",\"pattern\":" << quoted(options.outputPattern) << ",\"requestedCount\":" << plan.targets.size() << ",\"committedCount\":" << committedCount << ",\"writeSeconds\":" << writeSeconds << ",\"integrateIncludesScheduledWrites\":true,\"records\":[";
    for (std::size_t index = 0; index < plan.targets.size(); ++index) {
        if (index != 0) output << ',';
        const auto& target = plan.targets[index];
        output << "{\"ordinal\":" << index+1 << ",\"requestedTime\":" << target.requestedTime << ",\"path\":" << quoted(target.destination);
        if (records != nullptr && index < records->size()) {
            const auto& record = (*records)[index];
            output << ",\"emittedTime\":" << record.emittedTime << ",\"eventKind\":" << quoted(eventKindName(record.eventKind)) << ",\"writeSeconds\":" << record.writeSeconds << ",\"status\":" << quoted(record.committed ? "committed" : "failed") << ",\"failure\":" << quoted(record.failure);
        } else {
            output << ",\"emittedTime\":null,\"eventKind\":null,\"writeSeconds\":0,\"status\":\"pending\",\"failure\":\"\"";
        }
        output << '}';
    }
    output << "]}";
    return output.str();
}

std::string failureJSON(ExitCode code, const std::string& stage, const std::string& message, const std::string& location = {}, const std::string& scheduledOutput = {}) {
    std::ostringstream output;
    output << "{\"schemaVersion\":\"wave-vortex-run-v1\",\"status\":\"failed\",\"exitCode\":" << static_cast<int>(code) << ",\"failure\":{\"stage\":" << quoted(stage) << ",\"message\":" << quoted(message) << ",\"location\":" << quoted(location) << '}';
    if (!scheduledOutput.empty()) output << ",\"scheduledOutput\":" << scheduledOutput;
    output << '}';
    return output.str();
}

void emit(const std::string& report, const std::string& path, std::ostream& stream) {
    stream << report << '\n';
    if (!path.empty()) {
        std::ofstream output(path,std::ios::binary|std::ios::trunc);
        if (output) output << report << '\n';
    }
}

std::unique_ptr<WVFFTEngine> provider(const Options& options, std::string& version, std::string& baseLibrary, std::string& threadLibrary, std::string& error) {
    if (options.provider == "reference") return std::make_unique<WVReferenceFFTEngine>();
#if WV_RUNTIME_HAS_NATIVE_FFTW
    std::unique_ptr<WVFFTEngine> result;
    const auto status = WVFFTWEngine::create(options.threads,result);
    if (!status) {
        error = status.message +
            " Rebuild the source runner with "
            "PortableRuntime/buildWaveVortexRun.sh on a supported "
            "Apple-silicon host.";
        return {};
    }
    const auto identity = WVFFTWEngine::linkedLibraries();
    version = identity.version;
    baseLibrary = identity.baseLibrary;
    threadLibrary = identity.threadLibrary;
    const std::filesystem::path expected = std::filesystem::weakly_canonical(WV_RUNTIME_EXPECTED_FFTW_ROOT);
    const auto base = std::filesystem::weakly_canonical(baseLibrary);
    const auto threads = std::filesystem::weakly_canonical(threadLibrary);
    const auto underExpected = [&expected](const std::filesystem::path& value) { return value.string().rfind(expected.string()+std::string(1,std::filesystem::path::preferred_separator),0) == 0; };
    if (version.find("3.3.11") == std::string::npos || !underExpected(base) || !underExpected(threads)) {
        error = "Loaded FFTW identity does not match the configured pinned "
                "3.3.11 provider. Rebuild the source runner with "
                "PortableRuntime/buildWaveVortexRun.sh.";
        return {};
    }
    return result;
#else
    (void)version; (void)baseLibrary; (void)threadLibrary;
    error = "This runner was built without the native FFTW provider. Build "
            "the source runner with PortableRuntime/buildWaveVortexRun.sh "
            "on a supported Apple-silicon host; no reference-provider "
            "fallback is performed.";
    return {};
#endif
}

} // namespace

int wavevortex::runtime::runWaveVortex(
    int argc, char** argv,
    std::shared_ptr<const WVExtensionCatalog> catalog) {
    const auto totalStart = Clock::now();
    Options options;
    std::string error;
    if (!catalog) {
        std::cerr << "The WaveVortex runner requires a frozen extension catalog.\n";
        return static_cast<int>(ExitCode::usage);
    }
    if (!parseOptions(argc,argv,options,error)) {
        emit(failureJSON(ExitCode::usage,"arguments",error),options.report,std::cerr);
        return static_cast<int>(ExitCode::usage);
    }

    Timings timings;
    timings.parse = seconds(totalStart);
    const auto preflightStart = Clock::now();
    phase(options,"inspect");
    WVCheckpointInspection inspection;
    WVModelOutputNetCDFInspection modelInspection;
    auto start = Clock::now();
    WVCheckpointStatus checkpointStatus;
    if (options.restartMode == "model") {
        checkpointStatus = WVModelOutputNetCDFSink::inspect(
            options.modelFiles, *catalog, modelInspection);
        if (checkpointStatus)
            inspection = modelInspection.latestRestart;
    } else {
        checkpointStatus = WVCheckpointReader::inspect(options.input,*catalog,inspection);
    }
    timings.inspect = seconds(start);
    if (!checkpointStatus) {
        emit(failureJSON(ExitCode::checkpoint,
                         options.restartMode == "model"
                             ? "model-graph-preflight" : "inspect",
                         checkpointStatus.message,checkpointStatus.location),
             options.report,std::cerr);
        return static_cast<int>(ExitCode::checkpoint);
    }
    if (options.restartMode == "coefficients" && !options.scheduledOutput()) {
        const auto inputPath = normalizedPath(options.input);
        const auto outputPath = normalizedPath(options.output);
        if (inputPath == outputPath) {
            emit(failureJSON(ExitCode::usage,"output-preflight",
                             "Coefficient output must not alias its input."),
                 options.report,std::cerr);
            return static_cast<int>(ExitCode::usage);
        }
        std::error_code destinationError;
        const auto destinationStatus =
            std::filesystem::symlink_status(outputPath,destinationError);
        if (destinationError == std::errc::no_such_file_or_directory)
            destinationError.clear();
        if (destinationError) {
            emit(failureJSON(ExitCode::output,"output-preflight",
                             "Unable to inspect output destination: " +
                                 destinationError.message(),
                             outputPath.string()),options.report,std::cerr);
            return static_cast<int>(ExitCode::output);
        }
        if (options.outputPolicy == "create" &&
            destinationStatus.type() != std::filesystem::file_type::not_found) {
            emit(failureJSON(ExitCode::output,"output-preflight",
                             "Output destination already exists; use "
                             "--output-policy replace to authorize replacement.",
                             outputPath.string()),options.report,std::cerr);
            return static_cast<int>(ExitCode::output);
        }
    }
    auto activeForcingSchedule = inspection.forcingSchedule;
    if (options.restartMode == "model" && modelInspection.isDynamicsLinear)
        activeForcingSchedule.entries.clear();
    const auto forcingStatus = WVModel::validateCheckpointForcingSchedule(
        inspection, activeForcingSchedule, *catalog);
    if (!forcingStatus) {
        emit(failureJSON(ExitCode::checkpoint,"forcing-preflight",forcingStatus.message),options.report,std::cerr);
        return static_cast<int>(ExitCode::checkpoint);
    }
    if (options.hasFinalTime && options.finalTime < inspection.t) {
        emit(failureJSON(ExitCode::usage,"arguments","--final-time precedes the selected checkpoint state."),options.report,std::cerr);
        return static_cast<int>(ExitCode::usage);
    }
    if (options.stepPolicy == cli::WVRunRequestStepPolicy::cflSelected &&
        !inspection.stateDescription.supportsFixedTimeStepSelection) {
        emit(failureJSON(ExitCode::usage,"integration-preflight",
                         "The selected numerical system does not provide "
                         "fixed-step CFL candidates."),
             options.report,std::cerr);
        return static_cast<int>(ExitCode::usage);
    }
    const auto requestedInterval = options.hasFinalTime
        ? options.finalTime - inspection.t
        : static_cast<double>(options.steps) * options.deltaT;
    if (options.stepPolicy == cli::WVRunRequestStepPolicy::adaptive &&
        !options.hasMaximumStep && !(requestedInterval > 0.0)) {
        emit(failureJSON(
                 ExitCode::usage,"integration-preflight",
                 "The default adaptive maximum step requires a positive "
                 "continuation interval."),
             options.report,std::cerr);
        return static_cast<int>(ExitCode::usage);
    }
    const auto defaultMaximumStep = options.hasFinalTime && requestedInterval > 0.0
        ? 0.1 * requestedInterval
        : options.deltaT;
    const auto requestedMaximumStep =
        options.hasMaximumStep ? options.maximumStep : defaultMaximumStep;
    const auto effectiveMaximumStep = options.hasFinalTime && requestedInterval > 0.0
        ? std::min(requestedMaximumStep, requestedInterval)
        : requestedMaximumStep;
    ScheduledOutputPlan scheduledPlan;
    if (options.scheduledOutput()) {
        ExitCode failureCode = ExitCode::usage;
        if (!prepareScheduledOutput(options,inspection.t,scheduledPlan,failureCode,error)) {
            emit(failureJSON(failureCode,"output-preflight",error,options.outputDirectory,scheduledOutputJSON(options,scheduledPlan,nullptr)),options.report,std::cerr);
            return static_cast<int>(failureCode);
        }
    }

    WVModelOutputConfiguration preparedOutputConfiguration;
    if (options.restartMode == "model") {
        WVModelOutputRequest outputRequest;
        outputRequest.policy = options.outputPolicy == "create"
            ? WVModelOutputPolicy::create
            : options.outputPolicy == "replace"
                ? WVModelOutputPolicy::replace
                : WVModelOutputPolicy::append;
        outputRequest.finalTime = options.finalTime;
        outputRequest.destinations = options.outputDestinations;
        const auto outputStatus = WVModel::prepareModelOutput(
            catalog,modelInspection,outputRequest,preparedOutputConfiguration);
        if (!outputStatus) {
            emit(failureJSON(ExitCode::output,"model-output-preflight",
                             outputStatus.message),options.report,std::cerr);
            return static_cast<int>(ExitCode::output);
        }
    }
    timings.preflight = seconds(preflightStart);

    std::string providerVersion;
    std::string baseLibrary;
    std::string threadLibrary;
    start = Clock::now();
    auto fftEngine = provider(options,providerVersion,baseLibrary,threadLibrary,error);
    timings.provider = seconds(start);
    if (!fftEngine) {
        emit(failureJSON(ExitCode::provider,"provider",error,{},options.scheduledOutput() ? scheduledOutputJSON(options,scheduledPlan,nullptr) : std::string{}),options.report,std::cerr);
        return static_cast<int>(ExitCode::provider);
    }

    WVCheckpoint sourceCheckpoint;
    phase(options,"read");
    start = Clock::now();
    checkpointStatus = options.restartMode == "model"
                           ? WVCheckpointStatus::ok()
                           : WVCheckpointReader::read(options.input,*catalog,sourceCheckpoint);
    timings.read = seconds(start);
    if (!checkpointStatus) {
        emit(failureJSON(ExitCode::checkpoint,"read",checkpointStatus.message,checkpointStatus.location,options.scheduledOutput() ? scheduledOutputJSON(options,scheduledPlan,nullptr) : std::string{}),options.report,std::cerr);
        return static_cast<int>(ExitCode::checkpoint);
    }

    WVKernelStatus kernelStatus;
    phase(options,"construct");
    start = Clock::now();
    WVModelIntegratorConfiguration integratorConfiguration;
    const auto retainDenseOutput =
        options.benchmarkDenseOutputsPerStep != 0 || options.scheduledOutput() ||
        options.benchmarkOutputCount != 0 || options.benchmarkModelDenseOutput ||
        (options.requestMode && options.restartMode == "model");
    if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK23) {
        integratorConfiguration.kind = WVModelIntegratorKind::adaptiveRK23;
        integratorConfiguration.adaptive.relativeTolerance = options.relativeTolerance;
        integratorConfiguration.adaptive.absoluteToleranceScale = options.absoluteTolerance;
        integratorConfiguration.adaptive.maximumStepSize = effectiveMaximumStep;
        integratorConfiguration.adaptive.retainDenseOutput = retainDenseOutput;
    } else if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK45) {
        integratorConfiguration.kind = WVModelIntegratorKind::adaptiveRK45;
        integratorConfiguration.adaptiveRK45.relativeTolerance = options.relativeTolerance;
        integratorConfiguration.adaptiveRK45.absoluteToleranceScale = options.absoluteTolerance;
        integratorConfiguration.adaptiveRK45.maximumStepSize = effectiveMaximumStep;
        integratorConfiguration.adaptiveRK45.retainDenseOutput = retainDenseOutput;
    } else if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK78) {
        integratorConfiguration.kind = WVModelIntegratorKind::adaptiveRK78;
        integratorConfiguration.adaptiveRK78.relativeTolerance = options.relativeTolerance;
        integratorConfiguration.adaptiveRK78.absoluteToleranceScale = options.absoluteTolerance;
        integratorConfiguration.adaptiveRK78.maximumStepSize = effectiveMaximumStep;
        integratorConfiguration.adaptiveRK78.retainDenseOutput = retainDenseOutput;
    } else {
        integratorConfiguration.kind = WVModelIntegratorKind::fixedRK4;
        integratorConfiguration.fixed.retainDenseOutput = retainDenseOutput;
    }
    WVModel model;
    WVModelState modelState;
    if (options.restartMode == "model") {
        kernelStatus = WVModel::createFromModelOutputInspection(
            catalog,std::move(modelInspection),
            std::move(preparedOutputConfiguration),std::move(fftEngine),
            integratorConfiguration,model,modelState);
    } else {
        sourceCheckpoint.forcingSchedule = activeForcingSchedule;
        kernelStatus = WVModel::createFromCheckpoint(
            catalog, std::move(sourceCheckpoint), std::move(fftEngine),
            integratorConfiguration, model, modelState);
    }
    if (!kernelStatus) {
        emit(failureJSON(ExitCode::provider,"construct",kernelStatus.message,{},options.scheduledOutput() ? scheduledOutputJSON(options,scheduledPlan,nullptr) : std::string{}),options.report,std::cerr);
        return static_cast<int>(ExitCode::provider);
    }
    auto &checkpoint = modelState.checkpoint();
    auto &integrator = detail::WVModelInternalAccess::integrator(model);
#if !WV_RUNTIME_HAS_DENSE_OUTPUT
    if (options.benchmarkDenseOutputsPerStep != 0) {
        emit(failureJSON(ExitCode::usage,"arguments","This archived baseline does not implement dense output."),options.report,std::cerr);
        return static_cast<int>(ExitCode::usage);
    }
#endif
    WVFixedStepRK4* fixedIntegrator = nullptr;
    WVAdaptiveRK23* adaptiveRK23Integrator = nullptr;
    WVAdaptiveRK45* adaptiveRK45Integrator = nullptr;
    WVAdaptiveRK78* adaptiveRK78Integrator = nullptr;
    if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK23) {
        adaptiveRK23Integrator = &static_cast<WVAdaptiveRK23&>(integrator);
    } else if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK45) {
        adaptiveRK45Integrator = &static_cast<WVAdaptiveRK45&>(integrator);
    } else if (options.integrator == cli::WVRunRequestIntegrationMethod::adaptiveRK78) {
        adaptiveRK78Integrator = &static_cast<WVAdaptiveRK78&>(integrator);
    } else {
        fixedIntegrator = &static_cast<WVFixedStepRK4&>(integrator);
    }
    double integrationInitialStep = fixedIntegrator == nullptr
        ? options.hasInitialStep
            ? std::min(options.initialStep,effectiveMaximumStep)
            : options.requestSchemaVersion ==
                      cli::WVRunRequest::schemaV2Version
                ? 0.0
                : std::min(options.deltaT,effectiveMaximumStep)
        : options.deltaT;
    double selectedFixedStep = integrationInitialStep;
    WVFixedTimeStepCandidates fixedTimeStepCandidates;
    bool evaluatedFixedTimeStepCandidates = false;
    const auto &coefficientDimensions =
        model.stateLayout().coefficientFamilies().front().spectralDimensions;
    auto state = modelState.mutableView();
    timings.construct = seconds(start);
    phase(options,"prepare");
    start = Clock::now();
    kernelStatus = model.prepareStateAfterRestart(modelState);
    state = modelState.mutableView();
    timings.prepare = seconds(start);
    if (!kernelStatus) {
        emit(failureJSON(ExitCode::integration,"prepare",kernelStatus.message,{},options.scheduledOutput() ? scheduledOutputJSON(options,scheduledPlan,nullptr) : std::string{}),options.report,std::cerr);
        return static_cast<int>(ExitCode::integration);
    }
    const bool resolveDefaultAdaptiveInitialStep =
        options.requestSchemaVersion == cli::WVRunRequest::schemaV2Version &&
        fixedIntegrator == nullptr && !options.hasInitialStep;
    if (options.stepPolicy == cli::WVRunRequestStepPolicy::cflSelected ||
        resolveDefaultAdaptiveInitialStep) {
        start = Clock::now();
        if (!model.supportsFixedTimeStepSelection()) {
            kernelStatus = {WVKernelStatusCode::unsupportedOperation,
                            "The resolved numerical system does not provide "
                            "fixed-step CFL candidates."};
        } else {
            kernelStatus = model.evaluateFixedTimeStepCandidates(
                modelState,
                resolveDefaultAdaptiveInitialStep ? 0.5 : options.cfl,
                fixedTimeStepCandidates);
        }
        timings.cfl = seconds(start);
        if (!kernelStatus) {
            emit(failureJSON(ExitCode::integration,"cfl-selection",
                             kernelStatus.message),options.report,std::cerr);
            return static_cast<int>(ExitCode::integration);
        }
        evaluatedFixedTimeStepCandidates = true;
        if (!resolveDefaultAdaptiveInitialStep &&
            options.timeStepConstraint ==
            cli::WVRunRequestTimeStepConstraint::advective)
            selectedFixedStep = fixedTimeStepCandidates.advective;
        else if (!resolveDefaultAdaptiveInitialStep &&
                 options.timeStepConstraint ==
                 cli::WVRunRequestTimeStepConstraint::oscillatory)
            selectedFixedStep = fixedTimeStepCandidates.oscillatory;
        else
            selectedFixedStep =
                std::min(fixedTimeStepCandidates.advective,
                         fixedTimeStepCandidates.oscillatory);
        if (!(selectedFixedStep > 0.0) || std::isnan(selectedFixedStep)) {
            emit(failureJSON(ExitCode::integration,"cfl-selection",
                             "The selected CFL step is not positive."),
                 options.report,std::cerr);
            return static_cast<int>(ExitCode::integration);
        }
        integrationInitialStep = std::isfinite(selectedFixedStep)
            ? selectedFixedStep
            : requestedInterval > 0.0 ? requestedInterval : 1.0;
        if (resolveDefaultAdaptiveInitialStep) {
            options.usedDefaultInitialStep = true;
            integrationInitialStep =
                std::min(integrationInitialStep,effectiveMaximumStep);
        }
        options.deltaT = integrationInitialStep;
    }
    timings.startup = timings.read + timings.construct + timings.prepare +
                      timings.cfl;

    phasePlateau(options,"steady-retained");
    const auto integrationBaselineRSS = currentRSSBytes();
#if WV_RUNTIME_HAS_DENSE_OUTPUT
    BenchmarkSink benchmarkSink;
    WVOutputDriverMetrics outputDriverMetrics;
    std::size_t outputPlanMaximumLiveBytes = 0;
    std::size_t outputOrchestrationMaximumLiveBytes = 0;
    std::unique_ptr<PhaseReportingCheckpointSink> scheduledSink;
    if (options.scheduledOutput()) {
        scheduledSink = std::make_unique<PhaseReportingCheckpointSink>(catalog,checkpoint,options);
    }
    const auto runOutput = [&](const std::vector<WVExplicitOutputTarget>& targets,
                               double finalTime,
                               WVOutputSink& sink) -> WVKernelStatus {
        WVOutputPlan plan;
        auto status = WVOutputPlan::createExplicit(model.stateLayout(),catalog,state.waveVortex.t,finalTime,targets,plan);
        if (!status) return status;
        const auto planBytes = plan.persistentBytes();
        outputPlanMaximumLiveBytes =
            std::max(outputPlanMaximumLiveBytes, planBytes);
        status = model.advanceToTime(modelState,finalTime,integrationInitialStep,
                                     plan,sink);
        state = modelState.mutableView();
        const auto metrics = model.metrics(&modelState).outputDriver;
        outputDriverMetrics.acceptedStepCount += metrics.acceptedStepCount;
        outputDriverMetrics.outputStateEvaluationCount += metrics.outputStateEvaluationCount;
        outputDriverMetrics.initialStateEventCount += metrics.initialStateEventCount;
        outputDriverMetrics.interpolatedStateEvaluationCount += metrics.interpolatedStateEvaluationCount;
        outputDriverMetrics.acceptedEndpointStateEventCount += metrics.acceptedEndpointStateEventCount;
        outputDriverMetrics.interpolationSeconds += metrics.interpolationSeconds;
        outputDriverMetrics.interpolationBufferCapacityBytes = std::max(outputDriverMetrics.interpolationBufferCapacityBytes,metrics.interpolationBufferCapacityBytes);
        outputDriverMetrics.interpolationBufferMaximumLiveBytes = std::max(outputDriverMetrics.interpolationBufferMaximumLiveBytes,metrics.interpolationBufferMaximumLiveBytes);
        outputDriverMetrics.routeStagingCapacityBytes = std::max(outputDriverMetrics.routeStagingCapacityBytes,metrics.routeStagingCapacityBytes);
        outputDriverMetrics.routeStagingMaximumLiveBytes = std::max(outputDriverMetrics.routeStagingMaximumLiveBytes,metrics.routeStagingMaximumLiveBytes);
        outputDriverMetrics.retainedStorageBytes = std::max(outputDriverMetrics.retainedStorageBytes,metrics.retainedStorageBytes);
        outputOrchestrationMaximumLiveBytes = std::max(
            outputOrchestrationMaximumLiveBytes,
            planBytes + metrics.retainedStorageBytes);
        outputDriverMetrics.generatedSemanticOccurrenceCount += metrics.generatedSemanticOccurrenceCount;
        return status;
    };
    const auto benchmarkTargets = [](const std::vector<double>& times) {
        std::vector<WVExplicitOutputTarget> targets;
        targets.reserve(times.size());
        for (std::size_t index = 0; index < times.size(); ++index)
            targets.push_back({times[index],"benchmark-output-"+std::to_string(index+1)});
        return targets;
    };
#endif
    double proposedStepSize = integrationInitialStep;
    const auto advanceBenchmarkSteps = [&](std::size_t count) -> WVKernelStatus {
        if (count == 0) return WVKernelStatus::ok();
        if (options.benchmarkDenseOutputsPerStep == 0) {
            for (std::size_t step = 0; step < count; ++step) {
                const auto stepStatus = model.step(modelState,proposedStepSize);
                if (!stepStatus) return stepStatus;
                state = modelState.mutableView();
                proposedStepSize = model.nextStepSize();
            }
            return WVKernelStatus::ok();
        }
#if WV_RUNTIME_HAS_DENSE_OUTPUT
        const auto finalTime = state.waveVortex.t+static_cast<double>(count)*options.deltaT;
        return runOutput(benchmarkTargets(interiorOutputTimes(state.waveVortex.t,options.deltaT,count,options.benchmarkDenseOutputsPerStep)),finalTime,benchmarkSink);
#else
        return {WVKernelStatusCode::unsupportedOperation,"This archived baseline does not implement dense output."};
#endif
    };
    if (options.benchmarkWarmupSteps != 0) {
        kernelStatus = advanceBenchmarkSteps(options.benchmarkWarmupSteps);
        if (!kernelStatus) {
            emit(failureJSON(ExitCode::integration,"warmup",kernelStatus.message),options.report,std::cerr);
            return static_cast<int>(ExitCode::integration);
        }
    }
    phase(options,"integrate");
    if (!waitForPhaseSample("integrate")) {
        emit(failureJSON(ExitCode::integration,"memory-phase","The external RSS sampler did not acknowledge the integration phase."),options.report,std::cerr);
        return static_cast<int>(ExitCode::integration);
    }
    start = Clock::now();
    if (options.benchmarkWarmupSteps != 0 || options.benchmarkDenseOutputsPerStep != 0) {
        kernelStatus = advanceBenchmarkSteps(options.steps);
    } else if (options.scheduledOutput()) {
#if WV_RUNTIME_HAS_DENSE_OUTPUT
        kernelStatus = runOutput(scheduledPlan.targets,options.finalTime,*scheduledSink);
#else
        kernelStatus = {WVKernelStatusCode::unsupportedOperation,"This archived baseline does not implement scheduled output."};
#endif
    } else if (options.restartMode == "model") {
        kernelStatus = model.advanceToTime(modelState,options.finalTime,
                                           integrationInitialStep);
        state = modelState.mutableView();
        outputDriverMetrics = model.metrics(&modelState).outputDriver;
        outputOrchestrationMaximumLiveBytes =
            outputDriverMetrics.retainedStorageBytes;
    } else if (options.benchmarkOutputCount != 0) {
#if WV_RUNTIME_HAS_DENSE_OUTPUT
        kernelStatus = runOutput(benchmarkTargets(uniformInteriorOutputTimes(state.waveVortex.t,options.finalTime,options.benchmarkOutputCount)),options.finalTime,benchmarkSink);
#else
        kernelStatus = {WVKernelStatusCode::unsupportedOperation,"This archived baseline does not implement dense output."};
#endif
    } else if (options.hasSteps) {
        for (std::size_t step = 0; step < options.steps && kernelStatus; ++step) {
            kernelStatus = model.step(modelState,proposedStepSize);
            if (kernelStatus) {
                state = modelState.mutableView();
                proposedStepSize = model.nextStepSize();
            }
        }
    } else {
        kernelStatus = model.advanceToTime(modelState,options.finalTime,
                                           integrationInitialStep);
        state = modelState.mutableView();
    }
    timings.integrate = seconds(start);
    const auto integrationPeakRSS = peakRSSBytes();
    if (!kernelStatus) {
        emit(failureJSON(options.scheduledOutput() ? ExitCode::output : ExitCode::integration,options.scheduledOutput() ? "scheduled-output" : "integrate",kernelStatus.message,{},options.scheduledOutput() ? scheduledOutputJSON(options,scheduledPlan,&scheduledSink->sink()) : std::string{}),options.report,std::cerr);
        return static_cast<int>(options.scheduledOutput() ? ExitCode::output : ExitCode::integration);
    }
    if (options.scheduledOutput() && scheduledSink->sink().metrics().checkpointWriteCount != scheduledPlan.targets.size()) {
        emit(failureJSON(ExitCode::output,"scheduled-output","Integration completed before every requested checkpoint was written.",{},scheduledOutputJSON(options,scheduledPlan,&scheduledSink->sink())),options.report,std::cerr);
        return static_cast<int>(ExitCode::output);
    }
    phasePlateau(options,"outputs-held");

    if (options.restartMode == "model") {
        start = Clock::now();
        checkpointStatus = model.closeOutput();
        timings.write = seconds(start);
    } else if (options.scheduledOutput()) {
        timings.write = scheduledSink->sink().metrics().checkpointWriteSeconds;
    } else {
        checkpoint.state.t = state.waveVortex.t;
        checkpoint.state.t0 = state.waveVortex.t0;
        start = Clock::now();
        phase(options,"write");
        checkpointStatus = WVCheckpointWriter::write(
            options.output,*catalog,checkpoint,
            options.outputPolicy == "replace"
                ? WVCheckpointCommitPolicy::replaceExisting
                : WVCheckpointCommitPolicy::createNew);
        timings.write = seconds(start);
    }
    timings.total = seconds(totalStart);
    if (!options.scheduledOutput() && !checkpointStatus) {
        emit(failureJSON(ExitCode::output,"write",checkpointStatus.message,checkpointStatus.location),options.report,std::cerr);
        return static_cast<int>(ExitCode::output);
    }

    const auto modelMetrics = model.metrics(&modelState);
    const auto& kernelMetrics = modelMetrics.kernel;
    const auto& forcingMetrics = modelMetrics.forcing;
    WVFixedStepRK4Metrics fixedMetrics;
    WVAdaptiveRK23Metrics adaptiveMetrics;
    if (fixedIntegrator != nullptr) {
        fixedMetrics = fixedIntegrator->metrics();
        adaptiveMetrics = fixedMetrics;
    }
    const std::vector<WVAdaptiveRKStepDiagnostic>* adaptiveDiagnostics = nullptr;
    const std::vector<std::uint64_t>* adaptiveToleranceComponentHashes = nullptr;
    const WVAdaptiveRKStageBufferLastUse* adaptiveStageBufferLastUse = nullptr;
    std::size_t adaptiveStageBufferLastUseCount = 0;
    std::uint64_t adaptiveToleranceHash = 0;
    const char* adaptiveController = "fixed-rk4";
    bool adaptiveDiagnosticsComplete = true;
    if (adaptiveRK23Integrator != nullptr) {
        adaptiveMetrics = adaptiveRK23Integrator->metrics();
        adaptiveDiagnostics = &adaptiveRK23Integrator->stepDiagnostics();
        adaptiveToleranceHash = adaptiveRK23Integrator->toleranceHash();
        adaptiveToleranceComponentHashes = &adaptiveRK23Integrator->toleranceComponentHashes();
        adaptiveController = WVAdaptiveRK23::controllerIdentifier();
        adaptiveDiagnosticsComplete = adaptiveRK23Integrator->stepDiagnosticsComplete();
        adaptiveStageBufferLastUse = WVAdaptiveRK23::stageBufferLastUseRecords();
        adaptiveStageBufferLastUseCount = WVAdaptiveRK23::stageBufferLastUseRecordCount();
    } else if (adaptiveRK45Integrator != nullptr) {
        adaptiveMetrics = adaptiveRK45Integrator->metrics();
        adaptiveDiagnostics = &adaptiveRK45Integrator->stepDiagnostics();
        adaptiveToleranceHash = adaptiveRK45Integrator->toleranceHash();
        adaptiveToleranceComponentHashes = &adaptiveRK45Integrator->toleranceComponentHashes();
        adaptiveController = WVAdaptiveRK45::controllerIdentifier();
        adaptiveDiagnosticsComplete = adaptiveRK45Integrator->stepDiagnosticsComplete();
        adaptiveStageBufferLastUse = WVAdaptiveRK45::stageBufferLastUseRecords();
        adaptiveStageBufferLastUseCount = WVAdaptiveRK45::stageBufferLastUseRecordCount();
    } else if (adaptiveRK78Integrator != nullptr) {
        adaptiveMetrics = adaptiveRK78Integrator->metrics();
        adaptiveDiagnostics = &adaptiveRK78Integrator->stepDiagnostics();
        adaptiveToleranceHash = adaptiveRK78Integrator->toleranceHash();
        adaptiveToleranceComponentHashes = &adaptiveRK78Integrator->toleranceComponentHashes();
        adaptiveController = WVAdaptiveRK78::controllerIdentifier();
        adaptiveDiagnosticsComplete = adaptiveRK78Integrator->stepDiagnosticsComplete();
        adaptiveStageBufferLastUse = WVAdaptiveRK78::stageBufferLastUseRecords();
        adaptiveStageBufferLastUseCount = WVAdaptiveRK78::stageBufferLastUseRecordCount();
    }
    if (fixedIntegrator != nullptr) {
        adaptiveStageBufferLastUse = WVFixedStepRK4::stageBufferLastUseRecords();
        adaptiveStageBufferLastUseCount = WVFixedStepRK4::stageBufferLastUseRecordCount();
    }
    const auto hasAdaptiveIntegrator = adaptiveDiagnostics != nullptr;
    const auto& integratorMetrics = fixedIntegrator != nullptr ? fixedMetrics : adaptiveMetrics;
    const auto& integratedObserverMetrics = modelMetrics.integratedObservers;
    const auto outputEvaluationMetrics = modelMetrics.outputEvaluation;
    const auto modelOutputMetrics = modelMetrics.output;
    const auto stepCount = fixedIntegrator != nullptr ? fixedMetrics.stepCount : adaptiveMetrics.acceptedStepCount;
    const auto rejectedStepCount = hasAdaptiveIntegrator ? adaptiveMetrics.rejectedStepCount : 0;
    const auto rightHandSideEvaluationCount = fixedIntegrator != nullptr ? fixedMetrics.rightHandSideEvaluationCount : adaptiveMetrics.rightHandSideEvaluationCount;
    const auto integratorWorkspaceCapacityBytes = integratorMetrics.workspaceCapacityBytes;
    const auto integratorWorkspaceLiveBytes = integratorMetrics.workspaceLiveBytes;
    const auto integratorWorkspaceMaximumLiveBytes = integratorMetrics.workspaceMaximumLiveBytes;
    const auto integratorStateSizedWorkspaceBytes =
        integratorMetrics.stateCapacityBytes *
        integratorMetrics.workspaceStateEquivalentCount;
    const auto integratorWorkspaceNonStateBytes =
        integratorWorkspaceCapacityBytes >= integratorStateSizedWorkspaceBytes
            ? integratorWorkspaceCapacityBytes -
                  integratorStateSizedWorkspaceBytes
            : 0;
    const auto integratorByteLedgerAgrees =
        integratorWorkspaceCapacityBytes >= integratorStateSizedWorkspaceBytes &&
        integratorWorkspaceCapacityBytes ==
            integratorStateSizedWorkspaceBytes +
                integratorWorkspaceNonStateBytes;
    const auto integratorDiagnosticBytes = adaptiveMetrics.diagnosticCapacityBytes;
    const auto continuousExtensionBytes =
        adaptiveMetrics.continuousExtensionWorkspaceCapacityBytes;
#if WV_RUNTIME_HAS_DENSE_OUTPUT
    const auto denseHistoryBytes = fixedIntegrator != nullptr ? fixedMetrics.denseHistoryCapacityBytes : adaptiveMetrics.denseHistoryCapacityBytes;
    const auto driverInterpolationMaximumLiveBytes = outputDriverMetrics.interpolationBufferMaximumLiveBytes;
    const auto outputDriverMaximumLiveBytes = outputDriverMetrics.retainedStorageBytes;
    const auto driverInterpolationSeconds = outputDriverMetrics.interpolationSeconds;
    const auto interpolatedOutputCount = options.scheduledOutput() ? outputDriverMetrics.interpolatedStateEvaluationCount : benchmarkSink.interpolatedCount;
    const auto scheduledOutputBytes =
        scheduledSink == nullptr
            ? 0
            : sizeof(PhaseReportingCheckpointSink) -
                  sizeof(WVCheckpointOutputSink) +
                  scheduledSink->sink().persistentBytes();
#else
    const std::size_t denseHistoryBytes = 0;
    const std::size_t driverInterpolationMaximumLiveBytes = 0;
    const std::size_t outputDriverMaximumLiveBytes = 0;
    const std::size_t outputPlanMaximumLiveBytes = 0;
    const std::size_t outputOrchestrationMaximumLiveBytes = 0;
    const double driverInterpolationSeconds = 0.0;
    const std::size_t interpolatedOutputCount = 0;
    const std::size_t scheduledOutputBytes = 0;
#endif
    const auto checkpointStateBytes =
        checkpointCoefficientStorageBytes(checkpoint);
    const auto knownPersistentBytes =
        modelMetrics.modelPersistentBytes + modelMetrics.catalogPersistentBytes +
        checkpointStateBytes +
        modelMetrics.integrationSystemPersistentBytes +
        modelMetrics.integratorPersistentBytes +
        modelMetrics.outputConfigurationPersistentBytes + scheduledOutputBytes;
    const auto fullModelPersistentBytes =
        knownPersistentBytes +
        (modelMetrics.statePersistentBytes > checkpointStateBytes
             ? modelMetrics.statePersistentBytes - checkpointStateBytes
             : 0) +
        modelMetrics.outputEvaluationPersistentBytes +
        modelMetrics.outputSinkPersistentBytes;
    const auto occurrenceWorkspaceRetainedBytes =
        modelOutputMetrics.occurrenceWorkspaceRetainedBytes;
    const auto occurrenceWorkspaceMaximumLiveBytes =
        modelOutputMetrics.occurrenceWorkspaceMaximumLiveBytes;
    const auto staticFullModelPersistentBytes =
        fullModelPersistentBytes > occurrenceWorkspaceRetainedBytes
            ? fullModelPersistentBytes - occurrenceWorkspaceRetainedBytes
            : 0;
    const auto fullModelRetainedBytes =
        staticFullModelPersistentBytes + occurrenceWorkspaceRetainedBytes;
    const auto fullModelMaximumLiveBytes =
        staticFullModelPersistentBytes +
        std::max(occurrenceWorkspaceRetainedBytes,
                 occurrenceWorkspaceMaximumLiveBytes) +
        outputOrchestrationMaximumLiveBytes;
    const auto integratorElementReads = fixedMetrics.stageStateConstructionElementReads+fixedMetrics.weightedFluxInitializationElementReads+fixedMetrics.weightedAccumulationElementReads+fixedMetrics.finalStateUpdateElementReads+fixedMetrics.acceptedStateCommitElementReads;
    const auto integratorElementWrites = fixedMetrics.stageStateConstructionElementWrites+fixedMetrics.stageFluxClearElementWrites+fixedMetrics.weightedFluxClearElementWrites+fixedMetrics.weightedFluxInitializationElementWrites+fixedMetrics.weightedAccumulationElementWrites+fixedMetrics.finalStateUpdateElementWrites+fixedMetrics.acceptedStateCommitElementWrites;
    const auto forcingElementReads = forcingMetrics.temporaryAccumulationElementReads+forcingMetrics.outputCopyElementReads;
    const auto forcingElementWrites = forcingMetrics.accumulatorClearElementWrites+forcingMetrics.temporaryFluxClearElementWrites+forcingMetrics.kernelOutputInitializationElementWrites+forcingMetrics.temporaryAccumulationElementWrites+forcingMetrics.outputCopyElementWrites+forcingMetrics.stateConstraintElementWrites;
    const auto trafficElementReads = integratorElementReads+forcingElementReads;
    const auto trafficElementWrites = integratorElementWrites+forcingElementWrites;
    const auto requestedStepJSON = [&](bool present, double value) {
        if (!present) return std::string{"null"};
        return jsonNumberOrInfinity(value);
    };
    std::ostringstream providerReport;
    providerReport << "{";
    if (options.requestSchemaVersion == cli::WVRunRequest::schemaV2Version) {
        const std::string requestedProviderJSON = options.hasRequestedProvider
            ? std::string{"\""}+jsonEscape(options.provider)+"\""
            : "null";
        const std::string requestedThreadsJSON = options.hasRequestedThreads
            ? std::to_string(options.threads)
            : "null";
        providerReport << "\"requestedId\":" << requestedProviderJSON << ',';
        providerReport << "\"requestedThreads\":" << requestedThreadsJSON << ',';
    }
    providerReport << "\"id\":" << quoted(options.provider)
                   << ",\"version\":" << quoted(providerVersion)
                   << ",\"threads\":" << options.threads
                   << ",\"baseLibrary\":" << quoted(baseLibrary)
                   << ",\"threadLibrary\":" << quoted(threadLibrary) << '}';
    std::ostringstream report;
    report << std::setprecision(17)
           << "{\"schemaVersion\":\"wave-vortex-run-v1\",\"status\":\"complete\",\"source\":{\"commit\":" << quoted(WV_RUNTIME_SOURCE_COMMIT) << "},"
           << "\"input\":" << quoted(options.input) << ",\"output\":" << quoted(options.output) << ",\"restartMode\":" << quoted(options.restartMode) << ",\"outputPolicy\":" << quoted(options.outputPolicy) << ",\"provider\":" << providerReport.str() << ','
           << "\"request\":{\"active\":" << (options.requestMode ? "true" : "false") << ",\"path\":" << quoted(options.requestPath) << ",\"schemaIdentifier\":" << quoted(options.requestSchemaIdentifier) << ",\"schemaVersion\":" << options.requestSchemaVersion << ",\"modelFiles\":" << stringArrayJSON(options.modelFiles) << ",\"destinations\":" << destinationMapJSON(options.outputDestinations) << "},"
           << "\"state\":{\"initialTime\":" << inspection.t << ",\"finalTime\":" << state.waveVortex.t << ",\"deltaT\":" << options.deltaT << ",\"stepCount\":" << stepCount << ",\"rejectedStepCount\":" << rejectedStepCount << ",\"rhsEvaluationCount\":" << rightHandSideEvaluationCount << ",\"shape\":" << sizeArrayJSON(coefficientDimensions) << "},"
           << "\"integrator\":{\"id\":" << quoted(options.integrator) << ",\"controller\":" << quoted(adaptiveController) << ",\"relativeTolerance\":" << options.relativeTolerance << ",\"absoluteTolerance\":" << options.absoluteTolerance << ",\"requestedInitialStep\":" << requestedStepJSON(options.hasInitialStep,options.initialStep) << ",\"effectiveInitialStep\":" << integrationInitialStep << ",\"requestedMaximumStep\":" << requestedStepJSON(options.hasMaximumStep,options.maximumStep) << ",\"effectiveMaximumStep\":" << (hasAdaptiveIntegrator ? effectiveMaximumStep : options.deltaT) << ",\"toleranceHash\":" << quoted(hasAdaptiveIntegrator ? std::to_string(adaptiveToleranceHash) : "") << ",\"toleranceHashClearedMantissaBits\":20,\"toleranceComponentHashes\":" << (hasAdaptiveIntegrator ? unsignedIntegerArrayJSON(*adaptiveToleranceComponentHashes) : "[]") << ",\"lastNormalizedError\":" << adaptiveMetrics.normalizedError << ",\"lastProposedStepSize\":" << adaptiveMetrics.lastProposedStepSize << ",\"lastAcceptedStepSize\":" << adaptiveMetrics.lastAcceptedStepSize << ",\"nextStepSize\":" << integrator.nextStepSize() << ",\"fsalReuseCount\":" << adaptiveMetrics.fsalReuseCount << ",\"fsalInvalidationCount\":" << adaptiveMetrics.fsalInvalidationCount << ",\"rejectedInitialDerivativeReuseCount\":" << adaptiveMetrics.rejectedInitialDerivativeReuseCount << ",\"constraintModifiedCoefficientCount\":" << adaptiveMetrics.constraintModifiedCoefficientCount << ",\"baseRightHandSideEvaluationCount\":" << adaptiveMetrics.baseRightHandSideEvaluationCount << ",\"continuousExtensionRightHandSideEvaluationCount\":" << adaptiveMetrics.continuousExtensionRightHandSideEvaluationCount << ",\"denseOutputEvaluationCount\":" << (fixedIntegrator != nullptr ? fixedMetrics.denseOutputEvaluationCount : adaptiveMetrics.denseOutputEvaluationCount) << ",\"denseOutputCacheBuildCount\":" << adaptiveMetrics.denseOutputCacheBuildCount << ",\"denseOutputCacheReuseCount\":" << adaptiveMetrics.denseOutputCacheReuseCount << ",\"denseOutputElementReads\":" << (fixedIntegrator != nullptr ? fixedMetrics.denseOutputElementReads : adaptiveMetrics.denseOutputElementReads) << ",\"denseOutputElementWrites\":" << (fixedIntegrator != nullptr ? fixedMetrics.denseOutputElementWrites : adaptiveMetrics.denseOutputElementWrites) << ",\"denseOutputSeconds\":" << (fixedIntegrator != nullptr ? fixedMetrics.denseOutputSeconds : adaptiveMetrics.denseOutputSeconds) << ",\"continuousExtensionSeconds\":" << adaptiveMetrics.continuousExtensionSeconds << ",\"stateCapacityBytes\":" << adaptiveMetrics.stateCapacityBytes << ",\"workspaceStateEquivalentCount\":" << adaptiveMetrics.workspaceStateEquivalentCount << ",\"workspaceMaximumLiveStateEquivalentCount\":" << adaptiveMetrics.workspaceMaximumLiveStateEquivalentCount << ",\"denseHistoryStateEquivalentCount\":" << adaptiveMetrics.denseHistoryStateEquivalentCount << ",\"retainedBaseStageCapacityBytes\":" << adaptiveMetrics.retainedBaseStageCapacityBytes << ",\"retainedBaseStageStateEquivalentCount\":" << adaptiveMetrics.retainedBaseStageStateEquivalentCount << ",\"continuousExtensionWorkspaceCapacityBytes\":" << adaptiveMetrics.continuousExtensionWorkspaceCapacityBytes << ",\"continuousExtensionWorkspaceMaximumLiveBytes\":" << adaptiveMetrics.continuousExtensionWorkspaceMaximumLiveBytes << ",\"continuousExtensionWorkspaceStateEquivalentCount\":" << adaptiveMetrics.continuousExtensionWorkspaceStateEquivalentCount << ",\"continuousExtensionWorkspaceMaximumLiveStateEquivalentCount\":" << adaptiveMetrics.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount << ",\"errorPolicyBytes\":" << adaptiveMetrics.errorPolicyBytes << ",\"diagnosticBytes\":" << adaptiveMetrics.diagnosticCapacityBytes << ",\"stageBufferLastUse\":" << (hasAdaptiveIntegrator ? adaptiveStageBufferLastUseJSON(adaptiveStageBufferLastUse,adaptiveStageBufferLastUseCount) : "[]") << ",\"acceptedStepDiagnosticsComplete\":" << (adaptiveDiagnosticsComplete ? "true" : "false") << ",\"acceptedSteps\":" << (hasAdaptiveIntegrator ? adaptiveStepDiagnosticsJSON(*adaptiveDiagnostics) : "[]") << "},"
           << "\"integratorStateSizedBuffers\":" << adaptiveStageBufferLastUseJSON(adaptiveStageBufferLastUse,adaptiveStageBufferLastUseCount) << ','
           ;
    if (options.requestSchemaVersion == cli::WVRunRequest::schemaV2Version)
        report << "\"integrationRequest\":"
               << integrationRequestV2JSON(
                      options,fixedTimeStepCandidates,
                      evaluatedFixedTimeStepCandidates,adaptiveController,
                      selectedFixedStep,integrationInitialStep,
                      effectiveMaximumStep)
               << ',';
    report << "\"integratorStorageLedger\":{\"scope\":\"exact portable integrator-owned storage\",\"exact\":true,\"persistentBytes\":"
           << modelMetrics.integratorPersistentBytes
           << ",\"workspaceCapacityBytes\":"
           << integratorWorkspaceCapacityBytes
           << ",\"workspaceMaximumLiveBytes\":"
           << integratorWorkspaceMaximumLiveBytes
           << ",\"stateCapacityBytes\":" << integratorMetrics.stateCapacityBytes
           << ",\"stateSizedWorkspaceBytes\":"
           << integratorStateSizedWorkspaceBytes
           << ",\"workspaceNonStateBytes\":"
           << integratorWorkspaceNonStateBytes
           << ",\"denseHistoryBytes\":" << denseHistoryBytes
           << ",\"continuousExtensionBytes\":"
           << continuousExtensionBytes
           << ",\"sharedAbstractionStateSizedCopyCount\":0"
           << ",\"byteLedgerAgreement\":"
           << (integratorByteLedgerAgrees ? "true" : "false") << "},";
    report << "\"authorBenchmark\":{\"warmupStepCount\":" << options.benchmarkWarmupSteps << ",\"sampleStepCount\":" << (options.hasSteps ? options.steps : 0) << ",\"denseOutputsPerStep\":" << options.benchmarkDenseOutputsPerStep << ",\"requestedOutputCount\":" << options.benchmarkOutputCount << ",\"modelDenseOutputEnabled\":" << (options.benchmarkModelDenseOutput ? "true" : "false") << ",\"interpolatedOutputCount\":" << interpolatedOutputCount << ",\"interpolationSeconds\":" << driverInterpolationSeconds << "},"
           << "\"forcing\":" << forcingJSON(activeForcingSchedule) << ','
           << "\"timingSeconds\":" << timingJSON(options,timings) << ','
           << "\"storageAccounting\":{\"scope\":\"application-and-provider-reported-cpp-storage\",\"excludes\":[\"allocator-metadata\",\"opaque-fftw-plan-internals\",\"opaque-netcdf-library-internals\"]},"
           << "\"storageBytes\":{\"checkpointState\":" << checkpointStateBytes << ",\"modelFacade\":" << modelMetrics.modelPersistentBytes << ",\"modelState\":" << modelMetrics.statePersistentBytes << ",\"extensionCatalog\":" << modelMetrics.catalogPersistentBytes << ",\"integrationSystem\":" << modelMetrics.integrationSystemPersistentBytes << ",\"integratorPersistent\":" << modelMetrics.integratorPersistentBytes << ",\"modelOutputConfiguration\":" << modelMetrics.outputConfigurationPersistentBytes << ",\"modelOutputEvaluation\":" << modelMetrics.outputEvaluationPersistentBytes << ",\"modelOutputSink\":" << modelMetrics.outputSinkPersistentBytes << ",\"modelOutput\":" << modelMetrics.outputPersistentBytes << ",\"descriptor\":" << kernelMetrics.descriptorBytes << ",\"planWrapper\":" << kernelMetrics.planBytes << ",\"kernelEngine\":" << kernelMetrics.engineBytes << ",\"kernelManagement\":" << kernelMetrics.kernelManagementBytes << ",\"kernelScratch\":" << kernelMetrics.scratchCapacityBytes << ",\"forcingSchedule\":" << forcingMetrics.scheduleBytes << ",\"forcingDerivedOperators\":" << forcingMetrics.derivedOperatorBytes << ",\"forcingWorkspace\":" << forcingMetrics.workspaceCapacityBytes << ",\"integratorWorkspace\":" << integratorWorkspaceCapacityBytes << ",\"integratorDiagnostics\":" << integratorDiagnosticBytes << ",\"denseHistory\":" << denseHistoryBytes << ",\"driverInterpolationMaximumLive\":" << driverInterpolationMaximumLiveBytes << ",\"occurrenceWorkspace\":" << occurrenceWorkspaceRetainedBytes << ",\"scheduledOutput\":" << scheduledOutputBytes << ",\"knownPersistent\":" << knownPersistentBytes << ",\"fullModelPersistent\":" << fullModelRetainedBytes << ",\"persistentFullHermitian\":0},"
           << "\"arrayTraffic\":{\"scope\":\"exact fixed-RK4 integration-boundary arrays; adaptive traffic is reported through method metrics\",\"elementBytes\":" << sizeof(WVComplex64) << ",\"integrator\":{\"stageStateConstructionReads\":" << fixedMetrics.stageStateConstructionElementReads << ",\"stageStateConstructionWrites\":" << fixedMetrics.stageStateConstructionElementWrites << ",\"stageFluxClearWrites\":" << fixedMetrics.stageFluxClearElementWrites << ",\"weightedFluxClearWrites\":" << fixedMetrics.weightedFluxClearElementWrites << ",\"weightedFluxInitializationReads\":" << fixedMetrics.weightedFluxInitializationElementReads << ",\"weightedFluxInitializationWrites\":" << fixedMetrics.weightedFluxInitializationElementWrites << ",\"weightedAccumulationReads\":" << fixedMetrics.weightedAccumulationElementReads << ",\"weightedAccumulationWrites\":" << fixedMetrics.weightedAccumulationElementWrites << ",\"finalStateUpdateReads\":" << fixedMetrics.finalStateUpdateElementReads << ",\"finalStateUpdateWrites\":" << fixedMetrics.finalStateUpdateElementWrites << ",\"acceptedStateCommitReads\":" << fixedMetrics.acceptedStateCommitElementReads << ",\"acceptedStateCommitWrites\":" << fixedMetrics.acceptedStateCommitElementWrites << "},\"forcing\":{\"accumulatorClearWrites\":" << forcingMetrics.accumulatorClearElementWrites << ",\"temporaryFluxClearWrites\":" << forcingMetrics.temporaryFluxClearElementWrites << ",\"kernelOutputInitializationWrites\":" << forcingMetrics.kernelOutputInitializationElementWrites << ",\"temporaryAccumulationReads\":" << forcingMetrics.temporaryAccumulationElementReads << ",\"temporaryAccumulationWrites\":" << forcingMetrics.temporaryAccumulationElementWrites << ",\"outputCopyReads\":" << forcingMetrics.outputCopyElementReads << ",\"outputCopyWrites\":" << forcingMetrics.outputCopyElementWrites << ",\"stateConstraintWrites\":" << forcingMetrics.stateConstraintElementWrites << "},\"totals\":{\"elementReads\":" << trafficElementReads << ",\"elementWrites\":" << trafficElementWrites << ",\"bytesRead\":" << trafficElementReads*sizeof(WVComplex64) << ",\"bytesWritten\":" << trafficElementWrites*sizeof(WVComplex64) << "}},"
           << "\"forcingOperations\":{\"evaluationCount\":" << forcingMetrics.evaluationCount << ",\"physicalFieldReconstructionCount\":" << forcingMetrics.physicalFieldReconstructionCount << ",\"physicalFieldReuseCount\":" << forcingMetrics.physicalFieldReuseCount << ",\"spatialTendencyProjectionCount\":" << forcingMetrics.spatialTendencyProjectionCount << ",\"spatialTendencyClearElementWrites\":" << forcingMetrics.spatialTendencyClearElementWrites << "},"
           << "\"livenessBytes\":{\"integratorWorkspaceLive\":" << integratorWorkspaceLiveBytes << ",\"integratorWorkspaceMaximumLive\":" << integratorWorkspaceMaximumLiveBytes << ",\"forcingWorkspaceLive\":" << forcingMetrics.workspaceLiveBytes << ",\"forcingWorkspaceMaximumLive\":" << forcingMetrics.workspaceMaximumLiveBytes << ",\"denseHistoryRetainedWithinWorkspace\":" << denseHistoryBytes << ",\"acceptedStepAdditionalArrayStorage\":0,\"contractAbstractionAdditionalArrayStorage\":0,\"contractAbstractionMaximumLiveArrayStorage\":" << driverInterpolationMaximumLiveBytes << ",\"outputDriverRetained\":0,\"outputDriverMaximumLive\":" << outputDriverMaximumLiveBytes << ",\"outputPlanMaximumLive\":" << outputPlanMaximumLiveBytes << ",\"outputOrchestrationMaximumLive\":" << outputOrchestrationMaximumLiveBytes << ",\"occurrenceWorkspaceRetained\":" << occurrenceWorkspaceRetainedBytes << ",\"occurrenceWorkspaceMaximumLive\":" << occurrenceWorkspaceMaximumLiveBytes << ",\"knownRetained\":" << knownPersistentBytes << ",\"knownMaximumLive\":" << knownPersistentBytes+outputOrchestrationMaximumLiveBytes << ",\"fullModelRetained\":" << fullModelRetainedBytes << ",\"fullModelMaximumLive\":" << fullModelMaximumLiveBytes << "},"
           << "\"rssBytes\":{\"integrationBaseline\":" << integrationBaselineRSS << ",\"processPeak\":" << integrationPeakRSS << ",\"peakIncrementLowerBound\":" << (integrationPeakRSS > integrationBaselineRSS ? integrationPeakRSS-integrationBaselineRSS : 0) << "},"
           << "\"integrationBreakdownSeconds\":{\"rightHandSide\":" << integratedObserverMetrics.rightHandSideSeconds << ",\"waveVortexFlux\":" << integratedObserverMetrics.waveVortexFluxSeconds << ",\"additionalStateClear\":" << integratedObserverMetrics.additionalStateClearSeconds << ",\"tracerAdvection\":" << integratedObserverMetrics.tracerAdvectionSeconds << ",\"tracerForward\":" << kernelMetrics.scalarForwardSeconds << ",\"tracerDerivativeAssembly\":" << kernelMetrics.scalarDerivativeAssemblySeconds << ",\"tracerVerticalDerivative\":" << kernelMetrics.scalarVerticalDerivativeSeconds << ",\"tracerInverse\":" << kernelMetrics.scalarInverseSeconds << ",\"tracerProduct\":" << kernelMetrics.scalarProductSeconds << ",\"tracerAntialias\":" << kernelMetrics.scalarAntialiasSeconds << ",\"particleAdvection\":" << integratedObserverMetrics.particleAdvectionSeconds << ",\"denseInterpolation\":" << driverInterpolationSeconds << ",\"observerEvaluation\":" << outputEvaluationMetrics.evaluationSeconds << ",\"outputPayloadWrite\":" << modelOutputMetrics.payloadWriteSeconds << ",\"outputSynchronization\":" << modelOutputMetrics.synchronizationSeconds << "},"
           << "\"execution\":{\"engine\":" << quoted(model.kernelProviderIdentifier()) << ",\"library\":" << quoted(model.kernelProviderLibraryIdentity()) << ",\"schedule\":" << quoted(model.forcingScheduleIdentifier()) << ",\"planCount\":" << kernelMetrics.planCount << ",\"noFallback\":true}";
    if (options.scheduledOutput()) report << ",\"scheduledOutput\":" << scheduledOutputJSON(options,scheduledPlan,&scheduledSink->sink());
    report << '}';
    emit(report.str(),options.report,std::cout);
    phase(options,"complete");
    return static_cast<int>(ExitCode::success);
}
