#include "WVTestExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <chrono>
#include <csignal>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/wait.h>
#include <unistd.h>
#endif

#ifndef WV_TEST_COMPACT_VARIABLE_POLICY
#define WV_TEST_COMPACT_VARIABLE_POLICY 0
#endif

namespace {

using namespace wavevortex::runtime;

void require(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }

std::string quote(const std::filesystem::path& path) { return "\""+path.string()+"\""; }

std::string explicitRestartMode(std::string arguments) {
    if (arguments.find("--restart-mode") == std::string::npos)
        arguments += " --restart-mode coefficients";
    if (arguments.find("--output-policy") == std::string::npos)
        arguments += arguments.find("--restart-mode model") == std::string::npos
                         ? " --output-policy create"
                         : " --output-policy append";
    return arguments;
}

int run(const std::string& arguments) { return std::system((quote(WV_RUNTIME_RUNNER)+" "+explicitRestartMode(arguments)).c_str()); }

int runWithPid(const std::string& arguments, const std::filesystem::path& pidFile) {
    return std::system((quote(WV_RUNTIME_RUNNER)+" "+explicitRestartMode(arguments)+" >/dev/null 2>&1 & echo $! > "+quote(pidFile)+"; wait $!").c_str());
}

std::vector<char> bytes(const std::filesystem::path& path) {
    std::ifstream input(path,std::ios::binary);
    return {std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
}

std::string text(const std::filesystem::path& path) {
    std::ifstream input(path,std::ios::binary);
    return {std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
}

void waitForText(const std::filesystem::path& path, const std::string& expected) {
    for (std::size_t attempt = 0; attempt < 1000; ++attempt) {
        if (std::filesystem::exists(path) && text(path).find(expected) != std::string::npos) return;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    throw std::runtime_error("timed out waiting for runner phase "+expected);
}

void waitForNonemptyFile(const std::filesystem::path& path) {
    for (std::size_t attempt = 0; attempt < 1000; ++attempt) {
        if (std::filesystem::exists(path) && !text(path).empty()) return;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    throw std::runtime_error("timed out waiting for "+path.string());
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

double relativeDifference(const WVCheckpoint& first, const WVCheckpoint& second) {
    double difference = 0.0;
    double scale = 0.0;
    const std::vector<wavevortex::WVComplex64>* firstValues[] = {&first.state.coefficients.Ap,&first.state.coefficients.Am,&first.state.coefficients.A0};
    const std::vector<wavevortex::WVComplex64>* secondValues[] = {&second.state.coefficients.Ap,&second.state.coefficients.Am,&second.state.coefficients.A0};
    for (std::size_t component = 0; component < 3; ++component) {
        require(firstValues[component]->size() == secondValues[component]->size(),"restart comparison shapes differ");
        for (std::size_t index = 0; index < firstValues[component]->size(); ++index) {
            const auto left = (*firstValues[component])[index];
            const auto right = (*secondValues[component])[index];
            difference = std::max(difference,std::hypot(left.real-right.real,left.imag-right.imag));
            scale = std::max(scale,std::max(std::hypot(left.real,left.imag),std::hypot(right.real,right.imag)));
        }
    }
    return difference/std::max(scale,std::numeric_limits<double>::min());
}

bool hasFiniteState(const WVCheckpoint& checkpoint) {
    const std::vector<wavevortex::WVComplex64>* values[] = {
        &checkpoint.state.coefficients.Ap,&checkpoint.state.coefficients.Am,&checkpoint.state.coefficients.A0};
    for (const auto* component : values) {
        for (const auto value : *component) {
            if (!std::isfinite(value.real) || !std::isfinite(value.imag)) return false;
        }
    }
    return true;
}

double jsonNumber(const std::string& json, const std::string& key) {
    const auto position = json.find("\""+key+"\":");
    require(position != std::string::npos,"JSON field is missing: "+key);
    const char* begin = json.c_str()+position+key.size()+3;
    char* end = nullptr;
    const double value = std::strtod(begin,&end);
    require(end != begin,"JSON field is not numeric: "+key);
    return value;
}

void requirePortableRunnerPolicy(const std::string& report,const char* route) {
    require(report.find("\"variableKernelPolicy\":{")!=std::string::npos,
        std::string(route)+" runner report omitted variable-kernel policy");
    require(report.find("\"matrixBackend\":\"scalar\"")!=std::string::npos &&
        report.find("\"transformKind\":\"constant-stratification\"")!=std::string::npos &&
        report.find("\"spectralSchedule\":\"established-interleaved\"")!=std::string::npos &&
        report.find("\"horizontalSchedule\":\"full-fft\"")!=std::string::npos &&
        report.find("\"effectiveFFTThreads\":1")!=std::string::npos &&
        report.find("\"noFallback\":true")!=std::string::npos,
        std::string(route)+" reference policy changed or was incompletely reported");
#if WV_TEST_COMPACT_VARIABLE_POLICY
    require(report.find("\"buildEnabled\":true")!=std::string::npos &&
        report.find("\"selection\":\"non-variable-transform\"")!=std::string::npos,
        std::string(route)+" did not report its inspected non-variable family");
#else
    require(report.find("\"buildEnabled\":false")!=std::string::npos &&
        report.find("\"selection\":\"build-disabled\"")!=std::string::npos,
        std::string(route)+" changed the default-off runner policy");
#endif
}

double createModelOutputFixture(const WVCheckpoint& checkpoint,const std::filesystem::path& path) {
    WVPortableObserverRecord record;
    const std::vector<std::size_t> shape={checkpoint.state.coefficients.shape.rows,
        checkpoint.state.coefficients.shape.columns};
    for (const auto* name:{"Ap","Am","A0"})
        record.stateBlocks.push_back({name,WVStateScalarType::complex64,shape,
            WVToleranceKind::coefficientEnergyScaled,1e-6,WVStateOwnership::integratorOwned,
            WVRestartRequirement::requiredDynamicState});
    WVObserverRecord coefficients;
    coefficients.identifier="coefficients";
    coefficients.name="wave-vortex coefficients";
    coefficients.typeIdentifier="WVCoefficients";
    coefficients.stateBlockIdentifiers={"Ap","Am","A0"};
    record.observers.push_back(std::move(coefficients));
    const double interval=1e-7,finalTime=checkpoint.state.t+2*interval;
    record.outputFiles={{"primary",path.string(),{{"restart","wave-vortex",
        {interval,checkpoint.state.t,finalTime},{"coefficients"},true}}}};
    WVPortableObserverDescriptor descriptor;
    auto status=WVPortableObserverDescriptor::create(record,test::extensionCatalog(),descriptor);
    require(bool(status),status.message);
    WVIntegrationStateLayout layout;
    status=WVIntegrationStateLayout::create(checkpoint.state.coefficients.shape,descriptor,layout);
    require(bool(status),status.message);
    WVOutputPlan plan;
    status=WVOutputPlan::create(descriptor,test::extensionCatalog(),checkpoint.state.t,finalTime,{},plan);
    require(bool(status) && plan.eventCount()>0,status.message);
    WVModelOutputNetCDFConfiguration configuration{test::extensionCatalog(),checkpoint,false};
    WVModelOutputNetCDFSink sink;
    auto persistence=WVModelOutputNetCDFSink::createNew(configuration,descriptor,plan,layout,nullptr,sink);
    require(bool(persistence),persistence.message);
    status=sink.preflight(plan); require(bool(status),status.message);
    for (std::size_t index=0;index<plan.eventCount();++index) {
        const auto planned=plan.event(index);
        require(planned.routeCount==1,"Model-output fixture route count differs");
        auto state=checkpoint.state.view(); state.t=planned.scheduledTime;
        WVOutputEvent event;
        event.eventOrdinal=planned.eventOrdinal;
        event.scheduledTime=planned.scheduledTime;
        event.state={state,nullptr,0};
        event.routes=planned.routes;
        event.routeCount=planned.routeCount;
        WVOutputDeliveryResult delivery;
        status=sink.deliver(event,planned.routes[0],delivery);
        require(bool(status),status.message);
    }
    persistence=sink.close(); require(bool(persistence),persistence.message);
    return plan.event(plan.eventCount()-1).scheduledTime;
}

} // namespace

int main() {
    try {
        const auto directory = std::filesystem::temp_directory_path()/"wave-vortex-run-contract";
        std::filesystem::remove_all(directory);
        std::filesystem::create_directories(directory);
        const auto input = std::filesystem::path(WV_RUNTIME_FIXTURE_DIR)/"forcing-mixed-hydrostatic.nc";
        {
          const auto callbackReport = directory / "callback-failure.json";
          std::vector<std::string> arguments{
              "wave-vortex-run",
              input.string(),
              (directory / "callback-failure.nc").string(),
              "--restart-mode",
              "coefficients",
              "--output-policy",
              "create",
              "--delta-t",
              "1e-5",
              "--steps",
              "2",
              "--fft-provider",
              "reference",
              "--report",
              callbackReport.string()};
          std::vector<char *> argv;
          for (auto &argument : arguments)
            argv.push_back(argument.data());
          const int code = runWaveVortex(
              static_cast<int>(argv.size()), argv.data(),
              test::extensionCatalog(), {[](const auto &) -> bool {
                throw std::runtime_error("injected callback error");
              }});
          const auto failure = text(callbackReport);
          require(code == 7 &&
                      failure.find("\"stage\":\"stop-callback\"") !=
                          std::string::npos &&
                      failure.find("\"reason\":\"callback-failure\"") !=
                          std::string::npos &&
                      jsonNumber(failure, "acceptedStepCount") == 0.0,
                  "reusable runner must report callback failure separately "
                  "from integration and output failure");
        }
        // Warmup and measured phases share one invocation-level stop result.
        for (const auto denseOutputs : {0, 1}) {
            for (const auto scenario : {"stop", "continue", "throw"}) {
                const auto stem = std::string("warmup-") + scenario + "-" +
                                  std::to_string(denseOutputs);
                const auto warmupReport = directory / (stem + ".json");
                const auto warmupOutput = directory / (stem + ".nc");
                std::vector<std::string> arguments{
                    "wave-vortex-run", input.string(), warmupOutput.string(),
                    "--restart-mode", "coefficients", "--output-policy", "create",
                    "--delta-t", "1e-5", "--steps", "2", "--fft-provider",
                    "reference", "--benchmark-warmup-steps", "1",
                    "--report", warmupReport.string()};
                if (denseOutputs != 0) {
                    arguments.push_back("--benchmark-dense-outputs-per-step");
                    arguments.push_back(std::to_string(denseOutputs));
                }
                std::vector<char *> argv;
                for (auto &argument : arguments) argv.push_back(argument.data());
                std::size_t callbackCount = 0;
                bool requested = false;
                double requestTime = 0.0;
                const int code = runWaveVortex(
                    static_cast<int>(argv.size()), argv.data(),
                    test::extensionCatalog(), {[&](const auto &progress) {
                        ++callbackCount;
                        if (progress.boundary != WVIntegrationBoundary::acceptedStep ||
                            requested || std::string(scenario) == "continue") return false;
                        requested = true;
                        requestTime = progress.acceptedTime;
                        if (std::string(scenario) == "throw")
                            throw std::runtime_error("injected warmup callback failure");
                        return true; // A one-shot request must survive phase transitions.
                    }});
                const auto reportText = text(warmupReport);
                require(jsonNumber(reportText, "callbackEvaluationCount") == callbackCount,
                        "runner must aggregate callback counts across warmup and measurement");
                if (std::string(scenario) == "throw") {
                    require(code == 7 &&
                                reportText.find("\"stage\":\"stop-callback\"") != std::string::npos &&
                                reportText.find("\"reason\":\"callback-failure\"") != std::string::npos &&
                                jsonNumber(reportText, "acceptedStepCount") == 1.0 &&
                                jsonNumber(reportText, "finalAcceptedTime") == requestTime &&
                                !std::filesystem::exists(warmupOutput),
                            "warmup callback failure must retain structured classification and state metrics");
                } else {
                    require(code == 0, "controlled benchmark warmup run failed");
                    WVCheckpoint stopped;
                    const auto status = WVCheckpointReader::read(
                        warmupOutput.string(), *test::extensionCatalog(), stopped);
                    require(static_cast<bool>(status), status.message);
                    if (std::string(scenario) == "stop") {
                        require(reportText.find("\"status\":\"stopped\"") != std::string::npos &&
                                    jsonNumber(reportText, "stepCount") == 1.0 &&
                                    jsonNumber(reportText, "requestedAtAcceptedTime") == requestTime &&
                                    stopped.state.t == requestTime,
                                "a warmup stop must preserve its request and skip all measured steps");
                    } else {
                        require(jsonNumber(reportText, "stepCount") == 3.0,
                                "continuing warmup must complete warmup and measured steps");
                    }
                }
            }
        }
        const auto output = directory/"output.nc";
        const auto report = directory/"report.json";
        require(run(quote(input)+" "+quote(output)+" --delta-t 1e-5 --steps 2 --fft-provider reference --report "+quote(report)) == 0,"runner step execution failed");
        WVCheckpoint checkpoint;
        auto status = WVCheckpointReader::read(output.string(), *test::extensionCatalog(),checkpoint);
        require(static_cast<bool>(status),status.message);
        require(checkpoint.state.t > 8.949 && checkpoint.state.t < 9.025 &&
                    hasFiniteState(checkpoint),
                "runner did not produce a finite advanced state");
        require(std::filesystem::file_size(report) > 0,"runner did not write its report");
        const auto reportText = text(report);
        requirePortableRunnerPolicy(reportText,"checkpoint route");
        const auto checkpointStateBytes = jsonNumber(reportText,"checkpointState");
        const auto modelStateBytes = jsonNumber(reportText,"modelState");
        const auto extensionCatalogBytes = jsonNumber(reportText,"extensionCatalog");
        const auto integrationSystemBytes = jsonNumber(reportText,"integrationSystem");
        const auto integratorBytes = jsonNumber(reportText,"integratorPersistent");
        const auto knownPersistentBytes = jsonNumber(reportText,"knownPersistent");
        const auto fullModelPersistentBytes = jsonNumber(reportText,"fullModelPersistent");
        require(reportText.find(
                    "\"scope\":\"application-and-provider-reported-cpp-storage\"") !=
                    std::string::npos &&
                    reportText.find("opaque-fftw-plan-internals") !=
                        std::string::npos &&
                    reportText.find("opaque-netcdf-library-internals") !=
                        std::string::npos,
                "runner omitted the exact storage-accounting scope");
        require(extensionCatalogBytes > 0.0,
                "runner omitted the retained extension catalog");
        require(knownPersistentBytes ==
                    jsonNumber(reportText,"modelFacade")+
                        extensionCatalogBytes+checkpointStateBytes+
                        integrationSystemBytes+integratorBytes+
                        jsonNumber(reportText,"modelOutputConfiguration")+
                        jsonNumber(reportText,"scheduledOutput"),
                "runner known-persistent ownership formula is not exact");
        require(fullModelPersistentBytes ==
                    knownPersistentBytes+
                        std::max(0.0,modelStateBytes-checkpointStateBytes)+
                        jsonNumber(reportText,"modelOutputEvaluation")+
                        jsonNumber(reportText,"modelOutputSink"),
                "runner full persistent accounting is not exact");
        const auto occurrenceRetained =
            jsonNumber(reportText,"occurrenceWorkspaceRetained");
        const auto occurrenceMaximumLive =
            jsonNumber(reportText,"occurrenceWorkspaceMaximumLive");
        const auto orchestrationMaximumLive =
            jsonNumber(reportText,"outputOrchestrationMaximumLive");
        require(jsonNumber(reportText,"outputDriverRetained") == 0.0 &&
                    jsonNumber(reportText,"knownRetained") ==
                        knownPersistentBytes &&
                    jsonNumber(reportText,"knownMaximumLive") ==
                        knownPersistentBytes + orchestrationMaximumLive &&
                    jsonNumber(reportText,"fullModelRetained") ==
                        fullModelPersistentBytes &&
                    jsonNumber(reportText,"fullModelMaximumLive") ==
                        fullModelPersistentBytes - occurrenceRetained +
                            std::max(occurrenceRetained,
                                     occurrenceMaximumLive) +
                            orchestrationMaximumLive,
                "runner full retained/liveness accounting is inconsistent");
        require(reportText.find("\"arrayTraffic\"") != std::string::npos && reportText.find("\"stageStateConstructionReads\":1296") != std::string::npos,"runner omitted exact RK4 traffic diagnostics");
        require(reportText.find("\"stageFluxClearWrites\":0") != std::string::npos && reportText.find("\"weightedFluxInitializationReads\":216") != std::string::npos,"runner omitted eliminated-clear and first-stage initialization diagnostics");
        require(reportText.find("\"contractAbstractionAdditionalArrayStorage\":0") != std::string::npos,"runner reported array-sized contract workspace");
        require(reportText.find("\"integrationBreakdownSeconds\"") != std::string::npos && reportText.find("\"waveVortexFlux\"") != std::string::npos && reportText.find("\"tracerAdvection\"") != std::string::npos && reportText.find("\"tracerForward\"") != std::string::npos && reportText.find("\"tracerAntialias\"") != std::string::npos && reportText.find("\"observerEvaluation\"") != std::string::npos,"runner omitted the integration cost decomposition");
        const auto modelInput=directory/"model-route.nc";
        const auto modelReport=directory/"model-route-report.json";
        const auto modelTime=createModelOutputFixture(checkpoint,modelInput);
        require(run(quote(modelInput)+" --restart-mode model --output-policy append --delta-t 1e-7 --final-time "+number(modelTime+1e-7)+" --fft-provider reference --report "+quote(modelReport))==0,
            "runner model-output route failed");
        requirePortableRunnerPolicy(text(modelReport),"model-output route");
        const auto protectedOutput = directory/"protected-output.nc";
        {
            std::ofstream stream(protectedOutput,std::ios::binary);
            stream << "protected destination";
        }
        const auto protectedBytes = bytes(protectedOutput);
        require(run(quote(input)+" "+quote(protectedOutput)+" --restart-mode coefficients --output-policy create --delta-t 1e-7 --steps 1 --fft-provider reference >/dev/null 2>&1") != 0,
                "create policy replaced an existing destination");
        require(bytes(protectedOutput) == protectedBytes,
                "failed create changed the existing destination");
        require(run(quote(input)+" "+quote(protectedOutput)+" --restart-mode coefficients --output-policy replace --delta-t 1e-7 --steps 1 --fft-provider reference >/dev/null 2>&1") == 0,
                "explicit replacement did not succeed");
        WVCheckpoint replacedCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(
                    protectedOutput.string(), *test::extensionCatalog(),
                    replacedCheckpoint)),
                "explicit replacement did not produce a checkpoint");
        const auto inputBytes = bytes(input);
        require(run(quote(input)+" "+quote(input)+" --restart-mode coefficients --output-policy replace --delta-t 1e-7 --steps 1 --fft-provider reference >/dev/null 2>&1") != 0,
                "coefficient replacement accepted an input/output alias");
        require(bytes(input) == inputBytes,
                "input/output alias failure changed the source");
        WVCheckpoint initialCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(input.string(), *test::extensionCatalog(),initialCheckpoint)),"scheduled-output input is unreadable");
        const double scheduledMidpoint = initialCheckpoint.state.t+5e-6;
        const double scheduledEndpoint = initialCheckpoint.state.t+1e-5;
        const double scheduledFinalTime = initialCheckpoint.state.t+2e-5;
        const auto scheduledDirectory = directory/"scheduled-fixed";
        const auto scheduledReport = directory/"scheduled-fixed-report.json";
        const auto fixedScheduledArguments = std::string(" --delta-t 1e-5 --final-time ")+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(initialCheckpoint.state.t)+" --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledEndpoint)+" --output-time "+number(scheduledFinalTime)+" --output-directory "+quote(scheduledDirectory)+" --report "+quote(scheduledReport);
        require(run(quote(input)+fixedScheduledArguments) == 0,"fixed runner scheduled-output execution failed");
        const std::vector<std::filesystem::path> scheduledPaths{
            scheduledDirectory/"checkpoint-000001.nc",scheduledDirectory/"checkpoint-000002.nc",
            scheduledDirectory/"checkpoint-000003.nc",scheduledDirectory/"checkpoint-000004.nc"};
        for (std::size_t index = 0; index < scheduledPaths.size(); ++index) {
            WVCheckpoint scheduledCheckpoint;
            require(static_cast<bool>(WVCheckpointReader::read(scheduledPaths[index].string(), *test::extensionCatalog(),scheduledCheckpoint)),"fixed scheduled checkpoint is not readable");
            const double expectedTimes[] = {initialCheckpoint.state.t,scheduledMidpoint,scheduledEndpoint,scheduledFinalTime};
            require(std::abs(scheduledCheckpoint.state.t-expectedTimes[index]) <= 1e-14,"fixed scheduled checkpoint has the wrong time");
        }
        const auto scheduledReportText = text(scheduledReport);
        require(scheduledReportText.find("\"requestedCount\":4") != std::string::npos && scheduledReportText.find("\"committedCount\":4") != std::string::npos,"fixed scheduled-output report omitted counts");
        require(scheduledReportText.find("\"eventKind\":\"initial\"") != std::string::npos && scheduledReportText.find("\"eventKind\":\"interpolated\"") != std::string::npos && scheduledReportText.find("\"eventKind\":\"accepted-endpoint\"") != std::string::npos,"fixed scheduled-output report omitted event kinds");
        const auto scheduledDriverMaximumLive =
            jsonNumber(scheduledReportText,"outputDriverMaximumLive");
        const auto scheduledPlanMaximumLive =
            jsonNumber(scheduledReportText,"outputPlanMaximumLive");
        const auto scheduledOrchestrationMaximumLive =
            jsonNumber(scheduledReportText,
                       "outputOrchestrationMaximumLive");
        require(jsonNumber(scheduledReportText,"outputDriverRetained") == 0.0 &&
                    scheduledDriverMaximumLive > 0.0 &&
                    scheduledPlanMaximumLive > 0.0 &&
                    scheduledOrchestrationMaximumLive ==
                        scheduledDriverMaximumLive +
                            scheduledPlanMaximumLive &&
                    jsonNumber(scheduledReportText,"knownMaximumLive") ==
                        jsonNumber(scheduledReportText,"knownRetained") +
                            scheduledOrchestrationMaximumLive,
                "scheduled-output driver/plan liveness is not exact");
        const auto scheduledSentinel = bytes(scheduledPaths.front());
        require(run(quote(input)+fixedScheduledArguments+" >/dev/null 2>&1") != 0,"scheduled output replaced existing checkpoints");
        require(bytes(scheduledPaths.front()) == scheduledSentinel,"scheduled output collision changed an existing checkpoint");

        const auto lateFailureDirectory = directory/"late-failure";
        const auto lateFailurePhase = directory/"late-failure.phase";
        const auto lateFailureReport = directory/"late-failure.json";
        const auto lateFailureArguments = quote(input)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledFinalTime)+" --output-directory "+quote(lateFailureDirectory)+" --phase-file "+quote(lateFailurePhase)+" --report "+quote(lateFailureReport)+" >/dev/null 2>&1";
        int lateFailureStatus = 0;
        std::thread lateFailureWorker([&] { lateFailureStatus = run(lateFailureArguments); });
        waitForText(lateFailurePhase,"output-committed:1");
        const auto racedDestination = lateFailureDirectory/"checkpoint-000002.nc";
        {
            std::ofstream stream(racedDestination,std::ios::binary);
            stream << "protected";
        }
        lateFailureWorker.join();
        require(lateFailureStatus != 0,"later scheduled-output collision unexpectedly succeeded");
        WVCheckpoint earlierCommitted;
        require(static_cast<bool>(WVCheckpointReader::read((lateFailureDirectory/"checkpoint-000001.nc").string(), *test::extensionCatalog(),earlierCommitted)),"later scheduled-output failure corrupted the earlier checkpoint");
        require(text(racedDestination) == "protected","later scheduled-output failure replaced the raced destination");
        const auto lateFailureReportText = text(lateFailureReport);
        require(lateFailureReportText.find("\"committedCount\":1") != std::string::npos && lateFailureReportText.find("\"status\":\"failed\"") != std::string::npos,"later scheduled-output failure report omitted partial results");

#if defined(__APPLE__) || defined(__linux__)
        for (const int scenario : {0, 1, 2}) {
          const auto stem = "graceful-" + std::to_string(scenario);
          const auto stopDirectory = directory / stem;
          const auto stopPhase = directory / (stem + ".phase");
          const auto stopPid = directory / (stem + ".pid");
          const auto stopReport = directory / (stem + ".json");
          const auto stopArguments =
              quote(input) + " --delta-t 1e-5 --final-time " +
              number(initialCheckpoint.state.t + 1e-3) +
              " --fft-provider reference --output-time " +
              number(initialCheckpoint.state.t) + " --output-time " +
              number(scheduledMidpoint) + " --output-time " +
              number(initialCheckpoint.state.t + 1e-3) +
              " --output-directory " + quote(stopDirectory) + " --phase-file " +
              quote(stopPhase) + " --report " + quote(stopReport);
          int stopStatus = -1;
          std::thread stopWorker(
              [&] { stopStatus = runWithPid(stopArguments, stopPid); });
          waitForText(stopPhase,
                      scenario == 1 ? "output-committed:2" : "steady-retained");
          waitForNonemptyFile(stopPid);
          std::ifstream processInput(stopPid);
          int process = 0;
          processInput >> process;
          require(process > 0 && ::kill(process, SIGINT) == 0,
                  "unable to request graceful CLI stop");
          int repeatedSignalResult = 0;
          int repeatedSignalError = 0;
          if (scenario == 2) {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            repeatedSignalResult = ::kill(process, SIGINT);
            repeatedSignalError = errno;
            require(repeatedSignalResult == 0 || repeatedSignalError == ESRCH,
                    "unable to issue repeated interrupt");
          }
          stopWorker.join();
          if (scenario == 2 && stopStatus != 0) {
            // A successful kill can target a zombie whose graceful exit status
            // is already fixed, so only an actual SIGINT status proves the
            // repeated signal won the race.
            const bool directSignal = WIFSIGNALED(stopStatus) &&
                WTERMSIG(stopStatus) == SIGINT;
            const bool shellSignal = WIFEXITED(stopStatus) &&
                WEXITSTATUS(stopStatus) == 128 + SIGINT;
            require(repeatedSignalResult == 0 &&
                        (directSignal || shellSignal),
                    "repeated SIGINT produced an unrelated process failure");
            continue;
          }
          require(stopStatus == 0, "graceful SIGINT must return success");
          const auto stoppedReport = text(stopReport);
          require(stoppedReport.find("\"status\":\"stopped\"") !=
                          std::string::npos &&
                      stoppedReport.find("\"reason\":\"stop-requested\"") !=
                          std::string::npos &&
                      jsonNumber(stoppedReport, "finalAcceptedTime") <
                          initialCheckpoint.state.t + 1e-3,
                  "graceful stop report must distinguish actual accepted time "
                  "from target");
          WVCheckpoint stoppedCheckpoint;
          require(bool(WVCheckpointReader::read(
                      (stopDirectory / "checkpoint-000001.nc").string(),
                      *test::extensionCatalog(), stoppedCheckpoint)),
                  "graceful stop corrupted committed checkpoint");
          require(
              !std::filesystem::exists(stopDirectory / "checkpoint-000003.nc"),
              "graceful stop wrote an occurrence after its accepted boundary");
        }

        const auto interruptedDirectory = directory/"interrupted";
        const auto interruptedPhase = directory/"interrupted.phase";
        const auto interruptedPid = directory/"interrupted.pid";
        const auto interruptedArguments = quote(input)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledFinalTime)+" --output-directory "+quote(interruptedDirectory)+" --phase-file "+quote(interruptedPhase);
        int interruptedStatus = 0;
        std::thread interruptedWorker([&] { interruptedStatus = runWithPid(interruptedArguments,interruptedPid); });
        waitForText(interruptedPhase,"output-committed:1");
        waitForNonemptyFile(interruptedPid);
        std::ifstream pidInput(interruptedPid);
        int processIdentifier = 0;
        pidInput >> processIdentifier;
        require(processIdentifier > 0 && ::kill(processIdentifier,SIGTERM) == 0,"unable to interrupt scheduled-output runner");
        interruptedWorker.join();
        require(interruptedStatus != 0,"interrupted scheduled-output runner reported success");
        WVCheckpoint interruptedCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read((interruptedDirectory/"checkpoint-000001.nc").string(), *test::extensionCatalog(),interruptedCheckpoint)),"interruption corrupted the committed checkpoint");
        require(!std::filesystem::exists(interruptedDirectory/"checkpoint-000002.nc"),"interruption created a later checkpoint");
#endif

        const auto fixedControl = directory/"fixed-scheduled-control.nc";
        require(run(quote(input)+" "+quote(fixedControl)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference >/dev/null 2>&1") == 0,"fixed scheduled-output control failed");
        WVCheckpoint fixedControlCheckpoint;
        WVCheckpoint fixedScheduledFinal;
        require(static_cast<bool>(WVCheckpointReader::read(fixedControl.string(), *test::extensionCatalog(),fixedControlCheckpoint)) && static_cast<bool>(WVCheckpointReader::read(scheduledPaths.back().string(), *test::extensionCatalog(),fixedScheduledFinal)),"fixed scheduled-output comparison is unreadable");
        require(relativeDifference(fixedControlCheckpoint,fixedScheduledFinal) == 0.0,"fixed scheduled output changed the accepted trajectory");
        const auto fixedInterpolatedRestart = directory/"fixed-interpolated-restart.nc";
        require(run(quote(scheduledPaths[1])+" "+quote(fixedInterpolatedRestart)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference >/dev/null 2>&1") == 0,"fixed interpolated checkpoint did not restart");
        WVCheckpoint fixedInterpolatedRestartCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(fixedInterpolatedRestart.string(), *test::extensionCatalog(),fixedInterpolatedRestartCheckpoint)),"fixed interpolated restart is unreadable");
        require(std::abs(fixedInterpolatedRestartCheckpoint.state.t-scheduledFinalTime) <= 1e-12 && hasFiniteState(fixedInterpolatedRestartCheckpoint),"fixed interpolated checkpoint did not produce a finite restart at the requested final time");
        const auto adaptiveOutput = directory/"adaptive-output.nc";
        const auto adaptiveReport = directory/"adaptive-report.json";
        require(run(quote(input)+" "+quote(adaptiveOutput)+" --delta-t 0.037 --initial-step 0.02 --maximum-step 0.01 --steps 2 --integrator adaptive-rk23 --relative-tolerance 1e-3 --absolute-tolerance 1e-6 --fft-provider reference --report "+quote(adaptiveReport)) == 0,"adaptive runner execution failed");
        WVCheckpoint adaptiveCheckpoint;
        status = WVCheckpointReader::read(adaptiveOutput.string(), *test::extensionCatalog(),adaptiveCheckpoint);
        require(static_cast<bool>(status) && adaptiveCheckpoint.state.t > checkpoint.state.t-0.074,"adaptive runner output is not readable or did not advance");
        const auto adaptiveReportText = text(adaptiveReport);
        require(adaptiveReportText.find("\"id\":\"adaptive-rk23\"") != std::string::npos && adaptiveReportText.find("\"controller\":\"matlab-ode23-v1\"") != std::string::npos && adaptiveReportText.find("\"requestedInitialStep\":0.02") != std::string::npos && adaptiveReportText.find("\"effectiveInitialStep\":0.01") != std::string::npos && adaptiveReportText.find("\"effectiveMaximumStep\":0.01") != std::string::npos && adaptiveReportText.find("\"toleranceHash\":") != std::string::npos && adaptiveReportText.find("\"toleranceHashClearedMantissaBits\":20") != std::string::npos && adaptiveReportText.find("\"acceptedSteps\":[") != std::string::npos && adaptiveReportText.find("\"rejectedStepCount\":") != std::string::npos && adaptiveReportText.find("\"nextStepSize\":") != std::string::npos && adaptiveReportText.find("\"denseOutputEvaluationCount\":") != std::string::npos && adaptiveReportText.find("\"denseOutputSeconds\":") != std::string::npos,"adaptive runner report omitted method diagnostics");
        const auto rk45Output = directory/"adaptive-rk45-output.nc";
        const auto rk45Report = directory/"adaptive-rk45-report.json";
        require(run(quote(input)+" "+quote(rk45Output)+" --delta-t 0.037 --initial-step 0.02 --maximum-step 0.01 --steps 2 --integrator adaptive-rk45 --relative-tolerance 1e-3 --absolute-tolerance 1e-6 --fft-provider reference --report "+quote(rk45Report)) == 0,"adaptive RK45 runner execution failed");
        WVCheckpoint rk45Checkpoint;
        status = WVCheckpointReader::read(rk45Output.string(), *test::extensionCatalog(),rk45Checkpoint);
        require(static_cast<bool>(status) && rk45Checkpoint.state.t > checkpoint.state.t-0.074,"adaptive RK45 runner output is not readable or did not advance");
        const auto rk45ReportText = text(rk45Report);
        require(rk45ReportText.find("\"id\":\"adaptive-rk45\"") != std::string::npos && rk45ReportText.find("\"controller\":\"matlab-ode45-v1\"") != std::string::npos && rk45ReportText.find("\"workspaceStateEquivalentCount\":7") != std::string::npos && rk45ReportText.find("\"denseHistoryStateEquivalentCount\":0") != std::string::npos && rk45ReportText.find("\"buffer\":\"k2/k7\"") != std::string::npos && rk45ReportText.find("\"producer\":\"second-stage then endpoint right-hand side\"") != std::string::npos && rk45ReportText.find("\"integratorStorageLedger\":") != std::string::npos && rk45ReportText.find("\"byteLedgerAgreement\":true") != std::string::npos && rk45ReportText.find("\"sharedAbstractionStateSizedCopyCount\":0") != std::string::npos && rk45ReportText.find("\"diagnosticBytes\":") != std::string::npos,"adaptive RK45 runner report omitted method or exact memory diagnostics");
        const auto rk78Output = directory/"adaptive-rk78-output.nc";
        const auto rk78Report = directory/"adaptive-rk78-report.json";
        require(run(quote(input)+" "+quote(rk78Output)+" --delta-t 0.037 --initial-step 1e-5 --maximum-step 1e-5 --steps 2 --integrator adaptive-rk78 --relative-tolerance 1e-3 --absolute-tolerance 1e-6 --fft-provider reference --report "+quote(rk78Report)) == 0,"adaptive RK78 runner execution failed");
        WVCheckpoint rk78Checkpoint;
        status = WVCheckpointReader::read(rk78Output.string(), *test::extensionCatalog(),rk78Checkpoint);
        require(static_cast<bool>(status) && rk78Checkpoint.state.t > checkpoint.state.t-0.074,"adaptive RK78 runner output is not readable or did not advance");
        const auto rk78ReportText = text(rk78Report);
        require(rk78ReportText.find("\"id\":\"adaptive-rk78\"") != std::string::npos && rk78ReportText.find("\"controller\":\"matlab-ode78-v1\"") != std::string::npos && rk78ReportText.find("\"workspaceStateEquivalentCount\":11") != std::string::npos && rk78ReportText.find("\"workspaceMaximumLiveStateEquivalentCount\":11") != std::string::npos && rk78ReportText.find("\"denseHistoryStateEquivalentCount\":0") != std::string::npos && rk78ReportText.find("\"retainedBaseStageStateEquivalentCount\":0") != std::string::npos && rk78ReportText.find("\"continuousExtensionRightHandSideEvaluationCount\":0") != std::string::npos && rk78ReportText.find("\"continuousExtensionWorkspaceMaximumLiveBytes\":0") != std::string::npos && rk78ReportText.find("\"continuousExtensionBytes\":0") != std::string::npos && rk78ReportText.find("\"buffer\":\"k2/k3/k5\"") != std::string::npos && rk78ReportText.find("\"errorPolicyBytes\":") != std::string::npos,"adaptive RK78 runner report omitted method or exact endpoint-only memory diagnostics");
        const auto adaptiveSeriesDirectory = directory/"scheduled-adaptive";
        const auto adaptiveSeriesReport = directory/"scheduled-adaptive-report.json";
        const auto adaptiveSeriesArguments = std::string(" --delta-t 1e-5 --final-time ")+number(scheduledFinalTime)+" --integrator adaptive-rk23 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledFinalTime)+" --output-directory "+quote(adaptiveSeriesDirectory)+" --output-pattern 'state-{index}-{time}.nc' --report "+quote(adaptiveSeriesReport);
        require(run(quote(input)+adaptiveSeriesArguments) == 0,"adaptive runner public scheduled-output execution failed");
        const auto adaptiveSeriesText = text(adaptiveSeriesReport);
        require(adaptiveSeriesText.find("\"requestedCount\":2") != std::string::npos && adaptiveSeriesText.find("\"committedCount\":2") != std::string::npos && adaptiveSeriesText.find("\"status\":\"committed\"") != std::string::npos,"adaptive scheduled-output report is incomplete");
        std::vector<std::filesystem::path> adaptiveSeriesPaths;
        for (const auto& entry : std::filesystem::directory_iterator(adaptiveSeriesDirectory)) adaptiveSeriesPaths.push_back(entry.path());
        std::sort(adaptiveSeriesPaths.begin(),adaptiveSeriesPaths.end());
        require(adaptiveSeriesPaths.size() == 2,"adaptive scheduled output wrote the wrong number of checkpoints");
        for (const auto& path : adaptiveSeriesPaths) {
            WVCheckpoint emitted;
            require(static_cast<bool>(WVCheckpointReader::read(path.string(), *test::extensionCatalog(),emitted)),"adaptive scheduled checkpoint is not readable");
        }
        const auto adaptiveSeriesControl = directory/"scheduled-adaptive-control.nc";
        const auto adaptiveSeriesControlReport = directory/"scheduled-adaptive-control.json";
        require(run(quote(input)+" "+quote(adaptiveSeriesControl)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --integrator adaptive-rk23 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference --report "+quote(adaptiveSeriesControlReport)) == 0,"adaptive scheduled-output control failed");
        WVCheckpoint adaptiveSeriesFinal;
        WVCheckpoint adaptiveSeriesControlCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(adaptiveSeriesPaths.back().string(), *test::extensionCatalog(),adaptiveSeriesFinal)) && static_cast<bool>(WVCheckpointReader::read(adaptiveSeriesControl.string(), *test::extensionCatalog(),adaptiveSeriesControlCheckpoint)),"adaptive scheduled-output comparison is unreadable");
        require(relativeDifference(adaptiveSeriesFinal,adaptiveSeriesControlCheckpoint) == 0.0,"adaptive scheduled output changed the accepted trajectory");
        const auto adaptiveSeriesControlText = text(adaptiveSeriesControlReport);
        for (const auto& field : {"stepCount","rejectedStepCount","rhsEvaluationCount"}) require(jsonNumber(adaptiveSeriesText,field) == jsonNumber(adaptiveSeriesControlText,field),"adaptive scheduled output changed "+std::string(field));
        const auto adaptiveInterpolatedRestart = directory/"adaptive-interpolated-restart.nc";
        require(run(quote(adaptiveSeriesPaths.front())+" "+quote(adaptiveInterpolatedRestart)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --integrator adaptive-rk23 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference >/dev/null 2>&1") == 0,"adaptive interpolated checkpoint did not restart");
        WVCheckpoint adaptiveInterpolatedRestartCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(adaptiveInterpolatedRestart.string(), *test::extensionCatalog(),adaptiveInterpolatedRestartCheckpoint)),"adaptive interpolated restart is unreadable");
        require(relativeDifference(adaptiveSeriesControlCheckpoint,adaptiveInterpolatedRestartCheckpoint) <= 1e-6,"adaptive interpolated restart exceeded tolerance-equivalent continuation");
        const auto rk78SeriesDirectory = directory/"scheduled-rk78";
        const auto rk78SeriesReport = directory/"scheduled-rk78-report.json";
        const auto rk78SeriesArguments = std::string(" --delta-t 2e-5 --initial-step 2e-5 --maximum-step 2e-5 --final-time ")+number(scheduledFinalTime)+" --integrator adaptive-rk78 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledFinalTime)+" --output-directory "+quote(rk78SeriesDirectory)+" --output-pattern 'state-{index}-{time}.nc' --report "+quote(rk78SeriesReport);
        require(run(quote(input)+rk78SeriesArguments) == 0,"adaptive RK78 public scheduled-output execution failed");
        std::vector<std::filesystem::path> rk78SeriesPaths;
        for (const auto& entry : std::filesystem::directory_iterator(rk78SeriesDirectory)) rk78SeriesPaths.push_back(entry.path());
        std::sort(rk78SeriesPaths.begin(),rk78SeriesPaths.end());
        require(rk78SeriesPaths.size() == 2,"adaptive RK78 scheduled output wrote the wrong number of checkpoints");
        const auto rk78SeriesControl = directory/"scheduled-rk78-control.nc";
        const auto rk78SeriesControlReport = directory/"scheduled-rk78-control.json";
        const auto rk78ControlArguments = std::string(" --delta-t 2e-5 --initial-step 2e-5 --maximum-step 2e-5 --final-time ")+number(scheduledFinalTime)+" --integrator adaptive-rk78 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference --report "+quote(rk78SeriesControlReport);
        require(run(quote(input)+" "+quote(rk78SeriesControl)+rk78ControlArguments) == 0,"adaptive RK78 scheduled-output control failed");
        WVCheckpoint rk78SeriesFinal;
        WVCheckpoint rk78SeriesControlCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(rk78SeriesPaths.back().string(), *test::extensionCatalog(),rk78SeriesFinal)) && static_cast<bool>(WVCheckpointReader::read(rk78SeriesControl.string(), *test::extensionCatalog(),rk78SeriesControlCheckpoint)),"adaptive RK78 scheduled-output comparison is unreadable");
        require(relativeDifference(rk78SeriesFinal,rk78SeriesControlCheckpoint) == 0.0,"adaptive RK78 scheduled output changed the endpoint trajectory");
        const auto rk78SeriesText = text(rk78SeriesReport);
        const auto rk78SeriesControlText = text(rk78SeriesControlReport);
        for (const auto& field : {"stepCount","rejectedStepCount","baseRightHandSideEvaluationCount"}) require(jsonNumber(rk78SeriesText,field) == jsonNumber(rk78SeriesControlText,field),"adaptive RK78 scheduled output changed "+std::string(field));
        require(jsonNumber(rk78SeriesText,"continuousExtensionRightHandSideEvaluationCount") == 4.0 && jsonNumber(rk78SeriesControlText,"continuousExtensionRightHandSideEvaluationCount") == 0.0,"adaptive RK78 scheduled output did not isolate extension RHS work");
        const auto rk78InterpolatedRestart = directory/"adaptive-rk78-interpolated-restart.nc";
        require(run(quote(rk78SeriesPaths.front())+" "+quote(rk78InterpolatedRestart)+" --delta-t 1e-5 --initial-step 1e-5 --maximum-step 1e-5 --final-time "+number(scheduledFinalTime)+" --integrator adaptive-rk78 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference >/dev/null 2>&1") == 0,"adaptive RK78 interpolated checkpoint did not restart");
        WVCheckpoint rk78InterpolatedRestartCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(rk78InterpolatedRestart.string(), *test::extensionCatalog(),rk78InterpolatedRestartCheckpoint)) && std::abs(rk78InterpolatedRestartCheckpoint.state.t-scheduledFinalTime) <= 1e-12 && hasFiniteState(rk78InterpolatedRestartCheckpoint),"adaptive RK78 interpolated checkpoint restart is invalid");
        const auto adaptiveDenseOutput = directory/"adaptive-dense-output.nc";
        const auto adaptiveDenseReport = directory/"adaptive-dense-report.json";
        require(run(quote(input)+" "+quote(adaptiveDenseOutput)+" --delta-t 1e-7 --steps 2 --integrator adaptive-rk23 --fft-provider reference --benchmark-dense-outputs-per-step 1 --report "+quote(adaptiveDenseReport)) == 0,"adaptive runner dense-output execution failed");
        const auto adaptiveDenseReportText = text(adaptiveDenseReport);
        require(adaptiveDenseReportText.find("\"denseOutputEvaluationCount\":2") != std::string::npos && adaptiveDenseReportText.find("\"interpolatedOutputCount\":2") != std::string::npos,"adaptive runner did not report its dense-output work");
        const auto rk45DenseOutput = directory/"adaptive-rk45-dense-output.nc";
        const auto rk45DenseReport = directory/"adaptive-rk45-dense-report.json";
        require(run(quote(input)+" "+quote(rk45DenseOutput)+" --delta-t 1e-7 --steps 2 --integrator adaptive-rk45 --fft-provider reference --benchmark-dense-outputs-per-step 1 --report "+quote(rk45DenseReport)) == 0,"adaptive RK45 runner dense-output execution failed");
        const auto rk45DenseReportText = text(rk45DenseReport);
        require(rk45DenseReportText.find("\"denseOutputEvaluationCount\":2") != std::string::npos && rk45DenseReportText.find("\"interpolatedOutputCount\":2") != std::string::npos && rk45DenseReportText.find("\"denseHistoryStateEquivalentCount\":6") != std::string::npos && rk45DenseReportText.find("\"denseHistoryRetainedWithinWorkspace\":") != std::string::npos && rk45DenseReportText.find("\"acceptedStepAdditionalArrayStorage\":0") != std::string::npos,"adaptive RK45 runner did not report its dense-output work or history");
        const auto rk78DenseOutput = directory/"adaptive-rk78-dense-output.nc";
        const auto rk78DenseReport = directory/"adaptive-rk78-dense-report.json";
        require(run(quote(input)+" "+quote(rk78DenseOutput)+" --delta-t 1e-7 --steps 2 --integrator adaptive-rk78 --fft-provider reference --benchmark-dense-outputs-per-step 1 --report "+quote(rk78DenseReport)) == 0,"adaptive RK78 runner dense-output execution failed");
        const auto rk78DenseReportText = text(rk78DenseReport);
        require(rk78DenseReportText.find("\"denseOutputEvaluationCount\":2") != std::string::npos && rk78DenseReportText.find("\"interpolatedOutputCount\":2") != std::string::npos && rk78DenseReportText.find("\"continuousExtensionRightHandSideEvaluationCount\":8") != std::string::npos && rk78DenseReportText.find("\"denseOutputCacheBuildCount\":2") != std::string::npos && rk78DenseReportText.find("\"denseOutputCacheReuseCount\":0") != std::string::npos && rk78DenseReportText.find("\"retainedBaseStageStateEquivalentCount\":8") != std::string::npos && rk78DenseReportText.find("\"continuousExtensionWorkspaceStateEquivalentCount\":4") != std::string::npos && rk78DenseReportText.find("\"continuousExtensionWorkspaceMaximumLiveStateEquivalentCount\":4") != std::string::npos && rk78DenseReportText.find("\"workspaceMaximumLiveStateEquivalentCount\":15") != std::string::npos && jsonNumber(rk78DenseReportText,"continuousExtensionBytes") > 0.0,"adaptive RK78 runner did not report lazy dense-output work or storage");
        const auto adaptiveScheduledOutput = directory/"adaptive-scheduled-output.nc";
        const auto adaptiveScheduledReport = directory/"adaptive-scheduled-report.json";
        require(run(quote(input)+" "+quote(adaptiveScheduledOutput)+" --delta-t 1e-7 --final-time 8.9510002 --integrator adaptive-rk23 --fft-provider reference --benchmark-output-count 4 --report "+quote(adaptiveScheduledReport)) == 0,"adaptive runner scheduled-output execution failed");
        const auto adaptiveScheduledReportText = text(adaptiveScheduledReport);
        require(adaptiveScheduledReportText.find("\"requestedOutputCount\":4") != std::string::npos && adaptiveScheduledReportText.find("\"interpolatedOutputCount\":4") != std::string::npos && adaptiveScheduledReportText.find("\"interpolationSeconds\":") != std::string::npos,"adaptive runner did not report scheduled-output work");
        const std::vector<std::string> forcingFixtures{
            "forcing-nonlinear.nc","forcing-adaptive-damping.nc","forcing-fixed-amplitude.nc","forcing-quadratic-bottom-friction.nc",
            "forcing-pseudo-topographic.nc","forcing-beta-plane.nc","forcing-mixed-hydrostatic.nc","forcing-mixed-nonhydrostatic.nc"};
        for (const auto& fixture : forcingFixtures) {
            const auto fixtureInput = std::filesystem::path(WV_RUNTIME_FIXTURE_DIR)/fixture;
            const auto fixtureOutput = directory/("adaptive-"+fixture);
            require(run(quote(fixtureInput)+" "+quote(fixtureOutput)+" --delta-t 1e-5 --steps 1 --integrator adaptive-rk23 --fft-provider reference >/dev/null 2>&1") == 0,"adaptive runner failed forcing fixture "+fixture);
        }
        const auto rk78FixedAmplitudeOutput = directory/"adaptive-rk78-forcing-fixed-amplitude.nc";
        require(run(quote(std::filesystem::path(WV_RUNTIME_FIXTURE_DIR)/"forcing-fixed-amplitude.nc")+" "+quote(rk78FixedAmplitudeOutput)+" --delta-t 1e-5 --steps 1 --integrator adaptive-rk78 --fft-provider reference --benchmark-dense-outputs-per-step 4 >/dev/null 2>&1") == 0,"adaptive RK78 dense output failed the fixed-amplitude constraint fixture");
        require(static_cast<bool>(WVCheckpointReader::read(input.string(), *test::extensionCatalog(),initialCheckpoint)),"adaptive restart input is unreadable");
        const auto continuous = directory/"adaptive-continuous.nc";
        const auto midpoint = directory/"adaptive-midpoint.nc";
        const auto restarted = directory/"adaptive-restarted.nc";
        const double midpointTime = initialCheckpoint.state.t+1e-5;
        const double finalTime = initialCheckpoint.state.t+2e-5;
        const auto adaptiveArguments = std::string(" --delta-t 1e-5 --integrator adaptive-rk23 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference");
        require(run(quote(input)+" "+quote(continuous)+adaptiveArguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"continuous adaptive restart control failed");
        require(run(quote(input)+" "+quote(midpoint)+adaptiveArguments+" --final-time "+number(midpointTime)+" >/dev/null 2>&1") == 0,"adaptive midpoint checkpoint failed");
        require(run(quote(midpoint)+" "+quote(restarted)+adaptiveArguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"adaptive restart continuation failed");
        WVCheckpoint continuousCheckpoint;
        WVCheckpoint restartedCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(continuous.string(), *test::extensionCatalog(),continuousCheckpoint)) && static_cast<bool>(WVCheckpointReader::read(restarted.string(), *test::extensionCatalog(),restartedCheckpoint)),"adaptive restart outputs are unreadable");
        require(relativeDifference(continuousCheckpoint,restartedCheckpoint) <= 1e-6,"adaptive restart exceeded tolerance-based trajectory equivalence");
        const auto rk45Continuous = directory/"adaptive-rk45-continuous.nc";
        const auto rk45Midpoint = directory/"adaptive-rk45-midpoint.nc";
        const auto rk45Restarted = directory/"adaptive-rk45-restarted.nc";
        const auto rk45Arguments = std::string(" --delta-t 1e-5 --integrator adaptive-rk45 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference");
        require(run(quote(input)+" "+quote(rk45Continuous)+rk45Arguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"continuous adaptive RK45 restart control failed");
        require(run(quote(input)+" "+quote(rk45Midpoint)+rk45Arguments+" --final-time "+number(midpointTime)+" >/dev/null 2>&1") == 0,"adaptive RK45 midpoint checkpoint failed");
        require(run(quote(rk45Midpoint)+" "+quote(rk45Restarted)+rk45Arguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"adaptive RK45 restart continuation failed");
        WVCheckpoint rk45ContinuousCheckpoint;
        WVCheckpoint rk45RestartedCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(rk45Continuous.string(), *test::extensionCatalog(),rk45ContinuousCheckpoint)) && static_cast<bool>(WVCheckpointReader::read(rk45Restarted.string(), *test::extensionCatalog(),rk45RestartedCheckpoint)),"adaptive RK45 restart outputs are unreadable");
        require(relativeDifference(rk45ContinuousCheckpoint,rk45RestartedCheckpoint) <= 1e-6,"adaptive RK45 restart exceeded tolerance-based trajectory equivalence");
        const auto rk78Continuous = directory/"adaptive-rk78-continuous.nc";
        const auto rk78Midpoint = directory/"adaptive-rk78-midpoint.nc";
        const auto rk78Restarted = directory/"adaptive-rk78-restarted.nc";
        const auto rk78Arguments = std::string(" --delta-t 1e-5 --initial-step 1e-5 --maximum-step 1e-5 --integrator adaptive-rk78 --relative-tolerance 1e-6 --absolute-tolerance 1e-8 --fft-provider reference");
        require(run(quote(input)+" "+quote(rk78Continuous)+rk78Arguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"continuous adaptive RK78 restart control failed");
        require(run(quote(input)+" "+quote(rk78Midpoint)+rk78Arguments+" --final-time "+number(midpointTime)+" >/dev/null 2>&1") == 0,"adaptive RK78 midpoint checkpoint failed");
        require(run(quote(rk78Midpoint)+" "+quote(rk78Restarted)+rk78Arguments+" --final-time "+number(finalTime)+" >/dev/null 2>&1") == 0,"adaptive RK78 restart continuation failed");
        WVCheckpoint rk78ContinuousCheckpoint;
        WVCheckpoint rk78RestartedCheckpoint;
        require(static_cast<bool>(WVCheckpointReader::read(rk78Continuous.string(), *test::extensionCatalog(),rk78ContinuousCheckpoint)) && static_cast<bool>(WVCheckpointReader::read(rk78Restarted.string(), *test::extensionCatalog(),rk78RestartedCheckpoint)),"adaptive RK78 restart outputs are unreadable");
        const auto rk78RestartDifference = relativeDifference(rk78ContinuousCheckpoint,rk78RestartedCheckpoint);
        require(rk78RestartDifference <= 1e-10,"adaptive RK78 checkpoint reconstruction exceeded its tolerance-equivalent trajectory: "+number(rk78RestartDifference));
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --steps 1 --final-time 9 --fft-provider reference >/dev/null 2>&1") != 0,"runner accepted both endpoint modes");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --steps 1 --fft-provider native-fftw >/dev/null 2>&1") != 0,"reference-only build silently substituted a provider");
        const auto sentinel = bytes(output);
        require(run(quote(input)+" "+quote(output)+" --delta-t bad --steps 1 --fft-provider reference >/dev/null 2>&1") != 0,"runner accepted malformed deltaT");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --steps 1 --integrator fixed-rk4 --relative-tolerance 1e-3 --fft-provider reference >/dev/null 2>&1") != 0,"fixed runner accepted adaptive tolerance options");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --initial-step 0.01 --steps 1 --integrator fixed-rk4 --fft-provider reference >/dev/null 2>&1") != 0,"fixed runner accepted adaptive step-control options");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --maximum-step bad --steps 1 --integrator adaptive-rk23 --fft-provider reference >/dev/null 2>&1") != 0,"adaptive runner accepted malformed maximum step");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --steps 1 --integrator adaptive-rk23 --benchmark-output-count 2 --fft-provider reference >/dev/null 2>&1") != 0,"adaptive runner accepted scheduled outputs without final-time integration");
        require(run(quote(input)+" "+quote(output)+" --delta-t 0.037 --steps 1 --integrator unknown --fft-provider reference >/dev/null 2>&1") != 0,"runner accepted an unknown integrator");
        const auto invalidSeriesDirectory = directory/"invalid-series";
        require(run(quote(input)+" "+quote(output)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-directory "+quote(invalidSeriesDirectory)+" >/dev/null 2>&1") != 0,"runner accepted positional OUTPUT with scheduled output");
        require(run(quote(input)+" --delta-t 1e-5 --steps 2 --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-directory "+quote(invalidSeriesDirectory)+" >/dev/null 2>&1") != 0,"runner accepted --steps with scheduled output");
        require(run(quote(input)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-time "+number(scheduledMidpoint)+" --output-directory "+quote(invalidSeriesDirectory)+" >/dev/null 2>&1") != 0,"runner accepted duplicate scheduled times");
        require(run(quote(input)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-directory "+quote(invalidSeriesDirectory)+" --output-pattern '../escape-{index}.nc' >/dev/null 2>&1") != 0,"runner accepted a scheduled-output path escape");
        const auto outputDirectoryFile = directory/"not-a-directory";
        {
            std::ofstream stream(outputDirectoryFile);
            stream << "file";
        }
        require(run(quote(input)+" --delta-t 1e-5 --final-time "+number(scheduledFinalTime)+" --fft-provider reference --output-time "+number(scheduledMidpoint)+" --output-directory "+quote(outputDirectoryFile)+" >/dev/null 2>&1") != 0,"runner accepted a file as the scheduled output directory");
        require(bytes(output) == sentinel,"argument failure changed an existing output");
        std::filesystem::remove_all(directory);
        std::cout << "Standalone runner tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
