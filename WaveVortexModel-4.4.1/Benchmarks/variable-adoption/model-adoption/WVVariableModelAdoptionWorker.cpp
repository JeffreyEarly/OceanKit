#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include "WVNativeFFTWEngine.hpp"

#if WV_MODEL_VARIABLE_SERVICES
#include "WaveVortexRuntime/WVVariableKernelServices.hpp"
#include "WaveVortexKernel/WVVariableExecutionOptions.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#endif

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
#include <sstream>
#include <string>
#include <vector>

#if defined(__APPLE__)
#include <mach/mach.h>
#include <sys/resource.h>
#endif

namespace {
using namespace wavevortex;
using namespace wavevortex::runtime;
using Clock = std::chrono::steady_clock;

std::string escape(const std::string &value) {
  std::ostringstream output;
  for (const unsigned char character : value) {
    switch (character) {
    case '\\':
    case '"': output << '\\' << character; break;
    case '\n': output << "\\n"; break;
    case '\r': output << "\\r"; break;
    case '\t': output << "\\t"; break;
    default: output << character;
    }
  }
  return output.str();
}

std::size_t currentRSSBytes() {
#if defined(__APPLE__)
  mach_task_basic_info information{};
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&information), &count) ==
      KERN_SUCCESS)
    return static_cast<std::size_t>(information.resident_size);
#endif
  return 0;
}

std::size_t peakRSSBytes() {
#if defined(__APPLE__)
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) == 0)
    return static_cast<std::size_t>(usage.ru_maxrss);
#endif
  return 0;
}

bool parseCount(const char *text, std::size_t &value) {
  if (text == nullptr || *text == '\0' || *text == '-') return false;
  char *end = nullptr;
  const auto parsed = std::strtoull(text, &end, 10);
  if (end == text || *end != '\0' ||
      parsed > std::numeric_limits<std::size_t>::max())
    return false;
  value = static_cast<std::size_t>(parsed);
  return true;
}

bool parseFinite(const char *text, double &value) {
  if (text == nullptr || *text == '\0') return false;
  char *end = nullptr;
  value = std::strtod(text, &end);
  return end != text && *end == '\0' && std::isfinite(value);
}

std::string transformIdentifier(WVPersistedTransformKind kind) {
  switch (kind) {
  case WVPersistedTransformKind::stratifiedQG: return "stratified-qg";
  case WVPersistedTransformKind::hydrostatic: return "hydrostatic";
  case WVPersistedTransformKind::boussinesq: return "boussinesq";
  case WVPersistedTransformKind::barotropicQG: return "barotropic-qg";
  case WVPersistedTransformKind::constantStratification:
    return "constant-stratification";
  }
  return "unknown";
}

void emitFailure(const std::string &stage, const std::string &message) {
  std::cerr << "{\"schemaVersion\":\"wvm-variable-model-worker-v1\","
            << "\"status\":\"failed\",\"failure\":{\"stage\":\""
            << escape(stage) << "\",\"message\":\"" << escape(message)
            << "\"}}\n";
}

bool writeReport(const std::string &path, const std::string &contents,
                 std::string &message) {
  if (std::filesystem::exists(path)) {
    message = "Refusing to replace an existing report.";
    return false;
  }
  const auto temporary = path + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) {
      message = "Unable to create the temporary report.";
      return false;
    }
    output << contents << '\n';
    if (!output) {
      message = "Unable to write the temporary report.";
      return false;
    }
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  if (error) {
    std::filesystem::remove(temporary);
    message = "Unable to publish the report: " + error.message();
    return false;
  }
  return true;
}

std::size_t retainedBytes(const WVModelMetrics &metrics) {
  return metrics.modelPersistentBytes + metrics.catalogPersistentBytes +
         metrics.statePersistentBytes + metrics.integrationSystemPersistentBytes +
         metrics.integratorPersistentBytes + metrics.outputPersistentBytes;
}

#if WV_MODEL_VARIABLE_SERVICES
struct BackendAudit {
  std::size_t instanceCount = 0;
  std::vector<std::string> identifiers;
};

WVVariableKernelServices servicesFor(const std::string &selection,
                                     std::size_t workers,
                                     std::size_t pointwiseWorkers,
                                     const std::shared_ptr<BackendAudit> &audit) {
  WVVariableKernelServices services;
  const bool scalar = selection == "frozen";
  services.matrixBackendFactory = [scalar, audit](
                                      std::unique_ptr<WVVerticalMatrixBackend> &backend) {
    const auto status = scalar ? WVCreateScalarMatrixBackend(backend)
                               : WVCreateAccelerateMatrixBackend(backend);
    if (status && backend) {
      ++audit->instanceCount;
      audit->identifiers.emplace_back(backend->identifier());
    }
    return status;
  };
  if (selection == "frozen") {
    services.execution.horizontalSchedule = WVRetainedHorizontalSchedule::fullFFT;
    services.execution.horizontalWorkers = 1;
    services.execution.streamedNonlinear = false;
    services.execution.spectralSchedule =
        WVVariableSpectralSchedule::establishedInterleaved;
  } else {
    services.execution.horizontalSchedule =
        WVRetainedHorizontalSchedule::streamingPrunedTile16;
    services.execution.horizontalWorkers = workers;
    services.execution.streamedNonlinear = true;
    services.execution.pointwiseWorkers = pointwiseWorkers;
    services.execution.spectralSchedule =
        selection == "compact"
            ? WVVariableSpectralSchedule::compactSplitFusedViews
            : WVVariableSpectralSchedule::establishedInterleaved;
  }
  return services;
}
#endif
} // namespace

int main(int argc, char **argv) {
  if (argc != 7 && argc != 8) {
    emitFailure("arguments", "Usage: worker OUTPUT_NC SELECTION HORIZONTAL_WORKERS FINAL_TIME STEP REPORT_JSON [POINTWISE_WORKERS]");
    return 2;
  }
  const auto processStarted = Clock::now();
  const std::string outputPath = argv[1];
  const std::string selection = argv[2];
  const std::string reportPath = argv[6];
  std::size_t workers = 0, requestedPointwiseWorkers = 1;
  double finalTime = 0.0, step = 0.0;
  if ((selection != "frozen" && selection != "interleaved" &&
       selection != "compact") ||
      !parseCount(argv[3], workers) || workers == 0 ||
      (argc == 8 && (!parseCount(argv[7], requestedPointwiseWorkers) ||
                     requestedPointwiseWorkers == 0)) ||
      !parseFinite(argv[4], finalTime) || !parseFinite(argv[5], step) ||
      step <= 0.0) {
    emitFailure("arguments", "Selection, worker count, final time, or step is invalid.");
    return 2;
  }
#if !WV_MODEL_VARIABLE_SERVICES
  if (selection != "frozen") {
    emitFailure("arguments", "This source build supports only the frozen selection.");
    return 2;
  }
#endif
  if (selection == "frozen" && requestedPointwiseWorkers != 1) {
    emitFailure("arguments", "The frozen selection requires one pointwise worker.");
    return 2;
  }
  if (!std::filesystem::is_regular_file(outputPath)) {
    emitFailure("arguments", "The copied model-output file does not exist.");
    return 2;
  }
  if (std::filesystem::exists(reportPath)) {
    emitFailure("arguments", "Refusing to replace an existing report.");
    return 2;
  }

  std::shared_ptr<const WVExtensionCatalog> catalog;
  auto status = makeBuiltInExtensionCatalog(catalog);
  if (!status) {
    emitFailure("catalog", status.message);
    return 3;
  }

  const auto inspectStarted = Clock::now();
  WVModelOutputNetCDFInspection inspection;
  auto checkpointStatus =
      WVModelOutputNetCDFSink::inspect({outputPath}, *catalog, inspection);
  if (!checkpointStatus) {
    emitFailure("inspect", checkpointStatus.message);
    return 3;
  }
  const double inspectSeconds =
      std::chrono::duration<double>(Clock::now() - inspectStarted).count();
  const auto transform = transformIdentifier(inspection.latestRestart.transformKind);
  if (transform != "stratified-qg" && transform != "hydrostatic" &&
      transform != "boussinesq") {
    emitFailure("inspect", "The harness admits only variable-stratification fixtures.");
    return 3;
  }
  if (inspection.latestRestart.forcingSchedule.entries.size() != 1 ||
      inspection.latestRestart.forcingSchedule.entries.front().typeIdentifier !=
          "WVNonlinearAdvection") {
    emitFailure("inspect", "The harness requires exactly nonlinear advection forcing.");
    return 3;
  }
  const double initialTime = inspection.latestRestart.t;
  if (!(finalTime > initialTime)) {
    emitFailure("arguments", "Final time must be later than the selected restart.");
    return 2;
  }

  WVModelOutputRequest request;
  request.policy = WVModelOutputPolicy::append;
  request.finalTime = finalTime;
  WVModelOutputConfiguration outputConfiguration;
  const auto preparationStarted = Clock::now();
  status = WVModel::prepareModelOutput(catalog, inspection, request,
                                       outputConfiguration);
  if (!status) {
    emitFailure("prepare-output", status.message);
    return 4;
  }
  const double preparationSeconds =
      std::chrono::duration<double>(Clock::now() - preparationStarted).count();

  const auto providerStarted = Clock::now();
  std::unique_ptr<WVFFTEngine> fft;
  status = WVFFTWEngine::create(1, fft);
  if (!status) {
    emitFailure("provider", status.message);
    return 4;
  }
  const auto providerIdentity = WVFFTWEngine::linkedLibraries();
  const auto expected = std::filesystem::weakly_canonical(
      std::filesystem::path(WV_RUNTIME_EXPECTED_FFTW_ROOT));
  const auto base =
      std::filesystem::weakly_canonical(providerIdentity.baseLibrary);
  const auto thread =
      std::filesystem::weakly_canonical(providerIdentity.threadLibrary);
  const auto underExpected = [&expected](const std::filesystem::path &value) {
    return value.string().rfind(
               expected.string() +
                   std::string(1, std::filesystem::path::preferred_separator),
               0) == 0;
  };
  if (providerIdentity.version.find("3.3.11") == std::string::npos ||
      !underExpected(base) || !underExpected(thread)) {
    emitFailure("provider", "Loaded FFTW does not match the pinned provider.");
    return 4;
  }
  const double providerSeconds =
      std::chrono::duration<double>(Clock::now() - providerStarted).count();

  WVModelIntegratorConfiguration integrator;
  integrator.kind = WVModelIntegratorKind::fixedRK4;
  integrator.fixed.retainDenseOutput = true;
  WVModel model;
  WVModelState state;
  std::string matrixBackend;
  std::size_t matrixBackendInstances = 0;
  const auto constructionStarted = Clock::now();
#if WV_MODEL_VARIABLE_SERVICES
  auto audit = std::make_shared<BackendAudit>();
  const auto services = servicesFor(selection, workers,
                                    requestedPointwiseWorkers, audit);
  status = WVModel::createFromModelOutputInspection(
      catalog, std::move(inspection), std::move(outputConfiguration),
      std::move(fft), integrator, model, state, {}, services);
  matrixBackendInstances = audit->instanceCount;
  if (!audit->identifiers.empty()) {
    matrixBackend = audit->identifiers.front();
    for (const auto &identifier : audit->identifiers) {
      if (identifier != matrixBackend) {
        emitFailure("construct", "Vertical operators received different matrix backends.");
        return 4;
      }
    }
  }
#else
  status = WVModel::createFromModelOutputInspection(
      catalog, std::move(inspection), std::move(outputConfiguration),
      std::move(fft), integrator, model, state, {});
  std::unique_ptr<WVVerticalMatrixBackend> probedBackend;
  if (status) {
    const auto backendStatus = WVCreateScalarMatrixBackend(probedBackend);
    if (!backendStatus || !probedBackend) {
      emitFailure("backend-probe", backendStatus.message);
      return 4;
    }
    matrixBackend = probedBackend->identifier();
  }
#endif
  if (!status) {
    emitFailure("construct", status.message);
    return 4;
  }
  const double constructionSeconds =
      std::chrono::duration<double>(Clock::now() - constructionStarted).count();
  if (matrixBackend.empty()) {
    emitFailure("construct", "No vertical matrix backend was constructed.");
    return 4;
  }

  const auto steadyRSS = currentRSSBytes();
  const auto integrationStarted = Clock::now();
  status = model.advanceToTime(state, finalTime, step);
  const double integrationSeconds =
      std::chrono::duration<double>(Clock::now() - integrationStarted).count();
  if (!status) {
    emitFailure("integrate", status.message);
    return 5;
  }
  const auto metricsBeforeClose = model.metrics(&state);
  const auto closeStarted = Clock::now();
  checkpointStatus = model.closeOutput();
  const double closeSeconds =
      std::chrono::duration<double>(Clock::now() - closeStarted).count();
  if (!checkpointStatus) {
    emitFailure("close-output", checkpointStatus.message);
    return 6;
  }
  const auto metrics = model.metrics(&state);
  const auto finalState = state.constView();
  if (finalState.waveVortex.t != finalTime) {
    emitFailure("integrate", "Final accepted time differs from the request.");
    return 5;
  }

  const bool compact = selection == "compact";
  const bool streamed = selection != "frozen";
  const char *horizontalSchedule =
      streamed ? "streaming-pruned-tile16" : "full-fft";
  const char *spectralSchedule = compact ? "compact-split-fused-views"
                                         : "established-interleaved";
  const auto effectiveWorkers = streamed ? workers : std::size_t{1};
  const auto effectivePointwiseWorkers =
      streamed ? requestedPointwiseWorkers : std::size_t{1};
  const auto peakRSS = peakRSSBytes();
  const double completeSeconds =
      std::chrono::duration<double>(Clock::now() - processStarted).count();

  std::ostringstream output;
  output << std::setprecision(17)
         << "{\"schemaVersion\":\"wvm-variable-model-worker-v1\","
         << "\"status\":\"complete\",\"interface\":\"standalone-cpp-wvmodel\","
         << "\"source\":{\"commit\":\"" << WV_MODEL_SOURCE_COMMIT
         << "\",\"root\":\"" << escape(WV_MODEL_SOURCE_ROOT)
         << "\",\"compiledSourceDirty\":"
         << (WV_MODEL_COMPILED_SOURCE_DIRTY ? "true" : "false")
         << ",\"workerSourceSHA256\":\"" << WV_MODEL_WORKER_SOURCE_SHA256
         << "\",\"buildType\":\"" << escape(WV_MODEL_BUILD_TYPE)
         << "\",\"variableServicesAvailable\":"
         << (WV_MODEL_VARIABLE_SERVICES ? "true" : "false")
         << ",\"compiler\":\"" << escape(__VERSION__) << "\"},"
         << "\"workload\":{\"transform\":\"" << transform
         << "\",\"initialTime\":" << initialTime
         << ",\"finalTime\":" << finalTime << ",\"step\":" << step
         << ",\"forcing\":\"WVNonlinearAdvection\","
         << "\"outputPath\":\"" << escape(outputPath) << "\"},"
         << "\"selection\":{\"id\":\"" << selection
         << "\",\"matrixBackendIdentifier\":\""
         << escape(matrixBackend) << "\",\"matrixBackendInstanceCount\":"
         << matrixBackendInstances
         << ",\"matrixBackendEvidence\":\""
#if WV_MODEL_VARIABLE_SERVICES
         << "injected-factory-observed-construction"
#else
         << "frozen-model-default-plus-matched-scalar-factory-probe"
#endif
         << "\",\"horizontalSchedule\":\"" << horizontalSchedule
         << "\",\"horizontalWorkers\":" << effectiveWorkers
         << ",\"streamedNonlinear\":" << (streamed ? "true" : "false")
         << ",\"spectralSchedule\":\"" << spectralSchedule
         << "\",\"compactSplitViews\":" << (compact ? "true" : "false")
         << ",\"pointwiseWorkers\":" << effectivePointwiseWorkers
         << ",\"generalVerticalWorkers\":1,\"fftwInternalThreads\":1,"
         << "\"pointwiseWorkersExposed\":"
         << (WV_MODEL_VARIABLE_SERVICES ? "true" : "false") << "},"
         << "\"provider\":{\"modelIdentifier\":\""
         << escape(model.kernelProviderIdentifier()) << "\",\"modelLibrary\":\""
         << escape(model.kernelProviderLibraryIdentity())
         << "\",\"version\":\"" << escape(providerIdentity.version)
         << "\",\"baseLibrary\":\"" << escape(providerIdentity.baseLibrary)
         << "\",\"threadLibrary\":\"" << escape(providerIdentity.threadLibrary)
         << "\",\"forcingSchedule\":\""
         << escape(model.forcingScheduleIdentifier()) << "\",\"noFallback\":true},"
         << "\"timingSeconds\":{\"inspect\":" << inspectSeconds
         << ",\"outputPreparation\":" << preparationSeconds
         << ",\"providerCreation\":" << providerSeconds
         << ",\"modelConstruction\":" << constructionSeconds
         << ",\"integrate\":" << integrationSeconds << ",\"close\":"
         << closeSeconds << ",\"workerCompleteBeforeReport\":" << completeSeconds
         << "},\"work\":{\"integratorSteps\":" << metrics.integrator.stepCount
         << ",\"acceptedSteps\":" << metrics.integrator.acceptedStepCount
         << ",\"rejectedSteps\":" << metrics.integrator.rejectedStepCount
         << ",\"rhsEvaluations\":" << metrics.integrator.rightHandSideEvaluationCount
         << ",\"forcingEvaluations\":" << metrics.forcing.evaluationCount
         << ",\"scalarAdvections\":"
         << metrics.integratedObservers.tracerEvaluationCount
         << ",\"kernelScalarMetricAvailable\":false"
         << ",\"tracerEvaluations\":" << metrics.integratedObservers.tracerEvaluationCount
         << ",\"particleVelocityEvaluations\":"
         << metrics.integratedObservers.velocityFieldEvaluationCount
         << ",\"denseOutputEvaluations\":" << metrics.integrator.denseOutputEvaluationCount
         << ",\"outputInterpolations\":"
         << metrics.outputDriver.interpolatedStateEvaluationCount
         << ",\"committedDeliveries\":" << metrics.outputDriver.committedDeliveryCount
         << ",\"committedRecords\":" << metrics.output.committedRecordCount
         << ",\"writtenBytes\":" << metrics.output.writtenBytes << "},"
         << "\"storageBytes\":{\"model\":" << metrics.modelPersistentBytes
         << ",\"catalog\":" << metrics.catalogPersistentBytes
         << ",\"state\":" << metrics.statePersistentBytes
         << ",\"integrationSystem\":" << metrics.integrationSystemPersistentBytes
         << ",\"integrator\":" << metrics.integratorPersistentBytes
         << ",\"output\":" << metrics.outputPersistentBytes
         << ",\"knownRetained\":" << retainedBytes(metrics)
         << ",\"knownRetainedBeforeClose\":" << retainedBytes(metricsBeforeClose)
         << ",\"kernelScratchCapacity\":" << metrics.kernel.scratchCapacityBytes
         << ",\"kernelScratchHighWater\":" << metrics.kernel.scratchHighWaterBytes
         << ",\"forcingWorkspaceCapacity\":" << metrics.forcing.workspaceCapacityBytes
         << ",\"forcingWorkspaceMaximumLive\":" << metrics.forcing.workspaceMaximumLiveBytes
         << ",\"integratorWorkspaceCapacity\":" << metrics.integrator.workspaceCapacityBytes
         << ",\"integratorWorkspaceMaximumLive\":"
         << metrics.integrator.workspaceMaximumLiveBytes
         << ",\"outputInterpolationCapacity\":"
         << metrics.outputDriver.interpolationBufferCapacityBytes
         << ",\"fullModelRetainedBeforeClose\":"
         << retainedBytes(metricsBeforeClose)
         << ",\"fullModelRetainedAfterClose\":" << retainedBytes(metrics)
         << ",\"steadyProcessRSS\":" << steadyRSS
         << ",\"peakProcessRSS\":" << peakRSS << "}}";

  std::string message;
  if (!writeReport(reportPath, output.str(), message)) {
    emitFailure("report", message);
    return 6;
  }
  std::cout << output.str() << '\n';
  return 0;
}
