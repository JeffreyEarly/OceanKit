#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef WV_BENCHMARK_SOURCE_REVISION
#define WV_BENCHMARK_SOURCE_REVISION "unknown"
#endif
#ifndef WV_BENCHMARK_SOURCE_DIRTY
#define WV_BENCHMARK_SOURCE_DIRTY 1
#endif
#ifndef WV_BENCHMARK_COMPILER
#define WV_BENCHMARK_COMPILER "unknown"
#endif
#ifndef WV_BENCHMARK_WORKER_SHA256
#define WV_BENCHMARK_WORKER_SHA256 "unknown"
#endif
#ifndef WV_BENCHMARK_CMAKE_SHA256
#define WV_BENCHMARK_CMAKE_SHA256 "unknown"
#endif
#ifndef WV_BENCHMARK_SOURCE_DIFF_SHA256
#define WV_BENCHMARK_SOURCE_DIFF_SHA256 "unknown"
#endif
#ifndef WV_BENCHMARK_SOURCE_STATUS_SHA256
#define WV_BENCHMARK_SOURCE_STATUS_SHA256 "unknown"
#endif

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;

namespace {

using Clock = std::chrono::steady_clock;

template <class Status> void require(Status status) {
  if (!status)
    throw std::runtime_error(status.message);
}

std::uint64_t hashBytes(const void *data, std::size_t bytes,
                        std::uint64_t hash = UINT64_C(1469598103934665603)) {
  const auto *values = static_cast<const unsigned char *>(data);
  for (std::size_t index = 0; index < bytes; ++index) {
    hash ^= values[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

std::uint64_t hashState(const WVState &state, std::size_t count) {
  auto hash = UINT64_C(1469598103934665603);
  for (const auto view : {state.coefficients.Ap, state.coefficients.Am,
                          state.coefficients.A0})
    hash = hashBytes(view.data, count * sizeof(WVComplex64), hash);
  return hash;
}

struct OutputIdentity {
  std::uint64_t flux = UINT64_C(1469598103934665603);
  std::uint64_t fields = UINT64_C(1469598103934665603);
  double maximumMagnitude = 0.0;
};

template <class Value>
double compareValue(Value value, Value reference) {
  const double difference = std::abs(value - reference);
  if (difference > 1e-12 + 1e-10 * std::abs(reference))
    throw std::runtime_error("Output changed beyond replay tolerance.");
  return difference;
}

OutputIdentity inspectOutputs(const std::vector<WVComplex64> &fp,
                              const std::vector<WVComplex64> &fm,
                              const std::vector<WVComplex64> &f0,
                              const std::vector<double> &fields) {
  OutputIdentity result;
  for (const auto *values : {&fp, &fm, &f0}) {
    for (const auto &value : *values) {
      if (!std::isfinite(value.real) || !std::isfinite(value.imag))
        throw std::runtime_error("Nonfinite flux output.");
      result.maximumMagnitude =
          std::max(result.maximumMagnitude, std::abs(value.real));
      result.maximumMagnitude =
          std::max(result.maximumMagnitude, std::abs(value.imag));
    }
    result.flux = hashBytes(values->data(),
                            values->size() * sizeof(WVComplex64), result.flux);
  }
  for (const auto value : fields) {
    if (!std::isfinite(value))
      throw std::runtime_error("Nonfinite physical-field output.");
    result.maximumMagnitude =
        std::max(result.maximumMagnitude, std::abs(value));
  }
  result.fields = hashBytes(fields.data(), fields.size() * sizeof(double),
                            result.fields);
  return result;
}

void writePayload(const std::string &path,
                  const std::vector<WVComplex64> &fp,
                  const std::vector<WVComplex64> &fm,
                  const std::vector<WVComplex64> &f0) {
  if (path.empty())
    return;
  if (std::filesystem::exists(path))
    throw std::runtime_error("Refusing to overwrite flux payload " + path);
  std::ofstream output(path, std::ios::binary);
  for (const auto *values : {&fp, &fm, &f0})
    output.write(reinterpret_cast<const char *>(values->data()),
                 static_cast<std::streamsize>(values->size() *
                                              sizeof(WVComplex64)));
  output.close();
  if (!output)
    throw std::runtime_error("Cannot write flux payload " + path);
}

json metrics(const WVHydrostaticKernelMetrics &value) {
  return {{"stateValidations", value.stateValidationCount},
          {"phasePreparations", value.phasePreparationCount},
          {"coefficientAssemblies", value.coefficientAssemblyCount},
          {"verticalPreparations", value.verticalPreparationCount},
          {"verticalOperatorExecutions", value.verticalOperatorExecutionCount},
          {"horizontalSpectrumReuses", value.horizontalSpectrumReuseCount},
          {"derivativeAdvectionConsumers",
           value.derivativeAdvectionConsumerCount},
          {"preparedVerticalDerivatives",
           value.preparedVerticalDerivativeCount},
          {"tiledNonlinearExecutions", value.tiledNonlinearCount},
          {"tiledColumnInverses", value.tiledColumnInverseCount},
          {"tiledRowInverses", value.tiledRowInverseCount},
          {"tiledReusedColumns", value.tiledReusedColumnCount}};
}

json metrics(const WVBoussinesqKernelMetrics &value) {
  auto result = json{{"stateValidations", value.stateValidationCount},
                     {"phasePreparations", value.phasePreparationCount},
                     {"coefficientAssemblies", value.coefficientAssemblyCount},
                     {"verticalPreparations", value.verticalPreparationCount},
                     {"verticalOperatorExecutions",
                      value.verticalOperatorExecutionCount},
                     {"verticalMatrixGroupExecutions",
                      value.verticalMatrixGroupExecutionCount},
                     {"horizontalSpectrumReuses",
                      value.horizontalSpectrumReuseCount},
                     {"derivativeAdvectionConsumers",
                      value.derivativeAdvectionConsumerCount},
                     {"preparedVerticalDerivatives",
                      value.preparedVerticalDerivativeCount},
                     {"tiledNonlinearExecutions", value.tiledNonlinearCount},
                     {"tiledColumnInverses", value.tiledColumnInverseCount},
                     {"tiledRowInverses", value.tiledRowInverseCount},
                     {"tiledReusedColumns", value.tiledReusedColumnCount}};
  return result;
}

template <class Kernel> class PreparedRunner {
public:
  PreparedRunner(std::unique_ptr<Kernel> kernel, const WVState &state,
                 std::size_t spectralCount, std::size_t spatialCount)
      : kernel_(std::move(kernel)), state_(state), spectralCount_(spectralCount),
        initialStateHash_(hashState(state_, spectralCount)), fp_(spectralCount),
        fm_(spectralCount), f0_(spectralCount), fields_(4 * spatialCount),
        flux_{{fp_.data(), kernel_->spectralShape()},
              {fm_.data(), kernel_->spectralShape()},
              {f0_.data(), kernel_->spectralShape()}},
        fieldView_{fields_.data(),
                   {kernel_->spatialShape().first, kernel_->spatialShape().second,
                    kernel_->spatialShape().third, 4}} {}

  json run(const json &command) {
    const auto warmups = command.at("warmups").get<std::size_t>();
    const auto samples = command.at("samples").get<std::size_t>();
    if (samples == 0 || warmups == 0)
      throw std::runtime_error("A run block needs at least one warmup and sample.");
    const auto capacity = kernel_->persistentBytes();
    maximumReplayDifference_ = 0.0;
    validationSeconds_ = 0.0;
    std::vector<double> warmupSeconds(warmups), sampleSeconds(samples);
    for (std::size_t index = 0; index < warmups; ++index)
      warmupSeconds[index] = execute(index + 1 == warmups);
    kernel_->resetMetrics();
    for (std::size_t index = 0; index < samples; ++index)
      sampleSeconds[index] = execute(index + 1 == samples);
    if (kernel_->persistentBytes() != capacity)
      throw std::runtime_error("Persistent capacity changed during run block.");
    if (hashState(state_, spectralCount_) != initialStateHash_)
      throw std::runtime_error("Input coefficients changed during run block.");
    const auto &storage = kernel_->storage();
    const auto producerMetrics = metrics(kernel_->metrics());
    if (producerMetrics.at("tiledNonlinearExecutions")
            .template get<std::size_t>() != samples)
      throw std::runtime_error("Not every timed sample used tiled nonlinear execution.");
    const auto payload = command.value("payload", std::string{});
    writePayload(payload, fp_, fm_, f0_);
    if (!payload.empty()) {
      const auto path = payload + ".fields.bin";
      if (std::filesystem::exists(path))
        throw std::runtime_error("Refusing to overwrite physical-field payload.");
      std::ofstream output(path, std::ios::binary);
      output.write(reinterpret_cast<const char*>(fields_.data()),
                   static_cast<std::streamsize>(fields_.size()*sizeof(double)));
      output.close();
      if (!output) throw std::runtime_error("Cannot write physical-field payload.");
    }
    return {{"event", "result"},
            {"id", command.value("id", std::string{})},
            {"warmupSeconds", warmupSeconds},
            {"sampleSeconds", sampleSeconds},
            {"warmups", warmups},
            {"samples", samples},
            {"persistentBytes", capacity},
            {"storage",
             {{"preparedBytes", storage.preparedBytes},
              {"workspaceBytes", storage.workspaceBytes},
              {"spectralScratchBytes", storage.spectralScratchBytes},
              {"realScratchBytes", storage.realScratchBytes},
              {"providerBytesLowerBound", storage.providerBytesLowerBound},
              {"planBytesLowerBound", storage.planBytesLowerBound}}},
            {"producerMetrics", producerMetrics},
            {"fluxHashFNV1a64", lastIdentity_.flux},
            {"fieldHashFNV1a64", lastIdentity_.fields},
            {"maximumOutputMagnitude", lastIdentity_.maximumMagnitude},
            {"maximumReplayDifference", maximumReplayDifference_},
            {"validationSeconds", validationSeconds_},
            {"validationPolicy", "last warmup and last sample per block"},
            {"payload", payload}};
  }

  Kernel &kernel() noexcept { return *kernel_; }

private:
  double execute(bool validate) {
    const auto start = Clock::now();
    auto status = kernel_->beginStateEvaluation(state_);
    const bool began = static_cast<bool>(status);
    if (began)
      status = kernel_->nonlinearFluxAndFields(state_, flux_, fieldView_);
    const auto endStatus = began ? kernel_->endStateEvaluation()
                                 : WVKernelStatus::ok();
    const auto seconds = std::chrono::duration<double>(Clock::now() - start).count();
    require(status);
    require(endStatus);
    if (!validate)
      return seconds;
    const auto validationStart = Clock::now();
    const auto now = inspectOutputs(fp_, fm_, f0_, fields_);
    lastIdentity_ = now;
    if (identity_) {
      for (std::size_t index = 0; index < spectralCount_; ++index) {
        for (const auto pair :
             {std::pair{fp_[index], referenceFp_[index]},
              std::pair{fm_[index], referenceFm_[index]},
              std::pair{f0_[index], referenceF0_[index]}}) {
          maximumReplayDifference_ =
              std::max(maximumReplayDifference_,
                       compareValue(pair.first.real, pair.second.real));
          maximumReplayDifference_ =
              std::max(maximumReplayDifference_,
                       compareValue(pair.first.imag, pair.second.imag));
        }
      }
      for (std::size_t index = 0; index < fields_.size(); ++index)
        maximumReplayDifference_ =
            std::max(maximumReplayDifference_,
                     compareValue(fields_[index], referenceFields_[index]));
    } else {
      identity_ = std::make_unique<OutputIdentity>(now);
      referenceFp_ = fp_;
      referenceFm_ = fm_;
      referenceF0_ = f0_;
      referenceFields_ = fields_;
    }
    validationSeconds_ +=
        std::chrono::duration<double>(Clock::now() - validationStart).count();
    return seconds;
  }

  std::unique_ptr<Kernel> kernel_;
  WVState state_;
  std::size_t spectralCount_ = 0;
  std::uint64_t initialStateHash_ = 0;
  std::vector<WVComplex64> fp_, fm_, f0_;
  std::vector<double> fields_;
  std::vector<WVComplex64> referenceFp_, referenceFm_, referenceF0_;
  std::vector<double> referenceFields_;
  WVFlux flux_;
  WVRealFieldBundleView fieldView_;
  std::unique_ptr<OutputIdentity> identity_;
  OutputIdentity lastIdentity_;
  double maximumReplayDifference_ = 0.0;
  double validationSeconds_ = 0.0;
};

template <class Runner>
void serve(Runner &runner, json ready) {
  const auto &kernel = runner.kernel();
  const auto identity = WVFFTWEngine::linkedLibraries();
  ready["event"] = "ready";
  ready["schema"] = "wvm-prepared-pipeline-v1";
  ready["source"] = {{"revision", WV_BENCHMARK_SOURCE_REVISION},
                     {"dirty", bool(WV_BENCHMARK_SOURCE_DIRTY)},
                     {"trackedDiffSha256", WV_BENCHMARK_SOURCE_DIFF_SHA256},
                     {"statusSha256", WV_BENCHMARK_SOURCE_STATUS_SHA256},
                     {"workerSha256", WV_BENCHMARK_WORKER_SHA256},
                     {"cmakeSha256", WV_BENCHMARK_CMAKE_SHA256}};
  ready["compiler"] = WV_BENCHMARK_COMPILER;
  ready["provider"] = {{"version", identity.version},
                       {"baseLibrary", identity.baseLibrary},
                       {"threadLibrary", identity.threadLibrary},
                       {"fftThreads", 1}};
  ready["matrixBackend"] = kernel.matrixBackendIdentifier();
  ready["horizontalSchedule"] = kernel.horizontalScheduleIdentifier();
  ready["grid"] = {kernel.geometry().Nx, kernel.geometry().Ny,
                   kernel.geometry().Nz};
  ready["measurement"] = "begin scope + nonlinearFluxAndFields + end scope; persisted forcing schedule excluded";
  std::cout << ready.dump() << '\n' << std::flush;

  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty())
      continue;
    const auto command = json::parse(line);
    const auto name = command.at("command").get<std::string>();
    if (name == "quit")
      return;
    if (name != "run")
      throw std::runtime_error("Unknown resident-worker command.");
    std::cout << runner.run(command).dump() << '\n' << std::flush;
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("Usage: wv-prepared-pipeline INPUT");
    const auto loadStart = Clock::now();
    std::shared_ptr<const WVExtensionCatalog> catalog;
    require(makeBuiltInExtensionCatalog(catalog));
    WVCheckpoint checkpoint;
    require(WVCheckpointReader::read(argv[1], *catalog, checkpoint));
    const auto loadSeconds =
        std::chrono::duration<double>(Clock::now() - loadStart).count();
    if (!checkpoint.stratifiedModalSource)
      throw std::runtime_error("A variable-stratification source is required.");
    const auto &geometry = checkpoint.stratifiedModalSource->geometry();
    if (geometry.transformClass != "WVTransformHydrostatic" &&
        geometry.transformClass != "WVTransformBoussinesq")
      throw std::runtime_error("Prepared pipeline supports Hydrostatic or Boussinesq input.");
    const auto spectralCount = geometry.Nj * geometry.Nkl;
    const auto spatialCount = geometry.Nx * geometry.Ny * geometry.Nz;
    const WVShape2D shape{geometry.Nj, geometry.Nkl};
    const auto coefficient = [&](const std::string &name) -> WVComplexConstView {
      for (const auto &family : checkpoint.transformState.coefficientFamilies)
        if (family.identifier == name) {
          if (family.values.size() != spectralCount)
            throw std::runtime_error("Fixture coefficient shape differs from modal source.");
          return {family.values.data(), shape};
        }
      throw std::runtime_error("Fixture is missing coefficient family " + name);
    };
    const WVState state{checkpoint.transformState.t,
                        checkpoint.transformState.t0,
                        {coefficient("Ap"), coefficient("Am"), coefficient("A0")}};
    WVVariableExecutionOptions options;
    options.horizontalSchedule =
        WVRetainedHorizontalSchedule::streamingPrunedTile16;
    options.horizontalWorkers = 12;
    options.pointwiseWorkers = 8;
    options.streamedNonlinear = true;
    options.spectralSchedule =
        WVVariableSpectralSchedule::compactSplitFusedViews;
    options.sharedFieldGradients = true;
    options.fusedDerivativeAdvection = true;
    options.tiledNonlinear = true;
    options.verticalGroupWorkers =
        geometry.transformClass == "WVTransformBoussinesq" ? 8 : 1;
    const auto preparationStart = Clock::now();
    std::unique_ptr<WVFFTEngine> fft;
    require(WVFFTWEngine::create(1, fft));
    json ready{{"input", argv[1]},
               {"family", geometry.transformClass},
               {"loadSeconds", loadSeconds},
               {"forcingEntriesPresent", checkpoint.forcingSchedule.entries.size()},
               {"options",
                {{"horizontalWorkers", options.horizontalWorkers},
                 {"pointwiseWorkers", options.pointwiseWorkers},
                 {"verticalGroupWorkers", options.verticalGroupWorkers},
                 {"compactSplitViews", options.usesCompactSplitViews()},
                 {"streamedNonlinear", options.streamedNonlinear},
                 {"sharedFieldGradients", options.sharedFieldGradients},
                 {"fusedDerivativeAdvection", options.fusedDerivativeAdvection},
                 {"tiledNonlinear", options.tiledNonlinear}}}};
    if (geometry.transformClass == "WVTransformHydrostatic") {
      std::unique_ptr<WVTransformHydrostaticKernel> kernel;
      require(WVTransformHydrostaticKernel::create(
          checkpoint.stratifiedModalSource, std::move(fft), kernel,
          WVCreateAccelerateMatrixBackend, options));
      ready["preparationSeconds"] =
          std::chrono::duration<double>(Clock::now() - preparationStart).count();
      if (!kernel->supportsTiledNonlinear())
        throw std::runtime_error("Hydrostatic tiled nonlinear path was not prepared.");
      PreparedRunner runner(std::move(kernel), state, spectralCount, spatialCount);
      serve(runner, std::move(ready));
    } else {
      std::unique_ptr<WVTransformBoussinesqKernel> kernel;
      require(WVTransformBoussinesqKernel::create(
          checkpoint.stratifiedModalSource, std::move(fft), kernel,
          WVCreateAccelerateMatrixBackend, options));
      ready["preparationSeconds"] =
          std::chrono::duration<double>(Clock::now() - preparationStart).count();
      if (!kernel->supportsTiledNonlinear())
        throw std::runtime_error("Boussinesq tiled nonlinear path was not prepared.");
      PreparedRunner runner(std::move(kernel), state, spectralCount, spatialCount);
      serve(runner, std::move(ready));
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
