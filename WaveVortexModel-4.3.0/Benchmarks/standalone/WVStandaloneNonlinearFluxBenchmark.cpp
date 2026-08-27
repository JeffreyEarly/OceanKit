#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WVNativeFFTWEngine.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
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
    if (character == '\\' || character == '"') output << '\\';
    output << character;
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

void emitFailure(const std::string &stage, const std::string &message) {
  std::cerr << "{\"schemaVersion\":\"three-interface-worker-v1\","
            << "\"status\":\"failed\",\"failure\":{\"stage\":\""
            << escape(stage) << "\",\"message\":\"" << escape(message)
            << "\"}}\n";
}

void phase(const std::string &path, const std::string &value,
           bool plateau = false) {
  if (path.empty()) return;
  const auto temporary = path + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    output << value;
  }
  std::error_code ignored;
  std::filesystem::rename(temporary, path, ignored);
  if (plateau) std::this_thread::sleep_for(std::chrono::milliseconds(100));
}
}  // namespace

int main(int argc, char **argv) {
  if (argc != 6 && argc != 8) {
    emitFailure("arguments", "Usage: worker INPUT THREADS WARMUPS SAMPLES OUTPUT [--phase-file PATH]");
    return 2;
  }
  if (argc == 8 && std::string(argv[6]) != "--phase-file") {
    emitFailure("arguments", "The optional argument must be --phase-file PATH.");
    return 2;
  }
  const std::string phasePath = argc == 8 ? argv[7] : "";
  phase(phasePath, "startup");
  std::size_t threads = 0, warmups = 0, samples = 0;
  if (!parseCount(argv[2], threads) || threads == 0 ||
      !parseCount(argv[3], warmups) || !parseCount(argv[4], samples) ||
      samples == 0) {
    emitFailure("arguments", "Thread and sample counts are invalid.");
    return 2;
  }

  WVCheckpoint checkpoint;
  std::shared_ptr<const WVExtensionCatalog> catalog;
  auto status = makeBuiltInExtensionCatalog(catalog);
  if (!status) {
    emitFailure("catalog", status.message);
    return 3;
  }
  phase(phasePath, "read");
  auto checkpointStatus =
      WVCheckpointReader::read(argv[1], *catalog, checkpoint);
  if (!checkpointStatus) {
    emitFailure("read", checkpointStatus.message);
    return 3;
  }
  if (checkpoint.forcingSchedule.entries.size() != 1 ||
      checkpoint.forcingSchedule.entries.front().typeIdentifier !=
          "WVNonlinearAdvection") {
    emitFailure("forcing", "The matched kernel case requires only nonlinear advection.");
    return 3;
  }

  std::unique_ptr<WVFFTEngine> fft;
  phase(phasePath, "construct");
  status = WVFFTWEngine::create(threads, fft);
  if (!status) {
    emitFailure("provider", status.message);
    return 4;
  }
  const auto identity = WVFFTWEngine::linkedLibraries();
  const auto expected = std::filesystem::weakly_canonical(
      std::filesystem::path(WV_RUNTIME_EXPECTED_FFTW_ROOT));
  const auto base = std::filesystem::weakly_canonical(identity.baseLibrary);
  const auto thread = std::filesystem::weakly_canonical(identity.threadLibrary);
  const auto underExpected = [&expected](const std::filesystem::path &value) {
    return value.string().rfind(expected.string() +
                                    std::string(1, std::filesystem::path::preferred_separator),
                                0) == 0;
  };
  if (identity.version.find("3.3.11") == std::string::npos ||
      !underExpected(base) || !underExpected(thread)) {
    emitFailure("provider", "Loaded FFTW does not match the pinned provider.");
    return 4;
  }

  std::unique_ptr<WVConstantStratificationForcingEngine> forcing;
  status = WVConstantStratificationForcingEngine::create(
      checkpoint.configuration, checkpoint.forcingSchedule, catalog,
      std::move(fft), forcing);
  if (!status) {
    emitFailure("construct", status.message);
    return 4;
  }
  const auto shape = forcing->stateShape();
  const auto count = shape.rows * shape.columns;
  std::vector<WVComplex64> fp(count), fm(count), f0(count);
  const WVState state{
      checkpoint.state.t, checkpoint.state.t0,
      {{checkpoint.state.coefficients.Ap.data(), shape},
       {checkpoint.state.coefficients.Am.data(), shape},
       {checkpoint.state.coefficients.A0.data(), shape}}};
  WVFlux flux{{fp.data(), shape}, {fm.data(), shape}, {f0.data(), shape}};
  for (std::size_t index = 0; index < warmups; ++index) {
    status = forcing->nonlinearFlux(state, flux);
    if (!status) {
      emitFailure("warmup", status.message);
      return 5;
    }
  }
  phase(phasePath, "steady-retained", true);
  const auto baselineRSS = currentRSSBytes();
  std::vector<double> raw(samples);
  phase(phasePath, "integrate");
  for (std::size_t index = 0; index < samples; ++index) {
    const auto start = Clock::now();
    status = forcing->nonlinearFlux(state, flux);
    raw[index] = std::chrono::duration<double>(Clock::now() - start).count();
    if (!status) {
      emitFailure("execute", status.message);
      return 5;
    }
  }
  const auto peakRSS = peakRSSBytes();
  auto sorted = raw;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  const auto &metrics = forcing->kernel().metrics();
  double checksum = 0.0;
  for (std::size_t index = 0; index < count; ++index)
    checksum += std::abs(fp[index].real) + std::abs(fp[index].imag) +
                std::abs(fm[index].real) + std::abs(fm[index].imag) +
                std::abs(f0[index].real) + std::abs(f0[index].imag);
  {
    std::ofstream output(argv[5], std::ios::binary | std::ios::trunc);
    if (!output) {
      emitFailure("output", "Unable to create the flux comparison file.");
      return 6;
    }
    output.write(reinterpret_cast<const char *>(fp.data()),
                 static_cast<std::streamsize>(fp.size() * sizeof(WVComplex64)));
    output.write(reinterpret_cast<const char *>(fm.data()),
                 static_cast<std::streamsize>(fm.size() * sizeof(WVComplex64)));
    output.write(reinterpret_cast<const char *>(f0.data()),
                 static_cast<std::streamsize>(f0.size() * sizeof(WVComplex64)));
    if (!output) {
      emitFailure("output", "Unable to write the flux comparison file.");
      return 6;
    }
  }
  phase(phasePath, "outputs-held", true);

  std::cout << std::setprecision(17)
            << "{\"schemaVersion\":\"three-interface-worker-v1\","
            << "\"status\":\"complete\",\"interface\":\"standalone-compiled\","
            << "\"operation\":\"nonlinearFlux\",\"sourceCommit\":\""
            << WV_RUNTIME_SOURCE_COMMIT << "\",\"provider\":{\"id\":\"native-neon-pthreads\","
            << "\"version\":\"" << escape(identity.version) << "\",\"threads\":"
            << threads << ",\"baseLibrary\":\"" << escape(identity.baseLibrary)
            << "\",\"threadLibrary\":\"" << escape(identity.threadLibrary)
            << "\",\"noFallback\":true},\"timing\":{\"medianSeconds\":" << median
            << ",\"samplesSeconds\":[";
  for (std::size_t index = 0; index < raw.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << raw[index];
  }
  std::cout << "]},\"memory\":{\"baselineProcessBytes\":" << baselineRSS
            << ",\"peakProcessBytes\":" << peakRSS
            << ",\"peakIncrementBytes\":"
            << (peakRSS > baselineRSS ? peakRSS - baselineRSS : 0)
            << "},\"execution\":{\"engine\":\"" << escape(forcing->kernel().engineIdentifier())
            << "\",\"library\":\"" << escape(forcing->kernel().engineLibraryIdentity())
            << "\",\"planCount\":" << metrics.planCount
            << ",\"checksum\":" << checksum << "}}\n";
  phase(phasePath, "complete");
  return 0;
}
