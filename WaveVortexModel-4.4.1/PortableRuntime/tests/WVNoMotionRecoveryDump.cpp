#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"
#include "nlohmann/json.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;

namespace {
std::size_t count(const json &value) {
  if (!value.is_number_integer() || value.get<double>() < 0 ||
      value.get<double>() >= static_cast<double>(std::numeric_limits<std::size_t>::max()))
    throw std::runtime_error("Expected a nonnegative representable integer.");
  return value.get<std::size_t>();
}
json reportJSON(const WVNoMotionRecoveryReport &r) {
  return {{"algorithm",r.algorithm},{"reason",r.reason},{"exitFlag",r.exitFlag},
          {"qualified",r.qualified},{"iterations",r.iterations},{"evaluations",r.evaluations},
          {"acceptedSteps",r.acceptedSteps},{"rejectedSteps",r.rejectedSteps},
          {"workspaceBytes",r.workspaceBytes},{"initialCost",r.initialCost},
          {"finalCost",r.finalCost},{"maximumResidual",r.maximumResidual},
          {"gradientNorm",r.gradientNorm},{"damping",r.damping}};
}
}

int main(int argc, char **argv) {
  try {
    if (argc != 3) throw std::runtime_error("Usage: WVNoMotionRecoveryDump request.json result.json");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open recovery request.");
    json request; input >> request;
    const auto mode = request.at("mode").get<std::string>();
    const auto weights = request.at("zInt").get<std::vector<double>>();
    const double depth = request.at("Lz").get<double>();
    const auto reference = request.value("rhoNM0",std::vector<double>{});
    WVNoMotionRecoveryOptions options;
    if (request.contains("options")) {
      const auto &o = request.at("options");
      if (o.contains("maximumIterations")) options.maximumIterations = count(o.at("maximumIterations"));
      if (o.contains("maximumEvaluations")) options.maximumEvaluations = count(o.at("maximumEvaluations"));
      options.gradientTolerance = o.value("gradientTolerance",options.gradientTolerance);
      options.stepTolerance = o.value("stepTolerance",options.stepTolerance);
      options.relativeCostTolerance = o.value("relativeCostTolerance",options.relativeCostTolerance);
    }
    WVRealVolumeConstView view;
    std::vector<double> density;
    if (mode == "moments" || mode == "recover") {
      const auto &shape = request.at("shape");
      if (!shape.is_array() || shape.size() != 3) throw std::runtime_error("shape requires Nx,Ny,Nz.");
      view.shape = {count(shape.at(0)),count(shape.at(1)),count(shape.at(2))};
      const auto nx = view.shape.first, ny = view.shape.second, nz = view.shape.third;
      const auto maximum = std::numeric_limits<std::size_t>::max();
      if (!nx || !ny || !nz || nx > maximum/ny || nx*ny > maximum/nz ||
          nx*ny*nz > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())/sizeof(double))
        throw std::runtime_error("Invalid or overflowing density shape.");
      const auto elements = nx*ny*nz;
      if (request.contains("densityFile")) {
        std::ifstream binary(request.at("densityFile").get<std::string>(),std::ios::binary|std::ios::ate);
        if (!binary || binary.tellg() != static_cast<std::streamoff>(elements*sizeof(double)))
          throw std::runtime_error("Binary density file length does not match shape.");
        density.resize(elements); binary.seekg(0);
        binary.read(reinterpret_cast<char *>(density.data()),static_cast<std::streamsize>(elements*sizeof(double)));
        if (!binary) throw std::runtime_error("Cannot read binary density volume.");
      } else density = request.at("rhoTotal").get<std::vector<double>>();
      if (density.size() != elements) throw std::runtime_error("rhoTotal length does not match shape.");
      view.data = density.data();
    }
    // Parsing, input ownership and filesystem IO are outside the measured call.
    WVDensityDistribution distribution;
    WVNoMotionRecoveryReport report;
    std::vector<double> profile;
    const auto target = request.value("targetMoments",std::vector<double>{});
    const double minimum = request.value("minimumDensity",0.0), maximum = request.value("maximumDensity",1.0);
    const auto started = std::chrono::steady_clock::now();
    WVKernelStatus status;
    if (mode == "moments")
      status = WVNoMotionProfileRecovery::moments(view,weights,depth,distribution);
    else if (mode == "recover")
      status = WVNoMotionProfileRecovery::recover(view,weights,depth,reference,profile,report,options);
    else if (mode == "fitMoments")
      status = WVNoMotionProfileRecovery::fitMoments(weights,depth,reference,minimum,maximum,target,profile,report,options);
    else throw std::runtime_error("Unknown recovery mode.");
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    json result{{"schema","wv-no-motion-recovery-dump-v1"},{"status",status ? "complete":"failed"},
                {"code",static_cast<unsigned>(status.code)},{"message",status.message},
                {"profile",profile},{"report",reportJSON(report)},{"elapsedSeconds",elapsed}};
    if (mode == "moments")
      result["distribution"] = {{"minimumDensity",distribution.minimumDensity},
          {"maximumDensity",distribution.maximumDensity},{"moments",distribution.moments},
          {"stableProfile",distribution.stableProfile},{"sampleCount",distribution.sampleCount},
          {"workspaceBytes",distribution.workspaceBytes}};
    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error("Cannot open recovery result.");
    output << result.dump() << '\n';
    if (!output) throw std::runtime_error("Cannot write recovery result.");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
