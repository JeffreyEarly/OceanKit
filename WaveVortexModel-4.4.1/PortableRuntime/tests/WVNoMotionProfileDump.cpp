#include "WaveVortexRuntime/WVNoMotionProfile.hpp"
#include "nlohmann/json.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace wavevortex::runtime;
using nlohmann::json;

namespace {
void require(const wavevortex::WVKernelStatus &status) {
  if (!status) throw std::runtime_error(status.message);
}
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
}

int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("Usage: WVNoMotionProfileDump request.json result.json");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open no-motion profile request.");
    json request;
    input >> request;
    const auto z = request.at("z").get<std::vector<double>>();
    const auto rho = request.at("rho").get<std::vector<double>>();
    const auto queries = request.at("queryHeights").get<std::vector<double>>();
    const auto targets = request.value("targetDensity", std::vector<double>{});
    const auto material = request.value("materialHeights", std::vector<double>{});
    const bool evaluateInverse = request.contains("targetDensity");
    const bool evaluateAPE = request.contains("materialHeights");
    if (evaluateAPE && material.size() != queries.size())
      throw std::runtime_error("materialHeights must have one entry per queryHeights entry.");
    const double gravity = evaluateAPE ? request.at("g").get<double>() : 1;
    const double referenceDensity = evaluateAPE ? request.at("rho0").get<double>() : 1;
    WVNoMotionProfile profile;
    auto start = Clock::now();
    require(WVNoMotionProfile::create(z, rho, profile));
    const double constructionSeconds = elapsed(start);
    // Flat caller-owned vectors preserve MATLAB(:) order. The primitive
    // allocates no per-query storage; JSON parsing is outside these timings.
    std::vector<double> density(queries.size()), inverse(targets.size()), ape(material.size());
    start = Clock::now();
    for (std::size_t index = 0; index < queries.size(); ++index)
      require(profile.density(queries[index], density[index]));
    const double densitySeconds = elapsed(start);
    start = Clock::now();
    if (evaluateInverse)
      for (std::size_t index = 0; index < targets.size(); ++index)
        require(profile.inverse(targets[index], inverse[index]));
    const double inverseSeconds = elapsed(start);
    start = Clock::now();
    if (evaluateAPE)
      for (std::size_t index = 0; index < queries.size(); ++index)
        require(profile.availablePotentialEnergy(queries[index], material[index],
                                                 gravity, referenceDensity, ape[index]));
    const double apeSeconds = elapsed(start);
    json result{{"schema", "wv-no-motion-profile-dump-v1"}, {"status", "complete"},
                {"density", density}, {"persistentBytes", profile.retainedBytes()},
                {"knotCount", profile.knotCount()},
                {"seconds", {{"construction", constructionSeconds}, {"density", densitySeconds},
                             {"inverse", inverseSeconds}, {"ape", apeSeconds}}}};
    if (evaluateInverse) result["inverse"] = inverse;
    if (evaluateAPE) result["ape"] = ape;
    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error("Cannot open no-motion profile result.");
    output << result.dump() << '\n';
    if (!output) throw std::runtime_error("Cannot write no-motion profile result.");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
