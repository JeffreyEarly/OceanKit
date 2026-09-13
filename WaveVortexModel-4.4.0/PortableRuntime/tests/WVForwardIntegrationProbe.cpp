#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"

#include <charconv>
#include <iostream>
#include <string>
#include <vector>

using namespace wavevortex::runtime;

// Qualification-only application: stop through the public source API while
// leaving ordinary runner options, reports and destination policy unchanged.
int main(int argc, char **argv) {
  if (argc < 6 || std::string(argv[1]) != "--stop-boundary" ||
      std::string(argv[3]) != "--stop-after") {
    std::cerr << "Usage: WVForwardIntegrationProbe --stop-boundary "
                 "accepted-step|output-occurrence --stop-after N "
                 "<ordinary runner arguments>\n";
    return 2;
  }
  const std::string name = argv[2];
  if (name != "accepted-step" && name != "output-occurrence") {
    std::cerr << "Unknown stop boundary.\n";
    return 2;
  }
  const std::string count = argv[4];
  std::size_t stopAfter = 0;
  const auto parsed = std::from_chars(count.data(), count.data() + count.size(),
                                      stopAfter);
  if (parsed.ec != std::errc{} || parsed.ptr != count.data() + count.size() ||
      stopAfter == 0) {
    std::cerr << "Stop count must be a positive integer.\n";
    return 2;
  }
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  std::shared_ptr<const WVExtensionCatalog> catalog;
  if (status) status = builder.freeze(catalog);
  if (!status) {
    std::cerr << status.message << '\n';
    return 3;
  }
  std::vector<char *> arguments{argv[0]};
  arguments.insert(arguments.end(), argv + 5, argv + argc);
  const auto boundary = name == "accepted-step"
                            ? WVIntegrationBoundary::acceptedStep
                            : WVIntegrationBoundary::outputOccurrence;
  std::size_t matchingCalls = 0;
  return runWaveVortex(
      static_cast<int>(arguments.size()), arguments.data(), std::move(catalog),
      {[&](const WVIntegrationProgress &progress) {
        return progress.boundary == boundary && ++matchingCalls == stopAfter;
      }});
}
