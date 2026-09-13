#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"
#include "WVTestQuadraticSchedule.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>

using namespace wavevortex::runtime;

int main(int argc, char **argv) {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (status)
    status = wavevortex::runtime::test::registerQuadraticSchedule(builder);
  std::shared_ptr<const WVExtensionCatalog> catalog;
  if (status)
    status = builder.freeze(catalog);
  if (!status || !catalog ||
      catalog->outputSchedules().registration(
          wavevortex::runtime::test::quadraticScheduleType, 1) == nullptr) {
    std::cerr << "Unable to construct an extended runner catalog.\n";
    return EXIT_FAILURE;
  }

  if (argc > 1)
    return runWaveVortex(argc, argv, std::move(catalog));

  // Exercise the exact public entry point through its deterministic usage
  // path. The standalone executable remains a thin built-in-only caller.
  constexpr int usageExitCode = 2;
  if (runWaveVortex(argc, argv, std::move(catalog)) != usageExitCode) {
    std::cerr << "The extended runner entry did not return its usage code.\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: source-linked extended runner composition\n";
  return EXIT_SUCCESS;
}
