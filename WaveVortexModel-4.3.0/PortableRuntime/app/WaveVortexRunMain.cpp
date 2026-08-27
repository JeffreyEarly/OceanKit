#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"

#include <iostream>
#include <memory>

int main(int argc, char **argv) {
  std::shared_ptr<const wavevortex::runtime::WVExtensionCatalog> catalog;
  const auto status = wavevortex::runtime::makeBuiltInExtensionCatalog(catalog);
  if (!status) {
    std::cerr << "Unable to construct built-in extension catalog: "
              << status.message << '\n';
    return 2;
  }
  return wavevortex::runtime::runWaveVortex(argc, argv, std::move(catalog));
}
