#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVRunner.hpp"

#include <csignal>
#include <iostream>
#include <memory>

namespace {
volatile std::sig_atomic_t interrupted = 0;
void requestStop(int) {
  interrupted = 1;
  std::signal(SIGINT, SIG_DFL);
}
} // namespace

int main(int argc, char **argv) {
  std::shared_ptr<const wavevortex::runtime::WVExtensionCatalog> catalog;
  const auto status = wavevortex::runtime::makeBuiltInExtensionCatalog(catalog);
  if (!status) {
    std::cerr << "Unable to construct built-in extension catalog: "
              << status.message << '\n';
    return 2;
  }
  interrupted = 0;
  const auto previous = std::signal(SIGINT, requestStop);
  wavevortex::runtime::WVIntegrationControl control;
  control.shouldStop = [](const auto &) { return interrupted != 0; };
  const int result = wavevortex::runtime::runWaveVortex(
      argc, argv, std::move(catalog), control);
  std::signal(SIGINT, previous);
  return result;
}
