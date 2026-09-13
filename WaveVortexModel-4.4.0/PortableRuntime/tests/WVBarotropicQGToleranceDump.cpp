#include "WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

double argument(char **values, int index) {
  return std::stod(values[index]);
}

void require(const WVKernelStatus &status) {
  if (!status)
    throw std::runtime_error(status.message);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 13) {
    std::cerr << "usage: tolerance-dump Nx Ny Lx Ly h j g rotationRate "
                 "latitude shouldAntialias planetaryRadius "
                 "absoluteToleranceScale\n";
    return 2;
  }
  try {
    WVTransformBarotropicQGConfiguration configuration;
    configuration.Nx = static_cast<std::size_t>(argument(argv, 1));
    configuration.Ny = static_cast<std::size_t>(argument(argv, 2));
    configuration.Lx = argument(argv, 3);
    configuration.Ly = argument(argv, 4);
    configuration.h = argument(argv, 5);
    configuration.j = static_cast<std::uint32_t>(argument(argv, 6));
    configuration.g = argument(argv, 7);
    configuration.rotationRate = argument(argv, 8);
    configuration.latitude = argument(argv, 9);
    configuration.shouldAntialias = argument(argv, 10) != 0.0;
    configuration.planetaryRadius = argument(argv, 11);
    const double absoluteToleranceScale = argument(argv, 12);

    std::unique_ptr<WVBarotropicQGIntegrationSystem> system;
    require(WVBarotropicQGIntegrationSystem::create(
        configuration, std::make_unique<WVReferenceFFTEngine>(), system));
    std::unique_ptr<WVIntegrationErrorPolicy> policy;
    require(system->createErrorPolicy(absoluteToleranceScale, policy));

    std::cout << std::setprecision(17)
              << "{\"absoluteToleranceScale\":"
              << absoluteToleranceScale << ",\"absoluteTolerance\":[";
    for (std::size_t index = 0; index < policy->elementCount(0); ++index) {
      if (index != 0)
        std::cout << ',';
      std::cout << policy->absoluteTolerance(0, index);
    }
    std::cout << "]}\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 3;
  }
}
