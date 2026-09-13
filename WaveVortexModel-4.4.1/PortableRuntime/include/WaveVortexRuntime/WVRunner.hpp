#pragma once

#include "WaveVortexRuntime/WVIntegrationContracts.hpp"
#include <memory>

namespace wavevortex::runtime {

class WVExtensionCatalog;

// Stable source API v1 reusable runner entry point. catalog must be a non-null
// catalog frozen by the application after all built-in and source-linked
// registrations are added. The by-value shared pointer transfers retained
// ownership into the invocation; no process-global registry or plug-in
// discovery is performed.
int runWaveVortex(int argc, char **argv,
                  std::shared_ptr<const WVExtensionCatalog> catalog);

// The reusable runner never installs process signal handlers. Applications
// may supply their own stop control; the standalone main maps SIGINT here.
int runWaveVortex(int argc, char **argv,
                  std::shared_ptr<const WVExtensionCatalog> catalog,
                  const WVIntegrationControl &control);

} // namespace wavevortex::runtime
