#pragma once

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

} // namespace wavevortex::runtime
