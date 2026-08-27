#pragma once

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <cstddef>

namespace wavevortex::runtime {

class WVConstantStratificationIntegrationSystem;

// Resolved MATLAB-compatible three-dimensional tracer. Its numerical
// differentiation remains owned by the shared constant-stratification kernel.
class WVTracer final {
public:
  const WVObserverRecord &record() const noexcept { return record_; }
  std::size_t stateBlock() const noexcept { return stateBlock_; }
  bool shouldAntialias() const noexcept { return record_.shouldAntialias; }

private:
  WVObserverRecord record_;
  std::size_t stateBlock_ = 0;
  friend class WVConstantStratificationIntegrationSystem;
};

} // namespace wavevortex::runtime
