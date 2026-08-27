#pragma once

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <cstddef>
#include <limits>

namespace wavevortex::runtime {

class WVConstantStratificationIntegrationSystem;

// Resolved MATLAB-compatible particle component. Accepted x and y state is
// deliberately unwrapped; periodicity is applied only while sampling fields.
class WVLagrangianParticles final {
public:
  const WVObserverRecord &record() const noexcept { return record_; }
  std::size_t particleCount() const noexcept { return particleCount_; }
  bool isXYOnly() const noexcept { return record_.isXYOnly; }
  std::size_t positionOffset() const noexcept { return positionOffset_; }

private:
  WVObserverRecord record_;
  std::size_t xBlock_ = 0;
  std::size_t yBlock_ = 0;
  std::size_t zBlock_ = std::numeric_limits<std::size_t>::max();
  std::size_t particleCount_ = 0;
  std::size_t positionOffset_ = 0;
  std::size_t uOutput_ = 0;
  std::size_t vOutput_ = 0;
  std::size_t wOutput_ = std::numeric_limits<std::size_t>::max();
  friend class WVConstantStratificationIntegrationSystem;
};

} // namespace wavevortex::runtime
