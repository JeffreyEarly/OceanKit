#pragma once

#include "WaveVortexRuntime/WVDensityDiagnosticContract.hpp"
#include "WaveVortexRuntime/WVNoMotionProfile.hpp"
#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace wavevortex::runtime::detail {

struct WVDensityEventGeometry {
  const std::vector<double> *heights = nullptr;
  const std::vector<double> *integrationWeights = nullptr;
  const std::vector<double> *initialProfile = nullptr;
  double depth = 0.0;
  double gravity = 0.0;
  double referenceDensity = 0.0;
};

enum class WVDensityEventField : std::uint8_t { rhoNm = 1, etaTrue = 2, ape = 4 };

struct WVDensityEventOutput {
  WVDensityEventField field = WVDensityEventField::rhoNm;
  double *data = nullptr;
  std::size_t elementCount = 0;
};

struct WVDensityEventView {
  const double *data = nullptr;
  std::size_t elementCount = 0;
};

struct WVDensityEventMetrics {
  std::size_t evaluationCount = 0;
  std::size_t reuseCount = 0;
  std::size_t recoveryAttemptCount = 0;
  std::size_t recoveryCount = 0;
  std::size_t profileConstructionAttemptCount = 0;
  std::size_t profileConstructionCount = 0;
  std::size_t inverseAttemptCount = 0;
  std::size_t inversePassCount = 0;
  std::size_t inverseSampleCount = 0;
  std::size_t apeAttemptCount = 0;
  std::size_t apePassCount = 0;
  std::size_t apeSampleCount = 0;
  // Owned vector capacities only, excluding borrowed volume/geometry, caller
  // outputs and object/allocator metadata. High water includes nested recovery
  // and profile construction peaks while the already prepared state is live.
  std::size_t liveBytes = 0;
  std::size_t highWaterBytes = 0;
  std::size_t recoveryWorkspaceHighWaterBytes = 0;
  std::size_t profileConstructionWorkspaceHighWaterBytes = 0;
};

// One immutable event binding, with demand-driven derived state. The volume
// and geometry are borrowed and must remain valid AND unchanged until release.
// The owning event must create a new binding for every later state, even if
// times and coefficient/storage addresses are identical. No pointer-key cache
// or state mutation detection is provided here.
class WVDensityEventEvaluation final {
public:
  static constexpr std::uint8_t rhoNmDemand = 1;
  static constexpr std::uint8_t etaTrueDemand = 2;
  static constexpr std::uint8_t apeDemand = 4;

  static WVKernelStatus create(
      WVRealVolumeConstView density, WVDensityEventGeometry geometry,
      WVDensityDiagnosticContract contract, WVDensityEventEvaluation &output,
      const WVNoMotionRecoveryOptions &options = {});

  // Complete missing stages; existing completed stages survive a later-stage
  // failure. rhoNm always means actual recovery, independently of the selected
  // reference used by etaTrue/ape. A zero demand performs no numerical work.
  WVKernelStatus prepare(std::uint8_t demands);
  WVKernelStatus reserveStorage(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference);
  // Reuse-policy event extension. Existing completed results and their backing
  // storage remain unchanged; only storage for missing demands may grow.
  WVKernelStatus reserveAdditionalStorage(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference);
  // Prepare the minimum reusable volume storage for low-memory execution.
  // One derived field uses one volume; simultaneous eta/APE uses two because
  // both published values must coexist. Material heights are overwritten in
  // place after their final consumer.
  WVKernelStatus reserveLowMemoryStorage(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference);

  // Validate all views and aliasing first, prepare every demand, then publish
  // together. Failure never changes any caller buffer. Duplicate field requests
  // with disjoint outputs are supported without repeating numerical work.
  WVKernelStatus evaluate(const WVDensityEventOutput *outputs, std::size_t count);

  // Borrowed immutable results; unprepared fields return an empty view. Views
  // remain valid through further demand extensions, until release/recreation.
  WVDensityEventView view(WVDensityEventField field) const noexcept;
  WVDensityEventView materialHeights() const noexcept;
  const WVDensityEventMetrics &metrics() const noexcept { return metrics_; }
  const WVNoMotionRecoveryReport &recoveryReport() const noexcept { return recoveryReport_; }
  bool initialized() const noexcept { return initialized_; }

  // Frees all owned vector storage and clears the borrowed binding. Counters
  // and high-water values remain available; liveBytes becomes zero.
  void release() noexcept;
  // Clear this event binding and all ready flags while preserving vector
  // capacities prepared by the owning output arena.
  void resetRetainingCapacity() noexcept;
  // Drop only logical derived values while retaining this immutable binding
  // and any bounded low-memory capacities prepared for repeated calls.
  void discardDerived() noexcept;

private:
  WVKernelStatus recoverActual();
  WVKernelStatus prepareProfile();
  WVKernelStatus invert();
  WVKernelStatus formEta();
  WVKernelStatus formAPE();
  void account() noexcept;
  WVKernelStatus validateOutputs(const WVDensityEventOutput *, std::size_t,
                                 std::uint8_t &) const;

  WVRealVolumeConstView density_;
  WVDensityEventGeometry geometry_;
  WVDensityDiagnosticContract contract_;
  WVNoMotionRecoveryOptions options_;
  WVNoMotionRecoveryReport recoveryReport_;
  WVDensityEventMetrics metrics_;
  WVNoMotionProfile profile_;
  std::vector<double> actualProfile_, materialHeights_, etaTrue_, ape_;
  std::vector<double> lowPrimary_, lowSecondary_;
  std::size_t sampleCount_ = 0;
  std::uint8_t lowMemoryDemands_ = 0;
  int etaLowSlot_ = -1;
  int apeLowSlot_ = -1;
  bool lowMemoryStorage_ = false;
  bool initialized_ = false;
  bool actualReady_ = false;
  bool profileReady_ = false;
  bool inverseReady_ = false;
  bool etaReady_ = false;
  bool apeReady_ = false;
};

} // namespace wavevortex::runtime::detail
