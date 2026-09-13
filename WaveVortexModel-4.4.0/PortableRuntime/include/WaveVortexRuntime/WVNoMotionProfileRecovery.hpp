#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cstddef>
#include <limits>
#include <vector>

namespace wavevortex::runtime {

struct WVNoMotionRecoveryOptions {
  std::size_t maximumIterations = 2000;
  std::size_t maximumEvaluations = 5000;
  double gradientTolerance = 1e-12;
  double stepTolerance = 1e-12;
  double relativeCostTolerance = 1e-12;
};

struct WVNoMotionRecoveryReport {
  const char *algorithm = "damped-least-squares";
  const char *reason = "not-started";
  int exitFlag = 0;
  std::size_t iterations = 0;
  std::size_t evaluations = 0;
  std::size_t acceptedSteps = 0;
  std::size_t rejectedSteps = 0;
  double initialCost = std::numeric_limits<double>::quiet_NaN();
  double finalCost = std::numeric_limits<double>::quiet_NaN();
  double maximumResidual = std::numeric_limits<double>::quiet_NaN();
  double gradientNorm = std::numeric_limits<double>::quiet_NaN();
  double damping = std::numeric_limits<double>::quiet_NaN();
  // Peak owned vector capacities, including unpublished profile staging;
  // excludes borrowed inputs, caller output, object/allocator metadata.
  std::size_t workspaceBytes = 0;
  // A positive solver exit alone does not qualify the physical profile.
  bool qualified = false;
};

struct WVDensityDistribution {
  double minimumDensity = std::numeric_limits<double>::quiet_NaN();
  double maximumDensity = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> moments;
  // Nonempty only for an exactly horizontally uniform, strictly stable field;
  // in that case moments is empty because no moment fit is needed.
  std::vector<double> stableProfile;
  std::size_t sampleCount = 0;
  std::size_t workspaceBytes = 0;
};

class WVNoMotionProfileRecovery final {
public:
  // Fit normalized raw moments of orders 1..n using n-2 log-gap parameters.
  // Initial shape comes from the initial profile; physical endpoints are the
  // supplied current extrema. Only a positive exit, maximum residual <=1e-8
  // and finite strictly decreasing result publish output. All failures leave
  // output unchanged; report preserves the solver termination and diagnostics.
  // Owned O(n^2) workspace is allocated once, with no iteration allocations.
  static WVKernelStatus fitMoments(
      const std::vector<double> &integrationWeights, double depth,
      const std::vector<double> &initialProfile, double minimumDensity,
      double maximumDensity, const std::vector<double> &targetMoments,
      std::vector<double> &output, WVNoMotionRecoveryReport &report,
      const WVNoMotionRecoveryOptions &options = {});

  static WVKernelStatus moments(
      WVRealVolumeConstView density,
      const std::vector<double> &integrationWeights, double depth,
      WVDensityDistribution &output,
      // Optional peak owned vector capacity, including partial failure work.
      // Excludes borrowed inputs and the prior caller output.
      std::size_t *workspaceBytes = nullptr);

  static WVKernelStatus recover(
      WVRealVolumeConstView density,
      const std::vector<double> &integrationWeights, double depth,
      const std::vector<double> &initialProfile, std::vector<double> &output,
      WVNoMotionRecoveryReport &report,
      const WVNoMotionRecoveryOptions &options = {});
};

} // namespace wavevortex::runtime
