#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime {
namespace {
constexpr std::size_t tileCapacity = 4096;

WVKernelStatus invalid(const char *message) {
  return {WVKernelStatusCode::invalidConfiguration, message};
}
WVKernelStatus numerical(const char *message) {
  return {WVKernelStatusCode::numericalFailure, message};
}
std::size_t retainedBytes(const WVDensityDistribution &distribution) {
  return sizeof(double) *
      (distribution.moments.capacity() + distribution.stableProfile.capacity());
}
} // namespace

WVKernelStatus WVNoMotionProfileRecovery::moments(
    WVRealVolumeConstView density, const std::vector<double> &integrationWeights,
    double depth, WVDensityDistribution &output, std::size_t *workspaceBytes) {
  if(workspaceBytes) *workspaceBytes=0;
  const auto nx = density.shape.first, ny = density.shape.second;
  const auto nz = density.shape.third;
  if (!nx || !ny || nz < 2 || integrationWeights.size() != nz)
    return {WVKernelStatusCode::invalidShape,
            "Density moments require a nonempty volume and matching vertical weights."};
  const auto maximum = std::numeric_limits<std::size_t>::max();
  if (nx > maximum / ny || nx * ny > maximum / nz ||
      nx * ny * nz > maximum / sizeof(double))
    return {WVKernelStatusCode::sizeOverflow, "Density volume size overflowed."};
  if (!density.data)
    return {WVKernelStatusCode::invalidPointer, "Density volume is null."};
  if (!std::isfinite(depth) || !(depth > 0.0))
    return invalid("Density moment depth must be finite and positive.");
  for (double weight : integrationWeights) {
    if (!std::isfinite(weight) || !(weight > 0.0))
      return invalid("Density moment weights must be finite and positive.");
  }
  try {
    WVDensityDistribution candidate;
    candidate.sampleCount = nx * ny * nz;
    candidate.stableProfile.resize(nz);
    if(workspaceBytes) *workspaceBytes=retainedBytes(candidate);
    candidate.minimumDensity = std::numeric_limits<double>::infinity();
    candidate.maximumDensity = -std::numeric_limits<double>::infinity();
    const auto plane = nx * ny;
    bool stable = true;
    for (std::size_t z = 0; z < nz; ++z) {
      const double first = density.data[z * plane];
      candidate.stableProfile[z] = first;
      if (z && !(first < candidate.stableProfile[z - 1])) stable = false;
      for (std::size_t point = 0; point < plane; ++point) {
        const double value = density.data[z * plane + point];
        if (!std::isfinite(value))
          return invalid("Density moment samples must be finite.");
        stable = stable && value == first;
        candidate.minimumDensity = std::min(candidate.minimumDensity, value);
        candidate.maximumDensity = std::max(candidate.maximumDensity, value);
      }
    }
    if (stable) {
      candidate.workspaceBytes = retainedBytes(candidate);
      output = std::move(candidate);
      return WVKernelStatus::ok();
    }
    candidate.stableProfile.clear();
    for (double weight : integrationWeights)
      if (!std::isfinite(weight / depth) || !(weight / depth > 0.0))
        return numerical("Normalized density moment weight is not representable.");
    const double range = candidate.maximumDensity - candidate.minimumDensity;
    if (!(range > 0.0) || !std::isfinite(range))
      return numerical("Current density distribution has no representable invertible range.");
    candidate.moments.assign(nz, 0.0);
    if(workspaceBytes) *workspaceBytes=retainedBytes(candidate);
    std::vector<double> compensation, planeMoments, planeCompensation, tile;
    const auto account=[&]() {
      candidate.workspaceBytes=retainedBytes(candidate)+sizeof(double)*
          (compensation.capacity()+planeMoments.capacity()+planeCompensation.capacity()+tile.capacity());
      if(workspaceBytes) *workspaceBytes=candidate.workspaceBytes;
    };
    compensation.resize(nz); account();
    planeMoments.resize(nz); account();
    planeCompensation.resize(nz); account();
    const auto tileSize = std::min(plane, tileCapacity);
    tile.resize(2 * tileSize); account();
    double *normalized = tile.data(), *power = tile.data() + tileSize;
    candidate.workspaceBytes = retainedBytes(candidate) +
        sizeof(double) * (compensation.capacity() + planeMoments.capacity() +
                          planeCompensation.capacity() + tile.capacity());
    // Borrow the volume. Only one bounded tile of normalized values and powers
    // is retained; no full-volume copy or parcel-by-moment matrix is formed.
    for (std::size_t z = 0; z < nz; ++z) {
      const double weight = integrationWeights[z] / depth;
      std::fill(planeMoments.begin(), planeMoments.end(), 0.0);
      std::fill(planeCompensation.begin(), planeCompensation.end(), 0.0);
      for (std::size_t begin = 0; begin < plane; begin += tileSize) {
        const auto count = std::min(tileSize, plane - begin);
        for (std::size_t point = 0; point < count; ++point) {
          normalized[point] = (density.data[z * plane + begin + point] -
                               candidate.minimumDensity) / range;
          power[point] = normalized[point];
        }
        for (std::size_t k = 0; k < nz; ++k) {
          // Independent short sums allow SIMD without fast-math or changing
          // the objective. A single serial reduction stalls each power pass.
          std::array<double, 8> lanes{};
          std::size_t point = 0;
          for (; point + lanes.size() <= count; point += lanes.size())
            for (std::size_t lane = 0; lane < lanes.size(); ++lane)
              lanes[lane] += power[point + lane];
          double sum = (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]) +
                       (lanes[4] + lanes[5]) + (lanes[6] + lanes[7]);
          for (; point < count; ++point) sum += power[point];
          const double increment = sum - planeCompensation[k];
          const double updated = planeMoments[k] + increment;
          planeCompensation[k] = (updated - planeMoments[k]) - increment;
          planeMoments[k] = updated;
          if (k + 1 < nz)
            for (std::size_t point = 0; point < count; ++point)
              power[point] *= normalized[point];
        }
      }
      // Apply the vertical weight after the full horizontal mean: weighting
      // partial tiles can underflow even when the plane's moment is representable.
      for (std::size_t k = 0; k < nz; ++k) {
        const double increment = weight * (planeMoments[k] / static_cast<double>(plane)) - compensation[k];
        const double updated = candidate.moments[k] + increment;
        compensation[k] = (updated - candidate.moments[k]) - increment;
        candidate.moments[k] = updated;
      }
    }
    for (double moment : candidate.moments)
      if (!std::isfinite(moment))
        return numerical("Density moment accumulation overflowed.");
    output = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure, "Density moment allocation failed."};
  } catch (const std::length_error &) {
    return {WVKernelStatusCode::sizeOverflow, "Density moment storage size overflowed."};
  }
}

WVKernelStatus WVNoMotionProfileRecovery::recover(
    WVRealVolumeConstView density, const std::vector<double> &integrationWeights,
    double depth, const std::vector<double> &initialProfile,
    std::vector<double> &output, WVNoMotionRecoveryReport &report,
    const WVNoMotionRecoveryOptions &options) {
  report = WVNoMotionRecoveryReport{};
  if (!options.maximumEvaluations || !std::isfinite(options.gradientTolerance) ||
      !(options.gradientTolerance > 0.0) || !std::isfinite(options.stepTolerance) ||
      !(options.stepTolerance > 0.0) || !std::isfinite(options.relativeCostTolerance) ||
      !(options.relativeCostTolerance > 0.0)) {
    report.reason = "invalid-options";
    return invalid("Density recovery requires positive finite tolerances and an evaluation budget.");
  }
  if (initialProfile.size() != density.shape.third || initialProfile.size() < 2) {
    report.reason = "invalid-reference";
    return {WVKernelStatusCode::invalidShape, "Initial density profile has the wrong vertical size."};
  }
  for (std::size_t i = 0; i < initialProfile.size(); ++i) {
    if (!std::isfinite(initialProfile[i]) || initialProfile[i] < 0.0 ||
        (i && !(initialProfile[i] < initialProfile[i - 1]))) {
      report.reason = "invalid-reference";
      return invalid("Initial density profile must be finite, nonnegative and strictly decreasing.");
    }
  }
  WVDensityDistribution distribution;
  auto status = moments(density, integrationWeights, depth, distribution,&report.workspaceBytes);
  if (!status) {
    report.reason = "invalid-density-distribution";
    return status;
  }
  if (!distribution.stableProfile.empty()) {
    // Match MATLAB's exact stable-rest path before any ill-conditioned fitting.
    report.algorithm = "stable-rest-profile";
    report.reason = "stable-rest-profile";
    report.exitFlag = 1;
    report.qualified = true;
    report.initialCost = report.finalCost = report.maximumResidual =
        report.gradientNorm = report.damping = 0.0;
    report.workspaceBytes = distribution.workspaceBytes;
    output = std::move(distribution.stableProfile);
    return WVKernelStatus::ok();
  }
  status = fitMoments(integrationWeights, depth, initialProfile,
                     distribution.minimumDensity, distribution.maximumDensity,
                     distribution.moments, output, report, options);
  report.workspaceBytes = std::max(distribution.workspaceBytes,
                                  retainedBytes(distribution) + report.workspaceBytes);
  return status;
}
} // namespace wavevortex::runtime
