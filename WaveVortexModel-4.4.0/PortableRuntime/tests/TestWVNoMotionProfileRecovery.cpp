#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"
#include "WVAllocationProbe.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {
void require(bool value, const char *message) {
  if (!value) throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance, const char *message) {
  require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
std::vector<double> uniformVolume(const std::vector<double> &profile, std::size_t plane) {
  std::vector<double> values(profile.size() * plane);
  for (std::size_t z = 0; z < profile.size(); ++z)
    std::fill_n(values.begin() + z * plane, plane, profile[z]);
  return values;
}
// Direct long-double parcel quadrature is independent of the production
// normalized-power recurrence, tiling, and log-gap parameterization.
std::vector<double> profileMoments(const std::vector<double> &profile,
                                   const std::vector<double> &weights) {
  std::vector<double> result(profile.size());
  const long double low = profile.back(), span = profile.front() - low;
  for (std::size_t order = 1; order <= profile.size(); ++order) {
    long double sum = 0;
    for (std::size_t z = 0; z < profile.size(); ++z)
      sum += weights[z] * std::pow((profile[z] - low) / span, static_cast<int>(order));
    result[order - 1] = static_cast<double>(sum);
  }
  return result;
}
void qualified(const WVNoMotionRecoveryReport &report) {
  require(report.qualified && report.exitFlag > 0 &&
              std::isfinite(report.maximumResidual) && report.maximumResidual <= 1e-8,
          "successful recovery lacks physical qualification");
}

void verifyIndependentWeightedMoments() {
  const std::vector<double> weights{.2, .3, .5};
  const std::vector<double> values{10, 14, 12, 12, 14, 10};
  WVDensityDistribution distribution;
  require(bool(WVNoMotionProfileRecovery::moments({values.data(), {2, 1, 3}}, weights, 1, distribution)),
          "unequal-volume moment extraction failed");
  require(distribution.minimumDensity == 10 && distribution.maximumDensity == 14 &&
              distribution.sampleCount == 6 && distribution.stableProfile.empty(),
          "moment distribution metadata changed");
  require(distribution.moments.size() == 3, "wrong moment order count");
  for (std::size_t order = 1; order <= 3; ++order)
    close(distribution.moments[order - 1], .35 + .3 * std::pow(.5, static_cast<int>(order)),
          2e-16, "moments lost unequal vertical parcel volumes");
  // Each occupied plane contributes a representable weight, although dividing
  // that weight between its two parcels before summing would underflow to zero.
  const double tiny = std::numeric_limits<double>::denorm_min();
  const std::vector<double> tinyWeights(3, tiny);
  for (const std::size_t plane : {2U, 32771U}) {
    const auto binaryValues = uniformVolume({1, 0, 1}, plane);
    require(bool(WVNoMotionProfileRecovery::moments({binaryValues.data(), {plane, 1, 3}},
                                                   tinyWeights, 1, distribution)),
            "representable subnormal plane moments were rejected");
    require(distribution.stableProfile.empty() && distribution.moments.size() == 3,
            "unstable binary distribution took stable-rest shortcut");
    for (double moment : distribution.moments)
      require(moment == 2 * tiny, "parcel or partial-tile underflow erased a representable plane moment");
  }
}

void verifyExactChangedStableRest() {
  for (const auto &changed : {std::vector<double>{-1, -2},
                             std::vector<double>{12, 11.2, 10.6, 10.2, 10}}) {
    const std::size_t n = changed.size();
    const auto values = uniformVolume(changed, 6);
    std::vector<double> initial(n), weights(n, 1.0 / n), output{719};
    for (std::size_t z = 0; z < n; ++z) initial[z] = 100 - static_cast<double>(z);
    WVDensityDistribution distribution;
    require(bool(WVNoMotionProfileRecovery::moments({values.data(), {3, 2, n}}, weights, 1, distribution)),
            "stable rest moment scan failed");
    require(distribution.stableProfile == changed && distribution.moments.empty(),
            "stable shortcut reconstructed rather than preserved the profile");
    WVNoMotionRecoveryReport report;
    require(bool(WVNoMotionProfileRecovery::recover({values.data(), {3, 2, n}}, weights, 1,
                                                    initial, output, report)), "stable rest recovery failed");
    qualified(report);
    require(output == changed && report.iterations == 0 && report.evaluations == 0 &&
                report.maximumResidual == 0 && std::strcmp(report.algorithm, "stable-rest-profile") == 0,
            "stable-rest shortcut changed density or invoked fitting");
    require(report.workspaceBytes == n * sizeof(double), "stable shortcut workspace includes borrowed inputs");
    // Exact stable rest does not need normalized weights. Valid finite raw
    // weights must not reject that shortcut when their quotient overflows.
    weights.assign(n, 1e308);
    output = {719};
    require(bool(WVNoMotionProfileRecovery::recover({values.data(), {3, 2, n}}, weights, 1e-308,
                                                    initial, output, report)),
            "unused overflowing weight normalization rejected exact stable rest");
    qualified(report);
    require(output == changed && report.evaluations == 0 && report.maximumResidual == 0,
            "extreme finite weights changed exact stable-rest recovery");
  }
}

void verifyVolumePermutationRecovery() {
  const std::vector<double> initial{20, 19.75, 19.5, 19.25, 19};
  const std::vector<double> weights{.125, .25, .25, .25, .125};
  for (const auto &changed : {std::vector<double>{12, 11.5, 11, 10.5, 10},
                             std::vector<double>{11.7, 11.33125, 10.925, 10.48125, 10}}) {
    auto values = uniformVolume(changed, 6);
    // Exchange equal-volume parcels across levels: the true sorted profile
    // and all moments are unchanged, while the stable-rest shortcut is false.
    std::swap(values[6], values[18]);
    std::vector<double> output{719};
    WVNoMotionRecoveryReport report;
    require(bool(WVNoMotionProfileRecovery::recover({values.data(), {3, 2, 5}}, weights, 1,
                                                    initial, output, report)), "permuted profile recovery failed");
    qualified(report);
    require(output.size() == changed.size() && output.front() == changed.front() &&
                output.back() == changed.back(), "current extrema were replaced by initial extrema");
    for (std::size_t z = 0; z < changed.size(); ++z)
      close(output[z], changed[z], 2e-7, "recovery differs from known parcel permutation truth");
    const auto expected = profileMoments(changed, weights), actual = profileMoments(output, weights);
    for (std::size_t k = 0; k < expected.size(); ++k)
      close(actual[k], expected[k], 1e-8, "qualified output fails independent moment closure");
  }
}

void verifyInvalidDistributionsPreserveOutput() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> valid{3, 1, 2, 2, 1, 3};
  const std::vector<double> weights{.25, .5, .25};
  WVDensityDistribution output;
  output.minimumDensity = 719; output.maximumDensity = 720;
  output.moments = {721}; output.stableProfile = {722}; output.sampleCount = 723;
  output.workspaceBytes = 724;
  const auto unchanged = [&]() {
    require(output.minimumDensity == 719 && output.maximumDensity == 720 &&
                output.moments == std::vector<double>{721} && output.stableProfile == std::vector<double>{722} &&
                output.sampleCount == 723 && output.workspaceBytes == 724,
            "failed moment extraction changed caller output");
  };
  for (const auto shape : {WVShape3D{0, 1, 3}, WVShape3D{2, 1, 1},
                           WVShape3D{std::numeric_limits<std::size_t>::max(), 2, 3}}) {
    require(!WVNoMotionProfileRecovery::moments({valid.data(), shape}, weights, 1, output),
            "invalid moment shape accepted"); unchanged();
  }
  require(!WVNoMotionProfileRecovery::moments({nullptr, {2, 1, 3}}, weights, 1, output),
          "null volume accepted"); unchanged();
  for (const auto &bad : {std::vector<double>{1}, std::vector<double>{0, .5, .5},
                          std::vector<double>{-.1, .5, .6}, std::vector<double>{nan, .5, .5}}) {
    require(!WVNoMotionProfileRecovery::moments({valid.data(), {2, 1, 3}}, bad, 1, output),
            "invalid vertical weights accepted"); unchanged();
  }
  for (const double depth : {0.0, -1.0, nan}) {
    require(!WVNoMotionProfileRecovery::moments({valid.data(), {2, 1, 3}}, weights, depth, output),
            "invalid depth accepted"); unchanged();
  }
  for (const auto &bad : {std::vector<double>(6, 2), std::vector<double>{3, 1, 2, nan, 1, 3},
                          std::vector<double>{3, 1, 2, 2, 1, std::numeric_limits<double>::infinity()}}) {
    require(!WVNoMotionProfileRecovery::moments({bad.data(), {2, 1, 3}}, weights, 1, output),
            "constant or nonfinite density accepted"); unchanged();
  }
  allocationProbe::failAfter = 0;
  const auto status = WVNoMotionProfileRecovery::moments({valid.data(), {2, 1, 3}}, weights, 1, output);
  allocationProbe::failAfter = -1;
  require(status.code == WVKernelStatusCode::allocationFailure, "moment allocation failure not reported");
  unchanged();
}

void verifySolverFailuresAndBudgets() {
  const std::vector<double> weights{.25, .5, .25}, initial{3, 2, 1}, target{.55, .4, .33};
  const auto stableValues = uniformVolume(initial, 2);
  const WVRealVolumeConstView stableView{stableValues.data(), {2, 1, 3}};
  std::vector<double> output{719, 720};
  WVNoMotionRecoveryReport report;
  const auto failed = [&](const WVKernelStatus &status) {
    require(!status && !report.qualified && output == std::vector<double>({719, 720}),
            "failed solver published output or qualified");
  };
  for (const auto &bad : {std::vector<double>{3, 1}, std::vector<double>{3, 3, 1},
                          std::vector<double>{3, 2, -1},
                          std::vector<double>{3, std::numeric_limits<double>::quiet_NaN(), 1}}) {
    failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, bad, 1, 3, target, output, report));
    // Validation must precede even the otherwise successful stable shortcut.
    failed(WVNoMotionProfileRecovery::recover(stableView, weights, 1, bad, output, report));
    require(std::strcmp(report.reason, "invalid-reference") == 0 && report.evaluations == 0,
            "recover bypassed initial-reference validation");
  }
  failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 3, 3, target, output, report));
  failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, {}, output, report));
  auto badTarget = target; badTarget[1] = std::numeric_limits<double>::quiet_NaN();
  failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, badTarget, output, report));
  for (int which = 0; which < 4; ++which) {
    WVNoMotionRecoveryOptions options;
    if (which == 0) options.maximumEvaluations = 0;
    if (which == 1) options.gradientTolerance = 0;
    if (which == 2) options.stepTolerance = -1;
    if (which == 3) options.relativeCostTolerance = std::numeric_limits<double>::infinity();
    failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, target, output, report, options));
    failed(WVNoMotionProfileRecovery::recover(stableView, weights, 1, initial, output, report, options));
    require(std::strcmp(report.reason, "invalid-options") == 0 && report.evaluations == 0,
            "recover bypassed option validation");
  }
  for (bool iterations : {false, true}) {
    WVNoMotionRecoveryOptions options;
    if (iterations) options.maximumIterations = 0; else options.maximumEvaluations = 1;
    failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, target, output, report, options));
    require(report.exitFlag == 0 && report.evaluations == 1, "budget exhaustion was treated as convergence");
  }
  WVNoMotionRecoveryOptions loose;
  loose.gradientTolerance = 1e9;
  failed(WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, target, output, report, loose));
  require(report.exitFlag > 0 && report.maximumResidual > 1e-8,
          "test failed to exercise positive but physically unqualified solver termination");
  allocationProbe::failAfter = 0;
  const auto status = WVNoMotionProfileRecovery::fitMoments(weights, 1, initial, 1, 3, target, output, report);
  allocationProbe::failAfter = -1;
  failed(status);
  require(status.code == WVKernelStatusCode::allocationFailure, "fit allocation failure not reported");
}

void verifyBoundedMomentStorage() {
  const std::vector<double> profile{8, 7, 6, 5, 4}, weights(5, .2);
  std::size_t priorBytes = 0;
  for (const std::size_t plane : {4097U, 32771U}) {
    auto values = uniformVolume(profile, plane);
    std::swap(values.front(), values[4 * plane]);
    WVDensityDistribution distribution;
    allocationProbe::calls = 0; allocationProbe::counting = true;
    const auto status = WVNoMotionProfileRecovery::moments({values.data(), {plane, 1, 5}}, weights, 1, distribution);
    allocationProbe::counting = false;
    require(bool(status) && allocationProbe::calls <= 6, "moment scan allocated per tile or sample");
    require(distribution.workspaceBytes <= 65536 + 5 * profile.size() * sizeof(double),
            "moment scan retained more than its bounded tile and vertical vectors");
    if (priorBytes) require(distribution.workspaceBytes == priorBytes, "moment workspace grew with volume size");
    priorBytes = distribution.workspaceBytes;
    const auto expected = profileMoments(profile, weights);
    for (std::size_t order = 0; order < expected.size(); ++order)
      close(distribution.moments[order], expected[order], 2e-15, "tile boundary/tail changed parcel moments");
  }
}

void verifyNoIterationAllocations() {
  const std::vector<double> weights{.125, .25, .25, .25, .125};
  const std::vector<double> initial{12, 11.5, 11, 10.5, 10};
  const std::vector<double> truth{12, 11.2, 10.6, 10.2, 10};
  const auto target = profileMoments(truth, weights);
  std::size_t initialAllocations = 0, initialBytes = 0;
  for (bool iterate : {false, true}) {
    std::vector<double> output;
    WVNoMotionRecoveryReport report;
    allocationProbe::calls = 0; allocationProbe::counting = true;
    const auto status = WVNoMotionProfileRecovery::fitMoments(weights, 1, iterate ? initial : truth,
                                                            10, 12, target, output, report);
    allocationProbe::counting = false;
    const std::size_t calls = allocationProbe::calls;
    require(calls > 0 && calls <= 32 && report.workspaceBytes <= sizeof(double) * (16 * 25 + 64 * 5),
            "fit workspace is not bounded by vertical problem size");
    if (!iterate) {
      require(bool(status) && report.evaluations == 1, "known solution required an optimization step");
      initialAllocations = calls; initialBytes = report.workspaceBytes;
    } else {
      require(bool(status) && report.iterations > 1, "allocation test did not exercise repeated iterations");
      qualified(report);
      require(calls == initialAllocations && report.workspaceBytes == initialBytes,
              "solver allocated additional workspace during iterations");
    }
  }
}
} // namespace

int main() {
  try {
    verifyIndependentWeightedMoments();
    verifyExactChangedStableRest();
    verifyVolumePermutationRecovery();
    verifyInvalidDistributionsPreserveOutput();
    verifySolverFailuresAndBudgets();
    verifyBoundedMomentStorage();
    verifyNoIterationAllocations();
    std::cout << "No-motion recovery contracts passed: independent volume moments, exact stable rest, "
                 "known parcel permutations, qualification, transactional failures, budgets and bounded workspace.\n";
    return 0;
  } catch (const std::exception &error) {
    allocationProbe::counting = false; allocationProbe::failAfter = -1;
    std::cerr << error.what() << '\n';
    return 1;
  }
}
