#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"

#include <algorithm>
#include <cmath>
#include <new>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime {
namespace {

struct Workspace {
  std::vector<double> weights, parameters, trialParameters, step;
  std::vector<double> gaps, logNodes, residual, trialResidual;
  std::vector<double> jacobian, trialJacobian, augmented, rightHandSide;
  std::vector<double> profile;

  std::size_t bytes() const noexcept {
    return sizeof(double) *
           (weights.capacity() + parameters.capacity() +
            trialParameters.capacity() + step.capacity() + gaps.capacity() +
            logNodes.capacity() + residual.capacity() + trialResidual.capacity() +
            jacobian.capacity() + trialJacobian.capacity() + augmented.capacity() +
            rightHandSide.capacity() + profile.capacity());
  }
};

bool finitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

// x_2 = -sum(exp(u)); successive log-density nodes add those positive gaps.
// Keep this order consistent with the MATLAB parameterization and Jacobian.
bool logDensityNodes(const std::vector<double> &parameters, Workspace &work) {
  double sum = 0.0;
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    if (!std::isfinite(parameters[i]))
      return false;
    work.gaps[i] = std::exp(parameters[i]);
    sum += work.gaps[i];
  }
  if (!std::isfinite(sum))
    return false;
  double prefix = 0.0;
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    work.logNodes[i] = -sum + prefix;
    prefix += work.gaps[i];
  }
  return true;
}

// The matrices use column-major storage. Only two O(n*p) objective arrays are
// needed; weighted exponential terms and their cumulative Jacobian sums are
// formed directly without another full matrix.
bool evaluate(const std::vector<double> &parameters,
              const std::vector<double> &targets, Workspace &work,
              std::vector<double> &residual, std::vector<double> &jacobian) {
  if (!logDensityNodes(parameters, work))
    return false;
  const auto n = targets.size();
  const auto p = parameters.size();
  for (std::size_t row = 0; row < n; ++row) {
    const double order = static_cast<double>(row + 1);
    double interior = 0.0;
    double derivativePrefix = 0.0;
    for (std::size_t column = 0; column < p; ++column) {
      const double weightedPower =
          std::exp(order * work.logNodes[column]) * work.weights[column + 1];
      interior += weightedPower;
      derivativePrefix += order * weightedPower;
      const double entry = -work.gaps[column] * derivativePrefix;
      if (!std::isfinite(entry))
        return false;
      jacobian[row + n * column] = entry;
    }
    residual[row] = interior + work.weights.back() - targets[row];
    if (!std::isfinite(residual[row]))
      return false;
  }
  return true;
}

double costOf(const std::vector<double> &residual) {
  double sum = 0.0;
  for (const double value : residual)
    sum += value * value;
  return sum / 2.0;
}

double gradientNorm(const std::vector<double> &jacobian,
                    const std::vector<double> &residual, std::size_t p) {
  double maximum = 0.0;
  const auto n = residual.size();
  for (std::size_t column = 0; column < p; ++column) {
    double entry = 0.0;
    for (std::size_t row = 0; row < n; ++row)
      entry += jacobian[row + n * column] * residual[row];
    if (!std::isfinite(entry))
      return std::numeric_limits<double>::infinity();
    maximum = std::max(maximum, std::abs(entry));
  }
  return maximum;
}

double vectorNorm(const std::vector<double> &values) {
  double norm = 0.0;
  for (const double value : values)
    norm = std::hypot(norm, value);
  return norm;
}

// Householder QR of [J; sqrt(lambda)*I], applied to [-r;0] in place.
// No normal equations, external linear-algebra dependency, or factorization
// workspace allocation is used. The strictly positive damping regularizes
// every parameter column, even when the raw-moment Jacobian is rank deficient.
bool solveAugmented(Workspace &work, double damping) {
  const auto n = work.residual.size();
  const auto p = work.parameters.size();
  const auto rows = n + p;
  const double diagonal = std::sqrt(damping);
  if (!finitePositive(diagonal))
    return false;
  std::fill(work.augmented.begin(), work.augmented.end(), 0.0);
  std::fill(work.rightHandSide.begin(), work.rightHandSide.end(), 0.0);
  for (std::size_t column = 0; column < p; ++column) {
    for (std::size_t row = 0; row < n; ++row)
      work.augmented[row + rows * column] = work.jacobian[row + n * column];
    work.augmented[n + column + rows * column] = diagonal;
  }
  for (std::size_t row = 0; row < n; ++row)
    work.rightHandSide[row] = -work.residual[row];

  for (std::size_t column = 0; column < p; ++column) {
    double norm = 0.0;
    for (std::size_t row = column; row < rows; ++row)
      norm = std::hypot(norm, work.augmented[row + rows * column]);
    if (!finitePositive(norm))
      return false;
    const double alpha = work.augmented[column + rows * column];
    const double beta = -std::copysign(norm, alpha);
    const double ratio = alpha / beta;
    const double tau = 1.0 - ratio;
    // Scaling by beta avoids overflow in the equivalent alpha-beta divisor.
    for (std::size_t row = column + 1; row < rows; ++row)
      work.augmented[row + rows * column] =
          (work.augmented[row + rows * column] / beta) / (ratio - 1.0);
    work.augmented[column + rows * column] = beta;
    for (std::size_t next = column + 1; next < p; ++next) {
      double dot = work.augmented[column + rows * next];
      for (std::size_t row = column + 1; row < rows; ++row)
        dot += work.augmented[row + rows * column] *
               work.augmented[row + rows * next];
      dot *= tau;
      work.augmented[column + rows * next] -= dot;
      for (std::size_t row = column + 1; row < rows; ++row)
        work.augmented[row + rows * next] -=
            work.augmented[row + rows * column] * dot;
    }
    double dot = work.rightHandSide[column];
    for (std::size_t row = column + 1; row < rows; ++row)
      dot += work.augmented[row + rows * column] * work.rightHandSide[row];
    dot *= tau;
    work.rightHandSide[column] -= dot;
    for (std::size_t row = column + 1; row < rows; ++row)
      work.rightHandSide[row] -= work.augmented[row + rows * column] * dot;
  }
  for (std::size_t reverse = p; reverse > 0; --reverse) {
    const auto row = reverse - 1;
    double value = work.rightHandSide[row];
    for (std::size_t column = row + 1; column < p; ++column)
      value -= work.augmented[row + rows * column] * work.step[column];
    work.step[row] = value / work.augmented[row + rows * row];
    if (!std::isfinite(work.step[row]))
      return false;
  }
  return true;
}

void updateReport(const Workspace &work, double cost, double damping,
                  WVNoMotionRecoveryReport &report) {
  report.finalCost = cost;
  report.damping = damping;
  report.maximumResidual = 0.0;
  for (const double value : work.residual)
    report.maximumResidual = std::max(report.maximumResidual, std::abs(value));
  report.gradientNorm = gradientNorm(work.jacobian, work.residual,
                                    work.parameters.size());
}

} // namespace

WVKernelStatus WVNoMotionProfileRecovery::fitMoments(
    const std::vector<double> &integrationWeights, double depth,
    const std::vector<double> &initialProfile, double minimumDensity,
    double maximumDensity, const std::vector<double> &targetMoments,
    std::vector<double> &output, WVNoMotionRecoveryReport &report,
    const WVNoMotionRecoveryOptions &options) {
  report = {};
  report.reason = "invalid-input";
  const auto n = integrationWeights.size();
  if (n < 3 || initialProfile.size() != n || targetMoments.size() != n)
    return {WVKernelStatusCode::invalidShape,
            "No-motion recovery requires n>=3 weights, profile knots and moments."};
  if (!finitePositive(depth) || !std::isfinite(minimumDensity) ||
      !std::isfinite(maximumDensity) || !(maximumDensity > minimumDensity) ||
      options.maximumEvaluations == 0 || !finitePositive(options.gradientTolerance) ||
      !finitePositive(options.stepTolerance) ||
      !finitePositive(options.relativeCostTolerance))
    return {WVKernelStatusCode::invalidConfiguration,
            "No-motion recovery requires ordered finite extrema and positive controls."};
  for (std::size_t i = 0; i < n; ++i) {
    if (!finitePositive(integrationWeights[i]) ||
        !std::isfinite(initialProfile[i]) || initialProfile[i] < 0.0 ||
        (i && !(initialProfile[i - 1] > initialProfile[i])) ||
        !std::isfinite(targetMoments[i]))
      return {WVKernelStatusCode::invalidConfiguration,
              "No-motion recovery inputs require finite moments and a strictly stable reference."};
  }
  const double densitySpan = maximumDensity - minimumDensity;
  const double initialSpan = initialProfile.front() - initialProfile.back();
  if (!finitePositive(densitySpan) || !finitePositive(initialSpan)) {
    report.reason = "nonfinite-normalization";
    return {WVKernelStatusCode::numericalFailure,
            "No-motion recovery density normalization is not representable."};
  }
  const auto p = n - 2;
  const auto limit = std::numeric_limits<std::size_t>::max();
  // Bound every matrix product and total owned byte count before allocation.
  if (n > limit - p || n > limit / p || n + p > limit / p) {
    report.reason = "size-overflow";
    return {WVKernelStatusCode::sizeOverflow, "No-motion recovery matrix size overflowed."};
  }
  const auto rows = n + p;
  const auto np = n * p;
  const auto augmentedCount = rows * p;
  const std::size_t counts[] = {n, p, p, p, p, p, n, n, np, np,
                                augmentedCount, rows, n};
  std::size_t total = 0;
  for (const auto count : counts) {
    if (count > limit / sizeof(double) - total) {
      report.reason = "size-overflow";
      return {WVKernelStatusCode::sizeOverflow, "No-motion recovery workspace size overflowed."};
    }
    total += count;
  }
  try {
    Workspace work;
    const auto allocate = [&](std::vector<double> &values, std::size_t count) {
      values.resize(count);
      report.workspaceBytes = work.bytes();
    };
    allocate(work.weights, n);
    allocate(work.parameters, p);
    allocate(work.trialParameters, p);
    allocate(work.step, p);
    allocate(work.gaps, p);
    allocate(work.logNodes, p);
    allocate(work.residual, n);
    allocate(work.trialResidual, n);
    allocate(work.jacobian, np);
    allocate(work.trialJacobian, np);
    allocate(work.augmented, augmentedCount);
    allocate(work.rightHandSide, rows);
    allocate(work.profile, n);

    report.reason = "nonfinite-normalization";
    for (std::size_t i = 0; i < n; ++i) {
      work.weights[i] = integrationWeights[n - 1 - i] / depth;
      if (!finitePositive(work.weights[i]))
        return {WVKernelStatusCode::numericalFailure,
                "No-motion recovery normalized weights are not representable."};
    }
    double previous = std::log((initialProfile[n - 2] - initialProfile.back()) /
                               initialSpan);
    for (std::size_t i = 0; i < p; ++i) {
      const double next = i + 1 == p ? 0.0 :
          std::log((initialProfile[n - 3 - i] - initialProfile.back()) / initialSpan);
      work.parameters[i] = std::log(next - previous);
      if (!std::isfinite(previous) || !std::isfinite(next) ||
          !std::isfinite(work.parameters[i]))
        return {WVKernelStatusCode::numericalFailure,
                "No-motion recovery initial log gaps are not representable."};
      previous = next;
    }
    report.reason = "nonfinite-initial-state";
    report.evaluations = 1;
    if (!evaluate(work.parameters, targetMoments, work, work.residual, work.jacobian)) {
      report.exitFlag = -1;
      return {WVKernelStatusCode::numericalFailure,
              "No-motion recovery initial moments or Jacobian are not finite."};
    }
    double cost = costOf(work.residual);
    report.initialCost = cost;
    double scale = 0.0;
    for (std::size_t column = 0; column < p; ++column) {
      double sum = 0.0;
      for (std::size_t row = 0; row < n; ++row) {
        const double entry = work.jacobian[row + n * column];
        sum += entry * entry;
      }
      scale = std::max(scale, sum);
    }
    const double minimum = std::numeric_limits<double>::min();
    double damping = std::max(1e-3 * scale, minimum);
    updateReport(work, cost, damping, report);
    if (!std::isfinite(cost) || !std::isfinite(scale) ||
        !std::isfinite(report.gradientNorm)) {
      report.exitFlag = -1;
      return {WVKernelStatusCode::numericalFailure,
              "No-motion recovery initial objective arithmetic overflowed."};
    }
    double multiplier = 2.0;
    report.reason = "iteration-limit";
    for (std::size_t iteration = 0; iteration < options.maximumIterations; ++iteration) {
      report.iterations = iteration + 1;
      const double gradient = gradientNorm(work.jacobian, work.residual, p);
      if (!std::isfinite(gradient)) {
        report.exitFlag = -1;
        report.reason = "nonfinite-gradient";
        break;
      }
      if (gradient <= options.gradientTolerance) {
        report.exitFlag = 1;
        report.reason = "gradient-tolerance";
        break;
      }
      if (report.evaluations >= options.maximumEvaluations) {
        report.reason = "evaluation-limit";
        break;
      }
      if (!solveAugmented(work, damping)) {
        report.exitFlag = -1;
        report.reason = "nonfinite-step";
        break;
      }
      const double stepNorm = vectorNorm(work.step);
      const double parameterNorm = vectorNorm(work.parameters);
      const double stepBound = options.stepTolerance * (parameterNorm + options.stepTolerance);
      if (!std::isfinite(stepNorm) || !std::isfinite(parameterNorm) ||
          !std::isfinite(stepBound)) {
        report.exitFlag = -1;
        report.reason = "nonfinite-step";
        break;
      }
      if (stepNorm <= stepBound) {
        report.exitFlag = 2;
        report.reason = "step-tolerance";
        break;
      }
      for (std::size_t i = 0; i < p; ++i)
        work.trialParameters[i] = work.parameters[i] + work.step[i];
      const bool trialFinite = evaluate(work.trialParameters, targetMoments, work,
                                        work.trialResidual, work.trialJacobian);
      ++report.evaluations;
      const double trialCost = trialFinite ? costOf(work.trialResidual) :
                                           std::numeric_limits<double>::infinity();
      double residualChange = 0.0;
      double squaredChange = 0.0;
      for (std::size_t row = 0; row < n; ++row) {
        double value = 0.0;
        for (std::size_t column = 0; column < p; ++column)
          value += work.jacobian[row + n * column] * work.step[column];
        residualChange += work.residual[row] * value;
        squaredChange += value * value;
      }
      const double predictedDecrease = -residualChange - squaredChange / 2.0;
      const double actualDecrease = cost - trialCost;
      if (trialFinite && std::isfinite(trialCost) && predictedDecrease > 0.0 &&
          actualDecrease > 0.0) {
        const double gain = actualDecrease / predictedDecrease;
        const double oldCost = cost;
        work.parameters.swap(work.trialParameters);
        work.residual.swap(work.trialResidual);
        work.jacobian.swap(work.trialJacobian);
        cost = trialCost;
        ++report.acceptedSteps;
        const double ratio = 2.0 * gain - 1.0;
        damping = std::max(minimum, damping * std::max(1.0 / 3.0,
                                                       1.0 - ratio * ratio * ratio));
        multiplier = 2.0;
        if (actualDecrease <= options.relativeCostTolerance * oldCost) {
          report.exitFlag = 3;
          report.reason = "relative-cost-tolerance";
          break;
        }
      } else {
        ++report.rejectedSteps;
        damping *= multiplier;
        multiplier *= 2.0;
        if (!std::isfinite(damping) || damping > 1e30 * std::max(scale, minimum)) {
          report.exitFlag = -2;
          report.reason = "damping-limit";
          break;
        }
      }
      if (report.evaluations >= options.maximumEvaluations) {
        report.reason = "evaluation-limit";
        break;
      }
    }
    updateReport(work, cost, damping, report);
    if (report.exitFlag <= 0 || !std::isfinite(report.maximumResidual) ||
        report.maximumResidual > 1e-8)
      return {WVKernelStatusCode::numericalFailure,
              "No-motion recovery did not qualify a normalized moment residual <=1e-8."};
    if (!logDensityNodes(work.parameters, work))
      return {WVKernelStatusCode::numericalFailure,
              "No-motion recovery final profile parameterization is not finite."};
    work.profile.front() = maximumDensity;
    work.profile.back() = minimumDensity;
    for (std::size_t i = 0; i < p; ++i)
      work.profile[n - 2 - i] = densitySpan * std::exp(work.logNodes[i]) + minimumDensity;
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(work.profile[i]) ||
          (i && !(work.profile[i - 1] > work.profile[i])))
        return {WVKernelStatusCode::numericalFailure,
                "No-motion recovery profile is not finite and strictly decreasing."};
    }
    output = std::move(work.profile);
    report.qualified = true;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    report.reason = "allocation-failure";
    return {WVKernelStatusCode::allocationFailure,
            "No-motion recovery workspace allocation failed."};
  } catch (const std::length_error &) {
    report.reason = "size-overflow";
    return {WVKernelStatusCode::sizeOverflow,
            "No-motion recovery workspace exceeds vector capacity."};
  }
}

} // namespace wavevortex::runtime
