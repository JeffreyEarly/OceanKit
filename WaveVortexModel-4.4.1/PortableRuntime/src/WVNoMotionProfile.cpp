#include "WaveVortexRuntime/WVNoMotionProfile.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime {
namespace {
WVKernelStatus invalid(const char *message) {
  return {WVKernelStatusCode::invalidConfiguration, message};
}
WVKernelStatus numerical(const char *message) {
  return {WVKernelStatusCode::numericalFailure, message};
}

// MATLAB eps(x) is spacing at the magnitude of x, including subnormals and
// zero. frexp also handles the largest finite double without nextafter -> inf.
double spacing(double magnitude) {
  if (magnitude < std::numeric_limits<double>::min())
    return std::numeric_limits<double>::denorm_min();
  int exponent = 0;
  std::frexp(magnitude, &exponent);
  return std::scalbn(1.0, exponent - std::numeric_limits<double>::digits);
}

double polynomial(const std::array<double, 4> &c, double x) {
  return ((c[0] * x + c[1]) * x + c[2]) * x + c[3];
}
} // namespace

WVKernelStatus WVNoMotionProfile::create(const std::vector<double> &heights,
                                         const std::vector<double> &densities,
                                         WVNoMotionProfile &output,
                                         std::size_t *creationWorkspaceBytes) {
  if (creationWorkspaceBytes)
    *creationWorkspaceBytes = 0;
  const auto count = heights.size();
  if (count < 2 || count != densities.size())
    return {WVKernelStatusCode::invalidShape,
            "No-motion profile requires matching vectors with at least two knots."};
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(heights[i]) || !std::isfinite(densities[i]))
      return invalid("No-motion profile knots must be finite.");
    if (i && (!(heights[i] > heights[i - 1]) ||
              !(densities[i] < densities[i - 1])))
      return invalid("No-motion heights must increase and densities decrease strictly.");
  }
  try {
    WVNoMotionProfile candidate;
    std::vector<double> widths, secants, slopes;
    const auto account = [&]() {
      if (creationWorkspaceBytes)
        *creationWorkspaceBytes = std::max(*creationWorkspaceBytes,
            candidate.retainedBytes() - sizeof(candidate) + sizeof(double) *
                (widths.capacity() + secants.capacity() + slopes.capacity()));
    };
    candidate.densityOffset_ = densities.back();
    candidate.densityScale_ = densities.front() - densities.back();
    const double span = heights.back() - heights.front();
    const double heightMagnitude =
        std::max(std::abs(heights.front()), std::abs(heights.back())) + span;
    if (!std::isfinite(candidate.densityScale_) ||
        !(candidate.densityScale_ > 0.0) || !std::isfinite(heightMagnitude))
      return numerical("No-motion profile normalization or height range overflowed.");
    candidate.densityTolerance_ =
        8.0 * spacing(std::max(std::abs(densities.front()),
                               std::abs(densities.back())));
    candidate.heightTolerance_ = 8.0 * spacing(heightMagnitude);
    candidate.heights_ = heights;
    account();
    candidate.densities_ = densities;
    account();
    candidate.normalizedDensity_.resize(count);
    account();
    candidate.coefficients_.resize(count - 1);
    account();
    widths.resize(count - 1);
    account();
    secants.resize(count - 1);
    account();
    slopes.resize(count);
    account();
    for (std::size_t i = 0; i < count; ++i) {
      const double normalized =
          (densities[i] - candidate.densityOffset_) / candidate.densityScale_;
      if (!std::isfinite(normalized) ||
          (i && !(normalized < candidate.normalizedDensity_[i - 1])))
        return numerical("No-motion density normalization lost a distinct knot.");
      candidate.normalizedDensity_[i] = normalized;
      if (i) {
        widths[i - 1] = heights[i] - heights[i - 1];
        secants[i - 1] =
            (normalized - candidate.normalizedDensity_[i - 1]) / widths[i - 1];
        if (!std::isfinite(secants[i - 1]) || !(secants[i - 1] < 0.0))
          return numerical("No-motion profile secant is not representable.");
      }
    }
    if (count == 2) {
      slopes[0] = slopes[1] = secants[0];
    } else {
      // Shape-preserving PCHIP uses a weighted harmonic mean of neighboring
      // secants. Strictly decreasing knots guarantee equal, nonzero signs.
      for (std::size_t i = 1; i + 1 < count; ++i) {
        const double leftWeight = 2.0 * widths[i] + widths[i - 1];
        const double rightWeight = widths[i] + 2.0 * widths[i - 1];
        const double numerator = leftWeight + rightWeight;
        const double denominator =
            leftWeight / secants[i - 1] + rightWeight / secants[i];
        if (!std::isfinite(numerator) || !std::isfinite(denominator))
          return numerical("No-motion PCHIP slope arithmetic overflowed.");
        slopes[i] = numerator / denominator;
      }
      const auto endpointSlope = [&](std::size_t near, std::size_t next) {
        return ((2.0 * widths[near] + widths[next]) * secants[near] -
                widths[near] * secants[next]) /
               (widths[near] + widths[next]);
      };
      slopes.front() = endpointSlope(0, 1);
      slopes.back() = endpointSlope(count - 2, count - 3);
      if (!std::isfinite(slopes.front()) || !std::isfinite(slopes.back()))
        return numerical("No-motion PCHIP endpoint arithmetic overflowed.");
      // Endpoints may have zero derivative while the profile remains strictly
      // monotone. The inverse must then use its bracket safeguard.
      slopes.front() = std::min(0.0, slopes.front());
      slopes.back() = std::min(0.0, slopes.back());
    }
    for (std::size_t i = 0; i + 1 < count; ++i) {
      const double leftCurvature = (secants[i] - slopes[i]) / widths[i];
      const double rightCurvature = (slopes[i + 1] - secants[i]) / widths[i];
      auto &c = candidate.coefficients_[i];
      c = {{(rightCurvature - leftCurvature) / widths[i],
            2.0 * leftCurvature - rightCurvature, slopes[i],
            candidate.normalizedDensity_[i]}};
      if (!std::all_of(c.begin(), c.end(),
                       [](double value) { return std::isfinite(value); }))
        return numerical("No-motion PCHIP coefficients are not representable.");
    }
    output = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "No-motion profile allocation failed."};
  } catch (const std::length_error &) {
    return {WVKernelStatusCode::sizeOverflow,
            "No-motion profile storage size overflowed."};
  }
}

bool WVNoMotionProfile::containsHeight(double height) const noexcept {
  return heights_.size() >= 2 && std::isfinite(height) &&
         height >= heights_.front() && height <= heights_.back();
}

std::size_t WVNoMotionProfile::heightInterval(double height) const noexcept {
  const auto upper = std::upper_bound(heights_.begin(), heights_.end(), height);
  return std::min(static_cast<std::size_t>(upper - heights_.begin() - 1),
                  coefficients_.size() - 1);
}

WVKernelStatus WVNoMotionProfile::density(double height, double &value) const {
  if (!containsHeight(height))
    return invalid("Density query height is outside the initialized profile.");
  const auto interval = heightInterval(height);
  const double result = densityOffset_ +
                        densityScale_ * polynomial(coefficients_[interval],
                                                   height - heights_[interval]);
  if (!std::isfinite(result))
    return numerical("No-motion density evaluation overflowed.");
  value = result;
  return WVKernelStatus::ok();
}

WVKernelStatus WVNoMotionProfile::inverse(double value, double &height) const {
  if (heights_.size() < 2 || !std::isfinite(value))
    return invalid("Density inversion requires a finite value and initialized profile.");
  // Compare distances rather than expanded endpoints: adding tolerance to a
  // near-maximum finite density must not turn its bound into infinity.
  if ((value < densities_.back() && densities_.back() - value > densityTolerance_) ||
      (value > densities_.front() && value - densities_.front() > densityTolerance_))
    return invalid("Density is outside the no-motion profile range.");
  const double target =
      (std::clamp(value, densities_.back(), densities_.front()) - densityOffset_) /
      densityScale_;
  const auto upperKnot = std::upper_bound(normalizedDensity_.begin(),
                                         normalizedDensity_.end(), target,
                                         std::greater<double>());
  const auto interval =
      std::min(static_cast<std::size_t>(upperKnot - normalizedDensity_.begin() - 1),
               coefficients_.size() - 1);
  const double base = heights_[interval];
  const auto &c = coefficients_[interval];
  double lower = 0.0;
  double upper = heights_[interval + 1] - base;
  double local = upper * (normalizedDensity_[interval] - target) /
                 (normalizedDensity_[interval] - normalizedDensity_[interval + 1]);
  for (unsigned iteration = 0; iteration < 64; ++iteration) {
    const double residual = polynomial(c, local) - target;
    if (!std::isfinite(residual))
      return numerical("No-motion density inversion arithmetic overflowed.");
    if (residual > 0.0)
      lower = local;
    if (residual < 0.0)
      upper = local;
    if (residual == 0.0 || upper - lower <= heightTolerance_) {
      const double result = base + local;
      if (!std::isfinite(result))
        return numerical("No-motion inverse height overflowed.");
      height = result;
      return WVKernelStatus::ok();
    }
    const double derivative = (3.0 * c[0] * local + 2.0 * c[1]) * local + c[2];
    if (!std::isfinite(derivative))
      return numerical("No-motion inverse derivative overflowed.");
    double next = local - residual / derivative;
    if (!std::isfinite(next) || next <= lower || next >= upper)
      next = (lower + upper) / 2.0;
    if (!std::isfinite(next))
      return numerical("No-motion inverse bracket overflowed.");
    local = next;
  }
  return numerical("No-motion density inversion did not converge in 64 iterations.");
}

WVKernelStatus WVNoMotionProfile::availablePotentialEnergy(
    double height, double materialHeight, double gravity, double referenceDensity,
    double &value) const {
  if (!containsHeight(height) || !containsHeight(materialHeight))
    return invalid("APE heights must be inside the initialized no-motion profile.");
  if (!std::isfinite(gravity) || !(gravity > 0.0) ||
      !std::isfinite(referenceDensity) || !(referenceDensity > 0.0))
    return invalid("APE gravity and reference density must be positive and finite.");
  const double factor = gravity * densityScale_ / referenceDensity;
  if (!std::isfinite(factor))
    return numerical("No-motion APE normalization overflowed.");
  double low = std::min(height, materialHeight);
  const double high = std::max(height, materialHeight);
  auto interval = heightInterval(low);
  double integral = 0.0;
  while (low < high) {
    const double end = std::min(high, heights_[interval + 1]);
    const double t = low - heights_[interval];
    const double d = end - low;
    const auto &c = coefficients_[interval];
    const double a = (3.0 * c[0] * t + 2.0 * c[1]) * t + c[2];
    const double b = 6.0 * c[0] * t + 2.0 * c[1];
    const double cubic = 3.0 * c[0];
    const double displacement = height - low;
    // Integrate (r-height)*rho'(r) directly over this crossed interval.
    // Density subtraction would erase the energy of tiny displacements.
    const double piece = -displacement * a * d +
                         (a - displacement * b) * (d * d) / 2.0 +
                         (b - displacement * cubic) * (d * d * d) / 3.0 +
                         cubic * (d * d * d * d) / 4.0;
    integral += piece;
    if (!std::isfinite(integral))
      return numerical("No-motion APE integral overflowed.");
    low = end;
    ++interval;
  }
  const double sign = height > materialHeight ? 1.0 :
                      height < materialHeight ? -1.0 : 0.0;
  const double result = factor * sign * integral;
  if (!std::isfinite(result))
    return numerical("No-motion APE evaluation overflowed.");
  value = result;
  return WVKernelStatus::ok();
}

std::size_t WVNoMotionProfile::retainedBytes() const noexcept {
  return sizeof(*this) + sizeof(double) *
                             (heights_.capacity() + densities_.capacity() +
                              normalizedDensity_.capacity()) +
         sizeof(std::array<double, 4>) * coefficients_.capacity();
}

} // namespace wavevortex::runtime
