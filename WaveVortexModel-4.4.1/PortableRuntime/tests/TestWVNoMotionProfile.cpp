#include "WaveVortexRuntime/WVNoMotionProfile.hpp"
#include "WVAllocationProbe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double absoluteTolerance,
           double relativeTolerance, const char *message) {
  require(std::isfinite(actual) && std::abs(actual - expected) <=
              absoluteTolerance + relativeTolerance * std::abs(expected), message);
}
WVNoMotionProfile profileFor(const std::vector<double> &z,
                             const std::vector<double> &rho) {
  WVNoMotionProfile profile;
  require(bool(WVNoMotionProfile::create(z, rho, profile)), "valid profile creation failed");
  return profile;
}
double density(const WVNoMotionProfile &profile, double z) {
  double result = -719;
  require(bool(profile.density(z, result)), "valid density query failed");
  return result;
}
double inverse(const WVNoMotionProfile &profile, double rho) {
  double result = -719;
  require(bool(profile.inverse(rho, result)), "valid inverse query failed");
  return result;
}
double ape(const WVNoMotionProfile &profile, double z, double s,
           double g = 1, double rho0 = 1) {
  double result = -719;
  require(bool(profile.availablePotentialEnergy(z, s, g, rho0, result)), "valid APE query failed");
  return result;
}

// Exact PCHIP pieces for heights [-1,0,1], densities [4,3,0]. The first
// endpoint derivative is zero; these formulas do not read production state.
double sharpDensity(double z) {
  if (z <= 0) {
    const double t = z + 1;
    return 4 - 1.5 * t * t + .5 * t * t * t;
  }
  return 3 - 1.5 * z - 2 * z * z + .5 * z * z * z;
}
double sharpDerivative(double z) {
  if (z <= 0) {
    const double t = z + 1;
    return -3 * t + 1.5 * t * t;
  }
  return -1.5 - 4 * z + 1.5 * z * z;
}

// Hermite values on unequal intervals [-1,0,2]. Known PCHIP nodal slopes
// follow weighted harmonic secants: -5/6, -27/23, -11/6.
double irregularDensity(double z) {
  const bool first = z <= 0;
  const double h = first ? 1 : 2;
  const double t = first ? z + 1 : z / 2;
  const double y0 = first ? 4 : 3, y1 = first ? 3 : 0;
  const double d0 = first ? -5.0 / 6 : -27.0 / 23;
  const double d1 = first ? -27.0 / 23 : -11.0 / 6;
  return (2 * t * t * t - 3 * t * t + 1) * y0 +
         (t * t * t - 2 * t * t + t) * h * d0 +
         (-2 * t * t * t + 3 * t * t) * y1 +
         (t * t * t - t * t) * h * d1;
}

double sharpAPEByQuadrature(double z, double s) {
  if (z == s) return 0;
  const double lower = std::min(z, s), upper = std::max(z, s);
  const std::array<double, 3> nodes{-std::sqrt(3.0 / 5), 0, std::sqrt(3.0 / 5)};
  const std::array<double, 3> weights{5.0 / 9, 8.0 / 9, 5.0 / 9};
  const std::array<double, 3> cuts{lower, std::clamp(0.0, lower, upper), upper};
  double result = 0;
  for (std::size_t interval = 0; interval < 2; ++interval) {
    const double halfWidth = (cuts[interval + 1] - cuts[interval]) / 2;
    const double center = (cuts[interval + 1] + cuts[interval]) / 2;
    for (std::size_t point = 0; point < nodes.size(); ++point) {
      const double r = center + halfWidth * nodes[point];
      result += halfWidth * weights[point] * (r - z) * sharpDerivative(r);
    }
  }
  return z > s ? result : -result;
}

void verifyLinearAndTinyDisplacements() {
  constexpr double N2 = .04;
  for (const auto &knots : {std::vector<double>{-1, 1},
                            std::vector<double>{-1, -.4, 0, .4, 1}}) {
    std::vector<double> rho;
    for (const double z : knots) rho.push_back(1 - N2 * z);
    const auto profile = profileFor(knots, rho);
    for (const double z : {-1.0, -.83, -.4, 0.0, .13, .4, .89, 1.0}) {
      close(density(profile, z), 1 - N2 * z, 3e-16, 0, "linear density changed");
      close(inverse(profile, 1 - N2 * z), z, 1e-14, 0, "linear inverse changed");
    }
    const double eps = std::numeric_limits<double>::epsilon();
    const std::vector<std::array<double, 2>> pairs{
        {-1, 1}, {1, -1}, {0, 0}, {-.7, .3}, {.4, -.4},
        {-.25 + 1e-12, -.25}, {-.25 - 1e-12, -.25},
        {std::nextafter(-.25, 0.0), -.25},
        {std::nextafter(-.25, -1.0), -.25},
        {eps, 0}, {-eps, 0}, {eps, -eps}, {-eps, eps},
        {std::nextafter(.4, 1.0), .4}, {std::nextafter(-.4, -1.0), -.4}};
    for (const auto &pair : pairs) {
      const double eta = pair[0] - pair[1];
      close(ape(profile, pair[0], pair[1]), N2 * eta * eta / 2,
            0, 1e-13, "linear/tiny/cross-knot APE changed");
    }
  }
}

void verifyNonlinearAndIrregularProfiles() {
  const auto sharp = profileFor({-1, 0, 1}, {4, 3, 0});
  for (const double z : {-1.0, -1 + 1e-4, -.91, -.37, 0.0, .17, .83, 1.0}) {
    close(density(sharp, z), sharpDensity(z), 2e-15, 0, "sharp PCHIP density differs from polynomial");
    close(inverse(sharp, sharpDensity(z)), z, 2e-11, 0, "zero-slope endpoint inverse changed");
  }
  for (const double z : {-1.0, -.71, 0.0, .13, 1.0})
    for (const double s : {-1.0, -.21, 0.0, .89, 1.0}) {
      const double value = ape(sharp, z, s);
      require(value >= 0, "stable nonlinear profile produced negative APE");
      close(value, sharpAPEByQuadrature(z, s), 3e-15, 2e-13,
            "nonlinear multi-interval APE differs from independent quadrature");
    }
  const auto irregular = profileFor({-1, 0, 2}, {4, 3, 0});
  for (const double z : {-1.0, -.87, -.1, 0.0, .13, .9, 1.7, 2.0}) {
    close(density(irregular, z), irregularDensity(z), 3e-15, 0, "irregular-knot PCHIP density changed");
    close(inverse(irregular, irregularDensity(z)), z, 1e-14, 0, "irregular-knot inverse changed");
  }
}

void verifyWeakAndSteepGradients() {
  constexpr double rho0 = 1025, g = 9.81, N2 = 1e-10;
  const std::vector<double> knots{-1000, -731, -207, 0};
  std::vector<double> rho;
  for (const double z : knots) rho.push_back(rho0 * (1 - N2 * z / g));
  const auto weak = profileFor(knots, rho);
  const double densityUlp = std::nextafter(rho0, std::numeric_limits<double>::infinity()) - rho0;
  for (const double z : {-850.0, -500.0, -100.0}) {
    close(inverse(weak, rho0 * (1 - N2 * z / g)), z,
          4 * densityUlp * g / (rho0 * N2), 0, "weak gradient inverse exceeds density resolution");
    close(ape(weak, z + 2, z, g, rho0), 2 * N2, 0, 2e-6, "weak gradient APE changed");
  }
  const auto steep = profileFor({-1e-6, 0, 1e-6}, {4, 3, 0});
  for (const double scaled : {-1.0, -.8, -.1, 0.0, .37, 1.0}) {
    close(density(steep, scaled * 1e-6), sharpDensity(scaled), 3e-15, 0, "steep density changed");
    close(inverse(steep, sharpDensity(scaled)), scaled * 1e-6, 1e-20, 0, "steep gradient inverse changed");
    close(ape(steep, scaled * 1e-6, -.3e-6), sharpAPEByQuadrature(scaled, -.3) * 1e-6,
          2e-21, 2e-13, "steep gradient APE changed");
  }
}

void verifyValidationAndTransactionalOutputs() {
  WVNoMotionProfile profile;
  double output = 719;
  require(!profile.density(0, output) && output == 719, "uninitialized density query changed output");
  require(!profile.inverse(0, output) && output == 719, "uninitialized inverse query changed output");
  require(!profile.availablePotentialEnergy(0, 0, 1, 1, output) && output == 719,
          "uninitialized APE query changed output");
  profile = profileFor({-1, 0, 1}, {4, 3, 0});
  const auto retained = profile.retainedBytes();
  const double value = density(profile, .2);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<std::pair<std::vector<double>, std::vector<double>>> invalid{
      {{}, {}}, {{0}, {1}}, {{-1, 0, 1}, {3, 2}},
      {{-1, 0, 0}, {3, 2, 1}}, {{-1, 0, 1}, {3, 2, 2}},
      {{-1, 0, 1}, {1, 2, 3}}, {{-1, nan, 1}, {3, 2, 1}},
      {{-1, 0, 1}, {3, nan, 1}}};
  for (const auto &inputs : invalid) {
    require(!WVNoMotionProfile::create(inputs.first, inputs.second, profile), "invalid profile accepted");
    require(profile.knotCount() == 3 && profile.retainedBytes() == retained &&
                density(profile, .2) == value, "failed creation changed the existing profile");
  }
  const double largest = std::numeric_limits<double>::max();
  require(WVNoMotionProfile::create({-1, 0, 1}, {largest, 0, -largest}, profile).code ==
              WVKernelStatusCode::numericalFailure,
          "unrepresentable density normalization was not rejected");
  require(WVNoMotionProfile::create({-1, 0, 1, 2}, {1e308, 2e-300, 1e-300, 0}, profile).code ==
              WVKernelStatusCode::numericalFailure,
          "collapsed normalized density intervals were not rejected");
  require(profile.retainedBytes() == retained && density(profile, .2) == value,
          "failed numerical construction changed the existing profile");
  const std::vector<double> replacementHeights{-1, 1}, replacementDensities{5, 0};
  allocationProbe::failAfter = 0;
  const auto allocationStatus = WVNoMotionProfile::create(replacementHeights, replacementDensities, profile);
  allocationProbe::failAfter = -1;
  require(allocationStatus.code == WVKernelStatusCode::allocationFailure &&
              profile.retainedBytes() == retained && density(profile, .2) == value,
          "allocation failure changed the existing profile");
  for (const double z : {-1.01, 1.01, nan, std::numeric_limits<double>::infinity()}) {
    require(!profile.density(z, output) && output == 719, "invalid height changed density output");
    require(!profile.availablePotentialEnergy(z, 0, 1, 1, output) && output == 719,
            "invalid current height changed APE output");
    require(!profile.availablePotentialEnergy(0, z, 1, 1, output) && output == 719,
            "invalid material height changed APE output");
  }
  const double ulp = std::nextafter(4.0, 5.0) - 4;
  close(inverse(profile, 4 + 8 * ulp), -1, 0, 0, "upper density roundoff was not clamped");
  close(inverse(profile, -8 * ulp), 1, 0, 0, "lower density roundoff was not clamped");
  for (const double rho : {4 + 9 * ulp, -9 * ulp, nan})
    require(!profile.inverse(rho, output) && output == 719, "invalid density changed inverse output");
  for (const double parameter : {0.0, -1.0, nan}) {
    require(!profile.availablePotentialEnergy(.2, -.3, parameter, 1, output) && output == 719,
            "invalid gravity changed APE output");
    require(!profile.availablePotentialEnergy(.2, -.3, 1, parameter, output) && output == 719,
            "invalid reference density changed APE output");
  }
  require(!profile.availablePotentialEnergy(1, -1, 1e308, 1, output) && output == 719,
          "overflowing APE changed caller output");
  const double tiny = std::numeric_limits<double>::denorm_min();
  const auto subnormal = profileFor({-1, 1}, {2 * tiny, tiny});
  close(inverse(subnormal, 2 * tiny), -1, 0, 0, "subnormal upper endpoint inverse changed");
  close(inverse(subnormal, tiny), 1, 0, 0, "subnormal lower endpoint inverse changed");
  close(inverse(subnormal, 10 * tiny), -1, 0, 0, "subnormal density upper clamp changed");
  close(inverse(subnormal, -7 * tiny), 1, 0, 0, "subnormal density lower clamp changed");
  require(!subnormal.inverse(11 * tiny, output) && output == 719 &&
              !subnormal.inverse(-8 * tiny, output) && output == 719,
          "subnormal density outside roundoff band changed inverse output");
}

void verifyBoundedStorageAndAllocationFreeQueries() {
  std::vector<double> heights(65), densities(65);
  for (std::size_t index = 0; index < heights.size(); ++index) {
    heights[index] = -1 + 2 * static_cast<double>(index) / (heights.size() - 1);
    densities[index] = 3 - heights[index] - .1 * heights[index] * heights[index];
  }
  const auto profile = profileFor(heights, densities);
  const auto bytes = profile.retainedBytes();
  require(bytes >= sizeof(profile) && bytes <= sizeof(profile) + 16 * heights.size() * sizeof(double),
          "profile storage is not bounded by its knot count");
  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  bool success = true;
  double checksum = 0;
  for (std::size_t index = 0; index < 70000; ++index) {
    const double z = -.99 + 1.98 * static_cast<double>(index % 997) / 996;
    double rho = 0, s = 0, value = 0;
    success = bool(profile.density(z, rho)) && success;
    success = bool(profile.inverse(rho, s)) && success;
    success = bool(profile.availablePotentialEnergy(-z, s, 9.81, 1025, value)) && success;
    checksum += value;
  }
  allocationProbe::counting = false;
  require(success && std::isfinite(checksum) && checksum > 0, "repeated scalar queries failed");
  require(allocationProbe::calls == 0, "successful profile queries allocated heap storage");
  require(profile.retainedBytes() == bytes && profile.knotCount() == heights.size(),
          "query history changed immutable profile storage");
}
}

int main() {
  try {
    verifyLinearAndTinyDisplacements();
    verifyNonlinearAndIrregularProfiles();
    verifyWeakAndSteepGradients();
    verifyValidationAndTransactionalOutputs();
    verifyBoundedStorageAndAllocationFreeQueries();
    std::cout << "WVNoMotionProfile contracts passed: analytic density/inverse, exact knots, "
                 "weak/steep gradients, sub-ulp APE, quadrature, validation, immutable storage, zero query allocations.\n";
    return 0;
  } catch (const std::exception &error) {
    allocationProbe::counting = false;
    std::cerr << error.what() << '\n';
    return 1;
  }
}
