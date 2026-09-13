#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace wavevortex::runtime {

// Immutable shape-preserving interpolant of an explicitly supplied no-motion
// density profile. This primitive does not recover a profile from model state.
class WVNoMotionProfile final {
public:
  // Heights must increase strictly and densities decrease strictly. Creation
  // owns the inputs and derived cubic coefficients; failure leaves output
  // unchanged. Finite inputs whose normalization or interpolation arithmetic
  // is not representable are rejected with numericalFailure.
  // Optional workspace reporting records the peak candidate and temporary
  // vector capacities, including partial allocation before failure. It excludes
  // borrowed inputs, the previous output object, and allocator/object metadata.
  static WVKernelStatus create(const std::vector<double> &heights,
                               const std::vector<double> &densities,
                               WVNoMotionProfile &output,
                               std::size_t *creationWorkspaceBytes = nullptr);

  // Queries allocate no workspace and leave their output unchanged on failure.
  // Heights are restricted to the closed profile domain. Inverse density may
  // exceed an endpoint by 8*eps(max(abs(densities))); it is then clamped, matching
  // MATLAB WVNoMotionProfile. An uninitialized profile rejects all queries.
  WVKernelStatus density(double height, double &value) const;
  WVKernelStatus inverse(double value, double &height) const;
  WVKernelStatus availablePotentialEnergy(double height, double materialHeight,
                                          double gravity, double referenceDensity,
                                          double &value) const;

  std::size_t knotCount() const noexcept { return heights_.size(); }
  // Includes the object and owned vector capacities, excluding allocator metadata.
  std::size_t retainedBytes() const noexcept;

private:
  std::size_t heightInterval(double height) const noexcept;
  bool containsHeight(double height) const noexcept;

  std::vector<double> heights_;
  std::vector<double> densities_;
  std::vector<double> normalizedDensity_;
  // Descending powers of height relative to each interval's lower endpoint.
  std::vector<std::array<double, 4>> coefficients_;
  double densityOffset_ = 0.0;
  double densityScale_ = 0.0;
  double densityTolerance_ = 0.0;
  double heightTolerance_ = 0.0;
};

} // namespace wavevortex::runtime
