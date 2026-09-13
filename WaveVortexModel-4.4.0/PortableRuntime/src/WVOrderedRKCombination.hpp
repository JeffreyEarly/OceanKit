#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <array>
#include <cstddef>

namespace wavevortex::runtime::rk_detail {

template <std::size_t N>
inline void orderedWeightedAffine(
    double *destination, const double *base,
    const std::array<const double *, N> &inputs,
    const std::array<double, N> &weights, std::size_t count,
    double step) noexcept {
  static_assert(N > 0, "An RK combination needs at least one derivative.");
  // The destination is private RK stage storage, disjoint from the base and
  // derivative spans. This first pass preserves setScaled's stored rounding
  // boundary before the ordered addScaled operations are fused below.
  for (std::size_t index = 0; index < count; ++index)
    destination[index] = weights[0] * inputs[0][index];
  for (std::size_t index = 0; index < count; ++index) {
    double sum = destination[index];
    for (std::size_t term = 1; term < N; ++term)
      sum = sum + weights[term] * inputs[term][index];
    destination[index] = base[index] + step * sum;
  }
}

template <std::size_t N>
inline void orderedWeightedAffine(
    WVComplex64 *destination, const WVComplex64 *base,
    const std::array<const WVComplex64 *, N> &inputs,
    const std::array<double, N> &weights, std::size_t count,
    double step) noexcept {
  static_assert(N > 0, "An RK combination needs at least one derivative.");
  // See the real overload for the disjoint-span and rounding contract.
  for (std::size_t index = 0; index < count; ++index) {
    destination[index].real = weights[0] * inputs[0][index].real;
    destination[index].imag = weights[0] * inputs[0][index].imag;
  }
  for (std::size_t index = 0; index < count; ++index) {
    WVComplex64 sum = destination[index];
    for (std::size_t term = 1; term < N; ++term) {
      sum.real = sum.real + weights[term] * inputs[term][index].real;
      sum.imag = sum.imag + weights[term] * inputs[term][index].imag;
    }
    destination[index] = {base[index].real + step * sum.real,
                          base[index].imag + step * sum.imag};
  }
}

} // namespace wavevortex::runtime::rk_detail
