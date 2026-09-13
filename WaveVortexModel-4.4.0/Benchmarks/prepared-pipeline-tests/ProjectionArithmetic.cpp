#include "WVCoefficientFormulas.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using wavevortex::WVComplex64;
using wavevortex::detail::add;
using wavevortex::detail::conjugate;
using wavevortex::detail::multiply;
using wavevortex::detail::subtract;

#if defined(__clang__) || defined(__GNUC__)
#define WV_NOINLINE __attribute__((noinline))
#else
#define WV_NOINLINE
#endif

namespace {

struct Inputs {
  std::size_t Nj = 0;
  std::vector<double> k, l, ApmN;
  std::vector<unsigned char> inertial;
  std::vector<WVComplex64> divergence, verticalVelocity, density, u, v, phase;
};

struct Outputs {
  std::vector<WVComplex64> Ap, Am;
};

WVComplex64 scale(WVComplex64 value, double factor) noexcept {
  return multiply(value, factor);
}

WV_NOINLINE void legacyStoredDivergence(const Inputs& input, bool hasW,
                                        WVComplex64* stored) {
  const auto count = input.divergence.size();
  for (std::size_t i = 0; i < count; ++i) {
    if (!hasW) {
      stored[i] = input.divergence[i];
      continue;
    }
    const auto mode = i / input.Nj;
    const double halfK = std::hypot(input.k[mode], input.l[mode]) / 2;
    stored[i] = add(input.divergence[i],
                    multiply(input.verticalVelocity[i], {0, halfK}));
  }
}

WV_NOINLINE void legacyFinalProjection(const Inputs& input,
                                       const WVComplex64* stored,
                                       Outputs& output) {
  const auto count = input.divergence.size();
  for (std::size_t i = 0; i < count; ++i) {
    const auto n = scale(input.density[i], input.ApmN[i]);
    auto ap = add(stored[i], n);
    auto am = subtract(stored[i], n);
    if (input.inertial[i]) {
      ap = scale(subtract(input.u[i], multiply(input.v[i], {0, 1})), .5);
      am = conjugate(ap);
    }
    output.Ap[i] = multiply(ap, conjugate(input.phase[i]));
    output.Am[i] = multiply(am, input.phase[i]);
  }
}

WV_NOINLINE void legacyProjection(const Inputs& input, bool hasW,
                                  Outputs& output) {
  std::vector<WVComplex64> stored(input.divergence.size());
  legacyStoredDivergence(input, hasW, stored.data());
  legacyFinalProjection(input, stored.data(), output);
}

WV_NOINLINE void fusedProjection(const Inputs& input, bool hasW,
                                 Outputs& output) {
  const auto modeCount = input.k.size();
  for (std::size_t mode = 0; mode < modeCount; ++mode) {
    const double halfK =
        hasW ? std::hypot(input.k[mode], input.l[mode]) / 2 : 0;
    for (std::size_t j = 0; j < input.Nj; ++j) {
      const auto i = j + input.Nj * mode;
      const auto n = scale(input.density[i], input.ApmN[i]);
      const auto d = hasW
                         ? add(input.divergence[i],
                               multiply(input.verticalVelocity[i], {0, halfK}))
                         : input.divergence[i];
      auto ap = add(d, n);
      auto am = subtract(d, n);
      if (input.inertial[i]) {
        ap = scale(subtract(input.u[i], multiply(input.v[i], {0, 1})), .5);
        am = conjugate(ap);
      }
      output.Ap[i] = multiply(ap, conjugate(input.phase[i]));
      output.Am[i] = multiply(am, input.phase[i]);
    }
  }
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result = 0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void requireSame(const std::vector<WVComplex64>& expected,
                 const std::vector<WVComplex64>& actual,
                 const char* channel, bool hasW) {
  if (expected.size() != actual.size())
    throw std::runtime_error("Projection output size changed");
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (bits(expected[i].real) == bits(actual[i].real) &&
        bits(expected[i].imag) == bits(actual[i].imag))
      continue;
    std::cerr << std::hex << std::showbase << "hasW=" << hasW
              << " channel=" << channel << " index=" << std::dec << i
              << " expected=(" << std::hex << bits(expected[i].real) << ','
              << bits(expected[i].imag) << ") actual=("
              << bits(actual[i].real) << ',' << bits(actual[i].imag) << ")\n";
    throw std::runtime_error("Fused projection changed stored-pass arithmetic");
  }
}

Inputs fixture() {
  Inputs input;
  input.Nj = 11;
  input.k = {0.0, 1.0, -2.0, 0x1.8p-20, -3.25, 9.0, -0.0};
  input.l = {-0.0, 0.0, 0.75, -0x1.2p-18, -4.5, 12.0, 2.0};
  const auto count = input.Nj * input.k.size();
  input.ApmN.resize(count);
  input.inertial.resize(count);
  input.divergence.resize(count);
  input.verticalVelocity.resize(count);
  input.density.resize(count);
  input.u.resize(count);
  input.v.resize(count);
  input.phase.resize(count);

  for (std::size_t i = 0; i < count; ++i) {
    const auto mode = i / input.Nj;
    const double sign = i % 2 == 0 ? 1.0 : -1.0;
    input.ApmN[i] = sign * (0x1.123456789abcp-7 + i * 0x1p-18);
    input.inertial[i] = static_cast<unsigned char>(mode == 0);
    input.divergence[i] = {
        sign * (0x1.3456789abcdep+4 + i * 0x1p-12),
        -sign * (0x1.bcdef01234567p-3 + i * 0x1p-16)};
    input.verticalVelocity[i] = {
        -sign * (0x1.76543210fedcbp-4 + i * 0x1p-17),
        sign * (0x1.23456789abcdep+3 + i * 0x1p-13)};
    input.density[i] = {
        sign * (0x1.0fedcba987654p+2 + i * 0x1p-15),
        sign * (0x1.89abcdef01234p-5 + i * 0x1p-19)};
    input.u[i] = {sign * (1.0 + i * 0x1p-10),
                  -sign * (0.25 + i * 0x1p-12)};
    input.v[i] = {-sign * (0.5 + i * 0x1p-11),
                  sign * (0.75 + i * 0x1p-13)};
    const double angle = 0.013 * static_cast<double>(i + 1);
    input.phase[i] = {std::cos(angle), std::sin(angle)};
  }

  // Signed-zero modes exercise the exact complex helper behavior. The first
  // mode is also the inertial branch, whose result must remain independent of
  // the discarded divergence and density values.
  input.divergence[0] = {0.0, -0.0};
  input.verticalVelocity[0] = {-0.0, 0.0};
  input.density[0] = {0.0, -0.0};
  input.u[0] = {-0.0, 0.0};
  input.v[0] = {0.0, -0.0};
  input.phase[0] = {1.0, -0.0};

  const std::size_t signedZero = input.Nj * 6 + 3;
  input.divergence[signedZero] = {-0.0, 0.0};
  input.verticalVelocity[signedZero] = {0.0, -0.0};
  input.density[signedZero] = {-0.0, 0.0};
  input.ApmN[signedZero] = -0.0;
  input.phase[signedZero] = {1.0, 0.0};

  // A large finite cancellation leaves a small stored residual before the
  // density product is added. It detects contraction across that old pass.
  const std::size_t cancellation = input.Nj + 5;
  const double product = std::ldexp(1.0, 499);
  input.verticalVelocity[cancellation] = {0.0, std::ldexp(1.0, 500)};
  input.divergence[cancellation] = {
      std::nextafter(product, std::numeric_limits<double>::infinity()),
      -std::ldexp(1.0, -500)};
  input.density[cancellation] = {std::ldexp(1.0, 500),
                                 -std::ldexp(1.0, 498)};
  input.ApmN[cancellation] = std::ldexp(1.0, -500);
  input.phase[cancellation] = {0x1.fffffffffffffp-1, 0x1p-27};

  // Keep extreme finite operands in a non-inertial non-cancelling lane too.
  const std::size_t extreme = 4 * input.Nj + 7;
  input.divergence[extreme] = {std::ldexp(1.0, 700),
                               -std::ldexp(1.0, 699)};
  input.verticalVelocity[extreme] = {std::ldexp(1.0, 690),
                                     std::ldexp(1.0, 688)};
  input.density[extreme] = {-std::ldexp(1.0, 710),
                            std::ldexp(1.0, 709)};
  input.ApmN[extreme] = std::ldexp(1.0, -40);
  input.phase[extreme] = {0.75, -0.5};
  return input;
}

void run(bool hasW) {
  const auto input = fixture();
  Outputs expected{std::vector<WVComplex64>(input.divergence.size()),
                   std::vector<WVComplex64>(input.divergence.size())};
  Outputs actual{std::vector<WVComplex64>(input.divergence.size()),
                 std::vector<WVComplex64>(input.divergence.size())};
  legacyProjection(input, hasW, expected);
  fusedProjection(input, hasW, actual);
  requireSame(expected.Ap, actual.Ap, "Ap", hasW);
  requireSame(expected.Am, actual.Am, "Am", hasW);
}

} // namespace

int main() try {
  run(false);
  run(true);
  std::cout << "Projection arithmetic preserves legacy stored-pass bits\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}

#undef WV_NOINLINE
