#pragma once

namespace wavevortex::runtime {

// Runtime-only numerical selection. MATLAB restart files do not encode the
// reference flag or a legacy solver choice; omission selects the corrected
// actual-density contract rather than inferring a choice from the checkpoint.
enum class WVNoMotionReference { actual, initial };

struct WVDensityDiagnosticContract {
  static constexpr const char *identifier = "wave-vortex-density-diagnostics-v1";
  WVNoMotionReference reference = WVNoMotionReference::actual;

  const char *referenceIdentifier() const noexcept {
    switch (reference) {
    case WVNoMotionReference::actual: return "actual";
    case WVNoMotionReference::initial: return "initial";
    }
    return "invalid";
  }
  const char *recoveryIdentifier() const noexcept {
    switch (reference) {
    case WVNoMotionReference::actual: return "dampedLeastSquares";
    case WVNoMotionReference::initial: return "not-required";
    }
    return "invalid";
  }
};

} // namespace wavevortex::runtime
