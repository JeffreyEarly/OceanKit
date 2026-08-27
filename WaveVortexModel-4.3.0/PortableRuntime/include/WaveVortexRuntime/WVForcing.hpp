#pragma once

#include "WaveVortexRuntime/WVForcingSchedule.hpp"
#include "WaveVortexRuntime/WVIntegrationContracts.hpp"

#include <cstddef>
#include <cstdint>
#include <array>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVConstantStratificationForcingEngine;
class WVConstantStratificationRightHandSideContext;

struct WVFixedAmplitudeConfiguration {
  std::vector<std::size_t> ApIndices;
  std::vector<WVComplex64> ApValues;
  std::vector<std::size_t> AmIndices;
  std::vector<WVComplex64> AmValues;
  std::vector<std::size_t> A0Indices;
  std::vector<WVComplex64> A0Values;
};

struct WVPseudoTopographicConfiguration {
  WVShape2D topographicShape;
  std::vector<double> topographicHeight;
  std::array<WVComplex64, 2> barotropicVelocityAmplitude{};
  double frequency = 0.0;
  std::string darwinSymbol;
  double rampDuration = 0.0;
  double startTime = 0.0;
  bool shouldAvoidAdaptiveDamping = true;
  double maximumForcedHorizontalWavenumber = 0.0;
  double maximumForcedVerticalMode = 0.0;
};

struct WVPseudoTopographicOperators {
  WVPseudoTopographicConfiguration configuration;
  std::vector<WVComplex64> responsePlusX;
  std::vector<WVComplex64> responsePlusY;
  std::vector<WVComplex64> responseMinusX;
  std::vector<WVComplex64> responseMinusY;
};

// Borrowed per-call services exposed to one resolved forcing. These operations
// are coarse: each call performs a complete transform or contiguous array loop.
class WVForcingExecutionContext final {
public:
  WVKernelStatus nonlinearAdvection();
  WVKernelStatus physicalFields(WVRealFieldBundleConstView &fields);
  WVRealFieldBundleView clearedSpatialTendency();
  WVKernelStatus projectSpatialTendency(
      const WVRealFieldBundleConstView &tendency);
  WVKernelStatus adaptiveDamping(const std::vector<double> &operatorValues);
  WVKernelStatus pseudoTopographicGeneration(
      const WVPseudoTopographicOperators &operators);
  WVKernelStatus betaPlaneAdvection(
      const std::vector<WVComplex64> &operatorValues);
  void zeroSelectedTendencies(
      const WVFixedAmplitudeConfiguration &configuration);
  WVKernelStatus linearCoefficientTendency(double rate);

private:
  WVConstantStratificationForcingEngine *engine_ = nullptr;
  const WVState *state_ = nullptr;
  WVFlux *flux_ = nullptr;
  bool *outputInitialized_ = nullptr;
  WVRealFieldBundleView *externalFields_ = nullptr;
  bool *externalFieldsPrepared_ = nullptr;
  friend class WVConstantStratificationForcingEngine;
};

// Stable source API v1 boundary for one exact MATLAB identity/version forcing.
// A factory-created implementation owns immutable typed configuration and
// derived operators and is retained by its resolved model. Calls occur once
// per forcing stage or constraint pass; implementations run their complete
// contiguous coefficient/grid loops below that boundary.
class WVForcing {
public:
  virtual ~WVForcing() = default;
  virtual const std::string &typeIdentifier() const noexcept = 0;
  virtual std::uint32_t contractVersion() const noexcept = 0;
  virtual const std::string &name() const noexcept = 0;
  virtual WVForcingStage stage() const noexcept = 0;
  virtual std::uint8_t priority() const noexcept = 0;
  virtual std::size_t ordinal() const noexcept = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
  virtual bool requiresPhysicalFields() const noexcept { return false; }
  virtual bool requiresForcingFields() const noexcept { return false; }
  virtual bool producesCompleteFlux() const noexcept { return false; }
  virtual std::size_t constraintWriteCount() const noexcept { return 0; }

  virtual WVKernelStatus addRightHandSide(
      WVForcingExecutionContext &context) const = 0;

  virtual WVStateConstraintResult
  applyConstraint(WVMutableCoefficients &) const {
    return {WVKernelStatus::ok(), 0, true};
  }
};

} // namespace wavevortex::runtime
