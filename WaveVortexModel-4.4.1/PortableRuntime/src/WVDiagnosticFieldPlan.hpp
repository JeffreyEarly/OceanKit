#pragma once
#include <array>

#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVPortableVariablePlan.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WaveVortexRuntime/WVDensityDiagnosticContract.hpp"

namespace wavevortex::runtime::detail {

// Immutable output/dependency plan. All coefficient copies and intermediate
// fields are owned by one evaluate() invocation and released before it returns.
class WVDiagnosticFieldPlan final {
public:
  static std::string configurationIdentifier(const WVFieldEvaluationService&);
  static bool required(const std::vector<WVFieldRequest>&, bool stratified = false) noexcept;
  static WVKernelStatus create(const WVFieldEvaluationService&,
      const std::vector<WVFieldRequest>&, WVFieldEvaluationPlan&, WVDensityDiagnosticContract = {},
      bool prepareScientificDependencies = true);
  // Private implementation qualification seam. Public field-plan construction
  // continues to enforce the catalog's unavailable density-output contracts.
  static WVKernelStatus createDensityQualification(const WVFieldEvaluationService&,
      const std::vector<WVFieldRequest>&, WVDensityDiagnosticContract,
      WVFieldEvaluationPlan&);
  static WVKernelStatus createDensityQualificationForActiveEvaluation(
      const WVFieldEvaluationService&,const std::vector<WVFieldRequest>&,
      WVDensityDiagnosticContract,WVFieldEvaluationPlan&);
  WVKernelStatus rebind(const WVFieldEvaluationService&, WVFieldEvaluationPlan&) const;
  WVKernelStatus evaluate(WVFieldEvaluationService&, const WVIntegrationState&,
      WVFieldOutputView*, std::size_t, const std::uint8_t* activeOutputs = nullptr) const;
  std::size_t persistentBytes() const noexcept;
  WVKernelStatus prepareEventArena(const WVFieldEvaluationService&) const;
  WVKernelStatus prepareEventArenaForActiveEvaluation(
      const WVFieldEvaluationService&) const;
  bool hasForcingDiagnostics() const noexcept {return !forcingIndices_.empty();}
  bool hasDensityDiagnostics() const noexcept {return hasDensity_;}

private:
  struct Output {
    WVPortableVariable variable = WVPortableVariable::invalid;
    WVPortableVariablePlan execution;
    std::size_t group = 0, dependency = 0;
    bool surface = false, extrema = false, verticalMean = false, forcing = false, density = false;
    bool sampled = false;
    std::size_t forcingSlot=0,forcingChannel=0,forcingPhysicalChannels=0;
    std::array<std::size_t,4> auxiliaries{};
    WVFieldOutputSpecification specification;
    WVFieldEvaluationPlan sampler;
    WVFieldSamplingRequest sampling;
  };
  struct Group {
    WVFieldEvaluationPlan fields;
    std::vector<WVFieldRequest> requests;
  };
  WVKernelStatus configure(const WVFieldEvaluationService&);
  WVKernelStatus prepareEventArenaImpl(
      const WVFieldEvaluationService&,bool active) const;
  static WVKernelStatus createImpl(const WVFieldEvaluationService&,
      const std::vector<WVFieldRequest>&, WVFieldEvaluationPlan&,
      bool densityQualification, WVDensityDiagnosticContract,
      bool prepareScientificDependencies);
  double omega(std::size_t) const noexcept;
  bool keep(std::size_t group, std::size_t family, std::size_t coefficient) const noexcept;
  double weight(std::size_t z) const noexcept;
  double stratification(std::size_t z) const noexcept;
  const WVFieldEvaluationService* owner_ = nullptr;
  const WVTransformConstantStratificationKernel* constant_ = nullptr;
  const WVTransformHydrostaticKernel* hydrostatic_ = nullptr;
  const WVTransformBoussinesqKernel* boussinesq_ = nullptr;
  const WVStratifiedModalGeometry* modal_ = nullptr;
  std::string configuration_;
  WVShape3D spatial_{};
  WVShape2D spectral_{};
  bool isHydrostatic_ = true, isBarotropic_ = false, isQG_ = false;
  double Lz_ = 0, constantN2_ = 0, barotropicG_ = 0;
  bool densityQualification_ = false, hasDensity_ = false;
  std::size_t densityDependency_=0;
  WVDensityDiagnosticContract densityContract_;
  std::vector<double> densityHeights_, densityWeights_, densityInitial_;
  double densityGravity_ = 0, densityReference_ = 0;
  std::array<Group,5> groups_;
  std::vector<Output> outputs_;
  std::vector<std::size_t> forcingIndices_;
  std::size_t forcingPhysicalChannels_=0;
  std::uint64_t scratchSignature_=1469598103934665603ULL;
  std::array<std::size_t,4> forcingPhysicalDependencies_{};
};

} // namespace wavevortex::runtime::detail
