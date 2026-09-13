#pragma once

#include "WaveVortexRuntime/generated/WVPortableVariableContracts.hpp"
#include <string>

namespace wavevortex::runtime {

inline bool isPortableForcingVariableName(std::string_view name) noexcept {
  return name.rfind("Fu_",0)==0 || name.rfind("Fv_",0)==0 || name.rfind("Fw_",0)==0 ||
      name.rfind("Feta_",0)==0 || name.rfind("Fqgpv_",0)==0;
}

inline std::string portableVariableConfigurationIdentifier(std::string_view transform,
    bool isHydrostatic,bool shouldAntialias) {
  std::string prefix;
  if(transform=="WVTransformConstantStratification") prefix=isHydrostatic ? "constant-hydrostatic" : "constant-nonhydrostatic";
  else if(transform=="WVTransformHydrostatic") prefix="hydrostatic";
  else if(transform=="WVTransformBoussinesq") prefix="boussinesq";
  else if(transform=="WVTransformBarotropicQG") prefix="barotropic";
  else if(transform=="WVTransformStratifiedQG") prefix="stratified-qg";
  else return {};
  return prefix+(shouldAntialias ? "-aa1" : "-aa0");
}

// Construction-time contracts only. This plan neither allocates fields nor
// evaluates operations. Numerical consumers retain their existing evaluators.
enum class WVPortableVariableStatus : std::uint8_t {
  supported, unknownVariable, unknownConfiguration, notApplicable,
  unsupportedSampling, customOperation, intentionalIncompatibility,
  requiresForcingBinding, requiresNoMotionSolver, pendingEvaluation,
  invalidContract, ambiguousForcingBinding
};

enum class WVPortableNoMotionSolver : std::uint8_t {
  unspecified, lsqnonlin, fminsearch, dampedLeastSquares
};

enum class WVPortableOperationSource : std::uint8_t {
  unverified, builtIn, custom
};

struct WVPortableVariableOptions {
  // Set builtIn only after validating the operation identity/provenance during
  // graph ingestion. Matching a user-supplied name does not establish identity.
  WVPortableOperationSource source = WVPortableOperationSource::unverified;
  bool shouldUseTrueNoMotionProfile = true;
  WVPortableNoMotionSolver noMotionSolver = WVPortableNoMotionSolver::unspecified;
  // The existing forcing registry must qualify identity, version, configuration,
  // stage and unique sanitized instance name before this can be true.
  bool hasQualifiedForcingBinding = false;
  bool requireEvaluator = false;
};

struct WVPortableVariablePlan {
  const WVPortableVariableContract *output = nullptr;
  std::string_view resolvedName;
  std::uint32_t forcingInstanceOrdinal = 0;
  std::array<WVPortableVariable, WVPortableVariableCatalog.size()> order{};
  std::size_t count = 0;
  std::uint64_t primitiveMask = 0;
  std::uint64_t intermediateMask = 0;
  WVPortableNoMotionSolver noMotionSolver = WVPortableNoMotionSolver::unspecified;
};

inline WVPortableVariableStatus resolvePortableVariablePlan(
    std::string_view name, std::string_view configuration,
    std::uint8_t sampling, const WVPortableVariableOptions &options,
    WVPortableVariablePlan &output) noexcept {
  output = {};
  if (options.source != WVPortableOperationSource::builtIn)
    return WVPortableVariableStatus::customOperation;
  if (options.noMotionSolver != WVPortableNoMotionSolver::unspecified &&
      options.noMotionSolver != WVPortableNoMotionSolver::lsqnonlin &&
      options.noMotionSolver != WVPortableNoMotionSolver::fminsearch &&
      options.noMotionSolver != WVPortableNoMotionSolver::dampedLeastSquares)
    return WVPortableVariableStatus::invalidContract;
  bool knownConfiguration = false;
  for (const auto candidate : WVPortableVariableConfigurations)
    knownConfiguration |= configuration == candidate;
  if (!knownConfiguration)
    return WVPortableVariableStatus::unknownConfiguration;
  for (const auto &entry : WVPortableVariableExclusions)
    if (name == entry.name)
      return WVPortableVariableStatus::intentionalIncompatibility;
  const auto *metadata = findPortableVariable(name);
  if (!metadata) return WVPortableVariableStatus::unknownVariable;
  const auto *contract = portableVariableContract(metadata->identifier, configuration);
  if (!contract) return WVPortableVariableStatus::notApplicable;
  if (sampling == 0 || (sampling & (sampling - 1)) != 0 ||
      (contract->metadata.samplingMask & sampling) == 0)
    return WVPortableVariableStatus::unsupportedSampling;
  if (std::string_view(contract->authority) == "forcing-instance-template" &&
      !options.hasQualifiedForcingBinding)
    return WVPortableVariableStatus::requiresForcingBinding;
  if (options.requireEvaluator && contract->runtime == WVPortableDiagnosticRuntime::intentionalIncompatibility)
    return WVPortableVariableStatus::intentionalIncompatibility;
  if (options.requireEvaluator &&
      contract->runtime != WVPortableDiagnosticRuntime::implemented)
    return WVPortableVariableStatus::pendingEvaluation;

  WVPortableVariablePlan candidate;
  candidate.output = contract;
  candidate.resolvedName = name;
  candidate.noMotionSolver = options.noMotionSolver;
  std::array<std::uint8_t, WVPortableVariableCatalog.size()> state{};
  bool requiresSolver = false;
  const auto visit = [&](auto &&self, WVPortableVariable id) noexcept -> bool {
    const auto ordinal = static_cast<std::size_t>(id);
    if (ordinal >= state.size() || state[ordinal] == 1) return false;
    if (state[ordinal] == 2) return true;
    const auto *row = portableVariableContract(id, configuration);
    if (!row || row->dependencyCount > row->dependencies.size()) return false;
    state[ordinal] = 1;
    for (std::size_t index = 0; index < row->dependencyCount; ++index)
      if (!self(self, row->dependencies[index])) return false;
    if (options.shouldUseTrueNoMotionProfile &&
        row->trueProfileDependency != WVPortableVariable::invalid &&
        !self(self, row->trueProfileDependency)) return false;
    requiresSolver |= id == WVPortableVariable::rho_nm;
    candidate.primitiveMask |= row->metadata.primitiveDependencyMask;
    candidate.intermediateMask |= row->intermediateMask;
    candidate.order[candidate.count++] = id;
    state[ordinal] = 2;
    return true;
  };
  if (!visit(visit, metadata->identifier))
    return WVPortableVariableStatus::invalidContract;
  if (requiresSolver && options.noMotionSolver == WVPortableNoMotionSolver::unspecified)
    return WVPortableVariableStatus::requiresNoMotionSolver;
  if (requiresSolver && options.requireEvaluator &&
      options.noMotionSolver != WVPortableNoMotionSolver::dampedLeastSquares)
    return WVPortableVariableStatus::invalidContract;
  output = candidate;
  return WVPortableVariableStatus::supported;
}

// A binding refers to a forcing graph node already qualified by the existing
// forcing registry. Views remain valid for the lifetime of the observation plan.
struct WVPortableForcingVariableBinding {
  std::string_view instanceName;
  std::uint32_t instanceOrdinal;
  bool qualified;
};

inline WVPortableVariableStatus resolvePortableForcingVariablePlan(
    std::string_view outputName, std::string_view configuration,
    std::uint8_t sampling, const WVPortableForcingVariableBinding *bindings,
    std::size_t bindingCount, const WVPortableVariableOptions &options,
    WVPortableVariablePlan &output) noexcept {
  output = {};
  if (!bindings && bindingCount != 0)
    return WVPortableVariableStatus::invalidContract;
  constexpr std::string_view exemplar = "portable_catalog_forcing";
  const WVPortableVariableContract *match = nullptr;
  const WVPortableForcingVariableBinding *binding = nullptr;
  for (const auto &row : WVPortableVariableContracts) {
    if (configuration != row.configuration ||
        std::string_view(row.authority) != "forcing-instance-template") continue;
    const std::string_view templateName = row.metadata.name;
    const auto prefix = templateName.substr(0, templateName.size() - exemplar.size());
    if (outputName.substr(0, prefix.size()) != prefix) continue;
    for (std::size_t index = 0; index < bindingCount; ++index) {
      const auto &candidate = bindings[index];
      if (candidate.instanceName.empty() ||
          outputName.size() != prefix.size() + candidate.instanceName.size()) continue;
      bool same = true;
      for (std::size_t character = 0; character < candidate.instanceName.size(); ++character) {
        auto value = candidate.instanceName[character];
        if (value == ' ' || value == '-') value = '_';
        same &= outputName[prefix.size() + character] == value;
      }
      if (!same) continue;
      if (match) return WVPortableVariableStatus::ambiguousForcingBinding;
      match = &row;
      binding = &candidate;
    }
  }
  if (!match || !binding || !binding->qualified)
    return WVPortableVariableStatus::requiresForcingBinding;
  auto boundOptions = options;
  boundOptions.hasQualifiedForcingBinding = true;
  const auto status = resolvePortableVariablePlan(match->metadata.name, configuration,
      sampling, boundOptions, output);
  if (status == WVPortableVariableStatus::supported) {
    output.resolvedName = outputName;
    output.forcingInstanceOrdinal = binding->instanceOrdinal;
  }
  return status;
}

} // namespace wavevortex::runtime
