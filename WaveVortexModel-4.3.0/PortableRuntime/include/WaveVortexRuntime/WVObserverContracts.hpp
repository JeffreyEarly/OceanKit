#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"
#include "WaveVortexRuntime/WVObservingSystem.hpp"
#include "WaveVortexRuntime/WVPortableImplementationContract.hpp"
#include "WaveVortexRuntime/WVPortableTypedRecord.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVExtensionCatalog;

inline constexpr std::uint32_t WVPortableObserverContractVersion = 1;
inline constexpr const char *WVPortableObserverContractIdentifier =
    "portable-observers-v1";

enum class WVStateScalarType : std::uint8_t { real64, complex64 };
enum class WVToleranceKind : std::uint8_t {
  coefficientEnergyScaled,
  uniformAbsolute
};
enum class WVRestartRequirement : std::uint8_t {
  requiredDynamicState,
  derivedState
};
enum class WVStateOwnership : std::uint8_t { integratorOwned, observerDerived };
enum class WVPositionInterpolation : std::uint8_t { linear, spline };

struct WVStateBlockRecord {
  std::string identifier;
  WVStateScalarType scalarType = WVStateScalarType::real64;
  std::vector<std::size_t> dimensions;
  WVToleranceKind toleranceKind = WVToleranceKind::uniformAbsolute;
  double absoluteTolerance = 0.0;
  WVStateOwnership ownership = WVStateOwnership::integratorOwned;
  WVRestartRequirement restartRequirement =
      WVRestartRequirement::requiredDynamicState;
};

struct WVObserverRecord {
  std::string identifier;
  std::string name;
  std::string typeIdentifier;
  std::uint32_t contractVersion = WVPortablePairContractVersion;
  // Construction-only, per-record configuration resolved before integration.
  // Legacy fixed-schema readers may infer this record from the compatibility
  // fields below; new observer providers should populate it directly.
  WVPortableTypedRecord configuration;
  std::vector<std::string> stateBlockIdentifiers;
  std::vector<std::string> fieldNames;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  bool isXYOnly = false;
  bool shouldAntialias = true;
  WVPositionInterpolation advectionInterpolation =
      WVPositionInterpolation::linear;
  WVPositionInterpolation trackedFieldInterpolation =
      WVPositionInterpolation::linear;
  double horizontalAbsoluteTolerance = 0.0;
  double verticalAbsoluteTolerance = 0.0;
  double outputScale = 1.0;
  double outputOffset = 0.0;
};

class WVResolvedObserver final {
public:
  const std::string &identifier() const noexcept { return identifier_; }
  const WVPortableTypedRecord &configuration() const noexcept {
    return configuration_;
  }
  const WVObservingSystem &implementation() const noexcept {
    return *implementation_;
  }
  std::shared_ptr<const WVObservingSystem>
  implementationHandle() const noexcept {
    return implementation_;
  }
  const WVObserverExecutionPlan &executionPlan() const noexcept {
    return executionPlan_;
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &observer,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const;
  WVKernelStatus observationBatch(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &batch) const;
  WVKernelStatus prepareOccurrence(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOccurrencePreparationContext &context,
      WVObserverOccurrenceWorkspace &workspace) const;
  std::size_t persistentBytes() const noexcept;

private:
  friend class WVPortableObserverDescriptor;
  std::string identifier_;
  WVPortableTypedRecord configuration_;
  WVObserverExecutionPlan executionPlan_;
  std::shared_ptr<const WVObservingSystem> implementation_;
};

struct WVOutputScheduleRecord {
  WVOutputScheduleRecord() = default;
  WVOutputScheduleRecord(double interval, double initial, double final)
      : outputInterval(interval), initialTime(initial), finalTime(final) {}
  double outputInterval = 0.0;
  double initialTime = 0.0;
  double finalTime = 0.0;
  // Empty identifies the legacy evenly-spaced schedule. New source-linked
  // providers use an exact identity/version and a construction-only record.
  std::string typeIdentifier;
  std::uint32_t contractVersion = 1;
  WVPortableTypedRecord configuration;
};

struct WVOutputGroupRecord {
  std::string identifier;
  std::string name;
  WVOutputScheduleRecord schedule;
  std::vector<std::string> observerIdentifiers;
  bool containsCompleteCoefficientRestart = false;
};

struct WVOutputFileRecord {
  std::string identifier;
  std::string destination;
  std::vector<WVOutputGroupRecord> groups;
};

struct WVPortableObserverRecord {
  std::string schemaIdentifier = WVPortableObserverContractIdentifier;
  std::uint32_t schemaVersion = WVPortableObserverContractVersion;
  std::vector<WVStateBlockRecord> stateBlocks;
  std::vector<WVObserverRecord> observers;
  std::vector<WVOutputFileRecord> outputFiles;
};

class WVPortableObserverDescriptor final {
public:
  static WVKernelStatus create(const WVPortableObserverRecord &record,
                               std::shared_ptr<const WVExtensionCatalog> catalog,
                               WVPortableObserverDescriptor &descriptor);

  const WVPortableObserverRecord &record() const noexcept;
  const std::vector<WVStateBlockRecord> &stateBlocks() const noexcept;
  const std::vector<WVObserverRecord> &observers() const noexcept;
  const std::vector<WVOutputFileRecord> &outputFiles() const noexcept;
  const WVObservingSystem *
  implementation(const WVObserverRecord &observer) const noexcept;
  const WVResolvedObserver *
  resolvedObserver(const WVObserverRecord &observer) const noexcept;
  const std::shared_ptr<const WVExtensionCatalog> &catalog() const noexcept;
  std::size_t persistentBytes() const noexcept;

private:
  class Impl;
  std::shared_ptr<const Impl> impl_;
};

} // namespace wavevortex::runtime
