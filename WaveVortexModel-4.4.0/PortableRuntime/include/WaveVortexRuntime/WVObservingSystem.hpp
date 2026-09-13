#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace wavevortex::runtime {

struct WVObserverRecord;
struct WVStateBlockRecord;
struct WVObserverOutputPlanningContext;
struct WVObserverOutputPlan;
struct WVObserverOccurrencePreparationContext;
struct WVObserverOccurrenceWorkspace;
class WVObserverOutputEvaluationContext;
struct WVObservationBatch;
enum class WVObservationBatchKind : std::uint8_t;

struct WVObserverExecutionPlan {
  std::string fieldListAttribute;
  // Empty omits the legacy NetCDF name attribute for encodings that did not
  // historically persist one.
  std::string persistedName;
  std::vector<std::string> outputFields;
  // A record contributes a complete coefficient restart when this list
  // contains Ap, Am, and A0. Other restart representations can expose the
  // same generic family identities without a class-specific special case.
  std::vector<std::string> coefficientRestartFamilies;
};

// Construction-time binding surface for coarse integrated observer
// operations. Implementations invoke at most one callback while the runtime
// graph is compiled; numerical element loops remain statically specialized.
struct WVObserverIntegrationBinder {
  std::function<WVKernelStatus(const WVObserverRecord &)> advectedPositions;
  std::function<WVKernelStatus(const WVObserverRecord &)> advectedScalar;
};

// Stable source API v1 boundary for one exact MATLAB identity/version observing
// system. A factory returns an immutable shared implementation for one observer
// record; the resolved descriptor and its runtime owners retain that shared
// lifetime. Calls occur once per observer, RHS stage, or output event;
// implementations must not introduce virtual dispatch inside element loops.
class WVObservingSystem {
public:
  virtual ~WVObservingSystem() = default;
  virtual const std::string &typeIdentifier() const noexcept = 0;
  virtual std::uint32_t contractVersion() const noexcept = 0;

  virtual WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &blocks,
      std::map<std::string, std::size_t> &integratedBlockOwnerCounts) const = 0;
  virtual WVKernelStatus executionPlan(
      const WVObserverRecord &observer,
      WVObserverExecutionPlan &plan) const = 0;
  virtual WVKernelStatus bindIntegration(
      const WVObserverRecord &observer,
      const WVObserverIntegrationBinder &binder) const;
  virtual WVKernelStatus outputPlan(
      const WVObserverRecord &observer,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const;
  virtual WVKernelStatus prepareOccurrence(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOccurrencePreparationContext &context,
      WVObserverOccurrenceWorkspace &workspace) const;
  virtual WVKernelStatus observationBatch(
      const WVObserverRecord &observer, const WVObserverOutputPlan &plan,
      const WVObserverOutputEvaluationContext &context,
      WVObservationBatchKind kind, WVObservationBatch &batch) const;
  virtual std::size_t persistentBytes() const noexcept = 0;
};

} // namespace wavevortex::runtime
