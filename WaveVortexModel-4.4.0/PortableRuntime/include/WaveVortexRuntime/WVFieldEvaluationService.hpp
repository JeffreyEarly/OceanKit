#pragma once

#include "WaveVortexRuntime/WVDensityDiagnosticContract.hpp"
#include "WaveVortexRuntime/WVNoMotionProfileRecovery.hpp"

#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"
#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WaveVortexRuntime/WVObserverContracts.hpp"
#include "WaveVortexRuntime/generated/WVPortableVariableCatalog.hpp"
#include "WaveVortexRuntime/WVPortableVariablePlan.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

struct WVIntegrationState;
class WVIntegrationStateLayout;
class WVFieldEvaluationService;
class WVFieldEvaluationSession;
class WVConstantStratificationForcingEngine;
class WVBarotropicQGForcingEngine;
class WVStratifiedQGForcingEngine;
class WVHydrostaticForcingEngine;
class WVBoussinesqForcingEngine;


namespace detail {
class WVBarotropicQGFieldEvaluationAdapter;
class WVStratifiedFieldEvaluationAdapter;
class WVDiagnosticFieldPlan;
class WVForcingDiagnosticBinding;
class WVFieldEvaluationEventWorkspace;
class WVFieldEvaluationArena;
class WVFieldEvaluationEventScope;
class WVSampledMovingFieldPlan;
}

enum class WVFieldSamplingKind : std::uint8_t {
  fullGrid,
  fixedVerticalProfiles,
  positions
};

struct WVFieldSamplingRequest {
  WVFieldSamplingKind kind = WVFieldSamplingKind::fullGrid;
  // Fixed-profile indices use MATLAB's one-based grid-index convention.
  std::vector<std::size_t> xIndices;
  std::vector<std::size_t> yIndices;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
};

struct WVFieldRequest {
  std::string identifier;
  std::string fieldName;
  WVFieldSamplingRequest sampling;
};

struct WVFieldOutputSpecification {
  std::string identifier;
  std::string fieldName;
  WVFieldSamplingKind samplingKind = WVFieldSamplingKind::fullGrid;
  std::vector<std::size_t> dimensions;
  std::size_t elementCount = 0;
  bool isComplex = false;
};

struct WVFieldOutputView {
  double *data = nullptr;
  std::size_t elementCount = 0;
  WVComplex64 *complexData = nullptr;
};

// One field sampled from a caller-supplied moving-position array. Offsets and
// counts refer to the shared coordinate views passed to evaluateMoving().
struct WVMovingFieldRequest {
  std::string identifier;
  std::string fieldName;
  std::size_t positionOffset = 0;
  std::size_t positionCount = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
};

struct WVMovingPositionView {
  const double *x = nullptr;
  const double *y = nullptr;
  const double *z = nullptr;
  std::size_t positionCount = 0;
};

// Construction-time request for a field whose positions and output extents
// are supplied by each observation occurrence. The position-set slot is an
// already-resolved ordinal; event evaluation performs no name lookup.
struct WVEventFieldRequest {
  std::string identifier;
  std::string fieldName;
  std::size_t positionSetSlot = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
};

struct WVEventFieldOutputSpecification {
  std::string identifier;
  std::string fieldName;
  WVPortableVariable fieldIdentifier = WVPortableVariable::invalid;
  WVPortableNaturalRank naturalRank = WVPortableNaturalRank::volume;
  std::uint64_t dependencyMask = 0;
  std::size_t positionSetSlot = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
};

// Borrowed event-workspace coordinates. Extents use logical MATLAB order and
// must have a product equal to positionCount. An omitted extent view implies a
// one-dimensional {positionCount} result. Horizontal fields may omit z;
// volume fields require it when positionCount is nonzero.
struct WVEventPositionSetView {
  const double *x = nullptr;
  const double *y = nullptr;
  const double *z = nullptr;
  std::size_t positionCount = 0;
  const std::size_t *extents = nullptr;
  std::size_t extentCount = 0;
};

class WVEventFieldEvaluationPlan final {
public:
  const std::vector<WVEventFieldOutputSpecification> &outputs() const noexcept {
    return outputs_;
  }
  std::size_t outputCount() const noexcept { return outputs_.size(); }
  std::size_t positionSetCount() const noexcept { return positionSetCount_; }
  std::uint64_t requestedFieldMask() const noexcept {
    return requestedFieldMask_;
  }
  std::uint64_t dependencyMask() const noexcept { return dependencyMask_; }
  std::uint64_t fieldPlanFingerprint() const noexcept { return fingerprint_; }
  std::size_t persistentBytes() const noexcept;

private:
  struct ResolvedRequest {
    WVPortableVariable field = WVPortableVariable::invalid;
    WVPortableNaturalRank nativeRank = WVPortableNaturalRank::volume;
    std::uint64_t dependencyMask = 0;
    std::size_t positionSetSlot = 0;
    WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
    std::size_t outputIndex = 0;
  };

  WVTransformConstantStratificationConfiguration configuration_;
  std::vector<ResolvedRequest> requests_;
  std::vector<WVEventFieldOutputSpecification> outputs_;
  std::vector<std::uint8_t> requiresZByPositionSet_;
  std::size_t positionSetCount_ = 0;
  std::uint64_t requestedFieldMask_ = 0;
  std::uint64_t dependencyMask_ = 0;
  std::uint64_t fingerprint_ = 0;
  std::shared_ptr<const void> transformPlan_;
  std::size_t transformPlanBytes_ = 0;
  bool genericSampling_ = false;
  std::string configurationIdentifier_;
  WVDensityDiagnosticContract densityContract_;
  const WVFieldEvaluationService *owner_ = nullptr;
  std::shared_ptr<const std::uint8_t> planIdentity_;

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
  friend class detail::WVStratifiedFieldEvaluationAdapter;

};

class WVMovingFieldEvaluationPlan final {
public:
  const std::vector<WVFieldOutputSpecification> &outputs() const noexcept {
    return outputs_;
  }
  std::size_t outputCount() const noexcept { return outputs_.size(); }
  std::size_t positionCount() const noexcept { return positionCount_; }
  std::size_t persistentBytes() const noexcept;

private:
  struct ResolvedRequest {
    std::size_t primitiveChannel = 0;
    std::size_t positionOffset = 0;
    std::size_t positionCount = 0;
    WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
    std::size_t outputIndex = 0;
  };
  WVTransformConstantStratificationConfiguration configuration_;
  std::vector<ResolvedRequest> requests_;
  std::vector<WVFieldOutputSpecification> outputs_;
  std::size_t positionCount_ = 0;
  std::shared_ptr<const void> transformPlan_;
  std::size_t transformPlanBytes_ = 0;
  std::shared_ptr<const detail::WVSampledMovingFieldPlan> sampledPlan_;
  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
  friend class detail::WVStratifiedFieldEvaluationAdapter;

};

class WVFieldEvaluationPlan final {
public:
  const std::vector<WVFieldOutputSpecification> &outputs() const noexcept {
    return outputs_;
  }
  std::size_t outputCount() const noexcept { return outputs_.size(); }
  std::size_t persistentBytes() const noexcept;
  bool hasDensityDiagnostics() const noexcept;

private:
  using Field = WVPortableVariable;
  using NativeRank = WVPortableNaturalRank;

  struct PositionWeights {
    bool outsideInterpolationDomain = false;
    std::array<std::size_t, 2> xLinearIndices{};
    std::array<std::size_t, 2> yLinearIndices{};
    std::array<std::size_t, 2> zLinearIndices{};
    std::array<double, 2> xLinearWeights{};
    std::array<double, 2> yLinearWeights{};
    std::array<double, 2> zLinearWeights{};
    std::vector<double> xSplineWeights;
    std::vector<double> ySplineWeights;
    std::vector<double> zSplineWeights;
    std::size_t persistentBytes() const noexcept;
  };

  struct ResolvedRequest {
    Field field = Field::u;
    std::uint64_t dependencyMask = 0;
    NativeRank nativeRank = NativeRank::volume;
    WVFieldSamplingKind samplingKind = WVFieldSamplingKind::fullGrid;
    WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
    std::vector<std::size_t> profileXIndices;
    std::vector<std::size_t> profileYIndices;
    std::vector<PositionWeights> positionWeights;
    std::size_t outputIndex = 0;
    std::size_t persistentBytes() const noexcept;
  };

  friend class WVObserverOutputEvaluationService;
  std::shared_ptr<const detail::WVDiagnosticFieldPlan> diagnosticPlan_;
  WVTransformConstantStratificationConfiguration configuration_;
  std::vector<ResolvedRequest> requests_;
  std::vector<WVFieldOutputSpecification> outputs_;
  std::uint64_t requestedFieldMask_ = 0;
  std::uint64_t dependencyMask_ = 0;
  std::shared_ptr<const void> transformPlan_;
  std::size_t transformPlanBytes_ = 0;

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
  friend class detail::WVStratifiedFieldEvaluationAdapter;
  friend class detail::WVDiagnosticFieldPlan;
};

struct WVPreparedFieldOutputSpecification {
  std::size_t planOutputIndex = 0;
  std::size_t positionSetSlot = 0;
  std::vector<std::size_t> dimensions;
  std::size_t elementCount = 0;
};

struct WVPreparedFieldGeometryMetrics {
  std::size_t positionSetCount = 0;
  std::size_t positionCount = 0;
  std::size_t retainedBytes = 0;
  std::size_t liveBytes = 0;
};

// Event-scoped, retry-stable interpolation geometry. Coordinate storage is
// borrowed from the occurrence workspace; resolved interpolation weights and
// output extents are owned here until every route for that occurrence commits.
class WVPreparedFieldGeometry final {
public:
  const std::vector<WVPreparedFieldOutputSpecification> &
  outputs() const noexcept {
    return outputs_;
  }
  std::size_t outputCount() const noexcept { return outputs_.size(); }
  std::size_t positionSetCount() const noexcept { return positionSets_.size(); }
  std::size_t positionCount() const noexcept { return positionCount_; }
  WVEventPositionSetView positionSet(std::size_t slot) const noexcept;
  std::uint64_t fieldPlanFingerprint() const noexcept {
    return fieldPlanFingerprint_;
  }
  std::uint64_t geometryFingerprint() const noexcept {
    return geometryFingerprint_;
  }
  bool sameGeometry(const WVPreparedFieldGeometry &other) const noexcept;
  std::size_t retainedBytes() const noexcept;
  std::size_t liveBytes() const noexcept;
  std::size_t borrowedCoordinateBytes() const noexcept {
    return borrowedCoordinateBytes_;
  }
  WVPreparedFieldGeometryMetrics metrics() const noexcept;

private:
  struct PositionSet {
    const double *x = nullptr;
    const double *y = nullptr;
    const double *z = nullptr;
    std::size_t positionCount = 0;
    std::vector<std::size_t> extents;
  };

  WVFieldEvaluationPlan evaluationPlan_;
  std::vector<PositionSet> positionSets_;
  std::vector<WVPreparedFieldOutputSpecification> outputs_;
  std::size_t positionCount_ = 0;
  std::size_t borrowedCoordinateBytes_ = 0;
  std::uint64_t fieldPlanFingerprint_ = 0;
  std::shared_ptr<const std::uint8_t> planIdentity_;
  std::uint64_t geometryFingerprint_ = 0;
  std::shared_ptr<const void> transformGeometry_;
  std::size_t transformGeometryBytes_ = 0;

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
  friend class detail::WVStratifiedFieldEvaluationAdapter;

};

// One independently keyed occurrence in a coarse same-state evaluation.
// Plans and prepared geometries remain separately owned by their occurrence;
// the service unions their reconstruction dependencies without copying their
// coordinate weights.
struct WVEventFieldEvaluationBatchEntry {
  const WVEventFieldEvaluationPlan *plan = nullptr;
  const WVPreparedFieldGeometry *geometry = nullptr;
  WVFieldOutputView *outputs = nullptr;
  std::size_t outputCount = 0;
};

struct WVFieldEvaluationMetrics {
  std::size_t evaluationCount = 0;
  std::size_t coincidentBatchCount = 0;
  std::size_t lastPlanBytes = 0;
  std::size_t maximumPlanBytes = 0;
  std::size_t servicePersistentBytes = 0;
  std::size_t transformPersistentBytes = 0;
  std::size_t scratchCapacityBytes = 0;
  std::size_t scratchHighWaterBytes = 0;
  std::size_t movingInterpolationWorkspaceBytes = 0;
  std::size_t transformCount = 0;
  std::size_t fftExecutionCount = 0;
  std::size_t primitiveFieldEvaluationCount = 0;
  std::size_t primitiveFieldReuseCount = 0;
  std::size_t fullGridWriteCount = 0;
  std::size_t profileWriteCount = 0;
  std::size_t linearInterpolationCount = 0;
  std::size_t splineInterpolationCount = 0;
  std::size_t outputElementWriteCount = 0;
  std::size_t movingEvaluationCount = 0;
  std::size_t movingPositionCount = 0;
  std::size_t movingPrimitiveTransformCount = 0;
  std::size_t eventPlanCreationCount = 0;
  std::size_t eventPlanFieldResolutionCount = 0;
  std::size_t eventGeometryPreparationCount = 0;
  std::size_t eventEvaluationCount = 0;
  std::size_t eventBatchEvaluationCount = 0;
  std::size_t eventBatchOccurrenceCount = 0;
  std::size_t eventBatchOutputCount = 0;
  std::size_t eventFieldReuseCount = 0;
  std::size_t eventFieldWorkspaceLiveBytes = 0;
  std::size_t eventFieldWorkspaceHighWaterBytes = 0;
  std::size_t eventFieldArenaPlannedBytes = 0;
  std::size_t eventFieldArenaPeakBytes = 0;
  std::size_t eventBatchInvocationWorkspaceBytes = 0;
  std::size_t eventPositionSetCount = 0;
  std::size_t eventPositionCount = 0;
  std::size_t lastEventPlanBytes = 0;
  std::size_t maximumEventPlanBytes = 0;
  std::size_t lastPreparedGeometryRetainedBytes = 0;
  std::size_t maximumPreparedGeometryRetainedBytes = 0;
  std::size_t lastPreparedGeometryLiveBytes = 0;
  std::size_t maximumPreparedGeometryLiveBytes = 0;
  std::size_t diagnosticEvaluationCount = 0;
  // Unique underlying field requests, and repeated references served by them.
  std::size_t diagnosticPrimitiveOutputCount = 0;
  std::size_t diagnosticIntermediateReuseCount = 0;
  std::size_t diagnosticWorkspaceLiveBytes = 0;
  std::size_t diagnosticWorkspaceHighWaterBytes = 0;
  // Peak allocation beyond servicePersistentBytes. Prepared arena capacities
  // checked out during an event are excluded from this additive value.
  std::size_t additionalTransientHighWaterBytes = 0;
  std::size_t densityRecoveryCount = 0;
  std::size_t densityProfileConstructionCount = 0;
  std::size_t densityInversePassCount = 0;
  std::size_t densityAPEPassCount = 0;
  std::size_t densityAPVPassCount = 0;
  std::size_t densityAPVReuseCount = 0;
  std::size_t densityReuseCount = 0;
  std::size_t densityWorkspaceLiveBytes = 0;
  std::size_t densityWorkspaceHighWaterBytes = 0;
  std::size_t catalogBytes = portableVariableCatalogBytes();
  WVVariableEvaluationMetrics variableEvaluation;
};

// Explicit immutable-state lifetime for callers that issue several field
// queries for one output event. The session is tied to one service and one
// borrowed state; ending it invalidates all event-scoped result views.
class WVFieldEvaluationSession final {
public:
  WVFieldEvaluationSession();
  ~WVFieldEvaluationSession();
  WVFieldEvaluationSession(const WVFieldEvaluationSession &) = delete;
  WVFieldEvaluationSession &operator=(const WVFieldEvaluationSession &) = delete;
  WVFieldEvaluationSession(WVFieldEvaluationSession &&) noexcept;
  WVFieldEvaluationSession &operator=(WVFieldEvaluationSession &&) noexcept;
  bool active() const noexcept;
  const WVKernelStatus &status() const noexcept;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
  friend class WVFieldEvaluationService;
};

class WVFieldEvaluationService final {
public:
  static WVKernelStatus
  create(const WVTransformConstantStratificationConfiguration &configuration,
         std::unique_ptr<WVFFTEngine> engine,
         std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  createBorrowing(WVTransformConstantStratificationKernel &transform,
                  std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  create(const WVTransformBarotropicQGConfiguration &configuration,
         std::unique_ptr<WVFFTEngine> engine,
         std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  createBorrowing(WVTransformBarotropicQGKernel &transform,
                  std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  create(std::shared_ptr<const WVStratifiedModalSource> source,
         std::unique_ptr<WVFFTEngine> engine,
         std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  createBorrowing(WVTransformHydrostaticKernel &transform,std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus createBorrowing(WVTransformBoussinesqKernel &transform,std::unique_ptr<WVFieldEvaluationService> &service);
  static WVKernelStatus
  createBorrowing(WVTransformStratifiedQGKernel &transform,
                  std::unique_ptr<WVFieldEvaluationService> &service);

  // The borrowed engine and its resolved schedule must outlive the service.
  static WVKernelStatus createBorrowing(WVConstantStratificationForcingEngine&, std::unique_ptr<WVFieldEvaluationService>&);
  static WVKernelStatus createBorrowing(WVBarotropicQGForcingEngine&, std::unique_ptr<WVFieldEvaluationService>&);
  static WVKernelStatus createBorrowing(WVStratifiedQGForcingEngine&, std::unique_ptr<WVFieldEvaluationService>&);
  static WVKernelStatus createBorrowing(WVHydrostaticForcingEngine&, std::unique_ptr<WVFieldEvaluationService>&);
  static WVKernelStatus createBorrowing(WVBoussinesqForcingEngine&, std::unique_ptr<WVFieldEvaluationService>&);

  WVFieldEvaluationService(const WVFieldEvaluationService &) = delete;
  WVFieldEvaluationService &
  operator=(const WVFieldEvaluationService &) = delete;
  WVFieldEvaluationService(WVFieldEvaluationService &&) = delete;
  WVFieldEvaluationService &operator=(WVFieldEvaluationService &&) = delete;
  ~WVFieldEvaluationService();

  WVKernelStatus setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy policy);
  WVKernelStatus validateVariableEvaluationPolicyChange(
      WVVariableEvaluationPolicy policy) const noexcept;
  WVVariableEvaluationPolicy variableEvaluationPolicy() const noexcept {
    return variableEvaluationPolicy_;
  }
  bool evaluationSessionActive() const noexcept { return eventWorkspace_ != nullptr; }
  WVKernelStatus beginEvaluationSession(const WVIntegrationState &state,
                                        WVFieldEvaluationSession &session);
  // Prepare the event-owned physical, raw-tendency, and projected-coefficient
  // cache nodes used by the exact built-in nonlinear-advection fast path.
  // This is a setup operation and performs no scientific computation.
  WVKernelStatus prepareBuiltinNonlinearCoefficientEvaluation();
  // Evaluate the exact resolved built-in nonlinear spatial forcing in the
  // active immutable-state event, publishing its physical fields, raw
  // tendency, and projected coefficients through the ordinary event caches.
  WVKernelStatus evaluateBuiltinNonlinearCoefficients(
      const WVIntegrationState &state, WVFlux &flux);

  static std::vector<std::string> supportedFieldNames();
  const std::vector<WVPortableForcingVariableBinding>& forcingVariableBindings() const noexcept;
  std::string portableVariableConfiguration() const;
  WVKernelStatus createPlan(const std::vector<WVFieldRequest> &requests,
                            WVFieldEvaluationPlan &plan,
                            WVDensityDiagnosticContract densityContract = {}) const;
  // Extend the prepared workload of the active immutable-state event. This is
  // intentionally distinct from createPlan: it is valid only between producer
  // calls in an explicit reuse-policy evaluation session.
  WVKernelStatus createPlanForActiveEvaluation(
      const std::vector<WVFieldRequest> &requests,
      WVFieldEvaluationPlan &plan,
      WVDensityDiagnosticContract densityContract = {}) const;
  // A null selection evaluates every output. Otherwise one byte per output
  // selects its dependencies and writes; inactive output views are untouched.
  WVKernelStatus evaluate(const WVFieldEvaluationPlan &plan,
                          const WVState &state, WVFieldOutputView *outputs,
                          std::size_t outputCount,
                          const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus evaluate(const WVFieldEvaluationPlan &plan,
                          const WVIntegrationState &state,
                          WVFieldOutputView *outputs,
                          std::size_t outputCount,
                          const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus
  createMovingPlan(const std::vector<WVMovingFieldRequest> &requests,
                   WVMovingFieldEvaluationPlan &plan,
                   WVDensityDiagnosticContract densityContract = {}) const;
  // Selection also skips coordinate values belonging only to inactive requests.
  WVKernelStatus evaluateMoving(const WVMovingFieldEvaluationPlan &plan,
                                const WVState &state,
                                WVMovingPositionView positions,
                                WVFieldOutputView *outputs,
                                std::size_t outputCount, const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus evaluateMoving(const WVMovingFieldEvaluationPlan &plan,
                                const WVIntegrationState &state,
                                WVMovingPositionView positions,
                                WVFieldOutputView *outputs,
                                std::size_t outputCount, const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus evaluateMovingFromAdvectionFields(
      const WVMovingFieldEvaluationPlan &plan, const WVState &state,
      const WVRealFieldBundleConstView &advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount, const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus evaluateMovingFromAdvectionFields(
      const WVMovingFieldEvaluationPlan &plan,
      const WVIntegrationState &state,
      const WVRealFieldBundleConstView &advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount, const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus
  createEventPlan(const std::vector<WVEventFieldRequest> &requests,
                  WVEventFieldEvaluationPlan &plan,
                  WVDensityDiagnosticContract densityContract = {});
  WVKernelStatus
  prepareEventGeometry(const WVEventFieldEvaluationPlan &plan,
                       const WVEventPositionSetView *positionSets,
                       std::size_t positionSetCount,
                       WVPreparedFieldGeometry &geometry);
  WVKernelStatus evaluateEvent(const WVEventFieldEvaluationPlan &plan,
                               const WVPreparedFieldGeometry &geometry,
                               const WVState &state, WVFieldOutputView *outputs,
                               std::size_t outputCount);
  WVKernelStatus evaluateEvent(const WVEventFieldEvaluationPlan &plan,
                               const WVPreparedFieldGeometry &geometry,
                               const WVIntegrationState &state,
                               WVFieldOutputView *outputs,
                               std::size_t outputCount);
  WVKernelStatus
  evaluateEventBatch(const WVState &state,
                     const WVEventFieldEvaluationBatchEntry *entries,
                     std::size_t entryCount);
  WVKernelStatus
  evaluateEventBatch(const WVIntegrationState &state,
                     const WVEventFieldEvaluationBatchEntry *entries,
                     std::size_t entryCount);
  WVRealFieldBundleView advectionFieldStorage() noexcept;

  const WVTransformConstantStratificationConfiguration &
  configuration() const noexcept;
  const WVStratifiedModalGeometry* stratifiedGeometry() const noexcept;
  bool hasLegacyConfiguration() const noexcept { return transform_ != nullptr; }
  WVKernelStatus createStateLayout(
      const WVPortableObserverDescriptor &descriptor,
      WVIntegrationStateLayout &layout) const;
  bool isCompatibleWith(const WVIntegrationStateLayout &layout) const noexcept;
  bool isCompatibleWith(
      const WVFieldEvaluationService &other) const noexcept;
  const WVFieldEvaluationMetrics &metrics() const noexcept;
  WVVariableEvaluationMetrics activeVariableEvaluationMetrics() const noexcept;
  bool activeDensityRecoveryReport(
      WVNoMotionRecoveryReport& report) const noexcept;
  WVVariableProducerMetrics producerMetrics() const noexcept;
  std::size_t persistentBytes() const noexcept;

private:
  struct PlanInvocation {
    const WVFieldEvaluationPlan *plan = nullptr;
    WVFieldOutputView *outputs = nullptr;
    std::size_t outputCount = 0;
    const std::uint8_t *activeOutputs = nullptr;
  };
  friend class detail::WVDiagnosticFieldPlan;
  friend class detail::WVFieldEvaluationEventScope;
  detail::WVFieldEvaluationEventWorkspace* eventWorkspace_ = nullptr;
  WVVariableEvaluationPolicy variableEvaluationPolicy_ =
      WVVariableEvaluationPolicy::reuse;
  WVVariableEvaluationContext variableEvaluationContext_;
  WVKernelStatus variableEvaluationPreparationStatus_ = WVKernelStatus::ok();
  mutable std::unique_ptr<detail::WVFieldEvaluationArena> eventArena_;
  mutable bool eventArenaForcingPrepared_ = false;
  WVFieldEvaluationService();
  WVKernelStatus beginStateEvaluation(const WVIntegrationState&,const void* owner);
  WVKernelStatus addStateEvaluationView(const WVIntegrationState&,
      std::size_t componentIdentity);
  WVKernelStatus removeStateEvaluationView(const WVIntegrationState&,
      std::size_t componentIdentity);
  void endStateEvaluation() noexcept;
  WVKernelStatus prepareForcingEvaluationContext();
  WVKernelStatus prepareEventArena(std::size_t requestCount) const;
  WVKernelStatus prepareEventArena(const WVFieldEvaluationPlan&,
      std::uint32_t componentIdentity=0) const;
  WVKernelStatus prepareEventArena(const WVEventFieldEvaluationPlan&) const;
  WVKernelStatus prepareEventField(const WVVariableEvaluationKey&,
      std::size_t elements,bool complex) const;
  WVKernelStatus prepareDensityEventArena(std::size_t sampleCount,
      std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference,bool apvNeeded) const;
  WVKernelStatus prepareEventArenaForActiveEvaluation(
      const WVFieldEvaluationPlan&,std::uint32_t componentIdentity=0) const;
  WVKernelStatus prepareEventArenaImpl(const WVFieldEvaluationPlan&,
      std::uint32_t componentIdentity,bool activePreparation) const;
  WVKernelStatus prepareEventFieldForActiveEvaluation(
      const WVVariableEvaluationKey&,std::size_t elements,bool complex) const;
  WVKernelStatus prepareDensityEventArenaForActiveEvaluation(
      std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
      WVNoMotionReference reference,bool apvNeeded) const;
  WVKernelStatus createPlanImpl(const std::vector<WVFieldRequest>& requests,
      WVFieldEvaluationPlan& plan,WVDensityDiagnosticContract densityContract,
      bool prepareScientificDependencies) const;
  WVKernelStatus initializeScratch();
  WVKernelStatus evaluatePlanBatch(const PlanInvocation *invocations,
                                   std::size_t invocationCount,
                                   const WVState &state);
  WVKernelStatus evaluateMovingImpl(
      const WVMovingFieldEvaluationPlan &plan, const WVState &state,
      const WVRealFieldBundleConstView *advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount, const std::uint8_t *activeOutputs = nullptr);
  WVKernelStatus evaluateSampledMovingImpl(
      const WVMovingFieldEvaluationPlan &plan,
      const WVIntegrationState &state, WVMovingPositionView positions,
      WVFieldOutputView *outputs, std::size_t outputCount,
      const std::uint8_t *activeOutputs);
  WVKernelStatus samplePreparedField(const WVFieldEvaluationPlan &plan,
                                     const double *source,
                                     WVFieldOutputView output);
  WVKernelStatus evaluateSampledEventBatch(
      const WVIntegrationState &state,
      const WVEventFieldEvaluationBatchEntry *entries,
      std::size_t entryCount);
  WVFieldEvaluationMetrics &mutableMetrics() noexcept;
  class MovingWorkspace;
  class SampledMovingWorkspace;
  std::unique_ptr<WVTransformConstantStratificationKernel> ownedTransform_;
  WVTransformConstantStratificationKernel *transform_ = nullptr;
  std::unique_ptr<detail::WVBarotropicQGFieldEvaluationAdapter>
      barotropicQG_;
  std::unique_ptr<detail::WVStratifiedFieldEvaluationAdapter> stratified_;
  std::unique_ptr<detail::WVForcingDiagnosticBinding> forcing_;
  WVFieldEvaluationPlan builtinNonlinearFieldsPlan_;
  std::vector<double> builtinNonlinearPhysical_;
  std::vector<WVComplex64> builtinNonlinearCoefficients_;
  bool builtinNonlinearPrepared_ = false;
  std::unique_ptr<MovingWorkspace> movingWorkspace_;
  mutable std::unique_ptr<SampledMovingWorkspace> sampledMovingWorkspace_;
  std::vector<double> realScratch_;
  std::vector<WVComplex64> complexScratch_;
  std::vector<PlanInvocation> eventBatchInvocations_;
  mutable WVFieldEvaluationMetrics metrics_;
  mutable WVVariableProducerMetrics outputProducerMetrics_;
  bool executing_ = false;
  bool stateEvaluationActive_ = false;
};

} // namespace wavevortex::runtime
