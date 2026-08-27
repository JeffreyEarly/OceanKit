#pragma once

#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"
#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WaveVortexRuntime/WVObserverContracts.hpp"
#include "WaveVortexRuntime/generated/WVPortableVariableCatalog.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

struct WVIntegrationState;
class WVIntegrationStateLayout;

namespace detail {
class WVBarotropicQGFieldEvaluationAdapter;
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
};

struct WVFieldOutputView {
  double *data = nullptr;
  std::size_t elementCount = 0;
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

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
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
  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
};

class WVFieldEvaluationPlan final {
public:
  const std::vector<WVFieldOutputSpecification> &outputs() const noexcept {
    return outputs_;
  }
  std::size_t outputCount() const noexcept { return outputs_.size(); }
  std::size_t persistentBytes() const noexcept;

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
    NativeRank nativeRank = NativeRank::volume;
    WVFieldSamplingKind samplingKind = WVFieldSamplingKind::fullGrid;
    WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
    std::vector<std::size_t> profileXIndices;
    std::vector<std::size_t> profileYIndices;
    std::vector<PositionWeights> positionWeights;
    std::size_t outputIndex = 0;
    std::size_t persistentBytes() const noexcept;
  };

  WVTransformConstantStratificationConfiguration configuration_;
  std::vector<ResolvedRequest> requests_;
  std::vector<WVFieldOutputSpecification> outputs_;
  std::uint64_t requestedFieldMask_ = 0;
  std::uint64_t dependencyMask_ = 0;
  std::shared_ptr<const void> transformPlan_;
  std::size_t transformPlanBytes_ = 0;

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
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
  std::uint64_t geometryFingerprint_ = 0;
  std::shared_ptr<const void> transformGeometry_;
  std::size_t transformGeometryBytes_ = 0;

  friend class WVFieldEvaluationService;
  friend class detail::WVBarotropicQGFieldEvaluationAdapter;
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
  std::size_t eventBatchInvocationWorkspaceBytes = 0;
  std::size_t eventPositionSetCount = 0;
  std::size_t eventPositionCount = 0;
  std::size_t lastEventPlanBytes = 0;
  std::size_t maximumEventPlanBytes = 0;
  std::size_t lastPreparedGeometryRetainedBytes = 0;
  std::size_t maximumPreparedGeometryRetainedBytes = 0;
  std::size_t lastPreparedGeometryLiveBytes = 0;
  std::size_t maximumPreparedGeometryLiveBytes = 0;
  std::size_t catalogBytes = portableVariableCatalogBytes();
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

  WVFieldEvaluationService(const WVFieldEvaluationService &) = delete;
  WVFieldEvaluationService &
  operator=(const WVFieldEvaluationService &) = delete;
  WVFieldEvaluationService(WVFieldEvaluationService &&) = delete;
  WVFieldEvaluationService &operator=(WVFieldEvaluationService &&) = delete;
  ~WVFieldEvaluationService();

  static std::vector<std::string> supportedFieldNames();
  WVKernelStatus createPlan(const std::vector<WVFieldRequest> &requests,
                            WVFieldEvaluationPlan &plan) const;
  WVKernelStatus evaluate(const WVFieldEvaluationPlan &plan,
                          const WVState &state, WVFieldOutputView *outputs,
                          std::size_t outputCount);
  WVKernelStatus evaluate(const WVFieldEvaluationPlan &plan,
                          const WVIntegrationState &state,
                          WVFieldOutputView *outputs,
                          std::size_t outputCount);
  WVKernelStatus
  createMovingPlan(const std::vector<WVMovingFieldRequest> &requests,
                   WVMovingFieldEvaluationPlan &plan) const;
  WVKernelStatus evaluateMoving(const WVMovingFieldEvaluationPlan &plan,
                                const WVState &state,
                                WVMovingPositionView positions,
                                WVFieldOutputView *outputs,
                                std::size_t outputCount);
  WVKernelStatus evaluateMoving(const WVMovingFieldEvaluationPlan &plan,
                                const WVIntegrationState &state,
                                WVMovingPositionView positions,
                                WVFieldOutputView *outputs,
                                std::size_t outputCount);
  WVKernelStatus evaluateMovingFromAdvectionFields(
      const WVMovingFieldEvaluationPlan &plan, const WVState &state,
      const WVRealFieldBundleConstView &advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount);
  WVKernelStatus evaluateMovingFromAdvectionFields(
      const WVMovingFieldEvaluationPlan &plan,
      const WVIntegrationState &state,
      const WVRealFieldBundleConstView &advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount);
  WVKernelStatus
  createEventPlan(const std::vector<WVEventFieldRequest> &requests,
                  WVEventFieldEvaluationPlan &plan);
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
  bool hasLegacyConfiguration() const noexcept { return transform_ != nullptr; }
  WVKernelStatus createStateLayout(
      const WVPortableObserverDescriptor &descriptor,
      WVIntegrationStateLayout &layout) const;
  bool isCompatibleWith(const WVIntegrationStateLayout &layout) const noexcept;
  bool isCompatibleWith(
      const WVFieldEvaluationService &other) const noexcept;
  const WVFieldEvaluationMetrics &metrics() const noexcept;
  std::size_t persistentBytes() const noexcept;

private:
  struct PlanInvocation {
    const WVFieldEvaluationPlan *plan = nullptr;
    WVFieldOutputView *outputs = nullptr;
    std::size_t outputCount = 0;
  };
  WVFieldEvaluationService() = default;
  WVKernelStatus initializeScratch();
  WVKernelStatus evaluatePlanBatch(const PlanInvocation *invocations,
                                   std::size_t invocationCount,
                                   const WVState &state);
  WVKernelStatus evaluateMovingImpl(
      const WVMovingFieldEvaluationPlan &plan, const WVState &state,
      const WVRealFieldBundleConstView *advectionFields,
      WVMovingPositionView positions, WVFieldOutputView *outputs,
      std::size_t outputCount);
  class MovingWorkspace;
  std::unique_ptr<WVTransformConstantStratificationKernel> ownedTransform_;
  WVTransformConstantStratificationKernel *transform_ = nullptr;
  std::unique_ptr<detail::WVBarotropicQGFieldEvaluationAdapter>
      barotropicQG_;
  std::unique_ptr<MovingWorkspace> movingWorkspace_;
  std::vector<double> realScratch_;
  std::vector<WVComplex64> complexScratch_;
  std::vector<PlanInvocation> eventBatchInvocations_;
  WVFieldEvaluationMetrics metrics_;
  bool executing_ = false;
};

} // namespace wavevortex::runtime
