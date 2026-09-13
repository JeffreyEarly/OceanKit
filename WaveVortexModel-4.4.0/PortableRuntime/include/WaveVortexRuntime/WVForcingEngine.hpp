#pragma once
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"
#include "WVForcingTendency.hpp"

#include "WaveVortexRuntime/WVIntegrationContracts.hpp"
#include "WaveVortexRuntime/WVForcing.hpp"
#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"

#include <memory>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {
namespace detail { class WVForcingDiagnosticBinding; }

class WVExtensionCatalog;
struct WVForcingEvaluationDependencies;

class WVConstantStratificationRightHandSideContext final {
public:
    bool hasAdvectionFields() const noexcept { return evaluation_ && evaluation_->active() && evaluation_->generation()==generation_ && advectionFields_.data; }
    WVRealFieldBundleConstView advectionFields() const noexcept { return hasAdvectionFields() ? advectionFields_ : WVRealFieldBundleConstView{}; }
    std::uint64_t generation() const noexcept { return generation_; }

private:
    const void* owner_ = nullptr;
    const WVVariableEvaluationContext* evaluation_ = nullptr;
    WVRealFieldBundleConstView advectionFields_;
    std::uint64_t generation_ = 0;
    friend class WVConstantStratificationForcingEngine;
};

struct WVForcingEngineMetrics {
    std::size_t scheduleBytes = 0;
    std::size_t derivedOperatorBytes = 0;
    std::size_t workspaceCapacityBytes = 0;
    std::size_t workspaceHighWaterBytes = 0;
    std::size_t evaluationCount = 0;
    std::size_t restoredCoefficientCount = 0;
    std::size_t resolvedSpatialCount = 0;
    std::size_t resolvedSpectralCount = 0;
    std::size_t resolvedAmplitudeCount = 0;
    std::size_t physicalFieldReconstructionCount = 0;
    std::size_t physicalFieldReuseCount = 0;
    std::size_t nonlinearProducerCount = 0;
    std::size_t horizontalSpeedReductionCount = 0;
    std::size_t verticalSpeedReductionCount = 0;
    std::array<std::array<std::size_t,4>,4> gridCalculusProducerCount{};
    std::array<std::size_t,2> constantLaplacianProducerCount{};
    std::size_t spatialTendencyProjectionCount = 0;
    std::size_t accumulatorClearElementWrites = 0;
    std::size_t spatialTendencyClearElementWrites = 0;
    std::size_t temporaryFluxClearElementWrites = 0;
    std::size_t kernelOutputInitializationElementWrites = 0;
    std::size_t temporaryAccumulationElementReads = 0;
    std::size_t temporaryAccumulationElementWrites = 0;
    std::size_t outputCopyElementReads = 0;
    std::size_t outputCopyElementWrites = 0;
    std::size_t stateConstraintElementWrites = 0;
    std::size_t workspaceLiveBytes = 0;
    std::size_t workspaceMaximumLiveBytes = 0;
};

// Evaluate a validated, immutable portable forcing schedule around the shared
// constant-stratification numerical kernel. Schedule construction resolves
// stage and priority order; evaluate() performs no forcing dispatch discovery.
class WVConstantStratificationForcingEngine final {
public:
    // Indices refer to the already resolved stage/priority order.
    std::size_t forcingCount() const noexcept { return forcing_.size(); }
    const WVForcing* forcingInstance(std::size_t index) const noexcept {
        return index<forcing_.size() ? forcing_[index].get() : nullptr;
    }
    const WVForcingEvaluationDependencies*
    forcingEvaluationDependencies(std::size_t index) const noexcept;
    // Optional u/v/w fields must describe this exact state and time.
    // They are borrowed for this invocation and must not alias state or outputs.
    WVKernelStatus evaluateForcingTendencies(const WVState&,
        const WVForcingTendencyOutput*,std::size_t,
        const WVRealFieldBundleConstView* preparedPhysical = nullptr,
        detail::WVForcingDiagnosticWorkspace* session = nullptr);
    const WVForcingTendencyMetrics& tendencyMetrics() const noexcept { return tendencyMetrics_; }

    // Validate the frozen schedule without constructing transforms, plans, or
    // array-sized derived operators and workspaces.
    static WVKernelStatus validateSchedule(
        const WVTransformConstantStratificationConfiguration& configuration,
        const WVFrozenForcingSchedule& schedule,
        WVShape2D coefficientShape,
        const WVExtensionCatalog& catalog);

    static WVKernelStatus create(
        const WVTransformConstantStratificationConfiguration& configuration,
        const WVFrozenForcingSchedule& schedule,
        std::shared_ptr<const WVExtensionCatalog> catalog,
        std::unique_ptr<WVFFTEngine> fftEngine,
        std::unique_ptr<WVConstantStratificationForcingEngine>& forcingEngine);

    ~WVConstantStratificationForcingEngine();
    WVConstantStratificationForcingEngine(const WVConstantStratificationForcingEngine&) = delete;
    WVConstantStratificationForcingEngine& operator=(const WVConstantStratificationForcingEngine&) = delete;

    WVKernelStatus beginStateEvaluation(const WVState&);
    void endStateEvaluation() noexcept;
    bool stateEvaluationActive() const noexcept { return evaluation_.active(); }
    WVKernelStatus validateStateEvaluation(const WVState&) const;
    WVKernelStatus nonlinearFlux(const WVState& state, WVFlux& flux);
    WVKernelStatus evaluateRightHandSideWithContext(
        const WVState& state, WVFlux& rightHandSide,
        WVRealFieldBundleView& advectionFieldStorage,
        WVConstantStratificationRightHandSideContext& context);
    WVKernelStatus advectFGridScalar(
        const WVConstantStratificationRightHandSideContext& context,
        const WVRealVolumeConstView& scalar, bool shouldAntialias,
        WVRealVolumeView& rightHandSide);
    WVStateConstraintResult restoreForcingAmplitudes(WVMutableCoefficients& coefficients);
    const WVForcingPreparation& preparation() const noexcept { return preparation_; }
    WVShape2D stateShape() const noexcept { return kernel_->descriptor().spectralShape(); }
    WVKernelStatus createErrorPolicy(double absoluteToleranceScale, std::unique_ptr<WVIntegrationErrorPolicy>& policy) const;

    const WVTransformConstantStratificationKernel& kernel() const noexcept { return *kernel_; }
    WVTransformConstantStratificationKernel& kernel() noexcept { return *kernel_; }
    const WVForcingEngineMetrics& metrics() const noexcept { return metrics_; }
    const std::string& scheduleIdentifier() const noexcept { return scheduleIdentifier_; }
    std::size_t persistentBytes() const noexcept;

    // Linear evolution retains instances for diagnostics and amplitude constraints.
    // Only their ordinary coefficient RHS contributions are disabled.
  WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy);
  WVKernelStatus validateVariableEvaluationPolicyChange(
      WVVariableEvaluationPolicy policy) const noexcept;
  const WVVariableEvaluationMetrics& variableEvaluationMetrics() const noexcept { return evaluation_.metrics(); }
    void setLinearDynamics(bool linear) noexcept { linearDynamics_ = linear; }

private:
    friend class detail::WVForcingDiagnosticBinding;
    WVKernelStatus evaluateForcingTendenciesImpl(const WVState&,
        const WVForcingTendencyOutput*,std::size_t,
        const WVRealFieldBundleConstView*,detail::WVForcingDiagnosticWorkspace*,
        WVFlux*);
    bool linearDynamics_ = false;
  WVVariableEvaluationContext evaluation_;
  WVVariableEvaluationPolicy evaluationPolicy_=WVVariableEvaluationPolicy::reuse;
    WVConstantStratificationForcingEngine() = default;

    WVKernelStatus initialize(const WVFrozenForcingSchedule& schedule);
    WVKernelStatus nonlinearFluxImpl(const WVState& state, WVFlux& flux, WVRealFieldBundleView* advectionFields, WVConstantStratificationRightHandSideContext* context);
    WVKernelStatus ensurePhysicalFields(const WVState& state, WVRealFieldBundleConstView& fields, WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared);
    WVRealFieldBundleView clearedSpatialTendency();
    WVKernelStatus projectSpatialTendency(const WVState& state, const WVRealFieldBundleConstView& tendency, WVFlux& flux);
    WVKernelStatus addProjectedSpatialTendency(const WVState& state, const WVRealFieldBundleConstView& tendency, WVFlux& flux, bool& outputInitialized);
    WVKernelStatus addAdaptiveDamping(const WVState& state, const std::vector<double>& damping, WVFlux& flux, WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared);
    WVKernelStatus addPseudoTopographicGeneration(const WVState& state, const WVPseudoTopographicOperators& operators, WVFlux& flux);
    void addBetaPlaneAdvection(const WVState& state, const std::vector<WVComplex64>& betaA0, WVFlux& flux) const;
    WVKernelStatus addLinearCoefficientTendency(const WVState& state, double rate, WVFlux& flux) const;
    void initializeOutputWithZeros(WVFlux& flux, bool& outputInitialized);
    WVKernelStatus addNonlinearFlux(const WVState& state, WVFlux& flux, bool& outputInitialized, WVRealFieldBundleView* externalFields, bool& externalFieldsPrepared);
    WVKernelStatus diagnosticLaplacian(const WVState&,double,double,WVLaplacianDirection,WVFlux&,bool&);

    std::unique_ptr<WVTransformConstantStratificationKernel> kernel_;
    std::shared_ptr<const WVExtensionCatalog> catalog_;
    std::vector<std::unique_ptr<WVForcing>> forcing_;
    std::vector<double> physicalFields_;
    std::vector<double> forcingFields_;
    std::vector<WVComplex64> temporaryFlux_, nonlinearCache_;
    WVForcingEngineMetrics metrics_;
    WVForcingPreparation preparation_;
    std::string scheduleIdentifier_;
    WVState evaluationState_{};
    WVRealFieldBundleView evaluationFields_{};
    std::size_t physicalFieldCount_=3;
    std::array<std::size_t,2> constantLaplacianUseCount_{};
    double horizontalMaximum_=0;
    bool executing_ = false;
    detail::WVForcingDiagnosticWorkspace* diagnosticWorkspace_=nullptr;
    WVForcingTendencyMetrics tendencyMetrics_;
    friend class WVForcingExecutionContext;
};

} // namespace wavevortex::runtime
