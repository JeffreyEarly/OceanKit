#pragma once

#include "WaveVortexRuntime/WVIntegrationContracts.hpp"
#include "WaveVortexRuntime/WVForcing.hpp"
#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"

#include <memory>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVExtensionCatalog;

class WVConstantStratificationRightHandSideContext final {
public:
    bool hasAdvectionFields() const noexcept { return advectionFields_.data != nullptr; }
    WVRealFieldBundleConstView advectionFields() const noexcept { return advectionFields_; }
    std::uint64_t generation() const noexcept { return generation_; }

private:
    const void* owner_ = nullptr;
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
    WVShape2D stateShape() const noexcept { return kernel_->descriptor().spectralShape(); }
    WVKernelStatus createErrorPolicy(double absoluteToleranceScale, std::unique_ptr<WVIntegrationErrorPolicy>& policy) const;

    const WVTransformConstantStratificationKernel& kernel() const noexcept { return *kernel_; }
    WVTransformConstantStratificationKernel& kernel() noexcept { return *kernel_; }
    const WVForcingEngineMetrics& metrics() const noexcept { return metrics_; }
    const std::string& scheduleIdentifier() const noexcept { return scheduleIdentifier_; }
    std::size_t persistentBytes() const noexcept;

private:
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
    void clearEvaluationWorkspace() noexcept;

    std::unique_ptr<WVTransformConstantStratificationKernel> kernel_;
    std::shared_ptr<const WVExtensionCatalog> catalog_;
    std::vector<std::unique_ptr<WVForcing>> forcing_;
    std::vector<double> physicalFields_;
    std::vector<double> forcingFields_;
    std::vector<WVComplex64> temporaryFlux_;
    WVForcingEngineMetrics metrics_;
    std::string scheduleIdentifier_;
    bool physicalFieldsValid_ = false;
    bool executing_ = false;
    std::uint64_t evaluationGeneration_ = 0;
    friend class WVForcingExecutionContext;
};

} // namespace wavevortex::runtime
