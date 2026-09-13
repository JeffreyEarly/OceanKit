#pragma once

#include "WVFFTEngine.hpp"
#include "WVStratifiedModalSource.hpp"
#include "WVVariableExecutionOptions.hpp"

#include <array>
#include <memory>
#include <vector>

namespace wavevortex {
namespace kernel_detail { class WVPreparedModeExecutor; struct WVCompactConstantSchedule; }

// Schedule selection is an explicit construction policy, never
// mutable process state or part of the scientific checkpoint identity.
enum class WVConstantNonlinearFluxSchedule { frozenStreamed, compactCandidate };
struct WVConstantKernelExecutionOptions {
#if !defined(WV_KERNEL_COMPACT_CONSTANT_CANDIDATE) || WV_KERNEL_COMPACT_CONSTANT_CANDIDATE
    WVConstantNonlinearFluxSchedule schedule = WVConstantNonlinearFluxSchedule::compactCandidate;
#else
    WVConstantNonlinearFluxSchedule schedule = WVConstantNonlinearFluxSchedule::frozenStreamed;
#endif
#if defined(WV_KERNEL_COMPACT_HORIZONTAL_WORKERS)
    std::size_t horizontalOuterWorkers = WV_KERNEL_COMPACT_HORIZONTAL_WORKERS;
#else
    std::size_t horizontalOuterWorkers = 12;
#endif
#if defined(WV_KERNEL_COMPACT_POINTWISE_WORKERS)
    std::size_t pointwiseWorkers = WV_KERNEL_COMPACT_POINTWISE_WORKERS;
#else
    std::size_t pointwiseWorkers = 12;
#endif
};

enum class WVLaplacianDirection : std::uint8_t { horizontal, vertical };

enum class WVDynamicalField : std::uint8_t { u, v, w, eta };
enum class WVConstantFField : std::uint8_t { pi, psi, qgpv };

struct WVKernelMetrics {
    std::size_t descriptorBytes = 0;
    std::size_t planBytes = 0;
    std::size_t engineBytes = 0;
    std::size_t kernelManagementBytes = 0;
    std::size_t scratchCapacityBytes = 0;
    std::size_t scratchHighWaterBytes = 0;
    std::size_t halfSpectrumScratchCapacityBytes = 0;
    std::size_t realScratchCapacityBytes = 0;
    std::size_t planCount = 0;
    std::size_t executionCount = 0;
    std::size_t horizontalExecutionCount = 0;
    std::size_t verticalExecutionCount = 0;
    std::size_t nonlinearFluxCallCount = 0;
    std::size_t nonlinearFluxPhaseEvaluationCount = 0;
    std::size_t advectionVelocityReconstructionCount = 0;
    std::size_t scalarAdvectionCount = 0;
    std::size_t scalarAntialiasCount = 0;
    std::size_t bytesCopied = 0;
    std::size_t stateValidationCount = 0;
    std::size_t derivedValidationCount = 0;
    std::size_t phasePreparationCount = 0;
    std::array<std::size_t,4> tendencyReconstructionCount{};
    // Field indices align with WVHydrostaticField so producer metrics have one
    // identity across constant and variable stratification families.
    std::array<std::array<std::array<std::size_t,5>,4>,16> reconstructionCount{};
    double phaseSeconds = 0.0;
    double reconstructionSeconds = 0.0;
    double derivativeReconstructionSeconds = 0.0;
    double productSeconds = 0.0;
    double projectionSeconds = 0.0;
    double coefficientAssemblySeconds = 0.0;
    double derivativeCoefficientAssemblySeconds = 0.0;
    double coefficientProjectionSeconds = 0.0;
    double scalarForwardSeconds = 0.0;
    double scalarDerivativeAssemblySeconds = 0.0;
    double scalarVerticalDerivativeSeconds = 0.0;
    double scalarInverseSeconds = 0.0;
    double scalarProductSeconds = 0.0;
    double scalarAntialiasSeconds = 0.0;
};

class WVTransformConstantStratificationKernel {
public:
    static WVKernelStatus create(const WVTransformConstantStratificationConfiguration& configuration, std::unique_ptr<WVFFTEngine> engine, std::unique_ptr<WVTransformConstantStratificationKernel>& kernel, WVConstantKernelExecutionOptions options = {});

    ~WVTransformConstantStratificationKernel();
    WVTransformConstantStratificationKernel(const WVTransformConstantStratificationKernel&) = delete;
    WVTransformConstantStratificationKernel& operator=(const WVTransformConstantStratificationKernel&) = delete;
    WVTransformConstantStratificationKernel(WVTransformConstantStratificationKernel&&) = delete;
    WVTransformConstantStratificationKernel& operator=(WVTransformConstantStratificationKernel&&) = delete;

    const WVTransformConstantStratificationDescriptor& descriptor() const noexcept { return descriptor_; }
    const WVKernelMetrics& metrics() const noexcept;
    const std::string& engineIdentifier() const noexcept { return engineIdentifier_; }
    const std::string& engineLibraryIdentity() const noexcept { return engineLibraryIdentity_; }
    const char* nonlinearFluxScheduleIdentifier() const noexcept;
    const char* phaseImplementationIdentifier() const noexcept;
    const char* coefficientArithmeticModeIdentifier() const noexcept;
    const char* inverseNormalizationPlacementIdentifier() const noexcept;
    const char* optimizationImplementationIdentifier() const noexcept;
    std::size_t coefficientWorkerCount() const noexcept;
    std::size_t pointwiseWorkerCount() const noexcept;
    std::size_t retainedHorizontalWorkerCount() const noexcept;
    const char* retainedHorizontalScheduleIdentifier() const noexcept;
    std::size_t verticalExecutionRowCount() const noexcept;
    std::size_t phaseReservationBytes() const noexcept;
    void setStageInstrumentation(bool enabled) noexcept;
    std::size_t persistentBytes() const noexcept;
    std::size_t scratchBytes() const noexcept { return (halfSpectrumScratch_.size() + realScratch_.size()) * sizeof(double) + preparedPhase_.size()*sizeof(WVComplex64); }

    // The borrowed coefficient arrays remain immutable until the matching end.
    WVKernelStatus beginStateEvaluation(const WVState&);
    WVKernelStatus beginStateEvaluation(const WVState&, const void* evaluationOwner);
    WVKernelStatus addStateEvaluationView(const WVState&, const void* evaluationOwner,
        std::size_t componentIdentity = 0);
    WVKernelStatus removeStateEvaluationView(const WVState&, const void* evaluationOwner,
        std::size_t componentIdentity);
    WVKernelStatus endStateEvaluation();
    bool stateEvaluationActive() const noexcept { return stateEvaluationActive_; }
    WVKernelStatus validateStateEvaluation(const WVState&) const noexcept;
    // Borrowed until the next kernel operation; scoped calls reuse the active phase.
    WVKernelStatus preparedPhase(const WVState&, WVComplexConstView&);
    // Prepare the optional two-channel inverse during model setup. The
    // transform itself never allocates plans in an evaluation scope.
    WVKernelStatus prepareHorizontalVelocityTransform();
    // Prepare raw MATLAB vertical/Fourier/calculus plans explicitly. Kernel
    // construction and existing constant-model storage remain unchanged until
    // this setup operation is requested.
    WVKernelStatus prepareMatlabPrimitives();

    WVKernelStatus applyVertical(WVStratifiedModalOperator, WVComplexConstView,
        WVComplexView);
    WVKernelStatus applyVerticalColumn(WVStratifiedModalOperator,
        std::size_t retainedColumn, WVComplexConstView inputColumn,
        WVComplexView outputColumn);
    WVKernelStatus applyVerticalCalculus(WVRealConstView input,
        bool inputIsF, unsigned order, bool integral, WVRealView output);
    WVKernelStatus horizontalForward(WVRealVolumeConstView, WVComplexView);
    WVKernelStatus horizontalInverse(WVComplexConstView, WVRealVolumeView);
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView,
        bool xDerivative, WVRealVolumeView);
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView,
        bool xDerivative, unsigned order, WVRealVolumeView);

    WVKernelStatus transformUVEtaToWaveVortex(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients);
    WVKernelStatus transformUVWEtaToWaveVortex(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients);
    WVKernelStatus transformUVEtaToWaveVortex(WVRealVolumeConstView u,
        WVRealVolumeConstView v, WVRealVolumeConstView eta, double t, double t0,
        WVMutableCoefficients& coefficients);
    WVKernelStatus transformUVWEtaToWaveVortex(WVRealVolumeConstView u,
        WVRealVolumeConstView v, WVRealVolumeConstView w,
        WVRealVolumeConstView eta, double t, double t0,
        WVMutableCoefficients& coefficients);
    WVKernelStatus transformWaveVortexToUVWEta(const WVState& state, WVRealFieldBundleView& fields);
    // Reconstruct a validated derived coefficient tendency without admitting it
    // as a primary state view in an active immutable evaluation.
    WVKernelStatus transformCoefficientTendencyToUVWEta(
        const WVState& tendency, WVRealFieldBundleView& fields);
    WVKernelStatus transformWaveVortexToUVW(const WVState& state, WVRealFieldBundleView& fields);
    WVKernelStatus transformWaveVortexToUV(const WVState& state, WVRealFieldBundleView& fields);
    WVKernelStatus transformStateFieldDerivatives(const WVState& state, WVDynamicalField field, WVRealFieldBundleView& derivatives);
    WVKernelStatus transformToSpatialDomainWithFAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields);
    WVKernelStatus transformToSpatialDomainWithFAllDerivatives(WVConstantFField field,
        const WVComplexConstView& Apm, const WVComplexConstView& A0,
        WVRealFieldBundleView& fields, std::size_t componentIdentity = 0);
    WVKernelStatus transformToSpatialDomainWithGAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields);
    // Add physical velocity/displacement Laplacian forcing in coefficient space.
    // Reuses the descriptor's field/projection factors and caller-owned flux.
    WVKernelStatus addLaplacianDamping(const WVState& state, double nu, double kappa,
        WVLaplacianDirection direction, WVFlux& flux);
    WVKernelStatus nonlinearFlux(const WVState& state, WVFlux& flux,
        WVStateDerivativeAccess* derivativeAccess = nullptr);
    WVKernelStatus nonlinearFluxWithAdvectionFields(const WVState& state, WVFlux& flux,
        WVRealFieldBundleView& advectionFields, WVStateDerivativeAccess* derivativeAccess = nullptr);
    // Optional observation output receives raw spatial tendencies before modal
    // projection, in [u,v,eta] or [u,v,w,eta] order. Storage is caller owned.
    // With projectFlux=false, spatialTendency is required and flux is untouched.
    WVKernelStatus nonlinearFluxUsingAdvectionFields(const WVState& state, WVFlux& flux, const WVRealFieldBundleConstView& advectionFields,
        WVRealFieldBundleView* spatialTendency = nullptr, bool projectFlux = true,
        WVStateDerivativeAccess* derivativeAccess = nullptr);
    // Call at setup when scalar advection is configured, before repeated RHS calls.
    WVKernelStatus prepareScalarAdvection();
    // MATLAB diffX/diffY/diffZG of an arbitrary full-grid G scalar. Horizontal
    // derivatives retain all non-Nyquist modes; vertical calculus retains Nj.
    // Reuses existing plans/scratch and publishes [dx,dy,dz] only on success.
    WVKernelStatus transformGGridScalarDerivatives(const WVRealVolumeConstView& scalar, WVRealFieldBundleView& derivatives);
    WVKernelStatus advectFGridScalar(const WVRealVolumeConstView& scalar, const WVRealFieldBundleConstView& advectionFields, bool shouldAntialias, WVRealVolumeView& rightHandSide);

private:
    WVTransformConstantStratificationKernel();
    WVKernelStatus validateStateContents(const WVState&,
        bool allowPreparedExecutor = false);
    WVKernelStatus validateState(const WVState&,
        bool allowPreparedExecutor = false);
    WVKernelStatus validateStateForCall(const WVState&);
    WVKernelStatus validateStateAndFluxForCall(const WVState&,const WVFlux&);
    WVKernelStatus validateMutableOutputOutsidePreparedState(const WVMutableCoefficients&) const;
    WVKernelStatus validateMutableOutputOutsidePreparedState(const void*,
        std::size_t) const;
    bool matchesStateEvaluation(const WVState&) const noexcept;
    std::size_t stateEvaluationComponent(const WVState&) const noexcept;
    WVKernelStatus prepareStatePhase(const WVState&);
    WVComplexConstView phaseForPreparedState() const noexcept;
    WVKernelStatus preparePlans();
    WVKernelStatus transformUVEtaToWaveVortexImpl(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues = {});
    WVKernelStatus transformUVWEtaToWaveVortexImpl(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues = {});
    WVKernelStatus transformWaveVortexToUVWEtaImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients = nullptr, WVComplexConstView phaseValues = {}, WVComplex64* generatedPhase = nullptr, bool countPrimaryReconstruction = true);
    WVKernelStatus transformWaveVortexToUVWImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients, WVComplexConstView phaseValues = {}, WVComplex64* generatedPhase = nullptr);
    WVKernelStatus transformWaveVortexToUVImpl(const WVState& state, WVRealFieldBundleView& fields, WVComplexConstView phaseValues = {}, WVComplex64* generatedPhase = nullptr);
    WVKernelStatus transformToSpatialDomainWithDerivativesImpl(const WVCoefficients& evolvedCoefficients, std::size_t target, WVRealFieldBundleView& derivatives);
    WVKernelStatus transformToSpatialDomainWithDerivativesFromStateImpl(const WVState& state, WVComplexConstView phaseValues, std::size_t target, WVRealFieldBundleView& derivatives, WVComplex64* generatedPhase = nullptr, const std::array<bool,3>* derivativeMask = nullptr);
    WVKernelStatus projectSingleFluxTargetImpl(const WVRealFieldBundleConstView& field, std::size_t target, WVComplexConstView phaseValues, WVFlux& flux);
    WVKernelStatus nonlinearFluxImpl(const WVState& state, WVFlux& flux, WVRealFieldBundleView* advectionFields, bool advectionFieldsPrepared = false,
        WVRealFieldBundleView* spatialTendency = nullptr, bool projectFlux = true,
        WVComplexConstView phaseValues = {}, WVStateDerivativeAccess* derivativeAccess = nullptr);
    WVKernelStatus ensureScalarInversePlan();
    WVKernelStatus ensureScalarPlans();
    WVKernelStatus prepareCompactScalarDerivatives(const WVRealVolumeConstView& scalar, bool sine);
    WVFFTPlan* scalarPlan(std::size_t index) const noexcept;
    const WVHalfSpectrumMappings& executionMapping() const noexcept;
    std::size_t executionRowCount() const noexcept;
    WVKernelStatus antialiasScalarInPlace(WVRealVolumeView& scalar);
    WVTransformConstantStratificationDescriptor descriptor_;
    std::shared_ptr<WVFFTEngine> engine_;
    std::string engineIdentifier_;
    std::string engineLibraryIdentity_;
    std::vector<std::unique_ptr<WVFFTPlan>> plans_;
    std::unique_ptr<WVFFTPlan> scalarInversePlan_;
    // Retained-row DCT/DST and one-column DCT/DST plans.
    std::array<std::unique_ptr<WVFFTPlan>,4> matlabPlans_;
    std::unique_ptr<kernel_detail::WVCompactConstantSchedule> compact_;
    std::vector<std::uint8_t> scalarAntialiasRows_;
    std::vector<double> halfSpectrumScratch_;
    std::vector<double> realScratch_;
    std::vector<WVComplex64> preparedPhase_;
    WVState preparedState_{};
    std::array<WVState,5> preparedStateViews_{};
    std::array<std::size_t,5> preparedStateComponents_{};
    std::size_t preparedStateViewCount_ = 0;
    const void* preparedStateOwner_ = nullptr;
    bool stateEvaluationActive_ = false;
    bool preparedPhaseReady_ = false;
    bool matlabPrimitivesPrepared_ = false;
    mutable WVKernelMetrics metrics_;
    std::unique_ptr<kernel_detail::WVPreparedModeExecutor> coefficientExecutor_;
    bool executing_ = false;
    bool stageInstrumentationEnabled_ = false;
};

} // namespace wavevortex
