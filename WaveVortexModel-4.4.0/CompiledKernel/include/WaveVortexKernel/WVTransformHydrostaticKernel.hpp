#pragma once
#include "WVStratifiedModalSource.hpp"
#include "WVVariableExecutionOptions.hpp"
#include <array>
#include <atomic>
#include <functional>

namespace wavevortex {
namespace spectral_detail { class WVVariableComplexBuffer; }
namespace kernel_detail { class WVPreparedModeExecutor; class WVPreparedFieldCache; }
inline constexpr const char* WVHydrostaticKernelContract = "wave-vortex-hydrostatic-kernel-v1";
enum class WVHydrostaticField { u, v, w, eta, pi, p, psi, qgpv, rhoE, rhoTotal, zetaX, zetaY, zetaZ, ssh, ssu, ssv };
enum class WVHydrostaticDerivative { value, x, y, z };
enum class WVHydrostaticComponent { all, wave, inertial, geostrophic, meanDensityAnomaly };
enum class WVHydrostaticFamily { F, G };

struct WVHydrostaticModeFactors {
    double omega = 0;
    WVComplex64 UAp{}, VAp{}, WAp{}, UA0{}, VA0{}, ApmD{};
    double NAp = 0, NA0 = 0, PA0 = 0, ApmN = 0, A0Z = 0, A0N = 0;
    double waveEnergy = 0, balancedEnergy = 0, psi = 0, qgpv = 0, enstrophy = 0;
    bool wave = false, inertial = false, geostrophic = false, meanDensityAnomaly = false;
};
struct WVHydrostaticStorage {
    std::size_t sharedScientificBytes = 0, preparedBytes = 0, workspaceBytes = 0;
    std::size_t spectralScratchBytes = 0, realScratchBytes = 0, factorBytes = 0;
    std::size_t providerBytesLowerBound = 0, planBytesLowerBound = 0;
};
struct WVHydrostaticKernelMetrics {
    std::size_t stateValidationCount = 0;
    std::size_t derivedValidationCount = 0;
    std::size_t phasePreparationCount = 0;
    std::size_t coefficientAssemblyCount = 0, verticalPreparationCount = 0;
    std::size_t verticalOperatorExecutionCount = 0;
    std::size_t derivativeAdvectionConsumerCount = 0;
    std::size_t tiledNonlinearCount = 0, tiledColumnInverseCount = 0;
    std::size_t tiledRowInverseCount = 0, tiledReusedColumnCount = 0;
    std::size_t horizontalSpectrumReuseCount = 0, preparedVerticalDerivativeCount = 0;
    std::array<std::size_t,4> tendencyReconstructionCount{};
    std::array<std::size_t,16> fieldReconstructionCount{};
    std::array<std::array<std::array<std::size_t,5>,4>,16> reconstructionCount{};
};

// Immutable scientific source plus one owned mutable workspace. Coefficients
// are interleaved [Nj,Nkl]; fields are column-major [Nx,Ny,Nz]. Separate kernel
// instances are required for concurrent calls. No eigensolver or persistence I/O.
class WVTransformHydrostaticKernel final {
public:
    ~WVTransformHydrostaticKernel();
    using MatrixBackendFactory = std::function<WVKernelStatus(std::unique_ptr<WVVerticalMatrixBackend>&)>;
    static WVKernelStatus create(std::shared_ptr<const WVStratifiedModalSource>,
        std::unique_ptr<WVFFTEngine>, std::unique_ptr<WVTransformHydrostaticKernel>&,
        MatrixBackendFactory = WVCreateScalarMatrixBackend, WVVariableExecutionOptions = {});
    const WVVariableExecutionOptions& executionOptions() const noexcept { return executionOptions_; }
    const char* horizontalScheduleIdentifier() const noexcept { return horizontalWorkspace_->scheduleIdentifier(); }
    const WVStratifiedModalGeometry& geometry() const noexcept { return source_->geometry(); }
    const std::vector<WVHydrostaticModeFactors>& factors() const noexcept { return factors_; }
    const WVHydrostaticStorage& storage() const noexcept;
    const WVHydrostaticKernelMetrics& metrics() const noexcept { return metrics_; }
    void resetMetrics() noexcept { metrics_ = {}; }
    const std::string& engineIdentifier() const noexcept { return engineIdentifier_; }
    const std::string& engineLibraryIdentity() const noexcept { return engineLibraryIdentity_; }
    const char* matrixBackendIdentifier() const noexcept { return vertical_[0]->backendIdentifier(); }
    WVShape2D spectralShape() const noexcept { return {geometry().Nj,geometry().Nkl}; }
    WVShape3D spatialShape() const noexcept { return {geometry().Nx,geometry().Ny,geometry().Nz}; }
    std::size_t persistentBytes() const noexcept;

    // Begin an explicit evaluation of one immutable borrowed state. State
    // validation and phase preparation are reused until endStateEvaluation().
    // The coefficient arrays and state identity must remain unchanged.
    WVKernelStatus beginStateEvaluation(const WVState&);
    WVKernelStatus beginStateEvaluation(const WVState&, const void* evaluationOwner);
    // Register another immutable coefficient view at the same evaluation
    // times, retaining the already prepared phase factors.
    WVKernelStatus addStateEvaluationView(const WVState&, const void* evaluationOwner,
        std::size_t componentIdentity = 0);
    WVKernelStatus removeStateEvaluationView(const WVState&, const void* evaluationOwner,
        std::size_t componentIdentity);
    WVKernelStatus endStateEvaluation();
    bool stateEvaluationActive() const noexcept { return stateEvaluationActive_; }
    WVKernelStatus validateStateEvaluation(const WVState&) const noexcept;
    // Validate a flux target against the active (or standalone) state without
    // executing coefficient evolution or projection.
    WVKernelStatus validateFluxOutput(const WVState&, const WVFlux&);
    // Borrowed until the next kernel operation; scoped calls reuse the active phase.
    WVKernelStatus preparedPhase(const WVState&, WVComplexConstView&);

    WVKernelStatus transformToSpatial(WVComplexConstView, WVHydrostaticFamily, WVRealVolumeView);
    WVKernelStatus transformFromSpatial(WVRealVolumeConstView, WVHydrostaticFamily, WVComplexView);
    WVKernelStatus applyVertical(WVStratifiedModalOperator, WVComplexConstView, WVComplexView);
    WVKernelStatus applyVerticalColumn(WVStratifiedModalOperator, std::size_t retainedColumn,
        WVComplexConstView inputColumn, WVComplexView outputColumn);
    WVKernelStatus horizontalForward(WVRealVolumeConstView, WVComplexView);
    WVKernelStatus horizontalInverse(WVComplexConstView, WVRealVolumeView);
    WVKernelStatus transformUVEtaToWaveVortex(WVRealVolumeConstView u, WVRealVolumeConstView v,
        WVRealVolumeConstView eta, double t, double t0, WVMutableCoefficients);
    // Surface fields use [Nx,Ny,1]. Horizontal vorticity supports value only;
    // other fields support first derivatives. rhoTotal belongs to the full flow.
    WVKernelStatus transformStateField(const WVState&, WVHydrostaticField, WVRealVolumeView,
        WVHydrostaticDerivative = WVHydrostaticDerivative::value,
        WVHydrostaticComponent = WVHydrostaticComponent::all);
    // Reconstruct derived [u,v,eta] coefficient tendencies while retaining the
    // active primary state's matching phase and foreign-state protection.
    WVKernelStatus transformCoefficientTendencyToUVEta(
        const WVState& tendency, WVRealFieldBundleView& fields);
    // Combine exact prepared derivative operands without another reconstruction.
    // zetaX expects (w_y,v_z); zetaY expects (u_z,w_x). Either input may be
    // the exact output view; partial overlap remains invalid.
    WVKernelStatus combinePreparedHorizontalVorticity(WVHydrostaticField,
        WVRealVolumeConstView firstDerivative, WVRealVolumeConstView secondDerivative,
        WVRealVolumeView output,
        WVHydrostaticComponent = WVHydrostaticComponent::all);
    // Convert exact prepared eta_z and eta operands to rho_e,z or rho_total,z.
    WVKernelStatus combinePreparedDensityZDerivative(WVHydrostaticField,
        WVRealVolumeConstView etaZ, WVRealVolumeConstView eta,
        WVRealVolumeView output,
        WVHydrostaticComponent = WVHydrostaticComponent::all);
    // Linear phase evolution gives current-time coefficients; stored amplitudes
    // are stationary under the f-plane linear dynamics. Exact in-place allowed.
    WVKernelStatus evolveCoefficients(const WVState&, WVMutableCoefficients);
    // Remove inactive modes and enforce real A0 means / conjugate inertial Am.
    WVKernelStatus constrainCoefficients(WVMutableCoefficients) const;
    // Caller-owned observation output captures the raw spatial contribution.
    // Prepared [u,v,w,eta] fields may be shared within one observation event.
    // With projectFlux=false, spatialTendency is required and flux is untouched.
    // Produces the complete physical bundle and projected flux together. Caller
    // publishes evaluator nodes only after success; requires prepared native support.
    bool supportsTiledNonlinear() const noexcept { return tiledNonlinearPrepared_; }
    WVKernelStatus nonlinearFluxAndFields(const WVState&, WVFlux&, WVRealFieldBundleView,
        WVRealFieldBundleView* spatialTendency = nullptr);
    WVKernelStatus nonlinearFlux(const WVState&, WVFlux&,
        WVRealFieldBundleView* spatialTendency = nullptr,
        const WVRealFieldBundleConstView* preparedFields = nullptr, bool projectFlux = true,
        WVStateDerivativeAccess* derivativeAccess = nullptr);
    // Reduce supplied physical fields with the prepared workers. The caller's
    // evaluation owns caching; this operation retains no field or result.
    WVKernelStatus reduceHorizontalSpeedMaximum(WVRealVolumeConstView u,
        WVRealVolumeConstView v, double& maximum);
    WVKernelStatus totalEnergy(const WVCoefficients&, double&,
        WVHydrostaticComponent = WVHydrostaticComponent::all) const;
    WVKernelStatus totalEnstrophy(const WVCoefficients&, double&) const;
    WVKernelStatus totalEnergySpatiallyIntegrated(const WVState&, double&,
        WVHydrostaticComponent = WVHydrostaticComponent::all);
    WVKernelStatus applyVerticalCalculus(WVRealConstView input,bool inputIsF,
        unsigned order,bool integral,WVRealView output);
    // MATLAB diffZF/diffZG orders 1..4 and intZF/intZG order 1. These preserve
    // all horizontal grid columns rather than truncating through a Fourier map.
    WVKernelStatus differentiateVertical(WVRealVolumeConstView, WVHydrostaticFamily,
        unsigned order, WVRealVolumeView);
    // Full-grid horizontal derivatives and three-dimensional passive advection.
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView, bool xDerivative, WVRealVolumeView);
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView, bool xDerivative,
        unsigned order, WVRealVolumeView);
    WVKernelStatus advectScalarWithAdvectionFields(WVRealVolumeConstView, WVRealFieldBundleConstView, bool antialias, WVRealVolumeView, bool xyOnly = false);
    WVKernelStatus integrateVertical(WVRealVolumeConstView, WVHydrostaticFamily, WVRealVolumeView);
private:
    WVTransformHydrostaticKernel() = default;
    WVKernelStatus spectral(WVComplexConstView) const;
    WVKernelStatus volume(WVRealVolumeConstView, bool surface = false) const;
    WVKernelStatus coefficients(const WVCoefficients&) const;
    WVKernelStatus outputs(WVMutableCoefficients) const;
    WVKernelStatus mutableOutputOutsidePreparedState(WVMutableCoefficients) const;
    WVKernelStatus mutableOutputOutsidePreparedState(const void*,std::size_t) const;
    WVKernelStatus stateContents(const WVState&) const;
    WVKernelStatus state(const WVState&);
    bool matchesStateEvaluation(const WVState&) const noexcept;
    std::size_t stateEvaluationComponent(const WVState&) const noexcept;
    std::size_t fieldPreparationView(const WVCoefficients&) const noexcept;
    WVKernelStatus validateStateForCall(const WVState&);
    WVKernelStatus preparePhaseForCall(const WVState&);
    WVKernelStatus prepareProjectionPhaseForCall(double t,double t0);
    WVKernelStatus disjoint(const void*,std::size_t,const void*,std::size_t) const;
    WVKernelStatus preparePhase(double t,double t0,const WVState* validatedState = nullptr);
    WVComplexOutput modalView(std::size_t slot = 0);
    WVComplexOutput gridView(std::size_t slot = 0);
    WVKernelStatus vertical(std::size_t,WVComplexInput,WVComplexOutput);
    WVKernelStatus verticalColumn(std::size_t,WVComplexInput,WVComplexOutput,std::size_t);
    WVKernelStatus project(const double*,WVComplexOutput,WVHydrostaticFamily);
    WVKernelStatus reconstruct(const WVCoefficients&,WVHydrostaticField,
        WVHydrostaticDerivative,WVHydrostaticComponent,double*,bool countPrimary = true,
        std::size_t metricComponent = 5,const WVRealOutputConsumer* consumer = nullptr,
        WVComplexInput* preparedSpectrum = nullptr,const WVComplexOutput* destination = nullptr);
    WVKernelStatus projectFields(const double*,const double*,const double*,WVMutableCoefficients);
    WVKernelStatus projectedFieldsToCoefficients(
        WVComplexInput,WVComplexInput,WVComplexInput,WVMutableCoefficients);
    WVKernelStatus verticalCalculus(const double*,WVHydrostaticFamily,unsigned,bool,double*,
        std::size_t columns,bool verticalFirst);
    bool tiledNonlinearPrepared_ = false;
    WVVariableExecutionOptions executionOptions_;
    std::shared_ptr<const WVStratifiedModalSource> source_;
    std::vector<WVHydrostaticModeFactors> factors_;
    mutable WVHydrostaticStorage storage_;
    std::size_t baseSpectralScratchBytes_ = 0;
    WVHydrostaticKernelMetrics metrics_;
    std::string engineIdentifier_,engineLibraryIdentity_;
    std::unique_ptr<WVRetainedHorizontalOperator> horizontal_;
    std::unique_ptr<WVRetainedHorizontalWorkspace> horizontalWorkspace_;
    std::array<std::unique_ptr<WVPreparedVerticalOperator>,4> vertical_;
    std::array<std::unique_ptr<WVVerticalWorkspace>,4> verticalWorkspace_;
    std::unique_ptr<spectral_detail::WVVariableComplexBuffer> spectralStorage_;
    std::unique_ptr<kernel_detail::WVPreparedModeExecutor> pointwise_;
    std::unique_ptr<kernel_detail::WVPreparedFieldCache> fieldCache_;
    std::vector<WVComplex64> phase_;
    std::vector<double> real_;
    WVState preparedState_{};
    std::array<WVState,5> preparedStateViews_{};
    std::array<std::size_t,5> preparedStateComponents_{};
    std::size_t preparedStateViewCount_ = 0;
    std::array<std::size_t,5> preparedStateViewIds_{};
    std::size_t nextStateViewId_ = 1;
    const void* preparedStateOwner_ = nullptr;
    bool stateEvaluationActive_ = false;
    std::size_t S_ = 0,R_ = 0,H_ = 0;
    std::atomic<bool> active_{false};
};
} // namespace wavevortex
