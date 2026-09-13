#pragma once
#include "WVStratifiedModalSource.hpp"
#include "WVVariableExecutionOptions.hpp"
#include <array>
#include <atomic>
#include <functional>

namespace wavevortex {
namespace spectral_detail { class WVVariableComplexBuffer; }
namespace kernel_detail { class WVPreparedModeExecutor; }
inline constexpr const char* WVStratifiedQGKernelContract = "wave-vortex-stratified-qg-kernel-v1";
enum class WVStratifiedQGField { u, v, w, eta, pi, p, psi, qgpv, rhoE, rhoTotal, zetaZ, ssh, ssu, ssv };
enum class WVStratifiedQGDerivative { value, x, y, z };
struct WVStratifiedQGModeFactors {
    double f = 0, beta = 0;
    std::vector<WVComplex64> u, v;
    std::vector<double> eta, pi, psi, qgpv, zetaZ, energy, enstrophy, kineticEnergy, potentialEnergy;
};
struct WVStratifiedQGStorage {
    std::size_t sharedScientificBytes = 0, preparedBytes = 0, workspaceBytes = 0;
    std::size_t spectralScratchBytes = 0, realScratchBytes = 0, factorBytes = 0;
    std::size_t providerBytesLowerBound = 0, planBytesLowerBound = 0;
};
struct WVStratifiedQGKernelMetrics {
    std::size_t stateValidationCount = 0;
    std::size_t horizontalSpeedMaximumReductionCount = 0;
    std::array<std::size_t,14> fieldReconstructionCount{};
    std::array<std::array<std::size_t,4>,14> reconstructionCount{};
    std::array<std::array<std::array<std::size_t,5>,4>,14> componentReconstructionCount{};
};

// One owned mutable workspace, canonical interleaved [Nj,Nkl] coefficients and
// column-major [Nx,Ny,Nz] fields. No wave coefficient arrays or eigenproblem.
// Calls on one kernel must not overlap; create separate kernels for concurrency.
class WVTransformStratifiedQGKernel final {
public:
    using MatrixBackendFactory = std::function<WVKernelStatus(std::unique_ptr<WVVerticalMatrixBackend>&)>;
    static WVKernelStatus create(std::shared_ptr<const WVStratifiedModalSource>,
        std::unique_ptr<WVFFTEngine>, std::unique_ptr<WVTransformStratifiedQGKernel>&,
        MatrixBackendFactory = WVCreateScalarMatrixBackend, WVVariableExecutionOptions = {});
    ~WVTransformStratifiedQGKernel();
    const WVVariableExecutionOptions& executionOptions() const noexcept { return executionOptions_; }
    const char* horizontalScheduleIdentifier() const noexcept { return horizontalWorkspace_->scheduleIdentifier(); }
    const WVStratifiedModalGeometry& geometry() const noexcept { return source_->geometry(); }
    const WVStratifiedQGModeFactors& factors() const noexcept { return factors_; }
    const WVStratifiedQGStorage& storage() const noexcept { return storage_; }
    const WVStratifiedQGKernelMetrics& metrics() const noexcept { return metrics_; }
    void resetMetrics() noexcept { metrics_ = {}; }
    const char* matrixBackendIdentifier() const noexcept { return vertical_[0]->backendIdentifier(); }
    std::size_t persistentBytes() const noexcept {
        return sizeof(*this)+storage_.sharedScientificBytes+storage_.preparedBytes+storage_.workspaceBytes+storage_.spectralScratchBytes+storage_.realScratchBytes+storage_.factorBytes;
    }
    const std::string& engineIdentifier() const noexcept { return engineIdentifier_; }
    WVShape2D spectralShape() const noexcept { return {geometry().Nj,geometry().Nkl}; }
    WVShape3D spatialShape() const noexcept { return {geometry().Nx,geometry().Ny,geometry().Nz}; }
    // Prepare the additional spectral slot used only by MATLAB vertical calculus.
    WVKernelStatus prepareMatlabPrimitives();

    // The borrowed A0 array remains immutable until endStateEvaluation().
    WVKernelStatus beginStateEvaluation(WVComplexConstView);
    WVKernelStatus beginStateEvaluation(WVComplexConstView, const void* evaluationOwner);
    WVKernelStatus addStateEvaluationView(WVComplexConstView,
        const void* evaluationOwner, std::size_t componentIdentity = 0);
    WVKernelStatus removeStateEvaluationView(WVComplexConstView,
        const void* evaluationOwner, std::size_t componentIdentity);
    WVKernelStatus endStateEvaluation();
    bool stateEvaluationActive() const noexcept { return stateEvaluationActive_; }
    WVKernelStatus validateStateEvaluation(WVComplexConstView) const noexcept;

    // Projection preserves the horizontal mean just as MATLAB's raw transform
    // does. Reconstructed QG fields mask all horizontal means as MATLAB does.
    WVKernelStatus transformQGPVToA0(WVRealVolumeConstView, WVComplexView);
    WVKernelStatus applyVertical(WVStratifiedModalOperator, WVComplexConstView, WVComplexView);
    WVKernelStatus applyVerticalColumn(WVStratifiedModalOperator, std::size_t retainedColumn,
        WVComplexConstView inputColumn, WVComplexView outputColumn);
    WVKernelStatus horizontalForward(WVRealVolumeConstView, WVComplexView);
    WVKernelStatus horizontalInverse(WVComplexConstView, WVRealVolumeView);
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView, bool xDerivative,
        WVRealVolumeView);
    WVKernelStatus differentiateHorizontal(WVRealVolumeConstView, bool xDerivative,
        unsigned order, WVRealVolumeView);
    // Diagnostic modal inverse preserves the mean and projects self-conjugate
    // Fourier values to their real part; ordinary field masks are unchanged.
    WVKernelStatus transformSpectralTendencyToSpatial(WVComplexConstView, WVRealVolumeView);
    WVKernelStatus transformUVEtaToA0(WVRealVolumeConstView u, WVRealVolumeConstView v,
        WVRealVolumeConstView eta, WVComplexView);
    WVKernelStatus transformA0ToField(WVComplexConstView, WVStratifiedQGField,
        WVRealVolumeView, WVStratifiedQGDerivative = WVStratifiedQGDerivative::value);
    // Convert exact prepared eta_z and eta operands to rho_e,z or rho_total,z.
    WVKernelStatus combinePreparedDensityZDerivative(WVStratifiedQGField,
        WVRealVolumeConstView etaZ, WVRealVolumeConstView eta,
        WVRealVolumeView output, std::size_t componentIdentity = 0);
    // Surface fields use [Nx,Ny,1]; other fields use [Nx,Ny,Nz].
    // Optional raw output captures the spatial contribution instead of projecting
    // it. Prepared [u,v] fields may be borrowed for one diagnostic invocation.
    WVKernelStatus nonlinearFlux(WVComplexConstView, WVComplexView, double beta = 0,
        WVRealVolumeView* raw = nullptr, const WVRealFieldBundleConstView* preparedUV = nullptr);
    // Configured physical closures normally overwrite A0 output. With raw
    // output they capture before projection and leave A0 output untouched.
    // Both paths reuse the prepared operators and mutable workspace.
    WVKernelStatus verticalDiffusivityFlux(WVComplexConstView, double kappaZ, WVComplexView, WVRealVolumeView* raw = nullptr);
    WVKernelStatus linearBottomFrictionFlux(WVComplexConstView, double rate, WVComplexView, WVRealVolumeView* raw = nullptr);
    WVKernelStatus quadraticBottomFrictionFlux(WVComplexConstView, double dragCoefficient, WVComplexView, WVRealVolumeView* raw = nullptr,
        const WVRealFieldBundleConstView* preparedUV = nullptr);
    WVKernelStatus linearFlux(WVComplexConstView, WVComplexView, double beta = 0) const;
    // F-plane A0 is stationary. An explicit beta gives exact linear Rossby evolution.
    WVKernelStatus evolveA0(WVComplexConstView, double elapsedTime, WVComplexView, double beta = 0) const;
    WVKernelStatus totalEnergy(WVComplexConstView, double&) const;
    WVKernelStatus totalEnstrophy(WVComplexConstView, double&) const;
    WVKernelStatus totalEnergySpatiallyIntegrated(WVComplexConstView, double&);
    WVKernelStatus totalEnstrophySpatiallyIntegrated(WVComplexConstView, double&);
    WVKernelStatus applyVerticalCalculus(WVRealConstView input,bool inputIsF,
        unsigned order,bool integral,WVRealView output);
    const std::string& engineLibraryIdentity() const noexcept { return engineLibraryIdentity_; }
    WVKernelStatus advectScalarWithAdvectionFields(WVRealVolumeConstView, WVRealFieldBundleConstView, bool antialias, WVRealVolumeView);
    WVKernelStatus uvMax(WVComplexConstView, double&);
private:
    WVTransformStratifiedQGKernel() = default;
    WVKernelStatus spectral(WVComplexConstView) const;
    WVKernelStatus validateState(WVComplexConstView) const;
    WVKernelStatus validateStateForCall(WVComplexConstView) const;
    WVKernelStatus mutableOutputOutsidePreparedState(WVComplexView) const;
    WVKernelStatus mutableOutputOutsidePreparedState(const void*,std::size_t) const;
    bool matchesStateEvaluation(WVComplexConstView) const noexcept;
    std::size_t stateEvaluationComponent(WVComplexConstView) const noexcept;
    WVKernelStatus volume(WVRealVolumeConstView, bool surface = false) const;
    WVKernelStatus disjoint(const void*, std::size_t, const void*, std::size_t) const;
    WVKernelStatus validateDiagnosticBuffers(WVComplexConstView, WVComplexView,
        const WVRealVolumeView*, const WVRealFieldBundleConstView* = nullptr) const;
    WVKernelStatus project(const double*, WVComplexOutput, std::size_t operation = 1);
    WVKernelStatus reconstruct(WVComplexConstView, WVStratifiedQGField, WVStratifiedQGDerivative, double*);
    WVKernelStatus vertical(std::size_t operation, WVComplexInput, WVComplexOutput);
    WVKernelStatus verticalColumn(std::size_t operation, WVComplexInput, WVComplexOutput,
        std::size_t retainedColumn);
    WVKernelStatus verticalCalculus(const double*,bool inputIsF,unsigned,bool,double*,
        std::size_t columns,bool verticalFirst);
    WVVariableExecutionOptions executionOptions_;
    std::shared_ptr<const WVStratifiedModalSource> source_;
    WVStratifiedQGModeFactors factors_;
    WVStratifiedQGStorage storage_;
    mutable WVStratifiedQGKernelMetrics metrics_;
    std::string engineIdentifier_, engineLibraryIdentity_;
    std::unique_ptr<WVRetainedHorizontalOperator> horizontal_;
    std::unique_ptr<WVRetainedHorizontalWorkspace> horizontalWorkspace_;
    std::array<std::unique_ptr<WVPreparedVerticalOperator>,4> vertical_;
    std::array<std::unique_ptr<WVVerticalWorkspace>,4> verticalWorkspace_;
    std::unique_ptr<spectral_detail::WVVariableComplexBuffer> modalSpectral_, gridSpectral_;
    std::unique_ptr<kernel_detail::WVPreparedModeExecutor> pointwise_;
    std::vector<double> real_;
    WVComplexConstView preparedState_{};
    std::array<WVComplexConstView,5> preparedStateViews_{};
    std::array<std::size_t,5> preparedStateComponents_{};
    std::size_t preparedStateViewCount_ = 0;
    const void* preparedStateOwner_ = nullptr;
    bool stateEvaluationActive_ = false;
    std::size_t S_ = 0, R_ = 0, H_ = 0;
    std::atomic<bool> active_{false};
};
} // namespace wavevortex
