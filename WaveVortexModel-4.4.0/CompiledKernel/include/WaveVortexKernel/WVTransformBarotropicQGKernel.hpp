#pragma once

#include "WVFFTEngine.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex {

struct WVTransformBarotropicQGConfiguration {
    std::uint32_t contractVersion = WVKernelContractVersion;
    std::size_t Nx = 0;
    std::size_t Ny = 0;
    double Lx = 0.0;
    double Ly = 0.0;
    double h = 0.0;
    std::uint32_t j = 1;
    double g = 0.0;
    double planetaryRadius = 0.0;
    double rotationRate = 0.0;
    double latitude = 0.0;
    bool shouldAntialias = true;
};

bool sameTransformConfiguration(
    const WVTransformBarotropicQGConfiguration& first,
    const WVTransformBarotropicQGConfiguration& second) noexcept;

struct WVBarotropicQGModes {
    double coriolisFrequency = 0.0;
    double deformationWavenumberSquared = 0.0;
    std::vector<WVComplex64> uFactor;
    std::vector<WVComplex64> vFactor;
    std::vector<double> etaFactor;
    std::vector<double> piFactor;
    std::vector<double> psiFactor;
    std::vector<double> qgpvFactor;
    std::vector<double> zetaZFactor;
    std::vector<double> energyFactor;
    std::vector<double> enstrophyFactor;
    std::size_t persistentBytes() const noexcept;
};

class WVTransformBarotropicQGDescriptor final {
public:
    static WVKernelStatus create(
        const WVTransformBarotropicQGConfiguration& configuration,
        WVTransformBarotropicQGDescriptor& descriptor);

    const WVTransformBarotropicQGConfiguration& configuration() const noexcept {
        return configuration_;
    }
    const std::vector<WVFourierMode>& fourierModes() const noexcept {
        return fourierModes_;
    }
    const WVHalfSpectrumMappings& halfSpectrumMappings() const noexcept {
        return halfSpectrumMappings_;
    }
    const WVBarotropicQGModes& modes() const noexcept { return modes_; }
    std::size_t Nkl() const noexcept { return fourierModes_.size(); }
    WVShape2D spectralShape() const noexcept { return {1, Nkl()}; }
    WVShape2D spatialShape() const noexcept {
        return {configuration_.Nx, configuration_.Ny};
    }
    std::size_t persistentBytes() const noexcept;

private:
    WVTransformBarotropicQGConfiguration configuration_;
    std::vector<WVFourierMode> fourierModes_;
    WVHalfSpectrumMappings halfSpectrumMappings_;
    WVBarotropicQGModes modes_;
};

enum class WVBarotropicQGField : std::uint8_t {
    u,
    v,
    eta,
    pi,
    psi,
    qgpv,
    zetaZ,
    ssh
};

struct WVBarotropicQGKernelMetrics {
    std::size_t descriptorBytes = 0;
    std::size_t planBytes = 0;
    std::size_t engineBytes = 0;
    std::size_t kernelManagementBytes = 0;
    std::size_t halfSpectrumScratchCapacityBytes = 0;
    std::size_t realScratchCapacityBytes = 0;
    std::size_t scratchCapacityBytes = 0;
    std::size_t scratchHighWaterBytes = 0;
    std::size_t persistentFullHermitianBytes = 0;
    std::size_t planCount = 0;
    std::size_t executionCount = 0;
    std::size_t forwardExecutionCount = 0;
    std::size_t inverseExecutionCount = 0;
    std::size_t fieldEvaluationCount = 0;
    std::size_t derivativeEvaluationCount = 0;
    std::size_t nonlinearFluxCallCount = 0;
    std::size_t antialiasedNonlinearFluxCallCount = 0;
    std::size_t forcingFieldReconstructionCount = 0;
    std::size_t forcingFieldReuseCount = 0;
    std::size_t forcingSpatialProjectionCount = 0;
    std::size_t horizontalSpeedMaximumReductionCount = 0;
    std::size_t bytesCopied = 0;
    std::size_t stateValidationCount = 0;
    std::array<std::array<std::size_t,3>,8> reconstructionCount{};
    std::array<std::array<std::array<std::size_t,5>,3>,8> componentReconstructionCount{};
};

// RHS-scoped transform workspace state. The numerical storage remains owned
// by the kernel's bounded 4H+5R scratch. This record tracks valid fields
// between coarse operations and can borrow caller-owned diagnostic output.
struct WVBarotropicQGOperationWorkspace {
    // Optional invocation-owned raw output. Spatial operations capture here
    // instead of projecting; the caller projects the cumulative stage once.
    WVRealView* spatialTendency = nullptr;
    bool spatialTendencyCaptured = false;
    const WVRealFieldBundleConstView* preparedVelocity = nullptr;
    bool physicalFieldsPrepared = false;
    bool qgpvDerivativesPrepared = false;
    bool horizontalSpeedMaximumPrepared = false;
    double horizontalSpeedMaximum = 0.0;
    std::size_t physicalFieldReconstructionCount = 0;
    std::size_t physicalFieldReuseCount = 0;
    std::size_t spatialTendencyProjectionCount = 0;
};

class WVTransformBarotropicQGKernel final {
public:
    static WVKernelStatus create(
        const WVTransformBarotropicQGConfiguration& configuration,
        std::unique_ptr<WVFFTEngine> engine,
        std::unique_ptr<WVTransformBarotropicQGKernel>& kernel);

    ~WVTransformBarotropicQGKernel() = default;
    WVTransformBarotropicQGKernel(const WVTransformBarotropicQGKernel&) = delete;
    WVTransformBarotropicQGKernel& operator=(const WVTransformBarotropicQGKernel&) = delete;
    WVTransformBarotropicQGKernel(WVTransformBarotropicQGKernel&&) = delete;
    WVTransformBarotropicQGKernel& operator=(WVTransformBarotropicQGKernel&&) = delete;

    const WVTransformBarotropicQGDescriptor& descriptor() const noexcept {
        return descriptor_;
    }
    const WVBarotropicQGKernelMetrics& metrics() const noexcept;
    const std::string& engineIdentifier() const noexcept {
        return engineIdentifier_;
    }
    const std::string& engineLibraryIdentity() const noexcept {
        return engineLibraryIdentity_;
    }
    const char* coefficientOrderingIdentifier() const noexcept {
        return "matlab-kl-radial-k-l";
    }
    const char* normalizationIdentifier() const noexcept {
        return "A0-is-qgpv";
    }
    const char* antialiasImplementationIdentifier() const noexcept {
        return "transform-level-radial-two-thirds";
    }
    std::size_t persistentBytes() const noexcept;
    std::size_t scratchBytes() const noexcept;

    // The borrowed A0 array remains immutable until endStateEvaluation().
    WVKernelStatus beginStateEvaluation(const WVComplexConstView& A0);
    WVKernelStatus beginStateEvaluation(const WVComplexConstView& A0, const void* evaluationOwner);
    WVKernelStatus addStateEvaluationView(const WVComplexConstView& A0,
        const void* evaluationOwner, std::size_t componentIdentity = 0);
    WVKernelStatus removeStateEvaluationView(const WVComplexConstView& A0,
        const void* evaluationOwner, std::size_t componentIdentity);
    WVKernelStatus endStateEvaluation();
    bool stateEvaluationActive() const noexcept { return stateEvaluationActive_; }
    WVKernelStatus validateStateEvaluation(const WVComplexConstView& A0) const noexcept;

    WVKernelStatus transformQGPVToA0(const WVRealConstView& qgpv,
                                     WVComplexView& A0);
    WVKernelStatus transformA0ToQGPV(const WVComplexConstView& A0,
                                     WVRealView& qgpv);
    // Raw retained Fourier transforms preserve the horizontal mean.
    WVKernelStatus horizontalForward(const WVRealConstView&, WVComplexView&);
    WVKernelStatus horizontalInverse(const WVComplexConstView&, WVRealView&);
    WVKernelStatus differentiateHorizontal(const WVRealConstView&,
        bool xDerivative, WVRealView&);
    WVKernelStatus differentiateHorizontal(const WVRealConstView&,
        bool xDerivative, unsigned order, WVRealView&);
    // Raw Fourier inverse for diagnostic contributions, including the mean.
    // Ordinary model fields retain their existing mean mask.
    WVKernelStatus transformSpectralTendencyToSpatial(const WVComplexConstView&, WVRealView&);
    WVKernelStatus transformA0ToField(const WVComplexConstView& A0,
                                     WVBarotropicQGField field,
                                     WVRealView& output);
    // Output channels are field, x derivative, and y derivative in that order.
    WVKernelStatus transformA0ToFieldWithDerivatives(
        const WVComplexConstView& A0, WVBarotropicQGField field,
        WVRealFieldBundleView& output);
    WVKernelStatus evolveA0(const WVComplexConstView& A0, double elapsedTime,
                            WVComplexView& evolvedA0) const;
    WVKernelStatus nonlinearFlux(const WVComplexConstView& A0,
                                 WVComplexView& F0);
    WVKernelStatus addPotentialVorticityAdvection(
        const WVComplexConstView& A0, WVComplexView& F0, bool accumulate,
        WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus addAdaptiveDamping(
        const WVComplexConstView& A0,
        const std::vector<double>& dampingOperator, WVComplexView& F0,
        bool accumulate, WVBarotropicQGOperationWorkspace& workspace);
    // Reuse the operation workspace's prepared u/v fields and retain the
    // reduction for subsequent adaptive forcing in the same RHS evaluation.
    WVKernelStatus horizontalSpeedMaximum(
        const WVComplexConstView& A0, double& maximumSpeed,
        WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus addLinearBottomFriction(
        const WVComplexConstView& A0, double rate, WVComplexView& F0,
        bool accumulate, WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus addQuadraticBottomFriction(
        const WVComplexConstView& A0, double drag, WVComplexView& F0,
        bool accumulate, WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus addBetaPlanePVAdvection(
        const WVComplexConstView& A0, double beta, WVComplexView& F0,
        bool accumulate, WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus totalEnergy(const WVComplexConstView& A0,
                               double& energy) const;
    WVKernelStatus totalEnstrophy(const WVComplexConstView& A0,
                                  double& enstrophy) const;
    WVKernelStatus totalEnergySpatiallyIntegrated(
        const WVComplexConstView& A0, double& energy);
    WVKernelStatus totalEnstrophySpatiallyIntegrated(
        const WVComplexConstView& A0, double& enstrophy);
    WVKernelStatus uvMax(const WVComplexConstView& A0,
                         double& maximumSpeed);
    WVKernelStatus advectScalar(const WVComplexConstView& A0,
                                const WVRealConstView& scalar,
                                bool shouldAntialias,
                                WVRealView& rightHandSide);
    WVKernelStatus prepareAdvectionFields(
        const WVComplexConstView& A0,
        WVBarotropicQGOperationWorkspace& workspace,
        WVRealFieldBundleConstView& fields);
    WVKernelStatus advectScalarWithAdvectionFields(
        const WVRealConstView& scalar,
        const WVRealFieldBundleConstView& advectionFields,
        bool shouldAntialias, WVRealView& rightHandSide);
    std::size_t enforceReality(WVComplexView& A0) const noexcept;
    WVKernelStatus constrainA0(WVComplexView& A0, std::size_t& modified) const;

private:
    WVTransformBarotropicQGKernel() = default;
    WVKernelStatus validateState(const WVComplexConstView&) const;
    WVKernelStatus validateStateForCall(const WVComplexConstView&) const;
    WVKernelStatus mutableOutputOutsidePreparedState(const WVComplexView&) const;
    WVKernelStatus mutableOutputOutsidePreparedState(const void*,std::size_t) const;
    bool matchesStateEvaluation(const WVComplexConstView&) const noexcept;
    std::size_t stateEvaluationComponent(const WVComplexConstView&) const noexcept;
    void recordReconstruction(const WVComplexConstView&,WVBarotropicQGField,
        std::size_t derivative) const noexcept;
    WVKernelStatus preparePlans();
    WVKernelStatus forward(const WVRealConstView& input, WVComplexView& output,
                           bool accumulate = false);
    WVKernelStatus inverse(const WVComplexConstView& input,
                           const WVComplex64* factors, WVRealView& output);
    WVKernelStatus inverse(const WVComplexConstView& input,
                           const double* factors, WVRealView& output);
    WVKernelStatus inverseNonlinearFields(const WVComplexConstView& A0);
    WVKernelStatus ensureForcingFields(
        const WVComplexConstView& A0, bool requireQGPVDerivatives,
        WVBarotropicQGOperationWorkspace& workspace);
    WVKernelStatus ensureHorizontalSpeedMaximum(
        const WVComplexConstView& A0,
        WVBarotropicQGOperationWorkspace& workspace);
    double reduceHorizontalSpeedMaximum() const noexcept;
    WVKernelStatus inverseQGPVDerivative(const WVComplexConstView& A0,
                                         bool xDerivative,
                                         WVRealView& output);
    WVKernelStatus spatialDerivative(const WVRealConstView& input,
                                     bool xDerivative,
                                     WVRealView& output);
    WVKernelStatus spatialDerivative(const WVRealConstView& input,
                                     bool xDerivative, unsigned order,
                                     WVRealView& output);
    WVKernelStatus antialiasScalarInPlace(WVRealView& scalar);
    WVKernelStatus validateForcingOperation(
        const WVComplexConstView& A0, const WVComplexView& F0,
        const WVRealView* spatialTendency = nullptr,
        const WVRealFieldBundleConstView* preparedVelocity = nullptr) const;
    WVKernelStatus fillHalfSpectrum(const WVComplexConstView& input,
                                    const WVComplex64* factors,
                                    std::size_t field, std::size_t fields);
    WVKernelStatus fillHalfSpectrum(const WVComplexConstView& input,
                                    const double* factors,
                                    std::size_t field, std::size_t fields);
    const WVComplex64* complexFactors(WVBarotropicQGField field) const noexcept;
    const double* realFactors(WVBarotropicQGField field) const noexcept;

    WVTransformBarotropicQGDescriptor descriptor_;
    std::unique_ptr<WVFFTEngine> engine_;
    std::string engineIdentifier_;
    std::string engineLibraryIdentity_;
    std::vector<std::unique_ptr<WVFFTPlan>> plans_;
    std::vector<WVComplex64> halfSpectrumScratch_;
    std::vector<double> realScratch_;
    mutable WVBarotropicQGKernelMetrics metrics_;
    WVComplexConstView preparedState_{};
    std::array<WVComplexConstView,5> preparedStateViews_{};
    std::array<std::size_t,5> preparedStateComponents_{};
    std::size_t preparedStateViewCount_ = 0;
    const void* preparedStateOwner_ = nullptr;
    bool stateEvaluationActive_ = false;
    bool executing_ = false;
};

} // namespace wavevortex
