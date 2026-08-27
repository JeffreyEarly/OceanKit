#pragma once

#include "WVFFTEngine.hpp"

#include <memory>
#include <vector>

namespace wavevortex {

enum class WVDynamicalField : std::uint8_t { u, v, w, eta };

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
    static WVKernelStatus create(const WVTransformConstantStratificationConfiguration& configuration, std::unique_ptr<WVFFTEngine> engine, std::unique_ptr<WVTransformConstantStratificationKernel>& kernel);

    ~WVTransformConstantStratificationKernel() = default;
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
    std::size_t phaseReservationBytes() const noexcept;
    void setStageInstrumentation(bool enabled) noexcept;
    std::size_t persistentBytes() const noexcept;
    std::size_t scratchBytes() const noexcept { return (halfSpectrumScratch_.size() + realScratch_.size()) * sizeof(double); }

    WVKernelStatus transformUVEtaToWaveVortex(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients);
    WVKernelStatus transformUVWEtaToWaveVortex(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients);
    WVKernelStatus transformWaveVortexToUVWEta(const WVState& state, WVRealFieldBundleView& fields);
    WVKernelStatus transformWaveVortexToUVW(const WVState& state, WVRealFieldBundleView& fields);
    WVKernelStatus transformStateFieldDerivatives(const WVState& state, WVDynamicalField field, WVRealFieldBundleView& derivatives);
    WVKernelStatus transformToSpatialDomainWithFAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields);
    WVKernelStatus transformToSpatialDomainWithGAllDerivatives(const WVComplexConstView& Apm, const WVComplexConstView& A0, WVRealFieldBundleView& fields);
    WVKernelStatus nonlinearFlux(const WVState& state, WVFlux& flux);
    WVKernelStatus nonlinearFluxWithAdvectionFields(const WVState& state, WVFlux& flux, WVRealFieldBundleView& advectionFields);
    WVKernelStatus nonlinearFluxUsingAdvectionFields(const WVState& state, WVFlux& flux, const WVRealFieldBundleConstView& advectionFields);
    WVKernelStatus prepareScalarAdvection();
    WVKernelStatus advectFGridScalar(const WVRealVolumeConstView& scalar, const WVRealFieldBundleConstView& advectionFields, bool shouldAntialias, WVRealVolumeView& rightHandSide);

private:
    WVTransformConstantStratificationKernel() = default;
    WVKernelStatus preparePlans();
    WVKernelStatus transformUVEtaToWaveVortexImpl(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues = {});
    WVKernelStatus transformUVWEtaToWaveVortexImpl(const WVRealFieldBundleConstView& fields, double t, double t0, WVMutableCoefficients& coefficients, WVComplexConstView phaseValues = {});
    WVKernelStatus transformWaveVortexToUVWEtaImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients = nullptr);
    WVKernelStatus transformWaveVortexToUVWImpl(const WVState& state, WVRealFieldBundleView& fields, const WVCoefficients* evolvedCoefficients, WVComplexConstView phaseValues = {});
    WVKernelStatus transformToSpatialDomainWithDerivativesImpl(const WVCoefficients& evolvedCoefficients, std::size_t target, WVRealFieldBundleView& derivatives);
    WVKernelStatus transformToSpatialDomainWithDerivativesFromStateImpl(const WVState& state, WVComplexConstView phaseValues, std::size_t target, WVRealFieldBundleView& derivatives);
    WVKernelStatus projectSingleFluxTargetImpl(const WVRealFieldBundleConstView& field, std::size_t target, WVComplexConstView phaseValues, WVFlux& flux);
    WVKernelStatus nonlinearFluxImpl(const WVState& state, WVFlux& flux, WVRealFieldBundleView* advectionFields, bool advectionFieldsPrepared = false);
    WVKernelStatus ensureScalarInversePlan();
    WVKernelStatus antialiasScalarInPlace(WVRealVolumeView& scalar);
    WVTransformConstantStratificationDescriptor descriptor_;
    std::unique_ptr<WVFFTEngine> engine_;
    std::string engineIdentifier_;
    std::string engineLibraryIdentity_;
    std::vector<std::unique_ptr<WVFFTPlan>> plans_;
    std::unique_ptr<WVFFTPlan> scalarInversePlan_;
    std::vector<std::uint8_t> scalarAntialiasRows_;
    std::vector<double> halfSpectrumScratch_;
    std::vector<double> realScratch_;
    mutable WVKernelMetrics metrics_;
    bool executing_ = false;
    bool stageInstrumentationEnabled_ = false;
};

} // namespace wavevortex
