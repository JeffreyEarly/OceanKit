#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

namespace wavevortex {

inline constexpr std::uint32_t WVKernelContractVersion = 4;

struct WVComplex64 {
    double real = 0.0;
    double imag = 0.0;
};

static_assert(std::is_standard_layout<WVComplex64>::value, "WVComplex64 must be standard layout");
static_assert(std::is_trivially_copyable<WVComplex64>::value, "WVComplex64 must be trivially copyable");
static_assert(sizeof(WVComplex64) == 2 * sizeof(double), "WVComplex64 must contain two doubles without padding");

enum class WVKernelStatusCode : std::uint32_t {
    success = 0,
    invalidConfiguration,
    invalidShape,
    invalidPointer,
    overlappingArrays,
    sizeOverflow,
    allocationFailure,
    fftPlanFailure,
    fftExecutionFailure,
    numericalFailure,
    unsupportedOperation,
    reentrantExecution
};

struct WVKernelStatus {
    WVKernelStatusCode code = WVKernelStatusCode::success;
    std::string message;

    explicit operator bool() const noexcept { return code == WVKernelStatusCode::success; }
    static WVKernelStatus ok() { return {}; }
};

struct WVShape2D {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::size_t elementCount() const;
};

struct WVShape3D {
    std::size_t first = 0;
    std::size_t second = 0;
    std::size_t third = 0;
    std::size_t elementCount() const;
};

struct WVShape4D {
    std::size_t first = 0;
    std::size_t second = 0;
    std::size_t third = 0;
    std::size_t fourth = 0;
    std::size_t elementCount() const;
};

template <typename T>
struct WVMatrixView {
    T* data = nullptr;
    WVShape2D shape;
    bool empty() const noexcept { return shape.rows == 0 || shape.columns == 0; }
};

template <typename T>
struct WVFieldBundleView {
    T* data = nullptr;
    WVShape4D shape;
    bool empty() const noexcept { return shape.first == 0 || shape.second == 0 || shape.third == 0 || shape.fourth == 0; }
};

template <typename T>
struct WVVolumeView {
    T* data = nullptr;
    WVShape3D shape;
    bool empty() const noexcept { return shape.first == 0 || shape.second == 0 || shape.third == 0; }
};

using WVRealConstView = WVMatrixView<const double>;
using WVRealView = WVMatrixView<double>;
using WVRealVolumeConstView = WVVolumeView<const double>;
using WVRealVolumeView = WVVolumeView<double>;
using WVComplexConstView = WVMatrixView<const WVComplex64>;
using WVComplexView = WVMatrixView<WVComplex64>;
using WVRealFieldBundleConstView = WVFieldBundleView<const double>;
using WVRealFieldBundleView = WVFieldBundleView<double>;

struct WVCoefficients {
    WVComplexConstView Ap;
    WVComplexConstView Am;
    WVComplexConstView A0;
};

struct WVMutableCoefficients {
    WVComplexView Ap;
    WVComplexView Am;
    WVComplexView A0;
};

struct WVState {
    double t = 0.0;
    double t0 = 0.0;
    WVCoefficients coefficients;
};

struct WVMutableState {
    double t = 0.0;
    double t0 = 0.0;
    WVMutableCoefficients coefficients;

    WVState view() const noexcept {
        return {t,t0,{{coefficients.Ap.data,coefficients.Ap.shape},{coefficients.Am.data,coefficients.Am.shape},{coefficients.A0.data,coefficients.A0.shape}}};
    }
};

struct WVFlux {
    WVComplexView Fp;
    WVComplexView Fm;
    WVComplexView F0;
};

struct WVTransformConstantStratificationConfiguration {
    std::uint32_t contractVersion = WVKernelContractVersion;
    std::size_t Nx = 0;
    std::size_t Ny = 0;
    std::size_t Nz = 0;
    std::size_t Nj = 0;
    double Lx = 0.0;
    double Ly = 0.0;
    double Lz = 0.0;
    double N0 = 0.0;
    double rho0 = 0.0;
    double g = 0.0;
    double planetaryRadius = 0.0;
    double rotationRate = 0.0;
    double latitude = 0.0;
    bool isHydrostatic = false;
    bool shouldAntialias = true;
};

// Exact identity for services and plans derived from a frozen transform
// configuration. Every listed value participates in the derived state.
bool sameTransformConfiguration(
    const WVTransformConstantStratificationConfiguration& first,
    const WVTransformConstantStratificationConfiguration& second) noexcept;

struct WVFourierMode {
    std::int64_t kMode = 0;
    std::int64_t lMode = 0;
    double k = 0.0;
    double l = 0.0;
    double Kh = 0.0;
    double cosAlpha = 0.0;
    double sinAlpha = 0.0;
    std::size_t dftPrimaryIndex = 0;
    std::size_t dftConjugateIndex = 0;
};

struct WVHalfSpectrumMappings {
    std::size_t NxHalf = 0;
    std::vector<std::size_t> directRows;
    std::vector<std::size_t> directWVIndices;
    std::vector<std::size_t> conjugatedRows;
    std::vector<std::size_t> conjugatedWVIndices;
    std::vector<std::size_t> storageRowsByWVIndex;
    std::vector<std::uint8_t> conjugatesStoredValueByWVIndex;
    std::vector<std::size_t> hermitianCompletionRows;
    std::vector<std::size_t> hermitianSourceRows;
    std::vector<std::size_t> selfConjugateRows;
    std::size_t persistentBytes() const noexcept;
};

struct WVConstantStratificationModes {
    double coriolisFrequency = 0.0;
    // Trapezoidal quadrature weight at z=-Lz. Bottom forcings divide by
    // this descriptor-owned weight before projecting their spatial tendency.
    double bottomQuadratureWeight = 0.0;
    std::vector<double> z;
    std::vector<double> j;
    std::vector<double> h0;
    // Vertical-only arrays use [Nj] ordering.
    std::vector<double> verticalWavenumber;
    std::vector<double> Fg;
    std::vector<double> Gg;
    std::vector<double> inertialScale;
    // G-family wave scaling is vertical-only for constant stratification.
    std::vector<double> gWaveScale;
    // Retains the #126 pre-scaling used to form ApmW at coefficient production.
    std::vector<double> apmWProjectionPrefactor;
    // Coefficient arrays use column-major [Nj,Nkl] ordering.
    std::vector<double> omega;
    std::vector<double> fWaveScale;
    std::vector<WVComplex64> UApField;
    std::vector<WVComplex64> VApField;
    std::vector<WVComplex64> WApField;
    std::vector<double> NApField;
    std::vector<WVComplex64> UA0Field;
    std::vector<WVComplex64> VA0Field;
    std::vector<double> NA0Field;
    std::vector<double> A0FromVorticity;
    std::vector<double> A0FromBuoyancy;
    std::vector<WVComplex64> ApmDProjection;
    std::vector<double> ApmNProjection;
    std::vector<double> ApmDScaled;
};

class WVTransformConstantStratificationDescriptor {
public:
    static WVKernelStatus create(const WVTransformConstantStratificationConfiguration& configuration, WVTransformConstantStratificationDescriptor& descriptor);

    const WVTransformConstantStratificationConfiguration& configuration() const noexcept { return configuration_; }
    const std::vector<WVFourierMode>& fourierModes() const noexcept { return fourierModes_; }
    const WVHalfSpectrumMappings& halfSpectrumMappings() const noexcept { return halfSpectrumMappings_; }
    const WVConstantStratificationModes& verticalModes() const noexcept { return verticalModes_; }
    std::size_t Nkl() const noexcept { return fourierModes_.size(); }
    WVShape2D spectralShape() const noexcept { return {configuration_.Nj, Nkl()}; }
    WVShape3D spatialShape() const noexcept { return {configuration_.Nx, configuration_.Ny, configuration_.Nz}; }
    double bottomQuadratureWeight() const noexcept { return verticalModes_.bottomQuadratureWeight; }
    std::size_t persistentBytes() const noexcept;

private:
    WVTransformConstantStratificationConfiguration configuration_;
    std::vector<WVFourierMode> fourierModes_;
    WVHalfSpectrumMappings halfSpectrumMappings_;
    WVConstantStratificationModes verticalModes_;
};

WVKernelStatus validateStateAndFlux(const WVTransformConstantStratificationDescriptor& descriptor, const WVState& state, const WVFlux& flux);

} // namespace wavevortex
