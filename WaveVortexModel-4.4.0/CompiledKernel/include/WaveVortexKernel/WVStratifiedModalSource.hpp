#pragma once
#include "WVSpectralOperators.hpp"

namespace wavevortex {

struct WVStratifiedModalGeometry {
    std::string transformClass, modelVersion;
    std::size_t Nx = 0, Ny = 0, Nz = 0, Nj = 0, Nkl = 0;
    double Lx = 0, Ly = 0, Lz = 0, g = 0, rho0 = 0;
    double latitude = 0, rotationRate = 0, planetaryRadius = 0;
    bool shouldAntialias = false;
    std::vector<double> x, y, z, j, k, l, N2, rho_nm0, dLnN2, P0, Q0, h_0, z_int;
    // Boussinesq-only data. Group membership is persisted, never inferred by
    // rounding radii. h_pm is [Nj,Nkl]; Ppm/Qpm are [Nj,K2unique].
    std::vector<double> K2unique, h_pm, Ppm, Qpm;
    std::vector<std::size_t> waveGroup;
    std::vector<WVRetainedModeKey> modes;
};

struct WVStratifiedModalInspection {
    WVStratifiedModalGeometry geometry;
    std::size_t scientificMatrixBytes = 0;
    // Matrix payloads are scanned through a fixed buffer, never retained here.
    static constexpr std::size_t matrixScanBytes = 4096*sizeof(double);
};

// Names refer to unpreconditioned scientific F/G operations. Matrix storage is
// column-major. Reconstruction maps modal to vertical grid; projection reverses it.
enum class WVStratifiedModalOperator { reconstructF, projectF, reconstructG, projectG, GToF, FToG,
    reconstructFw, projectFw, reconstructGw, projectGw, balancedGToWaveG,
    projectWaveDivergence, projectWaveVerticalVelocity };

enum class WVStratifiedScientificMatrix { PF0inv, QG0inv, PF0, QG0 };
struct WVScientificModalGroup {
    std::uint64_t identity = 0;
    std::vector<std::size_t> columns;
};

// Immutable scientific source boundary. Persistence adapters own validation and
// authoritative arrays; kernels consume data and prepare execution separately.
// Implementations must keep these identities and values fixed for their lifetime.
class WVStratifiedModalSource {
public:
    virtual ~WVStratifiedModalSource() = default;
    virtual const WVStratifiedModalGeometry& geometry() const noexcept = 0;
    virtual const std::string& sourceIdentity() const noexcept = 0;
    virtual const std::string& modeSetIdentity() const noexcept = 0;
    virtual std::size_t persistentBytes() const noexcept = 0;
    virtual WVKernelStatus prepareVertical(WVStratifiedModalOperator,
        WVComplexLayout input, WVComplexLayout output,
        std::unique_ptr<WVVerticalMatrixBackend>,
        std::unique_ptr<WVPreparedVerticalOperator>&,
        WVAccumulation = WVAccumulation::overwrite) const = 0;
    virtual WVKernelStatus horizontalSpecification(std::size_t planes,
        WVComplexRepresentation, const std::string& family,
        WVRetainedHorizontalSpecification&) const = 0;
};

} // namespace wavevortex
