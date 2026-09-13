#pragma once

#include "WVFFTEngine.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex {
namespace spectral_detail { struct VerticalGroupExecutorData; }

// Internal numerical service contract, independent of kernel contract 4 and
// portable extension source API v1. Scientific checkpoint records do not store it.
inline constexpr const char* WVSpectralOperatorContract = "wave-vortex-spectral-operators-v1";
enum class WVComplexRepresentation { interleaved, split };
enum class WVOperatorPlacement { outOfPlace, inPlace };
enum class WVFourierNormalization { forwardUnit, unitary, inverseUnit };
enum class WVMatrixAction { reconstruction, projection, crossFamily };
enum class WVAccumulation { overwrite, add };

// Strides are positive element strides, in complex values (interleaved) or
// doubles (each split component). Logical indices are [vertical, retained mode].
// modeSet names the exact ordered set, not an inferred radius or approximate key.
struct WVComplexLayout {
    std::size_t rows = 0, columns = 0, rowStride = 1, columnStride = 0;
    WVComplexRepresentation representation = WVComplexRepresentation::interleaved;
    std::string family, modeSet;
};
struct WVComplexInput {
    const WVComplex64* interleaved = nullptr;
    const double* real = nullptr;
    const double* imag = nullptr;
    // Capacity of the interleaved buffer, or minimum capacity of the two split buffers.
    std::size_t bytes = 0;
};
struct WVComplexOutput {
    WVComplex64* interleaved = nullptr;
    double* real = nullptr;
    double* imag = nullptr;
    std::size_t bytes = 0;
    WVComplexInput input() const noexcept { return {interleaved,real,imag,bytes}; }
};
struct WVRealGridLayout {
    std::size_t Nx = 0, Ny = 0, planes = 0;
    std::size_t xStride = 1, yStride = 0, planeStride = 0;
    std::string family;
};
struct WVRealInput { const double* data = nullptr; std::size_t bytes = 0; };
struct WVRealOutput { double* data = nullptr; std::size_t bytes = 0; };
struct WVRetainedAdvectionWork {
    // Prepared value spectra ordered u, v, w, eta; execution preserves them.
    std::array<WVComplexInput,4> base;
    // Ordered u, v, eta, w. Each active target begins as its prepared vertical
    // derivative spectrum and is replaced by the projected nonlinear flux.
    std::array<WVComplexOutput,4> targetSpectra;
    // Four contiguous [Nx,Ny,Nz] volumes ordered u, v, w, eta.
    WVRealOutput fields;
    // Optional completed physical tendencies. Three-target work is ordered
    // u, v, eta; four-target work is ordered u, v, w, eta.
    WVRealOutput tendencies;
    WVRealInput densityCorrection;
    std::size_t targets = 0;
};
struct WVRetainedAdvectionCounts {
    std::size_t columnInverses = 0;
    std::size_t rowInverses = 0;
    std::size_t reusedColumns = 0;
};
inline WVKernelStatus WVRetainedHorizontalPlan::inverseAndConsume(WVComplexInput,
    WVRealOutput,const WVRealOutputConsumer&) {
    return {WVKernelStatusCode::unsupportedOperation,"Provider has no inverse consumer."};
}
inline WVKernelStatus WVRetainedHorizontalPlan::prepareAdvection(std::size_t) {
    return {WVKernelStatusCode::unsupportedOperation,"Provider has no retained advection schedule."};
}
inline WVKernelStatus WVRetainedHorizontalPlan::advection(
    const WVRetainedAdvectionWork&,WVRetainedAdvectionCounts&) {
    return {WVKernelStatusCode::unsupportedOperation,"Provider has no retained advection schedule."};
}
struct WVRetainedModeKey { std::int64_t k = 0, l = 0; };
enum class WVRetainedHorizontalSchedule { fullFFT, streamingPrunedTile16 };
struct WVRetainedHorizontalSpecification {
    WVRealGridLayout grid;
    WVComplexLayout retained;
    double Lx = 0, Ly = 0;
    std::vector<WVRetainedModeKey> modes; // One representative per Hermitian orbit, in caller order.
    WVFourierNormalization normalization = WVFourierNormalization::forwardUnit;
    WVOperatorPlacement placement = WVOperatorPlacement::outOfPlace;
    WVRetainedHorizontalSchedule schedule = WVRetainedHorizontalSchedule::fullFFT;
    std::size_t outerWorkers = 1;
};

namespace spectral_detail {
struct HorizontalData;
struct HorizontalWorkspaceData;
struct VerticalData;
struct VerticalWorkspaceData;
}
class WVRetainedHorizontalWorkspace {
public:
    ~WVRetainedHorizontalWorkspace();
    WVRetainedHorizontalWorkspace(const WVRetainedHorizontalWorkspace&) = delete;
    WVRetainedHorizontalWorkspace& operator=(const WVRetainedHorizontalWorkspace&) = delete;
    std::size_t persistentBytes() const noexcept;
    std::size_t planBytesLowerBound() const noexcept;
    const char* scheduleIdentifier() const noexcept;
    const void* sharedResourceIdentity() const noexcept;
    std::size_t sharedResourceBytes() const noexcept;
    std::size_t workerCount() const noexcept;
private:
    friend class WVRetainedHorizontalOperator;
    WVRetainedHorizontalWorkspace();
    std::unique_ptr<spectral_detail::HorizontalWorkspaceData> data_;
};
class WVRetainedHorizontalOperator {
public:
    // Full FFT is the default; a pruned schedule requires explicit opt-in.
    static WVKernelStatus create(const WVRetainedHorizontalSpecification&, std::unique_ptr<WVFFTEngine>, std::unique_ptr<WVRetainedHorizontalOperator>&);
    static WVKernelStatus createShared(const WVRetainedHorizontalSpecification&, std::shared_ptr<WVFFTEngine>, std::unique_ptr<WVRetainedHorizontalOperator>&);
    // false omits full-grid derivative preparation and bounds fallback scratch
    // to one horizontal plane. spatialDerivative then returns unsupported.
    // With a retained provider, prepared full-grid derivatives also stream
    // through one plane; their frequency coverage is still the entire grid.
    WVKernelStatus createWorkspace(std::unique_ptr<WVRetainedHorizontalWorkspace>&, bool prepareSpatialDerivative = true) const;
    WVKernelStatus forward(WVRetainedHorizontalWorkspace&, WVRealInput, WVComplexOutput) const;
    WVKernelStatus inverse(WVRetainedHorizontalWorkspace&, WVComplexInput, WVRealOutput) const;
    // Consumer ranges refer to a contiguous [Nx,Ny,planes] output. The complete
    // output is retained, including when consumers run inside provider workers.
    WVKernelStatus inverseAndConsume(WVRetainedHorizontalWorkspace&, WVComplexInput,
        WVRealOutput,const WVRealOutputConsumer&) const;
    WVKernelStatus prepareAdvection(WVRetainedHorizontalWorkspace&,std::size_t targets) const;
    bool supportsAdvection(const WVRetainedHorizontalWorkspace&,std::size_t targets) const noexcept;
    WVKernelStatus advection(WVRetainedHorizontalWorkspace&,const WVRetainedAdvectionWork&,
        WVRetainedAdvectionCounts&) const;
    std::size_t persistentBytes() const noexcept;
    std::size_t providerBytesLowerBound() const noexcept;
    // Full-grid derivative before retained-mode projection (e.g. passive tracers).
    // Uses all resolved Fourier modes, with the selected even-grid Nyquist axis
    // zeroed for odd orders and retained for even orders, independent of the
    // retained spectral subset.
    WVKernelStatus spatialDerivative(WVRetainedHorizontalWorkspace&, WVRealInput, WVRealOutput, bool xDerivative) const;
    WVKernelStatus spatialDerivative(WVRetainedHorizontalWorkspace&, WVRealInput, WVRealOutput,
        bool xDerivative, unsigned order) const;
private:
    WVRetainedHorizontalOperator() = default;
    std::shared_ptr<const spectral_detail::HorizontalData> data_;
};

// Exact authoritative source identity. source must change with geometry,
// stratification, retained basis or normalization; name identifies one operator
// in that source. Equal numerical matrices with different identities stay distinct.
struct WVMatrixIdentity { std::string source, name; };
struct WVRealMatrixInput {
    const double* data = nullptr;
    std::size_t rows = 0, columns = 0, rowStride = 1, columnStride = 0, bytes = 0;
};
struct WVVerticalMatrixRecord {
    WVMatrixIdentity identity;
    WVRealMatrixInput values;
    WVMatrixAction action = WVMatrixAction::reconstruction;
    std::string inputFamily, outputFamily;
};
struct WVVerticalGroup {
    std::uint64_t identity = 0; // Exact source group identity; never derived from integer radius.
    std::size_t matrix = 0;
    std::vector<std::size_t> modes; // May be discontiguous; each retained column occurs exactly once.
};
struct WVVerticalSpecification {
    std::string sourceIdentity;
    WVMatrixAction action = WVMatrixAction::reconstruction;
    std::size_t Nz = 0, Nj = 0;
    WVComplexLayout input, output;
    std::vector<WVVerticalMatrixRecord> matrices;
    std::vector<WVVerticalGroup> groups;
    WVAccumulation accumulation = WVAccumulation::overwrite;
    WVOperatorPlacement placement = WVOperatorPlacement::outOfPlace;
};

// Adapter boundary: column-major A[m,k], B[k,n], C[m,n]. The prepared service
// validates dimensions/placement and supplies immutable matrices. Implementations
// must perform C=A*B+beta*C without application allocation or launching workers.
class WVVerticalMatrixBackend {
public:
    virtual ~WVVerticalMatrixBackend() = default;
    virtual const char* identifier() const noexcept = 0;
    virtual std::size_t maximumDimension() const noexcept = 0;
    virtual std::size_t persistentBytes() const noexcept = 0;
    // True only when independent calls on this immutable backend may execute
    // concurrently. Scheduling policy remains outside the backend.
    virtual bool supportsConcurrentCalls() const noexcept { return false; }
    virtual void split(std::size_t m, std::size_t k, std::size_t n, const double* a,
        const double* br, const double* bi, std::size_t ldb, double* cr, double* ci, std::size_t ldc, double beta) const noexcept = 0;
    virtual void interleaved(std::size_t m, std::size_t k, std::size_t n, const WVComplex64* a,
        const WVComplex64* b, std::size_t ldb, WVComplex64* c, std::size_t ldc, double beta) const noexcept = 0;
};
WVKernelStatus WVCreateScalarMatrixBackend(std::unique_ptr<WVVerticalMatrixBackend>&);

// One persistent synchronous executor may be shared by every prepared vertical
// operator owned by a kernel. It statically partitions the prepared group list.
class WVVerticalGroupExecutor {
public:
    static WVKernelStatus create(std::size_t workers,std::unique_ptr<WVVerticalGroupExecutor>&);
    ~WVVerticalGroupExecutor();
    WVVerticalGroupExecutor(const WVVerticalGroupExecutor&) = delete;
    WVVerticalGroupExecutor& operator=(const WVVerticalGroupExecutor&) = delete;
    std::size_t workerCount() const noexcept;
    std::size_t persistentBytes() const noexcept;
private:
    friend class WVPreparedVerticalOperator;
    WVVerticalGroupExecutor();
    std::unique_ptr<spectral_detail::VerticalGroupExecutorData> data_;
};

class WVVerticalWorkspace {
public:
    ~WVVerticalWorkspace();
    WVVerticalWorkspace(const WVVerticalWorkspace&) = delete;
    WVVerticalWorkspace& operator=(const WVVerticalWorkspace&) = delete;
    std::size_t persistentBytes() const noexcept;
private:
    friend class WVPreparedVerticalOperator;
    WVVerticalWorkspace();
    std::unique_ptr<spectral_detail::VerticalWorkspaceData> data_;
};
class WVPreparedVerticalOperator {
public:
    static WVKernelStatus create(const WVVerticalSpecification&, std::unique_ptr<WVVerticalMatrixBackend>, std::unique_ptr<WVPreparedVerticalOperator>&);
    WVKernelStatus createWorkspace(std::unique_ptr<WVVerticalWorkspace>&) const;
    WVKernelStatus createWorkspace(std::size_t groupWorkers,std::unique_ptr<WVVerticalWorkspace>&) const;
    WVKernelStatus execute(WVVerticalWorkspace&, WVComplexInput, WVComplexOutput) const;
    WVKernelStatus execute(WVVerticalWorkspace&, WVVerticalGroupExecutor&,
        WVComplexInput,WVComplexOutput) const;
    // Execute the exact prepared scientific matrix assigned to one retained
    // column. Other output columns are untouched.
    WVKernelStatus executeColumn(WVVerticalWorkspace&, WVComplexInput, WVComplexOutput,
        std::size_t retainedColumn) const;
    std::size_t persistentBytes() const noexcept;
    std::size_t matrixBytes() const noexcept;
    std::size_t uniqueMatrixCount() const noexcept;
    std::size_t preparedGroupCount() const noexcept;
    bool supportsConcurrentCalls() const noexcept;
    const char* backendIdentifier() const noexcept;
private:
    WVPreparedVerticalOperator() = default;
    std::shared_ptr<const spectral_detail::VerticalData> data_;
};

} // namespace wavevortex
