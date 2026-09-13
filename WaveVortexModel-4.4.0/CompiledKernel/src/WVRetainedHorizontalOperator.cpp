#include "WVSpectralValidation.hpp"
#include "WVHorizontalDerivativeMultiplier.hpp"
#include <new>

namespace wavevortex {
namespace spectral_detail {
struct HorizontalMode {
    std::size_t row, partner;
    bool conjugate, self;
};
struct HorizontalData {
    WVRetainedHorizontalSpecification spec;
    std::vector<HorizontalMode> mapping;
    std::shared_ptr<WVFFTEngine> engine;
    std::size_t realSpan = 0, complexSpan = 0, planeSize = 0, halfSize = 0;
    double forwardScale = 1, inverseScale = 1;
};
struct HorizontalWorkspaceData {
    std::shared_ptr<const HorizontalData> owner;
    std::unique_ptr<WVFFTPlan> forward, inverse;
    std::unique_ptr<WVRetainedHorizontalPlan> retained;
    bool derivativePrepared = true, batchedFullFFT = true;
    std::vector<double> real;
    std::vector<WVComplex64> half;
    std::atomic<bool> active{false};
};

std::size_t gridSpan(const WVRealGridLayout& grid) {
    if (grid.family.empty() || !grid.Nx || !grid.Ny || !grid.planes || !grid.xStride || !grid.yStride || !grid.planeStride)
        throw std::invalid_argument("Grid dimensions, strides and family must be specified.");
    // A permuted, padded nested grid is supported. Interleaving dimensions with
    // overlapping bounding intervals is rejected rather than guessed safe.
    std::array<std::pair<std::size_t,std::size_t>,3> dims{{{grid.xStride,grid.Nx},{grid.yStride,grid.Ny},{grid.planeStride,grid.planes}}};
    std::sort(dims.begin(),dims.end());
    std::size_t extent = 1;
    for (const auto& dim : dims) {
        if (dim.second > 1 && dim.first < extent) throw std::invalid_argument("Grid strides overlap.");
        const auto extra = product(dim.first,dim.second-1);
        if (extra > static_cast<std::size_t>(PTRDIFF_MAX)-extent) throw std::overflow_error("Grid span overflow.");
        extent += extra;
    }
    return product(extent,sizeof(double));
}
std::size_t index(std::int64_t mode, std::size_t count) {
    const auto half = static_cast<std::int64_t>(count/2);
    if (mode < -half || mode > half) throw std::invalid_argument("Retained Fourier mode lies outside the grid.");
    return mode < 0 ? static_cast<std::size_t>(static_cast<std::int64_t>(count)+mode) : static_cast<std::size_t>(mode);
}
WVFFTPlanSpecification fftSpecification(const HorizontalData& data, bool inverse, bool full = true) {
    const auto& g = data.spec.grid;
    const auto planes = full ? g.planes : 1;
    WVFFTPlanSpecification s;
    s.kind = inverse ? WVFFTPlanKind::horizontalComplexToReal2D : WVFFTPlanKind::horizontalRealToComplex2D;
    s.transformDimensions = {{g.Ny,static_cast<std::ptrdiff_t>(g.Nx),static_cast<std::ptrdiff_t>(g.Nx/2+1)},{g.Nx,1,1}};
    s.batchDimensions = {{planes,static_cast<std::ptrdiff_t>(data.planeSize),static_cast<std::ptrdiff_t>(data.halfSize)}};
    s.inputBytes = product(product(data.planeSize,planes),sizeof(double));
    s.outputBytes = product(product(data.halfSize,planes),sizeof(WVComplex64));
    s.destroysInput = inverse; // Only disposable workspace is passed to the full inverse.
    if (inverse) {
        for (auto& d : s.transformDimensions) std::swap(d.inputStride,d.outputStride);
        for (auto& d : s.batchDimensions) std::swap(d.inputStride,d.outputStride);
        std::swap(s.inputBytes,s.outputBytes);
    }
    return s;
}
WVKernelStatus validateHorizontalBuffers(const HorizontalData& d, WVRealInput real, WVComplexInput complex) {
    if (real.bytes < d.realSpan) return {WVKernelStatusCode::invalidShape,"Real grid buffer capacity is too small."};
    if (!addressFits(real.data,d.realSpan,alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid real grid storage."};
    auto status = validateStorage(d.spec.retained.representation,d.complexSpan,complex);
    if (!status) return status;
    if (realOverlap(real.data,d.realSpan,complex,d.complexSpan)) return {WVKernelStatusCode::overlappingArrays,"Horizontal input and output overlap."};
    return WVKernelStatus::ok();
}
} // namespace spectral_detail
using namespace spectral_detail;

WVRetainedHorizontalWorkspace::WVRetainedHorizontalWorkspace() = default;
WVRetainedHorizontalWorkspace::~WVRetainedHorizontalWorkspace() = default;
std::size_t WVRetainedHorizontalWorkspace::persistentBytes() const noexcept {
    return sizeof(*this)+sizeof(*data_)+data_->real.capacity()*sizeof(double)+data_->half.capacity()*sizeof(WVComplex64)+(data_->retained ? data_->retained->persistentBytes() : 0);
}
std::size_t WVRetainedHorizontalWorkspace::planBytesLowerBound() const noexcept {
    return (data_->forward ? data_->forward->persistentBytes()+data_->inverse->persistentBytes() : 0)+
        (data_->retained ? data_->retained->planBytesLowerBound() : 0);
}
const char* WVRetainedHorizontalWorkspace::scheduleIdentifier() const noexcept {
    return data_->retained ? data_->retained->identifier() : data_->batchedFullFFT ? "full-fft-gather" : "plane-streamed-full-fft-gather";
}
const void* WVRetainedHorizontalWorkspace::sharedResourceIdentity() const noexcept {
    return data_->retained ? data_->retained->sharedResourceIdentity() : nullptr;
}
std::size_t WVRetainedHorizontalWorkspace::sharedResourceBytes() const noexcept {
    return data_->retained ? data_->retained->sharedResourceBytes() : 0;
}
std::size_t WVRetainedHorizontalWorkspace::workerCount() const noexcept {
    return data_->retained ? data_->retained->workerCount() : 1;
}
std::size_t WVRetainedHorizontalOperator::persistentBytes() const noexcept {
    return sizeof(*this)+sizeof(*data_)+identityBytes(data_->spec.retained)+data_->spec.grid.family.capacity()+
        data_->spec.modes.capacity()*sizeof(WVRetainedModeKey)+data_->mapping.capacity()*sizeof(HorizontalMode);
}
std::size_t WVRetainedHorizontalOperator::providerBytesLowerBound() const noexcept { return data_->engine->persistentBytes(); }
WVKernelStatus WVRetainedHorizontalOperator::create(const WVRetainedHorizontalSpecification& spec,
    std::unique_ptr<WVFFTEngine> engine, std::unique_ptr<WVRetainedHorizontalOperator>& result) {
    try { return createShared(spec,std::shared_ptr<WVFFTEngine>(std::move(engine)),result); }
    catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Horizontal provider ownership allocation failed."}; }
}
WVKernelStatus WVRetainedHorizontalOperator::createShared(const WVRetainedHorizontalSpecification& spec,
    std::shared_ptr<WVFFTEngine> engine, std::unique_ptr<WVRetainedHorizontalOperator>& result) {
    try {
        if (!engine) return {WVKernelStatusCode::invalidPointer,"A full FFT provider is required."};
        if (spec.placement != WVOperatorPlacement::outOfPlace) return {WVKernelStatusCode::unsupportedOperation,"Retained horizontal execution is out-of-place."};
        if (!std::isfinite(spec.Lx) || !std::isfinite(spec.Ly) || spec.Lx <= 0 || spec.Ly <= 0)
            return {WVKernelStatusCode::invalidConfiguration,"Domain lengths must be finite and positive."};
        if (spec.grid.Nx > static_cast<std::size_t>(INT64_MAX) || spec.grid.Ny > static_cast<std::size_t>(INT64_MAX))
            return {WVKernelStatusCode::sizeOverflow,"Fourier grid exceeds signed mode-key range."};
        if (!spec.outerWorkers || (spec.schedule != WVRetainedHorizontalSchedule::fullFFT &&
            spec.schedule != WVRetainedHorizontalSchedule::streamingPrunedTile16))
            return {WVKernelStatusCode::invalidConfiguration,"Invalid horizontal schedule or worker count."};
        const auto rs = gridSpan(spec.grid), cs = validate(spec.retained);
        if (spec.retained.rows != spec.grid.planes || spec.retained.columns != spec.modes.size() || spec.grid.family != spec.retained.family)
            return {WVKernelStatusCode::invalidShape,"Horizontal families and extents must match the retained set."};
        auto d = std::make_shared<HorizontalData>();
        d->spec = spec; d->engine = std::move(engine); d->realSpan = rs; d->complexSpan = cs;
        d->planeSize = product(spec.grid.Nx,spec.grid.Ny);
        d->halfSize = product(spec.grid.Nx/2+1,spec.grid.Ny);
        // Check all execution allocations at preparation, before allocating the visited set.
        product(product(d->planeSize,spec.grid.planes),sizeof(double));
        product(product(d->halfSize,spec.grid.planes),sizeof(WVComplex64));
        switch (spec.normalization) {
            case WVFourierNormalization::forwardUnit: d->forwardScale = 1.0/d->planeSize; break;
            case WVFourierNormalization::unitary: d->forwardScale = d->inverseScale = 1.0/std::sqrt(static_cast<double>(d->planeSize)); break;
            case WVFourierNormalization::inverseUnit: d->inverseScale = 1.0/d->planeSize; break;
            default: return {WVKernelStatusCode::invalidConfiguration,"Unknown Fourier normalization."};
        }
        std::vector<bool> visited(d->planeSize,false);
        const auto nx = spec.grid.Nx, ny = spec.grid.Ny, half = nx/2+1;
        for (const auto& key : spec.modes) {
            const auto x = index(key.k,nx), y = index(key.l,ny);
            const auto cx = (nx-x)%nx, cy = (ny-y)%ny;
            if (visited[x+nx*y] || visited[cx+nx*cy]) return {WVKernelStatusCode::invalidConfiguration,"Retained set repeats a Hermitian orbit."};
            visited[x+nx*y] = visited[cx+nx*cy] = true;
            const bool conjugate = x > nx/2;
            const auto sx = conjugate ? cx : x, sy = conjugate ? cy : y;
            const bool boundary = sx == 0 || (nx%2 == 0 && sx == nx/2);
            const auto row = sx+half*sy;
            d->mapping.push_back({row,boundary ? sx+half*((ny-sy)%ny) : row,conjugate,x == cx && y == cy});
        }
        auto op = std::unique_ptr<WVRetainedHorizontalOperator>(new WVRetainedHorizontalOperator);
        op->data_ = std::move(d); result = std::move(op);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Horizontal setup allocation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
      catch (const std::invalid_argument& e) { return {WVKernelStatusCode::invalidConfiguration,e.what()}; }
}
WVKernelStatus WVRetainedHorizontalOperator::createWorkspace(std::unique_ptr<WVRetainedHorizontalWorkspace>& result, bool prepareSpatialDerivative) const {
    try {
        auto w = std::unique_ptr<WVRetainedHorizontalWorkspace>(new WVRetainedHorizontalWorkspace);
        w->data_ = std::make_unique<HorizontalWorkspaceData>();
        auto& d = *w->data_; d.owner = data_; d.derivativePrepared = prepareSpatialDerivative;
        if (data_->spec.schedule == WVRetainedHorizontalSchedule::streamingPrunedTile16) {
            auto status = data_->engine->createRetainedHorizontalPlan(data_->spec,d.retained);
            if (!status && status.code != WVKernelStatusCode::unsupportedOperation) return status;
            if (status && !d.retained) return {WVKernelStatusCode::fftPlanFailure,"Provider returned an empty retained plan."};
        }
        d.batchedFullFFT = prepareSpatialDerivative && !d.retained;
        if (prepareSpatialDerivative || !d.retained) {
            const auto planes = d.batchedFullFFT ? data_->spec.grid.planes : 1;
            d.real.resize(product(data_->planeSize,planes));
            d.half.resize(product(data_->halfSize,planes));
            auto status = data_->engine->createPlan(fftSpecification(*data_,false,d.batchedFullFFT),d.forward);
            if (!status) return status;
            status = data_->engine->createPlan(fftSpecification(*data_,true,d.batchedFullFFT),d.inverse);
            if (!status) return status;
            if (!d.forward || !d.inverse) return {WVKernelStatusCode::fftPlanFailure,"Provider returned an empty plan."};
        }
        result = std::move(w); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Horizontal workspace allocation failed."}; }
}
WVKernelStatus WVRetainedHorizontalOperator::forward(WVRetainedHorizontalWorkspace& workspace, WVRealInput input, WVComplexOutput output) const {
    auto& w = *workspace.data_; const auto& d = *data_;
    if (w.owner != data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another horizontal operator."};
    auto status = validateHorizontalBuffers(d,input,output.input()); if (!status) return status;
    ActiveCall guard(w.active); if (!guard.entered) return {WVKernelStatusCode::reentrantExecution,"Horizontal workspace is active."};
    if (w.retained) return w.retained->forward(input,output);
    const auto& g = d.spec.grid; const auto& l = d.spec.retained;
    if (!w.batchedFullFFT) {
        for (std::size_t p = 0; p < g.planes; ++p) {
            for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
                w.real[y*g.Nx+x] = input.data[p*g.planeStride+y*g.yStride+x*g.xStride];
            status = w.forward->execute(w.real.data(),w.half.data()); if (!status) return status;
            for (std::size_t mode = 0; mode < d.mapping.size(); ++mode) {
                const auto& map = d.mapping[mode]; auto value = w.half[map.row];
                value.real *= d.forwardScale; value.imag *= map.conjugate ? -d.forwardScale : d.forwardScale;
                if (map.self) value.imag = 0;
                write(output,p*l.rowStride+mode*l.columnStride,value);
            }
        }
        return WVKernelStatus::ok();
    }
    for (std::size_t p = 0; p < g.planes; ++p) for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
        w.real[p*d.planeSize+y*g.Nx+x] = input.data[p*g.planeStride+y*g.yStride+x*g.xStride];
    status = w.forward->execute(w.real.data(),w.half.data()); if (!status) return status;
    for (std::size_t mode = 0; mode < d.mapping.size(); ++mode) for (std::size_t p = 0; p < g.planes; ++p) {
        const auto& map = d.mapping[mode]; auto value = w.half[p*d.halfSize+map.row];
        value.real *= d.forwardScale; value.imag *= map.conjugate ? -d.forwardScale : d.forwardScale;
        if (map.self) value.imag = 0;
        write(output,p*l.rowStride+mode*l.columnStride,value);
    }
    return WVKernelStatus::ok();
}
WVKernelStatus WVRetainedHorizontalOperator::inverse(WVRetainedHorizontalWorkspace& workspace, WVComplexInput input, WVRealOutput output) const {
    return inverseAndConsume(workspace,input,output,{});
}
WVKernelStatus WVRetainedHorizontalOperator::inverseAndConsume(WVRetainedHorizontalWorkspace& workspace,
    WVComplexInput input,WVRealOutput output,const WVRealOutputConsumer& consumer) const {
    auto& w = *workspace.data_; const auto& d = *data_;
    if (w.owner != data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another horizontal operator."};
    auto status = validateHorizontalBuffers(d,{output.data,output.bytes},input); if (!status) return status;
    ActiveCall guard(w.active); if (!guard.entered) return {WVKernelStatusCode::reentrantExecution,"Horizontal workspace is active."};
    const auto& g = d.spec.grid; const auto& l = d.spec.retained;
    if (consumer.consume && (g.xStride!=1 || g.yStride!=g.Nx || g.planeStride!=d.planeSize))
        return {WVKernelStatusCode::invalidConfiguration,"Inverse consumers require contiguous physical output."};
    const auto consumePlane=[&](std::size_t p) {
        if (consumer.consume) consumer.consume(consumer.context,p*d.planeSize,(p+1)*d.planeSize,output.data);
    };
    // Validate all self-conjugate values before writing any caller output.
    for (std::size_t mode = 0; mode < d.mapping.size(); ++mode) if (d.mapping[mode].self)
        for (std::size_t p = 0; p < g.planes; ++p) if (read(input,p*l.rowStride+mode*l.columnStride).imag != 0)
            return {WVKernelStatusCode::invalidConfiguration,"Self-conjugate Fourier values must be real."};
    if (w.retained) {
        if (consumer.consume && w.retained->supportsInverseConsumer())
            return w.retained->inverseAndConsume(input,output,consumer);
        status=w.retained->inverse(input,output);
        if (status) for (std::size_t p=0;p<g.planes;++p) consumePlane(p);
        return status;
    }
    if (!w.batchedFullFFT) {
        for (std::size_t p = 0; p < g.planes; ++p) {
            std::fill(w.half.begin(),w.half.end(),WVComplex64{});
            for (std::size_t mode = 0; mode < d.mapping.size(); ++mode) {
                const auto& map = d.mapping[mode]; auto value = read(input,p*l.rowStride+mode*l.columnStride);
                if (map.conjugate) value.imag = -value.imag;
                w.half[map.row] = value;
                if (map.partner != map.row) w.half[map.partner] = {value.real,-value.imag};
            }
            status = w.inverse->execute(w.half.data(),w.real.data()); if (!status) return status;
            for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
                output.data[p*g.planeStride+y*g.yStride+x*g.xStride] = d.inverseScale*w.real[y*g.Nx+x];
            consumePlane(p);
        }
        return WVKernelStatus::ok();
    }
    std::fill(w.half.begin(),w.half.end(),WVComplex64{});
    for (std::size_t mode = 0; mode < d.mapping.size(); ++mode) for (std::size_t p = 0; p < g.planes; ++p) {
        const auto& map = d.mapping[mode]; auto value = read(input,p*l.rowStride+mode*l.columnStride);
        if (map.conjugate) value.imag = -value.imag;
        w.half[p*d.halfSize+map.row] = value;
        if (map.partner != map.row) w.half[p*d.halfSize+map.partner] = {value.real,-value.imag};
    }
    status = w.inverse->execute(w.half.data(),w.real.data()); if (!status) return status;
    for (std::size_t p = 0; p < g.planes; ++p) for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
        output.data[p*g.planeStride+y*g.yStride+x*g.xStride] = d.inverseScale*w.real[p*d.planeSize+y*g.Nx+x];
    for (std::size_t p=0;p<g.planes;++p) consumePlane(p);
    return WVKernelStatus::ok();
}
WVKernelStatus WVRetainedHorizontalOperator::prepareAdvection(
    WVRetainedHorizontalWorkspace& workspace,std::size_t targets) const {
    auto& w=*workspace.data_;
    if (w.owner!=data_)
        return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another horizontal operator."};
    if (!w.retained || !w.retained->supportsAdvection(targets))
        return {WVKernelStatusCode::unsupportedOperation,"Workspace has no compatible retained advection schedule."};
    ActiveCall guard(w.active);
    if (!guard.entered)
        return {WVKernelStatusCode::reentrantExecution,"Horizontal workspace is active."};
    return w.retained->prepareAdvection(targets);
}
bool WVRetainedHorizontalOperator::supportsAdvection(
    const WVRetainedHorizontalWorkspace& workspace,std::size_t targets) const noexcept {
    const auto& w=*workspace.data_;
    return w.owner==data_ && w.retained && w.retained->supportsAdvection(targets);
}
WVKernelStatus WVRetainedHorizontalOperator::advection(
    WVRetainedHorizontalWorkspace& workspace,const WVRetainedAdvectionWork& work,
    WVRetainedAdvectionCounts& counts) const {
    auto& w=*workspace.data_; const auto& d=*data_;
    if (w.owner!=data_)
        return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another horizontal operator."};
    if (work.targets!=3 && work.targets!=4)
        return {WVKernelStatusCode::invalidConfiguration,"Retained advection requires three or four targets."};
    if (!w.retained || !w.retained->supportsAdvection(work.targets))
        return {WVKernelStatusCode::unsupportedOperation,"Workspace has no compatible retained advection schedule."};
    const auto volumeBytes=product(product(d.planeSize,d.spec.grid.planes),sizeof(double));
    if (volumeBytes>static_cast<std::size_t>(PTRDIFF_MAX)/4)
        return {WVKernelStatusCode::sizeOverflow,"Advection physical-field span overflows ptrdiff_t."};
    const auto fieldBytes=4*volumeBytes;
    if (work.fields.bytes<fieldBytes)
        return {WVKernelStatusCode::invalidShape,"Advection physical-field capacity is too small."};
    if (!addressFits(work.fields.data,fieldBytes,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid advection physical-field storage."};
    const auto correctionBytes=product(d.spec.grid.planes,sizeof(double));
    if (work.densityCorrection.bytes<correctionBytes)
        return {WVKernelStatusCode::invalidShape,"Advection density-correction capacity is too small."};
    if (!addressFits(work.densityCorrection.data,correctionBytes,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid advection density-correction storage."};
    if (overlap(work.fields.data,fieldBytes,work.densityCorrection.data,correctionBytes))
        return {WVKernelStatusCode::overlappingArrays,"Advection fields overlap density correction."};
    const auto tendencyBytes=work.targets*volumeBytes;
    if ((work.tendencies.data==nullptr)!=(work.tendencies.bytes==0))
        return {WVKernelStatusCode::invalidPointer,"Incomplete advection tendency storage."};
    if (work.tendencies.data) {
        if (work.tendencies.bytes<tendencyBytes)
            return {WVKernelStatusCode::invalidShape,"Advection tendency capacity is too small."};
        if (!addressFits(work.tendencies.data,tendencyBytes,alignof(double)))
            return {WVKernelStatusCode::invalidPointer,"Invalid advection tendency storage."};
        if (overlap(work.tendencies.data,tendencyBytes,work.fields.data,fieldBytes) ||
            overlap(work.tendencies.data,tendencyBytes,work.densityCorrection.data,correctionBytes))
            return {WVKernelStatusCode::overlappingArrays,"Advection tendencies overlap real input storage."};
    }
    for (std::size_t field=0;field<4;++field) {
        auto status=validateStorage(d.spec.retained.representation,d.complexSpan,work.base[field]);
        if (!status) return status;
        if (realOverlap(work.fields.data,fieldBytes,work.base[field],d.complexSpan) ||
            realOverlap(work.densityCorrection.data,correctionBytes,work.base[field],d.complexSpan) ||
            (work.tendencies.data && realOverlap(work.tendencies.data,tendencyBytes,
                work.base[field],d.complexSpan)))
            return {WVKernelStatusCode::overlappingArrays,"Advection base spectra overlap real storage."};
    }
    for (std::size_t target=0;target<work.targets;++target) {
        const auto output=work.targetSpectra[target].input();
        auto status=validateStorage(d.spec.retained.representation,d.complexSpan,output);
        if (!status) return status;
        if (realOverlap(work.fields.data,fieldBytes,output,d.complexSpan) ||
            realOverlap(work.densityCorrection.data,correctionBytes,output,d.complexSpan) ||
            (work.tendencies.data && realOverlap(work.tendencies.data,tendencyBytes,
                output,d.complexSpan)))
            return {WVKernelStatusCode::overlappingArrays,"Advection target spectra overlap real storage."};
        for (std::size_t field=0;field<4;++field)
            if (storageOverlap(output,d.complexSpan,work.base[field],d.complexSpan))
                return {WVKernelStatusCode::overlappingArrays,"Advection target spectra overlap base spectra."};
        for (std::size_t other=0;other<target;++other)
            if (storageOverlap(output,d.complexSpan,work.targetSpectra[other].input(),d.complexSpan))
                return {WVKernelStatusCode::overlappingArrays,"Advection target spectra overlap each other."};
    }
    // Every input is checked before the first target or physical field can be written.
    const auto& layout=d.spec.retained;
    const auto hermitian=[&](WVComplexInput input) {
        for (std::size_t mode=0;mode<d.mapping.size();++mode) if (d.mapping[mode].self)
            for (std::size_t plane=0;plane<d.spec.grid.planes;++plane)
                if (read(input,plane*layout.rowStride+mode*layout.columnStride).imag!=0) return false;
        return true;
    };
    for (const auto& input:work.base) if (!hermitian(input))
        return {WVKernelStatusCode::invalidConfiguration,"Advection base self-conjugate values must be real."};
    for (std::size_t target=0;target<work.targets;++target)
        if (!hermitian(work.targetSpectra[target].input()))
            return {WVKernelStatusCode::invalidConfiguration,"Advection target self-conjugate values must be real."};
    for (std::size_t plane=0;plane<d.spec.grid.planes;++plane)
        if (!std::isfinite(work.densityCorrection.data[plane]))
            return {WVKernelStatusCode::invalidConfiguration,"Advection density correction must be finite."};
    ActiveCall guard(w.active);
    if (!guard.entered)
        return {WVKernelStatusCode::reentrantExecution,"Horizontal workspace is active."};
    return w.retained->advection(work,counts);
}
WVKernelStatus WVRetainedHorizontalOperator::spatialDerivative(WVRetainedHorizontalWorkspace& workspace, WVRealInput input, WVRealOutput output, bool xDerivative) const {
    return spatialDerivative(workspace,input,output,xDerivative,1);
}
WVKernelStatus WVRetainedHorizontalOperator::spatialDerivative(WVRetainedHorizontalWorkspace& workspace,
    WVRealInput input,WVRealOutput output,bool xDerivative,unsigned order) const {
    auto& w = *workspace.data_; const auto& d = *data_;
    if (!order) return {WVKernelStatusCode::invalidConfiguration,"Horizontal derivative order must be positive."};
    if (!w.derivativePrepared) return {WVKernelStatusCode::unsupportedOperation,"Workspace omitted full-grid derivative preparation."};
    if (w.owner != data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another horizontal operator."};
    if (input.bytes < d.realSpan || output.bytes < d.realSpan) return {WVKernelStatusCode::invalidShape,"Derivative grid capacity is too small."};
    if (!addressFits(input.data,d.realSpan,alignof(double)) || !addressFits(output.data,d.realSpan,alignof(double)))
        return {WVKernelStatusCode::invalidPointer,"Invalid derivative grid storage."};
    if (overlap(input.data,d.realSpan,output.data,d.realSpan)) return {WVKernelStatusCode::overlappingArrays,"Derivative input and output overlap."};
    ActiveCall guard(w.active); if (!guard.entered) return {WVKernelStatusCode::reentrantExecution,"Horizontal workspace is active."};
    const auto& g = d.spec.grid;
    const auto differentiate = [&](WVComplex64* values) {
        return kernel_detail::applyHorizontalDerivativeMultiplier(values,g.Nx,g.Ny,1,
            0,1,d.spec.Lx,d.spec.Ly,xDerivative,order);
    };
    if (!w.batchedFullFFT) {
        for (std::size_t p = 0; p < g.planes; ++p) {
            for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
                w.real[y*g.Nx+x] = input.data[p*g.planeStride+y*g.yStride+x*g.xStride];
            auto status = w.forward->execute(w.real.data(),w.half.data()); if (!status) return status;
            status=differentiate(w.half.data()); if (!status) return status;
            status = w.inverse->execute(w.half.data(),w.real.data()); if (!status) return status;
            for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
                output.data[p*g.planeStride+y*g.yStride+x*g.xStride] = w.real[y*g.Nx+x];
        }
        return WVKernelStatus::ok();
    }
    for (std::size_t p = 0; p < g.planes; ++p) for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
        w.real[p*d.planeSize+y*g.Nx+x] = input.data[p*g.planeStride+y*g.yStride+x*g.xStride];
    auto status = w.forward->execute(w.real.data(),w.half.data()); if (!status) return status;
    for (std::size_t p = 0; p < g.planes; ++p) {
        status=differentiate(w.half.data()+p*d.halfSize); if (!status) return status;
    }
    status = w.inverse->execute(w.half.data(),w.real.data()); if (!status) return status;
    for (std::size_t p = 0; p < g.planes; ++p) for (std::size_t y = 0; y < g.Ny; ++y) for (std::size_t x = 0; x < g.Nx; ++x)
        output.data[p*g.planeStride+y*g.yStride+x*g.xStride] = w.real[p*d.planeSize+y*g.Nx+x];
    return WVKernelStatus::ok();
}
} // namespace wavevortex
