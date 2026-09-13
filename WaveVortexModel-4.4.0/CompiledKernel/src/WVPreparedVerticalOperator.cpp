#include "WVSpectralValidation.hpp"
#include "WVPreparedModeExecutor.hpp"
#include <cstring>
#include <new>
#include <set>
#include <system_error>

namespace wavevortex {
namespace spectral_detail {
struct PreparedMatrix {
    WVMatrixIdentity identity;
    std::vector<double> real;
    std::vector<WVComplex64> complex;
};
struct PreparedGroup {
    std::uint64_t identity;
    std::size_t matrix;
    std::vector<std::size_t> modes;
    bool direct;
};
struct VerticalData {
    WVComplexLayout input, output;
    WVMatrixAction action;
    WVAccumulation accumulation;
    std::size_t Nz, Nj, inputSpan, outputSpan, packedColumns = 0;
    std::vector<PreparedMatrix> matrices;
    std::vector<PreparedGroup> groups;
    std::vector<std::size_t> columnMatrices;
    std::unique_ptr<WVVerticalMatrixBackend> backend;
};
struct VerticalWorkspaceData {
    std::shared_ptr<const VerticalData> owner;
    std::vector<double> br, bi, cr, ci;
    std::vector<WVComplex64> b, c;
    std::size_t groupWorkers = 1;
    std::atomic<bool> active{false};
};
struct VerticalGroupExecutorData {
    explicit VerticalGroupExecutorData(std::size_t workers):executor(workers) {}
    kernel_detail::WVPreparedModeExecutor executor;
    std::atomic<bool> active{false};
};
class ScalarBackend final : public WVVerticalMatrixBackend {
public:
    const char* identifier() const noexcept override { return "scalar-reference-v1"; }
    std::size_t maximumDimension() const noexcept override { return PTRDIFF_MAX; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    bool supportsConcurrentCalls() const noexcept override { return true; }
    void split(std::size_t m, std::size_t k, std::size_t n, const double* a, const double* br, const double* bi,
        std::size_t ldb, double* cr, double* ci, std::size_t ldc, double beta) const noexcept override {
        for (std::size_t col = 0; col < n; ++col) for (std::size_t row = 0; row < m; ++row) {
            double real = 0, imag = 0;
            for (std::size_t j = 0; j < k; ++j) { real += a[row+m*j]*br[j+ldb*col]; imag += a[row+m*j]*bi[j+ldb*col]; }
            const auto i = row+ldc*col;
            // beta=0 must not read uninitialized/NaN destination values.
            cr[i] = beta == 0 ? real : real+beta*cr[i]; ci[i] = beta == 0 ? imag : imag+beta*ci[i];
        }
    }
    void interleaved(std::size_t m, std::size_t k, std::size_t n, const WVComplex64* a, const WVComplex64* b,
        std::size_t ldb, WVComplex64* c, std::size_t ldc, double beta) const noexcept override {
        for (std::size_t col = 0; col < n; ++col) for (std::size_t row = 0; row < m; ++row) {
            WVComplex64 value{};
            for (std::size_t j = 0; j < k; ++j) {
                const auto aa = a[row+m*j], bb = b[j+ldb*col];
                value.real += aa.real*bb.real-aa.imag*bb.imag;
                value.imag += aa.real*bb.imag+aa.imag*bb.real;
            }
            const auto i = row+ldc*col;
            if (beta != 0) { value.real += beta*c[i].real; value.imag += beta*c[i].imag; }
            c[i] = value;
        }
    }
};
} // namespace spectral_detail
using namespace spectral_detail;
WVKernelStatus WVCreateScalarMatrixBackend(std::unique_ptr<WVVerticalMatrixBackend>& result) {
    try { result = std::make_unique<ScalarBackend>(); return WVKernelStatus::ok(); }
    catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Scalar backend allocation failed."}; }
}
WVVerticalGroupExecutor::WVVerticalGroupExecutor() = default;
WVVerticalGroupExecutor::~WVVerticalGroupExecutor() = default;
WVKernelStatus WVVerticalGroupExecutor::create(std::size_t workers,
    std::unique_ptr<WVVerticalGroupExecutor>& result) {
    try {
        if (!workers)
            return {WVKernelStatusCode::invalidConfiguration,"A vertical group executor requires a worker."};
        auto executor=std::unique_ptr<WVVerticalGroupExecutor>(new WVVerticalGroupExecutor);
        executor->data_=std::make_unique<VerticalGroupExecutorData>(workers);
        result=std::move(executor);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Vertical group executor allocation failed."};
    } catch (const std::system_error& e) {
        return {WVKernelStatusCode::allocationFailure,e.what()};
    } catch (const std::length_error& e) {
        return {WVKernelStatusCode::sizeOverflow,e.what()};
    }
}
std::size_t WVVerticalGroupExecutor::workerCount() const noexcept {
    return data_->executor.workerCount();
}
std::size_t WVVerticalGroupExecutor::persistentBytes() const noexcept {
    return sizeof(*this)+sizeof(*data_)+data_->executor.persistentBytes()-sizeof(data_->executor);
}
WVVerticalWorkspace::WVVerticalWorkspace() = default;
WVVerticalWorkspace::~WVVerticalWorkspace() = default;
std::size_t WVVerticalWorkspace::persistentBytes() const noexcept {
    return sizeof(*this)+sizeof(*data_)+(data_->br.capacity()+data_->bi.capacity()+data_->cr.capacity()+data_->ci.capacity())*sizeof(double)+
        (data_->b.capacity()+data_->c.capacity())*sizeof(WVComplex64);
}
std::size_t WVPreparedVerticalOperator::matrixBytes() const noexcept {
    std::size_t bytes = 0;
    for (const auto& m : data_->matrices) bytes += m.real.capacity()*sizeof(double)+m.complex.capacity()*sizeof(WVComplex64);
    return bytes;
}
std::size_t WVPreparedVerticalOperator::persistentBytes() const noexcept {
    std::size_t bytes = sizeof(*this)+sizeof(*data_)+identityBytes(data_->input)+identityBytes(data_->output)+matrixBytes()+
        data_->matrices.capacity()*sizeof(PreparedMatrix)+data_->groups.capacity()*sizeof(PreparedGroup)+
        data_->columnMatrices.capacity()*sizeof(std::size_t)+data_->backend->persistentBytes();
    for (const auto& m : data_->matrices) bytes += m.identity.source.capacity()+m.identity.name.capacity();
    for (const auto& g : data_->groups) bytes += g.modes.capacity()*sizeof(std::size_t);
    return bytes;
}
std::size_t WVPreparedVerticalOperator::uniqueMatrixCount() const noexcept { return data_->matrices.size(); }
std::size_t WVPreparedVerticalOperator::preparedGroupCount() const noexcept { return data_->groups.size(); }
bool WVPreparedVerticalOperator::supportsConcurrentCalls() const noexcept {
    return data_->backend->supportsConcurrentCalls();
}
const char* WVPreparedVerticalOperator::backendIdentifier() const noexcept { return data_->backend->identifier(); }
WVKernelStatus WVPreparedVerticalOperator::create(const WVVerticalSpecification& spec,
    std::unique_ptr<WVVerticalMatrixBackend> backend, std::unique_ptr<WVPreparedVerticalOperator>& result) {
    try {
        if (!backend) return {WVKernelStatusCode::invalidPointer,"A matrix backend is required."};
        if (spec.placement != WVOperatorPlacement::outOfPlace) return {WVKernelStatusCode::unsupportedOperation,"Vertical execution is out-of-place."};
        if (spec.sourceIdentity.empty() || !spec.Nz || !spec.Nj || spec.matrices.empty() || spec.groups.empty())
            return {WVKernelStatusCode::invalidConfiguration,"Vertical source, basis, matrices and groups are required."};
        const auto inputSpan = validate(spec.input), outputSpan = validate(spec.output);
        if (spec.input.columns != spec.output.columns || spec.input.modeSet != spec.output.modeSet || spec.input.representation != spec.output.representation)
            return {WVKernelStatusCode::invalidConfiguration,"Vertical views must share representation and exact ordered modes."};
        if (spec.accumulation != WVAccumulation::overwrite && spec.accumulation != WVAccumulation::add)
            return {WVKernelStatusCode::invalidConfiguration,"Unknown accumulation policy."};
        std::size_t inRows = 0, outRows = 0;
        switch (spec.action) {
            case WVMatrixAction::reconstruction: inRows = spec.Nj; outRows = spec.Nz; break;
            case WVMatrixAction::projection: inRows = spec.Nz; outRows = spec.Nj; break;
            case WVMatrixAction::crossFamily: inRows = outRows = spec.Nj; break;
            default: return {WVKernelStatusCode::invalidConfiguration,"Unknown matrix action."};
        }
        if (spec.input.rows != inRows || spec.output.rows != outRows)
            return {WVKernelStatusCode::invalidShape,"Vertical matrix action and basis extents disagree."};
        const auto limit = backend->maximumDimension();
        if (inRows > limit || outRows > limit)
            return {WVKernelStatusCode::sizeOverflow,"Matrix dimensions exceed the provider integer range."};
        product(product(inRows,outRows),sizeof(WVComplex64));
        auto d = std::make_shared<VerticalData>();
        d->input = spec.input; d->output = spec.output; d->action = spec.action; d->accumulation = spec.accumulation;
        d->Nz = spec.Nz; d->Nj = spec.Nj; d->inputSpan = inputSpan; d->outputSpan = outputSpan; d->backend = std::move(backend);
        d->columnMatrices.resize(spec.input.columns);
        std::vector<std::size_t> matrixIndices;
        for (const auto& source : spec.matrices) {
            const auto& m = source.values;
            if (source.action != spec.action || source.inputFamily != spec.input.family || source.outputFamily != spec.output.family)
                return {WVKernelStatusCode::invalidConfiguration,"Matrix action or family identity disagrees with the requested operation."};
            if (source.identity.source != spec.sourceIdentity || source.identity.name.empty())
                return {WVKernelStatusCode::invalidConfiguration,"Matrix source identity is inconsistent."};
            if (m.rows != outRows || m.columns != inRows) return {WVKernelStatusCode::invalidShape,"Matrix shape disagrees with field layouts."};
            const auto bytes = span(m.rows,m.columns,m.rowStride,m.columnStride,sizeof(double));
            if (!addressFits(m.data,bytes,alignof(double))) return {WVKernelStatusCode::invalidPointer,"Invalid matrix storage."};
            if (m.bytes < bytes) return {WVKernelStatusCode::invalidShape,"Matrix buffer capacity is too small."};
            std::size_t existing = 0;
            while (existing < d->matrices.size() && d->matrices[existing].identity.name != source.identity.name) ++existing;
            // Identity equality is exact. Compare every value (including signed
            // zero) before reusing an identity, without approximate deduplication.
            for (std::size_t c = 0; c < inRows; ++c) for (std::size_t r = 0; r < outRows; ++r) {
                const auto value = m.data[r*m.rowStride+c*m.columnStride];
                if (!std::isfinite(value)) return {WVKernelStatusCode::invalidConfiguration,"Matrix contains a nonfinite value."};
                if (existing < d->matrices.size()) {
                    const auto& old = d->matrices[existing];
                    const auto saved = old.real.empty() ? old.complex[r+outRows*c].real : old.real[r+outRows*c];
                    if (std::memcmp(&value,&saved,sizeof(double)) != 0)
                        return {WVKernelStatusCode::invalidConfiguration,"One exact matrix identity has different values."};
                }
            }
            if (existing == d->matrices.size()) {
                PreparedMatrix matrix; matrix.identity = source.identity;
                if (spec.input.representation == WVComplexRepresentation::split) matrix.real.resize(product(inRows,outRows));
                else matrix.complex.resize(product(inRows,outRows));
                for (std::size_t c = 0; c < inRows; ++c) for (std::size_t r = 0; r < outRows; ++r) {
                    const auto value = m.data[r*m.rowStride+c*m.columnStride];
                    if (matrix.real.empty()) matrix.complex[r+outRows*c] = {value,0}; else matrix.real[r+outRows*c] = value;
                }
                d->matrices.push_back(std::move(matrix));
            }
            matrixIndices.push_back(existing);
        }
        std::vector<bool> visited(spec.input.columns,false), usedMatrices(spec.matrices.size(),false);
        std::set<std::uint64_t> groups;
        const bool directLayout = spec.input.rowStride == 1 && spec.output.rowStride == 1 &&
            spec.input.columnStride >= inRows && spec.output.columnStride >= outRows &&
            spec.input.columnStride <= limit && spec.output.columnStride <= limit;
        for (const auto& source : spec.groups) {
            if (!groups.insert(source.identity).second || source.matrix >= matrixIndices.size() || source.modes.empty())
                return {WVKernelStatusCode::invalidConfiguration,"Invalid or duplicate vertical group."};
            if (source.modes.size() > limit) return {WVKernelStatusCode::sizeOverflow,"Group width exceeds the provider integer range."};
            for (std::size_t i = 0; i < source.modes.size(); ++i) {
                const auto mode = source.modes[i];
                if (mode >= visited.size() || visited[mode]) return {WVKernelStatusCode::invalidConfiguration,"Vertical group membership is invalid or repeated."};
                visited[mode] = true;
                d->columnMatrices[mode] = matrixIndices[source.matrix];
            }
            usedMatrices[source.matrix] = true;
            if (!directLayout) {
                d->packedColumns = std::max(d->packedColumns,source.modes.size());
                d->groups.push_back({source.identity,matrixIndices[source.matrix],source.modes,false});
                continue;
            }
            // Preserve the exact public group and matrix-source mapping while
            // presenting each maximal contiguous run as a direct matrix view.
            // Membership order remains authoritative; no sorting or inferred
            // grouping key participates in preparation.
            for (std::size_t first = 0; first < source.modes.size();) {
                std::size_t end = first+1;
                while (end < source.modes.size() && source.modes[end] == source.modes[end-1]+1) ++end;
                d->groups.push_back({source.identity,matrixIndices[source.matrix],
                    {source.modes.begin()+static_cast<std::ptrdiff_t>(first),source.modes.begin()+static_cast<std::ptrdiff_t>(end)},true});
                first = end;
            }
        }
        if (std::find(visited.begin(),visited.end(),false) != visited.end()) return {WVKernelStatusCode::invalidConfiguration,"Vertical groups do not cover the retained set."};
        if (std::find(usedMatrices.begin(),usedMatrices.end(),false) != usedMatrices.end()) return {WVKernelStatusCode::invalidConfiguration,"Unused matrix records are not permitted."};
        product(product(d->packedColumns,inRows),sizeof(WVComplex64)); product(product(d->packedColumns,outRows),sizeof(WVComplex64));
        auto op = std::unique_ptr<WVPreparedVerticalOperator>(new WVPreparedVerticalOperator);
        op->data_ = std::move(d); result = std::move(op); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Vertical setup allocation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
      catch (const std::invalid_argument& e) { return {WVKernelStatusCode::invalidConfiguration,e.what()}; }
}
WVKernelStatus WVPreparedVerticalOperator::createWorkspace(std::unique_ptr<WVVerticalWorkspace>& result) const {
    return createWorkspace(1,result);
}
WVKernelStatus WVPreparedVerticalOperator::createWorkspace(std::size_t groupWorkers,
    std::unique_ptr<WVVerticalWorkspace>& result) const {
    try {
        if (!groupWorkers)
            return {WVKernelStatusCode::invalidConfiguration,"Vertical group worker count must be positive."};
        auto w = std::unique_ptr<WVVerticalWorkspace>(new WVVerticalWorkspace);
        w->data_ = std::make_unique<VerticalWorkspaceData>(); auto& ws = *w->data_; ws.owner = data_; ws.groupWorkers=groupWorkers;
        const auto b = product(product(data_->input.rows,data_->packedColumns),groupWorkers);
        const auto c = product(product(data_->output.rows,data_->packedColumns),groupWorkers);
        if (data_->input.representation == WVComplexRepresentation::interleaved) { ws.b.resize(b); ws.c.resize(c); }
        else { ws.br.resize(b); ws.bi.resize(b); ws.cr.resize(c); ws.ci.resize(c); }
        result = std::move(w); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Vertical workspace allocation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
}
static void executePreparedGroups(const VerticalData& d,VerticalWorkspaceData& w,
    WVComplexInput input,WVComplexOutput output,double beta,std::size_t worker,
    std::size_t begin,std::size_t end) {
    const auto m = d.output.rows, k = d.input.rows;
    const auto packedBStride=k*d.packedColumns,packedCStride=m*d.packedColumns;
    for (std::size_t groupIndex=begin;groupIndex<end;++groupIndex) {
        const auto& group=d.groups[groupIndex];
        const auto& a = d.matrices[group.matrix]; const auto n = group.modes.size();
        const auto bOffset = group.modes.front()*d.input.columnStride, cOffset = group.modes.front()*d.output.columnStride;
        if (group.direct) {
            if (d.input.representation == WVComplexRepresentation::interleaved)
                d.backend->interleaved(m,k,n,a.complex.data(),input.interleaved+bOffset,d.input.columnStride,output.interleaved+cOffset,d.output.columnStride,beta);
            else d.backend->split(m,k,n,a.real.data(),input.real+bOffset,input.imag+bOffset,d.input.columnStride,output.real+cOffset,output.imag+cOffset,d.output.columnStride,beta);
        } else {
            const bool split = d.input.representation == WVComplexRepresentation::split;
            WVComplexOutput packedB{split ? nullptr : w.b.data()+worker*packedBStride,
                split ? w.br.data()+worker*packedBStride : nullptr,
                split ? w.bi.data()+worker*packedBStride : nullptr,0};
            WVComplexOutput packedC{split ? nullptr : w.c.data()+worker*packedCStride,
                split ? w.cr.data()+worker*packedCStride : nullptr,
                split ? w.ci.data()+worker*packedCStride : nullptr,0};
            for (std::size_t col = 0; col < n; ++col) for (std::size_t row = 0; row < k; ++row)
                write(packedB,row+k*col,read(input,row*d.input.rowStride+group.modes[col]*d.input.columnStride));
            if (split) d.backend->split(m,k,n,a.real.data(),packedB.real,packedB.imag,k,
                packedC.real,packedC.imag,m,0);
            else d.backend->interleaved(m,k,n,a.complex.data(),packedB.interleaved,k,
                packedC.interleaved,m,0);
            for (std::size_t col = 0; col < n; ++col) for (std::size_t row = 0; row < m; ++row) {
                const auto target = row*d.output.rowStride+group.modes[col]*d.output.columnStride;
                auto value = read(packedC.input(),row+m*col);
                if (beta != 0) { const auto old = read(output.input(),target); value.real += old.real; value.imag += old.imag; }
                write(output,target,value);
            }
        }
    }
}
WVKernelStatus WVPreparedVerticalOperator::execute(WVVerticalWorkspace& workspace, WVComplexInput input, WVComplexOutput output) const {
    auto& w = *workspace.data_; const auto& d = *data_;
    if (w.owner != data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another vertical operator."};
    auto status = validateStorage(d.input.representation,d.inputSpan,input); if (!status) return status;
    status = validateStorage(d.output.representation,d.outputSpan,output.input()); if (!status) return status;
    if (storageOverlap(input,d.inputSpan,output.input(),d.outputSpan)) return {WVKernelStatusCode::overlappingArrays,"Vertical input and output overlap."};
    ActiveCall guard(w.active); if (!guard.entered) return {WVKernelStatusCode::reentrantExecution,"Vertical workspace is active."};
    executePreparedGroups(d,w,input,output,d.accumulation==WVAccumulation::add ? 1 : 0,
        0,0,d.groups.size());
    return WVKernelStatus::ok();
}
WVKernelStatus WVPreparedVerticalOperator::execute(WVVerticalWorkspace& workspace,
    WVVerticalGroupExecutor& executor,WVComplexInput input,WVComplexOutput output) const {
    auto& w=*workspace.data_; const auto& d=*data_; auto& e=*executor.data_;
    if (w.owner!=data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another vertical operator."};
    if (!d.backend->supportsConcurrentCalls())
        return {WVKernelStatusCode::unsupportedOperation,"Vertical matrix backend does not permit concurrent calls."};
    if (w.groupWorkers!=e.executor.workerCount())
        return {WVKernelStatusCode::invalidConfiguration,"Vertical workspace and group executor worker counts differ."};
    auto status=validateStorage(d.input.representation,d.inputSpan,input); if (!status) return status;
    status=validateStorage(d.output.representation,d.outputSpan,output.input()); if (!status) return status;
    if (storageOverlap(input,d.inputSpan,output.input(),d.outputSpan))
        return {WVKernelStatusCode::overlappingArrays,"Vertical input and output overlap."};
    ActiveCall workspaceGuard(w.active); if (!workspaceGuard.entered)
        return {WVKernelStatusCode::reentrantExecution,"Vertical workspace is active."};
    ActiveCall executorGuard(e.active); if (!executorGuard.entered)
        return {WVKernelStatusCode::reentrantExecution,"Vertical group executor is active."};
    const auto workers=e.executor.workerCount();
    const auto block=d.groups.size()/workers+(d.groups.size()%workers!=0);
    const double beta=d.accumulation==WVAccumulation::add ? 1 : 0;
    e.executor.execute(d.groups.size(),[&](std::size_t begin,std::size_t end) {
        if (begin==end) return;
        executePreparedGroups(d,w,input,output,beta,begin/block,begin,end);
    });
    return WVKernelStatus::ok();
}
WVKernelStatus WVPreparedVerticalOperator::executeColumn(WVVerticalWorkspace& workspace,
    WVComplexInput input,WVComplexOutput output,std::size_t retainedColumn) const {
    auto& w=*workspace.data_; const auto& d=*data_;
    if (w.owner!=data_) return {WVKernelStatusCode::invalidConfiguration,"Workspace belongs to another vertical operator."};
    auto status=validateStorage(d.input.representation,d.inputSpan,input); if (!status) return status;
    status=validateStorage(d.output.representation,d.outputSpan,output.input()); if (!status) return status;
    if (storageOverlap(input,d.inputSpan,output.input(),d.outputSpan))
        return {WVKernelStatusCode::overlappingArrays,"Vertical input and output overlap."};
    if (retainedColumn>=d.input.columns)
        return {WVKernelStatusCode::invalidShape,"Vertical retained column is out of range."};
    ActiveCall guard(w.active); if (!guard.entered)
        return {WVKernelStatusCode::reentrantExecution,"Vertical workspace is active."};
    const auto m=d.output.rows,k=d.input.rows;
    const double beta=d.accumulation==WVAccumulation::add ? 1 : 0;
    const auto& a=d.matrices[d.columnMatrices[retainedColumn]];
    const auto bOffset=retainedColumn*d.input.columnStride;
    const auto cOffset=retainedColumn*d.output.columnStride;
    const auto limit=d.backend->maximumDimension();
    const bool direct=d.input.rowStride==1 && d.output.rowStride==1 &&
        d.input.columnStride>=k && d.output.columnStride>=m &&
        d.input.columnStride<=limit && d.output.columnStride<=limit;
    if (direct) {
        if (d.input.representation==WVComplexRepresentation::interleaved)
            d.backend->interleaved(m,k,1,a.complex.data(),input.interleaved+bOffset,
                d.input.columnStride,output.interleaved+cOffset,d.output.columnStride,beta);
        else d.backend->split(m,k,1,a.real.data(),input.real+bOffset,input.imag+bOffset,
            d.input.columnStride,output.real+cOffset,output.imag+cOffset,d.output.columnStride,beta);
        return WVKernelStatus::ok();
    }
    const bool split=d.input.representation==WVComplexRepresentation::split;
    WVComplexOutput packedB{split ? nullptr : w.b.data(),split ? w.br.data() : nullptr,
        split ? w.bi.data() : nullptr,0};
    WVComplexInput packedC{split ? nullptr : w.c.data(),split ? w.cr.data() : nullptr,
        split ? w.ci.data() : nullptr,0};
    for (std::size_t row=0;row<k;++row)
        write(packedB,row,read(input,row*d.input.rowStride+bOffset));
    if (split) d.backend->split(m,k,1,a.real.data(),w.br.data(),w.bi.data(),k,
        w.cr.data(),w.ci.data(),m,0);
    else d.backend->interleaved(m,k,1,a.complex.data(),w.b.data(),k,w.c.data(),m,0);
    for (std::size_t row=0;row<m;++row) {
        const auto target=row*d.output.rowStride+cOffset;
        auto value=read(packedC,row);
        if (beta!=0) { const auto old=read(output.input(),target); value.real+=old.real; value.imag+=old.imag; }
        write(output,target,value);
    }
    return WVKernelStatus::ok();
}
} // namespace wavevortex
