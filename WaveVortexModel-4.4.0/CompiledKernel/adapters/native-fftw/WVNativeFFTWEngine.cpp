#include "WVNativeFFTWEngine.hpp"
#include "WVNativeAdvectionArithmetic.hpp"
#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <condition_variable>
#include <stdexcept>
#include <system_error>
#include <thread>

#include <fftw3.h>

#include <climits>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <limits>
#include <new>
#include <utility>
#include <vector>

#if defined(__APPLE__) || defined(__linux__)
#include <dlfcn.h>
#endif

namespace wavevortex {
namespace {

std::mutex planningMutex;
std::atomic<std::size_t> activePlans{0};
std::atomic<std::size_t> totalPlansCreated{0};
std::atomic<std::size_t> totalPlansDestroyed{0};
std::atomic<std::size_t> outstandingPlanningBytes{0};
std::atomic<std::uint64_t> totalPlanningNanoseconds{0};

WVKernelStatus dimensions(const std::vector<WVFFTDimension>& source, std::vector<fftw_iodim64>& destination) {
    destination.clear();
    destination.reserve(source.size());
    for (const auto& dimension : source) {
        if (dimension.count == 0 || dimension.count > static_cast<std::size_t>(PTRDIFF_MAX) ||
            dimension.inputStride <= 0 || dimension.outputStride <= 0) {
            return {WVKernelStatusCode::invalidConfiguration, "FFTW dimensions must be nonzero, fit ptrdiff_t, and have positive strides."};
        }
        destination.push_back({static_cast<ptrdiff_t>(dimension.count),dimension.inputStride,dimension.outputStride});
    }
    return WVKernelStatus::ok();
}

std::string libraryContaining(const void* symbol) {
#if defined(__APPLE__) || defined(__linux__)
    Dl_info information{};
    if (symbol != nullptr && dladdr(symbol,&information) != 0 && information.dli_fname != nullptr) return information.dli_fname;
#endif
    return {};
}

std::string loadedLibrary() { return libraryContaining(reinterpret_cast<const void*>(&fftw_execute)); }

// Validate actual addressed spans, not arena capacities: a vertical sine plan
// starts inside a larger shared arena and intentionally leaves boundary holes.
WVKernelStatus addressedBytes(const WVFFTPlanSpecification& spec, bool input, std::size_t& bytes) {
    const bool complex = input ? spec.kind == WVFFTPlanKind::horizontalComplexToReal2D
                               : spec.kind == WVFFTPlanKind::horizontalRealToComplex2D;
    const auto elementBytes = complex ? sizeof(fftw_complex) : sizeof(double);
    const auto maximum = static_cast<std::size_t>(PTRDIFF_MAX) / elementBytes;
    std::size_t elements = 1;
    auto add = [&](const WVFFTDimension& dimension, bool lastTransform) {
        const auto count = complex && lastTransform ? dimension.count / 2 + 1 : dimension.count;
        const auto stride = static_cast<std::size_t>(input ? dimension.inputStride : dimension.outputStride);
        if (count - 1 > (maximum - elements) / stride) return false;
        elements += (count - 1) * stride;
        return true;
    };
    for (std::size_t i = 0; i < spec.transformDimensions.size(); ++i)
        if (!add(spec.transformDimensions[i], i + 1 == spec.transformDimensions.size()))
            return {WVKernelStatusCode::sizeOverflow,"FFTW transform span overflows ptrdiff_t."};
    for (const auto& dimension : spec.batchDimensions)
        if (!add(dimension, false)) return {WVKernelStatusCode::sizeOverflow,"FFTW batch span overflows ptrdiff_t."};
    bytes = elements * elementBytes;
    if (bytes > (input ? spec.inputBytes : spec.outputBytes))
        return {WVKernelStatusCode::invalidShape,"FFTW addressed span exceeds the declared buffer size."};
    return WVKernelStatus::ok();
}

class PlanningBuffer final {
public:
    explicit PlanningBuffer(std::size_t bytes) : data_(fftw_malloc(bytes)), bytes_(bytes) {
        if (!data_) throw std::bad_alloc();
        outstandingPlanningBytes += bytes_;
    }
    ~PlanningBuffer() { fftw_free(data_); outstandingPlanningBytes -= bytes_; }
    PlanningBuffer(const PlanningBuffer&) = delete;
    PlanningBuffer& operator=(const PlanningBuffer&) = delete;
    void* get() const noexcept { return data_; }
private:
    void* data_;
    std::size_t bytes_;
};

struct PlanDeleter {
    void operator()(fftw_plan plan) const noexcept {
        std::lock_guard<std::mutex> lock(planningMutex);
        fftw_destroy_plan(plan);
        --activePlans;
        ++totalPlansDestroyed;
    }
};
using OwnedPlan = std::unique_ptr<fftw_plan_s, PlanDeleter>;

class FFTWPlan final : public WVFFTPlan {
public:
    FFTWPlan(OwnedPlan plan, WVFFTPlanKind kind, bool inPlace, std::size_t inputSpan, std::size_t outputSpan)
        : plan_(std::move(plan)), kind_(kind), inPlace_(inPlace), inputSpan_(inputSpan), outputSpan_(outputSpan) {}
    WVKernelStatus execute(const void* input, void* output) override {
        if (input == nullptr || output == nullptr) return {WVKernelStatusCode::invalidPointer,"FFTW received a null execution pointer."};
        const auto in = reinterpret_cast<std::uintptr_t>(input);
        const auto out = reinterpret_cast<std::uintptr_t>(output);
        if (inPlace_) {
            if (in != out) return {WVKernelStatusCode::invalidPointer,"An in-place FFTW plan requires identical pointers."};
        } else if (in <= out ? out - in < inputSpan_ : in - out < outputSpan_) {
            return {WVKernelStatusCode::overlappingArrays,"An out-of-place FFTW plan requires disjoint addressed spans."};
        }
        switch (kind_) {
            case WVFFTPlanKind::horizontalRealToComplex2D:
                fftw_execute_dft_r2c(plan_.get(),const_cast<double*>(static_cast<const double*>(input)),static_cast<fftw_complex*>(output));
                break;
            case WVFFTPlanKind::horizontalComplexToReal2D:
                fftw_execute_dft_c2r(plan_.get(),const_cast<fftw_complex*>(static_cast<const fftw_complex*>(input)),static_cast<double*>(output));
                break;
            case WVFFTPlanKind::verticalDCTI:
            case WVFFTPlanKind::verticalDSTI:
                fftw_execute_r2r(plan_.get(),const_cast<double*>(static_cast<const double*>(input)),static_cast<double*>(output));
                break;
        }
        return WVKernelStatus::ok();
    }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
private:
    OwnedPlan plan_;
    WVFFTPlanKind kind_;
    bool inPlace_;
    std::size_t inputSpan_, outputSpan_;
};

// Prepared outer parallelism; execution never creates threads or allocates.
class RetainedWorkers {
public:
    explicit RetainedWorkers(std::size_t count) : count_(count) {
        threads_.reserve(count-1);
        try {
            for (std::size_t i=1; i<count; ++i) threads_.emplace_back([this,i] { loop(i); });
        } catch (...) { stop(); throw; }
    }
    ~RetainedWorkers() { stop(); }
    void run(void (*task)(void*,std::size_t) noexcept, void* context) noexcept {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            task_=task; context_=context; remaining_=count_-1; ++generation_;
        }
        ready_.notify_all();
        task(context,0);
        std::unique_lock<std::mutex> lock(mutex_);
        complete_.wait(lock,[this] { return remaining_==0; });
    }
    std::size_t bytes() const noexcept { return sizeof(*this)+threads_.capacity()*sizeof(std::thread); }
private:
    void stop() noexcept {
        { std::lock_guard<std::mutex> lock(mutex_); stopping_=true; }
        ready_.notify_all();
        for (auto& thread:threads_) if (thread.joinable()) thread.join();
    }
    void loop(std::size_t worker) noexcept {
        std::size_t seen=0;
        for (;;) {
            std::unique_lock<std::mutex> lock(mutex_);
            ready_.wait(lock,[&] { return stopping_ || generation_!=seen; });
            if (stopping_) return;
            seen=generation_; auto task=task_; auto context=context_;
            lock.unlock(); task(context,worker); lock.lock();
            if (--remaining_==0) complete_.notify_one();
        }
    }
    std::size_t count_,generation_=0,remaining_=0;
    bool stopping_=false;
    void (*task_)(void*,std::size_t) noexcept = nullptr;
    void* context_=nullptr;
    std::mutex mutex_;
    std::condition_variable ready_,complete_;
    std::vector<std::thread> threads_;
};

// FFT dimensions and the contiguous selected-column interval completely
// determine this shared resource. Caller ordering/strides remain plan-local.
class RetainedFFTWResources final {
public:
    RetainedFFTWResources(std::size_t nx,std::size_t ny,std::size_t columns,std::size_t count,std::size_t workers)
        : Nx(nx),Ny(ny),activeColumns(columns),modeCount(count),workerCount(workers),pool(workers) {
        const auto scratchCount=checked(workers,checked(nx/2+1,ny));
        const auto tileCount=checked(checked(workers,16),count);
        checked(scratchCount,sizeof(WVComplex64)); checked(tileCount,sizeof(WVComplex64));
        scratch.resize(scratchCount); tile.resize(tileCount);
    }
    WVKernelStatus prepare() {
        const auto nx=Nx,ny=Ny,half=nx/2+1;
        PlanningBuffer real(checked(checked(nx,ny),sizeof(double)));
        auto* fftScratch=reinterpret_cast<fftw_complex*>(scratch.data());
        fftw_iodim64 row{static_cast<ptrdiff_t>(nx),1,1};
        fftw_iodim64 rowBatch{static_cast<ptrdiff_t>(ny),static_cast<ptrdiff_t>(nx),static_cast<ptrdiff_t>(half)};
        fftw_iodim64 column{static_cast<ptrdiff_t>(ny),static_cast<ptrdiff_t>(half),static_cast<ptrdiff_t>(half)};
        fftw_iodim64 columnBatch{static_cast<ptrdiff_t>(activeColumns),1,1};
        // Owned plans are installed while locked but destroyed only after the
        // lock is released, including partial planning failure.
        bool success=true;
        {
            std::lock_guard<std::mutex> lock(planningMutex);
            const auto start=std::chrono::steady_clock::now();
            fftw_plan_with_nthreads(1);
            auto take=[&](OwnedPlan& owner,fftw_plan raw) {
                if (raw) { ++activePlans; ++totalPlansCreated; owner.reset(raw); }
                else success=false;
            };
            const unsigned flags=FFTW_MEASURE|FFTW_UNALIGNED;
            take(rowForward_,fftw_plan_guru64_dft_r2c(1,&row,1,&rowBatch,static_cast<double*>(real.get()),fftScratch,flags|FFTW_PRESERVE_INPUT));
            take(columnForward_,fftw_plan_guru64_dft(1,&column,1,&columnBatch,fftScratch,fftScratch,FFTW_FORWARD,flags));
            take(columnInverse_,fftw_plan_guru64_dft(1,&column,1,&columnBatch,fftScratch,fftScratch,FFTW_BACKWARD,flags));
            std::swap(rowBatch.is,rowBatch.os);
            take(rowInverse_,fftw_plan_guru64_dft_c2r(1,&row,1,&rowBatch,fftScratch,static_cast<double*>(real.get()),flags));
            totalPlanningNanoseconds+=static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count());
        }
        return success ? WVKernelStatus::ok() : WVKernelStatus{WVKernelStatusCode::fftPlanFailure,"Unable to prepare streaming retained FFTW plans."};
    }
    std::size_t bytes() const noexcept {
        return sizeof(*this)+pool.bytes()-sizeof(pool)+(scratch.capacity()+tile.capacity())*sizeof(WVComplex64);
    }
    static std::size_t checked(std::size_t a,std::size_t b) {
        if (a && b>static_cast<std::size_t>(PTRDIFF_MAX)/a) throw std::overflow_error("Retained FFTW extent overflow.");
        return a*b;
    }
    const std::size_t Nx,Ny,activeColumns,modeCount,workerCount;
    RetainedWorkers pool;
    std::vector<WVComplex64> scratch,tile;
    OwnedPlan rowForward_,columnForward_,columnInverse_,rowInverse_;
    std::atomic<bool> active{false};
};

class RetainedFFTWPlan final : public WVRetainedHorizontalPlan {
    struct Mode {
        std::size_t row,partner;
        bool conjugate,self;
        double k,l;
    };
    static constexpr std::size_t tileWidth=16;
public:
    explicit RetainedFFTWPlan(const WVRetainedHorizontalSpecification& spec)
        : spec_(spec), workers_(std::min(spec.outerWorkers,spec.grid.planes)),
          halfSize_((spec.grid.Nx/2+1)*spec.grid.Ny),
          xWavenumberScale_(2*std::acos(-1.0)/spec.Lx),
          yWavenumberScale_(2*std::acos(-1.0)/spec.Ly) {
        const auto nx=spec.grid.Nx,ny=spec.grid.Ny,half=nx/2+1;
        const auto index=[](std::int64_t k,std::size_t n) { return k<0 ? static_cast<std::size_t>(static_cast<std::int64_t>(n)+k) : static_cast<std::size_t>(k); };
        modes_.reserve(spec.modes.size());
        for (const auto& key:spec.modes) {
            const auto x=index(key.k,nx),y=index(key.l,ny),cx=(nx-x)%nx,cy=(ny-y)%ny;
            const bool conjugate=x>nx/2;
            const auto sx=conjugate?cx:x,sy=conjugate?cy:y,row=sx+half*sy;
            const bool boundary=sx==0 || (nx%2==0 && sx==nx/2);
            modes_.push_back({row,boundary?sx+half*((ny-sy)%ny):row,conjugate,x==cx&&y==cy,
                xWavenumberScale_*static_cast<double>(key.k),
                yWavenumberScale_*static_cast<double>(key.l)});
            activeColumns_=std::max(activeColumns_,sx+1);
            if ((nx%2==0 && sx==nx/2) || (ny%2==0 && y==ny/2)) hasNyquist_=true;
        }
        const double n=static_cast<double>(nx)*ny;
        if (spec.normalization==WVFourierNormalization::forwardUnit) forwardScale_=1/n;
        else if (spec.normalization==WVFourierNormalization::inverseUnit) inverseScale_=1/n;
        else forwardScale_=inverseScale_=1/std::sqrt(n);
    }
    WVKernelStatus prepare(std::vector<std::weak_ptr<RetainedFFTWResources>>& cache) {
        cache.erase(std::remove_if(cache.begin(),cache.end(),[](const auto& item) { return item.expired(); }),cache.end());
        for (const auto& item:cache) if (auto existing=item.lock()) {
            if (existing->Nx==spec_.grid.Nx && existing->Ny==spec_.grid.Ny &&
                existing->activeColumns==activeColumns_ && existing->modeCount==modes_.size() && existing->workerCount==workers_) {
                resources_=std::move(existing); return WVKernelStatus::ok();
            }
        }
        auto resources=std::make_shared<RetainedFFTWResources>(spec_.grid.Nx,spec_.grid.Ny,activeColumns_,modes_.size(),workers_);
        auto status=resources->prepare(); if (!status) return status;
        cache.push_back(resources); resources_=std::move(resources);
        return WVKernelStatus::ok();
    }
    WVKernelStatus forward(WVRealInput input,WVComplexOutput output) override {
        if (resources_->active.exchange(true)) return {WVKernelStatusCode::reentrantExecution,"Shared retained FFTW resource is active."};
        Context context{this,input,{}, {},output}; resources_->pool.run(forwardTask,&context);
        resources_->active.store(false); return WVKernelStatus::ok();
    }
    WVKernelStatus inverse(WVComplexInput input,WVRealOutput output) override {
        return inverseAndConsume(input,output,{});
    }
    bool supportsInverseConsumer() const noexcept override { return true; }
    WVKernelStatus inverseAndConsume(WVComplexInput input,WVRealOutput output,
        const WVRealOutputConsumer& consumer) override {
        if (resources_->active.exchange(true)) return {WVKernelStatusCode::reentrantExecution,"Shared retained FFTW resource is active."};
        Context context{this,{},output,input,{},&consumer}; resources_->pool.run(inverseTask,&context);
        resources_->active.store(false); return WVKernelStatus::ok();
    }
    bool supportsAdvection(std::size_t targets) const noexcept override {
        return (targets==3 || targets==4) && !hasNyquist_ &&
            spec_.grid.xStride==1 && spec_.grid.yStride==spec_.grid.Nx &&
            spec_.grid.planeStride==spec_.grid.Nx*spec_.grid.Ny;
    }
    WVKernelStatus prepareAdvection(std::size_t targets) override {
        if (!supportsAdvection(targets))
            return {WVKernelStatusCode::unsupportedOperation,"Retained FFTW advection requires contiguous planes without Nyquist modes."};
        if (resources_->active.exchange(true))
            return {WVKernelStatusCode::reentrantExecution,"Shared retained FFTW resource is active."};
        WVKernelStatus status=WVKernelStatus::ok();
        try {
            const auto maximumPlanes=(spec_.grid.planes/workers_)+
                (spec_.grid.planes%workers_!=0 ? 1 : 0);
            const auto depth=std::min(std::size_t{4},maximumPlanes);
            const auto tile=RetainedFFTWResources::checked(depth,modes_.size());
            const auto packedPerWorker=RetainedFFTWResources::checked(targets,tile);
            const auto stagePerWorker=RetainedFFTWResources::checked(
                targets,RetainedFFTWResources::checked(activeColumns_,spec_.grid.Ny));
            const auto realPerWorker=RetainedFFTWResources::checked(2,spec_.grid.Nx*spec_.grid.Ny);
            std::vector<WVComplex64> packed(
                RetainedFFTWResources::checked(workers_,packedPerWorker));
            std::vector<WVComplex64> stages(
                RetainedFFTWResources::checked(workers_,stagePerWorker));
            std::vector<double> real(RetainedFFTWResources::checked(workers_,realPerWorker));
            std::vector<WVRetainedAdvectionCounts> workerCounts(workers_);
            advectionPacked_.swap(packed); advectionStages_.swap(stages);
            advectionReal_.swap(real); advectionWorkerCounts_.swap(workerCounts);
            advectionTargets_=targets; advectionDepth_=depth;
        } catch (const std::bad_alloc&) {
            status={WVKernelStatusCode::allocationFailure,"Retained FFTW advection workspace allocation failed."};
        } catch (const std::length_error&) {
            status={WVKernelStatusCode::sizeOverflow,"Retained FFTW advection workspace exceeds container limits."};
        } catch (const std::overflow_error& e) {
            status={WVKernelStatusCode::sizeOverflow,e.what()};
        }
        resources_->active.store(false);
        return status;
    }
    WVKernelStatus advection(const WVRetainedAdvectionWork& work,
        WVRetainedAdvectionCounts& counts) override {
        if (!supportsAdvection(work.targets))
            return {WVKernelStatusCode::unsupportedOperation,"Retained FFTW advection is not supported for this plan."};
        if (advectionTargets_!=work.targets || !advectionDepth_)
            return {WVKernelStatusCode::invalidConfiguration,"Retained FFTW advection was not prepared for this target count."};
        if (resources_->active.exchange(true))
            return {WVKernelStatusCode::reentrantExecution,"Shared retained FFTW resource is active."};
        std::fill(advectionWorkerCounts_.begin(),advectionWorkerCounts_.end(),WVRetainedAdvectionCounts{});
        AdvectionContext context{this,&work};
        resources_->pool.run(advectionTask,&context);
        WVRetainedAdvectionCounts result;
        for (const auto& local:advectionWorkerCounts_) {
            result.columnInverses+=local.columnInverses;
            result.rowInverses+=local.rowInverses;
            result.reusedColumns+=local.reusedColumns;
        }
        resources_->active.store(false); counts=result;
        return WVKernelStatus::ok();
    }
    const char* identifier() const noexcept override { return "fftw-streaming-pruned-tile16"; }
    std::size_t workerCount() const noexcept override { return workers_; }
    std::size_t persistentBytes() const noexcept override {
        return sizeof(*this)+resources_->bytes()+modes_.capacity()*sizeof(Mode)+
            spec_.modes.capacity()*sizeof(WVRetainedModeKey)+spec_.grid.family.capacity()+
            spec_.retained.family.capacity()+spec_.retained.modeSet.capacity()+
            (advectionPacked_.capacity()+advectionStages_.capacity())*sizeof(WVComplex64)+
            advectionReal_.capacity()*sizeof(double)+
            advectionWorkerCounts_.capacity()*sizeof(WVRetainedAdvectionCounts);
    }
    // Shared resources include the four handles; opaque FFTW allocations and
    // OS worker stacks remain excluded, as for the full FFT adapter.
    std::size_t planBytesLowerBound() const noexcept override { return 0; }
    const void* sharedResourceIdentity() const noexcept override { return resources_.get(); }
    std::size_t sharedResourceBytes() const noexcept override { return resources_->bytes(); }
private:
    struct Context { RetainedFFTWPlan* plan; WVRealInput realInput; WVRealOutput realOutput; WVComplexInput complexInput; WVComplexOutput complexOutput; const WVRealOutputConsumer* consumer=nullptr; };
    struct AdvectionContext { RetainedFFTWPlan* plan; const WVRetainedAdvectionWork* work; };
    static WVComplex64 read(WVComplexInput input,std::size_t i) noexcept {
        return input.interleaved ? input.interleaved[i] : WVComplex64{input.real[i],input.imag[i]};
    }
    static void write(WVComplexOutput output,std::size_t i,WVComplex64 value) noexcept {
        if (output.interleaved) output.interleaved[i]=value;
        else { output.real[i]=value.real; output.imag[i]=value.imag; }
    }
    static WVComplex64 multiplyByI(WVComplex64 value,double factor) noexcept {
        return {-factor*value.imag,factor*value.real};
    }
    std::size_t partition(std::size_t worker) const noexcept {
        return (spec_.grid.planes/workers_)*worker+std::min(worker,spec_.grid.planes%workers_);
    }
    std::size_t targetField(std::size_t target) const noexcept {
        return target<2 ? target : target==2 ? 3 : 2;
    }
    void gather(WVComplexInput input,WVComplex64* tile,std::size_t first,
        std::size_t count) const noexcept {
        const auto& layout=spec_.retained;
        std::array<WVComplex64,4*32> block;
        for (std::size_t firstMode=0;firstMode<modes_.size();firstMode+=32) {
            const auto countModes=std::min(std::size_t{32},modes_.size()-firstMode);
            for (std::size_t mode=0;mode<countModes;++mode)
                for (std::size_t lane=0;lane<count;++lane)
                    block[mode*4+lane]=read(input,(first+lane)*layout.rowStride+
                        (firstMode+mode)*layout.columnStride);
            for (std::size_t lane=0;lane<count;++lane)
                for (std::size_t mode=0;mode<countModes;++mode)
                    tile[lane*modes_.size()+firstMode+mode]=block[mode*4+lane];
        }
    }
    void scatter(const WVComplex64* values,WVComplex64* scratch,unsigned derivative) const noexcept {
        std::fill_n(scratch,halfSize_,WVComplex64{});
        for (std::size_t mode=0;mode<modes_.size();++mode) {
            const auto& map=modes_[mode]; auto value=values[mode];
            if (derivative) value=multiplyByI(value,derivative==1 ? map.k : map.l);
            if (map.conjugate) value.imag=-value.imag;
            scratch[map.row]=value;
            if (map.partner!=map.row) scratch[map.partner]={value.real,-value.imag};
        }
    }
    void inverseY(WVComplex64* scratch,WVRetainedAdvectionCounts& counts) const noexcept {
        ++counts.columnInverses;
        fftw_execute_dft(resources_->columnInverse_.get(),reinterpret_cast<fftw_complex*>(scratch),
            reinterpret_cast<fftw_complex*>(scratch));
    }
    void inverseX(WVComplex64* scratch,double* output,
        WVRetainedAdvectionCounts& counts) const noexcept {
        ++counts.rowInverses;
        fftw_execute_dft_c2r(resources_->rowInverse_.get(),
            reinterpret_cast<fftw_complex*>(scratch),output);
        if (inverseScale_!=1)
            for (std::size_t i=0;i<spec_.grid.Nx*spec_.grid.Ny;++i) output[i]*=inverseScale_;
    }
    static void forwardTask(void* pointer,std::size_t worker) noexcept {
        auto& c=*static_cast<Context*>(pointer); auto& p=*c.plan;
        const auto& g=p.spec_.grid; const auto& layout=p.spec_.retained;
        auto* scratch=p.resources_->scratch.data()+worker*p.halfSize_;
        auto* tile=p.resources_->tile.data()+worker*tileWidth*p.modes_.size();
        const auto partition=[&](std::size_t w) { return (g.planes/p.workers_)*w+std::min(w,g.planes%p.workers_); };
        const auto begin=partition(worker),end=partition(worker+1);
        for (auto base=begin;base<end;base+=tileWidth) {
            const auto count=std::min(tileWidth,end-base);
            for (std::size_t lane=0;lane<count;++lane) {
                fftw_execute_dft_r2c(p.resources_->rowForward_.get(),const_cast<double*>(c.realInput.data+(base+lane)*g.planeStride),reinterpret_cast<fftw_complex*>(scratch));
                fftw_execute_dft(p.resources_->columnForward_.get(),reinterpret_cast<fftw_complex*>(scratch),reinterpret_cast<fftw_complex*>(scratch));
                for (std::size_t m=0;m<p.modes_.size();++m) {
                    const auto& map=p.modes_[m]; auto value=scratch[map.row];
                    value.real*=p.forwardScale_; value.imag*=map.conjugate?-p.forwardScale_:p.forwardScale_;
                    if (map.self) value.imag=0;
                    tile[lane*p.modes_.size()+m]=value;
                }
            }
            std::array<WVComplex64,tileWidth*32> block;
            for (std::size_t first=0;first<p.modes_.size();first+=32) {
                const auto modes=std::min(std::size_t{32},p.modes_.size()-first);
                for (std::size_t lane=0;lane<count;++lane) for (std::size_t m=0;m<modes;++m)
                    block[m*tileWidth+lane]=tile[lane*p.modes_.size()+first+m];
                for (std::size_t m=0;m<modes;++m) for (std::size_t lane=0;lane<count;++lane)
                    write(c.complexOutput,(base+lane)*layout.rowStride+(first+m)*layout.columnStride,block[m*tileWidth+lane]);
            }
        }
    }
    static void inverseTask(void* pointer,std::size_t worker) noexcept {
        auto& c=*static_cast<Context*>(pointer); auto& p=*c.plan;
        const auto& g=p.spec_.grid; const auto& layout=p.spec_.retained;
        auto* scratch=p.resources_->scratch.data()+worker*p.halfSize_;
        auto* tile=p.resources_->tile.data()+worker*tileWidth*p.modes_.size();
        const auto partition=[&](std::size_t w) { return (g.planes/p.workers_)*w+std::min(w,g.planes%p.workers_); };
        const auto begin=partition(worker),end=partition(worker+1);
        for (auto base=begin;base<end;base+=tileWidth) {
            const auto count=std::min(tileWidth,end-base);
            std::array<WVComplex64,tileWidth*32> block;
            for (std::size_t first=0;first<p.modes_.size();first+=32) {
                const auto modes=std::min(std::size_t{32},p.modes_.size()-first);
                for (std::size_t m=0;m<modes;++m) for (std::size_t lane=0;lane<count;++lane)
                    block[m*tileWidth+lane]=read(c.complexInput,(base+lane)*layout.rowStride+(first+m)*layout.columnStride);
                for (std::size_t lane=0;lane<count;++lane) for (std::size_t m=0;m<modes;++m)
                    tile[lane*p.modes_.size()+first+m]=block[m*tileWidth+lane];
            }
            for (std::size_t lane=0;lane<count;++lane) {
                std::fill(scratch,scratch+p.halfSize_,WVComplex64{});
                for (std::size_t m=0;m<p.modes_.size();++m) {
                    const auto& map=p.modes_[m]; auto value=tile[lane*p.modes_.size()+m];
                    if (map.conjugate) value.imag=-value.imag;
                    scratch[map.row]=value;
                    if (map.partner!=map.row) scratch[map.partner]={value.real,-value.imag};
                }
                fftw_execute_dft(p.resources_->columnInverse_.get(),reinterpret_cast<fftw_complex*>(scratch),reinterpret_cast<fftw_complex*>(scratch));
                auto* output=c.realOutput.data+(base+lane)*g.planeStride;
                fftw_execute_dft_c2r(p.resources_->rowInverse_.get(),reinterpret_cast<fftw_complex*>(scratch),output);
                if (p.inverseScale_!=1) for (std::size_t i=0;i<g.Nx*g.Ny;++i) output[i]*=p.inverseScale_;
                if (c.consumer && c.consumer->consume) {
                    const auto begin=(base+lane)*g.planeStride;
                    c.consumer->consume(c.consumer->context,begin,begin+g.Nx*g.Ny,c.realOutput.data);
                }
            }
        }
    }
    static void advectionTask(void* pointer,std::size_t worker) noexcept {
        auto& context=*static_cast<AdvectionContext*>(pointer);
        auto& plan=*context.plan; const auto& work=*context.work;
        const auto& grid=plan.spec_.grid; const auto& layout=plan.spec_.retained;
        const auto plane=grid.Nx*grid.Ny,volume=plane*grid.planes;
        const auto tileSize=plan.advectionDepth_*plan.modes_.size();
        const auto packedPerWorker=work.targets*tileSize;
        const auto stageSize=plan.activeColumns_*grid.Ny;
        auto* packed=plan.advectionPacked_.data()+worker*packedPerWorker;
        // The shared retained tile owns 16 planes per worker. Four base fields
        // require at most 4*depth planes, and prepared advection caps depth at 4.
        auto* baseTiles=plan.resources_->tile.data()+worker*tileWidth*plan.modes_.size();
        auto* stages=plan.advectionStages_.data()+worker*work.targets*stageSize;
        auto* scratch=plan.resources_->scratch.data()+worker*plan.halfSize_;
        auto* flux=plan.advectionReal_.data()+worker*2*plane;
        auto* derivative=flux+plane;
        auto& counts=plan.advectionWorkerCounts_[worker];
        const auto begin=plan.partition(worker),end=plan.partition(worker+1);
        const auto halfWidth=grid.Nx/2+1;
        for (auto first=begin;first<end;first+=plan.advectionDepth_) {
            const auto count=std::min(plan.advectionDepth_,end-first);
            for (std::size_t field=0;field<4;++field)
                plan.gather(work.base[field],baseTiles+field*tileSize,first,count);
            for (std::size_t target=0;target<work.targets;++target)
                plan.gather(work.targetSpectra[target].input(),
                    packed+target*tileSize,first,count);
            for (std::size_t lane=0;lane<count;++lane) {
                const auto z=first+lane,physicalOffset=z*plane;
                for (std::size_t field=0;field<4;++field) {
                    plan.scatter(baseTiles+field*tileSize+lane*plan.modes_.size(),scratch,0);
                    plan.inverseY(scratch,counts);
                    const auto target=field<2 ? field : field==3 ? 2 : 3;
                    if (target<work.targets)
                        for (std::size_t y=0;y<grid.Ny;++y)
                            std::copy_n(scratch+y*halfWidth,plan.activeColumns_,
                                stages+target*stageSize+y*plan.activeColumns_);
                    plan.inverseX(scratch,work.fields.data+field*volume+physicalOffset,counts);
                }
                for (std::size_t target=0;target<work.targets;++target) {
                    const auto field=plan.targetField(target);
                    std::fill_n(flux,plane,0.0);
                    for (std::size_t axis=0;axis<3;++axis) {
                        if (axis==0) {
                            std::fill_n(scratch,plan.halfSize_,WVComplex64{});
                            const auto* stage=stages+target*stageSize;
                            for (std::size_t y=0;y<grid.Ny;++y)
                                for (std::size_t x=0;x<plan.activeColumns_;++x)
                                    scratch[x+y*halfWidth]=multiplyByI(
                                        stage[x+y*plan.activeColumns_],
                                        plan.xWavenumberScale_*static_cast<double>(x));
                            ++counts.reusedColumns;
                        } else {
                            const auto* source=(axis==1 ? baseTiles+field*tileSize :
                                packed+target*tileSize)+lane*plan.modes_.size();
                            plan.scatter(source,scratch,axis==1 ? 2U : 0U);
                            plan.inverseY(scratch,counts);
                        }
                        plan.inverseX(scratch,derivative,counts);
                        const auto* velocity=work.fields.data+axis*volume+physicalOffset;
                        const auto* etaField=work.fields.data+3*volume+physicalOffset;
                        if (field==3 && axis==2)
                            native_detail::accumulateDensityAdvection(plane,flux,velocity,
                                derivative,etaField,work.densityCorrection.data[z]);
                        else
                            native_detail::accumulateAdvection(plane,flux,velocity,derivative);
                    }
                    if (work.tendencies.data) {
                        const auto outputChannel=work.targets==3 ? target : field;
                        std::copy_n(flux,plane,work.tendencies.data+
                            outputChannel*volume+physicalOffset);
                    }
                    fftw_execute_dft_r2c(plan.resources_->rowForward_.get(),flux,
                        reinterpret_cast<fftw_complex*>(scratch));
                    fftw_execute_dft(plan.resources_->columnForward_.get(),
                        reinterpret_cast<fftw_complex*>(scratch),reinterpret_cast<fftw_complex*>(scratch));
                    // Axis z was the target tile's last read for this lane.
                    auto* projected=packed+target*tileSize+lane*plan.modes_.size();
                    for (std::size_t mode=0;mode<plan.modes_.size();++mode) {
                        const auto& map=plan.modes_[mode]; auto value=scratch[map.row];
                        value.real*=plan.forwardScale_;
                        value.imag*=map.conjugate ? -plan.forwardScale_ : plan.forwardScale_;
                        if (map.self) value.imag=0;
                        projected[mode]=value;
                    }
                }
            }
            std::array<WVComplex64,4*32> block;
            for (std::size_t target=0;target<work.targets;++target)
                for (std::size_t firstMode=0;firstMode<plan.modes_.size();firstMode+=32) {
                    const auto countModes=std::min(std::size_t{32},plan.modes_.size()-firstMode);
                    const auto* projected=packed+target*tileSize;
                    for (std::size_t lane=0;lane<count;++lane)
                        for (std::size_t mode=0;mode<countModes;++mode)
                            block[mode*4+lane]=projected[lane*plan.modes_.size()+firstMode+mode];
                    for (std::size_t mode=0;mode<countModes;++mode)
                        for (std::size_t lane=0;lane<count;++lane)
                            write(work.targetSpectra[target],(first+lane)*layout.rowStride+
                                (firstMode+mode)*layout.columnStride,block[mode*4+lane]);
                }
        }
    }
    WVRetainedHorizontalSpecification spec_;
    std::size_t workers_,halfSize_,activeColumns_=0;
    std::vector<Mode> modes_;
    double forwardScale_=1,inverseScale_=1;
    double xWavenumberScale_=0,yWavenumberScale_=0;
    std::shared_ptr<RetainedFFTWResources> resources_;
    bool hasNyquist_=false;
    std::size_t advectionTargets_=0,advectionDepth_=0;
    std::vector<WVComplex64> advectionPacked_,advectionStages_;
    std::vector<double> advectionReal_;
    std::vector<WVRetainedAdvectionCounts> advectionWorkerCounts_;
};

} // namespace

struct WVFFTWRetainedCache {
    std::mutex mutex;
    std::vector<std::weak_ptr<RetainedFFTWResources>> resources;
};
WVFFTWEngine::WVFFTWEngine(std::size_t threadCount, std::string loadedLibraryPath)
    : threadCount_(threadCount), loadedLibraryPath_(std::move(loadedLibraryPath)), retainedCache_(std::make_unique<WVFFTWRetainedCache>()) {}
WVFFTWEngine::~WVFFTWEngine() = default;
std::size_t WVFFTWEngine::persistentBytes() const noexcept {
    std::lock_guard<std::mutex> lock(retainedCache_->mutex);
    return sizeof(*this)+loadedLibraryPath_.capacity()+sizeof(*retainedCache_)+
        retainedCache_->resources.capacity()*sizeof(std::weak_ptr<RetainedFFTWResources>);
}

WVKernelStatus WVFFTWEngine::create(std::size_t threadCount, std::unique_ptr<WVFFTEngine>& engine) {
    if (threadCount == 0 || threadCount > static_cast<std::size_t>(INT_MAX)) return {WVKernelStatusCode::invalidConfiguration,"FFTW thread count must lie in [1,INT_MAX]."};
    std::lock_guard<std::mutex> lock(planningMutex);
    if (fftw_init_threads() == 0) return {WVKernelStatusCode::fftPlanFailure,"FFTW thread initialization failed."};
    try {
        engine = std::unique_ptr<WVFFTEngine>(new WVFFTWEngine(threadCount,loadedLibrary()));
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"Unable to allocate the FFTW engine."};
    }
}

WVFFTWLifetimeMetrics WVFFTWEngine::lifetimeMetrics() noexcept {
    return {activePlans.load(),totalPlansCreated.load(),totalPlansDestroyed.load(),outstandingPlanningBytes.load(),1e-9*static_cast<double>(totalPlanningNanoseconds.load())};
}

WVFFTWLibraryIdentity WVFFTWEngine::linkedLibraries(const std::string& expectedOpenMPRuntime) {
    WVFFTWLibraryIdentity identity;
    identity.version = fftw_version;
    identity.baseLibrary = loadedLibrary();
    identity.threadLibrary = libraryContaining(reinterpret_cast<const void*>(&fftw_init_threads));
#if defined(__APPLE__) || defined(__linux__)
    if (!expectedOpenMPRuntime.empty()) {
        void* handle = dlopen(expectedOpenMPRuntime.c_str(),RTLD_LAZY | RTLD_NOLOAD);
        if (handle != nullptr) {
            identity.openMPRuntimeLibrary = libraryContaining(dlsym(handle,"omp_get_max_threads"));
            dlclose(handle);
        }
    }
#endif
    return identity;
}

std::string WVFFTWEngine::identifier() const { return "fftw"; }

WVKernelStatus WVFFTWEngine::createRetainedHorizontalPlan(const WVRetainedHorizontalSpecification& spec,
    std::unique_ptr<WVRetainedHorizontalPlan>& result) {
    // This hook receives a specification already validated by the operator.
    if (spec.schedule != WVRetainedHorizontalSchedule::streamingPrunedTile16 ||
        spec.grid.xStride != 1 || spec.grid.yStride != spec.grid.Nx)
        return {WVKernelStatusCode::unsupportedOperation,"Streaming retained FFTW requires contiguous horizontal planes."};
    try {
        auto plan=std::make_unique<RetainedFFTWPlan>(spec);
        std::lock_guard<std::mutex> lock(retainedCache_->mutex);
        auto status=plan->prepare(retainedCache_->resources); if (!status) return status;
        result=std::move(plan); return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Retained FFTW workspace allocation failed."}; }
      catch (const std::system_error&) { return {WVKernelStatusCode::allocationFailure,"Retained FFTW worker preparation failed."}; }
      catch (const std::overflow_error& e) { return {WVKernelStatusCode::sizeOverflow,e.what()}; }
}

WVKernelStatus WVFFTWEngine::createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) {
    try {
        const bool horizontal = specification.kind == WVFFTPlanKind::horizontalRealToComplex2D ||
                                specification.kind == WVFFTPlanKind::horizontalComplexToReal2D;
        const bool vertical = specification.kind == WVFFTPlanKind::verticalDCTI ||
                              specification.kind == WVFFTPlanKind::verticalDSTI;
        if ((!horizontal && !vertical) || specification.transformDimensions.size() != (horizontal ? 2u : 1u) ||
            specification.batchDimensions.size() > static_cast<std::size_t>(INT_MAX))
            return {WVKernelStatusCode::invalidConfiguration,"FFTW requires rank-two horizontal or rank-one vertical transforms."};
        if (horizontal && specification.inPlace)
            return {WVKernelStatusCode::unsupportedOperation,"Padded in-place horizontal FFTW storage is not supported."};
        if (specification.kind == WVFFTPlanKind::horizontalComplexToReal2D && !specification.destroysInput)
            return {WVKernelStatusCode::unsupportedOperation,"FFTW cannot preserve multidimensional c2r input; prepare destructive scratch explicitly."};
        if (specification.kind == WVFFTPlanKind::verticalDCTI && specification.transformDimensions.front().count < 2)
            return {WVKernelStatusCode::invalidConfiguration,"DCT-I requires at least two elements."};
        if (specification.inPlace) {
            for (const auto& dimension : specification.transformDimensions)
                if (dimension.inputStride != dimension.outputStride)
                    return {WVKernelStatusCode::unsupportedOperation,"In-place FFTW strides must match."};
            for (const auto& dimension : specification.batchDimensions)
                if (dimension.inputStride != dimension.outputStride)
                    return {WVKernelStatusCode::unsupportedOperation,"In-place FFTW batch strides must match."};
        }
        std::vector<fftw_iodim64> transformDimensions, batchDimensions;
        auto status = dimensions(specification.transformDimensions,transformDimensions);
        if (!status) return status;
        status = dimensions(specification.batchDimensions,batchDimensions);
        if (!status) return status;
        std::size_t inputSpan = 0, outputSpan = 0;
        status = addressedBytes(specification,true,inputSpan);
        if (!status) return status;
        status = addressedBytes(specification,false,outputSpan);
        if (!status) return status;
        std::vector<fftw_r2r_kind> kinds(vertical ? 1 : 0,
            specification.kind == WVFFTPlanKind::verticalDCTI ? FFTW_REDFT00 : FFTW_RODFT00);
        PlanningBuffer input(specification.inputBytes);
        // Avoid a second allocation for in-place transforms.
        std::unique_ptr<PlanningBuffer> output;
        if (!specification.inPlace) output = std::make_unique<PlanningBuffer>(specification.outputBytes);
        void* outputData = output ? output->get() : input.get();
        OwnedPlan rawPlan;
        const auto planningStart = std::chrono::steady_clock::now();
        {
            std::lock_guard<std::mutex> lock(planningMutex);
            fftw_plan_with_nthreads(static_cast<int>(threadCount_));
            unsigned flags = FFTW_MEASURE | FFTW_UNALIGNED;
            if (!specification.inPlace && !specification.destroysInput) flags |= FFTW_PRESERVE_INPUT;
            const int rank = static_cast<int>(transformDimensions.size());
            const int howManyRank = static_cast<int>(batchDimensions.size());
            fftw_plan created = nullptr;
            switch (specification.kind) {
                case WVFFTPlanKind::horizontalRealToComplex2D:
                    created = fftw_plan_guru64_dft_r2c(rank,transformDimensions.data(),howManyRank,batchDimensions.data(),static_cast<double*>(input.get()),static_cast<fftw_complex*>(outputData),flags);
                    break;
                case WVFFTPlanKind::horizontalComplexToReal2D:
                    created = fftw_plan_guru64_dft_c2r(rank,transformDimensions.data(),howManyRank,batchDimensions.data(),static_cast<fftw_complex*>(input.get()),static_cast<double*>(outputData),flags);
                    break;
                case WVFFTPlanKind::verticalDCTI:
                case WVFFTPlanKind::verticalDSTI:
                    created = fftw_plan_guru64_r2r(rank,transformDimensions.data(),howManyRank,batchDimensions.data(),static_cast<double*>(input.get()),static_cast<double*>(outputData),kinds.data(),flags);
                    break;
            }
            // Ownership and accounting begin immediately, before wrapper allocation.
            rawPlan.reset(created);
            if (created) { ++activePlans; ++totalPlansCreated; }
        }
        totalPlanningNanoseconds.fetch_add(static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-planningStart).count()));
        if (!rawPlan) return {WVKernelStatusCode::fftPlanFailure,"FFTW was unable to create the requested guru plan."};
        plan = std::make_unique<FFTWPlan>(std::move(rawPlan),specification.kind,specification.inPlace,inputSpan,outputSpan);
        return WVKernelStatus::ok();
    } catch (const std::bad_alloc&) {
        return {WVKernelStatusCode::allocationFailure,"FFTW allocation failed."};
    }
}

} // namespace wavevortex
