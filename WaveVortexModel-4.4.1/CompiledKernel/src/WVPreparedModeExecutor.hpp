#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <mutex>
#include <memory>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

namespace wavevortex::kernel_detail {

// One synchronous caller owns the workspace. The operation is borrowed until
// every worker finishes; dispatch requires neither std::function nor allocation.
class WVPreparedModeExecutor final {
public:
    struct LaunchThread {
        template <typename Work>
        std::thread operator()(Work work) const { return std::thread(work); }
    };

    explicit WVPreparedModeExecutor(std::size_t count)
        : WVPreparedModeExecutor(count, LaunchThread{}) {}

    // The launch seam lets tests fail after an actual worker has started.
    template <typename Launch>
    WVPreparedModeExecutor(std::size_t count, Launch launch) {
        if (count == 0) throw std::invalid_argument("A mode executor requires a worker.");
        workers_.reserve(count - 1);
        try {
            for (std::size_t i = 1; i < count; ++i)
                workers_.emplace_back(launch([this, i] { worker(i); }));
        } catch (...) {
            stop();
            throw;
        }
    }

    ~WVPreparedModeExecutor() { stop(); }
    WVPreparedModeExecutor(const WVPreparedModeExecutor&) = delete;
    WVPreparedModeExecutor& operator=(const WVPreparedModeExecutor&) = delete;

    std::size_t workerCount() const noexcept { return workers_.size() + 1; }
    std::size_t persistentBytes() const noexcept {
        // Thread-runtime heap storage, stacks and OS storage are not observable.
        return sizeof(*this) + workers_.capacity() * sizeof(std::thread);
    }

    template <typename Operation>
    void execute(std::size_t count, Operation&& operation) {
        if (workers_.empty()) {
            operation(std::size_t{0}, count);
            return;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            count_ = count;
            blockSize_ = count / workerCount() + (count % workerCount() != 0);
            operation_ = const_cast<void*>(static_cast<const void*>(std::addressof(operation)));
            invoke_ = [](void* value, std::size_t begin, std::size_t end) {
                (*static_cast<std::remove_reference_t<Operation>*>(value))(begin, end);
            };
            failure_ = nullptr;
            remaining_ = workers_.size();
            generation_ = !generation_;
        }
        ready_.notify_all();
        try {
            operation(std::size_t{0}, std::min(blockSize_, count));
        } catch (...) {
            recordFailure();
        }
        std::unique_lock<std::mutex> lock(mutex_);
        finished_.wait(lock, [this] { return remaining_ == 0; });
        operation_ = nullptr;
        invoke_ = nullptr;
        auto failure = failure_;
        failure_ = nullptr;
        lock.unlock();
        if (failure) std::rethrow_exception(failure);
    }

private:
    void recordFailure() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!failure_) failure_ = std::current_exception();
    }

    void worker(std::size_t index) {
        bool observedGeneration = false;
        std::unique_lock<std::mutex> lock(mutex_);
        while (true) {
            ready_.wait(lock, [&] { return stopping_ || generation_ != observedGeneration; });
            if (stopping_) return;
            observedGeneration = generation_;
            const auto begin = blockSize_ == 0 || index > count_ / blockSize_
                ? count_ : index * blockSize_;
            const auto end = begin + std::min(blockSize_, count_ - begin);
            auto invoke = invoke_;
            auto operation = operation_;
            lock.unlock();
            try {
                invoke(operation, begin, end);
            } catch (...) {
                recordFailure();
            }
            lock.lock();
            if (--remaining_ == 0) finished_.notify_one();
        }
    }

    void stop() noexcept {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (auto& worker : workers_) if (worker.joinable()) worker.join();
    }

    std::mutex mutex_;
    std::condition_variable ready_, finished_;
    std::vector<std::thread> workers_;
    void* operation_ = nullptr;
    void (*invoke_)(void*, std::size_t, std::size_t) = nullptr;
    std::exception_ptr failure_;
    std::size_t count_ = 0, blockSize_ = 0, remaining_ = 0;
    bool generation_ = false, stopping_ = false;
};

} // namespace wavevortex::kernel_detail
