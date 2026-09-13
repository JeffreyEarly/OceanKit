#pragma once

#include "WaveVortexKernel/WVFFTEngine.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace wavevortex {
struct WVFFTWRetainedCache;

struct WVFFTWLifetimeMetrics {
    std::size_t activePlans = 0;
    std::size_t totalPlansCreated = 0;
    std::size_t totalPlansDestroyed = 0;
    std::size_t outstandingPlanningBytes = 0;
    double totalPlanningSeconds = 0.0;
};

struct WVFFTWLibraryIdentity {
    std::string version;
    std::string baseLibrary;
    std::string threadLibrary;
    std::string openMPRuntimeLibrary;
};

// Native FFTW adapter for the portable engine contract.
// The build layer supplies the pinned external FFTW provider.
class WVFFTWEngine final : public WVFFTEngine {
public:
    static WVKernelStatus create(std::size_t threadCount, std::unique_ptr<WVFFTEngine>& engine);
    ~WVFFTWEngine() override;
    static WVFFTWLifetimeMetrics lifetimeMetrics() noexcept;
    static WVFFTWLibraryIdentity linkedLibraries(const std::string& expectedOpenMPRuntime = {});

    std::string identifier() const override;
    std::string libraryIdentity() const override { return loadedLibraryPath_; }
    std::size_t persistentBytes() const noexcept override;
    WVKernelStatus createRetainedHorizontalPlan(const WVRetainedHorizontalSpecification&,
        std::unique_ptr<WVRetainedHorizontalPlan>&) override;
    WVKernelStatus createPlan(const WVFFTPlanSpecification& specification, std::unique_ptr<WVFFTPlan>& plan) override;
    std::size_t threadCount() const noexcept { return threadCount_; }
    const std::string& loadedLibraryPath() const noexcept { return loadedLibraryPath_; }

private:
    WVFFTWEngine(std::size_t threadCount, std::string loadedLibraryPath);
    std::size_t threadCount_ = 1;
    std::string loadedLibraryPath_;
    std::unique_ptr<WVFFTWRetainedCache> retainedCache_;
};

} // namespace wavevortex
