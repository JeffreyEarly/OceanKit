#include "WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVObserverAdapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <utility>
#include <vector>

namespace wavevortex::runtime {
namespace {

template <typename T>
std::size_t vectorBytes(const std::vector<T> &values) noexcept {
  return values.capacity() * sizeof(T);
}

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

bool validCoordinate(const std::vector<double> &coordinate,
                     double length) noexcept {
  if (coordinate.size() < 2 || !std::isfinite(length) || length <= 0.0)
    return false;
  const double spacing = length / static_cast<double>(coordinate.size());
  const double scale = std::max(1.0, std::abs(length));
  const double tolerance =
      64.0 * std::numeric_limits<double>::epsilon() * scale;
  for (std::size_t index = 0; index < coordinate.size(); ++index) {
    if (!std::isfinite(coordinate[index]) ||
        (index > 0 && !(coordinate[index] > coordinate[index - 1])) ||
        std::abs(coordinate[index] - spacing * static_cast<double>(index)) >
            tolerance)
      return false;
  }
  return true;
}

std::size_t blockIndex(const WVIntegrationStateLayout &layout,
                       const std::string &identifier) noexcept {
  const auto &blocks = layout.additionalBlocks();
  const auto found = std::find_if(
      blocks.begin(), blocks.end(), [&](const auto &block) {
        return block.identifier == identifier;
      });
  return found == blocks.end()
             ? std::numeric_limits<std::size_t>::max()
             : static_cast<std::size_t>(found - blocks.begin());
}

class QGUnifiedErrorPolicy final : public WVIntegrationErrorPolicy {
public:
  QGUnifiedErrorPolicy(std::unique_ptr<WVIntegrationErrorPolicy> coefficients,
                       const WVIntegrationStateLayout &layout)
      : coefficients_(std::move(coefficients)),
        coefficientCount_(layout.coefficientFamilyCount()) {
    counts_.reserve(coefficientCount_ + layout.additionalBlocks().size());
    tolerances_.reserve(counts_.capacity());
    for (std::size_t family = 0; family < coefficientCount_; ++family) {
      counts_.push_back(coefficients_->elementCount(family));
      tolerances_.push_back(0.0);
    }
    for (const auto &block : layout.additionalBlocks()) {
      counts_.push_back(block.elementCount);
      tolerances_.push_back(block.absoluteTolerance);
    }
  }
  std::size_t componentCount() const noexcept override {
    return counts_.size();
  }
  std::size_t elementCount(std::size_t component) const noexcept override {
    return component < counts_.size() ? counts_[component] : 0;
  }
  double absoluteTolerance(std::size_t component,
                           std::size_t index) const noexcept override {
    return component < coefficientCount_
               ? coefficients_->absoluteTolerance(component, index)
               : component < tolerances_.size() ? tolerances_[component]
                                                 : 0.0;
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + coefficients_->persistentBytes() +
           vectorBytes(counts_) + vectorBytes(tolerances_);
  }

private:
  std::unique_ptr<WVIntegrationErrorPolicy> coefficients_;
  std::size_t coefficientCount_ = 0;
  std::vector<std::size_t> counts_;
  std::vector<double> tolerances_;
};

class BarotropicQGErrorPolicy final : public WVIntegrationErrorPolicy {
public:
  static WVKernelStatus create(
      const WVTransformBarotropicQGDescriptor &descriptor,
      double absoluteToleranceScale,
      std::unique_ptr<WVIntegrationErrorPolicy> &result) {
    result.reset();
    if (!std::isfinite(absoluteToleranceScale) ||
        absoluteToleranceScale <= 0.0)
      return invalid(
          "Adaptive absolute-tolerance scale must be finite and positive.");
    try {
      auto policy = std::unique_ptr<BarotropicQGErrorPolicy>(
          new BarotropicQGErrorPolicy());
      const auto &horizontal = descriptor.fourierModes();
      policy->tolerances_.assign(horizontal.size(), 1.0);
      std::vector<double> uniqueKh;
      uniqueKh.reserve(horizontal.size());
      for (const auto &mode : horizontal)
        uniqueKh.push_back(std::abs(mode.Kh));
      std::sort(uniqueKh.begin(), uniqueKh.end());
      uniqueKh.erase(std::unique(uniqueKh.begin(), uniqueKh.end()),
                     uniqueKh.end());
      if (uniqueKh.size() < 2)
        return invalid("Adaptive Barotropic QG tolerances require at least "
                       "two horizontal radial wavenumbers.");
      double deltaK = 0.0;
      for (std::size_t index = 1; index < uniqueKh.size(); ++index)
        deltaK = std::max(deltaK, uniqueKh[index] - uniqueKh[index - 1]);
      if (!(deltaK > 0.0) || !std::isfinite(deltaK))
        return invalid("Adaptive Barotropic QG radial spacing is invalid.");
      std::vector<double> centers;
      for (double center = 0.0;
           center <= uniqueKh.back() + 0.5 * deltaK; center += deltaK)
        centers.push_back(center);
      std::vector<std::size_t> bins(horizontal.size(), centers.size());
      std::vector<std::size_t> counts(centers.size(), 0);
      for (std::size_t mode = 0; mode < horizontal.size(); ++mode) {
        for (std::size_t bin = 0; bin < centers.size(); ++bin) {
          if (centers[bin] - 0.5 * deltaK < horizontal[mode].Kh &&
              horizontal[mode].Kh <= centers[bin] + 0.5 * deltaK) {
            bins[mode] = bin;
            ++counts[bin];
            break;
          }
        }
        if (bins[mode] == centers.size())
          return invalid("A Barotropic QG mode was not assigned to an "
                         "adaptive radial bin.");
      }
      const auto &energy = descriptor.modes().energyFactor;
      for (std::size_t mode = 0; mode < horizontal.size(); ++mode) {
        if (!(energy[mode] > 0.0) || !std::isfinite(energy[mode]))
          continue;
        const auto bin = bins[mode];
        const double radialEnergy =
            centers[bin] + 0.5 * deltaK -
            std::max(centers[bin] - 0.5 * deltaK, 0.0);
        const double energyPerCoefficient =
            radialEnergy / static_cast<double>(counts[bin]);
        policy->tolerances_[mode] =
            absoluteToleranceScale *
            std::sqrt(energyPerCoefficient / energy[mode]);
      }
      result = std::move(policy);
      return WVKernelStatus::ok();
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate Barotropic QG adaptive tolerances."};
    }
  }

  std::size_t componentCount() const noexcept override { return 1; }
  std::size_t elementCount(std::size_t component) const noexcept override {
    return component == 0 ? tolerances_.size() : 0;
  }
  double absoluteTolerance(std::size_t component,
                           std::size_t index) const noexcept override {
    return component == 0 && index < tolerances_.size()
               ? tolerances_[index]
               : std::numeric_limits<double>::quiet_NaN();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + vectorBytes(tolerances_);
  }

private:
  std::vector<double> tolerances_;
};

} // namespace

std::size_t
WVBarotropicQGNumericalConfiguration::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      stateDescription.transformIdentifier.capacity() +
                      vectorBytes(stateDescription.spatialDimensions) +
                      stateDescription.coefficientFamilies.capacity() *
                          sizeof(WVCoefficientFamilyDescription);
  for (const auto &family : stateDescription.coefficientFamilies)
    bytes += family.identifier.capacity() +
             vectorBytes(family.spectralDimensions);
  return bytes;
}

WVKernelStatus decodeBarotropicQGNumericalConfiguration(
    const WVBarotropicQGPersistedNumericalRecord &record,
    WVBarotropicQGNumericalConfiguration &configuration) {
  if (!validCoordinate(record.x, record.Lx) ||
      !validCoordinate(record.y, record.Ly))
    return invalid("Persisted Barotropic QG x and y coordinates must be "
                   "finite, strictly increasing, uniformly periodic axes.");
  if (!std::isfinite(record.t) || !std::isfinite(record.t0))
    return invalid("Persisted Barotropic QG times must be finite.");
  WVTransformBarotropicQGConfiguration transform;
  transform.Nx = record.x.size();
  transform.Ny = record.y.size();
  transform.Lx = record.Lx;
  transform.Ly = record.Ly;
  transform.h = record.h;
  transform.j = record.j;
  transform.g = record.g;
  transform.planetaryRadius = record.planetaryRadius;
  transform.rotationRate = record.rotationRate;
  transform.latitude = record.latitude;
  transform.shouldAntialias = record.shouldAntialias;
  WVTransformBarotropicQGDescriptor descriptor;
  auto status = WVTransformBarotropicQGDescriptor::create(transform,
                                                           descriptor);
  if (!status)
    return status;
  try {
    WVBarotropicQGNumericalConfiguration candidate;
    candidate.transform = transform;
    candidate.t = record.t;
    candidate.t0 = record.t0;
    candidate.stateDescription = {
        "WVTransformBarotropicQG", {transform.Nx, transform.Ny},
        {{"A0", {descriptor.Nkl()},
          WVToleranceKind::coefficientEnergyScaled}},
        true};
    configuration = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to retain decoded Barotropic QG configuration."};
  }
}

WVKernelStatus WVBarotropicQGIntegrationSystem::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVBarotropicQGIntegrationSystem> &system) {
  std::shared_ptr<const WVExtensionCatalog> catalog;
  auto status = makeBuiltInExtensionCatalog(catalog);
  if (!status)
    return status;
  return create(configuration, defaultNonlinearAdvectionSchedule(),
                std::move(catalog), std::move(engine), system);
}

WVKernelStatus WVBarotropicQGIntegrationSystem::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVBarotropicQGIntegrationSystem> &system) {
  return createImpl(configuration, schedule, nullptr, std::move(catalog),
                    std::move(engine), system);
}

WVKernelStatus WVBarotropicQGIntegrationSystem::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor &descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVBarotropicQGIntegrationSystem> &system) {
  return createImpl(configuration, schedule, &descriptor, std::move(catalog),
                    std::move(engine), system);
}

WVKernelStatus WVBarotropicQGIntegrationSystem::createImpl(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVBarotropicQGIntegrationSystem> &system) {
  system.reset();
  if (!catalog)
    return invalid("A Barotropic QG integration system requires an extension "
                   "catalog.");
  if (descriptor != nullptr && descriptor->catalog() != catalog)
    return invalid("The Barotropic QG observer descriptor and numerical "
                   "system require the same extension catalog.");
  try {
    auto candidate = std::unique_ptr<WVBarotropicQGIntegrationSystem>(
        new WVBarotropicQGIntegrationSystem());
    auto status = WVBarotropicQGForcingEngine::create(
        configuration, schedule, catalog, std::move(engine),
        candidate->forcingEngine_);
    if (!status)
      return status;
    WVTransformStateDescription stateDescription{
        "WVTransformBarotropicQG", {configuration.Nx, configuration.Ny},
        {{"A0", {candidate->kernel().descriptor().Nkl()},
          WVToleranceKind::coefficientEnergyScaled}},
        true};
    status = descriptor == nullptr
                 ? WVIntegrationStateLayout::createCoefficientOnly(
                       std::move(stateDescription), candidate->layout_)
                 : WVIntegrationStateLayout::create(
                       std::move(stateDescription), *descriptor,
                       candidate->layout_);
    if (!status)
      return status;

    std::vector<std::size_t> ownerCounts(
        candidate->layout_.additionalBlocks().size(), 0);
    std::vector<WVMovingFieldRequest> velocityRequests;
    std::size_t positionOffset = 0;
    if (descriptor != nullptr) {
      WVObserverIntegrationBinder binder;
      binder.advectedPositions =
          [&](const WVObserverRecord &observer) -> WVKernelStatus {
        if (!observer.isXYOnly)
          return {WVKernelStatusCode::unsupportedOperation,
                  "WVTransformBarotropicQG supports XY particles only."};
        if (!observer.z.empty())
          return {WVKernelStatusCode::unsupportedOperation,
                  "Barotropic QG XY particles do not persist a vertical "
                  "coordinate."};
        if (observer.stateBlockIdentifiers.size() != 2)
          return invalid("Barotropic QG particles require x and y state "
                         "blocks only.");
        Particle particle;
        particle.record = observer;
        particle.particleCount = observer.x.size();
        particle.positionOffset = positionOffset;
        particle.xBlock = blockIndex(
            candidate->layout_, observer.stateBlockIdentifiers[0]);
        particle.yBlock = blockIndex(
            candidate->layout_, observer.stateBlockIdentifiers[1]);
        if (particle.xBlock == std::numeric_limits<std::size_t>::max() ||
            particle.yBlock == std::numeric_limits<std::size_t>::max())
          return invalid("Barotropic QG particle state blocks are absent from "
                         "the integration layout.");
        ++ownerCounts[particle.xBlock];
        ++ownerCounts[particle.yBlock];
        particle.uOutput = velocityRequests.size();
        velocityRequests.push_back(
            {observer.identifier + "-u", "u", positionOffset,
             particle.particleCount, observer.advectionInterpolation});
        particle.vOutput = velocityRequests.size();
        velocityRequests.push_back(
            {observer.identifier + "-v", "v", positionOffset,
             particle.particleCount, observer.advectionInterpolation});
        positionOffset += particle.particleCount;
        candidate->particles_.push_back(std::move(particle));
        return WVKernelStatus::ok();
      };
      binder.advectedScalar =
          [&](const WVObserverRecord &observer) -> WVKernelStatus {
        if (!observer.isXYOnly)
          return {WVKernelStatusCode::unsupportedOperation,
                  "WVTransformBarotropicQG supports rank-two tracers only."};
        if (observer.stateBlockIdentifiers.size() != 1)
          return invalid("Barotropic QG tracer state is not uniquely "
                         "identified.");
        const auto block = blockIndex(
            candidate->layout_, observer.stateBlockIdentifiers.front());
        if (block == std::numeric_limits<std::size_t>::max())
          return invalid("Barotropic QG tracer state is absent from the "
                         "integration layout.");
        if (candidate->layout_.additionalBlocks()[block].dimensions !=
            std::vector<std::size_t>{configuration.Nx, configuration.Ny})
          return {WVKernelStatusCode::invalidShape,
                  "A Barotropic QG tracer must have shape [Nx,Ny]."};
        ++ownerCounts[block];
        candidate->tracers_.push_back({observer, block});
        return WVKernelStatus::ok();
      };
      for (const auto &observer : descriptor->observers()) {
        const auto *resolved = descriptor->resolvedObserver(observer);
        if (resolved == nullptr)
          return {WVKernelStatusCode::unsupportedOperation,
                  "Barotropic QG integration received an unresolved observer."};
        status = resolved->implementation().bindIntegration(observer, binder);
        if (!status)
          return status;
      }
      for (std::size_t block = 0; block < ownerCounts.size(); ++block)
        if (ownerCounts[block] != 1)
          return invalid("Integrated Barotropic QG state block " +
                         candidate->layout_.additionalBlocks()[block]
                             .identifier +
                         " must resolve exactly once.");
      status = WVFieldEvaluationService::createBorrowing(
          candidate->kernel(), candidate->fields_);
      if (!status)
        return status;
      status = candidate->fields_->createMovingPlan(
          velocityRequests, candidate->velocityPlan_);
      if (!status)
        return status;
    }
    candidate->x_.resize(positionOffset);
    candidate->y_.resize(positionOffset);
    candidate->velocityStorage_.resize(velocityRequests.size());
    candidate->velocityViews_.resize(velocityRequests.size());
    for (std::size_t index = 0; index < velocityRequests.size(); ++index) {
      candidate->velocityStorage_[index].resize(
          velocityRequests[index].positionCount);
      candidate->velocityViews_[index] = {
          candidate->velocityStorage_[index].data(),
          candidate->velocityStorage_[index].size()};
    }
    candidate->observerMetrics_.positionCapacityBytes =
        (candidate->x_.capacity() + candidate->y_.capacity()) *
        sizeof(double);
    for (const auto &values : candidate->velocityStorage_)
      candidate->observerMetrics_.velocityCapacityBytes +=
          values.capacity() * sizeof(double);
    candidate->observerMetrics_.persistentBytes =
        candidate->persistentBytes();
    system = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Barotropic QG integration system."};
  }
}

WVKernelStatus WVBarotropicQGIntegrationSystem::evaluateRightHandSide(
    const WVIntegrationState &state,
    WVIntegrationFlux &rightHandSide) {
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "The Barotropic QG integration system is not reentrant."};
  auto status = validateIntegrationState(layout_, state);
  if (!status)
    return status;
  if (rightHandSide.coefficientFamilyCount != 1 ||
      rightHandSide.coefficientFamilies == nullptr ||
      rightHandSide.coefficientFamilies[0].layout !=
          &layout_.coefficientFamilies()[0] ||
      rightHandSide.additionalBlockCount != layout_.additionalBlocks().size() ||
      (rightHandSide.additionalBlockCount != 0 &&
       rightHandSide.additionalBlocks == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG RHS storage does not match compact F0 and its "
            "resolved observer blocks."};
  executing_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{executing_};
  const auto A0 = coefficientFamilyView(layout_, state, 0);
  auto F0 = coefficientFamilyView(layout_, rightHandSide, 0);
  const WVComplexConstView A0View{
      A0.data, {1, A0.layout->elementCount}};
  WVComplexView F0View{F0.data, {1, F0.layout->elementCount}};
  const bool needsAdvectionFields = !tracers_.empty() || !particles_.empty();
  WVRealFieldBundleConstView advectionFields;
  const auto forcingStarted = std::chrono::steady_clock::now();
  status = forcingEngine_->evaluateRightHandSide(
      A0View, F0View,
      needsAdvectionFields ? &advectionFields : nullptr);
  if (!status)
    return status;
  observerMetrics_.waveVortexFluxSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - forcingStarted).count();
  const auto clearStarted = std::chrono::steady_clock::now();
  for (std::size_t block = 0; block < rightHandSide.additionalBlockCount;
       ++block) {
    const auto &metadata = *rightHandSide.additionalBlocks[block].layout;
    if (metadata.scalarType == WVStateScalarType::real64)
      std::fill_n(rightHandSide.additionalBlocks[block].realData,
                  metadata.elementCount, 0.0);
    else
      std::fill_n(rightHandSide.additionalBlocks[block].complexData,
                  metadata.elementCount, WVComplex64{});
  }
  observerMetrics_.additionalStateClearSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    clearStarted).count();

  const auto tracerStarted = std::chrono::steady_clock::now();
  const auto spatial = kernel().descriptor().spatialShape();
  for (const auto &tracer : tracers_) {
    const auto &input = state.additionalBlocks[tracer.stateBlock];
    auto &output = rightHandSide.additionalBlocks[tracer.stateBlock];
    const WVRealConstView scalar{input.realData, spatial};
    WVRealView scalarFlux{output.realData, spatial};
    status = kernel().advectScalarWithAdvectionFields(
        scalar, advectionFields, tracer.record.shouldAntialias, scalarFlux);
    if (!status)
      return status;
    ++observerMetrics_.tracerEvaluationCount;
    observerMetrics_.tracerValueWriteCount += spatial.elementCount();
    if (tracer.record.shouldAntialias)
      ++observerMetrics_.antialiasedTracerEvaluationCount;
  }
  observerMetrics_.tracerAdvectionSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - tracerStarted).count();

  const auto particleStarted = std::chrono::steady_clock::now();
  for (const auto &particle : particles_) {
    std::copy_n(state.additionalBlocks[particle.xBlock].realData,
                particle.particleCount,
                x_.data() + particle.positionOffset);
    std::copy_n(state.additionalBlocks[particle.yBlock].realData,
                particle.particleCount,
                y_.data() + particle.positionOffset);
  }
  if (!particles_.empty()) {
    status = fields_->evaluateMovingFromAdvectionFields(
        velocityPlan_, state, advectionFields,
        {x_.data(), y_.data(), nullptr, x_.size()}, velocityViews_.data(),
        velocityViews_.size());
    if (!status)
      return status;
    ++observerMetrics_.velocityFieldEvaluationCount;
  }
  for (const auto &particle : particles_) {
    std::copy_n(velocityStorage_[particle.uOutput].data(),
                particle.particleCount,
                rightHandSide.additionalBlocks[particle.xBlock].realData);
    std::copy_n(velocityStorage_[particle.vOutput].data(),
                particle.particleCount,
                rightHandSide.additionalBlocks[particle.yBlock].realData);
    observerMetrics_.particleValueWriteCount += 2 * particle.particleCount;
  }
  observerMetrics_.particleAdvectionSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - particleStarted).count();
  ++observerMetrics_.rightHandSideEvaluationCount;
  if (needsAdvectionFields)
    ++observerMetrics_.sharedRightHandSideContextCount;
  return WVKernelStatus::ok();
}

WVStateConstraintResult
WVBarotropicQGIntegrationSystem::enforceStateConstraints(
    WVMutableIntegrationState &state) {
  const auto status = validateMutableIntegrationState(layout_, state);
  if (!status)
    return {status, 0, false};
  auto A0 = coefficientFamilyView(layout_, state, 0);
  WVComplexView A0View{A0.data, {1, A0.layout->elementCount}};
  auto result = forcingEngine_->restoreForcingAmplitudes(A0View);
  if (!result)
    return result;
  const auto realityModified = kernel().enforceReality(A0View);
  result.modifiedCoefficientCount += realityModified;
  result.fsalCompatible = result.fsalCompatible && realityModified == 0;
  return result;
}

WVKernelStatus WVBarotropicQGIntegrationSystem::createErrorPolicy(
    double absoluteToleranceScale,
    std::unique_ptr<WVIntegrationErrorPolicy> &policy) const {
  std::unique_ptr<WVIntegrationErrorPolicy> coefficients;
  auto status = BarotropicQGErrorPolicy::create(
      kernel().descriptor(), absoluteToleranceScale, coefficients);
  if (!status)
    return status;
  if (layout_.additionalBlocks().empty()) {
    policy = std::move(coefficients);
    return WVKernelStatus::ok();
  }
  try {
    policy = std::make_unique<QGUnifiedErrorPolicy>(
        std::move(coefficients), layout_);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Barotropic QG observer error policy."};
  }
}

WVKernelStatus WVBarotropicQGIntegrationSystem::initializeParticleState(
    WVMutableIntegrationState &state) const {
  const auto status = validateMutableIntegrationState(layout_, state);
  if (!status)
    return status;
  for (const auto &particle : particles_) {
    std::copy(particle.record.x.begin(), particle.record.x.end(),
              state.additionalBlocks[particle.xBlock].realData);
    std::copy(particle.record.y.begin(), particle.record.y.end(),
              state.additionalBlocks[particle.yBlock].realData);
  }
  return WVKernelStatus::ok();
}

WVKernelStatus
WVBarotropicQGIntegrationSystem::evaluateFixedTimeStepCandidates(
    const WVIntegrationState &state, double cfl,
    WVFixedTimeStepCandidates &candidates) {
  if (!std::isfinite(cfl) || cfl <= 0.0)
    return invalid("CFL must be finite and positive.");
  auto status = validateIntegrationState(layout_, state);
  if (!status)
    return status;
  const auto started = std::chrono::steady_clock::now();
  const auto A0 = coefficientFamilyView(layout_, state, 0);
  WVComplexConstView A0View{A0.data, {1, A0.layout->elementCount}};
  WVFixedTimeStepCandidates result;
  status = kernel().uvMax(A0View, result.maximumHorizontalSpeed);
  if (!status)
    return status;
  if (!std::isfinite(result.maximumHorizontalSpeed) ||
      result.maximumHorizontalSpeed < 0.0)
    return invalid("Barotropic QG CFL velocity is nonfinite or negative.");
  double maximumHorizontalWavenumber = 0.0;
  for (const auto &mode : kernel().descriptor().fourierModes())
    maximumHorizontalWavenumber =
        std::max({maximumHorizontalWavenumber, std::abs(mode.k),
                  std::abs(mode.l)});
  if (!(maximumHorizontalWavenumber > 0.0) ||
      !std::isfinite(maximumHorizontalWavenumber))
    return invalid("The effective horizontal resolution is unavailable.");
  result.effectiveHorizontalGridResolution =
      std::acos(-1.0) / maximumHorizontalWavenumber;
  if (result.maximumHorizontalSpeed > 0.0)
    result.horizontalAdvective =
        cfl * result.effectiveHorizontalGridResolution /
        result.maximumHorizontalSpeed;
  result.advective = result.horizontalAdvective;
  result.evaluationSeconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  candidates = result;
  return WVKernelStatus::ok();
}

std::size_t
WVBarotropicQGIntegrationSystem::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) + layout_.persistentBytes() +
                      (forcingEngine_ == nullptr
                           ? 0
                           : forcingEngine_->persistentBytes()) +
                      (fields_ == nullptr ? 0 : fields_->persistentBytes()) +
                      velocityPlan_.persistentBytes() - sizeof(velocityPlan_) +
                      particles_.capacity() * sizeof(Particle) +
                      tracers_.capacity() * sizeof(Tracer) +
                      velocityStorage_.capacity() *
                          sizeof(std::vector<double>) +
                      velocityViews_.capacity() * sizeof(WVFieldOutputView) +
                      (x_.capacity() + y_.capacity()) * sizeof(double);
  for (const auto &values : velocityStorage_)
    bytes += values.capacity() * sizeof(double);
  return bytes;
}

} // namespace wavevortex::runtime
