#include "WaveVortexRuntime/WVConstantStratificationIntegrationSystem.hpp"
#include "WVObserverAdapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

std::size_t observerRecordDynamicBytes(
    const WVObserverRecord &record) noexcept {
  std::size_t bytes = record.identifier.capacity() + record.name.capacity() +
                      record.typeIdentifier.capacity() +
                      record.configuration.persistentBytes() -
                          sizeof(WVPortableTypedRecord) +
                      record.stateBlockIdentifiers.capacity() *
                          sizeof(std::string) +
                      record.fieldNames.capacity() * sizeof(std::string) +
                      (record.x.capacity() + record.y.capacity() +
                       record.z.capacity()) *
                          sizeof(double);
  for (const auto &identifier : record.stateBlockIdentifiers)
    bytes += identifier.capacity();
  for (const auto &field : record.fieldNames)
    bytes += field.capacity();
  return bytes;
}

std::size_t blockIndex(const WVIntegrationStateLayout &layout,
                       const std::string &identifier) {
  const auto &blocks = layout.additionalBlocks();
  const auto found = std::find_if(
      blocks.begin(), blocks.end(), [&](const auto &candidate) {
        return candidate.identifier == identifier;
      });
  return found == blocks.end()
             ? std::numeric_limits<std::size_t>::max()
             : static_cast<std::size_t>(found - blocks.begin());
}

class UnifiedErrorPolicy final : public WVIntegrationErrorPolicy {
public:
  UnifiedErrorPolicy(std::unique_ptr<WVIntegrationErrorPolicy> coefficients,
                     const WVIntegrationStateLayout &layout)
      : coefficients_(std::move(coefficients)),
        coefficientFamilyCount_(layout.coefficientFamilyCount()) {
    counts_.reserve(coefficientFamilyCount_ + layout.additionalBlocks().size());
    tolerances_.reserve(coefficientFamilyCount_ +
                        layout.additionalBlocks().size());
    for (std::size_t component = 0; component < coefficientFamilyCount_;
         ++component) {
      counts_.push_back(coefficients_->elementCount(component));
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
    if (component < coefficientFamilyCount_)
      return coefficients_->absoluteTolerance(component, index);
    return component < tolerances_.size() ? tolerances_[component] : 0.0;
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + coefficients_->persistentBytes() +
           counts_.capacity() * sizeof(std::size_t) +
           tolerances_.capacity() * sizeof(double);
  }

private:
  std::unique_ptr<WVIntegrationErrorPolicy> coefficients_;
  std::size_t coefficientFamilyCount_ = 0;
  std::vector<std::size_t> counts_;
  std::vector<double> tolerances_;
};

} // namespace

WVKernelStatus WVConstantStratificationIntegrationSystem::create(
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVConstantStratificationIntegrationSystem> &system) {
  return createImpl(configuration, schedule, nullptr, std::move(catalog),
                    std::move(engine), system);
}

WVKernelStatus WVConstantStratificationIntegrationSystem::create(
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor &descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVConstantStratificationIntegrationSystem> &system) {
  return createImpl(configuration, schedule, &descriptor, std::move(catalog),
                    std::move(engine), system);
}

WVKernelStatus WVConstantStratificationIntegrationSystem::createImpl(
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVConstantStratificationIntegrationSystem> &system) {
  system.reset();
  if (!catalog)
    return {WVKernelStatusCode::invalidConfiguration,
            "An integration system requires an extension catalog."};
  if (descriptor != nullptr && descriptor->catalog() != catalog)
    return {WVKernelStatusCode::invalidConfiguration,
            "The observer descriptor and integration system require the same "
            "extension catalog."};
  try {
    auto candidate =
        std::unique_ptr<WVConstantStratificationIntegrationSystem>(
            new WVConstantStratificationIntegrationSystem());
    auto status = WVConstantStratificationForcingEngine::create(
        configuration, schedule, std::move(catalog), std::move(engine),
        candidate->forcing_);
    if (!status)
      return status;
    const auto coefficientShape =
        candidate->forcing_->kernel().descriptor().spectralShape();
    WVTransformStateDescription stateDescription{
        "WVTransformConstantStratification",
        {configuration.Nx, configuration.Ny, configuration.Nz},
        {{"Ap", {coefficientShape.rows, coefficientShape.columns},
          WVToleranceKind::coefficientEnergyScaled},
         {"Am", {coefficientShape.rows, coefficientShape.columns},
          WVToleranceKind::coefficientEnergyScaled},
         {"A0", {coefficientShape.rows, coefficientShape.columns},
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
    static const std::vector<WVObserverRecord> noObservers;
    const auto &observers =
        descriptor == nullptr ? noObservers : descriptor->observers();
    WVObserverIntegrationBinder integrationBinder;
    integrationBinder.advectedScalar =
        [&](const WVObserverRecord &observer) -> WVKernelStatus {
      const auto block = blockIndex(
          candidate->layout_, observer.stateBlockIdentifiers.front());
      if (block == std::numeric_limits<std::size_t>::max())
        return invalid(
            "Tracer state block is absent from the integration layout.");
      const auto &dimensions =
          candidate->layout_.additionalBlocks()[block].dimensions;
      if (observer.isXYOnly) {
        if (dimensions != std::vector<std::size_t>(
                              {configuration.Nx, configuration.Ny}))
          return {WVKernelStatusCode::invalidShape,
                  "A two-dimensional tracer must have shape [Nx,Ny]."};
        return {WVKernelStatusCode::unsupportedOperation,
                "Two-dimensional tracer integration requires a future barotropic runtime; the constant-stratification runtime supports three-dimensional tracers only."};
      }
      if (dimensions != std::vector<std::size_t>(
                            {configuration.Nx, configuration.Ny,
                             configuration.Nz}))
        return {WVKernelStatusCode::invalidShape,
                "A constant-stratification tracer must have shape [Nx,Ny,Nz]."};
      ++ownerCounts[block];
      WVTracer tracer;
      tracer.record_ = observer;
      tracer.stateBlock_ = block;
      candidate->tracers_.push_back(std::move(tracer));
      return WVKernelStatus::ok();
    };
    integrationBinder.advectedPositions =
        [&](const WVObserverRecord &observer) -> WVKernelStatus {
      if (observer.isXYOnly && observer.z.size() != observer.x.size())
        return invalid("Constant-stratification XY particles require one fixed "
                       "z coordinate per particle.");
      WVLagrangianParticles particles;
      particles.record_ = observer;
      particles.particleCount_ = observer.x.size();
      particles.positionOffset_ = positionOffset;
      particles.xBlock_ = blockIndex(candidate->layout_,
                                     observer.stateBlockIdentifiers[0]);
      particles.yBlock_ = blockIndex(candidate->layout_,
                                     observer.stateBlockIdentifiers[1]);
      particles.zBlock_ = observer.isXYOnly
                              ? std::numeric_limits<std::size_t>::max()
                              : blockIndex(candidate->layout_,
                                           observer.stateBlockIdentifiers[2]);
      if (particles.xBlock_ == std::numeric_limits<std::size_t>::max() ||
          particles.yBlock_ == std::numeric_limits<std::size_t>::max() ||
          (!observer.isXYOnly &&
           particles.zBlock_ == std::numeric_limits<std::size_t>::max()))
        return invalid("Particle state blocks are absent from the integration "
                       "layout or appear in an incompatible order.");
      ++ownerCounts[particles.xBlock_];
      ++ownerCounts[particles.yBlock_];
      if (!observer.isXYOnly)
        ++ownerCounts[particles.zBlock_];
      const auto addVelocity = [&](const char *field) {
        const auto output = velocityRequests.size();
        velocityRequests.push_back(
            {observer.identifier + '-' + field, field, positionOffset,
             particles.particleCount_, observer.advectionInterpolation});
        return output;
      };
      particles.uOutput_ = addVelocity("u");
      particles.vOutput_ = addVelocity("v");
      particles.wOutput_ = observer.isXYOnly
                               ? std::numeric_limits<std::size_t>::max()
                               : addVelocity("w");
      positionOffset += particles.particleCount_;
      candidate->particles_.push_back(std::move(particles));
      return WVKernelStatus::ok();
    };
    for (const auto &observer : observers) {
      const auto *resolved = descriptor->resolvedObserver(observer);
      if (resolved == nullptr)
        return {WVKernelStatusCode::unsupportedOperation,
                "Integration system received an unresolved observer."};
      status = resolved->implementation().bindIntegration(
          observer, integrationBinder);
      if (!status)
        return status;
    }
    for (std::size_t block = 0; block < ownerCounts.size(); ++block) {
      if (ownerCounts[block] != 1)
        return invalid("Integration state block " +
                       candidate->layout_.additionalBlocks()[block].identifier +
                       " must resolve to exactly one integrated observer.");
    }

    const bool requiresScalarAntialiasPlan = std::any_of(
        candidate->tracers_.begin(), candidate->tracers_.end(),
        [](const auto &tracer) { return tracer.shouldAntialias(); });
    if (requiresScalarAntialiasPlan) {
      status = candidate->forcing_->kernel().prepareScalarAdvection();
      if (!status)
        return status;
    }

    // Coefficient-only integration is the zero-additional-block case. It must
    // not retain the field-evaluation scratch used by observing systems.
    if (!observers.empty()) {
      status = WVFieldEvaluationService::createBorrowing(
          candidate->forcing_->kernel(), candidate->fields_);
      if (!status)
        return status;
      status = candidate->fields_->createMovingPlan(
          velocityRequests, candidate->velocityPlan_);
      if (!status)
        return status;
    }
    candidate->x_.resize(positionOffset);
    candidate->y_.resize(positionOffset);
    candidate->z_.resize(positionOffset);
    candidate->velocityStorage_.resize(velocityRequests.size());
    candidate->velocityViews_.resize(velocityRequests.size());
    for (std::size_t index = 0; index < velocityRequests.size(); ++index) {
      candidate->velocityStorage_[index].resize(
          velocityRequests[index].positionCount);
      candidate->velocityViews_[index] = {
          candidate->velocityStorage_[index].data(),
          candidate->velocityStorage_[index].size()};
    }
    candidate->metrics_.positionCapacityBytes =
        (candidate->x_.capacity() + candidate->y_.capacity() +
         candidate->z_.capacity()) *
        sizeof(double);
    for (const auto &values : candidate->velocityStorage_)
      candidate->metrics_.velocityCapacityBytes +=
          values.capacity() * sizeof(double);
    candidate->metrics_.persistentBytes = candidate->persistentBytes();
    system = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the constant-stratification particle system."};
  }
}

WVConstantStratificationIntegrationSystem::~WVConstantStratificationIntegrationSystem() =
    default;

WVKernelStatus
WVConstantStratificationIntegrationSystem::evaluateRightHandSide(
    const WVIntegrationState &state, WVIntegrationFlux &rightHandSide) {
  const auto rightHandSideStarted = std::chrono::steady_clock::now();
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "The constant-stratification integration system is not reentrant."};
  auto status = validateIntegrationState(layout_, state);
  if (!status)
    return status;
  if (rightHandSide.additionalBlockCount != layout_.additionalBlocks().size() ||
      (rightHandSide.additionalBlockCount != 0 &&
       rightHandSide.additionalBlocks == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Integration RHS storage does not match the frozen layout."};
  executing_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{executing_};
  const bool needsAdvectionContext = !particles_.empty() || !tracers_.empty();
  WVConstantStratificationRightHandSideContext context;
  auto advectionStorage = needsAdvectionContext
                              ? fields_->advectionFieldStorage()
                              : WVRealFieldBundleView{};
  const auto waveVortexFluxStarted = std::chrono::steady_clock::now();
  status = !needsAdvectionContext
               ? forcing_->nonlinearFlux(state.waveVortex,
                                         rightHandSide.waveVortex)
               : forcing_->evaluateRightHandSideWithContext(
                     state.waveVortex, rightHandSide.waveVortex,
                     advectionStorage, context);
  if (!status)
    return status;
  metrics_.waveVortexFluxSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - waveVortexFluxStarted).count();
  if (needsAdvectionContext) ++metrics_.sharedRightHandSideContextCount;
  const auto additionalStateClearStarted = std::chrono::steady_clock::now();
  for (std::size_t block = 0; block < rightHandSide.additionalBlockCount;
       ++block) {
    const auto &layout = *rightHandSide.additionalBlocks[block].layout;
    if (layout.scalarType == WVStateScalarType::real64)
      std::fill_n(rightHandSide.additionalBlocks[block].realData,
                  layout.elementCount, 0.0);
    else
      std::fill_n(rightHandSide.additionalBlocks[block].complexData,
                  layout.elementCount, WVComplex64{});
  }
  metrics_.additionalStateClearSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - additionalStateClearStarted).count();
  for (const auto &particles : particles_) {
    const auto offset = particles.positionOffset_;
    const auto count = particles.particleCount_;
    std::copy_n(state.additionalBlocks[particles.xBlock_].realData, count,
                x_.data() + offset);
    std::copy_n(state.additionalBlocks[particles.yBlock_].realData, count,
                y_.data() + offset);
    if (particles.isXYOnly())
      std::copy(particles.record_.z.begin(), particles.record_.z.end(),
                z_.data() + offset);
    else
      std::copy_n(state.additionalBlocks[particles.zBlock_].realData, count,
                  z_.data() + offset);
  }
  const auto spatial = forcing_->kernel().descriptor().spatialShape();
  const auto tracerStarted = std::chrono::steady_clock::now();
  for (const auto &tracer : tracers_) {
    const auto &input = state.additionalBlocks[tracer.stateBlock_];
    auto &output = rightHandSide.additionalBlocks[tracer.stateBlock_];
    const WVRealVolumeConstView scalar{
        input.realData, {spatial.first, spatial.second, spatial.third}};
    WVRealVolumeView tracerFlux{
        output.realData, {spatial.first, spatial.second, spatial.third}};
    status = forcing_->advectFGridScalar(
        context, scalar, tracer.shouldAntialias(), tracerFlux);
    if (!status)
      return status;
    ++metrics_.tracerEvaluationCount;
    metrics_.tracerValueWriteCount += spatial.elementCount();
    if (tracer.shouldAntialias())
      ++metrics_.antialiasedTracerEvaluationCount;
  }
  metrics_.tracerAdvectionSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - tracerStarted).count();
  const auto particleStarted = std::chrono::steady_clock::now();
  if (!particles_.empty()) {
    status = fields_->evaluateMovingFromAdvectionFields(
        velocityPlan_, state.waveVortex, context.advectionFields(),
        {x_.data(), y_.data(), z_.data(), x_.size()}, velocityViews_.data(),
        velocityViews_.size());
    if (!status)
      return status;
    ++metrics_.velocityFieldEvaluationCount;
  }
  for (const auto &particles : particles_) {
    const auto count = particles.particleCount_;
    std::copy_n(velocityStorage_[particles.uOutput_].data(), count,
                rightHandSide.additionalBlocks[particles.xBlock_].realData);
    std::copy_n(velocityStorage_[particles.vOutput_].data(), count,
                rightHandSide.additionalBlocks[particles.yBlock_].realData);
    if (!particles.isXYOnly())
      std::copy_n(velocityStorage_[particles.wOutput_].data(), count,
                  rightHandSide.additionalBlocks[particles.zBlock_].realData);
    metrics_.particleValueWriteCount +=
        count * (particles.isXYOnly() ? 2 : 3);
  }
  metrics_.particleAdvectionSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - particleStarted).count();
  ++metrics_.rightHandSideEvaluationCount;
  metrics_.rightHandSideSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - rightHandSideStarted).count();
  return WVKernelStatus::ok();
}

WVStateConstraintResult
WVConstantStratificationIntegrationSystem::enforceStateConstraints(
    WVMutableIntegrationState &state) {
  const auto status = validateMutableIntegrationState(layout_, state);
  if (!status)
    return {status, 0, false};
  return forcing_->restoreForcingAmplitudes(state.waveVortex.coefficients);
}

WVKernelStatus
WVConstantStratificationIntegrationSystem::initializeParticleState(
    WVMutableIntegrationState &state) const {
  const auto status = validateMutableIntegrationState(layout_, state);
  if (!status)
    return status;
  for (const auto &particles : particles_) {
    std::copy(particles.record_.x.begin(), particles.record_.x.end(),
              state.additionalBlocks[particles.xBlock_].realData);
    std::copy(particles.record_.y.begin(), particles.record_.y.end(),
              state.additionalBlocks[particles.yBlock_].realData);
    if (!particles.isXYOnly())
      std::copy(particles.record_.z.begin(), particles.record_.z.end(),
                state.additionalBlocks[particles.zBlock_].realData);
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVConstantStratificationIntegrationSystem::createErrorPolicy(
    double absoluteToleranceScale,
    std::unique_ptr<WVIntegrationErrorPolicy> &policy) const {
  policy.reset();
  std::unique_ptr<WVIntegrationErrorPolicy> coefficients;
  auto status = forcing_->createErrorPolicy(absoluteToleranceScale, coefficients);
  if (!status)
    return status;
  try {
    policy = std::make_unique<UnifiedErrorPolicy>(
        std::move(coefficients), layout_);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Integration error-policy allocation failed."};
  }
}

WVKernelStatus
WVConstantStratificationIntegrationSystem::evaluateFixedTimeStepCandidates(
    const WVIntegrationState &state, double cfl,
    WVFixedTimeStepCandidates &candidates) {
  if (!std::isfinite(cfl) || cfl <= 0.0)
    return invalid("CFL must be finite and positive.");
  auto status = validateIntegrationState(layout_, state);
  if (!status)
    return status;
  const auto started = std::chrono::steady_clock::now();
  const auto &descriptor = forcing_->kernel().descriptor();
  const auto &configuration = descriptor.configuration();
  const auto spatial = descriptor.spatialShape();
  const auto fieldElements = spatial.elementCount();
  try {
    const auto velocityElements = 3 * fieldElements;
    std::unique_ptr<double[]> velocity(new double[velocityElements]);
    WVRealFieldBundleView fields{
        velocity.get(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 3}};
    status = forcing_->kernel().transformWaveVortexToUVW(state.waveVortex,
                                                         fields);
    if (!status)
      return status;

    WVFixedTimeStepCandidates result;
    result.transientWorkspaceMaximumLiveBytes =
        velocityElements * sizeof(double);
    double maximumHorizontalWavenumber = 0.0;
    for (const auto &mode : descriptor.fourierModes())
      maximumHorizontalWavenumber =
          std::max({maximumHorizontalWavenumber, std::abs(mode.k),
                    std::abs(mode.l)});
    if (!(maximumHorizontalWavenumber > 0.0) ||
        !std::isfinite(maximumHorizontalWavenumber))
      return invalid("The effective horizontal resolution is unavailable.");
    result.effectiveHorizontalGridResolution =
        std::acos(-1.0) / maximumHorizontalWavenumber;

    const double *u = velocity.get();
    const double *v = u + fieldElements;
    const double *w = v + fieldElements;
    for (std::size_t index = 0; index < fieldElements; ++index) {
      if (!std::isfinite(u[index]) || !std::isfinite(v[index]) ||
          !std::isfinite(w[index]))
        return invalid("CFL velocity reconstruction produced a nonfinite "
                       "value.");
      result.maximumHorizontalSpeed =
          std::max(result.maximumHorizontalSpeed,
                   std::sqrt(u[index] * u[index] + v[index] * v[index]));
    }
    if (result.maximumHorizontalSpeed > 0.0)
      result.horizontalAdvective =
          cfl * result.effectiveHorizontalGridResolution /
          result.maximumHorizontalSpeed;

    const auto horizontalElements = configuration.Nx * configuration.Ny;
    const auto &z = descriptor.verticalModes().z;
    double maximumVerticalRate = 0.0;
    for (std::size_t iz = 1; iz < configuration.Nz; ++iz) {
      const double denominator = 1.5 * (z[iz] - z[iz - 1]);
      if (!(denominator > 0.0) || !std::isfinite(denominator))
        return invalid("The vertical CFL spacing is invalid.");
      for (std::size_t horizontal = 0;
           horizontal < horizontalElements; ++horizontal) {
        const auto index = horizontal + horizontalElements * iz;
        maximumVerticalRate =
            std::max(maximumVerticalRate,
                     std::abs(w[index] / denominator));
      }
    }
    if (maximumVerticalRate > 0.0)
      result.verticalAdvective = cfl / maximumVerticalRate;
    result.advective =
        std::min(result.horizontalAdvective, result.verticalAdvective);

    const auto &omega = descriptor.verticalModes().omega;
    for (std::size_t mode = 0; mode < descriptor.Nkl(); ++mode)
      for (std::size_t j = 1; j < configuration.Nj; ++j)
        result.highestActiveWaveFrequency = std::max(
            result.highestActiveWaveFrequency,
            std::abs(omega[j + configuration.Nj * mode]));
    if (result.highestActiveWaveFrequency > 0.0)
      result.oscillatory =
          cfl * 2.0 * std::acos(-1.0) /
          result.highestActiveWaveFrequency;
    result.evaluationSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    candidates = result;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate transient CFL velocity workspace."};
  }
}

std::size_t
WVConstantStratificationIntegrationSystem::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) + layout_.persistentBytes() +
                      forcing_->persistentBytes() +
                      (fields_ ? fields_->persistentBytes() : 0) +
                      velocityPlan_.persistentBytes() - sizeof(velocityPlan_) +
                      particles_.capacity() * sizeof(WVLagrangianParticles) +
                      tracers_.capacity() * sizeof(WVTracer) +
                      velocityStorage_.capacity() *
                          sizeof(std::vector<double>) +
                      metrics_.positionCapacityBytes +
                      metrics_.velocityCapacityBytes +
                      velocityViews_.capacity() * sizeof(WVFieldOutputView);
  for (const auto &particles : particles_)
    bytes += observerRecordDynamicBytes(particles.record_);
  for (const auto &tracer : tracers_)
    bytes += observerRecordDynamicBytes(tracer.record_);
  return bytes;
}

} // namespace wavevortex::runtime
