#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <set>
#include <sstream>
#include <utility>

namespace wavevortex::runtime {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

template <typename T>
std::size_t vectorBytes(const std::vector<T> &values) noexcept {
  return values.capacity() * sizeof(T);
}

std::size_t stageRank(WVForcingStage stage) noexcept {
  return static_cast<std::size_t>(stage);
}

const std::vector<double> *realValues(const WVPortableTypedRecord &record,
                                      const char *name) {
  const auto *value = record.value(name);
  return value == nullptr
             ? nullptr
             : std::get_if<std::vector<double>>(&value->storage);
}

const std::vector<std::int64_t> *
integerValues(const WVPortableTypedRecord &record, const char *name) {
  const auto *value = record.value(name);
  return value == nullptr
             ? nullptr
             : std::get_if<std::vector<std::int64_t>>(&value->storage);
}

bool emptyConfiguration(const WVFrozenForcingEntry &entry) {
  return entry.configuration.values.empty();
}

double vanishingFilter(double value, double cutoff, double maximum) noexcept {
  value = std::abs(value);
  if (value < cutoff)
    return 0.0;
  if (value > maximum)
    return 1.0;
  if (maximum == cutoff)
    return value >= maximum ? 1.0 : 0.0;
  const double ratio = (value - maximum) / (value - cutoff);
  return std::exp(-(ratio * ratio));
}

bool containsForcingType(const WVForcingFactoryRegistration &registration,
                         const char *type) {
  return std::find(registration.forcingTypes.begin(),
                   registration.forcingTypes.end(), type) !=
         registration.forcingTypes.end();
}

const char *qgForcingType(WVForcingStage stage) noexcept {
  switch (stage) {
  case WVForcingStage::spatial:
    return "PVSpatial";
  case WVForcingStage::spectral:
    return "PVSpectral";
  case WVForcingStage::spectralAmplitude:
    return "PVSpectralAmplitude";
  }
  return "";
}

bool finitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool isSelfConjugate(std::int64_t mode, std::size_t count) noexcept {
  return mode == 0 ||
         (count % 2 == 0 &&
          mode == -static_cast<std::int64_t>(count / 2));
}

bool isPrimary(std::int64_t kMode, std::int64_t lMode, std::size_t Nx,
               std::size_t Ny) noexcept {
  const bool kSelf = isSelfConjugate(kMode, Nx);
  const bool lSelf = isSelfConjugate(lMode, Ny);
  return lMode > 0 || (lSelf && (kMode > 0 || kSelf));
}

bool isNyquist(std::int64_t kMode, std::int64_t lMode, std::size_t Nx,
                std::size_t Ny) noexcept {
  return (Nx % 2 == 0 &&
          kMode == -static_cast<std::int64_t>(Nx / 2)) ||
         (Ny % 2 == 0 &&
          lMode == -static_cast<std::int64_t>(Ny / 2));
}

WVKernelStatus allocationLightCoefficientCount(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::size_t &count) {
  count = 0;
  if (configuration.contractVersion != WVKernelContractVersion ||
      configuration.Nx < 2 || configuration.Ny < 2 || configuration.j > 1 ||
      !finitePositive(configuration.Lx) ||
      !finitePositive(configuration.Ly) ||
      !finitePositive(configuration.h) || !finitePositive(configuration.g) ||
      !finitePositive(configuration.planetaryRadius) ||
      !finitePositive(configuration.rotationRate) ||
      !std::isfinite(configuration.latitude) ||
      std::abs(configuration.latitude) < 5.0 ||
      std::abs(configuration.latitude) > 85.0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Invalid Barotropic QG configuration during forcing preflight."};
  const double maximumK = 2.0 * pi *
      static_cast<double>(configuration.Nx / 2) / configuration.Lx;
  const auto positiveK = (configuration.Nx + 1) / 2;
  const auto positiveL = (configuration.Ny + 1) / 2;
  for (std::size_t iL = 0; iL < configuration.Ny; ++iL) {
    const auto lMode = iL < positiveL
                           ? static_cast<std::int64_t>(iL)
                           : static_cast<std::int64_t>(iL) -
                                 static_cast<std::int64_t>(configuration.Ny);
    for (std::size_t iK = 0; iK < configuration.Nx; ++iK) {
      const auto kMode = iK < positiveK
                             ? static_cast<std::int64_t>(iK)
                             : static_cast<std::int64_t>(iK) -
                                   static_cast<std::int64_t>(configuration.Nx);
      if (!isPrimary(kMode, lMode, configuration.Nx, configuration.Ny) ||
          isNyquist(kMode, lMode, configuration.Nx, configuration.Ny))
        continue;
      const double k = 2.0 * pi * static_cast<double>(kMode) /
                       configuration.Lx;
      const double l = 2.0 * pi * static_cast<double>(lMode) /
                       configuration.Ly;
      if (configuration.shouldAntialias &&
          std::sqrt(k * k + l * l) > 2.0 * maximumK / 3.0)
        continue;
      if (count == std::numeric_limits<std::size_t>::max())
        return {WVKernelStatusCode::sizeOverflow,
                "Barotropic QG coefficient count overflowed."};
      ++count;
    }
  }
  if (count == 0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG forcing preflight retained no Fourier modes."};
  return WVKernelStatus::ok();
}

class ResolvedBarotropicQGForcing : public WVBarotropicQGForcing {
public:
  explicit ResolvedBarotropicQGForcing(const WVFrozenForcingEntry &entry)
      : typeIdentifier_(entry.typeIdentifier), name_(entry.name),
        contractVersion_(entry.contractVersion), stage_(entry.stage),
        priority_(entry.priority), ordinal_(entry.ordinal) {}
  const std::string &typeIdentifier() const noexcept override {
    return typeIdentifier_;
  }
  std::uint32_t contractVersion() const noexcept override {
    return contractVersion_;
  }
  const std::string &name() const noexcept override { return name_; }
  WVForcingStage stage() const noexcept override { return stage_; }
  std::uint8_t priority() const noexcept override { return priority_; }
  std::size_t ordinal() const noexcept override { return ordinal_; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + metadataDynamicBytes();
  }

protected:
  std::size_t metadataDynamicBytes() const noexcept {
    return typeIdentifier_.capacity() + name_.capacity();
  }

private:
  std::string typeIdentifier_;
  std::string name_;
  std::uint32_t contractVersion_ = 0;
  WVForcingStage stage_ = WVForcingStage::spatial;
  std::uint8_t priority_ = 255;
  std::size_t ordinal_ = 0;
};

class QGNonlinearAdvection final : public ResolvedBarotropicQGForcing {
public:
  using ResolvedBarotropicQGForcing::ResolvedBarotropicQGForcing;
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.nonlinearAdvection();
  }
};

class QGAdaptiveDamping final : public ResolvedBarotropicQGForcing {
public:
  QGAdaptiveDamping(WVFrozenForcingEntry entry, std::vector<double> damping)
      : ResolvedBarotropicQGForcing(entry), damping_(std::move(damping)) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.adaptiveDamping(damping_);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + metadataDynamicBytes() + vectorBytes(damping_);
  }

private:
  std::vector<double> damping_;
};

class QGLinearBottomFriction final : public ResolvedBarotropicQGForcing {
public:
  QGLinearBottomFriction(WVFrozenForcingEntry entry, double rate)
      : ResolvedBarotropicQGForcing(entry), rate_(rate) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.linearBottomFriction(rate_);
  }

private:
  double rate_ = 0.0;
};

class QGQuadraticBottomFriction final : public ResolvedBarotropicQGForcing {
public:
  QGQuadraticBottomFriction(WVFrozenForcingEntry entry, double drag)
      : ResolvedBarotropicQGForcing(entry), drag_(drag) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.quadraticBottomFriction(drag_);
  }

private:
  double drag_ = 0.0;
};

class QGBetaPlanePVAdvection final : public ResolvedBarotropicQGForcing {
public:
  QGBetaPlanePVAdvection(WVFrozenForcingEntry entry, double beta)
      : ResolvedBarotropicQGForcing(entry), beta_(beta) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.betaPlanePVAdvection(beta_);
  }

private:
  double beta_ = 0.0;
};

class QGFixedAmplitude final : public ResolvedBarotropicQGForcing {
public:
  QGFixedAmplitude(WVFrozenForcingEntry entry,
                   WVBarotropicQGFixedAmplitudeConfiguration values)
      : ResolvedBarotropicQGForcing(entry), values_(std::move(values)) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    context.zeroSelectedTendencies(values_);
    return WVKernelStatus::ok();
  }
  WVStateConstraintResult applyConstraint(WVComplexView &A0) const override {
    std::size_t modified = 0;
    for (std::size_t index = 0; index < values_.A0Indices.size(); ++index) {
      const auto destination = values_.A0Indices[index];
      const auto previous = A0.data[destination];
      const auto value = values_.A0Values[index];
      if (previous.real != value.real || previous.imag != value.imag)
        ++modified;
      A0.data[destination] = value;
    }
    return {WVKernelStatus::ok(), modified, false};
  }
  std::size_t constraintWriteCount() const noexcept override {
    return values_.A0Indices.size();
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + metadataDynamicBytes() +
           vectorBytes(values_.A0Indices) + vectorBytes(values_.A0Values);
  }

private:
  WVBarotropicQGFixedAmplitudeConfiguration values_;
};

WVKernelStatus decodeFixedAmplitude(
    const WVFrozenForcingEntry &entry, std::size_t coefficientCount,
    WVBarotropicQGFixedAmplitudeConfiguration *configuration) {
  const auto rejectNonQGFamily = [&](const char *indexName,
                                     const char *realName,
                                     const char *imagName) {
    const auto *indices = integerValues(entry.configuration, indexName);
    const auto *real = realValues(entry.configuration, realName);
    const auto *imag = realValues(entry.configuration, imagName);
    return (indices != nullptr && !indices->empty()) ||
           (real != nullptr && !real->empty()) ||
           (imag != nullptr && !imag->empty());
  };
  if (rejectNonQGFamily("ApIndices", "ApValuesReal", "ApValuesImag") ||
      rejectNonQGFamily("AmIndices", "AmValuesReal", "AmValuesImag"))
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG fixed-amplitude forcing accepts only A0."};
  const auto *indices = integerValues(entry.configuration, "A0Indices");
  const auto *real = realValues(entry.configuration, "A0ValuesReal");
  const auto *imag = realValues(entry.configuration, "A0ValuesImag");
  if (indices == nullptr && real == nullptr && imag == nullptr)
    return WVKernelStatus::ok();
  if (indices == nullptr || real == nullptr || imag == nullptr ||
      indices->size() != real->size() || real->size() != imag->size())
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG fixed-amplitude arrays must be complete and equal-length."};
  std::set<std::size_t> unique;
  for (std::size_t index = 0; index < indices->size(); ++index) {
    if ((*indices)[index] < 0 ||
        static_cast<std::size_t>((*indices)[index]) >= coefficientCount ||
        !std::isfinite((*real)[index]) || !std::isfinite((*imag)[index]))
      return {WVKernelStatusCode::invalidConfiguration,
              "Barotropic QG fixed-amplitude state is outside compact A0 or nonfinite."};
    const auto converted = static_cast<std::size_t>((*indices)[index]);
    if (!unique.insert(converted).second)
      return {WVKernelStatusCode::invalidConfiguration,
              "Barotropic QG fixed-amplitude A0 indices must be unique."};
    if (configuration != nullptr) {
      configuration->A0Indices.push_back(converted);
      configuration->A0Values.push_back(
          {(*real)[index], (*imag)[index]});
    }
  }
  return WVKernelStatus::ok();
}

std::vector<double> adaptiveDampingOperator(
    const WVTransformBarotropicQGDescriptor &descriptor,
    WVKernelStatus &status) {
  const auto &configuration = descriptor.configuration();
  double maximumComponent = 0.0;
  for (const auto &horizontal : descriptor.fourierModes())
    maximumComponent = std::max(
        maximumComponent,
        std::max(std::abs(horizontal.k), std::abs(horizontal.l)));
  if (!(maximumComponent > 0.0)) {
    status = {WVKernelStatusCode::invalidConfiguration,
              "Adaptive damping requires resolved horizontal wavenumbers."};
    return {};
  }
  const double effectiveResolution = pi / maximumComponent;
  const double dklMinimum =
      std::min(2.0 * pi / configuration.Lx,
               2.0 * pi / configuration.Ly);
  const double klCutoff =
      dklMinimum * std::pow(maximumComponent / dklMinimum, 0.75);
  const double prefactor = effectiveResolution / (pi * pi);
  std::vector<double> damping(descriptor.Nkl());
  for (std::size_t index = 0; index < descriptor.Nkl(); ++index) {
    const auto &horizontal = descriptor.fourierModes()[index];
    const double filter =
        vanishingFilter(horizontal.Kh, klCutoff, maximumComponent);
    damping[index] =
        -prefactor * filter *
        (horizontal.k * horizontal.k + horizontal.l * horizontal.l);
  }
  status = WVKernelStatus::ok();
  return damping;
}

} // namespace

namespace detail {

WVKernelStatus preflightBarotropicQGEmptyForcing(
    const WVFrozenForcingEntry &entry, std::size_t) {
  return emptyConfiguration(entry)
             ? WVKernelStatus::ok()
             : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                              "This Barotropic QG forcing accepts no configuration values."};
}

WVKernelStatus preflightBarotropicQGFixedAmplitude(
    const WVFrozenForcingEntry &entry, std::size_t coefficientCount) {
  return decodeFixedAmplitude(entry, coefficientCount, nullptr);
}

WVKernelStatus preflightBarotropicQGScalarForcing(
    const WVFrozenForcingEntry &entry, std::size_t) {
  if (entry.configuration.values.size() != 1 ||
      !entry.configuration.values.front().dimensions.empty())
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG scalar forcing requires exactly one scalar value."};
  const auto *values = std::get_if<std::vector<double>>(
      &entry.configuration.values.front().storage);
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()))
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG scalar forcing requires one finite real value."};
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGNonlinearAdvection(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  forcing = std::make_unique<QGNonlinearAdvection>(entry);
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGAdaptiveDamping(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &descriptor, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  WVKernelStatus status;
  auto damping = adaptiveDampingOperator(descriptor, status);
  if (!status)
    return status;
  forcing =
      std::make_unique<QGAdaptiveDamping>(entry, std::move(damping));
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGFixedAmplitude(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &descriptor, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  WVBarotropicQGFixedAmplitudeConfiguration values;
  auto status = decodeFixedAmplitude(entry, descriptor.Nkl(), &values);
  if (!status)
    return status;
  forcing =
      std::make_unique<QGFixedAmplitude>(entry, std::move(values));
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGQuadraticBottomFriction(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  const auto *values = realValues(entry.configuration, "Cd");
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()) || values->front() < 0.0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic quadratic drag requires one finite nonnegative Cd."};
  forcing = std::make_unique<QGQuadraticBottomFriction>(
      entry, values->front() / 4000.0);
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGLinearBottomFriction(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  const auto *values = realValues(entry.configuration, "r");
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()) || values->front() < 0.0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic linear drag requires one finite nonnegative r."};
  forcing =
      std::make_unique<QGLinearBottomFriction>(entry, values->front());
  return WVKernelStatus::ok();
}

WVKernelStatus createBarotropicQGBetaPlanePVAdvection(
    const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &descriptor, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  const auto &configuration = descriptor.configuration();
  const double beta = 2.0 * configuration.rotationRate *
                      std::cos(configuration.latitude * pi / 180.0) /
                      configuration.planetaryRadius;
  forcing = std::make_unique<QGBetaPlanePVAdvection>(entry, beta);
  return WVKernelStatus::ok();
}

} // namespace detail

WVKernelStatus WVBarotropicQGForcingExecutionContext::nonlinearAdvection() {
  const auto status = engine_->kernel().addPotentialVorticityAdvection(
      A0_, F0_, outputInitialized_, workspace_);
  if (status)
    outputInitialized_ = true;
  return status;
}

WVKernelStatus WVBarotropicQGForcingExecutionContext::adaptiveDamping(
    const std::vector<double> &dampingOperator) {
  const auto status = engine_->kernel().addAdaptiveDamping(
      A0_, dampingOperator, F0_, outputInitialized_, workspace_);
  if (status)
    outputInitialized_ = true;
  return status;
}

WVKernelStatus WVBarotropicQGForcingExecutionContext::linearBottomFriction(
    double rate) {
  const auto status = engine_->kernel().addLinearBottomFriction(
      A0_, rate, F0_, outputInitialized_, workspace_);
  if (status)
    outputInitialized_ = true;
  return status;
}

WVKernelStatus
WVBarotropicQGForcingExecutionContext::quadraticBottomFriction(double drag) {
  const auto status = engine_->kernel().addQuadraticBottomFriction(
      A0_, drag, F0_, outputInitialized_, workspace_);
  if (status)
    outputInitialized_ = true;
  return status;
}

WVKernelStatus WVBarotropicQGForcingExecutionContext::betaPlanePVAdvection(
    double beta) {
  const auto status = engine_->kernel().addBetaPlanePVAdvection(
      A0_, beta, F0_, outputInitialized_, workspace_);
  if (status)
    outputInitialized_ = true;
  return status;
}

void WVBarotropicQGForcingExecutionContext::zeroSelectedTendencies(
    const WVBarotropicQGFixedAmplitudeConfiguration &configuration) {
  if (!outputInitialized_) {
    engine_->initializeOutputWithZeros(F0_);
    outputInitialized_ = true;
  }
  for (const auto index : configuration.A0Indices)
    F0_.data[index] = {};
}

WVBarotropicQGForcingEngine::~WVBarotropicQGForcingEngine() = default;

WVKernelStatus WVBarotropicQGForcingEngine::validateSchedule(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule, std::size_t coefficientCount,
    const WVExtensionCatalog &catalog) {
  if (schedule.profileIdentifier != WVForcingScheduleProfileIdentifier ||
      schedule.profileVersion != WVForcingScheduleProfileVersion)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported frozen forcing schedule profile."};
  if (coefficientCount == 0)
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG forcing requires nonempty compact A0."};
  std::set<std::string> names;
  for (const auto &entry : schedule.entries) {
    const auto *registration = catalog.forcings().registration(
        entry.typeIdentifier, entry.contractVersion);
    if (registration == nullptr || !registration->isSupported ||
        !registration->barotropicQGFactory ||
        registration->contractVersion != entry.contractVersion)
      return {WVKernelStatusCode::unsupportedOperation,
              "The frozen schedule has no matching Barotropic QG forcing implementation."};
    if (entry.stage != registration->barotropicQGStage ||
        !containsForcingType(*registration, qgForcingType(entry.stage)))
      return {WVKernelStatusCode::invalidConfiguration,
              "A Barotropic QG forcing record is assigned to an incompatible stage."};
    if (entry.name.empty() || !names.insert(entry.name).second)
      return {WVKernelStatusCode::invalidConfiguration,
              "Forcing names must be nonempty and unique."};
    auto status = catalog.forcings().validateConfiguration(entry);
    if (!status)
      return status;
    if (!registration->barotropicQGPreflight)
      return {WVKernelStatusCode::invalidConfiguration,
              "A Barotropic QG forcing registration lacks allocation-light preflight."};
    status = registration->barotropicQGPreflight(entry, coefficientCount);
    if (!status)
      return status;
  }
  (void)configuration;
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> fftEngine,
    std::unique_ptr<WVBarotropicQGForcingEngine> &forcingEngine) {
  forcingEngine.reset();
  if (!catalog)
    return {WVKernelStatusCode::invalidPointer,
            "Barotropic QG forcing construction requires an extension catalog."};
  if (!fftEngine)
    return {WVKernelStatusCode::invalidPointer,
            "Barotropic QG forcing construction requires an FFT engine."};
  std::size_t coefficientCount = 0;
  auto status = allocationLightCoefficientCount(configuration,
                                                 coefficientCount);
  if (!status)
    return status;
  status = validateSchedule(configuration, schedule, coefficientCount,
                            *catalog);
  if (!status)
    return status;
  try {
    auto candidate = std::unique_ptr<WVBarotropicQGForcingEngine>(
        new WVBarotropicQGForcingEngine());
    candidate->catalog_ = std::move(catalog);
    status = WVTransformBarotropicQGKernel::create(
        configuration, std::move(fftEngine), candidate->kernel_);
    if (!status)
      return status;
    status = candidate->initialize(schedule);
    if (!status)
      return status;
    forcingEngine = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Barotropic QG forcing-engine allocation failed."};
  } catch (const std::exception &error) {
    return {WVKernelStatusCode::invalidConfiguration, error.what()};
  }
}

WVKernelStatus WVBarotropicQGForcingEngine::initialize(
    const WVFrozenForcingSchedule &schedule) {
  std::vector<const WVFrozenForcingEntry *> entries;
  entries.reserve(schedule.entries.size());
  for (const auto &entry : schedule.entries)
    entries.push_back(&entry);
  std::stable_sort(entries.begin(), entries.end(),
                   [](const auto *left, const auto *right) {
                     if (stageRank(left->stage) != stageRank(right->stage))
                       return stageRank(left->stage) <
                              stageRank(right->stage);
                     if (left->priority != right->priority)
                       return left->priority < right->priority;
                     return left->ordinal < right->ordinal;
                   });
  const bool hasAdaptiveDamping = std::any_of(
      entries.begin(), entries.end(), [this](const auto *entry) {
        const auto *registration = catalog_->forcings().registration(
            entry->typeIdentifier, entry->contractVersion);
        return registration != nullptr &&
               registration->providesAdaptiveDamping;
      });
  for (const auto *entry : entries) {
    std::unique_ptr<WVBarotropicQGForcing> resolved;
    auto status = catalog_->forcings().createBarotropicQG(
        *entry, kernel_->descriptor(), hasAdaptiveDamping, resolved);
    if (!status)
      return status;
    if (!resolved || resolved->stage() != entry->stage)
      return {WVKernelStatusCode::invalidConfiguration,
              "A Barotropic QG forcing factory returned an incompatible implementation."};
    if (entry->stage == WVForcingStage::spatial)
      ++metrics_.resolvedSpatialCount;
    else if (entry->stage == WVForcingStage::spectral)
      ++metrics_.resolvedSpectralCount;
    else
      ++metrics_.resolvedAmplitudeCount;
    metrics_.derivedOperatorBytes += resolved->persistentBytes();
    forcing_.push_back(std::move(resolved));
  }
  std::ostringstream identifier;
  identifier << WVForcingScheduleProfileIdentifier << ':';
  for (std::size_t index = 0; index < entries.size(); ++index) {
    if (index != 0)
      identifier << ',';
    identifier << entries[index]->typeIdentifier;
  }
  scheduleIdentifier_ = identifier.str();
  metrics_.scheduleBytes =
      scheduleIdentifier_.capacity() +
      forcing_.capacity() * sizeof(std::unique_ptr<WVBarotropicQGForcing>);
  metrics_.workspaceCapacityBytes = 0;
  return WVKernelStatus::ok();
}

void WVBarotropicQGForcingEngine::initializeOutputWithZeros(
    WVComplexView &F0) {
  std::fill_n(F0.data, kernel_->descriptor().Nkl(), WVComplex64{});
}

WVKernelStatus WVBarotropicQGForcingEngine::evaluateRightHandSide(
    const WVComplexConstView &A0, WVComplexView &F0,
    WVRealFieldBundleConstView *advectionFields) {
  if (advectionFields != nullptr)
    *advectionFields = {};
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Barotropic QG forcing-engine execution is not reentrant."};
  const auto expected = kernel_->descriptor().spectralShape();
  if (A0.shape.rows != expected.rows || A0.shape.columns != expected.columns ||
      F0.shape.rows != expected.rows || F0.shape.columns != expected.columns)
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG forcing requires compact [1,Nkl] A0 and F0."};
  if (A0.data == nullptr || F0.data == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Barotropic QG forcing received null compact state storage."};
  if (A0.data == F0.data)
    return {WVKernelStatusCode::overlappingArrays,
            "Barotropic QG A0 and F0 must not overlap."};
  executing_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{executing_};
  WVBarotropicQGForcingExecutionContext context;
  context.engine_ = this;
  context.A0_ = A0;
  context.F0_ = F0;
  for (const auto &forcing : forcing_) {
    ++metrics_.forcingCallCount;
    const auto status = forcing->addRightHandSide(context);
    if (!status)
      return status;
  }
  if (!context.outputInitialized_)
    initializeOutputWithZeros(F0);
  if (advectionFields != nullptr) {
    const auto status = kernel_->prepareAdvectionFields(
        A0, context.workspace_, *advectionFields);
    if (!status)
      return status;
  }
  metrics_.physicalFieldReconstructionCount +=
      context.workspace_.physicalFieldReconstructionCount;
  metrics_.physicalFieldReuseCount +=
      context.workspace_.physicalFieldReuseCount;
  metrics_.spatialTendencyProjectionCount +=
      context.workspace_.spatialTendencyProjectionCount;
  ++metrics_.evaluationCount;
  return WVKernelStatus::ok();
}

WVStateConstraintResult
WVBarotropicQGForcingEngine::restoreForcingAmplitudes(WVComplexView &A0) {
  const auto expected = kernel_->descriptor().spectralShape();
  if (A0.shape.rows != expected.rows || A0.shape.columns != expected.columns)
    return {{WVKernelStatusCode::invalidShape,
             "Barotropic QG constraints require compact [1,Nkl] A0."},
            0, false};
  if (A0.data == nullptr)
    return {{WVKernelStatusCode::invalidPointer,
             "Barotropic QG constraint A0 storage is null."},
            0, false};
  ++metrics_.constraintOperationCount;
  std::size_t modified = 0;
  bool fsalCompatible = true;
  for (const auto &forcing : forcing_) {
    const auto result = forcing->applyConstraint(A0);
    if (!result.status)
      return result;
    modified += result.modifiedCoefficientCount;
    fsalCompatible = fsalCompatible && result.fsalCompatible;
    metrics_.restoredCoefficientCount += forcing->constraintWriteCount();
    metrics_.stateConstraintElementWrites += forcing->constraintWriteCount();
  }
  return {WVKernelStatus::ok(), modified, fsalCompatible};
}

std::size_t WVBarotropicQGForcingEngine::persistentBytes() const noexcept {
  return sizeof(*this) +
         (kernel_ == nullptr ? 0 : kernel_->persistentBytes()) +
         metrics_.scheduleBytes + metrics_.derivedOperatorBytes;
}

} // namespace wavevortex::runtime
