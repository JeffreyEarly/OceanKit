#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"
#include "WVForcingDiagnosticWorkspace.hpp"

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

class ScopedBarotropicQGEvaluation final {
public:
  ScopedBarotropicQGEvaluation(WVBarotropicQGForcingEngine &engine,
                               const WVComplexConstView &state)
      : engine_(engine) {
    if (engine_.stateEvaluationActive())
      status_ = engine_.validateStateEvaluation(state);
    else {
      status_ = engine_.beginStateEvaluation(state);
      owns_ = static_cast<bool>(status_);
    }
  }
  ~ScopedBarotropicQGEvaluation() {
    if (owns_)
      (void)engine_.endStateEvaluation();
  }
  const WVKernelStatus &status() const noexcept { return status_; }

private:
  WVBarotropicQGForcingEngine &engine_;
  bool owns_ = false;
  WVKernelStatus status_ = WVKernelStatus::ok();
};

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
  bool supportsTendencyDiagnostics() const noexcept override { return true; }
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
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  using ResolvedBarotropicQGForcing::ResolvedBarotropicQGForcing;
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.nonlinearAdvection();
  }
};

class QGAdaptiveDamping final : public ResolvedBarotropicQGForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
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
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
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
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
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
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  QGBetaPlanePVAdvection(WVFrozenForcingEntry entry, double beta)
      : ResolvedBarotropicQGForcing(entry), beta_(beta) {}
  WVKernelStatus addRightHandSide(
      WVBarotropicQGForcingExecutionContext &context) const override {
    return context.betaPlanePVAdvection(beta_);
  }

private:
  double beta_ = 0.0;
};

class QGExplicitAntialiasing final : public ResolvedBarotropicQGForcing {
public:
  QGExplicitAntialiasing(const WVFrozenForcingEntry &entry, std::vector<std::size_t> indices)
      : ResolvedBarotropicQGForcing(entry), indices_(std::move(indices)) {}
  WVKernelStatus addRightHandSide(WVBarotropicQGForcingExecutionContext &context) const override {
    context.filterTendency(indices_);
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes()+vectorBytes(indices_); }
private:
  std::vector<std::size_t> indices_;
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

WVKernelStatus preflightBarotropicQGExplicitAntialiasing(const WVFrozenForcingEntry &entry, std::size_t) {
  return preflightExplicitAntialiasing(entry,false);
}
WVKernelStatus createBarotropicQGExplicitAntialiasing(const WVFrozenForcingEntry &entry,
    const WVTransformBarotropicQGDescriptor &descriptor, bool,
    std::unique_ptr<WVBarotropicQGForcing> &forcing) {
  const auto &c = descriptor.configuration();
  auto status = preflightExplicitAntialiasing(entry,c.shouldAntialias);
  if (!status) return status;
  const double Nj = realValues(entry.configuration,"Nj")->front();
  const double maximumK = 2.0*pi*static_cast<double>(c.Nx/2)/c.Lx;
  std::vector<std::size_t> indices;
  for (std::size_t kl = 0; kl < descriptor.Nkl(); ++kl)
    if (descriptor.fourierModes()[kl].Kh > 2.0*maximumK/3.0 || c.j > Nj-1.0)
      indices.push_back(kl);
  forcing = std::make_unique<QGExplicitAntialiasing>(entry,std::move(indices));
  return WVKernelStatus::ok();
}

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

void WVBarotropicQGForcingExecutionContext::filterTendency(const std::vector<std::size_t> &indices) {
  if (!outputInitialized_) {
    engine_->initializeOutputWithZeros(F0_);
    outputInitialized_ = true;
  }
  for (const auto index : indices) F0_.data[index] = {};
}

WVKernelStatus WVBarotropicQGForcingExecutionContext::nonlinearAdvection() {
  if (workspace_.spatialTendency != nullptr) {
    const auto status = engine_->kernel().addPotentialVorticityAdvection(
        A0_, F0_, outputInitialized_, workspace_);
    if (status)
      outputInitialized_ = true;
    return status;
  }
  const WVVariableEvaluationKey nonlinear{
      WVVariableEvaluationNode::forcingTendency, 0};
  if (engine_->evaluationPolicy_ == WVVariableEvaluationPolicy::lowMemory) {
    const auto status = engine_->evaluation_.evaluate(
        nonlinear, engine_->nonlinearScratch_.size() * sizeof(WVComplex64),
        [&] {
          return engine_->kernel().addPotentialVorticityAdvection(
              A0_, F0_, outputInitialized_, workspace_);
        });
    if (status)
      outputInitialized_ = true;
    (void)engine_->evaluation_.evict(nonlinear);
    return status;
  }
  const auto status = engine_->evaluation_.evaluate(
      nonlinear, engine_->nonlinearScratch_.size() * sizeof(WVComplex64), [&] {
        WVComplexView output{engine_->nonlinearScratch_.data(), F0_.shape};
        return engine_->kernel().addPotentialVorticityAdvection(
            A0_, output, false, workspace_);
      });
  if (!status)
    return status;
  if (!outputInitialized_) {
    std::copy(engine_->nonlinearScratch_.begin(),
              engine_->nonlinearScratch_.end(), F0_.data);
    outputInitialized_ = true;
  } else {
    for (std::size_t index = 0; index < engine_->nonlinearScratch_.size();
         ++index) {
      F0_.data[index].real += engine_->nonlinearScratch_[index].real;
      F0_.data[index].imag += engine_->nonlinearScratch_[index].imag;
    }
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingExecutionContext::adaptiveDamping(
    const std::vector<double> &dampingOperator) {
  auto producer = [&](double& maximum) {
    const auto before=engine_->kernel().metrics().horizontalSpeedMaximumReductionCount;
    const auto status=engine_->kernel().horizontalSpeedMaximum(A0_,maximum,workspace_);
    engine_->metrics_.horizontalSpeedReductionCount+=
        engine_->kernel().metrics().horizontalSpeedMaximumReductionCount-before;
    return status;
  };
  const auto status = diagnosticWorkspace_ ?
      diagnosticWorkspace_->evaluateHorizontalMaximum(
          engine_->horizontalSpeedMaximum_,producer) :
      engine_->evaluation_.evaluate(
          {WVVariableEvaluationNode::reduction, 0}, sizeof(double), [&] {
            return producer(engine_->horizontalSpeedMaximum_);
          });
  if (!status)
    return status;
  workspace_.horizontalSpeedMaximum=engine_->horizontalSpeedMaximum_;
  workspace_.horizontalSpeedMaximumPrepared=true;
  const auto dampingStatus = engine_->kernel().addAdaptiveDamping(
      A0_, dampingOperator, F0_, outputInitialized_, workspace_);
  if (dampingStatus)
    outputInitialized_ = true;
  return dampingStatus;
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

const WVForcingEvaluationDependencies*
WVBarotropicQGForcingEngine::forcingEvaluationDependencies(
    std::size_t index) const noexcept {
  const auto* forcing=forcingInstance(index);
  if(!forcing || !catalog_) return nullptr;
  const auto* registration=catalog_->forcings().registration(
      forcing->typeIdentifier(),forcing->contractVersion());
  return registration ? &registration->evaluationDependencies : nullptr;
}

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
    if (registration->modelPreflight) {
      status = registration->modelPreflight(entry,configuration.shouldAntialias);
      if (!status) return status;
    }
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
  nonlinearScratch_.resize(kernel_->descriptor().Nkl());
  auto status = evaluation_.prepare(
      {{WVVariableEvaluationNode::physicalField, 0},
       {WVVariableEvaluationNode::reduction, 0},
       {WVVariableEvaluationNode::forcingTendency, 0}});
  if (!status)
    return status;
  metrics_.scheduleBytes =
      scheduleIdentifier_.capacity() +
      forcing_.capacity() * sizeof(std::unique_ptr<WVBarotropicQGForcing>);
  metrics_.workspaceCapacityBytes = vectorBytes(nonlinearScratch_);
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::beginStateEvaluation(
    const WVComplexConstView &A0) {
  if (evaluation_.active() || executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Barotropic QG state evaluation is already active."};
  const bool ownsKernelScope=!kernel_->stateEvaluationActive();
  auto status = ownsKernelScope ? kernel_->beginStateEvaluation(A0) :
                                  kernel_->validateStateEvaluation(A0);
  if (!status)
    return status;
  status = evaluation_.begin(this, evaluationPolicy_);
  if (!status) {
    if (ownsKernelScope) (void)kernel_->endStateEvaluation();
    return status;
  }
  evaluationOwnsKernelScope_=ownsKernelScope;
  evaluationState_ = A0;
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::validateVariableEvaluationPolicyChange(
    WVVariableEvaluationPolicy policy) const noexcept {
  if (executing_ || evaluation_.active())
    return {WVKernelStatusCode::reentrantExecution,
            "Cannot change an active evaluation policy."};
  if (policy != WVVariableEvaluationPolicy::reuse &&
      policy != WVVariableEvaluationPolicy::lowMemory)
    return {WVKernelStatusCode::invalidConfiguration,
            "Unknown variable evaluation policy."};
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::setVariableEvaluationPolicy(
    WVVariableEvaluationPolicy policy) {
  auto status=validateVariableEvaluationPolicyChange(policy);
  if(!status) return status;
  if (policy == WVVariableEvaluationPolicy::reuse &&
      nonlinearScratch_.size() != kernel_->descriptor().Nkl()) {
    std::vector<WVComplex64> prepared;
    try {
      prepared.resize(kernel_->descriptor().Nkl());
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate the Barotropic QG nonlinear cache."};
    }
    nonlinearScratch_.swap(prepared);
  } else if (policy == WVVariableEvaluationPolicy::lowMemory) {
    std::vector<WVComplex64>().swap(nonlinearScratch_);
  }
  evaluationPolicy_ = policy;
  metrics_.workspaceCapacityBytes = vectorBytes(nonlinearScratch_);
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::endStateEvaluation() {
  if (!evaluation_.active())
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG state evaluation is not active."};
  evaluation_.end();
  evaluationState_ = {};
  const bool ownsKernelScope=evaluationOwnsKernelScope_;
  evaluationOwnsKernelScope_=false;
  return ownsKernelScope ? kernel_->endStateEvaluation() : WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::validateStateEvaluation(
    const WVComplexConstView &A0) const noexcept {
  if (!evaluation_.active() || A0.data != evaluationState_.data ||
      A0.shape.rows != evaluationState_.shape.rows ||
      A0.shape.columns != evaluationState_.shape.columns)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic QG coefficients do not belong to the active evaluation."};
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGForcingEngine::horizontalSpeedMaximum(
    const WVComplexConstView &A0, double &maximum) {
  ScopedBarotropicQGEvaluation evaluation(*this, A0);
  if (!evaluation.status())
    return evaluation.status();
  auto status = evaluation_.evaluate(
      {WVVariableEvaluationNode::reduction, 0}, sizeof(double), [&] {
        WVBarotropicQGOperationWorkspace workspace;
        const auto before=kernel_->metrics().horizontalSpeedMaximumReductionCount;
        const auto status=kernel_->horizontalSpeedMaximum(
            A0, horizontalSpeedMaximum_, workspace);
        metrics_.horizontalSpeedReductionCount+=
            kernel_->metrics().horizontalSpeedMaximumReductionCount-before;
        return status;
      });
  if (status)
    maximum = horizontalSpeedMaximum_;
  return status;
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
  ScopedBarotropicQGEvaluation evaluation(*this, A0);
  if (!evaluation.status())
    return evaluation.status();
  executing_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{executing_};
  WVBarotropicQGForcingExecutionContext context;
  context.engine_ = this;
  context.A0_ = A0;
  context.F0_ = F0;
  if (!linearDynamics_) for (const auto &forcing : forcing_) {
    ++metrics_.forcingCallCount;
    const auto status = forcing->addRightHandSide(context);
    if (!status)
      return status;
  }
  if (!context.outputInitialized_)
    initializeOutputWithZeros(F0);
  if (advectionFields != nullptr) {
    const auto status = evaluation_.evaluate(
        {WVVariableEvaluationNode::physicalField, 0},
        2 * kernel_->descriptor().spatialShape().elementCount() *
            sizeof(double),
        [&] {
          return kernel_->prepareAdvectionFields(
              A0, context.workspace_, *advectionFields);
        });
    if (!status)
      return status;
    if (advectionFields->data == nullptr) {
      const auto status = kernel_->prepareAdvectionFields(
          A0, context.workspace_, *advectionFields);
      if (!status)
        return status;
    }
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

WVKernelStatus WVBarotropicQGForcingEngine::evaluateForcingTendencies(
    const WVComplexConstView& A0,const WVForcingTendencyOutput* outputs,
    std::size_t count,const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session) {
  if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing diagnostics require an idle engine."};
  tendencyMetrics_.workspaceLastPeakBytes=0;
  const auto spectral=kernel().descriptor().spectralShape();
  const auto plane=kernel().descriptor().spatialShape();
  const WVShape4D spatial{plane.rows,plane.columns,1,1};
  const WVState state{0,0,{{},{},A0}};
  auto status=detail::validateForcingTendencyOutputs(forcing_,spectral,spatial,state,outputs,count);
  if (!status || !count) return status;
  status=detail::validatePreparedDiagnosticFields(preparedPhysical,spatial,2,state,outputs,count);
  if (!status) return status;
  if (A0.shape.rows!=spectral.rows || A0.shape.columns!=spectral.columns)
    return {WVKernelStatusCode::invalidShape,"Expected compact QG diagnostic state."};
  const auto address=reinterpret_cast<std::uintptr_t>(A0.data);
  if (!address || address%alignof(WVComplex64) || spectral.elementCount()*sizeof(WVComplex64)>UINTPTR_MAX-address)
    return {WVKernelStatusCode::invalidPointer,"Invalid QG diagnostic state storage."};
  ScopedBarotropicQGEvaluation evaluation(*this, A0);
  if (!evaluation.status())
    return evaluation.status().code==WVKernelStatusCode::numericalFailure ?
        WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
            "QG diagnostic state must be finite."} : evaluation.status();
  try {
    detail::WVForcingDiagnosticLedger localLedger(tendencyMetrics_);
    std::unique_ptr<detail::WVForcingDiagnosticWorkspace> local;
    if (!session) {
      local=std::make_unique<detail::WVForcingDiagnosticWorkspace>(spectral,spatial,1,2);
      std::vector<WVForcingStage> stages;
      local->nonlinearUseCount=0;
      for(std::size_t index=0;index<forcing_.size();++index) {
        const auto& forcing=forcing_[index];
        stages.push_back(forcing->stage());
        const auto* dependencies=forcingEvaluationDependencies(index);
        if(!dependencies) return {WVKernelStatusCode::invalidConfiguration,
            "Forcing evaluation dependencies are unavailable."};
        local->nonlinearUseCount+=dependencies->nonlinearUseCount;
      }
      status=localLedger.context.prepare(
          detail::WVForcingDiagnosticWorkspace::dependencyKeys(forcing_.size()));
      if(!status) return status;
      status=localLedger.context.begin(this,evaluationPolicy_); if(!status) return status;
      status=local->beginScopedEvaluation(localLedger.context,stages); if(!status) return status;
      session=local.get();
    }
    auto& work=*session;
    status=work.bind(this,state); if (!status) return status;
    const auto S=spectral.elementCount();
    const auto R=plane.elementCount();
    if (work.spectral.rows!=spectral.rows || work.spectral.columns!=spectral.columns ||
        work.spatial.first!=spatial.first || work.spatial.second!=spatial.second ||
        work.spatial.third!=spatial.third || work.spatial.fourth!=spatial.fourth ||
        work.flux.size()!=S || work.previous.size()!=S || work.temporary.size()!=S ||
        work.cumulative.size()!=spatial.elementCount() || work.raw.size()!=spatial.elementCount() ||
        work.physical.size()!=2*R)
      return {WVKernelStatusCode::invalidShape,
              "Forcing diagnostic session has incompatible Barotropic QG storage."};
    auto flux=work.fluxView();
    if (!work.initialized()) {
      if (preparedPhysical) {
        std::copy_n(preparedPhysical->data,preparedPhysical->shape.elementCount(),
                    work.physical.data());
        work.physicalPrepared=true;
      }
      status=kernel().evolveA0(A0,0,flux.F0); if (!status) return status;
      std::fill(work.flux.begin(),work.flux.end(),WVComplex64{});
      work.markInitialized();
    }
    WVBarotropicQGForcingExecutionContext context;
    context.engine_=this; context.A0_=A0; context.outputInitialized_=true;
    context.diagnosticWorkspace_=&work;
    WVRealFieldBundleConstView persistentVelocity{
        work.physical.data(),{plane.rows,plane.columns,1,2}};
    if (work.physicalPrepared)
      context.workspace_.preparedVelocity=&persistentVelocity;
    executing_=true;
    struct Guard {
      WVBarotropicQGForcingEngine& engine;
      detail::WVForcingDiagnosticWorkspace& work;
      WVBarotropicQGOperationWorkspace& physical;
      ~Guard() {
        engine.tendencyMetrics_.workspaceLastPeakBytes=work.bytes();
        engine.tendencyMetrics_.workspaceHighWaterBytes=std::max(engine.tendencyMetrics_.workspaceHighWaterBytes,work.bytes());
        engine.tendencyMetrics_.workspaceLiveBytes=0;
        engine.metrics_.physicalFieldReconstructionCount+=physical.physicalFieldReconstructionCount;
        engine.metrics_.physicalFieldReuseCount+=physical.physicalFieldReuseCount;
        engine.executing_=false;
      }
    } guard{*this,work,context.workspace_};
    tendencyMetrics_.workspaceLiveBytes=work.bytes();
    WVRealView raw{work.raw.data(),plane};
    return detail::evaluateForcingTendencySequence(forcing_,work,outputs,count,tendencyMetrics_,
      [&](const WVBarotropicQGForcing& forcing,WVFlux& destination) {
        context.F0_=destination.F0;
        context.workspace_.spatialTendency=forcing.stage()==WVForcingStage::spatial ? &raw : nullptr;
        context.workspace_.spatialTendencyCaptured=false;
        if (forcing.requiresDiagnosticPhysicalFields() &&
            !work.physicalPrepared) {
          WVRealFieldBundleConstView fields;
          const auto prepared=kernel().prepareAdvectionFields(A0,context.workspace_,fields);
          if (!prepared) return prepared;
          std::copy_n(fields.data,2*R,work.physical.data());
          work.physicalPrepared=true;
          context.workspace_.preparedVelocity=&persistentVelocity;
        }
        const auto result=forcing.addRightHandSide(context);
        work.spatialCaptured=context.workspace_.spatialTendencyCaptured;
        return result;
      },
      [&](WVRealFieldBundleConstView fields,WVFlux& destination) {
        return kernel().transformQGPVToA0({fields.data,plane},destination.F0);
      },
      [&](const std::vector<WVComplex64>& difference,WVRealFieldBundleView destination) {
        WVRealView output{destination.data,plane};
        return kernel().transformSpectralTendencyToSpatial({difference.data(),spectral},output);
      });
  } catch (const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,"Unable to allocate event forcing diagnostic workspace."};
  }
}

WVStateConstraintResult
WVBarotropicQGForcingEngine::restoreForcingAmplitudes(WVComplexView &A0) {
  if (evaluation_.active() || executing_)
    return {{WVKernelStatusCode::reentrantExecution,
             "Barotropic QG constraints cannot mutate an active evaluation."},
            0, false};
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
         metrics_.scheduleBytes + metrics_.derivedOperatorBytes +
         metrics_.workspaceCapacityBytes +
         evaluation_.persistentBytes();
}

} // namespace wavevortex::runtime
