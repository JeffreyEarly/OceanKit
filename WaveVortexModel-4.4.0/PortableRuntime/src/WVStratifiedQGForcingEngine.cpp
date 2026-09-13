#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"

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

class ScopedStratifiedQGEvaluation final {
public:
  ScopedStratifiedQGEvaluation(WVStratifiedQGForcingEngine &engine,
                              WVComplexConstView state)
      : engine_(engine) {
    if (engine_.stateEvaluationActive())
      status_ = engine_.validateStateEvaluation(state);
    else {
      status_ = engine_.beginStateEvaluation(state);
      owns_ = static_cast<bool>(status_);
    }
  }
  ~ScopedStratifiedQGEvaluation() {
    if (owns_)
      (void)engine_.endStateEvaluation();
  }
  const WVKernelStatus &status() const noexcept { return status_; }

private:
  WVStratifiedQGForcingEngine &engine_;
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

class ResolvedStratifiedQGForcing : public WVStratifiedQGForcing {
public:
  explicit ResolvedStratifiedQGForcing(const WVFrozenForcingEntry &entry)
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

class QGNonlinearAdvection final : public ResolvedStratifiedQGForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  using ResolvedStratifiedQGForcing::ResolvedStratifiedQGForcing;
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
    return context.nonlinearAdvection();
  }
};

class QGAdaptiveDamping final : public ResolvedStratifiedQGForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  QGAdaptiveDamping(WVFrozenForcingEntry entry, std::vector<double> damping)
      : ResolvedStratifiedQGForcing(entry), damping_(std::move(damping)) {}
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
    return context.adaptiveDamping(damping_);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + metadataDynamicBytes() + vectorBytes(damping_);
  }

private:
  std::vector<double> damping_;
};

class QGLinearBottomFriction final : public ResolvedStratifiedQGForcing {
public:
  QGLinearBottomFriction(WVFrozenForcingEntry entry, double rate)
      : ResolvedStratifiedQGForcing(entry), rate_(rate) {}
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
    return context.linearBottomFriction(rate_);
  }

private:
  double rate_ = 0.0;
};

class QGQuadraticBottomFriction final : public ResolvedStratifiedQGForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  QGQuadraticBottomFriction(WVFrozenForcingEntry entry, double drag)
      : ResolvedStratifiedQGForcing(entry), drag_(drag) {}
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
    return context.quadraticBottomFriction(drag_);
  }

private:
  double drag_ = 0.0;
};

class QGBetaPlanePVAdvection final : public ResolvedStratifiedQGForcing {
public:
  bool requiresDiagnosticPhysicalFields() const noexcept override { return true; }
  QGBetaPlanePVAdvection(WVFrozenForcingEntry entry, double beta)
      : ResolvedStratifiedQGForcing(entry), beta_(beta) {}
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
    return context.betaPlanePVAdvection(beta_);
  }

private:
  double beta_ = 0.0;
};

class QGVerticalDiffusivity final : public ResolvedStratifiedQGForcing {
public:
  QGVerticalDiffusivity(const WVFrozenForcingEntry &entry,double kappaZ):ResolvedStratifiedQGForcing(entry),kappaZ_(kappaZ) {}
  WVKernelStatus addRightHandSide(WVStratifiedQGForcingExecutionContext &context) const override { return context.verticalDiffusivity(kappaZ_); }
private:
  double kappaZ_;
};

class QGExplicitAntialiasing final : public ResolvedStratifiedQGForcing {
public:
  QGExplicitAntialiasing(const WVFrozenForcingEntry &entry, std::vector<std::size_t> indices)
      : ResolvedStratifiedQGForcing(entry), indices_(std::move(indices)) {}
  WVKernelStatus addRightHandSide(WVStratifiedQGForcingExecutionContext &context) const override {
    context.filterTendency(indices_);
    return WVKernelStatus::ok();
  }
  std::size_t persistentBytes() const noexcept override { return sizeof(*this)+metadataDynamicBytes()+vectorBytes(indices_); }
private:
  std::vector<std::size_t> indices_;
};

class QGFixedAmplitude final : public ResolvedStratifiedQGForcing {
public:
  QGFixedAmplitude(WVFrozenForcingEntry entry,
                   WVStratifiedQGFixedAmplitudeConfiguration values)
      : ResolvedStratifiedQGForcing(entry), values_(std::move(values)) {}
  WVKernelStatus addRightHandSide(
      WVStratifiedQGForcingExecutionContext &context) const override {
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
  WVStratifiedQGFixedAmplitudeConfiguration values_;
};

WVKernelStatus decodeFixedAmplitude(
    const WVFrozenForcingEntry &entry, std::size_t coefficientCount,
    WVStratifiedQGFixedAmplitudeConfiguration *configuration) {
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
            "Stratified QG fixed-amplitude forcing accepts only A0."};
  const auto *indices = integerValues(entry.configuration, "A0Indices");
  const auto *real = realValues(entry.configuration, "A0ValuesReal");
  const auto *imag = realValues(entry.configuration, "A0ValuesImag");
  if (indices == nullptr && real == nullptr && imag == nullptr)
    return WVKernelStatus::ok();
  if (indices == nullptr || real == nullptr || imag == nullptr ||
      indices->size() != real->size() || real->size() != imag->size())
    return {WVKernelStatusCode::invalidConfiguration,
            "Stratified QG fixed-amplitude arrays must be complete and equal-length."};
  std::set<std::size_t> unique;
  for (std::size_t index = 0; index < indices->size(); ++index) {
    if ((*indices)[index] < 0 ||
        static_cast<std::size_t>((*indices)[index]) >= coefficientCount ||
        !std::isfinite((*real)[index]) || !std::isfinite((*imag)[index]))
      return {WVKernelStatusCode::invalidConfiguration,
              "Stratified QG fixed-amplitude state is outside compact A0 or nonfinite."};
    const auto converted = static_cast<std::size_t>((*indices)[index]);
    if (!unique.insert(converted).second)
      return {WVKernelStatusCode::invalidConfiguration,
              "Stratified QG fixed-amplitude A0 indices must be unique."};
    if (configuration != nullptr) {
      configuration->A0Indices.push_back(converted);
      configuration->A0Values.push_back(
          {(*real)[index], (*imag)[index]});
    }
  }
  return WVKernelStatus::ok();
}

std::vector<double> adaptiveDampingOperator(const WVStratifiedModalGeometry &g,WVKernelStatus &status) {
  double maximumComponent=0;
  for (std::size_t kl=0;kl<g.Nkl;++kl) maximumComponent=std::max(maximumComponent,std::max(std::abs(g.k[kl]),std::abs(g.l[kl])));
  if (!(maximumComponent>0)) { status={WVKernelStatusCode::invalidConfiguration,"Adaptive damping requires resolved horizontal wavenumbers."}; return {}; }
  const double jMax=*std::max_element(g.j.begin(),g.j.end());
  const auto jIndex=static_cast<std::size_t>(std::max_element(g.j.begin(),g.j.end())-g.j.begin());
  const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
  const double maxRadiusSquared=g.g*g.h_0[jIndex]/(f*f);
  const double delta=pi/maximumComponent, dk=std::min(2*pi/g.Lx,2*pi/g.Ly);
  const double cutoff=dk*std::pow(maximumComponent/dk,.75);
  const double prefactor=delta/(pi*pi), verticalPrefactor=maxRadiusSquared/delta;
  std::vector<double> damping(g.Nj*g.Nkl);
  for (std::size_t kl=0;kl<g.Nkl;++kl) for (std::size_t j=0;j<g.Nj;++j) {
    const double kh=std::hypot(g.k[kl],g.l[kl]);
    double qj=1;
    if (g.Nj>2) { const double dj=g.j[1]-g.j[0]; qj=vanishingFilter(g.j[j],dj*std::pow(jMax/dj,.75),jMax); }
    damping[j+g.Nj*kl]=-prefactor*vanishingFilter(kh,cutoff,maximumComponent)*kh*kh-verticalPrefactor*qj*f*f/(g.g*g.h_0[j]);
  }
  status=WVKernelStatus::ok(); return damping;
}

} // namespace

namespace detail {

WVKernelStatus preflightStratifiedQGExplicitAntialiasing(const WVFrozenForcingEntry &entry, std::size_t) {
  return preflightExplicitAntialiasing(entry,false);
}
WVKernelStatus createStratifiedQGExplicitAntialiasing(const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &descriptor, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  const auto &c = descriptor;
  auto status = preflightExplicitAntialiasing(entry,c.shouldAntialias);
  if (!status) return status;
  const double Nj = realValues(entry.configuration,"Nj")->front();
  const double maximumK = 2.0*pi*static_cast<double>(c.Nx/2)/c.Lx;
  std::vector<std::size_t> indices;
  for (std::size_t kl=0;kl<c.Nkl;++kl) for (std::size_t j=0;j<c.Nj;++j)
    if (std::hypot(c.k[kl],c.l[kl])>2.0*maximumK/3.0 || c.j[j]>Nj-1.0)
      indices.push_back(j+c.Nj*kl);
  forcing = std::make_unique<QGExplicitAntialiasing>(entry,std::move(indices));
  return WVKernelStatus::ok();
}

WVKernelStatus preflightStratifiedQGEmptyForcing(
    const WVFrozenForcingEntry &entry, std::size_t) {
  return emptyConfiguration(entry)
             ? WVKernelStatus::ok()
             : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                              "This Stratified QG forcing accepts no configuration values."};
}

WVKernelStatus preflightStratifiedQGFixedAmplitude(
    const WVFrozenForcingEntry &entry, std::size_t coefficientCount) {
  return decodeFixedAmplitude(entry, coefficientCount, nullptr);
}

WVKernelStatus preflightStratifiedQGScalarForcing(
    const WVFrozenForcingEntry &entry, std::size_t) {
  if (entry.configuration.values.size() != 1 ||
      !entry.configuration.values.front().dimensions.empty())
    return {WVKernelStatusCode::invalidConfiguration,
            "Stratified QG scalar forcing requires exactly one scalar value."};
  const auto *values = std::get_if<std::vector<double>>(
      &entry.configuration.values.front().storage);
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()))
    return {WVKernelStatusCode::invalidConfiguration,
            "Stratified QG scalar forcing requires one finite real value."};
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGNonlinearAdvection(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  forcing = std::make_unique<QGNonlinearAdvection>(entry);
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGAdaptiveDamping(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &descriptor, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  WVKernelStatus status;
  auto damping = adaptiveDampingOperator(descriptor, status);
  if (!status)
    return status;
  forcing =
      std::make_unique<QGAdaptiveDamping>(entry, std::move(damping));
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGFixedAmplitude(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &descriptor, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  WVStratifiedQGFixedAmplitudeConfiguration values;
  auto status = decodeFixedAmplitude(entry, descriptor.Nj*descriptor.Nkl, &values);
  if (!status)
    return status;
  forcing =
      std::make_unique<QGFixedAmplitude>(entry, std::move(values));
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGQuadraticBottomFriction(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  const auto *values = realValues(entry.configuration, "Cd");
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()) || values->front() < 0.0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic quadratic drag requires one finite nonnegative Cd."};
  forcing = std::make_unique<QGQuadraticBottomFriction>(
      entry, values->front());
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGLinearBottomFriction(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  const auto *values = realValues(entry.configuration, "r");
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()) || values->front() < 0.0)
    return {WVKernelStatusCode::invalidConfiguration,
            "Barotropic linear drag requires one finite nonnegative r."};
  forcing =
      std::make_unique<QGLinearBottomFriction>(entry, values->front());
  return WVKernelStatus::ok();
}

WVKernelStatus createStratifiedQGBetaPlanePVAdvection(
    const WVFrozenForcingEntry &entry,
    const WVStratifiedModalGeometry &descriptor, bool,
    std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  const auto &configuration = descriptor;
  const double beta = 2.0 * configuration.rotationRate *
                      std::cos(configuration.latitude * pi / 180.0) /
                      configuration.planetaryRadius;
  forcing = std::make_unique<QGBetaPlanePVAdvection>(entry, beta);
  return WVKernelStatus::ok();
}

WVKernelStatus preflightStratifiedQGVerticalDiffusivity(const WVFrozenForcingEntry &entry,std::size_t) {
  return preflightVerticalDiffusivity(entry,false);
}
WVKernelStatus createStratifiedQGVerticalDiffusivity(const WVFrozenForcingEntry &entry,const WVStratifiedModalGeometry &,bool,std::unique_ptr<WVStratifiedQGForcing> &forcing) {
  forcing=std::make_unique<QGVerticalDiffusivity>(entry,realValues(entry.configuration,"kappa_z")->front());
  return WVKernelStatus::ok();
}

} // namespace detail

void WVStratifiedQGForcingExecutionContext::filterTendency(const std::vector<std::size_t> &indices) {
  if (!outputInitialized_) {
    engine_->initializeOutputWithZeros(F0_);
    outputInitialized_ = true;
  }
  for (const auto index : indices) F0_.data[index] = {};
}

WVKernelStatus WVStratifiedQGForcingExecutionContext::accumulate(WVKernelStatus status) {
  if (!status) return status;
  if (!outputInitialized_) { engine_->initializeOutputWithZeros(F0_); outputInitialized_=true; }
  for (std::size_t i=0;i<engine_->tendencyScratch_.size();++i) {
    F0_.data[i].real+=engine_->tendencyScratch_[i].real;
    F0_.data[i].imag+=engine_->tendencyScratch_[i].imag;
  }
  return WVKernelStatus::ok();
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::nonlinearAdvection() {
  if (engine_->diagnosticWorkspace_) {
    auto& work=*engine_->diagnosticWorkspace_;
    WVRealVolumeView raw{work.raw.data(),engine_->kernel().spatialShape()};
    WVRealFieldBundleConstView fields;
    auto prepared=engine_->diagnosticVelocity(A0_,fields); if (!prepared) return prepared;
    const auto result=engine_->kernel().nonlinearFlux(A0_,{engine_->tendencyScratch_.data(),F0_.shape},0,&raw,&fields);
    work.spatialCaptured=bool(result);
    if (result) engine_->metrics_.physicalFieldReconstructionCount+=2;
    return result;
  }
  if (engine_->evaluationPolicy_ == WVVariableEvaluationPolicy::lowMemory) {
    WVRealFieldBundleConstView fields;
    auto status = engine_->evaluationVelocity(A0_, fields);
    if (!status)
      return status;
    status = engine_->kernel().nonlinearFlux(
        A0_, {engine_->tendencyScratch_.data(), F0_.shape}, 0, nullptr,
        &fields);
    if (status) {
      engine_->metrics_.physicalFieldReconstructionCount += 2;
      engine_->metrics_.spatialTendencyProjectionCount += 1;
    }
    (void)engine_->evaluation_.evict(
        {WVVariableEvaluationNode::physicalField, 0});
    (void)engine_->evaluation_.evict(
        {WVVariableEvaluationNode::physicalField, 1});
    return accumulate(status);
  }
  const WVVariableEvaluationKey nonlinear{
      WVVariableEvaluationNode::forcingTendency, 0};
  const bool ready = engine_->evaluation_.ready(nonlinear);
  auto status = engine_->evaluation_.evaluate(
      nonlinear, engine_->nonlinearScratch_.size() * sizeof(WVComplex64), [&] {
        WVRealFieldBundleConstView fields;
        auto prepared = engine_->evaluationVelocity(A0_, fields);
        if (!prepared)
          return prepared;
        return engine_->kernel().nonlinearFlux(
            A0_, {engine_->nonlinearScratch_.data(), F0_.shape}, 0, nullptr,
            &fields);
      });
  if (!status)
    return status;
  if (!ready) {
    engine_->metrics_.physicalFieldReconstructionCount += 2;
    engine_->metrics_.spatialTendencyProjectionCount += 1;
  }
  std::copy(engine_->nonlinearScratch_.begin(),
            engine_->nonlinearScratch_.end(), engine_->tendencyScratch_.begin());
  return accumulate(WVKernelStatus::ok());
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::adaptiveDamping(const std::vector<double> &damping) {
  double speed=0;
  if (engine_->diagnosticWorkspace_) {
    auto status=engine_->diagnosticWorkspace_->evaluateHorizontalMaximum(
        speed,[&](double& maximum) {
          return engine_->computeHorizontalSpeedMaximum(A0_,maximum);
        });
    if (!status) return status;
  } else {
    if (engine_->evaluationPolicy_ == WVVariableEvaluationPolicy::lowMemory) {
      auto status=engine_->horizontalSpeedMaximum(A0_,speed);
      if (!status) return status;
      (void)engine_->evaluation_.evict(
          {WVVariableEvaluationNode::reduction, 0});
      (void)engine_->evaluation_.evict(
          {WVVariableEvaluationNode::physicalField, 0});
      (void)engine_->evaluation_.evict(
          {WVVariableEvaluationNode::physicalField, 1});
    } else {
      auto status=engine_->horizontalSpeedMaximum(A0_,speed); if (!status) return status;
    }
  }
  if (!outputInitialized_) { engine_->initializeOutputWithZeros(F0_); outputInitialized_=true; }
  for (std::size_t i=0;i<damping.size();++i) {
    F0_.data[i].real+=speed*damping[i]*A0_.data[i].real;
    F0_.data[i].imag+=speed*damping[i]*A0_.data[i].imag;
  }
  return WVKernelStatus::ok();
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::linearBottomFriction(double rate) {
  if (engine_->diagnosticWorkspace_) {
    auto& work=*engine_->diagnosticWorkspace_;
    WVRealVolumeView raw{work.raw.data(),engine_->kernel().spatialShape()};
    const auto result=engine_->kernel().linearBottomFrictionFlux(A0_,rate,{engine_->tendencyScratch_.data(),F0_.shape},&raw);
    work.spatialCaptured=bool(result);
    if (result) ++engine_->metrics_.physicalFieldReconstructionCount;
    return result;
  }
  const auto status = engine_->kernel().linearBottomFrictionFlux(A0_,rate,{engine_->tendencyScratch_.data(),F0_.shape});
  if (status) {
    engine_->metrics_.physicalFieldReconstructionCount += 1;
    engine_->metrics_.spatialTendencyProjectionCount += 1;
  }
  return accumulate(status);
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::quadraticBottomFriction(double drag) {
  if (engine_->diagnosticWorkspace_) {
    auto& work=*engine_->diagnosticWorkspace_;
    WVRealVolumeView raw{work.raw.data(),engine_->kernel().spatialShape()};
    WVRealFieldBundleConstView fields;
    auto prepared=engine_->diagnosticVelocity(A0_,fields); if (!prepared) return prepared;
    const auto result=engine_->kernel().quadraticBottomFrictionFlux(A0_,drag,{engine_->tendencyScratch_.data(),F0_.shape},&raw,&fields);
    work.spatialCaptured=bool(result);
    return result;
  }
  if (engine_->evaluationPolicy_ == WVVariableEvaluationPolicy::lowMemory) {
    WVRealFieldBundleConstView fields;
    auto status = engine_->evaluationVelocity(A0_, fields);
    if (!status)
      return status;
    status = engine_->kernel().quadraticBottomFrictionFlux(
        A0_, drag, {engine_->tendencyScratch_.data(), F0_.shape}, nullptr,
        &fields);
    if (status) {
      engine_->metrics_.spatialTendencyProjectionCount += 2;
    }
    (void)engine_->evaluation_.evict(
        {WVVariableEvaluationNode::physicalField, 0});
    (void)engine_->evaluation_.evict(
        {WVVariableEvaluationNode::physicalField, 1});
    return accumulate(status);
  }
  WVRealFieldBundleConstView fields;
  auto status = engine_->evaluationVelocity(A0_, fields);
  if (!status)
    return status;
  status = engine_->kernel().quadraticBottomFrictionFlux(
      A0_, drag, {engine_->tendencyScratch_.data(), F0_.shape}, nullptr,
      &fields);
  if (status) {
    engine_->metrics_.spatialTendencyProjectionCount += 2;
  }
  return accumulate(status);
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::verticalDiffusivity(double kappaZ) {
  if (engine_->diagnosticWorkspace_) {
    auto& work=*engine_->diagnosticWorkspace_;
    WVRealVolumeView raw{work.raw.data(),engine_->kernel().spatialShape()};
    const auto result=engine_->kernel().verticalDiffusivityFlux(A0_,kappaZ,{engine_->tendencyScratch_.data(),F0_.shape},&raw);
    work.spatialCaptured=bool(result);
    return result;
  }
  return accumulate(engine_->kernel().verticalDiffusivityFlux(A0_,kappaZ,{engine_->tendencyScratch_.data(),F0_.shape}));
}
WVKernelStatus WVStratifiedQGForcingExecutionContext::betaPlanePVAdvection(double beta) {
  if (engine_->diagnosticWorkspace_) {
    WVRealFieldBundleConstView fields;
    auto status=engine_->diagnosticVelocity(A0_,fields); if (!status) return status;
    auto& work=*engine_->diagnosticWorkspace_;
    const auto R=engine_->kernel().spatialShape().elementCount();
    for (std::size_t i=0;i<R;++i) work.raw[i]=-beta*fields.data[R+i];
    work.spatialCaptured=true;
    return WVKernelStatus::ok();
  }
  return accumulate(engine_->kernel().linearFlux(A0_,{engine_->tendencyScratch_.data(),F0_.shape},beta));
}

void WVStratifiedQGForcingExecutionContext::zeroSelectedTendencies(
    const WVStratifiedQGFixedAmplitudeConfiguration &configuration) {
  if (!outputInitialized_) {
    engine_->initializeOutputWithZeros(F0_);
    outputInitialized_ = true;
  }
  for (const auto index : configuration.A0Indices)
    F0_.data[index] = {};
}

WVStratifiedQGForcingEngine::~WVStratifiedQGForcingEngine() = default;

const WVForcingEvaluationDependencies*
WVStratifiedQGForcingEngine::forcingEvaluationDependencies(
    std::size_t index) const noexcept {
  const auto* forcing=forcingInstance(index);
  if(!forcing || !catalog_) return nullptr;
  const auto* registration=catalog_->forcings().registration(
      forcing->typeIdentifier(),forcing->contractVersion());
  return registration ? &registration->evaluationDependencies : nullptr;
}

WVKernelStatus WVStratifiedQGForcingEngine::validateSchedule(
    const WVStratifiedModalGeometry &configuration,
    const WVFrozenForcingSchedule &schedule, std::size_t coefficientCount,
    const WVExtensionCatalog &catalog) {
  if (schedule.profileIdentifier != WVForcingScheduleProfileIdentifier ||
      schedule.profileVersion != WVForcingScheduleProfileVersion)
    return {WVKernelStatusCode::unsupportedOperation,
            "Unsupported frozen forcing schedule profile."};
  if (coefficientCount == 0)
    return {WVKernelStatusCode::invalidShape,
            "Stratified QG forcing requires nonempty compact A0."};
  std::set<std::string> names;
  for (const auto &entry : schedule.entries) {
    const auto *registration = catalog.forcings().registration(
        entry.typeIdentifier, entry.contractVersion);
    if (registration == nullptr || !registration->isSupported ||
        !registration->stratifiedQGFactory ||
        registration->contractVersion != entry.contractVersion)
      return {WVKernelStatusCode::unsupportedOperation,
              "The frozen schedule has no matching Stratified QG forcing implementation."};
    if (entry.stage != registration->stratifiedQGStage ||
        !containsForcingType(*registration, qgForcingType(entry.stage)))
      return {WVKernelStatusCode::invalidConfiguration,
              "A Stratified QG forcing record is assigned to an incompatible stage."};
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
    if (!registration->stratifiedQGPreflight)
      return {WVKernelStatusCode::invalidConfiguration,
              "A Stratified QG forcing registration lacks allocation-light preflight."};
    status = registration->stratifiedQGPreflight(entry, coefficientCount);
    if (!status)
      return status;
  }
  (void)configuration;
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::create(
    std::shared_ptr<const WVStratifiedModalSource> source,
    const WVFrozenForcingSchedule &schedule,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> fftEngine,
    std::unique_ptr<WVStratifiedQGForcingEngine> &forcingEngine,
    const WVVariableKernelServices &services) {
  if (!catalog)
    return {WVKernelStatusCode::invalidPointer,
            "Stratified QG forcing construction requires an extension catalog."};
  if (!fftEngine)
    return {WVKernelStatusCode::invalidPointer,
            "Stratified QG forcing construction requires an FFT engine."};
  if (!source) return {WVKernelStatusCode::invalidPointer,"Stratified QG requires an immutable scientific source."};
  const auto& configuration=source->geometry();
  if (!configuration.Nj || configuration.Nkl>std::numeric_limits<std::size_t>::max()/configuration.Nj)
    return {WVKernelStatusCode::sizeOverflow,"Invalid Stratified QG coefficient count."};
  const auto coefficientCount=configuration.Nj*configuration.Nkl;
  auto status=validateSchedule(configuration,schedule,coefficientCount,*catalog);
  if (!status)
    return status;
  try {
    auto candidate = std::unique_ptr<WVStratifiedQGForcingEngine>(
        new WVStratifiedQGForcingEngine());
    candidate->catalog_ = std::move(catalog);
    status = WVTransformStratifiedQGKernel::create(
        std::move(source), std::move(fftEngine), candidate->kernel_,
        services.matrixBackendFactory, services.execution);
    if (!status)
      return status;
    candidate->tendencyScratch_.resize(coefficientCount);
    status = candidate->initialize(schedule);
    if (!status)
      return status;
    forcingEngine = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Stratified QG forcing-engine allocation failed."};
  } catch (const std::exception &error) {
    return {WVKernelStatusCode::invalidConfiguration, error.what()};
  }
}

WVKernelStatus WVStratifiedQGForcingEngine::initialize(
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
    std::unique_ptr<WVStratifiedQGForcing> resolved;
    auto status = catalog_->forcings().createStratifiedQG(
        *entry, kernel_->geometry(), hasAdaptiveDamping, resolved);
    if (!status)
      return status;
    if (!resolved || resolved->stage() != entry->stage)
      return {WVKernelStatusCode::invalidConfiguration,
              "A Stratified QG forcing factory returned an incompatible implementation."};
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
  const auto spatialCount = kernel_->spatialShape().elementCount();
  velocityScratch_.resize(3 * spatialCount);
  nonlinearScratch_.resize(tendencyScratch_.size());
  auto status = evaluation_.prepare(
      {{WVVariableEvaluationNode::physicalField, 0},
       {WVVariableEvaluationNode::physicalField, 1},
       {WVVariableEvaluationNode::reduction, 0},
       {WVVariableEvaluationNode::forcingTendency, 0}});
  if (!status)
    return status;
  metrics_.scheduleBytes =
      scheduleIdentifier_.capacity() +
      forcing_.capacity() * sizeof(std::unique_ptr<WVStratifiedQGForcing>);
  metrics_.workspaceCapacityBytes = vectorBytes(tendencyScratch_) +
                                    vectorBytes(nonlinearScratch_) +
                                    vectorBytes(velocityScratch_);
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::beginStateEvaluation(
    const WVComplexConstView &A0) {
  if (evaluation_.active() || executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Stratified QG state evaluation is already active."};
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
  evaluationVelocity_ = {};
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::validateVariableEvaluationPolicyChange(
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

WVKernelStatus WVStratifiedQGForcingEngine::setVariableEvaluationPolicy(
    WVVariableEvaluationPolicy policy) {
  auto status=validateVariableEvaluationPolicyChange(policy);
  if(!status) return status;
  if (policy == WVVariableEvaluationPolicy::reuse &&
      (velocityScratch_.size() !=
           3 * kernel_->spatialShape().elementCount() ||
       nonlinearScratch_.size() != tendencyScratch_.size())) {
    std::vector<double> preparedVelocity;
    std::vector<WVComplex64> preparedNonlinear;
    try {
      if(velocityScratch_.size()!=3*kernel_->spatialShape().elementCount())
        preparedVelocity.resize(3*kernel_->spatialShape().elementCount());
      if(nonlinearScratch_.size()!=tendencyScratch_.size())
        preparedNonlinear.resize(tendencyScratch_.size());
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate the Stratified QG velocity cache."};
    }
    if(!preparedVelocity.empty()) velocityScratch_.swap(preparedVelocity);
    if(!preparedNonlinear.empty()) nonlinearScratch_.swap(preparedNonlinear);
  } else if (policy == WVVariableEvaluationPolicy::lowMemory) {
    std::vector<WVComplex64>().swap(nonlinearScratch_);
  }
  evaluationPolicy_ = policy;
  metrics_.workspaceCapacityBytes = vectorBytes(tendencyScratch_) +
                                    vectorBytes(nonlinearScratch_) +
                                    vectorBytes(velocityScratch_);
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::endStateEvaluation() {
  if (!evaluation_.active())
    return {WVKernelStatusCode::invalidConfiguration,
            "Stratified QG state evaluation is not active."};
  evaluation_.end();
  evaluationState_ = {};
  evaluationVelocity_ = {};
  const bool ownsKernelScope=evaluationOwnsKernelScope_;
  evaluationOwnsKernelScope_=false;
  return ownsKernelScope ? kernel_->endStateEvaluation() : WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::validateStateEvaluation(
    const WVComplexConstView &A0) const noexcept {
  if (!evaluation_.active() || A0.data != evaluationState_.data ||
      A0.shape.rows != evaluationState_.shape.rows ||
      A0.shape.columns != evaluationState_.shape.columns)
    return {WVKernelStatusCode::invalidConfiguration,
            "Stratified QG coefficients do not belong to the active evaluation."};
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::evaluationVelocity(
    WVComplexConstView A0, WVRealFieldBundleConstView &fields) {
  auto status = validateStateEvaluation(A0);
  if (!status)
    return status;
  const auto shape = kernel_->spatialShape();
  const auto count = shape.elementCount();
  const auto required = 3 * count;
  if (velocityScratch_.size() != required) {
    try {
      velocityScratch_.resize(required);
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Unable to allocate the Stratified QG velocity cache."};
    }
    metrics_.workspaceCapacityBytes = vectorBytes(tendencyScratch_) +
                                      vectorBytes(nonlinearScratch_) +
                                      vectorBytes(velocityScratch_);
  }
  const WVStratifiedQGField names[] = {WVStratifiedQGField::u,
                                      WVStratifiedQGField::v};
  for (std::size_t channel = 0; channel < 2; ++channel) {
    const WVVariableEvaluationKey key{WVVariableEvaluationNode::physicalField,
                                      static_cast<std::uint32_t>(channel)};
    const bool ready = evaluation_.ready(key);
    status = evaluation_.evaluate(key, count * sizeof(double), [&] {
      return kernel_->transformA0ToField(
          A0, names[channel],
          {velocityScratch_.data() + channel * count, shape});
    });
    if (!status)
      return status;
    if (ready)
      ++metrics_.physicalFieldReuseCount;
    else
      ++metrics_.physicalFieldReconstructionCount;
  }
  evaluationVelocity_ = {
      velocityScratch_.data(),
      {shape.first, shape.second, shape.third, 2}};
  fields = evaluationVelocity_;
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::horizontalSpeedMaximum(
    WVComplexConstView A0, double &maximum) {
  ScopedStratifiedQGEvaluation evaluation(*this, A0);
  if (!evaluation.status())
    return evaluation.status();
  auto status = evaluation_.evaluate(
      {WVVariableEvaluationNode::reduction, 0}, sizeof(double), [&] {
        return computeHorizontalSpeedMaximum(A0,horizontalSpeedMaximum_);
      });
  if (status)
    maximum = horizontalSpeedMaximum_;
  return status;
}

WVKernelStatus WVStratifiedQGForcingEngine::computeHorizontalSpeedMaximum(
    WVComplexConstView A0,double& maximum) {
  WVRealFieldBundleConstView fields;
  auto status=diagnosticWorkspace_ ? diagnosticVelocity(A0,fields) :
      evaluationVelocity(A0,fields);
  if (!status) return status;
  const auto count=kernel_->spatialShape().elementCount();
  maximum=0.0;
  for(std::size_t index=0;index<count;++index)
    maximum=std::max(maximum,
        std::hypot(fields.data[index],fields.data[count+index]));
  ++metrics_.horizontalSpeedMaximumReductionCount;
  return WVKernelStatus::ok();
}

void WVStratifiedQGForcingEngine::initializeOutputWithZeros(
    WVComplexView &F0) {
  std::fill_n(F0.data, kernel_->geometry().Nj*kernel_->geometry().Nkl, WVComplex64{});
}

WVKernelStatus WVStratifiedQGForcingEngine::evaluateRightHandSide(
    const WVComplexConstView &A0, WVComplexView &F0,
    WVRealFieldBundleConstView *advectionFields) {
  if (advectionFields != nullptr)
    *advectionFields = {};
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Stratified QG forcing-engine execution is not reentrant."};
  const auto expected = kernel_->spectralShape();
  if (A0.shape.rows != expected.rows || A0.shape.columns != expected.columns ||
      F0.shape.rows != expected.rows || F0.shape.columns != expected.columns)
    return {WVKernelStatusCode::invalidShape,
            "Stratified QG forcing requires compact [Nj,Nkl] A0 and F0."};
  if (A0.data == nullptr || F0.data == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Stratified QG forcing received null compact state storage."};
  const auto inputAddress = reinterpret_cast<std::uintptr_t>(A0.data);
  const auto outputAddress = reinterpret_cast<std::uintptr_t>(F0.data);
  const auto bytes = expected.elementCount()*sizeof(WVComplex64);
  if (inputAddress > std::numeric_limits<std::uintptr_t>::max()-bytes || outputAddress > std::numeric_limits<std::uintptr_t>::max()-bytes)
    return {WVKernelStatusCode::invalidPointer,"Stratified QG state storage exceeds the address range."};
  if (inputAddress < outputAddress+bytes && outputAddress < inputAddress+bytes)
    return {WVKernelStatusCode::overlappingArrays,
            "Stratified QG A0 and F0 must not overlap."};
  ScopedStratifiedQGEvaluation evaluation(*this, A0);
  if (!evaluation.status())
    return evaluation.status();
  executing_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{executing_};
  WVStratifiedQGForcingExecutionContext context;
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
    WVRealFieldBundleConstView velocity;
    const auto status = evaluationVelocity(A0, velocity);
    if (!status)
      return status;
    *advectionFields = velocity;
    advectionFields->shape.fourth = 3;
  }
  ++metrics_.evaluationCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedQGForcingEngine::diagnosticVelocity(WVComplexConstView a,WVRealFieldBundleConstView& fields) {
  auto& work=*diagnosticWorkspace_;
  const auto shape=kernel().spatialShape(); const auto R=shape.elementCount();
  if (!work.physicalPrepared) {
    auto status=kernel().transformA0ToField(a,WVStratifiedQGField::u,{work.physical.data(),shape}); if (!status) return status;
    status=kernel().transformA0ToField(a,WVStratifiedQGField::v,{work.physical.data()+R,shape}); if (!status) return status;
    work.physicalPrepared=true; metrics_.physicalFieldReconstructionCount+=2;
  } else ++metrics_.physicalFieldReuseCount;
  fields={work.physical.data(),{shape.first,shape.second,shape.third,2}};
  return WVKernelStatus::ok();
}
WVKernelStatus WVStratifiedQGForcingEngine::evaluateForcingTendencies(
    const WVComplexConstView& A0,const WVForcingTendencyOutput* outputs,
    std::size_t count,const WVRealFieldBundleConstView* preparedPhysical,
    detail::WVForcingDiagnosticWorkspace* session) {
  if (executing_) return {WVKernelStatusCode::reentrantExecution,"Forcing diagnostics require an idle engine."};
  tendencyMetrics_.workspaceLastPeakBytes=0;
  const auto spectral=kernel().spectralShape(); const auto volume=kernel().spatialShape();
  const WVShape4D spatial{volume.first,volume.second,volume.third,1};
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
  ScopedStratifiedQGEvaluation evaluation(*this, A0);
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
    const auto R=volume.elementCount();
    if (work.spectral.rows!=spectral.rows || work.spectral.columns!=spectral.columns ||
        work.spatial.first!=spatial.first || work.spatial.second!=spatial.second ||
        work.spatial.third!=spatial.third || work.spatial.fourth!=spatial.fourth ||
        work.flux.size()!=S || work.previous.size()!=S || work.temporary.size()!=S ||
        work.cumulative.size()!=spatial.elementCount() || work.raw.size()!=spatial.elementCount() ||
        work.physical.size()!=2*R)
      return {WVKernelStatusCode::invalidShape,
              "Forcing diagnostic session has incompatible Stratified QG storage."};
    if (!work.initialized()) {
      if (preparedPhysical) {
        std::copy_n(preparedPhysical->data,preparedPhysical->shape.elementCount(),
                    work.physical.data());
        work.physicalPrepared=true;
      }
      work.markInitialized();
    }
    WVStratifiedQGForcingExecutionContext context;
    context.engine_=this; context.A0_=A0; context.outputInitialized_=true;
    executing_=true; diagnosticWorkspace_=&work;
    struct Guard {
      WVStratifiedQGForcingEngine& engine;
      ~Guard() {
        engine.tendencyMetrics_.workspaceLastPeakBytes=engine.diagnosticWorkspace_->bytes();
        engine.tendencyMetrics_.workspaceHighWaterBytes=std::max(engine.tendencyMetrics_.workspaceHighWaterBytes,engine.diagnosticWorkspace_->bytes());
        engine.tendencyMetrics_.workspaceLiveBytes=0;
        engine.diagnosticWorkspace_=nullptr; engine.executing_=false;
      }
    } guard{*this};
    tendencyMetrics_.workspaceLiveBytes=work.bytes();
    return detail::evaluateForcingTendencySequence(forcing_,work,outputs,count,tendencyMetrics_,
      [&](const WVStratifiedQGForcing& forcing,WVFlux& destination) {
        context.F0_=destination.F0;
        return forcing.addRightHandSide(context);
      },
      [&](WVRealFieldBundleConstView fields,WVFlux& destination) {
        return kernel().transformQGPVToA0({fields.data,volume},destination.F0);
      },
      [&](const std::vector<WVComplex64>& difference,WVRealFieldBundleView destination) {
        return kernel().transformSpectralTendencyToSpatial({difference.data(),spectral},{destination.data,volume});
      });
  } catch (const std::bad_alloc&) {
    return {WVKernelStatusCode::allocationFailure,"Unable to allocate event forcing diagnostic workspace."};
  }
}

WVStateConstraintResult
WVStratifiedQGForcingEngine::restoreForcingAmplitudes(WVComplexView &A0) {
  if (evaluation_.active() || executing_)
    return {{WVKernelStatusCode::reentrantExecution,
             "Stratified QG constraints cannot mutate an active evaluation."},
            0, false};
  const auto expected = kernel_->spectralShape();
  if (A0.shape.rows != expected.rows || A0.shape.columns != expected.columns)
    return {{WVKernelStatusCode::invalidShape,
             "Stratified QG constraints require compact [Nj,Nkl] A0."},
            0, false};
  if (A0.data == nullptr)
    return {{WVKernelStatusCode::invalidPointer,
             "Stratified QG constraint A0 storage is null."},
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

std::size_t WVStratifiedQGForcingEngine::persistentBytes() const noexcept {
  return sizeof(*this) +
         (kernel_ == nullptr ? 0 : kernel_->persistentBytes()) +
         metrics_.scheduleBytes + metrics_.derivedOperatorBytes +
         metrics_.workspaceCapacityBytes + evaluation_.persistentBytes();
}

} // namespace wavevortex::runtime
