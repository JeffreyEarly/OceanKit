#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"

#include "WVReferenceFFTEngine.hpp"
#include "WVStratifiedModalTestFixture.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

void requireClose(double actual, double expected, const std::string &message) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  require(std::abs(actual - expected) <= 2.0e-12 * scale, message);
}

void requireWithinOneE12(double actual, double expected,
                         const std::string &message) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  require(std::abs(actual - expected) <= 1.0e-12 * scale, message);
}

WVTransformConstantStratificationConfiguration
configuration(std::size_t nx, std::size_t ny, bool hydrostatic,
              bool antialias) {
  WVTransformConstantStratificationConfiguration value;
  value.Nx = nx;
  value.Ny = ny;
  value.Nz = 7;
  value.Nj = 6;
  value.Lx = 15000.0;
  value.Ly = 12000.0;
  value.Lz = 1300.0;
  value.N0 = 5.2e-3;
  value.rho0 = 1025.0;
  value.g = 9.81;
  value.planetaryRadius = 6.371e6;
  value.rotationRate = 7.2921e-5;
  value.latitude = 33.0;
  value.isHydrostatic = hydrostatic;
  value.shouldAntialias = antialias;
  return value;
}

struct OwnedState {
  WVShape2D shape;
  std::vector<WVComplex64> Ap;
  std::vector<WVComplex64> Am;
  std::vector<WVComplex64> A0;

  WVState view(double t = 37.25, double t0 = -3.5) const noexcept {
    return {t,
            t0,
            {{Ap.data(), shape}, {Am.data(), shape}, {A0.data(), shape}}};
  }
};

OwnedState stateFor(
    const WVTransformConstantStratificationConfiguration &configuration) {
  WVTransformConstantStratificationDescriptor descriptor;
  const auto status = WVTransformConstantStratificationDescriptor::create(
      configuration, descriptor);
  require(static_cast<bool>(status), "unable to create state descriptor");
  OwnedState state;
  state.shape = descriptor.spectralShape();
  const auto count = state.shape.elementCount();
  state.Ap.resize(count);
  state.Am.resize(count);
  state.A0.resize(count);
  for (std::size_t index = 0; index < count; ++index) {
    const double value = static_cast<double>(index + 1);
    state.Ap[index] = {2.3e-3 * std::sin(0.37 * value),
                       -1.7e-3 * std::cos(0.19 * value)};
    state.Am[index] = {-1.1e-3 * std::cos(0.23 * value),
                       1.9e-3 * std::sin(0.41 * value)};
    state.A0[index] = {1.3e-3 * std::sin(0.29 * value),
                       0.7e-3 * std::cos(0.31 * value)};
  }
  return state;
}

WVFieldRequest full(std::string name) {
  return {"full_" + name, std::move(name), {}};
}

std::size_t outputIndex(const WVFieldEvaluationPlan &plan,
                        const std::string &identifier) {
  const auto &outputs = plan.outputs();
  const auto iterator =
      std::find_if(outputs.begin(), outputs.end(), [&](const auto &output) {
        return output.identifier == identifier;
      });
  require(iterator != outputs.end(), "missing output " + identifier);
  return static_cast<std::size_t>(iterator - outputs.begin());
}

class WrappedPlan final : public WVFFTPlan {
public:
  WrappedPlan(std::unique_ptr<WVFFTPlan> plan,
              std::shared_ptr<std::size_t> active)
      : plan_(std::move(plan)), active_(std::move(active)) {
    ++*active_;
  }
  ~WrappedPlan() override { --*active_; }
  WVKernelStatus execute(const void *input, void *output) override {
    return plan_->execute(input, output);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + plan_->persistentBytes();
  }

private:
  std::unique_ptr<WVFFTPlan> plan_;
  std::shared_ptr<std::size_t> active_;
};

class CountingEngine final : public WVFFTEngine {
public:
  explicit CountingEngine(std::shared_ptr<std::size_t> active)
      : active_(std::move(active)) {}
  std::string identifier() const override { return "counting-reference"; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + reference_.persistentBytes() - sizeof(reference_) +
           sizeof(std::size_t);
  }
  WVKernelStatus createPlan(const WVFFTPlanSpecification &specification,
                            std::unique_ptr<WVFFTPlan> &plan) override {
    std::unique_ptr<WVFFTPlan> referencePlan;
    auto status = reference_.createPlan(specification, referencePlan);
    if (!status)
      return status;
    plan = std::make_unique<WrappedPlan>(std::move(referencePlan), active_);
    return WVKernelStatus::ok();
  }

private:
  WVReferenceFFTEngine reference_;
  std::shared_ptr<std::size_t> active_;
};

struct ActivePlanProbe {
  WVFieldEvaluationService* service=nullptr;
  bool armed=false,attempted=false;
  WVKernelStatus status;
};
class ActivePlanProbePlan final : public WVFFTPlan {
public:
  ActivePlanProbePlan(std::unique_ptr<WVFFTPlan> plan,
      std::shared_ptr<ActivePlanProbe> probe)
      :plan_(std::move(plan)),probe_(std::move(probe)) {}
  WVKernelStatus execute(const void* input,void* output) override {
    if(probe_->armed && !probe_->attempted) {
      probe_->attempted=true;
      WVFieldEvaluationPlan plan;
      probe_->status=probe_->service->createPlanForActiveEvaluation(
          {full("v")},plan);
    }
    return plan_->execute(input,output);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this)+plan_->persistentBytes();
  }
private:
  std::unique_ptr<WVFFTPlan> plan_;
  std::shared_ptr<ActivePlanProbe> probe_;
};
class ActivePlanProbeEngine final : public WVFFTEngine {
public:
  explicit ActivePlanProbeEngine(std::shared_ptr<ActivePlanProbe> probe)
      :probe_(std::move(probe)) {}
  std::string identifier() const override {return "active-plan-probe";}
  std::size_t persistentBytes() const noexcept override {return sizeof(*this);}
  WVKernelStatus createPlan(const WVFFTPlanSpecification& specification,
      std::unique_ptr<WVFFTPlan>& plan) override {
    std::unique_ptr<WVFFTPlan> inner;
    auto status=reference_.createPlan(specification,inner);
    if(status) plan=std::make_unique<ActivePlanProbePlan>(
        std::move(inner),probe_);
    return status;
  }
private:
  WVReferenceFFTEngine reference_;
  std::shared_ptr<ActivePlanProbe> probe_;
};

class FailurePlan final : public WVFFTPlan {
public:
  WVKernelStatus execute(const void *, void *) override {
    return {WVKernelStatusCode::fftExecutionFailure,
            "injected field-evaluation FFT failure"};
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
};

class FailureEngine final : public WVFFTEngine {
public:
  explicit FailureEngine(WVKernelStatusCode creationCode)
      : creationCode_(creationCode) {}
  std::string identifier() const override { return "failure-injection"; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }
  WVKernelStatus createPlan(const WVFFTPlanSpecification &,
                            std::unique_ptr<WVFFTPlan> &plan) override {
    if (creationCode_ != WVKernelStatusCode::success)
      return {creationCode_, "injected field-evaluation planning failure"};
    plan = std::make_unique<FailurePlan>();
    return WVKernelStatus::ok();
  }

private:
  WVKernelStatusCode creationCode_;
};

struct FailOnceControl {
  bool armed = false;
};

class FailOncePlan final : public WVFFTPlan {
public:
  FailOncePlan(std::unique_ptr<WVFFTPlan> plan,
               std::shared_ptr<FailOnceControl> control)
      : plan_(std::move(plan)), control_(std::move(control)) {}
  WVKernelStatus execute(const void *input, void *output) override {
    if (control_->armed) {
      control_->armed = false;
      return {WVKernelStatusCode::fftExecutionFailure,
              "injected one-shot field-evaluation FFT failure"};
    }
    return plan_->execute(input, output);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + plan_->persistentBytes();
  }

private:
  std::unique_ptr<WVFFTPlan> plan_;
  std::shared_ptr<FailOnceControl> control_;
};

class FailOnceEngine final : public WVFFTEngine {
public:
  explicit FailOnceEngine(std::shared_ptr<FailOnceControl> control)
      : control_(std::move(control)) {}
  std::string identifier() const override { return "reference-fail-once"; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + reference_.persistentBytes() - sizeof(reference_);
  }
  WVKernelStatus createPlan(const WVFFTPlanSpecification &specification,
                            std::unique_ptr<WVFFTPlan> &plan) override {
    std::unique_ptr<WVFFTPlan> referencePlan;
    auto status = reference_.createPlan(specification, referencePlan);
    if (!status)
      return status;
    plan = std::make_unique<FailOncePlan>(std::move(referencePlan), control_);
    return WVKernelStatus::ok();
  }

private:
  std::shared_ptr<FailOnceControl> control_;
  WVReferenceFFTEngine reference_;
};

void verifyCatalog() {
  std::vector<std::string> expected = {
      "u",       "v",         "w",       "eta",    "pi",
      "p",       "psi",       "qgpv",    "rho_e",  "rho_total",
      "rho_bar", "zeta_x",    "zeta_y",  "zeta_z", "ssu",
      "ssv",     "ssh",       "energy",  "uvMax",  "wMax"};
  for(const auto& metadata:WVPortableVariableCatalog)
    if(metadata.ordinal>=23 && findExecutablePortableVariable(metadata.name)) expected.emplace_back(metadata.name);
  require(WVFieldEvaluationService::supportedFieldNames() == expected,
          "field catalog changed");
}

void verifyPlanValidation() {
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      configuration(6, 5, true, true),
      std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status), "validation service creation failed");
  WVFieldEvaluationPlan plan;
  status = service->createPlan({{"bad", "not_a_field", {}}}, plan);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "unknown field was accepted");
  status = service->createPlan({{"pending", "w_w", {}}}, plan);
  require(static_cast<bool>(status), "qualified diagnostic plan was rejected");
  status = service->createPlan({full("u"), {"full_u", "v", {}}}, plan);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "duplicate identifier was accepted");
  WVFieldSamplingRequest profiles;
  profiles.kind = WVFieldSamplingKind::fixedVerticalProfiles;
  profiles.xIndices = {0};
  profiles.yIndices = {1};
  status = service->createPlan({{"profile", "u", profiles}}, plan);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "zero-based profile index was accepted");
  WVFieldSamplingRequest positions;
  positions.kind = WVFieldSamplingKind::positions;
  positions.x = {0.0};
  positions.y = {0.0};
  positions.z = {0.0};
  status = service->createPlan({{"energy_at_position", "energy", positions}},
                               plan);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "scalar position sampling was accepted");

  auto equatorial = configuration(6, 5, true, true);
  equatorial.latitude = 0.0;
  status = WVFieldEvaluationService::create(
      equatorial, std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status), "equatorial service creation failed");
  status = service->createPlan({full("psi")}, plan);
  require(status.code == WVKernelStatusCode::unsupportedOperation,
          "undefined equatorial streamfunction was not rejected by planning");
}

void verifyPhaseDiagnostics(bool hydrostatic, bool antialias) {
  const auto config = configuration(8, 6, hydrostatic, antialias);
  auto state = stateFor(config);
  const auto original = state;
  WVTransformConstantStratificationDescriptor descriptor;
  require(bool(WVTransformConstantStratificationDescriptor::create(config, descriptor)),
          "phase descriptor creation failed");
  std::unique_ptr<WVFieldEvaluationService> service;
  require(bool(WVFieldEvaluationService::create(config,
      std::make_unique<WVReferenceFFTEngine>(), service)), "phase service creation failed");
  WVFieldEvaluationPlan plan;
  require(bool(service->createPlan({full("phase"), full("conjPhase"), full("Apt"),
      full("Amt"), full("A0t")}, plan)), "phase diagnostic plan rejected");
  const auto count = state.shape.elementCount();
  std::array<std::vector<WVComplex64>, 5> storage;
  std::array<WVFieldOutputView, 5> views{};
  for (std::size_t index = 0; index < storage.size(); ++index) {
    const auto &spec = plan.outputs()[index];
    require(spec.isComplex && spec.dimensions ==
                std::vector<std::size_t>{state.shape.rows, state.shape.columns} &&
                spec.elementCount == count,
            "phase output lost its complex [j,kl] shape");
    storage[index].assign(count, {713, -719});
    views[index] = {nullptr, count, storage[index].data()};
  }
  const auto retained = service->persistentBytes() + plan.persistentBytes();
  const double t0 = 17, t = 371.25;
  const std::array<std::uint8_t, 5> inactive{};
  require(bool(service->evaluate(plan, state.view(t, t0), views.data(), views.size(),
                                 inactive.data())), "inactive phase plan failed");
  for (const auto &output : storage)
    for (const auto value : output)
      require(value.real == 713 && value.imag == -719,
              "inactive phase output was modified");
  require(service->metrics().diagnosticWorkspaceHighWaterBytes == 0 &&
              service->metrics().diagnosticIntermediateReuseCount == 0,
          "inactive phase outputs allocated a phase table");
  const std::array<std::uint8_t, 5> onlyA0{0, 0, 0, 0, 1};
  require(bool(service->evaluate(plan, state.view(t, t0), views.data(), views.size(),
                                 onlyA0.data())), "A0-only diagnostic selection failed");
  require(service->metrics().diagnosticWorkspaceHighWaterBytes == 0 &&
              service->metrics().diagnosticIntermediateReuseCount == 0,
          "A0-only output allocated or reused an unrequested phase table");

  const auto evaluateAndCompare = [&](double time, double referenceTime) {
    const auto previousReuse = service->metrics().diagnosticIntermediateReuseCount;
    require(bool(service->evaluate(plan, state.view(time, referenceTime),
                                   views.data(), views.size())), "phase evaluation failed");
    const double pi = std::acos(-1.0);
    const double f = 2 * config.rotationRate * std::sin(config.latitude * pi / 180);
    for (std::size_t horizontal = 0; horizontal < state.shape.columns; ++horizontal) {
      const double kh = descriptor.fourierModes()[horizontal].Kh;
      for (std::size_t vertical = 0; vertical < state.shape.rows; ++vertical) {
        const auto index = vertical + state.shape.rows * horizontal;
        const double m = static_cast<double>(vertical) * pi / config.Lz;
        // Independent constant-N dispersion relation; j=0 follows the
        // model's h=1 convention, including otherwise unused wave slots.
        const double gravityTerm = vertical == 0 ? config.g * kh * kh :
            hydrostatic ? config.N0 * config.N0 * kh * kh / (m * m) :
            (config.N0 * config.N0 - f * f) * kh * kh / (m * m + kh * kh);
        const double angle = std::sqrt(f * f + gravityTerm) * (time - referenceTime);
        const double real = std::cos(angle), imag = std::sin(angle);
        requireClose(storage[0][index].real, real, "phase real differs from analytic dispersion");
        requireClose(storage[0][index].imag, imag, "phase does not use t-t0");
        requireClose(storage[1][index].real, real, "conjPhase real differs");
        requireClose(storage[1][index].imag, -imag, "conjPhase is not the conjugate");
        requireClose(std::hypot(storage[0][index].real, storage[0][index].imag), 1,
                     "phase is not unit magnitude");
        const auto ap = state.Ap[index], am = state.Am[index];
        requireClose(storage[2][index].real, ap.real * real - ap.imag * imag, "Apt phase product real");
        requireClose(storage[2][index].imag, ap.real * imag + ap.imag * real, "Apt phase product imaginary");
        requireClose(storage[3][index].real, am.real * real + am.imag * imag, "Amt conjugate product real");
        requireClose(storage[3][index].imag, am.imag * real - am.real * imag, "Amt conjugate product imaginary");
        require(storage[4][index].real == state.A0[index].real &&
                    storage[4][index].imag == state.A0[index].imag, "phase output changed A0t");
      }
    }
    require(service->metrics().diagnosticIntermediateReuseCount == previousReuse + 3,
            "phase/conjPhase/Apt/Amt did not share their phase table");
    require(service->metrics().diagnosticWorkspaceHighWaterBytes <= count * sizeof(WVComplex64) &&
                service->metrics().diagnosticWorkspaceLiveBytes == 0 &&
                service->metrics().diagnosticPrimitiveOutputCount == 0 &&
                service->metrics().fftExecutionCount == 0,
            "phase diagnostics retained scratch or reconstructed physical fields");
    require(service->persistentBytes() + plan.persistentBytes() == retained,
            "phase evaluation changed retained plan/service storage");
  };
  evaluateAndCompare(t0, t0);
  evaluateAndCompare(t, t0);
  const auto phaseAtTime = storage[0];
  evaluateAndCompare(t + 1000, t0 + 1000);
  require(std::memcmp(phaseAtTime.data(), storage[0].data(), count * sizeof(WVComplex64)) == 0,
          "phase depends on absolute time instead of elapsed time");
  evaluateAndCompare(t + 31.5, t0);
  require(std::memcmp(phaseAtTime.data(), storage[0].data(), count * sizeof(WVComplex64)) != 0,
          "phase reused a stale event time");
  require(std::memcmp(original.Ap.data(), state.Ap.data(), count * sizeof(WVComplex64)) == 0 &&
              std::memcmp(original.Am.data(), state.Am.data(), count * sizeof(WVComplex64)) == 0 &&
              std::memcmp(original.A0.data(), state.A0.data(), count * sizeof(WVComplex64)) == 0,
          "phase diagnostics modified coefficients");
  const auto phaseBeforeZeroing = storage[0];
  std::fill(state.Ap.begin(), state.Ap.end(), WVComplex64{});
  std::fill(state.Am.begin(), state.Am.end(), WVComplex64{});
  std::fill(state.A0.begin(), state.A0.end(), WVComplex64{});
  evaluateAndCompare(t + 31.5, t0);
  require(std::memcmp(phaseBeforeZeroing.data(), storage[0].data(), count * sizeof(WVComplex64)) == 0,
          "phase depends on wave-vortex amplitudes");

  WVFieldSamplingRequest profiles;
  profiles.kind = WVFieldSamplingKind::fixedVerticalProfiles;
  profiles.xIndices = {1}; profiles.yIndices = {1};
  WVFieldSamplingRequest positions;
  positions.kind = WVFieldSamplingKind::positions;
  positions.x = {0}; positions.y = {0}; positions.z = {0};
  for (const auto *name : {"phase", "conjPhase"}) {
    WVFieldEvaluationPlan invalid;
    require(!service->createPlan({{"sample", name, profiles}}, invalid), "spectral phase accepted profile sampling");
    require(!service->createPlan({{"sample", name, positions}}, invalid), "spectral phase accepted position sampling");
  }
  for (const auto *name : {"rho_nm"}) {
    WVFieldEvaluationPlan invalid;
    require(bool(service->createPlan({full(name)}, invalid)), "qualified density full-grid plan rejected");
    require(!service->createPlan({{"sample",name,positions}}, invalid), "density accepted unsupported position sampling");
  }
}

void verifyDerivedMovingSampling() {
  const auto config = configuration(6, 5, true, false);
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(bool(status), "derived moving service creation failed");
  const auto state = stateFor(config);
  const double dx = config.Lx / config.Nx;
  const double dy = config.Ly / config.Ny;
  const double dz = config.Lz / (config.Nz - 1);
  WVFieldSamplingRequest sampling;
  sampling.kind = WVFieldSamplingKind::positions;
  sampling.interpolation = WVPositionInterpolation::spline;
  sampling.x = {0.0, 2.0 * dx, config.Lx + 1.25 * dx, -0.4 * dx};
  sampling.y = {0.0, 3.0 * dy, -0.75 * dy, config.Ly + 1.1 * dy};
  sampling.z = {-config.Lz, -config.Lz + 2.0 * dz,
                -config.Lz + 4.0 * dz, 0.0};
  const std::array<const char *, 10> names{
      "p", "pi", "psi", "qgpv", "ssh", "ssu", "ssv",
      "zeta_x", "zeta_y", "zeta_z"};
  std::vector<WVFieldRequest> fixedRequests;
  std::vector<WVMovingFieldRequest> movingRequests;
  for (std::size_t index = 0; index < names.size(); ++index) {
    fixedRequests.push_back(
        {"fixed-" + std::to_string(index), names[index], sampling});
    movingRequests.push_back(
        {"moving-" + std::to_string(index), names[index], 0,
         sampling.x.size(), sampling.interpolation});
  }
  WVFieldEvaluationPlan fixed;
  status = service->createPlan(fixedRequests, fixed);
  require(bool(status), "derived fixed-position plan failed: " + status.message);
  std::vector<std::vector<double>> fixedStorage(names.size());
  std::vector<WVFieldOutputView> fixedViews;
  for (std::size_t index = 0; index < names.size(); ++index) {
    fixedStorage[index].resize(sampling.x.size());
    fixedViews.push_back({fixedStorage[index].data(), fixedStorage[index].size()});
  }
  status = service->evaluate(fixed, state.view(), fixedViews.data(), fixedViews.size());
  require(bool(status), "derived fixed-position evaluation failed");

  WVMovingFieldEvaluationPlan moving;
  status = service->createMovingPlan(movingRequests, moving);
  require(bool(status), "derived moving plan failed: " + status.message);
  const auto retained = moving.persistentBytes();
  const auto serviceRetained=service->persistentBytes();
  std::vector<std::vector<double>> movingStorage(
      names.size(), std::vector<double>(sampling.x.size(), -919.0));
  std::vector<WVFieldOutputView> movingViews;
  for (auto &output : movingStorage)
    movingViews.push_back({output.data(), output.size()});
  std::vector<std::uint8_t> active(names.size());
  for (std::size_t index = 0; index < active.size(); ++index)
    active[index] = index % 2;
  status = service->evaluateMoving(
      moving, state.view(),
      {sampling.x.data(), sampling.y.data(), sampling.z.data(),
       sampling.x.size()},
      movingViews.data(), movingViews.size(), active.data());
  require(bool(status), "selected derived moving evaluation failed");
  for (std::size_t index = 0; index < names.size(); ++index)
    require(active[index] ? movingStorage[index] == fixedStorage[index]
                          : std::all_of(movingStorage[index].begin(),
                                        movingStorage[index].end(),
                                        [](double value) { return value == -919.0; }),
            "derived moving active-output selection changed values");
  for (auto &output : movingStorage)
    std::fill(output.begin(), output.end(), -919.0);
  status = service->evaluateMoving(
      moving, state.view(),
      {sampling.x.data(), sampling.y.data(), sampling.z.data(),
       sampling.x.size()},
      movingViews.data(), movingViews.size());
  require(bool(status), "derived moving evaluation failed");
  require(movingStorage == fixedStorage,
          "derived moving fields differ from fixed-position interpolation");
  require(moving.persistentBytes() == retained &&
              service->persistentBytes()==serviceRetained &&
              service->metrics().diagnosticWorkspaceLiveBytes == 0,
          "derived moving evaluation retained temporary workspace");

  auto invalidZ = sampling.z;
  invalidZ[0] = std::numeric_limits<double>::quiet_NaN();
  for (auto &output : movingStorage)
    std::fill(output.begin(), output.end(), -919.0);
  status = service->evaluateMoving(
      moving, state.view(),
      {sampling.x.data(), sampling.y.data(), invalidZ.data(),
       sampling.x.size()},
      movingViews.data(), movingViews.size());
  require(!status && std::all_of(movingStorage.begin(), movingStorage.end(),
      [](const auto &output) { return std::all_of(output.begin(), output.end(),
          [](double value) { return value == -919.0; }); }),
      "invalid derived moving coordinates published partial output");
  require(service->metrics().diagnosticWorkspaceLiveBytes == 0,
          "failed derived moving evaluation retained workspace");

  std::unique_ptr<WVFieldEvaluationService> other;
  require(bool(WVFieldEvaluationService::create(
              config, std::make_unique<WVReferenceFFTEngine>(), other)),
          "foreign derived moving service creation failed");
  require(!other->evaluateMoving(
              moving, state.view(),
              {sampling.x.data(), sampling.y.data(), sampling.z.data(),
               sampling.x.size()},
              movingViews.data(), movingViews.size()),
          "derived moving plan crossed service ownership");
}

void verifyDiagnosticSamplingRoutes() {
  const auto config = configuration(6, 5, false, true);
  std::unique_ptr<WVFieldEvaluationService> service;
  require(bool(WVFieldEvaluationService::create(
              config, std::make_unique<WVReferenceFFTEngine>(), service)),
          "diagnostic sampling service creation failed");
  const auto state = stateFor(config);
  WVFieldSamplingRequest sampling;
  sampling.kind = WVFieldSamplingKind::positions;
  sampling.interpolation = WVPositionInterpolation::spline;
  sampling.x = {0.4 * config.Lx / config.Nx,
                config.Lx + 2.2 * config.Lx / config.Nx};
  sampling.y = {1.3 * config.Ly / config.Ny,
                -0.7 * config.Ly / config.Ny};
  sampling.z = {-0.65 * config.Lz, -0.2 * config.Lz};
  const std::array<const char *, 3> names{"u_g", "p_w", "ssh_io"};
  std::vector<WVFieldRequest> fixedRequests;
  std::vector<WVMovingFieldRequest> movingRequests;
  std::vector<WVEventFieldRequest> eventRequests;
  for (std::size_t index = 0; index < names.size(); ++index) {
    fixedRequests.push_back(
        {"fixed-diagnostic-" + std::to_string(index), names[index], sampling});
    movingRequests.push_back(
        {"moving-diagnostic-" + std::to_string(index), names[index], 0,
         sampling.x.size(), sampling.interpolation});
    eventRequests.push_back(
        {"event-diagnostic-" + std::to_string(index), names[index], 0,
         sampling.interpolation});
  }
  WVFieldEvaluationPlan fixed;
  require(bool(service->createPlan(fixedRequests, fixed)),
          "diagnostic fixed-position plan failed");
  std::vector<std::vector<double>> fixedStorage(
      names.size(), std::vector<double>(sampling.x.size()));
  std::vector<WVFieldOutputView> fixedViews;
  for (auto &field : fixedStorage)
    fixedViews.push_back({field.data(), field.size()});
  require(bool(service->evaluate(fixed, state.view(), fixedViews.data(),
                                 fixedViews.size())),
          "diagnostic fixed-position evaluation failed");

  WVFieldSamplingRequest profiles;
  profiles.kind = WVFieldSamplingKind::fixedVerticalProfiles;
  profiles.xIndices = {1, config.Nx};
  profiles.yIndices = {1, config.Ny};
  WVFieldEvaluationPlan fullAndProfiles;
  require(bool(service->createPlan(
              {{"full-u-g", "u_g", {}}, {"profile-u-g", "u_g", profiles},
               {"full-p-w", "p_w", {}}, {"profile-p-w", "p_w", profiles}},
              fullAndProfiles)),
          "diagnostic fixed-profile plan failed");
  std::vector<std::vector<double>> profileStorage;
  std::vector<WVFieldOutputView> profileViews;
  profileStorage.reserve(fullAndProfiles.outputCount());
  profileViews.reserve(fullAndProfiles.outputCount());
  for (const auto &output : fullAndProfiles.outputs()) {
    profileStorage.emplace_back(output.elementCount);
    profileViews.push_back(
        {profileStorage.back().data(), profileStorage.back().size()});
  }
  require(bool(service->evaluate(fullAndProfiles, state.view(),
                                 profileViews.data(), profileViews.size())),
          "diagnostic fixed-profile evaluation failed");
  const auto plane = config.Nx * config.Ny;
  for (const auto &pair : {std::pair<std::size_t, std::size_t>{0, 1},
                           std::pair<std::size_t, std::size_t>{2, 3}})
    for (std::size_t z = 0; z < config.Nz; ++z) {
      requireClose(profileStorage[pair.second][z],
                   profileStorage[pair.first][plane * z],
                   "diagnostic first profile differs from full field");
      requireClose(profileStorage[pair.second][z + config.Nz],
                   profileStorage[pair.first][config.Nx - 1 +
                       config.Nx * (config.Ny - 1) + plane * z],
                   "diagnostic last profile differs from full field");
    }

  WVMovingFieldEvaluationPlan moving;
  require(bool(service->createMovingPlan(movingRequests, moving)),
          "diagnostic moving plan failed");
  std::vector<std::vector<double>> movingStorage(
      names.size(), std::vector<double>(sampling.x.size()));
  std::vector<WVFieldOutputView> movingViews;
  for (auto &field : movingStorage)
    movingViews.push_back({field.data(), field.size()});
  require(bool(service->evaluateMoving(
              moving, state.view(),
              {sampling.x.data(), sampling.y.data(), sampling.z.data(),
               sampling.x.size()},
              movingViews.data(), movingViews.size())),
          "diagnostic moving evaluation failed");
  require(movingStorage == fixedStorage,
          "diagnostic moving fields differ from fixed positions");

  WVEventFieldEvaluationPlan event;
  require(bool(service->createEventPlan(eventRequests, event)),
          "diagnostic event plan failed");
  const std::size_t extents[]{2};
  WVEventPositionSetView positionSet{
      sampling.x.data(), sampling.y.data(), sampling.z.data(),
      sampling.x.size(), extents, 1};
  WVPreparedFieldGeometry geometry;
  require(bool(service->prepareEventGeometry(event, &positionSet, 1, geometry)),
          "diagnostic event geometry failed");
  std::vector<std::vector<double>> eventStorage(
      names.size(), std::vector<double>(sampling.x.size()));
  std::vector<WVFieldOutputView> eventViews;
  for (auto &field : eventStorage)
    eventViews.push_back({field.data(), field.size()});
  require(bool(service->evaluateEvent(event, geometry, state.view(),
                                      eventViews.data(), eventViews.size())),
          "diagnostic event evaluation failed");
  require(eventStorage == fixedStorage,
          "diagnostic event fields differ from fixed positions");

  auto movedEvent = std::move(event);
  for (auto &field : eventStorage)
    std::fill(field.begin(), field.end(), -919.0);
  require(bool(service->evaluateEvent(movedEvent, geometry, state.view(),
                                      eventViews.data(), eventViews.size())) &&
              eventStorage == fixedStorage,
          "prepared diagnostic geometry did not survive a plan move");
  WVEventFieldEvaluationPlan identicalEvent;
  require(bool(service->createEventPlan(eventRequests, identicalEvent)),
          "identical diagnostic event plan failed");
  for (auto &field : eventStorage)
    std::fill(field.begin(), field.end(), -919.0);
  require(!service->evaluateEvent(identicalEvent, geometry, state.view(),
                                  eventViews.data(), eventViews.size()) &&
              std::all_of(eventStorage.begin(), eventStorage.end(),
                          [](const auto &field) {
                            return std::all_of(
                                field.begin(), field.end(),
                                [](double value) { return value == -919.0; });
                          }),
          "prepared diagnostic geometry accepted an identical replacement plan");
  const auto secondEventView = eventViews[1];
  eventViews[1] = eventViews[0];
  require(service->evaluateEvent(movedEvent, geometry, state.view(),
                                 eventViews.data(), eventViews.size())
                  .code == WVKernelStatusCode::overlappingArrays &&
              std::all_of(eventStorage.front().begin(),
                          eventStorage.front().end(),
                          [](double value) { return value == -919.0; }),
          "sampled event output overlap was accepted or published output");
  eventViews[1] = secondEventView;
  event = std::move(movedEvent);
  for (std::size_t index = 0; index < eventStorage.size(); ++index)
    std::copy(fixedStorage[index].begin(), fixedStorage[index].end(),
              eventStorage[index].begin());

  WVEventFieldEvaluationPlan primitiveEvent;
  require(bool(service->createEventPlan(
              {{"primitive-u", "u", 0, WVPositionInterpolation::spline}},
              primitiveEvent)),
          "primitive sibling event plan failed");
  WVPreparedFieldGeometry primitiveGeometry;
  require(bool(service->prepareEventGeometry(
              primitiveEvent, &positionSet, 1, primitiveGeometry)),
          "primitive sibling event geometry failed");
  std::vector<double> primitive(sampling.x.size());
  WVFieldOutputView primitiveView{primitive.data(), primitive.size()};
  WVEventFieldEvaluationBatchEntry mixedEntries[]{
      {&primitiveEvent, &primitiveGeometry, &primitiveView, 1},
      {&event, &geometry, eventViews.data(), eventViews.size()}};
  require(bool(service->evaluateEventBatch(state.view(), mixedEntries, 2)),
          "mixed primitive/diagnostic event batch failed");
  require(eventStorage == fixedStorage,
          "mixed event batch changed diagnostic values");

  auto otherConfig = config;
  otherConfig.Nx = 7;
  std::unique_ptr<WVFieldEvaluationService> other;
  require(bool(WVFieldEvaluationService::create(
              otherConfig, std::make_unique<WVReferenceFFTEngine>(), other)),
          "foreign diagnostic event service creation failed");
  WVPreparedFieldGeometry foreignGeometry;
  require(!other->prepareEventGeometry(event, &positionSet, 1, foreignGeometry),
          "diagnostic event plan crossed same-family grid ownership");
}

void verifyFailureAndLifecycleContracts() {
  const auto config = configuration(6, 5, false, false);
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config,
      std::make_unique<FailureEngine>(WVKernelStatusCode::allocationFailure),
      service);
  require(status.code == WVKernelStatusCode::allocationFailure,
          "allocation failure was not propagated");
  status = WVFieldEvaluationService::create(
      config, std::make_unique<FailureEngine>(WVKernelStatusCode::success),
      service);
  require(static_cast<bool>(status), "execution-failure service creation failed");
  WVFieldEvaluationPlan plan;
  status = service->createPlan({full("u")}, plan);
  require(static_cast<bool>(status), "execution-failure plan creation failed");
  const auto state = stateFor(config);
  std::vector<double> output(plan.outputs()[0].elementCount, 19.0);
  WVFieldOutputView view{output.data(), output.size()};
  status = service->evaluate(plan, state.view(), &view, 1);
  require(status.code == WVKernelStatusCode::fftExecutionFailure,
          "FFT execution failure was not propagated");

  const auto activePlans = std::make_shared<std::size_t>(0);
  status = WVFieldEvaluationService::create(
      config, std::make_unique<CountingEngine>(activePlans), service);
  require(static_cast<bool>(status), "counting service creation failed");
  require(*activePlans > 0, "service did not retain its private FFT plans");
  service.reset();
  require(*activePlans == 0, "service destruction leaked FFT plans");
}

void verifyLowMemoryComponentFailureCleanup() {
  const auto config = configuration(6, 5, true, true);
  const auto control = std::make_shared<FailOnceControl>();
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<FailOnceEngine>(control), service);
  require(bool(status), "low-memory component failure service creation failed");
  require(bool(service->setVariableEvaluationPolicy(
              WVVariableEvaluationPolicy::lowMemory)),
          "low-memory component failure policy was rejected");
  WVFieldEvaluationPlan plan;
  require(bool(service->createPlan({full("u_g"), full("v_g")}, plan)),
          "low-memory component failure plan creation failed");
  std::array<std::vector<double>, 2> values;
  std::array<WVFieldOutputView, 2> views;
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index].resize(plan.outputs()[index].elementCount);
    views[index] = {values[index].data(), values[index].size()};
  }
  const auto stateStorage = stateFor(config);
  const WVIntegrationState state{stateStorage.view()};
  const auto retained = service->persistentBytes();
  const auto ledgerBefore = service->metrics().variableEvaluation;
  const auto producersBefore = service->producerMetrics();
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state, session)),
            "low-memory component failure session did not start");
    control->armed = true;
    status = service->evaluate(plan, state, views.data(), views.size());
    require(status.code == WVKernelStatusCode::fftExecutionFailure,
            "low-memory component query did not expose the injected failure");
    require(bool(service->evaluate(plan, state, views.data(), views.size())),
            "low-memory component query did not recover in the same session");
  }
  const auto ledgerAfter = service->metrics().variableEvaluation;
  const auto producersAfter = service->producerMetrics();
  require(ledgerAfter.producerExecutions == ledgerBefore.producerExecutions + 7 &&
              ledgerAfter.evictions == ledgerBefore.evictions + 7 &&
              ledgerAfter.recomputations == ledgerBefore.recomputations + 3 &&
              ledgerAfter.liveBytes == 0,
          "low-memory component failure leaked pins or skipped explicit recomputation");
  require(producersAfter.stateValidations == producersBefore.stateValidations + 3,
          "low-memory component retry reused a stale registered coefficient view");
  require(service->persistentBytes() == retained,
          "low-memory component failure discarded prepared complex capacity");
}

void verifyEvaluation(std::size_t nx, std::size_t ny, bool hydrostatic,
                      bool antialias) {
  const auto config = configuration(nx, ny, hydrostatic, antialias);
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status), "field service creation failed");

  std::vector<WVFieldRequest> requests;
  for (const auto &metadata : WVPortableVariableCatalog)
    if(metadata.ordinal>=3 && metadata.ordinal<23) requests.push_back(full(metadata.name));
  WVFieldSamplingRequest profiles;
  profiles.kind = WVFieldSamplingKind::fixedVerticalProfiles;
  profiles.xIndices = {1, nx};
  profiles.yIndices = {1, ny};
  requests.push_back({"u_profiles", "u", profiles});

  const double dx = config.Lx / static_cast<double>(nx);
  const double dy = config.Ly / static_cast<double>(ny);
  const double dz = config.Lz / static_cast<double>(config.Nz - 1);
  WVFieldSamplingRequest linear;
  linear.kind = WVFieldSamplingKind::positions;
  linear.interpolation = WVPositionInterpolation::linear;
  linear.x = {0.35 * dx, config.Lx + 0.35 * dx, 2.0 * dx};
  linear.y = {0.6 * dy, -config.Ly + 0.6 * dy, 3.0 * dy};
  linear.z = {-config.Lz + 1.4 * dz, -config.Lz + 1.4 * dz,
              -config.Lz + 2.0 * dz};
  requests.push_back({"u_linear", "u", linear});
  WVFieldSamplingRequest spline = linear;
  spline.interpolation = WVPositionInterpolation::spline;
  requests.push_back({"u_spline", "u", spline});

  WVFieldEvaluationPlan plan;
  status = service->createPlan(requests, plan);
  require(static_cast<bool>(status), "field plan creation failed: " + status.message);
  require(plan.outputCount() == requests.size(), "wrong plan output count");
  const auto planBytes = plan.persistentBytes();
  require(planBytes > sizeof(plan), "plan storage metric omitted payloads");

  std::vector<std::vector<double>> storage;
  std::vector<WVFieldOutputView> views;
  storage.reserve(plan.outputCount());
  views.reserve(plan.outputCount());
  std::size_t expectedWrites = 0;
  for (const auto &output : plan.outputs()) {
    storage.emplace_back(output.elementCount, -9.87654321e250);
    expectedWrites += output.elementCount;
  }
  for (auto &output : storage)
    views.push_back({output.data(), output.size()});
  const auto state = stateFor(config);
  status = service->evaluate(plan, state.view(), views.data(), views.size());
  require(static_cast<bool>(status), "field evaluation failed: " + status.message);
  for (const auto &output : storage)
    require(std::none_of(output.begin(), output.end(), [](double value) {
              return value == -9.87654321e250;
            }),
            "field evaluation left caller-owned elements unwritten");

  const auto &u = storage[outputIndex(plan, "full_u")];
  const auto &v = storage[outputIndex(plan, "full_v")];
  const auto &w = storage[outputIndex(plan, "full_w")];
  const auto &eta = storage[outputIndex(plan, "full_eta")];
  const auto &piField = storage[outputIndex(plan, "full_pi")];
  const auto &pressure = storage[outputIndex(plan, "full_p")];
  const auto &rhoE = storage[outputIndex(plan, "full_rho_e")];
  const auto &rhoTotal = storage[outputIndex(plan, "full_rho_total")];
  const auto &rhoBar = storage[outputIndex(plan, "full_rho_bar")];
  const auto &ssu = storage[outputIndex(plan, "full_ssu")];
  const auto &ssv = storage[outputIndex(plan, "full_ssv")];
  const auto &ssh = storage[outputIndex(plan, "full_ssh")];
  const auto horizontalCount = nx * ny;
  const auto fieldCount = horizontalCount * config.Nz;
  const double densityScale = config.rho0 * config.N0 * config.N0 / config.g;
  for (std::size_t index = 0; index < fieldCount; ++index) {
    requireClose(pressure[index], config.rho0 * config.g * piField[index],
                 "pressure scaling mismatch");
    requireClose(rhoE[index], densityScale * eta[index],
                 "density anomaly scaling mismatch");
  }
  for (std::size_t z = 0; z < config.Nz; ++z) {
    const double zCoordinate = -config.Lz + dz * static_cast<double>(z);
    double meanEta = 0.0;
    for (std::size_t horizontal = 0; horizontal < horizontalCount;
         ++horizontal) {
      const auto index = horizontal + horizontalCount * z;
      meanEta += eta[index];
      requireClose(rhoTotal[index],
                   config.rho0 - densityScale * zCoordinate +
                       densityScale * eta[index],
                   "total density mismatch");
    }
    meanEta /= static_cast<double>(horizontalCount);
    requireClose(rhoBar[z],
                 config.rho0 - densityScale * zCoordinate +
                     densityScale * meanEta,
                 "mean density mismatch");
  }
  const auto surfaceOffset = fieldCount - horizontalCount;
  for (std::size_t index = 0; index < horizontalCount; ++index) {
    requireClose(ssu[index], u[surfaceOffset + index],
                 "surface u mismatch");
    requireClose(ssv[index], v[surfaceOffset + index],
                 "surface v mismatch");
    requireClose(ssh[index], piField[surfaceOffset + index],
                 "surface height mismatch");
  }
  double expectedUVMax = 0.0;
  double expectedWMax = 0.0;
  for (std::size_t index = 0; index < fieldCount; ++index) {
    expectedUVMax = std::max(
        expectedUVMax, std::sqrt(u[index] * u[index] + v[index] * v[index]));
    expectedWMax = std::max(expectedWMax, std::abs(w[index]));
  }
  requireClose(storage[outputIndex(plan, "full_uvMax")][0], expectedUVMax,
               "uvMax mismatch");
  requireClose(storage[outputIndex(plan, "full_wMax")][0], expectedWMax,
               "wMax mismatch");

  const auto &profile = storage[outputIndex(plan, "u_profiles")];
  for (std::size_t z = 0; z < config.Nz; ++z) {
    requireClose(profile[z], u[horizontalCount * z],
                 "first one-based profile mismatch");
    requireClose(profile[z + config.Nz],
                 u[(nx - 1) + nx * (ny - 1) + horizontalCount * z],
                 "last one-based profile mismatch");
  }
  const auto &linearOutput = storage[outputIndex(plan, "u_linear")];
  const auto &splineOutput = storage[outputIndex(plan, "u_spline")];
  requireClose(linearOutput[0], linearOutput[1],
               "linear horizontal wrapping mismatch");
  requireClose(splineOutput[0], splineOutput[1],
               "spline horizontal wrapping mismatch");
  requireClose(linearOutput[2],
               u[2 + nx * 3 + horizontalCount * 2],
               "linear interpolation did not reproduce a grid knot");
  requireClose(splineOutput[2],
               u[2 + nx * 3 + horizontalCount * 2],
               "spline interpolation did not reproduce a grid knot");

  std::unique_ptr<WVTransformConstantStratificationKernel> derivativeKernel;
  status = WVTransformConstantStratificationKernel::create(
      config, std::make_unique<WVReferenceFFTEngine>(), derivativeKernel);
  require(static_cast<bool>(status), "derivative reference creation failed");
  std::vector<double> derivatives(3 * fieldCount);
  std::vector<double> expectedZetaX(fieldCount, 0.0);
  std::vector<double> expectedZetaY(fieldCount, 0.0);
  std::vector<double> expectedZetaZ(fieldCount, 0.0);
  WVRealFieldBundleView derivativeView{
      derivatives.data(), {nx, ny, config.Nz, 3}};
  status = derivativeKernel->transformStateFieldDerivatives(
      state.view(), WVDynamicalField::u, derivativeView);
  require(static_cast<bool>(status), "u derivative reference failed");
  for (std::size_t index = 0; index < fieldCount; ++index) {
    expectedZetaY[index] += derivatives[2 * fieldCount + index];
    expectedZetaZ[index] -= derivatives[fieldCount + index];
  }
  status = derivativeKernel->transformStateFieldDerivatives(
      state.view(), WVDynamicalField::v, derivativeView);
  require(static_cast<bool>(status), "v derivative reference failed");
  for (std::size_t index = 0; index < fieldCount; ++index) {
    expectedZetaX[index] -= derivatives[2 * fieldCount + index];
    expectedZetaZ[index] += derivatives[index];
  }
  status = derivativeKernel->transformStateFieldDerivatives(
      state.view(), WVDynamicalField::w, derivativeView);
  require(static_cast<bool>(status), "w derivative reference failed");
  for (std::size_t index = 0; index < fieldCount; ++index) {
    expectedZetaX[index] += derivatives[fieldCount + index];
    expectedZetaY[index] -= derivatives[index];
  }
  const auto &zetaX = storage[outputIndex(plan, "full_zeta_x")];
  const auto &zetaY = storage[outputIndex(plan, "full_zeta_y")];
  const auto &zetaZ = storage[outputIndex(plan, "full_zeta_z")];
  for (std::size_t index = 0; index < fieldCount; ++index) {
    requireClose(zetaX[index], expectedZetaX[index], "zeta_x mismatch");
    requireClose(zetaY[index], expectedZetaY[index], "zeta_y mismatch");
    requireClose(zetaZ[index], expectedZetaZ[index], "zeta_z mismatch");
  }

  const auto firstOutputs = storage;
  status = service->evaluate(plan, state.view(), views.data(), views.size());
  require(static_cast<bool>(status), "repeated evaluation failed");
  require(storage == firstOutputs, "repeated evaluation changed exact outputs");
  require(plan.persistentBytes() == planBytes,
          "evaluation mutated immutable plan storage");
  const auto &metrics = service->metrics();
  const auto coefficientCount = state.shape.elementCount();
  const auto expectedCapacity =
      6 * fieldCount * sizeof(double) +
      2 * coefficientCount * sizeof(WVComplex64);
  const auto expectedHighWater = std::max(
      6 * fieldCount * sizeof(double),
      5 * fieldCount * sizeof(double) +
          2 * coefficientCount * sizeof(WVComplex64));
  require(metrics.evaluationCount == 2, "wrong repeated evaluation count");
  require(metrics.coincidentBatchCount == 2, "wrong coincident batch count");
  require(metrics.transformCount == 14, "coincident transforms were not reused");
  require(metrics.primitiveFieldEvaluationCount == 20,
          "wrong primitive evaluation count");
  require(metrics.primitiveFieldReuseCount >= 6,
          "repeated field consumers did not register reuse");
  require(metrics.outputElementWriteCount == 2 * expectedWrites,
          "output write metric is not exact");
  require(metrics.scratchCapacityBytes == expectedCapacity,
          "scratch capacity metric is not exact");
  require(metrics.scratchHighWaterBytes == expectedHighWater,
          "scratch high-water metric is not exact");
  require(metrics.lastPlanBytes == planBytes &&
              metrics.maximumPlanBytes == planBytes,
          "plan storage metrics are not exact");
  require(metrics.servicePersistentBytes == service->persistentBytes(),
          "service persistent storage metric is not exact");
  require(metrics.catalogBytes == portableVariableCatalogBytes(),
          "static catalog storage metric is not exact");
  require(metrics.scratchHighWaterBytes <= metrics.scratchCapacityBytes,
          "scratch high-water exceeds bounded capacity");

  WVMovingFieldEvaluationPlan movingPlan;
  status = service->createMovingPlan(
      {{"moving-linear", "u", 0, linear.x.size(),
        WVPositionInterpolation::linear},
       {"moving-spline", "u", 0, spline.x.size(),
        WVPositionInterpolation::spline}},
      movingPlan);
  require(static_cast<bool>(status), "moving plan creation failed");
  std::vector<double> movingLinear(linear.x.size());
  std::vector<double> movingSpline(spline.x.size());
  std::array<WVFieldOutputView, 2> movingViews{{
      {movingLinear.data(), movingLinear.size()},
      {movingSpline.data(), movingSpline.size()}}};
  status = service->evaluateMoving(
      movingPlan, state.view(),
      {linear.x.data(), linear.y.data(), linear.z.data(), linear.x.size()},
      movingViews.data(), movingViews.size());
  require(static_cast<bool>(status), "moving evaluation failed");
  for (std::size_t index = 0; index < linear.x.size(); ++index) {
    requireClose(movingLinear[index], linearOutput[index],
                 "moving linear interpolation differs from fixed positions");
    requireClose(movingSpline[index], splineOutput[index],
                 "moving spline interpolation differs from fixed positions");
  }
  const auto movingWorkspaceBytes =
      service->metrics().movingInterpolationWorkspaceBytes;
  const auto splineFactorBytes =
      (config.Nx * config.Nx + config.Ny * config.Ny +
       config.Nz * config.Nz) *
          sizeof(double) +
      (config.Nx + config.Ny + config.Nz) * sizeof(std::size_t);
  require(movingWorkspaceBytes >= splineFactorBytes,
          "moving workspace ledger omitted spline factor storage");
  status = service->evaluateMoving(
      movingPlan, state.view(),
      {linear.x.data(), linear.y.data(), linear.z.data(), linear.x.size()},
      movingViews.data(), movingViews.size());
  require(static_cast<bool>(status) &&
              service->metrics().movingInterpolationWorkspaceBytes ==
                  movingWorkspaceBytes,
          "repeated moving evaluation changed bounded interpolation storage");

  WVFieldOutputView badShape{views[0].data, views[0].elementCount - 1};
  status = service->evaluate(plan, state.view(), &badShape, 1);
  require(status.code == WVKernelStatusCode::invalidShape,
          "wrong output-view count/shape was accepted");
  std::vector<WVFieldOutputView> overlapping = views;
  overlapping[1].data = overlapping[0].data;
  status = service->evaluate(plan, state.view(), overlapping.data(),
                             overlapping.size());
  require(status.code == WVKernelStatusCode::overlappingArrays,
          "overlapping outputs were accepted");
}

void verifyEventFieldEvaluation() {
  const auto config = configuration(6, 5, true, true);
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status), "event field service creation failed");
  const auto state = stateFor(config);
  const double dx = config.Lx / static_cast<double>(config.Nx);
  const double dy = config.Ly / static_cast<double>(config.Ny);
  const double dz = config.Lz / static_cast<double>(config.Nz - 1);
  const std::vector<double> x{0.35 * dx, config.Lx + 0.35 * dx,
                              2.0 * dx, 3.4 * dx};
  const std::vector<double> y{0.6 * dy, -config.Ly + 0.6 * dy,
                              3.0 * dy, 1.25 * dy};
  const std::vector<double> z{-config.Lz + 1.4 * dz,
                              -config.Lz + 1.4 * dz,
                              -config.Lz + 2.0 * dz, 25.0};
  const std::vector<std::size_t> extents{2, 2};

  std::vector<const WVPortableVariableMetadata *> positionFields;
  for (const auto &metadata : WVPortableVariableCatalog)
    if (metadata.kind == WVPortableVariableKind::field &&
        (metadata.samplingMask & portablePositionSampling) != 0)
      positionFields.push_back(&metadata);
  require(positionFields.size() == 16,
          "position-sampleable catalog coverage changed");

  auto compareWithFixed = [&](WVPositionInterpolation interpolation) {
    std::vector<WVEventFieldRequest> eventRequests;
    std::vector<WVFieldRequest> fixedRequests;
    eventRequests.reserve(positionFields.size());
    fixedRequests.reserve(positionFields.size());
    for (const auto *metadata : positionFields) {
      const std::string name = metadata->name;
      eventRequests.push_back(
          {"event-" + name, name, 0, interpolation});
      WVFieldSamplingRequest sampling;
      sampling.kind = WVFieldSamplingKind::positions;
      sampling.x = x;
      sampling.y = y;
      sampling.z = z;
      sampling.interpolation = interpolation;
      fixedRequests.push_back({"fixed-" + name, name, std::move(sampling)});
    }

    WVEventFieldEvaluationPlan eventPlan;
    const auto resolutionBefore =
        service->metrics().eventPlanFieldResolutionCount;
    status = service->createEventPlan(eventRequests, eventPlan);
    require(static_cast<bool>(status),
            "event field plan creation failed: " + status.message);
    require(eventPlan.outputCount() == positionFields.size() &&
                eventPlan.positionSetCount() == 1 &&
                eventPlan.requestedFieldMask() != 0 &&
                eventPlan.dependencyMask() != 0 &&
                eventPlan.fieldPlanFingerprint() != 0 &&
                eventPlan.persistentBytes() > sizeof(eventPlan),
            "event field plan omitted resolved construction metadata");
    require(service->metrics().eventPlanFieldResolutionCount ==
                resolutionBefore + positionFields.size(),
            "event fields were not resolved exactly once at construction");
    for (std::size_t index = 0; index < positionFields.size(); ++index) {
      const auto &output = eventPlan.outputs()[index];
      require(output.fieldIdentifier == positionFields[index]->identifier &&
                  output.naturalRank == positionFields[index]->naturalRank &&
                  output.dependencyMask ==
                      positionFields[index]->primitiveDependencyMask &&
                  output.positionSetSlot == 0 &&
                  output.interpolation == interpolation,
              "event plan changed a resolved field operation");
    }

    WVEventFieldEvaluationPlan equivalentPlan;
    status = service->createEventPlan(eventRequests, equivalentPlan);
    require(static_cast<bool>(status) &&
                equivalentPlan.fieldPlanFingerprint() ==
                    eventPlan.fieldPlanFingerprint(),
            "equivalent event field plans have different identities");
    auto changedRequests = eventRequests;
    changedRequests.front().interpolation =
        interpolation == WVPositionInterpolation::linear
            ? WVPositionInterpolation::spline
            : WVPositionInterpolation::linear;
    WVEventFieldEvaluationPlan changedPlan;
    status = service->createEventPlan(changedRequests, changedPlan);
    require(static_cast<bool>(status) &&
                changedPlan.fieldPlanFingerprint() !=
                    eventPlan.fieldPlanFingerprint(),
            "field-plan identity omitted its interpolation operation");

    WVEventPositionSetView set{x.data(), y.data(), z.data(), x.size(),
                               extents.data(), extents.size()};
    const auto xBefore = x;
    const auto yBefore = y;
    const auto zBefore = z;
    WVPreparedFieldGeometry prepared;
    status = service->prepareEventGeometry(eventPlan, &set, 1, prepared);
    require(static_cast<bool>(status),
            "event geometry preparation failed: " + status.message);
    require(prepared.fieldPlanFingerprint() ==
                    eventPlan.fieldPlanFingerprint() &&
                prepared.geometryFingerprint() != 0 &&
                prepared.positionSetCount() == 1 &&
                prepared.positionCount() == x.size() &&
                prepared.outputCount() == eventPlan.outputCount(),
            "prepared event geometry lost its resolved identity");
    require(prepared.positionSet(0).positionCount == x.size() &&
                prepared.positionSet(0).extentCount == extents.size(),
            "prepared geometry lost its event position set");
    for (const auto &output : prepared.outputs())
      require(output.dimensions == extents &&
                  output.elementCount == x.size(),
              "event output did not preserve its supplied logical extents");
    const auto preparedMetrics = prepared.metrics();
    require(preparedMetrics.positionSetCount == 1 &&
                preparedMetrics.positionCount == x.size() &&
                preparedMetrics.retainedBytes == prepared.retainedBytes() &&
                preparedMetrics.liveBytes == prepared.liveBytes() &&
                prepared.liveBytes() >= prepared.retainedBytes(),
            "prepared geometry storage metrics are incomplete");
    require(x == xBefore && y == yBefore && z == zBefore,
            "event geometry preparation mutated source coordinates");

    WVPreparedFieldGeometry samePrepared;
    status = service->prepareEventGeometry(eventPlan, &set, 1,
                                           samePrepared);
    require(static_cast<bool>(status) &&
                prepared.sameGeometry(samePrepared),
            "identical event geometry did not preserve cache identity");
    auto differentX = x;
    differentX[0] += 0.125 * dx;
    WVEventPositionSetView differentSet{
        differentX.data(), y.data(), z.data(), differentX.size(),
        extents.data(), extents.size()};
    WVPreparedFieldGeometry differentPrepared;
    status = service->prepareEventGeometry(eventPlan, &differentSet, 1,
                                           differentPrepared);
    require(static_cast<bool>(status) &&
                !prepared.sameGeometry(differentPrepared),
            "different event coordinates shared a geometry identity");

    WVFieldEvaluationPlan fixedPlan;
    status = service->createPlan(fixedRequests, fixedPlan);
    require(static_cast<bool>(status),
            "fixed comparison plan creation failed: " + status.message);
    std::vector<std::vector<double>> fixedStorage;
    std::vector<std::vector<double>> eventStorage;
    std::vector<WVFieldOutputView> fixedViews;
    std::vector<WVFieldOutputView> eventViews;
    fixedStorage.reserve(fixedPlan.outputCount());
    eventStorage.reserve(eventPlan.outputCount());
    fixedViews.reserve(fixedPlan.outputCount());
    eventViews.reserve(eventPlan.outputCount());
    for (std::size_t index = 0; index < fixedPlan.outputCount(); ++index) {
      fixedStorage.emplace_back(x.size(), -19.0);
      eventStorage.emplace_back(x.size(), -23.0);
    }
    for (auto &values : fixedStorage)
      fixedViews.push_back({values.data(), values.size()});
    for (auto &values : eventStorage)
      eventViews.push_back({values.data(), values.size()});
    status = service->evaluate(fixedPlan, state.view(), fixedViews.data(),
                               fixedViews.size());
    require(static_cast<bool>(status),
            "fixed comparison evaluation failed: " + status.message);
    const auto resolutionsAfterPreparation =
        service->metrics().eventPlanFieldResolutionCount;
    status = service->evaluateEvent(eventPlan, prepared, state.view(),
                                    eventViews.data(), eventViews.size());
    require(static_cast<bool>(status),
            "event field evaluation failed: " + status.message);
    require(service->metrics().eventPlanFieldResolutionCount ==
                resolutionsAfterPreparation,
            "event evaluation repeated construction-time field resolution");
    for (std::size_t field = 0; field < positionFields.size(); ++field)
      for (std::size_t position = 0; position < x.size(); ++position)
        requireWithinOneE12(
            eventStorage[field][position], fixedStorage[field][position],
            std::string("event/fixed mismatch for ") +
                positionFields[field]->name);

    std::vector<std::vector<double>> untouched(
        eventPlan.outputCount(), std::vector<double>(x.size(), 713.0));
    std::vector<WVFieldOutputView> badViews;
    badViews.reserve(untouched.size());
    for (auto &values : untouched)
      badViews.push_back({values.data(), values.size()});
    --badViews.back().elementCount;
    status = service->evaluateEvent(eventPlan, prepared, state.view(),
                                    badViews.data(), badViews.size());
    require(status.code == WVKernelStatusCode::invalidShape &&
                std::all_of(untouched.begin(), untouched.end(),
                            [](const auto &values) {
                              return std::all_of(
                                  values.begin(), values.end(),
                                  [](double value) { return value == 713.0; });
                            }),
            "invalid event outputs mutated caller storage");
  };

  compareWithFixed(WVPositionInterpolation::linear);
  compareWithFixed(WVPositionInterpolation::spline);

  WVEventFieldEvaluationPlan volumePlan;
  status = service->createEventPlan(
      {{"volume", "u", 0, WVPositionInterpolation::linear}}, volumePlan);
  require(static_cast<bool>(status), "volume event plan creation failed");
  WVEventPositionSetView missingZ{x.data(), y.data(), nullptr, x.size(),
                                  extents.data(), extents.size()};
  WVPreparedFieldGeometry validationGeometry;
  status = service->prepareEventGeometry(volumePlan, &missingZ, 1,
                                         validationGeometry);
  require(status.code == WVKernelStatusCode::invalidPointer,
          "volume event geometry accepted missing z coordinates");
  const std::vector<std::size_t> wrongExtents{3, 2};
  WVEventPositionSetView badExtent{x.data(), y.data(), z.data(), x.size(),
                                   wrongExtents.data(), wrongExtents.size()};
  status = service->prepareEventGeometry(volumePlan, &badExtent, 1,
                                         validationGeometry);
  require(status.code == WVKernelStatusCode::invalidShape,
          "event geometry accepted incompatible logical extents");
  auto nonfiniteX = x;
  nonfiniteX[1] = std::numeric_limits<double>::quiet_NaN();
  WVEventPositionSetView nonfiniteSet{
      nonfiniteX.data(), y.data(), z.data(), nonfiniteX.size(),
      extents.data(), extents.size()};
  status = service->prepareEventGeometry(volumePlan, &nonfiniteSet, 1,
                                         validationGeometry);
  require(status.code == WVKernelStatusCode::invalidConfiguration,
          "event geometry accepted a nonfinite coordinate");

  WVEventFieldEvaluationPlan horizontalPlan;
  status = service->createEventPlan(
      {{"surface-u", "ssu", 0, WVPositionInterpolation::linear},
       {"surface-v", "ssv", 0, WVPositionInterpolation::linear},
       {"surface-height", "ssh", 0, WVPositionInterpolation::linear}},
      horizontalPlan);
  require(static_cast<bool>(status),
          "derived horizontal event plan creation failed");
  WVPreparedFieldGeometry horizontalGeometry;
  status = service->prepareEventGeometry(horizontalPlan, &missingZ, 1,
                                         horizontalGeometry);
  require(static_cast<bool>(status),
          "horizontal fields required irrelevant z coordinates");
  std::vector<std::vector<double>> horizontalStorage(
      horizontalPlan.outputCount(), std::vector<double>(x.size()));
  std::vector<WVFieldOutputView> horizontalViews;
  for (auto &values : horizontalStorage)
    horizontalViews.push_back({values.data(), values.size()});
  status = service->evaluateEvent(horizontalPlan, horizontalGeometry,
                                  state.view(), horizontalViews.data(),
                                  horizontalViews.size());
  require(static_cast<bool>(status),
          "derived horizontal event evaluation failed");

  WVEventFieldEvaluationPlan unsupportedPlan;
  status = service->createEventPlan(
      {{"scalar", "energy", 0, WVPositionInterpolation::linear}},
      unsupportedPlan);
  require(status.code == WVKernelStatusCode::unsupportedOperation,
          "event plan accepted a non-position-sampleable field");

  const std::vector<std::size_t> zeroExtents{0};
  WVEventPositionSetView emptySet{nullptr, nullptr, nullptr, 0,
                                  zeroExtents.data(), zeroExtents.size()};
  WVPreparedFieldGeometry emptyGeometry;
  status = service->prepareEventGeometry(volumePlan, &emptySet, 1,
                                         emptyGeometry);
  require(static_cast<bool>(status) && emptyGeometry.positionCount() == 0 &&
              emptyGeometry.outputs()[0].dimensions == zeroExtents &&
              emptyGeometry.outputs()[0].elementCount == 0,
          "zero-length event geometry was not preserved");
  WVFieldOutputView emptyOutput{nullptr, 0};
  const auto transformsBeforeEmpty = service->metrics().transformCount;
  status = service->evaluateEvent(volumePlan, emptyGeometry, state.view(),
                                  &emptyOutput, 1);
  require(static_cast<bool>(status) &&
              service->metrics().transformCount == transformsBeforeEmpty,
          "zero-length event performed unnecessary field reconstruction");

  const std::vector<double> longerX{0.0, dx, 2.0 * dx, 3.0 * dx,
                                     4.0 * dx, 5.0 * dx};
  const std::vector<double> longerY{0.0,      dy,       2.0 * dy,
                                    3.0 * dy, 4.0 * dy, 0.5 * dy};
  const std::vector<double> longerZ(longerX.size(), -0.5 * config.Lz);
  const std::vector<std::size_t> longerExtents{3, 2};
  WVEventPositionSetView longerSet{longerX.data(),       longerY.data(),
                                   longerZ.data(),       longerX.size(),
                                   longerExtents.data(), longerExtents.size()};
  WVPreparedFieldGeometry longerGeometry;
  status =
      service->prepareEventGeometry(volumePlan, &longerSet, 1, longerGeometry);
  require(static_cast<bool>(status) && longerGeometry.positionCount() == 6 &&
              longerGeometry.outputs()[0].dimensions == longerExtents &&
              service->metrics().maximumPreparedGeometryRetainedBytes >=
                  longerGeometry.retainedBytes() &&
              service->metrics().maximumPreparedGeometryLiveBytes >=
                  longerGeometry.liveBytes(),
          "variable event geometry or high-water storage metrics are wrong");

  WVEventFieldEvaluationPlan batchPlan;
  status = service->createEventPlan(
      {{"batch-u", "u", 0, WVPositionInterpolation::spline},
       {"batch-eta", "eta", 0, WVPositionInterpolation::spline}},
      batchPlan);
  require(static_cast<bool>(status), "event batch plan creation failed");
  WVPreparedFieldGeometry firstBatchGeometry;
  WVPreparedFieldGeometry secondBatchGeometry;
  WVEventPositionSetView firstBatchSet{
      x.data(), y.data(), z.data(), x.size(), extents.data(), extents.size()};
  status = service->prepareEventGeometry(batchPlan, &firstBatchSet, 1,
                                         firstBatchGeometry);
  require(static_cast<bool>(status), "first event batch geometry failed");
  status = service->prepareEventGeometry(batchPlan, &longerSet, 1,
                                         secondBatchGeometry);
  require(static_cast<bool>(status), "second event batch geometry failed");
  std::array<std::vector<double>, 2> firstBatchStorage{
      std::vector<double>(x.size()), std::vector<double>(x.size())};
  std::array<std::vector<double>, 2> secondBatchStorage{
      std::vector<double>(longerX.size()), std::vector<double>(longerX.size())};
  std::array<WVFieldOutputView, 2> firstBatchViews{
      {{firstBatchStorage[0].data(), firstBatchStorage[0].size()},
       {firstBatchStorage[1].data(), firstBatchStorage[1].size()}}};
  std::array<WVFieldOutputView, 2> secondBatchViews{
      {{secondBatchStorage[0].data(), secondBatchStorage[0].size()},
       {secondBatchStorage[1].data(), secondBatchStorage[1].size()}}};
  const std::array<WVEventFieldEvaluationBatchEntry, 2> batchEntries{
      {{&batchPlan, &firstBatchGeometry, firstBatchViews.data(),
        firstBatchViews.size()},
       {&batchPlan, &secondBatchGeometry, secondBatchViews.data(),
        secondBatchViews.size()}}};
  const auto metricsBeforeBatch = service->metrics();
  status = service->evaluateEventBatch(state.view(), batchEntries.data(),
                                       batchEntries.size());
  require(static_cast<bool>(status),
          "coincident event field batch evaluation failed");
  const auto &metricsAfterBatch = service->metrics();
  require(metricsAfterBatch.transformCount ==
              metricsBeforeBatch.transformCount + 1,
          "event batch repeated primitive reconstruction");
  require(metricsAfterBatch.primitiveFieldEvaluationCount ==
              metricsBeforeBatch.primitiveFieldEvaluationCount + 2,
          "event batch reported the wrong primitive field count");
  require(metricsAfterBatch.splineInterpolationCount ==
              metricsBeforeBatch.splineInterpolationCount +
                  2 * (x.size() + longerX.size()),
          "event batch skipped an occurrence interpolation");
  require(metricsAfterBatch.eventEvaluationCount ==
                  metricsBeforeBatch.eventEvaluationCount + 2 &&
              metricsAfterBatch.eventBatchEvaluationCount ==
                  metricsBeforeBatch.eventBatchEvaluationCount + 1 &&
              metricsAfterBatch.eventBatchOccurrenceCount ==
                  metricsBeforeBatch.eventBatchOccurrenceCount + 2 &&
              metricsAfterBatch.eventBatchOutputCount ==
                  metricsBeforeBatch.eventBatchOutputCount + 4,
          "event batch evaluation counters are not occurrence-exact");
  require(metricsAfterBatch.eventBatchInvocationWorkspaceBytes >=
                  batchEntries.size() *
                      (2 * sizeof(void *) + sizeof(std::size_t)) &&
              metricsAfterBatch.servicePersistentBytes ==
                  service->persistentBytes(),
          "event batch omitted physical invocation workspace storage");
}

// MATLAB's three-dimensional transform explicitly supplies extrapval=0 to
// interpn after circular shifts. On a six-point y grid the four-cell shift
// can still leave a query beyond the final grid knot.
void verifyShiftedSplineZero(WVFieldEvaluationService &service,
                             const WVIntegrationState &state,
                             double lx, double ly, double depth) {
  WVFieldSamplingRequest sampling;
  sampling.kind = WVFieldSamplingKind::positions;
  sampling.interpolation = WVPositionInterpolation::spline;
  sampling.x = {0.0, 0.0, lx, 0.0};
  sampling.y = {0.0, 1.75 * ly / 6.0, 0.0, 1.75 * ly / 6.0 + ly};
  sampling.z.assign(4, -depth);
  sampling.z[2] = 1.0;
  WVFieldEvaluationPlan fixed;
  auto status = service.createPlan({{"density", "rho_total", sampling}}, fixed);
  require(bool(status), status.message);
  std::array<double, 4> values{};
  WVFieldOutputView output{values.data(), values.size()};
  auto verify = [&] {
    require(values[0] > 1000.0, "in-grid spline control lost background density");
    for (std::size_t i = 1; i < values.size(); ++i)
      require(values[i] == 0.0, "shifted spline extrapolation must be zero");
  };
  status = service.evaluate(fixed, state, &output, 1);
  require(bool(status), status.message);
  verify();
  WVMovingFieldEvaluationPlan moving;
  status = service.createMovingPlan({{"density", "rho_total", 0, 4,
                                      WVPositionInterpolation::spline}}, moving);
  require(bool(status), status.message);
  status = service.evaluateMoving(moving, state,
      {sampling.x.data(), sampling.y.data(), sampling.z.data(), 4}, &output, 1);
  require(bool(status), status.message);
  verify();
  WVEventFieldEvaluationPlan event;
  status = service.createEventPlan({{"density", "rho_total", 0,
                                     WVPositionInterpolation::spline}}, event);
  require(bool(status), status.message);
  const std::size_t extent = 4;
  WVEventPositionSetView positions{sampling.x.data(), sampling.y.data(),
                                  sampling.z.data(), 4, &extent, 1};
  WVPreparedFieldGeometry geometry;
  status = service.prepareEventGeometry(event, &positions, 1, geometry);
  require(bool(status), status.message);
  status = service.evaluateEvent(event, geometry, state, &output, 1);
  require(bool(status), status.message);
  verify();
}

void verifyBarotropicSplineExtrapolation() {
  WVTransformBarotropicQGConfiguration config;
  config.Nx = 8; config.Ny = 6; config.Lx = 8; config.Ly = 6;
  config.h = 0.8; config.g = 9.81; config.planetaryRadius = 6.371e6;
  config.rotationRate = 7.2921e-5; config.latitude = 33;
  WVTransformBarotropicQGDescriptor descriptor;
  auto status = WVTransformBarotropicQGDescriptor::create(config, descriptor);
  require(bool(status), status.message);
  std::unique_ptr<WVFieldEvaluationService> service;
  status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(bool(status), status.message);
  std::vector<WVComplex64> a0(descriptor.Nkl());
  WVCoefficientFamilyLayout family{"A0", {descriptor.Nkl()}};
  family.elementCount = descriptor.Nkl();
  WVCoefficientFamilyConstView coefficients{&family, a0.data()};
  WVIntegrationState state;
  state.coefficientFamilies = &coefficients;
  state.coefficientFamilyCount = 1;
  WVMovingFieldEvaluationPlan plan;
  status = service->createMovingPlan({{"u", "u", 0, 3,
                                      WVPositionInterpolation::spline}}, plan);
  require(bool(status), status.message);
  // Inverse circular shift of X^3+2Y^3, so the actual four-cell shifted
  // interpolation grid is an exact cubic with an independent analytic oracle.
  std::vector<double> velocity(8 * 6 * 2);
  for (std::size_t y = 0; y < 6; ++y)
    for (std::size_t x = 0; x < 8; ++x) {
      const double xx = (x + 4) % 8, yy = (y + 4) % 6;
      velocity[x + 8 * y] = xx * xx * xx + 2 * yy * yy * yy;
    }
  const std::array<double, 3> x{6.5, 6.5, 6.5}, y{1.75, 1.0, 7.75};
  std::array<double, 3> values{};
  WVFieldOutputView output{values.data(), values.size()};
  status = service->evaluateMovingFromAdvectionFields(plan, state,
      {velocity.data(), {8, 6, 1, 2}}, {x.data(), y.data(), nullptr, 3},
      &output, 1);
  require(bool(status), status.message);
  // MATLAB BQG omits extrapval: interpn(...,'spline') extends the end cubic.
  requireClose(values[0], 395.84375, "BQG spline end-piece extrapolation");
  requireClose(values[1], 265.625, "BQG spline final-knot control");
  requireClose(values[2], values[0], "BQG extrapolation periodic wrapping");
}

void verifySmallGridSplineBoundaries() {
  auto config = configuration(8, 6, true, true);
  auto constantState = stateFor(config);
  std::fill(constantState.Ap.begin(), constantState.Ap.end(), WVComplex64{});
  std::fill(constantState.Am.begin(), constantState.Am.end(), WVComplex64{});
  std::fill(constantState.A0.begin(), constantState.A0.end(), WVComplex64{});
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(bool(status), status.message);
  verifyShiftedSplineZero(*service, {constantState.view()}, config.Lx,
                          config.Ly, config.Lz);

  wavevortex::test_fixture::Temporary file;
  wavevortex::test_fixture::fixture(file.path);
  std::shared_ptr<const WVStratifiedModalRecord> source;
  const auto readStatus = WVStratifiedModalReader::read(file.path.string(), source);
  require(bool(readStatus), readStatus.message);
  status = WVFieldEvaluationService::create(
      source, std::make_unique<WVReferenceFFTEngine>(), service);
  require(bool(status), status.message);
  WVEventFieldEvaluationPlan unavailableDensityEvent;
  require(!service->createEventPlan(
              {{"unavailable-density", "eta_true", 0,
                WVPositionInterpolation::linear}},
              unavailableDensityEvent),
          "stratified-QG event plan accepted an unavailable density diagnostic");
  const auto &g = source->geometry();
  std::vector<WVComplex64> a0(g.Nj * g.Nkl);
  WVCoefficientFamilyLayout family{"A0", {g.Nj, g.Nkl}};
  WVCoefficientFamilyConstView coefficients{&family, a0.data()};
  WVIntegrationState state;
  state.coefficientFamilies = &coefficients;
  state.coefficientFamilyCount = 1;
  verifyShiftedSplineZero(*service, state, g.Lx, g.Ly, g.Lz);
}

void verifyVariableEvaluationSessions() {
  const auto config=configuration(6,5,true,true);
  const auto owned=stateFor(config);
  WVIntegrationState state{owned.view()};
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status=WVFieldEvaluationService::create(
      config,std::make_unique<WVReferenceFFTEngine>(),service);
  require(bool(status),status.message);
  WVFieldEvaluationPlan first,second,derivedFull,derivedPoints;
  WVMovingFieldEvaluationPlan moving;
  status=service->createPlan({full("u")},first);
  require(bool(status),status.message);
  WVFieldSamplingRequest points;
  points.kind=WVFieldSamplingKind::positions;
  points.x={0.25}; points.y={0.5}; points.z={-0.75};
  status=service->createPlan({{"point-u","u",points}},second);
  require(bool(status),status.message);
  status=service->createPlan({full("p"),full("rho_e")},derivedFull);
  require(bool(status),status.message);
  status=service->createPlan({{"point-p","p",points},
      {"point-rho","rho_e",points}},derivedPoints);
  require(bool(status),status.message);
  status=service->createMovingPlan({{"moving-u","u",0,1},
      {"moving-v","v",0,1}},moving);
  require(bool(status),status.message);
  std::vector<double> fullValues(first.outputs()[0].elementCount);
  std::vector<double> pointValues(second.outputs()[0].elementCount);
  WVFieldOutputView fullView{fullValues.data(),fullValues.size()};
  WVFieldOutputView pointView{pointValues.data(),pointValues.size()};
  std::array<std::vector<double>,2> derivedFullValues,
      derivedPointValues;
  std::array<WVFieldOutputView,2> derivedFullViews,derivedPointViews;
  for(std::size_t index=0;index<2;++index) {
    derivedFullValues[index].resize(
        derivedFull.outputs()[index].elementCount);
    derivedPointValues[index].resize(
        derivedPoints.outputs()[index].elementCount);
    derivedFullViews[index]={derivedFullValues[index].data(),
        derivedFullValues[index].size()};
    derivedPointViews[index]={derivedPointValues[index].data(),
        derivedPointValues[index].size()};
  }

  const auto standaloneBefore=service->metrics().variableEvaluation;
  const auto standaloneProducerBefore=service->producerMetrics();
  std::array<double,2> movingValues{};
  WVFieldOutputView movingViews[]={{&movingValues[0],1},{&movingValues[1],1}};
  double movingX=points.x[0],movingY=points.y[0],movingZ=points.z[0];
  const WVMovingPositionView movingPosition{&movingX,&movingY,&movingZ,1};
  require(bool(service->evaluate(first,state,&fullView,1)) &&
      !service->evaluationSessionActive() &&
      bool(service->evaluate(first,state,&fullView,1)) &&
      !service->evaluationSessionActive(),
      "standalone field evaluation did not close its event scope");
  require(bool(service->evaluateMoving(moving,state,movingPosition,
          movingViews,2)) && !service->evaluationSessionActive() &&
      bool(service->evaluateMoving(moving,state,movingPosition,
          movingViews,2)) && !service->evaluationSessionActive(),
      "standalone moving evaluation did not close its event scope");
  const auto standaloneAfter=service->metrics().variableEvaluation;
  const auto standaloneProducerAfter=service->producerMetrics();
  require(standaloneAfter.contexts==standaloneBefore.contexts+4 &&
      standaloneAfter.producerExecutions==
          standaloneBefore.producerExecutions+4 &&
      standaloneProducerAfter.reconstructions[0][0][0]==
          standaloneProducerBefore.reconstructions[0][0][0]+4,
      "standalone field calls did not use one fresh producer scope each");
  const std::uint8_t inactiveField=0;
  const std::array<std::uint8_t,2> inactiveMoving{};
  require(bool(service->evaluate(first,state,&fullView,1,&inactiveField)) &&
      bool(service->evaluateMoving(moving,state,movingPosition,movingViews,2,
          inactiveMoving.data())) &&
      service->metrics().variableEvaluation.contexts==standaloneAfter.contexts &&
      !service->evaluationSessionActive(),
      "inactive standalone calls opened an event scope");

  const auto preparedArenaBytes=service->metrics().eventFieldArenaPlannedBytes;
  const auto preparedPersistentBytes=service->persistentBytes();
  auto before=service->metrics().variableEvaluation;
  const auto producerBefore=service->producerMetrics();
  {
    WVFieldEvaluationSession session;
    status=service->beginEvaluationSession(state,session);
    require(bool(status) && session.active(),"reuse evaluation session did not start");
    require(!service->setVariableEvaluationPolicy(WVVariableEvaluationPolicy::lowMemory),
            "active session allowed its policy to change");
    require(bool(service->evaluate(first,state,&fullView,1)),"reuse session first query failed");
    require(bool(service->evaluate(second,state,&pointView,1)),"reuse session second query failed");
    require(bool(service->evaluate(derivedFull,state,derivedFullViews.data(),2)) &&
            bool(service->evaluate(derivedPoints,state,derivedPointViews.data(),2)),
        "reuse session derived queries failed");
    auto foreign=state;
    foreign.waveVortex.t+=1;
    require(!service->evaluate(first,foreign,&fullView,1),
            "evaluation session accepted a foreign immutable state");
  }
  auto after=service->metrics().variableEvaluation;
  const auto producerAfterReuse=service->producerMetrics();
  require(after.contexts==before.contexts+1 &&
              after.producerExecutions==before.producerExecutions+4 &&
              after.cacheHits>=before.cacheHits+6 &&
              after.duplicateExecutions==before.duplicateExecutions,
          "reuse session did not execute one shared natural-grid producer");
  require(producerAfterReuse.phasePreparations==producerBefore.phasePreparations+1 &&
              producerAfterReuse.reconstructions[0][0][0]==
                  producerBefore.reconstructions[0][0][0]+1 &&
              producerAfterReuse.reconstructions[5][0][0]==
                  producerBefore.reconstructions[5][0][0]+1 &&
              producerAfterReuse.reconstructions[8][0][0]==
                  producerBefore.reconstructions[8][0][0]+1,
          "reuse session did not suppress the duplicate kernel producer");
  require(service->persistentBytes()==preparedPersistentBytes &&
      service->metrics().eventFieldArenaPlannedBytes==preparedArenaBytes &&
      service->metrics().eventFieldArenaPeakBytes<=preparedArenaBytes,
      "prepared repeated field queries grew the output-event arena");

  require(bool(service->setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy::lowMemory)),"low-memory policy was rejected");
  before=after;
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
            "low-memory evaluation session did not start");
    require(bool(service->evaluate(first,state,&fullView,1)),
            "low-memory first query failed");
    require(bool(service->evaluate(second,state,&pointView,1)),
            "low-memory second query failed");
  }
  after=service->metrics().variableEvaluation;
  const auto producerAfterLowMemory=service->producerMetrics();
  require(after.contexts==before.contexts+1 &&
              after.producerExecutions==before.producerExecutions+2 &&
              after.recomputations==before.recomputations+1 &&
              after.evictions==before.evictions+2 &&
              after.duplicateExecutions==before.duplicateExecutions,
          "low-memory session did not explicitly evict and recompute");
  require(producerAfterLowMemory.reconstructions[0][0][0]==
              producerAfterReuse.reconstructions[0][0][0]+2,
          "low-memory session did not recompute the kernel producer");

  require(bool(service->setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy::reuse)),"reuse policy restore failed");
  WVFieldEvaluationPlan components;
  require(bool(service->createPlan({full("u_g"),full("v_g")},components)),
      "component session plan creation failed");
  std::vector<double> componentU(components.outputs()[0].elementCount);
  std::vector<double> componentV(components.outputs()[1].elementCount);
  WVFieldOutputView componentViews[]={{componentU.data(),componentU.size()},
      {componentV.data(),componentV.size()}};
  const auto componentBefore=service->producerMetrics();
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "component evaluation session did not start");
    require(bool(service->evaluate(components,state,componentViews,2)) &&
        bool(service->evaluate(components,state,componentViews,2)),
        "component evaluation session failed");
  }
  const auto componentAfter=service->producerMetrics();
  require(componentAfter.phasePreparations==componentBefore.phasePreparations+1 &&
      componentAfter.reconstructions[0][0][1]==
          componentBefore.reconstructions[0][0][1]+1 &&
      componentAfter.reconstructions[1][0][1]==
          componentBefore.reconstructions[1][0][1]+1,
      "component outputs repeated phase or reconstruction producers");

  const std::array<const char*,6> energyNames{
      "energy","geostrophicEnergy","energy_g","energy_w","energy_io",
      "energy_mda"};
  std::array<WVFieldEvaluationPlan,6> energySingles;
  std::vector<WVFieldRequest> energyRequests;
  for(std::size_t index=0;index<energyNames.size();++index) {
    require(bool(service->createPlan({full(energyNames[index])},
            energySingles[index])),
        "component energy single plan creation failed");
    energyRequests.push_back(full(energyNames[index]));
  }
  WVFieldEvaluationPlan energyBatch;
  require(bool(service->createPlan(energyRequests,energyBatch)),
      "component energy batch plan creation failed");
  std::reverse(energyRequests.begin(),energyRequests.end());
  WVFieldEvaluationPlan reverseEnergyBatch;
  require(bool(service->createPlan(energyRequests,reverseEnergyBatch)),
      "reverse component energy batch plan creation failed");
  std::array<double,6> energyReference{};
  for(std::size_t index=0;index<energyNames.size();++index) {
    WVFieldOutputView output{&energyReference[index],1};
    const auto energyStatus=service->evaluate(
        energySingles[index],state,&output,1);
    require(bool(energyStatus),std::string("component energy reference ")+
        energyNames[index]+" failed: "+energyStatus.message);
  }
  require(energyReference[1]==energyReference[2],
      "geostrophic energy aliases produced different references");
  for(const auto policy:{WVVariableEvaluationPolicy::reuse,
          WVVariableEvaluationPolicy::lowMemory}) {
    require(bool(service->setVariableEvaluationPolicy(policy)),
        "component energy evaluation policy change failed");
    std::array<double,6> batched{},reversed{};
    std::array<WVFieldOutputView,6> batchViews{},reverseViews{};
    for(std::size_t index=0;index<energyNames.size();++index) {
      batchViews[index]={&batched[index],1};
      reverseViews[index]={&reversed[index],1};
    }
    require(bool(service->evaluate(energyBatch,state,batchViews.data(),
                batchViews.size())) &&
            bool(service->evaluate(reverseEnergyBatch,state,reverseViews.data(),
                reverseViews.size())),
        "component energy batched evaluation failed");
    for(std::size_t index=0;index<energyNames.size();++index) {
      require(batched[index]==energyReference[index],
          "component energy batching changed a value");
      require(reversed[energyNames.size()-1-index]==energyReference[index],
          "reverse component energy batching changed a value");
    }
    for(const bool reverse:{false,true}) {
      WVFieldEvaluationSession session;
      require(bool(service->beginEvaluationSession(state,session)),
          "component energy ordered session did not start");
      for(std::size_t ordinal=0;ordinal<energyNames.size();++ordinal) {
        const auto index=reverse ? energyNames.size()-1-ordinal : ordinal;
        double value=0;
        WVFieldOutputView output{&value,1};
        require(bool(service->evaluate(energySingles[index],state,&output,1)) &&
                value==energyReference[index],
            "component energy session order changed a value");
      }
    }
  }
  require(bool(service->setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy::reuse)),
      "reuse policy restore after component energy checks failed");

  WVFieldEvaluationPlan fullPi,pointPi;
  require(bool(service->createPlan({full("pi")},fullPi)) &&
      bool(service->createPlan({{"point-pi","pi",points}},pointPi)),
      "F-bundle session plan creation failed");
  std::vector<double> fullPiValues(fullPi.outputs()[0].elementCount);
  std::vector<double> pointPiValues(pointPi.outputs()[0].elementCount);
  WVFieldOutputView fullPiView{fullPiValues.data(),fullPiValues.size()};
  WVFieldOutputView pointPiView{pointPiValues.data(),pointPiValues.size()};
  const auto fBefore=service->producerMetrics();
  const auto fLedgerBefore=service->metrics().variableEvaluation;
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "F-bundle evaluation session did not start");
    require(bool(service->evaluate(fullPi,state,&fullPiView,1)) &&
        bool(service->evaluate(first,state,&fullView,1)) &&
        bool(service->evaluate(pointPi,state,&pointPiView,1)),
        "pi/u/pi F-bundle evaluation session failed");
  }
  const auto fAfter=service->producerMetrics();
  const auto fLedgerAfter=service->metrics().variableEvaluation;
  require(fAfter.phasePreparations==fBefore.phasePreparations+1 &&
      fAfter.reconstructions[0][0][0]==fBefore.reconstructions[0][0][0]+1,
      "pi/u/pi session repeated phase or velocity reconstruction");
  for(std::size_t derivative=0;derivative<4;++derivative)
    require(fAfter.reconstructions[4][derivative][0]==
        fBefore.reconstructions[4][derivative][0]+1,
        "pi/u/pi session repeated an identified F-bundle producer");
  require(fLedgerAfter.producerExecutions==fLedgerBefore.producerExecutions+2 &&
      fLedgerAfter.cacheHits==fLedgerBefore.cacheHits+1 &&
      fLedgerAfter.duplicateExecutions==fLedgerBefore.duplicateExecutions,
      "pi/u/pi session did not reuse one complete fused F bundle");
}

void verifyActiveDemandPlanExtension() {
  const auto config=configuration(6,5,false,true);
  const auto owned=stateFor(config);
  const WVIntegrationState state{owned.view()};
  std::unique_ptr<WVFieldEvaluationService> service,reference;
  require(bool(WVFieldEvaluationService::create(config,
      std::make_unique<WVReferenceFFTEngine>(),service)),
      "active demand service creation failed");
  require(bool(WVFieldEvaluationService::create(config,
      std::make_unique<WVReferenceFFTEngine>(),reference)),
      "active demand reference creation failed");
  WVFieldEvaluationPlan uPlan,unionPlan;
  require(bool(service->createPlan({full("u")},uPlan)) &&
      bool(reference->createPlan({full("u"),full("zeta_x")},unionPlan)),
      "active demand initial plans failed");
  WVFieldEvaluationPlan unchanged=uPlan;
  auto status=service->createPlanForActiveEvaluation(
      {full("zeta_x")},unchanged);
  require(status.code==WVKernelStatusCode::invalidConfiguration &&
      unchanged.outputCount()==uPlan.outputCount() &&
      unchanged.outputs()[0].fieldName==uPlan.outputs()[0].fieldName,
      "active demand plan was accepted without a session or replaced its output");

  std::array<std::vector<double>,2> expected;
  std::array<WVFieldOutputView,2> expectedViews;
  for(std::size_t i=0;i<expected.size();++i) {
    expected[i].resize(unionPlan.outputs()[i].elementCount);
    expectedViews[i]={expected[i].data(),expected[i].size()};
  }
  require(bool(reference->evaluate(unionPlan,state,expectedViews.data(),2)),
      "active demand union reference failed");

  std::vector<double> u(uPlan.outputs()[0].elementCount);
  WVFieldOutputView uView{u.data(),u.size()};
  WVFieldEvaluationPlan zetaPlan;
  const auto producersBefore=service->producerMetrics();
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "active demand session failed");
    require(bool(service->evaluate(uPlan,state,&uView,1)),
        "active demand initial field failed");
    const auto retainedU=u;
    WVFieldEvaluationPlan failed=uPlan;
    status=service->createPlanForActiveEvaluation(
        {{"invalid-active","not_a_registered_variable",{}}},failed);
    require(!status && failed.outputCount()==uPlan.outputCount() && u==retainedU,
        "failed active demand changed its plan or published field copy");
    status=service->createPlanForActiveEvaluation({full("zeta_x")},zetaPlan);
    require(bool(status),status.message);
    require(u==retainedU,"active demand preparation changed an earlier output copy");
    std::vector<double> zeta(zetaPlan.outputs()[0].elementCount);
    WVFieldOutputView zetaView{zeta.data(),zeta.size()};
    require(bool(service->evaluate(zetaPlan,state,&zetaView,1)),
        "active demand derivative field failed");
    for(std::size_t i=0;i<u.size();++i) requireClose(u[i],expected[0][i],
        "active demand changed the initial field");
    for(std::size_t i=0;i<zeta.size();++i) requireClose(zeta[i],expected[1][i],
        "active demand derivative differs from preplanned union");
  }
  const auto producersAfter=service->producerMetrics();
  require(producersAfter.phasePreparations==producersBefore.phasePreparations+1 &&
      producersAfter.reconstructions[0][0][0]==
          producersBefore.reconstructions[0][0][0]+1,
      "active demand extension repeated the retained velocity producer");
  const auto retained=service->persistentBytes();
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "active demand second event failed");
    std::vector<double> zeta(zetaPlan.outputs()[0].elementCount);
    WVFieldOutputView zetaView{zeta.data(),zeta.size()};
    require(bool(service->evaluate(zetaPlan,state,&zetaView,1)),
        "active demand retained plan failed in the next event");
  }
  require(service->persistentBytes()==retained,
      "active demand retained capacity grew in the next event");

  require(bool(service->setVariableEvaluationPolicy(
      WVVariableEvaluationPolicy::lowMemory)),
      "active demand low-memory setup failed");
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "active demand low-memory session failed");
    WVFieldEvaluationPlan rejected=uPlan;
    status=service->createPlanForActiveEvaluation({full("v")},rejected);
    require(status.code==WVKernelStatusCode::invalidConfiguration &&
        rejected.outputCount()==uPlan.outputCount(),
        "active demand accepted low-memory preparation");
  }
}

void verifyActiveDemandReentrantRejection() {
  const auto config=configuration(6,5,true,true);
  const auto owned=stateFor(config);
  const WVIntegrationState state{owned.view()};
  const auto probe=std::make_shared<ActivePlanProbe>();
  std::unique_ptr<WVFieldEvaluationService> service;
  require(bool(WVFieldEvaluationService::create(config,
      std::make_unique<ActivePlanProbeEngine>(probe),service)),
      "active demand probe service failed");
  probe->service=service.get();
  WVFieldEvaluationPlan plan;
  require(bool(service->createPlan({full("u")},plan)),
      "active demand probe plan failed");
  std::vector<double> values(plan.outputs()[0].elementCount);
  WVFieldOutputView output{values.data(),values.size()};
  {
    WVFieldEvaluationSession session;
    require(bool(service->beginEvaluationSession(state,session)),
        "active demand probe session failed");
    probe->armed=true;
    require(bool(service->evaluate(plan,state,&output,1)),
        "active demand probe outer evaluation failed");
  }
  require(probe->attempted &&
      probe->status.code==WVKernelStatusCode::reentrantExecution,
      "active demand plan was not rejected inside a producer");
}

void verifyFusedVariableEvaluationLedger() {
  const WVVariableEvaluationKey a{WVVariableEvaluationNode::forcingTendency,1};
  const WVVariableEvaluationKey b{WVVariableEvaluationNode::forcingTendency,2};
  const WVVariableEvaluationKey c{WVVariableEvaluationNode::forcingTendency,3};
  WVVariableEvaluationContext context;
  require(bool(context.prepare({a,b,c})),"fused ledger preparation failed");
  require(bool(context.begin(&context)),"fused ledger scope failed");
  const std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> ab{{a,8},{b,16}};
  const std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> bc{{b,16},{c,8}};
  require(!context.evaluateGroup(ab,[&]() {
      return context.evaluate(a,8,[](){return WVKernelStatus::ok();});
    }) && !context.ready(a) && !context.ready(b),
      "fused ledger cycle did not roll back atomically");
  require(bool(context.evaluateGroup(ab,[](){return WVKernelStatus::ok();})) &&
      context.metrics().producerExecutions==1,
      "fused ledger did not count one actual producer");
  require(!context.evaluateGroup(bc,[](){return WVKernelStatus::ok();}),
      "fused ledger accepted overlapping ready nodes");
  context.end();
  require(bool(context.begin(&context,WVVariableEvaluationPolicy::lowMemory)),
      "low-memory fused ledger scope failed");
  require(bool(context.evaluateGroup(ab,[](){return WVKernelStatus::ok();})) &&
      context.evict(a) && !context.evaluateGroup(ab,[](){return WVKernelStatus::ok();}) &&
      context.evict(b) && bool(context.evaluateGroup(ab,[](){return WVKernelStatus::ok();})) &&
      context.metrics().recomputations==1,
      "low-memory fused ledger did not require complete eviction before recomputation");
  context.end();
}

} // namespace

int main() {
  try {
    for (const bool hydrostatic : {false, true})
      for (const bool antialias : {false, true})
        verifyPhaseDiagnostics(hydrostatic, antialias);
    verifyBarotropicSplineExtrapolation();
    verifySmallGridSplineBoundaries();
    verifyVariableEvaluationSessions();
    verifyActiveDemandPlanExtension();
    verifyActiveDemandReentrantRejection();
    verifyFusedVariableEvaluationLedger();
    verifyCatalog();
    verifyPlanValidation();
    verifyFailureAndLifecycleContracts();
    verifyLowMemoryComponentFailureCleanup();
    verifyEvaluation(6, 5, true, true);
    verifyEvaluation(7, 6, false, false);
    verifyDerivedMovingSampling();
    verifyDiagnosticSamplingRoutes();
    verifyEventFieldEvaluation();
    std::cout << "WVFieldEvaluationService portable contracts passed: "
                 "hydrostatic/nonhydrostatic, odd/even, antialiasing, "
                 "zero/Nyquist, wrapping, profiles, linear/spline, reuse, "
                 "event-variable geometry, every position field, storage, "
                 "lifecycle, and failures.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
