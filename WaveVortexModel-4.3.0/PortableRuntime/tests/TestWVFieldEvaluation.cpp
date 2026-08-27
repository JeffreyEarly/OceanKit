#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"

#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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

void verifyCatalog() {
  const std::vector<std::string> expected = {
      "u",       "v",         "w",       "eta",    "pi",
      "p",       "psi",       "qgpv",    "rho_e",  "rho_total",
      "rho_bar", "zeta_x",    "zeta_y",  "zeta_z", "ssu",
      "ssv",     "ssh",       "energy",  "uvMax",  "wMax"};
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

void verifyEvaluation(std::size_t nx, std::size_t ny, bool hydrostatic,
                      bool antialias) {
  const auto config = configuration(nx, ny, hydrostatic, antialias);
  std::unique_ptr<WVFieldEvaluationService> service;
  auto status = WVFieldEvaluationService::create(
      config, std::make_unique<WVReferenceFFTEngine>(), service);
  require(static_cast<bool>(status), "field service creation failed");

  std::vector<WVFieldRequest> requests;
  for (const auto &name : WVFieldEvaluationService::supportedFieldNames())
    requests.push_back(full(name));
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

} // namespace

int main() {
  try {
    verifyCatalog();
    verifyPlanValidation();
    verifyFailureAndLifecycleContracts();
    verifyEvaluation(6, 5, true, true);
    verifyEvaluation(7, 6, false, false);
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
