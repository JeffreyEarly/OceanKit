#include "WVBarotropicQGFieldEvaluationAdapter.hpp"

#include "WaveVortexRuntime/WVIntegrationState.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

std::size_t checkedProduct(std::size_t first, std::size_t second) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first)
    throw std::overflow_error("Barotropic QG field size overflow");
  return first * second;
}

double wrapped(double coordinate, double length) noexcept {
  double result = std::fmod(coordinate, length);
  if (result < 0.0)
    result += length;
  return result >= length ? 0.0 : result;
}

std::uint64_t configurationFingerprint(
    const WVTransformBarotropicQGConfiguration &configuration) noexcept {
  std::uint64_t result = 1469598103934665603ULL;
  const auto append = [&](const auto &value) {
    const auto *bytes = reinterpret_cast<const std::uint8_t *>(&value);
    for (std::size_t index = 0; index < sizeof(value); ++index) {
      result ^= bytes[index];
      result *= 1099511628211ULL;
    }
  };
  append(configuration.contractVersion);
  append(configuration.Nx);
  append(configuration.Ny);
  append(configuration.Lx);
  append(configuration.Ly);
  append(configuration.h);
  append(configuration.j);
  append(configuration.g);
  append(configuration.planetaryRadius);
  append(configuration.rotationRate);
  append(configuration.latitude);
  append(configuration.shouldAntialias);
  return result;
}

class SplineSystem final {
public:
  explicit SplineSystem(std::size_t count) : count_(count) {
    if (count < 2)
      throw std::invalid_argument("Spline interpolation requires two points.");
    if (count < 4)
      return;
    lu_.assign(checkedProduct(count, count), 0.0);
    pivots_.resize(count);
    const auto set = [&](std::size_t row, std::size_t column, double value) {
      lu_[column + count * row] = value;
    };
    set(0, 0, -1.0);
    set(0, 1, 2.0);
    set(0, 2, -1.0);
    for (std::size_t row = 1; row + 1 < count; ++row) {
      set(row, row - 1, 1.0);
      set(row, row, 4.0);
      set(row, row + 1, 1.0);
    }
    set(count - 1, count - 3, -1.0);
    set(count - 1, count - 2, 2.0);
    set(count - 1, count - 1, -1.0);
    std::vector<double> transpose(checkedProduct(count, count));
    for (std::size_t row = 0; row < count; ++row)
      for (std::size_t column = 0; column < count; ++column)
        transpose[column + count * row] = lu_[row + count * column];
    lu_.swap(transpose);
    factor();
  }

  void weightsInto(double spacing, double query, std::size_t circularShift,
                   std::vector<double> &result,
                   std::vector<double> &shifted,
                   std::vector<double> &rhs) const {
    const double normalized = query / spacing;
    std::fill(shifted.begin(), shifted.end(), 0.0);
    std::fill(rhs.begin(), rhs.end(), 0.0);
    std::fill(result.begin(), result.end(), 0.0);
    if (count_ == 2) {
      shifted[0] = 1.0 - normalized;
      shifted[1] = normalized;
    } else if (count_ == 3) {
      shifted[0] =
          (normalized - 1.0) * (normalized - 2.0) / 2.0;
      shifted[1] = -normalized * (normalized - 2.0);
      shifted[2] = normalized * (normalized - 1.0) / 2.0;
    } else {
      std::size_t interval = normalized <= 0.0
                                 ? 0
                                 : static_cast<std::size_t>(
                                       std::floor(normalized));
      interval = std::min(interval, count_ - 2);
      const double fraction = std::clamp(
          normalized - static_cast<double>(interval), 0.0, 1.0);
      const double first = 1.0 - fraction;
      const double second = fraction;
      rhs[interval] =
          (first * first * first - first) * spacing * spacing / 6.0;
      rhs[interval + 1] =
          (second * second * second - second) * spacing * spacing / 6.0;
      solve(rhs);
      shifted[interval] += first;
      shifted[interval + 1] += second;
      const double scale = 6.0 / (spacing * spacing);
      for (std::size_t row = 1; row + 1 < count_; ++row) {
        shifted[row - 1] += scale * rhs[row];
        shifted[row] -= 2.0 * scale * rhs[row];
        shifted[row + 1] += scale * rhs[row];
      }
    }
    for (std::size_t shiftedIndex = 0; shiftedIndex < count_; ++shiftedIndex) {
      const auto original =
          (shiftedIndex + count_ - circularShift % count_) % count_;
      result[original] += shifted[shiftedIndex];
    }
  }

  std::vector<double> weights(double spacing, double query,
                              std::size_t circularShift) const {
    std::vector<double> result(count_, 0.0);
    std::vector<double> shifted(count_, 0.0);
    std::vector<double> rhs(count_, 0.0);
    weightsInto(spacing, query, circularShift, result, shifted, rhs);
    return result;
  }

  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + lu_.capacity() * sizeof(double) +
           pivots_.capacity() * sizeof(std::size_t);
  }

private:
  void factor() {
    for (std::size_t column = 0; column < count_; ++column) {
      std::size_t pivot = column;
      double maximum = std::abs(lu_[column + count_ * column]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        const double value = std::abs(lu_[column + count_ * row]);
        if (value > maximum) {
          maximum = value;
          pivot = row;
        }
      }
      if (maximum == 0.0)
        throw std::invalid_argument("Spline interpolation matrix is singular.");
      pivots_[column] = pivot;
      if (pivot != column)
        for (std::size_t entry = 0; entry < count_; ++entry)
          std::swap(lu_[entry + count_ * column],
                    lu_[entry + count_ * pivot]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        lu_[column + count_ * row] /= lu_[column + count_ * column];
        const double multiplier = lu_[column + count_ * row];
        for (std::size_t entry = column + 1; entry < count_; ++entry)
          lu_[entry + count_ * row] -=
              multiplier * lu_[entry + count_ * column];
      }
    }
  }

  void solve(std::vector<double> &rhs) const {
    for (std::size_t column = 0; column < count_; ++column) {
      if (pivots_[column] != column)
        std::swap(rhs[column], rhs[pivots_[column]]);
      for (std::size_t row = column + 1; row < count_; ++row)
        rhs[row] -= lu_[column + count_ * row] * rhs[column];
    }
    for (std::size_t reverse = 0; reverse < count_; ++reverse) {
      const auto row = count_ - 1 - reverse;
      for (std::size_t column = row + 1; column < count_; ++column)
        rhs[row] -= lu_[column + count_ * row] * rhs[column];
      rhs[row] /= lu_[row + count_ * row];
    }
  }

  std::size_t count_ = 0;
  std::vector<double> lu_;
  std::vector<std::size_t> pivots_;
};

enum class ScalarField : std::uint8_t { none, energy, uvMax };

struct Weight {
  std::array<std::size_t, 2> xIndices{};
  std::array<std::size_t, 2> yIndices{};
  std::array<double, 2> xWeights{};
  std::array<double, 2> yWeights{};
  std::vector<double> xSpline;
  std::vector<double> ySpline;
};

struct Request {
  WVBarotropicQGField field = WVBarotropicQGField::u;
  ScalarField scalar = ScalarField::none;
  WVFieldSamplingKind sampling = WVFieldSamplingKind::fullGrid;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::vector<Weight> weights;
  std::size_t output = 0;
};

struct Plan {
  std::uint64_t fingerprint = 0;
  std::vector<Request> requests;
  std::size_t persistentBytes() const noexcept {
    std::size_t bytes = sizeof(*this) + requests.capacity() * sizeof(Request);
    for (const auto &request : requests) {
      bytes += request.weights.capacity() * sizeof(Weight);
      for (const auto &weight : request.weights)
        bytes += (weight.xSpline.capacity() + weight.ySpline.capacity()) *
                 sizeof(double);
    }
    return bytes;
  }
};

struct MovingRequest {
  WVBarotropicQGField field = WVBarotropicQGField::u;
  std::size_t offset = 0;
  std::size_t count = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::size_t output = 0;
};

struct MovingPlan {
  std::uint64_t fingerprint = 0;
  std::vector<MovingRequest> requests;
  std::size_t positionCount = 0;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + requests.capacity() * sizeof(MovingRequest);
  }
};

struct EventRequest {
  WVBarotropicQGField field = WVBarotropicQGField::u;
  std::size_t positionSet = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::size_t output = 0;
};

struct EventPlan {
  std::uint64_t fingerprint = 0;
  std::vector<EventRequest> requests;
  std::size_t positionSetCount = 0;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + requests.capacity() * sizeof(EventRequest);
  }
};

struct EventGeometry {
  std::shared_ptr<const Plan> plan;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + (plan ? plan->persistentBytes() : 0);
  }
};

WVKernelStatus resolveField(const std::string &name,
                            WVBarotropicQGField &field,
                            ScalarField &scalar) {
  scalar = ScalarField::none;
  if (name == "u")
    field = WVBarotropicQGField::u;
  else if (name == "v")
    field = WVBarotropicQGField::v;
  else if (name == "eta")
    field = WVBarotropicQGField::eta;
  else if (name == "pi")
    field = WVBarotropicQGField::pi;
  else if (name == "psi")
    field = WVBarotropicQGField::psi;
  else if (name == "qgpv")
    field = WVBarotropicQGField::qgpv;
  else if (name == "zeta_z")
    field = WVBarotropicQGField::zetaZ;
  else if (name == "ssh")
    field = WVBarotropicQGField::ssh;
  else if (name == "energy")
    scalar = ScalarField::energy;
  else if (name == "uvMax")
    scalar = ScalarField::uvMax;
  else
    return {WVKernelStatusCode::unsupportedOperation,
            "WVTransformBarotropicQG does not support field " + name + "."};
  return WVKernelStatus::ok();
}

Weight weightsFor(double x, double y, WVPositionInterpolation interpolation,
                  const WVTransformBarotropicQGConfiguration &configuration,
                  const SplineSystem *xSpline,
                  const SplineSystem *ySpline) {
  const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
  const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
  const double xWrapped = wrapped(x, configuration.Lx);
  const double yWrapped = wrapped(y, configuration.Ly);
  const auto xLower = std::min(
      static_cast<std::size_t>(std::floor(xWrapped / dx)),
      configuration.Nx - 1);
  const auto yLower = std::min(
      static_cast<std::size_t>(std::floor(yWrapped / dy)),
      configuration.Ny - 1);
  Weight result;
  if (interpolation == WVPositionInterpolation::linear) {
    result.xIndices = {xLower, (xLower + 1) % configuration.Nx};
    result.yIndices = {yLower, (yLower + 1) % configuration.Ny};
    const double xFraction =
        (xWrapped - static_cast<double>(xLower) * dx) / dx;
    const double yFraction =
        (yWrapped - static_cast<double>(yLower) * dy) / dy;
    result.xWeights = {1.0 - xFraction, xFraction};
    result.yWeights = {1.0 - yFraction, yFraction};
  } else {
    const bool xBoundary = configuration.Nx < 4 || xLower < 3 ||
                           xLower > configuration.Nx - 4;
    const bool yBoundary = configuration.Ny < 4 || yLower < 3 ||
                           yLower > configuration.Ny - 4;
    const std::size_t xShift = xBoundary ? 4 : 0;
    const std::size_t yShift = yBoundary ? 4 : 0;
    const double xQuery = xBoundary
                              ? wrapped(x + 4.0 * dx, configuration.Lx)
                              : xWrapped;
    const double yQuery = yBoundary
                              ? wrapped(y + 4.0 * dy, configuration.Ly)
                              : yWrapped;
    result.xSpline = xSpline->weights(dx, xQuery, xShift);
    result.ySpline = ySpline->weights(dy, yQuery, yShift);
  }
  return result;
}

double interpolate(const double *field, const Weight &weight,
                   WVPositionInterpolation interpolation,
                   std::size_t Nx, std::size_t Ny) noexcept {
  if (interpolation == WVPositionInterpolation::linear) {
    double value = 0.0;
    for (std::size_t y = 0; y < 2; ++y)
      for (std::size_t x = 0; x < 2; ++x)
        value += weight.xWeights[x] * weight.yWeights[y] *
                 field[weight.xIndices[x] + Nx * weight.yIndices[y]];
    return value;
  }
  double value = 0.0;
  for (std::size_t y = 0; y < Ny; ++y)
    for (std::size_t x = 0; x < Nx; ++x)
      value += weight.xSpline[x] * weight.ySpline[y] * field[x + Nx * y];
  return value;
}

WVKernelStatus coefficientView(const WVIntegrationState &state,
                               const WVTransformBarotropicQGKernel &kernel,
                               WVComplexConstView &A0) {
  if (state.coefficientFamilyCount != 1 ||
      state.coefficientFamilies == nullptr ||
      state.coefficientFamilies[0].layout == nullptr ||
      state.coefficientFamilies[0].layout->identifier != "A0" ||
      state.coefficientFamilies[0].layout->elementCount !=
          kernel.descriptor().Nkl() ||
      state.coefficientFamilies[0].data == nullptr)
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG field evaluation requires compact A0 state."};
  A0 = {state.coefficientFamilies[0].data,
        {1, kernel.descriptor().Nkl()}};
  return WVKernelStatus::ok();
}

} // namespace

struct WVBarotropicQGFieldEvaluationAdapter::MovingInterpolationWorkspace {
  MovingInterpolationWorkspace(std::size_t Nx, std::size_t Ny)
      : xSpline(Nx), ySpline(Ny), xWeights(Nx), yWeights(Ny),
        xShifted(Nx), yShifted(Ny), xRhs(Nx), yRhs(Ny) {}

  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) +
           xSpline.persistentBytes() + ySpline.persistentBytes() +
           (xWeights.capacity() + yWeights.capacity() +
            xShifted.capacity() + yShifted.capacity() + xRhs.capacity() +
            yRhs.capacity()) *
               sizeof(double);
  }

  std::size_t scratchBytes() const noexcept {
    return (xWeights.capacity() + yWeights.capacity() +
            xShifted.capacity() + yShifted.capacity() + xRhs.capacity() +
            yRhs.capacity()) *
           sizeof(double);
  }

  SplineSystem xSpline;
  SplineSystem ySpline;
  std::vector<double> xWeights;
  std::vector<double> yWeights;
  std::vector<double> xShifted;
  std::vector<double> yShifted;
  std::vector<double> xRhs;
  std::vector<double> yRhs;
};

WVBarotropicQGFieldEvaluationAdapter::~WVBarotropicQGFieldEvaluationAdapter() =
    default;

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter>(
        new WVBarotropicQGFieldEvaluationAdapter());
    auto status = WVTransformBarotropicQGKernel::create(
        configuration, std::move(engine), candidate->ownedKernel_);
    if (!status)
      return status;
    candidate->kernel_ = candidate->ownedKernel_.get();
    candidate->fieldScratch_.resize(
        candidate->kernel_->descriptor().spatialShape().elementCount());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            configuration.Nx, configuration.Ny);
    candidate->metrics_.transformPersistentBytes =
        candidate->kernel_->persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        candidate->fieldScratch_.capacity() * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate Barotropic QG field service."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::createBorrowing(
    WVTransformBarotropicQGKernel &kernel,
    std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVBarotropicQGFieldEvaluationAdapter>(
        new WVBarotropicQGFieldEvaluationAdapter());
    candidate->kernel_ = &kernel;
    candidate->fieldScratch_.resize(
        kernel.descriptor().spatialShape().elementCount());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            kernel.descriptor().configuration().Nx,
            kernel.descriptor().configuration().Ny);
    candidate->metrics_.transformPersistentBytes = kernel.persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        candidate->fieldScratch_.capacity() * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed Barotropic QG field service."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::createPlan(
    const std::vector<WVFieldRequest> &requests,
    WVFieldEvaluationPlan &plan) const {
  try {
    auto implementation = std::make_shared<Plan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    implementation->requests.reserve(requests.size());
    WVFieldEvaluationPlan candidate;
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    const bool needsSpline = std::any_of(
        requests.begin(), requests.end(), [](const auto &request) {
          return request.sampling.kind == WVFieldSamplingKind::positions &&
                 request.sampling.interpolation ==
                     WVPositionInterpolation::spline;
        });
    std::unique_ptr<SplineSystem> xSpline;
    std::unique_ptr<SplineSystem> ySpline;
    if (needsSpline) {
      xSpline = std::make_unique<SplineSystem>(configuration().Nx);
      ySpline = std::make_unique<SplineSystem>(configuration().Ny);
    }
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &input = requests[index];
      if (input.identifier.empty() || !identifiers.insert(input.identifier).second)
        return invalid("Barotropic QG field identifiers must be nonempty and unique.");
      Request request;
      auto status = resolveField(input.fieldName, request.field, request.scalar);
      if (!status)
        return status;
      request.sampling = input.sampling.kind;
      request.interpolation = input.sampling.interpolation;
      request.output = index;
      WVFieldOutputSpecification output{
          input.identifier, input.fieldName, input.sampling.kind, {}, 0};
      if (request.scalar != ScalarField::none) {
        if (input.sampling.kind != WVFieldSamplingKind::fullGrid)
          return invalid("Scalar Barotropic QG fields do not support position sampling.");
        output.elementCount = 1;
      } else if (input.sampling.kind == WVFieldSamplingKind::fullGrid) {
        output.dimensions = {configuration().Nx, configuration().Ny};
        output.elementCount = checkedProduct(configuration().Nx,
                                             configuration().Ny);
      } else if (input.sampling.kind == WVFieldSamplingKind::positions) {
        if (input.sampling.x.empty() ||
            input.sampling.x.size() != input.sampling.y.size())
          return invalid("Barotropic QG position arrays must be equal and nonempty.");
        if (input.sampling.interpolation != WVPositionInterpolation::linear &&
            input.sampling.interpolation != WVPositionInterpolation::spline)
          return invalid("Barotropic QG interpolation method is unsupported.");
        request.weights.reserve(input.sampling.x.size());
        for (std::size_t position = 0; position < input.sampling.x.size();
             ++position) {
          if (!std::isfinite(input.sampling.x[position]) ||
              !std::isfinite(input.sampling.y[position]))
            return invalid("Barotropic QG positions must be finite.");
          request.weights.push_back(weightsFor(
              input.sampling.x[position], input.sampling.y[position],
              input.sampling.interpolation, configuration(), xSpline.get(),
              ySpline.get()));
        }
        output.dimensions = {input.sampling.x.size()};
        output.elementCount = input.sampling.x.size();
      } else {
        return {WVKernelStatusCode::unsupportedOperation,
                "Barotropic QG does not support fixed vertical profiles."};
      }
      implementation->requests.push_back(std::move(request));
      candidate.outputs_.push_back(std::move(output));
    }
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Barotropic QG field plan."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::evaluate(
    const WVFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  const auto plan =
      std::static_pointer_cast<const Plan>(publicPlan.transformPlan_);
  if (!plan || plan->fingerprint != configurationFingerprint(configuration()))
    return invalid("Barotropic QG field plan belongs to another transform.");
  if (outputCount != publicPlan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG outputs do not match the resolved plan."};
  for (std::size_t output = 0; output < outputCount; ++output)
    if (outputs[output].data == nullptr ||
        outputs[output].elementCount !=
            publicPlan.outputs_[output].elementCount)
      return {WVKernelStatusCode::invalidShape,
              "A Barotropic QG output has the wrong shape."};
  WVComplexConstView A0;
  auto status = coefficientView(state, *kernel_, A0);
  if (!status)
    return status;
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Barotropic QG field evaluation is not reentrant."};
  executing_ = true;
  struct Guard {
    bool &executing;
    ~Guard() { executing = false; }
  } guard{executing_};

  std::array<bool, 8> evaluated{};
  const auto spatial = kernel_->descriptor().spatialShape();
  for (const auto &request : plan->requests) {
    if (request.scalar != ScalarField::none) {
      double value = 0.0;
      status = request.scalar == ScalarField::energy
                   ? kernel_->totalEnergy(A0, value)
                   : kernel_->uvMax(A0, value);
      if (!status)
        return status;
      outputs[request.output].data[0] = value;
      ++metrics_.outputElementWriteCount;
      continue;
    }
    const auto fieldIndex = static_cast<std::size_t>(request.field);
    if (!evaluated[fieldIndex]) {
      const auto before = kernel_->metrics().executionCount;
      WVRealView view{fieldScratch_.data(), spatial};
      status = kernel_->transformA0ToField(A0, request.field, view);
      if (!status)
        return status;
      metrics_.fftExecutionCount +=
          kernel_->metrics().executionCount - before;
      ++metrics_.transformCount;
      ++metrics_.primitiveFieldEvaluationCount;
      evaluated[fieldIndex] = true;
      for (const auto &destination : plan->requests) {
        if (destination.scalar != ScalarField::none ||
            destination.field != request.field)
          continue;
        auto &output = outputs[destination.output];
        if (destination.sampling == WVFieldSamplingKind::fullGrid) {
          std::copy(fieldScratch_.begin(), fieldScratch_.end(), output.data);
          ++metrics_.fullGridWriteCount;
        } else {
          for (std::size_t position = 0;
               position < destination.weights.size(); ++position)
            output.data[position] = interpolate(
                fieldScratch_.data(), destination.weights[position],
                destination.interpolation, configuration().Nx,
                configuration().Ny);
          if (destination.interpolation == WVPositionInterpolation::linear)
            metrics_.linearInterpolationCount += destination.weights.size();
          else
            metrics_.splineInterpolationCount += destination.weights.size();
        }
        metrics_.outputElementWriteCount += output.elementCount;
        if (&destination != &request)
          ++metrics_.primitiveFieldReuseCount;
      }
    }
  }
  ++metrics_.evaluationCount;
  metrics_.scratchHighWaterBytes = std::max(
      metrics_.scratchHighWaterBytes,
      fieldScratch_.size() * sizeof(double));
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::createMovingPlan(
    const std::vector<WVMovingFieldRequest> &requests,
    WVMovingFieldEvaluationPlan &plan) const {
  try {
    auto implementation = std::make_shared<MovingPlan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    WVMovingFieldEvaluationPlan candidate;
    std::set<std::string> identifiers;
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &input = requests[index];
      WVBarotropicQGField field;
      ScalarField scalar;
      auto status = resolveField(input.fieldName, field, scalar);
      if (!status || scalar != ScalarField::none)
        return status ? invalid("A scalar field cannot be sampled at moving positions.") : status;
      if (input.identifier.empty() || !identifiers.insert(input.identifier).second ||
          input.positionCount == 0 ||
          (input.interpolation != WVPositionInterpolation::linear &&
           input.interpolation != WVPositionInterpolation::spline) ||
          input.positionOffset > std::numeric_limits<std::size_t>::max() -
                                     input.positionCount)
        return invalid("Barotropic QG moving-field request is invalid.");
      implementation->positionCount =
          std::max(implementation->positionCount,
                   input.positionOffset + input.positionCount);
      implementation->requests.push_back(
          {field, input.positionOffset, input.positionCount,
           input.interpolation, index});
      candidate.outputs_.push_back(
          {input.identifier, input.fieldName, WVFieldSamplingKind::positions,
           {input.positionCount}, input.positionCount});
    }
    candidate.positionCount_ = implementation->positionCount;
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Barotropic QG moving-field plan."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::evaluateMoving(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state, WVMovingPositionView positions,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  return evaluateMovingImpl(publicPlan, state, nullptr, positions, outputs,
                            outputCount);
}

WVKernelStatus
WVBarotropicQGFieldEvaluationAdapter::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  return evaluateMovingImpl(publicPlan, state, &advectionFields, positions,
                            outputs, outputCount);
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::evaluateMovingImpl(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView *advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  const auto plan =
      std::static_pointer_cast<const MovingPlan>(publicPlan.transformPlan_);
  if (!plan || plan->fingerprint != configurationFingerprint(configuration()))
    return invalid("Barotropic QG moving plan belongs to another transform.");
  if (positions.positionCount != plan->positionCount ||
      (positions.positionCount != 0 &&
       (positions.x == nullptr || positions.y == nullptr)) ||
      outputCount != publicPlan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG moving positions or outputs have the wrong shape."};
  for (std::size_t output = 0; output < outputCount; ++output)
    if (outputs[output].data == nullptr ||
        outputs[output].elementCount !=
            publicPlan.outputs_[output].elementCount)
      return {WVKernelStatusCode::invalidShape,
              "A Barotropic QG moving-field output has the wrong shape."};
  for (std::size_t position = 0; position < positions.positionCount;
       ++position)
    if (!std::isfinite(positions.x[position]) ||
        !std::isfinite(positions.y[position]))
      return invalid("Barotropic QG moving positions must be finite.");
  WVComplexConstView A0;
  auto status = coefficientView(state, *kernel_, A0);
  if (!status)
    return status;
  const auto &configuration = this->configuration();
  if (advectionFields != nullptr) {
    if (advectionFields->data == nullptr ||
        advectionFields->shape.first != configuration.Nx ||
        advectionFields->shape.second != configuration.Ny ||
        advectionFields->shape.third != 1 ||
        advectionFields->shape.fourth != 2)
      return {WVKernelStatusCode::invalidShape,
              "Prepared Barotropic QG advection fields must have shape "
              "[Nx,Ny,1,2]."};
    if (std::any_of(plan->requests.begin(), plan->requests.end(),
                    [](const auto &request) {
                      return request.field != WVBarotropicQGField::u &&
                             request.field != WVBarotropicQGField::v;
                    }))
      return {WVKernelStatusCode::unsupportedOperation,
              "Prepared Barotropic QG advection fields support only u and v."};
  }
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Barotropic QG moving-field evaluation is not reentrant."};
  executing_ = true;
  struct Guard {
    bool &executing;
    ~Guard() { executing = false; }
  } guard{executing_};

  const auto spatial = kernel_->descriptor().spatialShape();
  const auto R = spatial.elementCount();
  const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
  const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
  auto &workspace = *movingInterpolation_;
  const double *sampleField = fieldScratch_.data();
  const auto sample = [&](double x, double y,
                          WVPositionInterpolation interpolation) {
    const double xWrapped = wrapped(x, configuration.Lx);
    const double yWrapped = wrapped(y, configuration.Ly);
    const auto xLower = std::min(
        static_cast<std::size_t>(std::floor(xWrapped / dx)),
        configuration.Nx - 1);
    const auto yLower = std::min(
        static_cast<std::size_t>(std::floor(yWrapped / dy)),
        configuration.Ny - 1);
    if (interpolation == WVPositionInterpolation::linear) {
      const double xFraction =
          (xWrapped - static_cast<double>(xLower) * dx) / dx;
      const double yFraction =
          (yWrapped - static_cast<double>(yLower) * dy) / dy;
      const auto xUpper = (xLower + 1) % configuration.Nx;
      const auto yUpper = (yLower + 1) % configuration.Ny;
      return (1.0 - yFraction) *
                 ((1.0 - xFraction) *
                      sampleField[xLower + configuration.Nx * yLower] +
                  xFraction *
                      sampleField[xUpper + configuration.Nx * yLower]) +
             yFraction *
                 ((1.0 - xFraction) *
                      sampleField[xLower + configuration.Nx * yUpper] +
                  xFraction *
                      sampleField[xUpper + configuration.Nx * yUpper]);
    }
    const bool xBoundary = configuration.Nx < 4 || xLower < 3 ||
                           xLower > configuration.Nx - 4;
    const bool yBoundary = configuration.Ny < 4 || yLower < 3 ||
                           yLower > configuration.Ny - 4;
    workspace.xSpline.weightsInto(
        dx, xBoundary ? wrapped(x + 4.0 * dx, configuration.Lx) : xWrapped,
        xBoundary ? 4 : 0, workspace.xWeights, workspace.xShifted,
        workspace.xRhs);
    workspace.ySpline.weightsInto(
        dy, yBoundary ? wrapped(y + 4.0 * dy, configuration.Ly) : yWrapped,
        yBoundary ? 4 : 0, workspace.yWeights, workspace.yShifted,
        workspace.yRhs);
    double value = 0.0;
    for (std::size_t iy = 0; iy < configuration.Ny; ++iy)
      for (std::size_t ix = 0; ix < configuration.Nx; ++ix)
        value += workspace.xWeights[ix] * workspace.yWeights[iy] *
                 sampleField[ix + configuration.Nx * iy];
    return value;
  };

  std::array<bool, 8> evaluated{};
  for (const auto &request : plan->requests) {
    const auto fieldIndex = static_cast<std::size_t>(request.field);
    if (evaluated[fieldIndex])
      continue;
    if (advectionFields == nullptr) {
      const auto before = kernel_->metrics().executionCount;
      WVRealView field{fieldScratch_.data(), spatial};
      status = kernel_->transformA0ToField(A0, request.field, field);
      if (!status)
        return status;
      metrics_.fftExecutionCount += kernel_->metrics().executionCount - before;
      ++metrics_.transformCount;
      ++metrics_.movingPrimitiveTransformCount;
      sampleField = fieldScratch_.data();
    } else {
      const auto channel = request.field == WVBarotropicQGField::u ? 0 : 1;
      sampleField = advectionFields->data + channel * R;
      ++metrics_.primitiveFieldReuseCount;
    }
    ++metrics_.primitiveFieldEvaluationCount;
    evaluated[fieldIndex] = true;
    bool firstDestination = true;
    for (const auto &destination : plan->requests) {
      if (destination.field != request.field)
        continue;
      auto &output = outputs[destination.output];
      for (std::size_t position = 0; position < destination.count; ++position)
        output.data[position] = sample(
            positions.x[destination.offset + position],
            positions.y[destination.offset + position],
            destination.interpolation);
      if (destination.interpolation == WVPositionInterpolation::linear)
        metrics_.linearInterpolationCount += destination.count;
      else
        metrics_.splineInterpolationCount += destination.count;
      metrics_.outputElementWriteCount += destination.count;
      if (!firstDestination)
        ++metrics_.primitiveFieldReuseCount;
      firstDestination = false;
    }
  }
  ++metrics_.evaluationCount;
  ++metrics_.movingEvaluationCount;
  metrics_.movingPositionCount += positions.positionCount;
  metrics_.scratchHighWaterBytes = std::max(
      metrics_.scratchHighWaterBytes,
      fieldScratch_.size() * sizeof(double) + workspace.scratchBytes());
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::createEventPlan(
    const std::vector<WVEventFieldRequest> &requests,
    WVEventFieldEvaluationPlan &plan) {
  try {
    auto implementation = std::make_shared<EventPlan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    WVEventFieldEvaluationPlan candidate;
    candidate.outputs_.reserve(requests.size());
    for (std::size_t index = 0; index < requests.size(); ++index) {
      WVBarotropicQGField field;
      ScalarField scalar;
      auto status = resolveField(requests[index].fieldName, field, scalar);
      if (!status)
        return status;
      if (scalar != ScalarField::none)
        return invalid("Scalar fields cannot use event position sampling.");
      implementation->positionSetCount = std::max(
          implementation->positionSetCount, requests[index].positionSetSlot + 1);
      implementation->requests.push_back(
          {field, requests[index].positionSetSlot,
           requests[index].interpolation, index});
      candidate.outputs_.push_back(
          {requests[index].identifier, requests[index].fieldName,
           WVPortableVariable::invalid, WVPortableNaturalRank::horizontal, 0,
           requests[index].positionSetSlot, requests[index].interpolation});
    }
    candidate.positionSetCount_ = implementation->positionSetCount;
    candidate.fingerprint_ = implementation->fingerprint;
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    ++metrics_.eventPlanCreationCount;
    metrics_.eventPlanFieldResolutionCount += requests.size();
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Barotropic QG event-field plan."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::prepareEventGeometry(
    const WVEventFieldEvaluationPlan &publicPlan,
    const WVEventPositionSetView *positionSets,
    std::size_t positionSetCount, WVPreparedFieldGeometry &geometry) {
  const auto eventPlan =
      std::static_pointer_cast<const EventPlan>(publicPlan.transformPlan_);
  if (!eventPlan || positionSetCount != eventPlan->positionSetCount ||
      (positionSetCount != 0 && positionSets == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Barotropic QG event geometry does not match its plan."};
  try {
    std::vector<WVFieldRequest> requests;
    requests.reserve(eventPlan->requests.size());
    WVPreparedFieldGeometry candidate;
    candidate.positionSets_.reserve(positionSetCount);
    candidate.positionCount_ = 0;
    for (std::size_t slot = 0; slot < positionSetCount; ++slot) {
      const auto &view = positionSets[slot];
      if (view.positionCount != 0 && (view.x == nullptr || view.y == nullptr))
        return {WVKernelStatusCode::invalidPointer,
                "Barotropic QG event coordinates are null."};
      WVPreparedFieldGeometry::PositionSet stored;
      stored.x = view.x;
      stored.y = view.y;
      stored.positionCount = view.positionCount;
      if (view.extentCount == 0)
        stored.extents = {view.positionCount};
      else
        stored.extents.assign(view.extents, view.extents + view.extentCount);
      std::size_t count = 1;
      for (const auto extent : stored.extents)
        count = checkedProduct(count, extent);
      if (count != view.positionCount)
        return {WVKernelStatusCode::invalidShape,
                "Barotropic QG event extents do not match coordinates."};
      candidate.positionCount_ += view.positionCount;
      candidate.borrowedCoordinateBytes_ +=
          2 * view.positionCount * sizeof(double);
      candidate.positionSets_.push_back(std::move(stored));
    }
    for (const auto &request : eventPlan->requests) {
      const auto &set = positionSets[request.positionSet];
      WVFieldSamplingRequest sampling;
      sampling.kind = WVFieldSamplingKind::positions;
      sampling.interpolation = request.interpolation;
      sampling.x.assign(set.x, set.x + set.positionCount);
      sampling.y.assign(set.y, set.y + set.positionCount);
      requests.push_back({publicPlan.outputs_[request.output].identifier,
                          publicPlan.outputs_[request.output].fieldName,
                          std::move(sampling)});
      candidate.outputs_.push_back(
          {request.output, request.positionSet,
           candidate.positionSets_[request.positionSet].extents,
           set.positionCount});
    }
    WVFieldEvaluationPlan evaluationPlan;
    auto status = createPlan(requests, evaluationPlan);
    if (!status)
      return status;
    auto geometryImpl = std::make_shared<EventGeometry>();
    geometryImpl->plan =
        std::static_pointer_cast<const Plan>(evaluationPlan.transformPlan_);
    candidate.fieldPlanFingerprint_ = publicPlan.fingerprint_;
    candidate.geometryFingerprint_ = publicPlan.fingerprint_;
    for (const auto &set : candidate.positionSets_) {
      const auto byteCount = set.positionCount * sizeof(double);
      for (const auto *data : {set.x, set.y}) {
        const auto *bytes = reinterpret_cast<const std::uint8_t *>(data);
        for (std::size_t byte = 0; byte < byteCount; ++byte) {
          candidate.geometryFingerprint_ ^= bytes[byte];
          candidate.geometryFingerprint_ *= 1099511628211ULL;
        }
      }
    }
    candidate.transformGeometry_ = geometryImpl;
    candidate.transformGeometryBytes_ = geometryImpl->persistentBytes();
    candidate.evaluationPlan_ = std::move(evaluationPlan);
    geometry = std::move(candidate);
    ++metrics_.eventGeometryPreparationCount;
    metrics_.eventPositionSetCount += positionSetCount;
    metrics_.eventPositionCount += geometry.positionCount_;
    metrics_.lastPreparedGeometryRetainedBytes = geometry.retainedBytes();
    metrics_.maximumPreparedGeometryRetainedBytes = std::max(
        metrics_.maximumPreparedGeometryRetainedBytes,
        metrics_.lastPreparedGeometryRetainedBytes);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate Barotropic QG event geometry."};
  }
}

WVKernelStatus WVBarotropicQGFieldEvaluationAdapter::evaluateEventBatch(
    const WVIntegrationState &state,
    const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entryCount != 0 && entries == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Barotropic QG event batch entries are null."};
  for (std::size_t entry = 0; entry < entryCount; ++entry) {
    if (entries[entry].plan == nullptr || entries[entry].geometry == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Barotropic QG event entry is incomplete."};
    const auto status = evaluate(
        entries[entry].geometry->evaluationPlan_, state,
        entries[entry].outputs, entries[entry].outputCount);
    if (!status)
      return status;
  }
  ++metrics_.eventBatchEvaluationCount;
  metrics_.eventBatchOccurrenceCount += entryCount;
  metrics_.eventEvaluationCount += entryCount;
  return WVKernelStatus::ok();
}

bool WVBarotropicQGFieldEvaluationAdapter::isCompatibleWith(
    const WVIntegrationStateLayout &layout) const noexcept {
  const auto &spatial = layout.spatialDimensions();
  return spatial == std::vector<std::size_t>{configuration().Nx,
                                             configuration().Ny} &&
         layout.coefficientFamilyCount() == 1 &&
         layout.coefficientFamilies()[0].identifier == "A0" &&
         layout.coefficientFamilies()[0].elementCount ==
             kernel_->descriptor().Nkl();
}

const WVTransformBarotropicQGConfiguration &
WVBarotropicQGFieldEvaluationAdapter::configuration() const noexcept {
  return kernel_->descriptor().configuration();
}

std::size_t
WVBarotropicQGFieldEvaluationAdapter::persistentBytes() const noexcept {
  return sizeof(*this) +
         (ownedKernel_ ? kernel_->persistentBytes() : 0) +
         fieldScratch_.capacity() * sizeof(double) +
         (movingInterpolation_ ? movingInterpolation_->persistentBytes() : 0);
}

} // namespace wavevortex::runtime::detail
