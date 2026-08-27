#include "WaveVortexRuntime/WVFieldEvaluationService.hpp"
#include "WaveVortexRuntime/WVIntegrationState.hpp"
#include "WVBarotropicQGFieldEvaluationAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime {
namespace {

enum Dependency : std::uint64_t {
  primitiveValues = 1ULL << 0,
  pressureHeight = 1ULL << 1,
  streamfunction = 1ULL << 2,
  potentialVorticity = 1ULL << 3,
  uDerivatives = 1ULL << 4,
  vDerivatives = 1ULL << 5,
  wDerivatives = 1ULL << 6,
  spectralEnergy = 1ULL << 7
};

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

constexpr std::uint64_t fingerprintOffset = 1469598103934665603ULL;
constexpr std::uint64_t fingerprintPrime = 1099511628211ULL;

void appendFingerprint(std::uint64_t &fingerprint, const void *data,
                       std::size_t byteCount) noexcept {
  const auto *bytes = static_cast<const std::uint8_t *>(data);
  for (std::size_t index = 0; index < byteCount; ++index) {
    fingerprint ^= bytes[index];
    fingerprint *= fingerprintPrime;
  }
}

template <typename Value>
void appendFingerprint(std::uint64_t &fingerprint,
                       const Value &value) noexcept {
  appendFingerprint(fingerprint, &value, sizeof(value));
}

void appendConfigurationFingerprint(
    std::uint64_t &fingerprint,
    const WVTransformConstantStratificationConfiguration &configuration)
    noexcept {
  appendFingerprint(fingerprint, configuration.contractVersion);
  appendFingerprint(fingerprint, configuration.Nx);
  appendFingerprint(fingerprint, configuration.Ny);
  appendFingerprint(fingerprint, configuration.Nz);
  appendFingerprint(fingerprint, configuration.Nj);
  appendFingerprint(fingerprint, configuration.Lx);
  appendFingerprint(fingerprint, configuration.Ly);
  appendFingerprint(fingerprint, configuration.Lz);
  appendFingerprint(fingerprint, configuration.N0);
  appendFingerprint(fingerprint, configuration.rho0);
  appendFingerprint(fingerprint, configuration.g);
  appendFingerprint(fingerprint, configuration.planetaryRadius);
  appendFingerprint(fingerprint, configuration.rotationRate);
  appendFingerprint(fingerprint, configuration.latitude);
  appendFingerprint(fingerprint, configuration.isHydrostatic);
  appendFingerprint(fingerprint, configuration.shouldAntialias);
}

std::size_t checkedProduct(std::size_t first, std::size_t second) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first)
    throw std::overflow_error("field-evaluation size overflow");
  return first * second;
}

bool memoryOverlaps(const void *first, std::size_t firstBytes,
                    const void *second, std::size_t secondBytes) {
  if (!first || !second || firstBytes == 0 || secondBytes == 0)
    return false;
  const auto firstAddress = reinterpret_cast<std::uintptr_t>(first);
  const auto secondAddress = reinterpret_cast<std::uintptr_t>(second);
  return firstAddress < secondAddress + secondBytes &&
         secondAddress < firstAddress + firstBytes;
}

WVComplex64 subtract(WVComplex64 first, WVComplex64 second) noexcept {
  return {first.real - second.real, first.imag - second.imag};
}

WVComplex64 multiply(WVComplex64 first, WVComplex64 second) noexcept {
  return {first.real * second.real - first.imag * second.imag,
          first.real * second.imag + first.imag * second.real};
}

WVComplex64 multiply(WVComplex64 value, double scale) noexcept {
  return {value.real * scale, value.imag * scale};
}

WVComplex64 conjugate(WVComplex64 value) noexcept {
  return {value.real, -value.imag};
}

double squaredMagnitude(WVComplex64 value) noexcept {
  return value.real * value.real + value.imag * value.imag;
}

double wrapped(double coordinate, double length) noexcept {
  double value = std::fmod(coordinate, length);
  if (value < 0.0)
    value += length;
  if (value >= length)
    value = 0.0;
  return value;
}

class SplineSystem final {
public:
  explicit SplineSystem(std::size_t count) : count_(count) {
    if (count < 2)
      throw std::invalid_argument(
          "Spline interpolation requires at least two points per axis.");
    if (count < 4)
      return;
    lu_.assign(checkedProduct(count, count), 0.0);
    pivots_.resize(count);
    // Factor the transpose of the uniform-grid not-a-knot system. The
    // common grid-spacing factor cancels between the system and its RHS.
    auto setSystem = [&](std::size_t row, std::size_t column, double value) {
      lu_[column + count * row] = value;
    };
    setSystem(0, 0, -1.0);
    setSystem(0, 1, 2.0);
    setSystem(0, 2, -1.0);
    for (std::size_t row = 1; row + 1 < count; ++row) {
      setSystem(row, row - 1, 1.0);
      setSystem(row, row, 4.0);
      setSystem(row, row + 1, 1.0);
    }
    setSystem(count - 1, count - 3, -1.0);
    setSystem(count - 1, count - 2, 2.0);
    setSystem(count - 1, count - 1, -1.0);
    std::vector<double> transpose(checkedProduct(count, count));
    for (std::size_t row = 0; row < count; ++row)
      for (std::size_t column = 0; column < count; ++column)
        transpose[column + count * row] = lu_[row + count * column];
    lu_.swap(transpose);
    factor();
  }

  std::vector<double> weights(double origin, double spacing, double query,
                              std::size_t circularShift = 0) const {
    std::vector<double> result;
    weightsInto(origin, spacing, query, result, circularShift);
    return result;
  }

  void weightsInto(double origin, double spacing, double query,
                   std::vector<double> &result,
                   std::size_t circularShift = 0,
                   std::vector<double> *rightHandSideWorkspace = nullptr,
                   std::vector<double> *shiftedWorkspace = nullptr) const {
    const double normalized = (query - origin) / spacing;
    auto restoreOriginalOrdering = [&](const std::vector<double> &shifted) {
      result.assign(count_, 0.0);
      for (std::size_t shiftedIndex = 0; shiftedIndex < count_;
           ++shiftedIndex) {
        const auto originalIndex =
            (shiftedIndex + count_ - circularShift % count_) % count_;
        result[originalIndex] += shifted[shiftedIndex];
      }
    };
    if (count_ == 2) {
      const std::vector<double> linear{1.0 - normalized, normalized};
      restoreOriginalOrdering(linear);
      return;
    }
    if (count_ == 3) {
      std::vector<double> quadratic = {
          (normalized - 1.0) * (normalized - 2.0) / 2.0,
          -normalized * (normalized - 2.0),
          normalized * (normalized - 1.0) / 2.0};
      restoreOriginalOrdering(quadratic);
      return;
    }
    std::size_t interval = normalized <= 0.0
                               ? 0
                               : static_cast<std::size_t>(std::floor(normalized));
    interval = std::min(interval, count_ - 2);
    const double fraction = std::clamp(normalized - static_cast<double>(interval),
                                       0.0, 1.0);
    const double first = 1.0 - fraction;
    const double second = fraction;
    std::vector<double> localRightHandSide;
    auto &rhs = rightHandSideWorkspace == nullptr ? localRightHandSide
                                                   : *rightHandSideWorkspace;
    rhs.assign(count_, 0.0);
    rhs[interval] = (first * first * first - first) * spacing * spacing / 6.0;
    rhs[interval + 1] =
        (second * second * second - second) * spacing * spacing / 6.0;
    solve(rhs);
    std::vector<double> localShifted;
    auto &shifted = shiftedWorkspace == nullptr ? localShifted
                                                : *shiftedWorkspace;
    shifted.assign(count_, 0.0);
    shifted[interval] += first;
    shifted[interval + 1] += second;
    const double rhsScale = 6.0 / (spacing * spacing);
    for (std::size_t row = 1; row + 1 < count_; ++row) {
      shifted[row - 1] += rhsScale * rhs[row];
      shifted[row] -= 2.0 * rhsScale * rhs[row];
      shifted[row + 1] += rhsScale * rhs[row];
    }
    restoreOriginalOrdering(shifted);
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
        const double candidate = std::abs(lu_[column + count_ * row]);
        if (candidate > maximum) {
          maximum = candidate;
          pivot = row;
        }
      }
      if (maximum == 0.0)
        throw std::invalid_argument("Spline interpolation system is singular.");
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

  void solve(std::vector<double> &rightHandSide) const {
    for (std::size_t column = 0; column < count_; ++column) {
      if (pivots_[column] != column)
        std::swap(rightHandSide[column], rightHandSide[pivots_[column]]);
      for (std::size_t row = column + 1; row < count_; ++row)
        rightHandSide[row] -=
            lu_[column + count_ * row] * rightHandSide[column];
    }
    for (std::size_t reverse = 0; reverse < count_; ++reverse) {
      const std::size_t row = count_ - 1 - reverse;
      for (std::size_t column = row + 1; column < count_; ++column)
        rightHandSide[row] -=
            lu_[column + count_ * row] * rightHandSide[column];
      rightHandSide[row] /= lu_[row + count_ * row];
    }
  }

  std::size_t count_ = 0;
  std::vector<double> lu_;
  std::vector<std::size_t> pivots_;
};

class ExecutionGuard final {
public:
  explicit ExecutionGuard(bool &executing)
      : executing_(executing), entered_(!executing) {
    if (entered_)
      executing_ = true;
  }
  ~ExecutionGuard() {
    if (entered_)
      executing_ = false;
  }
  bool entered() const noexcept { return entered_; }

private:
  bool &executing_;
  bool entered_ = false;
};

} // namespace

class WVFieldEvaluationService::MovingWorkspace final {
public:
  explicit MovingWorkspace(
      const WVTransformConstantStratificationConfiguration &configuration)
      : xSpline(configuration.Nx), ySpline(configuration.Ny),
        zSpline(configuration.Nz) {
    xWeights.reserve(configuration.Nx);
    yWeights.reserve(configuration.Ny);
    zWeights.reserve(configuration.Nz);
    xRightHandSide.reserve(configuration.Nx);
    yRightHandSide.reserve(configuration.Ny);
    zRightHandSide.reserve(configuration.Nz);
    xShifted.reserve(configuration.Nx);
    yShifted.reserve(configuration.Ny);
    zShifted.reserve(configuration.Nz);
  }
  SplineSystem xSpline;
  SplineSystem ySpline;
  SplineSystem zSpline;
  std::vector<double> xWeights;
  std::vector<double> yWeights;
  std::vector<double> zWeights;
  std::vector<double> xRightHandSide;
  std::vector<double> yRightHandSide;
  std::vector<double> zRightHandSide;
  std::vector<double> xShifted;
  std::vector<double> yShifted;
  std::vector<double> zShifted;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) +
           xSpline.persistentBytes() - sizeof(SplineSystem) +
           ySpline.persistentBytes() - sizeof(SplineSystem) +
           zSpline.persistentBytes() - sizeof(SplineSystem) +
           (xWeights.capacity() + yWeights.capacity() + zWeights.capacity() +
            xRightHandSide.capacity() + yRightHandSide.capacity() +
            zRightHandSide.capacity() + xShifted.capacity() +
            yShifted.capacity() + zShifted.capacity()) *
               sizeof(double);
  }
};

std::size_t WVMovingFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) +
                      requests_.capacity() * sizeof(ResolvedRequest) +
                      outputs_.capacity() * sizeof(WVFieldOutputSpecification);
  for (const auto &output : outputs_)
    bytes += output.identifier.capacity() + output.fieldName.capacity() +
             output.dimensions.capacity() * sizeof(std::size_t);
  return bytes + transformPlanBytes_;
}

std::size_t WVEventFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t bytes =
      sizeof(*this) + requests_.capacity() * sizeof(ResolvedRequest) +
      outputs_.capacity() * sizeof(WVEventFieldOutputSpecification) +
      requiresZByPositionSet_.capacity() * sizeof(std::uint8_t);
  for (const auto &output : outputs_)
    bytes += output.identifier.capacity() + output.fieldName.capacity();
  return bytes + transformPlanBytes_;
}

std::size_t WVFieldEvaluationPlan::PositionWeights::persistentBytes() const
    noexcept {
  return xSplineWeights.capacity() * sizeof(double) +
         ySplineWeights.capacity() * sizeof(double) +
         zSplineWeights.capacity() * sizeof(double);
}

std::size_t WVFieldEvaluationPlan::ResolvedRequest::persistentBytes() const
    noexcept {
  std::size_t value = profileXIndices.capacity() * sizeof(std::size_t) +
                      profileYIndices.capacity() * sizeof(std::size_t) +
                      positionWeights.capacity() * sizeof(PositionWeights);
  for (const auto &weights : positionWeights)
    value += weights.persistentBytes();
  return value;
}

std::size_t WVFieldEvaluationPlan::persistentBytes() const noexcept {
  std::size_t value = sizeof(*this) +
                      requests_.capacity() * sizeof(ResolvedRequest) +
                      outputs_.capacity() * sizeof(WVFieldOutputSpecification);
  for (const auto &request : requests_)
    value += request.persistentBytes();
  for (const auto &output : outputs_) {
    value += output.identifier.capacity() + output.fieldName.capacity() +
             output.dimensions.capacity() * sizeof(std::size_t);
  }
  return value + transformPlanBytes_;
}

WVEventPositionSetView
WVPreparedFieldGeometry::positionSet(std::size_t slot) const noexcept {
  if (slot >= positionSets_.size())
    return {};
  const auto &set = positionSets_[slot];
  return {set.x, set.y, set.z, set.positionCount, set.extents.data(),
          set.extents.size()};
}

bool WVPreparedFieldGeometry::sameGeometry(
    const WVPreparedFieldGeometry &other) const noexcept {
  if (fieldPlanFingerprint_ != other.fieldPlanFingerprint_ ||
      geometryFingerprint_ != other.geometryFingerprint_ ||
      positionSets_.size() != other.positionSets_.size())
    return false;
  for (std::size_t slot = 0; slot < positionSets_.size(); ++slot) {
    const auto &first = positionSets_[slot];
    const auto &second = other.positionSets_[slot];
    if (first.positionCount != second.positionCount ||
        first.extents != second.extents ||
        (first.z == nullptr) != (second.z == nullptr))
      return false;
    const auto bytes = first.positionCount * sizeof(double);
    if (bytes != 0 &&
        (std::memcmp(first.x, second.x, bytes) != 0 ||
         std::memcmp(first.y, second.y, bytes) != 0 ||
         (first.z != nullptr &&
          std::memcmp(first.z, second.z, bytes) != 0)))
      return false;
  }
  return true;
}

std::size_t WVPreparedFieldGeometry::retainedBytes() const noexcept {
  std::size_t bytes =
      sizeof(*this) + positionSets_.capacity() * sizeof(PositionSet) +
      outputs_.capacity() * sizeof(WVPreparedFieldOutputSpecification);
  for (const auto &set : positionSets_)
    bytes += set.extents.capacity() * sizeof(std::size_t);
  for (const auto &output : outputs_)
    bytes += output.dimensions.capacity() * sizeof(std::size_t);
  const auto evaluationBytes = evaluationPlan_.persistentBytes();
  if (evaluationBytes >= sizeof(evaluationPlan_))
    bytes += evaluationBytes - sizeof(evaluationPlan_);
  return bytes + transformGeometryBytes_;
}

std::size_t WVPreparedFieldGeometry::liveBytes() const noexcept {
  return retainedBytes() + borrowedCoordinateBytes_;
}

WVPreparedFieldGeometryMetrics
WVPreparedFieldGeometry::metrics() const noexcept {
  return {positionSets_.size(), positionCount_, retainedBytes(), liveBytes()};
}

WVKernelStatus WVFieldEvaluationService::create(
    const WVTransformConstantStratificationConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status = WVTransformConstantStratificationKernel::create(
        configuration, std::move(engine), candidate->ownedTransform_);
    if (!status)
      return status;
    candidate->transform_ = candidate->ownedTransform_.get();
    status = candidate->initializeScratch();
    if (!status)
      return status;
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate bounded field-evaluation storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformConstantStratificationKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    candidate->transform_ = &transform;
    const auto status = candidate->initializeScratch();
    if (!status)
      return status;
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed field-evaluation storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  }
}

WVKernelStatus WVFieldEvaluationService::create(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status = detail::WVBarotropicQGFieldEvaluationAdapter::create(
        configuration, std::move(engine), candidate->barotropicQG_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->barotropicQG_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the Barotropic QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::createBorrowing(
    WVTransformBarotropicQGKernel &transform,
    std::unique_ptr<WVFieldEvaluationService> &service) {
  service.reset();
  try {
    auto candidate = std::unique_ptr<WVFieldEvaluationService>(
        new WVFieldEvaluationService());
    auto status =
        detail::WVBarotropicQGFieldEvaluationAdapter::createBorrowing(
            transform, candidate->barotropicQG_);
    if (!status)
      return status;
    candidate->metrics_ = candidate->barotropicQG_->metrics();
    service = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate the borrowed Barotropic QG field boundary."};
  }
}

WVKernelStatus WVFieldEvaluationService::initializeScratch() {
  const auto fieldElements = transform_->descriptor().spatialShape().elementCount();
  const auto coefficientElements =
      transform_->descriptor().spectralShape().elementCount();
  realScratch_.resize(checkedProduct(6, fieldElements));
  complexScratch_.resize(checkedProduct(2, coefficientElements));
  movingWorkspace_ = std::make_unique<MovingWorkspace>(
      transform_->descriptor().configuration());
  metrics_.transformPersistentBytes = transform_->persistentBytes();
  metrics_.scratchCapacityBytes =
      realScratch_.capacity() * sizeof(double) +
      complexScratch_.capacity() * sizeof(WVComplex64);
  metrics_.movingInterpolationWorkspaceBytes =
      movingWorkspace_->persistentBytes();
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVFieldEvaluationService::~WVFieldEvaluationService() = default;

std::vector<std::string> WVFieldEvaluationService::supportedFieldNames() {
  std::vector<std::string> result;
  result.reserve(WVPortableVariableCatalog.size());
  for (const auto &variable : WVPortableVariableCatalog)
    if (variable.kind == WVPortableVariableKind::field)
      result.emplace_back(variable.name);
  return result;
}

WVKernelStatus WVFieldEvaluationService::createPlan(
    const std::vector<WVFieldRequest> &requests,
    WVFieldEvaluationPlan &plan) const {
  if (barotropicQG_)
    return barotropicQG_->createPlan(requests, plan);
  try {
    WVFieldEvaluationPlan candidate;
    const auto &configuration = transform_->descriptor().configuration();
    candidate.configuration_ = configuration;
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    bool needsSpline = false;
    for (const auto &request : requests)
      needsSpline = needsSpline ||
                    (request.sampling.kind == WVFieldSamplingKind::positions &&
                     request.sampling.interpolation ==
                         WVPositionInterpolation::spline);
    std::unique_ptr<SplineSystem> xSpline;
    std::unique_ptr<SplineSystem> ySpline;
    std::unique_ptr<SplineSystem> zSpline;
    if (needsSpline) {
      xSpline = std::make_unique<SplineSystem>(configuration.Nx);
      ySpline = std::make_unique<SplineSystem>(configuration.Ny);
    }
    const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
    const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
    const double dz = configuration.Lz /
                      static_cast<double>(configuration.Nz - 1);
    for (std::size_t outputIndex = 0; outputIndex < requests.size();
         ++outputIndex) {
      const auto &request = requests[outputIndex];
      if (request.identifier.empty())
        return invalid("Every field request must have a nonempty identifier.");
      if (!identifiers.insert(request.identifier).second)
        return invalid("Field request identifiers must be unique: " +
                       request.identifier + ".");
      const auto *metadata = findPortableVariable(request.fieldName);
      if (metadata == nullptr || metadata->kind != WVPortableVariableKind::field)
        return invalid("Unknown or unsupported field: " + request.fieldName + ".");
      const auto fieldIndex = static_cast<std::size_t>(metadata->ordinal);
      const auto field = metadata->identifier;
      const auto rank = metadata->naturalRank;

      WVFieldEvaluationPlan::ResolvedRequest resolved;
      if (request.sampling.kind != WVFieldSamplingKind::fullGrid &&
          request.sampling.kind !=
              WVFieldSamplingKind::fixedVerticalProfiles &&
          request.sampling.kind != WVFieldSamplingKind::positions)
        return invalid("Unknown field-sampling kind.");
      const std::uint8_t requestedSampling =
          request.sampling.kind == WVFieldSamplingKind::fullGrid
              ? portableFullGridSampling
              : request.sampling.kind ==
                        WVFieldSamplingKind::fixedVerticalProfiles
                    ? portableFixedVerticalProfileSampling
                    : portablePositionSampling;
      if ((metadata->samplingMask & requestedSampling) == 0)
        return invalid("Sampling mode is unsupported for field " +
                       request.fieldName + ".");
      if (request.sampling.kind == WVFieldSamplingKind::positions &&
          request.sampling.interpolation != WVPositionInterpolation::linear &&
          request.sampling.interpolation != WVPositionInterpolation::spline)
        return invalid("Unknown position-interpolation method.");
      resolved.field = field;
      resolved.nativeRank = rank;
      resolved.samplingKind = request.sampling.kind;
      resolved.interpolation = request.sampling.interpolation;
      resolved.outputIndex = outputIndex;
      WVFieldOutputSpecification output;
      output.identifier = request.identifier;
      output.fieldName = request.fieldName;
      output.samplingKind = request.sampling.kind;
      if (request.sampling.kind == WVFieldSamplingKind::fullGrid) {
        if (rank == WVFieldEvaluationPlan::NativeRank::volume)
          output.dimensions = {configuration.Nx, configuration.Ny,
                               configuration.Nz};
        else if (rank == WVFieldEvaluationPlan::NativeRank::horizontal)
          output.dimensions = {configuration.Nx, configuration.Ny};
        else if (rank == WVFieldEvaluationPlan::NativeRank::vertical)
          output.dimensions = {configuration.Nz};
      } else if (request.sampling.kind ==
                 WVFieldSamplingKind::fixedVerticalProfiles) {
        if (rank != WVFieldEvaluationPlan::NativeRank::volume)
          return invalid("Fixed vertical profiles require a three-dimensional field: " +
                         request.fieldName + ".");
        if (request.sampling.xIndices.empty() ||
            request.sampling.xIndices.size() !=
                request.sampling.yIndices.size())
          return invalid("Fixed-profile xIndices and yIndices must be nonempty and have equal length.");
        resolved.profileXIndices.reserve(request.sampling.xIndices.size());
        resolved.profileYIndices.reserve(request.sampling.yIndices.size());
        for (std::size_t index = 0; index < request.sampling.xIndices.size();
             ++index) {
          const auto xIndex = request.sampling.xIndices[index];
          const auto yIndex = request.sampling.yIndices[index];
          if (xIndex == 0 || xIndex > configuration.Nx || yIndex == 0 ||
              yIndex > configuration.Ny)
            return invalid("Fixed-profile indices must use MATLAB one-based values within the model grid.");
          resolved.profileXIndices.push_back(xIndex - 1);
          resolved.profileYIndices.push_back(yIndex - 1);
        }
        output.dimensions = {configuration.Nz,
                             request.sampling.xIndices.size()};
      } else {
        if (rank == WVFieldEvaluationPlan::NativeRank::scalar ||
            rank == WVFieldEvaluationPlan::NativeRank::vertical)
          return invalid("Position sampling is unsupported for field " +
                         request.fieldName + ".");
        const auto positionCount = request.sampling.x.size();
        if (positionCount == 0 || request.sampling.y.size() != positionCount)
          return invalid("Position x and y arrays must be nonempty and have equal length.");
        if (rank == WVFieldEvaluationPlan::NativeRank::volume &&
            request.sampling.z.size() != positionCount)
          return invalid("Three-dimensional position sampling requires equal-length x, y, and z arrays.");
        if (rank == WVFieldEvaluationPlan::NativeRank::horizontal &&
            !request.sampling.z.empty() &&
            request.sampling.z.size() != positionCount)
          return invalid("Optional z positions must be empty or match x and y.");
        resolved.positionWeights.reserve(positionCount);
        for (std::size_t position = 0; position < positionCount; ++position) {
          const double x = request.sampling.x[position];
          const double y = request.sampling.y[position];
          const double z = rank == WVFieldEvaluationPlan::NativeRank::volume
                               ? request.sampling.z[position]
                               : 0.0;
          if (!std::isfinite(x) || !std::isfinite(y) ||
              (rank == WVFieldEvaluationPlan::NativeRank::volume &&
               !std::isfinite(z)))
            return invalid("Position coordinates must be finite.");
          WVFieldEvaluationPlan::PositionWeights weights;
          const double xWrapped = wrapped(x, configuration.Lx);
          const double yWrapped = wrapped(y, configuration.Ly);
          const auto xLower = std::min(
              static_cast<std::size_t>(std::floor(xWrapped / dx)),
              configuration.Nx - 1);
          const auto yLower = std::min(
              static_cast<std::size_t>(std::floor(yWrapped / dy)),
              configuration.Ny - 1);
          if (request.sampling.interpolation ==
              WVPositionInterpolation::linear) {
            weights.xLinearIndices = {xLower, (xLower + 1) % configuration.Nx};
            weights.yLinearIndices = {yLower, (yLower + 1) % configuration.Ny};
            const double xFraction =
                (xWrapped - static_cast<double>(xLower) * dx) / dx;
            const double yFraction =
                (yWrapped - static_cast<double>(yLower) * dy) / dy;
            weights.xLinearWeights = {1.0 - xFraction, xFraction};
            weights.yLinearWeights = {1.0 - yFraction, yFraction};
            if (rank == WVFieldEvaluationPlan::NativeRank::volume) {
              weights.outsideInterpolationDomain =
                  z < -configuration.Lz || z > 0.0;
              if (!weights.outsideInterpolationDomain) {
                const double normalizedZ = (z + configuration.Lz) / dz;
                const auto zLower = std::min(
                    static_cast<std::size_t>(std::max(0.0, std::floor(normalizedZ))),
                    configuration.Nz - 2);
                const double zFraction =
                    std::clamp(normalizedZ - static_cast<double>(zLower), 0.0,
                               1.0);
                weights.zLinearIndices = {zLower, zLower + 1};
                weights.zLinearWeights = {1.0 - zFraction, zFraction};
              }
            }
          } else {
            const bool xBoundary = xLower < 3 || xLower > configuration.Nx - 4;
            const bool yBoundary = yLower < 3 || yLower > configuration.Ny - 4;
            const std::size_t xShift = xBoundary ? 4 : 0;
            const std::size_t yShift = yBoundary ? 4 : 0;
            const double xQuery = xBoundary
                                      ? wrapped(x + 4.0 * dx, configuration.Lx)
                                      : xWrapped;
            const double yQuery = yBoundary
                                      ? wrapped(y + 4.0 * dy, configuration.Ly)
                                      : yWrapped;
            // MATLAB's interpn(...,"spline",0) applies the zero fill value
            // when a boundary-shifted query lies beyond the last stored
            // periodic grid point (the duplicated x=L/y=L points are absent).
            weights.outsideInterpolationDomain =
                xQuery > static_cast<double>(configuration.Nx - 1) * dx ||
                yQuery > static_cast<double>(configuration.Ny - 1) * dy;
            if (!weights.outsideInterpolationDomain) {
              weights.xSplineWeights =
                  xSpline->weights(0.0, dx, xQuery, xShift);
              weights.ySplineWeights =
                  ySpline->weights(0.0, dy, yQuery, yShift);
            }
            if (rank == WVFieldEvaluationPlan::NativeRank::volume) {
              weights.outsideInterpolationDomain =
                  weights.outsideInterpolationDomain ||
                  z < -configuration.Lz || z > 0.0;
              if (!weights.outsideInterpolationDomain) {
                if (!zSpline)
                  zSpline = std::make_unique<SplineSystem>(configuration.Nz);
                weights.zSplineWeights =
                    zSpline->weights(-configuration.Lz, dz, z);
              }
            }
          }
          resolved.positionWeights.push_back(std::move(weights));
        }
        output.dimensions = {positionCount};
      }
      output.elementCount = 1;
      for (const auto dimension : output.dimensions)
        output.elementCount = checkedProduct(output.elementCount, dimension);
      candidate.requestedFieldMask_ |= 1ULL << fieldIndex;
      if (field == WVFieldEvaluationPlan::Field::psi &&
          transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Streamfunction evaluation is undefined when the Coriolis "
                "frequency is zero."};
      candidate.dependencyMask_ |= metadata->primitiveDependencyMask;
      candidate.requests_.push_back(std::move(resolved));
      candidate.outputs_.push_back(std::move(output));
    }
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate immutable field-evaluation plan storage."};
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::invalid_argument &error) {
    return invalid(error.what());
  }
}

WVKernelStatus WVFieldEvaluationService::evaluate(
    const WVFieldEvaluationPlan &plan, const WVState &state,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  const PlanInvocation invocation{&plan, outputs, outputCount};
  return evaluatePlanBatch(&invocation, 1, state);
}

WVKernelStatus WVFieldEvaluationService::evaluate(
    const WVFieldEvaluationPlan &plan, const WVIntegrationState &state,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  if (barotropicQG_)
    return barotropicQG_->evaluate(plan, state, outputs, outputCount);
  return evaluate(plan, state.waveVortex, outputs, outputCount);
}

WVKernelStatus
WVFieldEvaluationService::evaluatePlanBatch(const PlanInvocation *invocations,
                                            std::size_t invocationCount,
                                            const WVState &state) {
  if (invocationCount == 0)
    return WVKernelStatus::ok();
  if (invocations == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Field-evaluation plan invocations have a null pointer."};

  std::size_t totalOutputCount = 0;
  std::size_t planBytes = 0;
  std::uint64_t requestedFieldMask = 0;
  std::uint64_t dependencyMask = 0;
  bool allOutputsEmpty = true;
  for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
       ++invocationIndex) {
    const auto &invocation = invocations[invocationIndex];
    if (invocation.plan == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation invocation has a null plan."};
    const auto &plan = *invocation.plan;
    if (!sameTransformConfiguration(plan.configuration_,
                                    transform_->descriptor().configuration()))
      return invalid("The evaluation plan was created for a different "
                     "transform configuration.");
    if (invocation.outputCount != plan.outputs_.size())
      return {WVKernelStatusCode::invalidShape,
              "The output-view count must match the evaluation plan."};
    if (invocation.outputCount != 0 && invocation.outputs == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation outputs have a null view pointer."};
    if (totalOutputCount >
        std::numeric_limits<std::size_t>::max() - invocation.outputCount)
      return {WVKernelStatusCode::sizeOverflow,
              "Field-evaluation batch output count overflows size_t."};
    totalOutputCount += invocation.outputCount;
    const auto invocationPlanBytes = plan.persistentBytes();
    planBytes = planBytes > std::numeric_limits<std::size_t>::max() -
                                invocationPlanBytes
                    ? std::numeric_limits<std::size_t>::max()
                    : planBytes + invocationPlanBytes;
    requestedFieldMask |= plan.requestedFieldMask_;
    dependencyMask |= plan.dependencyMask_;
    for (const auto &output : plan.outputs_)
      allOutputsEmpty = allOutputsEmpty && output.elementCount == 0;
  }
  if (!std::isfinite(state.t) || !std::isfinite(state.t0))
    return invalid("Field-evaluation state times must be finite.");
  const auto spectral = transform_->descriptor().spectralShape();
  const auto coefficientBytes = spectral.elementCount() * sizeof(WVComplex64);
  const void *coefficientInputs[] = {state.coefficients.Ap.data,
                                     state.coefficients.Am.data,
                                     state.coefficients.A0.data};
  const WVComplexConstView coefficientViews[] = {
      state.coefficients.Ap, state.coefficients.Am, state.coefficients.A0};
  for (const auto &view : coefficientViews) {
    if (view.shape.rows != spectral.rows ||
        view.shape.columns != spectral.columns)
      return {WVKernelStatusCode::invalidShape,
              "Field-evaluation coefficients must have shape [Nj,Nkl]."};
    if (view.data == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Field-evaluation coefficients have a null pointer."};
  }
  for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
       ++invocationIndex) {
    const auto &invocation = invocations[invocationIndex];
    const auto &plan = *invocation.plan;
    for (std::size_t outputIndex = 0; outputIndex < invocation.outputCount;
         ++outputIndex) {
      const auto &output = invocation.outputs[outputIndex];
      if (output.elementCount != plan.outputs_[outputIndex].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Caller-owned output has the wrong element count for request " +
                    plan.outputs_[outputIndex].identifier + "."};
      if (output.elementCount != 0 && output.data == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Caller-owned output has a null pointer for request " +
                    plan.outputs_[outputIndex].identifier + "."};
      const auto bytes = output.elementCount * sizeof(double);
      for (const auto *input : coefficientInputs)
        if (memoryOverlaps(output.data, bytes, input, coefficientBytes))
          return {WVKernelStatusCode::overlappingArrays,
                  "Field outputs must not overlap coefficient inputs."};
      for (std::size_t otherInvocation = invocationIndex;
           otherInvocation < invocationCount; ++otherInvocation) {
        const auto &other = invocations[otherInvocation];
        const std::size_t firstOther =
            otherInvocation == invocationIndex ? outputIndex + 1 : 0;
        for (std::size_t otherOutput = firstOther;
             otherOutput < other.outputCount; ++otherOutput)
          if (memoryOverlaps(
                  output.data, bytes, other.outputs[otherOutput].data,
                  other.outputs[otherOutput].elementCount * sizeof(double)))
            return {WVKernelStatusCode::overlappingArrays,
                    "Caller-owned field outputs must not overlap each other."};
      }
    }
  }
  ExecutionGuard guard(executing_);
  if (!guard.entered())
    return {WVKernelStatusCode::reentrantExecution,
            "Field evaluation is not reentrant."};

  ++metrics_.evaluationCount;
  if (totalOutputCount > 1)
    ++metrics_.coincidentBatchCount;
  metrics_.lastPlanBytes = planBytes;
  metrics_.maximumPlanBytes =
      std::max(metrics_.maximumPlanBytes, metrics_.lastPlanBytes);

  if (allOutputsEmpty)
    return WVKernelStatus::ok();

  const auto &configuration = transform_->descriptor().configuration();
  const auto spatial = transform_->descriptor().spatialShape();
  const auto fieldElements = spatial.elementCount();
  const auto horizontalElements = configuration.Nx * configuration.Ny;
  const auto coefficientElements = spectral.elementCount();
  auto updateScratchHighWater = [&](std::size_t realElements,
                                    std::size_t complexElements) {
    const auto bytes =
        realElements * sizeof(double) + complexElements * sizeof(WVComplex64);
    metrics_.scratchHighWaterBytes =
        std::max(metrics_.scratchHighWaterBytes, bytes);
  };
  auto invokeTransform = [&](auto &&operation) {
    const auto before = transform_->metrics().executionCount;
    ++metrics_.transformCount;
    const auto status = operation();
    metrics_.fftExecutionCount += transform_->metrics().executionCount - before;
    return status;
  };

  auto fieldRequested = [&](WVFieldEvaluationPlan::Field field) {
    return (requestedFieldMask & (1ULL << static_cast<std::size_t>(field))) !=
           0;
  };

  auto writeField = [&](WVFieldEvaluationPlan::Field field,
                        WVFieldEvaluationPlan::NativeRank rank,
                        const double *source) -> WVKernelStatus {
    std::size_t sourceElements = 1;
    if (rank == WVFieldEvaluationPlan::NativeRank::volume)
      sourceElements = fieldElements;
    else if (rank == WVFieldEvaluationPlan::NativeRank::horizontal)
      sourceElements = horizontalElements;
    else if (rank == WVFieldEvaluationPlan::NativeRank::vertical)
      sourceElements = configuration.Nz;
    const double *samplingSource = source;
    const WVFieldEvaluationPlan::ResolvedRequest *firstFullRequest = nullptr;
    std::size_t firstFullInvocation = invocationCount;
    std::size_t consumer = 0;
    for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
         ++invocationIndex) {
      const auto &invocation = invocations[invocationIndex];
      for (const auto &request : invocation.plan->requests_) {
        if (request.field != field)
          continue;
        if (request.samplingKind == WVFieldSamplingKind::fullGrid &&
            firstFullRequest == nullptr) {
          std::copy(source, source + sourceElements,
                    invocation.outputs[request.outputIndex].data);
          samplingSource = invocation.outputs[request.outputIndex].data;
          firstFullRequest = &request;
          firstFullInvocation = invocationIndex;
          ++metrics_.fullGridWriteCount;
          metrics_.outputElementWriteCount += sourceElements;
        }
        ++consumer;
      }
    }
    if (consumer > 1)
      metrics_.primitiveFieldReuseCount += consumer - 1;
    for (std::size_t invocationIndex = 0; invocationIndex < invocationCount;
         ++invocationIndex) {
      const auto &invocation = invocations[invocationIndex];
      for (const auto &request : invocation.plan->requests_) {
        if (request.field != field)
          continue;
        auto &output = invocation.outputs[request.outputIndex];
        if (request.samplingKind == WVFieldSamplingKind::fullGrid) {
          if (invocationIndex != firstFullInvocation ||
              &request != firstFullRequest) {
            std::copy(samplingSource, samplingSource + sourceElements,
                      output.data);
            ++metrics_.fullGridWriteCount;
            metrics_.outputElementWriteCount += sourceElements;
          }
          continue;
        }
        if (request.samplingKind ==
            WVFieldSamplingKind::fixedVerticalProfiles) {
          for (std::size_t profile = 0;
               profile < request.profileXIndices.size(); ++profile) {
            const auto horizontalIndex =
                request.profileXIndices[profile] +
                configuration.Nx * request.profileYIndices[profile];
            for (std::size_t z = 0; z < configuration.Nz; ++z)
              output.data[z + configuration.Nz * profile] =
                  samplingSource[horizontalIndex + horizontalElements * z];
          }
          ++metrics_.profileWriteCount;
          metrics_.outputElementWriteCount += output.elementCount;
          continue;
        }
        for (std::size_t position = 0;
             position < request.positionWeights.size(); ++position) {
          const auto &weights = request.positionWeights[position];
          double value = 0.0;
          if (!weights.outsideInterpolationDomain) {
            if (request.interpolation == WVPositionInterpolation::linear) {
              const std::size_t zCount =
                  rank == WVFieldEvaluationPlan::NativeRank::volume ? 2 : 1;
              for (std::size_t iz = 0; iz < zCount; ++iz)
                for (std::size_t iy = 0; iy < 2; ++iy)
                  for (std::size_t ix = 0; ix < 2; ++ix) {
                    const auto zIndex =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zLinearIndices[iz]
                            : 0;
                    const auto sourceIndex =
                        weights.xLinearIndices[ix] +
                        configuration.Nx * weights.yLinearIndices[iy] +
                        horizontalElements * zIndex;
                    const double zWeight =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zLinearWeights[iz]
                            : 1.0;
                    value += samplingSource[sourceIndex] *
                             weights.xLinearWeights[ix] *
                             weights.yLinearWeights[iy] * zWeight;
                  }
              ++metrics_.linearInterpolationCount;
            } else {
              const std::size_t zCount =
                  rank == WVFieldEvaluationPlan::NativeRank::volume
                      ? configuration.Nz
                      : 1;
              for (std::size_t iz = 0; iz < zCount; ++iz)
                for (std::size_t iy = 0; iy < configuration.Ny; ++iy)
                  for (std::size_t ix = 0; ix < configuration.Nx; ++ix) {
                    const auto sourceIndex =
                        ix + configuration.Nx * iy + horizontalElements * iz;
                    const double zWeight =
                        rank == WVFieldEvaluationPlan::NativeRank::volume
                            ? weights.zSplineWeights[iz]
                            : 1.0;
                    value += samplingSource[sourceIndex] *
                             weights.xSplineWeights[ix] *
                             weights.ySplineWeights[iy] * zWeight;
                  }
              ++metrics_.splineInterpolationCount;
            }
          } else if (request.interpolation == WVPositionInterpolation::linear) {
            ++metrics_.linearInterpolationCount;
          } else {
            ++metrics_.splineInterpolationCount;
          }
          output.data[position] = value;
        }
        metrics_.outputElementWriteCount += output.elementCount;
      }
    }
    return WVKernelStatus::ok();
  };

  if ((dependencyMask & primitiveValues) != 0) {
    updateScratchHighWater(4 * fieldElements, 0);
    WVRealFieldBundleView primitiveBundle{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    auto status = invokeTransform(
        [&]() { return transform_->transformWaveVortexToUVWEta(state, primitiveBundle); });
    if (!status)
      return status;
    const bool primitiveNeeded[] = {
        fieldRequested(WVFieldEvaluationPlan::Field::u) ||
            fieldRequested(WVFieldEvaluationPlan::Field::ssu) ||
            fieldRequested(WVFieldEvaluationPlan::Field::uvMax),
        fieldRequested(WVFieldEvaluationPlan::Field::v) ||
            fieldRequested(WVFieldEvaluationPlan::Field::ssv) ||
            fieldRequested(WVFieldEvaluationPlan::Field::uvMax),
        fieldRequested(WVFieldEvaluationPlan::Field::w) ||
            fieldRequested(WVFieldEvaluationPlan::Field::wMax),
        fieldRequested(WVFieldEvaluationPlan::Field::eta) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoE) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal) ||
            fieldRequested(WVFieldEvaluationPlan::Field::rhoBar)};
    metrics_.primitiveFieldEvaluationCount +=
        static_cast<std::size_t>(primitiveNeeded[0]) +
        static_cast<std::size_t>(primitiveNeeded[1]) +
        static_cast<std::size_t>(primitiveNeeded[2]) +
        static_cast<std::size_t>(primitiveNeeded[3]);
    const double *u = realScratch_.data();
    const double *v = u + fieldElements;
    const double *w = v + fieldElements;
    const double *eta = w + fieldElements;
    if (fieldRequested(WVFieldEvaluationPlan::Field::u))
      writeField(WVFieldEvaluationPlan::Field::u,
                 WVFieldEvaluationPlan::NativeRank::volume, u);
    if (fieldRequested(WVFieldEvaluationPlan::Field::v))
      writeField(WVFieldEvaluationPlan::Field::v,
                 WVFieldEvaluationPlan::NativeRank::volume, v);
    if (fieldRequested(WVFieldEvaluationPlan::Field::w))
      writeField(WVFieldEvaluationPlan::Field::w,
                 WVFieldEvaluationPlan::NativeRank::volume, w);
    if (fieldRequested(WVFieldEvaluationPlan::Field::eta))
      writeField(WVFieldEvaluationPlan::Field::eta,
                 WVFieldEvaluationPlan::NativeRank::volume, eta);
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssu))
      writeField(WVFieldEvaluationPlan::Field::ssu,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 u + fieldElements - horizontalElements);
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssv))
      writeField(WVFieldEvaluationPlan::Field::ssv,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 v + fieldElements - horizontalElements);
    if (fieldRequested(WVFieldEvaluationPlan::Field::uvMax)) {
      double maximum = 0.0;
      for (std::size_t index = 0; index < fieldElements; ++index)
        maximum = std::max(maximum, std::sqrt(u[index] * u[index] +
                                              v[index] * v[index]));
      writeField(WVFieldEvaluationPlan::Field::uvMax,
                 WVFieldEvaluationPlan::NativeRank::scalar, &maximum);
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::wMax)) {
      double maximum = 0.0;
      for (std::size_t index = 0; index < fieldElements; ++index)
        maximum = std::max(maximum, std::abs(w[index]));
      writeField(WVFieldEvaluationPlan::Field::wMax,
                 WVFieldEvaluationPlan::NativeRank::scalar, &maximum);
    }
    double *derived = realScratch_.data() + 4 * fieldElements;
    const double densityScale =
        configuration.rho0 * configuration.N0 * configuration.N0 /
        configuration.g;
    if (fieldRequested(WVFieldEvaluationPlan::Field::rhoE) ||
        fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal)) {
      updateScratchHighWater(5 * fieldElements, 0);
      for (std::size_t index = 0; index < fieldElements; ++index)
        derived[index] = densityScale * eta[index];
      if (fieldRequested(WVFieldEvaluationPlan::Field::rhoE))
        writeField(WVFieldEvaluationPlan::Field::rhoE,
                   WVFieldEvaluationPlan::NativeRank::volume, derived);
      if (fieldRequested(WVFieldEvaluationPlan::Field::rhoTotal)) {
        for (std::size_t z = 0; z < configuration.Nz; ++z) {
          const double rhoNoMotion =
              configuration.rho0 - densityScale *
                                       transform_->descriptor().verticalModes().z[z];
          for (std::size_t horizontal = 0; horizontal < horizontalElements;
               ++horizontal) {
            const auto index = horizontal + horizontalElements * z;
            derived[index] = rhoNoMotion + densityScale * eta[index];
          }
        }
        writeField(WVFieldEvaluationPlan::Field::rhoTotal,
                   WVFieldEvaluationPlan::NativeRank::volume, derived);
      }
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::rhoBar)) {
      updateScratchHighWater(4 * fieldElements + configuration.Nz, 0);
      for (std::size_t z = 0; z < configuration.Nz; ++z) {
        double meanEta = 0.0;
        for (std::size_t horizontal = 0; horizontal < horizontalElements;
             ++horizontal)
          meanEta += eta[horizontal + horizontalElements * z];
        meanEta /= static_cast<double>(horizontalElements);
        derived[z] =
            configuration.rho0 -
            densityScale * transform_->descriptor().verticalModes().z[z] +
            densityScale * meanEta;
      }
      writeField(WVFieldEvaluationPlan::Field::rhoBar,
                 WVFieldEvaluationPlan::NativeRank::vertical, derived);
    }
  }

  auto fillFFieldCoefficients = [&](WVFieldEvaluationPlan::Field field) {
    const auto &modes = transform_->descriptor().verticalModes();
    const auto &horizontalModes = transform_->descriptor().fourierModes();
    auto *wave = complexScratch_.data();
    auto *zeroFrequency = wave + coefficientElements;
    const double elapsed = state.t - state.t0;
    const double f = modes.coriolisFrequency;
    for (std::size_t mode = 0; mode < horizontalModes.size(); ++mode) {
      const double Kh2 = horizontalModes[mode].Kh * horizontalModes[mode].Kh;
      for (std::size_t j = 0; j < configuration.Nj; ++j) {
        const auto index = j + configuration.Nj * mode;
        wave[index] = {};
        zeroFrequency[index] = {};
        if (field == WVFieldEvaluationPlan::Field::pi) {
          const WVComplex64 phase{std::cos(modes.omega[index] * elapsed),
                                  std::sin(modes.omega[index] * elapsed)};
          const auto positive = multiply(state.coefficients.Ap.data[index], phase);
          const auto negative =
              multiply(state.coefficients.Am.data[index], conjugate(phase));
          const double rawNAp =
              modes.NApField[index] / modes.gWaveScale[j];
          wave[index] = multiply(subtract(positive, negative), rawNAp);
        }
        if (horizontalModes[mode].Kh > 0.0) {
          const double inverseRossbySquared =
              j == 0 ? 0.0 : f * f / (configuration.g * modes.h0[j]);
          const double denominator = Kh2 + inverseRossbySquared;
          if (field == WVFieldEvaluationPlan::Field::pi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index],
                         -(f / configuration.g) / denominator);
          else if (field == WVFieldEvaluationPlan::Field::psi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], -1.0 / denominator);
          else if (field == WVFieldEvaluationPlan::Field::qgpv)
            zeroFrequency[index] = state.coefficients.A0.data[index];
        } else if (j > 0) {
          if (field == WVFieldEvaluationPlan::Field::pi)
            zeroFrequency[index] = state.coefficients.A0.data[index];
          else if (field == WVFieldEvaluationPlan::Field::psi)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], configuration.g / f);
          else if (field == WVFieldEvaluationPlan::Field::qgpv)
            zeroFrequency[index] =
                multiply(state.coefficients.A0.data[index], -f / modes.h0[j]);
        }
      }
    }
  };

  auto evaluateFField = [&](WVFieldEvaluationPlan::Field field) {
    fillFFieldCoefficients(field);
    updateScratchHighWater(4 * fieldElements, 2 * coefficientElements);
    WVRealFieldBundleView fieldAndDerivatives{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    const WVComplexConstView wave{complexScratch_.data(), spectral};
    const WVComplexConstView zeroFrequency{
        complexScratch_.data() + coefficientElements, spectral};
    return invokeTransform([&]() {
      return transform_->transformToSpatialDomainWithFAllDerivatives(
          wave, zeroFrequency, fieldAndDerivatives);
    });
  };

  if ((dependencyMask & pressureHeight) != 0) {
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::pi);
    if (!status)
      return status;
    ++metrics_.primitiveFieldEvaluationCount;
    const double *pressureHeightField = realScratch_.data();
    if (fieldRequested(WVFieldEvaluationPlan::Field::pi))
      writeField(WVFieldEvaluationPlan::Field::pi,
                 WVFieldEvaluationPlan::NativeRank::volume,
                 pressureHeightField);
    if (fieldRequested(WVFieldEvaluationPlan::Field::p)) {
      double *pressure = realScratch_.data() + 4 * fieldElements;
      updateScratchHighWater(5 * fieldElements, 2 * coefficientElements);
      const double pressureScale = configuration.rho0 * configuration.g;
      for (std::size_t index = 0; index < fieldElements; ++index)
        pressure[index] = pressureScale * pressureHeightField[index];
      writeField(WVFieldEvaluationPlan::Field::p,
                 WVFieldEvaluationPlan::NativeRank::volume, pressure);
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::ssh))
      writeField(WVFieldEvaluationPlan::Field::ssh,
                 WVFieldEvaluationPlan::NativeRank::horizontal,
                 pressureHeightField + fieldElements - horizontalElements);
  }

  if ((dependencyMask & streamfunction) != 0) {
    if (transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
      return {WVKernelStatusCode::unsupportedOperation,
              "Streamfunction evaluation is undefined when the Coriolis frequency is zero."};
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::psi);
    if (!status)
      return status;
    ++metrics_.primitiveFieldEvaluationCount;
    writeField(WVFieldEvaluationPlan::Field::psi,
               WVFieldEvaluationPlan::NativeRank::volume,
               realScratch_.data());
  }

  if ((dependencyMask & potentialVorticity) != 0) {
    auto status = evaluateFField(WVFieldEvaluationPlan::Field::qgpv);
    if (!status)
      return status;
    ++metrics_.primitiveFieldEvaluationCount;
    writeField(WVFieldEvaluationPlan::Field::qgpv,
               WVFieldEvaluationPlan::NativeRank::volume,
               realScratch_.data());
  }

  if ((dependencyMask & spectralEnergy) != 0) {
    const auto &modes = transform_->descriptor().verticalModes();
    const auto &horizontalModes = transform_->descriptor().fourierModes();
    const double f = modes.coriolisFrequency;
    double energy = 0.0;
    for (std::size_t mode = 0; mode < horizontalModes.size(); ++mode) {
      const double Kh = horizontalModes[mode].Kh;
      const double Kh2 = Kh * Kh;
      for (std::size_t j = 0; j < configuration.Nj; ++j) {
        const auto index = j + configuration.Nj * mode;
        const double M = modes.verticalWavenumber[j];
        double hWave = 1.0;
        if (j > 0)
          hWave = configuration.isHydrostatic
                      ? configuration.N0 * configuration.N0 /
                            (configuration.g * M * M)
                      : (configuration.N0 * configuration.N0 - f * f) /
                            (configuration.g * (M * M + Kh2));
        double waveFactor = 0.0;
        if (Kh > 0.0 && j > 0)
          waveFactor = 2.0 * hWave;
        else if (Kh == 0.0)
          waveFactor = j == 0 ? configuration.Lz : hWave;
        double zeroFactor = 0.0;
        if (Kh > 0.0) {
          if (j == 0)
            zeroFactor = configuration.Lz / Kh2;
          else {
            const double denominator =
                Kh2 + f * f / (configuration.g * modes.h0[j]);
            zeroFactor = modes.h0[j] / denominator;
          }
        } else if (j > 0) {
          zeroFactor = configuration.g / 2.0;
        }
        energy += waveFactor *
                      (squaredMagnitude(state.coefficients.Ap.data[index]) +
                       squaredMagnitude(state.coefficients.Am.data[index])) +
                  zeroFactor *
                      squaredMagnitude(state.coefficients.A0.data[index]);
      }
    }
    writeField(WVFieldEvaluationPlan::Field::energy,
               WVFieldEvaluationPlan::NativeRank::scalar, &energy);
  }

  const auto derivativeDependencies =
      uDerivatives | vDerivatives | wDerivatives;
  if ((dependencyMask & derivativeDependencies) != 0) {
    updateScratchHighWater(6 * fieldElements, 0);
    double *zetaX = realScratch_.data();
    double *zetaY = zetaX + fieldElements;
    double *zetaZ = zetaY + fieldElements;
    double *derivatives = zetaZ + fieldElements;
    std::fill(zetaX, zetaX + 3 * fieldElements, 0.0);
    WVRealFieldBundleView derivativeBundle{
        derivatives,
        {configuration.Nx, configuration.Ny, configuration.Nz, 3}};
    if ((dependencyMask & uDerivatives) != 0) {
      auto status = invokeTransform([&]() {
        return transform_->transformStateFieldDerivatives(
            state, WVDynamicalField::u, derivativeBundle);
      });
      if (!status)
        return status;
      ++metrics_.primitiveFieldEvaluationCount;
      const double *uy = derivatives + fieldElements;
      const double *uz = uy + fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        zetaY[index] += uz[index];
        zetaZ[index] -= uy[index];
      }
    }
    if ((dependencyMask & vDerivatives) != 0) {
      auto status = invokeTransform([&]() {
        return transform_->transformStateFieldDerivatives(
            state, WVDynamicalField::v, derivativeBundle);
      });
      if (!status)
        return status;
      ++metrics_.primitiveFieldEvaluationCount;
      const double *vx = derivatives;
      const double *vz = derivatives + 2 * fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        zetaX[index] -= vz[index];
        zetaZ[index] += vx[index];
      }
    }
    if ((dependencyMask & wDerivatives) != 0) {
      auto status = invokeTransform([&]() {
        return transform_->transformStateFieldDerivatives(
            state, WVDynamicalField::w, derivativeBundle);
      });
      if (!status)
        return status;
      ++metrics_.primitiveFieldEvaluationCount;
      const double *wx = derivatives;
      const double *wy = derivatives + fieldElements;
      for (std::size_t index = 0; index < fieldElements; ++index) {
        zetaX[index] += wy[index];
        zetaY[index] -= wx[index];
      }
    }
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaX))
      writeField(WVFieldEvaluationPlan::Field::zetaX,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaX);
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaY))
      writeField(WVFieldEvaluationPlan::Field::zetaY,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaY);
    if (fieldRequested(WVFieldEvaluationPlan::Field::zetaZ))
      writeField(WVFieldEvaluationPlan::Field::zetaZ,
                 WVFieldEvaluationPlan::NativeRank::volume, zetaZ);
  }

  return WVKernelStatus::ok();
}

WVKernelStatus WVFieldEvaluationService::createEventPlan(
    const std::vector<WVEventFieldRequest> &requests,
    WVEventFieldEvaluationPlan &plan) {
  if (barotropicQG_)
    return barotropicQG_->createEventPlan(requests, plan);
  try {
    WVEventFieldEvaluationPlan candidate;
    const auto &configuration = transform_->descriptor().configuration();
    candidate.configuration_ = configuration;
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    std::size_t positionSetCount = 0;
    for (std::size_t outputIndex = 0; outputIndex < requests.size();
         ++outputIndex) {
      const auto &request = requests[outputIndex];
      if (request.identifier.empty() ||
          !identifiers.insert(request.identifier).second)
        return invalid(
            "Event-field request identifiers must be nonempty and unique.");
      if (request.positionSetSlot ==
          std::numeric_limits<std::size_t>::max())
        return {WVKernelStatusCode::sizeOverflow,
                "An event-field position-set slot overflows its plan."};
      if (request.interpolation != WVPositionInterpolation::linear &&
          request.interpolation != WVPositionInterpolation::spline)
        return invalid("Event-field interpolation method is invalid.");
      const auto *metadata = findPortableVariable(request.fieldName);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field ||
          (metadata->samplingMask & portablePositionSampling) == 0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Event-position sampling does not support field " +
                    request.fieldName + "."};
      if (metadata->naturalRank != WVPortableNaturalRank::volume &&
          metadata->naturalRank != WVPortableNaturalRank::horizontal)
        return {WVKernelStatusCode::unsupportedOperation,
                "Event-position sampling requires a volume or horizontal "
                "field: " +
                    request.fieldName + "."};
      if (metadata->identifier == WVPortableVariable::psi &&
          transform_->descriptor().verticalModes().coriolisFrequency == 0.0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Streamfunction evaluation is undefined when the Coriolis "
                "frequency is zero."};

      candidate.requests_.push_back(
          {metadata->identifier, metadata->naturalRank,
           metadata->primitiveDependencyMask, request.positionSetSlot,
           request.interpolation, outputIndex});
      candidate.outputs_.push_back(
          {request.identifier, request.fieldName, metadata->identifier,
           metadata->naturalRank, metadata->primitiveDependencyMask,
           request.positionSetSlot, request.interpolation});
      candidate.requestedFieldMask_ |=
          1ULL << static_cast<std::size_t>(metadata->ordinal);
      candidate.dependencyMask_ |= metadata->primitiveDependencyMask;
      positionSetCount =
          std::max(positionSetCount, request.positionSetSlot + 1);
    }
    candidate.positionSetCount_ = positionSetCount;
    candidate.requiresZByPositionSet_.assign(positionSetCount, 0);
    for (const auto &request : candidate.requests_)
      if (request.nativeRank == WVPortableNaturalRank::volume)
        candidate.requiresZByPositionSet_[request.positionSetSlot] = 1;

    std::uint64_t fingerprint = fingerprintOffset;
    appendConfigurationFingerprint(fingerprint, configuration);
    appendFingerprint(fingerprint, candidate.positionSetCount_);
    appendFingerprint(fingerprint, candidate.requestedFieldMask_);
    appendFingerprint(fingerprint, candidate.dependencyMask_);
    const auto requestCount = candidate.requests_.size();
    appendFingerprint(fingerprint, requestCount);
    for (std::size_t index = 0; index < candidate.requests_.size(); ++index) {
      const auto &resolved = candidate.requests_[index];
      const auto &output = candidate.outputs_[index];
      appendFingerprint(fingerprint, resolved.field);
      appendFingerprint(fingerprint, resolved.nativeRank);
      appendFingerprint(fingerprint, resolved.dependencyMask);
      appendFingerprint(fingerprint, resolved.positionSetSlot);
      appendFingerprint(fingerprint, resolved.interpolation);
      appendFingerprint(fingerprint, output.identifier.data(),
                        output.identifier.size());
      appendFingerprint(fingerprint, output.fieldName.data(),
                        output.fieldName.size());
    }
    candidate.fingerprint_ = fingerprint;
    const auto planBytes = candidate.persistentBytes();
    plan = std::move(candidate);
    ++metrics_.eventPlanCreationCount;
    metrics_.eventPlanFieldResolutionCount += requests.size();
    metrics_.lastEventPlanBytes = planBytes;
    metrics_.maximumEventPlanBytes =
        std::max(metrics_.maximumEventPlanBytes, planBytes);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate an event-field evaluation plan."};
  }
}

WVKernelStatus WVFieldEvaluationService::prepareEventGeometry(
    const WVEventFieldEvaluationPlan &plan,
    const WVEventPositionSetView *positionSets,
    std::size_t positionSetCount, WVPreparedFieldGeometry &geometry) {
  if (barotropicQG_)
    return barotropicQG_->prepareEventGeometry(
        plan, positionSets, positionSetCount, geometry);
  if (!sameTransformConfiguration(
          plan.configuration_, transform_->descriptor().configuration()))
    return invalid(
        "The event-field plan belongs to a different transform configuration.");
  if (positionSetCount != plan.positionSetCount_)
    return {WVKernelStatusCode::invalidShape,
            "Event position-set count must match the resolved field plan."};
  if (positionSetCount != 0 && positionSets == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Event position sets have a null view pointer."};

  try {
    WVPreparedFieldGeometry candidate;
    candidate.fieldPlanFingerprint_ = plan.fingerprint_;
    candidate.positionSets_.reserve(positionSetCount);
    std::uint64_t geometryFingerprint = fingerprintOffset;
    appendFingerprint(geometryFingerprint, plan.fingerprint_);
    appendFingerprint(geometryFingerprint, positionSetCount);
    for (std::size_t slot = 0; slot < positionSetCount; ++slot) {
      const auto &view = positionSets[slot];
      if (view.extentCount != 0 && view.extents == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Event position-set extents have a null pointer."};
      if (view.positionCount != 0 &&
          (view.x == nullptr || view.y == nullptr ||
           (plan.requiresZByPositionSet_[slot] != 0 && view.z == nullptr)))
        return {WVKernelStatusCode::invalidPointer,
                "Event coordinates required by a resolved field are missing."};

      WVPreparedFieldGeometry::PositionSet prepared;
      prepared.x = view.x;
      prepared.y = view.y;
      prepared.z = view.z;
      prepared.positionCount = view.positionCount;
      if (view.extentCount == 0)
        prepared.extents.push_back(view.positionCount);
      else
        prepared.extents.assign(view.extents,
                                view.extents + view.extentCount);
      std::size_t extentProduct = 1;
      for (const auto extent : prepared.extents)
        extentProduct = checkedProduct(extentProduct, extent);
      if (extentProduct != view.positionCount)
        return {WVKernelStatusCode::invalidShape,
                "Event position-set extents do not match its sample count."};
      for (std::size_t position = 0; position < view.positionCount;
           ++position) {
        if (!std::isfinite(view.x[position]) ||
            !std::isfinite(view.y[position]) ||
            (view.z != nullptr && !std::isfinite(view.z[position])))
          return invalid("Event position coordinates must be finite.");
      }
      if (candidate.positionCount_ >
          std::numeric_limits<std::size_t>::max() - view.positionCount)
        return {WVKernelStatusCode::sizeOverflow,
                "Event position count overflows its prepared geometry."};
      candidate.positionCount_ += view.positionCount;
      const std::size_t coordinateArrayCount = view.z == nullptr ? 2 : 3;
      const auto coordinateBytes = checkedProduct(
          checkedProduct(view.positionCount, coordinateArrayCount),
          sizeof(double));
      if (candidate.borrowedCoordinateBytes_ >
          std::numeric_limits<std::size_t>::max() - coordinateBytes)
        return {WVKernelStatusCode::sizeOverflow,
                "Event coordinate storage overflows its metrics."};
      candidate.borrowedCoordinateBytes_ += coordinateBytes;

      appendFingerprint(geometryFingerprint, slot);
      appendFingerprint(geometryFingerprint, view.positionCount);
      const auto extentCount = prepared.extents.size();
      appendFingerprint(geometryFingerprint, extentCount);
      for (const auto extent : prepared.extents)
        appendFingerprint(geometryFingerprint, extent);
      const bool hasZ = view.z != nullptr;
      appendFingerprint(geometryFingerprint, hasZ);
      const auto coordinateByteCount = view.positionCount * sizeof(double);
      appendFingerprint(geometryFingerprint, view.x, coordinateByteCount);
      appendFingerprint(geometryFingerprint, view.y, coordinateByteCount);
      if (hasZ)
        appendFingerprint(geometryFingerprint, view.z, coordinateByteCount);
      candidate.positionSets_.push_back(std::move(prepared));
    }
    candidate.geometryFingerprint_ = geometryFingerprint;

    const auto &configuration = transform_->descriptor().configuration();
    candidate.evaluationPlan_.configuration_ = configuration;
    candidate.evaluationPlan_.requestedFieldMask_ = plan.requestedFieldMask_;
    candidate.evaluationPlan_.dependencyMask_ = plan.dependencyMask_;
    candidate.evaluationPlan_.requests_.reserve(plan.requests_.size());
    candidate.evaluationPlan_.outputs_.reserve(plan.outputs_.size());
    candidate.outputs_.reserve(plan.outputs_.size());
    const double dx =
        configuration.Lx / static_cast<double>(configuration.Nx);
    const double dy =
        configuration.Ly / static_cast<double>(configuration.Ny);
    const double dz =
        configuration.Lz / static_cast<double>(configuration.Nz - 1);

    for (const auto &eventRequest : plan.requests_) {
      const auto &set =
          candidate.positionSets_[eventRequest.positionSetSlot];
      WVFieldEvaluationPlan::ResolvedRequest request;
      request.field = eventRequest.field;
      request.nativeRank = eventRequest.nativeRank;
      request.samplingKind = WVFieldSamplingKind::positions;
      request.interpolation = eventRequest.interpolation;
      request.outputIndex = eventRequest.outputIndex;
      request.positionWeights.reserve(set.positionCount);
      for (std::size_t position = 0; position < set.positionCount;
           ++position) {
        const double x = set.x[position];
        const double y = set.y[position];
        const double z = eventRequest.nativeRank ==
                                 WVPortableNaturalRank::volume
                             ? set.z[position]
                             : 0.0;
        const double xWrapped = wrapped(x, configuration.Lx);
        const double yWrapped = wrapped(y, configuration.Ly);
        const auto xLower = std::min(
            static_cast<std::size_t>(std::floor(xWrapped / dx)),
            configuration.Nx - 1);
        const auto yLower = std::min(
            static_cast<std::size_t>(std::floor(yWrapped / dy)),
            configuration.Ny - 1);
        WVFieldEvaluationPlan::PositionWeights weights;
        if (eventRequest.interpolation == WVPositionInterpolation::linear) {
          weights.xLinearIndices = {xLower,
                                    (xLower + 1) % configuration.Nx};
          weights.yLinearIndices = {yLower,
                                    (yLower + 1) % configuration.Ny};
          const double xFraction =
              (xWrapped - static_cast<double>(xLower) * dx) / dx;
          const double yFraction =
              (yWrapped - static_cast<double>(yLower) * dy) / dy;
          weights.xLinearWeights = {1.0 - xFraction, xFraction};
          weights.yLinearWeights = {1.0 - yFraction, yFraction};
          if (eventRequest.nativeRank == WVPortableNaturalRank::volume) {
            weights.outsideInterpolationDomain =
                z < -configuration.Lz || z > 0.0;
            if (!weights.outsideInterpolationDomain) {
              const double normalizedZ = (z + configuration.Lz) / dz;
              const auto zLower = std::min(
                  static_cast<std::size_t>(
                      std::max(0.0, std::floor(normalizedZ))),
                  configuration.Nz - 2);
              const double zFraction = std::clamp(
                  normalizedZ - static_cast<double>(zLower), 0.0, 1.0);
              weights.zLinearIndices = {zLower, zLower + 1};
              weights.zLinearWeights = {1.0 - zFraction, zFraction};
            }
          }
        } else {
          const bool xBoundary =
              xLower < 3 || xLower > configuration.Nx - 4;
          const bool yBoundary =
              yLower < 3 || yLower > configuration.Ny - 4;
          const std::size_t xShift = xBoundary ? 4 : 0;
          const std::size_t yShift = yBoundary ? 4 : 0;
          const double xQuery =
              xBoundary ? wrapped(x + 4.0 * dx, configuration.Lx) : xWrapped;
          const double yQuery =
              yBoundary ? wrapped(y + 4.0 * dy, configuration.Ly) : yWrapped;
          weights.outsideInterpolationDomain =
              xQuery > static_cast<double>(configuration.Nx - 1) * dx ||
              yQuery > static_cast<double>(configuration.Ny - 1) * dy;
          if (!weights.outsideInterpolationDomain) {
            movingWorkspace_->xSpline.weightsInto(
                0.0, dx, xQuery, weights.xSplineWeights, xShift,
                &movingWorkspace_->xRightHandSide,
                &movingWorkspace_->xShifted);
            movingWorkspace_->ySpline.weightsInto(
                0.0, dy, yQuery, weights.ySplineWeights, yShift,
                &movingWorkspace_->yRightHandSide,
                &movingWorkspace_->yShifted);
          }
          if (eventRequest.nativeRank == WVPortableNaturalRank::volume) {
            weights.outsideInterpolationDomain =
                weights.outsideInterpolationDomain ||
                z < -configuration.Lz || z > 0.0;
            if (!weights.outsideInterpolationDomain)
              movingWorkspace_->zSpline.weightsInto(
                  -configuration.Lz, dz, z, weights.zSplineWeights, 0,
                  &movingWorkspace_->zRightHandSide,
                  &movingWorkspace_->zShifted);
          }
        }
        request.positionWeights.push_back(std::move(weights));
      }
      candidate.evaluationPlan_.requests_.push_back(std::move(request));
      candidate.evaluationPlan_.outputs_.push_back(
          {std::string{}, std::string{}, WVFieldSamplingKind::positions, {},
           set.positionCount});
      candidate.outputs_.push_back(
          {eventRequest.outputIndex, eventRequest.positionSetSlot,
           set.extents, set.positionCount});
    }

    const auto retainedBytes = candidate.retainedBytes();
    const auto liveBytes = candidate.liveBytes();
    geometry = std::move(candidate);
    ++metrics_.eventGeometryPreparationCount;
    metrics_.eventPositionSetCount += positionSetCount;
    metrics_.eventPositionCount += geometry.positionCount_;
    metrics_.lastPreparedGeometryRetainedBytes = retainedBytes;
    metrics_.maximumPreparedGeometryRetainedBytes =
        std::max(metrics_.maximumPreparedGeometryRetainedBytes,
                 retainedBytes);
    metrics_.lastPreparedGeometryLiveBytes = liveBytes;
    metrics_.maximumPreparedGeometryLiveBytes =
        std::max(metrics_.maximumPreparedGeometryLiveBytes, liveBytes);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &) {
    return {WVKernelStatusCode::sizeOverflow,
            "Event geometry extents or storage overflow size_t."};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate prepared event geometry."};
  }
}

WVKernelStatus WVFieldEvaluationService::evaluateEvent(
    const WVEventFieldEvaluationPlan &plan,
    const WVPreparedFieldGeometry &geometry, const WVState &state,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  const WVEventFieldEvaluationBatchEntry entry{&plan, &geometry, outputs,
                                               outputCount};
  return evaluateEventBatch(state, &entry, 1);
}

WVKernelStatus WVFieldEvaluationService::evaluateEvent(
    const WVEventFieldEvaluationPlan &plan,
    const WVPreparedFieldGeometry &geometry,
    const WVIntegrationState &state, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  const WVEventFieldEvaluationBatchEntry entry{&plan, &geometry, outputs,
                                               outputCount};
  return evaluateEventBatch(state, &entry, 1);
}

WVKernelStatus WVFieldEvaluationService::evaluateEventBatch(
    const WVState &state, const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entryCount != 0 && entries == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Event field batch has a null entry pointer."};
  std::size_t batchOutputCount = 0;
  for (std::size_t entryIndex = 0; entryIndex < entryCount; ++entryIndex) {
    const auto &entry = entries[entryIndex];
    if (entry.plan == nullptr || entry.geometry == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Event field batch entry has a null plan or geometry."};
    const auto &plan = *entry.plan;
    const auto &geometry = *entry.geometry;
    if (!sameTransformConfiguration(plan.configuration_,
                                    transform_->descriptor().configuration()))
      return invalid("The event-field plan belongs to a different transform "
                     "configuration.");
    if (geometry.fieldPlanFingerprint_ != plan.fingerprint_ ||
        geometry.outputCount() != plan.outputCount())
      return invalid(
          "Prepared event geometry does not match the resolved field plan.");
    if (entry.outputCount != geometry.outputs_.size())
      return {WVKernelStatusCode::invalidShape,
              "Event output-view count must match the prepared geometry."};
    if (entry.outputCount != 0 && entry.outputs == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Event outputs have a null view pointer."};
    for (std::size_t outputIndex = 0; outputIndex < entry.outputCount;
         ++outputIndex)
      if (entry.outputs[outputIndex].elementCount !=
          geometry.outputs_[outputIndex].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Event output has the wrong element count for request " +
                    plan.outputs_[outputIndex].identifier + "."};
    batchOutputCount =
        batchOutputCount >
                std::numeric_limits<std::size_t>::max() - entry.outputCount
            ? std::numeric_limits<std::size_t>::max()
            : batchOutputCount + entry.outputCount;
  }

  try {
    eventBatchInvocations_.clear();
    eventBatchInvocations_.reserve(entryCount);
    for (std::size_t entryIndex = 0; entryIndex < entryCount; ++entryIndex) {
      const auto &entry = entries[entryIndex];
      eventBatchInvocations_.push_back(
          {&entry.geometry->evaluationPlan_, entry.outputs, entry.outputCount});
    }
  } catch (const std::bad_alloc &) {
    eventBatchInvocations_.clear();
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate event field batch invocation storage."};
  }
  metrics_.eventBatchInvocationWorkspaceBytes =
      eventBatchInvocations_.capacity() * sizeof(PlanInvocation);
  metrics_.servicePersistentBytes = persistentBytes();
  const auto status = evaluatePlanBatch(eventBatchInvocations_.data(),
                                        eventBatchInvocations_.size(), state);
  eventBatchInvocations_.clear();
  if (status) {
    metrics_.eventEvaluationCount += entryCount;
    ++metrics_.eventBatchEvaluationCount;
    metrics_.eventBatchOccurrenceCount += entryCount;
    metrics_.eventBatchOutputCount += batchOutputCount;
  }
  return status;
}

WVKernelStatus WVFieldEvaluationService::evaluateEventBatch(
    const WVIntegrationState &state,
    const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (barotropicQG_)
    return barotropicQG_->evaluateEventBatch(state, entries, entryCount);
  return evaluateEventBatch(state.waveVortex, entries, entryCount);
}

WVKernelStatus WVFieldEvaluationService::createMovingPlan(
    const std::vector<WVMovingFieldRequest> &requests,
    WVMovingFieldEvaluationPlan &plan) const {
  if (barotropicQG_)
    return barotropicQG_->createMovingPlan(requests, plan);
  try {
    WVMovingFieldEvaluationPlan candidate;
    candidate.configuration_ = transform_->descriptor().configuration();
    candidate.requests_.reserve(requests.size());
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &request = requests[index];
      if (request.identifier.empty() ||
          !identifiers.insert(request.identifier).second)
        return invalid(
            "Moving-field request identifiers must be nonempty and unique.");
      if (request.positionCount == 0 ||
          request.positionOffset >
              std::numeric_limits<std::size_t>::max() - request.positionCount)
        return invalid("Moving-field request position range is invalid.");
      if (request.interpolation != WVPositionInterpolation::linear &&
          request.interpolation != WVPositionInterpolation::spline)
        return invalid("Moving-field interpolation method is invalid.");
      const auto *metadata = findPortableVariable(request.fieldName);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field ||
          metadata->movingPrimitiveChannel < 0 ||
          (metadata->samplingMask & portablePositionSampling) == 0)
        return {WVKernelStatusCode::unsupportedOperation,
                "Moving-position sampling does not support field " +
                    request.fieldName + "."};
      const auto channel =
          static_cast<std::size_t>(metadata->movingPrimitiveChannel);
      candidate.positionCount_ =
          std::max(candidate.positionCount_,
                   request.positionOffset + request.positionCount);
      candidate.requests_.push_back({channel, request.positionOffset,
                                     request.positionCount,
                                     request.interpolation, index});
      candidate.outputs_.push_back(
          {request.identifier, request.fieldName,
           WVFieldSamplingKind::positions, {request.positionCount},
           request.positionCount});
    }
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a moving-field evaluation plan."};
  }
}

WVKernelStatus WVFieldEvaluationService::evaluateMoving(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  return evaluateMovingImpl(plan,state,nullptr,positions,outputs,outputCount);
}

WVKernelStatus WVFieldEvaluationService::evaluateMoving(
    const WVMovingFieldEvaluationPlan &plan,
    const WVIntegrationState &state, WVMovingPositionView positions,
    WVFieldOutputView *outputs, std::size_t outputCount) {
  if (barotropicQG_)
    return barotropicQG_->evaluateMoving(plan, state, positions, outputs,
                                         outputCount);
  return evaluateMoving(plan, state.waveVortex, positions, outputs,
                        outputCount);
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  return evaluateMovingImpl(plan,state,&advectionFields,positions,outputs,outputCount);
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &plan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  if (barotropicQG_)
    return barotropicQG_->evaluateMovingFromAdvectionFields(
        plan, state, advectionFields, positions, outputs, outputCount);
  return evaluateMovingFromAdvectionFields(
      plan, state.waveVortex, advectionFields, positions, outputs,
      outputCount);
}

WVRealFieldBundleView WVFieldEvaluationService::advectionFieldStorage() noexcept {
  if (barotropicQG_)
    return {};
  const auto spatial = transform_->descriptor().spatialShape();
  return {realScratch_.data(),{spatial.first,spatial.second,spatial.third,3}};
}

WVKernelStatus WVFieldEvaluationService::evaluateMovingImpl(
    const WVMovingFieldEvaluationPlan &plan, const WVState &state,
    const WVRealFieldBundleConstView *preparedAdvectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount) {
  if (!sameTransformConfiguration(
          plan.configuration_, transform_->descriptor().configuration()))
    return invalid(
        "The moving-field plan belongs to a different transform configuration.");
  if (positions.positionCount != plan.positionCount_ ||
      (positions.positionCount != 0 &&
       (positions.x == nullptr || positions.y == nullptr ||
        positions.z == nullptr)))
    return {WVKernelStatusCode::invalidShape,
            "Moving coordinates must match the plan's shared position count."};
  if (outputCount != plan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Moving-field outputs must match the plan."};
  const auto spectral = transform_->descriptor().spectralShape();
  for (const auto view : {state.coefficients.Ap, state.coefficients.Am,
                          state.coefficients.A0})
    if (view.data == nullptr || view.shape.rows != spectral.rows ||
        view.shape.columns != spectral.columns)
      return {WVKernelStatusCode::invalidShape,
              "Moving-field coefficients must have shape [Nj,Nkl]."};
  for (std::size_t index = 0; index < outputCount; ++index)
    if (outputs[index].data == nullptr ||
        outputs[index].elementCount != plan.outputs_[index].elementCount)
      return {WVKernelStatusCode::invalidShape,
              "Moving-field output shape does not match its request."};
  for (std::size_t index = 0; index < positions.positionCount; ++index)
    if (!std::isfinite(positions.x[index]) ||
        !std::isfinite(positions.y[index]) ||
        !std::isfinite(positions.z[index]))
      return invalid("Moving coordinates must be finite.");
  ExecutionGuard guard(executing_);
  if (!guard.entered())
    return {WVKernelStatusCode::reentrantExecution,
            "Field evaluation is not reentrant."};

  const auto &configuration = transform_->descriptor().configuration();
  const auto spatial = transform_->descriptor().spatialShape();
  const auto R = spatial.elementCount();
  const auto horizontalCount = configuration.Nx * configuration.Ny;
  const double *primitiveFields = realScratch_.data();
  if (preparedAdvectionFields == nullptr) {
    WVRealFieldBundleView fields{
        realScratch_.data(),
        {configuration.Nx, configuration.Ny, configuration.Nz, 4}};
    const auto before = transform_->metrics().executionCount;
    const auto status = transform_->transformWaveVortexToUVWEta(state, fields);
    if (!status)
      return status;
    metrics_.fftExecutionCount += transform_->metrics().executionCount - before;
    ++metrics_.movingPrimitiveTransformCount;
  } else {
    if (preparedAdvectionFields->data == nullptr ||
        preparedAdvectionFields->shape.first != configuration.Nx ||
        preparedAdvectionFields->shape.second != configuration.Ny ||
        preparedAdvectionFields->shape.third != configuration.Nz ||
        preparedAdvectionFields->shape.fourth != 3)
      return {WVKernelStatusCode::invalidShape,
              "Prepared advection fields must have shape [Nx,Ny,Nz,3]."};
    if (std::any_of(plan.requests_.begin(),plan.requests_.end(),[](const auto &request) {
          return request.primitiveChannel > 2;
        }))
      return {WVKernelStatusCode::unsupportedOperation,
              "Prepared advection fields support only u, v, and w requests."};
    primitiveFields = preparedAdvectionFields->data;
    ++metrics_.primitiveFieldReuseCount;
  }
  ++metrics_.evaluationCount;
  ++metrics_.movingEvaluationCount;
  metrics_.movingPositionCount += positions.positionCount;
  ++metrics_.transformCount;
  metrics_.primitiveFieldEvaluationCount += preparedAdvectionFields == nullptr ? 4 : 3;
  metrics_.scratchHighWaterBytes =
      std::max(metrics_.scratchHighWaterBytes, 4 * R * sizeof(double));

  const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
  const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
  const double dz = configuration.Lz / static_cast<double>(configuration.Nz - 1);
  const double densityScale =
      configuration.rho0 * configuration.N0 * configuration.N0 /
      configuration.g;
  const bool needsTotalDensity = std::any_of(
      plan.requests_.begin(), plan.requests_.end(), [](const auto &request) {
        return request.primitiveChannel == 5;
      });
  double *totalDensity = nullptr;
  if (needsTotalDensity) {
    totalDensity = realScratch_.data() + 4 * R;
    const auto &vertical = transform_->descriptor().verticalModes().z;
    const double *eta = primitiveFields + 3 * R;
    for (std::size_t iz = 0; iz < configuration.Nz; ++iz) {
      const double densityNoMotion =
          configuration.rho0 - densityScale * vertical[iz];
      for (std::size_t horizontal = 0; horizontal < horizontalCount;
           ++horizontal) {
        const auto index = horizontal + horizontalCount * iz;
        totalDensity[index] = densityNoMotion + densityScale * eta[index];
      }
    }
    metrics_.scratchHighWaterBytes =
        std::max(metrics_.scratchHighWaterBytes, 5 * R * sizeof(double));
  }
  for (const auto &request : plan.requests_) {
    auto &output = outputs[request.outputIndex];
    for (std::size_t local = 0; local < request.positionCount; ++local) {
      const auto position = request.positionOffset + local;
      const double x = wrapped(positions.x[position], configuration.Lx);
      const double y = wrapped(positions.y[position], configuration.Ly);
      const double z = positions.z[position];
      double value = 0.0;
      if (z >= -configuration.Lz && z <= 0.0) {
        const auto channel = std::min<std::size_t>(request.primitiveChannel, 3);
        const double *source = request.primitiveChannel == 5
                                   ? totalDensity
                                   : primitiveFields + channel * R;
        if (request.interpolation == WVPositionInterpolation::linear) {
          const auto x0 = std::min(
              static_cast<std::size_t>(std::floor(x / dx)),
              configuration.Nx - 1);
          const auto y0 = std::min(
              static_cast<std::size_t>(std::floor(y / dy)),
              configuration.Ny - 1);
          const double normalizedZ = (z + configuration.Lz) / dz;
          const auto z0 = std::min(
              static_cast<std::size_t>(
                  std::max(0.0, std::floor(normalizedZ))),
              configuration.Nz - 2);
          const double wx1 = (x - static_cast<double>(x0) * dx) / dx;
          const double wy1 = (y - static_cast<double>(y0) * dy) / dy;
          const double wz1 = std::clamp(normalizedZ - static_cast<double>(z0),
                                        0.0, 1.0);
          const std::array<std::size_t, 2> xi{{x0,
                                               (x0 + 1) % configuration.Nx}};
          const std::array<std::size_t, 2> yi{{y0,
                                               (y0 + 1) % configuration.Ny}};
          const std::array<double, 2> wx{{1.0 - wx1, wx1}};
          const std::array<double, 2> wy{{1.0 - wy1, wy1}};
          const std::array<double, 2> wz{{1.0 - wz1, wz1}};
          for (std::size_t iz = 0; iz < 2; ++iz)
            for (std::size_t iy = 0; iy < 2; ++iy)
              for (std::size_t ix = 0; ix < 2; ++ix)
                value += source[xi[ix] + configuration.Nx * yi[iy] +
                                horizontalCount * (z0 + iz)] *
                         wx[ix] * wy[iy] * wz[iz];
          ++metrics_.linearInterpolationCount;
        } else {
          const auto xLower = std::min(
              static_cast<std::size_t>(std::floor(x / dx)),
              configuration.Nx - 1);
          const auto yLower = std::min(
              static_cast<std::size_t>(std::floor(y / dy)),
              configuration.Ny - 1);
          const bool xBoundary =
              xLower < 3 || xLower > configuration.Nx - 4;
          const bool yBoundary =
              yLower < 3 || yLower > configuration.Ny - 4;
          const std::size_t xShift = xBoundary ? 4 : 0;
          const std::size_t yShift = yBoundary ? 4 : 0;
          const double xQuery =
              xBoundary ? wrapped(x + 4.0 * dx, configuration.Lx) : x;
          const double yQuery =
              yBoundary ? wrapped(y + 4.0 * dy, configuration.Ly) : y;
          if (xQuery <= static_cast<double>(configuration.Nx - 1) * dx &&
              yQuery <= static_cast<double>(configuration.Ny - 1) * dy) {
            movingWorkspace_->xSpline.weightsInto(
                0.0, dx, xQuery, movingWorkspace_->xWeights, xShift,
                &movingWorkspace_->xRightHandSide,
                &movingWorkspace_->xShifted);
            movingWorkspace_->ySpline.weightsInto(
                0.0, dy, yQuery, movingWorkspace_->yWeights, yShift,
                &movingWorkspace_->yRightHandSide,
                &movingWorkspace_->yShifted);
            movingWorkspace_->zSpline.weightsInto(
                -configuration.Lz, dz, z, movingWorkspace_->zWeights, 0,
                &movingWorkspace_->zRightHandSide,
                &movingWorkspace_->zShifted);
            for (std::size_t iz = 0; iz < configuration.Nz; ++iz)
              for (std::size_t iy = 0; iy < configuration.Ny; ++iy)
                for (std::size_t ix = 0; ix < configuration.Nx; ++ix)
                  value += source[ix + configuration.Nx * iy +
                                  horizontalCount * iz] *
                           movingWorkspace_->xWeights[ix] *
                           movingWorkspace_->yWeights[iy] *
                           movingWorkspace_->zWeights[iz];
          }
          ++metrics_.splineInterpolationCount;
        }
        if (request.primitiveChannel == 4)
          value *= densityScale;
      } else if (request.interpolation == WVPositionInterpolation::linear) {
        ++metrics_.linearInterpolationCount;
      } else {
        ++metrics_.splineInterpolationCount;
      }
      output.data[local] = value;
    }
    metrics_.outputElementWriteCount += output.elementCount;
  }
  return WVKernelStatus::ok();
}

const WVTransformConstantStratificationConfiguration &
WVFieldEvaluationService::configuration() const noexcept {
  return transform_->descriptor().configuration();
}

WVKernelStatus WVFieldEvaluationService::createStateLayout(
    const WVPortableObserverDescriptor &descriptor,
    WVIntegrationStateLayout &layout) const {
  (void)descriptor;
  if (barotropicQG_) {
    const auto &configuration = barotropicQG_->configuration();
    WVTransformBarotropicQGDescriptor transform;
    auto status = WVTransformBarotropicQGDescriptor::create(configuration,
                                                             transform);
    if (!status)
      return status;
    return WVIntegrationStateLayout::createCoefficientOnly(
        {"WVTransformBarotropicQG",
         {configuration.Nx, configuration.Ny},
         {{"A0", {transform.Nkl()},
           WVToleranceKind::coefficientEnergyScaled}},
         true},
        layout);
  }
  const auto &configuration = transform_->descriptor().configuration();
  const auto shape = transform_->descriptor().spectralShape();
  return WVIntegrationStateLayout::createCoefficientOnly(
      {"WVTransformConstantStratification",
       {configuration.Nx, configuration.Ny, configuration.Nz},
       {{"Ap", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled},
        {"Am", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled},
        {"A0", {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled}},
       true},
      layout);
}

bool WVFieldEvaluationService::isCompatibleWith(
    const WVIntegrationStateLayout &layout) const noexcept {
  if (barotropicQG_)
    return barotropicQG_->isCompatibleWith(layout);
  if (!transform_ || !layout.hasLegacyCoefficientTriple())
    return false;
  const auto shape = transform_->descriptor().spectralShape();
  const auto &configuration = transform_->descriptor().configuration();
  return layout.coefficientShape().rows == shape.rows &&
         layout.coefficientShape().columns == shape.columns &&
         layout.spatialDimensions() ==
             std::vector<std::size_t>{configuration.Nx, configuration.Ny,
                                      configuration.Nz};
}

bool WVFieldEvaluationService::isCompatibleWith(
    const WVFieldEvaluationService &other) const noexcept {
  if (static_cast<bool>(barotropicQG_) !=
      static_cast<bool>(other.barotropicQG_))
    return false;
  if (barotropicQG_)
    return sameTransformConfiguration(barotropicQG_->configuration(),
                                      other.barotropicQG_->configuration());
  return transform_ != nullptr && other.transform_ != nullptr &&
         sameTransformConfiguration(transform_->descriptor().configuration(),
                                    other.transform_->descriptor()
                                        .configuration());
}

const WVFieldEvaluationMetrics &
WVFieldEvaluationService::metrics() const noexcept {
  return barotropicQG_ ? barotropicQG_->metrics() : metrics_;
}

std::size_t WVFieldEvaluationService::persistentBytes() const noexcept {
  if (barotropicQG_)
    return sizeof(*this) + barotropicQG_->persistentBytes();
  return sizeof(*this) +
         (ownedTransform_ ? transform_->persistentBytes() : 0) +
         realScratch_.capacity() * sizeof(double) +
         complexScratch_.capacity() * sizeof(WVComplex64) +
         eventBatchInvocations_.capacity() * sizeof(PlanInvocation) +
         (movingWorkspace_ ? movingWorkspace_->persistentBytes() : 0);
}

} // namespace wavevortex::runtime
