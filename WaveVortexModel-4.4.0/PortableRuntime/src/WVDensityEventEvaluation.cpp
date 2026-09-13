#include "WVDensityEventEvaluation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {
WVKernelStatus invalid(const char *message) {
  return {WVKernelStatusCode::invalidConfiguration, message};
}
bool positive(double value) {
  return std::isfinite(value) && value > 0.0;
}
bool overlaps(const void *a, std::size_t aBytes, const void *b, std::size_t bBytes) {
  if (!a || !b || !aBytes || !bBytes)
    return false;
  const auto x = reinterpret_cast<std::uintptr_t>(a);
  const auto y = reinterpret_cast<std::uintptr_t>(b);
  return x <= y ? y - x < aBytes : x - y < bBytes;
}
bool validField(WVDensityEventField field) {
  return field == WVDensityEventField::rhoNm ||
         field == WVDensityEventField::etaTrue || field == WVDensityEventField::ape;
}
} // namespace

WVKernelStatus WVDensityEventEvaluation::create(
    WVRealVolumeConstView density, WVDensityEventGeometry geometry,
    WVDensityDiagnosticContract contract, WVDensityEventEvaluation &output,
    const WVNoMotionRecoveryOptions &options) {
  const auto nx = density.shape.first, ny = density.shape.second;
  const auto nz = density.shape.third;
  if (!nx || !ny || nz < 2)
    return {WVKernelStatusCode::invalidShape,
            "Density event requires a nonempty volume with at least two heights."};
  const auto maximum = std::numeric_limits<std::size_t>::max();
  if (nx > maximum / ny || nx * ny > maximum / nz ||
      nx * ny * nz > maximum / sizeof(double))
    return {WVKernelStatusCode::sizeOverflow, "Density event volume size overflowed."};
  if (!density.data || !geometry.heights || !geometry.integrationWeights ||
      !geometry.initialProfile)
    return {WVKernelStatusCode::invalidPointer, "Density event binding is incomplete."};
  if (geometry.heights->size() != nz || geometry.integrationWeights->size() != nz ||
      geometry.initialProfile->size() != nz)
    return {WVKernelStatusCode::invalidShape, "Density event geometry has incompatible heights."};
  if (!positive(geometry.depth) || !positive(geometry.gravity) ||
      !positive(geometry.referenceDensity))
    return invalid("Density event depth, gravity and reference density must be positive and finite.");
  if ((contract.reference != WVNoMotionReference::actual &&
       contract.reference != WVNoMotionReference::initial) ||
      options.maximumEvaluations == 0 || !positive(options.gradientTolerance) ||
      !positive(options.stepTolerance) || !positive(options.relativeCostTolerance))
    return invalid("Density event reference or recovery controls are invalid.");
  for (std::size_t i = 0; i < nz; ++i) {
    const double height = (*geometry.heights)[i];
    const double initial = (*geometry.initialProfile)[i];
    if (!std::isfinite(height) || !std::isfinite(initial) || initial < 0.0 ||
        !positive((*geometry.integrationWeights)[i]) ||
        (i && (!(height > (*geometry.heights)[i - 1]) ||
               !(initial < (*geometry.initialProfile)[i - 1]))))
      return invalid("Density event geometry requires finite ordered heights and stable initial density.");
  }
  if(output.lowMemoryStorage_) {
    const auto demands=output.lowMemoryDemands_;
    const bool derived=(demands&(etaTrueDemand|apeDemand))!=0;
    const bool combined=(demands&(etaTrueDemand|apeDemand))==
        (etaTrueDemand|apeDemand);
    const bool actual=(demands&rhoNmDemand) ||
        (contract.reference==WVNoMotionReference::actual && derived);
    if((derived && output.lowPrimary_.capacity()<nx*ny*nz) ||
        (combined && output.lowSecondary_.capacity()<nx*ny*nz) ||
        (actual && output.actualProfile_.capacity()<nz))
      return invalid("Prepared low-memory density storage is incompatible with the event.");
  }
  // No current-field scan or scientific allocation until a demand is active.
  // Rebinding preserves capacities prepared by the service-owned event arena.
  output.resetRetainingCapacity();
  output.density_ = density;
  output.geometry_ = geometry;
  output.contract_ = contract;
  output.options_ = options;
  output.sampleCount_ = nx * ny * nz;
  output.initialized_ = true;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::reserveStorage(
    std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
    WVNoMotionReference reference) {
  if((demands&~(rhoNmDemand|etaTrueDemand|apeDemand))!=0)
    return invalid("Density event storage demand is invalid.");
  try {
    lowMemoryStorage_=false;
    lowMemoryDemands_=0;
    std::vector<double>{}.swap(lowPrimary_);
    std::vector<double>{}.swap(lowSecondary_);
    if((demands&rhoNmDemand) ||
        (reference==WVNoMotionReference::actual &&
          (demands&(etaTrueDemand|apeDemand))))
      actualProfile_.reserve(profileCount);
    if(demands&(etaTrueDemand|apeDemand))
      materialHeights_.reserve(sampleCount);
    if(demands&etaTrueDemand) etaTrue_.reserve(sampleCount);
    if(demands&apeDemand) ape_.reserve(sampleCount);
    account();
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    account();
    return {WVKernelStatusCode::allocationFailure,
        "Unable to prepare density event storage."};
  } catch(const std::length_error&) {
    account();
    return {WVKernelStatusCode::sizeOverflow,
        "Prepared density event storage exceeds vector capacity."};
  }
}

WVKernelStatus WVDensityEventEvaluation::reserveAdditionalStorage(
    std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
    WVNoMotionReference reference) {
  if((demands&~(rhoNmDemand|etaTrueDemand|apeDemand))!=0)
    return invalid("Density event storage demand is invalid.");
  if(!initialized_)
    return reserveStorage(sampleCount,profileCount,demands,reference);
  if(lowMemoryStorage_ || reference!=contract_.reference ||
      sampleCount!=sampleCount_ || !geometry_.heights ||
      profileCount!=geometry_.heights->size())
    return invalid("Active density storage extension does not match its immutable event binding.");
  try {
    if(((demands&rhoNmDemand) ||
        (reference==WVNoMotionReference::actual &&
         (demands&(etaTrueDemand|apeDemand)))) && !actualReady_)
      actualProfile_.reserve(profileCount);
    if((demands&etaTrueDemand) && !etaReady_) {
      materialHeights_.reserve(sampleCount);
      etaTrue_.reserve(sampleCount);
    }
    if((demands&apeDemand) && !apeReady_) {
      materialHeights_.reserve(sampleCount);
      ape_.reserve(sampleCount);
    }
    account();
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    account();
    return {WVKernelStatusCode::allocationFailure,
        "Unable to extend active density event storage."};
  } catch(const std::length_error&) {
    account();
    return {WVKernelStatusCode::sizeOverflow,
        "Active density event storage exceeds vector capacity."};
  }
}

WVKernelStatus WVDensityEventEvaluation::reserveLowMemoryStorage(
    std::size_t sampleCount,std::size_t profileCount,std::uint8_t demands,
    WVNoMotionReference reference) {
  if((demands&~(rhoNmDemand|etaTrueDemand|apeDemand))!=0)
    return invalid("Density event low-memory storage demand is invalid.");
  if(initialized_)
    return invalid("Cannot resize low-memory density storage during an event.");
  try {
    lowMemoryStorage_=true;
    lowMemoryDemands_=demands;
    std::vector<double>{}.swap(materialHeights_);
    std::vector<double>{}.swap(etaTrue_);
    std::vector<double>{}.swap(ape_);
    const bool derived=(demands&(etaTrueDemand|apeDemand))!=0;
    const bool actual=(demands&rhoNmDemand) ||
        (reference==WVNoMotionReference::actual && derived);
    if(actual)
      actualProfile_.reserve(profileCount);
    else std::vector<double>{}.swap(actualProfile_);
    if(derived)
      lowPrimary_.reserve(sampleCount);
    else std::vector<double>{}.swap(lowPrimary_);
    if((demands&(etaTrueDemand|apeDemand))==(etaTrueDemand|apeDemand))
      lowSecondary_.reserve(sampleCount);
    else std::vector<double>{}.swap(lowSecondary_);
    account();
    return WVKernelStatus::ok();
  } catch(const std::bad_alloc&) {
    account();
    return {WVKernelStatusCode::allocationFailure,
        "Unable to prepare low-memory density event storage."};
  } catch(const std::length_error&) {
    account();
    return {WVKernelStatusCode::sizeOverflow,
        "Prepared low-memory density storage exceeds vector capacity."};
  }
}

void WVDensityEventEvaluation::account() noexcept {
  metrics_.liveBytes = sizeof(double) *
      (actualProfile_.capacity() + materialHeights_.capacity() +
       etaTrue_.capacity() + ape_.capacity() + lowPrimary_.capacity() +
       lowSecondary_.capacity()) +
      profile_.retainedBytes() - sizeof(profile_);
  metrics_.highWaterBytes = std::max(metrics_.highWaterBytes, metrics_.liveBytes);
}

WVKernelStatus WVDensityEventEvaluation::recoverActual() {
  if (actualReady_)
    return WVKernelStatus::ok();
  ++metrics_.recoveryAttemptCount;
  const auto before = metrics_.liveBytes;
  const auto status = WVNoMotionProfileRecovery::recover(
      density_, *geometry_.integrationWeights, geometry_.depth,
      *geometry_.initialProfile, actualProfile_, recoveryReport_, options_);
  metrics_.recoveryWorkspaceHighWaterBytes =
      std::max(metrics_.recoveryWorkspaceHighWaterBytes, recoveryReport_.workspaceBytes);
  metrics_.highWaterBytes =
      std::max(metrics_.highWaterBytes, before + recoveryReport_.workspaceBytes);
  account();
  if (!status)
    return status;
  actualReady_ = true;
  ++metrics_.recoveryCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::prepareProfile() {
  if (profileReady_)
    return WVKernelStatus::ok();
  if (contract_.reference == WVNoMotionReference::actual) {
    const auto status = recoverActual();
    if (!status)
      return status;
  }
  ++metrics_.profileConstructionAttemptCount;
  const auto before = metrics_.liveBytes;
  std::size_t constructionBytes = 0;
  const auto status = WVNoMotionProfile::create(
      *geometry_.heights,
      contract_.reference == WVNoMotionReference::actual ? actualProfile_ :
                                                          *geometry_.initialProfile,
      profile_, &constructionBytes);
  metrics_.profileConstructionWorkspaceHighWaterBytes =
      std::max(metrics_.profileConstructionWorkspaceHighWaterBytes, constructionBytes);
  metrics_.highWaterBytes = std::max(metrics_.highWaterBytes, before + constructionBytes);
  account();
  if (!status)
    return status;
  profileReady_ = true;
  ++metrics_.profileConstructionCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::invert() {
  if (inverseReady_)
    return WVKernelStatus::ok();
  ++metrics_.inverseAttemptCount;
  auto& material=lowMemoryStorage_ ? lowPrimary_ : materialHeights_;
  material.resize(sampleCount_);
  account();
  for (std::size_t i = 0; i < sampleCount_; ++i) {
    ++metrics_.inverseSampleCount;
    const auto status = profile_.inverse(density_.data[i], material[i]);
    if (!status)
      return status;
  }
  inverseReady_ = true;
  ++metrics_.inversePassCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::formEta() {
  if (etaReady_)
    return WVKernelStatus::ok();
  auto* material=lowMemoryStorage_ ? lowPrimary_.data() : materialHeights_.data();
  double* destination=nullptr;
  if(lowMemoryStorage_) {
    const bool preserveMaterial=(lowMemoryDemands_&apeDemand)!=0 && !apeReady_;
    auto& values=preserveMaterial ? lowSecondary_ : lowPrimary_;
    values.resize(sampleCount_);
    if(!preserveMaterial) material=values.data();
    destination=values.data();
    etaLowSlot_=preserveMaterial ? 1 : 0;
    if(!preserveMaterial) inverseReady_=false;
  } else {
    etaTrue_.resize(sampleCount_);
    destination=etaTrue_.data();
  }
  account();
  const auto plane = density_.shape.first * density_.shape.second;
  for (std::size_t z = 0; z < density_.shape.third; ++z) {
    for (std::size_t point = 0; point < plane; ++point) {
      const auto i = z * plane + point;
      const double displacement = (*geometry_.heights)[z] - material[i];
      if (!std::isfinite(displacement))
        return {WVKernelStatusCode::numericalFailure, "Density event displacement overflowed."};
      destination[i] = displacement;
    }
  }
  etaReady_ = true;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::formAPE() {
  if (apeReady_)
    return WVKernelStatus::ok();
  ++metrics_.apeAttemptCount;
  auto* material=lowMemoryStorage_ ? lowPrimary_.data() : materialHeights_.data();
  double* destination=nullptr;
  if(lowMemoryStorage_) {
    lowPrimary_.resize(sampleCount_);
    material=lowPrimary_.data();
    destination=lowPrimary_.data();
    apeLowSlot_=0;
    inverseReady_=false;
  } else {
    ape_.resize(sampleCount_);
    destination=ape_.data();
  }
  account();
  const auto plane = density_.shape.first * density_.shape.second;
  for (std::size_t z = 0; z < density_.shape.third; ++z) {
    for (std::size_t point = 0; point < plane; ++point) {
      const auto i = z * plane + point;
      ++metrics_.apeSampleCount;
      const auto status = profile_.availablePotentialEnergy(
          (*geometry_.heights)[z], material[i], geometry_.gravity,
          geometry_.referenceDensity, destination[i]);
      if (!status)
        return status;
    }
  }
  apeReady_ = true;
  ++metrics_.apePassCount;
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::prepare(std::uint8_t demands) {
  if (demands == 0)
    return WVKernelStatus::ok();
  if ((demands & ~(rhoNmDemand | etaTrueDemand | apeDemand)) != 0)
    return invalid("Density event demand contains an unknown diagnostic.");
  if (!initialized_)
    return invalid("Density event is not initialized.");
  if(lowMemoryStorage_ && (demands&~lowMemoryDemands_)!=0)
    return invalid("Density event demand was not prepared for low-memory execution.");
  const std::size_t reused =
      ((demands & rhoNmDemand) && actualReady_ ? 1u : 0u) +
      ((demands & etaTrueDemand) && etaReady_ ? 1u : 0u) +
      ((demands & apeDemand) && apeReady_ ? 1u : 0u);
  try {
    if (demands & rhoNmDemand) {
      const auto status = recoverActual();
      if (!status)
        return status;
    }
    if (demands & (etaTrueDemand | apeDemand)) {
      auto status = prepareProfile();
      if (!status)
        return status;
      const bool needsEta=(demands&etaTrueDemand) && !etaReady_;
      const bool needsAPE=(demands&apeDemand) && !apeReady_;
      if(lowMemoryStorage_ && !inverseReady_ && (needsEta || needsAPE)) {
        // Preserve an earlier result that occupies the primary slot before it
        // is reused for a new inversion during a demand extension.
        if((etaReady_ && etaLowSlot_==0) || (apeReady_ && apeLowSlot_==0)) {
          lowSecondary_.resize(sampleCount_);
          std::copy_n(lowPrimary_.data(),sampleCount_,lowSecondary_.data());
          if(etaReady_ && etaLowSlot_==0) etaLowSlot_=1;
          if(apeReady_ && apeLowSlot_==0) apeLowSlot_=1;
          account();
        }
      }
      if(needsEta || needsAPE) {
        status = invert();
        if (!status)
          return status;
      }
      if (needsEta) {
        status = formEta();
        if (!status)
          return status;
      }
      if (needsAPE) {
        status = formAPE();
        if (!status)
          return status;
      }
    }
    ++metrics_.evaluationCount;
    metrics_.reuseCount += reused;
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    account();
    return {WVKernelStatusCode::allocationFailure, "Density event workspace allocation failed."};
  } catch (const std::length_error &) {
    account();
    return {WVKernelStatusCode::sizeOverflow, "Density event workspace exceeds vector capacity."};
  }
}

WVDensityEventView WVDensityEventEvaluation::view(WVDensityEventField field) const noexcept {
  switch (field) {
  case WVDensityEventField::rhoNm:
    return actualReady_ ? WVDensityEventView{actualProfile_.data(), actualProfile_.size()} :
                          WVDensityEventView{};
  case WVDensityEventField::etaTrue:
    if(!etaReady_) return {};
    if(lowMemoryStorage_) {
      const auto& values=etaLowSlot_==0 ? lowPrimary_ : lowSecondary_;
      return {values.data(),values.size()};
    }
    return {etaTrue_.data(), etaTrue_.size()};
  case WVDensityEventField::ape:
    if(!apeReady_) return {};
    if(lowMemoryStorage_) {
      const auto& values=apeLowSlot_==0 ? lowPrimary_ : lowSecondary_;
      return {values.data(),values.size()};
    }
    return {ape_.data(), ape_.size()};
  }
  return {};
}

WVDensityEventView WVDensityEventEvaluation::materialHeights() const noexcept {
  if(!inverseReady_) return {};
  const auto& values=lowMemoryStorage_ ? lowPrimary_ : materialHeights_;
  return {values.data(),values.size()};
}

WVKernelStatus WVDensityEventEvaluation::validateOutputs(
    const WVDensityEventOutput *outputs, std::size_t count, std::uint8_t &demands) const {
  demands = 0;
  if (!initialized_)
    return invalid("Density event is not initialized.");
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(*outputs))
    return {WVKernelStatusCode::sizeOverflow, "Density output descriptor size overflowed."};
  if (count && !outputs)
    return {WVKernelStatusCode::invalidPointer, "Density output descriptors are missing."};
  const auto overlapStatus = []() {
    return WVKernelStatus{WVKernelStatusCode::overlappingArrays,
                          "Density outputs overlap event inputs, prepared storage or one another."};
  };
  for (std::size_t i = 0; i < count; ++i) {
    const auto &output = outputs[i];
    if (!validField(output.field))
      return invalid("Density output field is unknown.");
    const auto expected = output.field == WVDensityEventField::rhoNm ?
                           density_.shape.third : sampleCount_;
    if (output.elementCount != expected)
      return {WVKernelStatusCode::invalidShape, "Density output has the wrong element count."};
    if (!output.data)
      return {WVKernelStatusCode::invalidPointer, "Density output storage is missing."};
    const auto bytes = expected * sizeof(double);
    if (overlaps(output.data, bytes, density_.data, sampleCount_ * sizeof(double)) ||
        overlaps(output.data, bytes, this, sizeof(*this)) ||
        overlaps(output.data, bytes, outputs, count * sizeof(*outputs)))
      return overlapStatus();
    const std::vector<double> *vectors[] = {
        geometry_.heights, geometry_.integrationWeights, geometry_.initialProfile,
        &actualProfile_, &materialHeights_, &etaTrue_, &ape_, &lowPrimary_,
        &lowSecondary_};
    for (const auto *values : vectors)
      if (overlaps(output.data, bytes, values->data(), values->capacity() * sizeof(double)))
        return overlapStatus();
    for (std::size_t previous = 0; previous < i; ++previous)
      if (overlaps(output.data, bytes, outputs[previous].data,
                   outputs[previous].elementCount * sizeof(double)))
        return overlapStatus();
    demands |= static_cast<std::uint8_t>(output.field);
  }
  return WVKernelStatus::ok();
}

WVKernelStatus WVDensityEventEvaluation::evaluate(
    const WVDensityEventOutput *outputs, std::size_t count) {
  if (count == 0)
    return WVKernelStatus::ok();
  std::uint8_t demands = 0;
  auto status = validateOutputs(outputs, count, demands);
  if (!status)
    return status;
  status = prepare(demands);
  if (!status)
    return status;
  // All validation and fallible calculations precede the first caller write.
  for (std::size_t i = 0; i < count; ++i) {
    const auto values = view(outputs[i].field);
    std::copy_n(values.data, values.elementCount, outputs[i].data);
  }
  return WVKernelStatus::ok();
}

void WVDensityEventEvaluation::release() noexcept {
  std::vector<double>{}.swap(actualProfile_);
  std::vector<double>{}.swap(materialHeights_);
  std::vector<double>{}.swap(etaTrue_);
  std::vector<double>{}.swap(ape_);
  std::vector<double>{}.swap(lowPrimary_);
  std::vector<double>{}.swap(lowSecondary_);
  profile_=WVNoMotionProfile{};
  density_ = {};
  geometry_ = {};
  sampleCount_ = 0;
  lowMemoryDemands_=0;
  etaLowSlot_=apeLowSlot_=-1;
  lowMemoryStorage_=false;
  initialized_ = false;
  actualReady_=profileReady_=inverseReady_=etaReady_=apeReady_=false;
  account();
}

void WVDensityEventEvaluation::resetRetainingCapacity() noexcept {
  actualProfile_.clear();
  materialHeights_.clear();
  etaTrue_.clear();
  ape_.clear();
  lowPrimary_.clear();
  lowSecondary_.clear();
  profile_=WVNoMotionProfile{};
  density_={};
  geometry_={};
  sampleCount_=0;
  initialized_=false;
  etaLowSlot_=apeLowSlot_=-1;
  actualReady_=profileReady_=inverseReady_=etaReady_=apeReady_=false;
  account();
}

void WVDensityEventEvaluation::discardDerived() noexcept {
  if(lowMemoryStorage_) {
    actualProfile_.clear();
    lowPrimary_.clear();
    lowSecondary_.clear();
  } else {
    std::vector<double>{}.swap(actualProfile_);
    std::vector<double>{}.swap(materialHeights_);
    std::vector<double>{}.swap(etaTrue_);
    std::vector<double>{}.swap(ape_);
  }
  profile_ = WVNoMotionProfile{};
  etaLowSlot_=apeLowSlot_=-1;
  actualReady_ = profileReady_ = inverseReady_ = etaReady_ = apeReady_ = false;
  account();
}

} // namespace wavevortex::runtime::detail
