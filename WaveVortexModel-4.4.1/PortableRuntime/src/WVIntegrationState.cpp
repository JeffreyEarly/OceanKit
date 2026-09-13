#include "WaveVortexRuntime/WVIntegrationState.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <set>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}
bool sameShape(WVShape2D first, WVShape2D second) noexcept {
  return first.rows == second.rows && first.columns == second.columns;
}

std::size_t checkedCount(const std::vector<std::size_t> &dimensions,
                         bool &invalidValue) noexcept {
  std::size_t result = 1;
  invalidValue = dimensions.empty();
  for (const auto dimension : dimensions) {
    if (dimension == 0 ||
        result > std::numeric_limits<std::size_t>::max() / dimension) {
      invalidValue = true;
      return 0;
    }
    result *= dimension;
  }
  return result;
}

WVTransformStateDescription legacyDescription(WVShape2D shape) {
  WVTransformStateDescription result;
  result.transformIdentifier = "WVTransformConstantStratification";
  for (const char *identifier : {"Ap", "Am", "A0"})
    result.coefficientFamilies.push_back(
        {identifier, {shape.rows, shape.columns},
         WVToleranceKind::coefficientEnergyScaled});
  return result;
}

bool legacyTriple(
    const std::vector<WVCoefficientFamilyLayout> &families) noexcept {
  if (families.size() != 3)
    return false;
  const char *identifiers[] = {"Ap", "Am", "A0"};
  for (std::size_t index = 0; index < 3; ++index)
    if (families[index].identifier != identifiers[index] ||
        families[index].spectralDimensions.size() != 2 ||
        families[index].spectralDimensions !=
            families.front().spectralDimensions)
      return false;
  return true;
}

} // namespace

WVKernelStatus WVIntegrationStateLayout::initializeDescription(
    WVTransformStateDescription description,
    const WVPortableObserverRecord *source,
    WVIntegrationStateLayout &candidate) {
  if (description.transformIdentifier.empty() ||
      description.coefficientFamilies.empty())
    return invalid("A transform state description requires a nonempty "
                   "identity and at least one coefficient family.");
  if (!description.spatialDimensions.empty()) {
    bool badSpatialCount = false;
    (void)checkedCount(description.spatialDimensions, badSpatialCount);
    if (badSpatialCount)
      return {WVKernelStatusCode::invalidShape,
              "Transform spatial dimensions must be nonzero."};
  }

  std::set<std::string> identifiers;
  candidate.transformIdentifier_ = std::move(description.transformIdentifier);
  candidate.spatialDimensions_ = std::move(description.spatialDimensions);
  for (auto &descriptionFamily : description.coefficientFamilies) {
    if (descriptionFamily.identifier.empty() ||
        !identifiers.insert(descriptionFamily.identifier).second)
      return invalid("Coefficient-family identifiers must be unique and "
                     "nonempty.");
    bool badCount = false;
    const auto count =
        checkedCount(descriptionFamily.spectralDimensions, badCount);
    if (badCount)
      return {WVKernelStatusCode::sizeOverflow,
              "Coefficient-family dimensions are invalid."};
    if (candidate.coefficientElementCount_ >
        std::numeric_limits<std::size_t>::max() - count)
      return {WVKernelStatusCode::sizeOverflow,
              "Coefficient-family storage overflows size_t."};
    candidate.coefficientFamilies_.push_back(
        {std::move(descriptionFamily.identifier),
         std::move(descriptionFamily.spectralDimensions),
         descriptionFamily.toleranceKind, count,
         candidate.coefficientElementCount_});
    candidate.coefficientElementCount_ += count;
  }
  candidate.hasLegacyCoefficientTriple_ =
      legacyTriple(candidate.coefficientFamilies_);
  if (candidate.hasLegacyCoefficientTriple_)
    candidate.coefficientShape_ =
        {candidate.coefficientFamilies_.front().spectralDimensions[0],
         candidate.coefficientFamilies_.front().spectralDimensions[1]};

  if (source == nullptr) {
    for (const auto &family : candidate.coefficientFamilies_) {
      WVStateBlockRecord record;
      record.identifier = family.identifier;
      record.scalarType = WVStateScalarType::complex64;
      record.dimensions = family.spectralDimensions;
      record.toleranceKind = family.toleranceKind;
      record.ownership = WVStateOwnership::integratorOwned;
      record.restartRequirement =
          WVRestartRequirement::requiredDynamicState;
      candidate.stateBlockRecords_.push_back(std::move(record));
    }
    return WVKernelStatus::ok();
  }

  candidate.stateBlockRecords_ = source->stateBlocks;
  candidate.observerRecords_ = source->observers;
  std::set<std::string> observedFamilies;
  for (const auto &record : source->stateBlocks) {
    bool badCount = false;
    const auto count = checkedCount(record.dimensions, badCount);
    if (badCount)
      return {WVKernelStatusCode::sizeOverflow,
              "Integration state-block element count is invalid."};
    const auto family = std::find_if(
        candidate.coefficientFamilies_.begin(),
        candidate.coefficientFamilies_.end(), [&](const auto &value) {
          return value.identifier == record.identifier;
        });
    if (family != candidate.coefficientFamilies_.end()) {
      if (!observedFamilies.insert(record.identifier).second ||
          record.scalarType != WVStateScalarType::complex64 ||
          record.ownership != WVStateOwnership::integratorOwned ||
          record.dimensions != family->spectralDimensions ||
          record.toleranceKind != family->toleranceKind)
        return {WVKernelStatusCode::invalidShape,
                "Canonical coefficient blocks must match the resolved "
                "transform-family contract."};
      continue;
    }
    if (record.ownership != WVStateOwnership::integratorOwned)
      continue;
    WVAdditionalStateBlockLayout block{
        record.identifier, record.scalarType, record.dimensions,
        record.toleranceKind, record.absoluteTolerance, record.ownership,
        record.restartRequirement, count, 0};
    auto &total = record.scalarType == WVStateScalarType::real64
                      ? candidate.realElementCount_
                      : candidate.complexElementCount_;
    if (total > std::numeric_limits<std::size_t>::max() - count)
      return {WVKernelStatusCode::sizeOverflow,
              "Integration state storage overflows size_t."};
    block.scalarOffset = total;
    total += count;
    candidate.additionalBlocks_.push_back(std::move(block));
  }
  if (observedFamilies != identifiers)
    return invalid("Integration layout is missing one or more resolved "
                   "coefficient families.");
  return WVKernelStatus::ok();
}

WVKernelStatus WVIntegrationStateLayout::createCoefficientOnly(
    WVShape2D coefficientShape, WVIntegrationStateLayout &layout) {
  if (coefficientShape.rows == 0 || coefficientShape.columns == 0)
    return invalid("Integration coefficient shape must be nonzero.");
  auto description = legacyDescription(coefficientShape);
  return createCoefficientOnly(std::move(description), layout);
}

WVKernelStatus
WVIntegrationStateLayout::create(WVShape2D coefficientShape,
                               const WVPortableObserverDescriptor &descriptor,
                               WVIntegrationStateLayout &layout) {
  return create(coefficientShape, descriptor.record(), layout);
}

WVKernelStatus
WVIntegrationStateLayout::create(WVShape2D coefficientShape,
                               const WVPortableObserverRecord &source,
                               WVIntegrationStateLayout &layout) {
  if (coefficientShape.rows == 0 || coefficientShape.columns == 0)
    return invalid("Integration coefficient shape must be nonzero.");
  auto description = legacyDescription(coefficientShape);
  return create(std::move(description), source, layout);
}

WVKernelStatus WVIntegrationStateLayout::createCoefficientOnly(
    WVTransformStateDescription description,
    WVIntegrationStateLayout &layout) {
  try {
    WVIntegrationStateLayout candidate;
    auto result =
        initializeDescription(std::move(description), nullptr, candidate);
    if (!result)
      return result;
    layout = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Integration state-layout allocation failed."};
  }
}

WVKernelStatus WVIntegrationStateLayout::create(
    WVTransformStateDescription description,
    const WVPortableObserverDescriptor &descriptor,
    WVIntegrationStateLayout &layout) {
  return create(std::move(description), descriptor.record(), layout);
}

WVKernelStatus WVIntegrationStateLayout::create(
    WVTransformStateDescription description,
    const WVPortableObserverRecord &source,
    WVIntegrationStateLayout &layout) {
  try {
    WVIntegrationStateLayout candidate;
    auto result =
        initializeDescription(std::move(description), &source, candidate);
    if (!result)
      return result;
    layout = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Integration state-layout allocation failed."};
  }
}

std::size_t WVIntegrationStateLayout::persistentBytes() const noexcept {
  std::size_t bytes =
      transformIdentifier_.capacity() +
      spatialDimensions_.capacity() * sizeof(std::size_t) +
      coefficientFamilies_.capacity() * sizeof(WVCoefficientFamilyLayout) +
      additionalBlocks_.capacity() * sizeof(WVAdditionalStateBlockLayout) +
      stateBlockRecords_.capacity() * sizeof(WVStateBlockRecord) +
      observerRecords_.capacity() * sizeof(WVObserverRecord);
  for (const auto &family : coefficientFamilies_)
    bytes += family.identifier.capacity() +
             family.spectralDimensions.capacity() * sizeof(std::size_t);
  for (const auto &block : additionalBlocks_)
    bytes += block.identifier.capacity() +
             block.dimensions.capacity() * sizeof(std::size_t);
  for (const auto &block : stateBlockRecords_)
    bytes += block.identifier.capacity() +
             block.dimensions.capacity() * sizeof(std::size_t);
  for (const auto &observer : observerRecords_) {
    bytes += observer.identifier.capacity() + observer.name.capacity() +
             observer.typeIdentifier.capacity() +
             observer.configuration.persistentBytes() -
                 sizeof(WVPortableTypedRecord) +
             observer.stateBlockIdentifiers.capacity() * sizeof(std::string) +
             observer.fieldNames.capacity() * sizeof(std::string) +
             (observer.x.capacity() + observer.y.capacity() +
              observer.z.capacity()) *
                 sizeof(double);
    for (const auto &identifier : observer.stateBlockIdentifiers)
      bytes += identifier.capacity();
    for (const auto &name : observer.fieldNames)
      bytes += name.capacity();
  }
  return bytes;
}

WVKernelStatus
WVCoefficientStateStorage::initialize(const WVIntegrationStateLayout &layout) {
  try {
    storage_.assign(layout.coefficientElementCount(), WVComplex64{});
    mutableViews_.clear();
    constViews_.clear();
    mutableViews_.reserve(layout.coefficientFamilyCount());
    constViews_.reserve(layout.coefficientFamilyCount());
    for (const auto &family : layout.coefficientFamilies()) {
      auto *data = storage_.data() + family.scalarOffset;
      mutableViews_.push_back({&family, data});
      constViews_.push_back({&family, data});
    }
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    clear();
    return {WVKernelStatusCode::allocationFailure,
            "Coefficient state allocation failed."};
  }
}

std::size_t WVCoefficientStateStorage::capacityBytes() const noexcept {
  return storage_.capacity() * sizeof(WVComplex64) +
         mutableViews_.capacity() * sizeof(WVCoefficientFamilyView) +
         constViews_.capacity() * sizeof(WVCoefficientFamilyConstView);
}

void WVCoefficientStateStorage::clear() noexcept {
  storage_.clear();
  mutableViews_.clear();
  constViews_.clear();
}

std::size_t WVTransformStateCheckpoint::persistentBytes() const noexcept {
  std::size_t bytes = sizeof(*this) + transformIdentifier.capacity() +
                      spatialDimensions.capacity() * sizeof(std::size_t) +
                      coefficientFamilies.capacity() *
                          sizeof(WVCoefficientFamilyCheckpoint);
  for (const auto &family : coefficientFamilies)
    bytes += family.identifier.capacity() +
             family.spectralDimensions.capacity() * sizeof(std::size_t) +
             family.values.capacity() * sizeof(WVComplex64);
  return bytes;
}

WVKernelStatus captureTransformStateCheckpoint(
    const WVIntegrationStateLayout &layout,
    const WVIntegrationState &state,
    WVTransformStateCheckpoint &checkpoint) {
  auto result = validateIntegrationState(layout, state);
  if (!result)
    return result;
  try {
    WVTransformStateCheckpoint candidate;
    candidate.transformIdentifier = layout.transformIdentifier();
    candidate.spatialDimensions = layout.spatialDimensions();
    candidate.t = state.waveVortex.t;
    candidate.t0 = state.waveVortex.t0;
    candidate.coefficientFamilies.reserve(layout.coefficientFamilyCount());
    for (std::size_t index = 0; index < layout.coefficientFamilyCount();
         ++index) {
      const auto &metadata = layout.coefficientFamilies()[index];
      const auto view = coefficientFamilyView(layout, state, index);
      WVCoefficientFamilyCheckpoint family;
      family.identifier = metadata.identifier;
      family.spectralDimensions = metadata.spectralDimensions;
      family.values.assign(view.data, view.data + metadata.elementCount);
      candidate.coefficientFamilies.push_back(std::move(family));
    }
    checkpoint = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Transform checkpoint allocation failed."};
  }
}

WVKernelStatus restoreTransformStateCheckpoint(
    const WVTransformStateCheckpoint &checkpoint,
    const WVIntegrationStateLayout &layout,
    WVCoefficientStateStorage &storage,
    WVMutableIntegrationState &state) {
  if (checkpoint.transformIdentifier != layout.transformIdentifier() ||
      (!checkpoint.spatialDimensions.empty() &&
       !layout.spatialDimensions().empty() &&
       checkpoint.spatialDimensions != layout.spatialDimensions()) ||
      checkpoint.coefficientFamilies.size() !=
          layout.coefficientFamilyCount() ||
      !std::isfinite(checkpoint.t) || !std::isfinite(checkpoint.t0))
    return {WVKernelStatusCode::invalidShape,
            "Transform checkpoint identity, rank, or time does not match its "
            "resolved layout."};
  for (std::size_t index = 0; index < layout.coefficientFamilyCount();
       ++index) {
    const auto &expected = layout.coefficientFamilies()[index];
    const auto &actual = checkpoint.coefficientFamilies[index];
    if (actual.identifier != expected.identifier ||
        actual.spectralDimensions != expected.spectralDimensions ||
        actual.values.size() != expected.elementCount)
      return {WVKernelStatusCode::invalidShape,
              "Transform checkpoint coefficient-family identity or rank "
              "does not match its resolved layout."};
  }
  WVCoefficientStateStorage candidate;
  auto result = candidate.initialize(layout);
  if (!result)
    return result;
  for (std::size_t index = 0; index < layout.coefficientFamilyCount();
       ++index)
    std::copy(checkpoint.coefficientFamilies[index].values.begin(),
              checkpoint.coefficientFamilies[index].values.end(),
              candidate.mutableFamilies()[index].data);
  storage = std::move(candidate);
  state = {};
  state.waveVortex.t = checkpoint.t;
  state.waveVortex.t0 = checkpoint.t0;
  state.coefficientFamilies = storage.mutableFamilies();
  state.coefficientFamilyCount = storage.familyCount();
  if (layout.hasLegacyCoefficientTriple()) {
    const auto shape = layout.coefficientShape();
    auto *families = storage.mutableFamilies();
    state.waveVortex.coefficients =
        {{families[0].data, shape}, {families[1].data, shape},
         {families[2].data, shape}};
  }
  return WVKernelStatus::ok();
}

WVKernelStatus
WVAdditionalStateStorage::initialize(const WVIntegrationStateLayout &layout) {
  try {
    realStorage_.assign(layout.realElementCount(), 0.0);
    complexStorage_.assign(layout.complexElementCount(), WVComplex64{});
    mutableViews_.clear();
    constViews_.clear();
    mutableViews_.reserve(layout.additionalBlocks().size());
    constViews_.reserve(layout.additionalBlocks().size());
    for (const auto &block : layout.additionalBlocks()) {
      auto *real = block.scalarType == WVStateScalarType::real64
                       ? realStorage_.data() + block.scalarOffset
                       : nullptr;
      auto *complex = block.scalarType == WVStateScalarType::complex64
                          ? complexStorage_.data() + block.scalarOffset
                          : nullptr;
      mutableViews_.push_back({&block, real, complex});
      constViews_.push_back({&block, real, complex});
    }
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    clear();
    return {WVKernelStatusCode::allocationFailure,
            "Additional state storage allocation failed."};
  }
}

std::size_t WVAdditionalStateStorage::capacityBytes() const noexcept {
  return realStorage_.capacity() * sizeof(double) +
         complexStorage_.capacity() * sizeof(WVComplex64) +
         mutableViews_.capacity() * sizeof(WVAdditionalStateBlockView) +
         constViews_.capacity() * sizeof(WVAdditionalStateBlockConstView);
}

void WVAdditionalStateStorage::clear() noexcept {
  realStorage_.clear();
  complexStorage_.clear();
  mutableViews_.clear();
  constViews_.clear();
}

WVIntegrationState
integrationConstView(const WVMutableIntegrationState &state,
                   std::vector<WVAdditionalStateBlockConstView> &blockViews) {
  blockViews.clear();
  blockViews.reserve(state.additionalBlockCount);
  for (std::size_t index = 0; index < state.additionalBlockCount; ++index)
    blockViews.push_back({state.additionalBlocks[index].layout,
                          state.additionalBlocks[index].realData,
                          state.additionalBlocks[index].complexData});
  return {state.waveVortex.view(), blockViews.data(), blockViews.size()};
}

WVIntegrationState integrationConstView(
    const WVMutableIntegrationState &state,
    std::vector<WVCoefficientFamilyConstView> &coefficientViews,
    std::vector<WVAdditionalStateBlockConstView> &blockViews) {
  coefficientViews.clear();
  coefficientViews.reserve(state.coefficientFamilyCount);
  for (std::size_t index = 0; index < state.coefficientFamilyCount; ++index)
    coefficientViews.push_back({state.coefficientFamilies[index].layout,
                                state.coefficientFamilies[index].data});
  blockViews.clear();
  blockViews.reserve(state.additionalBlockCount);
  for (std::size_t index = 0; index < state.additionalBlockCount; ++index)
    blockViews.push_back({state.additionalBlocks[index].layout,
                          state.additionalBlocks[index].realData,
                          state.additionalBlocks[index].complexData});
  return {state.waveVortex.view(), blockViews.data(), blockViews.size(),
          coefficientViews.data(), coefficientViews.size()};
}

WVCoefficientFamilyConstView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    const WVIntegrationState &state, std::size_t family) noexcept {
  if (family >= layout.coefficientFamilyCount())
    return {};
  if (state.coefficientFamilies != nullptr &&
      state.coefficientFamilyCount == layout.coefficientFamilyCount())
    return state.coefficientFamilies[family];
  if (!layout.hasLegacyCoefficientTriple())
    return {};
  const WVComplexConstView views[] = {state.waveVortex.coefficients.Ap,
                                      state.waveVortex.coefficients.Am,
                                      state.waveVortex.coefficients.A0};
  return {&layout.coefficientFamilies()[family], views[family].data};
}

WVCoefficientFamilyView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    WVMutableIntegrationState &state, std::size_t family) noexcept {
  if (family >= layout.coefficientFamilyCount())
    return {};
  if (state.coefficientFamilies != nullptr &&
      state.coefficientFamilyCount == layout.coefficientFamilyCount())
    return state.coefficientFamilies[family];
  if (!layout.hasLegacyCoefficientTriple())
    return {};
  WVComplexView views[] = {state.waveVortex.coefficients.Ap,
                           state.waveVortex.coefficients.Am,
                           state.waveVortex.coefficients.A0};
  return {&layout.coefficientFamilies()[family], views[family].data};
}

WVCoefficientFamilyView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    WVIntegrationFlux &flux, std::size_t family) noexcept {
  if (family >= layout.coefficientFamilyCount())
    return {};
  if (flux.coefficientFamilies != nullptr &&
      flux.coefficientFamilyCount == layout.coefficientFamilyCount())
    return flux.coefficientFamilies[family];
  if (!layout.hasLegacyCoefficientTriple())
    return {};
  WVComplexView views[] = {flux.waveVortex.Fp, flux.waveVortex.Fm,
                           flux.waveVortex.F0};
  return {&layout.coefficientFamilies()[family], views[family].data};
}

namespace {

template <typename State>
WVKernelStatus validateCoefficientFamilies(
    const WVIntegrationStateLayout &layout, State &state) {
  if (state.coefficientFamilyCount != 0) {
    if (state.coefficientFamilyCount != layout.coefficientFamilyCount() ||
        state.coefficientFamilies == nullptr)
      return {WVKernelStatusCode::invalidShape,
              "Coefficient-family count does not match its transform "
              "layout."};
    for (std::size_t index = 0; index < state.coefficientFamilyCount; ++index)
      if (state.coefficientFamilies[index].layout !=
              &layout.coefficientFamilies()[index] ||
          state.coefficientFamilies[index].data == nullptr)
        return {WVKernelStatusCode::invalidPointer,
                "Coefficient-family order or storage does not match its "
                "frozen transform layout."};
    return WVKernelStatus::ok();
  }
  if (!layout.hasLegacyCoefficientTriple())
    return {WVKernelStatusCode::invalidShape,
            "Transform-neutral coefficient-family views are required for "
            "this transform layout."};
  const auto shape = layout.coefficientShape();
  const auto &c = state.waveVortex.coefficients;
  if (c.Ap.data == nullptr || c.Am.data == nullptr || c.A0.data == nullptr ||
      !sameShape(c.Ap.shape, shape) || !sameShape(c.Am.shape, shape) ||
      !sameShape(c.A0.shape, shape))
    return {WVKernelStatusCode::invalidShape,
            "Legacy integration coefficients must use the resolved "
            "[Nj,Nkl] shape."};
  return WVKernelStatus::ok();
}

template <typename State>
WVKernelStatus validateAdditionalBlocks(
    const WVIntegrationStateLayout &layout, const State &state) {
  if (state.additionalBlockCount != layout.additionalBlocks().size() ||
      (state.additionalBlockCount && state.additionalBlocks == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Integration state-block count does not match its layout."};
  for (std::size_t index = 0; index < state.additionalBlockCount; ++index) {
    const auto &expected = layout.additionalBlocks()[index];
    const auto &actual = state.additionalBlocks[index];
    if (actual.layout != &expected)
      return invalid(
          "Integration state-block order does not match its frozen layout.");
    if ((expected.scalarType == WVStateScalarType::real64 &&
         actual.realData == nullptr) ||
        (expected.scalarType == WVStateScalarType::complex64 &&
         actual.complexData == nullptr))
      return {WVKernelStatusCode::invalidPointer,
              "Integration state block has a null data pointer."};
  }
  return WVKernelStatus::ok();
}

} // namespace

WVKernelStatus validateIntegrationState(const WVIntegrationStateLayout &layout,
                                      const WVIntegrationState &state) {
  if (!std::isfinite(state.waveVortex.t) || !std::isfinite(state.waveVortex.t0))
    return invalid("Integration state times must be finite.");
  auto result = validateCoefficientFamilies(layout, state);
  return result ? validateAdditionalBlocks(layout, state) : result;
}

WVKernelStatus
validateMutableIntegrationState(const WVIntegrationStateLayout &layout,
                              const WVMutableIntegrationState &state) {
  if (!std::isfinite(state.waveVortex.t) ||
      !std::isfinite(state.waveVortex.t0))
    return invalid("Integration state times must be finite.");
  auto result = validateCoefficientFamilies(layout, state);
  return result ? validateAdditionalBlocks(layout, state) : result;
}

bool sameIntegrationStateLayout(const WVIntegrationStateLayout &first,
                                const WVIntegrationStateLayout &second) noexcept {
  if (first.transformIdentifier() != second.transformIdentifier() ||
      (!first.spatialDimensions().empty() &&
       !second.spatialDimensions().empty() &&
       first.spatialDimensions() != second.spatialDimensions()) ||
      first.coefficientFamilies().size() !=
          second.coefficientFamilies().size() ||
      first.stateBlockRecords().size() != second.stateBlockRecords().size() ||
      first.observerRecords().size() != second.observerRecords().size())
    return false;
  for (std::size_t index = 0; index < first.coefficientFamilies().size();
       ++index) {
    const auto &a = first.coefficientFamilies()[index];
    const auto &b = second.coefficientFamilies()[index];
    if (a.identifier != b.identifier ||
        a.spectralDimensions != b.spectralDimensions ||
        a.toleranceKind != b.toleranceKind)
      return false;
  }
  for (std::size_t i = 0; i < first.stateBlockRecords().size(); ++i) {
    const auto &a = first.stateBlockRecords()[i];
    const auto &b = second.stateBlockRecords()[i];
    if (a.identifier != b.identifier || a.scalarType != b.scalarType ||
        a.dimensions != b.dimensions || a.toleranceKind != b.toleranceKind ||
        a.absoluteTolerance != b.absoluteTolerance ||
        a.ownership != b.ownership ||
        a.restartRequirement != b.restartRequirement)
      return false;
  }
  const auto &aObservers = first.observerRecords();
  const auto &bObservers = second.observerRecords();
  for (std::size_t i = 0; i < aObservers.size(); ++i) {
    const auto &a = aObservers[i];
    const auto &b = bObservers[i];
    if (a.identifier != b.identifier || a.name != b.name ||
        a.typeIdentifier != b.typeIdentifier ||
        a.contractVersion != b.contractVersion ||
        a.stateBlockIdentifiers != b.stateBlockIdentifiers ||
        a.fieldNames != b.fieldNames || a.x != b.x || a.y != b.y ||
        a.z != b.z || a.isXYOnly != b.isXYOnly ||
        a.shouldAntialias != b.shouldAntialias ||
        a.advectionInterpolation != b.advectionInterpolation ||
        a.trackedFieldInterpolation != b.trackedFieldInterpolation ||
        a.horizontalAbsoluteTolerance != b.horizontalAbsoluteTolerance ||
        a.verticalAbsoluteTolerance != b.verticalAbsoluteTolerance ||
        a.outputScale != b.outputScale || a.outputOffset != b.outputOffset)
      return false;
  }
  return true;
}

} // namespace wavevortex::runtime
