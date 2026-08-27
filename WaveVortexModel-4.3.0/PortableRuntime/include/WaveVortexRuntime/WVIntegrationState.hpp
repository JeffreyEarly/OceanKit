#pragma once

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace wavevortex::runtime {

struct WVAdditionalStateBlockLayout {
  std::string identifier;
  WVStateScalarType scalarType = WVStateScalarType::real64;
  std::vector<std::size_t> dimensions;
  WVToleranceKind toleranceKind = WVToleranceKind::uniformAbsolute;
  double absoluteTolerance = 0.0;
  WVStateOwnership ownership = WVStateOwnership::integratorOwned;
  WVRestartRequirement restartRequirement =
      WVRestartRequirement::requiredDynamicState;
  std::size_t elementCount = 0;
  std::size_t scalarOffset = 0;
};

// Immutable transform-owned description of one canonical coefficient family.
// Spectral dimensions use the transform's natural MATLAB order; families need
// not share a rank or shape. scalarOffset addresses the packed complex
// integration workspace and is resolved once during preflight.
struct WVCoefficientFamilyLayout {
  std::string identifier;
  std::vector<std::size_t> spectralDimensions;
  WVToleranceKind toleranceKind = WVToleranceKind::coefficientEnergyScaled;
  std::size_t elementCount = 0;
  std::size_t scalarOffset = 0;
};

struct WVCoefficientFamilyDescription {
  std::string identifier;
  std::vector<std::size_t> spectralDimensions;
  WVToleranceKind toleranceKind = WVToleranceKind::coefficientEnergyScaled;
};

// Allocation-light transform identity and rank contract resolved before any
// state-sized storage is loaded or allocated.
struct WVTransformStateDescription {
  std::string transformIdentifier;
  std::vector<std::size_t> spatialDimensions;
  std::vector<WVCoefficientFamilyDescription> coefficientFamilies;
  bool supportsFixedTimeStepSelection = false;
};

// Resolves identifiers and offsets once so numerical stages use ordered views
// without virtual dispatch or repeated name lookup. The originating state and
// observer records remain available for exact cross-boundary compatibility
// checks.
class WVIntegrationStateLayout final {
public:
  static WVKernelStatus createCoefficientOnly(WVShape2D coefficientShape,
                                              WVIntegrationStateLayout &layout);
  static WVKernelStatus create(WVShape2D coefficientShape,
                               const WVPortableObserverDescriptor &descriptor,
                               WVIntegrationStateLayout &layout);
  static WVKernelStatus create(WVShape2D coefficientShape,
                               const WVPortableObserverRecord &record,
                               WVIntegrationStateLayout &layout);
  static WVKernelStatus createCoefficientOnly(
      WVTransformStateDescription description,
      WVIntegrationStateLayout &layout);
  static WVKernelStatus create(
      WVTransformStateDescription description,
      const WVPortableObserverDescriptor &descriptor,
      WVIntegrationStateLayout &layout);
  static WVKernelStatus create(
      WVTransformStateDescription description,
      const WVPortableObserverRecord &record,
      WVIntegrationStateLayout &layout);
  const std::string &transformIdentifier() const noexcept {
    return transformIdentifier_;
  }
  const std::vector<std::size_t> &spatialDimensions() const noexcept {
    return spatialDimensions_;
  }
  const std::vector<WVCoefficientFamilyLayout> &
  coefficientFamilies() const noexcept {
    return coefficientFamilies_;
  }
  std::size_t coefficientFamilyCount() const noexcept {
    return coefficientFamilies_.size();
  }
  std::size_t coefficientElementCount() const noexcept {
    return coefficientElementCount_;
  }
  bool hasLegacyCoefficientTriple() const noexcept {
    return hasLegacyCoefficientTriple_;
  }
  // Stable v1 compatibility view. For the legacy Ap/Am/A0 contract this is
  // their common [Nj,Nkl] shape. Transform-neutral code uses
  // coefficientFamilies() instead.
  WVShape2D coefficientShape() const noexcept { return coefficientShape_; }
  const std::vector<WVAdditionalStateBlockLayout> &
  additionalBlocks() const noexcept {
    return additionalBlocks_;
  }
  const std::vector<WVStateBlockRecord> &stateBlockRecords() const noexcept {
    return stateBlockRecords_;
  }
  const std::vector<WVObserverRecord> &observerRecords() const noexcept {
    return observerRecords_;
  }
  std::size_t realElementCount() const noexcept { return realElementCount_; }
  std::size_t complexElementCount() const noexcept {
    return complexElementCount_;
  }
  std::size_t integratedScalarCount() const noexcept {
    return 2 * coefficientElementCount_ + realElementCount_ +
           2 * complexElementCount_;
  }
  std::size_t persistentBytes() const noexcept;

private:
  static WVKernelStatus initializeDescription(
      WVTransformStateDescription description,
      const WVPortableObserverRecord *source,
      WVIntegrationStateLayout &candidate);
  std::string transformIdentifier_;
  std::vector<std::size_t> spatialDimensions_;
  std::vector<WVCoefficientFamilyLayout> coefficientFamilies_;
  std::size_t coefficientElementCount_ = 0;
  bool hasLegacyCoefficientTriple_ = false;
  WVShape2D coefficientShape_;
  std::vector<WVAdditionalStateBlockLayout> additionalBlocks_;
  std::vector<WVStateBlockRecord> stateBlockRecords_;
  std::vector<WVObserverRecord> observerRecords_;
  std::size_t realElementCount_ = 0;
  std::size_t complexElementCount_ = 0;
};

struct WVCoefficientFamilyConstView {
  const WVCoefficientFamilyLayout *layout = nullptr;
  const WVComplex64 *data = nullptr;
};

struct WVCoefficientFamilyView {
  const WVCoefficientFamilyLayout *layout = nullptr;
  WVComplex64 *data = nullptr;
};

struct WVAdditionalStateBlockConstView {
  const WVAdditionalStateBlockLayout *layout = nullptr;
  const double *realData = nullptr;
  const WVComplex64 *complexData = nullptr;
};

struct WVAdditionalStateBlockView {
  const WVAdditionalStateBlockLayout *layout = nullptr;
  double *realData = nullptr;
  WVComplex64 *complexData = nullptr;
};

struct WVIntegrationState {
  WVState waveVortex;
  const WVAdditionalStateBlockConstView *additionalBlocks = nullptr;
  std::size_t additionalBlockCount = 0;
  const WVCoefficientFamilyConstView *coefficientFamilies = nullptr;
  std::size_t coefficientFamilyCount = 0;
};

struct WVMutableIntegrationState {
  WVMutableState waveVortex;
  WVAdditionalStateBlockView *additionalBlocks = nullptr;
  std::size_t additionalBlockCount = 0;
  WVCoefficientFamilyView *coefficientFamilies = nullptr;
  std::size_t coefficientFamilyCount = 0;
};

struct WVIntegrationFlux {
  WVFlux waveVortex;
  WVAdditionalStateBlockView *additionalBlocks = nullptr;
  std::size_t additionalBlockCount = 0;
  WVCoefficientFamilyView *coefficientFamilies = nullptr;
  std::size_t coefficientFamilyCount = 0;
};

// Owns exactly the coefficient families declared by one transform layout.
// An A0-only transform therefore owns one compact allocation and no dummy Ap
// or Am storage.
class WVCoefficientStateStorage final {
public:
  WVKernelStatus initialize(const WVIntegrationStateLayout &layout);
  WVCoefficientFamilyView *mutableFamilies() noexcept {
    return mutableViews_.data();
  }
  const WVCoefficientFamilyConstView *constFamilies() const noexcept {
    return constViews_.data();
  }
  std::size_t familyCount() const noexcept { return mutableViews_.size(); }
  std::size_t capacityBytes() const noexcept;
  void clear() noexcept;

private:
  std::vector<WVComplex64> storage_;
  std::vector<WVCoefficientFamilyView> mutableViews_;
  std::vector<WVCoefficientFamilyConstView> constViews_;
};

struct WVCoefficientFamilyCheckpoint {
  std::string identifier;
  std::vector<std::size_t> spectralDimensions;
  std::vector<WVComplex64> values;
};

// Transform-neutral owning checkpoint payload used by transform-specific
// persistence adapters. NetCDF encoding remains adapter-owned; this record
// names and sizes only the families present in the resolved transform.
struct WVTransformStateCheckpoint {
  std::string transformIdentifier;
  std::vector<std::size_t> spatialDimensions;
  double t = 0.0;
  double t0 = 0.0;
  std::vector<WVCoefficientFamilyCheckpoint> coefficientFamilies;
  std::size_t persistentBytes() const noexcept;
};

WVKernelStatus captureTransformStateCheckpoint(
    const WVIntegrationStateLayout &layout,
    const WVIntegrationState &state,
    WVTransformStateCheckpoint &checkpoint);
WVKernelStatus restoreTransformStateCheckpoint(
    const WVTransformStateCheckpoint &checkpoint,
    const WVIntegrationStateLayout &layout,
    WVCoefficientStateStorage &storage,
    WVMutableIntegrationState &state);

// Owns only additional state; canonical Ap, Am, and A0 remain caller-owned.
class WVAdditionalStateStorage final {
public:
  WVKernelStatus initialize(const WVIntegrationStateLayout &layout);
  WVAdditionalStateBlockView *mutableBlocks() noexcept {
    return mutableViews_.data();
  }
  const WVAdditionalStateBlockConstView *constBlocks() const noexcept {
    return constViews_.data();
  }
  std::size_t blockCount() const noexcept { return mutableViews_.size(); }
  std::size_t capacityBytes() const noexcept;
  void clear() noexcept;

private:
  std::vector<double> realStorage_;
  std::vector<WVComplex64> complexStorage_;
  std::vector<WVAdditionalStateBlockView> mutableViews_;
  std::vector<WVAdditionalStateBlockConstView> constViews_;
};

WVIntegrationState
integrationConstView(const WVMutableIntegrationState &state,
                   std::vector<WVAdditionalStateBlockConstView> &blockViews);
WVIntegrationState integrationConstView(
    const WVMutableIntegrationState &state,
    std::vector<WVCoefficientFamilyConstView> &coefficientViews,
    std::vector<WVAdditionalStateBlockConstView> &blockViews);
WVCoefficientFamilyConstView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    const WVIntegrationState &state, std::size_t family) noexcept;
WVCoefficientFamilyView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    WVMutableIntegrationState &state, std::size_t family) noexcept;
WVCoefficientFamilyView coefficientFamilyView(
    const WVIntegrationStateLayout &layout,
    WVIntegrationFlux &flux, std::size_t family) noexcept;
WVKernelStatus validateIntegrationState(const WVIntegrationStateLayout &layout,
                                      const WVIntegrationState &state);
WVKernelStatus
validateMutableIntegrationState(const WVIntegrationStateLayout &layout,
                              const WVMutableIntegrationState &state);
bool sameIntegrationStateLayout(const WVIntegrationStateLayout &first,
                                const WVIntegrationStateLayout &second) noexcept;

} // namespace wavevortex::runtime
