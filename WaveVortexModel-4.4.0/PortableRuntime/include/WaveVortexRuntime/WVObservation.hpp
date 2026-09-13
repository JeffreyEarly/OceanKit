#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace wavevortex::runtime {

inline constexpr std::uint32_t WVObservationSchemaContractVersion = 1;
inline constexpr std::size_t WVNoResolvedObservationVariable =
    static_cast<std::size_t>(-1);

enum class WVObservationScalarType : std::uint8_t {
  real64,
  complex64,
  integer64,
  boolean8,
  text
};

enum class WVObservationAxisKind : std::uint8_t { fixed, unlimited };

enum class WVObservationCoordinateRole : std::uint8_t {
  none,
  recordTime,
  sampleTime,
  x,
  y,
  z,
  identifier,
  depth,
  pass,
  profile
};

enum class WVObservationValueLayout : std::uint8_t {
  staticValue,
  initialValue,
  record,
  flat
};

enum class WVObservationRaggedRole : std::uint8_t {
  none,
  rowCount,
  rowOffset
};

enum class WVObservationBufferOwnership : std::uint8_t { borrowed, owned };
enum class WVObservationBatchKind : std::uint8_t { initial, event };

struct WVObservationAttribute {
  std::string name;
  std::string value;
};

struct WVObservationStringListAttribute {
  std::string name;
  std::vector<std::string> values;
};

struct WVObservationAxis {
  std::string identifier;
  std::string name;
  WVObservationAxisKind kind = WVObservationAxisKind::fixed;
  // Fixed axes have a positive extent. Unlimited axes use zero here and take
  // their occurrence extent from each batch.
  std::size_t extent = 0;
  WVObservationCoordinateRole coordinateRole =
      WVObservationCoordinateRole::none;
};

struct WVObservationVariable {
  std::string identifier;
  std::string name;
  WVObservationScalarType scalarType = WVObservationScalarType::real64;
  // Dimension identifiers use MATLAB logical order. Persistence adapters may
  // reverse them for NetCDF without changing the schema contract.
  std::vector<std::string> dimensionIdentifiers;
  WVObservationValueLayout layout = WVObservationValueLayout::record;
  std::string units;
  std::string description;
  std::vector<WVObservationAttribute> attributes;
  WVObservationCoordinateRole coordinateRole =
      WVObservationCoordinateRole::none;
  WVObservationRaggedRole raggedRole = WVObservationRaggedRole::none;
  // A row-count or row-offset variable is a rank-one integer over one
  // unlimited parent axis. It relates those parent rows to this child
  // unlimited axis. Each child has one relationship parent and the complete
  // relationship graph is acyclic.
  std::string raggedChildAxisIdentifier;
};

// One explicitly owned or borrowed contiguous typed value. Owned storage is
// held only in the matching owned* member. Borrowed storage is held only in
// the matching borrowed* pointer and must outlive the batch consumer.
struct WVObservationValue {
  std::string variableIdentifier;
  WVObservationScalarType scalarType = WVObservationScalarType::real64;
  WVObservationBufferOwnership ownership =
      WVObservationBufferOwnership::borrowed;
  std::vector<std::size_t> extents;
  // Construction-resolved schema ordinal used by the event persistence path.
  // Identifiers remain serialization metadata and are not searched at runtime.
  std::size_t resolvedVariableIndex = WVNoResolvedObservationVariable;

  const double *borrowedReal64 = nullptr;
  const WVComplex64 *borrowedComplex64 = nullptr;
  const std::int64_t *borrowedInteger64 = nullptr;
  const std::uint8_t *borrowedBoolean8 = nullptr;
  const std::string *borrowedText = nullptr;

  std::vector<double> ownedReal64;
  std::vector<WVComplex64> ownedComplex64;
  std::vector<std::int64_t> ownedInteger64;
  std::vector<std::uint8_t> ownedBoolean8;
  std::vector<std::string> ownedText;

  static WVObservationValue
  borrowReal(std::string identifier, std::vector<std::size_t> extents,
             const double *values);
  static WVObservationValue
  borrowComplex(std::string identifier, std::vector<std::size_t> extents,
                const WVComplex64 *values);
  static WVObservationValue
  borrowInteger(std::string identifier, std::vector<std::size_t> extents,
                const std::int64_t *values);
  static WVObservationValue
  borrowBoolean(std::string identifier, std::vector<std::size_t> extents,
                const std::uint8_t *values);
  static WVObservationValue
  borrowText(std::string identifier, std::vector<std::size_t> extents,
             const std::string *values);

  static WVObservationValue ownReal(std::string identifier,
                                    std::vector<std::size_t> extents,
                                    std::vector<double> values);
  static WVObservationValue ownComplex(std::string identifier,
                                       std::vector<std::size_t> extents,
                                       std::vector<WVComplex64> values);
  static WVObservationValue ownInteger(std::string identifier,
                                       std::vector<std::size_t> extents,
                                       std::vector<std::int64_t> values);
  static WVObservationValue ownBoolean(std::string identifier,
                                       std::vector<std::size_t> extents,
                                       std::vector<std::uint8_t> values);
  static WVObservationValue ownText(std::string identifier,
                                    std::vector<std::size_t> extents,
                                    std::vector<std::string> values);

  std::size_t elementCount() const noexcept;
  std::size_t liveBytes() const noexcept;
  std::size_t retainedBytes() const noexcept;
  const double *real64Data() const noexcept;
  const WVComplex64 *complex64Data() const noexcept;
  const std::int64_t *integer64Data() const noexcept;
  const std::uint8_t *boolean8Data() const noexcept;
  const std::string *textData() const noexcept;
};

struct WVObservationMetadataVariable {
  std::string name;
  WVObservationValue value;
  bool isLogicalType = false;
};

// Data-only legacy metadata declarations. The generic persistence adapter
// owns all NetCDF calls; observing systems only populate these values.
struct WVObservationMetadata {
  std::vector<WVObservationAttribute> attributes;
  std::vector<WVObservationStringListAttribute> stringListAttributes;
  std::vector<WVObservationMetadataVariable> variables;
};

// Data-only persistence contract with its own exact identifier and version.
// Its version is independent of the C++ source API and paired implementation
// versions. Implementations declare schemas but never own a NetCDF writer.
struct WVObservationSchema {
  std::string identifier;
  std::uint32_t version = WVObservationSchemaContractVersion;
  // Legacy MATLAB-compatible schemas suppress portable schema attributes so
  // existing fixed output remains byte/schema compatible.
  bool preservesLegacyEncoding = false;
  WVObservationMetadata metadata;
  std::vector<WVObservationAxis> axes;
  std::vector<WVObservationVariable> variables;
};

struct WVObservationBatchMetrics {
  std::size_t liveBytes = 0;
  std::size_t retainedStorageBytes = 0;
};

// One initial or event value set for an exact schema identifier/version.
// Borrowed buffers must outlive batch consumption; owned buffers move with the
// batch. Persistence validates the complete batch before output mutation.
struct WVObservationBatch {
  WVObservationBatch() = default;
  WVObservationBatch(WVObservationBatch &&) noexcept = default;
  WVObservationBatch &operator=(WVObservationBatch &&) noexcept = default;
  WVObservationBatch(const WVObservationBatch &) = delete;
  WVObservationBatch &operator=(const WVObservationBatch &) = delete;

  std::string schemaIdentifier;
  std::uint32_t schemaVersion = WVObservationSchemaContractVersion;
  WVObservationBatchKind kind = WVObservationBatchKind::event;
  std::vector<WVObservationValue> values;

  WVObservationBatchMetrics metrics() const noexcept;
};

WVKernelStatus validateObservationSchema(const WVObservationSchema &schema);
WVKernelStatus validateObservationBatch(const WVObservationSchema &schema,
                                        const WVObservationBatch &batch);
std::size_t
observationSchemaRetainedBytes(const WVObservationSchema &schema) noexcept;
WVKernelStatus encodeObservationSchemaManifest(
    const WVObservationSchema &schema, std::vector<std::uint8_t> &bytes);
WVKernelStatus decodeObservationSchemaManifest(
    const std::vector<std::uint8_t> &bytes, WVObservationSchema &schema);

} // namespace wavevortex::runtime
