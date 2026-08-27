#include "WaveVortexRuntime/WVObservation.hpp"
#include "WaveVortexRuntime/WVPortableTypedRecord.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

bool validIdentifier(const std::string &value) {
  if (value.empty())
    return false;
  return std::all_of(value.begin(), value.end(), [](unsigned char character) {
    return (character >= 'a' && character <= 'z') ||
           (character >= 'A' && character <= 'Z') ||
           (character >= '0' && character <= '9') || character == '-' ||
           character == '_' || character == '.';
  });
}

bool multiply(std::size_t left, std::size_t right,
              std::size_t &product) noexcept {
  if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right)
    return false;
  product = left * right;
  return true;
}

std::size_t scalarBytes(WVObservationScalarType type) noexcept {
  switch (type) {
  case WVObservationScalarType::real64:
    return sizeof(double);
  case WVObservationScalarType::complex64:
    return sizeof(WVComplex64);
  case WVObservationScalarType::integer64:
    return sizeof(std::int64_t);
  case WVObservationScalarType::boolean8:
    return sizeof(std::uint8_t);
  case WVObservationScalarType::text:
    return sizeof(std::string);
  }
  return 0;
}

template <typename T>
const T *storageData(WVObservationBufferOwnership ownership,
                     const T *borrowed,
                     const std::vector<T> &owned) noexcept {
  return ownership == WVObservationBufferOwnership::borrowed
             ? borrowed
             : owned.data();
}

WVObservationValue borrowedValue(std::string identifier,
                                 WVObservationScalarType type,
                                 std::vector<std::size_t> extents) {
  WVObservationValue value;
  value.variableIdentifier = std::move(identifier);
  value.scalarType = type;
  value.ownership = WVObservationBufferOwnership::borrowed;
  value.extents = std::move(extents);
  return value;
}

WVObservationValue ownedValue(std::string identifier,
                              WVObservationScalarType type,
                              std::vector<std::size_t> extents) {
  WVObservationValue value;
  value.variableIdentifier = std::move(identifier);
  value.scalarType = type;
  value.ownership = WVObservationBufferOwnership::owned;
  value.extents = std::move(extents);
  return value;
}

bool applicable(const WVObservationVariable &variable,
                WVObservationBatchKind kind) noexcept {
  const bool initial = variable.layout == WVObservationValueLayout::staticValue ||
                       variable.layout == WVObservationValueLayout::initialValue;
  return kind == WVObservationBatchKind::initial ? initial : !initial;
}

} // namespace

WVObservationValue
WVObservationValue::borrowReal(std::string identifier,
                               std::vector<std::size_t> extents,
                               const double *values) {
  auto value = borrowedValue(std::move(identifier),
                             WVObservationScalarType::real64,
                             std::move(extents));
  value.borrowedReal64 = values;
  return value;
}

WVObservationValue
WVObservationValue::borrowComplex(std::string identifier,
                                  std::vector<std::size_t> extents,
                                  const WVComplex64 *values) {
  auto value = borrowedValue(std::move(identifier),
                             WVObservationScalarType::complex64,
                             std::move(extents));
  value.borrowedComplex64 = values;
  return value;
}

WVObservationValue
WVObservationValue::borrowInteger(std::string identifier,
                                  std::vector<std::size_t> extents,
                                  const std::int64_t *values) {
  auto value = borrowedValue(std::move(identifier),
                             WVObservationScalarType::integer64,
                             std::move(extents));
  value.borrowedInteger64 = values;
  return value;
}

WVObservationValue
WVObservationValue::borrowBoolean(std::string identifier,
                                  std::vector<std::size_t> extents,
                                  const std::uint8_t *values) {
  auto value = borrowedValue(std::move(identifier),
                             WVObservationScalarType::boolean8,
                             std::move(extents));
  value.borrowedBoolean8 = values;
  return value;
}

WVObservationValue
WVObservationValue::borrowText(std::string identifier,
                               std::vector<std::size_t> extents,
                               const std::string *values) {
  auto value = borrowedValue(std::move(identifier), WVObservationScalarType::text,
                             std::move(extents));
  value.borrowedText = values;
  return value;
}

WVObservationValue
WVObservationValue::ownReal(std::string identifier,
                            std::vector<std::size_t> extents,
                            std::vector<double> values) {
  auto value = ownedValue(std::move(identifier), WVObservationScalarType::real64,
                          std::move(extents));
  value.ownedReal64 = std::move(values);
  return value;
}

WVObservationValue
WVObservationValue::ownComplex(std::string identifier,
                               std::vector<std::size_t> extents,
                               std::vector<WVComplex64> values) {
  auto value = ownedValue(std::move(identifier),
                          WVObservationScalarType::complex64,
                          std::move(extents));
  value.ownedComplex64 = std::move(values);
  return value;
}

WVObservationValue
WVObservationValue::ownInteger(std::string identifier,
                               std::vector<std::size_t> extents,
                               std::vector<std::int64_t> values) {
  auto value = ownedValue(std::move(identifier),
                          WVObservationScalarType::integer64,
                          std::move(extents));
  value.ownedInteger64 = std::move(values);
  return value;
}

WVObservationValue
WVObservationValue::ownBoolean(std::string identifier,
                               std::vector<std::size_t> extents,
                               std::vector<std::uint8_t> values) {
  auto value = ownedValue(std::move(identifier),
                          WVObservationScalarType::boolean8,
                          std::move(extents));
  value.ownedBoolean8 = std::move(values);
  return value;
}

WVObservationValue
WVObservationValue::ownText(std::string identifier,
                            std::vector<std::size_t> extents,
                            std::vector<std::string> values) {
  auto value = ownedValue(std::move(identifier), WVObservationScalarType::text,
                          std::move(extents));
  value.ownedText = std::move(values);
  return value;
}

std::size_t WVObservationValue::elementCount() const noexcept {
  std::size_t count = 1;
  for (const auto extent : extents)
    if (!multiply(count, extent, count))
      return std::numeric_limits<std::size_t>::max();
  return count;
}

std::size_t WVObservationValue::liveBytes() const noexcept {
  const auto count = elementCount();
  if (count == std::numeric_limits<std::size_t>::max())
    return count;
  std::size_t bytes = 0;
  if (!multiply(count, scalarBytes(scalarType), bytes))
    return std::numeric_limits<std::size_t>::max();
  if (scalarType == WVObservationScalarType::text) {
    const auto *values = textData();
    for (std::size_t index = 0; values != nullptr && index < count; ++index) {
      if (bytes > std::numeric_limits<std::size_t>::max() -
                      values[index].size())
        return std::numeric_limits<std::size_t>::max();
      bytes += values[index].size();
    }
  }
  return bytes;
}

std::size_t WVObservationValue::retainedBytes() const noexcept {
  if (ownership == WVObservationBufferOwnership::borrowed)
    return variableIdentifier.capacity() +
           extents.capacity() * sizeof(std::size_t);
  std::size_t bytes = variableIdentifier.capacity() +
                      extents.capacity() * sizeof(std::size_t) +
                      ownedReal64.capacity() * sizeof(double) +
                      ownedComplex64.capacity() * sizeof(WVComplex64) +
                      ownedInteger64.capacity() * sizeof(std::int64_t) +
                      ownedBoolean8.capacity() * sizeof(std::uint8_t) +
                      ownedText.capacity() * sizeof(std::string);
  for (const auto &value : ownedText)
    bytes += value.capacity();
  return bytes;
}

const double *WVObservationValue::real64Data() const noexcept {
  return scalarType == WVObservationScalarType::real64
             ? storageData(ownership, borrowedReal64, ownedReal64)
             : nullptr;
}

const WVComplex64 *WVObservationValue::complex64Data() const noexcept {
  return scalarType == WVObservationScalarType::complex64
             ? storageData(ownership, borrowedComplex64, ownedComplex64)
             : nullptr;
}

const std::int64_t *WVObservationValue::integer64Data() const noexcept {
  return scalarType == WVObservationScalarType::integer64
             ? storageData(ownership, borrowedInteger64, ownedInteger64)
             : nullptr;
}

const std::uint8_t *WVObservationValue::boolean8Data() const noexcept {
  return scalarType == WVObservationScalarType::boolean8
             ? storageData(ownership, borrowedBoolean8, ownedBoolean8)
             : nullptr;
}

const std::string *WVObservationValue::textData() const noexcept {
  return scalarType == WVObservationScalarType::text
             ? storageData(ownership, borrowedText, ownedText)
             : nullptr;
}

std::size_t
observationSchemaRetainedBytes(const WVObservationSchema &schema) noexcept {
  std::size_t bytes = schema.identifier.capacity() +
                      schema.axes.capacity() * sizeof(WVObservationAxis) +
                      schema.variables.capacity() *
                          sizeof(WVObservationVariable) +
                      schema.metadata.attributes.capacity() *
                          sizeof(WVObservationAttribute) +
                      schema.metadata.stringListAttributes.capacity() *
                          sizeof(WVObservationStringListAttribute) +
                      schema.metadata.variables.capacity() *
                          sizeof(WVObservationMetadataVariable);
  for (const auto &axis : schema.axes)
    bytes += axis.identifier.capacity() + axis.name.capacity();
  for (const auto &variable : schema.variables) {
    bytes += variable.identifier.capacity() + variable.name.capacity() +
             variable.units.capacity() + variable.description.capacity() +
             variable.raggedChildAxisIdentifier.capacity() +
             variable.dimensionIdentifiers.capacity() * sizeof(std::string) +
             variable.attributes.capacity() * sizeof(WVObservationAttribute);
    for (const auto &identifier : variable.dimensionIdentifiers)
      bytes += identifier.capacity();
    for (const auto &attribute : variable.attributes)
      bytes += attribute.name.capacity() + attribute.value.capacity();
  }
  for (const auto &attribute : schema.metadata.attributes)
    bytes += attribute.name.capacity() + attribute.value.capacity();
  for (const auto &attribute : schema.metadata.stringListAttributes) {
    bytes += attribute.name.capacity() +
             attribute.values.capacity() * sizeof(std::string);
    for (const auto &value : attribute.values)
      bytes += value.capacity();
  }
  for (const auto &variable : schema.metadata.variables)
    bytes += variable.name.capacity() + variable.value.retainedBytes();
  return bytes;
}

WVObservationBatchMetrics WVObservationBatch::metrics() const noexcept {
  WVObservationBatchMetrics result;
  for (const auto &value : values) {
    const auto live = value.liveBytes();
    const auto retained = value.retainedBytes();
    if (result.liveBytes > std::numeric_limits<std::size_t>::max() - live)
      result.liveBytes = std::numeric_limits<std::size_t>::max();
    else
      result.liveBytes += live;
    if (result.retainedStorageBytes >
        std::numeric_limits<std::size_t>::max() - retained)
      result.retainedStorageBytes = std::numeric_limits<std::size_t>::max();
    else
      result.retainedStorageBytes += retained;
  }
  result.retainedStorageBytes +=
      schemaIdentifier.capacity() + values.capacity() * sizeof(WVObservationValue);
  return result;
}

WVKernelStatus validateObservationSchema(const WVObservationSchema &schema) {
  if (!validIdentifier(schema.identifier) || schema.version == 0)
    return invalid("Observation schema identity and version must be valid.");
  std::map<std::string, const WVObservationAxis *> axes;
  std::set<std::string> axisNames;
  for (const auto &axis : schema.axes) {
    if (!validIdentifier(axis.identifier) || axis.name.empty() ||
        !axisNames.insert(axis.name).second ||
        !axes.emplace(axis.identifier, &axis).second)
      return invalid("Observation axes require unique identities and names.");
    if ((axis.kind == WVObservationAxisKind::fixed && axis.extent == 0) ||
        (axis.kind == WVObservationAxisKind::unlimited && axis.extent != 0))
      return invalid("Observation axis extent is inconsistent with its kind.");
  }
  std::set<std::string> variableIdentifiers;
  std::set<std::string> variableNames;
  std::map<std::string, std::string> raggedParentByChild;
  for (const auto &variable : schema.variables) {
    if (!validIdentifier(variable.identifier) || variable.name.empty() ||
        !variableIdentifiers.insert(variable.identifier).second ||
        !variableNames.insert(variable.name).second)
      return invalid(
          "Observation variables require unique identities and names.");
    std::set<std::string> dimensions;
    bool hasUnlimitedAxis = false;
    for (const auto &identifier : variable.dimensionIdentifiers) {
      const auto found = axes.find(identifier);
      if (found == axes.end() || !dimensions.insert(identifier).second)
        return invalid("Observation variable dimensions are undeclared or repeated.");
      hasUnlimitedAxis =
          hasUnlimitedAxis ||
          found->second->kind == WVObservationAxisKind::unlimited;
    }
    if (variable.layout == WVObservationValueLayout::flat &&
        !hasUnlimitedAxis)
      return invalid("Flat observation variables require an unlimited axis.");
    if (variable.layout != WVObservationValueLayout::flat && hasUnlimitedAxis)
      return invalid(
          "Only flat observation variables may use unlimited axes.");
    if (variable.raggedRole != WVObservationRaggedRole::none) {
      const auto child = axes.find(variable.raggedChildAxisIdentifier);
      if (variable.scalarType != WVObservationScalarType::integer64 ||
          child == axes.end() ||
          child->second->kind != WVObservationAxisKind::unlimited)
        return invalid(
            "Ragged relationships require integer values and a declared child unlimited axis.");
      if (variable.dimensionIdentifiers.size() != 1)
        return invalid(
            "Ragged relationship variables require exactly one parent axis.");
      const auto parent = axes.find(variable.dimensionIdentifiers.front());
      if (parent == axes.end() ||
          parent->second->kind != WVObservationAxisKind::unlimited)
        return invalid(
            "Ragged relationship parent axes must be declared and unlimited.");
      if (parent->first == child->first)
        return invalid("Ragged relationships cannot contain self edges.");
      if (!raggedParentByChild.emplace(child->first, parent->first).second)
        return invalid(
            "A ragged child axis cannot have multiple relationship parents.");
    } else if (!variable.raggedChildAxisIdentifier.empty()) {
      return invalid("Non-ragged variables cannot name a ragged child axis.");
    }
    std::set<std::string> attributeNames;
    for (const auto &attribute : variable.attributes)
      if (attribute.name.empty() ||
          !attributeNames.insert(attribute.name).second)
        return invalid("Observation variable attributes must be unique.");
  }
  for (const auto &relationship : raggedParentByChild) {
    std::set<std::string> visited;
    auto axisIdentifier = relationship.first;
    while (true) {
      if (!visited.insert(axisIdentifier).second)
        return invalid("Ragged relationships must form an acyclic graph.");
      const auto parent = raggedParentByChild.find(axisIdentifier);
      if (parent == raggedParentByChild.end())
        break;
      axisIdentifier = parent->second;
    }
  }
  return WVKernelStatus::ok();
}

WVKernelStatus validateObservationBatch(const WVObservationSchema &schema,
                                        const WVObservationBatch &batch) {
  const auto schemaStatus = validateObservationSchema(schema);
  if (!schemaStatus)
    return schemaStatus;
  if (batch.schemaIdentifier != schema.identifier ||
      batch.schemaVersion != schema.version)
    return invalid("Observation batch schema identity or version drifted.");
  std::map<std::string, const WVObservationAxis *> axes;
  for (const auto &axis : schema.axes)
    axes.emplace(axis.identifier, &axis);
  std::map<std::string, const WVObservationVariable *> variables;
  std::size_t expectedValueCount = 0;
  for (const auto &variable : schema.variables) {
    variables.emplace(variable.identifier, &variable);
    expectedValueCount += applicable(variable, batch.kind) ? 1 : 0;
  }
  if (batch.values.size() != expectedValueCount)
    return invalid("Observation batch does not contain every applicable declared variable.");

  std::set<std::string> observedVariables;
  std::map<std::string, std::size_t> unlimitedExtents;
  std::map<std::string, const WVObservationValue *> valuesByIdentifier;
  for (const auto &value : batch.values) {
    const auto found = variables.find(value.variableIdentifier);
    if (found == variables.end())
      return invalid("Observation batch contains an undeclared variable.");
    const auto &variable = *found->second;
    const auto variableIndex = static_cast<std::size_t>(
        found->second - schema.variables.data());
    if (value.resolvedVariableIndex != WVNoResolvedObservationVariable &&
        value.resolvedVariableIndex != variableIndex)
      return invalid("Observation batch resolved-variable slot drifted.");
    if (!applicable(variable, batch.kind) ||
        !observedVariables.insert(value.variableIdentifier).second)
      return invalid("Observation batch variable is duplicated or belongs to another phase.");
    if (value.scalarType != variable.scalarType ||
        value.extents.size() != variable.dimensionIdentifiers.size())
      return invalid("Observation batch value type or rank differs from its schema.");
    std::size_t expectedElements = 1;
    for (std::size_t index = 0; index < value.extents.size(); ++index) {
      const auto &axis = *axes.at(variable.dimensionIdentifiers[index]);
      const auto extent = value.extents[index];
      if (axis.kind == WVObservationAxisKind::fixed) {
        if (extent != axis.extent)
          return invalid("Observation batch fixed-axis extent changed.");
      } else {
        const auto prior = unlimitedExtents.find(axis.identifier);
        if (prior == unlimitedExtents.end())
          unlimitedExtents.emplace(axis.identifier, extent);
        else if (prior->second != extent)
          return invalid("Observation batch unlimited-axis extents are inconsistent.");
      }
      if (!multiply(expectedElements, extent, expectedElements))
        return {WVKernelStatusCode::sizeOverflow,
                "Observation batch element count overflows size_t."};
    }
    if (value.elementCount() != expectedElements)
      return {WVKernelStatusCode::sizeOverflow,
              "Observation batch element count overflows size_t."};
    const bool hasData =
        expectedElements == 0 ||
        (value.scalarType == WVObservationScalarType::real64 &&
         value.real64Data() != nullptr) ||
        (value.scalarType == WVObservationScalarType::complex64 &&
         value.complex64Data() != nullptr) ||
        (value.scalarType == WVObservationScalarType::integer64 &&
         value.integer64Data() != nullptr) ||
        (value.scalarType == WVObservationScalarType::boolean8 &&
         value.boolean8Data() != nullptr) ||
        (value.scalarType == WVObservationScalarType::text &&
         value.textData() != nullptr);
    if (!hasData)
      return invalid("Observation batch value has no compatible contiguous buffer.");
    if (value.ownership == WVObservationBufferOwnership::owned) {
      const std::size_t ownedCount =
          value.scalarType == WVObservationScalarType::real64
              ? value.ownedReal64.size()
              : value.scalarType == WVObservationScalarType::complex64
                    ? value.ownedComplex64.size()
                    : value.scalarType == WVObservationScalarType::integer64
                          ? value.ownedInteger64.size()
                          : value.scalarType == WVObservationScalarType::boolean8
                                ? value.ownedBoolean8.size()
                                : value.ownedText.size();
      if (ownedCount != expectedElements)
        return invalid("Owned observation buffer size differs from its extents.");
    }
    if (value.scalarType == WVObservationScalarType::boolean8) {
      const auto *booleans = value.boolean8Data();
      for (std::size_t index = 0; index < expectedElements; ++index)
        if (booleans[index] > 1)
          return invalid("Boolean observation values must contain zero or one.");
    }
    valuesByIdentifier.emplace(value.variableIdentifier, &value);
  }

  for (const auto &variable : schema.variables) {
    if (!applicable(variable, batch.kind) ||
        variable.raggedRole == WVObservationRaggedRole::none)
      continue;
    const auto &value = *valuesByIdentifier.at(variable.identifier);
    const auto childExtent = unlimitedExtents.find(
        variable.raggedChildAxisIdentifier);
    if (childExtent == unlimitedExtents.end())
      return invalid("Ragged child axis is absent from the observation batch.");
    const auto *integers = value.integer64Data();
    const auto count = value.elementCount();
    if (variable.raggedRole == WVObservationRaggedRole::rowCount) {
      std::size_t total = 0;
      for (std::size_t index = 0; index < count; ++index) {
        if (integers[index] < 0 ||
            static_cast<std::uint64_t>(integers[index]) >
                std::numeric_limits<std::size_t>::max())
          return invalid("Ragged row counts must be nonnegative and bounded.");
        const auto row = static_cast<std::size_t>(integers[index]);
        if (total > std::numeric_limits<std::size_t>::max() - row)
          return {WVKernelStatusCode::sizeOverflow,
                  "Ragged row-count sum overflows size_t."};
        total += row;
      }
      if (total != childExtent->second)
        return invalid("Ragged row counts do not span the child axis.");
    } else {
      if (count == 0) {
        if (childExtent->second != 0)
          return invalid(
              "An empty ragged parent cannot reference a nonempty child axis.");
        continue;
      }
      if (integers[0] != 0)
        return invalid("Ragged row offsets must begin at zero in each batch.");
      for (std::size_t index = 0; index < count; ++index)
        if (integers[index] < 0 ||
            static_cast<std::uint64_t>(integers[index]) > childExtent->second ||
            (index > 0 && integers[index] < integers[index - 1]))
          return invalid("Ragged row offsets are malformed.");
    }
  }
  return WVKernelStatus::ok();
}

WVKernelStatus encodeObservationSchemaManifest(
    const WVObservationSchema &schema, std::vector<std::uint8_t> &bytes) {
  const auto schemaStatus = validateObservationSchema(schema);
  if (!schemaStatus)
    return schemaStatus;
  WVPortableTypedRecord manifest;
  manifest.schemaIdentifier = "portable-observation-schema-manifest-v1";
  manifest.schemaVersion = 1;
  const auto addText = [&](std::string name,
                           std::vector<std::string> values) {
    const auto count = values.size();
    manifest.values.push_back(
        {std::move(name), {count}, std::move(values)});
  };
  const auto addInteger = [&](std::string name,
                              std::vector<std::int64_t> values) {
    const auto count = values.size();
    manifest.values.push_back(
        {std::move(name), {count}, std::move(values)});
  };
  manifest.values.push_back(
      {"schemaIdentifier", {}, std::vector<std::string>{schema.identifier}});
  manifest.values.push_back(
      {"schemaVersion", {},
       std::vector<std::int64_t>{static_cast<std::int64_t>(schema.version)}});
  manifest.values.push_back(
      {"preservesLegacyEncoding", {},
       std::vector<std::uint8_t>{
           static_cast<std::uint8_t>(schema.preservesLegacyEncoding)}});

  std::vector<std::string> metadataAttributeNames;
  std::vector<std::string> metadataAttributeValues;
  for (const auto &attribute : schema.metadata.attributes) {
    metadataAttributeNames.push_back(attribute.name);
    metadataAttributeValues.push_back(attribute.value);
  }
  addText("metadataAttributeNames", std::move(metadataAttributeNames));
  addText("metadataAttributeValues", std::move(metadataAttributeValues));

  std::vector<std::string> metadataStringListNames;
  std::vector<std::int64_t> metadataStringListCounts;
  std::vector<std::string> metadataStringListValues;
  for (const auto &attribute : schema.metadata.stringListAttributes) {
    metadataStringListNames.push_back(attribute.name);
    metadataStringListCounts.push_back(
        static_cast<std::int64_t>(attribute.values.size()));
    metadataStringListValues.insert(metadataStringListValues.end(),
                                    attribute.values.begin(),
                                    attribute.values.end());
  }
  addText("metadataStringListNames", std::move(metadataStringListNames));
  addInteger("metadataStringListCounts",
             std::move(metadataStringListCounts));
  addText("metadataStringListValues", std::move(metadataStringListValues));

  std::vector<std::string> metadataVariableNames;
  std::vector<std::string> metadataVariableIdentifiers;
  std::vector<std::int64_t> metadataVariableTypes;
  std::vector<std::uint8_t> metadataVariableLogical;
  std::vector<std::int64_t> metadataVariableDimensionCounts;
  std::vector<std::int64_t> metadataVariableExtents;
  for (std::size_t index = 0; index < schema.metadata.variables.size();
       ++index) {
    const auto &variable = schema.metadata.variables[index];
    const auto &value = variable.value;
    const auto elementCount = value.elementCount();
    if (elementCount == std::numeric_limits<std::size_t>::max())
      return {WVKernelStatusCode::sizeOverflow,
              "Observation metadata-variable extent overflows size_t."};
    metadataVariableNames.push_back(variable.name);
    metadataVariableIdentifiers.push_back(value.variableIdentifier);
    metadataVariableTypes.push_back(
        static_cast<std::int64_t>(value.scalarType));
    metadataVariableLogical.push_back(
        static_cast<std::uint8_t>(variable.isLogicalType));
    metadataVariableDimensionCounts.push_back(
        static_cast<std::int64_t>(value.extents.size()));
    for (const auto extent : value.extents) {
      if (extent >
          static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        return {WVKernelStatusCode::sizeOverflow,
                "Observation metadata-variable extent exceeds int64."};
      metadataVariableExtents.push_back(static_cast<std::int64_t>(extent));
    }
    const auto storageName =
        "metadataVariable." + std::to_string(index);
    switch (value.scalarType) {
    case WVObservationScalarType::real64: {
      const auto *data = value.real64Data();
      if (elementCount > 0 && data == nullptr)
        return invalid("Observation metadata variable has no real storage.");
      std::vector<double> values;
      if (elementCount > 0)
        values.assign(data, data + elementCount);
      manifest.values.push_back(
          {storageName, value.extents, std::move(values)});
      break;
    }
    case WVObservationScalarType::complex64: {
      const auto *data = value.complex64Data();
      if (elementCount > 0 && data == nullptr)
        return invalid(
            "Observation metadata variable has no complex storage.");
      std::vector<double> real(elementCount);
      std::vector<double> imaginary(elementCount);
      for (std::size_t element = 0; element < elementCount; ++element) {
        real[element] = data[element].real;
        imaginary[element] = data[element].imag;
      }
      manifest.values.push_back(
          {storageName + ".real", value.extents, std::move(real)});
      manifest.values.push_back(
          {storageName + ".imag", value.extents, std::move(imaginary)});
      break;
    }
    case WVObservationScalarType::integer64: {
      const auto *data = value.integer64Data();
      if (elementCount > 0 && data == nullptr)
        return invalid(
            "Observation metadata variable has no integer storage.");
      std::vector<std::int64_t> values;
      if (elementCount > 0)
        values.assign(data, data + elementCount);
      manifest.values.push_back(
          {storageName, value.extents, std::move(values)});
      break;
    }
    case WVObservationScalarType::boolean8: {
      const auto *data = value.boolean8Data();
      if (elementCount > 0 && data == nullptr)
        return invalid(
            "Observation metadata variable has no Boolean storage.");
      std::vector<std::uint8_t> values;
      if (elementCount > 0)
        values.assign(data, data + elementCount);
      manifest.values.push_back(
          {storageName, value.extents, std::move(values)});
      break;
    }
    case WVObservationScalarType::text: {
      const auto *data = value.textData();
      if (elementCount > 0 && data == nullptr)
        return invalid("Observation metadata variable has no text storage.");
      std::vector<std::string> values;
      if (elementCount > 0)
        values.assign(data, data + elementCount);
      manifest.values.push_back(
          {storageName, value.extents, std::move(values)});
      break;
    }
    }
  }
  addText("metadataVariableNames", std::move(metadataVariableNames));
  addText("metadataVariableIdentifiers",
          std::move(metadataVariableIdentifiers));
  addInteger("metadataVariableTypes", std::move(metadataVariableTypes));
  manifest.values.push_back(
      {"metadataVariableLogical", {metadataVariableLogical.size()},
       std::move(metadataVariableLogical)});
  addInteger("metadataVariableDimensionCounts",
             std::move(metadataVariableDimensionCounts));
  addInteger("metadataVariableExtents", std::move(metadataVariableExtents));

  std::vector<std::string> axisIdentifiers;
  std::vector<std::string> axisNames;
  std::vector<std::int64_t> axisKinds;
  std::vector<std::int64_t> axisExtents;
  std::vector<std::int64_t> axisRoles;
  for (const auto &axis : schema.axes) {
    if (axis.extent >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
      return {WVKernelStatusCode::sizeOverflow,
              "Observation schema axis extent exceeds int64."};
    axisIdentifiers.push_back(axis.identifier);
    axisNames.push_back(axis.name);
    axisKinds.push_back(static_cast<std::int64_t>(axis.kind));
    axisExtents.push_back(static_cast<std::int64_t>(axis.extent));
    axisRoles.push_back(static_cast<std::int64_t>(axis.coordinateRole));
  }
  addText("axisIdentifiers", std::move(axisIdentifiers));
  addText("axisNames", std::move(axisNames));
  addInteger("axisKinds", std::move(axisKinds));
  addInteger("axisExtents", std::move(axisExtents));
  addInteger("axisRoles", std::move(axisRoles));

  std::vector<std::string> variableIdentifiers;
  std::vector<std::string> variableNames;
  std::vector<std::int64_t> variableTypes;
  std::vector<std::int64_t> variableLayouts;
  std::vector<std::string> variableUnits;
  std::vector<std::string> variableDescriptions;
  std::vector<std::int64_t> variableRoles;
  std::vector<std::int64_t> variableRaggedRoles;
  std::vector<std::string> variableRaggedChildren;
  std::vector<std::int64_t> dimensionCounts;
  std::vector<std::string> dimensions;
  std::vector<std::int64_t> attributeCounts;
  std::vector<std::string> attributeNames;
  std::vector<std::string> attributeValues;
  for (const auto &variable : schema.variables) {
    variableIdentifiers.push_back(variable.identifier);
    variableNames.push_back(variable.name);
    variableTypes.push_back(static_cast<std::int64_t>(variable.scalarType));
    variableLayouts.push_back(static_cast<std::int64_t>(variable.layout));
    variableUnits.push_back(variable.units);
    variableDescriptions.push_back(variable.description);
    variableRoles.push_back(
        static_cast<std::int64_t>(variable.coordinateRole));
    variableRaggedRoles.push_back(
        static_cast<std::int64_t>(variable.raggedRole));
    variableRaggedChildren.push_back(variable.raggedChildAxisIdentifier);
    dimensionCounts.push_back(
        static_cast<std::int64_t>(variable.dimensionIdentifiers.size()));
    dimensions.insert(dimensions.end(), variable.dimensionIdentifiers.begin(),
                      variable.dimensionIdentifiers.end());
    attributeCounts.push_back(
        static_cast<std::int64_t>(variable.attributes.size()));
    for (const auto &attribute : variable.attributes) {
      attributeNames.push_back(attribute.name);
      attributeValues.push_back(attribute.value);
    }
  }
  addText("variableIdentifiers", std::move(variableIdentifiers));
  addText("variableNames", std::move(variableNames));
  addInteger("variableTypes", std::move(variableTypes));
  addInteger("variableLayouts", std::move(variableLayouts));
  addText("variableUnits", std::move(variableUnits));
  addText("variableDescriptions", std::move(variableDescriptions));
  addInteger("variableRoles", std::move(variableRoles));
  addInteger("variableRaggedRoles", std::move(variableRaggedRoles));
  addText("variableRaggedChildren", std::move(variableRaggedChildren));
  addInteger("dimensionCounts", std::move(dimensionCounts));
  addText("dimensions", std::move(dimensions));
  addInteger("attributeCounts", std::move(attributeCounts));
  addText("attributeNames", std::move(attributeNames));
  addText("attributeValues", std::move(attributeValues));
  return encodePortableTypedRecord(manifest, bytes);
}

WVKernelStatus decodeObservationSchemaManifest(
    const std::vector<std::uint8_t> &bytes, WVObservationSchema &schema) {
  WVPortableTypedRecord manifest;
  auto status = decodePortableTypedRecord(
      bytes, manifest, {4 * 1024 * 1024, false, true});
  if (!status)
    return status;
  if (manifest.schemaIdentifier !=
          "portable-observation-schema-manifest-v1" ||
      manifest.schemaVersion != 1)
    return invalid("Unsupported observation-schema manifest.");
  const auto texts = [&](const char *name)
      -> const std::vector<std::string> * {
    const auto *value = manifest.value(name);
    return value != nullptr &&
                   std::holds_alternative<std::vector<std::string>>(
                       value->storage)
               ? &std::get<std::vector<std::string>>(value->storage)
               : nullptr;
  };
  const auto integers = [&](const char *name)
      -> const std::vector<std::int64_t> * {
    const auto *value = manifest.value(name);
    return value != nullptr &&
                   std::holds_alternative<std::vector<std::int64_t>>(
                       value->storage)
               ? &std::get<std::vector<std::int64_t>>(value->storage)
               : nullptr;
  };
  const auto booleans = [&](const char *name)
      -> const std::vector<std::uint8_t> * {
    const auto *value = manifest.value(name);
    return value != nullptr &&
                   std::holds_alternative<std::vector<std::uint8_t>>(
                       value->storage)
               ? &std::get<std::vector<std::uint8_t>>(value->storage)
               : nullptr;
  };
  const auto *schemaIdentifiers = texts("schemaIdentifier");
  const auto *schemaVersions = integers("schemaVersion");
  const auto *preservesLegacyEncoding =
      booleans("preservesLegacyEncoding");
  const auto *metadataAttributeNames = texts("metadataAttributeNames");
  const auto *metadataAttributeValues = texts("metadataAttributeValues");
  const auto *metadataStringListNames = texts("metadataStringListNames");
  const auto *metadataStringListCounts =
      integers("metadataStringListCounts");
  const auto *metadataStringListValues = texts("metadataStringListValues");
  const auto *metadataVariableNames = texts("metadataVariableNames");
  const auto *metadataVariableIdentifiers =
      texts("metadataVariableIdentifiers");
  const auto *metadataVariableTypes = integers("metadataVariableTypes");
  const auto *metadataVariableLogical =
      booleans("metadataVariableLogical");
  const auto *metadataVariableDimensionCounts =
      integers("metadataVariableDimensionCounts");
  const auto *metadataVariableExtents =
      integers("metadataVariableExtents");
  const auto *axisIdentifiers = texts("axisIdentifiers");
  const auto *axisNames = texts("axisNames");
  const auto *axisKinds = integers("axisKinds");
  const auto *axisExtents = integers("axisExtents");
  const auto *axisRoles = integers("axisRoles");
  const auto *variableIdentifiers = texts("variableIdentifiers");
  const auto *variableNames = texts("variableNames");
  const auto *variableTypes = integers("variableTypes");
  const auto *variableLayouts = integers("variableLayouts");
  const auto *variableUnits = texts("variableUnits");
  const auto *variableDescriptions = texts("variableDescriptions");
  const auto *variableRoles = integers("variableRoles");
  const auto *variableRaggedRoles = integers("variableRaggedRoles");
  const auto *variableRaggedChildren = texts("variableRaggedChildren");
  const auto *dimensionCounts = integers("dimensionCounts");
  const auto *dimensions = texts("dimensions");
  const auto *attributeCounts = integers("attributeCounts");
  const auto *attributeNames = texts("attributeNames");
  const auto *attributeValues = texts("attributeValues");
  if (schemaIdentifiers == nullptr || schemaIdentifiers->size() != 1 ||
      schemaVersions == nullptr || schemaVersions->size() != 1 ||
      preservesLegacyEncoding == nullptr ||
      preservesLegacyEncoding->size() != 1 ||
      metadataAttributeNames == nullptr ||
      metadataAttributeValues == nullptr ||
      metadataStringListNames == nullptr ||
      metadataStringListCounts == nullptr ||
      metadataStringListValues == nullptr ||
      metadataVariableNames == nullptr ||
      metadataVariableIdentifiers == nullptr ||
      metadataVariableTypes == nullptr ||
      metadataVariableLogical == nullptr ||
      metadataVariableDimensionCounts == nullptr ||
      metadataVariableExtents == nullptr ||
      axisIdentifiers == nullptr || axisNames == nullptr ||
      axisKinds == nullptr || axisExtents == nullptr || axisRoles == nullptr ||
      variableIdentifiers == nullptr || variableNames == nullptr ||
      variableTypes == nullptr || variableLayouts == nullptr ||
      variableUnits == nullptr || variableDescriptions == nullptr ||
      variableRoles == nullptr || variableRaggedRoles == nullptr ||
      variableRaggedChildren == nullptr || dimensionCounts == nullptr ||
      dimensions == nullptr || attributeCounts == nullptr ||
      attributeNames == nullptr || attributeValues == nullptr)
    return invalid("Observation-schema manifest is incomplete.");
  const auto axisCount = axisIdentifiers->size();
  if (axisNames->size() != axisCount || axisKinds->size() != axisCount ||
      axisExtents->size() != axisCount || axisRoles->size() != axisCount)
    return invalid("Observation-schema manifest axes are inconsistent.");
  const auto variableCount = variableIdentifiers->size();
  if (variableNames->size() != variableCount ||
      variableTypes->size() != variableCount ||
      variableLayouts->size() != variableCount ||
      variableUnits->size() != variableCount ||
      variableDescriptions->size() != variableCount ||
      variableRoles->size() != variableCount ||
      variableRaggedRoles->size() != variableCount ||
      variableRaggedChildren->size() != variableCount ||
      dimensionCounts->size() != variableCount ||
      attributeCounts->size() != variableCount ||
      attributeNames->size() != attributeValues->size())
    return invalid("Observation-schema manifest variables are inconsistent.");
  if (metadataAttributeNames->size() != metadataAttributeValues->size() ||
      metadataStringListNames->size() != metadataStringListCounts->size())
    return invalid(
        "Observation-schema manifest metadata attributes are inconsistent.");
  std::size_t metadataStringValueCount = 0;
  for (const auto count : *metadataStringListCounts) {
    if (count < 0 ||
        static_cast<std::uint64_t>(count) >
            std::numeric_limits<std::size_t>::max() -
                metadataStringValueCount)
      return invalid(
          "Observation-schema manifest metadata string-list size is invalid.");
    metadataStringValueCount += static_cast<std::size_t>(count);
  }
  if (metadataStringValueCount != metadataStringListValues->size())
    return invalid(
        "Observation-schema manifest metadata string lists are inconsistent.");
  const auto metadataVariableCount = metadataVariableNames->size();
  if (metadataVariableIdentifiers->size() != metadataVariableCount ||
      metadataVariableTypes->size() != metadataVariableCount ||
      metadataVariableLogical->size() != metadataVariableCount ||
      metadataVariableDimensionCounts->size() != metadataVariableCount)
    return invalid(
        "Observation-schema manifest metadata variables are inconsistent.");
  std::size_t metadataExtentCount = 0;
  for (std::size_t index = 0; index < metadataVariableCount; ++index) {
    const auto type = (*metadataVariableTypes)[index];
    const auto dimensionCount = (*metadataVariableDimensionCounts)[index];
    if (type < 0 || type > 4 || (*metadataVariableLogical)[index] > 1 ||
        dimensionCount < 0 ||
        static_cast<std::uint64_t>(dimensionCount) >
            std::numeric_limits<std::size_t>::max() - metadataExtentCount)
      return invalid(
          "Observation-schema manifest metadata-variable value is invalid.");
    metadataExtentCount += static_cast<std::size_t>(dimensionCount);
  }
  if (metadataExtentCount != metadataVariableExtents->size())
    return invalid(
        "Observation-schema manifest metadata-variable extents are inconsistent.");
  if ((*schemaVersions)[0] <= 0 ||
      static_cast<std::uint64_t>((*schemaVersions)[0]) >
          std::numeric_limits<std::uint32_t>::max())
    return invalid("Observation-schema manifest version is invalid.");

  WVObservationSchema candidate;
  candidate.identifier = (*schemaIdentifiers)[0];
  candidate.version = static_cast<std::uint32_t>((*schemaVersions)[0]);
  candidate.preservesLegacyEncoding = (*preservesLegacyEncoding)[0] != 0;
  for (std::size_t index = 0; index < metadataAttributeNames->size(); ++index)
    candidate.metadata.attributes.push_back(
        {(*metadataAttributeNames)[index],
         (*metadataAttributeValues)[index]});
  std::size_t metadataStringOffset = 0;
  for (std::size_t index = 0; index < metadataStringListNames->size();
       ++index) {
    const auto count =
        static_cast<std::size_t>((*metadataStringListCounts)[index]);
    WVObservationStringListAttribute attribute;
    attribute.name = (*metadataStringListNames)[index];
    attribute.values.insert(
        attribute.values.end(),
        metadataStringListValues->begin() +
            static_cast<std::ptrdiff_t>(metadataStringOffset),
        metadataStringListValues->begin() + static_cast<std::ptrdiff_t>(
                                                metadataStringOffset + count));
    metadataStringOffset += count;
    candidate.metadata.stringListAttributes.push_back(std::move(attribute));
  }
  std::size_t metadataExtentOffset = 0;
  for (std::size_t index = 0; index < metadataVariableCount; ++index) {
    const auto dimensionCount = static_cast<std::size_t>(
        (*metadataVariableDimensionCounts)[index]);
    std::vector<std::size_t> extents;
    extents.reserve(dimensionCount);
    for (std::size_t dimension = 0; dimension < dimensionCount; ++dimension) {
      const auto extent =
          (*metadataVariableExtents)[metadataExtentOffset + dimension];
      if (extent < 0 ||
          static_cast<std::uint64_t>(extent) >
              std::numeric_limits<std::size_t>::max())
        return invalid(
            "Observation-schema manifest metadata-variable extent is invalid.");
      extents.push_back(static_cast<std::size_t>(extent));
    }
    metadataExtentOffset += dimensionCount;
    const auto storageName =
        "metadataVariable." + std::to_string(index);
    const auto type = static_cast<WVObservationScalarType>(
        (*metadataVariableTypes)[index]);
    WVObservationValue value;
    if (type == WVObservationScalarType::complex64) {
      const auto *real = manifest.value(storageName + ".real");
      const auto *imaginary = manifest.value(storageName + ".imag");
      if (real == nullptr || imaginary == nullptr ||
          real->dimensions != extents || imaginary->dimensions != extents ||
          !std::holds_alternative<std::vector<double>>(real->storage) ||
          !std::holds_alternative<std::vector<double>>(imaginary->storage))
        return invalid(
            "Observation-schema manifest complex metadata storage is invalid.");
      const auto &realValues = std::get<std::vector<double>>(real->storage);
      const auto &imaginaryValues =
          std::get<std::vector<double>>(imaginary->storage);
      if (realValues.size() != imaginaryValues.size())
        return invalid(
            "Observation-schema manifest complex metadata storage differs.");
      std::vector<WVComplex64> values(realValues.size());
      for (std::size_t element = 0; element < values.size(); ++element)
        values[element] = {realValues[element], imaginaryValues[element]};
      value = WVObservationValue::ownComplex(
          (*metadataVariableIdentifiers)[index], std::move(extents),
          std::move(values));
    } else {
      const auto *storage = manifest.value(storageName);
      if (storage == nullptr || storage->dimensions != extents)
        return invalid(
            "Observation-schema manifest metadata storage is invalid.");
      switch (type) {
      case WVObservationScalarType::real64:
        if (!std::holds_alternative<std::vector<double>>(storage->storage))
          return invalid(
              "Observation-schema manifest real metadata storage is invalid.");
        value = WVObservationValue::ownReal(
            (*metadataVariableIdentifiers)[index], std::move(extents),
            std::get<std::vector<double>>(storage->storage));
        break;
      case WVObservationScalarType::integer64:
        if (!std::holds_alternative<std::vector<std::int64_t>>(
                storage->storage))
          return invalid(
              "Observation-schema manifest integer metadata storage is invalid.");
        value = WVObservationValue::ownInteger(
            (*metadataVariableIdentifiers)[index], std::move(extents),
            std::get<std::vector<std::int64_t>>(storage->storage));
        break;
      case WVObservationScalarType::boolean8:
        if (!std::holds_alternative<std::vector<std::uint8_t>>(
                storage->storage))
          return invalid(
              "Observation-schema manifest Boolean metadata storage is invalid.");
        value = WVObservationValue::ownBoolean(
            (*metadataVariableIdentifiers)[index], std::move(extents),
            std::get<std::vector<std::uint8_t>>(storage->storage));
        break;
      case WVObservationScalarType::text:
        if (!std::holds_alternative<std::vector<std::string>>(
                storage->storage))
          return invalid(
              "Observation-schema manifest text metadata storage is invalid.");
        value = WVObservationValue::ownText(
            (*metadataVariableIdentifiers)[index], std::move(extents),
            std::get<std::vector<std::string>>(storage->storage));
        break;
      case WVObservationScalarType::complex64:
        return invalid(
            "Observation-schema manifest complex metadata storage is invalid.");
      }
    }
    candidate.metadata.variables.push_back(
        {(*metadataVariableNames)[index], std::move(value),
         (*metadataVariableLogical)[index] != 0});
  }
  for (std::size_t index = 0; index < axisCount; ++index) {
    if ((*axisKinds)[index] < 0 || (*axisKinds)[index] > 1 ||
        (*axisExtents)[index] < 0 || (*axisRoles)[index] < 0 ||
        (*axisRoles)[index] >
            static_cast<std::int64_t>(WVObservationCoordinateRole::profile))
      return invalid("Observation-schema manifest axis value is invalid.");
    candidate.axes.push_back(
        {(*axisIdentifiers)[index], (*axisNames)[index],
         static_cast<WVObservationAxisKind>((*axisKinds)[index]),
         static_cast<std::size_t>((*axisExtents)[index]),
         static_cast<WVObservationCoordinateRole>((*axisRoles)[index])});
  }
  std::size_t dimensionOffset = 0;
  std::size_t attributeOffset = 0;
  for (std::size_t index = 0; index < variableCount; ++index) {
    if ((*variableTypes)[index] < 0 || (*variableTypes)[index] > 4 ||
        (*variableLayouts)[index] < 0 || (*variableLayouts)[index] > 3 ||
        (*variableRoles)[index] < 0 ||
        (*variableRoles)[index] >
            static_cast<std::int64_t>(WVObservationCoordinateRole::profile) ||
        (*variableRaggedRoles)[index] < 0 ||
        (*variableRaggedRoles)[index] > 2 || (*dimensionCounts)[index] < 0 ||
        (*attributeCounts)[index] < 0)
      return invalid("Observation-schema manifest variable value is invalid.");
    const auto dimensionCount =
        static_cast<std::size_t>((*dimensionCounts)[index]);
    const auto attributeCount =
        static_cast<std::size_t>((*attributeCounts)[index]);
    if (dimensionCount > dimensions->size() - dimensionOffset ||
        attributeCount > attributeNames->size() - attributeOffset)
      return invalid("Observation-schema manifest offsets are invalid.");
    WVObservationVariable variable;
    variable.identifier = (*variableIdentifiers)[index];
    variable.name = (*variableNames)[index];
    variable.scalarType =
        static_cast<WVObservationScalarType>((*variableTypes)[index]);
    variable.layout =
        static_cast<WVObservationValueLayout>((*variableLayouts)[index]);
    variable.units = (*variableUnits)[index];
    variable.description = (*variableDescriptions)[index];
    variable.coordinateRole =
        static_cast<WVObservationCoordinateRole>((*variableRoles)[index]);
    variable.raggedRole =
        static_cast<WVObservationRaggedRole>((*variableRaggedRoles)[index]);
    variable.raggedChildAxisIdentifier = (*variableRaggedChildren)[index];
    variable.dimensionIdentifiers.insert(
        variable.dimensionIdentifiers.end(),
        dimensions->begin() + static_cast<std::ptrdiff_t>(dimensionOffset),
        dimensions->begin() +
            static_cast<std::ptrdiff_t>(dimensionOffset + dimensionCount));
    for (std::size_t attribute = 0; attribute < attributeCount; ++attribute)
      variable.attributes.push_back(
          {(*attributeNames)[attributeOffset + attribute],
           (*attributeValues)[attributeOffset + attribute]});
    dimensionOffset += dimensionCount;
    attributeOffset += attributeCount;
    candidate.variables.push_back(std::move(variable));
  }
  if (dimensionOffset != dimensions->size() ||
      attributeOffset != attributeNames->size())
    return invalid("Observation-schema manifest has trailing values.");
  status = validateObservationSchema(candidate);
  if (!status)
    return status;
  schema = std::move(candidate);
  return WVKernelStatus::ok();
}

} // namespace wavevortex::runtime
