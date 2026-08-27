#include "WVLegacyObservationNetCDFAdapter.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVLegacyObserverCompatibility.hpp"
#include "WVModelOutputTransformAdapters.hpp"
#include "WVNetCDF.hpp"
#include "WVObserverAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <netcdf.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime {
namespace {

bool hexDecode(const std::string &text, std::vector<std::uint8_t> &bytes) {
  if (text.size() % 2 != 0)
    return false;
  const auto digit = [](char value) {
    if (value >= '0' && value <= '9')
      return value - '0';
    if (value >= 'a' && value <= 'f')
      return value - 'a' + 10;
    if (value >= 'A' && value <= 'F')
      return value - 'A' + 10;
    return -1;
  };
  bytes.resize(text.size() / 2);
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    const int high = digit(text[2 * index]);
    const int low = digit(text[2 * index + 1]);
    if (high < 0 || low < 0)
      return false;
    bytes[index] = static_cast<std::uint8_t>((high << 4) | low);
  }
  return true;
}

WVCheckpointStatus failure(WVCheckpointStatusCode code, std::string message,
                           std::string location) {
  return {code, std::move(message), std::move(location)};
}

WVCheckpointStatus optionalTextAttribute(int group, const char *name,
                                         std::string &value, bool &present,
                                         const std::string &path) {
  nc_type type = NC_NAT;
  std::size_t length = 0;
  const int inquiry = nc_inq_att(group, NC_GLOBAL, name, &type, &length);
  if (inquiry == NC_ENOTATT) {
    present = false;
    value.clear();
    return WVCheckpointStatus::ok();
  }
  if (inquiry != NC_NOERR)
    return detail::netcdfFailure(inquiry, "Attribute inspection",
                                 path + "/@" + name);
  present = true;
  return detail::readTextAttribute(group, name, value, path);
}

WVCheckpointStatus stringListAttribute(int group, const char *name,
                                       std::vector<std::string> &values,
                                       bool &present, const std::string &path) {
  nc_type type = NC_NAT;
  std::size_t length = 0;
  const int inquiry = nc_inq_att(group, NC_GLOBAL, name, &type, &length);
  if (inquiry == NC_ENOTATT) {
    present = false;
    values.clear();
    return WVCheckpointStatus::ok();
  }
  if (inquiry != NC_NOERR)
    return detail::netcdfFailure(inquiry, "String-list inspection",
                                 path + "/@" + name);
  present = true;
  values.clear();
  if (type == NC_STRING) {
    std::vector<char *> raw(length, nullptr);
    const int result = nc_get_att_string(group, NC_GLOBAL, name, raw.data());
    if (result != NC_NOERR)
      return detail::netcdfFailure(result, "String-list read",
                                   path + "/@" + name);
    values.reserve(length);
    for (const auto *entry : raw)
      values.emplace_back(entry == nullptr ? "" : entry);
    nc_free_string(length, raw.data());
    return WVCheckpointStatus::ok();
  }
  if (type == NC_CHAR) {
    std::string raw(length, '\0');
    const int result = nc_get_att_text(group, NC_GLOBAL, name, raw.data());
    if (result != NC_NOERR)
      return detail::netcdfFailure(result, "String-list read",
                                   path + "/@" + name);
    values.push_back(std::move(raw));
    return WVCheckpointStatus::ok();
  }
  return failure(WVCheckpointStatusCode::typeMismatch,
                 "Observer field-name attributes must contain strings.",
                 path + "/@" + name);
}

std::string portableIdentifier(std::string value) {
  for (auto &character : value) {
    const bool valid = (character >= 'a' && character <= 'z') ||
                       (character >= 'A' && character <= 'Z') ||
                       (character >= '0' && character <= '9') ||
                       character == '-' || character == '_' || character == '.';
    if (!valid)
      character = '-';
  }
  if (value.empty())
    value = "unnamed";
  return value;
}

bool sameObserverConfiguration(const WVObserverRecord &left,
                               const WVObserverRecord &right) {
  const auto sameTypedConfiguration = [&]() {
    if (left.configuration.schemaIdentifier.empty() &&
        right.configuration.schemaIdentifier.empty())
      return true;
    std::vector<std::uint8_t> leftBytes;
    std::vector<std::uint8_t> rightBytes;
    return static_cast<bool>(
               encodePortableTypedRecord(left.configuration, leftBytes)) &&
           static_cast<bool>(
               encodePortableTypedRecord(right.configuration, rightBytes)) &&
           leftBytes == rightBytes;
  };
  return left.name == right.name &&
         left.typeIdentifier == right.typeIdentifier &&
         left.contractVersion == right.contractVersion &&
         sameTypedConfiguration() &&
         left.stateBlockIdentifiers == right.stateBlockIdentifiers &&
         left.fieldNames == right.fieldNames && left.x == right.x &&
         left.y == right.y && left.z == right.z &&
         left.isXYOnly == right.isXYOnly &&
         left.shouldAntialias == right.shouldAntialias &&
         left.advectionInterpolation == right.advectionInterpolation &&
         left.trackedFieldInterpolation == right.trackedFieldInterpolation &&
         left.horizontalAbsoluteTolerance ==
             right.horizontalAbsoluteTolerance &&
         left.verticalAbsoluteTolerance == right.verticalAbsoluteTolerance &&
         left.outputScale == right.outputScale &&
         left.outputOffset == right.outputOffset;
}

bool sameRestartValues(const std::vector<double> &left,
                       const std::vector<double> &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    const double scale =
        std::max({1.0, std::abs(left[index]), std::abs(right[index])});
    if (!std::isfinite(left[index]) || !std::isfinite(right[index]) ||
        std::abs(left[index] - right[index]) >
            32 * std::numeric_limits<double>::epsilon() * scale)
      return false;
  }
  return true;
}

WVCheckpointStatus variableShape(int group, const std::string &name,
                                 std::vector<std::string> &dimensionNames,
                                 std::vector<std::size_t> &dimensionLengths,
                                 const std::string &path) {
  int variable = -1;
  auto result =
      detail::checkedNetCDF(nc_inq_varid(group, name.c_str(), &variable),
                            "Variable lookup", path + "/" + name);
  if (!result)
    return result;
  int count = 0;
  result = detail::checkedNetCDF(nc_inq_varndims(group, variable, &count),
                                 "Variable-rank inspection", path + "/" + name);
  if (!result)
    return result;
  std::vector<int> dimensions(static_cast<std::size_t>(count));
  result =
      detail::checkedNetCDF(nc_inq_vardimid(group, variable, dimensions.data()),
                            "Variable-dimension inspection", path + "/" + name);
  if (!result)
    return result;
  dimensionNames.clear();
  dimensionLengths.clear();
  for (const int dimension : dimensions) {
    char rawName[NC_MAX_NAME + 1] = {};
    std::size_t length = 0;
    result = detail::checkedNetCDF(
        nc_inq_dim(group, dimension, rawName, &length),
        "Variable-dimension inspection", path + "/" + name);
    if (!result)
      return result;
    dimensionNames.emplace_back(rawName);
    dimensionLengths.push_back(length);
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus readLatestRealSlab(int group, const std::string &name,
                                      std::size_t recordIndex,
                                      std::vector<double> &values,
                                      const std::string &path) {
  int variable = -1;
  auto result =
      detail::checkedNetCDF(nc_inq_varid(group, name.c_str(), &variable),
                            "Variable lookup", path + "/" + name);
  if (!result)
    return result;
  std::vector<std::string> dimensionNames;
  std::vector<std::size_t> dimensions;
  result = variableShape(group, name, dimensionNames, dimensions, path);
  if (!result)
    return result;
  if (dimensions.empty() || dimensionNames.front() != "t")
    return failure(WVCheckpointStatusCode::shapeMismatch,
                   "Dynamic observer values must begin with time.",
                   path + "/" + name);
  std::vector<std::size_t> start(dimensions.size(), 0);
  std::vector<std::size_t> count = dimensions;
  start.front() = recordIndex;
  count.front() = 1;
  std::size_t elementCount = 1;
  for (std::size_t index = 1; index < count.size(); ++index)
    elementCount *= count[index];
  values.resize(elementCount);
  return detail::checkedNetCDF(nc_get_vara_double(group, variable, start.data(),
                                                  count.data(), values.data()),
                               "Dynamic observer-state read",
                               path + "/" + name);
}

WVCheckpointStatus recordAtTime(int group, double selectedTime,
                                std::size_t &recordIndex, bool &found,
                                const std::string &path) {
  std::size_t count = 0;
  auto result = detail::dimensionLength(group, "t", count, path);
  if (!result)
    return result;
  int variable = -1;
  result = detail::checkedNetCDF(nc_inq_varid(group, "t", &variable),
                                 "Output-time lookup", path + "/t");
  if (!result)
    return result;
  std::vector<double> times(count);
  result =
      detail::checkedNetCDF(nc_get_var_double(group, variable, times.data()),
                            "Output-time read", path + "/t");
  if (!result)
    return result;
  found = false;
  for (std::size_t index = 0; index < times.size(); ++index) {
    const double scale =
        std::max({1.0, std::abs(times[index]), std::abs(selectedTime)});
    if (std::abs(times[index] - selectedTime) <=
        8 * std::numeric_limits<double>::epsilon() * scale) {
      if (found)
        return failure(WVCheckpointStatusCode::ambiguousState,
                       "An output group contains the selected restart time "
                       "more than once.",
                       path + "/t");
      found = true;
      recordIndex = index;
    }
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus
readObserverStateAtTime(const std::string &filePath,
                        const WVOutputGroupRecord &groupRecord,
                        const WVObserverRecord &observer,
                        const WVObserverFactoryRegistration &registration,
                        double selectedTime,
                        std::vector<std::vector<double>> &values, bool &found) {
  detail::WVNetCDFFile file;
  auto result = detail::WVNetCDFFile::openReadOnly(filePath, file);
  if (!result)
    return result;
  int group = -1;
  const std::string groupPath = "/" + groupRecord.name;
  result = detail::checkedNetCDF(
      nc_inq_ncid(file.id(), groupRecord.name.c_str(), &group),
      "Observer restart-group lookup", groupPath);
  if (!result)
    return result;
  std::size_t recordIndex = 0;
  result = recordAtTime(group, selectedTime, recordIndex, found, groupPath);
  if (!result || !found)
    return result;

  const auto *compatibility =
      detail::WVObserverFactoryRegistrationAccess::legacyCompatibility(
          registration);
  if (compatibility == nullptr || !compatibility->operationResolver)
    return failure(WVCheckpointStatusCode::unsupportedObserver,
                   "Dynamic observer state uses an unsupported observer.",
                   groupPath);
  values.clear();
  detail::WVLegacyObserverOperationBinder operations;
  operations.movingPositions = [&]() -> WVKernelStatus {
    const auto channels = detail::particlePositionChannels(observer.isXYOnly);
    values.resize(channels.size());
    for (std::size_t index = 0; index < channels.size(); ++index) {
      result = readLatestRealSlab(
          group,
          observer.name + "_" + detail::movingFieldChannelName(channels[index]),
          recordIndex, values[index], groupPath);
      if (!result)
        return {WVKernelStatusCode::invalidConfiguration, result.message};
    }
    if (observer.isXYOnly && !observer.z.empty()) {
      values.emplace_back();
      result = readLatestRealSlab(group, observer.name + "_z", recordIndex,
                                  values.back(), groupPath);
      if (!result)
        return {WVKernelStatusCode::invalidConfiguration, result.message};
    }
    return WVKernelStatus::ok();
  };
  operations.integratedState = [&]() -> WVKernelStatus {
    values.resize(1);
    result = readLatestRealSlab(group, observer.name, recordIndex, values[0],
                                groupPath);
    if (!result)
      return {WVKernelStatusCode::invalidConfiguration, result.message};
    return WVKernelStatus::ok();
  };
  const auto noDynamicState = [&]() {
    found = false;
    return WVKernelStatus::ok();
  };
  operations.fullField = noDynamicState;
  operations.fixedVerticalProfiles = noDynamicState;
  operations.fixedPositions = noDynamicState;
  const auto operationStatus =
      compatibility->operationResolver(observer, operations);
  if (!operationStatus && result)
    return failure(WVCheckpointStatusCode::unsupportedObserver,
                   operationStatus.message, groupPath);
  if (!result)
    return result;
  return WVCheckpointStatus::ok();
}

struct RestartStateLocation {
  const WVOutputFileRecord *file = nullptr;
  const WVOutputGroupRecord *group = nullptr;
  std::size_t recordIndex = 0;
};

std::string restartStateLocation(const RestartStateLocation &location) {
  return location.file->destination + "/" + location.group->name;
}

WVCheckpointStatus inspectRestartSlab(const RestartStateLocation &location,
                                      const std::string &name) {
  detail::WVNetCDFFile file;
  auto result = detail::WVNetCDFFile::openReadOnly(
      location.file->destination, file);
  if (!result)
    return result;
  int group = -1;
  const auto path = restartStateLocation(location);
  result = detail::checkedNetCDF(
      nc_inq_ncid(file.id(), location.group->name.c_str(), &group),
      "Observer restart-group lookup", path);
  if (!result)
    return result;
  std::vector<std::string> dimensionNames;
  std::vector<std::size_t> dimensionLengths;
  result = variableShape(group, name, dimensionNames, dimensionLengths, path);
  if (!result)
    return result;
  if (dimensionNames.empty() || dimensionNames.front() != "t" ||
      location.recordIndex >= dimensionLengths.front())
    return failure(WVCheckpointStatusCode::shapeMismatch,
                   "Dynamic observer values must begin with time.",
                   path + "/" + name);
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus compareRestartSlabs(
    const RestartStateLocation &left, const RestartStateLocation &right,
    const std::string &name) {
  detail::WVNetCDFFile leftFile;
  detail::WVNetCDFFile rightFile;
  auto result = detail::WVNetCDFFile::openReadOnly(left.file->destination,
                                                   leftFile);
  if (!result)
    return result;
  result = detail::WVNetCDFFile::openReadOnly(right.file->destination,
                                              rightFile);
  if (!result)
    return result;
  int leftGroup = -1;
  int rightGroup = -1;
  const auto leftPath = restartStateLocation(left);
  const auto rightPath = restartStateLocation(right);
  result = detail::checkedNetCDF(
      nc_inq_ncid(leftFile.id(), left.group->name.c_str(), &leftGroup),
      "Observer restart-group lookup", leftPath);
  if (!result)
    return result;
  result = detail::checkedNetCDF(
      nc_inq_ncid(rightFile.id(), right.group->name.c_str(), &rightGroup),
      "Observer restart-group lookup", rightPath);
  if (!result)
    return result;
  int leftVariable = -1;
  int rightVariable = -1;
  result = detail::checkedNetCDF(
      nc_inq_varid(leftGroup, name.c_str(), &leftVariable),
      "Dynamic observer-state lookup", leftPath + "/" + name);
  if (!result)
    return result;
  result = detail::checkedNetCDF(
      nc_inq_varid(rightGroup, name.c_str(), &rightVariable),
      "Dynamic observer-state lookup", rightPath + "/" + name);
  if (!result)
    return result;
  std::vector<std::string> leftNames;
  std::vector<std::string> rightNames;
  std::vector<std::size_t> leftDimensions;
  std::vector<std::size_t> rightDimensions;
  result = variableShape(leftGroup, name, leftNames, leftDimensions, leftPath);
  if (!result)
    return result;
  result =
      variableShape(rightGroup, name, rightNames, rightDimensions, rightPath);
  if (!result)
    return result;
  if (leftNames != rightNames || leftNames.empty() ||
      leftNames.front() != "t" || leftDimensions.size() != rightDimensions.size() ||
      !std::equal(leftDimensions.begin() + 1, leftDimensions.end(),
                  rightDimensions.begin() + 1) ||
      left.recordIndex >= leftDimensions.front() ||
      right.recordIndex >= rightDimensions.front())
    return failure(WVCheckpointStatusCode::ambiguousState,
                   "Dynamic observer state has incompatible shapes across "
                   "output groups.",
                   leftPath + " and " + rightPath);

  if (leftDimensions.size() == 1) {
    const std::size_t leftIndex[] = {left.recordIndex};
    const std::size_t rightIndex[] = {right.recordIndex};
    double leftValue = 0.0;
    double rightValue = 0.0;
    result = detail::checkedNetCDF(
        nc_get_var1_double(leftGroup, leftVariable, leftIndex, &leftValue),
        "Dynamic observer-state read", leftPath + "/" + name);
    if (!result)
      return result;
    result = detail::checkedNetCDF(
        nc_get_var1_double(rightGroup, rightVariable, rightIndex, &rightValue),
        "Dynamic observer-state read", rightPath + "/" + name);
    if (!result)
      return result;
    const double scale =
        std::max({1.0, std::abs(leftValue), std::abs(rightValue)});
    if (!std::isfinite(leftValue) || !std::isfinite(rightValue) ||
        std::abs(leftValue - rightValue) >
            32 * std::numeric_limits<double>::epsilon() * scale)
      return failure(WVCheckpointStatusCode::ambiguousState,
                     "Dynamic observer state conflicts across output groups.",
                     leftPath + " and " + rightPath);
    return WVCheckpointStatus::ok();
  }

  constexpr std::size_t maximumChunkElements = 4096;
  const std::size_t lastDimension = leftDimensions.size() - 1;
  const std::size_t innerCount = leftDimensions[lastDimension];
  std::size_t outerCount = 1;
  for (std::size_t dimension = 1; dimension < lastDimension; ++dimension)
    outerCount *= leftDimensions[dimension];
  std::vector<double> leftValues(std::min(innerCount, maximumChunkElements));
  std::vector<double> rightValues(leftValues.size());
  std::vector<std::size_t> leftStart(leftDimensions.size(), 0);
  std::vector<std::size_t> rightStart(rightDimensions.size(), 0);
  std::vector<std::size_t> count(leftDimensions.size(), 1);
  leftStart.front() = left.recordIndex;
  rightStart.front() = right.recordIndex;
  for (std::size_t outer = 0; outer < outerCount; ++outer) {
    std::size_t remainder = outer;
    for (std::size_t dimension = lastDimension; dimension-- > 1;) {
      leftStart[dimension] = remainder % leftDimensions[dimension];
      rightStart[dimension] = leftStart[dimension];
      remainder /= leftDimensions[dimension];
    }
    for (std::size_t offset = 0; offset < innerCount;
         offset += maximumChunkElements) {
      const auto chunk =
          std::min(maximumChunkElements, innerCount - offset);
      leftStart[lastDimension] = offset;
      rightStart[lastDimension] = offset;
      count[lastDimension] = chunk;
      result = detail::checkedNetCDF(
          nc_get_vara_double(leftGroup, leftVariable, leftStart.data(),
                             count.data(), leftValues.data()),
          "Dynamic observer-state read", leftPath + "/" + name);
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_get_vara_double(rightGroup, rightVariable, rightStart.data(),
                             count.data(), rightValues.data()),
          "Dynamic observer-state read", rightPath + "/" + name);
      if (!result)
        return result;
      for (std::size_t index = 0; index < chunk; ++index) {
        const double scale =
            std::max({1.0, std::abs(leftValues[index]),
                      std::abs(rightValues[index])});
        if (!std::isfinite(leftValues[index]) ||
            !std::isfinite(rightValues[index]) ||
            std::abs(leftValues[index] - rightValues[index]) >
                32 * std::numeric_limits<double>::epsilon() * scale)
          return failure(
              WVCheckpointStatusCode::ambiguousState,
              "Dynamic observer state conflicts across output groups.",
              leftPath + " and " + rightPath);
      }
    }
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus readRestartConfigurationSlab(
    const RestartStateLocation &location, const std::string &name,
    std::vector<double> &values) {
  detail::WVNetCDFFile file;
  auto result = detail::WVNetCDFFile::openReadOnly(
      location.file->destination, file);
  if (!result)
    return result;
  int group = -1;
  const auto path = restartStateLocation(location);
  result = detail::checkedNetCDF(
      nc_inq_ncid(file.id(), location.group->name.c_str(), &group),
      "Observer restart-group lookup", path);
  if (!result)
    return result;
  return readLatestRealSlab(group, name, location.recordIndex, values, path);
}

WVCheckpointStatus readWholeDoubleVariable(int group, const std::string &name,
                                           std::vector<double> &values,
                                           const std::string &path) {
  int variable = -1;
  auto result =
      detail::checkedNetCDF(nc_inq_varid(group, name.c_str(), &variable),
                            "Variable lookup", path + "/" + name);
  if (!result)
    return result;
  std::vector<std::string> names;
  std::vector<std::size_t> dimensions;
  result = variableShape(group, name, names, dimensions, path);
  if (!result)
    return result;
  std::size_t count = 1;
  for (const auto dimension : dimensions)
    count *= dimension;
  values.resize(count);
  return detail::checkedNetCDF(
      nc_get_var_double(group, variable, values.data()), "Variable read",
      path + "/" + name);
}

WVCheckpointStatus mergeStateBlock(WVPortableObserverRecord &record,
                                   const WVStateBlockRecord &block) {
  const auto found =
      std::find_if(record.stateBlocks.begin(), record.stateBlocks.end(),
                   [&](const auto &candidate) {
                     return candidate.identifier == block.identifier;
                   });
  if (found == record.stateBlocks.end()) {
    record.stateBlocks.push_back(block);
    return WVCheckpointStatus::ok();
  }
  if (found->scalarType != block.scalarType ||
      found->dimensions != block.dimensions ||
      found->toleranceKind != block.toleranceKind ||
      found->absoluteTolerance != block.absoluteTolerance ||
      found->ownership != block.ownership ||
      found->restartRequirement != block.restartRequirement)
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "Shared observer state-block definitions conflict.",
                   "/observingSystems");
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus mergeObserver(WVPortableObserverRecord &record,
                                 WVObserverRecord observer) {
  const auto found =
      std::find_if(record.observers.begin(), record.observers.end(),
                   [&](const auto &candidate) {
                     return candidate.identifier == observer.identifier;
                   });
  if (found == record.observers.end()) {
    record.observers.push_back(std::move(observer));
    return WVCheckpointStatus::ok();
  }
  if (!sameObserverConfiguration(*found, observer))
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "Shared observer metadata conflicts across output groups.",
                   "/observingSystems");
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus parseObserver(int outputGroup, int metadataGroup,
                                 const std::string &outputPath,
                                 const WVExtensionCatalog &catalog,
                                 WVPortableObserverRecord &portable,
                                 std::string &identifier) {
  std::string className;
  auto result = detail::readTextAttribute(metadataGroup, "AnnotatedClass",
                                          className, outputPath);
  if (!result)
    return result;
  WVObserverRecord observer;
  observer.typeIdentifier = className;
  bool present = false;
  std::string schemaIdentifier;
  result = optionalTextAttribute(metadataGroup,
                                 "portableObservationSchemaIdentifier",
                                 schemaIdentifier, present, outputPath);
  if (!result)
    return result;
  const bool hasPortableSchema = present;
  if (hasPortableSchema) {
    unsigned int schemaVersion = 0;
    result = detail::checkedNetCDF(
        nc_get_att_uint(metadataGroup, NC_GLOBAL,
                        "portableObservationSchemaVersion", &schemaVersion),
        "Observation-schema version read", outputPath);
    if (!result)
      return result;
    unsigned int observerVersion = 0;
    result = detail::checkedNetCDF(
        nc_get_att_uint(metadataGroup, NC_GLOBAL,
                        "portableObserverContractVersion", &observerVersion),
        "Observer contract-version read", outputPath);
    if (!result)
      return result;
    observer.contractVersion = observerVersion;
    std::string encodedConfiguration;
    result = detail::readTextAttribute(metadataGroup,
                                       "portableObserverConfiguration",
                                       encodedConfiguration, outputPath);
    if (!result)
      return result;
    std::vector<std::uint8_t> configurationBytes;
    if (!hexDecode(encodedConfiguration, configurationBytes))
      return failure(WVCheckpointStatusCode::invalidValue,
                     "Portable observer configuration is malformed.",
                     outputPath);
    const auto decoded = decodePortableTypedRecord(
        configurationBytes, observer.configuration,
        {1024 * 1024, true, true});
    if (!decoded)
      return failure(WVCheckpointStatusCode::invalidValue, decoded.message,
                     outputPath);
    result = detail::readTextAttribute(metadataGroup, "name", observer.name,
                                       outputPath);
    if (!result)
      return result;
    result = detail::readTextAttribute(metadataGroup, "portableIdentifier",
                                       observer.identifier, outputPath);
    if (!result)
      return result;
    identifier = observer.identifier;
    return mergeObserver(portable, std::move(observer));
  }
  observer.contractVersion = WVPortablePairContractVersion;
  const auto *registration = catalog.observers().registration(
      observer.typeIdentifier, observer.contractVersion);
  const auto *compatibility =
      registration == nullptr
          ? nullptr
          : detail::WVObserverFactoryRegistrationAccess::legacyCompatibility(
                *registration);
  if (compatibility == nullptr || !compatibility->operationResolver)
    return failure(WVCheckpointStatusCode::unsupportedObserver,
                   "Unsupported MATLAB observing-system class '" + className +
                       "'.",
                   outputPath + "/@AnnotatedClass");
  WVObserverExecutionPlan execution;
  execution.fieldListAttribute =
      compatibility->persistence.fieldListAttribute;
  execution.coefficientRestartFamilies =
      compatibility->persistence.coefficientRestartFamilies;
  result = optionalTextAttribute(metadataGroup, "name", observer.name, present,
                                 outputPath);
  if (!result)
    return result;
  if (!present)
    observer.name = className;
  result = optionalTextAttribute(metadataGroup, "portableIdentifier",
                                 observer.identifier, present, outputPath);
  if (!result)
    return result;
  const bool hadPortableIdentifier = present;
  if (!hadPortableIdentifier) {
    if (!compatibility->persistence.defaultIdentifier.empty())
      observer.identifier = compatibility->persistence.defaultIdentifier;
    else
      observer.identifier = portableIdentifier(className + "-" + observer.name);
  }
  identifier = observer.identifier;

  if (!execution.fieldListAttribute.empty()) {
    result = stringListAttribute(metadataGroup,
                                 execution.fieldListAttribute.c_str(),
                                 observer.fieldNames, present, outputPath);
    if (!result)
      return result;
    if (!present)
      observer.fieldNames.clear();
    else
      observer.fieldNames.erase(std::remove(observer.fieldNames.begin(),
                                            observer.fieldNames.end(),
                                            std::string{}),
                                observer.fieldNames.end());
    if (!hadPortableIdentifier &&
        compatibility->persistence.appendFieldsToDefaultIdentifier) {
      for (const auto &field : observer.fieldNames)
        observer.identifier += "-" + portableIdentifier(field);
      identifier = observer.identifier;
    }
  }

  const bool canonicalCoefficientProvider =
      execution.coefficientRestartFamilies ==
          std::vector<std::string>({"Ap", "Am", "A0"}) &&
      execution.fieldListAttribute.empty();
  WVCheckpointStatus operationResult = WVCheckpointStatus::ok();
  const auto adaptOperation = [&](auto action) {
    operationResult = action();
    return operationResult
               ? WVKernelStatus::ok()
               : WVKernelStatus{WVKernelStatusCode::invalidConfiguration,
                                operationResult.message};
  };
  detail::WVLegacyObserverOperationBinder operations;
  operations.fullField = [&]() {
    return adaptOperation([&]() -> WVCheckpointStatus {
    if (canonicalCoefficientProvider) {
    double tolerance = 0.0;
    result = detail::readDoubleScalar(metadataGroup, "absTolerance", tolerance,
                                      outputPath);
    if (!result)
      return result;
    int apVariable = -1;
    bool hasAp = false;
    result = detail::variableIdIfPresent(outputGroup, "Ap_real", apVariable,
                                         hasAp, outputPath);
    if (!result)
      return result;
    observer.stateBlockIdentifiers =
        hasAp ? std::vector<std::string>{"Ap", "Am", "A0"}
              : std::vector<std::string>{"A0"};
    std::vector<std::string> names;
    std::vector<std::size_t> dimensions;
    result = variableShape(outputGroup, hasAp ? "Ap_real" : "A0_real",
                           names, dimensions, outputPath);
    if (!result)
      return result;
    if ((hasAp && (names.size() != 3 || names[0] != "t" ||
                   names[1] != "kl" || names[2] != "j")) ||
        (!hasAp && (names.size() != 2 || names[0] != "t" ||
                    names[1] != "kl")))
      return failure(
          WVCheckpointStatusCode::shapeMismatch,
          hasAp ? "Coefficient output must use [t,kl,j] NetCDF order."
                : "Compact A0 output must use [t,kl] NetCDF order.",
          outputPath + (hasAp ? "/Ap_real" : "/A0_real"));
    const std::vector<std::size_t> logical =
        hasAp ? std::vector<std::size_t>{dimensions[2], dimensions[1]}
              : std::vector<std::size_t>{dimensions[1]};
    for (const auto &name : observer.stateBlockIdentifiers) {
      result = mergeStateBlock(portable,
                               {name, WVStateScalarType::complex64, logical,
                                WVToleranceKind::coefficientEnergyScaled,
                                tolerance, WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
      if (!result)
        return result;
    }
    } else {
    const bool hasLegacyCoefficients =
        std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  "Ap") != observer.fieldNames.end() &&
        std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  "Am") != observer.fieldNames.end() &&
        std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  "A0") != observer.fieldNames.end();
    const bool hasCompactA0 = !hasLegacyCoefficients &&
        std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  "A0") != observer.fieldNames.end();
    if (hasLegacyCoefficients || hasCompactA0) {
      std::vector<std::string> names;
      std::vector<std::size_t> dimensions;
      int coefficient = -1;
      bool paired = false;
      const char *baseName = hasLegacyCoefficients ? "Ap" : "A0";
      result = detail::variableIdIfPresent(
          outputGroup, std::string(baseName) + "_real", coefficient, paired,
          outputPath);
      if (!result)
        return result;
      result = variableShape(outputGroup,
                             paired ? std::string(baseName) + "_real"
                                    : std::string(baseName),
                             names, dimensions, outputPath);
      if (!result)
        return result;
      const bool legacyShape =
          hasLegacyCoefficients &&
          (names.size() == 2 || names.size() == 3) &&
          names[names.size() - 2] == "kl" && names.back() == "j";
      const bool compactShape = hasCompactA0 &&
          (names.size() == 1 || names.size() == 2) &&
          names.back() == "kl";
      if (!legacyShape && !compactShape)
        return failure(WVCheckpointStatusCode::shapeMismatch,
                       "Eulerian coefficients use an incompatible NetCDF "
                       "order.",
                       outputPath + "/" + baseName);
      const std::vector<std::size_t> logical =
          hasLegacyCoefficients
              ? std::vector<std::size_t>{
                    dimensions.back(), dimensions[dimensions.size() - 2]}
              : std::vector<std::size_t>{dimensions.back()};
      const std::vector<std::string> families =
          hasLegacyCoefficients
              ? std::vector<std::string>{"Ap", "Am", "A0"}
              : std::vector<std::string>{"A0"};
      for (const auto &name : families) {
        const auto existing = std::find_if(
            portable.stateBlocks.begin(), portable.stateBlocks.end(),
            [&](const auto &block) { return block.identifier == name; });
        const double tolerance = existing == portable.stateBlocks.end()
                                     ? 1e-6
                                     : existing->absoluteTolerance;
        result = mergeStateBlock(portable,
                                 {name, WVStateScalarType::complex64, logical,
                                  WVToleranceKind::coefficientEnergyScaled,
                                  tolerance, WVStateOwnership::integratorOwned,
                                  WVRestartRequirement::requiredDynamicState});
        if (!result)
          return result;
      }
    }
    }
    return WVCheckpointStatus::ok();
    });
  };
  operations.movingPositions = [&]() {
    return adaptOperation([&]() -> WVCheckpointStatus {
    result = detail::readLogicalScalar(metadataGroup, "isXYOnly",
                                       observer.isXYOnly, outputPath);
    if (!result)
      return result;
    result = detail::readDoubleScalar(metadataGroup, "absToleranceXY",
                                      observer.horizontalAbsoluteTolerance,
                                      outputPath);
    if (!result)
      return result;
    result = detail::readDoubleScalar(metadataGroup, "absToleranceZ",
                                      observer.verticalAbsoluteTolerance,
                                      outputPath);
    if (!result)
      return result;
    std::string interpolation;
    result = optionalTextAttribute(metadataGroup, "advectionInterpolation",
                                   interpolation, present, outputPath);
    if (!result)
      return result;
    observer.advectionInterpolation = interpolation == "spline"
                                          ? WVPositionInterpolation::spline
                                          : WVPositionInterpolation::linear;
    result = optionalTextAttribute(metadataGroup, "trackedVarInterpolation",
                                   interpolation, present, outputPath);
    if (!result)
      return result;
    observer.trackedFieldInterpolation = interpolation == "spline"
                                             ? WVPositionInterpolation::spline
                                             : WVPositionInterpolation::linear;
    std::size_t count = 0;
    result = detail::dimensionLength(outputGroup, observer.name + "_id", count,
                                     outputPath);
    if (!result)
      return result;
    observer.x.assign(count, 0.0);
    observer.y.assign(count, 0.0);
    int fixedZVariable = -1;
    if (!observer.isXYOnly ||
        nc_inq_varid(outputGroup, (observer.name + "_z").c_str(),
                     &fixedZVariable) == NC_NOERR)
      observer.z.assign(count, 0.0);
    const auto channels = detail::particlePositionChannels(observer.isXYOnly);
    for (std::size_t index = 0; index < channels.size(); ++index) {
      const std::string block = observer.identifier + "-" +
                                detail::movingFieldChannelName(channels[index]);
      observer.stateBlockIdentifiers.push_back(block);
      result = mergeStateBlock(portable,
                               {block,
                                WVStateScalarType::real64,
                                {count},
                                WVToleranceKind::uniformAbsolute,
                                index < 2 ? observer.horizontalAbsoluteTolerance
                                          : observer.verticalAbsoluteTolerance,
                                WVStateOwnership::integratorOwned,
                                WVRestartRequirement::requiredDynamicState});
      if (!result)
        return result;
    }
    return WVCheckpointStatus::ok();
    });
  };
  operations.integratedState = [&]() {
    return adaptOperation([&]() -> WVCheckpointStatus {
    result = detail::readLogicalScalar(metadataGroup, "isXYOnly",
                                       observer.isXYOnly, outputPath);
    if (!result)
      return result;
    result = detail::readLogicalScalar(metadataGroup, "shouldAntialias",
                                       observer.shouldAntialias, outputPath);
    if (!result)
      return result;
    double tolerance = 0.0;
    result = detail::readDoubleScalar(metadataGroup, "absTolerance", tolerance,
                                      outputPath);
    if (!result)
      return result;
    std::vector<std::string> names;
    std::vector<std::size_t> dimensions;
    result = variableShape(outputGroup, observer.name, names, dimensions,
                           outputPath);
    if (!result)
      return result;
    if (names.empty() || names.front() != "t" ||
        (dimensions.size() != 3 && dimensions.size() != 4))
      return failure(WVCheckpointStatusCode::shapeMismatch,
                     "Tracer output has incompatible dimensions.",
                     outputPath + "/" + observer.name);
    std::vector<std::size_t> logical(dimensions.rbegin(),
                                     dimensions.rend() - 1);
    const std::string block = observer.identifier + "-state";
    observer.stateBlockIdentifiers = {block};
    result =
        mergeStateBlock(portable, {block, WVStateScalarType::real64, logical,
                                   WVToleranceKind::uniformAbsolute, tolerance,
                                   WVStateOwnership::integratorOwned,
                                   WVRestartRequirement::requiredDynamicState});
    if (!result)
      return result;
    return WVCheckpointStatus::ok();
    });
  };
  operations.fixedPositions = [&]() {
    return adaptOperation([&]() -> WVCheckpointStatus {
    result = detail::readDoubleScalar(metadataGroup, "outputScale",
                                      observer.outputScale, outputPath);
    if (!result)
      return result;
    result = detail::readDoubleScalar(metadataGroup, "outputOffset",
                                      observer.outputOffset, outputPath);
    if (!result)
      return result;
    std::string interpolation;
    result = optionalTextAttribute(metadataGroup, "trackedVarInterpolation",
                                   interpolation, present, outputPath);
    if (!result)
      return result;
    observer.trackedFieldInterpolation = interpolation == "spline"
                                             ? WVPositionInterpolation::spline
                                             : WVPositionInterpolation::linear;
    result = readWholeDoubleVariable(outputGroup, observer.name + "_x",
                                     observer.x, outputPath);
    if (!result)
      return result;
    result = readWholeDoubleVariable(outputGroup, observer.name + "_y",
                                     observer.y, outputPath);
    if (!result)
      return result;
    result = readWholeDoubleVariable(outputGroup, observer.name + "_z",
                                     observer.z, outputPath);
    if (!result)
      return result;
    return WVCheckpointStatus::ok();
    });
  };
  operations.fixedVerticalProfiles = [&]() {
    return adaptOperation([&]() -> WVCheckpointStatus {
    result = readWholeDoubleVariable(outputGroup, observer.name + "_x",
                                     observer.x, outputPath);
    if (!result)
      return result;
    result = readWholeDoubleVariable(outputGroup, observer.name + "_y",
                                     observer.y, outputPath);
    if (!result)
      return result;
    result = readWholeDoubleVariable(outputGroup, observer.name + "_z",
                                     observer.z, outputPath);
    if (!result)
      return result;
    return WVCheckpointStatus::ok();
    });
  };
  const auto operationStatus =
      compatibility->operationResolver(observer, operations);
  if (!operationStatus && operationResult)
    return failure(WVCheckpointStatusCode::unsupportedObserver,
                   operationStatus.message, outputPath);
  if (!operationResult)
    return operationResult;
  return mergeObserver(portable, std::move(observer));
}


} // namespace

namespace detail {

WVCheckpointStatus parsePersistedObserver(
    int outputGroup, int metadataGroup, const std::string &outputPath,
    const WVExtensionCatalog &catalog,
    WVPortableObserverRecord &portable, std::string &identifier) {
  return parseObserver(outputGroup, metadataGroup, outputPath, catalog, portable,
                       identifier);
}

bool persistedObserverCarriesCoefficientState(
    int metadataGroup, const WVExtensionCatalog &catalog) noexcept {
  std::string className;
  const auto status = readTextAttribute(metadataGroup, "AnnotatedClass",
                                        className, "/observingSystems");
  if (!status)
    return false;
  WVObserverRecord observer;
  observer.typeIdentifier = className;
  const auto *registration = catalog.observers().registration(
      observer.typeIdentifier, observer.contractVersion);
  const auto *compatibility =
      registration == nullptr
          ? nullptr
          : WVObserverFactoryRegistrationAccess::legacyCompatibility(
                *registration);
  return compatibility != nullptr &&
         compatibility->persistence.coefficientRestartFamilies ==
             std::vector<std::string>({"Ap", "Am", "A0"}) &&
         compatibility->persistence.fieldListAttribute.empty();
}

bool legacyObservationAttributeMatches(std::string_view name,
                                       std::string_view expected,
                                       std::string_view observed) noexcept {
  if (expected == observed)
    return true;
  if (name != "long_name")
    return false;
  return (expected == "x position of mooring" &&
          observed == "x coordinate position of mooring") ||
         (expected == "y position of mooring" &&
          observed == "y coordinate position of mooring") ||
         (expected == "x position of particle" &&
          observed ==
              "x coordinate, recorded along the particle trajectory") ||
         (expected == "y position of particle" &&
          observed ==
              "y coordinate, recorded along the particle trajectory") ||
         (expected == "z position of particle" &&
          observed ==
              "z coordinate, recorded along the particle trajectory");
}

WVCheckpointStatus inspectPersistedObserverRestartState(
    const std::vector<WVOutputFileRecord> &files,
    std::vector<WVObserverRecord> &observers, double selectedTime,
    const WVExtensionCatalog &catalog) {
  for (auto &observer : observers) {
    const auto *registration = catalog.observers().registration(
        observer.typeIdentifier, observer.contractVersion);
    const auto *compatibility =
        registration == nullptr
            ? nullptr
            : WVObserverFactoryRegistrationAccess::legacyCompatibility(
                  *registration);
    if (compatibility == nullptr || !compatibility->operationResolver ||
        observer.stateBlockIdentifiers.empty())
      continue;
    std::vector<std::string> stateVariables;
    std::string fixedConfigurationVariable;
    WVLegacyObserverOperationBinder operations;
    const auto noDynamicState = [] { return WVKernelStatus::ok(); };
    operations.fullField = noDynamicState;
    operations.fixedVerticalProfiles = noDynamicState;
    operations.fixedPositions = noDynamicState;
    operations.movingPositions = [&] {
      const auto channels = particlePositionChannels(observer.isXYOnly);
      stateVariables.reserve(channels.size());
      for (const auto channel : channels)
        stateVariables.push_back(
            observer.name + "_" + movingFieldChannelName(channel));
      if (observer.isXYOnly && !observer.z.empty())
        fixedConfigurationVariable = observer.name + "_z";
      return WVKernelStatus::ok();
    };
    operations.integratedState = [&] {
      stateVariables = {observer.name};
      return WVKernelStatus::ok();
    };
    const auto operationStatus =
        compatibility->operationResolver(observer, operations);
    if (!operationStatus)
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     operationStatus.message,
                     "/observingSystems/" + observer.identifier);
    if (stateVariables.empty())
      continue;

    RestartStateLocation reference;
    bool foundReference = false;
    for (const auto &fileRecord : files) {
      for (const auto &groupRecord : fileRecord.groups) {
        if (std::find(groupRecord.observerIdentifiers.begin(),
                      groupRecord.observerIdentifiers.end(),
                      observer.identifier) ==
            groupRecord.observerIdentifiers.end())
          continue;
        WVNetCDFFile file;
        auto result = WVNetCDFFile::openReadOnly(fileRecord.destination, file);
        if (!result)
          return result;
        int group = -1;
        const std::string groupPath = "/" + groupRecord.name;
        result = checkedNetCDF(
            nc_inq_ncid(file.id(), groupRecord.name.c_str(), &group),
            "Observer restart-group lookup", groupPath);
        if (!result)
          return result;
        RestartStateLocation candidate{&fileRecord, &groupRecord, 0};
        bool found = false;
        result = recordAtTime(group, selectedTime, candidate.recordIndex, found,
                              groupPath);
        if (!result)
          return result;
        if (!found)
          continue;
        if (!foundReference) {
          for (const auto &name : stateVariables) {
            result = inspectRestartSlab(candidate, name);
            if (!result)
              return result;
          }
          if (!fixedConfigurationVariable.empty()) {
            result = readRestartConfigurationSlab(
                candidate, fixedConfigurationVariable, observer.z);
            if (!result)
              return result;
          }
          reference = candidate;
          foundReference = true;
          continue;
        }
        for (const auto &name : stateVariables) {
          result = compareRestartSlabs(reference, candidate, name);
          if (!result)
            return result;
        }
        if (!fixedConfigurationVariable.empty()) {
          result = compareRestartSlabs(reference, candidate,
                                       fixedConfigurationVariable);
          if (!result)
            return result;
        }
      }
    }
    if (!foundReference)
      return failure(WVCheckpointStatusCode::missingVariable,
                     "No output group contains required dynamic observer "
                     "state at the selected restart time.",
                     "/observingSystems/" + observer.identifier);
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus resolvePersistedObserverRestartState(
    const std::vector<WVOutputFileRecord> &files,
    std::vector<WVObserverRecord> &observers, double selectedTime,
    const WVExtensionCatalog &catalog,
    std::map<std::string, std::vector<std::vector<double>>> &resolvedState) {
  resolvedState.clear();
  for (auto &observer : observers) {
    const auto *registration = catalog.observers().registration(
        observer.typeIdentifier, observer.contractVersion);
    const auto *compatibility =
        registration == nullptr
            ? nullptr
            : WVObserverFactoryRegistrationAccess::legacyCompatibility(
                  *registration);
    if (compatibility == nullptr || !compatibility->operationResolver ||
        observer.stateBlockIdentifiers.empty())
      continue;
    bool hasDynamicState = false;
    bool hasAdvectedPositions = false;
    WVLegacyObserverOperationBinder operations;
    const auto noDynamicState = [] { return WVKernelStatus::ok(); };
    operations.fullField = noDynamicState;
    operations.fixedVerticalProfiles = noDynamicState;
    operations.fixedPositions = noDynamicState;
    operations.movingPositions = [&] {
      hasDynamicState = true;
      hasAdvectedPositions = true;
      return WVKernelStatus::ok();
    };
    operations.integratedState = [&] {
      hasDynamicState = true;
      return WVKernelStatus::ok();
    };
    const auto operationStatus =
        compatibility->operationResolver(observer, operations);
    if (!operationStatus)
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     operationStatus.message,
                     "/observingSystems/" + observer.identifier);
    if (!hasDynamicState)
      continue;
    std::vector<std::vector<double>> selectedValues;
    std::string selectedLocation;
    bool hasSelectedValues = false;
    for (const auto &file : files) {
      for (const auto &group : file.groups) {
        if (std::find(group.observerIdentifiers.begin(),
                      group.observerIdentifiers.end(), observer.identifier) ==
            group.observerIdentifiers.end())
          continue;
        std::vector<std::vector<double>> values;
        bool found = false;
        auto result = readObserverStateAtTime(
            file.destination, group, observer, *registration, selectedTime,
            values, found);
        if (!result)
          return result;
        if (!found)
          continue;
        const std::string location = file.destination + "/" + group.name;
        if (!hasSelectedValues) {
          selectedValues = std::move(values);
          selectedLocation = location;
          hasSelectedValues = true;
          continue;
        }
        const bool compatible =
            values.size() == selectedValues.size() &&
            std::equal(values.begin(), values.end(), selectedValues.begin(),
                       [](const auto &left, const auto &right) {
                         return sameRestartValues(left, right);
                       });
        if (!compatible)
          return failure(
              WVCheckpointStatusCode::ambiguousState,
              "Dynamic observer state conflicts across output groups.",
              selectedLocation + " and " + location);
      }
    }
    if (!hasSelectedValues)
      return failure(WVCheckpointStatusCode::missingVariable,
                     "No output group contains required dynamic observer "
                     "state at the selected restart time.",
                     "/observingSystems/" + observer.identifier);
    const bool hasFixedParticleZ =
        hasAdvectedPositions && observer.isXYOnly && !observer.z.empty();
    if (selectedValues.size() !=
        observer.stateBlockIdentifiers.size() + (hasFixedParticleZ ? 1 : 0))
      return failure(WVCheckpointStatusCode::shapeMismatch,
                     "Dynamic observer state has an incompatible block "
                     "count.",
                     selectedLocation);
    if (hasAdvectedPositions) {
      observer.x = selectedValues[0];
      observer.y = selectedValues[1];
      if (hasFixedParticleZ)
        observer.z = selectedValues.back();
      else if (selectedValues.size() == 3)
        observer.z = selectedValues[2];
    }
    if (hasFixedParticleZ)
      selectedValues.pop_back();
    resolvedState.emplace(observer.identifier, std::move(selectedValues));
  }
  return WVCheckpointStatus::ok();
}

} // namespace detail
} // namespace wavevortex::runtime
