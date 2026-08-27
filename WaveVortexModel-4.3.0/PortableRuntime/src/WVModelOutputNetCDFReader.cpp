#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include "WVLegacyObservationNetCDFAdapter.hpp"
#include "WVModelOutputTransformAdapters.hpp"
#include "WVNetCDF.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <netcdf.h>
#include <new>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace wavevortex::runtime {
namespace {

bool hexDecode(const std::string &text, std::vector<std::uint8_t> &bytes) {
  if (text.size() % 2 != 0)
    return false;
  auto digit = [](char value) -> int {
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
    const auto high = digit(text[2 * index]);
    const auto low = digit(text[2 * index + 1]);
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

WVCheckpointStatus optionalLogicalAttribute(int group, const char *name,
                                            bool &value, bool &present,
                                            const std::string &path) {
  nc_type type = NC_NAT;
  std::size_t length = 0;
  const int inquiry = nc_inq_att(group, NC_GLOBAL, name, &type, &length);
  if (inquiry == NC_ENOTATT) {
    present = false;
    value = false;
    return WVCheckpointStatus::ok();
  }
  if (inquiry != NC_NOERR)
    return detail::netcdfFailure(inquiry, "Attribute inspection",
                                 path + "/@" + name);
  if ((type != NC_BYTE && type != NC_UBYTE) || length != 1)
    return failure(WVCheckpointStatusCode::typeMismatch,
                   "Logical attributes must contain one byte.",
                   path + "/@" + name);
  unsigned char raw = 0;
  const int read = nc_get_att_uchar(group, NC_GLOBAL, name, &raw);
  if (read != NC_NOERR)
    return detail::netcdfFailure(read, "Attribute read", path + "/@" + name);
  if (raw > 1)
    return failure(WVCheckpointStatusCode::invalidValue,
                   "Logical attributes must contain zero or one.",
                   path + "/@" + name);
  present = true;
  value = raw != 0;
  return WVCheckpointStatus::ok();
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

WVCheckpointStatus optionalVariableTextAttribute(
    int group, int variable, const char *name, std::string &value,
    bool &present, const std::string &path) {
  nc_type type = NC_NAT;
  std::size_t length = 0;
  const int inquiry = nc_inq_att(group, variable, name, &type, &length);
  if (inquiry == NC_ENOTATT) {
    present = false;
    value.clear();
    return WVCheckpointStatus::ok();
  }
  if (inquiry != NC_NOERR)
    return detail::netcdfFailure(inquiry, "Variable-attribute inspection",
                                 path + "/@" + name);
  present = true;
  if (type == NC_CHAR) {
    value.assign(length, '\0');
    return detail::checkedNetCDF(
        nc_get_att_text(group, variable, name, value.data()),
        "Variable-attribute read", path + "/@" + name);
  }
  if (type == NC_STRING && length == 1) {
    char *raw = nullptr;
    const auto result = detail::checkedNetCDF(
        nc_get_att_string(group, variable, name, &raw),
        "Variable-attribute read", path + "/@" + name);
    if (!result)
      return result;
    value = raw == nullptr ? "" : raw;
    if (raw != nullptr)
      nc_free_string(1, &raw);
    return WVCheckpointStatus::ok();
  }
  return failure(WVCheckpointStatusCode::typeMismatch,
                 "Portable observation attributes must contain one string.",
                 path + "/@" + name);
}

WVCheckpointStatus variableStringListAttribute(
    int group, int variable, const char *name, std::vector<std::string> &values,
    bool &present, const std::string &path) {
  nc_type type = NC_NAT;
  std::size_t length = 0;
  const int inquiry = nc_inq_att(group, variable, name, &type, &length);
  if (inquiry == NC_ENOTATT) {
    present = false;
    values.clear();
    return WVCheckpointStatus::ok();
  }
  if (inquiry != NC_NOERR)
    return detail::netcdfFailure(inquiry, "Variable-attribute inspection",
                                 path + "/@" + name);
  if (type != NC_STRING)
    return failure(WVCheckpointStatusCode::typeMismatch,
                   "Observation dimension identifiers must be strings.",
                   path + "/@" + name);
  std::vector<char *> raw(length, nullptr);
  const auto result = detail::checkedNetCDF(
      nc_get_att_string(group, variable, name, raw.data()),
      "Variable-attribute read", path + "/@" + name);
  if (!result)
    return result;
  values.clear();
  values.reserve(length);
  for (const auto *entry : raw)
    values.emplace_back(entry == nullptr ? "" : entry);
  nc_free_string(length, raw.data());
  present = true;
  return WVCheckpointStatus::ok();
}

bool coordinateRole(const std::string &value,
                    WVObservationCoordinateRole &role) noexcept {
  if (value == "none") {
    role = WVObservationCoordinateRole::none;
    return true;
  }
  if (value == "record-time")
    role = WVObservationCoordinateRole::recordTime;
  else if (value == "sample-time")
    role = WVObservationCoordinateRole::sampleTime;
  else if (value == "x")
    role = WVObservationCoordinateRole::x;
  else if (value == "y")
    role = WVObservationCoordinateRole::y;
  else if (value == "z")
    role = WVObservationCoordinateRole::z;
  else if (value == "identifier")
    role = WVObservationCoordinateRole::identifier;
  else if (value == "depth")
    role = WVObservationCoordinateRole::depth;
  else if (value == "pass")
    role = WVObservationCoordinateRole::pass;
  else if (value == "profile")
    role = WVObservationCoordinateRole::profile;
  else
    return false;
  return true;
}

bool valueLayout(const std::string &value,
                 WVObservationValueLayout &layout) noexcept {
  if (value == "static")
    layout = WVObservationValueLayout::staticValue;
  else if (value == "initial")
    layout = WVObservationValueLayout::initialValue;
  else if (value == "record")
    layout = WVObservationValueLayout::record;
  else if (value == "flat")
    layout = WVObservationValueLayout::flat;
  else
    return false;
  return true;
}

bool raggedRole(const std::string &value,
                WVObservationRaggedRole &role) noexcept {
  if (value == "none")
    role = WVObservationRaggedRole::none;
  else if (value == "row-count")
    role = WVObservationRaggedRole::rowCount;
  else if (value == "row-offset")
    role = WVObservationRaggedRole::rowOffset;
  else
    return false;
  return true;
}

bool sameAttributes(const std::vector<WVObservationAttribute> &left,
                    const std::vector<WVObservationAttribute> &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index)
    if (left[index].name != right[index].name ||
        left[index].value != right[index].value)
      return false;
  return true;
}

bool sameObservationSchemaContract(const WVObservationSchema &left,
                                   const WVObservationSchema &right) {
  if (left.identifier != right.identifier || left.version != right.version ||
      left.axes.size() != right.axes.size() ||
      left.variables.size() != right.variables.size())
    return false;
  for (std::size_t index = 0; index < left.axes.size(); ++index) {
    const auto &a = left.axes[index];
    const auto &b = right.axes[index];
    if (a.identifier != b.identifier || a.name != b.name || a.kind != b.kind ||
        a.extent != b.extent || a.coordinateRole != b.coordinateRole)
      return false;
  }
  for (std::size_t index = 0; index < left.variables.size(); ++index) {
    const auto &a = left.variables[index];
    const auto &b = right.variables[index];
    if (a.identifier != b.identifier || a.name != b.name ||
        a.scalarType != b.scalarType ||
        a.dimensionIdentifiers != b.dimensionIdentifiers ||
        a.layout != b.layout || a.units != b.units ||
        a.description != b.description ||
        !sameAttributes(a.attributes, b.attributes) ||
        a.coordinateRole != b.coordinateRole ||
        a.raggedRole != b.raggedRole ||
        a.raggedChildAxisIdentifier != b.raggedChildAxisIdentifier)
      return false;
  }
  return true;
}

bool sameCompleteObservationSchema(const WVObservationSchema &left,
                                   const WVObservationSchema &right) {
  std::vector<std::uint8_t> leftBytes;
  std::vector<std::uint8_t> rightBytes;
  return encodeObservationSchemaManifest(left, leftBytes) &&
         encodeObservationSchemaManifest(right, rightBytes) &&
         leftBytes == rightBytes;
}

WVCheckpointStatus readProvisionalObservationSchema(
    int outputGroup, int metadataGroup, const std::string &observerIdentifier,
    const std::string &path,
    std::vector<WVInspectedObservationSchema> &schemas) {
  WVObservationSchema schema;
  bool present = false;
  auto result = optionalTextAttribute(
      metadataGroup, "portableObservationSchemaIdentifier", schema.identifier,
      present, path);
  if (!result || !present)
    return result;
  unsigned int version = 0;
  result = detail::checkedNetCDF(
      nc_get_att_uint(metadataGroup, NC_GLOBAL,
                      "portableObservationSchemaVersion", &version),
      "Observation-schema version read", path);
  if (!result)
    return result;
  schema.version = version;
  WVObservationSchema declaredSchema;
  std::string encodedManifest;
  bool hasManifest = false;
  result = optionalTextAttribute(
      metadataGroup, "portableObservationSchemaManifest", encodedManifest,
      hasManifest, path);
  if (!result)
    return result;
  if (hasManifest) {
    std::vector<std::uint8_t> manifestBytes;
    if (!hexDecode(encodedManifest, manifestBytes))
      return failure(WVCheckpointStatusCode::invalidValue,
                     "Observation-schema manifest is malformed.", path);
    const auto manifestStatus =
        decodeObservationSchemaManifest(manifestBytes, declaredSchema);
    if (!manifestStatus)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     manifestStatus.message, path);
    if (declaredSchema.identifier != schema.identifier ||
        declaredSchema.version != schema.version)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Observation-schema manifest identity differs.", path);
  }

  int variableCount = 0;
  result = detail::checkedNetCDF(
      nc_inq_varids(outputGroup, &variableCount, nullptr),
      "Observation-variable enumeration", path);
  if (!result)
    return result;
  std::vector<int> variables(static_cast<std::size_t>(variableCount));
  result = detail::checkedNetCDF(
      nc_inq_varids(outputGroup, &variableCount, variables.data()),
      "Observation-variable enumeration", path);
  if (!result)
    return result;
  int unlimitedCount = 0;
  result = detail::checkedNetCDF(
      nc_inq_unlimdims(outputGroup, &unlimitedCount, nullptr),
      "Observation unlimited-axis enumeration", path);
  if (!result)
    return result;
  std::vector<int> unlimited(static_cast<std::size_t>(unlimitedCount));
  if (unlimitedCount > 0) {
    result = detail::checkedNetCDF(
        nc_inq_unlimdims(outputGroup, &unlimitedCount, unlimited.data()),
        "Observation unlimited-axis enumeration", path);
    if (!result)
      return result;
  }

  for (const int variable : variables) {
    char rawName[NC_MAX_NAME + 1] = {};
    result = detail::checkedNetCDF(
        nc_inq_varname(outputGroup, variable, rawName),
        "Observation-variable name inspection", path);
    if (!result)
      return result;
    const std::string variablePath = path + "/../" + rawName;
    std::string observedSchema;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableObservationSchemaIdentifier",
        observedSchema, present, variablePath);
    if (!result)
      return result;
    if (!present || observedSchema != schema.identifier)
      continue;
    std::string observedOwner;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableObservationObserverIdentifier",
        observedOwner, present, variablePath);
    if (!result)
      return result;
    if (!present)
      return failure(WVCheckpointStatusCode::missingAttribute,
                     "Portable observation variable lacks an owning observer.",
                     variablePath);
    if (observedOwner != observerIdentifier)
      continue;
    std::string identifier;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableObservationVariableIdentifier",
        identifier, present, variablePath);
    if (!result || !present)
      return result ? failure(WVCheckpointStatusCode::missingAttribute,
                              "Portable observation variable lacks an identity.",
                              variablePath)
                    : result;
    if (std::any_of(schema.variables.begin(), schema.variables.end(),
                    [&](const auto &candidate) {
                      return candidate.identifier == identifier;
                    }))
      continue;

    WVObservationVariable specification;
    specification.identifier = std::move(identifier);
    specification.name = rawName;
    nc_type type = NC_NAT;
    result = detail::checkedNetCDF(
        nc_inq_vartype(outputGroup, variable, &type),
        "Observation-variable type inspection", variablePath);
    if (!result)
      return result;
    unsigned char complex = 0;
    if (nc_get_att_uchar(outputGroup, variable, "isComplex", &complex) ==
            NC_NOERR &&
        complex != 0) {
      specification.scalarType = WVObservationScalarType::complex64;
      const auto suffix = specification.name.rfind("_real");
      if (suffix == specification.name.size() - 5)
        specification.name.resize(suffix);
    } else if (type == NC_DOUBLE)
      specification.scalarType = WVObservationScalarType::real64;
    else if (type == NC_INT64)
      specification.scalarType = WVObservationScalarType::integer64;
    else if (type == NC_UBYTE)
      specification.scalarType = WVObservationScalarType::boolean8;
    else if (type == NC_STRING)
      specification.scalarType = WVObservationScalarType::text;
    else
      return failure(WVCheckpointStatusCode::typeMismatch,
                     "Portable observation variable has an unsupported type.",
                     variablePath);

    std::string layout;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableObservationValueLayout", layout,
        present, variablePath);
    if (!result || !present)
      return result ? failure(WVCheckpointStatusCode::missingAttribute,
                              "Portable observation variable lacks a layout.",
                              variablePath)
                    : result;
    if (!valueLayout(layout, specification.layout))
      return failure(WVCheckpointStatusCode::invalidValue,
                     "Portable observation variable has an unknown layout.",
                     variablePath);
    result = variableStringListAttribute(
        outputGroup, variable, "portableObservationDimensionIdentifiers",
        specification.dimensionIdentifiers, present, variablePath);
    if (!result)
      return result;
    if (!present)
      specification.dimensionIdentifiers.clear();
    std::vector<std::string> axisRoleNames;
    result = variableStringListAttribute(
        outputGroup, variable, "portableObservationAxisCoordinateRoles",
        axisRoleNames, present, variablePath);
    if (!result)
      return result;
    if (!specification.dimensionIdentifiers.empty() &&
        (!present || axisRoleNames.size() !=
                         specification.dimensionIdentifiers.size()))
      return failure(
          WVCheckpointStatusCode::missingAttribute,
          "Portable observation variable lacks complete axis roles.",
          variablePath);
    std::vector<WVObservationCoordinateRole> axisRoles(axisRoleNames.size());
    for (std::size_t index = 0; index < axisRoleNames.size(); ++index)
      if (!coordinateRole(axisRoleNames[index], axisRoles[index]))
        return failure(
            WVCheckpointStatusCode::invalidValue,
            "Portable observation axis has an unknown coordinate role.",
            variablePath);
    int rank = 0;
    result = detail::checkedNetCDF(
        nc_inq_varndims(outputGroup, variable, &rank),
        "Observation-variable rank inspection", variablePath);
    if (!result)
      return result;
    std::vector<int> dimensions(static_cast<std::size_t>(rank));
    if (rank > 0) {
      result = detail::checkedNetCDF(
          nc_inq_vardimid(outputGroup, variable, dimensions.data()),
          "Observation-variable dimension inspection", variablePath);
      if (!result)
        return result;
    }
    const std::size_t expectedRank =
        specification.dimensionIdentifiers.size() +
        (specification.layout == WVObservationValueLayout::record ? 1 : 0);
    if (expectedRank != dimensions.size())
      return failure(WVCheckpointStatusCode::shapeMismatch,
                     "Portable observation dimensions and layout disagree.",
                     variablePath);

    std::string attribute;
    result = optionalVariableTextAttribute(outputGroup, variable, "units",
                                           specification.units, present,
                                           variablePath);
    if (!result)
      return result;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "long_name", specification.description, present,
        variablePath);
    if (!result)
      return result;
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableCoordinateRole", attribute, present,
        variablePath);
    if (!result)
      return result;
    if (present && !coordinateRole(attribute, specification.coordinateRole))
      return failure(
          WVCheckpointStatusCode::invalidValue,
          "Portable observation variable has an unknown coordinate role.",
          variablePath);
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableRaggedRole", attribute, present,
        variablePath);
    if (!result)
      return result;
    if (present && !raggedRole(attribute, specification.raggedRole))
      return failure(WVCheckpointStatusCode::invalidValue,
                     "Portable observation variable has an unknown ragged role.",
                     variablePath);
    result = optionalVariableTextAttribute(
        outputGroup, variable, "portableRaggedChildAxis",
        specification.raggedChildAxisIdentifier, present, variablePath);
    if (!result)
      return result;

    if (hasManifest) {
      const auto declared = std::find_if(
          declaredSchema.variables.begin(), declaredSchema.variables.end(),
          [&](const auto &candidate) {
            return candidate.identifier == specification.identifier;
          });
      if (declared == declaredSchema.variables.end())
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Portable observation variable is absent from its manifest.",
                       variablePath);
      specification.attributes = declared->attributes;
      for (const auto &custom : specification.attributes) {
        std::string observed;
        bool hasCustom = false;
        result = optionalVariableTextAttribute(
            outputGroup, variable, custom.name.c_str(), observed, hasCustom,
            variablePath);
        if (!result)
          return result;
        if (!hasCustom || observed != custom.value)
          return failure(
              WVCheckpointStatusCode::schemaMismatch,
              "Portable observation custom attributes differ from their manifest.",
              variablePath);
      }
    }

    for (std::size_t logical = 0;
         logical < specification.dimensionIdentifiers.size(); ++logical) {
      const auto &axisIdentifier = specification.dimensionIdentifiers[logical];
      const std::size_t netcdfIndex = dimensions.size() - 1 - logical;
      char axisName[NC_MAX_NAME + 1] = {};
      std::size_t extent = 0;
      result = detail::checkedNetCDF(
          nc_inq_dim(outputGroup, dimensions[netcdfIndex], axisName, &extent),
          "Observation-axis inspection", variablePath);
      if (!result)
        return result;
      const bool isUnlimited =
          std::find(unlimited.begin(), unlimited.end(),
                    dimensions[netcdfIndex]) != unlimited.end();
      const auto existing = std::find_if(
          schema.axes.begin(), schema.axes.end(), [&](const auto &candidate) {
            return candidate.identifier == axisIdentifier;
          });
      const auto role = axisRoles[logical];
      WVObservationAxis observed{
          axisIdentifier, axisName,
          isUnlimited ? WVObservationAxisKind::unlimited
                      : WVObservationAxisKind::fixed,
          isUnlimited ? 0 : extent, role};
      if (existing == schema.axes.end())
        schema.axes.push_back(std::move(observed));
      else if (existing->name != observed.name ||
               existing->kind != observed.kind ||
               existing->extent != observed.extent)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Portable observation axes conflict.", variablePath);
      else if (existing->coordinateRole == WVObservationCoordinateRole::none &&
               role != WVObservationCoordinateRole::none)
        existing->coordinateRole = role;
      else if (existing->coordinateRole != role)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Portable observation axis roles conflict.",
                       variablePath);
    }
    schema.variables.push_back(std::move(specification));
  }
  const auto schemaStatus = validateObservationSchema(schema);
  if (!schemaStatus)
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   schemaStatus.message, path);
  if (hasManifest && !sameObservationSchemaContract(schema, declaredSchema))
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "Persisted observation variables differ from their schema "
                   "manifest.",
                   path);
  if (hasManifest)
    schema = std::move(declaredSchema);
  const auto existing = std::find_if(
      schemas.begin(), schemas.end(), [&](const auto &candidate) {
        return candidate.observerIdentifier == observerIdentifier;
      });
  if (existing == schemas.end())
    schemas.push_back({observerIdentifier, std::move(schema)});
  else if (!sameCompleteObservationSchema(existing->schema, schema))
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "Shared portable observation schemas conflict.", path);
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus
validateReadableHistory(int group, const WVOutputScheduleRecord &schedule,
                        double continuationTime, const std::string &path,
                        WVOutputScheduleCursor &continuation,
                        WVOutputDestinationProgress &destinationProgress) {
  std::size_t timeCount = 0;
  auto result = detail::dimensionLength(group, "t", timeCount, path);
  if (!result)
    return result;
  continuation = {};
  destinationProgress = {};
  destinationProgress.recordCount = timeCount;
  const auto atOrBeforeContinuation = [&](double time) {
    const double tolerance =
        8 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(time), std::abs(continuationTime)});
    return time <= continuationTime + tolerance;
  };
  int timeVariable = -1;
  result = detail::checkedNetCDF(nc_inq_varid(group, "t", &timeVariable),
                                 "Time-variable lookup", path + "/t");
  if (!result)
    return result;
  std::vector<double> times(timeCount);
  if (timeCount > 0) {
    result = detail::checkedNetCDF(
        nc_get_var_double(group, timeVariable, times.data()),
        "Output-time read", path + "/t");
    if (!result)
      return result;
  }
  if (schedule.typeIdentifier.empty()) {
    for (std::size_t index = 0; index < times.size(); ++index) {
      const double raw =
          (times[index] - schedule.initialTime) / schedule.outputInterval;
      const auto ordinal =
          static_cast<WVOutputScheduleOrdinal>(std::llround(raw));
      const double expected =
          schedule.initialTime +
          static_cast<double>(ordinal) * schedule.outputInterval;
      const double tolerance =
          8 * std::numeric_limits<double>::epsilon() *
          std::max({1.0, std::abs(times[index]), std::abs(expected)});
      if (!std::isfinite(times[index]) || ordinal < 0 ||
          std::abs(times[index] - expected) > tolerance ||
          (index > 0 && !(times[index] > times[index - 1])))
        return failure(WVCheckpointStatusCode::appendConflict,
                       "Output history is not a strictly increasing subset of "
                       "its configured schedule lattice.",
                       path + "/t");
      destinationProgress.committedScheduleCursor.committedOrdinal = ordinal;
      if (atOrBeforeContinuation(times[index]))
        continuation.committedOrdinal = ordinal;
    }
  } else {
    int ordinalVariable = -1, cursorVariable = -1;
    result = detail::checkedNetCDF(
        nc_inq_varid(group, "portableScheduleOrdinal", &ordinalVariable),
        "Schedule-ordinal lookup", path);
    if (!result)
      return result;
    result = detail::checkedNetCDF(
        nc_inq_varid(group, "portableScheduleCursor", &cursorVariable),
        "Schedule-cursor lookup", path);
    if (!result)
      return result;
    std::vector<long long> ordinals(timeCount);
    result = detail::checkedNetCDF(
        nc_get_var_longlong(group, ordinalVariable, ordinals.data()),
        "Schedule-ordinal read", path);
    if (!result)
      return result;
    for (std::size_t index = 0; index < timeCount; ++index) {
      if (!std::isfinite(times[index]) || ordinals[index] < 0 ||
          (index > 0 && (!(times[index] > times[index - 1]) ||
                         ordinals[index] <= ordinals[index - 1])))
        return failure(WVCheckpointStatusCode::appendConflict,
                       "Algorithmic output history is not strictly monotone.",
                       path + "/t");
      const std::size_t position[] = {index};
      char *encodedCursor = nullptr;
      result = detail::checkedNetCDF(
          nc_get_var1_string(group, cursorVariable, position, &encodedCursor),
          "Schedule-cursor read", path);
      if (!result)
        return result;
      const std::string encoded = encodedCursor == nullptr ? "" : encodedCursor;
      if (encodedCursor != nullptr)
        nc_free_string(1, &encodedCursor);
      std::vector<std::uint8_t> cursorBytes;
      if (!hexDecode(encoded, cursorBytes) ||
          cursorBytes.size() > WVMaximumOutputScheduleCursorBytes)
        return failure(WVCheckpointStatusCode::appendConflict,
                       "Algorithmic output cursor is malformed or oversized.",
                       path);
      WVPortableTypedRecord cursor;
      if (!cursorBytes.empty()) {
        const auto decoded = decodePortableTypedRecord(
            cursorBytes, cursor,
            {WVMaximumOutputScheduleCursorBytes, true, false});
        if (!decoded)
          return failure(WVCheckpointStatusCode::appendConflict,
                         decoded.message, path);
      }
      WVOutputScheduleCursor current{
          static_cast<WVOutputScheduleOrdinal>(ordinals[index]),
          std::move(cursor)};
      destinationProgress.committedScheduleCursor = current;
      if (atOrBeforeContinuation(times[index]))
        continuation = std::move(current);
    }
  }
  if (!times.empty()) {
    destinationProgress.hasCommittedTime = true;
    destinationProgress.lastCommittedTime = times.back();
  }
  int unlimitedCount = 0;
  result = detail::checkedNetCDF(
      nc_inq_unlimdims(group, &unlimitedCount, nullptr),
      "Unlimited-axis enumeration", path);
  if (!result)
    return result;
  std::vector<int> unlimitedDimensions(
      static_cast<std::size_t>(unlimitedCount));
  if (unlimitedCount > 0) {
    result = detail::checkedNetCDF(
        nc_inq_unlimdims(group, &unlimitedCount, unlimitedDimensions.data()),
        "Unlimited-axis enumeration", path);
    if (!result)
      return result;
  }
  int variableCount = 0;
  result = detail::checkedNetCDF(nc_inq_varids(group, &variableCount, nullptr),
                                 "Output-variable enumeration", path);
  if (!result)
    return result;
  std::vector<int> variables(static_cast<std::size_t>(variableCount));
  result = detail::checkedNetCDF(
      nc_inq_varids(group, &variableCount, variables.data()),
      "Output-variable enumeration", path);
  if (!result)
    return result;
  for (const int variable : variables) {
    nc_type type = NC_NAT;
    int dimensionCount = 0;
    if (nc_inq_vartype(group, variable, &type) != NC_NOERR ||
        nc_inq_varndims(group, variable, &dimensionCount) != NC_NOERR)
      return failure(WVCheckpointStatusCode::netcdfFailure,
                     "Unable to inspect output variable.", path);
    if ((type != NC_DOUBLE && type != NC_INT64 && type != NC_UBYTE &&
         type != NC_STRING) ||
        dimensionCount == 0)
      continue;
    std::vector<int> dimensions(static_cast<std::size_t>(dimensionCount));
    result = detail::checkedNetCDF(
        nc_inq_vardimid(group, variable, dimensions.data()),
        "Output-variable dimension inspection", path);
    if (!result)
      return result;
    char firstName[NC_MAX_NAME + 1] = {};
    result = detail::checkedNetCDF(
        nc_inq_dimname(group, dimensions.front(), firstName),
        "Output-variable dimension-name inspection", path);
    if (!result)
      return result;
    const bool timeLeading = std::string(firstName) == "t";
    const bool hasUnlimitedDimension = std::any_of(
        dimensions.begin(), dimensions.end(), [&](const int dimension) {
          return std::find(unlimitedDimensions.begin(),
                           unlimitedDimensions.end(),
                           dimension) != unlimitedDimensions.end();
        });
    if (!timeLeading && !hasUnlimitedDimension)
      continue;
    std::vector<std::size_t> lengths(dimensions.size(), 1);
    const std::size_t firstDataDimension = timeLeading ? 1 : 0;
    for (std::size_t dimension = firstDataDimension;
         dimension < dimensions.size(); ++dimension) {
      result = detail::checkedNetCDF(
          nc_inq_dimlen(group, dimensions[dimension], &lengths[dimension]),
          "Output-variable dimension-length inspection", path);
      if (!result)
        return result;
    }
    if (std::find(lengths.begin() + firstDataDimension, lengths.end(), 0) !=
        lengths.end())
      continue;

    constexpr std::size_t maximumValidationChunkElements = 4096;
    std::vector<double> realValues(maximumValidationChunkElements);
    std::vector<long long> integerValues(maximumValidationChunkElements);
    std::vector<unsigned char> booleanValues(maximumValidationChunkElements);
    std::vector<char *> textValues(maximumValidationChunkElements, nullptr);
    const auto readChunk = [&](const std::size_t *start,
                               const std::size_t *count,
                               std::size_t elementCount,
                               bool &containsFill) -> WVCheckpointStatus {
      containsFill = false;
      int readStatus = NC_NOERR;
      if (type == NC_DOUBLE) {
        readStatus = nc_get_vara_double(group, variable, start, count,
                                        realValues.data());
        containsFill =
            std::find(realValues.begin(), realValues.begin() + elementCount,
                      NC_FILL_DOUBLE) != realValues.begin() + elementCount;
      } else if (type == NC_INT64) {
        readStatus = nc_get_vara_longlong(group, variable, start, count,
                                          integerValues.data());
        containsFill =
            std::find(integerValues.begin(),
                      integerValues.begin() + elementCount, NC_FILL_INT64) !=
            integerValues.begin() + elementCount;
      } else if (type == NC_UBYTE) {
        readStatus = nc_get_vara_uchar(group, variable, start, count,
                                       booleanValues.data());
        containsFill =
            std::find(booleanValues.begin(),
                      booleanValues.begin() + elementCount, NC_FILL_UBYTE) !=
            booleanValues.begin() + elementCount;
      } else {
        std::fill(textValues.begin(), textValues.begin() + elementCount,
                  nullptr);
        readStatus = nc_get_vara_string(group, variable, start, count,
                                        textValues.data());
        if (readStatus == NC_NOERR) {
          containsFill = std::any_of(
              textValues.begin(), textValues.begin() + elementCount,
              [](const char *value) { return value == nullptr; });
          nc_free_string(elementCount, textValues.data());
        }
      }
      return detail::checkedNetCDF(readStatus, "Committed output-record read",
                                   path);
    };

    if (timeLeading && dimensions.size() == 1) {
      const std::size_t count[] = {1};
      for (std::size_t record = 0; record < timeCount; ++record) {
        const std::size_t start[] = {record};
        bool containsFill = false;
        result = readChunk(start, count, 1, containsFill);
        if (!result)
          return result;
        if (containsFill) {
          char name[NC_MAX_NAME + 1] = {};
          nc_inq_varname(group, variable, name);
          return failure(
              WVCheckpointStatusCode::incompleteRecord,
              "A committed output record contains unwritten payload values.",
              path + "/" + name);
        }
      }
      continue;
    }

    std::size_t chunkDimension = dimensions.size() - 1;
    std::size_t trailingElements = 1;
    while (chunkDimension > firstDataDimension &&
           lengths[chunkDimension] <=
               maximumValidationChunkElements / trailingElements) {
      trailingElements *= lengths[chunkDimension];
      --chunkDimension;
    }
    const auto chunkExtent = std::max<std::size_t>(
        1, std::min(lengths[chunkDimension],
                    maximumValidationChunkElements / trailingElements));
    std::vector<std::size_t> start(dimensions.size(), 0);
    std::vector<std::size_t> count(dimensions.size(), 1);
    for (std::size_t dimension = chunkDimension + 1;
         dimension < dimensions.size(); ++dimension)
      count[dimension] = lengths[dimension];

    std::size_t prefixCount = 1;
    for (std::size_t dimension = firstDataDimension;
         dimension < chunkDimension; ++dimension) {
      if (lengths[dimension] != 0 &&
          prefixCount >
              std::numeric_limits<std::size_t>::max() / lengths[dimension])
        return failure(WVCheckpointStatusCode::shapeMismatch,
                       "Output-variable dimensions are oversized.", path);
      prefixCount *= lengths[dimension];
    }
    const auto outerRecordCount = timeLeading ? timeCount : 1;
    for (std::size_t record = 0; record < outerRecordCount; ++record) {
      if (timeLeading)
        start.front() = record;
      for (std::size_t prefix = 0; prefix < prefixCount; ++prefix) {
        std::size_t remainder = prefix;
        for (std::size_t dimension = chunkDimension;
             dimension-- > firstDataDimension;) {
          start[dimension] = remainder % lengths[dimension];
          remainder /= lengths[dimension];
        }
        for (std::size_t offset = 0; offset < lengths[chunkDimension];
             offset += chunkExtent) {
          start[chunkDimension] = offset;
          count[chunkDimension] =
              std::min(chunkExtent, lengths[chunkDimension] - offset);
          const auto elementCount = count[chunkDimension] * trailingElements;
          bool containsFill = false;
          result = readChunk(start.data(), count.data(), elementCount,
                             containsFill);
          if (!result)
            return result;
          if (containsFill) {
            char name[NC_MAX_NAME + 1] = {};
            nc_inq_varname(group, variable, name);
            return failure(
                WVCheckpointStatusCode::incompleteRecord,
                "A committed output record contains unwritten payload values.",
                path + "/" + name);
          }
        }
      }
    }
  }

  for (const int dimension : unlimitedDimensions) {
    char rawName[NC_MAX_NAME + 1] = {};
    std::size_t physicalCount = 0;
    result = detail::checkedNetCDF(
        nc_inq_dim(group, dimension, rawName, &physicalCount),
        "Unlimited-axis inspection", path);
    if (!result)
      return result;
    const std::string name = rawName;
    if (name == "t")
      continue;
    WVOutputDestinationAxisProgress axis;
    axis.axisIdentifier = name;
    axis.physicalCount = physicalCount;
    if (timeCount > 0) {
      int progressVariable = -1;
      result = detail::checkedNetCDF(
          nc_inq_varid(group, ("portableCommitted_" + name).c_str(),
                       &progressVariable),
          "Unlimited-axis progress lookup", path + "/" + name);
      if (!result)
        return result;
      const std::size_t position[] = {timeCount - 1};
      long long committedCount = -1;
      result = detail::checkedNetCDF(
          nc_get_var1_longlong(group, progressVariable, position,
                               &committedCount),
          "Unlimited-axis progress read", path + "/" + name);
      if (!result)
        return result;
      if (committedCount < 0 ||
          static_cast<std::uint64_t>(committedCount) != physicalCount)
        return failure(WVCheckpointStatusCode::incompleteRecord,
                       "Unlimited-axis physical and committed offsets differ.",
                       path + "/" + name);
      axis.committedCount = static_cast<std::size_t>(committedCount);
    } else if (physicalCount != 0) {
      return failure(WVCheckpointStatusCode::incompleteRecord,
                     "An empty output group has committed ragged storage.",
                     path + "/" + name);
    }
    destinationProgress.unlimitedAxes.push_back(std::move(axis));
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus parseOutputFile(const std::string &path,
                                   const WVExtensionCatalog &catalog,
                                   double continuationTime,
                                   const WVTransformStateDescription
                                       &stateDescription,
                                   WVPortableObserverRecord &portable,
                                   std::vector<WVInspectedObservationSchema>
                                       &observationSchemas,
                                   std::vector<WVOutputScheduleContinuation>
                                       &continuations,
                                   std::vector<WVOutputDestinationProgress>
                                       &destinationProgress,
                                   WVOutputFileRecord &fileRecord,
                                   bool &isDynamicsLinear) {
  detail::WVNetCDFFile file;
  auto result = detail::WVNetCDFFile::openReadOnly(path, file);
  if (!result)
    return result;
  bool present = false;
  result = optionalLogicalAttribute(file.id(), "WVModelIsDynamicsLinear",
                                    isDynamicsLinear, present, "/");
  if (!result)
    return result;
  if (!present)
    isDynamicsLinear = false;
  result = optionalTextAttribute(file.id(), "portableFileIdentifier",
                                 fileRecord.identifier, present, "/");
  if (!result)
    return result;
  if (!present)
    fileRecord.identifier = portableIdentifier(
        std::filesystem::absolute(path).lexically_normal().string());
  fileRecord.destination = path;

  std::vector<int> children;
  result = detail::childGroups(file.id(), children, "/");
  if (!result)
    return result;
  std::size_t restartGroupCount = 0;
  for (const int groupId : children) {
    std::string groupName;
    result = detail::groupName(groupId, groupName, "/");
    if (!result)
      return result;
    const std::string groupPath = "/" + groupName;
    std::string className;
    result = optionalTextAttribute(groupId, "AnnotatedClass", className,
                                   present, groupPath);
    if (!result)
      return result;
    if (!present || className != "WVModelOutputGroupEvenlySpaced")
      continue;
    WVOutputGroupRecord group;
    group.name = groupName;
    result = optionalTextAttribute(groupId, "portableIdentifier",
                                   group.identifier, present, groupPath);
    if (!result)
      return result;
    if (!present)
      group.identifier = portableIdentifier(groupName);
    result = optionalTextAttribute(groupId, "portableScheduleTypeIdentifier",
                                   group.schedule.typeIdentifier, present,
                                   groupPath);
    if (!result)
      return result;
    if (present) {
      unsigned int version = 0;
      result = detail::checkedNetCDF(
          nc_get_att_uint(groupId, NC_GLOBAL, "portableScheduleContractVersion",
                          &version),
          "Schedule-version attribute read", groupPath);
      if (!result)
        return result;
      group.schedule.contractVersion = version;
      std::string encodedConfiguration;
      bool hasConfiguration = false;
      result = optionalTextAttribute(groupId, "portableScheduleConfiguration",
                                     encodedConfiguration, hasConfiguration,
                                     groupPath);
      if (!result)
        return result;
      std::vector<std::uint8_t> configurationBytes;
      if (!hasConfiguration ||
          !hexDecode(encodedConfiguration, configurationBytes))
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Algorithmic output schedule configuration is missing "
                       "or malformed.",
                       groupPath);
      const auto decoded = decodePortableTypedRecord(
          configurationBytes, group.schedule.configuration);
      if (!decoded)
        return failure(WVCheckpointStatusCode::invalidValue, decoded.message,
                       groupPath);
    }
    result = detail::readDoubleScalar(groupId, "outputInterval",
                                      group.schedule.outputInterval, groupPath);
    if (!result)
      return result;
    result = detail::readDoubleScalar(groupId, "initialTime",
                                      group.schedule.initialTime, groupPath);
    if (!result)
      return result;
    result = detail::readDoubleScalar(groupId, "finalTime",
                                      group.schedule.finalTime, groupPath);
    if (!result)
      return result;
    bool hasCoefficients = !stateDescription.coefficientFamilies.empty();
    for (const auto &family : stateDescription.coefficientFamilies) {
      int coefficientVariable = -1;
      bool hasFamily = false;
      result = detail::variableIdIfPresent(
          groupId, family.identifier + "_real", coefficientVariable,
          hasFamily, groupPath);
      if (!result)
        return result;
      if (!hasFamily) {
        result = detail::variableIdIfPresent(
            groupId, family.identifier, coefficientVariable, hasFamily,
            groupPath);
        if (!result)
          return result;
      }
      hasCoefficients = hasCoefficients && hasFamily;
    }
    int metadataRoot = -1;
    const int metadataLookup =
        nc_inq_ncid(groupId, "observingSystems", &metadataRoot);
    if (metadataLookup == NC_ENOGRP)
      return failure(WVCheckpointStatusCode::missingVariable,
                     "Output group lacks observing-system metadata.",
                     groupPath + "/observingSystems");
    result =
        detail::checkedNetCDF(metadataLookup, "Observer metadata-root lookup",
                              groupPath + "/observingSystems");
    if (!result)
      return result;
    std::vector<int> observers;
    result = detail::childGroups(metadataRoot, observers,
                                 groupPath + "/observingSystems");
    if (!result)
      return result;
    std::string rootObserverClass;
    result =
        optionalTextAttribute(metadataRoot, "AnnotatedClass", rootObserverClass,
                              present, groupPath + "/observingSystems");
    if (!result)
      return result;
    if (present) {
      std::string observerIdentifier;
      result = detail::parsePersistedObserver(
          groupId, metadataRoot, groupPath + "/observingSystems", catalog,
          portable,
          observerIdentifier);
      if (!result)
        return result;
      result = readProvisionalObservationSchema(
          groupId, metadataRoot, observerIdentifier,
          groupPath + "/observingSystems", observationSchemas);
      if (!result)
        return result;
      group.observerIdentifiers.push_back(std::move(observerIdentifier));
    }
    std::map<int, std::string> parsedIdentifiers;
    const auto parseOne = [&](int observerGroup) {
      std::string observerIdentifier;
      auto observerResult = detail::parsePersistedObserver(
          groupId, observerGroup, groupPath + "/observingSystems", catalog,
          portable,
          observerIdentifier);
      if (observerResult)
        observerResult = readProvisionalObservationSchema(
            groupId, observerGroup, observerIdentifier,
            groupPath + "/observingSystems", observationSchemas);
      if (observerResult)
        parsedIdentifiers.emplace(observerGroup, std::move(observerIdentifier));
      return observerResult;
    };
    std::set<int> parsedObservers;
    bool hasCoefficientObserver =
        detail::persistedObserverCarriesCoefficientState(metadataRoot,
                                                         catalog);
    for (const int observerGroup : observers) {
      std::string observerClass;
      result =
          optionalTextAttribute(observerGroup, "AnnotatedClass", observerClass,
                                present, groupPath + "/observingSystems");
      if (!result)
        return result;
      if (present &&
          detail::persistedObserverCarriesCoefficientState(observerGroup,
                                                           catalog)) {
        hasCoefficientObserver = true;
        result = parseOne(observerGroup);
        if (!result)
          return result;
        parsedObservers.insert(observerGroup);
      }
    }
    for (const int observerGroup : observers) {
      if (parsedObservers.find(observerGroup) != parsedObservers.end())
        continue;
      result = parseOne(observerGroup);
      if (!result)
        return result;
    }
    for (const int observerGroup : observers)
      group.observerIdentifiers.push_back(parsedIdentifiers.at(observerGroup));

    group.containsCompleteCoefficientRestart =
        detail::modelOutputGroupCarriesCompleteCoefficientRestart(
            stateDescription, hasCoefficients, hasCoefficientObserver);
    restartGroupCount +=
        group.containsCompleteCoefficientRestart ? 1 : 0;

    WVOutputScheduleCursor continuation;
    WVOutputDestinationProgress committed;
    result = validateReadableHistory(groupId, group.schedule, continuationTime,
                                     groupPath, continuation, committed);
    if (!result)
      return result;
    continuations.push_back(
        {fileRecord.identifier, group.identifier, std::move(continuation)});
    committed.fileIdentifier = fileRecord.identifier;
    committed.groupIdentifier = group.identifier;
    destinationProgress.push_back(std::move(committed));
    fileRecord.groups.push_back(std::move(group));
  }
  if (fileRecord.groups.empty())
    return failure(WVCheckpointStatusCode::missingVariable,
                   "File contains no supported output group.", path);
  if (restartGroupCount != 1)
    return failure(WVCheckpointStatusCode::ambiguousState,
                   "Each restartable output file must contain exactly one "
                   "complete coefficient stream.",
                   path);
  return WVCheckpointStatus::ok();
}

} // namespace

WVCheckpointStatus
WVModelOutputNetCDFSink::inspect(const std::vector<std::string> &paths,
                                 const WVExtensionCatalog &catalog,
                                 WVModelOutputNetCDFInspection &inspection) {
  if (paths.empty())
    return failure(WVCheckpointStatusCode::openFailure,
                   "At least one output path is required.", "/");
  try {
    WVModelOutputNetCDFInspection candidate;
    candidate.paths = paths;
    std::string selectedPath;
    WVCheckpointInspection selectedCheckpoint;

    // Select and cross-check the restart using allocation-light checkpoint
    // inspection before reconstructing any output graph state.
    for (const auto &path : paths) {
      WVCheckpointInspection checkpoint;
      auto result = WVCheckpointReader::inspect(path, catalog, checkpoint);
      if (!result)
        return result;
      if (selectedPath.empty()) {
        selectedCheckpoint = std::move(checkpoint);
        selectedPath = path;
      } else {
        if (!detail::sameModelOutputTransformConfiguration(
                selectedCheckpoint, checkpoint))
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Output files do not describe one model.", path);
        if (checkpoint.t > selectedCheckpoint.t) {
          selectedCheckpoint = std::move(checkpoint);
          selectedPath = path;
        }
      }
    }

    bool firstFile = true;
    for (const auto &path : paths) {
      WVOutputFileRecord fileRecord;
      bool isDynamicsLinear = false;
      auto result = parseOutputFile(
          path, catalog, selectedCheckpoint.t,
          selectedCheckpoint.stateDescription, candidate.observerRecord,
          candidate.observationSchemas, candidate.scheduleContinuations,
          candidate.destinationProgress, fileRecord, isDynamicsLinear);
      if (!result)
        return result;
      if (firstFile) {
        candidate.isDynamicsLinear = isDynamicsLinear;
        firstFile = false;
      } else if (candidate.isDynamicsLinear != isDynamicsLinear) {
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Output files disagree on the model dynamics mode.",
                       path);
      }
      candidate.observerRecord.outputFiles.push_back(std::move(fileRecord));
    }

    auto selectedFile = std::find_if(
        candidate.observerRecord.outputFiles.begin(),
        candidate.observerRecord.outputFiles.end(), [&](const auto &file) {
          return std::filesystem::absolute(file.destination)
                     .lexically_normal() ==
                 std::filesystem::absolute(selectedPath).lexically_normal();
        });
    if (selectedFile == candidate.observerRecord.outputFiles.end())
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Selected restart file is absent from the output graph.",
                     selectedPath);
    const std::string selectedGroupName =
        selectedCheckpoint.metadata.stateGroupPath.size() > 1
            ? selectedCheckpoint.metadata.stateGroupPath.substr(1)
            : std::string{};
    auto selectedGroup = std::find_if(
        selectedFile->groups.begin(), selectedFile->groups.end(),
        [&](const auto &group) { return group.name == selectedGroupName; });
    if (selectedGroup == selectedFile->groups.end())
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Selected restart group is absent from the output graph.",
                     selectedCheckpoint.metadata.stateGroupPath);
    // Validate legacy dynamic observer state with bounded scratch storage.
    // Full state slabs are materialized only by restoreState(), after the
    // canonical graph and every destination have passed preflight.
    auto legacyStateStatus = detail::inspectPersistedObserverRestartState(
        candidate.observerRecord.outputFiles,
        candidate.observerRecord.observers, selectedCheckpoint.t, catalog);
    if (!legacyStateStatus)
      return legacyStateStatus;
    candidate.latestRestart = std::move(selectedCheckpoint);
    candidate.latestRestartPath = std::move(selectedPath);
    inspection = std::move(candidate);
    return WVCheckpointStatus::ok();
  } catch (const std::bad_alloc &) {
    return failure(WVCheckpointStatusCode::writeFailure,
                   "Output inspection allocation failed.", "/");
  }
}

WVCheckpointStatus WVModelOutputNetCDFSink::restoreState(
    const WVModelOutputNetCDFInspection &inspection,
    const WVExtensionCatalog &catalog,
    const WVIntegrationStateLayout &stateLayout, WVCheckpoint &checkpoint,
    WVAdditionalStateStorage &additionalState) {
  if (inspection.latestRestartPath.empty())
    return failure(WVCheckpointStatusCode::openFailure,
                   "Output inspection has no selected restart path.", "/");
  try {
    WVCheckpoint restoredCheckpoint;
    auto result = WVCheckpointReader::read(
        inspection.latestRestartPath, catalog, restoredCheckpoint,
        WVCheckpointStateSelection::atIndex(
            inspection.latestRestart.metadata.selectedStateIndex));
    if (!result)
      return result;
    if (restoredCheckpoint.transformKind !=
            inspection.latestRestart.transformKind ||
        restoredCheckpoint.state.t != inspection.latestRestart.t ||
        restoredCheckpoint.state.t0 != inspection.latestRestart.t0 ||
        restoredCheckpoint.stateDescription.transformIdentifier !=
            inspection.latestRestart.stateDescription.transformIdentifier ||
        restoredCheckpoint.stateDescription.coefficientFamilies.size() !=
            inspection.latestRestart.stateDescription.coefficientFamilies
                .size())
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Selected restart changed after graph inspection.",
                     inspection.latestRestartPath);

    std::map<std::string, std::vector<std::vector<double>>> resolvedState;
    auto restoredObservers = inspection.observerRecord.observers;
    auto legacyStateStatus = detail::resolvePersistedObserverRestartState(
        inspection.observerRecord.outputFiles,
        restoredObservers, inspection.latestRestart.t, catalog, resolvedState);
    if (!legacyStateStatus)
      return legacyStateStatus;

    WVAdditionalStateStorage restoredAdditionalState;
    const auto storageStatus = restoredAdditionalState.initialize(stateLayout);
    if (!storageStatus)
      return failure(WVCheckpointStatusCode::descriptorFailure,
                     storageStatus.message, "/observingSystems");
    for (const auto &observer : restoredObservers) {
      if (observer.stateBlockIdentifiers.empty())
        continue;
      const auto values = resolvedState.find(observer.identifier);
      if (values == resolvedState.end())
        continue;
      for (std::size_t valueIndex = 0;
           valueIndex < observer.stateBlockIdentifiers.size(); ++valueIndex) {
        const auto block =
            std::find_if(restoredAdditionalState.mutableBlocks(),
                         restoredAdditionalState.mutableBlocks() +
                             restoredAdditionalState.blockCount(),
                         [&](const auto &view) {
                           return view.layout->identifier ==
                                  observer.stateBlockIdentifiers[valueIndex];
                         });
        if (block == restoredAdditionalState.mutableBlocks() +
                         restoredAdditionalState.blockCount())
          return failure(WVCheckpointStatusCode::descriptorFailure,
                         "Resolved dynamic observer state is absent from the "
                         "integration layout.",
                         "/observingSystems/" + observer.identifier);
        if (block->layout->elementCount != values->second[valueIndex].size())
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "Dynamic observer restart state has an incompatible "
                         "shape.",
                         "/observingSystems/" + observer.identifier);
        std::copy(values->second[valueIndex].begin(),
                  values->second[valueIndex].end(), block->realData);
      }
    }
    checkpoint = std::move(restoredCheckpoint);
    additionalState = std::move(restoredAdditionalState);
    return WVCheckpointStatus::ok();
  } catch (const std::bad_alloc &) {
    return failure(WVCheckpointStatusCode::writeFailure,
                   "Output restart-state allocation failed.", "/");
  }
}

} // namespace wavevortex::runtime
