#include "WVModelOutputTransformAdapters.hpp"

#include "WVNetCDF.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <utility>
#include <variant>

namespace wavevortex::runtime::detail {
namespace {

WVCheckpointStatus failed(WVCheckpointStatusCode code, std::string message,
                          std::string location) {
  return {code, std::move(message), std::move(location)};
}

WVCheckpointStatus defineLogical(int group, const std::string &name,
                                 int &variable, const std::string &path) {
  auto result = checkedNetCDF(
      nc_def_var(group, name.c_str(), NC_UBYTE, 0, nullptr, &variable),
      "Variable definition", path + "/" + name);
  if (!result)
    return result;
  return putByteAttribute(group, variable, "isLogicalType", 1,
                          path + "/" + name);
}

WVCheckpointStatus defineForcingEntry(int group,
                                      const WVFrozenForcingEntry &entry,
                                      const WVForcingCatalog &catalog,
                                      const std::array<int, 5> &rootDimensions,
                                      const std::string &path) {
  auto result = putTextAttribute(group, NC_GLOBAL, "AnnotatedClass",
                                 entry.typeIdentifier, path);
  if (!result)
    return result;
  const auto *registration =
      catalog.registration(entry.typeIdentifier, entry.contractVersion);
  if (registration == nullptr || !registration->isSupported)
    return failed(WVCheckpointStatusCode::unsupportedForcing,
                  "Unsupported forcing reached output definition.", path);
  if (registration->persistence.writesNameAttribute) {
    result = putTextAttribute(group, NC_GLOBAL, "name", entry.name, path);
    if (!result)
      return result;
  }
  for (const auto &field : registration->persistence.fields) {
    const auto *value = entry.configuration.value(field.recordName);
    if (value == nullptr && field.optional)
      continue;
    if (value == nullptr)
      return failed(WVCheckpointStatusCode::malformedForcing,
                    "Required forcing configuration value is missing.", path);
    if (field.encoding == WVForcingPersistenceEncoding::textAttribute) {
      const auto *values =
          std::get_if<std::vector<std::string>>(&value->storage);
      if (values == nullptr || values->size() != 1)
        return failed(WVCheckpointStatusCode::malformedForcing,
                      "Forcing text attribute is malformed.", path);
      if (!values->front().empty()) {
        result = putTextAttribute(group, NC_GLOBAL, field.netcdfName.c_str(),
                                  values->front(), path);
        if (!result)
          return result;
      }
      continue;
    }
    std::vector<int> dimensions;
    if (field.dimensions == WVForcingDimensionRule::ownLength) {
      if (value->valueCount() == 0)
        continue;
      int dimension = -1;
      result = checkedNetCDF(
          nc_def_dim(group, field.netcdfName.c_str(), value->valueCount(),
                     &dimension),
          "Forcing dimension definition", path + "/" + field.netcdfName);
      if (!result)
        return result;
      dimensions = {dimension};
    } else if (field.dimensions ==
               WVForcingDimensionRule::referencedLength) {
      int dimension = -1;
      const int inquiry =
          nc_inq_dimid(group, field.dimensionReference.c_str(), &dimension);
      if (inquiry == NC_EBADDIM && field.optional)
        continue;
      result = checkedNetCDF(inquiry, "Forcing dimension lookup",
                             path + "/" + field.dimensionReference);
      if (!result)
        return result;
      dimensions = {dimension};
    } else if (field.dimensions == WVForcingDimensionRule::horizontalYX) {
      dimensions = {rootDimensions[3], rootDimensions[2]};
    } else if (field.dimensions == WVForcingDimensionRule::componentPair) {
      int dimension = -1;
      result = checkedNetCDF(
          nc_def_dim(group, "barotropicVelocityComponent", 2, &dimension),
          "Forcing dimension definition",
          path + "/barotropicVelocityComponent");
      if (!result)
        return result;
      int coordinate = -1;
      result = defineDoubleVariable(group, "barotropicVelocityComponent",
                                    {dimension}, coordinate, path);
      if (!result)
        return result;
      dimensions = {dimension};
    }
    int variable = -1;
    if (field.encoding == WVForcingPersistenceEncoding::realVariable)
      result = defineDoubleVariable(group, field.netcdfName, dimensions,
                                    variable, path);
    else if (field.encoding == WVForcingPersistenceEncoding::logicalVariable)
      result = defineLogical(group, field.netcdfName, variable, path);
    else if (field.encoding ==
             WVForcingPersistenceEncoding::zeroBasedIndexVariable)
      result = checkedNetCDF(
          nc_def_var(group, field.netcdfName.c_str(), NC_UINT64,
                     static_cast<int>(dimensions.size()), dimensions.data(),
                     &variable),
          "Forcing variable definition", path + "/" + field.netcdfName);
    else if (field.encoding == WVForcingPersistenceEncoding::complexVariable)
      result = defineComplexVariable(group, field.netcdfName, dimensions, path);
    if (!result)
      return result;
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus variableId(int group, const std::string &name, int &variable,
                              const std::string &path) {
  return checkedNetCDF(nc_inq_varid(group, name.c_str(), &variable),
                       "Variable lookup", path + "/" + name);
}

WVCheckpointStatus writeDouble(int group, const std::string &name, double value,
                               const std::string &path) {
  int variable = -1;
  auto result = variableId(group, name, variable, path);
  if (!result)
    return result;
  return checkedNetCDF(nc_put_var_double(group, variable, &value),
                       "Variable write", path + "/" + name);
}

WVCheckpointStatus writeLogical(int group, const std::string &name, bool value,
                                const std::string &path) {
  int variable = -1;
  auto result = variableId(group, name, variable, path);
  if (!result)
    return result;
  const unsigned char raw = value ? 1 : 0;
  return checkedNetCDF(nc_put_var_uchar(group, variable, &raw),
                       "Variable write", path + "/" + name);
}

WVCheckpointStatus writeDoubles(int group, const std::string &name,
                                const std::vector<double> &values,
                                const std::string &path) {
  int variable = -1;
  auto result = variableId(group, name, variable, path);
  if (!result)
    return result;
  return checkedNetCDF(nc_put_var_double(group, variable, values.data()),
                       "Variable write", path + "/" + name);
}

WVCheckpointStatus writeForcingEntry(int group,
                                     const WVFrozenForcingEntry &entry,
                                     const WVForcingCatalog &catalog,
                                     const std::string &path) {
  const auto *registration =
      catalog.registration(entry.typeIdentifier, entry.contractVersion);
  if (registration == nullptr || !registration->isSupported)
    return failed(WVCheckpointStatusCode::unsupportedForcing,
                  "Unsupported forcing reached output writing.", path);
  for (const auto &field : registration->persistence.fields) {
    const auto *value = entry.configuration.value(field.recordName);
    if (value == nullptr && field.optional)
      continue;
    if (value == nullptr)
      return failed(WVCheckpointStatusCode::malformedForcing,
                    "Required forcing value is missing.", path);
    WVCheckpointStatus result = WVCheckpointStatus::ok();
    if (field.encoding == WVForcingPersistenceEncoding::realVariable) {
      const auto *values = std::get_if<std::vector<double>>(&value->storage);
      if (values == nullptr)
        return failed(WVCheckpointStatusCode::malformedForcing,
                      "Forcing real value has the wrong type.", path);
      result = writeDoubles(group, field.netcdfName, *values, path);
    } else if (field.encoding == WVForcingPersistenceEncoding::logicalVariable) {
      const auto *values =
          std::get_if<std::vector<std::uint8_t>>(&value->storage);
      if (values == nullptr || values->size() != 1)
        return failed(WVCheckpointStatusCode::malformedForcing,
                      "Forcing logical value is malformed.", path);
      result = writeLogical(group, field.netcdfName, values->front() != 0, path);
    } else if (field.encoding ==
               WVForcingPersistenceEncoding::zeroBasedIndexVariable) {
      const auto *values =
          std::get_if<std::vector<std::int64_t>>(&value->storage);
      if (values == nullptr)
        return failed(WVCheckpointStatusCode::malformedForcing,
                      "Forcing index value has the wrong type.", path);
      if (values->empty())
        continue;
      std::vector<unsigned long long> oneBased(values->size());
      for (std::size_t index = 0; index < values->size(); ++index)
        oneBased[index] = static_cast<unsigned long long>((*values)[index] + 1);
      int variable = -1;
      result = variableId(group, field.netcdfName, variable, path);
      if (result)
        result = checkedNetCDF(
            nc_put_var_ulonglong(group, variable, oneBased.data()),
            "Forcing index write", path + "/" + field.netcdfName);
    } else if (field.encoding == WVForcingPersistenceEncoding::complexVariable) {
      const auto *real = std::get_if<std::vector<double>>(&value->storage);
      const auto *imaginaryValue =
          entry.configuration.value(field.imaginaryRecordName);
      const auto *imag = imaginaryValue == nullptr
                             ? nullptr
                             : std::get_if<std::vector<double>>(
                                   &imaginaryValue->storage);
      if (real == nullptr || imag == nullptr || real->size() != imag->size())
        return failed(WVCheckpointStatusCode::malformedForcing,
                      "Forcing complex value is malformed.", path);
      if (field.dimensions == WVForcingDimensionRule::referencedLength &&
          real->empty())
        continue;
      if (field.dimensions == WVForcingDimensionRule::componentPair) {
        result = writeDoubles(group, "barotropicVelocityComponent", {1.0, 2.0},
                              path);
        if (!result)
          return result;
      }
      result = writeDoubles(group, field.netcdfName + "_real", *real, path);
      if (result)
        result = writeDoubles(group, field.netcdfName + "_imag", *imag, path);
    }
    if (!result)
      return result;
  }
  return WVCheckpointStatus::ok();
}

} // namespace

WVCheckpointStatus checkedNetCDF(int code, const std::string &operation,
                                 const std::string &location) {
  return code == NC_NOERR ? WVCheckpointStatus::ok()
                          : netcdfFailure(code, operation, location);
}

WVCheckpointStatus putTextAttribute(int group, int variable, const char *name,
                                    const std::string &value,
                                    const std::string &location) {
  return checkedNetCDF(
      nc_put_att_text(group, variable, name, value.size(), value.data()),
      "Text-attribute definition", location + "/@" + name);
}

WVCheckpointStatus putByteAttribute(int group, int variable, const char *name,
                                    unsigned char value,
                                    const std::string &location) {
  return checkedNetCDF(
      nc_put_att_uchar(group, variable, name, NC_UBYTE, 1, &value),
      "Byte-attribute definition", location + "/@" + name);
}

WVCheckpointStatus defineDoubleVariable(int group, const std::string &name,
                                        const std::vector<int> &dimensions,
                                        int &variable,
                                        const std::string &path) {
  return checkedNetCDF(
      nc_def_var(group, name.c_str(), NC_DOUBLE,
                 static_cast<int>(dimensions.size()),
                 dimensions.empty() ? nullptr : dimensions.data(), &variable),
      "Variable definition", path + "/" + name);
}

WVCheckpointStatus defineComplexVariable(int group, const std::string &name,
                                         const std::vector<int> &dimensions,
                                         const std::string &path) {
  int real = -1;
  int imag = -1;
  auto result =
      defineDoubleVariable(group, name + "_real", dimensions, real, path);
  if (!result)
    return result;
  result = defineDoubleVariable(group, name + "_imag", dimensions, imag, path);
  if (!result)
    return result;
  for (const auto &marker :
       std::array<std::pair<const char *, unsigned char>, 3>{
           {{"isComplex", 1}, {"isRealPart", 1}, {"isImaginaryPart", 0}}}) {
    result = putByteAttribute(group, real, marker.first, marker.second,
                              path + "/" + name + "_real");
    if (!result)
      return result;
  }
  for (const auto &marker :
       std::array<std::pair<const char *, unsigned char>, 3>{
           {{"isComplex", 1}, {"isRealPart", 0}, {"isImaginaryPart", 1}}}) {
    result = putByteAttribute(group, imag, marker.first, marker.second,
                              path + "/" + name + "_imag");
    if (!result)
      return result;
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus defineModelOutputRoot(
    int root, const WVCheckpoint &checkpoint, const WVForcingCatalog &catalog,
    bool isDynamicsLinear,
    std::array<int, 5> &dimensions, std::vector<int> &forcingGroups,
    std::vector<const WVFrozenForcingEntry *> &forcingEntries) {
  const bool isBarotropicQG =
      checkpoint.transformKind == WVPersistedTransformKind::barotropicQG;
  const auto &legacy = checkpoint.configuration;
  const auto &qg = checkpoint.barotropicQGConfiguration;
  const std::size_t Nkl = isBarotropicQG
                              ? checkpoint.transformState.coefficientFamilies
                                    .front().values.size()
                              : checkpoint.state.coefficients.shape.columns;
  const std::array<std::pair<const char *, std::size_t>, 5> definitions =
      isBarotropicQG
          ? std::array<std::pair<const char *, std::size_t>, 5>{
                {{nullptr, 0}, {"kl", Nkl}, {"x", qg.Nx}, {"y", qg.Ny},
                 {nullptr, 0}}}
          : std::array<std::pair<const char *, std::size_t>, 5>{
                {{"j", legacy.Nj}, {"kl", Nkl}, {"x", legacy.Nx},
                 {"y", legacy.Ny}, {"z", legacy.Nz}}};
  dimensions.fill(-1);
  for (std::size_t index = 0; index < definitions.size(); ++index) {
    if (definitions[index].first == nullptr)
      continue;
    auto result = checkedNetCDF(
        nc_def_dim(root, definitions[index].first, definitions[index].second,
                   &dimensions[index]),
        "Dimension definition", "/" + std::string(definitions[index].first));
    if (!result)
      return result;
    int variable = -1;
    result = defineDoubleVariable(root, definitions[index].first,
                                  {dimensions[index]}, variable, "/");
    if (!result)
      return result;
  }
  const std::vector<const char *> scalars =
      isBarotropicQG
          ? std::vector<const char *>{"Lx", "Ly", "g", "h", "j",
                                      "latitude", "planetaryRadius",
                                      "rotationRate", "t0"}
          : std::vector<const char *>{"Lx", "Ly", "Lz", "N0", "g",
                                      "latitude", "planetaryRadius", "rho0",
                                      "rotationRate", "t0"};
  for (const char *name : scalars) {
    int variable = -1;
    const auto result = defineDoubleVariable(root, name, {}, variable, "/");
    if (!result)
      return result;
  }
  const std::vector<const char *> logicals =
      isBarotropicQG ? std::vector<const char *>{"shouldAntialias"}
                     : std::vector<const char *>{"isHydrostatic",
                                                 "shouldAntialias"};
  for (const char *name : logicals) {
    int variable = -1;
    const auto result = defineLogical(root, name, variable, "/");
    if (!result)
      return result;
  }
  const std::string transformClass =
      isBarotropicQG ? "WVTransformBarotropicQG"
                     : "WVTransformConstantStratification";
  auto result = putTextAttribute(root, NC_GLOBAL, "AnnotatedClass",
                                 transformClass, "/");
  if (!result)
    return result;
  result = putTextAttribute(root, NC_GLOBAL, "WVTransform", transformClass,
                            "/");
  if (!result)
    return result;
  result = putTextAttribute(root, NC_GLOBAL, "model_version",
                            checkpoint.metadata.modelVersion, "/");
  if (!result)
    return result;
  result = putByteAttribute(root, NC_GLOBAL, "WVModelIsDynamicsLinear",
                            isDynamicsLinear ? 1 : 0, "/");
  if (!result)
    return result;
  result =
      putTextAttribute(root, NC_GLOBAL, "source",
                       "Created with the WaveVortex portable runtime", "/");
  if (!result)
    return result;
  result = putTextAttribute(root, NC_GLOBAL, "history",
                            "MATLAB-compatible multi-group output initialized "
                            "by the WaveVortex portable runtime.",
                            "/");
  if (!result)
    return result;
  const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm utc{};
#if defined(_WIN32)
  gmtime_s(&utc, &now);
#else
  gmtime_r(&now, &utc);
#endif
  std::ostringstream date;
  date << std::put_time(&utc, "%Y-%m-%d %H:%M:%S UTC");
  result = putTextAttribute(root, NC_GLOBAL, "date_created", date.str(), "/");
  if (!result)
    return result;

  forcingEntries.clear();
  for (const auto &entry : checkpoint.forcingSchedule.entries)
    forcingEntries.push_back(&entry);
  std::sort(forcingEntries.begin(), forcingEntries.end(),
            [](const auto *left, const auto *right) {
              return left->ordinal < right->ordinal;
            });
  if (forcingEntries.empty())
    return WVCheckpointStatus::ok();
  int forcingRoot = -1;
  result = checkedNetCDF(nc_def_grp(root, "forcing", &forcingRoot),
                         "Forcing-group definition", "/forcing");
  if (!result)
    return result;
  forcingGroups.reserve(forcingEntries.size());
  if (forcingEntries.size() == 1) {
    forcingGroups.push_back(forcingRoot);
    return defineForcingEntry(forcingRoot, *forcingEntries.front(), catalog, dimensions,
                              "/forcing");
  }
  for (std::size_t index = 0; index < forcingEntries.size(); ++index) {
    int group = -1;
    const std::string name = "forcing-" + std::to_string(index + 1);
    result = checkedNetCDF(nc_def_grp(forcingRoot, name.c_str(), &group),
                           "Forcing-group definition", "/forcing/" + name);
    if (!result)
      return result;
    forcingGroups.push_back(group);
    result = defineForcingEntry(group, *forcingEntries[index], catalog, dimensions,
                                "/forcing/" + name);
    if (!result)
      return result;
  }
  return WVCheckpointStatus::ok();
}

WVCheckpointStatus writeModelOutputRoot(
    int root, const WVCheckpoint &checkpoint,
    const WVForcingCatalog &catalog,
    const std::vector<int> &forcingGroups,
    const std::vector<const WVFrozenForcingEntry *> &forcingEntries) {
  const bool isBarotropicQG =
      checkpoint.transformKind == WVPersistedTransformKind::barotropicQG;
  const auto &configuration = checkpoint.configuration;
  const auto &qg = checkpoint.barotropicQGConfiguration;
  if (isBarotropicQG) {
    const std::size_t Nkl =
        checkpoint.transformState.coefficientFamilies.front().values.size();
    std::vector<double> kl(Nkl), x(qg.Nx), y(qg.Ny);
    for (std::size_t index = 0; index < kl.size(); ++index)
      kl[index] = static_cast<double>(index);
    for (std::size_t index = 0; index < x.size(); ++index)
      x[index] = static_cast<double>(index) * qg.Lx /
                 static_cast<double>(qg.Nx);
    for (std::size_t index = 0; index < y.size(); ++index)
      y[index] = static_cast<double>(index) * qg.Ly /
                 static_cast<double>(qg.Ny);
    const std::array<std::pair<const char *, const std::vector<double> *>, 3>
        coordinates{{{"kl", &kl}, {"x", &x}, {"y", &y}}};
    for (const auto &coordinate : coordinates) {
      const auto result =
          writeDoubles(root, coordinate.first, *coordinate.second, "/");
      if (!result)
        return result;
    }
    const std::array<std::pair<const char *, double>, 9> scalars{{
        {"Lx", qg.Lx},
        {"Ly", qg.Ly},
        {"g", qg.g},
        {"h", qg.h},
        {"j", static_cast<double>(qg.j)},
        {"latitude", qg.latitude},
        {"planetaryRadius", qg.planetaryRadius},
        {"rotationRate", qg.rotationRate},
        {"t0", checkpoint.state.t0},
    }};
    for (const auto &scalar : scalars) {
      const auto result = writeDouble(root, scalar.first, scalar.second, "/");
      if (!result)
        return result;
    }
    auto result =
        writeLogical(root, "shouldAntialias", qg.shouldAntialias, "/");
    if (!result)
      return result;
    for (std::size_t index = 0; index < forcingEntries.size(); ++index) {
      const std::string path =
          forcingEntries.size() == 1
              ? "/forcing"
              : "/forcing/forcing-" + std::to_string(index + 1);
      result = writeForcingEntry(forcingGroups[index], *forcingEntries[index],
                                 catalog, path);
      if (!result)
        return result;
    }
    return WVCheckpointStatus::ok();
  }
  std::vector<double> j(configuration.Nj),
      kl(checkpoint.state.coefficients.shape.columns), x(configuration.Nx),
      y(configuration.Ny), z(configuration.Nz);
  for (std::size_t index = 0; index < j.size(); ++index)
    j[index] = static_cast<double>(index);
  for (std::size_t index = 0; index < kl.size(); ++index)
    kl[index] = static_cast<double>(index);
  for (std::size_t index = 0; index < x.size(); ++index)
    x[index] = static_cast<double>(index) * configuration.Lx /
               static_cast<double>(configuration.Nx);
  for (std::size_t index = 0; index < y.size(); ++index)
    y[index] = static_cast<double>(index) * configuration.Ly /
               static_cast<double>(configuration.Ny);
  for (std::size_t index = 0; index < z.size(); ++index)
    z[index] =
        -configuration.Lz + static_cast<double>(index) * configuration.Lz /
                                static_cast<double>(configuration.Nz - 1);
  for (const auto &coordinate :
       std::array<std::pair<const char *, const std::vector<double> *>, 5>{
           {{"j", &j}, {"kl", &kl}, {"x", &x}, {"y", &y}, {"z", &z}}}) {
    const auto result =
        writeDoubles(root, coordinate.first, *coordinate.second, "/");
    if (!result)
      return result;
  }
  for (const auto &scalar : std::array<std::pair<const char *, double>, 10>{
           {{"Lx", configuration.Lx},
            {"Ly", configuration.Ly},
            {"Lz", configuration.Lz},
            {"N0", configuration.N0},
            {"g", configuration.g},
            {"latitude", configuration.latitude},
            {"planetaryRadius", configuration.planetaryRadius},
            {"rho0", configuration.rho0},
            {"rotationRate", configuration.rotationRate},
            {"t0", checkpoint.state.t0}}}) {
    const auto result = writeDouble(root, scalar.first, scalar.second, "/");
    if (!result)
      return result;
  }
  auto result =
      writeLogical(root, "isHydrostatic", configuration.isHydrostatic, "/");
  if (!result)
    return result;
  result =
      writeLogical(root, "shouldAntialias", configuration.shouldAntialias, "/");
  if (!result)
    return result;
  for (std::size_t index = 0; index < forcingEntries.size(); ++index) {
    const std::string path =
        forcingEntries.size() == 1
            ? "/forcing"
            : "/forcing/forcing-" + std::to_string(index + 1);
    result =
        writeForcingEntry(forcingGroups[index], *forcingEntries[index], catalog, path);
    if (!result)
      return result;
  }
  return WVCheckpointStatus::ok();
}

namespace {

class ConstantStratificationOutputTransformAdapter final
    : public WVModelOutputTransformAdapter {
public:
  explicit ConstantStratificationOutputTransformAdapter(
      WVTransformConstantStratificationConfiguration configuration)
      : configuration_(std::move(configuration)) {}

  WVCheckpointStatus validate(
      const WVCheckpoint &checkpoint,
      const WVIntegrationStateLayout &layout) const override {
    if (checkpoint.transformKind !=
            WVPersistedTransformKind::constantStratification ||
        (!checkpoint.stateDescription.transformIdentifier.empty() &&
         checkpoint.stateDescription.transformIdentifier !=
             layout.transformIdentifier()))
      return failed(WVCheckpointStatusCode::shapeMismatch,
                    "Output checkpoint and state-layout transforms differ.",
                    "/");
    WVTransformConstantStratificationDescriptor transform;
    const auto status = WVTransformConstantStratificationDescriptor::create(
        configuration_, transform);
    if (!status)
      return failed(WVCheckpointStatusCode::descriptorFailure,
                    status.message, "/");
    const auto shape = layout.coefficientShape();
    if (!layout.hasLegacyCoefficientTriple() ||
        shape.rows != configuration_.Nj ||
        shape.columns != transform.Nkl())
      return failed(WVCheckpointStatusCode::shapeMismatch,
                    "Output state layout and transform configuration differ.",
                    "/");
    return WVCheckpointStatus::ok();
  }

  const WVComplex64 *coefficientData(
      const WVCheckpoint &checkpoint, const WVIntegrationStateLayout &,
      std::size_t family) const noexcept override {
    const WVComplex64 *values[] = {
        checkpoint.state.coefficients.Ap.data(),
        checkpoint.state.coefficients.Am.data(),
        checkpoint.state.coefficients.A0.data()};
    return family < 3 ? values[family] : nullptr;
  }

  void bindConstructionState(
      const WVCheckpoint &checkpoint,
      WVIntegrationState &state) const noexcept override {
    state.waveVortex = checkpoint.state.view();
  }

  WVCheckpointStatus defineRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog, bool isDynamicsLinear,
      std::array<int, 5> &dimensions, std::vector<int> &forcingGroups,
      std::vector<const WVFrozenForcingEntry *> &forcingEntries) const override {
    return defineModelOutputRoot(root, checkpoint, catalog, isDynamicsLinear,
                                 dimensions, forcingGroups, forcingEntries);
  }

  WVCheckpointStatus writeRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog,
      const std::vector<int> &forcingGroups,
      const std::vector<const WVFrozenForcingEntry *> &forcingEntries)
      const override {
    return writeModelOutputRoot(root, checkpoint, catalog, forcingGroups,
                                forcingEntries);
  }

  bool sameConfiguration(
      const WVCheckpointInspection &inspection) const noexcept override {
    return inspection.transformKind ==
               WVPersistedTransformKind::constantStratification &&
           sameTransformConfiguration(configuration_, inspection.configuration);
  }

  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }

private:
  WVTransformConstantStratificationConfiguration configuration_;
};

class BarotropicQGOutputTransformAdapter final
    : public WVModelOutputTransformAdapter {
public:
  explicit BarotropicQGOutputTransformAdapter(
      WVTransformBarotropicQGConfiguration configuration)
      : configuration_(std::move(configuration)) {}

  WVCheckpointStatus validate(
      const WVCheckpoint &checkpoint,
      const WVIntegrationStateLayout &layout) const override {
    if (checkpoint.transformKind != WVPersistedTransformKind::barotropicQG ||
        checkpoint.stateDescription.transformIdentifier !=
            layout.transformIdentifier())
      return failed(WVCheckpointStatusCode::shapeMismatch,
                    "Output checkpoint and state-layout transforms differ.",
                    "/");
    WVTransformBarotropicQGDescriptor transform;
    const auto status = WVTransformBarotropicQGDescriptor::create(
        configuration_, transform);
    if (!status)
      return failed(WVCheckpointStatusCode::descriptorFailure,
                    status.message, "/");
    if (layout.coefficientFamilyCount() != 1 ||
        layout.coefficientFamilies()[0].identifier != "A0" ||
        layout.coefficientFamilies()[0].spectralDimensions !=
            std::vector<std::size_t>{transform.Nkl()} ||
        checkpoint.transformState.coefficientFamilies.size() != 1 ||
        checkpoint.transformState.coefficientFamilies[0].values.size() !=
            transform.Nkl())
      return failed(WVCheckpointStatusCode::shapeMismatch,
                    "Compact Barotropic QG output state does not match its "
                    "transform.",
                    "/");
    return WVCheckpointStatus::ok();
  }

  const WVComplex64 *coefficientData(
      const WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout,
      std::size_t family) const noexcept override {
    return family < layout.coefficientFamilyCount() &&
                   family < checkpoint.transformState.coefficientFamilies.size()
               ? checkpoint.transformState.coefficientFamilies[family]
                     .values.data()
               : nullptr;
  }

  void bindConstructionState(
      const WVCheckpoint &checkpoint,
      WVIntegrationState &state) const noexcept override {
    state.waveVortex.t = checkpoint.transformState.t;
    state.waveVortex.t0 = checkpoint.transformState.t0;
  }

  WVCheckpointStatus defineRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog, bool isDynamicsLinear,
      std::array<int, 5> &dimensions, std::vector<int> &forcingGroups,
      std::vector<const WVFrozenForcingEntry *> &forcingEntries) const override {
    return defineModelOutputRoot(root, checkpoint, catalog, isDynamicsLinear,
                                 dimensions, forcingGroups, forcingEntries);
  }

  WVCheckpointStatus writeRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog,
      const std::vector<int> &forcingGroups,
      const std::vector<const WVFrozenForcingEntry *> &forcingEntries)
      const override {
    return writeModelOutputRoot(root, checkpoint, catalog, forcingGroups,
                                forcingEntries);
  }

  bool sameConfiguration(
      const WVCheckpointInspection &inspection) const noexcept override {
    return inspection.transformKind ==
               WVPersistedTransformKind::barotropicQG &&
           sameTransformConfiguration(
               configuration_, inspection.barotropicQGConfiguration);
  }

  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this);
  }

private:
  WVTransformBarotropicQGConfiguration configuration_;
};

} // namespace

WVCheckpointStatus createModelOutputTransformAdapter(
    const WVCheckpoint &checkpoint,
    std::unique_ptr<WVModelOutputTransformAdapter> &adapter) {
  adapter.reset();
  try {
    if (checkpoint.transformKind == WVPersistedTransformKind::barotropicQG)
      adapter = std::make_unique<BarotropicQGOutputTransformAdapter>(
          checkpoint.barotropicQGConfiguration);
    else
      adapter =
          std::make_unique<ConstantStratificationOutputTransformAdapter>(
              checkpoint.configuration);
    return WVCheckpointStatus::ok();
  } catch (const std::bad_alloc &) {
    return failed(WVCheckpointStatusCode::writeFailure,
                  "Unable to allocate the model-output transform adapter.",
                  "/");
  }
}

bool sameModelOutputTransformConfiguration(
    const WVCheckpointInspection &left,
    const WVCheckpointInspection &right) noexcept {
  if (left.transformKind != right.transformKind)
    return false;
  return left.transformKind == WVPersistedTransformKind::barotropicQG
             ? sameTransformConfiguration(left.barotropicQGConfiguration,
                                          right.barotropicQGConfiguration)
             : sameTransformConfiguration(left.configuration,
                                          right.configuration);
}

bool modelOutputGroupCarriesCompleteCoefficientRestart(
    const WVTransformStateDescription &description,
    bool hasDeclaredCoefficientFamilies,
    bool hasCoefficientObserver) noexcept {
  if (!hasDeclaredCoefficientFamilies)
    return false;
  // In the compact QG contract A0 is also a valid Eulerian field name, so its
  // variable alone cannot identify restart ownership. The legacy constant-
  // stratification contract historically treats a complete Ap/Am/A0 field
  // triple as restart state, including linear passive-field output.
  return description.transformIdentifier != "WVTransformBarotropicQG" ||
         hasCoefficientObserver;
}

} // namespace wavevortex::runtime::detail
