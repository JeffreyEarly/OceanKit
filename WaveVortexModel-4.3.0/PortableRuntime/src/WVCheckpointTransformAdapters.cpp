#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include "WVForcingScheduleDecoder.hpp"
#include "WVLegacyObservationNetCDFAdapter.hpp"
#include "WVNetCDF.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace wavevortex::runtime {

namespace {

using detail::WVNetCDFFile;

struct GroupRecord {
    int id = -1;
    std::string path;
};

struct StateGroupRecord {
    int id = -1;
    std::string path;
};

WVCheckpointStatus status(WVCheckpointStatusCode code, std::string message, std::string location) {
    return {code, std::move(message), std::move(location)};
}

WVCheckpointStatus readOptionalTextAttribute(int groupId, const std::string& name, std::string& value, bool& present, const std::string& groupPath) {
    nc_type type = NC_NAT;
    std::size_t length = 0;
    const int result = nc_inq_att(groupId, NC_GLOBAL, name.c_str(), &type, &length);
    if (result == NC_ENOTATT) {
        present = false;
        value.clear();
        return WVCheckpointStatus::ok();
    }
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Attribute inspection", groupPath + "/@" + name);
    present = true;
    return detail::readTextAttribute(groupId, name, value, groupPath);
}

WVCheckpointStatus inspectGroupTree(int groupId, const std::string& path, std::vector<GroupRecord>& records) {
    records.push_back({groupId, path});
    std::vector<int> children;
    auto result = detail::childGroups(groupId, children, path);
    if (!result) return result;
    for (const int child : children) {
        std::string name;
        result = detail::groupName(child, name, path);
        if (!result) return result;
        const std::string childPath = path == "/" ? "/" + name : path + "/" + name;
        result = inspectGroupTree(child, childPath, records);
        if (!result) return result;
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus validateVersion(const std::string& version) {
    const auto separator = version.find('.');
    const std::string majorText = separator == std::string::npos ? version : version.substr(0, separator);
    int major = -1;
    const auto conversion = std::from_chars(majorText.data(), majorText.data() + majorText.size(), major);
    if (majorText.empty() || conversion.ec != std::errc{} || conversion.ptr != majorText.data() + majorText.size() || major != 4) {
        return status(WVCheckpointStatusCode::unsupportedModelVersion,
            "The structural checkpoint profile accepts WaveVortexModel 4.x files; found model_version '" + version + "'.",
            "/@model_version");
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus validateCoordinate(const std::string& name, const std::vector<double>& values) {
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (!std::isfinite(values[index])) {
            return status(WVCheckpointStatusCode::invalidValue, "Coordinate '" + name + "' contains a non-finite value.", "/" + name);
        }
        if (index > 0 && !(values[index] > values[index - 1])) {
            return status(WVCheckpointStatusCode::invalidValue, "Coordinate '" + name + "' must be strictly increasing.", "/" + name);
        }
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readConfiguration(int rootId, WVTransformConstantStratificationConfiguration& configuration) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    for (auto pair : {std::pair<const char*, std::vector<double>*>("x", &x), {"y", &y}, {"z", &z}}) {
        auto result = detail::readDoubleCoordinate(rootId, pair.first, *pair.second, "/");
        if (!result) return result;
        result = validateCoordinate(pair.first, *pair.second);
        if (!result) return result;
    }

    configuration.Nx = x.size();
    configuration.Ny = y.size();
    configuration.Nz = z.size();
    struct ScalarField { const char* name; double* value; };
    for (const ScalarField field : {
            ScalarField{"Lx", &configuration.Lx}, {"Ly", &configuration.Ly}, {"Lz", &configuration.Lz},
            {"N0", &configuration.N0}, {"g", &configuration.g}, {"rho0", &configuration.rho0},
            {"planetaryRadius", &configuration.planetaryRadius}, {"rotationRate", &configuration.rotationRate},
            {"latitude", &configuration.latitude}}) {
        const auto result = detail::readDoubleScalar(rootId, field.name, *field.value, "/");
        if (!result) return result;
        if (!std::isfinite(*field.value)) return status(WVCheckpointStatusCode::invalidValue, "Configuration variable '" + std::string(field.name) + "' must be finite.", "/" + std::string(field.name));
    }
    auto result = detail::readLogicalScalar(rootId, "isHydrostatic", configuration.isHydrostatic, "/");
    if (!result) return result;
    result = detail::readLogicalScalar(rootId, "shouldAntialias", configuration.shouldAntialias, "/");
    if (!result) return result;
    if (configuration.Nx < 2 || configuration.Ny < 2 || configuration.Nz < 2 ||
        configuration.Lx <= 0.0 || configuration.Ly <= 0.0 || configuration.Lz <= 0.0 ||
        configuration.N0 <= 0.0 || configuration.g <= 0.0 || configuration.rho0 <= 0.0 ||
        configuration.planetaryRadius <= 0.0 || configuration.rotationRate <= 0.0 ||
        configuration.latitude < -90.0 || configuration.latitude > 90.0) {
        return status(WVCheckpointStatusCode::invalidValue, "The constant-stratification configuration contains an invalid size or physical scalar.", "/");
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readBarotropicQGConfiguration(
    int rootId, WVTransformBarotropicQGConfiguration& configuration) {
    std::vector<double> x;
    std::vector<double> y;
    for (auto pair : {std::pair<const char*, std::vector<double>*>{"x", &x},
                      {"y", &y}}) {
        auto result = detail::readDoubleCoordinate(rootId, pair.first,
                                                   *pair.second, "/");
        if (!result) return result;
        result = validateCoordinate(pair.first, *pair.second);
        if (!result) return result;
    }
    configuration.Nx = x.size();
    configuration.Ny = y.size();
    struct ScalarField { const char* name; double* value; };
    for (const ScalarField field : {
            ScalarField{"Lx", &configuration.Lx},
            {"Ly", &configuration.Ly}, {"h", &configuration.h},
            {"g", &configuration.g},
            {"planetaryRadius", &configuration.planetaryRadius},
            {"rotationRate", &configuration.rotationRate},
            {"latitude", &configuration.latitude}}) {
        auto result = detail::readDoubleScalar(rootId, field.name,
                                                *field.value, "/");
        if (!result) return result;
        if (!std::isfinite(*field.value))
            return status(WVCheckpointStatusCode::invalidValue,
                          "Barotropic QG configuration variable '" +
                              std::string(field.name) + "' must be finite.",
                          "/" + std::string(field.name));
    }
    double j = 0.0;
    auto result = detail::readDoubleScalar(rootId, "j", j, "/");
    if (!result) return result;
    if (!std::isfinite(j) || (j != 0.0 && j != 1.0))
        return status(WVCheckpointStatusCode::invalidValue,
                      "Barotropic QG mode j must be exactly 0 or 1.", "/j");
    configuration.j = static_cast<std::uint32_t>(j);
    result = detail::readLogicalScalar(rootId, "shouldAntialias",
                                       configuration.shouldAntialias, "/");
    if (!result) return result;
    WVTransformBarotropicQGDescriptor descriptor;
    const auto descriptorStatus =
        WVTransformBarotropicQGDescriptor::create(configuration, descriptor);
    if (!descriptorStatus)
        return status(WVCheckpointStatusCode::descriptorFailure,
                      "Unable to rebuild the Barotropic QG descriptor: " +
                          descriptorStatus.message,
                      "/");
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus groupCarriesCoefficientRestart(
    const GroupRecord& group, const WVExtensionCatalog& catalog,
    bool& carriesCoefficientRestart) {
    carriesCoefficientRestart = true;
    int metadataRoot = -1;
    const int lookup =
        nc_inq_ncid(group.id, "observingSystems", &metadataRoot);
    if (lookup == NC_ENOGRP)
        return WVCheckpointStatus::ok();
    if (lookup != NC_NOERR)
        return detail::netcdfFailure(
            lookup, "Observer metadata-root lookup",
            group.path + "/observingSystems");
    carriesCoefficientRestart =
        detail::persistedObserverCarriesCoefficientState(metadataRoot,
                                                         catalog);
    std::vector<int> observers;
    auto result = detail::childGroups(metadataRoot, observers,
                                      group.path + "/observingSystems");
    if (!result) return result;
    for (const int observer : observers)
        carriesCoefficientRestart =
            carriesCoefficientRestart ||
            detail::persistedObserverCarriesCoefficientState(observer,
                                                               catalog);
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus findStateGroup(const std::vector<GroupRecord>& groups,
                                  StateGroupRecord& stateGroup) {
    std::vector<StateGroupRecord> candidates;
    for (const auto& group : groups) {
        std::size_t completeFamilies = 0;
        for (const char* family : {"Ap", "Am", "A0"}) {
            int plainId = -1;
            int realId = -1;
            int imagId = -1;
            bool plain = false;
            bool real = false;
            bool imag = false;
            auto result = detail::variableIdIfPresent(group.id, family, plainId, plain, group.path);
            if (!result) return result;
            result = detail::variableIdIfPresent(group.id, std::string(family) + "_real", realId, real, group.path);
            if (!result) return result;
            result = detail::variableIdIfPresent(group.id, std::string(family) + "_imag", imagId, imag, group.path);
            if (!result) return result;
            if (plain && (real || imag))
                return status(WVCheckpointStatusCode::ambiguousState, "A coefficient family cannot contain both a plain-real variable and a real/imaginary pair.", group.path + "/" + family);
            if (real != imag)
                return status(WVCheckpointStatusCode::missingComplexPartner, "A complex coefficient family is missing its real or imaginary partner.", group.path + "/" + family);
            if (plain || (real && imag)) ++completeFamilies;
        }
        if (completeFamilies == 0) continue;
        if (completeFamilies != 3) {
            return status(WVCheckpointStatusCode::missingComplexPartner, "A checkpoint state group must contain complete Ap, Am, and A0 real/imaginary pairs.", group.path);
        }
        int timeId = -1;
        bool hasTime = false;
        const auto result = detail::variableIdIfPresent(
            group.id, "t", timeId, hasTime, group.path);
        if (!result) return result;
        if (!hasTime) return status(WVCheckpointStatusCode::missingVariable, "The checkpoint state group does not contain time variable 't'.", group.path + "/t");
        candidates.push_back({group.id, group.path});
    }
    if (candidates.empty()) return status(WVCheckpointStatusCode::missingVariable, "No complete Ap, Am, A0 checkpoint state was found.", "/");
    if (candidates.size() != 1) return status(WVCheckpointStatusCode::ambiguousState, "A checkpoint must contain exactly one complete Ap, Am, A0 state group.", "/");
    stateGroup = candidates.front();
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus findBarotropicQGStateGroup(
    const std::vector<GroupRecord>& groups,
    const WVExtensionCatalog& catalog, StateGroupRecord& stateGroup) {
    std::vector<StateGroupRecord> candidates;
    for (const auto& group : groups) {
        int realId = -1;
        int imagId = -1;
        int plainId = -1;
        bool real = false;
        bool imag = false;
        bool plain = false;
        auto result = detail::variableIdIfPresent(
            group.id, "A0_real", realId, real, group.path);
        if (!result) return result;
        result = detail::variableIdIfPresent(
            group.id, "A0_imag", imagId, imag, group.path);
        if (!result) return result;
        result = detail::variableIdIfPresent(
            group.id, "A0", plainId, plain, group.path);
        if (!result) return result;
        if (!plain && !real && !imag) continue;
        bool carriesCoefficientRestart = false;
        result = groupCarriesCoefficientRestart(
            group, catalog, carriesCoefficientRestart);
        if (!result) return result;
        if (!carriesCoefficientRestart) continue;
        if (plain && (real || imag))
            return status(WVCheckpointStatusCode::ambiguousState,
                          "Barotropic QG A0 cannot have both plain and "
                          "complex-pair encodings.",
                          group.path + "/A0");
        if (real != imag)
            return status(WVCheckpointStatusCode::missingComplexPartner,
                          "Barotropic QG A0 is missing a complex partner.",
                          group.path + "/A0");
        for (const char* forbidden : {"Ap", "Ap_real", "Ap_imag", "Am",
                                      "Am_real", "Am_imag"}) {
            int id = -1;
            bool present = false;
            result = detail::variableIdIfPresent(group.id, forbidden, id,
                                                 present, group.path);
            if (!result) return result;
            if (present)
                return status(WVCheckpointStatusCode::schemaMismatch,
                              "Barotropic QG state must not contain dummy Ap "
                              "or Am variables.",
                              group.path + "/" + forbidden);
        }
        int timeId = -1;
        bool hasTime = false;
        result = detail::variableIdIfPresent(group.id, "t", timeId, hasTime,
                                             group.path);
        if (!result) return result;
        if (!hasTime)
            return status(WVCheckpointStatusCode::missingVariable,
                          "The Barotropic QG state group has no time variable.",
                          group.path + "/t");
        candidates.push_back({group.id, group.path});
    }
    if (candidates.empty())
        return status(WVCheckpointStatusCode::missingVariable,
                      "No compact Barotropic QG A0 state was found.", "/");
    if (candidates.size() != 1)
        return status(WVCheckpointStatusCode::ambiguousState,
                      "A Barotropic QG checkpoint must contain one compact "
                      "A0 state group.",
                      "/");
    stateGroup = candidates.front();
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus inquireVariable(int groupId, int variableId, nc_type& type, std::vector<int>& dimensionIds, const std::string& location) {
    int dimensionCount = 0;
    int result = nc_inq_varndims(groupId, variableId, &dimensionCount);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Variable-dimension inspection", location);
    result = nc_inq_vartype(groupId, variableId, &type);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Variable-type inspection", location);
    dimensionIds.resize(static_cast<std::size_t>(dimensionCount));
    if (dimensionCount > 0) {
        result = nc_inq_vardimid(groupId, variableId, dimensionIds.data());
        if (result != NC_NOERR) return detail::netcdfFailure(result, "Variable-dimension inspection", location);
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus dimensionName(int groupId, int dimensionId, std::string& name, std::size_t& length, const std::string& location) {
    char rawName[NC_MAX_NAME + 1] = {};
    const int result = nc_inq_dim(groupId, dimensionId, rawName, &length);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Dimension inspection", location);
    name = rawName;
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readByteVariableAttribute(int groupId, int variableId, const std::string& variableName, const std::string& attributeName, unsigned char expected, const std::string& groupPath) {
    nc_type type = NC_NAT;
    std::size_t length = 0;
    int result = nc_inq_att(groupId, variableId, attributeName.c_str(), &type, &length);
    if (result == NC_ENOTATT) return status(WVCheckpointStatusCode::missingAttribute, "Complex component '" + variableName + "' is missing attribute '" + attributeName + "'.", groupPath + "/" + variableName + "/@" + attributeName);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Variable-attribute inspection", groupPath + "/" + variableName + "/@" + attributeName);
    if ((type != NC_UBYTE && type != NC_BYTE) || length != 1) return status(WVCheckpointStatusCode::typeMismatch, "Complex marker attribute '" + attributeName + "' must be a scalar byte.", groupPath + "/" + variableName + "/@" + attributeName);
    unsigned char value = 0;
    result = nc_get_att_uchar(groupId, variableId, attributeName.c_str(), &value);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Variable-attribute read", groupPath + "/" + variableName + "/@" + attributeName);
    if (value != expected) return status(WVCheckpointStatusCode::invalidValue, "Complex marker attribute '" + attributeName + "' has an unexpected value.", groupPath + "/" + variableName + "/@" + attributeName);
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus inspectTime(int groupId, const std::string& groupPath, WVCheckpointStateSelection selection, std::size_t& selectedIndex, std::size_t& stateCount, double& time) {
    int variableId = -1;
    int result = nc_inq_varid(groupId, "t", &variableId);
    if (result == NC_ENOTVAR) return status(WVCheckpointStatusCode::missingVariable, "Missing state time variable 't'.", groupPath + "/t");
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Time-variable lookup", groupPath + "/t");
    nc_type type = NC_NAT;
    std::vector<int> dimensions;
    auto checkpointStatus = inquireVariable(groupId, variableId, type, dimensions, groupPath + "/t");
    if (!checkpointStatus) return checkpointStatus;
    if (type != NC_DOUBLE || dimensions.size() > 1) return status(WVCheckpointStatusCode::shapeMismatch, "State time must be a scalar or one-dimensional double.", groupPath + "/t");
    if (dimensions.empty()) {
        stateCount = 1;
        if (selection.kind == WVCheckpointStateSelectionKind::index && selection.index != 0) return status(WVCheckpointStatusCode::stateIndexOutOfRange, "State index is out of range for a scalar checkpoint.", groupPath + "/t");
        selectedIndex = 0;
        result = nc_get_var_double(groupId, variableId, &time);
    } else {
        std::string name;
        checkpointStatus = dimensionName(groupId, dimensions.front(), name, stateCount, groupPath + "/t");
        if (!checkpointStatus) return checkpointStatus;
        if (name != "t" || stateCount == 0) return status(WVCheckpointStatusCode::shapeMismatch, "State time must use a nonempty 't' dimension.", groupPath + "/t");
        selectedIndex = selection.kind == WVCheckpointStateSelectionKind::latest ? stateCount - 1 : selection.index;
        if (selectedIndex >= stateCount) return status(WVCheckpointStatusCode::stateIndexOutOfRange, "Requested checkpoint state index is out of range.", groupPath + "/t");
        const std::size_t start[] = {selectedIndex};
        const std::size_t count[] = {1};
        result = nc_get_vara_double(groupId, variableId, start, count, &time);
    }
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Time-variable read", groupPath + "/t");
    if (!std::isfinite(time)) return status(WVCheckpointStatusCode::invalidValue, "Checkpoint time must be finite.", groupPath + "/t");
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readComplexCoefficient(
    int groupId,
    const std::string& groupPath,
    const std::string& baseName,
    std::size_t selectedIndex,
    std::size_t stateCount,
    std::size_t Nj,
    std::size_t Nkl,
    std::vector<WVComplex64>* output) {
    const std::string realName = baseName + "_real";
    const std::string imagName = baseName + "_imag";
    int realId = -1;
    int imagId = -1;
    int result = nc_inq_varid(groupId, realName.c_str(), &realId);
    if (result == NC_ENOTVAR) {
        int plainId = -1;
        result = nc_inq_varid(groupId, baseName.c_str(), &plainId);
        if (result == NC_ENOTVAR) return status(WVCheckpointStatusCode::missingComplexPartner, "Missing coefficient '" + baseName + "' and real component '" + realName + "'.", groupPath + "/" + baseName);
        if (result != NC_NOERR) return detail::netcdfFailure(result, "Coefficient lookup", groupPath + "/" + baseName);
        nc_type plainType = NC_NAT;
        std::vector<int> plainDimensions;
        auto plainStatus = inquireVariable(groupId, plainId, plainType, plainDimensions, groupPath + "/" + baseName);
        if (!plainStatus) return plainStatus;
        if (plainType != NC_DOUBLE || plainDimensions.size() != 2)
            return status(WVCheckpointStatusCode::shapeMismatch, "A plain-real initial coefficient must use double [kl,j] storage.", groupPath + "/" + baseName);
        const std::array<const char*, 2> expectedNames = {"kl", "j"};
        const std::array<std::size_t, 2> expectedLengths = {Nkl, Nj};
        for (std::size_t index = 0; index < plainDimensions.size(); ++index) {
            std::string name;
            std::size_t length = 0;
            plainStatus = dimensionName(groupId, plainDimensions[index], name, length, groupPath + "/" + baseName);
            if (!plainStatus) return plainStatus;
            if (name != expectedNames[index] || length != expectedLengths[index])
                return status(WVCheckpointStatusCode::shapeMismatch, "A plain-real initial coefficient has incompatible dimensions.", groupPath + "/" + baseName);
        }
        if (output != nullptr) {
            std::vector<double> realValues(Nj * Nkl);
            result = nc_get_var_double(groupId, plainId, realValues.data());
            if (result != NC_NOERR) return detail::netcdfFailure(result, "Initial coefficient read", groupPath + "/" + baseName);
            output->resize(realValues.size());
            for (std::size_t index = 0; index < realValues.size(); ++index)
                (*output)[index] = {realValues[index], 0.0};
        }
        return WVCheckpointStatus::ok();
    }
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Complex-component lookup", groupPath + "/" + realName);
    result = nc_inq_varid(groupId, imagName.c_str(), &imagId);
    if (result == NC_ENOTVAR) return status(WVCheckpointStatusCode::missingComplexPartner, "Missing imaginary component '" + imagName + "'.", groupPath + "/" + imagName);
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Complex-component lookup", groupPath + "/" + imagName);

    nc_type realType = NC_NAT;
    nc_type imagType = NC_NAT;
    std::vector<int> realDimensions;
    std::vector<int> imagDimensions;
    auto checkpointStatus = inquireVariable(groupId, realId, realType, realDimensions, groupPath + "/" + realName);
    if (!checkpointStatus) return checkpointStatus;
    checkpointStatus = inquireVariable(groupId, imagId, imagType, imagDimensions, groupPath + "/" + imagName);
    if (!checkpointStatus) return checkpointStatus;
    if (realType != NC_DOUBLE || imagType != NC_DOUBLE) return status(WVCheckpointStatusCode::typeMismatch, "Complex coefficient components must be double precision.", groupPath + "/" + baseName);
    if (realDimensions != imagDimensions) return status(WVCheckpointStatusCode::shapeMismatch, "Complex coefficient components must have identical dimensions.", groupPath + "/" + baseName);
    // Linear WaveVortex files keep Ap/Am/A0 as initial-only [kl,j]
    // variables while the group's time axis continues to grow. Coefficient
    // cadence is therefore determined by the coefficient rank, not by the
    // number of output records in the group.
    const bool timeSeries = realDimensions.size() == 3;
    const std::size_t expectedDimensions = timeSeries ? 3 : 2;
    if (realDimensions.size() != expectedDimensions) return status(WVCheckpointStatusCode::shapeMismatch, "Coefficient '" + baseName + "' must use [kl,j] or [t,kl,j] NetCDF dimensions.", groupPath + "/" + baseName);

    const std::array<const char*, 3> expectedNames = {"t", "kl", "j"};
    const std::size_t expectedOffset = timeSeries ? 0 : 1;
    const std::array<std::size_t, 3> expectedLengths = {stateCount, Nkl, Nj};
    for (std::size_t index = 0; index < realDimensions.size(); ++index) {
        std::string name;
        std::size_t length = 0;
        checkpointStatus = dimensionName(groupId, realDimensions[index], name, length, groupPath + "/" + baseName);
        if (!checkpointStatus) return checkpointStatus;
        const std::size_t expectedIndex = index + expectedOffset;
        if (name != expectedNames[expectedIndex] || length != expectedLengths[expectedIndex]) return status(WVCheckpointStatusCode::shapeMismatch, "Coefficient '" + baseName + "' has incompatible NetCDF dimensions.", groupPath + "/" + baseName);
    }

    for (const auto& marker : std::array<std::tuple<int, std::string, unsigned char, unsigned char>, 2>{
            std::tuple<int, std::string, unsigned char, unsigned char>{realId, realName, 1, 0},
            {imagId, imagName, 0, 1}}) {
        checkpointStatus = readByteVariableAttribute(groupId, std::get<0>(marker), std::get<1>(marker), "isComplex", 1, groupPath);
        if (!checkpointStatus) return checkpointStatus;
        checkpointStatus = readByteVariableAttribute(groupId, std::get<0>(marker), std::get<1>(marker), "isRealPart", std::get<2>(marker), groupPath);
        if (!checkpointStatus) return checkpointStatus;
        checkpointStatus = readByteVariableAttribute(groupId, std::get<0>(marker), std::get<1>(marker), "isImaginaryPart", std::get<3>(marker), groupPath);
        if (!checkpointStatus) return checkpointStatus;
    }

    if (Nkl > std::numeric_limits<std::size_t>::max() / Nj) return status(WVCheckpointStatusCode::shapeMismatch, "Coefficient shape exceeds addressable storage.", groupPath + "/" + baseName);
    if (output == nullptr) return WVCheckpointStatus::ok();
    const std::size_t count = Nj * Nkl;
    std::vector<double> real(count);
    std::vector<double> imag(count);
    if (timeSeries) {
        const std::size_t start[] = {selectedIndex, 0, 0};
        const std::size_t slab[] = {1, Nkl, Nj};
        result = nc_get_vara_double(groupId, realId, start, slab, real.data());
        if (result == NC_NOERR) result = nc_get_vara_double(groupId, imagId, start, slab, imag.data());
    } else {
        result = nc_get_var_double(groupId, realId, real.data());
        if (result == NC_NOERR) result = nc_get_var_double(groupId, imagId, imag.data());
    }
    if (result != NC_NOERR) return detail::netcdfFailure(result, "Complex coefficient read", groupPath + "/" + baseName);
    output->resize(count);
    for (std::size_t index = 0; index < count; ++index) (*output)[index] = {real[index], imag[index]};
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readCompactComplexCoefficient(
    int groupId, const std::string& groupPath,
    const std::string& baseName, std::size_t selectedIndex,
    std::size_t stateCount, std::size_t elementCount,
    std::vector<WVComplex64>* output) {
    const std::string realName = baseName + "_real";
    const std::string imagName = baseName + "_imag";
    int realId = -1;
    int imagId = -1;
    int result = nc_inq_varid(groupId, realName.c_str(), &realId);
    if (result == NC_ENOTVAR) {
        int plainId = -1;
        result = nc_inq_varid(groupId, baseName.c_str(), &plainId);
        if (result == NC_ENOTVAR)
            return status(WVCheckpointStatusCode::missingComplexPartner,
                          "Missing compact coefficient '" + baseName + "'.",
                          groupPath + "/" + baseName);
        if (result != NC_NOERR)
            return detail::netcdfFailure(result, "Coefficient lookup",
                                         groupPath + "/" + baseName);
        nc_type type = NC_NAT;
        std::vector<int> dimensions;
        auto checkpointStatus = inquireVariable(
            groupId, plainId, type, dimensions, groupPath + "/" + baseName);
        if (!checkpointStatus) return checkpointStatus;
        if (type != NC_DOUBLE || dimensions.size() != 1)
            return status(WVCheckpointStatusCode::shapeMismatch,
                          "A compact initial coefficient must use [kl].",
                          groupPath + "/" + baseName);
        std::string name;
        std::size_t length = 0;
        checkpointStatus = dimensionName(groupId, dimensions[0], name, length,
                                         groupPath + "/" + baseName);
        if (!checkpointStatus) return checkpointStatus;
        if (name != "kl" || length != elementCount)
            return status(WVCheckpointStatusCode::shapeMismatch,
                          "Compact coefficient kl dimension is incompatible.",
                          groupPath + "/" + baseName);
        if (output != nullptr) {
            std::vector<double> values(elementCount);
            result = nc_get_var_double(groupId, plainId, values.data());
            if (result != NC_NOERR)
                return detail::netcdfFailure(
                    result, "Compact coefficient read",
                    groupPath + "/" + baseName);
            output->resize(elementCount);
            for (std::size_t index = 0; index < elementCount; ++index)
                (*output)[index] = {values[index], 0.0};
        }
        return WVCheckpointStatus::ok();
    }
    if (result != NC_NOERR)
        return detail::netcdfFailure(result, "Complex-component lookup",
                                     groupPath + "/" + realName);
    result = nc_inq_varid(groupId, imagName.c_str(), &imagId);
    if (result == NC_ENOTVAR)
        return status(WVCheckpointStatusCode::missingComplexPartner,
                      "Missing compact imaginary component '" + imagName +
                          "'.",
                      groupPath + "/" + imagName);
    if (result != NC_NOERR)
        return detail::netcdfFailure(result, "Complex-component lookup",
                                     groupPath + "/" + imagName);
    nc_type realType = NC_NAT;
    nc_type imagType = NC_NAT;
    std::vector<int> realDimensions;
    std::vector<int> imagDimensions;
    auto checkpointStatus = inquireVariable(
        groupId, realId, realType, realDimensions, groupPath + "/" + realName);
    if (!checkpointStatus) return checkpointStatus;
    checkpointStatus = inquireVariable(
        groupId, imagId, imagType, imagDimensions, groupPath + "/" + imagName);
    if (!checkpointStatus) return checkpointStatus;
    if (realType != NC_DOUBLE || imagType != NC_DOUBLE ||
        realDimensions != imagDimensions ||
        (realDimensions.size() != 1 && realDimensions.size() != 2))
        return status(WVCheckpointStatusCode::shapeMismatch,
                      "Compact complex coefficients must use [kl] or [t,kl].",
                      groupPath + "/" + baseName);
    const bool timeSeries = realDimensions.size() == 2;
    for (std::size_t dimension = 0; dimension < realDimensions.size();
         ++dimension) {
        std::string name;
        std::size_t length = 0;
        checkpointStatus = dimensionName(
            groupId, realDimensions[dimension], name, length,
            groupPath + "/" + baseName);
        if (!checkpointStatus) return checkpointStatus;
        const std::string expectedName =
            timeSeries && dimension == 0 ? "t" : "kl";
        const std::size_t expectedLength =
            timeSeries && dimension == 0 ? stateCount : elementCount;
        if (name != expectedName || length != expectedLength)
            return status(WVCheckpointStatusCode::shapeMismatch,
                          "Compact coefficient dimensions are incompatible.",
                          groupPath + "/" + baseName);
    }
    for (const auto& marker :
         std::array<std::tuple<int, std::string, unsigned char,
                               unsigned char>, 2>{
             std::tuple<int, std::string, unsigned char, unsigned char>{
                 realId, realName, 1, 0},
             {imagId, imagName, 0, 1}}) {
        checkpointStatus = readByteVariableAttribute(
            groupId, std::get<0>(marker), std::get<1>(marker), "isComplex", 1,
            groupPath);
        if (!checkpointStatus) return checkpointStatus;
        checkpointStatus = readByteVariableAttribute(
            groupId, std::get<0>(marker), std::get<1>(marker), "isRealPart",
            std::get<2>(marker), groupPath);
        if (!checkpointStatus) return checkpointStatus;
        checkpointStatus = readByteVariableAttribute(
            groupId, std::get<0>(marker), std::get<1>(marker),
            "isImaginaryPart", std::get<3>(marker), groupPath);
        if (!checkpointStatus) return checkpointStatus;
    }
    if (output == nullptr) return WVCheckpointStatus::ok();
    std::vector<double> real(elementCount);
    std::vector<double> imag(elementCount);
    if (timeSeries) {
        const std::size_t start[] = {selectedIndex, 0};
        const std::size_t count[] = {1, elementCount};
        result = nc_get_vara_double(groupId, realId, start, count, real.data());
        if (result == NC_NOERR)
            result = nc_get_vara_double(groupId, imagId, start, count,
                                        imag.data());
    } else {
        result = nc_get_var_double(groupId, realId, real.data());
        if (result == NC_NOERR)
            result = nc_get_var_double(groupId, imagId, imag.data());
    }
    if (result != NC_NOERR)
        return detail::netcdfFailure(result, "Compact coefficient read",
                                     groupPath + "/" + baseName);
    output->resize(elementCount);
    for (std::size_t index = 0; index < elementCount; ++index)
        (*output)[index] = {real[index], imag[index]};
    return WVCheckpointStatus::ok();
}

bool forcingOrdinal(const std::string& name, std::size_t& ordinal) {
    static constexpr const char* prefix = "forcing-";
    if (name.rfind(prefix, 0) != 0) return false;
    const std::string digits = name.substr(std::char_traits<char>::length(prefix));
    if (digits.empty()) return false;
    const auto conversion = std::from_chars(digits.data(), digits.data() + digits.size(), ordinal);
    return conversion.ec == std::errc{} && conversion.ptr == digits.data() + digits.size() && ordinal > 0;
}

WVCheckpointStatus readForcingHeaders(
    const std::vector<GroupRecord>& groups,
    std::vector<WVCheckpointForcingHeader>& headers,
    std::vector<detail::WVForcingGroupSource>& sources) {
    const auto forcing = std::find_if(groups.begin(), groups.end(), [](const GroupRecord& group) { return group.path == "/forcing"; });
    if (forcing == groups.end()) return WVCheckpointStatus::ok();
    std::string singletonClass;
    bool hasSingletonClass = false;
    auto result = readOptionalTextAttribute(forcing->id, "AnnotatedClass", singletonClass, hasSingletonClass, forcing->path);
    if (!result) return result;
    std::vector<int> children;
    result = detail::childGroups(forcing->id, children, forcing->path);
    if (!result) return result;
    if (hasSingletonClass) {
        if (!children.empty()) return status(WVCheckpointStatusCode::invalidValue, "A singleton forcing group cannot also contain forcing-N records.", forcing->path);
        headers.push_back({1, forcing->path, singletonClass});
        sources.push_back({forcing->id, 1, forcing->path, singletonClass});
        return WVCheckpointStatus::ok();
    }
    for (const int child : children) {
        std::string name;
        result = detail::groupName(child, name, forcing->path);
        if (!result) return result;
        std::size_t ordinal = 0;
        if (!forcingOrdinal(name, ordinal)) continue;
        std::string annotatedClass;
        result = detail::readTextAttribute(child, "AnnotatedClass", annotatedClass, forcing->path + "/" + name);
        if (!result) return result;
        headers.push_back({ordinal, forcing->path + "/" + name, annotatedClass});
        sources.push_back({child, ordinal, forcing->path + "/" + name, annotatedClass});
    }
    std::sort(headers.begin(), headers.end(), [](const auto& left, const auto& right) { return left.ordinal < right.ordinal; });
    std::sort(sources.begin(), sources.end(), [](const auto& left, const auto& right) { return left.ordinal < right.ordinal; });
    for (std::size_t index = 0; index < headers.size(); ++index) {
        if (headers[index].ordinal != index + 1) return status(WVCheckpointStatusCode::invalidValue, "Forcing group ordinals must be contiguous and one-based.", "/forcing");
    }
    return WVCheckpointStatus::ok();
}

} // namespace

WVCoefficients WVCheckpointCoefficients::view() const noexcept {
    return {{Ap.data(), shape}, {Am.data(), shape}, {A0.data(), shape}};
}

WVState WVCheckpointState::view() const noexcept {
    return {t, t0, coefficients.view()};
}

namespace {

WVCheckpointStatus inspectOpenFile(
    int rootId,
    WVCheckpointStateSelection selection,
    const WVExtensionCatalog& catalog,
    WVCheckpointInspection& inspection,
    StateGroupRecord& stateGroup) {
    WVCheckpointInspection candidate;
    auto result = detail::readTextAttribute(rootId, "model_version", candidate.metadata.modelVersion, "/");
    if (!result) return result;
    result = validateVersion(candidate.metadata.modelVersion);
    if (!result) return result;

    std::string wvTransform;
    std::string annotatedClass;
    bool hasWVTransform = false;
    bool hasAnnotatedClass = false;
    result = readOptionalTextAttribute(rootId, "WVTransform", wvTransform, hasWVTransform, "/");
    if (!result) return result;
    result = readOptionalTextAttribute(rootId, "AnnotatedClass", annotatedClass, hasAnnotatedClass, "/");
    if (!result) return result;
    if (!hasWVTransform && !hasAnnotatedClass) return status(WVCheckpointStatusCode::missingAttribute, "Checkpoint root has neither WVTransform nor AnnotatedClass metadata.", "/");
    if (hasWVTransform && hasAnnotatedClass && wvTransform != annotatedClass) return status(WVCheckpointStatusCode::unsupportedTransform, "WVTransform and AnnotatedClass root metadata disagree.", "/");
    candidate.metadata.transformClass = hasWVTransform ? wvTransform : annotatedClass;
    const bool isConstant = candidate.metadata.transformClass ==
                            "WVTransformConstantStratification";
    const bool isBarotropicQG = candidate.metadata.transformClass ==
                                "WVTransformBarotropicQG";
    if (!isConstant && !isBarotropicQG)
        return status(WVCheckpointStatusCode::unsupportedTransform,
                      "The portable runtime profile does not support transform '" +
                          candidate.metadata.transformClass + "'.", "/");
    candidate.transformKind = isBarotropicQG
                                  ? WVPersistedTransformKind::barotropicQG
                                  : WVPersistedTransformKind::constantStratification;

    result = isConstant
                 ? readConfiguration(rootId, candidate.configuration)
                 : readBarotropicQGConfiguration(
                       rootId, candidate.barotropicQGConfiguration);
    if (!result) return result;
    result = detail::readDoubleScalar(rootId, "t0", candidate.t0, "/");
    if (!result) return result;
    if (!std::isfinite(candidate.t0)) return status(WVCheckpointStatusCode::invalidValue, "Checkpoint reference time t0 must be finite.", "/t0");

    std::vector<GroupRecord> groups;
    result = inspectGroupTree(rootId, "/", groups);
    if (!result) return result;
    result = isConstant ? findStateGroup(groups, stateGroup)
                        : findBarotropicQGStateGroup(groups, catalog,
                                                    stateGroup);
    if (!result) return result;
    candidate.metadata.stateGroupPath = stateGroup.path;
    result = inspectTime(stateGroup.id, stateGroup.path, selection, candidate.metadata.selectedStateIndex, candidate.metadata.stateCount, candidate.t);
    if (!result) return result;

    std::size_t Nkl = 0;
    result = detail::dimensionLength(stateGroup.id, "kl", Nkl, stateGroup.path);
    if (!result) return result;
    if (Nkl == 0)
        return status(WVCheckpointStatusCode::invalidValue,
                      "Checkpoint spectral dimensions must be nonempty.",
                      stateGroup.path);
    std::size_t coefficientCount = Nkl;
    if (isConstant) {
        std::size_t Nj = 0;
        result = detail::dimensionLength(stateGroup.id, "j", Nj,
                                         stateGroup.path);
        if (!result) return result;
        if (Nj == 0)
            return status(WVCheckpointStatusCode::invalidValue,
                          "Checkpoint j dimension must be nonempty.",
                          stateGroup.path);
        candidate.configuration.Nj = Nj;
        candidate.coefficientShape = {Nj, Nkl};
        candidate.stateDescription = {
            candidate.metadata.transformClass,
            {candidate.configuration.Nx, candidate.configuration.Ny,
             candidate.configuration.Nz},
            {{"Ap", {Nj, Nkl}, WVToleranceKind::coefficientEnergyScaled},
             {"Am", {Nj, Nkl}, WVToleranceKind::coefficientEnergyScaled},
             {"A0", {Nj, Nkl}, WVToleranceKind::coefficientEnergyScaled}},
            true};
        coefficientCount = Nj * Nkl;
        for (const char* family : {"Ap", "Am", "A0"}) {
            result = readComplexCoefficient(
                stateGroup.id, stateGroup.path, family,
                candidate.metadata.selectedStateIndex,
                candidate.metadata.stateCount, Nj, Nkl, nullptr);
            if (!result) return result;
        }
    } else {
        candidate.coefficientShape = {1, Nkl};
        candidate.stateDescription = {
            candidate.metadata.transformClass,
            {candidate.barotropicQGConfiguration.Nx,
             candidate.barotropicQGConfiguration.Ny},
            {{"A0", {Nkl}, WVToleranceKind::coefficientEnergyScaled}}, true};
        result = readCompactComplexCoefficient(
            stateGroup.id, stateGroup.path, "A0",
            candidate.metadata.selectedStateIndex,
            candidate.metadata.stateCount, Nkl, nullptr);
        if (!result) return result;
    }

    std::vector<detail::WVForcingGroupSource> forcingSources;
    result = readForcingHeaders(groups, candidate.metadata.forcingHeaders, forcingSources);
    if (!result) return result;
    result = isConstant
                 ? detail::decodeForcingSchedule(
                       forcingSources, candidate.configuration,
                       coefficientCount, catalog, candidate.forcingSchedule)
                 : detail::decodeForcingSchedule(
                       forcingSources, candidate.barotropicQGConfiguration,
                       coefficientCount, catalog, candidate.forcingSchedule);
    if (!result) return result;
    if (isConstant) {
        WVTransformConstantStratificationDescriptor descriptor;
        const auto descriptorStatus =
            WVTransformConstantStratificationDescriptor::create(
                candidate.configuration, descriptor);
        if (!descriptorStatus)
            return status(WVCheckpointStatusCode::descriptorFailure,
                          "Unable to rebuild the constant-stratification descriptor: " +
                              descriptorStatus.message,
                          stateGroup.path);
        if (descriptor.Nkl() != Nkl ||
            descriptor.spectralShape().rows !=
                candidate.coefficientShape.rows)
            return status(WVCheckpointStatusCode::shapeMismatch,
                          "Stored coefficients do not match the rebuilt descriptor.",
                          stateGroup.path);
    } else {
        WVTransformBarotropicQGDescriptor descriptor;
        const auto descriptorStatus =
            WVTransformBarotropicQGDescriptor::create(
                candidate.barotropicQGConfiguration, descriptor);
        if (!descriptorStatus || descriptor.Nkl() != Nkl)
            return status(WVCheckpointStatusCode::shapeMismatch,
                          "Stored compact A0 does not match the rebuilt "
                          "Barotropic QG descriptor.",
                          stateGroup.path);
    }
    inspection = std::move(candidate);
    return WVCheckpointStatus::ok();
}

} // namespace

WVCheckpointStatus WVCheckpointReader::inspect(const std::string& path, const WVExtensionCatalog& catalog, WVCheckpointInspection& inspection, WVCheckpointStateSelection selection) {
    WVNetCDFFile file;
    auto result = WVNetCDFFile::openReadOnly(path, file);
    if (!result) return result;
    StateGroupRecord stateGroup;
    return inspectOpenFile(file.id(), selection, catalog, inspection, stateGroup);
}

WVCheckpointStatus WVCheckpointReader::read(const std::string& path, const WVExtensionCatalog& catalog, WVCheckpoint& checkpoint, WVCheckpointStateSelection selection) {
    WVNetCDFFile file;
    auto result = WVNetCDFFile::openReadOnly(path, file);
    if (!result) return result;
    const int rootId = file.id();
    WVCheckpointInspection inspection;
    StateGroupRecord stateGroup;
    result = inspectOpenFile(rootId, selection, catalog, inspection, stateGroup);
    if (!result) return result;
    WVCheckpoint candidate;
    candidate.transformKind = inspection.transformKind;
    candidate.configuration = inspection.configuration;
    candidate.barotropicQGConfiguration =
        inspection.barotropicQGConfiguration;
    candidate.stateDescription = inspection.stateDescription;
    candidate.metadata = inspection.metadata;
    candidate.forcingSchedule = inspection.forcingSchedule;
    candidate.state.t = inspection.t;
    candidate.state.t0 = inspection.t0;
    if (inspection.transformKind ==
        WVPersistedTransformKind::constantStratification) {
        candidate.state.coefficients.shape = inspection.coefficientShape;
        const auto Nj = inspection.coefficientShape.rows;
        const auto Nkl = inspection.coefficientShape.columns;
        result = readComplexCoefficient(stateGroup.id, stateGroup.path, "Ap", candidate.metadata.selectedStateIndex, candidate.metadata.stateCount, Nj, Nkl, &candidate.state.coefficients.Ap);
        if (!result) return result;
        result = readComplexCoefficient(stateGroup.id, stateGroup.path, "Am", candidate.metadata.selectedStateIndex, candidate.metadata.stateCount, Nj, Nkl, &candidate.state.coefficients.Am);
        if (!result) return result;
        result = readComplexCoefficient(stateGroup.id, stateGroup.path, "A0", candidate.metadata.selectedStateIndex, candidate.metadata.stateCount, Nj, Nkl, &candidate.state.coefficients.A0);
        if (!result) return result;
    } else {
        candidate.transformState.transformIdentifier =
            inspection.stateDescription.transformIdentifier;
        candidate.transformState.spatialDimensions =
            inspection.stateDescription.spatialDimensions;
        candidate.transformState.t = inspection.t;
        candidate.transformState.t0 = inspection.t0;
        candidate.transformState.coefficientFamilies.push_back(
            {"A0", inspection.stateDescription.coefficientFamilies[0]
                       .spectralDimensions,
             {}});
        result = readCompactComplexCoefficient(
            stateGroup.id, stateGroup.path, "A0",
            candidate.metadata.selectedStateIndex,
            candidate.metadata.stateCount,
            inspection.coefficientShape.columns,
            &candidate.transformState.coefficientFamilies[0].values);
        if (!result) return result;
    }

    checkpoint = std::move(candidate);
    return WVCheckpointStatus::ok();
}

std::size_t checkpointCoefficientStorageBytes(
    const WVCheckpoint& checkpoint) noexcept {
    if (checkpoint.transformKind == WVPersistedTransformKind::barotropicQG) {
        std::size_t bytes = 0;
        for (const auto& family : checkpoint.transformState.coefficientFamilies)
            bytes += family.values.capacity() * sizeof(WVComplex64);
        return bytes;
    }
    return (checkpoint.state.coefficients.Ap.capacity() +
            checkpoint.state.coefficients.Am.capacity() +
            checkpoint.state.coefficients.A0.capacity()) *
           sizeof(WVComplex64);
}

} // namespace wavevortex::runtime
