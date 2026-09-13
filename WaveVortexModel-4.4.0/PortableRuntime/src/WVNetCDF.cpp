#include "WVNetCDF.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace wavevortex::runtime::detail {

namespace {

WVCheckpointStatus missing(WVCheckpointStatusCode code, const std::string& kind, const std::string& name, const std::string& location) {
    return {code, "Missing " + kind + " '" + name + "' at " + location + ".", location};
}

} // namespace

WVNetCDFFile::~WVNetCDFFile() {
    if (id_ >= 0) {
        nc_close(id_);
        id_ = -1;
    }
}

WVNetCDFFile::WVNetCDFFile(WVNetCDFFile&& other) noexcept : id_(other.id_) {
    other.id_ = -1;
}

WVNetCDFFile& WVNetCDFFile::operator=(WVNetCDFFile&& other) noexcept {
    if (this != &other) {
        if (id_ >= 0) nc_close(id_);
        id_ = other.id_;
        other.id_ = -1;
    }
    return *this;
}

WVCheckpointStatus WVNetCDFFile::openReadOnly(const std::string& path, WVNetCDFFile& file) {
    int id = -1;
    const int status = nc_open(path.c_str(), NC_NOWRITE, &id);
    if (status != NC_NOERR) {
        return {WVCheckpointStatusCode::openFailure, "Unable to open NetCDF checkpoint '" + path + "': " + nc_strerror(status), path};
    }
    file = WVNetCDFFile(id);
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus netcdfFailure(int status, const std::string& operation, const std::string& location) {
    return {WVCheckpointStatusCode::netcdfFailure, operation + " failed at " + location + ": " + nc_strerror(status), location};
}

WVCheckpointStatus readTextAttribute(int groupId, const std::string& name, std::string& value, const std::string& groupPath) {
    nc_type type = NC_NAT;
    std::size_t length = 0;
    int status = nc_inq_att(groupId, NC_GLOBAL, name.c_str(), &type, &length);
    if (status == NC_ENOTATT) return missing(WVCheckpointStatusCode::missingAttribute, "attribute", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Attribute inspection", groupPath + "/@" + name);

    if (type == NC_CHAR) {
        std::string result(length, '\0');
        status = nc_get_att_text(groupId, NC_GLOBAL, name.c_str(), result.data());
        if (status != NC_NOERR) return netcdfFailure(status, "Attribute read", groupPath + "/@" + name);
        value = std::move(result);
        return WVCheckpointStatus::ok();
    }
    if (type == NC_STRING && length == 1) {
        char* result = nullptr;
        status = nc_get_att_string(groupId, NC_GLOBAL, name.c_str(), &result);
        if (status != NC_NOERR) return netcdfFailure(status, "Attribute read", groupPath + "/@" + name);
        value = result == nullptr ? std::string{} : std::string(result);
        nc_free_string(1, &result);
        return WVCheckpointStatus::ok();
    }
    return {WVCheckpointStatusCode::typeMismatch, "Attribute '" + name + "' at " + groupPath + " must be text.", groupPath + "/@" + name};
}

WVCheckpointStatus readDoubleScalar(int groupId, const std::string& name, double& value, const std::string& groupPath) {
    int variableId = -1;
    int status = nc_inq_varid(groupId, name.c_str(), &variableId);
    if (status == NC_ENOTVAR) return missing(WVCheckpointStatusCode::missingVariable, "variable", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable lookup", groupPath + "/" + name);
    nc_type type = NC_NAT;
    int dimensions = 0;
    status = nc_inq_var(groupId, variableId, nullptr, &type, &dimensions, nullptr, nullptr);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable inspection", groupPath + "/" + name);
    if (type != NC_DOUBLE || dimensions != 0) {
        return {WVCheckpointStatusCode::typeMismatch, "Variable '" + name + "' at " + groupPath + " must be a scalar double.", groupPath + "/" + name};
    }
    status = nc_get_var_double(groupId, variableId, &value);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable read", groupPath + "/" + name);
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readLogicalScalar(int groupId, const std::string& name, bool& value, const std::string& groupPath) {
    int variableId = -1;
    int status = nc_inq_varid(groupId, name.c_str(), &variableId);
    if (status == NC_ENOTVAR) return missing(WVCheckpointStatusCode::missingVariable, "variable", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable lookup", groupPath + "/" + name);
    nc_type type = NC_NAT;
    int dimensions = 0;
    status = nc_inq_var(groupId, variableId, nullptr, &type, &dimensions, nullptr, nullptr);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable inspection", groupPath + "/" + name);
    if ((type != NC_UBYTE && type != NC_BYTE) || dimensions != 0) {
        return {WVCheckpointStatusCode::typeMismatch, "Variable '" + name + "' at " + groupPath + " must be a scalar byte logical.", groupPath + "/" + name};
    }
    unsigned char raw = 0;
    status = nc_get_var_uchar(groupId, variableId, &raw);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable read", groupPath + "/" + name);
    if (raw > 1) return {WVCheckpointStatusCode::invalidValue, "Logical variable '" + name + "' at " + groupPath + " must contain 0 or 1.", groupPath + "/" + name};
    value = raw != 0;
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus readDoubleCoordinate(int groupId, const std::string& name, std::vector<double>& values, const std::string& groupPath) {
    int dimensionId = -1;
    int status = nc_inq_dimid(groupId, name.c_str(), &dimensionId);
    if (status == NC_EBADDIM) return missing(WVCheckpointStatusCode::missingDimension, "dimension", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Dimension lookup", groupPath + "/" + name);
    std::size_t length = 0;
    status = nc_inq_dimlen(groupId, dimensionId, &length);
    if (status != NC_NOERR) return netcdfFailure(status, "Dimension inspection", groupPath + "/" + name);
    if (length == 0) return {WVCheckpointStatusCode::invalidValue, "Coordinate '" + name + "' at " + groupPath + " must not be empty.", groupPath + "/" + name};

    int variableId = -1;
    status = nc_inq_varid(groupId, name.c_str(), &variableId);
    if (status == NC_ENOTVAR) return missing(WVCheckpointStatusCode::missingVariable, "coordinate variable", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable lookup", groupPath + "/" + name);
    nc_type type = NC_NAT;
    int dimensions = 0;
    status = nc_inq_vartype(groupId, variableId, &type);
    if (status == NC_NOERR) status = nc_inq_varndims(groupId, variableId, &dimensions);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable inspection", groupPath + "/" + name);
    if (type != NC_DOUBLE || dimensions != 1) {
        return {WVCheckpointStatusCode::shapeMismatch, "Coordinate variable '" + name + "' at " + groupPath + " must be a double over its matching dimension.", groupPath + "/" + name};
    }
    int variableDimensionId = -1;
    status = nc_inq_vardimid(groupId, variableId, &variableDimensionId);
    if (status != NC_NOERR) return netcdfFailure(status, "Variable inspection", groupPath + "/" + name);
    if (variableDimensionId != dimensionId) return {WVCheckpointStatusCode::shapeMismatch, "Coordinate variable '" + name + "' at " + groupPath + " must use its matching dimension.", groupPath + "/" + name};
    values.resize(length);
    status = nc_get_var_double(groupId, variableId, values.data());
    if (status != NC_NOERR) return netcdfFailure(status, "Coordinate read", groupPath + "/" + name);
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus childGroups(int groupId, std::vector<int>& groups, const std::string& groupPath) {
    int count = 0;
    int status = nc_inq_grps(groupId, &count, nullptr);
    if (status != NC_NOERR) return netcdfFailure(status, "Child-group inspection", groupPath);
    groups.resize(static_cast<std::size_t>(count));
    if (count > 0) {
        status = nc_inq_grps(groupId, nullptr, groups.data());
        if (status != NC_NOERR) return netcdfFailure(status, "Child-group inspection", groupPath);
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus groupName(int groupId, std::string& name, const std::string& parentPath) {
    char result[NC_MAX_NAME + 1] = {};
    const int status = nc_inq_grpname(groupId, result);
    if (status != NC_NOERR) return netcdfFailure(status, "Group-name read", parentPath);
    name = result;
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus variableIdIfPresent(int groupId, const std::string& name, int& variableId, bool& present, const std::string& groupPath) {
    const int status = nc_inq_varid(groupId, name.c_str(), &variableId);
    if (status == NC_ENOTVAR) {
        variableId = -1;
        present = false;
        return WVCheckpointStatus::ok();
    }
    if (status != NC_NOERR) return netcdfFailure(status, "Variable lookup", groupPath + "/" + name);
    present = true;
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus dimensionLength(int groupId, const std::string& name, std::size_t& length, const std::string& groupPath) {
    int dimensionId = -1;
    int status = nc_inq_dimid(groupId, name.c_str(), &dimensionId);
    if (status == NC_EBADDIM) return missing(WVCheckpointStatusCode::missingDimension, "dimension", name, groupPath);
    if (status != NC_NOERR) return netcdfFailure(status, "Dimension lookup", groupPath + "/" + name);
    status = nc_inq_dimlen(groupId, dimensionId, &length);
    if (status != NC_NOERR) return netcdfFailure(status, "Dimension inspection", groupPath + "/" + name);
    return WVCheckpointStatus::ok();
}

} // namespace wavevortex::runtime::detail
