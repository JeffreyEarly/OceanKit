#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"

#include <netcdf.h>

#include <cstddef>
#include <string>
#include <vector>

namespace wavevortex::runtime::detail {

class WVNetCDFFile final {
public:
    WVNetCDFFile() = default;
    ~WVNetCDFFile();
    WVNetCDFFile(const WVNetCDFFile&) = delete;
    WVNetCDFFile& operator=(const WVNetCDFFile&) = delete;
    WVNetCDFFile(WVNetCDFFile&& other) noexcept;
    WVNetCDFFile& operator=(WVNetCDFFile&& other) noexcept;

    static WVCheckpointStatus openReadOnly(const std::string& path, WVNetCDFFile& file);
    int id() const noexcept { return id_; }

private:
    explicit WVNetCDFFile(int id) noexcept : id_(id) {}
    int id_ = -1;
};

WVCheckpointStatus netcdfFailure(int status, const std::string& operation, const std::string& location);
WVCheckpointStatus readTextAttribute(int groupId, const std::string& name, std::string& value, const std::string& groupPath);
WVCheckpointStatus readDoubleScalar(int groupId, const std::string& name, double& value, const std::string& groupPath);
WVCheckpointStatus readLogicalScalar(int groupId, const std::string& name, bool& value, const std::string& groupPath);
WVCheckpointStatus readDoubleCoordinate(int groupId, const std::string& name, std::vector<double>& values, const std::string& groupPath);
WVCheckpointStatus childGroups(int groupId, std::vector<int>& groups, const std::string& groupPath);
WVCheckpointStatus groupName(int groupId, std::string& name, const std::string& parentPath);
WVCheckpointStatus variableIdIfPresent(int groupId, const std::string& name, int& variableId, bool& present, const std::string& groupPath);
WVCheckpointStatus dimensionLength(int groupId, const std::string& name, std::size_t& length, const std::string& groupPath);

} // namespace wavevortex::runtime::detail
