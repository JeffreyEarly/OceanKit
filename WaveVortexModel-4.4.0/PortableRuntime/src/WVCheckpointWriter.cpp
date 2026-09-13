#include "WaveVortexRuntime/WVCheckpointWriter.hpp"
#include "WaveVortexRuntime/WVForcingContracts.hpp"

#include "WVCheckpointWriterTestHooks.hpp"
#include "WVNetCDF.hpp"

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace wavevortex::runtime {

namespace {

using detail::WVCheckpointWriterFailurePoint;

std::atomic<WVCheckpointWriterFailurePoint> injectedFailure{WVCheckpointWriterFailurePoint::none};
std::atomic<std::uint64_t> temporarySequence{0};

WVCheckpointStatus status(WVCheckpointStatusCode code, std::string message, std::string location) {
    return {code, std::move(message), std::move(location)};
}

bool sameComplex(const std::vector<WVComplex64>& left, const std::vector<WVComplex64>& right) noexcept {
    if (left.size() != right.size()) return false;
    for (std::size_t index = 0; index < left.size(); ++index) {
        const bool sameReal = left[index].real == right[index].real || (std::isnan(left[index].real) && std::isnan(right[index].real));
        const bool sameImag = left[index].imag == right[index].imag || (std::isnan(left[index].imag) && std::isnan(right[index].imag));
        if (!sameReal || !sameImag) return false;
    }
    return true;
}

bool compatibleModelVersion(const std::string& version) {
    const auto separator = version.find('.');
    return !version.empty() && (separator == std::string::npos ? version : version.substr(0, separator)) == "4";
}

WVCheckpointStatus validateEntry(
    const WVFrozenForcingEntry& entry,
    const WVForcingCatalog& catalog,
    const WVTransformConstantStratificationConfiguration&,
    std::size_t) {
    const auto* registration = catalog.registration(entry.typeIdentifier, entry.contractVersion);
    if (registration == nullptr || !registration->isSupported || !registration->factory || registration->contractVersion != entry.contractVersion) {
        return status(WVCheckpointStatusCode::unsupportedForcing, "The checkpoint contains a forcing entry that portable runtime v1 cannot write.", "/forcing/@AnnotatedClass");
    }
    if (entry.name.empty()) return status(WVCheckpointStatusCode::malformedForcing, "Forcing names must not be empty.", "/forcing/@name");
    if (entry.stage != registration->stage) return status(WVCheckpointStatusCode::malformedForcing, "Forcing stage does not match its registered implementation.", "/forcing");
    const auto validation = catalog.validateConfiguration(entry);
    return validation ? WVCheckpointStatus::ok() : status(WVCheckpointStatusCode::malformedForcing,validation.message,"/forcing");
}

WVCheckpointStatus validateCheckpoint(const WVCheckpoint& checkpoint, const WVForcingCatalog& catalog, std::vector<const WVFrozenForcingEntry*>& sourceOrder) {
    if (checkpoint.metadata.profileIdentifier != WVCheckpointProfileIdentifier || checkpoint.metadata.profileVersion != WVCheckpointProfileVersion) return status(WVCheckpointStatusCode::invalidValue, "Unsupported checkpoint profile.", "/");
    if (!compatibleModelVersion(checkpoint.metadata.modelVersion)) return status(WVCheckpointStatusCode::unsupportedModelVersion, "The writer accepts only WaveVortexModel 4.x checkpoints.", "/@model_version");
    if (checkpoint.metadata.transformClass != "WVTransformConstantStratification") return status(WVCheckpointStatusCode::unsupportedTransform, "The writer supports only WVTransformConstantStratification.", "/");
    if (checkpoint.forcingSchedule.profileIdentifier != WVForcingScheduleProfileIdentifier || checkpoint.forcingSchedule.profileVersion != WVForcingScheduleProfileVersion) return status(WVCheckpointStatusCode::unsupportedForcing, "Unsupported frozen forcing schedule profile.", "/forcing");
    if (!std::isfinite(checkpoint.state.t) || !std::isfinite(checkpoint.state.t0)) return status(WVCheckpointStatusCode::invalidValue, "Checkpoint times must be finite.", "/t");

    WVTransformConstantStratificationDescriptor descriptor;
    const auto descriptorStatus = WVTransformConstantStratificationDescriptor::create(checkpoint.configuration, descriptor);
    if (!descriptorStatus) return status(WVCheckpointStatusCode::descriptorFailure, "Unable to validate the checkpoint descriptor: " + descriptorStatus.message, "/");
    const auto shape = checkpoint.state.coefficients.shape;
    if (shape.rows != checkpoint.configuration.Nj || shape.columns != descriptor.Nkl()) return status(WVCheckpointStatusCode::shapeMismatch, "Checkpoint coefficient shape does not match the rebuilt descriptor.", "/");
    const std::size_t coefficientCount = shape.elementCount();
    for (const auto* values : {&checkpoint.state.coefficients.Ap, &checkpoint.state.coefficients.Am, &checkpoint.state.coefficients.A0}) {
        if (values->size() != coefficientCount) return status(WVCheckpointStatusCode::shapeMismatch, "Checkpoint coefficient arrays have incompatible lengths.", "/");
    }

    sourceOrder.clear();
    sourceOrder.reserve(checkpoint.forcingSchedule.entries.size());
    for (const auto& entry : checkpoint.forcingSchedule.entries) sourceOrder.push_back(&entry);
    std::sort(sourceOrder.begin(), sourceOrder.end(), [](const auto* left, const auto* right) { return left->ordinal < right->ordinal; });
    std::unordered_set<std::string> names;
    for (std::size_t index = 0; index < sourceOrder.size(); ++index) {
        const auto& entry = *sourceOrder[index];
        if (entry.ordinal != index + 1) return status(WVCheckpointStatusCode::invalidValue, "Forcing ordinals must be contiguous and one-based.", "/forcing");
        if (!names.insert(entry.name).second) return status(WVCheckpointStatusCode::duplicateForcing, "Forcing names must be unique.", "/forcing/@name");
        const auto result = validateEntry(entry, catalog, checkpoint.configuration, coefficientCount);
        if (!result) return result;
    }
    return WVCheckpointStatus::ok();
}

class WritableNetCDFFile final {
public:
    WritableNetCDFFile() = default;
    ~WritableNetCDFFile() { if (id_ >= 0) nc_close(id_); }
    WritableNetCDFFile(const WritableNetCDFFile&) = delete;
    WritableNetCDFFile& operator=(const WritableNetCDFFile&) = delete;

    static WVCheckpointStatus create(const std::string& path, WritableNetCDFFile& output) {
        int id = -1;
        const int result = nc_create(path.c_str(), NC_NETCDF4 | NC_NOCLOBBER, &id);
        if (result != NC_NOERR) return status(WVCheckpointStatusCode::writeFailure, "Unable to create temporary checkpoint: " + std::string(nc_strerror(result)), path);
        output.id_ = id;
        return WVCheckpointStatus::ok();
    }

    int id() const noexcept { return id_; }
    WVCheckpointStatus sync(const std::string& path) {
        const int result = nc_sync(id_);
        return result == NC_NOERR ? WVCheckpointStatus::ok() : detail::netcdfFailure(result, "Checkpoint synchronization", path);
    }
    WVCheckpointStatus close(const std::string& path) {
        const int id = std::exchange(id_, -1);
        const int result = nc_close(id);
        return result == NC_NOERR ? WVCheckpointStatus::ok() : detail::netcdfFailure(result, "Checkpoint close", path);
    }

private:
    int id_ = -1;
};

WVCheckpointStatus checked(int netcdfStatus, const std::string& operation, const std::string& location) {
    return netcdfStatus == NC_NOERR ? WVCheckpointStatus::ok() : detail::netcdfFailure(netcdfStatus, operation, location);
}

WVCheckpointStatus textAttribute(int group, int variable, const char* name, const std::string& value, const std::string& location) {
    return checked(nc_put_att_text(group, variable, name, value.size(), value.data()), "Text-attribute definition", location + "/@" + name);
}

WVCheckpointStatus byteAttribute(int group, int variable, const char* name, unsigned char value, const std::string& location) {
    return checked(nc_put_att_uchar(group, variable, name, NC_UBYTE, 1, &value), "Byte-attribute definition", location + "/@" + name);
}

WVCheckpointStatus defineDouble(int group, const char* name, const std::vector<int>& dimensions, int& variable, const std::string& path) {
    return checked(nc_def_var(group, name, NC_DOUBLE, static_cast<int>(dimensions.size()), dimensions.empty() ? nullptr : dimensions.data(), &variable), "Variable definition", path + "/" + name);
}

WVCheckpointStatus defineLogical(int group, const char* name, int& variable, const std::string& path) {
    auto result = checked(nc_def_var(group, name, NC_UBYTE, 0, nullptr, &variable), "Variable definition", path + "/" + name);
    if (!result) return result;
    return byteAttribute(group, variable, "isLogicalType", 1, path + "/" + name);
}

WVCheckpointStatus markComplex(int group, int variable, bool realPart, const std::string& path) {
    for (const auto& marker : std::array<std::pair<const char*, unsigned char>, 3>{{
            {"isComplex", 1}, {"isRealPart", static_cast<unsigned char>(realPart)}, {"isImaginaryPart", static_cast<unsigned char>(!realPart)}}}) {
        const auto result = byteAttribute(group, variable, marker.first, marker.second, path);
        if (!result) return result;
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus defineComplexPair(int group, const std::string& base, const std::vector<int>& dimensions, const std::string& path) {
    int real = -1;
    int imag = -1;
    auto result = defineDouble(group, (base + "_real").c_str(), dimensions, real, path);
    if (!result) return result;
    result = defineDouble(group, (base + "_imag").c_str(), dimensions, imag, path);
    if (!result) return result;
    result = markComplex(group, real, true, path + "/" + base + "_real");
    if (!result) return result;
    return markComplex(group, imag, false, path + "/" + base + "_imag");
}

WVCheckpointStatus defineRoot(int root, const WVCheckpoint& checkpoint, std::array<int, 5>& dimensions) {
    const auto& configuration = checkpoint.configuration;
    const std::array<std::pair<const char*, std::size_t>, 5> definitions = {{{"j", configuration.Nj}, {"kl", checkpoint.state.coefficients.shape.columns}, {"x", configuration.Nx}, {"y", configuration.Ny}, {"z", configuration.Nz}}};
    for (std::size_t index = 0; index < definitions.size(); ++index) {
        auto result = checked(nc_def_dim(root, definitions[index].first, definitions[index].second, &dimensions[index]), "Dimension definition", "/" + std::string(definitions[index].first));
        if (!result) return result;
        int variable = -1;
        result = defineDouble(root, definitions[index].first, {dimensions[index]}, variable, "/");
        if (!result) return result;
    }
    for (const char* name : {"Lx", "Ly", "Lz", "N0", "g", "latitude", "planetaryRadius", "rho0", "rotationRate", "t", "t0"}) {
        int variable = -1;
        const auto result = defineDouble(root, name, {}, variable, "/");
        if (!result) return result;
    }
    for (const char* name : {"isHydrostatic", "shouldAntialias"}) {
        int variable = -1;
        const auto result = defineLogical(root, name, variable, "/");
        if (!result) return result;
    }
    for (const char* base : {"Ap", "Am", "A0"}) {
        const auto result = defineComplexPair(root, base, {dimensions[1], dimensions[0]}, "/");
        if (!result) return result;
    }
    auto result = textAttribute(root, NC_GLOBAL, "AnnotatedClass", "WVTransformConstantStratification", "/");
    if (!result) return result;
    result = textAttribute(root, NC_GLOBAL, "WVTransform", "WVTransformConstantStratification", "/");
    if (!result) return result;
    result = textAttribute(root, NC_GLOBAL, "model_version", checkpoint.metadata.modelVersion, "/");
    if (!result) return result;
    result = textAttribute(root, NC_GLOBAL, "source", "Created with the WaveVortex portable runtime", "/");
    if (!result) return result;
    result = textAttribute(root, NC_GLOBAL, "history", "Transactional root-level checkpoint written by the WaveVortex portable runtime.", "/");
    if (!result) return result;
    const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif
    std::ostringstream date;
    date << std::put_time(&utc, "%Y-%m-%d %H:%M:%S UTC");
    return textAttribute(root, NC_GLOBAL, "date_created", date.str(), "/");
}

WVCheckpointStatus defineForcingEntry(int group, const WVFrozenForcingEntry& entry, const WVForcingCatalog& catalog, const std::array<int, 5>& rootDimensions, const std::string& path) {
    auto result = textAttribute(group, NC_GLOBAL, "AnnotatedClass", entry.typeIdentifier, path);
    if (!result) return result;
    const auto* registration = catalog.registration(entry.typeIdentifier, entry.contractVersion);
    if (registration == nullptr || !registration->isSupported) return status(WVCheckpointStatusCode::unsupportedForcing,"Unsupported forcing reached checkpoint definition.",path);
    if (registration->persistence.writesNameAttribute) {
        result = textAttribute(group,NC_GLOBAL,"name",entry.name,path);
        if (!result) return result;
    }
    for (const auto& field : registration->persistence.fields) {
        const auto* value = entry.configuration.value(field.recordName);
        if (value == nullptr && field.optional) continue;
        if (value == nullptr) return status(WVCheckpointStatusCode::malformedForcing,"Required forcing configuration value is missing.",path+"/"+field.netcdfName);
        if (field.encoding == WVForcingPersistenceEncoding::textAttribute) {
            const auto* values = std::get_if<std::vector<std::string>>(&value->storage);
            if (values == nullptr || values->size() != 1) return status(WVCheckpointStatusCode::malformedForcing,"Forcing text attribute is malformed.",path);
            if (!values->front().empty()) { result = textAttribute(group,NC_GLOBAL,field.netcdfName.c_str(),values->front(),path); if (!result) return result; }
            continue;
        }
        std::vector<int> dimensions;
        if (field.dimensions == WVForcingDimensionRule::ownLength) {
            if (value->valueCount() == 0) continue;
            int dimension = -1;
            result = checked(nc_def_dim(group,field.netcdfName.c_str(),value->valueCount(),&dimension),"Forcing dimension definition",path+"/"+field.netcdfName);
            if (!result) return result;
            dimensions = {dimension};
        } else if (field.dimensions == WVForcingDimensionRule::referencedLength) {
            int dimension = -1;
            const int inquiry = nc_inq_dimid(group,field.dimensionReference.c_str(),&dimension);
            if (inquiry == NC_EBADDIM && field.optional) continue;
            result = checked(inquiry,"Forcing dimension lookup",path+"/"+field.dimensionReference);
            if (!result) return result;
            dimensions = {dimension};
        } else if (field.dimensions == WVForcingDimensionRule::horizontalYX) {
            dimensions = {rootDimensions[3],rootDimensions[2]};
        } else if (field.dimensions == WVForcingDimensionRule::componentPair) {
            int dimension = -1;
            result = checked(nc_def_dim(group,"barotropicVelocityComponent",2,&dimension),"Forcing dimension definition",path+"/barotropicVelocityComponent");
            if (!result) return result;
            int coordinate = -1;
            result = defineDouble(group,"barotropicVelocityComponent",{dimension},coordinate,path);
            if (!result) return result;
            dimensions = {dimension};
        }
        int variable = -1;
        if (field.encoding == WVForcingPersistenceEncoding::realVariable) result = defineDouble(group,field.netcdfName.c_str(),dimensions,variable,path);
        else if (field.encoding == WVForcingPersistenceEncoding::logicalVariable) result = defineLogical(group,field.netcdfName.c_str(),variable,path);
        else if (field.encoding == WVForcingPersistenceEncoding::zeroBasedIndexVariable) result = checked(nc_def_var(group,field.netcdfName.c_str(),NC_UINT64,static_cast<int>(dimensions.size()),dimensions.data(),&variable),"Forcing variable definition",path+"/"+field.netcdfName);
        else if (field.encoding == WVForcingPersistenceEncoding::complexVariable) result = defineComplexPair(group,field.netcdfName,dimensions,path);
        if (!result) return result;
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus variableId(int group, const std::string& name, int& variable, const std::string& path) {
    return checked(nc_inq_varid(group, name.c_str(), &variable), "Variable lookup", path + "/" + name);
}

WVCheckpointStatus writeDouble(int group, const std::string& name, double value, const std::string& path) {
    int variable = -1;
    auto result = variableId(group, name, variable, path);
    if (!result) return result;
    return checked(nc_put_var_double(group, variable, &value), "Variable write", path + "/" + name);
}

WVCheckpointStatus writeLogical(int group, const std::string& name, bool value, const std::string& path) {
    int variable = -1;
    auto result = variableId(group, name, variable, path);
    if (!result) return result;
    const unsigned char raw = value ? 1 : 0;
    return checked(nc_put_var_uchar(group, variable, &raw), "Variable write", path + "/" + name);
}

WVCheckpointStatus writeDoubles(int group, const std::string& name, const std::vector<double>& values, const std::string& path) {
    int variable = -1;
    auto result = variableId(group, name, variable, path);
    if (!result) return result;
    return checked(nc_put_var_double(group, variable, values.data()), "Variable write", path + "/" + name);
}

WVCheckpointStatus writeComplexPair(int group, const std::string& base, const std::vector<WVComplex64>& values, const std::string& path) {
    std::vector<double> real(values.size());
    std::vector<double> imag(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) {
        real[index] = values[index].real;
        imag[index] = values[index].imag;
    }
    auto result = writeDoubles(group, base + "_real", real, path);
    if (!result) return result;
    return writeDoubles(group, base + "_imag", imag, path);
}

WVCheckpointStatus writeRoot(int root, const WVCheckpoint& checkpoint) {
    const auto& configuration = checkpoint.configuration;
    std::vector<double> j(configuration.Nj), kl(checkpoint.state.coefficients.shape.columns), x(configuration.Nx), y(configuration.Ny), z(configuration.Nz);
    for (std::size_t index = 0; index < j.size(); ++index) j[index] = static_cast<double>(index);
    for (std::size_t index = 0; index < kl.size(); ++index) kl[index] = static_cast<double>(index);
    for (std::size_t index = 0; index < x.size(); ++index) x[index] = static_cast<double>(index) * configuration.Lx / static_cast<double>(configuration.Nx);
    for (std::size_t index = 0; index < y.size(); ++index) y[index] = static_cast<double>(index) * configuration.Ly / static_cast<double>(configuration.Ny);
    for (std::size_t index = 0; index < z.size(); ++index) z[index] = -configuration.Lz + static_cast<double>(index) * configuration.Lz / static_cast<double>(configuration.Nz - 1);
    for (const auto& coordinate : std::array<std::pair<const char*, const std::vector<double>*>, 5>{{{"j", &j}, {"kl", &kl}, {"x", &x}, {"y", &y}, {"z", &z}}}) {
        const auto result = writeDoubles(root, coordinate.first, *coordinate.second, "/");
        if (!result) return result;
    }
    for (const auto& scalar : std::array<std::pair<const char*, double>, 11>{{
            {"Lx", configuration.Lx}, {"Ly", configuration.Ly}, {"Lz", configuration.Lz}, {"N0", configuration.N0}, {"g", configuration.g},
            {"latitude", configuration.latitude}, {"planetaryRadius", configuration.planetaryRadius}, {"rho0", configuration.rho0}, {"rotationRate", configuration.rotationRate},
            {"t", checkpoint.state.t}, {"t0", checkpoint.state.t0}}}) {
        const auto result = writeDouble(root, scalar.first, scalar.second, "/");
        if (!result) return result;
    }
    auto result = writeLogical(root, "isHydrostatic", configuration.isHydrostatic, "/");
    if (!result) return result;
    result = writeLogical(root, "shouldAntialias", configuration.shouldAntialias, "/");
    if (!result) return result;
    result = writeComplexPair(root, "Ap", checkpoint.state.coefficients.Ap, "/");
    if (!result) return result;
    result = writeComplexPair(root, "Am", checkpoint.state.coefficients.Am, "/");
    if (!result) return result;
    return writeComplexPair(root, "A0", checkpoint.state.coefficients.A0, "/");
}

WVCheckpointStatus writeForcingEntry(int group, const WVFrozenForcingEntry& entry, const WVForcingCatalog& catalog, const std::string& path) {
    const auto* registration = catalog.registration(entry.typeIdentifier, entry.contractVersion);
    if (registration == nullptr || !registration->isSupported) return status(WVCheckpointStatusCode::unsupportedForcing,"Unsupported forcing reached checkpoint writing.",path);
    for (const auto& field : registration->persistence.fields) {
        const auto* value = entry.configuration.value(field.recordName);
        if (value == nullptr && field.optional) continue;
        if (value == nullptr) return status(WVCheckpointStatusCode::malformedForcing,"Required forcing value is missing.",path);
        WVCheckpointStatus result = WVCheckpointStatus::ok();
        if (field.encoding == WVForcingPersistenceEncoding::realVariable) {
            const auto* values = std::get_if<std::vector<double>>(&value->storage);
            if (values == nullptr) return status(WVCheckpointStatusCode::malformedForcing,"Forcing real value has the wrong type.",path);
            result = writeDoubles(group,field.netcdfName,*values,path);
        } else if (field.encoding == WVForcingPersistenceEncoding::logicalVariable) {
            const auto* values = std::get_if<std::vector<std::uint8_t>>(&value->storage);
            if (values == nullptr || values->size() != 1) return status(WVCheckpointStatusCode::malformedForcing,"Forcing logical value is malformed.",path);
            result = writeLogical(group,field.netcdfName,values->front() != 0,path);
        } else if (field.encoding == WVForcingPersistenceEncoding::zeroBasedIndexVariable) {
            const auto* values = std::get_if<std::vector<std::int64_t>>(&value->storage);
            if (values == nullptr) return status(WVCheckpointStatusCode::malformedForcing,"Forcing index value has the wrong type.",path);
            if (values->empty()) continue;
            std::vector<unsigned long long> oneBased(values->size());
            for (std::size_t index = 0; index < values->size(); ++index) oneBased[index] = static_cast<unsigned long long>((*values)[index]+1);
            int variable = -1;
            result = variableId(group,field.netcdfName,variable,path);
            if (result) result = checked(nc_put_var_ulonglong(group,variable,oneBased.data()),"Forcing index write",path+"/"+field.netcdfName);
        } else if (field.encoding == WVForcingPersistenceEncoding::complexVariable) {
            const auto* real = std::get_if<std::vector<double>>(&value->storage);
            const auto* imaginaryValue = entry.configuration.value(field.imaginaryRecordName);
            const auto* imag = imaginaryValue == nullptr ? nullptr : std::get_if<std::vector<double>>(&imaginaryValue->storage);
            if (real == nullptr || imag == nullptr || real->size() != imag->size()) return status(WVCheckpointStatusCode::malformedForcing,"Forcing complex value is malformed.",path);
            if (field.dimensions == WVForcingDimensionRule::referencedLength && real->empty()) continue;
            if (field.dimensions == WVForcingDimensionRule::componentPair) { result = writeDoubles(group,"barotropicVelocityComponent",{1.0,2.0},path); if (!result) return result; }
            result = writeDoubles(group,field.netcdfName+"_real",*real,path);
            if (result) result = writeDoubles(group,field.netcdfName+"_imag",*imag,path);
        }
        if (!result) return result;
    }
    return WVCheckpointStatus::ok();
}

WVCheckpointStatus verifyTemporary(const std::string& path, const WVExtensionCatalog& catalog, const WVCheckpoint& expected) {
    WVCheckpoint actual;
    auto result = WVCheckpointReader::read(path, catalog, actual);
    if (!result) return status(WVCheckpointStatusCode::writeFailure, "The temporary checkpoint failed structural validation: " + result.message, result.location);
    if (!sameTransformConfiguration(expected.configuration, actual.configuration) || expected.state.t != actual.state.t || expected.state.t0 != actual.state.t0 || expected.state.coefficients.shape.rows != actual.state.coefficients.shape.rows || expected.state.coefficients.shape.columns != actual.state.coefficients.shape.columns || !sameComplex(expected.state.coefficients.Ap, actual.state.coefficients.Ap) || !sameComplex(expected.state.coefficients.Am, actual.state.coefficients.Am) || !sameComplex(expected.state.coefficients.A0, actual.state.coefficients.A0) || expected.forcingSchedule.entries.size() != actual.forcingSchedule.entries.size()) {
        return status(WVCheckpointStatusCode::writeFailure, "The temporary checkpoint did not reproduce the requested checkpoint state.", path);
    }
    return WVCheckpointStatus::ok();
}

std::filesystem::path temporaryPathFor(const std::filesystem::path& destination) {
    const auto parent = destination.has_parent_path() ? destination.parent_path() : std::filesystem::current_path();
    const auto sequence = temporarySequence.fetch_add(1, std::memory_order_relaxed);
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return parent / (destination.filename().string() + ".tmp-" + std::to_string(stamp) + "-" + std::to_string(sequence));
}

class TemporaryPath final {
public:
    explicit TemporaryPath(std::filesystem::path value) : value_(std::move(value)) {}
    ~TemporaryPath() { if (!committed_) { std::error_code ignored; std::filesystem::remove(value_, ignored); } }
    const std::filesystem::path& value() const noexcept { return value_; }
    void committed() noexcept { committed_ = true; }
private:
    std::filesystem::path value_;
    bool committed_ = false;
};

WVCheckpointStatus writeTemporary(
    const std::filesystem::path& path,
    const WVForcingCatalog& catalog,
    const WVCheckpoint& checkpoint,
    const std::vector<const WVFrozenForcingEntry*>& sourceOrder) {
    WritableNetCDFFile file;
    auto result = WritableNetCDFFile::create(path.string(), file);
    if (!result) return result;
    std::array<int, 5> rootDimensions{};
    result = defineRoot(file.id(), checkpoint, rootDimensions);
    if (!result) return result;

    std::vector<int> forcingGroups;
    forcingGroups.reserve(sourceOrder.size());
    if (!sourceOrder.empty()) {
        int forcingRoot = -1;
        result = checked(nc_def_grp(file.id(), "forcing", &forcingRoot), "Forcing-group definition", "/forcing");
        if (!result) return result;
        if (sourceOrder.size() == 1) {
            forcingGroups.push_back(forcingRoot);
            result = defineForcingEntry(forcingRoot, *sourceOrder.front(), catalog, rootDimensions, "/forcing");
            if (!result) return result;
        } else {
            for (std::size_t index = 0; index < sourceOrder.size(); ++index) {
                int group = -1;
                const std::string name = "forcing-" + std::to_string(index + 1);
                result = checked(nc_def_grp(forcingRoot, name.c_str(), &group), "Forcing-group definition", "/forcing/" + name);
                if (!result) return result;
                forcingGroups.push_back(group);
                result = defineForcingEntry(group, *sourceOrder[index], catalog, rootDimensions, "/forcing/" + name);
                if (!result) return result;
            }
        }
    }
    if (injectedFailure.load(std::memory_order_relaxed) == WVCheckpointWriterFailurePoint::afterDefinition) return status(WVCheckpointStatusCode::writeFailure, "Injected checkpoint failure after definition.", path.string());
    result = checked(nc_enddef(file.id()), "End checkpoint definition", path.string());
    if (!result) return result;
    result = writeRoot(file.id(), checkpoint);
    if (!result) return result;
    for (std::size_t index = 0; index < sourceOrder.size(); ++index) {
        const std::string groupPath = sourceOrder.size() == 1 ? "/forcing" : "/forcing/forcing-" + std::to_string(index + 1);
        result = writeForcingEntry(forcingGroups[index], *sourceOrder[index], catalog, groupPath);
        if (!result) return result;
    }
    if (injectedFailure.load(std::memory_order_relaxed) == WVCheckpointWriterFailurePoint::afterWrite) return status(WVCheckpointStatusCode::writeFailure, "Injected checkpoint failure after writing.", path.string());
    result = file.sync(path.string());
    if (!result) return result;
    return file.close(path.string());
}

} // namespace

namespace detail {

void setCheckpointWriterFailurePoint(WVCheckpointWriterFailurePoint point) noexcept {
    injectedFailure.store(point, std::memory_order_relaxed);
}

} // namespace detail

WVCheckpointStatus WVCheckpointWriter::write(const std::string& path, const WVExtensionCatalog& catalog, const WVCheckpoint& checkpoint, WVCheckpointCommitPolicy commitPolicy) {
    if (path.empty()) return status(WVCheckpointStatusCode::writeFailure, "Checkpoint destination path must not be empty.", path);
    try {
        std::vector<const WVFrozenForcingEntry*> sourceOrder;
        auto result = validateCheckpoint(checkpoint, catalog.forcings(), sourceOrder);
        if (!result) return result;

        const std::filesystem::path destination = std::filesystem::absolute(path).lexically_normal();
        const auto parent = destination.parent_path();
        std::error_code filesystemError;
        if (!std::filesystem::exists(parent, filesystemError) || filesystemError) return status(WVCheckpointStatusCode::writeFailure, "Checkpoint destination parent does not exist.", parent.string());
        if (std::filesystem::is_directory(destination, filesystemError)) return status(WVCheckpointStatusCode::writeFailure, "Checkpoint destination is a directory.", destination.string());
        const auto destinationStatus = std::filesystem::symlink_status(destination, filesystemError);
        if (filesystemError == std::errc::no_such_file_or_directory) filesystemError.clear();
        if (commitPolicy == WVCheckpointCommitPolicy::createNew && filesystemError) return status(WVCheckpointStatusCode::commitFailure, "Unable to inspect create-new checkpoint destination: " + filesystemError.message(), destination.string());
        if (commitPolicy == WVCheckpointCommitPolicy::createNew && destinationStatus.type() != std::filesystem::file_type::not_found) return status(WVCheckpointStatusCode::commitFailure, "Checkpoint destination already exists.", destination.string());
        filesystemError.clear();

        TemporaryPath temporary(temporaryPathFor(destination));
        result = writeTemporary(temporary.value(), catalog.forcings(), checkpoint, sourceOrder);
        if (!result) return result;
        result = verifyTemporary(temporary.value().string(), catalog, checkpoint);
        if (!result) return result;
        if (injectedFailure.load(std::memory_order_relaxed) == WVCheckpointWriterFailurePoint::beforeCommit) return status(WVCheckpointStatusCode::commitFailure, "Injected checkpoint failure before commit.", destination.string());

        if (commitPolicy == WVCheckpointCommitPolicy::createNew) {
            std::filesystem::create_hard_link(temporary.value(), destination, filesystemError);
            if (filesystemError) return status(WVCheckpointStatusCode::commitFailure, "Unable to atomically create checkpoint without replacement: " + filesystemError.message(), destination.string());
            std::filesystem::remove(temporary.value(), filesystemError);
            if (!filesystemError) temporary.committed();
        } else {
            std::filesystem::rename(temporary.value(), destination, filesystemError);
            if (filesystemError) return status(WVCheckpointStatusCode::commitFailure, "Unable to atomically replace checkpoint: " + filesystemError.message(), destination.string());
            temporary.committed();
        }
        return WVCheckpointStatus::ok();
    } catch (const std::exception& exception) {
        return status(WVCheckpointStatusCode::writeFailure, "Checkpoint writing failed: " + std::string(exception.what()), path);
    }
}

} // namespace wavevortex::runtime
