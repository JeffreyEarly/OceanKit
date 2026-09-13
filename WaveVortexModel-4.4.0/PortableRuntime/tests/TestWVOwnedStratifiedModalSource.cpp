#include "WaveVortexKernel/WVOwnedStratifiedModalSource.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVBoussinesqModalTestFixture.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::test_fixture;

namespace {
WVStratifiedModalArrays arraysFrom(const WVStratifiedModalRecord& record) {
    WVStratifiedModalArrays arrays;
    arrays.geometry=record.geometry();
    arrays.PF0inv=record.PF0inv(); arrays.QG0inv=record.QG0inv();
    arrays.PF0=record.PF0(); arrays.QG0=record.QG0();
    arrays.PFpmInv=record.PFpmInv(); arrays.QGpmInv=record.QGpmInv();
    arrays.PFpm=record.PFpm(); arrays.QGpm=record.QGpm(); arrays.QGwg=record.QGwg();
    return arrays;
}

std::shared_ptr<const WVStratifiedModalRecord> read(const std::filesystem::path& path) {
    std::shared_ptr<const WVStratifiedModalRecord> result;
    const auto status=WVStratifiedModalReader::read(path.string(),result);
    require(static_cast<bool>(status),status.message.c_str());
    return result;
}

void sharedFamilies() {
    Temporary file; fixture(file.path);
    const auto record=read(file.path);
    for (const auto* family:{"WVTransformStratifiedQG","WVTransformHydrostatic"}) {
        auto arrays=arraysFrom(*record); arrays.geometry.transformClass=family;
        std::shared_ptr<const WVOwnedStratifiedModalSource> source;
        auto status=WVOwnedStratifiedModalSource::create(std::move(arrays),source);
        require(static_cast<bool>(status),status.message.c_str());
        require(source->geometry().transformClass==family,"Owned source lost its transform family.");
        require(source->PF0inv()==record->PF0inv() && source->groups().size()==1 &&
            source->groups()[0].columns==std::vector<std::size_t>({0,1,2}),
            "Owned shared modal state changed values or membership.");
        WVComplexLayout input{3,3,1,3,WVComplexRepresentation::interleaved,
            "F-modal",source->modeSetIdentity()};
        WVComplexLayout output{7,3,1,7,WVComplexRepresentation::interleaved,
            "F-grid",source->modeSetIdentity()};
        std::unique_ptr<WVVerticalMatrixBackend> backend;
        require(static_cast<bool>(WVCreateScalarMatrixBackend(backend)),"Matrix backend creation failed.");
        std::unique_ptr<WVPreparedVerticalOperator> prepared;
        status=source->prepareVertical(WVStratifiedModalOperator::reconstructF,
            input,output,std::move(backend),prepared);
        require(static_cast<bool>(status) && prepared->uniqueMatrixCount()==1,
            "Owned shared modal operator preparation failed.");
    }
}

void waveFamily() {
    Temporary file; boussinesqFixture(file.path);
    const auto record=read(file.path);
    auto arrays=arraysFrom(*record);
    std::shared_ptr<const WVOwnedStratifiedModalSource> source;
    auto status=WVOwnedStratifiedModalSource::create(std::move(arrays),source);
    require(static_cast<bool>(status),status.message.c_str());
    require(source->groups().size()==3 && source->groups()[2].columns==std::vector<std::size_t>({2}),
        "Owned Boussinesq source lost exact persisted groups.");
    WVComplexLayout input{3,3,1,3,WVComplexRepresentation::interleaved,
        "Fw-modal",source->modeSetIdentity()};
    WVComplexLayout output{7,3,1,7,WVComplexRepresentation::interleaved,
        "F-grid",source->modeSetIdentity()};
    std::unique_ptr<WVVerticalMatrixBackend> backend;
    require(static_cast<bool>(WVCreateScalarMatrixBackend(backend)),"Matrix backend creation failed.");
    std::unique_ptr<WVPreparedVerticalOperator> prepared;
    status=source->prepareVertical(WVStratifiedModalOperator::reconstructFw,
        input,output,std::move(backend),prepared);
    require(static_cast<bool>(status) && prepared->uniqueMatrixCount()==3 &&
        prepared->preparedGroupCount()==3,"Owned grouped wave preparation failed.");
}

void failurePreservesResult() {
    Temporary file; fixture(file.path);
    const auto record=read(file.path);
    auto valid=arraysFrom(*record);
    std::shared_ptr<const WVOwnedStratifiedModalSource> source;
    require(static_cast<bool>(WVOwnedStratifiedModalSource::create(std::move(valid),source)),
        "Valid owned source creation failed.");
    const auto retained=source;
    auto invalid=arraysFrom(*record); invalid.PF0inv.pop_back();
    const auto status=WVOwnedStratifiedModalSource::create(std::move(invalid),source);
    require(!status && source==retained,"Failed owned source creation replaced the published result.");
}
} // namespace

int main() {
    try {
        sharedFamilies(); waveFamily(); failurePreservesResult();
        std::cout << "owned stratified modal source tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
