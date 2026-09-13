#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVAllocationProbe.hpp"
#include "WVStratifiedModalTestFixture.hpp"
#include <netcdf.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;
namespace {
using namespace wavevortex::test_fixture;
void contracts(const std::filesystem::path& path) {
    WVStratifiedModalInspection inspection; auto status=WVStratifiedModalReader::inspect(path.string(),inspection); require(static_cast<bool>(status),status.message.c_str());
    std::shared_ptr<const WVStratifiedModalRecord> record;status=WVStratifiedModalReader::read(path.string(),record);require(static_cast<bool>(status),status.message.c_str());
    require(record->PF0inv().size()*4*sizeof(double)==inspection.scientificMatrixBytes,"Scientific matrix accounting mismatch");
    WVVerticalMatrixRecord raw; require(static_cast<bool>(record->matrixView(WVStratifiedScientificMatrix::PF0inv,raw)),"Scientific view failed");
    require(raw.values.data==record->PF0inv().data()&&raw.values.rows==7&&raw.values.columns==3&&raw.identity.source==record->sourceIdentity()&&raw.identity.name=="PF0inv","Scientific view lost exact identity or orientation");
    require(record->groups().size()==1&&record->groups()[0].identity==0&&record->groups()[0].columns==std::vector<std::size_t>{0,1,2},"Scientific group membership changed");
    const auto saved=record->PF0inv();
    WVComplexLayout input{3,3,1,3,WVComplexRepresentation::interleaved,"F-modal",record->modeSetIdentity()},output{7,3,1,7,WVComplexRepresentation::interleaved,"F-grid",record->modeSetIdentity()};
    std::unique_ptr<WVVerticalMatrixBackend> backend; require(static_cast<bool>(WVCreateScalarMatrixBackend(backend)),"Backend failed");
    std::unique_ptr<WVPreparedVerticalOperator> prepared;require(static_cast<bool>(record->prepareVertical(WVStratifiedModalOperator::reconstructF,input,output,std::move(backend),prepared)),"Preparation failed");
    std::unique_ptr<WVVerticalWorkspace> workspace;require(static_cast<bool>(prepared->createWorkspace(workspace)),"Workspace failed");
    std::vector<WVComplex64> in(9,{1,2}),out(21);
    allocationProbe::calls=0;allocationProbe::counting=true;
    for(int i=0;i<10;++i)require(static_cast<bool>(prepared->execute(*workspace,{in.data(),nullptr,nullptr,in.size()*sizeof(WVComplex64)},{out.data(),nullptr,nullptr,out.size()*sizeof(WVComplex64)})),"Execution failed");
    allocationProbe::counting=false;require(allocationProbe::calls==0,"Prepared modal execution allocated");
    require(record->PF0inv()==saved,"Preparation altered scientific matrices");
    change(path,"P0",7,1);
    require(record->geometry().P0[1]==2,"Scientific record borrowed mutable file data");
    std::shared_ptr<const WVStratifiedModalRecord> fresh;require(static_cast<bool>(WVStratifiedModalReader::read(path.string(),fresh)),"Rebuild read failed");
    require(fresh->sourceIdentity()!=record->sourceIdentity()&&fresh->geometry().P0[1]==7,"Changed scientific input reused identity");
    require(static_cast<bool>(WVCreateScalarMatrixBackend(backend)),"Backend failed");
    require(fresh->prepareVertical(WVStratifiedModalOperator::reconstructF,input,output,std::move(backend),prepared).code==WVKernelStatusCode::invalidConfiguration,"Old field identities accepted");
    record.reset();require(static_cast<bool>(prepared->execute(*workspace,{in.data(),nullptr,nullptr,in.size()*sizeof(WVComplex64)},{out.data(),nullptr,nullptr,out.size()*sizeof(WVComplex64)})),"Prepared operator depended on released scientific input");
}
void rejection(const std::filesystem::path& path) {
    struct Bad {const char* name;double value;std::size_t first,second;};
    for(const auto& bad:{Bad{"P0",0,1,0},{"Q0",0,1,0},{"h_0",0,1,0},{"N2",-1,3,0},{"rho_nm0",2000,6,0},{"z",-1100,2,0},
        {"j",1,2,0},{"j",6,2,0},{"k",1.123,1,0},{"k",0,1,0},{"l",-0.0001,2,0},{"QG0inv",1,1,0},{"QG0inv",1,0,3},{"QG0",1,0,1},
        {"PF0inv",2,0,3},{"PF0inv",2,1,2},{"PF0",NC_FILL_DOUBLE,5,1},{"QG0",std::numeric_limits<double>::quiet_NaN(),5,1},{"z_int",-1000,2,0}}) {
        fixture(path);std::shared_ptr<const WVStratifiedModalRecord> record;require(static_cast<bool>(WVStratifiedModalReader::read(path.string(),record)),"Base fixture invalid");
        const auto* saved=record.get();WVStratifiedModalInspection inspection;inspection.scientificMatrixBytes=17;
        change(path,bad.name,bad.value,bad.first,bad.second);
        require(!WVStratifiedModalReader::inspect(path.string(),inspection)&&inspection.scientificMatrixBytes==17,"Invalid record accepted or inspection output mutated");
        require(!WVStratifiedModalReader::read(path.string(),record)&&record.get()==saved,"Invalid record replaced scientific owner");
    }
    fixture(path);{File f(path);int v;nc(nc_inq_varid(f.id,"PF0",&v));nc(nc_rename_var(f.id,v,"missing_PF0"));}
    WVStratifiedModalInspection inspection;require(WVStratifiedModalReader::inspect(path.string(),inspection).code==WVCheckpointStatusCode::missingVariable,"Missing matrix accepted");
    fixture(path);{File f(path);int v;nc(nc_inq_varid(f.id,"PF0",&v));nc(nc_rename_var(f.id,v,"old_PF0"));int j,z;nc(nc_inq_dimid(f.id,"j",&j));nc(nc_inq_dimid(f.id,"z",&z));int dims[]={j,z};nc(nc_def_var(f.id,"PF0",NC_DOUBLE,2,dims,&v));}
    require(WVStratifiedModalReader::inspect(path.string(),inspection).code==WVCheckpointStatusCode::shapeMismatch,"Transposed matrix schema accepted");
    fixture(path,4101);change(path,"PF0",std::numeric_limits<double>::infinity(),4100,2);
    require(!WVStratifiedModalReader::inspect(path.string(),inspection),"Late payload outside scan buffer was not validated");
}
void allocationFailures(const std::filesystem::path& path) {
    fixture(path);const auto name=path.string();bool success=false;
    for(long i=0;i<256;++i){std::shared_ptr<const WVStratifiedModalRecord> record;allocationProbe::failAfter=i;const auto status=WVStratifiedModalReader::read(name,record);allocationProbe::failAfter=-1;if(status){success=true;break;}require(!record&&status.code==WVCheckpointStatusCode::allocationFailure,"Allocation failure escaped or published partial record");}
    require(success,"Allocation sweep did not reach successful read");
}
}
int main(){try{Temporary temporary;fixture(temporary.path);contracts(temporary.path);rejection(temporary.path);allocationFailures(temporary.path);std::cout<<"Modal records: validated scientific ownership, bounded payload scan, invalid records, rebuilds and prepared allocations passed.\n";return 0;}catch(const std::exception& e){allocationProbe::counting=false;allocationProbe::failAfter=-1;std::cerr<<e.what()<<'\n';return 1;}}
