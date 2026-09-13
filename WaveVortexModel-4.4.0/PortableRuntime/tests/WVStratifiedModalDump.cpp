#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {
template<class Status> void require(const Status& status) { if (!status) throw std::runtime_error(status.message); }
json complexValues(const std::vector<WVComplex64>& a) {
    std::vector<double> real,imag; for (auto v:a) { real.push_back(v.real); imag.push_back(v.imag); }
    return {{"real",real},{"imag",imag}};
}
}
int main(int argc,char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("Usage: WVStratifiedModalDump checkpoint.nc");
        WVStratifiedModalInspection inspection; require(WVStratifiedModalReader::inspect(argv[1],inspection));
        std::shared_ptr<const WVStratifiedModalRecord> record; require(WVStratifiedModalReader::read(argv[1],record));
        const auto& g=record->geometry();
        json report={{"contract",WVStratifiedModalRecordContract},{"matrixBytes",inspection.scientificMatrixBytes},{"scanBytes",inspection.matrixScanBytes},
            {"PF0inv",record->PF0inv()},{"QG0inv",record->QG0inv()},{"PF0",record->PF0()},{"QG0",record->QG0()},
            {"x",g.x},{"y",g.y},{"z",g.z},{"j",g.j},{"k",g.k},{"l",g.l},{"N2",g.N2},{"rho_nm0",g.rho_nm0},{"dLnN2",g.dLnN2},{"P0",g.P0},{"Q0",g.Q0},{"h_0",g.h_0},{"z_int",g.z_int}};
        struct Operation { WVStratifiedModalOperator kind; const char* name; const char* input; const char* output; bool reconstruction,projection; };
        for (const auto& op : {Operation{WVStratifiedModalOperator::reconstructF,"Finv","F-modal","F-grid",true,false},
            {WVStratifiedModalOperator::projectF,"F","F-grid","F-modal",false,true},
            {WVStratifiedModalOperator::reconstructG,"Ginv","G-modal","G-grid",true,false},
            {WVStratifiedModalOperator::projectG,"G","G-grid","G-modal",false,true},
            {WVStratifiedModalOperator::GToF,"GToF","G-modal","F-modal",false,false},
            {WVStratifiedModalOperator::FToG,"FToG","F-modal","G-modal",false,false}}) {
            const auto m=op.reconstruction ? g.Nz : g.Nj, k=op.projection ? g.Nz : g.Nj;
            for (auto representation : {WVComplexRepresentation::interleaved,WVComplexRepresentation::split}) {
                WVComplexLayout input{k,g.Nkl,1,k,representation,op.input,record->modeSetIdentity()}, output{m,g.Nkl,1,m,representation,op.output,record->modeSetIdentity()};
                std::unique_ptr<WVVerticalMatrixBackend> backend; require(WVCreateScalarMatrixBackend(backend));
                std::unique_ptr<WVPreparedVerticalOperator> prepared; require(record->prepareVertical(op.kind,input,output,std::move(backend),prepared));
                if (prepared->uniqueMatrixCount()!=1) throw std::runtime_error("Shared matrix expanded per mode.");
                std::unique_ptr<WVVerticalWorkspace> workspace; require(prepared->createWorkspace(workspace));
                std::vector<WVComplex64> b(k*g.Nkl),c(m*g.Nkl); std::vector<double> br(b.size()),bi(b.size()),cr(c.size()),ci(c.size());
                for (std::size_t col=0;col<g.Nkl;++col) for (std::size_t row=0;row<k;++row) {
                    const auto i=row+k*col; b[i]={0.1*(row+1)+0.03*col,-0.2*(row+1)+0.01*col}; br[i]=b[i].real; bi[i]=b[i].imag;
                }
                if (representation==WVComplexRepresentation::interleaved) require(prepared->execute(*workspace,{b.data(),nullptr,nullptr,b.size()*sizeof(WVComplex64)},{c.data(),nullptr,nullptr,c.size()*sizeof(WVComplex64)}));
                else {
                    require(prepared->execute(*workspace,{nullptr,br.data(),bi.data(),br.size()*sizeof(double)},{nullptr,cr.data(),ci.data(),cr.size()*sizeof(double)}));
                    for (std::size_t i=0;i<c.size();++i) c[i]={cr[i],ci[i]};
                }
                report[op.name][representation==WVComplexRepresentation::split ? "split" : "interleaved"]=complexValues(c);
            }
        }
        WVRetainedHorizontalSpecification h; require(record->horizontalSpecification(g.Nj,WVComplexRepresentation::interleaved,"F-modal",h));
        std::unique_ptr<WVRetainedHorizontalOperator> horizontal; require(WVRetainedHorizontalOperator::create(h,std::make_unique<WVReferenceFFTEngine>(),horizontal));
        std::unique_ptr<WVRetainedHorizontalWorkspace> workspace; require(horizontal->createWorkspace(workspace));
        std::vector<double> grid(g.Nx*g.Ny*g.Nj); std::vector<WVComplex64> modes(g.Nj*g.Nkl);
        for (std::size_t i=0;i<grid.size();++i) grid[i]=std::sin(0.07*i)+0.2*std::cos(0.19*i);
        require(horizontal->forward(*workspace,{grid.data(),grid.size()*sizeof(double)},{modes.data(),nullptr,nullptr,modes.size()*sizeof(WVComplex64)}));
        report["horizontal"]=complexValues(modes);
        std::cout << report.dump() << '\n'; return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
