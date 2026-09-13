#include "WaveVortexKernel/WVSpectralOperators.hpp"
#include "WVNativeFFTWEngine.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace wavevortex;
static void require(WVKernelStatus s) { if (!s) throw std::runtime_error(s.message); }
int main(int argc,char** argv) try {
    if (argc!=4) throw std::runtime_error("Usage: NativeStage mode-file targets output-prefix");
    std::ifstream input(argv[1]); std::size_t nx,ny,nz,m;
    double lx,ly; input>>nx>>ny>>nz>>m>>lx>>ly;
    if (!input || !m) throw std::runtime_error("Invalid mode file");
    const auto targets=static_cast<std::size_t>(std::stoul(argv[2]));
    WVRetainedHorizontalSpecification spec;
    spec.grid={nx,ny,nz,1,nx,nx*ny,"F"}; spec.Lx=lx; spec.Ly=ly;
    spec.schedule=WVRetainedHorizontalSchedule::streamingPrunedTile16;
    spec.outerWorkers=12; spec.normalization=WVFourierNormalization::forwardUnit;
    spec.retained.rows=nz; spec.retained.columns=m; spec.retained.rowStride=1;
    spec.retained.columnStride=nz; spec.retained.representation=WVComplexRepresentation::split;
    spec.retained.family="F"; spec.retained.modeSet="deterministic-stage";
    for (std::size_t i=0;i<m;++i) { std::int64_t k,l;input>>k>>l;spec.modes.push_back({k,l}); }
    if (!input) throw std::runtime_error("Truncated mode file");
    std::unique_ptr<WVFFTEngine> engine; require(WVFFTWEngine::create(1,engine));
    std::unique_ptr<WVRetainedHorizontalOperator> op;
    require(WVRetainedHorizontalOperator::create(spec,std::move(engine),op));
    std::unique_ptr<WVRetainedHorizontalWorkspace> workspace;
    require(op->createWorkspace(workspace,false));require(op->prepareAdvection(*workspace,targets));
    const auto n=nz*m,volume=nx*ny*nz;
    std::vector<double> br(4*n),bi(4*n),zr(targets*n),zi(targets*n),tr(targets*n),ti(targets*n);
    std::vector<double> fields(4*volume),correction(nz);
    for (std::size_t f=0;f<4;++f) for (std::size_t mode=0;mode<m;++mode)
        for (std::size_t z=0;z<nz;++z) {
            const auto i=f*n+mode*nz+z;
            br[i]=std::sin(.017*(1+f+mode+2*z));
            bi[i]=spec.modes[mode].k==0 && spec.modes[mode].l==0 ? 0.0 : std::cos(.031*(1+2*f+3*mode+z));
        }
    for (std::size_t t=0;t<targets;++t) for (std::size_t mode=0;mode<m;++mode)
        for (std::size_t z=0;z<nz;++z) {
            const auto i=t*n+mode*nz+z;
            zr[i]=.002*std::cos(.019*(1+t+mode+z));
            zi[i]=spec.modes[mode].k==0 && spec.modes[mode].l==0 ? 0.0 : .002*std::sin(.023*(1+2*t+mode+z));
        }
    for (std::size_t z=0;z<nz;++z) correction[z]=.001*std::cos(.04*z);
    WVRetainedAdvectionWork work; work.targets=targets;
    work.fields={fields.data(),fields.size()*sizeof(double)};
    work.densityCorrection={correction.data(),correction.size()*sizeof(double)};
    for (std::size_t f=0;f<4;++f) work.base[f]={nullptr,br.data()+f*n,bi.data()+f*n,n*sizeof(double)};
    for (std::size_t t=0;t<targets;++t) work.targetSpectra[t]={nullptr,tr.data()+t*n,ti.data()+t*n,n*sizeof(double)};
    WVRetainedAdvectionCounts counts;std::vector<double> seconds;
    for (int repeat=-2;repeat<9;++repeat) {
        std::copy(zr.begin(),zr.end(),tr.begin());std::copy(zi.begin(),zi.end(),ti.begin());
        const auto start=std::chrono::steady_clock::now();require(op->advection(*workspace,work,counts));
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if (repeat>=0) seconds.push_back(elapsed);
        if (counts.columnInverses!=nz*(4+2*targets) || counts.rowInverses!=nz*(4+3*targets) || counts.reusedColumns!=nz*targets)
            throw std::runtime_error("Producer count mismatch");
    }
    std::ofstream output(std::string(argv[3])+".bin",std::ios::binary);
    for (const auto* v:{&fields,&tr,&ti}) output.write(reinterpret_cast<const char*>(v->data()),static_cast<std::streamsize>(v->size()*sizeof(double)));
    if (!output) throw std::runtime_error("Cannot write stage output");
    std::cout<<std::setprecision(17)<<"{\"workspaceBytes\":"<<workspace->persistentBytes()<<",\"seconds\":[";
    for (std::size_t i=0;i<seconds.size();++i) std::cout<<(i?",":"")<<seconds[i];
    std::cout<<"],\"columnInverses\":"<<counts.columnInverses<<",\"rowInverses\":"<<counts.rowInverses<<",\"reusedColumns\":"<<counts.reusedColumns<<"}\n";
} catch (const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
