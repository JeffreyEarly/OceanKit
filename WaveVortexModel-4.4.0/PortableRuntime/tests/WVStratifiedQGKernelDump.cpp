#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"
#include "nlohmann/json.hpp"
#if WV_TEST_NATIVE_FFTW
#include "WVNativeFFTWEngine.hpp"
#endif
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;
namespace {
bool hasSuffix(const std::string& value,const std::string& suffix) { return value.size()>=suffix.size() && value.compare(value.size()-suffix.size(),suffix.size(),suffix)==0; }
template<class T> void require(T status) { if (!status) throw std::runtime_error(status.message); }
json complexValues(const std::vector<WVComplex64>& a) { std::vector<double> r,i; for (auto x:a) { r.push_back(x.real); i.push_back(x.imag); } return {{"real",r},{"imag",i}}; }
std::vector<WVComplex64> complexInput(const json& a) { auto r=a.at("real").get<std::vector<double>>(),i=a.at("imag").get<std::vector<double>>(); if (r.size()!=i.size()) throw std::runtime_error("Input lengths differ."); std::vector<WVComplex64> b(r.size()); for (std::size_t n=0;n<r.size();++n) b[n]={r[n],i[n]}; return b; }
}
int main(int argc,char** argv) {
    try {
        if (argc!=5) throw std::runtime_error("Usage: WVStratifiedQGKernelDump record.nc input.json output.json reference|native|native-accelerate[-pruned|-compact]");
        std::ifstream inputFile(argv[2]); json data; inputFile>>data;
        std::shared_ptr<const WVStratifiedModalRecord> record; require(WVStratifiedModalReader::read(argv[1],record));
        std::unique_ptr<WVFFTEngine> engine;
        const std::string requestedProvider=argv[4];
        const bool compact=hasSuffix(requestedProvider,"-compact");
        const bool pruned=compact || hasSuffix(requestedProvider,"-pruned");
        const auto suffixLength=compact ? 8u : (pruned ? 7u : 0u);
        const std::string provider=suffixLength ? requestedProvider.substr(0,requestedProvider.size()-suffixLength) : requestedProvider;
        if (provider=="reference") engine=std::make_unique<WVReferenceFFTEngine>();
        else {
#if WV_TEST_NATIVE_FFTW
            require(WVFFTWEngine::create(1,engine));
#else
            throw std::runtime_error("Native FFTW adapter was not built.");
#endif
        }
        WVVariableExecutionOptions options;
        if (pruned) options={WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true,
            compact ? WVVariableSpectralSchedule::compactSplitFusedViews :
                WVVariableSpectralSchedule::establishedInterleaved};
        if (compact) options.pointwiseWorkers=2;
        std::unique_ptr<WVTransformStratifiedQGKernel> kernel; require(WVTransformStratifiedQGKernel::create(record,std::move(engine),kernel,provider=="native-accelerate" ? WVCreateAccelerateMatrixBackend : WVCreateScalarMatrixBackend,options));
        const auto& g=kernel->geometry(); const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz;
        auto a=complexInput(data.at("A0")); if (a.size()!=S) throw std::runtime_error("Wrong coefficient size."); const auto original=a;
        WVComplexConstView A0{a.data(),kernel->spectralShape()}; std::vector<WVComplex64> c(S); WVComplexView out{c.data(),kernel->spectralShape()};
        json result; result["engine"]=kernel->engineIdentifier(); result["matrixBackend"]=kernel->matrixBackendIdentifier(); result["contract"]=WVStratifiedQGKernelContract; result["horizontalSchedule"]=kernel->horizontalScheduleIdentifier(); result["streamedNonlinear"]=kernel->executionOptions().streamedNonlinear;
        result["compactSplitViews"]=kernel->executionOptions().usesCompactSplitViews(); result["pointwiseWorkers"]=kernel->executionOptions().pointwiseWorkers;
        result["spectralSchedule"]=kernel->executionOptions().usesCompactSplitViews() ? "compact-split-fused-views" : "established-interleaved";
        const auto& f=kernel->factors(); result["factors"]={{"u",complexValues(f.u)},{"v",complexValues(f.v)},{"eta",f.eta},{"pi",f.pi},{"psi",f.psi},{"qgpv",f.qgpv},{"zetaZ",f.zetaZ},{"energy",f.energy},{"enstrophy",f.enstrophy},{"ke",f.kineticEnergy},{"pe",f.potentialEnergy}};
        struct Field { const char* name; WVStratifiedQGField id; bool surface; };
        const Field fields[]={{"u",WVStratifiedQGField::u,false},{"v",WVStratifiedQGField::v,false},{"w",WVStratifiedQGField::w,false},{"eta",WVStratifiedQGField::eta,false},{"pi",WVStratifiedQGField::pi,false},{"p",WVStratifiedQGField::p,false},{"psi",WVStratifiedQGField::psi,false},{"qgpv",WVStratifiedQGField::qgpv,false},{"rho_e",WVStratifiedQGField::rhoE,false},{"rho_total",WVStratifiedQGField::rhoTotal,false},{"zeta_z",WVStratifiedQGField::zetaZ,false},{"ssh",WVStratifiedQGField::ssh,true},{"ssu",WVStratifiedQGField::ssu,true},{"ssv",WVStratifiedQGField::ssv,true}};
        const char* names[]={"value","x","y","z"};
        for (const auto& field:fields) for (int d=0;d<4;++d) {
            std::vector<double> b(field.surface ? g.Nx*g.Ny : R);
            require(kernel->transformA0ToField(A0,field.id,{b.data(),{g.Nx,g.Ny,field.surface ? 1 : g.Nz}},static_cast<WVStratifiedQGDerivative>(d)));
            result["fields"][field.name][names[d]]=b;
        }
        auto reconstructedQ=result["fields"]["qgpv"]["value"].get<std::vector<double>>();
        auto reconstructedU=result["fields"]["u"]["value"].get<std::vector<double>>();
        auto reconstructedV=result["fields"]["v"]["value"].get<std::vector<double>>();
        auto reconstructedEta=result["fields"]["eta"]["value"].get<std::vector<double>>();
        require(kernel->transformQGPVToA0({reconstructedQ.data(),kernel->spatialShape()},out)); result["roundtripQ"]=complexValues(c);
        require(kernel->transformUVEtaToA0({reconstructedU.data(),kernel->spatialShape()},{reconstructedV.data(),kernel->spatialShape()},{reconstructedEta.data(),kernel->spatialShape()},out)); result["roundtripUVEta"]=complexValues(c);
        auto q=data.at("q").get<std::vector<double>>(),u=data.at("u").get<std::vector<double>>(),v=data.at("v").get<std::vector<double>>(),eta=data.at("eta").get<std::vector<double>>();
        for (const auto* x:{&q,&u,&v,&eta}) if(x->size()!=R) throw std::runtime_error("Wrong spatial input size.");
        require(kernel->transformQGPVToA0({q.data(),kernel->spatialShape()},out)); result["projectQ"]=complexValues(c);
        require(kernel->transformUVEtaToA0({u.data(),kernel->spatialShape()},{v.data(),kernel->spatialShape()},{eta.data(),kernel->spatialShape()},out)); result["projectUVEta"]=complexValues(c);
        require(kernel->nonlinearFlux(A0,out)); result["flux"]=complexValues(c);
        require(kernel->nonlinearFlux(A0,out,f.beta)); result["fluxBeta"]=complexValues(c);
        require(kernel->verticalDiffusivityFlux(A0,.002,out)); result["verticalDiffusivity"]=complexValues(c);
        require(kernel->linearBottomFrictionFlux(A0,1e-5,out)); result["linearBottomFriction"]=complexValues(c);
        require(kernel->quadraticBottomFrictionFlux(A0,.003,out)); result["quadraticBottomFriction"]=complexValues(c);
        require(kernel->linearFlux(A0,out,f.beta)); result["linearBeta"]=complexValues(c);
        require(kernel->evolveA0(A0,1234,out)); result["stationary"]=complexValues(c);
        require(kernel->evolveA0(A0,1234,out,f.beta)); result["evolvedBeta"]=complexValues(c);
        double scalar=0; require(kernel->totalEnergy(A0,scalar)); result["energy"]=scalar;
        require(kernel->totalEnstrophy(A0,scalar)); result["enstrophy"]=scalar;
        require(kernel->totalEnergySpatiallyIntegrated(A0,scalar)); result["spatialEnergy"]=scalar;
        require(kernel->totalEnstrophySpatiallyIntegrated(A0,scalar)); result["spatialEnstrophy"]=scalar;
        require(kernel->uvMax(A0,scalar)); result["uvMax"]=scalar;
        allocationProbe::calls=0; allocationProbe::counting=true;
        for (int i=0;i<3;++i) { require(kernel->nonlinearFlux(A0,out)); require(kernel->verticalDiffusivityFlux(A0,.002,out)); require(kernel->linearBottomFrictionFlux(A0,1e-5,out)); require(kernel->quadraticBottomFrictionFlux(A0,.003,out)); require(kernel->transformA0ToField(A0,WVStratifiedQGField::eta,{q.data(),kernel->spatialShape()})); }
        allocationProbe::counting=false; result["preparedAllocations"]=allocationProbe::calls.load();
        bool unchanged=true; for (std::size_t i=0;i<S;++i) unchanged=unchanged && a[i].real==original[i].real && a[i].imag==original[i].imag;
        result["inputPreserved"]=unchanged; result["scientificBytes"]=kernel->storage().sharedScientificBytes;
        result["spectralScratchBytes"]=kernel->storage().spectralScratchBytes; result["realScratchBytes"]=kernel->storage().realScratchBytes;
        std::ofstream outputFile(argv[3]); outputFile<<result.dump()<<'\n'; return 0;
    } catch (const std::exception& e) { allocationProbe::counting=false; allocationProbe::failAfter=-1; std::cerr<<e.what()<<'\n'; return 1; }
}
