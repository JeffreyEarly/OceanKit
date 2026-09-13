#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
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
        if (argc!=5) throw std::runtime_error("Usage: WVHydrostaticKernelDump record.nc input.json output.json reference|native|native-accelerate[-pruned|-compact]");
        std::ifstream inputFile(argv[2]); json data; inputFile>>data;
        std::shared_ptr<const WVStratifiedModalRecord> record; require(WVStratifiedModalReader::read(argv[1],record));
        std::weak_ptr<const WVStratifiedModalRecord> owner=record;
        std::unique_ptr<WVFFTEngine> engine; const std::string requestedProvider=argv[4]; const bool compact=hasSuffix(requestedProvider,"-compact"); const bool pruned=compact || hasSuffix(requestedProvider,"-pruned"); const std::string provider=pruned ? requestedProvider.substr(0,requestedProvider.size()-(compact ? 8 : 7)) : requestedProvider;
        if (provider=="reference") engine=std::make_unique<WVReferenceFFTEngine>();
        else if (provider=="native" || provider=="native-accelerate") {
#if WV_TEST_NATIVE_FFTW
            require(WVFFTWEngine::create(1,engine));
#else
            throw std::runtime_error("Native FFTW adapter was not built.");
#endif
        } else throw std::runtime_error("Unknown provider.");
        std::unique_ptr<WVTransformHydrostaticKernel> kernel;
        WVVariableExecutionOptions options; if (pruned) options={WVRetainedHorizontalSchedule::streamingPrunedTile16,2,true};
        if (compact) { options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews; options.pointwiseWorkers=2; options.fusedDerivativeAdvection=true; options.tiledNonlinear=true; }
        require(WVTransformHydrostaticKernel::create(record,std::move(engine),kernel,provider=="native-accelerate" ? WVCreateAccelerateMatrixBackend : WVCreateScalarMatrixBackend,options)); record.reset();
        const auto& g=kernel->geometry(); const auto S=g.Nj*g.Nkl,R=g.Nx*g.Ny*g.Nz; const auto shape=kernel->spectralShape(); const auto volume=kernel->spatialShape();
        std::array<std::vector<WVComplex64>,3> a,c;
        const char* coefficients[]={"Ap","Am","A0"};
        for (std::size_t k=0;k<3;++k) { a[k]=complexInput(data.at(coefficients[k])); if(a[k].size()!=S) throw std::runtime_error("Wrong coefficient size."); c[k].resize(S); }
        const auto original=a;
        WVState state{data.at("t"),data.at("t0"),{{a[0].data(),shape},{a[1].data(),shape},{a[2].data(),shape}}};
        WVMutableCoefficients out{{c[0].data(),shape},{c[1].data(),shape},{c[2].data(),shape}};
        WVFlux flux{out.Ap,out.Am,out.A0};
        std::vector<double> tiledFields(kernel->supportsTiledNonlinear() ? 4*R : 0);
        const auto nonlinear=[&](const WVState& sample,WVFlux& target) {
            return kernel->supportsTiledNonlinear() ? kernel->nonlinearFluxAndFields(sample,target,
                {tiledFields.data(),{g.Nx,g.Ny,g.Nz,4}}) : kernel->nonlinearFlux(sample,target);
        };
        auto packed=[&]() { return json{{"Ap",complexValues(c[0])},{"Am",complexValues(c[1])},{"A0",complexValues(c[2])}}; };
        json result; result["engine"]=kernel->engineIdentifier(); result["backend"]=kernel->matrixBackendIdentifier(); result["contract"]=WVHydrostaticKernelContract; result["horizontalSchedule"]=kernel->horizontalScheduleIdentifier(); result["streamedNonlinear"]=kernel->executionOptions().streamedNonlinear; result["compactSplitViews"]=kernel->executionOptions().usesCompactSplitViews(); result["pointwiseWorkers"]=kernel->executionOptions().pointwiseWorkers;
        const auto& factors=kernel->factors();
        const std::pair<const char*,double WVHydrostaticModeFactors::*> realFactors[]={
            {"Omega",&WVHydrostaticModeFactors::omega},{"NAp",&WVHydrostaticModeFactors::NAp},{"NA0",&WVHydrostaticModeFactors::NA0},{"PA0",&WVHydrostaticModeFactors::PA0},
            {"ApmN",&WVHydrostaticModeFactors::ApmN},{"A0Z",&WVHydrostaticModeFactors::A0Z},{"A0N",&WVHydrostaticModeFactors::A0N},{"Apm_TE_factor",&WVHydrostaticModeFactors::waveEnergy},
            {"A0_TE_factor",&WVHydrostaticModeFactors::balancedEnergy},{"A0_Psi_factor",&WVHydrostaticModeFactors::psi},{"A0_QGPV_factor",&WVHydrostaticModeFactors::qgpv},{"A0_TZ_factor",&WVHydrostaticModeFactors::enstrophy}};
        for (auto f:realFactors) { std::vector<double> b; for(auto x:factors) b.push_back(x.*f.second); result["factors"][f.first]=b; }
        const std::pair<const char*,WVComplex64 WVHydrostaticModeFactors::*> complexFactors[]={
            {"UAp",&WVHydrostaticModeFactors::UAp},{"VAp",&WVHydrostaticModeFactors::VAp},{"WAp",&WVHydrostaticModeFactors::WAp},{"UA0",&WVHydrostaticModeFactors::UA0},{"VA0",&WVHydrostaticModeFactors::VA0},{"ApmD",&WVHydrostaticModeFactors::ApmD}};
        for (auto f:complexFactors) { std::vector<WVComplex64> b; for(auto x:factors) b.push_back(x.*f.second); result["factors"][f.first]=complexValues(b); }
        struct Field { const char* name; WVHydrostaticField id; bool surface; };
        const Field fields[]={{"u",WVHydrostaticField::u,false},{"v",WVHydrostaticField::v,false},{"w",WVHydrostaticField::w,false},{"eta",WVHydrostaticField::eta,false},{"pi",WVHydrostaticField::pi,false},{"p",WVHydrostaticField::p,false},{"psi",WVHydrostaticField::psi,false},{"qgpv",WVHydrostaticField::qgpv,false},{"rho_e",WVHydrostaticField::rhoE,false},{"rho_total",WVHydrostaticField::rhoTotal,false},{"zeta_x",WVHydrostaticField::zetaX,false},{"zeta_y",WVHydrostaticField::zetaY,false},{"zeta_z",WVHydrostaticField::zetaZ,false},{"ssh",WVHydrostaticField::ssh,true},{"ssu",WVHydrostaticField::ssu,true},{"ssv",WVHydrostaticField::ssv,true}};
        const char* names[]={"value","x","y","z"};
        for (const auto& field:fields) for (int d=0;d<4;++d) {
            if ((field.id==WVHydrostaticField::zetaX || field.id==WVHydrostaticField::zetaY) && d) continue;
            std::vector<double> b(field.surface ? g.Nx*g.Ny : R);
            require(kernel->transformStateField(state,field.id,{b.data(),{g.Nx,g.Ny,field.surface ? 1 : g.Nz}},static_cast<WVHydrostaticDerivative>(d)));
            result["fields"][field.name][names[d]]=b;
        }
        auto U=result["fields"]["u"]["value"].get<std::vector<double>>(),V=result["fields"]["v"]["value"].get<std::vector<double>>(),N=result["fields"]["eta"]["value"].get<std::vector<double>>();
        require(kernel->transformUVEtaToWaveVortex({U.data(),volume},{V.data(),volume},{N.data(),volume},state.t,state.t0,out)); result["roundtrip"]=packed();
        auto u=data.at("u").get<std::vector<double>>(),v=data.at("v").get<std::vector<double>>(),eta=data.at("eta").get<std::vector<double>>(),q=data.at("q").get<std::vector<double>>();
        for (auto* x:{&u,&v,&eta,&q}) if (x->size()!=R) throw std::runtime_error("Wrong spatial input size.");
        require(kernel->transformUVEtaToWaveVortex({u.data(),volume},{v.data(),volume},{eta.data(),volume},state.t,state.t0,out)); result["project"]=packed();
        require(kernel->evolveCoefficients(state,out)); result["evolved"]=packed();
        require(nonlinear(state,flux)); result["flux"]=packed();
        double scalar=0; require(kernel->totalEnstrophy(state.coefficients,scalar)); result["enstrophy"]=scalar;
        const char* components[]={"all","wave","inertial","geostrophic","meanDensityAnomaly"};
        for (int component=0;component<5;++component) {
            const auto id=static_cast<WVHydrostaticComponent>(component);
            require(kernel->totalEnergy(state.coefficients,scalar,id)); result["components"][components[component]]["energy"]=scalar;
            require(kernel->totalEnergySpatiallyIntegrated(state,scalar,id)); result["components"][components[component]]["spatialEnergy"]=scalar;
            for (const auto& field:fields) if (!field.surface && field.id!=WVHydrostaticField::rhoTotal) {
                require(kernel->transformStateField(state,field.id,{U.data(),volume},WVHydrostaticDerivative::value,id)); result["components"][components[component]][field.name]=U;
            }
        }
        auto raw=a[0];
        for (std::size_t i=0;i<S;++i) if (factors[i].inertial) raw[i].imag=0;
        for (int family=0;family<2;++family) {
            const auto id=static_cast<WVHydrostaticFamily>(family); const char* name=family ? "G" : "F";
            require(kernel->transformToSpatial({raw.data(),shape},id,{U.data(),volume})); result["raw"][name]["inverse"]=U;
            require(kernel->transformFromSpatial({q.data(),volume},id,out.Ap)); result["raw"][name]["forward"]=complexValues(c[0]);
            for (unsigned order=1;order<=4;++order) { require(kernel->differentiateVertical({q.data(),volume},id,order,{U.data(),volume})); result["raw"][name]["d"+std::to_string(order)]=U; }
            require(kernel->integrateVertical({q.data(),volume},id,{U.data(),volume})); result["raw"][name]["integral"]=U;
        }
        for (auto& x:c) std::fill(x.begin(),x.end(),WVComplex64{.1,.2});
        require(kernel->constrainCoefficients(out)); result["constrained"]=packed();
        // A short fixed-step trajectory exercises changing phases and nonlinear
        // amplitudes without depending on the runtime integrator (#300).
        std::array<std::vector<WVComplex64>,3> trajectory=a, stage=a;
        std::array<std::array<std::vector<WVComplex64>,3>,4> slopes;
        for (auto& slope:slopes) for (auto& x:slope) x.resize(S);
        constexpr double dt=.5;
        for (int step=0;step<4;++step) {
            for (int k=0;k<4;++k) {
                const double fraction=k==0 ? 0 : k==3 ? 1 : .5;
                for (int j=0;j<3;++j) for (std::size_t i=0;i<S;++i) {
                    stage[j][i]=trajectory[j][i];
                    if (k) { stage[j][i].real+=dt*fraction*slopes[k-1][j][i].real; stage[j][i].imag+=dt*fraction*slopes[k-1][j][i].imag; }
                }
                WVState sample{state.t+dt*(step+fraction),state.t0,{{stage[0].data(),shape},{stage[1].data(),shape},{stage[2].data(),shape}}};
                WVFlux target{{slopes[k][0].data(),shape},{slopes[k][1].data(),shape},{slopes[k][2].data(),shape}};
                require(nonlinear(sample,target));
            }
            for (int j=0;j<3;++j) for (std::size_t i=0;i<S;++i) {
                trajectory[j][i].real+=dt/6*(slopes[0][j][i].real+2*slopes[1][j][i].real+2*slopes[2][j][i].real+slopes[3][j][i].real);
                trajectory[j][i].imag+=dt/6*(slopes[0][j][i].imag+2*slopes[1][j][i].imag+2*slopes[2][j][i].imag+slopes[3][j][i].imag);
            }
        }
        for (int j=0;j<3;++j) result["trajectory"][coefficients[j]]=complexValues(trajectory[j]);
        const auto bytes=kernel->persistentBytes();
        allocationProbe::calls=0; allocationProbe::counting=true;
        for (int i=0;i<3;++i) {
            require(nonlinear(state,flux)); require(kernel->transformStateField(state,WVHydrostaticField::rhoE,{U.data(),volume},WVHydrostaticDerivative::z));
            require(kernel->transformUVEtaToWaveVortex({u.data(),volume},{v.data(),volume},{eta.data(),volume},state.t,state.t0,out));
            require(kernel->differentiateVertical({q.data(),volume},WVHydrostaticFamily::F,4,{U.data(),volume}));
            require(kernel->integrateVertical({q.data(),volume},WVHydrostaticFamily::G,{U.data(),volume}));
        }
        allocationProbe::counting=false; result["tiledNonlinearExecutions"]=kernel->metrics().tiledNonlinearCount; result["preparedAllocations"]=allocationProbe::calls.load();
        bool unchanged=true; for (int j=0;j<3;++j) for (std::size_t i=0;i<S;++i) unchanged=unchanged && a[j][i].real==original[j][i].real && a[j][i].imag==original[j][i].imag;
        result["inputPreserved"]=unchanged; result["storageStable"]=bytes==kernel->persistentBytes(); result["realScratchBytes"]=kernel->storage().realScratchBytes;
        kernel.reset(); result["ownerReleased"]=owner.expired();
        std::ofstream outputFile(argv[3]); outputFile<<result.dump()<<'\n'; return 0;
    } catch (const std::exception& e) { allocationProbe::counting=false; allocationProbe::failAfter=-1; std::cerr<<e.what()<<'\n'; return 1; }
}
