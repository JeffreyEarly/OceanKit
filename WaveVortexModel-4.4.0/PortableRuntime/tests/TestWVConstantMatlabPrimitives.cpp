#include "WaveVortexKernel/WVTransformConstantStratificationKernel.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace wavevortex;
namespace {
constexpr double pi=3.141592653589793238462643383279502884;
void require(bool value,const std::string& message) { if (!value) { std::cerr<<"FAIL: "<<message<<'\n'; std::exit(1); } }
std::complex<double> value(WVComplex64 a) { return {a.real,a.imag}; }
WVComplex64 value(std::complex<double> a) { return {a.real(),a.imag()}; }
void near(WVComplex64 a,std::complex<double> b,double tolerance,const std::string& message) {
    if (std::abs(value(a)-b)>=tolerance*std::max(1.0,std::abs(b))) {
        std::cerr<<"FAIL: "<<message<<" actual="<<a.real<<","<<a.imag
            <<" expected="<<b.real()<<","<<b.imag()<<'\n'; std::exit(1);
    }
}
WVTransformConstantStratificationConfiguration configuration(bool hydrostatic) {
    WVTransformConstantStratificationConfiguration c;
    c.Nx=8;c.Ny=6;c.Nz=7;c.Nj=4;c.Lx=2*pi;c.Ly=2*pi;c.Lz=1300;
    c.N0=5.2e-3;c.rho0=1027;c.g=9.80665;c.planetaryRadius=6.3712e6;
    c.rotationRate=7.292115e-5;c.latitude=33;c.isHydrostatic=hydrostatic;c.shouldAntialias=true;
    return c;
}
std::unique_ptr<WVTransformConstantStratificationKernel> kernel(bool hydrostatic,
    WVConstantNonlinearFluxSchedule schedule) {
    auto c=configuration(hydrostatic);
    std::unique_ptr<WVTransformConstantStratificationKernel> result;
    WVConstantKernelExecutionOptions options; options.schedule=schedule; options.horizontalOuterWorkers=2; options.pointwiseWorkers=2;
    require(bool(WVTransformConstantStratificationKernel::create(c,
        std::make_unique<WVReferenceFFTEngine>(),result,options)),"constant kernel creation");
    return result;
}
std::unique_ptr<WVTransformConstantStratificationKernel> kernel(
    WVTransformConstantStratificationConfiguration c,
    WVConstantNonlinearFluxSchedule schedule) {
    std::unique_ptr<WVTransformConstantStratificationKernel> result;
    WVConstantKernelExecutionOptions options; options.schedule=schedule; options.horizontalOuterWorkers=2; options.pointwiseWorkers=2;
    require(bool(WVTransformConstantStratificationKernel::create(c,
        std::make_unique<WVReferenceFFTEngine>(),result,options)),"constant kernel creation");
    return result;
}
std::complex<double> dct(const std::vector<WVComplex64>& a,std::size_t column,
    std::size_t Nz,std::size_t j) {
    std::complex<double> sum=.5*value(a[Nz*column])+
        .5*value(a[Nz*column+Nz-1])*std::cos(pi*static_cast<double>(j));
    for (std::size_t z=1;z+1<Nz;++z)
        sum+=value(a[Nz*column+z])*std::cos(pi*static_cast<double>(j*z)/static_cast<double>(Nz-1));
    return 2.0*sum/static_cast<double>(Nz-1);
}
std::complex<double> dst(const std::vector<WVComplex64>& a,std::size_t column,
    std::size_t Nz,std::size_t j) {
    if (!j) return {};
    std::complex<double> sum;
    for (std::size_t z=1;z+1<Nz;++z)
        sum+=value(a[Nz*column+z])*std::sin(pi*static_cast<double>(j*z)/static_cast<double>(Nz-1));
    return 2.0*sum/static_cast<double>(Nz-1);
}
std::complex<double> inverseDct(const std::vector<WVComplex64>& a,std::size_t column,
    std::size_t Nj,std::size_t z) {
    std::complex<double> sum=.5*value(a[Nj*column]);
    for (std::size_t j=1;j<Nj;++j)
        sum+=value(a[Nj*column+j])*std::cos(pi*static_cast<double>(j*z)/6.0);
    return sum;
}
std::complex<double> inverseDst(const std::vector<WVComplex64>& a,std::size_t column,
    std::size_t Nj,std::size_t z) {
    std::complex<double> sum;
    for (std::size_t j=1;j<Nj;++j)
        sum+=value(a[Nj*column+j])*std::sin(pi*static_cast<double>(j*z)/6.0);
    return sum;
}
std::vector<double> verticalCalculusOracle(const std::vector<double>& input,
    const WVTransformConstantStratificationConfiguration& c,std::size_t columns,
    bool inputIsF,unsigned order,bool integral) {
    std::vector<double> result(c.Nz*columns),coefficients(c.Nz);
    const bool outputIsF=integral ? !inputIsF : ((order%2)==0 ? inputIsF : !inputIsF);
    const double signsF[]={-1,-1,1,1},signsG[]={1,-1,-1,1};
    for(std::size_t column=0;column<columns;++column) {
        for(std::size_t j=0;j<c.Nz;++j) {
            double sum=0;
            if(inputIsF) {
                sum=.5*input[c.Nz*column]+.5*input[c.Nz*column+c.Nz-1]*std::cos(pi*static_cast<double>(j));
                for(std::size_t z=1;z+1<c.Nz;++z)
                    sum+=input[z+c.Nz*column]*std::cos(pi*static_cast<double>(j*z)/static_cast<double>(c.Nz-1));
            } else if(j!=0 && j+1!=c.Nz) {
                for(std::size_t z=1;z+1<c.Nz;++z)
                    sum+=input[z+c.Nz*column]*std::sin(pi*static_cast<double>(j*z)/static_cast<double>(c.Nz-1));
            }
            coefficients[j]=2*sum/static_cast<double>(c.Nz-1);
            double multiplier=0;
            if(j<c.Nj && j!=0) {
                const double m=pi*static_cast<double>(j)/c.Lz;
                multiplier=integral ? (inputIsF?1.0:-1.0)/m :
                    (inputIsF?signsF[order-1]:signsG[order-1])*std::pow(m,static_cast<int>(order));
            }
            coefficients[j]*=multiplier;
        }
        for(std::size_t z=0;z<c.Nz;++z) {
            double sum=outputIsF ? .5*coefficients[0] : 0;
            for(std::size_t j=1;j<c.Nj;++j)
                sum+=coefficients[j]*(outputIsF ?
                    std::cos(pi*static_cast<double>(j*z)/static_cast<double>(c.Nz-1)) :
                    std::sin(pi*static_cast<double>(j*z)/static_cast<double>(c.Nz-1)));
            result[z+c.Nz*column]=sum;
        }
        if(integral && !inputIsF) {
            const double bottom=result[c.Nz*column];
            for(std::size_t z=0;z<c.Nz;++z) result[z+c.Nz*column]-=bottom;
        }
    }
    return result;
}
void testVerticalCalculus(bool hydrostatic,WVConstantNonlinearFluxSchedule schedule) {
    auto k=kernel(hydrostatic,schedule); const auto& c=k->descriptor().configuration();
    std::vector<double> unprepared(c.Nz),unpreparedOutput(c.Nz);
    require(k->applyVerticalCalculus({unprepared.data(),{c.Nz,1}},true,1,false,
        {unpreparedOutput.data(),{c.Nz,1}}).code==WVKernelStatusCode::unsupportedOperation,
        "vertical calculus ran before setup");
    require(bool(k->prepareMatlabPrimitives()),"vertical calculus setup");
    const auto stable=k->persistentBytes(),scratch=k->scratchBytes(),plans=k->metrics().planCount;
    for(const std::size_t columns:{std::size_t{1},std::size_t{5}}) {
        std::vector<double> input(c.Nz*columns),output(c.Nz*columns),original;
        for(std::size_t q=0;q<input.size();++q)
            input[q]=.37*std::sin(.41*static_cast<double>(q+1))+
                .19*std::cos(.23*static_cast<double>((q%9)+2));
        original=input;
        for(const bool inputIsF:{true,false}) {
            for(const unsigned order:{1U,2U,3U,4U}) {
                const auto expected=verticalCalculusOracle(input,c,columns,inputIsF,order,false);
                require(bool(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},inputIsF,order,false,
                    {output.data(),{c.Nz,columns}})),"constant vertical derivative");
                for(std::size_t q=0;q<output.size();++q)
                    require(std::abs(output[q]-expected[q])<2e-12*std::max(1.0,std::abs(expected[q])),
                        "constant vertical derivative oracle");
            }
            const auto expected=verticalCalculusOracle(input,c,columns,inputIsF,1,true);
            require(bool(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},inputIsF,1,true,
                {output.data(),{c.Nz,columns}})),"constant vertical integral");
            for(std::size_t q=0;q<output.size();++q)
                require(std::abs(output[q]-expected[q])<2e-11*std::max(1.0,std::abs(expected[q])),
                    "constant vertical integral oracle");
            if(!inputIsF) for(std::size_t column=0;column<columns;++column)
                require(std::abs(output[c.Nz*column])<2e-13,"G integral bottom was not zero");
        }
        require(input==original,"vertical calculus mutated its input");
        require(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},true,0,false,
            {output.data(),{c.Nz,columns}}).code==WVKernelStatusCode::invalidConfiguration,
            "zero vertical derivative accepted");
        require(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},true,5,false,
            {output.data(),{c.Nz,columns}}).code==WVKernelStatusCode::invalidConfiguration,
            "fifth vertical derivative accepted");
        require(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},true,2,true,
            {output.data(),{c.Nz,columns}}).code==WVKernelStatusCode::invalidConfiguration,
            "higher vertical integral accepted");
        require(k->applyVerticalCalculus({input.data(),{c.Nz,columns}},true,1,false,
            {input.data(),{c.Nz,columns}}).code==WVKernelStatusCode::overlappingArrays,
            "vertical calculus alias accepted");
        require(k->applyVerticalCalculus({input.data(),{c.Nz-1,columns}},true,1,false,
            {output.data(),{c.Nz,columns}}).code==WVKernelStatusCode::invalidShape,
            "vertical calculus shape accepted");
    }
    const auto S=c.Nj*k->descriptor().Nkl(); std::vector<WVComplex64> stateStorage(3*S);
    const WVShape2D stateShape{c.Nj,k->descriptor().Nkl()};
    WVState state{0,0,{{stateStorage.data(),stateShape},{stateStorage.data()+S,stateShape},{stateStorage.data()+2*S,stateShape}}};
    require(bool(k->beginStateEvaluation(state)),"begin active state for vertical calculus");
    std::vector<double> stateInput(c.Nz);
    require(k->applyVerticalCalculus({stateInput.data(),{c.Nz,1}},true,1,false,
        {reinterpret_cast<double*>(stateStorage.data()),{c.Nz,1}}).code==WVKernelStatusCode::overlappingArrays,
        "vertical calculus overwrote active state");
    require(bool(k->endStateEvaluation()),"end active state for vertical calculus");
    require(k->persistentBytes()==stable && k->scratchBytes()==scratch && k->metrics().planCount==plans,
        "vertical calculus execution changed storage");
}
void testVertical(bool hydrostatic,WVConstantNonlinearFluxSchedule schedule) {
    auto k=kernel(hydrostatic,schedule); const auto& d=k->descriptor(); const auto& c=d.configuration();
    const auto beforePlans=k->metrics().planCount,beforeBytes=k->persistentBytes(),beforeScratch=k->scratchBytes();
    require(k->applyVertical(WVStratifiedModalOperator::projectF,{},{}).code==WVKernelStatusCode::unsupportedOperation,"raw vertical ran before setup");
    require(k->persistentBytes()==beforeBytes && k->scratchBytes()==beforeScratch && k->metrics().planCount==beforePlans,"unprepared primitive changed kernel storage");
    require(bool(k->prepareMatlabPrimitives()),"prepare MATLAB primitives");
    const auto preparedBytes=k->persistentBytes(),preparedScratch=k->scratchBytes(),preparedPlans=k->metrics().planCount;
    require(preparedPlans>beforePlans && bool(k->prepareMatlabPrimitives()) && k->persistentBytes()==preparedBytes && k->scratchBytes()==preparedScratch && k->metrics().planCount==preparedPlans,"MATLAB setup was not idempotent");
    std::vector<WVComplex64> grid(c.Nz*d.Nkl());
    for (std::size_t q=0;q<grid.size();++q) grid[q]={.07*static_cast<double>(q+1),-.03*static_cast<double>((q%11)+1)};
    std::vector<WVComplex64> modal(c.Nj*d.Nkl()),reconstructed(c.Nz*d.Nkl());
    const auto& m=d.verticalModes();
    using Op=WVStratifiedModalOperator;
    for (const auto op:{Op::projectF,Op::projectG,Op::projectFw,Op::projectGw}) {
        require(bool(k->applyVertical(op,{grid.data(),{c.Nz,d.Nkl()}},{modal.data(),{c.Nj,d.Nkl()}})),"constant projection");
        const bool sine=op==Op::projectG || op==Op::projectGw,wave=op==Op::projectFw || op==Op::projectGw;
        for (std::size_t column=0;column<d.Nkl();++column) for (std::size_t j=0;j<c.Nj;++j) {
            const double scale=sine ? (wave?m.gWaveScale[j]:m.Gg[j]) : (wave?m.fWaveScale[j+c.Nj*column]:m.Fg[j]);
            near(modal[j+c.Nj*column],(sine?dst(grid,column,c.Nz,j):dct(grid,column,c.Nz,j))/scale,1e-8,"constant projection oracle");
        }
    }
    for (std::size_t q=0;q<modal.size();++q) modal[q]={.04*static_cast<double>(q+2),.025*static_cast<double>((q%7)+1)};
    for (const auto op:{Op::reconstructF,Op::reconstructG,Op::reconstructFw,Op::reconstructGw}) {
        require(bool(k->applyVertical(op,{modal.data(),{c.Nj,d.Nkl()}},{reconstructed.data(),{c.Nz,d.Nkl()}})),"constant reconstruction");
        const bool sine=op==Op::reconstructG || op==Op::reconstructGw,wave=op==Op::reconstructFw || op==Op::reconstructGw;
        std::vector<WVComplex64> scaled=modal;
        for (std::size_t column=0;column<d.Nkl();++column) for (std::size_t j=0;j<c.Nj;++j) {
            const double factor=sine ? (wave?m.gWaveScale[j]:m.Gg[j]) : (wave?m.fWaveScale[j+c.Nj*column]:m.Fg[j]);
            scaled[j+c.Nj*column]=value(value(scaled[j+c.Nj*column])*factor);
        }
        for (std::size_t column=0;column<d.Nkl();++column) for (std::size_t z=0;z<c.Nz;++z)
            near(reconstructed[z+c.Nz*column],sine?inverseDst(scaled,column,c.Nj,z):inverseDct(scaled,column,c.Nj,z),3e-12,"constant reconstruction oracle");
    }
    std::size_t zero=0; while (zero<d.Nkl() && d.fourierModes()[zero].Kh!=0) ++zero; require(zero<d.Nkl(),"zero mode absent");
    std::vector<WVComplex64> column(grid.begin()+zero*c.Nz,grid.begin()+(zero+1)*c.Nz),columnOut(c.Nj);
    require(bool(k->applyVerticalColumn(Op::projectFw,zero,{column.data(),{c.Nz,1}},{columnOut.data(),{c.Nj,1}})),"Fio column projection");
    for (std::size_t j=0;j<c.Nj;++j) near(columnOut[j],dct(column,0,c.Nz,j)/m.fWaveScale[j+c.Nj*zero],1e-8,"Fio oracle including j0");
    const auto selected=d.Nkl()-1; column.assign(grid.begin()+selected*c.Nz,grid.begin()+(selected+1)*c.Nz);
    require(bool(k->applyVerticalColumn(Op::projectG,selected,{column.data(),{c.Nz,1}},{columnOut.data(),{c.Nj,1}})),"selected G column projection");
    for (std::size_t j=0;j<c.Nj;++j) near(columnOut[j],dst(column,0,c.Nz,j)/m.Gg[j],1e-8,"selected G column oracle");
    std::vector<WVComplex64> columnGrid(c.Nz);
    require(bool(k->applyVerticalColumn(Op::reconstructG,selected,{columnOut.data(),{c.Nj,1}},{columnGrid.data(),{c.Nz,1}})),"selected G column reconstruction");
    near(columnGrid.front(),{},2e-13,"G reconstruction bottom endpoint");
    near(columnGrid.back(),{},2e-13,"G reconstruction top endpoint");
    std::vector<WVComplex64> balanced(c.Nj*d.Nkl());
    require(bool(k->applyVertical(Op::balancedGToWaveG,{modal.data(),{c.Nj,d.Nkl()}},{balanced.data(),{c.Nj,d.Nkl()}})),"balanced G-to-wave-G");
    for (std::size_t q=0;q<balanced.size();++q) { const auto j=q%c.Nj; near(balanced[q],value(modal[q])*(m.Gg[j]/m.gWaveScale[j]),2e-13,"balanced G-to-wave-G oracle"); }
    require(k->applyVertical(Op::projectF,{grid.data(),{c.Nz,d.Nkl()}},{grid.data(),{c.Nj,d.Nkl()}}).code==WVKernelStatusCode::overlappingArrays,"vertical alias accepted");
    std::vector<WVComplex64> stateStorage(3*c.Nj*d.Nkl()); const WVShape2D stateShape{c.Nj,d.Nkl()};
    WVState state{0,0,{{stateStorage.data(),stateShape},{stateStorage.data()+c.Nj*d.Nkl(),stateShape},{stateStorage.data()+2*c.Nj*d.Nkl(),stateShape}}};
    require(bool(k->beginStateEvaluation(state)),"begin constant active state");
    require(k->applyVertical(Op::projectF,{grid.data(),{c.Nz,d.Nkl()}},{stateStorage.data(),stateShape}).code==WVKernelStatusCode::overlappingArrays,"vertical primitive overwrote active state");
    require(bool(k->endStateEvaluation()),"end constant active state");
    require(bool(k->applyVerticalColumn(Op::projectFw,zero,{column.data(),{c.Nz,1}},{columnOut.data(),{c.Nj,1}})),"vertical primitive did not recover after rejection");
    require(k->persistentBytes()==preparedBytes && k->scratchBytes()==preparedScratch && k->metrics().planCount==preparedPlans,"vertical execution changed storage");
}
double coordinate(std::size_t i,std::size_t n) { return 2*pi*static_cast<double>(i)/static_cast<double>(n); }
std::complex<double> dft(const std::vector<double>& a,const WVTransformConstantStratificationConfiguration& c,std::size_t z,std::int64_t k,std::int64_t l) {
    std::complex<double> sum; const auto P=c.Nx*c.Ny;
    for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x)
        sum+=a[x+c.Nx*y+P*z]*std::exp(std::complex<double>(0,-k*coordinate(x,c.Nx)-l*coordinate(y,c.Ny)));
    return sum/static_cast<double>(P);
}
double derivativeTerm(double amplitude,unsigned order,double mode,double phase,bool sine) {
    const auto v=std::pow(std::complex<double>(0,mode),static_cast<int>(order))*std::exp(std::complex<double>(0,mode*phase));
    return amplitude*(sine?v.imag():v.real());
}
void testHorizontal(bool hydrostatic,WVConstantNonlinearFluxSchedule schedule) {
    auto k=kernel(hydrostatic,schedule); require(bool(k->prepareMatlabPrimitives()),"horizontal setup");
    const auto& d=k->descriptor(); const auto& c=d.configuration(); const auto R=c.Nx*c.Ny*c.Nz,H=c.Nz*d.Nkl();
    const auto stable=k->persistentBytes(),scratch=k->scratchBytes(),plans=k->metrics().planCount;
    std::vector<double> field(R);
    for (std::size_t z=0;z<c.Nz;++z) for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x)
        field[x+c.Nx*y+c.Nx*c.Ny*z]=1+.1*z+.7*std::cos(coordinate(x,c.Nx))-.4*std::sin(2*coordinate(y,c.Ny))+.11*std::cos(4*coordinate(x,c.Nx))+.09*std::cos(3*coordinate(y,c.Ny));
    std::vector<WVComplex64> spectrum(H); require(bool(k->horizontalForward({field.data(),{c.Nx,c.Ny,c.Nz}},{spectrum.data(),{c.Nz,d.Nkl()}})),"horizontal forward");
    for (std::size_t mode=0;mode<d.Nkl();++mode) for (std::size_t z=0;z<c.Nz;++z)
        near(spectrum[z+c.Nz*mode],dft(field,c,z,d.fourierModes()[mode].kMode,d.fourierModes()[mode].lMode),4e-12,"retained DFT oracle");
    std::vector<double> restored(R); require(bool(k->horizontalInverse({spectrum.data(),{c.Nz,d.Nkl()}},{restored.data(),{c.Nx,c.Ny,c.Nz}})),"horizontal inverse");
    std::vector<double> derivative(R),first(R);
    for (unsigned order:{1U,2U,4U}) for (bool xDirection:{true,false}) {
        require(bool(k->differentiateHorizontal({field.data(),{c.Nx,c.Ny,c.Nz}},xDirection,order,{derivative.data(),{c.Nx,c.Ny,c.Nz}})),"constant horizontal derivative");
        for (std::size_t z=0;z<c.Nz;++z) for (std::size_t y=0;y<c.Ny;++y) for (std::size_t x=0;x<c.Nx;++x) {
            double expected=0;
            if (xDirection) { expected+=derivativeTerm(.7,order,1,coordinate(x,c.Nx),false); if (!(order%2)) expected+=derivativeTerm(.11,order,-4,coordinate(x,c.Nx),false); }
            else { expected+=derivativeTerm(-.4,order,2,coordinate(y,c.Ny),true); if (!(order%2)) expected+=derivativeTerm(.09,order,-3,coordinate(y,c.Ny),false); }
            require(std::abs(derivative[x+c.Nx*y+c.Nx*c.Ny*z]-expected)<4e-10,"constant derivative oracle");
        }
        if (order==1) { require(bool(k->differentiateHorizontal({field.data(),{c.Nx,c.Ny,c.Nz}},xDirection,{first.data(),{c.Nx,c.Ny,c.Nz}})),"first derivative overload"); require(std::memcmp(first.data(),derivative.data(),R*sizeof(double))==0,"order-one path changed"); }
    }
    require(k->differentiateHorizontal({field.data(),{c.Nx,c.Ny,c.Nz}},true,0,{derivative.data(),{c.Nx,c.Ny,c.Nz}}).code==WVKernelStatusCode::invalidConfiguration,"zero derivative accepted");
    require(k->differentiateHorizontal({field.data(),{c.Nx,c.Ny,c.Nz}},true,2,{field.data(),{c.Nx,c.Ny,c.Nz}}).code==WVKernelStatusCode::overlappingArrays,"derivative alias accepted");
    require(k->persistentBytes()==stable && k->scratchBytes()==scratch && k->metrics().planCount==plans,"horizontal execution changed storage");
}
void testArbitraryHorizontal(std::size_t nx,std::size_t ny,std::size_t nz,
    bool hydrostatic,bool antialias,WVConstantNonlinearFluxSchedule schedule) {
    auto c=configuration(hydrostatic); c.Nx=nx; c.Ny=ny; c.Nz=nz;
    c.Nj=antialias && nz-1>3 ? 2*(nz-1)/3 : nz-1;
    c.shouldAntialias=antialias;
    auto k=kernel(c,schedule); require(bool(k->prepareMatlabPrimitives()),"arbitrary horizontal setup");
    const auto& d=k->descriptor(); const auto R=nx*ny*nz,H=nz*d.Nkl();
    std::vector<double> field(R);
    for(std::size_t q=0;q<R;++q)
        field[q]=std::sin(.173*static_cast<double>(q+1))+
            .31*std::cos(.071*static_cast<double>((q+3)*(q%7+1)));
    std::vector<WVComplex64> spectrum(H);
    require(bool(k->horizontalForward({field.data(),{nx,ny,nz}},{spectrum.data(),{nz,d.Nkl()}})),"arbitrary horizontal forward");
    for(std::size_t mode=0;mode<d.Nkl();++mode) for(std::size_t z=0;z<nz;++z)
        near(spectrum[z+nz*mode],dft(field,c,z,d.fourierModes()[mode].kMode,
            d.fourierModes()[mode].lMode),5e-13,"arbitrary retained DFT oracle");
    std::vector<double> restored(R);
    require(bool(k->horizontalInverse({spectrum.data(),{nz,d.Nkl()}},{restored.data(),{nx,ny,nz}})),"arbitrary horizontal inverse");
    for(std::size_t z=0;z<nz;++z) for(std::size_t y=0;y<ny;++y) for(std::size_t x=0;x<nx;++x) {
        double expected=0;
        for(std::size_t mode=0;mode<d.Nkl();++mode) {
            const auto key=d.fourierModes()[mode];
            const auto coefficient=value(spectrum[z+nz*mode]);
            const auto phase=std::exp(std::complex<double>(0,key.kMode*coordinate(x,nx)+key.lMode*coordinate(y,ny)));
            const bool self=(2*key.kMode)%static_cast<std::int64_t>(nx)==0 &&
                (2*key.lMode)%static_cast<std::int64_t>(ny)==0;
            expected+=(self ? 1.0 : 2.0)*std::real(coefficient*phase);
        }
        require(std::abs(restored[x+nx*y+nx*ny*z]-expected)<5e-12,"arbitrary retained inverse oracle");
    }
    std::vector<WVComplex64> general(H);
    for(std::size_t q=0;q<H;++q)
        general[q]={.019*static_cast<double>(q+1),-.013*static_cast<double>((q%17)+1)};
    require(bool(k->horizontalInverse({general.data(),{nz,d.Nkl()}},{restored.data(),{nx,ny,nz}})),"general complex horizontal inverse");
    for(std::size_t z=0;z<nz;++z) for(std::size_t y=0;y<ny;++y) for(std::size_t x=0;x<nx;++x) {
        double expected=0;
        for(std::size_t mode=0;mode<d.Nkl();++mode) {
            const auto key=d.fourierModes()[mode];
            auto coefficient=value(general[z+nz*mode]);
            const bool self=(2*key.kMode)%static_cast<std::int64_t>(nx)==0 &&
                (2*key.lMode)%static_cast<std::int64_t>(ny)==0;
            if(self) coefficient=std::complex<double>(coefficient.real(),0);
            const auto phase=std::exp(std::complex<double>(0,key.kMode*coordinate(x,nx)+key.lMode*coordinate(y,ny)));
            expected+=(self ? 1.0 : 2.0)*std::real(coefficient*phase);
        }
        require(std::abs(restored[x+nx*y+nx*ny*z]-expected)<5e-12,"general complex retained inverse oracle");
    }
}
void testSeparateProjection(bool hydrostatic,WVConstantNonlinearFluxSchedule schedule) {
    auto k=kernel(hydrostatic,schedule); const auto& d=k->descriptor(); const auto& c=d.configuration(); const auto R=c.Nx*c.Ny*c.Nz,S=c.Nj*d.Nkl(); const auto channels=hydrostatic?3U:4U;
    std::vector<double> fields(channels*R); for (std::size_t i=0;i<fields.size();++i) fields[i]=std::sin(.013*(i+1));
    std::vector<WVComplex64> ap(S),am(S),a0(S),bp(S),bm(S),b0(S);
    WVMutableCoefficients a{{ap.data(),{c.Nj,d.Nkl()}},{am.data(),{c.Nj,d.Nkl()}},{a0.data(),{c.Nj,d.Nkl()}}};
    WVMutableCoefficients b{{bp.data(),{c.Nj,d.Nkl()}},{bm.data(),{c.Nj,d.Nkl()}},{b0.data(),{c.Nj,d.Nkl()}}};
    WVRealFieldBundleConstView bundle{fields.data(),{c.Nx,c.Ny,c.Nz,channels}};
    if (hydrostatic) { require(bool(k->transformUVEtaToWaveVortex(bundle,2,0,a)),"bundled H projection"); require(bool(k->transformUVEtaToWaveVortex({fields.data(),{c.Nx,c.Ny,c.Nz}},{fields.data()+R,{c.Nx,c.Ny,c.Nz}},{fields.data()+2*R,{c.Nx,c.Ny,c.Nz}},2,0,b)),"separate H projection"); }
    else { require(bool(k->transformUVWEtaToWaveVortex(bundle,2,0,a)),"bundled NH projection"); require(bool(k->transformUVWEtaToWaveVortex({fields.data(),{c.Nx,c.Ny,c.Nz}},{fields.data()+R,{c.Nx,c.Ny,c.Nz}},{fields.data()+2*R,{c.Nx,c.Ny,c.Nz}},{fields.data()+3*R,{c.Nx,c.Ny,c.Nz}},2,0,b)),"separate NH projection"); }
    require(std::memcmp(ap.data(),bp.data(),S*sizeof(WVComplex64))==0 && std::memcmp(am.data(),bm.data(),S*sizeof(WVComplex64))==0 && std::memcmp(a0.data(),b0.data(),S*sizeof(WVComplex64))==0,"separate projection changed arithmetic");
}
}
int main() {
    for (const bool hydrostatic:{true,false}) for (const auto schedule:{WVConstantNonlinearFluxSchedule::frozenStreamed,WVConstantNonlinearFluxSchedule::compactCandidate}) {
        testVertical(hydrostatic,schedule); testVerticalCalculus(hydrostatic,schedule);
        testHorizontal(hydrostatic,schedule); testSeparateProjection(hydrostatic,schedule);
        for(const bool antialias:{false,true}) {
            testArbitraryHorizontal(6,6,5,hydrostatic,antialias,schedule);
            testArbitraryHorizontal(8,6,7,hydrostatic,antialias,schedule);
        }
    }
    std::cout<<"Constant MATLAB primitive contracts passed\n";
}
