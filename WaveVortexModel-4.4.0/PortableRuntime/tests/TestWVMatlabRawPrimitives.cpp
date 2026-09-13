#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WVReferenceFFTEngine.hpp"
#include "WVBoussinesqModalTestFixture.hpp"
#include "WVStratifiedModalTestFixture.hpp"

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using namespace wavevortex::test_fixture;

namespace {
constexpr double pi=3.141592653589793238462643383279502884;
void hydrostaticFixture(const std::filesystem::path& path) {
    fixture(path);
    File file(path);
    for (const auto* name : {"WVTransform","AnnotatedClass"})
        nc(nc_put_att_text(file.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformHydrostatic"),"WVTransformHydrostatic"));
}
std::shared_ptr<const WVStratifiedModalRecord> read(const std::filesystem::path& path) {
    std::shared_ptr<const WVStratifiedModalRecord> record;
    const auto status=WVStratifiedModalReader::read(path.string(),record);
    require(static_cast<bool>(status),status.message.c_str());
    return record;
}

template <typename Kernel>
std::unique_ptr<Kernel> makeKernel(const std::shared_ptr<const WVStratifiedModalRecord>& record,
    WVVariableExecutionOptions options = {}) {
    std::unique_ptr<Kernel> kernel;
    const auto status=Kernel::create(record,std::make_unique<WVReferenceFFTEngine>(),kernel,
        WVCreateScalarMatrixBackend,options);
    require(static_cast<bool>(status),status.message.c_str());
    return kernel;
}

template <typename Kernel>
void horizontalRoundTripAndDerivatives(Kernel& kernel,const WVStratifiedModalGeometry& geometry) {
    const auto volume=kernel.spatialShape();
    const WVShape2D retained{geometry.Nz,geometry.Nkl};
    const auto count=volume.elementCount();
    std::vector<double> field(count),restored(count),derivative(count);
    for (std::size_t z=0;z<geometry.Nz;++z) for (std::size_t y=0;y<geometry.Ny;++y) for (std::size_t x=0;x<geometry.Nx;++x) {
        const auto i=x+geometry.Nx*(y+geometry.Ny*z);
        field[i]=std::cos(2*pi*x/geometry.Nx)+.4*std::sin(2*pi*y/geometry.Ny);
    }
    std::vector<WVComplex64> spectrum(retained.elementCount());
    require(static_cast<bool>(kernel.horizontalForward({field.data(),volume},{spectrum.data(),retained})),"Raw horizontal forward failed.");
    require(static_cast<bool>(kernel.horizontalInverse({spectrum.data(),retained},{restored.data(),volume})),"Raw horizontal inverse failed.");
    for (std::size_t i=0;i<count;++i) require(std::abs(restored[i]-field[i])<1e-10,"Raw horizontal round trip changed values.");

    std::vector<WVComplex64> general(retained.elementCount()),projected(general.size());
    for (std::size_t mode=0;mode<geometry.Nkl;++mode) for (std::size_t z=0;z<geometry.Nz;++z) {
        const auto i=z+geometry.Nz*mode;
        general[i]={.13+.017*i,std::cos(.23*i+.4)};
        projected[i]=general[i];
        const auto& key=geometry.modes[mode];
        const auto selfK=key.k==0 || (geometry.Nx%2==0 &&
            (key.k==static_cast<std::int64_t>(geometry.Nx/2) || key.k==-static_cast<std::int64_t>(geometry.Nx/2)));
        const auto selfL=key.l==0 || (geometry.Ny%2==0 &&
            (key.l==static_cast<std::int64_t>(geometry.Ny/2) || key.l==-static_cast<std::int64_t>(geometry.Ny/2)));
        if (selfK && selfL) projected[i].imag=0;
    }
    const auto original=general;
    std::vector<double> projectedResult(count);
    require(static_cast<bool>(kernel.horizontalInverse({general.data(),retained},{restored.data(),volume})),
        "Raw horizontal inverse rejected general complex retained input.");
    require(static_cast<bool>(kernel.horizontalInverse({projected.data(),retained},{projectedResult.data(),volume})),
        "Raw horizontal inverse rejected projected retained input.");
    for (std::size_t i=0;i<general.size();++i)
        require(general[i].real==original[i].real && general[i].imag==original[i].imag,
            "Raw horizontal inverse changed its retained input.");
    for (std::size_t z=0;z<geometry.Nz;++z) for (std::size_t y=0;y<geometry.Ny;++y) for (std::size_t x=0;x<geometry.Nx;++x) {
        double expected=0;
        for (std::size_t mode=0;mode<geometry.Nkl;++mode) {
            const auto& key=geometry.modes[mode]; const auto value=projected[z+geometry.Nz*mode];
            const double phase=2*pi*(static_cast<double>(key.k)*x/geometry.Nx+
                static_cast<double>(key.l)*y/geometry.Ny);
            const auto selfK=key.k==0 || (geometry.Nx%2==0 &&
                (key.k==static_cast<std::int64_t>(geometry.Nx/2) || key.k==-static_cast<std::int64_t>(geometry.Nx/2)));
            const auto selfL=key.l==0 || (geometry.Ny%2==0 &&
                (key.l==static_cast<std::int64_t>(geometry.Ny/2) || key.l==-static_cast<std::int64_t>(geometry.Ny/2)));
            expected+=(selfK && selfL ? 1.0 : 2.0)*(value.real*std::cos(phase)-value.imag*std::sin(phase));
        }
        const auto i=x+geometry.Nx*(y+geometry.Ny*z);
        require(std::abs(restored[i]-expected)<1e-10,"Raw inverse differs from the direct Hermitian Fourier sum.");
        require(std::abs(restored[i]-projectedResult[i])<1e-12,
            "Raw inverse did not discard imaginary self-conjugate values like MATLAB symmetric inverse.");
    }
    const auto persistent=kernel.persistentBytes();
    for(std::size_t z=0;z<geometry.Nz;++z) for(std::size_t y=0;y<geometry.Ny;++y) for(std::size_t x=0;x<geometry.Nx;++x)
        field[x+geometry.Nx*(y+geometry.Ny*z)]=std::cos(2*pi*x/geometry.Nx)+std::cos(pi*x)+.4*std::sin(4*pi*y/geometry.Ny);
    for (unsigned order : {1U,2U,4U}) for (bool xDerivative : {true,false}) {
        require(static_cast<bool>(kernel.differentiateHorizontal({field.data(),volume},xDerivative,order,{derivative.data(),volume})),
            "Raw horizontal derivative failed.");
        for (std::size_t z=0;z<geometry.Nz;++z) for (std::size_t y=0;y<geometry.Ny;++y) for (std::size_t x=0;x<geometry.Nx;++x) {
            const auto i=x+geometry.Nx*(y+geometry.Ny*z);
            const double coordinate=2*pi*(xDerivative ? x : y)/(xDerivative ? geometry.Nx : geometry.Ny);
            const double fundamental=xDerivative ? 2*pi/geometry.Lx : 2*pi/geometry.Ly;
            const double nyquist=pi;
            const double wave=xDerivative ? std::real(std::pow(std::complex<double>(0,fundamental),static_cast<int>(order))*std::exp(std::complex<double>(0,coordinate))) : 0;
            const double nyq=xDerivative ? std::real(std::pow(std::complex<double>(0,nyquist*geometry.Nx/geometry.Lx),static_cast<int>(order))*std::exp(std::complex<double>(0,pi*x))) : 0;
            const double yWave=!xDerivative ? .4*std::imag(std::pow(std::complex<double>(0,2*fundamental),static_cast<int>(order))*std::exp(std::complex<double>(0,2*coordinate))) : 0;
            const double expected=wave+nyq+yWave;
            require(std::abs(derivative[i]-expected)<std::pow(xDerivative?pi*geometry.Nx/geometry.Lx:4*pi/geometry.Ly,order)*1e-10,"Raw horizontal derivative disagrees with analytic DFT.");
        }
    }
    require(kernel.persistentBytes()==persistent,"Horizontal primitives changed persistent storage.");
}

template<class Kernel>
void verticalOperators(Kernel& kernel,
    const std::shared_ptr<const WVStratifiedModalRecord>& source) {
    const auto& g=source->geometry();
    const std::array<WVStratifiedModalOperator,11> operations={
        WVStratifiedModalOperator::reconstructF,WVStratifiedModalOperator::projectF,
        WVStratifiedModalOperator::reconstructG,WVStratifiedModalOperator::projectG,
        WVStratifiedModalOperator::reconstructFw,WVStratifiedModalOperator::projectFw,
        WVStratifiedModalOperator::reconstructGw,WVStratifiedModalOperator::projectGw,
        WVStratifiedModalOperator::balancedGToWaveG,WVStratifiedModalOperator::projectWaveDivergence,
        WVStratifiedModalOperator::projectWaveVerticalVelocity};
    const std::array<std::size_t,11> inputRows={g.Nj,g.Nz,g.Nj,g.Nz,g.Nj,g.Nz,g.Nj,g.Nz,g.Nj,g.Nz,g.Nz};
    const std::array<std::size_t,11> outputRows={g.Nz,g.Nj,g.Nz,g.Nj,g.Nz,g.Nj,g.Nz,g.Nj,g.Nj,g.Nj,g.Nj};
    const std::array<const char*,11> inputFamilies={"F-modal","F-grid","G-modal","G-grid","Fw-modal","F-grid","Gw-modal","G-grid","G-modal","F-grid","G-grid"};
    const std::array<const char*,11> outputFamilies={"F-grid","F-modal","G-grid","G-modal","F-grid","Fw-modal","G-grid","Gw-modal","Gw-modal","Gw-modal","Gw-modal"};
    const auto persistent=kernel.persistentBytes();
    for (std::size_t index=0;index<(g.transformClass=="WVTransformBoussinesq"?operations.size():4);++index) {
        const WVShape2D inputShape{inputRows[index],g.Nkl},outputShape{outputRows[index],g.Nkl};
        std::vector<WVComplex64> input(inputShape.elementCount()),expected(outputShape.elementCount()),actual(outputShape.elementCount());
        for (std::size_t i=0;i<input.size();++i) input[i]={std::sin(.07*i+index),std::cos(.11*i-index)};
        WVComplexLayout inputLayout{inputShape.rows,inputShape.columns,1,inputShape.rows,WVComplexRepresentation::interleaved,inputFamilies[index],source->modeSetIdentity()};
        WVComplexLayout outputLayout{outputShape.rows,outputShape.columns,1,outputShape.rows,WVComplexRepresentation::interleaved,outputFamilies[index],source->modeSetIdentity()};
        std::unique_ptr<WVPreparedVerticalOperator> prepared;
        std::unique_ptr<WVVerticalMatrixBackend> backend;
        require(static_cast<bool>(WVCreateScalarMatrixBackend(backend)),"Matrix backend creation failed.");
        require(static_cast<bool>(source->prepareVertical(operations[index],inputLayout,outputLayout,std::move(backend),prepared)),"Boussinesq operator preparation failed.");
        std::unique_ptr<WVVerticalWorkspace> workspace;
        require(static_cast<bool>(prepared->createWorkspace(workspace)),"Boussinesq operator workspace failed.");
        require(static_cast<bool>(prepared->execute(*workspace,{input.data(),nullptr,nullptr,input.size()*sizeof(WVComplex64)},{expected.data(),nullptr,nullptr,expected.size()*sizeof(WVComplex64)})),"Boussinesq reference execution failed.");
        require(static_cast<bool>(kernel.applyVertical(operations[index],{input.data(),inputShape},{actual.data(),outputShape})),"Boussinesq raw operator failed.");
        for (std::size_t i=0;i<actual.size();++i) require(std::abs(actual[i].real-expected[i].real)<1e-12 && std::abs(actual[i].imag-expected[i].imag)<1e-12,"Boussinesq raw operator differs from reference.");
        std::vector<WVComplex64> alias(std::max(input.size(),actual.size()));
        require(kernel.applyVertical(operations[index],{alias.data(),inputShape},{alias.data(),outputShape}).code==WVKernelStatusCode::overlappingArrays,"Raw operator accepted overlapping input/output.");
        for(std::size_t column=0;column<g.Nkl;++column) {
        std::vector<WVComplex64> columnInput(inputRows[index]),columnOutput(outputRows[index]);
        for (std::size_t row=0;row<inputRows[index];++row) columnInput[row]=input[row+inputRows[index]*column];
        require(static_cast<bool>(kernel.applyVerticalColumn(operations[index],column,{columnInput.data(),{inputRows[index],1}},{columnOutput.data(),{outputRows[index],1}})),"Boussinesq raw column operator failed.");
        for (std::size_t row=0;row<outputRows[index];++row) {
            const auto reference=expected[row+outputRows[index]*column];
            require(std::abs(columnOutput[row].real-reference.real)<1e-12 && std::abs(columnOutput[row].imag-reference.imag)<1e-12,"Boussinesq column differs from full operator.");
        }
        }
    }
    require(kernel.persistentBytes()==persistent,"Boussinesq raw operators changed persistent storage.");
}

template<class Kernel>
std::vector<double> elementaryVertical(Kernel& kernel,const WVStratifiedModalGeometry& g,
    const std::vector<double>& input,unsigned operation) {
    std::vector<WVComplex64> grid(g.Nz),modal(g.Nj),result(g.Nz);
    for (std::size_t z=0;z<g.Nz;++z) grid[z]={input[z],0};
    const auto projection=operation==0 || operation==3 ? WVStratifiedModalOperator::projectF : WVStratifiedModalOperator::projectG;
    require(static_cast<bool>(kernel.applyVerticalColumn(projection,0,{grid.data(),{g.Nz,1}},{modal.data(),{g.Nj,1}})),
        "Vertical-calculus oracle projection failed.");
    for (std::size_t j=0;j<g.Nj;++j) {
        const double factor=operation==1 || operation==2 ? 1/g.h_0[j] : operation==3 ? g.h_0[j] : 1;
        modal[j].real*=factor; modal[j].imag*=factor;
    }
    const auto reconstruction=operation==1 || operation==4 ? WVStratifiedModalOperator::reconstructF : WVStratifiedModalOperator::reconstructG;
    require(static_cast<bool>(kernel.applyVerticalColumn(reconstruction,0,{modal.data(),{g.Nj,1}},{result.data(),{g.Nz,1}})),
        "Vertical-calculus oracle reconstruction failed.");
    std::vector<double> output(g.Nz);
    const auto bottom=result.front();
    for (std::size_t z=0;z<g.Nz;++z) {
        auto value=result[z].real;
        if (operation==0 || operation==2) value*=(-g.N2[z]/g.g);
        if (operation==4) value=(value-bottom.real)*(-g.g);
        output[z]=value;
    }
    return output;
}

template<class Kernel>
std::vector<double> verticalOracle(Kernel& kernel,const WVStratifiedModalGeometry& g,
    const std::vector<double>& input,std::size_t columns,bool inputIsF,unsigned order,bool integral) {
    std::vector<double> expected(input.size());
    for (std::size_t column=0;column<columns;++column) {
        std::vector<double> value(g.Nz);
        for (std::size_t z=0;z<g.Nz;++z) {
            value[z]=input[z+g.Nz*column];
            if (integral && !inputIsF) value[z]/=g.N2[z];
        }
        const auto apply=[&](unsigned operation) { value=elementaryVertical(kernel,g,value,operation); };
        if (integral) apply(inputIsF ? 3 : 4);
        else if (inputIsF) {
            apply(0); if (order>=3) apply(2); if (order==2 || order==4) apply(1);
        } else {
            if (order==1) apply(1); else { apply(2); if (order==3) apply(1); if (order==4) apply(2); }
        }
        for (std::size_t z=0;z<g.Nz;++z) expected[z+g.Nz*column]=value[z];
    }
    return expected;
}

template<class Kernel>
void verticalCalculus(Kernel& kernel,const WVStratifiedModalGeometry& g) {
    const auto persistent=kernel.persistentBytes();
    for (const auto columns : {std::size_t{1},g.Nkl+2}) {
        const WVShape2D shape{g.Nz,columns};
        std::vector<double> input(shape.elementCount()),actual(input.size());
        for (std::size_t i=0;i<input.size();++i) input[i]=std::sin(.19*i)+.03*i;
        const auto original=input;
        for (bool inputIsF : {true,false}) {
            for (unsigned order=1;order<=4;++order) {
                const auto expected=verticalOracle(kernel,g,input,columns,inputIsF,order,false);
                require(static_cast<bool>(kernel.applyVerticalCalculus({input.data(),shape},inputIsF,order,false,{actual.data(),shape})),
                    "Raw vertical derivative failed.");
                for (std::size_t i=0;i<actual.size();++i)
                    require(std::abs(actual[i]-expected[i])<1e-11*std::max(1.0,std::abs(expected[i])),
                        "Raw vertical derivative differs from the independent column oracle.");
            }
            const auto expected=verticalOracle(kernel,g,input,columns,inputIsF,1,true);
            require(static_cast<bool>(kernel.applyVerticalCalculus({input.data(),shape},inputIsF,1,true,{actual.data(),shape})),
                "Raw vertical integral failed.");
            for (std::size_t i=0;i<actual.size();++i)
                require(std::abs(actual[i]-expected[i])<1e-11*std::max(1.0,std::abs(expected[i])),
                    "Raw vertical integral differs from the independent column oracle.");
            if (!inputIsF) for (std::size_t column=0;column<columns;++column)
                require(std::abs(actual[g.Nz*column])<1e-12,"Raw G integral is not bottom-zero.");
        }
        for (std::size_t i=0;i<input.size();++i) require(input[i]==original[i],"Raw vertical calculus changed its input.");
    }
    std::vector<double> input(g.Nz*2),output(input.size()); const WVShape2D shape{g.Nz,2};
    require(kernel.applyVerticalCalculus({input.data(),{g.Nz-1,2}},true,1,false,{output.data(),shape}).code==WVKernelStatusCode::invalidShape,
        "Raw vertical calculus accepted incorrect input rows.");
    require(kernel.applyVerticalCalculus({input.data(),shape},true,1,false,{output.data(),{g.Nz,1}}).code==WVKernelStatusCode::invalidShape,
        "Raw vertical calculus accepted a mismatched output shape.");
    require(kernel.applyVerticalCalculus({input.data(),shape},true,0,false,{output.data(),shape}).code==WVKernelStatusCode::unsupportedOperation,
        "Raw vertical calculus accepted derivative order zero.");
    require(kernel.applyVerticalCalculus({input.data(),shape},true,2,true,{output.data(),shape}).code==WVKernelStatusCode::unsupportedOperation,
        "Raw vertical calculus accepted a higher-order integral.");
    require(kernel.applyVerticalCalculus({input.data(),shape},true,1,false,{input.data(),shape}).code==WVKernelStatusCode::overlappingArrays,
        "Raw vertical calculus accepted overlapping arrays.");
    require(kernel.applyVerticalCalculus({nullptr,shape},true,1,false,{output.data(),shape}).code==WVKernelStatusCode::invalidPointer,
        "Raw vertical calculus accepted null input storage.");
    const WVShape2D huge{g.Nz,std::numeric_limits<std::size_t>::max()};
    require(kernel.applyVerticalCalculus({input.data(),huge},true,1,false,{output.data(),huge}).code==WVKernelStatusCode::sizeOverflow,
        "Raw vertical calculus accepted an overflowing extent.");
    require(kernel.persistentBytes()==persistent,"Raw vertical calculus changed persistent storage.");
}

template<class Kernel>
void check(const std::shared_ptr<const WVStratifiedModalRecord>& record) {
    for(bool split:{false,true}) {
        WVVariableExecutionOptions options;
        if(split) { options.spectralSchedule=WVVariableSpectralSchedule::compactSplitFusedViews; options.streamedNonlinear=true; options.horizontalSchedule=WVRetainedHorizontalSchedule::streamingPrunedTile16; }
        auto kernel=makeKernel<Kernel>(record,options);
        if constexpr (std::is_same_v<Kernel,WVTransformStratifiedQGKernel>) {
            const auto& g=record->geometry(); const WVShape2D shape{g.Nz,1};
            std::vector<double> input(g.Nz),output(g.Nz); const auto baseline=kernel->persistentBytes();
            require(kernel->applyVerticalCalculus({input.data(),shape},true,1,false,{output.data(),shape}).code==WVKernelStatusCode::unsupportedOperation,
                "QG raw vertical calculus ran before explicit preparation.");
            require(kernel->persistentBytes()==baseline,"Unprepared QG vertical calculus changed storage.");
            require(static_cast<bool>(kernel->prepareMatlabPrimitives()),"QG MATLAB primitive preparation failed.");
            const auto prepared=kernel->persistentBytes();
            require(prepared>baseline,"QG MATLAB primitive preparation did not account for its spectral slot.");
            require(static_cast<bool>(kernel->prepareMatlabPrimitives()) && kernel->persistentBytes()==prepared,
                "QG MATLAB primitive preparation is not idempotent.");
        }
        verticalOperators(*kernel,record); verticalCalculus(*kernel,record->geometry());
        horizontalRoundTripAndDerivatives(*kernel,record->geometry());
    }
}
void families() {
    Temporary h; hydrostaticFixture(h.path); check<WVTransformHydrostaticKernel>(read(h.path));
    Temporary qg; fixture(qg.path); check<WVTransformStratifiedQGKernel>(read(qg.path));
    Temporary b; boussinesqFixture(b.path); check<WVTransformBoussinesqKernel>(read(b.path));
}
}
int main() {
    try { families(); std::cout << "raw MATLAB primitive tests passed\n"; return 0; }
    catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
