#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexKernel/WVTransformStratifiedQGKernel.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"
#include "WVNativeFFTWEngine.hpp"
#include "WVAccelerateMatrixBackend.hpp"
#include "nlohmann/json.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/resource.h>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;
using nlohmann::json;

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
template<class Status> void require(Status status) {
    if (!status) throw std::runtime_error(status.message);
}

std::size_t count(const char* value) {
    std::string text(value);
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Expected a nonnegative integer.");
    return std::stoull(text);
}

double scalarValue(const WVStratifiedModalGeometry& g, std::size_t i, std::size_t j, std::size_t k) {
    const double x = static_cast<double>(i)*g.Lx/g.Nx, y = static_cast<double>(j)*g.Ly/g.Ny, z = g.z.empty() ? (g.Nz==1 ? 0 : -static_cast<double>(k)*g.Lz/(g.Nz-1)) : g.z[k];
    const double px = 2.0 * pi / g.Lx, py = 2.0 * pi / g.Ly;
    const std::size_t mx = std::max<std::size_t>(1, g.Nx / 2 - 1);
    const std::size_t my = std::max<std::size_t>(1, g.Ny / 2 - 1);
    const double amplitude = 1.0 + .17 * z / g.Lz + .09 * std::sin(.7 * z / g.Lz);
    return amplitude * (std::sin(px*x) + .31*std::cos(2*px*x + py*y) +
        .23*std::sin(mx*px*x + my*py*y) + .11*std::cos((g.Nx/2)*px*x) +
        .19*std::cos(py*y) + .13*std::sin(px*x + (g.Ny/2)*py*y));
}

double scalarDerivative(const WVStratifiedModalGeometry& g, std::size_t i, std::size_t j, std::size_t k) {
    const double x = static_cast<double>(i)*g.Lx/g.Nx, y = static_cast<double>(j)*g.Ly/g.Ny, z = g.z.empty() ? (g.Nz==1 ? 0 : -static_cast<double>(k)*g.Lz/(g.Nz-1)) : g.z[k];
    const double px = 2.0 * pi / g.Lx, py = 2.0 * pi / g.Ly;
    const std::size_t mx = std::max<std::size_t>(1, g.Nx / 2 - 1);
    const std::size_t my = std::max<std::size_t>(1, g.Ny / 2 - 1);
    const double amplitude = 1.0 + .17 * z / g.Lz + .09 * std::sin(.7 * z / g.Lz);
    const auto derivativeMode = [](std::size_t mode, std::size_t n) { return (n%2==0 && mode==n/2) ? 0.0 : static_cast<double>(mode); };
    const double mx2=derivativeMode(2,g.Nx), mxx=derivativeMode(mx,g.Nx), myy=derivativeMode(my,g.Ny);
    const double nxq=derivativeMode(g.Nx/2,g.Nx), nyq=derivativeMode(g.Ny/2,g.Ny);
    const double qx = amplitude * (px*std::cos(px*x) - .31*mx2*px*std::sin(2*px*x + py*y) +
        .23*mxx*px*std::cos(mx*px*x + my*py*y) - .11*nxq*px*std::sin((g.Nx/2)*px*x) +
        .13*px*std::cos(px*x + (g.Ny/2)*py*y));
    const double qy = amplitude * (-.31*py*std::sin(2*px*x + py*y) +
        .23*myy*py*std::cos(mx*px*x + my*py*y) - .19*py*std::sin(py*y) + .13*nyq*py*std::cos(px*x + (g.Ny/2)*py*y));
    return -.7*qx + .4*qy;
}

template<class Kernel, class Execute>
json measure(Kernel& kernel, Execute execute, const WVStratifiedModalGeometry& g, const std::vector<double>& scalar,
    const std::vector<double>& velocity, std::vector<double>& output, std::size_t warmups, std::size_t samples,
    const std::string& schedule) {
    const auto persistentBytes=kernel.persistentBytes();
    const WVShape3D shape{g.Nx,g.Ny,g.Nz};
    const WVRealVolumeConstView scalarView{scalar.data(),shape};
    const WVRealFieldBundleConstView velocityView{velocity.data(),{g.Nx,g.Ny,g.Nz,3}};
    const WVRealVolumeView outputView{output.data(),shape};
    for (std::size_t i=0; i<warmups; ++i) require(execute(scalarView,velocityView,outputView));
    std::vector<double> times(samples);
    for (double& elapsed : times) {
        const auto start = std::chrono::steady_clock::now();
        const auto status = execute(scalarView,velocityView,outputView);
        elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cerr << json{{"sampleSeconds",elapsed},{"success",bool(status)}}.dump() << std::endl;
        require(status);
    }
    if (persistentBytes!=kernel.persistentBytes()) throw std::runtime_error("Persistent storage changed during execution.");
    double maxError=0, maxReference=0, squaredError=0, squaredReference=0;
    for (std::size_t k=0; k<g.Nz; ++k) for (std::size_t j=0; j<g.Ny; ++j) for (std::size_t i=0; i<g.Nx; ++i) {
        const std::size_t index=i+g.Nx*(j+g.Ny*k);
        const double expected=scalarDerivative(g,i,j,k), error=output[index]-expected;
        if (scalar[index]!=scalarValue(g,i,j,k) || velocity[index]!=.7 || velocity[scalar.size()+index]!=-.4 || velocity[2*scalar.size()+index]!=0)
            throw std::runtime_error("Scalar execution changed caller input.");
        if (!std::isfinite(output[index])) throw std::runtime_error("Nonfinite scalar derivative output.");
        maxError=std::max(maxError,std::abs(error)); maxReference=std::max(maxReference,std::abs(expected));
        squaredError+=error*error; squaredReference+=expected*expected;
    }
    const double maxScaled=maxError/std::max(1.0,maxReference);
    const double relativeL2=std::sqrt(squaredError/std::max(squaredReference,std::numeric_limits<double>::min()));
    if (maxScaled>1e-10 || relativeL2>1e-10)
        throw std::runtime_error("Scalar derivative differs from the analytic oracle.");
    const auto& storage=kernel.storage();
    return { {"workspaceBytes",storage.workspaceBytes}, {"preparedBytes",storage.preparedBytes},
        {"realScratchBytes",storage.realScratchBytes}, {"spectralScratchBytes",storage.spectralScratchBytes},
        {"planBytesLowerBound",storage.planBytesLowerBound}, {"providerBytesLowerBound",storage.providerBytesLowerBound},
        {"inputPreserved",true}, {"samplesSeconds",times}, {"horizontalSchedule",schedule},
        {"persistentBytes",persistentBytes}, {"maximumScaleNormalizedError",maxError/std::max(1.0,maxReference)},
        {"relativeL2Error",relativeL2}, {"analyticComparison",{{"passed",true},{"maximumScaleNormalizedError",maxScaled},{"relativeL2Error",relativeL2}}} };
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 7) throw std::runtime_error("Usage: worker INPUT frozen|pruned-streamed WORKERS WARMUPS SAMPLES OUTPUT");
        const std::string selection=argv[2];
        if (selection!="frozen" && selection!="pruned-streamed") throw std::runtime_error("Unknown derivative selection.");
        const auto workers=count(argv[3]), warmups=count(argv[4]), samples=count(argv[5]);
        if (!workers || !samples || std::filesystem::exists(argv[6])) throw std::runtime_error("Invalid counts or existing output.");
        std::shared_ptr<const WVExtensionCatalog> catalog; require(makeBuiltInExtensionCatalog(catalog));
        WVCheckpoint checkpoint; require(WVCheckpointReader::read(argv[1],*catalog,checkpoint));
        if (!checkpoint.stratifiedModalSource) throw std::runtime_error("A variable-stratification source is required.");
        const auto& g=checkpoint.stratifiedModalSource->geometry();
        if (g.Nx<4 || g.Ny<4 || g.Nz<1) throw std::runtime_error("Derivative workload requires dimensions >=4 in x/y.");
        std::unique_ptr<WVFFTEngine> fft; require(WVFFTWEngine::create(1,fft));
        WVVariableExecutionOptions options;
        if (selection=="pruned-streamed") { options.horizontalSchedule=WVRetainedHorizontalSchedule::streamingPrunedTile16; options.horizontalWorkers=workers; options.streamedNonlinear=true; }
        const std::size_t R=g.Nx*g.Ny*g.Nz;
        std::vector<double> scalar(R), velocity(3*R), output(R);
        for (std::size_t k=0;k<g.Nz;++k) for (std::size_t j=0;j<g.Ny;++j) for (std::size_t i=0;i<g.Nx;++i) {
            const std::size_t index=i+g.Nx*(j+g.Ny*k); scalar[index]=scalarValue(g,i,j,k);
            velocity[index]=.7; velocity[R+index]=-.4; velocity[2*R+index]=0;
        }
        json result;
        if (g.transformClass=="WVTransformStratifiedQG") {
            std::unique_ptr<WVTransformStratifiedQGKernel> kernel;
            require(WVTransformStratifiedQGKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&](auto q,auto u,auto rhs){return kernel->advectScalarWithAdvectionFields(q,u,false,rhs);},g,scalar,velocity,output,warmups,samples,kernel->horizontalScheduleIdentifier());
        } else if (g.transformClass=="WVTransformHydrostatic") {
            std::unique_ptr<WVTransformHydrostaticKernel> kernel;
            require(WVTransformHydrostaticKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&](auto q,auto u,auto rhs){return kernel->advectScalarWithAdvectionFields(q,u,false,rhs,true);},g,scalar,velocity,output,warmups,samples,kernel->horizontalScheduleIdentifier());
        } else if (g.transformClass=="WVTransformBoussinesq") {
            std::unique_ptr<WVTransformBoussinesqKernel> kernel;
            require(WVTransformBoussinesqKernel::create(checkpoint.stratifiedModalSource,std::move(fft),kernel,WVCreateAccelerateMatrixBackend,options));
            result=measure(*kernel,[&](auto q,auto u,auto rhs){return kernel->advectScalarWithAdvectionFields(q,u,false,rhs,true);},g,scalar,velocity,output,warmups,samples,kernel->horizontalScheduleIdentifier());
        } else throw std::runtime_error("Unsupported transform family.");
        std::ofstream file(argv[6],std::ios::binary); file.write(reinterpret_cast<const char*>(output.data()),static_cast<std::streamsize>(output.size()*sizeof(double))); file.close();
        if (!file) throw std::runtime_error("Cannot write scalar derivative payload.");
        rusage usage{}; if (getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("Cannot measure peak RSS.");
        const auto identity=WVFFTWEngine::linkedLibraries();
        result["schema"]="wvm-variable-derivative-v1"; result["family"]=g.transformClass; result["selection"]=selection; result["grid"]={g.Nx,g.Ny,g.Nz}; result["warmups"]=warmups; result["horizontalWorkers"]=options.horizontalWorkers; result["matrixBackend"]="accelerate"; result["workload"]="scalar-horizontal"; result["antialias"]=false; result["xyOnly"]=g.transformClass!="WVTransformStratifiedQG"; result["provider"]={{"version",identity.version},{"baseLibrary",identity.baseLibrary},{"threadLibrary",identity.threadLibrary},{"fftThreads",1}}; result["peakRSSBytes"]=static_cast<std::size_t>(usage.ru_maxrss);
        std::cout << result.dump() << '\n'; return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
