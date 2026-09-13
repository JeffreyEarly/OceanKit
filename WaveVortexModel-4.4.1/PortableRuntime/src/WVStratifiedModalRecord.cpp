#include "WaveVortexRuntime/WVStratifiedModalRecord.hpp"
#include "WVNetCDF.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>
#include <tuple>

namespace wavevortex::runtime {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
WVCheckpointStatus invalid(const std::string& message, const std::string& name = "") {
    return {WVCheckpointStatusCode::invalidValue,message,"/"+name};
}
std::size_t product(std::size_t a, std::size_t b) {
    if (a && b > static_cast<std::size_t>(PTRDIFF_MAX)/a) throw std::overflow_error("Modal extent exceeds the addressable range.");
    return a*b;
}
bool close(double a, double b, double scale = 1) {
    return std::abs(a-b) <= 128*std::numeric_limits<double>::epsilon()*std::max({scale,std::abs(a),std::abs(b)});
}
WVCheckpointStatus matrixVariable(int file, const char* name, const char* fast, const char* slow,
    std::size_t rows, std::size_t columns, int& variable) {
    int code = nc_inq_varid(file,name,&variable);
    if (code == NC_ENOTVAR) return {WVCheckpointStatusCode::missingVariable,"Required numeric modal input is missing.",std::string("/")+name};
    if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal variable lookup",name);
    nc_type type; int rank = 0;
    code = nc_inq_var(file,variable,nullptr,&type,&rank,nullptr,nullptr);
    if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal variable inspection",name);
    if (type != NC_DOUBLE || rank != (slow ? 2 : 1)) return {WVCheckpointStatusCode::shapeMismatch,"Modal input must be a double with its declared rank.",name};
    std::array<int,2> dims{}; code = nc_inq_vardimid(file,variable,dims.data());
    if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal dimensions",name);
    for (int i = 0; i < rank; ++i) {
        const char* expected = slow && i == 0 ? slow : fast;
        int id; code = nc_inq_dimid(file,expected,&id);
        if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal dimension lookup",expected);
        std::size_t length; code = nc_inq_dimlen(file,id,&length);
        if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal dimension length",expected);
        if (dims[i] != id || length != (slow && i == 0 ? columns : rows))
            return {WVCheckpointStatusCode::shapeMismatch,"Modal matrix orientation or dimension identity disagrees with its schema.",name};
    }
    product(product(rows,columns),sizeof(double));
    return WVCheckpointStatus::ok();
}
WVCheckpointStatus vectorInput(int file, const char* name, const char* dimension, std::size_t size, std::vector<double>& output) {
    int variable; auto status = matrixVariable(file,name,dimension,nullptr,size,1,variable); if (!status) return status;
    double fill = NC_FILL_DOUBLE; int noFill = 0;
    const int fillCode = nc_inq_var_fill(file,variable,&noFill,&fill);
    if (fillCode != NC_NOERR) return detail::netcdfFailure(fillCode,"Modal fill metadata",name);
    output.resize(size); const int code = nc_get_var_double(file,variable,output.data());
    if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal vector read",name);
    for (const double value : output) if (!std::isfinite(value) || value == fill || value == NC_FILL_DOUBLE) return invalid("Modal vector contains a nonfinite or missing value.",name);
    return WVCheckpointStatus::ok();
}
WVCheckpointStatus scanMatrix(int file, const char* name, const WVStratifiedModalGeometry& g, std::vector<double>* output) {
    const bool inverse = std::string(name).find("inv") != std::string::npos;
    const bool G = name[0] == 'Q';
    const auto rows = inverse ? g.Nz : g.Nj, columns = inverse ? g.Nj : g.Nz;
    int variable; auto status = matrixVariable(file,name,inverse ? "z" : "j",inverse ? "j" : "z",rows,columns,variable); if (!status) return status;
    double fill = NC_FILL_DOUBLE; int noFill = 0;
    const int fillCode = nc_inq_var_fill(file,variable,&noFill,&fill);
    if (fillCode != NC_NOERR) return detail::netcdfFailure(fillCode,"Modal fill metadata",name);
    if (output) output->resize(product(rows,columns));
    std::array<double,4096> buffer{};
    for (std::size_t c = 0; c < columns; ++c) {
        double maximum = 0;
        for (std::size_t first = 0; first < rows; first += std::min(buffer.size(),rows-first)) {
            const std::size_t count = std::min(buffer.size(),rows-first);
            const std::size_t start[] = {c,first}, extent[] = {1,count};
            const int code = nc_get_vara_double(file,variable,start,extent,buffer.data());
            if (code != NC_NOERR) return detail::netcdfFailure(code,"Modal payload scan",name);
            for (std::size_t i = 0; i < count; ++i) {
                const auto r = first+i; const double value = buffer[i];
                if (!std::isfinite(value) || value == fill || value == NC_FILL_DOUBLE) return invalid("Modal matrix contains a nonfinite or missing value.",name);
                const auto z = inverse ? r : c, j = inverse ? c : r;
                if (G && (g.j[j] == 0 || z == 0 || z+1 == g.Nz) && value != 0)
                    return invalid("Rigid-lid G matrices require zero boundaries and a zero barotropic mode.",name);
                if (inverse && !G && g.j[j] == 0 && !close(value,1))
                    return invalid("The preconditioned barotropic F column must equal one.",name);
                maximum = std::max(maximum,std::abs(value));
                if (output) (*output)[r+rows*c] = value;
            }
        }
        if (inverse && !(G && g.j[c] == 0) && !close(maximum,1))
            return invalid("Reconstruction columns must use the persisted unit-maximum preconditioning.",name);
    }
    return WVCheckpointStatus::ok();
}
// Scan each declared column through the same bounded buffer used by shared
// matrices. NetCDF dimensions are reversed relative to MATLAB column-major data.
WVCheckpointStatus scanWaveArray(int file,const char* name,const WVStratifiedModalGeometry& g,std::vector<double>* output) {
    const std::string key=name;
    const bool inverse=key=="PFpmInv" || key=="QGpmInv";
    const bool cross=key=="QGwg", scale=key=="Ppm" || key=="Qpm", depth=key=="h_pm";
    const bool G=key.front()=='Q';
    const std::size_t rows=inverse ? g.Nz : g.Nj,columns=scale || depth ? (depth ? g.Nkl : g.K2unique.size()) : (cross ? g.Nj : inverse ? g.Nj : g.Nz);
    const auto groups=scale || depth ? 1 : g.K2unique.size();
    const int rank=scale || depth ? 2 : 3;
    const char* dimensions[3]={scale ? "K2unique" : depth ? "kl" : "K2unique",inverse || cross ? "j" : "z",inverse ? "z" : "j"};
    if (rank==2) dimensions[1]="j";
    const std::size_t lengths[3]={rank==3 ? groups : columns,rank==3 ? columns : rows,rows};
    int variable=-1; auto code=nc_inq_varid(file,name,&variable);
    if (code==NC_ENOTVAR) return {WVCheckpointStatusCode::missingVariable,"Required Boussinesq modal array is missing.",name};
    if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave array lookup",name);
    int actualRank=0; nc_type type; code=nc_inq_var(file,variable,nullptr,&type,&actualRank,nullptr,nullptr);
    if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave array metadata",name);
    if (type!=NC_DOUBLE || actualRank!=rank) return {WVCheckpointStatusCode::shapeMismatch,"Wave array type or rank disagrees with the persisted schema.",name};
    int ids[3]; code=nc_inq_vardimid(file,variable,ids);
    if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave array dimensions",name);
    for (int i=0;i<rank;++i) {
        int expected=-1; std::size_t length=0;
        code=nc_inq_dimid(file,dimensions[i],&expected); if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave dimension lookup",name);
        code=nc_inq_dimlen(file,expected,&length); if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave dimension length",name);
        if (ids[i]!=expected || length!=lengths[i]) return {WVCheckpointStatusCode::shapeMismatch,"Wave array orientation or dimension identity is invalid.",name};
    }
    const auto count=product(product(rows,columns),groups); product(count,sizeof(double));
    if (output) output->resize(count);
    double fill=NC_FILL_DOUBLE; int noFill=0; code=nc_inq_var_fill(file,variable,&noFill,&fill);
    if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave fill metadata",name);
    std::array<double,4096> buffer{};
    for (std::size_t group=0;group<groups;++group) for (std::size_t column=0;column<columns;++column) {
        double maximum=0;
        for (std::size_t first=0;first<rows;first+=std::min(rows-first,buffer.size())) {
            const auto n=std::min(rows-first,buffer.size());
            const std::size_t start[]={group,column,first},extent[]={1,1,n};
            code=nc_get_vara_double(file,variable,start+(rank==2),extent+(rank==2),buffer.data());
            if (code!=NC_NOERR) return detail::netcdfFailure(code,"Wave payload scan",name);
            for (std::size_t i=0;i<n;++i) {
                const auto row=first+i; const double value=buffer[i];
                if (!std::isfinite(value) || value==fill || value==NC_FILL_DOUBLE) return invalid("Wave array contains a nonfinite or missing value.",name);
                if (scale || depth) {
                    if (value<=0 || (g.j[row]==0 && value!=1)) return invalid("Wave scales/depths require positive values and unit barotropic sentinels.",name);
                } else if (cross) {
                    if ((g.j[row]==0 || g.j[column]==0) && value!=0) return invalid("Wave G cross-projection must exclude the barotropic mode.",name);
                } else {
                    const auto z=inverse ? row : column,j=inverse ? column : row;
                    if (G && (g.j[j]==0 || z==0 || z+1==g.Nz) && value!=0) return invalid("Wave G matrices require rigid-lid boundaries and zero barotropic mode.",name);
                    if (inverse && !G && g.j[j]==0 && !close(value,1)) return invalid("Wave F barotropic column must equal one.",name);
                }
                maximum=std::max(maximum,std::abs(value));
                if (output) (*output)[row+rows*(column+columns*group)]=value;
            }
        }
        if (inverse && !(G && g.j[column]==0) && !close(maximum,1)) return invalid("Wave reconstruction must retain unit-maximum preconditioning.",name);
    }
    return WVCheckpointStatus::ok();
}
WVCheckpointStatus inspectFile(int file, WVStratifiedModalInspection& result) {
    WVStratifiedModalInspection candidate; auto& g = candidate.geometry;
    auto status = detail::readTextAttribute(file,"model_version",g.modelVersion,"/"); if (!status) return status;
    if (g.modelVersion.rfind("4.",0) != 0) return {WVCheckpointStatusCode::unsupportedModelVersion,"Modal reader requires a v4 scientific record.","/@model_version"};
    status = detail::readTextAttribute(file,"WVTransform",g.transformClass,"/"); if (!status) return status;
    std::string annotated; status = detail::readTextAttribute(file,"AnnotatedClass",annotated,"/"); if (!status) return status;
    if (annotated != g.transformClass) return {WVCheckpointStatusCode::schemaMismatch,"Root transform identities disagree.","/"};
    if (g.transformClass != "WVTransformStratifiedQG" && g.transformClass != "WVTransformHydrostatic" && g.transformClass != "WVTransformBoussinesq")
        return {WVCheckpointStatusCode::unsupportedTransform,"This modal record reader supports v4 Stratified QG, Hydrostatic and Boussinesq records.","/"};
    struct Dimension { const char* name; std::size_t* value; };
    for (const auto& d : {Dimension{"x",&g.Nx},{"y",&g.Ny},{"z",&g.Nz},{"j",&g.Nj},{"kl",&g.Nkl}}) {
        status = detail::dimensionLength(file,d.name,*d.value,"/"); if (!status) return status;
        if (!*d.value) return invalid("Modal dimensions must be nonempty.",d.name);
    }
    if (g.Nx < 2 || g.Ny < 2 || g.Nz < 3 || g.Nj >= g.Nz || g.Nx > static_cast<std::size_t>(INT64_MAX) || g.Ny > static_cast<std::size_t>(INT64_MAX))
        return invalid("Invalid horizontal or rigid-lid retained vertical dimensions.");
    product(product(product(g.Nx,g.Ny),g.Nz),sizeof(double));
    product(product(g.Nj,g.Nkl),sizeof(WVComplex64));
    candidate.scientificMatrixBytes = product(product(g.Nz,g.Nj),4*sizeof(double));
    struct Scalar { const char* name; double* value; };
    for (const auto& s : {Scalar{"Lx",&g.Lx},{"Ly",&g.Ly},{"Lz",&g.Lz},{"g",&g.g},{"rho0",&g.rho0},{"latitude",&g.latitude},{"rotationRate",&g.rotationRate},{"planetaryRadius",&g.planetaryRadius}}) {
        status = detail::readDoubleScalar(file,s.name,*s.value,"/"); if (!status) return status;
        if (!std::isfinite(*s.value) || (s.value != &g.latitude && *s.value <= 0)) return invalid("Invalid physical scalar.",s.name);
    }
    if (std::abs(g.latitude) < 5 || std::abs(g.latitude) > 85) return invalid("Absolute latitude must lie in [5,85].","latitude");
    status = detail::readLogicalScalar(file,"shouldAntialias",g.shouldAntialias,"/"); if (!status) return status;
    struct Vector { const char* name; const char* dimension; std::size_t size; std::vector<double>* values; };
    for (const auto& v : {Vector{"x","x",g.Nx,&g.x},{"y","y",g.Ny,&g.y},{"z","z",g.Nz,&g.z},{"j","j",g.Nj,&g.j},
        {"k","kl",g.Nkl,&g.k},{"l","kl",g.Nkl,&g.l},{"N2","z",g.Nz,&g.N2},{"rho_nm0","z",g.Nz,&g.rho_nm0},{"dLnN2","z",g.Nz,&g.dLnN2},
        {"P0","j",g.Nj,&g.P0},{"Q0","j",g.Nj,&g.Q0},{"h_0","j",g.Nj,&g.h_0},{"z_int","z",g.Nz,&g.z_int}}) {
        status = vectorInput(file,v.name,v.dimension,v.size,*v.values); if (!status) return status;
    }
    for (std::size_t i = 0; i < g.Nx; ++i) if (!close(g.x[i],g.Lx*static_cast<double>(i)/g.Nx,g.Lx)) return invalid("x must use the periodic origin and spacing.","x");
    for (std::size_t i = 0; i < g.Ny; ++i) if (!close(g.y[i],g.Ly*static_cast<double>(i)/g.Ny,g.Ly)) return invalid("y must use the periodic origin and spacing.","y");
    // MATLAB WKB quadrature inversion places endpoints to interpolation accuracy.
    if (std::abs(g.z.front()+g.Lz) > 1e-10*g.Lz || std::abs(g.z.back()) > 1e-10*g.Lz) return invalid("Vertical grid must include bottom -Lz and surface zero.","z");
    long double integral = 0;
    for (std::size_t i = 0; i < g.Nz; ++i) {
        if ((i && g.z[i] <= g.z[i-1]) || g.N2[i] <= 0 || g.rho_nm0[i] <= 0 || (i && g.rho_nm0[i] > g.rho_nm0[i-1]))
            return invalid("Vertical coordinates and stable stratification are inconsistent.","z/N2/rho_nm0");
        integral += g.z_int[i];
    }
    if (std::abs(integral-g.Lz) > 1e-10*g.Lz) return invalid("Vertical quadrature must integrate a constant to Lz.","z_int");
    for (std::size_t i = 0; i < g.Nj; ++i) {
        if (g.j[i] != std::floor(g.j[i]) || g.j[i] < 0 || g.j[i] >= static_cast<double>(g.Nz-1) || (i && g.j[i] <= g.j[i-1]))
            return invalid("Retained vertical keys must be unique increasing rigid-lid mode indices.","j");
        if (g.P0[i] <= 0 || g.Q0[i] <= 0 || g.h_0[i] <= 0) return invalid("Preconditioners and equivalent depths must be positive.","P0/Q0/h_0");
    }
    if (g.j.front() != 0 || !close(g.P0.front(),1) || !close(g.Q0.front(),1) || !close(g.h_0.front(),1)) return invalid("The barotropic mode requires its documented unit sentinel values.","j/P0/Q0/h_0");
    std::vector<double> kl; status = vectorInput(file,"kl","kl",g.Nkl,kl); if (!status) return status;
    std::set<std::pair<std::int64_t,std::int64_t>> seen;
    const double cutoff = (2.0/3.0)*(2*pi*static_cast<double>(g.Nx/2)/g.Lx);
    for (std::size_t i = 0; i < g.Nkl; ++i) {
        if (kl[i] != static_cast<double>(i)) return invalid("kl must name the ordered retained columns.","kl");
        const double k = g.k[i]*g.Lx/(2*pi), l = g.l[i]*g.Ly/(2*pi);
        if (!std::isfinite(k) || !std::isfinite(l) || std::abs(std::round(k)) > static_cast<double>((g.Nx-1)/2) || std::abs(std::round(l)) > static_cast<double>((g.Ny-1)/2) || !close(k,std::round(k)) || !close(l,std::round(l)))
            return invalid("Horizontal keys are outside the Nyquist-excluding integer grid.","k/l");
        const auto ki = static_cast<std::int64_t>(std::round(k)), li = static_cast<std::int64_t>(std::round(l));
        if (li < 0 || (li == 0 && ki < 0) || !seen.insert({ki,li}).second) return invalid("Horizontal keys must be unique primary Hermitian representatives.","k/l");
        if (g.shouldAntialias && std::hypot(g.k[i],g.l[i]) > cutoff && !close(std::hypot(g.k[i],g.l[i]),cutoff,cutoff)) return invalid("Retained key violates the declared antialias mask.","k/l");
        g.modes.push_back({ki,li});
    }
    for (const char* name : {"PF0inv","QG0inv","PF0","QG0"}) { status = scanMatrix(file,name,g,nullptr); if (!status) return status; }
    if (g.transformClass=="WVTransformBoussinesq") {
        std::size_t groups=0; status=detail::dimensionLength(file,"K2unique",groups,"/"); if (!status) return status;
        if (!groups || groups>g.Nkl) return invalid("Wave group count must be nonempty and bounded by retained columns.","K2unique");
        status=vectorInput(file,"K2unique","K2unique",groups,g.K2unique); if (!status) return status;
        std::vector<double> membership; status=vectorInput(file,"iK2unique","kl",g.Nkl,membership); if (!status) return status;
        for (std::size_t i=0;i<groups;++i) if (g.K2unique[i]<0 || (i && g.K2unique[i]<=g.K2unique[i-1])) return invalid("Wave group radii must be nonnegative and strictly increasing.","K2unique");
        if (g.K2unique.front()!=0) return invalid("Wave groups must include the inertial zero wavenumber.","K2unique");
        std::vector<std::size_t> first(groups,g.Nkl); g.waveGroup.resize(g.Nkl);
        for (std::size_t i=0;i<g.Nkl;++i) {
            const auto value=membership[i];
            if (value<1 || value>static_cast<double>(groups) || value!=std::floor(value)) return invalid("Wave membership must be a valid one-based persisted group index.","iK2unique");
            const auto group=static_cast<std::size_t>(value)-1; g.waveGroup[i]=group; first[group]=std::min(first[group],i);
            const double radius=g.k[i]*g.k[i]+g.l[i]*g.l[i];
            if (!close(radius,g.K2unique[group],std::max(radius,g.K2unique[group]))) return invalid("Persisted wave membership disagrees with its horizontal wavenumber.","iK2unique");
        }
        if (g.waveGroup.front()!=0 || std::find(first.begin(),first.end(),g.Nkl)!=first.end()) return invalid("Wave groups require complete membership and a first inertial column.","iK2unique");
        struct WaveScale { const char* name; std::vector<double>* values; };
        const WaveScale scales[]={{"Ppm",&g.Ppm},{"Qpm",&g.Qpm},{"h_pm",&g.h_pm}};
        for (const auto& entry:scales) {
            status=scanWaveArray(file,entry.name,g,entry.values); if (!status) return status;
        }
        const double f=2*g.rotationRate*std::sin(g.latitude*pi/180);
        for (double N2:g.N2) if (N2==f*f) return invalid("Boussinesq wave projection is singular where N2 equals f squared.","N2");
        for (std::size_t i=0;i<g.Nkl;++i) for (std::size_t j=0;j<g.Nj;++j)
            if (g.h_pm[j+g.Nj*i]!=g.h_pm[j+g.Nj*first[g.waveGroup[i]]]) return invalid("Equivalent depths differ within one persisted wave group.","h_pm");
        for (const char* name:{"PFpmInv","QGpmInv","PFpm","QGpm","QGwg"}) { status=scanWaveArray(file,name,g,nullptr); if (!status) return status; }
        const auto waveBytes=product(product(product(4,g.Nz)+g.Nj,g.Nj),product(groups,sizeof(double)));
        if (waveBytes>static_cast<std::size_t>(PTRDIFF_MAX)-candidate.scientificMatrixBytes) throw std::overflow_error("Scientific matrix byte count overflow.");
        candidate.scientificMatrixBytes+=waveBytes;
    }
    result = std::move(candidate); return WVCheckpointStatus::ok();
}
}
WVCheckpointStatus WVStratifiedModalReader::inspect(const std::string& path, WVStratifiedModalInspection& result) {
    try {
        detail::WVNetCDFFile file; auto status = detail::WVNetCDFFile::openReadOnly(path,file); if (!status) return status;
        return inspectFile(file.id(),result);
    } catch (const std::bad_alloc&) { return {WVCheckpointStatusCode::allocationFailure,"Modal inspection allocation failed.",path}; }
      catch (const std::overflow_error& e) { return invalid(e.what()); }
}
WVCheckpointStatus WVStratifiedModalReader::read(const std::string& path, std::shared_ptr<const WVStratifiedModalRecord>& result) {
    try {
        detail::WVNetCDFFile file; auto status = detail::WVNetCDFFile::openReadOnly(path,file); if (!status) return status;
        WVStratifiedModalInspection inspection; status = inspectFile(file.id(),inspection); if (!status) return status;
        auto candidate = std::shared_ptr<WVStratifiedModalRecord>(new WVStratifiedModalRecord);
        WVStratifiedModalArrays arrays; arrays.geometry=std::move(inspection.geometry);
        struct Matrix { const char* name; std::vector<double>* values; };
        for (const auto& m : {Matrix{"PF0inv",&arrays.PF0inv},{"QG0inv",&arrays.QG0inv},{"PF0",&arrays.PF0},{"QG0",&arrays.QG0}}) {
            status = scanMatrix(file.id(),m.name,arrays.geometry,m.values); if (!status) return status;
        }
        if (!arrays.geometry.K2unique.empty()) for (const auto& m : {Matrix{"PFpmInv",&arrays.PFpmInv},{"QGpmInv",&arrays.QGpmInv},{"PFpm",&arrays.PFpm},{"QGpm",&arrays.QGpm},{"QGwg",&arrays.QGwg}}) {
            status=scanWaveArray(file.id(),m.name,arrays.geometry,m.values); if (!status) return status;
        }
        int functionVariable = -1;
        const int functionCode = nc_inq_varid(file.id(),"N2Function",&functionVariable);
        if (functionCode != NC_ENOTVAR) {
            if (functionCode != NC_NOERR) return detail::netcdfFailure(functionCode,"Function payload lookup","N2Function");
            int rank=0, dimension=-1; nc_type type;
            int code=nc_inq_var(file.id(),functionVariable,nullptr,&type,&rank,nullptr,nullptr);
            if (code != NC_NOERR) return detail::netcdfFailure(code,"Function payload metadata","N2Function");
            if (rank != 1 || type != NC_UBYTE) return invalid("N2Function must be an opaque uint8 vector.","N2Function");
            nc_type markerType; std::size_t markerLength=0; unsigned char marker=0;
            code=nc_inq_att(file.id(),functionVariable,"isFunctionHandleType",&markerType,&markerLength);
            if (code != NC_NOERR || markerType != NC_UBYTE || markerLength != 1)
                return invalid("N2Function requires its scalar function-handle marker.","N2Function");
            code=nc_get_att_uchar(file.id(),functionVariable,"isFunctionHandleType",&marker);
            if (code != NC_NOERR || marker != 1) return invalid("Invalid N2Function function-handle marker.","N2Function");
            code=nc_inq_vardimid(file.id(),functionVariable,&dimension); if (code != NC_NOERR) return detail::netcdfFailure(code,"Function payload dimension","N2Function");
            std::size_t size=0; code=nc_inq_dimlen(file.id(),dimension,&size); if (code != NC_NOERR) return detail::netcdfFailure(code,"Function payload length","N2Function");
            candidate->N2FunctionPayload_.resize(size);
            code=nc_get_var_uchar(file.id(),functionVariable,candidate->N2FunctionPayload_.data()); if (code != NC_NOERR) return detail::netcdfFailure(code,"Function payload read","N2Function");
        }
        const auto sourceStatus=WVOwnedStratifiedModalSource::create(std::move(arrays),candidate->source_);
        if (!sourceStatus) {
            if (sourceStatus.code==WVKernelStatusCode::allocationFailure)
                return {WVCheckpointStatusCode::allocationFailure,sourceStatus.message,path};
            return invalid(sourceStatus.message);
        }
        result = std::move(candidate); return WVCheckpointStatus::ok();
    } catch (const std::bad_alloc&) { return {WVCheckpointStatusCode::allocationFailure,"Modal record allocation failed.",path}; }
      catch (const std::overflow_error& e) { return invalid(e.what()); }
}
std::size_t WVStratifiedModalRecord::persistentBytes() const noexcept {
    return sizeof(*this)+N2FunctionPayload_.capacity()+source_->persistentBytes();
}
WVKernelStatus WVStratifiedModalRecord::matrixView(WVStratifiedScientificMatrix matrix, WVVerticalMatrixRecord& result) const {
    return source_->matrixView(matrix,result);
}
WVKernelStatus WVStratifiedModalRecord::prepareVertical(WVStratifiedModalOperator operation, WVComplexLayout input, WVComplexLayout output,
    std::unique_ptr<WVVerticalMatrixBackend> backend, std::unique_ptr<WVPreparedVerticalOperator>& result, WVAccumulation accumulation) const {
    return source_->prepareVertical(operation,std::move(input),std::move(output),
        std::move(backend),result,accumulation);
}
WVKernelStatus WVStratifiedModalRecord::horizontalSpecification(std::size_t planes, WVComplexRepresentation representation,
    const std::string& family, WVRetainedHorizontalSpecification& result) const {
    return source_->horizontalSpecification(planes,representation,family,result);
}
} // namespace wavevortex::runtime
