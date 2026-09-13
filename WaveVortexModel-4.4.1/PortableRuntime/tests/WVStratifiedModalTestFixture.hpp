#pragma once
#include <netcdf.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace wavevortex::test_fixture {
inline void require(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
inline void nc(int code) { if (code!=NC_NOERR) throw std::runtime_error(nc_strerror(code)); }
struct File {
    int id=-1;
    explicit File(const std::filesystem::path& path,bool create=false) { if (create) nc(nc_create(path.c_str(),NC_NETCDF4|NC_CLOBBER,&id)); else nc(nc_open(path.c_str(),NC_WRITE,&id)); }
    ~File() { if(id>=0) nc_close(id); }
};
struct Temporary {
    std::filesystem::path path=std::filesystem::temp_directory_path()/("wvm-modal-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".nc");
    ~Temporary(){std::error_code ec;std::filesystem::remove(path,ec);}
};
inline void fixture(const std::filesystem::path& path,std::size_t Nz=7) {
    File f(path,true); int x,y,z,j,kl;
    nc(nc_def_dim(f.id,"x",8,&x)); nc(nc_def_dim(f.id,"y",6,&y)); nc(nc_def_dim(f.id,"z",Nz,&z)); nc(nc_def_dim(f.id,"j",3,&j)); nc(nc_def_dim(f.id,"kl",3,&kl));
    for(const auto* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformStratifiedQG"),"WVTransformStratifiedQG"));
    nc(nc_put_att_text(f.id,NC_GLOBAL,"model_version",5,"4.3.0"));
    auto scalar=[&](const char* name,double value) { int v; nc(nc_def_var(f.id,name,NC_DOUBLE,0,nullptr,&v)); nc(nc_put_var_double(f.id,v,&value)); };
    scalar("Lx",19000);scalar("Ly",11000);scalar("Lz",1000);scalar("g",9.81);scalar("rho0",1025);scalar("latitude",33);scalar("rotationRate",7.2921e-5);scalar("planetaryRadius",6.371e6);
    int mask; nc(nc_def_var(f.id,"shouldAntialias",NC_UBYTE,0,nullptr,&mask)); unsigned char yes=1; nc(nc_put_var_uchar(f.id,mask,&yes));
    auto vector=[&](const char* name,int dimension,const std::vector<double>& data) {int v;nc(nc_def_var(f.id,name,NC_DOUBLE,1,&dimension,&v));nc(nc_put_var_double(f.id,v,data.data()));};
    std::vector<double> xx(8),yy(6),zz(Nz),rho(Nz),weights(Nz),N2(Nz,1e-4),zero(Nz,0);
    for(std::size_t i=0;i<8;++i)xx[i]=19000.0*i/8;
    for(std::size_t i=0;i<6;++i)yy[i]=11000.0*i/6;
    for(std::size_t i=0;i<Nz;++i){zz[i]=-1000+1000.0*i/(Nz-1);rho[i]=1025-0.01*zz[i];weights[i]=1000.0/(Nz-1)*((i==0||i+1==Nz)?0.5:1);}
    vector("x",x,xx);vector("y",y,yy);vector("z",z,zz);vector("j",j,{0,1,2});vector("kl",kl,{0,1,2});
    const double pi=std::acos(-1.0);
    vector("k",kl,{0,2*pi/19000,0});vector("l",kl,{0,0,2*pi/11000});
    vector("N2",z,N2);vector("rho_nm0",z,rho);vector("dLnN2",z,zero);vector("z_int",z,weights);
    vector("P0",j,{1,2,3});vector("Q0",j,{1,4,5});vector("h_0",j,{1,0.5,0.2});
    std::vector<double> F(Nz*3),G(Nz*3),PF(Nz*3),QG(Nz*3);
    for(std::size_t c=0;c<3;++c){double maximum=0;for(std::size_t r=0;r<Nz;++r){F[r+Nz*c]=std::cos(c*pi*r/(Nz-1));G[r+Nz*c]=c==0||r==0||r+1==Nz?0:std::sin(c*pi*r/(Nz-1));maximum=std::max(maximum,std::abs(G[r+Nz*c]));PF[c+3*r]=F[r+Nz*c]/Nz;QG[c+3*r]=G[r+Nz*c]/Nz;}if(c)for(std::size_t r=0;r<Nz;++r)G[r+Nz*c]/=maximum;}
    auto matrix=[&](const char* name,int slow,int fast,const std::vector<double>& data){int v;const int dims[]={slow,fast};nc(nc_def_var(f.id,name,NC_DOUBLE,2,dims,&v));nc(nc_put_var_double(f.id,v,data.data()));};
    matrix("PF0inv",j,z,F);matrix("QG0inv",j,z,G);matrix("PF0",z,j,PF);matrix("QG0",z,j,QG);
}
inline void change(const std::filesystem::path& path,const char* name,double value,std::size_t first=0,std::size_t second=0) {
    File f(path);int v;nc(nc_inq_varid(f.id,name,&v));const std::size_t index[]={first,second};nc(nc_put_var1_double(f.id,v,index,&value));
}
} // namespace wavevortex::test_fixture
