#pragma once
#include "WVStratifiedModalTestFixture.hpp"

namespace wavevortex::test_fixture {
inline void boussinesqFixture(const std::filesystem::path& path) {
    fixture(path); File f(path);
    for (const auto* name:{"WVTransform","AnnotatedClass"}) nc(nc_put_att_text(f.id,NC_GLOBAL,name,std::char_traits<char>::length("WVTransformBoussinesq"),"WVTransformBoussinesq"));
    int z,j,kl,group; nc(nc_inq_dimid(f.id,"z",&z)); nc(nc_inq_dimid(f.id,"j",&j)); nc(nc_inq_dimid(f.id,"kl",&kl)); nc(nc_def_dim(f.id,"K2unique",3,&group));
    const double pi=std::acos(-1.0);
    const auto vector=[&](const char* name,int dimension,const std::vector<double>& a) { int v; nc(nc_def_var(f.id,name,NC_DOUBLE,1,&dimension,&v)); nc(nc_put_var_double(f.id,v,a.data())); };
    vector("K2unique",group,{0,std::pow(2*pi/19000,2),std::pow(2*pi/11000,2)});
    vector("iK2unique",kl,{1,2,3});
    const auto matrix=[&](const char* name,const char* source,bool inverse) {
        int old; nc(nc_inq_varid(f.id,source,&old)); std::vector<double> data(21); nc(nc_get_var_double(f.id,old,data.data()));
        std::vector<double> all; for (int n=0;n<3;++n) all.insert(all.end(),data.begin(),data.end());
        int v; int dimensions[]={group,inverse ? j : z,inverse ? z : j}; nc(nc_def_var(f.id,name,NC_DOUBLE,3,dimensions,&v)); nc(nc_put_var_double(f.id,v,all.data()));
    };
    matrix("PFpmInv","PF0inv",true); matrix("QGpmInv","QG0inv",true); matrix("PFpm","PF0",false); matrix("QGpm","QG0",false);
    struct Scale { const char* name; const char* source; };
    const Scale scales[]={{"Ppm","P0"},{"Qpm","Q0"},{"h_pm","h_0"}};
    for (const auto& entry:scales) {
        int old; nc(nc_inq_varid(f.id,entry.source,&old)); std::vector<double> values(3); nc(nc_get_var_double(f.id,old,values.data()));
        std::vector<double> data; for (int n=0;n<3;++n) data.insert(data.end(),values.begin(),values.end());
        int v; int dimensions[]={std::string(entry.name)=="h_pm" ? kl : group,j}; nc(nc_def_var(f.id,entry.name,NC_DOUBLE,2,dimensions,&v)); nc(nc_put_var_double(f.id,v,data.data()));
    }
    std::vector<double> cross(27,0); for (int n=0;n<3;++n) { cross[4+9*n]=1; cross[8+9*n]=1; }
    int v; int dimensions[]={group,j,j}; nc(nc_def_var(f.id,"QGwg",NC_DOUBLE,3,dimensions,&v)); nc(nc_put_var_double(f.id,v,cross.data()));
}
inline void changeWave(const std::filesystem::path& path,const char* name,double value,std::size_t group,std::size_t column,std::size_t row) {
    File f(path); int v; nc(nc_inq_varid(f.id,name,&v)); const std::size_t index[]={group,column,row}; nc(nc_put_var1_double(f.id,v,index,&value));
}
} // namespace wavevortex::test_fixture
