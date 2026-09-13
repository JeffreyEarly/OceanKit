#include "WVOrderedRKCombination.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <type_traits>
#include <vector>
using wavevortex::WVComplex64;
using wavevortex::runtime::rk_detail::orderedWeightedAffine;

template<class T> __attribute__((noinline)) void setScaled(T* out,const T* in,std::size_t n,double a) {
    for (std::size_t i=0;i<n;++i) {
        if constexpr(std::is_same_v<T,double>) out[i]=a*in[i];
        else out[i]={a*in[i].real,a*in[i].imag};
    }
}
template<class T> __attribute__((noinline)) void addScaled(T* out,const T* in,std::size_t n,double a) {
    for (std::size_t i=0;i<n;++i) {
        if constexpr(std::is_same_v<T,double>) out[i]=out[i]+a*in[i];
        else out[i]={out[i].real+a*in[i].real,out[i].imag+a*in[i].imag};
    }
}
template<class T> __attribute__((noinline)) void affine(T* out,const T* base,std::size_t n,double h) {
    for (std::size_t i=0;i<n;++i) {
        if constexpr(std::is_same_v<T,double>) out[i]=base[i]+h*out[i];
        else out[i]={base[i].real+h*out[i].real,base[i].imag+h*out[i].imag};
    }
}
template<std::size_t N,class T> __attribute__((noinline)) void legacy(T* out,const T* base,const std::array<const T*,N>& inputs,const std::array<double,N>& weights,std::size_t n,double h) {
    setScaled(out,inputs[0],n,weights[0]);
    for(std::size_t j=1;j<N;++j) addScaled(out,inputs[j],n,weights[j]);
    affine(out,base,n,h);
}
template<std::size_t N,class T> __attribute__((noinline)) void fused(T* out,const T* base,const std::array<const T*,N>& inputs,const std::array<double,N>& weights,std::size_t n,double h) {
    orderedWeightedAffine(out,base,inputs,weights,n,h);
}
template<std::size_t N,class T> void screen(std::size_t n) {
    std::vector<T> storage((N+3)*n);auto* base=storage.data();auto* old=base+n;auto* candidate=old+n;
    std::array<const T*,N> inputs;std::array<double,N> weights;
    for(std::size_t i=0;i<storage.size();++i) {
        const double v=static_cast<double>(static_cast<int>(i%101)-50)/67.0;
        if constexpr(std::is_same_v<T,double>) storage[i]=v;
        else storage[i]={v,-.17*v+.013};
    }
    for(std::size_t j=0;j<N;++j) {inputs[j]=storage.data()+(j+3)*n;weights[j]=(j%2 ? -.19 : .23)*(j+1);}
    constexpr double h=.17;legacy(old,base,inputs,weights,n,h);fused(candidate,base,inputs,weights,n,h);
    if(std::memcmp(old,candidate,n*sizeof(T))) throw std::runtime_error("Legacy/fused bits differ");
    std::array<std::vector<double>,2> seconds;
    for(int sample=-2;sample<9;++sample) for(int turn=0;turn<2;++turn) {
        const int role=(sample+2+turn)%2;auto* out=role ? candidate : old;
        const auto start=std::chrono::steady_clock::now();
        if(role) fused(out,base,inputs,weights,n,h);else legacy(out,base,inputs,weights,n,h);
        asm volatile("" : : "g"(out) : "memory");
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(sample>=0) seconds[role].push_back(elapsed);
    }
    if(std::memcmp(old,candidate,n*sizeof(T))) throw std::runtime_error("Timed output bits differ");
    std::cout<<std::setprecision(17)<<"{\"kind\":\""<<(std::is_same_v<T,double>?"real":"complex")<<"\",\"terms\":"<<N<<",\"count\":"<<n<<",\"workspaceBytes\":"<<storage.size()*sizeof(T)<<",\"bitwiseEqual\":true,\"seconds\":[";
    for(std::size_t role=0;role<2;++role) {std::cout<<(role?",[":"[");for(std::size_t i=0;i<seconds[role].size();++i)std::cout<<(i?",":"")<<seconds[role][i];std::cout<<"]";}
    std::cout<<"]}\n";
}
int main(int argc,char** argv) try {
    const std::size_t n=argc>1?std::stoull(argv[1]):617706;
    screen<2,double>(n);screen<4,double>(n);screen<8,double>(n);screen<9,double>(n);
    screen<2,WVComplex64>(n);screen<4,WVComplex64>(n);screen<8,WVComplex64>(n);screen<9,WVComplex64>(n);
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
