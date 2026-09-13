#include "WVFieldEvaluationEventWorkspace.hpp"
#include "WVStratifiedFieldEvaluationAdapter.hpp"

#include "WaveVortexRuntime/WVIntegrationState.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <stdexcept>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {
WVPortableVariable portableField(WVHydrostaticField field) noexcept {
  switch(field) {
  case WVHydrostaticField::u:return WVPortableVariable::u;
  case WVHydrostaticField::v:return WVPortableVariable::v;
  case WVHydrostaticField::w:return WVPortableVariable::w;
  case WVHydrostaticField::eta:return WVPortableVariable::eta;
  case WVHydrostaticField::pi:return WVPortableVariable::pi;
  case WVHydrostaticField::p:return WVPortableVariable::p;
  case WVHydrostaticField::psi:return WVPortableVariable::psi;
  case WVHydrostaticField::qgpv:return WVPortableVariable::qgpv;
  case WVHydrostaticField::rhoE:return WVPortableVariable::rhoE;
  case WVHydrostaticField::rhoTotal:return WVPortableVariable::rhoTotal;
  case WVHydrostaticField::zetaX:return WVPortableVariable::zetaX;
  case WVHydrostaticField::zetaY:return WVPortableVariable::zetaY;
  case WVHydrostaticField::zetaZ:return WVPortableVariable::zetaZ;
  case WVHydrostaticField::ssh:return WVPortableVariable::ssh;
  case WVHydrostaticField::ssu:return WVPortableVariable::ssu;
  case WVHydrostaticField::ssv:return WVPortableVariable::ssv;
  }
  return WVPortableVariable::invalid;
}

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

std::size_t checkedProduct(std::size_t first, std::size_t second) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first)
    throw std::overflow_error("Stratified QG field size overflow");
  return first * second;
}

double wrapped(double coordinate, double length) noexcept {
  double result = std::fmod(coordinate, length);
  if (result < 0.0)
    result += length;
  return result >= length ? 0.0 : result;
}

std::uint64_t configurationFingerprint(
    const WVStratifiedModalGeometry &configuration) noexcept {
  std::uint64_t result = 1469598103934665603ULL;
  const auto append = [&](const auto &value) {
    const auto *bytes = reinterpret_cast<const std::uint8_t *>(&value);
    for (std::size_t index = 0; index < sizeof(value); ++index) {
      result ^= bytes[index];
      result *= 1099511628211ULL;
    }
  };
  for (char c:configuration.transformClass) append(c);
  append(configuration.Nz);
  append(configuration.Nj);
  append(configuration.Nkl);
  append(configuration.Lz);
  for (double z:configuration.z) append(z);
  append(configuration.Nx);
  append(configuration.Ny);
  append(configuration.Lx);
  append(configuration.Ly);
  for (double j:configuration.j) append(j);
  append(configuration.g);
  append(configuration.planetaryRadius);
  append(configuration.rotationRate);
  append(configuration.latitude);
  append(configuration.shouldAntialias);
  return result;
}

class SplineSystem final {
public:
  explicit SplineSystem(std::size_t count) : count_(count) {
    if (count < 2)
      throw std::invalid_argument("Spline interpolation requires two points.");
    if (count < 4)
      return;
    lu_.assign(checkedProduct(count, count), 0.0);
    pivots_.resize(count);
    const auto set = [&](std::size_t row, std::size_t column, double value) {
      lu_[column + count * row] = value;
    };
    set(0, 0, -1.0);
    set(0, 1, 2.0);
    set(0, 2, -1.0);
    for (std::size_t row = 1; row + 1 < count; ++row) {
      set(row, row - 1, 1.0);
      set(row, row, 4.0);
      set(row, row + 1, 1.0);
    }
    set(count - 1, count - 3, -1.0);
    set(count - 1, count - 2, 2.0);
    set(count - 1, count - 1, -1.0);
    std::vector<double> transpose(checkedProduct(count, count));
    for (std::size_t row = 0; row < count; ++row)
      for (std::size_t column = 0; column < count; ++column)
        transpose[column + count * row] = lu_[row + count * column];
    lu_.swap(transpose);
    factor();
  }

  void weightsInto(double spacing, double query, std::size_t circularShift,
                   std::vector<double> &result,
                   std::vector<double> &shifted,
                   std::vector<double> &rhs) const {
    const double normalized = query / spacing;
    std::fill(shifted.begin(), shifted.end(), 0.0);
    std::fill(rhs.begin(), rhs.end(), 0.0);
    std::fill(result.begin(), result.end(), 0.0);
    // MATLAB supplies extrapval=0 even after periodic boundary shifts.
    // On small grids a shifted query can still lie past the final knot.
    if (query < 0.0 || query > static_cast<double>(count_ - 1) * spacing)
      return;
    if (count_ == 2) {
      shifted[0] = 1.0 - normalized;
      shifted[1] = normalized;
    } else if (count_ == 3) {
      shifted[0] =
          (normalized - 1.0) * (normalized - 2.0) / 2.0;
      shifted[1] = -normalized * (normalized - 2.0);
      shifted[2] = normalized * (normalized - 1.0) / 2.0;
    } else {
      std::size_t interval = normalized <= 0.0
                                 ? 0
                                 : static_cast<std::size_t>(
                                       std::floor(normalized));
      interval = std::min(interval, count_ - 2);
      const double fraction = std::clamp(
          normalized - static_cast<double>(interval), 0.0, 1.0);
      const double first = 1.0 - fraction;
      const double second = fraction;
      rhs[interval] =
          (first * first * first - first) * spacing * spacing / 6.0;
      rhs[interval + 1] =
          (second * second * second - second) * spacing * spacing / 6.0;
      solve(rhs);
      shifted[interval] += first;
      shifted[interval + 1] += second;
      const double scale = 6.0 / (spacing * spacing);
      for (std::size_t row = 1; row + 1 < count_; ++row) {
        shifted[row - 1] += scale * rhs[row];
        shifted[row] -= 2.0 * scale * rhs[row];
        shifted[row + 1] += scale * rhs[row];
      }
    }
    for (std::size_t shiftedIndex = 0; shiftedIndex < count_; ++shiftedIndex) {
      const auto original =
          (shiftedIndex + count_ - circularShift % count_) % count_;
      result[original] += shifted[shiftedIndex];
    }
  }

  std::vector<double> weights(double spacing, double query,
                              std::size_t circularShift) const {
    std::vector<double> result(count_, 0.0);
    std::vector<double> shifted(count_, 0.0);
    std::vector<double> rhs(count_, 0.0);
    weightsInto(spacing, query, circularShift, result, shifted, rhs);
    return result;
  }

  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + lu_.capacity() * sizeof(double) +
           pivots_.capacity() * sizeof(std::size_t);
  }

private:
  void factor() {
    for (std::size_t column = 0; column < count_; ++column) {
      std::size_t pivot = column;
      double maximum = std::abs(lu_[column + count_ * column]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        const double value = std::abs(lu_[column + count_ * row]);
        if (value > maximum) {
          maximum = value;
          pivot = row;
        }
      }
      if (maximum == 0.0)
        throw std::invalid_argument("Spline interpolation matrix is singular.");
      pivots_[column] = pivot;
      if (pivot != column)
        for (std::size_t entry = 0; entry < count_; ++entry)
          std::swap(lu_[entry + count_ * column],
                    lu_[entry + count_ * pivot]);
      for (std::size_t row = column + 1; row < count_; ++row) {
        lu_[column + count_ * row] /= lu_[column + count_ * column];
        const double multiplier = lu_[column + count_ * row];
        for (std::size_t entry = column + 1; entry < count_; ++entry)
          lu_[entry + count_ * row] -=
              multiplier * lu_[entry + count_ * column];
      }
    }
  }

  void solve(std::vector<double> &rhs) const {
    for (std::size_t column = 0; column < count_; ++column) {
      if (pivots_[column] != column)
        std::swap(rhs[column], rhs[pivots_[column]]);
      for (std::size_t row = column + 1; row < count_; ++row)
        rhs[row] -= lu_[column + count_ * row] * rhs[column];
    }
    for (std::size_t reverse = 0; reverse < count_; ++reverse) {
      const auto row = count_ - 1 - reverse;
      for (std::size_t column = row + 1; column < count_; ++column)
        rhs[row] -= lu_[column + count_ * row] * rhs[column];
      rhs[row] /= lu_[row + count_ * row];
    }
  }

  std::size_t count_ = 0;
  std::vector<double> lu_;
  std::vector<std::size_t> pivots_;
};

class VerticalSpline final {
public:
  explicit VerticalSpline(const std::vector<double>& z):z_(z),second_(z.size()*z.size(),0) {
    const auto n=z_.size(); if(n<4) return;
    std::vector<double> matrix(n*n,0);
    auto set=[&](std::size_t r,std::size_t c,double v){matrix[r*n+c]=v;};
    auto h=[&](std::size_t i){return z_[i+1]-z_[i];};
    set(0,0,-h(1)); set(0,1,h(0)+h(1)); set(0,2,-h(0));
    for(std::size_t r=1;r+1<n;++r) {
      set(r,r-1,h(r-1)); set(r,r,2*(h(r-1)+h(r))); set(r,r+1,h(r));
      second_[r*n+r-1]=6/h(r-1); second_[r*n+r]=-6*(1/h(r-1)+1/h(r)); second_[r*n+r+1]=6/h(r);
    }
    set(n-1,n-3,-h(n-2)); set(n-1,n-2,h(n-3)+h(n-2)); set(n-1,n-1,-h(n-3));
    for(std::size_t col=0;col<n;++col) {
      std::size_t pivot=col;
      for(std::size_t r=col+1;r<n;++r) if(std::abs(matrix[r*n+col])>std::abs(matrix[pivot*n+col])) pivot=r;
      if(matrix[pivot*n+col]==0) throw std::invalid_argument("Singular vertical interpolation grid.");
      for(std::size_t c=0;c<n;++c) {std::swap(matrix[col*n+c],matrix[pivot*n+c]);std::swap(second_[col*n+c],second_[pivot*n+c]);}
      const double diagonal=matrix[col*n+col];
      for(std::size_t c=0;c<n;++c) {matrix[col*n+c]/=diagonal;second_[col*n+c]/=diagonal;}
      for(std::size_t r=0;r<n;++r) if(r!=col) {
        const double factor=matrix[r*n+col];
        for(std::size_t c=0;c<n;++c) {matrix[r*n+c]-=factor*matrix[col*n+c];second_[r*n+c]-=factor*second_[col*n+c];}
      }
    }
  }
  void weightsInto(double z,std::vector<double>& weights) const {
    const auto n=z_.size(); std::fill(weights.begin(),weights.end(),0);
    if(z<z_.front() || z>z_.back()) return;
    if(n<4) {
      for(std::size_t i=0;i<n;++i) {weights[i]=1;for(std::size_t j=0;j<n;++j) if(i!=j) weights[i]*=(z-z_[j])/(z_[i]-z_[j]);}
      return;
    }
    const auto i=std::min(static_cast<std::size_t>(std::upper_bound(z_.begin(),z_.end(),z)-z_.begin()-1),n-2);
    const double h=z_[i+1]-z_[i],b=(z-z_[i])/h,a=1-b;
    for(std::size_t j=0;j<n;++j) weights[j]=(h*h/6)*((a*a*a-a)*second_[i*n+j]+(b*b*b-b)*second_[(i+1)*n+j]);
    weights[i]+=a; weights[i+1]+=b;
  }
  std::size_t persistentBytes() const noexcept {return sizeof(*this)+(z_.capacity()+second_.capacity())*sizeof(double);}
private:
  std::vector<double> z_,second_;
};

bool surface(WVHydrostaticField field) {return field==WVHydrostaticField::ssh || field==WVHydrostaticField::ssu || field==WVHydrostaticField::ssv;}

enum class ScalarField : std::uint8_t { none, energy, uvMax, wMax };

struct Weight {
  std::array<std::size_t, 2> xIndices{};
  std::array<std::size_t, 2> yIndices{};
  std::array<double, 2> xWeights{};
  std::array<double, 2> yWeights{};
  std::vector<double> xSpline;
  std::vector<double> ySpline;
  std::vector<double> zSpline;
  std::array<std::size_t,2> zIndices{};
  std::array<double,2> zWeights{1,0};
  bool outside = false;
};

struct Request {
  WVHydrostaticField field = WVHydrostaticField::u;
  ScalarField scalar = ScalarField::none;
  WVFieldSamplingKind sampling = WVFieldSamplingKind::fullGrid;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::vector<Weight> weights;
  std::vector<std::size_t> profileIndices;
  std::size_t output = 0;
};

struct Plan {
  std::uint64_t fingerprint = 0;
  std::vector<Request> requests;
  std::size_t persistentBytes() const noexcept {
    std::size_t bytes = sizeof(*this) + requests.capacity() * sizeof(Request);
    for (const auto &request : requests) {
      bytes += request.weights.capacity() * sizeof(Weight)+request.profileIndices.capacity()*sizeof(std::size_t);
      for (const auto &weight : request.weights)
        bytes += (weight.xSpline.capacity() + weight.ySpline.capacity() + weight.zSpline.capacity()) *
                 sizeof(double);
    }
    return bytes;
  }
};

struct MovingRequest {
  WVHydrostaticField field = WVHydrostaticField::u;
  std::size_t offset = 0;
  std::size_t count = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::size_t output = 0;
};

struct MovingPlan {
  std::uint64_t fingerprint = 0;
  std::vector<MovingRequest> requests;
  std::size_t positionCount = 0;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + requests.capacity() * sizeof(MovingRequest);
  }
};

struct EventRequest {
  WVHydrostaticField field = WVHydrostaticField::u;
  std::size_t positionSet = 0;
  WVPositionInterpolation interpolation = WVPositionInterpolation::linear;
  std::size_t output = 0;
};

struct EventPlan {
  std::uint64_t fingerprint = 0;
  std::vector<EventRequest> requests;
  std::size_t positionSetCount = 0;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + requests.capacity() * sizeof(EventRequest);
  }
};

struct EventGeometry {
  std::shared_ptr<const EventPlan> eventPlan;
  std::shared_ptr<const Plan> plan;
  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) + (plan ? plan->persistentBytes() : 0);
  }
};

WVKernelStatus resolveField(const std::string &name,
                            WVHydrostaticField &field,
                            ScalarField &scalar, bool hydrostatic) {
  scalar = ScalarField::none;
  if (name == "u")
    field = WVHydrostaticField::u;
  else if (name == "v")
    field = WVHydrostaticField::v;
  else if (name == "w" && hydrostatic) field=WVHydrostaticField::w;
  else if (name == "wMax" && hydrostatic) scalar=ScalarField::wMax;
  else if (name == "zeta_x" && hydrostatic) field=WVHydrostaticField::zetaX;
  else if (name == "zeta_y" && hydrostatic) field=WVHydrostaticField::zetaY;
  else if (name == "p") field = WVHydrostaticField::p;
  else if (name == "rho_e") field = WVHydrostaticField::rhoE;
  else if (name == "rho_total") field = WVHydrostaticField::rhoTotal;
  else if (name == "ssu") field = WVHydrostaticField::ssu;
  else if (name == "ssv") field = WVHydrostaticField::ssv;
  else if (name == "eta")
    field = WVHydrostaticField::eta;
  else if (name == "pi")
    field = WVHydrostaticField::pi;
  else if (name == "psi")
    field = WVHydrostaticField::psi;
  else if (name == "qgpv")
    field = WVHydrostaticField::qgpv;
  else if (name == "zeta_z")
    field = WVHydrostaticField::zetaZ;
  else if (name == "ssh")
    field = WVHydrostaticField::ssh;
  else if (name == "energy")
    scalar = ScalarField::energy;
  else if (name == "uvMax")
    scalar = ScalarField::uvMax;
  else
    return {WVKernelStatusCode::unsupportedOperation,
            "This stratified transform does not support field " + name + "."};
  return WVKernelStatus::ok();
}

Weight weightsFor(double x, double y, double z, bool horizontal, WVPositionInterpolation interpolation,
                  const WVStratifiedModalGeometry &configuration,
                  const SplineSystem *xSpline,
                  const SplineSystem *ySpline, const VerticalSpline *zSpline) {
  const double dx = configuration.Lx / static_cast<double>(configuration.Nx);
  const double dy = configuration.Ly / static_cast<double>(configuration.Ny);
  const double xWrapped = wrapped(x, configuration.Lx);
  const double yWrapped = wrapped(y, configuration.Ly);
  const auto xLower = std::min(
      static_cast<std::size_t>(std::floor(xWrapped / dx)),
      configuration.Nx - 1);
  const auto yLower = std::min(
      static_cast<std::size_t>(std::floor(yWrapped / dy)),
      configuration.Ny - 1);
  Weight result;
  if (interpolation == WVPositionInterpolation::linear) {
    result.xIndices = {xLower, (xLower + 1) % configuration.Nx};
    result.yIndices = {yLower, (yLower + 1) % configuration.Ny};
    const double xFraction =
        (xWrapped - static_cast<double>(xLower) * dx) / dx;
    const double yFraction =
        (yWrapped - static_cast<double>(yLower) * dy) / dy;
    result.xWeights = {1.0 - xFraction, xFraction};
    result.yWeights = {1.0 - yFraction, yFraction};
  } else {
    const bool xBoundary = configuration.Nx < 4 || xLower < 3 ||
                           xLower > configuration.Nx - 4;
    const bool yBoundary = configuration.Ny < 4 || yLower < 3 ||
                           yLower > configuration.Ny - 4;
    const std::size_t xShift = xBoundary ? 4 : 0;
    const std::size_t yShift = yBoundary ? 4 : 0;
    const double xQuery = xBoundary
                              ? wrapped(x + 4.0 * dx, configuration.Lx)
                              : xWrapped;
    const double yQuery = yBoundary
                              ? wrapped(y + 4.0 * dy, configuration.Ly)
                              : yWrapped;
    result.xSpline = xSpline->weights(dx, xQuery, xShift);
    result.ySpline = ySpline->weights(dy, yQuery, yShift);
  }
  if (!horizontal) {
    result.outside=z<configuration.z.front() || z>configuration.z.back();
    if (!result.outside) {
      auto iz=std::min(static_cast<std::size_t>(std::upper_bound(configuration.z.begin(),configuration.z.end(),z)-configuration.z.begin()-1),configuration.Nz-2);
      const double fraction=(z-configuration.z[iz])/(configuration.z[iz+1]-configuration.z[iz]);
      result.zIndices={iz,iz+1}; result.zWeights={1-fraction,fraction};
      if(interpolation==WVPositionInterpolation::spline) {result.zSpline.resize(configuration.Nz);zSpline->weightsInto(z,result.zSpline);}
    }
  } else if(interpolation==WVPositionInterpolation::spline) result.zSpline={1};
  return result;
}

double interpolate(const double *field,const Weight& weight,WVPositionInterpolation interpolation,std::size_t Nx,std::size_t Ny,std::size_t Nz) noexcept {
  if(weight.outside) return 0;
  double value=0;
  if(interpolation==WVPositionInterpolation::linear) {
    for(std::size_t z=0;z<(Nz==1?1:2);++z) for(std::size_t y=0;y<2;++y) for(std::size_t x=0;x<2;++x)
      value+=weight.xWeights[x]*weight.yWeights[y]*weight.zWeights[z]*field[weight.xIndices[x]+Nx*(weight.yIndices[y]+Ny*weight.zIndices[z])];
  } else {
    for(std::size_t z=0;z<Nz;++z) for(std::size_t y=0;y<Ny;++y) for(std::size_t x=0;x<Nx;++x)
      value+=weight.xSpline[x]*weight.ySpline[y]*weight.zSpline[z]*field[x+Nx*(y+Ny*z)];
  }
  return value;
}

WVKernelStatus coefficientView(const WVIntegrationState& state,const WVStratifiedModalGeometry& g,WVState& result) {
  const bool hydro=(g.transformClass=="WVTransformHydrostatic" || g.transformClass=="WVTransformBoussinesq");
  const auto count=hydro ? 3U : 1U;
  if (state.coefficientFamilyCount!=count || !state.coefficientFamilies) return invalid("Wrong stratified coefficient family count.");
  const char* names[]={"Ap","Am","A0"};
  WVComplexConstView views[3]{};
  for (std::size_t i=0;i<count;++i) {
    const auto& f=state.coefficientFamilies[i]; const auto family=hydro ? i : 2;
    if (!f.data || !f.layout || f.layout->identifier!=names[family] || (f.layout->spectralDimensions.size()!=2 || f.layout->spectralDimensions[0]!=g.Nj || f.layout->spectralDimensions[1]!=g.Nkl)) return invalid("Invalid stratified field coefficient storage.");
    views[family]={f.data,{g.Nj,g.Nkl}};
  }
  result={state.waveVortex.t,state.waveVortex.t0,{views[0],views[1],views[2]}};
  return WVKernelStatus::ok();
}

} // namespace

struct WVStratifiedFieldEvaluationAdapter::MovingInterpolationWorkspace {
  MovingInterpolationWorkspace(std::size_t Nx, std::size_t Ny,const std::vector<double>& z)
      : xSpline(Nx), ySpline(Ny), zSpline(z), xWeights(Nx), yWeights(Ny), zWeights(z.size()),
        xShifted(Nx), yShifted(Ny), xRhs(Nx), yRhs(Ny) {}

  std::size_t persistentBytes() const noexcept {
    return sizeof(*this) +
           xSpline.persistentBytes() + ySpline.persistentBytes() + zSpline.persistentBytes() +
           (xWeights.capacity() + yWeights.capacity() + zWeights.capacity() +
            xShifted.capacity() + yShifted.capacity() + xRhs.capacity() +
            yRhs.capacity()) *
               sizeof(double);
  }

  std::size_t scratchBytes() const noexcept {
    return (xWeights.capacity() + yWeights.capacity() + zWeights.capacity() +
            xShifted.capacity() + yShifted.capacity() + xRhs.capacity() +
            yRhs.capacity()) *
           sizeof(double);
  }

  SplineSystem xSpline;
  SplineSystem ySpline;
  VerticalSpline zSpline;
  std::vector<double> xWeights;
  std::vector<double> yWeights;
  std::vector<double> zWeights;
  std::vector<double> xShifted;
  std::vector<double> yShifted;
  std::vector<double> xRhs;
  std::vector<double> yRhs;
};

WVStratifiedFieldEvaluationAdapter::~WVStratifiedFieldEvaluationAdapter() =
    default;

WVKernelStatus WVStratifiedFieldEvaluationAdapter::beginStateEvaluation(
    const WVIntegrationState& state,const void* owner) {
  WVState coefficients;
  const auto status=coefficientView(state,configuration(),coefficients);
  if(!status) return status;
  if(hydrostaticKernel_) return hydrostaticKernel_->beginStateEvaluation(coefficients,owner);
  if(boussinesqKernel_) return boussinesqKernel_->beginStateEvaluation(coefficients,owner);
  return kernel_->beginStateEvaluation(coefficients.coefficients.A0,owner);
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::addStateEvaluationView(
    const WVIntegrationState& state,const void* owner,
    std::size_t componentIdentity) {
  WVState coefficients;
  const auto status=coefficientView(state,configuration(),coefficients);
  if(!status) return status;
  if(hydrostaticKernel_) return hydrostaticKernel_->addStateEvaluationView(
      coefficients,owner,componentIdentity);
  if(boussinesqKernel_) return boussinesqKernel_->addStateEvaluationView(
      coefficients,owner,componentIdentity);
  return kernel_->addStateEvaluationView(
      coefficients.coefficients.A0,owner,componentIdentity);
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::removeStateEvaluationView(
    const WVIntegrationState& state,const void* owner,
    std::size_t componentIdentity) {
  WVState coefficients;
  const auto status=coefficientView(state,configuration(),coefficients);
  if(!status) return status;
  if(hydrostaticKernel_) return hydrostaticKernel_->removeStateEvaluationView(
      coefficients,owner,componentIdentity);
  if(boussinesqKernel_) return boussinesqKernel_->removeStateEvaluationView(
      coefficients,owner,componentIdentity);
  return kernel_->removeStateEvaluationView(
      coefficients.coefficients.A0,owner,componentIdentity);
}

void WVStratifiedFieldEvaluationAdapter::endStateEvaluation() noexcept {
  if(hydrostaticKernel_) (void)hydrostaticKernel_->endStateEvaluation();
  else if(boussinesqKernel_) (void)boussinesqKernel_->endStateEvaluation();
  else if(kernel_) (void)kernel_->endStateEvaluation();
}

WVVariableProducerMetrics
WVStratifiedFieldEvaluationAdapter::producerMetrics() const noexcept {
  WVVariableProducerMetrics result;
  if(hydrostaticKernel_) {
    const auto& metrics=hydrostaticKernel_->metrics();
    result.stateValidations=metrics.stateValidationCount;
    result.phasePreparations=metrics.phasePreparationCount;
    result.derivedValidations=metrics.derivedValidationCount;
    result.coefficientAssemblies=metrics.coefficientAssemblyCount;
    result.verticalOperatorExecutions=metrics.verticalOperatorExecutionCount;
    result.verticalPreparations=metrics.verticalPreparationCount;
    result.horizontalSpectrumReuses=metrics.horizontalSpectrumReuseCount;
    result.derivativeAdvectionConsumers=metrics.derivativeAdvectionConsumerCount;
    result.tiledNonlinearExecutions=metrics.tiledNonlinearCount;
    result.tiledColumnInverses=metrics.tiledColumnInverseCount;
    result.tiledRowInverses=metrics.tiledRowInverseCount;
    result.tiledReusedColumns=metrics.tiledReusedColumnCount;

    result.preparedVerticalDerivatives=metrics.preparedVerticalDerivativeCount;
    result.tendencyReconstructions=metrics.tendencyReconstructionCount;
    result.reconstructions=metrics.reconstructionCount;
  } else if(boussinesqKernel_) {
    const auto& metrics=boussinesqKernel_->metrics();
    result.stateValidations=metrics.stateValidationCount;
    result.phasePreparations=metrics.phasePreparationCount;
    result.derivedValidations=metrics.derivedValidationCount;
    result.coefficientAssemblies=metrics.coefficientAssemblyCount;
    result.verticalOperatorExecutions=metrics.verticalOperatorExecutionCount;
    result.verticalMatrixGroupExecutions=metrics.verticalMatrixGroupExecutionCount;
    result.verticalPreparations=metrics.verticalPreparationCount;
    result.horizontalSpectrumReuses=metrics.horizontalSpectrumReuseCount;
    result.derivativeAdvectionConsumers=metrics.derivativeAdvectionConsumerCount;
    result.tiledNonlinearExecutions=metrics.tiledNonlinearCount;
    result.tiledColumnInverses=metrics.tiledColumnInverseCount;
    result.tiledRowInverses=metrics.tiledRowInverseCount;
    result.tiledReusedColumns=metrics.tiledReusedColumnCount;

    result.preparedVerticalDerivatives=metrics.preparedVerticalDerivativeCount;
    result.tendencyReconstructions=metrics.tendencyReconstructionCount;
    result.reconstructions=metrics.reconstructionCount;
  } else {
    const auto& metrics=kernel_->metrics();
    result.stateValidations=metrics.stateValidationCount;
    result.horizontalSpeedReductions=
        metrics.horizontalSpeedMaximumReductionCount;
    for(std::size_t field=0;field<metrics.componentReconstructionCount.size();++field)
      for(std::size_t derivative=0;
          derivative<metrics.componentReconstructionCount[field].size();++derivative)
        result.reconstructions[field][derivative]=
            metrics.componentReconstructionCount[field][derivative];
  }
  result.horizontalSpeedReductions+=
      outputProducerMetrics_.horizontalSpeedReductions;
  result.verticalSpeedReductions+=
      outputProducerMetrics_.verticalSpeedReductions;
  result.energyReductions+=outputProducerMetrics_.energyReductions;
  return result;
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::create(
    std::shared_ptr<const WVStratifiedModalSource> source,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVStratifiedFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVStratifiedFieldEvaluationAdapter>(
        new WVStratifiedFieldEvaluationAdapter());
    if (!source) return invalid("A stratified scientific source is required.");
    WVKernelStatus status;
    if (source->geometry().transformClass=="WVTransformHydrostatic") {
      status=WVTransformHydrostaticKernel::create(source,std::move(engine),candidate->ownedHydrostatic_);
      candidate->hydrostaticKernel_=candidate->ownedHydrostatic_.get();
    } else if (source->geometry().transformClass=="WVTransformBoussinesq") {
      status=WVTransformBoussinesqKernel::create(source,std::move(engine),candidate->ownedBoussinesq_);
      candidate->boussinesqKernel_=candidate->ownedBoussinesq_.get();
    } else {
      status=WVTransformStratifiedQGKernel::create(source,std::move(engine),candidate->ownedKernel_);
      candidate->kernel_=candidate->ownedKernel_.get();
    }
    if (!status) return status;
    const auto& g=source->geometry();
    candidate->fieldScratch_.resize(g.Nx*g.Ny*g.Nz);
    if (candidate->hydrostaticKernel_ || candidate->boussinesqKernel_) candidate->speedScratch_.resize(candidate->fieldScratch_.size());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            source->geometry().Nx, source->geometry().Ny, source->geometry().z);
    candidate->metrics_.transformPersistentBytes =
        candidate->boussinesqKernel_ ? candidate->boussinesqKernel_->persistentBytes() : candidate->hydrostaticKernel_ ? candidate->hydrostaticKernel_->persistentBytes() : candidate->kernel_->persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        (candidate->fieldScratch_.capacity()+candidate->speedScratch_.capacity()) * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate Stratified QG field service."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createBorrowing(
    WVTransformStratifiedQGKernel &kernel,
    std::unique_ptr<WVStratifiedFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVStratifiedFieldEvaluationAdapter>(
        new WVStratifiedFieldEvaluationAdapter());
    candidate->kernel_ = &kernel;
    candidate->fieldScratch_.resize(
        kernel.spatialShape().elementCount());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            kernel.geometry().Nx, kernel.geometry().Ny, kernel.geometry().z);
    candidate->metrics_.transformPersistentBytes = kernel.persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        (candidate->fieldScratch_.capacity()+candidate->speedScratch_.capacity()) * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed Stratified QG field service."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createBorrowing(
    WVTransformHydrostaticKernel &kernel,
    std::unique_ptr<WVStratifiedFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVStratifiedFieldEvaluationAdapter>(
        new WVStratifiedFieldEvaluationAdapter());
    candidate->hydrostaticKernel_ = &kernel;
    candidate->speedScratch_.resize(kernel.spatialShape().elementCount());
    candidate->fieldScratch_.resize(
        kernel.spatialShape().elementCount());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            kernel.geometry().Nx, kernel.geometry().Ny, kernel.geometry().z);
    candidate->metrics_.transformPersistentBytes = kernel.persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        (candidate->fieldScratch_.capacity()+candidate->speedScratch_.capacity()) * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed Stratified QG field service."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createBorrowing(
    WVTransformBoussinesqKernel &kernel,
    std::unique_ptr<WVStratifiedFieldEvaluationAdapter> &adapter) {
  adapter.reset();
  try {
    auto candidate = std::unique_ptr<WVStratifiedFieldEvaluationAdapter>(
        new WVStratifiedFieldEvaluationAdapter());
    candidate->boussinesqKernel_ = &kernel;
    candidate->speedScratch_.resize(kernel.spatialShape().elementCount());
    candidate->fieldScratch_.resize(
        kernel.spatialShape().elementCount());
    candidate->movingInterpolation_ =
        std::make_unique<MovingInterpolationWorkspace>(
            kernel.geometry().Nx, kernel.geometry().Ny, kernel.geometry().z);
    candidate->metrics_.transformPersistentBytes = kernel.persistentBytes();
    candidate->metrics_.scratchCapacityBytes =
        (candidate->fieldScratch_.capacity()+candidate->speedScratch_.capacity()) * sizeof(double) +
        candidate->movingInterpolation_->scratchBytes();
    candidate->metrics_.servicePersistentBytes = candidate->persistentBytes();
    adapter = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate borrowed Boussinesq field service."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createPlan(
    const std::vector<WVFieldRequest> &requests,
    WVFieldEvaluationPlan &plan) const {
  try {
    auto implementation = std::make_shared<Plan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    implementation->requests.reserve(requests.size());
    WVFieldEvaluationPlan candidate;
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    const bool needsSpline = std::any_of(
        requests.begin(), requests.end(), [](const auto &request) {
          return request.sampling.kind == WVFieldSamplingKind::positions &&
                 request.sampling.interpolation ==
                     WVPositionInterpolation::spline;
        });
    std::unique_ptr<SplineSystem> xSpline;
    std::unique_ptr<SplineSystem> ySpline;
    std::unique_ptr<VerticalSpline> zSpline;
    if (needsSpline) {
      xSpline = std::make_unique<SplineSystem>(configuration().Nx);
      ySpline = std::make_unique<SplineSystem>(configuration().Ny);
      zSpline = std::make_unique<VerticalSpline>(configuration().z);
    }
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &input = requests[index];
      if (input.identifier.empty() || !identifiers.insert(input.identifier).second)
        return invalid("Stratified QG field identifiers must be nonempty and unique.");
      Request request;
      auto status = resolveField(input.fieldName, request.field, request.scalar,(hydrostaticKernel_!=nullptr || boussinesqKernel_!=nullptr));
      if (!status)
        return status;
      request.sampling = input.sampling.kind;
      request.interpolation = input.sampling.interpolation;
      request.output = index;
      WVFieldOutputSpecification output{
          input.identifier, input.fieldName, input.sampling.kind, {}, 0};
      if (request.scalar != ScalarField::none) {
        if (input.sampling.kind != WVFieldSamplingKind::fullGrid)
          return invalid("Scalar Stratified QG fields do not support position sampling.");
        output.elementCount = 1;
      } else if (input.sampling.kind == WVFieldSamplingKind::fullGrid) {
        output.dimensions = {configuration().Nx, configuration().Ny};
        output.elementCount = checkedProduct(configuration().Nx,configuration().Ny);
        if (!surface(request.field)) {output.dimensions.push_back(configuration().Nz);output.elementCount=checkedProduct(output.elementCount,configuration().Nz);}
      } else if (input.sampling.kind == WVFieldSamplingKind::positions) {
        if (input.sampling.x.empty() ||
            input.sampling.x.size() != input.sampling.y.size() || (!surface(request.field) && input.sampling.z.size()!=input.sampling.x.size()))
          return invalid("Stratified QG position arrays must be equal and nonempty.");
        if (input.sampling.interpolation != WVPositionInterpolation::linear &&
            input.sampling.interpolation != WVPositionInterpolation::spline)
          return invalid("Stratified QG interpolation method is unsupported.");
        request.weights.reserve(input.sampling.x.size());
        for (std::size_t position = 0; position < input.sampling.x.size();
             ++position) {
          if (!std::isfinite(input.sampling.x[position]) ||
              !std::isfinite(input.sampling.y[position]) || (!surface(request.field) && !std::isfinite(input.sampling.z[position])))
            return invalid("Stratified QG positions must be finite.");
          request.weights.push_back(weightsFor(
              input.sampling.x[position], input.sampling.y[position],surface(request.field)?0:input.sampling.z[position],surface(request.field),
              input.sampling.interpolation, configuration(), xSpline.get(),ySpline.get(),zSpline.get()));
        }
        output.dimensions = {input.sampling.x.size()};
        output.elementCount = input.sampling.x.size();
      } else if(input.sampling.kind==WVFieldSamplingKind::fixedVerticalProfiles && !surface(request.field)) {
        if(input.sampling.xIndices.empty() || input.sampling.xIndices.size()!=input.sampling.yIndices.size()) return invalid("Profile grid indices must be paired and nonempty.");
        for(std::size_t i=0;i<input.sampling.xIndices.size();++i) {
          const auto x=input.sampling.xIndices[i],y=input.sampling.yIndices[i];
          if(!x || x>configuration().Nx || !y || y>configuration().Ny) return invalid("Profile index outside the grid.");
          request.profileIndices.push_back(x-1+configuration().Nx*(y-1));
        }
        output.dimensions={configuration().Nz,request.profileIndices.size()};output.elementCount=checkedProduct(configuration().Nz,request.profileIndices.size());
      } else return invalid("Unsupported SQG sampling request.");
      implementation->requests.push_back(std::move(request));
      candidate.outputs_.push_back(std::move(output));
    }
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Stratified QG field plan."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::evaluate(
    const WVFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const auto plan =
      std::static_pointer_cast<const Plan>(publicPlan.transformPlan_);
  if (!plan || plan->fingerprint != configurationFingerprint(configuration()))
    return invalid("Stratified QG field plan belongs to another transform.");
  if (outputCount != publicPlan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Stratified QG outputs do not match the resolved plan."};
  for (std::size_t output = 0; output < outputCount; ++output)
    if ((!activeOutputs || activeOutputs[output]) && (outputs[output].data == nullptr ||
        outputs[output].elementCount !=
            publicPlan.outputs_[output].elementCount))
      return {WVKernelStatusCode::invalidShape,
              "A Stratified QG output has the wrong shape."};
  WVState amplitudes;
  auto status = coefficientView(state, configuration(), amplitudes);
  if (!status)
    return status;
  if (executing_)
    return {WVKernelStatusCode::reentrantExecution,
            "Stratified QG field evaluation is not reentrant."};
  executing_ = true;
  struct Guard {
    bool &executing;
    ~Guard() { executing = false; }
  } guard{executing_};

  std::array<bool, 16> evaluated{};
  const WVShape3D spatial{configuration().Nx,configuration().Ny,configuration().Nz};
  for (const auto &request : plan->requests) {
    if (activeOutputs && !activeOutputs[request.output]) continue;
    if (request.scalar != ScalarField::none) {
      double value = 0.0;
      const auto variable=request.scalar==ScalarField::energy ?
          WVPortableVariable::energy : request.scalar==ScalarField::uvMax ?
          WVPortableVariable::uvMax : WVPortableVariable::wMax;
      const auto produce=[&]() {
        const auto produced=scalarValue(
            amplitudes,static_cast<unsigned>(request.scalar),value);
        if(produced) {
          if(request.scalar==ScalarField::energy)
            ++outputProducerMetrics_.energyReductions;
          else if(request.scalar==ScalarField::uvMax &&
              (hydrostaticKernel_ || boussinesqKernel_))
            ++outputProducerMetrics_.horizontalSpeedReductions;
          else ++outputProducerMetrics_.verticalSpeedReductions;
        }
        return produced;
      };
      bool reused=false;
      const auto component=request.scalar==ScalarField::energy && eventWorkspace_ ?
          eventWorkspace_->component() : 0u;
      status=eventWorkspace_ ? eventWorkspace_->evaluate(
          {WVVariableEvaluationNode::reduction,
            static_cast<std::uint32_t>(variable),component},&value,1,produce,reused) :
          produce();
      if (!status)
        return status;
      outputs[request.output].data[0] = value;
      ++metrics_.outputElementWriteCount;
      continue;
    }
    const auto fieldIndex = static_cast<std::size_t>(request.field);
    if (!evaluated[fieldIndex]) {
      WVRealVolumeView view{fieldScratch_.data(),{spatial.first,spatial.second,surface(request.field)?1:spatial.third}};
      bool reused=false;
      status = transformField(amplitudes, request.field, view,&reused);
      if (!status)
        return status;

      if(reused) ++metrics_.primitiveFieldReuseCount;
      else {++metrics_.transformCount; ++metrics_.primitiveFieldEvaluationCount;}
      evaluated[fieldIndex] = true;
      for (const auto &destination : plan->requests) {
        if (activeOutputs && !activeOutputs[destination.output]) continue;
        if (destination.scalar != ScalarField::none ||
            destination.field != request.field)
          continue;
        auto &output = outputs[destination.output];
        if (destination.sampling == WVFieldSamplingKind::fullGrid) {
          std::copy_n(fieldScratch_.data(),output.elementCount,output.data);
          ++metrics_.fullGridWriteCount;
        } else if(destination.sampling==WVFieldSamplingKind::fixedVerticalProfiles) {
          const auto plane=configuration().Nx*configuration().Ny;
          for(std::size_t p=0;p<destination.profileIndices.size();++p) for(std::size_t z=0;z<configuration().Nz;++z)
            output.data[z+configuration().Nz*p]=fieldScratch_[destination.profileIndices[p]+plane*z];
        } else {
          for (std::size_t position = 0;
               position < destination.weights.size(); ++position)
            output.data[position] = interpolate(
                fieldScratch_.data(), destination.weights[position],
                destination.interpolation, configuration().Nx,
                configuration().Ny,surface(destination.field)?1:configuration().Nz);
          if (destination.interpolation == WVPositionInterpolation::linear)
            metrics_.linearInterpolationCount += destination.weights.size();
          else
            metrics_.splineInterpolationCount += destination.weights.size();
        }
        metrics_.outputElementWriteCount += output.elementCount;
        if (&destination != &request)
          ++metrics_.primitiveFieldReuseCount;
      }
    }
  }
  ++metrics_.evaluationCount;
  metrics_.scratchHighWaterBytes = std::max(
      metrics_.scratchHighWaterBytes,
      fieldScratch_.size() * sizeof(double));
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createMovingPlan(
    const std::vector<WVMovingFieldRequest> &requests,
    WVMovingFieldEvaluationPlan &plan) const {
  try {
    auto implementation = std::make_shared<MovingPlan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    WVMovingFieldEvaluationPlan candidate;
    std::set<std::string> identifiers;
    for (std::size_t index = 0; index < requests.size(); ++index) {
      const auto &input = requests[index];
      WVHydrostaticField field;
      ScalarField scalar;
      auto status = resolveField(input.fieldName, field, scalar,(hydrostaticKernel_!=nullptr || boussinesqKernel_!=nullptr));
      if (!status || scalar != ScalarField::none)
        return status ? invalid("A scalar field cannot be sampled at moving positions.") : status;
      if (input.identifier.empty() || !identifiers.insert(input.identifier).second ||
          input.positionCount == 0 ||
          (input.interpolation != WVPositionInterpolation::linear &&
           input.interpolation != WVPositionInterpolation::spline) ||
          input.positionOffset > std::numeric_limits<std::size_t>::max() -
                                     input.positionCount)
        return invalid("Stratified QG moving-field request is invalid.");
      implementation->positionCount =
          std::max(implementation->positionCount,
                   input.positionOffset + input.positionCount);
      implementation->requests.push_back(
          {field, input.positionOffset, input.positionCount,
           input.interpolation, index});
      candidate.outputs_.push_back(
          {input.identifier, input.fieldName, WVFieldSamplingKind::positions,
           {input.positionCount}, input.positionCount});
    }
    candidate.positionCount_ = implementation->positionCount;
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Stratified QG moving-field plan."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::samplePreparedField(
    const WVFieldEvaluationPlan &publicPlan, const double *source,
    WVFieldOutputView output) {
  const auto plan = std::static_pointer_cast<const Plan>(publicPlan.transformPlan_);
  if (!plan || plan->fingerprint != configurationFingerprint(configuration()) ||
      plan->requests.size() != 1 || publicPlan.outputs_.size() != 1)
    return invalid("Prepared stratified sampler does not match this transform.");
  const auto &request = plan->requests.front();
  if (!source || !output.data ||
      output.elementCount != publicPlan.outputs_.front().elementCount)
    return {WVKernelStatusCode::invalidShape,
            "Prepared stratified sampler storage has the wrong shape."};
  const auto &g = configuration();
  const auto plane = g.Nx * g.Ny;
  if (request.sampling == WVFieldSamplingKind::fixedVerticalProfiles) {
    for (std::size_t profile = 0; profile < request.profileIndices.size(); ++profile)
      for (std::size_t z = 0; z < g.Nz; ++z)
        output.data[z + g.Nz * profile] =
            source[request.profileIndices[profile] + plane * z];
    ++metrics_.profileWriteCount;
  } else if (request.sampling == WVFieldSamplingKind::positions) {
    const auto depth = surface(request.field) ? 1 : g.Nz;
    for (std::size_t position = 0; position < request.weights.size(); ++position)
      output.data[position] = interpolate(source, request.weights[position],
                                          request.interpolation, g.Nx, g.Ny,
                                          depth);
    if (request.interpolation == WVPositionInterpolation::linear)
      metrics_.linearInterpolationCount += request.weights.size();
    else
      metrics_.splineInterpolationCount += request.weights.size();
  } else {
    std::copy_n(source, output.elementCount, output.data);
    ++metrics_.fullGridWriteCount;
  }
  metrics_.outputElementWriteCount += output.elementCount;
  return WVKernelStatus::ok();
}

void WVStratifiedFieldEvaluationAdapter::recordSampledMoving(
    std::size_t positionCount) noexcept {
  ++metrics_.movingEvaluationCount;
  metrics_.movingPositionCount += positionCount;
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::evaluateMoving(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state, WVMovingPositionView positions,
    WVFieldOutputView *outputs, std::size_t outputCount, const std::uint8_t *activeOutputs) {
  return evaluateMovingImpl(publicPlan, state, nullptr, positions, outputs,
                            outputCount, activeOutputs);
}

WVKernelStatus
WVStratifiedFieldEvaluationAdapter::evaluateMovingFromAdvectionFields(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView &advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  return evaluateMovingImpl(publicPlan, state, &advectionFields, positions,
                            outputs, outputCount, activeOutputs);
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::evaluateMovingImpl(
    const WVMovingFieldEvaluationPlan &publicPlan,
    const WVIntegrationState &state,
    const WVRealFieldBundleConstView *advectionFields,
    WVMovingPositionView positions, WVFieldOutputView *outputs,
    std::size_t outputCount, const std::uint8_t *activeOutputs) {
  const auto plan =
      std::static_pointer_cast<const MovingPlan>(publicPlan.transformPlan_);
  if (!plan || plan->fingerprint != configurationFingerprint(configuration()))
    return invalid("Stratified QG moving plan belongs to another transform.");
  if (positions.positionCount != plan->positionCount ||
      (positions.positionCount != 0 &&
       (positions.x == nullptr || positions.y == nullptr)) ||
      outputCount != publicPlan.outputs_.size() ||
      (outputCount != 0 && outputs == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Stratified QG moving positions or outputs have the wrong shape."};
  for (std::size_t output = 0; output < outputCount; ++output)
    if ((!activeOutputs || activeOutputs[output]) && (outputs[output].data == nullptr ||
        outputs[output].elementCount !=
            publicPlan.outputs_[output].elementCount))
      return {WVKernelStatusCode::invalidShape,
              "A Stratified QG moving-field output has the wrong shape."};
  const auto finitePosition=[&](std::size_t index) {return std::isfinite(positions.x[index]) && std::isfinite(positions.y[index]);};
  if(!activeOutputs) {
    for(std::size_t index=0;index<positions.positionCount;++index)
      if(!finitePosition(index)) return invalid("Stratified QG moving positions must be finite.");
  } else {
    bool anyActive=false;
    for(const auto& request:plan->requests) if(activeOutputs[request.output]) {
      anyActive=true;
      for(std::size_t index=request.offset;index<request.offset+request.count;++index)
        if(!finitePosition(index)) return invalid("Stratified QG moving positions must be finite.");
    }
    if(!anyActive) return WVKernelStatus::ok();
  }

  for(const auto& request:plan->requests) {
    if((activeOutputs && !activeOutputs[request.output]) || surface(request.field)) continue;
    if(!positions.z) return invalid("SQG volume samples require z coordinates.");
    for(std::size_t index=request.offset;index<request.offset+request.count;++index)
      if(!std::isfinite(positions.z[index])) return invalid("SQG positions must be finite.");
  }
  WVState amplitudes;
  auto status = coefficientView(state, configuration(), amplitudes);
  if (!status)
    return status;
  const auto& g=configuration();
  if(executing_) return {WVKernelStatusCode::reentrantExecution,"SQG field workspace is active."};
  executing_=true; struct Guard {bool& value;~Guard(){value=false;}} guard{executing_};
  const auto R=g.Nx*g.Ny*g.Nz; auto& workspace=*movingInterpolation_;
  if(advectionFields && (advectionFields->data==nullptr || advectionFields->shape.first!=g.Nx || advectionFields->shape.second!=g.Ny || advectionFields->shape.third!=g.Nz || advectionFields->shape.fourth!=3)) return invalid("SQG advection fields require [Nx,Ny,Nz,3].");
  for(const auto& request:plan->requests) if((!activeOutputs || activeOutputs[request.output]) && !surface(request.field) && !positions.z) return invalid("SQG volume samples require z coordinates.");
  std::array<bool,16> evaluated{};
  for(const auto& request:plan->requests) {
    if(activeOutputs && !activeOutputs[request.output]) continue;
    const auto fieldIndex=static_cast<std::size_t>(request.field);if(evaluated[fieldIndex]) continue;
    const double* values=fieldScratch_.data();
    if(advectionFields && (request.field==WVHydrostaticField::u || request.field==WVHydrostaticField::v || request.field==WVHydrostaticField::w)) {
      values=advectionFields->data+static_cast<std::size_t>(request.field)*R; ++metrics_.primitiveFieldReuseCount;
    } else {
      bool reused=false;
      status=transformField(amplitudes,request.field,{fieldScratch_.data(),{g.Nx,g.Ny,surface(request.field)?1:g.Nz}},&reused);if(!status) return status;
      if(reused) ++metrics_.primitiveFieldReuseCount;
      else {++metrics_.transformCount; ++metrics_.movingPrimitiveTransformCount;}
    }
    evaluated[fieldIndex]=true;
    for(const auto& destination:plan->requests) if((!activeOutputs || activeOutputs[destination.output]) && destination.field==request.field) {
      for(std::size_t p=0;p<destination.count;++p) {
        const auto index=p+destination.offset; const double x=positions.x[index],y=positions.y[index];
        const bool horizontal=surface(request.field); const double z=horizontal?g.z.front():positions.z[index];
        if(!std::isfinite(z)) return invalid("SQG positions must be finite.");
        double value=0;
        if(z>=g.z.front() && z<=g.z.back()) {
          const double dx=g.Lx/g.Nx,dy=g.Ly/g.Ny,xx=wrapped(x,g.Lx),yy=wrapped(y,g.Ly);
          const auto ix=std::min(static_cast<std::size_t>(xx/dx),g.Nx-1),iy=std::min(static_cast<std::size_t>(yy/dy),g.Ny-1);
          const auto iz=horizontal?0:std::min(static_cast<std::size_t>(std::upper_bound(g.z.begin(),g.z.end(),z)-g.z.begin()-1),g.Nz-2);
          if(destination.interpolation==WVPositionInterpolation::linear) {
            const double fx=xx/dx-ix,fy=yy/dy-iy,fz=horizontal?0:(z-g.z[iz])/(g.z[iz+1]-g.z[iz]);
            for(std::size_t k=0;k<(horizontal?1:2);++k) for(std::size_t j=0;j<2;++j) for(std::size_t i=0;i<2;++i)
              value+=(i?fx:1-fx)*(j?fy:1-fy)*(k?fz:1-fz)*values[(ix+i)%g.Nx+g.Nx*((iy+j)%g.Ny+g.Ny*(iz+k))];
            ++metrics_.linearInterpolationCount;
          } else {
            const bool bx=g.Nx<4 || ix<3 || ix>g.Nx-4,by=g.Ny<4 || iy<3 || iy>g.Ny-4;
            workspace.xSpline.weightsInto(dx,bx?wrapped(x+4*dx,g.Lx):xx,bx?4:0,workspace.xWeights,workspace.xShifted,workspace.xRhs);
            workspace.ySpline.weightsInto(dy,by?wrapped(y+4*dy,g.Ly):yy,by?4:0,workspace.yWeights,workspace.yShifted,workspace.yRhs);
            if(!horizontal) workspace.zSpline.weightsInto(z,workspace.zWeights);
            for(std::size_t k=0;k<(horizontal?1:g.Nz);++k) for(std::size_t j=0;j<g.Ny;++j) for(std::size_t i=0;i<g.Nx;++i)
              value+=workspace.xWeights[i]*workspace.yWeights[j]*(horizontal?1:workspace.zWeights[k])*values[i+g.Nx*(j+g.Ny*k)];
            ++metrics_.splineInterpolationCount;
          }
        }
        outputs[destination.output].data[p]=value;
      }
      metrics_.outputElementWriteCount+=destination.count;
    }
  }
  ++metrics_.evaluationCount;
  ++metrics_.movingEvaluationCount;
  metrics_.movingPositionCount += positions.positionCount;
  metrics_.scratchHighWaterBytes = std::max(
      metrics_.scratchHighWaterBytes,
      fieldScratch_.size() * sizeof(double) + workspace.scratchBytes());
  metrics_.servicePersistentBytes = persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::createEventPlan(
    const std::vector<WVEventFieldRequest> &requests,
    WVEventFieldEvaluationPlan &plan) {
  try {
    auto implementation = std::make_shared<EventPlan>();
    implementation->fingerprint = configurationFingerprint(configuration());
    WVEventFieldEvaluationPlan candidate;
    candidate.outputs_.reserve(requests.size());
    std::set<std::string> identifiers;
    for (std::size_t index = 0; index < requests.size(); ++index) {
      if (requests[index].identifier.empty() || !identifiers.insert(requests[index].identifier).second)
        return invalid("SQG event output identifiers must be nonempty and unique.");
      if (requests[index].positionSetSlot == std::numeric_limits<std::size_t>::max())
        return invalid("SQG event position-set slot exceeds the supported range.");
      WVHydrostaticField field;
      ScalarField scalar;
      auto status = resolveField(requests[index].fieldName, field, scalar,(hydrostaticKernel_!=nullptr || boussinesqKernel_!=nullptr));
      if (!status)
        return status;
      if (scalar != ScalarField::none)
        return invalid("Scalar fields cannot use event position sampling.");
      implementation->positionSetCount = std::max(
          implementation->positionSetCount, requests[index].positionSetSlot + 1);
      implementation->requests.push_back(
          {field, requests[index].positionSetSlot,
           requests[index].interpolation, index});
      candidate.outputs_.push_back(
          {requests[index].identifier, requests[index].fieldName,
           findExecutablePortableVariable(requests[index].fieldName)->identifier, surface(field)?WVPortableNaturalRank::horizontal:WVPortableNaturalRank::volume, 0,
           requests[index].positionSetSlot, requests[index].interpolation});
    }
    candidate.positionSetCount_ = implementation->positionSetCount;
    candidate.fingerprint_ = implementation->fingerprint;
    candidate.transformPlan_ = implementation;
    candidate.transformPlanBytes_ = implementation->persistentBytes();
    plan = std::move(candidate);
    ++metrics_.eventPlanCreationCount;
    metrics_.eventPlanFieldResolutionCount += requests.size();
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate a Stratified QG event-field plan."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::prepareEventGeometry(
    const WVEventFieldEvaluationPlan &publicPlan,
    const WVEventPositionSetView *positionSets,
    std::size_t positionSetCount, WVPreparedFieldGeometry &geometry) {
  const auto eventPlan =
      std::static_pointer_cast<const EventPlan>(publicPlan.transformPlan_);
  if (!eventPlan || eventPlan->fingerprint != configurationFingerprint(configuration()) ||
      positionSetCount != eventPlan->positionSetCount ||
      (positionSetCount != 0 && positionSets == nullptr))
    return {WVKernelStatusCode::invalidShape,
            "Stratified QG event geometry does not match its plan."};
  try {
    std::vector<WVFieldRequest> requests;
    requests.reserve(eventPlan->requests.size());
    WVPreparedFieldGeometry candidate;
    candidate.positionSets_.reserve(positionSetCount);
    candidate.positionCount_ = 0;
    for (std::size_t slot = 0; slot < positionSetCount; ++slot) {
      const auto &view = positionSets[slot];
      if (view.extentCount != 0 && view.extents == nullptr)
        return {WVKernelStatusCode::invalidPointer,"SQG event extents are null."};
      if (view.positionCount != 0 && (view.x == nullptr || view.y == nullptr))
        return {WVKernelStatusCode::invalidPointer,
                "Stratified QG event coordinates are null."};
      WVPreparedFieldGeometry::PositionSet stored;
      stored.x = view.x;
      stored.y = view.y;
      stored.z = view.z;
      stored.positionCount = view.positionCount;
      if (view.extentCount == 0)
        stored.extents = {view.positionCount};
      else
        stored.extents.assign(view.extents, view.extents + view.extentCount);
      std::size_t count = 1;
      for (const auto extent : stored.extents)
        count = checkedProduct(count, extent);
      if (count != view.positionCount)
        return {WVKernelStatusCode::invalidShape,
                "Stratified QG event extents do not match coordinates."};
      candidate.positionCount_ += view.positionCount;
      candidate.borrowedCoordinateBytes_ +=
          3 * view.positionCount * sizeof(double);
      candidate.positionSets_.push_back(std::move(stored));
    }
    for (const auto &request : eventPlan->requests) {
      const auto &set = positionSets[request.positionSet];
      WVFieldSamplingRequest sampling;
      sampling.kind = WVFieldSamplingKind::positions;
      sampling.interpolation = request.interpolation;
      if (set.positionCount != 0) {
        sampling.x.assign(set.x, set.x + set.positionCount);
        sampling.y.assign(set.y, set.y + set.positionCount);
      }
      if(set.z) sampling.z.assign(set.z,set.z+set.positionCount);
      requests.push_back({publicPlan.outputs_[request.output].identifier,
                          publicPlan.outputs_[request.output].fieldName,
                          std::move(sampling)});
      candidate.outputs_.push_back(
          {request.output, request.positionSet,
           candidate.positionSets_[request.positionSet].extents,
           set.positionCount});
    }
    WVFieldEvaluationPlan evaluationPlan;
    auto status = createPlan(requests, evaluationPlan);
    if (!status)
      return status;
    auto geometryImpl = std::make_shared<EventGeometry>();
    geometryImpl->eventPlan = eventPlan;
    geometryImpl->plan =
        std::static_pointer_cast<const Plan>(evaluationPlan.transformPlan_);
    candidate.fieldPlanFingerprint_ = publicPlan.fingerprint_;
    candidate.geometryFingerprint_ = publicPlan.fingerprint_;
    for (const auto &set : candidate.positionSets_) {
      const auto byteCount = set.positionCount * sizeof(double);
      for (const auto *data : {set.x, set.y, set.z}) {
        if(!data) continue;
        const auto *bytes = reinterpret_cast<const std::uint8_t *>(data);
        for (std::size_t byte = 0; byte < byteCount; ++byte) {
          candidate.geometryFingerprint_ ^= bytes[byte];
          candidate.geometryFingerprint_ *= 1099511628211ULL;
        }
      }
    }
    candidate.transformGeometry_ = geometryImpl;
    candidate.transformGeometryBytes_ = geometryImpl->persistentBytes();
    candidate.evaluationPlan_ = std::move(evaluationPlan);
    geometry = std::move(candidate);
    ++metrics_.eventGeometryPreparationCount;
    metrics_.eventPositionSetCount += positionSetCount;
    metrics_.eventPositionCount += geometry.positionCount_;
    metrics_.lastPreparedGeometryRetainedBytes = geometry.retainedBytes();
    metrics_.maximumPreparedGeometryRetainedBytes = std::max(
        metrics_.maximumPreparedGeometryRetainedBytes,
        metrics_.lastPreparedGeometryRetainedBytes);
    return WVKernelStatus::ok();
  } catch (const std::overflow_error &error) {
    return {WVKernelStatusCode::sizeOverflow, error.what()};
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Unable to allocate Stratified QG event geometry."};
  }
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::evaluateEventBatch(
    const WVIntegrationState &state,
    const WVEventFieldEvaluationBatchEntry *entries,
    std::size_t entryCount) {
  if (entryCount != 0 && entries == nullptr)
    return {WVKernelStatusCode::invalidPointer,
            "Stratified QG event batch entries are null."};
  for (std::size_t entry = 0; entry < entryCount; ++entry) {
    if (entries[entry].plan == nullptr || entries[entry].geometry == nullptr)
      return {WVKernelStatusCode::invalidPointer,
              "Stratified QG event entry is incomplete."};
    const auto geometry = std::static_pointer_cast<const EventGeometry>(entries[entry].geometry->transformGeometry_);
    if (!geometry || geometry->eventPlan.get() != entries[entry].plan->transformPlan_.get())
      return invalid("SQG event geometry was prepared for a different field plan.");
    const auto status = evaluate(
        entries[entry].geometry->evaluationPlan_, state,
        entries[entry].outputs, entries[entry].outputCount);
    if (!status)
      return status;
  }
  ++metrics_.eventBatchEvaluationCount;
  metrics_.eventBatchOccurrenceCount += entryCount;
  metrics_.eventEvaluationCount += entryCount;
  return WVKernelStatus::ok();
}

bool WVStratifiedFieldEvaluationAdapter::isCompatibleWith(const WVIntegrationStateLayout& layout) const noexcept {
  const auto& g=configuration(); const bool hydro=(hydrostaticKernel_!=nullptr || boussinesqKernel_!=nullptr);
  if (layout.transformIdentifier()!=g.transformClass || layout.spatialDimensions()!=std::vector<std::size_t>{g.Nx,g.Ny,g.Nz} || layout.coefficientFamilyCount()!=(hydro?3U:1U)) return false;
  const char* names[]={"Ap","Am","A0"};
  for (std::size_t i=0;i<layout.coefficientFamilyCount();++i) if (layout.coefficientFamilies()[i].identifier!=names[hydro?i:2] || layout.coefficientFamilies()[i].elementCount!=g.Nj*g.Nkl) return false;
  return true;
}

const WVStratifiedModalGeometry &
WVStratifiedFieldEvaluationAdapter::configuration() const noexcept {
  return boussinesqKernel_ ? boussinesqKernel_->geometry() : hydrostaticKernel_ ? hydrostaticKernel_->geometry() : kernel_->geometry();
}

std::size_t
WVStratifiedFieldEvaluationAdapter::persistentBytes() const noexcept {
  return sizeof(*this) +
         (ownedKernel_ ? kernel_->persistentBytes() : 0) + (ownedHydrostatic_ ? ownedHydrostatic_->persistentBytes() : 0) + (ownedBoussinesq_ ? ownedBoussinesq_->persistentBytes() : 0) + speedScratch_.capacity()*sizeof(double) +
         fieldScratch_.capacity() * sizeof(double) +
         (movingInterpolation_ ? movingInterpolation_->persistentBytes() : 0);
}

WVKernelStatus WVStratifiedFieldEvaluationAdapter::transformField(const WVState& state,WVHydrostaticField field,WVRealVolumeView out,bool* wasReused) {
  if(eventWorkspace_ && surface(field)) {
    // Surface operations select the upper plane of the same volume kernel.
    // Keep that expensive volume in the event workspace, not another alias.
    const auto base=field==WVHydrostaticField::ssu ? WVHydrostaticField::u :
        field==WVHydrostaticField::ssv ? WVHydrostaticField::v : WVHydrostaticField::pi;
    const auto& g=configuration();
    const auto status=transformField(state,base,{fieldScratch_.data(),{g.Nx,g.Ny,g.Nz}},wasReused);
    if(!status) return status;
    const auto plane=g.Nx*g.Ny;
    std::copy_n(fieldScratch_.data()+plane*(g.Nz-1),plane,out.data);
    return WVKernelStatus::ok();
  }
  bool reused=false;
  const auto operation=[&](){return transformUncachedField(state,field,out);};
  const WVVariableEvaluationKey key{WVVariableEvaluationNode::physicalField,
      static_cast<std::uint32_t>(portableField(field)),
      eventWorkspace_ ? eventWorkspace_->component() : 0u};
  const auto status=eventWorkspace_ ? eventWorkspace_->evaluate(key,out.data,
      out.shape.elementCount(),operation,reused) : operation();
  if(wasReused) *wasReused=reused;
  return status;
}
WVKernelStatus WVStratifiedFieldEvaluationAdapter::transformUncachedField(const WVState& state,WVHydrostaticField field,WVRealVolumeView out) {
  if(eventWorkspace_ && (field==WVHydrostaticField::zetaX ||
      field==WVHydrostaticField::zetaY) &&
      (hydrostaticKernel_ || boussinesqKernel_)) {
    const auto shape=out.shape;
    const auto component=eventWorkspace_->component();
    const auto evaluateDerivative=[&](WVHydrostaticField source,
        std::uint32_t derivative,std::vector<double>& storage) {
      const WVVariableEvaluationKey key{WVVariableEvaluationNode::derivative,
          static_cast<std::uint32_t>(portableField(source)),component,
          derivative};
      bool reused=false;
      WVRealVolumeView destination{storage.data(),shape};
      const auto operation=[&]() {
        if(hydrostaticKernel_)
          return hydrostaticKernel_->transformStateField(state,source,
              destination,static_cast<WVHydrostaticDerivative>(derivative));
        WVBoussinesqField mapped=source==WVHydrostaticField::u ?
            WVBoussinesqField::u : source==WVHydrostaticField::v ?
            WVBoussinesqField::v : WVBoussinesqField::w;
        return boussinesqKernel_->transformStateField(state,mapped,destination,
            static_cast<WVBoussinesqDerivative>(derivative));
      };
      return eventWorkspace_->evaluate(key,storage.data(),storage.size(),
          operation,reused);
    };
    const auto firstField=field==WVHydrostaticField::zetaX ?
        WVHydrostaticField::w : WVHydrostaticField::u;
    const auto firstDerivative=field==WVHydrostaticField::zetaX ?
        static_cast<std::uint32_t>(WVHydrostaticDerivative::y) :
        static_cast<std::uint32_t>(WVHydrostaticDerivative::z);
    const auto secondField=field==WVHydrostaticField::zetaX ?
        WVHydrostaticField::v : WVHydrostaticField::w;
    const auto secondDerivative=field==WVHydrostaticField::zetaX ?
        static_cast<std::uint32_t>(WVHydrostaticDerivative::z) :
        static_cast<std::uint32_t>(WVHydrostaticDerivative::x);
    auto status=evaluateDerivative(firstField,firstDerivative,fieldScratch_);
    if(!status) return status;
    status=evaluateDerivative(secondField,secondDerivative,speedScratch_);
    if(!status) return status;
    const WVRealVolumeConstView first{fieldScratch_.data(),shape};
    const WVRealVolumeConstView second{speedScratch_.data(),shape};
    if(hydrostaticKernel_)
      return hydrostaticKernel_->combinePreparedHorizontalVorticity(
          field,first,second,out,
          static_cast<WVHydrostaticComponent>(component));
    const auto target=field==WVHydrostaticField::zetaX ?
        WVBoussinesqField::zetaX : WVBoussinesqField::zetaY;
    return boussinesqKernel_->combinePreparedHorizontalVorticity(
        target,first,second,out,
        static_cast<WVBoussinesqComponent>(component));
  }
  if (hydrostaticKernel_) return hydrostaticKernel_->transformStateField(state,field,out);
  if (boussinesqKernel_) {
    WVBoussinesqField mapped;
    switch(field) {
      case WVHydrostaticField::u: mapped=WVBoussinesqField::u; break;
      case WVHydrostaticField::v: mapped=WVBoussinesqField::v; break;
      case WVHydrostaticField::w: mapped=WVBoussinesqField::w; break;
      case WVHydrostaticField::eta: mapped=WVBoussinesqField::eta; break;
      case WVHydrostaticField::pi: mapped=WVBoussinesqField::pi; break;
      case WVHydrostaticField::p: mapped=WVBoussinesqField::p; break;
      case WVHydrostaticField::psi: mapped=WVBoussinesqField::psi; break;
      case WVHydrostaticField::qgpv: mapped=WVBoussinesqField::qgpv; break;
      case WVHydrostaticField::rhoE: mapped=WVBoussinesqField::rhoE; break;
      case WVHydrostaticField::rhoTotal: mapped=WVBoussinesqField::rhoTotal; break;
      case WVHydrostaticField::zetaX: mapped=WVBoussinesqField::zetaX; break;
      case WVHydrostaticField::zetaY: mapped=WVBoussinesqField::zetaY; break;
      case WVHydrostaticField::zetaZ: mapped=WVBoussinesqField::zetaZ; break;
      case WVHydrostaticField::ssh: mapped=WVBoussinesqField::ssh; break;
      case WVHydrostaticField::ssu: mapped=WVBoussinesqField::ssu; break;
      case WVHydrostaticField::ssv: mapped=WVBoussinesqField::ssv; break;
      default: return invalid("Unsupported Boussinesq field.");
    }
    return boussinesqKernel_->transformStateField(state,mapped,out);
  }
  WVStratifiedQGField qg;
  switch(field) {
    case WVHydrostaticField::u: qg=WVStratifiedQGField::u; break;
    case WVHydrostaticField::v: qg=WVStratifiedQGField::v; break;
    case WVHydrostaticField::w: qg=WVStratifiedQGField::w; break;
    case WVHydrostaticField::eta: qg=WVStratifiedQGField::eta; break;
    case WVHydrostaticField::pi: qg=WVStratifiedQGField::pi; break;
    case WVHydrostaticField::p: qg=WVStratifiedQGField::p; break;
    case WVHydrostaticField::psi: qg=WVStratifiedQGField::psi; break;
    case WVHydrostaticField::qgpv: qg=WVStratifiedQGField::qgpv; break;
    case WVHydrostaticField::rhoE: qg=WVStratifiedQGField::rhoE; break;
    case WVHydrostaticField::rhoTotal: qg=WVStratifiedQGField::rhoTotal; break;
    case WVHydrostaticField::zetaZ: qg=WVStratifiedQGField::zetaZ; break;
    case WVHydrostaticField::ssh: qg=WVStratifiedQGField::ssh; break;
    case WVHydrostaticField::ssu: qg=WVStratifiedQGField::ssu; break;
    case WVHydrostaticField::ssv: qg=WVStratifiedQGField::ssv; break;
    default: return invalid("Unsupported Stratified QG field.");
  }
  return kernel_->transformA0ToField(state.coefficients.A0,qg,out);
}
WVKernelStatus WVStratifiedFieldEvaluationAdapter::scalarValue(const WVState& state,unsigned scalar,double& value) {
  if (!hydrostaticKernel_ && !boussinesqKernel_) return scalar==static_cast<unsigned>(ScalarField::energy) ? kernel_->totalEnergy(state.coefficients.A0,value) : kernel_->uvMax(state.coefficients.A0,value);
  if (scalar==static_cast<unsigned>(ScalarField::energy)) return boussinesqKernel_ ? boussinesqKernel_->totalEnergy(state.coefficients,value) : hydrostaticKernel_->totalEnergy(state.coefficients,value);
  const auto shape=boussinesqKernel_ ? boussinesqKernel_->spatialShape() : hydrostaticKernel_->spatialShape();
  auto s=transformField(state,scalar==static_cast<unsigned>(ScalarField::wMax) ? WVHydrostaticField::w : WVHydrostaticField::u,{fieldScratch_.data(),shape}); if (!s) return s;
  if (scalar==static_cast<unsigned>(ScalarField::uvMax)) { s=transformField(state,WVHydrostaticField::v,{speedScratch_.data(),shape}); if (!s) return s; }
  if (hydrostaticKernel_ && scalar==static_cast<unsigned>(ScalarField::uvMax))
    return hydrostaticKernel_->reduceHorizontalSpeedMaximum(
        {fieldScratch_.data(),shape},{speedScratch_.data(),shape},value);
  value=0; for (std::size_t i=0;i<fieldScratch_.size();++i) value=std::max(value,scalar==static_cast<unsigned>(ScalarField::wMax) ? std::abs(fieldScratch_[i]) : std::hypot(fieldScratch_[i],speedScratch_[i]));
  return WVKernelStatus::ok();
}
} // namespace wavevortex::runtime::detail
