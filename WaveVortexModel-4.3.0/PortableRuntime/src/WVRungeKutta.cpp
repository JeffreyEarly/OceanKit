#include "WaveVortexRuntime/WVRungeKutta.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVComplex64 scaledSum(WVComplex64 value, WVComplex64 increment,
                      double scale) noexcept {
  return {value.real + scale * increment.real,
          value.imag + scale * increment.imag};
}

double timeTolerance(double first, double second) noexcept {
  return 8.0 * std::numeric_limits<double>::epsilon() *
         std::max({1.0, std::abs(first), std::abs(second)});
}

std::uint64_t hashTolerance(std::uint64_t hash, double value) noexcept {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  // Discard inconsequential low mantissa noise so independently evaluated
  // MATLAB and C++ tolerance formulas have one reproducible audit identity.
  // Clearing 20 bits retains approximately 32 bits of significand precision;
  // the production-size formula audit differs by at most 6.3e-16 relatively.
  bits = (bits + UINT64_C(0x80000)) & ~UINT64_C(0xfffff);
  hash ^= bits;
  return hash * UINT64_C(1099511628211);
}

class IntegrationBuffer {
public:
  WVKernelStatus initialize(const WVIntegrationStateLayout &layout) {
    try {
      layout_ = &layout;
      coefficientCount_ = layout.coefficientElementCount();
      complex_.assign(coefficientCount_ + layout.complexElementCount(),
                      WVComplex64{});
      real_.assign(layout.realElementCount(), 0.0);
      mutableCoefficientFamilies_.clear();
      constCoefficientFamilies_.clear();
      mutableCoefficientFamilies_.reserve(layout.coefficientFamilyCount());
      constCoefficientFamilies_.reserve(layout.coefficientFamilyCount());
      for (const auto &family : layout.coefficientFamilies()) {
        auto *data = complex_.data() + family.scalarOffset;
        mutableCoefficientFamilies_.push_back({&family, data});
        constCoefficientFamilies_.push_back({&family, data});
      }
      mutableBlocks_.clear();
      constBlocks_.clear();
      mutableBlocks_.reserve(layout.additionalBlocks().size());
      constBlocks_.reserve(layout.additionalBlocks().size());
      for (const auto &block : layout.additionalBlocks()) {
        auto *real = block.scalarType == WVStateScalarType::real64
                         ? real_.data() + block.scalarOffset
                         : nullptr;
        auto *complex =
            block.scalarType == WVStateScalarType::complex64
                ? complex_.data() + coefficientCount_ + block.scalarOffset
                : nullptr;
        mutableBlocks_.push_back({&block, real, complex});
        constBlocks_.push_back({&block, real, complex});
      }
      return WVKernelStatus::ok();
    } catch (const std::bad_alloc &) {
      return {WVKernelStatusCode::allocationFailure,
              "Integration workspace allocation failed."};
    }
  }

  WVMutableIntegrationState mutableState(double t, double t0) noexcept {
    WVMutableState legacy;
    legacy.t = t;
    legacy.t0 = t0;
    if (layout_->hasLegacyCoefficientTriple()) {
      const auto shape = layout_->coefficientShape();
      legacy.coefficients =
          {{mutableCoefficientFamilies_[0].data, shape},
           {mutableCoefficientFamilies_[1].data, shape},
           {mutableCoefficientFamilies_[2].data, shape}};
    }
    return {legacy, mutableBlocks_.data(), mutableBlocks_.size(),
            mutableCoefficientFamilies_.data(),
            mutableCoefficientFamilies_.size()};
  }
  WVIntegrationState state(double t, double t0) const noexcept {
    WVState legacy;
    legacy.t = t;
    legacy.t0 = t0;
    if (layout_->hasLegacyCoefficientTriple()) {
      const auto shape = layout_->coefficientShape();
      legacy.coefficients =
          {{constCoefficientFamilies_[0].data, shape},
           {constCoefficientFamilies_[1].data, shape},
           {constCoefficientFamilies_[2].data, shape}};
    }
    return {legacy, constBlocks_.data(), constBlocks_.size(),
            constCoefficientFamilies_.data(),
            constCoefficientFamilies_.size()};
  }
  WVIntegrationFlux flux() noexcept {
    WVFlux legacy;
    if (layout_->hasLegacyCoefficientTriple()) {
      const auto shape = layout_->coefficientShape();
      legacy = {{mutableCoefficientFamilies_[0].data, shape},
                {mutableCoefficientFamilies_[1].data, shape},
                {mutableCoefficientFamilies_[2].data, shape}};
    }
    return {legacy, mutableBlocks_.data(), mutableBlocks_.size(),
            mutableCoefficientFamilies_.data(),
            mutableCoefficientFamilies_.size()};
  }
  void copyFrom(const WVIntegrationState &source) noexcept {
    for (std::size_t family = 0; family < layout_->coefficientFamilyCount();
         ++family) {
      const auto view = coefficientFamilyView(*layout_, source, family);
      const auto &metadata = layout_->coefficientFamilies()[family];
      std::copy_n(view.data, metadata.elementCount,
                  complex_.data() + metadata.scalarOffset);
    }
    for (std::size_t block = 0; block < source.additionalBlockCount; ++block) {
      const auto &layout = *source.additionalBlocks[block].layout;
      if (layout.scalarType == WVStateScalarType::real64)
        std::copy_n(source.additionalBlocks[block].realData,
                    layout.elementCount, real_.data() + layout.scalarOffset);
      else
        std::copy_n(
            source.additionalBlocks[block].complexData, layout.elementCount,
            complex_.data() + coefficientCount_ + layout.scalarOffset);
    }
  }
  void copyTo(WVMutableIntegrationState &destination) const noexcept {
    for (std::size_t family = 0; family < layout_->coefficientFamilyCount();
         ++family) {
      const auto view = coefficientFamilyView(*layout_, destination, family);
      const auto &metadata = layout_->coefficientFamilies()[family];
      std::copy_n(complex_.data() + metadata.scalarOffset,
                  metadata.elementCount, view.data);
    }
    for (std::size_t block = 0; block < destination.additionalBlockCount;
         ++block) {
      const auto &layout = *destination.additionalBlocks[block].layout;
      if (layout.scalarType == WVStateScalarType::real64)
        std::copy_n(real_.data() + layout.scalarOffset, layout.elementCount,
                    destination.additionalBlocks[block].realData);
      else
        std::copy_n(complex_.data() + coefficientCount_ +
                        layout.scalarOffset,
                    layout.elementCount,
                    destination.additionalBlocks[block].complexData);
    }
  }
  void assign(const IntegrationBuffer &source) noexcept {
    complex_ = source.complex_;
    real_ = source.real_;
  }
  void setScaled(const IntegrationBuffer &source, double scale) noexcept {
    for (std::size_t i = 0; i < complex_.size(); ++i) {
      complex_[i].real = scale * source.complex_[i].real;
      complex_[i].imag = scale * source.complex_[i].imag;
    }
    for (std::size_t i = 0; i < real_.size(); ++i)
      real_[i] = scale * source.real_[i];
  }
  void setAffine(const WVIntegrationState &base, const IntegrationBuffer &increment,
                 double scale) noexcept {
    for (std::size_t family = 0; family < layout_->coefficientFamilyCount();
         ++family) {
      const auto coefficients = coefficientFamilyView(*layout_, base, family);
      const auto &metadata = layout_->coefficientFamilies()[family];
      for (std::size_t index = 0; index < metadata.elementCount; ++index) {
        const auto flatIndex = metadata.scalarOffset + index;
        complex_[flatIndex] = scaledSum(coefficients.data[index],
                                        increment.complex_[flatIndex], scale);
      }
    }
    for (std::size_t block = 0; block < base.additionalBlockCount; ++block) {
      const auto &layout = *base.additionalBlocks[block].layout;
      if (layout.scalarType == WVStateScalarType::real64) {
        for (std::size_t index = 0; index < layout.elementCount; ++index) {
          const auto flatIndex = layout.scalarOffset + index;
          real_[flatIndex] = base.additionalBlocks[block].realData[index] +
                             scale * increment.real_[flatIndex];
        }
      } else {
        for (std::size_t index = 0; index < layout.elementCount; ++index) {
          const auto flatIndex =
              coefficientCount_ + layout.scalarOffset + index;
          complex_[flatIndex] =
              scaledSum(base.additionalBlocks[block].complexData[index],
                        increment.complex_[flatIndex], scale);
        }
      }
    }
  }
  void addScaled(const IntegrationBuffer &source, double scale) noexcept {
    for (std::size_t i = 0; i < complex_.size(); ++i)
      complex_[i] = scaledSum(complex_[i], source.complex_[i], scale);
    for (std::size_t i = 0; i < real_.size(); ++i)
      real_[i] += scale * source.real_[i];
  }
  void setWeightedCandidate(const WVIntegrationState &base, double h,
                            const IntegrationBuffer &k1, double w1,
                            const IntegrationBuffer &k2, double w2,
                            const IntegrationBuffer &k3, double w3) noexcept {
    copyFrom(base);
    addScaled(k1, h * w1);
    addScaled(k2, h * w2);
    addScaled(k3, h * w3);
  }
  std::size_t capacityBytes() const noexcept {
    return complex_.capacity() * sizeof(WVComplex64) +
           real_.capacity() * sizeof(double) +
           mutableBlocks_.capacity() * sizeof(WVAdditionalStateBlockView) +
           constBlocks_.capacity() * sizeof(WVAdditionalStateBlockConstView);
  }
  std::size_t valueCapacityBytes() const noexcept {
    return complex_.capacity() * sizeof(WVComplex64) +
           real_.capacity() * sizeof(double);
  }
  std::size_t coefficientViewBytes() const noexcept {
    return mutableCoefficientFamilies_.capacity() *
               sizeof(WVCoefficientFamilyView) +
           constCoefficientFamilies_.capacity() *
               sizeof(WVCoefficientFamilyConstView);
  }
  const std::vector<WVComplex64> &complex() const noexcept { return complex_; }
  const std::vector<double> &real() const noexcept { return real_; }
  std::vector<WVComplex64> &complex() noexcept { return complex_; }
  std::vector<double> &real() noexcept { return real_; }
  std::size_t coefficientCount() const noexcept { return coefficientCount_; }

private:
  const WVIntegrationStateLayout *layout_ = nullptr;
  std::size_t coefficientCount_ = 0;
  std::vector<WVComplex64> complex_;
  std::vector<double> real_;
  std::vector<WVCoefficientFamilyView> mutableCoefficientFamilies_;
  std::vector<WVCoefficientFamilyConstView> constCoefficientFamilies_;
  std::vector<WVAdditionalStateBlockView> mutableBlocks_;
  std::vector<WVAdditionalStateBlockConstView> constBlocks_;
};

WVKernelStatus constrain(WVIntegrationSystem &system,
                         IntegrationBuffer &buffer, double t, double t0) {
  auto state = buffer.mutableState(t, t0);
  return system.enforceStateConstraints(state).status;
}

WVKernelStatus evaluate(WVIntegrationSystem &system,
                        const IntegrationBuffer &state, double t, double t0,
                        IntegrationBuffer &derivative,
                        WVIntegratorMetrics &metrics) {
  auto flux = derivative.flux();
  const auto status = system.evaluateRightHandSide(state.state(t, t0), flux);
  if (status)
    ++metrics.rightHandSideEvaluationCount;
  return status;
}

void makeExternalViews(
    const WVMutableIntegrationState &state,
    std::vector<WVCoefficientFamilyConstView> &coefficientViews,
    std::vector<WVAdditionalStateBlockConstView> &blockViews) {
  coefficientViews.clear();
  coefficientViews.reserve(state.coefficientFamilyCount);
  for (std::size_t index = 0; index < state.coefficientFamilyCount; ++index)
    coefficientViews.push_back({state.coefficientFamilies[index].layout,
                                state.coefficientFamilies[index].data});
  blockViews.clear();
  blockViews.reserve(state.additionalBlockCount);
  for (std::size_t index = 0; index < state.additionalBlockCount; ++index)
    blockViews.push_back({state.additionalBlocks[index].layout,
                          state.additionalBlocks[index].realData,
                          state.additionalBlocks[index].complexData});
}

struct RK23ControllerPolicy {
  static double acceptedErrorPower(double error) noexcept {
    return std::cbrt(error);
  }
  static double rejectedErrorPower(double error) noexcept {
    return std::pow(error, -1.0 / 3.0);
  }
};

struct RK45ControllerPolicy {
  static double acceptedErrorPower(double error) noexcept {
    return std::pow(error, 1.0 / 5.0);
  }
  static double rejectedErrorPower(double error) noexcept {
    return std::pow(error, -1.0 / 5.0);
  }
};

struct RK78MethodPolicy {
  static constexpr const char *controllerIdentifier = "matlab-ode78-v1";
  static constexpr const char *methodIdentifier = "adaptive-rk78";

  static constexpr double c2 = 0.05;
  static constexpr double c3 = 0.1065625;
  static constexpr double c4 = 0.15984375;
  static constexpr double c5 = 0.39;
  static constexpr double c6 = 0.465;
  static constexpr double c7 = 0.155;
  static constexpr double c8 = 0.943;
  static constexpr double c9 = 0.901802041735857;
  static constexpr double c10 = 0.909;
  static constexpr double c11 = 0.94;

  static constexpr double a21 = 0.05;
  static constexpr double a31 = -0.0069931640625;
  static constexpr double a32 = 0.1135556640625;
  static constexpr double a41 = 0.0399609375;
  static constexpr double a43 = 0.1198828125;
  static constexpr double a51 = 0.36139756280045754;
  static constexpr double a53 = -1.3415240667004928;
  static constexpr double a54 = 1.3701265039000352;
  static constexpr double a61 = 0.049047202797202795;
  static constexpr double a64 = 0.23509720422144048;
  static constexpr double a65 = 0.18085559298135673;
  static constexpr double a71 = 0.06169289044289044;
  static constexpr double a74 = 0.11236568314640277;
  static constexpr double a75 = -0.03885046071451367;
  static constexpr double a76 = 0.01979188712522046;
  static constexpr double a81 = -1.767630240222327;
  static constexpr double a84 = -62.5;
  static constexpr double a85 = -6.061889377376669;
  static constexpr double a86 = 5.6508231982227635;
  static constexpr double a87 = 65.62169641937624;
  static constexpr double a91 = -1.1809450665549708;
  static constexpr double a94 = -41.50473441114321;
  static constexpr double a95 = -4.434438319103725;
  static constexpr double a96 = 4.260408188586133;
  static constexpr double a97 = 43.75364022446172;
  static constexpr double a98 = 0.00787142548991231;
  static constexpr double a101 = -1.2814059994414884;
  static constexpr double a104 = -45.047139960139866;
  static constexpr double a105 = -4.731362069449577;
  static constexpr double a106 = 4.514967016593808;
  static constexpr double a107 = 47.44909557172985;
  static constexpr double a108 = 0.010592282971116612;
  static constexpr double a109 = -0.0057468422638446166;
  static constexpr double a111 = -1.7244701342624853;
  static constexpr double a114 = -60.92349008483054;
  static constexpr double a115 = -5.951518376222393;
  static constexpr double a116 = 5.556523730698456;
  static constexpr double a117 = 63.98301198033305;
  static constexpr double a118 = 0.014642028250414961;
  static constexpr double a119 = 0.06460408772358203;
  static constexpr double a1110 = -0.0793032316900888;
  static constexpr double a121 = -3.301622667747079;
  static constexpr double a124 = -118.01127235975251;
  static constexpr double a125 = -10.141422388456112;
  static constexpr double a126 = 9.139311332232058;
  static constexpr double a127 = 123.37594282840426;
  static constexpr double a128 = 4.62324437887458;
  static constexpr double a129 = -3.3832777380682018;
  static constexpr double a1210 = 4.527592100324618;
  static constexpr double a1211 = -5.828495485811623;
  static constexpr double a131 = -3.039515033766309;
  static constexpr double a134 = -109.26086808941763;
  static constexpr double a135 = -9.290642497400293;
  static constexpr double a136 = 8.43050498176491;
  static constexpr double a137 = 114.20100103783314;
  static constexpr double a138 = -0.9637271342145479;
  static constexpr double a139 = -5.0348840888021895;
  static constexpr double a1310 = 5.958130824002923;

  static constexpr double b1 = 0.04427989419007951;
  static constexpr double b6 = 0.3541049391724449;
  static constexpr double b7 = 0.2479692154956438;
  static constexpr double b8 = -15.694202038838084;
  static constexpr double b9 = 25.084064965558564;
  static constexpr double b10 = -31.738367786260277;
  static constexpr double b11 = 22.938283273988784;
  static constexpr double b12 = -0.2361324633071542;

  static constexpr double e1 = 3.272103901028776e-05;
  static constexpr double e6 = 0.0005046250618777735;
  static constexpr double e7 = -0.00012117235897844563;
  static constexpr double e8 = 20.142336771313868;
  static constexpr double e9 = -5.237178599439828;
  static constexpr double e10 = 8.156744408794658;
  static constexpr double e11 = -22.938283273988784;
  static constexpr double e12 = 0.2361324633071542;
  static constexpr double e13 = -0.36016794372897754;

  static constexpr double c15 = 0.3110177634953864;
  static constexpr double c16 = 0.1725;
  static constexpr double c17 = 0.7846;

  inline static constexpr double a15[] = {
      0.04620700646754963,     0.045039041608424805,
      0.23368166977134244,     37.83901368421068,
      -15.949113289454246,     23.028368351816102,
      -44.85578507769412,      -0.06379858768647444,
      -0.012595035543861663};
  inline static constexpr double a16[] = {
      0.05037946855482041,     0.041098361310460796,
      0.17180541533481958,     4.614105319981519,
      -1.7916678830853965,     2.531658930485041,
      -5.324977860205731,      -0.03065532595385635,
      -0.005254479979429613,   -0.08399194644224793};
  inline static constexpr double a17[] = {
      0.0408289713299708,      0.4244479514247632,
      0.23260915312752345,     2.677982520711806,
      0.7420826657338945,      0.1460377847941461,
      -3.579344509890565,      0.11388443896001738,
      0.012677906510331901,    -0.07443436349946675,
      0.047827480797578516};

  inline static constexpr double continuousWeights[6][12] = {
      {-7.238550783576432811855355839508646327161,
       11.15330887588935170976376962782446833855,
       2.34875229807309355640904629061136935335,
       -1027.321675339240679090464776362465090654,
       1568.546608927281956416687915664731868885,
       -2000.882061921041961546811133479107090218,
       1496.620400693446268810344884971434468267,
       -16.41320775560933621675902845723196069900,
       -4.29672443178246482824254064733546854251,
       -20.41628069294821485579834313809132051248,
       16.53007184264271512356106095760699278945,
       -18.63064171313429626683549958846959067803},
      {26.00913483254676138219215542805486438340,
       -91.7609656398961659890179437322816238711,
       -11.6724894172018429369093778842231443146,
       9198.71432360760879019681406218311101879,
       -13995.38852541600542155322174511897930298,
       17864.36380347691630038038755096765127729,
       -13397.55405171476021512904990709508924800,
       147.6097045407002371315249807692915435608,
       38.6444746111678092366406218271498656093,
       153.5213232524836445391962375168798263930,
       -96.6861433615782065041742809436987893361,
       164.1994112280183092456176460821337125030},
      {-50.23684777762566731759165474184543812128,
       291.7074241722059450113911477530513089255,
       -3.339139076505928386509206543237093540,
       -33189.78048157363822223641020734287802492,
       50256.2124698102445419491620666726469821,
       -64205.1907515562863000297926577113695108,
       48323.5602199437493999696912750109765015,
       -535.719963714732106447158760197417632645,
       -140.3503471762808981414524290552248895548,
       -436.5502610211220460266289847121377276100,
       268.959934219531723149495873437076657635,
       -579.272256249540441494196462569641132906},
      {52.12072084601022449485077581012685809554,
       -430.4096692910862817449451677633631387823,
       94.885262249720610030798242337479596095,
       57750.0831348887181073584126028277545727,
       -86974.5128036219909523950692144595063700,
       111224.8489930378077126420609392735999202,
       -84051.4283423393032636942266780744607468,
       938.286247077820650371318861625025573381,
       246.3954669697502467443139611011701827640,
       598.214644262650861959065070073603792110,
       -428.681909788964647271837835032326719249,
       980.198255708866731505258442280896479501},
      {-27.06472451211777193118825764262673140465,
       299.4531188198997479843407054776900024282,
       -143.071126583012024456409244370652716962,
       -47698.93315706261990169947144294597707756,
       71494.7977095997701213661747332399327008,
       -91509.3392102130338542605593697286718077,
       69399.8582111570893316100585838633124312,
       -779.438309639349328345148153897689081893,
       -205.8341686964167118696204191085878165880,
       -398.7823950071290897160364203878571043995,
       354.578231152433375494079868740183658991,
       -786.224179015513894176220583239056456901},
      {5.454547288952965694339504452480078562780,
       -79.78911199784015209705095616004766020335,
       61.0967097444217359754873031115590556707,
       14951.54365344033382142012769129774268946,
       -22324.57139433374168317029445568645401598,
       28594.46085938937782634638310955782423389,
       -21748.11815446623273761450332307272543593,
       245.4393970278627292916961100938952065362,
       65.44129872356201885836080588282812631205,
       104.0129692060648441002024406476025340187,
       -114.7001840640649599911246871588418008302,
       239.7294100413035911863764570341369884827}};

  static double acceptedErrorPower(double error) noexcept {
    return std::pow(error, 1.0 / 8.0);
  }
  static double rejectedErrorPower(double error) noexcept {
    return std::pow(error, -1.0 / 8.0);
  }

  inline static constexpr WVAdaptiveRKStageBufferLastUse stageBufferLastUse[] = {
      {"stage", "stage-state construction", "accepted-state commit or continuous-extension scratch", 17},
      {"k1", "initial right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k2/k3/k5", "stage right-hand side with lifecycle reuse", "stage-13 construction or accepted initial state", 17},
      {"k4/k13", "stage right-hand side with lifecycle reuse", "embedded error estimate", 13},
      {"k6", "sixth-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k7", "seventh-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k8", "eighth-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k9", "ninth-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k10", "tenth-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k11", "eleventh-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k12", "twelfth-stage right-hand side", "accepted solution, error estimate, or continuous extension", 17},
      {"k14", "lazy continuous-extension right-hand side", "continuous-extension polynomial", 17},
      {"k15", "lazy continuous-extension right-hand side", "continuous-extension polynomial", 17},
      {"k16", "lazy continuous-extension right-hand side", "continuous-extension polynomial", 17},
      {"k17", "lazy continuous-extension right-hand side", "continuous-extension polynomial", 17}};
};

// Shared adaptive machinery is intentionally compile-time composed with each
// concrete method. Method workspaces supply explicit stage formulas and error
// increments; this driver supplies the method-neutral tolerance, controller,
// retry bookkeeping, and interval behavior without storing state-sized data.
class AdaptiveRungeKuttaDriver {
public:
  template <typename Options>
  static WVKernelStatus validateOptions(const Options &options,
                                        const char *methodName) {
    if (!(options.relativeTolerance > 0.0) ||
        !(options.absoluteToleranceScale > 0.0) ||
        !(options.safetyFactor > 0.0 && options.safetyFactor <= 1.0) ||
        !(options.rejectionFloorFactor > 0.0 &&
          options.rejectionFloorFactor < 1.0) ||
        !(options.maximumStepFactor >= 1.0) ||
        !(options.maximumStepSize > 0.0))
      return {WVKernelStatusCode::invalidConfiguration,
              std::string(methodName) +
                  " integration tolerances and controller limits must be positive."};
    return WVKernelStatus::ok();
  }

  template <typename Options>
  static WVKernelStatus initializeErrorPolicy(
      WVIntegrationSystem &system, const Options &options,
      std::unique_ptr<WVIntegrationErrorPolicy> &errorPolicy,
      std::uint64_t &toleranceHash,
      std::vector<std::uint64_t> &toleranceComponentHashes) {
    auto status = system.createErrorPolicy(options.absoluteToleranceScale,
                                           errorPolicy);
    if (!status)
      return status;
    const auto &layout = system.stateLayout();
    const auto coefficientFamilyCount = layout.coefficientFamilyCount();
    if (!errorPolicy ||
        errorPolicy->componentCount() !=
            coefficientFamilyCount + layout.additionalBlocks().size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Adaptive error policy does not match the integration layout."};
    for (std::size_t component = 0; component < coefficientFamilyCount;
         ++component)
      if (errorPolicy->elementCount(component) !=
          layout.coefficientFamilies()[component].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Adaptive coefficient tolerance shape does not match the integration layout."};
    for (std::size_t block = 0; block < layout.additionalBlocks().size();
         ++block)
      if (errorPolicy->elementCount(coefficientFamilyCount + block) !=
          layout.additionalBlocks()[block].elementCount)
        return {WVKernelStatusCode::invalidShape,
                "Adaptive state-block tolerance shape does not match the integration layout."};
    toleranceHash = UINT64_C(1469598103934665603);
    toleranceComponentHashes.assign(errorPolicy->componentCount(),
                                    UINT64_C(1469598103934665603));
    for (std::size_t component = 0;
         component < errorPolicy->componentCount(); ++component)
      for (std::size_t index = 0; index < errorPolicy->elementCount(component);
           ++index) {
        const auto tolerance = errorPolicy->absoluteTolerance(component, index);
        toleranceHash = hashTolerance(toleranceHash, tolerance);
        toleranceComponentHashes[component] =
            hashTolerance(toleranceComponentHashes[component], tolerance);
      }
    return WVKernelStatus::ok();
  }

  template <typename ComplexError, typename RealError>
  static double normalizedError(
      const WVIntegrationErrorPolicy &errorPolicy,
      const WVIntegrationStateLayout &layout,
      const IntegrationBuffer &candidate, const WVIntegrationState &initial,
      double relativeTolerance, ComplexError complexError,
      RealError realError, bool includeCandidateInScale = true) noexcept {
    double error = 0.0;
    const auto &complexCandidate = candidate.complex();
    const auto accumulateComplex = [&](std::size_t flatIndex, double absTol,
                                       WVComplex64 initialValue) {
      const auto increment = complexError(flatIndex);
      const auto initialMagnitude =
          std::hypot(initialValue.real, initialValue.imag);
      const auto candidateMagnitude = includeCandidateInScale
                                          ? std::hypot(
                                                complexCandidate[flatIndex].real,
                                                complexCandidate[flatIndex].imag)
                                          : initialMagnitude;
      const auto valueScale =
          std::max(absTol, relativeTolerance *
                               std::max(initialMagnitude, candidateMagnitude));
      const auto ratio =
          std::hypot(increment.real, increment.imag) / valueScale;
      if (!std::isfinite(ratio))
        return false;
      error = std::max(error, ratio);
      return true;
    };
    bool complexErrorFinite = true;
    for (std::size_t family = 0;
         family < layout.coefficientFamilyCount() && complexErrorFinite;
         ++family) {
      const auto &span = layout.coefficientFamilies()[family];
      const auto initialFamily = coefficientFamilyView(layout, initial, family);
      for (std::size_t index = 0; index < span.elementCount; ++index) {
        const auto flatIndex = span.scalarOffset + index;
        complexErrorFinite = accumulateComplex(
            flatIndex, errorPolicy.absoluteTolerance(family, index),
            initialFamily.data[index]);
        if (!complexErrorFinite)
          break;
      }
    }
    const auto coefficientValues = candidate.coefficientCount();
    for (std::size_t blockIndex = 0;
         blockIndex < layout.additionalBlocks().size() && complexErrorFinite;
         ++blockIndex) {
      const auto &span = layout.additionalBlocks()[blockIndex];
      if (span.scalarType != WVStateScalarType::complex64)
        continue;
      for (std::size_t index = 0; index < span.elementCount; ++index) {
        const auto flatIndex =
            coefficientValues + span.scalarOffset + index;
        complexErrorFinite = accumulateComplex(
            flatIndex,
            errorPolicy.absoluteTolerance(layout.coefficientFamilyCount() +
                                              blockIndex,
                                          index),
            initial.additionalBlocks[blockIndex].complexData[index]);
        if (!complexErrorFinite)
          break;
      }
    }
    if (!complexErrorFinite)
      return std::numeric_limits<double>::infinity();
    const auto &realCandidate = candidate.real();
    for (std::size_t blockIndex = 0;
         blockIndex < layout.additionalBlocks().size(); ++blockIndex) {
      const auto &span = layout.additionalBlocks()[blockIndex];
      if (span.scalarType != WVStateScalarType::real64)
        continue;
      for (std::size_t index = 0; index < span.elementCount; ++index) {
        const auto flatIndex = span.scalarOffset + index;
        const auto initialMagnitude = std::abs(
            initial.additionalBlocks[blockIndex].realData[index]);
        const auto candidateMagnitude = includeCandidateInScale
                                            ? std::abs(realCandidate[flatIndex])
                                            : initialMagnitude;
        const auto scale = std::max(
            errorPolicy.absoluteTolerance(layout.coefficientFamilyCount() +
                                              blockIndex,
                                          index),
            relativeTolerance *
                std::max(initialMagnitude, candidateMagnitude));
        const auto ratio = std::abs(realError(flatIndex)) / scale;
        if (!std::isfinite(ratio))
          return std::numeric_limits<double>::infinity();
        error = std::max(error, ratio);
      }
    }
    return error;
  }

  template <typename ControllerPolicy>
  static double controllerFactor(double error,
                                 std::size_t rejectedAttemptCount,
                                 double safetyFactor,
                                 double rejectionFloorFactor,
                                 double repeatedRejectionFactor,
                                 double maximumStepFactor) noexcept {
    if (std::isfinite(error) && error <= 1.0) {
      if (rejectedAttemptCount != 0)
        return 1.0;
      const auto temp =
          1.25 * ControllerPolicy::acceptedErrorPower(error);
      return temp > 0.2 ? 1.0 / temp : maximumStepFactor;
    }
    if (rejectedAttemptCount == 0)
      return std::max(rejectionFloorFactor,
                      safetyFactor *
                          ControllerPolicy::rejectedErrorPower(error));
    return repeatedRejectionFactor;
  }

  template <typename Diagnostic, typename EnsureWorkspace>
  static WVKernelStatus prepareStateAfterRestart(
      WVIntegrationSystem &system, WVMutableIntegrationState &state,
      bool &hasAcceptedStep, bool &fsalAvailable, double &nextStepSize,
      std::vector<Diagnostic> &diagnostics, WVIntegratorMetrics &metrics,
      EnsureWorkspace ensureWorkspace) {
    hasAcceptedStep = false;
    fsalAvailable = false;
    nextStepSize = 0.0;
    diagnostics.clear();
    metrics.diagnosticCapacityBytes =
        diagnostics.capacity() * sizeof(Diagnostic);
    auto status = ensureWorkspace();
    if (!status)
      return status;
    const auto result = system.enforceStateConstraints(state);
    metrics.constraintModifiedCoefficientCount += result.modifiedCoefficientCount;
    return result.status;
  }

  static void recordRejectedAttempt(WVIntegratorMetrics &metrics,
                                    std::size_t &rejectedAttemptCount,
                                    double nextStepSize,
                                    double &stepSize) noexcept {
    ++rejectedAttemptCount;
    ++metrics.rejectedStepCount;
    ++metrics.rejectedInitialDerivativeReuseCount;
    stepSize = nextStepSize;
  }

  static double minimumStepSize(double time) noexcept {
    const auto magnitude = std::abs(time);
    const auto spacing =
        std::nextafter(magnitude, std::numeric_limits<double>::infinity()) -
        magnitude;
    return 16.0 * spacing;
  }

  template <typename Step, typename NextStepSize>
  static WVKernelStatus advanceToTime(WVMutableIntegrationState &state,
                                      double finalTime, double stepSize,
                                      bool stretchFinalStep,
                                      const char *methodName, Step step,
                                      NextStepSize nextStepSize) {
    if (finalTime < state.waveVortex.t || !std::isfinite(finalTime))
      return {WVKernelStatusCode::invalidConfiguration,
              std::string(methodName) +
                  " cannot advance backward or to a nonfinite time."};
    while (state.waveVortex.t < finalTime) {
      const auto remaining = finalTime - state.waveVortex.t;
      if (remaining <= timeTolerance(state.waveVortex.t, finalTime)) {
        state.waveVortex.t = finalTime;
        break;
      }
      const auto use = stretchFinalStep && 1.1 * stepSize >= remaining
                           ? remaining
                           : std::min(stepSize, remaining);
      const auto status = step(use);
      if (!status)
        return status;
      stepSize = nextStepSize();
    }
    state.waveVortex.t = finalTime;
    return WVKernelStatus::ok();
  }
};

} // namespace

class WVFixedStepRK4::Workspace {
public:
  IntegrationBuffer stage, derivative, weighted, initialDerivative;
  std::vector<WVCoefficientFamilyConstView> acceptedCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> acceptedBlockViews;
  std::size_t capacityBytes() const noexcept {
    return stage.capacityBytes() + derivative.capacityBytes() +
           weighted.capacityBytes() + initialDerivative.capacityBytes();
  }
  std::size_t coefficientViewBytes() const noexcept {
    return stage.coefficientViewBytes() + derivative.coefficientViewBytes() +
           weighted.coefficientViewBytes() +
           initialDerivative.coefficientViewBytes();
  }
};

WVFixedStepRK4::WVFixedStepRK4(
    WVIntegrationSystem &system, WVFixedStepRK4Options options)
    : system_(system), options_(options) {}
WVFixedStepRK4::~WVFixedStepRK4() { delete workspace_; }

std::size_t WVFixedStepRK4::persistentBytes() const noexcept {
  return sizeof(*this) +
         (workspace_ == nullptr
              ? 0
              : sizeof(Workspace) + workspace_->capacityBytes() +
                    workspace_->coefficientViewBytes() +
                    workspace_->acceptedCoefficientViews.capacity() *
                        sizeof(WVCoefficientFamilyConstView) +
                    workspace_->acceptedBlockViews.capacity() *
                        sizeof(WVAdditionalStateBlockConstView));
}

WVKernelStatus
WVFixedStepRK4::ensureWorkspace(const WVMutableIntegrationState &state) {
  auto status = validateMutableIntegrationState(system_.stateLayout(), state);
  if (!status)
    return status;
  if (workspace_ != nullptr)
    return WVKernelStatus::ok();
  try {
    workspace_ = new Workspace;
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "RK4 workspace allocation failed."};
  }
  IntegrationBuffer *buffers[] = {&workspace_->stage, &workspace_->derivative,
                                &workspace_->weighted};
  for (auto *buffer : buffers) {
    status = buffer->initialize(system_.stateLayout());
    if (!status) {
      delete workspace_;
      workspace_ = nullptr;
      return status;
    }
  }
  if (options_.retainDenseOutput) {
    status = workspace_->initialDerivative.initialize(system_.stateLayout());
    if (!status) {
      delete workspace_;
      workspace_ = nullptr;
      return status;
    }
  }
  metrics_.workspaceCapacityBytes = workspace_->capacityBytes();
  metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.stateCapacityBytes = workspace_->stage.valueCapacityBytes();
  metrics_.workspaceStateEquivalentCount = options_.retainDenseOutput ? 4 : 3;
  metrics_.workspaceMaximumLiveStateEquivalentCount =
      metrics_.workspaceStateEquivalentCount;
  metrics_.denseHistoryCapacityBytes = options_.retainDenseOutput
                                           ? workspace_->initialDerivative.capacityBytes()
                                           : 0;
  metrics_.denseHistoryStateEquivalentCount =
      options_.retainDenseOutput ? 1 : 0;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
  return WVKernelStatus::ok();
}

WVKernelStatus WVFixedStepRK4::prepareStateAfterRestart(
    WVMutableIntegrationState &state) {
  hasAcceptedStep_ = false;
  acceptedStateConstrained_ = false;
  nextStepSize_ = 0.0;
  auto status = ensureWorkspace(state);
  if (!status)
    return status;
  const auto result = system_.enforceStateConstraints(state);
  acceptedStateConstrained_ = static_cast<bool>(result);
  metrics_.constraintModifiedCoefficientCount += result.modifiedCoefficientCount;
  return result.status;
}

WVKernelStatus WVFixedStepRK4::step(WVMutableIntegrationState &state,
                                             double h) {
  if (stepping_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK4 stepping is not reentrant."};
  if (!(h > 0.0) || !std::isfinite(h))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK4 step size must be finite and positive."};
  auto status = ensureWorkspace(state);
  if (!status)
    return status;
  stepping_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{stepping_};
  hasAcceptedStep_ = false;
  std::vector<WVCoefficientFamilyConstView> coefficientViews;
  std::vector<WVAdditionalStateBlockConstView> stateViews;
  const auto baseView =
      integrationConstView(state, coefficientViews, stateViews);
  const double initialTime = state.waveVortex.t;
  if (acceptedStateConstrained_) {
    auto derivative = workspace_->derivative.flux();
    status = system_.evaluateRightHandSide(baseView, derivative);
  }
  else {
    workspace_->stage.copyFrom(baseView);
    auto initial = workspace_->stage.mutableState(initialTime, state.waveVortex.t0);
    const auto constraint = system_.enforceStateConstraints(initial);
    metrics_.constraintModifiedCoefficientCount += constraint.modifiedCoefficientCount;
    if (!constraint)
      return constraint.status;
    auto derivative = workspace_->derivative.flux();
    status = system_.evaluateRightHandSide(
        workspace_->stage.state(initialTime, state.waveVortex.t0), derivative);
  }
  if (!status)
    return status;
  ++metrics_.rightHandSideEvaluationCount;
  workspace_->weighted.assign(workspace_->derivative);
  metrics_.weightedFluxInitializationElementReads +=
      workspace_->derivative.complex().size();
  metrics_.weightedFluxInitializationElementWrites +=
      workspace_->weighted.complex().size();
  if (options_.retainDenseOutput)
    workspace_->initialDerivative.assign(workspace_->derivative);
  workspace_->stage.setAffine(baseView, workspace_->derivative, 0.5 * h);
  metrics_.stageStateConstructionElementReads +=
      2 * workspace_->stage.complex().size();
  metrics_.stageStateConstructionElementWrites += workspace_->stage.complex().size();
  status = constrain(system_, workspace_->stage, initialTime + 0.5 * h,
                     state.waveVortex.t0);
  if (!status)
    return status;
  status = evaluate(system_, workspace_->stage, initialTime + 0.5 * h,
                    state.waveVortex.t0, workspace_->derivative, metrics_);
  if (!status)
    return status;
  workspace_->weighted.addScaled(workspace_->derivative, 2.0);
  metrics_.weightedAccumulationElementReads +=
      2 * workspace_->derivative.complex().size();
  metrics_.weightedAccumulationElementWrites +=
      workspace_->weighted.complex().size();
  workspace_->stage.setAffine(baseView, workspace_->derivative, 0.5 * h);
  metrics_.stageStateConstructionElementReads +=
      2 * workspace_->stage.complex().size();
  metrics_.stageStateConstructionElementWrites += workspace_->stage.complex().size();
  status = constrain(system_, workspace_->stage, initialTime + 0.5 * h,
                     state.waveVortex.t0);
  if (!status)
    return status;
  status = evaluate(system_, workspace_->stage, initialTime + 0.5 * h,
                    state.waveVortex.t0, workspace_->derivative, metrics_);
  if (!status)
    return status;
  workspace_->weighted.addScaled(workspace_->derivative, 2.0);
  metrics_.weightedAccumulationElementReads +=
      2 * workspace_->derivative.complex().size();
  metrics_.weightedAccumulationElementWrites +=
      workspace_->weighted.complex().size();
  workspace_->stage.setAffine(baseView, workspace_->derivative, h);
  metrics_.stageStateConstructionElementReads +=
      2 * workspace_->stage.complex().size();
  metrics_.stageStateConstructionElementWrites += workspace_->stage.complex().size();
  status = constrain(system_, workspace_->stage, initialTime + h,
                     state.waveVortex.t0);
  if (!status)
    return status;
  status = evaluate(system_, workspace_->stage, initialTime + h,
                    state.waveVortex.t0, workspace_->derivative, metrics_);
  if (!status)
    return status;
  workspace_->weighted.addScaled(workspace_->derivative, 1.0);
  metrics_.weightedAccumulationElementReads +=
      2 * workspace_->derivative.complex().size();
  metrics_.weightedAccumulationElementWrites +=
      workspace_->weighted.complex().size();
  workspace_->stage.setAffine(baseView, workspace_->weighted, h / 6.0);
  metrics_.finalStateUpdateElementReads +=
      2 * workspace_->weighted.complex().size();
  metrics_.finalStateUpdateElementWrites += workspace_->stage.complex().size();
  status = constrain(system_, workspace_->stage, initialTime + h,
                     state.waveVortex.t0);
  if (!status)
    return status;
  if (options_.retainDenseOutput)
    workspace_->weighted.copyFrom(baseView);
  workspace_->stage.copyTo(state);
  metrics_.acceptedStateCommitElementReads += workspace_->stage.complex().size();
  metrics_.acceptedStateCommitElementWrites += workspace_->stage.complex().size();
  state.waveVortex.t = initialTime + h;
  acceptedStateConstrained_ = true;
  makeExternalViews(state, workspace_->acceptedCoefficientViews,
                    workspace_->acceptedBlockViews);
  acceptedStep_ = {
      initialTime,
      state.waveVortex.t,
      {state.waveVortex.view(), workspace_->acceptedBlockViews.data(),
       workspace_->acceptedBlockViews.size(),
       workspace_->acceptedCoefficientViews.data(),
       workspace_->acceptedCoefficientViews.size()},
      {metrics_.acceptedStepCount + 1, 0, 4U,
       h, h, h, 0.0},
      options_.retainDenseOutput ? this : nullptr};
  nextStepSize_ = h;
  metrics_.lastStepSize = h;
  metrics_.nextStepSize = h;
  ++metrics_.acceptedStepCount;
  ++metrics_.stepCount;
  hasAcceptedStep_ = true;
  return WVKernelStatus::ok();
}

WVKernelStatus
WVFixedStepRK4::advanceToTime(WVMutableIntegrationState &state,
                                       double finalTime, double h) {
  if (finalTime < state.waveVortex.t || !std::isfinite(finalTime))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK4 cannot advance backward or to a nonfinite time."};
  while (state.waveVortex.t < finalTime) {
    const auto stepSize = std::min(h, finalTime - state.waveVortex.t);
    if (!(stepSize > 0.0))
      break;
    const auto status = step(state, stepSize);
    if (!status)
      return status;
  }
  state.waveVortex.t = finalTime;
  return WVKernelStatus::ok();
}

WVKernelStatus WVFixedStepRK4::evaluateDenseOutput(
    double time, WVMutableIntegrationState &output) const {
  if (!options_.retainDenseOutput || !hasAcceptedStep_)
    return {WVKernelStatusCode::unsupportedOperation,
            "RK4 dense output is unavailable."};
  if (evaluatingDenseOutput_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK4 dense-output evaluation is not reentrant."};
  auto status = validateMutableIntegrationState(system_.stateLayout(), output);
  if (!status)
    return status;
  const auto tolerance =
      timeTolerance(acceptedStep_.initialTime, acceptedStep_.finalTime);
  if (time < acceptedStep_.initialTime - tolerance ||
      time > acceptedStep_.finalTime + tolerance)
    return {WVKernelStatusCode::invalidConfiguration,
            "Dense-output time is outside the accepted interval."};
  evaluatingDenseOutput_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{evaluatingDenseOutput_};
  const double h = acceptedStep_.finalTime - acceptedStep_.initialTime;
  double theta = h == 0.0 ? 0.0 : (time - acceptedStep_.initialTime) / h;
  theta = std::max(0.0, std::min(1.0, theta));
  const double theta2 = theta * theta, theta3 = theta2 * theta,
               ew = 3 * theta2 - 2 * theta3, iw = 1 - ew,
               isw = h * (theta - 2 * theta2 + theta3),
               esw = h * (theta3 - theta2);
  const auto started = std::chrono::steady_clock::now();
  workspace_->stage.copyFrom(acceptedStep_.endpoint);
  auto &c = workspace_->stage.complex();
  for (std::size_t i = 0; i < c.size(); ++i) {
    const auto a = workspace_->weighted.complex()[i], b = c[i],
               k1 = workspace_->initialDerivative.complex()[i],
               k4 = workspace_->derivative.complex()[i];
    c[i] = {iw * a.real + ew * b.real + isw * k1.real + esw * k4.real,
            iw * a.imag + ew * b.imag + isw * k1.imag + esw * k4.imag};
  }
  auto &r = workspace_->stage.real();
  for (std::size_t i = 0; i < r.size(); ++i) {
    const auto a = workspace_->weighted.real()[i], b = r[i];
    r[i] = iw * a + ew * b + isw * workspace_->initialDerivative.real()[i] +
           esw * workspace_->derivative.real()[i];
  }
  auto mutableStage = workspace_->stage.mutableState(
      time, acceptedStep_.endpoint.waveVortex.t0);
  status = system_.enforceStateConstraints(mutableStage).status;
  if (!status)
    return status;
  workspace_->stage.copyTo(output);
  output.waveVortex.t = time;
  output.waveVortex.t0 = acceptedStep_.endpoint.waveVortex.t0;
  ++metrics_.denseOutputEvaluationCount;
  metrics_.denseOutputElementReads += 4 * c.size();
  metrics_.denseOutputElementWrites += c.size();
  metrics_.denseOutputSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return WVKernelStatus::ok();
}

class WVAdaptiveRK23::Workspace {
public:
  IntegrationBuffer stage, k1, k2, k3, k4;
  std::vector<WVCoefficientFamilyConstView> baseCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> baseBlockViews;
  std::vector<WVCoefficientFamilyConstView> acceptedCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> acceptedBlockViews;
  std::size_t capacityBytes() const noexcept {
    return stage.capacityBytes() + k1.capacityBytes() + k2.capacityBytes() +
           k3.capacityBytes() + k4.capacityBytes();
  }
  std::size_t coefficientViewBytes() const noexcept {
    return stage.coefficientViewBytes() + k1.coefficientViewBytes() +
           k2.coefficientViewBytes() + k3.coefficientViewBytes() +
           k4.coefficientViewBytes();
  }
  std::size_t externalViewBytes() const noexcept {
    return (baseCoefficientViews.capacity() +
            acceptedCoefficientViews.capacity()) *
               sizeof(WVCoefficientFamilyConstView) +
           (baseBlockViews.capacity() + acceptedBlockViews.capacity()) *
               sizeof(WVAdditionalStateBlockConstView);
  }
};

WVAdaptiveRK23::WVAdaptiveRK23(
    WVIntegrationSystem &system,
    WVAdaptiveRK23Options options)
    : system_(system), options_(options) {}
WVAdaptiveRK23::~WVAdaptiveRK23() { delete workspace_; }

std::size_t WVAdaptiveRK23::persistentBytes() const noexcept {
  return sizeof(*this) +
         (workspace_ == nullptr
              ? 0
              : sizeof(Workspace) + workspace_->capacityBytes() +
                    workspace_->coefficientViewBytes() +
                    workspace_->externalViewBytes()) +
         (errorPolicy_ == nullptr ? 0 : errorPolicy_->persistentBytes()) +
         stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK23StepDiagnostic) +
         toleranceComponentHashes_.capacity() * sizeof(std::uint64_t);
}

WVKernelStatus
WVAdaptiveRK23::ensureWorkspace(const WVMutableIntegrationState &state) {
  auto status = validateMutableIntegrationState(system_.stateLayout(), state);
  if (!status)
    return status;
  status = AdaptiveRungeKuttaDriver::validateOptions(options_, "RK23");
  if (!status)
    return status;
  if (workspace_)
    return WVKernelStatus::ok();
  status = AdaptiveRungeKuttaDriver::initializeErrorPolicy(
      system_, options_, errorPolicy_, toleranceHash_,
      toleranceComponentHashes_);
  if (!status)
    return status;
  try {
    workspace_ = new Workspace;
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "RK23 workspace allocation failed."};
  }
  IntegrationBuffer *buffers[] = {&workspace_->stage, &workspace_->k1,
                                &workspace_->k2,    &workspace_->k3,
                                &workspace_->k4};
  for (auto *b : buffers) {
    status = b->initialize(system_.stateLayout());
    if (!status) {
      delete workspace_;
      workspace_ = nullptr;
      return status;
    }
  }
  try {
    const auto coefficientFamilyCount =
        system_.stateLayout().coefficientFamilyCount();
    const auto blockCount = system_.stateLayout().additionalBlocks().size();
    workspace_->baseCoefficientViews.reserve(coefficientFamilyCount);
    workspace_->baseBlockViews.reserve(blockCount);
    workspace_->acceptedCoefficientViews.reserve(coefficientFamilyCount);
    workspace_->acceptedBlockViews.reserve(blockCount);
  } catch (const std::bad_alloc &) {
    delete workspace_;
    workspace_ = nullptr;
    return {WVKernelStatusCode::allocationFailure,
            "RK23 state-view workspace allocation failed."};
  }
  metrics_.workspaceCapacityBytes =
      workspace_->capacityBytes() + workspace_->coefficientViewBytes() +
      workspace_->externalViewBytes();
  metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.stateCapacityBytes = workspace_->stage.valueCapacityBytes();
  metrics_.workspaceStateEquivalentCount = 5;
  metrics_.workspaceMaximumLiveStateEquivalentCount = 5;
  metrics_.denseHistoryCapacityBytes =
      options_.retainDenseOutput
          ? workspace_->k1.valueCapacityBytes() +
                workspace_->k2.valueCapacityBytes() +
                workspace_->k3.valueCapacityBytes() +
                workspace_->k4.valueCapacityBytes()
          : 0;
  metrics_.denseHistoryStateEquivalentCount =
      options_.retainDenseOutput ? 4 : 0;
  metrics_.errorPolicyBytes = errorPolicy_->persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVAdaptiveRK23::prepareStateAfterRestart(
    WVMutableIntegrationState &state) {
  return AdaptiveRungeKuttaDriver::prepareStateAfterRestart(
      system_, state, hasAcceptedStep_, fsalAvailable_, nextStepSize_,
      stepDiagnostics_, metrics_, [&]() { return ensureWorkspace(state); });
}

WVKernelStatus WVAdaptiveRK23::step(WVMutableIntegrationState &state,
                                    double proposedStepSize) {
  if (stepping_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK23 stepping is not reentrant."};
  if (!(proposedStepSize > 0.0) || !std::isfinite(proposedStepSize))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK23 step size must be finite and positive."};
  auto status = ensureWorkspace(state);
  if (!status)
    return status;
  stepping_ = true;
  struct Guard {
    bool &v;
    ~Guard() { v = false; }
  } guard{stepping_};
  hasAcceptedStep_ = false;
  const auto baseView = integrationConstView(
      state, workspace_->baseCoefficientViews, workspace_->baseBlockViews);
  const auto t = state.waveVortex.t, t0 = state.waveVortex.t0;
  double h = std::min(proposedStepSize, options_.maximumStepSize);
  bool initialDerivativeAvailable = false;
  bool reusedFSALDerivative = false;
  if (fsalAvailable_) {
    std::swap(workspace_->k1, workspace_->k4);
    initialDerivativeAvailable = true;
    reusedFSALDerivative = true;
    fsalAvailable_ = false;
    ++metrics_.fsalReuseCount;
  }
  std::size_t rejectedThisStep = 0;
  std::size_t evaluationsThisStep = 0;
  for (;;) {
    if (!std::isfinite(h) || !(t + h > t)) {
      fsalAvailable_ = false;
      return {WVKernelStatusCode::numericalFailure,
              "RK23 cannot advance time with the proposed step."};
    }
    if (!initialDerivativeAvailable) {
      const auto before = metrics_.rightHandSideEvaluationCount;
      auto derivative = workspace_->k1.flux();
      status = system_.evaluateRightHandSide(baseView, derivative);
      if (status)
        ++metrics_.rightHandSideEvaluationCount;
      evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
      if (!status) {
        fsalAvailable_ = false;
        return status;
      }
      initialDerivativeAvailable = true;
    }
    workspace_->stage.setAffine(baseView, workspace_->k1, 0.5 * h);
    status = constrain(system_, workspace_->stage, t + 0.5 * h, t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    auto before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + 0.5 * h, t0,
                      workspace_->k2, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    workspace_->stage.setAffine(baseView, workspace_->k2, 0.75 * h);
    status = constrain(system_, workspace_->stage, t + 0.75 * h, t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + 0.75 * h, t0,
                      workspace_->k3, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    workspace_->stage.setWeightedCandidate(
        baseView, h, workspace_->k1, 2.0 / 9.0, workspace_->k2, 1.0 / 3.0,
        workspace_->k3, 4.0 / 9.0);
    auto candidateState = workspace_->stage.mutableState(t + h, t0);
    const auto endpointConstraint = system_.enforceStateConstraints(candidateState);
    metrics_.constraintModifiedCoefficientCount +=
        endpointConstraint.modifiedCoefficientCount;
    if (!endpointConstraint) {
      fsalAvailable_ = false;
      return endpointConstraint.status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h, t0, workspace_->k4,
                      metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    const auto complexError = [&](std::size_t i) noexcept {
      return WVComplex64{
          h * ((2.0 / 9.0 - 7.0 / 24.0) * workspace_->k1.complex()[i].real +
               (1.0 / 3.0 - 0.25) * workspace_->k2.complex()[i].real +
               (4.0 / 9.0 - 1.0 / 3.0) * workspace_->k3.complex()[i].real -
               0.125 * workspace_->k4.complex()[i].real),
          h * ((2.0 / 9.0 - 7.0 / 24.0) * workspace_->k1.complex()[i].imag +
               (1.0 / 3.0 - 0.25) * workspace_->k2.complex()[i].imag +
               (4.0 / 9.0 - 1.0 / 3.0) * workspace_->k3.complex()[i].imag -
               0.125 * workspace_->k4.complex()[i].imag)};
    };
    const auto realError = [&](std::size_t i) noexcept {
      return h * ((2.0 / 9.0 - 7.0 / 24.0) * workspace_->k1.real()[i] +
                  (1.0 / 3.0 - 0.25) * workspace_->k2.real()[i] +
                  (4.0 / 9.0 - 1.0 / 3.0) * workspace_->k3.real()[i] -
                  0.125 * workspace_->k4.real()[i]);
    };
    const auto error = AdaptiveRungeKuttaDriver::normalizedError(
        *errorPolicy_, system_.stateLayout(), workspace_->stage, baseView,
        options_.relativeTolerance, complexError, realError);
    const bool accepted = std::isfinite(error) && error <= 1.0;
    const auto factor =
        AdaptiveRungeKuttaDriver::controllerFactor<RK23ControllerPolicy>(
            error, rejectedThisStep, options_.safetyFactor,
            options_.rejectionFloorFactor, options_.rejectionFloorFactor,
            options_.maximumStepFactor);
    nextStepSize_ =
        std::min(options_.maximumStepSize, h * factor);
    metrics_.lastProposedStepSize = proposedStepSize;
    metrics_.normalizedError = error;
    metrics_.nextStepSize = nextStepSize_;
    if (accepted) {
      workspace_->stage.copyTo(state);
      state.waveVortex.t = t + h;
      makeExternalViews(state, workspace_->acceptedCoefficientViews,
                        workspace_->acceptedBlockViews);
      acceptedStep_ = {
          t,
          state.waveVortex.t,
          {state.waveVortex.view(), workspace_->acceptedBlockViews.data(),
           workspace_->acceptedBlockViews.size(),
           workspace_->acceptedCoefficientViews.data(),
           workspace_->acceptedCoefficientViews.size()},
          {metrics_.acceptedStepCount + 1, rejectedThisStep,
           evaluationsThisStep, h, proposedStepSize, nextStepSize_, error},
          options_.retainDenseOutput ? this : nullptr};
      ++metrics_.acceptedStepCount;
      metrics_.lastAcceptedStepSize = h;
      fsalAvailable_ = endpointConstraint.modifiedCoefficientCount == 0 &&
                       endpointConstraint.fsalCompatible;
      if (!fsalAvailable_)
        ++metrics_.fsalInvalidationCount;
      hasAcceptedStep_ = true;
      if (stepDiagnostics_.size() < options_.maximumRecordedStepDiagnostics)
        stepDiagnostics_.push_back(
            {t, h, error, nextStepSize_, rejectedThisStep,
             evaluationsThisStep, reusedFSALDerivative});
      metrics_.diagnosticCapacityBytes =
          stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK23StepDiagnostic);
      return WVKernelStatus::ok();
    }
    AdaptiveRungeKuttaDriver::recordRejectedAttempt(
        metrics_, rejectedThisStep, nextStepSize_, h);
    if (!(h > 0.0) || t + h == t)
      return {WVKernelStatusCode::numericalFailure,
              "RK23 step size underflowed after rejection."};
  }
}

WVKernelStatus
WVAdaptiveRK23::advanceToTime(WVMutableIntegrationState &state,
                                       double finalTime, double h) {
  return AdaptiveRungeKuttaDriver::advanceToTime(
      state, finalTime, h, false, "RK23",
      [&](double use) { return step(state, use); },
      [&]() { return nextStepSize_; });
}

WVKernelStatus WVAdaptiveRK23::evaluateDenseOutput(
    double time, WVMutableIntegrationState &output) const {
  if (!options_.retainDenseOutput || !hasAcceptedStep_)
    return {WVKernelStatusCode::unsupportedOperation,
            "RK23 dense output is unavailable."};
  if (evaluatingDenseOutput_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK23 dense-output evaluation is not reentrant."};
  auto status = validateMutableIntegrationState(system_.stateLayout(), output);
  if (!status)
    return status;
  const double h = acceptedStep_.finalTime - acceptedStep_.initialTime;
  const auto tol =
      timeTolerance(acceptedStep_.initialTime, acceptedStep_.finalTime);
  if (time < acceptedStep_.initialTime - tol ||
      time > acceptedStep_.finalTime + tol)
    return {WVKernelStatusCode::invalidConfiguration,
            "Dense-output time is outside the accepted interval."};
  evaluatingDenseOutput_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{evaluatingDenseOutput_};
  double theta = h == 0 ? 0 : (time - acceptedStep_.initialTime) / h;
  if (std::abs(time - acceptedStep_.initialTime) <= tol)
    theta = 0.0;
  if (std::abs(time - acceptedStep_.finalTime) <= tol)
    theta = 1.0;
  theta = std::max(0.0, std::min(1.0, theta));
  const double t2 = theta * theta, t3 = t2 * theta;
  const double weights[] = {
      theta - (4.0 / 3.0) * t2 + (5.0 / 9.0) * t3 - 2.0 / 9.0,
      t2 - (2.0 / 3.0) * t3 - 1.0 / 3.0,
      (4.0 / 3.0) * t2 - (8.0 / 9.0) * t3 - 4.0 / 9.0,
      -t2 + t3};
  const auto started = std::chrono::steady_clock::now();
  workspace_->stage.copyFrom(acceptedStep_.endpoint);
  const IntegrationBuffer *derivatives[] = {&workspace_->k1, &workspace_->k2,
                                          &workspace_->k3, &workspace_->k4};
  auto &complex = workspace_->stage.complex();
  for (std::size_t i = 0; i < complex.size(); ++i) {
    for (std::size_t derivative = 0; derivative < 4; ++derivative) {
      complex[i].real +=
          h * weights[derivative] * derivatives[derivative]->complex()[i].real;
      complex[i].imag +=
          h * weights[derivative] * derivatives[derivative]->complex()[i].imag;
    }
  }
  auto &real = workspace_->stage.real();
  for (std::size_t i = 0; i < real.size(); ++i) {
    for (std::size_t derivative = 0; derivative < 4; ++derivative)
      real[i] += h * weights[derivative] * derivatives[derivative]->real()[i];
  }
  auto stage = workspace_->stage.mutableState(
      time, acceptedStep_.endpoint.waveVortex.t0);
  status = system_.enforceStateConstraints(stage).status;
  if (!status)
    return status;
  workspace_->stage.copyTo(output);
  output.waveVortex.t = time;
  output.waveVortex.t0 = acceptedStep_.endpoint.waveVortex.t0;
  ++metrics_.denseOutputEvaluationCount;
  metrics_.denseOutputElementReads += 5 * (complex.size() + real.size());
  metrics_.denseOutputElementWrites += complex.size() + real.size();
  metrics_.denseOutputSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return WVKernelStatus::ok();
}

namespace {

constexpr WVAdaptiveRKStageBufferLastUse rk4StageBufferLastUse[] = {
    {"stage", "stage-state construction", "accepted-state commit or dense-output scratch", 4},
    {"derivative", "right-hand-side evaluation", "weighted accumulation", 4},
    {"weighted", "weighted derivative initialization", "final-state update", 4},
    {"initialDerivative", "accepted-step initial derivative", "dense-output interpolation", 4}};

constexpr WVAdaptiveRKStageBufferLastUse rk23StageBufferLastUse[] = {
    {"stage", "stage-state construction", "accepted-state commit or dense-output scratch", 4},
    {"k1", "initial derivative or FSAL reuse", "continuous extension or next-step FSAL swap", 4},
    {"k2", "second-stage right-hand side", "continuous extension", 4},
    {"k3", "third-stage right-hand side", "continuous extension", 4},
    {"k4", "endpoint right-hand side", "continuous extension or next-step FSAL swap", 4}};

constexpr WVAdaptiveRKStageBufferLastUse rk45StageBufferLastUse[] = {
    {"stage", "stage-state construction", "accepted-state commit or dense-output scratch", 7},
    {"k1", "initial derivative or FSAL reuse", "continuous extension or next-step FSAL swap", 7},
    {"k2/k7", "second-stage then endpoint right-hand side", "continuous extension or next-step FSAL swap", 7},
    {"k3", "third-stage right-hand side", "continuous extension", 7},
    {"k4", "fourth-stage right-hand side", "continuous extension", 7},
    {"k5", "fifth-stage right-hand side", "continuous extension", 7},
    {"k6", "sixth-stage right-hand side", "continuous extension", 7}};

} // namespace

const WVAdaptiveRKStageBufferLastUse *
WVFixedStepRK4::stageBufferLastUseRecords() noexcept {
  return rk4StageBufferLastUse;
}

std::size_t WVFixedStepRK4::stageBufferLastUseRecordCount() noexcept {
  return sizeof(rk4StageBufferLastUse) / sizeof(rk4StageBufferLastUse[0]);
}

const WVAdaptiveRKStageBufferLastUse *
WVAdaptiveRK23::stageBufferLastUseRecords() noexcept {
  return rk23StageBufferLastUse;
}

std::size_t WVAdaptiveRK23::stageBufferLastUseRecordCount() noexcept {
  return sizeof(rk23StageBufferLastUse) / sizeof(rk23StageBufferLastUse[0]);
}

class WVAdaptiveRK45::Workspace {
public:
  IntegrationBuffer stage, k1, k2OrK7, k3, k4, k5, k6;
  std::vector<WVCoefficientFamilyConstView> baseCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> baseBlockViews;
  std::vector<WVCoefficientFamilyConstView> acceptedCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> acceptedBlockViews;
  std::size_t capacityBytes() const noexcept {
    return stage.capacityBytes() + k1.capacityBytes() +
           k2OrK7.capacityBytes() + k3.capacityBytes() + k4.capacityBytes() +
           k5.capacityBytes() + k6.capacityBytes();
  }
  std::size_t coefficientViewBytes() const noexcept {
    return stage.coefficientViewBytes() + k1.coefficientViewBytes() +
           k2OrK7.coefficientViewBytes() + k3.coefficientViewBytes() +
           k4.coefficientViewBytes() + k5.coefficientViewBytes() +
           k6.coefficientViewBytes();
  }
  std::size_t externalViewBytes() const noexcept {
    return (baseCoefficientViews.capacity() +
            acceptedCoefficientViews.capacity()) *
               sizeof(WVCoefficientFamilyConstView) +
           (baseBlockViews.capacity() + acceptedBlockViews.capacity()) *
           sizeof(WVAdditionalStateBlockConstView);
  }
};

WVAdaptiveRK45::WVAdaptiveRK45(WVIntegrationSystem &system,
                               WVAdaptiveRK45Options options)
    : system_(system), options_(options) {}

WVAdaptiveRK45::~WVAdaptiveRK45() { delete workspace_; }

std::size_t WVAdaptiveRK45::persistentBytes() const noexcept {
  return sizeof(*this) +
         (workspace_ == nullptr
              ? 0
              : sizeof(Workspace) + workspace_->capacityBytes() +
                    workspace_->coefficientViewBytes() +
                    workspace_->externalViewBytes()) +
         (errorPolicy_ == nullptr ? 0 : errorPolicy_->persistentBytes()) +
         stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK45StepDiagnostic) +
         toleranceComponentHashes_.capacity() * sizeof(std::uint64_t);
}

WVKernelStatus
WVAdaptiveRK45::ensureWorkspace(const WVMutableIntegrationState &state) {
  auto status = validateMutableIntegrationState(system_.stateLayout(), state);
  if (!status)
    return status;
  status = AdaptiveRungeKuttaDriver::validateOptions(options_, "RK45");
  if (!status)
    return status;
  if (!(options_.repeatedRejectionFactor > 0.0 &&
        options_.repeatedRejectionFactor < 1.0))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK45 repeated-rejection factor must be between zero and one."};
  if (workspace_ != nullptr)
    return WVKernelStatus::ok();
  status = AdaptiveRungeKuttaDriver::initializeErrorPolicy(
      system_, options_, errorPolicy_, toleranceHash_,
      toleranceComponentHashes_);
  if (!status)
    return status;
  try {
    workspace_ = new Workspace;
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "RK45 workspace allocation failed."};
  }
  IntegrationBuffer *buffers[] = {
      &workspace_->stage, &workspace_->k1, &workspace_->k2OrK7,
      &workspace_->k3,    &workspace_->k4, &workspace_->k5,
      &workspace_->k6};
  for (auto *buffer : buffers) {
    status = buffer->initialize(system_.stateLayout());
    if (!status) {
      delete workspace_;
      workspace_ = nullptr;
      return status;
    }
  }
  try {
    workspace_->baseCoefficientViews.reserve(
        system_.stateLayout().coefficientFamilyCount());
    workspace_->baseBlockViews.reserve(
        system_.stateLayout().additionalBlocks().size());
    workspace_->acceptedCoefficientViews.reserve(
        system_.stateLayout().coefficientFamilyCount());
    workspace_->acceptedBlockViews.reserve(
        system_.stateLayout().additionalBlocks().size());
  } catch (const std::bad_alloc &) {
    delete workspace_;
    workspace_ = nullptr;
    return {WVKernelStatusCode::allocationFailure,
            "RK45 state-view workspace allocation failed."};
  }
  metrics_.workspaceCapacityBytes =
      workspace_->capacityBytes() + workspace_->coefficientViewBytes() +
      workspace_->externalViewBytes();
  metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.stateCapacityBytes = workspace_->stage.valueCapacityBytes();
  metrics_.workspaceStateEquivalentCount = 7;
  metrics_.workspaceMaximumLiveStateEquivalentCount = 7;
  metrics_.denseHistoryCapacityBytes =
      options_.retainDenseOutput
          ? workspace_->k1.valueCapacityBytes() +
                workspace_->k2OrK7.valueCapacityBytes() +
                workspace_->k3.valueCapacityBytes() +
                workspace_->k4.valueCapacityBytes() +
                workspace_->k5.valueCapacityBytes() +
                workspace_->k6.valueCapacityBytes()
          : 0;
  metrics_.denseHistoryStateEquivalentCount =
      options_.retainDenseOutput ? 6 : 0;
  metrics_.errorPolicyBytes = errorPolicy_->persistentBytes();
  return WVKernelStatus::ok();
}

WVKernelStatus WVAdaptiveRK45::prepareStateAfterRestart(
    WVMutableIntegrationState &state) {
  return AdaptiveRungeKuttaDriver::prepareStateAfterRestart(
      system_, state, hasAcceptedStep_, fsalAvailable_, nextStepSize_,
      stepDiagnostics_, metrics_, [&]() { return ensureWorkspace(state); });
}

WVKernelStatus WVAdaptiveRK45::step(WVMutableIntegrationState &state,
                                    double proposedStepSize) {
  return stepImplementation(state, proposedStepSize, false);
}

WVKernelStatus WVAdaptiveRK45::stepImplementation(
    WVMutableIntegrationState &state, double proposedStepSize,
    bool allowFinalStepStretch) {
  if (stepping_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK45 stepping is not reentrant."};
  if (!(proposedStepSize > 0.0) || !std::isfinite(proposedStepSize))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK45 step size must be finite and positive."};
  auto status = ensureWorkspace(state);
  if (!status)
    return status;
  stepping_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{stepping_};
  hasAcceptedStep_ = false;
  const auto baseView = integrationConstView(
      state, workspace_->baseCoefficientViews, workspace_->baseBlockViews);
  const auto t = state.waveVortex.t;
  const auto t0 = state.waveVortex.t0;
  double h = allowFinalStepStretch
                 ? proposedStepSize
                 : std::min(proposedStepSize, options_.maximumStepSize);
  bool initialDerivativeAvailable = false;
  bool reusedFSALDerivative = false;
  if (fsalAvailable_) {
    std::swap(workspace_->k1, workspace_->k2OrK7);
    initialDerivativeAvailable = true;
    reusedFSALDerivative = true;
    fsalAvailable_ = false;
    ++metrics_.fsalReuseCount;
  }
  std::size_t rejectedThisStep = 0;
  std::size_t evaluationsThisStep = 0;
  for (;;) {
    if (!std::isfinite(h) || !(t + h > t)) {
      fsalAvailable_ = false;
      return {WVKernelStatusCode::numericalFailure,
              "RK45 cannot advance time with the proposed step."};
    }
    if (!initialDerivativeAvailable) {
      const auto before = metrics_.rightHandSideEvaluationCount;
      auto derivative = workspace_->k1.flux();
      status = system_.evaluateRightHandSide(baseView, derivative);
      if (status)
        ++metrics_.rightHandSideEvaluationCount;
      evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
      if (!status) {
        fsalAvailable_ = false;
        return status;
      }
      initialDerivativeAvailable = true;
    }

    workspace_->stage.setAffine(baseView, workspace_->k1, h * (1.0 / 5.0));
    status = constrain(system_, workspace_->stage, t + h * (1.0 / 5.0), t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    auto before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h * (1.0 / 5.0), t0,
                      workspace_->k2OrK7, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    workspace_->stage.copyFrom(baseView);
    workspace_->stage.addScaled(workspace_->k1, h * (3.0 / 40.0));
    workspace_->stage.addScaled(workspace_->k2OrK7, h * (9.0 / 40.0));
    status = constrain(system_, workspace_->stage, t + h * (3.0 / 10.0), t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h * (3.0 / 10.0), t0,
                      workspace_->k3, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    workspace_->stage.copyFrom(baseView);
    workspace_->stage.addScaled(workspace_->k1, h * (44.0 / 45.0));
    workspace_->stage.addScaled(workspace_->k2OrK7, h * (-56.0 / 15.0));
    workspace_->stage.addScaled(workspace_->k3, h * (32.0 / 9.0));
    status = constrain(system_, workspace_->stage, t + h * (4.0 / 5.0), t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h * (4.0 / 5.0), t0,
                      workspace_->k4, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    workspace_->stage.copyFrom(baseView);
    workspace_->stage.addScaled(workspace_->k1, h * (19372.0 / 6561.0));
    workspace_->stage.addScaled(workspace_->k2OrK7,
                                h * (-25360.0 / 2187.0));
    workspace_->stage.addScaled(workspace_->k3, h * (64448.0 / 6561.0));
    workspace_->stage.addScaled(workspace_->k4, h * (-212.0 / 729.0));
    status = constrain(system_, workspace_->stage, t + h * (8.0 / 9.0), t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h * (8.0 / 9.0), t0,
                      workspace_->k5, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    workspace_->stage.copyFrom(baseView);
    workspace_->stage.addScaled(workspace_->k1, h * (9017.0 / 3168.0));
    workspace_->stage.addScaled(workspace_->k2OrK7, h * (-355.0 / 33.0));
    workspace_->stage.addScaled(workspace_->k3, h * (46732.0 / 5247.0));
    workspace_->stage.addScaled(workspace_->k4, h * (49.0 / 176.0));
    workspace_->stage.addScaled(workspace_->k5, h * (-5103.0 / 18656.0));
    status = constrain(system_, workspace_->stage, t + h, t0);
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h, t0, workspace_->k6,
                      metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    workspace_->stage.copyFrom(baseView);
    workspace_->stage.addScaled(workspace_->k1, h * (35.0 / 384.0));
    workspace_->stage.addScaled(workspace_->k3, h * (500.0 / 1113.0));
    workspace_->stage.addScaled(workspace_->k4, h * (125.0 / 192.0));
    workspace_->stage.addScaled(workspace_->k5, h * (-2187.0 / 6784.0));
    workspace_->stage.addScaled(workspace_->k6, h * (11.0 / 84.0));
    auto candidateState = workspace_->stage.mutableState(t + h, t0);
    const auto endpointConstraint =
        system_.enforceStateConstraints(candidateState);
    metrics_.constraintModifiedCoefficientCount +=
        endpointConstraint.modifiedCoefficientCount;
    if (!endpointConstraint) {
      fsalAvailable_ = false;
      return endpointConstraint.status;
    }
    before = metrics_.rightHandSideEvaluationCount;
    status = evaluate(system_, workspace_->stage, t + h, t0,
                      workspace_->k2OrK7, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    if (!status) {
      fsalAvailable_ = false;
      return status;
    }

    const auto complexError = [&](std::size_t i) noexcept {
      return WVComplex64{
          h * ((71.0 / 57600.0) * workspace_->k1.complex()[i].real -
               (71.0 / 16695.0) * workspace_->k3.complex()[i].real +
               (71.0 / 1920.0) * workspace_->k4.complex()[i].real -
               (17253.0 / 339200.0) * workspace_->k5.complex()[i].real +
               (22.0 / 525.0) * workspace_->k6.complex()[i].real -
               (1.0 / 40.0) * workspace_->k2OrK7.complex()[i].real),
          h * ((71.0 / 57600.0) * workspace_->k1.complex()[i].imag -
               (71.0 / 16695.0) * workspace_->k3.complex()[i].imag +
               (71.0 / 1920.0) * workspace_->k4.complex()[i].imag -
               (17253.0 / 339200.0) * workspace_->k5.complex()[i].imag +
               (22.0 / 525.0) * workspace_->k6.complex()[i].imag -
               (1.0 / 40.0) * workspace_->k2OrK7.complex()[i].imag)};
    };
    const auto realError = [&](std::size_t i) noexcept {
      return h * ((71.0 / 57600.0) * workspace_->k1.real()[i] -
                  (71.0 / 16695.0) * workspace_->k3.real()[i] +
                  (71.0 / 1920.0) * workspace_->k4.real()[i] -
                  (17253.0 / 339200.0) * workspace_->k5.real()[i] +
                  (22.0 / 525.0) * workspace_->k6.real()[i] -
                  (1.0 / 40.0) * workspace_->k2OrK7.real()[i]);
    };
    const auto error = AdaptiveRungeKuttaDriver::normalizedError(
        *errorPolicy_, system_.stateLayout(), workspace_->stage, baseView,
        options_.relativeTolerance, complexError, realError);
    const auto accepted = std::isfinite(error) && error <= 1.0;
    const auto factor =
        AdaptiveRungeKuttaDriver::controllerFactor<RK45ControllerPolicy>(
            error, rejectedThisStep, options_.safetyFactor,
            options_.rejectionFloorFactor, options_.repeatedRejectionFactor,
            options_.maximumStepFactor);
    nextStepSize_ = std::min(options_.maximumStepSize, h * factor);
    metrics_.lastProposedStepSize = proposedStepSize;
    metrics_.normalizedError = error;
    metrics_.nextStepSize = nextStepSize_;
    if (accepted) {
      workspace_->stage.copyTo(state);
      state.waveVortex.t = t + h;
      makeExternalViews(state, workspace_->acceptedCoefficientViews,
                        workspace_->acceptedBlockViews);
      acceptedStep_ = {
          t,
          state.waveVortex.t,
          {state.waveVortex.view(), workspace_->acceptedBlockViews.data(),
           workspace_->acceptedBlockViews.size(),
           workspace_->acceptedCoefficientViews.data(),
           workspace_->acceptedCoefficientViews.size()},
          {metrics_.acceptedStepCount + 1, rejectedThisStep,
           evaluationsThisStep, h, proposedStepSize, nextStepSize_, error},
          options_.retainDenseOutput ? this : nullptr};
      ++metrics_.acceptedStepCount;
      ++metrics_.stepCount;
      metrics_.lastStepSize = h;
      metrics_.lastAcceptedStepSize = h;
      fsalAvailable_ = endpointConstraint.modifiedCoefficientCount == 0 &&
                       endpointConstraint.fsalCompatible;
      if (!fsalAvailable_)
        ++metrics_.fsalInvalidationCount;
      hasAcceptedStep_ = true;
      if (stepDiagnostics_.size() < options_.maximumRecordedStepDiagnostics)
        stepDiagnostics_.push_back(
            {t, h, error, nextStepSize_, rejectedThisStep,
             evaluationsThisStep, reusedFSALDerivative});
      metrics_.diagnosticCapacityBytes =
          stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK45StepDiagnostic);
      return WVKernelStatus::ok();
    }
    AdaptiveRungeKuttaDriver::recordRejectedAttempt(
        metrics_, rejectedThisStep, nextStepSize_, h);
    if (!(h > 0.0) || t + h == t)
      return {WVKernelStatusCode::numericalFailure,
              "RK45 step size underflowed after rejection."};
  }
}

WVKernelStatus WVAdaptiveRK45::advanceToTime(
    WVMutableIntegrationState &state, double finalTime, double h) {
  return AdaptiveRungeKuttaDriver::advanceToTime(
      state, finalTime, h, true, "RK45",
      [&](double use) {
        return stepImplementation(
            state, use,
            use > options_.maximumStepSize &&
                use <= 1.1 * options_.maximumStepSize);
      },
      [&]() { return nextStepSize_; });
}

WVKernelStatus WVAdaptiveRK45::evaluateDenseOutput(
    double time, WVMutableIntegrationState &output) const {
  if (!options_.retainDenseOutput || !hasAcceptedStep_)
    return {WVKernelStatusCode::unsupportedOperation,
            "RK45 dense output is unavailable."};
  if (evaluatingDenseOutput_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK45 dense-output evaluation is not reentrant."};
  auto status = validateMutableIntegrationState(system_.stateLayout(), output);
  if (!status)
    return status;
  const double h = acceptedStep_.finalTime - acceptedStep_.initialTime;
  const auto tolerance =
      timeTolerance(acceptedStep_.initialTime, acceptedStep_.finalTime);
  if (time < acceptedStep_.initialTime - tolerance ||
      time > acceptedStep_.finalTime + tolerance)
    return {WVKernelStatusCode::invalidConfiguration,
            "Dense-output time is outside the accepted interval."};
  evaluatingDenseOutput_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{evaluatingDenseOutput_};
  double theta = h == 0.0 ? 0.0 : (time - acceptedStep_.initialTime) / h;
  if (std::abs(time - acceptedStep_.initialTime) <= tolerance)
    theta = 0.0;
  if (std::abs(time - acceptedStep_.finalTime) <= tolerance)
    theta = 1.0;
  theta = std::max(0.0, std::min(1.0, theta));
  const double theta2 = theta * theta;
  const double weights[] = {
      theta + theta2 * (-183.0 / 64.0 +
                        theta * (37.0 / 12.0 - 145.0 / 128.0 * theta)) -
          35.0 / 384.0,
      theta2 * (1500.0 / 371.0 +
                theta * (-1000.0 / 159.0 + 1000.0 / 371.0 * theta)) -
          500.0 / 1113.0,
      theta2 * (-125.0 / 32.0 +
                theta * (125.0 / 12.0 - 375.0 / 64.0 * theta)) -
          125.0 / 192.0,
      theta2 * (9477.0 / 3392.0 +
                theta * (-729.0 / 106.0 + 25515.0 / 6784.0 * theta)) +
          2187.0 / 6784.0,
      theta2 * (-11.0 / 7.0 +
                theta * (11.0 / 3.0 - 55.0 / 28.0 * theta)) -
          11.0 / 84.0,
      theta2 * (3.0 / 2.0 + theta * (-4.0 + 5.0 / 2.0 * theta))};
  const auto started = std::chrono::steady_clock::now();
  workspace_->stage.copyFrom(acceptedStep_.endpoint);
  const IntegrationBuffer *derivatives[] = {
      &workspace_->k1, &workspace_->k3, &workspace_->k4,
      &workspace_->k5, &workspace_->k6, &workspace_->k2OrK7};
  auto &complex = workspace_->stage.complex();
  for (std::size_t i = 0; i < complex.size(); ++i)
    for (std::size_t derivative = 0; derivative < 6; ++derivative) {
      complex[i].real +=
          h * weights[derivative] * derivatives[derivative]->complex()[i].real;
      complex[i].imag +=
          h * weights[derivative] * derivatives[derivative]->complex()[i].imag;
    }
  auto &real = workspace_->stage.real();
  for (std::size_t i = 0; i < real.size(); ++i)
    for (std::size_t derivative = 0; derivative < 6; ++derivative)
      real[i] += h * weights[derivative] * derivatives[derivative]->real()[i];
  auto stage = workspace_->stage.mutableState(
      time, acceptedStep_.endpoint.waveVortex.t0);
  status = system_.enforceStateConstraints(stage).status;
  if (!status)
    return status;
  workspace_->stage.copyTo(output);
  output.waveVortex.t = time;
  output.waveVortex.t0 = acceptedStep_.endpoint.waveVortex.t0;
  ++metrics_.denseOutputEvaluationCount;
  metrics_.denseOutputElementReads += 7 * (complex.size() + real.size());
  metrics_.denseOutputElementWrites += complex.size() + real.size();
  metrics_.denseOutputSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return WVKernelStatus::ok();
}

double WVAdaptiveRK45::initialTime() const noexcept {
  return hasAcceptedStep_ ? acceptedStep_.initialTime : 0.0;
}

double WVAdaptiveRK45::finalTime() const noexcept {
  return hasAcceptedStep_ ? acceptedStep_.finalTime : 0.0;
}

const WVAcceptedStep *WVAdaptiveRK45::lastAcceptedStep() const noexcept {
  return hasAcceptedStep_ ? &acceptedStep_ : nullptr;
}

const WVIntegratorMetrics &WVAdaptiveRK45::metrics() const noexcept {
  return metrics_;
}

const std::vector<WVAdaptiveRK45StepDiagnostic> &
WVAdaptiveRK45::stepDiagnostics() const noexcept {
  return stepDiagnostics_;
}

std::uint64_t WVAdaptiveRK45::toleranceHash() const noexcept {
  return toleranceHash_;
}

const std::vector<std::uint64_t> &
WVAdaptiveRK45::toleranceComponentHashes() const noexcept {
  return toleranceComponentHashes_;
}

bool WVAdaptiveRK45::stepDiagnosticsComplete() const noexcept {
  return stepDiagnostics_.size() == metrics_.acceptedStepCount;
}

const WVAdaptiveRKStageBufferLastUse *
WVAdaptiveRK45::stageBufferLastUseRecords() noexcept {
  return rk45StageBufferLastUse;
}

std::size_t WVAdaptiveRK45::stageBufferLastUseRecordCount() noexcept {
  return sizeof(rk45StageBufferLastUse) / sizeof(rk45StageBufferLastUse[0]);
}

double WVAdaptiveRK45::nextStepSize() const noexcept { return nextStepSize_; }

class WVAdaptiveRK78::Workspace {
public:
  IntegrationBuffer stage, k1, k2OrK3OrK5, k4OrK13, k6, k7, k8, k9, k10,
      k11, k12;
  std::vector<WVCoefficientFamilyConstView> baseCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> baseBlockViews;
  std::vector<WVCoefficientFamilyConstView> acceptedCoefficientViews;
  std::vector<WVAdditionalStateBlockConstView> acceptedBlockViews;

  std::size_t capacityBytes() const noexcept {
    return stage.capacityBytes() + k1.capacityBytes() +
           k2OrK3OrK5.capacityBytes() + k4OrK13.capacityBytes() +
           k6.capacityBytes() + k7.capacityBytes() + k8.capacityBytes() +
           k9.capacityBytes() + k10.capacityBytes() + k11.capacityBytes() +
           k12.capacityBytes();
  }
  std::size_t coefficientViewBytes() const noexcept {
    return stage.coefficientViewBytes() + k1.coefficientViewBytes() +
           k2OrK3OrK5.coefficientViewBytes() +
           k4OrK13.coefficientViewBytes() + k6.coefficientViewBytes() +
           k7.coefficientViewBytes() + k8.coefficientViewBytes() +
           k9.coefficientViewBytes() + k10.coefficientViewBytes() +
           k11.coefficientViewBytes() + k12.coefficientViewBytes();
  }
  std::size_t externalViewBytes() const noexcept {
    return (baseCoefficientViews.capacity() +
            acceptedCoefficientViews.capacity()) *
               sizeof(WVCoefficientFamilyConstView) +
           (baseBlockViews.capacity() + acceptedBlockViews.capacity()) *
               sizeof(WVAdditionalStateBlockConstView);
  }
};

class WVAdaptiveRK78::ContinuousExtensionWorkspace {
public:
  IntegrationBuffer k14, k15, k16, k17;

  std::size_t capacityBytes() const noexcept {
    return k14.capacityBytes() + k15.capacityBytes() + k16.capacityBytes() +
           k17.capacityBytes();
  }
  std::size_t coefficientViewBytes() const noexcept {
    return k14.coefficientViewBytes() + k15.coefficientViewBytes() +
           k16.coefficientViewBytes() + k17.coefficientViewBytes();
  }
};

WVAdaptiveRK78::WVAdaptiveRK78(WVIntegrationSystem &system,
                               WVAdaptiveRK78Options options)
    : system_(system), options_(options) {}

WVAdaptiveRK78::~WVAdaptiveRK78() {
  delete continuousExtension_;
  delete workspace_;
}

std::size_t WVAdaptiveRK78::persistentBytes() const noexcept {
  return sizeof(*this) +
         (workspace_ == nullptr
              ? 0
              : sizeof(Workspace) + workspace_->capacityBytes() +
                    workspace_->coefficientViewBytes() +
                    workspace_->externalViewBytes()) +
         (continuousExtension_ == nullptr
              ? 0
              : sizeof(ContinuousExtensionWorkspace) +
                    continuousExtension_->capacityBytes() +
                    continuousExtension_->coefficientViewBytes()) +
         (errorPolicy_ == nullptr ? 0 : errorPolicy_->persistentBytes()) +
         stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK78StepDiagnostic) +
         toleranceComponentHashes_.capacity() * sizeof(std::uint64_t);
}

WVKernelStatus
WVAdaptiveRK78::ensureWorkspace(const WVMutableIntegrationState &state) {
  auto status = validateMutableIntegrationState(system_.stateLayout(), state);
  if (!status)
    return status;
  status = AdaptiveRungeKuttaDriver::validateOptions(options_, "RK78");
  if (!status)
    return status;
  if (!(options_.repeatedRejectionFactor > 0.0 &&
        options_.repeatedRejectionFactor < 1.0))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK78 repeated-rejection factor must be between zero and one."};
  if (workspace_ != nullptr)
    return WVKernelStatus::ok();
  status = AdaptiveRungeKuttaDriver::initializeErrorPolicy(
      system_, options_, errorPolicy_, toleranceHash_,
      toleranceComponentHashes_);
  if (!status)
    return status;
  try {
    workspace_ = new Workspace;
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "RK78 workspace allocation failed."};
  }
  IntegrationBuffer *buffers[] = {
      &workspace_->stage,       &workspace_->k1,
      &workspace_->k2OrK3OrK5, &workspace_->k4OrK13,
      &workspace_->k6,          &workspace_->k7,
      &workspace_->k8,          &workspace_->k9,
      &workspace_->k10,         &workspace_->k11,
      &workspace_->k12};
  for (auto *buffer : buffers) {
    status = buffer->initialize(system_.stateLayout());
    if (!status) {
      delete workspace_;
      workspace_ = nullptr;
      return status;
    }
  }
  try {
    workspace_->baseCoefficientViews.reserve(
        system_.stateLayout().coefficientFamilyCount());
    workspace_->baseBlockViews.reserve(
        system_.stateLayout().additionalBlocks().size());
    workspace_->acceptedCoefficientViews.reserve(
        system_.stateLayout().coefficientFamilyCount());
    workspace_->acceptedBlockViews.reserve(
        system_.stateLayout().additionalBlocks().size());
  } catch (const std::bad_alloc &) {
    delete workspace_;
    workspace_ = nullptr;
    return {WVKernelStatusCode::allocationFailure,
            "RK78 state-view workspace allocation failed."};
  }
  metrics_.workspaceCapacityBytes =
      workspace_->capacityBytes() + workspace_->coefficientViewBytes() +
      workspace_->externalViewBytes();
  metrics_.workspaceMaximumLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
  metrics_.stateCapacityBytes = workspace_->stage.valueCapacityBytes();
  metrics_.workspaceStateEquivalentCount = 11;
  metrics_.workspaceMaximumLiveStateEquivalentCount = 11;
  metrics_.retainedBaseStageCapacityBytes =
      options_.retainDenseOutput
          ? workspace_->k1.valueCapacityBytes() +
                workspace_->k6.valueCapacityBytes() +
                workspace_->k7.valueCapacityBytes() +
                workspace_->k8.valueCapacityBytes() +
                workspace_->k9.valueCapacityBytes() +
                workspace_->k10.valueCapacityBytes() +
                workspace_->k11.valueCapacityBytes() +
                workspace_->k12.valueCapacityBytes()
          : 0;
  metrics_.retainedBaseStageStateEquivalentCount =
      options_.retainDenseOutput ? 8 : 0;
  metrics_.denseHistoryCapacityBytes =
      metrics_.retainedBaseStageCapacityBytes;
  metrics_.denseHistoryStateEquivalentCount =
      metrics_.retainedBaseStageStateEquivalentCount;
  metrics_.errorPolicyBytes = errorPolicy_->persistentBytes();
  return WVKernelStatus::ok();
}

void WVAdaptiveRK78::releaseContinuousExtension() const noexcept {
  delete continuousExtension_;
  continuousExtension_ = nullptr;
  continuousExtensionReady_ = false;
  metrics_.continuousExtensionWorkspaceCapacityBytes = 0;
  metrics_.continuousExtensionWorkspaceStateEquivalentCount = 0;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes;
}

WVKernelStatus WVAdaptiveRK78::ensureContinuousExtension() const {
  if (continuousExtension_ != nullptr)
    return WVKernelStatus::ok();
  auto *extension = new (std::nothrow) ContinuousExtensionWorkspace;
  if (extension == nullptr)
    return {WVKernelStatusCode::allocationFailure,
            "RK78 continuous-extension workspace allocation failed."};
  IntegrationBuffer *buffers[] = {&extension->k14, &extension->k15,
                                  &extension->k16, &extension->k17};
  for (auto *buffer : buffers) {
    const auto status = buffer->initialize(system_.stateLayout());
    if (!status) {
      delete extension;
      return status;
    }
  }
  continuousExtension_ = extension;
  const auto bytes = extension->capacityBytes() +
                     extension->coefficientViewBytes();
  metrics_.continuousExtensionWorkspaceCapacityBytes = bytes;
  metrics_.continuousExtensionWorkspaceStateEquivalentCount = 4;
  metrics_.continuousExtensionWorkspaceMaximumLiveBytes =
      std::max(metrics_.continuousExtensionWorkspaceMaximumLiveBytes, bytes);
  metrics_.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount = 4;
  metrics_.workspaceLiveBytes = metrics_.workspaceCapacityBytes + bytes;
  metrics_.workspaceMaximumLiveBytes =
      std::max(metrics_.workspaceMaximumLiveBytes,
               metrics_.workspaceLiveBytes);
  metrics_.workspaceMaximumLiveStateEquivalentCount =
      std::max<std::size_t>(metrics_.workspaceMaximumLiveStateEquivalentCount,
                            15);
  return WVKernelStatus::ok();
}

WVKernelStatus WVAdaptiveRK78::prepareStateAfterRestart(
    WVMutableIntegrationState &state) {
  releaseContinuousExtension();
  return AdaptiveRungeKuttaDriver::prepareStateAfterRestart(
      system_, state, hasAcceptedStep_, derivativeReuseAvailable_,
      nextStepSize_, stepDiagnostics_, metrics_,
      [&]() { return ensureWorkspace(state); });
}

WVKernelStatus WVAdaptiveRK78::step(WVMutableIntegrationState &state,
                                    double proposedStepSize) {
  return stepImplementation(state, proposedStepSize, false);
}

WVKernelStatus WVAdaptiveRK78::stepImplementation(
    WVMutableIntegrationState &state, double proposedStepSize,
    bool allowFinalStepStretch) {
  if (stepping_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK78 stepping is not reentrant."};
  if (!(proposedStepSize > 0.0) || !std::isfinite(proposedStepSize))
    return {WVKernelStatusCode::invalidConfiguration,
            "RK78 step size must be finite and positive."};
  auto status = ensureWorkspace(state);
  if (!status)
    return status;
  releaseContinuousExtension();
  stepping_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{stepping_};
  hasAcceptedStep_ = false;
  derivativeReuseAvailable_ = false;
  const auto baseView = integrationConstView(
      state, workspace_->baseCoefficientViews, workspace_->baseBlockViews);
  const auto t = state.waveVortex.t;
  const auto t0 = state.waveVortex.t0;
  double h = allowFinalStepStretch
                 ? proposedStepSize
                 : std::min(proposedStepSize, options_.maximumStepSize);
  bool initialDerivativeAvailable = false;
  std::size_t rejectedThisStep = 0;
  std::size_t evaluationsThisStep = 0;

  const auto evaluateStage = [&](IntegrationBuffer &derivative,
                                 double stageTime) {
    const auto before = metrics_.rightHandSideEvaluationCount;
    const auto result = evaluate(system_, workspace_->stage, stageTime, t0,
                                 derivative, metrics_);
    evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
    metrics_.baseRightHandSideEvaluationCount +=
        metrics_.rightHandSideEvaluationCount - before;
    return result;
  };

  for (;;) {
    if (!std::isfinite(h) || !(t + h > t))
      return {WVKernelStatusCode::numericalFailure,
              "RK78 cannot advance time with the proposed step."};
    if (!initialDerivativeAvailable) {
      const auto before = metrics_.rightHandSideEvaluationCount;
      auto derivative = workspace_->k1.flux();
      status = system_.evaluateRightHandSide(baseView, derivative);
      if (status)
        ++metrics_.rightHandSideEvaluationCount;
      evaluationsThisStep += metrics_.rightHandSideEvaluationCount - before;
      metrics_.baseRightHandSideEvaluationCount +=
          metrics_.rightHandSideEvaluationCount - before;
      if (!status)
        return status;
      initialDerivativeAvailable = true;
    }

    workspace_->stage.setAffine(baseView, workspace_->k1,
                                h * RK78MethodPolicy::a21);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c2, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k2OrK3OrK5,
                           t + h * RK78MethodPolicy::c2);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a31);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a32);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c3, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k2OrK3OrK5,
                           t + h * RK78MethodPolicy::c3);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a41);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a43);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c4, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k4OrK13,
                           t + h * RK78MethodPolicy::c4);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a51);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a53);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a54);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c5, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k2OrK3OrK5,
                           t + h * RK78MethodPolicy::c5);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a61);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a64);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a65);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c6, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k6, t + h * RK78MethodPolicy::c6);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a71);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a74);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a75);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a76);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c7, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k7, t + h * RK78MethodPolicy::c7);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a81);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a84);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a85);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a86);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a87);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c8, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k8, t + h * RK78MethodPolicy::c8);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a91);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a94);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a95);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a96);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a97);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::a98);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c9, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k9, t + h * RK78MethodPolicy::c9);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a101);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a104);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a105);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a106);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a107);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::a108);
    workspace_->stage.addScaled(workspace_->k9, RK78MethodPolicy::a109);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c10, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k10, t + h * RK78MethodPolicy::c10);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a111);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a114);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a115);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a116);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a117);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::a118);
    workspace_->stage.addScaled(workspace_->k9, RK78MethodPolicy::a119);
    workspace_->stage.addScaled(workspace_->k10,
                                RK78MethodPolicy::a1110);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage,
                       t + h * RK78MethodPolicy::c11, t0);
    if (!status)
      return status;
    status =
        evaluateStage(workspace_->k11, t + h * RK78MethodPolicy::c11);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a121);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a124);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a125);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a126);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a127);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::a128);
    workspace_->stage.addScaled(workspace_->k9, RK78MethodPolicy::a129);
    workspace_->stage.addScaled(workspace_->k10,
                                RK78MethodPolicy::a1210);
    workspace_->stage.addScaled(workspace_->k11,
                                RK78MethodPolicy::a1211);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage, t + h, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k12, t + h);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::a131);
    workspace_->stage.addScaled(workspace_->k4OrK13,
                                RK78MethodPolicy::a134);
    workspace_->stage.addScaled(workspace_->k2OrK3OrK5,
                                RK78MethodPolicy::a135);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::a136);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::a137);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::a138);
    workspace_->stage.addScaled(workspace_->k9, RK78MethodPolicy::a139);
    workspace_->stage.addScaled(workspace_->k10,
                                RK78MethodPolicy::a1310);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    status = constrain(system_, workspace_->stage, t + h, t0);
    if (!status)
      return status;
    status = evaluateStage(workspace_->k4OrK13, t + h);
    if (!status)
      return status;

    workspace_->stage.setScaled(workspace_->k1, RK78MethodPolicy::b1);
    workspace_->stage.addScaled(workspace_->k6, RK78MethodPolicy::b6);
    workspace_->stage.addScaled(workspace_->k7, RK78MethodPolicy::b7);
    workspace_->stage.addScaled(workspace_->k8, RK78MethodPolicy::b8);
    workspace_->stage.addScaled(workspace_->k9, RK78MethodPolicy::b9);
    workspace_->stage.addScaled(workspace_->k10, RK78MethodPolicy::b10);
    workspace_->stage.addScaled(workspace_->k11, RK78MethodPolicy::b11);
    workspace_->stage.addScaled(workspace_->k12, RK78MethodPolicy::b12);
    workspace_->stage.setAffine(baseView, workspace_->stage, h);
    auto candidateState = workspace_->stage.mutableState(t + h, t0);
    const auto endpointConstraint =
        system_.enforceStateConstraints(candidateState);
    metrics_.constraintModifiedCoefficientCount +=
        endpointConstraint.modifiedCoefficientCount;
    if (!endpointConstraint)
      return endpointConstraint.status;

    const auto weightedComplexError = [&](std::size_t index,
                                          bool imaginary) noexcept {
      const auto component = [&](const IntegrationBuffer &buffer) {
        return imaginary ? buffer.complex()[index].imag
                         : buffer.complex()[index].real;
      };
      return h * (RK78MethodPolicy::e1 * component(workspace_->k1) +
                  RK78MethodPolicy::e6 * component(workspace_->k6) +
                  RK78MethodPolicy::e7 * component(workspace_->k7) +
                  RK78MethodPolicy::e8 * component(workspace_->k8) +
                  RK78MethodPolicy::e9 * component(workspace_->k9) +
                  RK78MethodPolicy::e10 * component(workspace_->k10) +
                  RK78MethodPolicy::e11 * component(workspace_->k11) +
                  RK78MethodPolicy::e12 * component(workspace_->k12) +
                  RK78MethodPolicy::e13 * component(workspace_->k4OrK13));
    };
    const auto complexError = [&](std::size_t index) noexcept {
      return WVComplex64{weightedComplexError(index, false),
                         weightedComplexError(index, true)};
    };
    const auto realError = [&](std::size_t index) noexcept {
      return h * (RK78MethodPolicy::e1 * workspace_->k1.real()[index] +
                  RK78MethodPolicy::e6 * workspace_->k6.real()[index] +
                  RK78MethodPolicy::e7 * workspace_->k7.real()[index] +
                  RK78MethodPolicy::e8 * workspace_->k8.real()[index] +
                  RK78MethodPolicy::e9 * workspace_->k9.real()[index] +
                  RK78MethodPolicy::e10 * workspace_->k10.real()[index] +
                  RK78MethodPolicy::e11 * workspace_->k11.real()[index] +
                  RK78MethodPolicy::e12 * workspace_->k12.real()[index] +
                  RK78MethodPolicy::e13 * workspace_->k4OrK13.real()[index]);
    };
    const auto error = AdaptiveRungeKuttaDriver::normalizedError(
        *errorPolicy_, system_.stateLayout(), workspace_->stage, baseView,
        options_.relativeTolerance, complexError, realError,
        rejectedThisStep == 0);
    const auto accepted = std::isfinite(error) && error <= 1.0;
    const auto factor =
        AdaptiveRungeKuttaDriver::controllerFactor<RK78MethodPolicy>(
            error, rejectedThisStep, options_.safetyFactor,
            options_.rejectionFloorFactor, options_.repeatedRejectionFactor,
            options_.maximumStepFactor);
    nextStepSize_ = std::min(options_.maximumStepSize, h * factor);
    metrics_.lastProposedStepSize = proposedStepSize;
    metrics_.normalizedError = error;
    metrics_.nextStepSize = nextStepSize_;
    if (accepted) {
      if (options_.retainDenseOutput)
        workspace_->k2OrK3OrK5.copyFrom(baseView);
      workspace_->stage.copyTo(state);
      state.waveVortex.t = t + h;
      makeExternalViews(state, workspace_->acceptedCoefficientViews,
                        workspace_->acceptedBlockViews);
      acceptedStep_ = {
          t,
          state.waveVortex.t,
          {state.waveVortex.view(), workspace_->acceptedBlockViews.data(),
           workspace_->acceptedBlockViews.size(),
           workspace_->acceptedCoefficientViews.data(),
           workspace_->acceptedCoefficientViews.size()},
          {metrics_.acceptedStepCount + 1, rejectedThisStep,
           evaluationsThisStep, h, proposedStepSize, nextStepSize_, error},
          options_.retainDenseOutput ? this : nullptr};
      ++metrics_.acceptedStepCount;
      ++metrics_.stepCount;
      metrics_.lastStepSize = h;
      metrics_.lastAcceptedStepSize = h;
      derivativeReuseAvailable_ = false;
      hasAcceptedStep_ = true;
      if (stepDiagnostics_.size() < options_.maximumRecordedStepDiagnostics)
        stepDiagnostics_.push_back({t, h, error, nextStepSize_,
                                    rejectedThisStep, evaluationsThisStep,
                                    false});
      metrics_.diagnosticCapacityBytes =
          stepDiagnostics_.capacity() * sizeof(WVAdaptiveRK78StepDiagnostic);
      return WVKernelStatus::ok();
    }
    const auto attemptedStepSize = h;
    AdaptiveRungeKuttaDriver::recordRejectedAttempt(
        metrics_, rejectedThisStep, nextStepSize_, h);
    if (attemptedStepSize <= AdaptiveRungeKuttaDriver::minimumStepSize(t))
      return {WVKernelStatusCode::numericalFailure,
              "RK78 cannot meet integration tolerances at the minimum step size."};
    if (!(h > 0.0) || t + h == t)
      return {WVKernelStatusCode::numericalFailure,
              "RK78 step size underflowed after rejection."};
  }
}

WVKernelStatus WVAdaptiveRK78::advanceToTime(
    WVMutableIntegrationState &state, double finalTime, double h) {
  return AdaptiveRungeKuttaDriver::advanceToTime(
      state, finalTime, h, true, "RK78",
      [&](double use) {
        return stepImplementation(
            state, use,
            use > options_.maximumStepSize &&
                use <= 1.1 * options_.maximumStepSize);
      },
      [&]() { return nextStepSize_; });
}

WVKernelStatus WVAdaptiveRK78::evaluateDenseOutput(
    double time, WVMutableIntegrationState &output) const {
  if (!options_.retainDenseOutput || !hasAcceptedStep_)
    return {WVKernelStatusCode::unsupportedOperation,
            "RK78 dense output is unavailable."};
  if (evaluatingDenseOutput_)
    return {WVKernelStatusCode::reentrantExecution,
            "RK78 dense-output evaluation is not reentrant."};
  auto status = validateMutableIntegrationState(system_.stateLayout(), output);
  if (!status)
    return status;
  for (std::size_t family = 0; family < output.coefficientFamilyCount;
       ++family)
    if (output.coefficientFamilies[family].data ==
        acceptedStep_.endpoint.coefficientFamilies[family].data)
      return {WVKernelStatusCode::invalidConfiguration,
              "RK78 dense output must not alias accepted integration state."};
  for (std::size_t block = 0; block < output.additionalBlockCount; ++block)
    if ((output.additionalBlocks[block].realData != nullptr &&
         output.additionalBlocks[block].realData ==
             acceptedStep_.endpoint.additionalBlocks[block].realData) ||
        (output.additionalBlocks[block].complexData != nullptr &&
         output.additionalBlocks[block].complexData ==
             acceptedStep_.endpoint.additionalBlocks[block].complexData))
      return {WVKernelStatusCode::invalidConfiguration,
              "RK78 dense output must not alias accepted integration state."};
  const auto tolerance =
      timeTolerance(acceptedStep_.initialTime, acceptedStep_.finalTime);
  if (time < acceptedStep_.initialTime - tolerance ||
      time > acceptedStep_.finalTime + tolerance)
    return {WVKernelStatusCode::invalidConfiguration,
            "Dense-output time is outside the accepted interval."};
  evaluatingDenseOutput_ = true;
  struct Guard {
    bool &value;
    ~Guard() { value = false; }
  } guard{evaluatingDenseOutput_};
  const auto started = std::chrono::steady_clock::now();
  const double h = acceptedStep_.finalTime - acceptedStep_.initialTime;
  double theta = h == 0.0 ? 0.0 : (time - acceptedStep_.initialTime) / h;
  if (std::abs(time - acceptedStep_.initialTime) <= tolerance)
    theta = 0.0;
  if (std::abs(time - acceptedStep_.finalTime) <= tolerance)
    theta = 1.0;
  theta = std::max(0.0, std::min(1.0, theta));
  std::size_t denseReadMultiplier = 1;

  if (theta == 0.0) {
    workspace_->k2OrK3OrK5.copyTo(output);
  } else if (theta == 1.0) {
    workspace_->stage.copyFrom(acceptedStep_.endpoint);
    workspace_->stage.copyTo(output);
  } else {
    denseReadMultiplier = 13;
    const bool reusedCache = continuousExtensionReady_;
    if (!continuousExtensionReady_) {
      status = ensureContinuousExtension();
      if (!status)
        return status;
      const auto extensionStarted = std::chrono::steady_clock::now();
      const IntegrationBuffer *baseStages[] = {
          &workspace_->k1,  &workspace_->k6,  &workspace_->k7,
          &workspace_->k8,  &workspace_->k9,  &workspace_->k10,
          &workspace_->k11, &workspace_->k12};
      const auto setBaseStage = [&](const double *weights) {
        workspace_->stage.copyFrom(workspace_->k2OrK3OrK5.state(
            acceptedStep_.initialTime, acceptedStep_.endpoint.waveVortex.t0));
        for (std::size_t index = 0; index < 8; ++index)
          workspace_->stage.addScaled(*baseStages[index],
                                      h * weights[index]);
      };
      const auto evaluateExtensionStage = [&](IntegrationBuffer &derivative,
                                              double stageTime) {
        const auto before = metrics_.rightHandSideEvaluationCount;
        const auto result = evaluate(system_, workspace_->stage, stageTime,
                                     acceptedStep_.endpoint.waveVortex.t0,
                                     derivative, metrics_);
        metrics_.continuousExtensionRightHandSideEvaluationCount +=
            metrics_.rightHandSideEvaluationCount - before;
        return result;
      };
      const double baseWeights[] = {
          RK78MethodPolicy::b1,  RK78MethodPolicy::b6,
          RK78MethodPolicy::b7,  RK78MethodPolicy::b8,
          RK78MethodPolicy::b9,  RK78MethodPolicy::b10,
          RK78MethodPolicy::b11, RK78MethodPolicy::b12};
      setBaseStage(baseWeights);
      status = constrain(system_, workspace_->stage,
                         acceptedStep_.finalTime,
                         acceptedStep_.endpoint.waveVortex.t0);
      if (status)
        status = evaluateExtensionStage(continuousExtension_->k14,
                                        acceptedStep_.finalTime);
      if (status) {
        setBaseStage(RK78MethodPolicy::a15);
        workspace_->stage.addScaled(continuousExtension_->k14,
                                    h * RK78MethodPolicy::a15[8]);
        status = constrain(
            system_, workspace_->stage,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c15,
            acceptedStep_.endpoint.waveVortex.t0);
      }
      if (status)
        status = evaluateExtensionStage(
            continuousExtension_->k15,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c15);
      if (status) {
        setBaseStage(RK78MethodPolicy::a16);
        workspace_->stage.addScaled(continuousExtension_->k14,
                                    h * RK78MethodPolicy::a16[8]);
        workspace_->stage.addScaled(continuousExtension_->k15,
                                    h * RK78MethodPolicy::a16[9]);
        status = constrain(
            system_, workspace_->stage,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c16,
            acceptedStep_.endpoint.waveVortex.t0);
      }
      if (status)
        status = evaluateExtensionStage(
            continuousExtension_->k16,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c16);
      if (status) {
        setBaseStage(RK78MethodPolicy::a17);
        workspace_->stage.addScaled(continuousExtension_->k14,
                                    h * RK78MethodPolicy::a17[8]);
        workspace_->stage.addScaled(continuousExtension_->k15,
                                    h * RK78MethodPolicy::a17[9]);
        workspace_->stage.addScaled(continuousExtension_->k16,
                                    h * RK78MethodPolicy::a17[10]);
        status = constrain(
            system_, workspace_->stage,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c17,
            acceptedStep_.endpoint.waveVortex.t0);
      }
      if (status)
        status = evaluateExtensionStage(
            continuousExtension_->k17,
            acceptedStep_.initialTime + h * RK78MethodPolicy::c17);
      metrics_.continuousExtensionSeconds +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        extensionStarted)
              .count();
      if (!status) {
        releaseContinuousExtension();
        return status;
      }
      continuousExtensionReady_ = true;
      ++metrics_.denseOutputCacheBuildCount;
    }
    if (reusedCache)
      ++metrics_.denseOutputCacheReuseCount;

    double weights[12];
    for (std::size_t index = 0; index < 12; ++index) {
      double value = RK78MethodPolicy::continuousWeights[5][index];
      for (std::size_t power = 5; power-- > 0;)
        value = value * theta +
                RK78MethodPolicy::continuousWeights[power][index];
      weights[index] = value * theta * theta;
    }
    weights[0] += theta;
    workspace_->stage.copyFrom(workspace_->k2OrK3OrK5.state(
        acceptedStep_.initialTime, acceptedStep_.endpoint.waveVortex.t0));
    const IntegrationBuffer *stages[] = {
        &workspace_->k1,             &workspace_->k6,
        &workspace_->k7,             &workspace_->k8,
        &workspace_->k9,             &workspace_->k10,
        &workspace_->k11,            &workspace_->k12,
        &continuousExtension_->k14,  &continuousExtension_->k15,
        &continuousExtension_->k16,  &continuousExtension_->k17};
    for (std::size_t index = 0; index < 12; ++index)
      workspace_->stage.addScaled(*stages[index], h * weights[index]);
    auto interpolated = workspace_->stage.mutableState(
        time, acceptedStep_.endpoint.waveVortex.t0);
    status = system_.enforceStateConstraints(interpolated).status;
    if (!status)
      return status;
    workspace_->stage.copyTo(output);
  }
  output.waveVortex.t = time;
  output.waveVortex.t0 = acceptedStep_.endpoint.waveVortex.t0;
  ++metrics_.denseOutputEvaluationCount;
  const auto elementCount = workspace_->stage.complex().size() +
                            workspace_->stage.real().size();
  metrics_.denseOutputElementReads += denseReadMultiplier * elementCount;
  metrics_.denseOutputElementWrites += elementCount;
  metrics_.denseOutputSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return WVKernelStatus::ok();
}

double WVAdaptiveRK78::initialTime() const noexcept {
  return hasAcceptedStep_ ? acceptedStep_.initialTime : 0.0;
}

double WVAdaptiveRK78::finalTime() const noexcept {
  return hasAcceptedStep_ ? acceptedStep_.finalTime : 0.0;
}

const WVAcceptedStep *WVAdaptiveRK78::lastAcceptedStep() const noexcept {
  return hasAcceptedStep_ ? &acceptedStep_ : nullptr;
}

const WVIntegratorMetrics &WVAdaptiveRK78::metrics() const noexcept {
  return metrics_;
}

const std::vector<WVAdaptiveRK78StepDiagnostic> &
WVAdaptiveRK78::stepDiagnostics() const noexcept {
  return stepDiagnostics_;
}

std::uint64_t WVAdaptiveRK78::toleranceHash() const noexcept {
  return toleranceHash_;
}

const std::vector<std::uint64_t> &
WVAdaptiveRK78::toleranceComponentHashes() const noexcept {
  return toleranceComponentHashes_;
}

bool WVAdaptiveRK78::stepDiagnosticsComplete() const noexcept {
  return stepDiagnostics_.size() == metrics_.acceptedStepCount;
}

const char *WVAdaptiveRK78::controllerIdentifier() noexcept {
  return RK78MethodPolicy::controllerIdentifier;
}

const char *WVAdaptiveRK78::methodIdentifier() noexcept {
  return RK78MethodPolicy::methodIdentifier;
}

const WVAdaptiveRKStageBufferLastUse *
WVAdaptiveRK78::stageBufferLastUseRecords() noexcept {
  return RK78MethodPolicy::stageBufferLastUse;
}

std::size_t WVAdaptiveRK78::stageBufferLastUseRecordCount() noexcept {
  return sizeof(RK78MethodPolicy::stageBufferLastUse) /
         sizeof(RK78MethodPolicy::stageBufferLastUse[0]);
}

double WVAdaptiveRK78::nextStepSize() const noexcept { return nextStepSize_; }

} // namespace wavevortex::runtime
