#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

#include <cmath>
#include <cstddef>

namespace wavevortex::detail {

inline WVComplex64 add(WVComplex64 a, WVComplex64 b) noexcept { return {a.real + b.real, a.imag + b.imag}; }
inline WVComplex64 subtract(WVComplex64 a, WVComplex64 b) noexcept { return {a.real - b.real, a.imag - b.imag}; }
inline WVComplex64 multiply(WVComplex64 a, WVComplex64 b) noexcept { return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real}; }
inline WVComplex64 multiply(WVComplex64 a, double b) noexcept { return {a.real * b, a.imag * b}; }
inline WVComplex64 conjugate(WVComplex64 a) noexcept { return {a.real, -a.imag}; }
inline WVComplex64 phase(double angle) noexcept { return {std::cos(angle), std::sin(angle)}; }

struct EvolvedWaveVortexCoefficients {
    WVComplex64 Ap;
    WVComplex64 Am;
    WVComplex64 A0;
};

inline EvolvedWaveVortexCoefficients evolveWaveVortexCoefficients(WVComplex64 Ap, WVComplex64 Am, WVComplex64 A0, WVComplex64 phaseValue) noexcept {
    return {multiply(Ap,phaseValue),multiply(Am,conjugate(phaseValue)),A0};
}

template <std::size_t Target>
inline WVComplex64 coefficientValueForField(const WVConstantStratificationModes& modes, std::size_t index, const EvolvedWaveVortexCoefficients& coefficients) noexcept {
    if constexpr (Target == 0) return add(add(multiply(modes.UApField[index],coefficients.Ap),multiply(conjugate(modes.UApField[index]),coefficients.Am)),multiply(modes.UA0Field[index],coefficients.A0));
    if constexpr (Target == 1) return add(add(multiply(modes.VApField[index],coefficients.Ap),multiply(conjugate(modes.VApField[index]),coefficients.Am)),multiply(modes.VA0Field[index],coefficients.A0));
    if constexpr (Target == 2) return add(multiply(modes.WApField[index],coefficients.Ap),multiply(modes.WApField[index],coefficients.Am));
    return add(add(multiply(coefficients.Ap,modes.NApField[index]),multiply(coefficients.Am,-modes.NApField[index])),multiply(coefficients.A0,modes.NA0Field[index]));
}

struct DerivativeSpectrum {
    WVComplex64 value;
    WVComplex64 x;
    WVComplex64 y;
    WVComplex64 z;
};

template <bool CosineFamily, bool Conjugated>
inline DerivativeSpectrum normalizedDerivativeSpectrum(WVComplex64 value, double k, double l, double verticalWavenumber, bool isVerticalEndpoint) noexcept {
    if constexpr (Conjugated) value = conjugate(value);
    auto x = multiply(value,WVComplex64{0.0,Conjugated ? -k : k});
    auto y = multiply(value,WVComplex64{0.0,Conjugated ? -l : l});
    auto z = multiply(value,CosineFamily ? -verticalWavenumber : verticalWavenumber);
    if constexpr (CosineFamily) {
        x = multiply(x,0.5);
        y = multiply(y,0.5);
        z = isVerticalEndpoint ? WVComplex64{} : multiply(z,0.5);
    } else {
        x = isVerticalEndpoint ? WVComplex64{} : multiply(x,0.5);
        y = isVerticalEndpoint ? WVComplex64{} : multiply(y,0.5);
        z = multiply(z,0.5);
    }
    return {value,x,y,z};
}

template <bool FFamily, bool Conjugated>
inline DerivativeSpectrum fieldFamilyDerivativeSpectrum(WVComplex64 Apm, WVComplex64 A0, double waveScale, double geostrophicScale, double k, double l, double verticalWavenumber) noexcept {
    auto value = add(multiply(Apm,waveScale),multiply(A0,geostrophicScale));
    auto x = multiply(value,WVComplex64{0.0,k});
    auto y = multiply(value,WVComplex64{0.0,l});
    auto z = multiply(value,FFamily ? -verticalWavenumber : verticalWavenumber);
    if constexpr (Conjugated) {
        value = conjugate(value);
        x = conjugate(x);
        y = conjugate(y);
        z = conjugate(z);
    }
    return {value,x,y,z};
}

inline WVComplex64 horizontalVorticity(WVComplex64 U, WVComplex64 V, double k, double l) noexcept {
    return subtract(multiply(V,WVComplex64{0.0,k}),multiply(U,WVComplex64{0.0,l}));
}

inline WVComplex64 geostrophicCoefficient(WVComplex64 U, WVComplex64 V, WVComplex64 N, double k, double l, double vorticityScale, double buoyancyScale) noexcept {
    return add(multiply(horizontalVorticity(U,V,k,l),vorticityScale),multiply(N,buoyancyScale));
}

inline WVComplex64 buoyancyProjection(WVComplex64 N, WVComplex64 A0, double geostrophicFieldScale, double waveProjectionScale) noexcept {
    return multiply(subtract(N,multiply(A0,geostrophicFieldScale)),waveProjectionScale);
}

struct WaveCoefficientPair {
    WVComplex64 Ap;
    WVComplex64 Am;
};

inline WaveCoefficientPair waveCoefficientPair(WVComplex64 waveContribution, WVComplex64 buoyancyContribution) noexcept {
    return {add(waveContribution,buoyancyContribution),subtract(waveContribution,buoyancyContribution)};
}

inline WaveCoefficientPair inertialCoefficientPair(WVComplex64 U, WVComplex64 V, double scale) noexcept {
    const auto Ap = multiply(subtract(U,multiply(V,WVComplex64{0.0,1.0})),scale);
    return {Ap,conjugate(Ap)};
}

template <std::size_t Target>
inline WaveCoefficientPair inertialCoefficientPair(WVComplex64 value, double scale) noexcept {
    WVComplex64 Ap{};
    if constexpr (Target == 0) Ap = multiply(value,scale);
    if constexpr (Target == 1) Ap = multiply(multiply(value,WVComplex64{0.0,-1.0}),scale);
    return {Ap,conjugate(Ap)};
}

inline WaveCoefficientPair atReferenceTime(const WaveCoefficientPair& coefficients, WVComplex64 phaseValue) noexcept {
    return {multiply(coefficients.Ap,conjugate(phaseValue)),multiply(coefficients.Am,phaseValue)};
}

inline void accumulateAtReferenceTime(WVComplex64& Fp, WVComplex64& Fm, const WaveCoefficientPair& contribution, WVComplex64 phaseValue) noexcept {
    Fp = add(Fp,multiply(contribution.Ap,conjugate(phaseValue)));
    Fm = add(Fm,multiply(contribution.Am,phaseValue));
}

} // namespace wavevortex::detail
