#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "WVForcingImplementations.hpp"

#include <utility>

namespace wavevortex::runtime {
namespace {

} // namespace

std::vector<WVForcingFactoryRegistration> builtInForcingFactories() {
  using Encoding = WVForcingPersistenceEncoding;
  using Dimensions = WVForcingDimensionRule;
  const auto field = [](Encoding encoding, std::string recordName,
                        std::string imaginaryRecordName,
                        std::string netcdfName, Dimensions dimensions,
                        std::string reference = {}, bool optional = false,
                        bool nonnegative = false, bool positive = false,
                        bool allowInfinity = false) {
    return WVForcingPersistenceField{encoding, std::move(recordName),
                                     std::move(imaginaryRecordName),
                                     std::move(netcdfName), dimensions,
                                     std::move(reference), optional,
                                     nonnegative, positive, allowInfinity};
  };
  return {
      {"WVNonlinearAdvection", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "nonlinear advection", WVForcingStage::spatial, 127, {},
       detail::createNonlinearAdvectionForcing, true, "", false,
       WVForcingStage::spatial,
       detail::preflightBarotropicQGEmptyForcing,
       detail::createBarotropicQGNonlinearAdvection},
      {"WVAntialiasing", WVPortablePairContractVersion,
       {"Spectral", "PVSpectral"}, "antialiasing", WVForcingStage::spectral,
       255, {}, {}, false,
       "Transform-level antialiasing is represented by shouldAntialias; the "
       "diagnostic WVAntialiasing closure is not supported.", false,
       WVForcingStage::spatial, {}, {}},
      {"WVAdaptiveDamping", WVPortablePairContractVersion,
       {"Spectral", "PVSpectral"}, "adaptive damping",
       WVForcingStage::spectral, 255, {}, detail::createAdaptiveDampingForcing,
       true, "", true, WVForcingStage::spectral,
       detail::preflightBarotropicQGEmptyForcing,
       detail::createBarotropicQGAdaptiveDamping},
      {"WVFixedAmplitudeForcing", WVPortablePairContractVersion,
       {"SpectralAmplitude", "PVSpectralAmplitude"}, "fixed amplitude",
       WVForcingStage::spectralAmplitude, 255,
       {{{field(Encoding::zeroBasedIndexVariable, "ApIndices", {},
                "Ap_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "ApValuesReal", "ApValuesImag",
                "Apbar", Dimensions::referencedLength, "Ap_indices", true),
          field(Encoding::zeroBasedIndexVariable, "AmIndices", {},
                "Am_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "AmValuesReal", "AmValuesImag",
                "Ambar", Dimensions::referencedLength, "Am_indices", true),
          field(Encoding::zeroBasedIndexVariable, "A0Indices", {},
                "A0_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "A0ValuesReal", "A0ValuesImag",
                "A0bar", Dimensions::referencedLength, "A0_indices", true)}},
        true},
       detail::createFixedAmplitudeForcing, true, "", false,
       WVForcingStage::spectralAmplitude,
       detail::preflightBarotropicQGFixedAmplitude,
       detail::createBarotropicQGFixedAmplitude},
      {"WVBottomFrictionQuadratic", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "quadratic bottom friction", WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "Cd", {}, "Cd", Dimensions::scalar,
                {}, false, true)}},
        false},
       detail::createQuadraticBottomFriction, true, "", false,
       WVForcingStage::spatial,
       detail::preflightBarotropicQGScalarForcing,
       detail::createBarotropicQGQuadraticBottomFriction},
      {"WVPseudoTopographicWaveGeneration", WVPortablePairContractVersion,
       {"Spectral"}, "pseudo-topographic wave generation",
       WVForcingStage::spectral, 255,
       {{{field(Encoding::realVariable, "topographicHeight", {},
                "topographicHeight", Dimensions::horizontalYX),
          field(Encoding::complexVariable, "barotropicVelocityAmplitudeReal",
                "barotropicVelocityAmplitudeImag",
                "barotropicVelocityAmplitude", Dimensions::componentPair),
          field(Encoding::realVariable, "frequency", {}, "frequency",
                Dimensions::scalar, {}, false, false, true),
          field(Encoding::textAttribute, "darwinSymbol", {}, "darwinSymbol",
                Dimensions::scalar, {}, true),
          field(Encoding::realVariable, "rampDuration", {}, "rampDuration",
                Dimensions::scalar, {}, false, true),
          field(Encoding::realVariable, "startTime", {}, "startTime",
                Dimensions::scalar),
          field(Encoding::logicalVariable, "shouldAvoidAdaptiveDamping", {},
                "shouldAvoidAdaptiveDamping", Dimensions::scalar),
          field(Encoding::realVariable, "maximumForcedHorizontalWavenumber",
                {}, "maximumForcedHorizontalWavenumber", Dimensions::scalar,
                {}, false, true, false, true),
          field(Encoding::realVariable, "maximumForcedVerticalMode", {},
                "maximumForcedVerticalMode", Dimensions::scalar, {}, false,
                true, false, true)}},
        true},
       detail::createPseudoTopographicForcing, true, "", false,
       WVForcingStage::spatial, {}, {}},
      {"WVBetaPlanePVAdvection", WVPortablePairContractVersion,
       {"Spectral", "PVSpatial"}, "beta-plane advection of qgpv",
       WVForcingStage::spectral, 255, {}, detail::createBetaPlaneForcing, true,
       "", false, WVForcingStage::spatial,
       detail::preflightBarotropicQGEmptyForcing,
       detail::createBarotropicQGBetaPlanePVAdvection},
      {"WVHorizontalDamping", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial"}, "horizontal damping",
       WVForcingStage::spectral, 255, {}, {}, false,
       "WVHorizontalDamping is not implemented by portable runtime v1.",
       false, WVForcingStage::spatial, {}, {}},
      {"WVVerticalDamping", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial"}, "vertical damping",
       WVForcingStage::spectral, 255, {}, {}, false,
       "WVVerticalDamping is not implemented by portable runtime v1.",
       false, WVForcingStage::spatial, {}, {}},
      {"WVThermalDamping", WVPortablePairContractVersion, {"PVSpatial"},
       "thermal damping", WVForcingStage::spectral, 255, {}, {}, false,
       "WVThermalDamping is not implemented by portable runtime v1.",
       false, WVForcingStage::spatial, {}, {}},
      {"WVBottomFrictionLinear", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "linear bottom friction", WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "r", {}, "r", Dimensions::scalar,
                {}, false, true)}},
        false},
       detail::createLinearBottomFriction, true, "", false,
       WVForcingStage::spatial,
       detail::preflightBarotropicQGScalarForcing,
       detail::createBarotropicQGLinearBottomFriction},
      {"WVVerticalDiffusivity", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "vertical diffusivity", WVForcingStage::spectral, 255, {}, {}, false,
       "WVVerticalDiffusivity is not implemented by portable runtime v1.",
       false, WVForcingStage::spatial, {}, {}},
      {"WVNarrowBandGeostrophicForcing", WVPortablePairContractVersion,
       {"PVSpectralAmplitude"}, "narrow-band geostrophic forcing",
       WVForcingStage::spectralAmplitude, 255,
       {{{field(Encoding::zeroBasedIndexVariable, "ApIndices", {},
                "Ap_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "ApValuesReal", "ApValuesImag",
                "Apbar", Dimensions::referencedLength, "Ap_indices", true),
          field(Encoding::zeroBasedIndexVariable, "AmIndices", {},
                "Am_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "AmValuesReal", "AmValuesImag",
                "Ambar", Dimensions::referencedLength, "Am_indices", true),
          field(Encoding::zeroBasedIndexVariable, "A0Indices", {},
                "A0_indices", Dimensions::ownLength, {}, true),
          field(Encoding::complexVariable, "A0ValuesReal", "A0ValuesImag",
                "A0bar", Dimensions::referencedLength, "A0_indices", true),
          field(Encoding::realVariable, "r", {}, "r", Dimensions::scalar,
                {}, true, true),
          field(Encoding::realVariable, "k_r", {}, "k_r",
                Dimensions::scalar, {}, true, true),
          field(Encoding::realVariable, "k_f", {}, "k_f",
                Dimensions::scalar, {}, true, true),
          field(Encoding::realVariable, "j_f", {}, "j_f",
                Dimensions::scalar, {}, true, true),
          field(Encoding::realVariable, "u_rms", {}, "u_rms",
                Dimensions::scalar, {}, true, true),
          field(Encoding::textAttribute, "initialPV", {}, "initialPV",
                Dimensions::scalar, {}, true)}},
        true},
       {}, true, "", false, WVForcingStage::spectralAmplitude,
       detail::preflightBarotropicQGFixedAmplitude,
       detail::createBarotropicQGFixedAmplitude}};
}

WVFrozenForcingSchedule defaultNonlinearAdvectionSchedule() {
  WVFrozenForcingSchedule schedule;
  WVFrozenForcingEntry entry;
  entry.typeIdentifier = "WVNonlinearAdvection";
  entry.contractVersion = WVPortablePairContractVersion;
  entry.name = "nonlinear advection";
  entry.stage = WVForcingStage::spatial;
  entry.priority = 127;
  entry.ordinal = 1;
  entry.configuration.schemaIdentifier =
      "wave-vortex-forcing-configuration-v1";
  entry.configuration.schemaVersion = 1;
  schedule.entries.push_back(std::move(entry));
  return schedule;
}

} // namespace wavevortex::runtime
