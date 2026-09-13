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
  std::vector<WVForcingFactoryRegistration> factories{
      {"WVNonlinearAdvection", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "nonlinear advection", WVForcingStage::spatial, 127, {},
       detail::createNonlinearAdvectionForcing, true, "", false,
       WVForcingStage::spatial,
       detail::preflightBarotropicQGEmptyForcing,
       detail::createBarotropicQGNonlinearAdvection, {}, {},
       WVForcingStage::spatial, detail::preflightStratifiedQGEmptyForcing,
       detail::createStratifiedQGNonlinearAdvection},
      {"WVAntialiasing", WVPortablePairContractVersion,
       {"Spectral", "PVSpectral"}, "antialias filter", WVForcingStage::spectral,
       127, {{{field(Encoding::realVariable, "Nj", {}, "Nj", Dimensions::scalar,
                     {}, false, true)}}, false},
       detail::createExplicitAntialiasing, true, "", false,
       WVForcingStage::spectral,
       detail::preflightBarotropicQGExplicitAntialiasing,
       detail::createBarotropicQGExplicitAntialiasing,
       detail::preflightExplicitAntialiasing, detail::prepareExplicitAntialiasing,
       WVForcingStage::spectral, detail::preflightStratifiedQGExplicitAntialiasing,
       detail::createStratifiedQGExplicitAntialiasing},
      {"WVAdaptiveDamping", WVPortablePairContractVersion,
       {"Spectral", "PVSpectral"}, "adaptive damping",
       WVForcingStage::spectral, 255, {}, detail::createAdaptiveDampingForcing,
       true, "", true, WVForcingStage::spectral,
       detail::preflightBarotropicQGEmptyForcing,
       detail::createBarotropicQGAdaptiveDamping, {}, {},
       WVForcingStage::spectral, detail::preflightStratifiedQGEmptyForcing,
       detail::createStratifiedQGAdaptiveDamping},
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
       detail::createBarotropicQGFixedAmplitude, {}, {},
       WVForcingStage::spectralAmplitude, detail::preflightStratifiedQGFixedAmplitude,
       detail::createStratifiedQGFixedAmplitude},
      {"WVBottomFrictionQuadratic", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"},
       "quadratic bottom friction", WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "Cd", {}, "Cd", Dimensions::scalar,
                {}, false, true)}},
        false},
       detail::createQuadraticBottomFriction, true, "", false,
       WVForcingStage::spatial,
       detail::preflightBarotropicQGScalarForcing,
       detail::createBarotropicQGQuadraticBottomFriction, {}, {},
       WVForcingStage::spatial, detail::preflightStratifiedQGScalarForcing,
       detail::createStratifiedQGQuadraticBottomFriction},
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
       detail::createBarotropicQGBetaPlanePVAdvection, {}, {},
       WVForcingStage::spatial, detail::preflightStratifiedQGEmptyForcing,
       detail::createStratifiedQGBetaPlanePVAdvection},
      {"WVHorizontalDamping", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial"}, "horizontal scalar diffusivity",
       WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "nu", {}, "nu", Dimensions::scalar, {}, false, true),
          field(Encoding::realVariable, "kappa", {}, "kappa", Dimensions::scalar, {}, false, true)}}, false},
       detail::createHorizontalDamping, true, "", false, WVForcingStage::spatial, {}, {},
       detail::preflightLaplacianDamping},
      {"WVVerticalDamping", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial"}, "vertical scalar diffusivity",
       WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "nu", {}, "nu", Dimensions::scalar, {}, false, true),
          field(Encoding::realVariable, "kappa", {}, "kappa", Dimensions::scalar, {}, false, true)}}, false},
       detail::createVerticalDamping, true, "", false, WVForcingStage::spatial, {}, {},
       detail::preflightLaplacianDamping},
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
       detail::createBarotropicQGLinearBottomFriction, {}, {},
       WVForcingStage::spatial, detail::preflightStratifiedQGScalarForcing,
       detail::createStratifiedQGLinearBottomFriction},
      {"WVVerticalDiffusivity", WVPortablePairContractVersion,
       {"HydrostaticSpatial", "NonhydrostaticSpatial", "PVSpatial"}, "vertical diffusivity",
       WVForcingStage::spatial, 255,
       {{{field(Encoding::realVariable, "kappa_z", {}, "kappa_z", Dimensions::scalar, {}, false, true),
          field(Encoding::logicalVariable, "shouldForceMeanDensityAnomaly", {}, "shouldForceMeanDensityAnomaly", Dimensions::scalar)}}, false},
       detail::createVerticalDiffusivity, true, "", false, WVForcingStage::spatial, {}, {},
       detail::preflightVerticalDiffusivity, {},
       WVForcingStage::spatial, detail::preflightStratifiedQGVerticalDiffusivity,
       detail::createStratifiedQGVerticalDiffusivity},
      {"WVNarrowBandGeostrophicForcing", WVPortablePairContractVersion,
       {"SpectralAmplitude", "PVSpectralAmplitude"}, "narrow-band geostrophic forcing",
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
       detail::createFixedAmplitudeForcing, true, "", false, WVForcingStage::spectralAmplitude,
       detail::preflightBarotropicQGFixedAmplitude,
       detail::createBarotropicQGFixedAmplitude, {}, {},
       WVForcingStage::spectralAmplitude, detail::preflightStratifiedQGFixedAmplitude,
       detail::createStratifiedQGFixedAmplitude}};
  for (auto& registration:factories) {
    registration.supportsTendencyDiagnostics=true;
    if (registration.matlabClassName=="WVNonlinearAdvection")
      registration.evaluationDependencies.nonlinearUseCount=1;
    if (registration.matlabClassName=="WVAdaptiveDamping")
      registration.evaluationDependencies.horizontalMaximumNeeded=true;
    if (registration.matlabClassName=="WVHorizontalDamping") {
      ++registration.evaluationDependencies.constantLaplacianUseCount[0];
      // u/v/eta require Dxx and Dyy in the hydrostatic family.
      for (const auto field:{0U,1U,3U}) for (std::size_t kind=0;kind<2;++kind)
        ++registration.evaluationDependencies.hydrostaticGridCalculusUseCount[4*field+kind];
      // u/v/w/eta require Dxx and Dyy in the Boussinesq family.
      for (std::size_t field=0;field<4;++field) for (std::size_t kind=0;kind<2;++kind)
        ++registration.evaluationDependencies.boussinesqGridCalculusUseCount[4*field+kind];
    }
    if (registration.matlabClassName=="WVVerticalDamping") {
      ++registration.evaluationDependencies.constantLaplacianUseCount[1];
      // Hydrostatic u/v use F-Dzz; eta uses G-Dzz.
      for (const auto field:{0U,1U})
        ++registration.evaluationDependencies.hydrostaticGridCalculusUseCount[4*field+2];
      ++registration.evaluationDependencies.hydrostaticGridCalculusUseCount[4*3+3];
      // Boussinesq u/v use F-Dzz; w/eta use G-Dzz.
      for (const auto field:{0U,1U})
        ++registration.evaluationDependencies.boussinesqGridCalculusUseCount[4*field+2];
      for (const auto field:{2U,3U})
        ++registration.evaluationDependencies.boussinesqGridCalculusUseCount[4*field+3];
    }
    if (registration.matlabClassName=="WVVerticalDiffusivity") {
      ++registration.evaluationDependencies.constantLaplacianUseCount[1];
      ++registration.evaluationDependencies.hydrostaticGridCalculusUseCount[4*3+3];
      ++registration.evaluationDependencies.boussinesqGridCalculusUseCount[4*3+3];
    }
    if (registration.matlabClassName=="WVNonlinearAdvection") registration.hydrostaticFactory=detail::createHydrostaticNonlinearAdvectionForcing;
    if (registration.matlabClassName=="WVNonlinearAdvection") registration.boussinesqFactory=detail::createBoussinesqNonlinearAdvectionForcing;
    if (registration.matlabClassName=="WVAntialiasing") { registration.hydrostaticFactory=detail::createHydrostaticExplicitAntialiasing; registration.prepareHydrostaticResolution=detail::prepareHydrostaticExplicitAntialiasing; }
    if (registration.matlabClassName=="WVAntialiasing") { registration.boussinesqFactory=detail::createBoussinesqExplicitAntialiasing; registration.prepareBoussinesqResolution=detail::prepareBoussinesqExplicitAntialiasing; }
    if (registration.matlabClassName=="WVAdaptiveDamping") registration.hydrostaticFactory=detail::createHydrostaticAdaptiveDampingForcing;
    if (registration.matlabClassName=="WVAdaptiveDamping") registration.boussinesqFactory=detail::createBoussinesqAdaptiveDampingForcing;
    if (registration.matlabClassName=="WVFixedAmplitudeForcing") registration.hydrostaticFactory=detail::createHydrostaticFixedAmplitudeForcing;
    if (registration.matlabClassName=="WVFixedAmplitudeForcing") registration.boussinesqFactory=detail::createBoussinesqFixedAmplitudeForcing;
    if (registration.matlabClassName=="WVNarrowBandGeostrophicForcing") registration.hydrostaticFactory=detail::createHydrostaticFixedAmplitudeForcing;
    if (registration.matlabClassName=="WVNarrowBandGeostrophicForcing") registration.boussinesqFactory=detail::createBoussinesqFixedAmplitudeForcing;
    if (registration.matlabClassName=="WVBottomFrictionQuadratic") registration.hydrostaticFactory=detail::createHydrostaticQuadraticBottomFriction;
    if (registration.matlabClassName=="WVBottomFrictionQuadratic") registration.boussinesqFactory=detail::createBoussinesqQuadraticBottomFriction;
    if (registration.matlabClassName=="WVBottomFrictionLinear") registration.hydrostaticFactory=detail::createHydrostaticLinearBottomFriction;
    if (registration.matlabClassName=="WVBottomFrictionLinear") registration.boussinesqFactory=detail::createBoussinesqLinearBottomFriction;
    if (registration.matlabClassName=="WVPseudoTopographicWaveGeneration") registration.hydrostaticFactory=detail::createHydrostaticPseudoTopographicForcing;
    if (registration.matlabClassName=="WVPseudoTopographicWaveGeneration") registration.boussinesqFactory=detail::createBoussinesqPseudoTopographicForcing;
    if (registration.matlabClassName=="WVBetaPlanePVAdvection") registration.hydrostaticFactory=detail::createHydrostaticBetaPlaneForcing;
    if (registration.matlabClassName=="WVBetaPlanePVAdvection") registration.boussinesqFactory=detail::createBoussinesqBetaPlaneForcing;
    if (registration.matlabClassName=="WVHorizontalDamping") registration.hydrostaticFactory=detail::createHydrostaticHorizontalDamping;
    if (registration.matlabClassName=="WVHorizontalDamping") registration.boussinesqFactory=detail::createBoussinesqHorizontalDamping;
    if (registration.matlabClassName=="WVVerticalDamping") registration.hydrostaticFactory=detail::createHydrostaticVerticalDamping;
    if (registration.matlabClassName=="WVVerticalDamping") registration.boussinesqFactory=detail::createBoussinesqVerticalDamping;
    if (registration.matlabClassName=="WVVerticalDiffusivity") registration.hydrostaticFactory=detail::createHydrostaticVerticalDiffusivity;
    if (registration.matlabClassName=="WVVerticalDiffusivity") registration.boussinesqFactory=detail::createBoussinesqVerticalDiffusivity;
  }
  return factories;
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
