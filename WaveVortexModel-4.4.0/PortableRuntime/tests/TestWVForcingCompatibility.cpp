#include "WaveVortexRuntime/WVForcingContracts.hpp"
#include "generated/WVForcingCompatibilityRows.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVForcingEngine.hpp"
#include "WaveVortexRuntime/WVHydrostaticForcingEngine.hpp"
#include "WaveVortexRuntime/WVBoussinesqForcingEngine.hpp"
#include "WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp"
#include "WaveVortexRuntime/WVStratifiedQGForcingEngine.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace wavevortex;
using namespace wavevortex::runtime;

void testForcingIncompatibilities() {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (!status) throw std::runtime_error(status.message);
  std::shared_ptr<const WVExtensionCatalog> catalog;
  status = builder.freeze(catalog);
  if (!status) throw std::runtime_error(status.message);
  std::size_t rejected = 0;
  for (const auto& row : test::forcingCompatibilityRows) {
    if (row.applicable) continue;
    const auto* registration = catalog->forcings().registration(row.identity,row.version);
    WVFrozenForcingEntry entry;
    entry.typeIdentifier=row.identity; entry.contractVersion=row.version;
    entry.name=registration->defaultName;
    entry.stage=row.stratifiedQG ? registration->stratifiedQGStage : row.barotropic ? registration->barotropicQGStage : registration->stage;
    entry.priority=registration->priority;
    entry.configuration={"wave-vortex-forcing-configuration-v1",1,{}};
    if (entry.typeIdentifier=="WVAntialiasing")
      entry.configuration.values.push_back({"Nj",{},std::vector<double>{1}});
    WVFrozenForcingSchedule schedule; schedule.entries.push_back(entry);
    if (row.boussinesq) {
      WVStratifiedModalGeometry c; c.transformClass="WVTransformBoussinesq"; c.Nj=4; c.Nkl=24; c.shouldAntialias=row.antialias;
      status=WVBoussinesqForcingEngine::validateSchedule(c,schedule,{4,24},*catalog);
    } else if (row.hydrostatic) {
      WVStratifiedModalGeometry c; c.transformClass="WVTransformHydrostatic"; c.Nj=4; c.Nkl=24; c.shouldAntialias=row.antialias;
      status=WVHydrostaticForcingEngine::validateSchedule(c,schedule,{4,24},*catalog);
    } else if (row.stratifiedQG) {
      WVStratifiedModalGeometry c; c.Nx=8; c.Ny=6; c.Nz=9; c.Nj=4; c.Nkl=24; c.shouldAntialias=row.antialias;
      status=WVStratifiedQGForcingEngine::validateSchedule(c,schedule,96,*catalog);
    } else if (row.barotropic) {
      WVTransformBarotropicQGConfiguration c;
      c.Nx=8; c.Ny=6; c.Lx=17000; c.Ly=11000; c.h=1; c.j=1; c.shouldAntialias=row.antialias;
      status=WVBarotropicQGForcingEngine::validateSchedule(c,schedule,24,*catalog);
    } else {
      WVTransformConstantStratificationConfiguration c;
      c.Nx=8; c.Ny=6; c.Nz=5; c.Nj=4; c.shouldAntialias=row.antialias;
      status=WVConstantStratificationForcingEngine::validateSchedule(c,schedule,{4,24},*catalog);
    }
    if (status || (status.code!=WVKernelStatusCode::unsupportedOperation && status.code!=WVKernelStatusCode::invalidConfiguration))
      throw std::runtime_error(std::string("MATLAB incompatibility accepted: ")+row.identity);
    ++rejected;
  }
  if (rejected!=20) throw std::runtime_error("Incomplete baseline rejection evidence.");
}

int main() {
  try {
    testForcingIncompatibilities();
    const auto registrations = builtInForcingFactories();
    for (const auto &row : test::forcingCompatibilityRows) {
      const auto found = std::find_if(registrations.begin(), registrations.end(),
          [&](const auto &entry) { return entry.matlabClassName == row.identity &&
                                        entry.contractVersion == row.version; });
      if (found == registrations.end())
        throw std::runtime_error(std::string("Missing current registry pair: ") + row.identity);
      const bool available = found->isSupported &&
          (row.boussinesq ? static_cast<bool>(found->boussinesqFactory) : row.hydrostatic ? static_cast<bool>(found->hydrostaticFactory) : row.stratifiedQG ? static_cast<bool>(found->stratifiedQGFactory) : row.barotropic ? static_cast<bool>(found->barotropicQGFactory)
                          : static_cast<bool>(found->factory));
      if (available != row.factoryAvailable)
        throw std::runtime_error(std::string("Stale factory availability: ") + row.identity);
      // Unavailable placeholders do not define authoritative scientific stages.
      if (!row.applicable || !available) continue;
      const std::string stage(row.stage);
      WVForcingStage expected;
      if (stage == "HydrostaticSpatial" || stage == "NonhydrostaticSpatial" ||
          stage == "PVSpatial") expected = WVForcingStage::spatial;
      else if (stage == "Spectral" || stage == "PVSpectral")
        expected = WVForcingStage::spectral;
      else if (stage == "SpectralAmplitude" || stage == "PVSpectralAmplitude")
        expected = WVForcingStage::spectralAmplitude;
      else throw std::runtime_error("Unknown generated scientific stage.");
      const auto actual = row.stratifiedQG ? found->stratifiedQGStage : row.barotropic ? found->barotropicQGStage : found->stage;
      if (actual != expected || found->priority != row.priority ||
          std::find(found->forcingTypes.begin(), found->forcingTypes.end(), stage) ==
              found->forcingTypes.end())
        throw std::runtime_error(std::string("MATLAB/registry stage or priority drift: ") + row.identity);
    }
    std::cout << "Checked " << test::forcingCompatibilityRows.size()
              << " MATLAB-derived rows against the current forcing registry.\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
