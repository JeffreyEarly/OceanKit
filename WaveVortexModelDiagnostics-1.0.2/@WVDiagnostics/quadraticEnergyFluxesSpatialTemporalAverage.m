function forcing_fluxes = quadraticEnergyFluxesSpatialTemporalAverage(self,options)
% Compute spatial-temporal average of forcing fluxes.
%
% Compute spatial-temporal average of forcing fluxes
% Returns the spatial-temporal average of energy fluxes from external forcing for each reservoir.
%
% - Topic: Diagnostics — Energy Fluxes — Spatial-temporal averages, [1 1]
% - Declaration: forcing_fluxes = quadraticEnergyFluxesSpatialTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns forcing_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
end

forcing_fluxes = self.filterFluxesForReservoir(self.quadraticEnergyFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices),filter=@(v) mean(v));
end