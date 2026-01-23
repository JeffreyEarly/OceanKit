function forcing_fluxes = quadraticEnergyFluxesTemporalAverage(self,options)
% Compute temporally averaged forcing fluxes.
%
% Compute temporally averaged forcing fluxes
% Returns the temporally averaged energy fluxes from external forcing for each reservoir.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, [j kRadial]
% - Declaration: forcing_fluxes = quadraticEnergyFluxesTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns forcing_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
end

if isinf(options.timeIndices)
    filter_space = @(v) mean(v,3);
else
    filter_space = @(v) mean(v(:,:,options.timeIndices),3);
end
forcing_fluxes = self.filterFluxesForReservoir(self.quadraticEnergyFluxes(energyReservoirs=options.energyReservoirs),filter=filter_space);
end