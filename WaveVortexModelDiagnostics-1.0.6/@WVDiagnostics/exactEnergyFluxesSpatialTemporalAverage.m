function forcing_fluxes = exactEnergyFluxesSpatialTemporalAverage(self,options)
% Compute spatial-temporal average of the exact forcing fluxes.
%
% Compute spatial-temporal average of the exact forcing fluxes
% Returns the spatial-temporal average of the exact energy fluxes from external forcing
%
% - Topic: Diagnostics — Energy Fluxes — Spatial-temporal averages, [1 1]
% - Declaration: forcing_fluxes = exactEnergyFluxesSpatialTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns forcing_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

forcing_fluxes = self.exactEnergyFluxesOverTime(timeIndices=options.timeIndices);
for iForce = 1:length(forcing_fluxes)
    forcing_fluxes(iForce).te = mean(forcing_fluxes(iForce).te);
end
end