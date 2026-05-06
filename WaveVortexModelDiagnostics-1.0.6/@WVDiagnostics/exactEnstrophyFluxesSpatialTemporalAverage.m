function enstrophy_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
% Compute spatial-temporal average of the exact enstrophy fluxes.
%
% Compute spatial-temporal average of the exact enstrophy fluxes
% Returns the spatial-temporal average of the exact enstrophy fluxes from external forcing
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — Spatial-temporal averages, [1 1]
% - Declaration: forcing_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns enstrophy_fluxes: diagnosed flux values
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

enstrophy_fluxes = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
for iForce = 1:length(enstrophy_fluxes)
    enstrophy_fluxes(iForce).Z0 = mean(enstrophy_fluxes(iForce).Z0);
end
end