function enstrophy_fluxes = exactEnstrophyFluxesTemporalAverage(self,options)
% Compute temporally averaged enstrophy fluxes.
%
% Compute temporally averaged enstrophy fluxes
% Returns the temporally averaged enstrophy fluxes from external forcing for each reservoir.
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — Temporal averages, [j kRadial]
% - Declaration: enstrophy_fluxes = exactEnstrophyFluxesTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns enstrophy_fluxes: struct array with averaged fluxes
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

if isinf(options.timeIndices)
    filter_space = @(v) mean(v,3);
else
    filter_space = @(v) mean(v(:,:,options.timeIndices),3);
end
enstrophy_fluxes = self.exactEnstrophyFluxes;
for iForce=1:length(enstrophy_fluxes)
    enstrophy_fluxes(iForce).Z0 = filter_space(enstrophy_fluxes(iForce).Z0);
end
end