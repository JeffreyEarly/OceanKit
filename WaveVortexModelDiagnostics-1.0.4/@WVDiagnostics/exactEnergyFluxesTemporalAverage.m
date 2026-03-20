function forcing_fluxes = exactEnergyFluxesTemporalAverage(self,options)
% Temporally averaged exact energy fluxes.
%
% Temporally averaged exact energy fluxes
% Returns the temporally averaged enstrophy fluxes from external forcing for each reservoir.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, [j kRadial]
% - Declaration: enstrophy_fluxes = exactEnergyFluxesTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns forcing_fluxes: diagnosed flux values
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

if isinf(options.timeIndices)
    filter_space = @(v) mean(v,3);
else
    filter_space = @(v) mean(v(:,:,options.timeIndices),3);
end
forcing_fluxes = self.exactEnergyFluxes;
for iForce=1:length(forcing_fluxes)
    forcing_fluxes(iForce).te = filter_space(forcing_fluxes(iForce).te);
end
end