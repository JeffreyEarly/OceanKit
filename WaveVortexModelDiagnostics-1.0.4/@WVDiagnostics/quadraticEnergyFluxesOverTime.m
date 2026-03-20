function [forcing_fluxes,t] = quadraticEnergyFluxesOverTime(self,options)
% Compute forcing fluxes over time.
%
% Compute forcing fluxes over time
% Returns the energy fluxes from external forcing for each reservoir as a function of time.
%
% - Topic: Diagnostics — Energy Fluxes — Time series, [t 1]
% - Declaration: forcing_fluxes = quadraticEnergyFluxesOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter timeIndices: (optional) indices for time selection (default: Inf)
% - Returns forcing_fluxes: struct array with fluxes over time
% - Returns t: Summary table of energy flux diagnostics
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
forcing_fluxes = self.quadraticEnergyFluxes(energyReservoirs=options.energyReservoirs);

forcing_fluxes = self.filterFluxesForReservoir(forcing_fluxes,filter=filter_space);
t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end