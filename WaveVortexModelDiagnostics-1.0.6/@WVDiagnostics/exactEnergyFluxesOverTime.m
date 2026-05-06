function [energy_fluxes, t] = exactEnergyFluxesOverTime(self,options)
% Exact energy fluxes over time.
%
% Exact energy fluxes over time
% Returns the exact energy fluxes from external forcing for each time step.
%
% - Topic: Diagnostics — Energy Fluxes — Time series, [t 1]
% - Declaration: forcing_fluxes = exactEnergyFluxesOverTime(self)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns energy_fluxes: diagnosed flux values
% - Returns t: Summary table of energy flux diagnostics
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
energy_fluxes = self.exactEnergyFluxes();
for iForce=1:length(energy_fluxes)
    energy_fluxes(iForce).te = filter_space(energy_fluxes(iForce).te);
end

t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end