function [forcing_fluxes, t] = exactEnstrophyFluxesOverTime(self,options)
% Compute exact enstrophy fluxes over time.
%
% Compute exact enstrophy fluxes over time
% Returns the exact enstrophy fluxes from external forcing for each time step.
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — Time series, [t 1]
% - Declaration: forcing_fluxes = exactEnstrophyFluxesOverTime(self)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns forcing_fluxes: struct array with exact fluxes
% - Returns t: Summary table of enstrophy flux diagnostics
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_space = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter_space = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
forcing_fluxes = self.exactEnstrophyFluxes();
for iForce=1:length(forcing_fluxes)
    forcing_fluxes(iForce).Z0 = filter_space(forcing_fluxes(iForce).Z0);
end

t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end