function [enstrophy_fluxes,t] = quadraticEnstrophyFluxesOverTime(self,options)
% Compute enstrophy fluxes over time.
%
% Compute enstrophy fluxes over time
% Returns the enstrophy fluxes from external forcing
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — Time series, [t 1]
% - Declaration: enstrophy_fluxes = quadraticEnstrophyFluxesOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns enstrophy_fluxes: struct array with averaged fluxes
% - Returns t: Summary table of enstrophy flux diagnostics
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

if isinf(options.timeIndices)
    filter = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
enstrophy_fluxes = self.quadraticEnstrophyFluxes;
for iForce=1:length(enstrophy_fluxes)
    enstrophy_fluxes(iForce).Z0 = filter(enstrophy_fluxes(iForce).Z0);
end
t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end