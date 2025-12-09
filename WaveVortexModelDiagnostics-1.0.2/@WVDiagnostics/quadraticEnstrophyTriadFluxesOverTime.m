function [enstrophy_fluxes,t] = quadraticEnstrophyTriadFluxesOverTime(self,options)
% Compute enstrophy inertial (aka, triad) fluxes over time.
%
% Compute enstrophy inertial (aka, triad) fluxes over time
% Returns the enstrophy fluxes from external forcing
%
% - Topic: Diagnostics — Enstrophy fluxes — Triad interactions — Fluxes over time, [t 1]
% - Declaration: enstrophy_fluxes = quadraticEnstrophyFluxesOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter triadComponents: (optional) input argument `triadComponents` (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Returns enstrophy_fluxes: struct array with averaged fluxes
% - Returns t: Summary table of enstrophy flux diagnostics
arguments
    self WVDiagnostics
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
end

if isinf(options.timeIndices)
    filter = @(v) reshape( sum(sum(v,1),2), [], 1);
else
    filter = @(v) reshape( sum(sum(v(:,:,options.timeIndices),1),2), [], 1);
end
enstrophy_fluxes = self.quadraticEnstrophyTriadFluxes(triadComponents=options.triadComponents);
for iForce=1:length(enstrophy_fluxes)
    enstrophy_fluxes(iForce).Z0 = filter(enstrophy_fluxes(iForce).Z0);
end
t = self.t_diag;
if ~isinf(options.timeIndices)
    t = t(options.timeIndices);
end
end