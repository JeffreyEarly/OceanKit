function fig = plotEnstrophyTriadFluxOverTime(self,options)
% Plot enstrophy triad fluxes over time.
%
% Plot enstrophy triad fluxes over time.
% Plots time series of enstrophy triad fluxes computed with the quadratic
% approximation for the specified triad components. An optional filter may
% be applied to each time series before plotting. Axis labels use the class
% scaling properties (tscale, tscale_units, z_flux_scale, z_flux_scale_units).
%
% - Topic: Figures — Potential Enstrophy
% - Declaration: fig = plotEnstrophyTriadFluxOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter triadComponents: (optional) vector of TriadFlowComponent objects to include (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Parameter timeIndices: (optional) Indices of model times to use (default: Inf -> all times)
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter filter: (optional) Function handle accepting (v,t) used to preprocess each flux series before plotting (default: @(v,t) v)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v,t) v;
end
[forcing_fluxes, t] = self.quadraticEnstrophyTriadFluxesOverTime(timeIndices=options.timeIndices,triadComponents=options.triadComponents);

fig = figure(Visible=options.visible);
tl = tiledlayout(1,1,TileSpacing="compact");

for iForce = 1:length(forcing_fluxes)
    plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale,t)), hold on
end
legend(forcing_fluxes.fancyName)

xlabel("time (" + self.tscale_units + ")")
ylabel("enstrophy flux (" + self.z_flux_scale_units + ")")
xlim([min(t) max(t)]/self.tscale);

end