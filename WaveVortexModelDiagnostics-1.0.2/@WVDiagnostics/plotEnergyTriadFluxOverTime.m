function fig = plotEnergyTriadFluxOverTime(self,options)
% Plot inertial flux for each reservoir over time
%
% Plots the energy flux between reservoirs due to inertial interactions as a function of time.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotInertialFluxOverTime(self,options)
% - Parameter energyReservoirs: vector of EnergyReservoir objects (default: [geostrophic, wave, total])
% - Parameter visible: figure visibility (default: "on")
% - Parameter filter: function handle to filter fluxes (default: @(v) v)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave];
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v) v;
end
[inertial_fluxes,t] = self.quadraticEnergyTriadFluxesOverTime(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

fig = figure(Visible=options.visible);
tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
for iReservoir = 1:length(options.energyReservoirs)
    nexttile(tl);
    for iForce = 1:length(inertial_fluxes)
        plot(t/self.tscale,options.filter(inertial_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/self.flux_scale)), hold on
    end
    legend(inertial_fluxes.fancyName)

    fancyName = options.energyReservoirs(iReservoir).fancyName;
    xlabel("time (" + self.tscale_units + ")")
    ylabel("flux into " + fancyName + " (" + self.flux_scale_units + ")")
    xlim([min(t) max(t)]/self.tscale);
end
end