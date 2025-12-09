function fig = plotEnergyOverTime(self,options)
% Plot energy for each reservoir over time.
%
% Plots the energy in each specified reservoir as a function of time.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotEnergyOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter shouldIncludeExactTotalEnergy: (optional) include exact total energy (default: true)
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Parameter visible: (optional) figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total];
    options.shouldIncludeExactTotalEnergy = true
    options.timeIndices = Inf;
    options.visible = "on"
end
[reservoirs, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,shouldIncludeExactTotalEnergy=options.shouldIncludeExactTotalEnergy,timeIndices=options.timeIndices);

fig = figure(Visible=options.visible);

for iReservoir = 1:length(reservoirs)
    switch reservoirs(iReservoir).name
        case "te"
            plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2, Color=[0 0 0]), hold on
        case "te_quadratic"
            plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2, Color=[0 0 0], LineStyle="-."), hold on
        otherwise
            plot(t/self.tscale,reservoirs(iReservoir).energy/self.escale,LineWidth=2), hold on
    end
end
legend(reservoirs.fancyName)

xlabel("time (" + self.tscale_units + ")")
ylabel("energy (" + self.escale_units + ")")
xlim([min(t) max(t)]/self.tscale);
end