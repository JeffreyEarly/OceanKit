function fig = plotEnergyFluxOverTime(self,options)
% Plot energy fluxes as a function of time.
%
% You can plot either the exact or quadratic approximations to the energy fluxes.
%
% If you plot the quadratic fluxes, you can specify which energy reservoirs to include, which will create a subplot for each reservoir.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotEnergyFluxOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter approximation: (optional) {'exact','quadratic'} which approximation to use (default: 'exact')
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects to include (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter timeIndices: (optional) Indices of model times to plot/average (default: Inf -> all times)
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter filter: (optional) Function handle to apply to flux series before plotting (default: @(v) v)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.approximation {mustBeMember(options.approximation,{'quadratic','exact'})} = 'exact'
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v) v;
end
if options.approximation == "exact"
    [forcing_fluxes, t] = self.exactEnergyFluxesOverTime(timeIndices=options.timeIndices);

    fig = figure(Visible=options.visible);
    tl = tiledlayout(1,1,TileSpacing="compact");
    for iForce = 1:length(forcing_fluxes)
        plot(t/self.tscale,options.filter(forcing_fluxes(iForce).te/self.flux_scale)), hold on
    end
    legend(forcing_fluxes.fancyName)

    xlabel("time (" + self.tscale_units + ")")
    ylabel("flux (" + self.flux_scale_units + ")")
    xlim([min(t) max(t)]/self.tscale);
else
    [forcing_fluxes, t] = self.quadraticEnergyFluxesOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

    fig = figure(Visible=options.visible);
    tl = tiledlayout(length(options.energyReservoirs),1,TileSpacing="compact");
    for iReservoir = 1:length(options.energyReservoirs)
        nexttile(tl);
        for iForce = 1:length(forcing_fluxes)
            plot(t/self.tscale,options.filter(forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name)/self.flux_scale)), hold on
        end
        legend(forcing_fluxes.fancyName)

        fancyName = options.energyReservoirs(iReservoir).fancyName;
        xlabel("time (" + self.tscale_units + ")")
        ylabel("flux into " + fancyName + " (" + self.flux_scale_units + ")")
        xlim([min(t) max(t)]/self.tscale);
    end
end
end