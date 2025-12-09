function fig = plotEnstrophyFluxOverTime(self,options)
% Plot enstrophy fluxes for reservoirs as a function of time.
%
% Draws time series of enstrophy fluxes (exact or quadratic approximation)
% into specified enstrophy reservoirs. Supports selecting approximation,
% time indices, simple filtering for visualisation, and toggles for showing
% nonlinear advection, total flux, and dZ/dt.
%
% Some useful filters:
% filter=@(v,t) movmean(v,21);
% filter=@(v,t) cumtrapz(t,v)./(t+1)
%
% - Topic: Figures — Potential Enstrophy
% - Declaration: fig = plotEnstrophyFluxOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter approximation: (optional) {'quadratic','exact'} which approximation to use (default: 'exact')
% - Parameter timeIndices: (optional) Indices of model times to plot/average (default: Inf -> all times)
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter filter: (optional) Function handle to filter flux time series before plotting (default: @(v) v)
% - Parameter shouldShowNonlinearAdvection: (optional, logical) Show nonlinear advection term (default: true)
% - Parameter shouldShowTotal: (optional, logical) Show summed total flux (default: true)
% - Parameter shouldShowDtEnstrophy: (optional, logical) Show dZ/dt series for comparison (default: true)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.approximation {mustBeMember(options.approximation,{'quadratic','exact'})} = 'exact'
    options.timeIndices = Inf;
    options.visible = "on"
    options.filter = @(v) v;
    options.shouldShowNonlinearAdvection = true
    options.shouldShowTotal = true
    options.shouldShowDtEnstrophy = true
end
if options.approximation == "exact"
    [forcing_fluxes, t] = self.exactEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
    Z = self.exactEnstrophyOverTime(timeIndices=options.timeIndices);
else
    [forcing_fluxes, t] = self.quadraticEnstrophyFluxesOverTime(timeIndices=options.timeIndices);
    Z = self.quadraticEnstrophyOverTime(timeIndices=options.timeIndices);
end
if ~options.shouldShowNonlinearAdvection
    forcing_fluxes(1) = [];
end

fig = figure(Visible=options.visible);
tl = tiledlayout(1,1,TileSpacing="compact");
total = zeros(size(forcing_fluxes(1).Z0));
for iForce = 1:length(forcing_fluxes)
    plot(t/self.tscale,options.filter(forcing_fluxes(iForce).Z0/self.z_flux_scale)), hold on
    total = total + forcing_fluxes(iForce).Z0;
end
if options.shouldShowTotal
    plot(t/self.tscale,options.filter(total/self.z_flux_scale),Color=0*[1 1 1],LineWidth=2), hold on
    legendValues = cat(1,forcing_fluxes.fancyName,"total");
else
    legendValues = forcing_fluxes.fancyName;
end

if options.shouldShowDtEnstrophy
    t2 = t(2:end) - (t(2)-t(1))/2;
    dZdt = diff(Z)./diff(t);
    plot(t2/self.tscale,options.filter(dZdt/self.z_flux_scale),Color=0*[1 1 1],LineWidth=2,LineStyle="--"), hold on
    legendValues = cat(1,legendValues,"$\frac{d Z}{dt}$");
end

legend(legendValues);

xlabel("time (" + self.tscale_units + ")")
ylabel("flux (" + self.z_flux_scale_units + ")")
xlim([min(t) max(t)]/self.tscale);

if options.approximation == "exact"
    title(tl,"available potential enstrophy flux")
else
    title(tl,"quasigeostrophic potential enstrophy flux")
end
% mean(total/self.z_flux_scale)
end