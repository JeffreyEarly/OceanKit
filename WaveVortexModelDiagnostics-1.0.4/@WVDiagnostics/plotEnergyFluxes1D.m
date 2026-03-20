function fig = plotEnergyFluxes1D(self,options)
% Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency.
%
% Plot 1D energy fluxes as a function of pseudo-wavelength / KE-PE / frequency
% Produces one-dimensional summaries of energy fluxes (forcing and inertial)
% projected onto a radial (pseudo-)wavenumber axis, KE/(KE+APE) axis, and
% frequency axis. Computes temporal averages over the requested timeIndices,
% prepares per-flux metadata (colors, styles, alpha), and draws a compact
% 2x2 multipanel figure that highlights triad and forcing contributions.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotEnergyFluxes1D(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) Indices of model times to average/plot (default: Inf -> all times)
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter triadLineWidth: (optional) Line width for triad/inertial flux lines (default: 2)
% - Parameter forcingLineWidth: (optional) Line width for forcing flux lines (default: 1.5)
% - Parameter fluxTolerance: (optional) Minimum absolute flux magnitude to display (default: 5e-2)
% - Parameter keFractions: (optional) Fractional ticks for KE/(KE+APE) axis (default: [.01,0.025,.1,.25,.5,.75,.9,0.975,.99])
% - Parameter forcingFluxAttributes: Array/struct to override forcing attributes (color, fancyName, alpha)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.timeIndices = Inf
    options.visible = "on"
    options.triadLineWidth = 2
    options.forcingLineWidth = 1.5
    options.fluxTolerance = 5e-2
    options.keFractions = [.01,0.025,.1,.25,.5,.75,.9,0.975,.99]
    options.forcingFluxAttributes
end

% yLimits = [-2.2,2.2]; % set automatically below


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data: inertial fluxes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[inertial_fluxes_g, inertial_fluxes_w, kp] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(timeIndices=options.timeIndices);
[inertial_fluxes_w_omega, omegaAxisSparse] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega(timeIndices=options.timeIndices);
[inertial_fluxes_g_kepe, kePeAxis] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe(timeIndices=options.timeIndices);

fluxNames = ["ggg", "ggw", "tx-wwg", "tx-ggw"];
fluxLineStyle = ["-", ":", "--", "-."];
for iName=1:length(fluxNames)
    inertial_fluxes_g([inertial_fluxes_g.name] == fluxNames(iName)).color=0*[1 1 1];
    inertial_fluxes_g([inertial_fluxes_g.name] == fluxNames(iName)).lineStyle=fluxLineStyle(iName);
    inertial_fluxes_g([inertial_fluxes_g.name] == fluxNames(iName)).flux_kepe = inertial_fluxes_g_kepe([inertial_fluxes_g.name] == fluxNames(iName)).flux;
end

fluxNames = ["www", "wwg", "tx-wwg", "tx-ggw"];
fluxLineStyle = ["-", ":", "--", "-."];
for iName=1:length(fluxNames)
    inertial_fluxes_w([inertial_fluxes_w.name] == fluxNames(iName)).color=0*[1 1 1];
    inertial_fluxes_w([inertial_fluxes_w.name] == fluxNames(iName)).lineStyle=fluxLineStyle(iName);
    inertial_fluxes_w([inertial_fluxes_w.name] == fluxNames(iName)).flux_omega = inertial_fluxes_w_omega([inertial_fluxes_w.name] == fluxNames(iName)).flux;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data: energy spectrum rate of change
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wvt = self.wvt;

% initial
self.iTime=options.timeIndices(1);
TE_A0_j_kR_initial = wvt.transformToRadialWavenumber(wvt.A0_TE_factor .* abs(self.wvt.A0).^2);
TE_Apm_j_kR_initial = wvt.transformToRadialWavenumber(wvt.Apm_TE_factor .* (abs(self.wvt.Ap).^2 + abs(self.wvt.Am).^2));

% final
self.iTime=options.timeIndices(end);
TE_A0_j_kR_final = wvt.transformToRadialWavenumber(wvt.A0_TE_factor .* abs(self.wvt.A0).^2);
TE_Apm_j_kR_final = wvt.transformToRadialWavenumber(wvt.Apm_TE_factor .* (abs(self.wvt.Ap).^2 + abs(self.wvt.Am).^2));

% change over time
ddt_TE_A0 = (TE_A0_j_kR_final - TE_A0_j_kR_initial)/(self.t_wv(options.timeIndices(end)) - self.t_wv(options.timeIndices(1)));
ddt_TE_Apm = (TE_Apm_j_kR_final - TE_Apm_j_kR_initial)/(self.t_wv(options.timeIndices(end)) - self.t_wv(options.timeIndices(1)));

% zero-pad ddt_TE_A* to match size self matrices with anti-aliased modes
ddt_TE_A0 = self.transformToPseudoRadialWavenumber(EnergyReservoir.geostrophic_mda,paddata(ddt_TE_A0,[length(self.j),length(self.kRadial)],FillValue=0,Side='trailing'));
ddt_TE_Apm = self.transformToPseudoRadialWavenumber(EnergyReservoir.wave,paddata(ddt_TE_Apm,[length(self.j),length(self.kRadial)],FillValue=0,Side='trailing'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data: forcing fluxes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

energy_fluxes = self.quadraticEnergyFluxesTemporalAverage(energyReservoirs=[EnergyReservoir.geostrophic_mda, EnergyReservoir.wave],timeIndices=options.timeIndices);
energy_fluxes(1) = []; % dump the nonlinear advection

% assign default colors, alpha, fancy name to the forcings
C = orderedcolors("gem"); 
for i=1:length(energy_fluxes)
    forcing_fluxes_g(i).name = energy_fluxes(i).name;
    forcing_fluxes_g(i).fancyName = energy_fluxes(i).fancyName;
    forcing_fluxes_g(i).color = C(mod(i,size(C,1))+1,:);
    forcing_fluxes_g(i).alpha = 0.3;
    forcing_fluxes_g(i).flux = self.transformToPseudoRadialWavenumber(EnergyReservoir.geostrophic_mda,energy_fluxes(i).te_gmda);
    forcing_fluxes_g(i).flux_kepe = self.transformToKePeAxis(energy_fluxes(i).te_gmda);

    forcing_fluxes_w(i).name = energy_fluxes(i).name;
    forcing_fluxes_w(i).fancyName = energy_fluxes(i).fancyName;
    forcing_fluxes_w(i).color = C(mod(i,size(C,1))+1,:);
    forcing_fluxes_w(i).alpha = 0.3;
    forcing_fluxes_w(i).flux = self.transformToPseudoRadialWavenumber(EnergyReservoir.wave,energy_fluxes(i).te_wave);
    forcing_fluxes_w(i).flux_omega = self.transformToOmegaAxis(energy_fluxes(i).te_wave);
end

% override the forcings, if the user specifies
if isfield(options,"forcingFluxAttributes")
    forcing_fluxes_g_old = forcing_fluxes_g;
    forcing_fluxes_w_old = forcing_fluxes_w;
    clear forcing_fluxes_g;
    clear forcing_fluxes_w;
    for i=1:length(options.forcingFluxAttributes)
        idx = [forcing_fluxes_w_old.name] == options.forcingFluxAttributes(i).name;
        if any(idx)
            forcing_fluxes_g(i) = forcing_fluxes_g_old(idx);
            forcing_fluxes_w(i) = forcing_fluxes_w_old(idx);
            fields = ["fancyName","color","alpha"];
            for iField=1:length(fields)
                if isfield(options.forcingFluxAttributes(i),fields(iField))
                    forcing_fluxes_g(i).(fields(iField)) = options.forcingFluxAttributes(i).(fields(iField));
                    forcing_fluxes_w(i).(fields(iField)) = options.forcingFluxAttributes(i).(fields(iField));
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Visualization: figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radialWavelengthSparse = 2*pi./kp/1000;
radialWavelengthSparse(1) = 1.5*radialWavelengthSparse(2);
radialWavelength = 2*pi./self.kPseudoRadial/1000;
radialWavelength(1) = radialWavelengthSparse(1); %1.5*radialWavelength(2);

filter = @(v) cumsum(v)/self.flux_scale;

omegaAxis = self.omegaAxis/wvt.f;
omegaAxis = [0.9,omegaAxis]; % add point below f

omegaAxisSparse = omegaAxisSparse/wvt.f;
omegaAxisSparse = [0.9;omegaAxisSparse]; % add point below f
filter_omega = @(v) cumsum([0;v])/self.flux_scale;

fig = figure(Visible=options.visible,Units="points",Position=[50 50 800 500]);
fig.PaperPositionMode = "auto";
fig.Color = "w";
tl = tiledlayout(fig,GridSize=[2 2],TileSpacing="compact",Padding="compact");

% max and min fluxes for ylimits
ffgMax = max(max(abs(filter([forcing_fluxes_g.flux]))));
ffwMax = max(max(abs(filter([forcing_fluxes_w.flux]))));
ifgMax = max(max(abs(filter([inertial_fluxes_g.flux]))));
ifwMax = max(max(abs(filter([inertial_fluxes_w.flux]))));
yLimits = [-max([ffgMax ffwMax ifgMax ifwMax]),max([ffgMax ffwMax ifgMax ifwMax])];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization: geostrophic fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = nexttile(tl,1);
fluxes = inertial_fluxes_g;
for i=1:length(fluxes)
    if isempty(filter(fluxes(i).flux))
        continue;
    elseif max(abs(filter(fluxes(i).flux))) < options.fluxTolerance
        continue;
    end
    plot(radialWavelengthSparse,filter(fluxes(i).flux),LineWidth=options.triadLineWidth,Color=fluxes(i).color,LineStyle=fluxes(i).lineStyle,DisplayName=fluxes(i).fancyName), hold on
end
fluxes = forcing_fluxes_g;
for i=1:length(fluxes)
    flux = filter(fluxes(i).flux);
    if isempty(flux)
        continue;
    elseif max(abs(flux)) < options.fluxTolerance
        continue;
    end
    plot(radialWavelength,flux,LineWidth=options.forcingLineWidth,Color=fluxes(i).color,DisplayName=fluxes(i).fancyName), hold on
    drawPatchForFluxWithColor(radialWavelength,flux,fluxes(i).color);
end
plot(ax,radialWavelength,filter(ddt_TE_A0),'c',LineWidth=options.forcingLineWidth,DisplayName='d/dt geostrophic energy')


ax.XScale = "log";
ax.XDir = "reverse";
% ax.XTickLabels = [];
ax.XLabel.String = "pseudo-wavelength (km)";
ax.XLim = [min(radialWavelengthSparse) max(radialWavelengthSparse)];
ax.YLim = yLimits;
ax.YLabel.String = "geostrophic energy flux (" + self.flux_scale_units + ")";

lgd1 = legend(Location="southwest",Interpreter="latex");
lgd1.NumColumns = 2;

text(ax,max(xlim)*.95,max(ylim),'a)','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization: geostrophic fluxes KePe axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kePeAxisSimple = reshape(1:length(kePeAxis),size(kePeAxis));

ax = nexttile(tl,2);
fluxes = inertial_fluxes_g;
for i=1:length(fluxes)
    if isempty(filter(fluxes(i).flux))
        continue;
    elseif max(abs(filter(fluxes(i).flux))) < options.fluxTolerance
        continue;
    end
    plot(kePeAxisSimple,filter(fluxes(i).flux_kepe),LineWidth=options.triadLineWidth,Color=fluxes(i).color,LineStyle=fluxes(i).lineStyle,DisplayName=fluxes(i).fancyName), hold on
end
fluxes = forcing_fluxes_g;
for i=1:length(fluxes)
    flux = filter(fluxes(i).flux_kepe);
    if isempty(flux)
        continue;
    elseif max(abs(flux)) < options.fluxTolerance
        continue;
    end
    plot(kePeAxisSimple,flux,LineWidth=options.forcingLineWidth,Color=fluxes(i).color,DisplayName=fluxes(i).fancyName), hold on
    drawPatchForFluxWithColor(kePeAxisSimple,flux,fluxes(i).color);
end
% plot(ax,radialWavelength,filter(ddt_TE_Apm),'c',LineWidth=options.forcingLineWidth,DisplayName='d/dt wave energy')

ax.XTick = interp1(kePeAxis,kePeAxisSimple,options.keFractions);
ax.XTickLabels = string(options.keFractions);
ax.XLim = [kePeAxisSimple(1) kePeAxisSimple(end)];
ax.YLim = yLimits;
ax.XLabel.String = "KE/(KE+APE)";
ax.YTickLabels = [];

text(ax,min(xlim)*1.05,max(ylim),'b)','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization: wave fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = nexttile(tl,3);
fluxes = inertial_fluxes_w;
for i=1:length(fluxes)
    if isempty(filter(fluxes(i).flux))
        continue;
    elseif max(abs(filter(fluxes(i).flux))) < options.fluxTolerance
        continue;
    end
    plot(radialWavelengthSparse,filter(fluxes(i).flux),LineWidth=options.triadLineWidth,Color=fluxes(i).color,LineStyle=fluxes(i).lineStyle,DisplayName=fluxes(i).fancyName), hold on
end
fluxes = forcing_fluxes_w;
for i=1:length(fluxes)
    flux = filter(fluxes(i).flux);
    if isempty(flux)
        continue;
    elseif max(abs(flux)) < options.fluxTolerance
        continue;
    end
    plot(radialWavelength,flux,LineWidth=options.forcingLineWidth,Color=fluxes(i).color,DisplayName=fluxes(i).fancyName), hold on
    drawPatchForFluxWithColor(radialWavelength,flux,fluxes(i).color);
end
plot(ax,radialWavelength,filter(ddt_TE_Apm),'c',LineWidth=options.forcingLineWidth,DisplayName='d/dt wave energy')

ax.XScale = "log";
ax.XDir = "reverse";
ax.XLim = [min(radialWavelengthSparse) max(radialWavelengthSparse)];
ax.YLim = yLimits;
ax.XLabel.String = "pseudo-wavelength (km)";
ax.YLabel.String = "wave energy flux (" + self.flux_scale_units + ")";

lgd2 = legend(Location="southwest",Interpreter="latex");
lgd2.NumColumns = 2;

text(ax,max(xlim)*.95,max(ylim),'c)','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization: wave fluxes frequency axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = nexttile(tl,4);
fluxes = inertial_fluxes_w;
for i=1:length(fluxes)
    if isempty(filter(fluxes(i).flux))
        continue;
    elseif max(abs(filter(fluxes(i).flux))) < options.fluxTolerance
        continue;
    end
    plot(omegaAxisSparse,filter_omega(fluxes(i).flux_omega),LineWidth=options.triadLineWidth,Color=fluxes(i).color,LineStyle=fluxes(i).lineStyle,DisplayName=fluxes(i).fancyName), hold on
end
fluxes = forcing_fluxes_w;
for i=1:length(fluxes)
    flux = filter_omega(fluxes(i).flux_omega);
    if isempty(flux)
        continue;
    elseif max(abs(flux)) < options.fluxTolerance
        continue;
    end
    plot(omegaAxis,flux,LineWidth=options.forcingLineWidth,Color=fluxes(i).color,DisplayName=fluxes(i).fancyName), hold on
    drawPatchForFluxWithColor(omegaAxis,flux,fluxes(i).color);
end
% plot(ax,radialWavelength,filter(ddt_TE_Apm),'c',LineWidth=options.forcingLineWidth,DisplayName='d/dt wave energy')

ax.XScale = "log";
ax.XLim = [min(omegaAxis) max(omegaAxis)];
ax.YLim = yLimits;
ax.XLabel.String = "frequency (f)";
ax.YTickLabels = [];

text(ax,min(xlim)*1.05,max(ylim),'d)','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top')

% setTileSpacingBetween(tl, 0.5)

% % get all YLim and set all to match largest
% % This would work, except drawPatchForFluxWithColor above also uses YLim
% axesHandles = findobj(fig.Children, 'Type', 'axes');
% % Get all Y data from each axes
% allYLim = [];
% for iax = 1:length(axesHandles)
%     allYLim = [allYLim axesHandles(iax).YLim];
% end
% % Compute common Y limits
% yLimits = [min(allYLim), max(allYLim)];
% % Set YLim for all axes
% for ax = axesHandles'
%     ax.YLim = yLimits;
% end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to highlight a flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function drawPatchForFluxWithColor(fluxAxis,flux,color)
        % find radialWavelength for 10%-90% forcing interval
        cv = abs(flux) / max(abs(flux));
        if all(isfinite(cv)) == false
            return
        end
        % Remove duplicates values
        uniqueInd = [false; diff(cv) ~= 0]; % keep only points where cv changes
        firstOne = find(uniqueInd ~= 0, 1, 'first');
        if ~isempty(firstOne) & firstOne~=1
            uniqueInd(firstOne-1) = 1;
        end
        cv_unique = cv(uniqueInd);
        r_unique = fluxAxis(uniqueInd);
        r01 = interp1(cv_unique, r_unique, 0.16);
        r09 = interp1(cv_unique, r_unique, 0.84);
        % patch coordinates
        patchX = [r09,r01,r01,r09];
        if flux(end)>0
            patchY = [0,0,yLimits(2),yLimits(2)];
        else
            patchY = [yLimits(1),yLimits(1),0,0];
        end
        % add forcing/damping patch
        patch(ax,patchX,patchY,color,'FaceAlpha',fluxes(i).alpha,'EdgeColor','none','HandleVisibility','off')
    end



end