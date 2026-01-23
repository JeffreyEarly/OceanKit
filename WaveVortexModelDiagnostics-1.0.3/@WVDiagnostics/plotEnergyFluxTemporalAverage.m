function fig = plotEnergyFluxTemporalAverage(self,options)
% Plot temporally averaged energy flux diagnostics.
%
% Compute and plot temporally averaged energy flux diagnostics (exact or
% quadratic approximation) on a chosen axes projection. Supports grouping
% fluxes, overlaying frequency/ratio contours, and customizing colormap,
% saturation, and figure visibility.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotEnergyFluxTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter approximation: (optional) {'quadratic','exact'} approximation to use (default: 'quadratic')
% - Parameter energyReservoir: (optional) EnergyReservoir to plot (default: EnergyReservoir.total)
% - Parameter triadComponents: (optional) TriadFlowComponent vector for inertial flux selection (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Parameter showForcingFluxes: (optional) (optional, logical) Include forcing fluxes when using quadratic approximation (default: true)
% - Parameter timeIndices: (optional) Time indices to average over (default: Inf -> all times)
% - Parameter axes: (optional) Plot projection, one of {'jk','j','jWavenumber','k','k-pseudo-isotropic','omega'} (default: 'jk')
% - Parameter filter: (optional) Function handle applied to plotted data (default: @(v) v)
% - Parameter shouldOverlayWaveFrequencies: (optional) (optional, logical) Overlay wave-frequency contours on 'jk' axes (default: false)
% - Parameter shouldOverlayGeostrophicKineticPotentialRatioContours: (optional) (optional, logical) Overlay geostrophic KE/PE fraction contours (default: true)
% - Parameter colormap: (optional) Colormap to use for 'jk' axes (default: WVDiagnostics.crameri('-bam'))
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter overSaturationFactor: (optional) Scalar or two-element colormap limits (default: 10)
% - Parameter fluxGroups: (optional) Cell array of index groups to combine for plotting (default: [])
% - Parameter simpleName: (optional) Cell array of simple names for fluxGroups (default: [])
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.approximation {mustBeMember(options.approximation,{'quadratic','exact'})} = 'quadratic'
    options.energyReservoir = EnergyReservoir.total;
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
    options.showForcingFluxes = true;
    options.timeIndices = Inf;
    options.axes {mustBeMember(options.axes,{'jk','j','jWavenumber','k','k-pseudo-isotropic','omega'})} = 'jk'
    options.filter = @(v) v;
    options.shouldOverlayWaveFrequencies = false
    options.shouldOverlayGeostrophicKineticPotentialRatioContours = true
    options.colormap = WVDiagnostics.crameri('-bam')
    options.visible = "on"
    options.overSaturationFactor = 10;
    options.fluxGroups = [];
    options.simpleName = [];
end

if ~any( options.energyReservoir == [EnergyReservoir.geostrophic_mda,EnergyReservoir.wave,EnergyReservoir.total] )
    error('The energy reservoir must either geostrophic_mda, wave, or total.');
end

wvt = self.wvt;
if options.approximation == "exact"
    fluxes = self.exactEnergyFluxesTemporalAverage(timeIndices=options.timeIndices);
else
    forcing_fluxes = self.quadraticEnergyFluxesTemporalAverage(energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
    inertial_fluxes = self.quadraticEnergyTriadFluxesTemporalAverage(triadComponents=options.triadComponents,energyReservoirs=options.energyReservoir,timeIndices=options.timeIndices);
    if options.showForcingFluxes
        fluxes = cat(2,forcing_fluxes,inertial_fluxes);
    else
        fluxes = inertial_fluxes;
    end
end

% Custom combination and ordering of fluxes for plotting
% Use fluxes.name to see all fluxes, or plot them all un-grouped. Specify
% indices of fluxes to group.
% General idea: all forcing, dominant triad, all other triads, all damping 
if ~isempty(options.fluxGroups)
    fluxesGrouped = struct;
    for g = 1:numel(options.fluxGroups)
        fluxesGrouped(g).name = cat(2,fluxes(options.fluxGroups{g}).name);
        fluxesGrouped(g).fancyName = cat(2,fluxes(options.fluxGroups{g}).fancyName);
        fluxesGrouped(g).simpleName = options.simpleName{g};
        fluxesGrouped(g).(options.energyReservoir.name) = sum( cat(3, fluxes(options.fluxGroups{g}).(options.energyReservoir.name)), 3);
    end
    % replace fluxes with fluxesGrouped
    fluxes = fluxesGrouped;
end

% set color limits 
if options.axes == "jk"
    if isscalar(options.overSaturationFactor)
        if options.approximation == "exact"
            colorLimits = max(arrayfun( @(v) max(abs(v.te(:))), fluxes))*[-1 1]/options.overSaturationFactor;            
        else
            colorLimits = max(arrayfun( @(v) max(abs(v.(options.energyReservoir.name)(:))), fluxes))*[-1 1]/options.overSaturationFactor;
        end
        colorLimits = colorLimits/self.flux_scale;
    elseif length(options.overSaturationFactor)==2
        colorLimits = options.overSaturationFactor;
    else
        error("options.overSaturationFactor should be length 1 (overSaturationFactor) or length 2 (colormap limits).")
    end
end

% create radial wavelength vector
switch options.axes
    case "k-pseudo-isotropic"
        radialWavelength = 2*pi./self.kPseudoRadial/1000;
    otherwise
        radialWavelength = 2*pi./self.kRadial/1000;
end
radialWavelength(1) = 1.5*radialWavelength(2);

% create j vector for log y-axis.
jForLogAxis = self.j;
jForLogAxis(1) = 0.75;

fig = figure(Visible=options.visible);
tl = tiledlayout(fig,"flow",TileSpacing='tight');

for iComponent = 1:length(fluxes)
    if options.approximation == "exact"
        val = fluxes(iComponent).te/self.flux_scale;
    else
        val = fluxes(iComponent).(options.energyReservoir.name)/self.flux_scale;
    end
    ax = nexttile(tl);
    switch options.axes
        case "jk"
            % % % pcolor(ax,options.energyReservoir.kFromKRadial(wvt.kRadial),wvt.j,val), shading flat
            % % % self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
            pcolor(ax,radialWavelength,jForLogAxis,options.filter(val)), shading flat
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            colormap(ax, options.colormap)
            if options.approximation == "exact"
            else
                if options.energyReservoir==EnergyReservoir.total
                    text(radialWavelength(1)*.95,jForLogAxis(1),{'MDA','Inertial'},'FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
                elseif options.energyReservoir==EnergyReservoir.geostrophic_mda
                    text(radialWavelength(1)*.95,jForLogAxis(1),'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
                elseif options.energyReservoir==EnergyReservoir.wave
                    text(radialWavelength(1)*.95,jForLogAxis(1),'Inertial','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
                end
            end
            line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)           
            if options.shouldOverlayWaveFrequencies
                self.overlayFrequencyContours;
            end
            if options.shouldOverlayGeostrophicKineticPotentialRatioContours
                self.overlayGeostrophicKineticPotentialFractionContours;
            end
            clim(ax,colorLimits);
            
        case "j"
            plot(self.j,zeros(size(self.j)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,self.j,options.filter(sum(val,2)))
        case "jWavenumber"
            jWavelength = 2*pi./self.jWavenumber/1000;
            jWavelength(1) = 1.5*jWavelength(2);
            plot(jWavelength,zeros(size(self.j)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,jWavelength,options.filter(sum(val,2)))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            xlim([jWavelength(end) jWavelength(1)])
        case "k"
            % % % plot(options.energyReservoir.kFromKRadial(wvt.kRadial),zeros(size(wvt.kRadial)),LineWidth=2,Color=0*[1 1 1]), hold on
            % % % plot(ax,wvt.kRadial,sum(val,1))
            % % % self.setLogWavelengthXAxis(num_ticks=6,roundToNearest=5)
            plot(radialWavelength,zeros(size(radialWavelength)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,radialWavelength,options.filter(sum(val,1)))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
        case "k-pseudo-isotropic"
            v = self.transformToPseudoRadialWavenumber(options.energyReservoir,val);
            plot(radialWavelength,zeros(size(radialWavelength)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,radialWavelength,options.filter(v))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
        case "omega"
            omegaAxis = self.omegaAxis/self.wvt.f;
            v = self.transformToOmegaAxis(val);
            plot(omegaAxis,zeros(size(omegaAxis)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,omegaAxis,options.filter(v))
            set(gca,'XScale','log')
    end
    if isempty(options.simpleName)
        title(ax,fluxes(iComponent).fancyName + " (" + string(sum(val(:))) + " " + self.flux_scale_units + ")" )
    else
        title(ax,fluxes(iComponent).simpleName + " (" + string(sum(val(:))) + " " + self.flux_scale_units + ")" )
    end
end

switch options.axes
    case "jk"
        % jkR plot labels
        xlabel(tl,'horizontal wavelength (km)')
        ylabel(tl,"vertical mode")
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = "energy flux (" + self.flux_scale_units + ")";
    case "j"
        % j plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl, 'vertical mode')
    case "jWavenumber"
        % kR plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl,'radius of deformation (km)')
    case "k"
        % kR plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl,'horizontal wavelength (km)')
    case "k-pseudo-isotropic"
        % kR plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl,'pseudo-wavelength (km)')
    case "omega"
        % kR plot options
        ylabel(tl,"energy flux (" + self.flux_scale_units + ")")
        xlabel(tl,'frequency (f)')
end

switch options.axes
    case "jk"
        for ii=1:tilenum(ax) %prod(tl_jkR.GridSize)
            if ii <= prod(tl.GridSize)-tl.GridSize(2)
                nexttile(ii)
                %xticklabels([])
                % set(gca,'XTickLabels',[],'FontSize',n_size)
                set(gca,'XTickLabels',[])
            end
            if ~(mod(ii,tl.GridSize(2))==1)
                nexttile(ii)
                %yticklabels([])
                % set(gca,'YTickLabels',[],'FontSize',n_size)
                set(gca,'YTickLabels',[])
            end
            if (mod(ii,tl.GridSize(2))==1)
                nexttile(ii)
                % set(gca,'YTick',[0 5 10 15])
                % set(gca,'YTickLabels',get(gca,'ytick'),'FontSize',n_size)
                set(gca,'YTickLabels',get(gca,'ytick'))
                % set(gca,'fontname','times')
            end
        end
end


if isinf(options.timeIndices)
    options.timeIndices = 1:length(self.t_diag);
end
minDay = string(round(self.t_diag(min(options.timeIndices))/86400));
maxDay = string(round(self.t_diag(max(options.timeIndices))/86400));
if options.approximation == "exact"
    title(tl,"exact energy flux, day " + minDay + "-" + maxDay)
else
    title(tl,"quadratic energy flux into " + options.energyReservoir.fancyName + " energy, day " + minDay + "-" + maxDay)
end

end