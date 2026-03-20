function fig = plotEnstrophyFluxTemporalAverage(self,options)
% Plot temporally averaged enstrophy flux diagnostics.
%
% Compute and plot temporally averaged enstrophy flux diagnostics (exact or
% quadratic approximation) on a chosen projection. Supports 'jk' (j vs
% radial wavelength) pcolor plots, 1D summaries along j or k axes, and a
% pseudo-isotropic radial projection. Allows applying a filter to the data,
% customizing the colormap and saturation, and controlling figure visibility.
%
% - Topic: Figures — Potential Enstrophy
% - Declaration: fig = plotEnstrophyFluxTemporalAverage(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) Time indices to average over (default: Inf -> all times)
% - Parameter approximation: (optional) {'quadratic','exact'} approximation to use (default: 'exact')
% - Parameter axes: (optional) Plot projection; one of {'jk','j','jWavenumber','k','k-pseudo-isotropic'} (default: 'jk')
% - Parameter filter: (optional) Function handle applied to plotted data (default: @(v) v)
% - Parameter colormap: (optional) Colormap for 'jk' axes (default: WVDiagnostics.crameri('-bam'))
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter overSaturationFactor: (optional) Scalar or two-element colormap limits (default: 10)
% - Parameter simpleName: (optional) Cell array of simple name strings for flux panels (default: [])
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
    options.approximation {mustBeMember(options.approximation,{'quadratic','exact'})} = 'exact'
    options.axes {mustBeMember(options.axes,{'jk','j','jWavenumber','k','k-pseudo-isotropic'})} = 'jk'
    options.filter = @(v) v;
    options.colormap = WVDiagnostics.crameri('-bam')
    options.visible = "on"
    options.overSaturationFactor = 10;
    options.simpleName = [];
end

if options.approximation == "exact"
    fluxes = self.exactEnstrophyFluxesTemporalAverage(timeIndices=options.timeIndices);
else
    fluxes = self.quadraticEnstrophyFluxesTemporalAverage(timeIndices=options.timeIndices);
end

% set color limits 
if options.axes == "jk"
    if isscalar(options.overSaturationFactor)
        colorLimits = max(arrayfun( @(v) max(abs(v.Z0(:))), fluxes))*[-1 1]/options.overSaturationFactor;
        colorLimits = colorLimits/self.z_flux_scale;
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
    val = fluxes(iComponent).Z0/self.z_flux_scale;
    ax = nexttile(tl);
    switch options.axes
        case "jk"
            pcolor(ax,radialWavelength,jForLogAxis,options.filter(val)), shading flat
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            colormap(ax, options.colormap)
            line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)           

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
            plot(radialWavelength,zeros(size(radialWavelength)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,radialWavelength,options.filter(sum(val,1)))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            xlim([radialWavelength(end) radialWavelength(1)])
        case "k-pseudo-isotropic"
            v = self.transformToPseudoRadialWavenumberA0(val);
            plot(radialWavelength,zeros(size(radialWavelength)),LineWidth=2,Color=0*[1 1 1]), hold on
            plot(ax,radialWavelength,options.filter(v))
            set(gca,'XDir','reverse')
            set(gca,'XScale','log')
            xlim([radialWavelength(end) radialWavelength(1)])
    end
    if isempty(options.simpleName)
        title(ax,fluxes(iComponent).fancyName + " (" + string(sum(val(:))) + " " + self.z_flux_scale_units + ")" )
    else
        title(ax,fluxes(iComponent).simpleName + " (" + string(sum(val(:))) + " " + self.z_flux_scale_units + ")" )
    end
end

switch options.axes
    case "jk"
        % jkR plot labels
        xlabel(tl,'wavelength (km)')
        ylabel(tl,"vertical mode")
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = "enstrophy flux (" + self.z_flux_scale_units + ")";
    case "j"
        % j plot options
        ylabel(tl,"enstrophy flux (" + self.z_flux_scale_units + ")")
        xlabel(tl, 'vertical mode')
    case "jWavenumber"
        % kR plot options
        ylabel(tl,"enstrophy flux (" + self.z_flux_scale_units + ")")
        xlabel(tl,'radius of deformation (km)')
    case "k"
        % kR plot options
        ylabel(tl,"enstrophy flux (" + self.z_flux_scale_units + ")")
        xlabel(tl,'wavelength (km)')
    case "k-pseudo-isotropic"
        % kR plot options
        ylabel(tl,"enstrophy flux (" + self.z_flux_scale_units + ")")
        xlabel(tl,'pseudo-wavelength (km)')
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
    title(tl,"available potential enstrophy flux, day " + minDay + "-" + maxDay)
else
    title(tl,"quasigeostrophic potential enstrophy flux, day " + minDay + "-" + maxDay)
end