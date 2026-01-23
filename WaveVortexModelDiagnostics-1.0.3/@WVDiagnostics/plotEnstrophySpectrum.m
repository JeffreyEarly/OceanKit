function fig = plotEnstrophySpectrum(self,options)
% Plot the enstrophy spectrum at a given time.
%
% Plot the enstrophy spectrum at a given time
% Produces a multiplanel figure of geostrophic and apv enstrophy diagnostics at a
% specified time index. Shows available potential enstrophy (APV), the
% quadratic (QGPV) approximation, their difference, plus modal and radial
% spectra. Supports a 2x2 "four-panel" layout or a single-panel
% pseudo-radial projection ("kPseudoRadial").
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotEnstrophySpectrum(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter iTime: time index in model output file (default: self.iTime)
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter title: Figure title (default: 'Enstrophy Spectrum')
% - Parameter clim: (optional) Color limits for log10 pcolor panels (default: [-13 -7])
% - Parameter figureHandle: Existing figure handle to draw into (default: create new figure)
% - Parameter style: (optional) Plot style, one of "four-panel" or "kPseudoRadial" (default: "four-panel")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.iTime
    options.visible = "on"
    options.title
    options.clim = [-13 -7]
    options.figureHandle
    options.style {mustBeMember(options.style,["four-panel","kPseudoRadial"])} = "four-panel"
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

if ~isfield(options,"title")
    options.title = 'Enstrophy Spectrum';
end

if self.diagnosticsHasExplicitAntialiasing
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% compute energy and enstrophy spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total enstrophy, A0
% TZ_A0_j_kl = wvt.A0_TZ_factor .* abs(wvt.A0).^2;
% prefactor = wvt.h_0/2; prefactor(1) = wvt.Lz/2;
% qgpv_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.qgpv));
% apv_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.apv));
% error_bar = apv_bar-qgpv_bar;
% TZ_A0_j_kl = prefactor.*abs(qgpv_bar).^2;
% TZ_APV_j_kl = prefactor.*abs(apv_bar).^2;
% TZ_Error_j_kl = prefactor.*abs(error_bar).^2;

qgpv_bar = mean(mean(wvt.qgpv,1),2);

TZ_A0_j_kl = self.spectrumWithFgTransform(wvt.qgpv - qgpv_bar,useExplicitAntialiasedWVT=self.diagnosticsHasExplicitAntialiasing);
TZ_APV_j_kl = self.spectrumWithFgTransform(wvt.apv,useExplicitAntialiasedWVT=self.diagnosticsHasExplicitAntialiasing);
TZ_Error_j_kl = self.spectrumWithFgTransform(wvt.apv - wvt.qgpv+ qgpv_bar,useExplicitAntialiasedWVT=self.diagnosticsHasExplicitAntialiasing);

TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);
TZ_A0_kR = sum(TZ_A0_j_kR,1);
TZ_A0_j = sum(TZ_A0_j_kR,2);
TZ_APV_j_kR = wvt.transformToRadialWavenumber(TZ_APV_j_kl);
TZ_APV_kR = sum(TZ_APV_j_kR,1);
TZ_APV_j = sum(TZ_APV_j_kR,2);
TZ_Error_j_kR = wvt.transformToRadialWavenumber(TZ_Error_j_kl);
TZ_Error_kR = sum(TZ_Error_j_kR,1);
TZ_Error_j = sum(TZ_Error_j_kR,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% enstrophy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~isfield(options,"figureHandle")
    fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Color', 'w');
else
    fig = options.figureHandle;
    clf(options.figureHandle)
    set(0, 'currentfigure', options.figureHandle);
end

if options.style == "kPseudoRadial"
    TZ_A0 = self.transformToPseudoRadialWavenumberA0(TZ_A0_j_kR);
    TZ_APV = self.transformToPseudoRadialWavenumberA0(TZ_APV_j_kR);
    TZ_Error = self.transformToPseudoRadialWavenumberA0(TZ_Error_j_kR);

    % create radial wavelength vector
    radialWavelength = 2*pi./self.kPseudoRadial/1000;
    radialWavelength(1) = 1.5*radialWavelength(2);

    plot(radialWavelength,TZ_APV), hold on
    plot(radialWavelength,TZ_A0),
    plot(radialWavelength,TZ_Error,Color=0*[1 1 1])
    set(gca,'XDir','reverse')
    xscale('log'); yscale('log')
    axis tight
    title('Potential Enstrophy Spectrum')
    ylabel('enstrophy (m s^{-2})');
    xlabel('pseudo-wavelength (km)')
    ylim(10.^options.clim);
        legend('apv','qgpv','error')
else

    % create radial wavelength vector
    radialWavelength = 2*pi./wvt.kRadial/1000;
    radialWavelength(1) = 1.5*radialWavelength(2);

    % create j vector for log y-axis.
    jForLogAxis = wvt.j;
    jForLogAxis(1) = 0.75;


    tl = tiledlayout(2,2,TileSpacing="compact");
    if options.title ~= "none"
        title(tl, options.title, 'Interpreter', 'none')
    end

    % plot the available potential enstrophy
    val = log10((TZ_APV_j_kR).');
    axIGW = nexttile;
    pcolor(radialWavelength,jForLogAxis,val.'), shading flat
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    ylabel('vertical mode')
    title('APV')
    xlabel('wavelength (km)')
    colormap(axIGW, self.cmocean('dense'));
    text(radialWavelength(1),max(jForLogAxis),'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
    line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1.5)
    clim(options.clim)
    self.overlayGeostrophicKineticPotentialFractionContours


    % plot the geostrophic enstrophy
    val = log10((TZ_A0_j_kR).');
    axGEO = nexttile;
    pcolor(radialWavelength,jForLogAxis,val.'), shading flat
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'YTickLabel',[]);
    title('QGPV')
    xlabel('wavelength (km)')
    colormap(axGEO, self.cmocean('dense'));
    text(radialWavelength(1),max(jForLogAxis),'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
    line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1.5)
    clim(options.clim)

    self.overlayGeostrophicKineticPotentialFractionContours
    self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

    % plot vertical mode spectrum
    axJ = nexttile;
    plot(wvt.j,TZ_APV_j), hold on
    plot(wvt.j,TZ_A0_j)
    plot(wvt.j,TZ_Error_j,Color=0*[1 1 1])
    yscale('log')
    ylim(10.^options.clim);
    axis tight
    ylabel('enstrophy (m s^{-2})');
    xlabel('vertical mode');
    title('Vertical Mode Spectrum')
    legend('apv','qgpv','error')

    % plot horizontal wavenumber spectrum
    axK = nexttile;
    plot(radialWavelength,TZ_APV_kR), hold on
    plot(radialWavelength,TZ_A0_kR),
    plot(radialWavelength,TZ_Error_kR,Color=0*[1 1 1])
    set(gca,'XDir','reverse')
    xscale('log'); yscale('log')
    axis tight
    title('Radial Wavenumber Spectrum')
    % ylabel('enstrophy (m s^{-2})');
    xlabel('wavelength (km)')
    ylim(10.^options.clim);

    % match limits
    % xlimK = get(axK,'xlim');
    % ylimK = get(axK,'ylim');
    % ylimJ = get(axJ,'ylim');
    % % set(axIGW,'xlim',xlimK);
    % set(axGEO,'xlim',xlimK);
    set(axJ,'xlim',[min(wvt.j) max(wvt.j)])
    set(axJ,'ylim',10.^options.clim)
    set(axK,'xlim',get(axGEO,'xlim'))
    set(axK,'ylim',10.^options.clim)


    % colorbar
    cb = colorbar(axGEO);
    cb.Layout.Tile = 'south';
    cb.Label.String = "log10(m s^{-2})";
end

end