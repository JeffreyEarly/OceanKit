function fig = plotEnergySpectrum(self,options)
% Plot the wave/geostrophic energy spectra at a given time.
%
% Plot the wave/geostrophic energy spectra at a given time
% Makes a nice multiplanel plot of the wave and geostrophic spectra at a
% given time.
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotEnergySpectrum(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter iTime: time index in model output file
% - Parameter visible: (optional) figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.iTime
    options.visible = "on"
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% compute energy and enstrophy spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total energy, A0
TE_A0_j_kl = wvt.A0_TE_factor .* abs(wvt.A0).^2; % m^2/s^3
TE_A0_j_kR = wvt.transformToRadialWavenumber(TE_A0_j_kl);
TE_A0_kR = sum(TE_A0_j_kR,1);
TE_A0_j = sum(TE_A0_j_kR,2);

% total energy, Apm
TE_Apm_j_kl = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2); % m^2/s^3
TE_Apm_j_kR = wvt.transformToRadialWavenumber(TE_Apm_j_kl);

% total energy, inertial
maskApInertial = wvt.inertialComponent.maskAp; % m^2/s^3
maskAmInertial = wvt.inertialComponent.maskAm;
TE_inertial_j_kl = wvt.Apm_TE_factor .* (abs(maskApInertial.*wvt.Ap).^2 + abs(maskAmInertial.*wvt.Am).^2);
TE_inertial_j_kR = wvt.transformToRadialWavenumber(TE_inertial_j_kl);
TE_inertial_kR = sum(TE_inertial_j_kR,1);
TE_inertial_j = sum(TE_inertial_j_kR,2);

% total energy, wave
maskApWave = wvt.waveComponent.maskAp; % m^2/s^3
maskAmWave = wvt.waveComponent.maskAm;
TE_wave_j_kl = wvt.Apm_TE_factor .* (abs(maskApWave.*wvt.Ap).^2 + abs(maskAmWave.*wvt.Am).^2);
TE_wave_j_kR = wvt.transformToRadialWavenumber(TE_wave_j_kl);
TE_wave_kR = sum(TE_wave_j_kR,1);
TE_wave_j = sum(TE_wave_j_kR,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Combined energy spectrum figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default line color order, for repeating
linesTemp = lines;

% create radial wavelength vector
radialWavelength = 2*pi./wvt.kRadial/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

% create j vector for log y-axis.
jForLogAxis = wvt.j;
jForLogAxis(1) = 0.75;

% create figure
fig = figure('Units', 'points', 'Position', [50 50 700 500],'Visible',options.visible);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
tl = tiledlayout(2,2,TileSpacing='tight');
title(tl,'Energy Spectrum')

% plot the wave energy
val = log10(TE_Apm_j_kR);
axIGW = nexttile;
pcolor(radialWavelength,jForLogAxis,val), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
set(gca,'YScale','log')
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
colormap(axIGW, self.cmocean('dense'));
ylabel('vertical mode')
title('Internal Gravity Wave')
xlabel('wavelength (km)')
text(radialWavelength(1),max(jForLogAxis),'inertial','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)

self.overlayFrequencyContours(frequencies = [1.01 1.05 1.2 1.5 2 4 8 16],textColor = [.5,.5,.5],labelSpacing = 400,lineWidth = 1)

% plot the geostrophic energy
val = log10(TE_A0_j_kR);
axGEO = nexttile;
pcolor(radialWavelength,jForLogAxis,val), shading flat
set(gca,'XDir','reverse')
set(gca,'XScale','log')
set(gca,'YScale','log')
clim([max(TE_Apm_j_kR(:))-6 max(TE_Apm_j_kR(:))])
colormap(axGEO, self.cmocean('dense'));
set(gca,'YTickLabel',[]);
title('Geostrophic')
xlabel('wavelength (km)')
text(radialWavelength(1),max(jForLogAxis),'MDA','FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left')
line([radialWavelength(2),radialWavelength(2)],[min(jForLogAxis),max(jForLogAxis)],'Color','k','LineWidth',1)

self.overlayGeostrophicKineticPotentialFractionContours
self.showRossbyRadiusYAxis(textColor=[.5,.5,.5])

% plot vertical mode spectrum
axJ = nexttile;
plot(wvt.j,TE_inertial_j+TE_wave_j+TE_A0_j,wvt.j,TE_A0_j,wvt.j,TE_inertial_j+TE_wave_j)
hold on
plot(wvt.j,TE_wave_j,'--','Color',linesTemp(3,:))
plot(wvt.j,TE_inertial_j,':','Color',linesTemp(3,:))
yscale('log')
ylabel('energy (m^3 s^{-2})');
xlabel('vertical mode');
axis tight
title('Vertical Mode Spectrum')
legend('Total','Geostrophic','IO+IGW','IGW','IO','Location','southwest')

% plot horizontal wavenumber spectrum
axK = nexttile;
plot(radialWavelength,TE_inertial_kR+TE_wave_kR+TE_A0_kR,radialWavelength,TE_A0_kR,radialWavelength,TE_inertial_kR+TE_wave_kR)
set(gca,'XDir','reverse')
xscale('log'); yscale('log')
axis tight
title('Radial Wavenumber Spectrum')
legend('Total','Geostrophic','IO+IGW','Location','southwest')
xlabel('wavelength (km)')
% yticklabels([])

% match limits
xlimK = get(axK,'xlim');
ylimK = get(axK,'ylim');
ylimJ = get(axJ,'ylim');
set(axIGW,'xlim',xlimK);
set(axGEO,'xlim',xlimK);
set(axJ,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])
set(axK,'ylim',[min([ylimK,ylimJ]),max([ylimK,ylimJ])])

% colorbar
cb = colorbar(axIGW);
cb.Layout.Tile = 'south';
cb.Label.String = "log10(m^3 s^{-2})";


end