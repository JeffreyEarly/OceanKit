function fig = plotFluidStateMultipanel(self,options)
% Plot multipanel summary of fluid state and spectra
%
% Create a compact multipanel figure showing horizontal (x-y) maps and
% vertical (x-z) sections of vertical vorticity for the total flow, the
% wave component, and the geostrophic component at a specified model time.
% Optionally includes log-energy spectra with KE/PE and frequency/wavelength
% contours. Axes are annotated in kilometers and depth in kilometers; color
% limits and annotation ticks are handled internally.
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotFluidStateMultipanel(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter visible: (optional) input argument `visible` (default: "on")
% - Parameter iTime: time-related parameter `iTime`
% - Parameter title: input argument `title`
% - Parameter shouldShowEnergySpectra: (optional) input argument `shouldShowEnergySpectra` (default: true)
% - Parameter shouldShowTotalFields: (optional) input argument `shouldShowTotalFields` (default: false)
% - Parameter figureHandle: input argument `figureHandle`
% - Parameter wavelengths: (optional) input argument `wavelengths` (default: [1,2,5,10,20,50,100,200,500])
% - Parameter wavelengthColor: (optional) input argument `wavelengthColor` (default: [.5,.5,.5])
% - Parameter frequencies: (optional) input argument `frequencies` (default: [1.01 1.05 1.2 2 4 8 16])
% - Parameter frequencyColor: (optional) input argument `frequencyColor` (default: [.7,.7,.7])
% - Parameter keFractions: (optional) input argument `keFractions` (default: [.01,.1,.25,.5,.75,.9,.99])
% - Parameter keFractionColor: (optional) input argument `keFractionColor` (default: [.7,.7,.7])
% - Parameter labelSpacing: (optional) input argument `labelSpacing` (default: 1000)
% - Parameter lineWidth: (optional) input argument `lineWidth` (default: 1)
% - Returns fig: Figure handle for the generated plot
arguments
    self WVDiagnostics
    options.visible = "on"
    options.iTime
    options.title
    options.shouldShowEnergySpectra = true
    options.shouldShowTotalFields = false
    options.figureHandle
    options.wavelengths = [1,2,5,10,20,50,100,200,500];
    options.wavelengthColor = [.5,.5,.5];
    options.frequencies = [1.01 1.05 1.2 2 4 8 16];
    options.frequencyColor = [.7,.7,.7];
    options.keFractions = [.01,.1,.25,.5,.75,.9,.99];
    options.keFractionColor = [.7,.7,.7];
    options.labelSpacing = 1000;
    options.lineWidth = 1;
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

if ~isfield(options,"title")
    options.title = sprintf('%d days',round(wvt.t/86400));
end

cmDivRWB = self.cmocean('balance'); % diverging positive-negative
cmLinMono = self.cmocean('dense') ; % linear monochrome for data>0

% set limits
zeta_limits = [-0.2 0.2];
energy_limits = [-8 0];
enstrophy_limits = [-16 -9];

% Create wavelength axes for wvt variables
% We will pretend the "0" wavenumber is actually evenly spaced
% from the nearest two wavenumbers
kPseudoLocation = wvt.kRadial;
kPseudoLocation(1) = exp(-log(kPseudoLocation(3)) + 2*log(kPseudoLocation(2)));
jPseudoLocation = self.jWavenumber(wvt.j+1);
jPseudoLocation(1) = exp(-log(jPseudoLocation(3)) + 2*log(jPseudoLocation(2)));
[KPseudoLocation,JPseudoLocation] = ndgrid(kPseudoLocation,jPseudoLocation);
kPseudoRadial = sqrt(JPseudoLocation.^2 + KPseudoLocation.^2);

% For interpolation to work correctly we need to repeat the
% first entry, but properly back at zero
% NOTE: these are from wvd, including the anti-aliased modes.
kPseudoLocationWVD = self.kRadial;
kPseudoLocationWVD(1) = exp(-log(kPseudoLocationWVD(3)) + 2*log(kPseudoLocationWVD(2)));
jPseudoLocationWVD = self.jWavenumber;
jPseudoLocationWVD(1) = exp(-log(jPseudoLocationWVD(3)) + 2*log(jPseudoLocationWVD(2)));
kPaddedWVD = cat(1,0,kPseudoLocationWVD);
jPaddedWVD = cat(1,0,jPseudoLocationWVD);
[KPaddedWVD,JPaddedWVD] = ndgrid(kPaddedWVD,jPaddedWVD);


% create some nice tick labels to show wavelength
% ticks_x = [500;300;200;100;50;30;20;10];
% ticks_x = [500;50;30;20;15;12;10];
% % % ticks_x = (2*pi./(1000.*linspace(wvt.kRadial(2),wvt.kRadial(end),6)));
% % % labels_x = cell(length(ticks_x),1);
% % % for i=1:length(ticks_x)
% % %     labels_x{i} = sprintf('%.0f',ticks_x(i));
% % % end
% % % ticks_x = 2*pi./(1e3*ticks_x);

% color for some annotations
textColor = '[.5,.5,.5]';

% location for x-z section
% iY = round(wvt.Nx/2);
iY = 1;

% % % % % create the lines of constant frequency
% % % % [omegaN,n] = wvt.transformToRadialWavenumber(abs(wvt.Omega),ones(size(wvt.Omega)));
% % % % omegaJK = (omegaN./n)/wvt.f;

% create the lines of constant deformation radius
deformationJK = repmat(sqrt(wvt.Lr2)./1000,1,length(wvt.kRadial));

nColumns = 2;
if options.shouldShowEnergySpectra
    nColumns = nColumns + 1;
end
if nColumns == 3
    figPos = [50 50 900 600];
else
    figPos = [50 50 600 615];
end

if ~isfield(options,"figureHandle")
    fig = figure(Units='points',Position=figPos,Visible = options.visible);
    set(gcf,'PaperPositionMode','auto')
else
    fig = options.figureHandle;
    clf(options.figureHandle)
    set(0, 'currentfigure', options.figureHandle);
end


% % % tl = tiledlayout(2,nColumns,TileSpacing="tight");
tl = tiledlayout(1,nColumns,TileSpacing="tight");

title(tl, options.title, 'Interpreter', 'none')

% compute some quantities
TE_A0_j_kl = wvt.A0_TE_factor .* abs(wvt.A0).^2; % m^2/s^3
TE_A0_j_kR = wvt.transformToRadialWavenumber(TE_A0_j_kl);
TE_Apm_j_kl = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2); % m^2/s^3
TE_Apm_j_kR = wvt.transformToRadialWavenumber(TE_Apm_j_kl);
TZ_A0_j_kl = wvt.A0_TZ_factor .* (wvt.A0.*conj(wvt.A0));
TZ_A0_j_kR = wvt.transformToRadialWavenumber(TZ_A0_j_kl);

% wave and geosgrophic vorticity v_x - u_y
zeta_z_g = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);  % geostrophic
zeta_z_w = wvt.diffX(wvt.v_w) - wvt.diffY(wvt.u_w);  % wave

% nested tiled layout allows common colorbar for subset of axes.
tl_inner = tiledlayout(tl,2,2,TileSpacing='tight');
tl_inner.Layout.TileSpan = [1,2];

% geostrophic surface vorticity
% % % ax = nexttile(tl,1);
ax = nexttile(tl_inner,1);
val = zeta_z_g(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:');hold off; % add line for x-z section
title("geostrophic vorticity \zeta_g")
% title("geostrophic surface vorticity")
axis square
xticklabels([])
% yticklabels([])
ylabel('y-distance (km)')
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
colormap(ax, cmDivRWB);
clim(ax, zeta_limits);

% geostrophic vorticity section
% % % ax = nexttile(tl,nColumns+1);
ax = nexttile(tl_inner,3);
val = squeeze(zeta_z_g(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
% title("geostrophic x-z vorticity")
% title("x-z vorticity")
xlabel('x-distance (km)')
ylabel('Depth (km)')
axis square
colormap(ax, cmDivRWB);
ylabel('Depth (km)')
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);
% % % cb = colorbar("southoutside");
% % % cb.Label.String = "$\zeta$/f";
% % % cb.Label.Interpreter = 'latex';

% wave surface vorticity
% % % ax = nexttile(tl,2);
ax = nexttile(tl_inner,2);
val = zeta_z_w(:,:,end)/wvt.f;
pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,'k:');hold off; % add line for x-z section
% title("surface vorticity")
title("wave vorticity \zeta_w")
axis square
colormap(ax, cmDivRWB);
yticklabels([])
xticklabels([])
% xlabel('x-distance (km)')
% ylabel('y-distance (km)')
set(gca,'YTick',xticks,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);

% wave vorticity section
% % % ax = nexttile(tl,nColumns+2);
ax = nexttile(tl_inner,4);
val = squeeze(zeta_z_w(:,iY,:)/wvt.f);
pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
% title("wave x-z vorticity")
% title("x-z vorticity")
axis square
colormap(ax, cmDivRWB);
xlabel('x-distance (km)')
% ylabel('Depth (km)')
yticklabels([])
set(gca,'Layer','top','TickLength',[0.015 0.015])
clim(ax, zeta_limits);
% % % cb = colorbar("southoutside");
% % % cb.Label.String = "$\zeta$/f";
% % % cb.Label.Interpreter = 'latex';

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.String = "$(f)$"; %"$\zeta$/f";
cb.Label.Interpreter = 'latex';

if options.shouldShowEnergySpectra

    % nested tiled layout allows common colorbar for subset of axes.
    tl_inner = tiledlayout(tl,2,1,TileSpacing='tight');
    tl_inner.Layout.Tile = 3;

    % geostrophic energy spectrum
    % % % ax = nexttile(tl,3);
    ax = nexttile(tl_inner,1);
    val = log10(TE_A0_j_kR);
    pcolor(ax,2*pi./kPseudoLocation/1000,2*pi./jPseudoLocation/1000,val), shading flat,
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    set(gca,'YScale','log')
    xticklabels([])
    title('geostrophic energy spectrum')
    % title('geostrophic energy spectrum','Units', 'normalized', 'Position', [0.5, 0.93, 0])
    axis square
    clim(ax,energy_limits);
    ylabel("radius of deformation (km)")
    set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
    colormap(ax, cmLinMono);
    set(gca,'layer','top'),
    hold on
    % add ke:pe ratio contours. 
    % Note: have to pad quantities to work right with our log-log axes. And here
    % have to "double pad" to work right with pcolor's box shift and
    % reversed axis direction. 
    fraction = self.geo_hke_jk./(self.geo_hke_jk+self.geo_pe_jk);
    % fractionPadded = cat(1,fraction(1,:),fraction);
    fractionPadded = cat(1,[fraction(1,:);fraction(1,:)],fraction(1:end-1,:));
    % fractionPadded = cat(2,fractionPadded(:,1),fractionPadded);
    fractionPadded = cat(2,[fractionPadded(:,1) fractionPadded(:,1)],fractionPadded(:,1:end-1));
    fractionJK = interpn(KPaddedWVD,JPaddedWVD,fractionPadded.',KPseudoLocation,JPseudoLocation,"linear");
    [C,h] = contour(ax,2*pi./KPseudoLocation/1000,2*pi./JPseudoLocation/1000,fractionJK,options.keFractions,'LineWidth',options.lineWidth,'Color',options.keFractionColor, DisplayName="KE/(KE+PE)", HandleVisibility='off');
    clabel(C,h,options.keFractions,'Color',options.keFractionColor,'LabelSpacing',options.labelSpacing)
    % add pseudoWavelength 
    [C,h] = contour(ax,2*pi./KPseudoLocation/1000,2*pi./JPseudoLocation/1000,2*pi./kPseudoRadial/1000,options.wavelengths,'LineWidth',options.lineWidth,'Color',options.wavelengthColor, DisplayName="pseudo-wavelength (km)");
    clabel(C,h,options.wavelengths,'Color',options.wavelengthColor,'LabelSpacing',options.labelSpacing)
    hold off

    % wave energy spectrum
    % % % ax = nexttile(tl,nColumns + 3);
    ax = nexttile(tl_inner,2);
    val = log10(TE_Apm_j_kR);
    pcolor(ax,2*pi./kPseudoLocation/1000,2*pi./jPseudoLocation/1000,val), shading flat,
    set(gca,'XDir','reverse')
    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    set(gca,'YScale','log')
    axis square
    cb = colorbar('southoutside');
    clim(ax,energy_limits);
    title('wave energy spectrum')
    % title('wave energy spectrum','Units', 'normalized', 'Position', [0.5, 0.93, 0])
    xlabel('horizontal wavelength (km)')
    % clim(ax,energy_limits);
    ylabel("radius of deformation (km)")
    set(gca, 'YAxisLocation', 'right','Layer','top','TickLength',[0.015 0.015]);
    colormap(ax, cmLinMono);
    cb.Label.String = "log10(m^3 s^{-2})";
    set(gca,'layer','top'),
    hold on
    % add frequency contours
    % omegaPadded = cat(1,self.omega_jk(1,:),self.omega_jk);
    omegaPadded = cat(1,[self.omega_jk(1,:);self.omega_jk(1,:)],self.omega_jk(1:end-1,:));
    % omegaPadded = cat(2,omegaPadded(:,1),omegaPadded);
    omegaPadded = cat(2,[omegaPadded(:,1) omegaPadded(:,1)],omegaPadded(:,1:end-1));
    omegaJK = interpn(KPaddedWVD,JPaddedWVD,omegaPadded.',KPseudoLocation,JPseudoLocation,"linear");
    [C,h] = contour(ax,2*pi./KPseudoLocation/1000,2*pi./JPseudoLocation/1000,omegaJK/self.wvt.f,options.frequencies,'LineWidth',options.lineWidth,'Color',options.frequencyColor, DisplayName="frequency (f)", HandleVisibility='off');
    clabel(C,h,options.frequencies,'Color',options.frequencyColor,'LabelSpacing',options.labelSpacing)
    % add pseudoWavelength 
    [C,h] = contour(ax,2*pi./KPseudoLocation/1000,2*pi./JPseudoLocation/1000,2*pi./kPseudoRadial/1000,options.wavelengths,'LineWidth',options.lineWidth,'Color',options.wavelengthColor, DisplayName="pseudo-wavelength (km)");
    clabel(C,h,options.wavelengths,'Color',options.wavelengthColor,'LabelSpacing',options.labelSpacing)
    hold off


end

% Add a large label to the left of the first row (row 1)
% textAxes = axes('Position', [0, 0, 1, 1], 'Visible', 'off'); % Dummy invisible axes
% text(0.07, 0.74, 'Geostrophic', 'FontSize', 12, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label
% text(0.07, 0.34, 'Wave', 'FontSize', 12, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Rotation', 90); % Rotated vertical label

end