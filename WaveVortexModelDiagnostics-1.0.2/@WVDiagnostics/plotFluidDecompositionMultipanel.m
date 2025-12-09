function fig = plotFluidDecompositionMultipanel(self,options)
% Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
%
% Plot a multipanel decomposition of the fluid state (wave vs geostrophic).
% Creates a compact multipanel figure showing horizontal (x-y) maps and
% vertical (x-z) sections of vertical vorticity for the total flow, the
% wave component, and the geostrophic component at a specified model time.
% Uses fields from the current WVTransform (self.wvt) and applies consistent
% colormap, plotting scales, and annotations. The x/y axes are presented in
% kilometers and the z (depth) axis in kilometers.
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotFluidDecompositionMultipanel(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter iTime: Time index to display; when provided sets self.iTime (default: self.iTime)
% - Parameter title: Figure title; default uses the model time in days
% - Parameter yForXZSlice: y coordinate (meters) at which to take the x-z slice; default uses center y
% - Returns fig: Handle to the generated figure
arguments
    self WVDiagnostics
    options.visible = "on"
    options.iTime
    options.title
    options.yForXZSlice
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

if ~isfield(options,"title")
    options.title = sprintf('%d days',round(wvt.t/86400));
end

cmDivRWB = self.cmocean('balance'); % diverging positive-negative

% set limits
zeta_limits = [-0.2 0.2];

% location for x-z section
if ~isfield(options,"yForXZSlice")
    iY = round(wvt.Nx/2);
else
    iY = round(options.yForXZSlice/(wvt.y(2)-wvt.y(1)));
end

nColumns = 3;
if nColumns == 3
    figPos = [50 50 600 300];
    % figPos = [50 50 900 600];
else
    figPos = [50 50 600 615];
end

fig = figure(Units='points',Position=figPos,Visible = options.visible);
set(gcf,'PaperPositionMode','auto')


% tl = tiledlayout(2,nColumns,TileSpacing="tight");
tl = tiledlayout(3,nColumns,'TileSpacing', 'tight', 'Padding', 'compact');

if options.title ~= "none"
    title(tl, options.title, 'Interpreter', 'none')
end

% wave and geosgrophic vorticity v_x - u_y
zeta_z_g = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);  % geostrophic
zeta_z_w = wvt.diffX(wvt.v_w) - wvt.diffY(wvt.u_w);  % wave

    function makeVorticityXZPlot(zeta_z)
        val = squeeze(zeta_z(:,iY,:)/wvt.f);
        pcolor(ax, wvt.x/1e3, wvt.z/1e3, val.'), shading interp,
        % axis square
        % axis equal
        % axis tight
        % daspect([1 .02 1]);
        % pbaspect([1 .75 1]);
        ylim([-wvt.Lz/1e3,0])
        colormap(ax, cmDivRWB);
        set(gca,'Layer','top','TickLength',[0.015 0.015])
        clim(ax, zeta_limits);
    end

    function makeVorticityXYPlot(zeta_z)
        val = circshift(zeta_z(:,:,end)/wvt.f,-iY,2);
        pcolor(ax, wvt.x/1e3, wvt.y/1e3, val.'), shading interp,
        % hold on; plot(wvt.x/1e3,ones(size(wvt.x))*wvt.y(iY)/1e3,Color=0*[1 1 1],LineWidth=2); % add line for x-z section
        % axis square
        axis equal
        % daspect([1 1 1]);
        % pbaspect([1 1 1]);
        xlim([0,wvt.Lx/1e3])
        ylim([0,wvt.Ly/1e3])
        colormap(ax, cmDivRWB);
        set(gca,'Layer','top','TickLength',[0.015 0.015])
        clim(ax, zeta_limits);
    end

% geostrophic vorticity section
ax = nexttile(tl,1,[2 1]);
makeVorticityXYPlot(wvt.zeta_z);
xticklabels([])
ylabel('y-distance (km)')
title(ax, "total")

ax = nexttile(tl,2,[2 1]);
makeVorticityXYPlot(zeta_z_w)
xticklabels([])
yticklabels([])
title(ax, "wave")

ax = nexttile(tl,3,[2 1]);
makeVorticityXYPlot(zeta_z_g);
xticklabels([])
yticklabels([])
title(ax, "geostrophic")
% cb = colorbar("eastoutside");
% cb.Label.String = "vertical vorticity (f)";
% cb.Label.Interpreter = 'latex';

ax = nexttile(tl,7);
makeVorticityXZPlot(wvt.zeta_z);
xlabel('x-distance (km)')
ylabel('depth (km)')

ax = nexttile(tl,8);
makeVorticityXZPlot(zeta_z_w)
xlabel('x-distance (km)')
yticklabels([])

ax = nexttile(tl,9);
makeVorticityXZPlot(zeta_z_g);
xlabel('x-distance (km)')
yticklabels([])
% cb = colorbar("eastoutside");
% cb.Label.String = "vertical vorticity (f)";
% cb.Label.Interpreter = 'latex';

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = "vertical vorticity (f)";
cb.Label.Interpreter = 'latex';

end