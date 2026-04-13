scriptDir = fileparts(mfilename("fullpath"));

SiteNumber = 2;
dof = 4;
iModel = 5;
bootstrapCount = 1000;
tailLength = 12;
distanceScale = 1000;
figurePosition = [50 50 1920/2 1080/2];

movieDir = fullfile(scriptDir, "Movies");
if exist(movieDir, "dir") == 0
    mkdir(movieDir);
end

moviePath = fullfile(movieDir, "Site2DecompositionMovie.mp4");
videoWriter = VideoWriter(moviePath, "MPEG-4");
videoWriter.FrameRate = 24;
videoWriter.Quality = 100;
open(videoWriter);

load(fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", sprintf("smoothedGriddedRho%dDrifters.mat", SiteNumber)));
x = x(:, 1:end-1);
y = y(:, 1:end-1);

load(fullfile(scriptDir, "BootstrapData", sprintf("Rho%dDrifterSplineFits%d_dof%d.mat", SiteNumber, bootstrapCount, dof)));
[~, mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood, 'descend');
indexBest = mostLikelyIndices(1);

p = bootstraps{iModel};
parameterEstimates = struct( ...
    'u0', p.u0(:, indexBest), ...
    'v0', p.v0(:, indexBest), ...
    'u1', p.u1(:, indexBest), ...
    'v1', p.v1(:, indexBest), ...
    'sigma_n', p.sigma_n(:, indexBest), ...
    'sigma_s', p.sigma_s(:, indexBest), ...
    'zeta', p.zeta(:, indexBest), ...
    'delta', p.delta(:, indexBest));

[u_meso, v_meso, u_bg, v_bg, u_sm, v_sm] = DecomposeTrajectories(x, y, t, parameterEstimates);
x_meso = x(1, :) + cumtrapz(t, u_meso);
y_meso = y(1, :) + cumtrapz(t, v_meso);
x_bg = cumtrapz(t, u_bg);
y_bg = cumtrapz(t, v_bg);
x_sm = cumtrapz(t, u_sm);
y_sm = cumtrapz(t, v_sm);

x_camera = mean(x_meso, 2);
y_camera = mean(y_meso, 2);
dx = x - x_camera;
dy = y - y_camera;
cameraXLimits = [min(dx(:)) max(dx(:))];
cameraYLimits = [min(dy(:)) max(dy(:))];

f0 = 2 * 7.2921E-5 * sin(lat0*pi/180);
scaleFactor = 1.5;
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

errorColor = 0.0*[1 1 1];
errorAlpha = 0.2;
stdEdgeColor = 'none';
meanLineWidth = 2.0;

sigma = sqrt(bootstraps{iModel}.sigma_n.^2 + bootstraps{iModel}.sigma_s.^2);
theta = atan2(bootstraps{iModel}.sigma_s, bootstraps{iModel}.sigma_n)/2;
zeta = bootstraps{iModel}.zeta;
nBootstraps = size(theta, 2);
indices90 = mostLikelyIndices(1:round(0.9*nBootstraps));
indices68 = mostLikelyIndices(1:round(0.68*nBootstraps));
error90 = @(value) [max(value(:, indices90), [], 2); flip(min(value(:, indices90), [], 2))];
error68 = @(value) [max(value(:, indices68), [], 2); flip(min(value(:, indices68), [], 2))];

for iBootstrap = 1:nBootstraps
    flipIndex = find(abs(diff(theta(:, iBootstrap)*180/pi)) > 90, 1, 'first');
    if ~isempty(flipIndex)
        theta(1:flipIndex, iBootstrap) = theta(1:flipIndex, iBootstrap)-pi;
    end
end

figureHandle = figure('Units', 'points', 'Position', figurePosition, 'Color', 'w');

for iTime = 1:length(t)
    clf(figureHandle);

    frameIndices = max(1, iTime-tailLength):iTime;

    sp1 = subplot(6, 3, [2 5 8 11 14 17]);
    plot(sp1, x/distanceScale, y/distanceScale, 'LineWidth', 1, 'Color', 0.8*[1 1 1]);
    hold(sp1, 'on');
    if iTime > 1
        plot(sp1, x_meso(frameIndices, :)/distanceScale, y_meso(frameIndices, :)/distanceScale, 'LineWidth', 2);
    end
    scatter(sp1, x_meso(iTime, :)/distanceScale, y_meso(iTime, :)/distanceScale, 8^2, 'k', 'filled');
    axis(sp1, 'equal');
    xlim(sp1, (x_camera(iTime)+cameraXLimits)/distanceScale);
    ylim(sp1, (y_camera(iTime)+cameraYLimits)/distanceScale);
    sp1.XGrid = 'on';
    sp1.YGrid = 'on';
    sp1.XMinorGrid = 'on';
    sp1.YMinorGrid = 'on';
    sp1.GridAlpha = 0.5;
    sp1.MinorGridAlpha = 0.5;
    sp1.XTick = -25:5:125;
    sp1.YTick = 0:5:230;
    sp1.XTickLabel = [];
    sp1.YTickLabel = [];
    title(sp1, sprintf('LatMix 2011 Site 2 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600, 24)), floor(mod(t(iTime)/60, 60))), 'FontSize', 24, 'FontName', 'Helvetica');

    sp2 = subplot(6, 3, [3 6 9]);
    hold(sp2, 'on');
    if iTime > 1
        plot(sp2, x_bg(frameIndices)/distanceScale, y_bg(frameIndices)/distanceScale, 'LineWidth', 2, 'Color', 0.4*[1 1 1]);
    end
    scatter(sp2, x_bg(iTime)/distanceScale, y_bg(iTime)/distanceScale, 8^2, 'k', 'filled');
    axis(sp2, 'equal');
    sp2.XGrid = 'on';
    sp2.YGrid = 'on';
    sp2.XMinorGrid = 'on';
    sp2.YMinorGrid = 'on';
    sp2.GridAlpha = 0.5;
    sp2.MinorGridAlpha = 0.5;
    sp2.XTick = -5:5:5;
    sp2.YTick = -5:5:5;
    sp2.XTickLabel = [];
    sp2.YTickLabel = [];
    sp2.Box = 'on';
    text(sp2, -2.75, 2.5, 'background', 'FontSize', 16, 'FontName', 'Helvetica');

    sp3 = subplot(6, 3, [12 15 18]);
    hold(sp3, 'on');
    if iTime > 1
        plot(sp3, x_sm(frameIndices, :)/distanceScale, y_sm(frameIndices, :)/distanceScale, 'LineWidth', 2);
    end
    scatter(sp3, x_sm(iTime, :)/distanceScale, y_sm(iTime, :)/distanceScale, 8^2, 'k', 'filled');
    axis(sp3, 'equal');
    sp3.XGrid = 'on';
    sp3.YGrid = 'on';
    sp3.XMinorGrid = 'on';
    sp3.YMinorGrid = 'on';
    sp3.GridAlpha = 0.5;
    sp3.MinorGridAlpha = 0.5;
    sp3.XTick = -5:5:5;
    sp3.YTick = -5:5:5;
    sp3.XTickLabel = [];
    sp3.YTickLabel = [];
    sp3.Box = 'on';
    text(sp3, -2.75, 2.5, 'submesoscale', 'FontSize', 16, 'FontName', 'Helvetica');

    sp2.XLim = 3.2*[-1 1];
    sp3.XLim = 3.2*[-1 1];
    sp2.YLim = 3.2*[-1 1];
    sp3.YLim = 3.2*[-1 1];

    sp4 = subplot(6, 3, [1 4]);
    f = fill([t; flip(t)]/86400, error90(sigma)/f0, errorColor, 'EdgeColor', stdEdgeColor); hold on
    f.FaceAlpha = errorAlpha;
    f = fill([t; flip(t)]/86400, error68(sigma)/f0, errorColor, 'EdgeColor', stdEdgeColor);
    f.FaceAlpha = errorAlpha;
    sp4.ColorOrderIndex = 1;
    plot(t/86400, sigma(:, indexBest)/f0, 'LineWidth', meanLineWidth*scaleFactor);
    ylabel('\sigma (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    xlim([min(t) max(t)]/86400);
    sp4.FontSize = figure_axis_tick_size;
    sp4.XTickLabel = [];
    ylimits = [0 0.8];
    ylim(ylimits);
    plot([t(iTime) t(iTime)]/86400, ylimits, 'LineWidth', 2, 'Color', 0*[1 1 1]);

    sp5 = subplot(6, 3, [7 10]);
    f = fill([t; flip(t)]/86400, error90(theta)*180/pi, errorColor, 'EdgeColor', stdEdgeColor); hold on
    f.FaceAlpha = errorAlpha;
    f = fill([t; flip(t)]/86400, error68(theta)*180/pi, errorColor, 'EdgeColor', stdEdgeColor);
    f.FaceAlpha = errorAlpha;
    sp5.ColorOrderIndex = 1;
    plot(t/86400, theta(:, indexBest)*180/pi, 'LineWidth', meanLineWidth*scaleFactor);
    ylabel('\theta (°)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    xlim([min(t) max(t)]/86400);
    ylimits = [-140 40];
    ylim(ylimits);
    plot([t(iTime) t(iTime)]/86400, ylimits, 'LineWidth', 2, 'Color', 0*[1 1 1]);
    yticks([-90 -45 0 45 90]);
    sp5.FontSize = figure_axis_tick_size;
    sp5.XTickLabel = [];

    sp6 = subplot(6, 3, [13 16]);
    f = fill([t; flip(t)]/86400, error90(zeta)/f0, errorColor, 'EdgeColor', stdEdgeColor); hold on
    f.FaceAlpha = errorAlpha;
    f = fill([t; flip(t)]/86400, error68(zeta)/f0, errorColor, 'EdgeColor', stdEdgeColor);
    f.FaceAlpha = errorAlpha;
    sp6.ColorOrderIndex = 7;
    plot(t/86400, zeta(:, indexBest)/f0, 'LineWidth', meanLineWidth*scaleFactor);
    ylabel('\zeta (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    xlim([min(t) max(t)]/86400);
    ylimits = [-1 1];
    ylim(ylimits);
    plot([t(iTime) t(iTime)]/86400, ylimits, 'LineWidth', 2, 'Color', 0*[1 1 1]);
    sp6.FontSize = figure_axis_tick_size;
    xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

    b = 0.07;
    dxPosition = 0.06;
    dyPosition = -0.005;
    width = 0.3;
    height = 0.27;
    sp4.Position = [dxPosition b+2*dyPosition+2*0.3 width height];
    sp5.Position = [dxPosition b+dyPosition+0.3 width height];
    sp6.Position = [dxPosition b width height];

    sp1.Position = [0.3 b 0.5 1-2*b];
    width = 0.42;
    height = width;
    sp2.Position = [0.64 0.51 width height];
    sp3.Position = [0.64 0.07 width height];

    banner = annotation('textarrow', 0.97*[1 1], 0.04*[1 1], 'String', {'Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories.', 'Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    banner.FontSize = 16;
    banner.FontName = 'Times';

    mesoLabel = annotation('textarrow', 0.48*[1 1], 0.92*[1 1], 'String', 'mesoscale', 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    mesoLabel.FontSize = 16;
    mesoLabel.FontName = 'Helvetica';

    drawnow;
    writeVideo(videoWriter, getframe(figureHandle));
end

close(videoWriter);
close(figureHandle);
