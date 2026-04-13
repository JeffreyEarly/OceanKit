scriptDir = fileparts(mfilename("fullpath"));

SiteNumber = 1;
dof = 4;
iModel = 3;
bootstrapCount = 1000;
tailLength = 12;
distanceScale = 1000;
figurePosition = [50 50 1920/2 1080/2];

movieDir = fullfile(scriptDir, "Movies");
if exist(movieDir, "dir") == 0
    mkdir(movieDir);
end

moviePath = fullfile(movieDir, "Site1DecompositionMovie.mp4");
videoWriter = VideoWriter(moviePath, "MPEG-4");
videoWriter.FrameRate = 24;
videoWriter.Quality = 100;
open(videoWriter);

load(fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", sprintf("smoothedGriddedRho%dDrifters.mat", SiteNumber)));
x = x(:, 1:end-1);
y = y(:, 1:end-1);

load(fullfile(scriptDir, "BootstrapData", sprintf("Rho%dDrifterSplineFits%d_dof%d.mat", SiteNumber, bootstrapCount, dof)));
[~, indexBest] = max(bootstraps{iModel}.jointlikelihood);

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
cameraXLimits = 1.1 * [min(dx(:)) max(dx(:))];
cameraYLimits = 1.1 * [min(dy(:)) max(dy(:))];

figureHandle = figure('Units', 'points', 'Position', figurePosition, 'Color', 'w');

for iTime = 1:length(t)
    clf(figureHandle);

    frameIndices = max(1, iTime-tailLength):iTime;

    sp1 = subplot(2, 3, [1 2 4 5]);
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
    sp1.XTick = -5:5:25;
    sp1.YTick = 0:5:35;
    sp1.XTickLabel = [];
    sp1.YTickLabel = [];
    title(sp1, sprintf('LatMix 2011 Site 1 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600, 24)), floor(mod(t(iTime)/60, 60))), 'FontSize', 24, 'FontName', 'Helvetica');

    sp2 = subplot(2, 3, 3);
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
    text(sp2, -1, 1.1, 'background', 'FontSize', 16, 'FontName', 'Helvetica');

    sp3 = subplot(2, 3, 6);
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
    text(sp3, -1, 1.1, 'submesoscale', 'FontSize', 16, 'FontName', 'Helvetica');

    sp2.XLim = [-1.25 1.25];
    sp3.XLim = [-1.25 1.25];
    sp2.YLim = [-1.25 1.25];
    sp3.YLim = [-1.25 1.25];

    b = 0.05;
    sp1.Position = [b/2 b 0.66 1-2*b];
    sp2.Position = [0.62 0.51 0.41 0.41];
    sp3.Position = [0.62 0.07 0.41 0.41];

    banner = annotation('textarrow', 0.814*[1 1], 0.045*[1 1], 'String', {' Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories. Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    banner.FontSize = 16;
    banner.FontName = 'Times';

    mesoLabel = annotation('textarrow', 0.13*[1 1], 0.9*[1 1], 'String', 'mesoscale', 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    mesoLabel.FontSize = 16;
    mesoLabel.FontName = 'Helvetica';

    drawnow;
    writeVideo(videoWriter, getframe(figureHandle));
end

close(videoWriter);
close(figureHandle);
