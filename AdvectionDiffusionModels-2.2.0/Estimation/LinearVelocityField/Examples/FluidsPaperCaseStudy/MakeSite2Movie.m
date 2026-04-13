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

moviePath = fullfile(movieDir, "Site2Movie.mp4");
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

[u_meso, v_meso] = DecomposeTrajectories(x, y, t, parameterEstimates);
x_meso = x(1, :) + cumtrapz(t, u_meso);
y_meso = y(1, :) + cumtrapz(t, v_meso);

x_camera = mean(x_meso, 2);
y_camera = mean(y_meso, 2);
dx = x - x_camera;
dy = y - y_camera;
cameraXLimits = 1.1 * [min(dx(:)) max(dx(:))];
cameraYLimits = 1.1 * [min(dy(:)) max(dy(:))];

dt = t(2)-t(1);
u_camera = gradient(x_camera, dt);
v_camera = gradient(y_camera, dt);
speed_camera = hypot(u_camera, v_camera);

figureHandle = figure('Units', 'points', 'Position', figurePosition, 'Color', 'w');

for iTime = 1:length(t)
    clf(figureHandle);

    frameIndices = max(1, iTime-tailLength):iTime;

    sp1 = subplot(1, 2, 1);
    plot(sp1, x/distanceScale, y/distanceScale);
    hold(sp1, 'on');
    rectangle('Parent', sp1, 'Position', [x_camera(iTime)+cameraXLimits(1) y_camera(iTime)+cameraYLimits(1) diff(cameraXLimits) diff(cameraYLimits)]/distanceScale, 'LineWidth', 2);
    scatter(sp1, x(iTime, :)/distanceScale, y(iTime, :)/distanceScale, 3^2, 'k', 'filled');
    axis(sp1, 'equal');
    sp1.XGrid = 'on';
    sp1.YGrid = 'on';
    sp1.XMinorGrid = 'on';
    sp1.YMinorGrid = 'on';
    sp1.XTick = -25:20:125;
    sp1.YTick = 0:20:230;
    xlim(sp1, [-30 130]);
    ylim(sp1, [-5 230]);

    sp2 = subplot(1, 2, 2);
    hold(sp2, 'on');
    if iTime > 1
        plot(sp2, x(frameIndices, :)/distanceScale, y(frameIndices, :)/distanceScale, 'LineWidth', 2);
    end
    scatter(sp2, x(iTime, :)/distanceScale, y(iTime, :)/distanceScale, 8^2, 'k', 'filled');
    axis(sp2, 'equal');
    xlim(sp2, (x_camera(iTime)+cameraXLimits)/distanceScale);
    ylim(sp2, (y_camera(iTime)+cameraYLimits)/distanceScale);
    sp2.Box = 'on';
    sp2.XGrid = 'on';
    sp2.YGrid = 'on';
    sp2.XMinorGrid = 'on';
    sp2.YMinorGrid = 'on';
    sp2.GridAlpha = 0.5;
    sp2.MinorGridAlpha = 0.5;
    sp2.XTick = -25:5:125;
    sp2.YTick = -25:5:225;
    sp2.XTickLabel = [];
    sp2.YTickLabel = [];
    title(sp2, sprintf('LatMix 2011 Site 2 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600, 24)), floor(mod(t(iTime)/60, 60))), 'FontSize', 21, 'FontName', 'Helvetica');

    sp1.Position = [0.07 0.04 0.38 0.93];
    p1 = sp1.Position;
    sp2.Position = [0.5 0.03 0.38 0.9];

    banner = annotation('textarrow', 0.96*[1 1], [1 1], 'String', {' Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories.', 'Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 90);
    banner.FontSize = 14;
    banner.FontName = 'Helvetica';

    speedLabel = annotation('textarrow', 0.85*[1 1], 0.9*[1 1], 'String', sprintf('%.1f cm/s', speed_camera(iTime)*100), 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    speedLabel.FontSize = 18;
    speedLabel.FontName = 'Helvetica';

    drawnow;
    writeVideo(videoWriter, getframe(figureHandle));
end

close(videoWriter);
close(figureHandle);
