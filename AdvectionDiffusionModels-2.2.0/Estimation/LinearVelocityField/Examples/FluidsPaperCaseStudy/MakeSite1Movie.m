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

moviePath = fullfile(movieDir, "Site1Movie.mp4");
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
    xlim(sp1, [-2.5 20]);
    ylim(sp1, [-2.5 37.5]);

    sp2 = subplot(1, 2, 2);
    hold(sp2, 'on');
    if iTime > 1
        plot(sp2, x(frameIndices, :)/distanceScale, y(frameIndices, :)/distanceScale, 'LineWidth', 2);
    end
    scatter(sp2, x(iTime, :)/distanceScale, y(iTime, :)/distanceScale, 8^2, 'k', 'filled');
    axis(sp2, 'equal');
    xlim(sp2, (x_camera(iTime)+cameraXLimits)/distanceScale);
    ylim(sp2, (y_camera(iTime)+cameraYLimits)/distanceScale);
    sp2.XGrid = 'on';
    sp2.YGrid = 'on';
    sp2.XMinorGrid = 'on';
    sp2.YMinorGrid = 'on';
    sp2.GridAlpha = 0.5;
    sp2.MinorGridAlpha = 0.5;
    sp2.XTick = -5:5:25;
    sp2.YTick = 0:5:35;
    sp2.XTickLabel = [];
    sp2.YTickLabel = [];
    title(sp2, sprintf('LatMix 2011 Site 1 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600, 24)), floor(mod(t(iTime)/60, 60))), 'FontSize', 24, 'FontName', 'Helvetica');

    b = 0.05;
    p1 = sp1.Position;
    sp1.Position = [b/2 b p1(3) 1-2*b];
    p1 = sp1.Position;
    sp2.Position = [p1(1)+p1(3) p1(2) 1-p1(3)-b 1-2*b];

    banner = annotation('textarrow', [p1(1)+p1(3) b], [0.075 b], 'String', {' Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories.', 'Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle', 'none', 'LineStyle', 'none', 'TextRotation', 0);
    banner.FontSize = 16;
    banner.FontName = 'Times';

    drawnow;
    writeVideo(videoWriter, getframe(figureHandle));
end

close(videoWriter);
close(figureHandle);
