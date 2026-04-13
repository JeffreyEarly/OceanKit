scriptDir = fileparts(mfilename("fullpath"));
siteNumber = 2;
siteLabel = "Site 2";
dataFilename = "smoothedGriddedRho2Drifters.mat";
nBootstraps = 100;
randomSeed = 0;
scoreStride = 6;
psiS = [3 3 4];
fastS = 3;
mesoscaleConstraint = "zeroVorticity";
shouldSaveFigures = 0;
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", dataFilename);

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);
queryTimes = t;

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:size(x, 2)
    trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3); %#ok<SAGROW>
end

cacheDirectory = fullfile(scriptDir, "BootstrapData");
if exist(cacheDirectory, "dir") == 0
    mkdir(cacheDirectory);
end

psiSTag = join(string(psiS), "-");
cacheFilename = "Rho" + siteNumber + ...
    "GriddedStreamfunctionBootstrapFits" + nBootstraps + ...
    "_seed" + randomSeed + ...
    "_stride" + scoreStride + ...
    "_fastS" + fastS + ...
    "_psiS" + psiSTag + ...
    "_mesoscaleConstraint-" + mesoscaleConstraint + ".nc";
cachePath = fullfile(cacheDirectory, cacheFilename);

if isfile(cachePath)
    bootstrap = GriddedStreamfunctionBootstrap.fromFile(cachePath);
else
    bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
        trajectories, ...
        nBootstraps=nBootstraps, ...
        randomSeed=randomSeed, ...
        queryTimes=queryTimes, ...
        scoreStride=scoreStride, ...
        psiS=psiS, ...
        fastS=fastS, ...
        mesoscaleConstraint=mesoscaleConstraint);
    bootstrap.writeToFile(cachePath, shouldOverwriteExisting=true);
end

trajectories = bootstrap.observedTrajectories;
nDrifters = numel(trajectories);
bestFit = bootstrap.bestFit();
decomposition = bestFit.decomposeTrajectories(trajectories);

scaleFactor = 1;
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

tBootstrap = bootstrap.queryTimes;
tBootstrapDays = tBootstrap/86400;
[~, rankedIndices] = sort(bootstrap.scores.joint, "descend");
indices90 = rankedIndices(1:max(1, round(0.90 * bootstrap.nBootstraps)));
indices68 = rankedIndices(1:max(1, round(0.68 * bootstrap.nBootstraps)));
indexBest = bootstrap.bestBootstrapIndex();

sigma = hypot(bootstrap.summary.sigma_n, bootstrap.summary.sigma_s);
thetaDegrees = GriddedStreamfunction.visualPrincipalStrainAngle(bootstrap.summary.sigma_n, bootstrap.summary.sigma_s);
zeta = bootstrap.summary.zeta;
uCenter = bootstrap.summary.uCenter;
vCenter = bootstrap.summary.vCenter;

errorColor = 0.0 * [1 1 1];
errorAlpha = 0.2;
stdEdgeColor = "none";
meanLineWidth = 2.0;

figure(Color="w", Units="points", Position=[50 50 figureWidthMedium 470]);
tl = tiledlayout(4, 1, TileSpacing="compact", Padding="compact");

ax1 = nexttile;
fillRankedEnvelope(ax1, tBootstrapDays, zeta(:, indices90)/f0, errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax1, tBootstrapDays, zeta(:, indices68)/f0, errorColor, errorAlpha, stdEdgeColor);
hold(ax1, "on");
plot(ax1, tBootstrapDays, zeta(:, indexBest)/f0, LineWidth=meanLineWidth);
ylabel(ax1, "\zeta (f_0)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax1, [tBootstrapDays(1), tBootstrapDays(end)]);
ylim(ax1, [-1 1]);
set(ax1, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax1, "on");
title(ax1, "Site 2 gridded streamfunction bootstrap fit", FontSize=figureTitleSize, FontName=figureFont);

ax2 = nexttile;
fillRankedEnvelope(ax2, tBootstrapDays, sigma(:, indices90)/f0, errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax2, tBootstrapDays, sigma(:, indices68)/f0, errorColor, errorAlpha, stdEdgeColor);
hold(ax2, "on");
plot(ax2, tBootstrapDays, sigma(:, indexBest)/f0, LineWidth=meanLineWidth);
ylabel(ax2, "\sigma (f_0)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax2, [tBootstrapDays(1), tBootstrapDays(end)]);
ylim(ax2, [0 0.8]);
set(ax2, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax2, "on");

ax3 = nexttile;
fillRankedEnvelope(ax3, tBootstrapDays, thetaDegrees(:, indices90), errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax3, tBootstrapDays, thetaDegrees(:, indices68), errorColor, errorAlpha, stdEdgeColor);
hold(ax3, "on");
plot(ax3, tBootstrapDays, thetaDegrees(:, indexBest), LineWidth=meanLineWidth);
ylabel(ax3, "\theta (deg)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax3, [tBootstrapDays(1), tBootstrapDays(end)]);
ylim(ax3, [-140 40]);
yticks(ax3, [-90 -45 0 45 90]);
set(ax3, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax3, "on");

ax4 = nexttile;
fillRankedEnvelope(ax4, tBootstrapDays, uCenter(:, indices90), errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax4, tBootstrapDays, uCenter(:, indices68), errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax4, tBootstrapDays, vCenter(:, indices90), errorColor, errorAlpha, stdEdgeColor);
fillRankedEnvelope(ax4, tBootstrapDays, vCenter(:, indices68), errorColor, errorAlpha, stdEdgeColor);
hold(ax4, "on");
uLine = plot(ax4, tBootstrapDays, uCenter(:, indexBest), LineWidth=meanLineWidth);
vLine = plot(ax4, tBootstrapDays, vCenter(:, indexBest), LineWidth=meanLineWidth);
xlabel(ax4, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(ax4, "u_c/v_c (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax4, [tBootstrapDays(1), tBootstrapDays(end)]);
ylim(ax4, [-5 5]);
legend(ax4, [uLine, vLine], "u_c", "v_c", Location="southwest", FontSize=figureLegendSize, NumColumns=2);
set(ax4, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax4, "on");

title(tl, siteLabel + " bootstrap parameter time series", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site2BootstrapParameters.eps"), "-depsc2");
end

backgroundTrajectory = decomposition.fixedFrame.background(1);
backgroundX = backgroundTrajectory.x(t);
backgroundY = backgroundTrajectory.y(t);
observedSpeedSquaredSum = zeros(size(t));
mesoscaleSpeedSquaredSum = zeros(size(t));
submesoscaleSpeedSquaredSum = zeros(size(t));
decompositionXLimits = [min(backgroundX), max(backgroundX)];
decompositionYLimits = [min(backgroundY), max(backgroundY)];

figure(Color="w", Units="points", Position=[50 50 figureWidthOneColumn 300]);
axMeso = axes;
hold(axMeso, "on");
for iDrifter = 1:nDrifters
    trajectory = trajectories(iDrifter);
    ti = trajectory.t;
    mesoscale = decomposition.fixedFrame.mesoscale(iDrifter);
    plot(axMeso, mesoscale.x(ti)/1000, mesoscale.y(ti)/1000, LineWidth=1.2);
end
axis(axMeso, "equal");
xlabel(axMeso, "x (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axMeso, "y (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
title(axMeso, "Site 2 fixed-frame mesoscale trajectories", FontSize=figureTitleSize, FontName=figureFont);
set(axMeso, FontSize=figureAxisTickSize, FontName=figureFont);
box(axMeso, "on");

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site2BootstrapDecompFigA.eps"), "-depsc2");
end

figure(Color="w", Units="points", Position=[50 50 figureWidthOneColumn 420]);
tlDecomp = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axBackground = nexttile;
plot(axBackground, backgroundX/1000, backgroundY/1000, "k", LineWidth=1.5);
axis(axBackground, "equal");
ylabel(axBackground, "y (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
title(axBackground, "Common background path", FontSize=figureTitleSize, FontName=figureFont);
set(axBackground, FontSize=figureAxisTickSize, FontName=figureFont);
box(axBackground, "on");

axSubmeso = nexttile;
hold(axSubmeso, "on");
for iDrifter = 1:nDrifters
    trajectory = trajectories(iDrifter);
    ti = trajectory.t;
    submesoscale = decomposition.fixedFrame.submesoscale(iDrifter);
    xSubmesoscale = submesoscale.x(ti);
    ySubmesoscale = submesoscale.y(ti);
    [decompositionXLimits, decompositionYLimits] = expandTrajectoryLimits( ...
        decompositionXLimits, decompositionYLimits, xSubmesoscale, ySubmesoscale);
    plot(axSubmeso, xSubmesoscale/1000, ySubmesoscale/1000, LineWidth=1.2);
end
axis(axSubmeso, "equal");
xlabel(axSubmeso, "x (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axSubmeso, "y (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
title(axSubmeso, "Fixed-frame submesoscale trajectories", FontSize=figureTitleSize, FontName=figureFont);
set(axSubmeso, FontSize=figureAxisTickSize, FontName=figureFont);
box(axSubmeso, "on");

xlim(axBackground, paddedLimits(decompositionXLimits/1000));
xlim(axSubmeso, paddedLimits(decompositionXLimits/1000));
ylim(axBackground, paddedLimits(decompositionYLimits/1000));
ylim(axSubmeso, paddedLimits(decompositionYLimits/1000));

title(tlDecomp, siteLabel + " fixed-frame decomposition", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site2BootstrapDecompFigB.eps"), "-depsc2");
end

figure(Color="w", Units="points", Position=[50 50 figureWidthOneColumn 420]);
tlCom = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axCom = nexttile;
hold(axCom, "on");
for iDrifter = 1:nDrifters
    trajectory = trajectories(iDrifter);
    ti = trajectory.t;
    xi = trajectory.x(ti);
    yi = trajectory.y(ti);
    [~, qObserved, rObserved] = bestFit.centeredCoordinates(ti, xi, yi);
    centeredMesoscale = decomposition.centeredFrame.mesoscale(iDrifter);
    fixedMesoscale = decomposition.fixedFrame.mesoscale(iDrifter);
    fixedSubmesoscale = decomposition.fixedFrame.submesoscale(iDrifter);

    if iDrifter == 1
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8 * [1 1 1], DisplayName="Observed");
        plot(axCom, centeredMesoscale.x(ti)/1000, centeredMesoscale.y(ti)/1000, LineWidth=1.2, DisplayName="Mesoscale fit");
    else
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8 * [1 1 1]);
        plot(axCom, centeredMesoscale.x(ti)/1000, centeredMesoscale.y(ti)/1000, LineWidth=1.2);
    end

    observedSpeedSquaredSum = observedSpeedSquaredSum + trajectory.u(ti).^2 + trajectory.v(ti).^2;
    mesoscaleSpeedSquaredSum = mesoscaleSpeedSquaredSum + fixedMesoscale.u(ti).^2 + fixedMesoscale.v(ti).^2;
    submesoscaleSpeedSquaredSum = submesoscaleSpeedSquaredSum + fixedSubmesoscale.u(ti).^2 + fixedSubmesoscale.v(ti).^2;
end
axis(axCom, "equal");
xlabel(axCom, "q (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axCom, "r (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
title(axCom, "Centered trajectories", FontSize=figureTitleSize, FontName=figureFont);
legend(axCom, "Observed", "Mesoscale fit", Location="best", FontSize=figureLegendSize);
set(axCom, FontSize=figureAxisTickSize, FontName=figureFont);
box(axCom, "on");

observedSpeedRms = sqrt(observedSpeedSquaredSum/nDrifters);
mesoscaleSpeedRms = sqrt(mesoscaleSpeedSquaredSum/nDrifters);
backgroundSpeedRms = hypot(backgroundTrajectory.u(t), backgroundTrajectory.v(t));
submesoscaleSpeedRms = sqrt(submesoscaleSpeedSquaredSum/nDrifters);

axRms = nexttile;
plot(axRms, t/86400, observedSpeedRms, Color=0.2 * [1 1 1], LineWidth=1.5);
hold(axRms, "on");
plot(axRms, t/86400, mesoscaleSpeedRms, LineWidth=1.5);
plot(axRms, t/86400, backgroundSpeedRms, LineWidth=1.5);
plot(axRms, t/86400, submesoscaleSpeedRms, LineWidth=1.5);
xlabel(axRms, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axRms, "RMS speed (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axRms, [t(1), t(end)]/86400);
legend(axRms, "Observed", "Mesoscale", "Background", "Submesoscale", Location="best", FontSize=figureLegendSize);
set(axRms, FontSize=figureAxisTickSize, FontName=figureFont);
box(axRms, "on");

title(tlCom, siteLabel + " COM-frame trajectories and RMS speeds", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site2BootstrapDecompFigC.eps"), "-depsc2");
end

function fillRankedEnvelope(ax, tDays, values, faceColor, faceAlpha, edgeColor)
envelope = rankedEnvelope(values);
hold(ax, "on");
fill(ax, [tDays; flipud(tDays)], envelope, faceColor, ...
    EdgeColor=edgeColor, FaceAlpha=faceAlpha);
end

function envelope = rankedEnvelope(values)
envelope = [max(values, [], 2); flipud(min(values, [], 2))];
end

function [xLimits, yLimits] = expandTrajectoryLimits(xLimits, yLimits, xValues, yValues)
xLimits = [min(xLimits(1), min(xValues)), max(xLimits(2), max(xValues))];
yLimits = [min(yLimits(1), min(yValues)), max(yLimits(2), max(yValues))];
end

function limits = paddedLimits(limits)
range = diff(limits);
if range <= 0
    range = max(abs(limits(1)), 1);
end
limits = limits + 0.05 * range * [-1 1];
end
