% Make simple Site 2 figures directly from a gridded streamfunction fit.
scriptDir = fileparts(mfilename("fullpath"));
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", "smoothedGriddedRho2Drifters.mat");
if ~isfile(dataPath)
    error("GriddedStreamfunction:MissingExampleData", "Expected local example data at %s.", dataPath);
end

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
nDrifters = size(x, 2);
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);

trajectoryCell = cell(nDrifters, 1);
for iDrifter = 1:nDrifters
    trajectoryCell{iDrifter} = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
end
trajectories = vertcat(trajectoryCell{:});

fit = GriddedStreamfunction.fromTrajectories(trajectories, psiS=[2 2 4]);

scaleFactor = 1;
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

tDays = t/86400;
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);
thetaDegrees = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);
uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
backgroundX = fit.backgroundTrajectory.x(t);
backgroundY = fit.backgroundTrajectory.y(t);

observedSpeedSquaredSum = zeros(size(t));
mesoscaleSpeedSquaredSum = zeros(size(t));
submesoscaleSpeedSquaredSum = zeros(size(t));
decompositionXLimits = [min(backgroundX), max(backgroundX)];
decompositionYLimits = [min(backgroundY), max(backgroundY)];

figure(Color="w", Units="points", Position=[50 50 figureWidthMedium 420]);
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

ax1 = nexttile;
plot(ax1, tDays, sigma/f0, LineWidth=1.5);
hold(ax1, "on");
plot(ax1, tDays, zeta/f0, LineWidth=1.5);
ylabel(ax1, "rate (f)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax1, [tDays(1), tDays(end)]);
legend(ax1, "\sigma", "\zeta", Location="best", FontSize=figureLegendSize);
set(ax1, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax1, "on");
title(ax1, "Site 2 gridded streamfunction fit", FontSize=figureTitleSize, FontName=figureFont);

ax2 = nexttile;
plot(ax2, tDays, thetaDegrees, LineWidth=1.5);
ylabel(ax2, "\theta (deg)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax2, [tDays(1), tDays(end)]);
set(ax2, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax2, "on");

ax3 = nexttile;
plot(ax3, tDays, uBackground, LineWidth=1.5);
hold(ax3, "on");
plot(ax3, tDays, vBackground, LineWidth=1.5);
xlabel(ax3, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(ax3, "u_{bg}, v_{bg} (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(ax3, [tDays(1), tDays(end)]);
legend(ax3, "u_{bg}", "v_{bg}", Location="best", FontSize=figureLegendSize);
set(ax3, FontSize=figureAxisTickSize, FontName=figureFont);
box(ax3, "on");

title(tl, "Site 2 parameter time series");

figure(Color="w", Units="points", Position=[50 50 figureWidthOneColumn 300]);
axMeso = axes;
hold(axMeso, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    mesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter);
    plot(axMeso, mesoscale.x(ti)/1000, mesoscale.y(ti)/1000, LineWidth=1.2);
end
axis(axMeso, "equal");
xlabel(axMeso, "x (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axMeso, "y (km)", FontSize=figureAxisLabelSize, FontName=figureFont);
title(axMeso, "Site 2 fixed-frame mesoscale trajectories", FontSize=figureTitleSize, FontName=figureFont);
set(axMeso, FontSize=figureAxisTickSize, FontName=figureFont);
box(axMeso, "on");

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
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    submesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter);
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

title(tlDecomp, "Site 2 fixed-frame decomposition");

figure(Color="w", Units="points", Position=[50 50 figureWidthOneColumn 420]);
tlCom = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axCom = nexttile;
hold(axCom, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xi = trajectory.x(ti);
    yi = trajectory.y(ti);
    [~, qObserved, rObserved] = fit.centeredCoordinates(ti, xi, yi);
    centeredMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter);
    fixedMesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter);
    fixedSubmesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter);

    if iDrifter == 1
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8*[1 1 1], DisplayName="Observed");
        plot(axCom, centeredMesoscale.x(ti)/1000, centeredMesoscale.y(ti)/1000, LineWidth=1.2, DisplayName="Mesoscale fit");
    else
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8*[1 1 1]);
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
backgroundSpeedRms = hypot(fit.backgroundTrajectory.u(t), fit.backgroundTrajectory.v(t));
submesoscaleSpeedRms = sqrt(submesoscaleSpeedSquaredSum/nDrifters);

axRms = nexttile;
plot(axRms, tDays, observedSpeedRms, Color=0.2*[1 1 1], LineWidth=1.5);
hold(axRms, "on");
plot(axRms, tDays, mesoscaleSpeedRms, LineWidth=1.5);
plot(axRms, tDays, backgroundSpeedRms, LineWidth=1.5);
plot(axRms, tDays, submesoscaleSpeedRms, LineWidth=1.5);
xlabel(axRms, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axRms, "RMS speed (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axRms, [tDays(1), tDays(end)]);
legend(axRms, "Observed", "Mesoscale", "Background", "Submesoscale", Location="best", FontSize=figureLegendSize);
set(axRms, FontSize=figureAxisTickSize, FontName=figureFont);
box(axRms, "on");

title(tlCom, "Site 2 COM-frame trajectories and RMS speeds");

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
