% Show a simple Site 1 fit using the current gridded-streamfunction estimator.
scriptDir = fileparts(mfilename("fullpath"));
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", "smoothedGriddedRho1Drifters.mat");
if ~isfile(dataPath)
    error("GriddedStreamfunction:MissingExampleData", "Expected local example data at %s.", dataPath);
end

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
nDrifters = size(x, 2);

trajectoryCell = cell(nDrifters, 1);
for iDrifter = 1:nDrifters
    trajectoryCell{iDrifter} = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
end
trajectories = vertcat(trajectoryCell{:});

fit = GriddedStreamfunction.fromTrajectories(trajectories, psiS=[2 2 3]);

tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);

% Use the first drifter as the default one-trajectory example.
iTrajectory = 1;
trajectory = fit.observedTrajectories(iTrajectory);
ti = trajectory.t;
xi = trajectory.x(ti);
yi = trajectory.y(ti);
uObserved = trajectory.u(ti);
vObserved = trajectory.v(ti);
uMesoscale = fit.uMesoscale(ti, xi, yi);
vMesoscale = fit.vMesoscale(ti, xi, yi);
uBackgroundTrajectory = fit.uBackground(ti);
vBackgroundTrajectory = fit.vBackground(ti);
uSubmesoscale = uObserved - uMesoscale - uBackgroundTrajectory;
vSubmesoscale = vObserved - vMesoscale - vBackgroundTrajectory;
uReconstruction = uMesoscale + uBackgroundTrajectory + uSubmesoscale;
vReconstruction = vMesoscale + vBackgroundTrajectory + vSubmesoscale;
uDecompositionError = max(abs(uObserved - uReconstruction));
vDecompositionError = max(abs(vObserved - vReconstruction));
if uDecompositionError > 1e-10 || vDecompositionError > 1e-10
    error("GriddedStreamfunction:DecompositionMismatch", ...
        "The raw decomposition failed to reconstruct the observed velocity for drifter %d.", iTrajectory);
end

fprintf("Site 1 gridded streamfunction fit\n");
fprintf("  drifters: %d\n", nDrifters);
fprintf("  samples:  %d\n", numel(t));
fprintf("  mean sigma/f0: %.4f\n", mean(sigma/f0));
fprintf("  mean zeta/f0:  %.4f\n", mean(zeta/f0));
fprintf("  mean u_bg:     %.4f m/s\n", mean(uBackground));
fprintf("  mean v_bg:     %.4f m/s\n", mean(vBackground));
fprintf("  drifter %d max |u-u_{recon}|: %.3e m/s\n", iTrajectory, uDecompositionError);
fprintf("  drifter %d max |v-v_{recon}|: %.3e m/s\n", iTrajectory, vDecompositionError);

figure(Color="w");
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

ax1 = nexttile;
plot(ax1, tDays, sigma/f0, LineWidth=2);
hold(ax1, "on");
plot(ax1, tDays, zeta/f0, LineWidth=2);
ylabel(ax1, "\sigma/f_0, \zeta/f_0");
xlim(ax1, [tDays(1), tDays(end)]);
legend(ax1, "\sigma/f_0", "\zeta/f_0", Location="best");
box(ax1, "on");
title(ax1, "Site 1 gridded streamfunction fit");

ax2 = nexttile;
plot(ax2, tDays, uBackground, LineWidth=2);
hold(ax2, "on");
plot(ax2, tDays, vBackground, LineWidth=2);
ylabel(ax2, "u_{bg}, v_{bg} (m/s)");
xlim(ax2, [tDays(1), tDays(end)]);
legend(ax2, "u_{bg}", "v_{bg}", Location="best");
box(ax2, "on");

ax3 = nexttile;
hold(ax3, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xi = trajectory.x(ti);
    yi = trajectory.y(ti);
    [~, qi, ri] = fit.centeredCoordinates(ti, xi, yi);
    qMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).x(ti);
    rMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).y(ti);

    if iDrifter == 1
        plot(ax3, qi/1000, ri/1000, Color=0.8*[1 1 1], DisplayName="Observed");
        plot(ax3, qMesoscale/1000, rMesoscale/1000, LineWidth=1.5, DisplayName="Mesoscale fit");
    else
        plot(ax3, qi/1000, ri/1000, Color=0.8*[1 1 1]);
        plot(ax3, qMesoscale/1000, rMesoscale/1000, LineWidth=1.5);
    end
end
axis(ax3, "equal");
xlabel(ax3, "q (km)");
ylabel(ax3, "r (km)");
legend(ax3, "Observed", "Mesoscale fit", Location="best");
box(ax3, "on");

title(tl, "GriddedStreamfunction Site 1");

figure(Color="w");
tlDecomposition = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axU = nexttile;
hold(axU, "on");
scatter(axU, ti/86400, uObserved, 5^2, "k", "filled", DisplayName="observed");
plot(axU, ti/86400, uReconstruction, "k", LineWidth=1.5, DisplayName="reconstruction");
plot(axU, ti/86400, uMesoscale, LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale");
plot(axU, ti/86400, uBackgroundTrajectory, LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background");
plot(axU, ti/86400, uSubmesoscale, LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale");
ylabel(axU, "u (m/s)");
xlim(axU, [tDays(1), tDays(end)]);
legend(axU, Location="best");
box(axU, "on");

axV = nexttile;
hold(axV, "on");
scatter(axV, ti/86400, vObserved, 5^2, "k", "filled", DisplayName="observed");
plot(axV, ti/86400, vReconstruction, "k", LineWidth=1.5, DisplayName="reconstruction");
plot(axV, ti/86400, vMesoscale, LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale");
plot(axV, ti/86400, vBackgroundTrajectory, LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background");
plot(axV, ti/86400, vSubmesoscale, LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale");
xlabel(axV, "time (days)");
ylabel(axV, "v (m/s)");
xlim(axV, [tDays(1), tDays(end)]);
legend(axV, Location="best");
box(axV, "on");

title(tlDecomposition, sprintf("Site 1 drifter %d raw velocity decomposition", iTrajectory));
