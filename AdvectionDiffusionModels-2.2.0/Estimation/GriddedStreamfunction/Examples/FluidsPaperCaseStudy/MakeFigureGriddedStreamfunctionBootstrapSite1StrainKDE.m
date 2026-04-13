scriptDir = fileparts(mfilename("fullpath"));
siteNumber = 1;
siteLabel = "Site " + siteNumber;
dataFilename = "smoothedGriddedRho" + siteNumber + "Drifters.mat";
nBootstraps = 500;
randomSeed = 0;
scoreStride = 6;
psiS = [2 2 1];
fastS = 3;
mesoscaleConstraint = "zeroVorticity";
timeDays = 0:6;
shouldSaveFigures = 0;
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", dataFilename);

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);

figureFont = "Helvetica";
figureWidth = 650;
figureHeight = 760;
titleSize = 10;
axisLabelSize = 10;
tickSize = 9;
annotationSize = 8;
pointColor = [0 0.4470 0.7410];
bandAlpha = 0.22;
selectionColor = 0.85 * [1 1 1];
ringColor = 0.75 * [1 1 1];
pctTarget = (0.8:-0.1:0.1).';
minimumHalfWidth = 1e-6;
wedgeMass = 0.68;
wedgeColor = 0.5 * [1 1 1];
wedgeAlpha = 0.35;
wedgeResolution = 100;

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
    trajectories = TrajectorySpline.empty(0, 1);
    for iDrifter = 1:size(x, 2)
        trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3); %#ok<SAGROW>
    end

    bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
        trajectories, ...
        nBootstraps=nBootstraps, ...
        randomSeed=randomSeed, ...
        queryTimes=t, ...
        scoreStride=scoreStride, ...
        psiS=psiS, ...
        fastS=fastS, ...
        mesoscaleConstraint=mesoscaleConstraint);
    bootstrap.writeToFile(cachePath, shouldOverwriteExisting=true);
end

timeDays = reshape(timeDays, [], 1);
tBootstrapDays = reshape(bootstrap.queryTimes/86400, [], 1);
displayDays = [timeDays; tBootstrapDays(end)];
timeIndices = zeros(size(displayDays));
for iTime = 1:numel(displayDays)
    [~, timeIndices(iTime)] = min(abs(tBootstrapDays - displayDays(iTime)));
end
matchedTimeDays = tBootstrapDays(timeIndices);
disp(table(displayDays, matchedTimeDays, timeIndices, VariableNames=["RequestedDay", "MatchedDay", "QueryIndex"]))

ringRadii = (2:2:10) * 1e-6 / f0;
sigmaN = bootstrap.summary.sigma_n / f0;
sigmaS = bootstrap.summary.sigma_s / f0;
fullSigmaN = reshape(bootstrap.fullSummary.sigma_n, [], 1) / f0;
fullSigmaS = reshape(bootstrap.fullSummary.sigma_s, [], 1) / f0;
allCoordinates = [sigmaN(:); sigmaS(:); fullSigmaN; fullSigmaS];
halfWidth = max(1.1 * max(abs(allCoordinates), [], "all", "omitnan"), minimumHalfWidth);
minimum = [-halfWidth, -halfWidth];
maximum = [halfWidth, halfWidth];

% The 2D panels stay in raw (\sigma_n,\sigma_s) space, but the 1D angle
% summary uses the physical extensional-axis half-angle.
thetaReference = 0.5 * unwrap(atan2(fullSigmaS, fullSigmaN));
statisticsSeries = cell(numel(tBootstrapDays), 1);
strainSummarySeries = cell(numel(tBootstrapDays), 1);
modeSigma = NaN(numel(tBootstrapDays), 1);
sigmaBounds = NaN(numel(tBootstrapDays), 2);
modeThetaDegrees = NaN(numel(tBootstrapDays), 1);
thetaBoundsDegrees = NaN(numel(tBootstrapDays), 2);

for iTime = 1:numel(tBootstrapDays)
    samples = [ ...
        reshape(sigmaN(iTime, :), [], 1), ...
        reshape(sigmaS(iTime, :), [], 1)];
    statisticsSeries{iTime} = KernelDensityEstimate.planarStatisticsFromData( ...
        samples, ...
        referencePoint=[fullSigmaN(iTime), fullSigmaS(iTime)], ...
        minimum=minimum, ...
        maximum=maximum, ...
        summaryMass=wedgeMass);
    strainSummarySeries{iTime} = KernelDensityEstimate.strainSummaryFromPlanarStatistics( ...
        statisticsSeries{iTime}, ...
        thetaReference=thetaReference(iTime));
    modeSigma(iTime) = strainSummarySeries{iTime}.modeSigma;
    sigmaBounds(iTime, :) = strainSummarySeries{iTime}.sigmaBounds;
    modeThetaDegrees(iTime) = strainSummarySeries{iTime}.modeTheta * 180/pi;
    thetaBoundsDegrees(iTime, :) = strainSummarySeries{iTime}.thetaBounds * 180/pi;
end

figure(Color="w", Units="points", Position=[50 50 figureWidth figureHeight]);
tl = tiledlayout(4, 4, TileSpacing="compact", Padding="compact");
cmap = parula(256);
cmap(1, :) = 1;
colormap(cmap);

for iTime = 1:numel(displayDays)
    ax = nexttile(tl, iTime);
    samples = [ ...
        reshape(sigmaN(timeIndices(iTime), :), [], 1), ...
        reshape(sigmaS(timeIndices(iTime), :), [], 1)];
    KernelDensityEstimate.plotPlanarStatistics(ax, statisticsSeries{timeIndices(iTime)}, ...
        samples=samples, ...
        contourMasses=pctTarget, ...
        showWedge=~strainSummarySeries{timeIndices(iTime)}.containsZero, ...
        scatterColor=pointColor, ...
        wedgeColor=wedgeColor, ...
        wedgeAlpha=wedgeAlpha, ...
        wedgeResolution=wedgeResolution, ...
        ringRadii=ringRadii, ...
        ringColor=ringColor, ...
        ringLineStyle=":", ...
        ringLineWidth=0.8);

    if iTime > 4
        xlabel(ax, "\sigma_n (f_0)", FontSize=axisLabelSize, FontName=figureFont);
    else
        ax.XTickLabel = [];
    end
    if mod(iTime - 1, 4) == 0
        ylabel(ax, "\sigma_s (f_0)", FontSize=axisLabelSize, FontName=figureFont);
    else
        ax.YTickLabel = [];
    end
    text(ax, 0.03, 0.97, sprintf("day %.1f", matchedTimeDays(iTime)), ...
        Units="normalized", HorizontalAlignment="left", VerticalAlignment="top", ...
        FontSize=annotationSize, FontName=figureFont, Color=[0.2 0.2 0.2]);
    set(ax, FontSize=tickSize, FontName=figureFont);
end

axMagnitude = nexttile(tl, 9, [1 4]);
bandHandle = fill(axMagnitude, [tBootstrapDays; flipud(tBootstrapDays)], [sigmaBounds(:, 2); flipud(sigmaBounds(:, 1))], ...
    pointColor, EdgeColor="none", FaceAlpha=bandAlpha);
hold(axMagnitude, "on");
lineHandle = plot(axMagnitude, tBootstrapDays, modeSigma, Color=pointColor, LineWidth=1.5);
for iTime = 1:numel(matchedTimeDays)
    xline(axMagnitude, matchedTimeDays(iTime), ":", Color=selectionColor, HandleVisibility="off");
end
ylabel(axMagnitude, "\sigma (f_0)", FontSize=axisLabelSize, FontName=figureFont);
xlim(axMagnitude, [tBootstrapDays(1), tBootstrapDays(end)]);
legend(axMagnitude, [lineHandle, bandHandle], "KDE mode", "68% contour band", Location="northwest", FontSize=tickSize);
set(axMagnitude, FontSize=tickSize, FontName=figureFont);
box(axMagnitude, "on");

axAngle = nexttile(tl, 13, [1 4]);
fill(axAngle, [tBootstrapDays; flipud(tBootstrapDays)], [thetaBoundsDegrees(:, 2); flipud(thetaBoundsDegrees(:, 1))], ...
    pointColor, EdgeColor="none", FaceAlpha=bandAlpha);
hold(axAngle, "on");
plot(axAngle, tBootstrapDays, modeThetaDegrees, Color=pointColor, LineWidth=1.5);
for iTime = 1:numel(matchedTimeDays)
    xline(axAngle, matchedTimeDays(iTime), ":", Color=selectionColor, HandleVisibility="off");
end
xlabel(axAngle, "time (days)", FontSize=axisLabelSize, FontName=figureFont);
ylabel(axAngle, "\theta (deg)", FontSize=axisLabelSize, FontName=figureFont);
xlim(axAngle, [tBootstrapDays(1), tBootstrapDays(end)]);
set(axAngle, FontSize=tickSize, FontName=figureFont);
box(axAngle, "on");

title(tl, siteLabel + " bootstrap strain KDE demo", FontSize=titleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site" + siteNumber + "BootstrapStrainKDE.eps"), "-depsc2");
end
