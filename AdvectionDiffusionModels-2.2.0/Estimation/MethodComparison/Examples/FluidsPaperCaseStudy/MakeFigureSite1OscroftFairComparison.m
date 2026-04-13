scriptDir = locateExampleScriptDirectory();
shouldSaveFigures = 0;
scaleFactor = 1;

requireMethodComparisonFunction("buildSite1OscroftFairComparison")

comparison = buildSite1OscroftFairComparison();
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

tDays = reshape(comparison.site.queryTimes / 86400, [], 1);
f0 = comparison.site.f0;
oldFair = comparison.old.shared;
newFair = comparison.new.shared;
oldStrain = kdeStrainSeriesFromSharedSummary(oldFair);
newStrain = kdeStrainSeriesFromSharedSummary(newFair);

oldColor = [0.8500 0.3250 0.0980];
newColor = [0 0.4470 0.7410];
lineWidth = 1.5;
bandAlpha = 0.22;

figure(Color="w", Units="points", Position=[50 50 figureWidthMedium 420]);
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

axSigma = nexttile;
plotBoundedSeries(axSigma, tDays, oldStrain.modeSigma/f0, oldStrain.sigmaBounds/f0, oldColor, "-", bandAlpha, lineWidth);
plotBoundedSeries(axSigma, tDays, newStrain.modeSigma/f0, newStrain.sigmaBounds/f0, newColor, "-", bandAlpha, lineWidth);
ylabel(axSigma, "\sigma (f_0)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axSigma, [tDays(1), tDays(end)]);
set(axSigma, FontSize=figureAxisTickSize, FontName=figureFont);
box(axSigma, "on");
title(axSigma, "Site 1 fair bootstrap comparison", FontSize=figureTitleSize, FontName=figureFont);

axTheta = nexttile;
plotBoundedSeries(axTheta, tDays, oldStrain.modeThetaDegrees, oldStrain.thetaBoundsDegrees, oldColor, "-", bandAlpha, lineWidth);
plotBoundedSeries(axTheta, tDays, newStrain.modeThetaDegrees, newStrain.thetaBoundsDegrees, newColor, "-", bandAlpha, lineWidth);
ylabel(axTheta, "\theta (deg)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axTheta, [tDays(1), tDays(end)]);
ylim(axTheta, [-100 100]);
yticks(axTheta, [-90 -45 0 45 90]);
set(axTheta, FontSize=figureAxisTickSize, FontName=figureFont);
box(axTheta, "on");

axVelocity = nexttile;
hold(axVelocity, "on");
oldU = plotComparisonSeries(axVelocity, tDays, oldFair.commonSummary.uCenter, oldFair.scores.common, oldColor, "-");
oldV = plotComparisonSeries(axVelocity, tDays, oldFair.commonSummary.vCenter, oldFair.scores.common, oldColor, "--");
newU = plotComparisonSeries(axVelocity, tDays, newFair.summary.uCenter, newFair.scores.common, newColor, "-");
newV = plotComparisonSeries(axVelocity, tDays, newFair.summary.vCenter, newFair.scores.common, newColor, "--");
xlabel(axVelocity, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axVelocity, "u_c/v_c (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axVelocity, [tDays(1), tDays(end)]);
ylim(axVelocity, [-0.25 0.25]);
legend(axVelocity, [oldU, oldV, newU, newV], ...
    "Oscroft u_c", "Oscroft v_c", "Gridded u_c", "Gridded v_c", ...
    Location="southwest", FontSize=figureLegendSize, NumColumns=2);
set(axVelocity, FontSize=figureAxisTickSize, FontName=figureFont);
box(axVelocity, "on");

title(tl, "Site 1 shared-bootstrap sanity check", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site1OscroftFairComparison.eps"), "-depsc2");
end

function lineHandle = plotBoundedSeries(ax, tDays, centerValues, bounds, baseColor, lineStyle, bandAlpha, lineWidth)
hold(ax, "on");
drawBoundedFill(ax, tDays, bounds(:, 1), bounds(:, 2), baseColor, bandAlpha);
lineHandle = plot(ax, tDays, centerValues, Color=baseColor, LineStyle=lineStyle, LineWidth=lineWidth);
end

function drawBoundedFill(ax, tDays, lowerValues, upperValues, baseColor, bandAlpha)
isFinite = isfinite(tDays) & isfinite(lowerValues) & isfinite(upperValues);
if ~any(isFinite)
    return
end

iStart = 1;
while iStart <= numel(tDays)
    while iStart <= numel(tDays) && ~isFinite(iStart)
        iStart = iStart + 1;
    end
    if iStart > numel(tDays)
        break
    end

    iEnd = iStart;
    while iEnd < numel(tDays) && isFinite(iEnd + 1)
        iEnd = iEnd + 1;
    end

    fill(ax, [tDays(iStart:iEnd); flipud(tDays(iStart:iEnd))], [upperValues(iStart:iEnd); flipud(lowerValues(iStart:iEnd))], ...
        baseColor, EdgeColor="none", FaceAlpha=bandAlpha, HandleVisibility="off");
    iStart = iEnd + 1;
end
end

function medianLine = plotComparisonSeries(ax, tDays, values, scores, baseColor, lineStyle)
rankFractions = [0.90 0.68];
rankAlphas = [0.08 0.14];
quantileAlpha = 0.22;

hold(ax, "on");
for iRank = 1:numel(rankFractions)
    indices = rankedIndices(scores, rankFractions(iRank));
    envelope = rankedEnvelope(values(:, indices));
    fill(ax, [tDays; flipud(tDays)], envelope, baseColor, ...
        EdgeColor="none", FaceAlpha=rankAlphas(iRank), HandleVisibility="off");
end

quantiles = quantile(values, [0.16 0.5 0.84], 2);
fill(ax, [tDays; flipud(tDays)], [quantiles(:, 3); flipud(quantiles(:, 1))], baseColor, ...
    EdgeColor="none", FaceAlpha=quantileAlpha, HandleVisibility="off");
medianLine = plot(ax, tDays, quantiles(:, 2), Color=baseColor, LineStyle=lineStyle, LineWidth=1.5);
end

function strainSeries = kdeStrainSeriesFromSharedSummary(sharedComparison)
summary = sharedComparison.commonSummary;
fullSummary = sharedComparison.fullCommonSummary;
fullSigmaN = reshape(fullSummary.sigma_n, [], 1);
fullSigmaS = reshape(fullSummary.sigma_s, [], 1);
% Use the physical extensional-axis half-angle rather than the older
% principal-axis display helper.
thetaReference = 0.5 * unwrap(atan2(fullSigmaS, fullSigmaN));
nTime = size(summary.sigma_n, 1);

modeSigma = NaN(nTime, 1);
sigmaBounds = NaN(nTime, 2);
modeThetaDegrees = NaN(nTime, 1);
thetaBoundsDegrees = NaN(nTime, 2);

for iTime = 1:nTime
    samples = [reshape(summary.sigma_n(iTime, :), [], 1), reshape(summary.sigma_s(iTime, :), [], 1)];
    referencePoint = [reshape(fullSummary.sigma_n(iTime), 1, 1), reshape(fullSummary.sigma_s(iTime), 1, 1)];
    statistics = KernelDensityEstimate.planarStatisticsFromData(samples, referencePoint=referencePoint);
    if isfinite(thetaReference(iTime))
        strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics, thetaReference=thetaReference(iTime));
    else
        strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);
    end

    modeSigma(iTime) = strainSummary.modeSigma;
    sigmaBounds(iTime, :) = strainSummary.sigmaBounds;
    modeThetaDegrees(iTime) = strainSummary.modeTheta * 180/pi;
    thetaBoundsDegrees(iTime, :) = strainSummary.thetaBounds * 180/pi;
end

strainSeries = struct( ...
    "modeSigma", modeSigma, ...
    "sigmaBounds", sigmaBounds, ...
    "modeThetaDegrees", modeThetaDegrees, ...
    "thetaBoundsDegrees", thetaBoundsDegrees);
end

function indices = rankedIndices(scores, fraction)
[~, rankedIndices] = sort(scores, "descend");
indices = rankedIndices(1:max(1, round(fraction * numel(scores))));
end

function envelope = rankedEnvelope(values)
envelope = [max(values, [], 2); flipud(min(values, [], 2))];
end

function scriptDir = locateExampleScriptDirectory()
scriptPath = string(mfilename("fullpath"));
if strlength(scriptPath) == 0
    scriptName = string(mfilename);
    if strlength(scriptName) > 0
        scriptPath = string(which(scriptName));
    end
end

if strlength(scriptPath) == 0 && usejava("desktop")
    scriptPath = string(matlab.desktop.editor.getActiveFilename);
end

if strlength(scriptPath) == 0
    error("GriddedStreamfunction:UnableToLocateExampleScript", ...
        "Could not determine the example script location. Run the saved script file or add Estimation/MethodComparison to the MATLAB path first.");
end

scriptDir = char(fileparts(scriptPath));
end

function requireMethodComparisonFunction(functionName)
if exist(functionName, "file") == 2
    return
end

error("GriddedStreamfunction:MissingMethodComparisonDependency", ...
    "%s is not on the MATLAB path. Reload the AdvectionDiffusionModels package path so Estimation/MethodComparison is available.", functionName);
end
