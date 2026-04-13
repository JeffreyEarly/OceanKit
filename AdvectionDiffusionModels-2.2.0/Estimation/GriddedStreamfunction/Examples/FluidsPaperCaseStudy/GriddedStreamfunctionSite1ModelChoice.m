scriptDir = fileparts(mfilename("fullpath"));
siteNumber = 2;
dataFilename = "smoothedGriddedRho" + siteNumber + "Drifters.mat";
nBootstraps = 100;
randomSeed = 0;
scoreStride = 6;
fastS = 3;
psiSTimeValues = 0:5;
mesoscaleConstraints = ["zeroVorticity", "none"];
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", dataFilename);
jlabPath = fullfile(scriptDir, "..", "..", "..", "..", "..", "jLab");

addpath(genpath(jlabPath));

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
queryTimes = t;

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:size(x, 2)
    trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3); %#ok<SAGROW>
end

cacheDirectory = fullfile(scriptDir, "BootstrapData");
if exist(cacheDirectory, "dir") == 0
    mkdir(cacheDirectory);
end

nCandidateModels = numel(psiSTimeValues) * numel(mesoscaleConstraints);
psiSColumn = strings(nCandidateModels, 1);
mesoscaleConstraintColumn = strings(nCandidateModels, 1);
mesoscaleDegreesOfFreedomColumn = NaN(nCandidateModels, 1);
fullFitKappaColumn = NaN(nCandidateModels, 1);
fullFitCoherenceColumn = NaN(nCandidateModels, 1);
bestFitKappaColumn = NaN(nCandidateModels, 1);
bestFitCoherenceColumn = NaN(nCandidateModels, 1);

iModel = 0;
for iConstraint = 1:numel(mesoscaleConstraints)
    mesoscaleConstraint = mesoscaleConstraints(iConstraint);
    for iPsiSTime = 1:numel(psiSTimeValues)
        iModel = iModel + 1;
        psiS = [3 3 psiSTimeValues(iPsiSTime)];
        psiSTag = join(string(psiS), "-");
        psiSLabel = "[" + join(string(psiS), " ") + "]";
        cacheFilename = "Rho" + siteNumber + ...
            "GriddedStreamfunctionBootstrapFits" + nBootstraps + ...
            "_seed" + randomSeed + ...
            "_stride" + scoreStride + ...
            "_fastS" + fastS + ...
            "_psiS" + psiSTag + ...
            "_mesoscaleConstraint-" + mesoscaleConstraint + ".nc";
        cachePath = fullfile(cacheDirectory, cacheFilename);

        fprintf("Site 1 model choice: psiS=%s, mesoscaleConstraint=%s\n", psiSLabel, mesoscaleConstraint);
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

        psiSColumn(iModel) = psiSLabel;
        mesoscaleConstraintColumn(iModel) = mesoscaleConstraint;
        mesoscaleDegreesOfFreedomColumn(iModel) = bootstrap.mesoscaleDegreesOfFreedom;
        fullFitKappaColumn(iModel) = bootstrap.fullFitKappa;
        fullFitCoherenceColumn(iModel) = bootstrap.fullFitCoherence;
        bestFitKappaColumn(iModel) = bootstrap.bestFitKappa;
        bestFitCoherenceColumn(iModel) = bootstrap.bestFitCoherence;
    end
end

diagnosticsTable = table( ...
    psiSColumn, ...
    mesoscaleConstraintColumn, ...
    mesoscaleDegreesOfFreedomColumn, ...
    fullFitKappaColumn, ...
    fullFitCoherenceColumn, ...
    bestFitKappaColumn, ...
    bestFitCoherenceColumn, ...
    VariableNames=["psiS", "mesoscaleConstraint", "mesoscaleDegreesOfFreedom", "fullFitKappa", "fullFitCoherence", "bestFitKappa", "bestFitCoherence"]);
diagnosticsTable = sortrows(diagnosticsTable, ["mesoscaleConstraint", "psiS"]);

disp(diagnosticsTable);
