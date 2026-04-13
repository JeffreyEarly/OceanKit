function results = diagnoseKernelDensityEstimateBandwidth2D()
% diagnoseKernelDensityEstimateBandwidth2D Compare ISE around the selected 2-D bandwidth.

repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
addPackageToPath(repoRoot);

rng(17)
truthCases = diagnosticTruthCases();
gridSize = [401 321];
bandwidthFactors = [1 0.75 1.25];
results = table();

for iCase = 1:numel(truthCases)
    selectedModel = KernelDensityEstimate.fromData(truthCases(iCase).data, minimum=truthCases(iCase).minimum, maximum=truthCases(iCase).maximum);

    for iFactor = 1:numel(bandwidthFactors)
        model = KernelDensityEstimate( ...
            data=truthCases(iCase).data, ...
            bandwidth=selectedModel.bandwidth .* bandwidthFactors(iFactor), ...
            minimum=truthCases(iCase).minimum, ...
            maximum=truthCases(iCase).maximum);
        [density, gridVectors] = model.densityOnGrid(gridSize=gridSize);
        truthDensity = truthCases(iCase).densityOnGrid(gridVectors{1}, gridVectors{2});
        squaredError = (density - truthDensity).^2;
        integratedSquaredError = trapz(gridVectors{2}, trapz(gridVectors{1}, squaredError, 1), 2);
        resultRow = table( ...
            string(truthCases(iCase).name), ...
            bandwidthLabel(bandwidthFactors(iFactor)), ...
            model.bandwidth(1), ...
            model.bandwidth(2), ...
            integratedSquaredError, ...
            'VariableNames', ["caseName", "bandwidthLabel", "bandwidthX", "bandwidthY", "integratedSquaredError"]);
        results = [results; resultRow]; %#ok<AGROW>
    end
end

if nargout == 0
    disp(results)
end
end

function truthCases = diagnosticTruthCases()
nSamples = 400;

singleMean = [0.25 -0.10];
singleSigma = [1.10 0.35];
singleData = [singleMean(1) + singleSigma(1) * randn(nSamples, 1), singleMean(2) + singleSigma(2) * randn(nSamples, 1)];
singleMinimum = singleMean - 4 * singleSigma;
singleMaximum = singleMean + 4 * singleSigma;

mixtureWeights = [0.55 0.45];
mixtureMeans = [-1.8 -0.5; 2.1 1.2];
mixtureSigmas = [0.55 0.30; 0.90 0.60];
nFirst = round(mixtureWeights(1) * nSamples);
nSecond = nSamples - nFirst;
mixtureData = [ ...
    mixtureMeans(1, 1) + mixtureSigmas(1, 1) * randn(nFirst, 1), mixtureMeans(1, 2) + mixtureSigmas(1, 2) * randn(nFirst, 1); ...
    mixtureMeans(2, 1) + mixtureSigmas(2, 1) * randn(nSecond, 1), mixtureMeans(2, 2) + mixtureSigmas(2, 2) * randn(nSecond, 1)];
mixtureMinimum = min(mixtureMeans - 4 * mixtureSigmas, [], 1);
mixtureMaximum = max(mixtureMeans + 4 * mixtureSigmas, [], 1);

truthCases = [ ...
    struct( ...
        "name", "Single anisotropic Gaussian", ...
        "data", singleData, ...
        "minimum", singleMinimum, ...
        "maximum", singleMaximum, ...
        "densityOnGrid", @(xVector, yVector) gaussianDensityOnGrid(singleMean, singleSigma, xVector, yVector)), ...
    struct( ...
        "name", "Axis-aligned Gaussian mixture", ...
        "data", mixtureData, ...
        "minimum", mixtureMinimum, ...
        "maximum", mixtureMaximum, ...
        "densityOnGrid", @(xVector, yVector) gaussianMixtureDensityOnGrid(mixtureWeights, mixtureMeans, mixtureSigmas, xVector, yVector))];
end

function label = bandwidthLabel(factor)
if factor == 1
    label = "selected";
else
    label = string(sprintf("%.2fx", factor));
end
end

function density = gaussianDensityOnGrid(meanValue, sigmaValue, xVector, yVector)
[X, Y] = ndgrid(xVector, yVector);
density = exp(-0.5 * (((X - meanValue(1))/sigmaValue(1)).^2 + ((Y - meanValue(2))/sigmaValue(2)).^2)) ...
    /(2*pi * sigmaValue(1) * sigmaValue(2));
end

function density = gaussianMixtureDensityOnGrid(weights, means, sigmas, xVector, yVector)
density = zeros(numel(xVector), numel(yVector));
for iComponent = 1:numel(weights)
    density = density + weights(iComponent) * gaussianDensityOnGrid(means(iComponent, :), sigmas(iComponent, :), xVector, yVector);
end
end

function addPackageToPath(repoRoot)
metadataPath = fullfile(repoRoot, "resources", "mpackage.json");
if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    for iFolder = 1:length(metadata.folders)
        addpath(fullfile(repoRoot, metadata.folders(iFolder).path));
    end
end

addpath(repoRoot);
end
