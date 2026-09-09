%% Hydrostatic mode accuracy fix comparison
% This diagnostic separates quadrature-node error from numerical-mode error
% in hydrostatic Parseval tests. Analytical constant and exponential modes
% are used only as references; the numerical production path remains the
% WKB or ordinary spectral solver.

nPointsList = [32 64];
nEVPList = [128 160 192 224 256 320];
latitude = 31;
g = 9.81;
shouldMakeFigures = true;
shouldRunNodeWeightOptimization = true;
nodeWeightOptimizationPointCounts = 32;
shouldRunRetainedCapacitySweep = true;
shouldRunPhaseShiftScan = true;
phaseShiftCount = 41;

profiles = [
    hydrostaticProfile("constant-limit",4000,3*2*pi/3600,Inf)
    hydrostaticProfile("exponential",4000,3*2*pi/3600,1300)
    ];
solverKinds = ["wkb", "spectral"];

for iProfile = 1:length(profiles)
    profile = profiles(iProfile);
    fprintf('\n%s profile\n', profile.name);

    for nPoints = nPointsList
        nEVP = max(nEVPList(nEVPList >= ceil(2.1*nPoints)));
        nEVP = nEVP(1);
        fprintf('\nN = %d quadrature points, nEVP = %d\n', nPoints, nEVP);

        referenceNodes = referenceQuadratureNodes(profile, nPoints, latitude, g);
        referenceData = modeData(profile, "analytical", referenceNodes, nEVP, latitude, g);
        rawNodes = numericalQuadratureNodes(profile, "wkb", nPoints, nEVP, latitude, g);
        polishedNodes = polishedNumericalQuadratureNodes(profile, "wkb", nPoints, nEVP, latitude, g);

        shouldOptimizeNodes = profile.name == "exponential" && shouldRunNodeWeightOptimization;
        shouldOptimizeNodes = shouldOptimizeNodes && any(nPoints == nodeWeightOptimizationPointCounts) && exist("fmincon", "file") == 2;
        optimizationRetainedBands = unique([nPoints - 2, max(1, round(nPoints/2) - 1)]);
        nOptimizationRows = double(shouldOptimizeNodes)*numel(optimizationRetainedBands);
        rows = repmat(emptyComparisonRow(), 1, 2 + 3*numel(solverKinds) + nOptimizationRows);
        rows(1) = comparisonRow(profile, "analytical + reference roots", "analytical", referenceNodes, referenceNodes, referenceData, nEVP, latitude, g);
        rows(2) = comparisonRow(profile, "analytical + raw numerical roots", "analytical", rawNodes, referenceNodes, referenceData, nEVP, latitude, g);
        iRow = 2;
        for solverKind = solverKinds
            label = solverKind + " modes + reference roots";
            rows(iRow+1) = comparisonRow(profile, label, solverKind, referenceNodes, referenceNodes, referenceData, nEVP, latitude, g);
            rows(iRow+2) = comparisonRow(profile, solverKind + " modes + raw roots", solverKind, rawNodes, referenceNodes, referenceData, nEVP, latitude, g);
            label = solverKind + " modes + polished roots";
            rows(iRow+3) = comparisonRow(profile, label, solverKind, polishedNodes, referenceNodes, referenceData, nEVP, latitude, g);
            iRow = iRow + 3;
        end
        if shouldOptimizeNodes
            for nRetained = optimizationRetainedBands
                rows(iRow+1) = nodeWeightOptimizationRow(profile, referenceNodes, referenceData, latitude, g, nRetained);
                iRow = iRow + 1;
            end
        elseif profile.name == "exponential" && shouldRunNodeWeightOptimization && exist("fmincon", "file") ~= 2
            fprintf('  Skipping node/weight optimization because fmincon is unavailable.\n');
        end

        printComparisonRows(rows);
        printRecommendationRows(rows);
        if profile.name == "exponential" && shouldRunPhaseShiftScan
            printPhaseShiftScan(profile, nPoints, referenceNodes, nEVP, latitude, g, phaseShiftCount);
        end
        if profile.name == "exponential" && shouldRunRetainedCapacitySweep
            printRetainedCapacitySweep(profile, referenceNodes, nEVP, latitude, g);
        end
    end

    fprintf('\nnEVP sweep, N = %d\n', nPointsList(end));
    sweepRows = nEVPSweep(profile, nPointsList(end), nEVPList, latitude, g);
    printSweepRows(sweepRows);

    if shouldMakeFigures
        figure(Name="Hydrostatic mode accuracy nEVP sweep");
        semilogy([sweepRows.nEVP], [sweepRows.maxESpec], "o-")
        grid on
        xlabel("nEVP")
        ylabel("max(E_{spec,F}, E_{spec,G})")
        title(sprintf("%s profile, N = %d", profile.name, nPointsList(end)), Interpreter="none")
    end
end

function profile = hydrostaticProfile(name, depth, N0, b)
profile.name = string(name);
profile.depth = depth;
profile.N0 = N0;
profile.b = b;
profile.zIn = [-depth 0];
if isfinite(b)
    profile.N2Function = @(z) N0*N0*exp(2*z/b);
else
    profile.N2Function = @(z) N0*N0*ones(size(z));
end
end

function z = referenceQuadratureNodes(profile, nPoints, latitude, g)
if profile.name == "constant-limit"
    z = linspace(-profile.depth, 0, nPoints).';
    return
end

rootMode = nPoints - 1;
zScan = linspace(-profile.depth, 0, max(4000, 80*nPoints)).';
reference = analyticalSolver(profile, zScan, rootMode, latitude, g);
[~, ~, h] = reference.modesAtFrequency(0);
[~, GFunction] = reference.ModeFunctionsForOmegaAndC(0, sqrt(g*h(rootMode)));
rootFunction = @(z) GFunction(z, 0, sqrt(g*h(rootMode)));
values = rootFunction(zScan);
brackets = find(values(1:end-1).*values(2:end) <= 0);
roots = zeros(numel(brackets), 1);
for iRoot = 1:numel(brackets)
    roots(iRoot) = fzero(rootFunction, zScan(brackets(iRoot):brackets(iRoot)+1));
end
roots(abs(roots + profile.depth) < 1e-9*profile.depth | abs(roots) < 1e-9*profile.depth) = [];
z = unique([profile.zIn(1); roots; profile.zIn(2)], "stable");
if numel(z) ~= nPoints
    error("HydrostaticModeAccuracy:ReferenceRoots", "Expected %d reference nodes, found %d.", nPoints, numel(z));
end
end

function z = numericalQuadratureNodes(profile, solverKind, nPoints, nEVP, latitude, g)
zReference = linspace(profile.zIn(1), profile.zIn(2), max(1024, 12*nPoints)).';
solver = numericalSolver(profile, solverKind, zReference, nEVP, nPoints - 1, latitude, g);
z = solver.GaussQuadraturePointsForModesAtFrequency(nPoints, 0);
end

function z = polishedNumericalQuadratureNodes(profile, solverKind, nPoints, nEVP, latitude, g)
zReference = linspace(profile.zIn(1), profile.zIn(2), max(1024, 12*nPoints)).';
solver = numericalSolver(profile, solverKind, zReference, nEVP, nPoints - 1, latitude, g);
[A, B] = solver.eigenmatricesForFrequency(0);
[V, D] = eig(A, B);
[h, order] = sort(real(solver.hFromLambda(diag(D))), "descend");
V = V(:, order);
resolvedModes = ceil(find(h > 0, 1, "last")/2);
if isempty(resolvedModes) || resolvedModes < nPoints
    error("HydrostaticModeAccuracy:NeedMoreModes", "Only %d modes resolved with nEVP = %d.", resolvedModes, nEVP);
end

rootMode = nPoints - 1;
rootCoefficients = V(:, rootMode);
interiorRoots = InternalModesSpectral.FindRootsFromChebyshevVector(rootCoefficients, solver.xDomain);
interiorRoots(interiorRoots <= solver.xMin | interiorRoots >= solver.xMax) = [];
xRoots = unique([solver.xMin; interiorRoots; solver.xMax], "stable");
while numel(xRoots) > nPoints
    xRoots = sort(xRoots);
    F = InternalModesSpectral.IntegrateChebyshevVector(rootCoefficients);
    value = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(xRoots, solver.xDomain, F);
    [~, minIndex] = min(abs(diff(value)));
    xRoots(minIndex+1) = [];
end

polishedXRoots = xRoots;
rootFunction = @(x) InternalModesSpectral.ValueOfFunctionAtPointOnGrid(x, solver.xDomain, rootCoefficients);
for iRoot = 2:(numel(xRoots)-1)
    left = 0.5*(xRoots(iRoot-1) + xRoots(iRoot));
    right = 0.5*(xRoots(iRoot) + xRoots(iRoot+1));
    try
        if rootFunction(left)*rootFunction(right) < 0
            polishedXRoots(iRoot) = fzero(rootFunction, [left right]);
        else
            polishedXRoots(iRoot) = fzero(rootFunction, xRoots(iRoot));
        end
    catch
        polishedXRoots(iRoot) = xRoots(iRoot);
    end
end
z = InternalModesSpectral.fInverseBisection(solver.x_function, polishedXRoots, min(profile.zIn), max(profile.zIn), 1e-12);
z = reshape(z, [], 1);
end

function row = comparisonRow(profile, label, solverKind, z, referenceNodes, referenceData, nEVP, latitude, g)
data = modeData(profile, solverKind, z, nEVP, latitude, g);
freeDeltaZ = leastSquaresWeights(data);
freeMetrics = parsevalMetrics(data, freeDeltaZ);
positiveMetrics = positiveParsevalMetrics(data, freeDeltaZ);
[modeErrorF, modeErrorG, hError] = referenceErrors(data, referenceData);

row = diagnosticRow(label, solverKind, freeDeltaZ, freeMetrics, positiveMetrics, z, referenceNodes, modeErrorF, modeErrorG, hError);
row.nRetained = size(data.PhiG, 2);
end

function row = nodeWeightOptimizationRow(profile, initialNodes, referenceData, latitude, g, nRetained)
fprintf('  Optimizing exponential analytical nodes and weights for r = %d from Bessel roots...\n', nRetained);
evaluator = exponentialModeEvaluator(profile, initialNodes, latitude, g);
initialData = retainModeBand(modeDataFromEvaluator(profile, evaluator, initialNodes, g), nRetained);
initialDeltaZ = leastSquaresWeights(initialData);
initialPositiveMetrics = positiveParsevalMetrics(initialData, initialDeltaZ);
if any(~isfinite(initialPositiveMetrics.deltaZ))
    error("HydrostaticModeAccuracy:PositiveWeightsUnavailable", "Node optimization requires an initial nonnegative weight vector.");
end
initialDeltaZ = initialPositiveMetrics.deltaZ;
referenceData = retainModeBand(referenceData, nRetained);
nInterior = numel(initialNodes) - 2;
initialVariables = [initialNodes(2:end-1); initialDeltaZ];

minimumSpacing = profile.depth/(100*numel(initialNodes)^2);
lowerBounds = [profile.zIn(1)*ones(nInterior, 1); zeros(numel(initialNodes), 1)];
upperBounds = [profile.zIn(2)*ones(nInterior, 1); profile.depth*ones(numel(initialNodes), 1)];
lowerBounds(1:nInterior) = lowerBounds(1:nInterior) + minimumSpacing;
upperBounds(1:nInterior) = upperBounds(1:nInterior) - minimumSpacing;
A = zeros(max(nInterior - 1, 0), numel(initialVariables));
b = -minimumSpacing*ones(size(A, 1), 1);
for iNode = 1:size(A, 1)
    A(iNode, iNode) = 1;
    A(iNode, iNode+1) = -1;
end

objective = @(variables) nodeWeightObjective(variables, evaluator, profile, g, nRetained);
options = optimoptions("fmincon", Algorithm="sqp", Display="off", MaxIterations=200, MaxFunctionEvaluations=8000, StepTolerance=1e-10);
[variables, objectiveValue, exitFlag, output] = fmincon(objective, initialVariables, A, b, [], [], lowerBounds, upperBounds, [], options);
if exitFlag <= 0
    fprintf('  Node/weight optimization ended with fmincon exit flag %d after %d evaluations.\n', exitFlag, output.funcCount);
end
fprintf('  Final normalized Gram residual objective: %.3e\n', objectiveValue);

[z, deltaZ] = unpackNodeWeightVariables(variables, profile);
data = retainModeBand(modeDataFromEvaluator(profile, evaluator, z, g), nRetained);
metrics = parsevalMetrics(data, deltaZ);
positiveMetrics = metricsWithWeights(metrics, deltaZ);
[modeErrorF, modeErrorG, hError] = referenceErrors(data, referenceData);
label = sprintf("analytical + optimized nodes/weights r=%d", nRetained);
row = diagnosticRow(label, "analytical", deltaZ, metrics, positiveMetrics, z, initialNodes, modeErrorF, modeErrorG, hError);
row.nRetained = nRetained;
end

function printPhaseShiftScan(profile, nPoints, referenceNodes, nEVP, latitude, g, phaseShiftCount)
retainedBands = unique([nPoints - 2, max(1, round(nPoints/2) - 1)]);
thetaValues = linspace(-pi/2, pi/2, phaseShiftCount);
fprintf('\n  Phase-shifted Bessel-root grids\n');
fprintf('  %8s %12s %12s %12s %12s %12s %12s\n', "retained", "best theta", "max E", "free E", "E_F", "E_G", "min dz");
for nRetained = retainedBands
    bestRow.maxESpec = Inf;
    bestRow.freeMaxESpec = NaN;
    bestRow.ESpecF = NaN;
    bestRow.ESpecG = NaN;
    bestRow.minDeltaZ = NaN;
    bestTheta = NaN;
    for theta = thetaValues
        try
            z = phaseShiftedReferenceNodes(profile, nPoints, theta, latitude, g);
            data = modeData(profile, "analytical", z, nEVP, latitude, g);
            data = retainModeBand(data, nRetained);
            deltaZ = leastSquaresWeights(data);
            freeMetrics = parsevalMetrics(data, deltaZ);
            positiveMetrics = positiveParsevalMetrics(data, deltaZ);
            if positiveMetrics.maxESpec < bestRow.maxESpec
                bestRow.maxESpec = positiveMetrics.maxESpec;
                bestRow.freeMaxESpec = max(freeMetrics.ESpecF, freeMetrics.ESpecG);
                bestRow.ESpecF = positiveMetrics.ESpecF;
                bestRow.ESpecG = positiveMetrics.ESpecG;
                bestRow.nodeError = max(abs(z(:) - referenceNodes(:)));
                bestRow.minDeltaZ = min(positiveMetrics.deltaZ);
                bestTheta = theta;
            end
        catch
        end
    end
    fprintf('  %8d %12.3f %12.3e %12.3e %12.3e %12.3e %12.3e\n', nRetained, bestTheta, bestRow.maxESpec, bestRow.freeMaxESpec, bestRow.ESpecF, bestRow.ESpecG, bestRow.minDeltaZ);
end
end

function printRetainedCapacitySweep(profile, referenceNodes, nEVP, latitude, g)
data = modeData(profile, "analytical", referenceNodes, nEVP, latitude, g);
nActive = size(data.PhiG, 2);
retainedModes = (1:nActive).';
freeErrors = zeros(nActive, 1);
positiveErrors = zeros(nActive, 1);
minimumPositiveWeights = zeros(nActive, 1);
for iRetained = 1:nActive
    retainedData = retainModeBand(data, retainedModes(iRetained));
    deltaZ = leastSquaresWeights(retainedData);
    freeMetrics = parsevalMetrics(retainedData, deltaZ);
    positiveMetrics = positiveParsevalMetrics(retainedData, deltaZ);
    freeErrors(iRetained) = max(freeMetrics.ESpecF, freeMetrics.ESpecG);
    positiveErrors(iRetained) = positiveMetrics.maxESpec;
    minimumPositiveWeights(iRetained) = min(positiveMetrics.deltaZ);
end

thresholds = [1e-2 1e-3 1e-4];
fprintf('\n  Retained-band capacity on canonical Bessel-root grid\n');
fprintf('  %10s %10s %10s %10s\n', "threshold", "max r", "r/N", "min dz");
for threshold = thresholds
    valid = retainedModes(positiveErrors <= threshold & isfinite(positiveErrors));
    if isempty(valid)
        fprintf('  %10.1e %10d %10.3f %10.3e\n', threshold, 0, 0, NaN);
    else
        maxRetained = max(valid);
        index = find(retainedModes == maxRetained, 1);
        fprintf('  %10.1e %10d %10.3f %10.3e\n', threshold, maxRetained, maxRetained/numel(referenceNodes), minimumPositiveWeights(index));
    end
end
halfRetained = max(1, round(numel(referenceNodes)/2) - 1);
halfIndex = find(retainedModes == halfRetained, 1);
fprintf('  half-band r = %d: positive max E = %.3e, free max E = %.3e\n', halfRetained, positiveErrors(halfIndex), freeErrors(halfIndex));
end

function z = phaseShiftedReferenceNodes(profile, nPoints, theta, latitude, g)
if abs(theta) < 100*eps
    z = referenceQuadratureNodes(profile, nPoints, latitude, g);
    return
end

zScan = linspace(-profile.depth, 0, max(5000, 100*nPoints)).';
reference = quietAnalyticalSolver(profile, zScan, nPoints, latitude, g);
[~, ~, ~, h] = evalc("reference.modesAtFrequency(0)");
[~, GRoot] = reference.ModeFunctionsForOmegaAndC(0, sqrt(g*h(nPoints - 1)));
[~, GNext] = reference.ModeFunctionsForOmegaAndC(0, sqrt(g*h(nPoints)));
rootFunction = @(z) cos(theta)*GRoot(z, 0, sqrt(g*h(nPoints-1))) + sin(theta)*GNext(z, 0, sqrt(g*h(nPoints)));
values = rootFunction(zScan);
brackets = find(values(1:end-1).*values(2:end) <= 0);
roots = zeros(numel(brackets), 1);
for iRoot = 1:numel(brackets)
    roots(iRoot) = fzero(rootFunction, zScan(brackets(iRoot):brackets(iRoot)+1));
end
roots(abs(roots + profile.depth) < 1e-9*profile.depth | abs(roots) < 1e-9*profile.depth) = [];
roots = unique(roots, "stable");
if numel(roots) < nPoints - 2
    error("HydrostaticModeAccuracy:PhaseShiftRoots", "Expected at least %d interior roots, found %d.", nPoints - 2, numel(roots));
end
if numel(roots) > nPoints - 2
    [~, order] = sort(abs(roots - mean(profile.zIn)), "ascend");
    keep = sort(order(1:(nPoints - 2)));
    roots = roots(keep);
end
z = unique([profile.zIn(1); roots; profile.zIn(2)], "stable");
if numel(z) ~= nPoints
    error("HydrostaticModeAccuracy:PhaseShiftNodes", "Expected %d nodes, found %d.", nPoints, numel(z));
end
end

function solver = quietAnalyticalSolver(profile, z, nModes, latitude, g) %#ok<INUSD>
[~, solver] = evalc("analyticalSolver(profile, z, nModes, latitude, g)");
end

function data = retainModeBand(data, nRetained)
data.PhiF = data.PhiF(:, 1:(nRetained+1));
data.PhiGFull = data.PhiGFull(:, 1:nRetained);
data.PhiG = data.PhiG(:, 1:nRetained);
data.GammaF0 = data.GammaF0(1:(nRetained+1), 1:(nRetained+1));
data.GammaG0 = eye(nRetained);
data.h = data.h(1:nRetained);
data.F = data.F(:, 1:nRetained);
data.G = data.G(:, 1:nRetained);
end

function row = diagnosticRow(label, solverKind, deltaZ, metrics, positiveMetrics, z, referenceNodes, modeErrorF, modeErrorG, hError)
row.label = string(label);
row.solverKind = string(solverKind);
row.freeMaxESpec = max(metrics.ESpecF, metrics.ESpecG);
row.freeESpecF = metrics.ESpecF;
row.freeESpecG = metrics.ESpecG;
row.maxESpec = positiveMetrics.maxESpec;
row.ESpecF = positiveMetrics.ESpecF;
row.ESpecG = positiveMetrics.ESpecG;
row.positiveMaxESpec = positiveMetrics.maxESpec;
row.maxNodeError = max(abs(z(:) - referenceNodes(:)));
row.modeErrorF = modeErrorF;
row.modeErrorG = modeErrorG;
row.hError = hError;
row.freeMinDeltaZ = min(deltaZ);
row.minDeltaZ = min(positiveMetrics.deltaZ);
row.maxDeltaZ = max(positiveMetrics.deltaZ);
row.sumDeltaZ = sum(positiveMetrics.deltaZ);
end

function value = nodeWeightObjective(variables, evaluator, profile, g, nRetained)
[z, deltaZ] = unpackNodeWeightVariables(variables, profile);
data = modeDataFromEvaluator(profile, evaluator, z, g);
data = retainModeBand(data, nRetained);
[A, b] = gramSystem(data);
residual = A*deltaZ - b;
value = residual.'*residual;
end

function [z, deltaZ] = unpackNodeWeightVariables(variables, profile)
nPoints = (numel(variables) + 2)/2;
nInterior = nPoints - 2;
z = [profile.zIn(1); variables(1:nInterior); profile.zIn(2)];
deltaZ = variables((nInterior+1):end);
end

function evaluator = exponentialModeEvaluator(profile, initialNodes, latitude, g)
nActive = numel(initialNodes) - 2;
solver = analyticalSolver(profile, initialNodes, nActive, latitude, g);
[~, ~, h] = solver.modesAtFrequency(0);
evaluator.h = reshape(h(1:nActive), [], 1);
evaluator.FFunctions = cell(nActive, 1);
evaluator.GFunctions = cell(nActive, 1);
evaluator.normalizations = zeros(nActive, 1);
for j = 1:nActive
    c = sqrt(g*evaluator.h(j));
    [FFunction, GFunction] = solver.ModeFunctionsForOmegaAndC(0, c);
    geostrophicNorm = integral(@(z) profile.N2Function(z).*GFunction(z, 0, c).^2/g, profile.zIn(1), profile.zIn(2));
    normalization = sqrt(geostrophicNorm);
    if FFunction(profile.zIn(2), 0, c) < 0
        normalization = -normalization;
    end
    evaluator.FFunctions{j} = FFunction;
    evaluator.GFunctions{j} = GFunction;
    evaluator.normalizations(j) = normalization;
end
end

function data = modeDataFromEvaluator(profile, evaluator, z, g)
nPoints = numel(z);
nActive = numel(evaluator.h);
F = zeros(nPoints, nActive);
G = zeros(nPoints, nActive);
for j = 1:nActive
    c = sqrt(g*evaluator.h(j));
    F(:, j) = evaluator.FFunctions{j}(z, 0, c)/evaluator.normalizations(j);
    G(:, j) = evaluator.GFunctions{j}(z, 0, c)/evaluator.normalizations(j);
end

data.z = z(:);
data.N2 = profile.N2Function(z(:));
data.PhiF = [ones(nPoints, 1), F];
data.PhiGFull = G;
data.PhiG = G(2:end-1, :);
data.GammaF0 = diag([profile.depth; evaluator.h]);
data.GammaG0 = eye(nActive);
data.h = evaluator.h;
data.F = F;
data.G = G;
data.depth = profile.depth;
data.g = g;
end

function row = emptyComparisonRow()
row.label = "";
row.solverKind = "";
row.nRetained = NaN;
row.maxESpec = NaN;
row.ESpecF = NaN;
row.ESpecG = NaN;
row.freeMaxESpec = NaN;
row.freeESpecF = NaN;
row.freeESpecG = NaN;
row.positiveMaxESpec = NaN;
row.maxNodeError = NaN;
row.modeErrorF = NaN;
row.modeErrorG = NaN;
row.hError = NaN;
row.freeMinDeltaZ = NaN;
row.minDeltaZ = NaN;
row.maxDeltaZ = NaN;
row.sumDeltaZ = NaN;
end

function rows = nEVPSweep(profile, nPoints, nEVPList, latitude, g)
referenceNodes = referenceQuadratureNodes(profile, nPoints, latitude, g);
referenceData = modeData(profile, "analytical", referenceNodes, max(nEVPList), latitude, g);
rows = struct([]);
for iEVP = 1:numel(nEVPList)
    nEVP = nEVPList(iEVP);
    if 2*nPoints >= nEVP
        continue
    end
    rawNodes = numericalQuadratureNodes(profile, "wkb", nPoints, nEVP, latitude, g);
    row = comparisonRow(profile, "wkb modes + raw roots", "wkb", rawNodes, referenceNodes, referenceData, nEVP, latitude, g);
    row.nEVP = nEVP;
    rows = [rows row]; %#ok<AGROW>
end
end

function data = modeData(profile, solverKind, z, nEVP, latitude, g)
nPoints = numel(z);
nActive = nPoints - 2;
nModes = nPoints - 1;
switch solverKind
    case "analytical"
        solver = analyticalSolver(profile, z, nModes, latitude, g);
    otherwise
        solver = numericalSolver(profile, solverKind, z, nEVP, nModes, latitude, g);
end
[F, G, h] = solver.modesAtFrequency(0);

data.z = z(:);
data.N2 = solver.N2(:);
data.PhiF = [ones(nPoints, 1), F(:, 1:nActive)];
data.PhiGFull = G(:, 1:nActive);
data.PhiG = G(2:end-1, 1:nActive);
data.GammaF0 = diag([profile.depth; reshape(h(1:nActive), [], 1)]);
data.GammaG0 = eye(nActive);
data.h = reshape(h(1:nActive), [], 1);
data.F = F(:, 1:nActive);
data.G = G(:, 1:nActive);
data.depth = profile.depth;
data.g = g;
end

function solver = analyticalSolver(profile, z, nModes, latitude, g)
if profile.name == "constant-limit"
    solver = InternalModesConstantStratification(N0=profile.N0, zIn=profile.zIn, zOut=z(:), latitude=latitude, nModes=nModes, g=g);
else
    solver = InternalModesExponentialStratification(N0=profile.N0, b=profile.b, zIn=profile.zIn, zOut=z(:), latitude=latitude, nModes=nModes, g=g);
end
solver.upperBoundary = UpperBoundary.rigidLid;
solver.normalization = Normalization.geostrophic;
end

function solver = numericalSolver(profile, solverKind, z, nEVP, nModes, latitude, g)
switch solverKind
    case "wkb"
        solver = InternalModesWKBSpectral(N2=profile.N2Function, zIn=profile.zIn, zOut=z(:), latitude=latitude, nEVP=nEVP, nModes=nModes, g=g);
    case "spectral"
        solver = InternalModesSpectral(N2=profile.N2Function, zIn=profile.zIn, zOut=z(:), latitude=latitude, nEVP=nEVP, nModes=nModes, g=g);
    otherwise
        error("HydrostaticModeAccuracy:UnknownSolver", "Unknown solver kind %s.", solverKind);
end
solver.upperBoundary = UpperBoundary.rigidLid;
solver.normalization = Normalization.geostrophic;
end

function deltaZ = leastSquaresWeights(data)
nF = size(data.PhiF, 2);
nG = size(data.PhiGFull, 2);
nPoints = numel(data.z);
scaleF = diag(1./sqrt(diag(data.GammaF0)));
AF = zeros(nF*nF, nPoints);
bF = zeros(nF*nF, 1);
row = 0;
for p = 1:nF
    for q = 1:nF
        row = row + 1;
        AF(row, :) = (scaleF(p, p)*data.PhiF(:, p).*data.PhiF(:, q)*scaleF(q, q)).';
        bF(row) = double(p == q);
    end
end

AG = zeros(nG*nG, nPoints);
bG = zeros(nG*nG, 1);
row = 0;
for p = 1:nG
    for q = 1:nG
        row = row + 1;
        AG(row, :) = (data.N2(:).*data.PhiGFull(:, p).*data.PhiGFull(:, q)/data.g).';
        bG(row) = double(p == q);
    end
end

A = [AF; AG];
b = [bF; bG];
if exist("lsqminnorm", "file") == 2
    deltaZ = lsqminnorm(A, b);
else
    warningState = warning("off", "MATLAB:rankDeficientMatrix");
    cleanup = onCleanup(@() warning(warningState));
    deltaZ = A\b;
end
deltaZ = real(deltaZ(:));
end

function metrics = parsevalMetrics(data, deltaZ)
GammaF = data.PhiF.'*diag(deltaZ)*data.PhiF;
scaleF = diag(1./sqrt(diag(data.GammaF0)));
GammaG = data.PhiG.'*diag(data.N2(2:end-1).*deltaZ(2:end-1)/data.g)*data.PhiG;
metrics.ESpecF = norm(scaleF*(GammaF - data.GammaF0)*scaleF, 2);
metrics.ESpecG = norm(GammaG - data.GammaG0, 2);
end

function metrics = positiveParsevalMetrics(data, unconstrainedDeltaZ)
metrics.maxESpec = NaN;
metrics.ESpecF = NaN;
metrics.ESpecG = NaN;
metrics.deltaZ = NaN(size(unconstrainedDeltaZ));
if min(unconstrainedDeltaZ) >= -sqrt(eps)*max(1, sum(abs(unconstrainedDeltaZ)))
    unconstrainedMetrics = parsevalMetrics(data, unconstrainedDeltaZ);
    metrics = metricsWithWeights(unconstrainedMetrics, unconstrainedDeltaZ);
    return
end

[A, b] = gramSystem(data);
if exist("lsqlin", "file") == 2
    options = optimoptions("lsqlin", Display="off");
    [deltaZPositive, ~, ~, exitFlag] = lsqlin(A, b, [], [], [], [], zeros(numel(data.z), 1), [], [], options);
    if exitFlag > 0
        positiveMetrics = parsevalMetrics(data, deltaZPositive);
        metrics = metricsWithWeights(positiveMetrics, deltaZPositive);
    end
elseif exist("lsqnonneg", "file") == 2
    deltaZPositive = lsqnonneg(A, b);
    positiveMetrics = parsevalMetrics(data, deltaZPositive);
    metrics = metricsWithWeights(positiveMetrics, deltaZPositive);
end
end

function metrics = metricsWithWeights(metrics, deltaZ)
metrics.maxESpec = max(metrics.ESpecF, metrics.ESpecG);
metrics.deltaZ = real(deltaZ(:));
end

function [A, b] = gramSystem(data)
nF = size(data.PhiF, 2);
nG = size(data.PhiGFull, 2);
nPoints = numel(data.z);
scaleF = diag(1./sqrt(diag(data.GammaF0)));
A = zeros(nF*nF + nG*nG, nPoints);
b = zeros(size(A, 1), 1);
row = 0;
for p = 1:nF
    for q = 1:nF
        row = row + 1;
        A(row, :) = (scaleF(p, p)*data.PhiF(:, p).*data.PhiF(:, q)*scaleF(q, q)).';
        b(row) = double(p == q);
    end
end
for p = 1:nG
    for q = 1:nG
        row = row + 1;
        A(row, :) = (data.N2(:).*data.PhiGFull(:, p).*data.PhiGFull(:, q)/data.g).';
        b(row) = double(p == q);
    end
end
end

function [modeErrorF, modeErrorG, hError] = referenceErrors(data, referenceData)
referenceF = interpolateColumns(referenceData.z, referenceData.F, data.z);
referenceG = interpolateColumns(referenceData.z, referenceData.G, data.z);
for j = 1:size(data.F, 2)
    if data.F(:, j).'*referenceF(:, j) < 0
        referenceF(:, j) = -referenceF(:, j);
        referenceG(:, j) = -referenceG(:, j);
    end
end
modeErrorF = norm(data.F - referenceF, "fro")/max(1, norm(referenceF, "fro"));
modeErrorG = norm(data.G - referenceG, "fro")/max(1, norm(referenceG, "fro"));
hError = norm(data.h - referenceData.h)/max(1, norm(referenceData.h));
end

function values = interpolateColumns(zReference, valuesReference, z)
values = zeros(numel(z), size(valuesReference, 2));
for j = 1:size(valuesReference, 2)
    values(:, j) = interp1(zReference, valuesReference(:, j), z, "pchip");
end
end

function printComparisonRows(rows)
fprintf(['  %-38s %5s %10s %10s %10s %10s %10s %10s %10s ' ...
    '%10s %10s %10s\n'], "candidate", "r", "max E", "free E", "E_F", "E_G", "node err", "mode F", "mode G", "min dz", "max dz", "sum dz");
for row = rows
    values = [row.maxESpec row.freeMaxESpec row.ESpecF row.ESpecG row.maxNodeError row.modeErrorF row.modeErrorG row.minDeltaZ row.maxDeltaZ row.sumDeltaZ];
    fprintf(['  %-38s %5d %10.3e %10.3e %10.3e %10.3e %10.3e ' ...
        '%10.3e %10.3e %10.3e %10.3e %10.3e\n'], row.label, row.nRetained, values);
end
end

function printRecommendationRows(rows)
fullRetained = rows(1).nRetained;
referenceError = rows(1).maxESpec;
rawRootPenalty = abs(rows(2).maxESpec - referenceError);
fullMask = [rows.nRetained] == fullRetained;
candidateErrors = [rows(fullMask).maxESpec];
candidateLabels = [rows(fullMask).label];
[bestError, bestProductionIndex] = min(candidateErrors);
fprintf('\n  Recommendations\n');
fprintf('    analytical raw-root penalty: %.3e\n', rawRootPenalty);
fprintf('    best full-band candidate: %s with max E = %.3e\n', candidateLabels(bestProductionIndex), bestError);
retainedMask = [rows.nRetained] < fullRetained;
if any(retainedMask)
    [bestRetainedError, bestRetainedIndex] = min([rows(retainedMask).maxESpec]);
    retainedRows = rows(retainedMask);
    fprintf('    best reduced-band candidate: %s with max E = %.3e\n', retainedRows(bestRetainedIndex).label, bestRetainedError);
end
if contains(candidateLabels(bestProductionIndex), "optimized")
    fprintf('    moving nodes helps; the remaining combined F/G Gram error is still significant.\n');
elseif referenceError > 0.1*bestError
    fprintf('    fixed-node Gram-fit error is the dominant remaining limit.\n');
elseif rawRootPenalty > 0.1*bestError
    fprintf('    root accuracy is a material part of the remaining limit.\n');
else
    fprintf('    mode/eigenvector error is the dominant remaining limit.\n');
end
end

function printSweepRows(rows)
fprintf('  %8s %10s %10s %10s %10s %10s\n', "nEVP", "max E", "E_F", "E_G", "node err", "mode G");
for row = rows
    fprintf('  %8d %10.3e %10.3e %10.3e %10.3e %10.3e\n', row.nEVP, row.maxESpec, row.ESpecF, row.ESpecG, row.maxNodeError, row.modeErrorG);
end
[~, bestIndex] = min([rows.maxESpec]);
fprintf('  best nEVP in sweep: %d\n', rows(bestIndex).nEVP);
end
