%% Scratch partial-depth quadrature investigation
% This script tests observation-window weights for partial-depth G-mode
% projections. The primary score is not full-depth Parseval accuracy, because
% a partial-depth window cannot see the whole water column. Instead, the
% first score asks how well the observation weights reproduce the physical
% G-mode inner product over the observed window. A separate Parseval score
% shows how far the same weighted Gram matrix is from a full-depth spectrum.

%% Scratch controls
nReference = 100;
rTarget = 20;
nCandidateModes = 48;
nObs = 26;
nWindowDense = 4096;
samplingScenarios = ["irregularPartialDepth","uniformUpper300m"];
weightPolicies = ["geometricWindow","aggregatedFullWindow","aggregatedTargetWindow","nonnegativeWindowGramFit"];
noRefitPolicies = ["geometricWindow","aggregatedFullWindow","aggregatedTargetWindow"];
selectorPolicies = ["normalizedQR","windowEnergyRank","windowEnergyPivot"];
rankTolerance = 5e-2;
maxConditionNumber = 10;
visibilityThreshold = 5e-2;

definition = hydrostaticScratchDefinition();
if nCandidateModes < rTarget
    error("PartialDepthQuadratureScratch:TooFewCandidateModes", ...
        "nCandidateModes must be at least rTarget.");
end

if ~hasNonnegativeLeastSquares()
    error("PartialDepthQuadratureScratch:MissingNonnegativeSolver", ...
        "This scratch script requires lsqlin or lsqnonneg so fitted quadrature weights remain nonnegative.");
end

fprintf('\nPartial-depth hydrostatic G-mode quadrature scratch study\n');
fprintf('Reference grid: nReference=%d; target prefix: rTarget=%d; observations=%d.\n', nReference, rTarget, nObs);
fprintf('Sampling scenarios: %s.\n', char(strjoin(samplingScenarios,", ")));
fprintf('Weight policies: %s.\n', char(strjoin(weightPolicies,", ")));
fprintf('Selector policies: %s.\n', char(strjoin(selectorPolicies,", ")));
fprintf('Window score compares to dense physical integration over the observed window.\n');
fprintf('Partial-window potential energy is the relevant Parseval quantity for partial-depth data.\n');
fprintf('Full-depth identity mismatch is reported only as a reminder that the window is incomplete.\n');
fprintf('Selection quality uses partial-window energy explained, conditioning, observed residual, and aliasing together.\n');

reference = referenceAtlas(definition,nReference,rTarget);
rows = struct([]);
selectionRows = struct([]);

for iScenario = 1:length(samplingScenarios)
    samplingScenario = samplingScenarios(iScenario);
    scenario = scenarioData(definition,samplingScenario,nObs,nWindowDense,nCandidateModes,rTarget);
    fprintf('\nScenario: %s (%s)\n', char(scenario.name), char(scenario.description));
    fprintf('Observed window: %.1f m to %.1f m, width %.1f m.\n', ...
        min(scenario.zObs),max(scenario.zObs),max(scenario.zObs) - min(scenario.zObs));
    fprintf('Full-depth identity mismatch floor ||Gamma_window-I||_2: %.3e.\n', scenario.parsevalFloor);

    dzGeometric = geometricIncrements(scenario.zObs);
    row = diagnosticRow(scenario,"geometricWindow",dzGeometric,rankTolerance,maxConditionNumber);
    rows = appendRow(rows,row);
    selectionRows = appendRows(selectionRows,selectorDiagnosticsForPolicy(scenario,row,selectorPolicies,visibilityThreshold));

    dzAggregatedFull = aggregateReferenceIncrementsOverWindow(reference.z,reference.dzFullActive,scenario.zObs);
    row = diagnosticRow(scenario,"aggregatedFullWindow",dzAggregatedFull,rankTolerance,maxConditionNumber);
    rows = appendRow(rows,row);
    selectionRows = appendRows(selectionRows,selectorDiagnosticsForPolicy(scenario,row,selectorPolicies,visibilityThreshold));

    dzAggregatedTarget = aggregateReferenceIncrementsOverWindow(reference.z,reference.dzTarget,scenario.zObs);
    row = diagnosticRow(scenario,"aggregatedTargetWindow",dzAggregatedTarget,rankTolerance,maxConditionNumber);
    rows = appendRow(rows,row);
    selectionRows = appendRows(selectionRows,selectorDiagnosticsForPolicy(scenario,row,selectorPolicies,visibilityThreshold));

    [dzWindowFit,exitflag] = nonnegativeWindowGramFit(scenario);
    row = diagnosticRow(scenario,"nonnegativeWindowGramFit",dzWindowFit,rankTolerance,maxConditionNumber);
    row.exitflag = exitflag;
    rows = appendRow(rows,row);
    selectionRows = appendRows(selectionRows,selectorDiagnosticsForPolicy(scenario,row,selectorPolicies,visibilityThreshold));

    printScenarioTables(rows,samplingScenario,noRefitPolicies);
    printSelectionTables(selectionRows,samplingScenario);
end

printNoRefitSummary(rows,noRefitPolicies);
printOverallSummary(rows,weightPolicies);
printSelectionSummary(selectionRows,weightPolicies);
plotPolicyDiagnostics(rows,samplingScenarios,weightPolicies);
plotWindowWeights(rows,samplingScenarios,weightPolicies);
plotBestHeatmaps(rows,samplingScenarios,noRefitPolicies);

%% Local helpers
function definition = hydrostaticScratchDefinition()
definition.depth = 4000;
definition.N0 = 3*2*pi/3600;
definition.L_gm = 1300;
definition.N2 = @(z) definition.N0*definition.N0*exp(2*z/definition.L_gm);
definition.latitude = 31;
definition.g = 9.81;
definition.omega = 0;
definition.zDomain = [-definition.depth 0];
end

function reference = referenceAtlas(definition,nReference,rTarget)
z = gRootGrid(definition,nReference);
fullActiveData = fullDepthModeData(definition,z,nReference - 2);
targetData = fullDepthModeData(definition,z,rTarget);

[dzFullActive,fullExitflag] = positiveJointLeastSquaresIncrements(fullActiveData);
[dzTarget,targetExitflag] = positiveJointLeastSquaresIncrements(targetData);

reference.z = z;
reference.dzFullActive = dzFullActive;
reference.dzTarget = dzTarget;
reference.fullExitflag = fullExitflag;
reference.targetExitflag = targetExitflag;
end

function scenario = scenarioData(definition,samplingScenario,nObs,nWindowDense,nCandidateModes,rTarget)
[zObs,description] = observationGridForScenario(samplingScenario,nObs);
zWindow = linspace(min(zObs),max(zObs),nWindowDense).';
dzWindow = geometricIncrements(zWindow);

GObs = gModesAtDepths(definition,zObs,nCandidateModes);
GWindow = gModesAtDepths(definition,zWindow,nCandidateModes);
windowWeights = definition.N2(zWindow).*dzWindow/definition.g;
gammaWindowReference = GWindow(:,1:rTarget).'*(windowWeights.*GWindow(:,1:rTarget));
candidateGammaWindowReference = GWindow.'*(windowWeights.*GWindow);

scenario.name = string(samplingScenario);
scenario.description = description;
scenario.zObs = zObs;
scenario.N2Obs = definition.N2(zObs);
scenario.g = definition.g;
scenario.GObs = GObs;
scenario.GWindow = GWindow;
scenario.windowWeights = windowWeights;
scenario.rTarget = rTarget;
scenario.nCandidateModes = nCandidateModes;
scenario.gammaWindowReference = symmetrize(gammaWindowReference);
scenario.candidateGammaWindowReference = symmetrize(candidateGammaWindowReference);
scenario.parsevalFloor = norm(scenario.gammaWindowReference - eye(rTarget),2);
end

function [zObs,description] = observationGridForScenario(samplingScenario,nObs)
switch string(samplingScenario)
    case "irregularPartialDepth"
        zObsUniform = linspace(-1600,-80,nObs).';
        zObs = zObsUniform + 35*sin((1:length(zObsUniform)).');
        zObs = sort(zObs);
        description = string(sprintf('irregular samples from %.0f m to %.0f m', min(zObs), max(zObs)));
    case "uniformUpper300m"
        zObs = linspace(-380,-80,nObs).';
        description = "uniform 300 m window from -380 m to -80 m";
    otherwise
        error("PartialDepthQuadratureScratch:UnknownSamplingScenario", ...
            "Unknown samplingScenario ""%s"".", char(string(samplingScenario)));
end
end

function z = gRootGrid(definition,nPoints)
nEVP = max(256,ceil(2.1*(nPoints + 1)));
zReference = linspace(definition.zDomain(1),definition.zDomain(2),1024).';
imReference = InternalModesWKBSpectral(N2=definition.N2,zIn=definition.zDomain,zOut=zReference,latitude=definition.latitude,nEVP=nEVP,g=definition.g);
imReference.normalization = Normalization.geostrophic;
imReference.upperBoundary = UpperBoundary.rigidLid;
z = imReference.GaussQuadraturePointsForModesAtFrequency(nPoints,definition.omega);
end

function G = gModesAtDepths(definition,z,nModes)
nEVP = max(256,ceil(2.1*(nModes + 1)));
im = InternalModesWKBSpectral(N2=definition.N2,zIn=definition.zDomain,zOut=z,latitude=definition.latitude,nEVP=nEVP,nModes=nModes,g=definition.g);
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;
[~,G] = im.modesAtFrequency(definition.omega);
G = G(:,1:nModes);
end

function data = fullDepthModeData(definition,z,r)
nEVP = max(256,ceil(2.1*(max(length(z),r + 2) + 1)));
im = InternalModesWKBSpectral(N2=definition.N2,zIn=definition.zDomain,zOut=z,latitude=definition.latitude,nEVP=nEVP,nModes=r,g=definition.g);
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;
[F,G,h] = im.modesAtFrequency(definition.omega);
h = h(:);

data.z = z(:);
data.nPoints = length(z);
data.r = r;
data.depth = definition.depth;
data.g = definition.g;
data.N2 = definition.N2(data.z);
data.N2Interior = data.N2(2:(end - 1));
data.PhiF = cat(2,ones(data.nPoints,1),F(:,1:r));
data.PhiG = G(2:(end - 1),1:r);
data.GammaF0 = diag([definition.depth; h(1:r)]);
data.GammaG0 = eye(r);
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2) - z(1))/2,(z(3:end) - z(1:(end - 2)))/2,(z(end) - z(end - 1))/2);
end

function dz = aggregateReferenceIncrementsOverWindow(zReference,dzReference,zObs)
isInWindow = zReference >= min(zObs) & zReference <= max(zObs);
zReferenceWindow = zReference(isInWindow);
dzReferenceWindow = dzReference(isInWindow);
[~,assignment] = min(abs(zReferenceWindow - zObs.'),[],2);
dz = accumarray(assignment,dzReferenceWindow,[length(zObs) 1],@sum,0);
end

function [dz,exitflag] = nonnegativeWindowGramFit(scenario)
[A,rhs] = windowGramFitSystem(scenario.GObs(:,1:scenario.rTarget),scenario.N2Obs/scenario.g,scenario.gammaWindowReference);
[dz,exitflag] = positiveLeastSquares(A,rhs,length(scenario.zObs));
end

function [A,rhs] = windowGramFitSystem(G,weightFactor,gammaTarget)
nObs = size(G,1);
r = size(G,2);
A = zeros(r*r,nObs);
rhs = zeros(r*r,1);
row = 0;
for iMode = 1:r
    for jMode = 1:r
        row = row + 1;
        A(row,:) = (weightFactor.*G(:,iMode).*G(:,jMode)).';
        rhs(row) = gammaTarget(iMode,jMode);
    end
end
end

function [dz,exitflag] = positiveJointLeastSquaresIncrements(data)
[A,rhs] = jointLeastSquaresSystem(data);
[dz,exitflag] = positiveLeastSquares(A,rhs,data.nPoints);
end

function [A,rhs] = jointLeastSquaresSystem(data)
nF = size(data.PhiF,2);
nG = size(data.PhiG,2);
nRows = nF*nF + nG*nG;
A = zeros(nRows,data.nPoints);
rhs = zeros(nRows,1);
row = 0;

for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        scale = sqrt(data.GammaF0(iMode,iMode)*data.GammaF0(jMode,jMode));
        A(row,:) = (data.PhiF(:,iMode).*data.PhiF(:,jMode)).'/scale;
        rhs(row) = double(iMode == jMode);
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(data.nPoints,1);
        weights(2:(end - 1)) = (data.N2Interior/data.g).*data.PhiG(:,iMode).*data.PhiG(:,jMode);
        A(row,:) = weights.';
        rhs(row) = double(iMode == jMode);
    end
end
end

function [dz,exitflag] = positiveLeastSquares(A,rhs,nWeights)
dz = [];
exitflag = NaN;
try
    if exist("lsqlin","file") == 2
        options = optimoptions("lsqlin",Display="off");
        [dz,~,~,exitflag] = lsqlin(A,rhs,[],[],[],[],zeros(nWeights,1),[],[],options);
    elseif exist("lsqnonneg","file") == 2
        [dz,~,~,exitflag] = lsqnonneg(A,rhs);
    end
catch
    dz = [];
    exitflag = NaN;
end
end

function value = hasNonnegativeLeastSquares()
value = exist("lsqlin","file") == 2 || exist("lsqnonneg","file") == 2;
end

function row = diagnosticRow(scenario,weightPolicy,dz,rankTolerance,maxConditionNumber)
GPrefix = scenario.GObs(:,1:scenario.rTarget);
weightsG = scenario.N2Obs.*dz/scenario.g;
gammaObs = GPrefix.'*(weightsG.*GPrefix);
gammaObs = symmetrize(gammaObs);

[retainedModes,selectionDiagnostics] = selectResolvableModes(scenario.GObs,weightsG,rankTolerance,maxConditionNumber);
projection = projectionDiagnostics(scenario.GObs,weightsG,retainedModes);
if isempty(retainedModes)
    qrWindowError = NaN;
    qrParsevalError = NaN;
else
    gammaQR = projection.gramMatrix;
    gammaWindowQR = scenario.candidateGammaWindowReference(retainedModes,retainedModes);
    qrWindowError = norm(symmetrize(gammaQR - gammaWindowQR),2);
    qrParsevalError = norm(symmetrize(gammaQR - eye(length(retainedModes))),2);
end

row.scenario = scenario.name;
row.description = scenario.description;
row.weightPolicy = string(weightPolicy);
row.zObs = scenario.zObs;
row.dz = dz;
row.weightsG = weightsG;
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.negativeDzCount = sum(dz < -1e-10);
row.gammaObs = gammaObs;
row.gammaWindowReference = scenario.gammaWindowReference;
row.windowError = norm(symmetrize(gammaObs - scenario.gammaWindowReference),2);
row.parsevalError = norm(symmetrize(gammaObs - eye(scenario.rTarget)),2);
row.parsevalFloor = scenario.parsevalFloor;
row.frobeniusWindowError = norm(gammaObs - scenario.gammaWindowReference,"fro");
row.maxDiagWindowError = max(abs(diag(gammaObs - scenario.gammaWindowReference)));
row.qrRetainedModes = retainedModes;
row.qrRetainedCount = length(retainedModes);
row.qrFinalRelativePivot = selectionDiagnostics.finalRelativePivot;
row.qrConditionNumber = selectionDiagnostics.conditionNumber;
row.qrObservationResidualFraction = projection.observationResidualFraction;
row.qrWindowError = qrWindowError;
row.qrParsevalError = qrParsevalError;
row.exitflag = NaN;
end

function rows = selectorDiagnosticsForPolicy(scenario,weightRow,selectorPolicies,visibilityThreshold)
retainedCount = weightRow.qrRetainedCount;
rows = struct([]);
if retainedCount == 0
    return
end

for iPolicy = 1:length(selectorPolicies)
    selectorPolicy = selectorPolicies(iPolicy);
    retainedModes = retainedModesForSelector(scenario,weightRow,selectorPolicy,retainedCount);
    rows = appendRow(rows,selectionDiagnosticRow(scenario,weightRow,selectorPolicy,retainedModes,visibilityThreshold));
end
end

function retainedModes = retainedModesForSelector(scenario,weightRow,selectorPolicy,retainedCount)
switch string(selectorPolicy)
    case "normalizedQR"
        retainedModes = weightRow.qrRetainedModes;
    case "windowEnergyRank"
        windowEnergy = diag(scenario.candidateGammaWindowReference);
        [~,order] = sort(windowEnergy,"descend");
        retainedModes = sort(order(1:retainedCount)).';
    case "windowEnergyPivot"
        weightedWindowModes = sqrt(scenario.windowWeights).*scenario.GWindow;
        [~,~,pivotOrder] = qr(weightedWindowModes,0);
        retainedModes = sort(pivotOrder(1:retainedCount)).';
    otherwise
        error("PartialDepthQuadratureScratch:UnknownSelectorPolicy", ...
            "Unknown selector policy ""%s"".", char(string(selectorPolicy)));
end
end

function row = selectionDiagnosticRow(scenario,weightRow,selectorPolicy,retainedModes,visibilityThreshold)
projection = projectionDiagnostics(scenario.GObs,weightRow.weightsG,retainedModes);
gammaWindowSelected = scenario.candidateGammaWindowReference(retainedModes,retainedModes);
windowEigenvalues = eig(symmetrize(gammaWindowSelected));
windowEigenvalues = sort(real(windowEigenvalues));
partialResidualFraction = denseWindowResidualFraction(scenario,retainedModes);

row.scenario = scenario.name;
row.weightPolicy = weightRow.weightPolicy;
row.selectorPolicy = string(selectorPolicy);
row.retainedModes = retainedModes;
row.retainedCount = length(retainedModes);
row.partialEnergyExplainedFraction = 1 - partialResidualFraction;
row.partialEnergyResidualFraction = partialResidualFraction;
row.meanVisibleEnergy = trace(gammaWindowSelected)/length(retainedModes);
row.minVisibleEigenvalue = windowEigenvalues(1);
row.medianVisibleEigenvalue = median(windowEigenvalues);
row.maxVisibleEigenvalue = windowEigenvalues(end);
row.visibleRank = sum(windowEigenvalues > visibilityThreshold);
row.conditionNumber = cond(projection.gramMatrix);
row.observationResidualFraction = projection.observationResidualFraction;
row.maxAliasing = projection.maxAliasing;
end

function residualFraction = denseWindowResidualFraction(scenario,retainedModes)
GSelected = scenario.GWindow(:,retainedModes);
gramMatrix = GSelected.'*(scenario.windowWeights.*GSelected);
A = gramMatrix\(scenario.windowWeights.*GSelected).';
projectedModes = GSelected*A*scenario.GWindow;
residualModes = scenario.GWindow - projectedModes;
residualFraction = weightedObservationVariance(residualModes,scenario.windowWeights)/weightedObservationVariance(scenario.GWindow,scenario.windowWeights);
end

function [retainedModes,diagnostics] = selectResolvableModes(B,weights,rankTolerance,maxConditionNumber)
weightedB = sqrt(weights) .* B;
columnNorms = vecnorm(weightedB,2,1);
validColumns = columnNorms > 0;
normalizedB = zeros(size(weightedB));
normalizedB(:,validColumns) = weightedB(:,validColumns) ./ columnNorms(validColumns);

if ~any(validColumns)
    retainedModes = [];
    diagnostics.finalRelativePivot = NaN;
    diagnostics.conditionNumber = NaN;
    return
end

validModeNumbers = find(validColumns);
[~,R,pivotOrderValid] = qr(normalizedB(:,validColumns),0);
pivotOrder = validModeNumbers(pivotOrderValid);
relativePivots = abs(diag(R))/abs(R(1,1));
nRank = find(relativePivots >= rankTolerance,1,'last');
if isempty(nRank)
    nRank = 0;
end

selectedInPivotOrder = zeros(1,nRank);
nSelected = 0;
finalPivotIndex = NaN;
for iPivot = 1:nRank
    candidateModes = sort([selectedInPivotOrder(1:nSelected) pivotOrder(iPivot)]);
    candidateConditionNumber = cond(normalizedB(:,candidateModes).'*normalizedB(:,candidateModes));
    if candidateConditionNumber <= maxConditionNumber
        nSelected = nSelected + 1;
        selectedInPivotOrder(nSelected) = pivotOrder(iPivot);
        finalPivotIndex = iPivot;
    end
end

retainedModes = sort(selectedInPivotOrder(1:nSelected));
if isempty(retainedModes)
    diagnostics.finalRelativePivot = NaN;
    diagnostics.conditionNumber = NaN;
else
    diagnostics.finalRelativePivot = relativePivots(finalPivotIndex);
    diagnostics.conditionNumber = cond(normalizedB(:,retainedModes).'*normalizedB(:,retainedModes));
end
end

function diagnostics = projectionDiagnostics(B,weights,retainedModes)
nCandidateModes = size(B,2);
if isempty(retainedModes)
    diagnostics.gramMatrix = [];
    diagnostics.resolutionMatrix = [];
    diagnostics.observationResidualFraction = NaN;
    diagnostics.maxAliasing = NaN;
    return
end

[A,~,gramMatrix] = weightedProjection(B,weights,retainedModes);
observationProjectionMatrix = B(:,retainedModes)*A;
observationResidualMatrix = B - observationProjectionMatrix*B;
diagnostics.gramMatrix = gramMatrix;
diagnostics.resolutionMatrix = A*B;
diagnostics.observationResidualFraction = weightedObservationVariance(observationResidualMatrix,weights)/weightedObservationVariance(B,weights);
diagnostics.rejectedModes = setdiff(1:nCandidateModes,retainedModes);
diagnostics.aliasingMatrix = diagnostics.resolutionMatrix(:,diagnostics.rejectedModes);
diagnostics.aliasingColumnNorm = vecnorm(diagnostics.aliasingMatrix,2,1);
diagnostics.maxAliasing = max([0 diagnostics.aliasingColumnNorm]);
end

function [A,resolutionMatrix,gramMatrix] = weightedProjection(B,weights,retainedModes)
BS = B(:,retainedModes);
weightedBS = weights .* BS;
gramMatrix = BS.'*weightedBS;
A = gramMatrix\weightedBS.';
resolutionMatrix = A*B;
end

function variance = weightedObservationVariance(values,weights)
variance = sum(weights .* sum(values.^2,2));
end

function rows = appendRow(rows,row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row;
end
end

function rows = appendRows(rows,newRows)
if isempty(newRows)
    return
end

if isempty(rows)
    rows = newRows;
else
    rows = [rows newRows];
end
end

function value = symmetrize(value)
value = (value + value.')/2;
end

function printScenarioTables(rows,samplingScenario,noRefitPolicies)
scenarioRows = rows([rows.scenario] == samplingScenario);
[~,order] = sort([scenarioRows.windowError]);
scenarioRows = scenarioRows(order);

fprintf('\nPartial-window Parseval quadrature diagnostics\n');
fprintf('%-28s %11s %11s %11s %11s %11s %7s\n', ...
    'weight policy','E_partial','E_full','full floor','sum dz','min dz','neg dz');
for iRow = 1:length(scenarioRows)
    row = scenarioRows(iRow);
    fprintf('%-28s %11.3e %11.3e %11.3e %11.3e %11.3e %7d\n', ...
        char(row.weightPolicy),row.windowError,row.parsevalError,row.parsevalFloor, ...
        row.sumDz,row.minDz,row.negativeDzCount);
end

fprintf('\nObservational QR diagnostics\n');
fprintf('%-28s %8s %11s %11s %11s %11s\n', ...
    'weight policy','retained','cond','obs resid','QR E_win','QR E_pars');
for iRow = 1:length(scenarioRows)
    row = scenarioRows(iRow);
    fprintf('%-28s %8d %11.3e %11.3e %11.3e %11.3e\n', ...
        char(row.weightPolicy),row.qrRetainedCount,row.qrConditionNumber, ...
        row.qrObservationResidualFraction,row.qrWindowError,row.qrParsevalError);
end

noRefitRows = scenarioRows(ismember([scenarioRows.weightPolicy],noRefitPolicies));
[~,bestNoRefitIndex] = min([noRefitRows.windowError]);
bestNoRefit = noRefitRows(bestNoRefitIndex);
[~,bestOverallIndex] = min([scenarioRows.windowError]);
bestOverall = scenarioRows(bestOverallIndex);
fprintf('\nBest no-refit initial guess: %s (E_partial=%.3e, E_full=%.3e).\n', ...
    char(bestNoRefit.weightPolicy),bestNoRefit.windowError,bestNoRefit.parsevalError);
fprintf('Best overall policy: %s (E_partial=%.3e, E_full=%.3e).\n', ...
    char(bestOverall.weightPolicy),bestOverall.windowError,bestOverall.parsevalError);
end

function printSelectionTables(selectionRows,samplingScenario)
scenarioRows = selectionRows([selectionRows.scenario] == samplingScenario);
[~,order] = sort([scenarioRows.partialEnergyExplainedFraction],"descend");
scenarioRows = scenarioRows(order);

fprintf('\nRetained-mode selector diagnostics\n');
fprintf('%-24s %-18s %8s %11s %11s %8s %11s %11s %11s\n', ...
    'weight policy','selector','retained','PE explained','mean vis','vis rank','cond','obs resid','max alias');
for iRow = 1:length(scenarioRows)
    row = scenarioRows(iRow);
    fprintf('%-24s %-18s %8d %11.3e %11.3e %8d %11.3e %11.3e %11.3e\n', ...
        char(row.weightPolicy),char(row.selectorPolicy),row.retainedCount, ...
        row.partialEnergyExplainedFraction,row.meanVisibleEnergy,row.visibleRank, ...
        row.conditionNumber,row.observationResidualFraction,row.maxAliasing);
end
end

function printSelectionSummary(selectionRows,weightPolicies)
scenarios = unique([selectionRows.scenario],"stable");
fprintf('\nBest selector by partial-window energy explained\n');
fprintf('%-24s %-24s %-18s %8s %11s %11s %11s\n', ...
    'scenario','weight policy','selector','retained','PE explained','obs resid','max alias');
for iScenario = 1:length(scenarios)
    for iWeight = 1:length(weightPolicies)
        rows = selectionRows([selectionRows.scenario] == scenarios(iScenario) & [selectionRows.weightPolicy] == weightPolicies(iWeight));
        if isempty(rows)
            continue
        end

        [~,bestIndex] = max([rows.partialEnergyExplainedFraction]);
        row = rows(bestIndex);
        fprintf('%-24s %-24s %-18s %8d %11.3e %11.3e %11.3e\n', ...
            char(row.scenario),char(row.weightPolicy),char(row.selectorPolicy), ...
            row.retainedCount,row.partialEnergyExplainedFraction,row.observationResidualFraction,row.maxAliasing);
    end
end
end

function printNoRefitSummary(rows,noRefitPolicies)
scenarios = unique([rows.scenario],"stable");
fprintf('\nBest no-refit initial guess by scenario\n');
fprintf('%-24s %-28s %11s %11s %11s %11s\n', ...
    'scenario','weight policy','E_partial','E_full','full floor','sum dz');
for iScenario = 1:length(scenarios)
    scenarioRows = rows([rows.scenario] == scenarios(iScenario) & ismember([rows.weightPolicy],noRefitPolicies));
    [~,bestIndex] = min([scenarioRows.windowError]);
    row = scenarioRows(bestIndex);
    fprintf('%-24s %-28s %11.3e %11.3e %11.3e %11.3e\n', ...
        char(row.scenario),char(row.weightPolicy),row.windowError,row.parsevalError,row.parsevalFloor,row.sumDz);
end
end

function printOverallSummary(rows,weightPolicies)
scenarios = unique([rows.scenario],"stable");
fprintf('\nBest overall policy by scenario\n');
fprintf('%-24s %-28s %11s %11s %8s %11s\n', ...
    'scenario','weight policy','E_partial','E_full','retained','obs resid');
for iScenario = 1:length(scenarios)
    scenarioRows = rows([rows.scenario] == scenarios(iScenario) & ismember([rows.weightPolicy],weightPolicies));
    [~,bestIndex] = min([scenarioRows.windowError]);
    row = scenarioRows(bestIndex);
    fprintf('%-24s %-28s %11.3e %11.3e %8d %11.3e\n', ...
        char(row.scenario),char(row.weightPolicy),row.windowError,row.parsevalError, ...
        row.qrRetainedCount,row.qrObservationResidualFraction);
end
end

function plotPolicyDiagnostics(rows,samplingScenarios,weightPolicies)
windowErrors = metricMatrix(rows,samplingScenarios,weightPolicies,"windowError");
parsevalErrors = metricMatrix(rows,samplingScenarios,weightPolicies,"parsevalError");

figure('Name','Partial-depth quadrature policy diagnostics')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
bar(windowErrors)
set(gca,'XTickLabel',cellstr(samplingScenarios))
ylabel('E_{partial}')
title('Partial-window Parseval quadrature error')
legend(cellstr(weightPolicies),'Location','northwest','Interpreter','none')
grid on

nexttile
bar(parsevalErrors)
set(gca,'XTickLabel',cellstr(samplingScenarios))
ylabel('E_{full}')
title('Full-depth identity mismatch')
legend(cellstr(weightPolicies),'Location','northwest','Interpreter','none')
grid on
end

function values = metricMatrix(rows,samplingScenarios,weightPolicies,fieldName)
values = NaN(length(samplingScenarios),length(weightPolicies));
for iScenario = 1:length(samplingScenarios)
    for iPolicy = 1:length(weightPolicies)
        index = find([rows.scenario] == samplingScenarios(iScenario) & [rows.weightPolicy] == weightPolicies(iPolicy),1);
        if ~isempty(index)
            values(iScenario,iPolicy) = rows(index).(char(fieldName));
        end
    end
end
end

function plotWindowWeights(rows,samplingScenarios,weightPolicies)
figure('Name','Partial-depth observation-window increments')
tiledlayout(1,length(samplingScenarios),'TileSpacing','compact','Padding','compact')

for iScenario = 1:length(samplingScenarios)
    nexttile
    hold on
    for iPolicy = 1:length(weightPolicies)
        index = find([rows.scenario] == samplingScenarios(iScenario) & [rows.weightPolicy] == weightPolicies(iPolicy),1);
        row = rows(index);
        plot(row.dz,row.zObs,'-o','DisplayName',char(row.weightPolicy));
    end
    xlabel('dz (m)')
    ylabel('z (m)')
    title(char(samplingScenarios(iScenario)),'Interpreter','none')
    legend('Location','best','Interpreter','none')
    grid on
end
end

function plotBestHeatmaps(rows,samplingScenarios,noRefitPolicies)
figure('Name','Partial-depth prefix Gram-error heatmaps')
tiledlayout(length(samplingScenarios),2,'TileSpacing','compact','Padding','compact')

for iScenario = 1:length(samplingScenarios)
    scenarioRows = rows([rows.scenario] == samplingScenarios(iScenario));
    noRefitRows = scenarioRows(ismember([scenarioRows.weightPolicy],noRefitPolicies));
    [~,bestNoRefitIndex] = min([noRefitRows.windowError]);
    bestNoRefit = noRefitRows(bestNoRefitIndex);
    refitRow = scenarioRows([scenarioRows.weightPolicy] == "nonnegativeWindowGramFit");

    plotHeatmap(bestNoRefit.gammaObs - bestNoRefit.gammaWindowReference, ...
        samplingScenarios(iScenario) + " no-refit: " + bestNoRefit.weightPolicy)
    plotHeatmap(refitRow.gammaObs - refitRow.gammaWindowReference, ...
        samplingScenarios(iScenario) + " refit: " + refitRow.weightPolicy)
end
end

function plotHeatmap(matrix,titleText)
nexttile
imagesc(matrix)
axis image
maxValue = max(abs(matrix(:)));
if maxValue > 0
    caxis(maxValue*[-1 1])
end
colorbar
xlabel('mode index')
ylabel('mode index')
title(titleText,'Interpreter','none')
end
