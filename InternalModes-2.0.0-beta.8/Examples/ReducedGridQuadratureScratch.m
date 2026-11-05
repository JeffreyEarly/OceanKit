%% Scratch reduced-grid quadrature investigation
% This script asks how a high-resolution full-depth quadrature rule should be
% used, or not used, when building a reduced grid for a retained hydrostatic
% mode band. The 100-point rule is treated as a reference atlas. Aggregating
% its weights onto fewer nodes is one candidate strategy, not the assumed
% answer.

%% Scratch controls
nReference = 100;
rTarget = 20;
nReducedList = [22 26 32 48];
nodePolicies = ["reducedGRoot","uniform","chebyshevLobatto","referenceUniformSubset","referenceQRSubset"];
shouldIncludeUnconstrainedDiagnostics = false;

definition = hydrostaticScratchDefinition();
if any(nReducedList < rTarget + 2)
    error("ReducedGridQuadratureScratch:TooFewReducedPoints", ...
        "Each reduced grid needs at least rTarget+2 points for rigid-lid G-mode Parseval diagnostics.");
end

fprintf('\nReduced-grid hydrostatic quadrature scratch study\n');
fprintf('Reference grid: nReference=%d; retained target: rTarget=%d.\n', nReference, rTarget);
fprintf('Reduced grids: %s.\n', mat2str(nReducedList));
fprintf('Node policies: %s.\n', char(strjoin(nodePolicies,", ")));

hasNonnegativeFit = hasNonnegativeLeastSquares();
if hasNonnegativeFit
    fprintf('Nonnegative Gram-fit weights enabled.\n');
else
    error("ReducedGridQuadratureScratch:MissingNonnegativeSolver", ...
        "This scratch script requires lsqlin or lsqnonneg so fitted quadrature weights remain nonnegative.");
end
if shouldIncludeUnconstrainedDiagnostics
    fprintf('Unconstrained algebraic Gram-fit diagnostics enabled.\n');
end

reference = referenceAtlas(definition,nReference,rTarget);
printReferenceSummary(reference);

rows = struct([]);
for iReduced = 1:length(nReducedList)
    nReduced = nReducedList(iReduced);
    for iPolicy = 1:length(nodePolicies)
        nodePolicy = nodePolicies(iPolicy);
        zReduced = reducedNodesForPolicy(nodePolicy,definition,nReduced,reference);
        data = modeData(definition,zReduced,rTarget);

        dzGeometric = geometricIncrements(zReduced);
        rows = appendRow(rows,diagnosticRow(data,nodePolicy,"geometric",dzGeometric));

        dzAggregatedFull = aggregateReferenceIncrements(reference.z,reference.dzFullActive,zReduced);
        rows = appendRow(rows,diagnosticRow(data,nodePolicy,"aggregatedFull",dzAggregatedFull));

        dzAggregatedTarget = aggregateReferenceIncrements(reference.z,reference.dzTarget,zReduced);
        rows = appendRow(rows,diagnosticRow(data,nodePolicy,"aggregatedTarget",dzAggregatedTarget));

        if shouldIncludeUnconstrainedDiagnostics
            dzUnconstrainedGramFit = jointLeastSquaresIncrements(data);
            rows = appendRow(rows,diagnosticRow(data,nodePolicy,"unconstrainedGramFit",dzUnconstrainedGramFit));
        end

        [dzNonnegativeGramFit,exitflag] = positiveJointLeastSquaresIncrements(data);
        if ~isempty(dzNonnegativeGramFit)
            row = diagnosticRow(data,nodePolicy,"nonnegativeGramFit",dzNonnegativeGramFit);
            row.exitflag = exitflag;
            rows = appendRow(rows,row);
        end
    end

    printRankedRows(rows,nReduced);
end

printAggregationComparison(rows);
printWeightPolicySummary(rows);
plotSummaryRows(rows,nReducedList);
plotReducedWeights(rows,[22 32]);
plotBestHeatmaps(rows,22);

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

fullActiveData = modeData(definition,z,nReference - 2);
targetData = modeData(definition,z,rTarget);

[dzFullActive,fullActiveExitflag] = positiveJointLeastSquaresIncrements(fullActiveData);
[dzTarget,targetExitflag] = positiveJointLeastSquaresIncrements(targetData);

reference.z = z;
reference.fullActiveData = fullActiveData;
reference.targetData = targetData;
reference.dzFullActive = dzFullActive;
reference.dzTarget = dzTarget;
reference.fullActiveDiagnostics = diagnosticRow(fullActiveData,"referenceGRoot","fullActiveNonnegativeGramFit",dzFullActive);
reference.fullActiveDiagnostics.exitflag = fullActiveExitflag;
reference.fullOnTargetDiagnostics = diagnosticRow(targetData,"referenceGRoot","fullActiveNonnegativeGramFit",dzFullActive);
reference.fullOnTargetDiagnostics.exitflag = fullActiveExitflag;
reference.targetDiagnostics = diagnosticRow(targetData,"referenceGRoot","targetNonnegativeGramFit",dzTarget);
reference.targetDiagnostics.exitflag = targetExitflag;
end

function z = reducedNodesForPolicy(nodePolicy,definition,nReduced,reference)
switch string(nodePolicy)
    case "reducedGRoot"
        z = gRootGrid(definition,nReduced);
    case "uniform"
        z = linspace(definition.zDomain(1),definition.zDomain(2),nReduced).';
    case "chebyshevLobatto"
        theta = pi*(0:(nReduced - 1)).'/(nReduced - 1);
        z = -definition.depth*(1 + cos(theta))/2;
    case "referenceUniformSubset"
        indices = round(linspace(1,length(reference.z),nReduced)).';
        z = reference.z(indices);
    case "referenceQRSubset"
        z = qrSubsetReferenceNodes(reference,nReduced);
    otherwise
        error("ReducedGridQuadratureScratch:UnknownNodePolicy", ...
            "Unknown node policy ""%s"".", char(string(nodePolicy)));
end
z = reshape(z,[],1);
end

function z = gRootGrid(definition,nPoints)
nEVP = max(256,ceil(2.1*(nPoints + 1)));
zReference = linspace(definition.zDomain(1),definition.zDomain(2),1024).';
imReference = InternalModesWKBSpectral(N2=definition.N2,zIn=definition.zDomain,zOut=zReference,latitude=definition.latitude,nEVP=nEVP,g=definition.g);
imReference.normalization = Normalization.geostrophic;
imReference.upperBoundary = UpperBoundary.rigidLid;
z = imReference.GaussQuadraturePointsForModesAtFrequency(nPoints,definition.omega);
end

function z = qrSubsetReferenceNodes(reference,nReduced)
interiorCount = nReduced - 2;
[A,~] = jointLeastSquaresSystem(reference.targetData);
AInterior = A(:,2:(end - 1));
columnNorms = vecnorm(AInterior,2,1);
validColumns = columnNorms > 0;
ANormalized = AInterior(:,validColumns) ./ columnNorms(validColumns);
[~,~,pivotOrder] = qr(ANormalized,0);
validIndices = find(validColumns);
selectedInterior = sort(validIndices(pivotOrder(1:interiorCount)) + 1).';
indices = [1; selectedInterior; length(reference.z)];
z = reference.z(indices);
end

function data = modeData(definition,z,r)
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

function dz = aggregateReferenceIncrements(zReference,dzReference,zReduced)
[~,assignment] = min(abs(zReference - zReduced.'),[],2);
dz = accumarray(assignment,dzReference,[length(zReduced) 1],@sum,0);
end

function dz = jointLeastSquaresIncrements(data)
[A,rhs] = jointLeastSquaresSystem(data);
dz = leastSquaresSolution(A,rhs);
end

function [dz,exitflag] = positiveJointLeastSquaresIncrements(data)
dz = [];
exitflag = NaN;

[A,rhs] = jointLeastSquaresSystem(data);
try
    if exist("lsqlin","file") == 2
        options = optimoptions("lsqlin",Display="off");
        [dz,~,~,exitflag] = lsqlin(A,rhs,[],[],[],[],zeros(data.nPoints,1),[],[],options);
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

function row = diagnosticRow(data,nodePolicy,weightPolicy,dz)
row.nReduced = data.nPoints;
row.nodePolicy = string(nodePolicy);
row.weightPolicy = string(weightPolicy);
row.comboLabel = row.nodePolicy + " / " + row.weightPolicy;
row.z = data.z;
row.dz = dz;
row.exitflag = NaN;
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.negativeDzCount = sum(dz < -1e-10);

weightsF = dz;
weightsG = dz(2:(end - 1)).*data.N2Interior/data.g;
f = componentDiagnostics(data.PhiF,weightsF,data.GammaF0);
g = componentDiagnostics(data.PhiG,weightsG,data.GammaG0);

row.fSpecError = f.specError;
row.gSpecError = g.specError;
row.maxSpecError = max(row.fSpecError,row.gSpecError);
row.fFrobeniusError = f.frobeniusError;
row.gFrobeniusError = g.frobeniusError;
row.fMaxDiagError = f.maxDiagError;
row.gMaxDiagError = g.maxDiagError;
row.fOffdiagSpecError = f.offdiagSpecError;
row.gOffdiagSpecError = g.offdiagSpecError;
row.fConditionNumber = f.conditionNumber;
row.gConditionNumber = g.conditionNumber;
row.fNormalizedDifference = f.normalizedDifference;
row.gNormalizedDifference = g.normalizedDifference;
end

function diagnostics = componentDiagnostics(Phi,weights,Gamma0)
Gamma = Phi.'*(weights.*Phi);
normalizedDifference = normalizedMatrix(Gamma - Gamma0,Gamma0);
normalizedGram = normalizedMatrix(Gamma,Gamma0);
offdiag = normalizedDifference - diag(diag(normalizedDifference));

diagnostics.specError = norm(symmetrize(normalizedDifference),2);
diagnostics.frobeniusError = norm(normalizedDifference,"fro");
diagnostics.maxDiagError = max(abs(diag(normalizedGram) - 1));
diagnostics.offdiagSpecError = norm(symmetrize(offdiag),2);
diagnostics.conditionNumber = cond(Phi);
diagnostics.normalizedDifference = normalizedDifference;
end

function normalized = normalizedMatrix(matrix,Gamma0)
scales = 1./sqrt(diag(Gamma0));
normalized = scales.*matrix.*scales.';
end

function value = symmetrize(value)
value = (value + value.')/2;
end

function x = leastSquaresSolution(A,b)
if exist("lsqminnorm","file") == 2
    x = lsqminnorm(A,b);
else
    warningState = warning("off","MATLAB:rankDeficientMatrix");
    cleanup = onCleanup(@() warning(warningState));
    x = A\b;
end
end

function rows = appendRow(rows,row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row;
end
end

function printReferenceSummary(reference)
fprintf('\nReference 100-point diagnostics\n');
fprintf('  full-active nonnegative fitted weights evaluated on r=98: max(E_F,E_G)=%10.3e, neg dz=%d, min dz=%10.3e\n', ...
    reference.fullActiveDiagnostics.maxSpecError,reference.fullActiveDiagnostics.negativeDzCount,reference.fullActiveDiagnostics.minDz);
fprintf('  full-active nonnegative fitted weights evaluated on r=20: max(E_F,E_G)=%10.3e, neg dz=%d, min dz=%10.3e\n', ...
    reference.fullOnTargetDiagnostics.maxSpecError,reference.fullOnTargetDiagnostics.negativeDzCount,reference.fullOnTargetDiagnostics.minDz);
fprintf('  target-band nonnegative fitted weights evaluated on r=20: max(E_F,E_G)=%10.3e, neg dz=%d, min dz=%10.3e\n', ...
    reference.targetDiagnostics.maxSpecError,reference.targetDiagnostics.negativeDzCount,reference.targetDiagnostics.minDz);
end

function printRankedRows(rows,nReduced)
rows = rows([rows.nReduced] == nReduced);
[~,order] = sort([rows.maxSpecError]);
rows = rows(order);

fprintf('\nRanked reduced-grid diagnostics for nReduced=%d\n', nReduced);
fprintf('%-23s %-18s %11s %11s %11s %11s %11s %11s %7s\n', ...
    'node policy','weight policy','max E','E_F','E_G','sum dz','min dz','max dz','neg dz');
for iRow = 1:length(rows)
    row = rows(iRow);
    fprintf('%-23s %-18s %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %7d\n', ...
        char(row.nodePolicy),char(row.weightPolicy),row.maxSpecError,row.fSpecError,row.gSpecError, ...
        row.sumDz,row.minDz,row.maxDz,row.negativeDzCount);
end
end

function printAggregationComparison(rows)
fullRows = rows([rows.weightPolicy] == "aggregatedFull");
targetRows = rows([rows.weightPolicy] == "aggregatedTarget");
targetBetterCount = 0;
fullBetterCount = 0;
tieCount = 0;
ratios = zeros(0,1);

for iRow = 1:length(fullRows)
    match = find([targetRows.nReduced] == fullRows(iRow).nReduced & [targetRows.nodePolicy] == fullRows(iRow).nodePolicy,1);
    if isempty(match)
        continue
    end

    ratio = targetRows(match).maxSpecError/fullRows(iRow).maxSpecError;
    ratios(end+1,1) = ratio;
    if ratio < 1
        targetBetterCount = targetBetterCount + 1;
    elseif ratio > 1
        fullBetterCount = fullBetterCount + 1;
    else
        tieCount = tieCount + 1;
    end
end

fprintf('\nAggregated nonnegative reference-weight comparison\n');
fprintf('  target-band aggregate better in %d cases; full-active aggregate better in %d cases; ties=%d.\n', ...
    targetBetterCount,fullBetterCount,tieCount);
if ~isempty(ratios)
    fprintf('  median target/full max-error ratio: %.3e.\n', median(ratios));
end
end

function printWeightPolicySummary(rows)
weightPolicies = ["geometric","aggregatedFull","aggregatedTarget","nonnegativeGramFit"];
bestPoliciesByGrid = strings(1,length(unique([rows.nReduced])));
nReducedValues = unique([rows.nReduced]);
for iReduced = 1:length(nReducedValues)
    gridRows = rows([rows.nReduced] == nReducedValues(iReduced) & ismember([rows.weightPolicy],weightPolicies));
    [~,bestIndex] = min([gridRows.maxSpecError]);
    bestPoliciesByGrid(iReduced) = gridRows(bestIndex).weightPolicy;
end

fprintf('\nWeight-policy summary across node choices\n');
fprintf('%-20s %7s %11s %11s %11s %-23s\n', ...
    'weight policy','wins','best max E','median best','worst best','best node at best E');
for iPolicy = 1:length(weightPolicies)
    policy = weightPolicies(iPolicy);
    policyRows = rows([rows.weightPolicy] == policy);
    bestByGrid = zeros(length(nReducedValues),1);
    bestNodeByGrid = strings(length(nReducedValues),1);
    for iReduced = 1:length(nReducedValues)
        gridRows = policyRows([policyRows.nReduced] == nReducedValues(iReduced));
        [bestByGrid(iReduced),bestIndex] = min([gridRows.maxSpecError]);
        bestNodeByGrid(iReduced) = gridRows(bestIndex).nodePolicy;
    end

    [bestOverall,bestGridIndex] = min(bestByGrid);
    winCount = sum(bestPoliciesByGrid == policy);
    fprintf('%-20s %7d %11.3e %11.3e %11.3e %-23s\n', ...
        char(policy),winCount,bestOverall,median(bestByGrid),max(bestByGrid),char(bestNodeByGrid(bestGridIndex)));
end

fprintf('\nBest policy by reduced grid\n');
fprintf('%8s %-20s %-23s %11s %11s %11s\n', 'nReduced','weight policy','node policy','max E','E_F','E_G');
for iReduced = 1:length(nReducedValues)
    gridRows = rows([rows.nReduced] == nReducedValues(iReduced) & ismember([rows.weightPolicy],weightPolicies));
    [~,bestIndex] = min([gridRows.maxSpecError]);
    row = gridRows(bestIndex);
    fprintf('%8d %-20s %-23s %11.3e %11.3e %11.3e\n', ...
        row.nReduced,char(row.weightPolicy),char(row.nodePolicy),row.maxSpecError,row.fSpecError,row.gSpecError);
end
end

function plotSummaryRows(rows,nReducedList)
comboLabels = unique([rows.comboLabel],"stable");
figure('Name','Reduced-grid Parseval error summary')
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')

nexttile
hold on
for iLabel = 1:length(comboLabels)
    label = comboLabels(iLabel);
    comboRows = rows([rows.comboLabel] == label);
    [nValues,order] = sort([comboRows.nReduced]);
    values = [comboRows.maxSpecError];
    semilogy(nValues,max(values(order),eps),'-o','DisplayName',char(label));
end
xlabel('reduced grid points')
ylabel('max(E_F,E_G)')
title('Worst retained-band Parseval error')
legend('Location','eastoutside')
grid on

nexttile
hold on
for iLabel = 1:length(comboLabels)
    label = comboLabels(iLabel);
    comboRows = rows([rows.comboLabel] == label);
    [nValues,order] = sort([comboRows.nReduced]);
    values = [comboRows.fSpecError];
    semilogy(nValues,max(values(order),eps),'-o','DisplayName',char(label));
end
xlabel('reduced grid points')
ylabel('E_F')
title('F-mode Parseval error')
grid on

nexttile
hold on
for iLabel = 1:length(comboLabels)
    label = comboLabels(iLabel);
    comboRows = rows([rows.comboLabel] == label);
    [nValues,order] = sort([comboRows.nReduced]);
    values = [comboRows.gSpecError];
    semilogy(nValues,max(values(order),eps),'-o','DisplayName',char(label));
end
xlabel('reduced grid points')
ylabel('E_G')
title('G-mode Parseval error')
grid on

end

function plotReducedWeights(rows,nPlotList)
figure('Name','Reduced-grid increment comparisons')
tiledlayout(1,length(nPlotList),'TileSpacing','compact','Padding','compact')

for iPlot = 1:length(nPlotList)
    nPlot = nPlotList(iPlot);
    nexttile
    hold on
    plotRows = rows([rows.nReduced] == nPlot);
    weightPolicies = unique([plotRows.weightPolicy],"stable");
    for iPolicy = 1:length(weightPolicies)
        policy = weightPolicies(iPolicy);
        policyRows = plotRows([plotRows.weightPolicy] == policy);
        [~,bestIndex] = min([policyRows.maxSpecError]);
        row = policyRows(bestIndex);
        plot(row.dz,row.z,'-o','DisplayName',char(row.comboLabel));
    end
    xlabel('dz (m)')
    ylabel('z (m)')
    title(sprintf('Best rows by weight policy, n=%d', nPlot))
    legend('Location','best')
    grid on
end
end

function plotBestHeatmaps(rows,nReduced)
plotRows = rows([rows.nReduced] == nReduced);
aggregatedRows = plotRows(startsWith([plotRows.weightPolicy],"aggregated"));
refitRows = plotRows([plotRows.weightPolicy] == "nonnegativeGramFit");
if isempty(aggregatedRows) || isempty(refitRows)
    return
end

[~,aggregatedIndex] = min([aggregatedRows.maxSpecError]);
[~,refitIndex] = min([refitRows.maxSpecError]);
aggregatedRow = aggregatedRows(aggregatedIndex);
refitRow = refitRows(refitIndex);

figure('Name',sprintf('Normalized Gram-error heatmaps at nReduced=%d', nReduced))
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
plotHeatmap(aggregatedRow.fNormalizedDifference,"F aggregate: " + aggregatedRow.comboLabel)
plotHeatmap(aggregatedRow.gNormalizedDifference,"G aggregate: " + aggregatedRow.comboLabel)
plotHeatmap(refitRow.fNormalizedDifference,"F refit: " + refitRow.comboLabel)
plotHeatmap(refitRow.gNormalizedDifference,"G refit: " + refitRow.comboLabel)
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
