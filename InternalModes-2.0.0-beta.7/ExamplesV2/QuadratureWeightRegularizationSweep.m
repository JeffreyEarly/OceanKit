%% Investigation: geometric regularization of fixed-point quadrature weights
% This sweep compares geometric control-volume weights, the production
% normalized-Gram fit, and a dimension-normalized regularized Gram fit. It
% intentionally leaves the production objective unchanged. The resulting
% tables remain in the workspace for inspection. The default run is a
% representative developer smoke matrix suitable for the stacked F/G fit.
% Set `runExtendedSweep=true` in the workspace before running this script to
% restore the original 32/64/120-mode and sensitivity investigation.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

if ~exist("runExtendedSweep","var")
    runExtendedSweep = false;
end
if runExtendedSweep
    lambdaValues = [0 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2];
    nModesValues = [4 8 16 32];
    pointKinds = ["modeRoots" "uniform" "surfaceClustered"];
else
    lambdaValues = [0 1e-8 1e-4];
    nModesValues = [4 8];
    pointKinds = ["modeRoots" "surfaceClustered"];
end
configurations = regularizationConfigurations();

fprintf("Running the main quadrature-weight regularization sweep...\n");
sweepResults = runConfigurationSweep(configurations,nModesValues,pointKinds,lambdaValues);

if runExtendedSweep
    fprintf("Running the high-order exponential-WKB mode-root cases...\n");
    highOrderConfiguration = configurations([configurations.name] == "hydrostaticGExponentialWKB");
    highOrderConfiguration.nEVP = 384;
    highOrderResults64 = runConfigurationSweep(highOrderConfiguration,64,"modeRoots",lambdaValues);
    highOrderResults120 = runConfigurationSweep(highOrderConfiguration,120,"modeRoots",lambdaValues);
    highOrderResults = [highOrderResults64;highOrderResults120];

    fprintf("Running point-perturbation and solver-resolution sensitivity cases...\n");
    sensitivityResults = runSensitivitySweep(highOrderConfiguration,[16 32],[160 256 384],[0 1e-4 1e-2],lambdaValues);
    sweepResults = [sweepResults;highOrderResults;sensitivityResults];
end

[candidateSummary,familySummary,pointSummary,modeCountSummary] = summarizeCandidates(sweepResults,lambdaValues(2:end));

fprintf("\nUniversal regularization candidates\n");
disp(candidateSummary)
fprintf("\nFamily-level regularization candidates\n");
disp(familySummary)
if any(candidateSummary.robust)
    robustCandidates = candidateSummary(candidateSummary.robust,:);
    [~,bestIndex] = min(robustCandidates.weightDisplacementP95);
    fprintf("Recommended universal regularization strength: %g\n",robustCandidates.lambda(bestIndex));
else
    fprintf("No regularization strength satisfies the universal robustness criteria.\n");
end

successfulRegularized = sweepResults(sweepResults.succeeded & sweepResults.exitFlag > 0 & sweepResults.rule == "regularized" & ~sweepResults.isStressCase,:);
figure(Name="Quadrature-weight regularization sweep",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
loglog(candidateSummary.lambda,candidateSummary.maximumGramErrorRatio,"o-",LineWidth=1.2)
hold on
yline(1.10,"--",LineWidth=0.8)
hold off
grid on
xlabel("regularization strength \lambda")
ylabel("worst Gram-error ratio")
title("In-space fidelity")

nexttile
loglog(candidateSummary.lambda,candidateSummary.maximumHeldOutLeakageRatio,"o-",LineWidth=1.2)
hold on
yline(1.25,"--",LineWidth=0.8)
hold off
grid on
xlabel("regularization strength \lambda")
ylabel("worst held-out leakage ratio")
title("Out-of-space leakage")

nexttile
loglog(successfulRegularized.gramOperatorError,max(successfulRegularized.relativeWeightDisplacementRMS,eps),".",MarkerSize=9)
grid on
xlabel("normalized Gram operator error")
ylabel("RMS relative weight displacement")
title("Accuracy versus geometric displacement")

pointGroups = successfulRegularized.pointKind;
pointGroups(startsWith(pointGroups,"modeRootsPerturbed")) = "modeRootsPerturbed";
figure(Name="Grouped quadrature-weight Pareto plots",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");
groupedParetoPanel(successfulRegularized,successfulRegularized.problemFamily,"Problem family")
groupedParetoPanel(successfulRegularized,pointGroups,"Point construction")
groupedParetoPanel(successfulRegularized,"n_m = " + string(successfulRegularized.nModes),"Mode count")
groupedParetoPanel(successfulRegularized,"lambda = " + string(successfulRegularized.lambda),"Regularization strength")

%% Sweep helpers
function configurations = regularizationConfigurations()
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
constantN2 = @(z) N0*N0*ones(size(z));
exponentialN2 = @(z) N0*N0*exp(2*z/b);
unitDomain = [-1 0];

configurations = repmat(struct("name","","family","","zDomain",[0 1],"coordinateKind","z","nEVP",192,"isStressCase",false,"factory",[]),1,0);
configurations(end+1) = configuration("canonicalConstant","canonical",unitDomain,"z",192,false,@() IMEigenvalueProblem(zDomain=unitDomain,p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet()));
configurations(end+1) = configuration("canonicalVariableMetric","canonical",unitDomain,"z",192,false,@() IMEigenvalueProblem(zDomain=unitDomain,p=1,q=0,r=@(z) exp(2*z),surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet()));
configurations(end+1) = configuration("canonicalSurfaceEndpointMetric","endpointMetric",unitDomain,"z",192,false,@() IMEigenvalueProblem(zDomain=unitDomain,p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0),bottomBoundary=IMBoundaryCondition.dirichlet()));
configurations(end+1) = configuration("hydrostaticGConstant","hydrostaticG",zDomain,"z",192,false,@() IMInternalModes.hydrostaticGModes(N2=constantN2,zDomain=zDomain,g=g));
configurations(end+1) = configuration("hydrostaticGExponentialZ","hydrostaticG",zDomain,"z",256,false,@() IMInternalModes.hydrostaticGModes(N2=exponentialN2,zDomain=zDomain,g=g));
configurations(end+1) = configuration("hydrostaticGExponentialWKB","hydrostaticG",zDomain,"wkb",256,false,@() IMInternalModes.hydrostaticGModes(N2=exponentialN2,zDomain=zDomain,g=g));
configurations(end+1) = configuration("hydrostaticGExponentialDensity","hydrostaticG",zDomain,"density",256,false,@() IMInternalModes.hydrostaticGModes(N2=exponentialN2,zDomain=zDomain,g=g));
configurations(end+1) = configuration("hydrostaticFExponential","hydrostaticF",zDomain,"wkb",256,false,@() IMInternalModes.hydrostaticFModes(N2=exponentialN2,zDomain=zDomain,g=g));
configurations(end+1) = configuration("fixedWavenumberWaveSignedMetric","waveSignedMetric",zDomain,"wkb",256,true,@() IMInternalModes.waveModesAtWavenumber(N2=exponentialN2,zDomain=zDomain,k=2*pi/1e5,f0=8e-4,g=g));
end

function value = configuration(name,family,zDomain,coordinateKind,nEVP,isStressCase,factory)
value = struct("name",string(name),"family",string(family),"zDomain",zDomain,"coordinateKind",string(coordinateKind),"nEVP",nEVP,"isStressCase",isStressCase,"factory",factory);
end

function results = runConfigurationSweep(configurations,nModesValues,pointKinds,lambdaValues)
records = repmat(emptyRecord(),0,1);
for iConfiguration = 1:numel(configurations)
    configuration = configurations(iConfiguration);
    fprintf("  %s (%s coordinates)\n",configuration.name,configuration.coordinateKind);
    nAvailableModes = max(nModesValues) + 2;
    try
        solver = IMSolverSpectral(nEVP=configuration.nEVP,coordinateKind=configuration.coordinateKind);
        basisSet = solver.solveEVP(configuration.factory(),nModes=nAvailableModes);
    catch exception
        for nModes = nModesValues
            for pointKind = pointKinds
                records = [records; failureRecords(configuration,nModes,pointKind,lambdaValues,0,exception)];
            end
        end
        continue
    end

    for nModes = nModesValues
        fprintf("    %d modes\n",nModes);
        for pointKind = pointKinds
            try
                z = pointsForKind(basisSet,nModes,pointKind);
                records = [records; evaluatePointSet(configuration,basisSet,z,nModes,pointKind,lambdaValues,0)];
            catch exception
                records = [records; failureRecords(configuration,nModes,pointKind,lambdaValues,0,exception)];
            end
        end
    end
end
results = struct2table(records);
end

function results = runSensitivitySweep(configuration,nModesValues,nEVPValues,perturbationValues,lambdaValues)
records = repmat(emptyRecord(),0,1);
for nEVP = nEVPValues
    localConfiguration = configuration;
    localConfiguration.nEVP = nEVP;
    try
        solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=localConfiguration.coordinateKind);
        basisSet = solver.solveEVP(localConfiguration.factory(),nModes=max(nModesValues)+2);
    catch exception
        for nModes = nModesValues
            for perturbation = perturbationValues
                pointKind = "modeRootsPerturbed" + string(perturbation);
                records = [records; failureRecords(localConfiguration,nModes,pointKind,lambdaValues,perturbation,exception)];
            end
        end
        continue
    end

    for nModes = nModesValues
        zModeRoot = basisSet.pointsFromModeRoots(nModes=nModes);
        for perturbation = perturbationValues
            pointKind = "modeRootsPerturbed" + string(perturbation);
            z = perturbInteriorPoints(zModeRoot,perturbation);
            try
                records = [records; evaluatePointSet(localConfiguration,basisSet,z,nModes,pointKind,lambdaValues,perturbation)];
            catch exception
                records = [records; failureRecords(localConfiguration,nModes,pointKind,lambdaValues,perturbation,exception)];
            end
        end
    end
end
results = struct2table(records);
end

function z = pointsForKind(basisSet,nModes,pointKind)
zDomain = basisSet.zDomain;
depth = diff(zDomain);
switch pointKind
    case "modeRoots"
        z = basisSet.pointsFromModeRoots(nModes=nModes);
    case "uniform"
        z = linspace(zDomain(1),zDomain(2),2*nModes+1).';
    case "surfaceClustered"
        sigma = linspace(0,1,3*nModes).';
        z = zDomain(1) + depth*(1 - (1-sigma).^2);
    otherwise
        error("QuadratureWeightRegularizationSweep:UnknownPointKind","Unknown point kind %s.",pointKind)
end
end

function z = perturbInteriorPoints(z,relativeAmplitude)
if relativeAmplitude == 0 || length(z) <= 2
    return
end
spacing = diff(z);
localSpacing = min(spacing(1:end-1),spacing(2:end));
phase = sin((1:length(localSpacing)).'*sqrt(2));
z(2:end-1) = z(2:end-1) + relativeAmplitude*localSpacing.*phase;
if any(diff(z) <= 0)
    error("QuadratureWeightRegularizationSweep:InvalidPerturbation","The requested perturbation changed the point ordering.")
end
end

function records = evaluatePointSet(configuration,basisSet,z,nModes,pointKind,lambdaValues,perturbation)
startTime = tic;
[~,pureFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes);
pureElapsed = toc(startTime);
caseID = configuration.name + "|" + configuration.coordinateKind + "|" + pointKind + "|" + string(configuration.nEVP) + "|" + string(nModes);
records = diagnosticRecord(configuration,basisSet,pureFit.geometricWeights,pureFit.geometricTransform,pureFit,nModes,pointKind,"geometric",NaN,perturbation,caseID,0,NaN);
records(end+1,1) = diagnosticRecord(configuration,basisSet,pureFit.weights,pureFit.transform,pureFit,nModes,pointKind,"pure",0,perturbation,caseID,pureElapsed,pureFit.exitFlag);

for lambda = lambdaValues(lambdaValues > 0)
    objective = @(context) dimensionNormalizedRegularizedObjective(context,lambda,nModes);
    startTime = tic;
    try
        [~,fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,objective=objective);
        elapsedSeconds = toc(startTime);
        records(end+1,1) = diagnosticRecord(configuration,basisSet,fit.weights,fit.transform,pureFit,nModes,pointKind,"regularized",lambda,perturbation,caseID,elapsedSeconds,fit.exitFlag);
    catch exception
        elapsedSeconds = toc(startTime);
        records(end+1,1) = failureRecord(configuration,nModes,pointKind,"regularized",lambda,perturbation,caseID,elapsedSeconds,exception);
    end
end
end

function specification = dimensionNormalizedRegularizedObjective(context,lambda,nModes)
nSamples = length(context.z);
regularizationMatrix = sqrt(lambda/nSamples)*diag(1./context.geometricWeights);
regularizationTarget = sqrt(lambda/nSamples)*ones(nSamples,1);
specification = struct("A",[context.normalizedGramA/nModes; regularizationMatrix],"b",[context.normalizedGramB/nModes; regularizationTarget],"name","dimensionNormalizedRegularizedGram");
end

function record = diagnosticRecord(configuration,basisSet,weights,transform,pureFit,nModes,pointKind,rule,lambda,perturbation,caseID,elapsedSeconds,exitFlag)
relativeDisplacement = (weights-pureFit.geometricWeights)./pureFit.geometricWeights;
weightRatios = weights./pureFit.geometricWeights;

if isa(transform,"IMInternalModesDiscreteTransform")
    gramFrobeniusError = 0;
    gramOperatorError = 0;
    maximumDiagonalGramError = 0;
    maximumOffDiagonalGramError = 0;
    heldOutLeakage = 0;
    inverseMatrixConditionNumber = 0;
    gramConditionNumber = 0;
    roundTripError = 0;
    for variable = transform.availableVariables
        active = transform.activeModeMask(variable=variable);
        target = transform.targetGramMatrix(variable=variable);
        gram = transform.gramMatrix(variable=variable);
        scale = 1./sqrt(abs(diag(target(active,active))));
        gramDifference = scale.*(gram(active,active)-target(active,active)).*scale.';
        offDiagonalDifference = gramDifference-diag(diag(gramDifference));
        gramFrobeniusError = max(gramFrobeniusError,norm(gramDifference,"fro"));
        gramOperatorError = max(gramOperatorError,norm(gramDifference,2));
        maximumDiagonalGramError = max(maximumDiagonalGramError,max(abs(diag(gramDifference))));
        maximumOffDiagonalGramError = max(maximumOffDiagonalGramError,max(abs(offDiagonalDifference),[],"all"));
        heldOutValues = internalModeValues(basisSet,variable,transform.z);
        nHeldOut = min(2,size(heldOutValues,2)-nModes);
        if nHeldOut > 0
            heldOutCoefficients = transform.transformForward(heldOutValues(:,nModes+(1:nHeldOut)),variable=variable);
            heldOutLeakage = max(heldOutLeakage,max(vecnorm(heldOutCoefficients,2,1)));
        end
        inverseMatrixConditionNumber = max(inverseMatrixConditionNumber,transform.inverseMatrixConditionNumber(variable=variable));
        gramConditionNumber = max(gramConditionNumber,transform.gramConditionNumber(variable=variable));
        roundTripError = max(roundTripError,transform.roundTripError(variable=variable));
    end
else
    targetNorms = diag(transform.targetGramMatrix);
    scale = 1./sqrt(abs(targetNorms));
    gramDifference = scale.*(transform.gramMatrix-transform.targetGramMatrix).*scale.';
    offDiagonalDifference = gramDifference-diag(diag(gramDifference));
    gramFrobeniusError = norm(gramDifference,"fro");
    gramOperatorError = norm(gramDifference,2);
    maximumDiagonalGramError = max(abs(diag(gramDifference)));
    maximumOffDiagonalGramError = max(abs(offDiagonalDifference),[],"all");
    heldOutValues = basisSet.u(transform.z);
    nHeldOut = min(2,size(heldOutValues,2)-nModes);
    if nHeldOut > 0
        heldOutCoefficients = transform.forwardMatrix*heldOutValues(:,nModes+(1:nHeldOut));
        heldOutLeakage = max(vecnorm(heldOutCoefficients,2,1));
    else
        heldOutLeakage = NaN;
    end
    inverseMatrixConditionNumber = transform.inverseMatrixConditionNumber;
    gramConditionNumber = transform.gramConditionNumber;
    roundTripError = transform.roundTripError;
end

record = emptyRecord();
record.caseID = caseID;
record.problem = configuration.name;
record.problemFamily = configuration.family;
record.coordinateKind = configuration.coordinateKind;
record.pointKind = string(pointKind);
record.nEVP = configuration.nEVP;
record.nModes = nModes;
record.nSamples = length(transform.z);
record.lambda = lambda;
record.rule = string(rule);
record.isStressCase = configuration.isStressCase;
record.perturbation = perturbation;
record.succeeded = true;
record.gramFrobeniusError = gramFrobeniusError;
record.gramOperatorError = gramOperatorError;
record.maximumDiagonalGramError = maximumDiagonalGramError;
record.maximumOffDiagonalGramError = maximumOffDiagonalGramError;
record.heldOutLeakage = heldOutLeakage;
record.relativeWeightDisplacementRMS = norm(relativeDisplacement)/sqrt(length(weights));
record.relativeWeightDisplacementMaximum = max(abs(relativeDisplacement));
record.minimumWeight = min(weights);
record.maximumWeight = max(weights);
record.minimumWeightRatio = min(weightRatios);
record.maximumWeightRatio = max(weightRatios);
record.zeroWeightCount = sum(abs(weights) <= 1e-12*max(1,diff(configuration.zDomain)));
record.depthError = sum(weights)-diff(configuration.zDomain);
record.inverseMatrixConditionNumber = inverseMatrixConditionNumber;
record.gramConditionNumber = gramConditionNumber;
record.roundTripError = roundTripError;
record.exitFlag = exitFlag;
record.elapsedSeconds = elapsedSeconds;
end

function values = internalModeValues(basisSet,variable,z)
if variable == "F"
    values = basisSet.F(z);
else
    values = basisSet.G(z);
end
end

function records = failureRecords(configuration,nModes,pointKind,lambdaValues,perturbation,exception)
caseID = configuration.name + "|" + configuration.coordinateKind + "|" + pointKind + "|" + string(configuration.nEVP) + "|" + string(nModes);
records = failureRecord(configuration,nModes,pointKind,"geometric",NaN,perturbation,caseID,0,exception);
records(end+1,1) = failureRecord(configuration,nModes,pointKind,"pure",0,perturbation,caseID,0,exception);
for lambda = lambdaValues(lambdaValues > 0)
    records(end+1,1) = failureRecord(configuration,nModes,pointKind,"regularized",lambda,perturbation,caseID,0,exception);
end
end

function record = failureRecord(configuration,nModes,pointKind,rule,lambda,perturbation,caseID,elapsedSeconds,exception)
record = emptyRecord();
record.caseID = caseID;
record.problem = configuration.name;
record.problemFamily = configuration.family;
record.coordinateKind = configuration.coordinateKind;
record.pointKind = string(pointKind);
record.nEVP = configuration.nEVP;
record.nModes = nModes;
record.lambda = lambda;
record.rule = string(rule);
record.isStressCase = configuration.isStressCase;
record.perturbation = perturbation;
record.succeeded = false;
record.errorIdentifier = string(exception.identifier);
record.errorMessage = string(exception.message);
record.elapsedSeconds = elapsedSeconds;
end

function record = emptyRecord()
record = struct("caseID","","problem","","problemFamily","","coordinateKind","","pointKind","","nEVP",NaN,"nModes",NaN,"nSamples",NaN,"lambda",NaN,"rule","","isStressCase",false,"perturbation",NaN,"succeeded",false,"gramFrobeniusError",NaN,"gramOperatorError",NaN,"maximumDiagonalGramError",NaN,"maximumOffDiagonalGramError",NaN,"heldOutLeakage",NaN,"relativeWeightDisplacementRMS",NaN,"relativeWeightDisplacementMaximum",NaN,"minimumWeight",NaN,"maximumWeight",NaN,"minimumWeightRatio",NaN,"maximumWeightRatio",NaN,"zeroWeightCount",NaN,"depthError",NaN,"inverseMatrixConditionNumber",NaN,"gramConditionNumber",NaN,"roundTripError",NaN,"exitFlag",NaN,"elapsedSeconds",NaN,"errorIdentifier","","errorMessage","");
end

function [candidateSummary,familySummary,pointSummary,modeCountSummary] = summarizeCandidates(results,lambdaValues)
positiveResults = results(~results.isStressCase,:);
candidateSummary = candidateRows(positiveResults,lambdaValues,"allPositiveDefinite");
familySummary = groupedCandidateRows(positiveResults,lambdaValues,positiveResults.problemFamily);
pointGroups = positiveResults.pointKind;
pointGroups(startsWith(pointGroups,"modeRootsPerturbed")) = "modeRootsPerturbed";
pointSummary = groupedCandidateRows(positiveResults,lambdaValues,pointGroups);
modeCountSummary = groupedCandidateRows(positiveResults,lambdaValues,"nModes=" + string(positiveResults.nModes));
end

function summary = groupedCandidateRows(results,lambdaValues,groups)
groupNames = unique(groups,"stable");
summary = candidateRows(results(groups == groupNames(1),:),lambdaValues,groupNames(1));
for iGroup = 2:numel(groupNames)
    summary = [summary; candidateRows(results(groups == groupNames(iGroup),:),lambdaValues,groupNames(iGroup))];
end
end

function summary = candidateRows(results,lambdaValues,groupName)
records = repmat(struct("group","","lambda",NaN,"caseCount",0,"successfulCaseCount",0,"maximumGramErrorRatio",NaN,"maximumHeldOutLeakageRatio",NaN,"weightDisplacementP95",NaN,"robust",false),numel(lambdaValues),1);
pureResults = results(results.rule == "pure",:);
for iLambda = 1:numel(lambdaValues)
    lambda = lambdaValues(iLambda);
    regularizedResults = results(results.rule == "regularized" & results.lambda == lambda,:);
    [isMatched,pureIndex] = ismember(regularizedResults.caseID,pureResults.caseID);
    successful = false(height(regularizedResults),1);
    successful(isMatched) = regularizedResults.succeeded(isMatched) & regularizedResults.exitFlag(isMatched) > 0 & pureResults.succeeded(pureIndex(isMatched)) & pureResults.exitFlag(pureIndex(isMatched)) > 0;
    matchedPureIndex = pureIndex(successful);
    gramRatio = (regularizedResults.gramOperatorError(successful)+1e-12)./(pureResults.gramOperatorError(matchedPureIndex)+1e-12);
    leakageRatio = (regularizedResults.heldOutLeakage(successful)+1e-12)./(pureResults.heldOutLeakage(matchedPureIndex)+1e-12);
    records(iLambda).group = string(groupName);
    records(iLambda).lambda = lambda;
    records(iLambda).caseCount = height(regularizedResults);
    records(iLambda).successfulCaseCount = sum(successful);
    if any(successful)
        records(iLambda).maximumGramErrorRatio = max(gramRatio);
        records(iLambda).maximumHeldOutLeakageRatio = max(leakageRatio);
        records(iLambda).weightDisplacementP95 = percentile(regularizedResults.relativeWeightDisplacementRMS(successful),0.95);
        records(iLambda).robust = all(successful) && max(gramRatio) <= 1.10 && max(leakageRatio) <= 1.25;
    end
end
summary = struct2table(records);
end

function value = percentile(values,fraction)
values = sort(values(:));
if isempty(values)
    value = NaN;
    return
end
position = 1 + (length(values)-1)*fraction;
lowerIndex = floor(position);
upperIndex = ceil(position);
if lowerIndex == upperIndex
    value = values(lowerIndex);
else
    value = values(lowerIndex) + (position-lowerIndex)*(values(upperIndex)-values(lowerIndex));
end
end

function groupedParetoPanel(results,groups,titleText)
nexttile
groupNames = unique(groups,"stable");
colors = lines(numel(groupNames));
hold on
for iGroup = 1:numel(groupNames)
    isGroup = groups == groupNames(iGroup);
    loglog(max(results.gramOperatorError(isGroup),eps),max(results.relativeWeightDisplacementRMS(isGroup),eps),".",Color=colors(iGroup,:),MarkerSize=9)
end
hold off
set(gca,XScale="log",YScale="log")
grid on
xlabel("normalized Gram operator error")
ylabel("RMS relative weight displacement")
title(titleText)
legend(groupNames,Location="best",Interpreter="none")
end
