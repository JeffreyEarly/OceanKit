function [transform, assessment] = discreteTransform(self, options)
% Build an aligned internal-mode F/G transform through the compatibility API.
%
% The point rule is shared, but each requested variable retains its own
% sampled metric, active columns, target Gram matrix, and forward
% projection. Omitted `variables` selects every directly representable
% channel in canonical order `F`, `G`. Point-limited construction still
% uses roots of the next mode in the EVP's solved formulation.
%
% Supplying explicit `z` and `weights` keeps that quadrature rule unchanged.
% The method assesses every leading mode set on the one rule and returns the
% largest set that passes all enabled policies. The Gram policy must pass
% independently for every requested channel. Optional leakage uses
% same-variable rejected modes. Coupled quadratic aliasing assesses `FF->F`
% and `GG->F` when F is enabled, and `FG->G` when G is enabled.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = discreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested mode-root point count
% - Parameter options.z: explicit increasing sample locations
% - Parameter options.weights: optional explicit shared weights
% - Parameter options.nModes: optional strict family prefix for explicit z
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.gramTolerance: per-channel normalized-Gram tolerance
% - Parameter options.leakageTolerance: optional rejected-mode tolerance
% - Parameter options.quadraticAliasingTolerance: optional coupled-product tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: retained aligned transform
% - Returns assessment: family diagnostics and grid-selection provenance
arguments
    self IMInternalModesBasis
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.weights (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nModes double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.variables (1,:) string = strings(1,0)
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.allowRetainedPrefix (1,1) logical = false
end

IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nPoints,"nPoints");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nModes,"nModes");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nCheckModes,"nCheckModes");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.leakageTolerance,"leakageTolerance");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.quadraticAliasingTolerance,"quadraticAliasingTolerance");

if self.evp.modeFamily == "meanDensityAnomaly" && ~isempty(options.quadraticAliasingTolerance)
    error("IMMeanDensityAnomalyModesBasis:UnavailableQuadraticAliasingPolicy", ...
        "Coupled quadratic aliasing is not defined for mean-density-anomaly modes. Build a Gram-only transform, or use G leakage only when the retained target Gram matrix is positive definite.");
end

hasPointCount = ~isempty(options.nPoints);
hasExplicitPoints = ~isempty(options.z);
hasExplicitWeights = ~isempty(options.weights);
hasExplicitModeCount = ~isempty(options.nModes);
if hasPointCount == hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "Specify exactly one of nPoints or z.");
end
if hasExplicitWeights && ~hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "weights can be supplied only with explicit z.");
end
if hasPointCount && (hasExplicitWeights || hasExplicitModeCount)
    error("IMBasisSet:InvalidDiscretePointSpecification", "The nPoints workflow determines its mode-root candidate band and cannot be combined with weights or nModes.");
end
if options.allowRetainedPrefix && ~hasExplicitModeCount
    error("IMBasisSet:InvalidDiscreteTransformOption", "allowRetainedPrefix is an internal refinement option and requires an explicit nModes band.");
end
if ~hasExplicitWeights && ~hasExplicitModeCount
    [transform,assessment] = self.certifiedDiscreteTransform(nPoints=options.nPoints,z=options.z,variables=options.variables, ...
        gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
        quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes);
    return
end

nAvailableModes = size(self.nativeModes,2);
if hasPointCount
    requestedPointCount = options.nPoints;
    [z,candidateModeCount,gridDesign] = IMDiscreteTransformTools.pointsForExactCount(self,requestedPointCount,nAvailableModes);
else
    z = options.z(:);
    gridDesign = IMDiscreteTransformTools.explicitGridDesign(self,z);
    requestedPointCount = length(z);
    if hasExplicitModeCount
        candidateModeCount = options.nModes;
    else
        candidateModeCount = min(nAvailableModes,length(z));
    end
end
if candidateModeCount > nAvailableModes
    error("IMBasisSet:InvalidDiscreteModeCount", "The basis set contains %d modes, but nModes=%d was requested.",nAvailableModes,candidateModeCount);
end
if isempty(options.variables)
    variables = self.directlyRepresentableDiscreteVariables(z);
else
    variables = self.canonicalDiscreteVariables(options.variables);
end
if isempty(variables)
    error("IMInternalModesBasis:NoAvailableDiscreteTransformVariable", "Neither F nor G has a directly representable sampled metric on the supplied points.");
end

if hasExplicitWeights
    weightFit = [];
    candidateTransform = self.buildInternalModesDiscreteTransform(z,options.weights(:),candidateModeCount,variables);
else
    [~,weightFit] = self.quadratureWeightsForPoints(z=z,nModes=candidateModeCount,variables=variables);
    candidateTransform = weightFit.transform;
end

[prefixData,variableDiagnostics] = preparePrefixData(candidateTransform,variables,~isempty(options.leakageTolerance));
gramError = zeros(candidateModeCount,1);
gramLimitingVariable = strings(candidateModeCount,1);
for nPrefix = 1:candidateModeCount
    errors = arrayfun(@(variable) variableDiagnostics.(char(variable)).gramError(nPrefix),variables);
    [gramError(nPrefix),iLimiting] = max(errors);
    gramLimitingVariable(nPrefix) = variables(iLimiting);
end
gramAccepted = IMDiscreteTransformTools.cumulativeAcceptance(gramError <= options.gramTolerance);
gramPolicy = IMDiscreteTransformTools.policyResult("gram",true,options.gramTolerance,gramError,gramAccepted,candidateModeCount);
gramPolicy.limitingVariable = gramLimitingVariable;

leakageError = nan(candidateModeCount,1);
leakageLimitingVariable = strings(candidateModeCount,1);
leakageLimitingModeNumber = nan(candidateModeCount,1);
if isempty(options.leakageTolerance)
    leakageAccepted = true(candidateModeCount,1);
    leakagePolicy = IMDiscreteTransformTools.policyResult("leakage",false,[],leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = [];
else
    requirePositiveTargets(candidateTransform,variables,"Rejected-mode leakage");
    if isempty(options.nCheckModes)
        nCheckModes = 2*candidateModeCount;
    else
        nCheckModes = options.nCheckModes;
    end
    if nCheckModes <= candidateModeCount
        error("IMBasisSet:InvalidLeakageCheckModeCount", "nCheckModes must be greater than the candidate mode count when leakage assessment is enabled.");
    end
    checkBasis = auxiliaryBasis(self,nCheckModes,candidateModeCount,nAvailableModes);
    [leakageError,leakageLimitingVariable,leakageLimitingModeNumber] = coupledLeakage(prefixData,candidateTransform,checkBasis,variables,nCheckModes);
    leakageAccepted = IMDiscreteTransformTools.cumulativeAcceptance(leakageError <= options.leakageTolerance);
    leakagePolicy = IMDiscreteTransformTools.policyResult("leakage",true,options.leakageTolerance,leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = nCheckModes;
end
leakagePolicy.limitingVariable = leakageLimitingVariable;
leakagePolicy.limitingModeNumber = leakageLimitingModeNumber;

quadraticError = nan(candidateModeCount,1);
quadraticLimitingChannel = strings(candidateModeCount,1);
quadraticLimitingModeNumberI = nan(candidateModeCount,1);
quadraticLimitingModeNumberJ = nan(candidateModeCount,1);
if isempty(options.quadraticAliasingTolerance)
    quadraticAccepted = true(candidateModeCount,1);
    quadraticPolicy = IMDiscreteTransformTools.policyResult("quadraticAliasing",false,[],quadraticError,quadraticAccepted,candidateModeCount);
else
    [quadraticError,quadraticLimitingChannel,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ] = coupledQuadraticAliasing(self,prefixData,candidateTransform,variables);
    quadraticAccepted = IMDiscreteTransformTools.cumulativeAcceptance(quadraticError <= options.quadraticAliasingTolerance);
    quadraticPolicy = IMDiscreteTransformTools.policyResult("quadraticAliasing",true,options.quadraticAliasingTolerance,quadraticError,quadraticAccepted,candidateModeCount);
end
quadraticPolicy.errorNorm = "inducedHilbertMajorant";
quadraticPolicy.projectionPairing = "signedPontryagin";
quadraticPolicy.limitingChannel = quadraticLimitingChannel;
quadraticPolicy.limitingModeNumberI = quadraticLimitingModeNumberI;
quadraticPolicy.limitingModeNumberJ = quadraticLimitingModeNumberJ;

combinedAccepted = gramAccepted & leakageAccepted & quadraticAccepted;
retainedModeCount = find(combinedAccepted,1,"last");
if isempty(retainedModeCount)
    if hasExplicitModeCount && ~options.allowRetainedPrefix
        throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
    end
    error("IMBasisSet:NoAcceptableDiscreteTransformPrefix", "No candidate family mode passes every enabled retained-band policy on this shared rule.");
end
if hasExplicitModeCount && ~options.allowRetainedPrefix && retainedModeCount < candidateModeCount
    throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
end

[limitingPolicy,limitingVariable,retentionReason] = combinedPolicyResult(candidateModeCount,retainedModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
if retainedModeCount == candidateModeCount
    transform = candidateTransform;
else
    transform = candidateTransform.prefixTransform(retainedModeCount);
end
modeCount = (1:candidateModeCount).';
lastModeNumber = candidateTransform.modeNumber(:);
prefixDiagnostics = table(modeCount,lastModeNumber,gramError,gramLimitingVariable,leakageError,leakageLimitingVariable,leakageLimitingModeNumber, ...
    quadraticError,quadraticLimitingChannel,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ,gramAccepted,leakageAccepted,quadraticAccepted,combinedAccepted, ...
    VariableNames=["modeCount","lastModeNumber","gramError","gramLimitingVariable","leakageError","leakageLimitingVariable","leakageLimitingModeNumber", ...
    "quadraticAliasingError","quadraticLimitingChannel","quadraticLimitingModeNumberI","quadraticLimitingModeNumberJ","gramAccepted","leakageAccepted","quadraticAccepted","combinedAccepted"]);
assessment = IMInternalModesDiscreteTransformAssessment(transform=transform,candidateTransform=candidateTransform,weightFit=weightFit, ...
    requestedPointCount=requestedPointCount,prefixDiagnostics=prefixDiagnostics,variableDiagnostics=variableDiagnostics,gramPolicy=gramPolicy, ...
    leakagePolicy=leakagePolicy,quadraticAliasingPolicy=quadraticPolicy,limitingPolicy=limitingPolicy,limitingVariable=limitingVariable,retentionReason=retentionReason,gridDesign=gridDesign);
end

function diagnostics = emptyVariableDiagnostics(nModes,modeNumber)
modeCount = (1:nModes).';
lastModeNumber = modeNumber(:);
activeModeCount = zeros(nModes,1);
gramError = zeros(nModes,1);
roundTripError = zeros(nModes,1);
inverseMatrixConditionNumber = nan(nModes,1);
gramConditionNumber = nan(nModes,1);
sampledRank = zeros(nModes,1);
targetGramIsPositiveDefinite = false(nModes,1);
diagnostics = table(modeCount,lastModeNumber,activeModeCount,gramError,roundTripError,inverseMatrixConditionNumber,gramConditionNumber,sampledRank,targetGramIsPositiveDefinite);
end

function [prefixData,variableDiagnostics] = preparePrefixData(candidateTransform,variables,retainForwardMatrices)
nCandidate = length(candidateTransform.modeNumber);
nSamples = length(candidateTransform.z);
prefixData = struct();
variableDiagnostics = struct();
for variable = variables
    field = char(variable);
    inverseFull = candidateTransform.inverseMatrix(variable=variable);
    metric = candidateTransform.metricMatrix(variable=variable);
    sampledGramFull = candidateTransform.gramMatrix(variable=variable);
    targetGramFull = candidateTransform.targetGramMatrix(variable=variable);
    activeFull = candidateTransform.activeModeMask(variable=variable);
    pairingOperatorFull = inverseFull.'*metric;
    sampledGramRank = zeros(nCandidate,1);
    if retainForwardMatrices
        forwardMatrices = cell(nCandidate,1);
    else
        forwardMatrices = cell(0,1);
    end
    record = emptyVariableDiagnostics(nCandidate,candidateTransform.modeNumber);
    for nPrefix = 1:nCandidate
        inverse = inverseFull(:,1:nPrefix);
        active = activeFull(1:nPrefix);
        sampledGram = sampledGramFull(1:nPrefix,1:nPrefix);
        targetGram = targetGramFull(1:nPrefix,1:nPrefix);
        forward = zeros(nPrefix,nSamples);
        if any(active)
            sampledActive = sampledGram(active,active);
            singularValues = svd(sampledActive);
            rankTolerance = max(size(sampledActive))*eps(max(1,norm(sampledActive,2)));
            prefixSampledGramRank = sum(singularValues > rankTolerance);
            if prefixSampledGramRank < nnz(active)
                forward(active,:) = pinv(sampledActive,rankTolerance)*pairingOperatorFull(active,:);
            else
                forward(active,:) = sampledActive\pairingOperatorFull(active,:);
            end
            targetActive = targetGram(active,active);
            targetNorms = diag(targetActive);
            scale = 1./sqrt(abs(targetNorms));
            gramError = norm(scale.*(sampledActive-targetActive).*scale.',2);
            if prefixSampledGramRank < nnz(active)
                gramError = Inf;
            end
            projector = diag(double(active));
            roundTripError = norm(forward*inverse-projector,2);
            inverseSingularValues = svd(inverse(:,active));
            inverseConditionNumber = singularValueConditionNumber(inverseSingularValues);
            gramConditionNumber = singularValueConditionNumber(singularValues);
            targetTolerance = 100*eps(max(1,norm(targetActive,2)));
            targetIsPositiveDefinite = min(eig(targetActive)) > targetTolerance;
            inverseRankTolerance = max(size(inverse(:,active)))*eps(max(inverseSingularValues));
            sampledRank = sum(inverseSingularValues > inverseRankTolerance);
        else
            prefixSampledGramRank = 0;
            gramError = 0;
            roundTripError = 0;
            inverseConditionNumber = NaN;
            gramConditionNumber = NaN;
            targetIsPositiveDefinite = true;
            sampledRank = 0;
        end
        sampledGramRank(nPrefix) = prefixSampledGramRank;
        if retainForwardMatrices
            forwardMatrices{nPrefix} = forward;
        end
        record.activeModeCount(nPrefix) = nnz(active);
        record.gramError(nPrefix) = gramError;
        record.roundTripError(nPrefix) = roundTripError;
        record.inverseMatrixConditionNumber(nPrefix) = inverseConditionNumber;
        record.gramConditionNumber(nPrefix) = gramConditionNumber;
        record.sampledRank(nPrefix) = sampledRank;
        record.targetGramIsPositiveDefinite(nPrefix) = targetIsPositiveDefinite;
    end
    prefixData.(field) = struct(sampledGramRank=sampledGramRank,forwardMatrices={forwardMatrices});
    variableDiagnostics.(field) = record;
end
end

function value = singularValueConditionNumber(singularValues)
if singularValues(end) == 0
    value = Inf;
else
    value = singularValues(1)/singularValues(end);
end
end

function requirePositiveTargets(transform,variables,policyName)
for variable = variables
    data = transform.channelDiagnostics(variable=variable);
    if ~data.targetGramIsPositiveDefinite
        error("IMBasisSet:UnavailableDiscreteTransformPolicy", "%s requires a positive-definite active target for every participating channel; %s is signed or indefinite. Use the Gram policy alone.",policyName,variable);
    end
end
end

function basis = auxiliaryBasis(source,nCheckModes,candidateModeCount,nAvailableModes)
basis = source;
if nCheckModes > nAvailableModes
    try
        basis = source.solver.solveEVP(source.evp,nModes=nCheckModes);
        basis.normalization = source.normalization;
    catch cause
        exception = MException("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy could not obtain %d source modes.",nCheckModes);
        throw(addCause(exception,cause))
    end
end
if size(basis.nativeModes,2) < nCheckModes
    error("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy requested %d modes, but the solver returned %d.",nCheckModes,size(basis.nativeModes,2));
end
scale = max(1,max(abs([source.eigenvalues(1:candidateModeCount) basis.eigenvalues(1:candidateModeCount)]),[],"all"));
if ~isequal(basis.modeNumber(1:candidateModeCount),source.modeNumber(1:candidateModeCount)) ...
        || any(abs(basis.eigenvalues(1:candidateModeCount)-source.eigenvalues(1:candidateModeCount)) > 1e-8*scale)
    error("IMBasisSet:AuxiliaryModeMismatch", "The leakage-policy solve did not reproduce the candidate family prefix.");
end
end

function [errors,limitingVariable,limitingModeNumber] = coupledLeakage(prefixData,candidateTransform,checkBasis,variables,nCheckModes)
nCandidate = length(candidateTransform.modeNumber);
errors = zeros(nCandidate,1);
limitingVariable = strings(nCandidate,1);
limitingModeNumber = nan(nCandidate,1);
for variable = variables
    values = variableValues(checkBasis,variable,candidateTransform.z);
    values = values(:,1:nCheckModes);
    gram = checkBasis.gramMatrix(variable=variable);
    norms = diag(gram(1:nCheckModes,1:nCheckModes));
    for nPrefix = 1:nCandidate
        rejected = (nPrefix+1):nCheckModes;
        positive = norms(rejected) > 100*eps(max(1,max(abs(norms))));
        rejected = rejected(positive);
        if isempty(rejected)
            value = 0;
            label = NaN;
        else
            data = prefixData.(char(variable));
            coefficients = data.forwardMatrices{nPrefix}*values(:,rejected);
            active = candidateTransform.activeModeMask(variable=variable);
            active = active(1:nPrefix);
            target = candidateTransform.targetGramMatrix(variable=variable);
            target = target(1:nPrefix,1:nPrefix);
            numerator = sum(coefficients(active,:).*(target(active,active)*coefficients(active,:)),1);
            leakage = sqrt(max(0,numerator)./norms(rejected).');
            [value,index] = max(leakage);
            label = checkBasis.modeNumber(rejected(index));
        end
        if value >= errors(nPrefix)
            errors(nPrefix) = value;
            limitingVariable(nPrefix) = variable;
            limitingModeNumber(nPrefix) = label;
        end
    end
end
end

function [errors,limitingChannel,limitingLabelI,limitingLabelJ] = coupledQuadraticAliasing(basisSet,prefixData,candidateTransform,variables)
nCandidate = length(candidateTransform.modeNumber);
errors = zeros(nCandidate,1);
limitingChannel = strings(nCandidate,1);
limitingLabelI = nan(nCandidate,1);
limitingLabelJ = nan(nCandidate,1);
integrationGrid = basisSet.solver.innerProductGrid(basisSet.zDomain);
integrationWeights = basisSet.solver.innerProductWeights(integrationGrid,basisSet.zDomain);
sampleValues = struct(F=candidateTransform.inverseMatrix(variable="F"),G=candidateTransform.inverseMatrix(variable="G"));
normalizationFactors = basisSet.normalizationFactors(basisSet.normalization);
normalizationFactors = normalizationFactors(1:nCandidate);
integrationF = basisSet.rawVariable("F",integrationGrid);
integrationG = basisSet.rawVariable("G",integrationGrid);
integrationValues = struct(F=integrationF(:,1:nCandidate)./normalizationFactors,G=integrationG(:,1:nCandidate)./normalizationFactors);
endpointValues = struct(F=flipud(candidateTransform.endpointValues(variable="F")),G=flipud(candidateTransform.endpointValues(variable="G")));
channels = strings(0,3);
if ismember("F",variables)
    channels = [channels;"F" "F" "F";"G" "G" "F"];
end
if ismember("G",variables)
    channels = [channels;"F" "G" "G"];
end
for iChannel = 1:size(channels,1)
    sourceA = channels(iChannel,1);
    sourceB = channels(iChannel,2);
    target = channels(iChannel,3);
    channelName = sourceA+sourceB+"->"+target;
    context = productProjectionContext(basisSet,integrationGrid,integrationWeights,integrationValues.(char(target)),endpointValues.(char(target)),target,nCandidate);
    sourceAEndpoint = endpointValues.(char(sourceA));
    sourceBEndpoint = endpointValues.(char(sourceB));
    pairIndices = productPairIndices(nCandidate,sourceA == sourceB);
    nPairs = size(pairIndices,1);
    sampledProducts = zeros(length(candidateTransform.z),nPairs);
    productNorms = zeros(nPairs,1);
    pairings = zeros(nCandidate,nPairs);
    zeroTolerance = zeros(nPairs,1);
    maximumBatchEntries = 2^21;
    batchSize = max(1,floor(maximumBatchEntries/length(integrationGrid)));
    for iFirst = 1:batchSize:nPairs
        batch = iFirst:min(iFirst+batchSize-1,nPairs);
        iMode = pairIndices(batch,1);
        jMode = pairIndices(batch,2);
        sampledProducts(:,batch) = sampleValues.(char(sourceA))(:,iMode).*sampleValues.(char(sourceB))(:,jMode);
        continuousProducts = integrationValues.(char(sourceA))(:,iMode).*integrationValues.(char(sourceB))(:,jMode);
        endpointProducts = sourceAEndpoint(:,iMode).*sourceBEndpoint(:,jMode);
        weightedProducts = context.volumeWeights.*continuousProducts;
        productNorms(batch) = (sum(continuousProducts.*weightedProducts,1) ...
            + sum(context.majorantEndpointMetric.*endpointProducts.*endpointProducts,1)).';
        pairings(:,batch) = context.targetValues.'*weightedProducts ...
            + context.targetEndpoint.'*(context.endpointMetric.*endpointProducts);
        zeroTolerance(batch) = 1e3*eps(max(1,productNorms(batch)));
    end
    sampledNorms = vecnorm(sampledProducts,2,1).';
    usable = ~(abs(productNorms) <= zeroTolerance & sampledNorms <= sqrt(zeroTolerance));
    if any(productNorms(usable) <= zeroTolerance(usable))
        error("IMBasisSet:UnavailableDiscreteTransformPolicy", "The %s product channel does not define a positive continuous product norm.",channelName);
    end

    targetInverse = candidateTransform.inverseMatrix(variable=target);
    targetMetric = candidateTransform.metricMatrix(variable=target);
    sampledPairings = targetInverse.'*(targetMetric*sampledProducts);
    activeFull = candidateTransform.activeModeMask(variable=target);
    sampledGramFull = candidateTransform.gramMatrix(variable=target);
    targetGramFull = candidateTransform.targetGramMatrix(variable=target);
    targetMajorantGramFull = candidateTransform.targetMajorantGramMatrix(variable=target);
    data = prefixData.(char(target));
    sampledCoefficients = zeros(0,nPairs);
    continuousCoefficients = zeros(0,nPairs);
    previousSelected = false(nPairs,1);
    previousSampledFullRank = true;
    sampledFactor = zeros(0,0);
    targetFactor = zeros(0,0);
    for nPrefix = 1:nCandidate
        active = activeFull(1:nPrefix);
        selected = usable & max(pairIndices,[],2) <= nPrefix;
        sampledGram = sampledGramFull(1:nPrefix,1:nPrefix);
        targetGram = targetGramFull(1:nPrefix,1:nPrefix);
        targetMajorantGram = targetMajorantGramFull(1:nPrefix,1:nPrefix);
        sampledFullRank = data.sampledGramRank(nPrefix) == nnz(active);
        [sampledCoefficients,sampledFactor] = updateProjectionSolutions(sampledCoefficients,sampledFactor,previousSelected,selected,active, ...
            sampledGram,sampledPairings,sampledFullRank,previousSampledFullRank);
        [continuousCoefficients,targetFactor] = updateProjectionSolutions(continuousCoefficients,targetFactor,previousSelected,selected,active, ...
            targetGram,pairings,true,true);
        if ~any(selected)
            previousSelected = selected;
            previousSampledFullRank = sampledFullRank;
            continue
        end
        selectedPairs = find(selected);
        difference = sampledCoefficients(:,selected)-continuousCoefficients(:,selected);
        targetMajorantActive = targetMajorantGram(active,active);
        numerator = real(sum(conj(difference).*(targetMajorantActive*difference),1));
        values = sqrt(max(0,numerator)./productNorms(selected).');
        [value,iLimiting] = max(values);
        if value >= errors(nPrefix)
            pair = pairIndices(selectedPairs(iLimiting),:);
            errors(nPrefix) = value;
            limitingChannel(nPrefix) = channelName;
            limitingLabelI(nPrefix) = candidateTransform.modeNumber(pair(1));
            limitingLabelJ(nPrefix) = candidateTransform.modeNumber(pair(2));
        end
        previousSelected = selected;
        previousSampledFullRank = sampledFullRank;
    end
end
end

function [solutions,factor] = updateProjectionSolutions(solutions,factor,previousSelected,selected,active,matrix,pairings,isFullRank,previousWasFullRank)
nActive = nnz(active);
nPairs = length(selected);
if nActive == 0
    solutions = zeros(0,nPairs);
    factor = zeros(0,0);
    return
end
activeMatrix = matrix(active,active);
activePairings = pairings(active,:);
previousActiveCount = size(solutions,1);
if isFullRank
    factor = extendPositiveFactor(factor,activeMatrix,previousActiveCount);
else
    factor = zeros(0,0);
end
requiresFullSolve = ~isFullRank || (~previousWasFullRank && previousActiveCount > 0);
if requiresFullSolve
    solutions = zeros(nActive,nPairs);
    if isFullRank
        solutions(:,selected) = solvePositiveSystem(activeMatrix,factor,activePairings(:,selected));
    else
        rankTolerance = max(nActive,1)*eps(max(1,norm(activeMatrix,2)));
        solutions(:,selected) = pinv(activeMatrix,rankTolerance)*activePairings(:,selected);
    end
    return
end
if nActive == previousActiveCount+1
    solutions(nActive,nPairs) = 0;
    if previousActiveCount > 0 && any(previousSelected)
        coupling = activeMatrix(1:end-1,end);
        if isempty(factor)
            previousSolution = activeMatrix(1:end-1,1:end-1)\coupling;
        else
            previousFactor = factor(1:end-1,1:end-1);
            previousSolution = previousFactor\(previousFactor.'\coupling);
        end
        schurComplement = activeMatrix(end,end)-coupling.'*previousSolution;
        lastCoefficient = (activePairings(end,previousSelected)-coupling.'*solutions(1:end-1,previousSelected))/schurComplement;
        solutions(1:end-1,previousSelected) = solutions(1:end-1,previousSelected)-previousSolution*lastCoefficient;
        solutions(end,previousSelected) = lastCoefficient;
    end
elseif nActive ~= previousActiveCount
    error("IMBasisSet:InvalidDiscreteTransformState", "Active internal-mode columns must enter a family prefix one at a time.");
end
newlySelected = selected & ~previousSelected;
if any(newlySelected)
    solutions(:,newlySelected) = solvePositiveSystem(activeMatrix,factor,activePairings(:,newlySelected));
end
end

function factor = extendPositiveFactor(factor,matrix,previousSize)
currentSize = size(matrix,1);
if currentSize == previousSize
    return
end
if currentSize ~= previousSize+1 || (previousSize > 0 && isempty(factor))
    [factor,flag] = chol(matrix);
elseif previousSize == 0
    flag = matrix(1,1) <= 0;
    factor = sqrt(max(0,matrix(1,1)));
else
    column = factor.'\matrix(1:end-1,end);
    diagonalSquared = matrix(end,end)-column.'*column;
    tolerance = 100*eps(max(1,norm(matrix,2)));
    if diagonalSquared > tolerance
        factor = [factor column;zeros(1,previousSize) sqrt(diagonalSquared)];
        flag = 0;
    else
        flag = 1;
    end
end
if flag ~= 0
    [factor,flag] = chol(matrix);
    if flag ~= 0
        factor = zeros(0,0);
    end
end
end

function solutions = solvePositiveSystem(matrix,factor,rightHandSide)
if isempty(factor)
    solutions = matrix\rightHandSide;
else
    solutions = factor\(factor.'\rightHandSide);
end
end

function pairIndices = productPairIndices(nModes,isSymmetric)
if isSymmetric
    nPairs = nModes*(nModes+1)/2;
    pairIndices = zeros(nPairs,2);
    iFirst = 1;
    for iMode = 1:nModes
        jMode = (iMode:nModes).';
        rows = iFirst:(iFirst+length(jMode)-1);
        pairIndices(rows,:) = [repmat(iMode,length(jMode),1) jMode];
        iFirst = rows(end)+1;
    end
else
    pairIndices = [repelem((1:nModes).',nModes) repmat((1:nModes).',nModes,1)];
end
end

function context = productProjectionContext(basisSet,z,integrationWeights,targetValues,targetEndpoint,targetVariable,nCandidate)
spec = basisSet.evp.innerProduct(targetVariable);
coefficientContext = basisSet.evp.contextForSolver(basisSet.solver);
weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,coefficientContext);
if isscalar(weight)
    weight = weight*ones(size(z));
else
    weight = weight(:);
end
endpointZ = [basisSet.zDomain(1);basisSet.zDomain(2)];
[available,reason,~,~,endpointMetric] = basisSet.sampledInternalModesMetric(targetVariable,endpointZ,zeros(2,1));
if ~available
    error("IMBasisSet:UnavailableDiscreteTransformPolicy", "The %s product channel cannot construct its endpoint metric. %s",targetVariable,reason);
end
[majorantAvailable,majorantReason,~,~,majorantEndpointMetric] = basisSet.sampledInternalModesMetric(targetVariable,endpointZ,zeros(2,1),true);
if ~majorantAvailable
    error("IMBasisSet:UnavailableDiscreteTransformPolicy", "The %s product channel cannot construct its endpoint majorant. %s",targetVariable,majorantReason);
end
context = struct(volumeWeights=integrationWeights.*weight,endpointMetric=diag(endpointMetric), ...
    majorantEndpointMetric=diag(majorantEndpointMetric), ...
    targetEndpoint=targetEndpoint(:,1:nCandidate),targetValues=targetValues(:,1:nCandidate));
end

function values = variableValues(basisSet,variable,z)
if variable == "F"
    values = basisSet.F(z);
else
    values = basisSet.G(z);
end
end

function throwStrictModeFailure(nModes,gramPolicy,leakagePolicy,quadraticPolicy)
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
reasons = strings(0,1);
for iPolicy = 1:numel(policies)
    policy = policies{iPolicy};
    if policy.enabled && ~policy.accepted(nModes)
        reasons(end+1) = policy.name+"="+string(policy.error(nModes))+" (tolerance "+string(policy.tolerance)+")"; %#ok<AGROW>
    end
end
error("IMBasisSet:StrictDiscreteModeCountRejected", "The explicit nModes=%d family band is not acceptable: %s.",nModes,join(reasons,"; "));
end

function [policyName,variable,reason] = combinedPolicyResult(candidateCount,retainedCount,gramPolicy,leakagePolicy,quadraticPolicy)
if retainedCount == candidateCount
    policyName = "none";
    variable = "";
    reason = sprintf("All %d candidate family modes pass every enabled policy.",candidateCount);
    return;
end
next = retainedCount+1;
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
for iPolicy = 1:numel(policies)
    policy = policies{iPolicy};
    if policy.enabled && ~policy.accepted(next)
        policyName = policy.name;
        if isfield(policy,"limitingVariable")
            variable = policy.limitingVariable(next);
        elseif isfield(policy,"limitingChannel")
            channel = policy.limitingChannel(next);
            parts = split(channel,"->");
            variable = parts(end);
        else
            variable = "";
        end
        reason = sprintf("Retained %d of %d family modes; %s first rejects the next prefix (%s).",retainedCount,candidateCount,policy.name,policy.reason);
        return;
    end
end
policyName = "combined";
variable = "";
reason = sprintf("Retained %d of %d family modes under the combined policies.",retainedCount,candidateCount);
end
