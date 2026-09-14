function [transform, assessment] = certifiedDiscreteTransform(self, options)
% Select and independently fit a certified aligned F/G family band.
%
% Fresh family-specific weights are fitted for every count considered by
% the Gram search. This prevents a poorly resolved large candidate family
% from contaminating the assessment of a smaller MDA or APV band. Optional
% leakage and coupled quadratic policies can reduce that Gram-certified
% band; each reduction is refitted until the fitted and retained counts
% agree.
%
% The count search fits each candidate independently but evaluates only its
% full-band Gram matrix. For fixed candidate weights, every normalized
% prefix error is a principal submatrix of the full-band error, so its
% spectral norm cannot be larger. The normalized continuous target differs
% from its diagonal signature matrix by a fixed perturbation. When the Gram
% tolerance is smaller than the resulting nonsingularity margin, the same
% bound proves that every prefix Gram matrix has full rank. Larger loose
% tolerances retain a lightweight cumulative-prefix fallback.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = certifiedDiscreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested mode-root point count
% - Parameter options.z: explicit increasing physical sample points
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.gridDesign: optional provenance returned by `modeRootGrid`
% - Parameter options.gramTolerance: per-channel normalized-Gram tolerance
% - Parameter options.leakageTolerance: optional rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional coupled-product tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: independently fitted certified aligned transform
% - Returns assessment: final diagnostics, grid provenance, and count-search table
arguments
    self IMInternalModesBasis
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.variables (1,:) string = strings(1,0)
    options.gridDesign struct = struct.empty
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nPoints,"nPoints");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nCheckModes,"nCheckModes");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.leakageTolerance,"leakageTolerance");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.quadraticAliasingTolerance,"quadraticAliasingTolerance");
if self.evp.modeFamily == "meanDensityAnomaly" && ~isempty(options.quadraticAliasingTolerance)
    error("IMMeanDensityAnomalyModesBasis:UnavailableQuadraticAliasingPolicy", ...
        "Coupled quadratic aliasing is not defined for mean-density-anomaly modes.");
end
hasPointCount = ~isempty(options.nPoints);
hasExplicitPoints = ~isempty(options.z);
if hasPointCount == hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "Specify exactly one of nPoints or z.");
end
if hasPointCount && ~isempty(options.gridDesign)
    error("IMBasisSet:InvalidDiscretePointSpecification", "gridDesign accompanies explicit z; nPoints creates its own mode-root provenance.");
end

nAvailableModes = size(self.nativeModes,2);
if hasPointCount
    [z,maximumModeCount,gridDesign] = IMDiscreteTransformTools.pointsForExactCount(self,options.nPoints,nAvailableModes);
else
    z = options.z(:);
    maximumModeCount = min(nAvailableModes,length(z));
    gridDesign = IMDiscreteTransformTools.validatedGridDesign(self,z,options.gridDesign);
end

[searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = emptySearchColumns();
searchContext = prepareGramSearchContext(self,z,maximumModeCount,options.variables);
gramCertifiedModeCount = [];
for modeCount = maximumModeCount:-1:1
    try
        [candidateGramError,candidateRetainedModeCount,accepted] = fitGramSearchCandidate(searchContext,modeCount,options.gramTolerance);
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,true,accepted,candidateRetainedModeCount,candidateGramError,"","");
        if accepted
            gramCertifiedModeCount = modeCount;
            break
        end
    catch exception
        if ~isRecoverableSearchFailure(exception)
            rethrow(exception)
        end
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
    end
end
if isempty(gramCertifiedModeCount)
    error("IMBasisSet:NoCertifiedDiscreteTransform", "No independently fitted family count passed the Gram tolerance on the supplied points.");
end

if ~isempty(options.leakageTolerance) || ~isempty(options.quadraticAliasingTolerance)
    modeCount = gramCertifiedModeCount;
    while true
        try
            [transform,assessment] = self.discreteTransform(z=z,nModes=modeCount,variables=options.variables, ...
                gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
                quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes, ...
                allowRetainedPrefix=true);
            accepted = assessment.retainedModeCount == modeCount;
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,true,accepted,assessment.retainedModeCount,assessment.gramPolicy.error(end),"","");
            if accepted
                break
            end
            modeCount = assessment.retainedModeCount;
        catch exception
            if ~isRecoverableSearchFailure(exception) || modeCount == 1
                rethrow(exception)
            end
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
            modeCount = modeCount-1;
        end
    end
else
    [transform,assessment] = self.discreteTransform(z=z,nModes=gramCertifiedModeCount,variables=options.variables, ...
        gramTolerance=options.gramTolerance,allowRetainedPrefix=true);
end

certificationSearch = table(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
    VariableNames=["stage","modeCount","fitSucceeded","accepted","retainedModeCount","gramError","failureIdentifier","failureMessage"]);
assessment = assessment.withCertificationMetadata(gridDesign,certificationSearch);
if assessment.weightFitModeCount ~= assessment.retainedModeCount
    error("IMBasisSet:InconsistentCertifiedTransform", "Certified weights must be fitted to the final retained family count.");
end
end

function context = prepareGramSearchContext(basisSet,z,maximumModeCount,requestedVariables)
depth = diff(basisSet.zDomain);
edges = [basisSet.zDomain(1);0.5*(z(1:end-1)+z(2:end));basisSet.zDomain(2)];
geometricWeights = diff(edges);
if isempty(requestedVariables)
    variables = basisSet.directlyRepresentableDiscreteVariables(z);
else
    variables = basisSet.canonicalDiscreteVariables(requestedVariables);
end
if isempty(variables)
    error("IMInternalModesBasis:NoAvailableDiscreteTransformVariable", "Neither F nor G has a directly representable sampled metric on the supplied points.");
end
preparation = basisSet.prepareInternalModesDiscreteTransform(z,geometricWeights,maximumModeCount,variables);
channelObjectives = struct();
for variable = variables
    field = char(variable);
    data = preparation.channelData.(field);
    inverseMatrix = preparation.("inverse"+variable);
    endpointGramMatrix = inverseMatrix.'*data.endpointMetricMatrix*inverseMatrix;
    endpointGramMatrix = 0.5*(endpointGramMatrix+endpointGramMatrix.');
    [A,b,pairs] = normalizedGramSystem(inverseMatrix,data.interiorWeight,data.targetGramMatrix,endpointGramMatrix,data.activeModeMask);
    channelObjectives.(field) = struct(variable=variable,inverseMatrix=inverseMatrix,interiorWeight=data.interiorWeight, ...
        targetGramMatrix=data.targetGramMatrix,endpointGramMatrix=endpointGramMatrix,activeModeMask=data.activeModeMask, ...
        normalizedGramA=A,normalizedGramB=b,normalizedGramModePairs=pairs);
end
prefixRankMargin = targetPrefixRankMargin(channelObjectives,variables);
context = struct(depth=depth,geometricWeights=geometricWeights,variables=variables,channelObjectives=channelObjectives,prefixRankMargin=prefixRankMargin);
end

function [fullBandError,retainedModeCount,accepted] = fitGramSearchCandidate(context,modeCount,tolerance)
objectiveMatrix = zeros(0,length(context.geometricWeights));
objectiveTarget = zeros(0,1);
for variable = context.variables
    channel = context.channelObjectives.(char(variable));
    selected = max(channel.normalizedGramModePairs,[],2) <= modeCount;
    objectiveMatrix = [objectiveMatrix;channel.normalizedGramA(selected,:)]; %#ok<AGROW>
    objectiveTarget = [objectiveTarget;channel.normalizedGramB(selected)]; %#ok<AGROW>
end
[objectiveMatrix,objectiveTarget] = IMDiscreteTransformTools.validateObjectiveSystem(objectiveMatrix,objectiveTarget,length(context.geometricWeights),true);
weights = IMDiscreteTransformTools.fitQuadrature(objectiveMatrix,objectiveTarget,context.geometricWeights,context.depth,true,true);
channelData = candidateGramData(context,weights,modeCount);
fullBandError = gramPrefixError(channelData,context.variables,modeCount);
accepted = fullBandError <= tolerance;
retainedModeCount = NaN;
if accepted
    retainedModeCount = modeCount;
    if tolerance >= context.prefixRankMargin
        [retainedModeCount,accepted] = looseTolerancePrefixAssessment(channelData,context.variables,modeCount,tolerance);
    end
end
end

function margin = targetPrefixRankMargin(channelObjectives,variables)
margin = Inf;
for variable = variables
    channel = channelObjectives.(char(variable));
    active = channel.activeModeMask;
    if ~any(active)
        continue
    end
    targetGram = channel.targetGramMatrix(active,active);
    targetNorms = diag(targetGram);
    scale = 1./sqrt(abs(targetNorms));
    normalizedTarget = scale.*targetGram.*scale.';
    signature = diag(sign(targetNorms));
    margin = min(margin,1-norm(normalizedTarget-signature,2));
end
if ~isfinite(margin)
    margin = 1;
end
end

function channelData = candidateGramData(context,weights,modeCount)
channelData = struct();
for variable = context.variables
    field = char(variable);
    channel = context.channelObjectives.(field);
    inverseMatrix = channel.inverseMatrix(:,1:modeCount);
    interiorMetric = channel.interiorWeight.*weights;
    sampledGram = inverseMatrix.'*(interiorMetric.*inverseMatrix)+channel.endpointGramMatrix(1:modeCount,1:modeCount);
    sampledGram = 0.5*(sampledGram+sampledGram.');
    channelData.(field) = struct(sampledGram=sampledGram,targetGram=channel.targetGramMatrix(1:modeCount,1:modeCount), ...
        activeMask=channel.activeModeMask(1:modeCount));
end
end

function [retainedModeCount,accepted] = looseTolerancePrefixAssessment(channelData,variables,nModes,tolerance)
retainedModeCount = nModes;
accepted = true;
for nPrefix = 1:nModes
    prefixError = gramPrefixError(channelData,variables,nPrefix);
    if ~(prefixError <= tolerance)
        retainedModeCount = nPrefix-1;
        accepted = false;
        break
    end
end
end

function value = gramPrefixError(channelData,variables,nPrefix)
errors = zeros(length(variables),1);
for iVariable = 1:length(variables)
    data = channelData.(char(variables(iVariable)));
    active = data.activeMask(1:nPrefix);
    if ~any(active)
        continue
    end
    sampledGram = data.sampledGram(1:nPrefix,1:nPrefix);
    targetGram = data.targetGram(1:nPrefix,1:nPrefix);
    sampledActive = sampledGram(active,active);
    singularValues = svd(sampledActive);
    rankTolerance = max(size(sampledActive))*eps(max(1,singularValues(1)));
    targetActive = targetGram(active,active);
    targetNorms = diag(targetActive);
    scale = 1./sqrt(abs(targetNorms));
    errors(iVariable) = norm(scale.*(sampledActive-targetActive).*scale.',2);
    if sum(singularValues > rankTolerance) < nnz(active)
        errors(iVariable) = Inf;
    end
end
value = max(errors);
end

function [A,b,pairs] = normalizedGramSystem(inverseMatrix,interiorWeight,targetGramMatrix,endpointGramMatrix,activeMask)
active = find(activeMask);
nRows = length(active)*(length(active)+1)/2;
A = zeros(nRows,size(inverseMatrix,1));
b = zeros(nRows,1);
pairs = zeros(nRows,2);
iRow = 0;
targetNorms = diag(targetGramMatrix);
for iActive = 1:length(active)
    iMode = active(iActive);
    for jActive = iActive:length(active)
        jMode = active(jActive);
        iRow = iRow+1;
        rowFactor = 1;
        if iMode ~= jMode
            rowFactor = sqrt(2);
        end
        scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
        A(iRow,:) = rowFactor*(interiorWeight.*inverseMatrix(:,iMode).*inverseMatrix(:,jMode)).'/scale;
        b(iRow) = rowFactor*(targetGramMatrix(iMode,jMode)-endpointGramMatrix(iMode,jMode))/scale;
        pairs(iRow,:) = [iMode jMode];
    end
end
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = emptySearchColumns()
stage = strings(0,1);
modeCount = zeros(0,1);
fitSucceeded = false(0,1);
accepted = false(0,1);
retainedModeCount = zeros(0,1);
gramError = zeros(0,1);
failureIdentifier = strings(0,1);
failureMessage = strings(0,1);
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = appendSearch( ...
        stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage, ...
        newStage,newModeCount,newFitSucceeded,newAccepted,newRetainedModeCount,newGramError,newFailureIdentifier,newFailureMessage)
stage(end+1,1) = newStage;
modeCount(end+1,1) = newModeCount;
fitSucceeded(end+1,1) = newFitSucceeded;
accepted(end+1,1) = newAccepted;
retainedModeCount(end+1,1) = newRetainedModeCount;
gramError(end+1,1) = newGramError;
failureIdentifier(end+1,1) = newFailureIdentifier;
failureMessage(end+1,1) = newFailureMessage;
end

function tf = isRecoverableSearchFailure(exception)
tf = ismember(string(exception.identifier),["IMBasisSet:DiscreteTransformPolicyFailed", ...
    "IMBasisSet:NoAcceptableDiscreteTransformPrefix","IMBasisSet:QuadratureWeightFitFailed"]);
end
