function [transform, assessment] = discreteTransform(self, options)
% Build a scalar Galerkin transform using the compatibility entry point.
%
% A transform can be built from explicit physical points `z`, or from an
% exact requested point count `nPoints`. The point-count workflow searches
% the mode-root grids returned by `pointsFromModeRoots` and selects the
% largest candidate modal prefix whose grid has exactly the requested
% number of points. `nPoints` and `z` are mutually exclusive.
%
% With neither `weights` nor `nModes`, this method delegates to
% `certifiedDiscreteTransform`: candidate counts receive independent weight
% fits, so a poorly represented large band cannot contaminate a smaller
% band's rule. Prefer that named method in new code. Use
% `fitDiscreteTransform(z=z,modeCount=N)` for a strict exact-band fit, and
% `modeRootGrid` to design and name a shared physical grid explicitly.
%
% Supplying `weights`, or the legacy `nModes` name, retains the lower-level
% fixed-rule behavior. Every leading prefix is assessed on those same
% weights. This path remains useful for diagnosing a caller-owned rule and
% for compatibility, but it does not perform independently refitted count
% selection.
%
% If `nModes` is supplied with `z`, that band is strict: every requested
% mode must pass every enabled policy, or construction throws an error.
% The `nPoints` workflow designs its physical grid from mode roots and then
% uses certified count selection. It therefore does not accept `nModes`.
%
% For prefix $$N$$, define
%
% $$
% E_N=S_N(\Gamma_N-\Gamma_{0,N})S_N,\qquad
% S_N=\operatorname{diag}\!\left(
% \left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}
% \right).
% $$
%
% The always-enabled Gram policy accepts
% $$\lVert E_N\rVert_2\leq\texttt{gramTolerance}$$. The shipped default
% is $$10^{-2}$$. Regression sweeps over constant and representative
% exponential stratifications in physical, WKB, and density coordinates
% found constant cosine rules near roundoff and ordinary exponential rules
% below approximately $$5\times10^{-3}$$; the factor-two margin preserves
% those bands while still rejecting the order-one DCT-I Nyquist failure.
%
% `leakageTolerance` enables rejected-mode leakage
%
% $$
% \ell_N=\max_{N<j\leq N_\mathrm{check}}
% \frac{\lVert\Pi_N^\mathrm{discrete}u_j\rVert_\mu}
% {\lVert u_j\rVert_\mu}.
% $$
%
% Its default check count is twice the candidate count. An explicit
% `nCheckModes` must exceed the candidate band so every prefix has at least
% one rejected comparison mode. `quadraticAliasingTolerance` enables
%
% $$
% q_N=\max_{1\leq i\leq j\leq N}
% \frac{\left\lVert\Pi_N^\mathrm{discrete}(u_i u_j)-
% \Pi_N^\mathrm{continuous}(u_i u_j)\right\rVert_\mu}
% {\lVert u_i u_j\rVert_\mu}.
% $$
%
% The continuous quadratic projection is evaluated on the source solver's
% inner-product grid, independently of the fixed transform rule. Leakage
% and quadratic aliasing require a positive-definite target Gram matrix;
% signed targets support Gram assessment only. `IMDiscreteTransformAssessment`
% defines every table column and policy-result field.
%
% ```matlab
% [transform,assessment] = basisSet.discreteTransform(nPoints=32);
% [transform,assessment] = basisSet.discreteTransform(z=z,nModes=8);
% [transform,assessment] = basisSet.discreteTransform(z=z,weights=w,nModes=8);
% assessment.prefixDiagnostics
% ```
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = discreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested point count for a mode-root grid
% - Parameter options.z: increasing explicit physical sample points
% - Parameter options.weights: optional fixed quadrature weights aligned with `z`
% - Parameter options.nModes: optional strict number of leading candidate modes for explicit `z`
% - Parameter options.gramTolerance: nonnegative normalized-Gram operator-error tolerance
% - Parameter options.leakageTolerance: optional positive rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional positive scalar quadratic-aliasing tolerance
% - Parameter options.nCheckModes: optional number of source modes used by the leakage policy
% - Returns transform: retained production transform
% - Returns assessment: transform diagnostics and, for certified construction, grid and search provenance
arguments
    self IMBasisSet
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.weights (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nModes double {mustBeReal, mustBeFinite} = zeros(0,1)
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
    [transform,assessment] = self.certifiedDiscreteTransform(nPoints=options.nPoints,z=options.z, ...
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
    error("IMBasisSet:InvalidDiscreteModeCount", "The basis set contains %d modes, but nModes=%d was requested.", nAvailableModes, candidateModeCount);
end

if hasExplicitWeights
    weightFit = [];
    candidateTransform = self.buildDiscreteTransform(z,options.weights(:),candidateModeCount);
else
    [~,weightFit] = self.quadratureWeightsForPoints(z=z,nModes=candidateModeCount);
    candidateTransform = weightFit.transform;
end
if (~isempty(options.leakageTolerance) || ~isempty(options.quadraticAliasingTolerance)) ...
        && ~candidateTransform.targetGramIsPositiveDefinite
    error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
        "Rejected-mode leakage and quadratic-aliasing policies require a positive-definite target Gram matrix; use the Gram policy alone for signed or indefinite targets.");
end

prefixTransforms = cell(candidateModeCount,1);
gramError = zeros(candidateModeCount,1);
roundTripError = zeros(candidateModeCount,1);
inverseConditionNumber = zeros(candidateModeCount,1);
gramConditionNumber = zeros(candidateModeCount,1);
sampledGramRank = zeros(candidateModeCount,1);
for nPrefix = 1:candidateModeCount
    prefixTransforms{nPrefix} = prefixTransform(candidateTransform,nPrefix);
    gramError(nPrefix) = prefixTransforms{nPrefix}.relativeGramOperatorError;
    roundTripError(nPrefix) = prefixTransforms{nPrefix}.roundTripError;
    inverseConditionNumber(nPrefix) = prefixTransforms{nPrefix}.inverseMatrixConditionNumber;
    gramConditionNumber(nPrefix) = prefixTransforms{nPrefix}.gramConditionNumber;
    sampledGramRank(nPrefix) = prefixTransforms{nPrefix}.sampledGramRank;
end
gramAccepted = IMDiscreteTransformTools.cumulativeAcceptance(gramError <= options.gramTolerance);
gramPolicy = IMDiscreteTransformTools.policyResult("gram",true,options.gramTolerance,gramError,gramAccepted,candidateModeCount);

leakageError = nan(candidateModeCount,1);
leakageLimitingModeNumber = nan(candidateModeCount,1);
if isempty(options.leakageTolerance)
    leakageAccepted = true(candidateModeCount,1);
    leakagePolicy = IMDiscreteTransformTools.policyResult("leakage",false,[],leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = [];
    leakagePolicy.limitingModeNumber = leakageLimitingModeNumber;
else
    if isempty(options.nCheckModes)
        nCheckModes = 2*candidateModeCount;
    else
        nCheckModes = options.nCheckModes;
    end
    if nCheckModes <= candidateModeCount
        error("IMBasisSet:InvalidLeakageCheckModeCount", "nCheckModes must be greater than the candidate mode count when leakage assessment is enabled.");
    end
    checkBasis = self;
    if nCheckModes > nAvailableModes
        try
            checkBasis = self.solver.solveEVP(self.evp,nModes=nCheckModes);
        catch cause
            exception = MException("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy could not obtain %d source modes.", nCheckModes);
            throw(addCause(exception,cause))
        end
        if size(checkBasis.nativeModes,2) < nCheckModes
            error("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy requested %d source modes, but the solver returned %d.", nCheckModes, size(checkBasis.nativeModes,2));
        end
        referenceEigenvalues = self.eigenvalues(1:candidateModeCount);
        checkEigenvalues = checkBasis.eigenvalues(1:candidateModeCount);
        eigenvalueScale = max(1,max(abs([referenceEigenvalues checkEigenvalues]),[],"all"));
        if ~isequal(checkBasis.modeNumber(1:candidateModeCount),self.modeNumber(1:candidateModeCount)) ...
                || any(abs(checkEigenvalues-referenceEigenvalues) > 1e-8*eigenvalueScale)
            error("IMBasisSet:AuxiliaryModeMismatch", "The leakage-policy solve did not reproduce the candidate modal prefix.");
        end
    end
    checkValues = checkBasis.u(candidateTransform.z);
    checkValues = checkValues(:,1:nCheckModes);
    checkGram = checkBasis.gramMatrix();
    checkNorms = diag(checkGram(1:nCheckModes,1:nCheckModes));
    if any(checkNorms <= 0)
        error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
            "Rejected-mode leakage requires positive continuous norms for every checked source mode.");
    end
    checkModeNumber = checkBasis.modeNumber(1:nCheckModes);
    [leakageError,leakageLimitingModeNumber] = rejectedModeLeakage(prefixTransforms,checkValues,checkNorms,checkModeNumber);
    leakageAccepted = IMDiscreteTransformTools.cumulativeAcceptance(leakageError <= options.leakageTolerance);
    leakagePolicy = IMDiscreteTransformTools.policyResult("leakage",true,options.leakageTolerance,leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = nCheckModes;
    leakagePolicy.limitingModeNumber = leakageLimitingModeNumber;
end

quadraticError = nan(candidateModeCount,1);
quadraticLimitingModeNumberI = nan(candidateModeCount,1);
quadraticLimitingModeNumberJ = nan(candidateModeCount,1);
if isempty(options.quadraticAliasingTolerance)
    quadraticAccepted = true(candidateModeCount,1);
    quadraticPolicy = IMDiscreteTransformTools.policyResult("quadraticAliasing",false,[],quadraticError,quadraticAccepted,candidateModeCount);
    quadraticPolicy.limitingModeNumberI = quadraticLimitingModeNumberI;
    quadraticPolicy.limitingModeNumberJ = quadraticLimitingModeNumberJ;
else
    integrationGrid = self.solver.innerProductGrid(self.zDomain);
    integrationModes = self.u(integrationGrid);
    integrationModes = integrationModes(:,1:candidateModeCount);
    context = self.evp.contextForSolver(self.solver);
    spec = self.evp.innerProduct();
    interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,integrationGrid,context);
    if isscalar(interiorWeight)
        interiorWeight = interiorWeight*ones(size(integrationGrid));
    else
        interiorWeight = interiorWeight(:);
    end
    if numel(interiorWeight) ~= length(integrationGrid) || ~isreal(interiorWeight) || any(~isfinite(interiorWeight))
        error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must return one finite real value for each inner-product point.");
    end
    endpointModeValues = zeros(0,candidateModeCount);
    endpointCoefficients = zeros(0,1);
    endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
    for iWeight = 1:numel(endpointWeights)
        endpointWeight = endpointWeights(iWeight);
        if endpointWeight.d ~= 0
            error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Scalar quadratic-aliasing assessment requires endpoint terms that depend only on values.");
        end
        if endpointWeight.c == 0
            continue
        end
        if string(endpointWeight.location) == "surface"
            zEndpoint = self.zDomain(2);
        else
            zEndpoint = self.zDomain(1);
        end
        valuesAtEndpoint = self.u(zEndpoint);
        endpointModeValues(end+1,:) = valuesAtEndpoint(1,1:candidateModeCount); %#ok<AGROW>
        endpointCoefficients(end+1,1) = endpointWeight.coefficient*endpointWeight.c*endpointWeight.c; %#ok<AGROW>
    end
    integrate = @(integrand) self.solver.integrateInnerProduct(integrationGrid,integrand,self.zDomain);
    [quadraticError,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ] = scalarQuadraticAliasing( ...
        candidateTransform,prefixTransforms,integrationModes,interiorWeight,integrate,endpointModeValues,endpointCoefficients);
    quadraticAccepted = IMDiscreteTransformTools.cumulativeAcceptance(quadraticError <= options.quadraticAliasingTolerance);
    quadraticPolicy = IMDiscreteTransformTools.policyResult("quadraticAliasing",true,options.quadraticAliasingTolerance,quadraticError,quadraticAccepted,candidateModeCount);
    quadraticPolicy.limitingModeNumberI = quadraticLimitingModeNumberI;
    quadraticPolicy.limitingModeNumberJ = quadraticLimitingModeNumberJ;
end

combinedAccepted = gramAccepted & leakageAccepted & quadraticAccepted;
retainedModeCount = find(combinedAccepted,1,"last");
if isempty(retainedModeCount)
    if hasExplicitModeCount && ~options.allowRetainedPrefix
        throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
    end
    error("IMBasisSet:NoAcceptableDiscreteTransformPrefix", "No candidate mode passes every enabled retained-band policy on this fixed quadrature rule.");
end
if hasExplicitModeCount && ~options.allowRetainedPrefix && retainedModeCount < candidateModeCount
    throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
end

[limitingPolicy,retentionReason] = combinedPolicyResult(candidateModeCount,retainedModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
transform = prefixTransforms{retainedModeCount};
modeCount = (1:candidateModeCount).';
lastModeNumber = candidateTransform.modeNumber(:);
prefixDiagnostics = table(modeCount,lastModeNumber,gramError,roundTripError,inverseConditionNumber,gramConditionNumber,sampledGramRank, ...
    leakageError,leakageLimitingModeNumber,quadraticError,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ, ...
    gramAccepted,leakageAccepted,quadraticAccepted,combinedAccepted, ...
    VariableNames=["modeCount","lastModeNumber","gramError","roundTripError","inverseMatrixConditionNumber","gramConditionNumber","sampledGramRank", ...
    "leakageError","leakageLimitingModeNumber","quadraticAliasingError","quadraticLimitingModeNumberI","quadraticLimitingModeNumberJ", ...
    "gramAccepted","leakageAccepted","quadraticAccepted","combinedAccepted"]);
assessment = IMDiscreteTransformAssessment(transform=transform,candidateTransform=candidateTransform,weightFit=weightFit, ...
    requestedPointCount=requestedPointCount,prefixDiagnostics=prefixDiagnostics,gramPolicy=gramPolicy, ...
    leakagePolicy=leakagePolicy,quadraticAliasingPolicy=quadraticPolicy,limitingPolicy=limitingPolicy,retentionReason=retentionReason,gridDesign=gridDesign);
end

function transform = prefixTransform(candidateTransform,nModes)
transform = IMDiscreteTransform(z=candidateTransform.z,weights=candidateTransform.weights, ...
    modeNumber=candidateTransform.modeNumber(1:nModes),normalization=candidateTransform.normalization, ...
    inverseMatrix=candidateTransform.inverseMatrix(:,1:nModes),metricMatrix=candidateTransform.metricMatrix, ...
    targetGramMatrix=candidateTransform.targetGramMatrix(1:nModes,1:nModes));
end

function [errors,limitingModeNumber] = rejectedModeLeakage(prefixTransforms,checkValues,checkNorms,checkModeNumber)
nCandidate = length(prefixTransforms);
nCheckModes = length(checkModeNumber);
errors = zeros(nCandidate,1);
limitingModeNumber = nan(nCandidate,1);
for nPrefix = 1:nCandidate
    transform = prefixTransforms{nPrefix};
    rejectedIndices = (nPrefix+1):nCheckModes;
    coefficients = transform.forwardMatrix*checkValues(:,rejectedIndices);
    targetNorms = abs(diag(transform.targetGramMatrix));
    numeratorSquared = sum(targetNorms.*abs(coefficients).^2,1);
    denominatorSquared = abs(checkNorms(rejectedIndices)).';
    leakage = sqrt(numeratorSquared./denominatorSquared);
    leakage(~isfinite(leakage)) = Inf;
    [errors(nPrefix),iLimiting] = max(leakage);
    limitingModeNumber(nPrefix) = checkModeNumber(rejectedIndices(iLimiting));
end
end

function [errors,limitingModeNumberI,limitingModeNumberJ] = scalarQuadraticAliasing(candidateTransform,prefixTransforms, ...
        integrationModes,interiorWeight,integrate,endpointModeValues,endpointCoefficients)
nCandidate = length(prefixTransforms);
nPairs = nCandidate*(nCandidate+1)/2;
pairIndices = zeros(nPairs,2);
sampleProducts = zeros(length(candidateTransform.z),nPairs);
continuousPairings = zeros(nCandidate,nPairs);
productNorms = zeros(nPairs,1);
iPair = 0;
for iMode = 1:nCandidate
    for jMode = iMode:nCandidate
        iPair = iPair + 1;
        pairIndices(iPair,:) = [iMode jMode];
        sampleProducts(:,iPair) = candidateTransform.inverseMatrix(:,iMode).*candidateTransform.inverseMatrix(:,jMode);
        product = integrationModes(:,iMode).*integrationModes(:,jMode);
        productNorm = integrate(interiorWeight.*product.*product);
        for iEndpoint = 1:length(endpointCoefficients)
            productAtEndpoint = endpointModeValues(iEndpoint,iMode)*endpointModeValues(iEndpoint,jMode);
            productNorm = productNorm + endpointCoefficients(iEndpoint)*productAtEndpoint*productAtEndpoint;
        end
        if productNorm <= 100*eps(max(1,abs(productNorm)))
            error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
                "Scalar quadratic-aliasing assessment requires a positive continuous norm for every retained-mode product.");
        end
        productNorms(iPair) = productNorm;
        for kMode = 1:nCandidate
            pairing = integrate(interiorWeight.*integrationModes(:,kMode).*product);
            for iEndpoint = 1:length(endpointCoefficients)
                productAtEndpoint = endpointModeValues(iEndpoint,iMode)*endpointModeValues(iEndpoint,jMode);
                pairing = pairing + endpointCoefficients(iEndpoint)*endpointModeValues(iEndpoint,kMode)*productAtEndpoint;
            end
            continuousPairings(kMode,iPair) = pairing;
        end
    end
end

errors = zeros(nCandidate,1);
limitingModeNumberI = nan(nCandidate,1);
limitingModeNumberJ = nan(nCandidate,1);
for nPrefix = 1:nCandidate
    transform = prefixTransforms{nPrefix};
    includedPairs = find(pairIndices(:,1) <= nPrefix & pairIndices(:,2) <= nPrefix);
    targetNorms = diag(transform.targetGramMatrix);
    discreteCoefficients = transform.forwardMatrix*sampleProducts(:,includedPairs);
    continuousCoefficients = continuousPairings(1:nPrefix,includedPairs)./targetNorms;
    coefficientDifference = discreteCoefficients-continuousCoefficients;
    numeratorSquared = sum(abs(targetNorms).*abs(coefficientDifference).^2,1);
    denominatorSquared = productNorms(includedPairs).';
    pairErrors = sqrt(numeratorSquared./denominatorSquared);
    pairErrors(~isfinite(pairErrors)) = Inf;
    [errors(nPrefix),iLimiting] = max(pairErrors);
    limitingPair = pairIndices(includedPairs(iLimiting),:);
    limitingModeNumberI(nPrefix) = candidateTransform.modeNumber(limitingPair(1));
    limitingModeNumberJ(nPrefix) = candidateTransform.modeNumber(limitingPair(2));
end
end

function throwStrictModeFailure(requestedModeCount,gramPolicy,leakagePolicy,quadraticPolicy)
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
failureDetails = strings(0,1);
for iPolicy = 1:length(policies)
    policy = policies{iPolicy};
    if policy.enabled && policy.acceptedModeCount < requestedModeCount
        failureDetails(end+1,1) = sprintf("%s (accepted %d, first failing value %.6g, tolerance %.6g)", ...
            policy.name,policy.acceptedModeCount,policy.limitingValue,policy.tolerance); %#ok<AGROW>
    end
end
error("IMBasisSet:DiscreteTransformPolicyFailed", ...
    "The explicitly requested %d-mode band failed retained-band policy: %s.", requestedModeCount, join(failureDetails,"; "));
end

function [limitingPolicy,reason] = combinedPolicyResult(candidateModeCount,retainedModeCount,gramPolicy,leakagePolicy,quadraticPolicy)
if retainedModeCount == candidateModeCount
    limitingPolicy = "none";
    reason = sprintf("All %d candidate modes pass every enabled retained-band policy.", candidateModeCount);
    return
end
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
limitingNames = strings(0,1);
for iPolicy = 1:length(policies)
    policy = policies{iPolicy};
    if policy.enabled && policy.acceptedModeCount == retainedModeCount
        limitingNames(end+1,1) = policy.name; %#ok<AGROW>
    end
end
limitingPolicy = join(limitingNames," and ");
reason = sprintf("Retained %d of %d candidate modes; the %s policy rejected the next prefix.", ...
    retainedModeCount,candidateModeCount,limitingPolicy);
end
