function [weights, weightFit] = quadratureWeightsForPoints(self, options)
% Find quadrature weights for fixed physical sample points.
%
% The points `z` are fixed. This method solves for one weight $$w_k$$ per
% point so the sampled inner products of the first `nModes` retained modes
% approximate their continuous Gram matrix. Let $$\Gamma(w)$$ be the
% sampled Gram matrix, let $$\Gamma_0$$ be its continuous target, and define
%
% $$
% S=\operatorname{diag}\!\left(
% \left|\operatorname{diag}\Gamma_0\right|^{-1/2}
% \right),
% \qquad
% E(w)=S\left(\Gamma(w)-\Gamma_0\right)S.
% $$
%
% The default `"normalizedGramFrobenius"` objective minimizes
% $$\|E(w)\|_{\mathrm F}$$. It combines the mismatch from every retained
% mode pair: diagonal entries of $$E$$ measure modal norm errors, while
% off-diagonal entries measure lost orthogonality. For target modal norms
% $$C_i=(\Gamma_0)_{ii}$$ and fixed endpoint contribution
% $$\Gamma_{\mathrm{endpoint}}$$, the corresponding least-squares system is
%
% $$
% (A_{\mathrm{LS}})_{(i,j),k}
% =\rho_{ij}\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
% \qquad
% (b_{\mathrm{LS}})_{(i,j)}
% =\rho_{ij}\frac{(\Gamma_0-\Gamma_{\mathrm{endpoint}})_{ij}}
% {\sqrt{|C_iC_j|}}.
% $$
%
% Only pairs with $$1\leq i\leq j\leq n_m$$ are stored, in row order
% $$(1,1),(1,2),\ldots,(1,n_m),(2,2),\ldots,(n_m,n_m)$$, with
%
% $$
% \rho_{ij}=
% \begin{cases}
% 1, & i=j,\\
% \sqrt{2}, & i<j.
% \end{cases}
% $$
%
% Because the normalized Gram mismatch $$E(w)$$ is symmetric,
%
% $$
% \|E(w)\|_{\mathrm F}^2
% =\sum_i E_{ii}(w)^2+2\sum_{i<j}E_{ij}(w)^2.
% $$
%
% The upper-triangle system therefore gives exactly the full Frobenius
% objective with $$n_m(n_m+1)/2$$ rows instead of $$n_m^2$$ rows.
%
% The fitted weights solve
%
% $$
% \min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
% $$
%
% By default they also satisfy
%
% $$
% w_k\geq0,
% \qquad
% \sum_k w_k=z_\mathrm{surface}-z_\mathrm{bottom}.
% $$
%
% Geometric control-volume weights provide the initial guess and reference
% comparison. Unlike those geometric widths, fitted weights are algebraic
% quadrature coefficients and need not remain positive or sum to the full
% depth when the corresponding constraints are disabled. The optimizer
% works with dimensionless weights $$x_k=w_k/D$$, where $$D$$ is the full
% depth. Custom objectives continue to define $$Aw-b$$ in physical units.
% This method requires `lsqlin` from Optimization Toolbox.
%
% For the built-in objective, `weightFit.residualNorm` is the aggregate
% error $$\|E(w)\|_{\mathrm F}$$. The resulting transform separately reports
% `relativeGramOperatorError` as $$\|E(w)\|_2$$, the largest Gram distortion
% over any normalized combination of retained modes. `roundTripError`
% measures algebraic coefficient recovery and can be tiny even when either
% Gram error is appreciable. The quadrature-weight regression sweep supports
% retaining the unregularized Frobenius objective; geometric weights remain
% the initial guess and comparison baseline.
%
% A custom objective is a function handle accepting a context struct and
% returning a scalar struct with fields `A`, `b`, and optional `name`. The
% context contains `z`, `modeNumber`, `normalization`, `inverseMatrix`,
% `interiorWeight`, `targetGramMatrix`, `endpointGramMatrix`,
% `geometricWeights`, `normalizedGramA`, `normalizedGramB`, and
% `normalizedGramModePairs`. Row `q` of `normalizedGramModePairs` contains
% the retained basis-column indices `[iMode jMode]` represented by row `q`
% of the normalized Gram system. The diagonal or $$\sqrt{2}$$ row factor is
% already included in `normalizedGramA` and `normalizedGramB`.
%
% ```matlab
% [weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
% transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
% [weightFit.residualNorm weightFit.geometricResidualNorm]
% [transform.relativeGramOperatorError weightFit.geometricTransform.relativeGramOperatorError]
% ```
%
% - Topic: Build discrete transforms
% - Declaration: [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
% - Parameter options.z: increasing fixed physical sample points
% - Parameter options.nModes: number of leading retained modes
% - Parameter options.objective: `"normalizedGramFrobenius"` or a custom least-squares callback
% - Parameter options.nonnegative: whether fitted weights must be nonnegative
% - Parameter options.constrainDepth: whether fitted weights must sum to the full depth
% - Returns weights: fitted quadrature weights aligned with `z`
% - Returns weightFit: optional fitted-versus-geometric diagnostics
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
    options.objective = "normalizedGramFrobenius"
    options.nonnegative (1,1) logical = true
    options.constrainDepth (1,1) logical = true
end

z = options.z(:);
nModes = options.nModes;
depth = diff(self.zDomain);
edges = [self.zDomain(1); 0.5*(z(1:end-1) + z(2:end)); self.zDomain(2)];
geometricWeights = diff(edges);
geometricTransform = self.buildDiscreteTransform(z,geometricWeights,nModes);

evpContext = self.evp.contextForSolver(self.solver);
innerProduct = self.evp.innerProduct();
interiorWeight = IMEigenvalueProblem.evaluateCoefficient(innerProduct.interiorWeight,z,evpContext);
if isscalar(interiorWeight)
    interiorWeight = interiorWeight*ones(size(z));
else
    interiorWeight = interiorWeight(:);
end

inverseMatrix = geometricTransform.inverseMatrix;
targetGramMatrix = geometricTransform.targetGramMatrix;
interiorMetricMatrix = diag(interiorWeight.*geometricWeights);
endpointMetricMatrix = geometricTransform.metricMatrix - interiorMetricMatrix;
endpointGramMatrix = inverseMatrix.'*endpointMetricMatrix*inverseMatrix;
endpointGramMatrix = 0.5*(endpointGramMatrix + endpointGramMatrix.');
targetNorms = diag(targetGramMatrix);

nSamples = length(z);
nGramRows = nModes*(nModes + 1)/2;
normalizedGramA = zeros(nGramRows,nSamples);
normalizedGramB = zeros(nGramRows,1);
normalizedGramModePairs = zeros(nGramRows,2);
iRow = 0;
for iMode = 1:nModes
    for jMode = iMode:nModes
        iRow = iRow + 1;
        rowFactor = 1;
        if iMode ~= jMode
            rowFactor = sqrt(2);
        end
        scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
        normalizedGramA(iRow,:) = rowFactor*(interiorWeight.*inverseMatrix(:,iMode).*inverseMatrix(:,jMode)).'/scale;
        normalizedGramB(iRow) = rowFactor*(targetGramMatrix(iMode,jMode) - endpointGramMatrix(iMode,jMode))/scale;
        normalizedGramModePairs(iRow,:) = [iMode jMode];
    end
end

objectiveContext = struct();
objectiveContext.z = z;
objectiveContext.modeNumber = geometricTransform.modeNumber;
objectiveContext.normalization = geometricTransform.normalization;
objectiveContext.inverseMatrix = inverseMatrix;
objectiveContext.interiorWeight = interiorWeight;
objectiveContext.targetGramMatrix = targetGramMatrix;
objectiveContext.endpointGramMatrix = endpointGramMatrix;
objectiveContext.geometricWeights = geometricWeights;
objectiveContext.normalizedGramA = normalizedGramA;
objectiveContext.normalizedGramB = normalizedGramB;
objectiveContext.normalizedGramModePairs = normalizedGramModePairs;

objective = options.objective;
if ischar(objective) || (isstring(objective) && isscalar(objective))
    objectiveName = string(objective);
    if objectiveName ~= "normalizedGramFrobenius"
        error("IMBasisSet:UnknownQuadratureObjective", "Unknown quadrature objective ""%s"".", objectiveName);
    end
    objectiveMatrix = normalizedGramA;
    objectiveTarget = normalizedGramB;
elseif isa(objective,"function_handle")
    try
        specification = objective(objectiveContext);
    catch cause
        exception = MException("IMBasisSet:QuadratureObjectiveFailed", "The custom quadrature objective failed while building its least-squares system.");
        throw(addCause(exception,cause))
    end
    if ~isstruct(specification) || ~isscalar(specification) || ~isfield(specification,"A") || ~isfield(specification,"b")
        error("IMBasisSet:InvalidQuadratureObjective", "A custom quadrature objective must return a scalar struct with fields A and b.");
    end
    objectiveMatrix = specification.A;
    objectiveTarget = specification.b;
    objectiveName = "custom";
    if isfield(specification,"name")
        if ~(ischar(specification.name) || (isstring(specification.name) && isscalar(specification.name)))
            error("IMBasisSet:InvalidQuadratureObjective", "The optional custom objective name must be a text scalar.");
        end
        objectiveName = string(specification.name);
    end
else
    error("IMBasisSet:InvalidQuadratureObjective", "objective must be ""normalizedGramFrobenius"" or a function handle returning a least-squares system.");
end

[objectiveMatrix,objectiveTarget] = IMDiscreteTransformTools.validateObjectiveSystem(objectiveMatrix,objectiveTarget,nSamples,false);
[weights,exitFlag,solverOutput] = IMDiscreteTransformTools.fitQuadrature(objectiveMatrix,objectiveTarget,geometricWeights,depth,options.nonnegative,options.constrainDepth);

if nargout > 1
    transform = self.buildDiscreteTransform(z,weights,nModes);
    weightFit = IMQuadratureWeightFit(transform=transform,geometricTransform=geometricTransform,objectiveName=objectiveName,objectiveMatrix=objectiveMatrix,objectiveTarget=objectiveTarget,nonnegativeConstraint=options.nonnegative,depthConstraint=options.constrainDepth,depthTarget=depth,exitFlag=exitFlag,solverOutput=solverOutput);
end
end
