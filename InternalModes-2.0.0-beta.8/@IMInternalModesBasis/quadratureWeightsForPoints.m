function [weights, weightFit] = quadratureWeightsForPoints(self, options)
% Fit one quadrature rule to aligned internal-mode F/G channels.
%
% Every requested directly representable channel contributes its compressed
% normalized-Gram Frobenius system. The systems are stacked without an
% additional channel weight and solved for one nonnegative, full-depth
% constrained weight vector by default. Identically zero variable columns
% are omitted from that variable's objective while their family positions
% remain present in the returned transform.
%
% A custom objective receives the stacked system in
% `context.normalizedGramA` and `context.normalizedGramB`, row provenance in
% `context.normalizedGramVariables`, and scalar-style per-variable contexts
% in `context.variableContexts.F` and `.G`.
%
% - Topic: Build discrete transforms
% - Declaration: [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
% - Parameter options.z: increasing fixed physical sample points
% - Parameter options.nModes: number of leading aligned family modes
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.objective: built-in stacked objective or custom callback
% - Parameter options.nonnegative: impose nonnegative shared weights
% - Parameter options.constrainDepth: impose full-depth shared weights
% - Returns weights: fitted shared quadrature weights
% - Returns weightFit: specialized shared-fit diagnostics
arguments
    self IMInternalModesBasis
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
    options.variables (1,:) string = strings(1,0)
    options.objective = "stackedNormalizedGramFrobenius"
    options.nonnegative (1,1) logical = true
    options.constrainDepth (1,1) logical = true
end

z = options.z(:);
nModes = options.nModes;
depth = diff(self.zDomain);
edges = [self.zDomain(1);0.5*(z(1:end-1)+z(2:end));self.zDomain(2)];
geometricWeights = diff(edges);
if isempty(options.variables)
    variables = self.directlyRepresentableDiscreteVariables(z);
else
    variables = self.canonicalDiscreteVariables(options.variables);
end
if isempty(variables)
    error("IMInternalModesBasis:NoAvailableDiscreteTransformVariable", "Neither F nor G has a directly representable sampled metric on the supplied points.");
end
geometricTransform = self.buildInternalModesDiscreteTransform(z,geometricWeights,nModes,variables);

normalizedGramA = zeros(0,length(z));
normalizedGramB = zeros(0,1);
normalizedGramVariables = strings(0,1);
normalizedGramModePairs = zeros(0,2);
variableContexts = struct();
for variable = variables
    variableContext = geometricTransform.quadratureFitContext(variable=variable);
    variableA = variableContext.normalizedGramA;
    variableB = variableContext.normalizedGramB;
    variablePairs = variableContext.normalizedGramModePairs;
    rows = size(variableA,1);
    normalizedGramA = [normalizedGramA;variableA]; %#ok<AGROW>
    normalizedGramB = [normalizedGramB;variableB]; %#ok<AGROW>
    normalizedGramVariables = [normalizedGramVariables;repmat(variable,rows,1)]; %#ok<AGROW>
    normalizedGramModePairs = [normalizedGramModePairs;variablePairs]; %#ok<AGROW>
    variableContexts.(char(variable)) = variableContext;
end

objectiveContext = struct(z=z,modeNumber=geometricTransform.modeNumber,normalization=geometricTransform.normalization, ...
    geometricWeights=geometricWeights,availableVariables=variables,normalizedGramA=normalizedGramA,normalizedGramB=normalizedGramB, ...
    normalizedGramVariables=normalizedGramVariables,normalizedGramModePairs=normalizedGramModePairs,variableContexts=variableContexts);

objective = options.objective;
objectiveRowVariables = normalizedGramVariables;
objectiveModePairs = normalizedGramModePairs;
if ischar(objective) || (isstring(objective) && isscalar(objective))
    objectiveName = string(objective);
    if ~ismember(objectiveName,["stackedNormalizedGramFrobenius","normalizedGramFrobenius"])
        error("IMBasisSet:UnknownQuadratureObjective", "Unknown internal-mode quadrature objective ""%s"".", objectiveName);
    end
    objectiveName = "stackedNormalizedGramFrobenius";
    objectiveMatrix = normalizedGramA;
    objectiveTarget = normalizedGramB;
elseif isa(objective,"function_handle")
    try
        specification = objective(objectiveContext);
    catch cause
        exception = MException("IMBasisSet:QuadratureObjectiveFailed", "The custom internal-mode quadrature objective failed while building its least-squares system.");
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
    if isfield(specification,"rowVariables")
        objectiveRowVariables = string(specification.rowVariables(:));
    elseif size(objectiveMatrix,1) ~= length(normalizedGramVariables)
        nDefaultRows = length(normalizedGramVariables);
        nCustomRows = size(objectiveMatrix,1) - nDefaultRows;
        if nCustomRows >= 0
            objectiveRowVariables = [normalizedGramVariables;repmat("custom",nCustomRows,1)];
        else
            objectiveRowVariables = repmat("custom",size(objectiveMatrix,1),1);
        end
    end
    if isfield(specification,"modePairs")
        objectiveModePairs = specification.modePairs;
    elseif size(objectiveMatrix,1) ~= size(normalizedGramModePairs,1)
        nDefaultRows = size(normalizedGramModePairs,1);
        nCustomRows = size(objectiveMatrix,1) - nDefaultRows;
        if nCustomRows >= 0
            objectiveModePairs = [normalizedGramModePairs;zeros(nCustomRows,2)];
        else
            objectiveModePairs = zeros(size(objectiveMatrix,1),2);
        end
    end
else
    error("IMBasisSet:InvalidQuadratureObjective", "objective must name the stacked normalized-Gram objective or be a function handle.");
end

nSamples = length(z);
[objectiveMatrix,objectiveTarget] = IMDiscreteTransformTools.validateObjectiveSystem(objectiveMatrix,objectiveTarget,nSamples,true);
if length(objectiveRowVariables) ~= size(objectiveMatrix,1) || size(objectiveModePairs,1) ~= size(objectiveMatrix,1)
    error("IMBasisSet:InvalidQuadratureObjective", "Custom objective row provenance must align with its rows.");
end
[weights,exitFlag,solverOutput] = IMDiscreteTransformTools.fitQuadrature(objectiveMatrix,objectiveTarget,geometricWeights,depth,options.nonnegative,options.constrainDepth);

if nargout > 1
    transform = geometricTransform.transformWithWeights(weights);
    weightFit = IMInternalModesQuadratureWeightFit(transform=transform,geometricTransform=geometricTransform,objectiveName=objectiveName, ...
        objectiveMatrix=objectiveMatrix,objectiveTarget=objectiveTarget,objectiveRowVariables=objectiveRowVariables,objectiveModePairs=objectiveModePairs, ...
        nonnegativeConstraint=options.nonnegative,depthConstraint=options.constrainDepth,depthTarget=depth,exitFlag=exitFlag,solverOutput=solverOutput);
end
end
