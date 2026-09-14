classdef (Hidden, Sealed) IMDiscreteTransformTools
    % Share retained-band and quadrature operations across transform builders.

    methods (Static)
        function validateOptionalPositiveInteger(value,name)
            if ~isempty(value) && (~isscalar(value) || value <= 0 || value ~= floor(value))
                error("IMBasisSet:InvalidDiscreteTransformOption", "%s must be empty or a positive integer scalar.",name);
            end
        end

        function validateOptionalPositiveTolerance(value,name)
            if ~isempty(value) && (~isscalar(value) || value <= 0)
                error("IMBasisSet:InvalidDiscreteTransformOption", "%s must be empty or a positive finite scalar.",name);
            end
        end

        function [z,candidateModeCount,gridDesign] = pointsForExactCount(basisSet,nPoints,nAvailableModes)
            pointCounts = nan(nAvailableModes,1);
            pointGrids = cell(nAvailableModes,1);
            pointGridDesigns = cell(nAvailableModes,1);
            maximumPossibleModeCount = min(nAvailableModes,nPoints);
            for nModes = maximumPossibleModeCount:-1:1
                try
                    [pointGrids{nModes},pointGridDesigns{nModes}] = basisSet.pointsFromModeRoots(nModes=nModes);
                    pointCounts(nModes) = length(pointGrids{nModes});
                catch exception
                    if nModes == nAvailableModes && strcmp(exception.identifier,"IMBasisSet:AuxiliaryModeUnavailable")
                        continue
                    end
                    rethrow(exception)
                end
                if pointCounts(nModes) == nPoints
                    candidateModeCount = nModes;
                    z = pointGrids{nModes};
                    gridDesign = pointGridDesigns{nModes};
                    return
                end
            end
            for nModes = (maximumPossibleModeCount+1):nAvailableModes
                try
                    [pointGrids{nModes},pointGridDesigns{nModes}] = basisSet.pointsFromModeRoots(nModes=nModes);
                    pointCounts(nModes) = length(pointGrids{nModes});
                catch exception
                    if nModes == nAvailableModes && strcmp(exception.identifier,"IMBasisSet:AuxiliaryModeUnavailable")
                        continue
                    end
                    rethrow(exception)
                end
            end
            attainableCounts = unique(pointCounts(isfinite(pointCounts)));
            [~,order] = sort(abs(attainableCounts-nPoints));
            nearestCounts = attainableCounts(order(1:min(3,length(order))));
            error("IMBasisSet:UnattainableDiscretePointCount", "No available mode-root grid has exactly %d points. Nearest attainable counts are %s.",nPoints,join(string(nearestCounts(:).'),", "));
        end

        function gridDesign = explicitGridDesign(basisSet,z)
            gridDesign = struct( ...
                "kind","explicit", ...
                "pointCount",length(z), ...
                "sourceClass",string(class(basisSet)), ...
                "sourceEVP",string(basisSet.evp.name), ...
                "sourceFamily",IMDiscreteTransformTools.sourceFamily(basisSet), ...
                "generatingVariable","", ...
                "generatingModeIndex",NaN, ...
                "generatingModeNumber",NaN, ...
                "representedModeCount",NaN, ...
                "interpretationForG","Caller-supplied physical points.");
        end

        function gridDesign = modeRootGridDesign(basisSet,z,nModes,generatingModeNumber)
            generatingVariable = "u";
            interpretationForG = "Not an internal-mode basis.";
            if isa(basisSet,"IMInternalModesBasis")
                generatingVariable = string(basisSet.evp.formulation);
                if generatingVariable == "F"
                    interpretationForG = "Roots of F; for nonzero eigenvalues these are G-extrema-like because Gz=lambda*F.";
                else
                    interpretationForG = "Roots of G.";
                end
            end
            gridDesign = struct( ...
                "kind","modeRoot", ...
                "pointCount",length(z), ...
                "sourceClass",string(class(basisSet)), ...
                "sourceEVP",string(basisSet.evp.name), ...
                "sourceFamily",IMDiscreteTransformTools.sourceFamily(basisSet), ...
                "generatingVariable",generatingVariable, ...
                "generatingModeIndex",nModes+1, ...
                "generatingModeNumber",generatingModeNumber, ...
                "representedModeCount",nModes, ...
                "interpretationForG",interpretationForG);
        end

        function gridDesign = validatedGridDesign(basisSet,z,gridDesign)
            if isempty(gridDesign)
                gridDesign = IMDiscreteTransformTools.explicitGridDesign(basisSet,z);
                return
            end
            requiredFields = ["kind","pointCount","sourceClass","sourceEVP","sourceFamily", ...
                "generatingVariable","generatingModeIndex","generatingModeNumber", ...
                "representedModeCount","interpretationForG"];
            if ~isscalar(gridDesign) || ~all(isfield(gridDesign,requiredFields))
                error("IMBasisSet:InvalidGridDesign", "gridDesign must be the scalar provenance struct returned by modeRootGrid.");
            end
            if ~isscalar(gridDesign.pointCount) || gridDesign.pointCount ~= length(z)
                error("IMBasisSet:InvalidGridDesign", "gridDesign.pointCount must equal length(z).");
            end
        end

        function accepted = cumulativeAcceptance(rawAccepted)
            accepted = cumprod(double(rawAccepted(:))) > 0;
        end

        function policy = policyResult(name,enabled,tolerance,errorValues,accepted,candidateModeCount)
            acceptedModeCount = find(accepted,1,"last");
            if isempty(acceptedModeCount)
                acceptedModeCount = 0;
            end
            if ~enabled
                limitingValue = NaN;
                reason = "The "+name+" policy is disabled.";
            elseif acceptedModeCount < candidateModeCount
                limitingValue = errorValues(acceptedModeCount+1);
                reason = sprintf("Accepted %d of %d candidate modes; the next prefix has value %.6g above tolerance %.6g.",acceptedModeCount,candidateModeCount,limitingValue,tolerance);
            else
                limitingValue = max(errorValues);
                reason = sprintf("All %d candidate modes pass tolerance %.6g; the largest value is %.6g.",candidateModeCount,tolerance,limitingValue);
            end
            policy = struct(name=string(name),enabled=logical(enabled),tolerance=tolerance,error=errorValues(:),accepted=accepted(:), ...
                acceptedModeCount=acceptedModeCount,maximumAcceptedModeCount=acceptedModeCount,limitingValue=limitingValue,reason=string(reason),limited=enabled && acceptedModeCount < candidateModeCount);
        end

        function [objectiveMatrix,objectiveTarget] = validateObjectiveSystem(objectiveMatrix,objectiveTarget,nSamples,allowEmpty)
            if allowEmpty && isempty(objectiveMatrix) && size(objectiveMatrix,2) == nSamples && isempty(objectiveTarget)
                objectiveMatrix = double(objectiveMatrix);
                objectiveTarget = zeros(0,1);
                return
            end
            matrixIsValid = isnumeric(objectiveMatrix) && ismatrix(objectiveMatrix) && isreal(objectiveMatrix) ...
                && size(objectiveMatrix,2) == nSamples && all(isfinite(objectiveMatrix),"all");
            if ~allowEmpty
                matrixIsValid = matrixIsValid && ~isempty(objectiveMatrix);
            end
            if ~matrixIsValid
                error("IMBasisSet:InvalidQuadratureObjective", "The objective matrix A must be finite, real, and have one column per sample point.");
            end
            if ~isnumeric(objectiveTarget) || ~isvector(objectiveTarget) || ~isreal(objectiveTarget) ...
                    || length(objectiveTarget) ~= size(objectiveMatrix,1) || any(~isfinite(objectiveTarget),"all")
                error("IMBasisSet:InvalidQuadratureObjective", "The objective target b must be finite, real, and have one entry per objective row.");
            end
            objectiveMatrix = double(objectiveMatrix);
            objectiveTarget = double(objectiveTarget(:));
        end

        function [weights,exitFlag,solverOutput] = fitQuadrature(objectiveMatrix,objectiveTarget,geometricWeights,depth,nonnegative,constrainDepth)
            nSamples = length(geometricWeights);
            if isempty(objectiveMatrix)
                weights = geometricWeights;
                exitFlag = 1;
                solverOutput = struct("message","No active modal Gram rows; returned the geometric rule satisfying the requested constraints.");
                return
            end
            if isempty(which("lsqlin"))
                error("IMBasisSet:MissingQuadratureOptimizer", "Fitting quadrature weights requires lsqlin from Optimization Toolbox. Supply weights explicitly when it is unavailable.");
            end
            if nonnegative
                lowerBounds = zeros(nSamples,1);
            else
                lowerBounds = [];
            end
            if constrainDepth
                equalityMatrix = ones(1,nSamples);
                equalityTarget = 1;
            else
                equalityMatrix = [];
                equalityTarget = [];
            end
            try
                [scaledWeights,~,~,exitFlag,solverOutput] = lsqlin(depth*objectiveMatrix,objectiveTarget,[],[],equalityMatrix,equalityTarget,lowerBounds,[],geometricWeights/depth,IMDiscreteTransformTools.quadratureSolverOptions());
            catch cause
                exception = MException("IMBasisSet:QuadratureWeightFitFailed", "lsqlin failed while fitting quadrature weights.");
                throw(addCause(exception,cause))
            end
            if isempty(scaledWeights) || exitFlag <= 0 || any(~isfinite(scaledWeights))
                error("IMBasisSet:QuadratureWeightFitFailed", "lsqlin did not converge to finite quadrature weights (exit flag %d).",exitFlag);
            end
            weights = depth*scaledWeights;
            constraintTolerance = 1e-10*max(1,depth);
            if nonnegative && any(weights < -constraintTolerance)
                error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested nonnegative constraint.");
            end
            if constrainDepth && abs(sum(weights)-depth) > constraintTolerance
                error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested full-depth constraint.");
            end
            weights(weights < 0 & weights >= -constraintTolerance) = 0;
        end
    end

    methods (Static, Access = private)
        function family = sourceFamily(basisSet)
            family = "scalar";
            if isa(basisSet,"IMInternalModesBasis")
                family = string(basisSet.evp.modeFamily);
            end
        end

        function options = quadratureSolverOptions()
            persistent cachedOptions
            if isempty(cachedOptions)
                cachedOptions = optimoptions("lsqlin",Algorithm="active-set",Display="off",ConstraintTolerance=1e-12,OptimalityTolerance=1e-12,StepTolerance=1e-14,MaxIterations=5000);
            end
            options = cachedOptions;
        end
    end
end
