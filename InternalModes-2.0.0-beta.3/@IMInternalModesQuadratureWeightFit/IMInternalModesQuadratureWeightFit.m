classdef IMInternalModesQuadratureWeightFit
    % Store diagnostics for a shared F/G quadrature-weight fit.
    %
    % The fitted and geometric transforms use identical points, aligned
    % modal families, and requested variables.  The least-squares rows are
    % the unweighted stack of each variable's compressed normalized-Gram
    % Frobenius system. `objectiveRowVariables` records which physical
    % channel generated every row.
    % The built-in objective is
    %
    % $$\sum_{V\in\mathcal V}\left\lVert
    % S_V(\Gamma_V(w)-\Gamma_{0,V})S_V
    % \right\rVert_\mathrm{F}^{2}.$$
    %
    % `variableResidualNorm` and `variableGeometricResidualNorm` recover the
    % individual channel contributions to this stacked objective.
    %
    % - Topic: Inspect fitted weights
    % - Topic: Assess fit quality
    % - Topic: Inspect constraints
    % - Declaration: classdef IMInternalModesQuadratureWeightFit

    properties (SetAccess = private)
        % Specialized transform using fitted shared weights.
        transform
        % Fitted shared quadrature weights.
        weights
        % Specialized transform using geometric control-volume weights.
        geometricTransform
        % Geometric control-volume weights.
        geometricWeights
        % Direct channels included in the stacked fit.
        availableVariables
        % Name of the least-squares objective.
        objectiveName
        % Least-squares matrix with one column per sample.
        objectiveMatrix
        % Least-squares target vector.
        objectiveTarget
        % F/G/custom provenance for every objective row.
        objectiveRowVariables
        % Aligned family-column pair for every objective row; [0 0] means custom.
        objectiveModePairs
        % Fitted objective residual vector.
        residual
        % Geometric objective residual vector.
        geometricResidual
        % Norm of the fitted stacked residual.
        residualNorm
        % Norm of the geometric stacked residual.
        geometricResidualNorm
        % Whether nonnegative weights were imposed.
        nonnegativeConstraint
        % Whether exact full-depth weight sum was imposed.
        depthConstraint
        % Full physical depth targeted by the constraint.
        depthTarget
        % Fitted weight-sum error relative to depth.
        depthError
        % Geometric weight-sum error relative to depth.
        geometricDepthError
        % `lsqlin` exit flag.
        exitFlag
        % `lsqlin` optimizer diagnostics.
        solverOutput
    end

    methods
        function self = IMInternalModesQuadratureWeightFit(options)
            % Create diagnostics for one shared quadrature-weight fit.
            % - Topic: Developer topics — Construct fit results
            % - Declaration: weightFit = IMInternalModesQuadratureWeightFit(options)
            % - Developer: true
            arguments
                options.transform (1,1) IMInternalModesDiscreteTransform
                options.geometricTransform (1,1) IMInternalModesDiscreteTransform
                options.objectiveName {mustBeTextScalar}
                options.objectiveMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.objectiveTarget (:,1) double {mustBeReal, mustBeFinite}
                options.objectiveRowVariables (:,1) string
                options.objectiveModePairs (:,2) double {mustBeInteger, mustBeNonnegative}
                options.nonnegativeConstraint (1,1) logical
                options.depthConstraint (1,1) logical
                options.depthTarget (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.exitFlag (1,1) double {mustBeReal, mustBeFinite}
                options.solverOutput (1,1) struct
            end

            transform = options.transform;
            geometricTransform = options.geometricTransform;
            if ~isequal(transform.z,geometricTransform.z) || ~isequal(transform.modeNumber,geometricTransform.modeNumber) ...
                    || transform.normalization ~= geometricTransform.normalization || ~isequal(transform.availableVariables,geometricTransform.availableVariables)
                error("IMInternalModesQuadratureWeightFit:IncompatibleTransforms", "The fitted and geometric transforms must use the same points, modes, normalization, and variables.");
            end
            nRows = size(options.objectiveMatrix,1);
            if size(options.objectiveMatrix,2) ~= length(transform.z) || length(options.objectiveTarget) ~= nRows ...
                    || length(options.objectiveRowVariables) ~= nRows || size(options.objectiveModePairs,1) ~= nRows
                error("IMInternalModesQuadratureWeightFit:InvalidObjective", "Objective arrays and provenance must have compatible dimensions.");
            end
            if any(~ismember(options.objectiveRowVariables,["F","G","custom"]))
                error("IMInternalModesQuadratureWeightFit:InvalidObjective", "objectiveRowVariables may contain F, G, or custom.");
            end

            residual = options.objectiveMatrix*transform.weights-options.objectiveTarget;
            geometricResidual = options.objectiveMatrix*geometricTransform.weights-options.objectiveTarget;
            self.transform = transform;
            self.weights = transform.weights;
            self.geometricTransform = geometricTransform;
            self.geometricWeights = geometricTransform.weights;
            self.availableVariables = transform.availableVariables;
            self.objectiveName = string(options.objectiveName);
            self.objectiveMatrix = options.objectiveMatrix;
            self.objectiveTarget = options.objectiveTarget;
            self.objectiveRowVariables = options.objectiveRowVariables(:);
            self.objectiveModePairs = options.objectiveModePairs;
            self.residual = residual;
            self.geometricResidual = geometricResidual;
            self.residualNorm = norm(residual,2);
            self.geometricResidualNorm = norm(geometricResidual,2);
            self.nonnegativeConstraint = options.nonnegativeConstraint;
            self.depthConstraint = options.depthConstraint;
            self.depthTarget = options.depthTarget;
            self.depthError = sum(transform.weights)-options.depthTarget;
            self.geometricDepthError = sum(geometricTransform.weights)-options.depthTarget;
            self.exitFlag = options.exitFlag;
            self.solverOutput = options.solverOutput;
        end

        function value = variableResidualNorm(self, options)
            % Return the fitted normalized-Gram residual for one variable.
            arguments
                self IMInternalModesQuadratureWeightFit
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = "F"
            end
            rows = self.objectiveRowVariables == string(options.variable);
            value = norm(self.residual(rows),2);
        end

        function value = variableGeometricResidualNorm(self, options)
            % Return the geometric normalized-Gram residual for one variable.
            arguments
                self IMInternalModesQuadratureWeightFit
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = "F"
            end
            rows = self.objectiveRowVariables == string(options.variable);
            value = norm(self.geometricResidual(rows),2);
        end
    end
end
