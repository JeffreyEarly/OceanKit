classdef IMScalarDiscreteTransform < IMDiscreteTransform
    % Attach scalar-grid geometry and physical projection data to recovery.
    %
    % The generic superclass owns synthesis and recovery. This adapter adds
    % scalar samples, quadrature and the genuine physical IMProjection record.
    % Galerkin Gram evidence and prescribed coefficient systems remain
    % distinct from numerical sample rank. Check projection's assessment
    % capability flags before interpreting its coefficient system as a Gram. Failed diagnostic preparation retains synthesis/physical
    % pairings but does not expose an applicable recovery matrix.
    %
    % ```matlab
    % transform = basisSet.discreteTransform(z=z,weights=w,nModes=8);
    % physicalCoefficients = transform.projection.projectPhysical(values);
    % recoveredCoefficients = transform.transformForward(values);
    % ```
    %
    % - Topic: Create scalar transforms
    % - Topic: Inspect scalar transforms
    % - Declaration: classdef IMScalarDiscreteTransform < IMDiscreteTransform

    properties (SetAccess = private)
        % Physical sample points, one per scalar observation.
        % - Topic: Inspect scalar transforms
        z (:,1) double
        % Physical quadrature weights, distinct from observation-noise weights.
        % - Topic: Inspect scalar transforms
        weights (:,1) double
        % Optional legacy integer mode numbers; columnLabels carry scientific identity.
        % - Topic: Inspect scalar transforms
        modeNumber (1,:) double
        % Continuous normalization captured by the scalar family.
        % - Topic: Inspect scalar transforms
        normalization (1,1) string
        % Physical pairing, continuous forms and independent Gram diagnostics.
        % - Topic: Inspect scalar transforms
        projection (1,1) IMProjection
    end

    properties (Dependent, SetAccess = private)
        % Sampled scalar synthesis, also available as synthesisMatrix.
        % - Topic: Inspect scalar transforms
        inverseMatrix
        % Genuine scalar physical sample metric.
        % - Topic: Inspect scalar transforms
        metricMatrix
        % Sampled physical Gram or explicitly prescribed coefficient system.
        %
        % Check projection.supportsGramAssessment. A prescribed DD identity
        % system does not measure a sampled APV Gram or observation rank.
        % - Topic: Inspect scalar transforms
        gramMatrix
        % Continuous physical Gram target or prescribed coefficient target.
        % - Topic: Inspect scalar transforms
        targetGramMatrix
        % Physical Gram discrepancy, NaN without Gram assessment capability.
        % - Topic: Inspect scalar transforms
        relativeGramOperatorError
        % Numerical recovery round trip, infinity when recovery failed.
        % - Topic: Inspect scalar transforms
        roundTripError
        % Physical sampled-basis condition, NaN for a prescribed source dual.
        % - Topic: Inspect scalar transforms
        inverseMatrixConditionNumber
        % Condition number of the sampled Gram or prescribed coefficient system.
        % - Topic: Inspect scalar transforms
        gramConditionNumber
        % Rank of the projection system, independent of sampled observability.
        % - Topic: Inspect scalar transforms
        sampledGramRank
        % Whether the active physical or prescribed coefficient target is positive.
        % - Topic: Inspect scalar transforms
        targetGramIsPositiveDefinite
        % Whether any scalar physical quadrature weight is negative.
        % - Topic: Inspect scalar transforms
        hasNegativeWeights
    end

    methods
        function self = IMScalarDiscreteTransform(physical,recovery,options)
            % Bind separately prepared physical and numerical scalar operators.
            %
            % The shared engine prepares recovery once. This adapter validates
            % its alignment with the physical record without another solve or
            % factorization. General symmetric targets and string identities
            % are retained. A prescribed dual may lack physical Gram evidence
            % while still having independently prepared sampled recovery.
            %
            % - Topic: Create scalar transforms
            % - Parameter physical: genuine IMProjection, optionally physical-only
            % - Parameter recovery: shared-engine IMDiscreteTransform with aligned synthesis and physical dual
            % - Parameter options.z: scalar sample column
            % - Parameter options.weights: unchanged physical quadrature weights
            % - Parameter options.normalization: continuous normalization name
            % - Parameter options.modeNumber: optional aligned legacy integer labels
            % - Returns self: scalar adapter retaining physical projection and numerical recovery
            arguments (Input)
                physical (1,1) IMProjection
                recovery (1,1) IMDiscreteTransform
                options.z (:,1) double {mustBeReal,mustBeFinite}
                options.weights (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string
                options.modeNumber (1,:) double {mustBeReal,mustBeFinite,mustBeInteger} = []
            end
            nSamples = numel(options.z);
            if numel(options.weights) ~= nSamples || recovery.sampleCount ~= nSamples || physical.sampleCount ~= nSamples
                error("IMScalarDiscreteTransform:InvalidShape","Samples, weights, physical pairing and recovery must describe the same observations.");
            end
            if ~isequal(physical.columnLabels,recovery.columnLabels) || ~isequal(physical.activeColumnMask,recovery.activeColumnMask) ...
                    || ~isequal(physical.physicalDualMatrix,recovery.physicalDualMatrix)
                error("IMScalarDiscreteTransform:InvalidRecovery","Recovery must retain the physical record's labels, scientific activity and unchanged physical dual.");
            end
            if physical.projectionKind == "galerkin" && ~isequal(physical.sampledBasis,recovery.synthesisMatrix)
                error("IMScalarDiscreteTransform:InvalidRecovery","Scalar synthesis must exactly match the physical Galerkin basis.");
            end
            if ~isempty(options.modeNumber) && numel(options.modeNumber) ~= recovery.columnCount
                error("IMScalarDiscreteTransform:InvalidShape","Optional modeNumber must align with the retained columns.");
            end
            if (isfield(recovery.metadata,"z") && ~isequal(recovery.metadata.z,options.z)) ...
                    || (isfield(recovery.metadata,"normalization") && ~isequal(string(recovery.metadata.normalization),options.normalization))
                error("IMScalarDiscreteTransform:InvalidRecovery","Prepared recovery metadata must match the sample order and normalization.");
            end
            self@IMDiscreteTransform(recovery);
            self.z = options.z;
            self.weights = options.weights;
            self.modeNumber = options.modeNumber;
            self.normalization = options.normalization;
            self.projection = physical;
        end

        function value = get.inverseMatrix(self), value = self.synthesisMatrix; end
        function value = get.metricMatrix(self), value = self.projection.metricMatrix; end
        function value = get.gramMatrix(self), value = self.projection.gramMatrix; end
        function value = get.targetGramMatrix(self), value = self.projection.targetGramMatrix; end
        function value = get.relativeGramOperatorError(self), value = self.projection.gramError; end
        function value = get.roundTripError(self)
            if self.success, value = self.diagnostics.roundTripError; else, value = Inf; end
        end
        function value = get.inverseMatrixConditionNumber(self), value = self.projection.inverseMatrixConditionNumber; end
        function value = get.gramConditionNumber(self), value = self.projection.gramConditionNumber; end
        function value = get.sampledGramRank(self), value = self.projection.sampledGramRank; end
        function value = get.targetGramIsPositiveDefinite(self), value = self.projection.targetGramIsPositiveDefinite; end
        function value = get.hasNegativeWeights(self), value = any(self.weights < 0); end
    end
end
