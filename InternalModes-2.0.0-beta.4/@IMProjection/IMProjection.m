classdef IMProjection
    % Construct a fixed sampled projection and measure its numerical quality.
    %
    % The signed sampled pairing defines coefficients through
    % $$A = (\Phi^T W\Phi)^{-1}\Phi^T W.$$
    % Rank-deficient pairings use a pseudoinverse and report infinite Gram
    % error. The separately supplied positive majorant measures coefficient
    % errors; it never replaces the physical pairing. This object does not
    % choose a grid, fit weights, accept tolerances, or select model counts.
    % `fromPrescribedDual` instead accepts an explicit physical source dual and its
    % coefficient system. It has no synthesis basis and cannot measure a
    % sampled-basis Gram discrepancy or round trip.
    %
    % ```matlab
    % projection = IMProjection(sampledBasis,metricMatrix,targetGramMatrix);
    % coefficients = projection.project(values);
    % ```
    %
    % - Topic: Create projections
    % - Topic: Inspect projection data
    % - Topic: Measure projection quality
    % - Topic: Apply projections
    % - Declaration: classdef IMProjection

    properties (SetAccess = private)
        % Number of sampled source coordinates accepted by project.
        % - Topic: Inspect projection data
        sampleCount
        % Number of output coefficient coordinates.
        % - Topic: Inspect projection data
        columnCount
        % Galerkin construction or a caller-prescribed physical source dual.
        % - Topic: Inspect projection data
        projectionKind
        % Whether an actual sampled basis defines a Gram comparison.
        % - Topic: Inspect projection data
        supportsGramAssessment
        % Whether synthesis data defines an active-column round trip.
        % - Topic: Inspect projection data
        supportsRoundTripAssessment
        % Caller-supplied construction identity and numerical provenance.
        % - Topic: Inspect projection data
        provenance
        % Signed map from source sample values to coefficient pairings.
        % - Topic: Inspect projection data
        samplePairingMatrix
        % Sampled basis, or empty for a prescribed dual without synthesis.
        % - Topic: Inspect projection data
        sampledBasis
        % Signed observation metric, or empty for a prescribed dual.
        % - Topic: Inspect projection data
        metricMatrix
        % Continuous signed Gram matrix in the supplied column coordinates.
        % - Topic: Inspect projection data
        targetGramMatrix
        % Positive coefficient metric used for error magnitudes.
        % - Topic: Inspect projection data
        majorantGramMatrix
        % Logical row identifying columns with a direct projection.
        % - Topic: Inspect projection data
        activeColumnMask
        % Scientific column labels, distinct from array indices.
        % - Topic: Inspect projection data
        columnLabels
        % Matrix mapping sample columns to basis coefficients.
        % - Topic: Inspect projection data
        forwardMatrix
        % Sampled Gram or explicitly prescribed signed coefficient system.
        % - Topic: Inspect projection data
        gramMatrix
        % Magnitude-normalized active Gram operator discrepancy.
        %
        % Scaling uses the absolute diagonal of the target Gram matrix.
        % A rank-deficient active sampled Gram reports infinity. A prescribed
        % dual reports NaN because its system is not a component Gram.
        % - Topic: Measure projection quality
        gramError
        % Operator norm of the active round-trip error, or NaN without synthesis.
        % - Topic: Measure projection quality
        roundTripError
        % Numerical rank of the active sampled Gram matrix.
        % - Topic: Measure projection quality
        sampledGramRank
        % Condition number of the active sampled basis, or NaN without synthesis.
        % - Topic: Measure projection quality
        inverseMatrixConditionNumber
        % Condition number of the active sampled Gram matrix.
        % - Topic: Measure projection quality
        gramConditionNumber
        % Whether the continuous target is positive definite on active columns.
        % - Topic: Inspect projection data
        targetGramIsPositiveDefinite
    end

    methods (Static)
        projection = fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,options)
    end

    methods (Access = private)
        self = initializeProjection(self,samplePairingMatrix,sampleGram,targetGram,majorant,activeMask,labels)
    end

    methods (Static, Hidden)
        recipe = createRecipe(options)
        projection = fromRecipe(recipe,z,weights,options)
        errors = relativeCoefficientNorm(coefficients,majorantGramMatrix,normSquared)
    end

    methods
        function self = IMProjection(sampledBasis,metricMatrix,targetGramMatrix,options)
            % Prepare a projection without making acceptance decisions.
            %
            % A diagonal target defaults to its absolute diagonal majorant.
            % A general target requires an explicit positive majorant. Active
            % target diagonal entries must be nonzero for Gram normalization.
            % Inactive columns must be numerically zero in sampled and target
            % data; their forward rows remain zero. Empty active sets are valid.
            % With no arguments, construct the canonical empty Galerkin projection.
            %
            % - Topic: Create projections
            % - Declaration: projection = IMProjection(sampledBasis,metricMatrix,targetGramMatrix,options)
            % - Parameter sampledBasis: real finite sample-by-column matrix
            % - Parameter metricMatrix: real symmetric sample pairing
            % - Parameter targetGramMatrix: real symmetric continuous pairing
            % - Parameter options.majorantGramMatrix: positive active coefficient metric
            % - Parameter options.activeColumnMask: logical row, default all columns
            % - Parameter options.columnLabels: string row of scientific labels; default positional labels
            % - Parameter options.provenance: construction identity and numerical evidence
            % - Returns self: fixed projection and numerical diagnostics
            arguments (Input)
                sampledBasis (:,:) double {mustBeReal,mustBeFinite} = []
                metricMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                targetGramMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                options.majorantGramMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                options.activeColumnMask (1,:) logical = true(1,size(sampledBasis,2))
                options.columnLabels (1,:) string = string(1:size(sampledBasis,2))
                options.provenance (1,1) struct = struct(kind="sampledBasis")
            end
            arguments (Output)
                self IMProjection
            end
            [nSamples,nColumns] = size(sampledBasis);
            if ~isequal(size(metricMatrix),[nSamples nSamples]) || ~isequal(size(targetGramMatrix),[nColumns nColumns])
                error("IMProjection:InvalidShape","The sample pairing and target Gram dimensions must match sampledBasis.");
            end
            metricMatrix = symmetricMatrix(metricMatrix,"metricMatrix");
            self.sampledBasis = sampledBasis;
            self.metricMatrix = metricMatrix;
            self.projectionKind = "galerkin";
            self.supportsGramAssessment = true;
            self.supportsRoundTripAssessment = true;
            self.provenance = options.provenance;
            sampleGram = 0.5*(sampledBasis.'*metricMatrix*sampledBasis + sampledBasis.'*metricMatrix.'*sampledBasis);
            self = self.initializeProjection(sampledBasis.'*metricMatrix,sampleGram,targetGramMatrix,options.majorantGramMatrix,options.activeColumnMask,options.columnLabels);
            active = find(self.activeColumnMask);
            for iColumn = find(~self.activeColumnMask)
                values = sampledBasis(:,iColumn);
                if norm(values,2) > 1e3*eps(max(1,norm(values,2)))
                    error("IMProjection:InvalidInactiveColumn","Inactive Galerkin columns must have numerically zero sampled values.");
                end
            end
            if isempty(active)
                self.gramError = 0;
                self.roundTripError = 0;
                self.inverseMatrixConditionNumber = NaN;
                return
            end
            targetActive = self.targetGramMatrix(active,active);
            gramActive = self.gramMatrix(active,active);
            scale = 1./sqrt(abs(diag(targetActive)));
            self.gramError = norm(scale.*(gramActive-targetActive).*scale.',2);
            if self.sampledGramRank < numel(active)
                self.gramError = Inf;
            end
            self.roundTripError = norm(self.forwardMatrix*sampledBasis-diag(double(options.activeColumnMask)),2);
            self.inverseMatrixConditionNumber = cond(sampledBasis(:,active));

        end
    end
end

function matrix = symmetricMatrix(matrix,name)
if norm(matrix-matrix.',2) > 100*eps(max(1,norm(matrix,2)))
    error("IMProjection:NonSymmetricMatrix","%s must be symmetric.",name);
end
matrix = 0.5*(matrix+matrix.');
end
