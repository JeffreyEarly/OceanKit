classdef IMProjection
    % Construct a fixed sampled projection and measure its numerical quality.
    %
    % The physical dual uses the continuous target Gram matrix,
    % $$B = H^{-1}\Phi^T W.$$
    % Direct matrix construction prepares recovery by default; pass
    % prepareRecovery=false for physical preparation alone. The higher-level
    % projectionOnGrid factories default to false.
    % Optional sampled recovery is available through the separate recovery
    % object. A positive semidefinite metric uses the shared recovery engine;
    % signed pairings use its preconditioned physical-test recovery. Deficient
    % recovery is explicitly unavailable and never replaces the physical dual.
    % The separately supplied positive majorant measures coefficient
    % errors; it never replaces the physical pairing. This object does not
    % choose a grid, fit weights, apply physical acceptance tolerances, or
    % select model counts.
    % `fromPrescribedDual` instead accepts an explicit physical source dual and its
    % coefficient system. It has no synthesis basis and cannot measure a
    % sampled-basis Gram discrepancy or round trip.
    %
    % ```matlab
    % projection = IMProjection(sampledBasis,metricMatrix,targetGramMatrix);
    % coefficients = projection.projectPhysical(values);
    % ```
    %
    % - Topic: Create projections
    % - Topic: Inspect projection data
    % - Topic: Measure projection quality
    % - Topic: Apply projections
    % - Declaration: classdef IMProjection

    properties (SetAccess = private)
        % Number of sampled source coordinates accepted by projectPhysical.
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
        % Sampled Gram or explicitly prescribed signed coefficient system.
        % - Topic: Inspect projection data
        gramMatrix
        % Magnitude-normalized active Gram operator discrepancy.
        %
        % Scaling uses the absolute diagonal of the target Gram matrix.
        % This remains the actual discrepancy when the sampled Gram
        % is deficient; its rank is a separate physical diagnostic. A
        % prescribed dual reports NaN because its system is not a component Gram.
        % - Topic: Measure projection quality
        gramError
        % Operator norm of the active round-trip error, or NaN without synthesis.
        % - Topic: Measure projection quality
        roundTripError
        % Numerical rank of the active sampled Gram relative to the positive majorant.
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
        % Prepared sample recovery, separate from physical pairing data.
        % - Topic: Inspect projection data
        recovery
        % Continuous-normalization physical coefficient dual on these samples.
        %
        % A prescribed source uses its explicitly supplied coefficient
        % system. A Galerkin record instead uses the continuous target form,
        % without correcting its physical quadrature error.
        % - Topic: Inspect projection data
        physicalDualMatrix
        % Availability of the physical coefficient functional.
        % - Topic: Inspect projection data
        physicalDualStatus
    end

    properties (Access = private)
        coefficientScale_ = zeros(0)
        recoveryOptions_ = struct(prepareRecovery=true,failureMode="diagnostic",rankTolerance=1e-12)
    end

    methods (Static)
        projection = fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,options)
    end

    methods (Access = private)
        self = initializeProjection(self,samplePairingMatrix,sampleGram,targetGram,majorant,activeMask,labels)
        coefficients = applyAssessmentOperator(self,values)
    end

    methods (Static, Hidden)
        recipe = createRecipe(options)
        projection = fromRecipe(recipe,z,weights,options)
        errors = relativeCoefficientNorm(coefficients,majorantGramMatrix,normSquared)
        [factor,isPositive] = positiveFitFactor(metric)
    end

    methods
        function self = IMProjection(sampledBasis,metricMatrix,targetGramMatrix,options)
            % Prepare a projection without making acceptance decisions.
            %
            % A diagonal target defaults to its absolute diagonal majorant.
            % A general target requires an explicit positive majorant. Active
            % target diagonal entries must be nonzero for Gram normalization.
            % Inactive synthesis columns must be exactly zero. Their physical
            % forms must be numerically zero and their forward rows stay zero.
            % Empty active sets are valid.
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
            % - Parameter options.prepareRecovery: prepare sample recovery as well as the physical record, true by default
            % - Parameter options.failureMode: error or diagnostic for unresolved sample recovery
            % - Parameter options.rankTolerance: relative rank cutoff in the declared scaled recovery coordinates
            % - Returns self: fixed projection and numerical diagnostics
            arguments (Input)
                sampledBasis (:,:) double {mustBeReal,mustBeFinite} = []
                metricMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                targetGramMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                options.majorantGramMatrix (:,:) double {mustBeReal,mustBeFinite} = []
                options.activeColumnMask (1,:) logical = true(1,size(sampledBasis,2))
                options.columnLabels (1,:) string = string(1:size(sampledBasis,2))
                options.provenance (1,1) struct = struct(kind="sampledBasis")
                options.prepareRecovery (1,1) logical = true
                options.failureMode (1,1) string {mustBeMember(options.failureMode,["error","diagnostic"])} = "diagnostic"
                options.rankTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-12
            end
            arguments (Output)
                self IMProjection
            end
            [nSamples,nColumns] = size(sampledBasis);
            if options.rankTolerance >= 1
                error("IMProjection:InvalidRankTolerance","rankTolerance must be strictly between zero and one.");
            end
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
            self.recoveryOptions_ = struct(prepareRecovery=options.prepareRecovery, ...
                failureMode=options.failureMode,rankTolerance=options.rankTolerance);
            sampleGram = 0.5*(sampledBasis.'*metricMatrix*sampledBasis + sampledBasis.'*metricMatrix.'*sampledBasis);
            self = self.initializeProjection(sampledBasis.'*metricMatrix,sampleGram,targetGramMatrix,options.majorantGramMatrix,options.activeColumnMask,options.columnLabels);
            active = find(self.activeColumnMask);
            for iColumn = find(~self.activeColumnMask)
                values = sampledBasis(:,iColumn);
                if any(values ~= 0)
                    error("IMProjection:InvalidInactiveColumn","Inactive Galerkin synthesis columns must be exactly zero; declare activity from the scientific family rather than sample magnitude.");
                end
            end
            self.roundTripError = NaN;
            self.inverseMatrixConditionNumber = NaN;
            if isempty(active)
                self.gramError = 0;
            else
                targetActive = self.targetGramMatrix(active,active);
                gramActive = self.gramMatrix(active,active);
                scale = 1./sqrt(abs(diag(targetActive)));
                self.gramError = norm(scale.*(gramActive-targetActive).*scale.',2);
                self.inverseMatrixConditionNumber = cond(sampledBasis(:,active));
            end
            self.supportsRoundTripAssessment = options.prepareRecovery;
            if options.prepareRecovery
                [fitFactor,positiveFit] = IMProjection.positiveFitFactor(metricMatrix);
                if positiveFit
                    prepared = prepareIMRecovery(sampledBasis,fitFactor=fitFactor, ...
                        coefficientScale=self.coefficientScale_,physicalDual=self.physicalDualMatrix, ...
                        activeColumnMask=self.activeColumnMask,columnLabels=self.columnLabels, ...
                        rankTolerance=options.rankTolerance,failureMode=options.failureMode);
                else
                    if self.physicalDualStatus ~= "available"
                        error("IMProjection:UnavailablePhysicalDual","Corrected physical recovery requires a representable normalized physical dual (%s). Use prepareRecovery=false to inspect the physical record, or change the declared coefficient units.",self.physicalDualStatus);
                    end
                    prepared = prepareIMRecovery(sampledBasis,policy="correctedPhysicalDual", ...
                        physicalDual=self.physicalDualMatrix,coefficientScale=self.coefficientScale_, ...
                        activeColumnMask=self.activeColumnMask,columnLabels=self.columnLabels, ...
                        rankTolerance=options.rankTolerance,failureMode=options.failureMode);
                end
                self.recovery = IMDiscreteTransform.fromPreparedRecovery(sampledBasis,prepared, ...
                    metadata=struct(construction="IMProjection",projectionProvenance=self.provenance));
                if prepared.success
                    self.roundTripError = prepared.diagnostics.roundTripError;
                else
                    self.roundTripError = Inf;
                end
            end

        end
    end
end

function matrix = symmetricMatrix(matrix,name)
if norm(matrix-matrix.',2) > 100*eps(max(1,norm(matrix,2)))
    error("IMProjection:NonSymmetricMatrix","%s must be symmetric.",name);
end
matrix = 0.5*(matrix+matrix.');
end
