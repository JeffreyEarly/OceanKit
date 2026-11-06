classdef IMDiscreteTransform
    % Apply one prepared recovery and synthesis pair to arbitrary observations.
    %
    % For observations y and scientific coefficients a, synthesis is y=A*a
    % and recovery is a=R*y. Observation rows may contain several variables,
    % endpoint coordinates or traces. Physical quadrature, signed pairings and
    % energy forms belong to a separate IMProjection, not this generic object.
    %
    % Construct with fromObservationMatrices. The shared preparation kernel
    % makes rank and policy decisions once; application only multiplies cached
    % matrices and preserves trailing field dimensions.
    %
    % ```matlab
    % transform = IMDiscreteTransform.fromObservationMatrices(A,observationFactor=eye(size(A,1)),columnLabels=labels);
    % a = transform.transformForward(values);
    % values = transform.transformBack(a);
    % ```
    %
    % - Topic: Create transforms
    % - Topic: Inspect transforms
    % - Topic: Apply transforms
    % - Declaration: classdef IMDiscreteTransform

    properties (SetAccess = private)
        % Number of observation rows, independent of physical point count.
        % - Topic: Inspect transforms
        sampleCount (1,1) double
        % Number of aligned scientific coefficient coordinates.
        % - Topic: Inspect transforms
        columnCount (1,1) double
        % Complex observation-by-coefficient synthesis A.
        % - Topic: Inspect transforms
        synthesisMatrix (:,:) double
        % Complex coefficient-by-observation recovery R; empty on failure.
        % - Topic: Inspect transforms
        forwardMatrix (:,:) double
        % Unique scientific coefficient identities in retained order.
        % - Topic: Inspect transforms
        columnLabels (1,:) string
        % Scientifically declared active coordinates; never inferred from zeros.
        % - Topic: Inspect transforms
        activeColumnMask (1,:) logical
        % Uncorrected physical coefficient dual when explicitly supplied.
        % - Topic: Inspect transforms
        physicalDualMatrix (:,:) double
        % Explicit numerical recovery construction policy.
        % - Topic: Inspect transforms
        recoveryPolicy (1,1) string
        % Whether the requested numerical construction succeeded.
        % - Topic: Inspect transforms
        success (1,1) logical
        % Construction status, including explicit diagnostic-only failures.
        % - Topic: Inspect transforms
        status (1,1) string
        % Caller-supplied observation layout and scientific provenance.
        % - Topic: Inspect transforms
        metadata (1,1) struct
        % Immutable factors, scaling, thresholds and declared extended span.
        %
        % This is shared-kernel numerical provenance, separate from caller
        % metadata. Output views preserve the original preparation and
        % composed retained-column map; their retained coefficient scale
        % remains undeclared until explicitly prepared.
        % - Topic: Inspect transforms
        recoverySpecification (1,1) struct
        % Recovery diagnostics with independently declared norm capabilities.
        % - Topic: Inspect transforms
        diagnostics (1,1) struct
    end

    methods (Access = protected)
        function self = IMDiscreteTransform(prepared)
            % Initialize trusted prepared state without another factorization.
            if isa(prepared,"IMDiscreteTransform")
                self = copyPrepared(self,prepared);
                return
            end
            self.synthesisMatrix = prepared.synthesisMatrix;
            self.sampleCount = size(prepared.synthesisMatrix,1);
            self.columnCount = size(prepared.synthesisMatrix,2);
            self.forwardMatrix = prepared.forwardMatrix;
            self.columnLabels = prepared.columnLabels;
            self.activeColumnMask = prepared.activeColumnMask;
            self.physicalDualMatrix = prepared.physicalDual;
            self.recoveryPolicy = prepared.policy;
            self.success = prepared.success;
            self.status = prepared.status;
            self.metadata = prepared.metadata;
            self.recoverySpecification = prepared.specification;
            self.diagnostics = prepared.diagnostics;
        end
    end

    methods (Static)
        function transform = fromObservationMatrices(synthesisMatrix,options)
            % Prepare recovery with explicit fitting and observation norms.
            %
            % sampledLeastSquares fits F*A directly. correctedPhysicalDual
            % solves the scaled weak overlap while preserving the supplied B.
            % extendedSpan prepares a declared larger family before selecting
            % coefficient rows. No policy silently drops active directions.
            % Algebraic weak recovery may omit an observation factor; its
            % observation/noise geometry then remains unavailable.
            %
            % - Topic: Create transforms
            % - Parameter synthesisMatrix: finite real/complex sample-by-column A
            % - Parameter options.policy: sampledLeastSquares, correctedPhysicalDual or extendedSpan
            % - Parameter options.columnLabels: unique coefficient identities
            % - Parameter options.activeColumnMask: declared scientific activity
            % - Parameter options.coefficientScale: nonsingular active coefficient scaling C
            % - Parameter options.observationFactor: nonsingular positive observation whitening D
            % - Parameter options.fitFactor: optional rectangular least-squares fitting factor F
            % - Parameter options.physicalDual: uncorrected aligned physical coefficient dual B
            % - Parameter options.rankTolerance: dimensionless relative rank cutoff
            % - Parameter options.failureMode: error or diagnostic
            % - Parameter options.extendedSynthesis: declared larger synthesis span
            % - Parameter options.extendedColumnLabels: identities of extended columns
            % - Parameter options.extendedActiveColumnMask: extended scientific activity
            % - Parameter options.extendedCoefficientScale: extended coefficient scaling
            % - Parameter options.retainedColumnIndices: retained locations in extended span
            % - Parameter options.extendedPolicy: recovery policy on the extended span
            % - Parameter options.extendedPhysicalDual: physical dual for extended recovery
            % - Parameter options.metadata: caller observation-layout/provenance record
            % - Returns transform: immutable matrices and construction diagnostics
            arguments (Input)
                synthesisMatrix (:,:) double {mustBeFinite}
                options.policy (1,1) string = "sampledLeastSquares"
                options.columnLabels (1,:) string = string(1:size(synthesisMatrix,2))
                options.activeColumnMask (1,:) logical = true(1,size(synthesisMatrix,2))
                options.coefficientScale (:,:) double {mustBeFinite} = []
                options.observationFactor (:,:) double {mustBeFinite} = []
                options.fitFactor (:,:) double {mustBeFinite} = []
                options.physicalDual (:,:) double {mustBeFinite} = []
                options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
                options.failureMode (1,1) string = "error"
                options.extendedSynthesis (:,:) double {mustBeFinite} = []
                options.extendedColumnLabels (1,:) string = strings(1,0)
                options.extendedActiveColumnMask (1,:) logical = true(1,0)
                options.extendedCoefficientScale (:,:) double {mustBeFinite} = []
                options.retainedColumnIndices (1,:) double {mustBeInteger,mustBePositive} = []
                options.extendedPolicy (1,1) string = "sampledLeastSquares"
                options.extendedPhysicalDual (:,:) double {mustBeFinite} = []
                options.metadata (1,1) struct = struct()
            end
            args = rmfield(options,"metadata");
            pairs = namedargs2cell(args);
            prepared = prepareIMRecovery(synthesisMatrix,pairs{:});
            transform = IMDiscreteTransform.fromPreparedRecovery(synthesisMatrix,prepared,metadata=options.metadata);
        end
    end

    methods (Static, Hidden)
        function transform = fromPreparedRecovery(synthesisMatrix,prepared,options)
            % Bind the trusted shared-kernel result without preparing it twice.
            % - Topic: Create transforms
            % - Developer: true
            arguments (Input)
                synthesisMatrix (:,:) double {mustBeFinite}
                prepared (1,1) struct
                options.metadata (1,1) struct = struct()
            end
            nSamples = size(synthesisMatrix,1); nColumns = size(synthesisMatrix,2);
            required = ["forwardMatrix","physicalDual","success","status","policy","columnLabels","activeColumnMask","diagnostics","specification"];
            if ~all(isfield(prepared,required)) || numel(prepared.columnLabels) ~= nColumns || numel(prepared.activeColumnMask) ~= nColumns
                error("IMDiscreteTransform:InvalidPreparedRecovery","Use the complete shared-kernel result aligned with synthesisMatrix.");
            end
            if prepared.success && (~isequal(size(prepared.forwardMatrix),[nColumns nSamples]) || any(~isfinite(prepared.forwardMatrix),"all"))
                error("IMDiscreteTransform:InvalidPreparedRecovery","Successful recovery requires one finite row per coefficient and one column per observation.");
            end
            if ~prepared.success && ~isempty(prepared.forwardMatrix)
                error("IMDiscreteTransform:InvalidPreparedRecovery","A failed preparation must not expose an applicable recovery matrix.");
            end
            prepared.synthesisMatrix = synthesisMatrix;
            prepared.metadata = options.metadata;
            transform = IMDiscreteTransform(prepared);
        end
    end

    methods
        function coefficients = transformForward(self,values)
            % Recover coefficients along the leading observation dimension.
            % - Topic: Apply transforms
            % - Parameter values: observation-by-field array with arbitrary trailing dimensions
            % - Returns coefficients: coefficient-by-field array preserving trailing dimensions
            arguments (Input)
                self (1,1) IMDiscreteTransform
                values double
            end
            if size(values,1) ~= self.sampleCount
                error("IMDiscreteTransform:InvalidSampleCount","values must have sampleCount rows in the declared observation order.");
            end
            if ~self.success
                error("IMDiscreteTransform:UnavailableRecovery","The requested recovery is unavailable: %s.",self.status);
            end
            coefficients = applyLeading(self.forwardMatrix,values,self.columnCount);
        end

        function values = transformBack(self,coefficients)
            % Synthesize observations along the leading coefficient dimension.
            % - Topic: Apply transforms
            % - Parameter coefficients: coefficient-by-field array preserving scientific order
            % - Returns values: observation-by-field array preserving trailing dimensions
            arguments (Input)
                self (1,1) IMDiscreteTransform
                coefficients double
            end
            if size(coefficients,1) ~= self.columnCount
                error("IMDiscreteTransform:InvalidCoefficientCount","coefficients must have columnCount rows in the declared scientific order.");
            end
            values = applyLeading(self.synthesisMatrix,coefficients,self.sampleCount);
        end

        function selected = selectOutputs(self,columns)
            % Keep output rows of this recovery without refitting its span.
            %
            % This is an extended-span view, not least squares reconstructed
            % on a smaller family. The original construction diagnostics stay
            % labeled by their original span. Only the unscaled selected
            % round trip is recomputed; norm-dependent selected diagnostics
            % require a new explicit preparation.
            %
            % - Topic: Apply transforms
            % - Parameter columns: distinct retained array columns in caller order
            % - Returns selected: generic IMDiscreteTransform coefficient view
            arguments (Input)
                self (1,1) IMDiscreteTransform
                columns (1,:) double {mustBeInteger,mustBePositive}
            end
            if any(columns > self.columnCount) || numel(unique(columns)) ~= numel(columns)
                error("IMDiscreteTransform:InvalidColumns","Select distinct existing coefficient columns.");
            end
            selected = IMDiscreteTransform(self);
            selected.synthesisMatrix = self.synthesisMatrix(:,columns);
            selected.columnCount = numel(columns);
            selected.columnLabels = self.columnLabels(columns);
            selected.activeColumnMask = self.activeColumnMask(columns);
            if self.success, selected.forwardMatrix = self.forwardMatrix(columns,:); end
            if ~isempty(self.physicalDualMatrix), selected.physicalDualMatrix = self.physicalDualMatrix(columns,:); end
            selected.recoveryPolicy = "extendedSpan";
            specification = self.recoverySpecification;
            hasExtension = specification.policy == "extendedSpan" && isfield(specification,"extendedRecoverySpecification") && isfield(specification.extendedRecoverySpecification,"policy");
            if hasExtension
                origin = specification.extendedRecoverySpecification;
                originSynthesis = specification.extendedSynthesis;
                originDual = specification.extendedPhysicalDual;
                retained = specification.retainedColumnIndices(columns);
                originDiagnostics = self.diagnostics.extended;
            else
                origin = specification;
                originSynthesis = self.synthesisMatrix;
                originDual = self.physicalDualMatrix;
                retained = columns;
                originDiagnostics = self.diagnostics;
            end
            selected.recoverySpecification = struct(policy="extendedSpan",columnLabels=selected.columnLabels,activeColumnMask=selected.activeColumnMask,coefficientScale=[],observationFactor=origin.observationFactor,fitFactor=origin.fitFactor,rankTolerance=origin.rankTolerance,failureMode=origin.failureMode,extendedSynthesis=originSynthesis,extendedColumnLabels=origin.columnLabels,extendedActiveColumnMask=origin.activeColumnMask,extendedCoefficientScale=origin.coefficientScale,retainedColumnIndices=retained,extendedPolicy=origin.policy,extendedPhysicalDual=originDual,extendedRecoverySpecification=origin,retainedCoefficientScaleAvailable=false);
            unavailable = struct(available=false,status="unavailable",rank=NaN,columnCount=NaN,conditionNumber=NaN,singularValues=zeros(0,1),scaledSingularValues=zeros(0,1),singularValueScaleFactors=zeros(1,0),singularValuesRepresentable=false,relativeTolerance=origin.rankTolerance,absoluteTolerance=NaN);
            selected.diagnostics = struct(observation=unavailable,fit=unavailable,test=unavailable,overlap=unavailable,roundTripError=NaN,scaledRoundTripError=NaN,roundTripStatus="unavailable",definingEquationResidual=NaN,definingEquationResidualStatus="unavailable",scaledDefiningEquationResidual=NaN,observationNoiseGain=NaN,observationNoiseGainStatus="unavailable",noiseOptimalForObservationNorm=false,extendedRejectionError=NaN,extendedRejectionStatus="unavailable",principalCosines=zeros(0,1),supportsPrincipalGeometry=false,extended=originDiagnostics,constructionColumnLabels=origin.columnLabels,scope="fixed recovery output view");
            if self.success
                delta = selected.forwardMatrix*selected.synthesisMatrix-diag(double(selected.activeColumnMask));
                [selected.diagnostics.roundTripError,selected.diagnostics.roundTripStatus] = finiteNorm(delta);
                omitted = origin.activeColumnMask & ~ismember(1:size(originSynthesis,2),retained);
                [selected.diagnostics.extendedRejectionError,selected.diagnostics.extendedRejectionStatus] = finiteNorm(selected.forwardMatrix*originSynthesis(:,omitted));
            end
        end
    end

    methods (Access = private)
        function self = copyPrepared(self,prepared)
            names = ["sampleCount","columnCount","synthesisMatrix","forwardMatrix","columnLabels","activeColumnMask","physicalDualMatrix","recoveryPolicy","success","status","metadata","recoverySpecification","diagnostics"];
            for name = names, self.(name) = prepared.(name); end
        end
    end
end

function result = applyLeading(matrix,values,nRows)
shape = size(values);
nPages = prod(shape(2:end));
result = reshape(matrix*reshape(values,size(values,1),nPages),[nRows shape(2:end)]);
end

function [value,status] = finiteNorm(matrix)
value = NaN;
status = "unrepresentable";
if all(isfinite(matrix),"all")
    measured = norm(matrix,2);
    if isfinite(measured), value = measured; status = "measured"; end
end
end
