classdef IMInternalModesDiscreteTransformAssessment
    % Store retained-band diagnostics for aligned F/G transforms.
    %
    % One family-specific point-and-weight rule is assessed for every
    % leading family prefix. Certified construction records independently
    % refitted candidate counts in `certificationSearch`; its final
    % `weightFitModeCount` equals `retainedModeCount`. `prefixDiagnostics`
    % reports the worst requested channel for
    % each policy, while `variablePrefixDiagnostics` exposes the underlying
    % per-variable Gram, rank, round-trip, and conditioning records.
    % The combined table includes `modeCount`, `lastModeNumber`, the Gram
    % error and limiting variable, leakage error and limiting rejected mode,
    % coupled quadratic error and limiting product channel/source labels,
    % and cumulative per-policy and combined acceptance flags.
    % Policy structs expose `tolerance`, per-prefix `error` and `accepted`
    % arrays, `maximumAcceptedModeCount`, `limitingValue`, and `reason`.
    % Internal-mode quadratic products are projected with the signed
    % Pontryagin pairing and signed Gram solve. Their relative error is
    % measured with the induced positive Hilbert majorant, recorded as
    % `quadraticAliasingPolicy.errorNorm`. Signed rejected-mode leakage
    % and MDA coupled quadratic products remain unsupported.
    %
    % - Topic: Create transform assessments
    % - Topic: Inspect transform assessments
    % - Topic: Inspect retained-band policies
    % - Declaration: classdef IMInternalModesDiscreteTransformAssessment

    properties (SetAccess = private)
        % Production transform for the common accepted family prefix.
        transform
        % Full candidate-family transform on the fixed rule.
        candidateTransform
        % Shared weight fit, or empty for caller-supplied weights.
        weightFit
        % Exact requested physical point count.
        requestedPointCount
        % Actual physical point count.
        actualPointCount
        % Number of candidate family columns.
        candidateModeCount
        % Number of commonly retained family columns.
        retainedModeCount
        % Physical labels for candidate family columns.
        candidateModeNumber
        % Physical labels for retained family columns.
        retainedModeNumber
        % Direct forward channels assessed in canonical order.
        availableVariables
        % Variable responsible for the combined retained limit.
        limitingVariable
        % Worst-channel and coupled-product prefix diagnostics.
        prefixDiagnostics
        % Worst-channel normalized-Gram policy result.
        gramPolicy
        % Worst-channel rejected-mode leakage policy result.
        leakagePolicy
        % Coupled quadratic policy: signed projection, majorant error norm.
        quadraticAliasingPolicy
        % Policy responsible for the retained count.
        limitingPolicy
        % Readable retained-band explanation.
        retentionReason
        % Provenance of the physical sample grid.
        gridDesign
        % Independently refitted count-selection attempts.
        certificationSearch
        % Number of family modes used to fit the stored weights.
        weightFitModeCount
    end

    properties (Access = private)
        variableDiagnostics
    end

    methods
        function self = IMInternalModesDiscreteTransformAssessment(options)
            % Create an aligned internal-mode transform assessment.
            % - Topic: Developer topics — Construct assessments
            % - Declaration: assessment = IMInternalModesDiscreteTransformAssessment(options)
            % - Developer: true
            arguments
                options.transform (1,1) IMInternalModesDiscreteTransform
                options.candidateTransform (1,1) IMInternalModesDiscreteTransform
                options.weightFit = []
                options.requestedPointCount (1,1) double {mustBeInteger, mustBePositive}
                options.prefixDiagnostics table
                options.variableDiagnostics (1,1) struct
                options.gramPolicy (1,1) struct
                options.leakagePolicy (1,1) struct
                options.quadraticAliasingPolicy (1,1) struct
                options.limitingPolicy {mustBeTextScalar}
                options.limitingVariable {mustBeTextScalar}
                options.retentionReason {mustBeTextScalar}
                options.gridDesign struct = struct.empty
                options.certificationSearch table = table()
            end
            if ~isempty(options.weightFit) && ~isa(options.weightFit,"IMInternalModesQuadratureWeightFit")
                error("IMInternalModesDiscreteTransformAssessment:InvalidWeightFit", "weightFit must be empty or an IMInternalModesQuadratureWeightFit.");
            end
            candidateModeCount = length(options.candidateTransform.modeNumber);
            retainedModeCount = length(options.transform.modeNumber);
            if height(options.prefixDiagnostics) ~= candidateModeCount
                error("IMInternalModesDiscreteTransformAssessment:InvalidPrefixDiagnostics", "prefixDiagnostics must have one row per candidate family mode.");
            end
            if retainedModeCount > candidateModeCount || ~isequal(options.transform.modeNumber,options.candidateTransform.modeNumber(1:retainedModeCount))
                error("IMInternalModesDiscreteTransformAssessment:InvalidRetainedPrefix", "transform must be a leading prefix of candidateTransform.");
            end
            if ~isequal(options.transform.z,options.candidateTransform.z) || ~isequal(options.transform.weights,options.candidateTransform.weights)
                error("IMInternalModesDiscreteTransformAssessment:IncompatibleTransforms", "Retained and candidate transforms must share points and weights.");
            end
            for variable = options.candidateTransform.availableVariables
                field = char(variable);
                if ~isfield(options.variableDiagnostics,field) || height(options.variableDiagnostics.(field)) ~= candidateModeCount
                    error("IMInternalModesDiscreteTransformAssessment:InvalidVariableDiagnostics", "Each available variable needs one diagnostic row per candidate prefix.");
                end
            end
            if ~isempty(options.weightFit) && ~isequal(options.weightFit.transform.weights,options.candidateTransform.weights)
                error("IMInternalModesDiscreteTransformAssessment:IncompatibleWeightFit", "weightFit must describe the candidate transform's fixed rule.");
            end

            self.transform = options.transform;
            self.candidateTransform = options.candidateTransform;
            self.weightFit = options.weightFit;
            self.requestedPointCount = options.requestedPointCount;
            self.actualPointCount = length(options.candidateTransform.z);
            self.candidateModeCount = candidateModeCount;
            self.retainedModeCount = retainedModeCount;
            self.candidateModeNumber = options.candidateTransform.modeNumber;
            self.retainedModeNumber = options.transform.modeNumber;
            self.availableVariables = options.candidateTransform.availableVariables;
            self.limitingVariable = string(options.limitingVariable);
            self.prefixDiagnostics = options.prefixDiagnostics;
            self.variableDiagnostics = options.variableDiagnostics;
            self.gramPolicy = options.gramPolicy;
            self.leakagePolicy = options.leakagePolicy;
            self.quadraticAliasingPolicy = options.quadraticAliasingPolicy;
            self.limitingPolicy = string(options.limitingPolicy);
            self.retentionReason = string(options.retentionReason);
            if isempty(options.gridDesign)
                options.gridDesign = struct(kind="unknown",pointCount=length(options.candidateTransform.z));
            end
            self.gridDesign = options.gridDesign;
            self.certificationSearch = options.certificationSearch;
            if isempty(options.weightFit)
                self.weightFitModeCount = NaN;
            else
                self.weightFitModeCount = length(options.weightFit.transform.modeNumber);
            end
        end

        function diagnostics = variablePrefixDiagnostics(self, options)
            % Return the per-prefix diagnostic table for F or G.
            arguments
                self IMInternalModesDiscreteTransformAssessment
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = "F"
            end
            variable = string(options.variable);
            if ~ismember(variable,self.availableVariables)
                reason = self.candidateTransform.forwardTransformReason(variable=variable);
                error("IMInternalModesDiscreteTransformAssessment:UnavailableVariable", "%s",reason);
            end
            diagnostics = self.variableDiagnostics.(char(variable));
        end

        function assessment = withCertificationMetadata(self,gridDesign,certificationSearch)
            % Return this assessment with grid and count-search provenance.
            arguments
                self IMInternalModesDiscreteTransformAssessment
                gridDesign (1,1) struct
                certificationSearch table
            end
            assessment = self;
            assessment.gridDesign = gridDesign;
            assessment.certificationSearch = certificationSearch;
        end
    end
end
