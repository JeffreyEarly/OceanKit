classdef IMDiscreteTransformAssessment
    % Store fixed-rule retained-band diagnostics for a scalar transform.
    %
    % `IMDiscreteTransformAssessment` records one final point-and-weight rule
    % and the leading modal prefixes assessed on that rule. For certified
    % construction, `certificationSearch` records the independently refitted
    % counts considered before this exact final band was selected, and
    % `weightFitModeCount` equals `retainedModeCount`. Lower-level fixed-rule
    % calls can still fit a larger candidate once and return a shorter
    % prefix; their differing counts make that behavior explicit.
    %
    % For prefix $$N$$, let $$\Gamma_N$$ be its sampled Gram matrix,
    % $$\Gamma_{0,N}$$ its diagonal continuous target, and
    %
    % $$
    % S_N=\operatorname{diag}\!\left(
    % \left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}
    % \right),\qquad
    % E_N=S_N(\Gamma_N-\Gamma_{0,N})S_N.
    % $$
    %
    % The default Gram policy uses
    % $$\lVert E_N\rVert_2$$, stored in `gramError`. By contrast, the
    % built-in quadrature objective reports $$\lVert E_N\rVert_\mathrm F$$
    % as `IMQuadratureWeightFit.residualNorm`. `roundTripError` is
    % $$\lVert A_\mathrm fA_\mathrm i-I\rVert_2$$; it measures algebraic
    % consistency and can be near roundoff even when the Gram rule is poor.
    % The two condition-number columns describe algebraic sensitivity and
    % are diagnostics rather than acceptance policies.
    %
    % When enabled, rejected-mode leakage is
    %
    % $$
    % \ell_N=\max_{N<j\leq N_\mathrm{check}}
    % \frac{\lVert\Pi_N^\mathrm{discrete}u_j\rVert_\mu}
    % {\lVert u_j\rVert_\mu},
    % $$
    %
    % and scalar quadratic aliasing is
    %
    % $$
    % q_N=\max_{1\leq i\leq j\leq N}
    % \frac{\left\lVert
    % \Pi_N^\mathrm{discrete}(u_i u_j)-
    % \Pi_N^\mathrm{continuous}(u_i u_j)
    % \right\rVert_\mu}{\lVert u_i u_j\rVert_\mu}.
    % $$
    %
    % The continuous quadratic reference is integrated on the source
    % solver's inner-product grid, independently of the fitted points and
    % weights. Leakage and quadratic aliasing require a positive-definite
    % target because $$\lVert\cdot\rVert_\mu$$ must be a true norm. Signed
    % targets retain Gram diagnostics but cannot enable those policies.
    %
    % `prefixDiagnostics` has one row per candidate prefix. Its columns are:
    %
    % - `modeCount`: prefix column count $$N$$.
    % - `lastModeNumber`: physical label of column $$N$$.
    % - `gramError`: $$\lVert E_N\rVert_2$$.
    % - `roundTripError`: $$\lVert A_\mathrm fA_\mathrm i-I\rVert_2$$.
    % - `inverseMatrixConditionNumber`: $$\kappa_2(A_\mathrm i)$$.
    % - `gramConditionNumber`: $$\kappa_2(\Gamma_N)$$.
    % - `sampledGramRank`: numerical rank of $$\Gamma_N$$.
    % - `leakageError`: $$\ell_N$$, or `NaN` when disabled.
    % - `leakageLimitingModeNumber`: rejected physical mode label attaining $$\ell_N$$.
    % - `quadraticAliasingError`: $$q_N$$, or `NaN` when disabled.
    % - `quadraticLimitingModeNumberI`, `quadraticLimitingModeNumberJ`: physical labels of the product attaining $$q_N$$.
    % - `gramAccepted`, `leakageAccepted`, `quadraticAccepted`: cumulative per-policy prefix decisions.
    % - `combinedAccepted`: intersection of every enabled cumulative policy.
    %
    % Acceptance is consecutive: after a policy first rejects a prefix,
    % larger prefixes remain rejected even if a later raw metric happens to
    % fall below tolerance. Each policy struct stores `enabled`, `tolerance`,
    % the per-prefix `error` and `accepted` arrays, `maximumAcceptedModeCount`,
    % `limitingValue`, and a readable `reason`. Leakage additionally stores
    % `nCheckModes` and its limiting mode labels; quadratic aliasing stores
    % its limiting product labels.
    %
    % ```matlab
    % [transform,assessment] = basisSet.discreteTransform(nPoints=32);
    % assessment.prefixDiagnostics
    % assessment.gramPolicy
    % assessment.leakagePolicy
    % assessment.quadraticAliasingPolicy
    % ```
    %
    % - Topic: Create transform assessments
    % - Topic: Inspect transform assessments
    % - Topic: Inspect retained-band policies
    % - Declaration: classdef IMDiscreteTransformAssessment

    properties (SetAccess = private)
        % Production transform using the retained modal prefix.
        %
        % - Topic: Inspect transform assessments
        transform

        % Full candidate-band transform built with the fixed rule.
        %
        % - Topic: Inspect transform assessments
        candidateTransform

        % Quadrature-weight fit used to build the fixed rule.
        %
        % This is empty when the caller supplied explicit weights. When it
        % is present, `weightFit.transform` is the full candidate transform,
        % even when policies reduce the returned production prefix.
        %
        % - Topic: Inspect transform assessments
        weightFit

        % Requested physical sample count.
        %
        % For `discreteTransform(nPoints=Nz)`, this is the exact requested
        % count. For explicit `z`, it is `length(z)`.
        %
        % - Topic: Inspect transform assessments
        requestedPointCount

        % Actual physical sample count.
        %
        % Exact point-count construction guarantees this equals
        % `requestedPointCount`; it is recorded separately so the contract
        % remains directly inspectable.
        %
        % - Topic: Inspect transform assessments
        actualPointCount

        % Number of modes in the full candidate band.
        %
        % - Topic: Inspect transform assessments
        candidateModeCount

        % Number of modes retained by the combined policies.
        %
        % - Topic: Inspect transform assessments
        retainedModeCount

        % Physical labels for the candidate modes.
        %
        % - Topic: Inspect transform assessments
        candidateModeNumber

        % Physical labels for the retained modes.
        %
        % - Topic: Inspect transform assessments
        retainedModeNumber

        % Per-prefix transform and policy diagnostics.
        %
        % The table has one row for every leading candidate prefix and uses
        % physical mode labels rather than assuming labels equal indices.
        % See the class overview for every column definition.
        %
        % - Topic: Inspect transform assessments
        prefixDiagnostics

        % Normalized-Gram retained-band policy result.
        %
        % This policy is always enabled. Its shipped default tolerance is
        % $$10^{-2}$$, selected by constant/exponential coordinate sweeps;
        % callers can override it through `discreteTransform`.
        %
        % - Topic: Inspect retained-band policies
        gramPolicy

        % Rejected-mode leakage policy result.
        %
        % This policy is disabled when its tolerance is empty. With the
        % default check band it samples rejected modes through twice the
        % candidate count.
        %
        % - Topic: Inspect retained-band policies
        leakagePolicy

        % Scalar quadratic-aliasing policy result.
        %
        % This policy is disabled when its tolerance is empty. It applies to
        % products projected back into the same scalar basis; coupled `F/G`
        % product channels are outside this class.
        %
        % - Topic: Inspect retained-band policies
        quadraticAliasingPolicy

        % Policy that imposed the final retained count.
        %
        % - Topic: Inspect retained-band policies
        limitingPolicy

        % Readable explanation of the retained-band decision.
        %
        % - Topic: Inspect retained-band policies
        retentionReason

        % Provenance of the physical sample grid.
        gridDesign

        % Independently refitted count-selection attempts.
        %
        % Empty for a direct fixed-band assessment. Certified construction
        % records each attempted count and whether its freshly fitted rule
        % passed the relevant selection stage. Search-only Gram rows record
        % `retainedModeCount=NaN` when the full candidate was rejected,
        % because rejected candidate prefixes are intentionally not built.
        certificationSearch

        % Number of modes used when fitting the stored quadrature weights.
        %
        % This equals `retainedModeCount` for a certified transform. It can
        % exceed it for legacy fixed-rule prefix assessment.
        weightFitModeCount
    end

    methods
        function self = IMDiscreteTransformAssessment(options)
            % Create a scalar discrete-transform assessment.
            %
            % - Topic: Create transform assessments
            % - Declaration: assessment = IMDiscreteTransformAssessment(options)
            % - Parameter options.transform: retained production transform
            % - Parameter options.candidateTransform: full candidate transform
            % - Parameter options.weightFit: quadrature-weight fit or empty
            % - Parameter options.requestedPointCount: requested point count
            % - Parameter options.prefixDiagnostics: per-prefix diagnostic table
            % - Parameter options.gramPolicy: Gram policy result
            % - Parameter options.leakagePolicy: leakage policy result
            % - Parameter options.quadraticAliasingPolicy: quadratic policy result
            % - Parameter options.limitingPolicy: limiting policy name
            % - Parameter options.retentionReason: readable retention reason
            % - Returns assessment: initialized assessment
            arguments
                options.transform (1,1) IMDiscreteTransform
                options.candidateTransform (1,1) IMDiscreteTransform
                options.weightFit = []
                options.requestedPointCount (1,1) double {mustBeInteger, mustBePositive}
                options.prefixDiagnostics table
                options.gramPolicy (1,1) struct
                options.leakagePolicy (1,1) struct
                options.quadraticAliasingPolicy (1,1) struct
                options.limitingPolicy {mustBeTextScalar}
                options.retentionReason {mustBeTextScalar}
                options.gridDesign struct = struct.empty
                options.certificationSearch table = table()
            end

            if ~isempty(options.weightFit) && ~isa(options.weightFit,"IMQuadratureWeightFit")
                error("IMDiscreteTransformAssessment:InvalidWeightFit", "weightFit must be empty or an IMQuadratureWeightFit.");
            end
            candidateModeCount = length(options.candidateTransform.modeNumber);
            retainedModeCount = length(options.transform.modeNumber);
            if height(options.prefixDiagnostics) ~= candidateModeCount
                error("IMDiscreteTransformAssessment:InvalidPrefixDiagnostics", "prefixDiagnostics must have one row per candidate mode.");
            end
            if retainedModeCount > candidateModeCount || ~isequal(options.transform.modeNumber,options.candidateTransform.modeNumber(1:retainedModeCount))
                error("IMDiscreteTransformAssessment:InvalidRetainedPrefix", "transform must contain a leading prefix of candidateTransform modes.");
            end
            if ~isequal(options.transform.z,options.candidateTransform.z) || ~isequal(options.transform.weights,options.candidateTransform.weights)
                error("IMDiscreteTransformAssessment:IncompatibleTransforms", "The retained and candidate transforms must use identical points and weights.");
            end
            requiredColumns = ["modeCount" "lastModeNumber" "gramError" "roundTripError" ...
                "inverseMatrixConditionNumber" "gramConditionNumber" "sampledGramRank" "leakageError" "leakageLimitingModeNumber" ...
                "quadraticAliasingError" "quadraticLimitingModeNumberI" "quadraticLimitingModeNumberJ" ...
                "gramAccepted" "leakageAccepted" "quadraticAccepted" "combinedAccepted"];
            if ~all(ismember(requiredColumns,string(options.prefixDiagnostics.Properties.VariableNames)))
                error("IMDiscreteTransformAssessment:InvalidPrefixDiagnostics", "prefixDiagnostics does not contain the complete scalar retained-band diagnostic schema.");
            end
            policies = {options.gramPolicy,options.leakagePolicy,options.quadraticAliasingPolicy};
            requiredPolicyFields = ["enabled" "tolerance" "error" "accepted" "maximumAcceptedModeCount" "limitingValue" "reason"];
            for iPolicy = 1:length(policies)
                policy = policies{iPolicy};
                if ~all(isfield(policy,requiredPolicyFields)) || length(policy.error) ~= candidateModeCount || length(policy.accepted) ~= candidateModeCount
                    error("IMDiscreteTransformAssessment:InvalidPolicy", "Every policy must contain the complete result schema and one error and acceptance value per candidate prefix.");
                end
            end
            if ~isempty(options.weightFit)
                fitTransform = options.weightFit.transform;
                if ~isequal(fitTransform.z,options.candidateTransform.z) ...
                        || ~isequal(fitTransform.weights,options.candidateTransform.weights) ...
                        || ~isequal(fitTransform.modeNumber,options.candidateTransform.modeNumber) ...
                        || ~isequal(fitTransform.forwardMatrix,options.candidateTransform.forwardMatrix)
                    error("IMDiscreteTransformAssessment:IncompatibleWeightFit", "weightFit.transform must be the full candidate transform used by this assessment.");
                end
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
            self.prefixDiagnostics = options.prefixDiagnostics;
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

        function assessment = withCertificationMetadata(self,gridDesign,certificationSearch)
            % Return this assessment with grid and count-search provenance.
            arguments
                self IMDiscreteTransformAssessment
                gridDesign (1,1) struct
                certificationSearch table
            end
            assessment = self;
            assessment.gridDesign = gridDesign;
            assessment.certificationSearch = certificationSearch;
        end
    end
end
