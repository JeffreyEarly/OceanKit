classdef IMEigenvalueProblem
    % Describe a canonical scalar eigenvalue problem.
    %
    % `IMEigenvalueProblem` describes the canonical Sturm-Liouville EVP
    %
    % $$
    % -\frac{\partial}{\partial z}
    % \left(p(z)\frac{\partial u}{\partial z}\right)
    % +q(z)u(z)=\lambda r(z)u(z).
    % $$
    %
    % Boundary conditions are written in the same eigenvalue-dependent
    % canonical form:
    %
    % $$
    % -\left[a_i u(z_i)
    % -b_i p(z_i)\frac{\partial u}{\partial z}(z_i)\right]
    % =\lambda\left[c_i u(z_i)
    % -d_i p(z_i)\frac{\partial u}{\partial z}(z_i)\right],
    % \qquad i\in\{\mathrm{bottom},\mathrm{surface}\}.
    % $$
    %
    % `IMEigenvalueProblem` stores this continuous scalar problem; solvers
    % discretize it and return an `IMBasisSet` of retained modes.
    %
    % For example, the Dirichlet problem
    %
    % $$
    % -\frac{\partial^2 u}{\partial z^2}=\lambda u,\qquad u(-1)=u(0)=0
    % $$
    %
    % is created with:
    %
    % ```matlab
    % evp = IMEigenvalueProblem(zDomain=[-1 0], ...
    %     p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
    %     r=@(z,~) ones(size(z)), ...
    %     surfaceBoundary=IMBoundaryCondition.dirichlet(), ...
    %     bottomBoundary=IMBoundaryCondition.dirichlet());
    % ```
    %
    % A solver discretizes this EVP and returns the retained basis set:
    %
    % ```matlab
    % solver = IMSolverSpectral(nEVP=64);
    % basisSet = solver.solveEVP(evp,nModes=4);
    % ```
    %
    % - Topic: Create EVPs
    % - Topic: Define canonical coefficients
    % - Topic: Summarize EVPs
    % - Topic: Inspect EVP configuration
    % - Topic: Inspect inner products
    % - Topic: Inspect definiteness diagnostics
    % - Topic: Inspect mode selection
    % - Topic: Developer topics
    % - Declaration: classdef IMEigenvalueProblem

    properties (SetAccess = private)
        % Short EVP name.
        %
        % `name` appears in solver diagnostics and generated error messages.
        %
        % - Topic: Inspect EVP configuration
        name = "canonical"

        % Physical vertical domain.
        %
        % `zDomain` is sorted as `[bottom surface]` and defines the
        % interval on which the canonical EVP is posed.
        %
        % - Topic: Define canonical coefficients
        zDomain = [-1 0]

        % Coefficient multiplying the derivative flux.
        %
        % `p` defines the flux term in $$-(p u')'$$. It may be a scalar,
        % a vector on the solver grid, or a function handle with signature
        % `values = p(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        p = @(z,~) ones(size(z))

        % Coefficient multiplying the solved variable on the left side.
        %
        % `q` defines the multiplication term in $$q u$$. It may be a
        % scalar, a vector on the solver grid, or a function handle with
        % signature `values = q(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        q = @(z,~) zeros(size(z))

        % Metric coefficient multiplying the eigenvalue side.
        %
        % `r` defines the interior metric in $$\lambda r u$$ and in the
        % default scalar inner product. It may be a scalar, a vector on the
        % solver grid, or a function handle with signature
        % `values = r(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        r = @(z,~) ones(size(z))

        % Surface boundary condition.
        %
        % `surfaceBoundary` stores the scalar boundary condition at
        % `zDomain(2)` in canonical `IMBoundaryCondition` form.
        %
        % - Topic: Define canonical coefficients
        surfaceBoundary = IMBoundaryCondition.dirichlet()

        % Bottom boundary condition.
        %
        % `bottomBoundary` stores the scalar boundary condition at
        % `zDomain(1)` in canonical `IMBoundaryCondition` form.
        %
        % - Topic: Define canonical coefficients
        bottomBoundary = IMBoundaryCondition.dirichlet()

        % Additional coefficient parameters.
        %
        % Fields are copied into the coefficient context so custom
        % coefficient functions can read named values without new public
        % properties. `parameters` is not the whole coefficient context; it
        % is merged into the context provided by the solver and EVP, such
        % as `zDomain`. Subclasses and factories may add or reserve
        % additional fields; check their documentation before reusing field
        % names.
        %
        % ```matlab
        % evp = IMEigenvalueProblem( ...
        %     p=@(z,ctx) ctx.alpha*ones(size(z)), ...
        %     parameters=struct("alpha",2));
        % ```
        %
        % - Topic: Inspect EVP configuration
        parameters = struct()
    end

    methods
        function self = IMEigenvalueProblem(options)
            % Create a canonical scalar EVP.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.p: derivative-flux coefficient
            % - Parameter options.q: left-side value coefficient
            % - Parameter options.r: eigenvalue-side metric coefficient
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Parameter options.parameters: named coefficient parameters
            % - Returns evp: canonical EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "canonical"
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.parameters struct = struct()
            end

            self.name = string(options.name);
            self.zDomain = sort(options.zDomain);
            self.p = options.p;
            self.q = options.q;
            self.r = options.r;
            self.surfaceBoundary = options.surfaceBoundary;
            self.bottomBoundary = options.bottomBoundary;
            self.parameters = options.parameters;
        end

        summarize(self, solver)

        function [A, B] = assemble(self, solver)
            % Build the canonical matrix pair on a solver grid.
            %
            % Interior rows discretize
            %
            % $$
            % -\frac{\partial}{\partial z}\left(p(z)\frac{\partial u}{\partial z}\right)
            % +q(z)u(z)=\lambda r(z)u(z).
            % $$
            %
            % The surface and bottom rows are replaced by the endpoint
            % conditions using endpoint values of `p`, producing the matrix
            % pencil $$A q = \lambda B q$$. This method is mainly for solver
            % implementations, diagnostics, and external eigensolver
            % experiments; ordinary workflows call `solver.solveEVP`.
            %
            % - Topic: Developer topics
            % - Declaration: [A,B] = assemble(evp,solver)
            % - Parameter solver: canonical EVP solver
            % - Returns A: left matrix
            % - Returns B: right matrix
            % - Developer: true
            arguments
                self IMEigenvalueProblem
                solver IMSolver
            end

            solver = solver.configuredForEVP(self);
            context = self.contextForSolver(solver);
            z = solver.zNative(:);
            [pValues, qValues, rValues] = self.coefficientValues(z, context);
            pzValues = solver.differentiateGridValues(pValues, 1);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            D2 = solver.physicalDerivativeMatrix(2);

            A = -diag(pValues)*D2 - diag(pzValues)*D1 + diag(qValues)*D0;
            B = diag(rValues)*D0;
            [A, B] = self.applyBoundaryRow(A, B, solver, "surface", self.surfaceBoundary, pValues);
            [A, B] = self.applyBoundaryRow(A, B, solver, "bottom", self.bottomBoundary, pValues);
        end

        function context = contextForSolver(self, solver)
            % Return the coefficient context for this EVP and solver.
            %
            % The context begins with `solver.context()`, adds `zDomain`,
            % then copies each field of `parameters`. Coefficient handles such
            % as `p(z,ctx)`, `q(z,ctx)`, and `r(z,ctx)` receive this struct.
            %
            % - Topic: Developer topics
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            % - Developer: true
            arguments
                self IMEigenvalueProblem
                solver IMSolver
            end

            context = solver.context();
            context.zDomain = self.zDomain;
            parameterFields = fieldnames(self.parameters);
            for iField = 1:numel(parameterFields)
                fieldName = parameterFields{iField};
                context.(fieldName) = self.parameters.(fieldName);
            end
        end

        function spec = innerProduct(self)
            % Return the scalar inner-product recipe.
            %
            % The canonical solved-variable inner product is
            %
            % $$
            % \langle u_i,u_j\rangle_\mu=
            % \int_{z_b}^{z_s}r(z)u_i(z)u_j(z)\,dz+
            % \sum_{\ell\in S}D_\ell^{-1}L_\ell[u_i]L_\ell[u_j],
            % $$
            %
            % where $$S$$ is the set of active, nondegenerate endpoints,
            % $$D_\ell=\sigma_\ell(a_\ell d_\ell-b_\ell c_\ell)$$, and
            %
            % $$
            % L_\ell[u]=c_\ell u(z_\ell)-d_\ell p(z_\ell)\frac{\partial u}{\partial z}(z_\ell).
            % $$
            %
            % The returned struct has fields `variable`, `interiorWeight`,
            % `surfaceWeights`, and `bottomWeights`. `interiorWeight`
            % stores $$r(z)$$. The endpoint arrays are the same endpoint
            % metric terms returned by `endpointWeights("surface")` and
            % `endpointWeights("bottom")`; each endpoint weight stores
            % `coefficient`, equal to $$D_\ell^{-1}$$, plus `c`, `d`,
            % and `location`.
            %
            % - Topic: Inspect inner products
            % - Declaration: spec = innerProduct(evp)
            % - Returns spec: struct with interior and endpoint metric terms
            arguments
                self IMEigenvalueProblem
            end

            spec.variable = "u";
            spec.interiorWeight = self.r;
            spec.surfaceWeights = self.endpointWeights("surface");
            spec.bottomWeights = self.endpointWeights("bottom");
        end

        function weights = endpointWeights(self, location)
            % Return endpoint metric terms implied by active conditions.
            %
            % For an active boundary condition
            %
            % $$
            % -[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')].
            % $$
            %
            % Yassin's indexing uses `z_1` for the bottom and `z_2` for the
            % surface, so
            %
            % $$
            % D_i=(-1)^{i+1}(a_i d_i-b_i c_i).
            % $$
            %
            % Each returned struct has fields `location`, `coefficient`,
            % `c`, and `d`, representing the endpoint metric contribution
            %
            % $$
            % D_i^{-1}(c_i u-d_i p u_z)^2.
            % $$
            %
            % The `coefficient` field is the scalar returned by
            % `boundary.endpointWeightCoefficient(location)`, while `c`
            % and `d` are copied from the boundary condition properties.
            % In other words, `IMBoundaryCondition` derives the
            % boundary-level scalar coefficient $$D_i^{-1}$$, and this
            % method packages that coefficient with `location`, `c`, and
            % `d` so basis sets can assemble the endpoint Gram-matrix
            % term. These are the endpoint-only terms used by the full
            % `innerProduct()` recipe: `endpointWeights("surface")`
            % appears as `innerProduct().surfaceWeights`, and
            % `endpointWeights("bottom")` appears as
            % `innerProduct().bottomWeights`.
            %
            % - Topic: Inspect inner products
            % - Declaration: weights = endpointWeights(evp,location)
            % - Parameter location: `"surface"`, `"bottom"`, or omitted for both endpoints
            % - Returns weights: endpoint metric terms
            arguments
                self IMEigenvalueProblem
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom", "all"])} = "all"
            end

            switch string(location)
                case "surface"
                    weights = self.weightForBoundary("surface", self.surfaceBoundary);
                case "bottom"
                    weights = self.weightForBoundary("bottom", self.bottomBoundary);
                case "all"
                    weights = [self.weightForBoundary("surface", self.surfaceBoundary); ...
                        self.weightForBoundary("bottom", self.bottomBoundary)];
            end
        end

        function value = negativeEndpointWeightCount(self, options)
            % Count negative endpoint metric weights.
            %
            % `negativeEndpointWeightCount` is the number of active
            % endpoint metric weights with negative coefficient, after the
            % supplied sign tolerance. The canonical scalar EVP has two
            % endpoints with at most one active metric weight per endpoint,
            % so this value is always `0`, `1`, or `2`.
            %
            % - Topic: Inspect inner products
            % - Declaration: value = negativeEndpointWeightCount(evp,options)
            % - Parameter options.tolerance: scalar sign tolerance
            % - Returns value: number of active endpoint terms with negative metric weight
            arguments
                self IMEigenvalueProblem
                options.tolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 0
            end

            weights = self.endpointWeights();
            value = 0;
            for iWeight = 1:numel(weights)
                if weights(iWeight).coefficient < -options.tolerance
                    value = value + 1;
                end
            end
        end

        function diagnostics = definitenessDiagnostics(self, solver)
            % Assess the norm and energy signs that control negative modes.
            %
            % On the solver grid, the canonical EVP has a norm determined
            % by the interior weight `r` and endpoint weights, and an energy
            % determined by the interior `p`, `q` terms plus endpoint
            % contributions. This diagnostic checks the sampled coefficients
            % and endpoint terms that determine whether the norm is positive
            % and whether the energy is nonnegative. These are grid-level
            % checks for the discretized problem, not continuum guarantees
            % between grid points.
            %
            % The returned struct includes sampled minima `pMin`, `qMin`,
            % `rMin`; sign flags `pPositive`, `qNonnegative`, and
            % `rPositive`; norm fields `endpointWeights`,
            % `negativeEndpointWeightCount`, `metricPositive`, and
            % `hasDegenerateEndpointMetric`; energy fields
            % `endpointNumeratorNegativeDirections`,
            % `endpointNumeratorNonnegative`, `interiorNonnegative`, and
            % `quadraticFormNonnegative`; and status fields
            % `assessmentLevel` and `reason`.
            %
            % - Topic: Inspect definiteness diagnostics
            % - Declaration: diagnostics = definitenessDiagnostics(evp,solver)
            % - Parameter solver: canonical EVP solver
            % - Returns diagnostics: struct with norm and energy checks
            arguments
                self IMEigenvalueProblem
                solver IMSolver
            end

            solver = solver.configuredForEVP(self);
            context = self.contextForSolver(solver);
            z = solver.zNative(:);
            [pValues, qValues, rValues] = self.coefficientValues(z, context);
            tolerance = 100*eps;
            pTol = IMEigenvalueProblem.signTolerance(pValues, tolerance);
            qTol = IMEigenvalueProblem.signTolerance(qValues, tolerance);
            rTol = IMEigenvalueProblem.signTolerance(rValues, tolerance);

            diagnostics.pMin = min(pValues);
            diagnostics.qMin = min(qValues);
            diagnostics.rMin = min(rValues);
            diagnostics.pPositive = all(isfinite(pValues)) && diagnostics.pMin > pTol;
            diagnostics.qNonnegative = all(isfinite(qValues)) && diagnostics.qMin >= -qTol;
            diagnostics.rPositive = all(isfinite(rValues)) && diagnostics.rMin > rTol;
            diagnostics.endpointWeights = self.endpointWeights();
            diagnostics.negativeEndpointWeightCount = self.negativeEndpointWeightCount(tolerance=0);
            diagnostics.metricPositive = diagnostics.rPositive && diagnostics.negativeEndpointWeightCount == 0;
            diagnostics.hasDegenerateEndpointMetric = self.hasDegenerateEndpointMetric();
            diagnostics.endpointNumeratorNegativeDirections = self.endpointNegativeDirections();
            diagnostics.endpointNumeratorNonnegative = diagnostics.endpointNumeratorNegativeDirections == 0;
            diagnostics.interiorNonnegative = diagnostics.pPositive && diagnostics.qNonnegative;
            diagnostics.quadraticFormNonnegative = diagnostics.interiorNonnegative && diagnostics.endpointNumeratorNonnegative;
            diagnostics.assessmentLevel = "grid";
            diagnostics.reason = "Grid-level norm and energy signs were checked.";
            if diagnostics.hasDegenerateEndpointMetric
                diagnostics.assessmentLevel = "unknown";
                diagnostics.reason = "An active endpoint determinant is numerically degenerate.";
            elseif ~(diagnostics.pPositive && diagnostics.rPositive && all(isfinite(qValues)))
                diagnostics.assessmentLevel = "unknown";
                diagnostics.reason = "One or more coefficient samples are nonfinite or fail the required signs.";
            end
        end

        function bounds = negativeEigenvalueBounds(self, solver, A)
            % Bound negative eigenvalues using a grid-level assessment.
            %
            % The returned counts describe how many negative eigenvalues are
            % supported by the discretized canonical problem, rather than by
            % raw negative finite-real eigenvalues alone. Exact negative
            % counts require the zero mode to be absent. The returned struct
            % includes `assessmentLevel`, `negativeEndpointWeightCount`,
            % `zeroModeStatus`, `minNegativeEigenvalueCount`,
            % `maxNegativeEigenvalueCount`, and `reason`.
            % `maxNegativeEigenvalueCount` may be the string `"unknown"`
            % when coefficient signs or endpoint determinants cannot be
            % assessed on the grid.
            %
            % - Topic: Inspect mode selection
            % - Declaration: bounds = negativeEigenvalueBounds(evp,solver,A)
            % - Parameter solver: canonical EVP solver
            % - Parameter A: assembled left matrix, used for the zero-mode check
            % - Returns bounds: struct with min/max counts and a reason
            arguments
                self IMEigenvalueProblem
                solver IMSolver
                A (:,:) double = []
            end

            definiteness = self.definitenessDiagnostics(solver);
            zeroMode = self.zeroModeAssessment(A);
            bounds.assessmentLevel = definiteness.assessmentLevel;
            bounds.negativeEndpointWeightCount = definiteness.negativeEndpointWeightCount;
            bounds.zeroModeStatus = zeroMode.zeroModeStatus;
            bounds.minNegativeEigenvalueCount = 0;
            bounds.maxNegativeEigenvalueCount = "unknown";
            bounds.reason = definiteness.reason;

            if definiteness.assessmentLevel == "unknown"
                return;
            end

            if definiteness.metricPositive && definiteness.quadraticFormNonnegative
                bounds.maxNegativeEigenvalueCount = 0;
                bounds.reason = "The grid-level norm is positive and the energy is nonnegative.";
                return;
            end

            if definiteness.quadraticFormNonnegative && definiteness.negativeEndpointWeightCount > 0
                if bounds.zeroModeStatus == "absent"
                    bounds.minNegativeEigenvalueCount = definiteness.negativeEndpointWeightCount;
                    bounds.maxNegativeEigenvalueCount = definiteness.negativeEndpointWeightCount;
                    bounds.reason = "The energy is nonnegative, the endpoint norm has an assessed negative-weight count, and zero is absent.";
                else
                    bounds.maxNegativeEigenvalueCount = definiteness.negativeEndpointWeightCount;
                    bounds.reason = "The negative endpoint weight count bounds the search, but exact negative count requires zero to be absent.";
                end
                return;
            end

            if definiteness.metricPositive && definiteness.interiorNonnegative
                bounds.maxNegativeEigenvalueCount = definiteness.endpointNumeratorNegativeDirections;
                bounds.reason = "The norm is positive and only endpoint energy terms can make the energy negative.";
                return;
            end
        end

        function diagnostics = modeSelectionDiagnostics(self, solver, A)
            % Summarize negative and zero mode selection.
            %
            % Negative modes are bounded by `negativeEigenvalueBounds`.
            % Zero modes are inferred from the assembled left matrix `A`:
            % `zeroModeStatus` is `"present"` when the smallest singular
            % value is zero to the rank-style tolerance
            % $$\sigma_{\min}(A)\le \max(\mathrm{size}(A))\epsilon(\sigma_{\max}(A)),$$
            % `"absent"` when it is larger, and `"unchecked"` when `A` is
            % omitted. Mode labels are ordered as
            % $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
            % The returned struct includes `assessmentLevel`,
            % `negativeEndpointWeightCount`, `minNegativeEigenvalueCount`,
            % `maxNegativeEigenvalueCount`, `zeroModeStatus`,
            % `zeroModeCount`, `zeroModeSingularValue`,
            % `zeroModeTolerance`, and `reason`.
            %
            % - Topic: Inspect mode selection
            % - Declaration: diagnostics = modeSelectionDiagnostics(evp,solver,A)
            % - Parameter solver: canonical EVP solver
            % - Parameter A: assembled left matrix
            % - Returns diagnostics: struct with negative and zero mode selection fields
            arguments
                self IMEigenvalueProblem
                solver IMSolver
                A (:,:) double = []
            end

            bounds = self.negativeEigenvalueBounds(solver, A);
            zeroMode = self.zeroModeAssessment(A);
            diagnostics.assessmentLevel = bounds.assessmentLevel;
            diagnostics.negativeEndpointWeightCount = bounds.negativeEndpointWeightCount;
            diagnostics.minNegativeEigenvalueCount = bounds.minNegativeEigenvalueCount;
            diagnostics.maxNegativeEigenvalueCount = bounds.maxNegativeEigenvalueCount;
            diagnostics.zeroModeStatus = zeroMode.zeroModeStatus;
            diagnostics.zeroModeCount = zeroMode.zeroModeCount;
            diagnostics.zeroModeSingularValue = zeroMode.zeroModeSingularValue;
            diagnostics.zeroModeTolerance = zeroMode.zeroModeTolerance;
            diagnostics.reason = bounds.reason + " Zero mode status is " + zeroMode.zeroModeStatus + ".";
        end

        function selection = selectModes(self, eigenvalues, nModes, solver, A)
            % Select and label retained finite-real eigenmodes.
            %
            % Mode-selection diagnostics decide when raw negative discrete
            % eigenvalues should be retained and whether a zero mode should
            % be included. Retained modes are labeled in the order
            % $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
            % The full diagnostics struct is stored in
            % `selection.modeSelectionDiagnostics`.
            %
            % - Topic: Developer topics
            % - Declaration: selection = selectModes(evp,eigenvalues,nModes,solver,A)
            % - Parameter eigenvalues: finite real candidate eigenvalues
            % - Parameter nModes: number of retained modes
            % - Parameter solver: canonical solver
            % - Parameter A: assembled left matrix
            % - Returns selection: selected indices and mode numbers
            % - Developer: true
            arguments
                self IMEigenvalueProblem
                eigenvalues (:,1) double {mustBeReal, mustBeFinite}
                nModes (1,1) double {mustBeInteger, mustBePositive}
                solver IMSolver
                A (:,:) double
            end

            tolerance = self.eigenvalueTolerance(eigenvalues);
            diagnostics = self.modeSelectionDiagnostics(solver, A);
            negativeCount = nnz(eigenvalues < -tolerance);
            if isnumeric(diagnostics.maxNegativeEigenvalueCount)
                negativeCount = min(negativeCount, diagnostics.maxNegativeEigenvalueCount);
            end

            negativeCandidates = find(eigenvalues < -tolerance);
            [~, negativeSort] = sort(eigenvalues(negativeCandidates), "ascend");
            negativeIndex = negativeCandidates(negativeSort(1:min(negativeCount,numel(negativeSort))));

            zeroIndex = zeros(0,1);
            if diagnostics.zeroModeStatus == "present"
                zeroCandidates = find(abs(eigenvalues) <= tolerance);
                if ~isempty(zeroCandidates)
                    [~, zeroSort] = sort(abs(eigenvalues(zeroCandidates)), "ascend");
                    zeroIndex = zeroCandidates(zeroSort(1));
                end
            end

            positiveCandidates = find(eigenvalues > tolerance);
            positiveCandidates = setdiff(positiveCandidates, [negativeIndex(:); zeroIndex(:)], "stable");
            [~, positiveSort] = sort(eigenvalues(positiveCandidates), "ascend");
            positiveIndex = positiveCandidates(positiveSort);

            sortIndex = [negativeIndex(:); zeroIndex(:); positiveIndex(:)];
            sortIndex = sortIndex(1:min(nModes,numel(sortIndex)));
            modeNumber = zeros(1,numel(sortIndex));
            nNegative = nnz(ismember(sortIndex, negativeIndex));
            modeNumber(1:nNegative) = -1:-1:-nNegative;
            nextIndex = nNegative + 1;
            if ~isempty(zeroIndex) && nextIndex <= numel(modeNumber) && sortIndex(nextIndex) == zeroIndex
                modeNumber(nextIndex) = 0;
                nextIndex = nextIndex + 1;
            end
            modeNumber(nextIndex:end) = 1:(numel(modeNumber) - nextIndex + 1);

            selection.sortIndex = sortIndex;
            selection.modeNumber = modeNumber;
            selection.modeSelectionDiagnostics = diagnostics;
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, modeSelectionDiagnostics)
            % Create the solved scalar basis set for this EVP.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,modeSelectionDiagnostics)
            % - Returns basisSet: solved scalar basis set
            % - Developer: true
            arguments
                self IMEigenvalueProblem
                solver IMSolver
                nativeModes (:,:) double
                eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                modeNumber (1,:) double {mustBeInteger}
                modeSelectionDiagnostics struct
            end

            basisSet = IMBasisSet(solver=solver, evp=self, nativeModes=nativeModes, eigenvalues=eigenvalues, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics);
        end
    end

    methods (Hidden)
        function mask = finiteGeneralizedEigenpairMask(~,eigenvectors,metricMatrix)
            % Reject numerical representations of infinite pencil modes.
            %
            % A finite generalized eigenpair must have a resolvable action
            % under the pencil metric. Candidates satisfying
            %
            % $$
            % \frac{\lVert Bv\rVert_2}{\lVert B\rVert_2\lVert v\rVert_2}
            % \leq 10^3\epsilon
            % $$
            %
            % are numerical representations of the singular pencil's
            % infinite eigenmodes and are removed before physical sorting
            % and mode labeling. Regular pencils such as $$B=I$$ retain all
            % nonzero eigenvectors.
            %
            % - Topic: Developer topics
            % - Declaration: mask = finiteGeneralizedEigenpairMask(evp,eigenvectors,metricMatrix)
            % - Parameter eigenvectors: finite-real generalized eigenvectors
            % - Parameter metricMatrix: assembled generalized metric matrix
            % - Returns mask: logical row mask of family-valid candidates
            % - Developer: true
            metricScale = norm(metricMatrix,2);
            vectorScale = vecnorm(eigenvectors,2,1);
            relativeMetricAction = vecnorm(metricMatrix*eigenvectors,2,1) ...
                ./max(metricScale*vectorScale,realmin);
            mask = relativeMetricAction > 1e3*eps;
        end
    end

    methods (Access = protected)
        function [pValues, qValues, rValues] = coefficientValues(self, z, context)
            arguments
                self IMEigenvalueProblem
                z (:,1) double {mustBeReal}
                context struct
            end

            pValues = IMEigenvalueProblem.evaluateCoefficient(self.p, z, context);
            qValues = IMEigenvalueProblem.evaluateCoefficient(self.q, z, context);
            rValues = IMEigenvalueProblem.evaluateCoefficient(self.r, z, context);
            pValues = IMEigenvalueProblem.expandCoefficient(pValues, z, "p");
            qValues = IMEigenvalueProblem.expandCoefficient(qValues, z, "q");
            rValues = IMEigenvalueProblem.expandCoefficient(rValues, z, "r");
        end
    end

    methods (Access = private)
        function [A, B] = applyBoundaryRow(self, A, B, solver, location, boundary, pValues)
            arguments
                self IMEigenvalueProblem
                A (:,:) double
                B (:,:) double
                solver IMSolver
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
                boundary (1,1) IMBoundaryCondition
                pValues (:,1) double
            end

            index = solver.boundaryIndex(location);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            pEndpoint = pValues(index);
            A(index,:) = -boundary.a*D0(index,:) + boundary.b*pEndpoint*D1(index,:);
            B(index,:) = boundary.c*D0(index,:) - boundary.d*pEndpoint*D1(index,:);
        end

        function weight = weightForBoundary(self, location, boundary)
            arguments
                self IMEigenvalueProblem
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
                boundary (1,1) IMBoundaryCondition
            end

            weight = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
            if ~boundary.isEigenvalueDependent()
                return;
            end
            coefficient = boundary.endpointWeightCoefficient(location);
            if ~isfinite(coefficient)
                return;
            end
            weight(1,1).location = string(location);
            weight(1,1).coefficient = coefficient;
            weight(1,1).c = boundary.c;
            weight(1,1).d = boundary.d;
        end

        function count = endpointNegativeDirections(self)
            arguments
                self IMEigenvalueProblem
            end

            endpoints = [
                struct("location", "surface", "boundary", self.surfaceBoundary)
                struct("location", "bottom", "boundary", self.bottomBoundary)
            ];
            count = 0;
            for iEndpoint = 1:numel(endpoints)
                boundary = endpoints(iEndpoint).boundary;
                location = endpoints(iEndpoint).location;
                if boundary.isEigenvalueDependent()
                    if self.isDegenerateEndpointMetric(location, boundary)
                        continue;
                    end
                    H = boundary.endpointNumeratorMatrix(location);
                    count = count + nnz(eig(0.5*(H + H.')) < -sqrt(eps)*max(1,norm(H,"fro")));
                elseif boundary.b ~= 0
                    beta = boundary.robinEnergyCoefficient(location);
                    count = count + double(beta < -sqrt(eps)*max(1,abs(beta)));
                end
            end
        end

        function tf = hasDegenerateEndpointMetric(self)
            arguments
                self IMEigenvalueProblem
            end

            tf = self.isDegenerateEndpointMetric("surface", self.surfaceBoundary) ...
                || self.isDegenerateEndpointMetric("bottom", self.bottomBoundary);
        end

        function tf = isDegenerateEndpointMetric(self, location, boundary)
            arguments
                self IMEigenvalueProblem
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
                boundary (1,1) IMBoundaryCondition
            end

            tf = boundary.isEigenvalueDependent() && ~isfinite(boundary.endpointWeightCoefficient(location));
        end

        function zeroMode = zeroModeAssessment(self, A)
            arguments
                self IMEigenvalueProblem
                A (:,:) double = []
            end

            zeroMode.zeroModeStatus = "unchecked";
            zeroMode.zeroModeCount = 0;
            zeroMode.zeroModeSingularValue = NaN;
            zeroMode.zeroModeTolerance = NaN;
            if isempty(A)
                return;
            end
            singularValues = svd(A);
            zeroMode.zeroModeSingularValue = min(singularValues);
            zeroMode.zeroModeTolerance = max(size(A))*eps(max(singularValues));
            zeroMode.zeroModeCount = nnz(singularValues <= zeroMode.zeroModeTolerance);
            if zeroMode.zeroModeCount > 0
                zeroMode.zeroModeStatus = "present";
            else
                zeroMode.zeroModeStatus = "absent";
            end
        end

        function tolerance = eigenvalueTolerance(self, eigenvalues)
            arguments
                self IMEigenvalueProblem
                eigenvalues (:,1) double {mustBeReal, mustBeFinite}
            end

            finiteScale = abs(eigenvalues(isfinite(eigenvalues)));
            if isempty(finiteScale)
                scale = 1;
            else
                scale = max(1,min(finiteScale(finiteScale > 0), [], "omitnan"));
                if isempty(scale) || ~isfinite(scale)
                    scale = max(1,max(finiteScale));
                end
            end
            tolerance = 1e-8*scale;
        end
    end

    methods (Static)
        function values = evaluateCoefficient(coefficient, z, context)
            % Evaluate a scalar, vector, or coefficient function.
            %
            % - Topic: Developer topics
            % - Declaration: values = IMEigenvalueProblem.evaluateCoefficient(coefficient,z,context)
            % - Developer: true
            arguments
                coefficient
                z double
                context struct
            end

            if isa(coefficient, "function_handle")
                try
                    values = coefficient(z, context);
                catch
                    try
                        values = coefficient(z);
                    catch
                        values = coefficient(context);
                    end
                end
            else
                values = coefficient;
            end
        end
    end

    methods (Static, Access = private)
        function values = expandCoefficient(values, z, name)
            arguments
                values {mustBeNumeric}
                z (:,1) double
                name {mustBeTextScalar}
            end

            if isscalar(values)
                values = values*ones(size(z));
            end
            values = values(:);
            if length(values) ~= length(z)
                error("IMEigenvalueProblem:InvalidCoefficientSize", ...
                    "Coefficient %s must evaluate to a scalar or one value per grid point.", name);
            end
        end

        function tolerance = signTolerance(values, tau)
            arguments
                values {mustBeNumeric}
                tau (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
            end

            tolerance = tau*max(1,max(abs(values(:))));
        end

    end
end
