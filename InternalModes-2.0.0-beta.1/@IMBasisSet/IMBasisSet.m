classdef IMBasisSet
    % Store solved scalar canonical EVP modes.
    %
    % `IMBasisSet` stores the scalar modes selected from a canonical EVP and
    % evaluates the solved variable `u` and derivative `uz`. Modal
    % normalization is applied lazily when values, Gram matrices, or spectra
    % are requested. For each retained mode,
    % $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
    % where the scale factor $$s_j$$ is supplied by a basis-set
    % normalization rule. Custom rules are added after solving with
    % `addNormalization`.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet = basisSet.addNormalization("constantScaled", @(basisSet,j) C*basisSet.innerProductNormFactor(j));
    % basisSet.normalization = "unity";
    % u = basisSet.u(z);
    % factors = basisSet.normalizationFactors("unity");
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate modes
    % - Topic: Normalize modes
    % - Topic: Analyze modes
    % - Topic: Build discrete transforms
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSet

    properties
        % Active normalization rule name.
        %
        % This string selects a rule in the basis-set normalization
        % registry. Create custom rules with `addNormalization`.
        % The selected rule returns the scale factor $$s_j$$ used by `u`,
        % `uz`, and Gram-matrix methods. Passing `normalization=name` to an
        % evaluation method overrides this property for that call.
        %
        % - Topic: Normalize modes
        normalization = "unity"
    end

    properties (SetAccess = protected)
        % Solver that created the native modes.
        %
        % This is the discretization object used to interpolate native
        % coefficient columns. It is mainly useful for diagnostics and
        % developer workflows.
        %
        % - Topic: Developer topics — Native representation
        % - Developer: true
        solver

        % EVP descriptor that was solved.
        %
        % - Topic: Inspect basis sets
        evp

        % Native mode columns.
        %
        % These are the columns returned by the numerical solver before
        % interpolation onto physical coordinates and before modal
        % normalization. Most workflows should call `u(z)` or `uz(z)`
        % instead.
        %
        % - Topic: Developer topics — Native representation
        % - Developer: true
        nativeModes

        % Retained eigenvalues.
        %
        % The retained values are the discrete eigenvalues
        % $$\lambda_j$$ from the assembled canonical problem.
        %
        % - Topic: Inspect basis sets
        eigenvalues

        % Retained-mode labels.
        %
        % `modeNumber` is a row vector with one integer label per retained
        % mode column. `modeNumber(j)` labels `eigenvalues(j)` and column
        % `j` of `u(z)`, `uz(z)`, and, for internal-mode basis sets, `F(z)`
        % and `G(z)`.
        %
        % Numerical solves label retained modes in the order
        % $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
        % Negative labels are ordinal labels for retained negative
        % eigenvalues sorted by eigenvalue order; they do not identify
        % surface, bottom, or any other fixed physical branch. The label
        % `0` marks an inferred zero, barotropic, or null mode when one is
        % retained. Positive labels mark ordinary positive/interior modes.
        %
        % ```matlab
        % basisSet.modeNumber
        % % [-1 0 1 2]
        % ```
        %
        % See `modeSelectionDiagnostics` for the diagnostics that explain
        % why negative or zero modes were retained.
        %
        % - Topic: Inspect basis sets
        modeNumber

        % Mode-selection diagnostics.
        %
        % This is the diagnostics struct returned by
        % `evp.modeSelectionDiagnostics` when the solver selected and
        % labeled retained modes. Numerical solves use it to record
        % negative-mode bounds and zero-mode status.
        %
        % - Topic: Inspect basis sets
        modeSelectionDiagnostics

        % Additional metadata.
        %
        % `metadata` stores creation or diagnostic information associated
        % with this solved basis set. It is not interpreted by `IMBasisSet`;
        % EVP coefficient parameters remain available through
        % `basisSet.evp.parameters`.
        %
        % - Topic: Inspect basis sets
        metadata

        % Physical vertical domain.
        %
        % - Topic: Inspect basis sets
        zDomain
    end

    properties (Access = private)
        % Map from normalization names to rule handles.
        normalizationNameMap
        % Optional whole-family normalization rule handles.
        normalizationVectorRuleMap
    end

    methods
        function self = IMBasisSet(options)
            % Create a solved scalar basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: canonical EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.modeNumber: retained-mode labels
            % - Parameter options.modeSelectionDiagnostics: mode-selection diagnostics
            % - Parameter options.normalization: active normalization rule name
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: solved scalar basis set
            arguments
                options.solver IMSolver
                options.evp IMEigenvalueProblem
                options.nativeModes (:,:) double
                options.eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
            end

            nModes = size(options.nativeModes,2);
            if length(options.eigenvalues) ~= nModes
                error("IMBasisSet:InvalidEigenvalueCount", "eigenvalues must contain one value for each native mode column.");
            end
            if length(options.modeNumber) ~= nModes
                error("IMBasisSet:InvalidModeNumber", "modeNumber must contain one value for each native mode column.");
            end
            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.modeSelectionDiagnostics = options.modeSelectionDiagnostics;
            if isempty(options.normalization)
                self.normalization = "unity";
            else
                self.normalization = options.normalization;
            end
            self.metadata = options.metadata;
            self.zDomain = self.evp.zDomain;
            self.normalizationNameMap = configureDictionary("string","cell");
            self.normalizationVectorRuleMap = configureDictionary("string","cell");
            self = self.addVectorNormalization("unity",@(basisSet,iMode) basisSet.innerProductNormFactor(iMode),@(basisSet) basisSet.innerProductNormFactors());
        end

        function self = addNormalization(self, name, rule)
            % Add a named normalization rule.
            %
            % `addNormalization` registers or replaces one rule on this
            % basis set. The rule has signature
            % `scale = rule(basisSet,iMode)` and returns the raw scale
            % factor $$s_j$$ for retained mode `iMode`. Evaluated modes use
            % $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
            %
            % For a constant-scaled norm,
            % $$s_j=C\|u_j\|_\mu,$$
            % use:
            %
            % ```matlab
            % C = 2;
            % basisSet = basisSet.addNormalization("constantScaled", @(basisSet,j) C*basisSet.innerProductNormFactor(j));
            % basisSet.normalization = "constantScaled";
            % ```
            %
            % For an eigenvalue-scaled norm,
            % $$s_j=\sqrt{|\lambda_j|}\|u_j\|_\mu,$$
            % use:
            %
            % ```matlab
            % basisSet = basisSet.addNormalization("eigenvalueScaled", @(basisSet,j) sqrt(abs(basisSet.eigenvalues(j)))*basisSet.innerProductNormFactor(j));
            % ```
            %
            % If `name` already exists, the rule is overwritten.
            %
            % - Topic: Normalize modes
            % - Declaration: basisSet = addNormalization(basisSet,name,rule)
            % - Parameter name: normalization rule name
            % - Parameter rule: function handle with signature `scale = rule(basisSet,iMode)`
            % - Returns basisSet: basis set with the rule installed
            arguments
                self IMBasisSet
                name {mustBeTextScalar}
                rule (1,1) function_handle
            end

            name = self.normalizationName(name);
            self.normalizationNameMap{name} = rule;
            self.normalizationVectorRuleMap{name} = [];
        end

        function names = normalizationNames(self)
            % Return installed normalization rule names.
            %
            % `normalizationNames` reports the rules available to
            % `normalizationFactors` and selectable by
            % `basisSet.normalization`.
            %
            % - Topic: Normalize modes
            % - Declaration: names = normalizationNames(basisSet)
            % - Returns names: string array of installed normalization names
            arguments
                self IMBasisSet
            end

            names = sort(string(keys(self.normalizationNameMap)));
        end

        function values = u(self, z, options)
            % Evaluate solved scalar modes.
            %
            % Returned columns are normalized as
            % $$u_j(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
            % where $$s_j$$ comes from `normalizationFactors`.
            %
            % - Topic: Evaluate modes
            % - Declaration: values = u(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: scalar mode values
            arguments
                self IMBasisSet
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            values = self.rawU(z) ./ self.normalizationFactors(options.normalization);
        end

        function values = uz(self, z, options)
            % Evaluate solved scalar vertical derivatives.
            %
            % Derivatives are scaled by the same modal factors used for
            % `u`, so $$u'_j(z)=u_j^{\mathrm{raw}\prime}(z)/s_j$$.
            %
            % - Topic: Evaluate modes
            % - Declaration: values = uz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: derivative mode values
            arguments
                self IMBasisSet
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            values = self.rawUz(z) ./ self.normalizationFactors(options.normalization);
        end

        function factors = normalizationFactors(self, normalization)
            % Return scale factors for a normalization rule.
            %
            % For a requested rule name this returns the row vector
            % $$s_j$$ used by evaluation methods:
            % $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
            % The default scalar `"unity"` rule is installed on every
            % scalar basis set and uses
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
            % product. Custom rules are added to the basis set with
            % `addNormalization`, and the active rule is selected by
            % `basisSet.normalization`.
            %
            % ```matlab
            % factors = basisSet.normalizationFactors("unity");
            % uUnity = basisSet.u(z,normalization="unity");
            % ```
            %
            % - Topic: Normalize modes
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of positive scale factors
            arguments
                self IMBasisSet
                normalization = self.normalization
            end

            name = self.normalizationName(normalization);
            if ~isKey(self.normalizationNameMap, name)
                error("IMBasisSet:UnsupportedNormalization", "The basis set does not define a ""%s"" normalization.", name);
            end
            normalizeMode = self.normalizationNameMap{name};
            nModes = size(self.nativeModes,2);
            vectorRule = self.normalizationVectorRuleMap{name};
            if isempty(vectorRule)
                factors = zeros(1,nModes);
                for iMode = 1:nModes
                    factors(iMode) = normalizeMode(self,iMode);
                end
            else
                factors = vectorRule(self);
                if ~isnumeric(factors) || ~isreal(factors) || ~isvector(factors) || numel(factors) ~= nModes
                    error("IMBasisSet:InvalidNormalizationRule", "The vector normalization rule ""%s"" must return one real numeric factor per retained mode.",name);
                end
                factors = reshape(factors,1,[]);
            end
            factors = abs(factors);
            factors(factors == 0 | ~isfinite(factors)) = 1;
        end

        function gram = gramMatrix(self, options)
            % Return a scalar Gram matrix.
            %
            % With no arguments this uses the full basis-set domain. Passing
            % `zBounds=[zMin zMax]` restricts the interior integral to that
            % interval and includes endpoint terms only when the interval
            % contains the corresponding physical endpoint. For normalized
            % scalar modes,
            %
            % $$
            % M_{ij}(z_a,z_b)=
            % \int_{z_a}^{z_b}r(z)u_i(z)u_j(z)\,dz+
            % \sum_{\ell\in S,\ z_\ell\in[z_a,z_b]}D_\ell^{-1}L_\ell[u_i]L_\ell[u_j],
            % $$
            %
            % where included endpoint terms use
            %
            % $$
            % L_\ell[u_j]=c_\ell u_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial u_j}{\partial z}(z_\ell).
            % $$
            %
            % Use `endpointGramTerms` to inspect the prepared endpoint
            % vectors that generate the rank-one endpoint updates.
            %
            % - Topic: Analyze modes
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns gram: scalar Gram matrix
            arguments
                self IMBasisSet
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            z = self.solver.innerProductGrid(options.zBounds);
            context = self.evp.contextForSolver(self.solver);
            normalizationFactors = self.normalizationFactors(self.normalization);
            values = self.rawU(z) ./ normalizationFactors;
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.solver.integrateInnerProduct(z, integrand, options.zBounds);
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end

            endpointTerms = self.endpointGramTerms(zBounds=options.zBounds, normalization=self.normalization, useNormalized=true);
            for iTerm = 1:numel(endpointTerms)
                valuesAtEndpoint = endpointTerms(iTerm).values(:);
                gram = gram + endpointTerms(iTerm).coefficient*(valuesAtEndpoint*valuesAtEndpoint.');
            end
        end

        function terms = endpointGramTerms(self, options)
            % Prepare rank-one endpoint terms for scalar Gram matrices.
            %
            % `endpointGramTerms` returns the endpoint pieces that enter the
            % scalar Gram matrix. Each returned element contains a full
            % row vector over retained modes, not one pairwise matrix
            % contribution. Canonical endpoint weights use
            %
            % $$
            % L_\ell[u_j]=c_\ell u_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial u_j}{\partial z}(z_\ell),
            % $$
            %
            % and the Gram matrix applies the rank-one update
            %
            % $$
            % M \leftarrow M+D_\ell^{-1}L_\ell L_\ell^\mathsf{T}.
            % $$
            %
            % Endpoint terms are omitted when `zBounds` does not include
            % that endpoint.
            %
            % - Topic: Developer topics — Gram-matrix assembly
            % - Declaration: terms = endpointGramTerms(basisSet,options)
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Parameter options.normalization: normalization rule name
            % - Parameter options.useNormalized: whether returned values use the normalization
            % - Returns terms: struct array with `location`, `coefficient`, `values`, and `kind`
            % - Developer: true
            arguments
                self IMBasisSet
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
                options.normalization = self.normalization
                options.useNormalized (1,1) logical = true
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end

            terms = struct("location", {}, "coefficient", {}, "values", {}, "kind", {});
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct();
            endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
            normalizationFactors = [];
            for iWeight = 1:numel(endpointWeights)
                endpointWeight = endpointWeights(iWeight);
                if string(endpointWeight.location) == "surface"
                    zEndpoint = self.zDomain(2);
                else
                    zEndpoint = self.zDomain(1);
                end

                tolerance = 100*eps(max(1,max(abs([options.zBounds(:); zEndpoint]))));
                endpointIsIncluded = abs(min(options.zBounds) - zEndpoint) <= tolerance || abs(max(options.zBounds) - zEndpoint) <= tolerance;
                if ~endpointIsIncluded
                    continue;
                end

                uEndpoint = self.rawU(zEndpoint);
                uzEndpoint = self.rawUz(zEndpoint);
                if options.useNormalized
                    if isempty(normalizationFactors)
                        normalizationFactors = self.normalizationFactors(options.normalization);
                    end
                    uEndpoint = uEndpoint ./ normalizationFactors;
                    uzEndpoint = uzEndpoint ./ normalizationFactors;
                end

                pEndpoint = IMEigenvalueProblem.evaluateCoefficient(self.evp.p, zEndpoint, context);
                values = endpointWeight.c*uEndpoint - endpointWeight.d*pEndpoint*uzEndpoint;
                terms(end+1) = struct("location", string(endpointWeight.location), "coefficient", endpointWeight.coefficient, "values", values, "kind", "endpointWeight");
            end
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial scalar Gram matrix.
            %
            % This computes the eigendecomposition of the symmetric Gram
            % matrix on `zBounds` and sorts window-mode eigenvalues from
            % largest to smallest.
            %
            % - Topic: Analyze modes
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMBasisSet
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            gram = self.gramMatrix(zBounds=options.zBounds);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients)
            % Compute a scalar modal spectrum.
            %
            % For modal coefficients $$c_j$$ this returns
            % $$S_j=M_{jj}|c_j|^2,$$ where $$M$$ is the full-domain scalar
            % Gram matrix.
            %
            % - Topic: Analyze modes
            % - Declaration: spectrum = spectrum(basisSet,coefficients)
            % - Parameter coefficients: modal coefficients
            % - Returns spectrum: modal spectrum
            arguments
                self IMBasisSet
                coefficients (:,1) double
            end

            spectrum = self.crossSpectrum(coefficients, coefficients);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB)
            % Compute a scalar modal cross-spectrum.
            %
            % For modal coefficient vectors $$a_j$$ and $$b_j$$ this
            % returns $$S_j=M_{jj}\operatorname{Re}(a_j b_j^*)$$.
            %
            % - Topic: Analyze modes
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMBasisSet
                coefficientsA (:,1) double {mustBeFinite}
                coefficientsB (:,1) double {mustBeFinite}
            end

            nModes = size(self.nativeModes,2);
            if length(coefficientsA) ~= nModes || length(coefficientsB) ~= nModes
                error("IMBasisSet:InvalidCoefficientCount", "Coefficient vectors must contain one value for each retained mode (%d).", nModes);
            end
            gram = self.gramMatrix();
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end

        [transform, assessment] = discreteTransform(self, options)
        [transform, assessment] = certifiedDiscreteTransform(self, options)
        [transform, assessment] = fitDiscreteTransform(self, options)
        [weights, weightFit] = quadratureWeightsForPoints(self, options)
        [z, gridDesign] = pointsFromModeRoots(self, options)
        [z, gridDesign] = modeRootGrid(self, options)

    end

    methods
        function self = orientModeSigns(self)
            % Orient scalar modes with a deterministic sign convention.
            %
            % This developer utility flips native mode columns so each
            % scalar mode has a deterministic sign based on its largest
            % amplitude on the solver integration grid. `IMBasisSet` is a
            % value class, so callers must keep the returned basis set:
            %
            % ```matlab
            % basisSet = basisSet.orientModeSigns();
            % ```
            %
            % - Topic: Developer topics — Native representation
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with oriented native mode signs
            % - Developer: true
            if isempty(self.nativeModes)
                return;
            end
            z = self.solver.innerProductGrid(self.zDomain);
            values = self.rawU(z);
            signs = ones(1,size(values,2));
            for iMode = 1:size(values,2)
                [~, index] = max(abs(values(:,iMode)));
                if values(index,iMode) < 0
                    signs(iMode) = -1;
                end
            end
            self.nativeModes = self.nativeModes .* signs;
        end

        function factor = innerProductNormFactor(self, iMode)
            % Return the scalar inner-product norm factor.
            %
            % This developer utility returns the raw scale factor used by
            % the default `"unity"` normalization rule before any modal
            % normalization has been applied:
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle_\mu|}.$$
            % The scalar inner product includes the interior weight and
            % any canonical endpoint terms prepared by `endpointGramTerms`.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw scalar inner-product scale factor
            % - Developer: true
            arguments
                self IMBasisSet
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            z = self.solver.innerProductGrid(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            values = self.rawU(z);
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            integrand = weight(:).*values(:,iMode).*values(:,iMode);
            value = self.solver.integrateInnerProduct(z, integrand, self.zDomain);
            endpointTerms = self.endpointGramTerms(zBounds=self.zDomain, useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                valuesAtEndpoint = endpointTerms(iTerm).values;
                value = value + endpointTerms(iTerm).coefficient*valuesAtEndpoint(iMode)^2;
            end
            factor = sqrt(abs(value));
        end

        function factors = innerProductNormFactors(self)
            % Return raw scalar inner-product factors for the whole family.
            z = self.solver.innerProductGrid(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            values = self.rawU(z);
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            else
                weight = weight(:);
            end
            normSquared = zeros(1,size(values,2));
            for iMode = 1:size(values,2)
                normSquared(iMode) = self.solver.integrateInnerProduct(z,weight.*values(:,iMode).*values(:,iMode),self.zDomain);
            end
            endpointTerms = self.endpointGramTerms(zBounds=self.zDomain,useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                endpointValues = endpointTerms(iTerm).values;
                normSquared = normSquared+endpointTerms(iTerm).coefficient*endpointValues.^2;
            end
            factors = sqrt(abs(normSquared));
        end

        function factor = maxAmplitudeNormFactor(self, iMode)
            % Return the maximum scalar amplitude.
            %
            % This developer utility returns the raw maximum-amplitude
            % scale
            % $$s_j=\max_z |u_j^{\mathrm{raw}}(z)|$$
            % on the basis-set integration grid. Normalization rules can
            % use this factor to make the largest scalar mode amplitude
            % equal to one.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = maxAmplitudeNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: maximum raw scalar amplitude
            % - Developer: true
            arguments
                self IMBasisSet
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            z = self.solver.innerProductGrid(self.zDomain);
            values = self.rawU(z);
            factor = max(abs(values(:,iMode)));
        end

        function values = rawU(self, z)
            % Evaluate raw solved scalar modes.
            %
            % This developer utility evaluates the unnormalized native
            % mode columns on the physical grid `z`. The public `u` method
            % applies the active normalization rule after this step:
            % $$u_j(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
            % where $$s_j$$ comes from `normalizationFactors`.
            %
            % - Topic: Developer topics — Native representation
            % - Declaration: values = rawU(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized scalar mode values
            % - Developer: true
            arguments
                self IMBasisSet
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.solver.evaluateNativeModes(self.nativeModes, z);
        end

        function values = rawUz(self, z)
            % Evaluate raw solved scalar derivatives.
            %
            % This developer utility evaluates the unnormalized physical
            % $$z$$-derivative of the solved native modes. The public `uz`
            % method applies the active normalization rule after this
            % step:
            % $$\frac{\partial u_j}{\partial z}(z)=
            % \frac{\partial u_j^{\mathrm{raw}}}{\partial z}(z)/s_j,$$
            % where $$s_j$$ comes from `normalizationFactors`.
            %
            % - Topic: Developer topics — Native representation
            % - Declaration: values = rawUz(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized scalar derivative values
            % - Developer: true
            arguments
                self IMBasisSet
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
        end
    end

    methods (Access = protected)
        function self = addVectorNormalization(self,name,rule,vectorRule)
            arguments
                self IMBasisSet
                name {mustBeTextScalar}
                rule (1,1) function_handle
                vectorRule (1,1) function_handle
            end
            name = self.normalizationName(name);
            self.normalizationNameMap{name} = rule;
            self.normalizationVectorRuleMap{name} = vectorRule;
        end

        function transform = buildDiscreteTransform(self, z, weights, nModes)
            arguments
                self IMBasisSet
                z (:,1) double {mustBeReal, mustBeFinite}
                weights (:,1) double {mustBeReal, mustBeFinite}
                nModes (1,1) double {mustBeInteger, mustBePositive}
            end

            z = z(:);
            weights = weights(:);
            if length(z) ~= length(weights)
                error("IMBasisSet:InvalidDiscreteWeights", "weights must contain one value for each sample point in z.");
            end
            if nModes > size(self.nativeModes,2)
                error("IMBasisSet:InvalidDiscreteModeCount", "The basis set contains %d modes, but nModes=%d was requested.", size(self.nativeModes,2), nModes);
            end
            if length(z) < nModes
                error("IMBasisSet:InsufficientDiscreteSamples", "At least %d sample points are required for %d retained modes.", nModes, nModes);
            end
            if length(z) < 2 || any(diff(z) <= 0)
                error("IMBasisSet:InvalidDiscreteGrid", "z must contain at least two strictly increasing, unique sample points.");
            end
            domainTolerance = 100*eps(max(1,max(abs(self.zDomain))));
            if any(z < self.zDomain(1) - domainTolerance) || any(z > self.zDomain(2) + domainTolerance)
                error("IMBasisSet:InvalidDiscreteGrid", "All sample points must lie inside the basis-set zDomain.");
            end
            if ~any(weights ~= 0)
                error("IMBasisSet:InvalidDiscreteWeights", "weights must contain at least one nonzero value.");
            end

            inverseMatrix = self.u(z);
            inverseMatrix = inverseMatrix(:,1:nModes);
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct();
            interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(interiorWeight)
                interiorWeight = interiorWeight*ones(size(z));
            elseif numel(interiorWeight) ~= length(z)
                error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must return one value for each sample point.");
            else
                interiorWeight = interiorWeight(:);
            end
            if ~isreal(interiorWeight) || any(~isfinite(interiorWeight))
                error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must be finite and real at every sample point.");
            end

            metricMatrix = diag(interiorWeight.*weights);
            if isfield(spec, "endpointInnerProductTerms") && ~isempty(spec.endpointInnerProductTerms)
                error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Cross-variable endpoint inner products are not supported by scalar discrete transforms.");
            end
            endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
            for iWeight = 1:numel(endpointWeights)
                endpointWeight = endpointWeights(iWeight);
                if endpointWeight.d ~= 0
                    error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Endpoint derivative traces cannot be determined from arbitrary point samples; scalar discrete transforms require endpoint terms that depend only on sampled values.");
                end
                if endpointWeight.c == 0
                    continue;
                end
                if string(endpointWeight.location) == "surface"
                    zEndpoint = self.zDomain(2);
                else
                    zEndpoint = self.zDomain(1);
                end
                endpointIndex = find(abs(z - zEndpoint) <= domainTolerance, 1);
                if isempty(endpointIndex)
                    error("IMBasisSet:MissingDiscreteEndpointSample", "The %s endpoint must be included in z to represent its value-only inner-product term.", string(endpointWeight.location));
                end
                if ~isfinite(endpointWeight.coefficient)
                    error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "The %s endpoint inner-product coefficient is not finite.", string(endpointWeight.location));
                end
                metricMatrix(endpointIndex,endpointIndex) = metricMatrix(endpointIndex,endpointIndex) + endpointWeight.coefficient*endpointWeight.c*endpointWeight.c;
            end

            continuousGramMatrix = self.gramMatrix();
            targetGramMatrix = diag(diag(continuousGramMatrix(1:nModes,1:nModes)));
            transform = IMDiscreteTransform(z=z,weights=weights,modeNumber=self.modeNumber(1:nModes),normalization=self.normalizationName(self.normalization),inverseMatrix=inverseMatrix,metricMatrix=metricMatrix,targetGramMatrix=targetGramMatrix);
        end

        function name = normalizationName(~, normalization)
            parts = split(string(normalization), ".");
            name = parts(end);
        end

        function unsupported(~, operationName)
            error("IMBasisSet:UnsupportedOperation", "IMBasisSet does not support %s for this basis.", operationName);
        end
    end

end
