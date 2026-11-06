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
    % rmsFactors = basisSet.normalizationFactors("rms");
    % basisSet = basisSet.addNormalization("twiceRms",@(~,j) 2*rmsFactors(j));
    % basisSet.normalization = "twiceRms";
    % u = basisSet.u(z);
    % factors = basisSet.normalizationFactors();
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
        normalization = "rms"
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

    properties (Access = protected)
        % Map from normalization names to rule handles.
        normalizationNameMap
        % Optional whole-family normalization rule handles.
        normalizationVectorRuleMap
        normalizationCache_
        normalizationBuiltins_ (1,:) string = strings(1,0)
    end

    properties (Hidden, SetAccess = protected)
        % Immutable construction identities; derived normalization stays separate.
        basisContractState struct = struct()
    end

    properties (Dependent, SetAccess = private)
        % Supported continuous operations and conditional preparation requirements.
        % - Topic: Inspect basis sets
        capabilities
        % Requested/resolved scientific identity, without evaluating normalization.
        % - Topic: Inspect basis sets
        scientificSpecification
        % Original backend request and current selection, with unverified accuracy.
        % - Topic: Inspect basis sets
        constructionRecord
        % Fixed, per-column or independent horizontal identity.
        % - Topic: Inspect basis sets
        horizontalKind
        k
        kSquared
        % Signs of the active normalization form.
        % - Topic: Normalize modes
        normSignatures
        % Number of retained columns, including all special sectors.
        nColumns
        % Scientific family name and retained column identities.
        family
        columnLabels
        modeInfo
    end

    methods
        projection = projectionOnGrid(self,z,weights,options)
        transform = scalarTransform(self,z,weights,options)
    end

    methods (Hidden)
        function value = isBuiltinNormalization(self,name)
            value = ismember(string(name),self.normalizationBuiltins_);
        end
        function self = withConstructionRequest(self,request,solver,selection)
            % Attach original creation evidence without evaluating the modes.
            arguments
                self
                request (1,1) struct
                solver
                selection (1,1) struct = struct()
            end
            state = self.basisContractState;
            state.request = request;
            state.requestedSolver = IMBasisContracts.solverRecord(solver,false);
            if isfield(selection,"sortIndex"), state.sourceEigenpairIndices = reshape(selection.sortIndex,1,[]); end
            if isfield(selection,"modeSelectionDiagnostics"), state.selectionDiagnostics = selection.modeSelectionDiagnostics; end
            self.basisContractState = state;
        end
        recipe = projectionRecipe(self,options)
    end

    methods
        function value = get.scientificSpecification(self)
            value = IMBasisContracts.scientificSpecification(self);
        end
        function value = get.constructionRecord(self)
            value = IMBasisContracts.constructionRecord(self);
        end
        function value = get.horizontalKind(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.kind;
        end
        function value = get.k(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.k;
        end
        function value = get.kSquared(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.kSquared;
        end
        function value = get.capabilities(self)
            value = IMBasisContracts.capabilities(self);
        end
        function value = get.normSignatures(self)
            value = IMBasisContracts.normalizationSignatures(self);
        end

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
                self.normalization = "rms";
            else
                self.normalization = options.normalization;
            end
            self.metadata = options.metadata;
            self.zDomain = self.evp.zDomain;
            self.normalizationNameMap = configureDictionary("string","cell");
            self.normalizationVectorRuleMap = configureDictionary("string","cell");
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            if ~isa(self,"IMInternalModesBasis")
                self = self.addVectorNormalization("rms",@(basis,j) IMInternalModesNormalizationTools.selectedNormFactor(basis,"rms",j),@(basis) basis.rmsNormFactors());
            end
            self.normalizationBuiltins_ = self.normalizationNames();
            if ~isa(self,"IMInternalModesBasis"), IMBasisContracts.validateNormalization(self,self.normalization); end
            self.basisContractState = IMBasisContracts.initialState(self);
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
            % To set a dimensionless scalar profile RMS of one half, use
            % a positive multiplier of its existing RMS divisor:
            % $$s_j=2\sqrt{D^{-1}\int |u_j^{\mathrm{raw}}|^2\,dz}.$$
            %
            % ```matlab
            % rmsFactors = basisSet.normalizationFactors("rms");
            % basisSet = basisSet.addNormalization("twiceRms",@(~,j) 2*rmsFactors(j));
            % basisSet.normalization = "twiceRms";
            % ```
            %
            % Each custom rule must return one numeric scalar that is real,
            % positive, and finite. Validation occurs when the rule is used.
            %
            % Existing custom names may be replaced; built-in names are reserved.
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

            name = IMInternalModesNormalizationTools.customName(name,"IMBasisSet");
            self.normalizationNameMap{name} = rule;
            self.normalizationVectorRuleMap{name} = [];
            self = self.resetNormalizationCache();
        end

        function value = get.nColumns(self)
            value = size(self.nativeModes,2);
        end

        function value = get.family(self)
            value = string(self.evp.name);
            if isfield(self.metadata,"family"), value = string(self.metadata.family); end
        end

        function value = get.modeInfo(self)
            if isa(self.evp,"IMInternalModes")
                value = IMInternalModesNormalizationTools.modeInfo(self.evp,self.modeNumber,self.eigenvalues,self.metadata,self.modeSelectionDiagnostics);
            else
                value = table(string(self.modeNumber(:)),repmat("mode",self.nColumns,1),self.modeNumber(:),self.eigenvalues(:),VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            end
        end

        function value = get.columnLabels(self)
            value = self.modeInfo.columnLabel.';
        end

        function self = set.normalization(self,name)
            name = self.normalizationName(name);
            if ~isempty(self.normalizationNameMap) && ~isKey(self.normalizationNameMap,name)
                error("IMBasisSet:UnsupportedNormalization","The basis does not define normalization %s.",name);
            end
            self.normalization = name;
        end

        function self = withNormalization(self,name)
            % Select a registered normalization without solving or evaluating it.
            name = self.normalizationName(name);
            if ~isKey(self.normalizationNameMap,name)
                error("IMBasisSet:UnsupportedNormalization","The basis does not define normalization %s.",name);
            end
            self.normalization = name;
        end

        function self = selectColumns(self,columns,options)
            % Keep the concrete basis type, scientific identities and coupled fields.
            arguments
                self IMBasisSet
                columns (1,:) double {mustBeInteger,mustBePositive}
                options.modeInfo table = table()
            end
            if any(columns > self.nColumns) || numel(unique(columns)) ~= numel(columns)
                error("IMBasisSet:InvalidColumns","Select distinct existing columns.");
            end
            source = self;
            self.basisContractState = IMBasisContracts.selectState(self.basisContractState,columns);
            info = self.modeInfo(columns,:);
            self = self.selectAlignedColumns(columns);
            self.nativeModes = self.nativeModes(:,columns);
            self.eigenvalues = self.eigenvalues(columns);
            self.modeNumber = self.modeNumber(columns);
            if width(options.modeInfo) > 0, info = options.modeInfo; end
            if height(info) ~= numel(columns)
                error("IMBasisSet:InvalidModeInfo","modeInfo must contain one row per selected column.");
            end
            self.metadata.modeInfo = info;
            self.modeSelectionDiagnostics.modeInfo = info;
            if isfield(self.modeSelectionDiagnostics,"geostrophicInverseDepth")
                self.modeSelectionDiagnostics.geostrophicInverseDepth = self.modeSelectionDiagnostics.geostrophicInverseDepth(columns);
            end
            self = self.resetNormalizationCache();
            for name = reshape(self.normalizationNames(),1,[])
                cached = source.normalizationCache_(name);
                if ~isempty(cached), self.normalizationCache_(name,cached(columns)); end
                if isempty(self.normalizationVectorRuleMap{name})
                    rule = self.normalizationNameMap{name};
                    self.normalizationNameMap{name} = @(~,j) rule(source,columns(j));
                elseif ~ismember(name,self.normalizationBuiltins_)
                    self.normalizationNameMap{name} = @(~,j) IMInternalModesNormalizationTools.selectedNormFactor(source,name,columns(j));
                    self.normalizationVectorRuleMap{name} = @(~) IMInternalModesNormalizationTools.columns(source.normalizationFactors(name),columns);
                end
            end
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

        function values = evaluate(self,z,options)
            % Evaluate the solved scalar channel and its first derivative.
            %
            % Use variable="u" and derivativeOrder=0 or 1. Per-call
            % normalization leaves the active basis normalization unchanged.
            %
            % - Topic: Evaluate modes
            % - Declaration: values = evaluate(self,z,options)
            % - Parameter z: physical coordinate column within zDomain
            % - Parameter options.variable: required solved channel u
            % - Parameter options.derivativeOrder: zero for values or one for the vertical derivative
            % - Parameter options.normalization: normalization rule for this call
            % - Returns values: sample-by-mode scalar profiles or derivatives
            arguments (Input)
                self IMBasisSet
                z (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.normalization = self.normalization
            end
            arguments (Output)
                values double
            end
            if options.variable == ""
                error("IMBasisSet:MissingVariable","Specify the solved channel with variable=""u"".");
            end
            if options.variable ~= "u"
                error("IMBasisSet:UnsupportedVariable","The solved scalar channel is u.");
            end
            if any(z < self.zDomain(1)) || any(z > self.zDomain(2))
                error("IMBasisSet:OutsideDomain","Evaluation points must lie within zDomain.");
            end
            if options.derivativeOrder == 0
                values = self.u(z,normalization=options.normalization);
            elseif options.derivativeOrder == 1
                values = self.uz(z,normalization=options.normalization);
            else
                error("IMBasisSet:UnsupportedDerivative","The solved scalar channel supports derivative orders zero and one.");
            end
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
            % The default scalar `"rms"` rule is installed on every
            % scalar basis set and uses the dimensionless profile convention
            % $$s_j=\sqrt{D^{-1}\int |u_j|^2\,dz},\quad D=z_s-z_b.$$ Custom rules are added to the basis set with
            % `addNormalization`, and the active rule is selected by
            % `basisSet.normalization`.
            %
            % ```matlab
            % factors = basisSet.normalizationFactors("rms");
            % uRms = basisSet.u(z,normalization="rms");
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
            if ~isempty(vectorRule)
                cached = self.normalizationCache_(name);
                if ~isempty(cached) || nModes == 0
                    factors = cached;
                    if nModes == 0, factors = zeros(1,0); end
                    return
                end
            end
            if isempty(vectorRule)
                factors = zeros(1,nModes);
                for iMode = 1:nModes
                    value = normalizeMode(self,iMode);
                    IMInternalModesNormalizationTools.validateCustomFactor(value,"IMBasisSet");
                    factors(iMode) = value;
                end
            else
                factors = vectorRule(self);
                if ~isnumeric(factors) || ~isreal(factors) || ~isvector(factors) || numel(factors) ~= nModes
                    error("IMBasisSet:InvalidNormalizationRule", "The vector normalization rule ""%s"" must return one real numeric factor per retained mode.",name);
                end
                factors = reshape(factors,1,[]);
            end
            if ~isreal(factors) || any(factors <= 0 | ~isfinite(factors))
                error("IMBasisSet:NullNormalization","Every retained mode requires a finite positive normalization divisor.");
            end
            if ~isempty(vectorRule), self.normalizationCache_(name,factors); end
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
            gram = self.scalarGramMatrices(options.zBounds);
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
                % The built-in spectral value-only endpoint needs no flux.
                valueOnly = any(strcmp(class(self),{'IMBasisSet','IMInternalModesBasis'})) && strcmp(class(self.solver),'IMSolverSpectral') && endpointWeight.d == 0;
                if valueOnly
                    uzEndpoint = zeros(size(uEndpoint));
                else
                    uzEndpoint = self.rawUz(zEndpoint);
                end
                if options.useNormalized
                    if isempty(normalizationFactors)
                        normalizationFactors = self.normalizationFactors(options.normalization);
                    end
                    uEndpoint = uEndpoint ./ normalizationFactors;
                    uzEndpoint = uzEndpoint ./ normalizationFactors;
                end

                if endpointWeight.d == 0
                    % A value-only trace has no flux coefficient to sample.
                    % In particular, 0*(f0^2/N2) is undefined at neutral endpoints.
                    values = endpointWeight.c*uEndpoint;
                else
                    pEndpoint = IMEigenvalueProblem.evaluateCoefficient(self.evp.p, zEndpoint, context);
                    values = endpointWeight.c*uEndpoint - endpointWeight.d*pEndpoint*uzEndpoint;
                end
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
        [z,weights] = referenceQuadrature(self,zBounds,refinement)
        points = referenceBreakpoints(self,zBounds)

        function [gram,majorant,isPositiveMetric] = scalarGramMatrices(self,zBounds)
            % Integrate the represented scalar pairing and its positive majorant.
            [z,weights] = self.referenceQuadrature(zBounds);
            values = self.u(z);
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(weight), weight = repmat(weight,size(z)); end
            if numel(weight) ~= numel(z) || ~isreal(weight) || any(~isfinite(weight(:)))
                error("IMBasisSet:InvalidInnerProductWeight","The scalar interior weight must return one finite real value per reference point.");
            end
            gram = values.'*((weights.*weight(:)).*values);
            majorant = values.'*((weights.*abs(weight(:))).*values);
            isPositiveMetric = all(weight(:) >= 0);
            terms = self.endpointGramTerms(zBounds=zBounds);
            for iTerm = 1:numel(terms)
                endpointValues = terms(iTerm).values(:);
                endpointMatrix = endpointValues*endpointValues.';
                gram = gram+terms(iTerm).coefficient*endpointMatrix;
                majorant = majorant+abs(terms(iTerm).coefficient)*endpointMatrix;
                isPositiveMetric = isPositiveMetric && terms(iTerm).coefficient >= 0;
            end
            gram = 0.5*(gram+gram.');
            majorant = 0.5*(majorant+majorant.');
        end

        function factors = rmsNormFactors(self)
            % Normalize the scalar profile without changing its scientific pairing.
            [z,weights] = self.normalizationQuadrature();
            values = self.rawU(z);
            factors = sqrt(sum(weights.*abs(values).^2,1)/diff(self.zDomain));
        end

        function [z,weights] = normalizationQuadrature(self)
            % Integrate the profiles returned by the public evaluator.
            [z,weights] = self.referenceQuadrature(self.zDomain);
        end

        function self = resetNormalizationCache(self)
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
        end

        function self = selectAlignedColumns(self,~)
            % Family subclasses select companion coefficient arrays here.
        end

        function self = addVectorNormalization(self,name,rule,vectorRule)
            arguments
                self IMBasisSet
                name {mustBeTextScalar}
                rule (1,1) function_handle
                vectorRule (1,1) function_handle
            end
            name = self.normalizationName(name);
            self.normalizationBuiltins_ = setdiff(self.normalizationBuiltins_,string(name),"stable");
            self.normalizationNameMap{name} = rule;
            self.normalizationVectorRuleMap{name} = vectorRule;
            self = self.resetNormalizationCache();
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

            [continuousGramMatrix,continuousMajorant] = self.scalarGramMatrices(self.zDomain);
            targetGramMatrix = continuousGramMatrix(1:nModes,1:nModes);
            majorantGramMatrix = continuousMajorant(1:nModes,1:nModes);
            physical = IMProjection(inverseMatrix,metricMatrix,targetGramMatrix,majorantGramMatrix=majorantGramMatrix,columnLabels=string(self.modeNumber(1:nModes)));
            transform = IMScalarDiscreteTransform(physical,physical.recovery,z=z,weights=weights,modeNumber=self.modeNumber(1:nModes),normalization=self.normalizationName(self.normalization));
        end

        function name = normalizationName(~, normalization)
            name = IMInternalModesNormalizationTools.name(normalization,"IMBasisSet");
        end

        function unsupported(~, operationName)
            error("IMBasisSet:UnsupportedOperation", "IMBasisSet does not support %s for this basis.", operationName);
        end
    end

    methods (Hidden)
        function factor = innerProductNormFactor(self, iMode)
            % Return the square-root magnitude of the raw scalar pairing.
            %
            % This mathematical utility returns the raw diagonal pairing
            % magnitude before any modal normalization has been applied:
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle_\mu|}.$$
            % The scalar inner product includes the interior weight and
            % any canonical endpoint terms prepared by `endpointGramTerms`.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: square-root magnitude of the weighted scalar pairing
            % - Developer: true
            arguments
                self IMBasisSet
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            [z,weights] = self.referenceQuadrature(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            values = self.rawU(z);
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            integrand = weight(:).*values(:,iMode).*values(:,iMode);
            value = sum(weights.*integrand);
            endpointTerms = self.endpointGramTerms(zBounds=self.zDomain, useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                valuesAtEndpoint = endpointTerms(iTerm).values;
                value = value + endpointTerms(iTerm).coefficient*valuesAtEndpoint(iMode)^2;
            end
            factor = sqrt(abs(value));
        end

        function factors = innerProductNormFactors(self)
            % Return raw scalar inner-product factors for the whole family.
            [z,weights] = self.referenceQuadrature(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            values = self.rawU(z);
            spec = self.evp.innerProduct();
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            else
                weight = weight(:);
            end
            normSquared = sum((weights.*weight).*values.^2,1);
            endpointTerms = self.endpointGramTerms(zBounds=self.zDomain,useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                endpointValues = endpointTerms(iTerm).values;
                normSquared = normSquared+endpointTerms(iTerm).coefficient*endpointValues.^2;
            end
            factors = sqrt(abs(normSquared));
        end

    end

end
