classdef IMAnalyticalInternalModesBasis
    % Store exact internal-mode functions from an analytical solution family.
    %
    % `IMAnalyticalInternalModesBasis` has the same user-facing `F(z)`,
    % `G(z)`, normalization, and Gram-matrix idiom as numerical internal-mode
    % bases, but it evaluates closed-form functions supplied by an
    % `IMAnalyticalSolution`.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % basisSet = solution.internalModes(evp,nModes=4);
    % G = basisSet.G(linspace(-5000,0,128).');
    % ```
    %
    % - Topic: Evaluate analytical modes
    % - Topic: Analyze Gram matrices
    % - Topic: Inspect analytical modes
    % - Topic: Developer topics
    % - Declaration: classdef IMAnalyticalInternalModesBasis

    properties
        % Active normalization rule name.
        %
        % - Topic: Evaluate analytical modes
        normalization = "rmsF"
    end

    properties (SetAccess = private)
        % Analytical solution family that created this basis.
        %
        % - Topic: Inspect analytical modes
        solution

        % Internal-mode EVP represented by these exact functions.
        %
        % - Topic: Inspect analytical modes
        evp

        % Equivalent depths for retained modes.
        %
        % - Topic: Inspect analytical modes
        h

        % Retained eigenvalues.
        %
        % - Topic: Inspect analytical modes
        eigenvalues

        % Retained mode labels.
        %
        % - Topic: Inspect analytical modes
        modeNumber

        % Physical vertical domain.
        %
        % - Topic: Inspect analytical modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect analytical modes
        N2

        % Additional creation metadata.
        %
        % - Topic: Inspect analytical modes
        metadata

        % Mode-selection diagnostics.
        %
        % Analytical bases do not run a discrete mode-selection pass, so this
        % is usually empty.
        %
        % - Topic: Inspect analytical modes
        modeSelectionDiagnostics
    end

    properties (Access = private)
        rawVariableFunction
        rawUzFunction
        orientationSigns
        normalizationNameMap
        normalizationVectorRules_
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
        nColumns
        family
        columnLabels
        modeInfo
        frequency
        capabilities
        normSignatures
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
        function self = IMAnalyticalInternalModesBasis(options)
            % Create an exact internal-mode basis.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: basisSet = IMAnalyticalInternalModesBasis(options)
            % - Parameter options.solution: analytical solution family, or `[]` for profile-independent sectors
            % - Parameter options.evp: internal-mode EVP
            % - Parameter options.h: equivalent depths
            % - Parameter options.modeNumber: retained-mode labels
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.rawVariableFunction: exact `F`/`G` evaluator
            % - Parameter options.rawUzFunction: exact solved-variable derivative evaluator
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: creation metadata
            % - Returns basisSet: exact internal-mode basis
            arguments
                options.solution = []
                options.evp IMInternalModes
                options.h double
                options.modeNumber double
                options.N2 (1,1) function_handle
                options.rawVariableFunction (1,1) function_handle
                options.rawUzFunction (1,1) function_handle
                options.normalization = []
                options.metadata struct = struct()
                options.modeSelectionDiagnostics struct = struct()
            end

            if ~isempty(options.solution) && ~isa(options.solution,"IMAnalyticalSolution")
                error("IMAnalyticalInternalModesBasis:InvalidSolution", ...
                    "solution must be an IMAnalyticalSolution or [] for an exact profile-independent sector.");
            end
            self.solution = options.solution;
            self.evp = options.evp;
            self.h = reshape(options.h,1,[]);
            self.eigenvalues = 1 ./ self.h;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.zDomain = options.evp.zDomain;
            self.N2 = options.N2;
            self.metadata = options.metadata;
            self.modeSelectionDiagnostics = options.modeSelectionDiagnostics;
            self.rawVariableFunction = options.rawVariableFunction;
            self.rawUzFunction = options.rawUzFunction;
            self.orientationSigns = ones(1,length(self.h));
            zOrientation = linspace(self.zDomain(1),self.zDomain(2),max(1025,4*self.nColumns+1)).';
            GValues = self.rawVariable("G",zOrientation);
            FValues = self.rawVariable("F",zOrientation);
            if self.evp.formulation == "G"
                GzSurface = self.rawUzFunction(self,self.zDomain(2));
            else
                FSurface = self.rawVariableFunction(self,"F",self.zDomain(2));
                GzSurface = FSurface.*self.eigenvalues;
            end
            allowFFallback = self.evp.formulation == "F" & self.eigenvalues == 0;
            self.orientationSigns = IMModeOrientationTools.shallowInteriorGPositive( ...
                GValues=GValues,GzSurface=GzSurface,FValues=FValues,depth=diff(self.zDomain), ...
                surfaceIndex=size(GValues,1),allowFFallback=allowFFallback);
            self.metadata.modeOrientation = IMModeOrientationTools.convention;
            self.normalizationNameMap = configureDictionary("string","cell");
            self.normalizationVectorRules_ = configureDictionary("string","cell");
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            self = self.installNormalizationRules();
            if isempty(options.normalization)
                self.normalization = IMInternalModesNormalizationTools.default(self.evp);
            else
                self.normalization = options.normalization;
            end
            self.normalizationBuiltins_ = self.normalizationNames();
            self.basisContractState = IMBasisContracts.initialState(self);
        end

        function self = addNormalization(self, name, rule)
            % Add a named normalization rule.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: basisSet = addNormalization(basisSet,name,rule)
            % - Parameter name: normalization rule name
            % - Parameter rule: function handle with signature `scale = rule(basisSet,iMode)`
            % - Returns basisSet: basis set with the rule installed
            arguments
                self IMAnalyticalInternalModesBasis
                name {mustBeTextScalar}
                rule (1,1) function_handle
            end

            name = IMInternalModesNormalizationTools.customName(name,"IMAnalyticalInternalModesBasis");
            self.normalizationNameMap{name} = rule;
            self.normalizationVectorRules_{name} = [];
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
        end

        function value = get.nColumns(self)
            value = numel(self.h);
        end

        function value = get.family(self)
            value = string(self.evp.name);
            if isfield(self.metadata,"family"), value = string(self.metadata.family); end
        end

        function value = get.modeInfo(self)
            value = IMInternalModesNormalizationTools.modeInfo(self.evp,self.modeNumber,self.eigenvalues,self.metadata,self.modeSelectionDiagnostics);
        end

        function value = get.columnLabels(self)
            value = self.modeInfo.columnLabel.';
        end

        function value = get.frequency(self)
            value = IMInternalModesNormalizationTools.frequency(self.evp,self.h);
        end

        function value = get.capabilities(self)
            value = IMBasisContracts.capabilities(self);
        end

        function value = get.normSignatures(self)
            value = IMBasisContracts.normalizationSignatures(self);
        end

        function self = set.normalization(self,name)
            name = self.normalizationName(name);
            if ~isempty(self.normalizationNameMap) && ~isKey(self.normalizationNameMap,name)
                error("IMAnalyticalInternalModesBasis:UnsupportedNormalization","The basis does not define normalization %s.",name);
            end
            self.normalization = name;
        end

        function self = withNormalization(self,name)
            name = self.normalizationName(name);
            if ~isKey(self.normalizationNameMap,name)
                error("IMAnalyticalInternalModesBasis:UnsupportedNormalization","The analytical basis does not define normalization %s.",name);
            end
            self.normalization = name;
        end

        function self = selectColumns(self,columns,options)
            arguments
                self IMAnalyticalInternalModesBasis
                columns (1,:) double {mustBeInteger,mustBePositive}
                options.modeInfo table = table()
            end
            if any(columns > self.nColumns) || numel(unique(columns)) ~= numel(columns)
                error("IMAnalyticalInternalModesBasis:InvalidColumns","Select distinct existing columns.");
            end
            source = self;
            self.basisContractState = IMBasisContracts.selectState(self.basisContractState,columns);
            info = self.modeInfo(columns,:);
            self.rawVariableFunction = @(~,variable,z) IMInternalModesNormalizationTools.columns(source.rawVariable(variable,z),columns);
            self.rawUzFunction = @(~,z) IMInternalModesNormalizationTools.columns(source.rawUz(z),columns);
            self.orientationSigns = ones(1,numel(columns));
            self.h = self.h(columns);
            self.eigenvalues = self.eigenvalues(columns);
            self.modeNumber = self.modeNumber(columns);
            if width(options.modeInfo) > 0, info = options.modeInfo; end
            if height(info) ~= numel(columns)
                error("IMAnalyticalInternalModesBasis:InvalidModeInfo","modeInfo must have one row per selected column.");
            end
            self.metadata.modeInfo = info;
            self.modeSelectionDiagnostics.modeInfo = info;
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            for name = reshape(self.normalizationNames(),1,[])
                cached = source.normalizationCache_(name);
                if ~isempty(cached), self.normalizationCache_(name,cached(columns)); end
                if isempty(self.normalizationVectorRules_{name})
                    rule = self.normalizationNameMap{name};
                    self.normalizationNameMap{name} = @(~,j) rule(source,columns(j));
                elseif ~ismember(name,self.normalizationBuiltins_)
                    self.normalizationNameMap{name} = @(~,j) IMInternalModesNormalizationTools.selectedNormFactor(source,name,columns(j));
                    self.normalizationVectorRules_{name} = @(~) IMInternalModesNormalizationTools.columns(source.normalizationFactors(name),columns);
                end
            end
        end

        function values = evaluate(self,z,options)
            % Evaluate exact fields and supported physical derivatives.
            %
            % Custom callbacks retain their supplied solved-variable first
            % derivative. Further derivative recovery uses only the defining
            % equations of a known physical family.
            %
            % Use variable="F" or variable="G" and derivativeOrder to select
            % a physical channel. Per-call normalization leaves the basis's
            % active normalization unchanged.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: values = evaluate(self,z,options)
            % - Parameter z: physical coordinate column within zDomain
            % - Parameter options.variable: required physical channel F or G
            % - Parameter options.derivativeOrder: supported nonnegative vertical derivative order, zero by default
            % - Parameter options.normalization: normalization rule for this call
            % - Returns values: sample-by-mode physical profiles or derivatives
            arguments (Input)
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.normalization = self.normalization
            end
            arguments (Output)
                values double
            end
            variable = options.variable;
            derivativeOrder = options.derivativeOrder;
            if variable == ""
                error("IMAnalyticalInternalModesBasis:MissingVariable","Specify the channel with variable= followed by a supported channel name.");
            end
            if ~ismember(variable,["F","G"])
                error("IMAnalyticalInternalModesBasis:UnsupportedVariable","Unsupported channel %s. Inspect capabilities for supported channels and derivatives.",variable);
            end
            IMBasisContracts.validateEvaluation(self,variable,derivativeOrder);
            IMBasisContracts.validateNormalization(self,options.normalization);
            if any(z < self.zDomain(1)) || any(z > self.zDomain(2))
                error("IMAnalyticalInternalModesBasis:OutsideDomain","Evaluation points must lie within zDomain.");
            end
            values = self.rawVariableDerivative(variable,z,derivativeOrder)./self.normalizationFactors(options.normalization);
        end

        function values = Fz(self,z,options)
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            values = self.evaluate(z,variable="F",derivativeOrder=1,normalization=options.normalization);
        end

        function values = Gz(self,z,options)
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            values = self.evaluate(z,variable="G",derivativeOrder=1,normalization=options.normalization);
        end

        function values = Gzz(self,z,options)
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            values = self.evaluate(z,variable="G",derivativeOrder=2,normalization=options.normalization);
        end

        transform = discreteTransform(self,grid,options)

        function names = normalizationNames(self)
            % Return installed normalization rule names.
            %
            % - Topic: Inspect analytical modes
            % - Declaration: names = normalizationNames(basisSet)
            % - Returns names: installed rule names
            names = sort(string(keys(self.normalizationNameMap)));
        end

        function F = F(self, z, options)
            % Evaluate exact `F` modes.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns F: exact `F` values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            F = self.evaluate(z,variable="F",normalization=options.normalization);
        end

        function G = G(self, z, options)
            % Evaluate exact `G` modes.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns G: exact `G` values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            G = self.evaluate(z,variable="G",normalization=options.normalization);
        end

        function values = uz(self, z, options)
            % Evaluate exact solved-variable derivatives.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: values = uz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns values: exact derivative values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            values = self.rawUz(z) ./ self.normalizationFactors(options.normalization);
        end

        function factors = normalizationFactors(self, normalization)
            % Return scale factors for a normalization rule.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization rule name
            % - Returns factors: row vector of scale factors
            arguments
                self IMAnalyticalInternalModesBasis
                normalization = self.normalization
            end

            name = self.normalizationName(normalization);
            if ~isKey(self.normalizationNameMap, name)
                error("IMAnalyticalInternalModesBasis:UnsupportedNormalization", "The analytical basis does not define a ""%s"" normalization.", name);
            end
            rule = self.normalizationNameMap{name};
            vectorRule = self.normalizationVectorRules_{name};
            if ~isempty(vectorRule)
                factors = self.normalizationCache_(name);
                if isempty(factors) && self.nColumns > 0
                    factors = vectorRule(self);
                    if ~isreal(factors) || any(factors <= 0 | ~isfinite(factors))
                        error("IMAnalyticalInternalModesBasis:NullNormalization","Every retained mode requires a finite positive normalization divisor.");
                    end
                    self.normalizationCache_(name,factors);
                elseif self.nColumns == 0
                    factors = zeros(1,0);
                end
                return
            end
            factors = zeros(1,length(self.h));
            for iMode = 1:length(self.h)
                value = rule(self,iMode);
                IMInternalModesNormalizationTools.validateCustomFactor(value,"IMAnalyticalInternalModesBasis");
                factors(iMode) = value;
            end
            if ~isreal(factors) || any(factors <= 0 | ~isfinite(factors))
                error("IMAnalyticalInternalModesBasis:NullNormalization","Every retained mode requires a finite positive normalization divisor.");
            end
        end

        function spec = majorantInnerProduct(self, options)
            % Return the induced positive Hilbert-majorant recipe.
            %
            % The recipe retains the positive interior weight and replaces
            % every signed endpoint coefficient by its absolute value. Use
            % `evp.innerProduct` for the signed Pontryagin recipe used by
            % projection and signed invariants.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spec = majorantInnerProduct(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spec: positive interior and absolute-endpoint recipe
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end
            spec = self.evp.majorantInnerProduct(options.variable);
        end

        function gram = gramMatrix(self, options)
            % Return the signed Gram matrix for exact `F` or `G` modes.
            %
            % For variable $$V$$, endpoint terms are included only when
            % `zBounds` reaches the corresponding physical endpoint:
            % $$M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
            % \sum_\ell \gamma_\ell L_\ell[V_i]L_\ell[V_j]+
            % \sum_\ell \alpha_\ell V_i(z_\ell)V_j(z_\ell).$$
            % This is the Pontryagin pairing used for projection and signed
            % invariants and can be indefinite. Use `majorantGramMatrix`
            % for positive magnitudes and error measures.
            % Use `endpointGramTerms` to inspect the prepared endpoint
            % vectors that generate the rank-one endpoint updates.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns gram: Gram matrix
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateZBounds(options.zBounds);
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
        end

        function gram = majorantGramMatrix(self, options)
            % Return the positive Hilbert-majorant Gram matrix.
            %
            % The majorant retains the positive interior contribution and
            % replaces every signed endpoint coefficient by its absolute
            % value. It coincides with `gramMatrix` when all endpoint
            % coefficients are nonnegative.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = majorantGramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns gram: positive Hilbert-majorant Gram matrix
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateZBounds(options.zBounds);
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true, true);
        end

        function value = majorantNorm(self, coefficients, options)
            % Return the positive Hilbert-majorant norm of coefficients.
            %
            % For coefficient vector $$c$$ this returns
            % $$\sqrt{c^*M_+c}$$. The quantity
            % $$\sqrt{|c^*Mc|}$$ formed from the signed Gram matrix is not
            % a norm for arbitrary states.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: value = majorantNorm(basisSet,coefficients,options)
            % - Parameter coefficients: one coefficient per retained mode
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns value: positive scalar norm
            arguments
                self IMAnalyticalInternalModesBasis
                coefficients (:,1) double {mustBeFinite}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateCoefficientVector(coefficients, "coefficients");
            gram = self.majorantGramMatrix(variable=options.variable,zBounds=options.zBounds);
            valueSquared = real(coefficients(:)'*(gram*coefficients(:)));
            tolerance = 1e3*eps(max(1,norm(gram,2)*norm(coefficients,2)^2));
            if valueSquared < -tolerance
                error("IMAnalyticalInternalModesBasis:InvalidMajorantGramMatrix", "The computed majorant quadratic form is negative beyond roundoff.");
            end
            value = sqrt(max(0,valueSquared));
        end

        function terms = endpointGramTerms(self, options)
            % Prepare rank-one endpoint terms for exact Gram matrices.
            %
            % `endpointGramTerms` returns endpoint vectors over all retained
            % analytical modes. Solved-form endpoint weights use
            % $$L_\ell[V_j]=c_\ell V_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial V_j}{\partial z}(z_\ell),$$
            % and contribute
            % $$M \leftarrow M+\gamma_\ell L_\ell L_\ell^\mathsf{T}.$$
            % Catalog endpoint value terms use $$V_j(z_\ell)$$ and
            % contribute
            % $$M \leftarrow M+\alpha_\ell V_\ell V_\ell^\mathsf{T}.$$
            % Endpoint terms are omitted when `zBounds` does not include
            % that endpoint.
            %
            % - Topic: Developer topics
            % - Declaration: terms = endpointGramTerms(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Parameter options.normalization: normalization rule name
            % - Parameter options.useNormalized: whether returned values use the normalization
            % - Returns terms: struct array with `location`, `coefficient`, `values`, and `kind`
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
                options.normalization = self.normalization
                options.useNormalized (1,1) logical = true
            end

            self.validateZBounds(options.zBounds);
            variable = string(options.variable);
            spec = self.evp.innerProduct(variable);
            IMAnalyticalInternalModesBasis.assertInnerProductAvailable(spec);
            terms = struct("location", {}, "coefficient", {}, "values", {}, "kind", {});
            context = self.context();
            normalizationFactors = [];

            if variable == self.evp.formulation
                endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
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

                    uEndpoint = self.rawVariable(self.evp.formulation, zEndpoint);
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

            for iTerm = 1:numel(spec.endpointInnerProductTerms)
                endpointTerm = spec.endpointInnerProductTerms(iTerm);
                if string(endpointTerm.location) == "surface"
                    zEndpoint = self.zDomain(2);
                else
                    zEndpoint = self.zDomain(1);
                end

                tolerance = 100*eps(max(1,max(abs([options.zBounds(:); zEndpoint]))));
                endpointIsIncluded = abs(min(options.zBounds) - zEndpoint) <= tolerance || abs(max(options.zBounds) - zEndpoint) <= tolerance;
                if ~endpointIsIncluded
                    continue;
                end

                values = self.rawVariable(endpointTerm.variable, zEndpoint);
                if options.useNormalized
                    if isempty(normalizationFactors)
                        normalizationFactors = self.normalizationFactors(options.normalization);
                    end
                    values = values ./ normalizationFactors;
                end
                terms(end+1) = struct("location", string(endpointTerm.location), "coefficient", endpointTerm.coefficient, "values", values, "kind", "endpointInnerProductTerm");
            end
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial-depth Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            gram = self.gramMatrix(variable=options.variable, zBounds=options.zBounds);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute a modal spectrum.
            %
            % This is a signed spectrum and can contain negative entries in
            % negative Pontryagin directions. Use `majorantNorm` for a
            % positive total magnitude.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spectrum: modal spectrum
            arguments
                self IMAnalyticalInternalModesBasis
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a modal cross-spectrum.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMAnalyticalInternalModesBasis
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            self.validateCoefficientVector(coefficientsA, "coefficientsA");
            self.validateCoefficientVector(coefficientsB, "coefficientsB");
            gram = self.gramMatrix(variable=options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end
    end

    methods
        function factor = maxAmplitudeNormFactor(self, iMode, options)
            % Return the maximum amplitude of `F` or `G`.
            %
            % This developer utility returns the raw maximum-amplitude
            % scale
            % $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$
            % for exact analytical modes on the analytical integration
            % grid.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAmplitudeNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: maximum raw variable amplitude
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            values = self.rawVariable(string(options.variable), self.integrationGrid(self.zDomain));
            factor = max(abs(values(:,iMode)));
        end

    end

    methods (Access = private)
        function self = installNormalizationRules(self)
            for name = IMInternalModesNormalizationTools.available(self.evp)
                self.normalizationNameMap{name} = @(basis,j) basis.profileNormFactor(j,name);
                self.normalizationVectorRules_{name} = @(basis) basis.profileNormFactors(name);
            end
        end

        function values = rawVariableDerivative(self,variable,z,order)
            if order == 0
                values = self.rawVariable(variable,z);
            elseif variable == self.evp.formulation && order == 1
                values = self.rawUz(z);
            elseif variable == "G" && order == 1 && IMInternalModesNormalizationTools.hasHydrostaticDerivativeRelations(self.evp)
                values = self.rawVariable("F",z)./self.h;
            elseif variable == "G" && order == 2 && (IMInternalModesNormalizationTools.isWave(self.evp) || IMInternalModesNormalizationTools.hasHydrostaticDerivativeRelations(self.evp))
                potential = 0;
                weight = self.N2(z)/self.evp.g;
                if self.evp.name == "waveModesAtWavenumber"
                    potential = self.evp.parameters.k^2;
                    weight = (self.N2(z)-self.evp.f0^2)/self.evp.g;
                elseif self.evp.name == "waveModesAtFrequency"
                    weight = (self.N2(z)-self.evp.parameters.omega^2)/self.evp.g;
                end
                values = (potential-weight./self.h).*self.rawVariable("G",z);
            elseif variable == "F" && order == 1 && IMInternalModesNormalizationTools.isWave(self.evp)
                values = diff(self.zDomain)*self.rawVariableDerivative("G",z,2);
            elseif variable == "F" && order == 1 && IMInternalModesNormalizationTools.hasHydrostaticDerivativeRelations(self.evp)
                values = -self.N2(z)/self.evp.g.*self.rawVariable("G",z);
            else
                error("IMAnalyticalInternalModesBasis:UnsupportedDerivative","Derivative order %d of %s is unavailable.",order,variable);
            end
        end

        function factor = profileNormFactor(self,j,name)
            factors = self.profileNormFactors(name);
            factor = factors(j);
        end

        function factors = profileNormFactors(self,name)
            [z,w] = self.integrationGrid(self.zDomain);
            if startsWith(name,"maxAbs")
                z = [self.zDomain(1);z;self.zDomain(2)];
                w = [0;w;0];
            end
            D = diff(self.zDomain);
            switch name
                case {"meanDensityAnomalyEnergy","generalizedPotentialEnstrophy"}
                    variable = "G";
                    if name == "generalizedPotentialEnstrophy", variable = "F"; end
                    spec = self.evp.innerProduct(variable);
                    self.assertInnerProductAvailable(spec);
                    values = self.rawVariable(variable,z);
                    weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,self.context());
                    normSquared = sum(w.*weight.*values.^2,1);
                    terms = self.endpointGramTerms(variable=variable,useNormalized=false);
                    for term = terms(:).', normSquared = normSquared+term.coefficient*term.values.^2; end
                    factors = sqrt(abs(normSquared));
                    if name == "generalizedPotentialEnstrophy", factors = factors/sqrt(D); end
                case "surfaceF"
                    factors = abs(self.rawVariable("F",self.zDomain(2)));
                case "waveEnergy"
                    F = self.rawVariable("F",z); G = self.rawVariable("G",z);
                    factors = sqrt(sum(w.*(abs(F).^2+(self.evp.parameters.k*D)^2*abs(G).^2),1)/D);
                otherwise
                    if startsWith(name,"rms"), variable = extractAfter(name,"rms"); else, variable = extractAfter(name,"maxAbs"); end
                    values = self.rawVariable(variable,z);
                    if startsWith(name,"rms")
                        factors = sqrt(sum(w.*abs(values).^2,1)/D);
                    else
                        factors = IMInternalModesNormalizationTools.maximumAmplitudes(@(zz) self.rawVariable(variable,zz),z,values);
                    end
            end
            if any(factors <= 0 | ~isfinite(factors))
                error("IMAnalyticalInternalModesBasis:NullNormalization","A retained mode has zero or nonfinite %s norm.",name);
            end
        end

        function gram = variableGramMatrix(self, variable, zBounds, useNormalized, useMajorant)
            if nargin < 5
                useMajorant = false;
            end
            if useMajorant
                spec = self.evp.majorantInnerProduct(variable);
            else
                spec = self.evp.innerProduct(variable);
            end
            IMAnalyticalInternalModesBasis.assertInnerProductAvailable(spec);
            [z,integrationWeights] = self.integrationGrid(zBounds);
            context = self.context();
            normalizationFactors = [];
            if useNormalized
                normalizationFactors = self.normalizationFactors(self.normalization);
                values = self.rawVariable(variable, z) ./ normalizationFactors;
            else
                values = self.rawVariable(variable, z);
            end
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = values.'*((integrationWeights.*weight(:)).*values);
            gram = 0.5*(gram+gram.');

            endpointTerms = self.endpointGramTerms(variable=variable, zBounds=zBounds, useNormalized=useNormalized);
            for iTerm = 1:numel(endpointTerms)
                valuesAtEndpoint = endpointTerms(iTerm).values(:);
                coefficient = endpointTerms(iTerm).coefficient;
                if useMajorant
                    coefficient = abs(coefficient);
                end
                gram = gram + coefficient*(valuesAtEndpoint*valuesAtEndpoint.');
            end
        end

        function [z,weights] = integrationGrid(self, zBounds)
            originalCount = self.nColumns;
            if ~isempty(fieldnames(self.basisContractState)), originalCount = self.basisContractState.originalColumnCount; end
            [z,weights] = legpts(max(128,4*originalCount+32),sort(zBounds));
            z = z(:); weights = weights(:);
        end

        function context = context(self)
            context.zDomain = self.zDomain;
            context.coordinateKind = "analyticalSolution";
            context.N2 = @(z) self.N2(z);
            context.f0 = self.evp.f0;
            context.g = self.evp.g;
            context.formulation = self.evp.formulation;
            parameterFields = fieldnames(self.evp.parameters);
            for iField = 1:numel(parameterFields)
                fieldName = parameterFields{iField};
                context.(fieldName) = self.evp.parameters.(fieldName);
            end
        end

        function validateZBounds(~, zBounds)
            if zBounds(1) >= zBounds(2)
                error("IMAnalyticalInternalModesBasis:InvalidInterval", "zBounds must be increasing.");
            end
        end

        function validateCoefficientVector(self, coefficients, argumentName)
            if any(~isfinite(coefficients(:)))
                error("IMAnalyticalInternalModesBasis:InvalidCoefficients", "%s must contain finite values.", char(argumentName));
            end
            if length(coefficients) ~= length(self.h)
                error("IMAnalyticalInternalModesBasis:InvalidCoefficientCount", "%s must contain one value for each retained mode (%d).", char(argumentName), length(self.h));
            end
        end

        function name = normalizationName(~, normalization)
            name = IMInternalModesNormalizationTools.name(normalization,"IMAnalyticalInternalModesBasis");
        end
    end

    methods (Static, Access = private)
        function assertInnerProductAvailable(spec)
            if isfield(spec, "hasInnerProduct") && spec.hasInnerProduct
                return;
            end
            error("IMAnalyticalInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this analytical basis. %s", string(spec.variable), string(spec.reason));
        end
    end
    methods (Hidden)
        function values = rawVariable(self, variable, z)
            if variable == "F" && IMInternalModesNormalizationTools.isWave(self.evp)
                values = diff(self.zDomain)*self.rawUz(z);
            else
                values = self.rawVariableFunction(self, string(variable), z(:)).*self.orientationSigns;
            end
        end

        function values = rawUz(self, z)
            values = self.rawUzFunction(self, z(:)).*self.orientationSigns;
        end

        function factor = innerProductNormFactor(self, iMode, options)
            % Return the raw inner-product norm factor.
            %
            % This developer utility returns the raw scale factor
            % $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$
            % for exact analytical `F` or `G` modes before the active
            % normalization is applied. This diagonal magnitude is not a norm
            % for arbitrary modal combinations in an indefinite pairing. The
            % requested variable must have a known inner product.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: raw inner-product scale factor
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            gram = self.variableGramMatrix(string(options.variable), self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

    end

end
