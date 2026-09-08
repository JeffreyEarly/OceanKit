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
        normalization = "unity"
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
    end

    methods
        function self = IMAnalyticalInternalModesBasis(options)
            % Create an exact internal-mode basis.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: basisSet = IMAnalyticalInternalModesBasis(options)
            % - Parameter options.solution: analytical solution family
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
                options.solution IMAnalyticalSolution
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

            self.solution = options.solution;
            self.evp = options.evp;
            self.h = reshape(options.h,1,[]);
            self.eigenvalues = 1 ./ self.h;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.zDomain = options.solution.zDomain;
            self.N2 = options.N2;
            self.metadata = options.metadata;
            self.modeSelectionDiagnostics = options.modeSelectionDiagnostics;
            self.rawVariableFunction = options.rawVariableFunction;
            self.rawUzFunction = options.rawUzFunction;
            self.orientationSigns = ones(1,length(self.h));
            zOrientation = self.integrationGrid(self.zDomain);
            GValues = self.rawVariableFunction(self,"G",zOrientation);
            FValues = self.rawVariableFunction(self,"F",zOrientation);
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
            self = self.installNormalizationRules();
            if isempty(options.normalization)
                self.normalization = IMInternalModesNormalizationTools.default(self.evp);
            else
                self.normalization = options.normalization;
            end
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

            name = self.normalizationName(name);
            self.normalizationNameMap{name} = rule;
        end

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

            F = self.rawVariable("F", z) ./ self.normalizationFactors(options.normalization);
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

            G = self.rawVariable("G", z) ./ self.normalizationFactors(options.normalization);
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
            factors = zeros(1,length(self.h));
            for iMode = 1:length(self.h)
                factors(iMode) = rule(self, iMode);
            end
            factors = abs(factors);
            factors(factors == 0 | ~isfinite(factors)) = 1;
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
            % - Parameter options.normalization: normalization rule name or enum value
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
        function factor = innerProductNormFactor(self, iMode, options)
            % Return the raw inner-product norm factor.
            %
            % This developer utility returns the raw scale factor
            % $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$
            % for exact analytical `F` or `G` modes before the active
            % normalization is applied. This is a per-mode normalization
            % convention, not a norm for arbitrary modal combinations. The
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

        function factor = geostrophicNormFactor(self, iMode)
            % Return the geostrophic normalization factor.
            %
            % This developer utility implements the analytical
            % `modeFamily="hydrostatic"` normalization rule. Baroclinic
            % modes use one shared scale factor for the coupled `F`/`G`
            % pair based on the raw `G` inner product. A barotropic zero
            % mode uses the `F` norm divided by
            % $$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$.
            %
            % - Topic: Developer topics
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw geostrophic scale factor
            % - Developer: true
            if self.modeNumber(iMode) == 0
                gram = self.variableGramMatrix("F", self.zDomain, false);
                factor = sqrt(abs(gram(iMode,iMode)))/sqrt(diff(self.zDomain));
                return;
            end
            factor = self.innerProductNormFactor(iMode, variable="G");
        end

        function factor = depthNormFactor(self, iMode)
            % Return the volume-only depth normalization factor.
            %
            % This developer utility evaluates the positive factor
            % $$s_j=\sqrt{D^{-1}\int_{z_b}^{z_s}(F_j^{\mathrm{raw}})^2\,dz}$$
            % without endpoint terms. The same factor scales exact `F` and
            % `G` modes.
            %
            % - Topic: Developer topics
            % - Declaration: factor = depthNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: positive volume root-mean-square `F` factor
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            z = self.integrationGrid(self.zDomain);
            F = self.rawVariable("F", z);
            factor = sqrt(max(0, real(trapz(z, F(:,iMode).*F(:,iMode)))/diff(self.zDomain)));
        end

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

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the raw surface `F` value.
            %
            % This developer utility returns the exact raw surface value
            % $$s_j=F_j^{\mathrm{raw}}(z_\mathrm{surface}).$$
            % It gives unit surface `F` value when that value is finite and
            % nonzero. If the raw value is unavailable, zero, or nonfinite,
            % it returns `1` so normalization remains well defined.
            %
            % - Topic: Developer topics
            % - Declaration: factor = surfacePressureNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw surface-pressure scale factor
            % - Developer: true
            values = self.rawVariable("F", self.zDomain(2));
            factor = values(1,iMode);
            if ~isfinite(factor) || abs(factor) <= 1e-12
                factor = 1;
            end
        end
    end

    methods (Access = private)
        function self = installNormalizationRules(self)
            availableNormalizations = IMInternalModesNormalizationTools.available(self.evp);
            self = self.addNormalization("unity", @(basisSet,iMode) basisSet.innerProductNormFactor(iMode));
            self = self.addNormalization("uMax", @(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode, variable="F"));
            self = self.addNormalization("wMax", @(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode, variable="G"));
            if ismember("surfacePressure",availableNormalizations)
                self = self.addNormalization("surfacePressure", @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode));
            end
            if ismember("geostrophic",availableNormalizations)
                self = self.addNormalization("geostrophic", @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode));
            end
            if ismember("depth",availableNormalizations)
                self = self.addNormalization("depth", @(basisSet,iMode) basisSet.depthNormFactor(iMode));
            end
            if ismember("kConstant",availableNormalizations)
                self = self.addNormalization("kConstant", @(basisSet,iMode) basisSet.innerProductNormFactor(iMode, variable="G"));
            end
        end

        function values = rawVariable(self, variable, z)
            values = self.rawVariableFunction(self, string(variable), z(:)).*self.orientationSigns;
        end

        function values = rawUz(self, z)
            values = self.rawUzFunction(self, z(:)).*self.orientationSigns;
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
            z = self.integrationGrid(zBounds);
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
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    value = trapz(z, weight(:).*values(:,iMode).*values(:,jMode));
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end

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

        function z = integrationGrid(~, zBounds)
            z = linspace(min(zBounds), max(zBounds), 1024).';
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
            parts = split(string(normalization), ".");
            name = parts(end);
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
end
