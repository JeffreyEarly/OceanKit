classdef IMInternalModesBasis < IMBasisSet
    % Store solved internal-mode basis functions.
    %
    % `IMInternalModesBasis` evaluates the physical variables with explicit
    % `F(z)` and `G(z)` methods. If the EVP solves `G`, then `F(z)` is
    % recovered by `evp.FfromGz`, whose default relation is
    % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
    % If the EVP solves `F`, then `G(z)` is recovered by `evp.GfromFz`.
    % The default inverse is hydrostatic, but wave factories install
    % wave-specific inverse relations.
    % Normalization is shared across both variables: if a rule gives scale
    % $$s_j$$, then both diagnostic variables for mode $$j$$ are divided by
    % that same factor. Standard rules are installed on the basis set, and
    % custom rules are added after solving with `addNormalization`.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % F = basisSet.F(z);
    % G = basisSet.G(z);
    % ```
    %
    % - Topic: Create internal-mode bases
    % - Topic: Evaluate modes
    % - Topic: Analyze modes
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModesBasis < IMBasisSet

    properties (SetAccess = protected)
        % Equivalent depths for the retained internal modes.
        %
        % These are computed from the parent EVP as
        % $$h_j=\texttt{evp.hFromEigenvalue}(\lambda_j),$$
        % where $$\lambda_j$$ is `eigenvalues(j)`.
        %
        % - Topic: Inspect basis sets
        h

        % Buoyancy frequency squared profile.
        %
        % `N2` has signature `values = N2(z)` and is copied from the
        % internal-mode EVP.
        %
        % - Topic: Inspect basis sets
        N2
    end

    methods
        function self = IMInternalModesBasis(options)
            % Create an internal-mode basis set.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMInternalModesBasis(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: internal-mode EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.modeSelectionDiagnostics: mode-selection diagnostics
            % - Parameter options.normalization: active normalization rule name or enum value
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: internal-mode basis set
            arguments
                options.solver IMSolver
                options.evp IMInternalModes
                options.nativeModes (:,:) double
                options.eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
            end

            self@IMBasisSet(solver=options.solver, evp=options.evp, nativeModes=options.nativeModes, eigenvalues=options.eigenvalues, modeNumber=options.modeNumber, modeSelectionDiagnostics=options.modeSelectionDiagnostics, normalization=options.normalization, metadata=options.metadata);
            self.N2 = options.evp.N2;
            self.h = reshape(options.evp.hFromEigenvalue(self.eigenvalues),1,[]);
            if length(self.h) ~= size(self.nativeModes,2)
                error("IMInternalModesBasis:InvalidEquivalentDepthCount", "hFromEigenvalue must return one equivalent depth for each native mode column.");
            end

            availableNormalizations = IMInternalModesNormalizationTools.available(self.evp);
            self = self.addVectorNormalization("uMax",@(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode,variable="F"),@(basisSet) basisSet.maxAmplitudeNormFactors(variable="F"));
            self = self.addVectorNormalization("wMax",@(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode,variable="G"),@(basisSet) basisSet.maxAmplitudeNormFactors(variable="G"));
            if ismember("surfacePressure",availableNormalizations)
                self = self.addVectorNormalization("surfacePressure",@(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode),@(basisSet) basisSet.surfacePressureNormFactors());
            end
            if ismember("geostrophic",availableNormalizations)
                self = self.addVectorNormalization("geostrophic",@(basisSet,iMode) basisSet.geostrophicNormFactor(iMode),@(basisSet) basisSet.geostrophicNormFactors());
            end
            if ismember("depth",availableNormalizations)
                self = self.addVectorNormalization("depth",@(basisSet,iMode) basisSet.depthNormFactor(iMode),@(basisSet) basisSet.depthNormFactors());
            end
            if ismember("generalizedPotentialEnstrophy",availableNormalizations)
                self = self.addVectorNormalization("generalizedPotentialEnstrophy", ...
                    @(basisSet,iMode) basisSet.generalizedPotentialEnstrophyNormFactor(iMode), ...
                    @(basisSet) basisSet.generalizedPotentialEnstrophyNormFactors());
            end
            if ismember("kConstant",availableNormalizations)
                self = self.addVectorNormalization("kConstant",@(basisSet,iMode) basisSet.innerProductNormFactor(iMode,variable="G"),@(basisSet) basisSet.innerProductNormFactors(variable="G"));
            end
            if isempty(options.normalization)
                self.normalization = IMInternalModesNormalizationTools.default(self.evp);
            end
        end

        function G = G(self, z, options)
            % Evaluate `G` modes.
            %
            % If the EVP formulation is `G`, this evaluates the solved
            % canonical variable. If the formulation is `F`, `G` is recovered
            % by `evp.GfromFz`. The default relation is
            % $$G_j(z)=-\frac{g}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z),$$
            % but individual EVPs may supply a different diagnostic
            % relation.
            %
            % - Topic: Evaluate modes
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name or enum value
            % - Returns G: evaluated `G` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            G = self.rawVariable("G", z) ./ self.normalizationFactors(options.normalization);
        end

        function F = F(self, z, options)
            % Evaluate `F` modes.
            %
            % If the EVP formulation is `F`, this evaluates the solved
            % canonical variable. If the formulation is `G`, `F` is recovered
            % by `evp.FfromGz`; the default relation is
            % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
            %
            % - Topic: Evaluate modes
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name or enum value
            % - Returns F: evaluated `F` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            F = self.rawVariable("F", z) ./ self.normalizationFactors(options.normalization);
        end

        function spec = majorantInnerProduct(self, options)
            % Return the induced positive Hilbert-majorant recipe.
            %
            % The recipe retains the positive interior weight and replaces
            % every signed endpoint coefficient by its absolute value. Use
            % `evp.innerProduct` for the signed Pontryagin recipe used by
            % projection and signed invariants.
            %
            % - Topic: Analyze modes
            % - Declaration: spec = majorantInnerProduct(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spec: positive interior and absolute-endpoint recipe
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end
            spec = self.evp.majorantInnerProduct(options.variable);
        end

        function gram = gramMatrix(self, options)
            % Return the signed Gram matrix for `F` or `G`.
            %
            % With no arguments this uses the solved formulation over the full
            % basis-set domain. Use `variable="F"` or `variable="G"` to choose
            % a physical variable, and `zBounds=[zMin zMax]` to restrict the
            % interior integral. Endpoint terms are included only when the
            % interval contains the corresponding physical endpoint:
            %
            % $$
            % M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
            % \sum_\ell \gamma_\ell L_\ell[V_i]L_\ell[V_j]+
            % \sum_\ell \alpha_\ell V_i(z_\ell)V_j(z_\ell).
            % $$
            %
            % This is the Pontryagin pairing used for modal projection and
            % signed invariants. It can be indefinite. Use
            % `majorantGramMatrix` when a positive matrix is required for
            % magnitudes, error tolerances, or convergence tests.
            %
            % Use `endpointGramTerms` to inspect the prepared endpoint
            % vectors that generate the rank-one endpoint updates.
            % The requested variable must have a known inner product; if it
            % does not, this method throws
            % `IMInternalModesBasis:UnavailableInnerProduct` rather than
            % returning an incomplete Gram matrix.
            %
            % - Topic: Analyze modes
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns gram: Gram matrix
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
        end

        function gram = majorantGramMatrix(self, options)
            % Return the positive Hilbert-majorant Gram matrix.
            %
            % The majorant retains the positive interior contribution and
            % replaces every signed endpoint coefficient by its absolute
            % value. It is the positive product associated with the natural
            % $$L^2\oplus\mathbb C^s$$ decomposition of the Pontryagin
            % space. It coincides with `gramMatrix` when all endpoint
            % coefficients are nonnegative.
            %
            % - Topic: Analyze modes
            % - Declaration: gram = majorantGramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns gram: positive Hilbert-majorant Gram matrix
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true, true);
        end

        function value = majorantNorm(self, coefficients, options)
            % Return the positive Hilbert-majorant norm of modal coefficients.
            %
            % For coefficient vector $$c$$ this returns
            % $$\sqrt{c^*M_+c}$$, where $$M_+$$ is
            % `majorantGramMatrix`. Do not replace this with
            % $$\sqrt{|c^*Mc|}$$ for a signed Gram matrix $$M$$; the latter
            % can vanish for a nonzero state and is not a norm.
            %
            % - Topic: Analyze modes
            % - Declaration: value = majorantNorm(basisSet,coefficients,options)
            % - Parameter coefficients: one coefficient per retained mode
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns value: positive scalar norm
            arguments
                self IMInternalModesBasis
                coefficients (:,1) double {mustBeFinite}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            nModes = size(self.nativeModes,2);
            if length(coefficients) ~= nModes
                error("IMBasisSet:InvalidCoefficientCount", "Coefficient vectors must contain one value for each retained mode (%d).", nModes);
            end
            gram = self.majorantGramMatrix(variable=options.variable,zBounds=options.zBounds);
            valueSquared = real(coefficients(:)'*(gram*coefficients(:)));
            tolerance = 1e3*eps(max(1,norm(gram,2)*norm(coefficients,2)^2));
            if valueSquared < -tolerance
                error("IMInternalModesBasis:InvalidMajorantGramMatrix", "The computed majorant quadratic form is negative beyond roundoff.");
            end
            value = sqrt(max(0,valueSquared));
        end

        function terms = endpointGramTerms(self, options)
            % Prepare rank-one endpoint terms for `F` or `G` Gram matrices.
            %
            % `endpointGramTerms` returns endpoint vectors over all retained
            % modes. Solved-form endpoint weights use
            %
            % $$
            % L_\ell[V_j]=c_\ell V_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial V_j}{\partial z}(z_\ell),
            % $$
            %
            % and contribute
            %
            % $$
            % M \leftarrow M+\gamma_\ell L_\ell L_\ell^\mathsf{T}.
            % $$
            %
            % Catalog endpoint value terms use $$V_j(z_\ell)$$ and
            % contribute
            %
            % $$
            % M \leftarrow M+\alpha_\ell V_\ell V_\ell^\mathsf{T}.
            % $$
            %
            % Endpoint terms are omitted when `zBounds` does not include
            % that endpoint.
            %
            % - Topic: Developer topics — Gram-matrix assembly
            % - Declaration: terms = endpointGramTerms(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Parameter options.normalization: normalization rule name or enum value
            % - Parameter options.useNormalized: whether returned values use the normalization
            % - Returns terms: struct array with `location`, `coefficient`, `values`, and `kind`
            % - Developer: true
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
                options.normalization = self.normalization
                options.useNormalized (1,1) logical = true
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end

            variable = string(options.variable);
            spec = self.evp.innerProduct(variable);
            if ~isfield(spec, "hasInnerProduct") || ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot be used as a Gram matrix. %s", string(spec.variable), string(spec.reason));
            end
            terms = struct("location", {}, "coefficient", {}, "values", {}, "kind", {});
            if variable == self.evp.formulation
                terms = endpointGramTerms@IMBasisSet(self, zBounds=options.zBounds, normalization=options.normalization, useNormalized=options.useNormalized);
            end

            normalizationFactors = [];
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
                terms(end+1) = struct("location", string(endpointTerm.location), "coefficient", endpointTerm.coefficient, "values", values, "kind", "endpointInnerProductTerm"); %#ok<AGROW>
            end
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial-depth Gram matrix for `F` or `G`.
            %
            % This computes the eigendecomposition of the symmetric Gram
            % matrix on `zBounds`. If `variable` is omitted, the solved
            % formulation is used. The requested variable must have a known
            % inner product.
            %
            % - Topic: Analyze modes
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute an internal-mode modal spectrum.
            %
            % If `options.variable` is omitted, the solved formulation is
            % used. The requested variable must have a known inner product.
            % This is a signed spectrum: entries associated with negative
            % Pontryagin directions can be negative. Use `majorantNorm` for
            % a positive total magnitude; a generally additive per-mode
            % majorant spectrum does not exist because the majorant Gram
            % matrix need not be diagonal in this basis.
            %
            % - Topic: Analyze modes
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: optional variable name, `"F"` or `"G"`
            % - Returns spectrum: modal spectrum
            arguments
                self IMInternalModesBasis
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute an internal-mode modal cross-spectrum.
            %
            % If `options.variable` is omitted, the solved formulation is
            % used. The requested variable must have a known inner product.
            %
            % - Topic: Analyze modes
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: optional variable name, `"F"` or `"G"`
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMInternalModesBasis
                coefficientsA (:,1) double {mustBeFinite}
                coefficientsB (:,1) double {mustBeFinite}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            nModes = size(self.nativeModes,2);
            if length(coefficientsA) ~= nModes || length(coefficientsB) ~= nModes
                error("IMBasisSet:InvalidCoefficientCount", "Coefficient vectors must contain one value for each retained mode (%d).", nModes);
            end
            gram = self.gramMatrix(variable=options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end

        [transform, assessment] = discreteTransform(self, options)
        [weights, weightFit] = quadratureWeightsForPoints(self, options)
    end

    methods
        function [self, signs] = orientModeSigns(self)
            % Orient modes so `G` is positive immediately below the surface.
            %
            % A resolved nonzero surface `G` value sets the sign directly.
            % When `G` vanishes at the surface, the sign is chosen from
            % $$-G_z(z_s)$$, the leading one-sided Taylor coefficient into
            % the ocean. The rigid-lid barotropic mode has `G` identically
            % zero, so that known `F`-form zero mode uses `F` as a fallback.
            % The same sign flip is applied to the coupled `F`/`G` pair.
            % `IMInternalModesBasis` is a value class, so callers must keep
            % the returned basis set:
            %
            % ```matlab
            % basisSet = basisSet.orientModeSigns();
            % ```
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with oriented native mode signs
            % - Developer: true
            signs = ones(1,size(self.nativeModes,2));
            if isempty(self.nativeModes)
                return;
            end

            zNative = self.solver.zNative;
            surfaceIndex = self.solver.boundaryIndex("surface");
            GValues = self.rawVariable("G",zNative);
            GzValues = self.solver.differentiateGridValues(GValues,1);
            FValues = self.rawVariable("F",zNative);
            allowFFallback = self.evp.formulation == "F" & self.eigenvalues == 0;
            signs = IMModeOrientationTools.shallowInteriorGPositive( ...
                GValues=GValues,GzSurface=GzValues(surfaceIndex,:),FValues=FValues, ...
                depth=diff(self.zDomain),surfaceIndex=surfaceIndex,allowFFallback=allowFFallback);
            self.nativeModes = self.nativeModes .* signs;
            self.metadata.modeOrientation = IMModeOrientationTools.convention;
        end

        function factor = innerProductNormFactor(self, iMode, options)
            % Return the `F` or `G` inner-product norm factor.
            %
            % This is the raw factor
            % $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$ for `variable` equal
            % to `F` or `G`. If `variable` is omitted, the solved
            % formulation is used. This is a per-mode normalization
            % convention, not a norm for arbitrary modal combinations. The
            % requested variable must have a known inner product. Custom normalization rules
            % registered with `addNormalization` call this method.
            % The factor is computed from raw, unnormalized modes before
            % the active basis normalization is applied.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: raw inner-product norm factor
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            variable = string(options.variable);
            z = self.solver.innerProductGrid(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct(variable);
            if ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot define a normalization factor. %s",variable,string(spec.reason));
            end
            values = self.rawVariable(variable,z);
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            else
                weight = weight(:);
            end
            value = self.solver.integrateInnerProduct(z,weight.*values(:,iMode).*values(:,iMode),self.zDomain);
            endpointTerms = self.endpointGramTerms(variable=variable,zBounds=self.zDomain,useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                endpointValues = endpointTerms(iTerm).values;
                value = value+endpointTerms(iTerm).coefficient*endpointValues(iMode)^2;
            end
            factor = sqrt(abs(value));
        end

        function factor = geostrophicNormFactor(self, iMode)
            % Return the hydrostatic geostrophic normalization factor.
            %
            % This developer utility implements the normalization rule
            % installed for `modeFamily="hydrostatic"`. It chooses one
            % shared raw scale factor for each coupled `F`/`G` mode.
            % Baroclinic modes use
            % $$s_j^2=\langle G_j,G_j\rangle_G,$$
            % so normalized modes satisfy
            % $$\langle G_j,G_j\rangle_G=1,\qquad
            % \langle F_j,F_j\rangle_F=h_j.$$
            % A barotropic zero mode uses the `F` norm divided by
            % $$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$ and is a
            % separate null-mode convention.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw geostrophic scale factor
            % - Developer: true
            if self.modeNumber(iMode) == 0
                factor = self.innerProductNormFactor(iMode,variable="F")/sqrt(diff(self.zDomain));
                return;
            end
            factor = self.innerProductNormFactor(iMode, variable="G");
        end

        function factor = depthNormFactor(self, iMode)
            % Return the volume-only depth normalization factor.
            %
            % This developer utility evaluates the positive factor
            % $$s_j=\sqrt{D^{-1}\int_{z_b}^{z_s}(F_j^{\mathrm{raw}})^2\,dz}$$
            % without endpoint terms. The same factor scales `F` and `G`.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = depthNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: positive volume root-mean-square `F` factor
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            z = self.solver.innerProductGrid(self.zDomain);
            F = self.rawVariable("F", z);
            integral = self.solver.integrateInnerProduct(z, F(:,iMode).*F(:,iMode), self.zDomain);
            factor = sqrt(max(0, real(integral)/diff(self.zDomain)));
        end

        function factor = generalizedPotentialEnstrophyNormFactor(self, iMode)
            % Return the depth-mean generalized-potential-enstrophy factor.
            %
            % The solved `F` inner product contains the interior and active
            % endpoint terms of the generalized potential enstrophy. This
            % factor scales that complete inner product by the physical
            % depth, so normalized modes satisfy
            % $$D^{-1}\langle F_i,F_j\rangle_\alpha=\delta_{ij}.$$
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = generalizedPotentialEnstrophyNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: generalized-potential-enstrophy normalization factor
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            factor = self.innerProductNormFactor(iMode,variable="F")/sqrt(diff(self.zDomain));
        end

        function factor = maxAmplitudeNormFactor(self, iMode, options)
            % Return the maximum amplitude of `F` or `G`.
            %
            % This is $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$ for the
            % requested variable on the basis-set integration grid. If
            % `variable` is omitted, the solved formulation is used.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = maxAmplitudeNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: maximum absolute variable amplitude
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            variable = string(options.variable);
            z = self.solver.innerProductGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            factor = max(abs(values(:,iMode)));
        end

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the raw surface `F` value.
            %
            % This developer utility returns the raw surface value
            % $$s_j=F_j^{\mathrm{raw}}(z_\mathrm{surface}).$$
            % It gives unit surface `F` value when the raw surface value is
            % finite and nonzero. If that value is unavailable, zero, or
            % nonfinite, it returns `1` so normalization remains well
            % defined.
            %
            % - Topic: Developer topics — Normalization rules
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

        function values = rawVariable(self, variable, z)
            % Evaluate raw physical `F` or `G` modes.
            %
            % This developer utility evaluates unnormalized internal-mode
            % physical variables. The public `F` and `G` methods apply the
            % active basis normalization after this step:
            % $$F_j(z)=F_j^{\mathrm{raw}}(z)/s_j,\qquad
            % G_j(z)=G_j^{\mathrm{raw}}(z)/s_j.$$
            % If `variable` is the solved formulation, values come
            % directly from the native modes. Otherwise the diagnostic
            % variable is recovered through the EVP-owned diagnostic
            % relation.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: values = rawVariable(basisSet,variable,z)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized physical mode values
            % - Developer: true
            arguments
                self IMInternalModesBasis
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])}
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            variable = string(variable);
            context = self.evp.contextForSolver(self.solver);
            if variable == self.evp.formulation
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end

            switch self.evp.formulation
                case "G"
                    if variable ~= "F"
                        self.unsupported("evaluate " + variable);
                    end
                    dGdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.FfromGz(z(:), dGdz, self.h, context);
                case "F"
                    if variable ~= "G"
                        self.unsupported("evaluate " + variable);
                    end
                    dFdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.GfromFz(z(:), dFdz, self.h, context);
                otherwise
                    self.unsupported("evaluate " + variable);
            end
        end
    end

    methods (Hidden)
        function factors = innerProductNormFactors(self, options)
            % Return raw variable inner-product factors for the whole family.
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.evp.formulation
            end
            variable = string(options.variable);
            z = self.solver.innerProductGrid(self.zDomain);
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct(variable);
            if ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot define normalization factors. %s",variable,string(spec.reason));
            end
            values = self.rawVariable(variable,z);
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            else
                weight = weight(:);
            end
            integrationWeights = self.solver.innerProductWeights(z,self.zDomain);
            normSquared = sum((integrationWeights.*weight).*values.*values,1);
            endpointTerms = self.endpointGramTerms(variable=variable,zBounds=self.zDomain,useNormalized=false);
            for iTerm = 1:numel(endpointTerms)
                endpointValues = endpointTerms(iTerm).values;
                normSquared = normSquared+endpointTerms(iTerm).coefficient*endpointValues.^2;
            end
            factors = sqrt(abs(normSquared));
        end

        function factors = geostrophicNormFactors(self)
            % Return hydrostatic geostrophic factors for the whole family.
            factors = self.innerProductNormFactors(variable="G");
            zeroModes = self.modeNumber == 0;
            if any(zeroModes)
                FFactors = self.innerProductNormFactors(variable="F");
                factors(zeroModes) = FFactors(zeroModes)/sqrt(diff(self.zDomain));
            end
        end

        function factors = depthNormFactors(self)
            % Return volume-only F depth factors for the whole family.
            z = self.solver.innerProductGrid(self.zDomain);
            F = self.rawVariable("F",z);
            integrationWeights = self.solver.innerProductWeights(z,self.zDomain);
            normSquared = sum(integrationWeights.*F.*F,1)/diff(self.zDomain);
            factors = sqrt(max(0,real(normSquared)));
        end

        function factors = generalizedPotentialEnstrophyNormFactors(self)
            % Return generalized-potential-enstrophy factors for the family.
            factors = self.innerProductNormFactors(variable="F")/sqrt(diff(self.zDomain));
        end

        function factors = maxAmplitudeNormFactors(self, options)
            % Return raw maximum-amplitude factors for the whole family.
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.evp.formulation
            end
            z = self.solver.innerProductGrid(self.zDomain);
            values = self.rawVariable(string(options.variable),z);
            factors = max(abs(values),[],1);
        end

        function factors = surfacePressureNormFactors(self)
            % Return raw surface-F factors for the whole family.
            factors = self.rawVariable("F",self.zDomain(2));
            invalid = ~isfinite(factors) | abs(factors) <= 1e-12;
            factors(invalid) = 1;
        end
    end

    methods (Access = protected)
        function transform = buildInternalModesDiscreteTransform(self, z, weights, nModes, variables)
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                weights (:,1) double {mustBeReal, mustBeFinite}
                nModes (1,1) double {mustBeInteger, mustBePositive}
                variables (1,:) string
            end

            preparation = self.prepareInternalModesDiscreteTransform(z,weights,nModes,variables);
            transform = IMInternalModesDiscreteTransform(z=preparation.z,weights=preparation.weights,modeNumber=preparation.modeNumber,h=preparation.h, ...
                normalization=preparation.normalization,inverseF=preparation.inverseF,inverseG=preparation.inverseG,endpointF=preparation.endpointF,endpointG=preparation.endpointG, ...
                channelData=preparation.channelData,primaryVariable=preparation.primaryVariable,zDomain=preparation.zDomain,g=preparation.g, ...
                modeFamily=preparation.modeFamily,N2Values=preparation.N2Values,problemMetadata=preparation.problemMetadata);
        end

        function preparation = prepareInternalModesDiscreteTransform(self, z, weights, nModes, variables)
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                weights (:,1) double {mustBeReal, mustBeFinite}
                nModes (1,1) double {mustBeInteger, mustBePositive}
                variables (1,:) string
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
            if any(z < self.zDomain(1)-domainTolerance) || any(z > self.zDomain(2)+domainTolerance)
                error("IMBasisSet:InvalidDiscreteGrid", "All sample points must lie inside the basis-set zDomain.");
            end
            if ~any(weights ~= 0)
                error("IMBasisSet:InvalidDiscreteWeights", "weights must contain at least one nonzero value.");
            end

            variables = self.canonicalDiscreteVariables(variables);
            metricPreparation = struct();
            for variable = ["F","G"]
                [available,reason,metricMatrix,interiorWeight,endpointMetricMatrix] = self.sampledInternalModesMetric(variable,z,weights);
                metricPreparation.(char(variable)) = struct(available=available,reason=reason,metricMatrix=metricMatrix, ...
                    interiorWeight=interiorWeight,endpointMetricMatrix=endpointMetricMatrix);
                if ismember(variable,variables) && ~available
                    error("IMInternalModesBasis:UnavailableDiscreteTransformVariable", "The %s channel cannot provide a direct forward transform. %s", variable, reason);
                end
            end

            normalizationFactors = self.normalizationFactors(self.normalization);
            normalizationFactors = normalizationFactors(1:nModes);
            inverseF = self.rawVariable("F",z);
            inverseG = self.rawVariable("G",z);
            inverseF = inverseF(:,1:nModes)./normalizationFactors;
            inverseG = inverseG(:,1:nModes)./normalizationFactors;
            endpointZ = [self.zDomain(2);self.zDomain(1)];
            endpointF = self.rawVariable("F",endpointZ);
            endpointG = self.rawVariable("G",endpointZ);
            endpointF = endpointF(:,1:nModes)./normalizationFactors;
            endpointG = endpointG(:,1:nModes)./normalizationFactors;

            channelData = struct();
            for variable = ["F","G"]
                preparation = metricPreparation.(char(variable));
                requested = ismember(variable,variables);
                if requested
                    rawGram = self.variableGramMatrix(variable,self.zDomain,false);
                    rawMajorantGram = self.variableGramMatrix(variable,self.zDomain,false,true);
                    targetGramMatrix = rawGram(1:nModes,1:nModes)./(normalizationFactors(:)*normalizationFactors(:).');
                    targetMajorantGramMatrix = rawMajorantGram(1:nModes,1:nModes)./(normalizationFactors(:)*normalizationFactors(:).');
                    targetGramMatrix = 0.5*(targetGramMatrix+targetGramMatrix.');
                    targetMajorantGramMatrix = 0.5*(targetMajorantGramMatrix+targetMajorantGramMatrix.');
                    if variable == "F"
                        sampled = inverseF;
                    else
                        sampled = inverseG;
                    end
                    activeMask = self.internalModesTransformActiveMask(variable,sampled,targetGramMatrix,targetMajorantGramMatrix);
                    channelData.(char(variable)) = struct(available=true,reason="",activeModeMask=activeMask,metricMatrix=preparation.metricMatrix, ...
                        targetGramMatrix=targetGramMatrix,targetMajorantGramMatrix=targetMajorantGramMatrix, ...
                        interiorWeight=preparation.interiorWeight,endpointMetricMatrix=preparation.endpointMetricMatrix);
                else
                    reason = preparation.reason;
                    if preparation.available
                        reason = "The " + variable + " channel was not requested when this transform was built.";
                    end
                    metricMatrix = zeros(0,0);
                    targetGramMatrix = zeros(0,0);
                    targetMajorantGramMatrix = zeros(0,0);
                    activeMask = false(1,nModes);
                    channelData.(char(variable)) = struct(available=false,reason=string(reason),activeModeMask=activeMask,metricMatrix=metricMatrix, ...
                        targetGramMatrix=targetGramMatrix,targetMajorantGramMatrix=targetMajorantGramMatrix);
                end
            end

            metadata = self.metadata;
            primaryVariable = string(self.evp.formulation);
            if ~ismember(primaryVariable,variables)
                primaryVariable = variables(1);
            end
            preparation = struct(z=z,weights=weights,modeNumber=self.modeNumber(1:nModes),h=self.h(1:nModes), ...
                normalization=self.normalizationName(self.normalization),inverseF=inverseF,inverseG=inverseG,endpointF=endpointF,endpointG=endpointG, ...
                channelData=channelData,primaryVariable=primaryVariable,zDomain=self.zDomain,g=self.evp.g,modeFamily=self.evp.modeFamily,N2Values=self.N2(z),problemMetadata=metadata);
        end

        function variables = directlyRepresentableDiscreteVariables(self, z)
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end
            variables = strings(1,0);
            placeholderWeights = ones(size(z));
            for variable = ["F","G"]
                [available,~] = self.sampledInternalModesMetric(variable,z,placeholderWeights);
                if available
                    variables(end+1) = variable; %#ok<AGROW>
                end
            end
        end

        function variables = canonicalDiscreteVariables(~, variables)
            variables = string(variables(:).');
            if any(~ismember(variables,["F","G"]))
                invalid = variables(find(~ismember(variables,["F","G"]),1));
                error("IMInternalModesBasis:InvalidDiscreteTransformVariable", "Unknown internal-mode transform variable ""%s"". Choose F or G.", invalid);
            end
            canonical = ["F","G"];
            variables = canonical(ismember(canonical,variables));
        end

        function [available,reason,metricMatrix,targetGramMatrix,targetMajorantGramMatrix,activeMask,interiorWeight,endpointMetricMatrix] = prepareInternalModesTransformChannel(self, variable, z, weights, nModes)
            [available,reason,metricMatrix,interiorWeight,endpointMetricMatrix] = self.sampledInternalModesMetric(variable,z,weights);
            targetGramMatrix = zeros(0,0);
            targetMajorantGramMatrix = zeros(0,0);
            activeMask = false(1,nModes);
            if ~available
                return;
            end
            normalizationFactors = self.normalizationFactors(self.normalization);
            normalizationFactors = normalizationFactors(1:nModes);
            rawGram = self.variableGramMatrix(string(variable),self.zDomain,false);
            rawMajorantGram = self.variableGramMatrix(string(variable),self.zDomain,false,true);
            targetGramMatrix = rawGram(1:nModes,1:nModes)./(normalizationFactors(:)*normalizationFactors(:).');
            targetMajorantGramMatrix = rawMajorantGram(1:nModes,1:nModes)./(normalizationFactors(:)*normalizationFactors(:).');
            targetGramMatrix = 0.5*(targetGramMatrix+targetGramMatrix.');
            targetMajorantGramMatrix = 0.5*(targetMajorantGramMatrix+targetMajorantGramMatrix.');
            sampled = self.rawVariable(string(variable),z);
            sampled = sampled(:,1:nModes)./normalizationFactors;
            activeMask = self.internalModesTransformActiveMask(string(variable),sampled,targetGramMatrix,targetMajorantGramMatrix);
        end

        function [available,reason,metricMatrix,interiorWeight,endpointMetricMatrix] = sampledInternalModesMetric(self, variable, z, weights, useMajorant)
            if nargin < 5
                useMajorant = false;
            end
            variable = string(variable);
            z = z(:);
            weights = weights(:);
            metricMatrix = zeros(0,0);
            interiorWeight = zeros(size(z));
            endpointMetricMatrix = zeros(length(z));
            if useMajorant
                spec = self.evp.majorantInnerProduct(variable);
            else
                spec = self.evp.innerProduct(variable);
            end
            if ~spec.hasInnerProduct
                available = false;
                reason = "Its continuous inner product is unavailable. " + string(spec.reason);
                return;
            end

            context = self.evp.contextForSolver(self.solver);
            interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
            if isscalar(interiorWeight)
                interiorWeight = interiorWeight*ones(size(z));
            else
                interiorWeight = interiorWeight(:);
            end
            if numel(interiorWeight) ~= length(z) || ~isreal(interiorWeight) || any(~isfinite(interiorWeight))
                error("IMBasisSet:InvalidDiscreteMetricWeight", "The %s interior weight must return one finite real value for each sample point.", variable);
            end
            metricMatrix = diag(interiorWeight.*weights);
            domainTolerance = 100*eps(max(1,max(abs(self.zDomain))));

            endpointWeights = [spec.surfaceWeights;spec.bottomWeights];
            for iWeight = 1:numel(endpointWeights)
                endpointWeight = endpointWeights(iWeight);
                if endpointWeight.d ~= 0
                    available = false;
                    reason = "Its " + string(endpointWeight.location) + " endpoint functional contains a derivative trace; direct transforms require value-only sampled endpoint terms.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                if endpointWeight.c == 0
                    continue;
                end
                [endpointIndex,zEndpoint] = endpointSampleIndex(z,self.zDomain,string(endpointWeight.location),domainTolerance);
                if isempty(endpointIndex)
                    available = false;
                    reason = "The " + string(endpointWeight.location) + " endpoint at z=" + string(zEndpoint) + " must be included in z to represent its value-only metric term.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                if ~isfinite(endpointWeight.coefficient)
                    available = false;
                    reason = "Its " + string(endpointWeight.location) + " endpoint metric coefficient is not finite.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                endpointMetricMatrix(endpointIndex,endpointIndex) = endpointMetricMatrix(endpointIndex,endpointIndex) + endpointWeight.coefficient*endpointWeight.c*endpointWeight.c;
            end
            for iTerm = 1:numel(spec.endpointInnerProductTerms)
                term = spec.endpointInnerProductTerms(iTerm);
                if string(term.variable) ~= variable
                    available = false;
                    reason = "Its " + string(term.location) + " endpoint inner product uses companion-variable " + string(term.variable) + " values; paired admissible-state projection is not yet supported.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                [endpointIndex,zEndpoint] = endpointSampleIndex(z,self.zDomain,string(term.location),domainTolerance);
                if isempty(endpointIndex)
                    available = false;
                    reason = "The " + string(term.location) + " endpoint at z=" + string(zEndpoint) + " must be included in z to represent its value-only metric term.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                if ~isfinite(term.coefficient)
                    available = false;
                    reason = "Its " + string(term.location) + " endpoint metric coefficient is not finite.";
                    metricMatrix = zeros(0,0);
                    return;
                end
                endpointMetricMatrix(endpointIndex,endpointIndex) = endpointMetricMatrix(endpointIndex,endpointIndex) + term.coefficient;
            end
            metricMatrix = metricMatrix + endpointMetricMatrix;
            available = true;
            reason = "";
        end

        function activeMask = internalModesTransformActiveMask(self, variable, sampled, targetGramMatrix, targetMajorantGramMatrix)
            targetNorms = diag(targetGramMatrix).';
            majorantNorms = diag(targetMajorantGramMatrix).';
            targetTolerance = 1e3*eps(max(1,majorantNorms));
            columnNorms = vecnorm(sampled,2,1);
            columnTolerance = 1e3*eps(max(1,columnNorms));
            zeroTarget = abs(targetNorms) <= targetTolerance;
            zeroSample = columnNorms <= columnTolerance;
            unsupportedZeroNorm = zeroTarget & ~zeroSample;
            if any(unsupportedZeroNorm)
                labels = self.modeNumber(1:length(unsupportedZeroNorm));
                labels = labels(unsupportedZeroNorm);
                error("IMInternalModesBasis:UnsupportedZeroNormMode", "The %s channel has nonzero sampled columns with zero continuous norm for mode label(s) %s.", variable, join(string(labels),", "));
            end
            activeMask = ~(zeroTarget & zeroSample);
        end

        function gram = variableGramMatrix(self, variable, zBounds, useNormalized, useMajorant)
            if nargin < 5
                useMajorant = false;
            end
            z = self.solver.innerProductGrid(zBounds);
            context = self.evp.contextForSolver(self.solver);
            if useMajorant
                spec = self.evp.majorantInnerProduct(variable);
            else
                spec = self.evp.innerProduct(variable);
            end
            if ~isfield(spec, "hasInnerProduct") || ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot be used as a Gram matrix. %s", string(spec.variable), string(spec.reason));
            end
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
            integrationWeights = self.solver.innerProductWeights(z,zBounds);
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
    end

end

function [index,zEndpoint] = endpointSampleIndex(z,zDomain,location,tolerance)
if location == "surface"
    zEndpoint = zDomain(2);
else
    zEndpoint = zDomain(1);
end
index = find(abs(z-zEndpoint) <= tolerance,1);
end
