classdef IMInternalModes < IMEigenvalueProblem
    % Describe canonical EVPs with internal-mode interpretation.
    %
    % `IMInternalModes` translates standard `F` and `G` internal-mode
    % problems into the canonical scalar EVP. The solved scalar `u` is
    % either `F` or `G`; the other variable is recovered diagnostically by
    % the relation handles `FfromGz` and `GfromFz` on the resulting
    % `IMInternalModesBasis`.
    % Internal-mode EVPs own the stratification profile `N2` and physical
    % vertical domain used by solvers and basis sets.
    %
    % ```matlab
    % N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
    % solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.G(linspace(-4000,0,200).');
    % ```
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Summarize internal-mode EVPs
    % - Topic: Inspect internal-mode configuration
    % - Topic: Inspect internal-mode inner products
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModes < IMEigenvalueProblem

    properties (SetAccess = private)
        % Solved physical variable, `"F"` or `"G"`.
        %
        % `formulation` tells the canonical solver which variable is the
        % native unknown `u`. The complementary variable is evaluated
        % diagnostically by `IMInternalModesBasis`. Coefficient handles can
        % read this value as `ctx.formulation`.
        %
        % - Topic: Inspect internal-mode configuration
        formulation = "G"

        % Physical mode-family declaration.
        %
        % `modeFamily` tells internal-mode utilities which physical
        % catalog and coupled normalization rules are meaningful for this
        % EVP. The default `"none"` installs only generic internal-mode
        % behavior. The `"hydrostatic"` family declares the hydrostatic
        % `F`/`G` family, enabling the generalized boundary-condition
        % catalog and the coupled `geostrophic` normalization convention.
        % The `"meanDensityAnomaly"` family solves generalized-energy
        % `G` modes and carries a surface-referenced diagnostic `F`.
        %
        % - Topic: Inspect internal-mode configuration
        modeFamily = "none"

        % Buoyancy frequency squared function.
        %
        % `N2` has signature `values = N2(z)` and is owned by the EVP so
        % solvers can prepare `z`, `wkb`, or `density` coordinates from the
        % same continuous stratification.
        %
        % - Topic: Inspect internal-mode configuration
        N2

        % Coriolis parameter.
        %
        % `f0` is stored in `parameters.f0` and is available to coefficient
        % handles as `ctx.f0`.
        %
        % - Topic: Inspect internal-mode configuration
        f0 = 0

        % Gravitational acceleration.
        %
        % `g` is stored in `parameters.g` and is available to coefficient
        % handles as `ctx.g`.
        %
        % - Topic: Inspect internal-mode configuration
        g = 9.81

        % Equivalent-depth conversion function.
        %
        % `hFromEigenvalue` maps retained eigenvalues to equivalent depths
        % for internal-mode basis sets. The handle has signature
        % `h = hFromEigenvalue(lambda)`, so
        % $$h_j=\texttt{hFromEigenvalue}(\lambda_j).$$
        %
        % - Topic: Inspect internal-mode configuration
        hFromEigenvalue = @(lambda) 1 ./ lambda

        % Diagnostic relation from `G` derivative to `F`.
        %
        % `FfromGz` has signature `F = FfromGz(z,dGdz,h,ctx)`. The default
        % relation is
        % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
        %
        % - Topic: Inspect internal-mode configuration
        FfromGz = @(z,dGdz,h,ctx) dGdz .* reshape(h, 1, [])

        % Diagnostic relation from `F` derivative to `G`.
        %
        % `GfromFz` has signature `G = GfromFz(z,dFdz,h,ctx)`. The default
        % relation is the hydrostatic inverse
        % $$G_j(z)=-\frac{g}{N^2(z)}
        % \frac{\partial F_j}{\partial z}(z).$$
        % Wave-mode factories install relation handles with the
        % appropriate wave correction factors.
        %
        % - Topic: Inspect internal-mode configuration
        GfromFz = @(z,dFdz,h,ctx) -(ctx.g./ctx.N2(z(:))).*dFdz
    end

    methods
        function self = IMInternalModes(options)
            % Create an internal-mode canonical EVP.
            %
            % The constructor copies `options.parameters`, then writes the
            % internal-mode fields `f0`, `g`, and `formulation` into the
            % parameter struct. These constructor-owned fields override
            % same-named entries in `options.parameters`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.formulation: solved variable, `"F"` or `"G"`
            % - Parameter options.modeFamily: physical family, `"none"`, `"hydrostatic"`, or `"meanDensityAnomaly"`
            % - Parameter options.p: canonical derivative-flux coefficient
            % - Parameter options.q: canonical left-side value coefficient
            % - Parameter options.r: canonical metric coefficient
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.FfromGz: diagnostic relation from `G` derivative to `F`
            % - Parameter options.GfromFz: diagnostic relation from `F` derivative to `G`
            % - Parameter options.parameters: named coefficient parameters
            % - Returns evp: internal-mode EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "internalModes"
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.N2 function_handle
                options.formulation {mustBeTextScalar, mustBeMember(options.formulation, ["F", "G"])} = "G"
                options.modeFamily {mustBeTextScalar, mustBeMember(options.modeFamily, ["none", "hydrostatic", "meanDensityAnomaly"])} = "none"
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.hFromEigenvalue function_handle = @(lambda) 1 ./ lambda
                options.FfromGz function_handle = @(z,dGdz,h,ctx) dGdz .* reshape(h, 1, [])
                options.GfromFz function_handle = @(z,dFdz,h,ctx) -(ctx.g./ctx.N2(z(:))).*dFdz
                options.parameters struct = struct()
            end

            formulation = string(options.formulation);
            modeFamily = string(options.modeFamily);
            parameters = options.parameters;
            parameters.f0 = options.f0;
            parameters.g = options.g;
            parameters.formulation = formulation;

            self@IMEigenvalueProblem(name=options.name, p=options.p, q=options.q, r=options.r, ...
                zDomain=options.zDomain, surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                parameters=parameters);
            self.formulation = formulation;
            self.modeFamily = modeFamily;
            self.N2 = options.N2;
            self.f0 = options.f0;
            self.g = options.g;
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.FfromGz = options.FfromGz;
            self.GfromFz = options.GfromFz;
        end

        summarize(self, solver)

        function context = contextForSolver(self, solver)
            % Return the internal-mode coefficient context.
            %
            % The returned context extends the canonical context with `N2`,
            % `f0`, `g`, and `formulation`.
            %
            % - Topic: Developer topics
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            % - Developer: true
            arguments
                self IMInternalModes
                solver IMSolver
            end

            context = contextForSolver@IMEigenvalueProblem(self, solver);
            context.N2 = @(z) self.N2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.formulation = self.formulation;
        end

        function spec = innerProduct(self, variable)
            % Return the signed `F` or `G` inner-product recipe.
            %
            % The solved `G` variable uses the EVP's canonical weight `r`:
            % $$N^2/g$$ for hydrostatic modes, $$(N^2-f_0^2)/g$$ at fixed
            % wavenumber, and $$(N^2-\omega^2)/g$$ at fixed frequency.
            % Diagnostic hydrostatic `G` uses $$N^2/g$$; `F` retains unit
            % interior weight. The returned struct has fields
            % `variable`, `kind`, `interiorWeight`, `surfaceWeights`,
            % `bottomWeights`, `endpointInnerProductTerms`,
            % `hasInnerProduct`, and `reason`. `hasInnerProduct` is true
            % when the variable has a known inner product. When it is false, Gram
            % matrices, spectra, and inner-product normalization for that
            % variable throw `IMInternalModesBasis:UnavailableInnerProduct`.
            % Diagnostic variables use the value-only hydrostatic endpoint
            % catalog only when `modeFamily` is `"hydrostatic"` and a
            % catalog row is known; other diagnostic inner products are
            % unavailable until a family catalog is added. Endpoint
            % inner-product terms from the catalog have the form
            % $$\alpha_\ell F_i(z_\ell)F_j(z_\ell)$$ or
            % $$\alpha_\ell G_i(z_\ell)G_j(z_\ell),$$
            % where $$z_\ell$$ is the bottom or surface endpoint. The
            % variable used in the endpoint term is stored as
            % `term.variable`, so a `G` inner product can contain an
            % endpoint term involving `F`, and conversely.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Returns spec: struct with interior and endpoint inner-product terms
            arguments
                self IMInternalModes
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])} = self.formulation
            end

            variable = string(variable);
            spec.variable = variable;
            spec.kind = "signedPontryagin";
            if variable == "G"
                spec.interiorWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                spec.interiorWeight = @(z,~) ones(size(z));
            end
            spec.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();
            if variable == self.formulation
                if variable == "G"
                    spec.interiorWeight = self.r;
                end
                spec.surfaceWeights = self.endpointWeights("surface");
                spec.bottomWeights = self.endpointWeights("bottom");
                [spec.hasInnerProduct, spec.reason] = self.solvedInnerProductAvailability(spec.surfaceWeights, spec.bottomWeights);
            else
                spec.surfaceWeights = IMInternalModes.emptyEndpointWeights();
                spec.bottomWeights = IMInternalModes.emptyEndpointWeights();
                catalog = IMHydrostaticInnerProductCatalog.resolve(self, variable);
                spec.endpointInnerProductTerms = catalog.endpointInnerProductTerms;
                spec.hasInnerProduct = catalog.hasInnerProduct;
                spec.reason = catalog.reason;
            end
        end

        function spec = majorantInnerProduct(self, variable)
            % Return the induced positive Hilbert-majorant recipe.
            %
            % The signed Pontryagin product returned by `innerProduct` is
            % the physical pairing used for orthogonality and projection.
            % Its natural $$L^2\oplus\mathbb C^s$$ coordinate decomposition
            % induces a positive Hilbert product by taking the absolute
            % interior weight and the absolute value of every endpoint
            % coefficient:
            %
            % $$
            % (U,V)_+=\int |w|\,\overline{U}V\,dz+
            % \sum_\ell |\alpha_\ell|\,
            % \overline{L_\ell[U]}L_\ell[V].
            % $$
            %
            % Use this recipe for magnitudes, error tolerances, and
            % convergence diagnostics. Use `innerProduct` for signed
            % invariants, projection functionals, and modal coefficients.
            % The two recipes coincide when the interior weight and every
            % endpoint coefficient are nonnegative.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = majorantInnerProduct(evp,variable)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Returns spec: positive interior and absolute-endpoint recipe
            arguments
                self IMInternalModes
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])} = self.formulation
            end

            spec = self.innerProduct(variable);
            spec.kind = "inducedHilbertMajorant";
            signedWeight = spec.interiorWeight;
            spec.interiorWeight = @(z,ctx) abs(signedWeight(z,ctx));
            for iWeight = 1:numel(spec.surfaceWeights)
                spec.surfaceWeights(iWeight).coefficient = abs(spec.surfaceWeights(iWeight).coefficient);
            end
            for iWeight = 1:numel(spec.bottomWeights)
                spec.bottomWeights(iWeight).coefficient = abs(spec.bottomWeights(iWeight).coefficient);
            end
            for iTerm = 1:numel(spec.endpointInnerProductTerms)
                spec.endpointInnerProductTerms(iTerm).coefficient = abs(spec.endpointInnerProductTerms(iTerm).coefficient);
            end
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, modeSelectionDiagnostics)
            % Create an internal-mode basis set.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,modeSelectionDiagnostics)
            % - Returns basisSet: internal-mode basis set
            % - Developer: true
            arguments
                self IMInternalModes
                solver IMSolver
                nativeModes (:,:) double
                eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                modeNumber (1,:) double {mustBeInteger}
                modeSelectionDiagnostics struct
            end

            metadata = struct();
            if string(self.name) == "geostrophicAPVModes"
                metadata.g0 = self.parameters.g0;
                metadata.gd = self.parameters.gd;
                metadata.surfaceBoundary = string(self.parameters.surfaceBoundary);
            end
            basisSet = IMInternalModesBasis(solver=solver, evp=self, nativeModes=nativeModes, eigenvalues=eigenvalues, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics, metadata=metadata);
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G` internal-mode EVP.
            %
            % This factory creates the hydrostatic `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % =
            % \lambda_j\frac{N^2(z)}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid-lid and
            % rigid-bottom conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % Physical hydrostatic endpoint laws written in `F` and `G`
            % can be converted with `IMHydrostaticBoundaryCondition`
            % before they are passed to this factory:
            %
            % ```matlab
            % law = IMHydrostaticBoundaryCondition(a=A/g,b=1);
            % surfaceBoundary = law.canonicalBoundary(formulation="G",g=g);
            % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g,surfaceBoundary=surfaceBoundary);
            % ```
            %
            % After conversion, `innerProduct("F")` and
            % `innerProduct("G")` use the hydrostatic endpoint catalog to
            % report which bilinear forms are known.
            % Solved hydrostatic basis sets install the `geostrophic`
            % normalization rule and use it by default because they set
            % `modeFamily` to `"hydrostatic"`. This factory sets
            % `parameters.formulation`, `parameters.f0`, and `parameters.g`.
            %
            % ```matlab
            % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
            % solver = IMSolverSpectral(nEVP=128);
            % basisSet = solver.solveEVP(evp,nModes=4);
            % G = basisSet.G(z);
            % ```
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticGModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: hydrostatic `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            evp = IMInternalModes(name="hydrostaticGModes", formulation="G", modeFamily="hydrostatic", ...
                N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g, f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = hydrostaticFModes(options)
            % Create the hydrostatic `F` internal-mode EVP.
            %
            % This factory creates the hydrostatic `F`-form problem
            %
            % $$
            % -\frac{\partial}{\partial z}
            % \left(
            % \frac{1}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z)
            % \right)
            % =
            % \lambda_j\frac{F_j(z)}{g},
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell F_j(z_\ell)
            % -b_\ell\frac{1}{N^2(z_\ell)}
            % \frac{\partial F_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell F_j(z_\ell)
            % -d_\ell\frac{1}{N^2(z_\ell)}
            % \frac{\partial F_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.neumann()`, giving
            %
            % $$
            % \frac{1}{N^2(z_s)}
            % \frac{\partial F_j}{\partial z}(z_s)=0,\qquad
            % \frac{1}{N^2(z_b)}
            % \frac{\partial F_j}{\partial z}(z_b)=0.
            % $$
            %
            % Through the hydrostatic relation
            %
            % $$
            % G_j(z)=-\frac{g}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z),
            % $$
            %
            % these are the same rigid-lid and rigid-bottom conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % Physical hydrostatic endpoint laws written in `F` and `G`
            % can be converted with `IMHydrostaticBoundaryCondition`
            % before they are passed to this factory:
            %
            % ```matlab
            % law = IMHydrostaticBoundaryCondition(b=B,c=C);
            % surfaceBoundary = law.canonicalBoundary(formulation="F",g=g);
            % evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain,g=g,surfaceBoundary=surfaceBoundary);
            % ```
            %
            % After conversion, `innerProduct("F")` and
            % `innerProduct("G")` use the hydrostatic endpoint catalog to
            % report which bilinear forms are known.
            % The barotropic zero mode is inferred from the canonical left
            % problem during mode selection.
            % This factory sets `parameters.formulation` and `parameters.g`;
            % `parameters.f0` is supplied by the internal-mode constructor
            % default.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticFModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
            end

            evp = IMInternalModes(name="hydrostaticFModes", formulation="F", modeFamily="hydrostatic", ...
                N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = geostrophicAPVModes(options)
            % Create signed generalized-energy geostrophic APV modes.
            %
            % This factory creates the hydrostatic `F` problem
            % $$-\frac{\partial}{\partial z}\left(\frac{1}{N^2}\frac{\partial F_j}{\partial z}\right)=\frac{F_j}{g h_j}$$
            % with diagnostic relation
            % $$G_j=-\frac{g}{N^2}\frac{\partial F_j}{\partial z}.$$
            % The endpoint accelerations `g0` and `gd` are required and
            % may be signed finite values, zero, or positive infinity.
            % Zero selects Dirichlet; positive infinity selects the
            % corresponding reciprocal-zero endpoint limit.
            %
            % The surface convention maps to canonical ordinary boundary
            % conditions as follows:
            %
            % | Convention | `g0` | Canonical `(a,b)` | Physical condition |
            % | --- | --- | --- | --- |
            % | free surface | finite nonzero | $$(-(1/g+1/g_0),1)$$ | $$G_s=(1+g/g_0)F_s$$ |
            % | free surface | zero | Dirichlet | $$F_s=0$$ |
            % | free surface | `Inf` | $$(-1/g,1)$$ | $$G_s=F_s$$ |
            % | rigid lid | finite nonzero | $$(-1/g_0,1)$$ | $$G_s=(g/g_0)F_s$$ |
            % | rigid lid | zero | Dirichlet | $$F_s=0$$ |
            % | rigid lid | `Inf` | Neumann | $$G_s=0$$ |
            %
            % The bottom condition is:
            %
            % | `gd` | Canonical `(a,b)` | Physical condition |
            % | --- | --- | --- |
            % | finite nonzero | $$(1/g_d,1)$$ | $$G_b=-(g/g_d)F_b$$ |
            % | zero | Dirichlet | $$F_b=0$$ |
            % | `Inf` | Neumann | $$G_b=0$$ |
            %
            % The signed generalized-energy inner products are
            % $$\langle F_i,F_j\rangle_F=\int_{z_b}^{z_s}F_iF_j\,dz$$
            % and
            % $$\langle G_i,G_j\rangle_G=\int_{z_b}^{z_s}\frac{N^2}{g}G_iG_j\,dz+c_sG_i(z_s)G_j(z_s)+\frac{g_d}{g}G_i(z_b)G_j(z_b).$$
            % For finite active endpoints, $$c_s=g_0/(g+g_0)$$
            % under the free-surface convention and $$c_s=g_0/g$$
            % under the rigid-lid convention. Dirichlet, Neumann,
            % positive-infinite, and the free-surface `g0=-g` reductions
            % omit any coefficient made singular by the reduced endpoint
            % constraint. The free-surface `g0=Inf` limit instead has
            % $$c_s=1$$.
            %
            % Solved bases use `Normalization.depth` by default, so
            % $$D^{-1}\int_{z_b}^{z_s}F_j^2\,dz=1$$ without endpoint
            % terms. The same positive factor scales `F` and `G`, including
            % negative-eigendepth modes. Signed `h` values are retained;
            % an exact zero eigenvalue is represented by `h=Inf`.
            % `g0`, `gd`, and the surface convention are copied into
            % `evp.parameters` and basis metadata.
            % `surfaceBoundary` is the endpoint-convention contract used
            % to identify the matching zero-APV family. `N2(z)` must
            % return values with the shape of `z`; evaluated basis methods
            % return one column per retained mode.
            %
            % ```matlab
            % evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-4000 0],g0=-0.02,gd=Inf);
            % basisSet = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=4);
            % F = basisSet.F(z);
            % ```
            %
            % ```matlab
            % evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-4000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            % ```
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.geostrophicAPVModes(options)
            % - Parameter options.N2: buoyancy frequency squared function in radians squared per second squared
            % - Parameter options.zDomain: physical vertical domain in meters
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.g0: signed surface acceleration in meters per second squared
            % - Parameter options.gd: signed bottom acceleration in meters per second squared
            % - Parameter options.surfaceBoundary: `"freeSurface"` or `"rigidLid"`
            % - Returns evp: generalized-energy geostrophic APV EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
                options.surfaceBoundary {mustBeTextScalar, mustBeMember(options.surfaceBoundary, ["freeSurface", "rigidLid"])} = "freeSurface"
            end

            if isnan(options.g0) || options.g0 == -Inf
                error("IMInternalModes:InvalidSurfaceAcceleration", "g0 must be signed finite, zero, or positive Inf.");
            end
            if isnan(options.gd) || options.gd == -Inf
                error("IMInternalModes:InvalidBottomAcceleration", "gd must be signed finite, zero, or positive Inf.");
            end

            surfaceConvention = string(options.surfaceBoundary);
            if options.g0 == 0
                surfaceCondition = IMBoundaryCondition.dirichlet();
            elseif surfaceConvention == "freeSurface"
                if isinf(options.g0)
                    surfaceCondition = IMBoundaryCondition(a=-1/options.g, b=1);
                else
                    surfaceCondition = IMBoundaryCondition(a=-(1/options.g + 1/options.g0), b=1);
                end
            elseif isinf(options.g0)
                surfaceCondition = IMBoundaryCondition.neumann();
            else
                surfaceCondition = IMBoundaryCondition(a=-1/options.g0, b=1);
            end

            if options.gd == 0
                bottomCondition = IMBoundaryCondition.dirichlet();
            elseif isinf(options.gd)
                bottomCondition = IMBoundaryCondition.neumann();
            else
                bottomCondition = IMBoundaryCondition(a=1/options.gd, b=1);
            end

            parameters = struct("g0", options.g0, "gd", options.gd, "surfaceBoundary", surfaceConvention);
            evp = IMInternalModes(name="geostrophicAPVModes", formulation="F", modeFamily="hydrostatic", N2=options.N2, zDomain=options.zDomain, p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, surfaceBoundary=surfaceCondition, bottomBoundary=bottomCondition, parameters=parameters);
        end

        function evp = geostrophicGeneralizedPotentialEnstrophyModes(options)
            % Create free-surface generalized-potential-enstrophy modes.
            %
            % At fixed positive horizontal wavenumber `k`, this factory
            % creates the `F` problem
            % $$
            % -\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\frac{\partial F_j}{\partial z}\right)
            % +k^2F_j=\Lambda_jF_j.
            % $$
            % A finite surface weight applies
            % $$
            % \frac{f_0^2}{N_s^2}F_j'(z_s)+\frac{f_0^2}{g}F_j(z_s)
            % =\Lambda_j\frac{f_0^2}{\alpha_0}F_j(z_s),
            % $$
            % while a finite bottom weight applies
            % $$
            % \frac{f_0^2}{N_b^2}F_j'(z_b)
            % =-\Lambda_j\frac{f_0^2}{\alpha_d}F_j(z_b).
            % $$
            % Positive infinity removes the corresponding eigenvalue-side
            % endpoint term. The default bottom is inactive. When `alpha0`
            % is omitted, the surface weight is
            % $$
            % \alpha_0=\frac{f_0^2}{b_\mathrm{eff}},\qquad
            % b_\mathrm{eff}=\frac{(\int N\,dz)^2}{4\int N^2\,dz}.
            % $$
            %
            % Solved bases use
            % `Normalization.generalizedPotentialEnstrophy` by default, so
            % $$
            % \frac{1}{D}\left[\int F_iF_j\,dz
            % +\frac{f_0^2}{\alpha_0}F_i(z_s)F_j(z_s)
            % +\frac{f_0^2}{\alpha_d}F_i(z_b)F_j(z_b)\right]=\delta_{ij},
            % $$
            % with inactive terms omitted. The equivalent depth is the
            % diagnostic quantity
            % $$h_j=f_0^2/[g(\Lambda_j-k^2)].$$
            %
            % ```matlab
            % evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
            %     N2=N2,zDomain=[-4000 0],k=2*pi/100e3,f0=1e-4);
            % basisSet = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=8);
            % ```
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(options)
            % - Parameter options.N2: buoyancy frequency squared function in radians squared per second squared
            % - Parameter options.zDomain: physical vertical domain in meters
            % - Parameter options.k: positive horizontal wavenumber in radians per meter
            % - Parameter options.f0: nonzero Coriolis parameter in radians per second
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.alpha0: optional positive surface generalized-potential-enstrophy weight
            % - Parameter options.alphaD: positive bottom generalized-potential-enstrophy weight or `Inf`
            % - Returns evp: generalized-potential-enstrophy mode EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.k (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.alpha0 double {mustBeReal} = zeros(0,1)
                options.alphaD (1,1) double {mustBeReal} = Inf
            end

            if options.f0 == 0
                error("IMInternalModes:InvalidCoriolisParameter", "f0 must be nonzero.");
            end
            if ~isempty(options.alpha0) && (~isscalar(options.alpha0) || isnan(options.alpha0) || options.alpha0 <= 0)
                error("IMInternalModes:InvalidSurfaceEnstrophyWeight", "alpha0 must be empty, positive finite, or positive Inf.");
            end
            if isnan(options.alphaD) || options.alphaD <= 0
                error("IMInternalModes:InvalidBottomEnstrophyWeight", "alphaD must be positive finite or positive Inf.");
            end

            zDomain = sort(options.zDomain);
            try
                I1 = integral(@(z) sqrt(options.N2(z)),zDomain(1),zDomain(2));
                I2 = integral(options.N2,zDomain(1),zDomain(2));
            catch cause
                exception = MException("IMInternalModes:InvalidStratificationScale", "N2 could not be integrated to determine the generalized-potential-enstrophy scale.");
                throw(addCause(exception,cause))
            end
            if ~isscalar(I1) || ~isreal(I1) || ~isfinite(I1) || I1 <= 0 ...
                    || ~isscalar(I2) || ~isreal(I2) || ~isfinite(I2) || I2 <= 0
                error("IMInternalModes:InvalidStratificationScale", "The integrals of N and N2 must be finite positive real scalars.");
            end
            bEffective = I1*I1/(4*I2);
            if isempty(options.alpha0)
                alpha0 = options.f0*options.f0/bEffective;
                alpha0Source = "stratification";
            else
                alpha0 = options.alpha0;
                alpha0Source = "explicit";
            end

            f02 = options.f0*options.f0;
            if isinf(alpha0)
                surfaceCondition = IMBoundaryCondition(a=-f02/options.g,b=1);
            else
                surfaceCondition = IMBoundaryCondition(a=-f02/options.g,b=1,c=f02/alpha0);
            end
            if isinf(options.alphaD)
                bottomCondition = IMBoundaryCondition.neumann();
            else
                bottomCondition = IMBoundaryCondition(a=0,b=1,c=-f02/options.alphaD);
            end

            k = options.k;
            f0 = options.f0;
            g = options.g;
            parameters = struct("k",k,"alpha0",alpha0,"alphaD",options.alphaD, ...
                "bEffective",bEffective,"alpha0Source",alpha0Source);
            evp = IMInternalModes(name="geostrophicGeneralizedPotentialEnstrophyModes", ...
                formulation="F",modeFamily="none",N2=options.N2,zDomain=zDomain, ...
                p=@(z,ctx) ctx.f0*ctx.f0./ctx.N2(z), ...
                q=@(z,ctx) ctx.k*ctx.k*ones(size(z)),r=@(z,~) ones(size(z)), ...
                surfaceBoundary=surfaceCondition,bottomBoundary=bottomCondition, ...
                f0=f0,g=g,hFromEigenvalue=@(lambda) f0*f0./(g*(lambda-k*k)),parameters=parameters);
        end

        function evp = meanDensityAnomalyModes(options)
            % Create generalized-energy mean-density-anomaly modes.
            %
            % This factory creates the `G`-form problem
            %
            % $$
            % -G_j''(z)=\frac{N^2(z)}{g h_j}G_j(z)
            % $$
            %
            % with endpoint conditions
            %
            % $$
            % g h_jG_j'(z_s)=g_0G_j(z_s),\qquad
            % g h_jG_j'(z_b)=-g_dG_j(z_b).
            % $$
            %
            % `g0` and `gd` are required. A finite value, including zero,
            % keeps that endpoint active; zero is the Neumann limit.
            % Positive infinity imposes Dirichlet data and omits the
            % corresponding generalized-energy endpoint term. `NaN` and
            % negative infinity are rejected.
            %
            % Solved modes use the signed generalized-energy normalization
            %
            % $$
            % \frac{1}{g}\int_{z_b}^{z_s}N^2G_iG_j\,dz
            % +\frac{g_0}{g}G_i(z_s)G_j(z_s)
            % +\frac{g_d}{g}G_i(z_b)G_j(z_b)
            % =\epsilon_j\delta_{ij},
            % $$
            %
            % with inactive terms omitted. The basis exposes
            % `basisSet.signatures` as $$\epsilon_j\in\{-1,+1\}$$.
            % For the continuous projection functional
            % $$\mathcal G_j[X]=\langle G_j,X\rangle_G$$, normalized
            % coefficients obey $$A_j=\epsilon_j\mathcal G_j[X]$$.
            % Its aligned diagnostic pressure modes are computed by
            % surface-referenced integration,
            %
            % $$
            % F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz',
            % \qquad F_j(z_s)=0.
            % $$
            %
            % Discrete transforms directly project the `G` channel and
            % synthesize both `G` and `F`. The diagnostic `F` channel does
            % not define an independent coefficient projection metric.
            %
            % ```matlab
            % evp = IMInternalModes.meanDensityAnomalyModes( ...
            %     N2=N2,zDomain=[-4000 0],g0=0.02,gd=Inf);
            % basisSet = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=8);
            % F = basisSet.F(z);
            % G = basisSet.G(z);
            % ```
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.meanDensityAnomalyModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.g0: signed finite surface acceleration, zero, or positive infinity
            % - Parameter options.gd: signed finite bottom acceleration, zero, or positive infinity
            % - Returns evp: generalized-energy mean-density-anomaly EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
            end

            evp = IMMeanDensityAnomalyModes(N2=options.N2,zDomain=options.zDomain,g=options.g,g0=options.g0,gd=options.gd);
        end

        function evp = waveModesAtWavenumber(options)
            % Create the fixed-wavenumber wave-mode EVP.
            %
            % This factory creates the fixed-wavenumber `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % +k^2G_j(z)
            % =
            % \lambda_j\frac{N^2(z)-f_0^2}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid endpoint
            % conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % A linear free-surface condition at the surface can be written
            % as
            %
            % $$
            % G_j(z_s)=h_j\frac{\partial G_j}{\partial z}(z_s),
            % \qquad \lambda_j=\frac{1}{h_j},
            % $$
            %
            % equivalently
            %
            % $$
            % \frac{\partial G_j}{\partial z}(z_s)
            % =
            % \lambda_j G_j(z_s).
            % $$
            %
            % In canonical boundary-condition coefficients this is
            % `IMBoundaryCondition(a=0,b=1,c=1,d=0)` at the surface.
            % Solved fixed-wavenumber basis sets install the `kConstant`
            % normalization rule and use it by default.
            % This factory adds `parameters.k` and sets
            % `parameters.formulation`, `parameters.f0`, and `parameters.g`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtWavenumber(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: fixed-wavenumber `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.k (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            k = options.k;
            parameters = struct("k", k);
            GfromFz = @(z,dFdz,h,ctx) -(ctx.g./(ctx.N2(z(:)) - ctx.f0*ctx.f0 - ctx.g*reshape(h,1,[])*k*k)).*dFdz;
            evp = IMInternalModes(name="waveModesAtWavenumber", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) k*k*ones(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                GfromFz=GfromFz, parameters=parameters);
        end

        function evp = waveModesAtFrequency(options)
            % Create the fixed-frequency wave-mode EVP.
            %
            % This factory creates the fixed-frequency `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % =
            % \lambda_j\frac{N^2(z)-\omega^2}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid endpoint
            % conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % A linear free-surface condition at the surface can be written
            % as
            %
            % $$
            % G_j(z_s)=h_j\frac{\partial G_j}{\partial z}(z_s),
            % \qquad \lambda_j=\frac{1}{h_j},
            % $$
            %
            % equivalently
            %
            % $$
            % \frac{\partial G_j}{\partial z}(z_s)
            % =
            % \lambda_j G_j(z_s).
            % $$
            %
            % In canonical boundary-condition coefficients this is
            % `IMBoundaryCondition(a=0,b=1,c=1,d=0)` at the surface.
            % Solved fixed-frequency basis sets use the generic `unity`
            % normalization by default. A fixed-frequency diagnostic `F`
            % inner-product normalization is deferred until the wave
            % diagnostic inner-product catalog is derived. This factory
            % adds `parameters.omega` and sets `parameters.formulation`,
            % `parameters.f0`, and `parameters.g`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtFrequency(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.omega: wave frequency
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: fixed-frequency `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.omega (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            omega = options.omega;
            parameters = struct("omega", omega);
            GfromFz = @(z,dFdz,h,ctx) -(ctx.g./(ctx.N2(z(:)) - omega*omega)).*dFdz;
            evp = IMInternalModes(name="waveModesAtFrequency", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                GfromFz=GfromFz, parameters=parameters);
        end
    end

    methods (Access = private)
        function [hasInnerProduct, reason] = solvedInnerProductAvailability(self, surfaceWeights, bottomWeights)
            surfaceActive = self.surfaceBoundary.isEigenvalueDependent();
            bottomActive = self.bottomBoundary.isEigenvalueDependent();
            if (surfaceActive && isempty(surfaceWeights)) || (bottomActive && isempty(bottomWeights))
                hasInnerProduct = false;
                reason = "At least one active boundary condition has a degenerate or unavailable endpoint metric weight.";
            elseif isempty(surfaceWeights) && isempty(bottomWeights)
                hasInnerProduct = true;
                reason = "The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.";
            else
                hasInnerProduct = true;
                reason = "The solved formulation inner product follows the canonical scalar EVP endpoint weights.";
            end
        end
    end

    methods (Static, Access = private)
        function weights = emptyEndpointWeights()
            weights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
        end

    end
end
