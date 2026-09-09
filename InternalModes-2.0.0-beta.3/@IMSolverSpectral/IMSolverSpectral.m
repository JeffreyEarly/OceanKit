classdef IMSolverSpectral < IMSolver
    % Solve physical-coordinate EVPs with a Chebyshev spectral discretization.
    %
    % `IMSolverSpectral` owns the numerical coordinate choice, Chebyshev
    % resolution, derivative matrices, and physical-coordinate pullback
    % rules. It is configured against an EVP or geostrophic zero-APV problem
    % before solving. Stretched-coordinate values, inverse maps, Jacobians,
    % and second derivatives use one smooth Chebfun representation, so
    % reconstruction and sampled-field calculus share the same coordinate.
    %
    % ```matlab
    % evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], k=1e-4);
    % solver = IMSolverSpectral(nEVP=64);
    % basisSet = solver.solveEVP(evp);
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Inspect solvers
    % - Topic: Solve EVPs
    % - Topic: Solve geostrophic zero-APV modes
    % - Topic: Assemble EVPs
    % - Topic: Evaluate native modes
    % - Topic: Developer topics
    % - Declaration: classdef IMSolverSpectral

    properties (SetAccess = private)
        % Number of native EVP coefficients.
        %
        % - Topic: Inspect solvers
        nEVP

        % Physical vertical domain.
        %
        % - Topic: Inspect solvers
        zDomain = [NaN NaN]

        % Native coordinate kind.
        %
        % `coordinateKind` is `"z"`, `"wkb"`, or `"density"`. The `"wkb"`
        % and `"density"` coordinates require an `IMInternalModes` EVP or
        % `IMGeostrophicZeroAPVModes` problem because their coordinate maps
        % use the problem-owned `N2` profile.
        %
        % - Topic: Inspect solvers
        coordinateKind

        % Native Lobatto grid.
        %
        % - Topic: Inspect solvers
        xNative

        % Physical points corresponding to `xNative`.
        %
        % - Topic: Inspect solvers
        zNative

        % Chebyshev basis matrix on the native grid.
        %
        % - Topic: Developer topics
        % - Developer: true
        T

        % Native first-derivative basis matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Tx

        % Native second-derivative basis matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Txx
    end

    properties (SetAccess = private)
        % Reference physical grid for inspecting the coordinate map.
        %
        % - Topic: Developer topics
        % - Developer: true
        zReference

        % Reference native coordinate grid.
        %
        % - Topic: Developer topics
        % - Developer: true
        xReference

        % Reference coordinate derivative $$dx/dz$$.
        %
        % - Topic: Developer topics
        % - Developer: true
        qReference

        % Reference physical derivative of $$dx/dz$$.
        %
        % - Topic: Developer topics
        % - Developer: true
        qzReference
    end

    properties (Access = private)
        coordinateMap_ = []
        coordinateDerivative_ = []
        coordinateSecondDerivative_ = []
    end

    properties (Access = protected)
        % EVP-provided buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function = []
    end

    methods
        function self = IMSolverSpectral(options)
            % Create a coordinate-aware spectral solver.
            %
            % The `"z"` coordinate works for any canonical EVP. The `"wkb"`
            % and `"density"` coordinates use the problem stratification
            % `N2`, so they can only be configured with `IMInternalModes`
            % EVPs or `IMGeostrophicZeroAPVModes` problems.
            %
            % - Topic: Create solvers
            % - Declaration: solver = IMSolverSpectral(options)
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.coordinateKind: native coordinate kind
            % - Returns solver: initialized spectral solver
            arguments
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP, 4)} = 64
                options.coordinateKind {mustBeTextScalar, mustBeMember(options.coordinateKind, ["z", "wkb", "density"])} = "z"
            end

            self.nEVP = options.nEVP;
            self.coordinateKind = string(options.coordinateKind);
        end

        function solver = configuredForEVP(self, evp)
            % Return a spectral solver configured for an EVP.
            %
            % - Topic: Assemble EVPs
            % - Topic: Developer topics
            % - Declaration: solver = configuredForEVP(solver,evp)
            % - Parameter evp: canonical EVP descriptor
            % - Returns solver: solver with grid and coordinate matrices initialized
            % - Developer: true
            arguments
                self IMSolverSpectral
                evp IMEigenvalueProblem
            end

            solver = self;
            solver.zDomain = evp.zDomain;
            if solver.coordinateKind == "z"
                solver.N2Function = [];
            elseif isa(evp, "IMInternalModes")
                solver.N2Function = @(z) evp.N2(z);
            else
                error("IMEigenvalueProblem:UnsupportedCoordinateKind", ...
                    "Only internal-mode EVPs support coordinateKind=""%s"".", solver.coordinateKind);
            end
            solver = solver.setupCoordinate();
            solver = solver.setupNativeGrid();
        end

        function solver = configuredForGeostrophicZeroAPVModes(self, problem)
            % Return a spectral solver configured for zero-APV modes.
            %
            % - Topic: Solve geostrophic zero-APV modes
            % - Topic: Developer topics
            % - Declaration: solver = configuredForGeostrophicZeroAPVModes(solver,problem)
            % - Parameter problem: geostrophic zero-APV problem
            % - Returns solver: solver with zero-APV grid and coordinate matrices initialized
            % - Developer: true
            arguments
                self IMSolverSpectral
                problem IMGeostrophicZeroAPVModes
            end

            solver = self;
            solver.zDomain = problem.zDomain;
            solver.N2Function = @(z) problem.N2(z);
            solver = solver.setupCoordinate();
            solver = solver.setupNativeGrid();
        end

        function context = context(self)
            % Return the framework coefficient context.
            %
            % Solvers provide discretization fields. EVPs add the physical
            % domain, medium, and constants.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: framework coefficient context
            context.coordinateKind = self.coordinateKind;
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = N2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            if isempty(self.N2Function)
                error("IMSolverSpectral:UnsupportedOperation", ...
                    "This spectral solver is not configured with an N2 function.");
            end
            values = self.N2Function(z);
        end

        function values = differentiateGridValues(self, values, derivativeOrder)
            % Differentiate values sampled on the native grid.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: values = differentiateGridValues(solver,values,derivativeOrder)
            % - Parameter values: one value per native grid point
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: differentiated grid values
            values = values(:,:);
            if size(values,1) ~= self.nEVP
                error("IMSolverSpectral:InvalidGridValues", ...
                    "Grid values must have one row per native grid point.");
            end
            if derivativeOrder == 0
                return;
            end
            coefficients = self.T \ values;
            values = self.physicalDerivativeMatrix(derivativeOrder)*coefficients;
        end

        function D = physicalDerivativeMatrix(self, derivativeOrder)
            % Return a native matrix for a physical derivative.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: D = physicalDerivativeMatrix(solver,derivativeOrder)
            % - Parameter derivativeOrder: physical derivative order
            % - Returns D: matrix mapping coefficients to derivative values
            q = self.qAtZ(self.zNative);
            qz = self.qzAtZ(self.zNative);
            switch derivativeOrder
                case 0
                    D = self.T;
                case 1
                    D = diag(q)*self.Tx;
                case 2
                    D = diag(q.*q)*self.Txx + diag(qz)*self.Tx;
                otherwise
                    error("IMSolverSpectral:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function index = boundaryIndex(self, location)
            % Return the native row for a physical boundary.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: index = boundaryIndex(solver,location)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Returns index: native row index
            switch string(location)
                case "surface"
                    index = 1;
                case "bottom"
                    index = self.nEVP;
                otherwise
                    error("IMSolverSpectral:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"" or ""bottom"".");
            end
        end

        function values = evaluateNativeModes(self, nativeModes, z)
            % Evaluate native Chebyshev coefficient columns at physical points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluateNativeModes(solver,nativeModes,z)
            % - Parameter nativeModes: Chebyshev coefficient columns
            % - Parameter z: physical evaluation points
            % - Returns values: mode values at `z`
            z = z(:);
            nPolynomials = size(nativeModes,1);
            if nPolynomials == self.nEVP && self.isNativePhysicalGrid(z)
                TOut = self.T;
            else
                x = self.clampNativeCoordinate(self.xOfZ(z));
                TOut = self.chebyshevPolynomialsAtNativePoints(x,nPolynomials);
            end
            values = TOut*nativeModes;
        end

        function values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
            % Evaluate physical derivatives of native modes.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluatePhysicalDerivative(solver,nativeModes,z,derivativeOrder)
            % - Parameter nativeModes: Chebyshev coefficient columns
            % - Parameter z: physical evaluation points
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: derivative values at `z`
            z = z(:);
            nPolynomials = size(nativeModes,1);
            if nPolynomials == self.nEVP && self.isNativePhysicalGrid(z)
                TOut = self.T;
                TxOut = self.Tx;
                TxxOut = self.Txx;
            else
                x = self.clampNativeCoordinate(self.xOfZ(z));
                [TOut, TxOut, TxxOut] = self.chebyshevPolynomialsAtNativePoints(x,nPolynomials);
            end
            q = self.qAtZ(z);
            qz = self.qzAtZ(z);
            switch derivativeOrder
                case 0
                    values = TOut*nativeModes;
                case 1
                    values = diag(q)*TxOut*nativeModes;
                case 2
                    values = (diag(q.*q)*TxxOut + diag(qz)*TxOut)*nativeModes;
                otherwise
                    error("IMSolverSpectral:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function nativeIntegral = integrateGridValuesFromSurface(self, gridValues)
            % Integrate physical grid values downward from the surface.
            %
            % For physical integrands $$f(z)$$ sampled on `zNative`, this
            % method returns native Chebyshev coefficient columns for
            %
            % $$
            % I(z)=\int_z^{z_s}f(z')\,dz',
            % \qquad I(z_s)=0.
            % $$
            %
            % The Chebyshev series is integrated in the active native
            % coordinate $$x$$ after applying the Jacobian
            % $$dz/dx=1/q$$, where $$q=dx/dz$$. The degree-increasing final
            % coefficient is retained, and `evaluateNativeModes` evaluates
            % the resulting representation at its exact degree.
            %
            % - Topic: Developer topics
            % - Declaration: nativeIntegral = integrateGridValuesFromSurface(solver,gridValues)
            % - Parameter gridValues: physical integrand columns sampled on `zNative`
            % - Returns nativeIntegral: native Chebyshev coefficient columns for the surface-referenced integrals
            % - Developer: true
            arguments
                self IMSolverSpectral
                gridValues (:,:) double {mustBeReal, mustBeFinite}
            end

            if size(gridValues,1) ~= self.nEVP
                error("IMSolverSpectral:InvalidGridValues", "Grid values must have one row per native grid point.");
            end

            q = self.qAtZ(self.zNative);
            integrandCoefficients = self.T \ (gridValues./q(:));
            nativeIntegral = zeros(size(integrandCoefficients,1)+1,size(integrandCoefficients,2));
            coordinateScale = (max(self.xNative)-min(self.xNative))/2;
            for iColumn = 1:size(integrandCoefficients,2)
                primitive = coordinateScale*InternalModesSpectral.IntegrateChebyshevVector(integrandCoefficients(:,iColumn));
                surfaceValue = sum(primitive);
                primitive = -primitive;
                primitive(1) = primitive(1)+surfaceValue;
                nativeIntegral(:,iColumn) = primitive;
            end
        end

        function zRoots = rootsOfNativeMode(self, nativeMode)
            % Return physical roots of a native Chebyshev mode.
            %
            % The native mode is the Chebyshev series
            %
            % $$
            % u(x)=\sum_{k=0}^{K}c_kT_k(\widehat{x}),
            % $$
            %
            % where $$\widehat{x}\in[-1,1]$$ is the affine image of the
            % solver coordinate. Roots are eigenvalues of the Chebyshev
            % colleague matrix, filtered by their series residual, and then
            % mapped to physical depth by inverting `xOfZ`.
            %
            % - Topic: Evaluate native modes
            % - Declaration: zRoots = rootsOfNativeMode(solver,nativeMode)
            % - Parameter nativeMode: one native Chebyshev coefficient column
            % - Returns zRoots: sorted physical roots in the closed domain
            % - Developer: true
            arguments
                self IMSolverSpectral
                nativeMode (:,1) double {mustBeReal, mustBeFinite}
            end

            if length(nativeMode) ~= self.nEVP
                error("IMSolverSpectral:InvalidNativeMode", ...
                    "nativeMode must contain one coefficient for each of the %d native basis functions.", self.nEVP);
            end
            if any(~isfinite(self.zDomain)) || isempty(self.xNative)
                error("IMSolverSpectral:UnconfiguredSolver", ...
                    "Configure the spectral solver for an EVP before finding mode roots.");
            end

            coefficients = nativeMode(:);
            coefficientScale = norm(coefficients, Inf);
            if coefficientScale == 0
                error("IMSolverSpectral:InvalidNativeMode", "The native mode is identically zero.");
            end
            coefficients = coefficients/coefficientScale;
            coefficientTolerance = 100*length(coefficients)*eps;
            degree = find(abs(coefficients) > coefficientTolerance, 1, "last") - 1;
            if isempty(degree) || degree == 0
                zRoots = zeros(0,1);
                return
            end
            coefficients = coefficients(1:(degree+1));

            if degree == 1
                rootsNormalized = -coefficients(1)/coefficients(2);
            else
                colleague = zeros(degree);
                colleague(1,2) = 1;
                for iRow = 2:(degree-1)
                    colleague(iRow,iRow-1) = 0.5;
                    colleague(iRow,iRow+1) = 0.5;
                end
                colleague(end,:) = -coefficients(1:degree).'/(2*coefficients(end));
                colleague(end,end-1) = colleague(end,end-1) + 0.5;
                rootsNormalized = eig(colleague);
            end

            rootTolerance = max(1e-12, 100*degree*degree*eps);
            isReal = abs(imag(rootsNormalized)) <= rootTolerance.*max(1,abs(real(rootsNormalized)));
            rootsNormalized = real(rootsNormalized(isReal));
            rootsNormalized = rootsNormalized(rootsNormalized >= -1-rootTolerance & rootsNormalized <= 1+rootTolerance);
            rootsNormalized = min(max(rootsNormalized,-1),1);
            if isempty(rootsNormalized)
                zRoots = zeros(0,1);
                return
            end

            polynomialValues = zeros(size(rootsNormalized));
            theta = acos(rootsNormalized);
            for iCoefficient = 0:degree
                polynomialValues = polynomialValues + coefficients(iCoefficient+1)*cos(iCoefficient*theta);
            end
            residualTolerance = max(1e-11, 1000*degree*eps*norm(coefficients,1));
            rootsNormalized = rootsNormalized(abs(polynomialValues) <= residualTolerance);
            rootsNormalized = sort(rootsNormalized(:));
            if isempty(rootsNormalized)
                zRoots = zeros(0,1);
                return
            end

            duplicateTolerance = max(1e-11, 10*rootTolerance);
            rootsNormalized = rootsNormalized([true; diff(rootsNormalized) > duplicateTolerance]);
            xMin = min(self.xNative);
            xMax = max(self.xNative);
            xRoots = 0.5*(xMax-xMin)*rootsNormalized + 0.5*(xMax+xMin);
            zRoots = zeros(size(xRoots));
            for iRoot = 1:length(xRoots)
                zRoots(iRoot) = fzero(@(z) self.xOfZ(z)-xRoots(iRoot), self.zDomain);
            end
            zRoots = sort(zRoots);
            zTolerance = max(100*eps(max(1,max(abs(self.zDomain)))), 1e-10*max(1,diff(self.zDomain)));
            zRoots = min(max(zRoots,self.zDomain(1)),self.zDomain(2));
            zRoots = zRoots([true; diff(zRoots) > zTolerance]);
        end

        function z = innerProductGrid(self, zBounds)
            % Return the native grid used for spectral inner products.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: z = innerProductGrid(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: physical points corresponding to the native grid
            arguments
                self IMSolverSpectral
                zBounds (1,2) double
            end

            z = self.zNative;
        end

        function value = integrateInnerProduct(self, z, integrand, zBounds)
            % Integrate inner-product values in the native Chebyshev coordinate.
            %
            % The supplied `integrand` is a physical-coordinate integrand
            % sampled at `z`. The method integrates
            % $$f(z(x))\,q^{-1}(z(x))$$ in the native coordinate $$x$$,
            % where $$q=dx/dz$$.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: value = integrateInnerProduct(solver,z,integrand,zBounds)
            % - Parameter z: physical points from `innerProductGrid`
            % - Parameter integrand: physical-coordinate integrand values
            % - Parameter zBounds: physical integration bounds
            % - Returns value: definite integral over `zBounds`
            arguments
                self IMSolverSpectral
                z (:,1) double
                integrand (:,1) double
                zBounds (1,2) double
            end

            if length(z) ~= self.nEVP || max(abs(z(:) - self.zNative(:))) > self.gridTolerance()
                error("IMSolverSpectral:InvalidInnerProductGrid", ...
                    "Spectral inner products must be evaluated on the solver's native grid.");
            end
            if length(integrand) ~= self.nEVP
                error("IMSolverSpectral:InvalidIntegrandSize", ...
                    "The integrand must have one value per native grid point.");
            end

            q = self.qAtZ(self.zNative);
            integrandCheb = InternalModesSpectral.fct(integrand(:)./q(:));
            if self.boundsCoverDomain(zBounds)
                value = sum(self.chebyshevIntegrationWeights().*integrandCheb);
            else
                xBounds = self.xOfZ(sort(zBounds(:)));
                xBounds = min(max(xBounds, min(self.xNative)), max(self.xNative));
                value = InternalModesSpectral.IntegrateChebyshevVectorWithLimits( ...
                    integrandCheb, self.xNative, min(xBounds), max(xBounds));
            end
        end

        function x = xOfZ(self, z)
            % Map physical coordinate to native coordinate.
            %
            % - Topic: Evaluate native modes
            % - Declaration: x = xOfZ(solver,z)
            % - Parameter z: physical coordinate
            % - Returns x: native coordinate
            shape = size(z);
            z = z(:);
            x = NaN(size(z));
            inside = z >= self.zDomain(1) & z <= self.zDomain(2);
            if any(inside), x(inside) = feval(self.coordinateMap_,z(inside)); end
            x = reshape(x,shape);
        end

        function z = zOfX(self, x)
            % Invert the same smooth map used by physical differentiation.
            %
            % - Topic: Evaluate native modes
            % - Declaration: z = zOfX(solver,x)
            % - Parameter x: native coordinate
            % - Returns z: physical coordinate
            shape = size(x);
            x = x(:);
            z = NaN(size(x));
            inside = x >= self.xReference(1) & x <= self.xReference(end);
            if ~any(inside), z = reshape(z,shape); return; end
            target = x(inside);
            lower = repmat(self.zDomain(1),size(target));
            upper = repmat(self.zDomain(2),size(target));
            for iteration = 1:55
                middle = (lower+upper)/2;
                below = feval(self.coordinateMap_,middle) < target;
                lower(below) = middle(below);
                upper(~below) = middle(~below);
            end
            z(inside) = (lower+upper)/2;
            z(x == self.xReference(1)) = self.zDomain(1);
            z(x == self.xReference(end)) = self.zDomain(2);
            z = reshape(z,shape);
        end
    end

    methods (Hidden)
        function weights = innerProductWeights(self, z, zBounds)
            % Return physical-value weights for the spectral integral.
            arguments
                self IMSolverSpectral
                z (:,1) double
                zBounds (1,2) double
            end

            if length(z) ~= self.nEVP || max(abs(z(:)-self.zNative(:))) > self.gridTolerance()
                error("IMSolverSpectral:InvalidInnerProductGrid", ...
                    "Spectral inner products must be evaluated on the solver's native grid.");
            end
            if ~self.boundsCoverDomain(zBounds)
                weights = innerProductWeights@IMSolver(self,z,zBounds);
                return
            end

            nIntervals = self.nEVP-1;
            theta = pi*(0:nIntervals).'/nIntervals;
            interiorTheta = theta(2:nIntervals);
            valueWeights = zeros(self.nEVP,1);
            if mod(nIntervals,2) == 0
                valueWeights([1 end]) = 1/(nIntervals^2-1);
                harmonics = 1:(nIntervals/2-1);
                interior = 1-sum(2*cos(interiorTheta*(2*harmonics))./(4*harmonics.^2-1),2);
                interior = interior-cos(nIntervals*interiorTheta)/(nIntervals^2-1);
            else
                valueWeights([1 end]) = 1/nIntervals^2;
                harmonics = 1:((nIntervals-1)/2);
                interior = 1-sum(2*cos(interiorTheta*(2*harmonics))./(4*harmonics.^2-1),2);
            end
            valueWeights(2:nIntervals) = 2*interior/nIntervals;
            valueWeights = (max(self.xNative)-min(self.xNative))*valueWeights/2;
            weights = valueWeights./self.qAtZ(self.zNative);
        end
    end

    methods (Access = protected)
        function scale = eigenproblemColumnScale(~,nColumns)
            % Solve for degree-scaled Chebyshev coefficients. Suppressing
            % high-degree coefficient roundoff improves derivative accuracy
            % without removing columns or changing the physical pencil.
            scale = max(1,0:nColumns-1).^(-2);
        end

        function self = setupCoordinate(self)
            nReference = max(2001, 20*self.nEVP);
            self.zReference = linspace(self.zDomain(1), self.zDomain(2), nReference).';
            profileFunction = self.N2Function;
            switch self.coordinateKind
                case "z", derivativeFunction = @(z) ones(size(z));
                case "wkb", derivativeFunction = @(z) sqrt(profileFunction(z));
                case "density", derivativeFunction = profileFunction;
            end
            self.qReference = derivativeFunction(self.zReference);
            if ~isreal(self.qReference) || any(~isfinite(self.qReference)) || any(self.qReference <= 0)
                error("IMSolverSpectral:InvalidCoordinate", ...
                    "The native coordinate derivative dx/dz must be positive.");
            end
            % Derive x, dx/dz, and d2x/dz2 from one smooth representation.
            % Independent trapezoidal/PCHIP maps and endpoint finite
            % differences otherwise impose a floor on second derivatives.
            self.coordinateDerivative_ = chebfun(derivativeFunction,self.zDomain);
            if min(self.coordinateDerivative_) <= 0
                error("IMSolverSpectral:InvalidCoordinate","The represented coordinate derivative must remain positive throughout the domain.");
            end
            self.coordinateMap_ = cumsum(self.coordinateDerivative_);
            self.coordinateMap_ = self.coordinateMap_-feval(self.coordinateMap_,self.zDomain(1));
            self.coordinateSecondDerivative_ = diff(self.coordinateDerivative_);
            self.xReference = feval(self.coordinateMap_,self.zReference);
            self.qReference = feval(self.coordinateDerivative_,self.zReference);
            self.qzReference = feval(self.coordinateSecondDerivative_,self.zReference);
        end

        function self = setupNativeGrid(self)
            xMin = min(self.xReference);
            xMax = max(self.xReference);
            self.xNative = ((xMax - xMin)/2)*(cos(((0:self.nEVP-1).')*pi/(self.nEVP-1)) + 1) + xMin;
            self.zNative = self.zOfX(self.xNative);
            self.zNative(1) = self.zDomain(2);
            self.zNative(end) = self.zDomain(1);
            [self.T, self.Tx, self.Txx] = self.chebyshevPolynomialsAtNativePoints(self.xNative);
        end

        function q = qAtZ(self, z)
            q = feval(self.coordinateDerivative_,z);
        end

        function qz = qzAtZ(self, z)
            qz = feval(self.coordinateSecondDerivative_,z);
        end

        function x = clampNativeCoordinate(self, x)
            x = min(max(x, min(self.xNative)), max(self.xNative));
        end

        function value = isNativePhysicalGrid(self, z)
            value = length(z) == self.nEVP && max(abs(z(:) - self.zNative(:))) <= self.gridTolerance();
        end

        function [T, Tx, Txx] = chebyshevPolynomialsAtNativePoints(self, x, nPolynomials)
            x = x(:);
            if nargin < 3
                nPolynomials = self.nEVP;
            end
            xMin = min(self.xNative);
            xMax = max(self.xNative);
            L = xMax - xMin;
            xNorm = (2/L)*(x - xMin) - 1;
            xNorm = min(max(xNorm, -1), 1);
            theta = acos(xNorm);

            T = zeros(length(x),nPolynomials);
            for iPoly = 0:(nPolynomials-1)
                T(:,iPoly+1) = cos(iPoly*theta);
            end

            if nargout < 2
                return
            end

            Tx = self.differentiateChebyshevBasis(T, L);
            if nargout < 3
                return
            end

            Txx = self.differentiateChebyshevBasis(Tx, L);
        end

        function Tx = differentiateChebyshevBasis(~, T, L)
            nPolys = size(T,2);
            Tx = zeros(size(T));
            Tx(:,2) = T(:,1);
            Tx(:,3) = 4*T(:,2);
            for j = 4:nPolys
                m = j - 1;
                Tx(:,j) = (m/(m-2))*Tx(:,j-2) + 2*m*T(:,j-1);
            end
            Tx = (2/L)*Tx;
        end

        function weights = chebyshevIntegrationWeights(self)
            np = (0:(self.nEVP-1)).';
            weights = -(1+(-1).^np)./(np.*np - 1);
            weights(2) = 0;
            weights = (max(self.xNative) - min(self.xNative))*weights/2;
        end

        function value = boundsCoverDomain(self, zBounds)
            tolerance = self.gridTolerance();
            zBounds = sort(zBounds);
            value = abs(zBounds(1) - self.zDomain(1)) <= tolerance && ...
                abs(zBounds(2) - self.zDomain(2)) <= tolerance;
        end

        function tolerance = gridTolerance(self)
            tolerance = 100*eps(max(1,max(abs(self.zDomain))));
        end

    end
end
