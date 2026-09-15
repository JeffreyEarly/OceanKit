classdef IMSolverFiniteDifference < IMSolver
    % Solve physical-coordinate EVPs on a supplied finite-difference grid.
    %
    % `IMSolverFiniteDifference` uses the user's physical `z`
    % grid as its native basis. It shares the `solveEVP` interface with
    % the spectral solver but validates the supplied grid against the EVP
    % domain before solving.
    %
    % ```matlab
    % solver = IMSolverFiniteDifference(z=linspace(-1000,0,65).');
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Inspect solvers
    % - Topic: Solve EVPs
    % - Topic: Solve geostrophic zero-APV modes
    % - Topic: Assemble EVPs
    % - Topic: Evaluate native modes
    % - Topic: Developer topics
    % - Declaration: classdef IMSolverFiniteDifference

    properties (SetAccess = private)
        % Number of native EVP values.
        %
        % - Topic: Inspect solvers
        nEVP

        % Physical vertical domain.
        %
        % - Topic: Inspect solvers
        zDomain

        % Native finite-difference grid.
        %
        % - Topic: Inspect solvers
        zNative

        % Native coordinate alias for API consistency.
        %
        % - Topic: Inspect solvers
        xNative

        % Native value matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        T

        % Native first-derivative matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Tx

        % Native second-derivative matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Txx
    end

    methods
        function self = IMSolverFiniteDifference(options)
            % Create a finite-difference solver from a physical grid.
            %
            % - Topic: Create solvers
            % - Declaration: solver = IMSolverFiniteDifference(options)
            % - Parameter options.z: finite-difference grid
            % - Returns solver: initialized finite-difference solver
            arguments
                options.z (:,1) double
            end

            self.zNative = sort(options.z(:), "descend");
            self.xNative = self.zNative;
            self.zDomain = [min(self.zNative) max(self.zNative)];
            self.nEVP = length(self.zNative);
            self.T = eye(self.nEVP);
            self.Tx = self.finiteDifferenceMatrix(1);
            self.Txx = self.finiteDifferenceMatrix(2);
        end

        function solver = configuredForEVP(self, evp)
            % Return a finite-difference solver configured for an EVP.
            %
            % - Topic: Solve EVPs
            % - Topic: Developer topics
            % - Declaration: solver = configuredForEVP(solver,evp)
            % - Parameter evp: canonical EVP descriptor
            % - Returns solver: validated solver
            % - Developer: true
            arguments
                self IMSolverFiniteDifference
                evp IMEigenvalueProblem
            end

            solver = self;
            expectedDomain = sort(evp.zDomain);
            tolerance = 100*eps(max([1 abs(expectedDomain) abs(self.zDomain)]));
            if max(abs(self.zDomain - expectedDomain)) > tolerance
                error("IMSolverFiniteDifference:DomainMismatch", ...
                    "The finite-difference grid domain must match evp.zDomain.");
            end
        end

        function solver = configuredForGeostrophicZeroAPVModes(self, problem)
            % Return a finite-difference solver configured for zero-APV modes.
            %
            % - Topic: Solve geostrophic zero-APV modes
            % - Topic: Developer topics
            % - Declaration: solver = configuredForGeostrophicZeroAPVModes(solver,problem)
            % - Parameter problem: geostrophic zero-APV problem
            % - Returns solver: validated solver
            % - Developer: true
            arguments
                self IMSolverFiniteDifference
                problem IMGeostrophicZeroAPVModes
            end

            solver = self;
            expectedDomain = sort(problem.zDomain);
            tolerance = 100*eps(max([1 abs(expectedDomain) abs(self.zDomain)]));
            if max(abs(self.zDomain - expectedDomain)) > tolerance
                error("IMSolverFiniteDifference:DomainMismatch", "The finite-difference grid domain must match problem.zDomain.");
            end
        end

        function context = context(self)
            % Return the framework coefficient context.
            %
            % Solvers provide discretization fields. EVPs add the physical
            % domain, medium, and constants.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: framework coefficient context
            context.coordinateKind = "finiteDifference";
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = N2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            error("IMSolverFiniteDifference:UnsupportedOperation", ...
                "Finite-difference solvers do not own N2; use the EVP or basis set instead.");
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
                error("IMSolverFiniteDifference:InvalidGridValues", ...
                    "Grid values must have one row per native grid point.");
            end
            values = self.physicalDerivativeMatrix(derivativeOrder)*values;
        end

        function D = physicalDerivativeMatrix(self, derivativeOrder)
            % Return a native finite-difference matrix for a physical derivative.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: D = physicalDerivativeMatrix(solver,derivativeOrder)
            % - Parameter derivativeOrder: physical derivative order
            % - Returns D: derivative matrix
            switch derivativeOrder
                case 0
                    D = self.T;
                case 1
                    D = self.Tx;
                case 2
                    D = self.Txx;
                otherwise
                    error("IMSolverFiniteDifference:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function index = boundaryIndex(self, location)
            % Return the native row for a physical boundary.
            %
            % - Topic: Solve EVPs
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
                    error("IMSolverFiniteDifference:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"" or ""bottom"".");
            end
        end

        function values = evaluateNativeModes(self, nativeModes, z)
            % Interpolate native finite-difference modes to requested points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluateNativeModes(solver,nativeModes,z)
            % - Parameter nativeModes: native mode columns
            % - Parameter z: physical evaluation points
            % - Returns values: interpolated mode values
            values = interp1(flip(self.zNative), flip(nativeModes,1), z(:), "pchip");
        end

        function values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
            % Evaluate finite-difference physical derivatives at requested points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluatePhysicalDerivative(solver,nativeModes,z,derivativeOrder)
            % - Parameter nativeModes: native mode columns
            % - Parameter z: physical evaluation points
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: derivative values
            derivativeValues = self.physicalDerivativeMatrix(derivativeOrder)*nativeModes;
            values = interp1(flip(self.zNative), flip(derivativeValues,1), z(:), "pchip");
        end

        function nativeIntegral = integrateGridValuesFromSurface(self, gridValues)
            % Integrate physical grid values downward from the surface.
            %
            % For physical integrands $$f(z)$$ sampled on the descending
            % finite-difference grid, the returned native value columns are
            %
            % $$
            % I(z)=\int_z^{z_s}f(z')\,dz',
            % \qquad I(z_s)=0.
            % $$
            %
            % - Topic: Developer topics
            % - Declaration: nativeIntegral = integrateGridValuesFromSurface(solver,gridValues)
            % - Parameter gridValues: physical integrand columns sampled on `zNative`
            % - Returns nativeIntegral: native finite-difference values for the surface-referenced integrals
            % - Developer: true
            arguments
                self IMSolverFiniteDifference
                gridValues (:,:) double {mustBeReal, mustBeFinite}
            end

            if size(gridValues,1) ~= self.nEVP
                error("IMSolverFiniteDifference:InvalidGridValues", "Grid values must have one row per native grid point.");
            end
            nativeIntegral = -cumtrapz(self.zNative,gridValues,1);
            nativeIntegral(1,:) = 0;
        end

        function z = innerProductGrid(self, zBounds)
            % Return a bounded physical grid for finite-difference inner products.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: z = innerProductGrid(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: ascending physical integration grid
            arguments
                self IMSolverFiniteDifference
                zBounds (1,2) double
            end

            zBounds = sort(zBounds);
            zAscending = sort(self.zNative);
            interior = zAscending > zBounds(1) & zAscending < zBounds(2);
            z = unique([zBounds(1); zAscending(interior); zBounds(2)]);
        end

        function value = integrateInnerProduct(self, z, integrand, zBounds)
            % Integrate finite-difference inner-product values with trapezoids.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: value = integrateInnerProduct(solver,z,integrand,zBounds)
            % - Parameter z: physical integration grid
            % - Parameter integrand: physical-coordinate integrand values
            % - Parameter zBounds: physical integration bounds
            % - Returns value: definite integral over `zBounds`
            arguments
                self IMSolverFiniteDifference
                z (:,1) double
                integrand (:,1) double
                zBounds (1,2) double
            end

            [zSorted, sortIndex] = sort(z(:));
            value = trapz(zSorted, integrand(sortIndex));
        end

        function x = xOfZ(~, z)
            % Return the native coordinate for a physical finite-difference grid.
            %
            % - Topic: Evaluate native modes
            % - Declaration: x = xOfZ(solver,z)
            % - Parameter z: physical coordinate
            % - Returns x: native coordinate
            x = z;
        end

        function z = zOfX(~, x)
            % Return the physical coordinate for a native finite-difference grid.
            %
            % - Topic: Evaluate native modes
            % - Declaration: z = zOfX(solver,x)
            % - Parameter x: native coordinate
            % - Returns z: physical coordinate
            z = x;
        end
    end

    methods (Hidden)
        function weights = innerProductWeights(~, z, ~)
            % Return physical-value weights for trapezoidal integration.
            [zSorted,sortIndex] = sort(z(:));
            spacing = diff(zSorted);
            sortedWeights = [spacing(1)/2;(spacing(1:end-1)+spacing(2:end))/2;spacing(end)/2];
            weights = zeros(size(z));
            weights(sortIndex) = sortedWeights;
        end
    end

    methods (Access = private)
        function D = finiteDifferenceMatrix(self, derivativeOrder)
            if derivativeOrder == 0
                D = eye(self.nEVP);
                return;
            end

            nStencil = min(self.nEVP, max(derivativeOrder+1, 5));
            D = zeros(self.nEVP,self.nEVP);
            for iRow = 1:self.nEVP
                firstIndex = iRow - floor(nStencil/2);
                firstIndex = max(1, min(firstIndex, self.nEVP - nStencil + 1));
                stencilIndex = firstIndex:(firstIndex + nStencil - 1);
                weights = self.finiteDifferenceWeights(self.zNative(iRow), self.zNative(stencilIndex), derivativeOrder);
                D(iRow,stencilIndex) = weights;
            end
        end

        function weights = finiteDifferenceWeights(~, z0, z, derivativeOrder)
            powers = 0:(length(z)-1);
            A = (z(:).' - z0).^powers(:);
            b = zeros(length(z),1);
            b(derivativeOrder+1) = factorial(derivativeOrder);
            weights = (A\b).';
        end

    end
end
