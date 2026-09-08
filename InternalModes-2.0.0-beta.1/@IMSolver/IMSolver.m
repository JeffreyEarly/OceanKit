classdef (Abstract) IMSolver
    % Define the shared protocol for canonical EVP solvers.
    %
    % Concrete solvers own the grid, coordinate mapping, derivative
    % matrices, integration rule, and interpolation of native modes. The
    % base class owns the common generalized-eigenvalue workflow.
    %
    % - Topic: Solve EVPs
    % - Topic: Solve geostrophic zero-APV modes
    % - Topic: Developer topics
    % - Declaration: classdef (Abstract) IMSolver

    methods
        function [z,weights] = nativeQuadratureRule(self,zBounds)
            % Return the solver's native physical quadrature rule.
            %
            % The returned points are increasing in physical coordinate and
            % `weights` are reordered with them. For a WKB-configured
            % `IMSolverSpectral`, the points are Chebyshev--Lobatto points in
            % $$x(z)=\int N(z)\,dz$$ and the weights act directly on values
            % sampled in physical $$z$$.
            %
            % - Topic: Solve EVPs
            % - Declaration: [z,weights] = nativeQuadratureRule(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: increasing physical quadrature points
            % - Returns weights: physical-coordinate quadrature weights
            arguments
                self IMSolver
                zBounds (1,2) double
            end

            nativeZ = self.innerProductGrid(zBounds);
            nativeWeights = self.innerProductWeights(nativeZ,zBounds);
            [z,index] = sort(nativeZ(:));
            weights = nativeWeights(index);
        end

        function [z,weights,Dz] = nativeDifferentiationRule(self,zBounds,derivativeOrder)
            % Return one ordered grid for physical quadrature and differentiation.
            %
            % The returned points increase in physical coordinate. `weights`
            % integrate values in that order, and `Dz*values` differentiates
            % columns sampled in the same order. This keeps the grid, weights,
            % and differentiation matrix as one inseparable numerical rule.
            %
            % - Topic: Solve EVPs
            % - Declaration: [z,weights,Dz] = nativeDifferentiationRule(solver,zBounds,derivativeOrder)
            % - Parameter zBounds: physical integration bounds
            % - Parameter derivativeOrder: physical derivative order; default `1`
            % - Returns z: increasing physical grid points
            % - Returns weights: physical-coordinate quadrature weights
            % - Returns Dz: matrix mapping values on `z` to physical derivatives on `z`
            arguments
                self IMSolver
                zBounds (1,2) double
                derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 1
            end

            nativeZ = self.zNative(:);
            quadratureZ = self.innerProductGrid(zBounds);
            if length(quadratureZ) ~= length(nativeZ) || max(abs(sort(quadratureZ(:))-sort(nativeZ))) > 100*eps(max(1,max(abs(nativeZ))))
                error("IMSolver:UnsupportedDifferentiationBounds", ...
                    "nativeDifferentiationRule requires bounds whose quadrature uses the complete native solver grid.");
            end
            nativeWeights = self.innerProductWeights(nativeZ,zBounds);
            nativeDerivative = self.differentiateGridValues(eye(length(nativeZ)),derivativeOrder);
            [z,index] = sort(nativeZ(:));
            weights = nativeWeights(index);
            Dz = nativeDerivative(index,index);
        end

        function basisSet = solveEVP(self, evp, options)
            % Solve an EVP and return a basis set.
            %
            % If the assembled matrices produce no finite real eigenvalues,
            % `solveEVP` throws a matrix-level diagnostic before returning.
            % Before selection, the EVP may reject family-specific
            % numerical representations that are not finite physical modes.
            %
            % - Topic: Solve EVPs
            % - Declaration: basisSet = solveEVP(solver,evp,options)
            % - Parameter evp: canonical EVP descriptor
            % - Parameter options.nModes: number of modes to retain
            % - Returns basisSet: solved basis set
            arguments
                self IMSolver
                evp IMEigenvalueProblem
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 100
            end

            solver = self.configuredForEVP(evp);
            [A, B] = evp.assemble(solver);
            % Equilibrate equations before QZ: spectral derivative rows can
            % otherwise overwhelm boundary equations as resolution grows.
            % Left scaling preserves the generalized eigenproblem and its
            % native eigenvectors. Keep the physical matrices below for
            % zero-mode classification and matrix-level diagnostics.
            rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
            rowScale(rowScale == 0) = 1;
            solveB = B./rowScale;
            [V, D] = eig(A./rowScale, solveB);
            eigenvalues = diag(D);
            valid = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues)) ...
                & abs(imag(eigenvalues)) < 1e-8*max(1,abs(real(eigenvalues)));
            if ~any(valid)
                error("IMSolver:NoValidEigenvalues", "%s", ...
                    self.noValidEigenvalueMessage(evp, A, B, eigenvalues, valid));
            end

            V = real(V(:,valid));
            eigenvalues = real(eigenvalues(valid));
            % Numerical null-space rejection must use the metric in the
            % same equation scaling as the eigensolve.
            familyValid = evp.finiteGeneralizedEigenpairMask(V,solveB);
            V = V(:,familyValid);
            eigenvalues = eigenvalues(familyValid);
            if isempty(eigenvalues)
                error("IMSolver:NoValidEigenvalues", ...
                    "EVP ""%s"" produced no finite real eigenmode with nonzero metric action; every candidate lies in the numerical null space of B.",evp.name);
            end
            selection = evp.selectModes(eigenvalues(:), options.nModes, solver, A);
            eigenvalues = eigenvalues(selection.sortIndex);
            eigenvalues(selection.modeNumber == 0) = 0;
            V = V(:,selection.sortIndex);
            basisSet = evp.makeBasisSet(solver, V, eigenvalues(:).', ...
                selection.modeNumber, selection.modeSelectionDiagnostics);
            basisSet = basisSet.orientModeSigns();
        end

        function zRoots = rootsOfNativeMode(self, nativeMode)
            % Return physical roots of one native mode.
            %
            % Concrete solvers may implement this developer hook when their
            % native representation supports accurate root finding. The base
            % implementation reports that mode-root grids are unavailable.
            %
            % - Topic: Developer topics
            % - Declaration: zRoots = rootsOfNativeMode(solver,nativeMode)
            % - Parameter nativeMode: one native mode column
            % - Returns zRoots: physical roots in the solver domain
            % - Developer: true
            arguments
                self IMSolver
                nativeMode (:,1) double {mustBeReal, mustBeFinite}
            end

            zRoots = zeros(0,1); %#ok<NASGU>
            error("IMSolver:UnsupportedModeRoots", ...
                "%s does not implement rootsOfNativeMode for %d native values.", class(self), length(nativeMode));
        end

        function basisSet = solveGeostrophicZeroAPVModes(self, problem)
            % Solve canonical geostrophic zero-APV boundary modes.
            %
            % The operator is factored once for each distinct requested
            % horizontal wavenumber. Every requested unit endpoint response is solved
            % in the same multiple-right-hand-side operation. No
            % generalized-energy coefficient or rotation enters this solve.
            %
            % - Topic: Solve geostrophic zero-APV modes
            % - Declaration: basisSet = solveGeostrophicZeroAPVModes(solver,problem)
            % - Parameter problem: geostrophic zero-APV problem
            % - Returns basisSet: canonical boundary-normalized basis
            arguments
                self IMSolver
                problem IMGeostrophicZeroAPVModes
            end

            solver = self.configuredForGeostrophicZeroAPVModes(problem);
            z = solver.zNative;
            n = length(z);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            D2 = solver.physicalDerivativeMatrix(2);
            N2Values = problem.N2(z);
            N2Values = N2Values(:);
            if length(N2Values) ~= n
                error("IMGeostrophicZeroAPVModes:InvalidStratification", "N2 must return one value for each solver grid point.");
            end
            if any(~isfinite(N2Values)) || any(N2Values <= 0)
                error("IMGeostrophicZeroAPVModes:InvalidStratification", "N2 must be finite and positive on the solver grid.");
            end

            pValues = problem.f0^2 ./ N2Values;
            pzValues = solver.differentiateGridValues(pValues, 1);
            baseMatrix = diag(pValues)*D2 + diag(pzValues)*D1;
            surfaceIndex = solver.boundaryIndex("surface");
            bottomIndex = solver.boundaryIndex("bottom");
            N2Surface = N2Values(surfaceIndex);
            N2Bottom = N2Values(bottomIndex);

            surfaceRow = D1(surfaceIndex,:);
            if problem.surfaceBoundary == "freeSurface"
                surfaceRow = surfaceRow + (N2Surface/problem.g)*D0(surfaceIndex,:);
            end
            bottomRow = D1(bottomIndex,:);

            nEndpoints = numel(problem.endpoints);
            rightHandSides = zeros(n,nEndpoints);
            for iEndpoint = 1:nEndpoints
                switch problem.endpoints(iEndpoint)
                    case "surface"
                        rightHandSides(surfaceIndex,iEndpoint) = -N2Surface/problem.g;
                    case "bottom"
                        rightHandSides(bottomIndex,iEndpoint) = -N2Bottom/problem.g;
                end
            end

            nK = numel(problem.k);
            nativeModes = zeros(n,nEndpoints,nK);
            solvedWavenumbers = zeros(1,nK);
            solvedModes = cell(1,nK);
            nSolved = 0;
            for iK = 1:nK
                solvedIndex = find(solvedWavenumbers(1:nSolved) == problem.k(iK),1);
                if isempty(solvedIndex)
                    matrix = baseMatrix - problem.k(iK)^2*D0;
                    matrix(surfaceIndex,:) = surfaceRow;
                    matrix(bottomIndex,:) = bottomRow;
                    nSolved = nSolved + 1;
                    solvedWavenumbers(nSolved) = problem.k(iK);
                    solvedModes{nSolved} = solver.solveBoundaryValueSystems(matrix,rightHandSides);
                    solvedIndex = nSolved;
                end
                nativeModes(:,:,iK) = solvedModes{solvedIndex};
            end

            metadata = problem.metadata;
            metadata.solutionKind = "geostrophicZeroAPVModes";
            metadata.k = problem.k;
            metadata.endpoints = problem.endpoints;
            metadata.surfaceBoundary = problem.surfaceBoundary;
            metadata.modesPerWavenumber = nEndpoints;
            metadata.factorizations = nSolved;
            basisSet = IMGeostrophicZeroAPVModesBasis(problem=problem, solver=solver, nativeModes=nativeModes, metadata=metadata);
        end
    end

    methods (Hidden)
        function weights = innerProductWeights(self, z, zBounds)
            % Return the linear weights used by integrateInnerProduct.
            %
            % This generic implementation obtains the exact weights of the
            % configured solver's integration operator.  It lets callers
            % integrate many columns with one matrix multiplication while
            % preserving each concrete solver's existing quadrature rule.
            arguments
                self IMSolver
                z (:,1) double
                zBounds (1,2) double
            end

            nGrid = length(z);
            weights = zeros(nGrid,1);
            impulse = zeros(nGrid,1);
            for iGrid = 1:nGrid
                impulse(iGrid) = 1;
                weights(iGrid) = self.integrateInnerProduct(z,impulse,zBounds);
                impulse(iGrid) = 0;
            end
        end
    end

    methods (Access = protected)
        function values = solveBoundaryValueSystems(~, matrix, rightHandSides)
            % Solve one boundary-value matrix for multiple response columns.
            matrixFactorization = decomposition(matrix);
            values = matrixFactorization \ rightHandSides;
        end
    end

    methods (Access = private)
        function message = noValidEigenvalueMessage(~, evp, A, B, eigenvalues, valid)
            normA = norm(A, "fro");
            normB = norm(B, "fro");
            rankA = rank(A);
            rankB = rank(B);
            finite = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues));
            threshold = 1e-12*max([normA, normB, 1]);
            hint = "";
            if normA <= threshold || normB <= threshold
                hint = " One assembled matrix is nearly zero; check EVP coefficient parameters and physical constants such as f0.";
            end
            message = sprintf('EVP "%s" produced no finite real eigenvalues. total eigenvalues=%d, finite eigenvalues=%d, valid-real eigenvalues=%d, norm(A,"fro")=%.3g, norm(B,"fro")=%.3g, rank(A)=%d, rank(B)=%d.%s', ...
                evp.name, length(eigenvalues), nnz(finite), nnz(valid), normA, normB, rankA, rankB, hint);
        end
    end

    methods (Abstract)
        solver = configuredForEVP(self, evp)

        % Return a solver configured for geostrophic zero-APV modes.
        %
        % Concrete solvers prepare their native grid, coordinate mapping,
        % and derivative matrices for the supplied zero-APV problem.
        %
        % - Topic: Solve geostrophic zero-APV modes
        % - Topic: Developer topics
        % - Declaration: solver = configuredForGeostrophicZeroAPVModes(solver,problem)
        % - Parameter problem: geostrophic zero-APV problem
        % - Returns solver: configured solver
        % - Developer: true
        solver = configuredForGeostrophicZeroAPVModes(self, problem)

        context = context(self)
        values = N2(self, z)
        D = physicalDerivativeMatrix(self, derivativeOrder)
        values = differentiateGridValues(self, values, derivativeOrder)
        index = boundaryIndex(self, location)
        values = evaluateNativeModes(self, nativeModes, z)
        values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
        nativeIntegral = integrateGridValuesFromSurface(self, gridValues)
        z = innerProductGrid(self, zBounds)
        value = integrateInnerProduct(self, z, integrand, zBounds)
    end
end
