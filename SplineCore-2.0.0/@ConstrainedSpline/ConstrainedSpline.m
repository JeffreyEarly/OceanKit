classdef ConstrainedSpline < TensorSpline
    % Tensor-product spline fit through noisy data values.
    %
    % `ConstrainedSpline` is the noisy-data fitting counterpart to
    % `InterpolatingSpline`. It fits a tensor-product spline basis to
    % observations sampled on a one-dimensional grid or a rectilinear tensor
    % grid, with optional robust weighting, correlated observation errors,
    % local point constraints, and global shape constraints.
    %
    % At each iteratively reweighted least-squares step it solves
    %
    % $$
    % \min_{\xi}\ (y - \mathbf{B}\xi)^{T} W (y - \mathbf{B}\xi)
    % $$
    %
    % subject to
    %
    % $$
    % A_{\mathrm{eq}}\xi = b_{\mathrm{eq}}, \qquad
    % A_{\mathrm{ineq}}\xi \le b_{\mathrm{ineq}}.
    % $$
    %
    % When the distribution model provides correlated errors, the code
    % forms an observation covariance
    %
    % $$
    % \Sigma_{ij} = \sigma_i \rho(x_i,x_j)\sigma_j
    % $$
    %
    % and applies the corresponding weighted solve through a matrix
    % factorization rather than explicitly forming $$\Sigma^{-1}$$.
    %
    % ## Basic usage
    %
    % Use `ConstrainedSpline` when you want to fit a tensor-product
    % spline to noisy values on a one-dimensional grid or rectilinear grid.
    %
    % ```matlab
    % x = linspace(0,1,20)';
    % y = exp(-20*(x-0.5).^2) + 0.05*randn(size(x));
    % spline = ConstrainedSpline(x, y, S=3, constraints=GlobalConstraint.positive());
    % yFit = spline(x);
    % ```
    %
    % - Topic: Create a constrained tensor spline
    % - Topic: Inspect fit results
    % - Topic: Analyze the fit
    % - Topic: Choose constraint locations
    % - Topic: Prepare knot sequences
    % - Topic: Prepare fit inputs
    % - Topic: Compile constraints
    % - Topic: Solve fit systems
    % - Declaration: classdef ConstrainedSpline < TensorSpline

    properties (SetAccess = private)
        % Grid vectors used to define the fitted rectilinear lattice.
        %
        % These are the one-dimensional coordinate vectors that define the
        % observation lattice passed to the constructor. They are the
        % grid-aligned counterpart to
        % [`dataPoints`](/spline-core/classes/constrainedspline/datapoints.html).
        %
        % - Topic: Inspect fit results
        gridVectors

        % Error model used while fitting the tensor spline.
        %
        % The `distribution` object defines how residuals are converted into
        % per-observation variances and, optionally, a correlation model. It
        % therefore controls the weight matrix in the objective
        %
        % $$
        % (y - \mathbf{B}\xi)^T W (y - \mathbf{B}\xi).
        % $$
        %
        % - Topic: Inspect fit results
        distribution

        % Observation locations as an N-by-D point matrix.
        %
        % Each row is one observation location in physical coordinates. For
        % gridded inputs, `dataPoints` is the explicit point-matrix form of
        % [`gridVectors`](/spline-core/classes/constrainedspline/gridvectors.html).
        %
        % - Topic: Inspect fit results
        dataPoints
        
        % Observation values as an N-by-1 vector.
        %
        % This is the flattened data vector `y` that appears in the
        % weighted least-squares objective.
        %
        % - Topic: Inspect fit results
        dataValues

        % Local point constraints used during fitting.
        %
        % These constraints are compiled into rows of
        % `Aeq * xi = beq` or `Aineq * xi <= bineq` by evaluating the spline
        % basis or its derivatives at specified points.
        %
        % - Topic: Inspect fit results
        pointConstraints

        % Global shape constraints used during fitting.
        %
        % These are coefficient-level sufficient conditions for shapes such
        % as positivity and monotonicity. They are compiled into the
        % inequality system `Aineq * xi <= bineq`.
        %
        % - Topic: Inspect fit results
        globalConstraints
    end

    properties (SetAccess = private, Hidden)
        % Inverse coefficient covariance or normal-equation system matrix.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        CmInv
        % Design matrix for the observation locations.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        X
        % Weight matrix or weights used by the fit.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        W
        % Linear equality constraints applied to the coefficient solve.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        Aeq
        % Right-hand side for equality constraints.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        beq
        % Linear inequality constraints applied to the coefficient solve.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        Aineq
        % Right-hand side for inequality constraints.
        %
        % - Topic: Inspect fit results
        % - Developer: true
        bineq
    end

    methods
        function self = ConstrainedSpline(grid, values, options)
            % Create a tensor-product spline fit to noisy observations.
            %
            % Use this constructor with a numeric vector in 1-D or a cell
            % array of grid vectors in higher dimensions when fitting noisy
            % tensor-product data sampled on a rectilinear grid.
            %
            % If the design matrix is $$\mathbf{B}$$, the coefficient vector
            % is estimated by an iteratively reweighted least-squares solve
            % with optional linear equality and inequality constraints. The
            % weights are updated from the current residuals through the
            % supplied error `distribution`.
            %
            % In one dimension, `S=N-1` together with `splineDOF=N` gives the
            % same least-squares polynomial fit as `polyfit(t,x,N-1)`.
            %
            % ```matlab
            % x = linspace(0,1,20)';
            % y = exp(-20*(x-0.5).^2) + 0.05*randn(size(x));
            % spline = ConstrainedSpline(x, y, S=3, constraints=GlobalConstraint.positive());
            % yFit = spline(x);
            % ```
            %
            % - Topic: Create a constrained tensor spline
            % - Declaration: self = ConstrainedSpline(grid,values,options)
            % - Parameter grid: numeric vector in 1-D or cell array of grid vectors in higher dimensions
            % - Parameter values: observation values
            % - Parameter options.S: optional spline degree scalar or vector with one entry per dimension
            % - Parameter options.knotPoints: optional knot vector in 1-D or cell array of knot vectors
            % - Parameter options.splineDOF: optional target number of splines per dimension
            % - Parameter options.distribution: optional error model object for the fit
            % - Parameter options.constraints: optional mixed SplineConstraint array
            % - Returns self: ConstrainedSpline instance
            arguments
                grid
                values {mustBeNumeric,mustBeReal,mustBeFinite}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
                options.knotPoints = []
                options.splineDOF {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = []
                options.distribution = []
                options.constraints SplineConstraint = SplineConstraint.empty(0,1)
            end

            if iscell(grid)
                if isempty(grid)
                    error('ConstrainedSpline:InvalidGrid', 'grid must not be empty.');
                end

                gridVectors = reshape(grid, 1, []);
                for iDim = 1:numel(gridVectors)
                    validateattributes(gridVectors{iDim}, {'numeric'}, {'vector','real','finite','nonempty'});
                    gridVectors{iDim} = reshape(gridVectors{iDim}, [], 1);
                end
            else
                validateattributes(grid, {'numeric'}, {'vector','real','finite','nonempty'});
                gridVectors = {reshape(grid, [], 1)};
            end

            numDimensions = numel(gridVectors);
            expectedSize = cellfun(@numel, gridVectors);
            if numDimensions == 1
                if ~(isvector(values) && numel(values) == expectedSize(1))
                    error('ConstrainedSpline:SizeMismatch', 'values must have size matching the lengths of the supplied grid inputs.');
                end
            else
                actualSize = size(values);
                if numel(actualSize) < numDimensions
                    actualSize = [actualSize, ones(1, numDimensions - numel(actualSize))];
                end

                if ~isequal(actualSize(1:numDimensions), expectedSize)
                    error('ConstrainedSpline:SizeMismatch', 'values must have size matching the lengths of the supplied grid inputs.');
                end
            end

            K = options.S + 1;
            K = TensorSpline.normalizeOrders(K, numDimensions);
            S = K - 1;

            if isempty(options.knotPoints)
                tKnot = [];
            elseif isnumeric(options.knotPoints)
                validateattributes(options.knotPoints, {'numeric'}, {'vector','real','finite'});
                if numDimensions ~= 1
                    error('ConstrainedSpline:InvalidKnotCell', 'knotPoints must be a knot vector in 1-D or a cell array with one knot vector per dimension.');
                end
                tKnot = {reshape(options.knotPoints, [], 1)};
            else
                tKnot = TensorSpline.normalizeKnotCell(options.knotPoints, numDimensions);
            end

            if isempty(options.splineDOF)
                splineDOF = [];
            elseif isscalar(options.splineDOF)
                splineDOF = repmat(options.splineDOF, 1, numDimensions);
            else
                splineDOF = reshape(options.splineDOF, 1, []);
                if numel(splineDOF) ~= numDimensions
                    error('ConstrainedSpline:InvalidDegreeOfFreedomOption', 'splineDOF must be scalar or have one element per dimension.');
                end
            end

            if isempty(tKnot)
                tKnot = cell(1, numDimensions);
                for iDim = 1:numDimensions
                    uniqueValues = unique(gridVectors{iDim}, 'sorted');
                    if numel(uniqueValues) < K(iDim)
                        tKnot{iDim} = [repmat(uniqueValues(1), K(iDim), 1); repmat(uniqueValues(end), K(iDim), 1)];
                    elseif isempty(splineDOF)
                        tKnot{iDim} = BSpline.knotPointsForDataPoints(gridVectors{iDim}, S=S(iDim));
                    else
                        tKnot{iDim} = BSpline.knotPointsForDataPoints(gridVectors{iDim}, S=S(iDim), splineDOF=splineDOF(iDim));
                    end
                end
            end
            distribution = options.distribution;
            if isempty(distribution)
                distribution = NormalDistribution(1);
            end

            pointMatrix = TensorSpline.pointsFromGridVectors(gridVectors);
            observedValues = values(:);

            for iDim = 1:numDimensions
                tKnot{iDim} = ConstrainedSpline.terminatedKnotPoints(tKnot{iDim}, S(iDim));
            end

            constraints = reshape(options.constraints, [], 1);
            [pointConstraints, globalConstraints] = ConstrainedSpline.normalizeConstraintInputs(constraints, numDimensions);

            Xbasis = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=tKnot, S=S);
            rho_X = [];
            if ~isempty(distribution.rho)
                delta = permute(pointMatrix, [1 3 2]) - permute(pointMatrix, [3 1 2]);
                rho_X = distribution.rho(sqrt(sum(delta.^2, 3)));
            end

            basisSize = reshape(cellfun(@numel, tKnot), 1, []) - reshape(K, 1, []);
            numCoefficients = prod(basisSize);
            Aeq = sparse([], [], [], 0, numCoefficients);
            beq = zeros(0,1);
            Aineq = sparse([], [], [], 0, numCoefficients);
            bineq = zeros(0,1);

            [pointAeq, pointBeq, pointAineq, pointBineq] = ConstrainedSpline.compilePointConstraints(pointConstraints, tKnot, K);
            [globalAineq, globalBineq] = ConstrainedSpline.compileGlobalConstraints(globalConstraints, tKnot, K);
            Aeq = [Aeq; pointAeq];
            beq = [beq; pointBeq];
            Aineq = [Aineq; pointAineq; globalAineq];
            bineq = [bineq; pointBineq; globalBineq];

            [coefficients,CmInv,W] = ConstrainedSpline.tensorModelSolution(  observedValues, Xbasis, distribution, rho_X, Aeq, beq, Aineq, bineq);

            self@TensorSpline(S=S, knotPoints=tKnot, xi=coefficients(:));
            self.coefficientsAreReadOnly = true;
            self.gridVectors = gridVectors;
            self.distribution = distribution;
            self.dataPoints = pointMatrix;
            self.dataValues = observedValues;
            self.pointConstraints = pointConstraints;
            self.globalConstraints = globalConstraints;
            self.CmInv = CmInv;
            self.X = Xbasis;
            self.W = W;
            self.Aeq = Aeq;
            self.beq = beq;
            self.Aineq = Aineq;
            self.bineq = bineq;
        end

        S = smoothingMatrix(self)
    end

    methods (Static)
        knotPoints = terminatedKnotPoints(knotPoints, S)
        tc = minimumConstraintPoints(knotPoints, S, T)
        [xi,CmInv,W] = tensorModelSolution(values, designMatrix, distribution, rho_X, Aeq, beq, Aineq, bineq)

        % Constraint compilation.
        [Aeq, beq, Aineq, bineq] = compilePointConstraints(pointConstraints, tKnot, K)
        [Aineq, bineq] = compileGlobalConstraints(globalConstraints, tKnot, K)
        A = monotonicDifferenceMatrix(basisSize, dim, direction)

        % Linear-system helpers.
        [normalMatrix, rhs] = weightedNormalEquations(X, x, W)
        [xi, systemMatrix] = constrainedWeightedSolution(normalMatrix, rhs, Aeq, beq, Aineq, bineq)
        W = weightMatrixFromSigma2(sigma2, rho_X)
        x = leftSolve(A, b)
        [pointConstraints, globalConstraints] = normalizeConstraintInputs(constraints, numDimensions)
    end
end
