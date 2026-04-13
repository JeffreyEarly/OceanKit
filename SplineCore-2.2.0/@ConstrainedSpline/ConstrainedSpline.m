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
    % Use `ConstrainedSpline.fromData(...)` for ordinary one-dimensional
    % noisy-data fitting and `ConstrainedSpline.fromGriddedValues(...)`
    % when the observations lie on a rectilinear tensor grid.
    %
    % ```matlab
    % t = linspace(0,1,20)';
    % x = exp(-20*(t-0.5).^2) + 0.05*randn(size(t));
    % spline = ConstrainedSpline.fromData(t, x, S=3, constraints=GlobalConstraint.positive());
    % xFit = spline(t);
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
    % - Topic: Persist spline state
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

    properties (Access = private)
        CmInv_
        X_
        W_
        Aeq_
        beq_
        Aineq_
        bineq_
    end

    properties (SetAccess = private, Hidden, Dependent)
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

    properties (Dependent, SetAccess = private)
        % Grid-axis objects used to define the fitted rectilinear lattice.
        %
        % `gridAxes` is the axis-object representation of the original fit
        % grid and is the canonical constructor vocabulary for persisted or
        % otherwise pre-solved constrained splines.
        %
        % - Topic: Inspect fit results
        gridAxes
    end

    properties (Dependent, Hidden, SetAccess = private)
        dataPointIndex
        dataDimension
    end

    methods
        function self = ConstrainedSpline(options)
            % Create a constrained spline from canonical solved state.
            %
            % Use this low-level constructor when you already have the
            % solved spline coefficients, fit grid, observations, and
            % semantic constraints. For ordinary one-dimensional fitting,
            % use `ConstrainedSpline.fromData(...)`. For rectilinear-grid
            % fitting, use `ConstrainedSpline.fromGriddedValues(...)`.
            %
            % - Topic: Create a constrained tensor spline
            % - Declaration: self = ConstrainedSpline(options)
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Parameter options.knotAxes: ordered knot-axis objects defining the spline basis
            % - Parameter options.xi: fitted coefficient vector or array
            % - Parameter options.gridAxes: ordered fit-grid axis objects
            % - Parameter options.distribution: error model used during the fit
            % - Parameter options.dataPoints: observation locations as an N-by-D point matrix
            % - Parameter options.dataValues: observation values as an N-by-1 vector
            % - Parameter options.pointConstraints: optional PointConstraint array used during fitting
            % - Parameter options.globalConstraints: optional GlobalConstraint array used during fitting
            % - Parameter options.xMean: optional additive output offset
            % - Parameter options.xStd: optional multiplicative output scale
            % - Returns self: ConstrainedSpline instance
            % - nav_order: 1
            arguments
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
                options.knotAxes (:,1) SplineAxis {mustBeNonempty}
                options.xi {mustBeNumeric,mustBeReal,mustBeFinite}
                options.gridAxes (:,1) SplineAxis {mustBeNonempty}
                options.distribution (1,1) Distribution
                options.dataPoints (:,:) double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.dataValues (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.pointConstraints (:,1) PointConstraint = PointConstraint.empty(0,1)
                options.globalConstraints (:,1) GlobalConstraint = GlobalConstraint.empty(0,1)
                options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
                options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
            end

            numDimensions = numel(options.knotAxes);
            if numel(options.gridAxes) ~= numDimensions
                error('ConstrainedSpline:AxisCountMismatch', 'knotAxes and gridAxes must have the same number of dimensions.');
            end

            if size(options.dataPoints, 2) ~= numDimensions
                error('ConstrainedSpline:InvalidDataPointDimension', 'dataPoints must have one column per spline dimension.');
            end

            if size(options.dataPoints, 1) ~= numel(options.dataValues)
                error('ConstrainedSpline:InvalidObservationCount', 'dataValues must have one entry per row of dataPoints.');
            end

            for iConstraint = 1:numel(options.pointConstraints)
                if options.pointConstraints(iConstraint).numDimensions ~= numDimensions
                    error('ConstrainedSpline:InvalidPointConstraintDimension', 'Each PointConstraint must match the spline dimension.');
                end
            end

            for iConstraint = 1:numel(options.globalConstraints)
                if ~isempty(options.globalConstraints(iConstraint).dimension) && options.globalConstraints(iConstraint).dimension > numDimensions
                    error('ConstrainedSpline:InvalidGlobalConstraintDimension', 'Each GlobalConstraint dimension must not exceed the spline dimension.');
                end
            end

            self@TensorSpline(S=options.S, knotAxes=options.knotAxes, xi=options.xi, xMean=options.xMean, xStd=options.xStd);
            self.coefficientsAreReadOnly = true;
            self.initializePersistedFitState(options.gridAxes, options.distribution, options.dataPoints, options.dataValues, options.pointConstraints, options.globalConstraints);
            self.initializeFitDiagnostics();
        end

        function value = get.gridAxes(self)
            % Return the grid-axis objects used to define the fitted lattice.
            %
            % - Topic: Inspect fit results
            % - Declaration: value = get.gridAxes(self)
            % - Parameter self: ConstrainedSpline instance
            % - Returns value: ordered SplineAxis array
            value = SplineAxis.arrayFromVectors(self.gridVectors);
        end

        function value = get.dataPointIndex(self)
            % Return the observation-index coordinate for persisted fit data.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.dataPointIndex(self)
            % - Parameter self: ConstrainedSpline instance
            % - Returns value: observation index vector
            value = reshape(1:size(self.dataPoints, 1), [], 1);
        end

        function value = get.dataDimension(self)
            % Return the coordinate-dimension index for persisted point matrices.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.dataDimension(self)
            % - Parameter self: ConstrainedSpline instance
            % - Returns value: spatial-dimension index vector
            value = reshape(1:size(self.dataPoints, 2), [], 1);
        end

        function value = get.CmInv(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.CmInv_;
        end

        function value = get.X(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.X_;
        end

        function value = get.W(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.W_;
        end

        function value = get.Aeq(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.Aeq_;
        end

        function value = get.beq(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.beq_;
        end

        function value = get.Aineq(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.Aineq_;
        end

        function value = get.bineq(self)
            self.ensureFitDiagnosticsInitialized();
            value = self.bineq_;
        end

        S = smoothingMatrix(self)
    end

    methods (Access = private)
        function initializePersistedFitState(self,gridAxes,distribution,dataPoints,dataValues,pointConstraints,globalConstraints)
            arguments
                self (1,1) ConstrainedSpline
                gridAxes (:,1) SplineAxis
                distribution (1,1) Distribution
                dataPoints (:,:) double {mustBeReal,mustBeFinite}
                dataValues (:,1) double {mustBeReal,mustBeFinite}
                pointConstraints (:,1) PointConstraint = PointConstraint.empty(0,1)
                globalConstraints (:,1) GlobalConstraint = GlobalConstraint.empty(0,1)
            end

            gridVectors = SplineAxis.vectorsFromArray(gridAxes);
            if numel(gridVectors) == 1
                self.gridVectors = gridVectors{1};
            else
                self.gridVectors = gridVectors;
            end
            self.distribution = distribution;
            self.dataPoints = dataPoints;
            self.dataValues = dataValues;
            self.pointConstraints = pointConstraints;
            self.globalConstraints = globalConstraints;
        end

        function initializeFitDiagnostics(self,CmInv,X,W,Aeq,beq,Aineq,bineq)
            arguments
                self (1,1) ConstrainedSpline
                CmInv = []
                X = []
                W = []
                Aeq = []
                beq = []
                Aineq = []
                bineq = []
            end

            self.CmInv_ = CmInv;
            self.X_ = X;
            self.W_ = W;
            self.Aeq_ = Aeq;
            self.beq_ = beq;
            self.Aineq_ = Aineq;
            self.bineq_ = bineq;
        end

        function ensureFitDiagnosticsInitialized(self)
            if ~isempty(self.CmInv_)
                return;
            end

            tKnot = SplineAxis.vectorsFromArray(self.knotAxes);
            Xbasis = TensorSpline.matrixForPointMatrix(self.dataPoints, knotPoints=tKnot, S=self.S);
            [Aeq, beq, Aineq, bineq] = ConstrainedSpline.rebuildConstraintSystems(self.pointConstraints, self.globalConstraints, tKnot, self.K);
            rho_X = ConstrainedSpline.correlationMatrixForDataPoints(self.dataPoints, self.distribution);
            residual = self.dataValues - Xbasis*reshape(self.xi, [], 1);
            sigma2 = self.distribution.w(residual);
            W = ConstrainedSpline.weightMatrixFromSigma2(sigma2, rho_X);
            [normalMatrix, ~] = ConstrainedSpline.weightedNormalEquations(Xbasis, self.dataValues, W);
            CmInv = ConstrainedSpline.diagnosticSystemMatrix(normalMatrix, Aeq, Aineq);
            self.initializeFitDiagnostics(CmInv, Xbasis, W, Aeq, beq, Aineq, bineq);
        end
    end

    methods (Static)
        function self = fromData(t, x, options)
            % Create a constrained spline fit from one-dimensional samples.
            %
            % Use this factory for ordinary one-dimensional noisy-data
            % fitting. It validates the shared sample vectors, then
            % delegates to the general rectilinear-grid factory so the
            % one-dimensional behavior stays aligned with
            % `ConstrainedSpline.fromGriddedValues(...)`.
            %
            % - Topic: Create a constrained tensor spline
            % - Declaration: self = fromData(t,x,options)
            % - Parameter t: one-dimensional sample locations
            % - Parameter x: observation values sampled at `t`
            % - Parameter options.S: spline degree
            % - Parameter options.knotPoints: optional one-dimensional knot vector
            % - Parameter options.splineDOF: optional target number of splines
            % - Parameter options.distribution: optional error model object for the fit
            % - Parameter options.constraints: optional mixed SplineConstraint array
            % - Returns self: ConstrainedSpline instance
            % - nav_order: 2
            arguments
                t (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                x (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
                options.knotPoints = []
                options.splineDOF {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = []
                options.distribution (1,1) Distribution = NormalDistribution(sigma=1)
                options.constraints SplineConstraint = SplineConstraint.empty(0,1)
            end

            if numel(t) ~= numel(x)
                error('ConstrainedSpline:SizeMismatch', 'x must have the same number of elements as t.');
            end

            optionCell = namedargs2cell(options);
            self = ConstrainedSpline.fromGriddedValues(t, x, optionCell{:});
        end

        function self = fromGriddedValues(gridVectors, values, options)
            % Create a constrained spline fit from values on a rectilinear grid.
            %
            % Use this factory for fitting from rectilinear-grid samples,
            % especially in two or more dimensions. In one dimension,
            % prefer `ConstrainedSpline.fromData(...)`. This method
            % validates the gridded input, chooses the knot sequence,
            % normalizes the constraint objects, solves the fit system, and
            % then delegates to the cheap solved-state constructor.
            %
            % - Topic: Create a constrained tensor spline
            % - Declaration: self = fromGriddedValues(gridVectors,values,options)
            % - Parameter gridVectors: numeric vector in 1-D or cell array of grid vectors in higher dimensions
            % - Parameter values: array of observation values on the supplied grid
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Parameter options.knotPoints: optional knot vector in 1-D or cell array of knot vectors
            % - Parameter options.splineDOF: optional target number of splines per dimension
            % - Parameter options.distribution: optional error model object for the fit
            % - Parameter options.constraints: optional mixed SplineConstraint array
            % - Returns self: ConstrainedSpline instance
            % - nav_order: 3
            arguments
                gridVectors {mustBeNonempty}
                values {mustBeNumeric,mustBeReal,mustBeFinite}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
                options.knotPoints = []
                options.splineDOF {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = []
                options.distribution (1,1) Distribution = NormalDistribution(sigma=1)
                options.constraints SplineConstraint = SplineConstraint.empty(0,1)
            end

            if iscell(gridVectors)
                gridVectors = reshape(gridVectors, 1, []);
                for iDim = 1:numel(gridVectors)
                    validateattributes(gridVectors{iDim}, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
                    gridVectors{iDim} = reshape(gridVectors{iDim}, [], 1);
                end
            else
                validateattributes(gridVectors, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
                gridVectors = {reshape(gridVectors, [], 1)};
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

            K = TensorSpline.normalizeOrders(options.S + 1, numDimensions);
            S = K - 1;
            if isempty(options.knotPoints)
                tKnot = [];
            elseif isnumeric(options.knotPoints)
                validateattributes(options.knotPoints, {'numeric'}, {'vector', 'real', 'finite'});
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

            pointMatrix = TensorSpline.pointsFromGridVectors(gridVectors);
            observedValues = values(:);
            for iDim = 1:numDimensions
                tKnot{iDim} = ConstrainedSpline.terminatedKnotPoints(tKnot{iDim}, S(iDim));
            end

            constraints = reshape(options.constraints, [], 1);
            [pointConstraints, globalConstraints] = ConstrainedSpline.normalizeConstraintInputs(constraints, numDimensions);
            Xbasis = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=tKnot, S=S);
            rho_X = ConstrainedSpline.correlationMatrixForDataPoints(pointMatrix, options.distribution);

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

            [coefficients, CmInv, W] = ConstrainedSpline.tensorModelSolution(observedValues, Xbasis, options.distribution, rho_X, Aeq, beq, Aineq, bineq);
            self = ConstrainedSpline(S=S, knotAxes=SplineAxis.arrayFromVectors(tKnot), xi=coefficients(:), gridAxes=SplineAxis.arrayFromVectors(gridVectors), distribution=options.distribution, dataPoints=pointMatrix, dataValues=observedValues, pointConstraints=pointConstraints, globalConstraints=globalConstraints, xMean=0, xStd=1);
            self.initializeFitDiagnostics(CmInv, Xbasis, W, Aeq, beq, Aineq, bineq);
        end

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

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = TensorSpline.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CAObjectProperty('gridAxes', 'Ordered grid-axis objects used by the fit.');
            propertyAnnotations(end+1) = CAObjectProperty('distribution', 'Distribution used while fitting the spline.');
            propertyAnnotations(end+1) = CADimensionProperty('dataPointIndex', '', 'Observation index for persisted fit data.');
            propertyAnnotations(end+1) = CADimensionProperty('dataDimension', '', 'Coordinate dimension for persisted point matrices.');
            propertyAnnotations(end+1) = CANumericProperty('dataPoints', {'dataPointIndex', 'dataDimension'}, '', 'Observation locations as an N-by-D point matrix.');
            propertyAnnotations(end+1) = CANumericProperty('dataValues', {'dataPointIndex'}, '', 'Observation values as an N-by-1 vector.');
            propertyAnnotations(end+1) = CAObjectProperty('pointConstraints', 'PointConstraint array used while fitting the spline.');
            propertyAnnotations(end+1) = CAObjectProperty('globalConstraints', 'GlobalConstraint array used while fitting the spline.');
        end

        function names = classRequiredPropertyNames()
            names = union(TensorSpline.classRequiredPropertyNames(), {'gridAxes', 'distribution', 'dataPoints', 'dataValues', 'pointConstraints', 'globalConstraints'});
        end
    end

    methods (Static, Access = private)
        function rho_X = correlationMatrixForDataPoints(dataPoints,distribution)
            arguments
                dataPoints (:,:) double {mustBeReal,mustBeFinite}
                distribution (1,1) Distribution
            end

            rho_X = [];
            if isempty(distribution.rho)
                return;
            end

            delta = permute(dataPoints, [1 3 2]) - permute(dataPoints, [3 1 2]);
            rho_X = distribution.rho(sqrt(sum(delta.^2, 3)));
        end

        function [Aeq, beq, Aineq, bineq] = rebuildConstraintSystems(pointConstraints,globalConstraints,tKnot,K)
            arguments
                pointConstraints (:,1) PointConstraint
                globalConstraints (:,1) GlobalConstraint
                tKnot (1,:) cell
                K
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
        end

        function systemMatrix = diagnosticSystemMatrix(normalMatrix,Aeq,Aineq)
            arguments
                normalMatrix
                Aeq
                Aineq
            end

            if isempty(Aineq)
                if isempty(Aeq)
                    systemMatrix = normalMatrix;
                else
                    systemMatrix = [normalMatrix, Aeq'; Aeq, zeros(size(Aeq,1))];
                end
            else
                systemMatrix = normalMatrix;
            end
        end
    end
end
