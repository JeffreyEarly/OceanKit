classdef TensorSpline < CAAnnotatedClass
    % Tensor-product spline over multiple dimensions.
    %
    % `TensorSpline` is the multidimensional extension of `BSpline`. It
    % stores one spline degree and one knot vector per dimension, together
    % with a tensor-product coefficient array.
    %
    % The represented spline is
    %
    % $$
    % f(x_1,\ldots,x_d) = x_{\mathrm{Mean}} + x_{\mathrm{Std}}
    % \sum_{j_1=1}^{M_1} \cdots \sum_{j_d=1}^{M_d}
    % \xi_{j_1,\ldots,j_d}
    % \prod_{k=1}^{d} B_{j_k,S_k}(x_k;\tau_k),
    % $$
    %
    % where each $$\tau_k$$ is the knot vector for one coordinate
    % direction. Evaluation is pointwise on matching-size query arrays:
    % paired column vectors evaluate paired sample locations, while matching
    % `ndgrid` arrays evaluate a tensor-product lattice.
    %
    % ## Basic usage
    %
    % Use `TensorSpline.fromKnotPoints(...)` when you already have knot
    % vectors and tensor-product coefficients and want to evaluate the
    % resulting spline on matching-size query arrays. The low-level
    % `TensorSpline(...)` constructor is the cheap solved-state constructor
    % used directly by persistence and other internal bootstrap paths.
    %
    % ```matlab
    % knotPoints = {[0;0;0;0;1;1;1;1], [0;0;0;0;1;1;1;1]};
    % xi = randn(16,1);
    % spline = TensorSpline.fromKnotPoints(knotPoints, xi, S=[3 3]);
    %
    % [Xq,Yq] = ndgrid(linspace(0,1,40), linspace(0,1,40));
    % F = spline(Xq, Yq);
    % ```
    %
    % - Topic: Create a spline
    % - Topic: Inspect spline properties
    % - Topic: Evaluate the spline
    % - Topic: Transform the spline
    % - Topic: Build spline bases
    % - Declaration: classdef TensorSpline < CAAnnotatedClass

    properties (SetAccess = private)
        % Spline order in each tensor dimension.
        %
        % `K(k)` is the spline order along tensor dimension `k`. On each
        % fixed tensor cell, the spline is a polynomial of degree
        % `K(k)-1` in coordinate `x_k`, so `K=[4 4]` means a bicubic
        % tensor spline and `K=[2 3 4]` means linear, quadratic, and cubic
        % behavior across the three coordinates.
        %
        % The matching degree vector is
        % [`S`](/spline-core/classes/tensorspline/s.html), with
        % `S = K - 1`.
        %
        % ```matlab
        % spline = TensorSpline.fromKnotPoints(knotPoints, xi, S=[3 3]);
        % spline.K
        % % returns [4 4]
        % ```
        %
        % - Topic: Inspect spline properties
        K
    end

    properties (Access = private)
        % Internal tensor-product spline coefficients reshaped to basisSize.
        xi_
        % Internal knot vectors for each tensor dimension.
        tKnot_
    end

    properties (Access = protected, Hidden)
        coefficientsAreReadOnly (1,1) logical = false
    end

    properties (SetAccess = private)
        % Mean added back to zero-order evaluations.
        %
        % `xMean` is the additive term in the stored tensor-product model
        %
        % $$
        % f(x_1,\ldots,x_d) = x_{\mathrm{Mean}} + x_{\mathrm{Std}}
        % \sum_{j_1,\ldots,j_d} \xi_{j_1,\ldots,j_d}
        % \prod_{k=1}^{d} B_{j_k,S_k}(x_k;\tau_k).
        % $$
        %
        % It is mainly a numerical device: large offsets can be removed
        % before solving for `xi`, then added back only during zero-order
        % evaluation. As in the 1-D case, derivatives are unaffected by
        % `xMean`.
        %
        % ```matlab
        % spline = TensorSpline.fromKnotPoints(knotPoints, xi, S=[3 3], xMean=2.1, xStd=0.4);
        % values = spline(Xq, Yq);
        % ```
        %
        % - Topic: Inspect spline properties
        xMean (1,1) double {mustBeReal,mustBeFinite} = 0
        % Multiplicative scale applied to evaluations.
        %
        % `xStd` is the multiplicative scale factor in
        %
        % $$
        % f(x_1,\ldots,x_d) = x_{\mathrm{Mean}} + x_{\mathrm{Std}}
        % \sum_{j_1,\ldots,j_d} \xi_{j_1,\ldots,j_d}
        % \prod_{k=1}^{d} B_{j_k,S_k}(x_k;\tau_k).
        % $$
        %
        % It is useful when the fitted field has large or very small
        % amplitude: the stored coefficient array can remain close to order
        % one while evaluations and derivatives are rescaled back to
        % physical units.
        %
        % - Topic: Inspect spline properties
        xStd (1,1) double {mustBeReal,mustBeFinite} = 1
    end

    properties (Dependent)
        % Polynomial degree in each tensor dimension.
        %
        % `S(k)` is the polynomial degree along tensor dimension `k`. So
        % `S=[1 1]` gives bilinear pieces, `S=[3 3]` gives bicubic pieces,
        % and mixed values such as `S=[1 3]` are allowed when different
        % coordinates need different smoothness or complexity.
        %
        % The matching order vector is
        % [`K`](/spline-core/classes/tensorspline/k.html), with `K = S + 1`.
        %
        % - Topic: Inspect spline properties
        S
        % Number of tensor dimensions.
        %
        % This is the number of coordinate directions in the tensor-product
        % spline. In one dimension it is `1`; in higher dimensions it equals
        % the number of knot vectors in
        % [`knotPoints`](/spline-core/classes/tensorspline/knotpoints.html).
        %
        % ```matlab
        % spline.numDimensions
        % % returns 2 for a surface spline, 3 for a volume spline
        % ```
        %
        % - Topic: Inspect spline properties
        numDimensions
        % Number of basis functions in each dimension.
        %
        % If tensor dimension `k` has knot vector `tau_k` and order `K_k`,
        % then the number of one-dimensional basis functions in that
        % direction is
        %
        % $$
        % M_k = \mathrm{numel}(\tau_k) - K_k.
        % $$
        %
        % `basisSize` stores the row vector `[M_1 ... M_d]`. The total
        % coefficient count is `prod(basisSize)`.
        %
        % - Topic: Inspect spline properties
        basisSize
        % Tensor-product spline coefficients reshaped to basisSize.
        %
        % The coefficient array weights the tensor basis in
        %
        % $$
        % f(x_1,\ldots,x_d) = x_{\mathrm{Mean}} + x_{\mathrm{Std}}
        % \sum_{j_1,\ldots,j_d} \xi_{j_1,\ldots,j_d}
        % \prod_{k=1}^{d} B_{j_k,S_k}(x_k;\tau_k).
        % $$
        %
        % `xi` is stored reshaped to
        % [`basisSize`](/spline-core/classes/tensorspline/basissize.html),
        % so in 2-D it is a matrix, in 3-D it is an array, and so on.
        % The matrix-building helper is
        % [`matrixForPointMatrix`](/spline-core/classes/tensorspline/matrixforpointmatrix.html).
        %
        % - Topic: Inspect spline properties
        xi
        % Knot-axis objects defining the spline basis.
        %
        % `knotAxes` is the public axis-object representation of the
        % spline basis. Each axis stores one knot vector and is the
        % canonical constructor/persistence vocabulary for tensor splines.
        %
        % - Topic: Inspect spline properties
        knotAxes
        % Knot vectors defining the spline basis.
        %
        % Returns a numeric vector in 1-D and a cell array in higher dimensions.
        %
        % These are the per-dimension knot vectors
        % $$\tau_1, \ldots, \tau_d$$ that define the separable basis
        % functions $$B_{j_k,S_k}(x_k;\tau_k)$$.
        %
        % ```matlab
        % spline.knotPoints
        % % one vector in 1-D, one cell entry per dimension otherwise
        % ```
        %
        % - Topic: Inspect spline properties
        knotPoints
        % Minimum and maximum values of the spline domain in each dimension.
        %
        % `domain(k,:)` gives the lower and upper coordinate limits for
        % tensor dimension `k`:
        %
        % $$
        % \mathrm{domain}(k,:) = [\tau_{k,1},\ \tau_{k,\mathrm{end}}].
        % $$
        %
        % So a 2-D spline returns a `2 x 2` array of `[min max]` pairs, one
        % row for each coordinate direction.
        %
        % - Topic: Inspect spline properties
        domain
    end

    properties (Dependent, Hidden, SetAccess = private)
        xiPersisted
        coefficientIndex
        splineDimension
    end

    methods
        function self = TensorSpline(options)
            % Create a tensor-product spline from canonical solved state.
            %
            % Use this low-level constructor when you already have the
            % spline degree, knot-axis objects, and coefficient state. For
            % ordinary numeric knot-vector construction, use
            % `TensorSpline.fromKnotPoints(...)`.
            %
            % The stored spline has the tensor-product form
            %
            % $$
            % f(x_1,\ldots,x_d) = x_{\mathrm{Mean}} + x_{\mathrm{Std}}
            % \sum_{j_1,\ldots,j_d} \xi_{j_1,\ldots,j_d}
            % \prod_{k=1}^{d} B_{j_k,S_k}(x_k;\tau_k).
            % $$
            %
            % ```matlab
            % spline = TensorSpline(S=[3 3], knotAxes=SplineAxis.arrayFromVectors(knotPoints), xi=xi);
            % [Xq,Yq] = ndgrid(linspace(0,1,40), linspace(0,1,40));
            % values = spline(Xq, Yq);
            % ```
            %
            % - Topic: Create a spline
            % - Declaration: self = TensorSpline(options)
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Parameter options.knotAxes: ordered SplineAxis array defining the tensor-product basis
            % - Parameter options.xi: tensor-product coefficient array or vector
            % - Parameter options.xMean: optional additive output offset
            % - Parameter options.xStd: optional multiplicative output scale
            % - Returns self: TensorSpline instance
            arguments
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
                options.knotAxes (:,1) SplineAxis {mustBeNonempty}
                options.xi {mustBeNumeric,mustBeReal,mustBeFinite}
                options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
                options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
            end

            knotAxes = reshape(options.knotAxes, [], 1);
            numDimensions = numel(knotAxes);
            tKnot = SplineAxis.vectorsFromArray(knotAxes);
            K = TensorSpline.normalizeOrders(options.S + 1, numDimensions);

            self@CAAnnotatedClass();
            self.K = K;
            self.tKnot_ = tKnot;
            self.xi = options.xi;
            self.xMean = options.xMean;
            self.xStd = options.xStd;
        end

        function value = get.S(self)
            % Return the spline polynomial degree in each dimension.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.S(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: row vector equal to K - 1
            value = self.K - 1;
        end

        function value = get.numDimensions(self)
            % Return the number of tensor dimensions.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.numDimensions(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: number of dimensions
            value = numel(self.K);
        end

        function value = get.basisSize(self)
            % Return the number of basis functions per dimension.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.basisSize(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: row vector of basis sizes
            value = TensorSpline.basisSizeFromKnotCell(self.tKnot_, self.K);
        end

        function value = get.xi(self)
            % Return the tensor-product spline coefficients.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.xi(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: coefficient array reshaped to basisSize
            value = self.xi_;
        end

        function set.xi(self, value)
            % Update the tensor-product coefficients with canonical reshaping.
            %
            % - Topic: Inspect spline properties
            % - Declaration: set.xi(self,value)
            % - Parameter self: TensorSpline instance
            % - Parameter value: coefficient vector or array with prod(basisSize) entries
            arguments
                self (1,1) TensorSpline
                value {mustBeNumeric,mustBeReal,mustBeFinite}
            end

            if self.coefficientsAreReadOnly
                error('ConstrainedSpline:ReadOnlyCoefficients',  'Constrained spline coefficients are read-only after fitting.');
            end

            basisSize = self.basisSize;
            if numel(value) ~= prod(basisSize)
                error('TensorSpline:InvalidCoefficientCount',  'xi must contain exactly prod(basisSize) coefficients.');
            end

            if isscalar(basisSize)
                self.xi_ = reshape(value, basisSize, 1);
            else
                self.xi_ = reshape(value, basisSize);
            end
        end

        function value = get.knotPoints(self)
            % Return the knot vectors defining the tensor-product basis.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.knotPoints(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: knot vector in 1-D, cell array otherwise
            if self.numDimensions == 1
                value = self.tKnot_{1};
            else
                value = self.tKnot_;
            end
        end

        function value = get.knotAxes(self)
            % Return the knot-axis objects defining the tensor-product basis.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.knotAxes(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: ordered SplineAxis array
            value = SplineAxis.arrayFromVectors(self.tKnot_);
        end

        function value = get.xiPersisted(self)
            % Return the flattened coefficient vector used for persistence.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.xiPersisted(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: coefficient column vector
            value = reshape(self.xi_, [], 1);
        end

        function value = get.coefficientIndex(self)
            % Return the coefficient-index coordinate for persisted coefficients.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.coefficientIndex(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: coefficient index vector
            value = reshape(1:numel(self.xi_), [], 1);
        end

        function value = get.splineDimension(self)
            % Return the dimension-index coordinate used for vector-valued persisted state.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.splineDimension(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: dimension index vector
            value = reshape(1:self.numDimensions, [], 1);
        end

        function value = get.domain(self)
            % Return the domain limits for each dimension.
            %
            % - Topic: Inspect spline properties
            % - Declaration: value = get.domain(self)
            % - Parameter self: TensorSpline instance
            % - Returns value: numDimensions-by-2 array of [min max] domain limits
            value = cell2mat(cellfun(@(tk) [tk(1), tk(end)], self.tKnot_, 'UniformOutput', false)');
        end
        varargout = subsref(self, index)
        values = valueAtPoints(self, X, options)

    end

    methods (Static)
        function self = fromKnotPoints(knotPoints, xi, options)
            % Create a tensor-product spline from numeric knot vectors and coefficients.
            %
            % Use this factory for ordinary scientific construction when
            % you have numeric knot vectors or a knot-vector cell array.
            %
            % - Topic: Create a spline
            % - Declaration: self = fromKnotPoints(knotPoints,xi,options)
            % - Parameter knotPoints: numeric knot vector in 1-D or cell array of knot vectors
            % - Parameter xi: tensor-product coefficient array or vector
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Parameter options.xMean: optional additive output offset
            % - Parameter options.xStd: optional multiplicative output scale
            % - Returns self: TensorSpline instance
            arguments
                knotPoints {mustBeNonempty}
                xi {mustBeNumeric,mustBeReal,mustBeFinite}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
                options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
                options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
            end

            if isnumeric(knotPoints)
                validateattributes(knotPoints, {'numeric'}, {'vector','real','finite','nonempty'});
                tKnot = {reshape(knotPoints, [], 1)};
            else
                tKnot = TensorSpline.normalizeKnotCell(knotPoints, numel(knotPoints));
            end

            self = TensorSpline(S=options.S, knotAxes=SplineAxis.arrayFromVectors(tKnot), xi=xi, xMean=options.xMean, xStd=options.xStd);
        end

        B = matrixForPointMatrix(pointMatrix, options)
        [pointMatrix, gridSize] = pointsFromGridVectors(gridVectors)
        [pointMatrix, supportVectors] = pointsOfSupportFromKnotPoints(knotPoints, options)
    end

    methods (Static, Hidden)
        function self = annotatedClassFromFile(path)
            ncfile = NetCDFFile(path);
            if ~isKey(ncfile.attributes, 'AnnotatedClass')
                error('TensorSpline:MissingAnnotatedClass', 'Unable to find the AnnotatedClass attribute in %s.', path);
            end

            className = string(ncfile.attributes('AnnotatedClass'));
            if ncfile.hasGroupWithName(className)
                group = ncfile.groupWithName(className);
            else
                group = ncfile;
            end
            self = TensorSpline.annotatedClassFromGroup(group);
        end

        function self = annotatedClassFromGroup(group)
            className = string(group.attributes('AnnotatedClass'));
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, feval(strcat(className, '.classRequiredPropertyNames')));
            vars.xi = vars.xiPersisted;
            vars = rmfield(vars, 'xiPersisted');
            varCell = namedargs2cell(vars);
            self = feval(className, varCell{:});
            self.restoreOptionalPersistedPropertiesFromGroup(group);
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('splineDimension', '', 'Index over tensor dimensions.');
            propertyAnnotations(end+1) = CANumericProperty('S', {'splineDimension'}, '', 'Spline degree vector $$S$$.');
            propertyAnnotations(end+1) = CAObjectProperty('knotAxes', 'Ordered knot-axis objects.');
            propertyAnnotations(end+1) = CADimensionProperty('coefficientIndex', '', 'Index over persisted spline coefficients.');
            propertyAnnotations(end+1) = CANumericProperty('xiPersisted', {'coefficientIndex'}, '', 'Flattened tensor-product spline coefficients $$\\xi$$.');
            propertyAnnotations(end+1) = CANumericProperty('xMean', {}, '', 'Additive output offset $$x_{\\mathrm{Mean}}$$.');
            propertyAnnotations(end+1) = CANumericProperty('xStd', {}, '', 'Multiplicative output scale $$x_{\\mathrm{Std}}$$.');
        end

        function names = classRequiredPropertyNames()
            names = {'S', 'knotAxes', 'xiPersisted', 'xMean', 'xStd'};
        end

        K = normalizeOrders(K, numDimensions)
        tKnot = normalizeKnotCell(tKnot, numDimensions)
        basisSize = basisSizeFromKnotCell(tKnot, K)
        derivativeOrders = normalizeDerivativeOrders(derivativeOrders, numDimensions)
        [xi, tKnot, K] = differentiateAlongDimension(xi, tKnot, K, derivativeOrder, dim)
        [xi, tKnot, K] = integrateAlongDimension(xi, tKnot, K, dim, xMean, xStd)
        spline = zeroSplineForDomain(domain, numDimensions, options)
    end

end
