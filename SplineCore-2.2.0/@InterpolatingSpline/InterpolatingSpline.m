classdef InterpolatingSpline < TensorSpline
    % Interpolating spline on one-dimensional samples or rectilinear grids.
    %
    % `InterpolatingSpline` is the exact-fit counterpart to
    % `ConstrainedSpline`. Use
    % `InterpolatingSpline.fromGriddedValues(...)` when your samples live
    % on a one-dimensional grid or a rectilinear tensor grid and should be
    % reproduced exactly. The low-level `InterpolatingSpline(...)`
    % constructor is the cheap solved-state constructor used for persisted
    % restart and other direct bootstrap paths.
    %
    % Supported construction forms:
    %   spline = InterpolatingSpline.fromGriddedValues(x, V)
    %   spline = InterpolatingSpline.fromGriddedValues({x, y, ...}, V)
    %   spline = InterpolatingSpline(S=S, knotAxes=..., xi=..., gridAxes=...)
    %
    % If $$\mathbf{B}$$ is the tensor-product basis matrix on the supplied
    % grid and $$\tilde{y}$$ is the normalized data vector, the stored
    % coefficients are chosen so that
    %
    % $$
    % \mathbf{B}\xi = \tilde{y}, \qquad
    % \tilde{y} = \frac{y - \bar{y}}{s_y},
    % $$
    %
    % where $$xMean = \bar{y}$$ and $$xStd = s_y$$ are stored so later
    % evaluation returns values on the original scale.
    %
    % ## Basic usage
    %
    % ```matlab
    % x = linspace(0,1,8)';
    % y = linspace(-1,1,9)';
    % [X,Y] = ndgrid(x, y);
    % F = sin(2*pi*X).*cos(pi*Y);
    % spline = InterpolatingSpline.fromGriddedValues({x, y}, F, S=[3 3]);
    % Fq = spline(X, Y);
    % ```
    %
    % - Topic: Create an interpolating spline
    % - Topic: Inspect interpolation grids
    % - Declaration: classdef InterpolatingSpline < TensorSpline

    properties (SetAccess = private)
        % Grid vectors used to define the interpolation lattice.
        %
        % These are the original 1-D sample locations in each coordinate
        % direction. In one dimension `gridVectors` contains one column
        % vector; in higher dimensions it stores one grid vector per tensor
        % axis.
        %
        % If `[X1,...,Xd] = ndgrid(gridVectors{:})`, then the spline
        % interpolates the supplied value array exactly on that lattice.
        %
        % - Topic: Inspect interpolation grids
        gridVectors
    end

    properties (Dependent, SetAccess = private)
        % Grid-axis objects used to define the interpolation lattice.
        %
        % `gridAxes` is the axis-object representation of the original
        % interpolation grid and is the canonical constructor vocabulary
        % for persisted or otherwise pre-solved interpolation state.
        %
        % - Topic: Inspect interpolation grids
        gridAxes
    end

    methods
        function self = InterpolatingSpline(options)
            % Create an interpolating spline from canonical solved state.
            %
            % Use this low-level constructor when you already have the
            % interpolation knot axes, coefficient state, and grid axes.
            % For ordinary interpolation from gridded samples, use
            % `InterpolatingSpline.fromGriddedValues(...)`.
            %
            % ```matlab
            % spline = InterpolatingSpline( ...
            %     S=[3 3], ...
            %     knotAxes=SplineAxis.arrayFromVectors(knotPoints), ...
            %     xi=xi, ...
            %     gridAxes=SplineAxis.arrayFromVectors({x, y}));
            % ```
            %
            % - Topic: Create an interpolating spline
            % - Declaration: self = InterpolatingSpline(options)
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Parameter options.knotAxes: ordered knot-axis objects defining the spline basis
            % - Parameter options.xi: interpolating coefficient vector or array
            % - Parameter options.gridAxes: ordered interpolation-grid axis objects
            % - Parameter options.xMean: optional additive output offset
            % - Parameter options.xStd: optional multiplicative output scale
            % - Returns self: InterpolatingSpline instance
            arguments
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
                options.knotAxes (:,1) SplineAxis {mustBeNonempty}
                options.xi {mustBeNumeric,mustBeReal,mustBeFinite}
                options.gridAxes (:,1) SplineAxis {mustBeNonempty}
                options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
                options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
            end

            if numel(options.knotAxes) ~= numel(options.gridAxes)
                error('InterpolatingSpline:AxisCountMismatch', 'knotAxes and gridAxes must have the same number of dimensions.');
            end

            self@TensorSpline(S=options.S, knotAxes=options.knotAxes, xi=options.xi, xMean=options.xMean, xStd=options.xStd);
            gridVectors = SplineAxis.vectorsFromArray(options.gridAxes);
            if numel(gridVectors) == 1
                self.gridVectors = gridVectors{1};
            else
                self.gridVectors = gridVectors;
            end
        end

        function value = get.gridAxes(self)
            % Return the grid-axis objects defining the interpolation lattice.
            %
            % - Topic: Inspect interpolation grids
            % - Declaration: value = get.gridAxes(self)
            % - Parameter self: InterpolatingSpline instance
            % - Returns value: ordered SplineAxis array
            value = SplineAxis.arrayFromVectors(self.gridVectors);
        end
    end

    methods (Static)
        function self = fromGriddedValues(gridVectors, values, options)
            % Create an interpolating spline from values on a rectilinear grid.
            %
            % Use this factory for ordinary interpolation from gridded
            % samples. It validates the gridded input, builds the knot
            % sequence, normalizes the values, solves for the interpolation
            % coefficients, and then delegates to the cheap constructor.
            %
            % - Topic: Create an interpolating spline
            % - Declaration: self = fromGriddedValues(gridVectors,values,options)
            % - Parameter gridVectors: numeric vector in 1-D or cell array of grid vectors in higher dimensions
            % - Parameter values: array of sampled values on the supplied grid
            % - Parameter options.S: spline degree scalar or vector with one entry per dimension
            % - Returns self: InterpolatingSpline instance
            arguments
                gridVectors {mustBeNonempty}
                values {mustBeNumeric,mustBeReal,mustBeFinite}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
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
                    error('InterpolatingSpline:SizeMismatch', 'values must have size matching the lengths of the supplied grid inputs.');
                end
            else
                actualSize = size(values);
                if numel(actualSize) < numDimensions
                    actualSize = [actualSize, ones(1, numDimensions - numel(actualSize))];
                end

                if ~isequal(actualSize(1:numDimensions), expectedSize)
                    error('InterpolatingSpline:SizeMismatch', 'values must have size matching the lengths of the supplied grid inputs.');
                end
            end

            K = TensorSpline.normalizeOrders(options.S + 1, numDimensions);
            S = K - 1;
            tKnot = cell(1, numDimensions);
            for iDim = 1:numDimensions
                tKnot{iDim} = BSpline.knotPointsForDataPoints(gridVectors{iDim}, S=S(iDim));
            end

            xMean = mean(values(:));
            values = values - xMean;
            xStd = std(values(:));
            if xStd > 0
                values = values / xStd;
            else
                xStd = 1;
            end

            [gridPoints, ~] = TensorSpline.pointsFromGridVectors(gridVectors);
            X = TensorSpline.matrixForPointMatrix(gridPoints, knotPoints=tKnot, S=S);
            xi = X \ values(:);
            self = InterpolatingSpline(S=S, knotAxes=SplineAxis.arrayFromVectors(tKnot), xi=xi, xMean=xMean, xStd=xStd, gridAxes=SplineAxis.arrayFromVectors(gridVectors));
        end
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = TensorSpline.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CAObjectProperty('gridAxes', 'Ordered interpolation-grid axes.');
        end

        function names = classRequiredPropertyNames()
            names = union(TensorSpline.classRequiredPropertyNames(), {'gridAxes'});
        end
    end
end
