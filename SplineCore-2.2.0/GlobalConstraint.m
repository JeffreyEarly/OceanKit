classdef GlobalConstraint < SplineConstraint
    % Specify a global shape constraint for a constrained spline fit.
    %
    % Use `GlobalConstraint` to describe semantic whole-domain constraints
    % such as positivity or monotonicity. These objects are intended to be
    % compiled into linear coefficient inequalities by constrained fitting
    % classes.
    %
    % ## Basic usage
    %
    % ```matlab
    % c1 = GlobalConstraint.positive();
    % c2 = GlobalConstraint.monotonicIncreasing(dimension=2);
    % ```
    %
    % - Topic: Specify global constraints
    % - Declaration: classdef GlobalConstraint

    properties (Constant, Hidden)
        positiveShape = "positive"
        monotonicIncreasingShape = "monotonicIncreasing"
        monotonicDecreasingShape = "monotonicDecreasing"
    end

    properties (SetAccess = private)
        % Shape-constraint kind.
        %
        % - Topic: Inspect global constraint properties
        shape (1,1) string = ""

        % Tensor dimension associated with the constraint, when applicable.
        %
        % - Topic: Inspect global constraint properties
        dimension = double.empty(1,0)
    end

    properties (Dependent, Hidden, SetAccess = private)
        dimensionPersisted
    end

    methods
        function self = GlobalConstraint(varargin)
            % Create a global spline-constraint object.
            %
            % Use the static helper methods `positive`,
            % `monotonicIncreasing`, and `monotonicDecreasing` for the
            % intended public construction style.
            %
            % - Topic: Specify global constraints
            % - Declaration: self = GlobalConstraint(options)
            % - Parameter options.shape: global constraint shape identifier
            % - Parameter options.dimension: optional tensor dimension for monotonic constraints
            % - Returns self: GlobalConstraint instance
            self@SplineConstraint();
            if nargin == 0
                return;
            end

            options = GlobalConstraint.parseInputs(varargin{:});
            self.shape = options.shape;
            if ~isnan(options.dimensionPersisted)
                self.dimension = options.dimensionPersisted;
            elseif ~isnan(options.dimension)
                self.dimension = options.dimension;
            else
                self.dimension = double.empty(1,0);
            end
        end

        function value = get.dimensionPersisted(self)
            % Return the persisted scalar view of `dimension`.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.dimensionPersisted(self)
            % - Parameter self: GlobalConstraint instance
            % - Returns value: scalar dimension index or `NaN` when no dimension is associated with the constraint
            if isempty(self.dimension)
                value = NaN;
            else
                value = self.dimension;
            end
        end
    end

    methods (Static)
        function self = positive()
            % Create a positivity constraint.
            %
            % ```matlab
            % c = GlobalConstraint.positive();
            % ```
            %
            % - Topic: Specify global constraints
            % - Declaration: self = positive()
            % - Returns self: positivity GlobalConstraint
            self = GlobalConstraint(shape=GlobalConstraint.positiveShape);
        end

        function self = monotonicIncreasing(options)
            % Create a monotone-increasing constraint along one dimension.
            %
            % ```matlab
            % c = GlobalConstraint.monotonicIncreasing(dimension=1);
            % ```
            %
            % - Topic: Specify global constraints
            % - Declaration: self = monotonicIncreasing(options)
            % - Parameter options.dimension: tensor dimension, default 1
            % - Returns self: monotone-increasing GlobalConstraint
            arguments
                options.dimension (1,1) double {mustBeInteger,mustBePositive} = 1
            end

            self = GlobalConstraint(shape=GlobalConstraint.monotonicIncreasingShape, dimension=options.dimension);
        end

        function self = monotonicDecreasing(options)
            % Create a monotone-decreasing constraint along one dimension.
            %
            % ```matlab
            % c = GlobalConstraint.monotonicDecreasing(dimension=2);
            % ```
            %
            % - Topic: Specify global constraints
            % - Declaration: self = monotonicDecreasing(options)
            % - Parameter options.dimension: tensor dimension, default 1
            % - Returns self: monotone-decreasing GlobalConstraint
            arguments
                options.dimension (1,1) double {mustBeInteger,mustBePositive} = 1
            end

            self = GlobalConstraint(shape=GlobalConstraint.monotonicDecreasingShape, dimension=options.dimension);
        end
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = SplineConstraint.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CAPropertyAnnotation('shape', 'Global shape-constraint kind.');
            propertyAnnotations(end+1) = CANumericProperty('dimensionPersisted', {}, '', 'Persisted tensor dimension associated with the constraint.');
        end

        function names = classRequiredPropertyNames()
            names = {'shape', 'dimensionPersisted'};
        end
    end

    methods (Static, Access = private)
        function options = parseInputs(options)
            arguments
                options.shape (1,1) string
                options.dimension (1,1) double {mustBeReal} = NaN
                options.dimensionPersisted (1,1) double {mustBeReal} = NaN
            end

            mustBeMember(options.shape, [GlobalConstraint.positiveShape, GlobalConstraint.monotonicIncreasingShape, GlobalConstraint.monotonicDecreasingShape]);

            if ~isnan(options.dimension) && ~isnan(options.dimensionPersisted)
                error('GlobalConstraint:AmbiguousDimension', 'Specify either dimension or dimensionPersisted, not both.');
            end

            if ~isnan(options.dimension) && (options.dimension <= 0 || round(options.dimension) ~= options.dimension)
                error('GlobalConstraint:InvalidDimension', 'dimension must be a positive integer when it is supplied.');
            end

            if ~isnan(options.dimensionPersisted) && (options.dimensionPersisted <= 0 || round(options.dimensionPersisted) ~= options.dimensionPersisted)
                error('GlobalConstraint:InvalidPersistedDimension', 'dimensionPersisted must be a positive integer or NaN.');
            end
        end
    end
end
