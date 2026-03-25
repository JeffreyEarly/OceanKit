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

    methods (Access = private)
        function self = GlobalConstraint
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
            self = GlobalConstraint;
            self.shape = GlobalConstraint.positiveShape;
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

            self = GlobalConstraint;
            self.shape = GlobalConstraint.monotonicIncreasingShape;
            self.dimension = options.dimension;
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

            self = GlobalConstraint;
            self.shape = GlobalConstraint.monotonicDecreasingShape;
            self.dimension = options.dimension;
        end
    end
end
