classdef SplineAxis < CAAnnotatedClass
    % Store one ordered numeric axis for spline persistence.
    %
    % `SplineAxis` is the internal persisted representation of one knot or
    % grid axis used by the spline classes. It keeps the NetCDF schema
    % object-based while preserving the public cell-array vector APIs on
    % `TensorSpline`, `InterpolatingSpline`, and `ConstrainedSpline`.
    %
    % ```matlab
    % axis = SplineAxis(values=linspace(0, 1, 8)');
    % ```
    %
    % - Topic: Create a spline
    % - Topic: Inspect spline properties
    % - Topic: Persist spline state
    % - Declaration: classdef SplineAxis < CAAnnotatedClass

    properties (SetAccess = private)
        % Ordered numeric axis values.
        %
        % `values` stores the one-dimensional coordinate vector for one
        % spline grid or knot axis.
        %
        % - Topic: Inspect spline properties
        values (:,1) double {mustBeReal,mustBeFinite} = double.empty(0,1)
    end

    methods
        function self = SplineAxis(options)
            % Create a spline axis from a coordinate vector.
            %
            % - Topic: Create a spline
            % - Declaration: self = SplineAxis(options)
            % - Parameter options.values: ordered numeric axis vector
            % - Returns self: SplineAxis instance
            arguments
                options.values (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
            end

            self@CAAnnotatedClass();
            self.values = options.values;
        end
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CADimensionProperty('values', '', 'Ordered spline-axis values.');
        end

        function names = classRequiredPropertyNames()
            names = {'values'};
        end

        function axes = arrayFromVectors(vectors)
            arguments
                vectors
            end

            if isa(vectors, 'SplineAxis')
                axes = reshape(vectors, [], 1);
                return;
            end

            if isnumeric(vectors)
                vectors = {reshape(vectors, [], 1)};
            else
                vectors = reshape(vectors, 1, []);
            end

            axes = SplineAxis.empty(0,1);
            for iAxis = 1:numel(vectors)
                validateattributes(vectors{iAxis}, {'numeric'}, {'vector','real','finite','nonempty'});
                axes(end+1,1) = SplineAxis(values=reshape(vectors{iAxis}, [], 1)); %#ok<AGROW>
            end
        end

        function vectors = vectorsFromArray(axes)
            arguments
                axes (:,1) SplineAxis
            end

            vectors = cell(1, numel(axes));
            for iAxis = 1:numel(axes)
                vectors{iAxis} = axes(iAxis).values;
            end
        end
    end
end
