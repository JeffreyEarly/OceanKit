classdef PointConstraint < SplineConstraint
    % Specify local equality or bound constraints at one or more points.
    %
    % Use `PointConstraint` to declare values or derivative conditions at
    % specific points in one or more dimensions. A single constraint object
    % can represent many points at once, which makes it suitable for both
    % simple one-dimensional constraints and large masked regions in tensor
    % fits.
    %
    % ## Basic usage
    %
    % ```matlab
    % c1 = PointConstraint.equal((0:10)', D=2, value=0);
    % c2 = PointConstraint.lowerBound([X(mask), Y(mask)], D=[0 0], value=0);
    % c3 = PointConstraint.equalOnMask({x,y}, islandMask, D=[0 0], value=0);
    % ```
    %
    % - Topic: Specify point constraints
    % - Declaration: classdef PointConstraint

    properties (Constant, Hidden)
        allowedRelations = ["==", ">=", "<="]
        equalRelation = "=="
        lowerBoundRelation = ">="
        upperBoundRelation = "<="
    end

    properties (SetAccess = private)
        % Constraint locations as an N-by-D point matrix.
        %
        % - Topic: Inspect point constraint properties
        points = zeros(0,0)

        % Derivative orders as an N-by-D matrix.
        %
        % - Topic: Inspect point constraint properties
        D = zeros(0,0)

        % Constraint relation: "==", ">=", or "<=".
        %
        % - Topic: Inspect point constraint properties
        relation string = PointConstraint.equalRelation

        % Target values as an N-by-1 vector.
        %
        % - Topic: Inspect point constraint properties
        value = zeros(0,1)
    end

    properties (Dependent)
        % Number of constrained dimensions.
        %
        % - Topic: Inspect point constraint properties
        numDimensions

        % Number of constrained points.
        %
        % - Topic: Inspect point constraint properties
        numConstraints
    end

    properties (Dependent, Hidden, SetAccess = private)
        constraintIndex
        constraintDimension
    end

    methods
        function self = PointConstraint(points, relation, value, options)
            % Create a pointwise equality or bound constraint.
            %
            % Use the static helper methods `equal`, `lowerBound`, and
            % `upperBound` for the intended public construction style.
            %
            % - Topic: Specify point constraints
            % - Declaration: self = PointConstraint(points,relation,value,options)
            % - Parameter points: point locations as a vector or N-by-D matrix
            % - Parameter relation: one of "==", ">=", or "<="
            % - Parameter value: scalar or one target value per point
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Returns self: PointConstraint instance
            arguments
                points = []
                relation {mustBeTextScalar} = PointConstraint.equalRelation
                value = []
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
            end

            derivativeOrders = options.D;
            relation = string(relation);
            if ~any(relation == PointConstraint.allowedRelations)
                error('PointConstraint:InvalidRelation', 'relation must be one of "==", ">=", or "<=".');
            end

            self@SplineConstraint();

            if isempty(points)
                return;
            end

            pointMatrix = PointConstraint.normalizePoints(points);
            numPoints = size(pointMatrix,1);
            numDimensions = size(pointMatrix,2);
            derivativeOrders = PointConstraint.normalizeDerivativeOrders(derivativeOrders, numPoints, numDimensions);
            targetValues = PointConstraint.normalizeValues(value, numPoints);
            self.points = pointMatrix;
            self.D = derivativeOrders;
            self.relation = relation;
            self.value = targetValues;
        end

        function value = get.numDimensions(self)
            % Return the number of constrained dimensions.
            %
            % - Topic: Inspect point constraint properties
            % - Declaration: value = get.numDimensions(self)
            % - Parameter self: PointConstraint instance
            % - Returns value: number of dimensions in points
            value = size(self.points, 2);
        end

        function value = get.numConstraints(self)
            % Return the number of constrained points.
            %
            % - Topic: Inspect point constraint properties
            % - Declaration: value = get.numConstraints(self)
            % - Parameter self: PointConstraint instance
            % - Returns value: number of rows in points
            value = size(self.points, 1);
        end

        function value = get.constraintIndex(self)
            % Return the index coordinate for persisted point constraints.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.constraintIndex(self)
            % - Parameter self: PointConstraint instance
            % - Returns value: constraint index vector
            value = reshape(1:self.numConstraints, [], 1);
        end

        function value = get.constraintDimension(self)
            % Return the coordinate-dimension index for persisted point-constraint matrices.
            %
            % - Topic: Persist spline state
            % - Developer: true
            % - Declaration: value = get.constraintDimension(self)
            % - Parameter self: PointConstraint instance
            % - Returns value: point-dimension index vector
            value = reshape(1:self.numDimensions, [], 1);
        end
    end

    methods (Static)
        function self = equal(points, options)
            % Create a pointwise equality constraint.
            %
            % ```matlab
            % c = PointConstraint.equal(tc, D=2, value=0);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = equal(points,options)
            % - Parameter points: point locations as a vector or N-by-D matrix
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one target value per point
            % - Returns self: equality PointConstraint
            arguments
                points {mustBeNumeric,mustBeReal,mustBeFinite}
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(points, PointConstraint.equalRelation, options.value, D=options.D);
        end

        function self = lowerBound(points, options)
            % Create a pointwise lower-bound constraint.
            %
            % ```matlab
            % c = PointConstraint.lowerBound(P, D=[0 1], value=0);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = lowerBound(points,options)
            % - Parameter points: point locations as a vector or N-by-D matrix
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one bound value per point
            % - Returns self: lower-bound PointConstraint
            arguments
                points {mustBeNumeric,mustBeReal,mustBeFinite}
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(points, PointConstraint.lowerBoundRelation, options.value, D=options.D);
        end

        function self = upperBound(points, options)
            % Create a pointwise upper-bound constraint.
            %
            % ```matlab
            % c = PointConstraint.upperBound(P, D=[0 0], value=1);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = upperBound(points,options)
            % - Parameter points: point locations as a vector or N-by-D matrix
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one bound value per point
            % - Returns self: upper-bound PointConstraint
            arguments
                points {mustBeNumeric,mustBeReal,mustBeFinite}
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(points, PointConstraint.upperBoundRelation, options.value, D=options.D);
        end

        function self = equalOnMask(grid, mask, options)
            % Create a pointwise equality constraint from a logical mask.
            %
            % Use this when a constrained region is naturally described by
            % a logical mask on a rectilinear grid.
            %
            % ```matlab
            % c = PointConstraint.equalOnMask({x,y}, mask, D=[0 0], value=0);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = equalOnMask(grid,mask,options)
            % - Parameter grid: vector or cell array of grid vectors
            % - Parameter mask: logical mask selecting constrained locations
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one target value per selected point
            % - Returns self: equality PointConstraint
            arguments
                grid
                mask
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(PointConstraint.pointsFromMask(grid, mask), PointConstraint.equalRelation, options.value, D=options.D);
        end

        function self = lowerBoundOnMask(grid, mask, options)
            % Create a pointwise lower-bound constraint from a logical mask.
            %
            % ```matlab
            % c = PointConstraint.lowerBoundOnMask({x,y}, mask, D=[0 1], value=0);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = lowerBoundOnMask(grid,mask,options)
            % - Parameter grid: vector or cell array of grid vectors
            % - Parameter mask: logical mask selecting constrained locations
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one bound value per selected point
            % - Returns self: lower-bound PointConstraint
            arguments
                grid
                mask
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(PointConstraint.pointsFromMask(grid, mask), PointConstraint.lowerBoundRelation, options.value, D=options.D);
        end

        function self = upperBoundOnMask(grid, mask, options)
            % Create a pointwise upper-bound constraint from a logical mask.
            %
            % ```matlab
            % c = PointConstraint.upperBoundOnMask({x,y}, mask, D=[0 0], value=1);
            % ```
            %
            % - Topic: Specify point constraints
            % - Declaration: self = upperBoundOnMask(grid,mask,options)
            % - Parameter grid: vector or cell array of grid vectors
            % - Parameter mask: logical mask selecting constrained locations
            % - Parameter options.D: derivative orders as a scalar, row vector, or N-by-D matrix
            % - Parameter options.value: scalar or one bound value per selected point
            % - Returns self: upper-bound PointConstraint
            arguments
                grid
                mask
                options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 0
                options.value {mustBeNumeric,mustBeReal,mustBeFinite} = 0
            end

            self = PointConstraint(PointConstraint.pointsFromMask(grid, mask), PointConstraint.upperBoundRelation, options.value, D=options.D);
        end
    end

    methods (Static, Hidden)
        function self = annotatedClassFromFile(path)
            ncfile = NetCDFFile(path, shouldReadOnly=true);
            cleanup = onCleanup(@() ncfile.close()); %#ok<NASGU>
            if isKey(ncfile.attributes, 'AnnotatedClass')
                className = string(ncfile.attributes('AnnotatedClass'));
                if ncfile.hasGroupWithName(className)
                    group = ncfile.groupWithName(className);
                else
                    group = ncfile;
                end
            else
                error('PointConstraint:MissingAnnotatedClass', 'Unable to find the AnnotatedClass attribute in %s.', path);
            end
            self = PointConstraint.annotatedClassFromGroup(group);
        end

        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, PointConstraint.classRequiredPropertyNames());
            self = PointConstraint(vars.points, vars.relation, vars.value, D=vars.D);
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = SplineConstraint.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CADimensionProperty('constraintIndex', '', 'Index over constrained points.');
            propertyAnnotations(end+1) = CADimensionProperty('constraintDimension', '', 'Index over constrained dimensions.');
            propertyAnnotations(end+1) = CANumericProperty('points', {'constraintIndex', 'constraintDimension'}, '', 'Constraint locations as an N-by-D point matrix.');
            propertyAnnotations(end+1) = CANumericProperty('D', {'constraintIndex', 'constraintDimension'}, '', 'Derivative orders as an N-by-D matrix.');
            propertyAnnotations(end+1) = CAPropertyAnnotation('relation', 'Constraint relation.');
            propertyAnnotations(end+1) = CANumericProperty('value', {'constraintIndex'}, '', 'Target values for each constrained point.');
        end

        function names = classRequiredPropertyNames()
            names = {'points', 'relation', 'value', 'D'};
        end
    end

    methods (Static, Access = private)
        function pointMatrix = normalizePoints(points)
            % Normalize point input to an N-by-D matrix.
            validateattributes(points, {'numeric'}, {'2d','real','finite','nonempty'});
            if isvector(points)
                pointMatrix = reshape(points, [], 1);
            else
                pointMatrix = points;
            end
        end

        function derivativeOrders = normalizeDerivativeOrders(D, numPoints, numDimensions)
            % Normalize derivative orders to an N-by-D matrix.
            validateattributes(D, {'numeric'}, {'2d','real','finite','nonnegative','integer','nonempty'});

            if numDimensions == 1
                if isscalar(D)
                    derivativeOrders = repmat(double(D), numPoints, 1);
                    return;
                end

                if isvector(D) && numel(D) == numPoints
                    derivativeOrders = reshape(double(D), [], 1);
                    return;
                end

                error('PointConstraint:InvalidDerivativeOrders',  'For one-dimensional constraints, D must be scalar or have one entry per point.');
            end

            if isvector(D)
                D = reshape(D, 1, []);
                if numel(D) == numDimensions
                    derivativeOrders = repmat(double(D), numPoints, 1);
                    return;
                end

                if isscalar(D)
                    error('PointConstraint:AmbiguousDerivativeOrders',  'For multi-dimensional constraints, scalar D is ambiguous. Supply one derivative order per dimension.');
                end
            end

            if isequal(size(D), [numPoints, numDimensions])
                derivativeOrders = double(D);
                return;
            end

            error('PointConstraint:InvalidDerivativeOrders',  'D must be a 1-by-numDimensions vector or an N-by-numDimensions matrix.');
        end

        function values = normalizeValues(value, numPoints)
            % Normalize target values to an N-by-1 vector.
            validateattributes(value, {'numeric'}, {'vector','real','finite','nonempty'});
            if isscalar(value)
                values = repmat(double(value), numPoints, 1);
                return;
            end

            if numel(value) ~= numPoints
                error('PointConstraint:InvalidValueSize',  'value must be scalar or provide one entry per constrained point.');
            end

            values = reshape(double(value), [], 1);
        end

        function pointMatrix = pointsFromMask(grid, mask)
            % Convert a masked grid description into an explicit point matrix.
            if iscell(grid)
                if isempty(grid)
                    error('PointConstraint:InvalidGrid',  'grid must not be empty.');
                end

                if ~all(cellfun(@isvector, grid))
                    error('PointConstraint:InvalidGrid',  'grid must be a vector or a cell array of grid vectors.');
                end

                [allPoints, gridSize] = TensorSpline.pointsFromGridVectors(grid);
                normalizedMask = PointConstraint.normalizeMask(mask, gridSize);
                pointMatrix = allPoints(normalizedMask(:), :);
                return;
            end

            validateattributes(grid, {'numeric'}, {'vector','real','finite','nonempty'});
            normalizedMask = PointConstraint.normalizeMask(mask, size(grid));
            pointMatrix = reshape(grid(normalizedMask), [], 1);
        end

        function mask = normalizeMask(mask, expectedSize)
            % Normalize a logical mask to the requested output size.
            validateattributes(mask, {'numeric','logical'}, {'nonempty'});

            if isvector(mask) && prod(expectedSize) == numel(mask)
                mask = reshape(logical(mask), expectedSize);
            elseif isequal(size(mask), expectedSize)
                mask = logical(mask);
            else
                error('PointConstraint:InvalidMaskSize',  'mask must match the supplied grid size.');
            end

            if ~any(mask, 'all')
                error('PointConstraint:EmptyMask',  'mask must select at least one constrained point.');
            end
        end
    end
end
