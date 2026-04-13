classdef TrajectorySpline < CAAnnotatedClass
    % Two-dimensional trajectory model parameterized by a shared scalar variable.
    %
    % `TrajectorySpline` stores a planar parametric trajectory as two
    % one-dimensional component splines,
    %
    % $$
    % x = x(t), \qquad y = y(t),
    % $$
    %
    % built from a shared parameter vector `t`. Use
    % `TrajectorySpline.fromData(...)` when raw coordinate samples should be
    % fit as one-dimensional `ConstrainedSpline` objects. The low-level
    % `TrajectorySpline(...)` constructor is the cheap canonical constructor
    % used for persisted restart and other direct bootstrap paths.
    %
    % ```matlab
    % t = linspace(0, 1, 20)';
    % x = cos(2*pi*t);
    % y = sin(2*pi*t);
    %
    % trajectory = TrajectorySpline.fromData(t, x, y, S=3);
    % xq = trajectory.x(t);
    % yq = trajectory.y(t);
    % ```
    %
    % - Topic: Create a trajectory spline
    % - Topic: Inspect trajectory properties
    % - Topic: Evaluate trajectory derivatives
    % - Declaration: classdef TrajectorySpline < CAAnnotatedClass

    properties (SetAccess = private)
        % Parameter samples shared by both component splines.
        %
        % `t` is stored as a column vector and provides the shared
        % parameterization for both coordinate splines `x(t)` and `y(t)`.
        %
        % - Topic: Inspect trajectory properties
        t (:,1) double {mustBeReal,mustBeFinite}

        % One-dimensional spline for the x-coordinate trajectory component.
        %
        % Evaluate the x-coordinate along the trajectory with
        % `trajectory.x(tq)`. When the trajectory is created through
        % `TrajectorySpline.fromData(...)`, this component is a
        % `ConstrainedSpline`.
        %
        % - Topic: Inspect trajectory properties
        x

        % One-dimensional spline for the y-coordinate trajectory component.
        %
        % Evaluate the y-coordinate along the trajectory with
        % `trajectory.y(tq)`. When the trajectory is created through
        % `TrajectorySpline.fromData(...)`, this component is a
        % `ConstrainedSpline`.
        %
        % - Topic: Inspect trajectory properties
        y
    end

    methods
        function self = TrajectorySpline(options)
            % Create a trajectory from canonical component-spline state.
            %
            % Use this low-level constructor when you already have the
            % shared trajectory parameter and the one-dimensional component
            % spline objects. For ordinary fitting from raw sample values,
            % use `TrajectorySpline.fromData(...)`.
            %
            % ```matlab
            % xSpline = ConstrainedSpline.fromData(t, x, S=3);
            % ySpline = ConstrainedSpline.fromData(t, y, S=3);
            % trajectory = TrajectorySpline(t=t, x=xSpline, y=ySpline);
            % ```
            %
            % - Topic: Create a trajectory spline
            % - Declaration: self = TrajectorySpline(options)
            % - Parameter options.t: strictly increasing shared trajectory parameter vector
            % - Parameter options.x: one-dimensional spline for the x-coordinate
            % - Parameter options.y: one-dimensional spline for the y-coordinate
            % - Returns self: TrajectorySpline instance
            arguments
                options.t (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.x (1,1) TensorSpline
                options.y (1,1) TensorSpline
            end

            if any(diff(options.t) <= 0)
                error('TrajectorySpline:NonmonotonicParameter', 't must be strictly increasing.');
            end

            if options.x.numDimensions ~= 1 || options.y.numDimensions ~= 1
                error('TrajectorySpline:InvalidComponentSpline', 'x and y must be one-dimensional splines.');
            end

            self@CAAnnotatedClass();
            self.t = options.t;
            self.x = options.x;
            self.y = options.y;
        end

        function values = u(self, t)
            % Evaluate the x-velocity $$u(t) = \dot{x}(t)$$ along the trajectory.
            %
            % Use this method when the x-component derivative should be
            % evaluated through the trajectory API rather than by reaching
            % into the component spline directly.
            %
            % - Topic: Evaluate trajectory derivatives
            % - Declaration: values = u(self,t)
            % - Parameter t: numeric query points with any shape
            % - Returns values: x-derivative values with the same shape as `t`
            arguments (Input)
                self (1,1) TrajectorySpline
                t {mustBeNumeric,mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values
            end

            values = self.x.valueAtPoints(t, D=1);
        end

        function values = v(self, t)
            % Evaluate the y-velocity $$v(t) = \dot{y}(t)$$ along the trajectory.
            %
            % Use this method when the y-component derivative should be
            % evaluated through the trajectory API rather than by reaching
            % into the component spline directly.
            %
            % - Topic: Evaluate trajectory derivatives
            % - Declaration: values = v(self,t)
            % - Parameter t: numeric query points with any shape
            % - Returns values: y-derivative values with the same shape as `t`
            arguments (Input)
                self (1,1) TrajectorySpline
                t {mustBeNumeric,mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values
            end

            values = self.y.valueAtPoints(t, D=1);
        end
    end

    methods (Static)
        function self = fromData(t, x, y, options)
            % Create a trajectory spline from raw x(t) and y(t) samples.
            %
            % Use this factory when the coordinate components of a planar
            % trajectory are known at the same parameter samples and should
            % each be fit as a one-dimensional `ConstrainedSpline`.
            %
            % The resulting component models satisfy
            %
            % $$
            % x(t_i) = x_i, \qquad y(t_i) = y_i,
            % $$
            %
            % for each supplied sample pair `x_i`, `y_i` at parameter value
            % `t_i`.
            %
            % ```matlab
            % t = linspace(0, 1, 20)';
            % x = cos(2*pi*t);
            % y = sin(2*pi*t);
            % trajectory = TrajectorySpline.fromData(t, x, y, S=3);
            % ```
            %
            % - Topic: Create a trajectory spline
            % - Declaration: self = fromData(t,x,y,options)
            % - Parameter t: strictly increasing shared trajectory parameter vector
            % - Parameter x: x-coordinate samples evaluated at `t`
            % - Parameter y: y-coordinate samples evaluated at `t`
            % - Parameter options.S: spline degree shared by both coordinate splines
            % - Returns self: TrajectorySpline instance
            arguments
                t (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                x (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                y (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.S (1,1) double {mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
            end

            if numel(t) ~= numel(x) || numel(t) ~= numel(y)
                error('TrajectorySpline:SizeMismatch', 't, x, and y must have the same number of elements.');
            end

            if any(diff(t) <= 0)
                error('TrajectorySpline:NonmonotonicParameter', 't must be strictly increasing.');
            end

            xSpline = ConstrainedSpline.fromData(t, x, S=options.S);
            ySpline = ConstrainedSpline.fromData(t, y, S=options.S);
            self = TrajectorySpline(t=t, x=xSpline, y=ySpline);
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
                error('TrajectorySpline:MissingAnnotatedClass', 'Unable to find the AnnotatedClass attribute in %s.', path);
            end
            self = TrajectorySpline.annotatedClassFromGroup(group);
        end

        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, {'t'});
            xSpline = TensorSpline.annotatedClassFromGroup(group.groupWithName('x'));
            ySpline = TensorSpline.annotatedClassFromGroup(group.groupWithName('y'));
            self = TrajectorySpline(t=vars.t, x=xSpline, y=ySpline);
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('t', '', 'Trajectory parameter samples.');
            propertyAnnotations(end+1) = CAObjectProperty('x', 'One-dimensional spline model for the x-coordinate trajectory component.');
            propertyAnnotations(end+1) = CAObjectProperty('y', 'One-dimensional spline model for the y-coordinate trajectory component.');
        end

        function names = classRequiredPropertyNames()
            names = {'t', 'x', 'y'};
        end
    end
end
