classdef (Abstract) SplineConstraint < matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Common superclass for local and global spline constraint objects.
    %
    % Use `SplineConstraint` when you want to pass a mixed array of
    % `PointConstraint` and `GlobalConstraint` objects through one API.
    %
    % ```matlab
    % constraints = [
    %     PointConstraint.equal(0, D=1, value=0)
    %     GlobalConstraint.positive()
    % ];
    % ```
    %
    % - Topic: Specify constraints
    % - Declaration: classdef SplineConstraint < matlab.mixin.Heterogeneous & CAAnnotatedClass

    methods
        function self = SplineConstraint()
            % Initialize the heterogeneous spline-constraint root.
            %
            % - Topic: Specify constraints
            % - Declaration: self = SplineConstraint()
            % - Returns self: SplineConstraint base-class instance
            self@CAAnnotatedClass();
        end
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end

        function names = classRequiredPropertyNames()
            names = {};
        end
    end
end
