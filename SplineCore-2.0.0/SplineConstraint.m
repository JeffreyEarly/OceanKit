classdef (Abstract) SplineConstraint < matlab.mixin.Heterogeneous
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
    % - Declaration: classdef SplineConstraint < matlab.mixin.Heterogeneous
end
