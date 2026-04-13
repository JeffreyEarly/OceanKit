function [pointConstraints, globalConstraints] = normalizeConstraintInputs(constraints, numDimensions)
% Split typed constraint inputs into point and global arrays.
%
% This helper partitions a mixed `SplineConstraint` array into its
% `PointConstraint` and `GlobalConstraint` components and checks that each
% constraint is dimensionally compatible with the target spline.
%
% - Topic: Prepare fit inputs
% - Developer: true
% - Declaration: [pointConstraints,globalConstraints] = normalizeConstraintInputs(constraints,numDimensions)
% - Parameter constraints: mixed SplineConstraint array
% - Parameter numDimensions: target spline dimensionality
% - Returns pointConstraints: PointConstraint array
% - Returns globalConstraints: GlobalConstraint array
if isempty(constraints)
    pointConstraints = PointConstraint.empty(0,1);
    globalConstraints = GlobalConstraint.empty(0,1);
    return;
end

constraints = reshape(constraints, [], 1);
isPointConstraint = arrayfun(@(constraint) isa(constraint, 'PointConstraint'), constraints);
if any(isPointConstraint)
    pointConstraints = reshape([constraints(isPointConstraint)], [], 1);
else
    pointConstraints = PointConstraint.empty(0,1);
end

if any(~isPointConstraint)
    globalConstraints = reshape([constraints(~isPointConstraint)], [], 1);
else
    globalConstraints = GlobalConstraint.empty(0,1);
end

for iConstraint = 1:numel(pointConstraints)
    if ~isempty(numDimensions) && pointConstraints(iConstraint).numDimensions ~= numDimensions
        error('ConstrainedSpline:PointConstraintDimensionMismatch',  'Each point constraint must match the spline dimensionality.');
    end
end

for iConstraint = 1:numel(globalConstraints)
    if ~isempty(numDimensions) && ~isempty(globalConstraints(iConstraint).dimension) && globalConstraints(iConstraint).dimension > numDimensions
        error('ConstrainedSpline:GlobalConstraintDimensionMismatch',  'Global constraint dimensions must not exceed the spline dimensionality.');
    end
end
