function [pointMatrix, supportVectors] = pointsOfSupportFromKnotPoints(knotPoints, options)
% Return representative support points for a tensor-product spline basis.
%
% Use these points when you need one representative location per tensor
% basis function, for example when constructing transformed splines from
% sampled values.
%
% The support vectors are computed one dimension at a time with
% `BSpline.pointsOfSupport`, then combined by a Cartesian product to form
% the returned point matrix.
%
% ```matlab
% [supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(knotPoints, S=[3 3]);
% values = spline(supportVectors{:});
% ```
%
% - Topic: Build spline bases
% - Declaration: [pointMatrix,supportVectors] = pointsOfSupportFromKnotPoints(knotPoints, options)
% - Parameter knotPoints: knot vector in 1-D or cell array of knot vectors
% - Parameter options.S: spline degree scalar or vector with one entry per dimension
% - Returns pointMatrix: matrix with one row per tensor support point
% - Returns supportVectors: cell array with one support vector per dimension
arguments
    knotPoints
    options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
end
S = options.S;

if isnumeric(knotPoints)
    validateattributes(knotPoints, {'numeric'}, {'vector','real','finite','nonempty'});
    numDimensions = 1;
    tKnot = {reshape(knotPoints, [], 1)};
else
    numDimensions = numel(knotPoints);
    tKnot = TensorSpline.normalizeKnotCell(knotPoints, numDimensions);
end
K = TensorSpline.normalizeOrders(S + 1, numDimensions);

supportVectors = cell(1, numDimensions);
for iDim = 1:numDimensions
    supportVectors{iDim} = BSpline.pointsOfSupportFromKnotPoints(tKnot{iDim}, S=K(iDim) - 1);
end

pointMatrix = TensorSpline.pointsFromGridVectors(supportVectors);
