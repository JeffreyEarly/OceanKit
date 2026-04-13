function knotPoints = terminatedKnotPoints(knotPoints, S)
% Ensure each knot vector has S+1 repeated knots at its boundaries.
%
% Use this helper when you want to terminate a manually supplied
% knot sequence before fitting. In 1-D it accepts a numeric knot
% vector; in higher dimensions it accepts a cell array with one
% knot vector per dimension.
%
% This operation increases the multiplicity of the first and last knot
% values to `S+1`, which makes the spline basis terminate at the grid
% endpoints.
%
% ```matlab
% knotPoints = ConstrainedSpline.terminatedKnotPoints(knotPoints, 3);
% ```
%
% - Topic: Prepare knot sequences
% - Declaration: knotPoints = terminatedKnotPoints(knotPoints, S)
% - Parameter knotPoints: knot vector in 1-D or cell array of knot vectors
% - Parameter S: spline degree scalar or vector with one entry per dimension
% - Returns knotPoints: terminated knot sequence
if isnumeric(knotPoints)
    validateattributes(knotPoints, {'numeric'}, {'vector','real','finite'});
    K = TensorSpline.normalizeOrders(S + 1, 1);
    knotPoints = reshape(knotPoints, [], 1);

    nLeft = find(knotPoints <= knotPoints(1), 1, 'last');
    nRight = numel(knotPoints) - find(knotPoints == knotPoints(end), 1, 'first') + 1;
    knotPoints = [repmat(knotPoints(1), K-nLeft, 1); knotPoints; repmat(knotPoints(end), K-nRight, 1)];
    return;
end

if ~iscell(knotPoints)
    error('ConstrainedSpline:InvalidKnotCell',  'knotPoints must be a knot vector in 1-D or a cell array with one knot vector per dimension.');
end

numDimensions = numel(knotPoints);
K = TensorSpline.normalizeOrders(S + 1, numDimensions);
knotPoints = reshape(knotPoints, 1, []);
for iDim = 1:numDimensions
    knotPoints{iDim} = ConstrainedSpline.terminatedKnotPoints(knotPoints{iDim}, K(iDim) - 1);
end
