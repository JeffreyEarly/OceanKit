function B = matrixForPointMatrix(pointMatrix, options)
% Evaluate the tensor-product basis matrix and optional derivatives.
%
% Use this to assemble a tensor-product design matrix for
% interpolation, regression, or basis inspection.
%
% If `pointMatrix` has rows `x_i`, then row `i` of the returned matrix is the
% Kronecker product of the one-dimensional basis rows evaluated at `x_i`.
% In other words,
%
% $$
% \mathbf{B}_{i,:} =
% B^{(1)}(x_{i,1}) \otimes \cdots \otimes B^{(d)}(x_{i,d}),
% $$
%
% where each factor is the one-dimensional basis row in one coordinate
% direction, optionally replaced by its derivative row.
%
% ```matlab
% [Xq, Yq] = ndgrid(xq, yq);
% B = TensorSpline.matrixForPointMatrix([Xq(:), Yq(:)], knotPoints=knotPoints, S=[3 3]);
% values = B * spline.xi(:);
% ```
%
% - Topic: Build spline bases
% - Declaration: B = matrixForPointMatrix(pointMatrix, options)
% - Parameter pointMatrix: query locations as a point matrix
% - Parameter options.knotPoints: knot vector in 1-D or cell array of knot vectors
% - Parameter options.S: spline degree scalar or vector with one entry per dimension
% - Parameter options.D: derivative order per dimension
% - Returns B: basis matrix with one row per query point
arguments
    pointMatrix {mustBeNumeric,mustBeReal}
    options.knotPoints
    options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
    options.D = 0
end
knotPoints = options.knotPoints;
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
derivativeOrders = TensorSpline.normalizeDerivativeOrders(options.D, numDimensions);
basisSize = TensorSpline.basisSizeFromKnotCell(tKnot, K);

if numDimensions == 1
    pointMatrix = reshape(pointMatrix, [], 1);
else
    if size(pointMatrix,2) ~= numDimensions
        error('TensorSpline:InvalidPointMatrix', 'Point matrix must have one column per dimension.');
    end
end

if any(derivativeOrders > K - 1)
    B = zeros(size(pointMatrix,1), prod(basisSize));
    return;
end

numPoints = size(pointMatrix,1);
dimensionMatrices = cell(1, numDimensions);
for iDim = 1:numDimensions
    Bi = BSpline.matrixForDataPoints(pointMatrix(:,iDim), knotPoints=tKnot{iDim}, S=K(iDim) - 1, D=derivativeOrders(iDim));
    dimensionMatrices{iDim} = reshape(Bi(:,:,derivativeOrders(iDim)+1), numPoints, []);
end

B = dimensionMatrices{1};
for iDim = 2:numDimensions
    previousBasis = B;
    currentBasis = dimensionMatrices{iDim};
    combinedBasis = zeros(numPoints, size(previousBasis,2) * size(currentBasis,2));
    for iPoint = 1:numPoints
        combinedBasis(iPoint,:) = kron(currentBasis(iPoint,:), previousBasis(iPoint,:));
    end
    B = combinedBasis;
end
