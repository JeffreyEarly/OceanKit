function basisSize = basisSizeFromKnotCell(tKnot, K)
% Compute basis sizes from knot vectors and spline orders.
basisSize = reshape(cellfun(@numel, tKnot), 1, []) - reshape(K, 1, []);
if any(basisSize <= 0)
    error('TensorSpline:InvalidBasisSize', 'Each knot vector must be longer than its spline order.');
end
