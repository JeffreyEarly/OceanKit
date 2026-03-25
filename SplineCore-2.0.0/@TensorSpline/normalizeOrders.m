function K = normalizeOrders(K, numDimensions)
% Normalize spline-order input to one order per dimension.
validateattributes(K, {'numeric'}, {'vector','real','finite','positive','integer'});
if isscalar(K)
    K = repmat(K, 1, numDimensions);
else
    K = reshape(K, 1, []);
    if numel(K) ~= numDimensions
        error('TensorSpline:InvalidOrderVector', 'K must be scalar or have one element per dimension.');
    end
end
