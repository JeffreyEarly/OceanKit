function derivativeOrders = normalizeDerivativeOrders(derivativeOrders, numDimensions)
% Normalize derivative-order input to one order per dimension.
%
% - Topic: Utility
% - Developer: true
% - Declaration: derivativeOrders = normalizeDerivativeOrders(derivativeOrders,numDimensions)
validateattributes(derivativeOrders, {'numeric'}, {'real','finite','nonnegative','integer'});
if isscalar(derivativeOrders)
    if numDimensions == 1
        derivativeOrders = derivativeOrders;
    elseif derivativeOrders == 0
        derivativeOrders = zeros(1, numDimensions);
    else
        error('TensorSpline:InvalidDerivativeOrders',  'Derivative orders must be a vector with one element per dimension.');
    end
else
    derivativeOrders = reshape(derivativeOrders, 1, []);
    if numel(derivativeOrders) ~= numDimensions
        error('TensorSpline:InvalidDerivativeOrders',  'Derivative orders must have one element per dimension.');
    end
end
