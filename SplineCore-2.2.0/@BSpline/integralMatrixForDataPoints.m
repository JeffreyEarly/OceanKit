function W = integralMatrixForDataPoints(dataPoints, queryPoints, options)
% Return the exact interpolation-to-antiderivative linear map.
%
% Use this expert utility when one-dimensional interpolation values are
% sampled on `dataPoints` and the zero-anchored antiderivative should be
% evaluated at `queryPoints` without constructing an intermediate spline
% object.
%
% If `y` is the column vector of interpolation values and `f` is the
% corresponding exact interpolating spline on the supplied terminated knot
% sequence, then the returned matrix satisfies
%
% $$
% W y = F(t_q), \qquad
% F(t) = \int_{t_1}^{t} f(s)\,ds.
% $$
%
% ```matlab
% W = BSpline.integralMatrixForDataPoints(t, tq, knotPoints=knotPoints, S=3);
% valuesInt = W * values;
% ```
%
% - Topic: Build spline bases
% - Declaration: W = integralMatrixForDataPoints(dataPoints, queryPoints, options)
% - Parameter dataPoints: interpolation sample locations
% - Parameter queryPoints: points at which to evaluate the antiderivative
% - Parameter options.knotPoints: terminated knot sequence for the interpolation basis
% - Parameter options.S: spline degree of the interpolation basis
% - Returns W: linear map from interpolation values on `dataPoints` to antiderivative values on `queryPoints`
arguments
    dataPoints (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
    queryPoints (:,1) double {mustBeReal,mustBeFinite}
    options.knotPoints (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
    options.S (1,1) double {mustBeInteger,mustBeNonnegative}
end

knotPoints = options.knotPoints;
S = options.S;
basisMatrix = BSpline.matrixForDataPoints(dataPoints, knotPoints=knotPoints, S=S);
numCoefficients = size(basisMatrix, 2);
if numCoefficients ~= numel(dataPoints)
    error('BSpline:InvalidInterpolationSystem', ...
        'knotPoints and S must define a one-dimensional interpolation basis with one coefficient per data point.');
end

numDataPoints = numel(dataPoints);
meanWeights = ones(1, numDataPoints) / numDataPoints;
constantSupportPoints = BSpline.pointsOfSupportFromKnotPoints(knotPoints, S=S);
constantBasisMatrix = BSpline.matrixForDataPoints(constantSupportPoints, knotPoints=knotPoints, S=S);
coefficientMap = basisMatrix \ (eye(numDataPoints) - ones(numDataPoints, 1) * meanWeights);
coefficientMap = coefficientMap + (constantBasisMatrix \ ones(numCoefficients, 1)) * meanWeights;

K = S + 1;
dt = (knotPoints(1 + K:numCoefficients + K) - knotPoints(1:numCoefficients)) / K;
integrationMatrix = tril(repmat(reshape(dt, 1, []), numCoefficients + 1, 1), -1);

knotPointsIntegrated = [knotPoints(1); knotPoints; knotPoints(end)];
integratedBasisMatrix = BSpline.matrixForDataPoints(queryPoints, knotPoints=knotPointsIntegrated, S=S + 1);
W = integratedBasisMatrix * (integrationMatrix * coefficientMap);
end
