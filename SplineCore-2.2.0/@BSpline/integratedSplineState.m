function [xiIntegrated, knotPointsIntegrated, SIntegrated] = integratedSplineState(xi, options)
% Return the canonical spline state of the zero-anchored antiderivative.
%
% Use this expert utility when B-spline coefficients are already known and
% the exact antiderivative state should be computed without constructing a
% temporary `BSpline` object first.
%
% For coefficient matrix `xi`, each column is integrated independently
% with the shared affine output normalization
%
% $$
% f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau),
% $$
%
% producing the canonical antiderivative spline state
%
% $$
% F(t) = \int_{\tau_1}^{t} f(s)\,ds.
% $$
%
% ```matlab
% [xiInt, knotPointsInt, SInt] = BSpline.integratedSplineState( ...
%     xi, knotPoints=knotPoints, S=3);
% intspline = BSpline(S=SInt, knotPoints=knotPointsInt, xi=xiInt);
% ```
%
% - Topic: Transform the spline
% - Declaration: [xiIntegrated,knotPointsIntegrated,SIntegrated] = integratedSplineState(xi, options)
% - Parameter xi: spline coefficient vector or matrix with one spline per column
% - Parameter options.knotPoints: terminated knot sequence for the input spline basis
% - Parameter options.S: spline degree of the input coefficients
% - Parameter options.xMean: optional additive output offset shared by every column
% - Parameter options.xStd: optional multiplicative output scale shared by every column
% - Returns xiIntegrated: antiderivative spline coefficient matrix with one extra row
% - Returns knotPointsIntegrated: terminated knot sequence for the antiderivative spline
% - Returns SIntegrated: spline degree of the antiderivative spline
arguments
    xi (:,:) double {mustBeReal,mustBeFinite,mustBeNonempty}
    options.knotPoints (:,1) double {mustBeReal,mustBeFinite,mustBeNonempty}
    options.S (1,1) double {mustBeInteger,mustBeNonnegative}
    options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
    options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
end

knotPoints = options.knotPoints;
S = options.S;
if any(diff(knotPoints) < 0)
    error('BSpline:InvalidKnotPoints', 'knotPoints must be non-decreasing.');
end

numCoefficients = size(xi, 1);
expectedCoefficientCount = numel(knotPoints) - S - 1;
if numCoefficients ~= expectedCoefficientCount
    error('BSpline:InvalidCoefficientCount', ...
        'xi must have %d rows for knotPoints with degree S=%d.', expectedCoefficientCount, S);
end

if abs(options.xMean) > 0 || abs(options.xStd - 1) > 0
    supportPoints = BSpline.pointsOfSupportFromKnotPoints(knotPoints, S=S);
    basisMatrix = BSpline.matrixForDataPoints(supportPoints, knotPoints=knotPoints, S=S);
    xi = options.xStd * xi + basisMatrix \ (options.xMean * ones(numel(supportPoints), 1));
end

K = S + 1;
dt = (knotPoints(1 + K:numCoefficients + K) - knotPoints(1:numCoefficients)) / K;
xiIntegrated = [zeros(1, size(xi, 2), 'like', xi); cumsum(xi .* reshape(dt, [], 1), 1)];
knotPointsIntegrated = [knotPoints(1); knotPoints; knotPoints(end)];
SIntegrated = S + 1;
end
