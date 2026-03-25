function f = evaluateFromPPCoefficients(options)
% Evaluate a cached piecewise-polynomial spline representation.
%
% On interval `i`, let `u = queryPoints - tpp(i)`. This function evaluates
%
% $$
% f_i^{(D)}(u) = \sum_{m=D}^{S} \frac{c_{i,m}}{(m-D)!} u^{m-D},
% $$
%
% where the interval coefficients $$c_{i,m}$$ are stored in `C`.
%
% ```matlab
% xq = BSpline.evaluateFromPPCoefficients(queryPoints=tQuery, C=C, tpp=tpp);
% dxq = BSpline.evaluateFromPPCoefficients(queryPoints=tQuery, C=C, tpp=tpp, D=1);
% ```
%
% - Topic: Represent piecewise polynomials
% - Developer: true
% - Declaration: f = evaluateFromPPCoefficients(options)
% - Parameter options.queryPoints: points at which to evaluate the splines
% - Parameter options.C: polynomial coefficients to be used in polyval, size(C) = [length(tpp)-1, K]
% - Parameter options.tpp: piece-wise polynomial intervals
% - Parameter options.D: number of derivatives
% - Returns f: array the same size as queryPoints
arguments
    options.queryPoints {mustBeNumeric,mustBeReal}
    options.C (:,:) double
    options.tpp (:,1) double {mustBeNumeric,mustBeReal}
    options.D (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end
queryPoints = options.queryPoints;
C = options.C;
tpp = options.tpp;
D = options.D;

K = size(C,2);
f = zeros(size(queryPoints), 'like', queryPoints);

if D > K-1
    % By construction the splines are zero for K or more derivs
    return;
end

scale = factorial((K-1-D):-1:0);
indices = 1:(K-D);
scaledC = C(:,indices)./scale;

% Evaluate on a sorted copy so interval bins are contiguous, then restore
% the original ordering and shape.
[tSorted, sortIndices] = sort(queryPoints(:), 'ascend');
t_pp_bin = discretize(tSorted, [-Inf; tpp(2:end-1); Inf]);
fSorted = zeros(size(tSorted), 'like', tSorted);

startIndex = 1;
while startIndex <= numel(tSorted)
    iBin = t_pp_bin(startIndex);
    endIndex = startIndex + find(t_pp_bin(startIndex:end) ~= iBin, 1, 'first') - 2;
    if isempty(endIndex)
        endIndex = numel(tSorted);
    end
    coeffs = scaledC(iBin,:);
    fSorted(startIndex:endIndex) = polyval(coeffs, tSorted(startIndex:endIndex) - tpp(iBin));
    startIndex = endIndex + 1;
end

fFlat = zeros(numel(queryPoints), 1, 'like', queryPoints);
fFlat(sortIndices) = fSorted;
f = reshape(fFlat, size(queryPoints));
end
