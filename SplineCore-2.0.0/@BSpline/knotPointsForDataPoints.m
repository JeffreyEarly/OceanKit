function knotPoints = knotPointsForDataPoints(dataPoints, options)
% Construct a terminated knot sequence from sample locations.
%
% Use this helper to choose a knot sequence directly from sample
% locations before interpolation or least-squares fitting.
%
% The implementation sorts the sample locations, optionally subsamples them
% to target `splineDOF`, constructs pseudo-sites on that reduced grid, then
% repeats the first and last knot values `S+1` times so the spline basis is
% terminated at the interval endpoints.
%
% When `splineDOF` is supplied, the current implementation chooses the
% sample stride as
%
% $$
% \Delta n = \left\lceil \frac{N}{\max(\mathrm{splineDOF}, S+1)} \right\rceil,
% $$
%
% where $$N = \mathrm{numel}(t)$$.
%
% To recover the old `dataDOF=d` behavior, convert it to
% `splineDOF = max(S+1, ceil(numel(t)/d))`:
%
% ```matlab
% oldDataDOF = 3;
% S = 3;
% splineDOF = max(S+1, ceil(numel(t)/oldDataDOF));
% knotPoints = BSpline.knotPointsForDataPoints(t, S=S, splineDOF=splineDOF);
% ```
%
% ```matlab
% knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
% X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3);
% xi = X \ x;
% spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
% ```
%
% - Topic: Build spline bases
% - Declaration: knotPoints = knotPointsForDataPoints(dataPoints, options)
% - Parameter dataPoints: observation times (N)
% - Parameter options.S: (optional) spline degree
% - Parameter options.splineDOF: (optional) approximate target number of splines
% - Returns knotPoints: vector of knot point locations
arguments
    dataPoints (:,1) double
    options.S (1,1) double {mustBeNonnegative,mustBeInteger} = 3
    options.splineDOF (1,1) double = NaN
end
K = options.S + 1;

if ~isnan(options.splineDOF)
    mustBePositive(options.splineDOF);
    mustBeInteger(options.splineDOF);
end

if ~isnan(options.splineDOF)
    targetSplines = max(options.splineDOF, K);
    sampleStride = ceil(numel(dataPoints)/targetSplines);
else
    sampleStride = 1;
end

tData = sort(dataPoints);
tData = [tData(1); tData(1+sampleStride:sampleStride:end-sampleStride); tData(end)];
M = numel(tData);
mustBeGreaterThanOrEqual(M, K);

N = length(tData);
t_pseudo = interp1((0:N-1)',tData,linspace(0,N-1,M).');

if mod(K,2) == 1
    % Odd spline order, so knots go in between points.
    dt = diff(t_pseudo);

    % This gives us M+1 knot points.
    knotPoints = [t_pseudo(1); t_pseudo(1:end-1)+dt/2; t_pseudo(end)];

    % Now remove start and end knots
    for i=1:((K-1)/2)
        knotPoints(2) = [];
        knotPoints(end-1) = [];
    end

else
    knotPoints = t_pseudo;

    % Now remove start and end knots
    for i=1:((K-2)/2)
        knotPoints(2) = [];
        knotPoints(end-1) = [];
    end

end
% Now we increase the multiplicity of the knot points at the beginning and
% the end of the interval so that the splines do not extend past the end
% points.
knotPoints = [repmat(knotPoints(1),K-1,1); knotPoints; repmat(knotPoints(end),K-1,1)];
end
