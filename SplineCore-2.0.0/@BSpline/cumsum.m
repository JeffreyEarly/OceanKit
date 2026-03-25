function intspline = cumsum(spline)
% Return the indefinite integral of a B-spline.
%
% Use this to construct an antiderivative spline that can be evaluated at
% arbitrary points after integration.
%
% For coefficient vector $$\xi$$, the integrated spline uses cumulative
% coefficients
%
% $$
% \beta_0 = 0, \qquad \beta_j = \sum_{m=1}^{j} \xi_m \frac{\tau_{m+K} - \tau_m}{K}.
% $$
%
% If the spline carries nontrivial `xMean` or `xStd`, the method first
% converts that affine output normalization into an equivalent coefficient
% representation before integrating.
%
% ```matlab
% F = cumsum(spline);
% values = F(tQuery);
% ```
%
% - Topic: Transform the spline
% - Declaration: intspline = cumsum(spline)
% - Parameter spline: BSpline instance to integrate
% - Returns intspline: BSpline representing the antiderivative
arguments
    spline (1,1) BSpline
end

xi = spline.xi;
K = spline.K;
knotPoints = spline.knotPoints;
M = length(xi);

if abs(spline.xMean) > 0 || abs(spline.xStd - 1) > 0
    t = BSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
    X = BSpline.matrixForDataPoints(t, knotPoints=spline.knotPoints, S=spline.S);
    xi = spline.xStd*spline.xi + X\(spline.xMean*ones(length(t),1));
end

dt = (knotPoints(1+K:M+K)-knotPoints(1:M))/K;
beta = [0; cumsum(xi.*dt)];

knotPoints = cat(1, spline.knotPoints(1), spline.knotPoints, spline.knotPoints(end));
intspline = BSpline(S=spline.S+1, knotPoints=knotPoints, xi=beta);
