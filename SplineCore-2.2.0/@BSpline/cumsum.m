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

[xiIntegrated, knotPointsIntegrated, SIntegrated] = BSpline.integratedSplineState( ...
    spline.xi, knotPoints=spline.knotPoints, S=spline.S, xMean=spline.xMean, xStd=spline.xStd);
intspline = BSpline(S=SIntegrated, knotPoints=knotPointsIntegrated, xi=xiIntegrated);
