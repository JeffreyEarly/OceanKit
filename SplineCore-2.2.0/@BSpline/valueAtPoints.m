function x_out = valueAtPoints(self, t, options)
% Evaluate the spline or one of its derivatives at arbitrary points.
%
% This is the primary explicit evaluation method. Pass
% `D=0` for spline values, `D=1` for the first
% derivative, and so on.
%
% The returned array represents
%
% $$
% f^{(D)}(t) = x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}^{(D)}(t;\tau),
% $$
%
% with `xMean` added back only when `D=0`.
%
% ```matlab
% x = spline.valueAtPoints(tQuery);
% d2x = spline.valueAtPoints(tQuery, D=2);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: x_out = valueAtPoints(self,t,options)
% - Parameter self: BSpline instance
% - Parameter t: evaluation points
% - Parameter options.D: derivative order to evaluate
% - Returns x_out: array matching the shape of t
arguments
    self (1,1) BSpline
    t {mustBeNumeric,mustBeReal}
    options.D (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

if options.D > self.K-1
    x_out = zeros(size(t), 'like', t);
    return;
end
x_out = BSpline.evaluateFromPPCoefficients(queryPoints=t, C=self.C, tpp=self.t_pp, D=options.D);
if ~isempty(self.xStd)
    x_out = self.xStd*x_out;
end
if ~isempty(self.xMean) && options.D == 0
    x_out = x_out + self.xMean;
end
