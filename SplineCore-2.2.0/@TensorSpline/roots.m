function values = roots(self)
% Return real roots of a one-dimensional tensor spline within its domain.
%
% Use this to locate zero crossings of a one-dimensional tensor spline
% over its support.
%
% In one dimension this simply delegates to the underlying `BSpline`
% piecewise-polynomial root finder after transferring the degree, knots,
% coefficients, and affine output normalization.
%
% ```matlab
% tZero = roots(spline);
% ```
%
% - Topic: Transform the spline
% - Declaration: values = roots(self)
% - Parameter self: TensorSpline instance
% - Returns values: sorted real roots in the spline domain
arguments
    self (1,1) TensorSpline
end

if self.numDimensions ~= 1
    error('TensorSpline:roots:UnsupportedDimension',  'roots is only defined for one-dimensional TensorSpline objects.');
end

spline1D = BSpline(S=self.S, knotPoints=self.knotPoints, xi=self.xi(:), xMean=self.xMean, xStd=self.xStd);
values = roots(spline1D);
