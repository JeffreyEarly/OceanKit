function poweredSpline = power(spline,exponent)
% Raise spline values to a real scalar power by refitting support values.
%
% This is useful for simple nonlinear transforms of a spline when an exact
% spline-space representation is not available.
%
% The implementation samples the spline at representative support points,
% forms `values.^exponent`, and fits a new spline to those transformed
% samples. When the transformed samples are nonnegative up to tolerance, it
% uses a positive constrained fit; otherwise it uses an unconstrained fit.
%
% ```matlab
% squaredSpline = spline.^2;
% rootSpline = spline.^(1/2);
% ```
%
% - Topic: Transform the spline
% - Declaration: poweredSpline = power(spline,exponent)
% - Parameter spline: BSpline instance
% - Parameter exponent: scalar exponent
% - Returns poweredSpline: BSpline approximating spline.^exponent
arguments
    spline (1,1) BSpline
    exponent (1,1) double {mustBeReal,mustBeFinite}
end
if exponent == 1
    poweredSpline = spline;
    return;
end

supportPoints = BSpline.pointsOfSupportFromKnotPoints(spline.knotPoints, S=spline.S);
values = spline.valueAtPoints(supportPoints);
values(abs(values) < 2*eps) = 0;

poweredK = ceil(exponent*spline.K);
knotPoints = BSpline.knotPointsForDataPoints(supportPoints, S=poweredK-1);
poweredValues = values.^exponent;
tolerance = 10*eps(max(1, max(abs(poweredValues))));
if all(isfinite(poweredValues)) && all(poweredValues >= -tolerance)
    poweredValues = max(poweredValues, 0);
    fittedSpline = ConstrainedSpline.fromGriddedValues(supportPoints, poweredValues, S=poweredK-1, knotPoints=knotPoints, constraints=GlobalConstraint.positive());
    poweredSpline = BSpline(S=poweredK-1, knotPoints=fittedSpline.knotPoints, xi=fittedSpline.xi(:));
else
    X = BSpline.matrixForDataPoints(supportPoints, knotPoints=knotPoints, S=poweredK-1);
    xi = X\poweredValues;
    poweredSpline = BSpline(S=poweredK-1, knotPoints=knotPoints, xi=xi);
end
