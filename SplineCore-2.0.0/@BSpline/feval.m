function values = feval(spline,t,options)
% Evaluate a B-spline at the supplied points.
%
% This is a thin wrapper around `valueAtPoints(...)`.
%
% ```matlab
% values = feval(spline, tQuery);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: values = feval(spline,t,options)
% - Parameter spline: BSpline instance
% - Parameter t: evaluation points
% - Parameter options.D: derivative order to evaluate
% - Returns values: spline values with the same shape as the query input
arguments
    spline (1,1) BSpline
    t {mustBeNumeric,mustBeReal}
    options.D (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

values = spline.valueAtPoints(t, D=options.D);
