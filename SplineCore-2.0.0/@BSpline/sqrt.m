function splinesqrt = sqrt(spline)
% Return a spline approximation to the square root of the spline output.
%
% This is a convenience wrapper around `spline.^(1/2)` and is most useful
% when the spline is nonnegative over its domain.
%
% ```matlab
% amplitudeSpline = sqrt(energySpline);
% ```
%
% - Topic: Transform the spline
% - Declaration: splinesqrt = sqrt(spline)
% - Parameter spline: BSpline instance
% - Returns splinesqrt: BSpline approximating sqrt(spline)
arguments
    spline (1,1) BSpline
end
splinesqrt = power(spline,1/2);
