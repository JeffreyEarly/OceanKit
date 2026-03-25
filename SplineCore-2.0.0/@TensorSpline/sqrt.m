function splinesqrt = sqrt(self)
% Return a tensor spline approximation to the square root of the spline output.
%
% This is a convenience wrapper around `spline.^(1/2)` and is most useful
% when the spline is nonnegative over its domain.
%
% ```matlab
% amplitudeSpline = sqrt(energySpline);
% ```
%
% - Topic: Transform the spline
% - Declaration: splinesqrt = sqrt(self)
% - Parameter self: TensorSpline instance
% - Returns splinesqrt: TensorSpline approximating sqrt(spline)
arguments
    self (1,1) TensorSpline
end

splinesqrt = power(self, 1/2);
