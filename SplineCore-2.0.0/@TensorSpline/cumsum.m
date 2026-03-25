function intspline = cumsum(self, dim)
% Return the indefinite integral along one tensor dimension.
%
% The integral is zero at the lower bound of the selected
% dimension while holding all other coordinates fixed.
%
% This applies the one-dimensional B-spline integration formula along the
% chosen tensor dimension and leaves the other dimensions unchanged.
%
% ```matlab
% F = cumsum(spline);
% Fy = cumsum(spline, 2);
% ```
%
% - Topic: Transform the spline
% - Declaration: intspline = cumsum(self,dim)
% - Parameter self: TensorSpline instance
% - Parameter dim: tensor dimension to integrate along
% - Returns intspline: TensorSpline representing the integral
arguments
    self (1,1) TensorSpline
    dim (1,1) double {mustBeInteger,mustBePositive} = 1
end

if dim > self.numDimensions
    error('TensorSpline:InvalidDimension',  'dim must not exceed the number of spline dimensions.');
end

[xi, tKnotDim, Kdim] = TensorSpline.integrateAlongDimension(  self.xi, self.tKnot_{dim}, self.K(dim), dim, self.xMean, self.xStd);

K = self.K;
K(dim) = Kdim;
tKnot = self.tKnot_;
tKnot{dim} = tKnotDim;

intspline = TensorSpline(S=K-1, knotPoints=tKnot, xi=xi);
