function poweredSpline = power(self, exponent)
% Raise tensor-spline values to a positive scalar power by refitting support values.
%
% This is useful for simple nonlinear transforms of a tensor spline when
% an exact spline-space representation is not available.
%
% The implementation samples the spline on representative tensor support
% points, applies the power to those sampled values, chooses a new tensor
% basis, and refits the transformed samples. This is an approximation, not
% exact spline-space algebra.
%
% ```matlab
% squaredSpline = spline.^2;
% amplitudeSpline = spline.^(1/2);
% ```
%
% - Topic: Transform the spline
% - Declaration: poweredSpline = power(self,exponent)
% - Parameter self: TensorSpline instance
% - Parameter exponent: positive scalar exponent
% - Returns poweredSpline: TensorSpline approximating spline.^exponent
arguments
    self (1,1) TensorSpline
    exponent (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
end

if exponent == 1
    poweredSpline = self;
    return;
end

[supportPoints, supportVectors] = TensorSpline.pointsOfSupportFromKnotPoints(self.knotPoints, S=self.S);
basisMatrix = TensorSpline.matrixForPointMatrix(supportPoints, knotPoints=self.knotPoints, S=self.S);
values = basisMatrix * self.xi(:);
if ~isempty(self.xStd)
    values = self.xStd * values;
end
if ~isempty(self.xMean)
    values = values + self.xMean;
end
values(abs(values) < 2*eps) = 0;

poweredValues = values.^exponent;
maxSupportedOrder = cellfun(@numel, supportVectors);
poweredK = min(maxSupportedOrder, max(1, ceil(exponent * self.K)));
poweredTKnot = cell(1, self.numDimensions);
for iDim = 1:self.numDimensions
    poweredTKnot{iDim} = BSpline.knotPointsForDataPoints(supportVectors{iDim}, S=poweredK(iDim)-1);
end

X = TensorSpline.matrixForPointMatrix(supportPoints, knotPoints=poweredTKnot, S=poweredK - 1);
xi = X \ poweredValues(:);
poweredSpline = TensorSpline(S=poweredK-1, knotPoints=poweredTKnot, xi=xi);
