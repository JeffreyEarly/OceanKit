function dspline = diff(self, derivativeOrders)
% Return a tensor spline representing mixed partial derivatives.
%
% Use a scalar derivative order in 1-D or a derivative-order
% vector with one entry per dimension.
%
% The implementation differentiates the tensor-product coefficients along
% each requested dimension and reduces the degree in those dimensions by
% the corresponding derivative orders.
%
% ```matlab
% dspline = diff(spline);
% dFdx = diff(spline, [1 0]);
% ```
%
% - Topic: Transform the spline
% - Declaration: dspline = diff(self,derivativeOrders)
% - Parameter self: TensorSpline instance
% - Parameter derivativeOrders: derivative order per dimension
% - Returns dspline: TensorSpline representing the derivative
arguments
    self (1,1) TensorSpline
    derivativeOrders {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 1
end

derivativeOrders = TensorSpline.normalizeDerivativeOrders(derivativeOrders, self.numDimensions);
if all(derivativeOrders == 0)
    dspline = self;
    return;
end

if any(derivativeOrders > self.K - 1)
    dspline = TensorSpline.zeroSplineForDomain(self.domain, self.numDimensions, xStd=self.xStd);
    return;
end

xi = self.xi;
K = self.K;
tKnot = self.tKnot_;
for iDim = 1:self.numDimensions
    if derivativeOrders(iDim) == 0
        continue;
    end

    [xi, tKnot{iDim}, K(iDim)] = TensorSpline.differentiateAlongDimension(  xi, tKnot{iDim}, K(iDim), derivativeOrders(iDim), iDim);
end

dspline = TensorSpline(S=K-1, knotPoints=tKnot, xi=xi, xMean=0, xStd=self.xStd);
