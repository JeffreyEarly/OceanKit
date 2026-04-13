function values = valueAtPoints(self, X, options)
% Evaluate the tensor spline or a mixed partial derivative.
%
% This is the primary explicit evaluation method. Supply one
% matching-size query array per tensor dimension. Paired column
% vectors give pointwise queries, while matching ndgrid arrays give
% gridded evaluation over a tensor-product lattice.
%
% For derivative order vector `D = [D_1 ... D_d]`, this evaluates
%
% $$
% \partial^{D} f(x_1,\ldots,x_d) =
% x_{\mathrm{Std}}
% \sum_{j_1,\ldots,j_d} \xi_{j_1,\ldots,j_d}
% \prod_{k=1}^{d} B_{j_k,S_k}^{(D_k)}(x_k;\tau_k),
% $$
%
% with `xMean` added back only when all entries of `D` are zero.
%
% ```matlab
% values = spline.valueAtPoints(xq, yq);
% dFdx = spline.valueAtPoints(xq, yq, D=[1 0]);
% [Xq,Yq] = ndgrid(linspace(-1,1,40), linspace(0,2,50));
% F = spline.valueAtPoints(Xq, Yq);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: values = valueAtPoints(self,X1,...,Xn,options)
% - Parameter self: TensorSpline instance
% - Parameter X1,...,Xn: matching-size query locations as one array per dimension
% - Parameter options.D: derivative order per dimension
% - Returns values: spline values reshaped to match the query input
arguments
    self (1,1) TensorSpline
end
arguments (Repeating)
    X {mustBeNumeric,mustBeReal}
end
arguments
    options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 0
end

queryInputs = X;
if isempty(queryInputs)
    error('TensorSpline:NotEnoughInputs', 'Specify one query input per spline dimension.');
end

if numel(queryInputs) ~= self.numDimensions
    error('TensorSpline:InvalidEvaluationInput', 'Use spline.valueAtPoints(X1,...,Xn) or spline.valueAtPoints(X1,...,Xn,D=...).');
end
derivativeOrders = TensorSpline.normalizeDerivativeOrders(options.D, self.numDimensions);

validateattributes(queryInputs{1}, {'numeric'}, {'real'});
outputSize = size(queryInputs{1});
numPoints = numel(queryInputs{1});
pointMatrix = zeros(numPoints, self.numDimensions);
for iDim = 1:self.numDimensions
    validateattributes(queryInputs{iDim}, {'numeric'}, {'real'});
    if ~isequal(size(queryInputs{iDim}), outputSize)
        error('TensorSpline:InvalidQueryArrays', 'All query inputs must have the same size.');
    end
    pointMatrix(:,iDim) = queryInputs{iDim}(:);
end

if any(derivativeOrders > self.K - 1)
    values = zeros(outputSize, 'like', pointMatrix);
    return;
end

basisMatrix = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=self.tKnot_, S=self.S, D=derivativeOrders);
values = basisMatrix * self.xi(:);

if ~isempty(self.xStd)
    values = self.xStd * values;
end

if ~isempty(self.xMean) && all(derivativeOrders == 0)
    values = values + self.xMean;
end

values = reshape(values, outputSize);
