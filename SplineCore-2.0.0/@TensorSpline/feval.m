function values = feval(spline, X, options)
% Evaluate a tensor spline at matching-size query arrays.
%
% This is a thin wrapper around `valueAtPoints(...)`. Paired
% column vectors give pointwise queries, while matching `ndgrid`
% arrays give gridded evaluation on a tensor-product lattice.
%
% ```matlab
% values = feval(spline, xq, yq);
% [Xq,Yq] = ndgrid(linspace(-1,1,40), linspace(0,2,50));
% F = feval(spline, Xq, Yq);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: values = feval(spline,X1,...,Xn,options)
% - Parameter spline: TensorSpline instance
% - Parameter X1,...,Xn: matching-size query locations as one array per dimension
% - Parameter options.D: derivative order per dimension
% - Returns values: spline values with the same shape as the query input
arguments
    spline (1,1) TensorSpline
end
arguments (Repeating)
    X {mustBeNumeric,mustBeReal}
end
arguments
    options.D {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 0
end

values = spline.valueAtPoints(X{:}, D=options.D);
