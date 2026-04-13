function f = mtimes(f, g)
% Multiply tensor-spline outputs by a scalar.
%
% ```matlab
% scaledSpline = 2.5 * spline;
% ```
%
% - Topic: Transform the spline
% - Declaration: f = mtimes(f,g)
% - Parameter f: TensorSpline instance or scalar
% - Parameter g: scalar or TensorSpline instance
% - Returns f: transformed TensorSpline or empty when multiplying by []
arguments
    f
    g
end

if ~isa(f, 'TensorSpline')
    [f, g] = deal(g, f);
end

if ~isa(f, 'TensorSpline')
    error('TensorSpline:mtimes:UnsupportedOperand',  'One operand must be a TensorSpline.');
elseif isempty(g)
    f = [];
elseif isnumeric(g) && isscalar(g)
    f = TensorSpline(S=f.S, knotAxes=f.knotAxes, xi=f.xi, xMean=g*f.xMean, xStd=g*f.xStd);
else
    error('TensorSpline:mtimes:UnsupportedOperand',  'Only scalar numeric multiplication is supported.');
end
