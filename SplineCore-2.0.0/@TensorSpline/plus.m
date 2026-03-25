function f = plus(f, g)
% Add a scalar offset to tensor-spline outputs.
%
% ```matlab
% shiftedSpline = spline + 3;
% ```
%
% - Topic: Transform the spline
% - Declaration: f = plus(f,g)
% - Parameter f: TensorSpline instance or scalar
% - Parameter g: scalar or TensorSpline instance
% - Returns f: transformed TensorSpline or empty when adding []
arguments
    f
    g
end

if ~isa(f, 'TensorSpline')
    [f, g] = deal(g, f);
end

if ~isa(f, 'TensorSpline')
    error('TensorSpline:plus:UnsupportedOperand',  'One operand must be a TensorSpline.');
elseif isempty(g)
    f = [];
elseif isnumeric(g) && isscalar(g)
    f = TensorSpline(S=f.S, knotPoints=f.knotPoints, xi=f.xi, xMean=f.xMean + g, xStd=f.xStd);
else
    error('TensorSpline:plus:UnsupportedOperand',  'Only scalar numeric offsets are supported.');
end
