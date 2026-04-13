function f = plus(f,g)
% Add a scalar offset to a spline output.
%
% This shifts the spline output without changing the knot sequence or
% spline coefficients.
%
% ```matlab
% shiftedSpline = spline + 3;
% ```
%
% - Topic: Transform the spline
% - Declaration: f = plus(f,g)
% - Parameter f: BSpline instance or scalar
% - Parameter g: scalar or BSpline instance
% - Returns f: transformed BSpline or empty when adding []
arguments
    f
    g
end

if ~isa(f,'BSpline')
    [f, g] = deal(g, f);
end

if ~isa(f,'BSpline')
    error('BSpline:plus:UnsupportedOperand', 'One operand must be a BSpline.');
elseif isempty(g)
    f = [];
elseif isnumeric(g) && isscalar(g)
    f = f.affineOutputTransform(1, g);
else
    error('BSpline:plus:UnsupportedOperand', 'Only scalar numeric offsets are supported.');
end
