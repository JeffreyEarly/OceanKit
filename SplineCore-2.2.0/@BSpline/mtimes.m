function f = mtimes(f,g)
% Multiply a spline output by a scalar.
%
% This rescales the spline output without refitting the spline
% coefficients.
%
% ```matlab
% scaledSpline = 2.5 * spline;
% ```
%
% - Topic: Transform the spline
% - Declaration: f = mtimes(f,g)
% - Parameter f: BSpline instance or scalar
% - Parameter g: scalar or BSpline instance
% - Returns f: transformed BSpline or empty when multiplying by []
arguments
    f
    g
end

if ~isa(f,'BSpline')
    [f, g] = deal(g, f);
end

if ~isa(f,'BSpline')
    error('BSpline:mtimes:UnsupportedOperand', 'One operand must be a BSpline.');
elseif isempty(g)
    f = [];
elseif isnumeric(g) && isscalar(g)
    f = f.affineOutputTransform(g, 0);
else
    error('BSpline:mtimes:UnsupportedOperand', 'Only scalar numeric multiplication is supported.');
end
