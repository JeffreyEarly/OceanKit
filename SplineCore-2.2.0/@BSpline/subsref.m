function varargout = subsref(self, index)
% Evaluate the spline with function-call syntax or defer to built-in indexing.
%
% Parentheses indexing `spline(t)` is redirected to
% `valueAtPoints(...)`, while dot indexing behaves like the default
% MATLAB handle-class implementation.
%
% Use `spline(t)` for values. For derivatives, use
% `valueAtPoints(t, D=...)`.
%
% ```matlab
% x = spline(tQuery);
% dxdt = spline.valueAtPoints(tQuery, D=1);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: varargout = subsref(self,index)
% - Parameter self: BSpline instance
% - Parameter index: MATLAB subscript structure
% - Returns varargout: indexed property access or spline values
idx = index(1).subs;
switch index(1).type
    case '()'
        if numel(idx) ~= 1
            error('BSpline:InvalidEvaluationInput', 'Use spline(t) for values and spline.valueAtPoints(t, D=...) for derivatives.');
        end
        varargout{1} = self.valueAtPoints(idx{1});
    case '.'
        [varargout{1:nargout}] = builtin('subsref',self,index);
    case '{}'
        error('The BSpline class does not know what to do with {}.');
    otherwise
        error('Unexpected syntax');
end
