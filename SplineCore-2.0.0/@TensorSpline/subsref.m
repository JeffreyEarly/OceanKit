function varargout = subsref(self, index)
% Evaluate the tensor spline with function-call syntax or defer to built-in indexing.
%
% Function-call syntax is a thin wrapper around `valueAtPoints(...)`.
% Use `spline(X1,...,Xn)` for pointwise values at matching-size
% query arrays. Paired column vectors evaluate paired sample locations,
% while matching `ndgrid` arrays evaluate a tensor-product lattice. For
% derivatives, use `valueAtPoints(X1,...,Xn,D=...)`.
%
% ```matlab
% values = spline(xq, yq);
% [Xq,Yq] = ndgrid(linspace(-1,1,40), linspace(0,2,50));
% F = spline(Xq, Yq);
% dFdx = spline.valueAtPoints(xq, yq, D=[1 0]);
% ```
%
% - Topic: Evaluate the spline
% - Declaration: varargout = subsref(self,index)
% - Parameter self: TensorSpline instance
% - Parameter index: MATLAB subscript structure
% - Returns varargout: indexed property access or spline values
idx = index(1).subs;
switch index(1).type
    case '()'
        if numel(idx) ~= self.numDimensions
            error('TensorSpline:InvalidEvaluationInput', 'Use spline(X1,...,Xn) for values and spline.valueAtPoints(X1,...,Xn,D=...) for derivatives.');
        end
        varargout{1} = self.valueAtPoints(idx{:});
    case '.'
        [varargout{1:nargout}] = builtin('subsref',self,index);
    case '{}'
        error('The TensorSpline class does not know what to do with {}.');
    otherwise
        error('Unexpected syntax');
end
