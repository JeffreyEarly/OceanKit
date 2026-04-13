function values = roots(spline)
% Return real roots of a spline within its domain.
%
% Use this to locate zero crossings of the spline over its support.
%
% The implementation works interval by interval on the cached
% piecewise-polynomial coefficients, then keeps only real roots that lie in
% the corresponding interval.
%
% ```matlab
% tZero = roots(spline);
% ```
%
% - Topic: Transform the spline
% - Declaration: values = roots(spline)
% - Parameter spline: BSpline instance
% - Returns values: sorted real roots in the spline domain
arguments
    spline (1,1) BSpline
end
values = [];
scale = factorial((spline.K-1):-1:0);
C = spline.xStd*spline.C;
C(:,end) = C(:,end) + spline.xMean;
t_pp = spline.t_pp;

for iBin=1:size(spline.C,1)
    localRoots = roots(C(iBin,:)./scale);
    intervalWidth = t_pp(iBin+1) - t_pp(iBin);
    isValidRoot = imag(localRoots) == 0 & real(localRoots) >= 0 & real(localRoots) <= intervalWidth;
    if any(isValidRoot)
        values = cat(1,values,real(localRoots(isValidRoot)) + t_pp(iBin));
    end
end

values = sort(values);
