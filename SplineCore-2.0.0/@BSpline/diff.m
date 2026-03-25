function dspline = diff(spline,n)
% Differentiate a B-spline representation.
%
% Use this when you want a new spline object representing the derivative,
% rather than only evaluating derivatives at a set of points.
%
% The implementation applies the standard B-spline coefficient-difference
% recursion and trims one knot from each end for each derivative taken.
%
% ```matlab
% dspline = diff(spline);
% curvatureSpline = diff(spline, 2);
% ```
%
% - Topic: Transform the spline
% - Declaration: dspline = diff(spline,n)
% - Parameter spline: BSpline instance to differentiate
% - Parameter n: derivative order
% - Returns dspline: BSpline representing the nth derivative
arguments
    spline (1,1) BSpline
    n (1,1) double {mustBeInteger,mustBeNonnegative} = 1
end

if n == 0
    dspline = spline;
elseif n >= spline.K
    dspline = BSpline(S=0, knotPoints=reshape(spline.domain,[],1), xi=0);
else
    D = n;
    xi = spline.xi;
    K = spline.K;
    knotPoints = spline.knotPoints;
    M = length(xi);
    
    alpha = zeros(length(xi),D+1);
    alpha(:,1) = xi; % first column is the existing coefficients
    
    for d=1:D
        dm = diff(alpha(:,d));
        dt = (knotPoints(1+K-d:M+K-d)-knotPoints(1:M))/(K-d);
        alpha(1:end-d,d+1) = dm./dt(d+1:end);
    end
    
    dspline = BSpline(S=K-D-1, knotPoints=knotPoints((1+D):(end-D)), xi=alpha(1:end-D,D+1));
    dspline.xStd = spline.xStd;
end
