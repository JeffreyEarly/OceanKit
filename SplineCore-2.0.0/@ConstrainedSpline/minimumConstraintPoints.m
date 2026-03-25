function tc = minimumConstraintPoints(knotPoints, S, T)
% Return a minimal set of one-dimensional locations for universal derivative constraints.
%
% For a terminated spline of order K, this chooses the smallest
% set of 1-D points needed to constrain all segments at
% polynomial degree T.
%
% If `D = S - T`, the returned locations provide enough one-dimensional
% sample points to constrain every piecewise-polynomial segment through
% derivative order `T` without oversampling all knots.
%
% ```matlab
% tc = ConstrainedSpline.minimumConstraintPoints(knotPoints, 3, 0);
% ```
%
% - Topic: Choose constraint locations
% - Declaration: tc = minimumConstraintPoints(knotPoints, S, T)
% - Parameter knotPoints: one-dimensional knot sequence
% - Parameter S: spline degree
% - Parameter T: constrained polynomial degree
% - Returns tc: one-dimensional constraint locations
arguments
    knotPoints (:,1) double {mustBeNumeric,mustBeReal,mustBeFinite}
    S (1,1) double {mustBeNonnegative,mustBeInteger}
    T (1,1) double {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger}
end

K = S + 1;
t = unique(knotPoints);
D = K - 1 - T;
if mod(D, 2) == 0
    ts = t(1) + (t(2)-t(1))/(D/2 + 2) * (1:(D/2 + 1)).';
    te = t(end-1) + (t(end)-t(end-1))/(D/2 + 2) * (1:(D/2 + 1)).';
    ti = t(2:end-2) + diff(t(2:end-1))/2;
    tc = cat(1, ts, ti, te);
else
    ts = t(1) + (t(2)-t(1))/((D-1)/2 + 1) * (0:((D-1)/2)).';
    te = t(end) - (t(end)-t(end-1))/((D-1)/2 + 1) * (0:((D-1)/2)).';
    ti = t(2:end-1);
    tc = cat(1, ts, ti, te);
end
