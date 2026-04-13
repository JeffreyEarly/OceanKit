function [doesIntersect, yIntersection, r2] = doesIntersectSegment(yInitial, yFinal, vInitial, vFinal, doesIntersect, yIntersection, r2)
% Compute vectorized segment-segment intersections.
%
% This developer helper retains the legacy argument pattern that allows
% MATLAB to reuse the output arrays `doesIntersect`, `yIntersection`, and
% `r2` during repeated intersection queries.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: [doesIntersect, yIntersection, r2] = doesIntersectSegment(yInitial,yFinal,vInitial,vFinal,doesIntersect,yIntersection,r2)
% - Parameter yInitial: initial segment points stored as `nParticles x 2`
% - Parameter yFinal: final segment points stored as `nParticles x 2`
% - Parameter vInitial: initial polygon-edge vertex stored as `1 x 2`
% - Parameter vFinal: final polygon-edge vertex stored as `1 x 2`
% - Parameter doesIntersect: reusable logical output buffer
% - Parameter yIntersection: reusable intersection-point output buffer
% - Parameter r2: reusable squared-distance output buffer
% - Returns doesIntersect: logical column vector marking intersecting segments
% - Returns yIntersection: intersection points backed off slightly from the segment start
% - Returns r2: squared distance from `yInitial` to `yIntersection`
x1x2 = yInitial(:,1) - yFinal(:,1);
y3y4 = vInitial(2) - vFinal(2);
y1y2 = yInitial(:,2) - yFinal(:,2);
x3x4 = vInitial(1) - vFinal(1);
x1x3 = yInitial(:,1) - vInitial(:,1);
y1y3 = yInitial(:,2) - vInitial(:,2);

denom = x1x2*y3y4 - y1y2*x3x4;
a = x1x3.*y3y4 - y1y3.*x3x4;
b = x1x2.*y1y3 - y1y2.*x1x3;

t = a./denom;
u = -b./denom;

doesIntersect = t >= 0 & t <= 1 & u >= 0 & u <= 1 & denom ~= 0;
if any(doesIntersect)
    dx = -0.999*t(doesIntersect).*x1x2(doesIntersect);
    dy = -0.999*t(doesIntersect).*y1y2(doesIntersect);
    r2(doesIntersect) = dx.^2 + dy.^2;
    yIntersection(doesIntersect, 1) = yInitial(doesIntersect, 1) + dx;
    yIntersection(doesIntersect, 2) = yInitial(doesIntersect, 2) + dy;
end
end
