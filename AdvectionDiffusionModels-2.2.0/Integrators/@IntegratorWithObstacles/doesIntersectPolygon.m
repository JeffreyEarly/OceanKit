function doesIntersect = doesIntersectPolygon(obstacle, yInitial, yFinal)
% Test whether a segment intersects any edge of a polygon.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: doesIntersect = doesIntersectPolygon(obstacle,yInitial,yFinal)
% - Parameter obstacle: cached polygon geometry
% - Parameter yInitial: initial segment points stored as `nParticles x 2`
% - Parameter yFinal: final segment points stored as `nParticles x 2`
% - Returns doesIntersect: logical column vector marking segments that intersect the polygon
doesIntersect = false(size(yInitial, 1), 1);
for iEdge = 1:(length(obstacle.Vertices) - 1)
    x1x2 = yInitial(:,1) - yFinal(:,1);
    y3y4 = -obstacle.dy(iEdge);
    y1y2 = yInitial(:,2) - yFinal(:,2);
    x3x4 = -obstacle.dx(iEdge);
    x1x3 = yInitial(:,1) - obstacle.Vertices(iEdge, 1);
    y1y3 = yInitial(:,2) - obstacle.Vertices(iEdge, 2);

    denom = x1x2*y3y4 - y1y2*x3x4;
    a = x1x3.*y3y4 - y1y3.*x3x4;
    b = x1x2.*y1y3 - y1y2.*x1x3;

    t = a./denom;
    u = -b./denom;
    doesIntersect = doesIntersect | (t >= 0 & t <= 1 & u >= 0 & u <= 1 & denom ~= 0);
end
end
