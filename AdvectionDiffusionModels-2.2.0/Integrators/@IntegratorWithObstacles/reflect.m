function [yFinal, yIntersection, didReflect] = reflect(obstacle, yInitial, yFinal)
% Reflect step segments off one polygonal obstacle boundary.
%
% This developer method evaluates the repository's existing crossing-count
% rules for segments that start outside, on, or across the polygon
% boundary. At most one reflection is applied per call, so callers repeat
% the method until no more reflections are needed.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: [yFinal, yIntersection, didReflect] = reflect(obstacle,yInitial,yFinal)
% - Parameter obstacle: cached polygon geometry for one reflecting obstacle
% - Parameter yInitial: initial segment points stored as `nParticles x 2`
% - Parameter yFinal: proposed final segment points stored as `nParticles x 2`
% - Returns yFinal: updated final points after at most one reflection
% - Returns yIntersection: point of reflection or zero row when no reflection occurs
% - Returns didReflect: logical column vector marking reflected segments
nPoints = size(yInitial, 1);

didReflect = false(nPoints, 1);
yIntersection = zeros(size(yFinal));

isYFinalInterior = IntegratorWithObstacles.isInterior(obstacle, yFinal(:,1), yFinal(:,2));
isYInitialOnBorder = false(nPoints, 1);
numCrossings = zeros(nPoints, 1);

firstNonzeroCrossing = zeros(size(yFinal));
distanceToFirstNonzeroCrossing = nan(nPoints, 1);
firstNonzeroCrossingEdge = nan(nPoints, 1);
zeroCrossingEdge = nan(nPoints, 1);

zeroThreshold = 0;
intersectionBackoff = 1e-6;
for iEdge = 1:(length(obstacle.Vertices) - 1)
    yIntersect = nan(size(yInitial));
    r2 = nan(nPoints, 1);
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

    doesIntersect = t >= 0 & t <= 1 & u >= 0 & u <= 1 & denom ~= 0;
    if any(doesIntersect)
        dx = -t(doesIntersect).*x1x2(doesIntersect);
        dy = -t(doesIntersect).*y1y2(doesIntersect);
        dr = sqrt(dx.^2 + dy.^2);
        nonzero = dr > 0;
        dx(nonzero) = max(dr(nonzero) - intersectionBackoff, 0).*dx(nonzero)./dr(nonzero);
        dy(nonzero) = max(dr(nonzero) - intersectionBackoff, 0).*dy(nonzero)./dr(nonzero);
        r2(doesIntersect) = dx.^2 + dy.^2;
        yIntersect(doesIntersect, 1) = yInitial(doesIntersect, 1) + dx;
        yIntersect(doesIntersect, 2) = yInitial(doesIntersect, 2) + dy;
    end

    numCrossings = numCrossings + doesIntersect;
    zeroCrossing = doesIntersect & r2 <= zeroThreshold;
    isYInitialOnBorder = isYInitialOnBorder | zeroCrossing;
    newZeroCrossing = zeroCrossing & isnan(zeroCrossingEdge);
    if any(newZeroCrossing)
        zeroCrossingEdge(newZeroCrossing) = iEdge;
    end

    newIntersection = doesIntersect & r2 > zeroThreshold & (isnan(distanceToFirstNonzeroCrossing) | r2 < distanceToFirstNonzeroCrossing);
    if any(newIntersection)
        distanceToFirstNonzeroCrossing(newIntersection) = r2(newIntersection);
        firstNonzeroCrossing(newIntersection, :) = yIntersect(newIntersection, :);
        firstNonzeroCrossingEdge(newIntersection) = iEdge;
    end
end

reflectAtInitialPoint = isYInitialOnBorder & ((isYFinalInterior & mod(numCrossings, 2) == 1) | (~isYFinalInterior & numCrossings > 0 & mod(numCrossings, 2) == 0));
reflectAtFirstNonzero = (~isYInitialOnBorder & numCrossings > 0) | (isYInitialOnBorder & ((isYFinalInterior & mod(numCrossings, 2) == 0 & ~isnan(firstNonzeroCrossingEdge)) | (~isYFinalInterior & numCrossings > 1 & mod(numCrossings, 2) == 1)));

didReflect = reflectAtInitialPoint | reflectAtFirstNonzero;
yIntersection(reflectAtFirstNonzero, :) = firstNonzeroCrossing(reflectAtFirstNonzero, :);
yIntersection(reflectAtInitialPoint, :) = yInitial(reflectAtInitialPoint, :);

reflectionEdge = nan(nPoints, 1);
reflectionEdge(reflectAtInitialPoint) = zeroCrossingEdge(reflectAtInitialPoint);
reflectionEdge(reflectAtFirstNonzero) = firstNonzeroCrossingEdge(reflectAtFirstNonzero);

reflecting = didReflect & ~isnan(reflectionEdge);
if any(reflecting)
    px = yFinal(reflecting, 1) - yIntersection(reflecting, 1);
    py = yFinal(reflecting, 2) - yIntersection(reflecting, 2);

    dx = obstacle.dx(reflectionEdge(reflecting));
    dy = obstacle.dy(reflectionEdge(reflecting));

    a = dx.^2 - dy.^2;
    b = 2*dx.*dy;
    d = dx.^2 + dy.^2;

    yFinal(reflecting, 1) = (a.*px + b.*py)./d + yIntersection(reflecting, 1);
    yFinal(reflecting, 2) = (b.*px - a.*py)./d + yIntersection(reflecting, 2);
end
end
