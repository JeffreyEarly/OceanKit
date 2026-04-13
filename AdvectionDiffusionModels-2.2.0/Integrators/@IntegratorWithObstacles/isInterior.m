function isInterior = isInterior(obstacle, x, y)
% Determine whether points lie inside one polygonal obstacle.
%
% This method evaluates the winding-number test using the cached polygon
% edge geometry instead of MATLAB's general-purpose polygon query.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: isInterior = isInterior(obstacle,x,y)
% - Parameter obstacle: cached polygon geometry
% - Parameter x: x-coordinate query values
% - Parameter y: y-coordinate query values
% - Returns isInterior: logical array marking points inside the polygon
windingNumber = zeros(size(x));
for iEdge = 1:(length(obstacle.Vertices) - 1)
    isBelow = obstacle.Vertices(iEdge, 2) <= y;
    upwardCrossing = isBelow & obstacle.Vertices(iEdge + 1, 2) > y;
    if any(upwardCrossing)
        windingNumber(upwardCrossing) = windingNumber(upwardCrossing) + (IntegratorWithObstacles.isLeft(obstacle, iEdge, x(upwardCrossing), y(upwardCrossing)) > 0);
    end

    downwardCrossing = ~isBelow & obstacle.Vertices(iEdge + 1, 2) <= y;
    if any(downwardCrossing)
        windingNumber(downwardCrossing) = windingNumber(downwardCrossing) - (IntegratorWithObstacles.isLeft(obstacle, iEdge, x(downwardCrossing), y(downwardCrossing)) < 0);
    end
end

isInterior = abs(windingNumber) > 0;
end
