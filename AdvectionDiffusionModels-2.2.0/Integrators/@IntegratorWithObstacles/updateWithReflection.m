function [y, dy, didReflect] = updateWithReflection(self, obstacle, y, yWrapped, dy)
% Reflect proposed step segments that enter one obstacle.
%
% The candidate segment from `yWrapped - dy` to `yWrapped` is first checked
% against the obstacle bounding box. Any segment that may intersect or end
% inside the obstacle is passed to `reflect` to compute a single
% reflection.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: [y, dy, didReflect] = updateWithReflection(self,obstacle,y,yWrapped,dy)
% - Parameter obstacle: cached polygon geometry for one reflecting obstacle
% - Parameter y: accepted state array to be updated for reflected particles
% - Parameter yWrapped: wrapped candidate endpoints prior to reflection
% - Parameter dy: candidate increments measured from the unwrapped state
% - Returns y: updated accepted states for reflected particles
% - Returns dy: updated increments after one reflection
% - Returns didReflect: logical column vector marking which particles reflected
yWrappedInitial = yWrapped - dy;

mayIntersect = IntegratorWithObstacles.doesIntersectPolygon(obstacle.boundingBox, yWrappedInitial, yWrapped);
mayIntersect = mayIntersect | IntegratorWithObstacles.isInterior(obstacle.boundingBox, yWrapped(:,1), yWrapped(:,2));

didReflect = false(size(y, 1), 1);
if ~any(mayIntersect)
    return
end

[yFixed, yIntersection, didReflectSubset] = IntegratorWithObstacles.reflect(obstacle, yWrappedInitial(mayIntersect, :), yWrapped(mayIntersect, :));
didReflect(mayIntersect) = didReflectSubset;

if ~any(didReflectSubset)
    return
end

allIndices = 1:size(yWrapped, 1);
subIndices = allIndices(mayIntersect);
subIndices = subIndices(didReflectSubset);

dIntersection = yIntersection(didReflectSubset, :) - yWrappedInitial(subIndices, :);
dReflection = yFixed(didReflectSubset, :) - yIntersection(didReflectSubset, :);

y(subIndices, :) = y(subIndices, :) + dIntersection;
dy(subIndices, :) = dReflection;
end
