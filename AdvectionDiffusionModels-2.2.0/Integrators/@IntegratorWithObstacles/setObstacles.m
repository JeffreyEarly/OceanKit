function setObstacles(self, obstacles)
% Cache polygon geometry used by the obstacle reflection algorithm.
%
% Each obstacle is stored together with a closed vertex sequence, segment
% differences, and a rectangular bounding box used for the fast pre-checks
% in `updateWithReflection`.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: setObstacles(self,obstacles)
% - Parameter obstacles: array of reflecting `polyshape` obstacles
self.obstacleData = repmat(emptyObstacleData(), length(obstacles), 1);

for iObstacle = 1:length(obstacles)
    obstacleData = emptyObstacleData();
    obstacleData.Vertices = obstacles(iObstacle).Vertices;
    obstacleData.Vertices(end + 1, :) = obstacles(iObstacle).Vertices(1, 1:2);
    obstacleData.polyshape = obstacles(iObstacle);
    obstacleData.dx = obstacleData.Vertices(2:end, 1) - obstacleData.Vertices(1:end-1, 1);
    obstacleData.dy = obstacleData.Vertices(2:end, 2) - obstacleData.Vertices(1:end-1, 2);

    [xbb, ybb] = boundingbox(obstacles(iObstacle));
    boundingBox.polyshape = polyshape([min(xbb); min(xbb); max(xbb); max(xbb)], [max(ybb); min(ybb); min(ybb); max(ybb)]);
    boundingBox.Vertices = boundingBox.polyshape.Vertices;
    boundingBox.Vertices(end + 1, :) = boundingBox.Vertices(1, :);
    boundingBox.dx = boundingBox.Vertices(2:end, 1) - boundingBox.Vertices(1:end-1, 1);
    boundingBox.dy = boundingBox.Vertices(2:end, 2) - boundingBox.Vertices(1:end-1, 2);

    obstacleData.boundingBox = boundingBox;
    self.obstacleData(iObstacle) = obstacleData;
end
end

function obstacleData = emptyObstacleData()
boundingBox = struct( ...
    "polyshape", polyshape.empty(0,1), ...
    "Vertices", zeros(0, 2), ...
    "dx", zeros(0, 1), ...
    "dy", zeros(0, 1));

obstacleData = struct( ...
    "Vertices", zeros(0, 2), ...
    "polyshape", polyshape.empty(0,1), ...
    "dx", zeros(0, 1), ...
    "dy", zeros(0, 1), ...
    "boundingBox", boundingBox);
end
