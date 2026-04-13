function isLeft = isLeft(obstacle, iEdge, x, y)
% Compute the signed left-of-edge test for one polygon edge.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: isLeft = isLeft(obstacle,iEdge,x,y)
% - Parameter obstacle: cached polygon geometry
% - Parameter iEdge: scalar polygon-edge index
% - Parameter x: x-coordinate values to test
% - Parameter y: y-coordinate values to test
% - Returns isLeft: signed edge-orientation value used by the winding-number test
isLeft = obstacle.dx(iEdge)*(y - obstacle.Vertices(iEdge, 2)) - (x - obstacle.Vertices(iEdge, 1))*obstacle.dy(iEdge);
end
