function regenerateObstacles(self)
% Rebuild the chamber walls and escape opening.
%
% The geometry consists of five rectangular wall sections whose union
% leaves an opening of width `W` in the left wall.
%
% - Topic: Rebuild model geometry
% - Declaration: regenerateObstacles(self)
a = self.L;
d = self.delta;
topWindow = (a - self.W) / 2 + self.W;
bottomWindow = (a - self.W) / 2;

bottom = polyshape([-d -d a + d a + d], [0 -d -d 0]);
right = polyshape([a a a + d a + d], [a + d -d -d a + d]);
top = polyshape([-d -d a + d a + d], [a + d a a a + d]);
leftTop = polyshape([-d -d 0 0], [a + d topWindow topWindow a + d]);
leftBottom = polyshape([-d -d 0 0], [bottomWindow -d -d bottomWindow]);

self.obstacles = union([bottom; right; top; leftTop; leftBottom]);
