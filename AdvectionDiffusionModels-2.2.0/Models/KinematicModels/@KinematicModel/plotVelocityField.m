function varargout = plotVelocityField(self, options)
% Plot the model velocity field with quivers.
%
% The plotting grid uses `2*numPoints` points in x and `numPoints` points
% in y to preserve the typical wide aspect ratio.
%
% - Topic: Plot model diagnostics
% - Declaration: plotVelocityField(self,options)
% - Parameter options.t: scalar plotting time in seconds
% - Parameter options.quiverScale: scale factor forwarded to `quiver`
% - Parameter options.numPoints: number of points in the y-direction
% - Returns X: optional x-grid in meters when two outputs are requested
% - Returns Y: optional y-grid in meters when two outputs are requested
arguments
    self (1,1) KinematicModel
    options.t (1,1) double = 0
    options.quiverScale (1,1) double = 1
    options.numPoints (1,1) double {mustBePositive, mustBeInteger} = 150
end

xg = linspace(min(self.xVisualLimits), max(self.xVisualLimits), 2 * options.numPoints).';
yg = linspace(min(self.yVisualLimits), max(self.yVisualLimits), options.numPoints).';
[X, Y] = meshgrid(xg, yg);

if ~isempty(self.obstacles)
    x = reshape(X, [], 1);
    y = reshape(Y, [], 1);
    mask = true(size(x));
    for iObstacle = 1:length(self.obstacles)
        mask = mask & ~isinterior(self.obstacles(iObstacle), x, y);
    end
    mask = reshape(mask, size(X));

    quiver(X / self.visualScale, Y / self.visualScale, mask .* self.u(options.t, X, Y), mask .* self.v(options.t, X, Y), options.quiverScale);
    hold on
    plot(scale(self.obstacles, 1 / self.visualScale))
else
    quiver(X / self.visualScale, Y / self.visualScale, self.u(options.t, X, Y), self.v(options.t, X, Y), options.quiverScale)
end

axis equal
xlim([min(self.xVisualLimits) max(self.xVisualLimits)] / self.visualScale)
ylim([min(self.yVisualLimits) max(self.yVisualLimits)] / self.visualScale)
xlabel('km')
ylabel('km')

if nargout == 2
    varargout{1} = X;
    varargout{2} = Y;
end
