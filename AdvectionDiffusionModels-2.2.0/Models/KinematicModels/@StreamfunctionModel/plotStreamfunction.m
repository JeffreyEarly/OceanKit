function varargout = plotStreamfunction(self, t)
% Plot streamfunction contours on the current axes.
%
% The contour grid spans `xVisualLimits` and `yVisualLimits` and excludes
% obstacle interiors when polygonal obstacles are present.
%
% - Topic: Plot streamfunctions
% - Declaration: plotStreamfunction(self,t)
% - Parameter t: scalar plotting time in seconds; the default is `0`
% - Returns X: optional x-grid in meters when two outputs are requested
% - Returns Y: optional y-grid in meters when two outputs are requested
arguments
    self (1,1) StreamfunctionModel
    t (1,1) double = 0
end

numPoints = 150;
xg = linspace(min(self.xVisualLimits), max(self.xVisualLimits), 2 * numPoints).';
yg = linspace(min(self.yVisualLimits), max(self.yVisualLimits), numPoints).';
[Xf, Yf] = meshgrid(xg, yg);

if ~isempty(self.obstacles)
    x = reshape(Xf, [], 1);
    y = reshape(Yf, [], 1);
    mask = ~isinterior(self.obstacles, x, y);
    mask = reshape(mask, size(Xf));
    psiValues = self.psi(t, Xf(mask), Yf(mask));
    levels = linspace(min(psiValues), max(psiValues), 20);

    contour(Xf / self.visualScale, Yf / self.visualScale, mask .* self.psi(t, Xf, Yf), levels)
    hold on
    plot(scale(self.obstacles, 1 / self.visualScale))
else
    contour(Xf / self.visualScale, Yf / self.visualScale, self.psi(t, Xf, Yf))
end

axis equal
xlim([min(self.xVisualLimits) max(self.xVisualLimits)] / self.visualScale)
ylim([min(self.yVisualLimits) max(self.yVisualLimits)] / self.visualScale)
xlabel('km')
ylabel('km')
title(sprintf('%s', self.name))

if nargout == 2
    varargout{1} = Xf;
    varargout{2} = Yf;
end
