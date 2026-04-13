function plotTrajectories(self, x, y, varargin)
% Plot particle trajectories on the current axes.
%
% When `isXPeriodic` is true, trajectories are redrawn in wrapped copies
% so the visible path remains continuous across the plotting window.
%
% - Topic: Plot model diagnostics
% - Declaration: plotTrajectories(self,x,y,varargin)
% - Parameter x: trajectory x positions in meters with shape `[nTimes nParticles]`
% - Parameter y: trajectory y positions in meters with shape `[nTimes nParticles]`
% - Parameter varargin: forwarded line styling arguments
arguments
    self (1,1) KinematicModel
    x double
    y double
end
arguments (Repeating)
    varargin
end

if self.isXPeriodic
    xMin = min(self.xlim);
    xMax = max(self.xlim);
    xWindowLength = xMax - xMin;
    nmin = floor(min(x(:) - xMin) / xWindowLength);
    nmax = floor(max(x(:) - xMin) / xWindowLength);
    axesHandle = gca;

    for n = nmin:nmax
        xShifted = x - xMin - n * xWindowLength;
        yShifted = y;
        mask = xShifted < 0 | xShifted > xWindowLength;
        xShifted(mask) = nan;
        yShifted(mask) = nan;
        axesHandle.ColorOrderIndex = 1;
        plot((xShifted + xMin) / self.visualScale, yShifted / self.visualScale, varargin{:})
    end

    nmin = floor(min(x(end,:) - xMin) / xWindowLength);
    nmax = floor(max(x(end,:) - xMin) / xWindowLength);
    for n = nmin:nmax
        xShifted = x(end,:) - xMin - n * xWindowLength;
        yShifted = y(end,:);
        mask = xShifted < 0 | xShifted > xWindowLength;
        xShifted(mask) = nan;
        yShifted(mask) = nan;
        scatter((xShifted + xMin) / self.visualScale, yShifted / self.visualScale, 8^2, 'k', 'filled')
    end
else
    if size(x,1) > 1
        plot(x / self.visualScale, y / self.visualScale, varargin{:})
    end
    hold on
    scatter(x(end,:) / self.visualScale, y(end,:) / self.visualScale, 8^2, 'k', 'filled')
    axis equal
end
