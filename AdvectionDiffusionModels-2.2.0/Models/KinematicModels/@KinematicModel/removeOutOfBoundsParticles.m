function [x0, y0] = removeOutOfBoundsParticles(self, x0, y0)
% Remove initial particles outside the valid domain.
%
% Particles are removed when they lie outside `xlim` and `ylim` or inside
% any polygonal obstacle. Periodic directions are wrapped before obstacle
% intersection tests.
%
% - Topic: Work with particle domains
% - Declaration: [x0, y0] = removeOutOfBoundsParticles(self,x0,y0)
% - Parameter x0: column vector of initial x positions in meters
% - Parameter y0: column vector of initial y positions in meters
% - Returns x0: filtered column vector of x positions in meters
% - Returns y0: filtered column vector of y positions in meters
arguments
    self (1,1) KinematicModel
    x0 (:,1) double
    y0 (:,1) double
end

outOfBounds = x0 < min(self.xlim) | x0 > max(self.xlim) | y0 < min(self.ylim) | y0 > max(self.ylim);
if ~isempty(self.obstacles)
    if self.isXPeriodic
        xWrapped = mod(x0 - min(self.xlim), max(self.xlim) - min(self.xlim)) + min(self.xlim);
    else
        xWrapped = x0;
    end

    if self.isYPeriodic
        yWrapped = mod(y0 - min(self.ylim), max(self.ylim) - min(self.ylim)) + min(self.ylim);
    else
        yWrapped = y0;
    end

    for iObstacle = 1:length(self.obstacles)
        outOfBounds = outOfBounds | isinterior(self.obstacles(iObstacle), xWrapped, yWrapped);
    end
end

if sum(outOfBounds) > 0
    fprintf('Removed %d particles because they were out of bounds.\n', sum(outOfBounds));
    x0(outOfBounds) = [];
    y0(outOfBounds) = [];
end
