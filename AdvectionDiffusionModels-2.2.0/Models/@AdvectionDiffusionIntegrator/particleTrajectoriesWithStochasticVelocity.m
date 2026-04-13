function [t, x, y] = particleTrajectoriesWithStochasticVelocity(self, x0, y0, T, dt, u, v)
% Integrate particles with extra stochastic velocity terms.
%
% This method integrates
%
% $$ d\mathbf{r}(t) = [\mathbf{u}_{model}(t,\mathbf{r}) + \mathbf{u}_{stochastic}(t)]\,dt + \sqrt{2\kappa}\,d\mathbf{W}_t, $$
%
% where `u(t)` and `v(t)` return column vectors with one value per active
% particle at the current output time.
%
% - Topic: Integrate particle trajectories
% - Declaration: [t, x, y] = particleTrajectoriesWithStochasticVelocity(self,x0,y0,T,dt,u,v)
% - Parameter x0: initial x positions in meters; any input shape is flattened to a column vector
% - Parameter y0: initial y positions in meters; any input shape is flattened to a column vector
% - Parameter T: total integration duration in seconds
% - Parameter dt: output time increment in seconds
% - Parameter u: function handle `u(t)` returning one stochastic x-velocity value per particle
% - Parameter v: function handle `v(t)` returning one stochastic y-velocity value per particle
% - Returns t: column vector of output times in seconds
% - Returns x: x positions in meters with shape `[length(t) nParticles]`
% - Returns y: y positions in meters with shape `[length(t) nParticles]`
arguments
    self (1,1) AdvectionDiffusionIntegrator
    x0 {mustBeNumeric}
    y0 {mustBeNumeric}
    T (1,1) double {mustBeNonnegative}
    dt (1,1) double {mustBePositive}
    u (1,1) function_handle
    v (1,1) function_handle
end

x0 = reshape(x0, [], 1);
y0 = reshape(y0, [], 1);

if ~self.kinematicModel.isSpherical
    flux = @(tValue, p) cat(2, u(tValue) + self.kinematicModel.u(tValue, p(:,1), p(:,2)), v(tValue) + self.kinematicModel.v(tValue, p(:,1), p(:,2)));
else
    flux = @(tValue, p) cat(2, (u(tValue) + self.kinematicModel.u(tValue, p(:,1), p(:,2))) ./ cos(p(:,2) / self.kinematicModel.R), v(tValue) + self.kinematicModel.v(tValue, p(:,1), p(:,2)));
end

[x0, y0] = self.kinematicModel.removeOutOfBoundsParticles(x0, y0);
if isempty(x0)
    fprintf('There were no particles left to advect. Aborting.\n');
end

nParticles = length(x0);
p0 = cat(2, x0, y0);
nTimes = round(T / dt);

xMin = min(self.kinematicModel.xlim);
xMax = max(self.kinematicModel.xlim);
yMin = min(self.kinematicModel.ylim);
yMax = max(self.kinematicModel.ylim);

if self.stepSize == 0
    singleIncrement = true;
    step = dt;
else
    singleIncrement = false;
    step = self.stepSize;
end

ymin = [xMin yMin];
ymax = [xMax yMax];
obstacles = self.kinematicModel.obstacles;
isPeriodic = [self.kinematicModel.isXPeriodic self.kinematicModel.isYPeriodic];
integrator = IntegratorWithObstacles(flux, p0, dt=step, kappa=self.kappa, ymin=ymin, ymax=ymax, obstacles=obstacles, isPeriodic=isPeriodic);

t = zeros(nTimes + 1, 1);
x = zeros(length(t), nParticles);
y = zeros(length(t), nParticles);

x(1,:) = x0;
y(1,:) = y0;
for iTime = 1:nTimes
    if singleIncrement
        integrator.advanceOneStep();
    else
        integrator.advanceToTime(iTime * dt);
    end

    p = integrator.currentY;
    x(iTime + 1,:) = p(:,1).';
    y(iTime + 1,:) = p(:,2).';
    t(iTime + 1) = integrator.currentTime;
end
