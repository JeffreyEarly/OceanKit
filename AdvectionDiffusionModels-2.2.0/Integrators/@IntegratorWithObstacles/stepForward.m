function y = stepForward(self, y, t, dt)
% Apply one obstacle-aware stochastic step.
%
% This developer method first advances the deterministic drift with the RK4
% update from `Integrator`, then adds the stored diffusion increment, wraps
% periodic coordinates, and repeatedly reflects any step segment that would
% enter a polygonal obstacle.
%
% - Topic: Handle obstacle reflections
% - Developer: true
% - Declaration: y = stepForward(self,y,t,dt)
% - Parameter y: accepted state array at time `t`
% - Parameter t: scalar time associated with the accepted state
% - Parameter dt: positive scalar timestep
% - Returns y: accepted state array after drift, diffusion, wrapping, and obstacle reflection
yAccepted = y;
y = stepForward@Integrator(self, y, t, dt);

dy = y - yAccepted;
for iDim = 1:self.nDims
    dy(:, iDim) = dy(:, iDim) + dt*self.diffusivityFlux{iDim}(y(:, iDim), dt);
end

particlesToCheck = transpose(1:size(y, 1));
nLoops = 0;
while ~isempty(particlesToCheck)
    nLoops = nLoops + 1;

    for iDim = 1:self.nDims
        y(particlesToCheck, iDim) = yAccepted(particlesToCheck, iDim) + dy(particlesToCheck, iDim);
        if self.isPeriodic(iDim)
            y(particlesToCheck, iDim) = mod(y(particlesToCheck, iDim) - self.ymin(iDim), self.ymax(iDim) - self.ymin(iDim)) + self.ymin(iDim);
        end
    end

    didReflectOnSomeObstacle = false(size(particlesToCheck));
    for iObstacle = 1:length(self.obstacleData)
        yCurrent = yAccepted(particlesToCheck, :) + dy(particlesToCheck, :);
        for iDim = 1:self.nDims
            if self.isPeriodic(iDim)
                yCurrent(:, iDim) = mod(yCurrent(:, iDim) - self.ymin(iDim), self.ymax(iDim) - self.ymin(iDim)) + self.ymin(iDim);
            end
        end

        [yNew, dyNew, didReflect] = self.updateWithReflection(self.obstacleData(iObstacle), yAccepted(particlesToCheck, :), yCurrent, dy(particlesToCheck, :));
        if any(didReflect)
            reflectedIndices = particlesToCheck(didReflect);
            yAccepted(reflectedIndices, :) = yNew(didReflect, :);
            dy(reflectedIndices, :) = dyNew(didReflect, :);
            didReflectOnSomeObstacle = didReflectOnSomeObstacle | didReflect;
        end
    end

    particlesToCheck(~didReflectOnSomeObstacle) = [];
    if nLoops > 20
        warning("IntegratorWithObstacles:TooManyReflectionLoops", "Too many reflection loops. Something bad must of happened.");
        particlesToCheck = [];
    end
end

y = yAccepted + dy;
end
