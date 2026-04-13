classdef IntegratorEulerMaruyama < Integrator
    % Integrate the Itô SDE $$dy = f(t,y)\,dt + g(t,y)\,dW_t$$.
    %
    % `IntegratorEulerMaruyama` advances the stochastic differential
    % equation
    %
    % $$
    % dy = f(t,y)\,dt + g(t,y)\,dW_t
    % $$
    %
    % with the Euler-Maruyama update
    %
    % $$
    % y_{n+1} = y_n + f(t_n,y_n)\,dt + g(t_n,y_n)\sqrt{dt}\,\eta_n,
    % $$
    %
    % where $$\eta_n$$ is standard normal with the same shape as `y_n`.
    % `advanceToTime` preserves the existing first-order linear
    % interpolation between accepted stochastic steps.
    %
    % ```matlab
    % f = @(t, y) -0.1 * y;
    % g = @(t, y) 0.05 + zeros(size(y));
    % y0 = [1 0; 0 1];
    % integrator = IntegratorEulerMaruyama(f, g, y0, dt=0.1);
    % y = integrator.advanceOneStep();
    % ```
    %
    % - Topic: Create the integrator
    % - Topic: Advance the integrator
    % - Declaration: classdef IntegratorEulerMaruyama < Integrator

    properties (Access = protected)
        diffusionFunction = []
        previousY double = []
        previousTime double = []
    end

    methods
        function self = IntegratorEulerMaruyama(f, g, y0, options)
            % Create an Euler-Maruyama integrator for $$dy = f\,dt + g\,dW_t$$.
            %
            % The constructor validates both `f(0,y0)` and `g(0,y0)` and
            % requires them to return arrays with the same shape as `y0`.
            %
            % - Topic: Create the integrator
            % - Declaration: self = IntegratorEulerMaruyama(f,g,y0,dt=...)
            % - Parameter f: deterministic drift function $$f(t,y)$$
            % - Parameter g: stochastic amplitude function $$g(t,y)$$
            % - Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
            % - Parameter dt: positive scalar timestep $$dt$$
            arguments
                f (1,1) function_handle
                g (1,1) function_handle
                y0 {mustBeNumeric, mustBeReal}
                options.dt (1,1) double {mustBePositive}
            end

            self@Integrator(f, y0, dt=options.dt);

            try
                g0 = g(0, y0);
            catch exception
                throwAsCaller(MException( ...
                    "IntegratorEulerMaruyama:DiffusionEvaluationFailed", ...
                    "Failed to evaluate g(0,y0): %s", ...
                    exception.message));
            end

            if ~isequal(size(y0), size(g0))
                throwAsCaller(MException( ...
                    "IntegratorEulerMaruyama:InvalidDiffusionOutputSize", ...
                    "g(0,y0) must return an array with size %s, but it returned size %s.", ...
                    mat2str(size(y0)), ...
                    mat2str(size(g0))));
            end

            self.diffusionFunction = g;
        end

        function y = advanceToTime(self, t)
            % Advance the SDE until the requested output time `t`.
            %
            % Accepted Euler-Maruyama steps are taken until
            % `currentTime >= t`. When `t` lies between the two most recent
            % accepted stochastic states, the returned value is the current
            % first-order linear interpolation between them.
            %
            % - Topic: Advance the integrator
            % - Declaration: y = advanceToTime(self,t)
            % - Parameter t: scalar target time in the same units as `dt`
            % - Returns y: state estimate $$y(t)$$ returned with first-order interpolation between accepted steps
            arguments
                self
                t (1,1) double {mustBeReal, mustBeFinite}
            end

            while self.currentTime < t
                self.previousY = self.currentY;
                self.previousTime = self.currentTime;
                self.currentY = self.stepForward(self.currentY, self.currentTime, self.stepSize);
                self.currentTime = self.currentTime + self.stepSize;
                self.totalIterations = self.totalIterations + 1;
            end

            if isempty(self.previousY)
                y = self.currentY;
                return
            end

            alpha = (self.currentTime - t)/(self.currentTime - self.previousTime);
            y = alpha*self.previousY + (1 - alpha)*self.currentY;
        end

        function y = advanceOneStep(self)
            % Advance the SDE by one fixed timestep `dt`.
            %
            % This method stores the previous accepted state so that a later
            % call to `advanceToTime` can apply the same linear
            % interpolation used by the legacy implementation.
            %
            % - Topic: Advance the integrator
            % - Declaration: y = advanceOneStep(self)
            % - Returns y: accepted state array $$y_{n+1}$$ after one stochastic timestep
            self.previousY = self.currentY;
            self.previousTime = self.currentTime;
            self.currentY = self.stepForward(self.currentY, self.currentTime, self.stepSize);
            self.currentTime = self.currentTime + self.stepSize;
            self.totalIterations = self.totalIterations + 1;
            y = self.currentY;
        end
    end

    methods (Access = protected)
        function y = stepForward(self, y, t, dt)
            y = y + dt*self.driftFunction(t, y) + sqrt(dt)*randn(size(y)).*self.diffusionFunction(t, y);
        end
    end
end
