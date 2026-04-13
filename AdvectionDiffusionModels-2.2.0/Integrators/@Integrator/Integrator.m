classdef Integrator < handle
    % Integrate the ODE $$\frac{dy}{dt} = f(t,y)$$ with fixed-step RK4.
    %
    % `Integrator` advances the state `y` with the classical fourth-order
    % Runge-Kutta update
    %
    % $$
    % \begin{aligned}
    % k_1 &= f(t_n,y_n), \\
    % k_2 &= f\!\left(t_n + \tfrac{1}{2}dt, y_n + \tfrac{1}{2}dt\,k_1\right), \\
    % k_3 &= f\!\left(t_n + \tfrac{1}{2}dt, y_n + \tfrac{1}{2}dt\,k_2\right), \\
    % k_4 &= f(t_n + dt, y_n + dt\,k_3), \\
    % y_{n+1} &= y_n + \tfrac{dt}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right).
    % \end{aligned}
    % $$
    %
    % The state array `y` is stored as `nParticles x nDims`, and every
    % evaluation of `f(t,y)` must return an array with the same shape.
    %
    % ```matlab
    % f = @(t, y) [y(:,2) -y(:,1)];
    % y0 = [1 0; 0 1];
    % integrator = Integrator(f, y0, dt=0.1);
    % y = integrator.advanceToTime(1);
    % ```
    %
    % - Topic: Create the integrator
    % - Topic: Inspect integrator state
    % - Topic: Advance the integrator
    % - Declaration: classdef Integrator < handle

    properties (SetAccess = protected)
        % Fixed timestep `dt` used by every accepted step.
        %
        % This property stores the scalar timestep that appears in the RK4
        % update formulas documented for `Integrator`.
        %
        % - Topic: Inspect integrator state
        % - Declaration: self.stepSize
        % - Returns stepSize: positive scalar timestep $$dt$$
        stepSize (1,1) double {mustBePositive} = 1

        % Current integration time `t`.
        %
        % `currentTime` is updated after every accepted step and has the
        % same physical units as the `t` inputs passed to the public
        % stepping methods.
        %
        % - Topic: Inspect integrator state
        % - Declaration: self.currentTime
        % - Returns currentTime: scalar time $$t$$ associated with `currentY`
        currentTime (1,1) double = 0

        % Number of accepted timesteps taken so far.
        %
        % `totalIterations` counts accepted fixed-size updates. It does not
        % count interpolated output points requested by `integrateToTime`.
        %
        % - Topic: Inspect integrator state
        % - Declaration: self.totalIterations
        % - Returns totalIterations: nonnegative integer count of accepted steps
        totalIterations (1,1) double {mustBeInteger, mustBeNonnegative} = 0

        % Current state `y`.
        %
        % `currentY` stores the most recent accepted state array with shape
        % `nParticles x nDims`.
        %
        % - Topic: Inspect integrator state
        % - Declaration: self.currentY
        % - Returns currentY: state array $$y$$ with shape `nParticles x nDims`
        currentY double = []
    end

    properties (Access = protected)
        driftFunction = []
        nReps (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        nDims (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        rk4Stages double = zeros(0,0,4)
    end

    methods
        function self = Integrator(f, y0, options)
            % Create an RK4 integrator for $$\frac{dy}{dt} = f(t,y)$$.
            %
            % The constructor validates the initial right-hand side
            % evaluation `f(0,y0)` and stores the fixed timestep `dt`. The
            % state array `y0` is used without reshaping, so its
            % `nParticles x nDims` layout becomes the integration contract
            % for all later calls to `f`.
            %
            % - Topic: Create the integrator
            % - Declaration: self = Integrator(f,y0,dt=...)
            % - Parameter f: function handle for the deterministic right-hand side $$f(t,y)$$
            % - Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
            % - Parameter dt: positive scalar timestep $$dt$$
            arguments
                f (1,1) function_handle
                y0 {mustBeNumeric, mustBeReal}
                options.dt (1,1) double {mustBePositive}
            end

            try
                f0 = f(0, y0);
            catch exception
                throwAsCaller(MException( ...
                    "Integrator:DriftEvaluationFailed", ...
                    "Failed to evaluate f(0,y0): %s", ...
                    exception.message));
            end

            if ~isequal(size(y0), size(f0))
                throwAsCaller(MException( ...
                    "Integrator:InvalidDriftOutputSize", ...
                    "f(0,y0) must return an array with size %s, but it returned size %s.", ...
                    mat2str(size(y0)), ...
                    mat2str(size(f0))));
            end

            [self.nReps, self.nDims] = size(y0);
            self.stepSize = options.dt;
            self.driftFunction = f;
            self.currentY = y0;
            self.rk4Stages = zeros(self.nReps, self.nDims, 4);
        end

        function [y, t] = integrateToTime(self, t)
            % Integrate the ODE to each requested output time.
            %
            % The returned array `y` has shape `nParticles x nDims x nTimes`
            % and stores the accepted state associated with each requested
            % time entry. The returned vector `t` records the integrator
            % time after each request.
            %
            % - Topic: Advance the integrator
            % - Declaration: [y, t] = integrateToTime(self,t)
            % - Parameter t: vector of requested output times with the same units as `dt`
            % - Returns y: state history $$y(t)$$ stored as `nParticles x nDims x nTimes`
            % - Returns t: vector of accepted output times corresponding to the returned states
            arguments
                self
                t {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
            end

            requestedTimes = t;
            y = zeros(self.nReps, self.nDims, numel(requestedTimes));
            t = zeros(size(requestedTimes));
            for iTime = 1:numel(requestedTimes)
                y(:, :, iTime) = self.advanceToTime(requestedTimes(iTime));
                t(iTime) = self.currentTime;
            end
        end

        function y = advanceToTime(self, t)
            % Advance the integrator until the accepted time reaches `t`.
            %
            % `advanceToTime` repeatedly applies the fixed-step update until
            % `currentTime >= t`, then returns the accepted state. For the
            % deterministic RK4 integrator there is no substep
            % interpolation; the returned state is the most recent
            % accepted state `y`.
            %
            % - Topic: Advance the integrator
            % - Declaration: y = advanceToTime(self,t)
            % - Parameter t: scalar target time in the same units as `dt`
            % - Returns y: accepted state array $$y$$ after advancing to the requested time
            arguments
                self
                t (1,1) double {mustBeReal, mustBeFinite}
            end

            while self.currentTime < t
                self.currentY = self.stepForward(self.currentY, self.currentTime, self.stepSize);
                self.currentTime = self.currentTime + self.stepSize;
                self.totalIterations = self.totalIterations + 1;
            end

            y = self.currentY;
        end

        function y = advanceOneStep(self)
            % Advance the ODE by one fixed timestep `dt`.
            %
            % This method applies one accepted RK4 update and then returns
            % the new state `y`.
            %
            % - Topic: Advance the integrator
            % - Declaration: y = advanceOneStep(self)
            % - Returns y: accepted state array $$y_{n+1}$$ after one timestep
            self.currentY = self.stepForward(self.currentY, self.currentTime, self.stepSize);
            self.currentTime = self.currentTime + self.stepSize;
            self.totalIterations = self.totalIterations + 1;
            y = self.currentY;
        end
    end

    methods (Access = protected)
        function y = stepForward(self, y, t, dt)
            self.rk4Stages(:, :, 1) = self.driftFunction(t, y);
            self.rk4Stages(:, :, 2) = self.driftFunction(t + 0.5*dt, y + 0.5*dt*self.rk4Stages(:, :, 1));
            self.rk4Stages(:, :, 3) = self.driftFunction(t + 0.5*dt, y + 0.5*dt*self.rk4Stages(:, :, 2));
            self.rk4Stages(:, :, 4) = self.driftFunction(t + dt, y + dt*self.rk4Stages(:, :, 3));

            y = y + (dt/6)*(self.rk4Stages(:, :, 1) + 2*self.rk4Stages(:, :, 2) + 2*self.rk4Stages(:, :, 3) + self.rk4Stages(:, :, 4));
        end
    end
end
