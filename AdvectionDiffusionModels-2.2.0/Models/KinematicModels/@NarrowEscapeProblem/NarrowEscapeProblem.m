classdef NarrowEscapeProblem < KinematicModel
    % Rectangular chamber with a narrow opening in the left wall.
    %
    % This model has zero deterministic velocity,
    %
    % $$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$
    %
    % and represents the geometry through polygonal obstacles that leave a
    % gap of width `W` in the left wall. The problem is useful for testing
    % stochastic escape through a small opening.
    %
    % ```matlab
    % model = NarrowEscapeProblem();
    % integrator = AdvectionDiffusionIntegrator(model, 20);
    % x0 = [0.25; 0.75] * model.L;
    % y0 = [0.4; 0.6] * model.L;
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 6 * 3600, 300);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Rebuild model geometry
    % - Declaration: classdef NarrowEscapeProblem < KinematicModel

    properties
        % Chamber length in meters.
        %
        % - Topic: Inspect model parameters
        L = 10e3

        % Wall thickness in meters.
        %
        % - Topic: Inspect model parameters
        delta = 1e3

        % Opening width in meters.
        %
        % - Topic: Inspect model parameters
        W = 3e3
    end

    methods
        function self = NarrowEscapeProblem()
            % Create the default narrow-escape geometry.
            %
            % - Topic: Create the model
            % - Declaration: self = NarrowEscapeProblem()
            % - Returns self: `NarrowEscapeProblem` instance
            self.xVisualLimits = 0.1 * (self.L + 2 * self.delta) * [-1 1] + [-self.delta self.L + self.delta];
            self.yVisualLimits = 0.1 * (self.L + 2 * self.delta) * [-1 1] + [-self.delta self.L + self.delta];
            self.visualScale = 1e3;
            self.regenerateObstacles();
            self.name = 'Narrow escape';
        end

        function set.W(self, W)
            self.W = W;
            self.regenerateObstacles();
        end

        function set.L(self, L)
            self.L = L;
            self.regenerateObstacles();
        end

        function set.delta(self, delta)
            self.delta = delta;
            self.regenerateObstacles();
        end
    end
end
