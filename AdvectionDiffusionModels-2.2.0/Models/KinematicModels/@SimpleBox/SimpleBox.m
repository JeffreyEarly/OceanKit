classdef SimpleBox < KinematicModel
    % Rectangular box with zero deterministic velocity.
    %
    % This model defines
    %
    % $$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$
    %
    % on a finite rectangular domain. Any motion comes from the external
    % integrator, typically through the diffusivity term.
    %
    % ```matlab
    % model = SimpleBox();
    % integrator = AdvectionDiffusionIntegrator(model, 20);
    % x0 = [0.25; 0.75] * model.Lx;
    % y0 = [0.25; 0.75] * model.Ly;
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 6 * 3600, 300);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate the velocity field
    % - Declaration: classdef SimpleBox < KinematicModel

    properties
        % Domain width in meters.
        %
        % - Topic: Inspect model parameters
        Lx = 10e3

        % Domain height in meters.
        %
        % - Topic: Inspect model parameters
        Ly = 5e3
    end

    methods
        function self = SimpleBox()
            % Create the default box model.
            %
            % - Topic: Create the model
            % - Declaration: self = SimpleBox()
            % - Returns self: `SimpleBox` instance
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];

            epsilonX = 0.05 * self.Lx;
            epsilonY = 0.05 * self.Ly;
            self.xVisualLimits = [0 - epsilonX self.Lx + epsilonX];
            self.yVisualLimits = [0 - epsilonY self.Ly + epsilonY];
            self.visualScale = 1e3;
            self.name = 'No-flow box';
        end

        function uValue = u(self, t, x, y)
            % Evaluate the zero x-velocity field.
            %
            % This model sets $$u(t,x,y) = 0$$ everywhere in the box.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x`
            uValue = zeros(size(x));
        end

        function vValue = v(self, t, x, y)
            % Evaluate the zero y-velocity field.
            %
            % This model sets $$v(t,x,y) = 0$$ everywhere in the box.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `y`
            vValue = zeros(size(y));
        end
    end
end
