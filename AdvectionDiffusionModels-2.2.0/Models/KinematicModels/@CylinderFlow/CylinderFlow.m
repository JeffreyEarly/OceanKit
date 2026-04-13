classdef CylinderFlow < StreamfunctionModel
    % Potential flow around a circular cylinder.
    %
    % The model streamfunction is
    %
    % $$ \psi(t,x,y) = U y \left(1 - \frac{R^2}{x^2 + y^2}\right), $$
    %
    % which gives the classical incompressible potential flow past a
    % cylinder of radius `R`. The cylinder interior is represented as a
    % polygonal obstacle.
    %
    % ```matlab
    % model = CylinderFlow();
    % integrator = AdvectionDiffusionIntegrator(model, 0);
    % x0 = [-2; -2] * model.R;
    % y0 = [-0.75; 0.75] * model.R;
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 12 * 3600, 600);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate the streamfunction
    % - Topic: Evaluate the velocity field
    % - Declaration: classdef CylinderFlow < StreamfunctionModel

    properties
        % Cylinder radius in meters.
        %
        % - Topic: Inspect model parameters
        R = 60e3

        % Far-field speed in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        U = 100 * (1e3 / 86400)
    end

    methods
        function self = CylinderFlow()
            % Create the default cylinder-flow model.
            %
            % - Topic: Create the model
            % - Declaration: self = CylinderFlow()
            % - Returns self: `CylinderFlow` instance
            self.xlim = 3 * self.R * [-1 1];
            self.ylim = 3 * self.R * [-1 1];
            self.isXPeriodic = true;

            self.xVisualLimits = 3 * self.R * [-1 1];
            self.yVisualLimits = 3 * self.R * [-1 1];
            self.visualScale = 1e3;

            theta = linspace(0, 2 * pi - 2 * pi / 30, 30);
            self.obstacles = polyshape(self.R * cos(theta), self.R * sin(theta));
            self.name = 'Cylinder flow';
        end

        function psiValue = psi(self, t, x, y)
            % Evaluate the potential-flow streamfunction.
            %
            % The implemented streamfunction is
            % $$\psi(t,x,y) = U y \left(1 - \frac{R^2}{x^2 + y^2}\right)$$.
            %
            % - Topic: Evaluate the streamfunction
            % - Declaration: psiValue = psi(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: streamfunction values with the same shape as `x` and `y`
            r2 = x.^2 + y.^2;
            psiValue = self.U * y .* (1 - self.R * self.R ./ r2);
        end

        function uValue = u(self, t, x, y)
            % Evaluate the x-velocity component.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x`
            r2 = x.^2 + y.^2;
            uValue = self.U - self.U * self.R * self.R * (x.^2 - y.^2) ./ (r2 .* r2);
        end

        function vValue = v(self, t, x, y)
            % Evaluate the y-velocity component.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `x`
            r2 = x.^2 + y.^2;
            vValue = -self.U * self.R * self.R * (2 * x .* y) ./ (r2 .* r2);
        end
    end
end
