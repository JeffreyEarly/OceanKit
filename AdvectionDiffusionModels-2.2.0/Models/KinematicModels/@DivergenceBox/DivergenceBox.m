classdef DivergenceBox < KinematicModel
    % Box with alternating Gaussian convergence and divergence cells.
    %
    % The model defines a steady velocity field by summing Gaussian cells
    % centered at `(x_c, y_c)`:
    %
    % $$ u(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{x - x_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$
    %
    % $$ v(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{y - y_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$
    %
    % where `s_i = (-1)^i` alternates the sign of neighboring cells.
    %
    % ```matlab
    % model = DivergenceBox();
    % integrator = AdvectionDiffusionIntegrator(model, 0);
    % x0 = [0.2; 0.8] * model.Lx;
    % y0 = [0.25; 0.75] * model.Ly;
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 3600, 300);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate the velocity field
    % - Declaration: classdef DivergenceBox < KinematicModel

    properties
        % Domain width in meters.
        %
        % - Topic: Inspect model parameters
        Lx = 10e3

        % Domain height in meters.
        %
        % - Topic: Inspect model parameters
        Ly = 5e3

        % Gaussian length scale in meters.
        %
        % - Topic: Inspect model parameters
        Lg = 0.5e3

        % Number of cells along x.
        %
        % - Topic: Inspect model parameters
        n = 2

        % Number of cells along y.
        %
        % - Topic: Inspect model parameters
        m = 1

        % Velocity scale in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        U = 1.0
    end

    methods
        function self = DivergenceBox()
            % Create the default divergence-box model.
            %
            % - Topic: Create the model
            % - Declaration: self = DivergenceBox()
            % - Returns self: `DivergenceBox` instance
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];

            epsilonX = 0.05 * self.Lx;
            epsilonY = 0.05 * self.Ly;
            self.xVisualLimits = [0 - epsilonX self.Lx + epsilonX];
            self.yVisualLimits = [0 - epsilonY self.Ly + epsilonY];
            self.visualScale = 1e3;
            self.name = 'Divergence box';
        end

        function uValue = u(self, t, x, y)
            % Evaluate the x-velocity component.
            %
            % Neighboring Gaussian cells alternate sign along the
            % x-direction.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x`
            uValue = zeros(size(x));
            for i = 1:self.n
                for j = 1:self.m
                    xc = (2 * i - 1) * self.Lx / (self.n + 2);
                    yc = j * self.Ly / (self.m + 1);
                    r2 = (x - xc).^2 + (y - yc).^2;
                    gaussian = exp(-r2 / (2 * self.Lg * self.Lg));
                    uValue = uValue + ((-1).^i) * exp(0.5) * self.U * ((x - xc) / self.Lg) .* gaussian;
                end
            end
        end

        function vValue = v(self, t, x, y)
            % Evaluate the y-velocity component.
            %
            % Neighboring Gaussian cells alternate sign along the
            % x-direction.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `x`
            vValue = zeros(size(x));
            for i = 1:self.n
                for j = 1:self.m
                    xc = (2 * i - 1) * self.Lx / (self.n + 2);
                    yc = j * self.Ly / (self.m + 1);
                    r2 = (x - xc).^2 + (y - yc).^2;
                    gaussian = exp(-r2 / (2 * self.Lg * self.Lg));
                    vValue = vValue + ((-1).^i) * exp(0.5) * self.U * ((y - yc) / self.Lg) .* gaussian;
                end
            end
        end
    end
end
