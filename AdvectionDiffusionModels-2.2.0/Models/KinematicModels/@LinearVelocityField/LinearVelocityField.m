classdef LinearVelocityField < StreamfunctionModel
    % Two-dimensional linear velocity field with analytical solutions.
    %
    % This model defines the linear velocity field
    %
    % $$ \dot{\mathbf{x}} = \mathbf{u}_0 + A\mathbf{x}, $$
    %
    % with
    %
    % $$ \mathbf{u}_0 = \begin{bmatrix} u_0 \\ v_0 \end{bmatrix}, \qquad
    % A = \frac{1}{2}
    % \begin{bmatrix}
    % \sigma_n & \sigma_s - \zeta \\
    % \sigma_s + \zeta & -\sigma_n
    % \end{bmatrix}. $$
    %
    % The strain components satisfy
    %
    % $$ \sigma_n = \sigma \cos(2\theta), \qquad \sigma_s = \sigma \sin(2\theta). $$
    %
    % ### Analytical solutions
    %
    % The method `momentTensorEvolution` integrates the corresponding
    % second-moment system
    %
    % $$ \dot{M} = AM + MA^{\top} + 2\kappa I. $$
    %
    % The method `particlePath` returns the analytical solution for the
    % particle paths, as computed in [Oscroft, et al.
    % 2020](https://doi.org/10.3390/fluids6010014).
    %
    % ```matlab
    % model = LinearVelocityField(sigma=1e-5, theta=pi / 8, zeta=0);
    % integrator = AdvectionDiffusionIntegrator(model, 0);
    % x0 = [-20e3; 20e3];
    % y0 = [-10e3; 10e3];
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 12 * 3600, 600);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate the velocity field
    % - Topic: Evaluate the streamfunction
    % - Topic: Analyze particle and moment evolution
    % - Topic: Convert flow parameters
    % - Declaration: classdef LinearVelocityField < StreamfunctionModel

    properties
        % Strain magnitude in $$s^{-1}$$.
        %
        % Together with `theta`, this parameter determines `sigma_n` and
        % `sigma_s`.
        %
        % - Topic: Inspect model parameters
        sigma (1,1) double = 0

        % Strain orientation in radians.
        %
        % The angle enters through the doubled-angle relations for
        % `sigma_n` and `sigma_s`.
        %
        % - Topic: Inspect model parameters
        theta (1,1) double = 0

        % Relative vorticity in $$s^{-1}$$.
        %
        % `zeta` controls the antisymmetric part of the velocity
        % gradient.
        %
        % - Topic: Inspect model parameters
        zeta (1,1) double = 0

        % Uniform background x-velocity in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        u0 (1,1) double = 0

        % Uniform background y-velocity in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        v0 (1,1) double = 0
    end

    properties (Dependent)
        % Normal strain component $$\sigma_n = \sigma \cos(2\theta)$$.
        %
        % - Topic: Inspect model parameters
        sigma_n

        % Shear strain component $$\sigma_s = \sigma \sin(2\theta)$$.
        %
        % - Topic: Inspect model parameters
        sigma_s
    end

    methods
        function self = LinearVelocityField(options)
            % Create a linear velocity field.
            %
            % Pass any subset of the model parameters by name. Omitted
            % parameters default to zero.
            %
            % - Topic: Create the model
            % - Declaration: self = LinearVelocityField(sigma=...,theta=...,zeta=...,u0=...,v0=...)
            % - Parameter sigma: optional strain magnitude in $$s^{-1}$$; the default is `0`
            % - Parameter theta: optional strain orientation in radians; the default is `0`
            % - Parameter zeta: optional relative vorticity in $$s^{-1}$$; the default is `0`
            % - Parameter u0: optional uniform background x-velocity in $$m s^{-1}$$; the default is `0`
            % - Parameter v0: optional uniform background y-velocity in $$m s^{-1}$$; the default is `0`
            % - Returns self: `LinearVelocityField` instance
            arguments
                options.sigma (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
                options.theta (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
                options.zeta (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
                options.u0 (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
                options.v0 (1,1) {mustBeNumeric, mustBeReal, mustBeFinite} = 0
            end

            self.xVisualLimits = 1e3 * [-1 1];
            self.yVisualLimits = 1e3 * [-1 1];
            self.u0 = options.u0;
            self.v0 = options.v0;
            self.zeta = options.zeta;
            self.sigma = options.sigma;
            self.theta = options.theta;
            self.updateName();
        end

        function set.zeta(self, value)
            self.zeta = value;
            self.updateName();
        end

        function set.sigma(self, value)
            self.sigma = value;
            self.updateName();
        end

        function set.theta(self, value)
            self.theta = value;
            self.updateName();
        end

        function value = get.sigma_n(self)
            value = self.sigma * cos(2 * self.theta);
        end

        function value = get.sigma_s(self)
            value = self.sigma * sin(2 * self.theta);
        end

        function psiValue = psi(self, t, x, y)
            % Evaluate the quadratic streamfunction.
            %
            % The implemented streamfunction is
            %
            % $$ \psi = -u_0 y + v_0 x + \frac{1}{4}(\sigma_s + \zeta)x^2 - \frac{1}{2}\sigma_n xy - \frac{1}{4}(\sigma_s - \zeta)y^2. $$
            %
            % - Topic: Evaluate the streamfunction
            % - Declaration: psiValue = psi(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: streamfunction values with the same shape as `x` and `y`
            psiValue = -self.u0 * y + self.v0 * x + 0.25 * (self.sigma_s + self.zeta) * x .* x - 0.5 * self.sigma_n * x .* y - 0.25 * (self.sigma_s - self.zeta) * y .* y;
        end

        function uValue = u(self, t, x, y)
            % Evaluate the x-velocity.
            %
            % The implemented velocity component is
            %
            % $$ u = u_0 + \frac{1}{2}\left(\sigma_n x + (\sigma_s - \zeta)y\right). $$
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`
            uValue = self.u0 + 0.5 * (self.sigma_n * x + (self.sigma_s - self.zeta) * y);
        end

        function vValue = v(self, t, x, y)
            % Evaluate the y-velocity.
            %
            % The implemented velocity component is
            %
            % $$ v = v_0 + \frac{1}{2}\left((\sigma_s + \zeta)x - \sigma_n y\right). $$
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`
            vValue = self.v0 + 0.5 * ((self.sigma_s + self.zeta) * x - self.sigma_n * y);
        end
    end

    methods (Static)
        function [sigma_n, sigma_s] = normalAndShearFromSigmaTheta(sigma, theta)
            % Convert strain magnitude and orientation to strain components.
            %
            % This helper evaluates $$\sigma_n = \sigma \cos(2\theta)$$ and
            % $$\sigma_s = \sigma \sin(2\theta)$$.
            %
            % - Topic: Convert flow parameters
            % - Declaration: [sigma_n, sigma_s] = normalAndShearFromSigmaTheta(sigma,theta)
            % - Parameter sigma: strain magnitude in $$s^{-1}$$
            % - Parameter theta: strain orientation in radians
            % - Returns sigma_n: normal strain component in $$s^{-1}$$
            % - Returns sigma_s: shear strain component in $$s^{-1}$$
            sigma_n = sigma * cos(2 * theta);
            sigma_s = sigma * sin(2 * theta);
        end

        function [sigma, theta] = sigmaThetaFromNormalAndShear(sigma_n, sigma_s)
            % Convert normal and shear strain to magnitude and orientation.
            %
            % This helper inverts the relations
            % $$\sigma_n = \sigma \cos(2\theta)$$ and
            % $$\sigma_s = \sigma \sin(2\theta)$$.
            %
            % - Topic: Convert flow parameters
            % - Declaration: [sigma, theta] = sigmaThetaFromNormalAndShear(sigma_n,sigma_s)
            % - Parameter sigma_n: normal strain component in $$s^{-1}$$
            % - Parameter sigma_s: shear strain component in $$s^{-1}$$
            % - Returns sigma: strain magnitude in $$s^{-1}$$
            % - Returns theta: strain orientation in radians
            sigma = sqrt(sigma_n .* sigma_n + sigma_s .* sigma_s);
            theta = atan2(sigma_s, sigma_n) / 2;
        end
    end

    methods (Access = private)
        function updateName(self)
            if self.sigma * self.sigma > self.zeta * self.zeta
                self.name = 'Strain-dominated linear velocity field';
            elseif self.sigma * self.sigma < self.zeta * self.zeta
                self.name = 'Vorticity-dominated linear velocity field';
            else
                self.name = 'Strain-vorticity matched linear velocity field';
            end
        end
    end
end
