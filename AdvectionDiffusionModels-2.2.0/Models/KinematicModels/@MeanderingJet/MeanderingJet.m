classdef MeanderingJet < StreamfunctionModel
    % Kinematic streamfunction model for a meandering jet.
    %
    % This model uses the Bower meandering-jet streamfunction
    %
    % $$ \psi(t,x,y) = UL\left[1 - \tanh(\gamma)\right], $$
    %
    % with
    %
    % $$ \gamma = \frac{y - A\cos(\theta)}{L\sqrt{1 + (kA\sin(\theta))^2}}, \qquad \theta = k(x - c_x t), \qquad k = \frac{2\pi}{L_x}. $$
    %
    % The resulting velocity field is incompressible and periodic in the
    % x-direction over the interval `[0 2*Lx]`.
    %
    % References:
    %   Amy Bower, "A simple kinematic mechanism for mixing fluid parcels
    %   across a meandering jet."
    %
    % ```matlab
    % model = MeanderingJet();
    % integrator = AdvectionDiffusionIntegrator(model, 0);
    % x0 = [0.25; 1.25] * model.Lx;
    % y0 = [-0.5; 0.5] * model.L;
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 3 * 86400, 1800);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate meander coordinates
    % - Topic: Evaluate the streamfunction
    % - Topic: Evaluate the velocity field
    % - Declaration: classdef MeanderingJet < StreamfunctionModel

    properties
        % Jet half-width in meters.
        %
        % - Topic: Inspect model parameters
        L = 40e3

        % Velocity scale in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        U = 100 * (1e3 / 86400)

        % Meander amplitude in meters.
        %
        % - Topic: Inspect model parameters
        A = 50e3

        % Meander wavelength in meters.
        %
        % - Topic: Inspect model parameters
        Lx = 350e3

        % Meander phase speed in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        cx = 0 * 10 * (1e3 / 86400)
    end

    properties (Dependent)
        % Meander wavenumber in $$m^{-1}$$.
        %
        % - Topic: Inspect model parameters
        k
    end

    methods
        function self = MeanderingJet()
            % Create the default meandering jet model.
            %
            % - Topic: Create the model
            % - Declaration: self = MeanderingJet()
            % - Returns self: `MeanderingJet` instance
            self.xlim = [0 2 * self.Lx];
            self.ylim = [-5 * self.L 5 * self.L];
            self.isXPeriodic = true;

            self.xVisualLimits = [0 2 * self.Lx];
            self.yVisualLimits = [-5 * self.L 5 * self.L];
            self.visualScale = 1e3;
            self.name = 'Meandering jet';
        end

        function value = get.k(self)
            value = 2 * pi / self.Lx;
        end

        function thetaValue = theta(self, t, x)
            % Evaluate the meander phase.
            %
            % The phase is $$\theta = k(x - c_x t)$$.
            %
            % - Topic: Evaluate meander coordinates
            % - Declaration: thetaValue = theta(self,t,x)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Returns thetaValue: meander phase in radians with the same shape as `x`
            thetaValue = self.k * (x - self.cx * t);
        end

        function gammaValue = gamma(self, t, x, y)
            % Evaluate the nondimensional cross-jet coordinate.
            %
            % `gamma` is the transformed cross-jet coordinate that appears
            % inside the Bower streamfunction.
            %
            % - Topic: Evaluate meander coordinates
            % - Declaration: gammaValue = gamma(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns gammaValue: nondimensional coordinate with the same shape as `x` and `y`
            gammaValue = (y - self.A * cos(self.theta(t, x))) / self.L ./ sqrt(1 + (self.k * self.A * sin(self.theta(t, x))).^2);
        end

        function psiValue = psi(self, t, x, y)
            % Evaluate the meandering-jet streamfunction.
            %
            % - Topic: Evaluate the streamfunction
            % - Declaration: psiValue = psi(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: streamfunction values with the same shape as `x` and `y`
            psiValue = self.U * self.L * (1 - tanh(self.gamma(t, x, y)));
        end

        function uValue = u(self, t, x, y)
            % Evaluate the along-jet velocity component.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x`
            uValue = self.U * (1 + (self.k * self.A * sin(self.theta(t, x))).^2).^(-1/2) .* sech(self.gamma(t, x, y)).^2;
        end

        function vValue = v(self, t, x, y)
            % Evaluate the cross-jet velocity component.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `x`
            vValue = -self.U * (self.A * self.k * (1 + self.A * self.A * self.k * self.k) * sin(self.theta(t, x))) .* (1 + (self.k * self.A * sin(self.theta(t, x))).^2).^(-3/2) .* sech(self.gamma(t, x, y)).^2;
        end
    end
end
