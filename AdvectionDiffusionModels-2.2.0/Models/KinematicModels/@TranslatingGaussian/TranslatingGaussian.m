classdef TranslatingGaussian < StreamfunctionModel
    % Translating Gaussian eddy streamfunction.
    %
    % The streamfunction is
    %
    % $$ \psi(t,x,y) = e^{1/2} U L \exp\left(-\frac{(x - c_x t)^2 + (y - c_y t)^2}{2L^2}\right), $$
    %
    % which produces a translating Gaussian vortex with peak velocity scale
    % `U` and translation velocity `(cx, cy)`.
    %
    % ```matlab
    % model = TranslatingGaussian();
    % integrator = AdvectionDiffusionIntegrator(model, 0);
    % x0 = [-model.L; model.L];
    % y0 = [0; 0];
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 5 * 86400, 1800);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Create the model
    % - Topic: Inspect model parameters
    % - Topic: Evaluate the streamfunction
    % - Topic: Evaluate the velocity field
    % - Declaration: classdef TranslatingGaussian < StreamfunctionModel

    properties
        % Eddy length scale in meters.
        %
        % - Topic: Inspect model parameters
        L = 60e3

        % Peak velocity scale in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        U = 0.12

        % Translation speed in x in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        cx = -0.0267

        % Translation speed in y in $$m s^{-1}$$.
        %
        % - Topic: Inspect model parameters
        cy = 0
    end

    methods
        function self = TranslatingGaussian()
            % Create the default translating Gaussian model.
            %
            % - Topic: Create the model
            % - Declaration: self = TranslatingGaussian()
            % - Returns self: `TranslatingGaussian` instance
            self.xVisualLimits = 3 * self.L * [-1 1];
            self.yVisualLimits = 3 * self.L * [-1 1];
            self.visualScale = 1e3;
            self.name = 'Translating Gaussian';
        end

        function psiValue = psi(self, t, x, y)
            % Evaluate the Gaussian streamfunction.
            %
            % The eddy center translates at the constant speed
            % `[(cx) (cy)]`.
            %
            % - Topic: Evaluate the streamfunction
            % - Declaration: psiValue = psi(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: streamfunction values with the same shape as `x` and `y`
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            psiValue = exp(0.5) * self.U * self.L * gaussian;
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
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            uValue = exp(0.5) * self.U * ((y - self.cy * t) / self.L) .* gaussian;
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
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            vValue = -exp(0.5) * self.U * ((x - self.cx * t) / self.L) .* gaussian;
        end
    end
end
