classdef KinematicModel < handle
    % Base class for two-dimensional deterministic velocity fields.
    %
    % A kinematic model defines particle motion through
    %
    % $$ \frac{dx}{dt} = u(t,x,y), \qquad \frac{dy}{dt} = v(t,x,y). $$
    %
    % Subclasses specify the velocity field itself and, optionally, finite
    % domain limits, periodic directions, polygonal obstacles, and plotting
    % limits. The base class provides common geometry filtering and
    % visualization helpers.
    %
    % ```matlab
    % model = KinematicModel();
    % integrator = AdvectionDiffusionIntegrator(model, 20);
    % x0 = [0; 5e3];
    % y0 = [0; 5e3];
    % [~, x, y] = integrator.particleTrajectories(x0, y0, 3600, 300);
    % figure
    % model.plotTrajectories(x, y)
    % ```
    %
    % - Topic: Configure model domains
    % - Topic: Inspect model settings
    % - Topic: Evaluate velocity fields
    % - Topic: Plot model diagnostics
    % - Topic: Work with particle domains
    % - Declaration: classdef KinematicModel < handle

    properties
        % Finite or infinite x-domain limits in meters.
        %
        % If both limits are finite and `isXPeriodic` is true, particles
        % wrap across this interval.
        %
        % - Topic: Configure model domains
        xlim (1,2) double = [-Inf Inf]

        % Finite or infinite y-domain limits in meters.
        %
        % If both limits are finite and `isYPeriodic` is true, particles
        % wrap across this interval.
        %
        % - Topic: Configure model domains
        ylim (1,2) double = [-Inf Inf]

        % Indicates whether the x-direction wraps periodically across `xlim`.
        %
        % - Topic: Configure model domains
        isXPeriodic (1,1) logical = false

        % Indicates whether the y-direction wraps periodically across `ylim`.
        %
        % - Topic: Configure model domains
        isYPeriodic (1,1) logical = false

        % Polygonal obstacle geometry as a `polyshape` array.
        %
        % Particles inside these polygons are removed before trajectory
        % integration.
        %
        % - Topic: Configure model domains
        obstacles = []

        % Indicates that the x-coordinate represents longitude on a sphere.
        %
        % When this flag is true, `AdvectionDiffusionIntegrator` applies
        % the spherical x-scaling expected by the model.
        %
        % - Topic: Configure model domains
        isSpherical (1,1) logical = false

        % Finite x-limits in meters used for plotting.
        %
        % These limits affect the plotting helpers but not the physical
        % integration domain.
        %
        % - Topic: Inspect model settings
        xVisualLimits (1,2) double = [-1 1]

        % Finite y-limits in meters used for plotting.
        %
        % These limits affect the plotting helpers but not the physical
        % integration domain.
        %
        % - Topic: Inspect model settings
        yVisualLimits (1,2) double = [-1 1]

        % Plotting scale factor.
        %
        % A value of `1e3` renders meters as km in the plotting helpers.
        %
        % - Topic: Inspect model settings
        visualScale (1,1) double {mustBePositive} = 1e3

        % Descriptive model name used in figures and generated docs.
        %
        % - Topic: Inspect model settings
        name = ''
    end

    properties (Access = private)
        uFunction = []
        vFunction = []
    end

    methods
        function self = KinematicModel()
            % Create an empty kinematic-model base instance.
            %
            % Configure velocity closures and domain properties after
            % construction when deriving custom models from `KinematicModel`.
            %
            % - Topic: Configure model domains
            % - Declaration: self = KinematicModel()
            % - Returns self: `KinematicModel` instance
        end

        function uValue = u(self, t, x, y)
            % Evaluate the x-velocity component.
            %
            % Output shape matches the input `x` and `y` arrays.
            %
            % - Topic: Evaluate velocity fields
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: x-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`
            if isempty(self.uFunction)
                uValue = zeros(size(x));
            else
                uValue = self.uFunction(t, x, y);
            end
        end

        function vValue = v(self, t, x, y)
            % Evaluate the y-velocity component.
            %
            % Output shape matches the input `x` and `y` arrays.
            %
            % - Topic: Evaluate velocity fields
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar evaluation time in seconds
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: y-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`
            if isempty(self.vFunction)
                vValue = zeros(size(x));
            else
                vValue = self.vFunction(t, x, y);
            end
        end
    end

    methods (Static)
        function bounds = boundsFromLimits(xlim, ylim)
            % Convert rectangular limits to polygon vertices.
            %
            % Use this helper when a rectangular plotting or integration
            % region should be represented as a closed polygon.
            %
            % - Topic: Work with particle domains
            % - Declaration: bounds = boundsFromLimits(xlim,ylim)
            % - Parameter xlim: two-element vector of x limits in meters
            % - Parameter ylim: two-element vector of y limits in meters
            % - Returns bounds: structure with fields `xv` and `yv`
            arguments
                xlim (1,2) double
                ylim (1,2) double
            end

            xv = [min(xlim) min(xlim) max(xlim) max(xlim)];
            yv = [max(ylim) min(ylim) min(ylim) max(ylim)];
            bounds = struct('xv', xv, 'yv', yv);
        end
    end
end
