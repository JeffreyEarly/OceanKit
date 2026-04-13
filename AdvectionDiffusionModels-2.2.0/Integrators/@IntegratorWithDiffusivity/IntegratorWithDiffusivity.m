classdef IntegratorWithDiffusivity < Integrator
    % Integrate $$dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t$$ in a box domain.
    %
    % `IntegratorWithDiffusivity` combines the deterministic RK4 drift
    % update from `Integrator` with additive stochastic diffusion,
    %
    % $$
    % dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t,
    % $$
    %
    % where `kappa` may be a scalar or a `1 x nDims` vector of componentwise
    % diffusivities $$\kappa_i$$. The bounds `ymin` and `ymax` define
    % independent reflecting or unbounded box conditions along each
    % coordinate, using the repository's existing reflected and wrapped
    % increment formulas.
    %
    % ```matlab
    % f = @(t, y) zeros(size(y));
    % y0 = [0.25 0.25; 0.75 0.75];
    % integrator = IntegratorWithDiffusivity(f, y0, dt=0.1, kappa=0.01, ymin=0, ymax=1);
    % y = integrator.advanceOneStep();
    % ```
    %
    % - Topic: Create the integrator
    % - Topic: Inspect integrator settings
    % - Declaration: classdef IntegratorWithDiffusivity < Integrator

    properties (SetAccess = protected)
        % Componentwise diffusivity `kappa`.
        %
        % `kappa` stores either the scalar diffusivity expanded across all
        % dimensions or the supplied `1 x nDims` diffusivity vector.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.kappa
        % - Returns kappa: scalar-expanded or componentwise diffusivity $$\kappa$$ with units of `y.^2 / t`
        kappa double = []

        % Lower box bounds `ymin`.
        %
        % `ymin` stores the lower coordinate bound for each state
        % dimension. Use `-Inf` for an unbounded lower side.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.ymin
        % - Returns ymin: row vector of lower coordinate bounds
        ymin double = []

        % Upper box bounds `ymax`.
        %
        % `ymax` stores the upper coordinate bound for each state
        % dimension. Use `Inf` for an unbounded upper side.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.ymax
        % - Returns ymax: row vector of upper coordinate bounds
        ymax double = []
    end

    properties (Access = protected)
        diffusivityFlux cell = {}
    end

    methods
        function self = IntegratorWithDiffusivity(f, y0, options)
            % Create an additive-diffusion integrator on a box domain.
            %
            % The constructor keeps the legacy scalar-expansion behavior
            % for `kappa`, `ymin`, and `ymax`, with defaults corresponding
            % to zero diffusivity on an unbounded domain.
            %
            % - Topic: Create the integrator
            % - Declaration: self = IntegratorWithDiffusivity(f,y0,dt=...,kappa=...,ymin=...,ymax=...)
            % - Parameter f: deterministic drift function $$f(t,y)$$
            % - Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
            % - Parameter dt: positive scalar timestep $$dt$$
            % - Parameter kappa: scalar or `1 x nDims` diffusivity vector $$\kappa$$; defaults to `0`
            % - Parameter ymin: scalar or `1 x nDims` lower box bounds; defaults to `-Inf`
            % - Parameter ymax: scalar or `1 x nDims` upper box bounds; defaults to `Inf`
            arguments
                f (1,1) function_handle
                y0 {mustBeNumeric, mustBeReal}
                options.dt (1,1) double {mustBePositive}
                options.kappa {mustBeNumeric, mustBeReal, mustBeNonnegative} = 0
                options.ymin {mustBeNumeric, mustBeReal} = -Inf
                options.ymax {mustBeNumeric, mustBeReal} = Inf
            end

            [kappa, ymin, ymax] = normalizeDiffusivityInputs(y0, options.kappa, options.ymin, options.ymax);

            self@Integrator(f, y0, dt=options.dt);

            self.kappa = kappa;
            self.ymin = ymin;
            self.ymax = ymax;
            self.diffusivityFlux = cell(self.nDims, 1);

            for iDim = 1:self.nDims
                if kappa(iDim) == 0
                    self.diffusivityFlux{iDim} = @(x, dt) zeros(size(x));
                elseif ymin(iDim) == -Inf && isfinite(ymax(iDim))
                    b = ymax(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) -abs((s/sqrt(dt))*randn(size(x)) - (b - x)/dt) + (b - x)/dt;
                elseif isfinite(ymin(iDim)) && ymax(iDim) == Inf
                    a = ymin(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) abs((s/sqrt(dt))*randn(size(x)) - (a - x)/dt) + (a - x)/dt;
                elseif isfinite(ymin(iDim)) && isfinite(ymax(iDim))
                    a = ymin(iDim);
                    b = ymax(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) abs(mod((s/sqrt(dt))*randn(size(x)) - (b - x)/dt, 2*(b - a)/dt) - (b - a)/dt) + (a - x)/dt;
                elseif ymin(iDim) == -Inf && ymax(iDim) == Inf
                    self.diffusivityFlux{iDim} = @(x, dt) sqrt(2*kappa(iDim)/dt)*randn(size(x));
                else
                    throwInputError("IntegratorWithDiffusivity:InvalidBoundarySpecification", "Invalid boundary specification for dimension %d.", iDim);
                end
            end
        end
    end

    methods (Access = protected)
        function y = stepForward(self, y, t, dt)
            y = stepForward@Integrator(self, y, t, dt);
            for iDim = 1:self.nDims
                y(:, iDim) = y(:, iDim) + dt*self.diffusivityFlux{iDim}(y(:, iDim), dt);
            end
        end
    end
end

function [kappa, ymin, ymax] = normalizeDiffusivityInputs(y0, kappa, ymin, ymax)
[~, nDims] = size(y0);

kappa = expandPerDimensionInput(kappa, nDims, "kappa", "IntegratorWithDiffusivity:InvalidKappaSize");
ymin = expandPerDimensionInput(ymin, nDims, "ymin", "IntegratorWithDiffusivity:InvalidYminSize");
ymax = expandPerDimensionInput(ymax, nDims, "ymax", "IntegratorWithDiffusivity:InvalidYmaxSize");
end

function value = expandPerDimensionInput(value, nDims, name, identifier)
if isscalar(value)
    value = repmat(value, 1, nDims);
elseif isvector(value) && numel(value) == nDims
    value = reshape(value, 1, []);
else
    throwInputError(identifier, "%s must be a scalar or have %d entries.", name, nDims);
end
end

function throwInputError(identifier, message, varargin)
throwAsCaller(MException(identifier, message, varargin{:}));
end
