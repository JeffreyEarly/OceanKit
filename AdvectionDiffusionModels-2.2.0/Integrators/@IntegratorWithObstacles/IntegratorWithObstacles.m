classdef IntegratorWithObstacles < Integrator
    % Integrate diffusion with reflecting polygonal obstacles in 2D.
    %
    % `IntegratorWithObstacles` advances the stochastic model
    %
    % $$
    % dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t
    % $$
    %
    % on a two-dimensional domain with coordinate bounds `ymin` and `ymax`,
    % optional periodic wrapping, and polygonal obstacles. Proposed step
    % segments that cross an obstacle are handled by repeated geometric
    % reflection against the obstacle boundary.
    %
    % ```matlab
    % f = @(t, y) 0.1 * ones(size(y));
    % y0 = [0.2 0.2; 0.2 0.8];
    % obstacles = polyshape([0.45 0.55 0.55 0.45], [0.35 0.35 0.65 0.65]);
    % integrator = IntegratorWithObstacles(f, y0, dt=0.1, kappa=0.01, ymin=0, ymax=1, obstacles=obstacles);
    % y = integrator.advanceOneStep();
    % ```
    %
    % - Topic: Create the integrator
    % - Topic: Inspect integrator settings
    % - Topic: Handle obstacle reflections
    % - Declaration: classdef IntegratorWithObstacles < Integrator

    properties (SetAccess = protected)
        % Componentwise diffusivity `kappa`.
        %
        % `kappa` stores either the scalar diffusivity expanded across both
        % state dimensions or the supplied `1 x 2` diffusivity vector.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.kappa
        % - Returns kappa: scalar-expanded or componentwise diffusivity $$\kappa$$ with units of `y.^2 / t`
        kappa double = []

        % Lower coordinate bounds `ymin`.
        %
        % `ymin` stores the lower bound for each coordinate. Use `-Inf` to
        % indicate an unbounded lower side.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.ymin
        % - Returns ymin: `1 x 2` row vector of lower coordinate bounds
        ymin double = []

        % Upper coordinate bounds `ymax`.
        %
        % `ymax` stores the upper bound for each coordinate. Use `Inf` to
        % indicate an unbounded upper side.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.ymax
        % - Returns ymax: `1 x 2` row vector of upper coordinate bounds
        ymax double = []

        % Periodicity flags for each coordinate direction.
        %
        % `isPeriodic` is normalized to a logical `1 x 2` row vector. A
        % true entry means the corresponding coordinate is wrapped into the
        % interval `[ymin(i), ymax(i))` before obstacle checks.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.isPeriodic
        % - Returns isPeriodic: logical `1 x 2` vector of periodic-direction flags
        isPeriodic (1,2) logical = [false false]

        % Polygonal obstacles used for reflecting boundaries.
        %
        % `obstacles` stores the original non-overlapping `polyshape`
        % objects supplied to the constructor.
        %
        % - Topic: Inspect integrator settings
        % - Declaration: self.obstacles
        % - Returns obstacles: array of reflecting `polyshape` obstacles
        obstacles polyshape = polyshape.empty(0,1)
    end

    properties (Access = protected)
        diffusivityFlux cell = {}
        obstacleData struct = struct.empty(0,1)
    end

    methods
        function self = IntegratorWithObstacles(f, y0, options)
            % Create a 2D obstacle-aware stochastic integrator.
            %
            % The constructor validates the 2D state layout, expands
            % scalar-valued `kappa`, `ymin`, `ymax`, and `isPeriodic`,
            % stores the domain bounds, and precomputes polygon and
            % bounding-box geometry used by the reflection algorithm.
            %
            % - Topic: Create the integrator
            % - Declaration: self = IntegratorWithObstacles(f,y0,dt=...,kappa=...,ymin=...,ymax=...,obstacles=polyshape.empty(0,1),isPeriodic=false)
            % - Parameter f: deterministic drift function $$f(t,y)$$
            % - Parameter y0: initial condition $$y_0$$ stored as `nParticles x 2`
            % - Parameter dt: positive scalar timestep $$dt$$
            % - Parameter kappa: scalar or `1 x 2` diffusivity vector $$\kappa$$; defaults to `0`
            % - Parameter ymin: scalar or `1 x 2` lower coordinate bounds; defaults to `-Inf`
            % - Parameter ymax: scalar or `1 x 2` upper coordinate bounds; defaults to `Inf`
            % - Parameter obstacles: non-overlapping reflecting `polyshape` obstacles
            % - Parameter isPeriodic: logical scalar or `1 x 2` periodicity flags; defaults to `false`
            arguments
                f (1,1) function_handle
                y0 {mustBeNumeric, mustBeReal}
                options.dt (1,1) double {mustBePositive}
                options.kappa {mustBeNumeric, mustBeReal, mustBeNonnegative} = 0
                options.ymin {mustBeNumeric, mustBeReal} = -Inf
                options.ymax {mustBeNumeric, mustBeReal} = Inf
                options.obstacles = polyshape.empty(0,1)
                options.isPeriodic = false
            end

            [~, nDims] = size(y0);
            if nDims ~= 2
                throwInputError("IntegratorWithObstacles:InvalidStateDimension", "y0 must have exactly 2 columns, but it has %d.", nDims);
            end

            [kappa, ymin, ymax] = normalizeObstacleInputs(options.kappa, options.ymin, options.ymax, nDims);
            obstacles = normalizeObstacles(options.obstacles);
            isPeriodic = normalizePeriodicity(options.isPeriodic, nDims);

            if ~isempty(obstacles) && any(any(overlaps(obstacles) - eye(length(obstacles))))
                throwInputError("IntegratorWithObstacles:OverlappingObstacles", "Obstacles must not overlap. Use union to combine overlapping polygons before construction.");
            end

            self@Integrator(f, y0, dt=options.dt);

            self.kappa = kappa;
            self.ymin = ymin;
            self.ymax = ymax;
            self.isPeriodic = isPeriodic;
            self.obstacles = obstacles;

            self.setObstacles(obstacles);
            for iObstacle = 1:length(obstacles)
                [isInside, isOnBoundary] = inpolygon(y0(:,1), y0(:,2), obstacles(iObstacle).Vertices(:,1), obstacles(iObstacle).Vertices(:,2));
                if any(isInside & ~isOnBoundary)
                    throwInputError("IntegratorWithObstacles:InitialStateInsideObstacle", "Initial particles inside obstacle are not allowed.");
                end
            end

            self.diffusivityFlux = cell(self.nDims, 1);
            for iDim = 1:self.nDims
                if kappa(iDim) == 0
                    self.diffusivityFlux{iDim} = @(x, dt) zeros(size(x));
                elseif ymin(iDim) == -Inf && isfinite(ymax(iDim)) && ~isPeriodic(iDim)
                    b = ymax(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) -abs((s/sqrt(dt))*randn(size(x)) - (b - x)/dt) + (b - x)/dt;
                elseif isfinite(ymin(iDim)) && ymax(iDim) == Inf && ~isPeriodic(iDim)
                    a = ymin(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) abs((s/sqrt(dt))*randn(size(x)) - (a - x)/dt) + (a - x)/dt;
                elseif isfinite(ymin(iDim)) && isfinite(ymax(iDim)) && ~isPeriodic(iDim)
                    a = ymin(iDim);
                    b = ymax(iDim);
                    s = sqrt(2*kappa(iDim));
                    self.diffusivityFlux{iDim} = @(x, dt) abs(mod((s/sqrt(dt))*randn(size(x)) - (b - x)/dt, 2*(b - a)/dt) - (b - a)/dt) + (a - x)/dt;
                else
                    self.diffusivityFlux{iDim} = @(x, dt) sqrt(2*kappa(iDim)/dt)*randn(size(x));
                end
            end
        end
    end

    methods (Access = protected)
        y = stepForward(self, y, t, dt)
        setObstacles(self, obstacles)
        [y, dy, didReflect] = updateWithReflection(self, obstacle, y, yWrapped, dy)
    end

    methods (Static, Access = protected)
        [yFinal, yIntersection, didReflect] = reflect(obstacle, yInitial, yFinal)
        isLeft = isLeft(obstacle, iEdge, x, y)
        isInterior = isInterior(obstacle, x, y)
        doesIntersect = doesIntersectPolygon(obstacle, yInitial, yFinal)
        [doesIntersect, yIntersection, r2] = doesIntersectSegment(yInitial, yFinal, vInitial, vFinal, doesIntersect, yIntersection, r2)
    end
end

function [kappa, ymin, ymax] = normalizeObstacleInputs(kappa, ymin, ymax, nDims)
kappa = expandPerDimensionInput(kappa, nDims, "kappa", "IntegratorWithObstacles:InvalidKappaSize");
ymin = expandPerDimensionInput(ymin, nDims, "ymin", "IntegratorWithObstacles:InvalidYminSize");
ymax = expandPerDimensionInput(ymax, nDims, "ymax", "IntegratorWithObstacles:InvalidYmaxSize");
end

function isPeriodic = normalizePeriodicity(isPeriodic, nDims)
if ~islogical(isPeriodic)
    throwInputError("IntegratorWithObstacles:InvalidPeriodicityType", "isPeriodic must be logical.");
end

isPeriodic = expandPerDimensionInput(isPeriodic, nDims, "isPeriodic", "IntegratorWithObstacles:InvalidPeriodicitySize");
end

function obstacles = normalizeObstacles(obstacles)
if isempty(obstacles)
    obstacles = polyshape.empty(0,1);
    return
end

if ~isa(obstacles, "polyshape")
    throwInputError("IntegratorWithObstacles:InvalidObstaclesType", "obstacles must be empty or a polyshape array.");
end
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
