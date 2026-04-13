classdef AdvectionDiffusionIntegrator
    % Integrate particles through a kinematic flow with diffusivity.
    %
    % This class advances particle positions according to
    %
    % $$ d\mathbf{r}(t) = \mathbf{u}(t,\mathbf{r})\,dt + \sqrt{2\kappa}\,d\mathbf{W}_t, $$
    %
    % where $$\mathbf{r}(t) = [x(t)\ y(t)]$$,
    % $$\mathbf{u}(t,\mathbf{r})$$ is supplied by a `KinematicModel`, and
    % $$\kappa$$ is the scalar tracer diffusivity in $$m^2 s^{-1}$$. Boundary
    % limits, periodic wrapping, and polygonal obstacles are delegated to
    % `IntegratorWithObstacles`.
    %
    % The output trajectory arrays follow the contract
    %
    % - `size(t) = [nTimes 1]`
    % - `size(x) = [nTimes nParticles]`
    % - `size(y) = [nTimes nParticles]`
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
    % - Topic: Create particle integrators
    % - Topic: Inspect integrator settings
    % - Topic: Integrate particle trajectories
    % - Declaration: classdef AdvectionDiffusionIntegrator

    properties
        % Underlying deterministic velocity model.
        %
        % `kinematicModel` supplies the advecting velocity field through
        % `u(t,x,y)` and `v(t,x,y)`.
        %
        % - Topic: Inspect integrator settings
        kinematicModel (1,1) KinematicModel

        % Scalar diffusivity in $$m^2 s^{-1}$$.
        %
        % Set `kappa = 0` to recover purely advective trajectories.
        %
        % - Topic: Inspect integrator settings
        kappa (1,1) double {mustBeNonnegative} = 0

        % Internal integrator step size in seconds.
        %
        % A value of `0` uses the requested output increment `dt` as the
        % integration step.
        %
        % - Topic: Inspect integrator settings
        stepSize (1,1) double {mustBeNonnegative} = 0
    end

    methods
        function self = AdvectionDiffusionIntegrator(kinematicModel, kappa)
            % Create an advection-diffusion integrator.
            %
            % `kinematicModel` supplies the deterministic velocity field,
            % and `kappa` sets the scalar tracer diffusivity.
            %
            % - Topic: Create particle integrators
            % - Declaration: self = AdvectionDiffusionIntegrator(kinematicModel,kappa)
            % - Parameter kinematicModel: `KinematicModel` instance defining `u(t,x,y)` and `v(t,x,y)`
            % - Parameter kappa: scalar diffusivity in $$m^2 s^{-1}$$
            % - Returns self: `AdvectionDiffusionIntegrator` instance
            arguments
                kinematicModel (1,1) KinematicModel
                kappa (1,1) double {mustBeNonnegative}
            end

            self.kinematicModel = kinematicModel;
            self.kappa = kappa;
        end
    end
end
