classdef GriddedStreamfunction < CAAnnotatedClass
    % Fit a COM-frame streamfunction estimator and trajectory decomposition.
    %
    % `GriddedStreamfunction` fits a mesoscale streamfunction
    % $$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
    % $$m_x(t), m_y(t)$$, and an anchored background trajectory
    % $$x^{\mathrm{bg}}(t), y^{\mathrm{bg}}(t)$$ from asynchronous drifter
    % trajectory splines using the coupled estimator described in
    % `asynchronous-com-fit-rev6.tex`.
    %
    % The estimator uses the fast temporal spline basis both to define
    % the COM smoothing operator and to represent the common background
    % velocity. It first fits the COM trajectory, then solves one
    % least-squares problem for the gauge-reduced mesoscale
    % streamfunction coefficients in the COM frame, and finally recovers
    % the background velocity from the COM-space residual in that same
    % fast basis. Optional hard `zeroVorticity` and `zeroStrain`
    % constraints are applied directly to the mesoscale coefficients.
    %
    % The fitted decomposition follows
    %
    % $$
    % \dot{x} = u^{\mathrm{meso}} + u^{\mathrm{bg}} + u^{\mathrm{sm}},
    % \qquad
    % \dot{y} = v^{\mathrm{meso}} + v^{\mathrm{bg}} + v^{\mathrm{sm}},
    % $$
    %
    % with centered coordinates
    % $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$
    %
    % The fixed-frame trajectory decomposition is stored as per-drifter
    % spline trajectories satisfying
    %
    % $$
    % x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
    % \qquad
    % y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}},
    % $$
    %
    % where the mesoscale trajectory carries the observed initial drifter
    % position and the background and submesoscale trajectories are
    % zero-anchored at the first drifter sample. In the centered frame,
    % the mesoscale trajectory carries the initial centered position and
    % the centered submesoscale trajectory is zero-anchored.
    %
    % ```matlab
    % fit = GriddedStreamfunction.fromTrajectories(trajectories);
    % tFit = fit.fitSupportTimes;
    % com = fit.centerOfMassTrajectory;
    % background = fit.backgroundTrajectory;
    % decomposition = fit.decomposition;
    % uMeso = fit.uMesoscale(tFit, com.x(tFit), com.y(tFit));
    % ```
    %
    % - Topic: Fit the estimator
    % - Topic: Persist and restart fits
    % - Topic: Inspect fit setup and structure
    % - Topic: Inspect fit setup and structure — Input trajectories
    % - Topic: Inspect fit setup and structure — Mesoscale basis
    % - Topic: Inspect fit setup and structure — Fast temporal basis
    % - Topic: Inspect primary outputs
    % - Topic: Apply fitted decomposition
    % - Topic: Evaluate derived fields
    % - Topic: Evaluate derived fields — Coordinate transform
    % - Topic: Evaluate derived fields — Mesoscale field evaluation
    % - Topic: Evaluate derived fields — Background evaluation
    % - Topic: Evaluate derived fields — Derived diagnostics
    % - Topic: Evaluate derived fields — Visualization helper
    % - Declaration: classdef GriddedStreamfunction < CAAnnotatedClass

    properties (SetAccess = private)
        % Fitted COM-frame mesoscale streamfunction spline.
        %
        % This spline evaluates $$\psi(\tilde{x},\tilde{y},t)$$ in the
        % centered coordinates
        % $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$ It is
        % the solved mesoscale basis state; methods such as
        % `psiMesoscale`, `uMesoscale`, and `vMesoscale` are derived
        % evaluations of this spline rather than separately fitted fields.
        %
        % - Topic: Inspect fit setup and structure — Mesoscale basis
        % - nav_order: 1
        streamfunctionSpline

        % Observed drifter trajectory splines used for the fit.
        %
        % `observedTrajectories` preserves the original
        % `TrajectorySpline` inputs in drifter order.
        %
        % - Topic: Inspect fit setup and structure — Input trajectories
        % - nav_order: 1
        observedTrajectories

        % Fitted center-of-mass trajectory.
        %
        % This is one of the estimator's primary solved outputs.
        % `centerOfMassTrajectory.x(t)` evaluates $$m_x(t)$$ and
        % `centerOfMassTrajectory.y(t)` evaluates $$m_y(t)$$ on the fit
        % support interval.
        %
        % ```matlab
        % tFit = fit.fitSupportTimes;
        % plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
        % axis equal
        % xlabel("x (m)")
        % ylabel("y (m)")
        % ```
        %
        % - Topic: Inspect primary outputs
        % - nav_order: 1
        centerOfMassTrajectory

        % Fitted common background trajectory.
        %
        % This is one of the estimator's primary solved outputs and stores
        % the single common anchored background path recovered by the fit.
        % `backgroundTrajectory.x(t)` evaluates $$x^{\mathrm{bg}}(t)$$
        % and `backgroundTrajectory.y(t)` evaluates
        % $$y^{\mathrm{bg}}(t)$$, with
        % $$x^{\mathrm{bg}}(t_0)=0$$ and $$y^{\mathrm{bg}}(t_0)=0$$ at
        % the global fit start time. The recovered background velocity is
        % obtained from `backgroundTrajectory.u(t)` and
        % `backgroundTrajectory.v(t)`. This is the single shared
        % background path for the fitted estimator; for drifter `k`,
        % `decomposition.fixedFrame.background(k)` is the same path
        % re-anchored so it starts at zero at that drifter's first sample
        % time.
        %
        % ```matlab
        % tFit = fit.fitSupportTimes;
        % plot(fit.backgroundTrajectory.x(tFit), fit.backgroundTrajectory.y(tFit))
        % axis equal
        % xlabel("x^{bg} (m)")
        % ylabel("y^{bg} (m)")
        % ```
        %
        % - Topic: Inspect primary outputs
        % - nav_order: 2
        backgroundTrajectory

        % Hard constraint applied to the fitted mesoscale streamfunction.
        %
        % `mesoscaleConstraint` is `"none"`, `"zeroVorticity"`, or
        % `"zeroStrain"`.
        %
        % - Topic: Inspect fit setup and structure — Mesoscale basis
        % - nav_order: 4
        mesoscaleConstraint string = "none"

        % Representative pooled times from the stride-rule fast basis.
        %
        % This is empty when `fastKnotPoints` are supplied directly.
        %
        % - Topic: Inspect fit setup and structure — Fast temporal basis
        % - nav_order: 4
        representativeTimes

        % Sorted unique observation times used as trajectory support.
        %
        % `fitSupportTimes` is the canonical solved support grid shared by
        % the fitted COM and background trajectories.
        %
        % - Topic: Inspect fit setup and structure — Fast temporal basis
        % - nav_order: 3
        fitSupportTimes

        % Identifiable mesoscale degrees of freedom after gauge and constraints.
        %
        % `mesoscaleDegreesOfFreedom` counts the solved mesoscale spline
        % degrees of freedom after removing the additive streamfunction
        % gauge and applying the selected hard mesoscale constraint. This
        % scalar is determined by the resolved mesoscale basis and the
        % chosen structural constraint, not by the fitted coefficient
        % values themselves.
        %
        % - Topic: Inspect fit setup and structure — Mesoscale basis
        % - nav_order: 5
        mesoscaleDegreesOfFreedom (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    end

    properties (Dependent)
        % Per-drifter decomposition trajectories in fixed and centered frames.
        %
        % This is one of the estimator's primary solved outputs. Use
        % `fit.decomposition` to inspect how the fitted estimator
        % reconstructs the same drifters that were used in the fit. Each
        % field is a `TrajectorySpline` column vector aligned one-for-one
        % with `observedTrajectories`, so
        % `decomposition.fixedFrame.mesoscale(k)` is the mesoscale
        % trajectory for `observedTrajectories(k)`. Call
        % `decomposeTrajectories` when the same fitted COM and mesoscale
        % fields should be applied to a different set of drifters.
        %
        % In the fixed frame,
        % `decomposition.fixedFrame.background`,
        % `decomposition.fixedFrame.mesoscale`, and
        % `decomposition.fixedFrame.submesoscale` satisfy
        %
        % $$
        % x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
        % \qquad
        % y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}}.
        % $$
        %
        % The canonical shared background path itself is stored in
        % `backgroundTrajectory`. The fixed-frame background component is
        % that same path re-anchored so that
        % $$x_k^{\mathrm{bg}}(t_{k,0}) = y_k^{\mathrm{bg}}(t_{k,0}) = 0$$
        % at the first sample time $$t_{k,0}$$ of drifter $$k$$. The
        % fixed-frame mesoscale trajectory carries the observed initial
        % drifter position, and the fixed-frame submesoscale residual is
        % zero-anchored.
        %
        % In the centered frame,
        % `decomposition.centeredFrame.mesoscale` and
        % `decomposition.centeredFrame.submesoscale` satisfy
        %
        % $$
        % \tilde{x}_k = \tilde{x}_k^{\mathrm{meso}} + \tilde{x}_k^{\mathrm{sm}},
        % \qquad
        % \tilde{y}_k = \tilde{y}_k^{\mathrm{meso}} + \tilde{y}_k^{\mathrm{sm}},
        % $$
        %
        % after subtracting `centerOfMassTrajectory` from the observed
        % drifter path. Use these centered-frame splines when you want to
        % compare the fitted mesoscale motion against the observed
        % COM-relative motion.
        %
        % ```matlab
        % iDrifter = 1;
        % trajectory = fit.observedTrajectories(iDrifter);
        % ti = trajectory.t;
        %
        % background = fit.decomposition.fixedFrame.background(iDrifter);
        % mesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter);
        % submesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter);
        %
        % plot(trajectory.x(ti), trajectory.y(ti), "k")
        % hold on
        % plot(background.x(ti), background.y(ti))
        % plot(mesoscale.x(ti), mesoscale.y(ti))
        % plot(submesoscale.x(ti), submesoscale.y(ti))
        % axis equal
        % ```
        %
        % - Topic: Inspect primary outputs
        % - nav_order: 3
        decomposition

        % Fast temporal knot vector used for COM and background fits.
        %
        % - Topic: Inspect fit setup and structure — Fast temporal basis
        % - nav_order: 1
        fastKnotPoints

        % Fast temporal spline degree for COM and background fits.
        %
        % - Topic: Inspect fit setup and structure — Fast temporal basis
        % - nav_order: 2
        fastS

        % Mesoscale tensor-product knot vectors `{qKnot, rKnot, tKnot}`.
        %
        % - Topic: Inspect fit setup and structure — Mesoscale basis
        % - nav_order: 2
        psiKnotPoints

        % Mesoscale spline degrees `[Sq Sr St]`.
        %
        % - Topic: Inspect fit setup and structure — Mesoscale basis
        % - nav_order: 3
        psiS
    end

    properties (Hidden, SetAccess = private)
        fixedFrameBackgroundTrajectories = TrajectorySpline.empty(0, 1)
        fixedFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        fixedFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        centeredFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        centeredFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
    end

    properties (Access = private)
        observedTrajectorySampleData struct = struct([])
    end

    methods
        function self = GriddedStreamfunction(options)
            % Create a fit from canonical solved-state properties.
            %
            % Use this low-level constructor when the fitted spline state
            % is already available, for example after reading a persisted
            % restart file. For fitting from observed drifter trajectories,
            % use `GriddedStreamfunction.fromTrajectories(...)`.
            %
            % ```matlab
            % fit = GriddedStreamfunction( ...
            %     streamfunctionSpline=psiSpline, ...
            %     observedTrajectories=trajectories, ...
            %     centerOfMassTrajectory=centerTrajectory, ...
            %     backgroundTrajectory=backgroundTrajectory, ...
            %     fixedFrameBackgroundTrajectories=backgroundTrajectories, ...
            %     fixedFrameMesoscaleTrajectories=mesoscaleTrajectories, ...
            %     fixedFrameSubmesoscaleTrajectories=submesoscaleTrajectories, ...
            %     centeredFrameMesoscaleTrajectories=centeredMesoscaleTrajectories, ...
            %     centeredFrameSubmesoscaleTrajectories=centeredSubmesoscaleTrajectories, ...
            %     mesoscaleConstraint="none", ...
            %     representativeTimes=representativeTimes, ...
            %     fitSupportTimes=fitSupportTimes);
            % ```
            %
            % - Topic: Fit the estimator
            % - Declaration: self = GriddedStreamfunction(options)
            % - Parameter options.streamfunctionSpline: fitted centered-frame mesoscale streamfunction spline
            % - Parameter options.observedTrajectories: observed drifter trajectory splines aligned with the fit
            % - Parameter options.centerOfMassTrajectory: fitted center-of-mass trajectory
            % - Parameter options.backgroundTrajectory: fitted common background trajectory
            % - Parameter options.fixedFrameBackgroundTrajectories: stored fixed-frame background decomposition trajectories
            % - Parameter options.fixedFrameMesoscaleTrajectories: stored fixed-frame mesoscale decomposition trajectories
            % - Parameter options.fixedFrameSubmesoscaleTrajectories: stored fixed-frame submesoscale decomposition trajectories
            % - Parameter options.centeredFrameMesoscaleTrajectories: stored centered-frame mesoscale decomposition trajectories
            % - Parameter options.centeredFrameSubmesoscaleTrajectories: stored centered-frame submesoscale decomposition trajectories
            % - Parameter options.mesoscaleConstraint: hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Parameter options.representativeTimes: pooled representative times from the stride-rule fast basis
            % - Parameter options.fitSupportTimes: sorted unique observation times used as trajectory support
            % - Parameter options.mesoscaleDegreesOfFreedom: stored identifiable mesoscale degrees of freedom after gauge reduction and hard constraints
            % - Returns self: canonical `GriddedStreamfunction` instance
            arguments
                options.streamfunctionSpline = TensorSpline.empty(0, 1)
                options.observedTrajectories = TrajectorySpline.empty(0, 1)
                options.centerOfMassTrajectory = TrajectorySpline.empty(0, 1)
                options.backgroundTrajectory = TrajectorySpline.empty(0, 1)
                options.fixedFrameBackgroundTrajectories = TrajectorySpline.empty(0, 1)
                options.fixedFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.fixedFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.centeredFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.centeredFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
                options.representativeTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fitSupportTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.mesoscaleDegreesOfFreedom (1,1) double {mustBeInteger,mustBeNonnegative} = 0
            end

            self@CAAnnotatedClass();

            if nargin == 0
                return
            end

            GriddedStreamfunction.validateCanonicalState(options);
            self.streamfunctionSpline = options.streamfunctionSpline;
            self.observedTrajectories = reshape(options.observedTrajectories, [], 1);
            self.centerOfMassTrajectory = options.centerOfMassTrajectory;
            self.backgroundTrajectory = reshape(options.backgroundTrajectory, [], 1);
            self.fixedFrameBackgroundTrajectories = reshape(options.fixedFrameBackgroundTrajectories, [], 1);
            self.fixedFrameMesoscaleTrajectories = reshape(options.fixedFrameMesoscaleTrajectories, [], 1);
            self.fixedFrameSubmesoscaleTrajectories = reshape(options.fixedFrameSubmesoscaleTrajectories, [], 1);
            self.centeredFrameMesoscaleTrajectories = reshape(options.centeredFrameMesoscaleTrajectories, [], 1);
            self.centeredFrameSubmesoscaleTrajectories = reshape(options.centeredFrameSubmesoscaleTrajectories, [], 1);
            self.mesoscaleConstraint = string(options.mesoscaleConstraint);
            self.representativeTimes = options.representativeTimes;
            self.fitSupportTimes = options.fitSupportTimes;
            self.mesoscaleDegreesOfFreedom = options.mesoscaleDegreesOfFreedom;
            self.refreshObservedTrajectorySampleData();
        end
    end

    methods (Static)
        function self = fromFile(path)
            % Read a fitted estimator from a NetCDF restart file.
            %
            % `fromFile` reconstructs the canonical solved state written by
            % `writeToFile` without rerunning the trajectory fit.
            %
            % ```matlab
            % fit = GriddedStreamfunction.fromFile("gridded-fit.nc");
            % decomposition = fit.decomposition;
            % tFit = fit.fitSupportTimes;
            % plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
            % axis equal
            % ```
            %
            % - Topic: Persist and restart fits
            % - nav_order: 1
            % - Declaration: self = fromFile(path)
            % - Parameter path: path to the NetCDF restart file
            % - Returns self: reconstructed `GriddedStreamfunction` estimator
            arguments (Input)
                path {mustBeTextScalar}
            end

            self = griddedStreamfunctionFromFile(char(path));
        end

        function self = fromTrajectories(trajectories, options)
            % Fit the estimator from drifter trajectory splines.
            %
            % Use this factory with one `TrajectorySpline` per drifter.
            % Positions are sampled with `trajectory.x(trajectory.t)` and
            % `trajectory.y(trajectory.t)`, and observed velocities are the
            % first derivatives of the same trajectory splines.
            %
            % The fast temporal basis is used for both the center-of-mass
            % trajectory and the recovered common background path. The
            % tensor basis defined by `psiS` and `psiKnotPoints` is used
            % only for the mesoscale streamfunction in the centered frame,
            % with only the additive streamfunction gauge removed. Set
            % `mesoscaleConstraint` to impose a hard zero-vorticity or
            % zero-strain mesoscale fit.
            %
            % ```matlab
            % fit = GriddedStreamfunction.fromTrajectories( ...
            %     trajectories, mesoscaleConstraint="zeroVorticity");
            % tFit = fit.fitSupportTimes;
            % plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
            % hold on
            % plot(fit.backgroundTrajectory.x(tFit), fit.backgroundTrajectory.y(tFit))
            % axis equal
            % ```
            %
            % - Topic: Fit the estimator
            % - nav_order: 1
            % - Declaration: self = fromTrajectories(trajectories,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter fastS: optional fast temporal spline degree, default `3`
            % - Parameter mesoscaleConstraint: optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Returns self: fitted `GriddedStreamfunction` estimator
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
                options.buildDecomposition (1,1) logical = true
                options.sampleData struct = struct([])
            end

            trajectories = reshape(trajectories, [], 1);
            psiS = reshape(options.psiS, 1, []);
            psiKnotPoints = options.psiKnotPoints;
            fastKnotPoints = options.fastKnotPoints;
            fastS = options.fastS;
            mesoscaleConstraint = string(options.mesoscaleConstraint);
            sampleData = options.sampleData;

            if ~isempty(fastKnotPoints)
                if ~(isnumeric(fastKnotPoints) && isvector(fastKnotPoints) && ...
                        isreal(fastKnotPoints) && all(isfinite(fastKnotPoints)))
                    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                        "fastKnotPoints must be a finite real vector.");
                end

                fastKnotPoints = reshape(fastKnotPoints, [], 1);
                if any(diff(fastKnotPoints) < 0)
                    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                        "fastKnotPoints must be nondecreasing.");
                end
            end

            if ~isempty(psiKnotPoints)
                if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
                    error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                        "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
                end

                psiKnotPoints = reshape(psiKnotPoints, 1, []);
                for iDim = 1:3
                    knotVector = psiKnotPoints{iDim};
                    if ~(isnumeric(knotVector) && isvector(knotVector) && ...
                            isreal(knotVector) && all(isfinite(knotVector)))
                        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                            "Each psi knot vector must be a finite real vector.");
                    end

                    knotVector = reshape(knotVector, [], 1);
                    if any(diff(knotVector) < 0)
                        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                            "Each psi knot vector must be nondecreasing.");
                    end

                    psiKnotPoints{iDim} = knotVector;
                end
            end

            self = GriddedStreamfunction();
            fitTrajectorySplines( ...
                self, ...
                trajectories, ...
                psiKnotPoints, ...
                psiS, ...
                fastKnotPoints, ...
                fastS, ...
                mesoscaleConstraint, ...
                options.buildDecomposition, ...
                sampleData);
        end
    end

    methods
        function decomposition = get.decomposition(self)
            if isempty(self.fixedFrameBackgroundTrajectories)
                decomposition = struct();
                return
            end

            decomposition = struct( ...
                "fixedFrame", struct( ...
                    "background", self.fixedFrameBackgroundTrajectories, ...
                    "mesoscale", self.fixedFrameMesoscaleTrajectories, ...
                    "submesoscale", self.fixedFrameSubmesoscaleTrajectories), ...
                "centeredFrame", struct( ...
                    "mesoscale", self.centeredFrameMesoscaleTrajectories, ...
                    "submesoscale", self.centeredFrameSubmesoscaleTrajectories));
        end

        function fastKnotPoints = get.fastKnotPoints(self)
            if isempty(self.centerOfMassTrajectory)
                fastKnotPoints = zeros(0, 1);
                return
            end

            fastKnotPoints = reshape(self.centerOfMassTrajectory.x.knotPoints, [], 1);
        end

        function fastS = get.fastS(self)
            if isempty(self.centerOfMassTrajectory)
                fastS = [];
                return
            end

            fastS = reshape(self.centerOfMassTrajectory.x.S, 1, []);
        end

        function psiKnotPoints = get.psiKnotPoints(self)
            if isempty(self.streamfunctionSpline)
                psiKnotPoints = cell(1, 0);
                return
            end

            psiKnotPoints = self.streamfunctionSpline.knotPoints;
        end

        function psiS = get.psiS(self)
            if isempty(self.streamfunctionSpline)
                psiS = [];
                return
            end

            psiS = reshape(self.streamfunctionSpline.S, 1, []);
        end
    end

    methods (Access = protected)
        function restoreOptionalPersistedPropertiesFromGroup(self, ~)
            self.refreshObservedTrajectorySampleData();
        end
    end

    methods
        function values = uBackground(self, t)
            % Evaluate the fitted background x-velocity.
            %
            % `uBackground` is a derived evaluation of the solved
            % `backgroundTrajectory`, not an additional fitted state
            % variable.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % plot(tFit, fit.uBackground(tFit))
            % xlabel("t (s)")
            % ylabel("u^{bg} (m/s)")
            % ```
            %
            % - Topic: Evaluate derived fields — Background evaluation
            % - nav_order: 1
            % - Declaration: values = uBackground(self,t)
            % - Parameter t: evaluation times in seconds
            % - Returns values: background x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.backgroundTrajectory.u(t);
        end

        function values = vBackground(self, t)
            % Evaluate the fitted background y-velocity.
            %
            % `vBackground` is a derived evaluation of the solved
            % `backgroundTrajectory`, not an additional fitted state
            % variable.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % plot(tFit, fit.vBackground(tFit))
            % xlabel("t (s)")
            % ylabel("v^{bg} (m/s)")
            % ```
            %
            % - Topic: Evaluate derived fields — Background evaluation
            % - nav_order: 2
            % - Declaration: values = vBackground(self,t)
            % - Parameter t: evaluation times in seconds
            % - Returns values: background y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.backgroundTrajectory.v(t);
        end

        function psiValue = psiMesoscale(self, t, x, y)
            % Evaluate the fitted mesoscale streamfunction.
            %
            % `psiMesoscale` evaluates the solved centered-frame spline at
            % fixed-frame query points by subtracting
            % `centerOfMassTrajectory` first. It is a derived evaluation
            % of `streamfunctionSpline`, not additional solved state.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % xCom = fit.centerOfMassTrajectory.x(tFit);
            % yCom = fit.centerOfMassTrajectory.y(tFit);
            % psiCenter = fit.psiMesoscale(tFit, xCom, yCom);
            % plot(tFit, psiCenter)
            % xlabel("t (s)")
            % ylabel("\psi")
            % ```
            %
            % - Topic: Evaluate derived fields — Mesoscale field evaluation
            % - nav_order: 1
            % - Declaration: psiValue = psiMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: mesoscale streamfunction values with the same shape as `x`
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            psiValue = self.streamfunctionSpline.valueAtPoints(q, r, tEval);
        end

        function uValue = uMesoscale(self, t, x, y)
            % Evaluate the mesoscale x-velocity $$-\psi_{\tilde{y}}$$.
            %
            % This is a derived velocity evaluation of the solved
            % mesoscale spline in centered coordinates.
            %
            % ```matlab
            % trajectory = fit.observedTrajectories(1);
            % ti = trajectory.t;
            % uMeso = fit.uMesoscale(ti, trajectory.x(ti), trajectory.y(ti));
            % plot(ti, uMeso)
            % xlabel("t (s)")
            % ylabel("u^{meso} (m/s)")
            % ```
            %
            % - Topic: Evaluate derived fields — Mesoscale field evaluation
            % - nav_order: 2
            % - Declaration: uValue = uMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: mesoscale x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            uValue = -self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 1 0]);
        end

        function vValue = vMesoscale(self, t, x, y)
            % Evaluate the mesoscale y-velocity $$\psi_{\tilde{x}}$$.
            %
            % This is a derived velocity evaluation of the solved
            % mesoscale spline in centered coordinates.
            %
            % ```matlab
            % trajectory = fit.observedTrajectories(1);
            % ti = trajectory.t;
            % vMeso = fit.vMesoscale(ti, trajectory.x(ti), trajectory.y(ti));
            % plot(ti, vMeso)
            % xlabel("t (s)")
            % ylabel("v^{meso} (m/s)")
            % ```
            %
            % - Topic: Evaluate derived fields — Mesoscale field evaluation
            % - nav_order: 3
            % - Declaration: vValue = vMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: mesoscale y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            vValue = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[1 0 0]);
        end

        function values = sigma_n(self, t, x, y)
            % Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.
            %
            % This is a derived diagnostic computed from the solved
            % mesoscale spline.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % xCom = fit.centerOfMassTrajectory.x(tFit);
            % yCom = fit.centerOfMassTrajectory.y(tFit);
            % plot(tFit, fit.sigma_n(tFit, xCom, yCom))
            % xlabel("t (s)")
            % ylabel("\sigma_n (s^{-1})")
            % ```
            %
            % - Topic: Evaluate derived fields — Derived diagnostics
            % - nav_order: 1
            % - Declaration: values = sigma_n(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: normal strain in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            values = -2 * self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[1 1 0]);
        end

        function values = sigma_s(self, t, x, y)
            % Evaluate the shear strain field $$\sigma_s = \psi_{\tilde{x}\tilde{x}} - \psi_{\tilde{y}\tilde{y}}$$.
            %
            % This is a derived diagnostic computed from the solved
            % mesoscale spline.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % xCom = fit.centerOfMassTrajectory.x(tFit);
            % yCom = fit.centerOfMassTrajectory.y(tFit);
            % plot(tFit, fit.sigma_s(tFit, xCom, yCom))
            % xlabel("t (s)")
            % ylabel("\sigma_s (s^{-1})")
            % ```
            %
            % - Topic: Evaluate derived fields — Derived diagnostics
            % - nav_order: 2
            % - Declaration: values = sigma_s(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: shear strain in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            dqq = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[2 0 0]);
            drr = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 2 0]);
            values = dqq - drr;
        end

        function values = zeta(self, t, x, y)
            % Evaluate the relative-vorticity field $$\zeta = \psi_{\tilde{x}\tilde{x}} + \psi_{\tilde{y}\tilde{y}}$$.
            %
            % This is a derived diagnostic computed from the solved
            % mesoscale spline.
            %
            % ```matlab
            % tFit = fit.fitSupportTimes;
            % xCom = fit.centerOfMassTrajectory.x(tFit);
            % yCom = fit.centerOfMassTrajectory.y(tFit);
            % plot(tFit, fit.zeta(tFit, xCom, yCom))
            % xlabel("t (s)")
            % ylabel("\zeta (s^{-1})")
            % ```
            %
            % - Topic: Evaluate derived fields — Derived diagnostics
            % - nav_order: 3
            % - Declaration: values = zeta(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: relative vorticity in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            dqq = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[2 0 0]);
            drr = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 2 0]);
            values = dqq + drr;
        end
    end

    methods (Access = private)
        fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS, mesoscaleConstraint, buildDecomposition, sampleData)
        [backgroundTrajectory, decomposition] = decomposeTrajectorySet(self, trajectories)
    end

    methods (Access = {?GriddedStreamfunction, ?GriddedStreamfunctionBootstrap, ?GriddedStreamfunctionUnitTests})
        sampleData = decompositionSampleData(self, trajectories)
    end

    methods (Static)
        theta = visualPrincipalStrainAngle(sigma_n, sigma_s, options)
    end

    methods (Static, Access = {?GriddedStreamfunction, ?GriddedStreamfunctionBootstrap, ?GriddedStreamfunctionUnitTests, ?GriddedStreamfunctionBootstrapUnitTests})
        sampleData = sampleTrajectoryData(trajectories)
        sampleData = resampledTrajectorySampleData(observedSampleData, sampledIndices)
    end

    methods (Static, Hidden)
        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, { ...
                'observedTrajectories', ...
                'fixedFrameBackgroundTrajectories', ...
                'fixedFrameMesoscaleTrajectories', ...
                'fixedFrameSubmesoscaleTrajectories', ...
                'centeredFrameMesoscaleTrajectories', ...
                'centeredFrameSubmesoscaleTrajectories', ...
                'mesoscaleConstraint', ...
                'representativeTimes', ...
                'fitSupportTimes', ...
                'mesoscaleDegreesOfFreedom'});
            vars.streamfunctionSpline = TensorSpline.annotatedClassFromGroup(group.groupWithName('streamfunctionSpline'));
            vars.centerOfMassTrajectory = TrajectorySpline.annotatedClassFromGroup(group.groupWithName('centerOfMassTrajectory'));
            vars.backgroundTrajectory = TrajectorySpline.annotatedClassFromGroup(group.groupWithName('backgroundTrajectory'));
            GriddedStreamfunction.validateCanonicalState(vars);

            self = GriddedStreamfunction();
            self.streamfunctionSpline = vars.streamfunctionSpline;
            self.observedTrajectories = reshape(vars.observedTrajectories, [], 1);
            self.centerOfMassTrajectory = vars.centerOfMassTrajectory;
            self.backgroundTrajectory = reshape(vars.backgroundTrajectory, [], 1);
            self.fixedFrameBackgroundTrajectories = reshape(vars.fixedFrameBackgroundTrajectories, [], 1);
            self.fixedFrameMesoscaleTrajectories = reshape(vars.fixedFrameMesoscaleTrajectories, [], 1);
            self.fixedFrameSubmesoscaleTrajectories = reshape(vars.fixedFrameSubmesoscaleTrajectories, [], 1);
            self.centeredFrameMesoscaleTrajectories = reshape(vars.centeredFrameMesoscaleTrajectories, [], 1);
            self.centeredFrameSubmesoscaleTrajectories = reshape(vars.centeredFrameSubmesoscaleTrajectories, [], 1);
            self.mesoscaleConstraint = string(vars.mesoscaleConstraint);
            self.representativeTimes = vars.representativeTimes;
            self.fitSupportTimes = vars.fitSupportTimes;
            self.mesoscaleDegreesOfFreedom = vars.mesoscaleDegreesOfFreedom;
            self.refreshObservedTrajectorySampleData();
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0, 0);
            propertyAnnotations(end+1) = CAObjectProperty('streamfunctionSpline', 'Fitted centered-frame mesoscale streamfunction spline.', className='TensorSpline', sizeText='(1,1)');
            propertyAnnotations(end+1) = CAObjectProperty('observedTrajectories', 'Observed drifter trajectory splines used for the fit.', className='TrajectorySpline', sizeText='nonempty vector');
            propertyAnnotations(end+1) = CAObjectProperty('centerOfMassTrajectory', 'Fitted center-of-mass trajectory.', className='TrajectorySpline', sizeText='(1,1)');
            propertyAnnotations(end+1) = CAObjectProperty('backgroundTrajectory', 'Fitted common background trajectory.', className='TrajectorySpline', sizeText='vector or empty array');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameBackgroundTrajectories', 'Stored fixed-frame background decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameMesoscaleTrajectories', 'Stored fixed-frame mesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameSubmesoscaleTrajectories', 'Stored fixed-frame submesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('centeredFrameMesoscaleTrajectories', 'Stored centered-frame mesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('centeredFrameSubmesoscaleTrajectories', 'Stored centered-frame submesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAPropertyAnnotation('mesoscaleConstraint', 'Hard mesoscale constraint applied to the fit.');
            propertyAnnotations(end+1) = CADimensionProperty('representativeTimes', '', 'Representative pooled times from the stride-rule fast basis.');
            propertyAnnotations(end+1) = CADimensionProperty('fitSupportTimes', '', 'Sorted unique observation times used as trajectory support.');
            propertyAnnotations(end+1) = CANumericProperty('mesoscaleDegreesOfFreedom', {}, '', 'Stored identifiable mesoscale degrees of freedom after gauge reduction and hard constraints.');
        end

        function names = classRequiredPropertyNames()
            names = { ...
                'streamfunctionSpline', ...
                'observedTrajectories', ...
                'centerOfMassTrajectory', ...
                'backgroundTrajectory', ...
                'fixedFrameBackgroundTrajectories', ...
                'fixedFrameMesoscaleTrajectories', ...
                'fixedFrameSubmesoscaleTrajectories', ...
                'centeredFrameMesoscaleTrajectories', ...
                'centeredFrameSubmesoscaleTrajectories', ...
                'mesoscaleConstraint', ...
                'representativeTimes', ...
                'fitSupportTimes', ...
                'mesoscaleDegreesOfFreedom'};
        end
    end

    methods (Static, Access = private)
        function validateCanonicalState(options)
            if isempty(options.streamfunctionSpline) && isempty(options.observedTrajectories) && ...
                    isempty(options.centerOfMassTrajectory) && isempty(options.backgroundTrajectory) && ...
                    isempty(options.fixedFrameBackgroundTrajectories) && isempty(options.fixedFrameMesoscaleTrajectories) && ...
                    isempty(options.fixedFrameSubmesoscaleTrajectories) && isempty(options.centeredFrameMesoscaleTrajectories) && ...
                    isempty(options.centeredFrameSubmesoscaleTrajectories) && isempty(options.representativeTimes) && ...
                    isempty(options.fitSupportTimes) && string(options.mesoscaleConstraint) == "none" && ...
                    options.mesoscaleDegreesOfFreedom == 0
                return
            end

            if ~(isa(options.streamfunctionSpline, 'TensorSpline') && isscalar(options.streamfunctionSpline))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'streamfunctionSpline must be a scalar TensorSpline.');
            end
            if ~(isa(options.observedTrajectories, 'TrajectorySpline') && isvector(options.observedTrajectories) && ~isempty(options.observedTrajectories))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'observedTrajectories must be a nonempty vector of TrajectorySpline objects.');
            end
            if ~(isa(options.centerOfMassTrajectory, 'TrajectorySpline') && isscalar(options.centerOfMassTrajectory))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'centerOfMassTrajectory must be a scalar TrajectorySpline.');
            end
            if ~(isa(options.backgroundTrajectory, 'TrajectorySpline') && isvector(options.backgroundTrajectory))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'backgroundTrajectory must be a TrajectorySpline or empty TrajectorySpline array.');
            end

            trajectoryArrays = { ...
                options.fixedFrameBackgroundTrajectories, ...
                options.fixedFrameMesoscaleTrajectories, ...
                options.fixedFrameSubmesoscaleTrajectories, ...
                options.centeredFrameMesoscaleTrajectories, ...
                options.centeredFrameSubmesoscaleTrajectories};
            for iArray = 1:numel(trajectoryArrays)
                if ~(isa(trajectoryArrays{iArray}, 'TrajectorySpline') && isvector(trajectoryArrays{iArray}))
                    error('GriddedStreamfunction:InvalidCanonicalState', ...
                        'Stored decomposition trajectories must be TrajectorySpline vectors.');
                end
            end

            if ~isempty(options.backgroundTrajectory) && ~isscalar(options.backgroundTrajectory)
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'backgroundTrajectory must be a scalar TrajectorySpline when it is present.');
            end
            if ~isscalar(options.mesoscaleDegreesOfFreedom) || ~isfinite(options.mesoscaleDegreesOfFreedom) || ...
                    options.mesoscaleDegreesOfFreedom < 0 || round(options.mesoscaleDegreesOfFreedom) ~= options.mesoscaleDegreesOfFreedom
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'mesoscaleDegreesOfFreedom must be a finite nonnegative integer scalar.');
            end
        end

        representativeTimes = representativeObservationTimes(tCell)
        psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
        Aeq = mesoscaleConstraintMatrix(psiKnotPoints, psiS, G, mesoscaleConstraint)
        fastState = fastBasisState(allT, fastKnotPoints, fastS, shouldComputeDerivative)
    end

    methods (Access = private)
        function refreshObservedTrajectorySampleData(self)
            if isempty(self.observedTrajectories)
                self.observedTrajectorySampleData = struct([]);
            else
                self.observedTrajectorySampleData = GriddedStreamfunction.sampleTrajectoryData(self.observedTrajectories);
            end
        end
    end
end
