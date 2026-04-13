classdef GriddedStreamfunctionBootstrap < CAAnnotatedClass
    % Bootstrap whole-drifter gridded-streamfunction fits and consensus scores.
    %
    % `GriddedStreamfunctionBootstrap` fits one deterministic
    % `GriddedStreamfunction` on the full drifter ensemble, then constructs
    % a whole-drifter bootstrap ensemble by resampling the input
    % `TrajectorySpline` objects with replacement. Each replicate is
    % summarized by COM-local mesoscale diagnostics
    % $$u_c(t), v_c(t), \sigma_n(t), \sigma_s(t), \zeta(t)$$ evaluated at
    % the fitted center-of-mass trajectory. The class also reports the
    % scalar diffusivity diagnostic
    % $$\kappa = \langle (x_{\mathrm{sm}}^2 + y_{\mathrm{sm}}^2)/(4\Delta t)\rangle,$$
    % together with a scalar mean coherence over the lower-frequency half
    % of the mean coherence spectrum and the full mean coherence spectrum
    % between the fitted centered-frame mesoscale and submesoscale
    % velocities when the jLab coherence tools are available.
    %
    % The bootstrap ensemble is ranked by a consensus score formed from
    % time-local `KernelDensityEstimate` fits to the bootstrap cloud. The
    % score uses a 2-D density fit for `(uCenter, vCenter)`, a 2-D
    % density fit for `(sigma_n, sigma_s)`, and a 1-D density fit for
    % `zeta`, omitting constrained blocks when `mesoscaleConstraint`
    % forces them to vanish.
    %
    % The structural mesoscale degrees of freedom are determined only by
    % the resolved mesoscale basis and the selected hard constraint, so
    % the full fit and every reconstructed bootstrap replicate share the
    % same scalar `mesoscaleDegreesOfFreedom`.
    %
    % ```matlab
    % bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
    %     trajectories, nBootstraps=100, randomSeed=7);
    % fullFit = bootstrap.fullFit;
    % bestFit = bootstrap.bestFit();
    % quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
    % plot(bootstrap.queryTimes, quantiles.zeta(:, 2))
    % ```
    %
    % - Topic: Create a bootstrap ensemble
    % - Topic: Persist and restart ensembles
    % - Topic: Inspect ensemble setup
    % - Topic: Inspect primary outputs
    % - Topic: Inspect primary outputs — Full-data fit
    % - Topic: Inspect primary outputs — Best bootstrap replicate
    % - Topic: Reconstruct bootstrap fits
    % - Topic: Evaluate derived diagnostics
    % - Topic: Evaluate derived diagnostics — Full-fit diagnostics
    % - Topic: Evaluate derived diagnostics — Best-fit diagnostics
    % - Topic: Evaluate derived diagnostics — Ensemble summaries
    % - Declaration: classdef GriddedStreamfunctionBootstrap < CAAnnotatedClass

    properties (SetAccess = private)
        % Full-data gridded-streamfunction fit used as the reference solution.
        %
        % `fullFit` is the deterministic `GriddedStreamfunction` fit on
        % the original input drifters before any resampling is applied.
        % This is one of the ensemble's primary outputs.
        %
        % ```matlab
        % fullFit = bootstrap.fullFit;
        % tFit = fullFit.fitSupportTimes;
        % plot(fullFit.centerOfMassTrajectory.x(tFit), fullFit.centerOfMassTrajectory.y(tFit))
        % axis equal
        % xlabel("x (m)")
        % ylabel("y (m)")
        % ```
        %
        % - Topic: Inspect primary outputs — Full-data fit
        % - nav_order: 1
        fullFit

        % Original drifter trajectories used to seed the bootstrap ensemble.
        %
        % `observedTrajectories` preserves the original
        % `TrajectorySpline` column vector in the user-supplied drifter
        % order.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 1
        observedTrajectories

        % Whole-drifter resampling indices for each bootstrap replicate.
        %
        % `bootstrapIndices(i,:)` contains the drifter indices drawn with
        % replacement for bootstrap replicate `i`.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 2
        bootstrapIndices

        % Times used to store COM-local bootstrap summaries.
        %
        % `queryTimes` is a column vector of times in seconds. When not
        % supplied explicitly, it is the sorted unique observation times
        % restricted to the common overlap interval of the original
        % drifters.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 5
        queryTimes

        % Times used to compute the consensus score.
        %
        % `scoreTimes` is always an exact subset of `queryTimes`.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 6
        scoreTimes

        % Number of bootstrap replicates in the ensemble.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 3
        nBootstraps (1,1) double {mustBeInteger,mustBeNonnegative} = 0

        % Random-number seed used for whole-drifter resampling.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 4
        randomSeed (1,1) double {mustBeInteger,mustBeFinite} = 0
    end

    properties (Dependent)
        % Full-data COM-local mesoscale summary evaluated on `queryTimes`.
        %
        % `fullSummary` contains the fields `uCenter`, `vCenter`,
        % `sigma_n`, `sigma_s`, and `zeta`, each stored as a column vector
        % aligned with `queryTimes`. These are the primary COM-local
        % diagnostic outputs of the deterministic `fullFit`, evaluated at
        % the fitted center-of-mass positions
        % `(fullFit.centerOfMassTrajectory.x(queryTimes),`
        % `fullFit.centerOfMassTrajectory.y(queryTimes))`.
        %
        % ```matlab
        % fullSummary = bootstrap.fullSummary;
        % plot(bootstrap.queryTimes, fullSummary.zeta)
        % xlabel("t (s)")
        % ylabel("\zeta_c (s^{-1})")
        % ```
        %
        % - Topic: Inspect primary outputs — Full-data fit
        % - nav_order: 2
        fullSummary

        % Bootstrap COM-local mesoscale summaries evaluated on `queryTimes`.
        %
        % `summary.uCenter`, `summary.vCenter`, `summary.sigma_n`,
        % `summary.sigma_s`, and `summary.zeta` are arrays of size
        % `[numel(queryTimes) nBootstraps]`. Row `i` corresponds to
        % `queryTimes(i)`, and column `j` corresponds to bootstrap
        % replicate `j`. These are derived ensemble summaries for
        % uncertainty analysis; use `fitForBootstrap` or `bestFit` when
        % you need a full reconstructed fit rather than the stored summary
        % arrays.
        %
        % ```matlab
        % summary = bootstrap.summary;
        % medianUCenter = median(summary.uCenter, 2);
        % plot(bootstrap.queryTimes, medianUCenter)
        % ```
        %
        % - Topic: Evaluate derived diagnostics — Ensemble summaries
        % - nav_order: 1
        summary

        % Full-fit scalar submesoscale diffusivity diagnostic.
        %
        % `fullFitKappa` is computed from the fitted full-data
        % submesoscale trajectories and is eager because it only requires
        % one fit-level diagnostic evaluation.
        %
        % - Topic: Evaluate derived diagnostics — Full-fit diagnostics
        % - nav_order: 1
        fullFitKappa

        % Full-fit scalar mean coherence between mesoscale and submesoscale velocities.
        %
        % `fullFitCoherence` is the mean of the finite lower-frequency
        % half of `fullFitCoherenceSpectrum.coherence`.
        %
        % - Topic: Evaluate derived diagnostics — Full-fit diagnostics
        % - nav_order: 2
        fullFitCoherence

        % Full-fit mean coherence spectrum on the common overlap grid.
        %
        % `fullFitCoherenceSpectrum.frequency` and
        % `fullFitCoherenceSpectrum.coherence` are column vectors. When
        % the coherence backend is unavailable or no usable common uniform
        % overlap grid exists, both vectors are empty.
        %
        % - Topic: Evaluate derived diagnostics — Full-fit diagnostics
        % - nav_order: 3
        fullFitCoherenceSpectrum

        % Structural mesoscale degrees of freedom shared by the ensemble.
        %
        % `mesoscaleDegreesOfFreedom` is a single scalar determined by the
        % resolved mesoscale basis after gauge reduction together with the
        % selected hard mesoscale constraint. It is the same value
        % reported by `fullFit.mesoscaleDegreesOfFreedom`,
        % `bestFit().mesoscaleDegreesOfFreedom`, and every
        % `fitForBootstrap(iBootstrap).mesoscaleDegreesOfFreedom`.
        %
        % - Topic: Inspect ensemble setup
        % - nav_order: 7
        mesoscaleDegreesOfFreedom

        % Best-bootstrap scalar submesoscale diffusivity diagnostic.
        %
        % `bestFitKappa` is eager and corresponds to the fit returned by
        % `bestFit()`.
        %
        % - Topic: Evaluate derived diagnostics — Best-fit diagnostics
        % - nav_order: 1
        bestFitKappa

        % Best-bootstrap scalar mean coherence.
        %
        % `bestFitCoherence` is the mean of the finite lower-frequency
        % half of `bestFitCoherenceSpectrum.coherence`.
        %
        % - Topic: Evaluate derived diagnostics — Best-fit diagnostics
        % - nav_order: 2
        bestFitCoherence

        % Best-bootstrap mean coherence spectrum on the common overlap grid.
        %
        % - Topic: Evaluate derived diagnostics — Best-fit diagnostics
        % - nav_order: 3
        bestFitCoherenceSpectrum

        % Lazy scalar diffusivity diagnostics for each bootstrap replicate.
        %
        % `bootstrapKappa` is a row vector of length `nBootstraps`.
        %
        % ```matlab
        % kappa = bootstrap.bootstrapKappa;
        % plot(kappa, ".")
        % xlabel("bootstrap replicate")
        % ylabel("\kappa")
        % ```
        %
        % - Topic: Evaluate derived diagnostics — Ensemble summaries
        % - nav_order: 3
        bootstrapKappa

        % Lazy scalar coherence diagnostics for each bootstrap replicate.
        %
        % `bootstrapCoherence` is a row vector of length `nBootstraps`.
        % Each value uses the finite lower-frequency half of the stored
        % mean coherence spectrum. When coherence is unavailable, the
        % vector is filled with `NaN`.
        %
        % ```matlab
        % coherence = bootstrap.bootstrapCoherence;
        % plot(coherence, ".")
        % xlabel("bootstrap replicate")
        % ylabel("coherence")
        % ```
        %
        % - Topic: Evaluate derived diagnostics — Ensemble summaries
        % - nav_order: 4
        bootstrapCoherence

        % Consensus-score components and joint score for each bootstrap fit.
        %
        % `scores.uv`, `scores.strain`, `scores.zeta`, and `scores.joint`
        % are row vectors of length `nBootstraps`. These are derived
        % KDE-based ranking diagnostics for the bootstrap ensemble rather
        % than additional fitted state.
        %
        % ```matlab
        % scores = bootstrap.scores;
        % [~, iBest] = max(scores.joint);
        % plot(scores.joint, ".")
        % hold on
        % plot(iBest, scores.joint(iBest), "o")
        % ```
        %
        % - Topic: Evaluate derived diagnostics — Ensemble summaries
        % - nav_order: 5
        scores

        % Exact reconstruction metadata for each bootstrap replicate.
        %
        % `bootstrapMetadata.fastKnotPoints{i}` stores the resolved fast
        % temporal knot vector used by bootstrap replicate `i`, and
        % `bootstrapMetadata.psiKnotPoints{i}` stores the corresponding
        % resolved mesoscale knot cell `{qKnot, rKnot, tKnot}`.
        %
        % - Topic: Reconstruct bootstrap fits
        % - nav_order: 1
        bootstrapMetadata
    end

    properties (Hidden, SetAccess = private)
        fullSummaryUCenter = zeros(0, 1)
        fullSummaryVCenter = zeros(0, 1)
        fullSummarySigmaN = zeros(0, 1)
        fullSummarySigmaS = zeros(0, 1)
        fullSummaryZeta = zeros(0, 1)
        summaryUCenter = zeros(0, 0)
        summaryVCenter = zeros(0, 0)
        summarySigmaN = zeros(0, 0)
        summarySigmaS = zeros(0, 0)
        summaryZeta = zeros(0, 0)
        fullFitKappaValue (1,1) double {mustBeReal} = NaN
        fullFitCoherenceValue (1,1) double {mustBeReal} = NaN
        fullFitCoherenceFrequency (:,1) double {mustBeReal} = zeros(0, 1)
        fullFitCoherenceValues (:,1) double {mustBeReal} = zeros(0, 1)
        bestFitKappaValue (1,1) double {mustBeReal} = NaN
        bestFitCoherenceValue (1,1) double {mustBeReal} = NaN
        bestFitCoherenceFrequency (:,1) double {mustBeReal} = zeros(0, 1)
        bestFitCoherenceValues (:,1) double {mustBeReal} = zeros(0, 1)
        scoreUv = zeros(1, 0)
        scoreStrain = zeros(1, 0)
        scoreZeta = zeros(1, 0)
        scoreJoint = zeros(1, 0)
        bootstrapFitMetadata = GriddedStreamfunctionBootstrapFitMetadata.empty(0, 1)
        bootstrapKappaCache = zeros(1, 0)
        bootstrapCoherenceCache = zeros(1, 0)
    end

    properties (Access = private)
        observedTrajectorySampleData struct = struct([])
    end

    properties (Dependent, Hidden, SetAccess = private)
        bootstrapIndex
        trajectoryIndex
        fullFitCoherenceFrequencyIndex
        bestFitCoherenceFrequencyIndex
    end

    methods
        function self = GriddedStreamfunctionBootstrap(options)
            % Create a bootstrap from canonical restart-state properties.
            %
            % Use this low-level constructor when the full fit, stored
            % summaries, scores, and reconstruction metadata already exist,
            % for example after reading a persisted restart file. For
            % ordinary bootstrap fitting from observed trajectories, use
            % `GriddedStreamfunctionBootstrap.fromTrajectories(...)`.
            %
            % - Topic: Create a bootstrap ensemble
            % - Declaration: self = GriddedStreamfunctionBootstrap(options)
            % - Parameter options.fullFit: full-data `GriddedStreamfunction` fit
            % - Parameter options.observedTrajectories: original drifter trajectories used for resampling
            % - Parameter options.bootstrapIndices: whole-drifter resampling indices for each replicate
            % - Parameter options.queryTimes: stored summary times
            % - Parameter options.scoreTimes: stored consensus-score times
            % - Parameter options.fullSummaryUCenter: persisted full-fit `uCenter` summary values
            % - Parameter options.fullSummaryVCenter: persisted full-fit `vCenter` summary values
            % - Parameter options.fullSummarySigmaN: persisted full-fit `sigma_n` summary values
            % - Parameter options.fullSummarySigmaS: persisted full-fit `sigma_s` summary values
            % - Parameter options.fullSummaryZeta: persisted full-fit `zeta` summary values
            % - Parameter options.summaryUCenter: persisted bootstrap `uCenter` summary matrix
            % - Parameter options.summaryVCenter: persisted bootstrap `vCenter` summary matrix
            % - Parameter options.summarySigmaN: persisted bootstrap `sigma_n` summary matrix
            % - Parameter options.summarySigmaS: persisted bootstrap `sigma_s` summary matrix
            % - Parameter options.summaryZeta: persisted bootstrap `zeta` summary matrix
            % - Parameter options.fullFitKappaValue: persisted full-fit scalar diffusivity diagnostic
            % - Parameter options.fullFitCoherenceValue: persisted full-fit scalar mean coherence
            % - Parameter options.fullFitCoherenceFrequency: persisted full-fit coherence-spectrum frequencies
            % - Parameter options.fullFitCoherenceValues: persisted full-fit coherence-spectrum values
            % - Parameter options.bestFitKappaValue: persisted best-fit scalar diffusivity diagnostic
            % - Parameter options.bestFitCoherenceValue: persisted best-fit scalar mean coherence
            % - Parameter options.bestFitCoherenceFrequency: persisted best-fit coherence-spectrum frequencies
            % - Parameter options.bestFitCoherenceValues: persisted best-fit coherence-spectrum values
            % - Parameter options.scoreUv: persisted bootstrap velocity consensus scores
            % - Parameter options.scoreStrain: persisted bootstrap strain consensus scores
            % - Parameter options.scoreZeta: persisted bootstrap vorticity consensus scores
            % - Parameter options.scoreJoint: persisted bootstrap joint consensus scores
            % - Parameter options.bootstrapFitMetadata: persisted exact reconstruction metadata for each bootstrap replicate
            % - Parameter options.nBootstraps: number of bootstrap replicates
            % - Parameter options.randomSeed: random seed used for resampling
            % - Returns self: canonical `GriddedStreamfunctionBootstrap` instance
            arguments
                options.fullFit = GriddedStreamfunction.empty(0, 1)
                options.observedTrajectories = TrajectorySpline.empty(0, 1)
                options.bootstrapIndices double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.queryTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.scoreTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fullSummaryUCenter (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fullSummaryVCenter (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fullSummarySigmaN (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fullSummarySigmaS (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fullSummaryZeta (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.summaryUCenter double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.summaryVCenter double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.summarySigmaN double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.summarySigmaS double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.summaryZeta double {mustBeReal,mustBeFinite} = zeros(0, 0)
                options.fullFitKappaValue (1,1) double {mustBeReal} = NaN
                options.fullFitCoherenceValue (1,1) double {mustBeReal} = NaN
                options.fullFitCoherenceFrequency (:,1) double {mustBeReal} = zeros(0, 1)
                options.fullFitCoherenceValues (:,1) double {mustBeReal} = zeros(0, 1)
                options.bestFitKappaValue (1,1) double {mustBeReal} = NaN
                options.bestFitCoherenceValue (1,1) double {mustBeReal} = NaN
                options.bestFitCoherenceFrequency (:,1) double {mustBeReal} = zeros(0, 1)
                options.bestFitCoherenceValues (:,1) double {mustBeReal} = zeros(0, 1)
                options.scoreUv double {mustBeReal,mustBeFinite} = zeros(1, 0)
                options.scoreStrain double {mustBeReal,mustBeFinite} = zeros(1, 0)
                options.scoreZeta double {mustBeReal,mustBeFinite} = zeros(1, 0)
                options.scoreJoint double {mustBeReal,mustBeFinite} = zeros(1, 0)
                options.bootstrapFitMetadata = GriddedStreamfunctionBootstrapFitMetadata.empty(0, 1)
                options.nBootstraps (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.randomSeed (1,1) double {mustBeInteger,mustBeFinite} = 0
            end

            self@CAAnnotatedClass();

            if nargin == 0
                return
            end

            GriddedStreamfunctionBootstrap.validateCanonicalState(options);
            self.fullFit = options.fullFit;
            self.observedTrajectories = reshape(options.observedTrajectories, [], 1);
            self.bootstrapIndices = options.bootstrapIndices;
            self.queryTimes = options.queryTimes;
            self.scoreTimes = options.scoreTimes;
            self.fullSummaryUCenter = options.fullSummaryUCenter;
            self.fullSummaryVCenter = options.fullSummaryVCenter;
            self.fullSummarySigmaN = options.fullSummarySigmaN;
            self.fullSummarySigmaS = options.fullSummarySigmaS;
            self.fullSummaryZeta = options.fullSummaryZeta;
            self.summaryUCenter = options.summaryUCenter;
            self.summaryVCenter = options.summaryVCenter;
            self.summarySigmaN = options.summarySigmaN;
            self.summarySigmaS = options.summarySigmaS;
            self.summaryZeta = options.summaryZeta;
            self.fullFitKappaValue = options.fullFitKappaValue;
            self.fullFitCoherenceValue = options.fullFitCoherenceValue;
            self.fullFitCoherenceFrequency = options.fullFitCoherenceFrequency;
            self.fullFitCoherenceValues = options.fullFitCoherenceValues;
            self.bestFitKappaValue = options.bestFitKappaValue;
            self.bestFitCoherenceValue = options.bestFitCoherenceValue;
            self.bestFitCoherenceFrequency = options.bestFitCoherenceFrequency;
            self.bestFitCoherenceValues = options.bestFitCoherenceValues;
            self.scoreUv = reshape(options.scoreUv, 1, []);
            self.scoreStrain = reshape(options.scoreStrain, 1, []);
            self.scoreZeta = reshape(options.scoreZeta, 1, []);
            self.scoreJoint = reshape(options.scoreJoint, 1, []);
            self.bootstrapFitMetadata = reshape(options.bootstrapFitMetadata, [], 1);
            self.nBootstraps = options.nBootstraps;
            self.randomSeed = options.randomSeed;
            self.refreshObservedTrajectorySampleData();
        end

        function fullSummary = get.fullSummary(self)
            fullSummary = struct( ...
                "uCenter", self.fullSummaryUCenter, ...
                "vCenter", self.fullSummaryVCenter, ...
                "sigma_n", self.fullSummarySigmaN, ...
                "sigma_s", self.fullSummarySigmaS, ...
                "zeta", self.fullSummaryZeta);
        end

        function summary = get.summary(self)
            summary = struct( ...
                "uCenter", self.summaryUCenter, ...
                "vCenter", self.summaryVCenter, ...
                "sigma_n", self.summarySigmaN, ...
                "sigma_s", self.summarySigmaS, ...
                "zeta", self.summaryZeta);
        end

        function kappa = get.fullFitKappa(self)
            kappa = self.fullFitKappaValue;
        end

        function coherence = get.fullFitCoherence(self)
            coherence = self.fullFitCoherenceValue;
        end

        function spectrum = get.fullFitCoherenceSpectrum(self)
            spectrum = GriddedStreamfunctionBootstrap.coherenceSpectrumFromStorage( ...
                self.fullFitCoherenceFrequency, self.fullFitCoherenceValues);
        end

        function kappa = get.bestFitKappa(self)
            kappa = self.bestFitKappaValue;
        end

        function coherence = get.bestFitCoherence(self)
            coherence = self.bestFitCoherenceValue;
        end

        function spectrum = get.bestFitCoherenceSpectrum(self)
            spectrum = GriddedStreamfunctionBootstrap.coherenceSpectrumFromStorage( ...
                self.bestFitCoherenceFrequency, self.bestFitCoherenceValues);
        end

        function kappa = get.bootstrapKappa(self)
            ensureBootstrapDiagnostics(self);
            kappa = self.bootstrapKappaCache;
        end

        function coherence = get.bootstrapCoherence(self)
            ensureBootstrapDiagnostics(self);
            coherence = self.bootstrapCoherenceCache;
        end

        function dof = get.mesoscaleDegreesOfFreedom(self)
            if isempty(self.fullFit)
                dof = 0;
                return
            end

            dof = self.fullFit.mesoscaleDegreesOfFreedom;
        end

        function scores = get.scores(self)
            scores = struct( ...
                "uv", self.scoreUv, ...
                "strain", self.scoreStrain, ...
                "zeta", self.scoreZeta, ...
                "joint", self.scoreJoint);
        end

        function metadata = get.bootstrapMetadata(self)
            nMetadata = numel(self.bootstrapFitMetadata);
            metadata = struct( ...
                "fastKnotPoints", {cell(nMetadata, 1)}, ...
                "psiKnotPoints", {cell(nMetadata, 1)});
            for iMetadata = 1:nMetadata
                metadata.fastKnotPoints{iMetadata} = self.bootstrapFitMetadata(iMetadata).fastKnotPoints;
                metadata.psiKnotPoints{iMetadata} = self.bootstrapFitMetadata(iMetadata).psiKnotPoints;
            end
        end

        function bootstrapIndex = get.bootstrapIndex(self)
            bootstrapIndex = reshape(1:self.nBootstraps, [], 1);
        end

        function trajectoryIndex = get.trajectoryIndex(self)
            trajectoryIndex = reshape(1:size(self.bootstrapIndices, 2), [], 1);
        end

        function frequencyIndex = get.fullFitCoherenceFrequencyIndex(self)
            frequencyIndex = reshape(1:numel(self.fullFitCoherenceFrequency), [], 1);
        end

        function frequencyIndex = get.bestFitCoherenceFrequencyIndex(self)
            frequencyIndex = reshape(1:numel(self.bestFitCoherenceFrequency), [], 1);
        end
    end

    methods
        function index = bestBootstrapIndex(self)
            % Return the bootstrap index with the highest joint consensus score.
            %
            % - Topic: Inspect primary outputs — Best bootstrap replicate
            % - nav_order: 1
            % - Declaration: index = bestBootstrapIndex(self)
            % - Returns index: 1-based bootstrap index of the top-ranked replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
            end
            arguments (Output)
                index (1,1) double
            end

            [~, index] = max(self.scoreJoint, [], 2);
        end

        function fit = fitForBootstrap(self, iBootstrap)
            % Reconstruct one bootstrap replicate exactly from saved metadata.
            %
            % Use this method to recover a full `GriddedStreamfunction`
            % object for a particular bootstrap replicate without storing
            % every replicate fit in memory.
            %
            % ```matlab
            % iBootstrap = bootstrap.bestBootstrapIndex();
            % fit = bootstrap.fitForBootstrap(iBootstrap);
            % tFit = fit.fitSupportTimes;
            % plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
            % axis equal
            % ```
            %
            % - Topic: Reconstruct bootstrap fits
            % - nav_order: 2
            % - Declaration: fit = fitForBootstrap(self,iBootstrap)
            % - Parameter iBootstrap: 1-based bootstrap replicate index
            % - Returns fit: reconstructed `GriddedStreamfunction` replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
                iBootstrap (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Output)
                fit (1,1) GriddedStreamfunction
            end

            if iBootstrap > self.nBootstraps
                error("GriddedStreamfunctionBootstrap:InvalidBootstrapIndex", ...
                    "iBootstrap must be between 1 and %d.", self.nBootstraps);
            end

            sampledTrajectories = reshape(self.observedTrajectories(self.bootstrapIndices(iBootstrap, :)), [], 1);
            metadata = self.bootstrapFitMetadata(iBootstrap);
            fitOptions = struct( ...
                "psiKnotPoints", {metadata.psiKnotPoints}, ...
                "psiS", self.fullFit.psiS, ...
                "fastKnotPoints", metadata.fastKnotPoints, ...
                "fastS", self.fullFit.fastS, ...
                "mesoscaleConstraint", self.fullFit.mesoscaleConstraint);
            fitOptions.sampleData = GriddedStreamfunction.resampledTrajectorySampleData( ...
                self.observedTrajectorySampleData, self.bootstrapIndices(iBootstrap, :));
            fitArguments = namedargs2cell(fitOptions);
            fit = GriddedStreamfunction.fromTrajectories(sampledTrajectories, fitArguments{:});
        end

        function fit = bestFit(self)
            % Reconstruct the top-ranked bootstrap replicate.
            %
            % The returned fit has the same structural
            % `mesoscaleDegreesOfFreedom` as `self.fullFit` and
            % `self.mesoscaleDegreesOfFreedom`; only the fitted spline
            % coefficients and diagnostics vary across replicates.
            %
            % ```matlab
            % fit = bootstrap.bestFit();
            % tFit = fit.fitSupportTimes;
            % plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
            % axis equal
            % ```
            %
            % - Topic: Inspect primary outputs — Best bootstrap replicate
            % - nav_order: 2
            % - Declaration: fit = bestFit(self)
            % - Returns fit: reconstructed `GriddedStreamfunction` for the best bootstrap replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
            end
            arguments (Output)
                fit (1,1) GriddedStreamfunction
            end

            fit = fitForBootstrap(self, bestBootstrapIndex(self));
        end

        function quantiles = summaryQuantiles(self, probabilities)
            % Return bootstrap quantiles for the stored summaries.
            %
            % The returned struct contains the time-varying summary fields
            % `uCenter`, `vCenter`, `sigma_n`, `sigma_s`, and `zeta` with
            % size `[numel(queryTimes) numel(probabilities)]`, together
            % with the scalar fields `kappa` and `coherence` of size
            % `[1 numel(probabilities)]`. The scalar coherence values use
            % the lower-frequency half of each bootstrap replicate's mean
            % coherence spectrum.
            %
            % ```matlab
            % quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
            % plot(bootstrap.queryTimes, quantiles.zeta(:, 2))
            % xlabel("t (s)")
            % ylabel("\zeta_c (s^{-1})")
            % ```
            %
            % - Topic: Evaluate derived diagnostics — Ensemble summaries
            % - nav_order: 2
            % - Declaration: quantiles = summaryQuantiles(self,probabilities)
            % - Parameter probabilities: vector of quantile probabilities on `[0,1]`
            % - Returns quantiles: struct of bootstrap quantiles for the stored summary fields
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
                probabilities {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
            end
            arguments (Output)
                quantiles struct
            end

            probabilities = reshape(probabilities, 1, []);
            if any(probabilities < 0 | probabilities > 1)
                error("GriddedStreamfunctionBootstrap:InvalidProbabilities", ...
                    "probabilities must lie on the closed interval [0, 1].");
            end

            quantiles = struct( ...
                "uCenter", quantile(self.summaryUCenter, probabilities, 2), ...
                "vCenter", quantile(self.summaryVCenter, probabilities, 2), ...
                "sigma_n", quantile(self.summarySigmaN, probabilities, 2), ...
                "sigma_s", quantile(self.summarySigmaS, probabilities, 2), ...
                "zeta", quantile(self.summaryZeta, probabilities, 2), ...
                "kappa", quantile(self.bootstrapKappa, probabilities, 2), ...
                "coherence", quantile(self.bootstrapCoherence, probabilities, 2));
        end
    end

    methods (Static)
        function self = fromFile(path)
            % Read a bootstrap ensemble from a NetCDF restart file.
            %
            % `fromFile` reconstructs the canonical bootstrap state written
            % by `writeToFile` without rerunning the whole-drifter
            % resampling workflow.
            %
            % ```matlab
            % bootstrap = GriddedStreamfunctionBootstrap.fromFile("gridded-bootstrap.nc");
            % bestFit = bootstrap.bestFit();
            % scores = bootstrap.scores;
            % plot(scores.joint, ".")
            % ```
            %
            % - Topic: Persist and restart ensembles
            % - nav_order: 1
            % - Declaration: self = fromFile(path)
            % - Parameter path: path to the NetCDF restart file
            % - Returns self: reconstructed `GriddedStreamfunctionBootstrap` ensemble
            arguments (Input)
                path {mustBeTextScalar}
            end

            self = griddedStreamfunctionBootstrapFromFile(char(path));
        end

        function self = fromTrajectories(trajectories, options)
            % Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.
            %
            % Use this factory when the uncertainty analysis should
            % reflect sensitivity to which drifters were observed, rather
            % than only the residual variability within one fitted
            % ensemble.
            %
            % The factory first fits the full-data
            % `GriddedStreamfunction`, then resamples the drifter
            % trajectories with replacement and fits one replicate per
            % bootstrap draw. Each replicate stores COM-local mesoscale
            % diagnostics on `queryTimes`, while scalar bootstrap
            % diagnostics are computed lazily on demand. The class also
            % stores the resolved spline knot vectors required to
            % reconstruct each replicate exactly later.
            %
            % ```matlab
            % bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
            %     trajectories, nBootstraps=100, randomSeed=7);
            % quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
            % plot(bootstrap.queryTimes, quantiles.uCenter(:, 2))
            % xlabel("t (s)")
            % ylabel("u_c (m/s)")
            % ```
            %
            % - Topic: Create a bootstrap ensemble
            % - nav_order: 1
            % - Declaration: self = fromTrajectories(trajectories,nBootstraps=...,randomSeed=...,queryTimes=...,scoreTimes=...,scoreStride=...,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter nBootstraps: number of whole-drifter bootstrap replicates, default `100`
            % - Parameter randomSeed: integer random seed used for resampling, default `0`
            % - Parameter queryTimes: optional strictly increasing query times for stored summaries
            % - Parameter scoreTimes: optional strictly increasing subset of `queryTimes` used for consensus scoring
            % - Parameter scoreStride: optional positive stride used to subsample `queryTimes` into `scoreTimes`, default `1`
            % - Parameter psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter fastS: optional fast temporal spline degree, default `3`
            % - Parameter mesoscaleConstraint: optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Returns self: fitted bootstrap ensemble
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.nBootstraps (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.nBootstraps, 3)} = 100
                options.randomSeed (1,1) double {mustBeInteger,mustBeFinite} = 0
                options.queryTimes = []
                options.scoreTimes = []
                options.scoreStride (1,1) double {mustBeInteger,mustBePositive} = 1
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger,mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger,mustBeNonnegative} = 3
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
            end

            trajectories = reshape(trajectories, [], 1);
            fitOptions = struct( ...
                "psiKnotPoints", {options.psiKnotPoints}, ...
                "psiS", options.psiS, ...
                "fastKnotPoints", {options.fastKnotPoints}, ...
                "fastS", options.fastS, ...
                "mesoscaleConstraint", string(options.mesoscaleConstraint));
            fitOptionsWithoutDecomposition = fitOptions;
            fitOptionsWithoutDecomposition.buildDecomposition = false;

            [queryTimes, scoreTimes, scoreIndices] = GriddedStreamfunctionBootstrap.resolveBootstrapTimes( ...
                trajectories, options.queryTimes, options.scoreTimes, options.scoreStride);

            self = GriddedStreamfunctionBootstrap();
            self.observedTrajectories = trajectories;
            self.observedTrajectorySampleData = GriddedStreamfunction.sampleTrajectoryData(trajectories);
            observedTotalDisplacementWeights = GriddedStreamfunctionBootstrap.totalDisplacementWeightsForSampleTimes( ...
                self.observedTrajectorySampleData.tCell);
            self.nBootstraps = options.nBootstraps;
            self.randomSeed = options.randomSeed;
            self.queryTimes = queryTimes;
            self.scoreTimes = scoreTimes;

            fitOptions.sampleData = self.observedTrajectorySampleData;
            fitArguments = namedargs2cell(fitOptions);
            self.fullFit = GriddedStreamfunction.fromTrajectories(trajectories, fitArguments{:});
            fullSummary = GriddedStreamfunctionBootstrap.extractSummaryFromFit(self.fullFit, queryTimes);
            fullDiagnostics = GriddedStreamfunctionBootstrap.diagnosticsFromFit(self.fullFit, observedTotalDisplacementWeights);
            self.fullSummaryUCenter = fullSummary.uCenter;
            self.fullSummaryVCenter = fullSummary.vCenter;
            self.fullSummarySigmaN = fullSummary.sigma_n;
            self.fullSummarySigmaS = fullSummary.sigma_s;
            self.fullSummaryZeta = fullSummary.zeta;
            self.fullFitKappaValue = fullDiagnostics.kappa;
            self.fullFitCoherenceValue = fullDiagnostics.coherence;
            fullCoherenceSpectrumStorage = GriddedStreamfunctionBootstrap.coherenceSpectrumStorage( ...
                fullDiagnostics.coherenceSpectrum);
            self.fullFitCoherenceFrequency = fullCoherenceSpectrumStorage.frequency;
            self.fullFitCoherenceValues = fullCoherenceSpectrumStorage.coherence;

            nQuery = numel(queryTimes);
            nTrajectories = numel(trajectories);
            self.bootstrapIndices = zeros(self.nBootstraps, nTrajectories);
            self.summaryUCenter = zeros(nQuery, self.nBootstraps);
            self.summaryVCenter = zeros(nQuery, self.nBootstraps);
            self.summarySigmaN = zeros(nQuery, self.nBootstraps);
            self.summarySigmaS = zeros(nQuery, self.nBootstraps);
            self.summaryZeta = zeros(nQuery, self.nBootstraps);
            self.bootstrapFitMetadata = GriddedStreamfunctionBootstrapFitMetadata.empty(0, 1);
            self.bootstrapKappaCache = zeros(1, 0);
            self.bootstrapCoherenceCache = zeros(1, 0);

            originalRng = rng;
            cleanupRng = onCleanup(@() rng(originalRng));
            rng(self.randomSeed);

            for iBootstrap = 1:self.nBootstraps
                sampledIndices = randi(nTrajectories, 1, nTrajectories);
                sampledTrajectories = reshape(trajectories(sampledIndices), [], 1);
                sampledData = GriddedStreamfunction.resampledTrajectorySampleData( ...
                    self.observedTrajectorySampleData, sampledIndices);
                fitOptionsWithoutDecomposition.sampleData = sampledData;
                fitArgumentsWithoutDecomposition = namedargs2cell(fitOptionsWithoutDecomposition);
                fit = GriddedStreamfunction.fromTrajectories(sampledTrajectories, fitArgumentsWithoutDecomposition{:});
                summary = GriddedStreamfunctionBootstrap.extractSummaryFromFit(fit, queryTimes);

                self.bootstrapIndices(iBootstrap, :) = sampledIndices;
                self.summaryUCenter(:, iBootstrap) = summary.uCenter;
                self.summaryVCenter(:, iBootstrap) = summary.vCenter;
                self.summarySigmaN(:, iBootstrap) = summary.sigma_n;
                self.summarySigmaS(:, iBootstrap) = summary.sigma_s;
                self.summaryZeta(:, iBootstrap) = summary.zeta;
                self.bootstrapFitMetadata(iBootstrap, 1) = GriddedStreamfunctionBootstrapFitMetadata.fromFit(fit);
            end

            scores = GriddedStreamfunctionBootstrap.computeConsensusScores(self.summary, scoreIndices, self.fullFit.mesoscaleConstraint);
            self.scoreUv = scores.uv;
            self.scoreStrain = scores.strain;
            self.scoreZeta = scores.zeta;
            self.scoreJoint = scores.joint;

            iBest = self.bestBootstrapIndex();
            bestFit = self.fitForBootstrap(iBest);
            bestDiagnostics = GriddedStreamfunctionBootstrap.diagnosticsFromFit( ...
                bestFit, observedTotalDisplacementWeights(self.bootstrapIndices(iBest, :)));
            self.bestFitKappaValue = bestDiagnostics.kappa;
            self.bestFitCoherenceValue = bestDiagnostics.coherence;
            bestCoherenceSpectrumStorage = GriddedStreamfunctionBootstrap.coherenceSpectrumStorage( ...
                bestDiagnostics.coherenceSpectrum);
            self.bestFitCoherenceFrequency = bestCoherenceSpectrumStorage.frequency;
            self.bestFitCoherenceValues = bestCoherenceSpectrumStorage.coherence;
        end
    end

    methods (Static, Hidden)
        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, { ...
                'observedTrajectories', ...
                'bootstrapIndices', ...
                'queryTimes', ...
                'scoreTimes', ...
                'fullSummaryUCenter', ...
                'fullSummaryVCenter', ...
                'fullSummarySigmaN', ...
                'fullSummarySigmaS', ...
                'fullSummaryZeta', ...
                'summaryUCenter', ...
                'summaryVCenter', ...
                'summarySigmaN', ...
                'summarySigmaS', ...
                'summaryZeta', ...
                'fullFitKappaValue', ...
                'fullFitCoherenceValue', ...
                'fullFitCoherenceFrequency', ...
                'fullFitCoherenceValues', ...
                'bestFitKappaValue', ...
                'bestFitCoherenceValue', ...
                'bestFitCoherenceFrequency', ...
                'bestFitCoherenceValues', ...
                'scoreUv', ...
                'scoreStrain', ...
                'scoreZeta', ...
                'scoreJoint', ...
                'bootstrapFitMetadata', ...
                'nBootstraps', ...
                'randomSeed'});
            vars.fullFit = GriddedStreamfunction.annotatedClassFromGroup(group.groupWithName('fullFit'));
            GriddedStreamfunctionBootstrap.validateCanonicalState(vars);

            self = GriddedStreamfunctionBootstrap();
            self.fullFit = vars.fullFit;
            self.observedTrajectories = reshape(vars.observedTrajectories, [], 1);
            self.bootstrapIndices = vars.bootstrapIndices;
            self.queryTimes = vars.queryTimes;
            self.scoreTimes = vars.scoreTimes;
            self.fullSummaryUCenter = vars.fullSummaryUCenter;
            self.fullSummaryVCenter = vars.fullSummaryVCenter;
            self.fullSummarySigmaN = vars.fullSummarySigmaN;
            self.fullSummarySigmaS = vars.fullSummarySigmaS;
            self.fullSummaryZeta = vars.fullSummaryZeta;
            self.summaryUCenter = vars.summaryUCenter;
            self.summaryVCenter = vars.summaryVCenter;
            self.summarySigmaN = vars.summarySigmaN;
            self.summarySigmaS = vars.summarySigmaS;
            self.summaryZeta = vars.summaryZeta;
            self.fullFitKappaValue = vars.fullFitKappaValue;
            self.fullFitCoherenceValue = vars.fullFitCoherenceValue;
            self.fullFitCoherenceFrequency = vars.fullFitCoherenceFrequency;
            self.fullFitCoherenceValues = vars.fullFitCoherenceValues;
            self.bestFitKappaValue = vars.bestFitKappaValue;
            self.bestFitCoherenceValue = vars.bestFitCoherenceValue;
            self.bestFitCoherenceFrequency = vars.bestFitCoherenceFrequency;
            self.bestFitCoherenceValues = vars.bestFitCoherenceValues;
            self.scoreUv = reshape(vars.scoreUv, 1, []);
            self.scoreStrain = reshape(vars.scoreStrain, 1, []);
            self.scoreZeta = reshape(vars.scoreZeta, 1, []);
            self.scoreJoint = reshape(vars.scoreJoint, 1, []);
            self.bootstrapFitMetadata = reshape(vars.bootstrapFitMetadata, [], 1);
            self.nBootstraps = vars.nBootstraps;
            self.randomSeed = vars.randomSeed;
            self.refreshObservedTrajectorySampleData();
            optionalVars = CAAnnotatedClass.propertyValuesFromGroup(group, ...
                {'bootstrapKappaCache', 'bootstrapCoherenceCache'}, shouldIgnoreMissingProperties=true);
            if isstruct(optionalVars)
                if isfield(optionalVars, 'bootstrapKappaCache')
                    self.bootstrapKappaCache = reshape(optionalVars.bootstrapKappaCache, 1, []);
                end
                if isfield(optionalVars, 'bootstrapCoherenceCache')
                    self.bootstrapCoherenceCache = reshape(optionalVars.bootstrapCoherenceCache, 1, []);
                end
            end
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0, 0);
            propertyAnnotations(end+1) = CADimensionProperty('queryTimes', '', 'Times used to store COM-local bootstrap summaries.');
            propertyAnnotations(end+1) = CADimensionProperty('scoreTimes', '', 'Times used to compute the consensus score.');
            propertyAnnotations(end+1) = CADimensionProperty('bootstrapIndex', '', 'Index over bootstrap replicates.');
            propertyAnnotations(end+1) = CADimensionProperty('trajectoryIndex', '', 'Index over observed drifter trajectories.');
            propertyAnnotations(end+1) = CADimensionProperty('fullFitCoherenceFrequencyIndex', '', 'Index over the full-fit coherence-spectrum frequencies.');
            propertyAnnotations(end+1) = CADimensionProperty('bestFitCoherenceFrequencyIndex', '', 'Index over the best-fit coherence-spectrum frequencies.');
            propertyAnnotations(end+1) = CAObjectProperty('fullFit', 'Full-data gridded-streamfunction fit used as the reference solution.', className='GriddedStreamfunction', sizeText='(1,1)');
            propertyAnnotations(end+1) = CAObjectProperty('observedTrajectories', 'Original drifter trajectories used to seed the bootstrap ensemble.', className='TrajectorySpline', sizeText='nonempty vector');
            propertyAnnotations(end+1) = CANumericProperty('bootstrapIndices', {'bootstrapIndex', 'trajectoryIndex'}, '', 'Whole-drifter resampling indices for each bootstrap replicate.');
            propertyAnnotations(end+1) = CANumericProperty('fullSummaryUCenter', {'queryTimes'}, '', 'Persisted full-fit `uCenter` summary values.');
            propertyAnnotations(end+1) = CANumericProperty('fullSummaryVCenter', {'queryTimes'}, '', 'Persisted full-fit `vCenter` summary values.');
            propertyAnnotations(end+1) = CANumericProperty('fullSummarySigmaN', {'queryTimes'}, '', 'Persisted full-fit `sigma_n` summary values.');
            propertyAnnotations(end+1) = CANumericProperty('fullSummarySigmaS', {'queryTimes'}, '', 'Persisted full-fit `sigma_s` summary values.');
            propertyAnnotations(end+1) = CANumericProperty('fullSummaryZeta', {'queryTimes'}, '', 'Persisted full-fit `zeta` summary values.');
            propertyAnnotations(end+1) = CANumericProperty('summaryUCenter', {'queryTimes', 'bootstrapIndex'}, '', 'Persisted bootstrap `uCenter` summary matrix.');
            propertyAnnotations(end+1) = CANumericProperty('summaryVCenter', {'queryTimes', 'bootstrapIndex'}, '', 'Persisted bootstrap `vCenter` summary matrix.');
            propertyAnnotations(end+1) = CANumericProperty('summarySigmaN', {'queryTimes', 'bootstrapIndex'}, '', 'Persisted bootstrap `sigma_n` summary matrix.');
            propertyAnnotations(end+1) = CANumericProperty('summarySigmaS', {'queryTimes', 'bootstrapIndex'}, '', 'Persisted bootstrap `sigma_s` summary matrix.');
            propertyAnnotations(end+1) = CANumericProperty('summaryZeta', {'queryTimes', 'bootstrapIndex'}, '', 'Persisted bootstrap `zeta` summary matrix.');
            propertyAnnotations(end+1) = CANumericProperty('fullFitKappaValue', {}, '', 'Persisted full-fit scalar diffusivity diagnostic.');
            propertyAnnotations(end+1) = CANumericProperty('fullFitCoherenceValue', {}, '', 'Persisted full-fit scalar mean coherence.');
            propertyAnnotations(end+1) = CANumericProperty('fullFitCoherenceFrequency', {'fullFitCoherenceFrequencyIndex'}, '', 'Persisted full-fit coherence-spectrum frequencies.');
            propertyAnnotations(end+1) = CANumericProperty('fullFitCoherenceValues', {'fullFitCoherenceFrequencyIndex'}, '', 'Persisted full-fit coherence-spectrum values.');
            propertyAnnotations(end+1) = CANumericProperty('bestFitKappaValue', {}, '', 'Persisted best-fit scalar diffusivity diagnostic.');
            propertyAnnotations(end+1) = CANumericProperty('bestFitCoherenceValue', {}, '', 'Persisted best-fit scalar mean coherence.');
            propertyAnnotations(end+1) = CANumericProperty('bestFitCoherenceFrequency', {'bestFitCoherenceFrequencyIndex'}, '', 'Persisted best-fit coherence-spectrum frequencies.');
            propertyAnnotations(end+1) = CANumericProperty('bestFitCoherenceValues', {'bestFitCoherenceFrequencyIndex'}, '', 'Persisted best-fit coherence-spectrum values.');
            propertyAnnotations(end+1) = CANumericProperty('scoreUv', {'bootstrapIndex'}, '', 'Persisted bootstrap velocity consensus scores.');
            propertyAnnotations(end+1) = CANumericProperty('scoreStrain', {'bootstrapIndex'}, '', 'Persisted bootstrap strain consensus scores.');
            propertyAnnotations(end+1) = CANumericProperty('scoreZeta', {'bootstrapIndex'}, '', 'Persisted bootstrap vorticity consensus scores.');
            propertyAnnotations(end+1) = CANumericProperty('scoreJoint', {'bootstrapIndex'}, '', 'Persisted bootstrap joint consensus scores.');
            propertyAnnotations(end+1) = CAObjectProperty('bootstrapFitMetadata', 'Persisted exact reconstruction metadata for each bootstrap replicate.');
            propertyAnnotations(end+1) = CANumericProperty('nBootstraps', {}, '', 'Number of bootstrap replicates in the ensemble.');
            propertyAnnotations(end+1) = CANumericProperty('randomSeed', {}, '', 'Random-number seed used for whole-drifter resampling.');
            propertyAnnotations(end+1) = CANumericProperty('bootstrapKappaCache', {'bootstrapIndex'}, '', 'Optional cached bootstrap scalar diffusivity diagnostics.');
            propertyAnnotations(end+1) = CANumericProperty('bootstrapCoherenceCache', {'bootstrapIndex'}, '', 'Optional cached bootstrap scalar coherence diagnostics.');
        end

        function names = classRequiredPropertyNames()
            names = { ...
                'fullFit', ...
                'observedTrajectories', ...
                'bootstrapIndices', ...
                'queryTimes', ...
                'scoreTimes', ...
                'fullSummaryUCenter', ...
                'fullSummaryVCenter', ...
                'fullSummarySigmaN', ...
                'fullSummarySigmaS', ...
                'fullSummaryZeta', ...
                'summaryUCenter', ...
                'summaryVCenter', ...
                'summarySigmaN', ...
                'summarySigmaS', ...
                'summaryZeta', ...
                'fullFitKappaValue', ...
                'fullFitCoherenceValue', ...
                'fullFitCoherenceFrequency', ...
                'fullFitCoherenceValues', ...
                'bestFitKappaValue', ...
                'bestFitCoherenceValue', ...
                'bestFitCoherenceFrequency', ...
                'bestFitCoherenceValues', ...
                'scoreUv', ...
                'scoreStrain', ...
                'scoreZeta', ...
                'scoreJoint', ...
                'bootstrapFitMetadata', ...
                'nBootstraps', ...
                'randomSeed'};
        end
    end

    methods (Static, Access = private)
        function validateCanonicalState(options)
            if isempty(options.fullFit) && isempty(options.observedTrajectories) && isempty(options.bootstrapIndices) && ...
                    isempty(options.queryTimes) && isempty(options.scoreTimes) && ...
                    isempty(options.fullSummaryUCenter) && isempty(options.fullSummaryVCenter) && ...
                    isempty(options.fullSummarySigmaN) && isempty(options.fullSummarySigmaS) && ...
                    isempty(options.fullSummaryZeta) && isempty(options.summaryUCenter) && ...
                    isempty(options.summaryVCenter) && isempty(options.summarySigmaN) && ...
                    isempty(options.summarySigmaS) && isempty(options.summaryZeta) && ...
                    isnan(options.fullFitKappaValue) && isnan(options.fullFitCoherenceValue) && ...
                    isempty(options.fullFitCoherenceFrequency) && isempty(options.fullFitCoherenceValues) && ...
                    isnan(options.bestFitKappaValue) && isnan(options.bestFitCoherenceValue) && ...
                    isempty(options.bestFitCoherenceFrequency) && isempty(options.bestFitCoherenceValues) && ...
                    isempty(options.scoreUv) && isempty(options.scoreStrain) && ...
                    isempty(options.scoreZeta) && isempty(options.scoreJoint) && ...
                    isempty(options.bootstrapFitMetadata) && ...
                    options.nBootstraps == 0 && options.randomSeed == 0
                return
            end

            if ~(isa(options.fullFit, 'GriddedStreamfunction') && isscalar(options.fullFit))
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'fullFit must be a scalar GriddedStreamfunction.');
            end
            if ~(isa(options.observedTrajectories, 'TrajectorySpline') && isvector(options.observedTrajectories) && ~isempty(options.observedTrajectories))
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'observedTrajectories must be a nonempty vector of TrajectorySpline objects.');
            end
            if ~(isa(options.bootstrapFitMetadata, 'GriddedStreamfunctionBootstrapFitMetadata') && isvector(options.bootstrapFitMetadata))
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'bootstrapFitMetadata must be a vector of GriddedStreamfunctionBootstrapFitMetadata objects.');
            end
            if any(round(options.bootstrapIndices) ~= options.bootstrapIndices, 'all') || any(options.bootstrapIndices < 1, 'all')
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'bootstrapIndices must contain positive integer drifter indices.');
            end

            nBootstraps = options.nBootstraps;
            nQuery = numel(options.queryTimes);
            if size(options.bootstrapIndices, 1) ~= nBootstraps
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'bootstrapIndices must have one row per bootstrap replicate.');
            end
            if size(options.summaryUCenter, 1) ~= nQuery || size(options.summaryUCenter, 2) ~= nBootstraps
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'summaryUCenter must have size [numel(queryTimes) nBootstraps].');
            end

            summaryMatrices = {options.summaryVCenter, options.summarySigmaN, options.summarySigmaS, options.summaryZeta};
            for iMatrix = 1:numel(summaryMatrices)
                if ~isequal(size(summaryMatrices{iMatrix}), size(options.summaryUCenter))
                    error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                        'All persisted bootstrap summary matrices must have the same size.');
                end
            end

            fullSummaryVectors = { ...
                options.fullSummaryUCenter, ...
                options.fullSummaryVCenter, ...
                options.fullSummarySigmaN, ...
                options.fullSummarySigmaS, ...
                options.fullSummaryZeta};
            for iVector = 1:numel(fullSummaryVectors)
                if numel(fullSummaryVectors{iVector}) ~= nQuery
                    error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                        'Each persisted full-fit summary vector must align with queryTimes.');
                end
            end

            if numel(options.fullFitCoherenceFrequency) ~= numel(options.fullFitCoherenceValues) || ...
                    numel(options.bestFitCoherenceFrequency) ~= numel(options.bestFitCoherenceValues)
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'Persisted coherence frequency and coherence arrays must have matching lengths.');
            end

            if numel(options.scoreUv) ~= nBootstraps || numel(options.scoreStrain) ~= nBootstraps || ...
                    numel(options.scoreZeta) ~= nBootstraps || numel(options.scoreJoint) ~= nBootstraps || ...
                    numel(options.bootstrapFitMetadata) ~= nBootstraps
                error('GriddedStreamfunctionBootstrap:InvalidCanonicalState', ...
                    'Persisted bootstrap scores and metadata must all have length nBootstraps.');
            end
        end

        [queryTimes, scoreTimes, scoreIndices] = resolveBootstrapTimes(trajectories, queryTimesOption, scoreTimesOption, scoreStride)
        summary = extractSummaryFromFit(fit, queryTimes)
        diagnostics = diagnosticsFromFit(fit, totalDisplacementWeightsByTrajectory)
        kappaEstimate = kappaEstimateFromFit(fit, totalDisplacementWeightsByTrajectory)
        spectrum = coherenceSpectrumFromStorage(frequency, coherence)
        spectrum = coherenceSpectrumStorage(spectrum)
        totalDisplacementWeightsByTrajectory = totalDisplacementWeightsForSampleTimes(tCell)
        [xTotalByTrajectory, yTotalByTrajectory, durationByTrajectory] = submesoscaleTotalDisplacements(sampleData, totalDisplacementWeightsByTrajectory)
        scores = computeConsensusScores(summary, scoreIndices, mesoscaleConstraint)
    end

    methods (Access = private)
        ensureBootstrapDiagnostics(self)
        function refreshObservedTrajectorySampleData(self)
            if isempty(self.observedTrajectories)
                self.observedTrajectorySampleData = struct([]);
            else
                self.observedTrajectorySampleData = GriddedStreamfunction.sampleTrajectoryData(self.observedTrajectories);
            end
        end
    end
end
