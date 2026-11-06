classdef IMHydrostaticRootTransformAssessment
    % Record how a hydrostatic root transform's retained band was qualified.
    %
    % `IMInternalModesDiscreteTransform.hydrostaticRoot` returns this immutable record with
    % its transform and grid. Each row of `search` describes one fitted
    % candidate band $$r$$ (aligned columns $$F_0,\ldots,F_r$$): the F and G
    % energy errors
    %
    % $$\varepsilon_V=\left\|H_V^{-1/2}(\Phi_V^\mathsf{T}W_V\Phi_V-H_V)H_V^{-1/2}\right\|_2,$$
    %
    % sampled Gram ranks and condition numbers, the weight range and the
    % number of zero weights, the largest normalized modal mean
    % $$|\sum_i w_iF_j(z_i)|/D$$, and the omitted-mode gains
    % $$\|H_r^{1/2}A_V\Phi_{V,\mathrm{om}}H_\mathrm{om}^{-1/2}\|_2$$ of the
    % selected forward operator $$A_V$$: the
    % largest retained energy amplitude produced per unit energy amplitude of
    % the next `omittedModeCount` internal modes. Gains depend on that finite
    % omitted set and are diagnostics, not acceptance criteria.
    %
    % ```matlab
    % [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
    % assessment.search
    % assessment.limitingCriterion
    % ```
    %
    % - Topic: Assess transform quality
    % - Declaration: classdef IMHydrostaticRootTransformAssessment

    properties (SetAccess = private)
        % `"strict"` for a requested band or `"automatic"` for a searched band.
        %
        % - Topic: Assess transform quality
        selection

        % Number of sample points, including both endpoints.
        %
        % - Topic: Assess transform quality
        nPoints

        % Retained internal count $$r$$; the aligned family has $$r+1$$ columns.
        %
        % - Topic: Assess transform quality
        retainedInternalModeCount

        % Acceptance tolerance applied to both $$\varepsilon_F$$ and $$\varepsilon_G$$.
        %
        % - Topic: Assess transform quality
        gramTolerance

        % Why a larger band was not returned.
        %
        % `"requested"` for strict selection; `"maximumRootGridBand"` when
        % the largest band $$N_z-2$$ was accepted; otherwise the rejection
        % reason of band $$r+1$$: `"energyErrorF"`, `"energyErrorG"`,
        % `"energyErrorFG"`, `"rankF"` or `"rankG"`.
        %
        % - Topic: Assess transform quality
        limitingCriterion

        % One row per fitted candidate band, in search order.
        %
        % - Topic: Assess transform quality
        search

        % The `search` row of the returned band.
        %
        % - Topic: Assess transform quality
        final

        % Forward operator used by the transform and the omitted-mode gains.
        %
        % `"truncatedCollocation"` or `"galerkin"`; see
        % `IMInternalModesDiscreteTransform.hydrostaticRoot`.
        %
        % - Topic: Assess transform quality
        forwardOperator

        % Square root-grid collocation diagnostics.
        %
        % Condition numbers of the column-equilibrated square F and interior G
        % systems and their spans. Computed for both forward operators.
        %
        % - Topic: Assess transform quality
        collocation

        % Independent refinement evidence for the solved modes and roots.
        %
        % Fields record both resolutions, the largest relative equivalent-depth
        % and shape errors, the largest root displacement relative to the
        % depth, the tolerance, and the full `assessModeConvergence` table.
        %
        % - Topic: Assess transform quality
        convergence
    end

    methods
        function self = IMHydrostaticRootTransformAssessment(options)
            % Create the assessment record.
            %
            % - Topic: Assess transform quality
            % - Declaration: assessment = IMHydrostaticRootTransformAssessment(options)
            % - Developer: true
            arguments
                options.selection (1,1) string {mustBeMember(options.selection,["strict","automatic"])}
                options.nPoints (1,1) double {mustBeInteger,mustBePositive}
                options.retainedInternalModeCount (1,1) double {mustBeInteger,mustBeNonnegative}
                options.gramTolerance (1,1) double {mustBePositive}
                options.limitingCriterion (1,1) string
                options.search table
                options.final table
                options.convergence (1,1) struct
                options.forwardOperator (1,1) string {mustBeMember(options.forwardOperator,["truncatedCollocation","galerkin"])}
                options.collocation (1,1) struct
            end
            self.selection = options.selection;
            self.nPoints = options.nPoints;
            self.retainedInternalModeCount = options.retainedInternalModeCount;
            self.gramTolerance = options.gramTolerance;
            self.limitingCriterion = options.limitingCriterion;
            self.search = options.search;
            self.final = options.final;
            self.convergence = options.convergence;
            self.forwardOperator = options.forwardOperator;
            self.collocation = options.collocation;
        end
    end
end
