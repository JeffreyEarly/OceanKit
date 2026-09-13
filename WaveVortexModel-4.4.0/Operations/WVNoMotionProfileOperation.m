classdef WVNoMotionProfileOperation < WVOperation
    properties (GetAccess=public, SetAccess=protected)
        solver
        % Most recent compute result, including solver identity and exit flag.
        %
        % - Topic: Internal
        % - Developer: true
        lastSolverOutput = struct()
    end

    methods

        function self = WVNoMotionProfileOperation(options)
            arguments
                options.solver (1,1) string {mustBeMember(options.solver,["auto","lsqnonlin","fminsearch","dampedLeastSquares"])} = "auto"
            end
            outputVariables(1) = WVVariableAnnotation('rho_nm',{'z'},'kg m-3', 'no-motion density profile');
            % outputVariables(1).isVariableWithLinearTimeStep = false;
            self@WVOperation('rho_nm',outputVariables,@disp);

            if options.solver ~= "auto"
                self.solver = options.solver;
            else
                self.solver = "dampedLeastSquares";
            end
        end

        function varargout = compute(self,wvt,varargin)
            self.lastSolverOutput = struct();
            if wvt.isCompiledBuiltinOperation(self) && self.solver == "dampedLeastSquares"
                try
                    varargout = wvt.compiledVariables({'rho_nm'});
                catch exception
                    self.lastSolverOutput = wvt.compiledDensityRecoveryReport();
                    rethrow(exception)
                end
                self.lastSolverOutput = wvt.compiledDensityRecoveryReport();
                return
            end
            [rho_nm,exitflag,output] = WVNoMotionProfileOperation.find_rho_nm(wvt.z_int, wvt.Lz, wvt.rho_total, wvt.rho_nm0,solver=self.solver);
            output.exitflag = exitflag;
            output.solver = self.solver;
            self.lastSolverOutput = output;
            if self.solver == "dampedLeastSquares" && (exitflag <= 0 || output.maximumResidual > 1e-8 || any(~isfinite(rho_nm)) || any(diff(rho_nm) >= 0))
                error('WVNoMotionProfileOperation:UnqualifiedFit','The no-motion fit did not reach a finite, strictly ordered profile with normalized moment residual at most 1e-8. Inspect lastSolverOutput before using these density diagnostics.');
            end
            varargout = {rho_nm};
        end

    end

    methods (Static)
        function tf = hasOptimizationToolboxSupport()
            tf = license("test","Optimization_Toolbox") ...
                && exist("lsqnonlin","file")==2 ...
                && exist("optimoptions","file")==2;
        end

        function [rho,exitflag,output] = find_rho_nm(z_int, Lz, rho_total, rho_nm0,options)
            % Fit a stably ordered profile to volume-weighted density moments.
            %
            % Horizontally uniform stable density is returned exactly. For
            % other states, the default solver uses the extrema of the current
            % density field. Explicit legacy solvers retain the initial
            % profile endpoints. Interior nodes use strictly ordered log-gap
            % parameters; residuals are normalized raw moments 1 through Nz.
            %
            % lsqnonlin requires Optimization Toolbox; fminsearch retains the
            % legacy toolbox-free fit. dampedLeastSquares uses augmented QR
            % with explicit bounded iteration/evaluation controls. Direct
            % callers must inspect exitflag and output: a small step alone
            % does not establish an accurate physical density distribution.
            %
            % - Topic: Internal
            % - Developer: true

            arguments
                z_int (:,1) double {mustBeFinite, mustBePositive}
                Lz (1,1) double {mustBeFinite, mustBePositive}
                rho_total (:,:,:) double {mustBeFinite}
                rho_nm0 (:,1) double {mustBeFinite, mustBeNonnegative}
                options.solver (1,1) string {mustBeMember(options.solver, ["lsqnonlin","fminsearch","dampedLeastSquares"])} = "dampedLeastSquares"
            end

            n = numel(z_int);
            if numel(rho_nm0) ~= n || size(rho_total,3) ~= n
                error('WVNoMotionProfileOperation:InvalidDimensions','Density, reference profile and integration weights must share the same vertical dimension.');
            end
            if any(diff(rho_nm0) >= 0)
                error('WVNoMotionProfileOperation:NonInvertibleReference','The initial reference density must strictly decrease with height.');
            end

            % A horizontally uniform stable state is already its own
            % no-motion profile, even after its extrema have changed.
            % Recover it exactly instead of solving an ill-conditioned
            % inverse moment problem against the old endpoints.
            stableProfile = reshape(rho_total(1,1,:),[],1);
            if all(diff(stableProfile) < 0) && all(rho_total == reshape(stableProfile,1,1,[]),'all')
                rho = stableProfile;
                exitflag = 1;
                output = struct(iterations=0,funcCount=0,maximumResidual=0,algorithm="stable-rest-profile",message="The density is already horizontally uniform and stably ordered.");
                return
            end

            z = flip(z_int/Lz);

            rho0 = rho_nm0(end);
            rhoD = rho_nm0(1);

            % Use the initial profile's normalized shape as the initial guess.
            u0 = WVNoMotionProfileOperation.rho_to_u(rho_nm0, rho0, rhoD);
            if options.solver == "dampedLeastSquares"
                rho0 = min(rho_total,[],"all");
                rhoD = max(rho_total,[],"all");
                if rhoD <= rho0
                    error('WVNoMotionProfileOperation:NonInvertibleDistribution','The current density distribution is constant and has no unique inverse material height.');
                end
            end
            rho_moment = WVNoMotionProfileOperation.moments_from_rho_tot(rho_total, rho0, rhoD, z_int, Lz);

            switch options.solver
                case "dampedLeastSquares"
                    [u,exitflag,output] = WVNoMotionProfileOperation.solveMoments(u0,z,rho_moment);
                case "lsqnonlin"
                    options = optimoptions("lsqnonlin", ...
                        "Display","none", ...
                        "Algorithm","trust-region-reflective", ...
                        "SpecifyObjectiveGradient",true, ...
                        "FunctionTolerance",1e-12, ...
                        "StepTolerance",1e-12, ...
                        "OptimalityTolerance",1e-12, ...
                        "MaxFunctionEvaluations",5e4, ...
                        "MaxIterations",2e3);

                    % Solve (overdetermined): n residuals, n-2 parameters
                    fun = @(u) WVNoMotionProfileOperation.residual_and_jacobian(u,z,rho_moment);
                    [u,~,~,exitflag,output] = lsqnonlin(fun, u0, [], [], options);
                case "fminsearch"
                    phi = @(u) sum(WVNoMotionProfileOperation.residual(u,z,rho_moment).^2);

                    fmopts = optimset( ...
                        "Display","notify", ...
                        "TolFun", 1e-12, ...
                        "TolX", 1e-5, ...
                        "MaxIter", 10e3, ...
                        "MaxFunEvals", 5e4);

                    [u,~,exitflag,output] = fminsearch(phi, u0, fmopts);
            end
            % Recover m
            rho = WVNoMotionProfileOperation.u_to_rho(u, rho0, rhoD);
        end

        % only issue now is that m0_j is decreasing, not increasing.

        function [u,exitflag,output] = solveMoments(u0,weights,target,options)
            % Fit normalized density moments with bounded damped least squares.
            % The augmented rectangular system is solved directly, avoiding normal
            % equations for the ill-conditioned raw-moment Jacobian.
            %
            % - Topic: Internal
            % - Developer: true
            arguments
                u0 (:,1) double {mustBeFinite,mustBeReal}
                weights (:,1) double {mustBeFinite,mustBeReal,mustBePositive}
                target (:,1) double {mustBeFinite,mustBeReal}
                options.maximumIterations (1,1) double {mustBeFinite,mustBeInteger,mustBeNonnegative} = 2000
                options.maximumEvaluations (1,1) double {mustBeFinite,mustBeInteger,mustBePositive} = 5000
                options.gradientTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-12
                options.stepTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-12
                options.relativeCostTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-12
            end
            if numel(weights) ~= numel(u0)+2 || numel(target) ~= numel(weights)
                error("WVNoMotionProfileOperation:InvalidDimensions","Use n weights and target moments with n-2 log-gap parameters.");
            end
            u = u0;
            [residual,J] = WVNoMotionProfileOperation.residual_and_jacobian(u,weights,target);
            cost = sum(residual.^2)/2;
            initialCost = cost;
            scale = max(sum(J.^2,1));
            damping = max(1e-3*scale,realmin);
            multiplier = 2;
            accepted = 0;
            rejected = 0;
            evaluations = 1;
            exitflag = 0;
            reason = "iteration-limit";
            iterations = 0;
            identity = eye(numel(u));
            for iteration = 1:options.maximumIterations
                iterations = iteration;
                gradientNorm = norm(J'*residual,Inf);
                if gradientNorm<=options.gradientTolerance
                    exitflag = 1;
                    reason = "gradient-tolerance";
                    break
                end
                if evaluations >= options.maximumEvaluations
                    reason = "evaluation-limit";
                    break
                end
                step = [J;sqrt(damping)*identity]\[-residual;zeros(numel(u),1)];
                if any(~isfinite(step))
                    reason = "nonfinite-step";
                    exitflag = -1;
                    break
                end
                if norm(step)<=options.stepTolerance*(norm(u)+options.stepTolerance)
                    exitflag = 2;
                    reason = "step-tolerance";
                    break
                end
                trial = u+step;
                [trialResidual,trialJ] = WVNoMotionProfileOperation.residual_and_jacobian(trial,weights,target);
                evaluations = evaluations+1;
                trialCost = sum(trialResidual.^2)/2;
                linearChange = J*step;
                predictedDecrease = -residual'*linearChange-sum(linearChange.^2)/2;
                actualDecrease = cost-trialCost;
                if isfinite(trialCost) && all(isfinite(trialJ),"all") && predictedDecrease>0 && actualDecrease>0
                    gain = actualDecrease/predictedDecrease;
                    oldCost = cost;
                    u = trial;
                    residual = trialResidual;
                    J = trialJ;
                    cost = trialCost;
                    accepted = accepted+1;
                    damping = max(realmin,damping*max(1/3,1-(2*gain-1)^3));
                    multiplier = 2;
                    if actualDecrease<=options.relativeCostTolerance*oldCost
                        exitflag = 3;
                        reason = "relative-cost-tolerance";
                        break
                    end
                else
                    rejected = rejected+1;
                    damping = damping*multiplier;
                    multiplier = 2*multiplier;
                    if ~isfinite(damping) || damping>1e30*max(scale,realmin)
                        exitflag = -2;
                        reason = "damping-limit";
                        break
                    end
                end
                if evaluations>=options.maximumEvaluations
                    reason = "evaluation-limit";
                    break
                end
            end
            output = struct(algorithm="damped-least-squares",reason=reason,iterations=iterations,evaluations=evaluations,acceptedSteps=accepted, ...
                rejectedSteps=rejected,initialCost=initialCost,finalCost=cost,maximumResidual=max(abs(residual)), ...
                gradientNorm=norm(J'*residual,Inf),damping=damping,parameters=options);
        end

        function rho_moment = moments_from_rho_tot(rho_total, rho0, rhoD, z_int, Lz)
            int_vol_avg = @(integrand) sum(mean(mean(shiftdim(z_int,-2).*integrand,1),2),3)/Lz;

            normalizedDensity = (rho_total-rho0)/(rhoD-rho0);
            power = normalizedDensity;
            rho_moment = zeros(size(rho_total,3),1);
            for i=1:size(rho_total,3)
                rho_moment(i) = int_vol_avg(power);
                if i < size(rho_total,3)
                    power = power.*normalizedDensity;
                end
            end
        end

        function u = rho_to_u(rho, rho0, rhoD)
            % Convert an initial guess of rho, create an intermediate with (m1=0, mn=M) and map to u=[u3..un].
            n = numel(rho);
            assert(n >= 3, 'Need n >= 3 to have at least one free interior node.');

            % map function to [0 1], monotonically increasing
            m = flip((rho - rho0)/(rhoD - rho0));

            x = log(m(2:end)); % x2..xn (length n-1), xn should be log(1)
            assert(abs(x(end)) == 0, 'the end point must be 1 exactly.');

            g = diff(x);               % gaps: g3..gn (length n-2), must be > 0
            assert(all(g>0), 'Interior nodes must be strictly increasing and < M.');
            u = log(g(:));             % unconstrained parameters
        end

        function rho = u_to_rho(u, rho0, rhoD)
            % Convert u=[u3..un] to full m with m1=0, mn=M.
            u = u(:);
            np = numel(u);          % np = n-2
            n = np + 2;

            g = exp(u);             % g3..gn > 0

            % Build x2..xn with xn fixed to U
            x2 = - sum(g);        % ensures xn = U
            x = zeros(n-1,1);       % x2..xn
            x(1) = x2;              % x2
            if n-1 > 1
                x(2:end) = x2 + cumsum(g);   % x3..xn
            end

            m = zeros(n,1);
            m(1) = 0;
            m(2:n-1) = exp(x(1:end-1));  % m2..m_{n-1}
            m(n) = 1;

            rho = flip((rhoD - rho0)*m + rho0);
        end

        function r = residual(u,z,c)
            % Residuals:
            %   r_k = sum_{j=2..n-1} z_j m_j^k + z_n M^k - c_k
            % since m1=0 contributes 0 for k>=1.
            %
            % u parameterizes log-gaps between x2..xn with xn fixed.

            u = u(:);
            z = z(:);
            c = c(:);

            n = numel(z);

            % Build x2..x_{n-1} from u
            g = exp(u);              % length np = gaps g3..gn
            x2 =- sum(g);

            % x2..xn in a vector of length n-1, but we only need x2..x_{n-1}
            x_pos = zeros(n-1,1);    % x2..xn
            x_pos(1) = x2;
            if n-1 > 1
                x_pos(2:end) = x2 + cumsum(g);   % x3..xn
            end
            x2_to_nminus1 = x_pos(1:end-1);      % x2..x_{n-1} length n-2

            % Compute residuals for k=1..n
            kvec = (1:n).';                        % (n,1)

            % Interior contribution: sum_{j=2..n-1} z_j exp(k*x_j)
            E = exp(kvec * (x2_to_nminus1.'));     % (n, n-2), columns correspond to j=2..n-1
            W = E .* (ones(n,1) * z(2:n-1).');     % (n, n-2)
            interior = sum(W,2);

            % Endpoint contribution from fixed mn = M
            r = interior + z(n) - c;
        end

        function [r,J] = residual_and_jacobian(u,z,c)
            % Residuals:
            %   r_k = sum_{j=2..n-1} z_j m_j^k + z_n M^k - c_k
            % since m1=0 contributes 0 for k>=1.
            %
            % u parameterizes log-gaps between x2..xn with xn fixed.

            u = u(:);
            z = z(:);
            c = c(:);

            n = numel(z);

            % Build x2..x_{n-1} from u
            g = exp(u);              % length np = gaps g3..gn
            x2 =- sum(g);

            % x2..xn in a vector of length n-1, but we only need x2..x_{n-1}
            x_pos = zeros(n-1,1);    % x2..xn
            x_pos(1) = x2;
            if n-1 > 1
                x_pos(2:end) = x2 + cumsum(g);   % x3..xn
            end
            x2_to_nminus1 = x_pos(1:end-1);      % x2..x_{n-1} length n-2

            % Compute residuals for k=1..n
            kvec = (1:n).';                        % (n,1)

            % Interior contribution: sum_{j=2..n-1} z_j exp(k*x_j)
            E = exp(kvec * (x2_to_nminus1.'));     % (n, n-2), columns correspond to j=2..n-1
            W = E .* (ones(n,1) * z(2:n-1).');     % (n, n-2)
            interior = sum(W,2);

            % Endpoint contribution from fixed mn = M
            r = interior + z(n) - c;

            % Jacobian J is n-by-(n-2): columns correspond to u_i = u_{i+3-?}
            % Here u index i=1..np corresponds to gap g_{j} with j = i+2? (since g3 is first)
            % Key derivative:
            %   dr_k/dx_j = k z_j exp(k x_j)
            % and due to the "keep xn fixed" construction:
            %   ∂x_j/∂u_i = -g_i   for j <= (i+1) in the x2..x_{n-1} indexing
            %              0       for later j
            %
            % More explicitly: u_i corresponds to gap g_{i+2} (i=1 -> g3).
            % It affects x2..x_{(i+2)-1} = x2..x_{i+1} (in original indices).
            %
            % We'll compute A(k,ell)=dr_k/dx_{(ell+1)} for ell=1..n-2 corresponding to x2..x_{n-1}.
            A = (kvec * ones(1,n-2)) .* W;         % (n, n-2)

            J = zeros(n, n-2);

            % For each i, sum A over the affected x's and multiply by (-g_i)
            % i affects the first i columns of A (since i=1 affects x2 only, i=2 affects x2..x3, etc.)
            prefix = cumsum(A, 2);                  % prefix(:,i) = sum_{ell=1..i} A(:,ell)

            for i = 1:(n-2)
                J(:,i) = -g(i) * prefix(:,i);
            end

        end
    end
end
