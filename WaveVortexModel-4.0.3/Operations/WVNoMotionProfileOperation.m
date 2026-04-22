classdef WVNoMotionProfileOperation < WVOperation
    properties (GetAccess=public, SetAccess=protected)
        solver
    end

    methods

        function self = WVNoMotionProfileOperation()
            arguments
            end
            outputVariables(1) = WVVariableAnnotation('rho_nm',{'z'},'kg m^{-3}', 'no-motion density profile');
            % outputVariables(1).isVariableWithLinearTimeStep = false;
            self@WVOperation('rho_nm',outputVariables,@disp);

            if license("test","Optimization_Toolbox") && exist("fsolve","file")==2
                self.solver = "fsolve";
            else
                self.solver = "fminsearch";
            end
        end

        function varargout = compute(self,wvt,varargin)
            rho_nm = WVNoMotionProfileOperation.find_rho_nm(wvt.z_int, wvt.Lz, wvt.rho_total, wvt.rho_nm0,solver=self.solver);
            varargout = {rho_nm};
        end

    end

    methods (Static)
        function [rho,exitflag,output] = find_rho_nm(z_int, Lz, rho_total, rho_nm0,options)
            %SOLVEMOMENTS_M1ZERO_MNFIXED
            % Solve moments:  sum_{j=1}^n z_j m_j^k = c_k,  k=1..n
            % with constraints:
            %   m(1) = 0 (fixed)
            %   m(n) = M (fixed, known, positive)
            %   0 < m(2) < ... < m(n-1) < M  (strict ordering enforced)
            %
            % Unknowns are u = [u3; ...; un] where gaps in log-space are
            %   g_j = exp(u_j) = x_j - x_{j-1} > 0,  j=3..n
            % and x_n = log(M) is fixed. Then x_2 is chosen so x_n stays fixed:
            %   x_2 = log(M) - sum_{j=3}^n g_j
            % and x_j = x_2 + sum_{t=3}^j g_t for j=3..n.
            % Finally m_j = exp(x_j) for j=2..n-1, and m_n = M.
            %
            % Inputs
            %   z       (n,1) positive weights
            %   c       (n,1) target moments for k=1..n
            %   M       (1,1) fixed value for m(n), must be > 0
            %   m0      (n,1) initial guess, must satisfy m0(1)=0, m0(n)=M, strictly increasing
            %   options (optional) optimoptions for lsqnonlin
            %
            % Outputs
            %   m       (n,1) solution with fixed endpoints
            %   u       (n-2,1) unconstrained parameters [u3..un]
            %   exitflag, output : from lsqnonlin
            %
            % Requires Optimization Toolbox (lsqnonlin).

            arguments
                z_int (:,1) double {mustBeFinite, mustBePositive}
                Lz (1,1) double {mustBeFinite, mustBePositive}
                rho_total (:,:,:) double {mustBeFinite}
                rho_nm0 (:,1) double {mustBeFinite, mustBeNonnegative}
                options.solver (1,1) string {mustBeMember(options.solver, ["fsolve","fminsearch"])} = "fminsearch"
            end

            n = numel(z_int);
            assert(numel(rho_nm0)==n, 'm0 must have length n.');
            assert(all(diff(rho_nm0)<0), 'm0 must be strictly decreasing.');

            z = flip(z_int/Lz);

            rho0 = rho_nm0(end);
            rhoD = rho_nm0(1);

            rho_moment = WVNoMotionProfileOperation.moments_from_rho_tot(rho_total, rho0, rhoD, z_int, Lz);
            % Initial u from m0 (log-gaps between nodes 2..n)
            u0 = WVNoMotionProfileOperation.rho_to_u(rho_nm0, rho0, rhoD);

            switch options.solver
                case "fsolve"
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

        function rho_moment = moments_from_rho_tot(rho_total, rho0, rhoD, z_int, Lz)
            int_vol_avg = @(integrand) sum(mean(mean(shiftdim(z_int,-2).*integrand,1),2),3)/Lz;

            rho_moment = zeros(size(rho_total,3),1);
            for i=1:size(rho_total,3)
                rho_moment(i) = int_vol_avg(((rho_total-rho0)/(rhoD - rho0)).^i);
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