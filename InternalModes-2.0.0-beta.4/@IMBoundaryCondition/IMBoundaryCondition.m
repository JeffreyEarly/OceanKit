classdef IMBoundaryCondition
    % Store one scalar canonical boundary condition.
    %
    % `IMBoundaryCondition` represents
    %
    % $$
    % -\left[
    % a\,u(z_\ell)
    % -b\,p(z_\ell)\frac{\partial u}{\partial z}(z_\ell)
    % \right]
    % =
    % \lambda\left[
    % c\,u(z_\ell)
    % -d\,p(z_\ell)\frac{\partial u}{\partial z}(z_\ell)
    % \right].
    % $$
    %
    % The coefficients define the boundary condition. The EVP supplies the
    % endpoint location when it assembles matrix rows or computes endpoint weight
    % coefficients. The stored properties `a`, `b`, `c`, and `d` define
    % the boundary condition; `determinant(location)` derives the signed
    % determinant $$D_i$$; and `endpointWeightCoefficient(location)`
    % derives the scalar $$D_i^{-1}$$ used by
    % `IMEigenvalueProblem.endpointWeights`.
    %
    % At endpoint $$z_\ell$$, define
    %
    % $$
    % P_\ell=p(z_\ell)\frac{\partial u}{\partial z}(z_\ell).
    % $$
    %
    % The signed determinant is
    %
    % $$
    % D_\ell=\sigma_\ell(ad-bc),\qquad
    % \sigma_\mathrm{bottom}=+1,\quad \sigma_\mathrm{surface}=-1.
    % $$
    %
    % These endpoint weights appear in the full solved-variable inner
    % product
    %
    % $$
    % \langle u_i,u_j\rangle_\mu
    % =
    % \int_{z_b}^{z_s}r(z)u_i(z)u_j(z)\,dz
    % +
    % \sum_{\ell\in S}D_\ell^{-1}L_\ell[u_i]L_\ell[u_j],
    % $$
    %
    % where $$S$$ is the active, nondegenerate endpoint set and
    %
    % $$
    % L_\ell[u]=c\,u(z_\ell)-d\,p(z_\ell)\frac{\partial u}{\partial z}(z_\ell).
    % $$
    %
    % The table's "Endpoint norm contribution" column is one summand in
    % this endpoint sum.
    %
    % Only active boundary conditions, where `(c,d) ~= (0,0)`, produce
    % endpoint norm weights:
    %
    % | Boundary condition | Endpoint equation | Endpoint norm contribution |
    % | --- | --- | --- |
    % | Dirichlet | $$u_\ell=0$$ | none |
    % | Neumann | $$P_\ell=0$$ | none |
    % | Robin | $$a u_\ell-bP_\ell=0$$ | none; `robinEnergyCoefficient` is an energy diagnostic |
    % | Active value condition | eigenvalue side uses $$c u_\ell$$ | $$D_\ell^{-1}(c u_\ell)^2$$ |
    % | Active flux condition | eigenvalue side uses $$-dP_\ell$$ | $$D_\ell^{-1}(-dP_\ell)^2$$ |
    % | General active condition | $$-[a u_\ell-bP_\ell]=\lambda[c u_\ell-dP_\ell]$$ | $$D_\ell^{-1}(c u_\ell-dP_\ell)^2$$ |
    % | Degenerate active condition | `(c,d) ~= (0,0)` and $$D_\ell=0$$ | unavailable |
    %
    % For example, `IMBoundaryCondition(a=0,b=1,c=1,d=0)` has
    % $$D_\mathrm{surface}=+1$$ and $$D_\mathrm{bottom}=-1$$, so the
    % same boundary condition contributes $$+u_s^2$$ at the surface and
    % $$-u_b^2$$ at the bottom.
    %
    % - Topic: Create boundary conditions
    % - Topic: Inspect boundary conditions
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundaryCondition

    properties (SetAccess = private)
        % Coefficient multiplying endpoint value on the left.
        %
        % - Topic: Inspect boundary conditions
        a = 1

        % Coefficient multiplying endpoint flux on the left.
        %
        % - Topic: Inspect boundary conditions
        b = 0

        % Coefficient multiplying endpoint value on the eigenvalue side.
        %
        % - Topic: Inspect boundary conditions
        c = 0

        % Coefficient multiplying endpoint flux on the eigenvalue side.
        %
        % - Topic: Inspect boundary conditions
        d = 0
    end

    methods
        function self = IMBoundaryCondition(options)
            % Create a scalar boundary condition.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundaryCondition(options)
            % - Parameter options.a: value coefficient on the left
            % - Parameter options.b: flux coefficient on the left
            % - Parameter options.c: value coefficient on the eigenvalue side
            % - Parameter options.d: flux coefficient on the eigenvalue side
            % - Returns boundary: boundary condition
            arguments
                options.a (1,1) double {mustBeReal, mustBeFinite} = 1
                options.b (1,1) double {mustBeReal, mustBeFinite} = 0
                options.c (1,1) double {mustBeReal, mustBeFinite} = 0
                options.d (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            coefficients = [options.a options.b options.c options.d];
            if all(coefficients == 0)
                error("IMBoundaryCondition:DegenerateCondition", ...
                    "At least one boundary coefficient must be nonzero.");
            end

            self.a = options.a;
            self.b = options.b;
            self.c = options.c;
            self.d = options.d;
        end

        function tf = isEigenvalueDependent(self, tolerance)
            % Return true when the eigenvalue side is active.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: tf = isEigenvalueDependent(boundary,tolerance)
            % - Parameter tolerance: scalar activity tolerance
            % - Returns tf: true for active endpoint metric terms
            arguments
                self IMBoundaryCondition
                tolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = self.activityTolerance()
            end

            tf = abs(self.c) > tolerance || abs(self.d) > tolerance;
        end

        function value = determinant(self, location)
            % Return the signed endpoint determinant.
            %
            % Yassin's endpoint indexing uses `z_1` for the bottom endpoint
            % and `z_2` for the surface endpoint, so the sign
            % `(-1)^(i+1)` is positive at the bottom and negative at the
            % surface.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: value = determinant(boundary,location)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Returns value: signed determinant
            arguments
                self IMBoundaryCondition
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
            end

            value = IMBoundaryCondition.endpointSign(location)*(self.a*self.d - self.b*self.c);
        end

        function value = endpointWeightCoefficient(self, location)
            % Return the endpoint weight coefficient.
            %
            % For the active boundary condition
            %
            % $$
            % -\left[
            % a_i u(z_i)
            % -b_i p(z_i)\frac{\partial u}{\partial z}(z_i)
            % \right]
            % =
            % \lambda\left[
            % c_i u(z_i)
            % -d_i p(z_i)\frac{\partial u}{\partial z}(z_i)
            % \right].
            % $$
            %
            % the stored properties `a`, `b`, `c`, and `d` define
            %
            % $$
            % D_i=(-1)^{i+1}(a_i d_i-b_i c_i),
            % $$
            %
            % where Yassin's endpoint indexing makes the sign positive at
            % the bottom and negative at the surface. This method returns
            % $$D_i^{-1}$$, the scalar coefficient multiplying
            %
            % $$
            % \left(
            % c_i u(z_i)
            % -d_i p(z_i)\frac{\partial u}{\partial z}(z_i)
            % \right)^2
            % $$
            %
            % in the endpoint part of the norm. `IMEigenvalueProblem`
            % copies this scalar into the `coefficient` field of each
            % `endpointWeights` struct.
            %
            % This method returns only the scalar coefficient
            % $$D_i^{-1}$$. The full endpoint norm term is assembled by
            % `IMEigenvalueProblem.endpointWeights` from this coefficient
            % and the boundary properties `c` and `d`.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: value = endpointWeightCoefficient(boundary,location)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Returns value: coefficient multiplying `(c*u-d*p*u_z)^2`
            arguments
                self IMBoundaryCondition
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
            end

            D = self.determinant(location);
            tolerance = 100*eps*max(1,abs(D));
            if abs(D) <= tolerance
                value = NaN;
                return;
            end
            value = 1/D;
        end

        function beta = robinEnergyCoefficient(self, location)
            % Return the ordinary Robin endpoint quadratic coefficient.
            %
            % For inactive boundary conditions, this is
            % `(-1)^(i+1)*a/b` with Yassin's `z_1` bottom and `z_2`
            % surface indexing.
            %
            % - Topic: Developer topics
            % - Declaration: beta = robinEnergyCoefficient(boundary,location)
            % - Developer: true
            arguments
                self IMBoundaryCondition
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
            end

            if self.b == 0
                beta = 0;
                return;
            end
            beta = IMBoundaryCondition.endpointSign(location)*self.a/self.b;
        end

        function H = endpointNumeratorMatrix(self, location)
            % Return the active-endpoint numerator matrix.
            %
            % The vector is `[u; p*u_z]`.
            %
            % - Topic: Developer topics
            % - Declaration: H = endpointNumeratorMatrix(boundary,location)
            % - Developer: true
            arguments
                self IMBoundaryCondition
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
            end

            D = self.determinant(location);
            H = -(1/D)*[self.a*self.c self.a*self.d; self.a*self.d self.b*self.d];
        end
    end

    methods (Access = private)
        function tolerance = activityTolerance(self)
            tolerance = 100*eps*max([1 abs(self.a) abs(self.b) abs(self.c) abs(self.d)]);
        end
    end

    methods (Static)
        function boundary = dirichlet()
            % Create `u=0`.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundaryCondition.dirichlet()
            boundary = IMBoundaryCondition(a=1, b=0, c=0, d=0);
        end

        function boundary = neumann()
            % Create `p*u_z=0`.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundaryCondition.neumann()
            boundary = IMBoundaryCondition(a=0, b=1, c=0, d=0);
        end

        function boundary = robin(a, b)
            % Create `a*u-b*p*u_z=0`.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundaryCondition.robin(a,b)
            % - Parameter a: endpoint value coefficient
            % - Parameter b: endpoint flux coefficient
            arguments
                a (1,1) double {mustBeReal, mustBeFinite}
                b (1,1) double {mustBeReal, mustBeFinite}
            end

            boundary = IMBoundaryCondition(a=a, b=b, c=0, d=0);
        end
    end

    methods (Static, Access = private)
        function value = endpointSign(location)
            arguments
                location {mustBeTextScalar, mustBeMember(location, ["surface", "bottom"])}
            end

            switch string(location)
                case "surface"
                    value = -1;
                case "bottom"
                    value = 1;
            end
        end
    end
end
