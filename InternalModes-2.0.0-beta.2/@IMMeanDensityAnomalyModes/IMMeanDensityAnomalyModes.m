classdef IMMeanDensityAnomalyModes < IMInternalModes
    % Describe generalized-energy mean-density-anomaly modes.
    %
    % `IMMeanDensityAnomalyModes` describes the `G`-form problem
    %
    % $$
    % -G_j''(z)=\frac{N^2(z)}{g h_j}G_j(z).
    % $$
    %
    % A finite surface acceleration imposes
    % $$g h_jG_j'(z_s)=g_0G_j(z_s)$$, and a finite bottom
    % acceleration imposes
    % $$g h_jG_j'(z_b)=-g_dG_j(z_b).$$ Positive infinity makes the
    % corresponding endpoint inactive and imposes Dirichlet data. Finite
    % zero therefore gives the Neumann limit.
    %
    % The solved `G` modes use the signed generalized-energy product
    %
    % $$
    % \langle G_i,G_j\rangle=
    % \frac{1}{g}\int_{z_b}^{z_s}N^2G_iG_j\,dz
    % +\frac{g_0}{g}G_i(z_s)G_j(z_s)
    % +\frac{g_d}{g}G_i(z_b)G_j(z_b),
    % $$
    %
    % with inactive endpoint terms omitted. The companion pressure mode is
    % the surface-referenced diagnostic
    %
    % $$
    % F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz',
    % \qquad F_j(z_s)=0.
    % $$
    %
    % Construct this descriptor through
    % `IMInternalModes.meanDensityAnomalyModes`.
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Inspect internal-mode configuration
    % - Declaration: classdef IMMeanDensityAnomalyModes < IMInternalModes

    properties (SetAccess = private)
        % Surface generalized-energy acceleration.
        %
        % Finite values, including zero, are active. Positive infinity
        % selects an inactive Dirichlet endpoint.
        %
        % - Topic: Inspect internal-mode configuration
        g0

        % Bottom generalized-energy acceleration.
        %
        % Finite values, including zero, are active. Positive infinity
        % selects an inactive Dirichlet endpoint.
        %
        % - Topic: Inspect internal-mode configuration
        gd

        % Active endpoint identities in canonical surface-bottom order.
        %
        % - Topic: Inspect internal-mode configuration
        activeEndpoints
    end

    methods
        function self = IMMeanDensityAnomalyModes(options)
            % Create a generalized-energy mean-density-anomaly EVP.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMMeanDensityAnomalyModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.g0: signed finite surface acceleration, zero, or positive infinity
            % - Parameter options.gd: signed finite bottom acceleration, zero, or positive infinity
            % - Returns evp: mean-density-anomaly EVP descriptor
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
            end

            if isnan(options.g0) || options.g0 == -Inf
                error("IMMeanDensityAnomalyModes:InvalidSurfaceAcceleration", "g0 must be signed finite, zero, or positive Inf.");
            end
            if isnan(options.gd) || options.gd == -Inf
                error("IMMeanDensityAnomalyModes:InvalidBottomAcceleration", "gd must be signed finite, zero, or positive Inf.");
            end

            activeEndpoints = strings(1,0);
            if isfinite(options.g0)
                surfaceCondition = IMBoundaryCondition(a=0,b=1,c=options.g0/options.g,d=0);
                activeEndpoints(end+1) = "surface";
            else
                surfaceCondition = IMBoundaryCondition.dirichlet();
            end
            if isfinite(options.gd)
                bottomCondition = IMBoundaryCondition(a=0,b=1,c=-options.gd/options.g,d=0);
                activeEndpoints(end+1) = "bottom";
            else
                bottomCondition = IMBoundaryCondition.dirichlet();
            end

            parameters = struct("g0",options.g0,"gd",options.gd,"activeEndpoints",activeEndpoints);
            self@IMInternalModes(name="meanDensityAnomalyModes",formulation="G",modeFamily="meanDensityAnomaly", ...
                N2=options.N2,zDomain=options.zDomain,p=@(z,~) ones(size(z)),q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g,g=options.g,surfaceBoundary=surfaceCondition, ...
                bottomBoundary=bottomCondition,parameters=parameters);
            self.g0 = options.g0;
            self.gd = options.gd;
            self.activeEndpoints = activeEndpoints;
        end

        function spec = innerProduct(self, variable)
            % Return the mean-density-anomaly `F` or `G` inner product.
            %
            % `G` has the signed generalized-energy product defined by the
            % canonical EVP. `F` is a diagnostic pressure shape and does
            % not define a direct projection metric.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: `"F"` or `"G"`
            % - Returns spec: continuous inner-product specification
            arguments
                self IMMeanDensityAnomalyModes
                variable {mustBeTextScalar, mustBeMember(variable,["F","G"])} = "G"
            end

            spec = innerProduct@IMInternalModes(self,variable);
            if string(variable) == "F"
                spec.hasInnerProduct = false;
                spec.reason = "F is a surface-referenced diagnostic pressure mode; mean-density-anomaly coefficient projection is defined through G.";
            end
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, modeSelectionDiagnostics)
            % Create a mean-density-anomaly basis set.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,modeSelectionDiagnostics)
            % - Returns basisSet: aligned mean-density-anomaly basis set
            % - Developer: true
            arguments
                self IMMeanDensityAnomalyModes
                solver IMSolver
                nativeModes (:,:) double
                eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                modeNumber (1,:) double {mustBeInteger}
                modeSelectionDiagnostics struct
            end

            metadata = struct("g0",self.g0,"gd",self.gd,"activeEndpoints",self.activeEndpoints, ...
                "diagnosticFConvention","surfaceReferencedN2GIntegral");
            basisSet = IMMeanDensityAnomalyModesBasis(solver=solver,evp=self,nativeModes=nativeModes, ...
                eigenvalues=eigenvalues,modeNumber=modeNumber,modeSelectionDiagnostics=modeSelectionDiagnostics,metadata=metadata);
        end
    end

end
