classdef IMGeostrophicZeroAPVModes
    % Describe canonical geostrophic zero-APV modes at fixed wavenumber.
    %
    % `IMGeostrophicZeroAPVModes` stores the coefficient-independent
    % boundary-value problem
    %
    % $$
    % F_z=-\frac{N^2}{g}G,\qquad G_z=-\frac{gk^2}{f_0^2}F.
    % $$
    %
    % This pair is regular for finite nonnegative N2, including neutral
    % intervals. The canonical response vector contains the requested subset of
    %
    % $$
    % B_{\mathrm s}[F]=
    % \begin{cases}
    % G(0)-F(0), & \texttt{freeSurface},\\
    % G(0), & \texttt{rigidLid},
    % \end{cases}
    % \qquad
    % B_{\mathrm d}[F]=G(z_b).
    % $$
    %
    % These responses are proportional to the physical endpoint anomaly
    % functionals:
    %
    % $$
    % \eta_0[F]=\frac{f_0}{g}B_{\mathrm s}[F],
    % \qquad
    % \eta_d[F]=\frac{f_0}{g}B_{\mathrm d}[F].
    % $$
    % For positive endpoint stratification, the same functionals can be written
    %
    % $$
    % \eta_0[F]=
    % \begin{cases}
    % -\dfrac{f_0}{N^2(0)}\dfrac{\partial F}{\partial z}(0)-\dfrac{f_0}{g}F(0), & \texttt{freeSurface},\\
    % -\dfrac{f_0}{N^2(0)}\dfrac{\partial F}{\partial z}(0), & \texttt{rigidLid},
    % \end{cases}
    % \qquad
    % \eta_d[F]= -\frac{f_0}{N^2(z_b)}\frac{\partial F}{\partial z}(z_b).
    % $$
    %
    % Solvers return one mode with unit response at each requested endpoint
    % and zero response at the other endpoint. Generalized-energy
    % coefficients are supplied later to basis orthogonalization methods; they do
    % not change this problem or decide which endpoint coordinates exist.
    % For two endpoints,
    %
    % $$
    % \mathbf B[F_0^{\mathrm{sur}}]=(1,0)^T,
    % \qquad
    % \mathbf B[F_0^{\mathrm{bot}}]=(0,1)^T.
    % $$
    %
    % ```matlab
    % problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-4000 0],rotationRate=7.2921e-5,latitude=33,k=1e-4,endpoints=["surface","bottom"],surfaceBoundary="freeSurface");
    % solver = IMSolverSpectral(nEVP=128);
    % boundaryModes = solver.solveGeostrophicZeroAPVModes(problem);
    % orthogonal = boundaryModes.orthogonalizeGeneralizedEnergy(g0=-0.035,gd=0);
    % ```
    %
    % - Topic: Create geostrophic zero-APV problems
    % - Topic: Summarize geostrophic zero-APV problems
    % - Topic: Inspect geostrophic zero-APV problems
    % - Declaration: classdef IMGeostrophicZeroAPVModes

    properties (SetAccess = private)
        % Authoritative source background, when constructed from one.
        % - Topic: Inspect geostrophic zero-APV problems
        background = []

        % Buoyancy frequency squared function.
        %
        % `N2` has signature `values = N2(z)` and defines the continuous
        % stratification used by numerical solvers.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        N2

        % Physical vertical domain.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        zDomain

        % Planetary rotation rate in radians per second.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        rotationRate

        % Latitude in degrees north.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        latitude

        % Gravitational acceleration.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        g

        % Horizontal wavenumbers.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        k

        % Canonically ordered endpoint coordinates.
        %
        % `endpoints` is a possibly empty subset of `["surface","bottom"]`.
        % Its order defines the second dimension of solved mode arrays.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        endpoints

        % Surface endpoint convention.
        %
        % `surfaceBoundary` is `"freeSurface"` or `"rigidLid"`.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        surfaceBoundary

        % Additional creation metadata.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        metadata
    end

    properties (Dependent, SetAccess = private)
        % Coriolis parameter in radians per second.
        %
        % $$f_0=2\Omega\sin\varphi$$ is derived from `rotationRate` and
        % `latitude` and must be nonzero for this problem.
        %
        % - Topic: Inspect geostrophic zero-APV problems
        f0
    end

    methods
        function self = IMGeostrophicZeroAPVModes(options)
            % Create a canonical geostrophic zero-APV problem.
            %
            % - Topic: Create geostrophic zero-APV problems
            % - Declaration: problem = IMGeostrophicZeroAPVModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.rotationRate: nonnegative planetary rotation rate in radians per second
            % - Parameter options.latitude: latitude in degrees north, between -90 and 90
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.endpoints: requested endpoint coordinates
            % - Parameter options.surfaceBoundary: surface endpoint convention
            % - Parameter options.metadata: additional metadata
            % - Returns problem: geostrophic zero-APV problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.rotationRate (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.latitude,-90), mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.endpoints {mustBeText} = ["surface", "bottom"]
                options.surfaceBoundary {mustBeTextScalar, mustBeMember(options.surfaceBoundary, ["freeSurface", "rigidLid"])} = "freeSurface"
                options.metadata struct = struct()
            end

            if 2*options.rotationRate*sind(options.latitude) == 0
                error("IMGeostrophicZeroAPVModes:InvalidCoriolis", "Rotation rate and latitude must give a nonzero f0 for geostrophic zero-APV modes.");
            end
            if isempty(options.k)
                error("IMGeostrophicZeroAPVModes:EmptyWavenumber", "k must contain at least one positive horizontal wavenumber.");
            end
            endpoints = reshape(string(options.endpoints), 1, []);
            canonicalEndpoints = ["surface", "bottom"];
            if any(~ismember(endpoints, canonicalEndpoints))
                error("IMGeostrophicZeroAPVModes:InvalidEndpoint", "endpoints must contain only ""surface"" and ""bottom"".");
            end
            if numel(unique(endpoints)) ~= numel(endpoints)
                error("IMGeostrophicZeroAPVModes:DuplicateEndpoint", "endpoints must not contain duplicate values.");
            end

            self.N2 = options.N2;
            self.zDomain = sort(options.zDomain);
            self.rotationRate = options.rotationRate;
            self.latitude = options.latitude;
            self.g = options.g;
            self.k = reshape(options.k, 1, []);
            self.endpoints = canonicalEndpoints(ismember(canonicalEndpoints, endpoints));
            self.surfaceBoundary = string(options.surfaceBoundary);
            self.metadata = options.metadata;
        end

        function value = get.f0(self)
            value = 2*self.rotationRate*sind(self.latitude);
        end

        function summarize(self)
            % Print a readable problem summary.
            %
            % - Topic: Summarize geostrophic zero-APV problems
            % - Declaration: summarize(problem)
            arguments
                self IMGeostrophicZeroAPVModes
            end

            fprintf("%s\n", class(self));
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  rotationRate: %g rad/s\n", self.rotationRate);
            fprintf("  latitude: %g degrees north\n", self.latitude);
            fprintf("  f0: %g\n", self.f0);
            fprintf("  g: %g\n", self.g);
            fprintf("  endpoints: %s\n", join(self.endpoints, ", "));
            fprintf("  surfaceBoundary: %s\n", self.surfaceBoundary);
            fprintf("  nWavenumbers: %d\n", numel(self.k));
            fprintf("  modesPerWavenumber: %d\n", self.modesPerWavenumber());
        end

        function count = modesPerWavenumber(self)
            % Return the number of canonical modes for each wavenumber.
            %
            % - Topic: Inspect geostrophic zero-APV problems
            % - Declaration: count = modesPerWavenumber(problem)
            % - Returns count: number of requested endpoint coordinates
            arguments
                self IMGeostrophicZeroAPVModes
            end

            count = numel(self.endpoints);
        end
    end

    methods (Static)
        function problem = atWavenumber(varargin)
            % Create canonical endpoint responses from a background or physical inputs.
            %
            % `atWavenumber(background,k,endpoints="surface")` takes all physical
            % state from the background. The explicit name-value form remains
            % available. Mixing a background with explicit physical inputs is
            % rejected so there is one authoritative physical column.
            %
            % - Topic: Create geostrophic zero-APV problems
            % - Declaration: problem = IMGeostrophicZeroAPVModes.atWavenumber(background,k,options)
            % - Parameter background: optional authoritative IMBackground
            % - Parameter k: positive horizontal wavenumbers
            % - Returns problem: canonical geostrophic zero-APV problem
            if ~isempty(varargin) && isa(varargin{1},"IMBackground")
                background = varargin{1};
                if numel(varargin) >= 2 && isnumeric(varargin{2})
                    k = varargin{2};
                    supplied = varargin(3:end);
                else
                    supplied = varargin(2:end);
                    if mod(numel(supplied),2) ~= 0
                        error("IMGeostrophicZeroAPVModes:InvalidOptions","Supply name-value options after the background.");
                    end
                    kIndex = find(string(supplied(1:2:end)) == "k");
                    if numel(kIndex) ~= 1
                        error("IMGeostrophicZeroAPVModes:MissingWavenumber","Supply one horizontal wavenumber input after the background, positionally or with k=.");
                    end
                    k = supplied{2*kIndex};
                    supplied(2*kIndex-1:2*kIndex) = [];
                end
                if mod(numel(supplied),2) ~= 0
                    error("IMGeostrophicZeroAPVModes:InvalidOptions","Supply name-value options after background and k.");
                end
                names = string(supplied(1:2:end));
                if any(ismember(names,["N2","zDomain","rotationRate","latitude","g","surfaceBoundary","k"]))
                    error("IMGeostrophicZeroAPVModes:ConflictingBackground","Physical inputs come from the background. Supply only endpoint or metadata options alongside it.");
                end
                problem = IMGeostrophicZeroAPVModes(supplied{:},N2=background.N2,zDomain=background.zDomain,rotationRate=background.rotationRate,latitude=background.latitude,g=background.g,k=k,surfaceBoundary=background.surfaceCondition);
                problem.background = background;
            else
                problem = IMGeostrophicZeroAPVModes(varargin{:});
            end
        end
    end
end
