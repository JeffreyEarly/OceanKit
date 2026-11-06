classdef IMBackground
    % Define the shared physical column used to create every mode family.
    %
    % A background derives density and buoyancy frequency from one supplied
    % profile in physical depth. The density at the upper endpoint
    % $$\rho(z_{\mathrm{surface}})=\rho_0$$ sets the Boussinesq reference
    % density and the integration constant.
    % Numerical coordinates and mode counts belong to the solver and mode
    % factories. Public state is read-only after construction.
    %
    % Source-specific factories also record analytical profile parameters.
    % The exponential profile is anchored at the upper domain endpoint:
    % $$N^2(z)=N_0^2\exp(2(z-z_{\mathrm{surface}})/b)$$.
    %
    % ```matlab
    % background = IMBackground(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
    % exactBackground = IMBackground.constantStratification(N0=1e-2,zDomain=[-1000 0],surfaceCondition="freeSurface");
    % ```
    %
    % - Topic: Create a mode background
    % - Topic: Inspect the physical column
    % - Topic: Create analytical solutions
    % - Declaration: classdef IMBackground

    properties (SetAccess = private)
        % Buoyancy-frequency-squared profile in radians squared per second squared.
        %
        % The function preserves the shape of its physical-depth input.
        %
        % - Topic: Inspect the physical column
        N2

        % Surface density and Boussinesq reference density in kg per cubic meter.
        % - Topic: Inspect the physical column
        rho0

        % Strictly increasing physical domain `[bottom surface]` in meters.
        %
        % - Topic: Inspect the physical column
        zDomain

        % Nonnegative planetary rotation rate in radians per second.
        %
        % - Topic: Inspect the physical column
        rotationRate

        % Latitude in degrees north, between -90 and 90.
        %
        % - Topic: Inspect the physical column
        latitude

        % Positive gravitational acceleration in meters per second squared.
        %
        % - Topic: Inspect the physical column
        g

        % Physical surface convention, `"freeSurface"` or `"rigidLid"`.
        %
        % - Topic: Inspect the physical column
        surfaceCondition

        % Source profile, `"arbitrary"`, `"constant"`, or `"exponential"`.
        %
        % - Topic: Inspect the physical column
        profileKind

        % Analytical source parameters recorded by a source-specific factory.
        %
        % Constant profiles record `N0`; exponential profiles record `N0`
        % and `b`, with `N0` defined at `zDomain(2)`.
        %
        % - Topic: Inspect the physical column
        profileParameters
    end

    properties (Access = private)
        % Supplies the density evaluator, preparing its primitive on first use.
        densityEvaluatorFactory_
    end

    properties (Dependent, SetAccess = private)
        % Shape-preserving density function in kg per cubic meter.
        % - Topic: Inspect the physical column
        rho

        % Signed Coriolis frequency $$f_0=2\Omega\sin(\mathrm{latitude})$$.
        %
        % - Topic: Inspect the physical column
        f0

        % Positive column depth $$D=z_{\mathrm{surface}}-z_{\mathrm{bottom}}$$.
        %
        % - Topic: Inspect the physical column
        D
    end

    methods
        function self = IMBackground(options)
            % Create a background from a declared stratification profile.
            %
            % Supply exactly one of `rho` or `N2`. A density source gives
            % $$N^2=-(g/\rho_0)\,\partial_z\rho$$ by Chebfun differentiation.
            % A stratification source is evaluated directly and supplies the
            % density through a surface-anchored Chebfun antiderivative.
            % Both conversions are independent of the solver coordinate and grid.
            % The supplied function must describe a fixed profile throughout
            % the lifetime of the background.
            %
            % - Topic: Create a mode background
            % - Declaration: background = IMBackground(options)
            % - Parameter options.N2: shape-preserving squared buoyancy frequency as a function of depth
            % - Parameter options.rho: shape-preserving density as a function of depth
            % - Parameter options.rho0: surface density, inferred from rho input or defaulting to 1025 for N2 input
            % - Parameter options.zDomain: strictly increasing `[bottom surface]` domain
            % - Parameter options.rotationRate: planetary rotation rate
            % - Parameter options.latitude: latitude in degrees north
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceCondition: required physical surface convention, freeSurface or rigidLid
            % - Returns self: read-only physical background
            arguments
                options.N2 (1,1) function_handle
                options.rho (1,1) function_handle
                options.rho0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.rotationRate (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.latitude,-90), mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceCondition {mustBeTextScalar}
            end
            if ~isfield(options,"surfaceCondition")
                error("IMBackground:MissingSurfaceCondition","Supply surfaceCondition=""freeSurface"" or surfaceCondition=""rigidLid"" explicitly.");
            end
            if options.zDomain(2) <= options.zDomain(1)
                error("IMBackground:InvalidDomain","zDomain must be strictly increasing, with the bottom before the surface.");
            end
            surfaceCondition = string(options.surfaceCondition);
            if ~ismember(surfaceCondition,["freeSurface", "rigidLid"])
                error("IMBackground:InvalidSurfaceCondition","surfaceCondition must be ""freeSurface"" or ""rigidLid"".");
            end
            hasN2 = isfield(options,"N2");
            hasDensity = isfield(options,"rho");
            if hasN2 == hasDensity
                error("IMBackground:InvalidProfileSource","Supply exactly one of rho or N2.");
            end
            if hasDensity
                zSurface = options.zDomain(2);
                rho = options.rho;
                rho0 = rho(zSurface);
                if ~isnumeric(rho0) || ~isscalar(rho0) || ~isreal(rho0) || ~isfinite(rho0) || rho0 <= 0
                    error("IMBackground:InvalidSurfaceDensity","rho must return one real, finite, positive density at the upper domain endpoint.");
                end
                rho0 = double(rho0);
                if isfield(options,"rho0") && options.rho0 ~= rho0
                    error("IMBackground:ConflictingSurfaceDensity","rho0 must equal rho at the upper domain endpoint. Omit rho0 to infer it from rho.");
                end
                densityProfile = chebfun(rho,options.zDomain,'splitting','on');
                representedSurface = densityProfile(zSurface);
                gradient = diff(densityProfile);
                buoyancyScale = -options.g/rho0;
                self.N2 = @(z) buoyancyScale*evaluateProfile(gradient,z);
                densityEvaluator = @(z) rho0+(evaluateProfile(densityProfile,z)-representedSurface);
                self.densityEvaluatorFactory_ = @() densityEvaluator;
            else
                rho0 = 1025;
                if isfield(options,"rho0"), rho0 = options.rho0; end
                self.N2 = options.N2;
                self.densityEvaluatorFactory_ = densityFromStratification(options.N2,options.zDomain,options.g,rho0);
            end
            self.rho0 = rho0;
            self.zDomain = options.zDomain;
            self.rotationRate = options.rotationRate;
            self.latitude = options.latitude;
            self.g = options.g;
            self.surfaceCondition = surfaceCondition;
            self.profileKind = "arbitrary";
            self.profileParameters = struct();
        end

        function value = get.rho(self)
            value = self.densityEvaluatorFactory_();
        end

        function value = get.f0(self)
            value = 2*self.rotationRate*sind(self.latitude);
        end

        function value = get.D(self)
            value = diff(self.zDomain);
        end

        function solution = analyticalSolution(self)
            % Create the exact solution supported by this background source.
            %
            % Parameters are taken from this background rather than a
            % separately supplied analytical stratification. The current
            % exponential analytical backend requires the surface at zero;
            % translated exponential backgrounds remain available to
            % numerical solvers.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = analyticalSolution(background)
            % - Returns solution: matching analytical solution family
            arguments
                self (1,1) IMBackground
            end
            switch self.profileKind
                case "constant"
                    solution = IMConstantStratificationSolution(N0=self.profileParameters.N0,zDomain=self.zDomain,rotationRate=self.rotationRate,latitude=self.latitude,g=self.g);
                case "exponential"
                    if self.zDomain(2) ~= 0
                        error("IMBackground:UnsupportedAnalyticalDomain","The exponential analytical solver requires zDomain(2)=0. Use a numerical solver for a translated domain.");
                    end
                    solution = IMExponentialStratificationSolution(N0=self.profileParameters.N0,b=self.profileParameters.b,zDomain=self.zDomain,rotationRate=self.rotationRate,latitude=self.latitude,g=self.g);
                otherwise
                    error("IMBackground:UnsupportedAnalyticalProfile","An arbitrary N2 profile has no declared analytical solution. Use a source-specific background factory or a numerical solver.");
            end
        end
    end

    methods (Static)
        function self = constantStratification(options)
            % Create a constant background with one authoritative $$N_0$$.
            %
            % - Topic: Create a mode background
            % - Declaration: background = IMBackground.constantStratification(options)
            % - Parameter options.N0: positive constant buoyancy frequency
            % - Parameter options.rho0: surface density, default 1025 kg per cubic meter
            % - Parameter options.zDomain: strictly increasing physical domain
            % - Parameter options.rotationRate: planetary rotation rate
            % - Parameter options.latitude: latitude in degrees north
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceCondition: required physical surface convention, freeSurface or rigidLid
            % - Returns self: constant-stratification background
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 5.2e-3
                options.rho0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.rotationRate (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.latitude,-90), mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceCondition {mustBeTextScalar}
            end
            if ~isfield(options,"surfaceCondition")
                error("IMBackground:MissingSurfaceCondition","Supply surfaceCondition=""freeSurface"" or surfaceCondition=""rigidLid"" explicitly.");
            end
            N0 = options.N0;
            options = rmfield(options,"N0");
            options.N2 = @(z) N0^2*ones(size(z));
            args = namedargs2cell(options);
            self = IMBackground(args{:});
            self.profileKind = "constant";
            self.profileParameters = struct("N0",N0);
        end

        function self = exponentialStratification(options)
            % Create an exponential background anchored at the surface.
            %
            % `N0` is the buoyancy frequency at `zDomain(2)`, so the profile
            % is unchanged by translating the physical coordinate origin.
            %
            % - Topic: Create a mode background
            % - Declaration: background = IMBackground.exponentialStratification(options)
            % - Parameter options.N0: positive surface buoyancy frequency
            % - Parameter options.rho0: surface density, default 1025 kg per cubic meter
            % - Parameter options.b: positive buoyancy-frequency e-folding depth
            % - Parameter options.zDomain: strictly increasing physical domain
            % - Parameter options.rotationRate: planetary rotation rate
            % - Parameter options.latitude: latitude in degrees north
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceCondition: required physical surface convention, freeSurface or rigidLid
            % - Returns self: exponential-stratification background
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 5.2e-3
                options.rho0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.b (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1300
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.rotationRate (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.latitude,-90), mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceCondition {mustBeTextScalar}
            end
            if ~isfield(options,"surfaceCondition")
                error("IMBackground:MissingSurfaceCondition","Supply surfaceCondition=""freeSurface"" or surfaceCondition=""rigidLid"" explicitly.");
            end
            N0 = options.N0;
            b = options.b;
            zSurface = options.zDomain(2);
            options = rmfield(options,["N0", "b"]);
            options.N2 = @(z) N0^2*exp(2*(z-zSurface)/b);
            args = namedargs2cell(options);
            self = IMBackground(args{:});
            self.profileKind = "exponential";
            self.profileParameters = struct("N0",N0,"b",b);
        end
    end
end

function factory = densityFromStratification(N2,zDomain,g,rho0)
% One primitive serves all density queries for the fixed source profile.
densityEvaluator = [];
factory = @getDensityEvaluator;

    function value = getDensityEvaluator()
        if isempty(densityEvaluator)
            primitive = cumsum(chebfun(N2,zDomain,'splitting','on'));
            surfaceIntegral = primitive(zDomain(2));
            densityEvaluator = @(z) rho0-(rho0/g)*(evaluateProfile(primitive,z)-surfaceIntegral);
        end
        value = densityEvaluator;
    end
end

function values = evaluateProfile(profile,z)
values = reshape(feval(profile,z(:)),size(z));
end
