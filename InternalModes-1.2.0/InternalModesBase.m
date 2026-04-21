classdef (Abstract) InternalModesBase < handle
    % Define the shared contract for all internal-mode solvers.
    %
    % `InternalModesBase` is the abstract base class for the concrete
    % solver implementations. It stores the common physical parameters,
    % normalization choices, and boundary-condition state, and it defines
    % the shared public interface for requesting modes at fixed wavenumber
    % or fixed frequency.
    %
    % Following Section 2.3 of Early, Lelong, and Smith (2020), the
    % package centers on the vertical structure functions $$F_j(z)$$ and
    % $$G_j(z)$$ with equivalent depths $$h_j$$, connected by
    %
    % $$
    % (N^2 - \omega^2) G_j = -g \, \partial_z F_j,
    % \qquad
    % F_j = h_j \, \partial_z G_j.
    % $$
    %
    % Concrete subclasses solve either the fixed-$$K$$ or fixed-$$\omega$$
    % eigenvalue problems using spectral, finite-difference, analytical,
    % or WKB-based formulations, but they all expose the same shared state
    % documented here.
    %
    % - Topic: Inspect grids and stratification
    % - Topic: Configure normalization and boundaries
    % - Topic: Compute modes
    % - Topic: Developer topics
    % - Declaration: classdef (Abstract) InternalModesBase < handle

    
    properties (Access = public)
        % Toggle diagnostic messages from the active solver.
        %
        % - Topic: Developer topics
        % - Developer: true
        shouldShowDiagnostics = 0 % flag to show diagnostic information, default = 0
        
        % Planetary rotation rate in radians per second.
        %
        % - Topic: Inspect grids and stratification
        rotationRate % rotation rate of the planetary body
        % Latitude in degrees used to compute `f0`.
        %
        % - Topic: Inspect grids and stratification
        latitude % Latitude for which the modes are being computed.
        % Coriolis parameter at the selected latitude.
        %
        % - Topic: Inspect grids and stratification
        f0 % Coriolis parameter at the above latitude.
        % Total water-column depth $$D = zMax - zMin$$.
        %
        % - Topic: Inspect grids and stratification
        Lz % Depth of the ocean.
        % Reference surface density.
        %
        % - Topic: Inspect grids and stratification
        rho0 % Density at the surface of the ocean.
        
        % Optional cap on the number of modes returned.
        %
        % - Topic: Compute modes
        nModes % Used to limit the number of modes to be returned.
        
        % Two-element depth domain `[zMin zMax]`.
        %
        % - Topic: Inspect grids and stratification
        zDomain % [zMin zMax]
        % Flag indicating whether the concrete solver requires monotonic density.
        %
        % - Topic: Developer topics
        % - Developer: true
        requiresMonotonicDensity
        
        % Last fixed frequency used to build an adaptive grid, when applicable.
        %
        % - Topic: Developer topics
        % - Developer: true
        gridFrequency = [] % last requested frequency from the user---set to f0 if a wavenumber was last requested
        % Mode normalization convention.
        %
        % - Topic: Configure normalization and boundaries
        normalization = Normalization.kConstant % Normalization used for the modes. Either Normalization.(kConstant, omegaConstant, uMax, or wMax).
        % Upper boundary condition at the ocean surface.
        %
        % - Topic: Configure normalization and boundaries
        upperBoundary = UpperBoundary.rigidLid  % Surface boundary condition. Either UpperBoundary.rigidLid (default) or UpperBoundary.freeSurface.
        % Lower boundary condition at the ocean bottom.
        %
        % - Topic: Configure normalization and boundaries
        lowerBoundary = LowerBoundary.freeSlip  % Lower boundary condition. Either LowerBoundary.freeSlip (default) or LowerBoundary.none.
    end
    
    properties (SetObservable, AbortSet, Access = public)
        % Output depth grid on which public profiles and modes are returned.
        %
        % - Topic: Inspect grids and stratification
        z % Depth coordinate grid used for all output (same as zOut).
    end    

    properties (Dependent)
        % Minimum depth in `zDomain`.
        %
        % - Topic: Inspect grids and stratification
        zMin
        % Maximum depth in `zDomain`.
        %
        % - Topic: Inspect grids and stratification
        zMax
    end
    
    properties (Abstract)
        % Background density sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho  % Density on the z grid.
        % Buoyancy frequency squared sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        % First depth derivative of the background density sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_z % First derivative of density on the z grid.
        % Second depth derivative of the background density sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_zz % Second derivative of density on the z grid.
    end
    
    properties (Access = protected)
        % Gravitational acceleration.
        %
        % - Topic: Developer topics
        % - Developer: true
        g % 9.81 meters per second.
        % Dispersion-relation map from `h` and `k` to `omega`.
        %
        % - Topic: Developer topics
        % - Developer: true
        omegaFromK % function handle to compute omega(h,k)
        % Dispersion-relation map from `h` and `omega` to `k`.
        %
        % - Topic: Developer topics
        % - Developer: true
        kFromOmega % function handle to compute k(h,omega)
    end
    
    methods (Abstract)
        % Return vertical modes for a fixed horizontal wavenumber.
        %
        % - Topic: Compute modes
        % - Declaration: [F,G,h,omega,varargout] = modesAtWavenumber(self,k,varargin)
        % - Parameter self: concrete InternalModesBase subclass instance
        % - Parameter k: horizontal wavenumber
        % - Parameter varargin: additional requests supported by the concrete solver
        % - Returns F: horizontal-velocity mode matrix on `zOut`
        % - Returns G: vertical-velocity mode matrix on `zOut`
        % - Returns h: equivalent-depth row vector
        % - Returns omega: frequency row vector implied by `h` and `k`
        % - Returns varargout: additional outputs forwarded from the concrete solver
        [F,G,h,omega,varargout] = modesAtWavenumber(self, k, varargin ) % Return the normal modes and eigenvalue at a given wavenumber.
        
        % Return vertical modes for a fixed frequency.
        %
        % - Topic: Compute modes
        % - Declaration: [F,G,h,k,varargout] = modesAtFrequency(self,omega,varargin)
        % - Parameter self: concrete InternalModesBase subclass instance
        % - Parameter omega: frequency in radians per second
        % - Parameter varargin: additional requests supported by the concrete solver
        % - Returns F: horizontal-velocity mode matrix on `zOut`
        % - Returns G: vertical-velocity mode matrix on `zOut`
        % - Returns h: equivalent-depth row vector
        % - Returns k: horizontal wavenumber row vector implied by `h` and `omega`
        % - Returns varargout: additional outputs forwarded from the concrete solver
        [F,G,h,k,varargout] = modesAtFrequency(self, omega, varargin ) % Return the normal modes and eigenvalue at a given frequency.
    end
    
    methods (Abstract, Access = protected)
        self = InitializeWithBSpline(self, rho) % Used internally by subclasses to intialize with a bspline.
        self = InitializeWithGrid(self, rho, z_in) % Used internally by subclasses to intialize with a density grid.
        self = InitializeWithFunction(self, rho, z_min, z_max, z_out) % Used internally by subclasses to intialize with a density function.
        self = InitializeWithN2Function(self, N2, z_min, z_max, z_out) % Used internally by subclasses to intialize with a density function.
    end
    
    methods (Hidden)
        function [F,G,h,omega,varargout] = ModesAtWavenumber(self, k, varargin)
            if isempty(varargin)
                [F,G,h,omega] = self.modesAtWavenumber(k);
            else
                varargout = cell(size(varargin));
                [F,G,h,omega,varargout{:}] = self.modesAtWavenumber(k, varargin{:});
            end
        end
        
        function [F,G,h,k,varargout] = ModesAtFrequency(self, omega, varargin)
            if isempty(varargin)
                [F,G,h,k] = self.modesAtFrequency(omega);
            else
                varargout = cell(size(varargin));
                [F,G,h,k,varargout{:}] = self.modesAtFrequency(omega, varargin{:});
            end
        end
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Set the normalization and upper boundary condition
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.normalization(obj,norm)
            if  (norm ~= Normalization.geostrophic && norm ~= Normalization.kConstant && norm ~= Normalization.omegaConstant && norm ~= Normalization.uMax && norm ~= Normalization.wMax && norm ~= Normalization.surfacePressure)
                error('Invalid normalization! Valid options: Normalization.kConstant, Normalization.omegaConstant, Normalization.uMax, Normalization.wMax')
            else
                obj.normalization = norm;
            end
        end
        
        function set.upperBoundary(self,upperBoundary)
            self.upperBoundary = upperBoundary;
            self.upperBoundaryDidChange();
        end
        
        function value = get.zMin(self)
            value = self.zDomain(1);
        end
        
        function value = get.zMax(self)
            value = self.zDomain(2);
        end
    end
    
    methods (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesBase(options)
            % Initialize the shared solver state used by concrete subclasses.
            %
            % - Topic: Developer topics
            % - Declaration: self = InternalModesBase(options)
            % - Parameter options.rho: density profile as gridded values, a spline, or a function handle
            % - Parameter options.N2: buoyancy-frequency function handle used instead of `rho`
            % - Parameter options.zIn: input depth grid or domain bounds
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns self: initialized InternalModesBase state for subclass construction
            % - Developer: true
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            if isempty(options.zIn)
                error('You must specify zIn');
            end

            self.requiresMonotonicDensity = self.requiresMonotonicDensitySetting();
            self.zDomain = [min(options.zIn) max(options.zIn)];
            self.Lz = self.zDomain(2)-self.zDomain(1);
            self.latitude = options.latitude;
            self.rotationRate = options.rotationRate;
            self.f0 = 2*(self.rotationRate)*sin(self.latitude*pi/180);
            self.rho0 = options.rho0;
            self.nModes = options.nModes;
            self.g = options.g;

            if isempty(options.zOut)
                self.z = options.zIn;
            else
                self.z = options.zOut;
            end
            
            if ~isequal(options.N2,@disp)
                self.InitializeWithN2Function(options.N2, min(options.zIn), max(options.zIn));
            elseif isa(options.rho,'function_handle') == true
                self.InitializeWithFunction(options.rho, min(options.zIn), max(options.zIn));
            elseif isa(options.rho,'BSpline') == true
                self.rho0 = rho(max(options.rho.domain));
                self.InitializeWithBSpline(options.rho);
            elseif isa(options.rho,'numeric') == true
                if length(options.rho) ~= length(options.zIn)
                    error('rho must be 1 dimensional and z must have the same length');
                end
                self.rho0 = min(options.rho);
                options.rho = reshape(options.rho,[],1);
                [zGrid,I] = sort(options.zIn,'ascend');
                rhoGrid = options.rho(I);
                self.InitializeWithGrid(rhoGrid,zGrid);
            else
                error('You must initialize InternalModes with rho, N2, rhoGrid, or rhoSpline');
            end
            
            self.kFromOmega = @(h,omega) sqrt((omega^2 - self.f0^2)./(self.g * h));
            self.omegaFromK = @(h,k) sqrt( self.g * h * k^2 + self.f0^2 );
        end
                
        function out=requiresMonotonicDensitySetting(~)
            out=0;
        end

        function self = upperBoundaryDidChange(self)
            % This function is called when the user changes the surface
            % boundary condition. By overriding this function, a subclass
            % can respond as necessary.
        end
   
    end
end
