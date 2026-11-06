classdef InternalModesDensitySpectral < InternalModesSpectral
    % Solve the vertical EVP on a density-stretched coordinate with Chebyshev collocation.
    %
    % `InternalModesDensitySpectral` implements the density-coordinate
    % formulation discussed in Section 4.2 of Early, Lelong, and Smith
    % (2020). It uses the stretched coordinate
    %
    % $$
    % s(z) = -\frac{g}{\rho_0}\rho(z) + g,
    % $$
    %
    % so that
    %
    % $$
    % F_j = h_j N^2 \partial_s G_j.
    % $$
    %
    % This coordinate concentrates points where the density varies most
    % rapidly and therefore requires a monotonic background density
    % profile.
    %
    % ```matlab
    % im = InternalModesDensitySpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
    % [F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
    % ```
    %
    % - Topic: Create and initialize modes
    % - Topic: Compute modes
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesDensitySpectral < InternalModesSpectral
    
    properties %(Access = private)
        % $$\partial_z N^2$$ sampled on the Lobatto grid in the density coordinate.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2z_xLobatto    	% (d/dz)N2 on the z_sLobatto grid   
    end

    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesDensitySpectral(options)
            % Initialize the density-coordinate spectral solver.
            %
            % - Topic: Create and initialize modes
            % - Declaration: im = InternalModesDensitySpectral(options)
            % - Parameter options.rho: density profile as gridded values, a spline, or a function handle
            % - Parameter options.N2: buoyancy-frequency function handle used instead of `rho`
            % - Parameter options.zIn: input depth grid or domain bounds
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.nEVP: number of collocation points in the density-coordinate EVP
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns im: density-coordinate spectral solver instance
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.nEVP = 512
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            parentArgs = namedargs2cell(options);
            self@InternalModesSpectral(parentArgs{:});
        end
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [A,B] = eigenmatricesForWavenumber(self, k )
            % Assemble the fixed-$$K$$ generalized EVP in density coordinates.
            %
            % - Topic: Compute modes
            % - Declaration: [A,B] = eigenmatricesForWavenumber(self,k)
            % - Parameter self: InternalModesDensitySpectral instance
            % - Parameter k: horizontal wavenumber
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = eigenmatricesForFrequency(self, omega )
            % Assemble the fixed-$$\omega$$ generalized EVP in density coordinates.
            %
            % - Topic: Compute modes
            % - Declaration: [A,B] = eigenmatricesForFrequency(self,omega)
            % - Parameter self: InternalModesDensitySpectral instance
            % - Parameter omega: frequency in radians per second
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = ApplyBoundaryConditions(self,A,B)
            % Apply the active boundary conditions in density coordinates.
            %
            % With a free surface, this enforces the manuscript relation
            % $$N^2 \partial_s G_j = G_j / h_j$$ at the surface.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = ApplyBoundaryConditions(self,A,B)
            % - Parameter self: InternalModesDensitySpectral instance
            % - Parameter A: left generalized-eigenproblem matrix
            % - Parameter B: right generalized-eigenproblem matrix
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            n = self.nEVP;
            
            iSurface = 1;
            iBottom = n;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip
                    A(iBottom,:) = T(iBottom,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.noSlip
                    A(iBottom,:) = Tz(iBottom,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    A(iSurface,:) = self.N2_xLobatto(iSurface) * Tz(iSurface,:);
                    B(iSurface,:) = T(iSurface,:);
                case UpperBoundary.rigidLid
                    A(iSurface,:) = T(iSurface,:);
                    B(iSurface,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
    end
    
    methods (Access = protected)
        function out = requiresMonotonicDensitySetting(~)
            out = 1;
        end
        
        function self = InitializeWithGrid(self, rho, zIn)
            InitializeWithGrid@InternalModesSpectral(self,rho,zIn);
            
            % This is our new (s)treched grid...we need to make s(z) be
            % monotonic!
            s = @(z) (-self.g/self.rho0)*self.rho_function(z) + self.g;
            
            self.xDomain = [s(self.zMin) s(self.zMax)];
            
            % if the data is noisy, let's make sure that there are no
            % variations below the grid scale. Specifically, each grid
            % point should uniquely map to a z-grid point.
            if any(diff(rho)./diff(zIn) > 0)
                n = self.nEVP;
                xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
                xIn = s(zIn);
                y = discretize(xLobatto,xIn);
                while (length(unique(y)) < length(xLobatto))
                    n = n-1;
                    xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
                    y = discretize(xLobatto,xIn);
                end
                
                K = 6; % quintic spline
                splineDegree = K - 1;
                z_knot = BSpline.knotPointsForDataPoints([zIn(1); zIn(unique(y)+1)], S=splineDegree);
                rho_interpolant = ConstrainedSpline.fromData(zIn, rho, ...
                    S=splineDegree, ...
                    knotPoints=z_knot, ...
                    distribution=NormalDistribution(sigma=1), ...
                    constraints=GlobalConstraint.monotonicDecreasing());
                
                self.rho_function = rho_interpolant;
                self.N2_function = (-self.g/self.rho0)*diff(self.rho_function);
                
                if self.shouldShowDiagnostics == 1
                    fprintf('Data was found to be noise. Creating a %d-order monotonic smoothing spline using %d knot points.\n', K, n);
                end
            end
        end
        
        function self = SetupEigenvalueProblem(self)     
            % Create a stretched grid from the density function
            self.x_function = @(z) (-self.g/self.rho0)*self.rho_function(z) + self.g;
            
            N2z_function = diff(self.N2_function);
            self.N2z_xLobatto = N2z_function(self.z_xLobatto);

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromVCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromVCheb = @(G_cheb,h) h * self.N2 .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.GFromVCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            self.FFromVCheb = @(G_cheb,h) h * self.N2_xLobatto .* InternalModesSpectral.ifct(self.Diff1_xCheb(G_cheb));
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (1 - self.f0*self.f0./self.N2_xLobatto) .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2)./self.N2_xLobatto)));
        end
    end
    
    methods (Hidden)
        function [A,B] = EigenmatricesForWavenumber(self, k)
            [A,B] = self.eigenmatricesForWavenumber(k);
        end
        
        function [A,B] = EigenmatricesForFrequency(self, omega)
            [A,B] = self.eigenmatricesForFrequency(omega);
        end
    end
    
end
