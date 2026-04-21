classdef InternalModesWKBSpectral < InternalModesSpectral
    % Solve the vertical EVP on a WKB stretched coordinate with Chebyshev collocation.
    %
    % `InternalModesWKBSpectral` implements the WKB-coordinate spectral
    % method described in Section 4.3 of Early, Lelong, and Smith (2020).
    % It introduces the stretched coordinate
    %
    % $$
    % s(z) = \int_z^D N(z') \, dz',
    % $$
    %
    % and solves the transformed fixed-$$K$$ and fixed-$$\omega$$
    % eigenproblems in $$s$$, with
    %
    % $$
    % F_j = h_j N \, \partial_s G_j.
    % $$
    %
    % Compared with `InternalModesSpectral`, this class concentrates grid
    % resolution where stratification is strong while preserving the
    % public constructor contract used by downstream packages.
    %
    % ```matlab
    % im = InternalModesWKBSpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
    % [F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
    % ```
    %
    % - Topic: Create and initialize modes
    % - Topic: Compute modes
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesWKBSpectral < InternalModesSpectral
    
    properties %(Access = private)
        % Derivative of $$N(z)$$ used when assembling the stretched-coordinate EVP.
        %
        % - Topic: Developer topics
        % - Developer: true
        Nz_function
        % $$\partial_z N$$ sampled on the Lobatto grid in the WKB coordinate.
        %
        % - Topic: Developer topics
        % - Developer: true
        Nz_xLobatto     	% (d/dz)N on the xiLobatto grid   
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKBSpectral(options)
            % Initialize the WKB-coordinate spectral solver.
            %
            % This constructor intentionally preserves the public
            % name-value API used by `wave-vortex-model`.
            %
            % - Topic: Create and initialize modes
            % - Declaration: im = InternalModesWKBSpectral(options)
            % - Parameter options.rho: density profile as gridded values, a spline, or a function handle
            % - Parameter options.N2: buoyancy-frequency function handle used instead of `rho`
            % - Parameter options.zIn: input depth grid or domain bounds
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.nEVP: number of collocation points in the stretched-coordinate EVP
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns im: WKB-coordinate spectral solver instance
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.nEVP = 512;
                options.rotationRate (1,1) double = 7.2921e-5;
                options.g (1,1) double = 9.81
            end
            self@InternalModesSpectral(rho=options.rho,N2=options.N2,zIn=options.zIn,zOut=options.zOut,latitude=options.latitude,rho0=options.rho0,nModes=options.nModes,nEVP=options.nEVP,rotationRate=options.rotationRate,g=options.g);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [A,B] = eigenmatricesForWavenumber(self, k )
            % Assemble the fixed-$$K$$ generalized EVP in the WKB coordinate.
            %
            % In manuscript notation this method discretizes
            %
            % $$
            % N^2 \partial_{ss} G_j + N_s \partial_s G_j - K^2 G_j
            % = \frac{f_0^2-N^2}{g h_j} G_j.
            % $$
            %
            % - Topic: Compute modes
            % - Declaration: [A,B] = eigenmatricesForWavenumber(self,k)
            % - Parameter self: InternalModesWKBSpectral instance
            % - Parameter k: horizontal wavenumber
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
                        
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end
        
        function [A,B] = eigenmatricesForFrequency(self, omega )
            % Assemble the fixed-$$\omega$$ generalized EVP in the WKB coordinate.
            %
            % - Topic: Compute modes
            % - Declaration: [A,B] = eigenmatricesForFrequency(self,omega)
            % - Parameter self: InternalModesWKBSpectral instance
            % - Parameter omega: frequency in radians per second
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = EigenmatricesForMDAModes(self )
            % Assemble the MDA generalized EVP used by the diagnostic helper.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = EigenmatricesForMDAModes(self)
            % - Parameter self: InternalModesWKBSpectral instance
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( - self.N2_xLobatto/self.g )*T;

            % upper-boundary
            A(1,:) = Tz(1,:); %-Tz(n,:);
            B(1,:) = 0 ;%1/self.Lz; %0*T(n,:);
            self.upperBoundary = UpperBoundary.mda;

            % lower-boundary
            A(n,:) = Tz(n,:); %self.Lz*Tz(n,:)-T(n,:);
            B(n,:) = 0; %1/self.Lz; %0*T(n,:);
            self.lowerBoundary = LowerBoundary.mda;
        end

        function [A,B] = EigenmatricesForGeostrophicGModes(self, k )
            % Assemble the geostrophic `G`-mode EVP in the WKB coordinate.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = EigenmatricesForGeostrophicGModes(self,k)
            % - Parameter self: InternalModesWKBSpectral instance
            % - Parameter k: horizontal wavenumber
            % - Returns A: left generalized-eigenproblem matrix
            % - Returns B: right generalized-eigenproblem matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( (- self.N2_xLobatto)/self.g )*T;

            A(n,:) = T(n,:);
            B(n,:) = 0;

            if (self.g/(self.f0*self.f0))*(k*k) < 1
                A(1,:) = sqrt(self.N2_xLobatto(1))*Tz(1,:) + (self.g/(self.f0*self.f0))*(k*k)*T(1,:);
            else
                A(1,:) = sqrt(self.N2_xLobatto(1))*(self.f0*self.f0)/(self.g *k*k)* Tz(1,:) + T(1,:);
            end
            B(1,:) = 0;
        end
                
        function [A,B] = ApplyBoundaryConditions(self,A,B)
            % Apply the active surface and bottom conditions in WKB coordinates.
            %
            % With a free surface, this enforces
            % $N \partial_s G_j = G_j / h_j$ at `z = 0`, while rigid-lid
            % and bottom conditions follow the same physical choices as
            % the base spectral solver.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = ApplyBoundaryConditions(self,A,B)
            % - Parameter self: InternalModesWKBSpectral instance
            % - Parameter A: left generalized-eigenproblem matrix
            % - Parameter B: right generalized-eigenproblem matrix
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            n = self.nEVP;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip
                    A(n,:) = T(n,:);
                    B(n,:) = 0;
                case LowerBoundary.noSlip
                    A(n,:) = Tz(n,:);
                    B(n,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % N*G_s = \frac{1}{h_j} G at the surface
                    A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
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
            N_function = sqrt(self.N2_function);
            s = cumsum(N_function);
            
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
                
                K = 4; % cubic spline
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
            % Create the stretched WKB grid
            N_function = sqrt(self.N2_function);
            self.x_function = cumsum(N_function);
            
            self.Nz_function = diff(N_function);
            self.Nz_xLobatto = self.Nz_function(self.z_xLobatto);

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromVCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromVCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.GFromVCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            self.FFromVCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            self.GeostrophicNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * self.N2_xLobatto .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
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
