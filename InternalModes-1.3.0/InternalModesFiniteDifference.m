classdef InternalModesFiniteDifference < InternalModesBase
    % Solve the vertical eigenvalue problems with finite-difference matrices.
    %
    % `InternalModesFiniteDifference` discretizes the fixed-$$K$$ and
    % fixed-$$\omega$$ eigenvalue problems directly on the supplied depth
    % grid. It uses arbitrary-order finite-difference stencils generated
    % from the Fornberg weights algorithm and optionally interpolates the
    % resulting modes onto a separate output grid.
    %
    % This class is mainly useful as a baseline against the spectral
    % methods discussed in Section 3 of Early, Lelong, and Smith (2020),
    % where the manuscript explains why finite differencing becomes less
    % attractive when many frequencies or wavenumbers are required.
    %
    % ```matlab
    % im = InternalModesFiniteDifference(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, orderOfAccuracy=4);
    % [F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
    % ```
    %
    % - Topic: Create and initialize modes
    % - Topic: Inspect grids and stratification
    % - Topic: Compute modes
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesFiniteDifference < InternalModesBase
    
    properties (Access = public)
        % Density profile sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho  % Density on the z grid.
        % Buoyancy-frequency profile sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        % Formal order of accuracy for the finite-difference stencils.
        %
        % - Topic: Create and initialize modes
        orderOfAccuracy = 4 % Order of accuracy of the finite difference matrices.
    end
    
    properties (Dependent)
        % First depth derivative of the background density on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_z % First derivative of density on the z grid.
        % Second depth derivative of the background density on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_zz % Second derivative of density on the z grid.
    end
    
    properties (Access = public)
        % Number of differentiation-grid points.
        %
        % - Topic: Developer topics
        % - Developer: true
        n                   % length of z_diff
        % Depth grid used for differentiation and the generalized EVP.
        %
        % - Topic: Developer topics
        % - Developer: true
        z_diff              % the z-grid used for differentiation
        % Background density sampled on `z_diff`.
        %
        % - Topic: Developer topics
        % - Developer: true
        rho_z_diff          % rho on the z_diff grid
        % Buoyancy frequency squared sampled on `z_diff`.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2_z_diff           % N2 on the z_diff grid
        % First-derivative finite-difference matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Diff1               % 1st derivative matrix, w/ 1st derivative boundaries
        % Second-derivative finite-difference matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Diff2               % 2nd derivative matrix, w/ BCs set by upperBoundary property
        
        % Transformation from `z_diff` functions to the public output grid.
        %
        % - Topic: Developer topics
        % - Developer: true
        T_out               % *function* handle that transforms from z_diff functions to z_out functions
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesFiniteDifference(options)
            % Initialize the finite-difference solver.
            %
            % - Topic: Create and initialize modes
            % - Declaration: im = InternalModesFiniteDifference(options)
            % - Parameter options.rho: density profile as gridded values, a spline, or a function handle
            % - Parameter options.N2: buoyancy-frequency function handle used instead of `rho`
            % - Parameter options.zIn: input depth grid or domain bounds
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.orderOfAccuracy: formal order of accuracy for the finite-difference stencils
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns im: finite-difference solver instance
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.orderOfAccuracy (1,1) double {mustBePositive} = 4
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            baseOptions = rmfield(options,'orderOfAccuracy');
            baseArgs = namedargs2cell(baseOptions);
            self@InternalModesBase(baseArgs{:});
            self.orderOfAccuracy = options.orderOfAccuracy;
            
            self.n = length(self.z_diff);
            self.Diff1 = InternalModesFiniteDifference.FiniteDifferenceMatrix(1, self.z_diff, 1, 1, self.orderOfAccuracy);
            self.Diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 2, 2, self.orderOfAccuracy);
            self.N2_z_diff = -(self.g/self.rho0) * self.Diff1 * self.rho_z_diff;
            
            self.InitializeOutputTransformation(self.z);
            self.rho = self.T_out(self.rho_z_diff);
            self.N2 = self.T_out(self.N2_z_diff);
            
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = self.n;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega,varargout] = modesAtWavenumber(self, k, varargin )
            % Return finite-difference modes for a fixed horizontal wavenumber.
            %
            % - Topic: Compute modes
            % - Declaration: [F,G,h,omega,varargout] = modesAtWavenumber(self,k,varargin)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter k: horizontal wavenumber
            % - Parameter varargin: additional requests forwarded through `ModesFromGEP`
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns omega: frequency row vector implied by `h` and `k`
            % - Returns varargout: additional outputs forwarded through `ModesFromGEP`
            
            self.gridFrequency = 0;
            
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \left( \partial_{zz} - K^2*I \right)
            % B = \frac{f_0^2 - N^2}{g}
            A = self.Diff2 - k*k*eye(self.n);
            B = diag(self.f0*self.f0 - self.N2_z_diff)/self.g;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
            
            h_func = @(lambda) 1.0 ./ lambda;
            varargout = cell(size(varargin));
            [F,G,h,varargout{:}] = ModesFromGEP(self,A,B,h_func, varargin{:});
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,varargout] = modesAtFrequency(self, omega, varargin )
            % Return finite-difference modes for a fixed frequency.
            %
            % - Topic: Compute modes
            % - Declaration: [F,G,h,k,varargout] = modesAtFrequency(self,omega,varargin)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter omega: frequency in radians per second
            % - Parameter varargin: additional requests forwarded through `ModesFromGEP`
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns k: horizontal wavenumber row vector implied by `h` and `omega`
            % - Returns varargout: additional outputs forwarded through `ModesFromGEP`
            
            self.gridFrequency = omega;
            
            A = self.Diff2;
            B = -diag(self.N2_z_diff - omega*omega)/self.g;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
                        
            h_func = @(lambda) 1.0 ./ lambda;
            varargout = cell(size(varargin));
            [F,G,h,varargout{:}] = ModesFromGEP(self,A,B,h_func, varargin{:});
            k = self.kFromOmega(h,omega);
        end
        
        function [A,B] = ApplyBoundaryConditions(self,A,B)
            % Apply the current lower and upper boundary conditions to an EVP pair.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = ApplyBoundaryConditions(self,A,B)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter A: left generalized-eigenproblem matrix
            % - Parameter B: right generalized-eigenproblem matrix
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            iSurface = size(A,1);
            iBottom = 1;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip % G = 0
                    A(iBottom,:) = 0;
                    A(iBottom,iBottom) = 1;
                    B(iBottom,:) = 0;
                case LowerBoundary.noSlip % G_z = 0
                    D = weights(self.z_diff(iBottom),self.z_diff,1);
                    A(iBottom,:) = D(2,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    range = (iSurface-(self.orderOfAccuracy+1-1)):iSurface;
                    D = InternalModesFiniteDifference.weights( self.z_diff(iSurface), self.z_diff(range), 1 );
                    A(iSurface,:) = 0;
                    A(iSurface,range) = D(2,:);
                    B(iSurface,:) = 0;
                    B(iSurface,iSurface) = 1;
                case UpperBoundary.rigidLid
                    A(iSurface,:) = 0;
                    A(iSurface,iSurface) = 1;
                    B(iSurface,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
        function psi = surfaceModesAtWavenumber(self, k)
            % Return the surface SQG mode at fixed horizontal wavenumber.
            %
            % - Topic: Compute modes
            % - Declaration: psi = surfaceModesAtWavenumber(self,k)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter k: horizontal wavenumber array
            % - Returns psi: surface SQG mode evaluated on `zOut`
            psi = self.boundaryModesAtWavenumber(k,01);
        end
        
        function psi = bottomModesAtWavenumber(self, k)
            % Return the bottom SQG mode at fixed horizontal wavenumber.
            %
            % - Topic: Compute modes
            % - Declaration: psi = bottomModesAtWavenumber(self,k)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter k: horizontal wavenumber array
            % - Returns psi: bottom SQG mode evaluated on `zOut`
            psi = self.boundaryModesAtWavenumber(k,0);
        end
        
        function psi = boundaryModesAtWavenumber(self, k, isSurface)
            % Return either the surface or bottom SQG mode at fixed horizontal wavenumber.
            %
            % - Topic: Compute modes
            % - Declaration: psi = boundaryModesAtWavenumber(self,k,isSurface)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter k: horizontal wavenumber array
            % - Parameter isSurface: logical flag selecting the surface mode when true and the bottom mode when false
            % - Returns psi: requested boundary mode evaluated on `zOut`
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            
            % f'' finite diff matrix with f' at the boundaries
            diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 1, 1, self.orderOfAccuracy);
            N2z_z_diff = -(self.g/self.rho0) * diff2 * self.rho_z_diff;
            A = self.N2_z_diff .* diff2 - N2z_z_diff .* self.Diff1;
            B = - (1/(self.f0*self.f0))* (self.N2_z_diff.*self.N2_z_diff) .* eye(self.n);
            
            b = zeros(self.n,1);
            if isSurface == 1
                b(end) = 1;
            else
                b(1) = 1;
            end

            psi = zeros(length(k),self.n);
            for ii = 1:length(k)
                M = A + k(ii)*k(ii)*B;
                M(1,:) = self.f0*diff2(1,:);
                M(end,:) = self.f0*diff2(end,:);
                psi(ii,:) = M\b;
            end
            
            sizeK(end+1) = self.n;
            psi = reshape(psi,sizeK);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.Diff1 * self.rho_z_diff;
        end
        
        function value = get.rho_zz(self)
            diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 2, 2, self.orderOfAccuracy);
            value = diff2 * self.rho_z_diff;
        end
    end
    
    methods (Hidden)
        function psi = SurfaceModesAtWavenumber(self, k)
            psi = self.surfaceModesAtWavenumber(k);
        end
        
        function psi = BottomModesAtWavenumber(self, k)
            psi = self.bottomModesAtWavenumber(k);
        end
        
        function psi = BoundaryModesAtWavenumber(self, k, isSurface)
            psi = self.boundaryModesAtWavenumber(k, isSurface);
        end
    end
    
    methods (Access = protected)     
        
        function self = InitializeWithBSpline(self, rho)
           error('Not yet implemented') 
        end
        
        function self = InitializeWithGrid(self, rho, z_in)
            % Used internally by subclasses to intialize with a density function.
            %
            % Superclass calls this method upon initialization when it
            % determines that the input is given in gridded form. The goal
            % is to initialize z_diff and rho_z_diff.
            self.z_diff = z_in;
            self.rho_z_diff = rho;
        end

        function self = InitializeWithFunction(self, rho, z_min, z_max)
            % Used internally by subclasses to intialize with a density grid.
            %
            % The superclass calls this method upon initialization when it
            % determines that the input is given in functional form. The
            % goal is to initialize z_diff and rho_z_diff.
            if length(self.z) < 5
                error('You need more than 5 point output points for finite differencing to work');
            end
            
            if (min(self.z) == z_min && max(self.z) == z_max)
                self.z_diff = self.z;
                self.rho_z_diff = rho(self.z_diff);
            else
                error('Other cases not yet implemented');
                % Eventually we may want to use stretched coordinates as a
                % default
            end
        end
        
        function self = InitializeWithN2Function(self, N2, zMin, zMax)
            error('InternalModesFiniteDifference:UnsupportedInitialization', ...
                'Initialization from N2 functions is not implemented for InternalModesFiniteDifference.');
        end
    end
    
    methods (Access = public)   
        function self = InitializeOutputTransformation(self, z_out)
            % Prepare the interpolation from the differentiation grid to the public output grid.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: self = InitializeOutputTransformation(self,z_out)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter z_out: requested output grid
            % - Returns self: updated finite-difference solver instance
            if isequal(self.z_diff,z_out)
                self.T_out = @(f_in) real(f_in);
            else % want to interpolate onto the output grid
                self.T_out = @(f_in) interp1(self.z_diff,real(f_in),z_out);
            end
        end
        

        function [F,G,h,varargout] = ModesFromGEP(self,A,B,h_func, varargin)
            % Solve a generalized EVP and map its modes onto the public grid.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [F,G,h,varargout] = ModesFromGEP(self,A,B,h_func,varargin)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter A: left generalized-eigenproblem matrix
            % - Parameter B: right generalized-eigenproblem matrix
            % - Parameter h_func: map from eigenvalue to equivalent depth
            % - Parameter varargin: optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns varargout: requested diagnostics from `NormalizeModes`
            [V,D] = eig( A, B );
            
            [h, permutation] = sort(real(h_func(diag(D))),'descend');
            G = V(:,permutation);
            
            F = zeros(self.n,self.n);
            for j=1:self.n
                F(:,j) = h(j) * self.Diff1 * G(:,j);
            end
            
            if isempty(varargin)
                [F_norm,G_norm] = self.NormalizeModes(F,G,self.N2_z_diff, self.z_diff);
            else
                varargout = cell(size(varargin));
                [F_norm,G_norm,varargout{:}] = self.NormalizeModes(F,G,self.N2_z_diff, self.z_diff, varargin{:});
            end
            
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            for iMode=1:self.nModes
                F(:,iMode) = self.T_out(F_norm(:,iMode));
                G(:,iMode) = self.T_out(G_norm(:,iMode));
            end
            h = reshape(h(1:self.nModes),1,[]);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Generical function to normalize
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,varargout] = NormalizeModes(self,F,G,N2,z,varargin)
            % Normalize finite-difference modes using the active convention.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [F,G,varargout] = NormalizeModes(self,F,G,N2,z,varargin)
            % - Parameter self: InternalModesFiniteDifference instance
            % - Parameter F: horizontal-velocity mode matrix
            % - Parameter G: vertical-velocity mode matrix
            % - Parameter N2: buoyancy-frequency profile associated with the modes
            % - Parameter z: depth grid associated with the modes
            % - Parameter varargin: optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`
            % - Returns F: normalized horizontal-velocity mode matrix
            % - Returns G: normalized vertical-velocity mode matrix
            % - Returns varargout: requested quadratic-integral and normalization diagnostics
            if z(2)-z(1) > 0
                direction = 'last';
            else
                direction = 'first';
            end
            
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                varargout{iArg} = zeros(1,length(G(1,:)));
            end
            
            [maxIndexZ] = find(N2-self.gridFrequency*self.gridFrequency>0,1,direction);  
            for j=1:length(G(1,:))
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs(F(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.wMax
                        A = max( abs(G(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.kConstant
                        if z(2)-z(1) > 0
                            G20 = G(end,j)^2;
                        else
                            G20 = G(1,j)^2;
                        end
                        A = abs(G20 + trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                    case Normalization.omegaConstant
                        A = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                end
                
                if F(maxIndexZ,j)< 0
                    F(:,j) = -F(:,j);
                    G(:,j) = -G(:,j);
                end
                
                for iArg=1:length(varargin)
                    if ( strcmp(varargin{iArg}, 'F2') )
                        varargout{iArg}(j) = abs(trapz( z, F(:,j) .^ 2));
                    elseif ( strcmp(varargin{iArg}, 'G2') )
                        varargout{iArg}(j) = abs(trapz(z, G(:,j).^2));
                    elseif ( strcmp(varargin{iArg}, 'N2G2') )
                        varargout{iArg}(j) = abs(trapz(z, N2.* (G(:,j).^2)));
                    elseif  ( strcmp(varargin{iArg}, 'uMax') )
                        B = max( abs(F(:,j)) );
                        varargout{iArg}(j) = abs(1/B);
                    elseif  ( strcmp(varargin{iArg}, 'wMax') )
                        B = max( abs(G(:,j)) );
                        varargout{iArg}(j) = abs(1/B);
                    elseif ( strcmp(varargin{iArg}, 'kConstant') )
                        if z(2)-z(1) > 0
                            G20 = G(end,j)^2;
                        else
                            G20 = G(1,j)^2;
                        end
                        B = abs(G20 + trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                        varargout{iArg}(j) = sqrt(abs(1/B));
                    elseif ( strcmp(varargin{iArg}, 'omegaConstant') )
                        B = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                        varargout{iArg}(j) = sqrt(abs(1/B));
                    else
                        error('Invalid option. You may request F2, G2, N2G2');
                    end
                end
            end
        end
    end
    
    methods (Static)
        function D = FiniteDifferenceMatrix(numDerivs, x, leftBCDerivs, rightBCDerivs, orderOfAccuracy)
            % Build a finite-difference differentiation matrix with boundary stencils.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: D = FiniteDifferenceMatrix(numDerivs,x,leftBCDerivs,rightBCDerivs,orderOfAccuracy)
            % - Parameter numDerivs: derivative order to approximate
            % - Parameter x: grid on which to build the stencil matrix
            % - Parameter leftBCDerivs: derivative order enforced at the left boundary
            % - Parameter rightBCDerivs: derivative order enforced at the right boundary
            % - Parameter orderOfAccuracy: formal order of accuracy
            % - Returns D: differentiation matrix
            
            n = length(x);
            D = zeros(n,n);
            
            % left boundary condition
            range = 1:(orderOfAccuracy+leftBCDerivs); % not +1 because we're computing inclusive
            c = InternalModesFiniteDifference.weights( x(1), x(range), leftBCDerivs );
            D(1,range) = c(leftBCDerivs+1,:);
            
            % central derivatives, including possible weird end points
            centralBandwidth = ceil(numDerivs/2)+ceil(orderOfAccuracy/2)-1;
            for i=2:(n-1)
                rangeLength = 2*centralBandwidth; % not +1 because we're computing inclusive
                startIndex = max(i-centralBandwidth, 1);
                endIndex = startIndex+rangeLength;
                if (endIndex > n)
                    endIndex = n;
                    startIndex = endIndex-rangeLength;
                end
                range = startIndex:endIndex;
                c = InternalModesFiniteDifference.weights( x(i), x(range), numDerivs );
                D(i,range) = c(numDerivs+1,:);
            end
            
            % right boundary condition
            range = (n-(orderOfAccuracy+rightBCDerivs-1)):n; % not +1 because we're computing inclusive
            c = InternalModesFiniteDifference.weights( x(n), x(range), rightBCDerivs );
            D(n,range) = c(rightBCDerivs+1,:);
        end
                
        function c = weights(z,x,m)
            % Return Fornberg finite-difference weights for one stencil location.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: c = weights(z,x,m)
            % - Parameter z: evaluation location
            % - Parameter x: stencil grid points
            % - Parameter m: highest derivative order requested
            % - Returns c: weight matrix whose rows correspond to derivative order
            %
            n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
            for i=2:n
                mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
                for j=1:i-1
                    c3=x(i)-x(j);  c2=c2*c3;
                    if j==i-1
                        c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
                        c(1,i)=-c1*c5*c(1,i-1)/c2;
                    end
                    c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
                    c(1,j)=c4*c(1,j)/c3;
                end
                c1=c2;
            end
            
        end
    end
    
end
