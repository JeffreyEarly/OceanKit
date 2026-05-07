classdef InternalModesConstantStratification < InternalModesBase
    % Solve the vertical mode problem for constant buoyancy frequency.
    %
    % `InternalModesConstantStratification` provides the closed-form
    % solutions used throughout the package for validation and smoke
    % testing. In this model,
    %
    % $$
    % N^2(z) = N_0^2
    % $$
    %
    % and the background density is the linear profile
    %
    % $$
    % \bar{\rho}(z) = \rho_0 - \frac{N_0^2 \rho_0}{g} z.
    % $$
    %
    % The fixed-$$K$$ and fixed-$$\omega$$ solutions reduce to trigonometric
    % or hyperbolic functions with analytically known equivalent depths,
    % making this class the main reference implementation for checking the
    % numerical solvers.
    %
    % ```matlab
    % im = InternalModesConstantStratification(N0=5.2e-3, zIn=[-5000 0], zOut=zOut, latitude=33);
    % [F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
    % ```
    %
    % - Topic: Create and initialize modes
    % - Topic: Inspect grids and stratification
    % - Topic: Compute modes
    % - Topic: Inspect analytical solutions
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesConstantStratification < InternalModesBase
    properties (Access = public)
        % Constant buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect grids and stratification
        N0
        % Density profile sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho
        % Constant buoyancy-frequency profile sampled on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        N2
        % First depth derivative of the background density on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_z
        % Second depth derivative of the background density on `zOut`.
        %
        % - Topic: Inspect grids and stratification
        rho_zz
    end
    
    methods
        function self = InternalModesConstantStratification(options) 
            % Initialize the constant-stratification analytical solver.
            %
            % - Topic: Create and initialize modes
            % - Declaration: im = InternalModesConstantStratification(options)
            % - Parameter options.N0: constant buoyancy frequency in radians per second
            % - Parameter options.zIn: two-element depth domain `[zBottom zSurface]`
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns im: constant-stratification solver instance
            arguments
                options.N0 (1,1) double = 5.2e-3
                options.zIn (1,2) double = [-1300 0]
                options.zOut (:,1) double = linspace(-1300,0,65).'
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            rho0 = options.rho0;
            N0 = options.N0;
            rhoFunction = @(z) -(N0*N0*rho0/options.g)*z + rho0;
            N2Function = @(z) N0*N0*ones(size(z));
            
            
            self@InternalModesBase(rho=rhoFunction,zIn=options.zIn,zOut=options.zOut,latitude=options.latitude,rho0=rho0,nModes=options.nModes,rotationRate=options.rotationRate,g=options.g);
            self.N0 = N0;
            self.rho = rhoFunction(self.z);
            self.N2 = N2Function(self.z);
            self.rho_z = -(self.N0*self.N0*self.rho0/self.g)*ones(size(self.z));
            self.rho_zz = zeros(size(self.z));
            
            % fprintf('Using the analytical form for constant stratification N0=%.7g\n',self.N0);
        end
                
        function [F,G,h,omega,varargout] = modesAtWavenumber(self, k, varargin )
            % Return the analytical modes for a fixed horizontal wavenumber.
            %
            % - Topic: Compute modes
            % - Declaration: [F,G,h,omega,varargout] = modesAtWavenumber(self,k,varargin)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter k: horizontal wavenumber
            % - Parameter varargin: optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns omega: frequency row vector implied by `h` and `k`
            % - Returns varargout: requested normalization and quadratic-integral diagnostics
            k_z = (1:self.nModes)*pi/self.Lz;
            if self.upperBoundary == UpperBoundary.freeSurface % add the free surface correction to the vertical wavenumber
                for i=1:self.nModes
                    f = @(xi) (xi+i*pi)*(self.N0*self.N0 - self.f0*self.f0)*self.Lz - self.g*(k*k*self.Lz*self.Lz+(xi+i*pi).*(xi+i*pi)).*tan(xi);
                    k_z(i) = k_z(i) + fzero(f,0)/self.Lz;
                end
            end
            h = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k+k_z.*k_z)); 
            
            [F,G,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = self.BaroclinicModesWithEigenvalue(k_z,h);
            
            if self.upperBoundary == UpperBoundary.freeSurface                
                [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = self.BarotropicModeAtWavenumber(k);
                h = cat(2,h0,h(1:end-1));
                F = cat(2,F0,F(:,1:end-1));
                G = cat(2,G0,G(:,1:end-1));
                F2 = cat(2,F20,F2(1:end-1));
                G2 = cat(2,G20,G2(1:end-1));
                N2G2 = cat(2,N2G20,N2G2(1:end-1));
                uMaxRatio = cat(2,uMaxRatio0,uMaxRatio(1:end-1));
                wMaxRatio = cat(2,wMaxRatio0,wMaxRatio(1:end-1));
                kConstantRatio = cat(2,kConstantRatio0,kConstantRatio(1:end-1));
                omegaConstantRatio = cat(2,omegaConstantRatio0,omegaConstantRatio(1:end-1));
            end

            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                if ( strcmp(varargin{iArg}, 'F2') )
                    varargout{iArg} = F2;
                elseif ( strcmp(varargin{iArg}, 'G2') )
                    varargout{iArg} = G2;
                elseif ( strcmp(varargin{iArg}, 'N2G2') )
                    varargout{iArg} = N2G2;
                elseif  ( strcmp(varargin{iArg}, 'uMax') )
                    varargout{iArg} = uMaxRatio;
                elseif  ( strcmp(varargin{iArg}, 'wMax') )
                    varargout{iArg} = wMaxRatio;
                elseif  ( strcmp(varargin{iArg}, 'kConstant') )
                    varargout{iArg} = kConstantRatio;
                elseif  ( strcmp(varargin{iArg}, 'omegaConstant') )
                    varargout{iArg} = omegaConstantRatio;
                else
                    error('Invalid option. You may request F2, G2, N2G2, uMax, wMax, kConstant, omegaConstant');
                end
            end

            omega = self.omegaFromK(h,k);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G,h,k,varargout] = modesAtFrequency(self, omega, varargin )
            % Return the analytical modes for a fixed frequency.
            %
            % - Topic: Compute modes
            % - Declaration: [F,G,h,k,varargout] = modesAtFrequency(self,omega,varargin)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter omega: frequency in radians per second
            % - Parameter varargin: optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns k: horizontal wavenumber row vector implied by `h` and `omega`
            % - Returns varargout: requested normalization and quadratic-integral diagnostics
            k_z = (1:self.nModes)*pi/self.Lz;
            if self.upperBoundary == UpperBoundary.freeSurface % add the free surface correction to the vertical wavenumber
                for i=1:self.nModes
                    f = @(xi) self.g*tan(xi)/(xi+i*pi) - (self.N0*self.N0 - omega*omega)*self.Lz/((xi+i*pi)^2);
                    k_z(i) = k_z(i) + fzero(f,0)/self.Lz;
                end
            end
            h = (self.N0*self.N0 - omega*omega)./(self.g * k_z.*k_z);
            if (omega >= self.N0)
               k_z = zeros(size(k_z));
               h = zeros(size(h));
            end
            
            [F,G,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = self.BaroclinicModesWithEigenvalue(k_z,h);
            
            if self.upperBoundary == UpperBoundary.freeSurface                
                [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = self.BarotropicModeAtFrequency(omega);
                h = cat(2,h0,h(1:end-1));
                F = cat(2,F0,F(:,1:end-1));
                G = cat(2,G0,G(:,1:end-1));
                F2 = cat(2,F20,F2(1:end-1));
                G2 = cat(2,G20,G2(1:end-1));
                N2G2 = cat(2,N2G20,N2G2(1:end-1));
                uMaxRatio = cat(2,uMaxRatio0,uMaxRatio(1:end-1));
                wMaxRatio = cat(2,wMaxRatio0,wMaxRatio(1:end-1));
                kConstantRatio = cat(2,kConstantRatio0,kConstantRatio(1:end-1));
                omegaConstantRatio = cat(2,omegaConstantRatio0,omegaConstantRatio(1:end-1));
            end

            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                if ( strcmp(varargin{iArg}, 'F2') )
                    varargout{iArg} = F2;
                elseif ( strcmp(varargin{iArg}, 'G2') )
                    varargout{iArg} = G2;
                elseif ( strcmp(varargin{iArg}, 'N2G2') )
                    varargout{iArg} = N2G2;
                elseif  ( strcmp(varargin{iArg}, 'uMax') )
                    varargout{iArg} = uMaxRatio;
                elseif  ( strcmp(varargin{iArg}, 'wMax') )
                    varargout{iArg} = wMaxRatio;
                elseif  ( strcmp(varargin{iArg}, 'kConstant') )
                    varargout{iArg} = kConstantRatio;
                elseif  ( strcmp(varargin{iArg}, 'omegaConstant') )
                    varargout{iArg} = omegaConstantRatio;
                else
                    error('Invalid option. You may request F2, G2, N2G2, uMax, wMax, kConstant, omegaConstant');
                end
            end

            k = self.kFromOmega(h,omega);
        end
        
        function [psi] = surfaceModesAtWavenumber(self, k)
            % Return the analytical surface SQG mode for constant stratification.
            %
            % - Topic: Compute modes
            % - Declaration: psi = surfaceModesAtWavenumber(self,k)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter k: horizontal wavenumber array
            % - Returns psi: surface SQG mode evaluated on `zOut`
            % size(psi) = [size(k); length(z)]
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            k = reshape(k,[],1);
            zCoord = reshape(self.z,1,[]);
            
            lambda = k * (self.N0/self.f0);
            % dividing by sinh is not well conditioned.
            % psi = (1/self.N0) * cosh( (zCoord + self.Lz) .* lambda ) ./ (k .* sinh(lambda*self.Lz));
            psi = (1/self.N0) * (exp( zCoord .* lambda ) + exp( -(zCoord + 2*self.Lz) .* lambda ) ) ./ (k .* (1-exp(-2*lambda*self.Lz)));
            
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        function [psi] = bottomModesAtWavenumber(self, k)
            % Return the analytical bottom SQG mode for constant stratification.
            %
            % - Topic: Compute modes
            % - Declaration: psi = bottomModesAtWavenumber(self,k)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter k: horizontal wavenumber array
            % - Returns psi: bottom SQG mode evaluated on `zOut`
            % size(psi) = [size(k); length(z)]
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            k = reshape(k,[],1);
            zCoord = reshape(self.z,1,[]);
            
            lambda = k * (self.N0/self.f0);           
            psi = -(1/self.N0) * (exp( (zCoord - self.Lz) .* lambda ) + exp( -(zCoord + self.Lz) .* lambda ) ) ./ (k .* (1-exp(-2*lambda*self.Lz)));
            
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BaroclinicModesWithEigenvalue(self, k_z, h)
            % Evaluate the analytical baroclinic mode shapes for given eigenvalues.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: [F,G,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BaroclinicModesWithEigenvalue(self,k_z,h)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter k_z: vertical wavenumber row vector
            % - Parameter h: equivalent-depth row vector
            % - Returns F: horizontal-velocity mode matrix
            % - Returns G: vertical-velocity mode matrix
            % - Returns F2: depth integrals of `F.^2`
            % - Returns G2: depth integrals of `G.^2`
            % - Returns N2G2: depth integrals of `N2 .* G.^2`
            % - Returns uMaxRatio: ratio from the active normalization to `uMax`
            % - Returns wMaxRatio: ratio from the active normalization to `wMax`
            % - Returns kConstantRatio: ratio from the active normalization to `kConstant`
            % - Returns omegaConstantRatio: ratio from the active normalization to `omegaConstant`
            N0_ = self.N0; % reference buoyancy frequency, radians/seconds
            A = self.BaroclinicModeNormalization(self.normalization,k_z);
            G = A .*  sin(k_z .* (self.z + self.Lz));
            F = A .*  repmat(h.*k_z,length(self.z),1) .* cos(k_z .* (self.z + self.Lz));
            F2 = A.*A.*h.*h.*k_z.*k_z*self.Lz/2;
            N2G2 = A.*A*N0_*N0_*self.Lz/2;
            G2 = A.*A*self.Lz/2;
            uMaxRatio = self.BaroclinicModeNormalization(Normalization.uMax, k_z, h)./A;
            wMaxRatio = self.BaroclinicModeNormalization(Normalization.wMax, k_z, h)./A;
            kConstantRatio = self.BaroclinicModeNormalization(Normalization.kConstant, k_z, h)./A;
            omegaConstantRatio = self.BaroclinicModeNormalization(Normalization.omegaConstant, k_z, h)./A;
        end

        function A = BaroclinicModeNormalization(self, norm, k_z, h)
            % Return the requested analytical normalization factor for baroclinic modes.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: A = BaroclinicModeNormalization(self,norm,k_z,h)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter norm: normalization convention from `Normalization`
            % - Parameter k_z: vertical wavenumber row vector
            % - Parameter h: equivalent-depth row vector
            % - Returns A: normalization factor row vector
            j = 1:self.nModes;
            switch norm
                case Normalization.kConstant
                    A = (-1).^j .* (sin(k_z*self.Lz).^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(1-sin(2*k_z*self.Lz)./(2*k_z*self.Lz))).^(-1/2);
                case Normalization.omegaConstant
                    A = (-1).^j./( h .* k_z .* sqrt(1/2 + sin(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                case Normalization.wMax
                    A = (-1).^j;
                case Normalization.uMax
                    A = (-1).^j./(h.*k_z);
                case Normalization.surfacePressure
                    A = 1./(h.*k_z.*cos(k_z*self.Lz));
            end
        end
        
        function [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = BarotropicModeAtWavenumber(self, k)
            % Return the analytical barotropic mode branch for fixed $$K$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = BarotropicModeAtWavenumber(self,k)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter k: horizontal wavenumber
            % - Returns F0: barotropic horizontal mode
            % - Returns G0: barotropic vertical mode
            % - Returns h0: barotropic equivalent depth
            % - Returns F20: depth integral of `F0.^2`
            % - Returns G20: depth integral of `G0.^2`
            % - Returns N2G20: depth integral of `N2 .* G0.^2`
            % - Returns uMaxRatio0: ratio from the active normalization to `uMax`
            % - Returns wMaxRatio0: ratio from the active normalization to `wMax`
            % - Returns kConstantRatio0: ratio from the active normalization to `kConstant`
            % - Returns omegaConstantRatio0: ratio from the active normalization to `omegaConstant`
            k_star = sqrt( (self.N0*self.N0 - self.f0*self.f0)/(self.g*self.Lz) );
                
            if (abs(k-k_star)/k_star < 1e-6) % transition (linear) solution
                solutionType = 'linear';
                h0 = self.Lz;
                k_z = 0;
            elseif k > k_star % hyperbolic solution
                solutionType = 'hyperbolic';
                f = @(q) self.Lz*(self.N0*self.N0 - self.f0*self.f0) - (1./q).*(self.g*(k*k*self.Lz*self.Lz-q.*q)).*tanh(q);
                k_initial = sqrt( k*k*self.Lz*self.Lz - self.Lz*(self.N0*self.N0 - self.f0*self.f0)/self.g);
                k_z = fzero(f, k_initial)/self.Lz;
                h0 = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k - k_z*k_z ));
            elseif k < k_star % trig solution
                solutionType = 'trig';
                f = @(q) self.Lz*(self.N0*self.N0 - self.f0*self.f0) - (1./q).*(self.g*(k*k*self.Lz*self.Lz+q.*q)).*tan(q);
                k_initial = sqrt( - k*k*self.Lz*self.Lz + self.Lz*(self.N0*self.N0 - self.f0*self.f0)/self.g);
                k_z = fzero(f, k_initial)/self.Lz;
                h0 = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k + k_z*k_z ));        
            end
            
            [F0,G0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = self.BarotropicMode(solutionType, k_z, h0);
        end
        
        function [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = BarotropicModeAtFrequency(self, omega)
            % Return the analytical barotropic mode branch for fixed $$\omega$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = BarotropicModeAtFrequency(self,omega)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter omega: frequency in radians per second
            % - Returns F0: barotropic horizontal mode
            % - Returns G0: barotropic vertical mode
            % - Returns h0: barotropic equivalent depth
            % - Returns F20: depth integral of `F0.^2`
            % - Returns G20: depth integral of `G0.^2`
            % - Returns N2G20: depth integral of `N2 .* G0.^2`
            % - Returns uMaxRatio0: ratio from the active normalization to `uMax`
            % - Returns wMaxRatio0: ratio from the active normalization to `wMax`
            % - Returns kConstantRatio0: ratio from the active normalization to `kConstant`
            % - Returns omegaConstantRatio0: ratio from the active normalization to `omegaConstant`
            if (abs(omega-self.N0)/self.N0 < 1e-6)
                solutionType = 'linear';
                h0 = self.Lz;
                k_z = 0;
            elseif omega > self.N0 % hyperbolic solution
                solutionType = 'hyperbolic';
                f = @(q) self.Lz*(omega*omega - self.N0*self.N0) - self.g * q .* tanh(q);
                k_z = fzero(f, sqrt(self.Lz*(omega*omega - self.N0*self.N0)/self.g) )/self.Lz;
                h0 = (omega*omega - self.N0*self.N0)./(self.g*k_z*k_z);
            elseif omega < self.N0 % trig solution
                solutionType = 'trig';
                f = @(q) self.Lz*(self.N0*self.N0 - omega*omega) - self.g * q .* tan(q);
                k_z = fzero(f, sqrt(self.Lz*(self.N0*self.N0 - omega*omega)/self.g) )/self.Lz;
                h0 = (self.N0*self.N0 - omega*omega)./(self.g * k_z.*k_z);
            end

            [F0,G0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = self.BarotropicMode(solutionType, k_z, h0);
        end

        function [F0,G0,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BarotropicMode(self, solutionType, k_z, h0)
            % Evaluate a chosen analytical barotropic branch.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [F0,G0,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BarotropicMode(self,solutionType,k_z,h0)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter solutionType: one of `linear`, `hyperbolic`, or `trig`
            % - Parameter k_z: barotropic vertical wavenumber
            % - Parameter h0: barotropic equivalent depth
            % - Returns F0: barotropic horizontal mode
            % - Returns G0: barotropic vertical mode
            % - Returns F2: depth integral of `F0.^2`
            % - Returns G2: depth integral of `G0.^2`
            % - Returns N2G2: depth integral of `N2 .* G0.^2`
            % - Returns uMaxRatio: ratio from the active normalization to `uMax`
            % - Returns wMaxRatio: ratio from the active normalization to `wMax`
            % - Returns kConstantRatio: ratio from the active normalization to `kConstant`
            % - Returns omegaConstantRatio: ratio from the active normalization to `omegaConstant`
            A = self.BarotropicModeNormalization(self.normalization, solutionType, k_z, h0);

            % It's safer to do a switch on solutionType, rather than check
            % that omega *or* k are equal to N0, k_star within tolerance.
            if strcmp(solutionType, 'linear')
                G0 = A*(self.z + self.Lz);
                F0 = A*self.Lz*ones(size(self.z));
                F2 = A*A*(self.Lz)^3;
                G2 = A*A*((self.Lz)^3)/3;
                N2G2 = A*A*self.N0*self.N0*((self.Lz)^3)/3;
            elseif strcmp(solutionType, 'hyperbolic')
                G0 = A*sinh(k_z*(self.z + self.Lz));
                F0 = A*h0*k_z*cosh(k_z*(self.z + self.Lz));
                F2 = A.*A.*h0.*h0.*k_z.*k_z*(sinh(2*k_z*self.Lz)./(4*k_z*self.Lz) - 1/2);
                G2 = A.*A*(sinh(2*k_z*self.Lz)/(4*k_z) - self.Lz/2);
                N2G2 = A.*A*self.N0*self.N0*(sinh(2*k_z*self.Lz)/(4*k_z) - self.Lz/2);
            elseif strcmp(solutionType, 'trig')
                G0 = A*sin(k_z*(self.z + self.Lz));
                F0 = A*h0*k_z*cos(k_z*(self.z + self.Lz));
                F2 = A.*A.*h0.*h0.*k_z.*k_z*(1/2 + sin(2*k_z*self.Lz)./(4*k_z*self.Lz));
                G2 = A.*A*(self.Lz/2 - sin(2*k_z*self.Lz)/(4*k_z));
                N2G2 = A.*A*self.N0*self.N0*(self.Lz/2 - sin(2*k_z*self.Lz)/(4*k_z));
            end
        
            uMaxRatio = self.BarotropicModeNormalization(Normalization.uMax, solutionType, k_z, h0)/A;
            wMaxRatio = self.BarotropicModeNormalization(Normalization.wMax, solutionType, k_z, h0)/A;
            kConstantRatio = self.BarotropicModeNormalization(Normalization.kConstant, solutionType, k_z, h0)/A;
            omegaConstantRatio = self.BarotropicModeNormalization(Normalization.omegaConstant, solutionType, k_z, h0)/A;
        end

        function A = BarotropicModeNormalization(self, norm, solutionType, k_z, h0)
            % Return the requested analytical normalization factor for the barotropic branch.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: A = BarotropicModeNormalization(self,norm,solutionType,k_z,h0)
            % - Parameter self: InternalModesConstantStratification instance
            % - Parameter norm: normalization convention from `Normalization`
            % - Parameter solutionType: one of `linear`, `hyperbolic`, or `trig`
            % - Parameter k_z: barotropic vertical wavenumber
            % - Parameter h0: barotropic equivalent depth
            % - Returns A: normalization factor
            if strcmp(solutionType, 'linear')
                switch norm
                    case Normalization.kConstant
                        A = 1/(self.Lz * sqrt(1 + (self.N0*self.N0 - self.f0*self.f0)*self.Lz/(2*self.g)));
                    case Normalization.omegaConstant
                        A = 1/self.Lz;
                    case Normalization.wMax
                        A = 1/self.Lz;
                    case Normalization.uMax
                        A = 1/self.Lz;
                    case Normalization.surfacePressure
                        A = 1/self.Lz;
                end
            elseif strcmp(solutionType, 'hyperbolic')
                switch norm
                    case Normalization.kConstant
                        A = (sinh(k_z*self.Lz)^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(sinh(2*k_z*self.Lz)/(2*k_z*self.Lz)-1)).^(-1/2);
                    case Normalization.omegaConstant
                        A = 1/( h0 * k_z * sqrt(1/2 + sinh(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                    case Normalization.wMax
                        A = 1/sinh(k_z*self.Lz);
                    case Normalization.uMax
                        A = 1/(h0*k_z*cosh(k_z*self.Lz));
                    case Normalization.surfacePressure
                        A = 1/(h0*k_z*cosh(k_z*self.Lz));
                end
            elseif strcmp(solutionType, 'trig')
                switch norm
                    case Normalization.kConstant
                        A = (sin(k_z*self.Lz)^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(1-sin(2*k_z*self.Lz)/(2*k_z*self.Lz))).^(-1/2);
                    case Normalization.omegaConstant
                        A = 1/( h0 * k_z * sqrt(1/2 + sin(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                    case Normalization.wMax
                        A = 1/sin(k_z*self.Lz);
                    case Normalization.uMax
                        A = 1/(h0*k_z);
                    case Normalization.surfacePressure
                        A = 1/(h0*k_z*cos(k_z*self.Lz));
                end
            end
        end
    end
    
    methods (Hidden)
        function psi = SurfaceModesAtWavenumber(self, k)
            psi = self.surfaceModesAtWavenumber(k);
        end
        
        function psi = BottomModesAtWavenumber(self, k)
            psi = self.bottomModesAtWavenumber(k);
        end
    end
    
    methods (Access = protected)
        function self = InitializeWithGrid(self, rho, zIn)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = floor(length(self.z));
            end
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = floor(length(self.z));
            end
        end
        
        function self = InitializeWithBSpline(self, rho, z_in)
            error('not yet implemented')
        end
        
                
        function self = InitializeWithN2Function(self, N2, zMin, zMax)
           error('Invalid initialization path'); 
        end
    end
    
    methods (Static)
        function flag = IsStratificationConstant(rho,z_in)
            % Test whether a supplied profile is close to constant stratification.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: flag = IsStratificationConstant(rho,z_in)
            % - Parameter rho: density profile as gridded values, a spline, or a function handle
            % - Parameter z_in: depth grid or domain bounds associated with `rho`
            % - Returns flag: logical scalar indicating whether the profile is approximately constant-stratification
            if isa(rho,'function_handle') || isa(rho,'BSpline') == true
%                 if numel(z_in) ~= 2
%                     error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
%                 end
                z=linspace(min(z_in),max(z_in),100)';
                drhodz = diff(rho(z))./diff(z);
            elseif isa(rho,'numeric') == true
                drhodz = diff(rho)./diff(z_in);
            end
            
            max_ddrho = std(drhodz)/abs(mean(drhodz));
            flag = max_ddrho < 3e-3;
        end
        
        function [N0, rho0] = BuoyancyFrequencyFromConstantStratification(rho,z_in)
            % Estimate `N0` and `rho0` from a constant-stratification profile.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: [N0,rho0] = BuoyancyFrequencyFromConstantStratification(rho,z_in)
            % - Parameter rho: density profile as gridded values, a spline, or a function handle
            % - Parameter z_in: depth grid or domain bounds associated with `rho`
            % - Returns N0: estimated constant buoyancy frequency
            % - Returns rho0: estimated reference surface density
            g = 9.81;
            if isa(rho,'function_handle') == true || isa(rho,'BSpline') == true
                rho0 = rho(max(z_in));
                drhodz = (rho(max(z_in)) - rho(min(z_in)))/( max(z_in) - min(z_in) );
            elseif isa(rho,'numeric') == true
                rho0 = min(rho);
                drhodz = (min(rho)-max(rho))/(max(z_in)-min(z_in));
            end
            N0 = sqrt(-g*drhodz/rho0);
        end
        
    end
end
