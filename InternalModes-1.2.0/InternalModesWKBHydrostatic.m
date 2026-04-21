classdef InternalModesWKBHydrostatic < InternalModesSpectral
    % Compute hydrostatic WKB mode approximations from a spectrally resolved stratification.
    %
    % `InternalModesWKBHydrostatic` uses the spectral initialization path of
    % [`InternalModesSpectral`](/internal-modes/classes/numerical-solvers/internalmodesspectral)
    % to resolve $$N^2(z)$$, then applies a hydrostatic WKB approximation for
    % the fixed-$$\omega$$ problem. This class is useful as an asymptotic
    % comparison tool rather than as the primary production solver.
    %
    % In this approximation, the vertical phase coordinate is built from
    % the positive part of $$N(z) - \omega$$, and the modal depth is
    % estimated from
    %
    % $$
    % h_j = \frac{1}{g}\left(\frac{d}{j\pi}\right)^2,
    % $$
    %
    % where $$d$$ is the accumulated hydrostatic WKB phase over the
    % oscillatory region.
    %
    % ```matlab
    % im = InternalModesWKBHydrostatic(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude);
    % [F, G, h, k] = im.modesAtFrequency(5*im.f0);
    % ```
    %
    % - Topic: Create and initialize modes
    % - Topic: Compute modes
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesWKBHydrostatic < InternalModesSpectral
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKBHydrostatic(options)
            % Initialize the hydrostatic WKB approximation.
            %
            % - Topic: Create and initialize modes
            % - Declaration: im = InternalModesWKBHydrostatic(options)
            % - Parameter options.rho: density profile as gridded values, a spline, or a function handle
            % - Parameter options.N2: buoyancy-frequency function handle used instead of `rho`
            % - Parameter options.zIn: input depth grid or domain bounds
            % - Parameter options.zOut: output depth grid
            % - Parameter options.latitude: latitude in degrees
            % - Parameter options.rho0: reference surface density
            % - Parameter options.nModes: optional cap on the number of modes returned
            % - Parameter options.nEVP: spectral initialization resolution
            % - Parameter options.rotationRate: planetary rotation rate in radians per second
            % - Parameter options.g: gravitational acceleration
            % - Returns im: hydrostatic WKB solver instance
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
        function [F,G,h,omega] = modesAtWavenumber(self, k )
            % Report that the hydrostatic WKB solver does not implement fixed-$$K$$ modes.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [F,G,h,omega] = modesAtWavenumber(self,k)
            % - Parameter self: InternalModesWKBHydrostatic instance
            % - Parameter k: horizontal wavenumber
            % - Returns F: not returned because this method throws
            % - Returns G: not returned because this method throws
            % - Returns h: not returned because this method throws
            % - Returns omega: not returned because this method throws
            error('Not yet implemented');
        end
        
        function [F,G,h,k] = modesAtFrequency(self, omega )
            % Return hydrostatic WKB modes at a fixed frequency.
            %
            % - Topic: Compute modes
            % - Declaration: [F,G,h,k] = modesAtFrequency(self,omega)
            % - Parameter self: InternalModesWKBHydrostatic instance
            % - Parameter omega: frequency in radians per second
            % - Returns F: horizontal-velocity mode matrix on `zOut`
            % - Returns G: vertical-velocity mode matrix on `zOut`
            % - Returns h: equivalent-depth row vector
            % - Returns k: horizontal wavenumber row vector implied by `h` and `omega`
            N = flip(sqrt(self.N2_xLobatto));
            z = flip(self.xLobatto);
            
            Nzeroed = N-omega;
            N(Nzeroed <= 0) = 0;
            
            xi = cumtrapz(z,N);
            xi_out = interp1(z,xi,self.z);
            d = xi(end);
            
            Nout = sqrt(self.N2);
            
            j = 1:self.nEVP;
            g = 9.81;
            G = sqrt(2*g/d) * (sin(xi_out*j*pi/d)./sqrt(Nout));
            h = ((1/g)*(d./(j*pi)).^2);
            
            zeroMask = xi_out > 0;
            F = sqrt(2*g/d) * (zeroMask .* h .* sqrt(Nout) .* j*pi/d .* cos(xi_out*j*pi/d));
            
            % Grab sign of F at the ocean surface
            Fsign = sign(h .* sqrt(N(end)) .* j*pi/d .* cos(xi(end)*j*pi/d));
            
            F = Fsign .* F;
            G = Fsign .* G;
            
            switch self.normalization
                case Normalization.kConstant
                    %no-op
                case Normalization.omegaConstant
                    A = sqrt( self.Lz ./ h);
                    F = A.*F;
                    G = G.*G;
                otherwise
                    error('This normalization is not available for the analytical WKB solution.');
            end
            
            if self.upperBoundary == UpperBoundary.freeSurface
                error('No analytical free surface solution.');
            end
            
            h = h';
            k = self.kFromOmega(h,omega);
        end
        
    end
        
end
