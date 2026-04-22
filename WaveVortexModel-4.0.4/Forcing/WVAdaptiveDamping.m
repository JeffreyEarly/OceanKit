classdef WVAdaptiveDamping < WVForcing
    % Adaptive small-scale damping
    %
    % This damping operator is a linear closure that dynamically changes
    % its amplitude to keep the Reynolds number at the grid scale equal to
    % one. This closure is ideal for a spin-up problem where the amplitude
    % of the flow is changing.
    % 
    % This closure has a number of noteworthy features:
    %
    % - It does not mix geostrophic and wave modes, which requires setting
    % the diffusivity equal to the viscosity.
    % - The properties `k_no_damp` and `j_no_damp` indicate the wavenumber and
    % mode below which there is zero damping, due to the spectral vanishing
    % viscosity filter.
    % - The properties `k_damp` and `j_damp` are *estimates* of the
    % wavenumber and mode above which significant damping will occur.
    %
    % The damping operator acts in the spectral domain, directly damping
    % the wave-vortex coefficients.
    %
    % $$
    % \begin{align}
    %     \partial_t A_\pm^{k\ell j} =& - \nu (k^2 + \ell^2 ) A_\pm^{k\ell j} - \nu_z \lambda_j^{-2} A_\pm^{k\ell j} \\
    %     \partial_t A_0^{k\ell j} =& - \nu (k^2 + \ell^2 ) A_0^{k\ell j} - \nu_z \lambda_j^{-2} A_0^{k\ell j}
    % \end{align}
    % $$
    %
    % where
    %
    % $$
    % \nu_z = \nu \lambda^2_\textrm{min} k^2_\textrm{max} = \nu \lambda^2_\textrm{min} \left( \frac{\pi}{\Delta} \right)^2
    % $$
    %
    % is chosen to make the damping isotropic. The notation here is that
    % $$\Delta$$ is the horizontal grid resolution and
    % $$\lambda^2_\textrm{min}$$ is the smallest resolved radius of
    % deformation. The value of $$\nu$$ is set as
    %
    % $$
    % \nu = \frac{U \Delta}{\pi^2}
    % $$
    %
    % where $$U$$ is the maximum fluid velocity.
    %
    % ### Usage
    %
    % Assuming there is a WVTransform instance wvt, to add this forcing,
    %
    % ```matlab
    % wvt.addForcing(WVAdaptiveDamping(wvt));
    % ```
    %
    % ### Notes
    %
    % This currently damps the non-hydrostatic wavemodes the same as the
    % hydrostatic geostrophic modes. The non-hydrostatic modes would have a
    % smaller deformation radius, and thus would be damped more strongly.
    % So arguably they're under-damped in a non-hydrostatic simulation.
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: Internal
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVAdaptiveDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % spectral matrix that multiplies Ap,Am,A0 to damp
        %
        % - Topic: Properties
        damp

        % wavenumber at which the significant scale damping starts.
        %
        % - Topic: Properties
        k_damp

        % wavenumber below which there is zero damping
        %
        % - Topic: Properties
        k_no_damp
        
        % wavenumber at which the significant scale damping starts.
        %
        % - Topic: Properties
        j_damp

        % wavenumber below which there is zero damping
        %
        % - Topic: Properties
        j_no_damp

        % effective resolution used in the damping calculation
        %
        % - Topic: Properties
        assumedEffectiveHorizontalGridResolution = Inf;
    end

    properties (Access = private, Hidden)
        forcingListener
    end

    methods (Access = private, Hidden)
        function forcingDidChangeNotification(self,~,~)
            if self.wvt.effectiveHorizontalGridResolution ~= self.assumedEffectiveHorizontalGridResolution
                self.buildDampingOperator();
            end
        end
    end

    methods
        function self = WVAdaptiveDamping(wvt)
            % initialize the WVAdaptiveDamping
            %
            % - Topic: Initialization
            % - Declaration: damp = WVAdaptiveDamping(wvt)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVAdaptiveDamping instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"adaptive damping",WVForcingType(["Spectral","PVSpectral"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.buildDampingOperator();
            self.forcingListener = addlistener(self.wvt,'forcingDidChange',@self.forcingDidChangeNotification);
        end

        function didGetRemovedFromTransform(self, wvt)
            delete(self.forcingListener);
            self.forcingListener = [];
        end

        function buildDampingOperator(self)
            % Builds the damping operator
            %
            % - Topic: Internal
            % - Declaration: buildDampingOperator(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: None
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
            end
            self.assumedEffectiveHorizontalGridResolution = self.wvt.effectiveHorizontalGridResolution;

            kl_max = pi/self.assumedEffectiveHorizontalGridResolution;
            j_max = self.wvt.effectiveJMax;
            j_index = find(self.wvt.j == self.wvt.effectiveJMax);
            [K,L,~] = self.wvt.kljGrid;
            [Qkl,Qj,self.k_no_damp,self.k_damp,self.j_no_damp,self.j_damp] = self.spectralVanishingViscosityFilter(kl_max, j_max);
            prefactor_xy = self.assumedEffectiveHorizontalGridResolution/(pi^2);
            prefactor_z = (pi*pi*self.wvt.Lr2(j_index)/(self.assumedEffectiveHorizontalGridResolution)^2)*prefactor_xy;

            Lr2inv = 1./self.wvt.Lr2;
            self.damp = -prefactor_xy*Qkl.*(K.^2 +L.^2) ;
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.damp = self.damp - prefactor_z*Qj.*Lr2inv;
            end
        end

        function [Qkl,Qj,kl_cutoff, kl_damp, j_cutoff, j_damp] = spectralVanishingViscosityFilter(self, kl_max, j_max)
            % Builds the spectral vanishing viscosity operator
            %
            % - Topic: Internal
            % - Declaration: spectralVanishingViscosityFilter(self, options)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter options: struct with field shouldAssumeAntialiasing
            % - Returns: Qkl, Qj, kl_cutoff, kl_damp
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                kl_max
                j_max
            end
            wvt_ = self.wvt;
            dkl_min = min(wvt_.dk, wvt_.dl);
            kl_cutoff = dkl_min*(kl_max/dkl_min)^(3/4);

            b = sqrt(-log(0.1));
            kl_damp = (kl_max+b*kl_cutoff)/(1+b); % approximately

            [K,L,J] = wvt_.kljGrid;
            Kh = sqrt(K.^2 + L.^2);

            Qkl = exp( - ((abs(Kh)-kl_max)./(abs(Kh)-kl_cutoff)).^2 );
            Qkl(abs(Kh)<kl_cutoff) = 0;
            Qkl(abs(Kh)>kl_max) = 1;

            if wvt_.Nj > 2
                dj = wvt_.j(2)-wvt_.j(1);
                j_cutoff = dj*(j_max/dj)^(3/4);
                j_damp = (j_max+b*j_cutoff)/(1+b); % approximately
                Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
                Qj(J<j_cutoff) = 0;
                Qj(J>j_max) = 1;
            else
                j_cutoff = 0;
                j_damp = 0;
                Qj = ones(size(J));
            end
        end

        % function [Qkl,Qj,kl_cutoff, kl_damp, j_cutoff, j_damp] = spectralVanishingViscosityFilter(self, options)
        %     % Builds the spectral vanishing viscosity operator
        %     %
        %     % - Declaration: spectralVanishingViscosityFilter(self, options)
        %     % - Parameter self: an instance of WVAdaptiveViscosity
        %     % - Parameter options: struct with field shouldAssumeAntialiasing
        %     % - Returns: Qkl, Qj, kl_cutoff, kl_damp
        %     arguments
        %         self WVAdaptiveDamping {mustBeNonempty}
        %         options.shouldAssumeAntialiasing logical = false
        %     end
        %     wvt_ = self.wvt;
        %     k_max = max(wvt_.k);
        %     l_max = max(wvt_.l);
        %     j_max = max(wvt_.j);
        %     if options.shouldAssumeAntialiasing == 1
        %         k_max = 2*k_max/3;
        %         l_max = 2*l_max/3;
        %         j_max = 2*j_max/3;
        %     end
        % 
        %     kl_max = min(k_max,l_max);
        %     dkl_min = min(wvt_.dk, wvt_.dl);
        %     kl_cutoff = dkl_min*(kl_max/dkl_min)^(3/4);
        % 
        %     b = sqrt(-log(0.1));
        %     kl_damp = (kl_max+b*kl_cutoff)/(1+b); % approximately
        % 
        %     [K,L,J] = wvt_.kljGrid;
        %     Kh = sqrt(K.^2 + L.^2);
        % 
        %     Qkl = exp( - ((abs(Kh)-kl_max)./(abs(Kh)-kl_cutoff)).^2 );
        %     Qkl(abs(Kh)<kl_cutoff) = 0;
        %     Qkl(abs(Kh)>kl_max) = 1;
        % 
        %     if wvt_.Nj > 2
        %         dj = wvt_.j(2)-wvt_.j(1);
        %         j_cutoff = dj*(j_max/dj)^(3/4);
        %         j_damp = (j_max+b*j_cutoff)/(1+b); % approximately
        %         Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
        %         Qj(J<j_cutoff) = 0;
        %         Qj(J>j_max) = 1;
        %     else
        %         j_cutoff = 0;
        %         j_damp = 0;
        %         Qj = ones(size(J));
        %     end
        % end
        % 
        function dampingTimeScale = dampingTimeScale(self)
            % Computes the minimum damping time scale
            %
            % - Topic: Properties
            % - Declaration: dampingTimeScale(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: dampingTimeScale
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
            end
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            % Adds spectral forcing
            %
            % - Declaration: addSpectralForcing(self, wvt, Fp, Fm, F0)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvt: a WVTransform instance
            % - Parameter Fp: positive frequency forcing
            % - Parameter Fm: negative frequency forcing
            % - Parameter F0: zero frequency forcing
            % - Returns: Fp, Fm, F0
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                Fp double {mustBeNonempty}
                Fm double {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            uvMax = wvt.uvMax;
            Fp = Fp + uvMax * self.damp .* wvt.Ap;
            Fm = Fm + uvMax * self.damp .* wvt.Am;
            F0 = F0 + uvMax * self.damp .* wvt.A0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            % Adds potential vorticity spectral forcing
            %
            % - Declaration: addPotentialVorticitySpectralForcing(self, wvt, F0)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvt: a WVTransform instance
            % - Parameter F0: zero frequency forcing
            % - Returns: F0
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            F0 = F0 + wvt.uvMax * self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVAdaptiveDamping(wvtX2);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            arguments
            end
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Returns the defined property annotations for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classDefinedPropertyAnnotations()
            % - Returns: propertyAnnotations
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end
end