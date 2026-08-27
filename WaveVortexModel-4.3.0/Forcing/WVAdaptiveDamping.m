classdef WVAdaptiveDamping < WVForcing
    % Adapt small-scale spectral damping to the current flow speed.
    %
    % This closure rebuilds its spectral shape when the transform's effective
    % resolution changes and scales its coefficient tendency by the current
    % maximum horizontal speed. It is useful when the flow amplitude evolves
    % substantially, such as during spin-up.
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
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
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
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Inspect forcing or damping scales
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    %
    % - Declaration: WVAdaptiveDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Unit-speed spectral damping operator in inverse meters.
        %
        % This array has `wvt.spectralMatrixSize`. The actual coefficient
        % damping rate is `wvt.uvMax*damp` in inverse seconds.
        %
        % - Topic: Properties
        damp

        % Estimated horizontal wavenumber for significant damping.
        %
        % Units are radians per meter. The filter is already nonzero below
        % this estimate; use `k_no_damp` for the exact zero-damping cutoff.
        %
        % - Topic: Properties
        k_damp

        % Horizontal wavenumber below which damping is exactly zero.
        %
        % Units are radians per meter.
        %
        % - Topic: Properties
        k_no_damp
        
        % Estimated vertical mode number for significant damping.
        %
        % This value is dimensionless. The filter is already nonzero below
        % this estimate; use `j_no_damp` for the exact zero-damping cutoff.
        %
        % - Topic: Properties
        j_damp

        % Vertical mode number below which damping is exactly zero.
        %
        % This value is dimensionless.
        %
        % - Topic: Properties
        j_no_damp

        % Effective horizontal resolution used to construct `damp`, in meters.
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
        function contract = portableImplementationContract(self)
            % Return the paired portable implementation contract.
            %
            % - Topic: Forcing internals
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only forcing contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"assumedEffectiveHorizontalGridResolution",double(self.assumedEffectiveHorizontalGridResolution),"kNoDamp",double(self.k_no_damp),"jNoDamp",double(self.j_no_damp));
            contract = self.supportedPortableImplementationContract("WVAdaptiveDamping",payload);
        end

        function self = WVAdaptiveDamping(wvt)
            % Create adaptive spectral damping for a transform.
            %
            % - Topic: Initialization
            % - Declaration: self = WVAdaptiveDamping(wvt)
            % - Parameter wvt: transform that owns and evaluates the closure
            % - Returns self: adaptive-damping closure owned by `wvt`
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
            % Build the unit-speed spectral damping operator.
            %
            % - Topic: Internal
            % - Declaration: buildDampingOperator(self)
            % - Parameter self: adaptive-damping instance to update
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
            % Build horizontal and vertical spectral-vanishing filters.
            %
            % - Topic: Internal
            % - Declaration: [Qkl,Qj,kl_cutoff,kl_damp,j_cutoff,j_damp] = spectralVanishingViscosityFilter(kl_max,j_max)
            % - Parameter kl_max: maximum resolved horizontal wavenumber in radians per meter
            % - Parameter j_max: maximum resolved vertical-mode number
            % - Returns Qkl: horizontal filter on the spectral grid
            % - Returns Qj: vertical filter on the spectral grid
            % - Returns kl_cutoff: exact horizontal zero-damping cutoff in radians per meter
            % - Returns kl_damp: estimated horizontal significant-damping wavenumber in radians per meter
            % - Returns j_cutoff: exact vertical zero-damping cutoff
            % - Returns j_damp: estimated vertical significant-damping mode
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
        %     % - Parameter self: an instance of WVAdaptiveDamping
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
            % Return the inverse maximum unit-speed damping coefficient.
            %
            % Despite the historical method name, this value has units of
            % meters because `damp` has units of inverse meters. For a
            % nonzero flow, divide this value by `wvt.uvMax` to obtain the
            % shortest instantaneous e-folding time in seconds.
            %
            % - Topic: Properties
            % - Declaration: dampingTimeScale = dampingTimeScale()
            % - Returns dampingTimeScale: inverse maximum absolute entry of `damp`, in meters
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
            end
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            % Add adaptive damping to wave-vortex coefficient tendencies.
            %
            % - Declaration: addSpectralForcing(self, wvt, Fp, Fm, F0)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Fp: accumulated `Ap` tendency
            % - Parameter Fm: accumulated `Am` tendency
            % - Parameter F0: accumulated `A0` tendency
            % - Returns Fp: damped `Ap` tendency
            % - Returns Fm: damped `Am` tendency
            % - Returns F0: damped `A0` tendency
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
            % Add adaptive damping to the QGPV coefficient tendency.
            %
            % - Declaration: addPotentialVorticitySpectralForcing(self, wvt, F0)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter F0: accumulated `A0` tendency
            % - Returns F0: damped `A0` tendency
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            F0 = F0 + wvt.uvMax * self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Create equivalent adaptive damping for another resolution.
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: adaptive damping owned by `wvtX2`
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
