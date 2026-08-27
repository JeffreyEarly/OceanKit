classdef WVNarrowBandGeostrophicForcing < WVFixedAmplitudeForcing
    % Initialize and hold a narrow band of geostrophic coefficients.
    %
    % `WVNarrowBandGeostrophicForcing` selects the radial band
    % $$k_f-\Delta k/2<K_h<k_f+\Delta k/2$$ at vertical mode `j_f` and
    % holds those `A0` coefficients fixed while they continue to participate
    % in nonlinear interactions. Active `WVAdaptiveDamping` support is
    % excluded by the inherited fixed-amplitude selection behavior.
    %
    % The diagnostic `modelSpectrum` is a function of radial wavenumber
    % $$k$$ in radians per meter. It returns the legacy geostrophic energy
    % spectrum in $$\mathrm{m^{3}\,s^{-2}}$$,
    %
    % $$
    % E(k)=\begin{cases}
    % \kappa_\epsilon k_r^{-5/3-m}k^m, & k<k_r,\\
    % \kappa_\epsilon k^{-5/3}, & k_r\le k\le k_f,\\
    % \kappa_\epsilon k_f^{4/3}k^{-3}, & k>k_f,
    % \end{cases}
    % \qquad m=3/2.
    % $$
    %
    % Supplying `r` makes it authoritative and derives the effective `k_r`.
    % Otherwise `k_r` is authoritative and derives the effective `r`. Both
    % effective values are stored. The hydrostatic/stratified and barotropic
    % branches retain their distinct surface scaling.
    %
    % Construction has an intentional transform side effect unless
    % `initialPV="none"`. The default `"narrow-band"` draws the configured
    % random spectrum and assigns only the selected band to `wvt.A0`;
    % `"full-spectrum"` assigns the complete draw. Both choices consume the
    % global random stream. `"none"` leaves `wvt.A0` unchanged and fixes its
    % current values in the selected band. Restart and resolution conversion
    % use persisted selected state, reconstruct `modelSpectrum` from the
    % scalar configuration, and do not initialize the target transform or
    % consume random numbers.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformBarotropicQG([40e3,40e3],[8,8],latitude=45,shouldAntialias=false);
    % force = WVNarrowBandGeostrophicForcing(wvt,k_f=2*wvt.dk,u_rms=0.1);
    % wvt.addForcing(force);
    % E = force.modelSpectrum(wvt.kRadial);
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Declaration: classdef WVNarrowBandGeostrophicForcing < WVFixedAmplitudeForcing

    properties (SetAccess = private)
        % Effective large-scale damping rate in inverse seconds.
        %
        % When supplied, `r` determines `k_r`. When omitted, it is derived
        % from `k_r`, `u_rms`, and the transform-specific surface scaling.
        %
        % - Topic: Inspect forcing configuration
        r (1,1) double

        % Effective arrest wavenumber $$k_r$$ in radians per meter.
        %
        % Supplying `r` takes precedence over the `k_r` option.
        %
        % - Topic: Inspect forcing configuration
        k_r (1,1) double

        % Center wavenumber $$k_f$$ of the forced band in radians per meter.
        %
        % The selected band has width `wvt.kRadial(2)-wvt.kRadial(1)`.
        %
        % - Topic: Inspect forcing configuration
        k_f (1,1) double

        % Forced vertical geostrophic-mode number.
        %
        % - Topic: Inspect forcing configuration
        j_f (1,1) double

        % Target surface root-mean-square speed in meters per second.
        %
        % This sets the total energy scale used by `modelSpectrum`.
        %
        % - Topic: Inspect forcing configuration
        u_rms (1,1) double

        % Potential-vorticity initialization choice.
        %
        % `"none"` preserves current transform coefficients,
        % `"narrow-band"` initializes only the selected band, and
        % `"full-spectrum"` initializes the complete geostrophic spectrum.
        %
        % - Topic: Inspect forcing configuration
        initialPV (1,1) string
    end

    properties (Dependent, SetAccess = private)
        % Configured radial geostrophic energy-spectrum function.
        %
        % The returned handle accepts radial wavenumber in radians per meter
        % and returns spectral density in $$\mathrm{m^{3}\,s^{-2}}$$. It is
        % reconstructed from persisted scalar configuration and is not
        % itself written to NetCDF.
        %
        % - Topic: Inspect forcing configuration
        modelSpectrum (1,1) function_handle
    end

    methods
        function contract = portableImplementationContract(self)
            % Return the inherited fixed-amplitude portable contract.
            %
            % The portable runtime needs only the canonical selected `A0`
            % state. Narrow-band construction parameters remain MATLAB-side
            % provenance and do not introduce a separate numerical forcing.
            %
            % - Topic: Forcing internals
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only fixed-amplitude contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"ApIndices",self.Ap_indices,"ApValues",self.Apbar,"AmIndices",self.Am_indices,"AmValues",self.Ambar,"A0Indices",self.A0_indices,"A0Values",self.A0bar);
            contract = WVInternal.portableImplementationContract( ...
                string(class(self)),"WVNarrowBandGeostrophicForcing", ...
                "supported","",payload);
        end

        function self = WVNarrowBandGeostrophicForcing(wvt,options)
            % Create narrow-band geostrophic fixed-amplitude forcing.
            %
            % Ordinary construction optionally initializes `wvt.A0` before
            % selecting the forced band. The `A0_indices` and `A0bar` options
            % form the canonical restart/conversion path: both must be
            % supplied together, and that path never initializes `wvt.A0`
            % or consumes random numbers.
            %
            % - Topic: Create the forcing
            % - Declaration: self = WVNarrowBandGeostrophicForcing(wvt,options)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Parameter options.name: unique forcing-registry name; default `"narrow-band geostrophic forcing"`
            % - Parameter options.r: optional authoritative damping rate in inverse seconds
            % - Parameter options.k_r: arrest wavenumber in radians per meter; default `2*wvt.dk`
            % - Parameter options.k_f: forced-band center in radians per meter; default `20*wvt.dk`
            % - Parameter options.j_f: forced vertical-mode number; default `1`
            % - Parameter options.u_rms: target surface root-mean-square speed in meters per second; default `0.2`
            % - Parameter options.initialPV: `"none"`, `"narrow-band"`, or `"full-spectrum"`; default `"narrow-band"`
            % - Parameter options.A0_indices: canonical persisted selected indices; requires `A0bar`
            % - Parameter options.A0bar: canonical persisted selected values; requires `A0_indices`
            % - Returns self: configured `WVNarrowBandGeostrophicForcing`
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.name (1,1) string = "narrow-band geostrophic forcing"
                options.r (1,1) double
                options.k_r (1,1) double = 2*wvt.dk
                options.k_f (1,1) double = 20*wvt.dk
                options.j_f (1,1) double = 1
                options.u_rms (1,1) double = 0.2
                options.initialPV (1,1) string {mustBeMember(options.initialPV,["none","narrow-band","full-spectrum"])} = "narrow-band"
                options.A0_indices (:,1) uint64
                options.A0bar (:,1) double
            end

            hasA0Indices = isfield(options,"A0_indices");
            hasA0bar = isfield(options,"A0bar");
            if xor(hasA0Indices,hasA0bar)
                error("WVNarrowBandGeostrophicForcing:IncompleteCanonicalState", "A0_indices and A0bar must be supplied together for canonical restoration.")
            end
            if hasA0Indices && numel(options.A0_indices) ~= numel(options.A0bar)
                error("WVNarrowBandGeostrophicForcing:InvalidCanonicalState", "A0_indices and A0bar must contain the same number of entries.")
            end

            parentOptions.name = options.name;
            if hasA0Indices
                parentOptions.A0_indices = options.A0_indices;
                parentOptions.A0bar = options.A0bar;
            end
            parentArguments = namedargs2cell(parentOptions);
            self@WVFixedAmplitudeForcing(wvt,parentArguments{:});

            self.k_f = options.k_f;
            self.j_f = options.j_f;
            self.u_rms = options.u_rms;
            self.initialPV = options.initialPV;

            [~,sbRatio,~,magicNumber] = WVNarrowBandGeostrophicForcing.scalingFactors(wvt,self.j_f);
            if hasA0Indices
                if ~isfield(options,"r")
                    error("WVNarrowBandGeostrophicForcing:IncompleteCanonicalConfiguration", "Canonical restoration requires both effective r and k_r.")
                end
                self.r = options.r;
                self.k_r = options.k_r;
                return
            elseif isfield(options,"r")
                self.r = options.r;
                self.k_r = options.r/(magicNumber*options.u_rms);
            else
                self.r = magicNumber*sbRatio*options.u_rms*options.k_r;
                self.k_r = options.k_r;
            end

            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            MA0 = zeros(wvt.spectralMatrixSize);
            MA0(wvt.Kh > self.k_f-deltaK/2 & wvt.Kh < self.k_f+deltaK/2 & wvt.J == self.j_f) = 1;

            if self.initialPV == "narrow-band" || self.initialPV == "full-spectrum"
                modelSpectrum = self.modelSpectrum;
                [~,~,wvt.A0] = wvt.geostrophicComponent.randomAmplitudesWithSpectrum(A0Spectrum=@(k,j) modelSpectrum(k),shouldOnlyRandomizeOrientations=1);

                if self.initialPV == "narrow-band"
                    wvt.A0 = MA0.*wvt.A0;
                else
                    [surfaceMag,~,h] = WVNarrowBandGeostrophicForcing.scalingFactors(wvt,self.j_f);
                    if isa(wvt,"WVGeometryDoublyPeriodicBarotropic")
                        u = wvt.u;
                        v = wvt.v;
                    else
                        u = wvt.ssu;
                        v = wvt.ssv;
                    end
                    zeta = wvt.ssh;
                    KE = mean(mean(0.5*(u.^2+v.^2)));
                    PE = mean(mean(0.5*(9.81*zeta.^2)/h));
                    u_rms_surface = mean(mean(sqrt(u.^2+v.^2)));
                    fprintf("surface u_rms: %.2g cm s-1\n",100*u_rms_surface);
                    fprintf("surface energy, %g.\n",KE+PE);
                    fprintf('desired energy: %g, actual energy %g\n',0.5*(surfaceMag*self.u_rms)^2,wvt.geostrophicEnergy/h);
                end
            end
            self.setGeostrophicForcingCoefficients(MA0.*wvt.A0,MA0=MA0);
        end

        function modelSpectrum = get.modelSpectrum(self)
            [surfaceMag,~,~,~] = WVNarrowBandGeostrophicForcing.scalingFactors(self.wvt,self.j_f);
            surfaceU_rms = surfaceMag*self.u_rms;
            m = 3/2;
            kappa_epsilon = 0.5*surfaceU_rms^2/(((3*m+5)/(2*m+2))*self.k_r^(-2/3)-self.k_f^(-2/3));
            model_viscous = @(k) kappa_epsilon*self.k_r^(-5/3-m)*k.^m;
            model_inverse = @(k) kappa_epsilon*k.^(-5/3);
            model_forward = @(k) kappa_epsilon*self.k_f^(4/3)*k.^(-3);
            modelSpectrum = @(k) model_viscous(k).*(k<self.k_r)+model_inverse(k).*(k>=self.k_r & k<=self.k_f)+model_forward(k).*(k>self.k_f);
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            % Convert selected state without initializing the target transform.
            %
            % The inherited fixed-amplitude mapping preserves common
            % spectral coefficients and truncates selections unavailable at
            % the target resolution. Scalar configuration is unchanged.
            %
            % - Topic: Convert forcing resolution
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: converted `WVNarrowBandGeostrophicForcing`
            [A0_indices,A0bar] = self.convertSelectedCoefficients(wvtX2,self.A0_indices,self.A0bar);
            force = WVNarrowBandGeostrophicForcing(wvtX2,name=string(self.name),r=self.r,k_r=self.k_r,k_f=self.k_f,j_f=self.j_f,u_rms=self.u_rms,initialPV=self.initialPV,A0_indices=uint64(A0_indices),A0bar=A0bar);
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            % Persist effective configuration and canonical selected state.
            %
            % `modelSpectrum` is omitted and reconstructed from the saved
            % scalar configuration during restoration.
            %
            % - Topic: Forcing persistence
            % - Declaration: writeToGroup(group,propertyAnnotations,attributes)
            % - Parameter group: transform-owned NetCDF forcing group
            % - Parameter propertyAnnotations: forcing properties to persist
            % - Parameter attributes: additional NetCDF attributes
            arguments (Input)
                self WVNarrowBandGeostrophicForcing {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end
            if isempty(propertyAnnotations)
                propertyAnnotations = self.propertyAnnotationWithName(self.requiredProperties());
            end
            writeToGroup@CAAnnotatedClass(self,group,propertyAnnotations,attributes);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            % Return configuration and selected state required for restart.
            %
            % - Topic: Forcing persistence
            % - Declaration: vars = classRequiredPropertyNames()
            % - Returns vars: required persisted property names
            vars = {'name','r','k_r','k_f','j_f','u_rms','initialPV','A0_indices','A0bar'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Return annotated NetCDF metadata for the forcing.
            %
            % - Topic: Forcing persistence
            % - Declaration: propertyAnnotations = classDefinedPropertyAnnotations()
            % - Returns propertyAnnotations: effective configuration and selected-state annotations
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            fixedAnnotations = WVFixedAmplitudeForcing.classDefinedPropertyAnnotations();
            fixedNames = string({fixedAnnotations.name});
            propertyAnnotations = fixedAnnotations(ismember(fixedNames,["name","A0_indices","A0bar"]));
            propertyAnnotations(end+1) = CANumericProperty('r',{},'s-1','effective large-scale damping rate');
            propertyAnnotations(end+1) = CANumericProperty('k_r',{},'rad m-1','effective arrest wavenumber');
            propertyAnnotations(end+1) = CANumericProperty('k_f',{},'rad m-1','forced-band center wavenumber');
            propertyAnnotations(end+1) = CANumericProperty('j_f',{},'1','forced vertical geostrophic-mode number');
            propertyAnnotations(end+1) = CANumericProperty('u_rms',{},'m s-1','target surface root-mean-square speed');
            propertyAnnotations(end+1) = CAPropertyAnnotation('initialPV','potential-vorticity initialization choice');
        end
    end

    methods (Static, Access = private)
        function [surfaceMag,sbRatio,h,magicNumber] = scalingFactors(wvt,j_f)
            if ~isa(wvt,"WVGeometryDoublyPeriodicBarotropic")
                F = wvt.FinvMatrix;
                surfaceMag = 1/F(end,j_f+1);
                sbRatio = abs(F(end,j_f+1)/F(1,j_f+1));
                h = wvt.h_0(j_f+1);
                magicNumber = 2.25;
            else
                surfaceMag = 1;
                sbRatio = 1;
                h = wvt.h;
                magicNumber = 0.0225;
            end
        end
    end
end
