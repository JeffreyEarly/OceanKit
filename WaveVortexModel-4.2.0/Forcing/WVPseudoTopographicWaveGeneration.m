classdef WVPseudoTopographicWaveGeneration < WVForcing
    % Generate internal waves from prescribed barotropic flow over topography.
    %
    % `WVPseudoTopographicWaveGeneration` projects the first-order bottom
    % velocity
    %
    % $$
    % g_b=\boldsymbol U_{\mathrm{bt}}(t)\boldsymbol{\cdot}\nabla_Hh
    % $$
    %
    % onto the rigid-lid wave modes using their bottom pressure. The
    % projection is precomputed, so ordinary forcing calls add spectral
    % wave tendencies without a pressure solve or spatial transform. By
    % default, generation is projected outside the exact support of an
    % active `WVAdaptiveDamping`. Optional horizontal-wavenumber and
    % vertical-mode bounds support other closures. The incoming balanced
    % tendency is left unchanged. Select a standard constituent with
    % `darwinSymbol`, or supply a custom angular `frequency`.
    %
    % ```matlab
    % forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=h,barotropicVelocityAmplitude=[0.05; 0],darwinSymbol="M2");
    % wvt.removeAllForcing();
    % wvt.addForcing(forcing);
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Generate topography
    % - Topic: Inspect the forcing
    % - Topic: Evaluate the forcing
    % - Topic: Restart persistence
    % - Topic: CAAnnotatedClass requirement
    % - Declaration: classdef WVPseudoTopographicWaveGeneration < WVForcing

    properties (SetAccess = private)
        % Upward-positive topographic height $$h(x,y)$$ in meters.
        %
        % The field is stationary and periodic on the transform's
        % horizontal grid.
        %
        % - Topic: Inspect the forcing
        topographicHeight (:,:) double

        % Complex barotropic velocity amplitude in meters per second.
        %
        % The two entries are the zonal and meridional amplitudes in
        % $$\boldsymbol U_{\mathrm{bt}}=R(t)\operatorname{Re}
        % \{\widehat{\boldsymbol U}_{\mathrm{bt}}e^{-i\omega(t-t_0)}\}$$.
        %
        % - Topic: Inspect the forcing
        barotropicVelocityAmplitude (2,1) double

        % Coordinate for the barotropic-velocity components.
        %
        % Values 1 and 2 identify the zonal and meridional components,
        % respectively. This coordinate is used for NetCDF persistence.
        %
        % - Topic: Restart persistence
        barotropicVelocityComponent (2,1) double = [1; 2]

        % Barotropic angular frequency $$\omega$$ in radians per second.
        %
        % - Topic: Inspect the forcing
        frequency (1,1) double

        % Darwin symbol used to select the tidal frequency.
        %
        % This is empty when `frequency` was supplied directly. On restart,
        % the persisted symbol is descriptive metadata and the persisted
        % frequency remains authoritative.
        %
        % - Topic: Inspect the forcing
        darwinSymbol (1,1) string = ""

        % Duration of the half-cosine startup ramp in seconds.
        %
        % - Topic: Inspect the forcing
        rampDuration (1,1) double

        % Time at which the prescribed barotropic forcing begins, in seconds.
        %
        % - Topic: Inspect the forcing
        startTime (1,1) double

        % Whether generation avoids active adaptive damping.
        %
        % When true, modes for which an active `WVAdaptiveDamping` has a
        % nonzero spectral operator are excluded from the generated wave
        % tendency.
        %
        % - Topic: Inspect the forcing
        shouldAvoidAdaptiveDamping (1,1) logical

        % Largest radial horizontal wavenumber forced, in radians per meter.
        %
        % The default `Inf` applies no manual horizontal restriction.
        %
        % - Topic: Inspect the forcing
        maximumForcedHorizontalWavenumber (1,1) double

        % Largest vertical wave-mode index forced.
        %
        % The default `Inf` applies no manual vertical-mode restriction.
        %
        % - Topic: Inspect the forcing
        maximumForcedVerticalMode (1,1) double
    end

    properties (Access = private)
        dHdx
        dHdy
        responsePlusX
        responsePlusY
        responseMinusX
        responseMinusY
        cachedSpectralGenerationMask
        shouldRefreshSpectralGenerationMask (1,1) logical = true
        forcingListener
    end

    methods
        function self = WVPseudoTopographicWaveGeneration(wvt,options)
            % Create a prescribed bottom wave-generation forcing.
            %
            % Supply either `frequency` or `darwinSymbol`, but not both.
            % Omitting both selects M2. Supported Darwin symbols are `M2`,
            % `S2`, `N2`, `K1`, and `O1`. A zero ramp duration activates
            % the harmonic current immediately at `startTime`. The
            % transform must contain a wave component and implement
            % `waveModeVerticalStructureAtIndex`. Generation avoids active
            % adaptive damping by default. Manual bounds use radial
            % horizontal wavenumber and vertical wave-mode index.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = WVPseudoTopographicWaveGeneration(wvt,options)
            % - Parameter wvt: supported wave-bearing `WVTransform` receiving the forcing
            % - Parameter options.topographicHeight: real stationary terrain of size $$N_x\times N_y$$ in meters
            % - Parameter options.barotropicVelocityAmplitude: finite complex two-component velocity amplitude in meters per second
            % - Parameter options.frequency: custom positive angular frequency in radians per second
            % - Parameter options.darwinSymbol: astronomical tidal constituent used to select the frequency
            % - Parameter options.rampDuration: nonnegative startup-ramp duration in seconds
            % - Parameter options.startTime: finite forcing start time in seconds
            % - Parameter options.shouldAvoidAdaptiveDamping: whether to exclude modes damped by `WVAdaptiveDamping`
            % - Parameter options.maximumForcedHorizontalWavenumber: largest forced radial horizontal wavenumber in radians per meter
            % - Parameter options.maximumForcedVerticalMode: largest forced vertical wave-mode index
            % - Parameter options.name: forcing name registered with the transform
            % - Returns forcing: configured `WVPseudoTopographicWaveGeneration`
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.topographicHeight double
                options.barotropicVelocityAmplitude double
                options.frequency double
                options.darwinSymbol (1,1) string
                options.rampDuration double = 0
                options.startTime double = wvt.t
                options.shouldAvoidAdaptiveDamping (1,1) logical = true
                options.maximumForcedHorizontalWavenumber double = Inf
                options.maximumForcedVerticalMode double = Inf
                options.name (1,1) string = "pseudo-topographic wave generation"
            end

            if ~WVPseudoTopographicWaveGeneration.isSupportedTransform(wvt)
                error("WVPseudoTopographicWaveGeneration:UnsupportedTransform", "WVPseudoTopographicWaveGeneration requires a wave-bearing transform that implements waveModeVerticalStructureAtIndex.")
            end
            if ~isequal(size(options.topographicHeight),[wvt.Nx wvt.Ny])
                error("WVPseudoTopographicWaveGeneration:InvalidTopographicHeightSize", "topographicHeight must have size [%d %d], matching the transform horizontal grid.", wvt.Nx, wvt.Ny)
            end
            if ~isreal(options.topographicHeight) || any(~isfinite(options.topographicHeight),"all")
                error("WVPseudoTopographicWaveGeneration:InvalidTopographicHeight", "topographicHeight must be real and finite.")
            end
            if ~isequal(size(options.barotropicVelocityAmplitude),[2 1]) || any(~isfinite(options.barotropicVelocityAmplitude),"all")
                error("WVPseudoTopographicWaveGeneration:InvalidBarotropicVelocityAmplitude", "barotropicVelocityAmplitude must be a finite complex 2-by-1 vector.")
            end
            hasFrequency = isfield(options,"frequency");
            hasDarwinSymbol = isfield(options,"darwinSymbol");
            if hasFrequency && hasDarwinSymbol
                error("WVPseudoTopographicWaveGeneration:ConflictingFrequencyOptions", "Specify either frequency or darwinSymbol, but not both.")
            elseif hasFrequency
                frequency = options.frequency;
                darwinSymbol = "";
            else
                if hasDarwinSymbol
                    darwinSymbol = options.darwinSymbol;
                else
                    darwinSymbol = "M2";
                end
                frequency = WVPseudoTopographicWaveGeneration.frequencyForDarwinSymbol(darwinSymbol);
            end
            if ~isscalar(frequency) || ~isreal(frequency) || ~isfinite(frequency) || frequency <= 0
                error("WVPseudoTopographicWaveGeneration:InvalidFrequency", "frequency must be a finite positive real scalar.")
            end
            if ~isscalar(options.rampDuration) || ~isreal(options.rampDuration) || ~isfinite(options.rampDuration) || options.rampDuration < 0
                error("WVPseudoTopographicWaveGeneration:InvalidRampDuration", "rampDuration must be a finite nonnegative real scalar.")
            end
            if ~isscalar(options.startTime) || ~isreal(options.startTime) || ~isfinite(options.startTime)
                error("WVPseudoTopographicWaveGeneration:InvalidStartTime", "startTime must be a finite real scalar.")
            end
            if ~isscalar(options.maximumForcedHorizontalWavenumber) || ~isreal(options.maximumForcedHorizontalWavenumber) || isnan(options.maximumForcedHorizontalWavenumber) || options.maximumForcedHorizontalWavenumber < 0
                error("WVPseudoTopographicWaveGeneration:InvalidMaximumForcedHorizontalWavenumber", "maximumForcedHorizontalWavenumber must be a nonnegative real scalar or Inf.")
            end
            if ~isscalar(options.maximumForcedVerticalMode) || ~isreal(options.maximumForcedVerticalMode) || isnan(options.maximumForcedVerticalMode) || options.maximumForcedVerticalMode < 0
                error("WVPseudoTopographicWaveGeneration:InvalidMaximumForcedVerticalMode", "maximumForcedVerticalMode must be a nonnegative real scalar or Inf.")
            end
            if strlength(options.name) == 0
                error("WVPseudoTopographicWaveGeneration:InvalidName", "name must be a nonempty string.")
            end

            self@WVForcing(wvt,options.name,WVForcingType("Spectral"));
            self.topographicHeight = options.topographicHeight;
            self.barotropicVelocityAmplitude = options.barotropicVelocityAmplitude;
            self.frequency = frequency;
            self.darwinSymbol = darwinSymbol;
            self.rampDuration = options.rampDuration;
            self.startTime = options.startTime;
            self.shouldAvoidAdaptiveDamping = options.shouldAvoidAdaptiveDamping;
            self.maximumForcedHorizontalWavenumber = options.maximumForcedHorizontalWavenumber;
            self.maximumForcedVerticalMode = options.maximumForcedVerticalMode;
            self.dHdx = wvt.diffX(self.topographicHeight);
            self.dHdy = wvt.diffY(self.topographicHeight);

            terrainFourier = wvt.transformFromSpatialDomainWithFourier(repmat(self.topographicHeight,1,1,wvt.Nz));
            [self.responsePlusX,self.responsePlusY,self.responseMinusX,self.responseMinusY] = self.buildResponses(wvt,terrainFourier(1,:));
            self.forcingListener = addlistener(self.wvt,'forcingDidChange',@self.forcingDidChangeNotification);
        end

        function velocity = barotropicVelocityAtTime(self,t)
            % Evaluate the prescribed horizontally uniform current.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: velocity = barotropicVelocityAtTime(t)
            % - Parameter t: finite scalar time in seconds
            % - Returns velocity: real two-component velocity in meters per second
            arguments (Input)
                self WVPseudoTopographicWaveGeneration
                t (1,1) double {mustBeFinite}
            end
            arguments (Output)
                velocity (2,1) double
            end

            elapsed = t-self.startTime;
            if elapsed < 0
                velocity = zeros(2,1);
                return
            end
            if self.rampDuration == 0 || elapsed >= self.rampDuration
                ramp = 1;
            else
                ramp = 0.5*(1-cos(pi*elapsed/self.rampDuration));
            end
            velocity = ramp*real(self.barotropicVelocityAmplitude*exp(-1i*self.frequency*elapsed));
        end

        function gBottom = bottomVelocityAtTime(self,t)
            % Evaluate $$g_b=\boldsymbol U_{\mathrm{bt}}\boldsymbol{\cdot}\nabla_Hh$$.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: gBottom = bottomVelocityAtTime(t)
            % - Parameter t: finite scalar time in seconds
            % - Returns gBottom: real bottom-normal velocity on the horizontal grid
            velocity = self.barotropicVelocityAtTime(t);
            gBottom = velocity(1)*self.dHdx+velocity(2)*self.dHdy;
        end

        function [Fp,Fm,F0] = addSpectralForcing(self,wvt,Fp,Fm,F0)
            % Add the precomputed wave-generation tendency.
            %
            % The physical modal tendencies are converted componentwise to
            % WaveVortexModel's stored interaction representation. `F0` is
            % returned without modification.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [Fp,Fm,F0] = addSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: transform at the current model time
            % - Parameter Fp: accumulated positive-wave tendency
            % - Parameter Fm: accumulated negative-wave tendency
            % - Parameter F0: accumulated balanced tendency
            % - Returns Fp: positive-wave tendency including this forcing
            % - Returns Fm: negative-wave tendency including this forcing
            % - Returns F0: unchanged incoming balanced tendency
            self.requireOriginatingTransform(wvt);
            velocity = self.barotropicVelocityAtTime(wvt.t);
            generationMask = self.spectralGenerationMask();
            Fpt = generationMask.*(velocity(1)*self.responsePlusX+velocity(2)*self.responsePlusY);
            Fmt = generationMask.*(velocity(1)*self.responseMinusX+velocity(2)*self.responseMinusY);
            Fp = Fp+Fpt.*wvt.conjPhase;
            Fm = Fm+Fmt.*wvt.phase;
        end

        function [mask,components] = spectralGenerationMask(self)
            % Return the spectral region eligible for bottom-wave generation.
            %
            % The common `mask` combines wave validity, the manual radial
            % horizontal-wavenumber and vertical-mode bounds, and the exact
            % zero-damping support of active `WVAdaptiveDamping` objects.
            % `components` reports those masks separately, including the
            % distinct positive- and negative-wave validity masks.
            %
            % - Topic: Inspect the forcing
            % - Declaration: [mask,components] = spectralGenerationMask()
            % - Returns mask: logical mask applied to both generated wave tendencies
            % - Returns components: structure containing each constituent and branch-specific effective mask
            arguments (Input)
                self WVPseudoTopographicWaveGeneration {mustBeNonempty}
            end
            arguments (Output)
                mask logical
                components struct
            end

            if self.shouldRefreshSpectralGenerationMask
                if nargout > 1
                    [self.cachedSpectralGenerationMask,components] = self.buildSpectralGenerationMask();
                else
                    self.cachedSpectralGenerationMask = self.buildSpectralGenerationMask();
                end
                self.shouldRefreshSpectralGenerationMask = false;
            elseif nargout > 1
                [~,components] = self.buildSpectralGenerationMask();
            end
            mask = self.cachedSpectralGenerationMask;
        end

        function forcing = forcingWithResolutionOfTransform(self,wvtX2)
            % Rebuild the forcing for a transform at another resolution.
            %
            % The terrain is transferred spectrally, preserving common
            % Fourier coefficients while truncating or zero-padding modes
            % that are not shared by the two transforms. All modal response
            % arrays are then rebuilt for `wvtX2`.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: compatible supported transform at the target resolution
            % - Returns forcing: equivalent forcing rebuilt for `wvtX2`
            arguments (Input)
                self WVPseudoTopographicWaveGeneration {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            arguments (Output)
                forcing WVPseudoTopographicWaveGeneration
            end

            if ~WVPseudoTopographicWaveGeneration.isSupportedTransform(wvtX2)
                error("WVPseudoTopographicWaveGeneration:UnsupportedTransform", "WVPseudoTopographicWaveGeneration requires a wave-bearing transform that implements waveModeVerticalStructureAtIndex.")
            end
            if ~isequal([self.wvt.Lx self.wvt.Ly self.wvt.Lz],[wvtX2.Lx wvtX2.Ly wvtX2.Lz])
                error("WVPseudoTopographicWaveGeneration:IncompatibleDomain", "Resolution conversion requires transforms with identical domain dimensions.")
            end

            terrainFourier = self.wvt.transformFromSpatialDomainWithFourier(repmat(self.topographicHeight,1,1,self.wvt.Nz));
            terrainFourierX2 = self.wvt.spectralVariableWithResolution(wvtX2,terrainFourier);
            terrainX2 = wvtX2.transformToSpatialDomainWithFourier(repmat(terrainFourierX2(1,:),wvtX2.Nz,1));
            forcing = WVPseudoTopographicWaveGeneration(wvtX2,topographicHeight=real(terrainX2(:,:,1)),barotropicVelocityAmplitude=self.barotropicVelocityAmplitude,frequency=self.frequency,rampDuration=self.rampDuration,startTime=self.startTime,shouldAvoidAdaptiveDamping=self.shouldAvoidAdaptiveDamping,maximumForcedHorizontalWavenumber=self.maximumForcedHorizontalWavenumber,maximumForcedVerticalMode=self.maximumForcedVerticalMode,name=string(self.name));
            forcing.darwinSymbol = self.darwinSymbol;
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            % Write the forcing to its transform-owned NetCDF group.
            %
            % Transform-derived gradients and response arrays are omitted.
            % The terrain is written separately because its `x` and `y`
            % dimensions are owned by the parent transform group.
            %
            % - Topic: Restart persistence
            % - Declaration: writeToGroup(group,propertyAnnotations,attributes)
            % - Parameter group: transform-owned NetCDF group for this forcing
            % - Parameter propertyAnnotations: forcing properties to persist
            % - Parameter attributes: additional NetCDF attributes
            arguments (Input)
                self WVPseudoTopographicWaveGeneration {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end

            isTopographicHeight = string({propertyAnnotations.name}) == "topographicHeight";
            writeToGroup@CAAnnotatedClass(self,group,propertyAnnotations(~isTopographicHeight),attributes);
            if any(isTopographicHeight)
                annotation = propertyAnnotations(find(isTopographicHeight,1));
                variableAttributes = annotation.attributes;
                variableAttributes('units') = annotation.units;
                variableAttributes('long_name') = annotation.description;
                group.addVariable(annotation.name,annotation.dimensions,self.topographicHeight,isComplex=annotation.isComplex,attributes=variableAttributes);
            end
        end

        function didGetRemovedFromTransform(self,~)
            % Release the transform forcing-change listener.
            %
            % - Topic: CAAnnotatedClass requirement
            % - Developer: true
            if ~isempty(self.forcingListener)
                delete(self.forcingListener);
                self.forcingListener = [];
            end
        end
    end

    methods (Access = private)
        function forcingDidChangeNotification(self,~,~)
            self.shouldRefreshSpectralGenerationMask = true;
        end

        function [mask,components] = buildSpectralGenerationMask(self)
            positiveWaveValidity = logical(self.wvt.waveComponent.maskAp);
            negativeWaveValidity = logical(self.wvt.waveComponent.maskAm);
            waveValidity = positiveWaveValidity | negativeWaveValidity;
            horizontalBound = self.wvt.Kh <= self.maximumForcedHorizontalWavenumber;
            verticalBound = self.wvt.J <= self.maximumForcedVerticalMode;
            adaptiveDamping = true(size(waveValidity));

            if self.shouldAvoidAdaptiveDamping
                forcing = self.wvt.forcing;
                for iForcing = 1:numel(forcing)
                    if isa(forcing(iForcing),"WVAdaptiveDamping")
                        adaptiveDamping = adaptiveDamping & forcing(iForcing).damp == 0;
                    end
                end
            end

            mask = waveValidity & horizontalBound & verticalBound & adaptiveDamping;
            if nargout > 1
                components = struct( ...
                    positiveWaveValidity=positiveWaveValidity, ...
                    negativeWaveValidity=negativeWaveValidity, ...
                    waveValidity=waveValidity, ...
                    horizontalBound=horizontalBound, ...
                    verticalBound=verticalBound, ...
                    adaptiveDamping=adaptiveDamping, ...
                    effectivePositive=mask & positiveWaveValidity, ...
                    effectiveNegative=mask & negativeWaveValidity);
            end
        end

        function requireOriginatingTransform(self,wvt)
            if wvt ~= self.wvt
                error("WVPseudoTopographicWaveGeneration:TransformMismatch", "The forcing can only be evaluated with the WVTransform instance used during construction.")
            end
        end
    end

    methods (Static, Access = private)
        function frequency = frequencyForDarwinSymbol(darwinSymbol)
            switch darwinSymbol
                case "M2"
                    periodHours = 12.420602;
                case "S2"
                    periodHours = 12.000000;
                case "N2"
                    periodHours = 12.65834751;
                case "K1"
                    periodHours = 23.93447213;
                case "O1"
                    periodHours = 25.81933871;
                otherwise
                    error("WVPseudoTopographicWaveGeneration:UnknownConstituent", "Unknown tidal constituent '%s'.",darwinSymbol)
            end
            frequency = 2*pi/(periodHours*3600);
        end

        function requiredPropertyNames = persistedRequiredPropertyNames()
            requiredPropertyNames = { ...
                'topographicHeight', ...
                'barotropicVelocityAmplitude', ...
                'frequency', ...
                'darwinSymbol', ...
                'rampDuration', ...
                'startTime', ...
                'shouldAvoidAdaptiveDamping', ...
                'maximumForcedHorizontalWavenumber', ...
                'maximumForcedVerticalMode', ...
                'name'};
        end

        function [responsePlusX,responsePlusY,responseMinusX,responseMinusY] = buildResponses(wvt,terrainFourier)
            [~,iBottom] = min(wvt.z);
            bottomF = wvt.waveModeVerticalStructureAtIndex(iBottom);
            piPlus = wvt.g*bottomF.*wvt.NAp;
            piMinus = wvt.g*bottomF.*wvt.NAm;

            dHdxFourier = 1i*wvt.K.*terrainFourier;
            dHdyFourier = 1i*wvt.L.*terrainFourier;
            responsePlusX = complex(zeros(size(piPlus)));
            responsePlusY = complex(zeros(size(piPlus)));
            responseMinusX = complex(zeros(size(piMinus)));
            responseMinusY = complex(zeros(size(piMinus)));
            maskPlus = logical(wvt.waveComponent.maskAp);
            maskMinus = logical(wvt.waveComponent.maskAm);
            responsePlusX(maskPlus) = conj(piPlus(maskPlus)).*dHdxFourier(maskPlus)./wvt.Apm_TE_factor(maskPlus);
            responsePlusY(maskPlus) = conj(piPlus(maskPlus)).*dHdyFourier(maskPlus)./wvt.Apm_TE_factor(maskPlus);
            responseMinusX(maskMinus) = conj(piMinus(maskMinus)).*dHdxFourier(maskMinus)./wvt.Apm_TE_factor(maskMinus);
            responseMinusY(maskMinus) = conj(piMinus(maskMinus)).*dHdyFourier(maskMinus)./wvt.Apm_TE_factor(maskMinus);
        end

        function tf = isSupportedTransform(wvt)
            tf = wvt.hasWaveComponent && ismethod(wvt,"waveModeVerticalStructureAtIndex");
        end

        function [radialWavenumber,radialPowerSpectrum,radialTargetPowerSpectrum,radialModeCount] = radialSpectrum(Kh,realizedPowerSpectrum,targetPowerSpectrum,binWidth)
            maximumWavenumber = max(Kh,[],"all");
            radialWavenumber = (0:ceil(maximumWavenumber/binWidth)).'*binWidth;
            radialPowerSpectrum = nan(size(radialWavenumber));
            radialTargetPowerSpectrum = nan(size(radialWavenumber));
            radialModeCount = zeros(size(radialWavenumber));
            for iBin = 1:numel(radialWavenumber)
                indices = Kh >= max(0,radialWavenumber(iBin)-binWidth/2) & Kh < radialWavenumber(iBin)+binWidth/2;
                radialModeCount(iBin) = nnz(indices);
                if radialModeCount(iBin) > 0
                    radialPowerSpectrum(iBin) = mean(realizedPowerSpectrum(indices));
                    radialTargetPowerSpectrum(iBin) = mean(targetPowerSpectrum(indices));
                end
            end
        end
    end

    methods (Static)
        function requiredPropertyNames = classRequiredPropertyNames()
            % Return the forcing properties required for restart.
            %
            % Derived gradients and modal responses are intentionally not
            % persisted; construction against the restored transform
            % rebuilds them.
            %
            % - Topic: Restart persistence
            % - Declaration: requiredPropertyNames = classRequiredPropertyNames()
            % - Returns requiredPropertyNames: property names required to reconstruct the forcing
            arguments (Output)
                requiredPropertyNames cell
            end

            requiredPropertyNames = { ...
                'topographicHeight', ...
                'barotropicVelocityAmplitude', ...
                'frequency', ...
                'darwinSymbol', ...
                'rampDuration', ...
                'startTime', ...
                'shouldAvoidAdaptiveDamping', ...
                'maximumForcedHorizontalWavenumber', ...
                'maximumForcedVerticalMode', ...
                'name'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Return metadata used to persist the forcing configuration.
            %
            % - Topic: Restart persistence
            % - Declaration: propertyAnnotations = classDefinedPropertyAnnotations()
            % - Returns propertyAnnotations: annotated forcing properties and dimensions
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end

            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('barotropicVelocityComponent','1','barotropic-velocity component; 1 is zonal and 2 is meridional');
            propertyAnnotations(end+1) = CANumericProperty('topographicHeight',{'x','y'},'m','upward-positive topographic height');
            propertyAnnotations(end+1) = CANumericProperty('barotropicVelocityAmplitude',{'barotropicVelocityComponent'},'m s^{-1}','complex barotropic-velocity amplitude',isComplex=true);
            propertyAnnotations(end+1) = CANumericProperty('frequency',{},'rad s^{-1}','barotropic angular frequency');
            propertyAnnotations(end+1) = CAPropertyAnnotation('darwinSymbol','Darwin symbol used to select the frequency, or empty for a custom frequency');
            propertyAnnotations(end+1) = CANumericProperty('rampDuration',{},'s','half-cosine startup-ramp duration');
            propertyAnnotations(end+1) = CANumericProperty('startTime',{},'s','model time at which the forcing begins');
            propertyAnnotations(end+1) = CANumericProperty('shouldAvoidAdaptiveDamping',{},'bool','whether generation avoids active adaptive damping');
            propertyAnnotations(end+1) = CANumericProperty('maximumForcedHorizontalWavenumber',{},'rad m^{-1}','largest radial horizontal wavenumber receiving generation');
            propertyAnnotations(end+1) = CANumericProperty('maximumForcedVerticalMode',{},'1','largest vertical wave-mode index receiving generation');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the forcing');
        end

        function forcing = forcingFromGroup(group,wvt)
            % Reconstruct a forcing from its annotated NetCDF group.
            %
            % - Topic: Restart persistence
            % - Declaration: forcing = forcingFromGroup(group,wvt)
            % - Parameter group: NetCDF group containing the forcing state
            % - Parameter wvt: restored transform receiving the forcing
            % - Returns forcing: reconstructed `WVPseudoTopographicWaveGeneration`
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            arguments (Output)
                forcing WVPseudoTopographicWaveGeneration
            end

            requiredPropertyNames = WVPseudoTopographicWaveGeneration.persistedRequiredPropertyNames();
            missingPropertyNames = string.empty(1,0);
            for iProperty = 1:numel(requiredPropertyNames)
                propertyName = requiredPropertyNames{iProperty};
                isPresent = group.hasVariableWithName(propertyName) || group.hasGroupWithName(propertyName) || isKey(group.attributes,propertyName);
                if ~isPresent
                    missingPropertyNames(end+1) = string(propertyName); %#ok<AGROW>
                end
            end
            if ~isempty(missingPropertyNames)
                error("WVPseudoTopographicWaveGeneration:IncompleteRestart", "The restart group is missing required bottom-wave-generation properties: %s.", join(missingPropertyNames,", "))
            end

            options = CAAnnotatedClass.propertyValuesFromGroup(group,requiredPropertyNames);
            options.barotropicVelocityAmplitude = reshape(options.barotropicVelocityAmplitude,2,1);
            options.shouldAvoidAdaptiveDamping = logical(options.shouldAvoidAdaptiveDamping);
            options.name = string(options.name);
            persistedDarwinSymbol = string(options.darwinSymbol);
            if strlength(persistedDarwinSymbol) > 0
                WVPseudoTopographicWaveGeneration.frequencyForDarwinSymbol(persistedDarwinSymbol);
            end
            options = rmfield(options,"darwinSymbol");
            optionArguments = namedargs2cell(options);
            forcing = WVPseudoTopographicWaveGeneration(wvt,optionArguments{:});
            forcing.darwinSymbol = persistedDarwinSymbol;
        end

        function [virtualDepth,topographicHeight,diagnostics] = goffAbyssalHillTopography(wvt,options)
            % Generate periodic Goff abyssal-hill topography.
            %
            % The isotropic two-dimensional power spectrum is
            %
            % $$P_h(K)=\frac{4\pi h_{\mathrm{rms}}^2}{K_c^2}
            % \left(1+\frac{K^2}{K_c^2}\right)^{-2}.$$
            %
            % The returned height is positive upward, zero mean, and
            % normalized to the requested post-filter RMS. Random phases
            % use a local stream and do not alter MATLAB's global random
            % state.
            %
            % ```matlab
            % [H,h,diagnostics] = WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,minimumWavelength=30e3);
            % ```
            %
            % - Topic: Generate topography
            % - Declaration: [virtualDepth,topographicHeight,diagnostics] = WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,options)
            % - Parameter wvt: periodic `WVTransform` defining the horizontal grid and mean depth
            % - Parameter options.rmsHeight: requested post-filter RMS topographic height in meters
            % - Parameter options.cornerWavenumber: Goff corner wavenumber in radians per meter
            % - Parameter options.minimumWavelength: shortest retained wavelength in meters
            % - Parameter options.randomSeed: nonnegative seed for a local Mersenne Twister stream
            % - Returns virtualDepth: positive water depth of size $$N_x\times N_y$$ in meters
            % - Returns topographicHeight: upward-positive zero-mean height of size $$N_x\times N_y$$ in meters
            % - Returns diagnostics: realized statistics and spectral diagnostics
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.rmsHeight (1,1) double {mustBePositive,mustBeFinite} = 100
                options.cornerWavenumber (1,1) double {mustBePositive,mustBeFinite} = 1e-4
                options.minimumWavelength (1,1) double {mustBePositive,mustBeFinite} = 15e3
                options.randomSeed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 2023
            end
            arguments (Output)
                virtualDepth (:,:) double
                topographicHeight (:,:) double
                diagnostics (1,1) struct
            end

            if options.randomSeed > double(intmax("uint32"))
                error("WVPseudoTopographicWaveGeneration:InvalidRandomSeed", "randomSeed must not exceed intmax('uint32').")
            end

            cutoffWavenumber = 2*pi/options.minimumWavelength;
            retainedAxialWavenumber = pi/wvt.effectiveHorizontalGridResolution();
            if cutoffWavenumber > retainedAxialWavenumber*(1+10*eps(retainedAxialWavenumber))
                minimumResolvedWavelength = 2*pi/retainedAxialWavenumber;
                error("WVPseudoTopographicWaveGeneration:UnresolvedCutoff", "minimumWavelength must be at least %.6g m for this transform; the requested %.6g m cutoff is unresolved.", minimumResolvedWavelength, options.minimumWavelength)
            end

            [K,L] = ndgrid(wvt.k_dft,wvt.l_dft);
            Kh = hypot(K,L);
            targetPowerSpectrum = 4*pi*options.rmsHeight^2/options.cornerWavenumber^2.*(1+(Kh/options.cornerWavenumber).^2).^(-2);
            retainedMask = Kh > 0 & Kh <= cutoffWavenumber;

            stream = RandStream("mt19937ar",Seed=options.randomSeed);
            phaseSource = fft2(randn(stream,wvt.Nx,wvt.Ny));
            phase = ones(size(phaseSource));
            nonzeroPhase = abs(phaseSource) > 0;
            phase(nonzeroPhase) = phaseSource(nonzeroPhase)./abs(phaseSource(nonzeroPhase));

            domainArea = wvt.Lx*wvt.Ly;
            normalizedFourierCoefficients = zeros(wvt.Nx,wvt.Ny);
            normalizedFourierCoefficients(retainedMask) = sqrt(targetPowerSpectrum(retainedMask)/domainArea).*phase(retainedMask);
            topographicHeight = real(ifft2(normalizedFourierCoefficients*wvt.Nx*wvt.Ny));
            topographicHeight = topographicHeight-mean(topographicHeight,"all");
            realizedRmsHeight = sqrt(mean(topographicHeight.^2,"all"));
            if realizedRmsHeight == 0
                error("WVPseudoTopographicWaveGeneration:EmptySpectrum", "The requested cutoff retains no nonzero Fourier modes.")
            end
            topographicHeight = options.rmsHeight*topographicHeight/realizedRmsHeight;

            virtualDepth = wvt.Lz-topographicHeight;
            if any(virtualDepth <= 0,"all")
                error("WVPseudoTopographicWaveGeneration:NonpositiveVirtualDepth", "The generated topography reaches or exceeds the mean model depth. Reduce rmsHeight or choose another randomSeed.")
            end

            dHdx = wvt.diffX(topographicHeight);
            dHdy = wvt.diffY(topographicHeight);
            slopeMagnitude = hypot(dHdx,dHdy);
            normalizedFourierCoefficients = fft2(topographicHeight)/(wvt.Nx*wvt.Ny);
            realizedPowerSpectrum = domainArea*abs(normalizedFourierCoefficients).^2;
            powerSpectrumScale = median(realizedPowerSpectrum(retainedMask)./targetPowerSpectrum(retainedMask));
            normalizedTargetPowerSpectrum = powerSpectrumScale*targetPowerSpectrum;
            [radialWavenumber,radialPowerSpectrum,radialTargetPowerSpectrum,radialModeCount] = WVPseudoTopographicWaveGeneration.radialSpectrum(Kh,realizedPowerSpectrum,normalizedTargetPowerSpectrum,min(wvt.dk,wvt.dl));

            diagnostics = struct( ...
                meanHeight=mean(topographicHeight,"all"), ...
                rmsHeight=sqrt(mean(topographicHeight.^2,"all")), ...
                rmsSlope=sqrt(mean(slopeMagnitude.^2,"all")), ...
                maximumSlope=max(slopeMagnitude,[],"all"), ...
                cutoffWavenumber=cutoffWavenumber, ...
                cornerWavenumber=options.cornerWavenumber, ...
                minimumWavelength=options.minimumWavelength, ...
                randomSeed=options.randomSeed, ...
                retainedAxialWavenumber=retainedAxialWavenumber, ...
                powerSpectrumScale=powerSpectrumScale, ...
                fourierCoefficients=normalizedFourierCoefficients, ...
                radialWavenumber=radialWavenumber, ...
                radialPowerSpectrum=radialPowerSpectrum, ...
                radialTargetPowerSpectrum=radialTargetPowerSpectrum, ...
                radialModeCount=radialModeCount);
        end
    end
end
