classdef WVTransform < matlab.mixin.indexing.RedefinesDot & CAAnnotatedClass
    % Represent a fluid state with orthogonal wave and geostrophic solutions.
    %
    % `WVTransform` is the abstract base class for wave-vortex transforms.
    % Each concrete transform represents the fluid state at time
    % `t` with spectral coefficients and exposes physical variables such as
    % velocity, isopycnal displacement, pressure, and potential vorticity.
    % Density is denoted by $$\rho$$. Energetic quantities are normalized per
    % unit reference density $$\rho_0$$ unless a method states otherwise.
    %
    % The distinguishing feature of a `WVTransform` is that an instantaneous
    % fluid state is represented as energetically orthogonal wave, inertial,
    % geostrophic, and mean-density-anomaly constituents. No temporal filter is
    % required to perform the decomposition. The same object can reconstruct
    % $$(u,v,w,\rho,p)$$, relative vorticity, potential vorticity, energetic
    % diagnostics, and custom registered variables from those coefficients.
    %
    % Choose one of five concrete transform classes:
    %
    % + `WVTransformConstantStratification`, in hydrostatic or nonhydrostatic mode
    % + `WVTransformHydrostatic`
    % + `WVTransformBoussinesq`
    % + `WVTransformStratifiedQG`
    % + `WVTransformBarotropicQG`
    %
    % Wave-bearing transforms store positive- and negative-frequency wave and
    % inertial coefficients in `Ap` and `Am`, and zero-frequency geostrophic
    % and mean-density-anomaly coefficients in `A0`. Quasigeostrophic
    % transforms use `A0` only. The `Apt`, `Amt`, and `A0t` variables are the
    % corresponding coefficients evaluated at the current transform time.
    %
    % - Topic: Create and restore a transform
    % - Topic: Inspect the domain
    % - Topic: Inspect the domain — Physical environment
    % - Topic: Inspect the domain — Physical environment — Planetary rotation
    % - Topic: Inspect the domain — Physical environment — Stratification and reference density
    % - Topic: Inspect the domain — Physical environment — Gravity
    % - Topic: Inspect the domain — Spatial grid
    % - Topic: Inspect the domain — Spatial grid — Coordinate axes
    % - Topic: Inspect the domain — Spatial grid — Coordinate arrays
    % - Topic: Inspect the domain — Spatial grid — Domain dimensions
    % - Topic: Inspect the domain — Spatial grid — Resolution and shape
    % - Topic: Inspect the domain — Spatial grid — Quadrature and integration
    % - Topic: Inspect the domain — Spectral grid
    % - Topic: Inspect the domain — Spectral grid — Compact grid vectors
    % - Topic: Inspect the domain — Spectral grid — Compact grid arrays
    % - Topic: Inspect the domain — Spectral grid — Wavenumber spacing
    % - Topic: Inspect the domain — Spectral grid — Horizontal wavenumber geometry
    % - Topic: Inspect the domain — Spectral grid — Resolution and shape
    % - Topic: Inspect the domain — Spectral grid — Vertical modes and scaling
    % - Topic: Inspect the domain — Spectral grid — Vertical-mode transformation matrices
    % - Topic: Inspect the domain — Transform configuration
    % - Topic: Initialize the flow
    % - Topic: Initialize the flow — General initialization
    % - Topic: Initialize the flow — Waves
    % - Topic: Initialize the flow — Waves — Individual modes
    % - Topic: Initialize the flow — Waves — Wave spectra
    % - Topic: Initialize the flow — Inertial oscillations
    % - Topic: Initialize the flow — Geostrophic motions
    % - Topic: Initialize the flow — Mean density anomalies
    % - Topic: Evaluate physical fields
    % - Topic: Evaluate physical fields — Registered variables
    % - Topic: Evaluate physical fields — On the model grid
    % - Topic: Evaluate physical fields — On the model grid — Velocity
    % - Topic: Evaluate physical fields — On the model grid — Density and displacement
    % - Topic: Evaluate physical fields — On the model grid — Pressure and surface fields
    % - Topic: Evaluate physical fields — On the model grid — Vorticity and geostrophic fields
    % - Topic: Evaluate physical fields — At arbitrary positions
    % - Topic: Evaluate physical fields — Isopycnal utilities
    % - Topic: Manage forcing and closures
    % - Topic: Manage forcing and closures — Configure forcing
    % - Topic: Manage forcing and closures — Inspect forcing and closures
    % - Topic: Manage forcing and closures — Summarize forcing
    % - Topic: Analyze the flow
    % - Topic: Analyze the flow — Flow diagnostics
    % - Topic: Analyze the flow — Density validity
    % - Topic: Analyze the flow — Potential vorticity and enstrophy
    % - Topic: Analyze the flow — Spectra
    % - Topic: Analyze the flow — Spectra — Spectral fields
    % - Topic: Analyze the flow — Spectra — Radial wavenumber
    % - Topic: Analyze the flow — Spectra — Frequency
    % - Topic: Analyze energy
    % - Topic: Analyze energy — Component energy
    % - Topic: Analyze energy — Total energy
    % - Topic: Analyze energy — Energy summaries
    % - Topic: Save transform state
    % - Topic: Convert representations
    % - Topic: Convert representations — Physical fields and coefficients
    % - Topic: Differentiate and integrate fields
    % - Topic: Inspect flow components
    % - Topic: Inspect flow components — Primary flow components
    % - Topic: Inspect flow components — Registered and combined components
    % - Topic: Inspect flow components — Summarize flow components
    % - Topic: Inspect wave-vortex coefficients
    % - Topic: Inspect wave-vortex coefficients — Stored coefficients
    % - Topic: Inspect wave-vortex coefficients — Coefficients at the current time
    % - Topic: Inspect wave-vortex coefficients — Coefficient evolution
    % - Topic: Create a related transform
    % - Topic: Extend a transform
    % - Topic: Extend a transform — Flow components
    % - Topic: Extend a transform — Operations and variables
    % - Topic: Get package information
    % - Topic: Projection and reconstruction coefficients
    % - Topic: Geometry and mode indexing
    % - Topic: Geometry and mode indexing — Mode numbers and validity
    % - Topic: Geometry and mode indexing — Linear-index conversion
    % - Topic: Geometry and mode indexing — DFT and WV layout metadata
    % - Topic: Geometry and mode indexing — Layout conversion
    % - Topic: Geometry and mode indexing — Masks and Hermitian bookkeeping
    % - Topic: Geometry and mode indexing — Additional geometry utilities
    % - Topic: Spectral transforms and operators
    % - Topic: Nonlinear flux and forcing internals
    % - Topic: Persistence internals
    % - Topic: Caches and registries
    % - Topic: Class internals
    % - Topic: Construction internals
    %
    % - Declaration: classdef WVTransform < matlab.mixin.indexing.RedefinesDot & CAAnnotatedClass
    
    % Public read and write properties
    properties (GetAccess=public, SetAccess=public)
        % Current transform time in seconds.
        %
        % - Topic: Domain attributes
        t = 0

        % Reference time for the stored wave phases, in seconds.
        %
        % - Topic: Domain attributes
        t0 = 0

        % Positive-frequency wave and inertial coefficients at reference time `t0`.
        %
        % `Ap` is a complex array with the transform's spectral layout. Only
        % locations selected by the primary wave and inertial component masks
        % are active. Together with `Am`, it reconstructs a real physical
        % state and has units of velocity.
        %
        % - Topic: Wave-vortex coefficients
        Ap = 0

        % Negative-frequency wave and inertial coefficients at reference time `t0`.
        %
        % `Am` is a complex array with the transform's spectral layout. Only
        % locations selected by the primary wave and inertial component masks
        % are active. Its conjugacy relations with `Ap` enforce a real
        % physical state, including `Am = conj(Ap)` on inertial modes.
        %
        % - Topic: Wave-vortex coefficients
        Am = 0

        % Zero-frequency geostrophic and mean-density-anomaly coefficients.
        %
        % `A0` is a complex spectral array with units of streamfunction,
        % $$\mathrm{m^{2}\,s^{-1}}$$. Geostrophic modes occupy nonzero
        % horizontal wavenumbers; mean-density-anomaly modes occupy the
        % horizontally averaged internal-mode locations. Quasigeostrophic
        % transforms use `A0` without `Ap` or `Am` wave content.
        %
        % - Topic: Wave-vortex coefficients
        A0 = 0
    end

    properties (GetAccess=public, SetAccess=public)
        Apm_TE_factor, A0_TE_factor, A0_TZ_factor, A0_QGPV_factor, A0_Psi_factor, A0_KE_factor, A0_PE_factor
    end
    % Public read-only properties
    properties (GetAccess=public, SetAccess=protected)
        forcingType
    end

    properties (GetAccess=public, SetAccess=private, Dependent)
        % Installed WaveVortexModel version.
        %
        % - Topic: Get package information
        version
    end

    properties (Abstract)
        totalEnergySpatiallyIntegrated
        totalEnergy
        isHydrostatic
    end

    properties (Dependent, SetAccess=private)
        hasClosure
        primaryFlowComponents
        nFluxedComponents
        forcing
    end

    properties %(Access=private)
        operationNameMap
        operationVariableNameMap
        timeDependentVariablesNameMap
        wvCoefficientDependentVariablesNameMap
        variableCache

        primaryFlowComponentNameMap
        flowComponentNameMap
        totalFlowComponent

        hasPVComponent logical
        hasWaveComponent logical

        forcingNameMap
        spatialFluxForcing WVForcing = WVForcing.empty(1,0)
        spectralFluxForcing WVForcing = WVForcing.empty(1,0)
        spectralAmplitudeForcing WVForcing = WVForcing.empty(1,0)
    end

    events
        forcingDidChange
    end
    
    methods (Abstract)
        % Construct the same transform family at a requested resolution.
        %
        % - Topic: Initialization
        wvtX2 = waveVortexTransformWithResolution(self,m)
        
        % Required for transformUVEtaToWaveVortex 
        % u_bar = transformFromSpatialDomainWithFio(self,u)
        u_bar = transformFromSpatialDomainWithFg(self, u)
        w_bar = transformFromSpatialDomainWithGg(self, w)
        % w_bar = transformWithG_wg(self, w_bar )

        % Required for transformWaveVortexToUVEta
        u = transformToSpatialDomainWithF(self, options)
        w = transformToSpatialDomainWithG(self, options )

        [Fp,Fm,F0] = nonlinearFlux(self)
    end
    
    methods (Static,Abstract)
        names = spatialDimensionNames()
        names = spectralDimensionNames()
    end

    methods (Access=protected)
        function varargout = dotReference(self,indexOp)
            % Typically the request will be directly for a WVOperation,
            % but sometimes it will be for a variable that can only be
            % produced as a bi-product of some operation.
            if isKey(self.operationVariableNameMap,indexOp(1).Name)
                varargout{1} = self.variableWithName(indexOp(1).Name);
                if length(indexOp) > 1
                    varargout{1} = varargout{1}.(indexOp(2:end));
                end
            elseif isKey(self.operationNameMap,indexOp(1).Name)
                op = self.operationNameMap{indexOp(1).Name};
                varargout = cell(1,op.nVarOut);
                [varargout{:}] = self.performOperation(op);
            else
                error("WVTransform:UnknownVariable","The variable %s does not exist",indexOp(1).Name);
            end
            
        end

        function self = dotAssign(self,indexOp,varargin)
            error("The WVVariableAnnotation %s is read-only.",indexOp(1).Name)
        end

        function n = dotListLength(self,indexOp,indexContext)
            if isKey(self.operationNameMap,indexOp(1).Name)
                modelOp = self.operationNameMap(indexOp(1).Name);
                n = modelOp{1}.nVarOut;
            else
                n=1;
            end
        end
    end

    methods
        function self = WVTransform(forcingType)
            % Initialize the internal WVTransform state for a concrete subclass.
            %
            % Concrete transform constructors call this method; users create
            % one of the supported subclasses instead.
            %
            % - Topic: Developer
            % - Developer: true
            arguments
                forcingType WVForcingType {mustBeNonempty}
            end
            
            self.variableCache = configureDictionary("string","cell");
            self.operationVariableNameMap = configureDictionary("string","WVVariableAnnotation"); %containers.Map(); % contains names of variables with associated operations
            self.operationNameMap = configureDictionary("string","cell");
            self.timeDependentVariablesNameMap = configureDictionary("string","cell");
            self.wvCoefficientDependentVariablesNameMap = configureDictionary("string","cell");
            
            self.updateDependentVariablesNameMap([],[]);
            addlistener(self,'propertyAnnotationsDidChange',@self.updateDependentVariablesNameMap);

            self.primaryFlowComponentNameMap = configureDictionary("string","cell");
            self.flowComponentNameMap = configureDictionary("string","cell");

            self.forcingNameMap = configureDictionary("string","cell");

            if length(intersect(WVForcing.spatialFluxTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spatial flux forcing type.")
            end
            if length(intersect(WVForcing.spectralFluxTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spectral flux forcing type.")
            end
            if length(intersect(WVForcing.spectralAmplitudeTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spectral amplitude forcing type.")
            end
            self.forcingType = forcingType;
        end

        function updateDependentVariablesNameMap(self,~,~)
            self.timeDependentVariablesNameMap = configureDictionary("string","cell");
            self.wvCoefficientDependentVariablesNameMap = configureDictionary("string","cell");
            annotations = self.propertyAnnotations;
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    if annotations(i).isVariableWithLinearTimeStep && ~isKey(self.timeDependentVariablesNameMap,annotations(i).name)
                        self.timeDependentVariablesNameMap{annotations(i).name} = annotations(i);
                    end
                    if annotations(i).isDependentOnApAmA0 && ~isKey(self.wvCoefficientDependentVariablesNameMap,annotations(i).name)
                        self.wvCoefficientDependentVariablesNameMap{annotations(i).name} = annotations(i);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wvtX2 = waveVortexTransformWithDoubleResolution(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Flow components
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Register a primary flow-component extension.
        %
        % - Topic: Flow components
        addPrimaryFlowComponent(self,primaryFlowComponent)
        % Return the registered primary flow-component names.
        %
        % - Topic: Flow components
        names = primaryFlowComponentNames(self)
        % Return a primary flow component by its short name.
        %
        % - Topic: Flow components
        val = primaryFlowComponentWithName(self,name)
        % Registered primary flow components.
        %
        % - Topic: Flow components
        function components = get.primaryFlowComponents(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                components WVPrimaryFlowComponent
            end
            components = [self.primaryFlowComponentNameMap{self.primaryFlowComponentNameMap.keys}];
        end

        % Register a diagnostic flow-component extension.
        %
        % - Topic: Flow components
        addFlowComponent(self,flowComponent)
        % Return all registered flow-component names.
        %
        % - Topic: Flow components
        names = flowComponentNames(self)
        % Return a registered flow component by short name.
        %
        % - Topic: Flow components
        val = flowComponentWithName(self,name)

        function n = get.nFluxedComponents(self)
            n = 2*self.hasWaveComponent + self.hasPVComponent;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Operations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Register an operation and its annotated output variables.
        %
        % - Topic: Operations
        addOperation(self,operation,options)
        % Remove a registered operation.
        %
        % - Topic: Operations
        removeOperation(self,transformOperation)
        % Return a registered operation by name.
        %
        % - Topic: Operations
        val = operationWithName(self,name)

        varargout = performOperation(self,modelOp)
        varargout = performOperationWithName(self,opName)

        % Evaluate one or more registered state variables.
        %
        % - Topic: State variables
        [varargout] = variableWithName(self, variableNames);
        
        % Access dynamical variables at arbitrary positions.
        %
        % The interpolation method may be `linear` or `spline`. Horizontal
        % coordinates are wrapped periodically before interpolation.
        % - Topic: State variables
        [varargout] = variableAtPositionWithName(self,x,y,z,variableNames,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variable cache
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addToVariableCache(self,name,var)
        clearVariableCacheOfApAmA0DependentVariables(self)
        clearVariableCacheOfTimeDependentVariables(self)
        removeFromVariableCache(self,name)
        varargout = fetchFromVariableCache(self,varargin)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Forcing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bool = get.hasClosure(self)
            bool = false;
            for iForce=1:length(self.forcing)
                bool = bool | self.forcing(iForce).isClosure;
            end
        end

        function version = get.version(self)
            arguments
                self WVTransform
            end
            version = WVTransform.cachedPackageVersion();
        end

        function removeAllForcing(self)
            % Remove every forcing and closure from this transform.
            %
            % - Topic: Forcing
            arguments
                self WVTransform {mustBeNonempty}
            end
            self.removeForcing(self.forcing);
            % self.forcingNameMap = configureDictionary("string","cell");
            % self.spatialFluxForcing = WVForcing.empty(1,0);
            % self.spectralFluxForcing = WVForcing.empty(1,0);
            % self.spectralAmplitudeForcing = WVForcing.empty(1,0);
            % notify(self,'forcingDidChange');
        end

        function setForcing(self,force)
            % Replace the complete forcing registry.
            %
            % - Topic: Forcing
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end
            [spatialForcing,spectralForcing,amplitudeForcing,nameMap] = self.stagedForcingRegistry(force);
            self.commitForcingRegistry(spatialForcing,spectralForcing,amplitudeForcing,nameMap);
        end

        function addForcing(self,force)
            % Add forcing or closure objects to this transform.
            %
            % - Topic: Forcing
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end

            requestedForcing = [self.forcing reshape(force,1,[])];
            [spatialForcing,spectralForcing,amplitudeForcing,nameMap] = self.stagedForcingRegistry(requestedForcing);
            self.commitForcingRegistry(spatialForcing,spectralForcing,amplitudeForcing,nameMap);
        end

        function removeForcing(self,force)
            % Remove the exact registered forcing objects.
            %
            % - Topic: Forcing
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end
            for iForce = 1:length(force)
                aForce = force(iForce);
                if ~isKey(self.forcingNameMap,aForce.name) || self.forcingNameMap{aForce.name} ~= aForce
                    error("The requested forcing is not registered with this transform.")
                end
            end

            retainedForcing = self.forcing;
            for iForce = 1:length(force)
                retainedForcing(retainedForcing == force(iForce)) = [];
            end
            [spatialForcing,spectralForcing,amplitudeForcing,nameMap] = self.stagedForcingRegistry(retainedForcing);
            self.commitForcingRegistry(spatialForcing,spectralForcing,amplitudeForcing,nameMap);
        end

        function forcing = get.forcing(self)
            % Registered forcing and closure objects in execution order.
            %
            % - Topic: Forcing
            forcing = cat(2,self.spatialFluxForcing,self.spectralFluxForcing,self.spectralAmplitudeForcing);
        end

        function bool = hasForcingWithName(self,name)
            % Test whether forcing objects are registered by name.
            %
            % - Topic: Forcing
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                bool logical
            end
            bool = isKey(self.forcingNameMap,string(name));
        end

        function names = forcingNames(self)
            % Return forcing and closure names in application order.
            %
            % Names follow the physical-space flux, spectral-flux, and
            % spectral-amplitude stages used by the nonlinear tendency.
            %
            % - Topic: Forcing
            arguments (Input)
                self WVTransform {mustBeNonempty}
            end
            arguments (Output)
                names (:,1) string
            end
            names = [self.spatialFluxForcing.name, self.spectralFluxForcing.name, self.spectralAmplitudeForcing.name];
        end

        function forcing = forcingWithName(self,name)
            % Return registered forcing objects by name.
            %
            % - Topic: Forcing
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                forcing WVForcing
            end
            forcing = WVForcing.empty(1,0);
            for iName = 1:length(name)
                requestedName = name{iName};
                if ~isKey(self.forcingNameMap,requestedName)
                    error("No forcing named '%s' is registered with this transform.",requestedName)
                end
                forcing(end+1) = self.forcingNameMap{requestedName};
            end
        end

        function restoreForcingAmplitudes(self)
            if any(ismember(WVForcingType("PVSpectralAmplitude"),self.forcingType))
                for i=1:length(self.spectralAmplitudeForcing)
                    self.A0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralAmplitude(self,self.A0);
                end
            end

            if any(ismember(WVForcingType("SpectralAmplitude"),self.forcingType))
                for i=1:length(self.spectralAmplitudeForcing)
                    [self.Ap, self.Am, self.A0] = self.spectralAmplitudeForcing(i).setSpectralAmplitude(self,self.Ap, self.Am, self.A0);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Domain attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set.t(self,value)
            self.t = value;
            self.clearVariableCacheOfTimeDependentVariables();
        end

        function set.Ap(self,value)
            self.Ap = value;
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        function set.Am(self,value)
            self.Am = value;
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        function set.A0(self,value)
            self.A0 = value;
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N)
        [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)

        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApUMask,AmUMask,A0UMask,ApUxMask,AmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)

        function [Fp,Fm,F0] = rk4NonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)
            function nlF = fluxAtTimeCellArray(t,y0,wvt)
                n = 0;
                wvt.t = t;
                if wvt.hasWaveComponent == true
                    n=n+1; wvt.Ap(:) = y0{n};
                    n=n+1; wvt.Am(:) = y0{n};
                end
                if wvt.hasPVComponent == true
                    n=n+1; wvt.A0(:) = y0{n};
                end

                nlF = cell(1,3);
                [nlF{:}] = wvt.nonlinearFluxForFlowComponents(uFlowComponent,gradUFlowComponent);
            end
            dt = 2*pi/max(abs(self.Omega(:)))/10;
            previousY = {self.Ap;self.Am;self.A0};
            rk4Integrator = WVArrayIntegrator(@(t,y0) fluxAtTimeCellArray(t,y0,self),[self.t self.t+dt],previousY,dt);
            self.Ap =previousY{1};
            self.Am =previousY{2};
            self.A0 =previousY{3};

            Fp = (rk4Integrator.currentY{1} - rk4Integrator.previousY{1})/dt;
            Fm = (rk4Integrator.currentY{2} - rk4Integrator.previousY{2})/dt;
            F0 = (rk4Integrator.currentY{3} - rk4Integrator.previousY{3})/dt;
        end

        function [Fp,Fm,F0,rk4Integrator] = rk4NonlinearFlux(self)
            function nlF = fluxAtTimeCellArray(t,y0,wvt)
                n = 0;
                wvt.t = t;
                if wvt.hasWaveComponent == true
                    n=n+1; wvt.Ap(:) = y0{n};
                    n=n+1; wvt.Am(:) = y0{n};
                end
                if wvt.hasPVComponent == true
                    n=n+1; wvt.A0(:) = y0{n};
                end

                nlF = cell(1,3);
                [nlF{:}] = wvt.nonlinearFlux;
            end
            dt = 2*pi/max(abs(self.Omega(:)))/10;
            previousY = {self.Ap;self.Am;self.A0};
            previousT = self.t;
            rk4Integrator = WVArrayIntegrator(@(t,y0) fluxAtTimeCellArray(t,y0,self),[self.t self.t+dt],previousY,dt);
            self.t = previousT;
            self.Ap =previousY{1};
            self.Am =previousY{2};
            self.A0 =previousY{3};

            Fp = (rk4Integrator.currentY{1} - rk4Integrator.previousY{1})/dt;
            Fm = (rk4Integrator.currentY{2} - rk4Integrator.previousY{2})/dt;
            F0 = (rk4Integrator.currentY{3} - rk4Integrator.previousY{3})/dt;
        end

        [Ep,Em,E0_A,E0_B] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options);
        Z0 = enstrophyFluxFromNonlinearFlux(self,F0,options);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % function energy = totalEnergySpatiallyIntegrated(self)
        %     if self.isHydrostatic == 1
        %         [u,v,eta] = self.variableWithName('u','v','eta');
        %         energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        %     else
        %         [u,v,w,eta] = self.variableWithName('u','v','w','eta');
        %         energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        %     end
        % end
        % 
        % function energy = totalEnergy(self)
        %     energy = sum( self.Apm_TE_factor(:).*( abs(self.Ap(:)).^2 + abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        % end

        function energy = totalEnergyOfFlowComponent(self,flowComponent)
            % Compute the energy carried by one flow component.
            %
            % The component masks select its active `Ap`, `Am`, and `A0`
            % coefficients before the transform's energy factors are summed.
            %
            % ```matlab
            % waveEnergy = wvt.totalEnergyOfFlowComponent(wvt.waveComponent);
            % ```
            %
            % - Topic: Energetics
            % - Declaration: energy = totalEnergyOfFlowComponent(flowComponent)
            % - Parameter flowComponent: component whose coefficient masks select the energy
            % - Returns energy: horizontally averaged, depth-integrated energy per unit reference density
            arguments (Input)
                self WVTransform
                flowComponent WVFlowComponent
            end
            arguments (Output)
                energy (1,1) double
            end
            energy = 0;
            if flowComponent.hasWaveComponent
                energy = energy + sum(self.Apm_TE_factor(:).*( flowComponent.maskAp(:).*abs(self.Ap(:)).^2 + flowComponent.maskAm(:).*abs(self.Am(:)).^2 ));
            end
            if flowComponent.hasPVComponent
                energy = energy +  sum(self.A0_TE_factor(:).*( flowComponent.maskA0(:).*abs(self.A0(:)).^2));
            end   
        end


        % function variable = dynamicalVariable(self,variableName,options)
        %     arguments(Input)
        %         self WVTransform {mustBeNonempty}
        %     end
        %     arguments (Input,Repeating)
        %         variableName char
        %     end
        %     arguments (Input)
        %         options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
        %     end
        %     arguments (Output)
        %         variable (:,1) cell
        %     end
        %     variable = cell(size(variableName));
        %     if ~isempty(options.flowComponent)
        %         for iVar=1:length(variableName)
        %             variableName{iVar} = append(variableName{iVar},'_',options.flowComponent.abbreviatedName);
        %         end
        %     end
        % end

        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithF(Apm=options.Apm,A0=options.A0);
            ux = self.diffX(u);
            uy = self.diffY(u);
            uz = self.diffZF(u);
        end

        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithG(Apm=options.Apm,A0=options.A0);
            wx = self.diffX(w);
            wy = self.diffY(w);
            wz = self.diffZG(w);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function names = variableNames(self)
            % Return the names of all registered state variables.
            %
            % - Topic: State variables
            annotations = self.propertyAnnotations;
            names = {};
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    names{end+1} = annotations(i).name;
                end
            end
        end

        function bool = hasVariableWithName(self,name)
            % Test whether state variables are registered by name.
            %
            % - Topic: State variables
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                bool logical
            end
            bool = ismember(string(name),self.variableNames);
        end

        summarizeEnergyContent(self)
        summarizeDegreesOfFreedom(self)
        summarizeModeEnergy(self,options)

        function summarizeVariables(self)
            % Print a table of registered state variables and cache status.
            %
            % - Topic: State variables
            annotations = self.propertyAnnotations;
            variableAnnotationNameMap = configureDictionary("string","cell");
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    variableAnnotationNameMap{annotations(i).name} = annotations(i);
                end
            end
            Dimension = cell(variableAnnotationNameMap.numEntries,1);
            Units = cell(variableAnnotationNameMap.numEntries,1);
            Description = cell(variableAnnotationNameMap.numEntries,1);
            Name = keys(variableAnnotationNameMap,'cell');
            Cached = cell(variableAnnotationNameMap.numEntries,1); %variableCache
            for iVar=1:length(Name)
                if isempty(variableAnnotationNameMap{Name{iVar}}.dimensions)
                    Dimension{iVar} = "()";
                else
                    Dimension{iVar} = join(["(",join(string(variableAnnotationNameMap{Name{iVar}}.dimensions),', '),")"]) ;
                end
                Units{iVar} = variableAnnotationNameMap{Name{iVar}}.units;
                Description{iVar} = variableAnnotationNameMap{Name{iVar}}.description;
                Cached{iVar} = isKey(self.variableCache,Name{iVar});
            end
            Name = string(Name);
            Dimension = string(Dimension);
            Units = string(Units);
            Cached = string(Cached);
            Description = string(Description);
            T = table(Name,Dimension,Units,Cached,Description);
            disp(T);
        end

        function summarizeFlowComponents(self)
            % Print a table of registered primary and diagnostic components.
            %
            % - Topic: Flow components
            Name = cell(self.flowComponentNameMap.numEntries,1);
            isPrimary = cell(self.flowComponentNameMap.numEntries,1);
            FullName = cell(self.flowComponentNameMap.numEntries,1);
            AbbreviatedName = cell(self.flowComponentNameMap.numEntries,1);
            flowComponentNames_ = self.flowComponentNames;
            for iVar = 1:length(flowComponentNames_)
                name = flowComponentNames_(iVar);
                Name{iVar} = name{1};
                isPrimary{iVar} = isKey(self.primaryFlowComponentNameMap,name{1});
                FullName{iVar} = self.flowComponentWithName(name{1}).name;
                AbbreviatedName{iVar} = self.flowComponentWithName(name{1}).abbreviatedName;
            end
            Name = string(Name);
            isPrimary = string(isPrimary);
            FullName = string(FullName);
            AbbreviatedName = string(AbbreviatedName);
            T = table(Name,isPrimary,FullName,AbbreviatedName);
            disp(T);
        end

        function summarizeForcing(self)
            % Print a table of registered forcing and closure objects.
            %
            % - Topic: Forcing
            Name = cell(length(self.forcing),1);
            IsClosure = cell(length(self.forcing),1);
            for iForce=1:length(self.forcing)
                Name{iForce} = self.forcing(iForce).name;
                IsClosure{iForce} = self.forcing(iForce).isClosure;
            end
            Name = string(Name);
            IsClosure = string(IsClosure);
            T = table(Name,IsClosure);
            disp(T);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initializing, adding and removing dynamical features
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   init — clears ALL variables Ap,Am,A0, then sets/adds
        %   set  - clears only the component requested, and sets with new value.
        %   add  - adds to existing component
        %   removeAll – remove all features of given type       

        initFromNetCDFFile(self,ncfile,options)

        function initForcingFromNetCDFFile(self,ncfile)
            arguments
                self WVTransform {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
            end
            % forcingGroupName = join( [string(class(self)),"forcing"],"-");
            % group = ncfile.groupWithName(class(self));
            group = ncfile;
            f = @(className,group) feval(strcat(className,'.forcingFromGroup'),group, self);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,{"forcing"},classConstructor=f,shouldIgnoreMissingProperties=true);
            if isfield(vars,"forcing")
                self.setForcing(vars.forcing);
            else
                self.removeAllForcing();
            end
        end

        initWithUVRho(self,u,v,rho,t)
        initWithUVEta(self,U,V,N)
        addUVEta(self,U,V,N)

        initWithRandomFlow(self,flowComponentNames)
        addRandomFlow(self,flowComponentNames)
        
        removeAll(self)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove internal waves from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ncfile = writeToFile(self,path,props,options)

        function flag = hasMeanPressureDifference(self)
            % Diagnose an MDA mean-pressure difference between the boundaries.
            %
            % Only the mean-density-anomaly (MDA) component can contribute
            % to a horizontally averaged pressure difference between the
            % top and bottom boundaries. The diagnostic evaluates
            %
            % $$
            % \Delta \bar p_{\mathrm{mda}} =
            % \left|\langle p_{\mathrm{mda,top}}\rangle_{xy}
            % -\langle p_{\mathrm{mda,bottom}}\rangle_{xy}\right|
            % $$
            %
            % relative to the maximum absolute MDA pressure. It returns
            % `true` when the relative difference is greater than `1e-5`.
            % Transforms without an MDA component return `false`. Other
            % flow components and common pressure-gauge offsets do not
            % enter the calculation.
            %
            % - Topic: Initial Conditions
            % - Declaration: flag = hasMeanPressureDifference()
            % - Returns flag: scalar logical indicating a resolved MDA boundary-pressure difference
            arguments (Input)
                self WVTransform {mustBeNonempty}
            end
            arguments (Output)
                flag (1,1) logical
            end

            if ~isa(self,"WVMeanDensityAnomalyMethods")
                flag = false;
                return
            end

            mdaCoefficients = self.mdaComponent.maskA0.*self.A0;
            if ~any(mdaCoefficients ~= 0,"all")
                flag = false;
                return
            end

            pressureOperation = self.operationForKnownVariable('p',flowComponent=self.mdaComponent);
            pMDA = pressureOperation.compute(self);
            pressureScale = max(abs(pMDA),[],"all");
            if pressureScale == 0
                flag = false;
                return
            end

            meanPressureDifference = abs(mean(pMDA(:,:,end),"all") - mean(pMDA(:,:,1),"all"));
            flag = meanPressureDifference > 1e-5*pressureScale;
        end

        [varargout] = spectralVariableWithResolution(self,wvtX2,varargin)
    end

    methods (Access=private)
        function [spatialForcing,spectralForcing,amplitudeForcing,nameMap] = stagedForcingRegistry(self,forcing)
            spatialForcing = WVForcing.empty(1,0);
            spectralForcing = WVForcing.empty(1,0);
            amplitudeForcing = WVForcing.empty(1,0);
            nameMap = configureDictionary("string","cell");

            for iForce = 1:length(forcing)
                aForce = forcing(iForce);
                if aForce.wvt ~= self
                    error("Every forcing must be initialized with the transform receiving it.")
                end

                compatibleTypes = intersect(aForce.forcingType,self.forcingType);
                isSpatial = any(ismember(compatibleTypes,WVForcing.spatialFluxTypes()));
                isSpectral = any(ismember(compatibleTypes,WVForcing.spectralFluxTypes()));
                isAmplitude = any(ismember(compatibleTypes,WVForcing.spectralAmplitudeTypes()));
                if sum([isSpatial isSpectral isAmplitude]) ~= 1
                    error("The transform does not support exactly one forcing category for '%s'.",aForce.name)
                end

                if isKey(nameMap,aForce.name)
                    existingForce = nameMap{aForce.name};
                    if existingForce == aForce
                        continue
                    end
                    spatialForcing(spatialForcing == existingForce) = [];
                    spectralForcing(spectralForcing == existingForce) = [];
                    amplitudeForcing(amplitudeForcing == existingForce) = [];
                end

                if isSpatial
                    spatialForcing(end+1) = aForce;
                elseif isSpectral
                    spectralForcing(end+1) = aForce;
                else
                    amplitudeForcing(end+1) = aForce;
                end
                nameMap{aForce.name} = aForce;
            end

            spatialForcing = WVTransform.sortForcingByPriority(spatialForcing);
            spectralForcing = WVTransform.sortForcingByPriority(spectralForcing);
            amplitudeForcing = WVTransform.sortForcingByPriority(amplitudeForcing);
        end

        function commitForcingRegistry(self,spatialForcing,spectralForcing,amplitudeForcing,nameMap)
            oldForcing = self.forcing;
            newForcing = [spatialForcing spectralForcing amplitudeForcing];
            if WVTransform.areIdenticalForcingArrays(oldForcing,newForcing)
                return
            end

            self.spatialFluxForcing = spatialForcing;
            self.spectralFluxForcing = spectralForcing;
            self.spectralAmplitudeForcing = amplitudeForcing;
            self.forcingNameMap = nameMap;

            for iForce = 1:length(oldForcing)
                if isempty(newForcing) || ~any(newForcing == oldForcing(iForce))
                    oldForcing(iForce).didGetRemovedFromTransform(self);
                end
            end
            notify(self,'forcingDidChange');
        end
    end

    methods (Access=protected)
        % protected — Access from methods in class or subclasses
        varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin);
    end

    methods (Static)
        % Initialize the a transform from file
        [wvt,ncfile] = waveVortexTransformFromFile(path,options)

        operations = classDefinedOperationForKnownVariable(variableName)
        propertyAnnotations = propertyAnnotationForKnownVariable(variableName,options)
        [transformToSpatialDomainWithF,transformToSpatialDomainWithG,mask,isMasked] = optimizedTransformsForFlowComponent(primaryFlowComponents,flowComponent)

        function [propertyAnnotations] = propertyAnnotationsForTransform(variableName,options)
            % return array of CAPropertyAnnotations for the WVTransform
            %
            % This function returns annotations for all properties defined
            % by the WVTransform. It selectively returns annotations for
            % the wave-vortex coefficients, as not all subclass will handle
            % these coefficients in the same way.
            %
            % - Topic: Developer
            % - Declaration: [propertyAnnotations,A0Prop,ApProp,AmProp] = WVTransform.propertyAnnotationsForTransform()
            % - Returns propertyAnnotations: array of CAPropertyAnnotation instances
            % - Returns A0Prop: CANumericProperty instance for A0
            % - Returns ApProp: CANumericProperty instance for Ap
            % - Returns AmProp: CANumericProperty instance for Am
            arguments (Input,Repeating)
                variableName char
            end
            arguments (Input)
                options.spectralDimensionNames = {'j','kl'}
            end
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CAObjectProperty('forcing','array of WVForcing objects');

            propertyAnnotations(end+1) = CANumericProperty('t',{}, 's', 'time of observations');
            propertyAnnotations(end).attributes('standard_name') = 'time';
            propertyAnnotations(end).attributes('axis') = 'T';

            propertyAnnotations(end+1) = CANumericProperty('t0',{},'s', 'reference time of Ap, Am, A0');

            annotation = WVVariableAnnotation('totalEnergy',{},'m3 s-2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            annotation = WVVariableAnnotation('totalEnergySpatiallyIntegrated',{},'m3 s-2', 'horizontally-averaged depth-integrated energy computed in the spatial domain');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            for iVar = 1:length(variableName)
                name = variableName{iVar};
                switch name
                    case 'A0'
                        prop = WVVariableAnnotation('A0',options.spectralDimensionNames,'m2 s-1', 'geostrophic coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'Ap'
                        prop = WVVariableAnnotation('Ap',options.spectralDimensionNames,'m s-1', 'positive wave coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'Am'
                        prop = WVVariableAnnotation('Am',options.spectralDimensionNames,'m s-1', 'negative wave coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'A0_TE_factor'
                        prop = CANumericProperty('A0_TE_factor',options.spectralDimensionNames,'m-1', 'multiplicative factor that multiplies $$A_0^2$$ to compute total energy.',isComplex=0);
                    case 'A0_KE_factor'
                        prop = CANumericProperty('A0_KE_factor',options.spectralDimensionNames,'m-1', 'multiplicative factor that multiplies $$A_0^2$$ to compute kinetic energy.',isComplex=0);
                    case 'A0_PE_factor'
                        prop = CANumericProperty('A0_PE_factor',options.spectralDimensionNames,'m-1', 'multiplicative factor that multiplies $$A_0^2$$ to compute potential energy.',isComplex=0);
                    case 'A0_QGPV_factor'
                        prop = CANumericProperty('A0_QGPV_factor',options.spectralDimensionNames,'m-2', 'multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).',isComplex=0);
                    case 'A0_Psi_factor'
                        prop = CANumericProperty('A0_Psi_factor',options.spectralDimensionNames,'', 'multiplicative factor that maps $$A_0$$ to geostrophic streamfunction. It has no single unit: geostrophic locations are dimensionless and MDA locations have units of $$\mathrm{s\,m^{-1}}$$.',isComplex=0);
                    case 'A0_TZ_factor'
                        prop = CANumericProperty('A0_TZ_factor',options.spectralDimensionNames,'m-3', 'multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.',isComplex=0);
                    case 'Apm_TE_factor'
                        prop = CANumericProperty('Apm_TE_factor',options.spectralDimensionNames,'m', 'multiplicative factor that multiplies $$A_\pm^2$$ to compute total energy.',isComplex=0);

                    otherwise
                        error('There is no variable named %s.',name)
                end
                propertyAnnotations(end+1) = prop;
            end
        end
    end

    methods (Static, Access=private)
        function forcing = sortForcingByPriority(forcing)
            if length(forcing) > 1
                [~,indices] = sort([forcing.priority],"ascend");
                forcing = forcing(indices);
            end
        end

        function tf = areIdenticalForcingArrays(first,second)
            tf = length(first) == length(second);
            if tf && ~isempty(first)
                tf = all(first == second);
            end
        end

        function version = cachedPackageVersion()
            persistent cachedVersion

            if ~isempty(cachedVersion)
                version = cachedVersion;
                return
            end

            classFilePath = which('WVTransform');
            packageRoot = fileparts(fileparts(classFilePath));
            manifestPath = fullfile(packageRoot,'resources','mpackage.json');

            if ~isfile(manifestPath)
                error("WVTransform:MissingPackageManifest", ...
                    "Could not find the WaveVortexModel package manifest at %s.", manifestPath);
            end

            try
                manifest = jsondecode(fileread(manifestPath));
            catch ME
                error("WVTransform:InvalidPackageManifest", ...
                    "Could not parse the WaveVortexModel package manifest at %s: %s", ...
                    manifestPath, ME.message);
            end

            if ~isstruct(manifest) || ~isfield(manifest,'version')
                error("WVTransform:InvalidPackageManifest", ...
                    "The WaveVortexModel package manifest at %s must define a version field.", manifestPath);
            end

            manifestVersion = manifest.version;
            if ~(ischar(manifestVersion) || (isstring(manifestVersion) && isscalar(manifestVersion)))
                error("WVTransform:InvalidPackageManifest", ...
                    "The WaveVortexModel package manifest at %s must define version as a text scalar.", manifestPath);
            end

            version = string(manifestVersion);
            if strlength(version) == 0
                error("WVTransform:InvalidPackageManifest", ...
                    "The WaveVortexModel package manifest at %s must define a nonempty version.", manifestPath);
            end

            cachedVersion = version;
        end
    end
end
