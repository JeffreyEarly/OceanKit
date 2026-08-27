classdef WVTransformStratifiedQG < WVGeometryDoublyPeriodicStratified & WVTransform & WVGeostrophicMethods
    % Represent stratified quasigeostrophic flow with variable stratification.
    %
    % To initialize an instance of the WVTransformStratifiedQG class you
    % must specify the domain size, the number of grid points, and either
    % the density profile or the stratification profile.
    %
    % ```matlab
    % N0 = 3*2*pi/3600;
    % L_gm = 1300;
    % N2 = @(z) N0*N0*exp(2*z/L_gm);
    % wvt = WVTransformStratifiedQG([100e3,100e3,4000],[64,64,65],N2Function=N2,latitude=30);
    % ```
    %
    % The quasigeostrophic state is stored in
    % [`A0`](/classes/transforms/wvtransform/a0.html), with current-time view
    % `A0t`. This transform has no active `Ap`, `Am`, `Apt`, or `Amt` content.
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
    % - Topic: Initialize the flow — Geostrophic motions
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
    %
    % - Declaration: classdef WVTransformStratifiedQG < [WVTransform](/classes/transforms/wvtransform/)
    properties (Dependent)
        totalEnergySpatiallyIntegrated
        totalEnergy
        isHydrostatic
    end

    properties (GetAccess=private,SetAccess=private)
        Fpv, F0
    end

    methods
        function self = WVTransformStratifiedQG(Lxyz, Nxyz, options)
            % Create a stratified quasigeostrophic transform.
            %
            % Creates a new instance of the WVTransformStratifiedQG class
            % for quasigeostrophic flow in variable stratification.
            %
            % Supply either `N2Function` or `rhoFunction`. This transform
            % represents the geostrophic `A0` coefficients and has no wave
            % `Ap` or `Am` content. Additional modal arrays accepted by the
            % constructor are reconstruction state used by persistence.
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformStratifiedQG(Lxyz,Nxyz,options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter options.N2Function: function returning squared buoyancy frequency on `[-Lz,0]`
            % - Parameter options.rhoFunction: function returning density on `[-Lz,0]`
            % - Parameter options.latitude: latitude in the supported domain; default `33`
            % - Parameter options.shouldAntialias: exclude quadratically aliased modes; default `true`
            % - Parameter options.rho0: reference density in kilograms per cubic meter; default `1025`
            % - Returns wvt: new `WVTransformStratifiedQG` instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.shouldAntialias (1,1) logical = true
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.j (:,1) double {mustBeNonempty}
                options.Nj (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double {mustBeSupportedLatitude} = 33
                options.g (1,1) double = 9.81
                
                options.dLnN2 (:,1) double
                options.PF0inv
                options.QG0inv
                options.PF0
                options.QG0
                options.h_0 (:,1) double
                options.P0 (:,1) double
                options.Q0 (:,1) double
                options.z_int (:,1) double
            end

            optionArgs = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, optionArgs{:})
            self@WVTransform(WVForcingType(["PVSpectral","PVSpatial","PVSpectralAmplitude"]));
            self@WVGeostrophicMethods();

            self.initializeGeostrophicComponent();

            % This is not good, I think this should go in the constructor.
            self.addForcing(WVNonlinearAdvection(self));

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = self.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            self.A0 = zeros(self.spectralMatrixSize);
            self.F0 = zeros(self.spectralMatrixSize);
            self.Fpv = zeros(self.spatialMatrixSize);
            if self.geostrophicComponent.normalization ~= "qgpv"
                error("This transform requires the geostrophic component to be normalized the the qgpv norm.");
                % self.A0PV = self.geostrophicComponent.multiplierForVariable(WVCoefficientMatrix.A0,"qgpv-inv");
            end
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformStratifiedQG([self.Lx self.Ly self.Lz],m,optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.A0);
        end

        function wvt = hydrostaticTransform(self)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvt = WVTransformHydrostatic([self.Lx self.Ly self.Lz],[self.Nx self.Ny self.Nz],optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvt);
            end
            wvt.setForcing(forcing);

            wvt.t0 = self.t0;
            wvt.t = self.t;
            wvt.A0 = self.spectralVariableWithResolution(wvt,self.A0);
        end

        function energy = get.totalEnergySpatiallyIntegrated(self)
            [u,v,eta] = self.variableWithName('u','v','eta');
            energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end
        
        function flag = get.isHydrostatic(self)
            flag = true;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear flux computation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function A0 = transformQGPVToWaveVortex(self,qgpv)
            A0 = self.transformFromSpatialDomainWithFg(self.transformFromSpatialDomainWithFourier(qgpv));
        end

        function F0 = nonlinearFlux(self)
            self.Fpv = 0*self.Fpv;
            for i=1:length(self.spatialFluxForcing)
                self.Fpv = self.spatialFluxForcing(i).addPotentialVorticitySpatialForcing(self,self.Fpv);
            end
            self.F0 = self.transformQGPVToWaveVortex(self.Fpv);
            for i=1:length(self.spectralFluxForcing)
                self.F0 = self.spectralFluxForcing(i).addPotentialVorticitySpectralForcing(self,self.F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                self.F0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralForcing(self,self.F0);
            end
            F0 = self.F0;
        end

        function F0 = fluxForForcing(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                F0 dictionary
            end
            F0 = configureDictionary("string","cell");
            self.Fpv = 0*self.Fpv;
            for i=1:length(self.spatialFluxForcing)
               Fpv0 = self.Fpv;
               self.Fpv = self.spatialFluxForcing(i).addPotentialVorticitySpatialForcing(self,self.Fpv);
               F0{self.spatialFluxForcing(i).name} = self.transformQGPVToWaveVortex(self.Fpv-Fpv0);
            end
            self.F0 = self.transformQGPVToWaveVortex(self.Fpv);
            for i=1:length(self.spectralFluxForcing)
                F0_i = self.F0;
                self.F0 = self.spectralFluxForcing(i).addPotentialVorticitySpectralForcing(self,self.F0);
                F0{self.spectralFluxForcing(i).name} = self.F0 - F0_i;
            end
            for i=1:length(self.spectralAmplitudeForcing)
                F0_i = self.F0;
                self.F0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralForcing(self,self.F0);
                F0{self.spectralAmplitudeForcing(i).name} = self.F0 - F0_i;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithFourier(self.PF0inv*(self.P0 .* options.A0));
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(self.QG0inv*(self.Q0 .* options.A0));
        end

    end

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformStratifiedQG.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformStratifiedQG.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformStratifiedQG.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'A0t','uvMax','zeta_z','ssh','ssu','ssv','u','v','eta','pi','p','psi','qgpv','rho_e','rho_total'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformStratifiedQG.spectralDimensionNames();
            spatialDimensionNames = WVTransformStratifiedQG.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformStratifiedQG.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

        function [Lxyz,Nxyz,options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratified.requiredPropertiesForGeometryFromGroup(group);
            % CAAnnotatedClass.throwErrorIfMissingProperties(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % newOptions = namedargs2cell(vars);
            % options = cat(2,geomOptions,newOptions);
            options = geomOptions;
        end

        function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
            % Restore a WVTransformStratifiedQG instance from an existing file
            %
            % This static method is called by WVTransform.waveVortexTransformFromFile
            % and should not need to be called directly.
            %
            % With one output, the temporary NetCDF file is closed before
            % returning. With two outputs, the caller owns the returned
            % NetCDFFile and must close it.
            %
            % - Topic: Initialization (Static)
            % - Declaration: wvt = waveVortexTransformFromFile(path,options)
            % - Parameter path: path to a NetCDF file
            % - Parameter iTime: (optional) time index to initialize from (default 1)
            % - Parameter shouldReadOnly: (optional) open the returned NetCDFFile read-only (default true)
            arguments (Input)
                path char {mustBeFile}
                options.iTime (1,1) double {mustBePositive} = 1
                options.shouldReadOnly logical = true
            end
            arguments (Output)
                wvt WVTransform
                ncfile NetCDFFile
            end
            ncfile = NetCDFFile(path,shouldReadOnly=options.shouldReadOnly);
            try
                wvt = WVTransformStratifiedQG.transformFromGroup(ncfile);
                wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);
                wvt.initForcingFromNetCDFFile(ncfile);
            catch exception
                if ~isempty(ncfile.id)
                    ncfile.close();
                end
                rethrow(exception)
            end
            if nargout < 2
                ncfile.close();
            end
        end


        function wvt = transformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                wvt WVTransform {mustBeNonempty}
            end  
            [Lxy, Nxy, options] = WVTransformStratifiedQG.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformStratifiedQG(Lxy,Nxy,options{:});
        end

    end

end
