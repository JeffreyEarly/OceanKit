classdef WVTransformConstantStratification < WVGeometryDoublyPeriodicStratifiedConstant & WVTransform & WVGeostrophicMethods & WVInternalGravityWaveMethods & WVInertialOscillationMethods & WVMeanDensityAnomalyMethods
    % Decompose constant-stratification flow into wave and geostrophic components.
    %
    % To initialize an instance of the WVTransformConstantStratification
    % class you must specify the domain size, the number of grid points,
    % and the constant buoyancy frequency.
    %
    % ```matlab
    % N0 = 3*2*pi/3600;
    % wvt = WVTransformConstantStratification([100e3,100e3,4000],[64,64,65],N0=N0,latitude=30);
    % wvtHydrostatic = WVTransformConstantStratification([100e3,100e3,4000],[64,64,65],N0=N0,latitude=30,isHydrostatic=true);
    % ```
    %
    % The transform state is stored in [`Ap`](/classes/transforms/wvtransform/ap.html),
    % [`Am`](/classes/transforms/wvtransform/am.html), and
    % [`A0`](/classes/transforms/wvtransform/a0.html). Their current-time
    % views are `Apt`, `Amt`, and `A0t`.
    %
    % - Topic: Create and restore a transform
    % - Topic: Inspect wave-vortex coefficients
    % - Topic: Inspect wave-vortex coefficients — Coefficients at the current time
    % - Topic: Inspect the domain
    % - Topic: Inspect the domain — Spatial grid
    % - Topic: Inspect the domain — Spectral grid
    % - Topic: Inspect the domain — Rotation and stratification
    % - Topic: Initialize the flow
    % - Topic: Initialize the flow — Waves
    % - Topic: Initialize the flow — Inertial oscillations
    % - Topic: Initialize the flow — Geostrophic motions
    % - Topic: Initialize the flow — Mean density anomalies
    % - Topic: Evaluate physical fields
    % - Topic: Evaluate physical fields — On the model grid
    % - Topic: Convert representations
    % - Topic: Convert representations — Physical fields and coefficients
    % - Topic: Differentiate and integrate fields
    % - Topic: Analyze the flow
    % - Topic: Analyze the flow — Energy and summaries
    % - Topic: Analyze the flow — Potential vorticity and enstrophy
    % - Topic: Analyze the flow — Spectra
    % - Topic: Projection and reconstruction coefficients
    % - Topic: Geometry and mode indexing
    % - Topic: Spectral transforms and operators
    % - Topic: Nonlinear flux and forcing internals
    % - Topic: Persistence internals
    % - Topic: Caches and registries
    % - Topic: Class internals
    %
    % - Declaration: classdef WVTransformConstantStratification < [WVTransform](/classes/transforms/wvtransform/)
    properties (Dependent)
        totalEnergySpatiallyIntegrated
        totalEnergy
    end
    properties (GetAccess=public, SetAccess=public)
        shouldUseTrueNoMotionProfile (1,1) logical = false
    end
    properties
        Fu, Fv, Feta
        nonlinearFluxFunction

        cos_alpha
        sin_alpha
        ApmD_scaled
        ApmW_scaled
    end

    methods
        function self = WVTransformConstantStratification(Lxyz, Nxyz, options)
            % Create a wave-vortex transform for constant stratification.
            %
            % Creates a new instance of the WVTransformConstantStratification
            % class appropriate for disentangling waves and vortices in
            % constant stratification.
            %
            % Set `isHydrostatic=true` for hydrostatic dynamics; the default
            % is the nonhydrostatic transform. `N0` is the constant buoyancy
            % frequency. The remaining geometry options configure the
            % rotating, doubly periodic domain.
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformConstantStratification(Lxyz,Nxyz,options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter options.N0: constant buoyancy frequency in radians per second; default `5.2e-3`
            % - Parameter options.isHydrostatic: use hydrostatic dynamics; default `false`
            % - Parameter options.latitude: latitude in the supported domain; default `33`
            % - Parameter options.shouldAntialias: exclude quadratically aliased modes; default `true`
            % - Parameter options.rho0: reference density in kilograms per cubic meter; default `1025`
            % - Returns wvt: new `WVTransformConstantStratification` instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.shouldAntialias (1,1) logical = true

                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double {mustBeSupportedLatitude} = 33
                options.g (1,1) double = 9.81
                
                options.isHydrostatic logical = false
            end

            optionArgs = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratifiedConstant(Lxyz, Nxyz, optionArgs{:})
            self@WVTransform(WVForcingType(["HydrostaticSpatial","Spectral","SpectralAmplitude"]));
            self@WVGeostrophicMethods();
            self@WVMeanDensityAnomalyMethods();
            self@WVInternalGravityWaveMethods();
            self@WVInertialOscillationMethods();

            self.initializeGeostrophicComponent();
            self.initializeMeanDensityAnomalyComponent();
            self.initializeInternalGravityWaveComponent();
            self.initializeInertialOscillationComponent();

            % This is not good, I think this should go in the constructor.
            self.addForcing(WVNonlinearAdvection(self));

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = self.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            self.addFlowComponent([self.geostrophicComponent,self.waveComponent,self.inertialComponent,self.mdaComponent]);
            self.addOperation(EtaTrueOperation(self));
            self.addOperation(APVOperation());
            self.addOperation(APEOperation(self));
            
            self.A0 = zeros(self.spectralMatrixSize);
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);

            self.Fu=zeros(self.spatialMatrixSize);
            self.Fv=zeros(self.spatialMatrixSize);
            self.Feta=zeros(self.spatialMatrixSize);

            k = shiftdim(self.k,-1);
            l = shiftdim(self.l,-1);
            kappa = sqrt(k.^2 + l.^2);
            self.cos_alpha = k./kappa;
            self.sin_alpha = l./kappa;
            self.cos_alpha(1) = 0;
            self.sin_alpha(1) = 0;

            signNorm = -2*(mod(self.j,2) == 1)+1; % equivalent to (-1)^j
            prefactor = signNorm * sqrt((self.g*self.Lz)/(2*(self.N0*self.N0 - self.f*self.f)));
            mj = (self.j*pi/self.Lz);
            self.ApmD_scaled = (mj/2) .* prefactor;
            self.ApmW_scaled = sqrt(-1) * (kappa/2) .* prefactor;

            if self.isHydrostatic
                self.nonlinearFluxFunction = @() self.nonlinearFluxHydrostatic();
            else
                self.nonlinearFluxFunction = @() self.nonlinearFluxNonhydrostatic();
            end
        end

        function set.shouldUseTrueNoMotionProfile(self,value)
            if self.shouldUseTrueNoMotionProfile == value
                self.shouldUseTrueNoMotionProfile = value;
                return
            end
            self.shouldUseTrueNoMotionProfile = value;
            self.removeFromVariableCache("rho_nm");
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','N0','rho0','planetaryRadius','rotationRate','latitude','g','isHydrostatic'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],m,optionArgs{:});
            wvtX2.shouldUseTrueNoMotionProfile = self.shouldUseTrueNoMotionProfile;
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.A0,wvtX2.Ap,wvtX2.Am] = self.spectralVariableWithResolution(wvtX2,self.A0,self.Ap,self.Am);
        end

        function wvt2 = waveVortexTransformWithExplicitAntialiasing(self)
            if self.shouldAntialias == false
                error("This function only applies to transforms that are dealiasing.")
            end
            names = {'shouldAntialias','N0','rho0','planetaryRadius','rotationRate','latitude','g','isHydrostatic'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
                if names{i} == "shouldAntialias"
                    optionArgs{2*i} = false;
                end
            end
            wvt2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],[self.Nx self.Ny self.Nz],optionArgs{:});
            wvt2.shouldUseTrueNoMotionProfile = self.shouldUseTrueNoMotionProfile;
            wvt2.removeAllForcing();
            wvt2.addForcing(WVAntialiasing(wvt2));

            for iForce=1:length(self.forcing)
                wvt2.addForcing(self.forcing(iForce).forcingWithResolutionOfTransform(wvt2));
            end

            wvt2.t0 = self.t0;
            wvt2.t = self.t;
            [wvt2.A0,wvt2.Ap,wvt2.Am] = self.spectralVariableWithResolution(wvt2,self.A0,self.Ap,self.Am);

        end

        function dx = effectiveHorizontalGridResolution(self)
            %returns the effective grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: flag = effectiveHorizontalGridResolution(other)
            % - Returns effectiveHorizontalGridResolution: double
            arguments
                self WVGeometryDoublyPeriodic
            end
            if self.hasForcingWithName("antialias filter")
                dx = self.forcingWithName("antialias filter").effectiveHorizontalGridResolution;
            else
                dx = effectiveHorizontalGridResolution@WVGeometryDoublyPeriodic(self);
            end
        end

        function j_max = effectiveJMax(self)
            if self.hasForcingWithName("antialias filter")
                j_max = self.forcingWithName("antialias filter").effectiveJMax;
            else
                j_max = effectiveJMax@WVGeometryDoublyPeriodicStratifiedConstant(self);
            end
        end

        function energy = get.totalEnergySpatiallyIntegrated(self)
            if self.isHydrostatic == 1
                [u,v,eta] = self.variableWithName('u','v','eta');
                energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
            else
                [u,v,w,eta] = self.variableWithName('u','v','w','eta');
                energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
            end
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.Apm_TE_factor(:).*( abs(self.Ap(:)).^2 + abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear flux computation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % function [Fp,Fm,F0] = nonlinearFlux(self)
        %     % self.Fu(:)=0;self.Fv(:)=0;self.Feta(:)=0;
        %     % self.Fu=0*self.Fu;self.Fv(:)=0*self.Fv;self.Feta(:)=0*self.Feta;
        %     self.Fu=zeros(self.spatialMatrixSize);self.Fv=zeros(self.spatialMatrixSize);self.Feta=zeros(self.spatialMatrixSize);
        %     for i=1:length(self.spatialFluxForcing)
        %        [self.Fu, self.Fv, self.Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, self.Fu, self.Fv, self.Feta);
        %     end
        %     [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(self.Fu, self.Fv, self.Feta);
        %     for i=1:length(self.spectralFluxForcing)
        %        [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
        %     end
        %     for i=1:length(self.spectralAmplitudeForcing)
        %        [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
        %     end  
        % end
        function [Fp,Fm,F0] = nonlinearFlux(self)
            [Fp,Fm,F0] = self.nonlinearFluxFunction();
        end
        function [Fp,Fm,F0] = nonlinearFluxHydrostatic(self)
            Fu=zeros(self.spatialMatrixSize);Fv=zeros(self.spatialMatrixSize);Feta=zeros(self.spatialMatrixSize); % this isn't good, need to cached
            for i=1:length(self.spatialFluxForcing)
                [Fu, Fv, Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, Fu, Fv, Feta);
            end
            [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu, Fv, Feta);
            for i=1:length(self.spectralFluxForcing)
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
            end
        end

        function [Fp,Fm,F0] = nonlinearFluxNonhydrostatic(self)
            Fu=zeros(self.spatialMatrixSize);Fv=zeros(self.spatialMatrixSize);Fw=zeros(self.spatialMatrixSize);Feta=zeros(self.spatialMatrixSize); % this isn't good, need to cached
            for i=1:length(self.spatialFluxForcing)
                [Fu, Fv, Fw, Feta] = self.spatialFluxForcing(i).addNonhydrostaticSpatialForcing(self, Fu, Fv, Fw, Feta);
            end
            [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu, Fv, Fw, Feta);
            for i=1:length(self.spectralFluxForcing)
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
            end
        end

        function F = fluxForForcing(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                F dictionary
            end
            F = configureDictionary("string","cell");
            if self.isHydrostatic
                Fu=0;Fv=0;Feta=0;
                for i=1:length(self.spatialFluxForcing)
                    Fu0=Fu;Fv0=Fv;Feta0=Feta;
                    [Fu, Fv, Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, Fu, Fv, Feta);
                    [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu-Fu0, Fv-Fv0, Feta-Feta0);
                    F{self.spatialFluxForcing(i).name} = struct("Fp",Fp,"Fm",Fm,"F0",F0);
                end
                [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu, Fv, Feta);
            else
                Fu=0;Fv=0;Fw=0;Feta=0;
                for i=1:length(self.spatialFluxForcing)
                    Fu0=Fu;Fv0=Fv;Fw0=Fw;Feta0=Feta;
                    [Fu, Fv, Fw, Feta] = self.spatialFluxForcing(i).addNonhydrostaticSpatialForcing(self, Fu, Fv, Fw, Feta);
                    [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu-Fu0, Fv-Fv0, Fw-Fw0, Feta-Feta0);
                    F{self.spatialFluxForcing(i).name} = struct("Fp",Fp,"Fm",Fm,"F0",F0);
                end
                [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu, Fv, Fw, Feta);
            end
            for i=1:length(self.spectralFluxForcing)
                Fp_i = Fp; Fm_i = Fm; F0_i = F0;
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
                F{self.spectralFluxForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                Fp_i = Fp; Fm_i = Fm; F0_i = F0;
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
                F{self.spectralAmplitudeForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
            end
        end

        function varargout = spatialFluxForForcingWithName(self,name)
            Fu_name = replace(replace(join( ["Fu_", string(name)],"")," ","_"),"-","_");
            Fv_name = replace(replace(join( ["Fv_", string(name)],"")," ","_"),"-","_");
            Feta_name = replace(replace(join( ["Feta_", string(name)],"")," ","_"),"-","_");
            [Fu_,Fv_,Feta_] = self.variableWithName(Fu_name,Fv_name,Feta_name);
            if self.isHydrostatic && nargout == 3
                varargout = {Fu_,Fv_,Feta_};
            elseif self.isHydrostatic == false && nargout == 4
                Fw_name = replace(replace(join( ["Fw_", string(name)],"")," ","_"),"-","_");
                F_w = self.variableWithName(Fw_name);
                varargout = {Fu_,Fv_,F_w,Feta_};
            end
        end


    end

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformConstantStratification.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformConstantStratification.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratifiedConstant.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformConstantStratification.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','Ap','Am','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'phase','conjPhase','A0t','Apt','Amt','uvMax','wMax','zeta_x','zeta_y','zeta_z','ssh','ssu','ssv','u','v','w','eta','pi','p','psi','qgpv','rho_e','rho_total','rho_nm'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformConstantStratification.spectralDimensionNames();
            spatialDimensionNames = WVTransformConstantStratification.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicStratifiedConstant.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','Ap','Am','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor','A0_Psi_factor','Apm_TE_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformConstantStratification.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

      function [Lxyz, Nxyz, options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratifiedConstant.requiredPropertiesForGeometryFromGroup(group);
            options = geomOptions;
        end

        function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
            % Initialize a WVTransformHydrostatic instance from an existing file
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
                wvt = WVTransformConstantStratification.transformFromGroup(ncfile);
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
            [Lxy, Nxy, options] = WVTransformConstantStratification.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformConstantStratification(Lxy,Nxy,options{:});
        end

    end

end
