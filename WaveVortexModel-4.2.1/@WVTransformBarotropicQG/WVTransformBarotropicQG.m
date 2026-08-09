classdef WVTransformBarotropicQG < WVGeometryDoublyPeriodicBarotropic & WVTransform & WVGeostrophicMethods
    % Represent two-dimensional equivalent-barotropic quasigeostrophic flow.
    %
    % This is a two-dimensional, single-layer transform. The `h` parameter
    % is the equivalent depth; `0.80` m is a representative first-baroclinic
    % value. The transform stores its state in `A0` and has no wave `Ap` or
    % `Am` content.
    %
    % ```matlab
    % Lxy = 50e3;
    % Nxy = 256;
    % latitude = 25;
    % wvt = WVTransformBarotropicQG([Lxy,Lxy],[Nxy,Nxy],h=0.8,latitude=latitude);
    % ```
    %
    % The quasigeostrophic state is stored in
    % [`A0`](/classes/transforms/wvtransform/a0.html), with current-time view
    % `A0t`. This transform has no active `Ap`, `Am`, `Apt`, or `Amt` content.
    %
    % - Topic: Create and restore a transform
    % - Topic: Inspect wave-vortex coefficients
    % - Topic: Inspect wave-vortex coefficients — Coefficients at the current time
    % - Topic: Inspect the domain
    % - Topic: Inspect the domain — Spatial grid
    % - Topic: Inspect the domain — Spectral grid
    % - Topic: Inspect the domain — Rotation and stratification
    % - Topic: Initialize the flow
    % - Topic: Initialize the flow — Geostrophic motions
    % - Topic: Evaluate physical fields
    % - Topic: Evaluate physical fields — On the model grid
    % - Topic: Convert representations
    % - Topic: Convert representations — Physical fields and coefficients
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
    % - Declaration: classdef WVTransformBarotropicQG < [WVTransform](/classes/transforms/wvtransform/)
    properties (Dependent)
        h_0
        totalEnergySpatiallyIntegrated
        totalEnergy
        isHydrostatic
    end

    properties %(GetAccess=private,SetAccess=private)
        Fpv, F0
    end

    methods
        function self = WVTransformBarotropicQG(Lxy, Nxy, options)
            % Create an equivalent-barotropic quasigeostrophic transform.
            %
            % ```matlab
            % Lxy = 50e3;
            % Nxy = 256;
            % wvt = WVTransformBarotropicQG([Lxy,Lxy],[Nxy,Nxy],h=0.8,latitude=30);
            % ```
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformBarotropicQG(Lxy,Nxy,options)
            % - Parameter Lxy: length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter shouldAntialias: (optional) whether or not to de-alias for quadratic multiplications
            % - Parameter options.h: equivalent depth in meters; default `0.8`
            % - Parameter options.latitude: latitude in the supported domain; default `33`
            % - Returns wvt: new `WVTransformBarotropicQG` instance
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.shouldAntialias (1,1) logical = true
                options.rotationRate (1,1) double = 7.2921E-5
                options.planetaryRadius (1,1) double = 6.371e6
                options.latitude (1,1) double {mustBeSupportedLatitude} = 33
                options.g (1,1) double = 9.81
                options.h (1,1) double = 0.8
                options.j (1,1) double {mustBeMember(options.j,[0 1])} = 1
            end
            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,optionCell{:});
            self@WVTransform(WVForcingType(["PVSpectral","PVSpatial","PVSpectralAmplitude"]));
            self@WVGeostrophicMethods();

            % This is not good, I think this should go in the constructor.
            self.addForcing(WVNonlinearAdvection(self));
            
            self.initializeGeostrophicComponent();

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = WVTransformBarotropicQG.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            self.A0 = zeros(self.spectralMatrixSize);
            self.F0 = zeros(self.spectralMatrixSize);
            self.Fpv = zeros(self.spatialMatrixSize);
            if self.geostrophicComponent.normalization ~= "qgpv"
                error("This transform requires the geostrophic component to be normalized the the qgpv norm.");
            end
        end

        function val = get.h_0(self)
            val = self.h;
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
            A0 = self.transformFromSpatialDomainWithFourier(qgpv);
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
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformFromSpatialDomainWithFg(~, u)
        end

        function w = transformFromSpatialDomainWithGg(~, w)
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
            u = self.transformToSpatialDomainWithFourier(options.A0);
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(options.A0);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','h','j','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformBarotropicQG([self.Lx self.Ly],m,optionArgs{:});

            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            wvtX2.A0 = self.spectralVariableWithResolution(wvtX2,self.A0);
        end

        function wvtX2 = waveVortexTransformWithDoubleResolution(self)
            % create a new WVTransform with double resolution
            %
            % - Topic: Initialization
            wvtX2 = self.waveVortexTransformWithResolution(2*[self.Nx self.Ny]);
        end

        function ratio = maxFg(self,k0, l0, j0)
            ratio = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function energy = get.totalEnergySpatiallyIntegrated(self)
            % Return horizontally averaged barotropic energy.
            %
            % The spatial invariant is
            % $$E = \frac{h}{2}\langle u^2+v^2\rangle + \frac{g}{2}\langle\eta^2\rangle.$$
            %
            % - Topic: Energetics
            % - Returns energy: horizontally averaged energy per unit density
            [u,v,eta] = self.variableWithName('u','v','eta');
            energy = self.h*mean(u.^2+v.^2,'all')/2 + self.g*mean(eta.^2,'all')/2;
        end

        function enstrophy = totalEnstrophySpatiallyIntegrated(self)
            % Return horizontally averaged barotropic potential enstrophy.
            %
            % The spatial invariant is
            % $$Z = \frac{h}{2}\langle q_{\mathrm{QG}}^2\rangle.$$
            %
            % - Topic: Energetics
            % - Declaration: enstrophy = totalEnstrophySpatiallyIntegrated()
            % - Returns enstrophy: horizontally averaged potential enstrophy
            enstrophy = self.h*mean(self.qgpv.^2,'all')/2;
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end

        function setSSH(self,ssh,options)
            arguments
                self WVTransformBarotropicQG
                ssh
                options.shouldRemoveMeanPressure double {mustBeMember(options.shouldRemoveMeanPressure,[0 1])} = 0
            end
            if options.shouldRemoveMeanPressure == 1
                sshbar = mean(mean(ssh(self.X,self.Y)));
            else
                sshbar = 0;
            end
            psi = @(X,Y,Z) (self.g/self.f)*(ssh(X,Y)-sshbar);

            self.setGeostrophicStreamfunction(psi);
        end
    end
    
    methods (Access=protected)
        % protected — Access from methods in class or subclasses
        varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin);
    end

    methods (Static)

        function names = spectralDimensionNames()
            % return a cell array of property names required by the class
            %
            % This function returns an array of property names required to be written
            % by the class, in order to restore its state.
            %
            % - Topic: Developer
            % - Declaration:  names = spectralDimensionNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'kl'};
        end

        function names = spatialDimensionNames()
            % return a cell array of the spatial dimension names
            %
            % This function returns an array of dimension names
            %
            % - Topic: Developer
            % - Declaration:  names = spatialDimensionNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'x','y'};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformBarotropicQG.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformBarotropicQG.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicBarotropic.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformBarotropicQG.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'A0t','uvMax','zeta_z','ssh','u','v','eta','pi','psi','qgpv'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformBarotropicQG.spectralDimensionNames();
            spatialDimensionNames = WVTransformBarotropicQG.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicBarotropic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformBarotropicQG.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

        function [Lxy,Nxy,options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options
            end
            [Lxy, Nxy, geomOptions] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryFromGroup(group);
            % CAAnnotatedClass.throwErrorIfMissingProperties(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % newOptions = namedargs2cell(vars);
            % options = cat(2,geomOptions,newOptions);
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
                wvt = WVTransformBarotropicQG.transformFromGroup(ncfile);
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
            [Lxy, Nxy, options] = WVTransformBarotropicQG.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformBarotropicQG(Lxy,Nxy,options{:});
        end

    end
end
