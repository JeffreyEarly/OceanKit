classdef WVFixedAmplitudeForcing < WVForcing
    % Hold selected wave-vortex coefficients at prescribed amplitudes.
    %
    % This forcing keeps selected wave-vortex coefficients at prescribed amplitudes while those modes continue to participate in nonlinear interactions.
    %
    % As a simple example, one can set an internal wave mode with amplitude $$1\ \mathrm{cm\,s^{-1}}$$, and that mode will continue to oscillate and maintain its amplitude. The wave will participate in all the nonlinear dynamics, but its amplitude will be maintained/restored at each time step.
    %
    %
    % There are several different ways to write this style of forcing mathematically. The equations of motion, written in the spectral domain, take the following form
    %
    % $$
    % \frac{\partial}{\partial t} A^{klj} = \sum_i F_i^{klj}
    % $$
    %
    % where $$F_i$$ are the contributions from the registered forcing
    % objects. The transform evaluates physical-space forcing, spectral
    % forcing, and then spectral-amplitude forcing. This forcing is evaluated
    % last: it zeros the tendency at selected indices and restores the
    % prescribed coefficient values after the integration step, giving
    % $$\partial_t A^{k\ell j}=0$$ for those modes.
    %
    % In practice, of course, we simply restore the amplitudes to their desired value at the last step, e.g.,
    %
    % ```matlab
    % A0(self.A0_indices) = self.A0bar
    % ```
    %
    % ### Notes
    %
    % - This approach is commonly used in forced-dissipative turbulence to maintain some fixed forcing.
    % - Every fixed mode removes a degree of freedom because it no longer
    % evolves freely. The setter methods therefore ignore coefficients below
    % $$10^{-6}$$ times the largest supplied magnitude unless an explicit
    % mask is provided.
    % - Avoid selecting modes in a closure's damping range. When
    % `WVAdaptiveDamping` is registered, the setter methods automatically
    % remove requested modes with $$K_h>k_\mathrm{damp}$$.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.setGeostrophicModes(kMode=1,lMode=0,j=1,u=0.01);
    % force = WVFixedAmplitudeForcing(wvt,name="geostrophic-mean-flow");
    % force.setGeostrophicForcingCoefficients(wvt.A0);
    % wvt.addForcing(force);
    % ```
    %
    % In practice you can initialize the flow in any way you want with any arbitrary structure, and then pass those coefficients to the forcing. The `WVFixedAmplitudeForcing` looks for coefficients that are small and ignores those.
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Configure the forcing
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    %
    % - Declaration: WVFixedAmplitudeForcing < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Linear indices of the selected `A0` coefficients.
        %
        % This column vector is empty by default and has one entry for each
        % value in `A0bar`.
        %
        % - Topic: Properties
        A0_indices (:,1) uint64 = []

        % Linear indices of the selected `Ap` coefficients.
        %
        % This column vector is empty by default and has one entry for each
        % value in `Apbar`.
        %
        % - Topic: Properties
        Ap_indices (:,1) uint64 = []

        % Linear indices of the selected `Am` coefficients.
        %
        % This column vector is empty by default and has one entry for each
        % value in `Ambar`.
        %
        % - Topic: Properties
        Am_indices (:,1) uint64 = []

        % Prescribed `A0` values in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % Values correspond element-by-element to `A0_indices`.
        %
        % - Topic: Properties
        A0bar (:,1) double = []

        % Prescribed `Ap` values in $$\mathrm{m\,s^{-1}}$$.
        %
        % Values correspond element-by-element to `Ap_indices`.
        %
        % - Topic: Properties
        Apbar (:,1) double = []

        % Prescribed `Am` values in $$\mathrm{m\,s^{-1}}$$.
        %
        % Values correspond element-by-element to `Am_indices`.
        %
        % - Topic: Properties
        Ambar (:,1) double = []
    end

    methods
        function contract = portableImplementationContract(self)
            % Return the paired portable implementation contract.
            %
            % - Topic: Forcing internals
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only forcing contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"ApIndices",self.Ap_indices,"ApValues",self.Apbar,"AmIndices",self.Am_indices,"AmValues",self.Ambar,"A0Indices",self.A0_indices,"A0Values",self.A0bar);
            contract = self.supportedPortableImplementationContract("WVFixedAmplitudeForcing",payload);
        end

        function self = WVFixedAmplitudeForcing(wvt,options)
            % Create fixed-amplitude forcing for selected coefficients.
            %
            % Supply a unique registry name. Coefficients may be selected
            % directly with paired value/index column vectors, or later with
            % `setWaveForcingCoefficients` and
            % `setGeostrophicForcingCoefficients`.
            %
            % See the [WVFixedAmplitudeForcing overview](/classes/forcing/wvfixedamplitudeforcing/)
            % for the tendency and restoration behavior, modeling cautions,
            % and a usage example.
            %
            % - Topic: Initialization
            % - Declaration: self = WVFixedAmplitudeForcing(wvt,options)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Parameter name: required unique forcing-registry name
            % - Parameter Apbar: optional column of prescribed `Ap` values; default empty
            % - Parameter Ambar: optional column of prescribed `Am` values; default empty
            % - Parameter A0bar: optional column of prescribed `A0` values; default empty
            % - Parameter Ap_indices: optional column of corresponding `Ap` linear indices; default empty
            % - Parameter Am_indices: optional column of corresponding `Am` linear indices; default empty
            % - Parameter A0_indices: optional column of corresponding `A0` linear indices; default empty
            % - Returns self: fixed-amplitude forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.name {mustBeText}

                options.Apbar (:,1) double = []
                options.Ambar (:,1) double = []
                options.A0bar (:,1) double = []
                options.A0_indices (:,1) uint64 = []
                options.Ap_indices (:,1) uint64 = []
                options.Am_indices (:,1) uint64 = []
            end
            self@WVForcing(wvt,options.name,WVForcingType(["SpectralAmplitude","PVSpectralAmplitude"]));

            if ~isfield(options,"name")
                error("You must specify a unique name for these spectral masks, e.g., geostrophic mean flow, or M2 tide.")
            end

            canInitializeDirectly = any(options.A0_indices(:)) | any(options.Ap_indices(:)) | any(options.Am_indices(:));
            if canInitializeDirectly == true
                self.Apbar= options.Apbar;
                self.Ambar= options.Ambar;
                self.A0bar= options.A0bar;
                self.A0_indices  = options.A0_indices;
                self.Ap_indices  = options.Ap_indices;
                self.Am_indices  = options.Am_indices;
            end
        end
        function setWaveForcingCoefficients(self,Apbar,Ambar,options)
            % Select positive- and negative-frequency coefficients to fix.
            %
            % Without explicit masks, coefficients whose magnitude is at
            % least $$10^{-6}$$ times the largest supplied magnitude are
            % selected. If adaptive damping is registered, selected modes
            % above its horizontal `k_damp` threshold are removed.
            %
            % - Topic: Setting the forcing
            % - Declaration: setWaveForcingCoefficients(Apbar,Ambar,options)
            % - Parameter Apbar: `Ap` values on the transform spectral grid
            % - Parameter Ambar: `Am` values on the transform spectral grid
            % - Parameter MAp: optional logical `Ap` selection mask; default `abs(Apbar) > 1e-6*max(abs(Apbar(:)))`
            % - Parameter MAm: optional logical `Am` selection mask; default `abs(Ambar) > 1e-6*max(abs(Ambar(:)))`
            arguments
                self WVFixedAmplitudeForcing {mustBeNonempty}
                Apbar (:,:) double {mustBeNonempty}
                Ambar (:,:) double {mustBeNonempty}
                options.MAp (:,:) logical = abs(Apbar) > 1e-6*max(abs(Apbar(:)))
                options.MAm (:,:) logical = abs(Ambar) > 1e-6*max(abs(Ambar(:)))
            end
            if self.wvt.hasForcingWithName("adaptive damping")
                svv = self.wvt.forcingWithName("adaptive damping");
                dampedIndicesAp = options.MAp(self.wvt.Kh > svv.k_damp);
                dampedIndicesAm = options.MAm(self.wvt.Kh > svv.k_damp);
                if any(dampedIndicesAp(:) | dampedIndicesAm(:))
                    warning('You have set %d forcing modes in the damping region. These will be removed.',sum(dampedIndicesAp(:))+sum(dampedIndicesAm(:)));
                    Apbar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MAp(self.wvt.Kh > svv.k_damp) = 0;
                    Ambar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MAm(self.wvt.Kh > svv.k_damp) = 0;
                end
            end

            self.Ap_indices = find(options.MAp);
            self.Apbar = Apbar(self.Ap_indices);

            self.Am_indices = find(options.MAm);
            self.Ambar = Ambar(self.Am_indices);

            fprintf('You are forcing at %d wave modes.\n',length(self.Ap_indices) + length(self.Am_indices));
        end

        function setGeostrophicForcingCoefficients(self,A0bar,options)
            % Select zero-frequency coefficients to fix.
            %
            % Without an explicit mask, coefficients whose magnitude is at
            % least $$10^{-6}$$ times the largest supplied magnitude are
            % selected. If adaptive damping is registered, selected modes
            % above its horizontal `k_damp` threshold are removed.
            %
            % - Topic: Setting the forcing
            % - Declaration: setGeostrophicForcingCoefficients(A0bar,options)
            % - Parameter A0bar: `A0` values on the transform spectral grid
            % - Parameter MA0: optional logical `A0` selection mask; default `abs(A0bar) > 1e-6*max(abs(A0bar(:)))`
            arguments
                self WVFixedAmplitudeForcing {mustBeNonempty}
                A0bar (:,:) double {mustBeNonempty}
                options.MA0 (:,:) logical = abs(A0bar) > 1e-6*max(abs(A0bar(:)))
            end

            if self.wvt.hasForcingWithName("adaptive damping")
                svv = self.wvt.forcingWithName("adaptive damping");
                dampedIndicesA0 = options.MA0(self.wvt.Kh > svv.k_damp);

                if any(dampedIndicesA0(:))
                    warning('You have set %d forcing modes in the damping region. These will be removed.',sum(dampedIndicesA0(:)));
                    A0bar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MA0(self.wvt.Kh > svv.k_damp) = 0;
                end
            end

            self.A0_indices = find(options.MA0);
            self.A0bar = A0bar(self.A0_indices);

            fprintf('You are forcing at %d geostrophic modes.\n',length(self.A0_indices));
        end

        function [model_spectrum, r] = setNarrowBandGeostrophicForcing(self, options)
            % Deprecated 4.x helper for narrow-band geostrophic forcing.
            %
            % New callers should construct `WVNarrowBandGeostrophicForcing`
            % directly. This compatibility entry point remains silent in
            % WaveVortexModel 4.x, initializes the same transform state, and
            % copies the subclass's selected coefficients into `self`.
            %
            % - Topic: Setting the forcing
            % - Declaration: [model_spectrum,r] = setNarrowBandGeostrophicForcing(options)
            % - Parameter r: optional large-scale damping rate in inverse seconds; when supplied, determines `k_r`
            % - Parameter k_r: optional arrest wavenumber in radians per meter; default `2*dk`
            % - Parameter k_f: optional forcing-band center in radians per meter; default `20*dk`
            % - Parameter j_f: optional forced vertical-mode number; default `1`
            % - Parameter u_rms: optional target surface root-mean-square speed in meters per second; default `0.2`
            % - Parameter initialPV: optional initialization choice `"none"`, `"narrow-band"`, or `"full-spectrum"`; default `"narrow-band"`
            % - Returns model_spectrum: configured radial spectrum function
            % - Returns r: effective damping rate in inverse seconds
            arguments
                self WVFixedAmplitudeForcing {mustBeNonempty}
                options.r (1,1) double
                options.k_r (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*2
                options.k_f (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*20
                options.j_f (1,1) double = 1
                options.u_rms (1,1) double = 0.2 % set the *total* energy (not just kinetic) equal to 0.5*u_rms^2
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end

            % TODO(v5): Remove this compatibility helper under Issue #164.
            optionArguments = namedargs2cell(options);
            force = WVNarrowBandGeostrophicForcing(self.wvt,optionArguments{:},name=string(self.name));
            self.A0_indices = force.A0_indices;
            self.A0bar = force.A0bar;
            model_spectrum = force.modelSpectrum;
            r = force.r;
        end
        
        function [Ap, Am, A0] = setSpectralAmplitude(self, wvt, Ap, Am, A0)
            Ap(self.Ap_indices) = self.Apbar;
            Am(self.Am_indices) = self.Ambar;
            A0(self.A0_indices) = self.A0bar;
        end

        function [Fp, Fm, F0] = setSpectralForcing(self, wvt, Fp, Fm, F0)
            Fp(self.Ap_indices) = 0;
            Fm(self.Am_indices) = 0;
            F0(self.A0_indices) = 0;
        end

        function A0 = setPotentialVorticitySpectralAmplitude(self, wvt, A0)
            arguments
                self WVForcing
                wvt WVTransform
                A0 
            end
            A0(self.A0_indices) = self.A0bar;
        end

        function F0 = setPotentialVorticitySpectralForcing(self, wvt, F0)
            arguments
                self WVForcing
                wvt WVTransform
                F0
            end
            F0(self.A0_indices) = 0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            options.name = self.name;
            [options.Ap_indices,options.Apbar] = self.convertSelectedCoefficients(wvtX2,self.Ap_indices,self.Apbar);
            [options.Am_indices,options.Ambar] = self.convertSelectedCoefficients(wvtX2,self.Am_indices,self.Ambar);
            [options.A0_indices,options.A0bar] = self.convertSelectedCoefficients(wvtX2,self.A0_indices,self.A0bar);

            optionArgs = namedargs2cell(options);
            force = WVFixedAmplitudeForcing(wvtX2,optionArgs{:});
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            % Writes this class to a NetCDF group
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: writeToGroup(group,propertyAnnotations,attributes)
            % - Parameter group: NetCDFGroup
            % - Parameter propertyAnnotations: CAPropertyAnnotation
            % - Parameter attributes: configureDictionary("string","string")
            arguments
                self CAAnnotatedClass
                group NetCDFGroup
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end
            % override the logic, and only pass non-zero coefficients.
            properties = {'name'};
            if ~isempty(self.Ap_indices)
                properties = union(properties,{'Ap_indices','Apbar'});
            end
            if ~isempty(self.Am_indices)
                properties = union(properties,{'Am_indices','Ambar'});
            end
            if ~isempty(self.A0_indices)
                properties = union(properties,{'A0_indices','A0bar'});
            end
            propertyAnnotations = self.propertyAnnotationWithName(properties);      
            writeToGroup@CAAnnotatedClass(self,group,propertyAnnotations,attributes);
        end
    end

    methods (Access=protected)
        function [targetIndices,targetValues] = convertSelectedCoefficients(self,wvtX2,sourceIndices,sourceValues)
            selection = zeros(self.wvt.spectralMatrixSize);
            selection(sourceIndices) = 1;
            targetSelection = self.wvt.spectralVariableWithResolution(wvtX2,selection);
            targetIndices = find(targetSelection);

            coefficientValues = complex(zeros(self.wvt.spectralMatrixSize));
            coefficientValues(sourceIndices) = sourceValues;
            targetCoefficientValues = self.wvt.spectralVariableWithResolution(wvtX2,coefficientValues);
            targetValues = targetCoefficientValues(targetIndices);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            vars = {'name','Ap_indices','Apbar','Am_indices','Ambar','A0_indices','A0bar'};
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
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the forcing');
            propertyAnnotations(end+1) = CADimensionProperty('Ap_indices', '1','indices into the Ap matrix');
            propertyAnnotations(end+1) = CANumericProperty('Apbar', {'Ap_indices'}, 'm s-1','prescribed Ap coefficient values',isComplex=true);
            propertyAnnotations(end+1) = CADimensionProperty('Am_indices', '1','indices into the Am matrix');
            propertyAnnotations(end+1) = CANumericProperty('Ambar', {'Am_indices'}, 'm s-1','prescribed Am coefficient values',isComplex=true);
            propertyAnnotations(end+1) = CADimensionProperty('A0_indices', '1','indices into the A0 matrix');
            propertyAnnotations(end+1) = CANumericProperty('A0bar', {'A0_indices'}, 'm2 s-1','prescribed A0 coefficient values',isComplex=true);
        end
    end

end
