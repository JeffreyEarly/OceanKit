classdef WVFixedAmplitudeForcing < WVForcing
    % Fixed amplitude forcing at the natural frequency of each mode
    %
    % The fixed-amplitude forcing maintains the amplitude of the wave or geostrophic features, while allowing energy and enstrophy to flux from the feature.
    %
    % As a simple example, one can set an internal wave mode with amplitude 1 cm/s, and that mode will continue to oscillate and maintain its amplitude. The wave will participate in all the nonlinear dynamics, but its amplitude will be maintained/restored at each time step.
    %
    %
    % There are several different ways to write this style of forcing mathematically. The equations of motion, written in the spectral domain, take the following form
    %
    % $$
    % \frac{\partial}{\partial t} A^{klj} = \sum_i F_i^{klj}
    % $$
    %
    % where $$F_i$$ are the different forces applied. The transform computes the spatial forcing (which includes nonlinear advection), the spectral forcing, followed by the spectral amplitude forcing. The `WVFixedAmplitudeForcing` is a spectral amplitude forcing and is thus comptued last. This forcing thus simply adds back the flux the from the spatial and spectral forcing, so that $$\frac{\partial}{\partial t} A^{klj} =0$$ for the modes in question.
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
    % - Every mode that is used in `WVFixedAmplitudeForcing` essentially removes a degree-of-freedom from the model, as that mode is no longer free to fully evolve. Thus when you pass the forcing wave-vortex coefficients, e.g. `A0`, it does not fix the amplitude of coefficients that are small to avoid removing degrees-of-freedom.
    % - One must also be careful not to forcing in the damping region. If you have some sort of small scale damping enabled, you probably do not want to be forcing at those smallest scales.
    %
    % ### Usage
    %
    % To setup a geostrophic mean flow,
    %
    % ```matlab
    % % initialize a transform
    % wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=@(z) N0*N0*exp(2*z/L_gm),latitude=33);
    %
    % % set a geostrophic mode, with no flow at the bottom boundary
    % wvt.setGeostrophicModes(k=0,l=5,j=1,phi=0,u=u0);
    % wvt.setGeostrophicModes(k=0,l=5,j=0,phi=0,u=max(max(wvt.u(:,:,1))));
    %
    % % pass the vortex coefficients to the forcing
    % force = WVFixedAmplitudeForcing(wvt,name="geostrophic-mean-flow");
    % force.setGeostrophicForcingCoefficients(wvt.A0);
    % wvt.addForcing(force);
    %
    % ```
    %
    % In practice you can initialize the flow in any way you want with any arbitrary structure, and then pass those coefficients to the forcing. The `WVFixedAmplitudeForcing` looks for coefficients that are small and ignores those.
    %
    % - Topic: Initialization
    % - Topic: Setting the forcing
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVFixedAmplitudeForcing < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % indices of modes in the `A0` matrix to fix
        %
        % - Topic: Properties
        A0_indices (:,1) uint64 = []

        % indices of modes in the `Ap` matrix to fix
        %
        % - Topic: Properties
        Ap_indices (:,1) uint64 = []

        % indices of modes in the `Am` matrix to fix
        %
        % - Topic: Properties
        Am_indices (:,1) uint64 = []

        % amplitudes of the fixed modes in the `A0` matrix
        %
        % - Topic: Properties
        A0bar (:,1) double = []

        % amplitudes of the fixed modes in the `Ap` matrix
        %
        % - Topic: Properties
        Apbar (:,1) double = []

        % amplitudes of the fixed modes in the `Am` matrix
        %
        % - Topic: Properties
        Ambar (:,1) double = []
    end

    methods
        function self = WVFixedAmplitudeForcing(wvt,options)
            % initialize the WVFixedAmplitudeForcing
            %
            % You must pass the instance of the WVTransform to be used and
            % you must also specify a unique name for the forcing.
            %
            % - Topic: Initialization
            % - Declaration: self = WVFixedAmplitudeForcing(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter name: (required) name of this forcing
            % - Parameter Apbar: (optional) amplitude of Ap matrix to fix
            % - Parameter Ambar: (optional) amplitude of Am matrix to fix
            % - Parameter A0bar: (optional) amplitude of A0 matrix to fix
            % - Parameter Ap_indices: (optional) index of coefficient in Ap matrix to fix
            % - Parameter Am_indices: (optional) index of coefficient in Am matrix to fix
            % - Parameter A0_indices: (optional) index of coefficient in A0 matrix to fix
            % - Returns self: a WVFixedAmplitudeForcing instance
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
            % set the amplitude to fix for the wave part of the flow
            %
            % This function will automatically remove modes set in the
            % damping region of the WVAdapativeDamping forcing, if present.
            %
            % - Topic: Setting the forcing
            % - Declaration:  setWaveForcingCoefficients(Apbar,Ambar,options)
            % - Parameter Apbar: Ap fixed amplitude
            % - Parameter Ambar: Am fixed amplitude
            % - Parameter MAp: (optional) forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes. Default is MAp = abs(Apbar) > 1e-6*max(abs(Apbar(:)))
            % - Parameter MAm: (optional) forcing mask, Am. 1s at the forced modes, 0s at the unforced modes. Default is MAm = abs(Ambar) > 1e-6*max(abs(Ambar(:)))
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
            % set amplitude to fix for the geostrophic part of the flow
            %
            % This function will automatically remove modes set in the
            % damping region of the WVAdapativeDamping forcing, if present.
            %
            % - Topic: Setting the forcing
            % - Declaration: setGeostrophicForcingCoefficients(A0bar,options)
            % - Parameter A0bar: A0 fixed amplitude
            % - Parameter MA0: (optional) forcing mask, A0. 1s at the forced modes, 0s at the unforced modes. Default is MA0 = abs(A0bar) > 1e-6*max(abs(A0bar(:)))
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
            % sets a narrow waveband of geostrophic forcing for forced-dissipative modeling
            %
            % to be moved to a subclass
            %
            % - Topic: Setting the forcing
            % - Declaration: setNarrowBandGeostrophicForcing(options)
            % - Parameter r: A0 fixed amplitude
            arguments
                self WVFixedAmplitudeForcing {mustBeNonempty}
                options.r (1,1) double
                options.k_r (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*2
                options.k_f (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*20
                options.j_f (1,1) double = 1
                options.u_rms (1,1) double = 0.2 % set the *total* energy (not just kinetic) equal to 0.5*u_rms^2
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end

            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                % the idea is to set the energy at the sea-surface and
                % so we need to know the relative amplitude of this
                % mode at the surface.
                F = self.wvt.FinvMatrix;
                surfaceMag = 1/F(end,options.j_f+1);
                sbRatio = abs(F(end,options.j_f+1)/F(1,options.j_f+1));
                % sbRatio = 1; % should we change the damping scale? Or no?
                h = self.wvt.h_0(options.j_f+1);
                magicNumber = 2.25;
            else
                surfaceMag = 1;
                sbRatio = 1;
                h = self.wvt.h;
                magicNumber = 0.0225;
            end


            if isfield(options,"r")
                k_r = options.r/(magicNumber*options.u_rms);
            else
                r = magicNumber*sbRatio*options.u_rms*options.k_r; % 1/s bracket [0.02 0.025]
                % fprintf('1/r is %.1f days, switching to %.1f days\n',1/(self.r*86400),1/(r*86400));
                k_r = options.k_r;
            end
            k_f = options.k_f;
            j_f = options.j_f;
            wvt = self.wvt;

            % smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(self.r),1,'first');
            % fprintf('(k_r=%.2f dk, k_f=%d dk, k_nu=%d dk.\n',k_r/wvt.dk,round(k_f/wvt.dk),round(self.k_damp/wvt.dk));
            % fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(k_f/(wvt.k(2)-wvt.k(1))));


            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            MA0 = zeros(wvt.spectralMatrixSize);
            MA0(wvt.Kh > k_f-deltaK/2 & wvt.Kh < k_f+deltaK/2 & wvt.J == j_f) = 1;

            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                u_rms = surfaceMag * options.u_rms;

                m = 3/2; % We don't really know what this number is.
                kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
                model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
                model_inverse = @(k) kappa_epsilon * k.^(-5/3);
                model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
                model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

                [~,~,wvt.A0] = wvt.geostrophicComponent.randomAmplitudesWithSpectrum(A0Spectrum= @(k,j) model_spectrum(k),shouldOnlyRandomizeOrientations=1);

                if strcmp(options.initialPV,'narrow-band')
                    wvt.A0 = MA0 .* wvt.A0;
                else
                    if isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
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
                    fprintf("surface u_rms: %.2g cm/s\n",100*u_rms_surface);
                    fprintf("surface energy, %g.\n",KE+PE);
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/h);
                end
            end
            self.setGeostrophicForcingCoefficients(MA0 .* wvt.A0,MA0=MA0);
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
            Abar = zeros(self.wvt.spectralMatrixSize);
            Abar(self.Ap_indices) = self.Apbar;
            [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
            options.Apbar = self.Apbar;
            options.Ap_indices = find(AbarX2);

            Abar = zeros(self.wvt.spectralMatrixSize);
            Abar(self.Am_indices) = self.Ambar;
            [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
            options.Ambar = self.Ambar;
            options.Am_indices = find(AbarX2);

            Abar = zeros(self.wvt.spectralMatrixSize);
            Abar(self.A0_indices) = self.A0bar;
            [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
            options.A0bar = self.A0bar;
            options.A0_indices = find(AbarX2);

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
            propertyAnnotations(end+1) = CADimensionProperty('Ap_indices', '','indices into the Ap matrix');
            propertyAnnotations(end+1) = CANumericProperty('Apbar', {'Ap_indices'}, '','Ap mean value',isComplex=true);
            propertyAnnotations(end+1) = CADimensionProperty('Am_indices', '','indices into the Am matrix');
            propertyAnnotations(end+1) = CANumericProperty('Ambar', {'Am_indices'}, '','Am mean value',isComplex=true);
            propertyAnnotations(end+1) = CADimensionProperty('A0_indices', '','indices into the A0 matrix');
            propertyAnnotations(end+1) = CANumericProperty('A0bar', {'A0_indices'}, '','A0 mean value',isComplex=true);
        end
    end

end