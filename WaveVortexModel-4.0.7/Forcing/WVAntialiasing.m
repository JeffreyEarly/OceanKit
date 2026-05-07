classdef WVAntialiasing < WVForcing
    % Antialiasing filter
    %
    % This forcing removes (sets to zero) energy in the largest 1/3 modes
    % to prevent quadratic aliasing. This is the correct de-aliasing for
    % the horizontal fourier modes, but is not exactly correct for the
    % vertical modes in variable stratification. You can thus manually set
    % `Nj`, otherwise it will default to `options.Nj = floor(2*wvt.Nj/3);`.
    %
    % **Very important** You should almost never use this forcing, as
    % de-aliasing is built-in at the transform level and enabled by
    % default. For performance reasons is far more optimal to simply never
    % compute the de-aliased modes. The purpose of this forcing to allow
    % direct measurement of the effect of the de-aliasing on energy and
    % potential enstrophy. It is quite slow and thus we recommend it be
    % used for diagnostic purposes only.
    %
    % ### Usage
    %
    % This is likely to change in the future, but at the moment several of
    % the transforms have a function that will make a new transform in the
    % identical state, but with the antialiasing filter explicitly added.
    %
    % ```matlab
    % wvtAA = wvt.waveVortexTransformWithExplicitAntialiasing();
    % ```
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Topic: Initializing
    % - Declaration: WVAntialiasing < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % spectral matrix that multiplies Ap,Am,A0 to zero out the aliased modes
        %
        % - Topic: Properties
        M
    end

    methods
        function self = WVAntialiasing(wvt,options)
            % initialize the WVAntialiasing
            %
            % - Declaration: nlFlux = WVNonlinearFlux(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter Nj: (optional) vertical mode above which energy will be set to zero.
            % - Returns self: a WVAntialiasing instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.Nj
            end
            self@WVForcing(wvt,"antialias filter",WVForcingType(["Spectral","PVSpectral"]));
            self.priority = 127;
            if wvt.shouldAntialias
                error("WVAntialiasing:AntialiasingNotSupported","Antialiasing is not supported for a transform that aliases at the transform level.");
            end
            self.wvt = wvt;
            self.isClosure = true;
            Aklz = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.k_dft,wvt.l_dft,wvt.Nj);
            self.M = wvt.transformFromDFTGridToWVGrid(Aklz);
            if ~isfield(options,"Nj")
                options.Nj = floor(2*wvt.Nj/3);
            end
            self.M(wvt.J > (options.Nj-1)) = 1;
        end

        function effectiveHorizontalGridResolution = effectiveHorizontalGridResolution(self)
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
                self WVAntialiasing
            end
            effectiveHorizontalGridResolution = pi/max(max(abs(self.wvt.L(~self.M)),abs(self.wvt.K(~self.M))));
        end

        function j_max = effectiveJMax(self)
            %returns the effective highest vertical mode
            %
            % The effective highest vertical modeis the highest fully resolved
            % mode in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: flag = effectiveJMax(other)
            % - Returns effectiveJMax: double
            arguments
                self WVAntialiasing
            end
            j_max = max(abs(self.wvt.J(~self.M)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            Fp = Fp - self.M .* Fp;
            Fm = Fm - self.M .* Fm;
            F0 = F0 - self.M .* F0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            F0 = F0 - self.M .* F0;
        end

        function [Ap, Am, A0] = setSpectralAmplitude(self, wvt, Ap, Am, A0)
            Ap = (~self.M) .* Ap;
            Am = (~self.M) .* Am;
            A0 = (~self.M) .* A0;
        end



        function A0 = setPotentialVorticitySpectralAmplitude(self, wvt, A0)
            arguments
                self WVForcing
                wvt WVTransform
                A0 
            end
            A0 = (~self.M) .* A0;
        end

        function F0 = setPotentialVorticitySpectralForcing(self, wvt, F0)
            arguments
                self WVForcing
                wvt WVTransform
                F0
            end
            F0 = F0 - self.M .* F0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVAntialiasing(wvtX2);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
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