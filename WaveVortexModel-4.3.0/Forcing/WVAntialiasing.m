classdef WVAntialiasing < WVForcing
    % Apply explicit spectral antialias filtering for diagnostics.
    %
    % This closure removes coefficient tendencies in aliased modes and sets
    % those coefficients to zero after every integration step. The
    % horizontal mask uses the transform's quadratic-aliasing rule. Vertical
    % modes with `j >= Nj` are discarded; `Nj` defaults to
    % `floor(2*wvt.Nj/3)`.
    %
    % Transform-level antialiasing is enabled by default and is more
    % efficient because discarded modes are never computed. Explicit
    % antialiasing is intended for measuring the filter's effect on energy
    % and potential enstrophy. Construct the transform with
    % `shouldAntialias=false` before adding this closure. It is compatible
    % with wave-bearing and QG transforms.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true,shouldAntialias=false);
    % wvt.addForcing(WVAntialiasing(wvt));
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Inspect forcing or damping scales
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    % - Declaration: WVAntialiasing < [WVForcing](/classes/forcing/wvforcing/)
    properties (GetAccess=public, SetAccess=protected)
        % Number of retained vertical modes.
        %
        % This value is preserved when the forcing is converted to another
        % resolution or restored from an annotated NetCDF file.
        %
        % - Topic: Properties
        Nj

        % Logical-shape spectral mask of discarded coefficients.
        %
        % This array has `wvt.spectralMatrixSize`; nonzero entries are
        % removed from coefficient tendencies and amplitudes.
        %
        % - Topic: Properties
        M
    end

    methods
        function self = WVAntialiasing(wvt,options)
            % Create explicit antialias filtering for a transform.
            %
            % The transform must have been constructed with
            % `shouldAntialias=false`.
            %
            % - Declaration: self = WVAntialiasing(wvt,options)
            % - Parameter wvt: transform that owns and evaluates the closure
            % - Parameter Nj: optional number of retained vertical modes; modes with `j >= Nj` are discarded; default `floor(2*wvt.Nj/3)`
            % - Returns self: explicit-antialiasing closure owned by `wvt`
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
            self.Nj = options.Nj;
            self.M(wvt.J > (options.Nj-1)) = 1;
        end

        function effectiveHorizontalGridResolution = effectiveHorizontalGridResolution(self)
            % Return the shortest fully retained horizontal wavelength.
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: effectiveHorizontalGridResolution = effectiveHorizontalGridResolution()
            % - Returns effectiveHorizontalGridResolution: effective horizontal resolution in meters
            arguments
                self WVAntialiasing
            end
            effectiveHorizontalGridResolution = pi/max(max(abs(self.wvt.L(~self.M)),abs(self.wvt.K(~self.M))));
        end

        function j_max = effectiveJMax(self)
            % Return the highest retained vertical-mode number.
            %
            % This dimensionless value accounts for the explicit vertical
            % mask and is appropriate when constructing damping operators.
            %
            % - Topic: Properties
            % - Declaration: j_max = effectiveJMax()
            % - Returns j_max: highest retained vertical-mode number
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
            force = WVAntialiasing(wvtX2,Nj=self.Nj);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            vars = {"Nj"};
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
            propertyAnnotations(end+1) = CANumericProperty('Nj', {}, '1','number of retained vertical modes used to construct the filter');
        end
    end
end
