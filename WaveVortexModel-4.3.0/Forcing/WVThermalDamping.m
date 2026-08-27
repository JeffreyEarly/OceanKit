classdef WVThermalDamping < WVForcing
    % Apply large-scale thermal damping to QGPV.
    %
    % For each deformation scale $$L_r$$, the implementation adds
    %
    % $$
    % \mathcal{S}_q=\frac{\alpha}{L_r^2}\psi.
    % $$
    %
    % This follows the large-scale thermal-damping formulation considered
    % by [Scott and Dritschel](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/halting-scale-and-energy-equilibration-in-twodimensional-quasigeostrophic-turbulence/BD0CAFC9019691ADC9B18A95D15445F9).
    %
    % ### Notes
    %
    % This forcing is compatible only with stratified and barotropic QG
    % transforms through their physical-space QGPV forcing stage.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformBarotropicQG([40e3,30e3],[8,6],h=0.8,latitude=45);
    % wvt.addForcing(WVThermalDamping(wvt,alpha=1/(200*86400)));
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Inspect forcing or damping scales
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Declaration: WVThermalDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Configured thermal-damping rate in $$\mathrm{s^{-1}}$$.
        %
        % The constructor default is `1/(200*86400)`.
        %
        % - Topic: Properties
        alpha

        % Deformation-scaled damping coefficient in $$\mathrm{s^{-1}\,m^{-2}}$$.
        %
        % This is `alpha/wvt.Lr2` and has the shape required by the QG
        % streamfunction field.
        %
        % - Topic: Properties
        alpha_scaled
    end

    methods
        function self = WVThermalDamping(wvt,options)
            % Create thermal damping for a QG transform.
            %
            % - Topic: Initialization
            % - Declaration: self = WVThermalDamping(wvt,options)
            % - Parameter wvt: stratified or barotropic QG transform that owns the forcing
            % - Parameter alpha: optional damping rate in inverse seconds; default `1/(200*86400)`
            % - Returns self: thermal-damping forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.alpha (1,1) double {mustBeNonnegative} = 1/(200*86400)
            end
            self@WVForcing(wvt,"thermal damping",WVForcingType("PVSpatial"));
            self.alpha = options.alpha;
            self.alpha_scaled = self.alpha/wvt.Lr2;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv = Fpv + self.alpha_scaled * wvt.psi;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVThermalDamping(wvtX2,alpha=self.alpha);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            vars = {'alpha'};
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
            propertyAnnotations(end+1) = CANumericProperty('alpha', {}, 's-1','thermal damping coefficient');
        end
    end

end
