classdef WVVerticalDiffusivity < WVForcing
    % Apply vertical diffusivity to the thermodynamic field.
    %
    % This forcing applies a vertical diffusivity with fixed $$\kappa_z$$
    % to the thermodynamic equation.
    % 
    % The specific form of the forcing is given by 
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= 0 \\
    % \mathcal{S}_v &= 0 \\
    % \mathcal{S}_w &= 0 \\
    % \mathcal{S}_\eta &= \kappa_z \frac{\partial^2 \eta}{\partial z^2} - \kappa_z \frac{\partial}{\partial z} \ln N^2
    % \end{align}
    % $$
    %
    % The horizontally uniform
    % $$-\kappa_z\partial_z\ln N^2$$ source projects onto the
    % mean-density-anomaly component. Set
    % `shouldForceMeanDensityAnomaly=false` to omit this source. The option
    % has no effect for constant stratification because the gradient is
    % zero, and it does not modify the wave modes.
    %
    % For `WVTransformStratifiedQG`, the QGPV definition
    %
    % $$
    % q=\partial_xv-\partial_yu-f\partial_z\eta
    % $$
    %
    % maps the displacement source
    % $$\mathcal{S}_\eta=\kappa_z\partial_{zz}\eta$$ to
    %
    % $$
    % \mathcal{S}_q=-f\partial_z\mathcal{S}_\eta
    % =-f\kappa_z\partial_{zzz}\eta.
    % $$
    %
    % Stratified QG contains only nonzero-horizontal-wavenumber geostrophic
    % modes, not a mean-density-anomaly component. Consequently,
    % `shouldForceMeanDensityAnomaly` does not alter its QGPV pathway.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.addForcing(WVVerticalDiffusivity(wvt,kappa_z=1e-6));
    % ```
    %
    % ### Notes
    %
    % This is currently implemented in the spatial domain. It applies to
    % wave-bearing three-dimensional transforms and has a separate QGPV
    % pathway for stratified QG. It is not compatible with barotropic QG,
    % which has no vertical structure.
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    %
    % - Declaration: WVVerticalDiffusivity < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Configured vertical diffusivity in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % The constructor default is `1e-5`.
        %
        % - Topic: Properties
        kappa_z

        % Whether to include the mean-density-anomaly source.
        %
        % The default is `true`. This controls the horizontally uniform
        % $$-\kappa_z\partial_z\ln N^2$$ source, which projects onto the
        % mean-density-anomaly component. It has no effect on wave modes,
        % for constant stratification, or for stratified QG.
        %
        % - Topic: Properties
        shouldForceMeanDensityAnomaly

        % Precomputed vertical logarithmic stratification gradient.
        %
        % This Internal value is zero when the correction is disabled or
        % stratification is constant, and otherwise has units of inverse
        % meters.
        %
        % - Topic: Properties
        dLnN2 = 0
    end

    methods
        function self = WVVerticalDiffusivity(wvt,options)
            % Create vertical diffusivity for a three-dimensional transform.
            %
            % - Topic: Initialization
            % - Declaration: self = WVVerticalDiffusivity(wvt,options)
            % - Parameter wvt: wave-bearing or stratified-QG transform that owns the forcing
            % - Parameter kappa_z: optional vertical diffusivity in square meters per second; default `1e-5`
            % - Parameter shouldForceMeanDensityAnomaly: optional flag controlling the horizontally uniform mean-density-anomaly source; default `true`
            % - Returns self: vertical-diffusivity forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.kappa_z double = 1e-5
                options.shouldForceMeanDensityAnomaly = true;
            end
            supportedTypes = ["HydrostaticSpatial","NonhydrostaticSpatial","PVSpatial"];
            if isa(wvt,"WVGeometryDoublyPeriodicBarotropic")
                supportedTypes = ["HydrostaticSpatial","NonhydrostaticSpatial"];
            end
            self@WVForcing(wvt,"vertical diffusivity",WVForcingType(supportedTypes));
            self.wvt = wvt;
            self.kappa_z = options.kappa_z;
            self.shouldForceMeanDensityAnomaly = options.shouldForceMeanDensityAnomaly;
            if self.shouldForceMeanDensityAnomaly && isprop(wvt,"dLnN2")
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Feta = Feta + self.kappa_z * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Feta = Feta + self.kappa_z * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            % Add the QGPV tendency induced by vertical displacement diffusion.
            %
            % This evaluates
            % $$-f\partial_z(\kappa_z\partial_{zz}\eta)$$.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: Fpv = addPotentialVorticitySpatialForcing(wvt,Fpv)
            % - Parameter wvt: stratified-QG transform evaluating the forcing
            % - Parameter Fpv: accumulated physical-space QGPV tendency
            % - Returns Fpv: QGPV tendency including vertical diffusivity
            % - Developer: true
            Fpv = Fpv - wvt.f * self.kappa_z * (wvt.diffZG(wvt.eta,n=3));
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Create equivalent vertical diffusivity for another resolution.
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: vertical diffusivity owned by `wvtX2`
            arguments
                self WVVerticalDiffusivity {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVVerticalDiffusivity(wvtX2,kappa_z=self.kappa_z,shouldForceMeanDensityAnomaly=self.shouldForceMeanDensityAnomaly);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            arguments
            end
            vars = {"kappa_z","shouldForceMeanDensityAnomaly"};
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
            propertyAnnotations(end+1) = CANumericProperty('kappa_z', {}, 'm2 s-1','vertical diffusivity');
            propertyAnnotations(end+1) = CANumericProperty('shouldForceMeanDensityAnomaly',{},'1', 'logical flag indicating whether to include the horizontally uniform mean-density-anomaly source');
        end
    end
end
