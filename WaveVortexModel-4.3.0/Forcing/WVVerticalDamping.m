classdef WVVerticalDamping < WVForcing
    % Apply vertical Laplacian viscosity and diffusivity.
    %
    % The damping is designed to mimic the VerticalScalarDiffusivity in
    % Oceananigans to allow for direct comparison between the models. This
    % applies to wave-bearing three-dimensional transforms and is intended
    % for use with
    % [`WVHorizontalDamping`](/classes/forcing/closures/wvhorizontaldamping/).
    % For an automatically scaled closure, use
    % [`WVAdaptiveDamping`](/classes/forcing/closures/wvadaptivedamping/).
    % 
    % The specific form of the forcing is given by 
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= \nu \frac{\partial^2 u}{\partial z^2} \\
    % \mathcal{S}_v &= \nu \frac{\partial^2 v }{\partial z^2} \\
    % \mathcal{S}_w &= \nu \frac{\partial^2 w}{\partial z^2} \\
    % \mathcal{S}_\eta &= \kappa \frac{\partial^2 \eta}{\partial z^2} - \kappa \frac{\partial}{\partial z} \ln N^2
    % \end{align}
    % $$
    %
    % Here $$\nu$$ is the vertical viscosity and $$\kappa$$ is the vertical
    % diffusivity. Combine this closure with
    % [`WVHorizontalDamping`](/classes/forcing/closures/wvhorizontaldamping/)
    % to damp horizontal gradients as well. For guidance on automatically
    % scaled coefficients, see
    % [`WVAdaptiveDamping`](/classes/forcing/closures/wvadaptivedamping/).
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.addForcing(WVVerticalDamping(wvt,nu=5e-4,kappa=1e-6));
    % ```
    %
    %
    % ### Notes
    %
    % This closure is evaluated in the spatial domain.
    %
    % For constant stratification, $$\partial_z \ln N^2=0$$ and the
    % stratification-gradient correction vanishes. The configured viscosity
    % and diffusivity are preserved when the forcing is copied to a
    % transform with a different resolution.
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    %
    % - Declaration: WVVerticalDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Vertical momentum viscosity in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % The constructor default is `5e-4`.
        %
        % - Topic: Properties
        nu

        % Vertical displacement diffusivity in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % The constructor default is `1e-6`.
        %
        % - Topic: Properties
        kappa

        % Precomputed vertical logarithmic stratification gradient.
        %
        % This Internal array is zero for constant stratification and has
        % units of inverse meters for variable stratification.
        %
        % - Topic: Properties
        dLnN2 = 0
    end

    methods
        function self = WVVerticalDamping(wvt,options)
            % Create vertical Laplacian damping for a transform.
            %
            % - Topic: Initialization
            % - Declaration: self = WVVerticalDamping(wvt,options)
            % - Parameter wvt: wave-bearing three-dimensional transform that owns the closure
            % - Parameter nu: optional vertical viscosity in square meters per second; default `5e-4`
            % - Parameter kappa: optional vertical diffusivity in square meters per second; default `1e-6`
            % - Returns self: vertical-damping closure owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.nu = 5e-4
                options.kappa = 1e-6
            end
            self@WVForcing(wvt,"vertical scalar diffusivity",WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.nu = options.nu;
            self.kappa = options.kappa;
            if isprop(wvt,"dLnN2")
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            end

            % construct the damping operator
            % [K,L,~] = self.wvt.kljGrid;
            % self.F_damp = -(K.^2 +L.^2);
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu = Fu +  self.nu*wvt.diffZF(wvt.u,n=2);
            Fv = Fv +  self.nu*wvt.diffZF(wvt.v,n=2);
            Feta = Feta + self.kappa * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu = Fu +  self.nu*wvt.diffZF(wvt.u,n=2);
            Fv = Fv +  self.nu*wvt.diffZF(wvt.v,n=2);
            Fw = Fw +  self.nu*wvt.diffZG(wvt.w,n=2);
            Feta = Feta + self.kappa * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Create equivalent vertical damping for another resolution.
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: vertical damping owned by `wvtX2`
            arguments
                self WVVerticalDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVVerticalDamping(wvtX2,nu=self.nu,kappa=self.kappa);
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
            vars = {"nu","kappa"};
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
            propertyAnnotations(end+1) = CANumericProperty('nu', {}, 'm2 s-1','viscosity');
            propertyAnnotations(end+1) = CANumericProperty('kappa', {}, 'm2 s-1','diffusivity');
        end
    end
end
