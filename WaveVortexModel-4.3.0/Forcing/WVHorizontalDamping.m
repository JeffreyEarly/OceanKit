classdef WVHorizontalDamping < WVForcing
    % Apply horizontal Laplacian viscosity and diffusivity.
    %
    % The damping is a simple horizontal Laplacian, designed to mimic the
    % [HorizontalScalarDiffusivity in
    % Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/appendix/library/#Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity)
    % to allow direct comparison between the models. It applies to
    % wave-bearing three-dimensional transforms and is intended for use with
    % [`WVVerticalDamping`](/classes/forcing/closures/wvverticaldamping/).
    % For an automatically scaled closure, use
    % [`WVAdaptiveDamping`](/classes/forcing/closures/wvadaptivedamping/).
    % 
    % The specific form of the forcing is given by 
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= \nu \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right) u \\
    % \mathcal{S}_v &= \nu \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)  v \\
    % \mathcal{S}_w &= \nu \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)  w \\
    % \mathcal{S}_\eta &= \kappa \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)  \eta
    % \end{align}
    % $$
    %
    % These are horizontal Laplacian viscosity, $$\nu$$, and diffusivity,
    % $$\kappa$$. Combine this closure with
    % [`WVVerticalDamping`](/classes/forcing/closures/wvverticaldamping/) to
    % damp vertical gradients as well. For guidance on automatically scaled
    % coefficients, see
    % [`WVAdaptiveDamping`](/classes/forcing/closures/wvadaptivedamping/).
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.addForcing(WVHorizontalDamping(wvt,nu=1e-4,kappa=1e-6));
    % ```
    %
    %
    % ### Notes
    %
    % This closure is evaluated in the spatial domain.
    %
    % The configured viscosity and diffusivity are preserved when the
    % forcing is copied to a transform with a different resolution.
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    %
    % - Declaration: WVHorizontalDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Horizontal momentum viscosity in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % The constructor default is `1e-4`.
        %
        % - Topic: Properties
        nu

        % Horizontal displacement diffusivity in $$\mathrm{m^{2}\,s^{-1}}$$.
        %
        % The constructor default is `1e-6`.
        %
        % - Topic: Properties
        kappa
    end

    methods
        function self = WVHorizontalDamping(wvt,options)
            % Create horizontal Laplacian damping for a transform.
            %
            % - Topic: Initialization
            % - Declaration: self = WVHorizontalDamping(wvt,options)
            % - Parameter wvt: wave-bearing three-dimensional transform that owns the closure
            % - Parameter nu: optional horizontal viscosity in square meters per second; default `1e-4`
            % - Parameter kappa: optional horizontal diffusivity in square meters per second; default `1e-6`
            % - Returns self: horizontal-damping closure owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.nu = 1e-4
                options.kappa = 1e-6
            end
            self@WVForcing(wvt,"horizontal scalar diffusivity",WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.nu = options.nu;
            self.kappa = options.kappa;

            % construct the damping operator
            % [K,L,~] = self.wvt.kljGrid;
            % self.F_damp = -(K.^2 +L.^2);
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu = Fu + self.nu*(wvt.diffX(wvt.u,n=2) + wvt.diffY(wvt.u,n=2));
            Fv = Fv +  self.nu*(wvt.diffX(wvt.v,n=2) + wvt.diffY(wvt.v,n=2));
            Feta = Feta + self.kappa*(wvt.diffX(wvt.eta,n=2) + wvt.diffY(wvt.eta,n=2));
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu = Fu + self.nu*(wvt.diffX(wvt.u,n=2) + wvt.diffY(wvt.u,n=2));
            Fv = Fv +  self.nu*(wvt.diffX(wvt.v,n=2) + wvt.diffY(wvt.v,n=2));
            Fw = Fw +  self.nu*(wvt.diffX(wvt.w,n=2) + wvt.diffY(wvt.w,n=2));
            Feta = Feta + self.kappa*(wvt.diffX(wvt.eta,n=2) + wvt.diffY(wvt.eta,n=2));
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Create equivalent horizontal damping for another resolution.
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: horizontal damping owned by `wvtX2`
            arguments
                self WVHorizontalDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVHorizontalDamping(wvtX2,nu=self.nu,kappa=self.kappa);
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
