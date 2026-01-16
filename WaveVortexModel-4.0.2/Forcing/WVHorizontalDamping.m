classdef WVHorizontalDamping < WVForcing
    % Horizontal laplacian damping with viscosity and diffusivity
    %
    % The damping is a simple horizontal Laplacian, designed to mimic the
    % [HorizontalScalarDiffusivity in
    % Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/appendix/library/#Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity)
    % to allow for direct comparison between the models. This is intended to be used in combination with WVVerticalScalarDiffusivity. In
    % general, you should be using the
    % [`WVAdaptiveDamping`](../WVAdaptiveDamping).
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
    % which is just your standard Laplacian viscosity, $$\nu$$, and diffusivity, $$\kappa$$, in
    % the horizontal. This should be combined with
    % [`WVVerticalDamping`](../WVVerticalDamping) for a complete closure. For help
    % choosing appropriate values, see the notes in
    % [`WVAdaptiveDamping`](../WVAdaptiveDamping).
    %
    % ### Usage
    %
    % Assuming there is a WVTransform instance wvt, to add this forcing,
    %
    % ```matlab
    % wvt.addForcing(WVHorizontalDamping(wvt,nu=1e-4, kappa=1e-6));
    % ```
    %
    %
    % ### Notes
    %
    % This is currently implemented in the spatial domain, an is
    % thus highly un-optimized.
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVHorizontalDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % horizontal viscosity
        %
        % - Topic: Properties
        nu

        % horizontal diffusivity
        %
        % - Topic: Properties
        kappa
    end

    methods
        function self = WVHorizontalDamping(wvt,options)
            % initialize the WVHorizontalDamping
            %
            % - Topic: Initialization
            % - Declaration: self = WVHorizontalDamping(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter nu: horizontal viscosity, default $$1 \cdot 10^{-4} \textrm{m}^2/\textrm{s}$$
            % - Parameter kappa: horizontal diffusivity, default $$1 \cdot 10^{-6} \textrm{m}^2/\textrm{s}$$
            % - Returns self: a WVHorizontalDamping instance
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
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVHorizontalDamping
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVHorizontalDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVHorizontalDamping(wvtX2);
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
            propertyAnnotations(end+1) = CANumericProperty('nu', {}, 'm^2 s^{-1}','viscosity');
            propertyAnnotations(end+1) = CANumericProperty('kappa', {}, 'm^2 s^{-1}','diffusivity');
        end
    end
end