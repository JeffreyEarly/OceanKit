classdef WVVerticalDamping < WVForcing
    % Vertical viscosity and diffusivity
    %
    % The damping is designed to mimic the VerticalScalarDiffusivity in
    % Oceananigans to allow for direct comparison between the models. This
    % is intended be used in combination with
    % WVHorizontalDamping. In general, you should be using the
    % [`WVAdaptiveDamping`](../WVAdaptiveDamping).
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
    % with viscosity, $$\nu$$, and diffusivity, $$\kappa$$. This should be combined with
    % [`WVHorizontalDamping`](../`WVHorizontalDamping`) for a complete closure. For help
    % choosing appropriate values, see the notes in
    % [`WVAdaptiveDamping`](../WVAdaptiveDamping).
    %
    % ### Usage
    %
    % Assuming there is a WVTransform instance wvt, to add this forcing,
    %
    % ```matlab
    % wvt.addForcing(WVVerticalDamping(wvt,nu=5e-4, kappa=1e-6));
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
    % - Declaration: WVVerticalDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % vertical viscosity
        %
        % - Topic: Properties
        nu

        % vertical diffusivity
        %
        % - Topic: Properties
        kappa

        % variable stratification factor
        %
        % - Topic: Properties
        dLnN2 = 0
    end

    methods
        function self = WVVerticalDamping(wvt,options)
            % initialize the WVVerticalDamping
            %
            % - Topic: Initialization
            % - Declaration: self = WVVerticalDamping(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter nu: vertical viscosity, default $$5 \cdot 10^{-4} \textrm{m}^2/\textrm{s}$$
            % - Parameter kappa: vertical diffusivity, default $$1 \cdot 10^{-6} \textrm{m}^2/\textrm{s}$$
            % - Returns self: a WVVerticalDamping instance
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
            self.dLnN2 = shiftdim(wvt.dLnN2,-2);

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
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVVerticalDamping
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVVerticalDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVVerticalDamping(wvtX2);
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