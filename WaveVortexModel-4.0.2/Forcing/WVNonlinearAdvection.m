classdef WVNonlinearAdvection < WVForcing
    % The advective flux, $$\mathbf{u}\cdot \nabla \mathbf{u}$$ and $$\mathbf{u}\cdot \nabla \eta$$
    %
    % The nonlinear advection forcing adds the nonlinear terms to the momentum and thermodynamic equation.
    %
    % The nonlinear terms are all computed in the spatial domain.
    %
    % For nonhydrostatic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= - \left( u \partial_x u + v \partial_y u + w \partial_z u \right) \\
    % \mathcal{S}_v &= - \left( u \partial_x v + v \partial_y v + w \partial_z v \right) \\
    % \mathcal{S}_w &= - \left(  u \partial_x w + v \partial_y w + w \partial_z w \right) \\
    % \mathcal{S}_\eta &= - \left( u \partial_x \eta + v \partial_y \eta  + w \left(\partial_z \eta +\eta \partial_z \ln N^2 \right) \right)
    % \end{align}
    % $$
    %
    % for hydrostatic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= - \left( u \partial_x u + v \partial_y u + w \partial_z u \right) \\
    % \mathcal{S}_v &= - \left( u \partial_x v + v \partial_y v + w \partial_z v \right) \\
    % \mathcal{S}_\eta &= - \left( u \partial_x \eta + v \partial_y \eta  + w \left(\partial_z \eta +\eta \partial_z \ln N^2 \right) \right)
    % \end{align}
    % $$
    %
    % and for quasigeostrophic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_\textrm{qgpv} &= - \left( u \partial_x q + v \partial_y q \right)
    % \end{align}
    % $$
    %
    % where $$q$$ is the qgpv.
    %
    % ### Notes
    %
    % This is the only forcing added to the transforms by default. You must explicitly remove it if you want to consider linear flows.
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVNonlinearAdvection < [WVForcing](/classes/wvforcing/)
    properties
        % variable stratification factor
        %
        % - Topic: Properties
        dLnN2 = 0
    end

    methods
        function self = WVNonlinearAdvection(wvt)
            % initialize the WVNonlinearAdvection nonlinear flux
            %
            % - Topic: Initialization
            % - Declaration: self = WVNonlinearAdvection(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVNonlinearAdvection instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"nonlinear advection",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.priority = 127;
            if isa(wvt,'WVStratification')
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            end
        end
        
        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu = Fu - (wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u));
            Fv = Fv - (wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v));
            Feta = Feta - (wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2));
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu = Fu - (wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u));
            Fv = Fv - (wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v));
            Fw = Fw - (wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w));
            Feta = Feta - (wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2));
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv = Fpv - (wvt.u.*wvt.diffX(wvt.qgpv) + wvt.v.*wvt.diffY(wvt.qgpv));
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVNonlinearAdvection(wvtX2);
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