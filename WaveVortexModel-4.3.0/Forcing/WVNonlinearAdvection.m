classdef WVNonlinearAdvection < WVForcing
    % Add nonlinear advection to the model equations.
    %
    % The nonlinear terms are evaluated in physical space and added to the
    % momentum, thermodynamic, or quasigeostrophic potential-vorticity
    % (QGPV) equation appropriate to the transform.
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
    % \mathcal{S}_\mathrm{qgpv} &= - \left( u \partial_x q + v \partial_y q \right)
    % \end{align}
    % $$
    %
    % where $$q$$ is QGPV.
    %
    % ### Notes
    %
    % Every supported transform installs this forcing by default. A
    % nonlinear `WVModel` evaluates it automatically. Analytical linear
    % evolution does not evaluate nonlinear forcing, so the object does not
    % need to be removed when using linear evolution.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % nonlinearAdvection = wvt.forcingWithName("nonlinear advection");
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    %
    % - Declaration: WVNonlinearAdvection < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Precomputed vertical logarithmic stratification gradient.
        %
        % This Internal array is zero for constant stratification and has
        % units of inverse meters for variable stratification.
        %
        % - Topic: Properties
        dLnN2 = 0
    end

    methods
        function contract = portableImplementationContract(self)
            % Return the paired portable implementation contract.
            %
            % - Topic: Forcing internals
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only forcing contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority);
            contract = self.supportedPortableImplementationContract("WVNonlinearAdvection",payload);
        end

        function self = WVNonlinearAdvection(wvt)
            % Create nonlinear advection for a transform.
            %
            % See the [WVNonlinearAdvection overview](/classes/forcing/wvnonlinearadvection/)
            % for the hydrostatic, nonhydrostatic, and QGPV equations and
            % its role as the default forcing.
            %
            % - Topic: Initialization
            % - Declaration: self = WVNonlinearAdvection(wvt)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Returns self: nonlinear-advection forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"nonlinear advection",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.priority = 127;
            if isa(wvt,'WVStratification') && isprop(wvt,'dLnN2')
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
