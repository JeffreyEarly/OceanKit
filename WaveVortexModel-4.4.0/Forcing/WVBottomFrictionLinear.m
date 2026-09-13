classdef WVBottomFrictionLinear < WVForcing
    % Apply linear drag at the bottom boundary.
    %
    % The parameter $$r$$ is an inverse time scale in $$\mathrm{s^{-1}}$$.
    % For a three-dimensional transform, the bottom tendency is scaled by
    % the bottom quadrature weight `z_int(1)` so its vertically integrated
    % effect does not change with vertical resolution:
    %
    % $$
    % r_\mathrm{scaled}=\frac{L_z}{z_\mathrm{int}(1)}r.
    % $$
    %
    % A barotropic transform has no vertical quadrature and uses
    % $$r_\mathrm{scaled}=r$$.
    %
    % Comparing this with quadratic drag gives the characteristic relation
    % $$L_z r=C_d\lvert\mathbf{u}\rvert$$.
    %
    % For both nonhydrostatic and hydrostatic transforms linear bottom drag
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= -r_\mathrm{scaled} u(x,y,-D) \\
    % \mathcal{S}_v &= -r_\mathrm{scaled} v(x,y,-D)  \\
    % \mathcal{S}_w &= 0 \\
    % \mathcal{S}_\eta &= 0
    % \end{align}
    % $$
    %
    % and for quasigeostrophic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_\mathrm{qgpv} &= -r_\mathrm{scaled}\zeta(x,y,-D)
    % \end{align}
    % $$
    %
    % where $$\zeta = \partial_x v - \partial_y u$$.
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.addForcing(WVBottomFrictionLinear(wvt,r=1/(200*86400)));
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Inspect forcing or damping scales
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    %
    % - Declaration: WVBottomFrictionLinear < [WVForcing](/classes/forcing/wvforcing/)
    properties (SetAccess = private)
        % Configured linear drag rate in $$\mathrm{s^{-1}}$$.
        %
        % The constructor default is `1/(200*86400)`, corresponding to a
        % 200-day time scale.
        %
        % - Topic: Properties
        r

        % Drag rate applied at the bottom grid point in $$\mathrm{s^{-1}}$$.
        %
        % This is `r*Lz/z_int(1)` for a three-dimensional transform and
        % `r` for a barotropic transform.
        %
        % - Topic: Properties
        r_scaled
    end

    methods
        function contract = portableImplementationContract(self)
            % Return the paired portable implementation contract.
            %
            % - Topic: Forcing persistence
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only forcing contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"r",double(self.r));
            contract = self.supportedPortableImplementationContract("WVBottomFrictionLinear",payload);
        end

        function self = WVBottomFrictionLinear(wvt,options)
            % Create linear bottom friction for a transform.
            %
            % See the [WVBottomFrictionLinear overview](/classes/forcing/wvbottomfrictionlinear/)
            % for the governing equations, resolution scaling, comparison
            % with quadratic drag, and a usage example.
            %
            % - Topic: Initialization
            % - Declaration: self = WVBottomFrictionLinear(wvt,options)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Parameter r: optional drag rate in inverse seconds; default `1/(200*86400)`
            % - Returns self: linear bottom-friction forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.r (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 1/(200*86400) % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"linear bottom friction",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.r = options.r;
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.r_scaled = self.r * wvt.Lz / wvt.z_int(1);
            else
                self.r_scaled = self.r;
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r_scaled*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r_scaled*wvt.v(:,:,1);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r_scaled*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r_scaled*wvt.v(:,:,1);
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv(:,:,1) = Fpv(:,:,1) - self.r_scaled * wvt.zeta_z(:,:,1);
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFrictionLinear(wvtX2,r=self.r);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            vars = {'r'};
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
            propertyAnnotations(end+1) = CANumericProperty('r', {}, 's-1','bottom friction');
        end
    end

end
