classdef WVBottomFrictionQuadratic < WVForcing
    % Apply quadratic drag at the bottom boundary.
    %
    % The dimensionless drag coefficient $$C_d$$ is divided by the bottom
    % quadrature weight for a three-dimensional transform:
    %
    % $$
    % c_d=\frac{C_d}{z_\mathrm{int}(1)}.
    % $$
    %
    % Barotropic QG uses a fixed 4000 m reference depth,
    % $$c_d=C_d/(4000\,\mathrm{m})$$. Comparing quadratic and linear drag
    % gives the characteristic relation $$L_z r=C_d\lvert\mathbf{u}\rvert$$.
    %
    % Using the notation that
    %
    % $$
    % |\mathbf{u}(x,y,-D)| = \sqrt{u^2(x,y,-D) + v^2(x,y,-D)}
    % $$
    %
    % is the magnitude of the total velocity at the bottom boundary. For
    % hydrostatic and nonhydrostatic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= -c_d |\mathbf{u}(x,y,-D)| u(x,y,-D) \\
    % \mathcal{S}_v &= -c_d |\mathbf{u}(x,y,-D)| v(x,y,-D)  \\
    % \mathcal{S}_w &= 0 \\
    % \mathcal{S}_\eta &= 0
    % \end{align}
    % $$
    %
    % and for quasigeostrophic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_\mathrm{qgpv} &= -c_d \left[ \partial_x \left( |\mathbf{u}|v \right) - \partial_y \left( |\mathbf{u}|u \right) \right]_{z=-D}
    % \end{align}
    % $$
    %
    % ### Example
    %
    % ```matlab
    % wvt = WVTransformConstantStratification([40e3,30e3,2e3],[8,6,5],N0=5.2e-3,latitude=45,isHydrostatic=true);
    % wvt.addForcing(WVBottomFrictionQuadratic(wvt,Cd=0.001));
    % ```
    %
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Inspect forcing or damping scales
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    %
    % - Declaration: WVBottomFrictionQuadratic < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Configured dimensionless quadratic drag coefficient.
        %
        % The constructor default is `1e-3`.
        %
        % - Topic: Properties
        Cd

        % Drag coefficient applied at the bottom in $$\mathrm{m^{-1}}$$.
        %
        % This is `Cd/z_int(1)` for a three-dimensional transform and
        % `Cd/4000` for a barotropic transform.
        %
        % - Topic: Properties
        cd
    end

    methods
        function contract = portableImplementationContract(self)
            % Return the paired portable implementation contract.
            %
            % - Topic: Forcing persistence
            % - Declaration: contract = portableImplementationContract(self)
            % - Returns contract: versioned data-only forcing contract
            % - Developer: true
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"Cd",double(self.Cd));
            contract = self.supportedPortableImplementationContract("WVBottomFrictionQuadratic",payload);
        end

        function self = WVBottomFrictionQuadratic(wvt,options)
            % Create quadratic bottom friction for a transform.
            %
            % See the [WVBottomFrictionQuadratic overview](/classes/forcing/wvbottomfrictionquadratic/)
            % for the governing equations, geometry-dependent scaling,
            % comparison with linear drag, and a usage example.
            %
            % - Topic: Initialization
            % - Declaration: self = WVBottomFrictionQuadratic(wvt,options)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Parameter Cd: optional dimensionless drag coefficient; default `1e-3`
            % - Returns self: quadratic bottom-friction forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
                options.Cd (1,1) double {mustBeNonnegative} = 1e-3 % https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"quadratic bottom friction",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.Cd = options.Cd;
            
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.cd = self.Cd/wvt.z_int(1);
            else
                % scaled by Lz, a typical ocean depth
                self.cd = self.Cd/4000;
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            ub = wvt.u(:,:,1);
            vb = wvt.v(:,:,1);
            cb = sqrt(ub.^2 + vb.^2);
            Fu(:,:,1) = Fu(:,:,1) - self.cd*ub.*cb;
            Fv(:,:,1) = Fv(:,:,1) - self.cd*vb.*cb;
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            ub = wvt.u(:,:,1);
            vb = wvt.v(:,:,1);
            cb = sqrt(ub.^2 + vb.^2);
            Fu(:,:,1) = Fu(:,:,1) - self.cd*ub.*cb;
            Fv(:,:,1) = Fv(:,:,1) - self.cd*vb.*cb;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            u_b = wvt.u(:,:,1);
            v_b = wvt.v(:,:,1);
            uv_mag = sqrt(u_b.^2 + v_b.^2);
            Fpv(:,:,1) = Fpv(:,:,1) - self.cd * (wvt.diffX(uv_mag.*v_b) - wvt.diffY(uv_mag.*u_b));
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFrictionQuadratic(wvtX2,Cd=self.Cd);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            vars = {'Cd'};
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
            propertyAnnotations(end+1) = CANumericProperty('Cd', {}, '1','dimensionless quadratic drag coefficient');
        end
    end

end
