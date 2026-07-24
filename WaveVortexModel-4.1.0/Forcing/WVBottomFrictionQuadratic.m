classdef WVBottomFrictionQuadratic < WVForcing
    % Quadratic bottom friction
    %
    % Applies quadratic bottom friction to the flow, i.e., $$\frac{du}{dt} = -\frac{Cd}{dz}\lvert \mathbf{u} \rvert\mathbf{u}(x,y,-D)$$. Cd is unitless, and dz is (approximately) the size of the grid spacing at the bottom boundary.
    %
    % To compare with linear bottom friction where $$\frac{du}{dt} = -r \cdot u(x,y,-D)$$, note that $$- \frac{L_z}{dz} r = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert$$ and you will find a characteristic velocity $$\lvert \mathbf{u} \rvert$$ of about 11.5 cm/s for Cd=0.002 and r=1/(200 days). If $$C_d=0.001$$, then the damping time scale has to double to 1/(400 days) for and equivalent characteristic velocity.
    %
    % For barotropic QG we want the scaling to work out similarly, but now have $$-r = -\frac{C_d}{D}\lvert \mathbf{u} \rvert$$ where $$D$$ works out to be $$L_z$$. Thus, we will scale the barotropic QG quadratic bottom drag by 4000 m, to match typical oceanic scales.
    %
    % Using the notation that
    %
    % $$
    % |\mathbf{u}(x,y,-D)| = \sqrt{u^2(x,y,-D) + v^2(x,y,-D)}
    % $$
    %
    % is the magnitude of the total velocity at the bottom boundary, both nonhydrostatic and hydrostatic transforms linear bottom drag have the form
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= -\frac{C_d}{dz} |\mathbf{u}(x,y,-D)| u(x,y,-D) \\
    % \mathcal{S}_v &= -\frac{C_d}{dz} |\mathbf{u}(x,y,-D)| v(x,y,-D)  \\
    % \mathcal{S}_w &= 0 \\
    % \mathcal{S}_\eta &= 0
    % \end{align}
    % $$
    %
    % and for quasigeostrophic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_\textrm{qgpv} &= -\frac{C_dd}{dz} \left( \frac{\partial}{\partial x} \left( |\mathbf{u}(x,y,-D)| v(x,y,-D) \right) - \frac{\partial}{\partial y} \left( |\mathbf{u}(x,y,-D)| u(x,y,-D) \right) \right)
    % \end{align}
    % $$
    %
    % ### Usage
    %
    % Assuming there is a WVTransform instance wvt, to add this forcing,
    %
    % ```matlab
    % wvt.addForcing(WVBottomFrictionQuadratic(Cd=0.001));
    % ```
    %
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVBottomFrictionQuadratic < [WVForcing](/classes/wvforcing/)
    properties
        % non-dimensional quadratic drag coefficient
        %
        % - Topic: Properties
        Cd

        % $$\frac{Cd}{dz}$$ scaled quadratic drag coefficient with units $$m^{-1}$$
        %
        % - Topic: Properties
        cd
    end

    methods
        function self = WVBottomFrictionQuadratic(wvt,options)
            % initialize the WVBottomFrictionQuadratic
            %
            % - Topic: Initialization
            % - Declaration: self = WVBottomFrictionQuadratic(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter Cd: (optional) non-dimensional quadratic damping coefficient. Default is 0.001
            % - Returns frictionalForce: a WVBottomFrictionQuadratic instance
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
            propertyAnnotations(end+1) = CANumericProperty('Cd', {}, '','non-dimensional quadratic drag coefficient');
        end
    end

end