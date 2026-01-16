classdef WVBottomFrictionLinear < WVForcing
    % Linear bottom friction
    %
    % Applies linear bottom friction to the flow, i.e., $$\frac{du}{dt} = -r \cdot u(x,y,-D)$$. The parameter $$r$$ has units of $$s^{-1}$$ and thus can be set as an inverse time scale.
    %
    % The linear bottom friction is scaled such that we actually apply, $$\frac{du}{dt} = -\frac{L_z}{dz} r \cdot u(x,y,-D)$$ and the volume integrated effect of friction remains the same regardless of resolution. $$L_z$$ is the total domain depth and $$dz$$ is the spacing at the bottom grid point.
    %
    % To compare with quadratic bottom friction where $$\frac{du}{dt} = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert $$, note that $$- \frac{L_z}{dz} r = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert$$ and you will find a characteristic velocity $$\lvert\mathbf{u}\rvert$$ of about 10 cm/s for $$C_d=0.002$$.
    %
    % For both nonhydrostatic and hydrostatic transforms linear bottom drag
    %
    % $$
    % \begin{align}
    % \mathcal{S}_u &= -\frac{L_z}{dz} r \cdot u(x,y,-D) \\
    % \mathcal{S}_v &= -\frac{L_z}{dz} r \cdot v(x,y,-D)  \\
    % \mathcal{S}_w &= 0 \\
    % \mathcal{S}_\eta &= 0
    % \end{align}
    % $$
    %
    % and for quasigeostrophic transforms,
    %
    % $$
    % \begin{align}
    % \mathcal{S}_\textrm{qgpv} &= -\frac{L_z}{dz} r \cdot \zeta(x,y,-D)
    % \end{align}
    % $$
    %
    % where $$\zeta = \partial_x v - \partial_y u$$.
    %
    % ### Usage
    %
    % Assuming there is a WVTransform instance wvt, to add this forcing,
    %
    % ```matlab
    % wvt.addForcing(WVBottomFrictionLinear(r=1/(200*86400)));
    % ```
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: CAAnnotatedClass requirement
    %
    % - Declaration: WVBottomFrictionLinear < [WVForcing](/classes/wvforcing/)
    properties
        % bottom friction, $$s^{-1}$$
        %
        % - Topic: Properties
        r

        % scaled bottom friction, $$\frac{Lz}{dz} r$$ with units $$s^{-1}$$
        %
        % - Topic: Properties
        r_scaled
    end

    methods
        function self = WVBottomFrictionLinear(wvt,options)
            % initialize the WVBottomFrictionLinear
            %
            % - Topic: Initialization
            % - Declaration: self = WVBottomFrictionLinear(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter r: (optional) linear bottom friction, try 1/(200*86400)
            % - Returns frictionalForce: a WVBottomFrictionLinear instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.r (1,1) double {mustBeNonnegative} = 1/(200*86400) % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
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
            propertyAnnotations(end+1) = CANumericProperty('r', {}, 's^{-1}','bottom friction');
        end
    end

end