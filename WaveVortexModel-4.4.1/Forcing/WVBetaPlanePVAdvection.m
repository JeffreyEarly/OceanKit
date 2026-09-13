classdef WVBetaPlanePVAdvection < WVForcing
    % Add beta-plane advection to the balanced QGPV tendency.
    %
    % On a beta plane, material conservation of total quasigeostrophic
    % potential vorticity gives
    %
    % $$
    % \frac{D}{Dt}(q+\beta y)=0,
    % $$
    %
    % so the right-hand-side tendency contributed by this forcing is
    %
    % $$
    % \left.\frac{\partial q}{\partial t}\right|_\beta=-\beta v_g.
    % $$
    %
    % QG transforms evaluate this expression directly in physical QGPV
    % space. Wave-bearing transforms apply the equivalent spectral tendency
    % only to the geostrophic `A0` coefficients; `Ap`, `Am`, inertial modes,
    % and mean-density-anomaly modes receive no direct beta tendency. This
    % retains beta advection for the balanced flow but is not a full
    % beta-plane treatment of internal-wave dynamics: wave frequencies and
    % structures continue to use the transform's constant Coriolis
    % parameter.
    %
    % ```matlab
    % wvt.addForcing(WVBetaPlanePVAdvection(wvt));
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    % - Declaration: WVBetaPlanePVAdvection < [WVForcing](/classes/forcing/wvforcing/)
    properties
        % Spectral multiplier mapping `A0` to its beta-plane tendency.
        %
        % This Internal array is $$-\beta V_{A0}$$. It is zero on modes
        % without geostrophic meridional velocity, including the horizontal
        % mean and mean-density-anomaly modes.
        %
        % - Topic: Forcing internals
        % - Developer: true
        betaA0
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
            contract = self.supportedPortableImplementationContract("WVBetaPlanePVAdvection",payload);
        end

        function self = WVBetaPlanePVAdvection(wvt)
            % Create beta-plane QGPV advection for a transform.
            %
            % - Topic: Create the forcing
            % - Declaration: self = WVBetaPlanePVAdvection(wvt)
            % - Parameter wvt: transform receiving the balanced beta-plane tendency
            % - Returns self: beta-plane QGPV-advection forcing owned by `wvt`
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"beta-plane advection of qgpv",WVForcingType(["Spectral" "PVSpatial"]));
            self.betaA0 = -wvt.beta * wvt.VA0;
            self.betaA0(1,1,1) = 0;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            % Add $$-\beta v_g$$ to a physical-space QGPV tendency.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: Fpv = addPotentialVorticitySpatialForcing(wvt,Fpv)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter Fpv: accumulated physical-space QGPV tendency
            % - Returns Fpv: QGPV tendency including beta-plane advection
            % - Developer: true
            Fpv = Fpv - wvt.beta * wvt.v;
        end

        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            % Add the balanced beta-plane tendency in spectral space.
            %
            % `Fp` and `Fm` are returned unchanged. The multiplier is zero
            % outside the geostrophic `A0` modes.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Fp,Fm,F0] = addSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: wave-bearing transform evaluating the forcing
            % - Parameter Fp: accumulated positive-frequency tendency
            % - Parameter Fm: accumulated negative-frequency tendency
            % - Parameter F0: accumulated zero-frequency tendency
            % - Returns Fp: unchanged positive-frequency tendency
            % - Returns Fm: unchanged negative-frequency tendency
            % - Returns F0: zero-frequency tendency including beta-plane advection
            % - Developer: true
            F0 = F0 + self.betaA0 .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            % Create beta-plane advection for another resolution.
            %
            % - Topic: Convert forcing resolution
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: beta-plane forcing owned by `wvtX2`
            % - Developer: true
            force = WVBetaPlanePVAdvection(wvtX2);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end

end
