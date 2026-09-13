classdef EtaTrueOperation < WVOperation
    % Compute displacement using the selected monotone no-motion profile.
    methods
        function self = EtaTrueOperation(wvt)
            arguments
                wvt WVTransform
            end
            outputVariables = WVVariableAnnotation('eta_true',{'x','y','z'},'m', 'true isopycnal deviation');
            self@WVOperation('eta_true',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            profile = self.profileForDiagnostics(wvt);
            materialHeight = profile.inverse(wvt.rho_total);
            varargout = {wvt.Z-materialHeight};
        end
    end

    methods (Static, Hidden)
        function profile = profileForDiagnostics(wvt)
            % Keep inversion and energetics on the same density representation.
            if wvt.shouldUseTrueNoMotionProfile
                density = wvt.rho_nm;
            else
                density = wvt.rho_nm0;
            end
            profile = WVNoMotionProfile(wvt.z,density);
        end
    end

    methods (Static)
        function y = fInverseBisection(f, x, yMin,yMax, tol)
            %FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
            % Taken from cumsum as part of chebfun.
            % chebfun/inv.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.

            a = yMin*ones(size(x));
            b = yMax*ones(size(x));
            c = (a + b)/2;

            while ( norm(b - a, inf) >= tol )
                vals = feval(f, c);
                % Bisection:
                I1 = ((vals-x) <= -tol);
                I2 = ((vals-x) >= tol);
                I3 = ~I1 & ~I2;
                a = I1.*c + I2.*a + I3.*c;
                b = I1.*b + I2.*c + I3.*c;
                c = (a+b)/2;
            end

            y = c;

        end
    end

end
