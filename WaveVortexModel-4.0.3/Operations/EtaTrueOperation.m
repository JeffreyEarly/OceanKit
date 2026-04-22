classdef EtaTrueOperation < WVOperation
% Computes the true vertical displacement, eta.
%
% 2026-01-15
%
% Using moments to compute rho_nm, and then a high order spline to compute
% the functional version. The following code gives,
%
% int_vol_avg = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3)/wvt.Lz;
% 
% data = wvt.rho_nm;
% spline_nm = BSpline(K,BSpline.knotPointsForDataPoints(wvt.z,K=K));
% spline_nm.x_mean = mean(data);
% spline_nm.x_std = std(data);
% Z = BSpline.Spline( wvt.z, spline_nm.t_knot, spline_nm.K );
% spline_nm.m = Z\((data - spline_nm.x_mean)/spline_nm.x_std);
% 
% rho_e = wvt.rho_total - shiftdim(wvt.rho_nm,-2);
% int_vol_avg(rho_e)/max(abs(rho_e(:)))
% int_vol_avg(wvt.eta_true)/max(abs(wvt.eta_true(:)))
%
% which produces relative errors of 3e-6 for rho_e and 1e-3 for eta. This
% corresponds to a volume average of 13cm for eta, about 10x better than
% reported in July 2025 below.
%
% 2025-07-14
% Note: this will recover rho_e with
%   rho_e = wvt.rhoFunction(wvt.Z-wvt.eta_true)-wvt.rhoFunction(wvt.Z);
% for which
%   drho = wvt.rho_e - rho_e;
%   sqrt(mean(drho(:).^2))
% gives 1e-13 for a realistic simulation.
%
% Using
%   int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);
% I find that
%   int_vol(eta_true)
% returns ~4000 (or an average of 1 m displacement).
%
% This should be compared with
%   int_vol(shiftdim(wvt.N2,-2).*wvt.eta)/wvt.N2(end)
% which returns ~-65 or O(100) m.
%
% This all assume an adiabatic simulation---this is really just a measure
% of the MDA.

    properties (GetAccess=public, SetAccess=protected)
        spline_nm
        Z
    end

    methods

        function self = EtaTrueOperation(wvt)
            arguments
                wvt 
            end
            outputVariables(1) = WVVariableAnnotation('eta_true',{'x','y','z'},'m', 'true isopycnal deviation');
            self@WVOperation('eta_true',outputVariables,@disp);

            K=min(wvt.Nz,8);
            data = wvt.rho0 - wvt.rho_nm0;
            self.spline_nm = BSpline(K,BSpline.knotPointsForDataPoints(wvt.z,K=K));
            self.spline_nm.x_mean = mean(data);
            self.spline_nm.x_std = std(data);
            self.Z = BSpline.Spline( wvt.z, self.spline_nm.t_knot, self.spline_nm.K );
            self.spline_nm.m = self.Z\((data - self.spline_nm.x_mean)/self.spline_nm.x_std);
        end

        function varargout = compute(self,wvt,varargin)
            data = wvt.rho0 - wvt.rho_nm;
            self.spline_nm.m = self.Z\((data - self.spline_nm.x_mean)/self.spline_nm.x_std);

            rho_total = (wvt.rhoFunction(wvt.Z) - wvt.rho0) + wvt.rho_e ;

            zMinusEta = EtaTrueOperation.fInverseBisection(self.spline_nm,-rho_total(:),-wvt.Lz,0,1e-12);
            zMinusEta = reshape(zMinusEta,size(wvt.X));
            eta_true = wvt.Z - zMinusEta;
            varargout = {eta_true};
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