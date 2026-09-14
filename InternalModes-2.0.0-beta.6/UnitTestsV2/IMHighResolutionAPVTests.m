classdef IMHighResolutionAPVTests < matlab.unittest.TestCase
    % Refine the solve without changing the retained physical mode band.
    methods (Test)
        function highResolutionRetainsAnalyticalEigendepths(testCase)
            N0 = 5.2e-3;
            b = 1300;
            D = 4000;
            g = 9.81;
            N2 = @(z) N0^2*exp(2*z/b);
            integratedN2 = N0^2*b*(1-exp(-2*D/b))/2;
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=N2,zDomain=[-D 0],g=g,g0=-integratedN2, ...
                gd=integratedN2,surfaceBoundary="freeSurface");
            exact = IMExponentialStratificationSolution( ...
                N0=N0,b=b,zDomain=[-D 0],g=g).internalModes(evp,nModes=84);
            z = linspace(-D,0,1601).';
            exactF = exact.F(z);
            for nEVP = [783 1551]
                solved = IMSolverSpectral(nEVP=nEVP).solveEVP(evp,nModes=84);
                testCase.verifyEqual(solved.modeNumber,exact.modeNumber)
                testCase.verifyLessThan(max(abs(solved.h-exact.h)./abs(exact.h)),1e-6)
                F = solved.F(z);
                orientation = sign(sum(F.*exactF,1));
                % Compare each scientific mode independently; no mixing or
                % fitting is allowed to conceal a changed resolved basis.
                testCase.verifyLessThan(max(vecnorm(F.*orientation-exactF)./vecnorm(exactF)),1e-5)
            end
        end
    end
end
