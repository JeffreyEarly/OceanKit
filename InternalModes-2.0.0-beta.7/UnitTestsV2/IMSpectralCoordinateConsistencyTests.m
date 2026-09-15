classdef IMSpectralCoordinateConsistencyTests < matlab.unittest.TestCase
    methods (Test)
        function stretchedMapsAndDerivativesAgreeWithExactExponential(testCase)
            D = 1000; b = 700; N0 = .01;
            N2 = @(z)N0^2*exp(2*z/b);
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-D 0],k=1e-3,f0=1e-4);
            for kind = ["wkb" "density"]
                if kind == "wkb", scale = N0; rate = 1/b; else, scale = N0^2; rate = 2/b; end
                solver = IMSolverSpectral(nEVP=65,coordinateKind=kind).configuredForEVP(evp);
                zQuery = linspace(-D,0,111);
                expectedX = scale/rate*(exp(rate*zQuery)-exp(-rate*D));
                actualX = solver.xOfZ(zQuery);
                testCase.verifyEqual(actualX,expectedX,AbsTol=5e-13)
                testCase.verifyEqual(solver.zOfX(actualX),zQuery,AbsTol=2e-11)
                testCase.verifyEqual(solver.qzReference,rate*scale*exp(rate*solver.zReference),RelTol=2e-11)
                [z,~,Dz] = solver.nativeDifferentiationRule([-D 0]);
                [~,~,Dzz] = solver.nativeDifferentiationRule([-D 0],2);
                value = exp(rate*z);
                testCase.verifyEqual(Dz*value,rate*value,RelTol=1e-8)
                testCase.verifyEqual(Dzz*value,rate^2*value,RelTol=1e-7)
                testCase.verifyTrue(all(isnan(solver.xOfZ([-D-1 1]))))
                testCase.verifyTrue(all(isnan(solver.zOfX([-1 2*expectedX(end)]))))
            end
        end

        function reconfigurationRebuildsMapsWithoutChangingTheOriginal(testCase)
            firstEVP = IMInternalModes.hydrostaticGModes(N2=@(z)1e-4*exp(2*z/700),zDomain=[-1000 0]);
            secondEVP = IMInternalModes.hydrostaticGModes(N2=@(z)4e-4*ones(size(z)),zDomain=[-1000 0]);
            first = IMSolverSpectral(nEVP=33,coordinateKind="wkb").configuredForEVP(firstEVP);
            z = [-1000;-500;0];
            original = first.xOfZ(z);
            second = first.configuredForEVP(secondEVP);
            testCase.verifyEqual(first.xOfZ(z),original)
            testCase.verifyEqual(second.xOfZ(z),.02*(z+1000),AbsTol=1e-12)
            testCase.verifyEqual(second.qzReference,zeros(size(second.qzReference)),AbsTol=1e-14)
            testCase.verifyGreaterThan(norm(first.zNative-second.zNative),1)
        end
    end
end
