classdef IMSpectralWaveMomentumTests < matlab.unittest.TestCase
    % Independent differentiation must resolve the physical momentum balance.
    properties (TestParameter)
        profile = {"constant","exponential"}
        coordinateKind = {"z","wkb","density"}
        nEVP = {64,128}
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            testCase.addTeardown(@() path(originalPath));
        end
    end

    methods (Test)
        function waveMomentumSurvivesSpectralRefinement(testCase,profile,coordinateKind,nEVP)
            D = 1000; f = 1e-4; g = 9.81; nModes = 6;
            if profile == "constant", N2 = @(z) 1e-4*ones(size(z)); else, N2 = @(z) 1e-4*exp(2*z/700); end
            kappa = [2*pi/1e5 2*pi/1e3];
            surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=coordinateKind);
            collection = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-D 0],f0=f,g=g,surfaceBoundary=surface,nModes=nModes);
            for page = 1:numel(kappa)
                basis = collection.bases{page};
                testCase.verifyEqual(basis.modeNumber,1:nModes)
                testCase.verifyGreaterThan(basis.h,0)
                calculus = IMSolverSpectral(nEVP=257,coordinateKind="wkb").configuredForEVP(basis.evp);
                [z,weights,Dz] = calculus.nativeDifferentiationRule([-D 0]);
                G = basis.G(z); F = basis.F(z);
                h = basis.h(:).'; omega2 = f^2+g*h*kappa(page)^2;
                % Differentiate native G twice; never infer F_z from the EVP.
                directPressure = g*h.*basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./basis.normalizationFactors(basis.normalization);
                inertia = -omega2.*G; buoyancy = N2(z).*G;
                scale = weightedNorm(directPressure,weights)+weightedNorm(inertia,weights)+weightedNorm(buoyancy,weights);
                testCase.verifyLessThan(weightedNorm(directPressure+inertia+buoyancy,weights)./scale,2e-8)
                % Preserve the WVM physical-field diagnostic's existing bound.
                testCase.verifyLessThan(weightedNorm(g*Dz*F+inertia+buoyancy,weights)./scale,2e-7)
                testCase.verifyLessThan(abs(G(end,:)-F(end,:))./(max(abs(G),[],1)+max(abs(F),[],1)),1e-8)
                testCase.verifyLessThan(abs(G(1,:))./max(abs(G),[],1),1e-8)
                if profile == "constant"
                    exact = IMConstantStratificationSolution(N0=.01,zDomain=[-D 0],f0=f,g=g).internalModes(basis.evp,nModes=nModes);
                    testCase.verifyEqual(h,exact.h,RelTol=1e-8)
                    exactG = exact.G(z); exactF = exact.F(z);
                    orientation = sign(sum(weights.*G.*exactG,1));
                    testCase.verifyLessThan(weightedNorm(G-exactG.*orientation,weights)./weightedNorm(exactG,weights),1e-8)
                    testCase.verifyLessThan(weightedNorm(F-exactF.*orientation,weights)./weightedNorm(exactF,weights),1e-8)
                end
            end
        end
    end
end

function value = weightedNorm(values,weights)
value = sqrt(sum(weights.*abs(values).^2,1));
end
