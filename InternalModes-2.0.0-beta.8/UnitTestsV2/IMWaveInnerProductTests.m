classdef IMWaveInnerProductTests < matlab.unittest.TestCase
    methods (Test)
        function waveEnergyPreservesSturmLiouvilleProduct(testCase)
            for surface = [IMBoundaryCondition.dirichlet(),IMBoundaryCondition(a=0,b=1,c=1,d=0)]
                for slope = [0 .5]
                    N2 = @(z)1e-4*(1+slope*z/1000);
                    f0 = .003; g = 9.81;
                    evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,surfaceBoundary=surface);
                    basis = IMSolverSpectral(nEVP=96,coordinateKind="wkb").solveEVP(evp,nModes=6);
                    basis = basis.withNormalization("waveEnergy");
                    [z,weights] = basis.solver.nativeQuadratureRule([-1000 0]);
                    G = basis.G(z);
                    % Unit wave energy gives D*lambda_j*<G_j,G_j> = 1.
                    gram = G'*(weights.*((N2(z)-f0^2)/g).*G);
                    if surface.c ~= 0, gram = gram+G(end,:)'*G(end,:); end
                    testCase.verifyEqual(gram,diag(basis.h/diff(basis.zDomain)),AbsTol=2e-8)
                    testCase.verifyEqual(basis.gramMatrix(variable="G"),gram,AbsTol=2e-10)
                    testCase.verifyEqual(basis.majorantGramMatrix(variable="G"),gram,AbsTol=2e-10)
                end
            end
        end

        function fixedFrequencyUsesItsOwnCanonicalWeight(testCase)
            N2 = @(z)1e-4*(1+.5*z/1000);
            omega = .003; g = 9.81;
            evp = IMInternalModes.waveModesAtFrequency(rotationRate=0,N2=N2,zDomain=[-1000 0],omega=omega,g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
            basis = IMSolverSpectral(nEVP=96,coordinateKind="wkb").solveEVP(evp,nModes=5);
            basis = basis.withNormalization("rmsF");
            [z,weights] = basis.solver.nativeQuadratureRule([-1000 0]);
            G = basis.G(z);
            gram = G'*(weights.*((N2(z)-omega^2)/g).*G);
            testCase.verifyEqual(gram,diag(basis.h/diff(basis.zDomain)),AbsTol=2e-8)
        end

        function analyticalRotatingModesUseTheSameCorrectedNormalization(testCase)
            N0 = .01; f0 = .003; g = 9.81;
            N2 = @(z)N0^2*ones(size(z));
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0,zDomain=[-1000 0]);
            basis = solution.internalModes(evp,nModes=4);
            basis = basis.withNormalization("waveEnergy");
            z = linspace(-1000,0,1001).';
            G = basis.G(z);
            testCase.verifyEqual(trapz(z,((N0^2-f0^2)/g)*G.^2),basis.h/diff(basis.zDomain),AbsTol=1e-10)
        end

        function solvedGUsesCanonicalRAndMajorantIsPositive(testCase)
            N2 = @(z)1e-4*ones(size(z));
            evp = IMInternalModes(rotationRate=0,N2=N2,zDomain=[-1000 0],r=@(z,~)z+500);
            signed = evp.innerProduct("G");
            positive = evp.majorantInnerProduct("G");
            z = [-1000;-600;-400;0];
            testCase.verifyEqual(signed.interiorWeight(z,struct()),z+500)
            testCase.verifyEqual(positive.interiorWeight(z,struct()),abs(z+500))
        end
    end
end
