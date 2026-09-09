classdef IMVariableWaveCountTests < matlab.unittest.TestCase
    methods (Test)
        function distinctCountsPreserveExactModePrefixes(testCase)
            solver = IMSolverSpectral(nEVP=64,coordinateKind="wkb");
            N2 = @(z) 1e-4*(1+.2*exp(z/300));
            kappa = [3e-4 0 1e-4 3e-4 0 2e-4];
            [collection,costs] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],f0=1e-4,nModes=[2 5 2 3],nInertialModes=4);
            testCase.verifyEqual(collection.basisIndex,[1 2 3 1 2 4])
            testCase.verifyEqual(costs.solves.requestedColumnCount,[2;4;5;3])
            testCase.verifyEqual(costs.preparationCount,1)
            testCase.verifyEqual(costs.solveCount,4)
            z = linspace(-1000,0,65).';
            for p = 1:4
                basis = collection.bases{p};
                scalar = solver.solveEVP(basis.evp,nModes=costs.solves.requestedColumnCount(p));
                testCase.verifyEqual(basis.modeNumber,scalar.modeNumber)
                testCase.verifyEqual(basis.h,scalar.h,RelTol=1e-8)
                testCase.verifyEqual(basis.F(z),scalar.F(z),RelTol=2e-7,AbsTol=1e-9)
                testCase.verifyEqual(basis.G(z),scalar.G(z),RelTol=2e-7,AbsTol=1e-9)
            end
            testCase.verifyEqual(collection.evaluate(z,variable="G",pages=[4 1]),repmat(collection.bases{1}.G(z),1,1,2),AbsTol=1e-10)
            testCase.verifyError(@()collection.evaluate(z,variable="G",pages=[1 3]),"IMBasisCollection:HeterogeneousColumns")
        end

        function countMapValidationPrecedesSolves(testCase)
            solver = IMSolverSpectral(nEVP=32);
            args = {"N2",@(z)1e-4*ones(size(z)),"zDomain",[-1000 0]};
            testCase.verifyError(@()solver.solveWaveModesAtWavenumbers([1e-4 2e-4 1e-4],args{:},nModes=[2 3]),"IMSolver:InvalidWaveModeCounts")
            testCase.verifyError(@()solver.solveWaveModesAtWavenumbers([1e-4 2e-4 1e-4],args{:},nModes=[2 3 4]),"IMSolver:ConflictingWaveModeCounts")
            [collection,costs] = solver.solveWaveModesAtWavenumbers([0 0],args{:},nModes=[],nInertialModes=3);
            testCase.verifyEqual(costs.solves.requestedColumnCount,3)
            testCase.verifyEqual(collection.basisIndex,[1 1])
        end
    end
end
