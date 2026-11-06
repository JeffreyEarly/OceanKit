classdef IMVariableWaveCountTests < matlab.unittest.TestCase
    methods (Test)
        function distinctCountsPreserveExactModePrefixes(testCase)
            solver = IMSolverSpectral(nEVP=64,coordinateKind="wkb");
            N2 = @(z) 1e-4*(1+.2*exp(z/300));
            kappa = [3e-4 0 1e-4 3e-4 0 2e-4];
            [collection,costs] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],rotationRate=1e-4/(2*sind(33)),latitude=33,nModes=[2 4 5 2 4 3],surfaceBoundary=IMBoundaryCondition.dirichlet());
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

        function repeatedCountsReuseTheLargestSolveAndPreserveRequests(testCase)
            background = IMBackground.exponentialStratification(N0=.01,b=700,zDomain=[-1000 0],surfaceCondition="freeSurface");
            kappa = [2e-4 0 2e-4 1e-4 2e-4 0 2e-4];
            counts = [2 0 5 3 2 4 0];
            z = linspace(-1000,0,67).';
            for solver = {IMSolverSpectral(nEVP=64),IMSolverFiniteDifference(z=z),IMSolverAnalytical()}
                [collection,costs] = solver{1}.solveWaveModesAtWavenumbers(background,kappa,nModes=counts);
                testCase.verifyEqual(costs.solveCount,3)
                testCase.verifyEqual(costs.solves.requestedColumnCount,[5;4;3])
                testCase.verifyEqual(collection.basisIndex,[1 2 3 4 1 5 6])
                testCase.verifyEqual(numel(collection.bases),6)
                for page = 1:numel(kappa)
                    basis = collection.bases{collection.basisIndex(page)};
                    record = basis.constructionRecord.selection;
                    testCase.verifyEqual(basis.nColumns,counts(page))
                    testCase.verifyEqual(record.request.nModes,counts(page))
                    largest = max(counts(kappa==kappa(page)));
                    testCase.verifyEqual(record.request.sharedSolve.nModes,largest)
                    testCase.verifyEqual(record.request.sharedSolve.request.nModes,largest)
                    testCase.verifyGreaterThanOrEqual(record.originalColumnCount,largest)
                    testCase.verifyEqual(record.activeSourceColumns,1:counts(page))
                    source = solver{1}.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=kappa(page)),nModes=largest);
                    expected = source.F(z);
                    testCase.verifyEqual(basis.F(z),expected(:,1:counts(page)),RelTol=2e-8,AbsTol=2e-10)
                end
                short = collection.bases{1}; full = collection.bases{3};
                fullFactors = full.normalizationFactors();
                testCase.verifyEqual(short.normalizationFactors(),fullFactors(1:2),RelTol=2e-12)
            end
        end

        function countMapValidationPrecedesSolves(testCase)
            solver = IMSolverSpectral(nEVP=32);
            args = {"N2",@(z)1e-4*ones(size(z)),"zDomain",[-1000 0]};
            testCase.verifyError(@()solver.solveWaveModesAtWavenumbers([1e-4 2e-4 1e-4],args{:},rotationRate=0,nModes=[2 3],surfaceBoundary=IMBoundaryCondition.dirichlet()),"IMSolver:InvalidWaveModeCounts")
            [collection,costs] = solver.solveWaveModesAtWavenumbers([0 0],args{:},rotationRate=0,nModes=3,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyEqual(costs.solves.requestedColumnCount,3)
            testCase.verifyEqual(collection.basisIndex,[1 1])
        end
    end
end
