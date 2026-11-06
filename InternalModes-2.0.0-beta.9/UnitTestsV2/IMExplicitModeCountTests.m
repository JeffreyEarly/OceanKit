classdef IMExplicitModeCountTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
        end
    end

    methods (Test)
        function scalarSolversRequireCountsForEveryInventory(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            wave = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            hydrostatic = IMInternalModes.hydrostaticFModes(background);
            for solver = IMExplicitModeCountTests.solvers()
                testCase.verifyError(@() solver{1}.solveEVP(wave),"IMSolver:MissingModeCount")
                testCase.verifyError(@() solver{1}.solveEVP(wave,nModes=[]),"IMSolver:MissingModeCount")
                testCase.verifyError(@() solver{1}.solveEVP(hydrostatic,branch="barotropic"),"IMSolver:MissingModeCount")
            end
        end

        function missingNumericalCountsFailBeforeSamplingStratification(testCase)
            sampleCalls = 0;
            problem = IMInternalModes.waveModesAtWavenumber(N2=@profile,zDomain=[-1000 0],k=2e-4,surfaceBoundary=IMBoundaryCondition.dirichlet());
            sampleCalls = 0;
            solvers = IMExplicitModeCountTests.solvers();
            for solver = solvers(1:2)
                testCase.verifyError(@() solver{1}.solveEVP(problem),"IMSolver:MissingModeCount")
                testCase.verifyError(@() solver{1}.solveWaveModesAtWavenumbers(2e-4,N2=@profile,zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.dirichlet()),"IMSolver:MissingModeCount")
            end
            testCase.verifyEqual(sampleCalls,0)
            function N2 = profile(z)
                sampleCalls = sampleCalls+1;
                N2 = 1e-4*ones(size(z));
            end
        end

        function explicitScalarZeroKeepsTheRequestedInventoryEmpty(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            wave = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            hydrostatic = IMInternalModes.hydrostaticFModes(background);
            for solver = IMExplicitModeCountTests.solvers()
                emptyWave = solver{1}.solveEVP(wave,nModes=0);
                emptyHydrostatic = solver{1}.solveEVP(hydrostatic,nModes=0);
                barotropic = solver{1}.solveEVP(hydrostatic,nModes=1,branch="barotropic");
                testCase.verifyEqual(emptyWave.nColumns,0)
                testCase.verifySize(emptyWave.G([-1000;-500;0]),[3 0])
                testCase.verifySize(emptyWave.Gz([-1000;-500;0]),[3 0])
                testCase.verifySize(emptyWave.Gzz([-1000;-500;0]),[3 0])
                testCase.verifyEqual(emptyHydrostatic.nColumns,0)
                testCase.verifySize(emptyHydrostatic.F([-1000;-500;0]),[3 0])
                testCase.verifyEqual(barotropic.nColumns,1)
                testCase.verifyEqual(barotropic.modeNumber,0)
            end
        end

        function bulkSolversRequireOneCountContractForEveryPage(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for solver = IMExplicitModeCountTests.solvers()
                for kappa = {[1e-4 2e-4],[1e-4 0],0}
                    testCase.verifyError(@() solver{1}.solveWaveModesAtWavenumbers(background,kappa{1}),"IMSolver:MissingModeCount")
                    testCase.verifyError(@() solver{1}.solveWaveModesAtWavenumbers(background,kappa{1},nModes=[]),"IMSolver:MissingModeCount")
                end
                testCase.verifyError(@() solver{1}.solveWaveModesAtWavenumbers(background,[],nModes=0),"IMSolver:EmptyWavenumberCollection")
                testCase.verifyError(@() solver{1}.solveWaveModesAtWavenumbers(background,[1e-4 0 2e-4],nModes=[2 3]),"IMSolver:InvalidWaveModeCounts")
                for count = [0 2]
                    [collection,diagnostics] = solver{1}.solveWaveModesAtWavenumbers(background,[0 0],nModes=count);
                    testCase.verifyEqual(collection.basisIndex,[1 1])
                    testCase.verifyEqual(numel(collection.bases),1)
                    testCase.verifyEqual(collection.bases{1}.nColumns,count)
                    testCase.verifyEqual(diagnostics.requestedCounts,[count count])
                end
                collection = solver{1}.solveWaveModesAtWavenumbers(background,[2e-4 0 1e-4 2e-4],nModes=[0 2 3 0]);
                testCase.verifyEqual(collection.basisIndex,[1 2 3 1])
                testCase.verifyEqual(cellfun(@(basis) basis.nColumns,collection.bases),[0 2 3])
                uniform = solver{1}.solveWaveModesAtWavenumbers(background,[1e-4 0],nModes=2);
                testCase.verifyEqual(cellfun(@(basis) basis.nColumns,uniform.bases),[2 2])
            end
        end

        function branchNamesAreValidatedEvenForEmptyRequests(testCase)
            sampleCalls = 0;
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            for solver = IMExplicitModeCountTests.solvers()
                testCase.verifyError(@() solver{1}.solveEVP(problem,branch="typo",nModes=0),"IMSolver:UnknownBranch")
                testCase.verifyError(@() solver{1}.solveEVP(problem,branch=["internal","typo"],nModes=[1 0]),"IMSolver:UnknownBranch")
                empty = solver{1}.solveEVP(problem,branch="surface",nModes=0);
                testCase.verifyEqual(empty.nColumns,0)
            end
            numericalProblem = IMInternalModes.waveModesAtWavenumber(N2=@profile,zDomain=[-1000 0],k=2e-4,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyError(@() IMSolverSpectral().solveEVP(numericalProblem,branch="typo",nModes=0),"IMSolver:UnknownBranch")
            testCase.verifyEqual(sampleCalls,0)
            function N2 = profile(z)
                sampleCalls = sampleCalls+1;
                N2 = 1e-4*ones(size(z));
            end
        end
    end

    methods (Static, Access = private)
        function values = solvers()
            values = {IMSolverSpectral(nEVP=24),IMSolverFiniteDifference(z=linspace(-1000,0,25).'),IMSolverAnalytical()};
        end
    end
end
