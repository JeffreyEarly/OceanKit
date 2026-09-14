classdef IMBulkWaveModesTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            originalPath = path;
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
            testCase.addTeardown(@() path(originalPath));
        end
    end
    methods (Test)
        function repeatedUnorderedWavenumbersAgreeWithScalarSolves(testCase)
            N2 = @(z) 1e-4*(1+0.2*exp(z/300));
            kappa = [3e-4 0 1e-4 3e-4 0];
            for coordinateKind = ["z","wkb","density"]
                solver = IMSolverSpectral(nEVP=64,coordinateKind=coordinateKind);
                [collection,diagnostics] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],f0=1e-4,nModes=5,nInertialModes=3);
                testCase.verifyEqual(collection.kappa,kappa)
                testCase.verifyEqual(collection.basisIndex,[1 2 3 1 2])
                testCase.verifyEqual(diagnostics.solvedKappa,[3e-4 0 1e-4])
                testCase.verifyEqual(diagnostics.solveCount,3)
                testCase.verifyEqual(diagnostics.reusedPageCount,2)
                testCase.verifyEqual(diagnostics.solves.requestedColumnCount,[5;3;5])
                testCase.verifyEqual(diagnostics.solveAccuracy,"unverified")
                testCase.verifyGreaterThanOrEqual(diagnostics.setupSeconds,0)
                testCase.verifyGreaterThanOrEqual(diagnostics.assemblySeconds,diagnostics.sharedAssemblySeconds)
                testCase.verifyGreaterThan(diagnostics.sharedAssemblyStorageBytes,0)
                z = linspace(-1000,0,31).';
                for i = 1:3
                    basis = collection.bases{i};
                    evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],f0=1e-4,k=diagnostics.solvedKappa(i));
                    scalar = solver.solveEVP(evp,nModes=diagnostics.solves.requestedColumnCount(i));
                    testCase.verifyEqual(basis.evp.parameters.k,diagnostics.solvedKappa(i))
                    testCase.verifyEqual(basis.modeNumber,scalar.modeNumber)
                    testCase.verifyEqual(basis.normalization,scalar.normalization)
                    testCase.verifyEqual(basis.eigenvalues,scalar.eigenvalues,RelTol=2e-8,AbsTol=1e-11)
                    testCase.verifyEqual(basis.G(z),scalar.G(z),RelTol=2e-7,AbsTol=1e-9)
                    testCase.verifyEqual(basis.F(z),scalar.F(z),RelTol=2e-7,AbsTol=1e-9)
                    testCase.verifyEqual(basis.uz(z),scalar.uz(z),RelTol=2e-7,AbsTol=1e-9)
                end
                testCase.verifyEqual(collection.evaluate(z,variable="G",pages=[4 1]),repmat(collection.bases{1}.G(z),1,1,2),AbsTol=1e-12)
                expectedF = cat(3,collection.bases{1}.F(z),collection.bases{3}.F(z),collection.bases{1}.F(z));
                testCase.verifyEqual(collection.evaluate(z,variable="F",pages=[4 3 1],sampleChunkSize=7,pageChunkSize=2),expectedF,AbsTol=1e-10,RelTol=1e-12)
            end
        end
        function freeSurfaceAndNondefaultBoundaryRowsRemainUnchanged(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            solver = IMSolverSpectral(nEVP=64);
            surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            bottom = IMBoundaryCondition(a=1,b=0.1,c=0,d=0);
            [collection,~] = solver.solveWaveModesAtWavenumbers([1e-4 4e-4],N2=N2,zDomain=[-1000 0],f0=1e-4,surfaceBoundary=surface,bottomBoundary=bottom,nModes=4);
            base = IMInternalModes.waveModesAtWavenumber(k=0,N2=N2,zDomain=[-1000 0],f0=1e-4,surfaceBoundary=surface,bottomBoundary=bottom);
            prepared = solver.configuredForEVP(base);
            [A0,B0,samples] = base.assembleConfigured(prepared);
            A2 = prepared.physicalDerivativeMatrix(0);
            endpoints = [prepared.boundaryIndex("surface"),prepared.boundaryIndex("bottom")];
            A2(endpoints,:) = 0;
            for i = 1:2
                basis = collection.bases{i};
                evp = basis.evp;
                [A,B] = evp.assemble(solver);
                expectedA = A0+evp.parameters.k^2*A2;
                testCase.verifyEqual(expectedA,A,AbsTol=2e-13)
                testCase.verifyEqual(expectedA(endpoints,:),A0(endpoints,:),AbsTol=0)
                testCase.verifyEqual(B,B0,AbsTol=0)
                samples.q = evp.parameters.k^2*ones(size(samples.q));
                testCase.verifyEqual(evp.preparedModeSelectionDiagnostics(samples,A),evp.modeSelectionDiagnostics(prepared,A))
                scalar = solver.solveEVP(evp,nModes=4);
                testCase.verifyEqual(basis.eigenvalues,scalar.eigenvalues,RelTol=2e-8,AbsTol=1e-11)
                z = [-1000;0];
                testCase.verifyEqual(basis.F(z),scalar.F(z),RelTol=2e-7,AbsTol=1e-8)
                testCase.verifyEqual(basis.G(z),scalar.G(z),RelTol=2e-7,AbsTol=1e-8)
            end
        end
        function sharedStratificationIsActuallySampledOnceForZCoordinates(testCase)
            sampleCalls = 0;
            solver = IMSolverSpectral(nEVP=40);
            [collection,diagnostics] = solver.solveWaveModesAtWavenumbers([2e-4 1e-4 2e-4 3e-4],N2=@profile,zDomain=[-1000 0],nModes=3);
            testCase.verifyEqual(sampleCalls,1)
            testCase.verifyEqual(numel(collection.bases),3)
            testCase.verifyEqual(diagnostics.preparationCount,1)
            function values = profile(z)
                sampleCalls = sampleCalls+1;
                values = 1e-4*ones(size(z));
            end
        end
        function customSolversCannotSilentlySharePreparation(testCase)
            solver = IMBulkUnsupportedSpectralSolver();
            testCase.verifyError(@() solver.solveWaveModesAtWavenumbers([1e-4 2e-4],N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],nModes=3),"IMSolver:UnsupportedBulkSolver")
        end
        function finiteDifferenceBulkMatchesItsOwnScalarDiscretization(testCase)
            solver = IMSolverFiniteDifference(z=linspace(-1000,0,41).');
            N2 = @(z) 1e-4*(1+0.1*exp(z/300));
            [collection,diagnostics] = solver.solveWaveModesAtWavenumbers([2e-4 1e-4 2e-4],N2=N2,zDomain=[-1000 0],nModes=4);
            testCase.verifyEqual(diagnostics.solveCount,2)
            for basis = collection.bases
                scalar = solver.solveEVP(basis{1}.evp,nModes=4);
                testCase.verifyEqual(basis{1}.modeNumber,scalar.modeNumber)
                testCase.verifyEqual(basis{1}.eigenvalues,scalar.eigenvalues,RelTol=1e-10,AbsTol=1e-12)
                testCase.verifyEqual(basis{1}.G(solver.zNative),scalar.G(solver.zNative),RelTol=1e-10,AbsTol=1e-10)
            end
        end
        function explicitCountsAndInputIdentityAreEnforced(testCase)
            solver = IMSolverSpectral(nEVP=16);
            N2 = @(z) 1e-4*ones(size(z));
            testCase.verifyError(@() solver.solveWaveModesAtWavenumbers([1e-4 0],N2=N2,zDomain=[-1000 0],nModes=3),"IMSolver:MissingInertialModeCount")
            testCase.verifyError(@() solver.solveWaveModesAtWavenumbers([],N2=N2,zDomain=[-1000 0],nModes=3),"IMSolver:EmptyWavenumberCollection")
            testCase.verifyError(@() solver.solveWaveModesAtWavenumbers(1e-4,N2=N2,zDomain=[-1000 0],nModes=30),"IMSolver:InsufficientWaveModes")
            testCase.verifyError(@() solver.solveWaveModesAtWavenumbers(0,N2=N2,zDomain=[-1000 0],nModes=3,nInertialModes=30),"IMSolver:InsufficientWaveModes")
            [collection,diagnostics] = solver.solveWaveModesAtWavenumbers(1e-4,N2=N2,zDomain=[-1000 0],nModes=3);
            testCase.verifyEqual(numel(collection.bases),1)
            testCase.verifyEqual(diagnostics.solveCount,1)
        end
        function negativeWaveLabelsAndScalarZeroLabelsArePreserved(testCase)
            solver = IMSolverSpectral(nEVP=48);
            N2 = @(z) ones(size(z));
            surface = IMBoundaryCondition(a=2,b=1);
            bottom = IMBoundaryCondition.neumann();
            [collection,~] = solver.solveWaveModesAtWavenumbers([0.2 0.1],N2=N2,zDomain=[-1 0],surfaceBoundary=surface,bottomBoundary=bottom,nModes=4);
            for i = 1:2
                basis = collection.bases{i};
                scalar = solver.solveEVP(basis.evp,nModes=4);
                testCase.verifyTrue(any(basis.modeNumber < 0))
                testCase.verifyEqual(basis.modeNumber,scalar.modeNumber)
                testCase.verifyEqual(basis.eigenvalues,scalar.eigenvalues,RelTol=1e-8,AbsTol=1e-10)
            end
            zeroEVP = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1 0]);
            scalar = solver.solveEVP(zeroEVP,nModes=4);
            testCase.verifyEqual(scalar.modeNumber,0:3)
            testCase.verifyEqual(scalar.eigenvalues(1),0)
            testCase.verifyTrue(all(isfinite(scalar.F([-1;0])),"all"))
        end
    end
end
