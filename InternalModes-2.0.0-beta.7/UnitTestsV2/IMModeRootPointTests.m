classdef IMModeRootPointTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function dirichletGridUsesNextSelectedMode(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisWithoutAuxiliary = solver.solveEVP(evp,nModes=4);
            basisWithAuxiliary = solver.solveEVP(evp,nModes=5);

            zAutomatic = basisWithoutAuxiliary.pointsFromModeRoots();
            zRetained = basisWithAuxiliary.pointsFromModeRoots(nModes=4);
            expected = linspace(-1,0,6).';

            testCase.verifyEqual(zAutomatic,expected,AbsTol=1e-10)
            testCase.verifyEqual(zRetained,zAutomatic,AbsTol=1e-12)
            testCase.verifyEqual(zAutomatic(1),-1,AbsTol=0)
            testCase.verifyEqual(zAutomatic(end),0,AbsTol=0)
        end

        function selectedColumnDefinesRootsWhenLabelsIncludeNullMode(testCase)
            solver = IMSolverSpectral(nEVP=96);
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1 0]);
            basisSet = solver.solveEVP(evp,nModes=5);
            z = basisSet.pointsFromModeRoots(nModes=4);
            F = basisSet.F(z(2:end-1));

            testCase.verifyEqual(basisSet.modeNumber,0:4,AbsTol=0)
            testCase.verifyEqual(basisSet.modeNumber(5),4,AbsTol=0)
            testCase.verifyEqual(length(z),6)
            testCase.verifyLessThan(max(abs(F(:,5))),1e-11)
        end

        function activeRobinGridKeepsPhysicalEndpoints(testCase)
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1, ...
                surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=5);
            z = basisSet.pointsFromModeRoots(nModes=4);
            values = basisSet.u(z);

            testCase.verifyEqual(z(1),-1,AbsTol=0)
            testCase.verifyEqual(z(end),0,AbsTol=0)
            testCase.verifyGreaterThan(abs(values(end,5)),1e-3)
            testCase.verifyLessThan(max(abs(values(2:end-1,5))),1e-10)
        end

        function allSpectralCoordinatesRecoverPhysicalRoots(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            coordinateKinds = ["z" "wkb" "density"];
            grids = zeros(6,length(coordinateKinds));
            for iKind = 1:length(coordinateKinds)
                solver = IMSolverSpectral(nEVP=128,coordinateKind=coordinateKinds(iKind));
                evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain);
                basisSet = solver.solveEVP(evp,nModes=5);
                grids(:,iKind) = basisSet.pointsFromModeRoots(nModes=4);
                G = basisSet.G(grids(2:end-1,iKind));
                testCase.verifyLessThan(max(abs(G(:,5))),1e-8)
            end

            testCase.verifyEqual(grids(:,2),grids(:,1),AbsTol=1e-3)
            testCase.verifyEqual(grids(:,3),grids(:,1),AbsTol=1e-3)
        end

        function highOrderExponentialRootsRemainUniqueAcrossCoordinates(testCase)
            [evp,zDomain] = testCase.exponentialProblem();
            coordinateKinds = ["z" "wkb" "density"];
            nModes = 120;
            nEVP = 384;
            grids = zeros(nModes+2,length(coordinateKinds));
            zReference = linspace(zDomain(1),zDomain(2),4001).';
            minimumSeparation = 1e-6*diff(zDomain);

            for iKind = 1:length(coordinateKinds)
                solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=coordinateKinds(iKind));
                basisSet = solver.solveEVP(evp,nModes=nModes+1);
                z = basisSet.pointsFromModeRoots(nModes=nModes);
                grids(:,iKind) = z;

                testCase.verifyEqual(length(z),nModes+2)
                testCase.verifyEqual(z(1),zDomain(1),AbsTol=0)
                testCase.verifyEqual(z(end),zDomain(2),AbsTol=0)
                testCase.verifyEqual(length(unique(z)),length(z))
                testCase.verifyGreaterThan(min(diff(z)),minimumSeparation)

                generatingMode = basisSet.nativeModes(:,nModes+1);
                rootValues = basisSet.solver.evaluatePhysicalDerivative(generatingMode,z(2:end-1),0);
                referenceValues = basisSet.solver.evaluatePhysicalDerivative(generatingMode,zReference,0);
                modeScale = max(abs(referenceValues));
                testCase.verifyGreaterThan(modeScale,0)
                testCase.verifyLessThan(max(abs(rootValues))/modeScale,1e-8)
            end

            testCase.verifyEqual(grids(:,2),grids(:,1),AbsTol=1e-3)
            testCase.verifyEqual(grids(:,3),grids(:,1),AbsTol=1e-3)
        end

        function highOrderAuxiliarySolveMatchesRetainedGeneratingMode(testCase)
            [evp,zDomain] = testCase.exponentialProblem();
            nModes = 120;
            solver = IMSolverSpectral(nEVP=384,coordinateKind="wkb");
            basisWithoutGeneratingMode = solver.solveEVP(evp,nModes=nModes);
            basisWithGeneratingMode = solver.solveEVP(evp,nModes=nModes+1);

            zAuxiliary = basisWithoutGeneratingMode.pointsFromModeRoots();
            zRetained = basisWithGeneratingMode.pointsFromModeRoots(nModes=nModes);

            testCase.verifyEqual(length(zAuxiliary),nModes+2)
            testCase.verifyEqual(length(unique(zAuxiliary)),length(zAuxiliary))
            testCase.verifyGreaterThan(min(diff(zAuxiliary)),1e-6*diff(zDomain))
            testCase.verifyEqual(zAuxiliary,zRetained,AbsTol=1e-6)
        end

        function highOrderWKBRootsConvergeWithResolution(testCase)
            [evp,zDomain] = testCase.exponentialProblem();
            nModes = 120;
            resolutions = [256 384];
            grids = zeros(nModes+2,length(resolutions));
            minimumSeparation = 1e-6*diff(zDomain);

            for iResolution = 1:length(resolutions)
                solver = IMSolverSpectral(nEVP=resolutions(iResolution),coordinateKind="wkb");
                basisSet = solver.solveEVP(evp,nModes=nModes+1);
                z = basisSet.pointsFromModeRoots(nModes=nModes);
                grids(:,iResolution) = z;

                testCase.verifyEqual(length(z),nModes+2)
                testCase.verifyEqual(length(unique(z)),length(z))
                testCase.verifyGreaterThan(min(diff(z)),minimumSeparation)
            end

            testCase.verifyEqual(grids(:,1),grids(:,2),AbsTol=1e-3)
        end

        function generatedGridBuildsFittedTransform(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisSet = solver.solveEVP(evp,nModes=4);
            z = basisSet.pointsFromModeRoots();
            [weights, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=4);
            transform = basisSet.discreteTransform(z=z,nModes=4);

            testCase.verifyGreaterThan(fit.exitFlag,0)
            testCase.verifyEqual(weights,fit.weights,AbsTol=0)
            testCase.verifyGreaterThanOrEqual(min(fit.weights),-1e-12)
            testCase.verifyEqual(sum(fit.weights),1,AbsTol=1e-10)
            testCase.verifyEqual(transform.z,z,AbsTol=0)
            testCase.verifyLessThan(transform.roundTripError,1e-11)
        end

        function unsupportedAndInvalidRequestsThrowStructuredErrors(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisSet = solver.solveEVP(evp,nModes=4);
            testCase.verifyError(@() basisSet.pointsFromModeRoots(nModes=5),"IMBasisSet:InvalidQuadratureModeCount")

            lowResolutionSolver = IMSolverSpectral(nEVP=4);
            lowResolutionBasis = lowResolutionSolver.solveEVP(evp,nModes=2);
            testCase.verifyError(@() lowResolutionBasis.pointsFromModeRoots(),"IMBasisSet:AuxiliaryModeUnavailable")

            inconsistentBasis = IMBasisSet(solver=basisSet.solver,evp=evp,nativeModes=basisSet.nativeModes, ...
                eigenvalues=basisSet.eigenvalues+1,modeNumber=basisSet.modeNumber);
            testCase.verifyError(@() inconsistentBasis.pointsFromModeRoots(),"IMBasisSet:AuxiliaryModeMismatch")

            finiteDifferenceSolver = IMSolverFiniteDifference(z=linspace(-1,0,65).');
            finiteDifferenceBasis = finiteDifferenceSolver.solveEVP(evp,nModes=3);
            testCase.verifyError(@() finiteDifferenceBasis.pointsFromModeRoots(nModes=2),"IMSolver:UnsupportedModeRoots")
        end

        function deficientRootSetIsRejected(testCase)
            [solver, evp] = testCase.regularProblem(8);
            solver = solver.configuredForEVP(evp);
            nativeModes = zeros(solver.nEVP,4);
            nativeModes(1,:) = 1;
            basisSet = IMBasisSet(solver=solver,evp=evp,nativeModes=nativeModes, ...
                eigenvalues=1:4,modeNumber=1:4);

            testCase.verifyError(@() basisSet.pointsFromModeRoots(nModes=3),"IMBasisSet:InsufficientQuadraturePoints")
        end

        function oldMethodNameIsAbsent(testCase)
            [solver, evp] = testCase.regularProblem(32);
            basisSet = solver.solveEVP(evp,nModes=2);

            testCase.verifyTrue(ismethod(basisSet,"pointsFromModeRoots"))
            testCase.verifyFalse(ismethod(basisSet,"quadraturePoints"))
        end
    end

    methods (Access = private)
        function [solver, evp] = regularProblem(~,nEVP)
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1, ...
                surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
        end

        function [evp,zDomain] = exponentialProblem(~)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain);
        end
    end
end
