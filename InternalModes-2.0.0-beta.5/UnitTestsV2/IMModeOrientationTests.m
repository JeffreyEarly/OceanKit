classdef IMModeOrientationTests < matlab.unittest.TestCase

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
        function numericalOrientationIsInvariantToRawEigenvectorSigns(testCase)
            z = linspace(-1,0,65).';
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1 0]);
            solver = IMSolverFiniteDifference(z=z).configuredForEVP(evp);
            GValues = [-ones(size(solver.zNative)), solver.zNative];

            first = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=GValues, ...
                eigenvalues=[1 2],modeNumber=[1 2]);
            second = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=-GValues, ...
                eigenvalues=[1 2],modeNumber=[1 2]);
            first = first.orientModeSigns();
            second = second.orientModeSigns();

            zCheck = linspace(-1,0,101).';
            testCase.verifyEqual(first.rawVariable("G",zCheck),second.rawVariable("G",zCheck),AbsTol=1e-12)
            testCase.verifyGreaterThan(first.rawVariable("G",-1e-6),zeros(1,2))
            testCase.verifyEqual(first.metadata.modeOrientation,"shallowInteriorGPositive-v1")
        end

        function FFormNullModeUsesDeterministicFallback(testCase)
            z = linspace(-1,0,65).';
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1 0]);
            solver = IMSolverFiniteDifference(z=z).configuredForEVP(evp);
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=-ones(size(z)), ...
                eigenvalues=0,modeNumber=0);

            basisSet = basisSet.orientModeSigns();

            testCase.verifyGreaterThan(basisSet.rawVariable("F",0),0)
            testCase.verifyEqual(basisSet.rawVariable("G",z),zeros(size(z)),AbsTol=1e-12)
        end

        function APVAndMDAFamiliesUseTheSameShallowInteriorConvention(testCase)
            D = 4000;
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            solver = IMSolverSpectral(nEVP=128);
            apvEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g0=0.02,gd=0.03);
            mdaEVP = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=[-D 0],g0=0.02,gd=0.03);

            apvBasis = solver.solveEVP(apvEVP,nModes=6);
            mdaBasis = solver.solveEVP(mdaEVP,nModes=6);
            zShallow = -1e-6*D;

            testCase.verifyGreaterThan(apvBasis.G(zShallow,normalization="wMax"),zeros(1,6))
            testCase.verifyGreaterThan(mdaBasis.G(zShallow,normalization="wMax"),zeros(1,6))
            testCase.verifyEqual(apvBasis.metadata.modeOrientation,IMModeOrientationTools.convention)
            testCase.verifyEqual(mdaBasis.metadata.modeOrientation,IMModeOrientationTools.convention)
        end

        function analyticalAndNumericalRigidSurfaceModesHaveMatchingPhase(testCase)
            D = 3000;
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-D 0]);
            analytical = IMConstantStratificationSolution(N0=N0,zDomain=[-D 0]).internalModes(evp,nModes=4);
            numerical = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=4);
            z = linspace(-D,0,2001).';
            analyticalG = analytical.G(z,normalization="wMax");
            numericalG = numerical.G(z,normalization="wMax");

            correlation = sum(analyticalG.*numericalG,1);
            testCase.verifyGreaterThan(correlation,zeros(1,4))
            testCase.verifyGreaterThan(analytical.G(-1e-6*D,normalization="wMax"),zeros(1,4))
            testCase.verifyEqual(analytical.metadata.modeOrientation,IMModeOrientationTools.convention)
        end
    end
end
