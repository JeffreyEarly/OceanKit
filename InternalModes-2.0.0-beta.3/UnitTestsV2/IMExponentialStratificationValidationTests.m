classdef IMExponentialStratificationValidationTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
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
        function rigidHydrostaticGModesEvaluateExponentialBasis(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 24).';

            testCase.verifyClass(solution, "IMExponentialStratificationSolution")
            testCase.verifyClass(basisSet, "IMAnalyticalInternalModesBasis")
            testCase.verifyTrue(all(isfinite(basisSet.h)))
            testCase.verifySize(basisSet.G(z), [24 3])
            testCase.verifySize(basisSet.F(z), [24 3])
            testCase.verifyEqual(basisSet.modeNumber, 1:3)
        end

        function freeSurfaceHydrostaticGModesIncludeBoundaryMode(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            freeSurface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, surfaceBoundary=freeSurface);

            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);

            testCase.verifyEqual(basisSet.modeNumber(1), -1)
            testCase.verifyTrue(all(isfinite(basisSet.metadata.roots)))
        end

        function hydrostaticFModesIncludeNullBranch(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);

            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            testCase.verifyEqual(basisSet.modeNumber(1), 0)
            F = basisSet.F(z, normalization=Normalization.unity);
            G = basisSet.G(z, normalization=Normalization.unity);
            factors = basisSet.normalizationFactors(Normalization.unity);
            testCase.verifyEqual(F(:,1), ones(size(z))/factors(1), AbsTol=1e-12)
            testCase.verifyEqual(G(:,1), zeros(size(z)), AbsTol=1e-12)
        end

        function analyticalGModesReturnSolvedDerivative(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            expected = basisSet.F(z, normalization="unity") ./ basisSet.h;

            testCase.verifyEqual(basisSet.uz(z, normalization="unity"), expected, RelTol=1e-12, AbsTol=1e-14)
        end

        function analyticalFModesReturnSolvedDerivative(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            expected = -(N2(z)/evp.g).*basisSet.G(z, normalization="unity");

            actual = basisSet.uz(z, normalization="unity");
            testCase.verifyEqual(actual, expected, RelTol=1e-12, AbsTol=1e-14)
            testCase.verifyEqual(actual(:,1), zeros(size(z)), AbsTol=1e-14)
            testCase.verifyTrue(all(isfinite(actual(:))))
        end

        function solutionCreatesExponentialBasis(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.waveModesAtFrequency(N2=N2, zDomain=zDomain, omega=1e-3);
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);

            basisSet = solution.internalModes(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(basisSet.evp.parameters.omega, 1e-3, AbsTol=0)
        end

        function providedEVPDomainMustMatchAnalyticalDomain(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=[-4000 0]);

            testCase.verifyError(@() solution.internalModes(evp, nModes=2), "IMExponentialStratificationSolution:DomainMismatch")
        end

        function unsupportedBoundaryIsRejected(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, bottomBoundary=IMBoundaryCondition.neumann());

            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain);

            testCase.verifyError(@() solution.internalModes(evp, nModes=2), "IMExponentialStratificationSolution:UnsupportedBoundary")
        end

        function exponentialZeroAPVModesUseCanonicalVariables(testCase)
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            z = linspace(zDomain(1), zDomain(2), 9).';
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain, f0=f0);

            exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary="rigidLid");

            testCase.verifyClass(exactModes,"IMAnalyticalGeostrophicZeroAPVModesBasis")
            testCase.verifySize(exactModes.F(z),[9 2 2])
            testCase.verifySize(exactModes.G(z),[9 2 2])
            testCase.verifyEqual(exactModes.endpointResponseMetric,repmat(eye(2),1,1,2),AbsTol=2e-13)
        end
    end
end
