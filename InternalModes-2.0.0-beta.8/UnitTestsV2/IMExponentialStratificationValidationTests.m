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
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());

            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);
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
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, surfaceBoundary=freeSurface);

            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);

            testCase.verifyEqual(basisSet.modeNumber(1), -1)
            testCase.verifyTrue(all(isfinite(basisSet.metadata.roots)))
        end

        function hydrostaticFModesIncludeNullBranch(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.neumann());

            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            testCase.verifyEqual(basisSet.modeNumber(1), 0)
            F = basisSet.F(z);
            G = basisSet.G(z);
            factors = basisSet.normalizationFactors();
            testCase.verifyEqual(F(:,1), ones(size(z))/factors(1), AbsTol=1e-12)
            testCase.verifyEqual(G(:,1), zeros(size(z)), AbsTol=1e-12)
        end

        function analyticalGModesReturnSolvedDerivative(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            expected = basisSet.F(z) ./ basisSet.h;

            testCase.verifyEqual(basisSet.uz(z), expected, RelTol=1e-12, AbsTol=1e-14)
        end

        function analyticalFModesReturnSolvedDerivative(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.neumann());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            expected = -(N2(z)/evp.g).*basisSet.G(z);

            actual = basisSet.uz(z);
            testCase.verifyEqual(actual, expected, RelTol=1e-12, AbsTol=1e-14)
            testCase.verifyEqual(actual(:,1), zeros(size(z)), AbsTol=1e-14)
            testCase.verifyTrue(all(isfinite(actual(:))))
        end

        function solutionCreatesExponentialBasis(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.waveModesAtFrequency(rotationRate=0,N2=N2, zDomain=zDomain, omega=1e-3,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);

            basisSet = solution.internalModes(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(basisSet.evp.parameters.omega, 1e-3, AbsTol=0)
        end

        function providedEVPAmplitudeMustMatchAnalyticalProfile(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) (2*N0)^2*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            testCase.verifyError(@() solution.internalModes(evp,nModes=2), "IMExponentialStratificationSolution:StratificationMismatch")
        end

        function providedEVPDecayScaleMustMatchAnalyticalProfile(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0^2*exp(2*z/(2*b));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            testCase.verifyError(@() solution.internalModes(evp,nModes=2), "IMExponentialStratificationSolution:StratificationMismatch")
        end

        function providedEVPInteriorMustMatchAnalyticalProfile(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            D = diff(zDomain);
            N2 = @(z) N0^2*exp(2*z/b).*(1 + 0.1*(z-zDomain(1)).*(zDomain(2)-z)/D^2);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            testCase.verifyEqual(evp.N2(zDomain),solution.N2(zDomain),AbsTol=0)
            testCase.verifyError(@() solution.internalModes(evp,nModes=2), "IMExponentialStratificationSolution:StratificationMismatch")
        end

        function deepWeakStratificationUsesItsLocalScale(testCase)
            N0 = 5.2e-3;
            b = 100;
            zDomain = [-5000 0];
            N2 = @(z) N0^2*exp(2*z/b).*(1+(z < -4000));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            testCase.verifyEqual(evp.N2(0),solution.N2(0),AbsTol=0)
            testCase.verifyError(@() solution.internalModes(evp,nModes=2),"IMExponentialStratificationSolution:StratificationMismatch")
        end

        function equivalentIndependentProfileIsAccepted(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) (N0*exp(z/b)).^2;
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2,zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            basisSet = solution.internalModes(evp,nModes=2);
            z = linspace(zDomain(1),zDomain(2),37).';

            testCase.verifyClass(basisSet,"IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(basisSet.N2(z),basisSet.evp.N2(z),RelTol=1e-12)
        end

        function malformedProvidedEVPProfilesAreRejected(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);
            profiles = {@(z) NaN(size(z)), @(z) N0^2*exp(2*z/b) + 1i*ones(size(z)), @(z) N0^2};

            for iProfile = 1:numel(profiles)
                evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=profiles{iProfile},zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());
                testCase.verifyError(@() solution.internalModes(evp,nModes=2), "IMExponentialStratificationSolution:InvalidStratification")
            end
        end

        function apvProfileMismatchIsRejectedBeforeCatalogSolve(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) (2*N0)^2*exp(2*z/b);
            evp = IMInternalModes.geostrophicAPVModes(rotationRate=0,N2=N2,zDomain=zDomain,g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            testCase.verifyError(@() solution.internalModes(evp,nModes=2), "IMExponentialStratificationSolution:StratificationMismatch")
        end

        function solutionCreatesMatchingDefaultEVP(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0,b=b,zDomain=zDomain);

            basisSet = solution.internalModes(nModes=2);
            z = linspace(zDomain(1),zDomain(2),37).';

            testCase.verifyClass(basisSet.evp,"IMInternalModes")
            testCase.verifyEqual(basisSet.evp.zDomain,zDomain)
            testCase.verifyEqual(basisSet.N2(z),basisSet.evp.N2(z),AbsTol=0)
        end

        function providedEVPDomainMustMatchAnalyticalDomain(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) (2*N0)^2*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.dirichlet());

            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=[-4000 0]);

            testCase.verifyError(@() solution.internalModes(evp, nModes=2), "IMExponentialStratificationSolution:DomainMismatch")
        end

        function unsupportedBoundaryIsRejected(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, bottomBoundary=IMBoundaryCondition.neumann(),surfaceBoundary=IMBoundaryCondition.dirichlet());

            solution = IMExponentialStratificationSolution(rotationRate=0,N0=N0, b=b, zDomain=zDomain);

            testCase.verifyError(@() solution.internalModes(evp, nModes=2), "IMExponentialStratificationSolution:UnsupportedBoundary")
        end

        function exponentialZeroAPVModesUseCanonicalVariables(testCase)
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            z = linspace(zDomain(1), zDomain(2), 9).';
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain, rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0));

            exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary="rigidLid");

            testCase.verifyClass(exactModes,"IMAnalyticalGeostrophicZeroAPVModesBasis")
            testCase.verifySize(exactModes.F(z),[9 2 2])
            testCase.verifySize(exactModes.G(z),[9 2 2])
            testCase.verifyEqual(exactModes.endpointResponseMetric,repmat(eye(2),1,1,2),AbsTol=2e-13)
        end
    end
end
