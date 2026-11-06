classdef IMConstantStratificationValidationTests < matlab.unittest.TestCase

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
        function rigidHydrostaticGModesMatchConstantSpectrum(testCase)
            N0 = 5.2e-3;
            g = 9.81;
            zDomain = [-5000 0];
            nModes = 4;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());

            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain, g=g);
            basisSet = solution.internalModes(evp, nModes=nModes);
            D = diff(zDomain);
            expectedH = (N0*N0)./(g*((1:nModes)*pi/D).^2);

            testCase.verifyClass(solution, "IMConstantStratificationSolution")
            testCase.verifyClass(basisSet, "IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(basisSet.h, expectedH, RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifySize(basisSet.G(linspace(zDomain(1),zDomain(2),12).'), [12 nModes])
        end

        function freeSurfaceHydrostaticGModesIncludeBoundaryMode(testCase)
            N0 = 5.2e-3;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*ones(size(z));
            freeSurface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, surfaceBoundary=freeSurface);

            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=3);

            testCase.verifyEqual(basisSet.modeNumber(1), -1)
            testCase.verifyTrue(all(isfinite(basisSet.h)))
            testCase.verifySize(basisSet.F(linspace(zDomain(1),zDomain(2),8).'), [8 3])
        end

        function hydrostaticFModesIncludeZeroMode(testCase)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.neumann());

            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=4);
            z = linspace(zDomain(1), zDomain(2), 10).';

            testCase.verifyEqual(basisSet.modeNumber(1), 0)
            testCase.verifyTrue(isinf(basisSet.h(1)))
            F = basisSet.F(z);
            G = basisSet.G(z);
            factors = basisSet.normalizationFactors();
            testCase.verifyEqual(F(:,1), ones(size(z))/factors(1), AbsTol=1e-12)
            testCase.verifyEqual(G(:,1), zeros(size(z)), AbsTol=1e-12)
        end

        function analyticalGModesReturnExactSolvedDerivative(testCase)
            N0 = 5.2e-3;
            g = 9.81;
            zDomain = [-5000 0];
            nModes = 3;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain, g=g);
            basisSet = solution.internalModes(evp, nModes=nModes);
            z = linspace(zDomain(1), zDomain(2), 11).';
            s = z - zDomain(1);

            expected = zeros(length(z), nModes);
            for iMode = 1:nModes
                k_z = basisSet.metadata.verticalWavenumbers(iMode);
                signValue = (-1)^(basisSet.modeNumber(iMode)+1);
                expected(:,iMode) = signValue*k_z*cos(k_z*s);
            end
            expected = expected ./ basisSet.normalizationFactors();

            testCase.verifyEqual(basisSet.uz(z), expected, RelTol=1e-12, AbsTol=1e-14)
        end

        function analyticalFModesReturnExactSolvedDerivative(testCase)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            nModes = 4;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(rotationRate=0,N2=N2, zDomain=zDomain,surfaceBoundary=IMBoundaryCondition.neumann());
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain);
            basisSet = solution.internalModes(evp, nModes=nModes);
            z = linspace(zDomain(1), zDomain(2), 11).';
            s = z - zDomain(1);

            expected = zeros(length(z), nModes);
            for iMode = 1:nModes
                k_z = basisSet.metadata.verticalWavenumbers(iMode);
                signValue = (-1)^(basisSet.modeNumber(iMode)+1);
                expected(:,iMode) = -signValue*k_z*sin(k_z*s);
            end
            expected = expected ./ basisSet.normalizationFactors();

            actual = basisSet.uz(z);
            testCase.verifyEqual(actual, expected, RelTol=1e-12, AbsTol=1e-14)
            testCase.verifyEqual(actual(:,1), zeros(size(z)), AbsTol=1e-14)
        end

        function solutionCreatesConstantBasis(testCase)
            zDomain = [-3000 0];
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.waveModesAtWavenumber(rotationRate=0,N2=N2, zDomain=zDomain, k=1e-4,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain);

            basisSet = solution.internalModes(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(basisSet.evp.parameters.k, 1e-4, AbsTol=0)
        end

        function solutionCreatesMatchingDefaultEVP(testCase)
            N0 = 5.2e-3;
            zDomain = [-3000 0];
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain);

            basisSet = solution.internalModes(nModes=2);

            testCase.verifyClass(basisSet.evp, "IMInternalModes")
            testCase.verifyEqual(basisSet.evp.zDomain, zDomain)
            testCase.verifyEqual(basisSet.evp.N2(zDomain(:)), N0*N0*ones(2,1), RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeSelectionDiagnostics, struct())
        end

        function suppliedEVPRejectsDifferentStratification(testCase)
            zDomain = [-1000 0];
            for N0 = [0.01 1e-10]
                solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,rotationRate=0);
                N2 = @(z) 4*N0^2*ones(size(z));
                waveEVP = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
                apvEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,rotationRate=0,g0=Inf,gd=Inf,surfaceBoundary="freeSurface");

                testCase.verifyError(@() solution.internalModes(waveEVP,nModes=2),"IMConstantStratificationSolution:StratificationMismatch")
                testCase.verifyError(@() solution.internalModes(apvEVP,nModes=2),"IMConstantStratificationSolution:StratificationMismatch")
            end
        end

        function suppliedEVPRejectsInteriorProfileMismatch(testCase)
            N0 = 0.01;
            zDomain = [-1000 0];
            N2 = @(z) N0^2*(1+0.1*sin(pi*(z-zDomain(1))/diff(zDomain)).^2);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,rotationRate=0);

            testCase.verifyEqual(N2(zDomain),solution.N2(zDomain),RelTol=1e-12)
            testCase.verifyError(@() solution.internalModes(evp,nModes=2),"IMConstantStratificationSolution:StratificationMismatch")
        end

        function suppliedEVPAcceptsEquivalentIndependentProfile(testCase)
            N0 = 0.01;
            zDomain = [-1000 0];
            N2 = @(z) exp(2*log(N0))*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,rotationRate=0);

            basisSet = solution.internalModes(evp,nModes=2);
            z = linspace(zDomain(1),zDomain(2),31).';

            testCase.verifyEqual(basisSet.N2(z),basisSet.evp.N2(z),RelTol=1e-12)
            expectedH = N0^2./(evp.g*((1:2)*pi/diff(zDomain)).^2);
            testCase.verifyEqual(basisSet.h,expectedH,RelTol=1e-12)
        end

        function suppliedEVPRejectsMalformedStratification(testCase)
            zDomain = [-1000 0];
            solution = IMConstantStratificationSolution(N0=0.01,zDomain=zDomain,rotationRate=0);
            profiles = {@(z) NaN(size(z)),@(z) 1e-4*(1+1i)*ones(size(z)),@(z) 1e-4,@(z) zeros(numel(z),2),@(z) error("IMTest:ProfileFailure","The profile could not be evaluated.")};
            for iProfile = 1:numel(profiles)
                evp = IMInternalModes.hydrostaticGModes(N2=profiles{iProfile},zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());

                testCase.verifyError(@() solution.internalModes(evp,nModes=2),"IMConstantStratificationSolution:InvalidStratification")
            end
        end

        function integerProfileUsesStructuredMismatchError(testCase)
            zDomain = [-1000 0];
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) zeros(size(z),"uint8"),zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(N0=0.01,zDomain=zDomain,rotationRate=0);

            testCase.verifyError(@() solution.internalModes(evp,nModes=2),"IMConstantStratificationSolution:StratificationMismatch")
        end

        function suppliedEVPDomainMismatchTakesPrecedence(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) 4e-4*ones(size(z)),zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(N0=0.01,zDomain=[-2000 0],rotationRate=0);

            testCase.verifyError(@() solution.internalModes(evp,nModes=2),"IMConstantStratificationSolution:DomainMismatch")
        end

        function suppliedEVPOwnsGravityAndRotation(testCase)
            N0 = 0.01;
            zDomain = [-1000 0];
            k = 1e-3;
            evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) N0^2*ones(size(z)),zDomain=zDomain,k=k,g=10,rotationRate=1e-4,latitude=30,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,g=9.81,rotationRate=0);

            basisSet = solution.internalModes(evp,nModes=2);
            expectedH = (N0^2-evp.f0^2)./(evp.g*(((1:2)*pi/diff(zDomain)).^2+k^2));

            testCase.verifyEqual(basisSet.evp,evp)
            testCase.verifyEqual(basisSet.h,expectedH,RelTol=1e-12)
        end

        function coordinateKindSolversAgreeWithRigidAnalyticalDepths(testCase)
            N0 = 5.2e-3;
            g = 9.81;
            zDomain = [-1000 0];
            nModes = 2;
            nEVP = 48;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(rotationRate=0,N2=N2, zDomain=zDomain, g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = IMConstantStratificationSolution(rotationRate=0,N0=N0, zDomain=zDomain, g=g);
            analytical = solution.internalModes(evp, nModes=nModes);

            coordinateKinds = ["z", "wkb", "density"];
            for coordinateKind = coordinateKinds
                solver = IMSolverSpectral(nEVP=nEVP, coordinateKind=coordinateKind);
                numerical = solver.solveEVP(evp, nModes=nModes);
                testCase.verifyEqual(numerical.h, analytical.h, RelTol=1e-4)
            end
        end

        function constantZeroAPVModesUseCanonicalVariables(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            z = linspace(zDomain(1), zDomain(2), 9).';
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0));

            exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary="rigidLid");

            testCase.verifyClass(exactModes,"IMAnalyticalGeostrophicZeroAPVModesBasis")
            testCase.verifySize(exactModes.F(z),[9 2 2])
            testCase.verifySize(exactModes.G(z),[9 2 2])
            testCase.verifyEqual(exactModes.endpointResponseMetric,repmat(eye(2),1,1,2),AbsTol=2e-13)
        end
    end
end
