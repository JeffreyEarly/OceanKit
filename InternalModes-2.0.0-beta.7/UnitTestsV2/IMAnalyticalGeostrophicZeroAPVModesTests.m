classdef IMAnalyticalGeostrophicZeroAPVModesTests < matlab.unittest.TestCase

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
        function analyticalHierarchyAdvertisesConcreteOperations(testCase)
            constant = IMConstantStratificationSolution(f0=1);
            exponential = IMExponentialStratificationSolution(f0=1);

            retiredMethods = ["internalMode"+"Availability", "sqg"+"Availability", "sqgModes"+"AtWavenumber"];
            for methodName = retiredMethods
                testCase.verifyFalse(ismethod(constant,methodName))
                testCase.verifyFalse(ismethod(exponential,methodName))
            end
            testCase.verifyTrue(ismethod(constant,"internalModes"))
            testCase.verifyTrue(ismethod(constant,"geostrophicZeroAPVModesAtWavenumber"))
            testCase.verifyTrue(ismethod(exponential,"internalModes"))
            testCase.verifyTrue(ismethod(exponential,"geostrophicZeroAPVModesAtWavenumber"))
            summaryText = evalc("constant.summarize()");
            testCase.verifySubstring(summaryText,"internalModes")
            testCase.verifySubstring(summaryText,"geostrophicZeroAPVModesAtWavenumber")
        end

        function constantResponsesCoverEndpointsAndSurfaceConventions(testCase)
            [N0,f0,g,zDomain] = testCase.constantParameters();
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
            endpointCases = {"surface","bottom",["surface" "bottom"]};
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                for iEndpoints = 1:numel(endpointCases)
                    exactModes = solution.geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4],endpoints=endpointCases{iEndpoints},surfaceBoundary=surfaceBoundary);
                    responses = testCase.endpointResponses(exactModes);
                    rows = testCase.endpointRows(exactModes.endpoints);
                    inactiveRows = setdiff(1:2,rows);
                    for iK = 1:numel(exactModes.k)
                        testCase.verifyEqual(responses(rows,:,iK),eye(numel(rows)),AbsTol=3e-13)
                        testCase.verifyEqual(responses(inactiveRows,:,iK),zeros(numel(inactiveRows),numel(rows)),AbsTol=3e-13)
                        testCase.verifyEqual(exactModes.endpointResponseMetric(:,:,iK),eye(numel(rows)),AbsTol=5e-13)
                    end
                end
            end
        end

        function exponentialResponsesCoverEndpointsAndSurfaceConventions(testCase)
            [N0,b,f0,g,zDomain] = testCase.exponentialParameters();
            solution = IMExponentialStratificationSolution(N0=N0,b=b,zDomain=zDomain,f0=f0,g=g);
            endpointCases = {"surface","bottom",["surface" "bottom"]};
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                for iEndpoints = 1:numel(endpointCases)
                    exactModes = solution.geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4],endpoints=endpointCases{iEndpoints},surfaceBoundary=surfaceBoundary);
                    responses = testCase.endpointResponses(exactModes);
                    rows = testCase.endpointRows(exactModes.endpoints);
                    inactiveRows = setdiff(1:2,rows);
                    for iK = 1:numel(exactModes.k)
                        testCase.verifyEqual(responses(rows,:,iK),eye(numel(rows)),AbsTol=5e-13)
                        testCase.verifyEqual(responses(inactiveRows,:,iK),zeros(numel(inactiveRows),numel(rows)),AbsTol=5e-13)
                    end
                end
            end
        end

        function constantModesMatchIndependentHyperbolicSolve(testCase)
            [N0,f0,g,zDomain] = testCase.constantParameters();
            k = [1e-4 2.3e-4];
            z = linspace(zDomain(1),zDomain(2),81).';
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary=surfaceBoundary);
                actualF = exactModes.F(z);
                actualG = exactModes.G(z);
                for iK = 1:numel(k)
                    [expectedF,expectedG] = testCase.independentConstantModes(z,k(iK),N0,f0,g,zDomain,surfaceBoundary);
                    testCase.verifyEqual(actualF(:,:,iK),expectedF,RelTol=2e-13,AbsTol=2e-13)
                    testCase.verifyEqual(actualG(:,:,iK),expectedG,RelTol=2e-13,AbsTol=2e-13)
                end
            end
        end

        function exponentialModesMatchIndependentBesselSolve(testCase)
            [N0,b,f0,g,zDomain] = testCase.exponentialParameters();
            k = [1e-4 1.8e-4];
            z = linspace(zDomain(1),zDomain(2),81).';
            solution = IMExponentialStratificationSolution(N0=N0,b=b,zDomain=zDomain,f0=f0,g=g);
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary=surfaceBoundary);
                actualF = exactModes.F(z);
                actualG = exactModes.G(z);
                for iK = 1:numel(k)
                    [expectedF,expectedG] = testCase.independentExponentialModes(z,k(iK),N0,b,f0,g,zDomain,surfaceBoundary);
                    testCase.verifyEqual(actualF(:,:,iK),expectedF,RelTol=5e-13,AbsTol=5e-13)
                    testCase.verifyEqual(actualG(:,:,iK),expectedG,RelTol=5e-13,AbsTol=5e-13)
                end
            end
        end

        function rigidLidLegacyDerivativeNormalizationIsRescaled(testCase)
            [N0,f0,g,zDomain] = testCase.constantParameters();
            k = 1.4e-4;
            z = linspace(zDomain(1),zDomain(2),41).';
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
            exactSurface = solution.geostrophicZeroAPVModesAtWavenumber(k,endpoints="surface",surfaceBoundary="rigidLid");

            depth = diff(zDomain);
            zRelative = z-zDomain(2);
            m = k*N0/f0;
            legacyProfile = (1/N0)*(exp(m*zRelative)+exp(-m*(zRelative+2*depth)))/(k*(1-exp(-2*m*depth)));
            legacyResponse = -g/(f0*N0*N0);
            testCase.verifyEqual(exactSurface.F(z),legacyProfile/legacyResponse,RelTol=3e-13,AbsTol=3e-13)
        end

        function exactQuadraticFormsMatchContinuousEnergy(testCase)
            [N0,f0,g,zDomain] = testCase.constantParameters();
            k = 1.7e-4;
            z = linspace(zDomain(1),zDomain(2),20001).';
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary=surfaceBoundary);
                F = exactModes.F(z);
                G = exactModes.G(z);
                F = F(:,:,1);
                G = G(:,:,1);
                integrand = zeros(numel(z),2,2);
                for iMode = 1:2
                    for jMode = 1:2
                        integrand(:,iMode,jMode) = k^2*F(:,iMode).*F(:,jMode) + (f0^2*N0^2/g^2)*G(:,iMode).*G(:,jMode);
                    end
                end
                expectedEnergy = squeeze(trapz(z,integrand,1))/(2*k^4);
                if surfaceBoundary == "freeSurface"
                    expectedEnergy = expectedEnergy + (f0^2/(2*g*k^4))*(F(end,:).'*F(end,:));
                end
                testCase.verifyEqual(exactModes.energyMatrix(:,:,1),expectedEnergy,RelTol=8e-8,AbsTol=2e-8)
                expectedGeneralized = expectedEnergy-0.03*exactModes.surfaceBuoyancyMatrix(:,:,1)+0.02*exactModes.bottomBuoyancyMatrix(:,:,1);
                testCase.verifyEqual(exactModes.generalizedEnergyMatrix(g0=-0.03,gd=0.02),expectedGeneralized,RelTol=8e-8,AbsTol=2e-8)
            end
        end

        function analyticalRotationsMatchNumericalContract(testCase)
            [N0,f0,g,zDomain] = testCase.constantParameters();
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*ones(size(z));
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
            exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,f0=f0,g=g,k=k);
            numericalModes = IMSolverSpectral(nEVP=80).solveGeostrophicZeroAPVModes(problem);
            g0 = -0.03;
            gd = 0.01;

            exactDepth = exactModes.rotateBoundaryDepth(g0=g0,gd=gd);
            numericalDepth = numericalModes.rotateBoundaryDepth(g0=g0,gd=gd);
            exactSurface = exactModes.rotateSurfaceBuoyancy(g0=g0,gd=gd);
            numericalSurface = numericalModes.rotateSurfaceBuoyancy(g0=g0,gd=gd);
            testCase.verifyEqual(exactDepth.h0,numericalDepth.h0,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEqual(exactDepth.rotationEigenvalues,numericalDepth.rotationEigenvalues,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEqual(exactSurface.rotationEigenvalues,numericalSurface.rotationEigenvalues,RelTol=2e-10,AbsTol=2e-10)
            testCase.verifyEmpty(exactSurface.h0)
            testCase.verifyEqual(exactDepth.signatures,sign(exactDepth.rotationEigenvalues))
            testCase.verifyLessThan(max(exactDepth.rotationResiduals.relativePencilResidual,[],"all"),1e-12)

            Hg = exactModes.generalizedEnergyMatrix(g0=g0,gd=gd);
            custom = exactModes.rotateWithPencil(name="custom",leftMatrix=Hg,rightMatrix=exactModes.endpointResponseMetric);
            testCase.verifyEmpty(custom.h0)
            testCase.verifyEqual(custom.rotationEigenvalues,exactDepth.rotationEigenvalues,RelTol=1e-12,AbsTol=1e-12)
        end

        function exponentialExactModesConvergeFromNumericalSolver(testCase)
            [N0,b,f0,g,zDomain] = testCase.exponentialParameters();
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*exp(2*z/b);
            solution = IMExponentialStratificationSolution(N0=N0,b=b,zDomain=zDomain,f0=f0,g=g);
            exactModes = solution.geostrophicZeroAPVModesAtWavenumber(k,surfaceBoundary="freeSurface");
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,f0=f0,g=g,k=k,surfaceBoundary="freeSurface");
            numericalModes = IMSolverSpectral(nEVP=128).solveGeostrophicZeroAPVModes(problem);
            z = linspace(zDomain(1),zDomain(2),101).';

            testCase.verifyEqual(numericalModes.F(z),exactModes.F(z),RelTol=2e-8,AbsTol=2e-9)
            testCase.verifyEqual(numericalModes.G(z),exactModes.G(z),RelTol=2e-8,AbsTol=2e-9)
            testCase.verifyEqual(numericalModes.energyMatrix,exactModes.energyMatrix,RelTol=3e-8,AbsTol=3e-8)
        end

        function signedCoriolisLeavesStructuresInvariant(testCase)
            [N0,~,g,zDomain] = testCase.constantParameters();
            positive = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=1e-4,g=g);
            negative = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=-1e-4,g=g);
            z = linspace(zDomain(1),zDomain(2),31).';
            positiveModes = positive.geostrophicZeroAPVModesAtWavenumber(1e-4);
            negativeModes = negative.geostrophicZeroAPVModesAtWavenumber(1e-4);
            testCase.verifyEqual(negativeModes.F(z),positiveModes.F(z),AbsTol=0)
            testCase.verifyEqual(negativeModes.G(z),positiveModes.G(z),AbsTol=0)
            testCase.verifyEqual(negativeModes.energyMatrix,positiveModes.energyMatrix,AbsTol=0)
        end

        function invalidRequestsThrowConcreteErrors(testCase)
            constant = IMConstantStratificationSolution(f0=0);
            testCase.verifyError(@() constant.geostrophicZeroAPVModesAtWavenumber(1),"IMGeostrophicZeroAPVModes:InvalidCoriolis")
            constant = IMConstantStratificationSolution(f0=1);
            testCase.verifyError(@() constant.geostrophicZeroAPVModesAtWavenumber(1,endpoints=string.empty(1,0)),"IMGeostrophicZeroAPVModes:NoEndpoint")
            testCase.verifyError(@() constant.geostrophicZeroAPVModesAtWavenumber(1,endpoints="middle"),"IMGeostrophicZeroAPVModes:InvalidEndpoint")
            testCase.verifyError(@() constant.geostrophicZeroAPVModesAtWavenumber(1,endpoints=["surface" "surface"]),"IMGeostrophicZeroAPVModes:DuplicateEndpoint")
            testCase.verifyError(@() constant.geostrophicZeroAPVModesAtWavenumber(1,surfaceBoundary="other"),"MATLAB:validators:mustBeMember")

            exponential = IMExponentialStratificationSolution(f0=0);
            testCase.verifyError(@() exponential.geostrophicZeroAPVModesAtWavenumber(1),"IMGeostrophicZeroAPVModes:InvalidCoriolis")
        end
    end

    methods (Static, Access = private)
        function [N0,f0,g,zDomain] = constantParameters()
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-1000 0];
        end

        function [N0,b,f0,g,zDomain] = exponentialParameters()
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-4000 0];
        end

        function responses = endpointResponses(basisSet)
            zEndpoints = [basisSet.zDomain(2); basisSet.zDomain(1)];
            F = basisSet.F(zEndpoints);
            G = basisSet.G(zEndpoints);
            nEndpoints = numel(basisSet.endpoints);
            nK = numel(basisSet.k);
            responses = zeros(2,nEndpoints,nK);
            for iK = 1:nK
                surface = G(1,:,iK);
                if basisSet.surfaceBoundary == "freeSurface"
                    surface = surface-F(1,:,iK);
                end
                responses(:,:,iK) = [surface; G(2,:,iK)];
            end
        end

        function rows = endpointRows(endpoints)
            rows = zeros(1,numel(endpoints));
            rows(endpoints == "surface") = 1;
            rows(endpoints == "bottom") = 2;
        end

        function [F,G] = independentConstantModes(z,k,N0,f0,g,zDomain,surfaceBoundary)
            m = k*N0/abs(f0);
            midpoint = mean(zDomain);
            y = z(:)-midpoint;
            rawF = [cosh(m*y) sinh(m*y)/m];
            rawFz = [m*sinh(m*y) cosh(m*y)];
            surfaceF = rawF(end,:);
            surfaceG = -(g/(N0*N0))*rawFz(end,:);
            bottomG = -(g/(N0*N0))*rawFz(1,:);
            if surfaceBoundary == "freeSurface"
                surfaceG = surfaceG-surfaceF;
            end
            coefficients = [surfaceG; bottomG]\eye(2);
            F = rawF*coefficients;
            G = -(g/(N0*N0))*rawFz*coefficients;
        end

        function [F,G] = independentExponentialModes(z,k,N0,b,f0,g,zDomain,surfaceBoundary)
            eta = b*k*N0/abs(f0);
            evaluate = @(zIn) IMAnalyticalGeostrophicZeroAPVModesTests.independentExponentialFundamental(zIn,eta,b);
            [rawF,rawFz] = evaluate(z(:));
            [surfaceF,surfaceFz] = evaluate(zDomain(2));
            [~,bottomFz] = evaluate(zDomain(1));
            surfaceG = -(g/(N0*N0))*surfaceFz;
            bottomN2 = N0*N0*exp(2*zDomain(1)/b);
            bottomG = -(g/bottomN2)*bottomFz;
            if surfaceBoundary == "freeSurface"
                surfaceG = surfaceG-surfaceF;
            end
            coefficients = [surfaceG; bottomG]\eye(2);
            F = rawF*coefficients;
            N2Values = N0*N0*exp(2*z(:)/b);
            G = -(g./N2Values).*rawFz*coefficients;
        end

        function [F,Fz] = independentExponentialFundamental(z,eta,b)
            z = z(:);
            exponential = exp(z/b);
            x = eta*exponential;
            F = [exponential.*besseli(1,x) exponential.*besselk(1,x)];
            derivativeScale = (x/b).*exponential;
            Fz = [derivativeScale.*besseli(0,x) -derivativeScale.*besselk(0,x)];
        end
    end
end
