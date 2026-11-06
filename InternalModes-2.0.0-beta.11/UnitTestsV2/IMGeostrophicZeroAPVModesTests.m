classdef IMGeostrophicZeroAPVModesTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.originalPath = path;
            addpath(repoRoot);
            addpath(fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function constructorOwnsEndpointsRatherThanEnergyWeights(testCase)
            N2 = @(z) ones(size(z));
            surface = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,endpoints="surface");
            bottom = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,endpoints="bottom",surfaceBoundary="rigidLid");
            both = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=[2 3]);

            testCase.verifyEqual(surface.endpoints,"surface")
            testCase.verifyEqual(bottom.endpoints,"bottom")
            testCase.verifyEqual(both.endpoints,["surface" "bottom"])
            testCase.verifyEqual(bottom.surfaceBoundary,"rigidLid")
            testCase.verifyEqual(both.modesPerWavenumber(),2)
            testCase.verifyFalse(isprop(surface,"g0"))
            testCase.verifyFalse(isprop(surface,"gd"))
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=0,latitude=33,k=2),"IMGeostrophicZeroAPVModes:InvalidCoriolis")
            empty = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,endpoints=strings(1,0));
            testCase.verifyEqual(empty.modesPerWavenumber(),0)
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,endpoints="middle"),"IMGeostrophicZeroAPVModes:InvalidEndpoint")
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,endpoints=["surface" "surface"]),"IMGeostrophicZeroAPVModes:DuplicateEndpoint")
            testCase.verifyTrue(testCase.throwsAny(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,k=2,surfaceBoundary="bad")))
        end

        function oldPublicNamesAreRemoved(testCase)
            testCase.verifyEqual(exist("IMSurfaceGeostrophicModes","class"),0)
            testCase.verifyEqual(exist("IMSurfaceGeostrophicModesBasis","class"),0)
            methodNames = string(methods("IMSolver"));
            testCase.verifyFalse(ismember("solveSurfaceGeostrophicModes",methodNames))
            testCase.verifyFalse(ismember("configuredForSurfaceGeostrophicModes",methodNames))
        end

        function canonicalResponsesCoverEveryEndpointAndSurfaceConvention(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-1200 0];
            N2 = @(z) N0^2*ones(size(z));
            endpointCases = {"surface","bottom",["surface" "bottom"]};
            surfaceCases = ["freeSurface" "rigidLid"];
            for iSurface = 1:numel(surfaceCases)
                for iEndpoints = 1:numel(endpointCases)
                    problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,k=[1e-4 2e-4],endpoints=endpointCases{iEndpoints},surfaceBoundary=surfaceCases(iSurface));
                    basisSet = IMSolverSpectral(nEVP=80).solveGeostrophicZeroAPVModes(problem);
                    responses = testCase.endpointResponses(basisSet);
                    selectedRows = testCase.endpointRows(problem.endpoints);
                    for iK = 1:numel(problem.k)
                        testCase.verifyEqual(responses(selectedRows,:,iK),eye(numel(problem.endpoints)),AbsTol=2e-10)
                        inactiveRows = setdiff(1:2,selectedRows);
                        testCase.verifyEqual(responses(inactiveRows,:,iK),zeros(numel(inactiveRows),numel(problem.endpoints)),AbsTol=2e-10)
                        testCase.verifyEqual(basisSet.endpointResponseMetric(:,:,iK),eye(numel(problem.endpoints)),AbsTol=5e-10)
                    end
                    testCase.verifySize(basisSet.F(linspace(zDomain(1),zDomain(2),7).'),[7 numel(problem.endpoints) numel(problem.k)])
                    testCase.verifySize(basisSet.G(linspace(zDomain(1),zDomain(2),7).'),[7 numel(problem.endpoints) numel(problem.k)])
                end
            end
        end

        function solverFactorsOncePerWavenumberAndUsesMultipleRightHandSides(testCase)
            N2 = @(z) ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,k=[1 2 1 3 2]);
            counter = IMGeostrophicZeroAPVTestCounter();
            solver = IMInstrumentedGeostrophicZeroAPVSolver(nEVP=32,solveCounter=counter);
            basisSet = solver.solveGeostrophicZeroAPVModes(problem);
            testCase.verifyEqual(counter.callCount,3)
            testCase.verifyEqual(counter.rightHandSideCounts,[2 2 2])
            testCase.verifyEqual(basisSet.metadata.factorizations,3)
            F = basisSet.F(linspace(-1,0,7).');
            testCase.verifyEqual(F(:,:,1),F(:,:,3),AbsTol=0)
            testCase.verifyEqual(F(:,:,2),F(:,:,5),AbsTol=0)
        end

        function canonicalBasisPreservesPhysicalEndpointCoordinates(testCase)
            boundaryModes = testCase.standardBoundaryModes([1e-4 2e-4]);
            testCase.verifyEqual(boundaryModes.normalization,"endpointResponse")
            testCase.verifyEmpty(boundaryModes.boundaryDepth)
            testCase.verifyEqual(boundaryModes.canonicalEndpointLabels,["surface" "bottom"])
            testCase.verifyEqual(boundaryModes.boundaryResponseMatrix,repmat(eye(2),1,1,2),AbsTol=0)
            testCase.verifyFalse(any(ismember(["rotationName","rotationMatrix","h0"],string(properties(boundaryModes)))))
        end

        function constantStratificationMatchesIndependentExactBoundarySolve(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-1000 0];
            k = [1e-4 2.5e-4];
            z = linspace(zDomain(1),zDomain(2),101).';
            N2 = @(zIn) N0^2*ones(size(zIn));
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,k=k,surfaceBoundary=surfaceBoundary);
                basisSet = IMSolverSpectral(nEVP=72).solveGeostrophicZeroAPVModes(problem);
                F = basisSet.F(z);
                G = basisSet.G(z);
                for iK = 1:numel(k)
                    [FExact,GExact] = testCase.constantExactModes(z,k(iK),N0,f0,g,zDomain,surfaceBoundary);
                    testCase.verifyEqual(F(:,:,iK),FExact,RelTol=2e-10,AbsTol=2e-10)
                    testCase.verifyEqual(G(:,:,iK),GExact,RelTol=2e-10,AbsTol=2e-10)
                end
            end
        end

        function quadraticFormsMatchIndependentContinuousEnergy(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-900 0];
            k = 1.7e-4;
            N2 = @(z) N0^2*ones(size(z));
            z = linspace(zDomain(1),zDomain(2),12001).';
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,k=k,surfaceBoundary=surfaceBoundary);
                basisSet = IMSolverSpectral(nEVP=80).solveGeostrophicZeroAPVModes(problem);
                F = basisSet.F(z);
                G = basisSet.G(z);
                F = F(:,:,1);
                G = G(:,:,1);
                integrand = zeros(length(z),2,2);
                for i = 1:2
                    for j = 1:2
                        integrand(:,i,j) = k^2*F(:,i).*F(:,j) + (f0^2*N0^2/g^2)*G(:,i).*G(:,j);
                    end
                end
                expectedEnergy = squeeze(trapz(z,integrand,1))/(2*k^4);
                if surfaceBoundary == "freeSurface"
                    expectedEnergy = expectedEnergy + (f0^2/(2*g*k^4))*(F(end,:).'*F(end,:));
                end
                responses = testCase.endpointResponses(basisSet);
                responseSurface = responses(1,:,1);
                responseBottom = responses(2,:,1);
                expectedSurface = (f0^2/(2*g^2*k^4))*(responseSurface.'*responseSurface);
                expectedBottom = (f0^2/(2*g^2*k^4))*(responseBottom.'*responseBottom);

                testCase.verifyEqual(basisSet.energyMatrix(:,:,1),expectedEnergy,RelTol=2e-7,AbsTol=3e-8)
                testCase.verifyEqual(basisSet.surfaceBuoyancyMatrix(:,:,1),expectedSurface,RelTol=1e-11,AbsTol=1e-11)
                testCase.verifyEqual(basisSet.bottomBuoyancyMatrix(:,:,1),expectedBottom,RelTol=1e-11,AbsTol=1e-11)
                testCase.verifyEqual(basisSet.energyMatrix(:,:,1),basisSet.energyMatrix(:,:,1).',AbsTol=0)
                testCase.verifyEqual(basisSet.surfaceBuoyancyMatrix(:,:,1),basisSet.surfaceBuoyancyMatrix(:,:,1).',AbsTol=0)
                testCase.verifyEqual(basisSet.bottomBuoyancyMatrix(:,:,1),basisSet.bottomBuoyancyMatrix(:,:,1).',AbsTol=0)
                expectedGeneralized = expectedEnergy-0.03*expectedSurface+0.02*expectedBottom;
                testCase.verifyEqual(basisSet.generalizedEnergyMatrix(g0=-0.03,gd=0.02),expectedGeneralized,RelTol=2e-7,AbsTol=3e-8)
            end
        end

        function generalizedEnergyAcceptsPositiveNegativeAndZeroWeights(testCase)
            boundaryModes = testCase.standardBoundaryModes(1e-4);
            coefficientCases = [0 0; 0.03 0; -0.03 0.02; 0 -0.02];
            for iCase = 1:size(coefficientCases,1)
                g0 = coefficientCases(iCase,1);
                gd = coefficientCases(iCase,2);
                expected = boundaryModes.energyMatrix+g0*boundaryModes.surfaceBuoyancyMatrix+gd*boundaryModes.bottomBuoyancyMatrix;
                actual = boundaryModes.generalizedEnergyMatrix(g0=g0,gd=gd);
                testCase.verifyEqual(actual,expected,AbsTol=0)
                testCase.verifyEqual(actual,permute(actual,[2 1 3]),AbsTol=0)
            end
            testCase.verifyTrue(testCase.throwsAny(@() boundaryModes.generalizedEnergyMatrix(g0=Inf)))
            testCase.verifyTrue(testCase.throwsAny(@() boundaryModes.generalizedEnergyMatrix(gd=NaN)))
        end

        function absentEndpointCoefficientIsIgnored(testCase)
            N2 = @(z) ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,k=1,endpoints="surface");
            basisSet = IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(problem);
            testCase.verifyEqual(basisSet.bottomBuoyancyMatrix,zeros(1,1),AbsTol=0)
            testCase.verifyEqual(basisSet.generalizedEnergyMatrix(g0=3,gd=-1e9),basisSet.generalizedEnergyMatrix(g0=3,gd=0),AbsTol=0)
        end

        function generalizedEnergyOrthogonalityHasSignedBoundaryDepth(testCase)
            boundaryModes = testCase.standardBoundaryModes([1e-4 2e-4]);
            g0 = -.035; gd = .018;
            orthogonal = boundaryModes.orthogonalizeGeneralizedEnergy(g0=g0,gd=gd);
            Hg = boundaryModes.generalizedEnergyMatrix(g0=g0,gd=gd);
            for page = 1:numel(boundaryModes.k)
                C = orthogonal.boundaryResponseMatrix(:,:,page);
                target = diag(orthogonal.boundaryDepth(:,page)/(2*boundaryModes.k(page)^2));
                testCase.verifyEqual(C.'*C,eye(2),AbsTol=2e-11)
                testCase.verifyEqual(C.'*Hg(:,:,page)*C,target,AbsTol=2e-8)
                testCase.verifyLessThanOrEqual(max(diff(orthogonal.boundaryDepth(:,page))),0)
                testCase.verifyGreaterThanOrEqual(testCase.referenceCoefficients(C),zeros(1,size(C,2)))
            end
            testCase.verifyLessThan(max(orthogonal.metadata.orthogonality.relativeEnergyResidual,[],"all"),1e-12)
            testCase.verifyLessThan(max(orthogonal.metadata.orthogonality.endpointOrthogonalityError),1e-12)
        end

        function publicEndpointChoicesFollowTheManuscript(testCase)
            boundary = testCase.standardBoundaryModes(1e-4);
            methodNames = string(methods(boundary));
            testCase.verifyFalse(any(ismember(["rotateSurfaceBuoyancy","rotateBoundaryDepth","rotateWithPencil"],methodNames)))
            testCase.verifyTrue(ismember("orthogonalizeGeneralizedEnergy",methodNames))
            testCase.verifyFalse(ismember("rotationRightMetric",boundary.normalizationNames()))
            testCase.verifyTrue(testCase.throwsAny(@() boundary.orthogonalizeGeneralizedEnergy()))
            testCase.verifyTrue(testCase.throwsAny(@() boundary.orthogonalizeGeneralizedEnergy(g0=0)))
        end

        function rotationsApplyIdenticalColumnMapToFAndG(testCase)
            boundaryModes = testCase.standardBoundaryModes(1e-4);
            rotatedModes = boundaryModes.orthogonalizeGeneralizedEnergy(g0=-0.035,gd=0.02);
            z = linspace(-1000,0,31).';
            C = rotatedModes.boundaryResponseMatrix(:,:,1);
            testCase.verifyEqual(rotatedModes.F(z),boundaryModes.F(z)*C,RelTol=1e-11,AbsTol=1e-11)
            testCase.verifyEqual(rotatedModes.G(z),boundaryModes.G(z)*C,RelTol=1e-11,AbsTol=1e-11)
        end

        function aSingleChosenEndpointIsACompletePhysicalFamily(testCase)
            N2 = @(z) ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,k=1,endpoints="bottom");
            basis = IMSolverSpectral(nEVP=32).solveGeostrophicZeroAPVModes(problem);
            orthogonal = basis.orthogonalizeGeneralizedEnergy(g0=Inf,gd=.4);
            testCase.verifyEqual(orthogonal.nColumns,1)
            testCase.verifyEqual(orthogonal.boundaryResponseMatrix,1,AbsTol=0)
            testCase.verifyEqual(orthogonal.boundaryDepth,2*basis.k^2*basis.generalizedEnergyMatrix(g0=Inf,gd=.4),AbsTol=1e-12)
        end

        function finiteDifferenceSolverPreservesCanonicalResponses(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-1000 0];
            N2 = @(z) N0^2*ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),k=1e-4);
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1),zDomain(2),129).');
            basisSet = solver.solveGeostrophicZeroAPVModes(problem);
            responses = testCase.endpointResponses(basisSet);
            testCase.verifyEqual(responses(:,:,1),eye(2),AbsTol=2e-8)
            mismatchedSolver = IMSolverFiniteDifference(z=linspace(-900,0,129).');
            testCase.verifyError(@() mismatchedSolver.solveGeostrophicZeroAPVModes(problem),"IMSolverFiniteDifference:DomainMismatch")
        end

        function exponentialStratificationConvergesAndSatisfiesIndependentResidual(testCase)
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            g = 9.81;
            zDomain = [-4000 0];
            k = 2*pi/1e6;
            N2 = @(z) N0^2*exp(2*z/b);
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,k=k);
            coarse = IMSolverSpectral(nEVP=96).solveGeostrophicZeroAPVModes(problem);
            refined = IMSolverSpectral(nEVP=192).solveGeostrophicZeroAPVModes(problem);
            z = linspace(zDomain(1),zDomain(2),4001).';
            FCoarse = coarse.F(z);
            FRefined = refined.F(z);
            testCase.verifyEqual(FCoarse,FRefined,RelTol=2e-8,AbsTol=2e-8)

            F = FRefined(:,:,1);
            dz = z(2)-z(1);
            p = f0^2./N2(z);
            for iMode = 1:2
                Fz = gradient(F(:,iMode),dz);
                residual = gradient(p.*Fz,dz)-k^2*F(:,iMode);
                interior = 4:(length(z)-3);
                scale = max(k^2*max(abs(F(interior,iMode))),1e-14);
                testCase.verifyLessThan(max(abs(residual(interior)))/scale,2e-3)
            end
        end

        function surfaceBoundaryMatchesAPVMetadataContract(testCase)
            N2 = @(z) ones(size(z));
            for surfaceBoundary = ["freeSurface" "rigidLid"]
                zeroProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,k=1,surfaceBoundary=surfaceBoundary);
                apvProblem = IMInternalModes.geostrophicAPVModes(rotationRate=0,N2=N2,zDomain=[-1 0],g=2,g0=3,gd=4,surfaceBoundary=surfaceBoundary);
                testCase.verifyEqual(zeroProblem.surfaceBoundary,apvProblem.parameters.surfaceBoundary)
            end
        end
    end

    methods (Access = private)
        function basisSet = standardBoundaryModes(~,k)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-1000 0];
            N2 = @(z) N0^2*ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=9.81,k=k);
            basisSet = IMSolverSpectral(nEVP=80).solveGeostrophicZeroAPVModes(problem);
        end

        function responses = endpointResponses(~,basisSet)
            nEndpoints = numel(basisSet.endpoints);
            nK = numel(basisSet.k);
            FSurface = reshape(basisSet.F(basisSet.zDomain(2)),nEndpoints,nK);
            GSurface = reshape(basisSet.G(basisSet.zDomain(2)),nEndpoints,nK);
            GBottom = reshape(basisSet.G(basisSet.zDomain(1)),nEndpoints,nK);
            surfaceResponse = GSurface;
            if basisSet.surfaceBoundary == "freeSurface"
                surfaceResponse = surfaceResponse-FSurface;
            end
            responses = zeros(2,nEndpoints,nK);
            responses(1,:,:) = reshape(surfaceResponse,1,nEndpoints,nK);
            responses(2,:,:) = reshape(GBottom,1,nEndpoints,nK);
        end

        function rows = endpointRows(~,endpoints)
            rows = zeros(1,numel(endpoints));
            rows(endpoints == "surface") = 1;
            rows(endpoints == "bottom") = 2;
        end

        function [F,G] = constantExactModes(~,z,k,N0,f0,g,zDomain,surfaceBoundary)
            depth = diff(zDomain);
            m = k*N0/f0;
            c = cosh(m*depth);
            s = sinh(m*depth);
            bottomRow = [0 -g*m/N0^2];
            surfaceRow = -(g/N0^2)*m*[s c];
            if surfaceBoundary == "freeSurface"
                surfaceRow = surfaceRow-[c s];
            end
            responseMatrix = [surfaceRow;bottomRow];
            coefficients = responseMatrix\eye(2);
            shiftedZ = z-zDomain(1);
            basisValues = [cosh(m*shiftedZ) sinh(m*shiftedZ)];
            basisDerivatives = m*[sinh(m*shiftedZ) cosh(m*shiftedZ)];
            F = basisValues*coefficients;
            G = -(g/N0^2)*basisDerivatives*coefficients;
        end

        function didThrow = throwsAny(~,functionHandle)
            didThrow = false;
            try
                functionHandle();
            catch
                didThrow = true;
            end
        end

        function coefficients = referenceCoefficients(~,matrix)
            coefficients = zeros(1,size(matrix,2));
            for iColumn = 1:size(matrix,2)
                [~,referenceIndex] = max(abs(matrix(:,iColumn)));
                coefficients(iColumn) = matrix(referenceIndex,iColumn);
            end
        end
    end
end
