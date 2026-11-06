classdef IMSmoothedBandwidthCoordinateTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addStudyPaths(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repoRoot,"ExamplesV2")));
        end
    end

    methods (Test)
        function constantControlKeepsTheGridAndPhysicalSpectrum(testCase)
            N2 = @(z) 2.5e-5*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            ordinary = IMSolverSpectral(nEVP=65,coordinateKind="wkb").configuredForEVP(evp);
            smoothed = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=65).configuredForEVP(evp);
            testCase.verifyEqual(smoothed.zNative,ordinary.zNative,AbsTol=2e-9);
            % Constant coordinates give the same physical grid while their
            % coordinate scales retain distinct N and bandwidth definitions.
            horizontalCoriolis = 2*evp.rotationRate*cosd(evp.latitude);
            expectedBandwidth = sqrt(hypot(N2(smoothed.zReference)+horizontalCoriolis^2-evp.f0^2,2*evp.f0*horizontalCoriolis));
            testCase.verifyEqual(ordinary.qReference,sqrt(N2(ordinary.zReference)),RelTol=2e-11);
            testCase.verifyEqual(smoothed.qReference,expectedBandwidth,RelTol=2e-11);
            testCase.verifyEqual(smoothed.N2(smoothed.zReference),N2(smoothed.zReference));
            basis = smoothed.solveEVP(evp,nModes=4);
            expectedH = 2.5e-5*1000^2./(evp.g*pi^2*(1:4).^2);
            testCase.verifyEqual(basis.h,expectedH,RelTol=2e-9);
        end

        function sinusoidalSourceMatchesTheExactNeumannFilter(testCase)
            zDomain = [-1000 0];
            smoothingLength = 150;
            k = pi/1000;
            meanN2 = 2e-5;
            amplitude = 1e-5;
            N2 = @(z) meanN2+amplitude*cos(k*(z-zDomain(1)));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=65,smoothingLength=smoothingLength).configuredForEVP(evp);
            z = solver.zReference;
            filteredAmplitude = amplitude/(1+(smoothingLength*k)^2);
            expectedQ = meanN2+filteredAmplitude*cos(k*(z-zDomain(1)));
            expectedQz = -k*filteredAmplitude*sin(k*(z-zDomain(1)));
            testCase.verifyEqual(feval(solver.smoothedBandwidthSquared,z),expectedQ,RelTol=2e-10);
            testCase.verifyEqual(solver.qReference.^2,expectedQ,RelTol=2e-10);
            testCase.verifyEqual(solver.qzReference,expectedQz./(2*sqrt(expectedQ)),AbsTol=2e-13);
            testCase.verifyEqual(solver.qzReference([1 end]),[0;0],AbsTol=2e-13);
        end

        function rotatingSourceSatisfiesIndependentPhysicalDepthResidual(testCase)
            zDomain = [-1000 0];
            N2 = @(z) 1e-8*(1+0.7*cos(pi*(z+1000)/1000));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,latitude=33,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=65).configuredForEVP(evp);
            horizontalCoriolis = 2*evp.rotationRate*cosd(evp.latitude);
            z = linspace(zDomain(1),zDomain(2),301).';
            B = hypot(N2(z)+horizontalCoriolis^2-evp.f0^2,2*evp.f0*horizontalCoriolis);
            Q = solver.smoothedBandwidthSquared;
            residual = feval(Q,z)-solver.smoothingLength^2*feval(diff(Q,2),z)-B;
            testCase.verifyLessThan(max(abs(residual))/max(B),2e-8);
            testCase.verifyEqual(feval(diff(Q),zDomain),[0 0],AbsTol=2e-15);
            testCase.verifyEqual(solver.N2(z),N2(z));
        end

        function partialNeutralLayerIsInvertibleWithoutRotation(testCase)
            N2 = @IMSmoothedBandwidthCoordinateTests.neutralProfile;
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=257).configuredForEVP(evp);
            z = linspace(-1000,0,151).';
            testCase.verifyGreaterThan(min(solver.qReference),0);
            testCase.verifyTrue(all(diff(solver.xReference) > 0));
            testCase.verifyEqual(solver.zOfX(solver.xOfZ(z)),z,AbsTol=3e-10);
            testCase.verifyEqual(solver.N2([-200;-100;0]),zeros(3,1));
            testCase.verifyEqual(solver.N2(z),N2(z));
            [nodes,weights,Dz] = solver.nativeDifferentiationRule(evp.zDomain);
            [~,~,Dzz] = solver.nativeDifferentiationRule(evp.zDomain,2);
            testCase.verifyEqual(Dz*nodes,ones(size(nodes)),AbsTol=1e-5);
            testCase.verifyEqual(Dzz*nodes,zeros(size(nodes)),AbsTol=2e-6);
            testCase.verifyTrue(all(weights > 0));
            testCase.verifyEqual(sum(weights),1000,AbsTol=1e-6);
            testCase.verifyEqual(weights.'*nodes,-500000,AbsTol=1e-3);
        end

        function invalidPhysicalProfilesAndDegenerateBandwidthAreRejected(testCase)
            solver = IMBandwidthComparison.SmoothedBandwidthSolver();
            neutral = IMInternalModes.hydrostaticGModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            negative = IMInternalModes.hydrostaticGModes(N2=@(z) -1e-8*ones(size(z)),zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyError(@() solver.configuredForEVP(neutral),"IMSolverSpectral:InvalidCoordinate");
            testCase.verifyError(@() solver.configuredForEVP(negative),"IMSolverSpectral:InvalidCoordinate");
            testCase.verifyError(@() IMBandwidthComparison.SmoothedBandwidthSolver(smoothingLength=0),"MATLAB:validators:mustBePositive");
            testCase.verifyError(@() IMBandwidthComparison.SmoothedBandwidthSolver(smoothingLength=Inf),"MATLAB:validators:mustBeFinite");
        end

        function canonicalReconfigurationRebuildsWithoutMutatingOriginal(testCase)
            N2 = @(z) 1e-8*exp(2*z/700);
            firstEVP = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            secondEVP = IMInternalModes.hydrostaticGModes(N2=@(z) 4e-8*ones(size(z)),zDomain=[-1500 0],latitude=-45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            first = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=33).configuredForEVP(firstEVP);
            originalX = first.xReference;
            originalQ = feval(first.smoothedBandwidthSquared,first.zReference);
            second = first.configuredForEVP(secondEVP);
            fresh = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=33).configuredForEVP(secondEVP);
            testCase.verifyEqual(first.xReference,originalX);
            testCase.verifyEqual(feval(first.smoothedBandwidthSquared,first.zReference),originalQ);
            testCase.verifyEqual(first.N2(first.zReference),N2(first.zReference));
            testCase.verifyEqual(second.zNative,fresh.zNative,AbsTol=1e-10);
            testCase.verifyEqual(second.xReference,fresh.xReference,AbsTol=1e-13);
            testCase.verifyEqual(second.N2(second.zReference),4e-8*ones(size(second.zReference)));
        end

        function zeroAPVConfigurationRebuildsAndRestoresPhysicalStratification(testCase)
            N2 = @(z) 1e-8*exp(2*z/700);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            original = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=33).configuredForEVP(evp);
            zeroN2 = @(z) 2e-8*ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=zeroN2,zDomain=[-1500 0],k=1e-3,latitude=-45);
            rebuilt = original.configuredForGeostrophicZeroAPVModes(problem);
            fresh = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=33).configuredForGeostrophicZeroAPVModes(problem);
            testCase.verifyEqual(rebuilt.zNative,fresh.zNative,AbsTol=1e-10);
            testCase.verifyEqual(rebuilt.N2(rebuilt.zReference),zeroN2(rebuilt.zReference));
            returned = rebuilt.configuredForEVP(evp);
            testCase.verifyEqual(returned.xReference,original.xReference,AbsTol=1e-13);
            testCase.verifyEqual(returned.N2(returned.zReference),N2(returned.zReference));
            testCase.verifyEqual(original.N2(original.zReference),N2(original.zReference));
        end
    end

    methods (Static, Access = private)
        function N2 = neutralProfile(z)
            t = min(max((-z-200)/200,0),1);
            N2 = 0.005^2*(10*t.^3-15*t.^4+6*t.^5);
        end
    end
end
