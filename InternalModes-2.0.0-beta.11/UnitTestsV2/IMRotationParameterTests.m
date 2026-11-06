classdef IMRotationParameterTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function defaultsAndContextUsePlanetaryInputs(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyEqual(evp.rotationRate,7.2921e-5);
            testCase.verifyEqual(evp.latitude,33);
            testCase.verifyEqual(evp.f0,2*7.2921e-5*sind(33));
            solver = IMSolverSpectral(nEVP=17).configuredForEVP(evp);
            context = evp.contextForSolver(solver);
            testCase.verifyEqual(context.rotationRate,evp.rotationRate);
            testCase.verifyEqual(context.latitude,evp.latitude);
            testCase.verifyEqual(context.f0,evp.f0);
            testCase.verifyEqual(evp.parameters.f0,evp.f0);
        end

        function hemisphereSignDoesNotChangeTheBandwidthGrid(testCase)
            N2 = @(z) 1e-8*exp(2*z/700);
            north = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            south = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,latitude=-45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyEqual(south.f0,-north.f0);
            first = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth").configuredForEVP(north);
            second = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth").configuredForEVP(south);
            testCase.verifyEqual(first.zNative,second.zNative,AbsTol=1e-11);
        end

        function planetaryInputsAreValidated(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            testCase.verifyError(@() IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=-1,surfaceBoundary=IMBoundaryCondition.dirichlet()),"MATLAB:validators:mustBeNonnegative");
            testCase.verifyError(@() IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],latitude=91,surfaceBoundary=IMBoundaryCondition.dirichlet()),"MATLAB:validators:mustBeLessThanOrEqual");
            testCase.verifyError(@() IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],latitude=-91,surfaceBoundary=IMBoundaryCondition.dirichlet()),"MATLAB:validators:mustBeGreaterThanOrEqual");
            testCase.verifyError(@() IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=NaN,surfaceBoundary=IMBoundaryCondition.dirichlet()),"MATLAB:validators:mustBeFinite");
        end

        function zeroRotationRetainsTheOriginalExponentialStretch(testCase)
            D = 1000; b = 700; N0 = .01;
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) N0^2*exp(2*z/b),zDomain=[-D 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyEqual(evp.f0,0);
            z = linspace(-D,0,111);
            for kind = ["wkb", "bandwidth"]
                solver = IMSolverSpectral(nEVP=65,coordinateKind=kind).configuredForEVP(evp);
                testCase.verifyEqual(solver.xOfZ(z),N0*b*(exp(z/b)-exp(-D/b)),AbsTol=5e-13);
                testCase.verifyEqual(solver.qReference,N0*exp(solver.zReference/b),RelTol=2e-12);
                testCase.verifyEqual(solver.qzReference,N0/b*exp(solver.zReference/b),RelTol=2e-11);
            end
        end

        function classicalWKBIgnoresRotationWithoutChangingPhysics(testCase)
            D = 1000; b = 700; N0 = 1e-4;
            N2 = @(z) N0^2*exp(2*z/b);
            nonrotating = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-D 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            rotating = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-D 0],latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            classical = IMSolverSpectral(nEVP=65,coordinateKind="wkb").configuredForEVP(rotating);
            reference = classical.configuredForEVP(nonrotating);
            bandwidth = IMSolverSpectral(nEVP=65,coordinateKind="bandwidth").configuredForEVP(rotating);
            z = linspace(-D,0,111);
            testCase.verifyEqual(classical.qReference,N0*exp(classical.zReference/b),RelTol=2e-12);
            testCase.verifyEqual(classical.qzReference,N0/b*exp(classical.zReference/b),RelTol=2e-11);
            testCase.verifyEqual(classical.xOfZ(z),N0*b*(exp(z/b)-exp(-D/b)),AbsTol=5e-15);
            testCase.verifyEqual(classical.zNative,reference.zNative,AbsTol=1e-11);
            testCase.verifyGreaterThan(norm(classical.zNative-bandwidth.zNative),1);
            classicalContext = classical.context();
            bandwidthContext = bandwidth.context();
            testCase.verifyEqual(classicalContext.coordinateKind,"wkb");
            testCase.verifyEqual(bandwidthContext.coordinateKind,"bandwidth");
            testCase.verifyEqual(rotating.f0,2*rotating.rotationRate*sind(45));
            testCase.verifyGreaterThan(abs(rotating.f0),0);
        end

        function classicalWKBRejectsNeutralColumnsWithRotation(testCase)
            neutral = IMInternalModes.hydrostaticGModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            classical = IMSolverSpectral(nEVP=33,coordinateKind="wkb");
            testCase.verifyError(@() classical.configuredForEVP(neutral),"IMSolverSpectral:InvalidCoordinate");
        end

        function neutralLayerHasInvertibleMapAndPhysicalCalculus(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth").configuredForEVP(evp);
            q = 2*evp.rotationRate;
            testCase.verifyEqual(solver.qReference,q*ones(size(solver.qReference)),AbsTol=1e-17);
            z = linspace(-1000,0,51);
            testCase.verifyEqual(solver.xOfZ(z),q*(z+1000),AbsTol=1e-13);
            testCase.verifyEqual(solver.zOfX(solver.xOfZ(z)),z,AbsTol=2e-11);
            [nodes,weights,Dz] = solver.nativeDifferentiationRule(evp.zDomain);
            [~,~,Dzz] = solver.nativeDifferentiationRule(evp.zDomain,2);
            testCase.verifyEqual(Dz*(nodes.^2),2*nodes,AbsTol=2e-8);
            testCase.verifyEqual(Dzz*(nodes.^2),2*ones(size(nodes)),AbsTol=2e-8);
            testCase.verifyTrue(all(weights > 0));
            testCase.verifyEqual(sum(weights),1000,AbsTol=1e-10);
        end

        function rotatingExponentialHasConsistentJacobianDerivative(testCase)
            b = 700; N2 = @(z) 1e-8*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],latitude=33,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=65,coordinateKind="bandwidth").configuredForEVP(evp);
            horizontal = 2*evp.rotationRate*cosd(evp.latitude);
            w = N2(solver.zReference)+horizontal^2-evp.f0^2;
            B = hypot(w,2*evp.f0*horizontal);
            expectedDerivative = w.*N2(solver.zReference)./(b*B.^(3/2));
            testCase.verifyEqual(solver.qzReference,expectedDerivative,AbsTol=1e-15);
            [z,~,Dz] = solver.nativeDifferentiationRule(evp.zDomain);
            testCase.verifyEqual(Dz*(z.^2),2*z,AbsTol=1e-7);
        end

        function rotationReconfigurationRebuildsWithoutMutatingOriginal(testCase)
            N2 = @(z) 1e-8*exp(2*z/700);
            firstEVP = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            secondEVP = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1000 0],latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            first = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth").configuredForEVP(firstEVP);
            original = first.xReference;
            second = first.configuredForEVP(secondEVP);
            testCase.verifyEqual(first.xReference,original);
            testCase.verifyGreaterThan(norm(first.zNative-second.zNative),1);
            testCase.verifyGreaterThan(min(second.qReference),0);
        end

        function zeroAPVConfigurationReceivesRotationMetadata(testCase)
            N2 = @(z) 1e-8*ones(size(z));
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,latitude=-45);
            solver = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth").configuredForGeostrophicZeroAPVModes(problem);
            horizontal = 2*problem.rotationRate*cosd(problem.latitude);
            expected = sqrt(hypot(1e-8+horizontal^2-problem.f0^2,2*problem.f0*horizontal));
            testCase.verifyEqual(solver.qReference,expected*ones(size(solver.qReference)),AbsTol=1e-17);
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-1000 0],k=1e-3,rotationRate=0),"IMGeostrophicZeroAPVModes:InvalidCoriolis");
        end

        function invalidAndDegenerateMonitorsRemainErrors(testCase)
            base = IMSolverSpectral(nEVP=33,coordinateKind="bandwidth");
            negative = IMInternalModes.hydrostaticGModes(N2=@(z) -1e-8*ones(size(z)),zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyError(@() base.configuredForEVP(negative),"IMSolverSpectral:InvalidCoordinate");
            neutral = IMInternalModes.hydrostaticGModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],rotationRate=0,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyError(@() base.configuredForEVP(neutral),"IMSolverSpectral:InvalidCoordinate");
            rotationRate = 7.2921e-5;
            degenerate = IMInternalModes.hydrostaticGModes(N2=@(z) (2*rotationRate)^2*ones(size(z)),zDomain=[-1000 0],rotationRate=rotationRate,latitude=90,surfaceBoundary=IMBoundaryCondition.dirichlet());
            testCase.verifyError(@() base.configuredForEVP(degenerate),"IMSolverSpectral:InvalidCoordinate");
        end
    end
end
