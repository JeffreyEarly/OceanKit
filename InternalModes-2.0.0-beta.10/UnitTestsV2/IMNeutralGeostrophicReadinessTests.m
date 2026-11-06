classdef IMNeutralGeostrophicReadinessTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function independentQuadratureIntegratesThePrescribedProfile(testCase)
            profile = neutralGeostrophicProfile();
            [z,weights] = neutralGeostrophicQuadrature(profile.joins,257);
            testCase.verifyGreaterThan(min(weights),0)
            testCase.verifyEqual(sum(weights),1000,AbsTol=2e-12)
            testCase.verifyEqual(sum(weights.*profile.N2(z)),profile.integratedN2,AbsTol=2e-16)
            testCase.verifyEqual(profile.N2([-200;-100;0]),zeros(3,1),AbsTol=0)
            testCase.verifyEqual(profile.N2([-1000;-500;-400]),(5e-3)^2*ones(3,1),AbsTol=0)
        end

        function independentConstantReferenceMatchesTheAnalyticalAPVFamily(testCase)
            N0 = 5e-3;
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) N0^2*ones(size(z)),zDomain=[-1000 0],g0=0,gd=-.01,surfaceBoundary="freeSurface");
            reference = neutralGeostrophicReference(evp,nModes=4,joins=[-1000 0]);
            analytical = IMConstantStratificationSolution(N0=N0,zDomain=[-1000 0]).internalModes(evp,nModes=4);
            [z,weights] = neutralGeostrophicQuadrature([-1000 0],257);
            values = reference.evaluate(z);
            % The analytical family's normalization uses its own integration
            % rule. Fit one shared amplitude and test both physical shapes.
            amplitudes = sum(weights.*values.F.*analytical.F(z),1)./sum(weights.*values.F.^2,1);
            testCase.verifyEqual(reference.h,analytical.h,RelTol=2e-9)
            testCase.verifyEqual(values.F.*amplitudes,analytical.F(z),RelTol=2e-8,AbsTol=2e-8)
            testCase.verifyEqual(values.G.*amplitudes,analytical.G(z),RelTol=2e-8,AbsTol=2e-8)
            testCase.verifyEqual(values.F.'*(weights.*values.F)/1000,eye(4),AbsTol=2e-10)
            testCase.verifyEqual(reference.scan.negativeFound,1)
        end

        function neutralAPVReferencesRefineIndependentlyAndHaveAnAffineCap(testCase)
            p = neutralGeostrophicProfile();
            evp = IMInternalModes.geostrophicAPVModes(N2=p.N2,zDomain=p.zDomain,g0=0,gd=-.01,surfaceBoundary="freeSurface");
            coarse = neutralGeostrophicReference(evp,nModes=4,tolerance=1e-11);
            tight = neutralGeostrophicReference(evp,nModes=4,tolerance=1e-13,rootBrackets=coarse.rootBrackets);
            [z,weights] = neutralGeostrophicQuadrature(p.joins,257);
            agreement = neutralGeostrophicComparison(tight,coarse,z,weights);
            testCase.verifyLessThan(agreement.maximumError,1e-7)
            cap = tight.evaluate((-200:20:0).');
            testCase.verifyEqual(cap.F,ones(size(cap.F,1),1)*cap.F(1,:),AbsTol=1e-10)
            expectedG = cap.G(1,:)+(0:20:200).'.*tight.inverseDepth.*cap.F(1,:);
            testCase.verifyEqual(cap.G,expectedG,RelTol=1e-12,AbsTol=1e-10)
        end

        function zeroAPVReferenceHasCanonicalEndpointResponses(testCase)
            p = neutralGeostrophicProfile();
            problem = IMGeostrophicZeroAPVModes(N2=p.N2,zDomain=p.zDomain,k=2*pi/1e3,latitude=45);
            reference = neutralGeostrophicReference(problem);
            values = reference.evaluate(p.zDomain.');
            response = [values.G(end,:)-values.F(end,:);values.G(1,:)];
            testCase.verifyEqual(response,eye(2),AbsTol=1e-10)
        end
    end
end
