classdef IMBackgroundProfileTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function profileSourceAndSurfaceDensityAreUnambiguous(testCase)
            zDomain = [-1000 0];
            N2 = @(z) 1e-4*ones(size(z));
            rho = @(z) 1030*ones(size(z));
            invalidRho = @(z) zeros(size(z));
            testCase.verifyError(@() IMBackground(zDomain=zDomain,surfaceCondition="rigidLid"), ...
                "IMBackground:InvalidProfileSource");
            testCase.verifyError(@() IMBackground(N2=N2,rho=rho,zDomain=zDomain,surfaceCondition="rigidLid"), ...
                "IMBackground:InvalidProfileSource");
            testCase.verifyError(@() IMBackground(rho=invalidRho,zDomain=zDomain,surfaceCondition="rigidLid"), ...
                "IMBackground:InvalidSurfaceDensity");
            testCase.verifyError(@() IMBackground(rho=rho,rho0=1025,zDomain=zDomain,surfaceCondition="rigidLid"), ...
                "IMBackground:ConflictingSurfaceDensity");

            background = IMBackground(rho=rho,rho0=1030,zDomain=zDomain,surfaceCondition="rigidLid");
            testCase.verifyEqual(background.rho0,1030);
            testCase.verifyEqual(background.N2(zDomain),zeros(size(zDomain)),AbsTol=1e-16);
        end

        function bothSourceFormsMatchTranslatedExponentialProfiles(testCase)
            zSurface = 125;
            zDomain = [-875 zSurface];
            rho0 = 1031.5;
            g = 10.4;
            surfaceN2 = 1.2e-4;
            scale = 600;
            N2 = @(z) surfaceN2*exp((z-zSurface)/scale);
            rho = @(z) rho0*(1-(surfaceN2*scale/g)*expm1((z-zSurface)/scale));
            z = reshape(zSurface-[0 950 200 350 1000 25 300 910 50 750 600 50],2,3,2);

            fromN2 = IMBackground(N2=N2,rho0=rho0,zDomain=zDomain,g=g,surfaceCondition="freeSurface");
            fromDensity = IMBackground(rho=rho,zDomain=zDomain,g=g,surfaceCondition="freeSurface");
            for caseBackground = {fromN2,fromDensity}
                background = caseBackground{1};
                testCase.verifyClass(background.N2,"function_handle");
                testCase.verifyClass(background.rho,"function_handle");
                testCase.verifyEqual(background.rho0,rho0,AbsTol=2e-12);
                testCase.verifyEqual(background.rho(zSurface),background.rho0,AbsTol=2e-12);
                testCase.verifyEqual(background.N2(z),N2(z),RelTol=3e-10,AbsTol=2e-15);
                testCase.verifyEqual(background.rho(z),rho(z),AbsTol=3e-11);
            end
        end

        function constantN2UsesDefaultSurfaceDensity(testCase)
            zDomain = [-400 0];
            N2 = 7e-5;
            background = IMBackground(N2=@(z) N2*ones(size(z)),zDomain=zDomain,surfaceCondition="rigidLid");
            z = linspace(zDomain(1),zDomain(2),7).';

            testCase.verifyEqual(background.rho0,1025);
            testCase.verifyEqual(background.rho(z),1025*(1-N2*z/background.g),AbsTol=2e-12);
        end

        function spectralDensityPreparationIsLazyAndReused(testCase)
            calls = 0;
            background = IMBackground(N2=@countedProfile,zDomain=[-400 0],surfaceCondition="rigidLid");
            testCase.verifyEqual(calls,0);
            density = background.rho;
            preparationCalls = calls;
            testCase.verifyGreaterThan(preparationCalls,0);
            testCase.verifyEqual(background.rho(-100),density(-100));
            testCase.verifyEqual(background.rho(0),background.rho0);
            testCase.verifyEqual(calls,preparationCalls);

            function values = countedProfile(z)
                calls = calls+1;
                values = 1e-4*ones(size(z));
            end
        end

        function analyticalFactoriesCarryTheDensityAnchor(testCase)
            zDomain = [-900 100];
            rho0 = 1035;
            N0 = .01;
            b = 500;
            constant = IMBackground.constantStratification(N0=N0,rho0=rho0,zDomain=zDomain,surfaceCondition="rigidLid");
            exponential = IMBackground.exponentialStratification( ...
                N0=N0,b=b,rho0=rho0,zDomain=zDomain,surfaceCondition="freeSurface");
            depth = diff(zDomain);

            testCase.verifyEqual(constant.profileParameters,struct("N0",N0));
            testCase.verifyEqual(exponential.profileParameters,struct("N0",N0,"b",b));
            testCase.verifyEqual(constant.rho(zDomain),rho0*[1+N0^2*depth/constant.g 1],AbsTol=3e-12);
            rhoBottom = rho0*(1-(N0^2*b/(2*exponential.g))*expm1(-2*depth/b));
            testCase.verifyEqual(exponential.rho(zDomain),[rhoBottom rho0],AbsTol=3e-12);
        end

        function neutralLayerRetainsItsSuppliedN2AndIntegratedDensity(testCase)
            zDomain = [-1000 0];
            zJoin = -250;
            scale = 1000;
            N2Amplitude = 1e-4;
            N2 = @(z) N2Amplitude*(max(zJoin-z,0)/scale).^2;
            background = IMBackground(N2=N2,zDomain=zDomain,surfaceCondition="rigidLid");
            z = [zDomain(1) -500 zJoin-1 zJoin zJoin+1 zDomain(2)];
            rho = background.rho0*(1+N2Amplitude*max(zJoin-z,0).^3/(3*background.g*scale^2));

            testCase.verifyEqual(background.N2(z(z>=zJoin)),zeros(1,3));
            testCase.verifyEqual(background.rho(z),rho,AbsTol=3e-10);
        end

        function densitySourceFeedsWaveSolversAtIndependentResolutions(testCase)
            zSurface = 100;
            D = 800;
            zDomain = zSurface+[-D 0];
            g = 10;
            rho0 = 1030;
            N2 = @(z) 1e-4*(1+.4*(z-zSurface)/D);
            rho = @(z) rho0*(1-(1e-4/g)*((z-zSurface)+.2*(z-zSurface).^2/D));
            background = IMBackground(rho=rho,zDomain=zDomain,g=g,latitude=-30,surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
            candidate = IMSolverSpectral(nEVP=48,coordinateKind="wkb").solveEVP(evp,nModes=3);
            reference = IMSolverSpectral(nEVP=72,coordinateKind="z").solveEVP(evp,nModes=3);

            testCase.verifyEqual(candidate.modeNumber,reference.modeNumber);
            testCase.verifyEqual(candidate.h,reference.h,RelTol=1e-6);
            for basis = {candidate,reference}
                solved = basis{1};
                z = solved.solver.zNative;
                testCase.verifyEqual(solved.evp.background.rho0,rho0,AbsTol=2e-12);
                testCase.verifyEqual(solved.evp.background.rho(z),rho(z),AbsTol=3e-11);
                testCase.verifyEqual(solved.evp.N2(z),N2(z),RelTol=3e-10,AbsTol=2e-15);
            end
        end
    end
end
