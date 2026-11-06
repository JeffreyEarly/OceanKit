classdef IMModeNeutralConstructionTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function fullyNeutralRigidHydrostaticFKeepsItsBarotropicMode(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="rigidLid");
            problem = IMInternalModes.hydrostaticFModes(background);
            basis = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=1);
            exact = IMSolverAnalytical().solveEVP(problem,nModes=1);
            testCase.verifyEqual(basis.nColumns,1);
            testCase.verifyEqual(basis.modeInfo.branch,"barotropic");
            testCase.verifyEqual(basis.eigenvalues,0);
            testCase.verifyEqual(basis.h,Inf);
            testCase.verifyEqual(basis.normalization,"rmsF");
            z = linspace(-1000,0,13).';
            testCase.verifyEqual(basis.F(z),ones(size(z)),AbsTol=2e-13);
            testCase.verifyEqual(basis.G(z),zeros(size(z)),AbsTol=2e-13);
            testCase.verifyEqual(exact.F(z),basis.F(z),AbsTol=2e-13);
            testCase.verifyEqual(exact.G(z),basis.G(z),AbsTol=2e-13);
            for backend = {IMSolverAnalytical(),IMSolverSpectral(nEVP=32)}
                explicit = backend{1}.solveEVP(problem,branch="barotropic",nModes=1);
                mixedEmpty = backend{1}.solveEVP(problem,branch=["barotropic","positive"],nModes=[1 0]);
                testCase.verifyEqual(explicit.F(z),ones(size(z)),AbsTol=2e-13);
                testCase.verifyEqual(mixedEmpty.F(z),ones(size(z)),AbsTol=2e-13);
            end
        end

        function rigidHydrostaticFUsesTheSameRegularAnalyticalAndNumericalPair(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            exact = IMSolverAnalytical().solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
            numerical = IMSolverSpectral(nEVP=64).solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
            testCase.verifyEqual(numerical.nColumns,4);
            testCase.verifyEqual(numerical.modeInfo.branch,["barotropic";"positive";"positive";"positive"]);
            testCase.verifyEqual(numerical.columnLabels,exact.columnLabels);
            testCase.verifyEqual(numerical.h(2:end),exact.h(2:end),RelTol=2e-8);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(numerical.F(z),exact.F(z),AbsTol=2e-9);
            testCase.verifyEqual(numerical.G(z),exact.G(z),AbsTol=2e-9);
        end

        function balancedExplicitCountsRetainRequestedSpecialSectors(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            inactive = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                apv = solver{1}.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=inactive.g0,gd=inactive.gd),nModes=1,branch="zeroEigenvalue");
                hydrostatic = solver{1}.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=1);
                testCase.verifyEqual(apv.nColumns,1);
                testCase.verifyEqual(hydrostatic.nColumns,1);
                testCase.verifyEqual(apv.eigenvalues,0);
                testCase.verifyEqual(hydrostatic.eigenvalues,0);
                testCase.verifyEqual(apv.modeInfo.branch,"zeroEigenvalue");
                testCase.verifyEqual(hydrostatic.modeInfo.branch,"barotropic");
            end
            neumann = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            mda = IMSolverSpectral(nEVP=48).solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=neumann.g0,gd=neumann.gd),nModes=1,branch="constant");
            testCase.verifyEqual(mda.nColumns,1);
            testCase.verifyEqual(mda.modeInfo.branch,"constant");
            testCase.verifyEqual(mda.eigenvalues,0);
            G = mda.G([-1000;-500;0]);
            testCase.verifyEqual(G,repmat(G(1,:),3,1),AbsTol=2e-12);
        end

        function zeroCountsKeepEmptyEVPSpectraAvailable(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                wave = solver{1}.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=0);
                inertial = solver{1}.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=0),nModes=0);
                hydrostatic = solver{1}.solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=0);
                for basis = {wave,inertial,hydrostatic}
                    testCase.verifyEqual(basis{1}.nColumns,0);
                    testCase.verifyEmpty(basis{1}.h);
                    testCase.verifyEmpty(basis{1}.eigenvalues);
                    testCase.verifySize(basis{1}.F([-1000;0]),[2 0]);
                    testCase.verifySize(basis{1}.G([-1000;0]),[2 0]);
                end
            end
        end

        function emptyEndpointInventoryHasNoFictitiousSpectrum(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=32)}
                problem = IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints=strings(1,0));
                basis = solver{1}.solveGeostrophicZeroAPVModes(problem);
                testCase.verifyEqual(basis.nColumns,0);
                testCase.verifyEmpty(basis.columnLabels);
                testCase.verifySize(basis.F([-1000;0]),[2 0]);
                testCase.verifySize(basis.G([-1000;0]),[2 0]);
                testCase.verifyFalse(isprop(basis,"h"));
                testCase.verifyFalse(isprop(basis,"eigenvalues"));
            end
        end

        function numericalAndExactEndpointsKeepTheSameExplicitResponses(testCase)
            background = IMBackground.constantStratification(N0=.005,zDomain=[-1000 0],surfaceCondition="freeSurface");
            z = linspace(-1000,0,41).';
            endpoints = {["surface","bottom"],"surface","bottom"};
            for iEndpoints = 1:numel(endpoints)
                exact = IMSolverAnalytical().solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints=endpoints{iEndpoints}));
                numerical = IMSolverSpectral(nEVP=64).solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints=endpoints{iEndpoints}));
                testCase.verifyEqual(numerical.columnLabels,exact.columnLabels);
                testCase.verifyEqual(numerical.nColumns,numel(endpoints{iEndpoints}));
                testCase.verifyEqual(numerical.F(z),exact.F(z),AbsTol=2e-9);
                testCase.verifyEqual(numerical.G(z),exact.G(z),AbsTol=2e-9);
            end
        end
    end
end
