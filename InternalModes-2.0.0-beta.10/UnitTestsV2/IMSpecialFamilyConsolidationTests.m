classdef IMSpecialFamilyConsolidationTests < matlab.unittest.TestCase
    % Verify direct special-family construction after facade removal.
    methods (Test)
        function neutralBackgroundSelectionAndNormalizedDual(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],nModesPerInterval=4);
            selected = basis.selectColumns([4 2]);
            testCase.verifyEqual(selected.localModeNumber,[4 2]);
            selected.normalization = "maxAbsG";
            z = linspace(-1000,0,1001).';
            weights = ones(size(z)); weights([1 end]) = 0.5;
            projection = selected.projectionOnGrid(z,weights,variable="APV");
            testCase.verifyEqual(projection.physicalDualMatrix*selected.APV(z),eye(2),AbsTol=1e-12);
            testCase.verifyEqual(selected.gramMatrix(),0.5*eye(2),AbsTol=1e-15);
            custom = selected.addNormalization("double",@(~,~) 2);
            testCase.verifyEqual(custom.G(z,normalization="double"),0.5*selected.G(z,normalization="rmsG"));
            testCase.verifyError(@() IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],nModesPerInterval=4,zDomain=[-1000 0]),"IMNeutralLabelModesBasis:ConflictingBackground");
            nonneutral = IMBackground(N2=@(z) 1e-30*ones(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            testCase.verifyError(@() IMNeutralLabelModesBasis(nonneutral,neutralIntervals=[-1000 0],nModesPerInterval=0),"IMNeutralLabelModesBasis:InconsistentInterval");
        end

        function zeroDepthWaveHasExactIdentitiesAndCounts(testCase)
            rotationRate = 7.2921e-5; latitude = 33;
            f0 = 2*rotationRate*sind(latitude);
            background = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],rotationRate=rotationRate,latitude=latitude,surfaceCondition="freeSurface");
            basis = IMZeroDepthWaveModesBasis(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=3);
            z = linspace(-1000,0,1001).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=1e-14);
            testCase.verifyEqual(basis.Fz(z),1000*basis.Gzz(z),AbsTol=1e-16);
            testCase.verifyEqual(basis.h,[0 0 0]);
            testCase.verifyEqual(basis.frequency,repmat(abs(f0),1,3));
            testCase.verifyEqual(basis.gramMatrix(),eye(3),AbsTol=1e-14);
            testCase.verifyEqual(basis.G([-1000;0]),zeros(2,3));
            selected = basis.selectColumns([3 1]);
            testCase.verifyEqual(selected.G(z),basis.G(z)*[0 1;0 0;1 0]);
            empty = IMZeroDepthWaveModesBasis(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifyEqual(empty.nColumns,0);
            testCase.verifySize(empty.F(z),[1001 0]);
            testCase.verifySize(empty.gramMatrix(),[0 0]);
        end

        function selectedSpecialColumnsKeepCustomNormalizationAndDuals(testCase)
            callbacks = 0;
            neutral = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            labels = IMNeutralLabelModesBasis(neutral,neutralIntervals=[-1000 0],nModesPerInterval=4);
            labels = labels.addNormalization("identityDependent",@countedFactor);
            labels.normalization = "identityDependent";
            selected = labels.selectColumns([4 2]);
            testCase.verifyEqual(callbacks,0);
            z = linspace(-1000,0,1001).';
            original = labels.G(z);
            testCase.verifyEqual(selected.G(z),original(:,[4 2]),AbsTol=1e-15);
            weights = ones(size(z)); weights([1 end]) = 0.5;
            projection = selected.projectionOnGrid(z,weights,variable="APV");
            testCase.verifyEqual(projection.physicalDualMatrix*selected.APV(z),eye(2),AbsTol=1e-12);
            referenceGram = labels.gramMatrix();
            testCase.verifyEqual(selected.gramMatrix(),referenceGram([4 2],[4 2]));

            f0 = neutral.f0;
            critical = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=neutral.zDomain,surfaceCondition="freeSurface");
            waves = IMZeroDepthWaveModesBasis(critical,1e-4,criticalIntervals=[-1000 0],nModesPerInterval=4);
            waves = waves.addNormalization("identityDependent",@countedFactor);
            waves.normalization = "identityDependent";
            countBeforeSelection = callbacks;
            selectedWaves = waves.selectColumns([4 2]);
            testCase.verifyEqual(callbacks,countBeforeSelection);
            originalF = waves.F(z); originalG = waves.G(z);
            testCase.verifyEqual(selectedWaves.F(z),originalF(:,[4 2]));
            testCase.verifyEqual(selectedWaves.G(z),originalG(:,[4 2]));
            testCase.verifyEqual(selectedWaves.F(z),1000*selectedWaves.Gz(z),AbsTol=1e-14);

            function factor = countedFactor(~,column)
                callbacks = callbacks+1;
                factor = 2*column;
            end
        end

        function analyticalBackendReturnsExactCountedFamily(testCase)
            background = IMBackground.constantStratification(N0=1e-2,zDomain=[-1000 0],surfaceCondition="freeSurface");
            solver = IMSolverAnalytical();
            problem = IMInternalModes.waveModesAtWavenumber(background,k=1e-3);
            basis = solver.solveEVP(problem,nModes=4);
            testCase.verifyClass(basis,"IMAnalyticalInternalModesBasis");
            testCase.verifyEqual(basis.nColumns,4);
            internal = solver.solveEVP(problem,nModes=2,branch="internal");
            testCase.verifyEqual(internal.modeInfo.branch,repmat("internal",2,1));
            empty = solver.solveEVP(problem,nModes=0);
            testCase.verifyEqual(empty.nColumns,0);
            testCase.verifyError(@() solver.solveEVP(problem,nModes=1,branch="subinertial"),"IMSolver:InsufficientModes");
        end

        function endpointBackgroundAndExactDerivatives(testCase)
            background = IMBackground.constantStratification(N0=1e-2,zDomain=[-1000 0],surfaceCondition="freeSurface");
            problem = IMGeostrophicZeroAPVModes.atWavenumber(background,[1e-4 2e-4],endpoints="surface");
            solver = IMSolverAnalytical();
            basis = solver.solveGeostrophicZeroAPVModes(problem);
            z = linspace(-1000,0,21).';
            F = basis.F(z,pages=[2 1 2]);
            G = basis.G(z,pages=[2 1 2]);
            testCase.verifyEqual(basis.evaluate(z,variable="F",derivativeOrder=1,pages=[2 1 2]),-background.N2(z)/background.g.*G,AbsTol=1e-14);
            testCase.verifyEqual(basis.Gz(z,pages=[2 1 2]),reshape(-background.g*problem.k([2 1 2]).^2/background.f0^2,1,1,[]).*F,AbsTol=1e-13);
            testCase.verifyEqual(basis.nColumns,1);
            testCase.verifyEqual(basis.columnLabels,"surface");
            testCase.verifyError(@() IMGeostrophicZeroAPVModes.atWavenumber(background,1e-4,g=9.81),"IMGeostrophicZeroAPVModes:ConflictingBackground");
        end
    end
end
