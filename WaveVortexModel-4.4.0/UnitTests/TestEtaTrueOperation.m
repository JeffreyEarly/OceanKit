classdef TestEtaTrueOperation < matlab.unittest.TestCase

    methods (Test, TestTags = "full")
        function testTransformsRegisterRhoNmWithoutPersistingFlag(testCase)
            testCase.verifyTrue(ismember('rho_nm',WVTransformConstantStratification.namesOfTransformVariables()));
            testCase.verifyTrue(ismember('rho_nm',WVTransformHydrostatic.namesOfTransformVariables()));
            testCase.verifyTrue(ismember('rho_nm',WVTransformBoussinesq.namesOfTransformVariables()));

            testCase.verifyFalse(ismember('shouldUseTrueNoMotionProfile', ...
                WVTransformConstantStratification.newRequiredPropertyNames()));
            testCase.verifyFalse(ismember('shouldUseTrueNoMotionProfile', ...
                WVTransformHydrostatic.newRequiredPropertyNames()));
            testCase.verifyFalse(ismember('shouldUseTrueNoMotionProfile', ...
                WVTransformBoussinesq.newRequiredPropertyNames()));

            testCase.verifyFalse(TestEtaTrueOperation.hasPropertyAnnotation( ...
                WVTransformConstantStratification.propertyAnnotationsForTransform(),"shouldUseTrueNoMotionProfile"));
            testCase.verifyFalse(TestEtaTrueOperation.hasPropertyAnnotation( ...
                WVTransformHydrostatic.propertyAnnotationsForTransform(),"shouldUseTrueNoMotionProfile"));
            testCase.verifyFalse(TestEtaTrueOperation.hasPropertyAnnotation( ...
                WVTransformBoussinesq.propertyAnnotationsForTransform(),"shouldUseTrueNoMotionProfile"));
        end

        function testShouldUseTrueNoMotionProfileIsPostInitializationOnly(testCase)
            testCase.verifyFalse(TestEtaTrueOperation.canConstructConstantWithShouldUseTrueNoMotionProfile());
            testCase.verifyFalse(TestEtaTrueOperation.canConstructHydrostaticWithShouldUseTrueNoMotionProfile());
            testCase.verifyFalse(TestEtaTrueOperation.canConstructBoussinesqWithShouldUseTrueNoMotionProfile());

            transforms = TestEtaTrueOperation.allTransforms();
            for index = 1:numel(transforms)
                wvt = transforms{index};
                testCase.verifyTrue(wvt.shouldUseTrueNoMotionProfile);
                wvt.shouldUseTrueNoMotionProfile = false;
                testCase.verifyFalse(wvt.shouldUseTrueNoMotionProfile);
                wvt.shouldUseTrueNoMotionProfile = true;
                testCase.verifyTrue(wvt.shouldUseTrueNoMotionProfile);
            end
        end

        function testExplicitFalseUsesOnlyTheReferenceProfile(testCase)
            transforms = TestEtaTrueOperation.allTransforms();
            for index = 1:numel(transforms)
                wvt = transforms{index};
                wvt.shouldUseTrueNoMotionProfile = false;
                wvt.addOperation(TestEtaTrueOperation.failingRhoNmOperation(),shouldOverwriteExisting=true,shouldSuppressWarning=true);
                expectedEta = .02*wvt.Lz*sin(2*pi*wvt.X/wvt.Lx).*sin(pi*(wvt.Z+wvt.Lz)/wvt.Lz);
                materialHeight = wvt.Z-expectedEta;
                slope = (wvt.rho_nm0(1)-wvt.rho_nm0(end))/wvt.Lz;
                wvt.addToVariableCache('rho_total',wvt.rho_nm0(end)-slope*materialHeight);
                testCase.verifyWarningFree(@()wvt.performOperationWithName('ape'));
                testCase.verifyEqual(wvt.eta_true,expectedEta,AbsTol=2e-7);
                testCase.verifyEqual(wvt.ape,(wvt.g*slope/wvt.rho0)*expectedEta.^2/2,AbsTol=1e-8);
                testCase.verifyFalse(isKey(wvt.variableCache,'rho_nm'));
            end
        end

        function testDefaultUsesTheActualChangedNonlinearRestProfile(testCase)
            transforms = TestEtaTrueOperation.allTransforms();
            for index = 1:numel(transforms)
                wvt = transforms{index};
                testCase.verifyTrue(wvt.shouldUseTrueNoMotionProfile);
                span = wvt.rho_nm0(1)-wvt.rho_nm0(end);
                q = (wvt.z+wvt.Lz)/wvt.Lz;
                changed = wvt.rho_nm0+.1*span*q.*(1-q);
                density = repmat(reshape(changed,1,1,[]),wvt.Nx,wvt.Ny,1);
                % Prescribe the physical density input, retaining the real
                % rho_nm, eta_true, APE, and APV registered operations.
                wvt.addToVariableCache('rho_total',density);
                testCase.verifyWarningFree(@()wvt.performOperationWithName('ape'));
                testCase.verifyEqual(wvt.rho_nm,changed,AbsTol=3e-13);
                testCase.verifyEqual(wvt.eta_true,zeros(wvt.spatialMatrixSize),AbsTol=2e-8);
                testCase.verifyEqual(wvt.ape,zeros(wvt.spatialMatrixSize),AbsTol=1e-18);
                testCase.verifyEqual(wvt.apv,zeros(wvt.spatialMatrixSize),AbsTol=1e-13);

                % Explicit false selects the original linear profile for
                % BOTH displacement and energy, even after true caches exist.
                wvt.shouldUseTrueNoMotionProfile = false;
                expectedEta = repmat(reshape(.1*wvt.Lz*q.*(1-q),1,1,[]),wvt.Nx,wvt.Ny,1);
                N2 = wvt.g*span/(wvt.rho0*wvt.Lz);
                testCase.verifyEqual(wvt.eta_true,expectedEta,AbsTol=2e-7);
                testCase.verifyEqual(wvt.ape,N2*expectedEta.^2/2,AbsTol=1e-8);
                testCase.verifyGreaterThan(max(abs(wvt.eta_true),[],"all"),10);
                wvt.shouldUseTrueNoMotionProfile = true;
                testCase.verifyEqual(wvt.apv,zeros(wvt.spatialMatrixSize),AbsTol=1e-13);
                testCase.verifyEqual(wvt.ape,zeros(wvt.spatialMatrixSize),AbsTol=1e-18);
            end
        end

        function testActualDefaultRecoversAnEqualVolumeParcelRearrangement(testCase)
            wvt = TestEtaTrueOperation.constantTransform();
            operation = WVNoMotionProfileOperation();
            testCase.verifyEqual(operation.solver,"dampedLeastSquares");
            wvt.addOperation(operation,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            density = repmat(reshape(wvt.rho_nm0,1,1,[]),wvt.Nx,wvt.Ny,1);
            testCase.assertEqual(wvt.z_int(2),wvt.z_int(4));
            density(1,1,[2 4]) = density(1,1,[4 2]);
            materialHeight = wvt.Z;
            materialHeight(1,1,[2 4]) = materialHeight(1,1,[4 2]);
            wvt.addToVariableCache('rho_total',density);
            expectedEta = wvt.Z-materialHeight;
            testCase.verifyWarningFree(@()wvt.performOperationWithName('eta_true'));
            testCase.verifyEqual(wvt.rho_nm,wvt.rho_nm0,AbsTol=2e-12);
            testCase.verifyEqual(wvt.eta_true,expectedEta,AbsTol=2e-7);
            testCase.verifyEqual(wvt.ape,wvt.N2(1)*expectedEta.^2/2,AbsTol=1e-7);
            testCase.verifyEqual(operation.lastSolverOutput.exitflag,1);
            testCase.verifyLessThan(operation.lastSolverOutput.maximumResidual,1e-12);
        end

        function testRejectedProfileDoesNotCacheDependentDiagnostics(testCase)
            wvt = TestEtaTrueOperation.constantTransform();
            rest = repmat(reshape(wvt.rho_nm0,1,1,[]),wvt.Nx,wvt.Ny,1);
            wvt.addToVariableCache('rho_total',wvt.rho0*ones(size(rest)));
            testCase.verifyError(@()wvt.performOperationWithName('ape'),'WVNoMotionProfileOperation:NonInvertibleDistribution');
            for name = ["rho_nm","eta_true","ape"]
                testCase.verifyFalse(isKey(wvt.variableCache,name));
            end
            wvt.addToVariableCache('rho_total',rest);
            testCase.verifyWarningFree(@()wvt.performOperationWithName('ape'));
            testCase.verifyEqual(wvt.eta_true,zeros(wvt.spatialMatrixSize),AbsTol=2e-8);
            testCase.verifyEqual(wvt.ape,zeros(wvt.spatialMatrixSize),AbsTol=1e-18);
        end

        function testActualDefaultAcceptsChangedRestExtrema(testCase)
            wvt = TestEtaTrueOperation.constantTransform();
            changed = wvt.rho0+2*(wvt.rho_nm0-wvt.rho0)+.5;
            wvt.addToVariableCache('rho_total',repmat(reshape(changed,1,1,[]),wvt.Nx,wvt.Ny,1));
            testCase.verifyWarningFree(@()wvt.performOperationWithName('ape'));
            testCase.verifyEqual(wvt.rho_nm,changed);
            testCase.verifyEqual(wvt.eta_true,zeros(wvt.spatialMatrixSize),AbsTol=2e-8);
            testCase.verifyEqual(wvt.ape,zeros(wvt.spatialMatrixSize),AbsTol=1e-18);
        end

        function testCoefficientDrivenDensityDiagnosticsSurviveOutputAndRestart(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            wvt = TestEtaTrueOperation.constantTransform();
            initialEta = .02*wvt.Lz*sin(pi*(wvt.Z+wvt.Lz)/wvt.Lz);
            expectedDensity = reshape(wvt.rho_nm0,1,1,[])+(wvt.rho0/wvt.g)*wvt.N2(1)*initialEta;
            wvt.initWithUVEta(zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),initialEta);
            testCase.verifyGreaterThan(max(abs(wvt.A0),[],"all"),0);
            testCase.verifyEqual(wvt.rho_total,expectedDensity,AbsTol=5e-13);
            expectedProfile = reshape(expectedDensity(1,1,:),[],1);
            testCase.verifyLessThan(max(diff(expectedProfile)),0);
            testCase.verifyGreaterThan(max(abs(expectedProfile-wvt.rho_nm0)),.01);
            testCase.verifyEqual(wvt.rho_nm,expectedProfile,AbsTol=5e-13);
            testCase.verifyEqual(wvt.eta_true,zeros(wvt.spatialMatrixSize),AbsTol=2e-8);
            testCase.verifyEqual(wvt.ape,zeros(wvt.spatialMatrixSize),AbsTol=1e-18);
            testCase.verifyEqual(wvt.apv,zeros(wvt.spatialMatrixSize),AbsTol=1e-13);

            model = WVModel(wvt,shouldUseLinearDynamics=true);
            modelCleanup = onCleanup(@()model.closeNetCDFFile());
            names = {'rho_nm','eta_true','ape','apv'};
            model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
            path = fullfile(fixture.Folder,'coefficient-density.nc');
            model.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true);
            model.setupIntegrator(integratorType="fixed",deltaT=.25);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
            restored = WVModel.modelFromFile(path);
            restoredCleanup = onCleanup(@()restored.closeNetCDFFile());
            testCase.verifyTrue(restored.wvt.shouldUseTrueNoMotionProfile);
            testCase.verifyEqual(restored.wvt.A0,wvt.A0);
            testCase.verifyEqual(restored.wvt.rho_total,expectedDensity,AbsTol=5e-13);
            testCase.verifyEqual(restored.wvt.rho_nm,expectedProfile,AbsTol=5e-13);
            restored.setupIntegrator(integratorType="fixed",deltaT=.25);
            restored.integrateToTime(1.5,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            restored.closeNetCDFFile();
            times = ncread(path,'/wave-vortex/t');
            testCase.verifyEqual(times,[0;.5;1;1.5]);
            expected = {expectedProfile,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize)};
            tolerances = [5e-13 2e-8 1e-18 1e-13];
            for index = 1:numel(names)
                actual = restored.wvt.variableWithName(names{index});
                testCase.verifyEqual(actual,expected{index},AbsTol=tolerances(index));
                saved = ncread(path,"/wave-vortex/"+names{index});
                testCase.verifyEqual(reshape(saved,[],numel(times)),repmat(expected{index}(:),1,numel(times)),AbsTol=tolerances(index));
            end
            clear restoredCleanup modelCleanup
        end

        function testShouldUseTrueNoMotionProfilePersistsThroughTransformCopies(testCase)
            wvtConstant = TestEtaTrueOperation.constantTransform(shouldAntialias=true);
            wvtConstant.shouldUseTrueNoMotionProfile = false;
            wvtConstantX2 = wvtConstant.waveVortexTransformWithResolution([12 12 7]);
            wvtConstantAntialias = wvtConstant.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyFalse(wvtConstantX2.shouldUseTrueNoMotionProfile);
            testCase.verifyFalse(wvtConstantAntialias.shouldUseTrueNoMotionProfile);

            wvtHydrostatic = TestEtaTrueOperation.hydrostaticTransform(shouldAntialias=true);
            wvtHydrostatic.shouldUseTrueNoMotionProfile = false;
            wvtHydrostaticX2 = wvtHydrostatic.waveVortexTransformWithResolution([12 12 7]);
            wvtHydrostaticAntialias = wvtHydrostatic.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyFalse(wvtHydrostaticX2.shouldUseTrueNoMotionProfile);
            testCase.verifyFalse(wvtHydrostaticAntialias.shouldUseTrueNoMotionProfile);
            testCase.verifyFalse(wvtHydrostatic.boussinesqTransform().shouldUseTrueNoMotionProfile);

            wvtBoussinesq = TestEtaTrueOperation.boussinesqTransform(shouldAntialias=true);
            wvtBoussinesq.shouldUseTrueNoMotionProfile = false;
            wvtBoussinesqX2 = wvtBoussinesq.waveVortexTransformWithResolution([12 12 7]);
            wvtBoussinesqAntialias = wvtBoussinesq.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyFalse(wvtBoussinesqX2.shouldUseTrueNoMotionProfile);
            testCase.verifyFalse(wvtBoussinesqAntialias.shouldUseTrueNoMotionProfile);
        end

        function testShouldUseTrueNoMotionProfileDoesNotPersistThroughRoundTrip(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            transforms = TestEtaTrueOperation.allTransforms();

            for iTransform = 1:numel(transforms)
                transforms{iTransform}.shouldUseTrueNoMotionProfile = false;
                path = fullfile(fixture.Folder,sprintf('eta-true-%d.nc',iTransform));
                writtenFile = transforms{iTransform}.writeToFile(path,shouldOverwriteExisting=true);
                writtenFile.close();
                [wvt2,ncfile] = WVTransform.waveVortexTransformFromFile(path);
                cleanup = onCleanup(@()TestEtaTrueOperation.closeIfOpen(ncfile));

                testCase.verifyTrue(wvt2.shouldUseTrueNoMotionProfile);

                ncfile.close();
                clear cleanup
            end
        end

        function testShouldUseTrueNoMotionProfileInvalidatesDensityDiagnosticCaches(testCase)
            transforms = {
                TestEtaTrueOperation.constantTransform()
                TestEtaTrueOperation.hydrostaticTransform()
                TestEtaTrueOperation.boussinesqTransform()
                };

            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                wvt.addToVariableCache("rho_nm",wvt.rho_nm0);
                wvt.addToVariableCache("eta",zeros(wvt.spatialMatrixSize));
                for name = ["eta_true","ape","apv"]
                    wvt.addToVariableCache(name,ones(wvt.spatialMatrixSize));
                end

                wvt.shouldUseTrueNoMotionProfile = true;
                testCase.verifyTrue(isKey(wvt.variableCache,"rho_nm"));
                testCase.verifyTrue(isKey(wvt.variableCache,"eta"));

                wvt.shouldUseTrueNoMotionProfile = false;
                testCase.verifyFalse(isKey(wvt.variableCache,"rho_nm"));
                for name = ["eta_true","ape","apv"]
                    testCase.verifyFalse(isKey(wvt.variableCache,name));
                end
                testCase.verifyTrue(isKey(wvt.variableCache,"eta"));

                wvt.addToVariableCache("rho_nm",wvt.rho_nm0);
                wvt.shouldUseTrueNoMotionProfile = false;
                testCase.verifyTrue(isKey(wvt.variableCache,"rho_nm"));
                for name = ["eta_true","ape","apv"]
                    wvt.addToVariableCache(name,ones(wvt.spatialMatrixSize));
                end
                wvt.shouldUseTrueNoMotionProfile = true;
                for name = ["rho_nm","eta_true","ape","apv"]
                    testCase.verifyFalse(isKey(wvt.variableCache,name));
                end
                testCase.verifyTrue(isKey(wvt.variableCache,"eta"));
            end
        end
    end

    methods (Static, Access=private)
        function transforms = allTransforms()
            transforms = {TestEtaTrueOperation.constantTransform(), ...
                TestEtaTrueOperation.hydrostaticTransform(),TestEtaTrueOperation.boussinesqTransform()};
        end

        function wvt = constantTransform(options)
            arguments
                options.shouldAntialias (1,1) logical = false
            end

            wvt = WVTransformConstantStratification([4e3, 4e3, 2e3], [8 8 5], ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias);
        end

        function wvt = hydrostaticTransform(options)
            arguments
                options.shouldAntialias (1,1) logical = false
            end

            N2 = @(z) (5.2e-3)^2 * ones(size(z));
            wvt = WVTransformHydrostatic([4e3, 4e3, 2e3], [8 8 5], ...
                N2Function=N2, ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias);
        end

        function wvt = boussinesqTransform(options)
            arguments
                options.shouldAntialias (1,1) logical = false
            end

            N2 = @(z) (5.2e-3)^2 * ones(size(z));
            wvt = WVTransformBoussinesq([4e3, 4e3, 2e3], [8 8 5], ...
                N2Function=N2, ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias);
        end

        function tf = hasPropertyAnnotation(propertyAnnotations,name)
            propertyNames = string({propertyAnnotations.name});
            tf = ismember(name,propertyNames);
        end

        function tf = canConstructConstantWithShouldUseTrueNoMotionProfile()
            try
                WVTransformConstantStratification([4e3, 4e3, 2e3], [8 8 5], ...
                    latitude=30, shouldAntialias=false, shouldUseTrueNoMotionProfile=true);
                tf = true;
            catch
                tf = false;
            end
        end

        function tf = canConstructHydrostaticWithShouldUseTrueNoMotionProfile()
            try
                N2 = @(z) (5.2e-3)^2 * ones(size(z));
                WVTransformHydrostatic([4e3, 4e3, 2e3], [8 8 5], ...
                    N2Function=N2, latitude=30, shouldAntialias=false, shouldUseTrueNoMotionProfile=true);
                tf = true;
            catch
                tf = false;
            end
        end

        function tf = canConstructBoussinesqWithShouldUseTrueNoMotionProfile()
            try
                N2 = @(z) (5.2e-3)^2 * ones(size(z));
                WVTransformBoussinesq([4e3, 4e3, 2e3], [8 8 5], ...
                    N2Function=N2, latitude=30, shouldAntialias=false, shouldUseTrueNoMotionProfile=true);
                tf = true;
            catch
                tf = false;
            end
        end

        function op = failingRhoNmOperation()
            outputVariables(1) = WVVariableAnnotation('rho_nm',{'z'},'kg m^{-3}', 'test no-motion density profile');
            op = WVOperation('rho_nm',outputVariables,@(~) error('TestEtaTrueOperation:RhoNmShouldNotBeComputed', ...
                'rho_nm should not be computed when shouldUseTrueNoMotionProfile is false.'));
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
