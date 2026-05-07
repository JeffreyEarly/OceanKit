classdef TestEtaTrueOperation < matlab.unittest.TestCase

    methods (Test)
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

            wvt = TestEtaTrueOperation.hydrostaticTransform();
            testCase.verifyFalse(wvt.shouldUseTrueNoMotionProfile);
            wvt.shouldUseTrueNoMotionProfile = true;
            testCase.verifyTrue(wvt.shouldUseTrueNoMotionProfile);
            wvt.shouldUseTrueNoMotionProfile = false;
            testCase.verifyFalse(wvt.shouldUseTrueNoMotionProfile);
        end

        function testEtaTrueUsesRhoNm0ByDefault(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.addOperation(TestEtaTrueOperation.failingRhoNmOperation(), ...
                shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));

            eta_true = wvt.eta_true;
            expectedEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,wvt.rho_nm0);
            testCase.verifyEqual(eta_true,expectedEtaTrue,AbsTol=1e-12);
        end

        function testEtaTrueUsesRegisteredRhoNmWhenRequested(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.shouldUseTrueNoMotionProfile = true;
            rho_nm = TestEtaTrueOperation.testRhoNmProfile(wvt);
            wvt.addOperation(TestEtaTrueOperation.fixedRhoNmOperation(rho_nm), ...
                shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));

            eta_true = wvt.eta_true;
            expectedEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,rho_nm);
            defaultEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,wvt.rho_nm0);

            testCase.verifyEqual(eta_true,expectedEtaTrue,AbsTol=1e-12);
            testCase.verifyGreaterThan(max(abs(eta_true(:) - defaultEtaTrue(:))),1e-8);
        end

        function testEtaTrueDoesNotWarnWhenUsingRhoNm0(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.addOperation(EtaTrueOperationToolboxUnavailable(wvt), ...
                shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));
        end

        function testEtaTrueWarnsOnceWhenOptimizationToolboxUnavailable(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.shouldUseTrueNoMotionProfile = true;
            rho_nm = TestEtaTrueOperation.testRhoNmProfile(wvt);
            wvt.addOperation(TestEtaTrueOperation.fixedRhoNmOperation(rho_nm), ...
                shouldOverwriteExisting=true,shouldSuppressWarning=true);
            wvt.addOperation(EtaTrueOperationToolboxUnavailable(wvt), ...
                shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarning(@() wvt.performOperationWithName('eta_true'), ...
                'EtaTrueOperation:OptimizationToolboxUnavailable');
            wvt.clearVariableCacheOfApAmA0DependentVariables();
            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));
        end

        function testShouldUseTrueNoMotionProfilePersistsThroughTransformCopies(testCase)
            wvtHydrostatic = TestEtaTrueOperation.hydrostaticTransform(shouldAntialias=true);
            wvtHydrostatic.shouldUseTrueNoMotionProfile = true;
            wvtHydrostaticX2 = wvtHydrostatic.waveVortexTransformWithResolution([12 12 7]);
            wvtHydrostaticAntialias = wvtHydrostatic.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyTrue(wvtHydrostaticX2.shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtHydrostaticAntialias.shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtHydrostatic.boussinesqTransform().shouldUseTrueNoMotionProfile);

            wvtBoussinesq = TestEtaTrueOperation.boussinesqTransform(shouldAntialias=true);
            wvtBoussinesq.shouldUseTrueNoMotionProfile = true;
            wvtBoussinesqX2 = wvtBoussinesq.waveVortexTransformWithResolution([12 12 7]);
            wvtBoussinesqAntialias = wvtBoussinesq.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyTrue(wvtBoussinesqX2.shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtBoussinesqAntialias.shouldUseTrueNoMotionProfile);
        end

        function testShouldUseTrueNoMotionProfileDoesNotPersistThroughRoundTrip(testCase)
            transforms = {
                TestEtaTrueOperation.hydrostaticTransform()
                TestEtaTrueOperation.boussinesqTransform()
                };

            for iTransform = 1:numel(transforms)
                transforms{iTransform}.shouldUseTrueNoMotionProfile = true;
                path = [tempname,'.nc'];
                cleanup = onCleanup(@() TestEtaTrueOperation.deleteIfExists(path));
                transforms{iTransform}.writeToFile(path,shouldOverwriteExisting=true);
                [wvt2,ncfile] = WVTransform.waveVortexTransformFromFile(path);

                testCase.verifyFalse(wvt2.shouldUseTrueNoMotionProfile);

                ncfile = [];
                wvt2 = [];
                clear cleanup
            end
        end

        function testShouldUseTrueNoMotionProfileInvalidatesOnlyRhoNmCache(testCase)
            transforms = {
                TestEtaTrueOperation.constantTransform()
                TestEtaTrueOperation.hydrostaticTransform()
                TestEtaTrueOperation.boussinesqTransform()
                };

            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                wvt.addToVariableCache("rho_nm",wvt.rho_nm0);
                wvt.addToVariableCache("eta",zeros(wvt.spatialMatrixSize));

                wvt.shouldUseTrueNoMotionProfile = false;
                testCase.verifyTrue(isKey(wvt.variableCache,"rho_nm"));
                testCase.verifyTrue(isKey(wvt.variableCache,"eta"));

                wvt.shouldUseTrueNoMotionProfile = true;
                testCase.verifyFalse(isKey(wvt.variableCache,"rho_nm"));
                testCase.verifyTrue(isKey(wvt.variableCache,"eta"));

                wvt.addToVariableCache("rho_nm",wvt.rho_nm0);
                wvt.shouldUseTrueNoMotionProfile = true;
                testCase.verifyTrue(isKey(wvt.variableCache,"rho_nm"));
            end
        end
    end

    methods (Static, Access=private)
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

        function op = fixedRhoNmOperation(rho_nm)
            outputVariables(1) = WVVariableAnnotation('rho_nm',{'z'},'kg m^{-3}', 'test no-motion density profile');
            op = WVOperation('rho_nm',outputVariables,@(~) rho_nm);
        end

        function op = failingRhoNmOperation()
            outputVariables(1) = WVVariableAnnotation('rho_nm',{'z'},'kg m^{-3}', 'test no-motion density profile');
            op = WVOperation('rho_nm',outputVariables,@(~) error('TestEtaTrueOperation:RhoNmShouldNotBeComputed', ...
                'rho_nm should not be computed when shouldUseTrueNoMotionProfile is false.'));
        end

        function rho_nm = testRhoNmProfile(wvt)
            delta = diff(wvt.rho_nm0);
            if isempty(delta)
                rho_nm = wvt.rho_nm0;
                return
            end

            weights = 1 + 0.2*cos(linspace(0,pi,numel(delta))');
            scaledDelta = delta .* weights;
            scaledDelta = scaledDelta * (sum(delta)/sum(scaledDelta));

            rho_nm = [wvt.rho_nm0(1); wvt.rho_nm0(1) + cumsum(scaledDelta)];
        end

        function eta_true = etaTrueForProfile(wvt,rho_nm)
            K = min(wvt.Nz,8);
            S = K - 1;
            data = wvt.rho0 - rho_nm;
            knotPoints = BSpline.knotPointsForDataPoints(wvt.z,S=S);
            Z = BSpline.matrixForDataPoints(wvt.z,knotPoints=knotPoints,S=S);
            xMean = mean(data);
            xStd = std(data);
            xi = Z\((data - xMean)/xStd);
            spline_nm = BSpline(S=S,knotPoints=knotPoints,xi=xi,xMean=xMean,xStd=xStd);

            rho_total = (wvt.rhoFunction(wvt.Z) - wvt.rho0) + wvt.rho_e;
            zMinusEta = EtaTrueOperation.fInverseBisection(spline_nm,-rho_total(:),-wvt.Lz,0,1e-12);
            zMinusEta = reshape(zMinusEta,size(wvt.X));
            eta_true = wvt.Z - zMinusEta;
        end

        function deleteIfExists(path)
            if isfile(path)
                try
                    delete(path);
                catch
                end
            end
        end
    end
end
