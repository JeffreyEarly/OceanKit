classdef TestEtaTrueOperation < matlab.unittest.TestCase

    methods (Test)
        function testConstantStratificationRegistersRhoNm(testCase)
            names = WVTransformConstantStratification.namesOfTransformVariables();
            requiredPropertyNames = WVTransformConstantStratification.newRequiredPropertyNames();
            propertyAnnotations = WVTransformConstantStratification.propertyAnnotationsForTransform();
            propertyNames = string({propertyAnnotations.name});

            testCase.verifyTrue(ismember('rho_nm',names));
            testCase.verifyTrue(ismember('shouldUseTrueNoMotionProfile',requiredPropertyNames));
            testCase.verifyTrue(ismember("shouldUseTrueNoMotionProfile",propertyNames));
        end

        function testEtaTrueUsesRhoNm0ByDefault(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.addOperation(TestEtaTrueOperation.failingRhoNmOperation(),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));

            eta_true = wvt.eta_true;
            expectedEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,wvt.rho_nm0);
            testCase.verifyEqual(eta_true,expectedEtaTrue,AbsTol=1e-12);
        end

        function testEtaTrueUsesRegisteredRhoNmWhenRequested(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform(shouldUseTrueNoMotionProfile=true);
            rho_nm = TestEtaTrueOperation.testRhoNmProfile(wvt);
            wvt.addOperation(TestEtaTrueOperation.fixedRhoNmOperation(rho_nm),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));

            eta_true = wvt.eta_true;
            expectedEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,rho_nm);
            defaultEtaTrue = TestEtaTrueOperation.etaTrueForProfile(wvt,wvt.rho_nm0);

            testCase.verifyEqual(eta_true,expectedEtaTrue,AbsTol=1e-12);
            testCase.verifyGreaterThan(max(abs(eta_true(:) - defaultEtaTrue(:))),1e-8);
        end

        function testEtaTrueDoesNotWarnWhenUsingRhoNm0(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform();
            wvt.addOperation(EtaTrueOperationToolboxUnavailable(wvt),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));
        end

        function testEtaTrueWarnsOnceWhenOptimizationToolboxUnavailable(testCase)
            wvt = TestEtaTrueOperation.hydrostaticTransform(shouldUseTrueNoMotionProfile=true);
            rho_nm = TestEtaTrueOperation.testRhoNmProfile(wvt);
            wvt.addOperation(TestEtaTrueOperation.fixedRhoNmOperation(rho_nm),shouldOverwriteExisting=true,shouldSuppressWarning=true);
            wvt.addOperation(EtaTrueOperationToolboxUnavailable(wvt),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyWarning(@() wvt.performOperationWithName('eta_true'),'EtaTrueOperation:OptimizationToolboxUnavailable');
            wvt.clearVariableCacheOfApAmA0DependentVariables();
            testCase.verifyWarningFree(@() wvt.performOperationWithName('eta_true'));
        end

        function testShouldUseTrueNoMotionProfilePersistsThroughTransformCopies(testCase)
            wvtHydrostatic = TestEtaTrueOperation.hydrostaticTransform(shouldUseTrueNoMotionProfile=true,shouldAntialias=true);
            testCase.verifyTrue(wvtHydrostatic.waveVortexTransformWithResolution([12 12 7]).shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtHydrostatic.waveVortexTransformWithExplicitAntialiasing().shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtHydrostatic.boussinesqTransform().shouldUseTrueNoMotionProfile);

            wvtBoussinesq = TestEtaTrueOperation.boussinesqTransform(shouldUseTrueNoMotionProfile=true,shouldAntialias=true);
            testCase.verifyTrue(wvtBoussinesq.waveVortexTransformWithResolution([12 12 7]).shouldUseTrueNoMotionProfile);
            testCase.verifyTrue(wvtBoussinesq.waveVortexTransformWithExplicitAntialiasing().shouldUseTrueNoMotionProfile);
        end

        function testShouldUseTrueNoMotionProfilePersistsThroughRoundTrip(testCase)
            transforms = {
                TestEtaTrueOperation.hydrostaticTransform(shouldUseTrueNoMotionProfile=true)
                TestEtaTrueOperation.boussinesqTransform(shouldUseTrueNoMotionProfile=true)
                };

            for iTransform = 1:numel(transforms)
                path = [tempname,'.nc'];
                cleanup = onCleanup(@() TestEtaTrueOperation.deleteIfExists(path));
                transforms{iTransform}.writeToFile(path,shouldOverwriteExisting=true);
                [wvt2,ncfile] = WVTransform.waveVortexTransformFromFile(path);

                testCase.verifyTrue(wvt2.shouldUseTrueNoMotionProfile);

                ncfile = [];
                wvt2 = [];
                clear cleanup
            end
        end
    end

    methods (Static, Access=private)
        function wvt = constantTransform(options)
            arguments
                options.shouldUseTrueNoMotionProfile (1,1) logical = false
                options.shouldAntialias (1,1) logical = false
            end

            wvt = WVTransformConstantStratification([4e3, 4e3, 2e3], [8 8 5], ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias, ...
                shouldUseTrueNoMotionProfile=options.shouldUseTrueNoMotionProfile);
        end

        function wvt = hydrostaticTransform(options)
            arguments
                options.shouldUseTrueNoMotionProfile (1,1) logical = false
                options.shouldAntialias (1,1) logical = false
            end

            N2 = @(z) (5.2e-3)^2 * ones(size(z));
            wvt = WVTransformHydrostatic([4e3, 4e3, 2e3], [8 8 5], ...
                N2Function=N2, ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias, ...
                shouldUseTrueNoMotionProfile=options.shouldUseTrueNoMotionProfile);
        end

        function wvt = boussinesqTransform(options)
            arguments
                options.shouldUseTrueNoMotionProfile (1,1) logical = false
                options.shouldAntialias (1,1) logical = false
            end

            N2 = @(z) (5.2e-3)^2 * ones(size(z));
            wvt = WVTransformBoussinesq([4e3, 4e3, 2e3], [8 8 5], ...
                N2Function=N2, ...
                latitude=30, ...
                shouldAntialias=options.shouldAntialias, ...
                shouldUseTrueNoMotionProfile=options.shouldUseTrueNoMotionProfile);
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
            data = wvt.rho0 - rho_nm;
            spline_nm = BSpline(K,BSpline.knotPointsForDataPoints(wvt.z,K=K));
            spline_nm.x_mean = mean(data);
            spline_nm.x_std = std(data);
            Z = BSpline.Spline(wvt.z, spline_nm.t_knot, spline_nm.K);
            spline_nm.m = Z\((data - spline_nm.x_mean)/spline_nm.x_std);

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
