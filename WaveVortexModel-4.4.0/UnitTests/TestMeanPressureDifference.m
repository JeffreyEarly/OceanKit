classdef TestMeanPressureDifference < matlab.unittest.TestCase
    properties
        transformType
        wvt
    end

    properties (ClassSetupParameter)
        transform = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG", ...
            barotropicQG="barotropicQG")
    end

    methods (TestClassSetup)
        function storeTransformType(testCase,transform)
            testCase.transformType = transform;
        end
    end

    methods (TestMethodSetup)
        function createFreshTransform(testCase)
            testCase.wvt = testCase.newTransform(testCase.transformType);
        end
    end

    methods (Test, TestTags = "full")
        function zeroStateReturnsScalarLogicalFalse(testCase)
            flag = testCase.wvt.hasMeanPressureDifference();
            testCase.verifyFalse(flag)
            testCase.verifyClass(flag,"logical")
            testCase.verifySize(flag,[1 1])
        end

        function nonMDAComponentsDoNotContribute(testCase)
            seedRandomNumberGenerator(testCase,2801);
            wvt = testCase.wvt;

            wvt.initWithRandomFlow('geostrophic',uvMax=0.01);
            testCase.verifyFalse(wvt.hasMeanPressureDifference())

            if wvt.hasWaveComponent
                wvt.initWithRandomFlow('wave',uvMax=0.01);
                testCase.verifyFalse(wvt.hasMeanPressureDifference())

                wvt.initWithRandomFlow('inertial',uvMax=0.01);
                testCase.verifyFalse(wvt.hasMeanPressureDifference())
            end
        end

        function internalMDAProducesPressureDifference(testCase)
            wvt = testCase.wvt;
            if ~isa(wvt,"WVMeanDensityAnomalyMethods")
                testCase.verifyFalse(wvt.hasMeanPressureDifference())
                return
            end

            mdaIndex = find(wvt.mdaComponent.maskA0,1);
            wvt.A0(mdaIndex) = 1;
            testCase.verifyTrue(wvt.hasMeanPressureDifference())
        end

        function nonMDASuperpositionDoesNotChangeResult(testCase)
            wvt = testCase.wvt;
            if ~isa(wvt,"WVMeanDensityAnomalyMethods")
                testCase.verifyFalse(wvt.hasMeanPressureDifference())
                return
            end

            mdaIndex = find(wvt.mdaComponent.maskA0,1);
            wvt.A0(mdaIndex) = 1;
            testCase.verifyTrue(wvt.hasMeanPressureDifference())

            seedRandomNumberGenerator(testCase,2802);
            wvt.addRandomFlow('geostrophic','wave','inertial',uvMax=0.01);
            testCase.verifyTrue(wvt.hasMeanPressureDifference())
        end

        function totalPressureGaugeOffsetDoesNotChangeResult(testCase)
            wvt = testCase.wvt;
            if ~isa(wvt,"WVMeanDensityAnomalyMethods")
                testCase.verifyFalse(wvt.hasMeanPressureDifference())
                return
            end

            mdaIndex = find(wvt.mdaComponent.maskA0,1);
            wvt.A0(mdaIndex) = 1;
            testCase.verifyTrue(wvt.hasMeanPressureDifference())

            totalPressure = wvt.variableWithName('p');
            gaugeOffset = 1e6*max(abs(totalPressure),[],"all");
            pressureAnnotation = WVVariableAnnotation('p',wvt.spatialDimensionNames,'kg/m/s2','test pressure with a common gauge offset');
            pressureOperation = WVOperation('p',pressureAnnotation,@(~) totalPressure+gaugeOffset);
            wvt.addOperation(pressureOperation,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            testCase.verifyEqual(wvt.variableWithName('p'),totalPressure+gaugeOffset)
            testCase.verifyTrue(wvt.hasMeanPressureDifference())
        end

        function relativeToleranceUsesAbsoluteDifference(testCase)
            if testCase.transformType ~= "constantHydrostatic"
                return
            end

            wvt = testCase.wvt;
            oddIndex = wvt.indexFromModeNumber(0,0,1);
            evenIndex = wvt.indexFromModeNumber(0,0,2);

            wvt.A0(evenIndex) = 1;
            evenPressure = wvt.variableWithName('p');
            pressureScale = max(abs(evenPressure),[],"all");
            testCase.verifyFalse(wvt.hasMeanPressureDifference())

            wvt.A0(:) = 0;
            wvt.A0(oddIndex) = 1;
            oddPressure = wvt.variableWithName('p');
            oddDifference = abs(mean(oddPressure(:,:,end),"all") - mean(oddPressure(:,:,1),"all"));

            belowToleranceAmplitude = 1e-7*pressureScale/oddDifference;
            wvt.A0(:) = 0;
            wvt.A0(evenIndex) = 1;
            wvt.A0(oddIndex) = belowToleranceAmplitude;
            testCase.verifyFalse(wvt.hasMeanPressureDifference())

            aboveToleranceAmplitude = 1e-3*pressureScale/oddDifference;
            wvt.A0(oddIndex) = aboveToleranceAmplitude;
            testCase.verifyTrue(wvt.hasMeanPressureDifference())

            wvt.A0(oddIndex) = -aboveToleranceAmplitude;
            testCase.verifyTrue(wvt.hasMeanPressureDifference())
        end
    end

    methods (Static, Access = private)
        function wvt = newTransform(transformType)
            Lxyz = [8e3 6e3 1e3];
            Nxyz = [8 6 9];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transformType
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "stratifiedQG"
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "barotropicQG"
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=false);
            end
        end
    end
end
