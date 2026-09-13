classdef TestInternalRuntimeArtifacts < matlab.unittest.TestCase
    methods (Test, TestTags="full")
        function correctedInternalNamesAreActive(testCase)
            testCase.verifyNotEmpty(meta.class.fromName('WVModelAdaptiveTimeStepMethods'))
            testCase.verifyNotEmpty(meta.class.fromName('WVModelAdaptiveTimeStepCellMethods'))
            testCase.verifyEmpty(meta.class.fromName('WVModelAdapativeTimeStepMethods'))
            testCase.verifyEmpty(meta.class.fromName('WVModelAdapativeTimeStepCellMethods'))

            modelMethods = string(methods('WVModel'));
            testCase.verifyTrue(any(modelMethods == "resetAdaptiveTimeStepIntegrator"))
            testCase.verifyFalse(any(modelMethods == "resetAdapativeTimeStepIntegrator"))

            wvt = WVTransformBarotropicQG([4000 3000],[8 6],latitude=45,shouldAntialias=false);
            testCase.verifyWarningFree(@()mustBeDoublyPeriodicFPlane(wvt))
            testCase.verifyError(@()mustBeDoublyPeriodicFPlane(1),'mustBeDoublyPeriodicFPlane:invalidClass')
            testCase.verifyEmpty(which('mustBeDoulbyPeriodicFPlane'))
        end

        function legacyArtifactsAreOffRuntimePath(testCase)
            testCase.verifyEmpty(which('WVOffGridTransform'))
            testCase.verifyEmpty(which('WVAdaptiveDiffusivity'))
            testCase.verifyEmpty(which('WVAdaptiveViscosity'))
            testCase.verifyEmpty(which('WVMeanFlowForcing'))
            testCase.verifyEmpty(which('WVSpectralVanishingViscosity'))
            stratificationMethods = string(methods('WVStratification'));
            transformMethods = string(methods('WVTransformConstantStratification'));
            testCase.verifyFalse(any(stratificationMethods == "verticalProjectionOperatorsWithFreeSurface"))
            testCase.verifyFalse(any(transformMethods == "speedTest"))
        end

        function retainedInternalTransformsRemainAvailable(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.verifyFalse(isfile(fullfile(repositoryRoot,"FastTransforms","@WVFastTransformDoublyPeriodicFFTW","WVFastTransformDoublyPeriodicFFTW.m")))
            testCase.verifyFalse(isfile(fullfile(repositoryRoot,"FastTransforms","fftw_dft2.m")))
            wvt = WVTransformBarotropicQG([4000 3000],[8 6],latitude=45,shouldAntialias=false);
            testCase.verifyTrue(ismethod(wvt,'transformToSpatialDomainFromDFTGridAtPosition'))
            testCase.verifyTrue(ismethod(wvt,'transformToSpatialDomainWithFourierAtPosition'))
            testCase.verifyError(@()wvt.variableAtPositionWithName(0,0,0,'u',interpolationMethod='finufft'),'MATLAB:validators:mustBeMember')
        end

        function finufftMatchesAnalyticalPeriodicFieldWhenAvailable(testCase)
            if exist('finufft_plan','class') ~= 8
                return
            end

            wvt = WVTransformBarotropicQG([4000 3000],[8 6],latitude=45,shouldAntialias=false);
            field = cos(2*pi*wvt.X/wvt.Lx) + 0.4*sin(2*pi*wvt.Y/wvt.Ly);
            x = [0 137 3999 4123];
            y = [0 211 2999 -77];
            expected = cos(2*pi*x/wvt.Lx) + 0.4*sin(2*pi*y/wvt.Ly);

            dftField = wvt.transformFromSpatialDomainToDFTGrid(field);
            fromDFT = wvt.transformToSpatialDomainFromDFTGridAtPosition(dftField,x,y);
            wvField = wvt.transformFromSpatialDomainWithFourier(field);
            fromWV = wvt.transformToSpatialDomainWithFourierAtPosition(wvField,x,y);

            testCase.verifyEqual(fromDFT,expected,AbsTol=2e-12)
            testCase.verifyEqual(fromWV,expected,AbsTol=2e-12)
        end
    end
end
