classdef TestCompiledPreview < matlab.unittest.TestCase
    methods (Test, TestTags="smoke")
        function matlabIsTheSideEffectFreeDefault(testCase)
            wvt = WVTransformConstantStratification([1200 900 100],[8 6 5],isHydrostatic=true);
            cleanup = onCleanup(@()delete(wvt));
            testCase.verifyEqual(wvt.computationalBackend,"matlab");
            metadata = wvt.computationalBackendMetadata;
            testCase.verifyEqual(metadata.schemaVersion,"1.0.0");
            testCase.verifyEqual(metadata.requestedBackend,"matlab");
            testCase.verifyEqual(metadata.activeBackend,"matlab");
            testCase.verifyEqual(metadata.provider.id,"matlab-builtin");
            testCase.verifyEqual(metadata.storage.status,"not-estimated");
            clear cleanup
        end

        function invalidBackendIsRejected(testCase)
            testCase.verifyError(@()WVTransformConstantStratification([1200 900 100],[8 6 5],computationalBackend="automatic"),"MATLAB:validators:mustBeMember");
        end
    end

    methods (Test, TestTags="full")
        function matlabResolutionAndExplicitAntialiasingRemainMatlab(testCase)
            wvt = WVTransformConstantStratification([1200 900 100],[8 6 5],isHydrostatic=true);
            cleanup = onCleanup(@()delete(wvt));
            resized = wvt.waveVortexTransformWithResolution([10 8 7]);
            resizedCleanup = onCleanup(@()delete(resized));
            explicit = wvt.waveVortexTransformWithExplicitAntialiasing();
            explicitCleanup = onCleanup(@()delete(explicit));
            testCase.verifyEqual(resized.computationalBackend,"matlab");
            testCase.verifyEqual(explicit.computationalBackend,"matlab");
            clear explicitCleanup resizedCleanup cleanup
        end
    end
end
