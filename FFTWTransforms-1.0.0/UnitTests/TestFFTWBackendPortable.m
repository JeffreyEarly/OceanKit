classdef TestFFTWBackendPortable < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRepositoryToPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(TestFFTWBackendPortable.repositoryRoot));
        end
    end

    methods (Test)
        function testLinuxR2025bStructuredUnavailability(testCase)
            isPortableEnvironment = string(version('-release')) == "2025b" && string(computer('arch')) == "glnxa64";
            testCase.assumeTrue(isPortableEnvironment,"This test targets the portable Linux/R2025b CI environment.");

            before = TestFFTWBackendPortable.generatedMexFiles;
            testCase.verifyWarningFree(@() FFTWBackend.capabilities());
            capabilities = FFTWBackend.capabilities();
            testCase.verifyEqual(capabilities.status,"unavailable");
            testCase.verifyFalse(capabilities.isAvailable);
            testCase.verifyFalse(capabilities.isComplete);
            testCase.verifyEqual(capabilities.reason.code,"unsupported-release");
            testCase.verifyFalse(capabilities.build.attempted);
            testCase.verifyFalse(capabilities.build.isPossible);
            TestFFTWBackendPortable.verifyUnavailableFeatures(testCase,capabilities,"unsupported-release");

            testCase.verifyWarningFree(@() FFTWBackend.build());
            buildResult = FFTWBackend.build();
            testCase.verifyEqual(buildResult.status,"unavailable");
            testCase.verifyTrue(buildResult.build.attempted);
            testCase.verifyFalse(buildResult.build.succeeded);
            testCase.verifyFalse(buildResult.build.installed);
            testCase.verifyEqual(buildResult.build.reason.code,"unsupported-platform");
            testCase.verifyEqual(TestFFTWBackendPortable.generatedMexFiles,before);
        end

        function testUnsupportedArchitectureIsActionable(testCase)
            context = TestableFFTWBackend.context();
            context.release = context.provider.supportedReleases;
            context.architecture = "glnxa64";
            result = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(result.status,"unavailable");
            testCase.verifyEqual(result.reason.code,"unsupported-architecture");
            testCase.verifyTrue(contains(result.reason.message,"maca64"));
            testCase.verifyFalse(result.build.isPossible);
            TestFFTWBackendPortable.verifyUnavailableFeatures(testCase,result,"unsupported-architecture");
        end
    end

    methods (Static, Access=private)
        function verifyUnavailableFeatures(testCase,capabilities,reasonCode)
            for featureName = ["r2c","c2r","dct1","dst1"]
                feature = capabilities.features.(featureName);
                testCase.verifyFalse(feature.isAvailable);
                testCase.verifyFalse(feature.selfTestPassed);
                testCase.verifyEqual(feature.reason.code,reasonCode);
                testCase.verifyNotEqual(feature.reason.message,"");
            end
        end

        function root = repositoryRoot()
            root = fileparts(fileparts(mfilename('fullpath')));
        end

        function names = generatedMexFiles()
            files = dir(fullfile(TestFFTWBackendPortable.repositoryRoot,"*."+mexext));
            names = sort(string({files.name}));
        end
    end
end
