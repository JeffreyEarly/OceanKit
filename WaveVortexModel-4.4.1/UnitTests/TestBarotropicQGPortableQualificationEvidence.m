classdef TestBarotropicQGPortableQualificationEvidence < matlab.unittest.TestCase
    methods (Test, TestTags="full")
        function evidenceNamesExecutableQualificationGates(testCase)
            repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            path = fullfile(repositoryRoot,"PortableRuntime","qualification", ...
                "barotropic-qg-v1.json");
            evidence = jsondecode(fileread(path));
            testCase.verifyEqual(string(evidence.schemaIdentifier), ...
                "wave-vortex-barotropic-qg-qualification-v1")
            testCase.verifyEqual(evidence.schemaVersion,1)
            testCase.verifyEqual(evidence.issue,286)
            testCase.verifyEqual(string(evidence.requestAuthor), ...
                "WVModel.writePortableRunRequest")
            testCase.verifyEqual(string(evidence.providers),["reference";"native-fftw"])
            testCase.verifyEqual(evidence.acceptance.coefficientRelativeErrorMaximum,1e-12)
            testCase.verifyEqual(evidence.acceptance.fieldRelativeErrorMaximum,1e-12)
            testCase.verifyEqual(string(evidence.acceptance.compactCoefficientFamilies),"A0")
            testCase.verifyEqual(string(evidence.acceptance.forbiddenCoefficientFamilies),["Ap";"Am"])
            testCase.verifyEqual(evidence.acceptance.persistentFullHermitianBytes,0)
            testCase.verifyFalse(evidence.acceptance.fallbackAllowed)
            testCase.verifyLessThanOrEqual(evidence.acceptance.maximumFFTThreads,16)
            requiredMatlab = string(evidence.requiredTests.matlabCppMatlab);
            testCase.verifyTrue(all(contains(requiredMatlab,"TestPortableRuntimeCompatibility/")))
            testCase.verifyNumElements(requiredMatlab,5)
            testCase.verifyEqual(numel(evidence.forcings),7)
            testCase.verifyEqual(numel(evidence.transactionalRejections),6)
        end
    end
end
