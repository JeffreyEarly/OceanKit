classdef TestPortableStratifiedQGQualificationEvidence < matlab.unittest.TestCase
    properties (SetAccess=private)
        report (1,1) struct
    end
    methods (TestClassSetup)
        function loadRecordedEvidence(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"tools")));
            testCase.report = jsondecode(fileread(fullfile(root,"PortableRuntime","qualification","stratified-qg-apple-silicon-v1.json")));
        end
    end
    methods (Test,TestTags="full")
        function measuredEvidenceIsComplete(testCase)
            validation = validatePortableStratifiedQGQualification(testCase.report,evidenceScope="historical");
            testCase.verifyEqual(validation.scope,"historical-workload");
            testCase.verifyFalse(validation.currentReadiness);
            testCase.verifyEqual(numel(testCase.report.rows),41);
            testCase.verifyEqual(numel(testCase.report.continuations),12);
            testCase.verifyEqual(numel(testCase.report.lifecycle),6);
        end
        function missingOrContradictoryEvidenceIsRejected(testCase)
            original = testCase.report;
            bad = original; bad.rows(end)=[]; testCase.verifyInvalid(bad);
            bad = original; bad.rows(end)=bad.rows(1); testCase.verifyInvalid(bad);
            bad = original; bad.rows(1).passed=false; testCase.verifyInvalid(bad);
            bad = original; bad.rows(1).testNames="TestPortableStratifiedQG/fieldSamplingMatchesMatlab"; testCase.verifyInvalid(bad);
            bad = original; bad.catalogSHA256="stale"; testCase.verifyInvalid(bad);
            bad = original; bad.tests(1).passed=false; testCase.verifyInvalid(bad);
            bad = original; bad.tests(end)=[]; testCase.verifyInvalid(bad);
            bad = original; bad.continuations(end)=[]; testCase.verifyInvalid(bad);
            bad = original; bad.continuations(end)=bad.continuations(1); testCase.verifyInvalid(bad);
            bad = original; bad.lifecycle(end)=bad.lifecycle(1); testCase.verifyInvalid(bad);
        end
        function numericalAndMemoryFailuresCannotClaimReadiness(testCase)
            bad = testCase.report; bad.continuations(1).whole.coefficients=NaN; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).whole.energy=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).definition.duration=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).whole.graph.dense.ordinal=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).runtime(1).execution.noFallback=false; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.lifecycle(1).retainedGrowthBytes=8; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.lifecycle(1).measurements(1).retainedBytesAfter=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.lifecycle(1).scientificOwnersReleased=false; testCase.verifyInvalid(bad);
            bad = testCase.report; index=find(string({bad.lifecycle.provider})=="native-fftw",1); bad.lifecycle(index).preparedStepAllocations=1; testCase.verifyInvalid(bad);
        end
    end
    methods (Access=private)
        function verifyInvalid(testCase,report)
            testCase.verifyError(@()validatePortableStratifiedQGQualification(report,evidenceScope="historical"),"WaveVortexModel:InvalidSQGQualification");
        end
    end
end
