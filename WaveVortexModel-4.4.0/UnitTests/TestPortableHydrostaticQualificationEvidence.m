classdef TestPortableHydrostaticQualificationEvidence < matlab.unittest.TestCase
    properties (SetAccess=private)
        report (1,1) struct
    end
    methods (TestClassSetup)
        function loadRecordedEvidence(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"tools")));
            testCase.report = jsondecode(fileread(fullfile(root,"PortableRuntime","qualification","hydrostatic-apple-silicon-v1.json")));
        end
    end
    methods (Test,TestTags="full")
        function measuredEvidenceIsComplete(testCase)
            validation = validatePortableHydrostaticQualification(testCase.report,evidenceScope="historical");
            testCase.verifyEqual(validation.scope,"historical-workload");
            testCase.verifyFalse(validation.currentReadiness);
            testCase.verifyEqual(numel(testCase.report.rows),47);
            testCase.verifyEqual(numel(testCase.report.continuations),12);
            testCase.verifyEqual(numel(testCase.report.lifecycle),7);
        end
        function referenceOnlyScopeStillRequiresEveryContinuation(testCase)
            referenceReport = testCase.report;
            referenceReport.providers = "reference";
            referenceReport.tests = referenceReport.tests(string({referenceReport.tests.name})~="TestPortableStratifiedQGQualification/lifecycleAndStorageRemainBounded");
            referenceReport.rows = referenceReport.rows(string({referenceReport.rows.provider})=="reference");
            referenceReport.continuations = referenceReport.continuations(string({referenceReport.continuations.provider})=="reference");
            grids = reshape([referenceReport.lifecycle.grid],3,[])';
            referenceReport.lifecycle = referenceReport.lifecycle(string({referenceReport.lifecycle.provider})=="reference" & ismember(grids,[8 6 9;12 10 13],"rows")');
            validatePortableHydrostaticQualification(referenceReport,evidenceScope="historical");
            testCase.verifyEqual(numel(referenceReport.continuations),6);
            testCase.verifyEqual(numel(referenceReport.lifecycle),2);
            bad = referenceReport; bad.continuations(end)=[]; testCase.verifyInvalid(bad);
            bad = referenceReport; bad.lifecycle(end)=[]; testCase.verifyInvalid(bad);
        end
        function contractsCannotClaimCompleteContinuationQualification(testCase)
            referenceReport = testCase.report;
            referenceReport.providers = "reference";
            referenceReport.tests = referenceReport.tests(string({referenceReport.tests.name})~="TestPortableStratifiedQGQualification/lifecycleAndStorageRemainBounded");
            referenceReport.rows = referenceReport.rows(string({referenceReport.rows.provider})=="reference");
            referenceReport.continuations = referenceReport.continuations(string({referenceReport.continuations.provider})=="reference");
            grids = reshape([referenceReport.lifecycle.grid],3,[])';
            referenceReport.lifecycle = referenceReport.lifecycle(string({referenceReport.lifecycle.provider})=="reference" & ismember(grids,[8 6 9;12 10 13],"rows")');
            referenceReport.schemaIdentifier = "wave-vortex-hydrostatic-contracts-v1";
            referenceReport.tests = referenceReport.tests(string({referenceReport.tests.name})~="TestPortableHydrostaticQualification/longerContinuationMatchesMatlab");
            completeContinuations = referenceReport.continuations;
            referenceReport.continuations = [];
            validatePortableHydrostaticQualification(referenceReport,evidenceScope="historical");
            bad = referenceReport; bad.schemaIdentifier="wave-vortex-hydrostatic-qualification-v1"; testCase.verifyInvalid(bad);
            bad = referenceReport; bad.continuations=completeContinuations; testCase.verifyInvalid(bad);
            bad = referenceReport; bad.tests(end)=[]; testCase.verifyInvalid(bad);
        end
        function missingOrContradictoryEvidenceIsRejected(testCase)
            original = testCase.report;
            bad = original; bad.rows(end)=[]; testCase.verifyInvalid(bad);
            bad = original; bad.rows(end)=bad.rows(1); testCase.verifyInvalid(bad);
            bad = original; bad.rows(1).passed=false; testCase.verifyInvalid(bad);
            bad = original; bad.rows(1).testNames="TestPortableHydrostatic/fieldSamplingMatchesMatlab"; testCase.verifyInvalid(bad);
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
            bad = testCase.report; bad.lifecycle(1).Nj=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.lifecycle(1).finalStateFinite=false; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.lifecycle(1).schemaIdentifier="wave-vortex-sqg-lifecycle-v1"; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).whole.coefficientFamilies.Ap=1; testCase.verifyInvalid(bad);
            bad = testCase.report; bad.continuations(1).runtime(1).integrationRequest.activeMethod="wrong"; testCase.verifyInvalid(bad);
            bad = testCase.report; index=find(string({bad.lifecycle.provider})=="native-fftw",1); bad.lifecycle(index).preparedStepAllocations=1; testCase.verifyInvalid(bad);
        end
    end
    methods (Access=private)
        function verifyInvalid(testCase,report)
            testCase.verifyError(@()validatePortableHydrostaticQualification(report,evidenceScope="historical"),"WaveVortexModel:InvalidHydrostaticQualification");
        end
    end
end
