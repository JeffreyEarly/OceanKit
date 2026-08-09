classdef TestWVFourierStorageLayoutIntegrationBenchmark < matlab.unittest.TestCase
    properties
        benchmarkFolder
        temporaryFolder
    end

    methods (TestClassSetup)
        function addBenchmarkPath(testCase)
            repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.benchmarkFolder = fullfile(repositoryRoot,"Benchmarks");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testCase.benchmarkFolder));
        end
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function reducedRunComparesProductionWithFrozenReference(testCase)
            caseId = "full-layout-64x48x17-antialias-0";
            result = runWVFourierStorageLayoutIntegrationBenchmark(caseIds=caseId,shouldWriteArtifacts=false,runId="storage-layout-integration-test");

            testCase.verifyEqual(result.status,"complete");
            testCase.verifyEqual(result.reference.issue,69);
            testCase.verifyEqual(strlength(result.reference.sha256),64);
            testCase.verifyTrue(result.readiness.passed);
            testCase.verifyFalse(result.readiness.completeGateSet);
            testCase.verifyEqual(result.readiness.gateCaseCount,0);
            testCase.verifyEqual(result.operationIds,["extract" "insert-primary" "insert-conjugate" "insert-complete" "forward-complete" "inverse-complete"]);

            benchmarkCase = result.cases;
            testCase.verifyEqual(benchmarkCase.status,"complete");
            testCase.verifyFalse(benchmarkCase.isGate);
            testCase.verifyEqual(benchmarkCase.mappingMethod,"two-dimensional-rows");
            testCase.verifyEqual(benchmarkCase.fourierStorageSize,[64 48]);
            testCase.verifyFalse(benchmarkCase.legacyMappingsAreMaterialized);
            testCase.verifyEqual(benchmarkCase.legacyMappingBytes,0);
            testCase.verifySize(benchmarkCase.warmupSchedules,[2 6]);
            testCase.verifySize(benchmarkCase.sampleSchedules,[7 6]);
            testCase.verifyEqual(numel(unique(benchmarkCase.sampleSchedules(1,:))),6);
            testCase.verifyNotEqual(benchmarkCase.sampleSchedules(1,1),benchmarkCase.sampleSchedules(2,1));
            for operation = benchmarkCase.operations
                testCase.verifyNumElements(operation.rawSeconds,7);
                testCase.verifyTrue(all(isfinite(operation.rawSeconds)));
                testCase.verifyEqual(operation.medianSeconds,median(operation.rawSeconds),AbsTol=eps(operation.medianSeconds));
                testCase.verifyGreaterThan(operation.referenceMedianSeconds,0);
                testCase.verifyEqual(operation.relativeToReference,operation.medianSeconds/operation.referenceMedianSeconds,AbsTol=10*eps(operation.relativeToReference));
                testCase.verifyLessThanOrEqual(operation.relativeError,1e-12);
                testCase.verifyTrue(operation.correctnessPassed);
                testCase.verifyFalse(operation.performanceGateApplies);
                testCase.verifyTrue(operation.performancePassed);
            end
        end

        function artifactsRoundTripAndStateIsRestored(testCase)
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            outputDirectory = fullfile(testCase.temporaryFolder,"artifacts");
            result = runWVFourierStorageLayoutIntegrationBenchmark(caseIds="full-layout-64x48x17-antialias-1",outputDirectory=outputDirectory,runId="storage-layout-integration-artifact");

            testCase.verifyTrue(isfile(fullfile(outputDirectory,"benchmark.json")));
            testCase.verifyTrue(isfile(fullfile(outputDirectory,"summary.md")));
            decoded = jsondecode(fileread(fullfile(outputDirectory,"benchmark.json")));
            testCase.verifyEqual(string(decoded.runId),result.runId);
            testCase.verifyEqual(string(decoded.cases.mappingMethod),"two-dimensional-rows");
            summary = string(fileread(fullfile(outputDirectory,"summary.md")));
            testCase.verifySubstring(summary,"## Timing and frozen-baseline comparison");
            testCase.verifySubstring(summary,"## Storage contract");
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
        end

        function releaseCandidateArtifactPassed(testCase)
            artifactPath = fullfile(testCase.benchmarkFolder,"results","reference", ...
                "transform-layout-v4.2.1-release-m5-max-r2026a-builtin","benchmark.json");
            artifact = jsondecode(fileread(artifactPath));

            testCase.verifyEqual(string(artifact.status),"complete");
            testCase.verifyEqual(string(artifact.runId),"transform-layout-v4.2.1-release-m5-max-r2026a-builtin");
            testCase.verifyEqual(strlength(string(artifact.environment.sourceCommit)),40);
            testCase.verifyEqual(strlength(string(artifact.environment.sourceTree)),40);
            testCase.verifyTrue(artifact.readiness.passed);
            testCase.verifyTrue(artifact.readiness.completeGateSet);
            testCase.verifyEqual(artifact.readiness.gateCaseCount,4);
            gateCases = artifact.cases([artifact.cases.isGate]);
            operations = [gateCases.operations];
            testCase.verifyTrue(all([operations.correctnessPassed]));
            testCase.verifyTrue(all([operations.performancePassed]));
        end

        function invalidCaseRestoresState(testCase)
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            testCase.verifyError(@()runWVFourierStorageLayoutIntegrationBenchmark(caseIds="missing-layout-case",shouldWriteArtifacts=false),"WaveVortexBenchmark:UnknownCase");
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
        end
    end
end
