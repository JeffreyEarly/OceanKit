classdef TestWaveVortexBenchmark < matlab.unittest.TestCase
    properties
        repositoryRoot
        benchmarkFolder
        temporaryFolder
    end

    methods (TestClassSetup)
        function addBenchmarkPath(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.benchmarkFolder = fullfile(testCase.repositoryRoot,"Benchmarks");
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
        function suiteRegistryIsVersionedAndTransformSpecific(testCase)
            suites = waveVortexBenchmarkSuites(["smoke-v1" "core-v1" "scaling-standard-v1" "scaling-large-v1"]);
            testCase.verifyEqual(string({suites.id}),["smoke-v1" "core-v1" "scaling-standard-v1" "scaling-large-v1"]);
            testCase.verifyEqual([suites.version],[1 1 1 1]);
            testCase.verifyEqual(numel(suites(1).cases),6);
            testCase.verifyEqual(numel(suites(2).cases),4);
            testCase.verifyEqual(numel(suites(3).cases),34);
            testCase.verifyEqual(numel(suites(4).cases),19);
            testCase.verifyTrue(all(ismember(["constant-nonhydrostatic" "constant-hydrostatic" "hydrostatic" "boussinesq" "stratified-qg" "barotropic-qg"],unique(string({suites(3).cases.transformId})))));
            testCase.verifyError(@()waveVortexBenchmarkSuites("missing-v1"),"WaveVortexBenchmark:UnknownSuite");
            testCase.verifyError(@()waveVortexBenchmarkBackends("missing"),"WaveVortexBenchmark:UnknownBackend");
        end

        function scoringUsesEqualFamilyWeights(testCase)
            scores = waveVortexBenchmarkScores([2;4;10],[1;2;10],["a";"a";"b"]);
            testCase.verifyEqual(scores.caseScores,[200;200;100],AbsTol=1e-12);
            testCase.verifyEqual(scores.familyScores,[200;100],AbsTol=1e-12);
            testCase.verifyEqual(scores.suiteScore,sqrt(20000),AbsTol=1e-12);
            testCase.verifyError(@()waveVortexBenchmarkScores([1;2],1,["a";"b"]),"WaveVortexBenchmark:ScoreSizeMismatch");
        end

        function smokeSuiteRunsWithProductionCacheSemantics(testCase)
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            results = runWaveVortexBenchmark(suites="smoke-v1",shouldMeasureMemory=false,shouldWriteArtifacts=false);
            testCase.verifyEqual(results.status,"complete");
            testCase.verifyEqual(string({results.suites.cases.status}),repmat("complete",1,6));
            for benchmarkCase = results.suites.cases
                backend = benchmarkCase.backends;
                testCase.verifyNumElements(backend.rawSeconds,benchmarkCase.sampleCount);
                testCase.verifyTrue(all(isfinite(backend.rawSeconds)));
                testCase.verifyGreaterThan(backend.medianSeconds,0);
                testCase.verifyGreaterThan(backend.firstCallSeconds,0);
                testCase.verifyGreaterThan(backend.sameStateCacheHitSeconds,0);
                testCase.verifyLessThanOrEqual(backend.relativeError,1e-12);
                testCase.verifyTrue(backend.correctnessPassed);
                testCase.verifyEqual(backend.memory.status,"not-requested");
            end
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
        end

        function artifactsAndReferenceRoundTrip(testCase)
            runFolder = fullfile(testCase.temporaryFolder,"run");
            referenceFolder = fullfile(testCase.temporaryFolder,"reference");
            results = runWaveVortexBenchmark(suites="smoke-v1",caseIds="smoke-constant-hydrostatic",outputDirectory=runFolder,referenceDirectory=referenceFolder,shouldMeasureMemory=false,shouldCreateReference=true,runId="test-run");
            testCase.verifyTrue(isfile(fullfile(runFolder,"benchmark.json")));
            testCase.verifyTrue(isfile(fullfile(runFolder,"summary.md")));
            referencePath = fullfile(referenceFolder,"smoke-v1-m5-max-r2026a-builtin");
            testCase.verifyTrue(isfile(fullfile(referencePath,"benchmark.json")));
            decoded = jsondecode(fileread(fullfile(runFolder,"benchmark.json")));
            testCase.verifyEqual(string(decoded.schemaVersion),"1.0.0");
            testCase.verifyEqual(decoded.suites.cases.backends.caseScore,100);
            summary = fileread(fullfile(runFolder,"summary.md"));
            testCase.verifySubstring(summary,"## Suite scores");
            testCase.verifySubstring(summary,"## Family scores");
            testCase.verifySubstring(summary,"## Timing and scores");
            testCase.verifySubstring(summary,"## Construction and cache diagnostics");
            testCase.verifySubstring(summary,"## Memory");
            testCase.verifyEqual(results.runId,"test-run");
            testCase.verifyEqual(results.status,"partial");
            testCase.verifyFalse(results.suites.selectionIsComplete);
        end

        function freshProcessMemoryIsRecorded(testCase)
            results = runWaveVortexBenchmark(suites="smoke-v1",caseIds="smoke-barotropic-qg",shouldMeasureMemory=true,shouldWriteArtifacts=false);
            memory = results.suites.cases.backends.memory;
            diagnostic = "Memory worker failed: " + string(memory.failure.identifier) + ": " + string(memory.failure.message);
            if string(memory.status) == "failed" && strcmp(getenv("GITHUB_ACTIONS"),"true")
                testCase.verifyEqual(string(memory.failure.identifier),"WaveVortexBenchmark:MemoryWorkerFailed",diagnostic);
                licenseMessages = ["License checkout failed" "Unable to find a license for MATLAB"];
                testCase.verifyTrue(any(contains(string(memory.failure.message),licenseMessages)),diagnostic);
                return
            end
            testCase.verifyEqual(string(memory.status),"complete",diagnostic);
            testCase.verifyGreaterThan(memory.baselineBytes,0);
            testCase.verifyGreaterThanOrEqual(memory.persistentIncrementBytes,0);
            testCase.verifyGreaterThanOrEqual(memory.peakIncrementBytes,memory.persistentIncrementBytes);
            testCase.verifyNotEmpty(memory.provider);
        end
    end
end
