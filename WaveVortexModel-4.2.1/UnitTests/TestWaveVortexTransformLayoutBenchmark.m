classdef TestWaveVortexTransformLayoutBenchmark < matlab.unittest.TestCase
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
        function registryDefinesEightDiagnosticCases(testCase)
            suite = waveVortexBenchmarkSuites("transform-layout-v1");
            testCase.verifyEqual(suite.kind,"transform-layout");
            testCase.verifyFalse(suite.isScored);
            testCase.verifyEqual(numel(suite.cases),8);
            testCase.verifyEqual(reshape([suite.cases.Nxyz],3,[]).',[64 48 17;64 48 17;65 63 17;65 63 17;256 256 65;256 256 65;512 512 129;512 512 129]);
            testCase.verifyEqual([suite.cases.shouldAntialias],logical([0 1 0 1 0 1 0 1]));
            testCase.verifyEqual([suite.cases.sampleCount],[7 7 7 7 7 7 3 3]);
            testCase.verifyEqual([suite.cases.warmupCount],2*ones(1,8));
        end

        function reducedOddEvenSuiteMatchesProduction(testCase)
            suite = reducedSuite();
            result = runWaveVortexTransformLayoutSuite(suite,1e-12,testCase.repositoryRoot);
            testCase.verifyEqual(result.status,"complete");
            testCase.verifyEqual(result.metadata.strategyIds,["wv-sorted-linear" "dft-sorted-linear" "two-dimensional-rows" "per-plane-linear"]);
            testCase.verifyEqual(result.metadata.operationIds,["extract" "insert-primary" "insert-conjugate" "insert-complete" "forward-complete" "inverse-complete"]);
            testCase.verifyEqual(result.metadata.assignmentOrder,["primary" "conjugate"]);
            testCase.verifyEqual(numel(result.metadata.productionSources),3);
            testCase.verifyTrue(all(strlength(string({result.metadata.productionSources.sha256})) == 64));

            expectedPairCount = 4*6;
            for benchmarkCase = result.cases
                testCase.verifyEqual(benchmarkCase.status,"complete");
                testCase.verifySize(benchmarkCase.warmupSchedules,[1 expectedPairCount]);
                testCase.verifySize(benchmarkCase.sampleSchedules,[2 expectedPairCount]);
                testCase.verifyEqual(numel(unique(benchmarkCase.sampleSchedules(1,:))),expectedPairCount);
                testCase.verifyNotEqual(benchmarkCase.sampleSchedules(1,1),benchmarkCase.sampleSchedules(2,1));
                testCase.verifyEqual(string({benchmarkCase.strategies.id}),result.metadata.strategyIds);
                testCase.verifyEqual(string({benchmarkCase.selections.operationId}),result.metadata.operationIds);
                for strategy = benchmarkCase.strategies
                    testCase.verifyEqual(strategy.assignmentSteps,["primary" "conjugate"]);
                    testCase.verifyTrue(strategy.persistentBufferReused);
                    testCase.verifyFalse(strategy.timedBufferClearing);
                    testCase.verifyTrue(strategy.sourceArraysUnchanged);
                    testCase.verifyEqual(strategy.copyObservation,"unavailable");
                    testCase.verifyEqual(strategy.mappingBytes,sum([strategy.mappingArrays.bytes]));
                    for mapping = strategy.mappingArrays
                        testCase.verifyEqual(mapping.bytes,8*prod(mapping.shape));
                    end
                    testCase.verifyTrue(all(string({strategy.mappingArrays.class}) == "uint64"));
                    testCase.verifyEqual(string({strategy.operations.id}),result.metadata.operationIds);
                    for operation = strategy.operations
                        testCase.verifyNumElements(operation.rawSeconds,benchmarkCase.sampleCount);
                        testCase.verifyTrue(all(isfinite(operation.rawSeconds)));
                        testCase.verifyEqual(operation.medianSeconds,median(operation.rawSeconds),AbsTol=eps(operation.medianSeconds));
                        testCase.verifyLessThanOrEqual(operation.relativeError,1e-12);
                        testCase.verifyTrue(operation.correctnessPassed);
                        testCase.verifyEqual(operation.resultBytes,operation.resultStorage.bytes);
                        if ismember(operation.id,["extract" "forward-complete"])
                            testCase.verifyEqual(operation.resultStorage.shape(1),benchmarkCase.Nxyz(3));
                            testCase.verifyEqual(operation.resultStorage.bytes,16*prod(operation.resultStorage.shape));
                            testCase.verifyTrue(operation.resultStorage.isComplex);
                        elseif operation.id == "inverse-complete"
                            testCase.verifyEqual(operation.resultStorage.shape,benchmarkCase.Nxyz);
                            testCase.verifyFalse(operation.resultStorage.isComplex);
                        else
                            testCase.verifyEqual(operation.resultStorage.bytes,0);
                        end
                    end
                end
            end
        end

        function currentPathPreferenceUsesThreePercentBand(testCase)
            ids = ["wv-sorted-linear" "candidate"];
            retained = waveVortexTransformLayoutPreference(ids,[1.029 1]);
            testCase.verifyEqual(retained.strictFastestStrategy,"candidate");
            testCase.verifyEqual(retained.preferredStrategy,"wv-sorted-linear");
            testCase.verifyTrue(retained.currentRetained);
            displaced = waveVortexTransformLayoutPreference(ids,[1.031 1]);
            testCase.verifyEqual(displaced.preferredStrategy,"candidate");
            testCase.verifyFalse(displaced.currentRetained);
            tied = waveVortexTransformLayoutPreference(ids,[1 1]);
            testCase.verifyEqual(tied.strictFastestStrategy,"wv-sorted-linear");
            testCase.verifyEqual(tied.preferredStrategy,"wv-sorted-linear");
            testCase.verifyError(@()waveVortexTransformLayoutPreference(ids,[NaN Inf]),"WaveVortexBenchmark:NoValidLayoutStrategy");
        end

        function sharedEntryPointWritesDiagnosticArtifactsAndRestoresState(testCase)
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            outputFolder = fullfile(testCase.temporaryFolder,"run");
            results = runWaveVortexBenchmark(suites="transform-layout-v1",caseIds="full-layout-64x48x17-antialias-0",outputDirectory=outputFolder,shouldMeasureMemory=true,runId="layout-test");
            testCase.verifyTrue(isfile(fullfile(outputFolder,"benchmark.json")));
            testCase.verifyTrue(isfile(fullfile(outputFolder,"summary.md")));
            decoded = jsondecode(fileread(fullfile(outputFolder,"benchmark.json")));
            testCase.verifyEqual(string(decoded.suites.kind),"transform-layout");
            testCase.verifyEqual(numel(decoded.suites.cases.strategies),4);
            summary = string(fileread(fullfile(outputFolder,"summary.md")));
            testCase.verifySubstring(summary,"## Transform-layout diagnostic");
            testCase.verifySubstring(summary,"### Extraction and complete-forward winners");
            testCase.verifySubstring(summary,"### Insertion and complete-inverse winners");
            testCase.verifySubstring(summary,"### Mapping-array and working-array storage");
            testCase.verifySubstring(summary,"### Correctness and observable copy semantics");
            testCase.verifyFalse(contains(summary,"fresh-process"));
            testCase.verifyEqual(results.configuration.shouldMeasureMemory,true);
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
        end

        function sharedEntryPointRestoresStateAfterFailure(testCase)
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            testCase.verifyError(@()runWaveVortexBenchmark(suites="transform-layout-v1",caseIds="missing-layout-case",shouldWriteArtifacts=false),"WaveVortexBenchmark:UnknownCase");
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
        end
    end
end

function suite = reducedSuite()
suite = waveVortexBenchmarkSuites("transform-layout-v1");
suite.cases = suite.cases(1:4);
sizes = [16 12 5;16 12 5;17 15 5;17 15 5];
for iCase = 1:numel(suite.cases)
    suite.cases(iCase).id = "reduced-layout-" + iCase;
    suite.cases(iCase).Nxyz = sizes(iCase,:);
    suite.cases(iCase).warmupCount = 1;
    suite.cases(iCase).sampleCount = 2;
end
end
