classdef TestCompiledPreviewBenchmark < matlab.unittest.TestCase
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
        function availabilityDecisionUsesSpeedCorrectnessAndLifecycle(testCase)
            comparison = comparisonFixture();
            decision = compiledPreviewBenchmarkDecision(comparison);
            testCase.verifyEqual(decision.status,"PREVIEW-AVAILABLE");
            testCase.verifyTrue(decision.available);

            comparison(2).compiledSpeedup = 1.249;
            testCase.verifyEqual(compiledPreviewBenchmarkDecision(comparison).status,"PREVIEW-NOT-AVAILABLE");
            comparison(2).compiledSpeedup = 1.5;
            comparison(1).maximumRelativeError = 1.01e-12;
            testCase.verifyEqual(compiledPreviewBenchmarkDecision(comparison).status,"PREVIEW-NOT-AVAILABLE");
            comparison(1).maximumRelativeError = 1e-15;
            comparison(1).exactRetainedRatio = 4;
            comparison(1).operationPeakRSSRatio = 4;
            testCase.verifyEqual(compiledPreviewBenchmarkDecision(comparison).status,"PREVIEW-AVAILABLE");
        end

        function publicNormalizationCarriesExactAndRSSMemory(testCase)
            raw = rawFixture();
            rawPath = fullfile(testCase.temporaryFolder,"compiled-preview-benchmark.json");
            writelines(jsonencode(raw),rawPath);
            [matlabDataset,compiledDataset] = publishedWaveVortexBenchmarksFromCompiledPreviewArtifact(rawPath,platformId="m5-max",platformName="Apple M5 Max",provenancePath="Benchmarks/results/reference/compiled-preview-v1-m5-max-r2026a/compiled-preview-benchmark.json");
            testCase.verifyEqual(matlabDataset.datasetId,"core-v1--matlab-builtin--m5-max--20260812T120000Z");
            testCase.verifyEqual(compiledDataset.datasetId,"core-v1--cpp-native-fftw--m5-max--20260812T120000Z");
            testCase.verifyEqual(compiledDataset.implementation.id,"cpp");
            testCase.verifyEqual(compiledDataset.toolchain.kind,"cpp");
            testCase.verifyEqual(compiledDataset.cases{1}.memory.exactRetainedBytes,240);
            testCase.verifyEqual(matlabDataset.cases{1}.memory.exactRetainedBytes,120);
            testCase.verifyEqual(compiledDataset.cases{1}.memory.peakIncrementBytes,100);
            testCase.verifyEqual(compiledDataset.cases{1}.timing.samplesSeconds,[0.4 0.5]);
        end

        function registryContainsCompiledOnlyForConstantStratification(testCase)
            backends = waveVortexBenchmarkBackends(["builtin" "compiled"]);
            testCase.verifyEqual(string({backends.id}),["builtin" "compiled"]);
            suite = waveVortexBenchmarkSuites("smoke-v1");
            nonconstant = suite.cases(find(~startsWith(string({suite.cases.transformId}),"constant-"),1));
            testCase.verifyError(@()createWaveVortexBenchmarkTransform(nonconstant,"compiled"),"WaveVortexBenchmark:UnsupportedBackend");
        end

    end

    methods (Test,TestTags="optional")
        function reducedPublicBenchmarkRunsBothBackends(testCase)
            if ~isCanonicalNativePlatform
                testCase.verifyError(@()runCompiledPreviewBenchmark(shouldWriteArtifacts=false),"WaveVortexBenchmark:CompiledPreviewUnsupportedPlatform")
                return
            end
            result = runCompiledPreviewBenchmark(sizes=[16 12 9],hydrostatic=[true false],processRunCount=1,warmupCount=0,mediumSampleCount=1,largeSampleCount=1,samplingIntervalSeconds=0.01,plateauSeconds=0.03,outputHoldSeconds=0.02,shouldWriteArtifacts=false);
            testCase.verifyEqual(result.status,"complete");
            testCase.verifyEqual(string({result.comparison.status}),["complete" "complete"]);
            testCase.verifyLessThanOrEqual(max([result.comparison.maximumRelativeError]),1e-12);
            testCase.verifyTrue(all([result.comparison.nativeExecutionPassed]));
            testCase.verifyTrue(all([result.comparison.lifecyclePassed]));
        end

        function injectedWorkerFailureWritesPartialArtifactAndRestoresState(testCase)
            if ~isCanonicalNativePlatform
                testCase.verifyError(@()runCompiledPreviewBenchmark(shouldWriteArtifacts=false),"WaveVortexBenchmark:CompiledPreviewUnsupportedPlatform")
                return
            end
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            outputDirectory = fullfile(testCase.temporaryFolder,"partial");
            testCase.verifyError(@()runCompiledPreviewBenchmark( ...
                sizes=[16 12 9],hydrostatic=true,processRunCount=1, ...
                warmupCount=0,mediumSampleCount=1,largeSampleCount=1, ...
                outputDirectory=outputDirectory,injectWorkerFailure=true), ...
                "WaveVortexBenchmark:CompiledPreviewWorkers");
            testCase.verifyEqual(string(pwd),string(originalDirectory));
            testCase.verifyEqual(string(path),string(originalPath));
            testCase.verifyEqual(rng,originalRng);
            artifact = jsondecode(fileread(fullfile(outputDirectory,"compiled-preview-benchmark.json")));
            testCase.verifyEqual(string(artifact.status),"failed");
            testCase.verifyEqual(string(artifact.failure.stage),"workers");
            testCase.verifyNotEqual(string(artifact.failure.identifier),"");
            testCase.verifyTrue(isfile(fullfile(outputDirectory,"summary.md")));
        end
    end
end

function value = isCanonicalNativePlatform
value = ismac && computer("arch") == "maca64" && startsWith(string(version("-release")),"2026a",IgnoreCase=true);
end

function value = comparisonFixture
record = struct("status","complete","compiledSpeedup",1.5,"maximumRelativeError",1e-15,"libraryIdentityPassed",true,"nativeExecutionPassed",true,"noFallback",true,"lifecyclePassed",true,"planCount",17,"persistentFullHermitianBytes",0,"exactRetainedRatio",2,"operationPeakRSSRatio",2);
value = repmat(record,2,1);
end

function raw = rawFixture
definition = struct("id","constant-hydrostatic-16x12x9","Nxyz",[16 12 9],"isHydrostatic",true,"shouldAntialias",true,"seed",17,"warmupCount",0,"sampleCount",2);
matlabRun = runFixture("matlab",definition,[0.8 1.0],120,1000,1050);
compiledRun = runFixture("compiled",definition,[0.4 0.5],240,1200,1300);
comparison = struct("id",definition.id,"matlabSeconds",0.9,"compiledSeconds",0.45,"maximumRelativeError",1e-15,"matlabExactRetainedBytes",120,"compiledExactRetainedBytes",240);
provider = struct("provider",struct("id","native-neon-pthreads","version","3.3.11","threadBackend","pthreads"),"compiler",struct("mexVersion","21.0.0"));
environment = struct("processor","Apple M5 Max","physicalMemoryBytes",64*2^30,"os","macOS","matlabVersion","26.1.0 (R2026a) Update 4","release","2026a","architecture","maca64","requestedThreads",18);
source = struct("commit",repmat('a',1,40),"isDirty",false);
configuration = struct("correctnessTolerance",1e-12,"exactScope","active application-owned arrays");
raw = struct("schemaVersion","compiled-preview-benchmark-v1","status","complete","runId","20260812T120000000Z","source",source,"environment",environment,"configuration",configuration,"provider",provider,"cases",definition,"runs",[matlabRun; compiledRun],"comparison",comparison);
end

function value = runFixture(implementation,definition,samples,exactBytes,steadyBytes,peakBytes)
rss = struct("steadyRetainedBytes",steadyBytes,"operationPeakBytes",peakBytes,"operationPeakIncrementBytes",peakBytes-steadyBytes);
value = struct("implementation",implementation,"case",definition,"rawSeconds",samples,"rss",rss,"ledger",struct("exactRetainedApplicationBytes",exactBytes));
end
