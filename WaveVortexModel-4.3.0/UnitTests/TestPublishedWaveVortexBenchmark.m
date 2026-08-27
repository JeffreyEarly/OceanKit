classdef TestPublishedWaveVortexBenchmark < matlab.unittest.TestCase
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
        function schemasAndCatalogDefineThePublishedBoundary(testCase)
            publishedSchema = jsondecode(fileread(fullfile(testCase.benchmarkFolder,"schemas","published-benchmark-v1.schema.json")));
            catalogSchema = jsondecode(fileread(fullfile(testCase.benchmarkFolder,"schemas","benchmark-catalog-v1.schema.json")));
            catalog = jsondecode(fileread(fullfile(testCase.benchmarkFolder,"results","catalog.json")));

            testCase.verifyEqual(string(publishedSchema.title),"WaveVortexModel published benchmark v1");
            testCase.verifyEqual(string(catalogSchema.title),"WaveVortexModel benchmark catalog v1");
            testCase.verifyEqual(string(catalog.schemaVersion),"benchmark-catalog-v1");
            testCase.verifyEqual(string({catalog.scoringReferences.suiteId}),["core-v1" "scaling-standard-v1" "scaling-large-v1"]);
            testCase.verifyEqual(string({catalog.scoringReferences.backendId}),repmat("builtin",1,3));
            testCase.verifyNumElements(catalog.publishedDatasets,7);
            for iReference = 1:numel(catalog.scoringReferences)
                reference = catalog.scoringReferences(iReference);
                relativePath = string(reference.rawArtifact);
                testCase.verifyTrue(isfile(fullfile(testCase.repositoryRoot,relativePath)));
                testCase.verifyTrue(startsWith(relativePath,"Benchmarks/results/reference/"));
                testCase.verifyFalse(contains(relativePath,["transform-layout" "retirement" "storage"]));
                testCase.verifyFalse(startsWith(relativePath,["/" "\\"]));
                testCase.verifyFalse(contains(relativePath,[".." "\\"]));
            end
            for iDataset = 1:numel(catalog.publishedDatasets)
                entry = catalog.publishedDatasets(iDataset);
                relativePath = string(entry.artifact);
                testCase.verifyTrue(isfile(fullfile(testCase.repositoryRoot,relativePath)));
                testCase.verifyTrue(startsWith(relativePath,"Benchmarks/results/published/"));
                testCase.verifyFalse(startsWith(relativePath,["/" "\\"]));
                testCase.verifyFalse(contains(relativePath,[".." "\\"]));
            end
        end

        function initialPlatformDatasetsMeetPublicationContract(testCase)
            catalog = jsondecode(fileread(fullfile(testCase.benchmarkFolder,"results","catalog.json")));
            standardPlatforms = strings(0,1);
            toolchainVersions = strings(0,1);
            for iDataset = 1:numel(catalog.publishedDatasets)
                entry = catalog.publishedDatasets(iDataset);
                dataset = jsondecode(fileread(fullfile(testCase.repositoryRoot,string(entry.artifact))));
                testCase.verifyEqual(string(dataset.schemaVersion),"published-benchmark-v1");
                testCase.verifyEqual(string(dataset.datasetId),string(entry.datasetId));
                if string(dataset.benchmark.suiteId) == "core-v1"
                    testCase.verifyEqual(string(dataset.implementation.version),"unreleased-preview");
                    testCase.verifyEqual(string(dataset.implementation.commit),"3b762e518bf5ef92681bab7ac48dfe54b10fc708");
                else
                    testCase.verifyEqual(string(dataset.implementation.version),"4.2.1");
                    testCase.verifyEqual(string(dataset.implementation.commit),"9652b116b3ffd4ee3372cc5cdeea9700cd6cbc32");
                end
                testCase.verifyFalse(logical(dataset.implementation.sourceDirty));
                testCase.verifySubstring(string(dataset.platform.processor),"Apple M");
                testCase.verifyGreaterThan(double(dataset.platform.physicalMemoryBytes),0);
                testCase.verifyGreaterThan(double(dataset.platform.threadCount),0);
                testCase.verifyTrue(isfile(fullfile(testCase.repositoryRoot,string(dataset.provenance.rawArtifact))));

                cases = dataset.cases;
                if ~iscell(cases)
                    cases = num2cell(cases);
                end
                statuses = string(cellfun(@(benchmarkCase)benchmarkCase.status,cases,"UniformOutput",false));
                if string(dataset.benchmark.suiteId) == "core-v1"
                    testCase.verifyNumElements(cases,4);
                    testCase.verifyEqual(statuses,repmat("complete",size(statuses)));
                elseif string(dataset.benchmark.suiteId) == "scaling-standard-v1"
                    testCase.verifyNumElements(cases,34);
                    testCase.verifyEqual(statuses,repmat("complete",size(statuses)));
                    standardPlatforms(end+1,1) = string(dataset.platform.id);
                else
                    testCase.verifyEqual(string(dataset.benchmark.suiteId),"scaling-large-v1");
                    testCase.verifyNumElements(cases,19);
                    if string(dataset.platform.id) == "m5-max"
                        testCase.verifyEqual(nnz(statuses == "complete"),16);
                        testCase.verifyEqual(nnz(statuses == "unavailable"),3);
                        reasons = string(cellfun(@(benchmarkCase)unavailableReason(benchmarkCase),cases(statuses == "unavailable"),"UniformOutput",false));
                        testCase.verifyTrue(all(contains(reasons,"48 GiB")));
                    else
                        testCase.verifyEqual(statuses,repmat("complete",size(statuses)));
                    end
                end
                for iCase = find(statuses == "complete")
                    benchmarkCase = cases{iCase};
                    testCase.verifyTrue(logical(benchmarkCase.correctness.passed));
                    testCase.verifyTrue(all(isfinite(benchmarkCase.timing.samplesSeconds)));
                    testCase.verifyGreaterThan(double(benchmarkCase.timing.medianSeconds),0);
                    testCase.verifyEqual(string(benchmarkCase.memory.status),"complete");
                end
                toolchainVersions(end+1,1) = string(dataset.toolchain.version);
            end
            testCase.verifyEqual(sort(standardPlatforms),sort(["lyra"; "m5-max"; "matilda"]));
            testCase.verifyTrue(any(contains(toolchainVersions,"R2025b") & contains(toolchainVersions,"Update 4")));
            testCase.verifyTrue(any(contains(toolchainVersions,"R2026a") & contains(toolchainVersions,"Update 4")));
        end

        function legacyArtifactNormalizesWithoutMutation(testCase)
            rawRelativePath = "Benchmarks/results/reference/scaling-standard-v1-m5-max-r2026a-builtin/benchmark.json";
            rawPath = fullfile(testCase.repositoryRoot,rawRelativePath);
            originalBytes = fileread(rawPath);
            dataset = publishedWaveVortexBenchmarkFromMatlabArtifact(rawPath,suiteId="scaling-standard-v1",platformId="m5-max",platformName="Apple M5 Max",provenancePath=rawRelativePath,processorName="Apple M5 Max",implementationVersion="4.2.1");

            expectedFields = ["schemaVersion" "datasetId" "collectedAt" "benchmark" "implementation" "platform" "toolchain" "provenance" "cases"];
            testCase.verifyEqual(string(fieldnames(dataset))',expectedFields);
            testCase.verifyEqual(dataset.datasetId,"scaling-standard-v1--matlab-builtin--m5-max--20260808T180708Z");
            testCase.verifyEqual(dataset.benchmark.operation,"nonlinearFlux");
            testCase.verifyEqual(dataset.implementation.version,"4.2.1");
            testCase.verifyEqual(dataset.platform.processor,"Apple M5 Max");
            testCase.verifyNumElements(dataset.cases,34);
            testCase.verifyEqual(fileread(rawPath),originalBytes);

            first = dataset.cases{1};
            testCase.verifyEqual(first.status,"complete");
            testCase.verifyEqual(first.timing.medianSeconds,median(first.timing.samplesSeconds),AbsTol=1e-15);
            testCase.verifyEqual(first.memory.peakIncrementBytes,first.memory.peakProcessBytes-first.memory.baselineProcessBytes);
        end

        function futureArtifactUsesAddedMetadataAndMarksMissingBackend(testCase)
            rawPath = fullfile(testCase.repositoryRoot,"Benchmarks","results","reference","scaling-standard-v1-m5-max-r2026a-builtin","benchmark.json");
            raw = jsondecode(fileread(rawPath));
            raw.schemaVersion = "1.1.0";
            raw.environment.packageName = "WaveVortexModel";
            raw.environment.packageVersion = "4.2.1";
            raw.environment.processorName = "Future Reference CPU";
            raw.suites.cases(1).backends = struct([]);
            fixturePath = fullfile(testCase.temporaryFolder,"raw.json");
            writeJson(fixturePath,raw);

            dataset = publishedWaveVortexBenchmarkFromMatlabArtifact(fixturePath,suiteId="scaling-standard-v1",platformId="future-host",platformName="Future host",provenancePath="Benchmarks/results/runs/future/benchmark.json");
            testCase.verifyEqual(dataset.implementation.displayName,"WaveVortexModel MATLAB");
            testCase.verifyEqual(dataset.implementation.version,"4.2.1");
            testCase.verifyEqual(dataset.platform.processor,"Future Reference CPU");
            testCase.verifyEqual(dataset.provenance.rawSchemaVersion,"1.1.0");
            testCase.verifyEqual(dataset.cases{1}.status,"unavailable");
            testCase.verifySubstring(dataset.cases{1}.unavailableReason,"builtin");
            testCase.verifyEqual(dataset.cases{2}.status,"complete");
        end

        function cppDocumentUsesTheSameLanguageNeutralVocabulary(testCase)
            cpp = exampleCppDataset;
            encoded = jsonencode(cpp);
            decoded = jsondecode(encoded);
            expectedFields = ["schemaVersion" "datasetId" "collectedAt" "benchmark" "implementation" "platform" "toolchain" "provenance" "cases"];

            testCase.verifyEqual(string(fieldnames(decoded))',expectedFields);
            testCase.verifyEqual(string(decoded.schemaVersion),"published-benchmark-v1");
            testCase.verifyEqual(string(decoded.implementation.id),"cpp");
            testCase.verifyEqual(string(decoded.toolchain.kind),"cpp");
            testCase.verifyFalse(isfield(decoded.toolchain.details,"matlabRelease"));
            testCase.verifyEqual(string(decoded.cases.status),"unavailable");
        end

        function runnerContainsNoMachineSpecificReferenceSelection(testCase)
            source = lower(string(fileread(fullfile(testCase.benchmarkFolder,"runWaveVortexBenchmark.m"))));
            testCase.verifySubstring(source,"scoringreferencefromcatalog");
            testCase.verifyFalse(contains(source,"m5-max"));
            testCase.verifyFalse(contains(source,"r2026a"));
            testCase.verifyFalse(contains(source,"referenceartifactpath"));
        end
    end
end

function dataset = exampleCppDataset
benchmark = struct("suiteId","scaling-standard-v1","suiteVersion",1,"operation","nonlinearFlux","correctnessTolerance",1e-12);
implementation = struct("id","cpp","displayName","WaveVortexModel C++","version","0.1.0","repository","https://github.com/JeffreyEarly/wave-vortex-model-cpp","commit",repmat('a',1,40),"backend","fftw","sourceDirty",false);
platform = struct("id","example-host","displayName","Example host","processor","Example CPU","physicalMemoryBytes",16*2^30,"os","Example OS","architecture","example64","threadCount",8);
toolchain = struct("kind","cpp","name","Clang","version","18.0.0","details",struct("buildType","Release"));
provenance = struct("rawArtifact","Benchmarks/results/cpp/example.json","rawSchemaVersion","cpp-v1");
configuration = struct("Lxyz",[15e3 15e3 1300],"Nxyz",[64 64 65],"isHydrostatic",false,"shouldAntialias",true,"seed",1,"warmupCount",2,"sampleCount",3);
benchmarkCase = struct("id","boussinesq-64x64x65","transformId","boussinesq","scoreFamily","boussinesq","configuration",configuration,"status","unavailable","unavailableReason","Transform not implemented by this backend.");
dataset = struct("schemaVersion","published-benchmark-v1","datasetId","scaling-standard-v1--cpp-fftw--example-host--20260810T123456Z","collectedAt","2026-08-10T12:34:56Z","benchmark",benchmark,"implementation",implementation,"platform",platform,"toolchain",toolchain,"provenance",provenance,"cases",{{benchmarkCase}});
end

function writeJson(pathname,value)
fileId = fopen(pathname,"w");
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s\n",jsonencode(value,PrettyPrint=true));
clear cleanup
end

function reason = unavailableReason(benchmarkCase)
reason = "";
if isfield(benchmarkCase,"unavailableReason")
    reason = string(benchmarkCase.unavailableReason);
end
end
