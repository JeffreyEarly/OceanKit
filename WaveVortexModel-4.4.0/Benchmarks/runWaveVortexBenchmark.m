function results = runWaveVortexBenchmark(options)
% Run versioned WaveVortex performance and memory benchmark suites.
%
% This authoring benchmark measures state-advanced nonlinear-advection
% evaluations while retaining production caches. Ordinary runs are written
% beneath Benchmarks/results/runs. Reference generation is explicit.
arguments
    options.suites (1,:) string = "core-v1"
    options.backends (1,:) string = "builtin"
    options.caseIds (1,:) string = strings(1,0)
    options.outputDirectory (1,1) string = ""
    options.referenceDirectory (1,1) string = ""
    options.catalogPath (1,1) string = ""
    options.shouldMeasureMemory (1,1) logical = true
    options.shouldWriteArtifacts (1,1) logical = true
    options.shouldCreateReference (1,1) logical = false
    options.correctnessTolerance (1,1) double {mustBePositive} = 1e-12
    options.runId (1,1) string = ""
end

benchmarkFolder = fileparts(mfilename("fullpath"));
repositoryRoot = fileparts(benchmarkFolder);
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
addRepositoryPaths(repositoryRoot,benchmarkFolder);

suites = waveVortexBenchmarkSuites(options.suites);
if ~isempty(options.caseIds)
    suites = filterSuiteCases(suites,options.caseIds);
end
backends = waveVortexBenchmarkBackends(options.backends);
if options.runId == ""
    options.runId = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmss'Z'"));
end
if options.catalogPath == ""
    options.catalogPath = fullfile(benchmarkFolder,"results","catalog.json");
end
if options.shouldCreateReference && options.referenceDirectory == ""
    error("WaveVortexBenchmark:ReferenceDirectoryRequired","referenceDirectory is required when shouldCreateReference is true.");
elseif ~options.shouldCreateReference && options.referenceDirectory ~= ""
    error("WaveVortexBenchmark:ReferenceDirectoryNotApplicable","referenceDirectory is used only when shouldCreateReference is true.");
end
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId + "-" + computer("arch") + "-" + version("-release"));
end

results = struct("schemaVersion","1.1.0","status","complete","runId",options.runId,"environment",benchmarkEnvironment(repositoryRoot),"configuration",struct("suiteIds",options.suites,"backendIds",options.backends,"caseIds",options.caseIds,"correctnessTolerance",options.correctnessTolerance,"shouldMeasureMemory",options.shouldMeasureMemory),"suites",emptySuiteResults());
for iSuite = 1:numel(suites)
    suiteResult = runSuite(suites(iSuite),backends,options,benchmarkFolder,repositoryRoot);
    results.suites(end+1) = suiteResult;
    if suiteResult.status ~= "complete"
        results.status = "partial";
    end
end

if options.shouldWriteArtifacts
    writeRunArtifacts(results,options.outputDirectory);
end
clear stateCleanup
end

function suiteResult = runSuite(suite,backends,options,benchmarkFolder,repositoryRoot)
[referenceFile,referenceArtifact,referenceOutputDirectory] = referenceLocations(suite,options,repositoryRoot);
if suite.kind == "transform-layout"
    suiteResult = runWaveVortexTransformLayoutSuite(suite,options.correctnessTolerance,repositoryRoot);
    suiteResult.referenceArtifact = referenceArtifact;
    if options.shouldCreateReference
        referenceResults = struct("schemaVersion","1.1.0","status",suiteResult.status,"runId",options.runId,"environment",benchmarkEnvironment(repositoryRoot),"configuration",struct("suiteIds",suite.id,"backendIds","builtin","correctnessTolerance",options.correctnessTolerance,"shouldMeasureMemory",false),"suites",suiteResult);
        writeRunArtifacts(referenceResults,referenceOutputDirectory);
    end
    return
end
suiteResult = struct("id",suite.id,"version",suite.version,"kind",suite.kind,"description",suite.description,"operation",suite.operation,"isScored",suite.isScored,"selectionIsComplete",suite.selectionIsComplete,"status","complete","cases",emptyCaseResults(),"familyScores",emptyScores(),"suiteScores",emptyScores(),"referenceArtifact","","metadata",struct);
for iCase = 1:numel(suite.cases)
    try
        caseResult = runCase(suite.cases(iCase),backends,options,benchmarkFolder,repositoryRoot);
    catch exception
        caseResult = failedCaseResult(suite.cases(iCase),exception);
        suiteResult.status = "partial";
    end
    suiteResult.cases(end+1) = caseResult;
end
if ~suite.selectionIsComplete
    suiteResult.status = "partial";
end

suiteResult.referenceArtifact = referenceArtifact;
if options.shouldCreateReference
    suiteResult = applySelfReferenceScores(suiteResult);
else
    suiteResult = applyReferenceScores(suiteResult,referenceFile);
end
suiteResult = calculateAggregateScores(suiteResult);

if options.shouldCreateReference
    referenceResults = struct("schemaVersion","1.1.0","status",suiteResult.status,"runId",options.runId,"environment",benchmarkEnvironment(repositoryRoot),"configuration",struct("suiteIds",suite.id,"backendIds","builtin","correctnessTolerance",options.correctnessTolerance,"shouldMeasureMemory",options.shouldMeasureMemory),"suites",suiteResult);
    writeRunArtifacts(referenceResults,referenceOutputDirectory);
end
end

function [referenceFile,referenceArtifact,referenceOutputDirectory] = referenceLocations(suite,options,repositoryRoot)
referenceFile = "";
referenceArtifact = "";
referenceOutputDirectory = "";
if options.shouldCreateReference
    referenceOutputDirectory = fullfile(options.referenceDirectory,suite.id);
    referenceFile = fullfile(referenceOutputDirectory,"benchmark.json");
    referenceArtifact = referenceOutputDirectory;
elseif suite.isScored
    [referenceFile,referenceArtifact] = scoringReferenceFromCatalog(options.catalogPath,suite.id,repositoryRoot);
end
end

function [absolutePath,relativePath] = scoringReferenceFromCatalog(catalogPath,suiteId,repositoryRoot)
if ~isfile(catalogPath)
    error("WaveVortexBenchmark:MissingCatalog","Benchmark catalog does not exist: %s.",catalogPath);
end
catalog = jsondecode(fileread(catalogPath));
if ~isfield(catalog,"schemaVersion") || string(catalog.schemaVersion) ~= "benchmark-catalog-v1" || ~isfield(catalog,"scoringReferences")
    error("WaveVortexBenchmark:InvalidCatalog","Benchmark catalog must use schema benchmark-catalog-v1 and contain scoringReferences.");
end
references = catalog.scoringReferences;
matches = string({references.suiteId}) == suiteId;
if nnz(matches) ~= 1
    error("WaveVortexBenchmark:MissingScoringReference","Benchmark catalog must contain exactly one scoring reference for suite %s.",suiteId);
end
reference = references(matches);
if string(reference.backendId) ~= "builtin"
    error("WaveVortexBenchmark:InvalidCatalog","Scoring reference for suite %s must use the builtin backend.",suiteId);
end
relativePath = string(reference.rawArtifact);
parts = split(relativePath,"/");
if startsWith(relativePath,["/" "\\"]) || ~isempty(regexp(relativePath,"^[A-Za-z]:","once")) || contains(relativePath,"\\") || any(parts == "..") || any(parts == "")
    error("WaveVortexBenchmark:InvalidCatalog","Scoring reference paths must be repository-relative and cannot contain traversal.");
end
absolutePath = fullfile(repositoryRoot,relativePath);
if ~isfile(absolutePath)
    error("WaveVortexBenchmark:MissingScoringReference","Scoring reference does not exist: %s.",relativePath);
end
end

function caseResult = runCase(benchmarkCase,backends,options,benchmarkFolder,repositoryRoot)
backendResults = emptyBackendResults();
models = cell(1,numel(backends));
states = cell(1,numel(backends));
for iBackend = 1:numel(backends)
    constructionTimer = tic;
    models{iBackend} = createWaveVortexBenchmarkTransform(benchmarkCase,backends(iBackend).id);
    constructionSeconds = toc(constructionTimer);
    states{iBackend} = initializeWaveVortexBenchmarkState(models{iBackend},benchmarkCase.seed);
    advanceWaveVortexBenchmarkState(models{iBackend},states{iBackend},0);
    firstTimer = tic;
    executeWaveVortexBenchmarkOperation(models{iBackend},benchmarkCase.operation);
    firstCallSeconds = toc(firstTimer);
    backendResults(end+1) = baseBackendResult(backends(iBackend).id,constructionSeconds,firstCallSeconds,benchmarkCase.sampleCount); %#ok<AGROW>
end

for iWarmup = 1:benchmarkCase.warmupCount
    for iBackend = 1:numel(backends)
        advanceWaveVortexBenchmarkState(models{iBackend},states{iBackend},iWarmup);
        executeWaveVortexBenchmarkOperation(models{iBackend},benchmarkCase.operation);
    end
end

for iSample = 1:benchmarkCase.sampleCount
    sampleOutputs = cell(1,numel(backends));
    startIndex = mod(iSample-1,numel(backends))+1;
    executionOrder = [startIndex:numel(backends),1:startIndex-1];
    for iOrder = executionOrder
        advanceWaveVortexBenchmarkState(models{iOrder},states{iOrder},benchmarkCase.warmupCount+iSample);
        sampleTimer = tic;
        sampleOutputs{iOrder} = executeWaveVortexBenchmarkOperation(models{iOrder},benchmarkCase.operation);
        backendResults(iOrder).rawSeconds(iSample) = toc(sampleTimer);
    end
    builtinIndex = find(string({backends.id}) == "builtin",1);
    if ~isempty(builtinIndex)
        for iBackend = 1:numel(backends)
            backendResults(iBackend).relativeError = max(backendResults(iBackend).relativeError,waveVortexBenchmarkRelativeError(sampleOutputs{builtinIndex},sampleOutputs{iBackend}));
        end
    end
end

for iBackend = 1:numel(backends)
    backendResults(iBackend).medianSeconds = median(backendResults(iBackend).rawSeconds);
    cacheTimer = tic;
    executeWaveVortexBenchmarkOperation(models{iBackend},benchmarkCase.operation);
    backendResults(iBackend).sameStateCacheHitSeconds = toc(cacheTimer);
    backendResults(iBackend).correctnessPassed = backendResults(iBackend).relativeError <= options.correctnessTolerance;
    if options.shouldMeasureMemory
        backendResults(iBackend).memory = measureCaseMemory(benchmarkCase,backends(iBackend).id,benchmarkFolder,repositoryRoot);
    end
end
caseResult = struct("id",benchmarkCase.id,"transformId",benchmarkCase.transformId,"scoreFamily",benchmarkCase.scoreFamily,"operation",benchmarkCase.operation,"Lxyz",benchmarkCase.Lxyz,"Nxyz",benchmarkCase.Nxyz,"isHydrostatic",benchmarkCase.isHydrostatic,"shouldAntialias",benchmarkCase.shouldAntialias,"seed",benchmarkCase.seed,"warmupCount",benchmarkCase.warmupCount,"sampleCount",benchmarkCase.sampleCount,"status","complete","failure",emptyFailure(),"backends",backendResults);
clear models
end

function result = baseBackendResult(backendId,constructionSeconds,firstCallSeconds,sampleCount)
result = struct("id",backendId,"status","complete","constructionSeconds",constructionSeconds,"firstCallSeconds",firstCallSeconds,"sameStateCacheHitSeconds",NaN,"rawSeconds",NaN(1,sampleCount),"medianSeconds",NaN,"relativeError",0,"correctnessPassed",false,"referenceMedianSeconds",NaN,"caseScore",NaN,"sameHostSpeedup",NaN,"memory",emptyMemory(),"failure",emptyFailure());
end

function caseResult = failedCaseResult(benchmarkCase,exception)
caseResult = struct("id",benchmarkCase.id,"transformId",benchmarkCase.transformId,"scoreFamily",benchmarkCase.scoreFamily,"operation",benchmarkCase.operation,"Lxyz",benchmarkCase.Lxyz,"Nxyz",benchmarkCase.Nxyz,"isHydrostatic",benchmarkCase.isHydrostatic,"shouldAntialias",benchmarkCase.shouldAntialias,"seed",benchmarkCase.seed,"warmupCount",benchmarkCase.warmupCount,"sampleCount",benchmarkCase.sampleCount,"status","failed","failure",exceptionFailure(exception),"backends",emptyBackendResults());
end

function suiteResult = applySelfReferenceScores(suiteResult)
for iCase = 1:numel(suiteResult.cases)
    for iBackend = 1:numel(suiteResult.cases(iCase).backends)
        backend = suiteResult.cases(iCase).backends(iBackend);
        if backend.id == "builtin" && backend.status == "complete"
            suiteResult.cases(iCase).backends(iBackend).referenceMedianSeconds = backend.medianSeconds;
            suiteResult.cases(iCase).backends(iBackend).caseScore = 100;
            suiteResult.cases(iCase).backends(iBackend).sameHostSpeedup = 1;
        end
    end
end
end

function suiteResult = applyReferenceScores(suiteResult,referenceFile)
if referenceFile == "" || ~isfile(referenceFile)
    return
end
reference = jsondecode(fileread(referenceFile));
referenceSuite = reference.suites(1);
for iCase = 1:numel(suiteResult.cases)
    referenceCaseIndex = find(string({referenceSuite.cases.id}) == suiteResult.cases(iCase).id,1);
    if isempty(referenceCaseIndex)
        continue
    end
    referenceBuiltinIndex = find(string({referenceSuite.cases(referenceCaseIndex).backends.id}) == "builtin",1);
    currentBuiltinIndex = find(string({suiteResult.cases(iCase).backends.id}) == "builtin",1);
    if isempty(referenceBuiltinIndex)
        continue
    end
    referenceMedian = referenceSuite.cases(referenceCaseIndex).backends(referenceBuiltinIndex).medianSeconds;
    currentBuiltinMedian = NaN;
    if ~isempty(currentBuiltinIndex)
        currentBuiltinMedian = suiteResult.cases(iCase).backends(currentBuiltinIndex).medianSeconds;
    end
    for iBackend = 1:numel(suiteResult.cases(iCase).backends)
        medianSeconds = suiteResult.cases(iCase).backends(iBackend).medianSeconds;
        suiteResult.cases(iCase).backends(iBackend).referenceMedianSeconds = referenceMedian;
        suiteResult.cases(iCase).backends(iBackend).caseScore = 100*referenceMedian/medianSeconds;
        suiteResult.cases(iCase).backends(iBackend).sameHostSpeedup = currentBuiltinMedian/medianSeconds;
    end
end
end

function suiteResult = calculateAggregateScores(suiteResult)
if ~suiteResult.isScored || ~suiteResult.selectionIsComplete
    return
end
backendIds = strings(1,0);
for iCase = 1:numel(suiteResult.cases)
    backendIds = union(backendIds,string({suiteResult.cases(iCase).backends.id}),"stable");
end
families = unique(string({suiteResult.cases.scoreFamily}),"stable");
for backendId = backendIds
    familyValues = NaN(1,numel(families));
    for iFamily = 1:numel(families)
        caseMask = string({suiteResult.cases.scoreFamily}) == families(iFamily);
        scores = NaN(1,nnz(caseMask));
        selectedCases = suiteResult.cases(caseMask);
        for iCase = 1:numel(selectedCases)
            backendIndex = find(string({selectedCases(iCase).backends.id}) == backendId,1);
            if ~isempty(backendIndex)
                scores(iCase) = selectedCases(iCase).backends(backendIndex).caseScore;
            end
        end
        if all(isfinite(scores))
            familyValues(iFamily) = exp(mean(log(scores)));
            suiteResult.familyScores(end+1) = struct("id",families(iFamily),"backendId",backendId,"score",familyValues(iFamily));
        end
    end
    if all(isfinite(familyValues))
        suiteResult.suiteScores(end+1) = struct("id",suiteResult.id,"backendId",backendId,"score",exp(mean(log(familyValues))));
    end
end
end

function memory = measureCaseMemory(benchmarkCase,backendId,benchmarkFolder,repositoryRoot)
memory = emptyMemory();
configPath = string(tempname) + ".json";
outputPath = string(tempname) + ".json";
config = struct("benchmarkCase",benchmarkCase,"backendId",backendId,"repositoryRoot",repositoryRoot,"benchmarkFolder",benchmarkFolder,"matlabPath",path);
writeText(configPath,jsonencode(config));
cleanup = onCleanup(@()deleteTemporaryFiles(configPath,outputPath));
matlabExecutable = fullfile(matlabroot,"bin","matlab");
statement = "addpath('" + replace(benchmarkFolder,"'","''") + "'); waveVortexBenchmarkMemoryWorker('" + replace(configPath,"'","''") + "','" + replace(outputPath,"'","''") + "')";
command = sprintf('"%s" -batch "%s"',matlabExecutable,replace(statement,'"','\"'));
[exitCode,commandOutput] = system(command);
if exitCode ~= 0 || ~isfile(outputPath)
    memory.status = "failed";
    memory.failure = struct("identifier","WaveVortexBenchmark:MemoryWorkerFailed","message",string(commandOutput),"stack",strings(0,1));
    return
end
memory = jsondecode(fileread(outputPath));
clear cleanup
end

function writeRunArtifacts(results,outputDirectory)
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writeText(fullfile(outputDirectory,"benchmark.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(outputDirectory,"summary.md"),benchmarkSummary(results));
end

function summary = benchmarkSummary(results)
summary = waveVortexBenchmarkSummary(results);
end

function environment = benchmarkEnvironment(repositoryRoot)
[~,commit] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
[~,dirty] = system(sprintf('git -C "%s" status --porcelain',repositoryRoot));
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
environment = struct("os",string(system_dependent("getos")),"processor",string(system_dependent("getcpu")),"processorName",processorName(),"physicalMemoryBytes",physicalMemoryBytes(),"matlabVersion",string(version),"matlabRelease",string(version("-release")),"architecture",string(computer("arch")),"requestedThreads",maxNumCompThreads,"sourceCommit",strtrim(string(commit)),"sourceDirty",strlength(strtrim(string(dirty))) > 0,"packageName",string(metadata.name),"packageVersion",string(metadata.version));
end

function name = processorName()
name = "";
if ismac
    [status,output] = system("sysctl -n machdep.cpu.brand_string");
    if status == 0
        name = strtrim(string(output));
    end
elseif isunix && isfile("/proc/cpuinfo")
    text = fileread("/proc/cpuinfo");
    token = regexp(text,"(?:model name|Hardware)\s*:\s*([^\r\n]+)","tokens","once");
    if ~isempty(token)
        name = strtrim(string(token{1}));
    end
elseif ispc
    name = strtrim(string(getenv("PROCESSOR_IDENTIFIER")));
end
if name == ""
    name = string(system_dependent("getcpu"));
end
end

function bytes = physicalMemoryBytes()
bytes = NaN;
if ispc
    [~,systemMemory] = memory;
    bytes = double(systemMemory.PhysicalMemory.Total);
elseif ismac
    [status,output] = system("sysctl -n hw.memsize");
    if status == 0
        bytes = str2double(strtrim(output));
    end
elseif isunix && isfile("/proc/meminfo")
    text = fileread("/proc/meminfo");
    token = regexp(text,"MemTotal:\s+(\d+)\s+kB","tokens","once");
    if ~isempty(token)
        bytes = 1024*str2double(token{1});
    end
end
end

function addRepositoryPaths(repositoryRoot,benchmarkFolder)
addpath(repositoryRoot,benchmarkFolder);
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
for iFolder = 1:numel(metadata.folders)
    folder = fullfile(repositoryRoot,metadata.folders(iFolder).path);
    if isfolder(folder)
        addpath(folder);
    end
end
end

function suites = filterSuiteCases(suites,caseIds)
foundIds = strings(1,0);
for iSuite = 1:numel(suites)
    suiteCaseIds = string({suites(iSuite).cases.id});
    keep = ismember(suiteCaseIds,caseIds);
    if nnz(keep) ~= numel(suiteCaseIds)
        suites(iSuite).selectionIsComplete = false;
    end
    suites(iSuite).cases = suites(iSuite).cases(keep);
    foundIds = union(foundIds,suiteCaseIds(keep),"stable");
end
unknownIds = setdiff(caseIds,foundIds);
if ~isempty(unknownIds)
    error("WaveVortexBenchmark:UnknownCase","Unknown benchmark case: %s.",strjoin(unknownIds,", "));
end
suites = suites(arrayfun(@(suite)~isempty(suite.cases),suites));
end

function restoreState(originalDirectory,originalPath,originalRng)
cd(originalDirectory);
path(originalPath);
rng(originalRng);
end

function writeText(pathname,contents)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to open %s for writing.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",contents);
clear cleanup
end

function deleteTemporaryFiles(varargin)
for iFile = 1:numel(varargin)
    if isfile(varargin{iFile})
        delete(varargin{iFile});
    end
end
end

function failure = exceptionFailure(exception)
failure = struct("identifier",string(exception.identifier),"message",string(exception.message),"stack",string({exception.stack.name}));
end

function failure = emptyFailure()
failure = struct("identifier","","message","","stack",strings(0,1));
end

function memory = emptyMemory()
memory = struct("status","not-requested","provider","","baselineBytes",NaN,"persistentBytes",NaN,"peakBytes",NaN,"persistentIncrementBytes",NaN,"peakIncrementBytes",NaN,"samplingIntervalSeconds",NaN,"failure",emptyFailure());
end

function results = emptySuiteResults()
results = struct("id",{},"version",{},"kind",{},"description",{},"operation",{},"isScored",{},"selectionIsComplete",{},"status",{},"cases",{},"familyScores",{},"suiteScores",{},"referenceArtifact",{},"metadata",{});
end

function results = emptyCaseResults()
results = struct("id",{},"transformId",{},"scoreFamily",{},"operation",{},"Lxyz",{},"Nxyz",{},"isHydrostatic",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"sampleCount",{},"status",{},"failure",{},"backends",{});
end

function results = emptyBackendResults()
results = struct("id",{},"status",{},"constructionSeconds",{},"firstCallSeconds",{},"sameStateCacheHitSeconds",{},"rawSeconds",{},"medianSeconds",{},"relativeError",{},"correctnessPassed",{},"referenceMedianSeconds",{},"caseScore",{},"sameHostSpeedup",{},"memory",{},"failure",{});
end

function results = emptyScores()
results = struct("id",{},"backendId",{},"score",{});
end
