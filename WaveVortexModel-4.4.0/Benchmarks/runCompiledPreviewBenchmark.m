function results = runCompiledPreviewBenchmark(options)
% Benchmark the exact public MATLAB and compiled nonlinear-flux entry points.
arguments
    options.sizes (:,3) double {mustBeInteger,mustBePositive} = [256 256 65; 512 512 129]
    options.hydrostatic (1,:) logical = [true false]
    options.processRunCount (1,1) double {mustBeInteger,mustBePositive} = 3
    options.warmupCount (1,1) double {mustBeInteger,mustBeNonnegative} = 2
    options.mediumSampleCount (1,1) double {mustBeInteger,mustBePositive} = 7
    options.largeSampleCount (1,1) double {mustBeInteger,mustBePositive} = 3
    options.samplingIntervalSeconds (1,1) double {mustBePositive} = 0.005
    options.plateauSeconds (1,1) double {mustBePositive} = 0.10
    options.outputHoldSeconds (1,1) double {mustBePositive} = 0.03
    options.outputDirectory (1,1) string = ""
    options.runId (1,1) string = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"))
    options.shouldWriteArtifacts (1,1) logical = true
    options.injectWorkerFailure (1,1) logical = false
end
if ~ismac || string(computer("arch")) ~= "maca64" || ~startsWith(string(version("-release")),"2026a",IgnoreCase=true)
    error("WaveVortexBenchmark:CompiledPreviewUnsupportedPlatform","The canonical compiled-preview benchmark targets MATLAB R2026a on macOS maca64.");
end
repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
benchmarkFolder = fullfile(repositoryRoot,"Benchmarks");
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
addRepositoryPaths(repositoryRoot,benchmarkFolder);
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId+"-compiled-preview-"+computer("arch")+"-"+version("-release"));
end
if options.shouldWriteArtifacts && isfolder(options.outputDirectory)
    error("WaveVortexBenchmark:CompiledPreviewOutputExists","Output already exists: %s",options.outputDirectory);
end
if options.shouldWriteArtifacts
    mkdir(options.outputDirectory);
end
results = initializeResult(options,repositoryRoot);
activeStage = "provider";
try
    capabilities = WVCompiledBackend.capabilities();
    if ~capabilities.isAvailable
        capabilities = WVCompiledBackend.build();
    end
    validateCapabilities(capabilities);
    results.provider = capabilities;
    moduleCleanup = onCleanup(@()clearModuleIfIdle(string(capabilities.module.name)));
    checkpoint(results,options);

    activeStage = "workers";
    cases = caseDefinitions(options);
    results.cases = cases;
    implementations = ["matlab" "compiled"];
    for iRun = 1:options.processRunCount
        caseOrder = mod((0:numel(cases)-1)+(iRun-1),numel(cases))+1;
        for iCase = caseOrder
            order = mod((0:1)+(iRun+iCase-2),2)+1;
            for iImplementation = order
                implementation = implementations(iImplementation);
                fprintf("Compiled preview: %s, %s, process %d/%d.\n",implementation,cases(iCase).id,iRun,options.processRunCount);
                results.runs(end+1,1) = runWorker(implementation,iRun,cases(iCase),capabilities,options,repositoryRoot,benchmarkFolder); %#ok<AGROW>
                checkpoint(results,options);
            end
        end
    end
    failed = results.runs(string({results.runs.status}) ~= "complete");
    if ~isempty(failed)
        messages = arrayfun(@(item)string(item.implementation)+"/"+string(item.case.id)+": "+string(item.failure.identifier)+" "+string(item.failure.message),failed);
        error("WaveVortexBenchmark:CompiledPreviewWorkers","One or more compiled-preview workers failed:%s%s",newline,strjoin(messages,newline));
    end

    activeStage = "correctness";
    results.correctness = correctnessRecords(cases);
    activeStage = "aggregation";
    results.comparison = comparisonRecords(results.runs,cases,results.correctness,capabilities);
    results.decision = compiledPreviewBenchmarkDecision(results.comparison);
    results.status = "complete";
    results.completedAtUTC = utcTimestamp;
    results.failure = emptyFailure;
    writeArtifacts(results,options);
catch exception
    results.status = "failed";
    results.completedAtUTC = utcTimestamp;
    results.failure = struct("stage",activeStage,"identifier",string(exception.identifier),"message",string(exception.message),"report",string(getReport(exception,"extended","hyperlinks","off")));
    writeArtifacts(results,options);
    rethrow(exception)
end
clear moduleCleanup
clear stateCleanup
end

function results = initializeResult(options,repositoryRoot)
[commit,tree,isDirty] = gitIdentity(repositoryRoot);
results = struct( ...
    "schemaVersion","compiled-preview-benchmark-v1", ...
    "status","running", ...
    "runId",options.runId, ...
    "generatedAtUTC",utcTimestamp, ...
    "completedAtUTC","", ...
    "environment",environmentRecord, ...
    "source",struct("repository","JeffreyEarly/wave-vortex-model","commit",commit,"tree",tree,"isDirty",isDirty), ...
    "configuration",struct("suiteId","core-v1","operation","public state-advanced nonlinearFlux","processRunCount",options.processRunCount,"warmupCount",options.warmupCount,"mediumSampleCount",options.mediumSampleCount,"largeSampleCount",options.largeSampleCount,"samplingIntervalSeconds",options.samplingIntervalSeconds,"rssBaseline","steady retained state of the active public backend","exactScope","active transform, populated cache, compiled-owned storage when active, and three outputs; canonical Ap/Am/A0 excluded","speedThreshold",1.25,"correctnessTolerance",1e-12), ...
    "provider",struct(), ...
    "cases",[], ...
    "runs",repmat(emptyRun,0,1), ...
    "correctness",[], ...
    "comparison",[], ...
    "decision",struct(), ...
    "failure",emptyFailure);
end

function validateCapabilities(capabilities)
if ~capabilities.isAvailable || string(capabilities.schemaVersion) ~= "1.0.0" || string(capabilities.provider.id) ~= "native-neon-pthreads" || ~capabilities.module.identityValidated || capabilities.libraries.openmp.detected || capabilities.contract.version ~= 4 || capabilities.contract.planCount ~= 17 || capabilities.featureValidation.maximumRelativeError > 1e-12
    error("WaveVortexBenchmark:CompiledPreviewCapability","The canonical public preview requires the validated native-neon-pthreads contract with no OpenMP and numerical error at most 1e-12.");
end
end

function cases = caseDefinitions(options)
suite = waveVortexBenchmarkSuites("core-v1");
cases = repmat(struct("id","","Nxyz",[],"isHydrostatic",false,"shouldAntialias",true,"seed",0,"warmupCount",options.warmupCount,"sampleCount",0),0,1);
for iSize = 1:size(options.sizes,1)
    for isHydrostatic = options.hydrostatic
        Nxyz = options.sizes(iSize,:);
        match = find(arrayfun(@(item)isequal(item.Nxyz,Nxyz)&&item.isHydrostatic==isHydrostatic,suite.cases),1);
        if isempty(match)
            seed = 175000+sum(Nxyz)+100*isHydrostatic;
        else
            seed = suite.cases(match).seed;
        end
        sampleCount = options.mediumSampleCount;
        if iSize == size(options.sizes,1) && size(options.sizes,1) > 1
            sampleCount = options.largeSampleCount;
        end
        identifier = sprintf("constant-%s-%dx%dx%d",conditional(isHydrostatic,"hydrostatic","nonhydrostatic"),Nxyz(1),Nxyz(2),Nxyz(3));
        cases(end+1,1) = struct("id",identifier,"Nxyz",Nxyz,"isHydrostatic",isHydrostatic,"shouldAntialias",true,"seed",seed,"warmupCount",options.warmupCount,"sampleCount",sampleCount); %#ok<AGROW>
    end
end
end

function run = runWorker(implementation,repeatIndex,definition,capabilities,options,repositoryRoot,benchmarkFolder)
config = struct( ...
    "implementation",implementation, ...
    "sourceCommit",gitValue(repositoryRoot,"rev-parse HEAD"), ...
    "repeatIndex",repeatIndex, ...
    "caseDefinition",definition, ...
    "expectedModuleHash",string(capabilities.module.sha256), ...
    "repositoryRoot",repositoryRoot, ...
    "benchmarkFolder",benchmarkFolder, ...
    "matlabPath",path, ...
    "samplerPath",fullfile(benchmarkFolder,"sampleProcessRSS.sh"), ...
    "samplingIntervalSeconds",options.samplingIntervalSeconds, ...
    "plateauSeconds",options.plateauSeconds, ...
    "outputHoldSeconds",options.outputHoldSeconds);
configPath = string(tempname)+".json";
outputPath = string(tempname)+".json";
cleanup = onCleanup(@()deleteTemporaryFiles(configPath,outputPath));
writeText(configPath,jsonencode(config));
if options.injectWorkerFailure
    configPath = configPath+".missing";
end
statement = "addpath('"+replace(benchmarkFolder,"'","''")+"'); compiledPreviewBenchmarkWorker('"+replace(configPath,"'","''")+"','"+replace(outputPath,"'","''")+"')";
command = sprintf('"%s" -batch "%s"',fullfile(matlabroot,"bin","matlab"),replace(statement,'"','\"'));
[exitCode,commandOutput] = system(command);
if exitCode ~= 0 || ~isfile(outputPath)
    run = emptyRun;
    run.implementation = implementation;
    run.sourceCommit = config.sourceCommit;
    run.repeatIndex = repeatIndex;
    run.case = definition;
    run.failure = struct("identifier","WaveVortexBenchmark:CompiledPreviewWorkerFailed","message",string(commandOutput),"report",string(commandOutput));
else
    run = normalizeRun(jsondecode(fileread(outputPath)));
end
clear cleanup
end

function run = normalizeRun(value)
run = struct("schemaVersion",string(value.schemaVersion),"status",string(value.status),"implementation",string(value.implementation),"sourceCommit",string(value.sourceCommit),"repeatIndex",value.repeatIndex,"case",value.case,"constructionSeconds",value.constructionSeconds,"rawSeconds",value.rawSeconds,"medianSeconds",value.medianSeconds,"ledger",value.ledger,"metadata",value.metadata,"rss",value.rss,"rssSamples",value.rssSamples,"lifecycle",value.lifecycle,"failure",value.failure);
end

function records = correctnessRecords(cases)
records = repmat(struct("id","","maximumRelativeError",NaN),numel(cases),1);
for iCase = 1:numel(cases)
    definition = cases(iCase);
    matlabWVT = WVTransformConstantStratification([15000 15000 1300],definition.Nxyz,isHydrostatic=definition.isHydrostatic,shouldAntialias=definition.shouldAntialias);
    compiledWVT = WVTransformConstantStratification([15000 15000 1300],definition.Nxyz,isHydrostatic=definition.isHydrostatic,shouldAntialias=definition.shouldAntialias,computationalBackend="compiled");
    cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
    state = initializeWaveVortexBenchmarkState(matlabWVT,definition.seed);
    advanceWaveVortexBenchmarkState(matlabWVT,state,definition.warmupCount+1);
    compiledWVT.Ap = matlabWVT.Ap;
    compiledWVT.Am = matlabWVT.Am;
    compiledWVT.A0 = matlabWVT.A0;
    compiledWVT.t = matlabWVT.t;
    expected = cell(1,3);
    actual = cell(1,3);
    [expected{:}] = matlabWVT.nonlinearFlux();
    [actual{:}] = compiledWVT.nonlinearFlux();
    errorValues = cellfun(@relativeError,actual,expected);
    records(iCase) = struct("id",definition.id,"maximumRelativeError",max(errorValues));
    clear cleanup
end
end

function comparison = comparisonRecords(runs,cases,correctness,capabilities)
comparison = repmat(struct("id","","Nxyz",[],"isHydrostatic",false,"status","failed","matlabSeconds",NaN,"compiledSeconds",NaN,"compiledSpeedup",NaN,"matlabProcessMedians",[],"compiledProcessMedians",[],"maximumRelativeError",NaN,"matlabExactRetainedBytes",NaN,"compiledExactRetainedBytes",NaN,"exactRetainedRatio",NaN,"matlabOperationPeakIncrementRSSBytes",NaN,"compiledOperationPeakIncrementRSSBytes",NaN,"operationPeakRSSRatio",NaN,"libraryIdentityPassed",false,"nativeExecutionPassed",false,"noFallback",false,"lifecyclePassed",false,"planCount",NaN,"persistentFullHermitianBytes",NaN),numel(cases),1);
for iCase = 1:numel(cases)
    selected = runs(string(arrayfun(@(item)item.case.id,runs,"UniformOutput",false))==cases(iCase).id);
    matlabRuns = selected(string({selected.implementation})=="matlab");
    compiledRuns = selected(string({selected.implementation})=="compiled");
    matlabTimes = [matlabRuns.medianSeconds];
    compiledTimes = [compiledRuns.medianSeconds];
    matlabExact = arrayfun(@(item)item.ledger.exactRetainedApplicationBytes,matlabRuns);
    compiledExact = arrayfun(@(item)item.ledger.exactRetainedApplicationBytes,compiledRuns);
    matlabRSS = arrayfun(@(item)item.rss.operationPeakIncrementBytes,matlabRuns);
    compiledRSS = arrayfun(@(item)item.rss.operationPeakIncrementBytes,compiledRuns);
    metadata = compiledRuns(1).metadata;
    identity = all(arrayfun(@(item)string(item.metadata.module.sha256)==string(capabilities.module.sha256) && item.metadata.module.identityValidated,compiledRuns));
    lifecycle = all(arrayfun(@(item)item.lifecycle.passed,compiledRuns));
    comparison(iCase) = struct( ...
        "id",string(cases(iCase).id), ...
        "Nxyz",cases(iCase).Nxyz, ...
        "isHydrostatic",cases(iCase).isHydrostatic, ...
        "status",conditional(all(string({selected.status})=="complete"),"complete","failed"), ...
        "matlabSeconds",median(matlabTimes), ...
        "compiledSeconds",median(compiledTimes), ...
        "compiledSpeedup",median(matlabTimes)/median(compiledTimes), ...
        "matlabProcessMedians",matlabTimes, ...
        "compiledProcessMedians",compiledTimes, ...
        "maximumRelativeError",correctness(iCase).maximumRelativeError, ...
        "matlabExactRetainedBytes",median(matlabExact), ...
        "compiledExactRetainedBytes",median(compiledExact), ...
        "exactRetainedRatio",median(compiledExact./matlabExact), ...
        "matlabOperationPeakIncrementRSSBytes",median(matlabRSS), ...
        "compiledOperationPeakIncrementRSSBytes",median(compiledRSS), ...
        "operationPeakRSSRatio",median(compiledRSS./matlabRSS), ...
        "libraryIdentityPassed",identity, ...
        "nativeExecutionPassed",string(metadata.provider.id)=="native-neon-pthreads", ...
        "noFallback",string(metadata.activeBackend)=="compiled", ...
        "lifecyclePassed",lifecycle, ...
        "planCount",metadata.runtimeMetrics.planCount, ...
        "persistentFullHermitianBytes",metadata.runtimeMetrics.persistentFullHermitianBytes);
end
end

function writeArtifacts(results,options)
if ~options.shouldWriteArtifacts
    return
end
writeText(fullfile(options.outputDirectory,"compiled-preview-benchmark.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(options.outputDirectory,"summary.md"),markdownSummary(results));
end

function checkpoint(results,options)
if options.shouldWriteArtifacts
    writeArtifacts(results,options);
end
end

function value = markdownSummary(results)
lines = ["# Compiled backend preview";"";"- Status: `"+results.status+"`";"- Decision: `"+fieldOr(results.decision,"status","pending")+"`";"- Provider: native FFTW 3.3.11 NEON/pthreads";"";"| Case | MATLAB (ms) | Compiled (ms) | Speedup | Error | Exact retained ratio | Operation RSS ratio |";"|---|---:|---:|---:|---:|---:|---:|"];
for item = results.comparison'
    lines(end+1) = "| "+item.id+" | "+sprintf('%.3f',1e3*item.matlabSeconds)+" | "+sprintf('%.3f',1e3*item.compiledSeconds)+" | "+sprintf('%.3fx',item.compiledSpeedup)+" | "+sprintf('%.3e',item.maximumRelativeError)+" | "+sprintf('%.3f',item.exactRetainedRatio)+" | "+sprintf('%.3f',item.operationPeakRSSRatio)+" |"; %#ok<AGROW>
end
lines = [lines;"";"Memory is descriptive and does not gate preview availability. Exact retained bytes count reachable application-owned arrays; allocator, MATLAB FFT, and FFTW plan-owned memory remain opaque and are represented by isolated operation RSS."];
if string(results.failure.identifier) ~= ""
    lines = [lines;"";"## Failure";"";"- Stage: `"+results.failure.stage+"`";"- `"+results.failure.identifier+"`: "+results.failure.message];
end
value = join(lines,newline)+newline;
end

function value = fieldOr(record,name,fallback)
if isstruct(record) && isfield(record,name)
    value = string(record.(name));
else
    value = fallback;
end
end

function deleteTransforms(varargin)
for iTransform = 1:numel(varargin)
    if ~isempty(varargin{iTransform}) && isvalid(varargin{iTransform})
        delete(varargin{iTransform});
    end
end
end

function value = relativeError(actual,expected)
value = max(abs(actual(:)-expected(:)),[],"omitmissing")/max(max(abs(expected(:)),[],"omitmissing"),realmin);
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

function [commit,tree,isDirty] = gitIdentity(root)
commit = gitValue(root,"rev-parse HEAD");
tree = gitValue(root,"rev-parse HEAD^{tree}");
[~,output] = system("git -C "+shellQuote(root)+" status --porcelain --untracked-files=no");
isDirty = strlength(strtrim(string(output))) > 0;
end

function value = gitValue(root,arguments)
[status,output] = system("git -C "+shellQuote(root)+" "+arguments);
if status ~= 0
    error("WaveVortexBenchmark:CompiledPreviewGit","Git command failed: %s",output);
end
value = string(strtrim(output));
end

function value = environmentRecord
value = struct("host",string(getenv("HOSTNAME")),"processor",processorName,"physicalMemoryBytes",physicalMemoryBytes,"os",string(system_dependent("getos")),"matlabVersion",string(version),"release",string(version("-release")),"architecture",string(computer("arch")),"platform",string(computer),"requestedThreads",maxNumCompThreads);
end

function value = processorName
[status,output] = system("sysctl -n machdep.cpu.brand_string");
value = strtrim(string(output));
if status ~= 0 || value == ""
    value = string(system_dependent("getcpu"));
end
end

function value = physicalMemoryBytes
[status,output] = system("sysctl -n hw.memsize");
value = str2double(strtrim(output));
if status ~= 0 || ~isfinite(value)
    value = NaN;
end
end

function quoted = shellQuote(value)
quoted = "'"+replace(string(value),"'","'""'""'")+"'";
end

function writeText(pathname,value)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to write %s",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",value);
clear cleanup
end

function deleteTemporaryFiles(varargin)
for iFile = 1:numel(varargin)
    if isfile(varargin{iFile})
        delete(varargin{iFile});
    end
end
end

function restoreState(directory,originalPath,originalRng)
cd(directory);
path(originalPath);
rng(originalRng);
end

function clearModuleIfIdle(moduleName)
if moduleName == ""
    return
end
try
    metrics = feval(char(moduleName),'moduleMetrics');
    if metrics.kernelCount == 0 && ~metrics.moduleLocked
        eval("clear "+moduleName);
    end
catch
end
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function value = utcTimestamp
value = string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function value = emptyFailure
value = struct("stage","","identifier","","message","","report","");
end

function value = emptyRun
value = struct("schemaVersion","1.0.0","status","failed","implementation","","sourceCommit","","repeatIndex",0,"case",struct(),"constructionSeconds",NaN,"rawSeconds",[],"medianSeconds",NaN,"ledger",struct(),"metadata",struct(),"rss",struct(),"rssSamples",struct(),"lifecycle",struct(),"failure",struct("identifier","","message","","report",""));
end
