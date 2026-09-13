function results = runWaveVortexBuiltinStorageBenchmark(options)
% Measure exact builtin transform storage and repeated process RSS.
arguments
    options.suiteId (1,1) string = "core-v1"
    options.caseIds (1,:) string = strings(1,0)
    options.processRunCount (1,1) double {mustBeInteger,mustBePositive} = 3
    options.samplingIntervalSeconds (1,1) double {mustBePositive} = 0.02
    options.plateauSeconds (1,1) double {mustBePositive} = 0.15
    options.outputDirectory (1,1) string = ""
    options.runId (1,1) string = ""
    options.shouldWriteArtifacts (1,1) logical = true
end

benchmarkFolder = string(fileparts(mfilename("fullpath")));
repositoryRoot = string(fileparts(benchmarkFolder));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
addRepositoryPaths(repositoryRoot,benchmarkFolder);

suite = waveVortexBenchmarkSuites(options.suiteId);
cases = suite.cases(startsWith(string({suite.cases.transformId}),"constant-"));
if ~isempty(options.caseIds)
    cases = cases(ismember(string({cases.id}),options.caseIds));
end
if isempty(cases)
    error("WaveVortexBenchmark:NoStorageCases","No constant-stratification cases matched the requested storage benchmark selection.");
end
if options.runId == ""
    options.runId = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmss'Z'"));
end
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId + "-builtin-storage-" + computer("arch") + "-" + version("-release"));
end

caseResults = emptyCases();
for iCase = 1:numel(cases)
    runs = emptyRuns();
    for iRepeat = 1:options.processRunCount
        runs(end+1) = runWorker(cases(iCase),iRepeat,options,benchmarkFolder,repositoryRoot); %#ok<AGROW>
    end
    caseResults(end+1) = aggregateCase(cases(iCase),runs); %#ok<AGROW>
end
status = "complete";
if any(string({caseResults.status}) ~= "complete")
    status = "partial";
end
[sourceCommit,sourceTree] = sourceIdentity(repositoryRoot);
results = struct( ...
    "schemaVersion","1.0.0", ...
    "status",status, ...
    "runId",options.runId, ...
    "environment",struct("matlabVersion",string(version),"matlabRelease",string(version("-release")),"architecture",string(computer("arch")),"os",string(system_dependent("getos"))), ...
    "source",struct("repository","JeffreyEarly/wave-vortex-model","commit",sourceCommit,"tree",sourceTree), ...
    "configuration",struct("suiteId",options.suiteId,"caseIds",options.caseIds,"processRunCount",options.processRunCount,"samplingIntervalSeconds",options.samplingIntervalSeconds,"plateauSeconds",options.plateauSeconds), ...
    "cases",caseResults);

if options.shouldWriteArtifacts
    writeArtifacts(results,options.outputDirectory);
end
clear stateCleanup
end

function run = runWorker(benchmarkCase,repeatIndex,options,benchmarkFolder,repositoryRoot)
configPath = string(tempname) + ".json";
outputPath = string(tempname) + ".json";
cleanup = onCleanup(@()deleteTemporaryFiles(configPath,outputPath));
config = struct("benchmarkCase",benchmarkCase,"repeatIndex",repeatIndex,"repositoryRoot",repositoryRoot,"benchmarkFolder",benchmarkFolder,"matlabPath",path,"samplerPath",fullfile(benchmarkFolder,"sampleProcessRSS.sh"),"samplingIntervalSeconds",options.samplingIntervalSeconds,"plateauSeconds",options.plateauSeconds);
writeText(configPath,jsonencode(config));
matlabExecutable = fullfile(matlabroot,"bin","matlab");
statement = "addpath('" + replace(benchmarkFolder,"'","''") + "'); waveVortexBuiltinStorageWorker('" + replace(configPath,"'","''") + "','" + replace(outputPath,"'","''") + "')";
command = sprintf('"%s" -batch "%s"',matlabExecutable,replace(statement,'"','\"'));
[exitCode,commandOutput] = system(command);
if exitCode ~= 0 || ~isfile(outputPath)
    run = struct("repeatIndex",repeatIndex,"status","failed","implementation","builtin","ledger",struct(),"rss",struct(),"metadata",struct(),"failure",struct("identifier","WaveVortexBenchmark:BuiltinStorageWorkerFailed","message",string(commandOutput)));
    return
end
decoded = jsondecode(fileread(outputPath));
run = struct("repeatIndex",repeatIndex,"status",string(decoded.status),"implementation",string(decoded.implementation),"ledger",decoded.ledger,"rss",decoded.rss,"metadata",decoded.metadata,"failure",decoded.failure);
clear cleanup
end

function result = aggregateCase(definition,runs)
complete = string({runs.status}) == "complete";
rssComplete = complete & arrayfun(@(run)isfield(run.rss,"status") && string(run.rss.status) == "complete",runs);
persistentBytes = NaN(1,numel(runs));
peak = NaN(1,numel(runs));
for iRun = 1:numel(runs)
    if rssComplete(iRun)
        persistentBytes(iRun) = runs(iRun).rss.persistentIncrementBytes;
        peak(iRun) = runs(iRun).rss.peakIncrementBytes;
    end
end
ledger = struct();
if any(complete)
    ledger = runs(find(complete,1)).ledger;
end
status = conditional(all(complete),"complete","partial");
rssStatus = conditional(all(rssComplete),"complete","unsupported");
rss = struct("status",rssStatus,"persistentIncrementBytes",persistentBytes,"peakIncrementBytes",peak,"medianPersistentIncrementBytes",median(persistentBytes,"omitnan"),"minimumPersistentIncrementBytes",min(persistentBytes,[],"omitnan"),"maximumPersistentIncrementBytes",max(persistentBytes,[],"omitnan"),"medianPeakIncrementBytes",median(peak,"omitnan"),"minimumPeakIncrementBytes",min(peak,[],"omitnan"),"maximumPeakIncrementBytes",max(peak,[],"omitnan"));
result = struct("id",definition.id,"transformId",definition.transformId,"Nxyz",definition.Nxyz,"isHydrostatic",definition.isHydrostatic,"shouldAntialias",definition.shouldAntialias,"seed",definition.seed,"warmupCount",definition.warmupCount,"status",status,"implementation","builtin","runs",runs,"ledger",ledger,"rss",rss);
end

function writeArtifacts(results,outputDirectory)
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writeText(fullfile(outputDirectory,"builtin-storage.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(outputDirectory,"summary.md"),summaryMarkdown(results));
end

function text = summaryMarkdown(results)
lines = [ ...
    "# Builtin transform storage benchmark"; ...
    ""; ...
    "- Status: `" + results.status + "`"; ...
    "- Source: `" + results.source.commit + "`"; ...
    "- MATLAB: `" + results.environment.matlabRelease + "`"; ...
    "- Architecture: `" + results.environment.architecture + "`"; ...
    ""; ...
    "| Case | Known persistent (MiB) | Maximum known live (MiB) | Opaque entries | Persistent RSS (MiB) | Peak RSS (MiB) |"; ...
    "|---|---:|---:|---:|---:|---:|"];
for benchmarkCase = results.cases
    if benchmarkCase.status == "complete"
        lines(end+1) = sprintf("| %s | %.3f | %.3f | %d | %.3f | %.3f |",benchmarkCase.id,benchmarkCase.ledger.knownPersistentBytes/2^20,benchmarkCase.ledger.knownMaximumLiveBytes/2^20,benchmarkCase.ledger.opaqueEntryCount,benchmarkCase.rss.medianPersistentIncrementBytes/2^20,benchmarkCase.rss.medianPeakIncrementBytes/2^20); %#ok<AGROW>
    else
        lines(end+1) = "| " + benchmarkCase.id + " | unavailable | unavailable | unavailable | unavailable | unavailable |"; %#ok<AGROW>
    end
end
text = join(lines,newline) + newline;
end

function [commit,tree] = sourceIdentity(repositoryRoot)
[commitStatus,commitOutput] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
[treeStatus,treeOutput] = system(sprintf('git -C "%s" rev-parse HEAD^{tree}',repositoryRoot));
commit = conditional(commitStatus == 0,string(strtrim(commitOutput)),"");
tree = conditional(treeStatus == 0,string(strtrim(treeOutput)),"");
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

function restoreState(originalDirectory,originalPath,originalRng)
cd(originalDirectory);
path(originalPath);
rng(originalRng);
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function value = emptyRuns()
value = struct("repeatIndex",{},"status",{},"implementation",{},"ledger",{},"rss",{},"metadata",{},"failure",{});
end

function value = emptyCases()
value = struct("id",{},"transformId",{},"Nxyz",{},"isHydrostatic",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"status",{},"implementation",{},"runs",{},"ledger",{},"rss",{});
end
