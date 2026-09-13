function results = runWaveVortexRetirementBenchmark(options)
% Compare the candidate builtin path with archived WaveVortex 4.2.1.
arguments
    options.baselineRef (1,1) string = "v4.2.1"
    options.candidateRef (1,1) string = "HEAD"
    options.baselineRoot (1,1) string = ""
    options.candidateRoot (1,1) string = ""
    options.caseIds (1,:) string = strings(1,0)
    options.processRunCount (1,1) double {mustBeInteger,mustBePositive} = 3
    options.correctnessTolerance (1,1) double {mustBePositive} = 1e-12
    options.performanceTolerance (1,1) double {mustBeNonnegative} = 0.03
    options.shouldRunStorageBenchmark (1,1) logical = true
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
addpath(benchmarkFolder,"-begin");

if options.runId == ""
    options.runId = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmss'Z'"));
end
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId + "-fftw-retirement-" + computer("arch") + "-" + version("-release"));
end

snapshotFolder = string(tempname);
mkdir(snapshotFolder);
snapshotCleanup = onCleanup(@()rmdir(snapshotFolder,"s"));
[baselineRoot,baselineIdentity] = prepareSource(options.baselineRoot,options.baselineRef,repositoryRoot,fullfile(snapshotFolder,"baseline"));
[candidateRoot,candidateIdentity] = prepareSource(options.candidateRoot,options.candidateRef,repositoryRoot,fullfile(snapshotFolder,"candidate"));
dependencyPath = sourceIndependentPath(originalPath,repositoryRoot);

suite = waveVortexBenchmarkSuites("core-v1");
cases = suite.cases;
if ~isempty(options.caseIds)
    cases = cases(ismember(string({cases.id}),options.caseIds));
end
if isempty(cases)
    error("WaveVortexBenchmark:NoRetirementCases","No core-v1 cases matched the retirement benchmark selection.");
end

caseResults = emptyCases();
for iCase = 1:numel(cases)
    runs = emptyRuns();
    for iRun = 1:options.processRunCount
        runs(end+1) = runWorker(cases(iCase),iRun,baselineRoot,candidateRoot,dependencyPath,benchmarkFolder); %#ok<AGROW>
    end
    caseResults(end+1) = aggregateCase(cases(iCase),runs,options.correctnessTolerance,options.performanceTolerance); %#ok<AGROW>
end

storage = struct("status","not-run");
if options.shouldRunStorageBenchmark
    storageDirectory = fullfile(options.outputDirectory,"storage");
    storage = runWaveVortexBuiltinStorageBenchmark(suiteId="core-v1",caseIds=string({cases.id}),processRunCount=options.processRunCount,outputDirectory=storageDirectory,shouldWriteArtifacts=options.shouldWriteArtifacts,runId=options.runId + "-storage");
end
passed = all([caseResults.passed]);
results = struct( ...
    "schemaVersion","1.0.0", ...
    "status",conditional(passed,"complete","failed"), ...
    "decision","RETIRE", ...
    "runId",options.runId, ...
    "environment",environmentRecord(), ...
    "configuration",struct("suiteId","core-v1","caseIds",string({cases.id}),"processRunCount",options.processRunCount,"correctnessTolerance",options.correctnessTolerance,"performanceTolerance",options.performanceTolerance,"warmupPolicy","core-v1","samplePolicy","7 medium / 3 large"), ...
    "baseline",baselineIdentity, ...
    "candidate",candidateIdentity, ...
    "cases",caseResults, ...
    "storage",storage, ...
    "gates",struct("correctnessPassed",all([caseResults.correctnessPassed]),"performancePassed",all([caseResults.performancePassed]),"builtinExecutionPassed",all([caseResults.builtinExecutionPassed]),"passed",passed));

if options.shouldWriteArtifacts
    writeArtifacts(results,options.outputDirectory);
end
clear snapshotCleanup stateCleanup
end

function [sourceRoot,identity] = prepareSource(suppliedRoot,sourceRef,repositoryRoot,destination)
if suppliedRoot ~= ""
    sourceRoot = suppliedRoot;
    identity = sourceIdentity(sourceRoot,sourceRef);
    return
end
mkdir(destination);
archivePath = destination + ".tar";
cleanup = onCleanup(@()deleteIfPresent(archivePath));
command = sprintf('git -C "%s" archive --format=tar --output="%s" "%s"',repositoryRoot,archivePath,sourceRef);
[status,output] = system(command);
if status ~= 0
    error("WaveVortexBenchmark:SourceArchiveFailed","Unable to archive %s: %s",sourceRef,output);
end
untar(archivePath,destination);
sourceRoot = destination;
identity = sourceIdentity(repositoryRoot,sourceRef);
clear cleanup
end

function run = runWorker(benchmarkCase,repeatIndex,baselineRoot,candidateRoot,dependencyPath,benchmarkFolder)
configPath = string(tempname) + ".json";
outputPath = string(tempname) + ".json";
cleanup = onCleanup(@()deleteTemporaryFiles(configPath,outputPath));
config = struct("benchmarkCase",benchmarkCase,"repeatIndex",repeatIndex,"baselineRoot",baselineRoot,"candidateRoot",candidateRoot,"dependencyPath",dependencyPath);
writeText(configPath,jsonencode(config));
matlabExecutable = fullfile(matlabroot,"bin","matlab");
statement = "addpath('" + replace(benchmarkFolder,"'","''") + "','-begin'); waveVortexRetirementWorker('" + replace(configPath,"'","''") + "','" + replace(outputPath,"'","''") + "')";
command = sprintf('"%s" -batch "%s"',matlabExecutable,replace(statement,'"','\"'));
[exitCode,commandOutput] = system(command);
if exitCode ~= 0 || ~isfile(outputPath)
    run = failedRun(repeatIndex,"WaveVortexBenchmark:RetirementWorkerFailed",string(commandOutput));
    return
end
run = jsondecode(fileread(outputPath));
run.status = string(run.status);
run.baseline.implementation = string(run.baseline.implementation);
run.candidate.implementation = string(run.candidate.implementation);
run.baseline.adapterClass = string(run.baseline.adapterClass);
run.candidate.adapterClass = string(run.candidate.adapterClass);
run.failure.identifier = string(run.failure.identifier);
run.failure.message = string(run.failure.message);
clear cleanup
end

function result = aggregateCase(definition,runs,correctnessTolerance,performanceTolerance)
complete = string({runs.status}) == "complete";
baselineMedians = NaN(1,numel(runs));
candidateMedians = NaN(1,numel(runs));
errors = Inf(1,numel(runs));
builtin = false(1,numel(runs));
for iRun = 1:numel(runs)
    if complete(iRun)
        baselineMedians(iRun) = runs(iRun).baseline.medianSeconds;
        candidateMedians(iRun) = runs(iRun).candidate.medianSeconds;
        errors(iRun) = runs(iRun).relativeError;
        builtin(iRun) = runs(iRun).baseline.implementation == "builtin" && runs(iRun).candidate.implementation == "builtin" && runs(iRun).baseline.adapterClass == "WVFastTransformDoublyPeriodicMatlab" && runs(iRun).candidate.adapterClass == "WVFastTransformDoublyPeriodicMatlab";
    end
end
baselineMedian = median(baselineMedians,"omitnan");
candidateMedian = median(candidateMedians,"omitnan");
ratio = candidateMedian/baselineMedian;
correctnessPassed = all(complete) && all(errors <= correctnessTolerance);
performancePassed = all(complete) && ratio <= 1+performanceTolerance;
builtinExecutionPassed = all(complete) && all(builtin);
result = struct("id",definition.id,"transformId",definition.transformId,"Nxyz",definition.Nxyz,"isHydrostatic",definition.isHydrostatic,"shouldAntialias",definition.shouldAntialias,"seed",definition.seed,"warmupCount",definition.warmupCount,"sampleCount",definition.sampleCount,"status",conditional(all(complete),"complete","failed"),"runs",runs,"baselineProcessMedians",baselineMedians,"candidateProcessMedians",candidateMedians,"baselineMedianSeconds",baselineMedian,"candidateMedianSeconds",candidateMedian,"candidateRelativeToBaseline",ratio,"maximumRelativeError",max(errors),"correctnessPassed",correctnessPassed,"performancePassed",performancePassed,"builtinExecutionPassed",builtinExecutionPassed,"passed",correctnessPassed && performancePassed && builtinExecutionPassed);
end

function identity = sourceIdentity(repositoryRoot,sourceRef)
[commitStatus,commitOutput] = system(sprintf('git -C "%s" rev-parse "%s^{commit}"',repositoryRoot,sourceRef));
[treeStatus,treeOutput] = system(sprintf('git -C "%s" rev-parse "%s^{tree}"',repositoryRoot,sourceRef));
identity = struct("ref",sourceRef,"commit",conditional(commitStatus == 0,string(strtrim(commitOutput)),"untracked-source"),"tree",conditional(treeStatus == 0,string(strtrim(treeOutput)),"untracked-source"));
end

function value = sourceIndependentPath(originalPath,repositoryRoot)
entries = split(string(originalPath),pathsep);
workspaceRoot = string(fileparts(repositoryRoot));
lowerEntries = lower(entries);
isWaveVortexSource = startsWith(entries,workspaceRoot) & (contains(lowerEntries,"wave-vortex-model") | contains(lowerEntries,"wavevortexmodel-"));
value = strjoin(entries(~isWaveVortexSource),pathsep);
end

function environment = environmentRecord()
environment = struct("matlabVersion",string(version),"matlabRelease",string(version("-release")),"architecture",string(computer("arch")),"os",string(system_dependent("getos")),"processor",processorName(),"memoryBytes",physicalMemoryBytes());
end

function value = processorName()
[status,output] = system("/usr/sbin/sysctl -n machdep.cpu.brand_string 2>/dev/null");
value = conditional(status == 0,string(strtrim(output)),"");
end

function value = physicalMemoryBytes()
[status,output] = system("/usr/sbin/sysctl -n hw.memsize 2>/dev/null");
value = conditional(status == 0,str2double(strtrim(output)),NaN);
end

function writeArtifacts(results,outputDirectory)
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writeText(fullfile(outputDirectory,"retirement-benchmark.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(outputDirectory,"summary.md"),summaryMarkdown(results));
end

function value = summaryMarkdown(results)
lines = ["# FFTW integration retirement gate";"";"- Decision: `" + results.decision + "`";"- Status: `" + results.status + "`";"- Baseline: `" + results.baseline.commit + "` (`v4.2.1`)";"- Candidate: `" + results.candidate.commit + "`";"- Gate: error at most `1e-12`; no candidate median over `1.03x` the same-host baseline";"";"| Case | v4.2.1 (ms) | Candidate (ms) | Candidate / baseline | Error | Builtin executed | Pass |";"|---|---:|---:|---:|---:|---|---|"];
for benchmarkCase = results.cases
    lines(end+1) = sprintf("| %s | %.3f | %.3f | %.3f | %.3g | %s | %s |",benchmarkCase.id,1e3*benchmarkCase.baselineMedianSeconds,1e3*benchmarkCase.candidateMedianSeconds,benchmarkCase.candidateRelativeToBaseline,benchmarkCase.maximumRelativeError,yesNo(benchmarkCase.builtinExecutionPassed),yesNo(benchmarkCase.passed)); %#ok<AGROW>
end
if isfield(results.storage,"cases")
    lines = [lines;"";"## Candidate builtin storage and RSS";"";"| Case | Known persistent (MiB) | Maximum known live (MiB) | Persistent RSS increment (MiB) | Peak RSS increment (MiB) |";"|---|---:|---:|---:|---:|"];
    for benchmarkCase = results.storage.cases
        lines(end+1) = sprintf("| %s | %.3f | %.3f | %.3f | %.3f |",benchmarkCase.id,benchmarkCase.ledger.knownPersistentBytes/2^20,benchmarkCase.ledger.knownMaximumLiveBytes/2^20,benchmarkCase.rss.medianPersistentIncrementBytes/2^20,benchmarkCase.rss.medianPeakIncrementBytes/2^20); %#ok<AGROW>
    end
end
value = join(lines,newline) + newline;
end

function value = yesNo(condition)
value = conditional(condition,"yes","no");
end

function run = failedRun(repeatIndex,identifier,message)
emptyImplementation = struct("implementation","","adapterClass","","rawSeconds",[],"medianSeconds",NaN);
run = struct("repeatIndex",repeatIndex,"status","failed","executionOrder",strings(0,1),"baseline",emptyImplementation,"candidate",emptyImplementation,"relativeError",Inf,"failure",struct("identifier",identifier,"message",message));
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
    deleteIfPresent(varargin{iFile});
end
end

function deleteIfPresent(pathname)
if isfile(pathname)
    delete(pathname);
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
value = struct("repeatIndex",{},"status",{},"executionOrder",{},"baseline",{},"candidate",{},"relativeError",{},"failure",{});
end

function value = emptyCases()
value = struct("id",{},"transformId",{},"Nxyz",{},"isHydrostatic",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"sampleCount",{},"status",{},"runs",{},"baselineProcessMedians",{},"candidateProcessMedians",{},"baselineMedianSeconds",{},"candidateMedianSeconds",{},"candidateRelativeToBaseline",{},"maximumRelativeError",{},"correctnessPassed",{},"performancePassed",{},"builtinExecutionPassed",{},"passed",{});
end
