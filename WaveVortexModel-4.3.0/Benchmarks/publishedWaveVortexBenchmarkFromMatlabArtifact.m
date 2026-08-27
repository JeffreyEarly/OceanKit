function dataset = publishedWaveVortexBenchmarkFromMatlabArtifact(rawArtifactPath,options)
% Normalize one MATLAB suite and backend for public benchmark use.
arguments
    rawArtifactPath (1,1) string
    options.suiteId (1,1) string
    options.platformId (1,1) string
    options.platformName (1,1) string
    options.provenancePath (1,1) string
    options.backendId (1,1) string = "builtin"
    options.processorName (1,1) string = ""
    options.implementationVersion (1,1) string = ""
    options.implementationRepository (1,1) string = "https://github.com/JeffreyEarly/wave-vortex-model"
end

if ~isfile(rawArtifactPath)
    error("WaveVortexBenchmark:MissingRawArtifact","Raw benchmark artifact does not exist: %s.",rawArtifactPath);
end
raw = jsondecode(fileread(rawArtifactPath));
if ~isfield(raw,"schemaVersion")
    error("WaveVortexBenchmark:InvalidRawArtifact","Raw benchmark artifact is missing schemaVersion.");
end
rawSchemaVersion = string(raw.schemaVersion);
if ~ismember(rawSchemaVersion,["1.0.0" "1.1.0"])
    error("WaveVortexBenchmark:UnsupportedRawSchema","Unsupported raw benchmark schema: %s.",rawSchemaVersion);
end
suite = selectedSuite(raw.suites,options.suiteId);
if ~isfield(suite,"operation") || ~ismember(string(suite.operation),["nonlinearAdvection" "nonlinearFlux"])
    error("WaveVortexBenchmark:UnsupportedPublishedOperation","Published benchmark v1 supports only the state-advanced nonlinearFlux contract.");
end
environment = raw.environment;
configuration = raw.configuration;
if ~isfield(environment,"sourceDirty") || logical(environment.sourceDirty)
    error("WaveVortexBenchmark:DirtyPublishedSource","Published benchmark artifacts require a clean source commit.");
end

processorName = options.processorName;
if processorName == "" && isfield(environment,"processorName")
    processorName = string(environment.processorName);
end
if processorName == "" && isfield(environment,"processor")
    processorName = string(environment.processor);
end
if processorName == ""
    error("WaveVortexBenchmark:MissingProcessorName","A human-readable processor name is required.");
end

implementationVersion = options.implementationVersion;
if implementationVersion == "" && isfield(environment,"packageVersion")
    implementationVersion = string(environment.packageVersion);
end
if implementationVersion == ""
    error("WaveVortexBenchmark:MissingImplementationVersion","An implementation version is required for this raw artifact.");
end

[collectedAt,timestamp] = collectionTime(raw.runId);
datasetId = options.suiteId + "--matlab-" + options.backendId + "--" + options.platformId + "--" + timestamp;
benchmark = struct("suiteId",string(suite.id),"suiteVersion",double(suite.version),"operation","nonlinearFlux","correctnessTolerance",double(configuration.correctnessTolerance));
displayName = "WaveVortexModel MATLAB";
if isfield(environment,"packageName")
    displayName = string(environment.packageName) + " MATLAB";
end
implementation = struct("id","matlab","displayName",displayName,"version",implementationVersion,"repository",options.implementationRepository,"commit",string(environment.sourceCommit),"backend",options.backendId,"sourceDirty",logical(environment.sourceDirty));
platform = struct("id",options.platformId,"displayName",options.platformName,"processor",processorName,"physicalMemoryBytes",double(environment.physicalMemoryBytes),"os",string(environment.os),"architecture",string(environment.architecture),"threadCount",double(environment.requestedThreads));
details = struct("matlabRelease",string(environment.matlabRelease));
toolchain = struct("kind","matlab","name","MATLAB","version",string(environment.matlabVersion),"details",details);
provenance = struct("rawArtifact",options.provenancePath,"rawSchemaVersion",rawSchemaVersion);

rawCases = suite.cases;
cases = cell(1,numel(rawCases));
for iCase = 1:numel(rawCases)
    if iscell(rawCases)
        rawCase = rawCases{iCase};
    else
        rawCase = rawCases(iCase);
    end
    cases{iCase} = normalizeCase(rawCase,options.backendId);
end
dataset = struct("schemaVersion","published-benchmark-v1","datasetId",datasetId,"collectedAt",collectedAt,"benchmark",benchmark,"implementation",implementation,"platform",platform,"toolchain",toolchain,"provenance",provenance,"cases",{cases});
end

function benchmarkCase = normalizeCase(rawCase,backendId)
if ~isfield(rawCase,"operation") || ~ismember(string(rawCase.operation),["nonlinearAdvection" "nonlinearFlux"])
    error("WaveVortexBenchmark:UnsupportedPublishedOperation","Raw cases must measure the state-advanced nonlinearFlux contract.");
end
configuration = struct("Lxyz",double(rawCase.Lxyz(:)'),"Nxyz",double(rawCase.Nxyz(:)'),"isHydrostatic",logical(rawCase.isHydrostatic),"shouldAntialias",logical(rawCase.shouldAntialias),"seed",double(rawCase.seed),"warmupCount",double(rawCase.warmupCount),"sampleCount",double(rawCase.sampleCount));
common = {"id",string(rawCase.id),"transformId",string(rawCase.transformId),"scoreFamily",string(rawCase.scoreFamily),"configuration",configuration};
if string(rawCase.status) ~= "complete"
    error("WaveVortexBenchmark:UnpublishableRawCase","Raw case %s did not complete and cannot be published.",string(rawCase.id));
end
backends = rawCase.backends;
backendIndex = 0;
for iBackend = 1:numel(backends)
    if iscell(backends)
        backend = backends{iBackend};
    else
        backend = backends(iBackend);
    end
    if string(backend.id) == backendId
        backendIndex = iBackend;
        break
    end
end
if backendIndex == 0
    benchmarkCase = struct(common{:},"status","unavailable","unavailableReason","Backend " + backendId + " was not recorded in the raw artifact.");
    return
end
if iscell(backends)
    backend = backends{backendIndex};
else
    backend = backends(backendIndex);
end
if string(backend.status) ~= "complete" || ~logical(backend.correctnessPassed)
    error("WaveVortexBenchmark:UnpublishableRawCase","Raw backend %s for case %s did not pass.",backendId,string(rawCase.id));
end
samples = double(backend.rawSeconds(:)');
timing = struct("medianSeconds",double(backend.medianSeconds),"samplesSeconds",samples);
correctness = struct("passed",logical(backend.correctnessPassed),"relativeError",double(backend.relativeError));
memory = normalizeMemory(backend.memory);
benchmarkCase = struct(common{:},"status","complete","timing",timing,"correctness",correctness,"memory",memory);
end

function memory = normalizeMemory(rawMemory)
if isfield(rawMemory,"status") && string(rawMemory.status) == "complete"
    memory = struct("status","complete","provider",string(rawMemory.provider),"baselineProcessBytes",double(rawMemory.baselineBytes),"peakProcessBytes",double(rawMemory.peakBytes),"peakIncrementBytes",double(rawMemory.peakIncrementBytes));
    return
end
reason = "Memory measurement was unavailable in the raw artifact.";
if isfield(rawMemory,"status") && string(rawMemory.status) == "not-requested"
    reason = "Memory measurement was not requested.";
elseif isfield(rawMemory,"failure") && isfield(rawMemory.failure,"message") && strlength(string(rawMemory.failure.message)) > 0
    reason = string(rawMemory.failure.message);
end
memory = struct("status","unavailable","reason",reason);
end

function suite = selectedSuite(suites,suiteId)
for iSuite = 1:numel(suites)
    if iscell(suites)
        candidate = suites{iSuite};
    else
        candidate = suites(iSuite);
    end
    if string(candidate.id) == suiteId
        suite = candidate;
        return
    end
end
error("WaveVortexBenchmark:MissingRawSuite","Raw benchmark artifact does not contain suite %s.",suiteId);
end

function [collectedAt,timestamp] = collectionTime(runId)
runId = string(runId);
if isempty(regexp(runId,"^\d{8}T\d{6}Z$","once"))
    error("WaveVortexBenchmark:InvalidCollectionTime","Raw benchmark runId must use UTC YYYYMMDDTHHmmssZ format.");
end
instant = datetime(runId,"InputFormat","yyyyMMdd'T'HHmmss'Z'","TimeZone","UTC");
collectedAt = string(instant,"yyyy-MM-dd'T'HH:mm:ss'Z'");
timestamp = string(instant,"yyyyMMdd'T'HHmmss'Z'");
end
