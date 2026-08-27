function [matlabDataset,compiledDataset] = publishedWaveVortexBenchmarksFromCompiledPreviewArtifact(rawArtifactPath,options)
% Normalize one canonical compiled-preview artifact into two public datasets.
arguments
    rawArtifactPath (1,1) string {mustBeFile}
    options.platformId (1,1) string = "m5-max"
    options.platformName (1,1) string = "Apple M5 Max"
    options.provenancePath (1,1) string
    options.implementationVersion (1,1) string = "unreleased-preview"
end
raw = jsondecode(fileread(rawArtifactPath));
if string(raw.schemaVersion) ~= "compiled-preview-benchmark-v1" || string(raw.status) ~= "complete"
    error("WaveVortexBenchmark:InvalidCompiledPreviewArtifact","The raw artifact must be a complete compiled-preview-benchmark-v1 result.");
end
[collectedAt,timestamp] = collectionTime(raw.runId);
matlabDataset = datasetFor("matlab",raw,collectedAt,timestamp,options);
compiledDataset = datasetFor("compiled",raw,collectedAt,timestamp,options);
end

function dataset = datasetFor(implementationId,raw,collectedAt,timestamp,options)
isCompiled = implementationId == "compiled";
publicImplementationId = conditional(isCompiled,"cpp","matlab");
backendId = conditional(isCompiled,"native-fftw","builtin");
displayName = conditional(isCompiled,"WaveVortexModel compiled preview","WaveVortexModel MATLAB");
datasetId = "core-v1--"+publicImplementationId+"-"+backendId+"--"+options.platformId+"--"+timestamp;
benchmark = struct("suiteId","core-v1","suiteVersion",1,"operation","nonlinearFlux","correctnessTolerance",double(raw.configuration.correctnessTolerance));
implementation = struct("id",publicImplementationId,"displayName",displayName,"version",options.implementationVersion,"repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",string(raw.source.commit),"backend",backendId,"sourceDirty",logical(raw.source.isDirty));
platform = struct("id",options.platformId,"displayName",options.platformName,"processor",string(raw.environment.processor),"physicalMemoryBytes",double(raw.environment.physicalMemoryBytes),"os",string(raw.environment.os),"architecture",string(raw.environment.architecture),"threadCount",double(raw.environment.requestedThreads));
if isCompiled
    details = struct("provider",string(raw.provider.provider.id),"fftwVersion",string(raw.provider.provider.version),"threadBackend",string(raw.provider.provider.threadBackend),"matlabHostRelease",string(raw.environment.release));
    toolchain = struct("kind","cpp","name","Apple Clang + FFTW","version",string(raw.provider.compiler.mexVersion),"details",details);
else
    details = struct("matlabRelease",string(raw.environment.release));
    toolchain = struct("kind","matlab","name","MATLAB","version",string(raw.environment.matlabVersion),"details",details);
end
provenance = struct("rawArtifact",options.provenancePath,"rawSchemaVersion",string(raw.schemaVersion));
cases = cell(1,numel(raw.comparison));
for iCase = 1:numel(raw.comparison)
    comparison = raw.comparison(iCase);
    definition = raw.cases(iCase);
    implementationMask = reshape(string({raw.runs.implementation}),[],1) == implementationId;
    caseMask = reshape(string(arrayfun(@(item)item.case.id,raw.runs,"UniformOutput",false)),[],1) == comparison.id;
    runs = raw.runs(implementationMask & caseMask);
    samples = [runs.rawSeconds];
    if isCompiled
        medianSeconds = comparison.compiledSeconds;
        exactRetainedBytes = comparison.compiledExactRetainedBytes;
    else
        medianSeconds = comparison.matlabSeconds;
        exactRetainedBytes = comparison.matlabExactRetainedBytes;
    end
    retainedValues = arrayfun(@(item)item.rss.steadyRetainedBytes,runs);
    peakValues = arrayfun(@(item)item.rss.operationPeakBytes,runs);
    baselineBytes = median(retainedValues);
    peakProcessBytes = median(peakValues);
    peakIncrementBytes = median(arrayfun(@(item)item.rss.operationPeakIncrementBytes,runs));
    configuration = struct("Lxyz",[15000 15000 1300],"Nxyz",double(definition.Nxyz(:)'),"isHydrostatic",logical(definition.isHydrostatic),"shouldAntialias",logical(definition.shouldAntialias),"seed",double(definition.seed),"warmupCount",double(definition.warmupCount),"sampleCount",double(definition.sampleCount));
    timing = struct("medianSeconds",double(medianSeconds),"samplesSeconds",double(samples(:)'));
    correctness = struct("passed",comparison.maximumRelativeError<=raw.configuration.correctnessTolerance,"relativeError",double(comparison.maximumRelativeError));
    memory = struct("status","complete","provider","macos-ps-rss-external","baselineProcessBytes",double(baselineBytes),"peakProcessBytes",double(peakProcessBytes),"peakIncrementBytes",double(peakIncrementBytes),"exactRetainedBytes",double(exactRetainedBytes),"exactRetainedDefinition",string(raw.configuration.exactScope));
    cases{iCase} = struct("id",string(comparison.id),"transformId",conditional(definition.isHydrostatic,"constant-hydrostatic","constant-nonhydrostatic"),"scoreFamily",conditional(definition.isHydrostatic,"constant-hydrostatic","constant-nonhydrostatic"),"configuration",configuration,"status","complete","timing",timing,"correctness",correctness,"memory",memory);
end
dataset = struct("schemaVersion","published-benchmark-v1","datasetId",datasetId,"collectedAt",collectedAt,"benchmark",benchmark,"implementation",implementation,"platform",platform,"toolchain",toolchain,"provenance",provenance,"cases",{cases});
end

function [collectedAt,timestamp] = collectionTime(runId)
token = regexp(string(runId),'^(\d{8}T\d{6})\d{0,3}Z$','tokens','once');
if isempty(token)
    error("WaveVortexBenchmark:InvalidCollectionTime","Compiled-preview runId must use a UTC YYYYMMDDTHHmmss[SSS]Z timestamp.");
end
timestamp = string(token{1})+"Z";
instant = datetime(timestamp,"InputFormat","yyyyMMdd'T'HHmmss'Z'","TimeZone","UTC");
collectedAt = string(instant,"yyyy-MM-dd'T'HH:mm:ss'Z'");
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
