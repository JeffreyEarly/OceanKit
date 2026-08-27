function waveVortexBenchmarkMemoryWorker(configPath,outputPath)
% Execute one benchmark case in a fresh MATLAB process for RSS accounting.
arguments
    configPath (1,1) string
    outputPath (1,1) string
end
config = jsondecode(fileread(configPath));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
cleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
result = emptyResult();
try
    path(config.matlabPath);
    addpath(config.repositoryRoot,config.benchmarkFolder);
    metadata = jsondecode(fileread(fullfile(config.repositoryRoot,"resources","mpackage.json")));
    for iFolder = 1:numel(metadata.folders)
        folder = fullfile(config.repositoryRoot,metadata.folders(iFolder).path);
        if isfolder(folder)
            addpath(folder);
        end
    end
    [baselineBytes,provider] = currentResidentBytes();
    wvt = createWaveVortexBenchmarkTransform(config.benchmarkCase,config.backendId);
    state = initializeWaveVortexBenchmarkState(wvt,config.benchmarkCase.seed);
    for iWarmup = 1:config.benchmarkCase.warmupCount
        advanceWaveVortexBenchmarkState(wvt,state,iWarmup);
        executeWaveVortexBenchmarkOperation(wvt,config.benchmarkCase.operation);
    end
    persistentBytes = currentResidentBytes();
    peakBytes = persistentBytes;
    nPeakSamples = max(3,config.benchmarkCase.sampleCount);
    for iSample = 1:nPeakSamples
        advanceWaveVortexBenchmarkState(wvt,state,config.benchmarkCase.warmupCount+iSample);
        outputs = executeWaveVortexBenchmarkOperation(wvt,config.benchmarkCase.operation); %#ok<NASGU>
        peakBytes = max(peakBytes,currentResidentBytes());
    end
    result = struct("status","complete","provider",provider,"baselineBytes",baselineBytes,"persistentBytes",persistentBytes,"peakBytes",peakBytes,"persistentIncrementBytes",max(0,persistentBytes-baselineBytes),"peakIncrementBytes",max(0,peakBytes-baselineBytes),"samplingIntervalSeconds",0,"failure",emptyFailure());
catch exception
    result.status = "failed";
    result.failure = struct("identifier",string(exception.identifier),"message",string(exception.message),"stack",string({exception.stack.name}));
end
writeText(outputPath,jsonencode(result,PrettyPrint=true));
clear cleanup
end

function [bytes,provider] = currentResidentBytes()
if ispc
    userMemory = memory;
    bytes = double(userMemory.MemUsedMATLAB);
    provider = "matlab-memory";
    return
end
pid = feature('getpid');
[status,output] = system(sprintf('ps -o rss= -p %d',pid));
if status ~= 0
    error("WaveVortexBenchmark:RSSUnavailable","Unable to query resident memory: %s",output);
end
kilobytes = str2double(strtrim(output));
if ~isfinite(kilobytes)
    error("WaveVortexBenchmark:RSSParseFailed","Unable to parse resident memory from: %s",output);
end
bytes = 1024*kilobytes;
if ismac
    provider = "macos-ps-rss";
else
    provider = "linux-ps-rss";
end
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

function result = emptyResult()
result = struct("status","failed","provider","","baselineBytes",NaN,"persistentBytes",NaN,"peakBytes",NaN,"persistentIncrementBytes",NaN,"peakIncrementBytes",NaN,"samplingIntervalSeconds",NaN,"failure",emptyFailure());
end

function failure = emptyFailure()
failure = struct("identifier","","message","","stack",strings(0,1));
end
