function waveVortexBuiltinStorageWorker(configPath,outputPath)
% Run one builtin storage and RSS measurement in an isolated MATLAB process.
arguments
    configPath (1,1) string
    outputPath (1,1) string
end

config = jsondecode(fileread(configPath));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
phasePath = string(tempname) + ".phase";
stopPath = string(tempname) + ".stop";
samplePath = string(tempname) + ".rss";
temporaryCleanup = onCleanup(@()deleteTemporaryFiles(phasePath,stopPath,samplePath));
sampler = emptySampler(config.samplingIntervalSeconds);
result = emptyResult(config);
wvt = [];

try
    path(config.matlabPath);
    addRepositoryPaths(config.repositoryRoot,config.benchmarkFolder);
    writePhase(phasePath,"baseline");
    sampler = startSampler(config.samplerPath,phasePath,stopPath,samplePath,config.samplingIntervalSeconds);
    pause(config.plateauSeconds);

    writePhase(phasePath,"construction");
    wvt = createWaveVortexBenchmarkTransform(config.benchmarkCase,"builtin");
    adapterClass = string(class(wvt.fastTransform));
    if adapterClass ~= "WVFastTransformDoublyPeriodicMatlab"
        error("WaveVortexBenchmark:UnexpectedStorageBackend","Expected the builtin adapter, but %s executed.",adapterClass);
    end
    state = initializeWaveVortexBenchmarkState(wvt,config.benchmarkCase.seed);

    writePhase(phasePath,"warmup");
    for iWarmup = 1:config.benchmarkCase.warmupCount
        advanceWaveVortexBenchmarkState(wvt,state,iWarmup);
        outputs = executeWaveVortexBenchmarkOperation(wvt,config.benchmarkCase.operation); %#ok<NASGU>
        clear outputs
    end

    writePhase(phasePath,"persistent");
    drawnow;
    pause(config.plateauSeconds);
    ledger = wvt.transformStorageLedger();

    writePhase(phasePath,"nonlinearFlux");
    advanceWaveVortexBenchmarkState(wvt,state,config.benchmarkCase.warmupCount+1);
    outputs = executeWaveVortexBenchmarkOperation(wvt,config.benchmarkCase.operation); %#ok<NASGU>
    drawnow;
    pause(max(0.03,3*config.samplingIntervalSeconds));
    clear outputs

    writePhase(phasePath,"cleanup");
    delete(wvt);
    clear wvt
    pause(max(0.05,2*config.samplingIntervalSeconds));
    sampler = stopSampler(sampler,stopPath,config.samplingIntervalSeconds);
    rss = samplerResult(sampler,samplePath,config.samplingIntervalSeconds);
    metadata = struct("implementation","builtin","adapterClass",adapterClass,"transformClass",string(config.benchmarkCase.transformId),"fourierStorageType",ledger.fourierStorageType);
    result = struct("status","complete","implementation","builtin","ledger",ledger,"rss",rss,"metadata",metadata,"failure",emptyFailure());
catch exception
    if ~isempty(wvt) && isvalid(wvt)
        delete(wvt);
    end
    sampler = stopSampler(sampler,stopPath,config.samplingIntervalSeconds);
    result.status = "failed";
    result.failure = exceptionFailure(exception);
    if isfile(samplePath)
        try
            result.rss = samplerResult(sampler,samplePath,config.samplingIntervalSeconds);
        catch
        end
    end
end
writeText(outputPath,jsonencode(result,PrettyPrint=true));
clear temporaryCleanup stateCleanup
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

function sampler = startSampler(samplerPath,phasePath,stopPath,samplePath,interval)
sampler = emptySampler(interval);
if ~(ismac || isunix) || ~isfile(samplerPath)
    sampler.reason = "External ps-based RSS sampling is unavailable on this platform.";
    return
end
pid = feature("getpid");
command = sprintf('"/bin/sh" "%s" "%d" "%s" "%s" "%s" "%.6f" >/dev/null 2>&1 & echo $!',samplerPath,pid,phasePath,stopPath,samplePath,interval);
[status,output] = system(command);
samplerPid = str2double(strtrim(output));
if status ~= 0 || ~isfinite(samplerPid)
    sampler.reason = "Unable to launch the external RSS sampler.";
    return
end
sampler.status = "running";
sampler.processId = samplerPid;
sampler.provider = conditional(ismac,"macos-ps-rss-external","linux-ps-rss-external");
end

function sampler = stopSampler(sampler,stopPath,interval)
if sampler.status ~= "running"
    return
end
writeText(stopPath,"stop");
pause(max(0.05,2*interval));
system(sprintf('kill %d >/dev/null 2>&1',sampler.processId));
sampler.status = "complete";
end

function rss = samplerResult(sampler,samplePath,interval)
rss = emptyRSS(interval);
rss.status = sampler.status;
rss.provider = sampler.provider;
rss.reason = sampler.reason;
if sampler.status ~= "complete" || ~isfile(samplePath)
    return
end
rss.samples = readSamples(samplePath,interval);
baseline = phaseValues(rss.samples,"baseline");
persistentValues = phaseValues(rss.samples,"persistent");
nonlinearValues = phaseValues(rss.samples,"nonlinearFlux");
construction = phaseValues(rss.samples,"construction");
warmup = phaseValues(rss.samples,"warmup");
if isempty(baseline) || isempty(persistentValues) || isempty(nonlinearValues)
    rss.status = "unsupported";
    rss.reason = "External sampling did not capture every required phase.";
    return
end
rss.baselineMedianBytes = median(baseline);
rss.constructionPeakBytes = maximumOrNaN(construction);
rss.warmupPeakBytes = maximumOrNaN(warmup);
rss.persistentMedianBytes = median(persistentValues);
rss.nonlinearFluxPeakBytes = max(nonlinearValues);
rss.persistentIncrementBytes = rss.persistentMedianBytes-rss.baselineMedianBytes;
rss.peakIncrementBytes = rss.nonlinearFluxPeakBytes-rss.baselineMedianBytes;
end

function samples = readSamples(samplePath,interval)
lines = splitlines(strtrim(string(fileread(samplePath))));
samples = repmat(struct("sampleIndex",0,"elapsedSeconds",0,"phase","","rssBytes",0),numel(lines),1);
for iLine = 1:numel(lines)
    fields = split(lines(iLine),sprintf('\t'));
    if numel(fields) ~= 3
        error("WaveVortexBenchmark:InvalidRSSSample","Invalid external RSS sample: %s",lines(iLine));
    end
    index = str2double(fields(1));
    rssKiB = str2double(fields(3));
    if ~isfinite(index) || ~isfinite(rssKiB)
        error("WaveVortexBenchmark:InvalidRSSSample","Invalid numeric external RSS sample: %s",lines(iLine));
    end
    samples(iLine) = struct("sampleIndex",index,"elapsedSeconds",index*interval,"phase",fields(2),"rssBytes",1024*rssKiB);
end
end

function values = phaseValues(samples,phase)
values = [samples(string({samples.phase}) == phase).rssBytes];
end

function value = maximumOrNaN(values)
if isempty(values)
    value = NaN;
else
    value = max(values);
end
end

function writePhase(pathname,phase)
temporaryPath = pathname + ".tmp";
writeText(temporaryPath,phase);
movefile(temporaryPath,pathname,"f");
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

function value = emptySampler(interval)
value = struct("status","unsupported","provider","","reason","","processId",NaN,"samplingIntervalSeconds",interval);
end

function value = emptyRSS(interval)
value = struct("status","unsupported","provider","","reason","","samplingIntervalSeconds",interval,"samples",repmat(struct("sampleIndex",0,"elapsedSeconds",0,"phase","","rssBytes",0),0,1),"baselineMedianBytes",NaN,"constructionPeakBytes",NaN,"warmupPeakBytes",NaN,"persistentMedianBytes",NaN,"nonlinearFluxPeakBytes",NaN,"persistentIncrementBytes",NaN,"peakIncrementBytes",NaN);
end

function value = emptyResult(config)
value = struct("status","failed","implementation","builtin","ledger",struct(),"rss",emptyRSS(config.samplingIntervalSeconds),"metadata",struct(),"failure",emptyFailure());
end

function value = emptyFailure()
value = struct("identifier","","message","","stack",strings(0,1));
end

function value = exceptionFailure(exception)
stack = repmat(struct("name","","file","","line",0),numel(exception.stack),1);
for iFrame = 1:numel(exception.stack)
    stack(iFrame) = struct("name",string(exception.stack(iFrame).name),"file",string(exception.stack(iFrame).file),"line",exception.stack(iFrame).line);
end
value = struct("identifier",string(exception.identifier),"message",string(exception.message),"stack",stack);
end
