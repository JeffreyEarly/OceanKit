function compiledPreviewBenchmarkWorker(configPath,outputPath)
% Run one public backend/case/repeat in a fresh MATLAB process.
config = jsondecode(fileread(configPath));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
phasePath = string(tempname)+".phase";
stopPath = string(tempname)+".stop";
samplePath = string(tempname)+".rss";
temporaryCleanup = onCleanup(@()deleteTemporaryFiles(phasePath,stopPath,samplePath));
sampler = emptySampler(config.samplingIntervalSeconds);
wvt = [];
result = emptyResult(config);
try
    path(config.matlabPath);
    addRepositoryPaths(config.repositoryRoot,config.benchmarkFolder);
    writePhase(phasePath,"startup");
    sampler = startSampler(config.samplerPath,phasePath,stopPath,samplePath,config.samplingIntervalSeconds);
    pause(config.plateauSeconds);

    definition = config.caseDefinition;
    backend = string(config.implementation);
    rng(definition.seed,"twister");
    writePhase(phasePath,"backend-construction");
    constructionTimer = tic;
    wvt = WVTransformConstantStratification([15000 15000 1300],definition.Nxyz,isHydrostatic=definition.isHydrostatic,shouldAntialias=definition.shouldAntialias,computationalBackend=backend);
    constructionSeconds = toc(constructionTimer);
    state = initializeWaveVortexBenchmarkState(wvt,definition.seed);
    advanceWaveVortexBenchmarkState(wvt,state,0);
    metadata = wvt.computationalBackendMetadata;
    validateActiveBackend(metadata,backend,string(config.expectedModuleHash));
    writePhase(phasePath,"backend-created");
    drawnow;
    pause(config.plateauSeconds);

    for iWarmup = 1:definition.warmupCount
        advanceWaveVortexBenchmarkState(wvt,state,iWarmup);
        writePhase(phasePath,"warmup");
        execute(wvt);
    end
    writePhase(phasePath,"steady-retained");
    drawnow;
    pause(config.plateauSeconds);

    rawSeconds = NaN(1,definition.sampleCount);
    outputs = cell(1,3);
    for iSample = 1:definition.sampleCount
        advanceWaveVortexBenchmarkState(wvt,state,definition.warmupCount+iSample);
        writePhase(phasePath,"operation");
        timer = tic;
        outputs = execute(wvt);
        rawSeconds(iSample) = toc(timer);
        writePhase(phasePath,"outputs-held");
        drawnow;
        pause(config.outputHoldSeconds);
        if iSample < definition.sampleCount
            clear outputs
            outputs = cell(1,3);
            writePhase(phasePath,"outputs-cleared");
        end
    end
    metadata = wvt.computationalBackendMetadata;
    ledger = compiledKernelMatlabRetainedLedger(wvt,outputs,metadata);
    clear outputs
    writePhase(phasePath,"outputs-cleared");
    drawnow;
    pause(config.plateauSeconds);

    writePhase(phasePath,"backend-destruction");
    delete(wvt);
    wvt = [];
    lifecycle = lifecycleRecord(backend);
    writePhase(phasePath,"backend-destroyed");
    drawnow;
    pause(config.plateauSeconds);
    sampler = stopSampler(sampler,stopPath,config.samplingIntervalSeconds);
    rssSamples = samplerResult(sampler,samplePath,config.samplingIntervalSeconds);
    rss = phaseRSS(rssSamples,config.samplingIntervalSeconds);
    completed = lifecycle.passed && rss.status == "complete";
    result = struct( ...
        "schemaVersion","1.0.0", ...
        "status",conditional(completed,"complete","failed"), ...
        "implementation",backend, ...
        "sourceCommit",string(config.sourceCommit), ...
        "repeatIndex",config.repeatIndex, ...
        "case",definition, ...
        "constructionSeconds",constructionSeconds, ...
        "rawSeconds",rawSeconds, ...
        "medianSeconds",median(rawSeconds), ...
        "ledger",ledger, ...
        "metadata",metadata, ...
        "rss",rss, ...
        "rssSamples",rssSamples, ...
        "lifecycle",lifecycle, ...
        "failure",emptyFailure);
catch exception
    if ~isempty(wvt) && isvalid(wvt)
        delete(wvt);
    end
    stopSampler(sampler,stopPath,config.samplingIntervalSeconds);
    result.failure = struct("identifier",string(exception.identifier),"message",string(exception.message),"report",string(getReport(exception,"extended","hyperlinks","off")));
end
writeText(outputPath,jsonencode(result,PrettyPrint=true));
clear temporaryCleanup stateCleanup
end

function outputs = execute(wvt)
outputs = cell(1,3);
[outputs{:}] = wvt.nonlinearFlux();
end

function validateActiveBackend(metadata,expected,expectedModuleHash)
if string(metadata.activeBackend) ~= expected
    error("WaveVortexBenchmark:CompiledPreviewFallback","Requested %s but %s executed.",expected,string(metadata.activeBackend));
end
if expected == "compiled"
    if string(metadata.provider.id) ~= "native-neon-pthreads" || ~metadata.module.identityValidated || string(metadata.module.sha256) ~= expectedModuleHash || metadata.contract.planCount ~= 17 || metadata.libraries.openmp.detected
        error("WaveVortexBenchmark:CompiledPreviewIdentity","The compiled worker did not execute the validated native provider contract.");
    end
end
end

function value = lifecycleRecord(backend)
value = struct("passed",true,"kernelCount",0,"activePlans",0,"outstandingPlanningBytes",0,"moduleLocked",false);
if backend ~= "compiled"
    return
end
metrics = wv_compiled_backend_mex('moduleMetrics');
value = struct( ...
    "passed",metrics.kernelCount == 0 && metrics.activePlans == 0 && metrics.outstandingPlanningBytes == 0 && ~metrics.moduleLocked, ...
    "kernelCount",metrics.kernelCount, ...
    "activePlans",metrics.activePlans, ...
    "outstandingPlanningBytes",metrics.outstandingPlanningBytes, ...
    "moduleLocked",metrics.moduleLocked);
end

function value = phaseRSS(rss,interval)
value = struct("status",rss.status,"provider",rss.provider,"samplingIntervalSeconds",interval,"steadyRetainedBytes",NaN,"operationPeakBytes",NaN,"operationPeakIncrementBytes",NaN,"backendDestroyedBytes",NaN,"constructionPeakBytes",NaN);
if rss.status ~= "complete" || isempty(rss.samples)
    return
end
phases = string({rss.samples.phase});
bytes = [rss.samples.rssBytes];
retained = bytes(phases=="steady-retained");
operation = bytes(phases=="operation" | phases=="outputs-held");
destroyed = bytes(phases=="backend-destroyed");
construction = bytes(phases=="backend-construction" | phases=="backend-created");
if isempty(retained) || isempty(operation) || isempty(destroyed)
    value.status = "unsupported";
    return
end
value.steadyRetainedBytes = median(retained);
value.operationPeakBytes = max(operation);
value.operationPeakIncrementBytes = max(0,value.operationPeakBytes-value.steadyRetainedBytes);
value.backendDestroyedBytes = median(destroyed);
if ~isempty(construction)
    value.constructionPeakBytes = max(construction);
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

function sampler = startSampler(samplerPath,phasePath,stopPath,samplePath,interval)
sampler = emptySampler(interval);
if ~(ismac||isunix) || ~isfile(samplerPath)
    sampler.reason = "External RSS sampler unavailable.";
    return
end
command = sprintf('"/bin/sh" "%s" "%d" "%s" "%s" "%s" "%.6f" >/dev/null 2>&1 & echo $!',samplerPath,feature('getpid'),phasePath,stopPath,samplePath,interval);
[status,output] = system(command);
samplerPid = str2double(strtrim(output));
if status ~= 0 || ~isfinite(samplerPid)
    sampler.reason = "Unable to launch RSS sampler.";
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
rss = struct("status",sampler.status,"provider",sampler.provider,"reason",sampler.reason,"samplingIntervalSeconds",interval,"samples",[]);
if sampler.status ~= "complete" || ~isfile(samplePath)
    return
end
lines = splitlines(strtrim(string(fileread(samplePath))));
samples = repmat(struct("sampleIndex",0,"elapsedSeconds",0,"phase","","rssBytes",0),numel(lines),1);
for iLine = 1:numel(lines)
    fields = split(lines(iLine),sprintf('\t'));
    index = str2double(fields(1));
    samples(iLine) = struct("sampleIndex",index,"elapsedSeconds",index*interval,"phase",fields(2),"rssBytes",1024*str2double(fields(3)));
end
rss.samples = samples;
end

function writePhase(pathname,phase)
temporary = pathname+".tmp";
writeText(temporary,phase);
movefile(temporary,pathname,"f");
end

function writeText(pathname,contents)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to open %s",pathname);
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

function restoreState(directory,originalPath,originalRng)
cd(directory);
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

function value = emptyFailure
value = struct("identifier","","message","","report","");
end

function value = emptyResult(config)
value = struct("schemaVersion","1.0.0","status","failed","implementation",string(config.implementation),"sourceCommit",string(config.sourceCommit),"repeatIndex",config.repeatIndex,"case",config.caseDefinition,"constructionSeconds",NaN,"rawSeconds",[],"medianSeconds",NaN,"ledger",struct(),"metadata",struct(),"rss",struct(),"rssSamples",struct(),"lifecycle",struct(),"failure",emptyFailure);
end
