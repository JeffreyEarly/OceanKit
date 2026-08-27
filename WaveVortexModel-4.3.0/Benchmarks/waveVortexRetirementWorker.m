function waveVortexRetirementWorker(configPath,outputPath)
% Compare two source snapshots in one fresh MATLAB process.
arguments
    configPath (1,1) string
    outputPath (1,1) string
end

config = jsondecode(fileread(configPath));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
comparisonFolder = string(tempname);
mkdir(comparisonFolder);
comparisonCleanup = onCleanup(@()rmdir(comparisonFolder,"s"));
result = failedResult(config.repeatIndex,"WaveVortexBenchmark:RetirementWorkerIncomplete","The retirement worker did not complete.");

try
    if mod(config.repeatIndex,2) == 1
        order = ["baseline" "candidate"];
    else
        order = ["candidate" "baseline"];
    end
    records = struct;
    referenceFiles = strings(0,1);
    for implementation = order
        sourceRoot = string(config.(implementation + "Root"));
        configureSource(sourceRoot,string(config.dependencyPath));
        [record,outputs] = runImplementation(config.benchmarkCase);
        records.(implementation) = record;
        if isempty(referenceFiles)
            referenceFiles = saveOutputs(outputs,comparisonFolder,implementation);
        else
            relativeError = compareOutputs(outputs,referenceFiles);
        end
        clear outputs
    end
    result = struct("repeatIndex",config.repeatIndex,"status","complete","executionOrder",order,"baseline",records.baseline,"candidate",records.candidate,"relativeError",relativeError,"failure",emptyFailure());
catch exception
    result = failedResult(config.repeatIndex,string(exception.identifier),string(exception.message));
end
writeText(outputPath,jsonencode(result,PrettyPrint=true));
clear comparisonCleanup stateCleanup
end

function [record,outputs] = runImplementation(benchmarkCase)
wvt = createWaveVortexBenchmarkTransform(benchmarkCase,"builtin");
cleanup = onCleanup(@()deleteTransform(wvt));
state = initializeWaveVortexBenchmarkState(wvt,benchmarkCase.seed);
for iWarmup = 1:benchmarkCase.warmupCount
    advanceWaveVortexBenchmarkState(wvt,state,iWarmup);
    outputs = executeWaveVortexBenchmarkOperation(wvt,benchmarkCase.operation); %#ok<NASGU>
end
rawSeconds = NaN(1,benchmarkCase.sampleCount);
for iSample = 1:benchmarkCase.sampleCount
    advanceWaveVortexBenchmarkState(wvt,state,benchmarkCase.warmupCount+iSample);
    timer = tic;
    outputs = executeWaveVortexBenchmarkOperation(wvt,benchmarkCase.operation);
    rawSeconds(iSample) = toc(timer);
end
adapterClass = string(class(wvt.fastTransform));
if adapterClass ~= "WVFastTransformDoublyPeriodicMatlab"
    error("WaveVortexBenchmark:UnexpectedRetirementBackend","Expected the builtin adapter, but %s executed.",adapterClass);
end
record = struct("implementation","builtin","adapterClass",adapterClass,"rawSeconds",rawSeconds,"medianSeconds",median(rawSeconds));
clear cleanup
end

function files = saveOutputs(outputs,folder,prefix)
files = strings(numel(outputs),1);
for iOutput = 1:numel(outputs)
    value = outputs{iOutput}; %#ok<NASGU>
    files(iOutput) = fullfile(folder,prefix + "-" + iOutput + ".mat");
    save(files(iOutput),"value","-v7.3");
end
end

function relativeError = compareOutputs(outputs,referenceFiles)
relativeError = 0;
for iOutput = 1:numel(outputs)
    reference = load(referenceFiles(iOutput),"value");
    numerator = max(abs(outputs{iOutput}(:)-reference.value(:)));
    denominator = max(max(abs(reference.value(:))),realmin("double"));
    relativeError = max(relativeError,numerator/denominator);
    clear reference
end
end

function configureSource(sourceRoot,dependencyPath)
path(dependencyPath);
metadata = jsondecode(fileread(fullfile(sourceRoot,"resources","mpackage.json")));
addpath(sourceRoot,"-begin");
for iFolder = numel(metadata.folders):-1:1
    addpath(fullfile(sourceRoot,metadata.folders(iFolder).path),"-begin");
end
addpath(fullfile(sourceRoot,"Benchmarks"),"-begin");
clear classes
end

function deleteTransform(wvt)
if ~isempty(wvt) && isvalid(wvt)
    delete(wvt);
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

function restoreState(originalDirectory,originalPath,originalRng)
cd(originalDirectory);
path(originalPath);
rng(originalRng);
end

function value = failedResult(repeatIndex,identifier,message)
emptyImplementation = struct("implementation","","adapterClass","","rawSeconds",[],"medianSeconds",NaN);
value = struct("repeatIndex",repeatIndex,"status","failed","executionOrder",strings(0,1),"baseline",emptyImplementation,"candidate",emptyImplementation,"relativeError",Inf,"failure",struct("identifier",identifier,"message",message));
end

function value = emptyFailure()
value = struct("identifier","","message","");
end
