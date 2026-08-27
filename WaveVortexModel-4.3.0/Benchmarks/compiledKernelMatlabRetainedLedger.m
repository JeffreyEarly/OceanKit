function ledger = compiledKernelMatlabRetainedLedger(wvt,outputs,backendMetadata)
% Count application-owned arrays retained by one public backend call.
arguments
    wvt WVTransformConstantStratification {mustBeNonempty}
    outputs (1,3) cell
    backendMetadata (1,1) struct = struct()
end

transform = wvt.transformStorageLedger();
cacheNames = string(wvt.variableCache.keys);
cacheEntries = repmat(struct("name","","className","","shape",[],"bytes",0),numel(cacheNames),1);
cacheBytes = 0;
for iName = 1:numel(cacheNames)
    value = wvt.variableCache{cacheNames(iName)};
    info = whos("value");
    cacheEntries(iName) = struct("name",cacheNames(iName),"className",string(class(value)),"shape",double(size(value)),"bytes",double(info.bytes));
    cacheBytes = cacheBytes + info.bytes;
end

outputEntries = repmat(struct("name","","className","","shape",[],"bytes",0),3,1);
outputBytes = 0;
names = ["Fp" "Fm" "F0"];
for iOutput = 1:3
    value = outputs{iOutput};
    info = whos("value");
    outputEntries(iOutput) = struct("name",names(iOutput),"className",string(class(value)),"shape",double(size(value)),"bytes",double(info.bytes));
    outputBytes = outputBytes + info.bytes;
end

compiledPersistentBytes = 0;
compiledDescriptorBytes = 0;
compiledScratchBytes = 0;
compiledPlanBytes = 0;
if isfield(backendMetadata,"runtimeMetrics") && isfield(backendMetadata.runtimeMetrics,"persistentBytes")
    metrics = backendMetadata.runtimeMetrics;
    compiledPersistentBytes = double(metrics.persistentBytes);
    compiledDescriptorBytes = double(metrics.descriptorBytes);
    compiledScratchBytes = double(metrics.scratchCapacityBytes);
    compiledPlanBytes = double(metrics.planBytes);
end

ledger = struct( ...
    "schemaVersion","compiled-preview-retained-memory-v1", ...
    "scope","active MATLAB transform storage, reachable variable-cache arrays, compiled backend-owned storage when active, and three returned flux arrays; canonical Ap/Am/A0 state excluded", ...
    "backend",string(wvt.computationalBackend), ...
    "transformRetainedBytes",double(transform.knownPersistentBytes), ...
    "cacheRetainedBytes",double(cacheBytes), ...
    "compiledPersistentBytes",compiledPersistentBytes, ...
    "compiledDescriptorBytes",compiledDescriptorBytes, ...
    "compiledScratchBytes",compiledScratchBytes, ...
    "compiledPlanWrapperLowerBoundBytes",compiledPlanBytes, ...
    "outputRetainedBytes",double(outputBytes), ...
    "exactRetainedApplicationBytes",double(transform.knownPersistentBytes+cacheBytes+compiledPersistentBytes+outputBytes), ...
    "cacheEntries",cacheEntries, ...
    "outputEntries",outputEntries, ...
    "opaqueMemory","MATLAB allocator, JIT temporaries, copy-on-write sharing, MATLAB FFT internals, and FFTW plan-owned allocations beyond the wrapper lower bound are represented only by isolated RSS");
end
