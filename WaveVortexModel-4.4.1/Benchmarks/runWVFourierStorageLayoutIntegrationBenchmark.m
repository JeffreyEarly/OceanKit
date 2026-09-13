function results = runWVFourierStorageLayoutIntegrationBenchmark(options)
% Compare the production Fourier-storage layout with issue #69.
arguments
    options.caseIds (1,:) string = strings(1,0)
    options.outputDirectory (1,1) string = ""
    options.shouldWriteArtifacts (1,1) logical = true
    options.runId (1,1) string = ""
    options.correctnessTolerance (1,1) double {mustBePositive} = 1e-12
    options.gateTolerance (1,1) double {mustBeNonnegative} = 0.03
end

benchmarkFolder = string(fileparts(mfilename("fullpath")));
repositoryRoot = string(fileparts(benchmarkFolder));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
addpath(repositoryRoot,benchmarkFolder);

suite = waveVortexBenchmarkSuites("transform-layout-v1");
if ~isempty(options.caseIds)
    selected = ismember(string({suite.cases.id}),options.caseIds);
    unknown = setdiff(options.caseIds,string({suite.cases.id}));
    if ~isempty(unknown)
        error("WaveVortexBenchmark:UnknownCase","Unknown transform-layout case: %s.",strjoin(unknown,", "));
    end
    suite.cases = suite.cases(selected);
end
referencePath = fullfile(benchmarkFolder,"results","reference","transform-layout-v1-m5-max-r2026a-builtin","benchmark.json");
reference = loadReference(referencePath);
if options.runId == ""
    options.runId = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmss'Z'"));
end
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId + "-storage-layout-integration-" + computer("arch") + "-" + version("-release"));
end

operationIds = ["extract" "insert-primary" "insert-conjugate" "insert-complete" "forward-complete" "inverse-complete"];
results = struct( ...
    "schemaVersion","1.0.0", ...
    "status","complete", ...
    "runId",options.runId, ...
    "environment",environmentRecord(repositoryRoot), ...
    "configuration",struct("caseIds",options.caseIds,"correctnessTolerance",options.correctnessTolerance,"gateTolerance",options.gateTolerance,"gateSizes",[256 256 65;512 512 129]), ...
    "reference",struct("issue",69,"path",relativePath(referencePath,repositoryRoot),"sha256",sha256File(referencePath),"strategy","wv-sorted-linear"), ...
    "operationIds",operationIds, ...
    "cases",emptyCaseResults(), ...
    "readiness",struct("passed",false,"completeGateSet",false,"gateCaseCount",0,"failedCriteria",strings(0,1)));

for iCase = 1:numel(suite.cases)
    try
        caseResult = runCase(suite.cases(iCase),reference,operationIds,options.correctnessTolerance,options.gateTolerance);
    catch exception
        caseResult = failedCase(suite.cases(iCase),exception);
        results.status = "partial";
    end
    results.cases(end+1) = caseResult;
end
results.readiness = readinessRecord(results.cases);
if ~results.readiness.passed
    results.status = "failed";
end
if options.shouldWriteArtifacts
    writeArtifacts(results,options.outputDirectory);
end
clear stateCleanup
end

function caseResult = runCase(benchmarkCase,reference,operationIds,tolerance,gateTolerance)
Nx = benchmarkCase.Nxyz(1);
Ny = benchmarkCase.Nxyz(2);
Nz = benchmarkCase.Nxyz(3);
rng(benchmarkCase.seed,"twister");
geometry = WVGeometryDoublyPeriodic(benchmarkCase.Lxyz(1:2),[Nx Ny],Nz=Nz,shouldAntialias=benchmarkCase.shouldAntialias,shouldExcludeNyquist=true,shouldExcludeConjugates=true,conjugateDimension=2);
transform = geometry.fastTransform;
layout = transform.fourierStorageLayout;
realInput = randn(Nx,Ny,Nz);
fullSpectrum = fft(fft(realInput,Nx,1),Ny,2)/(Nx*Ny);
fullSpectrumRows = reshape(fullSpectrum,layout.nFourierStorageRows,[]);
referenceWV = fullSpectrumRows(layout.fourierRowsForDirectWVIndices,:).';
referenceRows = layout.allocateFourierStorage(Nz);
referenceRows = layout.transformFromWVGridToFourierStorage(referenceRows,referenceWV);
referenceInverse = ifft(ifft(layout.reshapeFourierRowsToStorage(referenceRows),Nx,1),Ny,2,"symmetric")*(Nx*Ny);
buffer = WVTransformLayoutBenchmarkBuffer([Nx Ny Nz],"rows-2d");

errors = validateProduction(transform,layout,buffer,fullSpectrum,realInput,referenceWV,referenceRows,referenceInverse);
errorValues = cell2mat(struct2cell(errors));
if max(errorValues) > tolerance
    error("WaveVortexBenchmark:LayoutIntegrationCorrectnessFailed","Production layout exceeded the relative-error tolerance (maximum %.3g).",max(errorValues));
end

warmupSchedules = strings(benchmarkCase.warmupCount,numel(operationIds));
for iWarmup = 1:benchmarkCase.warmupCount
    order = rotatedOrder(numel(operationIds),iWarmup);
    warmupSchedules(iWarmup,:) = operationIds(order);
    for iOperation = order
        executeOperation(operationIds(iOperation),transform,layout,buffer,fullSpectrum,referenceWV,realInput);
    end
end

rawSeconds = NaN(numel(operationIds),benchmarkCase.sampleCount);
sampleSchedules = strings(benchmarkCase.sampleCount,numel(operationIds));
for iSample = 1:benchmarkCase.sampleCount
    order = rotatedOrder(numel(operationIds),iSample);
    sampleSchedules(iSample,:) = operationIds(order);
    for iOperation = order
        timer = tic;
        executeOperation(operationIds(iOperation),transform,layout,buffer,fullSpectrum,referenceWV,realInput);
        rawSeconds(iOperation,iSample) = toc(timer);
    end
end

referenceCase = referenceCaseFor(reference,benchmarkCase.id);
operations = repmat(struct("id","","rawSeconds",[],"medianSeconds",NaN,"referenceMedianSeconds",NaN,"relativeToReference",NaN,"relativeError",NaN,"correctnessPassed",false,"performanceGateApplies",false,"performancePassed",false),1,numel(operationIds));
isGate = isGateSize(benchmarkCase.Nxyz);
for iOperation = 1:numel(operationIds)
    samples = rawSeconds(iOperation,:);
    currentMedian = median(samples);
    baselineMedian = referenceMedian(referenceCase,operationIds(iOperation));
    ratio = currentMedian/baselineMedian;
    operations(iOperation) = struct( ...
        "id",operationIds(iOperation), ...
        "rawSeconds",samples, ...
        "medianSeconds",currentMedian, ...
        "referenceMedianSeconds",baselineMedian, ...
        "relativeToReference",ratio, ...
        "relativeError",errors.(errorField(operationIds(iOperation))), ...
        "correctnessPassed",errors.(errorField(operationIds(iOperation))) <= tolerance, ...
        "performanceGateApplies",isGate, ...
        "performancePassed",~isGate || ratio <= 1+gateTolerance);
end

caseResult = struct( ...
    "id",benchmarkCase.id, ...
    "Nxyz",benchmarkCase.Nxyz, ...
    "shouldAntialias",benchmarkCase.shouldAntialias, ...
    "seed",benchmarkCase.seed, ...
    "warmupCount",benchmarkCase.warmupCount, ...
    "sampleCount",benchmarkCase.sampleCount, ...
    "isGate",isGate, ...
    "status","complete", ...
    "failure",emptyFailure(), ...
    "mappingMethod",layout.mappingMethod, ...
    "fourierStorageSize",layout.fourierStorageSize, ...
    "mappingMemoryBytes",layout.mappingMemoryBytes, ...
    "persistentBufferBytes",complexArrayBytes(transform.complexBuffer), ...
    "warmupSchedules",warmupSchedules, ...
    "sampleSchedules",sampleSchedules, ...
    "operations",operations);
end

function errors = validateProduction(transform,layout,buffer,fullSpectrum,realInput,referenceWV,referenceRows,referenceInverse)
extracted = executeOperation("extract",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
buffer.reset();
executeOperation("insert-primary",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
primaryError = relativeError(buffer.value(layout.fourierRowsForDirectWVIndices,:),referenceRows(layout.fourierRowsForDirectWVIndices,:));
buffer.reset();
buffer.value(layout.fourierRowsForDirectWVIndices,:) = referenceRows(layout.fourierRowsForDirectWVIndices,:);
executeOperation("insert-conjugate",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
completionRows = [layout.hermitianCompletionRows;layout.selfConjugateFourierRows];
conjugateError = relativeError(buffer.value(completionRows,:),referenceRows(completionRows,:));
buffer.reset();
executeOperation("insert-complete",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
forward = executeOperation("forward-complete",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
inverse = executeOperation("inverse-complete",transform,layout,buffer,fullSpectrum,referenceWV,realInput);
errors = struct( ...
    "extract",relativeError(extracted,referenceWV), ...
    "insertPrimary",primaryError, ...
    "insertConjugate",conjugateError, ...
    "insertComplete",relativeError(buffer.value,referenceRows), ...
    "forwardComplete",relativeError(forward,referenceWV), ...
    "inverseComplete",relativeError(inverse,referenceInverse));
end

function output = executeOperation(operationId,transform,layout,buffer,fullSpectrum,wvInput,realInput)
switch operationId
    case "extract"
        output = layout.transformFromFourierStorageToWVGrid(reshape(fullSpectrum,layout.nFourierStorageRows,[]));
    case "insert-primary"
        buffer.value(layout.fourierRowsForDirectWVIndices,:) = wvInput(:,layout.directWVIndices).';
        output = [];
    case "insert-conjugate"
        if ~isempty(layout.hermitianCompletionRows)
            buffer.value(layout.hermitianCompletionRows,:) = conj(wvInput(:,layout.hermitianSourceWVIndices).');
        end
        if ~isempty(layout.selfConjugateFourierRows)
            buffer.value(layout.selfConjugateFourierRows,:) = real(buffer.value(layout.selfConjugateFourierRows,:));
        end
        output = [];
    case "insert-complete"
        buffer.value(layout.fourierRowsForDirectWVIndices,:) = wvInput(:,layout.directWVIndices).';
        if ~isempty(layout.hermitianCompletionRows)
            buffer.value(layout.hermitianCompletionRows,:) = conj(wvInput(:,layout.hermitianSourceWVIndices).');
        end
        if ~isempty(layout.selfConjugateFourierRows)
            buffer.value(layout.selfConjugateFourierRows,:) = real(buffer.value(layout.selfConjugateFourierRows,:));
        end
        output = [];
    case "forward-complete"
        output = transform.transformFromSpatialDomainWithFourier(realInput);
    case "inverse-complete"
        output = transform.transformToSpatialDomainWithFourier(wvInput);
    otherwise
        error("WaveVortexBenchmark:UnknownLayoutOperation","Unknown transform-layout operation %s.",operationId);
end
end

function reference = loadReference(pathname)
if ~isfile(pathname)
    error("WaveVortexBenchmark:MissingLayoutReference","The issue-69 reference artifact is missing: %s.",pathname);
end
reference = jsondecode(fileread(pathname));
end

function benchmarkCase = referenceCaseFor(reference,caseId)
cases = reference.suites.cases;
iCase = find(string({cases.id}) == caseId,1);
if isempty(iCase)
    error("WaveVortexBenchmark:MissingLayoutReferenceCase","The issue-69 artifact has no case %s.",caseId);
end
benchmarkCase = cases(iCase);
end

function value = referenceMedian(benchmarkCase,operationId)
iStrategy = find(string({benchmarkCase.strategies.id}) == "wv-sorted-linear",1);
operations = benchmarkCase.strategies(iStrategy).operations;
iOperation = find(string({operations.id}) == operationId,1);
value = operations(iOperation).medianSeconds;
end

function readiness = readinessRecord(cases)
gateCases = cases([cases.isGate]);
failed = strings(0,1);
for benchmarkCase = gateCases
    if benchmarkCase.status ~= "complete"
        failed(end+1,1) = benchmarkCase.id + ": " + benchmarkCase.failure.identifier; %#ok<AGROW>
        continue
    end
    for operation = benchmarkCase.operations
        if ~operation.correctnessPassed
            failed(end+1,1) = benchmarkCase.id + "/" + operation.id + ": correctness"; %#ok<AGROW>
        end
        if ~operation.performancePassed
            failed(end+1,1) = benchmarkCase.id + "/" + operation.id + ": " + sprintf("%.3fx reference",operation.relativeToReference); %#ok<AGROW>
        end
    end
end
readiness = struct("passed",isempty(failed),"completeGateSet",numel(gateCases) == 4,"gateCaseCount",numel(gateCases),"failedCriteria",failed);
end

function tf = isGateSize(Nxyz)
tf = isequal(Nxyz,[256 256 65]) || isequal(Nxyz,[512 512 129]);
end

function field = errorField(operationId)
switch operationId
    case "insert-primary"
        field = "insertPrimary";
    case "insert-conjugate"
        field = "insertConjugate";
    case "insert-complete"
        field = "insertComplete";
    case "forward-complete"
        field = "forwardComplete";
    case "inverse-complete"
        field = "inverseComplete";
    otherwise
        field = operationId;
end
end

function order = rotatedOrder(count,roundNumber)
start = mod(roundNumber-1,count)+1;
order = [start:count 1:start-1];
end

function value = relativeError(actual,expected)
if isempty(actual) && isempty(expected)
    value = 0;
    return
end
value = max(abs(actual(:)-expected(:)))/max(max(abs(expected(:))),realmin("double"));
end

function bytes = complexArrayBytes(value)
bytes = 16*numel(value);
end

function writeArtifacts(results,outputDirectory)
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
writeText(fullfile(outputDirectory,"benchmark.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(outputDirectory,"summary.md"),summaryText(results));
end

function summary = summaryText(results)
lines = ["# Fourier storage layout integration";"";"- Status: `" + results.status + "`";"- Gate passed: `" + string(results.readiness.passed) + "`";"- Reference: WaveVortex issue #69 (`wv-sorted-linear`)";"- Gate: no operation above `1.03x` the frozen median on the 256 and 512 workloads";"";"## Timing and frozen-baseline comparison";"";"| Case | Antialias | Gate | Operation | Production (ms) | Issue #69 (ms) | Relative | Error | Pass |";"|---|---:|---|---|---:|---:|---:|---:|---|"];
for benchmarkCase = results.cases
    if benchmarkCase.status ~= "complete"
        lines(end+1) = "| " + benchmarkCase.id + " | " + string(benchmarkCase.shouldAntialias) + " | " + yesNo(benchmarkCase.isGate) + " | failed | NaN | NaN | NaN | NaN | no |"; %#ok<AGROW>
        continue
    end
    for operation = benchmarkCase.operations
        lines(end+1) = sprintf("| %s | %d | %s | %s | %.3f | %.3f | %.3f | %.3g | %s |",benchmarkCase.id,benchmarkCase.shouldAntialias,yesNo(benchmarkCase.isGate),operation.id,1e3*operation.medianSeconds,1e3*operation.referenceMedianSeconds,operation.relativeToReference,operation.relativeError,yesNo(operation.correctnessPassed && operation.performancePassed)); %#ok<AGROW>
    end
end
lines = [lines;"";"## Storage contract";"";"| Case | Strategy | Mapping (MiB) | Persistent full buffer (MiB) |";"|---|---|---:|---:|"];
for benchmarkCase = results.cases
    if benchmarkCase.status == "complete"
        lines(end+1) = sprintf("| %s | %s | %.3f | %.3f |",benchmarkCase.id,benchmarkCase.mappingMethod,benchmarkCase.mappingMemoryBytes/2^20,benchmarkCase.persistentBufferBytes/2^20); %#ok<AGROW>
    end
end
if ~isempty(results.readiness.failedCriteria)
    lines = [lines;"";"## Failed criteria";""];
    for failure = string(results.readiness.failedCriteria(:)).'
        lines(end+1) = "- " + failure; %#ok<AGROW>
    end
end
summary = strjoin(lines,newline) + newline;
end

function value = yesNo(tf)
if tf
    value = "yes";
else
    value = "no";
end
end

function environment = environmentRecord(repositoryRoot)
[~,commit] = system("git -C " + quoted(repositoryRoot) + " rev-parse HEAD");
[~,tree] = system("git -C " + quoted(repositoryRoot) + " rev-parse HEAD^{tree}");
environment = struct("matlabVersion",string(version),"matlabRelease",string(version("-release")),"architecture",string(computer("arch")),"os",string(system_dependent("getos")),"sourceCommit",strtrim(string(commit)),"sourceTree",strtrim(string(tree)));
end

function value = quoted(value)
value = "'" + string(value) + "'";
end

function path = relativePath(path,root)
path = erase(string(path),string(root) + filesep);
end

function hash = sha256File(pathname)
fileId = fopen(pathname,"r");
if fileId < 0
    error("WaveVortexBenchmark:UnreadableArtifact","Unable to read %s.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
bytes = fread(fileId,Inf,"*uint8");
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(bytes);
hashBytes = typecast(digest.digest(),"uint8");
hash = lower(string(reshape(dec2hex(hashBytes,2).',1,[])));
clear cleanup
end

function writeText(pathname,text)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to write %s.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",text);
clear cleanup
end

function restoreState(originalDirectory,originalPath,originalRng)
cd(originalDirectory);
path(originalPath);
rng(originalRng);
end

function result = failedCase(benchmarkCase,exception)
result = struct("id",benchmarkCase.id,"Nxyz",benchmarkCase.Nxyz,"shouldAntialias",benchmarkCase.shouldAntialias,"seed",benchmarkCase.seed,"warmupCount",benchmarkCase.warmupCount,"sampleCount",benchmarkCase.sampleCount,"isGate",isGateSize(benchmarkCase.Nxyz),"status","failed","failure",struct("identifier",string(exception.identifier),"message",string(exception.message),"stack",string({exception.stack.name})),"mappingMethod","","fourierStorageSize",[],"mappingMemoryBytes",NaN,"persistentBufferBytes",NaN,"warmupSchedules",strings(0,0),"sampleSchedules",strings(0,0),"operations",emptyOperationResults());
end

function results = emptyCaseResults()
results = struct("id",{},"Nxyz",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"sampleCount",{},"isGate",{},"status",{},"failure",{},"mappingMethod",{},"fourierStorageSize",{},"mappingMemoryBytes",{},"persistentBufferBytes",{},"warmupSchedules",{},"sampleSchedules",{},"operations",{});
end

function results = emptyOperationResults()
results = struct("id",{},"rawSeconds",{},"medianSeconds",{},"referenceMedianSeconds",{},"relativeToReference",{},"relativeError",{},"correctnessPassed",{},"performanceGateApplies",{},"performancePassed",{});
end

function failure = emptyFailure()
failure = struct("identifier","","message","","stack",strings(0,1));
end
