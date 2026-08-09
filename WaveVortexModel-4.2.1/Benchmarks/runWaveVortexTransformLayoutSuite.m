function suiteResult = runWaveVortexTransformLayoutSuite(suite,correctnessTolerance,repositoryRoot)
% Benchmark full-spectrum mapping strategies without changing production code.
arguments
    suite (1,1) struct
    correctnessTolerance (1,1) double {mustBePositive} = 1e-12
    repositoryRoot (1,1) string = string(fileparts(fileparts(mfilename("fullpath"))))
end

strategyIds = ["wv-sorted-linear" "dft-sorted-linear" "two-dimensional-rows" "per-plane-linear"];
operationIds = ["extract" "insert-primary" "insert-conjugate" "insert-complete" "forward-complete" "inverse-complete"];
suiteResult = baseSuiteResult(suite,repositoryRoot,strategyIds,operationIds);
for iCase = 1:numel(suite.cases)
    try
        caseResult = runCase(suite.cases(iCase),strategyIds,operationIds,correctnessTolerance);
    catch exception
        caseResult = failedCaseResult(suite.cases(iCase),exception);
        suiteResult.status = "partial";
    end
    suiteResult.cases(end+1) = caseResult;
end
if ~suite.selectionIsComplete
    suiteResult.status = "partial";
end
end

function suiteResult = baseSuiteResult(suite,repositoryRoot,strategyIds,operationIds)
metadata = struct( ...
    "schema","transform-layout-v1", ...
    "currentStrategy","wv-sorted-linear", ...
    "tieTolerance",0.03, ...
    "strategyIds",strategyIds, ...
    "operationIds",operationIds, ...
    "assignmentOrder",["primary" "conjugate"], ...
    "copyObservation","unavailable", ...
    "copyObservationReason","No supported MATLAB API exposes copy-on-write or pointer state for these expressions.", ...
    "productionSources",productionSourceRecords(repositoryRoot));
suiteResult = struct("id",suite.id,"version",suite.version,"kind",suite.kind,"description",suite.description,"operation",suite.operation,"isScored",suite.isScored,"selectionIsComplete",suite.selectionIsComplete,"status","complete","cases",emptyCaseResults(),"familyScores",emptyScores(),"suiteScores",emptyScores(),"referenceArtifact","","metadata",metadata);
end

function caseResult = runCase(benchmarkCase,strategyIds,operationIds,correctnessTolerance)
Nx = benchmarkCase.Nxyz(1);
Ny = benchmarkCase.Nxyz(2);
Nz = benchmarkCase.Nxyz(3);
rng(benchmarkCase.seed,"twister");
geometry = WVGeometryDoublyPeriodic(benchmarkCase.Lxyz(1:2),[Nx Ny],Nz=Nz,shouldAntialias=benchmarkCase.shouldAntialias,shouldExcludeNyquist=true,shouldExcludeConjugates=true,conjugateDimension=2,fastTransform="builtin");
productionTransform = geometry.fastTransform;
realInput = randn(Nx,Ny,Nz);
realInputReference = realInput;
wvReference = productionTransform.transformFromSpatialDomainWithFourier(realInput);
wvInputReference = wvReference;
productionInverse = productionTransform.transformToSpatialDomainWithFourier(wvReference);
fullSpectrum = fft(fft(realInput,Nx,1),Ny,2)/(Nx*Ny);
fullBufferReference = complex(zeros(Nx,Ny,Nz));
fullBufferReference(geometry.dftPrimaryIndex) = wvReference;
fullBufferReference(geometry.dftConjugateIndex) = conj(wvReference(geometry.wvConjugateIndex));
runtimes = createStrategyRuntimes(geometry,Nz,strategyIds);
validations = repmat(emptyValidation(),1,numel(strategyIds));
for iStrategy = 1:numel(runtimes)
    validations(iStrategy) = validateStrategy(runtimes(iStrategy),fullSpectrum,wvReference,fullBufferReference,realInput,productionInverse,correctnessTolerance);
end

pairs = allOperationPairs(numel(strategyIds),numel(operationIds));
warmupSchedules = strings(benchmarkCase.warmupCount,size(pairs,1));
for iWarmup = 1:benchmarkCase.warmupCount
    order = rotatedOrder(size(pairs,1),iWarmup);
    warmupSchedules(iWarmup,:) = pairIds(pairs(order,:),strategyIds,operationIds);
    for iPair = order
        executeOperation(runtimes(pairs(iPair,1)),operationIds(pairs(iPair,2)),fullSpectrum,wvReference,realInput,Nx,Ny);
    end
end

rawSeconds = NaN(numel(strategyIds),numel(operationIds),benchmarkCase.sampleCount);
sampleSchedules = strings(benchmarkCase.sampleCount,size(pairs,1));
for iSample = 1:benchmarkCase.sampleCount
    order = rotatedOrder(size(pairs,1),iSample);
    sampleSchedules(iSample,:) = pairIds(pairs(order,:),strategyIds,operationIds);
    for iPair = order
        iStrategy = pairs(iPair,1);
        iOperation = pairs(iPair,2);
        sampleTimer = tic;
        executeOperation(runtimes(iStrategy),operationIds(iOperation),fullSpectrum,wvReference,realInput,Nx,Ny);
        rawSeconds(iStrategy,iOperation,iSample) = toc(sampleTimer);
    end
end

strategyResults = emptyStrategyResults();
for iStrategy = 1:numel(runtimes)
    operationResults = emptyOperationResults();
    for iOperation = 1:numel(operationIds)
        samples = reshape(rawSeconds(iStrategy,iOperation,:),1,[]);
        operationResults(end+1) = struct( ...
            "id",operationIds(iOperation), ...
            "rawSeconds",samples, ...
            "medianSeconds",median(samples), ...
            "relativeError",validationError(validations(iStrategy),operationIds(iOperation)), ...
            "correctnessPassed",validationPassed(validations(iStrategy),operationIds(iOperation)), ...
            "resultBytes",operationResultBytes(operationIds(iOperation),Nx,Ny,Nz,geometry.Nkl), ...
            "resultStorage",operationResultLedger(operationIds(iOperation),Nx,Ny,Nz,geometry.Nkl)); %#ok<AGROW>
    end
    strategyResults(end+1) = struct( ...
        "id",runtimes(iStrategy).id, ...
        "expressionId",runtimes(iStrategy).expressionId, ...
        "expressionDescription",runtimes(iStrategy).expressionDescription, ...
        "assignmentSteps",["primary" "conjugate"], ...
        "mappingArrays",runtimes(iStrategy).mappingLedger, ...
        "mappingBytes",sum([runtimes(iStrategy).mappingLedger.bytes]), ...
        "buffer",runtimes(iStrategy).bufferLedger, ...
        "persistentBufferReused",true, ...
        "timedBufferClearing",false, ...
        "sourceArraysUnchanged",isequaln(realInput,realInputReference) && isequaln(wvReference,wvInputReference), ...
        "copyObservation","unavailable", ...
        "copyObservationReason","No supported MATLAB API exposes copy-on-write or pointer state for these expressions.", ...
        "operations",operationResults); %#ok<AGROW>
end

selections = emptySelections();
for iOperation = 1:numel(operationIds)
    medians = arrayfun(@(strategy)strategy.operations(iOperation).medianSeconds,strategyResults);
    valid = arrayfun(@(strategy)strategy.operations(iOperation).correctnessPassed,strategyResults);
    selection = waveVortexTransformLayoutPreference(strategyIds(valid),medians(valid));
    selections(end+1) = struct("operationId",operationIds(iOperation),"strictFastestStrategy",selection.strictFastestStrategy,"strictFastestSeconds",selection.strictFastestMedianSeconds,"preferredStrategy",selection.preferredStrategy,"preferredSeconds",selection.preferredMedianSeconds,"currentSeconds",selection.currentMedianSeconds,"currentRelativeToFastest",selection.currentRelativeToFastest,"tieTolerance",selection.tieTolerance,"currentRetained",selection.currentRetained); %#ok<AGROW>
end

caseResult = struct( ...
    "id",benchmarkCase.id, ...
    "transformId",benchmarkCase.transformId, ...
    "scoreFamily",benchmarkCase.scoreFamily, ...
    "operation",benchmarkCase.operation, ...
    "Lxyz",benchmarkCase.Lxyz, ...
    "Nxyz",benchmarkCase.Nxyz, ...
    "isHydrostatic",benchmarkCase.isHydrostatic, ...
    "shouldAntialias",benchmarkCase.shouldAntialias, ...
    "seed",benchmarkCase.seed, ...
    "warmupCount",benchmarkCase.warmupCount, ...
    "sampleCount",benchmarkCase.sampleCount, ...
    "status","complete", ...
    "failure",emptyFailure(), ...
    "reference",struct("strategyId","wv-sorted-linear","forwardRelativeError",relativeInfinityError(wvReference,reshape(fullSpectrum(geometry.dftPrimaryIndex),Nz,geometry.Nkl)),"inputProjectionRelativeError",relativeInfinityError(productionInverse,realInput),"fullBufferRelativeError",relativeInfinityError(productionTransform.complexBuffer,fullBufferReference)), ...
    "storage",commonStorageLedger(realInput,fullSpectrum,wvReference,productionInverse), ...
    "warmupSchedules",warmupSchedules, ...
    "sampleSchedules",sampleSchedules, ...
    "strategies",strategyResults, ...
    "selections",selections);
end

function runtimes = createStrategyRuntimes(geometry,Nz,strategyIds)
Nx = geometry.Nx;
Ny = geometry.Ny;
currentPrimary = uint64(geometry.dftPrimaryIndex);
currentConjugate = uint64(geometry.dftConjugateIndex);
currentWvConjugate = uint64(geometry.wvConjugateIndex);
[dftPrimary,primaryOrder] = sort(currentPrimary);
[dftConjugate,conjugateOrder] = sort(currentConjugate);
wvPrimary = uint64(primaryOrder);
wvConjugate = currentWvConjugate(conjugateOrder);
primaryRows = uint64(geometry.dftPrimaryIndices2D);
wvConjugateRows = uint64(find(geometry.lMode_wv == 0 & geometry.kMode_wv ~= 0));
conjugateRows = uint64(geometry.dftConjugateIndices2D(wvConjugateRows));
planeOffsets = uint64((0:Nz-1)*Nx*Ny);

runtimeCells = cell(1,numel(strategyIds));
for iStrategy = 1:numel(strategyIds)
    id = strategyIds(iStrategy);
    switch id
        case "wv-sorted-linear"
            maps = struct("dftPrimary",currentPrimary,"dftConjugate",currentConjugate,"wvConjugate",currentWvConjugate);
            buffer = WVTransformLayoutBenchmarkBuffer([Nx Ny Nz],"full-3d");
            description = "Current WV-sorted linear gather and primary-then-conjugate insertion.";
            expressionId = "full-wv-linear-v1";
        case "dft-sorted-linear"
            maps = struct("dftPrimary",dftPrimary,"wvPrimary",wvPrimary,"dftConjugate",dftConjugate,"wvConjugate",wvConjugate);
            buffer = WVTransformLayoutBenchmarkBuffer([Nx Ny Nz],"full-3d");
            description = "DFT-sorted replicated linear indices with explicit WV permutations.";
            expressionId = "full-dft-linear-v1";
        case "two-dimensional-rows"
            maps = struct("primaryRows",primaryRows,"conjugateRows",conjugateRows,"wvConjugateRows",wvConjugateRows);
            buffer = WVTransformLayoutBenchmarkBuffer([Nx Ny Nz],"rows-2d");
            description = "Two-dimensional DFT rows with measured transpose between row and canonical layouts.";
            expressionId = "full-rows-2d-v1";
        case "per-plane-linear"
            maps = struct("primaryRows",primaryRows,"conjugateRows",conjugateRows,"wvConjugateRows",wvConjugateRows,"planeOffsets",planeOffsets);
            buffer = WVTransformLayoutBenchmarkBuffer([Nx Ny Nz],"full-3d");
            description = "Two-dimensional indices plus per-plane offsets and complete MATLAB loops.";
            expressionId = "full-per-plane-v1";
        otherwise
            error("WaveVortexBenchmark:UnknownLayoutStrategy","Unknown transform-layout strategy %s.",id);
    end
    runtimeCells{iStrategy} = struct("id",id,"expressionId",expressionId,"expressionDescription",description,"maps",maps,"mappingLedger",structLedger(maps),"buffer",buffer,"bufferLedger",arrayLedger("persistentFullBuffer",buffer.value));
end
runtimes = [runtimeCells{:}];
end

function validation = validateStrategy(runtime,fullSpectrum,wvReference,fullBufferReference,realInput,productionInverse,tolerance)
runtime.buffer.reset();
extracted = executeOperation(runtime,"extract",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
runtime.buffer.reset();
executeOperation(runtime,"insert-primary",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
primaryError = selectedBufferError(runtime,fullBufferReference,"primary");
runtime.buffer.reset();
executeOperation(runtime,"insert-conjugate",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
conjugateError = selectedBufferError(runtime,fullBufferReference,"conjugate");
runtime.buffer.reset();
executeOperation(runtime,"insert-complete",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
combinedError = relativeInfinityError(fullBuffer(runtime),fullBufferReference);
forwardOutput = executeOperation(runtime,"forward-complete",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
inverseOutput = executeOperation(runtime,"inverse-complete",fullSpectrum,wvReference,realInput,size(realInput,1),size(realInput,2));
errors = struct("extract",relativeInfinityError(extracted,wvReference),"insertPrimary",primaryError,"insertConjugate",conjugateError,"insertComplete",combinedError,"forwardComplete",relativeInfinityError(forwardOutput,wvReference),"inverseComplete",relativeInfinityError(inverseOutput,productionInverse));
values = cell2mat(struct2cell(errors));
if any(values > tolerance)
    error("WaveVortexBenchmark:LayoutCorrectnessFailed","Strategy %s exceeded the relative-error tolerance (maximum %.3g).",runtime.id,max(values));
end
validation = struct("errors",errors,"passed",true);
end

function output = executeOperation(runtime,operationId,fullSpectrum,wvInput,realInput,Nx,Ny)
switch operationId
    case "extract"
        output = extractWV(runtime,fullSpectrum);
    case "insert-primary"
        insertPrimary(runtime,wvInput);
        output = [];
    case "insert-conjugate"
        insertConjugate(runtime,wvInput);
        output = [];
    case "insert-complete"
        insertPrimary(runtime,wvInput);
        insertConjugate(runtime,wvInput);
        output = [];
    case "forward-complete"
        spectrum = fft(fft(realInput,Nx,1),Ny,2)/(Nx*Ny);
        output = extractWV(runtime,spectrum);
    case "inverse-complete"
        insertPrimary(runtime,wvInput);
        insertConjugate(runtime,wvInput);
        output = ifft(ifft(fullBuffer(runtime),Nx,1),Ny,2,"symmetric")*(Nx*Ny);
    otherwise
        error("WaveVortexBenchmark:UnknownLayoutOperation","Unknown transform-layout operation %s.",operationId);
end
end

function output = extractWV(runtime,fullSpectrum)
Nz = size(fullSpectrum,3);
switch runtime.id
    case "wv-sorted-linear"
        output = reshape(fullSpectrum(runtime.maps.dftPrimary),Nz,[]);
    case "dft-sorted-linear"
        output = complex(zeros(Nz,numel(runtime.maps.wvPrimary)/Nz));
        output(runtime.maps.wvPrimary) = fullSpectrum(runtime.maps.dftPrimary);
    case "two-dimensional-rows"
        rows = reshape(fullSpectrum,[],Nz);
        output = rows(runtime.maps.primaryRows,:).';
    case "per-plane-linear"
        output = complex(zeros(Nz,numel(runtime.maps.primaryRows)));
        for iZ = 1:Nz
            output(iZ,:) = fullSpectrum(runtime.maps.primaryRows+runtime.maps.planeOffsets(iZ));
        end
end
end

function insertPrimary(runtime,wvInput)
switch runtime.id
    case "wv-sorted-linear"
        runtime.buffer.value(runtime.maps.dftPrimary) = wvInput;
    case "dft-sorted-linear"
        runtime.buffer.value(runtime.maps.dftPrimary) = wvInput(runtime.maps.wvPrimary);
    case "two-dimensional-rows"
        runtime.buffer.value(runtime.maps.primaryRows,:) = wvInput.';
    case "per-plane-linear"
        for iZ = 1:size(wvInput,1)
            runtime.buffer.value(runtime.maps.primaryRows+runtime.maps.planeOffsets(iZ)) = wvInput(iZ,:);
        end
end
end

function insertConjugate(runtime,wvInput)
switch runtime.id
    case {"wv-sorted-linear","dft-sorted-linear"}
        runtime.buffer.value(runtime.maps.dftConjugate) = conj(wvInput(runtime.maps.wvConjugate));
    case "two-dimensional-rows"
        runtime.buffer.value(runtime.maps.conjugateRows,:) = conj(wvInput(:,runtime.maps.wvConjugateRows).');
    case "per-plane-linear"
        for iZ = 1:size(wvInput,1)
            runtime.buffer.value(runtime.maps.conjugateRows+runtime.maps.planeOffsets(iZ)) = conj(wvInput(iZ,runtime.maps.wvConjugateRows));
        end
end
end

function buffer = fullBuffer(runtime)
buffer = runtime.buffer.value;
if runtime.buffer.storageShape == "rows-2d"
    buffer = reshape(buffer,runtime.buffer.physicalSize);
end
end

function errorValue = selectedBufferError(runtime,reference,part)
buffer = fullBuffer(runtime);
switch part
    case "primary"
        if isfield(runtime.maps,"dftPrimary")
            indices = runtime.maps.dftPrimary;
        else
            indices = expandedIndices(runtime.maps.primaryRows,size(reference));
        end
    case "conjugate"
        if isfield(runtime.maps,"dftConjugate")
            indices = runtime.maps.dftConjugate;
        else
            indices = expandedIndices(runtime.maps.conjugateRows,size(reference));
        end
end
errorValue = relativeInfinityError(buffer(indices),reference(indices));
end

function indices = expandedIndices(rows,physicalSize)
offsets = uint64((0:physicalSize(3)-1)*physicalSize(1)*physicalSize(2));
indices = rows(:)+offsets;
indices = indices(:);
end

function pairs = allOperationPairs(strategyCount,operationCount)
[strategies,operations] = ndgrid(1:strategyCount,1:operationCount);
pairs = [strategies(:) operations(:)];
end

function order = rotatedOrder(count,roundNumber)
startIndex = mod(roundNumber-1,count)+1;
order = [startIndex:count 1:startIndex-1];
end

function ids = pairIds(pairs,strategyIds,operationIds)
ids = strings(1,size(pairs,1));
for iPair = 1:size(pairs,1)
    ids(iPair) = strategyIds(pairs(iPair,1)) + "/" + operationIds(pairs(iPair,2));
end
end

function errorValue = validationError(validation,operationId)
field = operationField(operationId);
errorValue = validation.errors.(field);
end

function passed = validationPassed(validation,operationId)
passed = validation.passed && isfinite(validationError(validation,operationId));
end

function field = operationField(operationId)
switch operationId
    case "extract"
        field = "extract";
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
end
end

function bytes = operationResultBytes(operationId,Nx,Ny,Nz,Nkl)
switch operationId
    case {"extract","forward-complete"}
        bytes = 16*Nz*Nkl;
    case "inverse-complete"
        bytes = 8*Nx*Ny*Nz;
    otherwise
        bytes = 0;
end
end

function ledger = operationResultLedger(operationId,Nx,Ny,Nz,Nkl)
switch operationId
    case {"extract","forward-complete"}
        ledger = struct("name","timedWVResult","class","double","shape",[Nz Nkl],"isComplex",true,"bytes",operationResultBytes(operationId,Nx,Ny,Nz,Nkl));
    case "inverse-complete"
        ledger = struct("name","timedRealResult","class","double","shape",[Nx Ny Nz],"isComplex",false,"bytes",operationResultBytes(operationId,Nx,Ny,Nz,Nkl));
    otherwise
        ledger = struct("name","none","class","","shape",[0 0],"isComplex",false,"bytes",0);
end
end

function ledger = commonStorageLedger(realInput,fullSpectrum,wvArray,realOutput)
ledger = [arrayLedger("realInput",realInput) arrayLedger("fullSpectrum",fullSpectrum) arrayLedger("wvArray",wvArray) arrayLedger("realOutput",realOutput)];
end

function ledger = structLedger(values)
names = string(fieldnames(values));
ledger = repmat(arrayLedger("",zeros(0,1)),1,numel(names));
for iValue = 1:numel(names)
    ledger(iValue) = arrayLedger(names(iValue),values.(names(iValue)));
end
end

function ledger = arrayLedger(name,value)
info = whos("value");
ledger = struct("name",string(name),"class",string(class(value)),"shape",double(size(value)),"bytes",double(info.bytes));
end

function errorValue = relativeInfinityError(actual,expected)
denominator = max(max(abs(expected(:))),realmin("double"));
errorValue = max(abs(actual(:)-expected(:)))/denominator;
end

function records = productionSourceRecords(repositoryRoot)
sources = [ ...
    struct("id","forward","path","FastTransforms/@WVFastTransformDoublyPeriodicMatlab/transformFromSpatialDomainWithFourier.m"), ...
    struct("id","inverse","path","FastTransforms/@WVFastTransformDoublyPeriodicMatlab/transformToSpatialDomainWithFourier.m"), ...
    struct("id","geometry-index","path","@WVGeometryDoublyPeriodic/WVGeometryDoublyPeriodic.m")];
records = repmat(struct("id","","path","","sha256",""),1,numel(sources));
for iSource = 1:numel(sources)
    pathname = fullfile(repositoryRoot,sources(iSource).path);
    records(iSource) = struct("id",sources(iSource).id,"path",sources(iSource).path,"sha256",sha256File(pathname));
end
end

function hash = sha256File(pathname)
fileId = fopen(pathname,"r");
if fileId < 0
    error("WaveVortexBenchmark:MissingProductionSource","Unable to read production source %s.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
bytes = fread(fileId,Inf,"*uint8");
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(bytes);
hashBytes = typecast(digest.digest(),"uint8");
hash = lower(string(reshape(dec2hex(hashBytes,2).',1,[])));
clear cleanup
end

function caseResult = failedCaseResult(benchmarkCase,exception)
caseResult = struct("id",benchmarkCase.id,"transformId",benchmarkCase.transformId,"scoreFamily",benchmarkCase.scoreFamily,"operation",benchmarkCase.operation,"Lxyz",benchmarkCase.Lxyz,"Nxyz",benchmarkCase.Nxyz,"isHydrostatic",benchmarkCase.isHydrostatic,"shouldAntialias",benchmarkCase.shouldAntialias,"seed",benchmarkCase.seed,"warmupCount",benchmarkCase.warmupCount,"sampleCount",benchmarkCase.sampleCount,"status","failed","failure",exceptionFailure(exception),"reference",struct,"storage",struct([]),"warmupSchedules",strings(0,0),"sampleSchedules",strings(0,0),"strategies",emptyStrategyResults(),"selections",emptySelections());
end

function validation = emptyValidation()
validation = struct("errors",struct("extract",NaN,"insertPrimary",NaN,"insertConjugate",NaN,"insertComplete",NaN,"forwardComplete",NaN,"inverseComplete",NaN),"passed",false);
end

function failure = exceptionFailure(exception)
failure = struct("identifier",string(exception.identifier),"message",string(exception.message),"stack",string({exception.stack.name}));
end

function failure = emptyFailure()
failure = struct("identifier","","message","","stack",strings(0,1));
end

function results = emptyCaseResults()
results = struct("id",{},"transformId",{},"scoreFamily",{},"operation",{},"Lxyz",{},"Nxyz",{},"isHydrostatic",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"sampleCount",{},"status",{},"failure",{},"reference",{},"storage",{},"warmupSchedules",{},"sampleSchedules",{},"strategies",{},"selections",{});
end

function results = emptyStrategyResults()
results = struct("id",{},"expressionId",{},"expressionDescription",{},"assignmentSteps",{},"mappingArrays",{},"mappingBytes",{},"buffer",{},"persistentBufferReused",{},"timedBufferClearing",{},"sourceArraysUnchanged",{},"copyObservation",{},"copyObservationReason",{},"operations",{});
end

function results = emptyOperationResults()
results = struct("id",{},"rawSeconds",{},"medianSeconds",{},"relativeError",{},"correctnessPassed",{},"resultBytes",{},"resultStorage",{});
end

function results = emptySelections()
results = struct("operationId",{},"strictFastestStrategy",{},"strictFastestSeconds",{},"preferredStrategy",{},"preferredSeconds",{},"currentSeconds",{},"currentRelativeToFastest",{},"tieTolerance",{},"currentRetained",{});
end

function results = emptyScores()
results = struct("id",{},"backendId",{},"score",{});
end
