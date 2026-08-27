function threeInterfaceMatlabWorker(configPath,outputPath)
% Execute one matched MATLAB-interface benchmark sample in a fresh process.
config = jsondecode(fileread(configPath));
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
phasePath = string(config.phasePath);
result = failedResult(config);
model = [];
wvt = [];
try
    path(config.matlabPath);
    addRepositoryPaths(config.repositoryRoot,config.benchmarkFolder);
    writePhase(phasePath,"startup");
    interfaceTimer = tic;
    backend = string(config.backend);
    benchmarkCase = string(config.case.id);
    requestedIntegrator = string(config.case.requestedIntegrator);
    if benchmarkCase == "nonlinear-flux"
        [wvt,reader] = WVTransform.waveVortexTransformFromFile(config.inputPath,iTime=Inf,shouldReadOnly=true,computationalBackend=backend);
        reader.close();
        validateBackend(wvt,backend);
        provider = providerRecord(wvt,backend,config);
        writePhase(phasePath,"steady-retained");
        pause(config.plateauSeconds);
        writePhase(phasePath,"integrate");
        waitForPhaseSample("integrate");
        operationTimer = tic;
        [Fp,Fm,F0] = wvt.nonlinearFlux();
        integrationSeconds = toc(operationTimer);
        writePhase(phasePath,"integration-complete");
        writeComplexBinary(config.comparisonPath,Fp,Fm,F0);
        finalState = stateRecord(wvt);
        outputAgreement = struct("kind","flux-binary","path",string(config.comparisonPath));
        actualIntegrator = "none";
        integrator = struct("requested",requestedIntegrator,"actual",actualIntegrator,"matched",requestedIntegrator==actualIntegrator);
    else
        model = benchmarkModelFromFile(config.inputPath,backend);
        wvt = model.wvt;
        validateBackend(wvt,backend);
        provider = providerRecord(wvt,backend,config);
        if requestedIntegrator == "fixed-rk4"
            model.setupIntegrator(integratorType="fixed",deltaT=config.case.deltaT);
        elseif startsWith(requestedIntegrator,"adaptive-")
            matlabIntegrator = matlabIntegratorFor(requestedIntegrator);
            model.setupIntegrator(integratorType="adaptive",integrator=matlabIntegrator,relTolerance=config.case.relativeTolerance,absTolerance=config.case.absoluteTolerance,shouldShowIntegrationStats=0);
            model.odeOptions = odeset(model.odeOptions,'InitialStep',config.case.initialStep,'Stats','on');
        else
            error("WaveVortexBenchmark:UnknownIntegrator","Unsupported requested integrator %s.",requestedIntegrator);
        end
        actualIntegrator = activeIntegrator(model);
        integrator = struct("requested",requestedIntegrator,"actual",actualIntegrator,"matched",requestedIntegrator==actualIntegrator);
        integrator.initialTime = double(model.t);
        integrator.finalTime = double(config.case.finalTime);
        if startsWith(requestedIntegrator,"adaptive-")
            integrator = adaptiveIntegratorRecordBeforeRun(model,config.case,integrator);
        end
        writePhase(phasePath,"steady-retained");
        pause(config.plateauSeconds);
        rhsEvaluationsBefore = double(model.nFluxComputations);
        writePhase(phasePath,"integrate");
        waitForPhaseSample("integrate");
        operationTimer = tic;
        if startsWith(requestedIntegrator,"adaptive-")
            statisticsText = evalc("runIntegration(model,config.case.finalTime)");
        else
            runIntegration(model,config.case.finalTime);
            statisticsText = "";
        end
        integrationSeconds = toc(operationTimer);
        writePhase(phasePath,"integration-complete");
        rhsEvaluationCount = double(model.nFluxComputations)-rhsEvaluationsBefore;
        if startsWith(requestedIntegrator,"adaptive-")
            integrator = completeAdaptiveIntegratorRecord(integrator,statisticsText,rhsEvaluationCount,requestedIntegrator);
        else
            integrator.acceptedStepCount = round((config.case.finalTime-integrator.initialTime)/config.case.deltaT);
            integrator.rejectedStepCount = 0;
            integrator.rhsEvaluationCount = rhsEvaluationCount;
            integrator.fsalReuseCount = 0;
            integrator.fsalInvalidationCount = 0;
        end
        finalState = stateRecord(wvt);
        model.closeNetCDFFile();
        integrator.outputRecordCounts = outputRecordCounts(config.inputPath);
        integrator.denseOutputEvaluationCount = scheduledInteriorOutputCount(integrator,config.case);
        integrator.storageAccounting = matlabIntegratorStorageAccounting;
        outputAgreement = struct("kind","model-output","path",string(config.inputPath));
    end
    interfaceTotalSeconds = toc(interfaceTimer);
    writePhase(phasePath,"outputs-held");
    pause(config.plateauSeconds);
    if ~isempty(model) && isvalid(model)
        delete(model);
        model = [];
    elseif ~isempty(wvt) && isvalid(wvt)
        delete(wvt);
        wvt = [];
    end
    writePhase(phasePath,"complete");
    result = struct( ...
        "schemaVersion","three-interface-worker-v1", ...
        "status","complete", ...
        "interface",string(config.interface), ...
        "case",config.case, ...
        "sourceCommit",string(config.sourceCommit), ...
        "timing",struct("interfaceTotalSeconds",interfaceTotalSeconds,"integrationSeconds",integrationSeconds), ...
        "memory",struct(), ...
        "provider",provider, ...
        "integrator",integrator, ...
        "finalState",finalState, ...
        "output",outputAgreement, ...
        "failure",emptyFailure);
catch exception
    if ~isempty(model) && isvalid(model)
        try
            model.closeNetCDFFile();
        catch
        end
        delete(model);
    elseif ~isempty(wvt) && isvalid(wvt)
        delete(wvt);
    end
    writePhase(phasePath,"failed");
    result.failure = struct("identifier",string(exception.identifier),"message",string(exception.message),"report",string(getReport(exception,"extended","hyperlinks","off")));
end
writeText(outputPath,jsonencode(result,PrettyPrint=true));
clear stateCleanup
end

function runIntegration(model,finalTime)
model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
end

function value = adaptiveIntegratorRecordBeforeRun(model,definition,value)
value.controller = "matlab-"+string(func2str(model.odeIntegrator))+"-v1";
value.relativeTolerance = double(odeget(model.odeOptions,'RelTol'));
tolerances = odeget(model.odeOptions,'AbsTol');
value.absoluteToleranceHash = threeInterfaceToleranceHash(tolerances);
value.absoluteToleranceHashClearedMantissaBits = 20;
value.absoluteToleranceComponentHashes = arrayfun(@(index)threeInterfaceToleranceHash(tolerances(model.arrayStartIndex(index):model.arrayEndIndex(index))),1:numel(model.arrayStartIndex));
value.requestedInitialStep = double(definition.initialStep);
value.effectiveInitialStep = double(odeget(model.odeOptions,'InitialStep'));
value.requestedMaximumStep = double(definition.maximumStep);
value.maximumStepPolicy = string(definition.maximumStepPolicy);
value.effectiveMaximumStep = 0.1*(double(definition.finalTime)-double(model.t));
value.initialTime = double(model.t);
value.finalTime = double(definition.finalTime);
end

function value = completeAdaptiveIntegratorRecord(value,statisticsText,rhsEvaluationCount,requestedIntegrator)
value.acceptedStepCount = statisticCount(statisticsText,"successful steps");
value.rejectedStepCount = statisticCount(statisticsText,"failed attempts");
reportedEvaluations = statisticCount(statisticsText,"function evaluations");
if reportedEvaluations ~= rhsEvaluationCount
    error("WaveVortexBenchmark:AdaptiveWorkCountMismatch","MATLAB %s reported %d function evaluations, but WVModel counted %d.",matlabIntegratorName(requestedIntegrator),reportedEvaluations,rhsEvaluationCount);
end
value.rhsEvaluationCount = rhsEvaluationCount;
value.fsalReuseCount = "not-exposed-by-matlab";
value.fsalInvalidationCount = "not-exposed-by-matlab";
end

function value = statisticCount(text,label)
token = regexp(text,'(?m)^\s*(\d+)\s+'+label,"tokens","once");
if isempty(token)
    error("WaveVortexBenchmark:MissingAdaptiveStatistics","MATLAB adaptive solver did not report %s.",label);
end
value = str2double(token{1});
end

function value = outputRecordCounts(pathname)
value = struct("waveVortex",recordCount(pathname,"/wave-vortex/t"),"dense",recordCount(pathname,"/dense/t"),"particles",recordCount(pathname,"/particles/t"),"tracers",recordCount(pathname,"/tracers/t"));
end

function identifier = activeIntegrator(model)
if string(model.integratorType) == "fixed"
    identifier = "fixed-rk4";
elseif string(model.integratorType) == "adaptive" && ismember(string(func2str(model.odeIntegrator)),["ode23" "ode45" "ode78"])
    identifier = "adaptive-rk"+extractAfter(string(func2str(model.odeIntegrator)),"ode");
else
    identifier = string(model.integratorType)+":"+string(func2str(model.odeIntegrator));
end
end

function integrator = matlabIntegratorFor(identifier)
switch string(identifier)
    case "adaptive-rk23"
        integrator = @ode23;
    case "adaptive-rk45"
        integrator = @ode45;
    case "adaptive-rk78"
        integrator = @ode78;
    otherwise
        error("WaveVortexBenchmark:UnknownIntegrator","Unsupported requested integrator %s.",identifier);
end
end

function name = matlabIntegratorName(identifier)
name = "ode"+extractAfter(string(identifier),"adaptive-rk");
end

function value = scheduledInteriorOutputCount(integrator,definition)
value = 0;
if ~isfield(definition,"workload") || string(definition.workload) ~= "composite-dense-output"
    return
end
if string(definition.requestedIntegrator) == "fixed-rk4"
    value = double(definition.denseOutputPointsPerStep);
else
    value = double(definition.denseOutputIntegrationRecordCount);
end
end

function value = matlabIntegratorStorageAccounting
value = struct("scope","MATLAB solver, allocator, and copy-on-write storage are opaque","exact",false,"persistentBytes","not-exposed-by-matlab","workspaceMaximumLiveBytes","not-exposed-by-matlab","denseHistoryBytes","not-exposed-by-matlab","stateEquivalentBufferCount","not-exposed-by-matlab");
end

function value = recordCount(pathname,variable)
try
    value = numel(ncread(pathname,variable));
catch exception
    if contains(string(exception.identifier),"netcdf") || contains(lower(string(exception.message)),"variable") || contains(lower(string(exception.message)),"group")
        value = 0;
    else
        rethrow(exception)
    end
end
end

function model = benchmarkModelFromFile(pathname,backend)
[wvt,reader] = WVTransform.waveVortexTransformFromFile(pathname,iTime=Inf,shouldReadOnly=true,computationalBackend=backend);
readerCleanup = onCleanup(@()closeReader(reader));
isDynamicsLinear = false;
if isKey(reader.attributes,'WVModelIsDynamicsLinear')
    isDynamicsLinear = logical(reader.attributes('WVModelIsDynamicsLinear'));
end
model = WVModel(wvt,shouldUseLinearDynamics=isDynamicsLinear);
reader.close();
clear readerCleanup
ncfile = NetCDFFile(char(pathname),shouldReadOnly=false);
try
    outputFile = WVModelOutputFile.modelOutputFileFromFile(ncfile,model);
    model.addOutputFile(outputFile);
catch exception
    if ~isempty(ncfile.id)
        ncfile.close();
    end
    rethrow(exception)
end
end

function validateBackend(wvt,backend)
metadata = wvt.computationalBackendMetadata;
if string(metadata.activeBackend) ~= backend
    error("WaveVortexBenchmark:InterfaceFallback","Requested %s but %s executed.",backend,string(metadata.activeBackend));
end
if backend == "compiled" && (string(metadata.provider.id) ~= "native-neon-pthreads" || ~metadata.module.identityValidated || metadata.libraries.openmp.detected || metadata.contract.planCount ~= 17)
    error("WaveVortexBenchmark:InterfaceProvider","The MATLAB compiled interface did not execute the validated native provider.");
end
end

function value = providerRecord(wvt,backend,config)
if backend == "compiled"
    metadata = wvt.computationalBackendMetadata;
    value = struct("id",string(metadata.provider.id),"version",string(metadata.provider.version),"threads",double(config.threadCount),"baseLibrary",string(metadata.libraries.base.path),"threadLibrary",string(metadata.libraries.thread.path),"noFallback",true);
else
    value = struct("id","matlab-builtin","version",string(version),"threads",double(config.threadCount),"baseLibrary","","threadLibrary","","noFallback",true);
end
end

function value = stateRecord(wvt)
value = struct("t",wvt.t,"shape",[wvt.Nj wvt.Nkl],"ApNorm",norm(wvt.Ap(:)),"AmNorm",norm(wvt.Am(:)),"A0Norm",norm(wvt.A0(:)));
end

function writeComplexBinary(pathname,varargin)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ComparisonWrite","Unable to open %s.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
for iValue = 1:numel(varargin)
    values = varargin{iValue}(:);
    interleaved = reshape([real(values).'; imag(values).'],[],1);
    fwrite(fileId,interleaved,"double");
end
clear cleanup
end

function writePhase(pathname,phase)
temporary = pathname+".tmp";
writeText(temporary,phase);
movefile(temporary,pathname,"f");
end

function waitForPhaseSample(expectedPhase)
acknowledgementPath = string(getenv("WV_RSS_PHASE_ACK"));
if acknowledgementPath == ""
    return
end
timer = tic;
while toc(timer) < 5
    if isfile(acknowledgementPath) && strtrim(string(fileread(acknowledgementPath))) == expectedPhase
        return
    end
    pause(0.001)
end
error("WaveVortexBenchmark:RSSPhaseHandshake","The external RSS sampler did not acknowledge phase %s.",expectedPhase)
end

function addRepositoryPaths(repositoryRoot,benchmarkFolder)
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
for iFolder = 1:numel(metadata.folders)
    folder = fullfile(repositoryRoot,metadata.folders(iFolder).path);
    if isfolder(folder)
        addpath(folder);
    end
end
addpath(repositoryRoot);
addpath(benchmarkFolder);
end

function closeReader(reader)
if ~isempty(reader) && isvalid(reader) && ~isempty(reader.id)
    reader.close();
end
end

function writeText(pathname,contents)
fileId = fopen(pathname,"w");
if fileId < 0
    error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to open %s.",pathname);
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",contents);
clear cleanup
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

function value = emptyFailure
value = struct("identifier","","message","","report","");
end

function value = failedResult(config)
value = struct("schemaVersion","three-interface-worker-v1","status","failed","interface",string(config.interface),"case",config.case,"sourceCommit",string(config.sourceCommit),"timing",struct(),"memory",struct(),"provider",struct(),"integrator",struct(),"finalState",struct(),"output",struct(),"failure",emptyFailure);
end
