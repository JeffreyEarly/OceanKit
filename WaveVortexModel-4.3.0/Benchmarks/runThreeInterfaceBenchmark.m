function results = runThreeInterfaceBenchmark(options)
% Compare matched MATLAB builtin, MATLAB compiled, and standalone compiled work.
arguments
    options.Nxyz (1,3) double {mustBeInteger,mustBePositive} = [256 256 129]
    options.Lxyz (1,3) double {mustBePositive} = [150e3 150e3 1300]
    options.processRunCount (1,1) double {mustBeInteger,mustBePositive} = 3
    options.deltaT (1,1) double = NaN
    options.adaptiveInitialStep (1,1) double = NaN
    options.relativeTolerance (1,1) double {mustBePositive} = 1e-3
    options.absoluteTolerance (1,1) double {mustBePositive} = 1e-6
    options.caseIds (1,:) string = strings(1,0)
    options.integrators (1,:) string {mustBeMember(options.integrators,["fixed-rk4" "adaptive-rk23" "adaptive-rk45" "adaptive-rk78"])} = ["fixed-rk4" "adaptive-rk23" "adaptive-rk45" "adaptive-rk78"]
    options.workloads (1,:) string {mustBeMember(options.workloads,["coefficient-endpoint" "composite-dense-output"])} = ["coefficient-endpoint" "composite-dense-output"]
    options.physicalConfigurations (1,:) string {mustBeMember(options.physicalConfigurations,["hydrostatic" "nonhydrostatic"])} = "nonhydrostatic"
    options.integrationStepCount (1,1) double {mustBeInteger,mustBePositive} = 56
    options.denseOutputPointsPerStep (1,1) double {mustBeInteger,mustBePositive} = 3
    options.adaptiveStepCount (1,1) double {mustBeInteger,mustBePositive} = 10
    options.adaptiveOutputCount (1,1) double {mustBeInteger,mustBePositive} = 2
    options.samplingIntervalSeconds (1,1) double {mustBePositive} = 0.005
    options.plateauSeconds (1,1) double {mustBePositive} = 0.05
    options.outputDirectory (1,1) string = ""
    options.archiveDirectory (1,1) string = ""
    options.runId (1,1) string = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"))
    options.shouldWriteArtifacts (1,1) logical = true
    options.injectWorkerFailure (1,1) logical = false
end
if ~ismac || string(computer("arch")) ~= "maca64" || isMATLABReleaseOlderThan("R2025b")
    error("WaveVortexBenchmark:ThreeInterfaceUnsupportedPlatform","The three-interface benchmark requires MATLAB R2025b or later on Apple silicon.");
end
legacyCaseIds = ["nonlinear-flux" "fixed-rk4-continuation" "adaptive-rk23-observer-output"];
if any(~ismember(options.caseIds,legacyCaseIds))
    error("WaveVortexBenchmark:UnknownCase","Unknown legacy three-interface case: %s",strjoin(options.caseIds(~ismember(options.caseIds,legacyCaseIds)),", "));
end
if numel(options.caseIds) > 1 && any(options.caseIds == "adaptive-rk23-observer-output")
    error("WaveVortexBenchmark:MixedOutputSchedules","The adaptive benchmark uses its own output schedule and must run separately from the frozen nonlinear-flux and fixed-RK4 cases.");
end
repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
benchmarkFolder = fullfile(repositoryRoot,"Benchmarks");
originalDirectory = pwd;
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng));
addRepositoryPaths(repositoryRoot,benchmarkFolder);
[options,physicalState,physicalEvidence,stepControlEvidence] = resolveBenchmarkContract(options);
if options.outputDirectory == ""
    options.outputDirectory = fullfile(benchmarkFolder,"results","runs",options.runId+"-three-interface");
end
if options.archiveDirectory == ""
    options.archiveDirectory = fullfile(fileparts(repositoryRoot),"wave-vortex-model-benchmark-artifacts","three-interface");
end
if options.shouldWriteArtifacts && isfolder(options.outputDirectory)
    error("WaveVortexBenchmark:ThreeInterfaceOutputExists","Output already exists: %s",options.outputDirectory);
end
if options.shouldWriteArtifacts
    mkdir(options.outputDirectory);
end
results = initializeResult(options,repositoryRoot,physicalEvidence,stepControlEvidence);
workFolder = string(tempname);
mkdir(workFolder);
workCleanup = onCleanup(@()removeFolder(workFolder));
verifyTemporaryCapacity(workFolder,options);
activeStage = "provider";
try
    capabilities = WVCompiledBackend.capabilities();
    if ~capabilities.isAvailable
        error("WaveVortexBenchmark:ValidatedProviderCacheRequired","The three-interface benchmark requires the prevalidated native FFTW provider cache; provider selection or construction is not performed by the benchmark.");
    end
    validateCapabilities(capabilities);
    results.provider = capabilities;
    activeStage = "build";
    executables = buildStandaloneWorkers(repositoryRoot,capabilities);
    activeStage = "fixture";
    definitions = caseDefinitions(options);
    [fixtures,fixtureRecords] = createMatchedFixtures(workFolder,options,definitions,physicalState);
    clear physicalState
    results.configuration.fixtures = fixtureRecords;
    if numel(fixtureRecords) == 1
        results.configuration.fixtureSHA256 = fixtureRecords.sha256;
    end
    checkpoint(results,options);

    activeStage = "workers";
    interfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
    results.cases = definitions;
    for iRepeat = 1:options.processRunCount
        caseOrder = mod((0:numel(definitions)-1)+(iRepeat-1),numel(definitions))+1;
        for iCase = caseOrder
            interfaceOrder = mod((0:2)+(iRepeat+iCase-2),3)+1;
            for iInterface = interfaceOrder
                interface = interfaces(iInterface);
                fprintf("Three-interface benchmark: %s, %s, repeat %d/%d.\n",interface,definitions(iCase).id,iRepeat,options.processRunCount);
                fixturePath = fixtures(fixtureKey(definitions(iCase)));
                results.runs(end+1,1) = runOne(interface,definitions(iCase),iRepeat,fixturePath,executables,capabilities,options,repositoryRoot,benchmarkFolder,workFolder);
                checkpoint(results,options);
                if string(results.runs(end).status) ~= "complete"
                    error("WaveVortexBenchmark:ThreeInterfaceWorker","%s/%s failed: %s",interface,definitions(iCase).id,results.runs(end).failure.message);
                end
            end
        end
        activeStage = "correctness-repeat";
        repeatMask = [results.runs.repeatIndex] == iRepeat;
        repeatComparison = aggregate(results.runs(repeatMask),definitions,1e-12);
        [results.runs(repeatMask),releasedBytes] = releaseValidatedOutputs(results.runs(repeatMask),workFolder);
        results.repeatComparisonEvidence(end+1,1) = struct("repeatIndex",iRepeat,"comparison",repeatComparison,"releasedRunCount",nnz(repeatMask),"releasedBytes",releasedBytes); %#ok<AGROW>
        checkpoint(results,options);
        activeStage = "workers";
    end
    failures = results.runs(string({results.runs.status})~="complete");
    if ~isempty(failures)
        messages = arrayfun(@(item)string(item.interface)+"/"+string(item.case.id)+": "+string(item.failure.message),failures);
        error("WaveVortexBenchmark:ThreeInterfaceWorkers","One or more workers failed:%s%s",newline,strjoin(messages,newline));
    end
    activeStage = "correctness";
    results.comparison = aggregate(results.runs,definitions,1e-12,results.repeatComparisonEvidence);
    validateThreeInterfaceBenchmarkContract(results);
    results.status = "complete";
    results.completedAtUTC = utcTimestamp;
    results.failure = emptyFailure;
    writeArtifacts(results,options);
    if options.shouldWriteArtifacts
        activeStage = "archive";
        results.externalArchive = archiveDetailedArtifact(options);
    end
catch exception
    results.status = "failed";
    results.completedAtUTC = utcTimestamp;
    results.failure = struct("stage",activeStage,"identifier",string(exception.identifier),"message",string(exception.message),"report",string(getReport(exception,"extended","hyperlinks","off")));
    writeArtifacts(results,options);
    rethrow(exception)
end
clear workCleanup stateCleanup
end

function [options,state,physicalEvidence,stepControlEvidence] = resolveBenchmarkContract(options)
state = struct();
physicalEvidence = struct();
stepControlEvidence = struct();
if ~isnan(options.deltaT) && options.deltaT <= 0
    error("WaveVortexBenchmark:InvalidStepControl","deltaT must be positive when supplied.")
end
if ~isnan(options.adaptiveInitialStep) && options.adaptiveInitialStep <= 0
    error("WaveVortexBenchmark:InvalidStepControl","adaptiveInitialStep must be positive when supplied.")
end
if ~isempty(options.caseIds)
    if isnan(options.deltaT), options.deltaT = 1/1024; end
    if isnan(options.adaptiveInitialStep), options.adaptiveInitialStep = options.deltaT; end
    return
end
if ~isequal(options.physicalConfigurations,"nonhydrostatic")
    error("WaveVortexBenchmark:PhysicalConfiguration","The canonical integrator study uses one nonhydrostatic physical configuration.")
end
wvt = WVTransformConstantStratification(options.Lxyz,options.Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=true);
[state,physicalEvidence] = initializeThreeInterfaceIntegratorState(wvt,4001);
model = WVModel(wvt);
cleanup = onCleanup(@()closeModels(model));
[fixedCFLCandidate,advectiveCFLCandidate,oscillatoryCFLCandidate] = model.timeStepForCFL(0.25);
[adaptiveCFLCandidate,adaptiveAdvectiveCFLCandidate,adaptiveOscillatoryCFLCandidate] = model.timeStepForCFL(0.5);
if isnan(options.deltaT)
    options.deltaT = 2^floor(log2(fixedCFLCandidate));
end
if options.deltaT > fixedCFLCandidate
    error("WaveVortexBenchmark:UnstableFixedStep","The requested RK4 step %.17g s exceeds the frozen-state CFL candidate %.17g s.",options.deltaT,fixedCFLCandidate)
end
if isnan(options.adaptiveInitialStep)
    options.adaptiveInitialStep = adaptiveCFLCandidate;
end
stepControlEvidence = struct( ...
    "fixedStepPolicy","largest power of two not exceeding the frozen-state CFL=0.25 candidate unless explicitly overridden", ...
    "fixedStep",options.deltaT, ...
    "fixedCFLCandidate",fixedCFLCandidate, ...
    "fixedAdvectiveCFLCandidate",advectiveCFLCandidate, ...
    "fixedOscillatoryCFLCandidate",oscillatoryCFLCandidate, ...
    "adaptiveInitialStepPolicy","frozen-state CFL=0.5 candidate unless explicitly overridden", ...
    "adaptiveInitialStep",options.adaptiveInitialStep, ...
    "adaptiveCFLCandidate",adaptiveCFLCandidate, ...
    "adaptiveAdvectiveCFLCandidate",adaptiveAdvectiveCFLCandidate, ...
    "adaptiveOscillatoryCFLCandidate",adaptiveOscillatoryCFLCandidate, ...
    "adaptiveMaximumStepPolicy","MATLAB default; no user maximum step", ...
    "finalTime",options.integrationStepCount*options.deltaT);
clear cleanup
end

function results = initializeResult(options,repositoryRoot,physicalEvidence,stepControlEvidence)
[commit,tree,isDirty] = gitIdentity(repositoryRoot);
isLegacy = ~isempty(options.caseIds);
results = struct( ...
    "schemaVersion",conditional(isLegacy,"three-interface-benchmark-v1","three-interface-benchmark-v2"), ...
    "status","running", ...
    "runId",options.runId, ...
    "generatedAtUTC",utcTimestamp, ...
    "completedAtUTC","", ...
    "source",struct("repository","JeffreyEarly/wave-vortex-model","commit",commit,"tree",tree,"isDirty",isDirty), ...
    "environment",environmentRecord, ...
    "configuration",struct("studyId",conditional(isLegacy,"legacy-three-interface-v1","integrator-runtime-memory-v1"),"Nxyz",options.Nxyz,"Lxyz",options.Lxyz,"processRunCount",options.processRunCount,"warmupCount",0,"samplesPerProcess",1,"deltaT",options.deltaT,"adaptiveInitialStep",options.adaptiveInitialStep,"maximumStepPolicy",conditional(isLegacy,"explicit","matlab-default"),"relativeTolerance",options.relativeTolerance,"absoluteTolerance",options.absoluteTolerance,"caseIds",options.caseIds,"integrators",options.integrators,"workloads",options.workloads,"physicalConfigurations",options.physicalConfigurations,"integrationStepCount",options.integrationStepCount,"denseOutputPointsPerStep",options.denseOutputPointsPerStep,"adaptiveStepCount",options.adaptiveStepCount,"adaptiveOutputCount",options.adaptiveOutputCount,"threadCount",min(18,maxNumCompThreads),"samplingIntervalSeconds",options.samplingIntervalSeconds,"fixtureSHA256","","fixtures",struct([]),"initialCondition",physicalEvidence,"stepControls",stepControlEvidence,"correctnessTolerance",1e-12,"outputAgreementTolerancePolicy",conditional(isLegacy,"legacy transform-relative tolerance","matched method RelTol plus base component AbsTol"),"timingBoundary","integration-only wall time begins immediately before integration and ends after required output delivery; startup, construction, provider creation, planning, and cleanup are excluded","rssBoundary","external process-tree RSS sampled only while the worker reports integrate or output-delivery phases; total live-process RSS peak is primary, steady-retained baseline, operation increment, process-lifetime peak, and final RSS are diagnostics","outputRetentionPolicy","complete all cross-interface and endpoint/dense correctness gates for one repeat, retain repeat evidence, then release validated raw NetCDF outputs"), ...
    "provider",struct(),"cases",[],"runs",repmat(emptyRun,0,1),"repeatComparisonEvidence",repmat(struct("repeatIndex",0,"comparison",[],"releasedRunCount",0,"releasedBytes",0),0,1),"comparison",[],"failure",emptyFailure);
end

function definitions = caseDefinitions(options)
if isempty(options.caseIds)
    finalTime = options.integrationStepCount*options.deltaT;
    denseOutputStartTime = 0;
    denseOutputInterval = options.deltaT/(options.denseOutputPointsPerStep+1);
    common = struct("Nxyz",options.Nxyz,"Lxyz",options.Lxyz,"deltaT",options.deltaT,"finalTime",finalTime,"relativeTolerance",options.relativeTolerance,"absoluteTolerance",options.absoluteTolerance,"outputRelativeTolerance",options.relativeTolerance,"outputAbsoluteTolerance",options.absoluteTolerance,"initialStep",options.adaptiveInitialStep,"maximumStep",[],"maximumStepPolicy","matlab-default","forcing","default WVNonlinearAdvection","shouldAntialias",true,"seed",4001,"operation","model-continuation","integrationStepCount",options.integrationStepCount,"denseOutputPointsPerStep",options.denseOutputPointsPerStep,"denseOutputStartTime",[],"denseOutputEndTime",[],"denseOutputRecordCount",0,"denseOutputIntegrationRecordCount",0);
    definitions = repmat(mergeStruct(common,struct("id","","requestedIntegrator","","physicalConfiguration","","isHydrostatic",false,"workload","","outputInterval",NaN,"observerGraph","")),0,1);
    for physicalConfiguration = options.physicalConfigurations
        for integrator = options.integrators
            for workload = options.workloads
                isDense = workload == "composite-dense-output";
                outputInterval = conditional(isDense,denseOutputInterval,finalTime);
                observerGraph = conditional(isDense,"one persisted restart record plus four first-step records containing u, two 3-D particles with tracked u, one 3-D tracer, and one source-linked mooring; coefficients at endpoints","coefficient-only state with endpoint-only output");
                identifier = physicalConfiguration+"--"+integrator+"--"+workload;
                schedule = struct("denseOutputStartTime",conditional(isDense,denseOutputStartTime,[]),"denseOutputEndTime",conditional(isDense,options.deltaT,[]),"denseOutputRecordCount",conditional(isDense,options.denseOutputPointsPerStep+2,0),"denseOutputIntegrationRecordCount",conditional(isDense,options.denseOutputPointsPerStep+1,0));
                definitions(end+1,1) = mergeStruct(common,mergeStruct(struct("id",identifier,"requestedIntegrator",integrator,"physicalConfiguration",physicalConfiguration,"isHydrostatic",physicalConfiguration=="hydrostatic","workload",workload,"outputInterval",outputInterval,"observerGraph",observerGraph),schedule)); %#ok<AGROW>
            end
        end
    end
    return
end
fixtureOutputInterval = options.deltaT/2;
if options.caseIds == "adaptive-rk23-observer-output"
    fixtureOutputInterval = options.adaptiveStepCount*options.deltaT/options.adaptiveOutputCount;
end
common = struct("Nxyz",options.Nxyz,"deltaT",options.deltaT,"finalTime",2*options.deltaT,"relativeTolerance",options.relativeTolerance,"absoluteTolerance",options.absoluteTolerance,"initialStep",options.deltaT,"maximumStep",options.deltaT,"maximumStepPolicy","explicit","outputInterval",fixtureOutputInterval,"observerGraph","cross-group Eulerian u, one mooring, two 3-D particles with tracked u, and one 3-D tracer","forcing","default WVNonlinearAdvection","shouldAntialias",true,"seed",4001);
available = [ ...
    mergeStruct(common,struct("id","nonlinear-flux","operation","nonlinearFlux","requestedIntegrator","none")); ...
    mergeStruct(common,struct("id","fixed-rk4-continuation","operation","model-continuation","requestedIntegrator","fixed-rk4")); ...
    mergeStruct(common,struct("id","adaptive-rk23-observer-output","operation","model-continuation","requestedIntegrator","adaptive-rk23","finalTime",(1+options.adaptiveStepCount)*options.deltaT))];
definitions = available(ismember(string({available.id}),options.caseIds));
end

function run = runOne(interface,definition,repeatIndex,fixturePath,executables,capabilities,options,repositoryRoot,benchmarkFolder,workFolder)
sampleFolder = fullfile(workFolder,sprintf('%s-%s-%d',interface,definition.id,repeatIndex));
mkdir(sampleFolder);
inputPath = fullfile(sampleFolder,"model.nc");
copyfile(fixturePath,inputPath);
comparisonPath = fullfile(sampleFolder,"comparison.bin");
phasePath = fullfile(sampleFolder,"phase.txt");
samplePath = fullfile(sampleFolder,"rss.tsv");
stdoutPath = fullfile(sampleFolder,"stdout.txt");
stderrPath = fullfile(sampleFolder,"stderr.txt");
writeText(phasePath,"launch");
processTimer = tic;
if startsWith(interface,"matlab-")
    backend = extractAfter(interface,"matlab-");
    if backend == "builtin", backend = "matlab"; end
    config = struct("interface",interface,"backend",backend,"sourceCommit",gitValue(repositoryRoot,"rev-parse HEAD"),"case",definition,"inputPath",inputPath,"comparisonPath",comparisonPath,"threadCount",min(18,maxNumCompThreads),"repositoryRoot",repositoryRoot,"benchmarkFolder",benchmarkFolder,"matlabPath",path,"phasePath",phasePath,"plateauSeconds",options.plateauSeconds);
    configPath = fullfile(sampleFolder,"config.json");
    outputPath = fullfile(sampleFolder,"worker.json");
    writeText(configPath,jsonencode(config));
    if options.injectWorkerFailure
        configPath = configPath+".missing";
    end
    statement = "addpath('"+replace(benchmarkFolder,"'","''")+"'); threeInterfaceMatlabWorker('"+replace(configPath,"'","''")+"','"+replace(outputPath,"'","''")+"')";
    workerCommand = sprintf('"%s" -batch "%s"',fullfile(matlabroot,"bin","matlab"),replace(statement,'"','\"'));
    command = sampledCommand(workerCommand,samplePath,phasePath,stdoutPath,stderrPath,options,benchmarkFolder);
    [exitCode,~] = system(command);
    processWallSeconds = toc(processTimer);
    commandOutput = readCommandOutput(stdoutPath,stderrPath);
    if exitCode ~= 0 || ~isfile(outputPath)
        run = emptyRun;
        run.interface = interface; run.case = definition; run.repeatIndex = repeatIndex;
        run.failure = struct("identifier","WaveVortexBenchmark:MatlabInterfaceWorkerFailed","message",string(commandOutput),"report",string(commandOutput));
        return
    end
    value = jsondecode(fileread(outputPath));
    run = normalizeMatlabRun(value,repeatIndex,processWallSeconds,processMemory(samplePath,options.samplingIntervalSeconds));
else
    if definition.id == "nonlinear-flux"
        workerCommand = shellQuote(executables.kernel)+" "+shellQuote(inputPath)+" "+string(min(18,maxNumCompThreads))+" 0 1 "+shellQuote(comparisonPath)+" --phase-file "+shellQuote(phasePath);
    else
        workerCommand = shellQuote(executables.runner)+" "+shellQuote(inputPath)+" --restart-mode model --output-policy append --delta-t "+numberText(definition.deltaT)+" --final-time "+numberText(definition.finalTime)+" --fft-provider native-fftw --threads "+string(min(18,maxNumCompThreads))+" --phase-file "+shellQuote(phasePath);
        workerCommand = workerCommand+" --integrator "+definition.requestedIntegrator;
        if isfield(definition,"workload") && string(definition.workload)=="composite-dense-output"
            workerCommand = workerCommand+" --benchmark-model-dense-output 1";
        end
        if startsWith(definition.requestedIntegrator,"adaptive-")
            workerCommand = workerCommand+" --relative-tolerance "+numberText(definition.relativeTolerance)+" --absolute-tolerance "+numberText(definition.absoluteTolerance)+" --initial-step "+numberText(definition.initialStep);
            if ~isempty(definition.maximumStep)
                workerCommand = workerCommand+" --maximum-step "+numberText(definition.maximumStep);
            end
        end
    end
    command = sampledCommand(workerCommand,samplePath,phasePath,stdoutPath,stderrPath,options,benchmarkFolder);
    [exitCode,~] = system(command);
    processWallSeconds = toc(processTimer);
    commandOutput = readCommandOutput(stdoutPath,stderrPath);
    if exitCode ~= 0
        run = emptyRun; run.interface = interface; run.case = definition; run.repeatIndex = repeatIndex;
        run.failure = struct("identifier","WaveVortexBenchmark:StandaloneInterfaceWorkerFailed","message",string(commandOutput),"report",string(commandOutput));
        return
    end
    value = jsondecode(fileread(stdoutPath));
    run = normalizeStandaloneRun(value,definition,repeatIndex,processWallSeconds,inputPath,comparisonPath,capabilities,processMemory(samplePath,options.samplingIntervalSeconds));
end
end

function run = normalizeMatlabRun(value,repeatIndex,processWallSeconds,memory)
run = emptyRun;
run.schemaVersion = string(value.schemaVersion); run.status = string(value.status); run.interface = string(value.interface); run.case = value.case; run.repeatIndex = repeatIndex; run.sourceCommit = string(value.sourceCommit); run.processWallSeconds = processWallSeconds; run.failure = value.failure;
if run.status ~= "complete"
    return
end
run.integrationSeconds = value.timing.integrationSeconds; run.interfaceTotalSeconds = value.timing.interfaceTotalSeconds; run.memory = memory; run.provider = value.provider; run.integrator = value.integrator; run.finalState = value.finalState; run.output = value.output;
end

function run = normalizeStandaloneRun(value,definition,repeatIndex,processWallSeconds,inputPath,comparisonPath,capabilities,memory)
run = emptyRun;
sourceCommit = "";
if isfield(value,"sourceCommit"), sourceCommit = string(value.sourceCommit); elseif isfield(value,"source") && isfield(value.source,"commit"), sourceCommit = string(value.source.commit); end
run.schemaVersion = string(value.schemaVersion); run.status = string(value.status); run.interface = "standalone-compiled"; run.case = definition; run.repeatIndex = repeatIndex; run.sourceCommit = sourceCommit; run.processWallSeconds = processWallSeconds;
if definition.id == "nonlinear-flux"
    run.integrationSeconds = value.timing.medianSeconds; run.interfaceTotalSeconds = processWallSeconds;
    run.memory = memory;
    run.integrator = struct("requested","none","actual","none","matched",true); run.finalState = struct(); run.output = struct("kind","flux-binary","path",comparisonPath);
else
    run.integrationSeconds = value.timingSeconds.integrate; run.interfaceTotalSeconds = value.timingSeconds.total;
    run.memory = memory;
    run.integrator = standaloneIntegratorRecord(value,definition,inputPath); run.finalState = value.state; run.output = struct("kind","model-output","path",inputPath);
    if isfield(value,"integrationBreakdownSeconds"), run.diagnostics = value.integrationBreakdownSeconds; end
end
noFallback = true;
if isfield(value.execution,"noFallback"), noFallback = logical(value.execution.noFallback); end
run.provider = struct("id","native-neon-pthreads","version",string(capabilities.provider.version),"threads",double(capabilities.contract.threadCount),"baseLibrary",string(value.provider.baseLibrary),"threadLibrary",string(value.provider.threadLibrary),"noFallback",noFallback);
run.failure = emptyFailure;
end

function comparison = aggregate(runs,definitions,tolerance,repeatEvidence)
if nargin < 4
    repeatEvidence = struct([]);
end
usesRepeatEvidence = ~isempty(repeatEvidence);
interfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
comparison = repmat(struct("id","","interfaces",[],"maximumRelativeError",NaN,"outputAgreementPassed",false,"outputGraph",emptyOutputGraph,"integratorAgreementPassed",false,"adaptiveWorkAgreementPassed",false,"absoluteToleranceFingerprintAgreementPassed",true,"memoryAgreementPassed",false,"endpointTrajectoryAgreementPassed",true,"matchedContractPassed",false),numel(definitions),1);
for iCase = 1:numel(definitions)
    selected = runs(string(arrayfun(@(item)item.case.id,runs,"UniformOutput",false))==definitions(iCase).id);
    records = repmat(struct("id","","processWallSeconds",NaN,"interfaceTotalSeconds",NaN,"integrationSeconds",NaN,"totalPeakRSSBytes",NaN,"incrementalPeakRSSBytes",NaN,"finalRSSBytes",NaN,"processWallRatio",NaN,"integrationRatio",NaN,"totalRSSRatio",NaN,"incrementalRSSRatio",NaN),numel(interfaces),1);
    builtin = selected(string({selected.interface})=="matlab-builtin");
    memoryPassed = all(arrayfun(@(item)string(item.memory.status)=="complete" && isfinite(item.memory.totalPeakRSSBytes) && item.memory.totalPeakRSSBytes>0,selected));
    builtinProcess = median([builtin.processWallSeconds]); builtinIntegration = median([builtin.integrationSeconds]); builtinPeak = median(arrayfun(@(item)item.memory.totalPeakRSSBytes,builtin)); builtinIncrement = median(arrayfun(@(item)item.memory.peakIncrementBytes,builtin));
    maximumError = 0; outputPassed = true; outputGraph = emptyOutputGraph; endpointTrajectoryPassed = true;
    if usesRepeatEvidence
        caseEvidence = arrayfun(@(item)item.comparison(iCase),repeatEvidence);
        maximumError = max([caseEvidence.maximumRelativeError]);
        outputPassed = all([caseEvidence.outputAgreementPassed]);
        endpointTrajectoryPassed = all([caseEvidence.endpointTrajectoryAgreementPassed]);
        for evidence = reshape(caseEvidence,1,[])
            outputGraph = mergeOutputGraph(outputGraph,evidence.outputGraph);
        end
    end
    for iInterface = 1:numel(interfaces)
        candidate = selected(string({selected.interface})==interfaces(iInterface));
        processValue = median([candidate.processWallSeconds]); integrationValue = median([candidate.integrationSeconds]); peakValue = median(arrayfun(@(item)item.memory.totalPeakRSSBytes,candidate)); incrementValue = median(arrayfun(@(item)item.memory.peakIncrementBytes,candidate)); finalValue = median(arrayfun(@(item)item.memory.finalRSSBytes,candidate));
        records(iInterface) = struct("id",interfaces(iInterface),"processWallSeconds",processValue,"interfaceTotalSeconds",median([candidate.interfaceTotalSeconds]),"integrationSeconds",integrationValue,"totalPeakRSSBytes",peakValue,"incrementalPeakRSSBytes",incrementValue,"finalRSSBytes",finalValue,"processWallRatio",processValue/builtinProcess,"integrationRatio",integrationValue/builtinIntegration,"totalRSSRatio",peakValue/builtinPeak,"incrementalRSSRatio",safeRatio(incrementValue,builtinIncrement));
        if iInterface > 1 && ~usesRepeatEvidence
            for iRepeat = 1:numel(candidate)
                reference = builtin([builtin.repeatIndex]==candidate(iRepeat).repeatIndex);
                [relativeTolerance,absoluteTolerance] = outputAgreementTolerances(definitions(iCase),tolerance);
                [errorValue,agreement,details] = compareOutputs(reference,candidate(iRepeat),relativeTolerance,absoluteTolerance);
                maximumError = max(maximumError,errorValue); outputPassed = outputPassed && agreement;
                outputGraph = mergeOutputGraph(outputGraph,details);
            end
        end
    end
    providersPassed = all(arrayfun(@(item)item.provider.noFallback,selected)) && all(arrayfun(@(item)item.interface=="matlab-builtin" || string(item.provider.id)=="native-neon-pthreads",selected));
    integratorsPassed = all(arrayfun(@(item)logical(item.integrator.matched) && string(item.integrator.requested)==definitions(iCase).requestedIntegrator && string(item.integrator.actual)==definitions(iCase).requestedIntegrator,selected));
    isAdaptive = startsWith(definitions(iCase).requestedIntegrator,"adaptive-");
    if isfield(definitions(iCase),"workload")
        adaptiveWorkPassed = ~isAdaptive || adaptiveControlsMatch(selected,definitions(iCase));
    else
        adaptiveWorkPassed = definitions(iCase).requestedIntegrator ~= "adaptive-rk23" || adaptiveWorkMatches(selected,definitions(iCase));
    end
    toleranceFingerprintPassed = ~isAdaptive || toleranceFingerprintMatches(selected);
    toleranceGatePassed = ~isfield(definitions(iCase),"workload") || toleranceFingerprintPassed;
    strictRelativeGate = ~isfield(definitions(iCase),"outputRelativeTolerance") && maximumError<=tolerance;
    methodToleranceGate = isfield(definitions(iCase),"outputRelativeTolerance") && outputPassed;
    numericalGatePassed = strictRelativeGate || methodToleranceGate;
    comparison(iCase) = struct("id",definitions(iCase).id,"interfaces",records,"maximumRelativeError",maximumError,"outputAgreementPassed",outputPassed,"outputGraph",outputGraph,"integratorAgreementPassed",integratorsPassed,"adaptiveWorkAgreementPassed",adaptiveWorkPassed,"absoluteToleranceFingerprintAgreementPassed",toleranceFingerprintPassed,"memoryAgreementPassed",memoryPassed,"endpointTrajectoryAgreementPassed",endpointTrajectoryPassed,"matchedContractPassed",providersPassed&&integratorsPassed&&adaptiveWorkPassed&&toleranceGatePassed&&memoryPassed&&endpointTrajectoryPassed&&numericalGatePassed);
end
if ~usesRepeatEvidence && all(arrayfun(@(definition)isfield(definition,"workload"),definitions))
    comparison = validateEndpointTrajectoryIndependence(comparison,runs,definitions,tolerance);
end
end

function [runs,releasedBytes] = releaseValidatedOutputs(runs,workFolder)
releasedBytes = 0;
for iRun = 1:numel(runs)
    pathname = string(runs(iRun).output.path);
    runs(iRun).output.retention = "released-after-repeat-correctness";
    runs(iRun).output.releasedBytes = 0;
    if pathname == "" || ~isfile(pathname)
        continue
    end
    if ~startsWith(pathname,workFolder+filesep)
        error("WaveVortexBenchmark:UnsafeOutputRelease","Refusing to release a validated output outside the benchmark temporary directory: %s",pathname);
    end
    information = dir(pathname);
    runs(iRun).output.releasedBytes = information.bytes;
    releasedBytes = releasedBytes+information.bytes;
    delete(pathname);
end
end

function passed = adaptiveControlsMatch(runs,definition)
required = ["controller" "relativeTolerance" "requestedInitialStep" "effectiveInitialStep" "requestedMaximumStep" "maximumStepPolicy" "effectiveMaximumStep" "initialTime" "finalTime" "acceptedStepCount" "rejectedStepCount" "rhsEvaluationCount" "denseOutputEvaluationCount" "fsalReuseCount" "fsalInvalidationCount" "outputRecordCounts"];
passed = all(arrayfun(@(run)all(isfield(run.integrator,required)),runs));
if ~passed, return, end
expectedController = "matlab-ode"+extractAfter(string(definition.requestedIntegrator),"adaptive-rk")+"-v1";
for run = reshape(runs,1,[])
    current = run.integrator;
    passed = passed && string(current.controller)==expectedController;
    passed = passed && numericControlMatches(current.relativeTolerance,definition.relativeTolerance);
    passed = passed && numericControlMatches(current.requestedInitialStep,definition.initialStep) && numericControlMatches(current.effectiveInitialStep,definition.initialStep);
    passed = passed && string(current.maximumStepPolicy)==string(definition.maximumStepPolicy);
    passed = passed && optionalNumericControlMatches(current.requestedMaximumStep,definition.maximumStep);
    passed = passed && numericControlMatches(current.effectiveMaximumStep,0.1*(double(current.finalTime)-double(current.initialTime)));
    passed = passed && numericControlMatches(current.finalTime,definition.finalTime);
end
end

function passed = numericControlMatches(actual,expected)
passed = isscalar(actual) && isscalar(expected) && abs(double(actual)-double(expected))<=8*eps(max([1 abs(double(expected))]));
end

function passed = optionalNumericControlMatches(actual,expected)
if isempty(actual) || isempty(expected)
    passed = isempty(actual) && isempty(expected);
else
    passed = numericControlMatches(actual,expected);
end
end

function comparison = validateEndpointTrajectoryIndependence(comparison,runs,definitions,tolerance)
physicalConfigurations = unique(string({definitions.physicalConfiguration}),"stable");
integrators = unique(string({definitions.requestedIntegrator}),"stable");
interfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
for physicalConfiguration = physicalConfigurations
    for integrator = integrators
        endpointIndex = find(string({definitions.physicalConfiguration})==physicalConfiguration & string({definitions.requestedIntegrator})==integrator & string({definitions.workload})=="coefficient-endpoint",1);
        denseIndex = find(string({definitions.physicalConfiguration})==physicalConfiguration & string({definitions.requestedIntegrator})==integrator & string({definitions.workload})=="composite-dense-output",1);
        if isempty(endpointIndex) || isempty(denseIndex), continue, end
        passed = true;
        for interface = interfaces
            runInterfaces = reshape(string({runs.interface}),[],1);
            runCaseIds = reshape(string(arrayfun(@(run)run.case.id,runs,"UniformOutput",false)),[],1);
            endpointRuns = runs(runInterfaces==interface & runCaseIds==definitions(endpointIndex).id);
            denseRuns = runs(runInterfaces==interface & runCaseIds==definitions(denseIndex).id);
            for iRepeat = 1:numel(endpointRuns)
                candidate = denseRuns([denseRuns.repeatIndex]==endpointRuns(iRepeat).repeatIndex);
                if numel(candidate)~=1 || string(endpointRuns(iRepeat).output.kind)~="model-output" || string(candidate.output.kind)~="model-output"
                    passed = false;
                    continue
                end
                [trajectoryPassed,~] = compareFinalCoefficientState(endpointRuns(iRepeat).output.path,candidate.output.path,tolerance);
                passed = passed && trajectoryPassed;
            end
        end
        comparison(endpointIndex).endpointTrajectoryAgreementPassed = passed;
        comparison(denseIndex).endpointTrajectoryAgreementPassed = passed;
        comparison(endpointIndex).matchedContractPassed = comparison(endpointIndex).matchedContractPassed && passed;
        comparison(denseIndex).matchedContractPassed = comparison(denseIndex).matchedContractPassed && passed;
    end
end
end

function [passed,maximumRelativeError] = compareFinalCoefficientState(referencePath,candidatePath,tolerance)
reference = ncinfo(referencePath,"/wave-vortex");
candidate = ncinfo(candidatePath,"/wave-vortex");
coefficientNames = ["Ap" "Am" "A0" "Ap_real" "Ap_imag" "Am_real" "Am_imag" "A0_real" "A0_imag"];
referenceNames = string({reference.Variables.Name});
candidateNames = string({candidate.Variables.Name});
names = intersect(referenceNames(ismember(referenceNames,coefficientNames)),candidateNames(ismember(candidateNames,coefficientNames)),"stable");
passed = ~isempty(names);
maximumRelativeError = 0;
for name = names
    referenceInfo = reference.Variables(referenceNames==name);
    candidateInfo = candidate.Variables(candidateNames==name);
    expected = readFinalRecord(referencePath,"/wave-vortex/"+name,referenceInfo);
    actual = readFinalRecord(candidatePath,"/wave-vortex/"+name,candidateInfo);
    if ~isequal(size(expected),size(actual))
        passed = false;
        maximumRelativeError = Inf;
        continue
    end
    difference = max(abs(double(actual)-double(expected)),[],"all");
    scale = max(abs(double(expected)),[],"all");
    relativeError = difference/max(scale,realmin("double"));
    maximumRelativeError = max(maximumRelativeError,relativeError);
    passed = passed && isfinite(relativeError) && difference<=1e-13+tolerance*scale;
end
referenceTime = ncread(referencePath,"/wave-vortex/t");
candidateTime = ncread(candidatePath,"/wave-vortex/t");
passed = passed && numericControlMatches(referenceTime(end),candidateTime(end));
end

function value = readFinalRecord(pathname,variablePath,information)
shape = double(information.Size);
dimensionNames = string({information.Dimensions.Name});
timeDimension = find(dimensionNames=="t",1);
if isempty(timeDimension)
    value = ncread(pathname,variablePath);
    return
end
start = ones(1,numel(shape));
count = shape;
start(timeDimension) = shape(timeDimension);
count(timeDimension) = 1;
value = ncread(pathname,variablePath,start,count);
end

function value = standaloneIntegratorRecord(report,definition,pathname)
value = struct("requested",definition.requestedIntegrator,"actual",string(report.integrator.id),"matched",definition.requestedIntegrator==string(report.integrator.id));
value.initialTime = double(report.state.initialTime);
value.finalTime = double(report.state.finalTime);
value.acceptedStepCount = double(report.state.stepCount);
value.rejectedStepCount = double(report.state.rejectedStepCount);
value.rhsEvaluationCount = double(report.state.rhsEvaluationCount);
value.denseOutputEvaluationCount = double(report.integrator.denseOutputEvaluationCount);
value.denseOutputCacheBuildCount = double(report.integrator.denseOutputCacheBuildCount);
value.denseOutputCacheReuseCount = double(report.integrator.denseOutputCacheReuseCount);
value.continuousExtensionRightHandSideEvaluationCount = double(report.integrator.continuousExtensionRightHandSideEvaluationCount);
value.fsalReuseCount = double(report.integrator.fsalReuseCount);
value.fsalInvalidationCount = double(report.integrator.fsalInvalidationCount);
value.constraintModifiedCoefficientCount = double(report.integrator.constraintModifiedCoefficientCount);
value.workspaceStateEquivalentCount = double(report.integrator.workspaceStateEquivalentCount);
value.workspaceMaximumLiveStateEquivalentCount = double(report.integrator.workspaceMaximumLiveStateEquivalentCount);
value.denseHistoryStateEquivalentCount = double(report.integrator.denseHistoryStateEquivalentCount);
value.continuousExtensionWorkspaceStateEquivalentCount = double(report.integrator.continuousExtensionWorkspaceStateEquivalentCount);
value.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount = double(report.integrator.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount);
value.stateSizedBuffers = report.integratorStateSizedBuffers;
value.storageAccounting = report.integratorStorageLedger;
value.sharedAbstractionStateSizedCopyCount = double(report.livenessBytes.contractAbstractionAdditionalArrayStorage>0);
if startsWith(definition.requestedIntegrator,"adaptive-")
    value.controller = string(report.integrator.controller);
    value.relativeTolerance = double(report.integrator.relativeTolerance);
    value.absoluteToleranceHash = string(report.integrator.toleranceHash);
    value.absoluteToleranceHashClearedMantissaBits = double(report.integrator.toleranceHashClearedMantissaBits);
    value.absoluteToleranceComponentHashes = string(report.integrator.toleranceComponentHashes(:)');
    value.requestedInitialStep = double(report.integrator.requestedInitialStep);
    value.effectiveInitialStep = double(report.integrator.effectiveInitialStep);
    value.requestedMaximumStep = double(report.integrator.requestedMaximumStep);
    value.maximumStepPolicy = string(definition.maximumStepPolicy);
    value.effectiveMaximumStep = double(report.integrator.effectiveMaximumStep);
end
value.outputRecordCounts = outputRecordCounts(pathname);
end

function passed = adaptiveWorkMatches(runs,definition)
required = ["controller" "relativeTolerance" "absoluteToleranceHash" "absoluteToleranceHashClearedMantissaBits" "absoluteToleranceComponentHashes" "requestedInitialStep" "effectiveInitialStep" "requestedMaximumStep" "effectiveMaximumStep" "initialTime" "finalTime" "acceptedStepCount" "rejectedStepCount" "rhsEvaluationCount" "denseOutputEvaluationCount" "outputRecordCounts"];
passed = all(arrayfun(@(run)all(isfield(run.integrator,required)),runs));
if ~passed
    return
end
reference = runs(1).integrator;
for run = reshape(runs,1,[])
    current = run.integrator;
    % The hashes remain useful diagnostics, but cannot be an equality gate:
    % independently evaluated MATLAB/C++ formulas can straddle an arbitrary
    % quantization boundary despite agreeing to roundoff. Controller inputs,
    % accepted work, output schedules, and complete numerical output are the
    % publication contract.
    exactFields = ["controller" "absoluteToleranceHashClearedMantissaBits" "acceptedStepCount" "rejectedStepCount" "rhsEvaluationCount" "denseOutputEvaluationCount"];
    passed = passed && all(arrayfun(@(name)string(current.(name))==string(reference.(name)),exactFields));
    numericFields = ["relativeTolerance" "requestedInitialStep" "effectiveInitialStep" "requestedMaximumStep" "effectiveMaximumStep" "initialTime" "finalTime"];
    for name = numericFields
        if name=="requestedMaximumStep"
            passed = passed && optionalNumericControlMatches(current.(name),reference.(name));
        else
            passed = passed && numericControlMatches(current.(name),reference.(name));
        end
    end
    passed = passed && isequal(orderfields(current.outputRecordCounts),orderfields(reference.outputRecordCounts));
end
passed = passed && reference.controller == "matlab-ode23-v1" && reference.relativeTolerance == definition.relativeTolerance && reference.requestedInitialStep == definition.initialStep && reference.effectiveInitialStep == definition.initialStep && reference.requestedMaximumStep == definition.maximumStep && reference.effectiveMaximumStep == definition.maximumStep && reference.finalTime == definition.finalTime;
end

function passed = toleranceFingerprintMatches(runs)
reference = runs(1).integrator;
passed = all(arrayfun(@(run)string(run.integrator.absoluteToleranceHash)==string(reference.absoluteToleranceHash),runs));
passed = passed && all(arrayfun(@(run)isequal(string(run.integrator.absoluteToleranceComponentHashes(:)),string(reference.absoluteToleranceComponentHashes(:))),runs));
end

function [errorValue,agreement,details] = compareOutputs(reference,candidate,relativeTolerance,absoluteTolerance)
if string(reference.output.kind) == "flux-binary"
    expected = readBinary(reference.output.path); actual = readBinary(candidate.output.path);
    errorValue = max(abs(actual-expected))/max(max(abs(expected)),realmin("double")); agreement = isfinite(errorValue) && errorValue <= relativeTolerance;
    details = struct("kind","flux-arrays","passed",agreement,"maximumRelativeError",errorValue,"maximumAbsoluteError",max(abs(actual-expected)),"variableCount",3,"recordCount",0,"categories",struct("name","coefficients","variableCount",3,"maximumAbsoluteError",max(abs(actual-expected)),"maximumRelativeError",errorValue,"passed",agreement),"differences",strings(0,1));
    return
end
details = compareWaveVortexOutputGraphs(reference.output.path,candidate.output.path,relativeTolerance=relativeTolerance,absoluteTolerance=absoluteTolerance);
details.kind = "complete-netcdf-output-graph";
errorValue = details.maximumRelativeError;
agreement = details.passed;
end

function [relativeTolerance,absoluteTolerance] = outputAgreementTolerances(definition,legacyTolerance)
relativeTolerance = legacyTolerance;
absoluteTolerance = 0;
if isfield(definition,"outputRelativeTolerance")
    relativeTolerance = double(definition.outputRelativeTolerance);
    absoluteTolerance = double(definition.outputAbsoluteTolerance);
end
end

function result = mergeOutputGraph(result,value)
if result.kind == "", result = value; return, end
result.passed = result.passed && value.passed;
result.maximumRelativeError = max(result.maximumRelativeError,value.maximumRelativeError);
result.maximumAbsoluteError = max(result.maximumAbsoluteError,value.maximumAbsoluteError);
result.variableCount = max(result.variableCount,value.variableCount);
result.recordCount = max(result.recordCount,value.recordCount);
result.differences = unique([result.differences(:); value.differences(:)],"stable");
for category = reshape(value.categories,1,[])
    index = find(string({result.categories.name})==string(category.name),1);
    if isempty(index)
        result.categories(end+1) = category;
    else
        result.categories(index).variableCount = max(result.categories(index).variableCount,category.variableCount);
        result.categories(index).maximumAbsoluteError = max(result.categories(index).maximumAbsoluteError,category.maximumAbsoluteError);
        result.categories(index).maximumRelativeError = max(result.categories(index).maximumRelativeError,category.maximumRelativeError);
        result.categories(index).passed = result.categories(index).passed && category.passed;
    end
end
end

function [fixtures,records] = createMatchedFixtures(workFolder,options,definitions,physicalState)
keys = unique(arrayfun(@fixtureKey,definitions),"stable");
fixtures = dictionary(string.empty(0,1),string.empty(0,1));
records = repmat(struct("key","","path","","sha256","","physicalConfiguration","","workload",""),numel(keys),1);
for iKey = 1:numel(keys)
    definition = definitions(find(arrayfun(@fixtureKey,definitions)==keys(iKey),1));
    pathname = fullfile(workFolder,"matched-"+replace(keys(iKey),"--","-")+"-model.nc");
    createMatchedFixture(pathname,options,definition,physicalState);
    fixtures(keys(iKey)) = pathname;
    physicalConfiguration = "legacy-nonhydrostatic";
    workload = "legacy-composite";
    if isfield(definition,"physicalConfiguration"), physicalConfiguration = string(definition.physicalConfiguration); end
    if isfield(definition,"workload"), workload = string(definition.workload); end
    records(iKey) = struct("key",keys(iKey),"path","external raw archive fixture: "+string(extractAfter(pathname,workFolder+filesep)),"sha256",sha256File(pathname),"physicalConfiguration",physicalConfiguration,"workload",workload);
end
end

function key = fixtureKey(definition)
if isfield(definition,"physicalConfiguration") && isfield(definition,"workload")
    key = string(definition.physicalConfiguration)+"--"+string(definition.workload);
else
    key = "legacy";
end
end

function createMatchedFixture(pathname,options,definition,physicalState)
isHydrostatic = false;
isComposite = true;
if isfield(definition,"isHydrostatic"), isHydrostatic = logical(definition.isHydrostatic); end
if isfield(definition,"workload"), isComposite = string(definition.workload)=="composite-dense-output"; end
outputInterval = definition.outputInterval;
isIntegratorStudy = isfield(definition,"workload");
if isIntegratorStudy
    domain = definition.Lxyz;
else
    domain = [15000 15000 1300];
end
wvt = WVTransformConstantStratification(domain,options.Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=isHydrostatic,shouldAntialias=true);
if isIntegratorStudy
    advanceWaveVortexBenchmarkState(wvt,physicalState,0);
else
    state = initializeWaveVortexBenchmarkState(wvt,4001);
    advanceWaveVortexBenchmarkState(wvt,state,0);
end
model = WVModel(wvt); cleanup = onCleanup(@()closeModels(model));
if isComposite
    model.setFloatPositions([10e3 70e3],[9e3 65e3],[-250 -850],'u');
    model.addTracer(sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly),"dye");
end
if ~isIntegratorStudy
    model.setupIntegrator(integratorType="fixed",deltaT=options.deltaT);
    model.integrateToTime(options.deltaT,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
end
endpointInterval = conditional(isIntegratorStudy,definition.finalTime,outputInterval);
outputFile = model.createNetCDFFileForModelOutput(pathname,outputInterval=endpointInterval,shouldOverwriteExisting=true);
if isComposite && isIntegratorStudy
    group = outputFile.outputGroupWithName(model.defaultOutputGroupName());
    particles = model.fluxedObservingSystemWithName("float");
    tracer = model.fluxedObservingSystemWithName("dye");
    group.removeObservingSystem([particles tracer]);
    denseGroup = outputFile.addNewEvenlySpacedOutputGroup("dense",outputInterval=outputInterval,initialTime=definition.denseOutputStartTime,finalTime=definition.denseOutputEndTime);
    denseGroup.addObservingSystem(WVEulerianFields(model,fieldNames={'u'}));
    denseGroup.addObservingSystem(WVMooring(model,name="mooring",x=[0 50e3],y=[0 40e3],trackedFieldNames={'u'}));
    particleGroup = outputFile.addNewEvenlySpacedOutputGroup("particles",outputInterval=outputInterval,initialTime=definition.denseOutputStartTime,finalTime=definition.denseOutputEndTime);
    particleGroup.addObservingSystem(particles);
    tracerGroup = outputFile.addNewEvenlySpacedOutputGroup("tracers",outputInterval=outputInterval,initialTime=definition.denseOutputStartTime,finalTime=definition.denseOutputEndTime);
    tracerGroup.addObservingSystem(tracer);
elseif isComposite
    model.eulerianObservingSystem.addNetCDFOutputVariables('u');
    group = outputFile.outputGroupWithName(model.defaultOutputGroupName());
    particles = model.fluxedObservingSystemWithName("float");
    tracer = model.fluxedObservingSystemWithName("dye");
    group.removeObservingSystem([particles tracer]);
    particleGroup = outputFile.addNewEvenlySpacedOutputGroup("particles",outputInterval=outputInterval);
    particleGroup.addObservingSystem(particles);
    tracerGroup = outputFile.addNewEvenlySpacedOutputGroup("tracers",outputInterval=outputInterval);
    tracerGroup.addObservingSystem(tracer);
    group.addObservingSystem(WVMooring(model,name="mooring",x=[0 5000],y=[0 4000],trackedFieldNames={'u'}));
end
model.outputTimesForIntegrationPeriod(model.t,model.t);
model.writeTimeStepToNetCDFFile(model.t);
model.recordNetCDFFileHistory();
model.closeNetCDFFile();
clear cleanup
end

function executables = buildStandaloneWorkers(repositoryRoot,capabilities)
buildDirectory = fullfile(tempdir,"wave-vortex-model-issue-312-three-interface-build");
providerRoot = fileparts(fileparts(string(capabilities.libraries.base.path)));
configure = "cmake -S "+shellQuote(fullfile(repositoryRoot,"PortableRuntime"))+" -B "+shellQuote(buildDirectory)+" -DCMAKE_BUILD_TYPE=Release -DWV_RUNTIME_ENABLE_NATIVE_FFTW=ON -DWV_RUNTIME_BUILD_BENCHMARKS=ON -DWV_RUNTIME_FFTW_ROOT="+shellQuote(providerRoot);
[status,output] = system(configure); if status~=0, error("WaveVortexBenchmark:StandaloneBuild","%s",output); end
[status,output] = system("cmake --build "+shellQuote(buildDirectory)+" --parallel --target wave-vortex-run wv-standalone-nonlinear-flux-benchmark"); if status~=0, error("WaveVortexBenchmark:StandaloneBuild","%s",output); end
executables = struct("runner",fullfile(buildDirectory,"wave-vortex-run"),"kernel",fullfile(buildDirectory,"wv-standalone-nonlinear-flux-benchmark"));
end

function validateCapabilities(value)
if ~value.isAvailable || string(value.provider.id)~="native-neon-pthreads" || ~value.module.identityValidated || value.libraries.openmp.detected || value.contract.planCount~=17 || value.featureValidation.maximumRelativeError>1e-12
    error("WaveVortexBenchmark:ThreeInterfaceCapability","The benchmark requires the validated native FFTW provider.");
end
end

function writeArtifacts(results,options)
if ~options.shouldWriteArtifacts, return, end
writeText(fullfile(options.outputDirectory,"three-interface-benchmark.json"),jsonencode(results,PrettyPrint=true));
writeText(fullfile(options.outputDirectory,"summary.md"),summaryMarkdown(results));
end

function checkpoint(results,options)
if options.shouldWriteArtifacts, writeArtifacts(results,options); end
end

function archive = archiveDetailedArtifact(options)
if ~isfolder(options.archiveDirectory)
    mkdir(options.archiveDirectory);
end
rawPath = fullfile(options.outputDirectory,"three-interface-benchmark.json");
archiveName = options.runId+"-three-interface-benchmark.json.gz";
generated = gzip(rawPath,options.archiveDirectory);
generatedPath = string(generated{1});
archivePath = fullfile(options.archiveDirectory,archiveName);
if generatedPath ~= archivePath
    movefile(generatedPath,archivePath,"f");
end
information = dir(archivePath);
archive = struct("fileName",archiveName,"sha256",sha256File(archivePath),"compressedBytes",information.bytes,"location","external author archive; not distributed with source");
end

function markdown = summaryMarkdown(results)
lines = ["# Matched three-interface benchmark"; ""; "Status: `"+results.status+"`."; ""; "| Case | Interface | Process wall (s) | Integration (s) | Peak RSS (GiB) | Increment RSS (GiB) |"; "|---|---|---:|---:|---:|---:|"];
if ~isempty(results.comparison)
    for benchmarkCase = reshape(results.comparison,1,[])
        for item = reshape(benchmarkCase.interfaces,1,[])
            lines(end+1) = sprintf('| %s | %s | %.6f | %.6f | %.3f | %.3f |',benchmarkCase.id,item.id,item.processWallSeconds,item.integrationSeconds,item.totalPeakRSSBytes/2^30,item.incrementalPeakRSSBytes/2^30); %#ok<AGROW>
        end
    end
end
markdown = strjoin(lines,newline)+newline;
end

function value = readBinary(pathname)
fileId=fopen(pathname,"r"); cleanup=onCleanup(@()fclose(fileId)); values=fread(fileId,Inf,"double"); value=complex(values(1:2:end),values(2:2:end)); clear cleanup
end

function command = sampledCommand(workerCommand,samplePath,phasePath,stdoutPath,stderrPath,options,benchmarkFolder)
sampler = fullfile(benchmarkFolder,"runProcessWithRSS.sh");
if ~isfile(sampler)
    error("WaveVortexBenchmark:MissingRSSSampler","The process-tree RSS sampler is missing.");
end
command = shellQuote(sampler)+" "+shellQuote(samplePath)+" "+shellQuote(phasePath)+" "+numberText(options.samplingIntervalSeconds)+" "+shellQuote(stdoutPath)+" "+shellQuote(stderrPath)+" -- /bin/sh -c "+shellQuote(workerCommand);
end

function output = readCommandOutput(stdoutPath,stderrPath)
output = "";
if isfile(stdoutPath), output = string(fileread(stdoutPath)); end
if isfile(stderrPath), output = output+newline+string(fileread(stderrPath)); end
end

function memory = processMemory(samplePath,interval)
memory = struct("status","failed","provider","macos-ps-process-tree","boundary","integration-phase-total-live-process-tree-rss","samplingIntervalSeconds",interval,"totalPeakRSSBytes",NaN,"processLifetimePeakRSSBytes",NaN,"baselineProcessBytes",NaN,"peakIncrementBytes",NaN,"finalRSSBytes",NaN,"integrationSampleCount",0,"maximumProcessCount",0,"samples",[]);
if ~isfile(samplePath)
    return
end
lines = splitlines(strtrim(string(fileread(samplePath))));
lines(lines=="") = [];
samples = repmat(struct("sampleIndex",0,"elapsedSeconds",0,"phase","","rssBytes",0,"processCount",0),0,1);
for iLine = 1:numel(lines)
    fields = split(lines(iLine),sprintf('\t'));
    if numel(fields) < 4
        continue
    end
    index = str2double(fields(1));
    rssBytes = 1024*str2double(fields(3));
    processCount = str2double(fields(4));
    if ~isfinite(index) || ~isfinite(rssBytes) || rssBytes<=0 || ~isfinite(processCount) || processCount<1
        continue
    end
    samples(end+1,1) = struct("sampleIndex",index,"elapsedSeconds",index*interval,"phase",fields(2),"rssBytes",rssBytes,"processCount",processCount); %#ok<AGROW>
end
bytes = [samples.rssBytes];
phases = string({samples.phase});
baseline = bytes(phases=="steady-retained");
integrationMask = phases=="integrate" | startsWith(phases,"output-committed:");
integrationBytes = bytes(integrationMask);
memory.samples = samples;
if isempty(samples) || isempty(integrationBytes)
    return
end
memory.status = "complete";
memory.totalPeakRSSBytes = max(integrationBytes);
memory.processLifetimePeakRSSBytes = max(bytes);
memory.integrationSampleCount = nnz(integrationMask);
if ~isempty(baseline)
    memory.baselineProcessBytes = median(baseline);
    memory.peakIncrementBytes = max(0,memory.totalPeakRSSBytes-memory.baselineProcessBytes);
end
memory.finalRSSBytes = bytes(end);
memory.maximumProcessCount = max([samples.processCount]);
end

function value = mergeStruct(first,second)
value=first; names=fieldnames(second); for iName=1:numel(names), value.(names{iName})=second.(names{iName}); end
end
function value=safeRatio(numerator,denominator), if denominator==0, value=conditional(numerator==0,1,Inf); else, value=numerator/denominator; end, end
function value=numberText(value), value=string(sprintf('%.17g',value)); end
function value=shellQuote(value), value="'"+replace(string(value),"'","'""'""'")+"'"; end
function closeModels(varargin), for i=1:numel(varargin), value=varargin{i}; if ~isempty(value)&&isvalid(value), try, value.closeNetCDFFile(); catch, end; delete(value); end, end, end
function removeFolder(pathname), if isfolder(pathname), rmdir(pathname,"s"); end, end
function verifyTemporaryCapacity(pathname,options)
realGridBytes = prod(double(options.Nxyz))*8;
requiredBytes = max(256*2^20,ceil(640*realGridBytes));
[status,output] = system("/bin/df -Pk "+shellQuote(pathname));
lines = splitlines(strtrim(string(output)));
if status~=0 || numel(lines)<2
    error("WaveVortexBenchmark:ThreeInterfaceDiskCapacity","Unable to determine free space for the three-interface benchmark temporary directory.");
end
fields = split(strtrim(lines(end)));
fields(fields=="") = [];
availableBytes = 1024*str2double(fields(4));
if ~isfinite(availableBytes) || availableBytes<requiredBytes
    error("WaveVortexBenchmark:ThreeInterfaceDiskCapacity","The canonical benchmark requires approximately %.1f GiB of temporary free space, but %.1f GiB is available.",requiredBytes/2^30,availableBytes/2^30);
end
end
function restoreState(directory,originalPath,originalRng), cd(directory); path(originalPath); rng(originalRng); end
function addRepositoryPaths(root,benchmark), metadata=jsondecode(fileread(fullfile(root,"resources","mpackage.json"))); for item=reshape(metadata.folders,1,[]), folder=fullfile(root,item.path); if isfolder(folder), addpath(folder); end, end, addpath(root); addpath(benchmark); end
function writeText(pathname,contents), parent=fileparts(pathname); if ~isfolder(parent), mkdir(parent); end, temporaryPath=pathname+".tmp"; if isfile(temporaryPath), delete(temporaryPath); end, [fileId,message]=fopen(temporaryPath,"w"); if fileId<0, error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to open %s: %s",temporaryPath,message); end, cleanup=onCleanup(@()fclose(fileId)); fprintf(fileId,"%s",contents); clear cleanup, [status,message]=movefile(temporaryPath,pathname,"f"); if ~status, error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to publish %s: %s",pathname,message); end, end
function value=sha256File(pathname), [status,output]=system(sprintf('/usr/bin/shasum -a 256 %s',shellQuote(pathname))); if status~=0, error("WaveVortexBenchmark:HashFailed","%s",output); end, value=string(extractBefore(strtrim(output),65)); end
function [commit,tree,isDirty]=gitIdentity(root), commit=gitValue(root,"rev-parse HEAD"); tree=gitValue(root,"rev-parse HEAD^{tree}"); [status,~]=system("git -C "+shellQuote(root)+" diff --quiet && git -C "+shellQuote(root)+" diff --cached --quiet"); isDirty=status~=0; end
function value=gitValue(root,args), [status,output]=system("git -C "+shellQuote(root)+" "+args); if status~=0,error("WaveVortexBenchmark:GitIdentity","%s",output);end,value=strtrim(string(output));end
function value=environmentRecord, [~,processor]=system("/usr/sbin/sysctl -n machdep.cpu.brand_string"); [~,memory]=system("/usr/sbin/sysctl -n hw.memsize"); value=struct("processor",strtrim(string(processor)),"physicalMemoryBytes",str2double(memory),"os",string(system_dependent("getos")),"matlabVersion",string(version),"release",string(version("-release")),"architecture",string(computer("arch")),"requestedThreads",min(18,maxNumCompThreads)); end
function value=utcTimestamp, value=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")); end
function value=conditional(condition,a,b), if condition,value=a;else,value=b;end,end
function value=emptyFailure, value=struct("stage","","identifier","","message","","report",""); end
function value=emptyOutputGraph, value=struct("kind","","passed",true,"maximumRelativeError",0,"maximumAbsoluteError",0,"variableCount",0,"recordCount",0,"categories",repmat(struct("name","","variableCount",0,"maximumAbsoluteError",0,"maximumRelativeError",0,"passed",true),0,1),"differences",strings(0,1)); end
function value=emptyRun, value=struct("schemaVersion","three-interface-worker-v1","status","failed","interface","","case",struct(),"repeatIndex",0,"sourceCommit","","processWallSeconds",NaN,"interfaceTotalSeconds",NaN,"integrationSeconds",NaN,"memory",struct(),"provider",struct(),"integrator",struct(),"finalState",struct(),"output",struct(),"diagnostics",struct(),"failure",struct("identifier","","message","","report","")); end
function value=outputRecordCounts(pathname), value=struct("waveVortex",netCDFRecordCount(pathname,"/wave-vortex/t"),"dense",netCDFRecordCount(pathname,"/dense/t"),"particles",netCDFRecordCount(pathname,"/particles/t"),"tracers",netCDFRecordCount(pathname,"/tracers/t")); end
function value=netCDFRecordCount(pathname,variable), try, value=numel(ncread(pathname,variable)); catch, value=0; end, end
