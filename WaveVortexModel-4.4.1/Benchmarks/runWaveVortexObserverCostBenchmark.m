function results = runWaveVortexObserverCostBenchmark(options)
% Separate integrated-observer state costs from dense-output delivery costs.
arguments
    options.Nxyz (1,3) double {mustBeInteger,mustBePositive} = [64 64 33]
    options.Lxyz (1,3) double {mustBePositive} = [150e3 150e3 1300]
    options.deltaT (1,1) double {mustBePositive} = 128
    options.integrationStepCount (1,1) double {mustBeInteger,mustBePositive} = 8
    options.denseOutputPointsPerStep (1,1) double {mustBeInteger,mustBePositive} = 3
    options.seed (1,1) double {mustBeInteger,mustBeNonnegative} = 4001
    options.caseIds (1,:) string = strings(1,0)
    options.shouldUseFreshProcess (1,1) logical = true
    options.samplingIntervalSeconds (1,1) double {mustBePositive} = 0.005
    options.plateauSeconds (1,1) double {mustBeNonnegative} = 0.05
    options.phasePath (1,1) string = ""
end

benchmarkFolder = string(fileparts(mfilename("fullpath")));
repositoryRoot = string(fileparts(benchmarkFolder));
originalPath = path;
originalRng = rng;
cleanup = onCleanup(@()restoreState(originalPath,originalRng));
addRepositoryPaths(repositoryRoot,benchmarkFolder);
definitions = caseDefinitions(options.caseIds,options.deltaT,options.integrationStepCount,options.denseOutputPointsPerStep);
configuration = struct("Nxyz",options.Nxyz,"Lxyz",options.Lxyz,"deltaT",options.deltaT,"integrationStepCount",options.integrationStepCount,"denseOutputPointsPerStep",options.denseOutputPointsPerStep,"seed",options.seed);
cases = cell(numel(definitions),1);
if options.shouldUseFreshProcess
    workFolder = string(tempname);
    mkdir(workFolder);
    workCleanup = onCleanup(@()rmdir(workFolder,"s"));
    for iCase = 1:numel(definitions)
        fprintf("Observer-cost benchmark: %s.\n",definitions(iCase).id);
        cases{iCase} = runFreshCase(definitions(iCase),configuration,benchmarkFolder,workFolder,options);
    end
    clear workCleanup
else
    plateauSeconds = options.plateauSeconds*(options.phasePath~="");
    for iCase = 1:numel(definitions)
        cases{iCase} = runCase(definitions(iCase),configuration,options.phasePath,plateauSeconds);
        cases{iCase}.memory = struct("provider","not-sampled","boundary","in-process profiling mode");
    end
end
cases = vertcat(cases{:});
results = struct("studyId","observer-dense-cost-decomposition-v1", ...
    "environment",struct("matlabRelease",string(version("-release")),"architecture",string(computer("arch"))), ...
    "configuration",configuration,"cases",cases);
clear cleanup
end

function value = runFreshCase(definition,configuration,benchmarkFolder,workFolder,options)
caseFolder = fullfile(workFolder,definition.id);
mkdir(caseFolder);
configPath = fullfile(caseFolder,"config.json");
outputPath = fullfile(caseFolder,"result.json");
phasePath = fullfile(caseFolder,"phase.txt");
samplePath = fullfile(caseFolder,"rss.tsv");
stdoutPath = fullfile(caseFolder,"stdout.txt");
stderrPath = fullfile(caseFolder,"stderr.txt");
writeText(configPath,jsonencode(struct("caseId",definition.id,"configuration",configuration,"phasePath",phasePath,"plateauSeconds",options.plateauSeconds)));
statement = "addpath('"+replace(benchmarkFolder,"'","''")+"'); waveVortexObserverCostBenchmarkWorker('"+replace(configPath,"'","''")+"','"+replace(outputPath,"'","''")+"')";
worker = sprintf('"%s" -batch "%s"',fullfile(matlabroot,"bin","matlab"),replace(statement,'"','\"'));
sampler = fullfile(benchmarkFolder,"runProcessWithRSS.sh");
command = shellQuote(sampler)+" "+shellQuote(samplePath)+" "+shellQuote(phasePath)+" "+string(sprintf('%.17g',options.samplingIntervalSeconds))+" "+shellQuote(stdoutPath)+" "+shellQuote(stderrPath)+" -- /bin/sh -c "+shellQuote(worker);
[exitCode,~] = system(command);
if exitCode~=0 || ~isfile(outputPath)
    error("WaveVortexBenchmark:ObserverCostWorkerFailed","Fresh MATLAB failed for %s:\n%s",definition.id,commandOutput(stdoutPath,stderrPath));
end
payload = jsondecode(fileread(outputPath));
if string(payload.status)~="complete"
    error("WaveVortexBenchmark:ObserverCostWorkerFailed","Fresh MATLAB failed for %s:\n%s",definition.id,payload.failure);
end
value = payload.case;
value.memory = processMemory(samplePath,options.samplingIntervalSeconds);
end

function value = runCase(definition,configuration,phasePath,plateauSeconds)
workFolder = string(tempname);
mkdir(workFolder);
outputPath = fullfile(workFolder,"output.nc");
fileCleanup = onCleanup(@()rmdir(workFolder,"s"));
writePhase(phasePath,"construction");
wvt = WVTransformConstantStratification(configuration.Lxyz,configuration.Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=true);
wvtCleanup = onCleanup(@()closeTransform(wvt));
initializeThreeInterfaceIntegratorState(wvt,configuration.seed);
model = WVModel(wvt);
modelCleanup = onCleanup(@()closeModel(model));
if definition.integratesObserverState
    model.setFloatPositions([10e3 70e3],[9e3 65e3],[-250 -850],'u');
    model.addTracer(sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly),"dye");
end
configureOutput(model,outputPath,definition);
model.setupIntegrator(integratorType="fixed",deltaT=configuration.deltaT);
stateElementCount = integratedStateElementCount(model);
writePhase(phasePath,"steady-retained");
pause(plateauSeconds);
rhsBefore = double(model.nFluxComputations);
writePhase(phasePath,"integrate");
waitForPhaseSample(phasePath);
timer = tic;
model.integrateToTime(definition.finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
elapsedSeconds = toc(timer);
writePhase(phasePath,"integration-complete");
acceptedStepCount = double(model.rk4Integrator.totalIterations);
rhsEvaluationCount = double(model.nFluxComputations)-rhsBefore;
finalState = struct("t",wvt.t,"ApNorm",norm(wvt.Ap(:)),"AmNorm",norm(wvt.Am(:)),"A0Norm",norm(wvt.A0(:)));
model.closeNetCDFFile();
recordCounts = [recordCount(outputPath,"/wave-vortex/t") recordCount(outputPath,"/dense/t") recordCount(outputPath,"/particles/t") recordCount(outputPath,"/tracers/t")];
if acceptedStepCount~=definition.integrationStepCount || rhsEvaluationCount~=1+4*definition.integrationStepCount
    error("WaveVortexBenchmark:ObserverCostWorkMismatch","Case %s did not execute its required fixed-RK4 work.",definition.id);
end
if ~isequal(recordCounts,definition.expectedRecordCounts)
    error("WaveVortexBenchmark:ObserverCostOutputMismatch","Case %s did not deliver its configured output graph.",definition.id);
end
value = rmfield(definition,"expectedRecordCounts");
value.elapsedSeconds = elapsedSeconds;
value.work = struct("acceptedStepCount",acceptedStepCount,"rhsEvaluationCount",rhsEvaluationCount,"scheduledInteriorOutputCount",definition.denseOutputPointsPerStep*definition.usesDenseOutput,"outputRecordCounts",recordCounts,"outputRecordNames",["waveVortex" "dense" "particles" "tracers"],"integratedStateElementCount",stateElementCount,"integratedStateNominalBytes",8*stateElementCount);
value.finalState = finalState;
writePhase(phasePath,"outputs-held");
pause(plateauSeconds);
writePhase(phasePath,"complete");
clear modelCleanup wvtCleanup fileCleanup
end

function configureOutput(model,outputPath,definition)
outputFile = model.createNetCDFFileForModelOutput(outputPath,outputInterval=definition.finalTime,shouldOverwriteExisting=true);
if definition.denseOutputTarget=="none", return, end
dense = outputFile.addNewEvenlySpacedOutputGroup("dense",outputInterval=definition.denseOutputInterval,initialTime=0,finalTime=definition.deltaT);
if definition.denseOutputTarget=="coefficients"
    dense.addObservingSystem(model.wvCoefficientFluxedObservingSystem);
    return
end
particles = model.fluxedObservingSystemWithName("float");
tracer = model.fluxedObservingSystemWithName("dye");
outputFile.outputGroupWithName(model.defaultOutputGroupName()).removeObservingSystem([particles tracer]);
dense.addObservingSystem(WVEulerianFields(model,fieldNames={'u'}));
dense.addObservingSystem(WVMooring(model,name="mooring",x=[0 50e3],y=[0 40e3],trackedFieldNames={'u'}));
particleGroup = outputFile.addNewEvenlySpacedOutputGroup("particles",outputInterval=definition.denseOutputInterval,initialTime=0,finalTime=definition.deltaT);
particleGroup.addObservingSystem(particles);
tracerGroup = outputFile.addNewEvenlySpacedOutputGroup("tracers",outputInterval=definition.denseOutputInterval,initialTime=0,finalTime=definition.deltaT);
tracerGroup.addObservingSystem(tracer);
end

function definitions = caseDefinitions(caseIds,deltaT,stepCount,densePointCount)
finalTime = stepCount*deltaT;
denseInterval = deltaT/(densePointCount+1);
denseRecords = densePointCount+2;
definitions = [ ...
    definition("coefficient-endpoint","Coefficients | endpoint output","coefficients only","endpoint output",false,"none",[2 0 0 0]); ...
    definition("coefficient-dense-output","Coefficients | dense coefficient output","coefficients only","first-step dense coefficient output",false,"coefficients",[2 denseRecords 0 0]); ...
    definition("integrated-observer-endpoint","Integrated tracer/particles | endpoint output","coefficients plus tracer and particles","endpoint output",true,"none",[2 0 0 0]); ...
    definition("composite-dense-output","Integrated tracer/particles | composite dense output","coefficients plus tracer and particles","first-step composite dense-output graph",true,"composite",[2 denseRecords denseRecords denseRecords])];
for iCase = 1:numel(definitions)
    definitions(iCase).deltaT = deltaT;
    definitions(iCase).finalTime = finalTime;
    definitions(iCase).integrationStepCount = stepCount;
    definitions(iCase).denseOutputPointsPerStep = densePointCount;
    definitions(iCase).denseOutputInterval = denseInterval;
end
knownIds = string({definitions.id});
if isempty(caseIds), return, end
if numel(unique(caseIds,"stable"))~=numel(caseIds)
    error("WaveVortexBenchmark:DuplicateObserverCostCase","Observer-cost benchmark case IDs must be unique.");
end
indices = arrayfun(@(id)find(knownIds==id,1),caseIds,UniformOutput=false);
if any(cellfun(@isempty,indices))
    error("WaveVortexBenchmark:UnknownObserverCostCase","Unknown observer-cost benchmark case: %s",strjoin(caseIds(cellfun(@isempty,indices)),", "));
end
definitions = definitions([indices{:}]);
end

function value = definition(id,label,integrationLabel,outputLabel,integratesObserverState,denseOutputTarget,expectedRecordCounts)
value = struct("id",id,"label",label,"integrationLabel",integrationLabel,"outputLabel",outputLabel,"integratesObserverState",integratesObserverState,"usesDenseOutput",denseOutputTarget~="none","denseOutputTarget",denseOutputTarget,"expectedRecordCounts",expectedRecordCounts);
end

function value = processMemory(samplePath,interval)
if ~isfile(samplePath), error("WaveVortexBenchmark:ObserverCostMemoryUnavailable","RSS samples are missing."); end
lines = splitlines(strtrim(string(fileread(samplePath))));
steady = []; integration = [];
for line = reshape(lines,1,[])
    fields = split(line,sprintf('\t'));
    if numel(fields)<4, continue, end
    bytes = 1024*str2double(fields(3));
    if fields(2)=="steady-retained", steady(end+1)=bytes; end %#ok<AGROW>
    if fields(2)=="integrate", integration(end+1)=bytes; end %#ok<AGROW>
end
if isempty(steady) || isempty(integration)
    error("WaveVortexBenchmark:ObserverCostMemoryUnavailable","RSS sampling did not capture steady-retained and integration phases.");
end
value = struct("provider","macos-ps-process-tree","boundary","total live-process-tree RSS during integration","samplingIntervalSeconds",interval,"peakRSSBytes",max(integration),"baselineRSSBytes",median(steady),"peakIncrementBytes",max(0,max(integration)-median(steady)),"integrationSampleCount",numel(integration));
end

function value = integratedStateElementCount(model)
value = 0;
for system = reshape(model.fluxedObservingSystems,1,[]), value=value+sum(system.lengthOfFluxComponents()); end
end

function value = recordCount(pathname,variable)
try
    value = numel(ncread(pathname,variable));
catch
    value = 0;
end
end

function writePhase(pathname,phase)
if pathname=="", return, end
writeText(pathname+".tmp",phase);
movefile(pathname+".tmp",pathname,"f");
end

function waitForPhaseSample(phasePath)
if phasePath=="", return, end
acknowledgementPath = string(getenv("WV_RSS_PHASE_ACK"));
timer = tic;
while acknowledgementPath~="" && toc(timer)<5
    if isfile(acknowledgementPath) && strtrim(string(fileread(acknowledgementPath)))=="integrate", return, end
    pause(0.001)
end
if acknowledgementPath~="", error("WaveVortexBenchmark:RSSPhaseHandshake","The RSS sampler did not acknowledge integration."); end
end

function addRepositoryPaths(repositoryRoot,benchmarkFolder)
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
for item=reshape(metadata.folders,1,[]), folder=fullfile(repositoryRoot,item.path); if isfolder(folder), addpath(folder); end, end
addpath(repositoryRoot);
addpath(benchmarkFolder);
end

function closeModel(model)
if isempty(model) || ~isvalid(model), return, end
try
    model.closeNetCDFFile();
catch
end
delete(model);
end

function closeTransform(wvt)
if ~isempty(wvt) && isvalid(wvt), delete(wvt); end
end

function value = commandOutput(stdoutPath,stderrPath)
value = "";
if isfile(stdoutPath), value=string(fileread(stdoutPath)); end
if isfile(stderrPath), value=value+newline+string(fileread(stderrPath)); end
end

function writeText(pathname,contents)
fileId = fopen(pathname,"w");
if fileId<0, error("WaveVortexBenchmark:ArtifactWriteFailed","Unable to open %s.",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",contents);
clear cleanup
end

function value = shellQuote(value), value="'"+replace(string(value),"'","'""'""'")+"'"; end
function restoreState(originalPath,originalRng), path(originalPath); rng(originalRng); end
