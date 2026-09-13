function report = hydroCompositeProfile(fixturePath,outputPath,options)
% Run the hydrostatic composite audit matched-model diagnostic from an immutable fixture.
arguments
    fixturePath (1,1) string {mustBeFile}
    outputPath (1,1) string
    options.backend (1,1) string {mustBeMember(options.backend,["compiled" "matlab"])} = "compiled"
    options.dependencyPaths (1,:) string = strings(1,0)
    options.shouldProfile (1,1) logical = false
    options.profilingToolsPath (1,1) string = ""
    options.reportPath (1,1) string = ""
    options.jsonReportPath (1,1) string = ""
    options.phasePath (1,1) string = ""
end

repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addRepositoryPaths(repositoryRoot,options.dependencyPaths);
if options.shouldProfile
    addProfilingTools(options.profilingToolsPath,repositoryRoot);
end

fixture = fixtureIdentity(fixturePath);
if isfile(outputPath)
    error("WaveVortexBenchmark:HydroAuditOutputExists","Refusing to overwrite the diagnostic output: %s",outputPath)
end
outputFolder = string(fileparts(outputPath));
if outputFolder ~= "" && ~isfolder(outputFolder), mkdir(outputFolder); end
copyfile(fixturePath,outputPath);

model = benchmarkModelFromFile(outputPath,options.backend);
wvt = model.wvt;
modelCleanup = onCleanup(@()closeModelAndTransform(model,wvt));
validateBackend(wvt,options.backend);
validateScientificFixture(wvt,outputPath,fixture.workload);
if double(model.t) ~= 0
    error("WaveVortexBenchmark:HydroAuditFixtureTime","The copied fixture must restart at t=0; found %.17g.",double(model.t))
end
model.setupIntegrator(integratorType="adaptive",integrator=@ode78,relTolerance=1e-3,absTolerance=1e-6,shouldShowIntegrationStats=0);
model.odeOptions = odeset(model.odeOptions,'InitialStep',269.39985873487035,'Stats','on');
tolerances = odeget(model.odeOptions,'AbsTol');
integrator = struct( ...
    "requested","adaptive-rk78", ...
    "actual","adaptive-rk78", ...
    "relativeTolerance",double(odeget(model.odeOptions,'RelTol')), ...
    "absoluteToleranceHash",threeInterfaceToleranceHash(tolerances), ...
    "absoluteToleranceHashClearedMantissaBits",20, ...
    "absoluteToleranceComponentHashes",arrayfun(@(index)threeInterfaceToleranceHash(tolerances(model.arrayStartIndex(index):model.arrayEndIndex(index))),1:numel(model.arrayStartIndex)), ...
    "requestedInitialStep",269.39985873487035, ...
    "effectiveInitialStep",double(odeget(model.odeOptions,'InitialStep')), ...
    "requestedMaximumStep",[], ...
    "maximumStepPolicy","matlab-default", ...
    "effectiveMaximumStep",716.8);

metricsBefore = runtimeMetrics(wvt);
rhsEvaluationsBefore = double(model.nFluxComputations);
profileAnalysis = [];
writePhase(options.phasePath,"integrate");
if options.shouldProfile
    profileLabel = "hydrostatic composite audit matched adaptive ode78";
    statisticsText = evalc("profileAnalysis = profileCodeHotspots(@()runIntegration(model),projectRoots=repositoryRoot,label=profileLabel,shouldPrintReport=false);");
    integrationSeconds = double(profileAnalysis.elapsedTime);
else
    operationTimer = tic;
    statisticsText = evalc("runIntegration(model)");
    integrationSeconds = toc(operationTimer);
end
writePhase(options.phasePath,"integration-complete");
integrator.acceptedStepCount = statisticCount(statisticsText,"successful steps");
integrator.rejectedStepCount = statisticCount(statisticsText,"failed attempts");
integrator.reportedRightHandSideEvaluationCount = statisticCount(statisticsText,"function evaluations");
rhsEvaluationCount = double(model.nFluxComputations)-rhsEvaluationsBefore;
if integrator.reportedRightHandSideEvaluationCount ~= rhsEvaluationCount
    error("WaveVortexBenchmark:HydroAuditWorkCount","MATLAB ode78 reported %d function evaluations, but WVModel counted %d.",integrator.reportedRightHandSideEvaluationCount,rhsEvaluationCount)
end
integrator.rhsEvaluationCount = rhsEvaluationCount;
metricsAfter = runtimeMetrics(wvt);
runtimeMetricDeltas = numericScalarDeltas(metricsBefore,metricsAfter);
counterMetrics = selectCounterMetrics(metricsAfter);
counterMetricDeltas = selectCounterMetrics(runtimeMetricDeltas);

model.closeNetCDFFile();
report = struct( ...
    "schemaVersion","wvm-hydroaudit-profile-v1", ...
    "fixture",fixture, ...
    "outputPath",outputPath, ...
    "backend",options.backend, ...
    "controls",struct("integrator","ode78","initialTime",0,"finalTime",7168,"initialStep",269.39985873487035,"relativeTolerance",1e-3,"absoluteTolerance",1e-6,"maximumStep",716.8,"maximumStepPolicy","MATLAB default: 0.1 integration span"), ...
    "integrationSeconds",integrationSeconds, ...
    "rhsEvaluationCount",rhsEvaluationCount, ...
    "integrator",integrator, ...
    "finalState",struct("t",double(wvt.t),"ApNorm",norm(wvt.Ap(:)),"AmNorm",norm(wvt.Am(:)),"A0Norm",norm(wvt.A0(:))), ...
    "outputRecordCounts",outputRecordCounts(outputPath), ...
    "backendMetadata",wvt.computationalBackendMetadata, ...
    "runtimeMetricsBefore",metricsBefore, ...
    "runtimeMetricsAfter",metricsAfter, ...
    "runtimeMetricDeltas",runtimeMetricDeltas, ...
    "producerPrimitiveCopyCounters",counterMetrics, ...
    "producerPrimitiveCopyCounterDeltas",counterMetricDeltas, ...
    "profileEnabled",options.shouldProfile);
if options.reportPath == "", options.reportPath = outputPath+".hydroaudit.mat"; end
save(options.reportPath,"report","profileAnalysis","-v7.3");
if options.jsonReportPath == "", options.jsonReportPath = outputPath+".hydroaudit.json"; end
writeText(options.jsonReportPath,jsonencode(report,PrettyPrint=true));
closeModelAndTransform(model,wvt);
model = [];
wvt = [];
clear modelCleanup pathCleanup
end

function runIntegration(model)
model.integrateToTime(7168,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
end

function value = statisticCount(text,label)
token = regexp(text,'(?m)^\s*(\d+)\s+'+label,"tokens","once");
if isempty(token), error("WaveVortexBenchmark:HydroAuditStatistics","MATLAB ode78 did not report %s.",label); end
value = str2double(token{1});
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
    model.addOutputFile(WVModelOutputFile.modelOutputFileFromFile(ncfile,model));
catch exception
    if ~isempty(ncfile.id), ncfile.close(); end
    rethrow(exception)
end
end

function validateBackend(wvt,backend)
metadata = wvt.computationalBackendMetadata;
if string(metadata.activeBackend) ~= backend
    error("WaveVortexBenchmark:HydroAuditBackendFallback","Requested %s but %s executed.",backend,string(metadata.activeBackend))
end
if backend == "compiled" && (string(metadata.provider.id) ~= "native-neon-pthreads" || ~metadata.module.identityValidated || metadata.libraries.openmp.detected)
    error("WaveVortexBenchmark:HydroAuditProvider","The compiled run did not use the validated native provider.")
end
end

function value = runtimeMetrics(wvt)
value = struct();
if isprop(wvt,"computationalBackendMetadata")
    metadata = wvt.computationalBackendMetadata;
    if isfield(metadata,"runtimeMetrics"), value = metadata.runtimeMetrics; end
end
end

function value = numericScalarDeltas(before,after)
value = struct();
names = intersect(string(fieldnames(before)),string(fieldnames(after)),"stable");
for name = reshape(names,1,[])
    first = before.(name);
    second = after.(name);
    if isnumeric(first) && isscalar(first) && isnumeric(second) && isscalar(second)
        value.(name) = double(second)-double(first);
    end
end
end

function value = selectCounterMetrics(metrics)
value = struct();
names = string(fieldnames(metrics));
selected = contains(lower(names),["producer" "primitive" "copy"]);
for name = reshape(names(any(selected,2)),1,[]), value.(name) = metrics.(name); end
end

function value = outputRecordCounts(pathname)
value = struct("waveVortex",recordCount(pathname,"/wave-vortex/t"),"dense",recordCount(pathname,"/dense/t"),"particles",recordCount(pathname,"/particles/t"),"tracers",recordCount(pathname,"/tracers/t"));
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

function value = fixtureIdentity(pathname)
hash = sha256File(pathname);
endpointHash = "23ecd022cb7094db8bf2cfdb7b1b224e63b5b209a2241692e91a39498cad4f16";
compositeHash = "ec97d7a1fa48aca4fcbe456f18b9991d84c886c8b95150d3551fa964364f8a4b";
if hash == endpointHash
    workload = "coefficient-endpoint";
elseif hash == compositeHash
    workload = "composite-dense-output";
else
    filename = string(pathname);
    if contains(filename,"coefficient-endpoint")
        workload = "coefficient-endpoint";
    elseif contains(filename,"composite-dense-output")
        workload = "composite-dense-output";
    else
        workload = "unknown";
    end
end
information = dir(pathname);
value = struct("path",pathname,"sha256",hash,"bytes",information.bytes,"workload",workload,"matchesQualifiedArchiveHash",ismember(hash,[endpointHash compositeHash]));
end

function value = sha256File(pathname)
[status,output] = system(sprintf('/usr/bin/shasum -a 256 %s',shellQuote(pathname)));
if status ~= 0, error("WaveVortexBenchmark:HydroAuditHashFailed","%s",output); end
value = string(extractBefore(strtrim(output),65));
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function addRepositoryPaths(repositoryRoot,dependencyPaths)
for dependencyPath = reshape(dependencyPaths,1,[])
    if ~isfolder(dependencyPath), error("WaveVortexBenchmark:HydroAuditDependency","Dependency folder does not exist: %s",dependencyPath); end
    addpath(dependencyPath);
end
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
for item = reshape(metadata.folders,1,[])
    folder = fullfile(repositoryRoot,item.path);
    if isfolder(folder), addpath(folder); end
end
addpath(repositoryRoot);
addpath(fullfile(repositoryRoot,"Benchmarks"));
end

function addProfilingTools(profilingToolsPath,repositoryRoot)
if profilingToolsPath == ""
    profilingToolsPath = fullfile(fileparts(repositoryRoot),"OceanKit","tools","profiling");
end
if isfolder(profilingToolsPath), addpath(profilingToolsPath); end
if exist("profileCodeHotspots","file") ~= 2
    error("WaveVortexBenchmark:HydroAuditProfilingTools","Add OceanKit/tools/profiling with profilingToolsPath when shouldProfile is true.")
end
end

function validateScientificFixture(wvt,pathname,workload)
valid = isa(wvt,"WVTransformHydrostatic") && isequal([wvt.Lx wvt.Ly wvt.Lz],[150e3 150e3 1300]);
valid = valid && isequal([wvt.Nx wvt.Ny wvt.Nz],[256 256 129]) && max(abs(wvt.N2-2e-5*exp(2*wvt.z/1300))) < 1e-17 && wvt.latitude == 45;
valid = valid && wvt.isHydrostatic && wvt.shouldAntialias;
if ~valid
    error("WaveVortexBenchmark:HydroAuditScientificFixture","The copied fixture does not have the frozen exponential hydrostatic 256-by-256-by-129 scientific configuration.")
end
counts = outputRecordCounts(pathname);
recordCounts = [counts.waveVortex counts.dense counts.particles counts.tracers];
if workload == "coefficient-endpoint" && ~isequal(recordCounts,[1 0 0 0])
    error("WaveVortexBenchmark:HydroAuditFixtureOutputs","The endpoint fixture does not contain exactly its one initial coefficient record.")
elseif workload == "composite-dense-output" && ~isequal(recordCounts,[1 1 1 1])
    error("WaveVortexBenchmark:HydroAuditFixtureOutputs","The composite fixture does not contain exactly one initial record in each output group.")
elseif workload == "unknown"
    error("WaveVortexBenchmark:HydroAuditFixtureWorkload","Name a regenerated fixture with coefficient-endpoint or composite-dense-output so its output contract can be checked.")
end
end

function writePhase(pathname,phase)
if pathname == "", return, end
fileId = fopen(pathname,"w");
if fileId < 0, error("WaveVortexBenchmark:HydroAuditPhase","Unable to write phase file: %s",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",phase);
clear cleanup
end

function writeText(pathname,contents)
parent = fileparts(pathname);
if parent ~= "" && ~isfolder(parent), mkdir(parent); end
fileId = fopen(pathname,"w");
if fileId < 0, error("WaveVortexBenchmark:HydroAuditReport","Unable to write report: %s",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",contents);
clear cleanup
end

function closeReader(reader)
if ~isempty(reader) && isvalid(reader) && ~isempty(reader.id), reader.close(); end
end

function closeModel(model)
if isempty(model) || ~isvalid(model), return, end
try
    model.closeNetCDFFile();
catch
end
delete(model);
end

function closeModelAndTransform(model,wvt)
closeModel(model);
if ~isempty(wvt) && isvalid(wvt), delete(wvt); end
end
