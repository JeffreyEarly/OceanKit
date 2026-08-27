function dataset = publishedThreeInterfaceBenchmarkFromArtifact(rawArtifactPath,options)
% Normalize a complete matched three-interface artifact for website publication.
arguments
    rawArtifactPath (1,1) string {mustBeFile}
    options.platformId (1,1) string = "m5-max"
    options.platformName (1,1) string = "Apple M5 Max"
    options.archiveFileName (1,1) string = ""
    options.archiveSHA256 (1,1) string = ""
    options.archiveCompressedBytes (1,1) double {mustBeNonnegative} = 0
    options.implementationVersion (1,1) string = "unreleased-preview"
end
raw = jsondecode(fileread(rawArtifactPath));
if string(raw.schemaVersion)=="three-interface-benchmark-v2"
    if string(raw.status)~="complete" || logical(raw.source.isDirty)
        error("WaveVortexBenchmark:InvalidThreeInterfaceArtifact","Publication requires a complete clean three-interface-benchmark-v2 artifact.");
    end
    validateThreeInterfaceBenchmarkContract(raw);
    dataset = normalizeIntegratorStudy(raw,rawArtifactPath,options);
    return
end
if string(raw.schemaVersion) ~= "three-interface-benchmark-v1" || string(raw.status) ~= "complete" || logical(raw.source.isDirty)
    error("WaveVortexBenchmark:InvalidThreeInterfaceArtifact","Publication requires a complete clean three-interface-benchmark-v1 artifact.");
end
dataset = normalizeLegacy(raw,options);
end

function dataset = normalizeIntegratorStudy(raw,rawArtifactPath,options)
[collectedAt,timestamp] = collectionTime(raw.runId);
datasetId = "three-interface--"+options.platformId+"--"+timestamp;
platform = struct("id",options.platformId,"displayName",options.platformName,"processor",string(raw.environment.processor),"physicalMemoryBytes",double(raw.environment.physicalMemoryBytes),"os",string(raw.environment.os),"architecture",string(raw.environment.architecture),"matlabVersion",string(raw.environment.matlabVersion),"threadCount",double(raw.configuration.threadCount));
provider = struct("id",string(raw.provider.provider.id),"version",string(raw.provider.provider.version),"threadBackend",string(raw.provider.provider.threadBackend),"scope","compiled-interfaces-only","moduleSHA256",string(raw.provider.module.sha256),"identityValidated",logical(raw.provider.module.identityValidated),"openMPDetected",logical(raw.provider.libraries.openmp.detected));
cases = cell(1,numel(raw.comparison));
for iCase = 1:numel(raw.comparison)
    comparison = raw.comparison(iCase);
    definition = raw.cases(iCase);
    caseRuns = raw.runs(string(arrayfun(@(run)run.case.id,raw.runs,"UniformOutput",false))==string(definition.id));
    interfaces = cell(1,numel(comparison.interfaces));
    for iInterface = 1:numel(comparison.interfaces)
        item = comparison.interfaces(iInterface);
        selected = caseRuns(string({caseRuns.interface})==string(item.id));
        interfaces{iInterface} = struct("id",string(item.id),"providerId",string(selected(1).provider.id),"integrationSeconds",double(item.integrationSeconds),"totalPeakRSSBytes",double(item.totalPeakRSSBytes),"integrationRatio",double(item.integrationRatio),"totalRSSRatio",double(item.totalRSSRatio),"integrationSamplesSeconds",double([selected.integrationSeconds]),"totalPeakRSSSamplesBytes",double(arrayfun(@(run)run.memory.totalPeakRSSBytes,selected)),"diagnostics",interfaceDiagnostics(selected));
    end
    contract = struct("Nxyz",double(definition.Nxyz(:)'),"Lxyz",double(raw.configuration.Lxyz(:)'),"physicalConfiguration",string(definition.physicalConfiguration),"isHydrostatic",logical(definition.isHydrostatic),"workload",string(definition.workload),"forcing",string(definition.forcing),"shouldAntialias",logical(definition.shouldAntialias),"integrator",string(definition.requestedIntegrator),"deltaT",double(definition.deltaT),"integrationStepCount",double(definition.integrationStepCount),"finalTime",double(definition.finalTime),"relativeTolerance",double(definition.relativeTolerance),"absoluteToleranceScale",double(definition.absoluteTolerance),"absoluteToleranceEvidence","per-component compact hashes in interface diagnostics","outputRelativeTolerance",double(definition.outputRelativeTolerance),"outputAbsoluteTolerance",double(definition.outputAbsoluteTolerance),"initialStep",double(definition.initialStep),"maximumStep",double(definition.maximumStep),"maximumStepPolicy",optionalString(definition,"maximumStepPolicy","explicit"),"outputInterval",double(definition.outputInterval),"denseOutputPointsPerStep",double(definition.denseOutputPointsPerStep),"denseOutputStartTime",optionalDouble(definition,"denseOutputStartTime"),"denseOutputEndTime",optionalDouble(definition,"denseOutputEndTime"),"denseOutputRecordCount",optionalDouble(definition,"denseOutputRecordCount"),"denseOutputIntegrationRecordCount",optionalDouble(definition,"denseOutputIntegrationRecordCount"),"observerGraph",string(definition.observerGraph),"processRunCount",double(raw.configuration.processRunCount),"warmupCount",double(raw.configuration.warmupCount),"samplesPerProcess",double(raw.configuration.samplesPerProcess));
    graph = comparison.outputGraph;
    graphSummary = struct("passed",logical(graph.passed),"variableCount",double(graph.variableCount),"recordCount",double(graph.recordCount),"maximumAbsoluteError",double(graph.maximumAbsoluteError),"maximumRelativeError",double(graph.maximumRelativeError),"categories",graph.categories);
    correctness = struct("passed",logical(comparison.matchedContractPassed),"maximumRelativeError",double(comparison.maximumRelativeError),"outputAgreementPassed",logical(comparison.outputAgreementPassed),"endpointTrajectoryAgreementPassed",logical(comparison.endpointTrajectoryAgreementPassed),"completeOutputGraph",graphSummary);
    cases{iCase} = struct("id",string(comparison.id),"physicalConfiguration",string(definition.physicalConfiguration),"workload",string(definition.workload),"integrator",string(definition.requestedIntegrator),"contract",contract,"interfaces",{interfaces},"correctness",correctness);
end
archive = struct("fileName",options.archiveFileName,"sha256",options.archiveSHA256,"compressedBytes",options.archiveCompressedBytes,"location","external sibling archive ../wave-vortex-model-benchmark-artifacts/three-interface; not distributed with source");
fixtures = raw.configuration.fixtures;
for iFixture = 1:numel(fixtures)
    if isfield(fixtures,"path"), fixtures(iFixture).path = string(fixtures(iFixture).path); end
end
provenance = struct("rawSchemaVersion",string(raw.schemaVersion),"rawArtifactSHA256",sha256File(rawArtifactPath),"externalArchive",archive,"fixtures",fixtures,"initialCondition",optionalStruct(raw.configuration,"initialCondition"),"stepControls",optionalStruct(raw.configuration,"stepControls"));
if isfield(raw,"revalidation"), provenance.revalidation = raw.revalidation; end
source = struct("repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",string(raw.source.commit),"tree",string(raw.source.tree),"sourceDirty",false,"version",options.implementationVersion);
dataset = struct("schemaVersion","published-three-interface-v3","datasetId",datasetId,"collectedAt",collectedAt,"studyId","integrator-runtime-memory-v1","source",source,"platform",platform,"provider",provider,"provenance",provenance,"cases",{cases});
end

function diagnostics = interfaceDiagnostics(runs)
integrators = [runs.integrator];
methodWork = struct("acceptedStepCounts",double([integrators.acceptedStepCount]),"rejectedStepCounts",double([integrators.rejectedStepCount]),"rhsEvaluationCounts",double([integrators.rhsEvaluationCount]),"denseOutputEvaluationCounts",double([integrators.denseOutputEvaluationCount]),"fsalReuse",sampleField(integrators,"fsalReuseCount"),"fsalInvalidations",sampleField(integrators,"fsalInvalidationCount"));
if isfield(integrators,"continuousExtensionRightHandSideEvaluationCount")
    methodWork.continuousExtensionRightHandSideEvaluationCounts = double([integrators.continuousExtensionRightHandSideEvaluationCount]);
end
controls = struct("requested",string(integrators(1).requested),"actual",string(integrators(1).actual),"matched",all([integrators.matched]));
if isfield(integrators,"controller"), controls.controller = string(integrators(1).controller); end
if isfield(integrators,"relativeTolerance"), controls.relativeTolerance = double(integrators(1).relativeTolerance); end
if isfield(integrators,"requestedInitialStep"), controls.requestedInitialStep = double(integrators(1).requestedInitialStep); controls.effectiveInitialStep = double(integrators(1).effectiveInitialStep); end
if isfield(integrators,"requestedMaximumStep"), controls.requestedMaximumStep = double(integrators(1).requestedMaximumStep); controls.effectiveMaximumStep = double(integrators(1).effectiveMaximumStep); end
if isfield(integrators,"maximumStepPolicy"), controls.maximumStepPolicy = string(integrators(1).maximumStepPolicy); end
if isfield(integrators,"absoluteToleranceHash"), controls.absoluteToleranceHashes = string({integrators.absoluteToleranceHash}); controls.absoluteToleranceComponentHashes = string(integrators(1).absoluteToleranceComponentHashes); end
memory = struct("boundary",string(runs(1).memory.boundary),"steadyRetainedRSSBytes",median(arrayfun(@(run)run.memory.baselineProcessBytes,runs)),"operationPeakIncrementRSSBytes",median(arrayfun(@(run)run.memory.peakIncrementBytes,runs)),"finalRSSBytes",median(arrayfun(@(run)run.memory.finalRSSBytes,runs)),"processLifetimePeakRSSBytes",median(arrayfun(@(run)run.memory.processLifetimePeakRSSBytes,runs)),"allocatorAndProviderAttribution","MATLAB allocator/COW and opaque FFT/provider storage are not exactly attributable from total RSS");
stateSizedBuffers = struct([]);
if isfield(integrators,"stateSizedBuffers"), stateSizedBuffers = integrators(1).stateSizedBuffers; end
storage = integrators(1).storageAccounting;
if isfield(integrators,"sharedAbstractionStateSizedCopyCount")
    storage.sharedAbstractionStateSizedCopyCount = double(integrators(1).sharedAbstractionStateSizedCopyCount);
else
    storage.sharedAbstractionStateSizedCopyCount = [];
end
for name = ["workspaceStateEquivalentCount" "workspaceMaximumLiveStateEquivalentCount" "denseHistoryStateEquivalentCount" "continuousExtensionWorkspaceStateEquivalentCount" "continuousExtensionWorkspaceMaximumLiveStateEquivalentCount"]
    if isfield(integrators,name), storage.(name) = integrators(1).(name); end
end
diagnostics = struct("controls",controls,"methodWork",methodWork,"integratorStorage",storage,"stateSizedBuffers",stateSizedBuffers,"memory",memory);
end

function value = sampleField(structures,name)
values = {structures.(name)};
if all(cellfun(@isnumeric,values)), value = double(cell2mat(values)); else, value = string(values); end
end

function value = optionalString(structure,name,defaultValue)
if isfield(structure,name), value = string(structure.(name)); else, value = string(defaultValue); end
end

function value = optionalDouble(structure,name)
if isfield(structure,name), value = double(structure.(name)); else, value = []; end
end

function value = optionalStruct(structure,name)
if isfield(structure,name), value = structure.(name); else, value = struct(); end
end

function value = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status~=0, error("WaveVortexBenchmark:HashFailed","%s",output); end
value = string(extractBefore(strtrim(output),65));
end

function value = shellQuote(value)
singleQuoteEscape = char([39 34 39 34 39]);
value = "'"+replace(string(value),"'",singleQuoteEscape)+"'";
end

function dataset = normalizeLegacy(raw,options)
validateThreeInterfaceBenchmarkContract(raw);
[collectedAt,timestamp] = collectionTime(raw.runId);
datasetId = "three-interface--"+options.platformId+"--"+timestamp;
platform = struct("id",options.platformId,"displayName",options.platformName,"processor",string(raw.environment.processor),"physicalMemoryBytes",double(raw.environment.physicalMemoryBytes),"os",string(raw.environment.os),"architecture",string(raw.environment.architecture),"matlabVersion",string(raw.environment.matlabVersion),"threadCount",double(raw.configuration.threadCount));
provider = struct("id",string(raw.provider.provider.id),"version",string(raw.provider.provider.version),"threadBackend",string(raw.provider.provider.threadBackend),"scope","compiled-interfaces-only","moduleSHA256",string(raw.provider.module.sha256),"identityValidated",logical(raw.provider.module.identityValidated),"openMPDetected",logical(raw.provider.libraries.openmp.detected));
cases = cell(1,numel(raw.comparison));
for iCase = 1:numel(raw.comparison)
    comparison = raw.comparison(iCase);
    definition = raw.cases(iCase);
    if ~isfield(comparison,"integratorAgreementPassed") || ~logical(comparison.integratorAgreementPassed)
        error("WaveVortexBenchmark:IntegratorMismatch","Publication requires the requested integrator to execute in every interface.");
    end
    if string(definition.requestedIntegrator) == "adaptive-rk23" && (~isfield(comparison,"adaptiveWorkAgreementPassed") || ~logical(comparison.adaptiveWorkAgreementPassed))
        error("WaveVortexBenchmark:AdaptiveWorkMismatch","Publication requires matched adaptive controller and work-count evidence.");
    end
    if ~isfield(comparison,"memoryAgreementPassed") || ~logical(comparison.memoryAgreementPassed)
        error("WaveVortexBenchmark:IncomparableMemory","Publication requires complete process-tree RSS measurements for every interface.");
    end
    if ~isfield(comparison,"outputGraph") || ~logical(comparison.outputGraph.passed)
        error("WaveVortexBenchmark:OutputGraphMismatch","Publication requires agreement across the complete saved output graph.");
    end
    requiredCategories = "coefficients";
    if string(definition.operation) == "model-continuation"
        requiredCategories = [requiredCategories "eulerianFields" "moorings" "particles" "tracers" "times"];
    end
    availableCategories = string({comparison.outputGraph.categories.name});
    if ~all(ismember(requiredCategories,availableCategories)) || any(~[comparison.outputGraph.categories.passed])
        error("WaveVortexBenchmark:OutputGraphMismatch","Publication requires every expected output category to be present and equivalent.");
    end
    caseRuns = raw.runs(reshape(string(arrayfun(@(run)run.case.id,raw.runs,"UniformOutput",false)),[],1)==string(comparison.id));
    if any(arrayfun(@(run)~logical(run.integrator.matched) || string(run.integrator.requested)~=string(definition.requestedIntegrator) || string(run.integrator.actual)~=string(definition.requestedIntegrator),caseRuns))
        error("WaveVortexBenchmark:IntegratorMismatch","Publication found a requested/actual integrator mismatch.");
    end
    interfaces = cell(1,numel(comparison.interfaces));
    for iInterface = 1:numel(comparison.interfaces)
        item = comparison.interfaces(iInterface);
        selected = raw.runs(reshape(string({raw.runs.interface}),[],1)==string(item.id) & reshape(string(arrayfun(@(run)run.case.id,raw.runs,"UniformOutput",false)),[],1)==string(comparison.id));
        if any(arrayfun(@(run)string(run.memory.status)~="complete" || string(run.memory.provider)~="macos-ps-process-tree" || ~isfinite(run.memory.totalPeakRSSBytes),selected))
            error("WaveVortexBenchmark:IncomparableMemory","Publication found a missing or incomparable primary memory measurement.");
        end
        interfaces{iInterface} = struct("id",string(item.id),"providerId",string(selected(1).provider.id),"processWallSeconds",double(item.processWallSeconds),"interfaceTotalSeconds",double(item.interfaceTotalSeconds),"integrationSeconds",double(item.integrationSeconds),"totalPeakRSSBytes",double(item.totalPeakRSSBytes),"incrementalPeakRSSBytes",double(item.incrementalPeakRSSBytes),"finalRSSBytes",double(item.finalRSSBytes),"processWallRatio",double(item.processWallRatio),"integrationRatio",double(item.integrationRatio),"totalRSSRatio",double(item.totalRSSRatio),"incrementalRSSRatio",double(item.incrementalRSSRatio),"processWallSamplesSeconds",double([selected.processWallSeconds]),"integrationSamplesSeconds",double([selected.integrationSeconds]));
    end
    contract = struct("Nxyz",double(definition.Nxyz(:)'),"Lxyz",double(raw.configuration.Lxyz(:)'),"forcing",string(definition.forcing),"shouldAntialias",logical(definition.shouldAntialias),"integrator",string(definition.requestedIntegrator),"deltaT",double(definition.deltaT),"finalTime",double(definition.finalTime),"relativeTolerance",double(definition.relativeTolerance),"absoluteTolerance",double(definition.absoluteTolerance),"outputInterval",double(definition.outputInterval),"observerGraph",string(definition.observerGraph),"processRunCount",double(raw.configuration.processRunCount),"warmupCount",double(raw.configuration.warmupCount),"samplesPerProcess",double(raw.configuration.samplesPerProcess));
    adaptiveWork = struct();
    if string(definition.requestedIntegrator) == "adaptive-rk23"
        exemplar = caseRuns(1).integrator;
        adaptiveWork = rmfield(exemplar,["requested" "actual" "matched"]);
        adaptiveWork.absoluteToleranceFingerprintAgreementPassed = toleranceFingerprintAgreement(comparison,caseRuns);
    end
    graph = comparison.outputGraph;
    graphSummary = struct("passed",logical(graph.passed),"variableCount",double(graph.variableCount),"recordCount",double(graph.recordCount),"maximumAbsoluteError",double(graph.maximumAbsoluteError),"maximumRelativeError",double(graph.maximumRelativeError),"categories",graph.categories);
    correctness = struct("passed",logical(comparison.matchedContractPassed),"maximumRelativeError",double(comparison.maximumRelativeError),"outputAgreementPassed",logical(comparison.outputAgreementPassed),"completeOutputGraph",graphSummary);
    cases{iCase} = struct("id",string(comparison.id),"operation",string(definition.operation),"contract",contract,"adaptiveWork",adaptiveWork,"interfaces",{interfaces},"correctness",correctness);
end
provenance = struct("rawSchemaVersion",string(raw.schemaVersion),"externalArchive",struct("fileName",options.archiveFileName,"sha256",options.archiveSHA256,"compressedBytes",options.archiveCompressedBytes));
source = struct("repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",string(raw.source.commit),"tree",string(raw.source.tree),"sourceDirty",false,"version",options.implementationVersion);
dataset = struct("schemaVersion","published-three-interface-v1","datasetId",datasetId,"collectedAt",collectedAt,"source",source,"platform",platform,"provider",provider,"provenance",provenance,"cases",{cases});
end

function value = toleranceFingerprintAgreement(comparison,runs)
if isfield(comparison,"absoluteToleranceFingerprintAgreementPassed")
    value = logical(comparison.absoluteToleranceFingerprintAgreementPassed);
    return
end
reference = runs(1).integrator;
value = all(arrayfun(@(run)string(run.integrator.absoluteToleranceHash)==string(reference.absoluteToleranceHash),runs));
value = value && all(arrayfun(@(run)isequal(string(run.integrator.absoluteToleranceComponentHashes(:)),string(reference.absoluteToleranceComponentHashes(:))),runs));
end

function [collectedAt,timestamp] = collectionTime(runId)
token = regexp(string(runId),'^(\d{8}T\d{6})\d{0,3}Z$','tokens','once');
if isempty(token), error("WaveVortexBenchmark:InvalidCollectionTime","runId must use UTC YYYYMMDDTHHmmss[SSS]Z."); end
timestamp = string(token{1})+"Z";
instant = datetime(timestamp,"InputFormat","yyyyMMdd'T'HHmmss'Z'","TimeZone","UTC");
collectedAt = string(instant,"yyyy-MM-dd'T'HH:mm:ss'Z'");
end
