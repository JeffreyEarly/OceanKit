function validateThreeInterfaceBenchmarkContract(raw)
% Validate the evidence required to publish a matched three-interface run.
requiredInterfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
schemaVersion = string(raw.schemaVersion);
isIntegratorStudy = ismember(schemaVersion,["three-interface-benchmark-v2" "three-interface-benchmark-v3"]);
isMatchedModelStudy = schemaVersion == "three-interface-benchmark-v3";
if ~ismember(schemaVersion,["three-interface-benchmark-v1" "three-interface-benchmark-v2" "three-interface-benchmark-v3"])
    error("WaveVortexBenchmark:MatchedContractFailed","The three-interface benchmark schema is not recognized.");
end
tolerance = double(raw.configuration.correctnessTolerance);
processRunCount = double(raw.configuration.processRunCount);
if isMatchedModelStudy
    validateMatchedModelStudy(raw);
end
validateCompiledProvider(raw);
definitions = raw.cases;
comparisons = raw.comparison;
caseIds = string({definitions.id});
validDefinitions = validCaseDefinitions(definitions,isIntegratorStudy,isMatchedModelStudy);
if isempty(definitions) || numel(definitions) ~= numel(comparisons) || numel(unique(caseIds)) ~= numel(caseIds) || ~validDefinitions
    error("WaveVortexBenchmark:MatchedContractFailed","Publication requires one or more unique recognized benchmark cases.");
end
for iCase = 1:numel(definitions)
    definition = definitions(iCase);
    comparisonIndex = find(string({comparisons.id}) == string(definition.id),1);
    if isempty(comparisonIndex)
        error("WaveVortexBenchmark:MatchedContractFailed","The comparison for %s is missing.",string(definition.id));
    end
    comparison = comparisons(comparisonIndex);
    caseRuns = raw.runs(string(arrayfun(@(run)run.case.id,raw.runs,"UniformOutput",false)) == string(definition.id));
    if numel(caseRuns) ~= processRunCount*numel(requiredInterfaces)
        error("WaveVortexBenchmark:MatchedContractFailed","Case %s does not contain the required complete fresh-process runs.",string(definition.id));
    end
    if isMatchedModelStudy
        validateMatchedAvailability(caseRuns,definition,processRunCount);
    elseif any(string({caseRuns.status}) ~= "complete")
        error("WaveVortexBenchmark:MatchedContractFailed","Case %s does not contain the required complete fresh-process runs.",string(definition.id));
    end
    measuredRuns = caseRuns(string({caseRuns.status})=="complete");
    validateIntegrators(measuredRuns,definition,comparison);
    validateAdaptiveWork(measuredRuns,definition,comparison,isIntegratorStudy);
    validateMemory(measuredRuns,comparison,isIntegratorStudy);
    validateRunProviders(caseRuns,requiredInterfaces,raw,isMatchedModelStudy);
    validateNumerics(comparison,definition,tolerance);
    if isIntegratorStudy
        validateIntegratorEvidence(measuredRuns,definition,comparison);
    end
    if ~logical(comparison.matchedContractPassed)
        error("WaveVortexBenchmark:MatchedContractFailed","Case %s did not pass its matched contract.",string(definition.id));
    end
end
end

function validateMatchedModelStudy(raw)
configuration = raw.configuration;
if string(configuration.studyId)~="matched-model-runtime-v1" || double(configuration.processRunCount)~=3 || ~logical(configuration.publicationEligible)
    error("WaveVortexBenchmark:MatchedModelPublicationContract","Publication requires the complete three-repeat matched-model-runtime-v1 contract.")
end
modelConfiguration = string(raw.modelConfiguration);
if ~isscalar(modelConfiguration) || ~ismember(modelConfiguration,["constant-nonhydrostatic" "hydrostatic-exponential" "boussinesq-exponential"])
    error("WaveVortexBenchmark:MatchedModelPublicationContract","Each artifact must identify exactly one recognized model configuration.")
end
if ~isequal(double(configuration.Nxyz(:)'),[256 256 129]) || ~isequal(double(configuration.Lxyz(:)'),[150e3 150e3 1300])
    error("WaveVortexBenchmark:MatchedModelPublicationContract","The matched model grid and domain do not match the publication contract.")
end
if ~isequal(reshape(string(configuration.integrators),1,[]),"adaptive-rk78") || ~isequal(reshape(string(configuration.workloads),1,[]),["coefficient-endpoint" "composite-dense-output"])
    error("WaveVortexBenchmark:MatchedModelPublicationContract","Publication requires RK78 and exactly the two matched workloads.")
end
model = configuration.model;
requiredModel = ["id" "transformClass" "physicalConfiguration" "isHydrostatic" "domainMeters" "grid" "latitudeDegrees" "shouldAntialias" "stratificationProfile" "N2ReferencePerSecondSquared" "exponentialScaleHeightMeters"];
if ~all(isfield(model,requiredModel)) || string(model.id)~=modelConfiguration || ~isequal(double(model.domainMeters(:)'),[150e3 150e3 1300]) || ~isequal(double(model.grid(:)'),[256 256 129]) || double(model.latitudeDegrees)~=45 || ~logical(model.shouldAntialias) || double(model.N2ReferencePerSecondSquared)~=2e-5
    error("WaveVortexBenchmark:PhysicalProvenance","The artifact omits required physical model provenance.")
end
isExponential = modelConfiguration~="constant-nonhydrostatic";
expectedClass = conditional(modelConfiguration=="constant-nonhydrostatic","WVTransformConstantStratification",conditional(modelConfiguration=="hydrostatic-exponential","WVTransformHydrostatic","WVTransformBoussinesq"));
expectedPhysical = conditional(modelConfiguration=="hydrostatic-exponential","hydrostatic","nonhydrostatic");
expectedProfile = conditional(isExponential,"N2(z) = 2e-5 exp(2 z / 1300) s^-2","N2 = 2e-5 s^-2");
validModel = string(model.transformClass)==expectedClass && string(model.physicalConfiguration)==expectedPhysical && logical(model.isHydrostatic)==(modelConfiguration=="hydrostatic-exponential") && string(model.stratificationProfile)==expectedProfile;
if isExponential
    validModel = validModel && isscalar(model.exponentialScaleHeightMeters) && double(model.exponentialScaleHeightMeters)==650;
else
    validModel = validModel && (isempty(model.exponentialScaleHeightMeters) || (isscalar(model.exponentialScaleHeightMeters) && isnan(model.exponentialScaleHeightMeters)));
end
initial = configuration.initialCondition;
validInitial = string(initial.id)=="gm0p5-red-geostrophic-j1-v1" && double(initial.seed)==4001 && double(initial.gmEnergyLevel)==0.5 && double(initial.geostrophicVerticalMode)==1 && double(initial.geostrophicMaximumSpeedTarget)==0.15;
if ~validModel || ~validInitial
    error("WaveVortexBenchmark:PhysicalProvenance","The artifact physical model or initialization provenance does not match the publication contract.")
end
fixtures = configuration.fixtures;
validFixtures = numel(fixtures)==2 && all(isfield(fixtures,["modelConfiguration" "physicalConfiguration" "stratificationProfile" "workload" "sha256" "bytes"]));
if validFixtures
    validFixtures = all(string({fixtures.modelConfiguration})==modelConfiguration) && isequal(sort(string({fixtures.workload})),sort(["coefficient-endpoint" "composite-dense-output"]));
    validFixtures = validFixtures && all(string({fixtures.physicalConfiguration})==expectedPhysical) && all(string({fixtures.stratificationProfile})==expectedProfile);
    validFixtures = validFixtures && all(strlength(string({fixtures.sha256}))==64) && all(double([fixtures.bytes])>0);
end
if ~validFixtures
    error("WaveVortexBenchmark:PhysicalProvenance","The two matched fixtures lack physical identity, hashes, or measured sizes.")
end
if ~isfield(configuration,"matlabWorker") || strlength(string(configuration.matlabWorker.sha256))~=64 || ~isfield(configuration,"standaloneWorkers") || strlength(string(configuration.standaloneWorkers.runner.sha256))~=64
    error("WaveVortexBenchmark:WorkerIdentity","The current MATLAB worker source fingerprint is missing.")
end
completeRuns = raw.runs(string({raw.runs.status})=="complete");
if strlength(string(raw.source.commit))~=40 || strlength(string(raw.source.tree))~=40 || any(string({completeRuns.sourceCommit})~=string(raw.source.commit))
    error("WaveVortexBenchmark:WorkerIdentity","A measured worker does not identify the artifact source commit.")
end
if ~isfield(configuration,"standaloneWorkers") || strlength(string(configuration.standaloneWorkers.runner.sha256))~=64
    error("WaveVortexBenchmark:WorkerIdentity","The frozen standalone runner fingerprint is missing.")
end
completeRuns = raw.runs(string({raw.runs.status})=="complete");
if strlength(string(raw.source.commit))~=40 || strlength(string(raw.source.tree))~=40 || any(string({completeRuns.sourceCommit})~=string(raw.source.commit))
    error("WaveVortexBenchmark:WorkerIdentity","A measured worker source identity does not match the artifact source.")
end
definitions = raw.cases;
if numel(definitions)~=2 || ~isequal(sort(string({definitions.workload})),sort(["coefficient-endpoint" "composite-dense-output"]))
    error("WaveVortexBenchmark:MatchedModelPublicationContract","Each matched-model artifact requires exactly the two benchmark workloads.")
end
for definition = reshape(definitions,1,[])
    expectedSchedule = conditional(string(definition.workload)=="composite-dense-output",[0 32 64 96 128],[0 7168]);
    valid = string(definition.modelConfiguration)==modelConfiguration && string(definition.requestedIntegrator)=="adaptive-rk78" && double(definition.finalTime)==7168;
    valid = valid && isequal(double(definition.outputScheduleSeconds(:)'),expectedSchedule) && double(definition.seed)==4001;
    if ~valid
        error("WaveVortexBenchmark:MatchedModelPublicationContract","A case changed the frozen model, RK78, duration, seed, or output schedule.")
    end
end
end

function validateMatchedAvailability(caseRuns,definition,processRunCount)
modelConfiguration = string(definition.modelConfiguration);
for interface = ["matlab-builtin" "matlab-compiled" "standalone-compiled"]
    selected = caseRuns(string({caseRuns.interface})==interface);
    if numel(selected)~=processRunCount
        error("WaveVortexBenchmark:InterfaceAvailability","Interface %s does not contain three explicit run records.",interface)
    end
    valid = all(string({selected.status})=="complete");
    if ~valid
        error("WaveVortexBenchmark:InterfaceAvailability","Interface %s has invalid completion or unavailability evidence for %s.",interface,modelConfiguration)
    end
end
end

function valid = validCaseDefinitions(definitions,isIntegratorStudy,isMatchedModelStudy)
if ~isIntegratorStudy
    valid = all(ismember(string({definitions.id}),["nonlinear-flux" "fixed-rk4-continuation" "adaptive-rk23-observer-output"]));
    return
end
required = ["physicalConfiguration" "isHydrostatic" "workload" "requestedIntegrator" "integrationStepCount" "denseOutputPointsPerStep" "outputRelativeTolerance" "outputAbsoluteTolerance"];
if isMatchedModelStudy, required(end+1) = "modelConfiguration"; end
valid = all(arrayfun(@(definition)all(isfield(definition,required)),definitions));
if ~valid, return, end
integrators = ["fixed-rk4" "adaptive-rk23" "adaptive-rk45" "adaptive-rk78"];
workloads = ["coefficient-endpoint" "composite-dense-output"];
physicalConfigurations = ["hydrostatic" "nonhydrostatic"];
valid = all(ismember(string({definitions.requestedIntegrator}),integrators));
valid = valid && all(ismember(string({definitions.workload}),workloads));
valid = valid && all(ismember(string({definitions.physicalConfiguration}),physicalConfigurations));
for definition = reshape(definitions,1,[])
    if isMatchedModelStudy
        expectedId = string(definition.modelConfiguration)+"--"+string(definition.requestedIntegrator)+"--"+string(definition.workload);
    else
        expectedId = string(definition.physicalConfiguration)+"--"+string(definition.requestedIntegrator)+"--"+string(definition.workload);
    end
    valid = valid && string(definition.id)==expectedId && logical(definition.isHydrostatic)==(string(definition.physicalConfiguration)=="hydrostatic");
    if isMatchedModelStudy
        modelConfigurations = ["constant-nonhydrostatic" "hydrostatic-exponential" "boussinesq-exponential"];
        valid = valid && ismember(string(definition.modelConfiguration),modelConfigurations);
        valid = valid && ((string(definition.modelConfiguration)=="hydrostatic-exponential" && string(definition.physicalConfiguration)=="hydrostatic") || (string(definition.modelConfiguration)~="hydrostatic-exponential" && string(definition.physicalConfiguration)=="nonhydrostatic"));
    end
end
end


function validateAdaptiveWork(caseRuns,definition,comparison,isIntegratorStudy)
if ~startsWith(string(definition.requestedIntegrator),"adaptive-")
    return
end
if ~isfield(comparison,"adaptiveWorkAgreementPassed") || ~logical(comparison.adaptiveWorkAgreementPassed)
    error("WaveVortexBenchmark:AdaptiveWorkMismatch","Publication requires matched adaptive controls and complete method-work evidence in every interface.");
end
if isIntegratorStudy && (~isfield(comparison,"absoluteToleranceFingerprintAgreementPassed") || ~logical(comparison.absoluteToleranceFingerprintAgreementPassed))
    error("WaveVortexBenchmark:AdaptiveToleranceMismatch","Publication requires identical quantized component absolute-tolerance fingerprints in every interface.");
end
% Legacy v1 evidence keeps fingerprint disagreement diagnostic because those
% immutable artifacts predate the strict integrator-study control contract.
required = ["controller" "relativeTolerance" "absoluteToleranceHash" "absoluteToleranceHashClearedMantissaBits" "absoluteToleranceComponentHashes" "requestedInitialStep" "effectiveInitialStep" "requestedMaximumStep" "effectiveMaximumStep" "initialTime" "finalTime" "acceptedStepCount" "rejectedStepCount" "rhsEvaluationCount" "denseOutputEvaluationCount" "fsalReuseCount" "fsalInvalidationCount" "outputRecordCounts"];
if any(arrayfun(@(run)~all(isfield(run.integrator,required)),caseRuns))
    error("WaveVortexBenchmark:AdaptiveWorkMismatch","An adaptive interface omitted required controller or work-count evidence.");
end
end

function validateCompiledProvider(raw)
provider = raw.provider;
valid = string(provider.status) == "available" && logical(provider.isAvailable);
valid = valid && string(provider.provider.id) == "native-neon-pthreads" && string(provider.provider.threadBackend) == "pthreads";
valid = valid && logical(provider.module.identityValidated) && ~logical(provider.libraries.openmp.detected);
valid = valid && double(provider.contract.threadCount) == double(raw.configuration.threadCount);
valid = valid && string(provider.libraries.base.path) ~= "" && string(provider.libraries.thread.path) ~= "";
valid = valid && isfinite(provider.featureValidation.maximumRelativeError) && provider.featureValidation.maximumRelativeError <= raw.configuration.correctnessTolerance;
if ~valid
    error("WaveVortexBenchmark:ProviderMismatch","Publication requires the validated native FFTW provider for both compiled interfaces.");
end
end

function validateIntegrators(caseRuns,definition,comparison)
expected = string(definition.requestedIntegrator);
valid = logical(comparison.integratorAgreementPassed);
valid = valid && all(arrayfun(@(run)logical(run.integrator.matched) && string(run.integrator.requested) == expected && string(run.integrator.actual) == expected,caseRuns));
if ~valid
    error("WaveVortexBenchmark:IntegratorMismatch","Publication requires the requested integrator to execute in every interface.");
end
end

function validateMemory(caseRuns,comparison,isIntegratorStudy)
valid = logical(comparison.memoryAgreementPassed);
valid = valid && all(arrayfun(@(run)validMemoryRecord(run,isIntegratorStudy),caseRuns));
if ~valid
    error("WaveVortexBenchmark:IncomparableMemory","Publication requires complete process-tree RSS measurements for every interface.");
end
end

function valid = validMemoryRecord(run,isIntegratorStudy)
memory = run.memory;
valid = string(memory.status) == "complete" && string(memory.provider) == "macos-ps-process-tree";
valid = valid && isfinite(memory.totalPeakRSSBytes) && memory.totalPeakRSSBytes > 0;
valid = valid && isfinite(memory.peakIncrementBytes) && memory.peakIncrementBytes >= 0;
valid = valid && isfinite(memory.finalRSSBytes) && memory.finalRSSBytes >= 0;
if isIntegratorStudy
    valid = valid && isfield(memory,"boundary") && string(memory.boundary)=="integration-phase-total-live-process-tree-rss";
    valid = valid && isfield(memory,"processLifetimePeakRSSBytes") && isfinite(memory.processLifetimePeakRSSBytes) && memory.processLifetimePeakRSSBytes>=memory.totalPeakRSSBytes;
    valid = valid && isfield(memory,"integrationSampleCount") && double(memory.integrationSampleCount)>0;
end
end

function validateRunProviders(caseRuns,requiredInterfaces,raw,isMatchedModelStudy)
for interface = requiredInterfaces
    selected = caseRuns(string({caseRuns.interface}) == interface);
    if numel(selected) ~= double(raw.configuration.processRunCount)
        error("WaveVortexBenchmark:ProviderMismatch","Interface %s does not contain the required provider records.",interface);
    end
    for iRun = 1:numel(selected)
        if string(selected(iRun).status)=="unavailable"
            continue
        end
        provider = selected(iRun).provider;
        valid = logical(provider.noFallback);
        if interface == "matlab-builtin"
            valid = valid && double(provider.threads) == double(raw.configuration.threadCount);
            valid = valid && string(provider.id) == "matlab-builtin" && string(provider.baseLibrary) == "" && string(provider.threadLibrary) == "";
        else
            valid = valid && string(provider.id) == string(raw.provider.provider.id) && string(provider.version) == string(raw.provider.provider.version);
            valid = valid && samePath(provider.baseLibrary,raw.provider.libraries.base.path) && samePath(provider.threadLibrary,raw.provider.libraries.thread.path);
            if isMatchedModelStudy
                valid = valid && validMatchedCompiledPolicy(provider,raw.configuration.threadCount,string(raw.modelConfiguration),interface);
            else
                valid = valid && double(provider.threads) == double(raw.configuration.threadCount);
            end
        end
        if ~valid
            error("WaveVortexBenchmark:ProviderMismatch","Interface %s did not execute its required transform provider without fallback.",interface);
        end
        if isMatchedModelStudy && startsWith(interface,"matlab-")
            validWorker = isfield(selected(iRun),"worker") && isfield(selected(iRun).worker,"sha256") && string(selected(iRun).worker.sha256)==string(raw.configuration.matlabWorker.sha256);
            if ~validWorker
                error("WaveVortexBenchmark:WorkerIdentity","A MATLAB run did not execute the frozen current worker source.")
            end
        end
    end
end
if isMatchedModelStudy && ismember(string(raw.modelConfiguration),["hydrostatic-exponential" "boussinesq-exponential"])
    matlabCompiled = caseRuns(string({caseRuns.interface})=="matlab-compiled");
    standaloneCompiled = caseRuns(string({caseRuns.interface})=="standalone-compiled");
    policyFields = ["policy" "matrixBackend" "compact" "effectiveFFTThreads" "horizontalWorkers" "pointwiseWorkers"];
    for field = policyFields
        matlabValues = arrayfun(@(run)string(run.provider.(field)),matlabCompiled);
        standaloneValues = arrayfun(@(run)string(run.provider.(field)),standaloneCompiled);
        if numel(unique(matlabValues))~=1 || numel(unique(standaloneValues))~=1 || matlabValues(1)~=standaloneValues(1)
            error("WaveVortexBenchmark:ProviderMismatch","The two compiled variable-model interfaces did not execute the same qualified native policy.")
        end
    end
end
end

function valid = validMatchedCompiledPolicy(provider,configuredThreadBudget,modelConfiguration,interface)
required = ["configuredThreadBudget" "requestedFFTThreads" "effectiveFFTThreads" "policy" "matrixBackend" "compact" "horizontalWorkers" "pointwiseWorkers"];
valid = all(isfield(provider,required)) && double(provider.configuredThreadBudget)==double(configuredThreadBudget);
if ~valid, return, end
valid = valid && double(provider.threads)==double(provider.effectiveFFTThreads);
if ismember(modelConfiguration,["hydrostatic-exponential" "boussinesq-exponential"])
    valid = valid && isempty(provider.requestedFFTThreads) && double(provider.effectiveFFTThreads)==1;
    valid = valid && string(provider.policy)=="compact-native-accelerate" && string(provider.matrixBackend)=="accelerate" && logical(provider.compact);
    valid = valid && double(provider.horizontalWorkers)>0 && double(provider.pointwiseWorkers)>0;
else
    valid = valid && double(provider.effectiveFFTThreads)==double(configuredThreadBudget);
    valid = valid && string(provider.policy)=="non-variable-transform";
    if interface=="matlab-compiled"
        valid = valid && isempty(provider.requestedFFTThreads);
    else
        valid = valid && double(provider.requestedFFTThreads)==double(configuredThreadBudget);
    end
end
end

function validateNumerics(comparison,definition,tolerance)
usesMethodTolerance = isfield(definition,"outputRelativeTolerance");
if ~isfinite(comparison.maximumRelativeError) || (~usesMethodTolerance && comparison.maximumRelativeError > tolerance)
    error("WaveVortexBenchmark:NumericalMismatch","Case %s exceeds the relative-error tolerance %.3e.",string(definition.id),tolerance);
end
if ~logical(comparison.outputAgreementPassed) || ~logical(comparison.outputGraph.passed)
    error("WaveVortexBenchmark:OutputGraphMismatch","Case %s does not agree across its complete output graph.",string(definition.id));
end
graph = comparison.outputGraph;
finiteGraphEvidence = isfinite(graph.maximumRelativeError) && isfinite(graph.maximumAbsoluteError) && all(isfinite([graph.categories.maximumRelativeError])) && all(isfinite([graph.categories.maximumAbsoluteError]));
if usesMethodTolerance
    validTolerances = isfinite(definition.outputRelativeTolerance) && definition.outputRelativeTolerance>0 && isfinite(definition.outputAbsoluteTolerance) && definition.outputAbsoluteTolerance>=0;
    numericalAgreement = validTolerances && finiteGraphEvidence && all([graph.categories.passed]);
else
    numericalAgreement = finiteGraphEvidence && graph.maximumRelativeError<=tolerance && all([graph.categories.passed]) && all([graph.categories.maximumRelativeError]<=tolerance);
end
if ~numericalAgreement
    error("WaveVortexBenchmark:OutputGraphMismatch","Case %s contains an output payload outside the relative-error tolerance.",string(definition.id));
end
requiredCategories = "coefficients";
if isfield(definition,"workload") && string(definition.workload)=="coefficient-endpoint"
    requiredCategories = [requiredCategories "times"];
elseif string(definition.operation) == "model-continuation"
    requiredCategories = [requiredCategories "eulerianFields" "moorings" "particles" "tracers" "times"];
end
if ~all(ismember(requiredCategories,string({graph.categories.name})))
    missingCategories = requiredCategories(~ismember(requiredCategories,string({graph.categories.name})));
    error("WaveVortexBenchmark:OutputGraphMismatch","Case %s is missing required output categories: %s.",string(definition.id),strjoin(missingCategories,", "));
end
end

function validateIntegratorEvidence(caseRuns,definition,comparison)
if ~isfield(comparison,"endpointTrajectoryAgreementPassed") || ~logical(comparison.endpointTrajectoryAgreementPassed)
    error("WaveVortexBenchmark:EndpointTrajectoryMismatch","Adding scheduled interior output changed an accepted endpoint trajectory.");
end
standalone = caseRuns(string({caseRuns.interface})=="standalone-compiled");
for run = reshape(standalone,1,[])
    integrator = run.integrator;
    required = ["acceptedStepCount" "rejectedStepCount" "rhsEvaluationCount" "denseOutputEvaluationCount" "fsalReuseCount" "fsalInvalidationCount" "workspaceStateEquivalentCount" "workspaceMaximumLiveStateEquivalentCount" "denseHistoryStateEquivalentCount" "continuousExtensionRightHandSideEvaluationCount" "continuousExtensionWorkspaceStateEquivalentCount" "continuousExtensionWorkspaceMaximumLiveStateEquivalentCount" "stateSizedBuffers" "storageAccounting" "sharedAbstractionStateSizedCopyCount"];
    if ~all(isfield(integrator,required))
        error("WaveVortexBenchmark:IncompleteIntegratorEvidence","The standalone integrator omitted required work or storage evidence.");
    end
    storage = integrator.storageAccounting;
    validStorage = isfield(storage,"exact") && logical(storage.exact) && isfield(storage,"byteLedgerAgreement") && logical(storage.byteLedgerAgreement);
    validStorage = validStorage && double(storage.persistentBytes)>0 && double(storage.workspaceCapacityBytes)>0 && double(storage.workspaceMaximumLiveBytes)>=double(storage.workspaceCapacityBytes);
    validStorage = validStorage && isfield(storage,"continuousExtensionBytes") && isfield(storage,"sharedAbstractionStateSizedCopyCount") && double(storage.sharedAbstractionStateSizedCopyCount)==0;
    if ~validStorage
        error("WaveVortexBenchmark:WorkspaceLedgerMismatch","The exact integrator byte ledger is incomplete or inconsistent.");
    end
    buffers = integrator.stateSizedBuffers;
    if isempty(buffers) || any(arrayfun(@(buffer)strlength(string(buffer.buffer))==0 || strlength(string(buffer.producer))==0 || strlength(string(buffer.lastUse))==0,buffers))
        error("WaveVortexBenchmark:WorkspaceLedgerMismatch","Every state-sized buffer requires a producer and last consumer.");
    end
    if double(integrator.sharedAbstractionStateSizedCopyCount)~=0
        error("WaveVortexBenchmark:AbstractionStateCopy","The shared integration abstraction introduced state-sized storage.");
    end
    isEndpoint = string(definition.workload)=="coefficient-endpoint";
    if isEndpoint
        endpointIsLazy = double(integrator.denseOutputEvaluationCount)==0 && double(integrator.denseHistoryStateEquivalentCount)==0 && double(integrator.continuousExtensionRightHandSideEvaluationCount)==0 && double(integrator.continuousExtensionWorkspaceStateEquivalentCount)==0 && double(integrator.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount)==0 && double(storage.denseHistoryBytes)==0;
        if ~endpointIsLazy
            error("WaveVortexBenchmark:UnexpectedDenseOutputWork","Endpoint-only execution allocated or evaluated dense-output state.");
        end
    else
        if double(integrator.denseOutputEvaluationCount)<=0 || double(integrator.denseHistoryStateEquivalentCount)<=0
            error("WaveVortexBenchmark:MissingDenseOutputWork","The composite workload did not exercise method-owned dense output.");
        end
        if string(definition.requestedIntegrator)=="adaptive-rk78" && (double(integrator.continuousExtensionRightHandSideEvaluationCount)<=0 || double(integrator.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount)~=4)
            error("WaveVortexBenchmark:MissingRK78LazyExtension","RK78 did not construct its four lazy continuous-extension buffers only when requested.");
        end
    end
end
for run = reshape(caseRuns(string({caseRuns.interface})~="standalone-compiled"),1,[])
    if ~isfield(run.integrator,"storageAccounting") || logical(run.integrator.storageAccounting.exact)
        error("WaveVortexBenchmark:MatlabStorageAttribution","MATLAB solver/allocator storage must remain explicitly opaque rather than be presented as exact.");
    end
    if string(definition.workload)=="coefficient-endpoint" && double(run.integrator.denseOutputEvaluationCount)~=0
        error("WaveVortexBenchmark:UnexpectedDenseOutputWork","MATLAB endpoint-only execution reported dense-output work.");
    end
end
end

function value = samePath(left,right)
value = string(left) == string(right);
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
