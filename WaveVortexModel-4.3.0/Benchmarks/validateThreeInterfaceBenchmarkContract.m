function validateThreeInterfaceBenchmarkContract(raw)
% Validate the evidence required to publish a matched three-interface run.
requiredInterfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
schemaVersion = string(raw.schemaVersion);
isIntegratorStudy = schemaVersion == "three-interface-benchmark-v2";
if ~ismember(schemaVersion,["three-interface-benchmark-v1" "three-interface-benchmark-v2"])
    error("WaveVortexBenchmark:MatchedContractFailed","The three-interface benchmark schema is not recognized.");
end
tolerance = double(raw.configuration.correctnessTolerance);
processRunCount = double(raw.configuration.processRunCount);
validateCompiledProvider(raw);
definitions = raw.cases;
comparisons = raw.comparison;
caseIds = string({definitions.id});
validDefinitions = validCaseDefinitions(definitions,isIntegratorStudy);
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
    if numel(caseRuns) ~= processRunCount*numel(requiredInterfaces) || any(string({caseRuns.status}) ~= "complete")
        error("WaveVortexBenchmark:MatchedContractFailed","Case %s does not contain the required complete fresh-process runs.",string(definition.id));
    end
    validateIntegrators(caseRuns,definition,comparison);
    validateAdaptiveWork(caseRuns,definition,comparison,isIntegratorStudy);
    validateMemory(caseRuns,comparison,isIntegratorStudy);
    validateRunProviders(caseRuns,requiredInterfaces,raw);
    validateNumerics(comparison,definition,tolerance);
    if isIntegratorStudy
        validateIntegratorEvidence(caseRuns,definition,comparison);
    end
    if ~logical(comparison.matchedContractPassed)
        error("WaveVortexBenchmark:MatchedContractFailed","Case %s did not pass its matched contract.",string(definition.id));
    end
end
end

function valid = validCaseDefinitions(definitions,isIntegratorStudy)
if ~isIntegratorStudy
    valid = all(ismember(string({definitions.id}),["nonlinear-flux" "fixed-rk4-continuation" "adaptive-rk23-observer-output"]));
    return
end
required = ["physicalConfiguration" "isHydrostatic" "workload" "requestedIntegrator" "integrationStepCount" "denseOutputPointsPerStep" "outputRelativeTolerance" "outputAbsoluteTolerance"];
valid = all(arrayfun(@(definition)all(isfield(definition,required)),definitions));
if ~valid, return, end
integrators = ["fixed-rk4" "adaptive-rk23" "adaptive-rk45" "adaptive-rk78"];
workloads = ["coefficient-endpoint" "composite-dense-output"];
physicalConfigurations = ["hydrostatic" "nonhydrostatic"];
valid = all(ismember(string({definitions.requestedIntegrator}),integrators));
valid = valid && all(ismember(string({definitions.workload}),workloads));
valid = valid && all(ismember(string({definitions.physicalConfiguration}),physicalConfigurations));
for definition = reshape(definitions,1,[])
    expectedId = string(definition.physicalConfiguration)+"--"+string(definition.requestedIntegrator)+"--"+string(definition.workload);
    valid = valid && string(definition.id)==expectedId && logical(definition.isHydrostatic)==(string(definition.physicalConfiguration)=="hydrostatic");
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

function validateRunProviders(caseRuns,requiredInterfaces,raw)
for interface = requiredInterfaces
    selected = caseRuns(string({caseRuns.interface}) == interface);
    if numel(selected) ~= double(raw.configuration.processRunCount)
        error("WaveVortexBenchmark:ProviderMismatch","Interface %s does not contain the required provider records.",interface);
    end
    for iRun = 1:numel(selected)
        provider = selected(iRun).provider;
        valid = logical(provider.noFallback) && double(provider.threads) == double(raw.configuration.threadCount);
        if interface == "matlab-builtin"
            valid = valid && string(provider.id) == "matlab-builtin" && string(provider.baseLibrary) == "" && string(provider.threadLibrary) == "";
        else
            valid = valid && string(provider.id) == string(raw.provider.provider.id) && string(provider.version) == string(raw.provider.provider.version);
            valid = valid && samePath(provider.baseLibrary,raw.provider.libraries.base.path) && samePath(provider.threadLibrary,raw.provider.libraries.thread.path);
        end
        if ~valid
            error("WaveVortexBenchmark:ProviderMismatch","Interface %s did not execute its required transform provider without fallback.",interface);
        end
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
