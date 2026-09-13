function artifact = measureHydrostaticAdapterBoundary(options)
% Measure the bounded MATLAB/compiled Hydrostatic adapter boundary.
% This is feasibility evidence, not a production performance qualification.
arguments
    options.outputPath (1,1) string = fullfile(tempdir,"wvm503-hydrostatic-adapter-boundary.json")
    options.repeatCount (1,1) double {mustBeInteger,mustBePositive} = 3
end
Lxyz = [150e3 150e3 1300];
Nxyz = [256 256 28];
Nj = 18; % Explicit vertical mode count for this EddyTide-sized feasibility case.
profile = @(z)2e-5*exp(2*z/1300);
rng(503,"twister");
setupTimer = tic;
wvt = WVTransformHydrostatic(Lxyz,Nxyz,Nj=Nj,N2Function=profile,latitude=45,shouldAntialias=true);
wvt.initWithRandomFlow(uvMax=0.01);
wvt.t0 = 37;
wvt.t = 113;
setupSeconds = toc(setupTimer);

createTimer = tic;
backend = WVCompiledTransformBackend.create(wvt);
createSeconds = toc(createTimer);
cleanup = onCleanup(@()deletePair(wvt,backend));

names = ["u" "v" "w" "eta"];
prepareTimer = tic;
backend.prepare(names);
prepareSeconds = toc(prepareTimer);

% Warm both paths before timed samples.
compiledMixed = backend.evaluate(wvt,names,shouldEvaluateFlux=true);
compiledFields = backend.evaluate(wvt,names);
warmMetrics = compiledFields.metrics;
expected = cellstr(names);
matlabFields = cell(1,numel(names));
% Warmup outputs are deliberately discarded to avoid retaining large copies.

mixedSeconds = zeros(options.repeatCount,1);
fieldsSeconds = zeros(options.repeatCount,1);
referenceMixedSeconds = zeros(options.repeatCount,1);
referenceFieldsSeconds = zeros(options.repeatCount,1);
for iRepeat = 1:options.repeatCount
    timer = tic; compiledMixed = backend.evaluate(wvt,names,shouldEvaluateFlux=true); mixedSeconds(iRepeat) = toc(timer);
    timer = tic; compiledFields = backend.evaluate(wvt,names); fieldsSeconds(iRepeat) = toc(timer);
    invalidateMATLABState(wvt);
    timer = tic; [matlabFp,matlabFm,matlabF0] = wvt.nonlinearFlux(); [matlabFields{:}] = wvt.variableWithName(expected{:}); referenceMixedSeconds(iRepeat) = toc(timer);
    invalidateMATLABState(wvt);
    timer = tic; [matlabFields{:}] = wvt.variableWithName(expected{:}); referenceFieldsSeconds(iRepeat) = toc(timer);
end

fieldErrors = zeros(1,numel(names));
for iName = 1:numel(names)
    fieldErrors(iName) = relativeError(compiledFields.values{iName},matlabFields{iName});
end
fluxErrors = [relativeError(compiledMixed.flux{1},matlabFp),relativeError(compiledMixed.flux{2},matlabFm),relativeError(compiledMixed.flux{3},matlabF0)];
metadata = compiledMixed.metrics;
finalMetrics = compiledFields.metrics;
shape = [wvt.Nx wvt.Ny wvt.Nz];
coefficientBytes = 3*wvt.Nj*wvt.Nkl*2*8;
fieldBytes = numel(names)*prod(shape)*8;
fluxBytes = 3*wvt.Nj*wvt.Nkl*2*8;
artifact = struct( ...
    "schemaVersion","wvm503-hydrostatic-adapter-boundary-v1", ...
    "classification","bounded feasibility measurement; not production timing qualification", ...
    "sourceCommit",git("rev-parse HEAD"),"sourceTreeDirty",git("status --short"),"matlabVersion",version,"moduleSHA256",moduleSHA256(backend.capabilities), ...
    "configuration",struct("Lxyz",Lxyz,"Nxyz",Nxyz,"Nj",Nj,"N2","2e-5*exp(2*z/1300)","latitude",45,"shouldAntialias",true,"seed",503,"t0",37,"t",113,"variableNames",names), ...
    "setup",struct("matlabConstructionAndInitializationSeconds",setupSeconds,"adapterCreateSeconds",createSeconds,"adapterPrepareSeconds",prepareSeconds), ...
    "timing",struct("repeatCount",options.repeatCount,"compiledMixedSeconds",mixedSeconds,"compiledFieldsOnlySeconds",fieldsSeconds,"matlabMixedSeconds",referenceMixedSeconds,"matlabFieldsOnlySeconds",referenceFieldsSeconds,"medianCompiledMixedSeconds",median(mixedSeconds),"medianCompiledFieldsOnlySeconds",median(fieldsSeconds),"medianMatlabMixedSeconds",median(referenceMixedSeconds),"medianMatlabFieldsOnlySeconds",median(referenceFieldsSeconds)), ...
    "transfer",struct("stateInputCopyBytes",metadata.stateInputCopyBytes,"coefficientInputBytesPerEvaluation",coefficientBytes,"fieldResultBytesPerEvaluation",fieldBytes,"fluxResultBytesPerEvaluation",fluxBytes), ...
    "storage",struct("sourceBytes",metadata.sourceBytes,"kernelBytes",metadata.kernelBytes,"engineBytes",metadata.engineBytes,"fieldServiceBytes",metadata.fieldServiceBytes,"preparedPlanBytes",metadata.preparedPlanBytes,"stablePreparedPlanBytes",compiledFields.metrics.preparedPlanBytes,"ownershipNotes",["kernelBytes includes sourceBytes" "engineBytes includes kernelBytes"],"combinedScientificExecutionBytesLowerBound",metadata.engineBytes+metadata.fieldServiceBytes+metadata.preparedPlanBytes), ...
    "evaluationAccounting",struct("mixedMetrics",compiledMixed.metrics,"fieldsMetrics",compiledFields.metrics,"producerExecutions",compiledMixed.metrics.producerExecutions,"stateValidations",compiledMixed.metrics.stateValidations,"phasePreparations",compiledMixed.metrics.phasePreparations,"duplicateExecutions",compiledMixed.metrics.duplicateExecutions), ...
    "parity",struct("fieldNames",names,"fieldRelativeErrors",fieldErrors,"fluxRelativeErrors",fluxErrors,"maximumRelativeError",max([fieldErrors fluxErrors]),"allFinite",all(isfinite([fieldErrors fluxErrors])) && all(cellfun(@(x)all(isfinite(x(:))),compiledFields.values)) && all(cellfun(@(x)all(isfinite(x(:))),compiledMixed.flux)),"tolerance",1e-12,"passed",all(isfinite([fieldErrors fluxErrors])) && max([fieldErrors fluxErrors])<=1e-12), ...
    "limitations",["No direct standalone kernel call was added; existing HydrostaticKernelDump remains the independent native-call validation route." "Host was not asserted idle; timing fields are feasibility observations only."]);
writeText(options.outputPath,jsonencode(artifact,PrettyPrint=true));
if ~artifact.parity.passed || metadata.stateInputCopyBytes ~= 0 || compiledFields.metrics.preparedPlanBytes ~= compiledMixed.metrics.preparedPlanBytes || compiledMixed.metrics.duplicateExecutions ~= 0 || finalMetrics.planPreparations ~= warmMetrics.planPreparations || finalMetrics.sourceBytes ~= warmMetrics.sourceBytes || finalMetrics.kernelBytes ~= warmMetrics.kernelBytes || finalMetrics.engineBytes ~= warmMetrics.engineBytes || finalMetrics.fieldServiceBytes ~= warmMetrics.fieldServiceBytes || finalMetrics.stateValidations-warmMetrics.stateValidations ~= 2*options.repeatCount || finalMetrics.phasePreparations-warmMetrics.phasePreparations ~= 2*options.repeatCount
    error("WaveVortexModel:CompiledTransformBoundaryMismatch","Hydrostatic adapter feasibility evidence failed a parity, copy, plan-stability, or duplicate-execution gate; the artifact was written before this error.");
end

function value = moduleSHA256(capabilities)
value = "";
if isfield(capabilities,"module") && isfield(capabilities.module,"sha256")
    value = string(capabilities.module.sha256);
end
end
clear cleanup
end

function invalidateMATLABState(wvt)
wvt.t = wvt.t + 1;
wvt.t = wvt.t - 1;
end

function value = relativeError(actual,expected)
if any(~isfinite(actual(:))) || any(~isfinite(expected(:)))
    value = Inf;
    return
end
value = max(abs(actual(:)-expected(:)))/max(max(abs(expected(:))),realmin);
end

function deletePair(wvt,backend)
if ~isempty(backend)
    delete(backend);
end
if ~isempty(wvt) && isvalid(wvt)
    delete(wvt);
end
end

function value = git(command)
[status,value] = system("git "+command);
if status ~= 0
    value = "unknown";
else
    value = strtrim(value);
end
end

function writeText(pathname,textValue)
folder = fileparts(pathname);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end
fid = fopen(pathname,"w");
cleanup = onCleanup(@()fclose(fid));
fwrite(fid,textValue);
clear cleanup
end
