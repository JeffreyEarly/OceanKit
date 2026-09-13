function report = issue513VerticalPrimitives(fixturePath,jsonPath,options)
% Time isolated public derivatives and layout permutations for issue 513.
arguments
    fixturePath (1,1) string {mustBeFile}
    jsonPath (1,1) string
    options.repositoryRoot (1,1) string = ""
    options.dependencyPaths (1,:) string = strings(1,0)
    options.repeatCount (1,1) double {mustBeInteger,mustBePositive} = 5
    options.phasePath (1,1) string = ""
end

repositoryRoot = options.repositoryRoot;
if repositoryRoot == "", repositoryRoot = string(fileparts(fileparts(mfilename("fullpath")))); end
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addRepositoryPaths(repositoryRoot,options.dependencyPaths);
if isfile(jsonPath)
    error("WaveVortexBenchmark:Issue513PrimitiveOutputExists","Refusing to overwrite the diagnostic report: %s",jsonPath)
end

[wvt,reader] = WVTransform.waveVortexTransformFromFile(fixturePath,iTime=Inf,shouldReadOnly=true,computationalBackend="compiled");
resourceCleanup = onCleanup(@()closeResources(reader,wvt));
reader.close();
metadata = validateAndRecordIdentity(wvt);
input = wvt.u;
if ~isreal(input) || ~isequal(size(input),wvt.spatialMatrixSize) || ~all(isfinite(input(:))) || norm(input(:)) == 0
    error("WaveVortexBenchmark:Issue513PrimitiveInput","The fixture did not produce a finite, nonzero real full-grid u field.")
end

metricsBeforeWarmup = metadata.runtimeMetrics;
writePhase(options.phasePath,"warmup");
warmDiffX = wvt.diffX(input);
warmDiffZF = wvt.diffZF(input);
warmInputLayout = reshape(permute(input,[3 1 2]),wvt.Nz,[]);
warmOutputSource = reshape(permute(warmDiffZF,[3 1 2]),wvt.Nz,[]);
warmOutputLayout = permute(reshape(warmOutputSource,wvt.Nz,wvt.Nx,wvt.Ny),[2 3 1]);
requireFinite(warmDiffX,"warm diffX");
requireFinite(warmDiffZF,"warm diffZF");
metricsAfterWarmup = wvt.computationalBackendMetadata.runtimeMetrics;

diffXSeconds = zeros(1,options.repeatCount);
diffZFSeconds = zeros(1,options.repeatCount);
inputPermuteSeconds = zeros(1,options.repeatCount);
outputPermuteSeconds = zeros(1,options.repeatCount);
diffZFRepeatRelativeError = zeros(1,options.repeatCount);

writePhase(options.phasePath,"measure-public-diffX");
for iRepeat = 1:options.repeatCount
    timer = tic;
    diffXOutput = wvt.diffX(input);
    diffXSeconds(iRepeat) = toc(timer);
    requireFinite(diffXOutput,"timed diffX");
end

writePhase(options.phasePath,"measure-public-diffZF");
for iRepeat = 1:options.repeatCount
    timer = tic;
    diffZFOutput = wvt.diffZF(input);
    diffZFSeconds(iRepeat) = toc(timer);
    requireFinite(diffZFOutput,"timed diffZF");
    diffZFRepeatRelativeError(iRepeat) = norm(diffZFOutput(:)-warmDiffZF(:))/max(norm(warmDiffZF(:)),realmin);
end

writePhase(options.phasePath,"measure-input-permute");
for iRepeat = 1:options.repeatCount
    timer = tic;
    inputLayout = reshape(permute(input,[3 1 2]),wvt.Nz,[]);
    inputPermuteSeconds(iRepeat) = toc(timer);
end

writePhase(options.phasePath,"measure-output-permute");
for iRepeat = 1:options.repeatCount
    timer = tic;
    outputLayout = permute(reshape(warmOutputSource,wvt.Nz,wvt.Nx,wvt.Ny),[2 3 1]);
    outputPermuteSeconds(iRepeat) = toc(timer);
end
if ~isequal(size(inputLayout),[wvt.Nz wvt.Nx*wvt.Ny]) || ~isequal(size(outputLayout),size(input)) || ~isequal(warmInputLayout,inputLayout) || ~isequal(warmOutputLayout,warmDiffZF)
    error("WaveVortexBenchmark:Issue513PrimitiveLayout","The isolated layout diagnostic produced an unexpected shape or ordering.")
end

finalMetadata = validateAndRecordIdentity(wvt);
metricsAfterMeasurement = finalMetadata.runtimeMetrics;
report = struct( ...
    "schemaVersion","wvm-issue513-vertical-primitives-v1", ...
    "fixture",struct("path",fixturePath,"sha256",sha256File(fixturePath)), ...
    "source",struct("repositoryRoot",repositoryRoot,"input","wvt.u","inputIsReal",isreal(input)), ...
    "dimensions",struct("Nxyz",[wvt.Nx wvt.Ny wvt.Nz],"spatialMatrixSize",size(input),"horizontalRowCount",wvt.Nx*wvt.Ny,"elementCount",numel(input)), ...
    "identity",rmfield(finalMetadata,"runtimeMetrics"), ...
    "repeatCount",options.repeatCount, ...
    "timingScope","isolated diagnostic timings; public derivative and permutation samples are separate and must not be summed as an integration runtime", ...
    "timings",struct("publicDiffXSeconds",diffXSeconds,"publicDiffZFSeconds",diffZFSeconds,"inputPermuteSeconds",inputPermuteSeconds,"outputPermuteSeconds",outputPermuteSeconds), ...
    "repeatComparison",struct("operator","diffZF","reference","one untimed warmed repeat","relativeErrors",diffZFRepeatRelativeError,"isNumericalCorrectnessGate",false), ...
    "metricsBeforeWarmup",metricsBeforeWarmup, ...
    "metricsAfterWarmup",metricsAfterWarmup, ...
    "metricsAfterMeasurement",metricsAfterMeasurement, ...
    "measurementMetricDeltas",numericScalarDeltas(metricsAfterWarmup,metricsAfterMeasurement), ...
    "measurementCounters",selectCounterMetrics(metricsAfterMeasurement), ...
    "measurementCounterDeltas",selectCounterMetrics(numericScalarDeltas(metricsAfterWarmup,metricsAfterMeasurement)));
writeText(jsonPath,jsonencode(report,PrettyPrint=true));
writePhase(options.phasePath,"complete");
closeResources(reader,wvt);
reader = [];
wvt = [];
clear resourceCleanup pathCleanup
end

function value = validateAndRecordIdentity(wvt)
metadata = wvt.computationalBackendMetadata;
if string(metadata.activeBackend) ~= "compiled" || string(metadata.provider.id) ~= "native-neon-pthreads" || ~metadata.module.identityValidated || metadata.libraries.openmp.detected
    error("WaveVortexBenchmark:Issue513PrimitiveProvider","The transform does not own an identity-validated native compiled backend.")
end
if string(metadata.module.sha256) ~= string(metadata.sourceIdentity.moduleSHA256) || sha256File(metadata.module.path) ~= string(metadata.module.sha256)
    error("WaveVortexBenchmark:Issue513PrimitiveIdentity","The live module, validated source identity, and module file hashes do not agree.")
end
[~,mexFiles] = inmem("-completenames");
loadedPaths = string(mexFiles);
loadedPaths = loadedPaths(endsWith(loadedPaths,filesep+string(metadata.module.name)+"."+mexext));
resolvedPath = string(which(metadata.module.name));
if numel(loadedPaths) ~= 1 || canonicalPath(loadedPaths) ~= canonicalPath(metadata.module.path) || canonicalPath(resolvedPath) ~= canonicalPath(metadata.module.path)
    error("WaveVortexBenchmark:Issue513PrimitiveModulePath","The live or resolved MEX path does not match the validated cached module.")
end
value = struct( ...
    "activeBackend",string(metadata.activeBackend), ...
    "providerId",string(metadata.provider.id), ...
    "sourceManifestSHA256",string(metadata.sourceIdentity.manifest.aggregateSHA256), ...
    "sourceModuleSHA256",string(metadata.sourceIdentity.moduleSHA256), ...
    "moduleName",string(metadata.module.name), ...
    "modulePath",string(metadata.module.path), ...
    "moduleSHA256",string(metadata.module.sha256), ...
    "moduleIdentityValidated",logical(metadata.module.identityValidated), ...
    "loadedModulePath",loadedPaths, ...
    "resolvedModulePath",resolvedPath, ...
    "baseLibraryPath",string(metadata.libraries.base.path), ...
    "baseLibrarySHA256",string(metadata.libraries.base.sha256), ...
    "threadLibraryPath",string(metadata.libraries.thread.path), ...
    "threadLibrarySHA256",string(metadata.libraries.thread.sha256), ...
    "openMPDetected",logical(metadata.libraries.openmp.detected), ...
    "runtimeMetrics",metadata.runtimeMetrics);
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
selected = contains(lower(names),["producer" "primitive" "copy" "execution"]);
for name = reshape(names(any(selected,2)),1,[]), value.(name) = metrics.(name); end
end

function requireFinite(value,label)
if ~all(isfinite(value(:)))
    error("WaveVortexBenchmark:Issue513PrimitiveFinite","%s produced nonfinite values.",label)
end
end

function addRepositoryPaths(repositoryRoot,dependencyPaths)
for dependencyPath = reshape(dependencyPaths,1,[])
    if ~isfolder(dependencyPath), error("WaveVortexBenchmark:Issue513PrimitiveDependency","Dependency folder does not exist: %s",dependencyPath); end
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

function writePhase(pathname,phase)
if pathname == "", return, end
writeText(pathname,phase,shouldOverwrite=true);
end

function writeText(pathname,contents,options)
arguments
    pathname (1,1) string
    contents (1,1) string
    options.shouldOverwrite (1,1) logical = false
end
if isfile(pathname) && ~options.shouldOverwrite
    error("WaveVortexBenchmark:Issue513PrimitiveOutputExists","Refusing to overwrite diagnostic output: %s",pathname)
end
parent = string(fileparts(pathname));
if parent ~= "" && ~isfolder(parent), mkdir(parent); end
fileId = fopen(pathname,"w");
if fileId < 0, error("WaveVortexBenchmark:Issue513PrimitiveWrite","Unable to write %s.",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",contents);
clear cleanup
end

function value = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status ~= 0, error("WaveVortexBenchmark:Issue513PrimitiveHash","Unable to hash %s.",pathname); end
value = extractBefore(string(strtrim(output))," ");
end

function pathname = canonicalPath(pathname)
[status,output] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0, error("WaveVortexBenchmark:Issue513PrimitivePath","Unable to resolve %s.",pathname); end
pathname = string(strtrim(output));
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function closeResources(reader,wvt)
if ~isempty(reader) && isvalid(reader) && ~isempty(reader.id), reader.close(); end
if ~isempty(wvt) && isvalid(wvt), delete(wvt); end
end
