function capabilities = wvCompiledBackendCapabilities(overrides)
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all");
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
emptyFailure = struct("stage","","identifier","","message","","stack",strings(0,1));
try
    constants = wvCompiledBackendConstants(overrides);
    support = wvCompiledBackendSupport(overrides);
    if ~pathContains(constants.packageRoot), addpath(constants.packageRoot); end
    attempt = readRecord(constants.attemptPath,emptyAttempt);
    validatedBuild = readRecord(constants.recordPath,struct());
    sourceIdentity = emptySourceIdentity;
    if isfield(validatedBuild,"sourceIdentity"), sourceIdentity = validatedBuild.sourceIdentity; end
    provider = struct( ...
        "id",constants.providerId, ...
        "version",constants.providerVersion, ...
        "sourceURL",constants.sourceURL, ...
        "sourceURLs",constants.sourceURLs, ...
        "sourceSHA256",constants.sourceSHA256, ...
        "threadBackend",constants.threadBackend, ...
        "configureFlags",constants.configureFlags, ...
        "compilerFlags",constants.compilerFlags, ...
        "deploymentTarget",constants.deploymentTarget, ...
        "simd","NEON", ...
        "sharedLibraries",true, ...
        "fortran",false, ...
        "openmp",false, ...
        "kernelScope","constant hydrostatic/nonhydrostatic, Hydrostatic, Boussinesq, Stratified QG and Barotropic QG");
    module = emptyModule(constants);
    libraries = emptyLibraries;
    validation = emptyValidation;
    storageEstimates = emptyStorageEstimates;
    compiler = struct();
    if ~isempty(fieldnames(validatedBuild))
        if isfield(validatedBuild,"compiler"), compiler = validatedBuild.compiler; end
        if isfield(validatedBuild,"validation"), validation = validatedBuild.validation; end
        if isfield(validatedBuild,"storageEstimates"), storageEstimates = validatedBuild.storageEstimates; end
    end
    status = "unsupported";
    failure = emptyFailure;
    if ~support.isSupported
        failure = expectedFailure("support","WaveVortexModel:CompiledBackendUnsupported",strjoin(support.reasons," "));
    elseif ~isfile(constants.installedModule)
        status = "not-built";
        failure = expectedFailure("availability","WaveVortexModel:CompiledBackendNotBuilt","The native compiled backend has not been built. Call WVCompiledBackend.build().");
    else
        [module,libraries,inspectionFailure] = inspectModule(constants,validatedBuild);
        if inspectionFailure.identifier == ""
            status = "available";
        else
            status = "invalid";
            failure = inspectionFailure;
        end
    end
    capabilities = assembleCapabilities(constants,support,provider,compiler,module,libraries,sourceIdentity,validation,storageEstimates,attempt,status,failure);
catch exception
    fallbackConstants = fallbackConstantsRecord;
    fallbackSupport = fallbackSupportRecord;
    try
        fallbackConstants = wvCompiledBackendConstants(overrides);
    catch
    end
    try
        fallbackSupport = wvCompiledBackendSupport(overrides);
    catch
    end
    capabilities = assembleCapabilities(fallbackConstants,fallbackSupport,struct(),struct(),emptyModule(fallbackConstants),emptyLibraries,emptySourceIdentity,emptyValidation,emptyStorageEstimates,emptyAttempt,"invalid",wvCompiledBackendFailure("inspection",exception));
end
clear pathCleanup warningCleanup
end

function capabilities = assembleCapabilities(constants,support,provider,compiler,module,libraries,sourceIdentity,validation,storageEstimates,attempt,status,failure)
capabilities = struct( ...
    "schemaVersion","1.0.0", ...
    "status",string(status), ...
    "isSupported",logical(support.isSupported), ...
    "isAvailable",string(status) == "available", ...
    "matlab",support.matlab, ...
    "platform",support.platform, ...
    "provider",provider, ...
    "compiler",compiler, ...
    "module",module, ...
    "libraries",libraries, ...
    "sourceIdentity",sourceIdentity, ...
    "contract",struct("version",double(constants.contractVersion),"threadCount",double(support.threadCount),"defaultThreadRule","min(18,maxNumCompThreads)","planCount",17,"planCountMeaning","logical-prepared-operation-slots"), ...
    "featureValidation",validation, ...
    "storageEstimates",storageEstimates, ...
    "buildAttempt",attempt, ...
    "cache",struct("root",string(constants.cacheRoot),"locationPolicy",string(constants.cacheLocationPolicy),"packageKey",string(constants.packageKey),"isIgnoredLocalState",true), ...
    "failure",failure);
end

function [module,libraries,failure] = inspectModule(constants,validatedBuild)
failure = struct("stage","","identifier","","message","","stack",strings(0,1));
module = emptyModule(constants);
libraries = emptyLibraries;
wasLoaded = moduleLoaded(constants.moduleName);
module.loadedBeforeInspection = wasLoaded;
module.path = string(constants.installedModule);
module.sha256 = sha256File(constants.installedModule);
try
    wvCompiledBackendResolveModule(constants.moduleName,constants.installedModule);
    info = feval(char(constants.moduleName),'moduleInfo');
    libraries.base = identityRecord(info.baseLibrary,info.version);
    libraries.thread = identityRecord(info.threadLibrary,info.version);
    libraries.openmp = identityRecord(info.openMPRuntimeLibrary,"");
    libraries.openmp.detected = string(info.openMPRuntimeLibrary) ~= "";
    module.engine = string(info.engine);
    % Inspect the installed binary; never infer its schedule from newer source
    % defaults. Source and binary identity must match the validated build record.
    for name = ["nonlinearFluxSchedule","executionScheduleVersion","workerPolicyIdentifier","requestedHorizontalWorkers","requestedPointwiseWorkers","planCountMeaning","matlabTransformBridgeVersion"]
        if isfield(info,name), module.(name) = info.(name); end
    end
    module.identityValidated = validateBinaryIdentities(constants,validatedBuild,module,libraries);
    if ~module.identityValidated
        error("WaveVortexModel:CompiledBackendIdentityMismatch","The installed module or native provider does not match the validated build record.");
    end
    module.isInstalled = true;
catch exception
    failure = wvCompiledBackendFailure("module-inspection",exception);
end
if ~wasLoaded
    clearModule(constants.moduleName);
end
module.loadedAfterInspection = moduleLoaded(constants.moduleName);
end

function tf = validateBinaryIdentities(constants,validatedBuild,module,libraries)
tf = libraries.base.path ~= "" && libraries.thread.path ~= "" && ~libraries.openmp.detected;
tf = tf && ~startsWith(libraries.base.path,string(matlabroot)) && ~startsWith(libraries.thread.path,string(matlabroot));
if isempty(fieldnames(validatedBuild)) || ~isfield(validatedBuild,"libraries") || ...
        ~isfield(validatedBuild,"module") || ~isfield(validatedBuild,"sourceIdentity")
    tf = false;
    return
end
expected = validatedBuild.libraries;
tf = tf && samePath(libraries.base.path,string(expected.base.path)) && samePath(libraries.thread.path,string(expected.thread.path));
tf = tf && string(validatedBuild.provider.id) == constants.providerId;
tf = tf && module.sha256 == string(validatedBuild.module.sha256);
tf = tf && libraries.base.sha256 == string(expected.base.sha256) && ...
    libraries.thread.sha256 == string(expected.thread.sha256);
identity = validatedBuild.sourceIdentity;
tf = tf && string(identity.moduleSHA256) == module.sha256 && ...
    string(identity.baseLibrarySHA256) == libraries.base.sha256 && ...
    string(identity.threadLibrarySHA256) == libraries.thread.sha256 && ...
    string(identity.providerId) == constants.providerId && ...
    string(identity.providerSourceSHA256) == constants.sourceSHA256;
end

function record = identityRecord(pathname,versionValue)
record = struct("path",string(pathname),"version",string(versionValue),"sha256","");
if record.path ~= "" && isfile(record.path)
    record.sha256 = sha256File(record.path);
end
end

function module = emptyModule(constants)
module = struct("name",string(constants.moduleName),"path",string(constants.installedModule),"sha256","","engine","","isInstalled",isfile(constants.installedModule),"identityValidated",false,"loadedBeforeInspection",false,"loadedAfterInspection",false);
end

function libraries = emptyLibraries
empty = struct("path","","version","","sha256","");
openmp = empty; openmp.detected = false;
libraries = struct("base",empty,"thread",empty,"openmp",openmp);
end

function identity = emptySourceIdentity
identity = struct("schemaVersion","","manifest",struct(),"moduleSHA256","", ...
    "baseLibrarySHA256","","threadLibrarySHA256","","providerId","", ...
    "providerSourceSHA256","");
end

function validation = emptyValidation
caseRecord = struct("status","not-run","isHydrostatic",false,"maximumRelativeError",NaN,"planCount",NaN,"lifecyclePassed",false);
validation = struct("status","not-run","maximumRelativeError",NaN,"hydrostatic",caseRecord,"nonhydrostatic",caseRecord,"deterministic",false,"tolerance",1e-12);
end

function estimates = emptyStorageEstimates
estimates = struct("status","not-run","formulaSource","compiled core metrics","cases",struct([]));
end

function attempt = emptyAttempt
attempt = struct("status","not-attempted","stage","","startedAtUTC","","completedAtUTC","","threadCount",NaN,"failure",struct("stage","","identifier","","message","","stack",strings(0,1)));
end

function failure = expectedFailure(stage,identifier,message)
failure = struct("stage",string(stage),"identifier",string(identifier),"message",string(message),"stack",strings(0,1));
end

function value = readRecord(pathname,defaultValue)
value = defaultValue;
if ~isfile(pathname), return, end
value = jsondecode(fileread(pathname));
end

function tf = pathContains(folder)
tf = any(string(strsplit(path,pathsep)) == string(folder));
end

function tf = moduleLoaded(moduleName)
[~,mexFiles] = inmem("-completenames");
tf = any(endsWith(string(mexFiles),filesep+moduleName+"."+mexext));
end

function clearModule(moduleName)
eval("clear "+moduleName);
end

function hash = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status ~= 0
    error("WaveVortexModel:CompiledBackendHash","Unable to hash %s.",pathname);
end
hash = extractBefore(string(strtrim(output))," ");
end

function tf = samePath(first,second)
tf = canonicalPath(first) == canonicalPath(second);
end

function pathname = canonicalPath(pathname)
[status,output] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0, pathname = string(pathname); else, pathname = string(strtrim(output)); end
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function constants = fallbackConstantsRecord
constants = struct("contractVersion",4,"cacheRoot","","cacheLocationPolicy","","packageKey","","moduleName","wv_compiled_backend_mex","installedModule","","packageRoot","");
end

function support = fallbackSupportRecord
support = struct("matlab",struct("release","","version","","minimumRelease","R2025b","isSupported",false),"platform",struct("architecture","","operatingSystem","","isSupported",false),"isSupported",false,"reasons","Inspection failed.","threadCount",NaN);
end
