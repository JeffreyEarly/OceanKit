function capabilities = wvCompiledBackendBuild(overrides)
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all");
originalDirectory = string(pwd);
directoryCleanup = onCleanup(@()cd(originalDirectory));
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
originalRng = rng;
rngCleanup = onCleanup(@()rng(originalRng));
constants = wvCompiledBackendConstants(overrides);
support = wvCompiledBackendSupport(overrides);
attempt = newAttempt(support.threadCount);
if ~support.isSupported
    capabilities = wvCompiledBackendCapabilities(overrides);
    capabilities.buildAttempt = attempt;
    clear rngCleanup pathCleanup directoryCleanup warningCleanup
    return
end

stage = "compiler";
cacheCreated = false;
try
    rejectLoadedModule(constants,overrides);
    injectFailure(overrides,stage);
    compiler = compilerRecord(overrides);
    ensureDirectory(constants.cacheRoot); cacheCreated = true;
    ensureDirectory(fullfile(constants.cacheRoot,"downloads"));
    ensureDirectory(fullfile(constants.cacheRoot,"source"));
    ensureDirectory(constants.buildDirectory);
    ensureDirectory(constants.installDirectory);
    ensureDirectory(constants.stageDirectory);
    ensureDirectory(constants.stateDirectory);
    attempt.compiler = compiler;
    writeJSON(constants.attemptPath,attempt);

    stage = "download";
    injectFailure(overrides,stage);
    downloadArchive(constants,overrides);

    stage = "checksum";
    injectFailure(overrides,stage);
    actualArchiveSHA256 = sha256File(constants.archivePath);
    if actualArchiveSHA256 ~= constants.sourceSHA256
        error("WaveVortexModel:CompiledBackendChecksumMismatch","FFTW archive SHA-256 was %s; expected %s.",actualArchiveSHA256,constants.sourceSHA256);
    end

    stage = "extract";
    injectFailure(overrides,stage);
    extractSource(constants,overrides);

    stage = "configure";
    injectFailure(overrides,stage);
    configureProvider(constants,compiler,overrides);

    stage = "provider-build";
    injectFailure(overrides,stage);
    buildProvider(constants,compiler,support.threadCount,overrides);

    stage = "provider-validation";
    injectFailure(overrides,stage);
    libraries = validateProviderLibraries(constants,overrides);

    stage = "mex-build";
    injectFailure(overrides,stage);
    stagedModule = buildStagedModule(constants,libraries,overrides);

    stage = "stage-validation";
    injectFailure(overrides,stage);
    [module,validation,storageEstimates] = validateModule(constants,stagedModule,libraries,support.threadCount);

    stage = "install";
    injectFailure(overrides,stage);
    transactionalInstall(constants,stagedModule,libraries,support.threadCount,overrides);

    stage = "record";
    record = struct( ...
        "schemaVersion","1.0.0", ...
        "provider",providerRecord(constants), ...
        "compiler",compiler, ...
        "module",module, ...
        "libraries",libraries, ...
        "validation",validation, ...
        "storageEstimates",storageEstimates, ...
        "threadCount",support.threadCount, ...
        "completedAtUTC",utcTimestamp);
    writeJSON(constants.recordPath,record);
    attempt.status = "complete";
    attempt.stage = "complete";
    attempt.completedAtUTC = utcTimestamp;
    attempt.failure = emptyFailure;
    writeJSON(constants.attemptPath,attempt);
    capabilities = wvCompiledBackendCapabilities(overrides);
catch exception
    attempt.status = "failed";
    attempt.stage = stage;
    attempt.completedAtUTC = utcTimestamp;
    attempt.failure = wvCompiledBackendFailure(stage,exception);
    if cacheCreated
        try
            writeJSON(constants.attemptPath,attempt);
        catch
        end
    end
    capabilities = wvCompiledBackendCapabilities(overrides);
    capabilities.status = "build-failed";
    capabilities.isAvailable = false;
    capabilities.buildAttempt = attempt;
    capabilities.failure = attempt.failure;
end
clear rngCleanup pathCleanup directoryCleanup warningCleanup
end

function attempt = newAttempt(threadCount)
attempt = struct("status","running","stage","preflight","startedAtUTC",utcTimestamp,"completedAtUTC","","threadCount",double(threadCount),"compiler",struct(),"failure",emptyFailure);
end

function failure = emptyFailure
failure = struct("stage","","identifier","","message","","stack",strings(0,1));
end

function provider = providerRecord(constants)
provider = struct("id",constants.providerId,"version",constants.providerVersion,"sourceURL",constants.sourceURL,"sourceURLs",constants.sourceURLs,"sourceSHA256",constants.sourceSHA256,"threadBackend",constants.threadBackend,"configureFlags",constants.configureFlags,"compilerFlags",constants.compilerFlags,"deploymentTarget",constants.deploymentTarget,"simd","NEON","openmp",false);
end

function rejectLoadedModule(constants,overrides)
loaded = moduleLoaded(constants.moduleName);
if isfield(overrides,"ModuleLoaded"), loaded = logical(overrides.ModuleLoaded); end
if loaded || mislocked(char(constants.moduleName))
    error("WaveVortexModel:CompiledBackendModuleLoaded","Clear all compiled-backend kernels and the %s MEX module before replacing it.",constants.moduleName);
end
end

function compiler = compilerRecord(overrides)
if isfield(overrides,"CompilerAvailable") && ~overrides.CompilerAvailable
    error("WaveVortexModel:CompiledBackendCompilerUnavailable","Apple Clang is unavailable.");
end
if isfield(overrides,"CompilerRecord")
    compiler = overrides.CompilerRecord;
    return
end
[status,clang] = system("/usr/bin/xcrun --find clang");
if status ~= 0, error("WaveVortexModel:CompiledBackendCompilerUnavailable","Apple Clang is unavailable."); end
[status,clangxx] = system("/usr/bin/xcrun --find clang++");
if status ~= 0, error("WaveVortexModel:CompiledBackendCompilerUnavailable","Apple Clang++ is unavailable."); end
[status,identity] = system(strtrim(string(clang))+" --version");
if status ~= 0 || ~contains(lower(string(identity)),"apple clang")
    error("WaveVortexModel:CompiledBackendCompilerIdentity","The native provider requires Apple Clang.");
end
[status,sdkPath] = system("/usr/bin/xcrun --show-sdk-path");
if status ~= 0 || ~isfolder(strtrim(string(sdkPath)))
    error("WaveVortexModel:CompiledBackendSDKUnavailable","The macOS SDK is unavailable.");
end
mexCompiler = mex.getCompilerConfigurations("C++","Selected");
if isempty(mexCompiler)
    error("WaveVortexModel:CompiledBackendMexCompilerUnavailable","No selected MATLAB C++ MEX compiler is available.");
end
compiler = struct("c",strtrim(string(clang)),"cxx",strtrim(string(clangxx)),"identity",strtrim(string(identity)),"sdkPath",strtrim(string(sdkPath)),"mexName",string(mexCompiler.Name),"mexManufacturer",string(mexCompiler.Manufacturer),"mexVersion",string(mexCompiler.Version));
end

function downloadArchive(constants,overrides)
if isfile(constants.archivePath), return, end
temporary = constants.archivePath+".partial";
deleteIfPresent(temporary);
cleanup = onCleanup(@()deleteIfPresent(temporary));
lastFailure = [];
for iURL = 1:numel(constants.sourceURLs)
    url = constants.sourceURLs(iURL);
    try
        if isfield(overrides,"DownloadFunction")
            overrides.DownloadFunction(url,temporary);
        else
            command = "/usr/bin/curl --fail --location --retry 2 --connect-timeout 30 --output "+shellQuote(temporary)+" "+shellQuote(url);
            runCommand(command,fullfile(constants.cacheRoot,"download-"+iURL+".log"),"download the pinned FFTW archive",overrides);
        end
        if isfile(temporary)
            actualSHA256 = sha256File(temporary);
            if actualSHA256 ~= constants.sourceSHA256
                error("WaveVortexModel:CompiledBackendChecksumMismatch","FFTW archive SHA-256 was %s; expected %s.",actualSHA256,constants.sourceSHA256);
            end
            break
        end
    catch exception
        lastFailure = exception;
        deleteIfPresent(temporary);
    end
end
if ~isfile(temporary) && ~isempty(lastFailure)
    rethrow(lastFailure)
end
if ~isfile(temporary)
    error("WaveVortexModel:CompiledBackendDownloadFailed","The FFTW download did not create %s.",temporary);
end
movefile(temporary,constants.archivePath,"f");
clear cleanup
end

function extractSource(constants,overrides)
if isfolder(constants.sourceDirectory), return, end
sourceParent = string(fileparts(constants.sourceDirectory));
ensureDirectory(sourceParent);
runCommand("/usr/bin/tar -xzf "+shellQuote(constants.archivePath)+" -C "+shellQuote(sourceParent),fullfile(constants.cacheRoot,"extract.log"),"extract FFTW",overrides);
if ~isfolder(constants.sourceDirectory)
    error("WaveVortexModel:CompiledBackendExtractFailed","The FFTW source directory was not created.");
end
end

function configureProvider(constants,compiler,overrides)
makefile = fullfile(constants.buildDirectory,"Makefile");
if isfile(makefile) && contains(string(fileread(makefile)),compiler.sdkPath) && contains(string(fileread(makefile)),"-mmacosx-version-min="+constants.deploymentTarget), return, end
if isfile(makefile)
    runCommand("cd "+shellQuote(constants.buildDirectory)+" && /usr/bin/make distclean",fullfile(constants.cacheRoot,"distclean.log"),"clean the incomplete FFTW build",overrides);
end
compileFlags = constants.compilerFlags+" -mmacosx-version-min="+constants.deploymentTarget+" -isysroot "+compiler.sdkPath;
linkFlags = "-mmacosx-version-min="+constants.deploymentTarget+" -isysroot "+compiler.sdkPath+" -Wl,-headerpad_max_install_names";
environment = "env SDKROOT="+shellQuote(compiler.sdkPath)+" MACOSX_DEPLOYMENT_TARGET="+constants.deploymentTarget+" CC="+shellQuote(compiler.c)+" CFLAGS="+shellQuote(compileFlags)+" LDFLAGS="+shellQuote(linkFlags);
command = "cd "+shellQuote(constants.buildDirectory)+" && "+environment+" "+shellQuote(fullfile(constants.sourceDirectory,"configure"))+" --prefix="+shellQuote(constants.installDirectory)+" "+constants.configureFlags;
runCommand(command,fullfile(constants.cacheRoot,"configure.log"),"configure FFTW",overrides);
configuration = string(fileread(fullfile(constants.cacheRoot,"configure.log")));
if ~contains(configuration,"checking whether a cycle counter is available... yes")
    error("WaveVortexModel:CompiledBackendCycleCounter","FFTW did not validate the Apple ARM cycle counter.");
end
end

function buildProvider(constants,compiler,threadCount,overrides)
baseLibrary = fullfile(constants.installDirectory,"lib","libfftw3.3.dylib");
threadLibrary = fullfile(constants.installDirectory,"lib","libfftw3_threads.3.dylib");
stampPath = fullfile(constants.installDirectory,"wv-provider-build-contract.txt");
stampValue = strjoin([constants.sourceSHA256 constants.compilerFlags constants.deploymentTarget compiler.identity],"|");
if isfile(baseLibrary) && isfile(threadLibrary) && isfile(stampPath) && string(fileread(stampPath)) == stampValue, return, end
if isfolder(constants.installDirectory), rmdir(constants.installDirectory,"s"); end
ensureDirectory(constants.installDirectory);
makeEnvironment = "env SDKROOT="+shellQuote(compiler.sdkPath)+" MACOSX_DEPLOYMENT_TARGET="+constants.deploymentTarget;
runCommand("cd "+shellQuote(constants.buildDirectory)+" && "+makeEnvironment+" /usr/bin/make -j"+threadCount,fullfile(constants.cacheRoot,"make.log"),"compile FFTW",overrides);
runCommand("cd "+shellQuote(constants.buildDirectory)+" && "+makeEnvironment+" /usr/bin/make install",fullfile(constants.cacheRoot,"install.log"),"install FFTW",overrides);
if ~isfile(baseLibrary) || ~isfile(threadLibrary)
    error("WaveVortexModel:CompiledBackendLibraryMissing","The FFTW build did not produce both shared pthread libraries.");
end
writeText(stampPath,stampValue);
end

function libraries = validateProviderLibraries(constants,overrides)
basePath = canonicalPath(fullfile(constants.installDirectory,"lib","libfftw3.3.dylib"));
threadPath = canonicalPath(fullfile(constants.installDirectory,"lib","libfftw3_threads.3.dylib"));
if ~isfile(basePath) || ~isfile(threadPath)
    error("WaveVortexModel:CompiledBackendLibraryMissing","The selected native provider libraries are missing.");
end
baseSymbols = runCommandForOutput("/usr/bin/nm -gU "+shellQuote(basePath),"inspect FFTW base symbols",overrides);
threadSymbols = runCommandForOutput("/usr/bin/nm -gU "+shellQuote(threadPath),"inspect FFTW thread symbols",overrides);
if ~contains(baseSymbols,"_fftw_execute") || ~contains(threadSymbols,"_fftw_init_threads")
    error("WaveVortexModel:CompiledBackendSymbolMissing","The native provider is missing required FFTW base or pthread symbols.");
end
dependencies = lower(runCommandForOutput("/usr/bin/otool -L "+shellQuote(basePath)+" && /usr/bin/otool -L "+shellQuote(threadPath),"inspect FFTW dependencies",overrides));
if containsOpenMP(dependencies)
    error("WaveVortexModel:CompiledBackendOpenMPRejected","The native provider links an OpenMP runtime.");
end
libraries = struct( ...
    "base",struct("path",basePath,"version",constants.providerVersion,"sha256",sha256File(basePath)), ...
    "thread",struct("path",threadPath,"version",constants.providerVersion,"sha256",sha256File(threadPath)), ...
    "openmp",struct("path","","version","","sha256","","detected",false));
end

function stagedModule = buildStagedModule(constants,libraries,overrides)
stagedModule = string(fullfile(constants.stageDirectory,constants.moduleName+"."+mexext));
deleteIfPresent(stagedModule);
if isfield(overrides,"MexBuildFunction")
    overrides.MexBuildFunction(constants,libraries,stagedModule);
    return
end
gateway = fullfile(constants.adapterDirectory,"wv_compiled_backend_mex.cpp");
engine = fullfile(constants.adapterDirectory,"WVNativeFFTWEngine.cpp");
runtimeSources = fullfile(constants.runtimeSourceDirectory,[ ...
    "WVExtensionCatalog.cpp"; ...
    "WVPortableImplementationContract.cpp"; ...
    "WVForcingContracts.cpp"; ...
    "WVForcingEngine.cpp"; ...
    "WVBarotropicQGForcingEngine.cpp"; ...
    "WVBarotropicQGIntegrationSystem.cpp"; ...
    "WVIntegrationState.cpp"; ...
    "WVRungeKutta.cpp"; ...
    "WVObserverAdapter.cpp"; ...
    "WVObserverContracts.cpp"; ...
    "WVObservation.cpp"; ...
    "WVPortableTypedRecord.cpp"; ...
    "WVOutputSchedule.cpp"; ...
    "WVFieldEvaluationService.cpp"; ...
    "WVBarotropicQGFieldEvaluationAdapter.cpp"; ...
    "WVConstantStratificationIntegrationSystem.cpp"; ...
    "WVModelTransformAdapters.cpp"; ...
    "WVModel.cpp"]);
runtimeHeaders = fullfile(constants.runtimeIncludeDirectory,"WaveVortexRuntime",[ ...
    "WVBarotropicQGForcingEngine.hpp"; ...
    "WVCheckpointReader.hpp"; ...
    "WVExtensionCatalog.hpp"; ...
    "WVFieldEvaluationService.hpp"; ...
    "WVForcing.hpp"; ...
    "WVForcingContracts.hpp"; ...
    "WVModel.hpp"; ...
    "WVModelOutputConfiguration.hpp"; ...
    "WVObservation.hpp"; ...
    "WVObserverContracts.hpp"; ...
    "WVObserverOutputProvider.hpp"; ...
    "WVObservingSystem.hpp"; ...
    "WVOutputSchedule.hpp"; ...
    "WVPortableImplementationContract.hpp"; ...
    "WVPortableTypedRecord.hpp"]);
runtimePrivateHeaders = fullfile(constants.runtimeSourceDirectory,[ ...
    "WVBarotropicQGFieldEvaluationAdapter.hpp"; ...
    "WVLegacyObserverCompatibility.hpp"; ...
    "WVModelTransformAdapters.hpp"; ...
    "WVModelInternalAccess.hpp"]);
coreSources = fullfile(constants.coreSourceDirectory,[ ...
    "WVKernelTypes.cpp"; ...
    "WVTransformBarotropicQGKernel.cpp"; ...
    "WVTransformConstantStratificationKernel.cpp"]);
requiredSources = [string(gateway);string(engine);fullfile(constants.adapterDirectory,"WVNativeFFTWEngine.hpp");coreSources;runtimeHeaders;runtimePrivateHeaders;runtimeSources];
if any(~isfile(requiredSources))
    error("WaveVortexModel:CompiledBackendSourceMissing","A tracked compiled-backend source file is missing: %s",strjoin(requiredSources(~isfile(requiredSources)),", "));
end
compilerFlags = "CXXFLAGS=$CXXFLAGS -std=c++17 -pthread -O3 -mcpu=native -mmacosx-version-min="+constants.deploymentTarget+" -DWV_KERNEL_NATIVE_OPTIMIZATION=1 -DWV_KERNEL_COEFFICIENT_WORKERS=2 -DWV_MODEL_ENABLE_OUTPUT=0";
linkerFlags = "LDFLAGS=$LDFLAGS -pthread -mmacosx-version-min="+constants.deploymentTarget+" -Wl,-rpath,"+fileparts(libraries.base.path);
coreArguments = cellstr(coreSources);
runtimeArguments = cellstr(runtimeSources);
mex("-R2018a",compilerFlags,gateway,coreArguments{:},runtimeArguments{:},engine,"-I"+constants.coreIncludeDirectory,"-I"+constants.runtimeIncludeDirectory,"-I"+constants.runtimeSourceDirectory,"-I"+constants.adapterDirectory,"-I"+fullfile(constants.installDirectory,"include"),linkerFlags,libraries.thread.path,libraries.base.path,"-outdir",constants.stageDirectory,"-output",constants.moduleName);
if ~isfile(stagedModule)
    error("WaveVortexModel:CompiledBackendMexMissing","The MEX build did not create the staged module.");
end
end

function [module,validation,storageEstimates] = validateModule(constants,modulePath,libraries,threadCount)
moduleDirectory = string(fileparts(modulePath));
originalDirectory = string(pwd);
directoryCleanup = onCleanup(@()cd(originalDirectory));
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(moduleDirectory,"-begin");
cd(moduleDirectory);
if moduleLoaded(constants.moduleName)
    error("WaveVortexModel:CompiledBackendModuleLoaded","The staged module name is already loaded.");
end
resolved = string(which(constants.moduleName));
if ~samePath(resolved,modulePath)
    error("WaveVortexModel:CompiledBackendStageIdentity","MATLAB resolved %s instead of the staged module %s.",resolved,modulePath);
end
info = feval(char(constants.moduleName),'moduleInfo');
if ~samePath(string(info.baseLibrary),libraries.base.path) || ~samePath(string(info.threadLibrary),libraries.thread.path)
    error("WaveVortexModel:CompiledBackendIdentityMismatch","dladdr did not resolve the exact selected FFTW base and thread libraries.");
end
if startsWith(string(info.baseLibrary),string(matlabroot)) || startsWith(string(info.threadLibrary),string(matlabroot))
    error("WaveVortexModel:CompiledBackendBundledFFTWRejected","The module resolved MATLAB-bundled FFTW.");
end
if string(info.openMPRuntimeLibrary) ~= ""
    error("WaveVortexModel:CompiledBackendOpenMPRejected","The module resolved an OpenMP runtime.");
end
dependencies = lower(runCommandForOutput("/usr/bin/otool -L "+shellQuote(modulePath),"inspect MEX dependencies",struct()));
if containsOpenMP(dependencies)
    error("WaveVortexModel:CompiledBackendOpenMPRejected","The MEX module links an OpenMP runtime.");
end
[validation,storageEstimates] = runSelfTests(constants,threadCount);
module = struct("name",constants.moduleName,"path",constants.installedModule,"stagedPath",string(modulePath),"sha256",sha256File(modulePath),"engine",string(info.engine));
clearModule(constants.moduleName);
if moduleLoaded(constants.moduleName)
    error("WaveVortexModel:CompiledBackendModuleCleanup","The staged MEX module remained loaded after validation.");
end
clear pathCleanup directoryCleanup
end

function [validation,storageEstimates] = runSelfTests(constants,threadCount)
definitions = struct("id",{"hydrostatic","nonhydrostatic"},"isHydrostatic",{true,false},"seed",{16901,16902});
caseResults = repmat(struct("status","failed","isHydrostatic",false,"maximumRelativeError",NaN,"errors",struct(),"planCount",NaN,"contractVersion",NaN,"lifecyclePassed",false,"threadCount",threadCount),2,1);
estimateCases = repmat(struct("id","","Nxyz",[],"isHydrostatic",false,"persistentBytes",NaN,"stateInputBytes",NaN,"fluxOutputBytes",NaN,"knownMaximumLiveOwnedBytes",NaN,"planBytes",NaN,"planMemoryAccounting",""),2,1);
for iCase = 1:2
    definition = definitions(iCase);
    rng(definition.seed,"twister");
    wvt = WVTransformConstantStratification([15000 12000 1300],[8 6 7],N0=5.2e-3,latitude=33,isHydrostatic=definition.isHydrostatic,shouldAntialias=true);
    wvt.initWithRandomFlow(uvMax=0.01);
    estimate = feval(char(constants.moduleName),'estimate',kernelConfiguration(wvt));
    before = feval(char(constants.moduleName),'moduleMetrics');
    handle = feval(char(constants.moduleName),'create',kernelConfiguration(wvt),threadCount);
    handleCleanup = onCleanup(@()deleteKernel(constants.moduleName,handle));
    [U,V,W,N] = wvt.transformWaveVortexToUVWEta(wvt.Ap,wvt.Am,wvt.A0,wvt.t);
    if definition.isHydrostatic, fields = cat(4,U,V,N); else, fields = cat(4,U,V,W,N); end
    [actualAp,actualAm,actualA0] = feval(char(constants.moduleName),'forward',handle,fields,wvt.t,wvt.t0);
    if definition.isHydrostatic
        [expectedAp,expectedAm,expectedA0] = wvt.transformUVEtaToWaveVortex(U,V,N);
    else
        [expectedAp,expectedAm,expectedA0] = wvt.transformUVWEtaToWaveVortex(U,V,W,N);
    end
    actualInverse = feval(char(constants.moduleName),'inverse',handle,wvt.Ap,wvt.Am,wvt.A0,wvt.t,wvt.t0);
    [actualFp,actualFm,actualF0] = feval(char(constants.moduleName),'nonlinearFlux',handle,wvt.Ap,wvt.Am,wvt.A0,wvt.t,wvt.t0);
    [expectedFp,expectedFm,expectedF0] = wvt.nonlinearFlux();
    errors = struct( ...
        "forward",max([relativeError(actualAp,expectedAp) relativeError(actualAm,expectedAm) relativeError(actualA0,expectedA0)]), ...
        "inverse",relativeError(actualInverse,cat(4,U,V,W,N)), ...
        "nonlinearFlux",max([relativeError(actualFp,expectedFp) relativeError(actualFm,expectedFm) relativeError(actualF0,expectedF0)]));
    metrics = feval(char(constants.moduleName),'metrics',handle);
    feval(char(constants.moduleName),'delete',handle); handle = []; clear handleCleanup
    after = feval(char(constants.moduleName),'moduleMetrics');
    lifecyclePassed = after.kernelCount == before.kernelCount && after.activePlans == before.activePlans && after.outstandingPlanningBytes == 0 && after.totalPlansCreated-after.totalPlansDestroyed == after.activePlans;
    maximumRelativeError = max(cell2mat(struct2cell(errors)));
    estimatePassed = estimate.descriptorBytes == metrics.descriptorBytes && estimate.scratchCapacityBytes == metrics.scratchCapacityBytes && estimate.stateInputBytes == metrics.stateInputBytes && estimate.fluxOutputBytes == metrics.fluxOutputBytes && estimate.knownMaximumLiveOwnedBytesLowerBound+metrics.planBytes+metrics.engineBytes+metrics.kernelManagementBytes == metrics.knownMaximumLiveOwnedBytes && metrics.persistentBytes == metrics.descriptorBytes+metrics.scratchCapacityBytes+metrics.planBytes+metrics.engineBytes+metrics.kernelManagementBytes;
    if maximumRelativeError > 1e-12 || metrics.planCount ~= 17 || metrics.contractVersion ~= constants.contractVersion || ~lifecyclePassed || ~estimatePassed
        error("WaveVortexModel:CompiledBackendSelfTest","The %s self-test failed: error %.3g, plans %d, contract %d, lifecycle %d.",definition.id,maximumRelativeError,metrics.planCount,metrics.contractVersion,lifecyclePassed);
    end
    caseResults(iCase) = struct("status","passed","isHydrostatic",definition.isHydrostatic,"maximumRelativeError",maximumRelativeError,"errors",errors,"planCount",metrics.planCount,"contractVersion",metrics.contractVersion,"lifecyclePassed",lifecyclePassed,"threadCount",threadCount);
    estimateCases(iCase) = struct("id",definition.id,"Nxyz",[estimate.Nx estimate.Ny estimate.Nz],"isHydrostatic",definition.isHydrostatic,"persistentBytes",metrics.persistentBytes,"stateInputBytes",metrics.stateInputBytes,"fluxOutputBytes",metrics.fluxOutputBytes,"knownMaximumLiveOwnedBytes",metrics.knownMaximumLiveOwnedBytes,"planBytes",metrics.planBytes,"planMemoryAccounting",string(metrics.planMemoryAccounting));
end
validation = struct("status","passed","maximumRelativeError",max([caseResults.maximumRelativeError]),"hydrostatic",caseResults(1),"nonhydrostatic",caseResults(2),"deterministic",true,"tolerance",1e-12);
storageEstimates = struct("status","measured-after-kernel-construction","formulaSource","Application-owned kernel, provider, descriptor, plan-wrapper, and bounded scratch storage; FFTW plan internals remain opaque","cases",estimateCases);
end

function transactionalInstall(constants,stagedModule,libraries,threadCount,overrides)
rejectLoadedModule(constants,overrides);
backup = fullfile(constants.stageDirectory,"previous-"+constants.moduleName+"."+mexext);
deleteIfPresent(backup);
hadPrevious = isfile(constants.installedModule);
if hadPrevious, movefile(constants.installedModule,backup,"f"); end
installed = false;
try
    copyfile(stagedModule,constants.installedModule,"f"); installed = true;
    injectFailure(overrides,"install-validation");
    [~,validation,~] = validateModule(constants,constants.installedModule,libraries,threadCount);
    if string(validation.status) ~= "passed"
        error("WaveVortexModel:CompiledBackendInstallValidation","The installed module did not pass validation.");
    end
    deleteIfPresent(backup);
catch exception
    if installed, deleteIfPresent(constants.installedModule); end
    if hadPrevious && isfile(backup), movefile(backup,constants.installedModule,"f"); end
    rethrow(exception)
end
end

function runCommand(command,logPath,description,overrides)
ensureDirectory(fileparts(logPath));
if isfield(overrides,"CommandFunction")
    [status,output] = overrides.CommandFunction(command,logPath);
else
    [status,output] = system(command+" > "+shellQuote(logPath)+" 2>&1");
end
if status ~= 0
    details = string(output);
    if isfile(logPath), details = string(fileread(logPath)); end
    error("WaveVortexModel:CompiledBackendCommandFailed","Unable to %s.%s%s",description,newline,details);
end
end

function output = runCommandForOutput(command,description,overrides)
if isfield(overrides,"CommandOutputFunction")
    [status,output] = overrides.CommandOutputFunction(command);
else
    [status,output] = system(command);
end
if status ~= 0
    error("WaveVortexModel:CompiledBackendCommandFailed","Unable to %s.%s%s",description,newline,output);
end
output = string(output);
end

function injectFailure(overrides,stage)
if isfield(overrides,"FailureStage") && string(overrides.FailureStage) == string(stage)
    error("WaveVortexModel:CompiledBackendInjectedFailure","Injected %s failure.",stage);
end
end

function configuration = kernelConfiguration(wvt)
configuration = struct("Nx",wvt.Nx,"Ny",wvt.Ny,"Nz",wvt.Nz,"Nj",wvt.Nj,"Lx",wvt.Lx,"Ly",wvt.Ly,"Lz",wvt.Lz,"N0",wvt.N0,"rho0",wvt.rho0,"g",wvt.g,"planetaryRadius",wvt.planetaryRadius,"rotationRate",wvt.rotationRate,"latitude",wvt.latitude,"isHydrostatic",wvt.isHydrostatic,"shouldAntialias",wvt.shouldAntialias);
end

function value = relativeError(actual,expected)
value = max(abs(actual(:)-expected(:)),[],"omitmissing")/max(max(abs(expected(:)),[],"omitmissing"),realmin);
end

function deleteKernel(moduleName,handle)
if isempty(handle), return, end
try
    feval(char(moduleName),'delete',handle);
catch
end
end

function tf = containsOpenMP(text)
tf = contains(lower(string(text)),["libomp" "libgomp" "libiomp" "openmp"]);
tf = any(tf);
end

function tf = moduleLoaded(moduleName)
[~,mexFiles] = inmem("-completenames");
tf = any(endsWith(string(mexFiles),filesep+moduleName+"."+mexext));
end

function clearModule(moduleName)
eval("clear "+moduleName);
end

function ensureDirectory(pathname)
if ~isfolder(pathname), mkdir(pathname); end
end

function deleteIfPresent(pathname)
if isfile(pathname), delete(pathname); end
end

function hash = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status ~= 0, error("WaveVortexModel:CompiledBackendHash","Unable to hash %s.",pathname); end
hash = extractBefore(string(strtrim(output))," ");
end

function pathname = canonicalPath(pathname)
[status,output] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0, error("WaveVortexModel:CompiledBackendPath","Unable to resolve %s.",pathname); end
pathname = string(strtrim(output));
end

function tf = samePath(first,second)
if first == "" || second == "", tf = false; return, end
tf = canonicalPath(first) == canonicalPath(second);
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function writeJSON(pathname,value)
ensureDirectory(fileparts(pathname));
temporary = string(pathname)+".tmp";
fileId = fopen(temporary,"w");
if fileId < 0, error("WaveVortexModel:CompiledBackendStateWrite","Unable to write %s.",temporary); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",jsonencode(value,PrettyPrint=true));
clear cleanup
movefile(temporary,pathname,"f");
end

function writeText(pathname,value)
fileId = fopen(pathname,"w");
if fileId < 0, error("WaveVortexModel:CompiledBackendStateWrite","Unable to write %s.",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",value);
clear cleanup
end

function value = utcTimestamp
value = string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end
