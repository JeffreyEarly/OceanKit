classdef TestWVCompiledBackend < matlab.unittest.TestCase
    properties
        repositoryRoot
    end

    methods (TestClassSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (TestMethodSetup)
        function unloadUnusedCompiledModule(testCase)
            % Other test classes may leave the owner-free MEX mapped after
            % deleting their transforms. Build tests need an unloaded module
            % so they exercise the requested build stage.
            [~,mexFiles] = inmem("-completenames");
            moduleIsLoaded = any(endsWith(string(mexFiles),filesep+"wv_compiled_backend_mex."+mexext));
            if moduleIsLoaded
                metrics = wv_compiled_backend_mex('moduleMetrics');
                testCase.assertEqual(metrics.kernelCount,0,"A prior test leaked a compiled kernel owner.");
                testCase.assertEqual(metrics.matlabTransformCount,0,"A prior test leaked a MATLAB transform owner.");
            end
            clear wv_compiled_backend_mex
        end
    end

    methods (Test, TestTags="smoke")
        function capabilitiesAreNonthrowingAndDoNotBuild(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = string(fixture.Folder);
            cacheRoot = fullfile(root,"cache");
            action = @()WVCompiledBackend.capabilitiesForTesting(struct("PackageRoot",root,"CacheRoot",cacheRoot,"Architecture","maca64","OperatingSystem","macOS","Release","R2025b","MaxThreads",16));
            capabilities = testCase.verifyWarningFree(action);
            testCase.verifyEqual(capabilities.schemaVersion,"1.0.0");
            testCase.verifyEqual(capabilities.status,"not-built");
            testCase.verifyEqual(capabilities.contract.threadCount,16);
            testCase.verifyEqual(capabilities.provider.id,"native-neon-pthreads");
            testCase.verifyEqual(capabilities.provider.sourceSHA256,"5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1");
            testCase.verifyEqual(capabilities.cache.root,string(cacheRoot));
            testCase.verifyFalse(isfolder(cacheRoot));
        end

        function defaultCacheIsExternalAndPackageSpecific(testCase)
            firstFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            secondFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            options = struct("Architecture","maca64","OperatingSystem","macOS","Release","R2025b");
            first = options; first.PackageRoot = string(firstFixture.Folder);
            second = options; second.PackageRoot = string(secondFixture.Folder);
            firstCapabilities = WVCompiledBackend.capabilitiesForTesting(first);
            withoutReleasePrefix = first; withoutReleasePrefix.Release = "2025b";
            repeatedCapabilities = WVCompiledBackend.capabilitiesForTesting(withoutReleasePrefix);
            secondCapabilities = WVCompiledBackend.capabilitiesForTesting(second);
            testCase.verifyEqual(firstCapabilities.cache.locationPolicy,"user-cache-package-path-v1");
            testCase.verifyEqual(firstCapabilities.cache.root,repeatedCapabilities.cache.root);
            testCase.verifyNotEqual(firstCapabilities.cache.root,secondCapabilities.cache.root);
            userCache = fullfile(string(java.lang.System.getProperty("user.home")),"Library","Caches","WaveVortexModel","compiled-backend-v1");
            testCase.verifyTrue(startsWith(firstCapabilities.cache.root,userCache+filesep));
            testCase.verifyFalse(startsWith(firstCapabilities.cache.root,string(firstFixture.Folder)+filesep));
        end

        function activationRejectsIncompleteOrMissingModuleIdentity(testCase)
            unavailable = struct("isAvailable",false,"module",struct());
            testCase.verifyError(@()WVCompiledBackend.activateModule(unavailable),"WaveVortexModel:CompiledBackendUnavailable");
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            missingPath = fullfile(fixture.Folder,"missing."+mexext);
            missing = struct("isAvailable",true,"module",struct("identityValidated",true,"name","wv_compiled_backend_mex","path",missingPath,"sha256",string(repmat('0',1,64))));
            testCase.verifyError(@()WVCompiledBackend.activateModule(missing),"WaveVortexModel:CompiledBackendNotBuilt");
        end

        function unsupportedSystemsProduceStructuredUnavailability(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = string(fixture.Folder);
            oldRelease = WVCompiledBackend.capabilitiesForTesting(struct("PackageRoot",root,"Architecture","maca64","OperatingSystem","macOS","Release","R2025a"));
            testCase.verifyEqual(oldRelease.status,"unsupported");
            testCase.verifyEqual(oldRelease.failure.identifier,"WaveVortexModel:CompiledBackendUnsupported");
            wrongArchitecture = WVCompiledBackend.capabilitiesForTesting(struct("PackageRoot",root,"Architecture","maci64","OperatingSystem","macOS","Release","R2025b"));
            testCase.verifyEqual(wrongArchitecture.status,"unsupported");
            testCase.verifyFalse(wrongArchitecture.platform.isSupported);
            testCase.verifyFalse(isfolder(fullfile(root,".compiled-backend-cache")));
        end

        function portableLinuxPathHasNoDownloadBuildOrOutput(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = string(fixture.Folder);
            capabilities = WVCompiledBackend.buildForTesting(struct("PackageRoot",root,"Architecture","glnxa64","OperatingSystem","Linux","Release","R2026a","DownloadFunction",@unexpectedDownload));
            testCase.verifyEqual(capabilities.status,"unsupported");
            testCase.verifyFalse(capabilities.isAvailable);
            testCase.verifyFalse(isfolder(fullfile(root,".compiled-backend-cache")));
            entries = dir(root); entries = entries(~ismember(string({entries.name}),["." ".."]));
            testCase.verifyEmpty(entries);
        end

        function capabilityAndUnsupportedBuildRestoreState(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            originalDirectory = string(pwd); originalPath = path; originalRng = rng; originalWarnings = warning;
            addpath(testCase.repositoryRoot,"-begin"); expectedPath = path;
            cd(fixture.Folder); rng(169,"twister"); expectedRng = rng;
            cleanup = onCleanup(@()restoreState(originalDirectory,originalPath,originalRng,originalWarnings));
            WVCompiledBackend.capabilitiesForTesting(struct("PackageRoot",string(fixture.Folder),"Architecture","glnxa64","OperatingSystem","Linux","Release","R2026a"));
            WVCompiledBackend.buildForTesting(struct("PackageRoot",string(fixture.Folder),"Architecture","glnxa64","OperatingSystem","Linux","Release","R2026a"));
            testCase.verifyEqual(string(pwd),string(fixture.Folder));
            testCase.verifyEqual(path,expectedPath);
            testCase.verifyEqual(rng,expectedRng);
            clear cleanup
        end
    end

    methods (Test, TestTags="full")
        function compilerDownloadChecksumAndBuildFailuresAreStructured(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = string(fixture.Folder);
            common = struct("PackageRoot",root,"CacheRoot",fullfile(root,"cache"),"Architecture","maca64","OperatingSystem","macOS","Release","R2025b","MaxThreads",2,"CompilerRecord",fakeCompilerRecord(root));
            compilerOptions = common; compilerOptions.CompilerAvailable = false;
            compilerFailure = WVCompiledBackend.buildForTesting(compilerOptions);
            testCase.verifyEqual(compilerFailure.status,"build-failed");
            testCase.verifyEqual(compilerFailure.failure.identifier,"WaveVortexModel:CompiledBackendCompilerUnavailable");
            testCase.verifyFalse(isfolder(fullfile(root,".compiled-backend-cache")));

            downloadOptions = common; downloadOptions.DownloadFunction = @failedDownload;
            downloadFailure = WVCompiledBackend.buildForTesting(downloadOptions);
            testCase.verifyEqual(downloadFailure.buildAttempt.stage,"download");
            testCase.verifyEqual(downloadFailure.failure.identifier,"WaveVortexModel:TestDownloadFailure");

            archiveFolder = fullfile(common.CacheRoot,"downloads");
            if ~isfolder(archiveFolder), mkdir(archiveFolder); end
            writelines("not fftw",fullfile(archiveFolder,"fftw-3.3.11.tar.gz"));
            checksumFailure = WVCompiledBackend.buildForTesting(common);
            testCase.verifyEqual(checksumFailure.buildAttempt.stage,"checksum");
            testCase.verifyEqual(checksumFailure.failure.identifier,"WaveVortexModel:CompiledBackendChecksumMismatch");

            delete(fullfile(archiveFolder,"fftw-3.3.11.tar.gz"));
            injected = common; injected.FailureStage = "download";
            buildFailure = WVCompiledBackend.buildForTesting(injected);
            testCase.verifyEqual(buildFailure.status,"build-failed");
            testCase.verifyEqual(buildFailure.failure.identifier,"WaveVortexModel:CompiledBackendInjectedFailure");
        end

        function loadedModuleReplacementIsRefused(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            options = struct("PackageRoot",string(fixture.Folder),"Architecture","maca64","OperatingSystem","macOS","Release","R2025b","ModuleLoaded",true);
            capabilities = WVCompiledBackend.buildForTesting(options);
            testCase.verifyEqual(capabilities.status,"build-failed");
            testCase.verifyEqual(capabilities.failure.identifier,"WaveVortexModel:CompiledBackendModuleLoaded");
            testCase.verifyFalse(isfolder(fullfile(fixture.Folder,".compiled-backend-cache")));
        end

        function unrelatedCwdAndTrackedExportResolveProductionSources(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            exportRoot = fullfile(fixture.Folder,"export"); mkdir(exportRoot);
            tracked = splitlines(strtrim(string(runGit(testCase.repositoryRoot,"ls-files"))));
            classFiles = tracked(startsWith(tracked,"@WVCompiledBackend/"));
            buildSource = string(fileread(fullfile(testCase.repositoryRoot,"@WVCompiledBackend","private","wvCompiledBackendBuild.m")));
            sourceTokens = regexp(buildSource,'"WV[^"/]+\.cpp"','match');
            sourceNames = unique(erase(string(sourceTokens(:)),'"'));
            coreFiles = tracked(startsWith(tracked,"CompiledKernel/src/") & endsWith(tracked,".cpp"));
            [~,coreNames] = fileparts(coreFiles);
            % The MEX gateway enters WVModel factories, so the complete kernel
            % library is a link dependency even for constant-model commands.
            testCase.verifyTrue(all(ismember(coreNames+".cpp",sourceNames)),"The production MEX source list omits a kernel translation unit.");
            runtimeFiles = "PortableRuntime/src/"+sourceNames;
            runtimeFiles = runtimeFiles(ismember(runtimeFiles,tracked));
            headers = tracked((startsWith(tracked,"CompiledKernel/") | startsWith(tracked,"PortableRuntime/")) & endsWith(tracked,".hpp"));
            required = [classFiles;"CompiledKernel/native-fftw-provider.env";"CompiledKernel/adapters/native-fftw/WVNativeFFTWEngine.cpp";"CompiledKernel/adapters/native-fftw/wv_compiled_backend_mex.cpp";coreFiles;runtimeFiles;headers];
            testCase.verifyTrue(all(ismember(required,tracked)));
            for relative = required'
                destination = fullfile(exportRoot,relative);
                if ~isfolder(fileparts(destination)), mkdir(fileparts(destination)); end
                copyfile(fullfile(testCase.repositoryRoot,relative),destination);
            end
            originalDirectory = pwd; originalPath = path; cleanup = onCleanup(@()restoreDirectoryAndPath(originalDirectory,originalPath));
            cd(fixture.Folder); addpath(exportRoot,"-begin"); clear WVCompiledBackend
            testCase.verifyTrue(startsWith(string(which("WVCompiledBackend")),string(exportRoot)));
            cacheRoot = fullfile(fixture.Folder,"cache");
            capabilities = WVCompiledBackend.capabilitiesForTesting(struct("PackageRoot",string(exportRoot),"CacheRoot",string(cacheRoot),"Architecture","glnxa64","OperatingSystem","Linux","Release","R2026a"));
            testCase.verifyEqual(capabilities.cache.root,string(cacheRoot));
            testCase.verifyFalse(isfolder(cacheRoot));
            clear cleanup
        end
    end

    methods (Test, TestTags="optional")
        function nativeBuildHasExactIdentitiesSelfTestsLifecycleAndRollback(testCase)
            if ~isCanonicalNativePlatform
                capabilities = WVCompiledBackend.build();
                testCase.verifyEqual(capabilities.status,"unsupported");
                testCase.verifyFalse(capabilities.isAvailable);
                return
            end
            capabilities = WVCompiledBackend.build();
            testCase.assertEqual(capabilities.status,"available",capabilities.failure.message);
            testCase.verifyEqual(capabilities.schemaVersion,"1.0.0");
            testCase.verifyTrue(capabilities.matlab.isSupported);
            testCase.verifyEqual(capabilities.matlab.minimumRelease,"R2025b");
            testCase.verifyEqual(capabilities.platform.architecture,"maca64");
            testCase.verifyEqual(capabilities.contract.threadCount,min(18,maxNumCompThreads));
            testCase.verifyEqual(capabilities.contract.planCountMeaning,"logical-prepared-operation-slots");
            testCase.verifyEqual(capabilities.module.executionScheduleVersion,1);
            testCase.verifyEqual(string(capabilities.module.workerPolicyIdentifier),"constant-stage-workers-v1");
            testCase.verifyEqual(string(capabilities.module.planCountMeaning),"logical-prepared-operation-slots");
            testCase.verifyEqual(string(capabilities.module.nonlinearFluxSchedule),"retained-compact-streamed-target-three-channel-v1");
            testCase.verifyEqual(capabilities.libraries.base.path,string(realpath(fullfile(capabilities.cache.root,"provider","native-neon-pthreads","lib","libfftw3.3.dylib"))));
            testCase.verifyEqual(capabilities.libraries.thread.path,string(realpath(fullfile(capabilities.cache.root,"provider","native-neon-pthreads","lib","libfftw3_threads.3.dylib"))));
            testCase.verifyFalse(capabilities.libraries.openmp.detected);
            testCase.verifyLessThanOrEqual(capabilities.featureValidation.maximumRelativeError,1e-12);
            testCase.verifyTrue(capabilities.featureValidation.hydrostatic.lifecyclePassed);
            testCase.verifyTrue(capabilities.featureValidation.nonhydrostatic.lifecyclePassed);
            testCase.verifyGreaterThan(capabilities.featureValidation.hydrostatic.planCount,0);
            testCase.verifyGreaterThan(capabilities.featureValidation.nonhydrostatic.planCount,0);

            missingSymbols = WVCompiledBackend.buildForTesting(struct("CommandOutputFunction",@missingSymbolOutput));
            testCase.verifyEqual(missingSymbols.status,"build-failed");
            testCase.verifyEqual(missingSymbols.buildAttempt.stage,"provider-validation");
            testCase.verifyEqual(missingSymbols.failure.identifier,"WaveVortexModel:CompiledBackendSymbolMissing");

            buildFailure = WVCompiledBackend.buildForTesting(struct("FailureStage","provider-build"));
            testCase.verifyEqual(buildFailure.status,"build-failed");
            testCase.verifyEqual(buildFailure.buildAttempt.stage,"provider-build");
            testCase.verifyEqual(buildFailure.failure.identifier,"WaveVortexModel:CompiledBackendInjectedFailure");

            recordPath = fullfile(capabilities.cache.root,"state","validated-build.json");
            originalRecord = string(fileread(recordPath));
            recordCleanup = onCleanup(@()writeText(recordPath,originalRecord));
            record = jsondecode(originalRecord); record.libraries.base.path = record.libraries.thread.path;
            writeText(recordPath,jsonencode(record,PrettyPrint=true));
            identityMismatch = WVCompiledBackend.capabilities();
            testCase.verifyEqual(identityMismatch.status,"invalid");
            testCase.verifyEqual(identityMismatch.failure.identifier,"WaveVortexModel:CompiledBackendIdentityMismatch");
            clear recordCleanup

            originalHash = capabilities.module.sha256;
            recordBeforeRollback = string(fileread(recordPath));
            rollback = WVCompiledBackend.buildForTesting(struct("FailureStage","install-validation"));
            testCase.verifyEqual(rollback.status,"build-failed");
            testCase.verifyEqual(string(fileread(recordPath)),recordBeforeRollback);
            afterRollback = WVCompiledBackend.capabilities();
            testCase.verifyEqual(afterRollback.status,"available");
            testCase.verifyEqual(afterRollback.module.sha256,originalHash);
            testCase.verifyFalse(afterRollback.module.loadedAfterInspection);
        end

        function compiledBackendExecutesWithoutFallback(testCase)
            if ~isCanonicalNativePlatform
                capabilities = WVCompiledBackend.capabilities();
                testCase.verifyEqual(capabilities.status,"unsupported");
                testCase.verifyFalse(capabilities.isAvailable);
                return
            end
            capabilities = WVCompiledBackend.capabilities();
            if ~capabilities.isAvailable
                capabilities = WVCompiledBackend.build();
            end
            testCase.assertTrue(capabilities.isAvailable,capabilities.failure.message);
            definitions = struct("isHydrostatic",{true,false},"seed",{17201,17202});
            for iDefinition = 1:numel(definitions)
                definition = definitions(iDefinition);
                matlabWVT = WVTransformConstantStratification([15000 12000 1300],[16 12 9],isHydrostatic=definition.isHydrostatic,shouldAntialias=true);
                compiledWVT = WVTransformConstantStratification([15000 12000 1300],[16 12 9],isHydrostatic=definition.isHydrostatic,shouldAntialias=true,computationalBackend="compiled");
                transformCleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                rng(definition.seed,"twister");
                matlabWVT.initWithRandomFlow(uvMax=0.01);
                compiledWVT.Ap = matlabWVT.Ap; compiledWVT.Am = matlabWVT.Am; compiledWVT.A0 = matlabWVT.A0;
                matlabWVT.t = 31; compiledWVT.t = 31;
                [expectedFp,expectedFm,expectedF0] = matlabWVT.nonlinearFlux();
                [actualFp,actualFm,actualF0] = compiledWVT.nonlinearFlux();
                maximumError = max([relativeError(actualFp,expectedFp) relativeError(actualFm,expectedFm) relativeError(actualF0,expectedF0)]);
                testCase.verifyLessThanOrEqual(maximumError,1e-12);
                metadata = compiledWVT.computationalBackendMetadata;
                testCase.verifyEqual(compiledWVT.computationalBackend,"compiled");
                testCase.verifyEqual(metadata.activeBackend,"compiled");
                testCase.verifyEqual(metadata.provider.id,"native-neon-pthreads");
                testCase.verifyTrue(metadata.module.identityValidated);
                testCase.verifyEqual(metadata.contract.version,capabilities.contract.version);
                testCase.verifyEqual(metadata.contract.matlabTransformBridgeVersion,capabilities.module.matlabTransformBridgeVersion);
                testCase.verifyGreaterThanOrEqual(metadata.contract.matlabTransformBridgeVersion,5);
                testCase.verifyGreaterThan(metadata.runtimeMetrics.engineBytes,0);
                testCase.verifyGreaterThan(metadata.runtimeMetrics.kernelBytes,0);
                testCase.verifyGreaterThan(metadata.runtimeMetrics.fieldServiceBytes,0);
                testCase.verifyEqual(metadata.runtimeMetrics.duplicateExecutions,0);

                if definition.isHydrostatic
                    resized = compiledWVT.waveVortexTransformWithResolution([18 14 11]);
                    resizedCleanup = onCleanup(@()delete(resized));
                    testCase.verifyEqual(resized.computationalBackend,"compiled");
                    explicit = compiledWVT.waveVortexTransformWithExplicitAntialiasing();
                    explicitCleanup = onCleanup(@()delete(explicit));
                    testCase.verifyEqual(explicit.computationalBackend,"compiled");
                    testCase.verifyFalse(explicit.shouldAntialias);
                    testCase.verifyTrue(explicit.hasForcingWithName("antialias filter"));
                    clear explicitCleanup
                    clear resizedCleanup

                    restartPath = string(tempname)+".nc";
                    restartCleanup = onCleanup(@()deleteIfPresent(restartPath));
                    ncfile = compiledWVT.writeToFile(char(restartPath),shouldOverwriteExisting=true);
                    ncfile.close();
                    restoredMatlab = WVTransform.waveVortexTransformFromFile(char(restartPath));
                    restoredCompiled = WVTransform.waveVortexTransformFromFile(char(restartPath),computationalBackend="compiled");
                    restorationCleanup = onCleanup(@()deleteTransforms(restoredMatlab,restoredCompiled));
                    testCase.verifyEqual(restoredMatlab.computationalBackend,"matlab");
                    testCase.verifyEqual(restoredCompiled.computationalBackend,"compiled");
                    clear restorationCleanup restartCleanup
                end

                matlabWVT.removeAllForcing();
                compiledWVT.removeAllForcing();
                [expectedFp,expectedFm,expectedF0] = matlabWVT.nonlinearFlux();
                [actualFp,actualFm,actualF0] = compiledWVT.nonlinearFlux();
                verifyExactComplexZero(testCase,actualFp,expectedFp);
                verifyExactComplexZero(testCase,actualFm,expectedFm);
                verifyExactComplexZero(testCase,actualF0,expectedF0);
                clear transformCleanup
            end
            metrics = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(metrics.kernelCount,0);
            testCase.verifyEqual(metrics.activePlans,0);
            testCase.verifyEqual(metrics.outstandingPlanningBytes,0);
        end

        function capabilityInspectionPreservesForeignLoadedModule(testCase)
            if ~isCanonicalNativePlatform
                return
            end
            capabilities = WVCompiledBackend.capabilities();
            if ~capabilities.isAvailable
                capabilities = WVCompiledBackend.build();
            end
            testCase.assertTrue(capabilities.isAvailable,capabilities.failure.message);
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            foreignModule = fullfile(fixture.Folder,capabilities.module.name+"."+mexext);
            copyfile(capabilities.module.path,foreignModule);
            originalPath = path;
            cleanup = onCleanup(@()clearCompiledModuleAndRestorePath(originalPath));
            addpath(fixture.Folder,"-begin");
            info = feval(char(capabilities.module.name),'moduleInfo');
            testCase.verifyTrue(isfield(info,"engine"));
            testCase.assertEqual(loadedCompiledModulePath(capabilities.module.name),string(realpath(foreignModule)));

            inspected = WVCompiledBackend.capabilities();
            testCase.verifyEqual(inspected.status,"invalid");
            testCase.verifyEqual(inspected.failure.identifier,"WaveVortexModel:CompiledBackendModuleConflict");
            testCase.verifyEqual(loadedCompiledModulePath(capabilities.module.name),string(realpath(foreignModule)));
            testCase.verifyEqual(string(which(capabilities.module.name)),string(foreignModule));
            clear cleanup
        end
    end
end

function value = isCanonicalNativePlatform
value = ismac && computer("arch") == "maca64";
end

function value = fakeCompilerRecord(sdkPath)
value = struct("c","/usr/bin/clang","cxx","/usr/bin/clang++","identity","Apple Clang test double","sdkPath",string(sdkPath),"mexName","test","mexManufacturer","test","mexVersion","test");
end

function deleteTransforms(varargin)
for iTransform = 1:numel(varargin)
    if ~isempty(varargin{iTransform}) && isvalid(varargin{iTransform})
        delete(varargin{iTransform});
    end
end
end

function value = relativeError(actual,expected)
value = max(abs(actual(:)-expected(:)),[],'omitmissing')/max(max(abs(expected(:)),[],'omitmissing'),realmin);
end

function verifyExactComplexZero(testCase,actual,expected)
testCase.verifyEqual(real(actual),expected);
testCase.verifyEqual(imag(actual),zeros(size(expected)));
end

function deleteIfPresent(pathname)
if isfile(pathname)
    delete(pathname);
end
end

function unexpectedDownload(~,~)
error("WaveVortexModel:UnexpectedDownload","Unsupported capability detection attempted a download.");
end

function failedDownload(~,~)
error("WaveVortexModel:TestDownloadFailure","Injected network failure.");
end

function [status,output] = missingSymbolOutput(~)
status = 0; output = "";
end

function output = runGit(repositoryRoot,arguments)
[status,output] = system("git -C "+shellQuote(repositoryRoot)+" "+arguments);
if status ~= 0, error("WaveVortexModel:TestGit","git %s failed.",arguments); end
end

function pathname = realpath(pathname)
[status,pathname] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0, error("WaveVortexModel:TestPath","Unable to resolve path."); end
pathname = strtrim(pathname);
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function restoreState(directory,originalPath,originalRng,originalWarnings)
cd(directory); path(originalPath); rng(originalRng); warning(originalWarnings);
end

function restoreDirectoryAndPath(directory,originalPath)
cd(directory); path(originalPath);
end

function writeText(pathname,value)
fileId = fopen(pathname,"w");
if fileId < 0, error("WaveVortexModel:TestWrite","Unable to write %s.",pathname); end
cleanup = onCleanup(@()fclose(fileId)); fprintf(fileId,"%s",value); clear cleanup
end

function pathname = loadedCompiledModulePath(moduleName)
[~,mexFiles] = inmem("-completenames");
paths = string(mexFiles);
paths = paths(endsWith(paths,filesep+moduleName+"."+mexext));
if numel(paths) ~= 1, pathname = ""; else, pathname = string(realpath(paths)); end
end

function clearCompiledModuleAndRestorePath(originalPath)
clear wv_compiled_backend_mex
path(originalPath);
end
