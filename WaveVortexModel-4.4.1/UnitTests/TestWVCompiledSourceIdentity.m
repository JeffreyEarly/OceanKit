classdef TestWVCompiledSourceIdentity < matlab.unittest.TestCase
    properties (SetAccess=private)
        repositoryRoot (1,1) string
    end

    methods (TestClassSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test,TestTags="smoke")
        function manifestIsPathIndependentAndTracksEditAddDelete(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            firstRoot = fullfile(fixture.Folder,"first");
            secondRoot = fullfile(fixture.Folder,"second");
            populateSourceTree(firstRoot);
            copyfile(firstRoot,secondRoot);
            first = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(firstRoot)));
            repeated = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(firstRoot)));
            relocated = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(secondRoot)));
            testCase.verifyEqual(first,repeated);
            testCase.verifyEqual(first.aggregateSHA256,relocated.aggregateSHA256);
            testCase.verifyEqual({first.files.path},{relocated.files.path});

            source = fullfile(firstRoot,"CompiledKernel","src","alpha.cpp");
            writelines("changed",source);
            edited = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(firstRoot)));
            testCase.verifyNotEqual(edited.aggregateSHA256,first.aggregateSHA256);
            addedPath = fullfile(firstRoot,"PortableRuntime","include","added.hpp");
            writelines("added",addedPath);
            added = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(firstRoot)));
            testCase.verifyNotEqual(added.aggregateSHA256,edited.aggregateSHA256);
            delete(addedPath);
            deleted = WVCompiledBackend.sourceManifestForTesting(struct("PackageRoot",string(firstRoot)));
            testCase.verifyEqual(deleted.aggregateSHA256,edited.aggregateSHA256);
        end
    end

    methods (Test,TestTags="optional")
        function validatedBinaryAndSourceIdentityQualify(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            qualification = WVCompiledBackend.qualifySourceIdentity();
            testCase.verifyTrue(qualification.matches);
            testCase.verifyEqual(string(qualification.currentManifest.aggregateSHA256),string(qualification.manifest.aggregateSHA256));
            testCase.verifyEqual(string(qualification.moduleSHA256),capabilities.module.sha256);
            testCase.verifyEqual(string(qualification.baseLibrarySHA256),capabilities.libraries.base.sha256);
            testCase.verifyEqual(string(qualification.threadLibrarySHA256),capabilities.libraries.thread.sha256);
        end

        function samePathHashReplacementIsRejected(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            recordPath = fullfile(capabilities.cache.root,"state","validated-build.json");
            original = string(fileread(recordPath));
            cleanup = onCleanup(@()writeText(recordPath,original));
            record = jsondecode(original);
            record.module.sha256 = string(repmat('0',1,64));
            writeText(recordPath,jsonencode(record,PrettyPrint=true));
            changedModule = WVCompiledBackend.capabilities();
            testCase.verifyEqual(changedModule.status,"invalid");
            testCase.verifyEqual(changedModule.failure.identifier,"WaveVortexModel:CompiledBackendIdentityMismatch");
            writeText(recordPath,original);
            record = jsondecode(original);
            record.libraries.base.sha256 = string(repmat('f',1,64));
            writeText(recordPath,jsonencode(record,PrettyPrint=true));
            changedProvider = WVCompiledBackend.capabilities();
            testCase.verifyEqual(changedProvider.status,"invalid");
            testCase.verifyEqual(changedProvider.failure.identifier,"WaveVortexModel:CompiledBackendIdentityMismatch");
            clear cleanup
        end

    end
end

function populateSourceTree(root)
paths = [ ...
    "CompiledKernel/src/alpha.cpp" ...
    "CompiledKernel/include/WaveVortexKernel/alpha.hpp" ...
    "CompiledKernel/adapters/native-fftw/adapter.cpp" ...
    "CompiledKernel/adapters/accelerate/matrix.cpp" ...
    "PortableRuntime/src/runtime.cpp" ...
    "PortableRuntime/include/WaveVortexRuntime/runtime.hpp" ...
    "@WVCompiledBackend/WVCompiledBackend.m" ...
    "@WVCompiledTransformBackend/WVCompiledTransformBackend.m" ...
    "CompiledKernel/native-fftw-provider.env"];
for index = 1:numel(paths)
    pathname = fullfile(root,paths(index));
    if ~isfolder(fileparts(pathname)), mkdir(fileparts(pathname)); end
    writelines("content "+index,pathname);
end
end

function writeText(pathname,value)
file = fopen(pathname,"w");
if file < 0, error("Unable to write test record."); end
cleanup = onCleanup(@()fclose(file));
fprintf(file,"%s",value);
clear cleanup
end
