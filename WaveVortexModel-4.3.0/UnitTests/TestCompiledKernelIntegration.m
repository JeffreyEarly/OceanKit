classdef TestCompiledKernelIntegration < matlab.unittest.TestCase
    properties
        repositoryRoot
        selection
        runtimeSelection
    end

    methods (TestClassSetup)
        function loadSelection(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.selection = jsondecode(fileread(fullfile(testCase.repositoryRoot,"CompiledKernel","source-selection.json")));
            testCase.runtimeSelection = jsondecode(fileread(fullfile(testCase.repositoryRoot,"PortableRuntime","source-selection.json")));
        end
    end

    methods (Test,TestTags="smoke")
        function selectionRecordsTheAcceptedContract(testCase)
            testCase.verifyEqual(string(testCase.selection.schemaVersion),"1.0.0");
            testCase.verifyEqual(string(testCase.selection.selectedFeatureSnapshot),"8d0b49236c703dfa7230a8875022fdb3e30283b0");
            testCase.verifyEqual(string(testCase.selection.preparation.sharedCoefficientFormulas),"97e28a892f55877a9267b2c4e3ca399b85596a34");
            testCase.verifyEqual(string(testCase.selection.preparation.nativeFFTWProvider),"1854d72ec07b0fb58c2b2aa972250d91629e669d");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.sharedRightHandSideEvaluation),"bcbf315dc46d2becc5964300ccf252452a9768d5");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.eagerScalarAdvectionPlanning),"issue-218");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.adaptiveRK23MatlabParity),"issue-240");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.wvModelRuntimeFacade),"issue-248");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.resolvedObservingSystems),"issue-258");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.observationBatches),"issue-260");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.resolvedForcings),"issue-261");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.linearBottomFriction),"issue-262");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.explicitExtensionCatalog),"issue-271");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.eventDependentObservationGeometry),"issue-272");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.losslessOutputGraphs),"issue-273");
            testCase.verifyEqual(string(testCase.selection.postSelectionExtensions.stableSourceAPI),"issue-249");
            testCase.verifyEqual(testCase.selection.contract.version,4);
            testCase.verifyEqual(string(testCase.selection.contract.coefficientShape),"[Nj,Nkl]");
            testCase.verifyEqual(testCase.selection.contract.planCount,17);
            testCase.verifyEqual(string(testCase.selection.contract.knownScratch),"4H+6R");
            testCase.verifyEqual(testCase.selection.contract.coefficientWorkers,2);
            testCase.verifyFalse(testCase.selection.contract.persistentFullHermitianSpectrum);
        end

        function runtimeSelectionRecordsStableSourceAPI(testCase)
            runtimeSelection = testCase.runtimeSelection;
            testCase.verifyEqual(string(runtimeSelection.schemaVersion),"1.0.0");
            testCase.verifyEqual(string(runtimeSelection.sourceAPI.identifier),"wave-vortex-portable-source-api-v1");
            testCase.verifyEqual(runtimeSelection.sourceAPI.major,1);
            testCase.verifyEqual(runtimeSelection.sourceAPI.minor,0);
            testCase.verifyEqual(string(runtimeSelection.sourceAPI.compatibility),"source-only");

            expectedInputs = [
                "PortableRuntime/CMakeLists.txt"
                "PortableRuntime/README.md"
                "PortableRuntime/buildWaveVortexRun.sh"
                "PortableRuntime/source-selection.json"
                "PortableRuntime/include"
                "PortableRuntime/src"
                "PortableRuntime/app"
                "PortableRuntime/contracts"
                "PortableRuntime/examples"
                "PortableRuntime/third_party/nlohmann"
                ];
            expectedHeaders = [
                "PortableRuntime/include/WaveVortexRuntime/WVPortableImplementationContract.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVPortableTypedRecord.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVExtensionCatalog.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObservingSystem.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObserverOutputProvider.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObservation.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVOutputSchedule.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVForcing.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVForcingContracts.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVModelOutputConfiguration.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVModel.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVRunner.hpp"
                ];
            actualInputs = reshape(string(runtimeSelection.retained.portableRuntimeInputs),[],1);
            actualHeaders = reshape(string(runtimeSelection.sourceAPI.publicHeaders),[],1);
            testCase.verifyEqual(sort(actualInputs),sort(expectedInputs));
            testCase.verifyEqual(sort(actualHeaders),sort(expectedHeaders));
            for relativePath = [actualInputs;actualHeaders]'
                retainedPath = fullfile(testCase.repositoryRoot,relativePath);
                testCase.verifyTrue(isfile(retainedPath) || isfolder(retainedPath),relativePath);
            end

            provenanceFields = [
                "sourceAPIVersion"
                "preFacadeBaseline"
                "preGeneralizationBaseline"
                "externalConsumerBaseline"
                "externalConsumerBaselineCommit"
                "externalConsumerBaselineWaveVortexCommit"
                ];
            expectedProvenance = [
                "1.0"
                "5940f4e4e206fbfb9a6cb4760af2ba2347e3571f"
                "900f98eb28629a56b7a2d26da38c597c7a709a96"
                "satmapkit/AlongTrackSimulator"
                "ba57981f336ad5bbbc0907dcd74fcd4fcd137708"
                "84dc1b53093650fdbc9be62b637cf298c46a95d3"
                ];
            runtimeExtensions = runtimeSelection.postSelectionExtensions;
            kernelExtensions = testCase.selection.postSelectionExtensions;
            testCase.verifyEqual(string(runtimeExtensions.losslessOutputGraphs),"issue-273");
            testCase.verifyEqual(string(runtimeExtensions.stableSourceAPI),"issue-249");
            testCase.verifyEqual(string(runtimeExtensions.barotropicQGKernel),"issue-280");
            testCase.verifyEqual(string(runtimeExtensions.barotropicQGForcing),"issue-281");
            testCase.verifyEqual(string(runtimeExtensions.barotropicQGModelOutputRestart),"issue-282");
            for iField = 1:numel(provenanceFields)
                field = provenanceFields(iField);
                testCase.verifyEqual(string(runtimeExtensions.(field)),expectedProvenance(iField),field);
                testCase.verifyEqual(string(kernelExtensions.(field)),expectedProvenance(iField),field);
            end
        end

        function keyRuntimeHashesMatchTheSelection(testCase)
            entries = {
                "CompiledKernel/src/WVCoefficientFormulas.hpp","coefficientFormulas"
                "CompiledKernel/src/WVTransformConstantStratificationKernel.cpp","kernel"
                "PortableRuntime/src/WVConstantStratificationIntegrationSystem.cpp","integrationSystem"
                "CompiledKernel/adapters/native-fftw/WVNativeFFTWEngine.cpp","nativeEngine"
                "CompiledKernel/adapters/native-fftw/wv_compiled_backend_mex.cpp","mexGateway"
                "@WVCompiledBackend/private/wvCompiledBackendBuild.m","buildOrchestrator"
                };
            for iEntry = 1:size(entries,1)
                actual = sha256File(fullfile(testCase.repositoryRoot,entries{iEntry,1}));
                testCase.verifyEqual(actual,string(testCase.selection.keySourceSHA256.(entries{iEntry,2})),entries{iEntry,1});
            end
        end
    end

    methods (Test,TestTags="full")
        function portableCoreHasNoHostOrVendorDependency(testCase)
            coreFiles = [
                filesUnder(fullfile(testCase.repositoryRoot,"CompiledKernel","include"),[".hpp" ".h"])
                filesUnder(fullfile(testCase.repositoryRoot,"CompiledKernel","src"),[".cpp" ".hpp" ".h"])
                ];
            forbidden = ["mex.h" "matrix.h" "MatlabDataArray" "fftw3.h" "netcdf.h" "Accelerate/Accelerate.h"];
            for file = coreFiles'
                source = string(fileread(file));
                testCase.verifyFalse(any(contains(source,forbidden)),"Portable dependency found in "+file);
            end
        end

        function repositoryTracksOnlySourceProducts(testCase)
            tracked = splitlines(strtrim(runGit(testCase.repositoryRoot,"ls-files")));
            required = [
                "CompiledKernel/CMakeLists.txt"
                "CompiledKernel/source-selection.json"
                "CompiledKernel/include/WaveVortexKernel/WVFFTEngine.hpp"
                "CompiledKernel/src/WVTransformConstantStratificationKernel.cpp"
                "CompiledKernel/adapters/native-fftw/wv_compiled_backend_mex.cpp"
                "CompiledKernel/native-fftw-provider.env"
                "CompiledKernel/adapters/reference/WVReferenceFFTEngine.cpp"
                "PortableRuntime/CMakeLists.txt"
                "PortableRuntime/README.md"
                "PortableRuntime/source-selection.json"
                "PortableRuntime/buildWaveVortexRun.sh"
                "PortableRuntime/include/WaveVortexRuntime/WVExtensionCatalog.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVForcing.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVForcingContracts.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVModel.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVModelOutputConfiguration.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObservation.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObserverOutputProvider.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVObservingSystem.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVOutputSchedule.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVPortableImplementationContract.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVPortableTypedRecord.hpp"
                "PortableRuntime/include/WaveVortexRuntime/WVRunner.hpp"
                "PortableRuntime/src/WVExtensionCatalog.cpp"
                "PortableRuntime/src/WVLegacyObserverCompatibility.hpp"
                "PortableRuntime/src/WVModelInternalAccess.hpp"
                "PortableRuntime/app/WaveVortexRun.cpp"
                "PortableRuntime/app/WaveVortexRunMain.cpp"
                "PortableRuntime/app/WVRunRequest.cpp"
                "PortableRuntime/app/WVRunRequest.hpp"
                "PortableRuntime/contracts/wave-vortex-run-request-v1.schema.json"
                "PortableRuntime/contracts/wave-vortex-run-request-v2.schema.json"
                "PortableRuntime/examples/portable-run-request-v1.json"
                "PortableRuntime/examples/portable-run-request-v2.json"
                "PortableRuntime/third_party/nlohmann/json.hpp"
                "PortableRuntime/third_party/nlohmann/LICENSE.MIT"
                "@WVCompiledBackend/WVCompiledBackend.m"
                ];
            testCase.verifyTrue(all(ismember(required,tracked)));
            forbidden = contains(tracked,".compiled-backend-cache/") | ...
                endsWith(tracked,[".a" ".o" ".dylib" ".tar.gz"]) | ...
                ~cellfun(@isempty,regexp(cellstr(tracked),'\.mex[a-z0-9]*$','once'));
            testCase.verifyFalse(any(forbidden),"Generated native products are tracked.");
            testCase.verifyFalse(any(startsWith(tracked,"Benchmarks/results/experiments/issue1")),"Engineering benchmark artifacts entered the integration branch.");
        end
    end
end

function files = filesUnder(root,extensions)
listing = dir(fullfile(root,"**","*"));
listing = listing(~[listing.isdir]);
paths = string(fullfile({listing.folder},{listing.name}))';
[~,~,extension] = fileparts(paths);
files = paths(ismember(extension,extensions));
end

function hash = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status ~= 0, error("WaveVortexModel:TestHash","Unable to hash %s.",pathname); end
hash = extractBefore(string(strtrim(output))," ");
end

function output = runGit(repositoryRoot,arguments)
[status,output] = system("git -C "+shellQuote(repositoryRoot)+" "+arguments);
if status ~= 0, error("WaveVortexModel:TestGit","git %s failed.",arguments); end
output = string(output);
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
