classdef TestReleaseVerification < matlab.unittest.TestCase
    properties
        repositoryRoot (1,1) string
    end

    methods (TestMethodSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test,TestTags="full")
        function releaseCallerUsesImmutablePilot(testCase)
            workflow = testCase.readFile(fullfile(".github","workflows","release-mpm.yml"));
            testCase.verifySubstring(workflow, ...
                "JeffreyEarly/OceanKit/.github/workflows/reusable-mpm-release.yml@mpm-release-v0.1.2");
            testCase.verifySubstring(workflow,"documentationPackageSpecifier: ClassDocumentation@1.3.2");
            testCase.verifySubstring(workflow,"shouldCheckWebsiteDocumentation: true");
            testCase.verifySubstring(workflow,"documentationCheckTask: docs:check");
            testCase.verifySubstring(workflow, ...
                "shouldPromoteUnreleased: ${{ github.event.inputs.bump != 'none' }}");
            testCase.verifySubstring(workflow,"shouldBuildWebsiteDocumentation:");
            testCase.verifySubstring(workflow,"shouldPackageForDistribution:");
        end

        function runtimeManifestExcludesAuthoringTests(testCase)
            manifest = jsondecode(testCase.readFile(fullfile("resources","mpackage.json")));
            folders = string({manifest.folders.path});
            testCase.verifyFalse(any(folders == "UnitTests"));
            testCase.verifyFalse(any(folders == "Benchmarks"));
            testCase.verifyFalse(any(folders == "DeveloperExperiments"));

            dependencies = string({manifest.dependencies.name});
            compatibleVersions = string({manifest.dependencies.compatibleVersions});
            expected = dictionary( ...
                ["ClassAnnotations" "InternalModes" "NetCDF" "SplineCore" "chebfun"], ...
                ["^1.2.1" "1.3.0" "^1.0.2" "^2.2.0" ""]);
            testCase.verifyEqual(sort(reshape(dependencies,[],1)),sort(reshape(keys(expected),[],1)));
            for iDependency = 1:numel(dependencies)
                testCase.verifyEqual(compatibleVersions(iDependency),expected(dependencies(iDependency)));
            end
            testCase.verifyEqual(string(manifest.releaseCompatibility),">=R2025b");
            testCase.verifyMatches(string(manifest.version),"^\d+\.\d+\.\d+$");
        end

        function packageWorkflowHasIsolatedR2025bGates(testCase)
            workflow = testCase.readFile(fullfile(".github","workflows","release-verification.yml"));
            tools = testCase.readFile(fullfile("tools","verifyWaveVortexModelPackage.m")) + ...
                testCase.readFile(fullfile("tools","prepareWaveVortexModelReleaseCandidate.m"));
            contract = workflow + tools;
            for required = [
                    "Clean install / MATLAB R2025b"
                    "Exported package / MATLAB R2025b"
                    "release: R2025b"
                    "contents: read"
                    "96f0b801c565406dd5a4ba2480334a3a481c3e2c"
                    "MATLAB_PREFDIR"
                    "Temporary=true"
                    "ClassDocumentation@1.3.2"
                    "prepareWaveVortexModelReleaseCandidate"
                    "verifyWaveVortexModelPackage"
                    "configureCIEnvironment(repositoryRoot,oceanKitRoot)"
                    "ReleaseCompatibility"
                    ">=R2025b"
                    "checkCIWorkspaceCleanliness.sh"
                    ]'
                testCase.verifySubstring(contract,required);
            end
            testCase.verifySubstring(workflow,"workflow_dispatch:");
            testCase.verifySubstring(workflow,"pull_request:");
            testCase.verifyFalse(contains(workflow,newline + "  push:"));
            testCase.verifyFalse(contains(workflow,'- "UnitTests/**"'));
            testCase.verifyFalse(contains(workflow,'- "Benchmarks/**"'));
            testCase.verifyFalse(contains(workflow,'- "DeveloperExperiments/**"'));
            testCase.verifyFalse(contains(workflow,"contents: write"));
            testCase.verifySubstring(workflow,"verifyWaveVortexModelPackage(sourceRoot,oceanKitRoot)");
            testCase.verifySubstring(workflow,"candidate.report.exportPath");
            testCase.verifySubstring(workflow,"candidate.report.version");
            testCase.verifyFalse(contains(workflow,newline + "    paths:"));
            testCase.verifySubstring(workflow,"final-integration");
            testCase.verifySubstring(workflow,"types: [opened, synchronize, reopened, labeled, unlabeled, closed]");
            testCase.verifySubstring(workflow,"Build and inspect the exported reference runtime");
            testCase.verifySubstring(workflow,"-DBUILD_TESTING=OFF");
            testCase.verifyFalse(contains(contract,'expectedVersion="4.2.0"'));
            testCase.verifyFalse(contains(contract,'WaveVortexModel-4.2.1'));
        end

        function workflowTopologySeparatesFocusedAndFinalGates(testCase)
            requiredWorkflow = testCase.readFile(fullfile(".github","workflows","ci.yml"));
            packageWorkflow = testCase.readFile(fullfile(".github","workflows","release-verification.yml"));
            extendedWorkflow = testCase.readFile(fullfile(".github","workflows","extended-ci.yml"));
            focusedJob = testCase.sectionBetween(requiredWorkflow,"  focused-matlab:","  smoke:");

            for required = [
                    "Focused MATLAB / R2025b"
                    "Required / WaveVortexModel"
                    "Portable C++ kernel contract"
                    "final-integration"
                    "buildtool(taskArguments{:})"
                    "WVM_RUN_ANALYZER"
                    "WVM_RUN_DOCUMENTATION"
                    "MATLAB provisioning incident"
                    "timeout-minutes: 10"
                    "cache restore/populate details only in its step log"
                    ]'
                testCase.verifySubstring(requiredWorkflow,required);
            end
            testCase.verifyEqual(count(focusedJob,"matlab-actions/setup-matlab@v3"),1);
            testCase.verifyEqual(count(focusedJob,"matlab-actions/run-command@v3"),1);
            testCase.verifyEqual(count(focusedJob,"buildtool test:smoke"),0);
            testCase.verifySubstring(focusedJob,'tasks = "test:smoke"');
            testCase.verifyEqual(count(focusedJob,"end;"),3);

            for workflow = [requiredWorkflow packageWorkflow extendedWorkflow]
                testCase.verifySubstring(workflow,"workflow_dispatch:");
                testCase.verifySubstring(workflow,"pull_request:");
                testCase.verifySubstring(workflow,"types: [opened, synchronize, reopened, labeled, unlabeled, closed]");
                testCase.verifySubstring(workflow,"cancel-in-progress: true");
                testCase.verifySubstring(workflow,"retain-final");
                testCase.verifySubstring(workflow,"github.base_ref == 'main'");
                testCase.verifySubstring(workflow,"final-integration");
            end
            testCase.verifyFalse(contains(packageWorkflow,newline + "    paths:"));
            testCase.verifyFalse(contains(extendedWorkflow,newline + "    paths:"));
            for context = [
                    "Smoke / MATLAB R2025b"
                    "Documentation / MATLAB R2025b"
                    "Code Analyzer / MATLAB R2025b"
                    "Clean install / MATLAB R2025b"
                    "Exported package / MATLAB R2025b"
                    "Full / MATLAB R2025b"
                    "Exhaustive / MATLAB R2025b"
                    "Optional / MATLAB R2025b"
                    ]'
                testCase.verifyTrue(contains(requiredWorkflow + packageWorkflow + extendedWorkflow,context));
            end
        end

        function portableRuntimeExportContractIsSourceOnly(testCase)
            helper = testCase.readFile("tools/prepareWaveVortexModelReleaseCandidate.m");
            buildScript = testCase.readFile(fullfile("PortableRuntime","buildWaveVortexRun.sh"));
            cmake = testCase.readFile(fullfile("PortableRuntime","CMakeLists.txt"));
            for required = [
                    "CompiledKernel/native-fftw-provider.env"
                    "CompiledKernel/include/WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
                    "CompiledKernel/src/WVTransformBarotropicQGKernel.cpp"
                    "PortableRuntime/README.md"
                    "PortableRuntime/buildWaveVortexRun.sh"
                    "PortableRuntime/source-selection.json"
                    "PortableRuntime/include/WaveVortexRuntime/WVCheckpointReader.hpp"
                    "PortableRuntime/include/WaveVortexRuntime/WVBarotropicQGIntegrationSystem.hpp"
                    "PortableRuntime/include/WaveVortexRuntime/WVExtensionCatalog.hpp"
                    "PortableRuntime/include/WaveVortexRuntime/WVFieldEvaluationService.hpp"
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
                    "PortableRuntime/src/WVCheckpointTransformAdapters.cpp"
                    "PortableRuntime/src/WVBarotropicQGFieldEvaluationAdapter.cpp"
                    "PortableRuntime/src/WVBarotropicQGFieldEvaluationAdapter.hpp"
                    "PortableRuntime/src/WVBarotropicQGIntegrationSystem.cpp"
                    "PortableRuntime/src/WVExtensionCatalog.cpp"
                    "PortableRuntime/src/WVFieldEvaluationService.cpp"
                    "PortableRuntime/src/WVLegacyObserverCompatibility.hpp"
                    "PortableRuntime/src/WVModelInternalAccess.hpp"
                    "PortableRuntime/src/WVModelOutputNetCDFReader.cpp"
                    "PortableRuntime/src/WVModelOutputTransformAdapters.cpp"
                    "PortableRuntime/src/WVModelTransformAdapters.cpp"
                    "PortableRuntime/src/WVObserverOutputEvaluationService.cpp"
                    "PortableRuntime/src/WVOutputTransformAdapters.cpp"
                    "PortableRuntime/src/WVOutputSchedule.cpp"
                    "PortableRuntime/app/WaveVortexRun.cpp"
                    "PortableRuntime/app/WaveVortexRunMain.cpp"
                    ]'
                testCase.verifySubstring(helper,required);
            end
            for required = [
                    "include/WaveVortexRuntime/WVExtensionCatalog.hpp"
                    "include/WaveVortexRuntime/WVRunner.hpp"
                    "src/WVExtensionCatalog.cpp"
                    "app/WaveVortexRunMain.cpp"
                    ]'
                testCase.verifySubstring(cmake,required);
            end
            testCase.verifySubstring(buildScript,'repository_root=$(CDPATH= cd -- "$script_directory/.." && pwd)');
            testCase.verifySubstring(cmake,'WaveVortexModel-${WV_RUNTIME_PACKAGE_VERSION}');
            testCase.verifySubstring(helper,"compiled product or downloaded archive");
            testCase.verifySubstring(helper,'exportedFiles == "wave-vortex-run"');
            testCase.verifySubstring(helper,'exportedFiles == "wave-vortex-run.exe"');
            testCase.verifyFalse(contains(helper,'startsWith(exportedFiles,"wave-vortex-run")'));
        end

        function routineWorkflowsUsePilotDependencySnapshot(testCase)
            requiredWorkflow = testCase.readFile(fullfile(".github","workflows","ci.yml"));
            extendedWorkflow = testCase.readFile(fullfile(".github","workflows","extended-ci.yml"));
            for workflow = [requiredWorkflow extendedWorkflow]
                testCase.verifySubstring(workflow,"96f0b801c565406dd5a4ba2480334a3a481c3e2c");
                testCase.verifySubstring(workflow,"release: R2025b");
                testCase.verifyFalse(contains(workflow,"eb71bc7b05e74776afd678b27c964cf53cd9d547"));
                testCase.verifyFalse(contains(workflow,"R2024b"));
            end
            for context = [
                    "Smoke / MATLAB R2025b"
                    "Documentation / MATLAB R2025b"
                    "Code Analyzer / MATLAB R2025b"
                    ]'
                testCase.verifySubstring(requiredWorkflow,context);
            end
            for context = [
                    "Full / MATLAB R2025b"
                    "Exhaustive / MATLAB R2025b"
                    "Optional / MATLAB R2025b"
                    ]'
                testCase.verifySubstring(extendedWorkflow,context);
            end
        end

        function releaseDocumentationMatchesCanonicalSources(testCase)
            canonicalCI = testCase.readFile(fullfile("Documentation","WebsiteDocumentation", ...
                "developers-guide","continuous-integration.md"));
            generatedCI = testCase.readFile(fullfile("docs","developers-guide","continuous-integration.md"));
            canonicalChecklist = testCase.readFile(fullfile("Documentation","WebsiteDocumentation", ...
                "developers-guide","release-checklist.md"));
            generatedChecklist = testCase.readFile(fullfile("docs","developers-guide","release-checklist.md"));
            readme = testCase.readFile(fullfile(testCase.repositoryRoot,"README.md"));
            installation = testCase.readFile(fullfile("Documentation","WebsiteDocumentation","installation.md"));
            testCase.verifyEqual(generatedCI,canonicalCI);
            testCase.verifyEqual(generatedChecklist,canonicalChecklist);
            for phrase = [
                    "Native-package CI: Clean install and Exported package on MATLAB R2025b"
                    "promote the complete `Unreleased` body"
                    "authoring commit first"
                    "OceanKit package snapshot second"
                    "immutable version tag third"
                    "GitHub release last"
                    ]'
                testCase.verifySubstring(generatedChecklist,phrase);
            end
            testCase.verifySubstring(canonicalCI,"supported runtime and native-package verification use MATLAB R2025b");
            testCase.verifyFalse(contains(canonicalCI," / MATLAB R2025a"));
            testCase.verifySubstring(readme,"requires MATLAB R2025b or newer");
            testCase.verifySubstring(installation,"requires MATLAB R2025b or newer");
            testCase.verifyFalse(contains(readme,"requires MATLAB R2025a"));
            testCase.verifyFalse(contains(installation,"requires MATLAB R2025a"));
        end

        function exportHelperRequiresPinnedDocumentationVersion(testCase)
            helper = testCase.readFile("tools/prepareWaveVortexModelReleaseCandidate.m");
            testCase.verifySubstring(helper,'installDocumentationPackage("ClassDocumentation@1.3.2"');
            testCase.verifySubstring(helper,'dependency.Version == "1.3.2"');
            testCase.verifyFalse(contains(helper,'dependency.Version == "1.3.0"'));
        end

        function generatedVersionHistoryIncludesCurrentReleaseRecords(testCase)
            changelog = testCase.readFile("CHANGELOG.md");
            versionHistory = testCase.readFile(fullfile("docs","version-history.md"));
            releaseGateRecord = "Pinned releases to the immutable OceanKit release pilot";
            releaseFloorRecord = "Raised the minimum MATLAB release to R2025b";
            compactMappingRecord = "Replaced vertically replicated horizontal Fourier mappings";
            for record = [releaseGateRecord releaseFloorRecord compactMappingRecord]
                testCase.verifySubstring(changelog,record);
                testCase.verifySubstring(versionHistory,record);
            end
        end
    end

    methods (Access=private)
        function section = sectionBetween(testCase,contents,startMarker,endMarker)
            startIndex = strfind(contents,startMarker);
            endIndex = strfind(contents,endMarker);
            testCase.assertNotEmpty(startIndex);
            testCase.assertNotEmpty(endIndex);
            testCase.assertLessThan(startIndex(1),endIndex(1));
            section = extractBetween(contents,startIndex(1),endIndex(1)-1);
        end

        function contents = readFile(testCase,path)
            if ~isfile(path)
                path = fullfile(testCase.repositoryRoot,path);
            end
            contents = string(fileread(path));
        end
    end
end
