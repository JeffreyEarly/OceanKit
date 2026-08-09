classdef TestReleaseWorkflow < matlab.unittest.TestCase
    properties
        TemporaryRoot
        OceanKitRoot
    end

    methods (TestMethodSetup)
        function createTemporaryRoot(testCase)
            testCase.OceanKitRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.TemporaryRoot = string(tempname);
            mkdir(testCase.TemporaryRoot);
            testCase.addTeardown(@()rmdir(testCase.TemporaryRoot,"s"));
            addpath(fullfile(testCase.OceanKitRoot,"tools"));
            testCase.addTeardown(@()rmpath(fullfile(testCase.OceanKitRoot,"tools")));
            warningState = warning("query","mpm:repository:UnableToConnectToRepository");
            warning("off","mpm:repository:UnableToConnectToRepository");
            testCase.addTeardown(@()warning(warningState));
            pathWarningState = warning("query","MATLAB:mpath:uninstalledPackagesOnPath");
            warning("off","MATLAB:mpath:uninstalledPackagesOnPath");
            testCase.addTeardown(@()warning(pathWarningState));
        end
    end

    methods (Test)
        function promoteUnreleasedPreservesMarkdown(testCase)
            changelogPath = fullfile(testCase.TemporaryRoot,"CHANGELOG.md");
            bodyPath = fullfile(testCase.TemporaryRoot,"release-body.md");
            originalBody = "- First `entry` with 100%." + newline + ...
                newline + "  Continuation with ""quotes""." + newline + ...
                newline + "```matlab" + newline + "value = ""unchanged"";" + newline + "```";
            source = "# Version History" + newline + newline + ...
                "## [Unreleased]" + newline + newline + originalBody + newline + newline + ...
                "## [1.2.3] - 2026-01-01" + newline + newline + "- Older." + newline;
            testCase.writeText(changelogPath,source);

            body = oceankitrelease.promoteUnreleased(changelogPath,"1.2.4","2026-08-09",bodyPath);

            testCase.verifyEqual(body,originalBody);
            testCase.verifyEqual(string(fileread(bodyPath)),originalBody + newline);
            promoted = string(fileread(changelogPath));
            testCase.verifyTrue(contains(promoted,"## [Unreleased]" + newline + newline + "## [1.2.4] - 2026-08-09"));
            testCase.verifyTrue(contains(promoted,originalBody));
            testCase.verifyTrue(contains(promoted,"## [1.2.3] - 2026-01-01"));
        end

        function invalidUnreleasedSectionsDoNotChangeFiles(testCase)
            cases = [
                "# Version History" + newline + newline + "## [1.0.0] - 2026-01-01" + newline
                "# Version History" + newline + newline + "## [Unreleased]" + newline + newline + "## [1.0.0] - 2026-01-01" + newline
                "# Version History" + newline + newline + "## [Unreleased]" + newline + newline + "- One" + newline + newline + "## [Unreleased]" + newline + newline + "- Two" + newline
                ];
            expectedIdentifiers = [
                "OceanKitRelease:InvalidUnreleasedSection"
                "OceanKitRelease:EmptyUnreleasedSection"
                "OceanKitRelease:InvalidUnreleasedSection"
                ];
            for iCase = 1:numel(cases)
                changelogPath = fullfile(testCase.TemporaryRoot,"CHANGELOG-" + iCase + ".md");
                bodyPath = fullfile(testCase.TemporaryRoot,"body-" + iCase + ".md");
                testCase.writeText(changelogPath,cases(iCase));
                testCase.verifyError(@()oceankitrelease.promoteUnreleased(changelogPath,"1.0.1","2026-08-09",bodyPath),expectedIdentifiers(iCase));
                testCase.verifyEqual(string(fileread(changelogPath)),cases(iCase));
                testCase.verifyFalse(isfile(bodyPath));
            end
        end

        function ciReleaseCalculatesVersionBumps(testCase)
            bumpTypes = ["none" "patch" "minor" "major"];
            expectedVersions = ["1.3.0" "1.3.1" "1.4.0" "2.0.0"];
            for iBump = 1:numel(bumpTypes)
                repositoryRoot = testCase.createPackageFixture("version-" + bumpTypes(iBump));
                ci_release(rootDir=repositoryRoot,bumpType=bumpTypes(iBump));
                manifest = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
                testCase.verifyEqual(string(manifest.version),expectedVersions(iBump));
            end
        end

        function legacyNotesBehaviorRemainsAvailable(testCase)
            repositoryRoot = testCase.createPackageFixture("legacy-notes");
            changelogPath = fullfile(repositoryRoot,"CHANGELOG.md");
            testCase.writeText(changelogPath,"# Version History" + newline + newline + "## [1.3.0] - 2026-01-01" + newline + "- Old" + newline);

            ci_release(rootDir=repositoryRoot,bumpType="patch",notes="First" + newline + "Second");

            changelog = string(fileread(changelogPath));
            testCase.verifyTrue(contains(changelog,"## [1.3.1]"));
            testCase.verifyTrue(contains(changelog,"- First" + newline + "- Second"));
            testCase.verifyFalse(contains(changelog,"## [Unreleased]"));
        end

        function fullPilotCandidatePromotesBuildsAndExports(testCase)
            repositoryRoot = testCase.createPackageFixture("full-pilot");
            changelogPath = fullfile(repositoryRoot,"CHANGELOG.md");
            websiteRoot = fullfile(repositoryRoot,"Documentation","WebsiteDocumentation");
            mkdir(websiteRoot);
            testCase.writeText(changelogPath,"# Version History" + newline + newline + ...
                "## [Unreleased]" + newline + newline + "- Multiline `release` entry." + newline + ...
                newline + "  With 100% and ""quotes""." + newline);
            testCase.writeText(fullfile(websiteRoot,"index.md"),"# Fixture documentation" + newline);
            testCase.git(repositoryRoot,"add CHANGELOG.md Documentation/WebsiteDocumentation/index.md");
            testCase.git(repositoryRoot,"commit -m documentation-baseline");
            releaseBodyPath = fullfile(testCase.TemporaryRoot,"pilot-release-body.md");
            exportRoot = fullfile(testCase.TemporaryRoot,"pilot-export");

            ci_release(rootDir=repositoryRoot,bumpType="patch",shouldBuildWebsiteDocumentation=true, ...
                shouldPackageForDistribution=true,shouldPromoteUnreleased=true, ...
                shouldRequireDocumentationBuilder=true,releaseDate="2026-08-09", ...
                releaseBodyPath=releaseBodyPath,outputRoot=exportRoot);

            manifest = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
            testCase.verifyEqual(string(manifest.version),"1.3.1");
            changelog = string(fileread(changelogPath));
            testCase.verifyTrue(contains(changelog,"## [Unreleased]" + newline + newline + "## [1.3.1] - 2026-08-09"));
            testCase.verifyEqual(string(fileread(releaseBodyPath)), ...
                "- Multiline `release` entry." + newline + newline + "  With 100% and ""quotes""." + newline);
            testCase.verifyTrue(contains(string(fileread(fullfile(repositoryRoot,"docs","version-history.md"))), ...
                "## [1.3.1] - 2026-08-09"));
            testCase.verifyTrue(isfolder(fullfile(exportRoot,"ClassDocumentation-1.3.1")));
            testCase.verifyFalse(isfolder(fullfile(repositoryRoot,"dist")));
        end

        function metadataCoversNonDistributionAndDistributionModes(testCase)
            repositoryRoot = testCase.createPackageFixture("metadata");
            metadataPath = fullfile(testCase.TemporaryRoot,"metadata.json");
            metadata = oceankitrelease.writeMetadata(metadataPath, ...
                mode="pilot",workflowRef="workflow@ref",workflowSha="abc",sourceStartSha="def", ...
                repositoryRoot=repositoryRoot,oldVersion="1.3.0",bumpType="none",releaseDate="2026-08-09", ...
                documentationChecked=true,documentationBuilt=false,unreleasedPromoted=false, ...
                shouldPackageForDistribution=false);
            testCase.verifyFalse(metadata.distributionCreated);
            testCase.verifyFalse(metadata.versionChanged);
            testCase.verifyEmpty(metadata.authoringPaths);
            decoded = jsondecode(fileread(metadataPath));
            testCase.verifyEqual(decoded.schemaVersion,1);
            testCase.verifyEmpty(decoded.snapshotFolder);
            testCase.verifyEmpty(decoded.exportPath);

            exportRoot = fullfile(testCase.TemporaryRoot,"export");
            exportFolder = fullfile(exportRoot,"ClassDocumentation-1.3.0");
            mkdir(fullfile(exportFolder,"resources"));
            copyfile(fullfile(repositoryRoot,"resources","mpackage.json"),fullfile(exportFolder,"resources","mpackage.json"));
            metadata = oceankitrelease.writeMetadata(metadataPath, ...
                mode="pilot",workflowRef="workflow@ref",workflowSha="abc",sourceStartSha="def", ...
                repositoryRoot=repositoryRoot,oldVersion="1.3.0",bumpType="none",releaseDate="2026-08-09", ...
                documentationChecked=false,documentationBuilt=false,unreleasedPromoted=false, ...
                shouldPackageForDistribution=true,exportRoot=exportRoot);
            testCase.verifyTrue(metadata.distributionCreated);
            testCase.verifyEqual(metadata.snapshotFolder,"ClassDocumentation-1.3.0");
            testCase.verifyEqual(metadata.exportPath,fullfile(exportRoot,"ClassDocumentation-1.3.0"));
        end

        function metadataCoversVersionOnlyAndFullReleaseModes(testCase)
            repositoryRoot = testCase.createPackageFixture("metadata-release-modes");
            metadataPath = fullfile(testCase.TemporaryRoot,"release-modes.json");
            package = matlab.mpm.Package(repositoryRoot);
            package.Version = matlab.mpm.Version(1,3,1);

            metadata = oceankitrelease.writeMetadata(metadataPath, ...
                mode="pilot",workflowRef="workflow@ref",workflowSha="abc",sourceStartSha="def", ...
                repositoryRoot=repositoryRoot,oldVersion="1.3.0",bumpType="patch",releaseDate="2026-08-09", ...
                documentationChecked=false,documentationBuilt=false,unreleasedPromoted=false, ...
                shouldPackageForDistribution=false);
            testCase.verifyTrue(metadata.versionChanged);
            testCase.verifyEqual(string(metadata.authoringPaths),"resources/mpackage.json");
            encoded = string(fileread(metadataPath));
            testCase.verifyNotEmpty(regexp(encoded,'"authoringPaths"\s*:\s*\[\s*"resources/mpackage.json"\s*\]','once'));

            changelogPath = fullfile(repositoryRoot,"CHANGELOG.md");
            versionHistoryPath = fullfile(repositoryRoot,"docs","version-history.md");
            mkdir(fullfile(repositoryRoot,"docs"));
            testCase.writeText(changelogPath,"# Version History" + newline);
            testCase.writeText(versionHistoryPath,"# Version History" + newline);
            testCase.git(repositoryRoot,"add CHANGELOG.md docs/version-history.md resources/mpackage.json");
            testCase.git(repositoryRoot,"commit -m release-baseline");
            testCase.writeText(changelogPath,"# Version History" + newline + newline + "## [1.3.1]" + newline);
            testCase.writeText(versionHistoryPath,"# Version History" + newline + newline + "## [1.3.1]" + newline);
            package = matlab.mpm.Package(repositoryRoot);
            package.Version = matlab.mpm.Version(1,3,2);
            releaseBodyPath = fullfile(testCase.TemporaryRoot,"full-release-body.md");
            testCase.writeText(releaseBodyPath,"- Released" + newline);
            exportRoot = fullfile(testCase.TemporaryRoot,"full-export");
            exportFolder = fullfile(exportRoot,"ClassDocumentation-1.3.2");
            mkdir(fullfile(exportFolder,"resources"));
            copyfile(fullfile(repositoryRoot,"resources","mpackage.json"),fullfile(exportFolder,"resources","mpackage.json"));

            metadata = oceankitrelease.writeMetadata(metadataPath, ...
                mode="pilot",workflowRef="workflow@ref",workflowSha="abc",sourceStartSha="def", ...
                repositoryRoot=repositoryRoot,oldVersion="1.3.1",bumpType="patch",releaseDate="2026-08-09", ...
                documentationChecked=true,documentationBuilt=true,unreleasedPromoted=true, ...
                shouldPackageForDistribution=true,exportRoot=exportRoot,releaseBodyPath=releaseBodyPath);
            testCase.verifyTrue(metadata.versionChanged);
            testCase.verifyTrue(metadata.documentationChecked);
            testCase.verifyTrue(metadata.documentationBuilt);
            testCase.verifyTrue(metadata.unreleasedPromoted);
            testCase.verifyTrue(metadata.distributionCreated);
            testCase.verifyEqual(sort(string(metadata.authoringPaths)), ...
                sort(["CHANGELOG.md";"docs/version-history.md";"resources/mpackage.json"]));
        end

        function requestedDocumentationTaskExecutesProgrammatically(testCase)
            repositoryRoot = testCase.createPackageFixture("documentation-task");
            markerPath = fullfile(testCase.TemporaryRoot,"documentation-task-ran.txt");
            oldMarker = string(getenv("OCEANKIT_RELEASE_TEST_MARKER"));
            setenv("OCEANKIT_RELEASE_TEST_MARKER",markerPath);
            testCase.addTeardown(@()setenv("OCEANKIT_RELEASE_TEST_MARKER",oldMarker));
            buildFile = [
                "function plan = buildfile"
                "import matlab.buildtool.Task"
                "plan = buildplan;"
                "plan(""docs:check"") = Task(Actions=@documentationCheckTask,DisableIncremental=true);"
                "end"
                ""
                "function documentationCheckTask(~)"
                "fileID = fopen(getenv(""OCEANKIT_RELEASE_TEST_MARKER""),""w"");"
                "assert(fileID >= 0)"
                "fprintf(fileID,""executed\n"");"
                "fclose(fileID);"
                "end"
                ];
            testCase.writeText(fullfile(repositoryRoot,"buildfile.m"),join(buildFile,newline) + newline);
            testCase.git(repositoryRoot,"add buildfile.m");
            testCase.git(repositoryRoot,"commit -m build-task");

            oceankitrelease.runDocumentationCheck(repositoryRoot,"docs:check");

            testCase.verifyEqual(string(fileread(markerPath)),"executed" + newline);
            didThrow = false;
            try
                oceankitrelease.runDocumentationCheck(repositoryRoot,"docs:missing");
            catch
                didThrow = true;
            end
            testCase.verifyTrue(didThrow,"A missing documentation task must fail.");
        end

        function pilotRejectsUnexpectedWorkingTreeDrift(testCase)
            repositoryRoot = testCase.createPackageFixture("drift");
            trackedPath = fullfile(repositoryRoot,"tracked.txt");
            testCase.writeText(trackedPath,"original" + newline);
            testCase.git(repositoryRoot,"add tracked.txt");
            testCase.git(repositoryRoot,"commit -m fixture-update");
            testCase.writeText(trackedPath,"changed" + newline);

            testCase.verifyError(@()oceankitrelease.changedAuthoringPaths(repositoryRoot,"pilot"), ...
                "OceanKitRelease:UnexpectedAuthoringDrift");
        end

        function validationHelpersFailClearly(testCase)
            testCase.verifyError(@()oceankitrelease.installDocumentationPackage("ClassDocumentation"), ...
                "OceanKitRelease:InvalidDocumentationPackageSpecifier");
            dependency = oceankitrelease.installDocumentationPackage("ClassDocumentation@1.3.0", ...
                repositoryRoot=testCase.OceanKitRoot);
            testCase.verifyEqual(dependency.Name,"ClassDocumentation");
            testCase.verifyEqual(dependency.Version,"1.3.0");
            testCase.verifyTrue(isfolder(dependency.Root));
            repositoryRoot = testCase.createPackageFixture("missing-buildfile");
            testCase.verifyError(@()oceankitrelease.runDocumentationCheck(repositoryRoot,"docs:check"), ...
                "OceanKitRelease:BuildFileNotFound");
        end
    end

    methods (Access = private)
        function repositoryRoot = createPackageFixture(testCase,name)
            repositoryRoot = fullfile(testCase.TemporaryRoot,name);
            mkdir(fullfile(repositoryRoot,"resources"));
            sourceManifest = fullfile(testCase.OceanKitRoot,"ClassDocumentation-1.3.0","resources","mpackage.json");
            copyfile(sourceManifest,fullfile(repositoryRoot,"resources","mpackage.json"));
            testCase.git(repositoryRoot,"init -b main");
            testCase.git(repositoryRoot,"config user.name Test");
            testCase.git(repositoryRoot,"config user.email test@example.com");
            testCase.git(repositoryRoot,"add resources/mpackage.json");
            testCase.git(repositoryRoot,"commit -m fixture");
        end

        function git(testCase,repositoryRoot,arguments)
            command = sprintf('git -C "%s" %s',repositoryRoot,arguments);
            [status,output] = system(command);
            testCase.assertEqual(status,0,output);
        end

        function writeText(testCase,path,text)
            fileID = fopen(path,"w");
            testCase.assertGreaterThanOrEqual(fileID,0,"Unable to create test fixture.");
            cleanup = onCleanup(@()fclose(fileID));
            fwrite(fileID,text);
            clear cleanup
        end
    end
end
