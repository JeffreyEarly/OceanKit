classdef TestProductionCodeAnalyzer < matlab.unittest.TestCase
    properties
        repositoryRoot
        productionReport
        temporaryFolder
    end

    methods (TestClassSetup)
        function analyzeProductionSource(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.repositoryRoot,"tools")));
            testCase.productionReport = analyzeProductionCode(testCase.repositoryRoot,ShouldPrint=false,ShouldFail=false);
        end
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function productionInventoryIsDeterministic(testCase)
            files = testCase.productionReport.Files;
            testCase.verifyNumElements(files,172);
            testCase.verifyEqual(files,sort(unique(files)));
            testCase.verifyTrue(all(isfile(fullfile(testCase.repositoryRoot,files))));

            expectedFiles = [
                "WVOperation.m"
                "@WVTransform/WVTransform.m"
                "FastTransforms/WVFourierStorageLayout.m"
                "FastTransforms/@WVFastTransformDoublyPeriodicFFTW/WVFastTransformDoublyPeriodicFFTW.m"
                "Forcing/WVNonlinearAdvection.m"
                "ObservingSystems/WVLagrangianParticles.m"
                "Integrators/WVModelAdaptiveTimeStepMethods.m"
                ];
            testCase.verifyTrue(all(ismember(expectedFiles,files)));
        end

        function authoringAndLegacyAreasAreExcluded(testCase)
            files = testCase.productionReport.Files;
            excludedPrefixes = [
                "UnitTests/"
                "tools/"
                "Documentation/"
                "docs/"
                "Benchmarks/"
                "DeveloperExperiments/"
                "OceanKit/"
                ".buildtool/"
                ];
            for prefix = excludedPrefixes'
                testCase.verifyFalse(any(startsWith(files,prefix)),"Unexpected analyzer input under " + prefix);
            end
            testCase.verifyFalse(any(files == "buildfile.m"));
        end

        function productionFindingsAreClassifiedAndNonblocking(testCase)
            findings = testCase.productionReport.Findings;
            testCase.verifyEmpty(testCase.productionReport.BlockingFindings);
            testCase.verifyFalse(any(startsWith(findings.Classification,"blocking")));
            testCase.verifyTrue(any(findings.CheckID == "AGROW" & findings.Classification == "performance"));
            testCase.verifyTrue(any(findings.CheckID == "INUSD" & findings.Classification == "style"));
            testCase.verifyTrue(any(findings.CheckID == "CTOINW" & findings.Classification == "accepted-false-positive"));
            testCase.verifyTrue(any(findings.CheckID == "MCNPR" & findings.Classification == "accepted-false-positive"));
            testCase.verifyFalse(any(findings.CheckID == "NASGU"));
        end

        function reportContainsReleaseLocationsAndDiagnostics(testCase)
            output = evalc("analyzeProductionCode(testCase.repositoryRoot,ShouldFail=false);");
            testCase.verifySubstring(output,"MATLAB Code Analyzer: release=R");
            testCase.verifySubstring(output,"files=172");
            testCase.verifySubstring(output,"[AGROW, performance]");
            testCase.verifySubstring(output,"Variable appears to change size");
            testCase.verifyFalse(contains(output,testCase.repositoryRoot));
        end

        function unreachableStatementBlocks(testCase)
            fixturePath = testCase.writeUnreachableFixture("unreachableCode.m",false);
            report = analyzeProductionCode(testCase.temporaryFolder,Files=fixturePath,ShouldPrint=false,ShouldFail=false);
            testCase.verifyEqual(report.BlockingFindings.CheckID,"UNRCH");
            testCase.verifyEqual(report.BlockingFindings.RelativeFile,"unreachableCode.m");
            testCase.verifyFalse(report.BlockingFindings.Suppressed);
            testCase.verifyError(@()analyzeProductionCode(testCase.temporaryFolder,Files=fixturePath,ShouldPrint=false),"WaveVortexModel:CodeAnalyzerFailed");
        end

        function suppressedUnreachableStatementStillBlocks(testCase)
            fixturePath = testCase.writeUnreachableFixture("suppressedUnreachableCode.m",true);
            report = analyzeProductionCode(testCase.temporaryFolder,Files=fixturePath,ShouldPrint=false,ShouldFail=false);
            testCase.verifyEqual(report.BlockingFindings.CheckID,"UNRCH");
            testCase.verifyTrue(report.BlockingFindings.Suppressed);
            testCase.verifyError(@()analyzeProductionCode(testCase.temporaryFolder,Files=fixturePath,ShouldPrint=false),"WaveVortexModel:CodeAnalyzerFailed");
        end
    end

    methods (Access=private)
        function fixturePath = writeUnreachableFixture(testCase,name,shouldSuppress)
            unreachableStatement = "y = x + 1;";
            if shouldSuppress
                unreachableStatement = unreachableStatement + " %#ok<UNRCH>";
            end
            fixturePath = fullfile(testCase.temporaryFolder,name);
            writelines([
                "function y = " + erase(name,".m") + "(x)"
                "y = x;"
                "return"
                unreachableStatement
                "end"
                ],fixturePath);
        end
    end
end
