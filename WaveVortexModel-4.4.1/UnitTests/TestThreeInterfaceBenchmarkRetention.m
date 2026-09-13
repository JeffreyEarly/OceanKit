classdef TestThreeInterfaceBenchmarkRetention < matlab.unittest.TestCase
    properties
        TemporaryFolder (1,1) string
    end
    methods (TestMethodSetup)
        function prepare(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"Benchmarks")));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryFolder = string(fixture.Folder);
        end
    end
    methods (Test,TestTags="full")
        function completedWorkersWithFailedComparisonKeepEveryPayload(testCase)
            [runs,comparison] = fixture(testCase.TemporaryFolder);
            ncwrite(runs(3).output.path,"u",2);
            outputComparison = compareWaveVortexOutputGraphs(runs(1).output.path,runs(3).output.path);
            testCase.assertFalse(outputComparison.passed);
            comparison.outputAgreementPassed = outputComparison.passed;
            testCase.verifyError(@()releaseThreeInterfaceBenchmarkOutputs(runs,comparison,testCase.TemporaryFolder),"WaveVortexBenchmark:ThreeInterfaceRepeatComparison");
            for run = runs
                testCase.verifyEqual(ncread(run.output.path,"u"),1+double(run.interface=="standalone-compiled"));
            end
        end

        function failedControlAndTrajectoryGatesKeepPayloads(testCase)
            [runs,comparison] = fixture(testCase.TemporaryFolder);
            for field = ["matchedContractPassed" "endpointTrajectoryAgreementPassed"]
                failed = comparison;
                failed.(field) = false;
                testCase.verifyError(@()releaseThreeInterfaceBenchmarkOutputs(runs,failed,testCase.TemporaryFolder),"WaveVortexBenchmark:ThreeInterfaceRepeatComparison");
                testCase.verifyTrue(all(arrayfun(@(run)isfile(run.output.path),runs)));
            end
        end

        function passingRepeatReleasesOnlyValidatedOutputs(testCase)
            [runs,comparison] = fixture(testCase.TemporaryFolder);
            configPath = fullfile(testCase.TemporaryFolder,"config.json");
            writelines("worker configuration",configPath);
            expectedBytes = sum(arrayfun(@(run)dir(run.output.path).bytes,runs));
            [released,bytes] = releaseThreeInterfaceBenchmarkOutputs(runs,comparison,testCase.TemporaryFolder);
            testCase.verifyEqual(bytes,expectedBytes);
            testCase.verifyFalse(any(arrayfun(@(run)isfile(run.output.path),runs)));
            testCase.verifyTrue(isfile(configPath));
            testCase.verifyEqual(string(arrayfun(@(run)run.output.retention,released,UniformOutput=false)),repmat("released-after-repeat-correctness",1,3));
        end

        function workerLogsRetainWarningsAndBoundLargeOutput(testCase)
            stdout = fullfile(testCase.TemporaryFolder,"stdout.txt");
            stderr = fullfile(testCase.TemporaryFolder,"stderr.txt");
            writelines("start "+string(repmat('x',1,100000))+" end",stdout);
            writelines("Warning: worker teardown diagnostic",stderr);
            diagnostics = readThreeInterfaceWorkerDiagnostics(stdout,stderr);
            testCase.verifyTrue(diagnostics.stdout.truncated);
            testCase.verifyGreaterThan(diagnostics.stdout.originalBytes,64*1024);
            testCase.verifyLessThan(strlength(diagnostics.stdout.text),64*1024+100);
            testCase.verifyTrue(startsWith(diagnostics.stdout.text,"start "));
            testCase.verifyTrue(contains(diagnostics.stdout.text," end"));
            testCase.verifySubstring(diagnostics.stderr.text,"Warning: worker teardown diagnostic");
            testCase.verifyFalse(diagnostics.stderr.truncated);
            testCase.verifyTrue(isfile(stdout));
        end
    end
    methods (Test,TestTags="optional")
        function failedWorkerRetainsEvidenceAndRestoresState(testCase)
            testCase.assumeTrue(ismac && string(computer("arch"))=="maca64" && ~isMATLABReleaseOlderThan("R2025b"));
            originalDirectory = pwd;
            originalPath = path;
            originalRng = rng;
            archive = fullfile(testCase.TemporaryFolder,"archive");
            testCase.verifyError(@()runThreeInterfaceBenchmark(Nxyz=[8 6 5],processRunCount=1,caseIds="nonlinear-flux",archiveDirectory=archive,shouldWriteArtifacts=false,injectWorkerFailure=true),"WaveVortexBenchmark:ThreeInterfaceWorker");
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(path,originalPath);
            testCase.verifyEqual(rng,originalRng);
            receipts = dir(fullfile(archive,"failures","*","three-interface-benchmark.json"));
            testCase.assertNumElements(receipts,1);
            receipt = jsondecode(fileread(fullfile(receipts.folder,receipts.name)));
            testCase.verifyEqual(string(receipt.status),"failed");
            testCase.verifyEqual(string(receipt.failure.stage),"workers");
            testCase.verifyTrue(isfolder(receipt.failureArtifacts.directory));
            testCase.verifyNotEmpty(receipt.source.commit);
            testCase.verifyTrue(receipt.provider.isAvailable);
            for pattern = ["config.json" "command.txt" "stdout.txt" "stderr.txt" "*.nc"]
                testCase.verifyNotEmpty(dir(fullfile(receipts.folder,"**",pattern)),"Missing retained evidence: "+pattern);
            end
        end
    end
end

function [runs,comparison] = fixture(folder)
interfaces = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
runs = repmat(struct("status","complete","interface","","output",struct("path","")),1,3);
for iRun = 1:3
    runs(iRun).interface = interfaces(iRun);
    runs(iRun).output.path = fullfile(folder,interfaces(iRun)+".nc");
    nccreate(runs(iRun).output.path,"u",Dimensions={"x",1});
    ncwrite(runs(iRun).output.path,"u",1);
end
comparison = struct("matchedContractPassed",true,"outputAgreementPassed",true,"endpointTrajectoryAgreementPassed",true);
end
