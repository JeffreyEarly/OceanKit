classdef TestThreeInterfaceBenchmark < matlab.unittest.TestCase
    properties
        RepositoryRoot (1,1) string
        TemporaryFolder (1,1) string
    end

    methods (TestClassSetup)
        function addBenchmarkPath(testCase)
            testCase.RepositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.RepositoryRoot,"Benchmarks")));
        end
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function normalizationPreservesThreeMatchedInterfaces(testCase)
            raw = rawFixture;
            rawPath = fullfile(testCase.TemporaryFolder,"three-interface-benchmark.json");
            writelines(jsonencode(raw),rawPath);
            dataset = publishedThreeInterfaceBenchmarkFromArtifact(rawPath,platformId="m5-max",platformName="Apple M5 Max");
            testCase.verifyEqual(dataset.schemaVersion,"published-three-interface-v1")
            testCase.verifyEqual(dataset.datasetId,"three-interface--m5-max--20260815T120000Z")
            testCase.verifyEqual(dataset.provider.scope,"compiled-interfaces-only")
            testCase.verifyNumElements(dataset.cases,3)
            for benchmarkCase = dataset.cases
                testCase.verifyNumElements(benchmarkCase{1}.interfaces,3)
                testCase.verifyTrue(benchmarkCase{1}.correctness.passed)
                testCase.verifyEqual(benchmarkCase{1}.interfaces{1}.providerId,"matlab-builtin")
                testCase.verifyEqual(benchmarkCase{1}.interfaces{2}.providerId,"native-neon-pthreads")
                testCase.verifyEqual(benchmarkCase{1}.interfaces{2}.integrationRatio,0.5)
                testCase.verifyEqual(benchmarkCase{1}.interfaces{3}.totalRSSRatio,0.25)
            end
        end

        function toleranceHashMatchesPortableReference(testCase)
            testCase.verifyEqual(threeInterfaceToleranceHash([1e-6 1e-5]),"11095511005457332475")
            testCase.verifyEqual(threeInterfaceToleranceHash([1e-8 1e-8]),"3713629807952050427")
            values = [1e-6 1e-8 1e-10];
            testCase.verifyEqual(threeInterfaceToleranceHash(values.*(1+eps)),threeInterfaceToleranceHash(values))
            testCase.verifyNotEqual(threeInterfaceToleranceHash(values.*(1+1e-8)),threeInterfaceToleranceHash(values))
        end

        function integratorStudyNormalizesPrimaryAndDiagnosticEvidence(testCase)
            raw = integratorStudyFixture;
            matlabRuns = find(string({raw.runs.interface})~="standalone-compiled");
            for iRun = matlabRuns
                raw.runs(iRun).integrator = rmfield(raw.runs(iRun).integrator,"sharedAbstractionStateSizedCopyCount");
            end
            rawPath = writeRaw(testCase,raw,"integrator-study.json");
            dataset = publishedThreeInterfaceBenchmarkFromArtifact(rawPath,archiveFileName="integrator-study.json.gz",archiveSHA256=repmat('d',1,64),archiveCompressedBytes=4096);
            testCase.verifyEqual(dataset.schemaVersion,"published-three-interface-v3")
            testCase.verifyEqual(dataset.studyId,"integrator-runtime-memory-v1")
            testCase.verifyEqual(dataset.provenance.rawArtifactSHA256,string(extractBefore(strtrim(systemOutput("/usr/bin/shasum -a 256 "+shellQuote(rawPath))),65)))
            testCase.verifyEqual(dataset.provenance.externalArchive.location,"external sibling archive ../wave-vortex-model-benchmark-artifacts/three-interface; not distributed with source")
            testCase.verifyNumElements(dataset.cases,2)
            endpoint = dataset.cases{1};
            dense = dataset.cases{2};
            testCase.verifyEqual(endpoint.workload,"coefficient-endpoint")
            testCase.verifyEqual(endpoint.interfaces{3}.integrationSeconds,0.25)
            testCase.verifyEqual(endpoint.interfaces{3}.totalPeakRSSBytes,2^30)
            testCase.verifyTrue(endpoint.interfaces{3}.diagnostics.integratorStorage.exact)
            testCase.verifyEqual(endpoint.interfaces{3}.diagnostics.integratorStorage.denseHistoryBytes,0)
            testCase.verifyEqual(dense.interfaces{3}.diagnostics.integratorStorage.continuousExtensionBytes,4*1024)
            testCase.verifyEqual(dense.interfaces{3}.diagnostics.methodWork.continuousExtensionRightHandSideEvaluationCounts,4)
            testCase.verifyTrue(dense.correctness.endpointTrajectoryAgreementPassed)
        end

        function integratorStudyRejectsWrongMemoryBoundary(testCase)
            raw = integratorStudyFixture;
            raw.runs(1).memory.boundary = "process-lifetime-rss";
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:IncomparableMemory")
        end

        function integratorStudyRejectsEndpointDenseStorage(testCase)
            raw = integratorStudyFixture;
            interfaces = reshape(string({raw.runs.interface}),[],1);
            workloads = reshape(string(arrayfun(@(run)run.case.workload,raw.runs,"UniformOutput",false)),[],1);
            standalone = find(interfaces=="standalone-compiled" & workloads=="coefficient-endpoint",1);
            testCase.assertNotEmpty(standalone)
            integrator = raw.runs(standalone).integrator;
            integrator.denseHistoryStateEquivalentCount = 1;
            integrator.storageAccounting.denseHistoryBytes = 1024;
            raw.runs(standalone).integrator = integrator;
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:UnexpectedDenseOutputWork")
        end

        function integratorStudyRejectsMissingRK78LazyExtension(testCase)
            raw = integratorStudyFixture;
            interfaces = reshape(string({raw.runs.interface}),[],1);
            workloads = reshape(string(arrayfun(@(run)run.case.workload,raw.runs,"UniformOutput",false)),[],1);
            standalone = find(interfaces=="standalone-compiled" & workloads=="composite-dense-output",1);
            testCase.assertNotEmpty(standalone)
            integrator = raw.runs(standalone).integrator;
            integrator.continuousExtensionWorkspaceMaximumLiveStateEquivalentCount = 0;
            raw.runs(standalone).integrator = integrator;
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:MissingRK78LazyExtension")
        end

        function integratorStudyRejectsChangedEndpointTrajectory(testCase)
            raw = integratorStudyFixture;
            raw.comparison(2).endpointTrajectoryAgreementPassed = false;
            raw.comparison(2).matchedContractPassed = false;
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:EndpointTrajectoryMismatch")
        end

        function integratorStudyRejectsComponentToleranceMismatch(testCase)
            raw = integratorStudyFixture;
            raw.comparison(1).absoluteToleranceFingerprintAgreementPassed = false;
            raw.comparison(1).matchedContractPassed = false;
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:AdaptiveToleranceMismatch")
        end

        function integratorStudyUsesMatchedMethodOutputTolerances(testCase)
            raw = integratorStudyFixture;
            tracerIndex = find(string({raw.comparison(2).outputGraph.categories.name})=="tracers",1);
            raw.comparison(2).maximumRelativeError = 1.25e-7;
            raw.comparison(2).outputGraph.maximumAbsoluteError = 1.25e-7;
            raw.comparison(2).outputGraph.maximumRelativeError = 1.25e-7;
            raw.comparison(2).outputGraph.categories(tracerIndex).maximumAbsoluteError = 1.25e-7;
            raw.comparison(2).outputGraph.categories(tracerIndex).maximumRelativeError = 1.25e-7;
            testCase.verifyWarningFree(@()validateThreeInterfaceBenchmarkContract(raw))
        end

        function storedPayloadEvidenceCanBeConservativelyReclassified(testCase)
            raw = integratorStudyFixture;
            comparison = raw.comparison(2);
            tracerIndex = find(string({comparison.outputGraph.categories.name})=="tracers",1);
            comparison.maximumRelativeError = 1.25e-7;
            comparison.outputAgreementPassed = false;
            comparison.outputGraph.passed = false;
            comparison.outputGraph.maximumAbsoluteError = 1.25e-7;
            comparison.outputGraph.maximumRelativeError = 1.25e-7;
            comparison.outputGraph.categories(tracerIndex).maximumAbsoluteError = 1.25e-7;
            comparison.outputGraph.categories(tracerIndex).maximumRelativeError = 1.25e-7;
            comparison.outputGraph.categories(tracerIndex).passed = false;
            comparison.outputGraph.differences = char("/tracers/dye payload differs");
            comparison.matchedContractPassed = false;
            [comparison,recovery] = reclassifyThreeInterfaceOutputEvidence(comparison,1e-3,1e-6);
            testCase.verifyTrue(comparison.outputGraph.passed)
            testCase.verifyTrue(comparison.matchedContractPassed)
            testCase.verifyEqual(recovery.reclassifiedCategories,"tracers")
            testCase.verifyEmpty(comparison.outputGraph.differences)
        end

        function storedStructuralDifferenceCannotBeReclassified(testCase)
            raw = integratorStudyFixture;
            comparison = raw.comparison(2);
            tracerIndex = find(string({comparison.outputGraph.categories.name})=="tracers",1);
            comparison.outputAgreementPassed = false;
            comparison.outputGraph.passed = false;
            comparison.outputGraph.categories(tracerIndex).passed = false;
            comparison.outputGraph.differences = "/tracers/dye shape differs";
            comparison.matchedContractPassed = false;
            [comparison,recovery] = reclassifyThreeInterfaceOutputEvidence(comparison,1e-3,1e-6);
            testCase.verifyFalse(comparison.outputGraph.passed)
            testCase.verifyFalse(comparison.matchedContractPassed)
            testCase.verifyFalse(recovery.passed)
            testCase.verifyEqual(comparison.outputGraph.differences,"/tracers/dye shape differs")
        end

        function rssSamplerAcknowledgesShortIntegrationPhase(testCase)
            phasePath = fullfile(testCase.TemporaryFolder,"phase.txt");
            samplePath = fullfile(testCase.TemporaryFolder,"rss.tsv");
            stdoutPath = fullfile(testCase.TemporaryFolder,"stdout.txt");
            stderrPath = fullfile(testCase.TemporaryFolder,"stderr.txt");
            workerPath = fullfile(testCase.TemporaryFolder,"worker.sh");
            writelines(["#!/bin/sh"; "printf '%s\n' integrate > ""$1"""; "while :; do"; "  if [ -f ""$WV_RSS_PHASE_ACK"" ] && [ ""$(sed -n '1p' ""$WV_RSS_PHASE_ACK"")"" = integrate ]; then exit 0; fi"; "  sleep 0.001"; "done"],workerPath);
            sampler = fullfile(testCase.RepositoryRoot,"Benchmarks","runProcessWithRSS.sh");
            command = strjoin([shellQuote(sampler) shellQuote(samplePath) shellQuote(phasePath) "0.005" shellQuote(stdoutPath) shellQuote(stderrPath) "-- /bin/sh" shellQuote(workerPath) shellQuote(phasePath)]," ");
            [status,output] = system(command);
            testCase.assertEqual(status,0,output)
            samples = string(fileread(samplePath));
            testCase.verifySubstring(samples,sprintf('\tintegrate\t'))
            testCase.verifyEqual(strtrim(string(fileread(phasePath+".sampled"))),"integrate")
        end

        function compositionPreservesFrozenCasesAndCorrectedAdaptiveEvidence(testCase)
            frozen = publishedThreeInterfaceBenchmarkFromArtifact(writeRaw(testCase,rawFixture,"frozen.json"),platformId="lyra",platformName="Apple M3 Max");
            adaptiveRaw = rawFixture;
            adaptiveRaw.runId = "20260816T120000000Z";
            adaptiveRaw.cases = adaptiveRaw.cases(3);
            adaptiveRaw.comparison = adaptiveRaw.comparison(3);
            adaptiveRaw.runs = adaptiveRaw.runs(arrayfun(@(run)string(run.case.id)=="adaptive-rk23-observer-output",adaptiveRaw.runs));
            adaptive = publishedThreeInterfaceBenchmarkFromArtifact(writeRaw(testCase,adaptiveRaw,"adaptive.json"),platformId="lyra",platformName="Apple M3 Max");
            adaptive.provider.moduleSHA256 = repmat('c',1,64);
            frozenPath = writeRaw(testCase,frozen,"frozen-published.json");
            adaptivePath = writeRaw(testCase,adaptive,"adaptive-published.json");
            composed = composePublishedThreeInterfaceBenchmark(frozenPath,adaptivePath);
            testCase.verifyEqual(composed.schemaVersion,"published-three-interface-v2")
            testCase.verifyEqual(arrayfun(@(index)string(composed.cases{index}.id),1:numel(composed.cases)),["nonlinear-flux" "fixed-rk4-continuation" "adaptive-rk23-observer-output"])
            testCase.verifyEqual(composed.cases{1}.evidence.datasetId,frozen.datasetId)
            testCase.verifyEqual(composed.cases{3}.evidence.datasetId,adaptive.datasetId)
            testCase.verifyEqual(string(composed.cases{1}.evidence.provider.moduleSHA256),string(frozen.provider.moduleSHA256))
            testCase.verifyEqual(string(composed.cases{3}.evidence.provider.moduleSHA256),string(adaptive.provider.moduleSHA256))
            testCase.verifyEqual(composed.provider.moduleSHA256,"per-case-evidence")
            testCase.verifyEqual(composed.provider.moduleIdentityScope,"case-evidence")
            testCase.verifyEqual(composed.provenance.composition,"frozen-valid-v1-plus-corrected-adaptive")
        end

        function compositionRejectsIncompatibleProviderConfiguration(testCase)
            frozen = publishedThreeInterfaceBenchmarkFromArtifact(writeRaw(testCase,rawFixture,"frozen.json"),platformId="lyra",platformName="Apple M3 Max");
            adaptiveRaw = rawFixture;
            adaptiveRaw.runId = "20260816T120000000Z";
            adaptiveRaw.cases = adaptiveRaw.cases(3);
            adaptiveRaw.comparison = adaptiveRaw.comparison(3);
            adaptiveRaw.runs = adaptiveRaw.runs(arrayfun(@(run)string(run.case.id)=="adaptive-rk23-observer-output",adaptiveRaw.runs));
            adaptive = publishedThreeInterfaceBenchmarkFromArtifact(writeRaw(testCase,adaptiveRaw,"adaptive.json"),platformId="lyra",platformName="Apple M3 Max");
            adaptive.provider.threadBackend = "openmp";
            frozenPath = writeRaw(testCase,frozen,"frozen-published.json");
            adaptivePath = writeRaw(testCase,adaptive,"adaptive-published.json");
            testCase.verifyError(@()composePublishedThreeInterfaceBenchmark(frozenPath,adaptivePath),"WaveVortexBenchmark:IncompatiblePublishedEvidence")
        end

        function dirtyOrIncompleteArtifactsCannotBePublished(testCase)
            raw = rawFixture;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            raw.source.isDirty = true;
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:InvalidThreeInterfaceArtifact")
            raw.source.isDirty = false; raw.status = "failed";
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:InvalidThreeInterfaceArtifact")
        end

        function integratorMismatchCannotBePublished(testCase)
            raw = rawFixture;
            raw.runs(2).integrator.actual = "adaptive-rk23";
            raw.runs(2).integrator.matched = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:IntegratorMismatch")
        end

        function adaptiveWorkMismatchCannotBePublished(testCase)
            raw = rawFixture;
            adaptiveRun = find(arrayfun(@(run)string(run.case.id)=="adaptive-rk23-observer-output",raw.runs),1);
            raw.runs(adaptiveRun).integrator.rhsEvaluationCount = raw.runs(adaptiveRun).integrator.rhsEvaluationCount+1;
            raw.comparison(3).adaptiveWorkAgreementPassed = false;
            raw.comparison(3).matchedContractPassed = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:AdaptiveWorkMismatch")
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:AdaptiveWorkMismatch")
        end

        function toleranceFingerprintMismatchIsDiagnostic(testCase)
            raw = rawFixture;
            adaptiveRuns = find(arrayfun(@(run)string(run.case.id)=="adaptive-rk23-observer-output",raw.runs));
            raw.runs(adaptiveRuns(end)).integrator.absoluteToleranceHash = "different";
            raw.runs(adaptiveRuns(end)).integrator.absoluteToleranceComponentHashes(1) = "different";
            raw.comparison(3).absoluteToleranceFingerprintAgreementPassed = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            validateThreeInterfaceBenchmarkContract(raw)
            publishedThreeInterfaceBenchmarkFromArtifact(rawPath);
        end

        function matchingLegacyFingerprintIsRecovered(testCase)
            raw = rawFixture;
            raw.comparison = rmfield(raw.comparison,"absoluteToleranceFingerprintAgreementPassed");
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            validateThreeInterfaceBenchmarkContract(raw)
            published = publishedThreeInterfaceBenchmarkFromArtifact(rawPath);
            testCase.verifyTrue(published.cases{3}.adaptiveWork.absoluteToleranceFingerprintAgreementPassed)
        end

        function mixedOutputSchedulesAreRejected(testCase)
            if ~isCanonicalNativePlatform
                return
            end
            testCase.verifyError(@()runThreeInterfaceBenchmark(caseIds=["fixed-rk4-continuation" "adaptive-rk23-observer-output"],shouldWriteArtifacts=false),"WaveVortexBenchmark:MixedOutputSchedules")
        end

        function incomparableMemoryCannotBePublished(testCase)
            raw = rawFixture;
            raw.runs(1).memory.status = "failed";
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:IncomparableMemory")
        end

        function outputGraphMismatchCannotBePublished(testCase)
            raw = rawFixture;
            raw.comparison(2).outputGraph.categories(4).passed = false;
            raw.comparison(2).outputGraph.passed = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:OutputGraphMismatch")
        end

        function failedMatchedContractCannotBePublished(testCase)
            raw = rawFixture;
            raw.comparison(1).matchedContractPassed = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:MatchedContractFailed")
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:MatchedContractFailed")
        end

        function finiteNonlinearFluxErrorAboveToleranceCannotBePublished(testCase)
            raw = rawFixture;
            raw.comparison(1).maximumRelativeError = 2e-12;
            raw.comparison(1).outputGraph.maximumRelativeError = 2e-12;
            raw.comparison(1).outputGraph.categories.maximumRelativeError = 2e-12;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()validateThreeInterfaceBenchmarkContract(raw),"WaveVortexBenchmark:NumericalMismatch")
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:NumericalMismatch")
        end

        function outputAgreementFlagCannotBePublished(testCase)
            raw = rawFixture;
            raw.comparison(2).outputAgreementPassed = false;
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:OutputGraphMismatch")
        end

        function providerMismatchCannotBePublished(testCase)
            raw = rawFixture;
            raw.runs(2).provider.id = "matlab-bundled";
            rawPath = fullfile(testCase.TemporaryFolder,"raw.json");
            writelines(jsonencode(raw),rawPath);
            testCase.verifyError(@()publishedThreeInterfaceBenchmarkFromArtifact(rawPath),"WaveVortexBenchmark:ProviderMismatch")
        end

        function benchmarkWorkerRemainsAuthorOnly(testCase)
            manifest = string(fileread(fullfile(testCase.RepositoryRoot,"resources","mpackage.json")));
            testCase.verifyFalse(contains(manifest,"Benchmarks"))
            cmake = string(fileread(fullfile(testCase.RepositoryRoot,"PortableRuntime","CMakeLists.txt")));
            testCase.verifySubstring(cmake,"WV_RUNTIME_BUILD_BENCHMARKS")
            testCase.verifySubstring(cmake,"wv-standalone-nonlinear-flux-benchmark")
        end

        function verboseThreeInterfaceArtifactsAreNotTracked(testCase)
            [status,output] = system("git -C "+shellQuote(testCase.RepositoryRoot)+" ls-files");
            testCase.assertEqual(status,0)
            tracked = splitlines(strtrim(string(output)));
            forbidden = endsWith(tracked,"/three-interface-benchmark.json") | endsWith(tracked,"-rss.tsv") | endsWith(tracked,"worker.json") | endsWith(tracked,"three-interface-benchmark.json.gz") | startsWith(tracked,"docs/benchmarks/raw/three-interface--");
            testCase.verifyFalse(any(forbidden),"Verbose three-interface results must remain in the external compressed archive.")
            compact = tracked(startsWith(tracked,"Benchmarks/results/published/three-interface--") | startsWith(tracked,"docs/benchmarks/data/three-interface--"));
            for path = reshape(compact,1,[])
                if ~isfile(fullfile(testCase.RepositoryRoot,path))
                    continue
                end
                information = dir(fullfile(testCase.RepositoryRoot,path));
                testCase.verifyLessThanOrEqual(information.bytes,512*1024,"Published three-interface records must remain compact.")
            end
        end

        function completeOutputGraphRejectsPayloadAndMetadataDifferences(testCase)
            reference = fullfile(testCase.TemporaryFolder,"reference.nc");
            candidate = fullfile(testCase.TemporaryFolder,"candidate.nc");
            createOutputGraphFixture(reference,false,false);
            copyfile(reference,candidate)
            ncwrite(candidate,"/particles/float_x",[1; 9]);
            comparison = compareWaveVortexOutputGraphs(reference,candidate);
            testCase.verifyFalse(comparison.passed)
            testCase.verifyTrue(any(contains(comparison.differences,"float_x payload")))
            copyfile(reference,candidate,"f")
            ncwriteatt(candidate,"/particles","outputInterval",2);
            comparison = compareWaveVortexOutputGraphs(reference,candidate);
            testCase.verifyFalse(comparison.passed)
            testCase.verifyTrue(any(contains(comparison.differences,"outputInterval")))
        end

        function completeOutputGraphRejectsShapeOrderAndExtraOutput(testCase)
            reference = fullfile(testCase.TemporaryFolder,"reference.nc");
            reordered = fullfile(testCase.TemporaryFolder,"reordered.nc");
            wrongShape = fullfile(testCase.TemporaryFolder,"wrong-shape.nc");
            createOutputGraphFixture(reference,false,false);
            createOutputGraphFixture(reordered,true,false);
            createOutputGraphFixture(wrongShape,false,true);
            testCase.verifyFalse(compareWaveVortexOutputGraphs(reference,reordered).passed)
            testCase.verifyFalse(compareWaveVortexOutputGraphs(reference,wrongShape).passed)
            copyfile(reference,reordered,"f")
            nccreate(reordered,"/particles/extra",Dimensions={"id",2});
            ncwrite(reordered,"/particles/extra",[1; 2]);
            testCase.verifyFalse(compareWaveVortexOutputGraphs(reference,reordered).passed)
        end

        function completeOutputGraphUsesBoundedMultidimensionalChunks(testCase)
            reference = fullfile(testCase.TemporaryFolder,"reference.nc");
            candidate = fullfile(testCase.TemporaryFolder,"candidate.nc");
            createOutputGraphFixture(reference,false,false);
            copyfile(reference,candidate)
            comparison = compareWaveVortexOutputGraphs(reference,candidate,maximumChunkBytes=16);
            testCase.verifyTrue(comparison.passed)
            testCase.verifyEqual(comparison.maximumAbsoluteError,0)
        end
    end

    methods (Test,TestTags="optional")
        function reducedDenseScheduleUsesAnchoredPhysicalContract(testCase)
            if ~isCanonicalNativePlatform
                return
            end
            result = runThreeInterfaceBenchmark(Nxyz=[8 6 5],processRunCount=1,deltaT=1/1024,adaptiveInitialStep=1/1024,integrators="fixed-rk4",workloads="composite-dense-output",physicalConfigurations="nonhydrostatic",integrationStepCount=20,samplingIntervalSeconds=0.005,plateauSeconds=0.2,shouldWriteArtifacts=false);
            testCase.verifyEqual(result.status,"complete")
            testCase.verifyEqual(result.cases.deltaT,1/1024)
            testCase.verifyEqual(result.cases.denseOutputPointsPerStep,3)
            testCase.verifyEqual(result.cases.denseOutputStartTime,0)
            testCase.verifyEqual(result.cases.denseOutputEndTime,1/1024)
            testCase.verifyEqual(result.cases.denseOutputRecordCount,5)
            testCase.verifyEqual(result.cases.denseOutputIntegrationRecordCount,4)
            testCase.verifyTrue(result.comparison.matchedContractPassed)
            testCase.verifyNumElements(result.repeatComparisonEvidence,1)
            testCase.verifyEqual(result.repeatComparisonEvidence.releasedRunCount,3)
            testCase.verifyGreaterThan(result.repeatComparisonEvidence.releasedBytes,0)
            testCase.verifyTrue(all(arrayfun(@(run)string(run.output.retention)=="released-after-repeat-correctness",result.runs)))
        end

        function reducedMatchedBenchmarkRunsAllInterfaces(testCase)
            if ~isCanonicalNativePlatform
                testCase.verifyError(@()runThreeInterfaceBenchmark(shouldWriteArtifacts=false),"WaveVortexBenchmark:ThreeInterfaceUnsupportedPlatform")
                return
            end
            outputDirectory = fullfile(testCase.TemporaryFolder,"results");
            archiveDirectory = fullfile(testCase.TemporaryFolder,"archive");
            result = runThreeInterfaceBenchmarkComparison(resolutions=[8 6 5],processRunCount=1,deltaT=1e-4,caseIds="adaptive-rk23-observer-output",adaptiveStepCount=10,adaptiveOutputCount=2,samplingIntervalSeconds=0.005,plateauSeconds=0.1,outputRoot=outputDirectory,archiveDirectory=archiveDirectory);
            testCase.verifySize(result,[1 1])
            testCase.verifyEqual(result.status,"complete")
            testCase.verifyEqual(string({result.comparison.id}),"adaptive-rk23-observer-output")
            testCase.verifyTrue(all([result.comparison.outputAgreementPassed]))
            testCase.verifyLessThanOrEqual(max([result.comparison.maximumRelativeError]),1e-12)
            testCase.verifyTrue(all([result.comparison.integratorAgreementPassed]))
            testCase.verifyTrue(all([result.comparison.memoryAgreementPassed]))
            testCase.verifyTrue(all(arrayfun(@(item)item.matchedContractPassed,result.comparison)))
            testCase.verifyTrue(isfile(fullfile(archiveDirectory,result.externalArchive.fileName)))
            testCase.verifyEqual(strlength(result.externalArchive.sha256),64)
            testCase.verifyGreaterThan(result.externalArchive.compressedBytes,0)
        end
    end
end

function createOutputGraphFixture(pathname,reordered,wrongShape)
if reordered
    variableNames = ["u" "t"];
else
    variableNames = ["t" "u"];
end
for name = variableNames
    if name == "t"
        nccreate(pathname,"/wave-vortex/t",Dimensions={"t",2});
        ncwrite(pathname,"/wave-vortex/t",[0; 1]);
    else
        nccreate(pathname,"/wave-vortex/u",Dimensions={"x",2,"t",2});
        ncwrite(pathname,"/wave-vortex/u",reshape(1:4,2,2));
    end
end
ncwriteatt(pathname,"/wave-vortex","outputInterval",1);
nccreate(pathname,"/wave-vortex/Ap_real",Dimensions={"kl",2,"j",2,"t",2});
ncwrite(pathname,"/wave-vortex/Ap_real",reshape(1:8,2,2,2));
nccreate(pathname,"/wave-vortex/mooring_u",Dimensions={"id",2,"t",2});
ncwrite(pathname,"/wave-vortex/mooring_u",reshape(1:4,2,2));
nccreate(pathname,"/particles/t",Dimensions={"t",2});
ncwrite(pathname,"/particles/t",[0; 1]);
particleCount = conditional(wrongShape,3,2);
nccreate(pathname,"/particles/float_x",Dimensions={"id",particleCount});
ncwrite(pathname,"/particles/float_x",(1:particleCount)');
ncwriteatt(pathname,"/particles/float_x","isParticle","1");
ncwriteatt(pathname,"/particles","outputInterval",1);
nccreate(pathname,"/tracers/t",Dimensions={"t",2});
ncwrite(pathname,"/tracers/t",[0; 1]);
nccreate(pathname,"/tracers/dye",Dimensions={"x",2,"y",2,"t",2});
ncwrite(pathname,"/tracers/dye",reshape(1:8,2,2,2));
ncwriteatt(pathname,"/tracers/dye","isTracer","1");
end

function value = isCanonicalNativePlatform
value = ismac && computer("arch") == "maca64";
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end


function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function pathname = writeRaw(testCase,value,name)
pathname = fullfile(testCase.TemporaryFolder,name);
writelines(jsonencode(value),pathname);
end

function raw = rawFixture
interfaces = [interfaceRecord("matlab-builtin",1,1,1,1); interfaceRecord("matlab-compiled",0.5,0.5,2,2); interfaceRecord("standalone-compiled",0.25,0.25,0.25,0.25)];
definitions = [caseDefinition("nonlinear-flux","nonlinearFlux","none"); caseDefinition("fixed-rk4-continuation","model-continuation","fixed-rk4"); caseDefinition("adaptive-rk23-observer-output","model-continuation","adaptive-rk23")];
comparison = repmat(struct("id","","interfaces",interfaces,"maximumRelativeError",1e-14,"outputAgreementPassed",true,"outputGraph",modelOutputGraph,"integratorAgreementPassed",true,"adaptiveWorkAgreementPassed",true,"absoluteToleranceFingerprintAgreementPassed",true,"memoryAgreementPassed",true,"matchedContractPassed",true),3,1);
runs = repmat(runRecord("matlab-builtin",definitions(1)),0,1);
for iCase = 1:3
    comparison(iCase).id = definitions(iCase).id;
    if definitions(iCase).operation == "nonlinearFlux"
        comparison(iCase).outputGraph = fluxOutputGraph;
    end
    for identifier = ["matlab-builtin" "matlab-compiled" "standalone-compiled"]
        runs(end+1,1) = runRecord(identifier,definitions(iCase)); %#ok<AGROW>
    end
end
provider = struct("status","available","isAvailable",true,"provider",struct("id","native-neon-pthreads","version","3.3.11","threadBackend","pthreads"),"module",struct("sha256",repmat('b',1,64),"identityValidated",true),"libraries",struct("base",struct("path","/tmp/libfftw3.3.dylib"),"thread",struct("path","/tmp/libfftw3_threads.3.dylib"),"openmp",struct("detected",false)),"contract",struct("threadCount",18),"featureValidation",struct("maximumRelativeError",1e-14));
environment = struct("processor","Apple M5 Max","physicalMemoryBytes",64*2^30,"os","macOS","architecture","maca64","matlabVersion","R2026a Update 4");
source = struct("commit",repmat('a',1,40),"tree",repmat('c',1,40),"isDirty",false);
configuration = struct("Lxyz",[150e3 150e3 1300],"processRunCount",1,"warmupCount",0,"samplesPerProcess",1,"threadCount",18,"correctnessTolerance",1e-12);
raw = struct("schemaVersion","three-interface-benchmark-v1","status","complete","runId","20260815T120000000Z","source",source,"environment",environment,"configuration",configuration,"provider",provider,"cases",definitions,"runs",runs,"comparison",comparison);
end

function raw = integratorStudyFixture
raw = rawFixture;
raw.schemaVersion = "three-interface-benchmark-v2";
raw.configuration.studyId = "integrator-runtime-memory-v1";
raw.configuration.fixtures = [struct("physicalConfiguration","nonhydrostatic","workload","coefficient-endpoint","sha256",repmat('1',1,64)); struct("physicalConfiguration","nonhydrostatic","workload","composite-dense-output","sha256",repmat('2',1,64))];
raw.configuration.initialCondition = struct("id","gm1-red-geostrophic-j1-v1","gmEnergyLevel",1,"geostrophicMaximumSpeedTarget",0.15);
raw.configuration.stepControls = struct("fixedStep",128,"adaptiveInitialStep",295.7935799274,"adaptiveMaximumStepPolicy","MATLAB default; no user maximum step","finalTime",7168);
raw.cases = [integratorStudyDefinition("coefficient-endpoint"); integratorStudyDefinition("composite-dense-output")];
interfaces = [interfaceRecord("matlab-builtin",1,1,1,1); interfaceRecord("matlab-compiled",0.5,0.5,1,1); interfaceRecord("standalone-compiled",0.25,0.25,1,1)];
raw.comparison = repmat(struct("id","","interfaces",interfaces,"maximumRelativeError",1e-14,"outputAgreementPassed",true,"outputGraph",modelOutputGraph,"integratorAgreementPassed",true,"adaptiveWorkAgreementPassed",true,"absoluteToleranceFingerprintAgreementPassed",true,"memoryAgreementPassed",true,"endpointTrajectoryAgreementPassed",true,"matchedContractPassed",true),2,1);
raw.runs = repmat(integratorStudyRun("matlab-builtin",raw.cases(1)),0,1);
for iCase = 1:2
    raw.comparison(iCase).id = raw.cases(iCase).id;
    if raw.cases(iCase).workload == "coefficient-endpoint"
        raw.comparison(iCase).outputGraph = endpointOutputGraph;
    end
    for interface = ["matlab-builtin" "matlab-compiled" "standalone-compiled"]
        raw.runs(end+1,1) = integratorStudyRun(interface,raw.cases(iCase)); %#ok<AGROW>
    end
end
end

function value = integratorStudyDefinition(workload)
value = caseDefinition("nonhydrostatic--adaptive-rk78--"+workload,"model-continuation","adaptive-rk78");
value.physicalConfiguration = "nonhydrostatic";
value.isHydrostatic = false;
value.workload = workload;
value.deltaT = 128;
value.integrationStepCount = 56;
value.denseOutputPointsPerStep = 3;
value.finalTime = 7168;
value.initialStep = 295.7935799274;
value.maximumStep = [];
value.maximumStepPolicy = "matlab-default";
value.outputRelativeTolerance = value.relativeTolerance;
value.outputAbsoluteTolerance = value.absoluteTolerance;
value.outputInterval = conditional(workload=="coefficient-endpoint",7168,32);
value.denseOutputStartTime = conditional(workload=="coefficient-endpoint",[],0);
value.denseOutputEndTime = conditional(workload=="coefficient-endpoint",[],128);
value.denseOutputRecordCount = conditional(workload=="coefficient-endpoint",0,5);
value.denseOutputIntegrationRecordCount = conditional(workload=="coefficient-endpoint",0,4);
value.observerGraph = conditional(workload=="coefficient-endpoint","coefficient-only state with endpoint-only output","one persisted restart record plus four first-step field, particle, tracer, and mooring records");
end

function value = integratorStudyRun(identifier,definition)
value = runRecord(identifier,caseDefinition("adaptive-rk23-observer-output","model-continuation","adaptive-rk23"));
value.case = definition;
value.interface = identifier;
value.integrationSeconds = conditional(identifier=="matlab-builtin",1,conditional(identifier=="matlab-compiled",0.5,0.25));
value.memory.boundary = "integration-phase-total-live-process-tree-rss";
value.memory.processLifetimePeakRSSBytes = 2^30;
value.memory.integrationSampleCount = 4;
value.memory.baselineProcessBytes = 2^28;
value.provider = runRecord(identifier,definition).provider;
isDense = definition.workload == "composite-dense-output";
integrator = struct("requested","adaptive-rk78","actual","adaptive-rk78","matched",true,"controller","matlab-ode78-v1","relativeTolerance",definition.relativeTolerance,"absoluteToleranceHash","12345","absoluteToleranceHashClearedMantissaBits",20,"absoluteToleranceComponentHashes",repmat("123",1,7),"requestedInitialStep",definition.initialStep,"effectiveInitialStep",definition.initialStep,"requestedMaximumStep",NaN,"maximumStepPolicy",definition.maximumStepPolicy,"effectiveMaximumStep",0.1*definition.finalTime,"initialTime",0,"finalTime",definition.finalTime,"acceptedStepCount",14,"rejectedStepCount",0,"rhsEvaluationCount",182,"denseOutputEvaluationCount",conditional(isDense,4,0),"fsalReuseCount",0,"fsalInvalidationCount",0,"outputRecordCounts",struct("waveVortex",2,"dense",conditional(isDense,5,0),"particles",conditional(isDense,5,0),"tracers",conditional(isDense,5,0)),"continuousExtensionRightHandSideEvaluationCount",conditional(isDense,4,0),"continuousExtensionWorkspaceStateEquivalentCount",0,"continuousExtensionWorkspaceMaximumLiveStateEquivalentCount",conditional(isDense,4,0),"workspaceStateEquivalentCount",conditional(isDense,15,11),"workspaceMaximumLiveStateEquivalentCount",conditional(isDense,15,11),"denseHistoryStateEquivalentCount",conditional(isDense,8,0),"sharedAbstractionStateSizedCopyCount",0);
integrator.requestedMaximumStep = definition.maximumStep;
if identifier == "standalone-compiled"
    count = integrator.workspaceStateEquivalentCount;
    integrator.stateSizedBuffers = arrayfun(@(index)struct("buffer","k"+index,"producer","stage "+index,"lastUse","accepted-step commit"),1:count);
    integrator.storageAccounting = struct("scope","integrator-owned exact allocation capacity","exact",true,"persistentBytes",2048,"workspaceCapacityBytes",count*1024,"workspaceMaximumLiveBytes",count*1024,"denseHistoryBytes",conditional(isDense,8*1024,0),"continuousExtensionBytes",conditional(isDense,4*1024,0),"sharedAbstractionStateSizedCopyCount",0,"nonStateBytes",256,"byteLedgerAgreement",true);
else
    integrator.storageAccounting = struct("scope","MATLAB solver and allocator storage is opaque","exact",false,"persistentBytes",NaN,"workspaceCapacityBytes",NaN,"workspaceMaximumLiveBytes",NaN,"denseHistoryBytes",NaN,"continuousExtensionBytes",NaN,"nonStateBytes",NaN,"byteLedgerAgreement",false);
    integrator.stateSizedBuffers = struct([]);
end
value.integrator = integrator;
end

function value = endpointOutputGraph
names = ["coefficients" "times"];
categories = arrayfun(@(name)struct("name",name,"variableCount",1,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"passed",true),names);
value = struct("kind","complete-netcdf-output-graph","passed",true,"maximumRelativeError",1e-14,"maximumAbsoluteError",0,"variableCount",2,"recordCount",2,"categories",categories,"differences",strings(0,1));
end

function value = modelOutputGraph
names = ["coefficients" "eulerianFields" "moorings" "particles" "tracers" "times"];
categories = arrayfun(@(name)struct("name",name,"variableCount",1,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"passed",true),names);
value = struct("kind","complete-netcdf-output-graph","passed",true,"maximumRelativeError",1e-14,"maximumAbsoluteError",0,"variableCount",6,"recordCount",12,"categories",categories,"differences",strings(0,1));
end

function value = fluxOutputGraph
category = struct("name","coefficients","variableCount",3,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"passed",true);
value = struct("kind","flux-arrays","passed",true,"maximumRelativeError",1e-14,"maximumAbsoluteError",0,"variableCount",3,"recordCount",0,"categories",category,"differences",strings(0,1));
end

function value = caseDefinition(identifier,operation,integrator)
value = struct("id",identifier,"operation",operation,"requestedIntegrator",integrator,"Nxyz",[256 256 129],"forcing","default WVNonlinearAdvection","shouldAntialias",true,"deltaT",1e-3,"finalTime",2e-3,"relativeTolerance",1e-3,"absoluteTolerance",1e-6,"initialStep",1e-3,"maximumStep",1e-3,"outputInterval",5e-4,"observerGraph","fields, particles, tracers");
end

function value = interfaceRecord(identifier,processRatio,integrationRatio,totalRatio,incrementRatio)
value = struct("id",identifier,"processWallSeconds",processRatio,"interfaceTotalSeconds",processRatio,"integrationSeconds",integrationRatio,"totalPeakRSSBytes",totalRatio*2^30,"incrementalPeakRSSBytes",incrementRatio*2^28,"finalRSSBytes",totalRatio*2^29,"processWallRatio",processRatio,"integrationRatio",integrationRatio,"totalRSSRatio",totalRatio,"incrementalRSSRatio",incrementRatio);
end

function value = runRecord(identifier,definition)
if identifier == "matlab-builtin"
    provider = struct("id","matlab-builtin","version","R2026a Update 4","threads",18,"baseLibrary","","threadLibrary","","noFallback",true);
else
    provider = struct("id","native-neon-pthreads","version","3.3.11","threads",18,"baseLibrary","/tmp/libfftw3.3.dylib","threadLibrary","/tmp/libfftw3_threads.3.dylib","noFallback",true);
end
memory = struct("status","complete","provider","macos-ps-process-tree","totalPeakRSSBytes",2^30,"peakIncrementBytes",2^28,"finalRSSBytes",2^29);
integrator = struct("requested",definition.requestedIntegrator,"actual",definition.requestedIntegrator,"matched",true);
if definition.requestedIntegrator == "adaptive-rk23"
    integrator.controller = "matlab-ode23-v1";
    integrator.relativeTolerance = definition.relativeTolerance;
    integrator.absoluteToleranceHash = "12345";
    integrator.absoluteToleranceHashClearedMantissaBits = 20;
    integrator.absoluteToleranceComponentHashes = repmat("123",1,7);
    integrator.requestedInitialStep = definition.initialStep;
    integrator.effectiveInitialStep = definition.initialStep;
    integrator.requestedMaximumStep = definition.maximumStep;
    integrator.effectiveMaximumStep = definition.maximumStep;
    integrator.initialTime = definition.deltaT;
    integrator.finalTime = definition.finalTime;
    integrator.acceptedStepCount = 1;
    integrator.rejectedStepCount = 0;
    integrator.rhsEvaluationCount = 4;
    integrator.denseOutputEvaluationCount = 0;
    integrator.outputRecordCounts = struct("waveVortex",3,"particles",3,"tracers",3);
    integrator.fsalReuseCount = "not-exposed-by-matlab";
    integrator.fsalInvalidationCount = "not-exposed-by-matlab";
end
value = struct("schemaVersion","three-interface-worker-v1","status","complete","interface",identifier,"case",definition,"processWallSeconds",1,"integrationSeconds",1,"memory",memory,"provider",provider,"integrator",integrator);
end

function output = systemOutput(command)
[status,output] = system(command);
if status ~= 0
    error("WaveVortexBenchmark:TestCommandFailed","%s",output)
end
output = string(output);
end
