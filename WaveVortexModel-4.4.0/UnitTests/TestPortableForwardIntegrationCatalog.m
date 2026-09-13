classdef TestPortableForwardIntegrationCatalog < matlab.unittest.TestCase
    properties (SetAccess=private)
        repositoryRoot (1,1) string
        catalog (1,1) struct
    end
    methods (TestClassSetup)
        function loadCatalog(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.repositoryRoot,"tools")));
            testCase.catalog = jsondecode(fileread(fullfile(testCase.repositoryRoot,"PortableRuntime","contracts","portable-forward-integration-v1.json")));
        end
    end
    methods (Test,TestTags="full")
        function regenerationMatchesCommittedArtifacts(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            result = generatePortableForwardIntegrationCatalog(outputRoot=string(fixture.Folder));
            testCase.verifyEqual(fileread(result.catalogPath),fileread(fullfile(testCase.repositoryRoot,"PortableRuntime","contracts","portable-forward-integration-v1.json")));
            testCase.verifyEqual(fileread(result.documentationPath),fileread(fullfile(testCase.repositoryRoot,"PortableRuntime","INTEGRATION.md")));
        end
        function committedReceiptsEstablishExecutedQualification(testCase)
            for provider = ["reference","native-fftw"]
                path = fullfile(testCase.repositoryRoot,"PortableRuntime","qualification","forward-integration-"+provider+"-v1.json");
                testCase.assertTrue(isfile(path),"Missing committed execution receipt: "+path);
                receipt = jsondecode(fileread(path));
                testCase.assertEqual(string(receipt.provider),provider);
                validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=testCase.repositoryRoot);
            end
        end
        function executableRevisionMustMatchDeclaredBuild(testCase)
            expected = string(repmat('a',1,40));
            stale = string(repmat('b',1,40));
            for executable = ["/fixture/wave-vortex-run","/fixture/WVForwardIntegrationProbe"]
                report = struct(source=struct(commit=expected));
                TestPortableForwardIntegration.validateExecutableSource(report,executable,expected);
                report.source.commit = stale;
                try
                    TestPortableForwardIntegration.validateExecutableSource(report,executable,expected);
                    testCase.assertFail("A stale executable must be rejected before its report can become a receipt.");
                catch exception
                    testCase.verifyEqual(string(exception.identifier),"WaveVortexModel:ExecutableSourceMismatch");
                    testCase.verifySubstring(string(exception.message),executable);
                    testCase.verifySubstring(string(exception.message),expected);
                    testCase.verifySubstring(string(exception.message),stale);
                end
            end
        end
        function executableRevisionCannotBeMissingOrMalformed(testCase)
            expected = string(repmat('a',1,40));
            reports = {struct(),struct(source=struct()),struct(source=[]), ...
                struct(source=struct(commit=[])),struct(source=struct(commit="")), ...
                struct(source=struct(commit=[expected,expected])),struct(source=struct(commit=missing))};
            for k = 1:numel(reports)
                testCase.verifyError(@()TestPortableForwardIntegration.validateExecutableSource(reports{k},"/fixture/wave-vortex-run",expected),"WaveVortexModel:ExecutableSourceMismatch");
            end
        end
        function collectorRejectsMissingOrMixedExecutionSources(testCase)
            [receipt,root] = testCase.receiptFixture();
            evidenceDirectory = fullfile(root,"fragments");
            mkdir(evidenceDirectory);
            testCase.verifyError(@()collectPortableForwardIntegrationQualification(evidenceDirectory,"reference",repositoryRoot=root),"WaveVortexModel:InvalidForwardIntegrationAssembly");
            receipt.buildSource = struct(root=root,commit=receipt.sourceCommit);
            for item = reshape(receipt.cases,1,[])
                fragment = receipt;
                fragment.cases = item;
                writeText(fullfile(evidenceDirectory,string(item.id)+"-reference.json"),jsonencode(fragment));
            end
            lastPath = fullfile(evidenceDirectory,string(receipt.cases(end).id)+"-reference.json");
            original = jsondecode(fileread(lastPath));
            mixed = original;
            mixed.sourceSHA256(1).sha256 = repmat('0',1,64);
            writeText(lastPath,jsonencode(mixed));
            testCase.verifyError(@()collectPortableForwardIntegrationQualification(evidenceDirectory,"reference",repositoryRoot=root),"WaveVortexModel:InvalidForwardIntegrationAssembly");
            mixed = original;
            mixed.executionBinaries(1).sha256 = repmat('0',1,64);
            writeText(lastPath,jsonencode(mixed));
            testCase.verifyError(@()collectPortableForwardIntegrationQualification(evidenceDirectory,"reference",repositoryRoot=root),"WaveVortexModel:InvalidForwardIntegrationAssembly");
            mixed = original;
            mixed.buildSource.commit = repmat('b',1,40);
            writeText(lastPath,jsonencode(mixed));
            testCase.verifyError(@()collectPortableForwardIntegrationQualification(evidenceDirectory,"reference",repositoryRoot=root),"WaveVortexModel:InvalidForwardIntegrationAssembly");
            writeText(lastPath,jsonencode(original));
            assembled = collectPortableForwardIntegrationQualification(evidenceDirectory,"reference",repositoryRoot=root);
            testCase.verifyEqual(string(assembled.sourceCommit),string(receipt.sourceCommit));
            testCase.verifyEqual(assembled.buildSource,jsondecode(jsonencode(receipt.buildSource)));
            testCase.verifyEqual(numel(assembled.artifacts),6);
            testCase.verifyTrue(isfile(fullfile(root,"PortableRuntime","qualification","forward-integration-reference-v1.json")));
            for artifact = reshape(assembled.artifacts,1,[])
                testCase.verifyEqual(portableForwardIntegrationSHA256(fullfile(root,string(artifact.path))),string(artifact.sha256));
            end
        end
        function supportRowsAndExclusionsAreComplete(testCase)
            validatePortableForwardIntegrationCatalog(testCase.catalog);
            testCase.verifyEqual(numel(testCase.catalog.rows),30);
            testCase.verifyEqual(numel(testCase.catalog.cases),6);
            bad = testCase.catalog; bad.rows(end) = [];
            testCase.verifyInvalidCatalog(bad);
            bad = testCase.catalog; bad.rows(end) = bad.rows(1);
            testCase.verifyInvalidCatalog(bad);
            bad = testCase.catalog; bad.rows(1).witnesses = {'unresolved'};
            testCase.verifyInvalidCatalog(bad);
            bad = testCase.catalog; bad.witnesses(1).symbol = 'missingFunction';
            testCase.verifyInvalidCatalog(bad);
            bad = testCase.catalog; bad.rows(1).support = 'pending';
            testCase.verifyInvalidCatalog(bad);
            bad = testCase.catalog; bad.exclusions(2).detail = 'Default MATLAB output is unsupported';
            testCase.verifyInvalidCatalog(bad);
        end
        function executionRequiresEveryCaseAndExactScope(testCase)
            [receipt,root] = testCase.receiptFixture();
            validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=root);
            bad = receipt; bad.cases(end) = [];
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(end) = bad.cases(1);
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).passed = false;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).provider = 'native-fftw';
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).scope = 'complete-model-readiness';
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).testName = 'TestPortableForwardIntegration/representativeLifecycleMatchesMatlab(configuration=boussinesq)';
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).checks(end) = [];
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).checks(1).passed = false;
            testCase.verifyInvalidReceipt(bad,root);
        end
        function sourceAndArtifactDigestsCannotGoStale(testCase)
            [receipt,root] = testCase.receiptFixture();
            bad = receipt; bad.catalogSHA256 = repmat('0',1,64);
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.sourceSHA256(1).sha256 = repmat('0',1,64);
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.sourceSHA256(end) = [];
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt;
            sampler = string({bad.sourceSHA256.path})=="PortableRuntime/src/WVFieldEvaluationService.cpp";
            bad.sourceSHA256(sampler) = [];
            testCase.verifyInvalidReceipt(bad,root);
            artifact = "PortableRuntime/qualification/fixture.json";
            writeText(fullfile(root,artifact),"{}");
            receipt.artifacts(end+1) = struct(path=artifact,sha256=portableForwardIntegrationSHA256(fullfile(root,artifact)));
            validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=root);
            writeText(fullfile(root,artifact),"{ }");
            testCase.verifyInvalidReceipt(receipt,root);
        end
        function numericalAndRetainedEvidenceCannotBePromoted(testCase)
            [receipt,root] = testCase.receiptFixture();
            bad = receipt; bad.cases(1).evolution.coefficients = 0;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).outcomes.whole(2).tracer = NaN;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).outcomes.stoppedMatlab(2).fields = 3e-7;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).outcomes.split(2) = [];
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).stop.requestedAtAcceptedTime = 138;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).runs.whole.integrationRequest.stepPolicy = 'explicit';
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).denseHeld.recordCount = 0;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.cases(1).runs.stoppedResume.state.initialTime = 17;
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt;
            for k = 1:numel(bad.cases)
                bad.cases(k).controls.stopBoundary = "accepted-step";
                bad.cases(k).stop.requestedBoundary = "accepted-step";
                bad.cases(k).stop.requestedAtOutputTime = [];
                bad.cases(k).runs.stopped.termination = bad.cases(k).stop;
            end
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.artifacts = struct(path={},sha256={});
            testCase.verifyInvalidReceipt(bad,root);
            bad = receipt; bad.executionBinaries(2).role = 'runner';
            testCase.verifyInvalidReceipt(bad,root);
            writeText(fullfile(root,string(receipt.artifacts.path)),"Unrelated passing tests");
            receipt.artifacts.sha256 = portableForwardIntegrationSHA256(fullfile(root,string(receipt.artifacts.path)));
            testCase.verifyInvalidReceipt(receipt,root);
        end
        function inheritedProofPreservesOriginalArtifactAndCorrections(testCase)
            [receipt,root] = testCase.receiptFixture();
            source = receipt.sourceSHA256(1);
            key = matlab.lang.makeValidName(char(source.path));
            original = struct(sourceSHA256=struct);
            original.sourceSHA256.(key) = source.sha256;
            original.sourceSHA256.PortableRuntime_README_md = repmat('0',1,64);
            original.sourceSHA256.PortableRuntime_CMakeLists_txt = repmat('0',1,64);
            original.sourceSHA256.PortableRuntime_tests_TestFixture_cpp = repmat('0',1,64);
            originalPath = "PortableRuntime/qualification/inherited-fixture.json";
            writeText(fullfile(root,originalPath),jsonencode(original));
            entry = struct(path=originalPath,sha256=portableForwardIntegrationSHA256(fullfile(root,originalPath)), ...
                classification="production-source-compatible",productionSources=source,verifiedCorrections=struct(path={},sha256={}));
            receipt.inheritedEvidence = entry;
            validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=root);
            bad = receipt; bad.inheritedEvidence.productionSources = struct(path={},sha256={});
            testCase.verifyInvalidReceipt(bad,root);
            % A correction must retain the original stale digest and record a
            % passing check for the new exact bytes; merely rebinding fails.
            original.sourceSHA256.(key) = repmat('0',1,64);
            original.ciCorrections = struct(file=source.path,correctedSHA256=source.sha256,focusedNativeTest=struct(passed=true));
            writeText(fullfile(root,originalPath),jsonencode(original));
            receipt.inheritedEvidence.sha256 = portableForwardIntegrationSHA256(fullfile(root,originalPath));
            testCase.verifyInvalidReceipt(receipt,root);
            receipt.inheritedEvidence.verifiedCorrections = source;
            validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=root);
            original.ciCorrections.focusedNativeTest.passed = false;
            writeText(fullfile(root,originalPath),jsonencode(original));
            receipt.inheritedEvidence.sha256 = portableForwardIntegrationSHA256(fullfile(root,originalPath));
            testCase.verifyInvalidReceipt(receipt,root);
        end
    end
    methods (Access=private)
        function verifyInvalidCatalog(testCase,catalog)
            testCase.verifyError(@()validatePortableForwardIntegrationCatalog(catalog),"WaveVortexModel:InvalidForwardIntegrationCatalog");
        end
        function verifyInvalidReceipt(testCase,receipt,root)
            testCase.verifyError(@()validatePortableForwardIntegrationQualification(receipt,testCase.catalog,repositoryRoot=root),"WaveVortexModel:InvalidForwardIntegrationQualification");
        end
        function [receipt,root] = receiptFixture(testCase)
            % Synthetic receipts exercise validation only and never qualify
            % scientific execution. Every copied source lives in a fixture.
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = string(fixture.Folder);
            paths = unique([string({testCase.catalog.witnesses.path}),reshape(string(testCase.catalog.relatedCatalogs),1,[]),string({testCase.catalog.historicalEvidence.path}), ...
                "PortableRuntime/contracts/portable-forward-integration-v1.json"]);
            for path = paths
                destination = fullfile(root,path);
                if ~isfolder(fileparts(destination)), mkdir(fileparts(destination)); end
                copyfile(fullfile(testCase.repositoryRoot,path),destination);
            end
            required = ["PortableRuntime/include/WaveVortexRuntime/WVRungeKutta.hpp", ...
                "PortableRuntime/src/WVRungeKutta.cpp","PortableRuntime/src/WVOrderedRKCombination.hpp","PortableRuntime/src/WVModel.cpp", ...
                "PortableRuntime/src/WVOutputOrchestration.cpp","PortableRuntime/src/WVModelOutputNetCDFWriter.cpp", ...
                "PortableRuntime/tests/WVForwardIntegrationProbe.cpp","UnitTests/TestPortableForwardIntegration.m", ...
                "PortableRuntime/app/WaveVortexRun.cpp","PortableRuntime/include/WaveVortexRuntime/WVIntegrationContracts.hpp", ...
                "CompiledKernel/src/WVTransformConstantStratificationKernel.cpp", ...
                "PortableRuntime/src/WVFieldEvaluationService.cpp","PortableRuntime/src/WVStratifiedFieldEvaluationAdapter.cpp", ...
                "PortableRuntime/src/WVBarotropicQGFieldEvaluationAdapter.cpp"];
            required = [required,"PortableRuntime/CMakeLists.txt","CompiledKernel/CMakeLists.txt", ...
                "PortableRuntime/app/WVRunnerVariablePolicy.cpp","PortableRuntime/app/WVRunnerVariablePolicy.hpp"];
            hashes = struct(path={},sha256={});
            for path = required
                content = "function representativeLifecycleMatchesMatlab(testCase,configuration)"+newline+"end"+newline;
                if path=="PortableRuntime/app/WaveVortexRun.cpp"
                    content = content+"% WVRunnerVariablePolicy.hpp"+newline;
                end
                writeText(fullfile(root,path),content);
                hashes(end+1) = struct(path=path,sha256=portableForwardIntegrationSHA256(fullfile(root,path))); %#ok<AGROW>
            end
            checkNames = ["nontrivialEvolution","heldAmplitudes","completeRestartState","allContinuationDirections","twoOutputDestinations","controlledStopResume","denseOutputParity"];
            checks = struct(name=cellstr(checkNames),passed=true);
            cases = struct(id={},configuration={},profile={},testName={},passed={},failed={},incomplete={},provider={},scope={},checks={},controls={},evolution={},outcomes={},denseHeld={},stop={},runs={});
            for item = reshape(testCase.catalog.cases,1,[])
                measured = measuredFixture(item);
                cases(end+1) = struct(id=item.id,configuration=item.configuration,profile=item.profile, ...
                    testName="TestPortableForwardIntegration/"+string(item.symbol)+"(configuration="+string(item.parameter)+")", ...
                    passed=true,failed=false,incomplete=false,provider="reference",scope="forward-integration",checks=checks,controls=measured.controls,evolution=measured.evolution, ...
                    outcomes=measured.outcomes,denseHeld=measured.denseHeld,stop=measured.stop,runs=measured.runs); %#ok<AGROW>
            end
            logPath = "PortableRuntime/qualification/fixture.log";
            lines = "FORWARD_INTEGRATION_PASS "+string({cases.configuration})+" provider=reference synthetic-validation-fixture";
            writeText(fullfile(root,logPath),join(lines,newline));
            artifacts = struct(path=logPath,sha256=portableForwardIntegrationSHA256(fullfile(root,logPath)));
            binaries = struct(role={"runner","probe"},sha256=repmat('a',1,64),sourceCommit=repmat('a',1,40));
            receipt = struct(schema="wave-vortex-forward-integration-qualification-v1",schemaVersion=1,scope="forward-integration",status="complete", ...
                sourceCommit=repmat('a',1,40),provider="reference",matlabRelease="fixture",platform="fixture", ...
                catalogSHA256=portableForwardIntegrationSHA256(fullfile(root,"PortableRuntime","contracts","portable-forward-integration-v1.json")), ...
                sourceSHA256=hashes,artifacts=artifacts,executionBinaries=binaries,cases=cases);
        end
    end
end

function writeText(path,text)
if ~isfolder(fileparts(path)), mkdir(fileparts(path)); end
fid = fopen(path,"w");
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,"%s",text);
end

function value = measuredFixture(definition)
method = string(definition.profile);
policy="adaptive";
if method=="rk4-cfl", method="fixed-rk4"; policy="cfl"; end
controls = struct(initialTime=17,splitTime=77,finalTime=137,initialStep=4,maximumStep=4,relativeTolerance=1e-9,absoluteTolerance=1e-10,method=method,stopBoundary="accepted-step");
stop = struct(reason="stop-requested",requestedBoundary="accepted-step",requestedAtAcceptedTime=21,requestedAtOutputTime=[],finalAcceptedTime=37,callbackEvaluationCount=2);
if ismember(string(definition.configuration),["constant-nonhydrostatic","stratified-qg","boussinesq"])
    controls.stopBoundary = "output-occurrence";
    stop.requestedBoundary = "output-occurrence";
    stop.requestedAtOutputTime = 20;
end
metrics = struct(coefficients=1e-9,fields=1e-9,tracer=1e-9,savedOutput=1e-9,particlePosition=1e-5,heldAmplitude=0,restoredStateExact=true);
outcomes = struct;
for name = ["whole","split","matlabToCpp","cppToMatlab","stoppedCpp","stoppedMatlab"]
    outcomes.(name) = [metrics,metrics];
end
runs = struct;
names = ["whole","firstSegment","secondSegment","matlabToCpp","stopped","stoppedResume"];
initial = [17,17,77,77,17,37]; final = [137,77,137,137,37,137];
for k = 1:numel(names)
    termination = struct(reason="reached-final-time",finalAcceptedTime=final(k));
    status = "complete";
    if names(k)=="stopped", termination=stop; status="stopped"; end
    runs.(names(k)) = struct(status=status,state=struct(initialTime=initial(k),finalTime=final(k),stepCount=5), ...
        termination=termination,integrationRequest=struct(activeMethod=method,stepPolicy=policy),acceptedStepCount=5,denseOutputEvaluationCount=1);
end
value = struct(controls=controls,evolution=struct(coefficients=1e-5,particlePosition=1e-3,tracer=1e-7),outcomes=outcomes, ...
    denseHeld=struct(recordCount=1,firstRecordTime=22,maximumHeldError=0),stop=stop,runs=runs);
end
