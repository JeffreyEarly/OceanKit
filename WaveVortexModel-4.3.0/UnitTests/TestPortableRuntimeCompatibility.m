classdef TestPortableRuntimeCompatibility < matlab.unittest.TestCase
    properties (SetAccess = private)
        RepositoryRoot (1,1) string
        Runner (1,1) string
        TemporaryFolder (1,1) string
    end

    methods (TestClassSetup)
        function buildReferenceRuntime(testCase)
            testCase.RepositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryFolder = string(fixture.Folder);
            buildDirectory = fullfile(testCase.TemporaryFolder,"build");
            configure = "cmake -S " + shellQuote(fullfile(testCase.RepositoryRoot,"PortableRuntime")) + ...
                " -B " + shellQuote(buildDirectory) + ...
                " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF";
            [status,output] = systemWithoutMatlabRuntime(configure);
            testCase.assertEqual(status,0,output)
            [status,output] = systemWithoutMatlabRuntime("cmake --build " + shellQuote(buildDirectory) + ...
                " --parallel --target wave-vortex-run");
            testCase.assertEqual(status,0,output)
            testCase.Runner = fullfile(buildDirectory,"wave-vortex-run");
            testCase.assertTrue(isfile(testCase.Runner))
        end
    end

    methods (Test, TestTags = "full")
        function matlabFixtureRunsAndRuntimeCheckpointRestores(testCase)
            inputPath = fullfile(testCase.RepositoryRoot,"PortableRuntime","tests", ...
                "fixtures","forcing-nonlinear.nc");
            outputPath = fullfile(testCase.TemporaryFolder,"runtime-output.nc");
            command = shellQuote(testCase.Runner) + " " + shellQuote(inputPath) + ...
                " " + shellQuote(outputPath) + ...
                " --restart-mode coefficients --output-policy create" + ...
                " --delta-t 0.01 --steps 1 --fft-provider reference";
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            report = jsondecode(output);
            testCase.verifyEqual(string(report.status),"complete")
            testCase.verifyEqual(string(report.execution.engine),"reference-direct")
            testCase.verifyTrue(report.execution.noFallback)

            [wvt,ncfile] = WVTransform.waveVortexTransformFromFile( ...
                outputPath,shouldReadOnly=true);
            cleanup = onCleanup(@()ncfile.close());
            testCase.verifyEqual(wvt.t,report.state.finalTime,AbsTol=4*eps(report.state.finalTime))
            testCase.verifyEqual([wvt.Nj wvt.Nkl],reshape(report.state.shape,1,[]))
            testCase.verifyFalse(any(~isfinite(wvt.Ap),"all"))
            testCase.verifyFalse(any(~isfinite(wvt.Am),"all"))
            testCase.verifyFalse(any(~isfinite(wvt.A0),"all"))
            clear cleanup
        end

        function matlabModelGraphRoundTripsThroughStandalone(testCase)
            sourcePath = fullfile(testCase.TemporaryFolder,"matlab-model.nc");
            controlPath = fullfile(testCase.TemporaryFolder,"matlab-control.nc");
            source = testCase.createNonlinearModel(sourcePath);
            source.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            source.closeNetCDFFile();
            copyfile(sourcePath,controlPath)

            command = shellQuote(testCase.Runner) + " " + shellQuote(sourcePath) + ...
                " --restart-mode model --output-policy append" + ...
                " --delta-t 1e-4 --final-time 2e-4 --fft-provider reference";
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            report = jsondecode(output);
            testCase.verifyEqual(string(report.status),"complete")
            testCase.verifyEqual(string(report.restartMode),"model")
            testCase.verifyEqual(string(report.outputPolicy),"append")
            testCase.verifySubstring(string(report.execution.schedule),"WVBottomFrictionLinear")
            testCase.verifyTrue(report.execution.noFallback)
            testCase.verifyGreaterThan(report.forcingOperations.physicalFieldReconstructionCount,0)
            testCase.verifyEqual(report.forcingOperations.physicalFieldReconstructionCount,report.forcingOperations.evaluationCount)
            testCase.verifyGreaterThanOrEqual(report.forcingOperations.physicalFieldReuseCount,2*report.forcingOperations.evaluationCount)
            testCase.verifyEqual(report.forcingOperations.spatialTendencyProjectionCount,2*report.forcingOperations.physicalFieldReconstructionCount)
            testCase.verifyGreaterThan(report.forcingOperations.spatialTendencyClearElementWrites,0)

            runtimeModel = WVModel.modelFromFile(char(sourcePath));
            runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
            controlModel = WVModel.modelFromFile(char(controlPath));
            controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
            controlModel.setupIntegrator(integratorType="fixed",deltaT=1e-4);
            controlModel.integrateToTime(2e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);

            testCase.verifyModelGraphsEqual(runtimeModel,controlModel)
            testCase.verifyEqual(runtimeModel.wvt.Ap,controlModel.wvt.Ap,AbsTol=1e-12)
            testCase.verifyEqual(runtimeModel.wvt.Am,controlModel.wvt.Am,AbsTol=1e-12)
            testCase.verifyEqual(runtimeModel.wvt.A0,controlModel.wvt.A0,AbsTol=1e-12)
            [runtimeX,runtimeY,runtimeZ] = runtimeModel.floatPositions();
            [controlX,controlY,controlZ] = controlModel.floatPositions();
            testCase.verifyEqual(runtimeX,controlX,AbsTol=1e-8)
            testCase.verifyEqual(runtimeY,controlY,AbsTol=1e-8)
            testCase.verifyEqual(runtimeZ,controlZ,AbsTol=1e-10)
            testCase.verifyEqual(runtimeModel.tracer("dye"),controlModel.tracer("dye"),AbsTol=1e-8)
            testCase.verifyEqual(runtimeModel.outputFileWithName("matlab-model.nc").outputGroupWithName("wave-vortex").incrementsWrittenToGroup,uint64(3))
            testCase.verifyEqual(runtimeModel.outputFileWithName("matlab-model.nc").outputGroupWithName("particles").incrementsWrittenToGroup,uint64(3))
            testCase.verifyEqual(runtimeModel.outputFileWithName("matlab-model.nc").outputGroupWithName("tracers").incrementsWrittenToGroup,uint64(3))
            clear runtimeCleanup controlCleanup
        end

        function linearAdaptiveGraphRoundTripsThroughStandalone(testCase)
            sourcePath = fullfile(testCase.TemporaryFolder,"linear-model.nc");
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5], ...
                N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
            wvt.Ap(2) = 2e-4 + 1e-4i;
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            outputFile = model.createNetCDFFileForModelOutput(sourcePath, ...
                outputInterval=1e-4,shouldOverwriteExisting=true);
            model.eulerianObservingSystem.addNetCDFOutputVariables('u');
            model.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            command = shellQuote(testCase.Runner) + " " + shellQuote(sourcePath) + ...
                " --restart-mode model --output-policy append" + ...
                " --integrator adaptive-rk23 --relative-tolerance 1e-6" + ...
                " --absolute-tolerance 1e-8 --delta-t 1e-4" + ...
                " --final-time 2e-4 --fft-provider reference";
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            restored = WVModel.modelFromFile(char(sourcePath));
            cleanup = onCleanup(@()restored.closeNetCDFFile());
            testCase.verifyTrue(restored.isDynamicsLinear)
            testCase.verifyEqual(restored.wvt.t,2e-4,AbsTol=1e-14)
            testCase.verifyEqual(restored.wvt.Ap(2),wvt.Ap(2),AbsTol=1e-14)
            testCase.verifyEqual(restored.outputFileWithName("linear-model.nc").outputGroupWithName("wave-vortex").incrementsWrittenToGroup,uint64(3))
            clear cleanup
        end

        function linearBottomFrictionMatchesMatlabFixedAndAdaptive(testCase)
            cases = {
                true, "fixed", 5, 0; ...
                false, "fixed", 7, 2.5e-7; ...
                true, "ode23", 7, 2.5e-7; ...
                false, "ode23", 5, 2.5e-7; ...
                true, "ode45", 5, 2.5e-7; ...
                false, "ode45", 7, 2.5e-7};
            for iCase = 1:size(cases,1)
                isHydrostatic = cases{iCase,1};
                integratorType = cases{iCase,2};
                Nz = cases{iCase,3};
                r = cases{iCase,4};
                caseName = sprintf("linear-%s-h%d-nz%d",integratorType,isHydrostatic,Nz);
                sourcePath = fullfile(testCase.TemporaryFolder,caseName+"-source.nc");
                controlPath = fullfile(testCase.TemporaryFolder,caseName+"-control.nc");
                model = testCase.createLinearBottomFrictionModel(sourcePath,isHydrostatic,Nz,r,integratorType);
                model.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                model.closeNetCDFFile();
                copyfile(sourcePath,controlPath)

                command = shellQuote(testCase.Runner) + " " + shellQuote(sourcePath) + ...
                    " --restart-mode model --output-policy append" + ...
                    " --delta-t 1e-5 --final-time 2e-4 --fft-provider reference";
                if ismember(integratorType,["ode23" "ode45"])
                    command = command + " --integrator adaptive-rk" + extractAfter(integratorType,"ode") + ...
                        " --relative-tolerance 1e-8 --absolute-tolerance 1e-10" + ...
                        " --maximum-step 1e-5";
                end
                [status,output] = systemWithoutMatlabRuntime(command);
                testCase.assertEqual(status,0,output)
                report = jsondecode(output);
                testCase.verifySubstring(string(report.execution.schedule),"WVBottomFrictionLinear")
                testCase.verifyTrue(report.execution.noFallback)
                testCase.verifyEqual(report.forcingOperations.physicalFieldReconstructionCount,report.forcingOperations.evaluationCount)

                runtimeModel = WVModel.modelFromFile(char(sourcePath));
                runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                testCase.configureLinearIntegrator(controlModel,integratorType);
                controlModel.integrateToTime(2e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.verifyEqual(runtimeModel.wvt.forcingWithName("linear bottom friction").r,r)
                testCase.verifyEqual(runtimeModel.wvt.forcingWithName("linear bottom friction").r_scaled,2*(Nz-1)*r,RelTol=10*eps)
                testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError(runtimeModel.wvt,controlModel.wvt),1e-12,caseName)
                clear runtimeCleanup controlCleanup
            end
        end

        function linearOde78EndpointMatchesMatlab(testCase)
            sourcePath = fullfile(testCase.TemporaryFolder,"linear-ode78-source.nc");
            controlPath = fullfile(testCase.TemporaryFolder,"linear-ode78-control.nc");
            runtimePath = fullfile(testCase.TemporaryFolder,"linear-ode78-runtime.nc");
            model = testCase.createLinearBottomFrictionModel(sourcePath,false,7,2.5e-7,"ode78");
            model.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
            copyfile(sourcePath,controlPath)

            command = shellQuote(testCase.Runner) + " " + shellQuote(sourcePath) + ...
                " " + shellQuote(runtimePath) + ...
                " --restart-mode coefficients --output-policy create" + ...
                " --integrator adaptive-rk78 --relative-tolerance 1e-8" + ...
                " --absolute-tolerance 1e-10 --delta-t 1e-5" + ...
                " --initial-step 1e-5 --maximum-step 1e-5" + ...
                " --final-time 2e-4 --fft-provider reference";
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            report = jsondecode(output);
            testCase.verifyEqual(string(report.integrator.id),"adaptive-rk78")
            testCase.verifyEqual(string(report.integrator.controller),"matlab-ode78-v1")
            testCase.verifyEqual(report.integrator.workspaceStateEquivalentCount,11)
            testCase.verifyEqual(report.integrator.denseHistoryStateEquivalentCount,0)

            [runtimeTransform,runtimeFile] = WVTransform.waveVortexTransformFromFile(runtimePath,shouldReadOnly=true);
            runtimeCleanup = onCleanup(@()runtimeFile.close());
            controlModel = WVModel.modelFromFile(char(controlPath));
            controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
            testCase.configureLinearIntegrator(controlModel,"ode78");
            controlModel.integrateToTime(2e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError(runtimeTransform,controlModel.wvt),1e-12)
            clear runtimeCleanup controlCleanup
        end

        function linearOde78ContinuousOutputMatchesMatlab(testCase)
            sourcePath = fullfile(testCase.TemporaryFolder,"linear-ode78-dense-source.nc");
            controlPath = fullfile(testCase.TemporaryFolder,"linear-ode78-dense-control.nc");
            outputDirectory = fullfile(testCase.TemporaryFolder,"linear-ode78-dense-output");
            model = testCase.createLinearBottomFrictionModel(sourcePath,false,7,2.5e-7,"ode78");
            model.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
            copyfile(sourcePath,controlPath)

            requestedTimes = [1.25e-4 1.5e-4 1.75e-4 2e-4];
            command = shellQuote(testCase.Runner) + " " + shellQuote(sourcePath) + ...
                " --restart-mode coefficients --output-policy create" + ...
                " --integrator adaptive-rk78 --relative-tolerance 1e-8" + ...
                " --absolute-tolerance 1e-10 --delta-t 1e-4" + ...
                " --initial-step 1e-4 --maximum-step 1e-4" + ...
                " --final-time 2e-4 --fft-provider reference";
            for requestedTime = requestedTimes
                command = command + " --output-time " + string(requestedTime);
            end
            command = command + " --output-directory " + shellQuote(outputDirectory);
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            report = jsondecode(output);
            testCase.verifyEqual(report.integrator.baseRightHandSideEvaluationCount,13)
            testCase.verifyEqual(report.integrator.continuousExtensionRightHandSideEvaluationCount,4)
            testCase.verifyEqual(report.integrator.denseOutputCacheBuildCount,1)
            testCase.verifyEqual(report.integrator.denseOutputCacheReuseCount,2)
            testCase.verifyEqual(report.integrator.retainedBaseStageStateEquivalentCount,8)
            testCase.verifyEqual(report.integrator.continuousExtensionWorkspaceStateEquivalentCount,4)
            testCase.verifyEqual(report.integrator.workspaceMaximumLiveStateEquivalentCount,15)

            runtimeFiles = dir(fullfile(outputDirectory,"*.nc"));
            testCase.assertNumElements(runtimeFiles,numel(requestedTimes))
            controlModel = WVModel.modelFromFile(char(controlPath));
            controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
            testCase.configureLinearIntegrator(controlModel,"ode78");
            for iTime = 1:numel(requestedTimes)
                controlModel.integrateToTime(requestedTimes(iTime), ...
                    shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                runtimePath = fullfile(runtimeFiles(iTime).folder,runtimeFiles(iTime).name);
                [runtimeTransform,runtimeFile] = WVTransform.waveVortexTransformFromFile( ...
                    runtimePath,shouldReadOnly=true);
                runtimeCleanup = onCleanup(@()runtimeFile.close());
                testCase.verifyEqual(runtimeTransform.t,requestedTimes(iTime),AbsTol=4*eps(requestedTimes(iTime)))
                testCase.verifyLessThanOrEqual( ...
                    testCase.normalizedCoefficientError(runtimeTransform,controlModel.wvt),1e-12)
                clear runtimeCleanup
            end
            clear controlCleanup
        end

        function matlabCFLReferenceFixturesMatchPortableConstants(testCase)
            cfl = 0.4;
            fixtures = {
                false, 0.1020286099012897, 23522.813868795667, 3.2273103327544224e8, 0.0013924156036910925, 1804.9741156372481; ...
                true, 0.10202859972462885, 23522.816215036815, 3.117452951358881e8, 0.0014450221120718154, 1739.2634354005847};
            for iFixture = 1:size(fixtures,1)
                isHydrostatic = fixtures{iFixture,1};
                wvt = WVTransformConstantStratification([15000 12000 1300],[6 5 7], ...
                    N0=5.2e-3,rho0=1027,g=9.80665,planetaryRadius=6.3712e6, ...
                    rotationRate=7.292115e-5,latitude=33, ...
                    isHydrostatic=isHydrostatic,shouldAntialias=true);
                testCase.setConstantCFLFixtureState(wvt)
                candidates = testCase.matlabCFLCandidates(wvt,cfl);
                testCase.verifyEqual(candidates.effectiveHorizontalGridResolution,6000.0000000000009,RelTol=1e-13)
                testCase.verifyEqual(candidates.maximumHorizontalSpeed,fixtures{iFixture,2},RelTol=1e-12)
                testCase.verifyEqual(candidates.horizontalAdvective,fixtures{iFixture,3},RelTol=1e-12)
                testCase.verifyEqual(candidates.verticalAdvective,fixtures{iFixture,4},RelTol=1e-12)
                testCase.verifyEqual(candidates.highestActiveWaveFrequency,fixtures{iFixture,5},RelTol=1e-12)
                testCase.verifyEqual(candidates.oscillatory,fixtures{iFixture,6},RelTol=1e-12)
            end

            qg = WVTransformBarotropicQG([15000 12000],[6 5],h=0.8,j=1, ...
                g=9.80665,planetaryRadius=6.3712e6,rotationRate=7.292115e-5, ...
                latitude=33,shouldAntialias=true);
            p = reshape(1:numel(qg.A0),[],1);
            qg.A0(:) = complex(1e-5*sin(0.19*(p+4)),1e-5*cos(0.05*(p+5)));
            testCase.verifyEqual(qg.effectiveHorizontalGridResolution,6000.0000000000009,RelTol=1e-13)
            testCase.verifyEqual(qg.uvMax,0.11956917129957827,RelTol=1e-12)
            testCase.verifyEqual(cfl*qg.effectiveHorizontalGridResolution/qg.uvMax,20072.063508635067,RelTol=1e-12)
            testCase.verifyFalse(qg.hasWaveComponent)
        end

        function runRequestV2MatchesCFLAndExactMatlabSolvers(testCase)
            sourcePath = fullfile(testCase.TemporaryFolder,"request-v2-source.nc");
            model = testCase.createLinearBottomFrictionModel(sourcePath,false,7,2.5e-7,"fixed");
            model.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            sourceTime = model.wvt.t;
            fixedFinalTime = sourceTime+1e-4;
            adaptiveFinalTime = sourceTime+2e-4;
            candidates = testCase.matlabCFLCandidates(model.wvt,0.25);
            selectedStep = min(candidates.advective,candidates.oscillatory);
            model.closeNetCDFFile();

            explicitPath = fullfile(testCase.TemporaryFolder,"request-v2-explicit.nc");
            cflPath = fullfile(testCase.TemporaryFolder,"request-v2-cfl.nc");
            copyfile(sourcePath,explicitPath)
            copyfile(sourcePath,cflPath)
            explicit = struct(method="fixed-rk4",finalTime=fixedFinalTime,initialStep=selectedStep);
            explicitReport = testCase.runV2Request(explicitPath,explicit,"request-v2-explicit-report.json");
            cflIntegration = struct(method="fixed-rk4",finalTime=fixedFinalTime,cfl=0.25,timeStepConstraint="min");
            cflReport = testCase.runV2Request(cflPath,cflIntegration,"request-v2-cfl-report.json");
            testCase.verifyEqual(string(cflReport.integrationRequest.requestedMethod),"fixed-rk4")
            testCase.verifyEqual(string(cflReport.integrationRequest.activeMethod),"fixed-rk4")
            testCase.verifyEqual(string(cflReport.integrationRequest.stepPolicy),"cfl")
            testCase.verifyTrue(cflReport.integrationRequest.noFallback)
            testCase.verifyEqual(cflReport.integrationRequest.candidates.horizontalAdvective,candidates.horizontalAdvective,RelTol=1e-12)
            if isinf(candidates.verticalAdvective)
                testCase.verifyEqual(string(cflReport.integrationRequest.candidates.verticalAdvective),"infinity")
            else
                testCase.verifyEqual(cflReport.integrationRequest.candidates.verticalAdvective,candidates.verticalAdvective,RelTol=1e-12)
            end
            testCase.verifyEqual(cflReport.integrationRequest.candidates.advective,candidates.advective,RelTol=1e-12)
            testCase.verifyEqual(cflReport.integrationRequest.candidates.oscillatory,candidates.oscillatory,RelTol=1e-12)
            testCase.verifyEqual(cflReport.integrationRequest.selectedStep,selectedStep,RelTol=1e-12)
            testCase.verifyGreaterThan(cflReport.integrationRequest.transientWorkspaceMaximumLiveBytes,0)
            testCase.verifyTrue(all(isfield(cflReport.timingSeconds,["parse" "preflight" "provider" "startup"])))
            explicitModel = WVModel.modelFromFile(char(explicitPath));
            explicitCleanup = onCleanup(@()explicitModel.closeNetCDFFile());
            cflModel = WVModel.modelFromFile(char(cflPath));
            cflCleanup = onCleanup(@()cflModel.closeNetCDFFile());
            testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError(explicitModel.wvt,cflModel.wvt),1e-12)
            testCase.verifyEqual(string(explicitReport.integrationRequest.stepPolicy),"explicit")
            testCase.verifyEqual(explicitReport.integrationRequest.selectedStep,selectedStep,RelTol=1e-12)
            testCase.verifyEqual(explicitReport.integrator.requestedInitialStep,selectedStep,RelTol=1e-12)
            clear explicitCleanup cflCleanup

            methods = ["ode23" "ode45" "ode78"];
            identifiers = ["adaptive-rk23" "adaptive-rk45" "adaptive-rk78"];
            controllers = ["matlab-ode23-v1" "matlab-ode45-v1" "matlab-ode78-v1"];
            for iMethod = 1:numel(methods)
                runtimePath = fullfile(testCase.TemporaryFolder,"request-v2-"+methods(iMethod)+"-runtime.nc");
                controlPath = fullfile(testCase.TemporaryFolder,"request-v2-"+methods(iMethod)+"-control.nc");
                copyfile(sourcePath,runtimePath)
                copyfile(sourcePath,controlPath)
                integration = struct(method=identifiers(iMethod),finalTime=adaptiveFinalTime, ...
                    initialStep=2e-4,maximumStep=2e-4, ...
                    relativeTolerance=1e-8,absoluteToleranceScale=1e-10);
                report = testCase.runV2Request(runtimePath,integration,"request-v2-"+methods(iMethod)+"-report.json");
                testCase.verifyEqual(string(report.integrationRequest.requestedMethod),identifiers(iMethod))
                testCase.verifyEqual(string(report.integrationRequest.activeMethod),identifiers(iMethod))
                testCase.verifyEqual(string(report.integrationRequest.controller),controllers(iMethod))
                testCase.verifyTrue(report.integrationRequest.noFallback)
                testCase.verifyEqual(report.integrator.requestedInitialStep,integration.initialStep,RelTol=1e-12)
                testCase.verifyEqual(report.integrator.requestedMaximumStep,integration.maximumStep,RelTol=1e-12)
                testCase.verifyEqual(report.integrator.relativeTolerance,integration.relativeTolerance)
                testCase.verifyEqual(report.integrator.absoluteTolerance,integration.absoluteToleranceScale)
                if methods(iMethod) == "ode78"
                    testCase.verifyEqual(report.integrator.continuousExtensionRightHandSideEvaluationCount,4)
                    testCase.verifyEqual(report.integrator.denseOutputCacheBuildCount,1)
                end
                runtimeModel = WVModel.modelFromFile(char(runtimePath));
                runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                testCase.configureLinearIntegrator(controlModel,methods(iMethod));
                controlModel.integrateToTime(adaptiveFinalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError(runtimeModel.wvt,controlModel.wvt),1e-12,methods(iMethod))
                clear runtimeCleanup controlCleanup
            end
        end

        function omittedDefaultsMatchExplicitMatlabBehavior(testCase)
            constantSource = fullfile(testCase.TemporaryFolder,"defaults-constant-source.nc");
            constantModel = testCase.createLinearBottomFrictionModel( ...
                constantSource,false,7,2.5e-7,"fixed");
            constantModel.integrateToTime(1e-4, ...
                shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            constantTime = constantModel.wvt.t;
            constantCandidates = testCase.matlabCFLCandidates(constantModel.wvt,.5);
            constantModel.closeNetCDFFile();

            qgSource = fullfile(testCase.TemporaryFolder,"defaults-qg-source.nc");
            qgModel = testCase.createBarotropicQGModel(qgSource,[5 4],1,true);
            qgModel.integrateToTime(.005, ...
                shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            qgTime = qgModel.wvt.t;
            qgCandidates = testCase.matlabCFLCandidates(qgModel.wvt,.5);
            qgModel.closeNetCDFFile();

            sourcePaths = [constantSource qgSource];
            sourceTimes = [constantTime qgTime];
            finalTimes = [constantTime+1e-4 qgTime+.005];
            candidates = {constantCandidates,qgCandidates};
            labels = ["constant" "qg"];
            for iCase = 1:numel(labels)
                label = labels(iCase);
                defaultPath = fullfile(testCase.TemporaryFolder, ...
                    "defaults-"+label+"-runtime.nc");
                explicitPath = fullfile(testCase.TemporaryFolder, ...
                    "defaults-"+label+"-explicit.nc");
                controlPath = fullfile(testCase.TemporaryFolder, ...
                    "defaults-"+label+"-matlab.nc");
                copyfile(sourcePaths(iCase),defaultPath)
                copyfile(sourcePaths(iCase),explicitPath)
                copyfile(sourcePaths(iCase),controlPath)
                selectedStep = min(candidates{iCase}.advective, ...
                    candidates{iCase}.oscillatory);
                maximumStep = .1*(finalTimes(iCase)-sourceTimes(iCase));
                effectiveInitialStep = min(selectedStep,maximumStep);

                defaultRequest = fullfile(testCase.TemporaryFolder, ...
                    "defaults-"+label+".json");
                defaultReportName = "defaults-"+label+"-report.json";
                WVModel.writePortableRunRequest(defaultRequest,defaultPath, ...
                    finalTime=finalTimes(iCase),fftProvider="reference",threads=1, ...
                    reportPath=defaultReportName);
                immutableDefaultRequest = testCase.fileBytes(defaultRequest);
                [status,output] = systemWithoutMatlabRuntime( ...
                    shellQuote(testCase.Runner)+" --request "+shellQuote(defaultRequest));
                testCase.assertEqual(status,0,output)
                testCase.verifyEqual(testCase.fileBytes(defaultRequest), ...
                    immutableDefaultRequest,label+" default request")
                defaultReport = jsondecode(fileread(fullfile( ...
                    testCase.TemporaryFolder,defaultReportName)));

                explicitRequest = fullfile(testCase.TemporaryFolder, ...
                    "defaults-"+label+"-explicit.json");
                explicitReportName = "defaults-"+label+"-explicit-report.json";
                WVModel.writePortableRunRequest(explicitRequest,explicitPath, ...
                    method="adaptive-rk78",finalTime=finalTimes(iCase), ...
                    initialStep=selectedStep,maximumStep=maximumStep, ...
                    relativeTolerance=1e-3,absoluteToleranceScale=1e-6, ...
                    fftProvider="reference",threads=1,reportPath=explicitReportName);
                [status,output] = systemWithoutMatlabRuntime( ...
                    shellQuote(testCase.Runner)+" --request "+shellQuote(explicitRequest));
                testCase.assertEqual(status,0,output)
                explicitReport = jsondecode(fileread(fullfile( ...
                    testCase.TemporaryFolder,explicitReportName)));

                testCase.verifyEmpty(defaultReport.integrationRequest.requestedMethod,label)
                testCase.verifyEqual(string(defaultReport.integrationRequest.activeMethod), ...
                    "adaptive-rk78",label)
                testCase.verifyEqual(string(defaultReport.integrationRequest.controller), ...
                    "matlab-ode78-v1",label)
                testCase.verifyEqual(string(defaultReport.integrationRequest.initialStepPolicy), ...
                    "cfl-0.5-default",label)
                testCase.verifyEqual(defaultReport.integrationRequest.cfl,.5,label)
                testCase.verifyEmpty(defaultReport.integrationRequest.requestedInitialStep,label)
                testCase.verifyEmpty(defaultReport.integrationRequest.requestedMaximumStep,label)
                testCase.verifyEmpty(defaultReport.integrationRequest.requestedRelativeTolerance,label)
                testCase.verifyEmpty(defaultReport.integrationRequest.requestedAbsoluteToleranceScale,label)
                testCase.verifyEqual(defaultReport.integrationRequest.activeRelativeTolerance,1e-3,label)
                testCase.verifyEqual(defaultReport.integrationRequest.activeAbsoluteToleranceScale,1e-6,label)
                testCase.verifyEqual(defaultReport.integrationRequest.selectedStep,selectedStep,RelTol=1e-12)
                testCase.verifyEqual(defaultReport.integrationRequest.effectiveIntegrationStep, ...
                    effectiveInitialStep,RelTol=1e-12)
                testCase.verifyEqual(defaultReport.integrationRequest.effectiveMaximumStep, ...
                    maximumStep,RelTol=1e-12)
                testCase.verifyEqual(string(defaultReport.provider.requestedId),"reference",label)
                testCase.verifyEqual(string(defaultReport.provider.id),"reference",label)
                testCase.verifyEqual(defaultReport.provider.requestedThreads,1,label)
                testCase.verifyEqual(defaultReport.provider.threads,1,label)
                testCase.verifyEqual(explicitReport.integrationRequest.effectiveIntegrationStep, ...
                    effectiveInitialStep,RelTol=1e-12)

                defaultModel = WVModel.modelFromFile(char(defaultPath));
                defaultCleanup = onCleanup(@()defaultModel.closeNetCDFFile());
                explicitModel = WVModel.modelFromFile(char(explicitPath));
                explicitCleanup = onCleanup(@()explicitModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                controlModel.setupIntegrator();
                controlModel.integrateToTime(finalTimes(iCase), ...
                    shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError( ...
                    defaultModel.wvt,explicitModel.wvt),1e-12,label+" explicit")
                testCase.verifyLessThanOrEqual(testCase.normalizedCoefficientError( ...
                    defaultModel.wvt,controlModel.wvt),1e-12,label+" MATLAB")
                clear defaultCleanup explicitCleanup controlCleanup
            end

            nativePath = fullfile(testCase.TemporaryFolder,"defaults-native-unavailable.nc");
            copyfile(constantSource,nativePath)
            nativeRequest = fullfile(testCase.TemporaryFolder,"defaults-native-unavailable.json");
            nativeReportName = "defaults-native-unavailable-report.json";
            WVModel.writePortableRunRequest(nativeRequest,nativePath, ...
                finalTime=constantTime+1e-4,reportPath=nativeReportName);
            sourceBytes = testCase.fileBytes(nativePath);
            requestBytes = testCase.fileBytes(nativeRequest);
            [status,~] = systemWithoutMatlabRuntime( ...
                shellQuote(testCase.Runner)+" --request "+shellQuote(nativeRequest));
            testCase.verifyNotEqual(status,0)
            testCase.verifyEqual(testCase.fileBytes(nativePath),sourceBytes)
            testCase.verifyEqual(testCase.fileBytes(nativeRequest),requestBytes)
            failure = jsondecode(fileread(fullfile( ...
                testCase.TemporaryFolder,nativeReportName)));
            testCase.verifyEqual(string(failure.status),"failed")
            testCase.verifyEqual(string(failure.failure.stage),"provider")
            testCase.verifySubstring(string(failure.failure.message), ...
                "PortableRuntime/buildWaveVortexRun.sh")
        end

        function barotropicQGModelOutputRoundTripsThroughStandalone(testCase)
            cases = {
                [5 4], 0, false; ...
                [6 5], 1, true; ...
                [5 6], 1, false; ...
                [6 4], 0, true};
            fieldNames = ["u" "v" "eta" "pi" "psi" "qgpv" "zeta_z" "ssh"];
            maximumCoefficientError = 0;
            maximumFieldError = 0;
            maximumParticleError = 0;
            maximumTracerError = 0;
            for iCase = 1:size(cases,1)
                Nxy = cases{iCase,1};
                j = cases{iCase,2};
                shouldAntialias = cases{iCase,3};
                caseName = sprintf("qg-%dx%d-j%d-a%d",Nxy(1),Nxy(2),j,shouldAntialias);
                runtimePath = fullfile(testCase.TemporaryFolder,caseName+"-runtime.nc");
                controlPath = fullfile(testCase.TemporaryFolder,caseName+"-control.nc");
                source = testCase.createBarotropicQGModel(runtimePath,Nxy,j,shouldAntialias);
                source.integrateToTime(0.01,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                source.closeNetCDFFile();
                copyfile(runtimePath,controlPath)

                requestPath = fullfile(testCase.TemporaryFolder,caseName+"-request.json");
                reportName = caseName+"-report.json";
                WVModel.writePortableRunRequest(requestPath,runtimePath, ...
                    schemaVersion=2,method="fixed-rk4",finalTime=.02,initialStep=.0025, ...
                    outputPolicy="append",fftProvider="reference",threads=1,reportPath=reportName);
                immutableRequest = fileread(requestPath);
                command = shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath);
                [status,output] = systemWithoutMatlabRuntime(command);
                testCase.assertEqual(status,0,output)
                testCase.verifyEqual(fileread(requestPath),immutableRequest)
                report = jsondecode(fileread(fullfile(testCase.TemporaryFolder,reportName)));
                testCase.verifyEqual(string(report.status),"complete",caseName)
                testCase.verifyEqual(string(report.execution.engine),"reference-direct",caseName)
                testCase.verifyTrue(report.execution.noFallback,caseName)
                testCase.verifyEqual(report.execution.planCount,3,caseName)
                testCase.verifyEqual(report.storageBytes.persistentFullHermitian,0,caseName)
                testCase.verifyGreaterThan(report.forcingOperations.physicalFieldReuseCount,0,caseName)
                testCase.verifyGreaterThan(report.integrationBreakdownSeconds.tracerAdvection,0,caseName)
                testCase.verifyGreaterThan(report.integrationBreakdownSeconds.particleAdvection,0,caseName)

                runtimeModel = WVModel.modelFromFile(char(runtimePath));
                runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                controlModel.setupIntegrator(integratorType="fixed",deltaT=0.0025);
                controlModel.integrateToTime(0.02,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                [coefficientError,fieldError,particleError,tracerError] = ...
                    testCase.verifyBarotropicQGModelsEqual(runtimeModel,controlModel,fieldNames,caseName);
                maximumCoefficientError = max(maximumCoefficientError,coefficientError);
                maximumFieldError = max(maximumFieldError,fieldError);
                maximumParticleError = max(maximumParticleError,particleError);
                maximumTracerError = max(maximumTracerError,tracerError);
                testCase.verifyModelGraphsEqual(runtimeModel,controlModel)
                runtimeForcingNames = runtimeModel.wvt.forcingNames;
                controlForcingNames = controlModel.wvt.forcingNames;
                testCase.verifyEqual(string(runtimeForcingNames(:)), ...
                    string(controlForcingNames(:)),caseName)
                testCase.verifyEqual(runtimeModel.outputFiles(1).outputGroups(1).incrementsWrittenToGroup,uint64(5),caseName)
                testCase.verifyNetCDFSchemasEqual(runtimePath,controlPath,caseName)

                for iGroup = 1:numel(runtimeModel.outputFiles(1).outputGroups)
                    runtimeModel.outputFiles(1).outputGroups(iGroup).finalTime = runtimeModel.wvt.t;
                    controlModel.outputFiles(1).outputGroups(iGroup).finalTime = controlModel.wvt.t;
                end
                runtimeModel.setupIntegrator(integratorType="fixed",deltaT=0.0025);
                runtimeFinalTime = runtimeModel.wvt.t + 4*0.0025;
                controlFinalTime = controlModel.wvt.t + 4*0.0025;
                runtimeModel.integrateToTime(runtimeFinalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                controlModel.integrateToTime(controlFinalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                [coefficientError,fieldError,particleError,tracerError] = ...
                    testCase.verifyBarotropicQGModelsEqual(runtimeModel,controlModel,fieldNames,caseName+" MATLAB continuation");
                maximumCoefficientError = max(maximumCoefficientError,coefficientError);
                maximumFieldError = max(maximumFieldError,fieldError);
                maximumParticleError = max(maximumParticleError,particleError);
                maximumTracerError = max(maximumTracerError,tracerError);
                clear runtimeCleanup controlCleanup
            end
            fprintf("Barotropic QG model-output MATLAB/C++ coefficient error: %.3e\n",maximumCoefficientError)
            fprintf("Barotropic QG model-output MATLAB/C++ field error: %.3e\n",maximumFieldError)
            fprintf("Barotropic QG model-output MATLAB/C++ particle error: %.3e\n",maximumParticleError)
            fprintf("Barotropic QG model-output MATLAB/C++ tracer error: %.3e\n",maximumTracerError)
        end

        function matlabWriterRequestsDecodeAndRunEveryIntegrationForm(testCase)
            sourceFixture = fullfile(testCase.TemporaryFolder,"writer-source.nc");
            sourceModel = testCase.createLinearBottomFrictionModel(sourceFixture,true,5,0,"fixed");
            sourceModel.integrateToTime(1e-4,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            sourceModel.closeNetCDFFile();
            forms = {
                1,"fixed-rk4","explicit";
                1,"adaptive-rk23","adaptive";
                2,"fixed-rk4","explicit";
                2,"fixed-rk4","advective";
                2,"fixed-rk4","oscillatory";
                2,"fixed-rk4","min";
                2,"adaptive-rk23","adaptive";
                2,"adaptive-rk45","adaptive";
                2,"adaptive-rk78","adaptive";
                };
            for iForm = 1:size(forms,1)
                schemaVersion = forms{iForm,1};
                method = forms{iForm,2};
                policy = forms{iForm,3};
                stem = "writer-v"+schemaVersion+"-"+method+"-"+policy;
                modelPath = fullfile(testCase.TemporaryFolder,stem+".nc");
                requestPath = fullfile(testCase.TemporaryFolder,stem+".json");
                reportPath = stem+"-report.json";
                copyfile(sourceFixture,modelPath)
                ncwriteatt(modelPath,"/","portableFileIdentifier",stem);
                times = ncread(modelPath,"/wave-vortex/t");
                finalTime = times(end);
                if policy == "explicit"
                    WVModel.writePortableRunRequest(requestPath,modelPath, ...
                        schemaVersion=schemaVersion,method=method,finalTime=finalTime,initialStep=.25, ...
                        fftProvider="reference",threads=1,reportPath=reportPath);
                elseif method == "fixed-rk4"
                    WVModel.writePortableRunRequest(requestPath,modelPath, ...
                        schemaVersion=schemaVersion,method=method,finalTime=finalTime,cfl=.25, ...
                        timeStepConstraint=policy,fftProvider="reference",threads=1,reportPath=reportPath);
                else
                    WVModel.writePortableRunRequest(requestPath,modelPath, ...
                        schemaVersion=schemaVersion,method=method,finalTime=finalTime, ...
                        initialStep=.25,maximumStep=1,relativeTolerance=1e-3, ...
                        absoluteToleranceScale=1e-6,fftProvider="reference",threads=1, ...
                        reportPath=reportPath);
                end
                immutableRequest = fileread(requestPath);
                command = shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath);
                [status,output] = systemWithoutMatlabRuntime(command);
                testCase.assertEqual(status,0,output)
                testCase.verifyEqual(fileread(requestPath),immutableRequest)
                report = jsondecode(fileread(fullfile(testCase.TemporaryFolder,reportPath)));
                testCase.verifyEqual(string(report.status),"complete")
                testCase.verifyEqual(string(report.request.schemaIdentifier),"wave-vortex-run-request-v"+schemaVersion)
                testCase.verifyEqual(report.request.schemaVersion,schemaVersion)
                testCase.verifyEqual(string(report.integrator.id),method)
                if schemaVersion == 2
                    testCase.verifyEqual(string(report.integrationRequest.requestedMethod),method)
                    testCase.verifyTrue(report.integrationRequest.noFallback)
                end
            end
        end

        function matlabWriterBarotropicQGRequestsDecodeAndRun(testCase)
            cases = {
                [5 4], 0, false; ...
                [6 5], 1, true; ...
                [5 6], 1, false; ...
                [6 4], 0, true};
            forms = {
                "fixed-rk4", "explicit"; ...
                "fixed-rk4", "cfl"; ...
                "adaptive-rk23", "adaptive"; ...
                "adaptive-rk45", "adaptive"; ...
                "adaptive-rk78", "adaptive"};
            fieldNames = ["u" "v" "eta" "pi" "psi" "qgpv" "zeta_z" "ssh"];
            for iCase = 1:size(cases,1)
                Nxy = cases{iCase,1};
                j = cases{iCase,2};
                shouldAntialias = cases{iCase,3};
                caseName = sprintf("qg-%dx%d-j%d-a%d",Nxy(1),Nxy(2),j,shouldAntialias);
                sourceFixture = fullfile(testCase.TemporaryFolder,"writer-"+caseName+"-source.nc");
                sourceModel = testCase.createBarotropicQGModel(sourceFixture,Nxy,j,shouldAntialias,.00125);
                sourceModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                sourceModel.closeNetCDFFile();
                for iForm = 1:size(forms,1)
                    method = forms{iForm,1};
                    stepPolicy = forms{iForm,2};
                    stem = "writer-"+caseName+"-"+method+"-"+stepPolicy;
                    modelPath = fullfile(testCase.TemporaryFolder,stem+"-runtime.nc");
                    controlPath = fullfile(testCase.TemporaryFolder,stem+"-control.nc");
                    requestPath = fullfile(testCase.TemporaryFolder,stem+".json");
                    reportName = stem+"-report.json";
                    copyfile(sourceFixture,modelPath)
                    copyfile(sourceFixture,controlPath)
                    ncwriteatt(modelPath,"/","portableFileIdentifier",stem);
                    ncwriteatt(controlPath,"/","portableFileIdentifier",stem);
                    if stepPolicy == "explicit"
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,initialStep=.0025, ...
                            outputPolicy="append",fftProvider="reference",threads=1,reportPath=reportName);
                    elseif stepPolicy == "cfl"
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,cfl=.25, ...
                            timeStepConstraint="advective",outputPolicy="append", ...
                            fftProvider="reference",threads=1,reportPath=reportName);
                    else
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,initialStep=.005,maximumStep=.005, ...
                            relativeTolerance=1e-8,absoluteToleranceScale=1e-10, ...
                            outputPolicy="append",fftProvider="reference",threads=1,reportPath=reportName);
                    end
                    immutableRequest = fileread(requestPath);
                    command = shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath);
                    [status,output] = systemWithoutMatlabRuntime(command);
                    testCase.assertEqual(status,0,output)
                    testCase.verifyEqual(fileread(requestPath),immutableRequest)
                    report = jsondecode(fileread(fullfile(testCase.TemporaryFolder,reportName)));
                    testCase.verifyEqual(string(report.status),"complete",stem)
                    testCase.verifyEqual(string(report.integrationRequest.requestedMethod),method,stem)
                    testCase.verifyEqual(string(report.integrationRequest.activeMethod),method,stem)
                    testCase.verifyEqual(string(report.integrationRequest.stepPolicy),stepPolicy,stem)
                    testCase.verifyTrue(report.integrationRequest.noFallback,stem)
                    testCase.verifyEqual(string(report.execution.engine),"reference-direct",stem)
                    testCase.verifyTrue(report.execution.noFallback,stem)
                    testCase.verifyEqual(report.execution.planCount,3,stem)
                    testCase.verifyEqual(report.storageBytes.persistentFullHermitian,0,stem)
                    testCase.verifyGreaterThan(report.storageBytes.fullModelPersistent,0,stem)
                    testCase.verifyGreaterThanOrEqual(report.livenessBytes.fullModelMaximumLive, ...
                        report.livenessBytes.fullModelRetained,stem)
                    testCase.verifyGreaterThan(report.integrator.denseOutputEvaluationCount,0,stem)
                    if method == "adaptive-rk78"
                        testCase.verifyGreaterThan( ...
                            report.integrator.continuousExtensionRightHandSideEvaluationCount,0,stem)
                    end

                    runtimeModel = WVModel.modelFromFile(char(modelPath));
                    runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                    controlModel = WVModel.modelFromFile(char(controlPath));
                    controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                    if stepPolicy == "explicit"
                        controlModel.setupIntegrator(integratorType="fixed",deltaT=.0025);
                    elseif stepPolicy == "cfl"
                        controlModel.setupIntegrator(integratorType="fixed", ...
                            deltaT=report.integrationRequest.selectedStep);
                    else
                        matlabIntegrator = str2func("ode"+extractAfter(method,"adaptive-rk"));
                        controlModel.setupIntegrator(integratorType="adaptive",integrator=matlabIntegrator, ...
                            absTolerance=1e-10,relTolerance=1e-8);
                        controlModel.odeOptions = odeset(controlModel.odeOptions, ...
                            'InitialStep',.005,'MaxStep',.005);
                    end
                    controlModel.integrateToTime(.01,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                    testCase.verifyBarotropicQGModelsEqual(runtimeModel,controlModel,fieldNames,stem);
                    testCase.verifyModelGraphsEqual(runtimeModel,controlModel)
                    testCase.verifyEqual(runtimeModel.outputFiles(1).outputGroups(1).incrementsWrittenToGroup, ...
                        uint64(9),stem)
                    outputInformation = ncinfo(modelPath,"/wave-vortex");
                    variableNames = string({outputInformation.Variables.Name});
                    testCase.verifyFalse(any(startsWith(variableNames,["Ap" "Am"])),stem)
                    runtimeModel.closeNetCDFFile();
                    controlModel.closeNetCDFFile();
                    clear runtimeCleanup controlCleanup
                    testCase.verifyBarotropicQGOutputSeriesEqual(modelPath,controlPath,fieldNames,stem)
                    testCase.verifyNetCDFSchemasEqual(modelPath,controlPath,stem)
                end
            end
        end

        function matlabWriterBarotropicQGForcingMatrixMatchesMatlab(testCase)
            forcingCases = {
                "nonlinear", "nonlinear"; ...
                "adaptive", "adaptive"; ...
                "fixed", "fixed"; ...
                "narrow", "narrow"; ...
                "linear", "linear"; ...
                "quadratic", "quadratic"; ...
                "beta", "beta"; ...
                "ordered-fixed", ["quadratic" "linear" "beta" "nonlinear" "adaptive" "fixed"]; ...
                "ordered-narrow", ["beta" "quadratic" "nonlinear" "adaptive" "narrow"]};
            fieldNames = ["u" "v" "eta" "pi" "psi" "qgpv" "zeta_z" "ssh"];
            for iCase = 1:size(forcingCases,1)
                caseName = "qg-forcing-"+forcingCases{iCase,1};
                sourcePath = fullfile(testCase.TemporaryFolder,caseName+"-source.nc");
                runtimePath = fullfile(testCase.TemporaryFolder,caseName+"-runtime.nc");
                controlPath = fullfile(testCase.TemporaryFolder,caseName+"-control.nc");
                requestPath = fullfile(testCase.TemporaryFolder,caseName+"-request.json");
                reportName = caseName+"-report.json";
                source = testCase.createBarotropicQGForcingModel( ...
                    sourcePath,forcingCases{iCase,2});
                source.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                source.closeNetCDFFile();
                copyfile(sourcePath,runtimePath)
                copyfile(sourcePath,controlPath)
                ncwriteatt(runtimePath,"/","portableFileIdentifier",caseName);
                ncwriteatt(controlPath,"/","portableFileIdentifier",caseName);
                WVModel.writePortableRunRequest(requestPath,runtimePath, ...
                    schemaVersion=2,method="fixed-rk4",finalTime=.01,initialStep=.0025, ...
                    outputPolicy="append",fftProvider="reference",threads=1,reportPath=reportName);
                immutableRequest = fileread(requestPath);
                [status,output] = systemWithoutMatlabRuntime( ...
                    shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath));
                testCase.assertEqual(status,0,output)
                testCase.verifyEqual(fileread(requestPath),immutableRequest)
                report = jsondecode(fileread(fullfile(testCase.TemporaryFolder,reportName)));
                testCase.verifyEqual(string(report.status),"complete",caseName)
                testCase.verifyTrue(report.execution.noFallback,caseName)
                testCase.verifyEqual(string(report.execution.engine),"reference-direct",caseName)
                testCase.verifyEqual(report.storageBytes.persistentFullHermitian,0,caseName)

                runtimeModel = WVModel.modelFromFile(char(runtimePath));
                runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                controlModel.setupIntegrator(integratorType="fixed",deltaT=.0025);
                controlModel.integrateToTime(.01,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.verifyBarotropicQGModelsEqual( ...
                    runtimeModel,controlModel,fieldNames,caseName);
                testCase.verifyForcingContractsEqual(runtimeModel.wvt,controlModel.wvt,caseName)
                testCase.verifyModelGraphsEqual(runtimeModel,controlModel)
                runtimeModel.closeNetCDFFile();
                controlModel.closeNetCDFFile();
                clear runtimeCleanup controlCleanup
                testCase.verifyNetCDFSchemasEqual(runtimePath,controlPath,caseName)
            end
        end

        function matlabWriterBarotropicQGRejectionsAreTransactional(testCase)
            mooringPath = fullfile(testCase.TemporaryFolder,"qg-reject-mooring.nc");
            mooringModel = testCase.createBarotropicQGModel(mooringPath,[6 5],1,true);
            mooringModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            mooringModel.closeNetCDFFile();
            testCase.verifyRejectedWriterRequest(mooringPath, ...
                @()ncwriteatt(mooringPath, ...
                "/wave-vortex/observingSystems/observingSystems-1", ...
                "AnnotatedClass","WVMooring"),"mooring")

            verticalPath = fullfile(testCase.TemporaryFolder,"qg-reject-vertical-particles.nc");
            verticalModel = testCase.createBarotropicQGModel(verticalPath,[6 5],1,true);
            verticalModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            verticalModel.closeNetCDFFile();
            particleMetadata = "/wave-vortex/observingSystems/observingSystems-3/isXYOnly";
            testCase.verifyRejectedWriterRequest(verticalPath, ...
                @()ncwrite(verticalPath,particleMetadata,uint8(0)),"vertical-particles")

            tracerPath = fullfile(testCase.TemporaryFolder,"qg-reject-rank3-tracer.nc");
            tracerModel = testCase.createBarotropicQGModel(tracerPath,[6 5],1,true);
            tracerModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            tracerModel.closeNetCDFFile();
            tracerMetadata = "/wave-vortex/observingSystems/observingSystems-4/isXYOnly";
            testCase.verifyRejectedWriterRequest(tracerPath, ...
                @()ncwrite(tracerPath,tracerMetadata,uint8(0)),"rank3-tracer")

            forcingPath = fullfile(testCase.TemporaryFolder,"qg-reject-forcing.nc");
            forcingModel = testCase.createBarotropicQGModel(forcingPath,[6 5],1,true);
            forcingModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            forcingModel.closeNetCDFFile();
            testCase.verifyRejectedWriterRequest(forcingPath, ...
                @()ncwriteatt(forcingPath,"/forcing","AnnotatedClass","WVThermalDamping"), ...
                "unsupported-forcing")

            contractPath = fullfile(testCase.TemporaryFolder,"qg-reject-contract.nc");
            contractModel = testCase.createBarotropicQGModel(contractPath,[6 5],1,true);
            contractModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            contractModel.closeNetCDFFile();
            testCase.verifyRejectedWriterRequest(contractPath, ...
                @()ncwriteatt(contractPath,"/","model_version","5.0.0"),"contract-version")

            primaryPath = fullfile(testCase.TemporaryFolder,"qg-reject-graph-primary.nc");
            secondaryPath = fullfile(testCase.TemporaryFolder,"qg-reject-graph-secondary.nc");
            graphModel = testCase.createBarotropicQGModel(primaryPath,[6 5],1,true);
            graphModel.createNetCDFFileForModelOutput(secondaryPath, ...
                outputInterval=.005,shouldOverwriteExisting=true);
            graphModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            graphModel.closeNetCDFFile();
            testCase.verifyRejectedWriterRequest([primaryPath secondaryPath], ...
                @()ncwrite(secondaryPath, ...
                "/wave-vortex/observingSystems/observingSystems-2/absTolerance",2e-6), ...
                "incompatible-graph")
        end

        function matlabWriterBarotropicQGMultiFileGroupsAndPolicies(testCase)
            primaryFixture = fullfile(testCase.TemporaryFolder,"qg-multifile-primary.nc");
            secondaryFixture = fullfile(testCase.TemporaryFolder,"qg-multifile-secondary.nc");
            sourceModel = testCase.createBarotropicQGModel(primaryFixture,[6 5],1,true);
            primaryFile = sourceModel.outputFiles(1);
            denseGroup = primaryFile.addNewEvenlySpacedOutputGroup( ...
                "dense",outputInterval=.0025,initialTime=0,finalTime=.01);
            denseGroup.addObservingSystem(WVEulerianFields(sourceModel, ...
                fieldNames={'u','qgpv'}));
            sourceModel.createNetCDFFileForModelOutput(secondaryFixture, ...
                outputInterval=.005,shouldOverwriteExisting=true);
            sourceModel.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            sourceModel.closeNetCDFFile();
            ncwriteatt(primaryFixture,"/","portableFileIdentifier","primary");
            ncwriteatt(secondaryFixture,"/","portableFileIdentifier","secondary");
            fieldNames = ["u" "v" "eta" "pi" "psi" "qgpv" "zeta_z" "ssh"];

            for policy = ["create" "replace" "append"]
                stem = "qg-multifile-"+policy;
                sourcePrimary = fullfile(testCase.TemporaryFolder,stem+"-source-primary.nc");
                sourceSecondary = fullfile(testCase.TemporaryFolder,stem+"-source-secondary.nc");
                controlPath = fullfile(testCase.TemporaryFolder,stem+"-control.nc");
                copyfile(primaryFixture,sourcePrimary)
                copyfile(secondaryFixture,sourceSecondary)
                copyfile(primaryFixture,controlPath)
                modelFiles = [sourcePrimary sourceSecondary];
                outputPaths = modelFiles;
                destinations = configureDictionary("string","string");
                if policy ~= "append"
                    outputPaths = [ ...
                        fullfile(testCase.TemporaryFolder,stem+"-output-primary.nc") ...
                        fullfile(testCase.TemporaryFolder,stem+"-output-secondary.nc")];
                    destinations(["primary","secondary"]) = outputPaths;
                    if policy == "replace"
                        copyfile(primaryFixture,outputPaths(1))
                        copyfile(secondaryFixture,outputPaths(2))
                    end
                end
                requestPath = fullfile(testCase.TemporaryFolder,stem+"-request.json");
                reportName = stem+"-report.json";
                WVModel.writePortableRunRequest(requestPath,modelFiles, ...
                    schemaVersion=2,method="fixed-rk4",finalTime=.01,initialStep=.0025, ...
                    outputPolicy=policy,destinations=destinations,fftProvider="reference", ...
                    threads=1,reportPath=reportName);
                immutableRequest = testCase.fileBytes(requestPath);
                [status,output] = systemWithoutMatlabRuntime( ...
                    shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath));
                testCase.assertEqual(status,0,output)
                testCase.verifyEqual(testCase.fileBytes(requestPath),immutableRequest,stem)
                report = jsondecode(fileread(fullfile(testCase.TemporaryFolder,reportName)));
                testCase.verifyEqual(string(report.status),"complete",stem)
                testCase.verifyEqual(string(report.outputPolicy),policy,stem)
                testCase.verifyTrue(report.execution.noFallback,stem)
                testCase.verifyTrue(all(isfile(outputPaths)),stem)
                primaryInformation = ncinfo(outputPaths(1));
                secondaryInformation = ncinfo(outputPaths(2));
                testCase.verifyEqual(sort(string({primaryInformation.Groups.Name})), ...
                    sort(["dense" "forcing" "wave-vortex"]),stem)
                testCase.verifyEqual(sort(string({secondaryInformation.Groups.Name})), ...
                    sort(["forcing" "wave-vortex"]),stem)

                runtimeModel = WVModel.modelFromFile(char(outputPaths(1)));
                runtimeCleanup = onCleanup(@()runtimeModel.closeNetCDFFile());
                controlModel = WVModel.modelFromFile(char(controlPath));
                controlCleanup = onCleanup(@()controlModel.closeNetCDFFile());
                controlModel.setupIntegrator(integratorType="fixed",deltaT=.0025);
                controlModel.integrateToTime(.01,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.verifyBarotropicQGModelsEqual(runtimeModel,controlModel,fieldNames,stem)
                runtimeModel.closeNetCDFFile();
                controlModel.closeNetCDFFile();
                clear runtimeCleanup controlCleanup
                testCase.verifyNetCDFSchemasEqual(outputPaths(1),controlPath,stem+" primary")
                testCase.verifyNetCDFSchemasEqual(outputPaths(2),secondaryFixture,stem+" secondary")
            end
        end
    end

    methods (Access = private)
        function model = createNonlinearModel(testCase,path)
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5], ...
                N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
            wvt.Ap(2) = 2e-6 + 1e-6i;
            wvt.Am(3) = -1e-6 + 0.5e-6i;
            wvt.A0(4) = 0.75e-6i;
            [x,y] = ndgrid(wvt.x,wvt.y);
            wvt.addForcing([ ...
                WVAdaptiveDamping(wvt) ...
                WVBottomFrictionLinear(wvt,r=2.5e-7) ...
                WVBottomFrictionQuadratic(wvt,Cd=1e-3) ...
                WVBetaPlanePVAdvection(wvt) ...
                WVPseudoTopographicWaveGeneration(wvt, ...
                    topographicHeight=2*cos(2*pi*x/wvt.Lx)+sin(2*pi*y/wvt.Ly), ...
                    barotropicVelocityAmplitude=[1e-4; -0.5e-4], ...
                    frequency=1.31e-4,name="terrain") ...
                WVFixedAmplitudeForcing(wvt,name="fixed", ...
                    A0_indices=uint64(4),A0bar=0.75e-6i)]);
            model = WVModel(wvt);
            outputFile = model.createNetCDFFileForModelOutput(path, ...
                outputInterval=1e-4,shouldOverwriteExisting=true);
            model.eulerianObservingSystem.addNetCDFOutputVariables('u');
            model.setFloatPositions([500 1500],[400 1200],[-250 -750],'u', ...
                absToleranceXY=1e-8,absToleranceZ=1e-8);
            tracer = sin(2*pi*wvt.X/wvt.Lx);
            model.addTracer(tracer,"dye");
            group = outputFile.outputGroupWithName(model.defaultOutputGroupName());
            particles = model.fluxedObservingSystemWithName("float");
            tracerObserver = model.fluxedObservingSystemWithName("dye");
            group.removeObservingSystem([particles tracerObserver]);
            particleGroup = outputFile.addNewEvenlySpacedOutputGroup("particles",outputInterval=1e-4);
            particleGroup.addObservingSystem(particles);
            tracerGroup = outputFile.addNewEvenlySpacedOutputGroup("tracers",outputInterval=1e-4);
            tracerGroup.addObservingSystem(tracerObserver);
            group.addObservingSystem(WVMooring(model,name="mooring", ...
                x=[0 1000],y=[0 900],trackedFieldNames={'u'}));
            model.setupIntegrator(integratorType="fixed",deltaT=1e-4);
            forcingNames = wvt.forcingNames;
            testCase.assertEqual(string(forcingNames(:)), ...
                ["nonlinear advection"; "linear bottom friction"; "quadratic bottom friction"; "adaptive damping"; ...
                 "beta-plane advection of qgpv"; "terrain"; "fixed"])
        end

        function model = createLinearBottomFrictionModel(testCase,path,isHydrostatic,Nz,r,integratorType)
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 Nz], ...
                N0=sqrt(2e-5),latitude=45,isHydrostatic=isHydrostatic,shouldAntialias=false);
            wvt.Ap(2) = 2e-6 + 1e-6i;
            wvt.Am(3) = -1e-6 + 0.5e-6i;
            wvt.A0(4) = 0.75e-6i;
            wvt.addForcing(WVBottomFrictionLinear(wvt,r=r));
            model = WVModel(wvt);
            model.createNetCDFFileForModelOutput(path,outputInterval=1e-4,shouldOverwriteExisting=true);
            testCase.configureLinearIntegrator(model,integratorType);
        end

        function model = createBarotropicQGModel(~,path,Nxy,j,shouldAntialias,outputInterval)
            if nargin < 6
                outputInterval = .005;
            end
            wvt = WVTransformBarotropicQG([15000 9000],Nxy,h=0.8,j=j, ...
                g=9.80665,planetaryRadius=6.3712e6,rotationRate=7.292115e-5, ...
                latitude=33,shouldAntialias=shouldAntialias);
            index = reshape(1:numel(wvt.A0),[],1);
            wvt.A0(:) = complex(2e-5*sin(0.37*index),1e-5*cos(0.23*(index+2)));
            wvt.A0(wvt.Kh == 0) = 0;
            selfConjugate = wvt.dftPrimaryIndices2D == wvt.dftConjugateIndices2D;
            wvt.A0(selfConjugate) = real(wvt.A0(selfConjugate));
            wvt.addForcing([WVAdaptiveDamping(wvt) ...
                WVBottomFrictionLinear(wvt,r=2.5e-7) ...
                WVBetaPlanePVAdvection(wvt)]);
            model = WVModel(wvt);
            model.eulerianObservingSystem.addNetCDFOutputVariables( ...
                'u','v','eta','pi','psi','qgpv','zeta_z','ssh');
            model.setDrifterPositions([0.13 0.77]*wvt.Lx,[0.21 0.63]*wvt.Ly, ...
                [],'u','qgpv',absToleranceXY=1e-5, ...
                advectionInterpolation="linear",trackedVarInterpolation="spline");
            tracer = 0.3 + 0.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly);
            model.addTracer(tracer,"dye");
            model.createNetCDFFileForModelOutput(char(path), ...
                outputInterval=outputInterval,shouldOverwriteExisting=true);
            model.setupIntegrator(integratorType="fixed",deltaT=0.0025);
        end

        function model = createBarotropicQGForcingModel(testCase,path,forcingKinds)
            wvt = WVTransformBarotropicQG([15000 9000],[6 5],h=.8,j=1, ...
                g=9.80665,planetaryRadius=6.3712e6,rotationRate=7.292115e-5, ...
                latitude=33,shouldAntialias=true);
            index = reshape(1:numel(wvt.A0),[],1);
            wvt.A0(:) = complex(2e-5*sin(.37*index),1e-5*cos(.23*(index+2)));
            wvt.A0(wvt.Kh == 0) = 0;
            selfConjugate = wvt.dftPrimaryIndices2D == wvt.dftConjugateIndices2D;
            wvt.A0(selfConjugate) = real(wvt.A0(selfConjugate));
            for forcingKind = forcingKinds
                switch forcingKind
                    case "nonlinear"
                        forcing = WVNonlinearAdvection(wvt);
                    case "adaptive"
                        forcing = WVAdaptiveDamping(wvt);
                    case "fixed"
                        selected = uint64([2;4]);
                        forcing = WVFixedAmplitudeForcing(wvt,name="fixed", ...
                            A0_indices=selected,A0bar=wvt.A0(selected));
                    case "narrow"
                        selected = uint64([2;4]);
                        forcing = WVNarrowBandGeostrophicForcing(wvt,name="narrow", ...
                            r=2.5e-7,k_r=2*wvt.dk,A0_indices=selected,A0bar=wvt.A0(selected));
                    case "linear"
                        forcing = WVBottomFrictionLinear(wvt,r=2.5e-7);
                    case "quadratic"
                        forcing = WVBottomFrictionQuadratic(wvt,Cd=1.5e-3);
                    case "beta"
                        forcing = WVBetaPlanePVAdvection(wvt);
                    otherwise
                        testCase.assertFail("Unknown QG qualification forcing "+forcingKind)
                end
                wvt.addForcing(forcing);
            end
            model = WVModel(wvt);
            model.eulerianObservingSystem.addNetCDFOutputVariables( ...
                'u','v','eta','pi','psi','qgpv','zeta_z','ssh');
            model.setDrifterPositions([.13 .77]*wvt.Lx,[.21 .63]*wvt.Ly, ...
                [],'u','qgpv',absToleranceXY=1e-5, ...
                advectionInterpolation="linear",trackedVarInterpolation="spline");
            tracer = .3+.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly);
            model.addTracer(tracer,"dye");
            model.createNetCDFFileForModelOutput(char(path), ...
                outputInterval=.005,shouldOverwriteExisting=true);
            model.setupIntegrator(integratorType="fixed",deltaT=.0025);
        end

        function configureLinearIntegrator(~,model,integratorType)
            if integratorType == "ode23"
                model.setupIntegrator(integratorType="adaptive",integrator=@ode23,absTolerance=1e-10,relTolerance=1e-8);
            elseif integratorType == "ode45"
                model.setupIntegrator(integratorType="adaptive",integrator=@ode45,absTolerance=1e-10,relTolerance=1e-8);
            elseif integratorType == "ode78"
                model.setupIntegrator(integratorType="adaptive",integrator=@ode78,absTolerance=1e-10,relTolerance=1e-8);
            else
                model.setupIntegrator(integratorType="fixed",deltaT=1e-5);
            end
        end

        function setConstantCFLFixtureState(~,wvt)
            p = reshape(1:numel(wvt.Ap),[],1);
            wvt.Ap(:) = complex(1e-5*sin(0.17*p),1e-5*cos(0.13*(p+1)));
            wvt.Am(:) = complex(1e-5*cos(0.11*(p+2)),1e-5*sin(0.07*(p+3)));
            wvt.A0(:) = complex(1e-5*sin(0.19*(p+4)),1e-5*cos(0.05*(p+5)));
            wvt.t0 = -0.25;
            wvt.t = 0.5;
        end

        function candidates = matlabCFLCandidates(~,wvt,cfl)
            candidates = struct;
            candidates.effectiveHorizontalGridResolution = wvt.effectiveHorizontalGridResolution;
            candidates.maximumHorizontalSpeed = wvt.uvMax;
            candidates.horizontalAdvective = cfl*candidates.effectiveHorizontalGridResolution/candidates.maximumHorizontalSpeed;
            candidates.verticalAdvective = Inf;
            if wvt.hasVariableWithName("w")
                W = wvt.w;
                W = W(:,:,2:end);
                ratio = abs(W./((3/2)*shiftdim(diff(wvt.z),-2)));
                candidates.verticalAdvective = cfl/max(ratio(:));
            end
            candidates.advective = min(candidates.horizontalAdvective,candidates.verticalAdvective);
            candidates.highestActiveWaveFrequency = 0;
            candidates.oscillatory = Inf;
            if wvt.hasWaveComponent
                omega = wvt.Omega;
                omega(1,:) = [];
                candidates.highestActiveWaveFrequency = max(abs(omega(:)));
                candidates.oscillatory = cfl*2*pi/candidates.highestActiveWaveFrequency;
            end
        end

        function report = runV2Request(testCase,modelPath,integration,reportName)
            reportPath = fullfile(testCase.TemporaryFolder,reportName);
            requestPath = reportPath+".request.json";
            request = struct;
            request.schemaIdentifier = "wave-vortex-run-request-v2";
            request.schemaVersion = 2;
            request.modelFiles = {char(modelPath)};
            request.integration = integration;
            request.output = struct(policy="append",destinations=struct);
            request.execution = struct(fftProvider="reference",threads=1);
            request.report = reportPath;
            writelines(jsonencode(request),requestPath)
            immutableRequest = fileread(requestPath);
            command = shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath);
            [status,output] = systemWithoutMatlabRuntime(command);
            testCase.assertEqual(status,0,output)
            testCase.verifyEqual(fileread(requestPath),immutableRequest)
            report = jsondecode(fileread(reportPath));
            testCase.verifyEqual(string(report.request.schemaIdentifier),"wave-vortex-run-request-v2")
            testCase.verifyEqual(report.request.schemaVersion,2)
            testCase.verifyTrue(report.execution.noFallback)
        end

        function verifyRejectedWriterRequest(testCase,sourcePaths,mutation,label)
            sourcePaths = reshape(string(sourcePaths),1,[]);
            destinations = configureDictionary("string","string");
            destinationPaths = strings(size(sourcePaths));
            for iSource = 1:numel(sourcePaths)
                identifier = label+"-"+iSource;
                ncwriteatt(sourcePaths(iSource),"/","portableFileIdentifier",identifier);
                destinationPaths(iSource) = fullfile(testCase.TemporaryFolder, ...
                    label+"-destination-"+iSource+".nc");
                writelines("qualification sentinel "+label+" "+iSource,destinationPaths(iSource))
                destinations(identifier) = destinationPaths(iSource);
            end
            requestPath = fullfile(testCase.TemporaryFolder,label+"-rejection-request.json");
            WVModel.writePortableRunRequest(requestPath,sourcePaths, ...
                schemaVersion=2,method="fixed-rk4",finalTime=.01,initialStep=.0025, ...
                outputPolicy="replace",destinations=destinations, ...
                fftProvider="reference",threads=1,reportPath=label+"-rejection-report.json");
            immutableRequest = testCase.fileBytes(requestPath);
            mutation();
            sourceBytes = cellfun(@(path)testCase.fileBytes(path),cellstr(sourcePaths), ...
                UniformOutput=false);
            destinationBytes = cellfun(@(path)testCase.fileBytes(path),cellstr(destinationPaths), ...
                UniformOutput=false);
            [status,output] = systemWithoutMatlabRuntime( ...
                shellQuote(testCase.Runner)+" --request "+shellQuote(requestPath));
            testCase.assertNotEqual(status,0,label)
            testCase.verifyEqual(testCase.fileBytes(requestPath),immutableRequest,label)
            for iSource = 1:numel(sourcePaths)
                testCase.verifyEqual(testCase.fileBytes(sourcePaths(iSource)), ...
                    sourceBytes{iSource},label+" source")
                testCase.verifyEqual(testCase.fileBytes(destinationPaths(iSource)), ...
                    destinationBytes{iSource},label+" destination")
            end
            failure = jsondecode(output);
            testCase.verifyEqual(string(failure.status),"failed",label)
            testCase.verifySubstring(string(failure.failure.stage),"preflight",label)
        end

        function bytes = fileBytes(~,path)
            file = fopen(path,"rb");
            assert(file >= 0,"Unable to open %s.",path)
            cleanup = onCleanup(@()fclose(file));
            bytes = fread(file,Inf,"*uint8");
            clear cleanup
        end

        function errorValue = normalizedCoefficientError(~,actual,reference)
            absoluteError = max([max(abs(actual.Ap-reference.Ap),[],"all"),max(abs(actual.Am-reference.Am),[],"all"),max(abs(actual.A0-reference.A0),[],"all")]);
            referenceScale = max([max(abs(reference.Ap),[],"all"),max(abs(reference.Am),[],"all"),max(abs(reference.A0),[],"all"),eps]);
            errorValue = absoluteError/referenceScale;
        end

        function [coefficientError,fieldError,particleError,tracerError] = ...
                verifyBarotropicQGModelsEqual(testCase,actual,expected,fieldNames,diagnostic)
            coefficientScale = max(max(abs(expected.wvt.A0),[],'all'),eps);
            coefficientError = max(abs(actual.wvt.A0-expected.wvt.A0),[],'all')/coefficientScale;
            testCase.verifyLessThanOrEqual(coefficientError,1e-12,diagnostic+" A0")
            fieldError = 0;
            for fieldName = fieldNames
                actualField = actual.wvt.(fieldName);
                expectedField = expected.wvt.(fieldName);
                scale = max(max(abs(expectedField),[],'all'),eps);
                fieldError = max(fieldError,max(abs(actualField-expectedField),[],'all')/scale);
            end
            testCase.verifyLessThanOrEqual(fieldError,1e-12,diagnostic+" fields")
            energyScale = max(abs(expected.wvt.totalEnergy),eps);
            testCase.verifyLessThanOrEqual(abs(actual.wvt.totalEnergy-expected.wvt.totalEnergy)/energyScale, ...
                1e-12,diagnostic+" energy")
            enstrophyScale = max(abs(expected.wvt.totalEnstrophy()),eps);
            testCase.verifyLessThanOrEqual(abs(actual.wvt.totalEnstrophy()-expected.wvt.totalEnstrophy())/enstrophyScale, ...
                1e-12,diagnostic+" enstrophy")
            [actualX,actualY,actualZ,actualTracked] = actual.drifterPositions();
            [expectedX,expectedY,expectedZ,expectedTracked] = expected.drifterPositions();
            particleError = max([max(abs(actualX-expectedX),[],'all'), ...
                max(abs(actualY-expectedY),[],'all'), ...
                max(abs(actualTracked.u-expectedTracked.u),[],'all'), ...
                max(abs(actualTracked.qgpv-expectedTracked.qgpv),[],'all')]);
            testCase.verifyEqual(actualZ,expectedZ,diagnostic+" particle z")
            testCase.verifyLessThanOrEqual(particleError,1e-8,diagnostic+" particles")
            actualTracer = actual.tracer("dye");
            expectedTracer = expected.tracer("dye");
            tracerError = max(abs(actualTracer-expectedTracer),[],'all');
            testCase.verifyLessThanOrEqual(tracerError,1e-8,diagnostic+" tracer")
            testCase.verifyLessThanOrEqual(abs(actual.wvt.t-expected.wvt.t), ...
                4*eps(expected.wvt.t),diagnostic)
        end

        function verifyModelGraphsEqual(testCase,actual,expected)
            testCase.verifyEqual(actual.isDynamicsLinear,expected.isDynamicsLinear)
            actualForcingNames = actual.wvt.forcingNames;
            expectedForcingNames = expected.wvt.forcingNames;
            testCase.verifyEqual(string(actualForcingNames(:)),string(expectedForcingNames(:)))
            actualFile = actual.outputFiles(1);
            expectedFile = expected.outputFiles(1);
            testCase.verifyEqual(sort(actualFile.outputGroupNames),sort(expectedFile.outputGroupNames))
            for groupName = reshape(expectedFile.outputGroupNames,1,[])
                actualGroup = actualFile.outputGroupWithName(groupName);
                expectedGroup = expectedFile.outputGroupWithName(groupName);
                testCase.verifyEqual(string(arrayfun(@class,actualGroup.observingSystems,UniformOutput=false)), ...
                    string(arrayfun(@class,expectedGroup.observingSystems,UniformOutput=false)))
                testCase.verifyEqual(string({actualGroup.observingSystems.name}), ...
                    string({expectedGroup.observingSystems.name}))
                testCase.verifyEqual(actualGroup.outputInterval,expectedGroup.outputInterval)
                testCase.verifyEqual(actualGroup.initialTime,expectedGroup.initialTime)
                testCase.verifyEqual(actualGroup.finalTime,expectedGroup.finalTime)
                testCase.verifyEqual(actualGroup.incrementsWrittenToGroup, ...
                    expectedGroup.incrementsWrittenToGroup)
            end
        end

        function verifyForcingContractsEqual(testCase,actual,expected,diagnostic)
            testCase.verifyEqual(numel(actual.forcing),numel(expected.forcing),diagnostic)
            for iForcing = 1:numel(expected.forcing)
                testCase.verifyEqual(class(actual.forcing(iForcing)), ...
                    class(expected.forcing(iForcing)),diagnostic)
                testCase.verifyEqual(actual.forcing(iForcing).portableImplementationContract(), ...
                    expected.forcing(iForcing).portableImplementationContract(),diagnostic)
                if isa(expected.forcing(iForcing),"WVNarrowBandGeostrophicForcing")
                    for propertyName = ["r" "k_r" "k_f" "j_f" "u_rms" "initialPV"]
                        testCase.verifyEqual(actual.forcing(iForcing).(propertyName), ...
                            expected.forcing(iForcing).(propertyName),diagnostic+" "+propertyName)
                    end
                end
            end
        end

        function verifyBarotropicQGOutputSeriesEqual(testCase,actualPath,expectedPath,fieldNames,diagnostic)
            spectralAndFields = ["A0_real" "A0_imag" fieldNames];
            for variableName = spectralAndFields
                actual = ncread(actualPath,"/wave-vortex/"+variableName);
                expected = ncread(expectedPath,"/wave-vortex/"+variableName);
                scale = max(max(abs(expected),[],"all"),eps);
                relativeError = max(abs(actual-expected),[],"all")/scale;
                testCase.verifyLessThanOrEqual(relativeError,1e-12, ...
                    diagnostic+" dense "+variableName)
            end
            for variableName = ["drifter_x" "drifter_y" "dye"]
                actual = ncread(actualPath,"/wave-vortex/"+variableName);
                expected = ncread(expectedPath,"/wave-vortex/"+variableName);
                testCase.verifyLessThanOrEqual(max(abs(actual-expected),[],"all"),1e-8, ...
                    diagnostic+" dense "+variableName)
            end
            actualTime = ncread(actualPath,"/wave-vortex/t");
            expectedTime = ncread(expectedPath,"/wave-vortex/t");
            testCase.verifyEqual(actualTime,expectedTime,diagnostic+" output times", ...
                AbsTol=4*eps(max(expectedTime)))
        end

        function verifyNetCDFSchemasEqual(testCase,actualPath,expectedPath,diagnostic)
            testCase.verifyEqual(netCDFSchemaSignature(actualPath), ...
                netCDFSchemaSignature(expectedPath),diagnostic+" NetCDF schema")
        end
    end
end

function value = shellQuote(value)
value = "'" + replace(string(value),"'","'""'""'") + "'";
end

function [status,output] = systemWithoutMatlabRuntime(command)
if isunix && ~ismac
    command = "env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH " + string(command);
end
[status,output] = system(command);
end

function signature = netCDFSchemaSignature(path)
signature = netCDFGroupSchema(ncinfo(path));
end

function signature = netCDFGroupSchema(information)
signature = struct;
signature.Name = string(information.Name);
signature.Dimensions = strings(0,1);
if ~isempty(information.Dimensions)
    dimensionNames = string({information.Dimensions.Name});
    [dimensionNames,order] = sort(dimensionNames);
    dimensionSignatures = strings(numel(order),1);
    for iDimension = 1:numel(order)
        dimension = information.Dimensions(order(iDimension));
        dimensionSignatures(iDimension) = dimensionNames(iDimension)+":"+string(dimension.Unlimited);
    end
    signature.Dimensions = dimensionSignatures;
end
signature.Attributes = netCDFAttributeSchema(information.Attributes);
signature.Variables = cell(0,1);
if ~isempty(information.Variables)
    [~,order] = sort(string({information.Variables.Name}));
    variableSignatures = cell(numel(order),1);
    for iVariable = 1:numel(order)
        variable = information.Variables(order(iVariable));
        value = struct;
        value.Name = string(variable.Name);
        value.Datatype = string(variable.Datatype);
        value.Dimensions = strings(1,0);
        if ~isempty(variable.Dimensions)
            value.Dimensions = string({variable.Dimensions.Name});
        end
        value.Attributes = netCDFAttributeSchema(variable.Attributes);
        variableSignatures{iVariable} = value;
    end
    signature.Variables = variableSignatures;
end
signature.Groups = cell(0,1);
if ~isempty(information.Groups)
    groupNames = string({information.Groups.Name});
    information.Groups(groupNames == "observingSystems" | ...
        endsWith(groupNames,"/observingSystems")) = [];
    [~,order] = sort(string({information.Groups.Name}));
    groupSignatures = cell(numel(order),1);
    for iGroup = 1:numel(order)
        groupSignatures{iGroup} = netCDFGroupSchema(information.Groups(order(iGroup)));
    end
    signature.Groups = groupSignatures;
end
end

function signature = netCDFAttributeSchema(attributes)
signature = strings(0,1);
if isempty(attributes)
    return
end
schemaAttributes = [ ...
    "AnnotatedClass" "WVModelIsDynamicsLinear" "WVTransform" "model_version" ...
    "name" "fieldNames" "trackedFieldNames" "advectionInterpolation" ...
    "trackedVarInterpolation" "initialPV" "darwinSymbol" ...
    "isComplex" "isRealPart" "isImaginaryPart" "isLogicalType" ...
    "isParticle" "particleName" "particleVariableName" "isTracer"];
attributes = attributes(ismember(string({attributes.Name}),schemaAttributes));
if isempty(attributes)
    return
end
[names,~] = sort(string({attributes.Name}));
signature = reshape(names,[],1);
end
