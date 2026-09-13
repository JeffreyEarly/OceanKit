classdef TestPortableStableForcing < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        runner (1,1) string
        executable (1,1) string
        providers (1,:) string
    end
    properties (TestParameter)
        outputDynamics = struct(linear=true,nonlinear=false)
        outputFamily = struct(constantHydrostatic="constant-hydrostatic",constantNonhydrostatic="constant-nonhydrostatic",barotropic="barotropic",stratifiedQG="stratified-qg",hydrostatic="hydrostatic",boussinesq="boussinesq")
    end
    methods (TestClassSetup)
        function buildProbe(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture(PreservingOnFailure=true));
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_STABLE_FORCING_DUMP"));
            testCase.runner = string(getenv("WV_STABLE_FORCING_RUNNER"));
            testCase.providers = "reference";
            if getenv("WV_STABLE_FORCING_NATIVE") == "1", testCase.providers = ["reference","native"]; end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 4 --target WVStableForcingDump wave-vortex-run");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVStableForcingDump");
                testCase.runner = fullfile(build,"wave-vortex-run");
            end
            testCase.assertTrue(isfile(testCase.executable));
            if testCase.runner == "", testCase.runner = fullfile(fileparts(testCase.executable),"wave-vortex-run"); end
            testCase.assertTrue(isfile(testCase.runner));
        end
    end
    methods (Test,TestTags="full")
        function forcingOutputContinuationMatchesMatlab(testCase,outputFamily,outputDynamics)
            family = string(outputFamily);
            if family == "barotropic"
                wvt = WVTransformBarotropicQG([17000 11000],[8 6],j=1,shouldAntialias=false);
            elseif family == "stratified-qg"
                wvt = WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
            else
                wvt = diagnosticWaveTransform(family,[8 6 9],false);
            end
            isQG = family == "barotropic" || family == "stratified-qg";
            if isQG
                n = reshape(1:numel(wvt.A0),size(wvt.A0));
                wvt.A0 = 1e-6*complex(sin(.17*n),cos(.23*n)).*(wvt.Kh>0);
                wvt.t = 37;
            end
            forces = [WVNonlinearAdvection(wvt),WVBottomFrictionLinear(wvt,r=2.5e-7),WVBottomFrictionQuadratic(wvt,Cd=.002)];
            fixed = testCase.forcing(wvt,"WVFixedAmplitudeForcing");
            if family == "stratified-qg"
                fixedIndex = find(wvt.Kh>0 & wvt.J==1,1);
                fixed = WVFixedAmplitudeForcing(wvt,name="catalog fixed",A0_indices=uint64(fixedIndex),A0bar=wvt.A0(fixedIndex));
                forces = [forces,WVVerticalDiffusivity(wvt,kappa_z=.002)];
            end
            forces = [forces,WVBetaPlanePVAdvection(wvt),WVAdaptiveDamping(wvt),testCase.forcing(wvt,"WVNarrowBandGeostrophicForcing"),fixed,WVAntialiasing(wvt,Nj=3)];
            if ~isQG
                terrain = 2*cos(2*pi*reshape(wvt.x,[],1)/wvt.Lx)+sin(2*pi*reshape(wvt.y,1,[])/wvt.Ly);
                tidal = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[.01;.005],frequency=1e-4,rampDuration=0);
                forces = [forces,WVHorizontalDamping(wvt,nu=.125,kappa=.002),WVVerticalDamping(wvt,nu=.125,kappa=.002),WVVerticalDiffusivity(wvt,kappa_z=.002),tidal];
            end
            wvt.setForcing(forces);
            operation = SpatialForcingOperation(wvt);
            wvt.addOperation(operation);
            names = string({operation.outputVariables.name});
            model = WVModel(wvt,shouldUseLinearDynamics=outputDynamics);
            model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
            source = fullfile(testCase.folder,"forcing-output-source.nc");
            file = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
            dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
            dense.addObservingSystem(WVEulerianFields(model,fieldNames=cellstr([names,"u","v"])));
            file.outputTimesForIntegrationPeriod(37,38);
            file.writeTimeStepToOutputFile(37);
            model.closeNetCDFFile();
            ncwriteatt(source,"/","portableFileIdentifier","forcing-output-primary");
            baselineSource = fullfile(testCase.folder,"forcing-baseline-source.nc");
            if ~outputDynamics
                baseline = WVModel(wvt,shouldUseLinearDynamics=false);
                baseline.eulerianObservingSystem.addNetCDFOutputVariables('u','v');
                baselineFile = baseline.createNetCDFFileForModelOutput(baselineSource,outputInterval=.5,shouldOverwriteExisting=true);
                baselineDense = baselineFile.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
                baselineDense.addObservingSystem(WVEulerianFields(baseline,fieldNames={'u','v'}));
                baselineFile.outputTimesForIntegrationPeriod(37,38);
                baselineFile.writeTimeStepToOutputFile(37);
                baseline.closeNetCDFFile();
                ncwriteatt(baselineSource,"/","portableFileIdentifier","forcing-baseline-primary");
            end
            controlPath = fullfile(testCase.folder,"forcing-output-matlab.nc");
            copyfile(source,controlPath);
            control = testCase.verifyWarningFree(@()WVModel.modelFromFile(controlPath));
            cleanup = onCleanup(@()control.closeNetCDFFile());
            control.setupIntegrator(integratorType="fixed",deltaT=.25);
            control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            control.closeNetCDFFile();
            clear cleanup
            for provider = testCase.providers
                for scenario = 1:4
                    outputPath = fullfile(testCase.folder,"forcing-output-runtime.nc");
                    copyfile(source,outputPath);
                    baselinePath = fullfile(testCase.folder,"forcing-baseline-runtime.nc");
                    if ~outputDynamics, copyfile(baselineSource,baselinePath); end
                    times = 38;
                    if scenario == 2, times = [37.5 38]; end
                    for finalTime = times
                        report = testCase.runForcingOutputContinuation(outputPath,provider,scenario,finalTime);
                        if ~outputDynamics
                            baselineReport = testCase.runForcingOutputContinuation(baselinePath,provider,scenario,finalTime);
                            testCase.verifyEqual(report.state,baselineReport.state,"Diagnostics changed integration progress or RHS counts.");
                            for metric = ["lastNormalizedError","lastAcceptedStepSize","nextStepSize","constraintModifiedCoefficientCount","baseRightHandSideEvaluationCount"]
                                testCase.verifyEqual(report.integrator.(metric),baselineReport.integrator.(metric),"Diagnostics changed "+metric);
                            end
                            testCase.verifyEqual(report.forcing,baselineReport.forcing,"Diagnostics changed the resolved forcing schedule.");
                            coefficientNames = ["A0_real","A0_imag"];
                            if ~isQG, coefficientNames = ["A0_real","A0_imag","Ap_real","Ap_imag","Am_real","Am_imag"]; end
                            for coefficient = coefficientNames
                                variablePath = "/wave-vortex/"+coefficient;
                                testCase.verifyEqual(ncread(outputPath,variablePath),ncread(baselinePath,variablePath),"Diagnostics changed saved "+coefficient);
                            end
                        end
                        testCase.verifyEqual(report.state.stepCount,4-2*double(scenario==2));
                        testCase.verifyEqual(report.diagnosticEvaluation.workspaceLiveBytes,0);
                        testCase.verifyGreaterThan(report.diagnosticEvaluation.evaluationCount,0);
                        testCase.verifyGreaterThan(report.integrator.denseOutputEvaluationCount,0);
                    end
                    comparisons = 0;
                    maximumError = 0;
                    for groupName = ["wave-vortex","dense"]
                        expectedTimes = ncread(controlPath,"/"+groupName+"/t");
                        actualTimes = ncread(outputPath,"/"+groupName+"/t");
                        testCase.verifyEqual(actualTimes,expectedTimes,family+" "+groupName+" times");
                        testCase.verifyNumElements(actualTimes,3+6*double(groupName=="dense"));
                        for name = names
                            variablePath = "/"+groupName+"/"+name;
                            expected = ncread(controlPath,variablePath);
                            actual = ncread(outputPath,variablePath);
                            testCase.verifySize(actual,size(expected),family+" "+name);
                            % Each record is compared independently so a
                            % large initial field cannot hide a bad dense value.
                            for record = 1:numel(actualTimes)
                                index = repmat({':'},1,ndims(expected)); index{end} = record;
                                reference = expected(index{:}); values = actual(index{:});
                                testCase.assertTrue(all(isfinite(reference(:))),"Nonfinite MATLAB saved values: "+family+" "+variablePath+" shape="+mat2str(size(expected)));
                                testCase.assertTrue(all(isfinite(values(:))),"Nonfinite C++ saved values: "+family+" "+variablePath+" shape="+mat2str(size(actual)));
                                scale = max(abs(reference(:)));
                                if startsWith(name,"Fw_pseudo_topographic")
                                    originalTime = wvt.t; wvt.t = actualTimes(record);
                                    zero = zeros(size(wvt.Ap));
                                    [Fp,Fm,~] = tidal.addSpectralForcing(wvt,zero,zero,zero);
                                    plus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*wvt.phase.*Fp);
                                    minus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAm.*wvt.conjPhase.*Fm);
                                    scale = max(abs([plus(:);minus(:)]));
                                    wvt.t = originalTime;
                                    testCase.assertGreaterThan(scale,0);
                                    testCase.verifyLessThanOrEqual(max(abs(plus(:)+minus(:))),1e-12*scale);
                                    testCase.verifyLessThanOrEqual(max(abs(reference(:))),1e-12*scale);
                                    testCase.verifyLessThanOrEqual(max(abs(values(:))),1e-12*scale);
                                end
                                error = max(abs(values(:)-reference(:)))/max(scale,realmin);
                                testCase.verifyLessThanOrEqual(error,1e-12,family+" "+provider+" scenario="+scenario+" "+groupName+" "+name+" t="+actualTimes(record));
                                maximumError = max(maximumError,error);
                                comparisons = comparisons+1;
                            end
                            for attribute = ["units","long_name"]
                                testCase.verifyEqual(string(ncreadatt(outputPath,variablePath,attribute)),string(ncreadatt(controlPath,variablePath,attribute)));
                            end
                        end
                    end
                    fprintf('FORCING_OUTPUT_DIAGNOSTICS %s provider=%s scenario=%d comparisons=%d max_relative=%.17g linear=%d\n',outputFamily,provider,scenario,comparisons,maximumError,outputDynamics);
                end
            end
        end
        function degenerateAdaptiveDampingMatchesMatlab(testCase)
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for retainedModes = [1 2 3]
                    grid = [8 6 9];
                    if retainedModes == 2, grid = [5 5 9]; end
                    wvt = diagnosticWaveTransform(family,grid,false);
                    damping = WVAdaptiveDamping(wvt);
                    wvt.setForcing([WVNonlinearAdvection(wvt),WVBetaPlanePVAdvection(wvt),WVAntialiasing(wvt,Nj=retainedModes),damping]);
                    testCase.verifyEqual(wvt.effectiveJMax,retainedModes-1);
                    testCase.verifyTrue(all(isfinite(damping.damp),"all"));
                    operation = SpatialForcingOperation(wvt);
                    expected = cell(1,operation.nVarOut);
                    [expected{:}] = operation.compute(wvt);
                    names = string({operation.outputVariables.name});
                    for value = expected, testCase.assertTrue(all(isfinite(value{1}),"all")); end
                    before = {wvt.Ap,wvt.Am,wvt.A0};
                    testCase.compareTendency(wvt);
                    wvt.addOperation(operation);
                    model = WVModel(wvt);
                    model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
                    source = fullfile(testCase.folder,"degenerate-source.nc");
                    file = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
                    file.outputTimesForIntegrationPeriod(wvt.t,wvt.t+1);
                    file.writeTimeStepToOutputFile(wvt.t);
                    model.closeNetCDFFile();
                    for index = 1:numel(names)
                        saved = ncread(source,"/wave-vortex/"+names(index));
                        testCase.verifyEqual(saved(:),expected{index}(:));
                    end
                    restored = testCase.verifyWarningFree(@()WVModel.modelFromFile(source));
                    cleanup = onCleanup(@()restored.closeNetCDFFile());
                    restoredDamping = restored.wvt.forcingWithName("adaptive damping");
                    testCase.verifyEqual(restoredDamping.damp,damping.damp);
                    restored.closeNetCDFFile(); clear cleanup
                    resultPath = fullfile(testCase.folder,"degenerate-tendencies.json");
                    for provider = testCase.providers
                        [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(resultPath)+" "+provider+" tendencies");
                        testCase.assertEqual(status,0,provider+": "+output);
                        actual = jsondecode(fileread(resultPath));
                        testCase.verifyEqual(actual.diagnosticWorkspaceLiveBytes,0);
                        testCase.verifyEqual(actual.diagnosticForcingEvaluationCount,numel(wvt.forcing));
                        comparisons = 0;
                        maximumError = 0;
                        for instance = reshape(actual.tendencies,1,[])
                            suffix = replace(string(instance.name),[" ","-"],"_");
                            for channel = string(fieldnames(instance.fields))'
                                reference = expected{names==channel+"_"+suffix};
                                values = instance.fields.(channel);
                                testCase.verifyTrue(all(isfinite(values(:))));
                                scale = max(abs(reference(:)));
                                testCase.verifyLessThanOrEqual(max(abs(values(:)-reference(:))),1e-12*max(scale,realmin),family+" modes="+retainedModes+" "+provider+" "+channel+" "+suffix);
                                maximumError = max(maximumError,max(abs(values(:)-reference(:)))/max(scale,realmin));
                                comparisons = comparisons+1;
                            end
                        end
                        testCase.verifyEqual(comparisons,operation.nVarOut);
                        fprintf('DEGENERATE_DAMPING %s modes=%d provider=%s comparisons=%d max_relative=%.17g\n',family,retainedModes,provider,comparisons,maximumError);
                    end
                    testCase.verifyEqual({wvt.Ap,wvt.Am,wvt.A0},before);
                end
            end
        end
        function isolatedWaveDiagnosticsHaveSpatialShapes(testCase)
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                wvt = diagnosticWaveTransform(family,[8 6 9],false);
                for force = [WVAdaptiveDamping(wvt),WVVerticalDiffusivity(wvt,kappa_z=.002)]
                    wvt.setForcing(force);
                    before = {wvt.Ap,wvt.Am,wvt.A0};
                    operation = SpatialForcingOperation(wvt);
                    values = cell(1,operation.nVarOut);
                    [values{:}] = operation.compute(wvt);
                    for value = values
                        testCase.verifySize(value{1},wvt.spatialMatrixSize);
                        testCase.verifyTrue(all(isfinite(value{1}),"all"));
                    end
                    testCase.verifyEqual({wvt.Ap,wvt.Am,wvt.A0},before);
                end
            end
        end
        function isolatedQGDiagnosticsPreserveOutputIdentity(testCase)
            wvt = WVTransformBarotropicQG([17000 11000],[8 6],shouldAntialias=false);
            wvt.setForcing(WVAdaptiveDamping(wvt));
            operation = SpatialForcingOperation(wvt);
            testCase.verifyEqual(string(operation.name),string(operation.outputVariables.name));
            testCase.verifyEqual(string(operation.name),"Fqgpv_"+replace(string(wvt.forcing.name),[" ","-"],"_"));
            wvt.addOperation(operation);
            testCase.verifyEqual(wvt.variableWithName(operation.name),zeros(wvt.spatialMatrixSize));
        end
        function isolatedTendenciesMatchMatlab(testCase)
            matrix = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
            for config = reshape(matrix.configurations,1,[])
                family = extractBefore(string(config.id),"-aa");
                if family == "barotropic"
                    wvt = WVTransformBarotropicQG([17000 11000],[8 6],j=1,shouldAntialias=config.shouldAntialias);
                elseif family == "stratified-qg"
                    wvt = WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=config.shouldAntialias);
                else
                    wvt = diagnosticWaveTransform(family,[8 6 9],config.shouldAntialias);
                end
                if family == "barotropic" || family == "stratified-qg"
                    n = reshape(1:numel(wvt.A0),size(wvt.A0));
                    wvt.A0 = 1e-6*complex(sin(.17*n),cos(.23*n)).*(wvt.Kh>0);
                    wvt.t = 37;
                end
                rows = matrix.rows(string({matrix.rows.configuration}) == string(config.id));
                for row = reshape(rows,1,[])
                    if string(row.matlab.applicability) ~= "applicable", continue; end
                    force = testCase.forcing(wvt,string(row.forcing));
                    if string(row.forcing) == "WVPseudoTopographicWaveGeneration"
                        % A two-dimensional terrain avoids an analytically
                        % vanishing Fu channel from one-dimensional symmetry.
                        terrain = 2*cos(2*pi*reshape(wvt.x,[],1)/wvt.Lx)+sin(2*pi*reshape(wvt.y,1,[])/wvt.Ly);
                        force = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[.01;.005],frequency=1e-4,rampDuration=0);
                    end
                    wvt.setForcing(force);
                    before = {wvt.Ap,wvt.Am,wvt.A0};
                    operation = SpatialForcingOperation(wvt);
                    expected = cell(1,operation.nVarOut);
                    [expected{:}] = operation.compute(wvt);
                    names = string({operation.outputVariables.name});
                    expectedChannels = extractBefore(names,"_");
                    cancellationScale = 0;
                    if ~wvt.isHydrostatic && string(row.forcing) == "WVPseudoTopographicWaveGeneration"
                        zero = zeros(size(wvt.Ap));
                        [Fp,Fm,~] = force.addSpectralForcing(wvt,zero,zero,zero);
                        plus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*wvt.phase.*Fp);
                        minus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAm.*wvt.conjPhase.*Fm);
                        cancellationScale = max(abs([plus(:);minus(:)]));
                        testCase.assertGreaterThan(cancellationScale,0);
                        testCase.verifyLessThanOrEqual(max(abs(plus(:)+minus(:))),1e-12*cancellationScale);
                    end
                    source = testCase.writeInitialModel(wvt);
                    resultPath = fullfile(testCase.folder,"isolated-tendencies.json");
                    for provider = testCase.providers
                        [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(resultPath)+" "+provider+" tendencies");
                        testCase.assertEqual(status,0,string(row.id)+" "+provider+": "+output);
                        actual = jsondecode(fileread(resultPath));
                        testCase.assertNumElements(actual.tendencies,1,string(row.id));
                        testCase.verifyEqual(actual.diagnosticWorkspaceLiveBytes,0);
                        testCase.verifyEqual(actual.diagnosticForcingEvaluationCount,1);
                        instance = actual.tendencies;
                        maximumError = 0;
                        channels = string(fieldnames(instance.fields))';
                        testCase.verifyEqual(numel(channels),numel(expectedChannels));
                        for channel = channels
                            name = channel+"_"+replace(string(instance.name),[" ","-"],"_");
                            index = find(names==name);
                            testCase.assertNumElements(index,1,name);
                            reference = expected{index};
                            values = instance.fields.(channel);
                            scale = max(abs(reference(:)));
                            if channel == "Fw" && cancellationScale > 0
                                testCase.verifyLessThanOrEqual(scale,1e-12*cancellationScale);
                                testCase.verifyLessThanOrEqual(max(abs(values(:))),1e-12*cancellationScale);
                                scale = cancellationScale;
                            end
                            error = max(abs(values(:)-reference(:)))/max(scale,realmin);
                            testCase.verifyLessThanOrEqual(error,1e-12,string(row.id)+" "+provider+" "+name);
                            maximumError = max(maximumError,error);
                        end
                        fprintf('ISOLATED_FORCING_DIAGNOSTICS %s provider=%s comparisons=%d max_relative=%.17g\n',row.id,provider,numel(channels),maximumError);
                    end
                    testCase.verifyEqual({wvt.Ap,wvt.Am,wvt.A0},before);
                end
            end
        end
        function fullGridWaveTendenciesMatchMatlab(testCase)
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for antialias = [false true]
                    for grid = {[8 6 9],[9 7 10]}
                        wvt = diagnosticWaveTransform(family,grid{1},antialias);
                        fixed = WVFixedAmplitudeForcing(wvt,name="held-coefficients",Ap_indices=uint64(2),Apbar=wvt.Ap(2),A0_indices=uint64(2),A0bar=wvt.A0(2));
                        terrain = 2*cos(2*pi*reshape(wvt.x,[],1)/wvt.Lx)+sin(2*pi*reshape(wvt.y,1,[])/wvt.Ly);
                        tidal = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[.01;.005],frequency=1e-4,rampDuration=0);
                        forces = [WVNonlinearAdvection(wvt),WVBottomFrictionLinear(wvt,r=2.5e-7),WVBottomFrictionQuadratic(wvt,Cd=.002),WVHorizontalDamping(wvt,nu=.125,kappa=.002),WVVerticalDamping(wvt,nu=.125,kappa=.002),WVVerticalDiffusivity(wvt,kappa_z=.002),WVAdaptiveDamping(wvt),WVBetaPlanePVAdvection(wvt),tidal,testCase.forcing(wvt,"WVNarrowBandGeostrophicForcing"),fixed];
                        if ~antialias, forces = [forces,WVAntialiasing(wvt,Nj=3)]; end %#ok<AGROW>
                        wvt.setForcing(forces);
                        before = {wvt.Ap,wvt.Am,wvt.A0};
                        verticalCancellationScale = 0;
                        if ~wvt.isHydrostatic
                            % Tidal Ap/Am contributions to Fw cancel analytically.
                            % Use their separate physical magnitudes to condition
                            % this zero-channel comparison, retaining the 1e-12
                            % bound without dividing by a roundoff-only residual.
                            zero = zeros(size(wvt.Ap));
                            [tidalFp,tidalFm,~] = tidal.addSpectralForcing(wvt,zero,zero,zero);
                            plus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*wvt.phase.*tidalFp);
                            minus = wvt.transformToSpatialDomainWithG(Apm=wvt.WAm.*wvt.conjPhase.*tidalFm);
                            verticalCancellationScale = max(abs([plus(:);minus(:)]));
                            testCase.assertGreaterThan(verticalCancellationScale,0);
                            testCase.verifyLessThanOrEqual(max(abs(plus(:)+minus(:))),1e-12*verticalCancellationScale);
                        end
                        operation = SpatialForcingOperation(wvt);
                        expected = cell(1,operation.nVarOut);
                        [expected{:}] = operation.compute(wvt);
                        names = string({operation.outputVariables.name});
                        source = testCase.writeInitialModel(wvt);
                        resultPath = fullfile(testCase.folder,"tendencies.json");
                        for provider = testCase.providers
                            [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(resultPath)+" "+provider+" tendencies");
                            testCase.assertEqual(status,0,family+" "+provider+": "+output);
                            actual = jsondecode(fileread(resultPath));
                            testCase.verifyEqual(actual.diagnosticWorkspaceLiveBytes,0);
                            testCase.verifyEqual(actual.diagnosticForcingEvaluationCount,numel(wvt.forcing));
                            comparisons = 0;
                            maximumError = 0;
                            for instance = reshape(actual.tendencies,1,[])
                                suffix = replace(string(instance.name),[" ","-"],"_");
                                for channel = string(fieldnames(instance.fields))'
                                    name = channel+"_"+suffix;
                                    index = find(names==name);
                                    testCase.assertNumElements(index,1,name);
                                    reference = expected{index};
                                    values = instance.fields.(channel);
                                    absoluteError = max(abs(values(:)-reference(:)));
                                    scale = max(abs(reference(:)));
                                    if channel=="Fw" && string(instance.type)=="WVPseudoTopographicWaveGeneration"
                                        testCase.verifyLessThanOrEqual(scale,1e-12*verticalCancellationScale,"MATLAB tidal Fw must remain a cancellation residual");
                                        testCase.verifyLessThanOrEqual(max(abs(values(:))),1e-12*verticalCancellationScale,"C++ tidal Fw must remain a cancellation residual");
                                        scale = verticalCancellationScale;
                                    end
                                    error = absoluteError/max(scale,realmin);
                                    if error > 1e-12
                                        fprintf('FORCING_DIAGNOSTIC_ERROR %s %s abs=%.17g reference=%.17g\n',family,name,absoluteError,scale);
                                    end
                                    testCase.verifyLessThanOrEqual(error,1e-12,family+" aa="+antialias+" "+provider+" "+name);
                                    maximumError = max(maximumError,error);
                                    comparisons = comparisons+1;
                                end
                            end
                            testCase.verifyEqual(comparisons,operation.nVarOut);
                            fprintf('FORCING_DIAGNOSTICS %s aa=%d grid=%s provider=%s comparisons=%d max_relative=%.17g\n',family,antialias,mat2str(grid{1}),provider,comparisons,maximumError);
                        end
                        testCase.verifyEqual({wvt.Ap,wvt.Am,wvt.A0},before);
                    end
                end
            end
        end
        function fullGridBarotropicTendenciesMatchMatlab(testCase)
            for j = [0 1]
                for antialias = [false true]
                    for grid = {[8 6],[9 7]}
                        wvt = WVTransformBarotropicQG([17000 11000],grid{1},j=j,shouldAntialias=antialias);
                        n = reshape(1:numel(wvt.A0),size(wvt.A0));
                        wvt.A0 = 1e-6*complex(sin(.17*n),cos(.23*n)).*(wvt.Kh>0);
                        wvt.t = 37;
                        fixed = WVFixedAmplitudeForcing(wvt,name="held-pv",A0_indices=uint64((1:numel(wvt.A0))'),A0bar=wvt.A0(:));
                        forces = [WVNonlinearAdvection(wvt),WVBottomFrictionLinear(wvt,r=2.5e-7),WVBottomFrictionQuadratic(wvt,Cd=.002),WVBetaPlanePVAdvection(wvt),WVAdaptiveDamping(wvt),WVNarrowBandGeostrophicForcing(wvt,initialPV="none",k_f=2*wvt.dk,j_f=j),fixed];
                        if ~antialias, forces = [forces,WVAntialiasing(wvt)]; end %#ok<AGROW>
                        wvt.setForcing(forces);
                        before = wvt.A0;
                        operation = SpatialForcingOperation(wvt);
                        expected = cell(1,operation.nVarOut);
                        [expected{:}] = operation.compute(wvt);
                        names = string({operation.outputVariables.name});
                        source = testCase.writeInitialModel(wvt);
                        resultPath = fullfile(testCase.folder,"qg-tendencies.json");
                        for provider = testCase.providers
                            [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(resultPath)+" "+provider+" tendencies");
                            testCase.assertEqual(status,0,provider+": "+output);
                            actual = jsondecode(fileread(resultPath));
                            testCase.verifyEqual(actual.diagnosticWorkspaceLiveBytes,0);
                            testCase.verifyEqual(actual.diagnosticForcingEvaluationCount,numel(wvt.forcing));
                            maximumError = 0;
                            for instance = reshape(actual.tendencies,1,[])
                                name = "Fqgpv_"+replace(string(instance.name),[" ","-"],"_");
                                index = find(names==name);
                                testCase.assertNumElements(index,1,name);
                                reference = expected{index};
                                values = instance.fields.Fqgpv;
                                error = max(abs(values(:)-reference(:)))/max(max(abs(reference(:))),realmin);
                                testCase.verifyLessThanOrEqual(error,1e-12,"j="+j+" aa="+antialias+" "+provider+" "+name);
                                maximumError = max(maximumError,error);
                            end
                            testCase.verifyEqual(numel(actual.tendencies),operation.nVarOut);
                            fprintf('QG_FORCING_DIAGNOSTICS j=%d aa=%d grid=%s provider=%s comparisons=%d max_relative=%.17g\n',j,antialias,mat2str(grid{1}),provider,operation.nVarOut,maximumError);
                        end
                        testCase.verifyEqual(wvt.A0,before);
                    end
                end
            end
        end
        function fullGridStratifiedTendenciesMatchMatlab(testCase)
            for antialias = [false true]
                for grid = {[8 6 9],[9 7 10]}
                    wvt = WVTransformStratifiedQG([17000 11000 1000],grid{1},Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
                    n = reshape(1:numel(wvt.A0),size(wvt.A0));
                    wvt.A0 = 1e-6*complex(sin(.17*n),cos(.23*n)).*(wvt.Kh>0);
                    wvt.t = 37;
                    fixed = WVFixedAmplitudeForcing(wvt,name="held-pv",A0_indices=uint64((1:numel(wvt.A0))'),A0bar=wvt.A0(:));
                    forces = [WVNonlinearAdvection(wvt),WVBottomFrictionLinear(wvt,r=2.5e-7),WVBottomFrictionQuadratic(wvt,Cd=.002),WVVerticalDiffusivity(wvt,kappa_z=.002),WVBetaPlanePVAdvection(wvt),WVAdaptiveDamping(wvt),testCase.forcing(wvt,"WVNarrowBandGeostrophicForcing"),fixed];
                    if ~antialias, forces = [forces,WVAntialiasing(wvt,Nj=3)]; end %#ok<AGROW>
                    wvt.setForcing(forces);
                    before = wvt.A0;
                    operation = SpatialForcingOperation(wvt);
                    expected = cell(1,operation.nVarOut);
                    [expected{:}] = operation.compute(wvt);
                    names = string({operation.outputVariables.name});
                    source = testCase.writeInitialModel(wvt);
                    resultPath = fullfile(testCase.folder,"stratified-tendencies.json");
                    for provider = testCase.providers
                        [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(resultPath)+" "+provider+" tendencies");
                        testCase.assertEqual(status,0,provider+": "+output);
                        actual = jsondecode(fileread(resultPath));
                        testCase.verifyEqual(actual.diagnosticWorkspaceLiveBytes,0);
                        testCase.verifyEqual(actual.diagnosticForcingEvaluationCount,numel(wvt.forcing));
                        maximumError = 0;
                        for instance = reshape(actual.tendencies,1,[])
                            name = "Fqgpv_"+replace(string(instance.name),[" ","-"],"_");
                            index = find(names==name);
                            testCase.assertNumElements(index,1,name);
                            reference = expected{index};
                            values = instance.fields.Fqgpv;
                            error = max(abs(values(:)-reference(:)))/max(max(abs(reference(:))),realmin);
                            testCase.verifyLessThanOrEqual(error,1e-12,"aa="+antialias+" "+provider+" "+name);
                            maximumError = max(maximumError,error);
                        end
                        testCase.verifyEqual(numel(actual.tendencies),operation.nVarOut);
                        fprintf('STRATIFIED_FORCING_DIAGNOSTICS aa=%d grid=%s provider=%s comparisons=%d max_relative=%.17g\n',antialias,mat2str(grid{1}),provider,operation.nVarOut,maximumError);
                    end
                    testCase.verifyEqual(wvt.A0,before);
                end
            end
        end
        function explicitAntialiasingMatchesMatlab(testCase)
            for family = ["hydrostatic","nonhydrostatic","barotropic"]
                for grid = {[8 6 5],[9 7 7]}
                    for retained = [0 1 3]
                        wvt = testCase.transform(family,grid{1},false);
                        nonlinear = WVNonlinearAdvection(wvt);
                        filter = WVAntialiasing(wvt,Nj=retained);
                        wvt.setForcing([nonlinear,filter]);
                        testCase.compareTendency(wvt);
                    end
                end
            end
        end
        function explicitFilterPreservesOrderingAndFixedAmplitudes(testCase)
            for family = ["hydrostatic","nonhydrostatic","barotropic"]
                for reverseOrder = [false true]
                    wvt = testCase.transform(family,[9 6 5],false);
                    filter = WVAntialiasing(wvt,Nj=2);
                    adaptive = WVAdaptiveDamping(wvt);
                    fixed = WVFixedAmplitudeForcing(wvt,name="fixed",A0_indices=uint64(2),A0bar=wvt.A0(2));
                    forces = [adaptive,WVNonlinearAdvection(wvt),filter,fixed];
                    if reverseOrder, forces = fliplr(forces); end
                    wvt.setForcing(forces);
                    before = wvt.A0;
                    testCase.compareTendency(wvt);
                    testCase.verifyEqual(wvt.A0,before);
                end
            end
        end
    end
    methods (Test,TestTags="full")
        function adaptiveClosuresUseMatlabEffectiveResolution(testCase)
            for family = ["hydrostatic","nonhydrostatic","barotropic"]
                wvt = testCase.transform(family,[8 6 5],false);
                forces = [WVNonlinearAdvection(wvt),WVAntialiasing(wvt,Nj=3),WVAdaptiveDamping(wvt),testCase.forcing(wvt,"WVFixedAmplitudeForcing")];
                if family ~= "barotropic"
                    forces = [forces,WVHorizontalDamping(wvt,nu=.125,kappa=.002),WVVerticalDamping(wvt,nu=.125,kappa=.002),WVVerticalDiffusivity(wvt,kappa_z=.002)]; %#ok<AGROW>
                end
                wvt.setForcing(forces);
                source = testCase.writeInitialModel(wvt);
                originalTime = wvt.t;
                control = WVModel(wvt);
                expectedStep = control.timeStepForCFL(.5);
                expectedResolution = wvt.effectiveHorizontalGridResolution;
                control.setupIntegrator();
                control.integrateToTime(originalTime+1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                for provider = testCase.providers
                    runtimePath = fullfile(testCase.folder,"adaptive.nc");
                    copyfile(source,runtimePath);
                    requestPath = fullfile(testCase.folder,"request.json");
                    reportPath = fullfile(testCase.folder,"report.json");
                    WVModel.writePortableRunRequest(requestPath,runtimePath,finalTime=originalTime+1,fftProvider=replace(provider,"native","native-fftw"),reportPath=reportPath);
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(requestPath));
                    testCase.assertEqual(status,0,output);
                    report = jsondecode(fileread(reportPath));
                    testCase.verifyEqual(report.integrationRequest.candidates.effectiveHorizontalGridResolution,expectedResolution,RelTol=1e-12);
                    testCase.verifyEqual(report.integrationRequest.selectedStep,expectedStep,RelTol=1e-12);
                    actual = WVModel.modelFromFile(char(runtimePath));
                    cleanup = onCleanup(@()actual.closeNetCDFFile());
                    families = "A0";
                    if family ~= "barotropic", families = ["Ap","Am","A0"]; end
                    for name = families
                        expected = control.wvt.(name);
                        testCase.verifyLessThanOrEqual(max(abs(actual.wvt.(name)-expected),[],"all"),2e-10*max(abs(expected),[],"all")+1e-16,provider+" adaptive "+name);
                    end
                    clear cleanup
                end
            end
        end
        function malformedClosuresLeaveAppendTargetsUnchanged(testCase)
            for identity = ["WVAntialiasing","WVHorizontalDamping","WVVerticalDamping","WVVerticalDiffusivity"]
                wvt = testCase.transform("nonhydrostatic",[8 6 5],false);
                wvt.setForcing(testCase.forcing(wvt,identity));
                source = testCase.writeInitialModel(wvt);
                group = forcingGroup(ncinfo(source),identity);
                testCase.assertNotEmpty(group);
                parameter = "nu";
                if identity == "WVAntialiasing", parameter = "Nj";
                elseif identity == "WVVerticalDiffusivity", parameter = "kappa_z"; end
                ncwrite(source,group+"/"+parameter,-1);
                before = readBytes(source);
                [status,~] = cleanSystem(shellQuote(testCase.runner)+" "+shellQuote(source)+" --restart-mode model --output-policy append --integrator fixed-rk4 --delta-t 0.5 --final-time 38 --fft-provider reference");
                testCase.verifyNotEqual(status,0);
                testCase.verifyEqual(readBytes(source),before,identity+" mutated an invalid append target.");
            end
        end
        function catalogPairsMatchMatlabAndAppend(testCase)
            matrix = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
            for row = reshape(matrix.rows,1,[])
                if string(row.matlab.applicability) ~= "applicable" || startsWith(string(row.configuration),"stratified-qg"), continue; end
                config = matrix.configurations(string({matrix.configurations.id}) == string(row.configuration));
                family = "hydrostatic";
                if string(config.transform) == "WVTransformBarotropicQG", family = "barotropic";
                elseif ~config.isHydrostatic, family = "nonhydrostatic"; end
                wvt = testCase.transform(family,[8 6 7],config.shouldAntialias);
                wvt.setForcing(testCase.forcing(wvt,string(row.forcing)));
                testCase.compareTendency(wvt);
                testCase.compareContinuation(wvt,false);
            end
        end
        function orderedClosureGraphsRestartAndConvertResolution(testCase)
            for family = ["hydrostatic","nonhydrostatic","barotropic"]
                wvt = testCase.transform(family,[8 6 5],false);
                forces = [WVAdaptiveDamping(wvt),WVNonlinearAdvection(wvt),WVAntialiasing(wvt,Nj=3),testCase.forcing(wvt,"WVFixedAmplitudeForcing")];
                if family ~= "barotropic"
                    forces = [forces,WVHorizontalDamping(wvt,nu=.125,kappa=.002),WVVerticalDamping(wvt,nu=.125,kappa=.002),WVVerticalDiffusivity(wvt,kappa_z=.002)]; %#ok<AGROW>
                end
                wvt.setForcing(forces);
                testCase.compareContinuation(wvt,true);
                targetGrid = [9 8 7];
                if family == "barotropic", targetGrid = targetGrid(1:2); end
                converted = wvt.waveVortexTransformWithResolution(targetGrid);
                testCase.compareTendency(converted);
                testCase.compareContinuation(converted,false);
            end
        end
        function laplacianDampingMatchesMatlab(testCase)
            for family = ["hydrostatic","nonhydrostatic"]
                for grid = {[8 6 5],[9 7 7]}
                    for antialias = [false true]
                        for rates = {[0 0],[.125 2.5e-5]}
                            for direction = ["horizontal","vertical"]
                                wvt = testCase.transform(family,grid{1},antialias);
                                values = rates{1};
                                if direction == "horizontal"
                                    force = WVHorizontalDamping(wvt,nu=values(1),kappa=values(2));
                                else
                                    force = WVVerticalDamping(wvt,nu=values(1),kappa=values(2));
                                end
                                wvt.setForcing(force);
                                testCase.compareTendency(wvt);
                            end
                        end
                    end
                end
            end
        end
        function verticalDiffusivityMatchesMatlab(testCase)
            for family = ["hydrostatic","nonhydrostatic"]
                for grid = {[8 6 5],[9 7 7]}
                    for antialias = [false true]
                        for flag = [false true]
                            wvt = testCase.transform(family,grid{1},antialias);
                            wvt.setForcing(WVVerticalDiffusivity(wvt,kappa_z=.0025,shouldForceMeanDensityAnomaly=flag));
                            testCase.compareTendency(wvt);
                        end
                    end
                end
            end
        end
        function constantNarrowBandAndOrderedClosuresMatchMatlab(testCase)
            for family = ["hydrostatic","nonhydrostatic"]
                for antialias = [false true]
                    wvt = testCase.transform(family,[9 6 5],antialias);
                    narrow = WVNarrowBandGeostrophicForcing(wvt,initialPV="none",k_f=2*wvt.dk,j_f=1);
                    wvt.setForcing(narrow);
                    testCase.compareTendency(wvt);
                    forces = [WVVerticalDiffusivity(wvt),WVHorizontalDamping(wvt),WVNonlinearAdvection(wvt),WVVerticalDamping(wvt),WVAdaptiveDamping(wvt),narrow];
                    if ~antialias, forces(end+1) = WVAntialiasing(wvt,Nj=2); end %#ok<AGROW>
                    wvt.setForcing(forces);
                    testCase.compareTendency(wvt);
                end
            end
        end
    end

    methods (Access=private)
        function report = runForcingOutputContinuation(testCase,outputPath,provider,scenario,finalTime)
            requestPath = fullfile(testCase.folder,"forcing-output-request.json");
            method = "fixed-rk4";
            if scenario == 3, method = "adaptive-rk45"; end
            if scenario == 4, method = "adaptive-rk78"; end
            if scenario <= 2
                WVModel.writePortableRunRequest(requestPath,outputPath,method=method,finalTime=finalTime,initialStep=.25,fftProvider=replace(provider,"native","native-fftw"),reportPath="forcing-output-report.json");
            else
                WVModel.writePortableRunRequest(requestPath,outputPath,method=method,finalTime=finalTime,initialStep=.25,maximumStep=.25,fftProvider=replace(provider,"native","native-fftw"),reportPath="forcing-output-report.json");
            end
            command = shellQuote(testCase.runner)+" --request "+shellQuote(requestPath);
            [status,output] = cleanSystem(command);
            testCase.assertEqual(status,0,provider+" scenario="+scenario+": "+output);
            report = jsondecode(fileread(fullfile(testCase.folder,"forcing-output-report.json")));
        end
        function force = forcing(~,wvt,identity)
            switch identity
                case "WVFixedAmplitudeForcing"
                    force = WVFixedAmplitudeForcing(wvt,name="catalog fixed",A0_indices=uint64(2),A0bar=wvt.A0(2));
                case "WVNarrowBandGeostrophicForcing"
                    force = WVNarrowBandGeostrophicForcing(wvt,initialPV="none",k_f=2*wvt.dk,j_f=1);
                case "WVPseudoTopographicWaveGeneration"
                    force = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=2*cos(2*pi*reshape(wvt.x,[],1)/wvt.Lx)*ones(1,wvt.Ny),barotropicVelocityAmplitude=[.01;.005],frequency=1e-4,rampDuration=0);
                case "WVAntialiasing"
                    force = WVAntialiasing(wvt,Nj=2);
                case {"WVHorizontalDamping","WVVerticalDamping"}
                    force = feval(identity,wvt,nu=.125,kappa=.002);
                case "WVVerticalDiffusivity"
                    force = WVVerticalDiffusivity(wvt,kappa_z=.002);
                otherwise
                    force = feval(identity,wvt);
            end
        end
        function source = writeInitialModel(testCase,wvt)
            model = WVModel(wvt);
            source = fullfile(testCase.folder,"source.nc");
            outputFile = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
            outputFile.outputTimesForIntegrationPeriod(wvt.t,wvt.t+2);
            outputFile.writeTimeStepToOutputFile(wvt.t);
            model.closeNetCDFFile();
        end
        function compareContinuation(testCase,wvt,segmented)
            source = testCase.writeInitialModel(wvt);
            originalTime = wvt.t;
            controlPath = fullfile(testCase.folder,"control.nc");
            copyfile(source,controlPath);
            control = WVModel.modelFromFile(char(controlPath));
            cleanup = onCleanup(@()control.closeNetCDFFile());
            control.setupIntegrator(integratorType="fixed",deltaT=.5);
            duration = 1+double(segmented);
            control.integrateToTime(originalTime+duration,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            control.closeNetCDFFile();
            clear cleanup
            for provider = testCase.providers
                runtimeProvider = replace(provider,"native","native-fftw");
                runtimePath = fullfile(testCase.folder,"runtime.nc");
                copyfile(source,runtimePath);
                for stop = 1:duration
                    command = shellQuote(testCase.runner)+" "+shellQuote(runtimePath)+" --restart-mode model --output-policy append --integrator fixed-rk4 --delta-t 0.5 --final-time "+string(originalTime+stop)+" --fft-provider "+runtimeProvider;
                    [status,output] = cleanSystem(command);
                    testCase.assertEqual(status,0,output);
                end
                actual = WVModel.modelFromFile(char(runtimePath));
                cleanup = onCleanup(@()actual.closeNetCDFFile());
                families = "A0";
                if ~isa(wvt,"WVTransformBarotropicQG"), families = ["Ap","Am","A0"]; end
                for name = families
                    expected = control.wvt.(name);
                    testCase.verifyLessThanOrEqual(max(abs(actual.wvt.(name)-expected),[],"all"),2e-10*max(abs(expected),[],"all")+1e-16,provider+" "+name);
                end
                testCase.verifyEqual(string({actual.wvt.forcing.name}),string({wvt.forcing.name}));
                for i = 1:numel(wvt.forcing)
                    testCase.verifyEqual(string(class(actual.wvt.forcing(i))),string(class(wvt.forcing(i))));
                end
                if segmented
                    uninterruptedPath = fullfile(testCase.folder,"uninterrupted.nc");
                    copyfile(source,uninterruptedPath);
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" "+shellQuote(uninterruptedPath)+" --restart-mode model --output-policy append --integrator fixed-rk4 --delta-t 0.5 --final-time "+string(originalTime+duration)+" --fft-provider "+runtimeProvider);
                    testCase.assertEqual(status,0,output);
                    uninterrupted = WVModel.modelFromFile(char(uninterruptedPath));
                    otherCleanup = onCleanup(@()uninterrupted.closeNetCDFFile());
                    for name = families
                        testCase.verifyEqual(actual.wvt.(name),uninterrupted.wvt.(name));
                    end
                    clear otherCleanup
                end
                clear cleanup
            end
        end
        function wvt = transform(~,family,grid,antialias)
            if family == "barotropic"
                wvt = WVTransformBarotropicQG([17000 11000],grid(1:2),j=1,shouldAntialias=antialias);
            else
                wvt = WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=family == "hydrostatic",shouldAntialias=antialias);
            end
            originalStream = rng;
            cleanup = onCleanup(@()rng(originalStream));
            rng(9321);
            wvt.initWithRandomFlow(uvMax=.01);
            wvt.t = 37;
            wvt.removeAllForcing();
        end
        function compareTendency(testCase,wvt)
            if isa(wvt,"WVTransformBarotropicQG")
                expected = struct(F0=wvt.nonlinearFlux());
            else
                [Fp,Fm,F0] = wvt.nonlinearFlux();
                expected = struct(Fp=Fp,Fm=Fm,F0=F0);
            end
            source = testCase.writeInitialModel(wvt);
            outputPath = fullfile(testCase.folder,"rhs.json");
            for provider = testCase.providers
                [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(outputPath)+" "+provider);
                testCase.assertEqual(status,0,output);
                actual = jsondecode(fileread(outputPath));
                % The direct reference FFT intentionally uses per-line scratch vectors.
                if provider == "native"
                    testCase.verifyEqual(actual.applicationAllocations,0,"Prepared native forcing execution allocated application memory.");
                end
                for name = string(fieldnames(expected))'
                    reference = expected.(name)(:);
                    result = complex(actual.(name).real,actual.(name).imag);
                    tolerance = 2e-11*max(abs(reference))+1e-18;
                    err = max(abs(result(:)-reference));
                    testCase.verifyLessThanOrEqual(err,tolerance,provider+" "+name);
                end
            end
        end
    end
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
function [status,output] = cleanSystem(command)
[status,output] = system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end

function result = forcingGroup(info,identity)
result = "";
for attribute = reshape(info.Attributes,1,[])
    if string(attribute.Name) == "AnnotatedClass" && string(attribute.Value) == identity
        result = string(info.Name);
        return
    end
end
for group = reshape(info.Groups,1,[])
    result = forcingGroup(group,identity);
    if result ~= "", return; end
end
end
function bytes = readBytes(path)
file = fopen(path,"r");
cleanup = onCleanup(@()fclose(file));
bytes = fread(file,Inf,"*uint8");
end

function wvt = diagnosticWaveTransform(family,grid,antialias)
switch family
    case "hydrostatic"
        wvt = WVTransformHydrostatic([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
    case "boussinesq"
        wvt = WVTransformBoussinesq([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
    otherwise
        wvt = WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=antialias);
end
n = reshape(1:numel(wvt.A0),size(wvt.A0));
wave = wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0;
inertial = (wvt.K.^2+wvt.L.^2)==0;
wvt.Ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
wvt.Am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
wvt.Am(inertial) = conj(wvt.Ap(inertial));
wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n));
wvt.A0(inertial) = .003*sin(n(inertial)).*(wvt.J(inertial)>0);
wvt.t0 = 17;
wvt.t = 37;
wvt.removeAllForcing();
end
