classdef TestPortableDiagnostics < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        executable (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function setup(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_DIAGNOSTIC_DUMP"));
            testCase.providers = "reference";
            if getenv("WV_DIAGNOSTIC_NATIVE") == "1", testCase.providers = ["reference","native"]; end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 4 --target WVDiagnosticFieldDump wave-vortex-run");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVDiagnosticFieldDump");
            end
        end
    end
    methods (Test,TestTags="full")
        function allTransformDiagnosticsMatchMatlab(testCase)
            catalog = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            families = ["constant-hydrostatic","constant-nonhydrostatic","barotropic","stratified-qg","hydrostatic","boussinesq"];
            evidenceRows = struct(configuration={},grid={},provider={},maximumRelativeError={},variableCount={}, ...
                retainedBytes={},scratchHighWaterBytes={},scratchLiveBytes={},primitiveOutputCount={});
            for family = families
                for antialias = [false true]
                    for odd = [false true]
                        grid = [8+double(odd),6+double(odd),9+double(odd)];
                        wvt = testCase.transform(family,grid,antialias);
                        configuration = family+"-aa"+double(antialias);
                        rows = catalog.contracts(string({catalog.contracts.configuration})==configuration);
                        names = reshape(arrayfun(@(row)string(row.metadata.name),rows),1,[]);
                        % Nonlinear density fits have dedicated analytic and bounded-parity coverage.
                        densityFields = intersect(names,["rho_nm","eta_true","ape","apv"],"stable");
                        names = names(~endsWith(names,"_portable_catalog_forcing") & ~ismember(names,["Ap","Am","A0",densityFields]));
                        for row = reshape(rows,1,[])
                            if string(row.authority)=="known-variable-factory"
                                if string(row.component)==""
                                    operation = wvt.operationForKnownVariable('energy');
                                else
                                    operation = wvt.operationForKnownVariable('energy',flowComponent=wvt.flowComponentWithName(string(row.component)));
                                end
                                wvt.addOperation(operation);
                            end
                        end
                        expected = struct;
                        for name = reshape(names,1,[])
                            value = wvt.(name);
                            expected.(name) = value(:);
                        end
                        model = WVModel(wvt,shouldUseLinearDynamics=family~="barotropic");
                        checkpoint = fullfile(testCase.folder,"source.nc");
                        file = model.createNetCDFFileForModelOutput(checkpoint,outputInterval=.5,shouldOverwriteExisting=true);
                        file.outputTimesForIntegrationPeriod(wvt.t,wvt.t+1);
                        file.writeTimeStepToOutputFile(wvt.t);
                        model.closeNetCDFFile();
                        request = fullfile(testCase.folder,"request.json");
                        writeText(request,jsonencode(struct(fields=names)));
                        for provider = testCase.providers
                            resultPath = fullfile(testCase.folder,"result.json");
                            [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(checkpoint)+" "+shellQuote(request)+" "+shellQuote(resultPath)+" "+provider);
                            testCase.assertEqual(status,0,configuration+" odd="+odd+" "+provider+": "+output);
                            actual = jsondecode(fileread(resultPath));
                            maximumError = 0;
                            for name = reshape(names,1,[])
                                value = actual.fields.(name).real;
                                if isfield(actual.fields.(name),'imag'), value = complex(value,actual.fields.(name).imag); end
                                scale = max(abs(expected.(name)));
                                error = max(abs(value(:)-expected.(name)))/max(scale,realmin);
                                maximumError = max(maximumError,error);
                                testCase.verifyLessThanOrEqual(error,1e-12,configuration+" odd="+odd+" "+provider+" "+name);
                            end
                            evidenceRows(end+1) = struct(configuration=configuration,grid=grid,provider=provider, ...
                                maximumRelativeError=maximumError,variableCount=numel(names), ...
                                retainedBytes=actual.retainedBytes,scratchHighWaterBytes=actual.scratchHighWaterBytes, ...
                                scratchLiveBytes=actual.scratchLiveBytes,primitiveOutputCount=actual.primitiveOutputs); %#ok<AGROW>
                            testCase.verifyEqual(actual.scratchLiveBytes,0);
                            testCase.verifyEqual(actual.diagnosticEvaluations,1);
                            testCase.verifyLessThan(actual.primitiveOutputs,numel(names));
                        end
                    end
                end
            end
            reportPath = string(getenv("WV_DIAGNOSTIC_REPORT"));
            if reportPath~=""
                [status,commit] = cleanSystem("git -C "+shellQuote(testCase.root)+" rev-parse HEAD");
                testCase.assertEqual(status,0,commit);
                report = struct(schema="portable-diagnostic-numerical-v1",issue=305,sourceCommit=string(strtrim(commit)), ...
                    matlabRelease=string(version('-release')),rows=evidenceRows, ...
                    passes=all([evidenceRows.maximumRelativeError]<=1e-12));
                writeText(reportPath,jsonencode(report,PrettyPrint=true));
            end
        end
        function outputContinuationMatchesMatlab(testCase)
            catalog = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            runner = fullfile(fileparts(testCase.executable),"wave-vortex-run");
            for family = ["constant-hydrostatic","constant-nonhydrostatic","barotropic","stratified-qg","hydrostatic","boussinesq"]
                wvt = testCase.transform(family,[8 6 9],true);
                rows = catalog.contracts(string({catalog.contracts.configuration})==family+"-aa1");
                % Forcing templates require actual instance bindings, covered
                % by TestPortableStableForcing's output continuation matrix.
                rows = rows([rows.ordinal]>=23 & string({rows.runtimeStatus})=="implemented" & string({rows.authority})~="forcing-instance-template");
                % Density recovery has its own prospective error bounds and
                % complete output/continuation coverage in TestPortableDensityOutput.
                densityNames = ["rho_nm","eta_true","ape","apv"];
                rows = rows(~ismember(arrayfun(@(row)string(row.metadata.name),rows),densityNames));
                names = reshape(arrayfun(@(row)string(row.metadata.name),rows),1,[]);
                for row = reshape(rows,1,[])
                    if string(row.authority)=="known-variable-factory"
                        wvt.addOperation(wvt.operationForKnownVariable('energy',flowComponent=wvt.flowComponentWithName(string(row.component))));
                    end
                end
                model = WVModel(wvt,shouldUseLinearDynamics=family~="barotropic");
                model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
                source = fullfile(testCase.folder,"lifecycle-source.nc");
                file = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
                dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
                dense.addObservingSystem(WVEulerianFields(model,fieldNames=cellstr(names)));
                secondSource = fullfile(testCase.folder,"lifecycle-second-source.nc");
                model.createNetCDFFileForModelOutput(secondSource,outputInterval=.5,shouldOverwriteExisting=true);
                for outputFile = reshape(model.outputFiles,1,[])
                    outputFile.outputTimesForIntegrationPeriod(37,38);
                    outputFile.writeTimeStepToOutputFile(37);
                end
                model.closeNetCDFFile();
                ncwriteatt(source,"/","portableFileIdentifier","diagnostic-primary");
                ncwriteatt(secondSource,"/","portableFileIdentifier","diagnostic-secondary");
                controlPath = fullfile(testCase.folder,"lifecycle-matlab.nc"); copyfile(source,controlPath);
                control = WVModel.modelFromFile(char(controlPath));
                control.setupIntegrator(integratorType="fixed",deltaT=.25);
                control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                control.closeNetCDFFile();
                for provider = testCase.providers
                    for scenario = 1:3
                        segmented = scenario==2;
                        outputPath = fullfile(testCase.folder,"lifecycle-runtime.nc"); copyfile(source,outputPath);
                        secondOutput = fullfile(testCase.folder,"lifecycle-second-runtime.nc"); copyfile(secondSource,secondOutput);
                        times = 38;
                        if segmented, times = [37.5 38]; end
                        for finalTime = times
                            requestPath = fullfile(testCase.folder,"lifecycle-request.json");
                            if scenario==3
                                WVModel.writePortableRunRequest(requestPath,[outputPath,secondOutput],method="adaptive-rk45",finalTime=finalTime,initialStep=.25,maximumStep=.25,fftProvider=replace(provider,"native","native-fftw"),reportPath="lifecycle-report.json");
                            else
                                WVModel.writePortableRunRequest(requestPath,[outputPath,secondOutput],method="fixed-rk4",finalTime=finalTime,initialStep=.25,fftProvider=replace(provider,"native","native-fftw"),reportPath="lifecycle-report.json");
                            end
                            [status,output] = cleanSystem(shellQuote(runner)+" --request "+shellQuote(requestPath));
                            testCase.assertEqual(status,0,family+" "+provider+" segmented="+segmented+": "+output);
                            report = jsondecode(fileread(fullfile(testCase.folder,"lifecycle-report.json")));
                            testCase.verifyEqual(report.state.stepCount,4-2*double(segmented));
                            testCase.verifyEqual(report.diagnosticEvaluation.workspaceLiveBytes,0);
                            testCase.verifyGreaterThan(report.diagnosticEvaluation.evaluationCount,0);
                            if family~="stratified-qg"
                                testCase.verifyGreaterThan(report.diagnosticEvaluation.workspaceHighWaterBytes,0);
                            end
                            % Prepared scientific buffers are already counted in retained storage.
                            testCase.verifyGreaterThanOrEqual(report.livenessBytes.fullModelMaximumLive, ...
                                report.livenessBytes.fullModelRetained+report.diagnosticEvaluation.additionalTransientHighWaterBytes);
                        end
                        for groupName = ["wave-vortex","dense"]
                            info = ncinfo(controlPath,"/"+groupName);
                            for variable = reshape(info.Variables,1,[])
                                name = string(variable.Name);
                                if ~ismember(erase(name,["_real","_imag"]),[names,"t","Ap","Am","A0"]), continue; end
                                expected = ncread(controlPath,"/"+groupName+"/"+name);
                                actual = ncread(outputPath,"/"+groupName+"/"+name);
                                testCase.verifySize(actual,size(expected),family+" "+name);
                                error = max(abs(actual(:)-expected(:)))/max(max(abs(expected(:))),realmin);
                                testCase.verifyLessThanOrEqual(error,1e-12,family+" "+provider+" scenario="+scenario+" "+groupName+" "+name);
                                if groupName=="wave-vortex"
                                    second = ncread(secondOutput,"/wave-vortex/"+name);
                                    error = max(abs(second(:)-expected(:)))/max(max(abs(expected(:))),realmin);
                                    testCase.verifyLessThanOrEqual(error,1e-12,family+" "+provider+" scenario="+scenario+" second file "+name);
                                end
                            end
                        end
                    end
                end
            end
        end

        function phaseOutputContinuationMatchesMatlab(testCase)
            runner = fullfile(fileparts(testCase.executable),"wave-vortex-run");
            testCase.assertTrue(isfile(runner));
            names = ["phase","conjPhase"];
            rows = struct(configuration={},grid={},provider={},maximumAbsoluteError={},scenarios={},ordinaryRecords={},denseRecords={});
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for antialias = [false true]
                    grid = [8 6 9]+double(~antialias);
                    wvt = testCase.transform(family,grid,antialias);
                    model = WVModel(wvt,shouldUseLinearDynamics=true);
                    model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
                    source = fullfile(testCase.folder,"phase-source.nc");
                    file = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
                    dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
                    dense.addObservingSystem(WVEulerianFields(model,fieldNames=cellstr(names)));
                    file.outputTimesForIntegrationPeriod(37,38);
                    file.writeTimeStepToOutputFile(37);
                    model.closeNetCDFFile();
                    control = fullfile(testCase.folder,"phase-matlab.nc"); copyfile(source,control);
                    testCase.continuePhaseInMatlab(control,38);
                    prefix = fullfile(testCase.folder,"phase-prefix.nc"); copyfile(source,prefix);
                    testCase.continuePhaseInMatlab(prefix,37.5);
                    for provider = testCase.providers
                        dumpRequest = fullfile(testCase.folder,"phase-fields.json");
                        dumpPath = fullfile(testCase.folder,"phase-fields-result.json");
                        writeText(dumpRequest,jsonencode(struct(fields=names)));
                        [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(dumpRequest)+" "+shellQuote(dumpPath)+" "+provider);
                        testCase.assertEqual(status,0,output);
                        dump = jsondecode(fileread(dumpPath));
                        maximumError = 0;
                        for name = names
                            actual = complex(dump.fields.(name).real,dump.fields.(name).imag);
                            expected = wvt.(name);
                            error = max(abs(actual(:)-expected(:)));
                            testCase.verifyLessThanOrEqual(error,1e-12,family+" "+provider+" "+name);
                            maximumError = max(maximumError,error);
                        end
                        scenarios = ["whole","cppRestart","matlabToCpp","cppToMatlab"];
                        for scenario = scenarios
                            destination = fullfile(testCase.folder,"phase-runtime.nc");
                            if scenario=="matlabToCpp", copyfile(prefix,destination);
                            else, copyfile(source,destination); end
                            if ismember(scenario,["cppRestart","cppToMatlab"])
                                testCase.continuePhaseInCpp(destination,provider,37.5,runner);
                                testCase.verifyPhaseRecords(destination,wvt,37.5);
                            end
                            if scenario=="cppToMatlab"
                                testCase.continuePhaseInMatlab(destination,38);
                            else
                                testCase.continuePhaseInCpp(destination,provider,38,runner);
                            end
                            maximumError = max(maximumError,testCase.verifyPhaseRecords(destination,wvt,38));
                            for group = ["wave-vortex","dense"]
                                expectedTimes = ncread(control,"/"+group+"/t");
                                testCase.verifyEqual(ncread(destination,"/"+group+"/t"),expectedTimes);
                                for name = names
                                    for part = ["real","imag"]
                                        variable = "/"+group+"/"+name+"_"+part;
                                        expected = ncread(control,variable);
                                        actual = ncread(destination,variable);
                                        error = max(abs(actual(:)-expected(:)));
                                        testCase.verifyLessThanOrEqual(error,1e-12,family+" "+provider+" "+scenario+" "+variable);
                                        maximumError = max(maximumError,error);
                                    end
                                end
                            end
                        end
                        rows(end+1) = struct(configuration=family+"-aa"+double(antialias),grid=grid,provider=provider, ...
                            maximumAbsoluteError=maximumError,scenarios=scenarios,ordinaryRecords=3,denseRecords=9); %#ok<AGROW>
                        fprintf("PHASE_OUTPUT_PASS %s aa=%d provider=%s maximumAbsoluteError=%.17g scenarios=4\n",family,antialias,provider,maximumError);
                    end
                end
            end
            reportPath = string(getenv("WV_PHASE_REPORT"));
            if reportPath~=""
                report = struct(schema="portable-phase-output-numerical-v1",rows=rows,matlabRelease=string(version('-release')), ...
                    initialTime=37,referenceTime=17,finalTime=38,step=.25,denseInterval=.125, ...
                    passes=all([rows.maximumAbsoluteError]<=1e-12));
                writeText(reportPath,jsonencode(report,PrettyPrint=true));
            end
        end

        function existingComponentOperationIsPreserved(testCase)
            wvt = testCase.transform("constant-hydrostatic",[8 6 9],false);
            operation = wvt.operationForKnownVariable('energy',flowComponent=wvt.flowComponentWithName("wave"));
            wvt.addOperation(operation);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            observer = WVEulerianFields(model,fieldNames={'energy_w'});
            testCase.verifyEqual(wvt.propertyAnnotationWithName('energy_w').modelOp,operation);
            testCase.verifyEqual(observer.fieldNames,"energy_w");
        end

        function incompatibilitiesPreserveMatlabContracts(testCase)
            wvt = testCase.transform("hydrostatic",[8 6 9],false);
            for name = ["phase","conjPhase"]
                testCase.verifyTrue(wvt.propertyAnnotationWithName(name).isComplex);
                testCase.verifyFalse(isreal(wvt.(name)));
            end
            testCase.verifyFalse(ismember('shouldUseTrueNoMotionProfile',wvt.annotatedPropertyNames));
            solver = WVNoMotionProfileOperation();
            testCase.verifyEqual(solver.solver,"dampedLeastSquares");
            for legacySolver = ["lsqnonlin","fminsearch"]
                legacy = WVNoMotionProfileOperation(solver=legacySolver);
                testCase.verifyEqual(legacy.solver,legacySolver);
            end
            testCase.verifyFalse(ismember('solver',wvt.annotatedPropertyNames));
            testCase.verifyTrue(wvt.shouldUseTrueNoMotionProfile);
            wvt.shouldUseTrueNoMotionProfile = false;
            filePath = fullfile(testCase.folder,"profile-contract.nc");
            file = wvt.writeToFile(filePath,shouldOverwriteExisting=true); file.close();
            [restored,file] = WVTransform.waveVortexTransformFromFile(filePath);
            cleanup = onCleanup(@()file.close());
            testCase.verifyTrue(restored.shouldUseTrueNoMotionProfile);
        end
    end
    methods (Access=private)
        function continuePhaseInMatlab(testCase,path,finalTime)
            model = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()model.closeNetCDFFile());
            testCase.assertEqual(model.wvt.t0,17);
            testCase.assertEqual(model.wvt.conjPhase,conj(model.wvt.phase));
            model.setupIntegrator(integratorType="fixed",deltaT=.25);
            model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.assertEqual(model.wvt.t,finalTime);
        end
        function continuePhaseInCpp(testCase,path,provider,finalTime,runner)
            request = fullfile(testCase.folder,"phase-run.json");
            reportPath = fullfile(testCase.folder,"phase-run-report.json");
            WVModel.writePortableRunRequest(request,path,method="fixed-rk4",finalTime=finalTime,initialStep=.25, ...
                fftProvider=replace(provider,"native","native-fftw"),reportPath=reportPath);
            [status,output] = cleanSystem(shellQuote(runner)+" --request "+shellQuote(request));
            testCase.assertEqual(status,0,output);
            report = jsondecode(fileread(reportPath));
            testCase.assertEqual(string(report.status),"complete");
            testCase.assertTrue(report.integrationRequest.noFallback);
            testCase.assertEqual(string(report.provider.id),replace(provider,"native","native-fftw"));
            testCase.assertGreaterThan(report.integrator.denseOutputEvaluationCount,0);
        end
        function maximumError = verifyPhaseRecords(testCase,path,wvt,finalTime)
            maximumError = 0;
            for group = ["wave-vortex","dense"]
                times = reshape(ncread(path,"/"+group+"/t"),1,[]);
                interval = .5;
                if group=="dense", interval=.125; end
                testCase.assertEqual(times,37:interval:finalTime);
                values = cell(1,2);
                for index = 1:2
                    names = ["phase","conjPhase"];
                    name = names(index);
                    realPath = "/"+group+"/"+name+"_real";
                    imagPath = "/"+group+"/"+name+"_imag";
                    realInfo = ncinfo(path,realPath);
                    imagInfo = ncinfo(path,imagPath);
                    testCase.assertEqual(realInfo.Size,[size(wvt.Ap),numel(times)]);
                    testCase.assertEqual(realInfo.Dimensions,imagInfo.Dimensions);
                    actual = complex(ncread(path,realPath),ncread(path,imagPath));
                    actual = reshape(actual,numel(wvt.Ap),[]);
                    expected = exp(wvt.iOmega(:).*(times-wvt.t0));
                    if index==2, expected=conj(expected); end
                    error = max(abs(actual-expected),[],"all");
                    testCase.verifyLessThanOrEqual(error,1e-12,group+" "+name);
                    testCase.verifyLessThanOrEqual(max(abs(abs(actual)-1),[],"all"),1e-12);
                    testCase.verifyGreaterThan(max(abs(actual(:,end)-actual(:,1))),1e-8);
                    maximumError = max(maximumError,error);
                    values{index} = actual;
                end
                testCase.verifyEqual(values{2},conj(values{1}));
            end
            restored = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()restored.closeNetCDFFile());
            testCase.assertEqual(restored.wvt.t,finalTime);
            testCase.assertEqual(restored.wvt.t0,wvt.t0);
            expected = exp(wvt.iOmega*(finalTime-wvt.t0));
            testCase.verifyEqual(restored.wvt.phase,expected,AbsTol=1e-12);
            testCase.verifyEqual(restored.wvt.conjPhase,conj(expected),AbsTol=1e-12);
        end
    end
    methods (Static,Access=private)
        function wvt = transform(family,grid,antialias)
            switch family
                case "barotropic"
                    wvt = WVTransformBarotropicQG([17000 11000],grid(1:2),j=1,shouldAntialias=antialias);
                case "stratified-qg"
                    wvt = WVTransformStratifiedQG([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
                case "boussinesq"
                    wvt = WVTransformBoussinesq([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
                otherwise
                    wvt = WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=antialias);
            end
            n = reshape(1:numel(wvt.A0),size(wvt.A0));
            wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n));
            if ~ismember(family,["barotropic","stratified-qg"])
                wave = wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0;
                inertial = (wvt.K.^2+wvt.L.^2)==0;
                wvt.Ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
                wvt.Am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
                wvt.Am(inertial) = conj(wvt.Ap(inertial));
                mda = inertial & wvt.J>0;
                wvt.A0(mda) = .003*sin(n(mda));
            end
            if family=="barotropic"
                wvt.A0(wvt.Kh==0) = 0;
                selfConjugate = wvt.dftPrimaryIndices2D==wvt.dftConjugateIndices2D;
                wvt.A0(selfConjugate) = real(wvt.A0(selfConjugate));
            elseif family=="stratified-qg"
                wvt.A0(:,wvt.k==0 & wvt.l==0) = 0;
            end
            wvt.t0 = 17;
            wvt.t = 37;
            wvt.removeAllForcing();
        end
    end
end
function writeText(path,text)
file = fopen(path,'w');
cleanup = onCleanup(@()fclose(file));
fprintf(file,'%s\n',text);
end
function [status,output] = cleanSystem(command)
[status,output] = system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
