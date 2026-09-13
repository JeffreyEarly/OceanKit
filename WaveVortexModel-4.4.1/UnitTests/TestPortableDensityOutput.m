classdef TestPortableDensityOutput < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        runner (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function prepareRunner(testCase)
            testCase.root=string(fileparts(fileparts(mfilename("fullpath"))));
            fixture=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture(PreservingOnFailure=true));
            testCase.folder=string(fixture.Folder);
            testCase.runner=string(getenv("WV_DENSITY_OUTPUT_RUNNER"));
            if testCase.runner=="" && getenv("WV_DIAGNOSTIC_DUMP")~=""
                testCase.runner=fullfile(fileparts(string(getenv("WV_DIAGNOSTIC_DUMP"))),"wave-vortex-run");
            end
            if testCase.runner==""
                build=fullfile(testCase.folder,"build");
                [status,output]=cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output]=cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target wave-vortex-run");
                testCase.assertEqual(status,0,output);
                testCase.runner=fullfile(build,"wave-vortex-run");
            end
            testCase.assertTrue(isfile(testCase.runner));
            testCase.providers="reference";
            if getenv("WV_DIAGNOSTIC_NATIVE")=="1", testCase.providers=["reference","native-fftw"]; end
        end
    end
    methods (Test,TestTags="full")
        function ordinaryAndDenseDensityOutputsMatchMatlab(testCase)
            rows={};
            for family=["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for antialias=[false true]
                    for selection=["actual","initial"]
                        [source,wvt]=testCase.authorSource(family,antialias,selection,1);
                        control=testCase.copyFiles(source,"control");
                        testCase.continueInMatlab(control,selection,38);
                        for provider=testCase.providers
                            actual=testCase.copyFiles(source,"runtime");
                            report=testCase.continueInCpp(actual,selection,provider,38);
                            errors=testCase.compareFiles(actual,control,wvt,selection,38);
                            rows{end+1}=struct(configuration=family,antialias=antialias,reference=selection,provider=provider, ...
                                errors=errors,ordinaryRecords=3,denseRecords=9,contract=report.densityDiagnosticContract); %#ok<AGROW>
                            fprintf("DENSITY_OUTPUT_PASS %s aa=%d reference=%s provider=%s\n",family,antialias,selection,provider);
                        end
                    end
                end
            end
            testCase.writeReport("ordinary-dense",rows);
        end
        function nonlinearDensityOutputPreservesIntegration(testCase)
            rows={};
            for family=["constant-nonhydrostatic","hydrostatic","boussinesq"]
                [source,wvt]=testCase.authorSource(family,true,"actual",1,false,true);
                source=testCase.copyFiles(source,"nonlinear-density-source");
                control=testCase.copyFiles(source,"nonlinear-control");
                testCase.continueInMatlab(control,"actual",38);
                [plainSource,~]=testCase.authorSource(family,true,"actual",1,false,false);
                for provider=testCase.providers
                    actual=testCase.copyFiles(source,"nonlinear-density");
                    plain=testCase.copyFiles(plainSource,"nonlinear-plain");
                    report=testCase.continueInCpp(actual,"actual",provider,38);
                    baseline=testCase.continueInCpp(plain,"actual",provider,38);
                    errors=testCase.compareFiles(actual,control,wvt,"actual",38);
                    profile=ncread(actual,"/wave-vortex/rho_nm");
                    profileEvolution=max(abs(profile(:,end)-profile(:,1)));
                    testCase.verifyGreaterThan(profileEvolution,8*eps(wvt.rho0),"Current distribution must evolve between nonlinear output events.");
                    for name=["Ap","Am","A0"]
                        for part=["real","imag"]
                            variable="/wave-vortex/"+name+"_"+part;
                            testCase.verifyEqual(ncread(actual,variable),ncread(plain,variable),"Density output must not modify accepted coefficients.");
                        end
                    end
                    testCase.verifyEqual(report.state.stepCount,baseline.state.stepCount);
                    testCase.verifyEqual(report.state.rhsEvaluationCount,baseline.state.rhsEvaluationCount);
                    if family=="constant-nonhydrostatic"
                        initialCheck=testCase.copyFiles(source,"initial-report-check");
                        initialReport=testCase.continueInCpp(initialCheck,"initial",provider,37.5);
                        testCase.verifyEqual(string(initialReport.densityDiagnosticContract.profileRecovery),"dampedLeastSquares");
                        testCase.verifyEqual(string(initialReport.densityDiagnosticContract.selectedProfileRecovery),"not-required");
                    end
                    metrics=report.densityEvaluation;
                    testCase.verifyGreaterThan(metrics.recoveryCount,0);
                    for field=["profileConstructionCount","inversePassCount","apePassCount","apvPassCount"]
                        testCase.verifyEqual(metrics.(field),metrics.recoveryCount,"Each density event shares its expensive stages.");
                    end
                    testCase.verifyEqual(metrics.workspaceLiveBytes,0);
                    testCase.verifyEqual(baseline.densityEvaluation.recoveryCount,0);
                    testCase.verifyEqual(baseline.densityEvaluation.workspaceHighWaterBytes,0);
                    rows{end+1}=struct(configuration=family,provider=provider,errors=errors,profileEvolution=profileEvolution, ...
                        density=metrics,withoutDensity=baseline.densityEvaluation,stepCount=report.state.stepCount,rhsEvaluationCount=report.state.rhsEvaluationCount,coefficientsExactlyEqual=true); %#ok<AGROW>
                    fprintf("DENSITY_NONLINEAR_PASS %s provider=%s profileEvolution=%.17g\n",family,provider,profileEvolution);
                end
            end
            testCase.writeReport("nonlinear",rows);
        end
        function adaptiveDensityOutputAndRestartPreserveIntegration(testCase)
            method="adaptive-rk45"; rows={};
            [source,wvt]=testCase.authorSource("hydrostatic",true,"actual",1,false,true);
            source=testCase.copyFiles(source,"adaptive-density-source");
            control=testCase.copyFiles(source,"adaptive-control");
            testCase.continueInMatlab(control,"actual",38,method);
            [plainSource,~]=testCase.authorSource("hydrostatic",true,"actual",1,false,false);
            for provider=testCase.providers
                actual=testCase.copyFiles(source,"adaptive-density");
                plain=testCase.copyFiles(plainSource,"adaptive-plain");
                report=testCase.continueInCpp(actual,"actual",provider,38,method);
                baseline=testCase.continueInCpp(plain,"actual",provider,38,method);
                errors=testCase.compareFiles(actual,control,wvt,"actual",38);
                for name=["Ap","Am","A0"]
                    for part=["real","imag"]
                        variable="/wave-vortex/"+name+"_"+part;
                        testCase.verifyEqual(ncread(actual,variable),ncread(plain,variable),"Adaptive density output must not modify accepted coefficients.");
                    end
                end
                for count=["stepCount","rejectedStepCount","rhsEvaluationCount"]
                    testCase.verifyEqual(report.state.(count),baseline.state.(count));
                end
                testCase.verifyEqual(report.densityEvaluation.workspaceLiveBytes,0);
                testCase.verifyGreaterThan(report.densityEvaluation.recoveryCount,0);
                testCase.verifyEqual(baseline.densityEvaluation.recoveryCount,0);
                restart=testCase.copyFiles(source,"adaptive-restart");
                testCase.continueInCpp(restart,"actual",provider,37.5,method);
                testCase.verifyRestart(restart,wvt,"actual",37.5);
                testCase.continueInCpp(restart,"actual",provider,38,method);
                restartErrors=testCase.compareFiles(restart,control,wvt,"actual",38);
                rows{end+1}=struct(configuration="hydrostatic",provider=provider,method=method,errors=errors,restartErrors=restartErrors, ...
                    density=report.densityEvaluation,coefficientsExactlyEqual=true,stepCount=report.state.stepCount,rejectedStepCount=report.state.rejectedStepCount,rhsEvaluationCount=report.state.rhsEvaluationCount); %#ok<AGROW>
                fprintf("DENSITY_ADAPTIVE_PASS hydrostatic provider=%s steps=%d rejected=%d\n",provider,report.state.stepCount,report.state.rejectedStepCount);
            end
            testCase.writeReport("adaptive",rows);
        end
        function twoDestinationsContinueInBothDirections(testCase)
            rows={};
            for family=["constant-hydrostatic","hydrostatic"]
                selection="actual";
                if family=="hydrostatic", selection="initial"; end
                [source,wvt]=testCase.authorSource(family,true,selection,2);
                control=testCase.copyFiles(source,"control");
                testCase.continueInMatlab(control,selection,38);
                prefix=testCase.copyFiles(source,"prefix");
                testCase.continueInMatlab(prefix,selection,37.5);
                for provider=testCase.providers
                    for scenario=["whole","cppRestart","matlabToCpp","cppToMatlab"]
                        if scenario=="matlabToCpp", actual=testCase.copyFiles(prefix,"runtime");
                        else, actual=testCase.copyFiles(source,"runtime"); end
                        if ismember(scenario,["cppRestart","cppToMatlab"])
                            testCase.continueInCpp(actual,selection,provider,37.5);
                            for path=actual, testCase.verifyRestart(path,wvt,selection,37.5); end
                        end
                        if scenario=="cppToMatlab", testCase.continueInMatlab(actual,selection,38);
                        else, testCase.continueInCpp(actual,selection,provider,38); end
                        errors=cell(1,2);
                        for index=1:2
                            errors{index}=testCase.compareFiles(actual(index),control(index),wvt,selection,38);
                        end
                        rows{end+1}=struct(configuration=family,reference=selection,provider=provider,scenario=scenario,destinations={errors}); %#ok<AGROW>
                    end
                end
            end
            testCase.writeReport("continuation",rows);
        end
    end
    methods (Access=private)
        function [paths,wvt]=authorSource(testCase,family,antialias,selection,count,isLinear,includeDensity)
            if nargin<6, isLinear=true; end
            if nargin<7, includeDensity=true; end
            grid=[9 7 10]; if antialias, grid=[8 6 9]; end
            wvt=makeTransform(family,grid,antialias);
            n=reshape(1:numel(wvt.A0),size(wvt.A0));
            wave=wvt.J>0 & wvt.Kh>0;
            wvt.Ap=.001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*wave;
            wvt.Am=.002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
            wvt.A0=1e-6*(sin(.7*n)+1i*cos(.3*n));
            meanDensity=wvt.Kh==0 & wvt.J>0;
            wvt.A0(meanDensity)=.003*sin(n(meanDensity));
            wvt.t0=17; wvt.t=37; wvt.removeAllForcing();
            wvt.shouldUseTrueNoMotionProfile=selection=="actual";
            if ~isLinear, wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]); end
            model=WVModel(wvt,shouldUseLinearDynamics=isLinear);
            names=cell(1,0);
            if ~isLinear, names={'u','eta'}; end
            if includeDensity, names=[names,{'rho_nm','eta_true','ape','apv'}]; end
            model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
            paths=fullfile(testCase.folder,"source-"+(1:count)+".nc");
            for index=1:count
                file=model.createNetCDFFileForModelOutput(paths(index),outputInterval=.5,shouldOverwriteExisting=true);
                dense=file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
                dense.addObservingSystem(WVEulerianFields(model,fieldNames=names));
                file.outputTimesForIntegrationPeriod(37,38); file.writeTimeStepToOutputFile(37);
            end
            model.closeNetCDFFile();
            for index=1:count, ncwriteatt(paths(index),"/","portableFileIdentifier","density-output-"+index); end
        end
        function paths=copyFiles(testCase,source,stem)
            paths=fullfile(testCase.folder,stem+"-"+(1:numel(source))+".nc");
            for index=1:numel(source), copyfile(source(index),paths(index)); end
        end
        function continueInMatlab(testCase,paths,selection,finalTime,method)
            if nargin<5, method="fixed-rk4"; end
            for path=paths
                model=WVModel.modelFromFile(char(path));
                cleanup=onCleanup(@()model.closeNetCDFFile());
                testCase.assertTrue(model.wvt.shouldUseTrueNoMotionProfile,"Legacy MATLAB restart intentionally restores the default actual profile.");
                model.wvt.shouldUseTrueNoMotionProfile=selection=="actual";
                if method=="fixed-rk4"
                    model.setupIntegrator(integratorType="fixed",deltaT=.25);
                else
                    model.setupIntegrator(integratorType="adaptive",integrator=@ode45,absTolerance=1e-10,relTolerance=1e-9);
                    model.odeOptions=odeset(model.odeOptions,'InitialStep',.25,'MaxStep',.25);
                end
                model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                testCase.assertEqual(model.wvt.t,finalTime); clear cleanup
            end
        end
        function report=continueInCpp(testCase,paths,selection,provider,finalTime,method)
            if nargin<6, method="fixed-rk4"; end
            requestPath=fullfile(testCase.folder,"request.json"); reportPath=fullfile(testCase.folder,"report.json");
            if method=="fixed-rk4"
                WVModel.writePortableRunRequest(requestPath,paths,method=method,finalTime=finalTime,initialStep=.25,fftProvider=provider,reportPath=reportPath);
            else
                WVModel.writePortableRunRequest(requestPath,paths,method=method,finalTime=finalTime,initialStep=.25,maximumStep=.25,relativeTolerance=1e-9,absoluteToleranceScale=1e-10,fftProvider=provider,reportPath=reportPath);
            end
            request=jsondecode(fileread(requestPath));
            testCase.assertFalse(isfield(request.execution,"densityDiagnostics"));
            if selection=="initial"
                request.execution.densityDiagnostics=struct(contract="wave-vortex-density-diagnostics-v1",reference="initial");
                writeJSON(requestPath,request);
            end
            [status,output]=cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(requestPath));
            testCase.assertEqual(status,0,output);
            report=jsondecode(fileread(reportPath));
            testCase.assertEqual(string(report.status),"complete");
            testCase.assertEqual(string(report.provider.id),provider);
            testCase.assertEqual(string(report.integrationRequest.activeMethod),method);
            testCase.assertTrue(report.integrationRequest.noFallback);
            contract=report.densityDiagnosticContract;
            testCase.assertEqual(string(contract.identifier),"wave-vortex-density-diagnostics-v1");
            testCase.assertEqual(string(contract.reference),selection);
            source="default"; if selection=="initial", source="request"; end
            testCase.assertEqual(string(contract.selectionSource),source);
            testCase.assertEqual(string(contract.outputEvaluation),"available");
            recovery="dampedLeastSquares"; if selection=="initial", recovery="not-required"; end
            testCase.assertEqual(string(contract.selectedProfileRecovery),recovery);
            if report.densityEvaluation.recoveryCount>0, recovery="dampedLeastSquares"; end
            testCase.assertEqual(string(contract.profileRecovery),recovery);
            testCase.assertGreaterThan(report.integrator.denseOutputEvaluationCount,0);
        end
        function errors=compareFiles(testCase,actual,expected,wvt,selection,finalTime)
            errors=struct(rho_nm=0,eta_true=0,ape=0,apv=0);
            for group=["wave-vortex","dense"]
                times=reshape(ncread(actual,"/"+group+"/t"),1,[]);
                interval=.5; if group=="dense", interval=.125; end
                testCase.assertEqual(times,37:interval:finalTime);
                testCase.assertEqual(times,reshape(ncread(expected,"/"+group+"/t"),1,[]));
                for name=["rho_nm","eta_true","ape","apv"]
                    path="/"+group+"/"+name;
                    reference=ncread(expected,path); values=ncread(actual,path);
                    if name=="eta_true"
                        referenceRecords=reshape(reference,numel(wvt.Z),[]);
                        testCase.assertGreaterThan(max(abs(referenceRecords),[],"all"),1e-5);
                        testCase.assertGreaterThan(max(abs(referenceRecords(:,end)-referenceRecords(:,1))),1e-8);
                    end
                    testCase.assertEqual(size(values),size(reference));
                    dimensions=[size(wvt.Z),numel(times)];
                    if name=="rho_nm", dimensions=[wvt.Nz,numel(times)]; end
                    info=ncinfo(actual,path); testCase.assertEqual(info.Size,dimensions);
                    testCase.assertTrue(all(isfinite(values),"all"));
                    tolerance=fieldTolerance(name,reference,wvt);
                    difference=max(abs(values-reference),[],"all");
                    testCase.verifyLessThanOrEqual(difference,tolerance,group+" "+name+" "+selection);
                    errors.(name)=max(errors.(name),difference);
                    if name=="ape", testCase.verifyGreaterThanOrEqual(min(values,[],"all"),0); end
                end
            end
            testCase.verifyRestart(actual,wvt,selection,finalTime);
        end
        function verifyRestart(testCase,path,wvt,selection,time)
            model=WVModel.modelFromFile(char(path)); cleanup=onCleanup(@()model.closeNetCDFFile());
            testCase.assertEqual(model.wvt.t,time); testCase.assertEqual(model.wvt.t0,17);
            testCase.assertTrue(model.wvt.shouldUseTrueNoMotionProfile);
            model.wvt.shouldUseTrueNoMotionProfile=selection=="actual";
            for name=["rho_nm","eta_true","ape","apv"]
                expected=model.wvt.(name); stored=ncread(path,"/wave-vortex/"+name);
                stored=reshape(stored,numel(expected),[]);
                testCase.verifyLessThanOrEqual(max(abs(stored(:,end)-expected(:))),fieldTolerance(name,expected,wvt),"restored "+name);
            end
        end
        function writeReport(~,name,rows)
            directory=string(getenv("WV_DENSITY_OUTPUT_EVIDENCE"));
            if directory~=""
                if ~isfolder(directory), mkdir(directory); end
                writeJSON(fullfile(directory,name+".json"),struct(schema="portable-density-output-v1",matlabRelease=string(version('-release')),rows=[rows{:}]));
            end
        end
    end
end
function tolerance=fieldTolerance(name,reference,wvt)
    switch name
        case "rho_nm", tolerance=1e-6*(max(wvt.rho_total,[],"all")-min(wvt.rho_total,[],"all"));
        case "eta_true", tolerance=1e-6*wvt.Lz;
        case "ape", tolerance=1e-5*max(abs(reference),[],"all")+1e-12;
        case "apv", tolerance=1e-6*max(abs(reference),[],"all")+1e-12;
    end
end
function wvt=makeTransform(family,grid,antialias)
    if family=="hydrostatic"
        wvt=WVTransformHydrostatic([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
    elseif family=="boussinesq"
        wvt=WVTransformBoussinesq([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
    else
        wvt=WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=antialias);
    end
end
function writeJSON(path,value)
    file=fopen(path,"w"); cleanup=onCleanup(@()fclose(file)); fwrite(file,jsonencode(value,PrettyPrint=true));
end
function value=shellQuote(value)
    value="'"+replace(string(value),"'","'""'""'")+"'";
end
function [status,output]=cleanSystem(command)
    [status,output]=system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
