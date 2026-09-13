classdef TestPortableSamplingOutput < matlab.unittest.TestCase
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
            testCase.runner=string(getenv("WV_SAMPLING_OUTPUT_RUNNER"));
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
        function ordinaryAndDenseSamplesMatchMatlab(testCase)
            definitions=[ ...
                struct(family="constant-nonhydrostatic",grid=[8 6 9],reference="actual",method="fixed-rk4"), ...
                struct(family="constant-nonhydrostatic",grid=[9 7 10],reference="initial",method="adaptive-rk45"), ...
                struct(family="hydrostatic",grid=[9 7 10],reference="actual",method="fixed-rk4"), ...
                struct(family="hydrostatic",grid=[8 6 9],reference="initial",method="adaptive-rk45")];
            rows={};
            for definition=definitions
                stem=definition.family+"-"+definition.reference+"-"+definition.method;
                [source,wvt]=testCase.authorSource(definition,1,true,"source-"+stem);
                control=testCase.copyFiles(source,"matlab-"+stem);
                testCase.continueInMatlab(control,definition.reference,38,definition.method);
                for provider=testCase.providers
                    actual=testCase.copyFiles(source,"cpp-"+stem+"-"+provider);
                    report=testCase.continueInCpp(actual,definition.reference,provider,38,definition.method);
                    errors=testCase.compareFiles(actual,control,wvt,definition.reference,38);
                    rows{end+1}=struct(configuration=definition.family,grid=definition.grid,reference=definition.reference,method=definition.method,provider=provider,errors=errors,ordinaryRecords=3,denseRecords=9,contract=report.densityDiagnosticContract); %#ok<AGROW>
                    fprintf("SAMPLING_OUTPUT_PASS %s grid=%s reference=%s method=%s provider=%s\n",definition.family,mat2str(definition.grid),definition.reference,definition.method,provider);
                end
            end
            testCase.writeReport("ordinary-dense",rows);
        end
        function sampledDiagnosticsPreserveIntegration(testCase)
            definition=struct(family="hydrostatic",grid=[9 7 10],reference="actual",method="adaptive-rk45");
            [sampledSource,wvt]=testCase.authorSource(definition,1,true,"sampled-control-source");
            [plainSource,~]=testCase.authorSource(definition,1,false,"plain-control-source");
            matlabControl=testCase.copyFiles(sampledSource,"sampled-matlab");
            testCase.continueInMatlab(matlabControl,definition.reference,38,definition.method);
            rows={};
            for provider=testCase.providers
                sampled=testCase.copyFiles(sampledSource,"sampled-cpp-"+provider);
                plain=testCase.copyFiles(plainSource,"plain-cpp-"+provider);
                report=testCase.continueInCpp(sampled,definition.reference,provider,38,definition.method);
                baseline=testCase.continueInCpp(plain,definition.reference,provider,38,definition.method);
                errors=testCase.compareFiles(sampled,matlabControl,wvt,definition.reference,38);
                for name=["Ap","Am","A0"]
                    for part=["real","imag"]
                        variable="/wave-vortex/"+name+"_"+part;
                        testCase.verifyEqual(ncread(sampled,variable),ncread(plain,variable),"Sampled diagnostics must not change accepted coefficients.");
                    end
                end
                for group=["wave-vortex","dense"]
                    for coordinate=["x","y","z"]
                        variable="/"+group+"/particles_"+coordinate;
                        testCase.verifyEqual(ncread(sampled,variable),ncread(plain,variable),"Sampled diagnostics must not change particle positions.");
                    end
                end
                for count=["stepCount","rejectedStepCount","rhsEvaluationCount"]
                    testCase.verifyEqual(report.state.(count),baseline.state.(count));
                end
                testCase.verifyGreaterThan(report.densityEvaluation.recoveryCount,0);
                testCase.verifyEqual(baseline.densityEvaluation.recoveryCount,0);
                testCase.verifyGreaterThan(report.diagnosticEvaluation.evaluationCount,0);
                testCase.verifyGreaterThan(report.diagnosticEvaluation.primitiveOutputCount,0);
                testCase.verifyEqual(baseline.diagnosticEvaluation.evaluationCount,0);
                rows{end+1}=struct(provider=provider,errors=errors,coefficientsExactlyEqual=true,particlePositionsExactlyEqual=true,stepCount=report.state.stepCount,rejectedStepCount=report.state.rejectedStepCount,rhsEvaluationCount=report.state.rhsEvaluationCount,density=report.densityEvaluation,diagnostics=report.diagnosticEvaluation); %#ok<AGROW>
            end
            testCase.writeReport("integration-control",rows);
        end
        function twoDestinationsContinueInBothDirections(testCase)
            definition=struct(family="hydrostatic",grid=[8 6 9],reference="initial",method="fixed-rk4");
            [source,wvt]=testCase.authorSource(definition,2,true,"restart-source");
            control=testCase.copyFiles(source,"restart-matlab-control");
            testCase.continueInMatlab(control,definition.reference,38,definition.method);
            prefix=testCase.copyFiles(source,"restart-matlab-prefix");
            testCase.continueInMatlab(prefix,definition.reference,37.5,definition.method);
            rows={};
            for provider=testCase.providers
                matlabToCpp=testCase.copyFiles(prefix,"restart-matlab-to-cpp-"+provider);
                report=testCase.continueInCpp(matlabToCpp,definition.reference,provider,38,definition.method);
                matlabToCppErrors=testCase.compareDestinations(matlabToCpp,control,wvt,definition.reference,38);

                cppToMatlab=testCase.copyFiles(source,"restart-cpp-to-matlab-"+provider);
                testCase.continueInCpp(cppToMatlab,definition.reference,provider,37.5,definition.method);
                for path=cppToMatlab, testCase.verifyRestart(path,wvt,definition.reference,37.5); end
                testCase.continueInMatlab(cppToMatlab,definition.reference,38,definition.method);
                cppToMatlabErrors=testCase.compareDestinations(cppToMatlab,control,wvt,definition.reference,38);
                rows{end+1}=struct(provider=provider,reference=definition.reference,method=definition.method,destinations=2,matlabToCpp={matlabToCppErrors},cppToMatlab={cppToMatlabErrors},contract=report.densityDiagnosticContract); %#ok<AGROW>
            end
            testCase.writeReport("continuation",rows);
        end
    end
    methods (Access=private)
        function [paths,wvt]=authorSource(testCase,definition,count,includeSampling,stem)
            wvt=makeTransform(definition.family,definition.grid);
            n=reshape(1:numel(wvt.A0),size(wvt.A0));
            horizontal=wvt.K.^2+wvt.L.^2>0;
            wave=wvt.J>0 & horizontal; inertial=~horizontal;
            wvt.Ap=.001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
            wvt.Am=.002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
            wvt.Am(inertial)=conj(wvt.Ap(inertial));
            wvt.A0=1e-6*(sin(.7*n)+1i*cos(.3*n));
            meanDensity=inertial & wvt.J>0;
            wvt.A0(meanDensity)=.003*sin(n(meanDensity));
            wvt.t0=17; wvt.t=37; wvt.removeAllForcing();
            wvt.shouldUseTrueNoMotionProfile=true; actualDisplacement=wvt.eta_true;
            wvt.shouldUseTrueNoMotionProfile=false; initialDisplacement=wvt.eta_true;
            testCase.assertGreaterThan(max(abs(actualDisplacement-initialDisplacement),[],"all"),1e-8,"Density reference fixture must distinguish actual from initial.");
            wvt.shouldUseTrueNoMotionProfile=definition.reference=="actual";
            wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]);
            model=WVModel(wvt);
            fields={'eta_true','ape','apv','p','zeta_z'};
            if ~includeSampling, fields={}; end
            x=wvt.Lx*[.13 .57 .91]; y=wvt.Ly*[.21 .63 .82]; z=-wvt.Lz*[.18 .47 .79];
            model.addParticles('particles',false,x,y,z,fields{:},trackedVarInterpolation='linear',advectionInterpolation='spline',absToleranceXY=1e-8,absToleranceZ=1e-8);
            paths=fullfile(testCase.folder,stem+"-"+(1:count)+".nc");
            for index=1:count
                file=model.createNetCDFFileForModelOutput(paths(index),outputInterval=.5,shouldOverwriteExisting=true);
                if includeSampling
                    file.outputGroups(1).addObservingSystem(newMooring(model));
                end
                dense=file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=37,finalTime=38);
                dense.addObservingSystem(model.fluxedObservingSystemWithName("particles"));
                if includeSampling, dense.addObservingSystem(newMooring(model)); end
                file.outputTimesForIntegrationPeriod(37,38); file.writeTimeStepToOutputFile(37);
            end
            model.closeNetCDFFile();
            for index=1:count, ncwriteatt(paths(index),"/","portableFileIdentifier","sampling-output-"+index); end
        end
        function paths=copyFiles(testCase,source,stem)
            paths=fullfile(testCase.folder,stem+"-"+(1:numel(source))+".nc");
            for index=1:numel(source), copyfile(source(index),paths(index)); end
        end
        function continueInMatlab(testCase,paths,reference,finalTime,method)
            for path=paths
                model=WVModel.modelFromFile(char(path)); cleanup=onCleanup(@()model.closeNetCDFFile());
                testCase.assertTrue(model.wvt.shouldUseTrueNoMotionProfile,"Legacy MATLAB restart intentionally restores the default actual profile.");
                model.wvt.shouldUseTrueNoMotionProfile=reference=="actual";
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
        function report=continueInCpp(testCase,paths,reference,provider,finalTime,method)
            requestPath=fullfile(testCase.folder,"request.json"); reportPath=fullfile(testCase.folder,"report.json");
            if method=="fixed-rk4"
                WVModel.writePortableRunRequest(requestPath,paths,method=method,finalTime=finalTime,initialStep=.25,fftProvider=provider,reportPath=reportPath);
            else
                WVModel.writePortableRunRequest(requestPath,paths,method=method,finalTime=finalTime,initialStep=.25,maximumStep=.25,relativeTolerance=1e-9,absoluteToleranceScale=1e-10,fftProvider=provider,reportPath=reportPath);
            end
            request=jsondecode(fileread(requestPath));
            testCase.assertFalse(isfield(request.execution,"densityDiagnostics"));
            if reference=="initial"
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
            testCase.assertEqual(string(contract.reference),reference);
            source="default"; if reference=="initial", source="request"; end
            testCase.assertEqual(string(contract.selectionSource),source);
            testCase.assertEqual(string(contract.outputEvaluation),"available");
            recovery="dampedLeastSquares"; if reference=="initial", recovery="not-required"; end
            testCase.assertEqual(string(contract.selectedProfileRecovery),recovery);
            testCase.assertGreaterThan(report.integrator.denseOutputEvaluationCount,0);
        end
        function errors=compareDestinations(testCase,actual,expected,wvt,reference,finalTime)
            errors=cell(1,numel(actual));
            for index=1:numel(actual)
                errors{index}=testCase.compareFiles(actual(index),expected(index),wvt,reference,finalTime);
            end
        end
        function errors=compareFiles(testCase,actual,expected,wvt,reference,finalTime)
            errors=struct;
            for group=["wave-vortex","dense"]
                interval=.5; if group=="dense", interval=.125; end
                times=reshape(ncread(actual,"/"+group+"/t"),1,[]);
                testCase.assertEqual(times,37:interval:finalTime);
                testCase.assertEqual(times,reshape(ncread(expected,"/"+group+"/t"),1,[]));
                for coordinate=["x","y","z"]
                    name="particles_"+coordinate;
                    expectedValues=ncread(expected,"/"+group+"/"+name);
                    values=ncread(actual,"/"+group+"/"+name);
                    difference=max(abs(values-expectedValues),[],"all");
                    tolerance=2e-10*max([1;abs(expectedValues(:))]);
                    testCase.verifyLessThanOrEqual(difference,tolerance,group+" "+name);
                    key=groupKey(group); errors.(key).(name)=difference;
                end
                for observer=["particles","mooring"]
                    for field=["eta_true","ape","apv","p","zeta_z"]
                        name=observer+"_"+field;
                        path="/"+group+"/"+name;
                        expectedValues=ncread(expected,path); values=ncread(actual,path);
                        testCase.assertEqual(size(values),size(expectedValues));
                        testCase.assertTrue(all(isfinite(values),"all"));
                        testCase.assertGreaterThan(max(abs(expectedValues),[],"all"),0,"Fixture field must be nontrivial: "+path);
                        info=ncinfo(actual,path);
                        if observer=="mooring", testCase.assertEqual(info.Size,[wvt.Nz 2 numel(times)]);
                        else, testCase.assertEqual(info.Size,[3 numel(times)]); end
                        tolerance=fieldTolerance(field,expectedValues,wvt);
                        difference=max(abs(values-expectedValues),[],"all");
                        testCase.verifyLessThanOrEqual(difference,tolerance,group+" "+name+" "+reference);
                        key=groupKey(group); errors.(key).(name)=difference;
                    end
                end
            end
            testCase.verifyRestart(actual,wvt,reference,finalTime);
        end
        function verifyRestart(testCase,path,wvt,reference,time)
            model=WVModel.modelFromFile(char(path)); cleanup=onCleanup(@()model.closeNetCDFFile());
            testCase.assertEqual(model.wvt.t,time); testCase.assertEqual(model.wvt.t0,17);
            testCase.assertTrue(model.wvt.shouldUseTrueNoMotionProfile);
            model.wvt.shouldUseTrueNoMotionProfile=reference=="actual";
            file=model.outputFileWithName(string(getFileName(path)));
            ordinary=file.outputGroupWithName("wave-vortex"); dense=file.outputGroupWithName("dense");
            ordinaryParticles=ordinary.observingSystemWithName("particles");
            testCase.verifyTrue(ordinaryParticles==dense.observingSystemWithName("particles"));
            testCase.verifyEqual(reshape(ordinaryParticles.trackedFieldNames,1,[]),["eta_true" "ape" "apv" "p" "zeta_z"]);
            testCase.verifyEqual(ordinaryParticles.advectionInterpolation,'spline');
            testCase.verifyEqual(ordinaryParticles.trackedVarInterpolation,'linear');
            ordinaryMooring=ordinary.observingSystemWithName("mooring"); denseMooring=dense.observingSystemWithName("mooring");
            testCase.verifyFalse(ordinaryMooring==denseMooring,"Each output group must restore its own sample-only mooring.");
            for mooring=[ordinaryMooring,denseMooring]
                testCase.verifyEqual(reshape(mooring.trackedFieldNames,1,[]),["eta_true" "ape" "apv" "p" "zeta_z"]);
                testCase.verifyEqual(mooring.x,model.wvt.Lx*[.17 .71]); testCase.verifyEqual(mooring.y,model.wvt.Ly*[.29 .68]);
            end
            [~,~,~,tracked]=model.particlePositions('particles');
            for field=["eta_true","ape","apv","p","zeta_z"]
                stored=ncread(path,"/wave-vortex/particles_"+field);
                stored=reshape(stored,size(stored,1),[]);
                trackedValues=tracked.(field);
                testCase.verifyLessThanOrEqual(max(abs(stored(:,end)-trackedValues(:))),fieldTolerance(field,trackedValues,wvt),"restored particle "+field);
                mooring=ordinaryMooring;
                values=model.wvt.(field); expected=zeros(model.wvt.Nz,numel(mooring.x));
                for index=1:numel(mooring.x), expected(:,index)=values(mooring.x_index(index),mooring.y_index(index),:); end
                stored=ncread(path,"/wave-vortex/mooring_"+field);
                stored=reshape(stored,size(stored,1),size(stored,2),[]);
                testCase.verifyLessThanOrEqual(max(abs(stored(:,:,end)-expected),[],"all"),fieldTolerance(field,expected,wvt),"restored mooring "+field);
            end
        end
        function writeReport(~,name,rows)
            directory=string(getenv("WV_SAMPLING_OUTPUT_EVIDENCE"));
            if directory~=""
                if ~isfolder(directory), mkdir(directory); end
                writeJSON(fullfile(directory,name+".json"),struct(schema="portable-sampling-output-v1",matlabRelease=string(version('-release')),rows=[rows{:}]));
            end
        end
    end
end
function mooring=newMooring(model)
mooring=WVMooring(model,name="mooring",x=model.wvt.Lx*[.17 .71],y=model.wvt.Ly*[.29 .68],trackedFieldNames={'eta_true','ape','apv','p','zeta_z'});
end
function tolerance=fieldTolerance(name,reference,wvt)
switch name
    case "eta_true", tolerance=1e-6*wvt.Lz;
    case "ape", tolerance=1e-5*max(abs(reference),[],"all")+1e-12;
    case "apv", tolerance=1e-6*max(abs(reference),[],"all")+1e-12;
    otherwise, tolerance=2e-11*max(abs(reference),[],"all")+1e-16;
end
end
function wvt=makeTransform(family,grid)
if family=="hydrostatic"
    wvt=WVTransformHydrostatic([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=true);
else
    wvt=WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=false,shouldAntialias=true);
end
end
function key=groupKey(group)
key=replace(group,"-","_");
end
function name=getFileName(path)
[~,base,extension]=fileparts(path); name=base+extension;
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
