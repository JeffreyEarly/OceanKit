classdef TestPortableStratifiedQG < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        executable (1,1) string
        runner (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function buildRuntime(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_STABLE_FORCING_DUMP"));
            testCase.runner = string(getenv("WV_STABLE_FORCING_RUNNER"));
            testCase.providers = "reference";
            if getenv("WV_STABLE_FORCING_NATIVE") == "1", testCase.providers = ["reference","native"]; end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 4 --target WVStableForcingDump WVStratifiedQGFieldDump wave-vortex-run");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVStableForcingDump");
                testCase.runner = fullfile(build,"wave-vortex-run");
            end
            testCase.assertTrue(isfile(testCase.executable));
        end
    end
    methods (Test,TestTags="full")
        function linearPassiveObserversMatchMatlab(testCase)
            wvt = testCase.transform([8 6 9],true);
            wvt.A0(:,wvt.k==0 & wvt.l==0) = 1e-6*(1+2i);
            source = testCase.writeInitialModel(wvt,true,true);
            for provider = testCase.providers
                path = fullfile(testCase.folder,"linear.nc"); copyfile(source,path);
                controlPath = fullfile(testCase.folder,"linear-control.nc"); copyfile(source,controlPath);
                request = fullfile(testCase.folder,"linear.json");
                WVModel.writePortableRunRequest(request,path,method="fixed-rk4",finalTime=38,initialStep=.25,fftProvider=replace(provider,"native","native-fftw"));
                [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(request)); testCase.assertEqual(status,0,output);
                control = WVModel.modelFromFile(char(controlPath)); cleanup = onCleanup(@()control.closeNetCDFFile());
                control.setupIntegrator(integratorType="fixed",deltaT=.25);
                control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                expectedTracer = control.tracer("dye"); clear cleanup
                actual = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()actual.closeNetCDFFile());
                testCase.verifyTrue(actual.isDynamicsLinear); testCase.verifyEqual(actual.wvt.A0,wvt.A0);
                testCase.verifyEqual(actual.tracer("dye"),expectedTracer,AbsTol=1e-12); clear cleanup
            end
        end
        function cflAndDefaultSteppingMatchMatlab(testCase)
            wvt = testCase.transform([9 7 9],false);
            wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt),testCase.forcing(wvt,"WVFixedAmplitudeForcing")]);
            source = testCase.writeInitialModel(wvt,true);
            model = WVModel(wvt); expectedStep = model.timeStepForCFL(.5);
            for provider = testCase.providers
                for mode = ["cfl","default-rk23","default"]
                    path = fullfile(testCase.folder,"stepping.nc"); copyfile(source,path);
                    controlPath = fullfile(testCase.folder,"stepping-control.nc"); copyfile(source,controlPath);
                    request = fullfile(testCase.folder,"stepping.json"); reportPath = fullfile(testCase.folder,"stepping-report.json");
                    if mode == "cfl"
                        WVModel.writePortableRunRequest(request,path,method="fixed-rk4",cfl=.5,timeStepConstraint="advective",finalTime=38,fftProvider=replace(provider,"native","native-fftw"),reportPath=reportPath);
                    elseif mode == "default-rk23"
                        WVModel.writePortableRunRequest(request,path,method="adaptive-rk23",finalTime=38,fftProvider=replace(provider,"native","native-fftw"),reportPath=reportPath);
                    else
                        WVModel.writePortableRunRequest(request,path,finalTime=38,fftProvider=replace(provider,"native","native-fftw"),reportPath=reportPath);
                    end
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(request)); testCase.assertEqual(status,0,output);
                    report = jsondecode(fileread(reportPath));
                    testCase.verifyEqual(report.integrationRequest.selectedStep,expectedStep,RelTol=1e-12);
                    testCase.verifyTrue(report.execution.noFallback);
                    control = WVModel.modelFromFile(char(controlPath)); cleanup = onCleanup(@()control.closeNetCDFFile());
                    if mode == "cfl", control.setupIntegrator(integratorType="fixed",deltaT=expectedStep);
                    elseif mode == "default-rk23", control.setupIntegrator(integratorType="adaptive",integrator=@ode23);
                    else, control.setupIntegrator(); end
                    control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                    expected = control.wvt.A0; expectedTracer = control.tracer("dye"); clear cleanup
                    actual = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()actual.closeNetCDFFile());
                    testCase.verifyEqual(actual.wvt.A0,expected,RelTol=2e-11,AbsTol=1e-16);
                    testCase.verifyEqual(actual.tracer("dye"),expectedTracer,AbsTol=1e-12); clear cleanup
                end
            end
        end
        function multiFileOutputPoliciesAndRestart(testCase)
            wvt = testCase.transform([8 6 9],true);
            wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]);
            [primary,secondary] = testCase.writeInitialModel(wvt,true);
            ncwriteatt(primary,"/","portableFileIdentifier","primary"); ncwriteatt(secondary,"/","portableFileIdentifier","secondary");
            controlPath = fullfile(testCase.folder,"policy-control.nc"); copyfile(primary,controlPath);
            control = WVModel.modelFromFile(char(controlPath)); cleanup = onCleanup(@()control.closeNetCDFFile());
            control.setupIntegrator(integratorType="fixed",deltaT=.25);
            control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            expected = control.wvt.A0; expectedTracer = control.tracer("dye"); clear cleanup
            for provider = testCase.providers
                for policy = ["create","replace","append"]
                    inputs = [fullfile(testCase.folder,"policy-primary.nc"),fullfile(testCase.folder,"policy-secondary.nc")];
                    copyfile(primary,inputs(1)); copyfile(secondary,inputs(2));
                    scientific = cell(2,2);
                    for index=1:2
                        scientific{index,1}=ncread(inputs(index),"PF0inv");
                        scientific{index,2}=ncread(inputs(index),"N2Function");
                    end
                    outputs = inputs; destinations = configureDictionary("string","string");
                    if policy ~= "append"
                        outputs = [fullfile(testCase.folder,"output-primary.nc"),fullfile(testCase.folder,"output-secondary.nc")];
                        for index=1:2
                            if isfile(outputs(index)), delete(outputs(index)); end
                            if policy=="replace", copyfile(inputs(index),outputs(index)); end
                        end
                        destinations(["primary","secondary"]) = outputs;
                    end
                    request = fullfile(testCase.folder,"policy.json");
                    WVModel.writePortableRunRequest(request,inputs,method="fixed-rk4",finalTime=37.5,initialStep=.25,outputPolicy=policy,destinations=destinations,fftProvider=replace(provider,"native","native-fftw"));
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(request)); testCase.assertEqual(status,0,output);
                    WVModel.writePortableRunRequest(request,outputs,method="fixed-rk4",finalTime=38,initialStep=.25,outputPolicy="append",fftProvider=replace(provider,"native","native-fftw"));
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(request)); testCase.assertEqual(status,0,output);
                    % Compare C++ output before writable MATLAB reload can
                    % serialize the opaque function handle again. Saved input
                    % bytes also make the append assertion non-tautological.
                    % New files use the selected restart's immutable source
                    % (primary at equal times). Append preserves each file's
                    % original payload; independent MATLAB serializations of
                    % the same function need not have identical opaque bytes.
                    for index=1:2
                        sourceIndex=1; if policy=="append", sourceIndex=index; end
                        testCase.verifyEqual(ncread(outputs(index),"PF0inv"),scientific{sourceIndex,1});
                        testCase.verifyEqual(ncread(outputs(index),"N2Function"),scientific{sourceIndex,2});
                    end
                    actual = WVModel.modelFromFile(char(outputs(1))); cleanup = onCleanup(@()actual.closeNetCDFFile());
                    testCase.verifyEqual(actual.wvt.A0,expected,RelTol=2e-11,AbsTol=1e-16);
                    testCase.verifyEqual(actual.tracer("dye"),expectedTracer,AbsTol=1e-12); clear cleanup
                    times = ncread(outputs(1),"/dense/t");
                    controlTimes = ncread(controlPath,"/dense/t");
                    expectedTimes = controlTimes;
                    if policy ~= "append", expectedTimes = expectedTimes(expectedTimes > 37); end
                    testCase.verifyEqual(times,expectedTimes);
                    [found,indices] = ismember(times,controlTimes); testCase.assertTrue(all(found));
                    expectedField = ncread(controlPath,"/dense/u"); observedField = ncread(outputs(1),"/dense/u");
                    testCase.verifyEqual(observedField,expectedField(:,:,:,indices),AbsTol=1e-12,RelTol=2e-11);
                    expectedTimes = (37:.25:38)';
                    if policy ~= "append", expectedTimes = expectedTimes(2:end); end
                    testCase.verifyEqual(ncread(outputs(2),"/wave-vortex/t"),expectedTimes);
                end
            end
        end
        function incompatibleModelGraphsAreTransactional(testCase)
            catalog = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
            rows = catalog.rows(startsWith(string({catalog.rows.configuration}),"stratified-qg") & string({catalog.rows.acceptance}) == "incompatible");
            testCase.assertNumElements(rows,7);
            for row = reshape(rows,1,[])
                wvt = testCase.transform([8 6 9],endsWith(string(row.configuration),"aa1"));
                wvt.setForcing(WVNonlinearAdvection(wvt));
                path = testCase.writeInitialModel(wvt);
                ncwriteatt(path,"/forcing","AnnotatedClass",row.forcing);
                before = fileBytes(path);
                [status,output] = cleanSystem(shellQuote(testCase.runner)+" "+shellQuote(path)+" --restart-mode model --output-policy append --integrator fixed-rk4 --delta-t .25 --final-time 38 --fft-provider reference");
                testCase.verifyNotEqual(status,0,output);
                testCase.verifyEqual(fileBytes(path),before,string(row.id));
            end
            wvt = testCase.transform([8 6 9],false);
            wvt.setForcing(WVNonlinearAdvection(wvt));
            path = testCase.writeInitialModel(wvt,true,false,false);
            before = fileBytes(path);
            [status,output] = cleanSystem(shellQuote(testCase.runner)+" "+shellQuote(path)+" --restart-mode model --output-policy append --integrator fixed-rk4 --delta-t .25 --final-time 38 --fft-provider reference");
            testCase.verifyNotEqual(status,0,output);
            testCase.verifyTrue(contains(output,"XY-only advection"),output);
            testCase.verifyEqual(fileBytes(path),before);
        end
        function modelIntegrationRoundTripsThroughMatlab(testCase)
            wvt = testCase.transform([9 7 9],false);
            wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt),WVVerticalDiffusivity(wvt,kappa_z=.002),WVBottomFrictionLinear(wvt,r=1e-5),WVBottomFrictionQuadratic(wvt,Cd=.003),WVBetaPlanePVAdvection(wvt)]);
            source = testCase.writeInitialModel(wvt,true);
            for provider = testCase.providers
                for method = ["fixed-rk4","adaptive-rk23","adaptive-rk78"]
                    runtimePath = fullfile(testCase.folder,"runtime.nc"); copyfile(source,runtimePath);
                    requestPath = fullfile(testCase.folder,"request.json");
                    if method == "fixed-rk4"
                        WVModel.writePortableRunRequest(requestPath,runtimePath,method=method,finalTime=38,initialStep=.25,fftProvider=replace(provider,"native","native-fftw"));
                    else
                        WVModel.writePortableRunRequest(requestPath,runtimePath,method=method,finalTime=38,initialStep=.5,maximumStep=.5,relativeTolerance=1e-8,absoluteToleranceScale=1e-10,fftProvider=replace(provider,"native","native-fftw"));
                    end
                    [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(requestPath));
                    testCase.assertEqual(status,0,output);
                    actual = WVModel.modelFromFile(char(runtimePath)); actualCleanup = onCleanup(@()actual.closeNetCDFFile());
                    controlPath = fullfile(testCase.folder,"control.nc"); copyfile(source,controlPath);
                    control = WVModel.modelFromFile(char(controlPath)); controlCleanup = onCleanup(@()control.closeNetCDFFile());
                    if method == "fixed-rk4"
                        control.setupIntegrator(integratorType="fixed",deltaT=.25);
                    else
                        control.setupIntegrator(integratorType="adaptive",integrator=str2func("ode"+extractAfter(method,"adaptive-rk")),absTolerance=1e-10,relTolerance=1e-8);
                        control.odeOptions = odeset(control.odeOptions,'InitialStep',.5,'MaxStep',.5);
                    end
                    control.integrateToTime(38,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                    testCase.verifyLessThanOrEqual(max(abs(actual.wvt.A0-control.wvt.A0),[],"all"),2e-11*max(abs(control.wvt.A0),[],"all")+1e-16,provider+" "+method);
                    testCase.verifyEqual(actual.tracer("dye"),control.tracer("dye"),AbsTol=1e-12);
                    [ax,ay,az,at] = actual.drifterPositions(); [cx,cy,cz,ct] = control.drifterPositions();
                    testCase.verifyEqual(ax,cx,AbsTol=1e-9); testCase.verifyEqual(ay,cy,AbsTol=1e-9); testCase.verifyEqual(az,cz);
                    testCase.verifyEqual(at.u,ct.u,AbsTol=1e-12); testCase.verifyEqual(at.qgpv,ct.qgpv,AbsTol=1e-14);
                    for name = ["u","v","eta","p","rho_e","qgpv","zeta_z","ssu","ssv","ssh","dye","drifter_x","drifter_y","mooring_u","mooring_qgpv"]
                        expected = ncread(controlPath,"/wave-vortex/"+name); observed = ncread(runtimePath,"/wave-vortex/"+name);
                        testCase.verifyEqual(observed,expected,AbsTol=1e-10,RelTol=2e-11);
                    end
                    clear actualCleanup controlCleanup
                end
            end
        end
        function fieldSamplingMatchesMatlab(testCase)
            z = -1000*linspace(1,0,9)'.^1.4;
            wvt = WVTransformStratifiedQG([17000 11000 1000],[9 7 9],Nj=4,z=z,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=false);
            n = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
            wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J);
            wvt.A0(:,wvt.k==0 & wvt.l==0) = 0;
            source = testCase.writeInitialModel(wvt);
            request = struct(x=[-1;17100;4300;8500;1000;3500],y=[40;-10;11100;800;5500;9000],z=[-1100;-1000;-510;-200;0;1],fields={{'u','v','eta','p','rho_e','rho_total','qgpv','zeta_z','ssu','ssv','ssh'}});
            requestPath = fullfile(testCase.folder,"positions.json");
            file = fopen(requestPath,"w"); cleanup = onCleanup(@()fclose(file)); fwrite(file,jsonencode(request)); clear cleanup
            executable = fullfile(fileparts(testCase.executable),"WVStratifiedQGFieldDump");
            outputPath = fullfile(testCase.folder,"fields.json");
            for provider = testCase.providers
                [status,output] = cleanSystem(shellQuote(executable)+" "+shellQuote(source)+" "+shellQuote(requestPath)+" "+shellQuote(outputPath)+" "+provider);
                testCase.assertEqual(status,0,output);
                actual = jsondecode(fileread(outputPath));
                for method = ["linear","spline"]
                    for name = string(request.fields)
                        expected = wvt.variableAtPositionWithName(request.x,request.y,request.z,char(name),interpolationMethod=char(method));
                        for path = ["fixed","moving","event"]
                            observed = actual.(method).(path).(name);
                            testCase.verifyLessThanOrEqual(max(abs(observed(:)-expected(:))),2e-11*max(abs(expected(:)))+1e-16,provider+" "+method+" "+path+" "+name);
                        end
                    end
                    if provider=="native", testCase.verifyEqual(actual.(method).allocations,0); end
                end
            end
        end
        function applicableForcingTendenciesMatchMatlab(testCase)
            catalog = jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-forcing-compatibility-v1.json")));
            for grid = {[8 6 9],[9 7 10]}
                for antialias = [false true]
                    wvt = testCase.transform(grid{1},antialias);
                    rows = catalog.rows(string({catalog.rows.configuration}) == "stratified-qg-aa"+double(antialias) & string({catalog.rows.acceptance}) == "supported");
                    for row = reshape(rows,1,[])
                        identity = string(row.forcing);
                        wvt.setForcing(testCase.forcing(wvt,identity));
                        testCase.compareTendency(wvt);
                        if grid{1}(1)==8, testCase.compareContinuation(wvt); end
                    end
                end
            end
        end
        function orderedForcingTendenciesMatchMatlab(testCase)
            for retained = [1 2 4]
                wvt = testCase.transform([9 7 9],false,Nj=retained);
                for forceMean = [false true]
                    forces = [WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt),WVAntialiasing(wvt,Nj=max(1,retained-1)),WVVerticalDiffusivity(wvt,kappa_z=.002,shouldForceMeanDensityAnomaly=forceMean),WVBottomFrictionLinear(wvt,r=1e-5),WVBottomFrictionQuadratic(wvt,Cd=.003),WVBetaPlanePVAdvection(wvt),testCase.forcing(wvt,"WVFixedAmplitudeForcing")];
                    for reverse = [false true]
                        if reverse, forces=fliplr(forces); end
                        wvt.setForcing(forces);
                        testCase.compareTendency(wvt);
                    end
                end
            end
        end
    end
    methods (Access=private)
        function wvt = transform(~,grid,antialias,options)
            arguments
                ~
                grid
                antialias
                options.Nj = 4
            end
            wvt = WVTransformStratifiedQG([17000 11000 1000],grid,Nj=options.Nj,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=antialias);
            n = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
            wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J);
            wvt.A0(:,wvt.k==0 & wvt.l==0) = 0;
            wvt.t = 37;
            wvt.removeAllForcing();
        end
        function force = forcing(~,wvt,identity)
            switch identity
                case "WVFixedAmplitudeForcing"
                    index = uint64(wvt.Nj+1);
                    force = WVFixedAmplitudeForcing(wvt,name="catalog fixed",A0_indices=index,A0bar=wvt.A0(index));
                case "WVNarrowBandGeostrophicForcing"
                    force = WVNarrowBandGeostrophicForcing(wvt,initialPV="none",k_f=2*wvt.dk,j_f=0);
                case "WVAntialiasing"
                    force = WVAntialiasing(wvt,Nj=2);
                case "WVVerticalDiffusivity"
                    force = WVVerticalDiffusivity(wvt,kappa_z=.002);
                otherwise
                    force = feval(identity,wvt);
            end
        end
        function [path,secondary] = writeInitialModel(testCase,wvt,withObservers,linear,tracerIsXYOnly)
            if nargin < 3, withObservers = false; end
            if nargin < 4, linear = false; end
            if nargin < 5, tracerIsXYOnly = true; end
            model = WVModel(wvt,shouldUseLinearDynamics=linear);
            if withObservers
                model.eulerianObservingSystem.addNetCDFOutputVariables('u','v','eta','p','rho_e','qgpv','zeta_z','ssu','ssv','ssh');
                model.setDrifterPositions([500 16500],[300 10900],[-750 -100],'u','qgpv',advectionInterpolation="linear",trackedVarInterpolation="spline");
                tracer = WVTracer(model,name="dye",phi=.3+.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly).*exp(wvt.Z/1000),isXYOnly=tracerIsXYOnly);
                model.addFluxedObservingSystem(tracer);
            end
            path = fullfile(testCase.folder,"source.nc");
            output = model.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true);
            if withObservers
                output.outputGroups(1).addObservingSystem(WVMooring(model,name="mooring",x=[0 4500],y=[0 3400],trackedFieldNames={'u','qgpv'}));
                dense = output.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=wvt.t,finalTime=wvt.t+2);
                dense.addObservingSystem(WVEulerianFields(model,fieldNames={'u','eta'}));
            end
            secondary = "";
            if nargout > 1
                secondary = fullfile(testCase.folder,"secondary.nc");
                model.createNetCDFFileForModelOutput(secondary,outputInterval=.25,shouldOverwriteExisting=true);
            end
            for file = reshape(model.outputFiles,1,[])
                file.outputTimesForIntegrationPeriod(wvt.t,wvt.t+2);
                file.writeTimeStepToOutputFile(wvt.t);
            end
            model.closeNetCDFFile();
        end
        function compareContinuation(testCase,wvt)
            source = testCase.writeInitialModel(wvt);
            controlPath = fullfile(testCase.folder,"pair-control.nc"); copyfile(source,controlPath);
            control = WVModel.modelFromFile(char(controlPath)); cleanup = onCleanup(@()control.closeNetCDFFile());
            control.setupIntegrator(integratorType="fixed",deltaT=.25);
            control.integrateToTime(wvt.t+.5,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            expected = control.wvt.A0; clear cleanup
            for provider = testCase.providers
                runtimePath = fullfile(testCase.folder,"pair-runtime.nc"); copyfile(source,runtimePath);
                request = fullfile(testCase.folder,"pair.json");
                WVModel.writePortableRunRequest(request,runtimePath,method="fixed-rk4",finalTime=wvt.t+.5,initialStep=.25,fftProvider=replace(provider,"native","native-fftw"));
                [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(request)); testCase.assertEqual(status,0,output);
                actual = WVModel.modelFromFile(char(runtimePath)); cleanup = onCleanup(@()actual.closeNetCDFFile());
                testCase.verifyEqual(actual.wvt.A0,expected,string(class(wvt.forcing)),RelTol=2e-11,AbsTol=1e-16); clear cleanup
            end
        end
        function compareTendency(testCase,wvt)
            expected = wvt.nonlinearFlux();
            source = testCase.writeInitialModel(wvt);
            outputPath = fullfile(testCase.folder,"rhs.json");
            for provider = testCase.providers
                [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(source)+" "+shellQuote(outputPath)+" "+provider);
                testCase.assertEqual(status,0,output);
                report = jsondecode(fileread(outputPath));
                actual = complex(report.F0.real,report.F0.imag);
                testCase.verifyLessThanOrEqual(max(abs(actual(:)-expected(:))),2e-11*max(abs(expected(:)))+1e-18,provider);
                if provider == "native", testCase.verifyEqual(report.applicationAllocations,0); end
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

function value = fileBytes(path)
file = fopen(path,"rb"); cleanup = onCleanup(@()fclose(file)); value = fread(file,Inf,"*uint8");
end
