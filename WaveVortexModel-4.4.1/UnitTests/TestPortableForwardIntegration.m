classdef TestPortableForwardIntegration < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        fixtureFolder (1,1) string
        runner (1,1) string
        probe (1,1) string
        providers (1,:) string
        provenance (1,1) struct
    end
    properties (TestParameter)
        configuration = struct(constantHydrostatic="constant-hydrostatic",constantNonhydrostatic="constant-nonhydrostatic", ...
            barotropic="barotropic-qg",stratifiedQG="stratified-qg",hydrostatic="hydrostatic",boussinesq="boussinesq")
    end
    methods (TestClassSetup)
        function prepareRunner(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture(PreservingOnFailure=true));
            testCase.folder = string(fixture.Folder);
            testCase.fixtureFolder = testCase.folder;
            testCase.runner = string(getenv("WV_FORWARD_INTEGRATION_RUNNER"));
            testCase.probe = string(getenv("WV_FORWARD_INTEGRATION_PROBE"));
            if testCase.runner=="", testCase.runner=string(getenv("WV_STABLE_FORCING_RUNNER")); end
            if testCase.probe=="" && testCase.runner~="", testCase.probe=fullfile(fileparts(testCase.runner),"WVForwardIntegrationProbe"); end
            testCase.providers = "reference";
            if getenv("WV_STABLE_FORCING_NATIVE") == "1", testCase.providers = ["reference","native-fftw"]; end
            if testCase.runner == ""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target wave-vortex-run WVForwardIntegrationProbe");
                testCase.assertEqual(status,0,output);
                testCase.runner = fullfile(build,"wave-vortex-run");
                testCase.probe = fullfile(build,"WVForwardIntegrationProbe");
            end
            testCase.assertTrue(isfile(testCase.runner));
            testCase.assertTrue(isfile(testCase.probe));
            sourceRoot=string(getenv("WV_FORWARD_INTEGRATION_SOURCE_ROOT"));
            if sourceRoot=="", sourceRoot=testCase.root; end
            [status,sourceCommit]=cleanSystem("git -C "+shellQuote(sourceRoot)+" rev-parse HEAD");
            testCase.assertEqual(status,0,sourceCommit);
            testCase.provenance=struct(sourceSHA256=sourceHashes(testCase.root),buildSource=struct(root=sourceRoot,commit=string(strtrim(sourceCommit))));
            testCase.provenance.executionBinaries=struct(role={"runner","probe"},path={testCase.runner,testCase.probe}, ...
                sha256={sha256(testCase.runner),sha256(testCase.probe)},sourceCommit=string(strtrim(sourceCommit)));
        end
    end
    methods (Test,TestTags="full")
        function representativeLifecycleMatchesMatlab(testCase,configuration)
            preparationTimer=tic;
            testCase.folder=fullfile(testCase.fixtureFolder,configuration);
            mkdir(testCase.folder);
            definition = profileFor(configuration);
            [source,initial,definition] = testCase.authorModel(definition);
            control = testCase.copyPair(source,"matlab-control");
            prefix = testCase.copyPair(source,"matlab-prefix");
            referenceModels=cell(1,2);
            for destination = 1:2
                testCase.runMatlab(control(destination),definition,definition.finalTime);
                testCase.runMatlab(prefix(destination),definition,definition.splitTime);
                referenceModels{destination}=WVModel.modelFromFile(char(control(destination)));
                referenceModels{destination}.closeNetCDFFile();
            end
            evolution = testCase.measureEvolution(control(1),initial,definition);
            preparationSeconds=toc(preparationTimer);
            for provider = testCase.providers
                providerTimer=tic;
                whole = testCase.copyPair(source,"whole");
                wholeReport = testCase.runCpp(whole,definition,provider,definition.finalTime);
                if definition.method=="fixed-rk4", testCase.assertEqual(wholeReport.integrationRequest.selectedStep,definition.initialStep,RelTol=1e-12); end
                denseHeld = testCase.verifyDenseHeldRecords(whole(1),wholeReport,definition);
                split = testCase.copyPair(source,"split");
                splitReport = testCase.runCpp(split,definition,provider,definition.splitTime);
                hybrid = testCase.copyPair(split,"cpp-matlab");
                for destination = 1:2
                    testCase.verifyRestoredState(hybrid(destination),definition,definition.splitTime);
                    testCase.runMatlab(hybrid(destination),definition,definition.finalTime);
                end
                resumedReport = testCase.runCpp(split,definition,provider,definition.finalTime);
                matlabPrefix = testCase.copyPair(prefix,"matlab-cpp");
                matlabPrefixReport = testCase.runCpp(matlabPrefix,definition,provider,definition.finalTime);
                stopped = testCase.copyPair(source,"stopped");
                stopReport = testCase.runCpp(stopped,definition,provider,definition.finalTime,true);
                testCase.assertEqual(string(stopReport.status),"stopped");
                stopTime = stopReport.termination.finalAcceptedTime;
                testCase.assertEqual(string(stopReport.termination.reason),"stop-requested");
                testCase.assertEqual(string(stopReport.termination.requestedBoundary),definition.stopBoundary);
                testCase.assertGreaterThan(stopTime,definition.initialTime);
                testCase.assertLessThan(stopTime,definition.finalTime);
                stoppedMatlab = testCase.copyPair(stopped,"stopped-matlab");
                for destination = 1:2
                    testCase.verifyRestoredState(stopped(destination),definition,stopTime);
                    testCase.runMatlab(stoppedMatlab(destination),definition,definition.finalTime);
                end
                stoppedResumeReport = testCase.runCpp(stopped,definition,provider,definition.finalTime);
                paths = {whole,split,matlabPrefix,hybrid,stopped,stoppedMatlab};
                names = ["whole","split","matlabToCpp","cppToMatlab","stoppedCpp","stoppedMatlab"];
                outcomes = struct;
                for scenario = 1:numel(names)
                    values = cell(1,2);
                    for destination = 1:2
                        values{destination} = testCase.compareOutput(paths{scenario}(destination),control(destination),definition,referenceModels{destination});
                    end
                    outcomes.(names(scenario)) = [values{:}];
                end
                testCase.assertGreaterThan(wholeReport.integrator.denseOutputEvaluationCount,0);
                checks = struct(name={"nontrivialEvolution","heldAmplitudes","completeRestartState","allContinuationDirections","twoOutputDestinations","controlledStopResume","denseOutputParity"},passed=true);
                row = struct(id=definition.id,configuration=configuration,profile=definition.profile, ...
                    testName="TestPortableForwardIntegration/representativeLifecycleMatchesMatlab(configuration="+definition.parameter+")", ...
                    passed=true,failed=false,incomplete=false,provider=provider,scope="forward-integration", ...
                    controls=definition,evolution=evolution,outcomes=outcomes,denseHeld=denseHeld,stop=stopReport.termination,checks=checks, ...
                    timingSeconds=struct(sharedMatlabPreparation=preparationSeconds,providerLifecycle=toc(providerTimer)), ...
                    runs=struct(whole=compactReport(wholeReport),firstSegment=compactReport(splitReport), ...
                    secondSegment=compactReport(resumedReport),matlabToCpp=compactReport(matlabPrefixReport), ...
                    stopped=compactReport(stopReport),stoppedResume=compactReport(stoppedResumeReport)));
                evidence = string(getenv("WV_FORWARD_INTEGRATION_EVIDENCE"));
                if evidence ~= ""
                    if ~isfolder(evidence), mkdir(evidence); end
                    [status,commit] = cleanSystem("git -C "+shellQuote(testCase.root)+" rev-parse HEAD");
                    testCase.assertEqual(status,0,commit);
                    receipt = struct(schema="wave-vortex-forward-integration-qualification-v1",schemaVersion=1, ...
                        scope="forward-integration",status="complete",sourceCommit=string(strtrim(commit)),provider=provider, ...
                        matlabRelease=string(version('-release')),platform=string(computer),cases=row);
                    receipt.sourceSHA256 = testCase.provenance.sourceSHA256;
                    receipt.buildSource=testCase.provenance.buildSource;
                    receipt.executionBinaries=testCase.provenance.executionBinaries;
                    receipt.artifacts = struct('path',{},'sha256',{});
                    catalog = fullfile(testCase.root,"PortableRuntime","contracts","portable-forward-integration-v1.json");
                    if isfile(catalog), receipt.catalogSHA256=sha256(catalog); end
                    writeJSON(fullfile(evidence,definition.id+"-"+provider+".json"),receipt);
                end
                fprintf("FORWARD_INTEGRATION_PASS %s provider=%s coefficientEvolution=%.17g particleEvolution=%.17g tracerEvolution=%.17g stopTime=%.17g\n", ...
                    configuration,provider,evolution.coefficients,evolution.particlePosition,evolution.tracer,stopTime);
            end
        end
    end
    methods (Static,Hidden)
        function validateExecutableSource(report,executablePath,expectedCommit)
            arguments
                report (1,1) struct
                executablePath (1,1) string
                expectedCommit (1,1) string
            end
            actualCommit = "<missing>";
            if isfield(report,"source") && isstruct(report.source) && isscalar(report.source) && isfield(report.source,"commit")
                value = report.source.commit;
                if (ischar(value) && isrow(value)) || (isstring(value) && isscalar(value))
                    actualCommit = string(value);
                end
            end
            if ismissing(actualCommit) || actualCommit~=expectedCommit || strlength(expectedCommit)==0
                error("WaveVortexModel:ExecutableSourceMismatch", ...
                    "Executable %s reports source revision %s; expected %s. Reconfigure and relink this executable against the declared build source before qualification.", ...
                    executablePath,actualCommit,expectedCommit);
            end
        end
    end
    methods (Access=private)
        function [paths,initial,definition] = authorModel(testCase,definition)
            grid = [8 6 9];
            domain = [17000 11000 1000];
            z = -1000*linspace(1,0,grid(3))'.^1.3;
            N2 = @(z) 1e-4*exp(z/700);
            switch definition.configuration
                case "barotropic-qg"
                    wvt = WVTransformBarotropicQG(domain(1:2),grid(1:2),j=1,shouldAntialias=true);
                case "stratified-qg"
                    wvt = WVTransformStratifiedQG(domain,grid,Nj=4,z=z,N2Function=N2,shouldAntialias=true);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(domain,grid,Nj=4,z=z,N2Function=N2,shouldAntialias=true);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(domain,grid,Nj=4,z=z,N2Function=N2,shouldAntialias=true);
                otherwise
                    wvt = WVTransformConstantStratification(domain,grid,N0=5.2e-3,isHydrostatic=definition.configuration=="constant-hydrostatic",shouldAntialias=true);
            end
            index = reshape(1:numel(wvt.A0),size(wvt.A0));
            wvt.A0 = 1e-6*(sin(.7*index)+1i*cos(.3*index))./(1+wvt.J).*(wvt.Kh>0);
            if wvt.hasWaveComponent
                inertial = wvt.Kh==0;
                wave = wvt.J>0 & ~inertial;
                wvt.Ap = .001*(sin(.7*index)+1i*cos(.3*index))./(1+wvt.J).*(wave|inertial);
                wvt.Am = .002*(sin(.3*index)+1i*cos(.7*index))./(1+wvt.J).*wave;
                wvt.Am(inertial) = conj(wvt.Ap(inertial));
                mda = inertial & wvt.J>0;
                wvt.A0(mda) = .003*sin(index(mda));
            elseif definition.configuration == "barotropic-qg"
                selfConjugate = wvt.dftPrimaryIndices2D==wvt.dftConjugateIndices2D;
                wvt.A0(selfConjugate) = real(wvt.A0(selfConjugate));
            end
            wvt.t0 = -3;
            wvt.t = definition.initialTime;
            wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt), ...
                WVBottomFrictionLinear(wvt,r=1e-5),WVBottomFrictionQuadratic(wvt,Cd=.003),WVBetaPlanePVAdvection(wvt)]);
            heldIndex = find(wvt.Kh>0 & wvt.J==wvt.j(min(2,numel(wvt.j))),1);
            definition.heldIndex = heldIndex;
            fixed = WVFixedAmplitudeForcing(wvt,name="qualification held",A0_indices=uint64(heldIndex),A0bar=wvt.A0(heldIndex));
            wvt.addForcing(fixed);
            wvt.restoreForcingAmplitudes();
            definition.heldReal = real(wvt.A0(heldIndex));
            definition.heldImag = imag(wvt.A0(heldIndex));
            model = WVModel(wvt);
            if definition.method == "fixed-rk4"
                definition.cfl = .5*definition.initialStep/model.timeStepForCFL(.5);
            end
            model.eulerianObservingSystem.addNetCDFOutputVariables('u','v','eta','qgpv');
            xyOnly = ~wvt.hasWaveComponent;
            particleZ = [-600 -100];
            if definition.configuration=="barotropic-qg", particleZ=[]; end
            model.addParticles('drifter',xyOnly,[500 4500],[300 3200],particleZ,'u',advectionInterpolation="spline",trackedVarInterpolation="linear");
            phi = .3+.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly);
            if definition.configuration~="barotropic-qg", phi=phi.*exp(wvt.Z/1000); end
            model.addFluxedObservingSystem(WVTracer(model,name="dye",phi=phi,isXYOnly=xyOnly));
            initial = stateSnapshot(model);
            paths = fullfile(testCase.folder,["source-primary.nc","source-secondary.nc"]);
            for destination = 1:2
                interval = 5;
                if destination==2, interval = 20; end
                file = model.createNetCDFFileForModelOutput(paths(destination),outputInterval=interval,shouldOverwriteExisting=true);
                if definition.configuration~="barotropic-qg"
                    file.outputGroups(1).addObservingSystem(WVMooring(model,name="mooring",x=[0 4000],y=[0 3200],trackedFieldNames={'u'}));
                end
                dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=5,initialTime=definition.initialTime,finalTime=definition.finalTime);
                dense.addObservingSystem(WVEulerianFields(model,fieldNames={'u','eta','qgpv'}));
                file.outputTimesForIntegrationPeriod(definition.initialTime,definition.finalTime);
                file.writeTimeStepToOutputFile(definition.initialTime);
            end
            model.closeNetCDFFile();
            for destination=1:2, ncwriteatt(paths(destination),"/","portableFileIdentifier","qualification-"+destination); end
        end
        function paths = copyPair(testCase,source,stem)
            paths = fullfile(testCase.folder,stem+["-primary.nc","-secondary.nc"]);
            for destination=1:2, copyfile(source(destination),paths(destination)); end
        end
        function runMatlab(testCase,path,definition,finalTime)
            model = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()model.closeNetCDFFile());
            if definition.method == "fixed-rk4"
                model.setupIntegrator(integratorType="fixed",deltaT=definition.initialStep);
            else
                method = str2func("ode"+extractAfter(definition.method,"adaptive-rk"));
                model.setupIntegrator(integratorType="adaptive",integrator=method,absTolerance=definition.absoluteTolerance,relTolerance=definition.relativeTolerance);
                model.odeOptions = odeset(model.odeOptions,'InitialStep',definition.initialStep,'MaxStep',definition.maximumStep);
            end
            model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.assertEqual(model.wvt.t,finalTime,AbsTol=4*eps(finalTime));
            testCase.assertEqual(model.wvt.A0(definition.heldIndex),complex(definition.heldReal,definition.heldImag),AbsTol=1e-18);
        end
        function report = runCpp(testCase,paths,definition,provider,finalTime,shouldStop)
            if nargin<6, shouldStop=false; end
            request = fullfile(testCase.folder,"run.json");
            reportPath = fullfile(testCase.folder,"run-report.json");
            if definition.method == "fixed-rk4"
                WVModel.writePortableRunRequest(request,paths,method=definition.method,finalTime=finalTime,cfl=definition.cfl,timeStepConstraint="min",fftProvider=provider,reportPath=reportPath);
            else
                WVModel.writePortableRunRequest(request,paths,method=definition.method,finalTime=finalTime,initialStep=definition.initialStep,maximumStep=definition.maximumStep, ...
                    relativeTolerance=definition.relativeTolerance,absoluteToleranceScale=definition.absoluteTolerance,fftProvider=provider,reportPath=reportPath);
            end
            executablePath = testCase.runner;
            executable = shellQuote(executablePath);
            if shouldStop
                after = 1;
                if definition.stopBoundary == "output-occurrence", after=2; end
                executablePath = testCase.probe;
                executable = shellQuote(executablePath)+" --stop-boundary "+definition.stopBoundary+" --stop-after "+after;
            end
            [status,output] = cleanSystem(executable+" --request "+shellQuote(request));
            testCase.assertEqual(status,0,output);
            report = jsondecode(fileread(reportPath));
            TestPortableForwardIntegration.validateExecutableSource(report,executablePath,testCase.provenance.buildSource.commit);
            testCase.assertTrue(report.integrationRequest.noFallback);
            testCase.assertTrue(report.execution.noFallback);
            testCase.assertEqual(string(report.provider.id),provider);
            testCase.assertEqual(string(report.integrationRequest.activeMethod),definition.method);
            if definition.method=="fixed-rk4", testCase.assertEqual(string(report.integrationRequest.stepPolicy),"cfl");
            else, testCase.assertEqual(string(report.integrationRequest.stepPolicy),"adaptive"); end
            if ~shouldStop, testCase.assertEqual(string(report.status),"complete"); end
        end
        function model = verifyRestoredState(testCase,path,definition,time)
            model = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()model.closeNetCDFFile());
            testCase.assertEqual(model.wvt.t,time,AbsTol=4*eps(time));
            for name = coefficientNames(model.wvt)
                actual = model.wvt.(name);
                values = complex(ncread(path,"/wave-vortex/"+name+"_real"),ncread(path,"/wave-vortex/"+name+"_imag"));
                values = reshape(values,numel(actual),[]);
                testCase.assertEqual(actual(:),values(:,end));
            end
            [x,y,z] = model.drifterPositions();
            for name=["x","y","z"]
                if name=="z" && definition.configuration=="barotropic-qg", testCase.assertEmpty(z); continue; end
                expected = ncread(path,"/wave-vortex/drifter_"+name);
                actual = x;
                if name=="y", actual=y; elseif name=="z", actual=z; end
                testCase.assertEqual(actual(:),expected(:,end));
            end
            phi = model.tracer("dye");
            saved = reshape(ncread(path,"/wave-vortex/dye"),numel(phi),[]);
            testCase.assertEqual(phi(:),saved(:,end));
            testCase.assertEqual(model.wvt.A0(definition.heldIndex),complex(definition.heldReal,definition.heldImag),AbsTol=1e-18);
        end
        function metrics = verifyDenseHeldRecords(testCase,path,report,definition)
            times=ncread(path,"/wave-vortex/t");
            if definition.method=="fixed-rk4"
                accepted=definition.initialTime+(0:report.state.stepCount)*report.integrationRequest.selectedStep;
                accepted=[accepted,definition.finalTime];
            else
                testCase.assertTrue(report.integrator.acceptedStepDiagnosticsComplete);
                steps=report.integrator.acceptedSteps;
                accepted=[definition.initialTime,[steps.initialTime]+[steps.acceptedStepSize]];
            end
            offStep=min(abs(times(:)-accepted(:)'),[],2)>64*eps(definition.finalTime);
            testCase.assertGreaterThan(nnz(offStep),0,"The coefficient stream must contain actual dense evaluations.");
            values=complex(ncread(path,"/wave-vortex/A0_real"),ncread(path,"/wave-vortex/A0_imag"));
            values=reshape(values,[],numel(times));
            error=max(abs(values(definition.heldIndex,offStep)-complex(definition.heldReal,definition.heldImag)));
            testCase.assertLessThanOrEqual(error,1e-18);
            metrics=struct(recordCount=nnz(offStep),firstRecordTime=times(find(offStep,1)),maximumHeldError=error);
        end
        function metrics = compareOutput(testCase,path,referencePath,definition,reference)
            actual = testCase.verifyRestoredState(path,definition,definition.finalTime);
            metrics = struct(coefficients=0,fields=0,tracer=relativeError(actual.tracer("dye"),reference.tracer("dye")),particlePosition=0,savedOutput=0,heldAmplitude=0,restoredStateExact=true);
            for name=coefficientNames(actual.wvt), metrics.coefficients=max(metrics.coefficients,relativeError(actual.wvt.(name),reference.wvt.(name))); end
            for name=["u","v","eta","qgpv"], metrics.fields=max(metrics.fields,relativeError(actual.wvt.(name),reference.wvt.(name))); end
            [ax,ay,az] = actual.drifterPositions(); [rx,ry,rz] = reference.drifterPositions();
            metrics.particlePosition=max(abs([ax(:)-rx(:);ay(:)-ry(:);az(:)-rz(:)]));
            testCase.assertEqual(actual.isDynamicsLinear,reference.isDynamicsLinear);
            testCase.assertEqual(actual.wvt.forcingNames,reference.wvt.forcingNames);
            for groupName=reshape(reference.outputFiles(1).outputGroupNames,1,[])
                aGroup=actual.outputFiles(1).outputGroupWithName(groupName);
                rGroup=reference.outputFiles(1).outputGroupWithName(groupName);
                testCase.assertEqual(string({aGroup.observingSystems.name}),string({rGroup.observingSystems.name}));
                testCase.assertEqual(string(arrayfun(@class,aGroup.observingSystems,UniformOutput=false)),string(arrayfun(@class,rGroup.observingSystems,UniformOutput=false)));
                testCase.assertEqual(aGroup.incrementsWrittenToGroup,rGroup.incrementsWrittenToGroup);
                testCase.assertEqual([aGroup.outputInterval,aGroup.initialTime,aGroup.finalTime],[rGroup.outputInterval,rGroup.initialTime,rGroup.finalTime]);
            end
            for group=["wave-vortex","dense"]
                testCase.assertEqual(ncread(path,"/"+group+"/t"),ncread(referencePath,"/"+group+"/t"));
                info=ncinfo(referencePath,"/"+group);
                for variable=reshape(info.Variables,1,[])
                    name=string(variable.Name);
                    if ~any(strcmp(name,["u","v","eta","qgpv","dye","drifter_x","drifter_y","drifter_z","drifter_u","drifter_qgpv","mooring_u","mooring_qgpv","Ap_real","Ap_imag","Am_real","Am_imag","A0_real","A0_imag"])), continue; end
                    observed=ncread(path,"/"+group+"/"+name); expected=ncread(referencePath,"/"+group+"/"+name);
                    testCase.assertSize(observed,size(expected));
                    testCase.assertTrue(all(isfinite(observed(:))));
                    metrics.savedOutput=max(metrics.savedOutput,relativeError(observed,expected));
                end
            end
            held=complex(ncread(path,"/wave-vortex/A0_real"),ncread(path,"/wave-vortex/A0_imag"));
            times=ncread(path,"/wave-vortex/t"); held=reshape(held,[],numel(times));
            metrics.heldAmplitude=max(abs(held(definition.heldIndex,:)-complex(definition.heldReal,definition.heldImag)));
            testCase.assertLessThanOrEqual(metrics.heldAmplitude,1e-18);
            testCase.assertLessThanOrEqual(max([metrics.coefficients,metrics.fields,metrics.tracer,metrics.savedOutput]),2e-7,path+" "+jsonencode(metrics));
            testCase.assertLessThanOrEqual(metrics.particlePosition,.02);
        end
        function evolution = measureEvolution(testCase,path,initial,definition)
            model=WVModel.modelFromFile(char(path)); cleanup=onCleanup(@()model.closeNetCDFFile());
            final=stateSnapshot(model);
            final.coefficients(definition.heldIndex)=0;
            initial.coefficients(definition.heldIndex)=0;
            evolution=struct(coefficients=relativeError(final.coefficients,initial.coefficients), ...
                particlePosition=max(abs(final.positions-initial.positions)),tracer=relativeError(final.tracer,initial.tracer));
            testCase.assertGreaterThan(evolution.coefficients,1e-6);
            testCase.assertGreaterThan(evolution.particlePosition,1e-4);
            testCase.assertGreaterThan(evolution.tracer,1e-8);
        end
    end
end

function definition=profileFor(configuration)
    configurations=["constant-hydrostatic","constant-nonhydrostatic","barotropic-qg","stratified-qg","hydrostatic","boussinesq"];
    parameters=["constantHydrostatic","constantNonhydrostatic","barotropic","stratifiedQG","hydrostatic","boussinesq"];
    method="adaptive-rk45"; profile="adaptive-rk45"; boundary="output-occurrence";
    if configuration=="constant-hydrostatic", method="fixed-rk4"; profile="rk4-cfl"; boundary="accepted-step"; end
    if configuration=="constant-nonhydrostatic", method="adaptive-rk23"; profile="adaptive-rk23"; end
    if ismember(configuration,["barotropic-qg","hydrostatic"]), boundary="accepted-step"; end
    suffix=extractAfter(method,"adaptive-");
    if method=="fixed-rk4", suffix="cfl-rk4"; end
    definition=struct(id="forward-"+configuration+"-"+suffix,configuration=configuration,parameter=parameters(configurations==configuration),profile=profile,method=method, ...
        initialTime=17,splitTime=77,finalTime=137,initialStep=4,maximumStep=4,relativeTolerance=1e-9,absoluteTolerance=1e-10,stopBoundary=boundary);
end
function snapshot=stateSnapshot(model)
    coefficients=model.wvt.A0(:);
    if model.wvt.hasWaveComponent, coefficients=[coefficients;model.wvt.Ap(:);model.wvt.Am(:)]; end
    [x,y,z]=model.drifterPositions();
    snapshot=struct(coefficients=coefficients,positions=[x(:);y(:);z(:)],tracer=model.tracer("dye"));
end
function names=coefficientNames(wvt)
    names="A0";
    if wvt.hasWaveComponent, names=["Ap","Am","A0"]; end
end
function report=compactReport(full)
    report=struct(status=full.status,integrationRequest=full.integrationRequest,termination=full.termination, ...
        state=full.state,acceptedStepCount=full.state.stepCount,rejectedStepCount=full.state.rejectedStepCount, ...
        denseOutputEvaluationCount=full.integrator.denseOutputEvaluationCount);
    if isfield(full,"variableKernelPolicy"), report.variableKernelPolicy=full.variableKernelPolicy; end
end
function value=relativeError(actual,expected)
    value=max(abs(actual(:)-expected(:)))/max(max(abs(expected(:))),realmin);
end
function value=sha256(path)
    [status,output]=system("shasum -a 256 "+shellQuote(path)); assert(status==0,output);
    value=extractBefore(string(output)," ");
end
function hashes=sourceHashes(root)
    paths=["PortableRuntime/include/WaveVortexRuntime/WVRungeKutta.hpp","PortableRuntime/include/WaveVortexRuntime/WVIntegrationContracts.hpp","PortableRuntime/src/WVRungeKutta.cpp","PortableRuntime/src/WVOrderedRKCombination.hpp", ...
        "PortableRuntime/src/WVModel.cpp","PortableRuntime/src/WVOutputOrchestration.cpp", ...
        "PortableRuntime/src/WVModelOutputNetCDFWriter.cpp","PortableRuntime/app/WaveVortexRun.cpp", ...
        "PortableRuntime/src/WVFieldEvaluationService.cpp","PortableRuntime/src/WVStratifiedFieldEvaluationAdapter.cpp", ...
        "PortableRuntime/src/WVBarotropicQGFieldEvaluationAdapter.cpp", ...
        "CompiledKernel/src/WVTransformConstantStratificationKernel.cpp", ...
        "PortableRuntime/tests/WVForwardIntegrationProbe.cpp","UnitTests/TestPortableForwardIntegration.m"];
    hashes=struct('path',{},'sha256',{});
    sourceRoot=string(getenv("WV_FORWARD_INTEGRATION_SOURCE_ROOT"));
    if sourceRoot=="", sourceRoot=root; end
    runnerSource=fullfile(sourceRoot,"PortableRuntime/app/WaveVortexRun.cpp");
    if isfile(runnerSource) && contains(fileread(runnerSource),"WVRunnerVariablePolicy.hpp")
        paths=[paths,"PortableRuntime/CMakeLists.txt","CompiledKernel/CMakeLists.txt", ...
            "PortableRuntime/app/WVRunnerVariablePolicy.cpp","PortableRuntime/app/WVRunnerVariablePolicy.hpp"];
    end
    for index=1:numel(paths)
        file=fullfile(sourceRoot,paths(index));
        if startsWith(paths(index),"UnitTests/"), file=fullfile(root,paths(index)); end
        hashes(index)=struct(path=paths(index),sha256=sha256(file));
    end
end
function writeJSON(path,value)
    file=fopen(path,"w"); cleanup=onCleanup(@()fclose(file)); fwrite(file,jsonencode(value,PrettyPrint=true));
end
function value=shellQuote(value)
    value="'"+replace(string(value),"'","'""'""'")+"'";
end
function [status,output]=cleanSystem(command)
    command="env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command;
    [status,output]=system(command);
end
