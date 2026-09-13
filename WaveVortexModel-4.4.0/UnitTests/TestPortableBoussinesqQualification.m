classdef TestPortableBoussinesqQualification < matlab.unittest.TestCase
    properties (SetAccess=private)
        folder (1,1) string
        runner (1,1) string
        providers (1,:) string
        evidenceFolder (1,1) string
    end
    methods (TestClassSetup)
        function prepareRuntime(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.evidenceFolder = string(getenv("WV_BOUSS_EVIDENCE_FOLDER"));
            if testCase.evidenceFolder == "", testCase.evidenceFolder = testCase.folder; end
            testCase.runner = string(getenv("WV_STABLE_FORCING_RUNNER"));
            testCase.providers = "reference";
            if getenv("WV_STABLE_FORCING_NATIVE") == "1", testCase.providers = ["reference","native-fftw"]; end
            if testCase.runner == ""
                root = string(fileparts(fileparts(mfilename("fullpath"))));
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 4 --target wave-vortex-run WVBoussinesqLifecycleProbe");
                testCase.assertEqual(status,0,output);
                testCase.runner = fullfile(build,"wave-vortex-run");
            end
        end
    end
    methods (Test,TestTags="full")
        function longerContinuationMatchesMatlab(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            manifest = jsondecode(fileread(fullfile(root,"PortableRuntime","contracts","boussinesq-qualification-cases-v1.json")));
            for definition = reshape(manifest.cases,1,[])
                definition.grid = reshape(definition.grid,1,[]);
                [source,initialState,verticalKineticFraction] = testCase.authorModel(definition);
                for provider = testCase.providers
                    paths = fullfile(testCase.folder,["matlab.nc","whole.nc","split.nc","hybrid.nc"]);
                    for path = paths, copyfile(source,path); end
                    reports = cell(1,3);
                    reports{1} = testCase.runSegment(paths(2),provider,definition,617,"whole");
                    reports{2} = testCase.runSegment(paths(3),provider,definition,317,"first");
                    copyfile(paths(3),paths(4));
                    testCase.runMatlab(paths(4),definition,reports{2}.integrationRequest.selectedStep);
                    reports{3} = testCase.runSegment(paths(3),provider,definition,617,"second");
                    testCase.runMatlab(paths(1),definition,reports{1}.integrationRequest.selectedStep);
                    whole = testCase.compareModels(paths(2),paths(1));
                    split = testCase.compareModels(paths(3),paths(1));
                    hybrid = testCase.compareModels(paths(4),paths(1));
                    % Segment boundaries can change adaptive steps. The bound
                    % includes MATLAB's default 1e-3 relative error control.
                    bound = 2e-7;
                    if definition.stepPolicy == "default", bound = 5e-4; end
                    testCase.verifyLessThanOrEqual(max([whole.coefficients whole.fields whole.energy whole.enstrophy whole.tracer whole.denseFields whole.mooring]),bound,definition.id+" "+provider);
                    testCase.verifyLessThanOrEqual(max([split.coefficients split.fields split.energy split.enstrophy split.tracer split.denseFields split.mooring]),bound,definition.id+" split "+provider);
                    testCase.verifyLessThanOrEqual(max([hybrid.coefficients hybrid.fields hybrid.energy hybrid.enstrophy hybrid.tracer hybrid.denseFields hybrid.mooring]),bound,definition.id+" MATLAB restart "+provider);
                    testCase.verifyLessThanOrEqual(max([whole.particlePosition split.particlePosition hybrid.particlePosition]),.02);
                    testCase.verifyEqual(whole.graph,hybrid.graph);
                    testCase.verifyEqual(whole.graph,split.graph);
                    actual = WVModel.modelFromFile(char(paths(2))); cleanup = onCleanup(@()actual.closeNetCDFFile());
                    change = relativeError([actual.wvt.Ap;actual.wvt.Am;actual.wvt.A0],initialState);
                    if ~definition.linear, testCase.verifyGreaterThan(change,1e-4,"The qualification trajectory must exercise nontrivial evolution."); end

                    clear cleanup
                    if definition.regime=="hydrostatic-limit", testCase.verifyLessThan(verticalKineticFraction,1e-3); else, testCase.verifyGreaterThan(verticalKineticFraction,.01); end
                    row = struct(definition=definition,provider=provider,whole=whole,segmented=split,matlabRestart=hybrid,verticalKineticFraction=verticalKineticFraction,coefficientChange=change,runtime={reports});
                    writeJSON(fullfile(testCase.evidenceFolder,definition.id+"-"+provider+".json"),row);
                end
            end
        end
        function lifecycleAndStorageRemainBounded(testCase)
            executable = fullfile(fileparts(testCase.runner),"WVBoussinesqLifecycleProbe");
            testCase.assertTrue(isfile(executable));
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            manifest = jsondecode(fileread(fullfile(root,"PortableRuntime","contracts","boussinesq-qualification-cases-v1.json")));
            for definition = reshape(manifest.lifecycleCases,1,[])
                if definition.nativeOnly && ~ismember("native-fftw",testCase.providers), continue; end
                definition.grid = reshape(definition.grid,1,[]);
                if isequal(testCase.providers,"reference") && ~ismember(definition.grid,manifest.referenceOnlyLifecycleGrids,"rows"), continue; end
                source = testCase.authorModel(definition);
                for provider = testCase.providers
                    if definition.nativeOnly && provider=="reference", continue; end
                    outputPath = fullfile(testCase.evidenceFolder,definition.id+"-"+provider+".json");
                    [status,output] = cleanSystem(shellQuote(executable)+" "+shellQuote(source)+" "+shellQuote(outputPath)+" "+provider);
                    testCase.assertEqual(status,0,output);
                    report = jsondecode(fileread(outputPath));
                    report.definition = definition;
                    writeJSON(outputPath,report);
                    testCase.verifyEqual(report.completedLifecycles,6);
                    testCase.verifyTrue(report.scientificOwnersReleased);
                    testCase.verifyEqual(report.retainedGrowthBytes,0);
                    if provider=="native-fftw", testCase.verifyEqual(report.preparedStepAllocations,0); end
                end
            end
        end
    end
    methods (Access=private)
        function [path,initialState,verticalKineticFraction] = authorModel(testCase,definition)
            z = -1000*linspace(1,0,definition.grid(3))'.^1.3;
            if definition.profile == "exponential"
                N2Function = @(z) 1e-4*exp(z/700);
            else
                N2Function = @(z) 1e-4*(1+.35*tanh((z+450)/160));
            end
            wvt = WVTransformBoussinesq(reshape(definition.Lxyz,1,[]),definition.grid,Nj=definition.Nj,z=z,N2Function=N2Function,shouldAntialias=definition.shouldAntialias);
            index = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
            wvt.A0 = 1e-6*(sin(.7*index)+1i*cos(.3*index))./(1+wvt.J);
            inertial = (wvt.K.^2+wvt.L.^2)==0;
            wave = wvt.J>0 & ~inertial;
            wvt.A0(inertial) = 0;
            wvt.Ap = .001*(sin(.7*index)+1i*cos(.3*index))./(1+wvt.J).*(wave|inertial);
            wvt.Am = .002*(sin(.3*index)+1i*cos(.7*index))./(1+wvt.J).*wave;
            if definition.regime=="nonhydrostatic"
                % Preserve all components while making propagating waves
                % substantial relative to inertial shear and balanced flow.
                wvt.A0(~inertial) = .1*wvt.A0(~inertial);
                wvt.Ap(inertial) = .01*wvt.Ap(inertial);
            end
            wvt.Am(inertial) = conj(wvt.Ap(inertial));
            mda = inertial & wvt.J>0;
            wvt.A0(mda) = .003*sin(index(mda));
            wvt.t0 = -3;
            wvt.t = 17;
            forces = [WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt),WVVerticalDiffusivity(wvt,kappa_z=.002),WVBottomFrictionLinear(wvt,r=1e-5),WVBottomFrictionQuadratic(wvt,Cd=.003),WVBetaPlanePVAdvection(wvt)];
            if ~definition.shouldAntialias && definition.Nj>1, forces(end+1)=WVAntialiasing(wvt,Nj=definition.Nj-1); end
            wvt.setForcing(forces);
            fixed = WVFixedAmplitudeForcing(wvt,name="qualified fixed",A0_indices=uint64(wvt.Nj+1),A0bar=wvt.A0(wvt.Nj+1));
            wvt.addForcing(fixed);
            % Start from constrained amplitudes so the evolution gate cannot
            % pass solely because antialiasing removes initial coefficients.
            wvt.restoreForcingAmplitudes();
            initialState = [wvt.Ap;wvt.Am;wvt.A0];
            verticalKineticFraction = sum(wvt.w.^2,"all")/sum(wvt.u.^2+wvt.v.^2+wvt.w.^2,"all");
            model = WVModel(wvt,shouldUseLinearDynamics=definition.linear);
            model.eulerianObservingSystem.addNetCDFOutputVariables('u','v','w','eta','p','rho_e','qgpv','zeta_x','zeta_y','zeta_z','ssh');
            model.addParticles('drifter',definition.xyOnly,wvt.Lx*[.999 .3],wvt.Ly*[.999 .1],[-600 -100],'u','qgpv',advectionInterpolation="spline",trackedVarInterpolation="linear");
            model.addFluxedObservingSystem(WVTracer(model,name="dye",phi=.3+.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly).*exp(wvt.Z/1000),isXYOnly=definition.xyOnly));
            path = fullfile(testCase.folder,"initial.nc");
            output = model.createNetCDFFileForModelOutput(path,outputInterval=100,shouldOverwriteExisting=true);
            output.outputGroups(1).addObservingSystem(WVMooring(model,name="mooring",x=wvt.Lx*[0 .25],y=wvt.Ly*[0 .3],trackedFieldNames={'u','qgpv'}));
            dense = output.addNewEvenlySpacedOutputGroup("dense",outputInterval=25,initialTime=17,finalTime=617);
            dense.addObservingSystem(WVEulerianFields(model,fieldNames={'u','eta'}));
            output.outputTimesForIntegrationPeriod(17,617); output.writeTimeStepToOutputFile(17);
            model.closeNetCDFFile();
        end
        function runMatlab(~,path,definition,selectedStep)
            matlabModel = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()matlabModel.closeNetCDFFile());
            if definition.method == "fixed-rk4"
                matlabModel.setupIntegrator(integratorType="fixed",deltaT=selectedStep);
            elseif definition.stepPolicy == "explicit"
                matlabModel.setupIntegrator(integratorType="adaptive",integrator=str2func("ode"+extractAfter(definition.method,"adaptive-rk")),absTolerance=1e-10,relTolerance=1e-9);
                matlabModel.odeOptions = odeset(matlabModel.odeOptions,'InitialStep',5,'MaxStep',10);
            else
                matlabModel.setupIntegrator(integratorType="adaptive",integrator=str2func("ode"+extractAfter(definition.method,"adaptive-rk")));
            end
            matlabModel.integrateToTime(617,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            clear cleanup
        end
        function report = runSegment(testCase,path,provider,definition,finalTime,label)
            requestPath = fullfile(testCase.folder,"request.json");
            reportPath = fullfile(testCase.folder,"report-"+label+".json");
            common = {"method",definition.method,"finalTime",finalTime,"fftProvider",provider,"reportPath",reportPath};
            if definition.stepPolicy == "cfl"
                model = WVModel.modelFromFile(char(path)); cleanup = onCleanup(@()model.closeNetCDFFile());
                cfl = .5*5/model.timeStepForCFL(.5);
                clear cleanup
                WVModel.writePortableRunRequest(requestPath,path,common{:},cfl=cfl,timeStepConstraint="min");
            elseif definition.stepPolicy == "default"
                WVModel.writePortableRunRequest(requestPath,path,common{:});
            elseif definition.method == "fixed-rk4"
                WVModel.writePortableRunRequest(requestPath,path,common{:},initialStep=5);
            else
                WVModel.writePortableRunRequest(requestPath,path,common{:},initialStep=5,maximumStep=10,relativeTolerance=1e-9,absoluteToleranceScale=1e-10);
            end
            scientificNames = ["K2unique","iK2unique","Ppm","Qpm","h_pm","PFpmInv","QGpmInv","PFpm","QGpm","QGwg","N2Function"];
            scientific = arrayfun(@(name)ncread(path,name),scientificNames,UniformOutput=false);
            [status,output] = cleanSystem(shellQuote(testCase.runner)+" --request "+shellQuote(requestPath));
            testCase.assertEqual(status,0,output);
            report = jsondecode(fileread(reportPath));
            preserved = true;
            for index=1:numel(scientificNames)
                values=ncread(path,scientificNames(index));
                preserved=preserved && isequal(values,scientific{index});
                testCase.verifyEqual(values,scientific{index});
            end
            report.qualificationScientificStatePreserved = preserved;
            testCase.verifyEqual(string(report.status),"complete");
            testCase.verifyTrue(report.execution.noFallback);
            testCase.verifyTrue(report.integratorStorageLedger.byteLedgerAgreement);
            testCase.verifyEqual(report.storageBytes.persistentFullHermitian,0);
        end
        function errors = compareModels(testCase,actualPath,expectedPath)
            actual = WVModel.modelFromFile(char(actualPath)); actualCleanup = onCleanup(@()actual.closeNetCDFFile());
            expected = WVModel.modelFromFile(char(expectedPath)); expectedCleanup = onCleanup(@()expected.closeNetCDFFile());
            errors = struct(coefficients=relativeError([actual.wvt.Ap;actual.wvt.Am;actual.wvt.A0],[expected.wvt.Ap;expected.wvt.Am;expected.wvt.A0]),fields=0,energy=relativeError(actual.wvt.totalEnergy,expected.wvt.totalEnergy),enstrophy=relativeError(actual.wvt.totalEnstrophy(),expected.wvt.totalEnstrophy()),tracer=relativeError(actual.tracer("dye"),expected.tracer("dye")),particlePosition=0,denseFields=0,mooring=0,graph=struct());
            errors.coefficientFamilies = struct();
            for name = ["Ap","Am","A0"]
                errors.coefficientFamilies.(name) = relativeError(actual.wvt.(name),expected.wvt.(name));
            end
            errors.coefficients = max([errors.coefficientFamilies.Ap errors.coefficientFamilies.Am errors.coefficientFamilies.A0]);
            for field = ["u","v","w","eta","p","rho_e","qgpv","zeta_x","zeta_y","zeta_z","ssh"]
                errors.fields = max(errors.fields,relativeError(actual.wvt.(field),expected.wvt.(field)));
            end
            [ax,ay,az,at] = actual.drifterPositions(); [ex,ey,ez,et] = expected.drifterPositions();
            errors.particlePosition = max([abs(ax-ex),abs(ay-ey),abs(az-ez)],[],"all");
            testCase.verifyEqual(at.u,et.u,AbsTol=1e-6); testCase.verifyEqual(at.qgpv,et.qgpv,AbsTol=1e-8);
            testCase.verifyEqual(string(actual.wvt.forcingNames),string(expected.wvt.forcingNames));
            testCase.verifyEqual(actual.isDynamicsLinear,expected.isDynamicsLinear);
            for groupName = ["wave-vortex","dense"]
                a = actual.outputFiles(1).outputGroupWithName(groupName); e = expected.outputFiles(1).outputGroupWithName(groupName);
                testCase.verifyEqual(string({a.observingSystems.name}),string({e.observingSystems.name}));
                testCase.verifyEqual(string(arrayfun(@class,a.observingSystems,UniformOutput=false)),string(arrayfun(@class,e.observingSystems,UniformOutput=false)));
                for property = ["outputInterval","initialTime","finalTime","incrementsWrittenToGroup"]
                    testCase.verifyEqual(a.(property),e.(property));
                end
                times = ncread(actualPath,"/"+groupName+"/t");
                testCase.verifyEqual(times,ncread(expectedPath,"/"+groupName+"/t"));
                actualSchema = ncinfo(actualPath,"/"+groupName);
                expectedSchema = ncinfo(expectedPath,"/"+groupName);
                testCase.verifyEqual(actualSchema.Variables,expectedSchema.Variables);
                errors.graph.(replace(groupName,"-","_")) = struct(times=times,ordinal=a.incrementsWrittenToGroup,observerNames=string({a.observingSystems.name}),variableNames=string({actualSchema.Variables.Name}));
            end
            clear actualCleanup expectedCleanup
            for field = ["mooring_u","mooring_qgpv"]
                errors.mooring = max(errors.mooring,relativeError(ncread(actualPath,"/wave-vortex/"+field),ncread(expectedPath,"/wave-vortex/"+field)));
            end
            for field = ["u","eta"]
                errors.denseFields = max(errors.denseFields,relativeError(ncread(actualPath,"/dense/"+field),ncread(expectedPath,"/dense/"+field)));
            end
        end
    end
end
function error = relativeError(actual,expected)
error = max(abs(actual-expected),[],"all")/max(1e-20,max(abs(expected),[],"all"));
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
function [status,output] = cleanSystem(command)
[status,output] = system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
function writeJSON(path,value)
file = fopen(path,"w"); assert(file>=0,"Unable to write qualification evidence."); cleanup = onCleanup(@()fclose(file)); fprintf(file,"%s\n",jsonencode(value,PrettyPrint=true));
end
