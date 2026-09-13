classdef TestPortableDensityEventEvaluation < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        executable (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function prepareProbe(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture(PreservingOnFailure=true));
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_DENSITY_EVENT_DUMP"));
            if testCase.executable=="" && getenv("WV_DIAGNOSTIC_DUMP")~=""
                testCase.executable = fullfile(fileparts(string(getenv("WV_DIAGNOSTIC_DUMP"))),"WVDensityEventDump");
            end
            if testCase.executable==""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target WVDensityEventDump");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVDensityEventDump");
            end
            testCase.assertTrue(isfile(testCase.executable));
            testCase.providers = "reference";
            if getenv("WV_DIAGNOSTIC_NATIVE")=="1", testCase.providers=["reference","native"]; end
        end
    end
    methods (Test,TestTags="full")
        function coefficientDrivenDensityMatchesMatlab(testCase)
            rows = {};
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for antialias = [false true]
                    grid = [9 7 10];
                    if antialias, grid=[8 6 9]; end
                    wvt = makeTransform(family,grid,antialias);
                    n = reshape(1:numel(wvt.A0),size(wvt.A0));
                    wave = wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0;
                    wvt.Ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*wave;
                    wvt.Am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
                    wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n));
                    mda = (wvt.K.^2+wvt.L.^2)==0 & wvt.J>0;
                    wvt.A0(mda) = .003*sin(n(mda));
                    wvt.t0=17; wvt.t=37;
                    checkpoint = testCase.writeSnapshot(wvt);
                    expectedProfile = wvt.rho_nm;
                    span = max(wvt.rho_total,[],"all")-min(wvt.rho_total,[],"all");
                    for selection = ["actual","initial"]
                        wvt.shouldUseTrueNoMotionProfile=selection=="actual";
                        expectedEta = wvt.eta_true; expectedAPE=wvt.ape; expectedAPV=wvt.apv;
                        testCase.assertGreaterThan(max(abs(expectedEta),[],"all"),1e-5);
                        for provider = testCase.providers
                            actual = testCase.runProbe(checkpoint,selection,provider);
                            label=family+" aa="+antialias+" "+selection+" "+provider;
                            profileError=max(abs(actual.fields.rho_nm.real(:)-expectedProfile(:)));
                            etaError=max(abs(actual.fields.eta_true.real(:)-expectedEta(:)));
                            apeError=max(abs(actual.fields.ape.real(:)-expectedAPE(:)));
                            apeScale=max(abs(expectedAPE),[],"all");
                            cppEta=reshape(actual.fields.eta_true.real,grid);
                            cppAPV=reshape(actual.fields.apv.real,grid);
                            fixedEtaAPV=wvt.zeta_z-wvt.zeta_x.*wvt.diffX(cppEta)-wvt.zeta_y.*wvt.diffY(cppEta)-(wvt.zeta_z+wvt.f).*wvt.diffZG(cppEta);
                            apvScale=max(abs(expectedAPV),[],"all");
                            apvError=max(abs(cppAPV-expectedAPV),[],"all");
                            apvOperatorError=max(abs(cppAPV-fixedEtaAPV),[],"all");
                            testCase.verifyLessThanOrEqual(apvError,1e-6*apvScale+1e-12,label+" apv");
                            testCase.verifyLessThanOrEqual(apvOperatorError,1e-10*apvScale+1e-13,label+" apv given displacement");
                            selectedProfile=wvt.rho_nm0;
                            if selection=="actual", selectedProfile=actual.fields.rho_nm.real(:); end
                            profile=WVNoMotionProfile(wvt.z,selectedProfile);
                            materialHeight=wvt.Z-reshape(actual.fields.eta_true.real,grid);
                            densityClosure=max(abs(profile.density(materialHeight)-wvt.rho_total),[],"all");
                            testCase.verifyLessThanOrEqual(densityClosure,16*eps(max(wvt.rho_total,[],"all")));
                            % The raw-moment fit is arithmetic-sensitive. These
                            % physical bounds match the qualified recovery contract.
                            testCase.verifyLessThanOrEqual(profileError,1e-6*span,label+" rho_nm");
                            testCase.verifyLessThanOrEqual(etaError,1e-6*wvt.Lz,label+" eta_true");
                            testCase.verifyLessThanOrEqual(apeError,1e-5*apeScale+1e-12,label+" ape");
                            testCase.verifyTrue(all(isfinite(actual.fields.eta_true.real)));
                            testCase.verifyGreaterThanOrEqual(min(actual.fields.ape.real),0);
                            testCase.verifyEqual(actual.fields.rho_nm.dimensions(:),wvt.Nz);
                            testCase.verifyEqual(actual.fields.eta_true.dimensions(:),grid(:));
                            testCase.verifyMetrics(actual);
                            rows{end+1}=struct(configuration=family,antialias=antialias,grid=grid,reference=selection,provider=provider, ...
                                profileMaximumError=profileError,profileTolerance=1e-6*span,etaMaximumError=etaError,etaTolerance=1e-6*wvt.Lz, ...
                                apeMaximumError=apeError,apeTolerance=1e-5*apeScale+1e-12,densityClosureMaximum=densityClosure, ...
                                apvMaximumError=apvError,apvTolerance=1e-6*apvScale+1e-12,apvOperatorMaximumError=apvOperatorError, ...
                                apvOperatorTolerance=1e-10*apvScale+1e-13,metrics=actual.metrics); %#ok<AGROW>
                        end
                    end
                end
            end
            reportPath=string(getenv("WV_DENSITY_EVENT_REPORT"));
            if reportPath~="", writeJSON(reportPath,struct(schema="density-event-matlab-parity-v1",rows=[rows{:}])); end
        end
        function stableRestHasIndependentZeroDisplacementAndEnergy(testCase)
            wvt=makeTransform("constant-hydrostatic",[8 6 9],true);
            checkpoint=testCase.writeSnapshot(wvt);
            for selection=["actual","initial"]
                for provider=testCase.providers
                    actual=testCase.runProbe(checkpoint,selection,provider);
                    testCase.verifyEqual(actual.fields.rho_nm.real(:),wvt.rho_nm0(:),AbsTol=4*eps(wvt.rho0));
                    testCase.verifyEqual(actual.fields.eta_true.real,zeros(numel(wvt.Z),1),AbsTol=1e-8);
                    testCase.verifyEqual(actual.fields.ape.real,zeros(numel(wvt.Z),1),AbsTol=1e-20);
                    testCase.verifyEqual(actual.fields.apv.real,zeros(numel(wvt.Z),1),AbsTol=1e-12);
                    testCase.verifyMetrics(actual);
                end
            end
        end
        function initialReferenceAPVMatchesResolvedAnalyticFlow(testCase)
            % A linear initial profile makes the imposed displacement exact.
            % This is an initial-reference calculus control, not a claim that
            % the arbitrary displaced profile preserves the density distribution.
            rows={};
            for antialias=[false true]
                wvt=WVTransformConstantStratification([2 2 2],[16 16 17],N0=.2,shouldAntialias=antialias);
                X=wvt.X; Y=wvt.Y; Z=wvt.Z; m=pi/2;
                u=.03*sin(pi*Y).*cos(m*Z); v=.02*sin(pi*X).*cos(m*Z);
                zetaX=.02*m*sin(pi*X).*sin(m*Z);
                zetaY=-.03*m*sin(pi*Y).*sin(m*Z);
                zetaZ=pi*(.02*cos(pi*X)-.03*cos(pi*Y)).*cos(m*Z);
                for axis=1:2
                    coordinate=X; horizontalVorticity=zetaX;
                    if axis==2, coordinate=Y; horizontalVorticity=zetaY; end
                    eta=.04*sin(pi*coordinate).*sin(m*(Z+2));
                    etaHorizontal=.04*pi*cos(pi*coordinate).*sin(m*(Z+2));
                    etaVertical=.04*m*sin(pi*coordinate).*cos(m*(Z+2));
                    expected=zetaZ-horizontalVorticity.*etaHorizontal-(zetaZ+wvt.f).*etaVertical;
                    wvt.initWithUVEta(u,v,eta);
                    wvt.shouldUseTrueNoMotionProfile=false;
                    testCase.assertGreaterThan(max(abs(expected-zetaZ),[],"all"),1e-4);
                    testCase.verifyEqual(wvt.apv,expected,AbsTol=1e-10);
                    checkpoint=testCase.writeSnapshot(wvt);
                    for provider=testCase.providers
                        actual=testCase.runProbe(checkpoint,"initial",provider,["eta_true","apv"]);
                        testCase.verifyEqual(reshape(actual.fields.apv.real,size(Z)),expected,AbsTol=1e-10);
                        testCase.verifyEqual(reshape(actual.fields.eta_true.real,size(Z)),eta,AbsTol=1e-10);
                        testCase.verifyEqual(actual.metrics.recoveryCount,0);
                        testCase.verifyEqual(actual.metrics.inversePassCount,1);
                        testCase.verifyEqual(actual.metrics.apvPassCount,1);
                        testCase.verifyEqual(actual.metrics.apePassCount,0);
                        testCase.verifyEqual(actual.metrics.liveBytes,0);
                        rows{end+1}=struct(axis=axis,antialias=antialias,provider=provider,reference="initial", ...
                            etaMaximumError=max(abs(actual.fields.eta_true.real(:)-eta(:))), ...
                            apvMaximumError=max(abs(actual.fields.apv.real(:)-expected(:))), ...
                            absoluteTolerance=1e-10,metrics=actual.metrics); %#ok<AGROW>
                    end
                end
            end
            reportPath=string(getenv("WV_DENSITY_EVENT_REPORT"));
            if reportPath~="", writeJSON(replace(reportPath,".json","-analytic.json"),struct(schema="density-apv-analytic-v1",rows=[rows{:}])); end
        end
    end
    methods (Access=private)
        function checkpoint=writeSnapshot(testCase,wvt)
            wvt.removeAllForcing();
            model=WVModel(wvt,shouldUseLinearDynamics=true);
            checkpoint=fullfile(testCase.folder,"source.nc");
            file=model.createNetCDFFileForModelOutput(checkpoint,outputInterval=.5,shouldOverwriteExisting=true);
            file.outputTimesForIntegrationPeriod(wvt.t,wvt.t+1);
            file.writeTimeStepToOutputFile(wvt.t);
            model.closeNetCDFFile();
        end
        function actual=runProbe(testCase,checkpoint,selection,provider,names)
            if nargin<5, names=["rho_nm","eta_true","ape","apv"]; end
            request=fullfile(testCase.folder,"request.json"); result=fullfile(testCase.folder,"result.json");
            writeJSON(request,struct(fields=names,reference=selection));
            [status,output]=cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(checkpoint)+" "+shellQuote(request)+" "+shellQuote(result)+" "+provider);
            testCase.assertEqual(status,0,selection+" "+provider+": "+output);
            actual=jsondecode(fileread(result));
        end
        function verifyMetrics(testCase,actual)
            testCase.verifyEqual(actual.metrics.recoveryCount,1);
            testCase.verifyEqual(actual.metrics.profileConstructionCount,1);
            testCase.verifyEqual(actual.metrics.inversePassCount,1);
            testCase.verifyEqual(actual.metrics.apePassCount,1);
            testCase.verifyEqual(actual.metrics.apvPassCount,1);
            testCase.verifyGreaterThan(actual.metrics.apvReuseCount,0);
            testCase.verifyGreaterThan(actual.metrics.reuseCount,0);
            testCase.verifyEqual(actual.metrics.liveBytes,0);
            testCase.verifyEqual(actual.metrics.eventLiveBytes,0);
            testCase.verifyEqual(actual.replayRecoveryCount,2);
        end
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
