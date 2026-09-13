classdef TestPortableNoMotionProfile < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        executable (1,1) string
    end
    methods (TestClassSetup)
        function prepareProbe(testCase)
            testCase.root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_NO_MOTION_PROFILE_DUMP"));
            if testCase.executable=="" && getenv("WV_DIAGNOSTIC_DUMP")~=""
                testCase.executable = fullfile(fileparts(string(getenv("WV_DIAGNOSTIC_DUMP"))),"WVNoMotionProfileDump");
            end
            if testCase.executable==""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target WVNoMotionProfileDump wave-vortex-run");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVNoMotionProfileDump");
            end
            testCase.assertTrue(isfile(testCase.executable));
        end
    end
    methods (Test,TestTags="full")
        function linearAndTinyDisplacementsMatchIndependentTruth(testCase)
            N2 = .04;
            for knots = {[-1;1],[-1;-.8;-.15;0;.7;1]}
                z = [-1;1;0;.7;-.8;1e-12;-1e-12;.5;.5+eps(.5);-.5-eps(.5)];
                material = [1;-1;0;-.8;.7;0;0;.5;.5;-.5];
                request = struct(z=knots{1},rho=1-N2*knots{1},queryHeights=z,materialHeights=material, ...
                    targetDensity=1-N2*material,g=1,rho0=1);
                actual = testCase.compareWithMatlab(request,"linear-"+numel(knots{1}));
                [~,expected] = DensityDiagnosticReference.linearReference(z,material,N2);
                testCase.verifyEqual(actual.density,1-N2*z,AbsTol=2e-15);
                testCase.verifyEqual(actual.inverse,material,AbsTol=1e-13);
                testCase.verifyEqual(actual.ape,expected,RelTol=2e-13);
                testCase.verifyGreaterThan(actual.ape(z~=material),zeros(nnz(z~=material),1));
            end
        end
        function nonlinearProfilesConvergeToIndependentTruth(testCase)
            references = [DensityDiagnosticReference.nonlinearRestProfile(N2=.04,curvature=.3), ...
                DensityDiagnosticReference.nonlinearRestProfile(N2=.07,curvature=-.4)];
            [x,z] = ndgrid(linspace(-.8,.8,19));
            [~,material] = DensityDiagnosticReference.inversePolarTwist(x,z,[.13 -.17],.7,1.3);
            for reference = references
                errors = zeros(3,2);
                counts = [17 65 257];
                for level = 1:numel(counts)
                    knots = linspace(-1,1,counts(level)).';
                    request = struct(z=knots,rho=reference.density(knots),queryHeights=z(:),materialHeights=material(:), ...
                        targetDensity=reference.density(material(:)),g=reference.gravity,rho0=reference.rho0);
                    actual = testCase.compareWithMatlab(request,"nonlinear-"+reference.N2+"-"+counts(level));
                    errors(level,1) = max(abs(actual.inverse-material(:)));
                    errors(level,2) = max(abs(actual.ape-reference.availablePotentialEnergy(z(:),material(:))));
                end
                testCase.verifyLessThan(errors(end,:),errors(1,:)/100);
                testCase.verifyLessThan(errors(end,:),[1e-8 1e-9]);
                fprintf("NO_MOTION_NONLINEAR_CONVERGENCE N2=%.3g inverse=%.17g ape=%.17g\n",reference.N2,errors(end,:));
            end
        end
        function irregularAndWeakProfilesPreserveCalculus(testCase)
            knots = [-2;-1.7;-.8;-.79;-.1;.25;1];
            densities = [9;8.99;8.7;6;5.9;3;1];
            material = [-2;-1.7;-.8;-.795;-.79;-.1;.25;1;-.25;.8;-.1;-.1];
            z = [1;-.1;.25;-.794;-2;.25;-.8;-2;.9;-.7;-.1;-.8];
            pp = pchip(knots,densities);
            request = struct(z=knots,rho=densities,queryHeights=z,materialHeights=material, ...
                targetDensity=ppval(pp,material),g=9.81,rho0=1025);
            actual = testCase.compareWithMatlab(request,"irregular");
            derivative = mkpp(pp.breaks,pp.coefs(:,1:end-1).*(pp.order-1:-1:1));
            expected = zeros(size(z));
            for index = 1:numel(z)
                lower = min(z(index),material(index));
                upper = max(z(index),material(index));
                if lower==upper, continue; end
                integrand = @(r) (r-z(index)).*ppval(derivative,r);
                expected(index) = (request.g/request.rho0)*sign(z(index)-material(index)) * ...
                    integral(integrand,lower,upper,Waypoints=knots(knots>lower & knots<upper),AbsTol=1e-12,RelTol=1e-12);
            end
            testCase.verifyEqual(actual.ape,expected,AbsTol=5e-13,RelTol=5e-12);
            rho0 = 1025; g = 9.81; N2 = 1e-10;
            knots = [-1000;-731;-207;0];
            density = @(z) rho0*(1-N2*z/g);
            material = [-850;-500;-100];
            z = material+[25;-10;2];
            request = struct(z=knots,rho=density(knots),queryHeights=z,materialHeights=material, ...
                targetDensity=density(material),g=g,rho0=rho0);
            actual = testCase.compareWithMatlab(request,"weak-gradient");
            testCase.verifyEqual(actual.inverse,material,AbsTol=4*eps(rho0)*g/(rho0*N2));
            testCase.verifyEqual(actual.ape,N2*(z-material).^2/2,RelTol=2e-6);
        end
        function capturedGivenProfilesMatchMatlab(testCase)
            for day = ["3000","3250"]
                path = fullfile(testCase.root,"UnitTests","fixtures","portable-no-motion-profile-day"+day+".json");
                fixture = jsondecode(fileread(path));
                testCase.assertEqual(string(fixture.schema),"wave-vortex-supplied-no-motion-profile-fixture-v1");
                for source = reshape(fixture.sourceSHA256,1,[])
                    testCase.assertEqual(sha256(fullfile(testCase.root,source.path)),string(source.sha256));
                end
                testCase.assertEqual(sha256(fullfile(testCase.root,fixture.generator.path)),string(fixture.generator.sha256));
                actual = testCase.compareWithMatlab(fixture.request,"captured-day"+day);
                testCase.verifyEqual(actual.density,fixture.expected.density,AbsTol=8*eps(max(abs(fixture.request.rho))));
                testCase.verifyEqual(actual.inverse,fixture.expected.inverse,AbsTol=64*eps(max(abs(fixture.request.z))+(max(fixture.request.z)-min(fixture.request.z))));
                testCase.verifyEqual(actual.ape,fixture.expected.ape,RelTol=2e-11);
            end
        end
        function runtimeReportsDensityReferenceWithoutChangingIntegration(testCase)
            runner = fullfile(fileparts(testCase.executable),"wave-vortex-run");
            testCase.assertTrue(isfile(runner));
            wvt = WVTransformConstantStratification([17000 11000 1000],[8 6 9],N0=5.2e-3);
            wvt.removeAllForcing();
            index = reshape(1:numel(wvt.A0),size(wvt.A0));
            wave = wvt.J>0 & wvt.Kh>0;
            wvt.Ap = .001*sin(.7*index).*wave;
            wvt.Am = .002*cos(.3*index).*wave;
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            source = fullfile(testCase.folder,"contract-source.nc");
            file = model.createNetCDFFileForModelOutput(source,outputInterval=1,shouldOverwriteExisting=true);
            file.outputTimesForIntegrationPeriod(0,1);
            file.writeTimeStepToOutputFile(0);
            model.closeNetCDFFile();
            baseline = struct;
            rows = struct(selection={},contract={},coefficientMaximumDifference={});
            for selection = ["absent","empty","initial"]
                destination = fullfile(testCase.folder,"contract-"+selection+".nc");
                copyfile(source,destination);
                requestPath = fullfile(testCase.folder,"contract-request.json");
                reportPath = fullfile(testCase.folder,"contract-report.json");
                WVModel.writePortableRunRequest(requestPath,destination,method="fixed-rk4",finalTime=1,initialStep=1,fftProvider="reference",reportPath=reportPath);
                request = jsondecode(fileread(requestPath));
                if selection=="empty"
                    request.execution.densityDiagnostics = struct;
                elseif selection=="initial"
                    request.execution.densityDiagnostics = struct(contract="wave-vortex-density-diagnostics-v1",reference="initial");
                end
                writeJSON(requestPath,request);
                [status,output] = cleanSystem(shellQuote(runner)+" --request "+shellQuote(requestPath));
                testCase.assertEqual(status,0,output);
                report = jsondecode(fileread(reportPath));
                testCase.assertEqual(string(report.status),"complete");
                testCase.assertTrue(report.integrationRequest.noFallback);
                testCase.assertEqual(report.state.stepCount,1);
                contract = report.densityDiagnosticContract;
                testCase.verifyEqual(string(contract.identifier),"wave-vortex-density-diagnostics-v1");
                testCase.verifyEqual(string(contract.outputEvaluation),"available");
                if selection=="initial"
                    testCase.verifyEqual(string(contract.reference),"initial");
                    testCase.verifyEqual(string(contract.profileRecovery),"not-required");
                else
                    testCase.verifyEqual(string(contract.reference),"actual");
                    testCase.verifyEqual(string(contract.profileRecovery),"dampedLeastSquares");
                end
                expectedSource = "request";
                if selection=="absent", expectedSource="default"; end
                testCase.verifyEqual(string(contract.selectionSource),expectedSource);
                maximumDifference = 0;
                info = ncinfo(destination,"/wave-vortex");
                variableNames = string({info.Variables.Name});
                for family = ["Ap","Am","A0"]
                    if ismember(family,variableNames)
                        values = ncread(destination,"/wave-vortex/"+family);
                    else
                        values = complex(ncread(destination,"/wave-vortex/"+family+"_real"),ncread(destination,"/wave-vortex/"+family+"_imag"));
                    end
                    if selection=="absent", baseline.(family)=values; end
                    testCase.verifyEqual(values,baseline.(family));
                    maximumDifference = max(maximumDifference,max(abs(values-baseline.(family)),[],"all"));
                end
                rows(end+1) = struct(selection=selection,contract=contract,coefficientMaximumDifference=maximumDifference); %#ok<AGROW>
            end
            directory = string(getenv("WV_NO_MOTION_PROFILE_REPORT_DIR"));
            if directory~=""
                if ~isfolder(directory), mkdir(directory); end
                writeJSON(fullfile(directory,"runtime-contract.json"),struct(schema="wave-vortex-density-reference-request-parity-v1",rows=rows));
            end
        end
    end
    methods (Access=private)
        function actual = compareWithMatlab(testCase,request,label)
            profile = WVNoMotionProfile(request.z,request.rho);
            expected = struct(density=profile.density(request.queryHeights),inverse=profile.inverse(request.targetDensity), ...
                ape=profile.availablePotentialEnergy(request.queryHeights,request.materialHeights,request.g,request.rho0));
            requestPath = fullfile(testCase.folder,"request.json");
            resultPath = fullfile(testCase.folder,"result.json");
            writeJSON(requestPath,request);
            [status,output] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(requestPath)+" "+shellQuote(resultPath));
            testCase.assertEqual(status,0,label+": "+output);
            actual = jsondecode(fileread(resultPath));
            testCase.assertEqual(string(actual.schema),"wv-no-motion-profile-dump-v1");
            testCase.assertEqual(string(actual.status),"complete");
            testCase.assertEqual(actual.knotCount,numel(request.z));
            for name = ["density","inverse","ape"]
                actual.(name) = actual.(name)(:);
                testCase.assertSize(actual.(name),size(expected.(name)));
                testCase.assertTrue(all(isfinite(actual.(name))),label+" "+name);
            end
            densityTolerance = 8*eps(max(abs(request.rho)));
            heightTolerance = 64*eps(max(abs(request.z))+(max(request.z)-min(request.z)));
            testCase.verifyEqual(actual.density,expected.density,AbsTol=densityTolerance);
            testCase.verifyEqual(actual.inverse,expected.inverse,AbsTol=heightTolerance);
            testCase.verifyEqual(actual.ape,expected.ape,RelTol=2e-11);
            testCase.verifyGreaterThanOrEqual(actual.ape,zeros(size(actual.ape)));
            testCase.verifyEqual(actual.ape(request.queryHeights==request.materialHeights),zeros(nnz(request.queryHeights==request.materialHeights),1));
            closure = max(abs(profile.density(actual.inverse)-request.targetDensity));
            testCase.verifyLessThanOrEqual(closure,densityTolerance);
            error = struct(density=max(abs(actual.density-expected.density)),inverse=max(abs(actual.inverse-expected.inverse)), ...
                ape=max(abs(actual.ape-expected.ape)),apeRelative=max(abs(actual.ape-expected.ape)./max(abs(expected.ape),realmin)),densityClosure=closure);
            fprintf("NO_MOTION_PROFILE_PARITY %s density=%.17g inverse=%.17g apeRelative=%.17g closure=%.17g\n",label,error.density,error.inverse,error.apeRelative,closure);
            directory = string(getenv("WV_NO_MOTION_PROFILE_REPORT_DIR"));
            if directory~=""
                if ~isfolder(directory), mkdir(directory); end
                report = struct(schema="wave-vortex-supplied-no-motion-profile-parity-v1",caseName=label, ...
                    knotCount=numel(request.z),queryCount=numel(request.queryHeights),matlabRelease=string(version('-release')), ...
                    maximumError=error,tolerance=struct(density=densityTolerance,inverse=heightTolerance,apeRelative=2e-11));
                writeJSON(fullfile(directory,label+".json"),report);
            end
        end
    end
end

function writeJSON(path,value)
stream = fopen(path,'w');
assert(stream>=0,"Cannot open the supplied-profile test artifact.");
cleanup = onCleanup(@()fclose(stream));
fprintf(stream,"%s\n",jsonencode(value,PrettyPrint=true));
end
function result = sha256(path)
[status,digest] = cleanSystem("shasum -a 256 "+shellQuote(path));
assert(status==0,"Cannot hash the supplied-profile fixture source.");
result = string(extractBefore(digest," "));
end
function [status,output] = cleanSystem(command)
[status,output] = system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
