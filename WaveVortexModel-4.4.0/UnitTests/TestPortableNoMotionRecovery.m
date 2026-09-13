classdef TestPortableNoMotionRecovery < matlab.unittest.TestCase
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
            testCase.executable = string(getenv("WV_NO_MOTION_RECOVERY_DUMP"));
            if testCase.executable=="" && getenv("WV_DIAGNOSTIC_DUMP")~=""
                testCase.executable = fullfile(fileparts(string(getenv("WV_DIAGNOSTIC_DUMP"))),"WVNoMotionRecoveryDump");
            end
            if testCase.executable==""
                build = fullfile(testCase.folder,"build");
                [status,output] = cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output] = cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target WVNoMotionRecoveryDump");
                testCase.assertEqual(status,0,output);
                testCase.executable = fullfile(build,"WVNoMotionRecoveryDump");
            end
            testCase.assertTrue(isfile(testCase.executable));
        end
    end
    methods (Test,TestTags="full")
        function weightedMomentsMatchIndependentParcelVolumes(testCase)
            density = reshape([1000 1004 1002 1002 1004 1000],2,1,3);
            request = struct(mode="moments",zInt=[.2;.3;.5],Lz=1,shape=[2 1 3],rhoTotal=density(:));
            actual = testCase.runProbe(request,"weighted-moments");
            testCase.assertEqual(string(actual.status),"complete");
            testCase.verifyEqual(actual.distribution.minimumDensity,1000);
            testCase.verifyEqual(actual.distribution.maximumDensity,1004);
            testCase.verifyEqual(actual.distribution.sampleCount,numel(density));
            expected = .35+.3*(.5.^(1:3)).';
            testCase.verifyEqual(actual.distribution.moments,expected,AbsTol=3e-16);
            matlab = WVNoMotionProfileOperation.moments_from_rho_tot(density,1000,1004,request.zInt,1);
            testCase.verifyEqual(actual.distribution.moments,matlab,AbsTol=3e-16);
        end
        function changedStableRestIsRecoveredExactly(testCase)
            z = linspace(-1,0,17).';
            initial = 1025-z;
            changed = 1025-2*z-.3*z.^2;
            weights = ones(size(z))/16; weights([1 end])=weights([1 end])/2;
            density = repmat(reshape(changed,1,1,[]),4,3,1);
            request = recoveryRequest(weights,density,initial);
            actual = testCase.runProbe(request,"stable-rest");
            testCase.verifyQualified(actual);
            testCase.verifyEqual(actual.profile,changed);
            testCase.verifyEqual(actual.report.iterations,0);
            testCase.verifyEqual(actual.report.maximumResidual,0);
            testCase.verifyEqual(string(actual.report.algorithm),"stable-rest-profile");
        end
        function equalVolumeRearrangementsRecoverIndependentProfiles(testCase)
            for curvature = [0 .3]
                z = linspace(-1,0,5).';
                initial = 1025-z;
                changed = 1026-2*z-curvature*z.^2;
                weights = [.125;.25;.25;.25;.125];
                density = repmat(reshape(changed,1,1,[]),4,3,1);
                density(1,1,[2 4]) = density(1,1,[4 2]);
                request = recoveryRequest(weights,density,initial);
                actual = testCase.runProbe(request,"rearranged-"+curvature);
                testCase.verifyQualified(actual);
                [matlab,exitFlag,output] = WVNoMotionProfileOperation.find_rho_nm(weights,1,density,initial);
                testCase.assertGreaterThan(exitFlag,0);
                testCase.assertLessThanOrEqual(output.maximumResidual,1e-8);
                testCase.verifyEqual(actual.profile,changed,AbsTol=2e-7);
                testCase.verifyEqual(actual.profile,matlab,AbsTol=2e-9);
                testCase.verifyEqual(actual.profile([1 end]),[max(density,[],"all");min(density,[],"all")]);
                expectedHeight = repmat(reshape(z,1,1,[]),4,3,1);
                expectedHeight(1,1,[2 4]) = expectedHeight(1,1,[4 2]);
                profile = WVNoMotionProfile(z,actual.profile);
                testCase.verifyEqual(profile.inverse(density),expectedHeight,AbsTol=2e-7);
                if curvature==0
                    currentHeight = repmat(reshape(z,1,1,[]),4,3,1);
                    expectedAPE = (9.81/1025)*(currentHeight-expectedHeight).^2;
                    testCase.verifyEqual(profile.availablePotentialEnergy(currentHeight,expectedHeight,9.81,1025),expectedAPE,AbsTol=1e-12);
                end
            end
        end
        function capturedMomentsQualifyAndPreservePhysicalProfile(testCase)
            for day = ["3000","3250"]
                fixturePath = fullfile(testCase.root,"UnitTests","fixtures","density-run18-day"+day+".json");
                fixture = jsondecode(fileread(fixturePath));
                minimum = fixture.rho_nm0(end); maximum = fixture.rho_nm0(1);
                request = struct(mode="fitMoments",zInt=fixture.z_int,Lz=fixture.Lz,rhoNM0=fixture.rho_nm0, ...
                    minimumDensity=minimum,maximumDensity=maximum,targetMoments=fixture.moments);
                actual = testCase.runProbe(request,"captured-"+day);
                testCase.verifyQualified(actual);
                testCase.verifyLessThan(actual.report.maximumResidual,1e-9);
                u0 = WVNoMotionProfileOperation.rho_to_u(fixture.rho_nm0,minimum,maximum);
                [u,exitFlag,output] = WVNoMotionProfileOperation.solveMoments(u0,flip(fixture.z_int/fixture.Lz),fixture.moments);
                testCase.assertGreaterThan(exitFlag,0);
                matlab = WVNoMotionProfileOperation.u_to_rho(u,minimum,maximum);
                profileDifference = max(abs(actual.profile-matlab));
                % Raw-moment Jacobians are numerically rank deficient here.
                % Changing floating-point contraction alone changes fitted
                % nodes while preserving residual accuracy. Compare physical
                % profiles at one ppm of density span and domain depth;
                % retain the independent moment and synthetic-truth gates.
                profileTolerance = 1e-6*(maximum-minimum);
                heightTolerance = 1e-6*fixture.Lz;
                testCase.verifyLessThanOrEqual(profileDifference,profileTolerance);
                cppProfile = WVNoMotionProfile(fixture.z,actual.profile);
                matlabProfile = WVNoMotionProfile(fixture.z,matlab);
                query = linspace(fixture.z(1),fixture.z(end),257).';
                density = matlabProfile.density(query);
                heightDifference = max(abs(cppProfile.inverse(density)-matlabProfile.inverse(density)));
                testCase.verifyLessThanOrEqual(heightDifference,heightTolerance);
                directMoments = sum((fixture.z_int/fixture.Lz).*((actual.profile-minimum)/(maximum-minimum)).^(1:numel(matlab)),1).';
                testCase.verifyLessThanOrEqual(max(abs(directMoments-fixture.moments)),1e-8);
                metrics = struct(profileMaximumDifference=profileDifference,inverseHeightMaximumDifference=heightDifference, ...
                    cppResidual=actual.report.maximumResidual,matlabResidual=output.maximumResidual, ...
                    directMomentResidual=max(abs(directMoments-fixture.moments)), ...
                    profileTolerance=profileTolerance,inverseHeightTolerance=heightTolerance);
                testCase.writeEvidence("captured-"+day+"-physical",metrics);
                fprintf("RECOVERY_CAPTURED day=%s rhoError=%.17g heightError=%.17g residual=%.17g\n",day,profileDifference,heightDifference,actual.report.maximumResidual);
            end
        end
        function unrepresentableDistributionsNeverQualify(testCase)
            weights = [.125;.25;.25;.25;.125];
            initial = linspace(1027,1025,5).';
            density = repmat(reshape(initial,1,1,[]),8,8,1);
            density(1,1,3) = initial(1)+10;
            actual = testCase.runProbe(recoveryRequest(weights,density,initial),"outlier");
            testCase.verifyEqual(string(actual.status),"failed");
            testCase.verifyFalse(actual.report.qualified);
            testCase.verifyGreaterThan(actual.report.maximumResidual,1e-8);
            actual = testCase.runProbe(recoveryRequest(weights,ones(2,2,5),initial),"constant-density");
            testCase.verifyEqual(string(actual.status),"failed");
        end
        function exhaustedBudgetsNeverQualify(testCase)
            request = struct(mode="fitMoments",zInt=[.25;.5;.25],Lz=1,rhoNM0=[3;2;1], ...
                minimumDensity=1,maximumDensity=3,targetMoments=[.55;.4;.33]);
            for budget = ["maximumIterations","maximumEvaluations"]
                request.options = struct;
                request.options.(budget) = double(budget=="maximumEvaluations");
                actual = testCase.runProbe(request,"exhausted-"+budget);
                testCase.verifyEqual(string(actual.status),"failed");
                testCase.verifyFalse(actual.report.qualified);
                testCase.verifyEqual(actual.report.exitFlag,0);
                testCase.verifyEqual(actual.report.evaluations,1);
            end
        end
    end
    methods (Access=private)
        function result = runProbe(testCase,request,label)
            input = fullfile(testCase.folder,"request.json");
            output = fullfile(testCase.folder,"result.json");
            writeJSON(input,request);
            [status,message] = cleanSystem(shellQuote(testCase.executable)+" "+shellQuote(input)+" "+shellQuote(output));
            testCase.assertEqual(status,0,message);
            result = jsondecode(fileread(output));
            testCase.writeEvidence(label,result);
        end
        function verifyQualified(testCase,result)
            testCase.assertEqual(string(result.status),"complete",string(result.message));
            testCase.assertTrue(result.report.qualified);
            testCase.verifyGreaterThan(result.report.exitFlag,0);
            testCase.verifyLessThanOrEqual(result.report.maximumResidual,1e-8);
            testCase.verifyTrue(all(isfinite(result.profile)));
            testCase.verifyLessThan(diff(result.profile),zeros(numel(result.profile)-1,1));
        end
        function writeEvidence(~,label,value)
            directory = string(getenv("WV_NO_MOTION_RECOVERY_REPORT_DIR"));
            if directory=="", return; end
            if ~isfolder(directory), mkdir(directory); end
            writeJSON(fullfile(directory,label+".json"),value);
        end
    end
end

function request = recoveryRequest(weights,density,initial)
request = struct(mode="recover",zInt=weights,Lz=1,rhoNM0=initial,shape=size(density),rhoTotal=density(:));
end
function writeJSON(path,value)
stream = fopen(path,'w');
assert(stream>=0,"Cannot open the no-motion recovery artifact.");
cleanup = onCleanup(@()fclose(stream));
fprintf(stream,"%s\n",jsonencode(value,PrettyPrint=true));
end
function [status,output] = cleanSystem(command)
[status,output] = system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
