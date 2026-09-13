classdef TestStratifiedQGCompiledKernel < matlab.unittest.TestCase
    properties (SetAccess = private)
        folder (1,1) string
        executable (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function buildKernel(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_QG_KERNEL_DUMP"));
            testCase.providers = "reference";
            testCase.providers = [testCase.providers,"reference-pruned","reference-compact"];
            if getenv("WV_QG_TEST_NATIVE") == "1"
                testCase.providers = [testCase.providers,"native","native-accelerate","native-pruned","native-accelerate-pruned","native-compact","native-accelerate-compact"];
            end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,text] = systemWithoutMatlabRuntime("cmake -S " + shellQuote(fullfile(root,"PortableRuntime")) + " -B " + shellQuote(build) + " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,text)
                [status,text] = systemWithoutMatlabRuntime("cmake --build " + shellQuote(build) + " --parallel 4 --target WVStratifiedQGKernelDump");
                testCase.assertEqual(status,0,text)
                testCase.executable = fullfile(build,"WVStratifiedQGKernelDump");
            end
            testCase.assertTrue(isfile(testCase.executable))
        end
    end
    methods (Test, TestTags = "full")
        function variableStratificationMatchesMatlab(testCase)
            profiles = {@(z) 1e-4*exp(z/700), @(z) 5e-5*(1+0.4*tanh((z+400)/180))};
            sizes = [8 6 9 4; 9 7 10 6; 8 7 7 1];
            for iProfile = 1:numel(profiles)
                for iSize = 1:size(sizes,1)
                    for antialias = [false true]
                        sz = sizes(iSize,:);
                        wvt = WVTransformStratifiedQG([19000 11000 1000],sz(1:3),Nj=sz(4),N2Function=profiles{iProfile},shouldAntialias=antialias);
                        testCase.compareKernel(wvt)
                    end
                end
            end
        end

        function retainedSubsetMatchesMatlab(testCase)
            base = WVTransformStratifiedQG([19000 11000 1000],[9 7 10],Nj=6,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=false);
            keep = [1 3 5];
            wvt = WVTransformStratifiedQG([base.Lx base.Ly base.Lz],[base.Nx base.Ny base.Nz],N2Function=base.N2Function,shouldAntialias=false,z=base.z,j=base.j(keep),Nj=numel(keep),dLnN2=base.dLnN2,PF0inv=base.PF0inv(:,keep),QG0inv=base.QG0inv(:,keep),PF0=base.PF0(keep,:),QG0=base.QG0(keep,:),P0=base.P0(keep),Q0=base.Q0(keep),h_0=base.h_0(keep),z_int=base.z_int);
            testCase.verifyEqual(wvt.j,[0;2;4])
            testCase.compareKernel(wvt)
        end

        function horizontalMeansAreNotGeostrophicFields(testCase)
            wvt = WVTransformStratifiedQG([19000 11000 1000],[8 6 9],Nj=4,N2Function=@(z) 1e-4*exp(z/700));
            wvt.A0 = complex(zeros(wvt.spectralMatrixSize));
            wvt.A0(:,1) = 1e-6*(1:wvt.Nj)';
            input = testCase.inputs(wvt);
            for provider = testCase.providers
                report = testCase.runKernel(wvt,input,provider);
                testCase.verifyEqual(report.fields.qgpv.value,zeros(wvt.Nx*wvt.Ny*wvt.Nz,1))
                testCase.verifyEqual(report.fields.u.value,zeros(wvt.Nx*wvt.Ny*wvt.Nz,1))
                testCase.verifyEqual(report.energy,0)
                testCase.verifyEqual(report.enstrophy,0)
                testCase.verifyEqual(complexOutput(report.stationary),wvt.A0(:))
            end
        end
    end
    methods (Access = private)
        function input = inputs(~,wvt)
            n = reshape(0:wvt.Nx*wvt.Ny*wvt.Nz-1,wvt.spatialMatrixSize);
            input.A0 = struct('real',real(wvt.A0(:)),'imag',imag(wvt.A0(:)));
            input.q = reshape(1e-6*(sin(0.07*n)+0.2*cos(0.19*n)),[],1);
            input.u = reshape(0.01*sin(0.11*n),[],1);
            input.v = reshape(0.02*cos(0.17*n),[],1);
            input.eta = reshape(0.3*sin(0.09*n),[],1);
        end

        function report = runKernel(testCase,wvt,input,provider)
            path = fullfile(testCase.folder,"qg.nc");
            file = wvt.writeToFile(char(path),shouldOverwriteExisting=true);
            file.close();
            inputPath = fullfile(testCase.folder,"input.json");
            outputPath = fullfile(testCase.folder,"output.json");
            fid = fopen(inputPath,'w');
            cleanup = onCleanup(@()fclose(fid));
            fwrite(fid,jsonencode(input),'char');
            clear cleanup
            [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.executable) + " " + shellQuote(path) + " " + shellQuote(inputPath) + " " + shellQuote(outputPath) + " " + provider);
            testCase.assertEqual(status,0,text)
            report = jsondecode(fileread(outputPath));
        end

        function compareKernel(testCase,wvt)
            n = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
            wvt.A0 = 1e-6*(sin(0.7*n)+1i*cos(0.3*n))./(1+wvt.J);
            wvt.A0(:,wvt.k==0 & wvt.l==0) = 0;
            wvt.t0 = 17;
            wvt.t = 83;
            input = testCase.inputs(wvt);
            fields = ["u","v","w","eta","pi","p","psi","qgpv","rho_e","rho_total","zeta_z","ssh","ssu","ssv"];
            for provider = testCase.providers
                report = testCase.runKernel(wvt,input,provider);
                testCase.verifyEqual(string(report.contract),"wave-vortex-stratified-qg-kernel-v1")
                schedule=string(report.horizontalSchedule);
                testCase.verifyEqual(schedule,expectedHorizontalSchedule(provider))
                testCase.verifyEqual(report.streamedNonlinear,endsWith(provider,["-pruned","-compact"]))
                testCase.verifyEqual(report.compactSplitViews,endsWith(provider,"-compact"))
                testCase.verifyEqual(report.pointwiseWorkers,1+double(endsWith(provider,"-compact")))
                testCase.verifyEqual(report.preparedAllocations,0)
                testCase.verifyTrue(report.inputPreserved)
                testCase.verifyEqual(report.realScratchBytes,4*wvt.Nx*wvt.Ny*wvt.Nz*8)
                factors = {"eta",wvt.NA0; "pi",wvt.PA0; "psi",wvt.A0_Psi_factor; "qgpv",wvt.A0_QGPV_factor; "energy",wvt.A0_TE_factor; "enstrophy",wvt.A0_TZ_factor; "ke",wvt.A0_KE_factor; "pe",wvt.A0_PE_factor};
                for i = 1:size(factors,1)
                    testCase.near(report.factors.(factors{i,1}),factors{i,2},"factor "+factors{i,1})
                end
                testCase.near(complexOutput(report.factors.u),wvt.UA0,"factor u")
                testCase.near(complexOutput(report.factors.v),wvt.VA0,"factor v")
                for field = fields
                    surface = any(field==["ssh","ssu","ssv"]);
                    base = field;
                    if field=="ssh", base="pi"; end
                    if field=="ssu", base="u"; end
                    if field=="ssv", base="v"; end
                    if field=="w"
                        value = zeros(wvt.spatialMatrixSize);
                    else
                        value = wvt.variableWithName(char(base));
                    end
                    dx = wvt.diffX(value); dy = wvt.diffY(value);
                    if base=="rho_total"
                        % The background is horizontally constant. Differentiate
                        % the anomaly to avoid cancellation against rho_nm0.
                        dx=wvt.diffX(wvt.rho_e); dy=wvt.diffY(wvt.rho_e);
                    end
                    if any(base==["eta","rho_e","rho_total"])
                        dz = wvt.diffZG(wvt.eta);
                        if base~="eta"
                            dz = (wvt.rho0/wvt.g)*reshape(wvt.N2,1,1,[]).*(dz+reshape(wvt.dLnN2,1,1,[]).*wvt.eta);
                            if base=="rho_total", dz=dz-(wvt.rho0/wvt.g)*reshape(wvt.N2,1,1,[]); end
                        end
                    else
                        dz = wvt.diffZF(value);
                    end
                    if surface
                        value=value(:,:,end); dx=dx(:,:,end); dy=dy(:,:,end); dz=dz(:,:,end);
                    end
                    testCase.near(report.fields.(field).value,value,field+" value")
                    testCase.near(report.fields.(field).x,dx,field+" dx")
                    testCase.near(report.fields.(field).y,dy,field+" dy")
                    testCase.near(report.fields.(field).z,dz,field+" dz")
                end
                testCase.near(complexOutput(report.projectQ),wvt.transformQGPVToWaveVortex(reshape(input.q,wvt.spatialMatrixSize)),"QGPV projection")
                [~,~,projected] = wvt.transformUVEtaToWaveVortex(reshape(input.u,wvt.spatialMatrixSize),reshape(input.v,wvt.spatialMatrixSize),reshape(input.eta,wvt.spatialMatrixSize));
                testCase.near(complexOutput(report.projectUVEta),projected,"UVEta projection")
                flux = wvt.nonlinearFlux();
                testCase.near(complexOutput(report.flux),flux,"nonlinear flux")
                betaFlux = wvt.transformQGPVToWaveVortex(-wvt.beta*wvt.v);
                testCase.near(complexOutput(report.linearBeta),betaFlux,"beta flux")
                betaForcing = WVBetaPlanePVAdvection(wvt);
                wvt.addForcing(betaForcing);
                testCase.near(complexOutput(report.fluxBeta),wvt.nonlinearFlux(),"combined flux")
                wvt.removeForcing(betaForcing);
                testCase.near(complexOutput(report.roundtripQ),wvt.A0,"QGPV round trip")
                testCase.near(complexOutput(report.roundtripUVEta),wvt.A0,"UVEta round trip")
                testCase.near(complexOutput(report.stationary),wvt.A0,"stationary A0")
                testCase.near(complexOutput(report.evolvedBeta),wvt.A0.*exp(-wvt.beta*wvt.VA0*1234),"Rossby evolution")
                testCase.near(report.energy,wvt.totalEnergy,"spectral energy")
                testCase.near(report.enstrophy,wvt.totalEnstrophy(),"spectral enstrophy")
                testCase.near(report.spatialEnergy,wvt.totalEnergySpatiallyIntegrated,"spatial energy")
                testCase.near(report.spatialEnstrophy,wvt.totalEnstrophySpatiallyIntegrated(),"spatial enstrophy")
                testCase.near(report.uvMax,wvt.uvMax,"maximum speed")
                closures = {WVVerticalDiffusivity(wvt,kappa_z=.002),WVBottomFrictionLinear(wvt,r=1e-5),WVBottomFrictionQuadratic(wvt,Cd=.003)};
                closureNames = ["verticalDiffusivity","linearBottomFriction","quadraticBottomFriction"];
                for iClosure = 1:numel(closures)
                    spatialFlux = closures{iClosure}.addPotentialVorticitySpatialForcing(wvt,zeros(wvt.spatialMatrixSize));
                    testCase.near(complexOutput(report.(closureNames(iClosure))),wvt.transformQGPVToWaveVortex(spatialFlux),closureNames(iClosure))
                end
            end
        end

        function near(testCase,actual,expected,label)
            testCase.assertEqual(numel(actual),numel(expected),label)
            testCase.verifyLessThanOrEqual(max(abs(actual(:)-expected(:))),2e-11*max(abs(expected),[],"all")+1e-18,label)
        end
    end
end

function value = complexOutput(record)
value = record.real + 1i*record.imag;
end
function value = shellQuote(value)
value = "'" + replace(string(value),"'","'""'""'") + "'";
end
function [status,output] = systemWithoutMatlabRuntime(command)
if isunix && ~ismac
    command = "env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH " + command;
end
[status,output] = system(command);
end

function schedule = expectedHorizontalSchedule(provider)
if endsWith(provider,["-pruned","-compact"]) && startsWith(provider,"native")
    schedule="fftw-streaming-pruned-tile16";
else
    schedule="full-fft-gather";
end
end
