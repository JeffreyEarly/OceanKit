classdef TestHydrostaticCompiledKernel < matlab.unittest.TestCase
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
            testCase.executable = string(getenv("WV_HYDRO_KERNEL_DUMP"));
            testCase.providers = ["reference","reference-pruned","reference-compact"];
            if getenv("WV_HYDRO_TEST_NATIVE") == "1"
                testCase.providers = [testCase.providers,"native","native-accelerate","native-pruned","native-accelerate-pruned","native-compact","native-accelerate-compact"];
            end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,text] = systemWithoutMatlabRuntime("cmake -S " + shellQuote(fullfile(root,"PortableRuntime")) + " -B " + shellQuote(build) + " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,text)
                [status,text] = systemWithoutMatlabRuntime("cmake --build " + shellQuote(build) + " --parallel 4 --target WVHydrostaticKernelDump");
                testCase.assertEqual(status,0,text)
                testCase.executable = fullfile(build,"WVHydrostaticKernelDump");
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
                        wvt = WVTransformHydrostatic([19000 11000 1000],sz(1:3),Nj=sz(4),N2Function=profiles{iProfile},shouldAntialias=antialias);
                        testCase.compareKernel(wvt)
                    end
                end
            end
        end
        function retainedSubsetAndSouthernLatitude(testCase)
            base = WVTransformHydrostatic([19000 11000 1000],[9 7 10],Nj=6,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=false,latitude=-33);
            keep = [1 3 5];
            wvt = WVTransformHydrostatic([base.Lx base.Ly base.Lz],[base.Nx base.Ny base.Nz],N2Function=base.N2Function,shouldAntialias=false,latitude=-33,z=base.z,j=base.j(keep),Nj=numel(keep),dLnN2=base.dLnN2,PF0inv=base.PF0inv(:,keep),QG0inv=base.QG0inv(:,keep),PF0=base.PF0(keep,:),QG0=base.QG0(keep,:),P0=base.P0(keep),Q0=base.Q0(keep),h_0=base.h_0(keep),z_int=base.z_int);
            testCase.verifyEqual(wvt.j,[0;2;4])
            testCase.compareKernel(wvt)
        end
    end
    methods (Access = private)
        function report = runKernel(testCase,wvt,input,provider)
            path = fullfile(testCase.folder,"hydrostatic.nc");
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
            wave = wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0;
            inertial = (wvt.K.^2+wvt.L.^2)==0;
            geo = ~inertial;
            mda = inertial & wvt.J>0;
            ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
            am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
            am(inertial) = conj(ap(inertial));
            a0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*geo;
            a0(mda) = .003*sin(n(mda));
            initial = {ap,am,a0};
            wvt.Ap=ap; wvt.Am=am; wvt.A0=a0; wvt.t0=17; wvt.t=83;
            input.t=wvt.t; input.t0=wvt.t0;
            names=["Ap","Am","A0"];
            for j=1:3
                input.(names(j))=struct('real',real(initial{j}(:)),'imag',imag(initial{j}(:)));
            end
            n=reshape(0:wvt.Nx*wvt.Ny*wvt.Nz-1,wvt.spatialMatrixSize);
            input.u=reshape(.01*sin(.11*n),[],1);
            input.v=reshape(.02*cos(.17*n),[],1);
            input.eta=reshape(.3*sin(.09*n),[],1);
            % Irregular full-grid values retain horizontal frequencies outside
            % the wave-vortex spectrum in the raw vertical calculus checks.
            input.q=reshape(sin(.07*n)+.2*cos(.19*n),[],1);
            fields=["u","v","w","eta","pi","p","psi","qgpv","rho_e","rho_total","zeta_x","zeta_y","zeta_z","ssh","ssu","ssv"];
            [fp,fm,f0]=wvt.nonlinearFlux(); flux={fp,fm,f0};
            [pp,pm,p0]=wvt.transformUVEtaToWaveVortex(reshape(input.u,wvt.spatialMatrixSize),reshape(input.v,wvt.spatialMatrixSize),reshape(input.eta,wvt.spatialMatrixSize)); projected={pp,pm,p0};
            trajectory=testCase.trajectory(wvt,initial);
            wvt.Ap=ap; wvt.Am=am; wvt.A0=a0; wvt.t=83;
            for provider=testCase.providers
                report=testCase.runKernel(wvt,input,provider);
                testCase.verifyEqual(string(report.contract),"wave-vortex-hydrostatic-kernel-v1")
                schedule=string(report.horizontalSchedule);
                testCase.verifyEqual(schedule,expectedHorizontalSchedule(provider))
                testCase.verifyEqual(report.streamedNonlinear,endsWith(provider,["-pruned","-compact"]))
                testCase.verifyEqual(report.compactSplitViews,endsWith(provider,"-compact"))
                testCase.verifyEqual(report.pointwiseWorkers,1+double(endsWith(provider,"-compact")))
                testCase.verifyEqual(report.preparedAllocations,0)
                testCase.verifyTrue(report.inputPreserved && report.storageStable && report.ownerReleased)
                expectedScratchBytes=10*wvt.Nx*wvt.Ny*wvt.Nz*8;
                if endsWith(provider,["-pruned","-compact"]), expectedScratchBytes=6*wvt.Nx*wvt.Ny*wvt.Nz*8; end
                testCase.verifyEqual(report.realScratchBytes,expectedScratchBytes)
                for name=["Omega","NAp","NA0","PA0","ApmN","A0Z","A0N","Apm_TE_factor","A0_TE_factor","A0_Psi_factor","A0_QGPV_factor","A0_TZ_factor"]
                    testCase.near(report.factors.(name),wvt.(name),"factor "+name)
                end
                for name=["UAp","VAp","WAp","UA0","VA0","ApmD"]
                    testCase.near(complexOutput(report.factors.(name)),wvt.(name),"factor "+name)
                end
                for field=fields
                    surface=any(field==["ssh","ssu","ssv"]);
                    base=field;
                    if field=="ssh", base="pi"; end
                    if field=="ssu", base="u"; end
                    if field=="ssv", base="v"; end
                    value=wvt.variableWithName(char(base));
                    testCase.near(report.fields.(field).value,surfaceValue(value,surface),field+" value")
                    if any(field==["zeta_x","zeta_y"]), continue; end
                    dx=wvt.diffX(value); dy=wvt.diffY(value);
                    if base=="rho_total"
                        dx=wvt.diffX(wvt.rho_e); dy=wvt.diffY(wvt.rho_e);
                    end
                    if any(base==["w","eta"])
                        dz=wvt.diffZG(value);
                    elseif any(base==["rho_e","rho_total"])
                        dz=(wvt.rho0/wvt.g)*reshape(wvt.N2,1,1,[]).*(wvt.diffZG(wvt.eta)+reshape(wvt.dLnN2,1,1,[]).*wvt.eta);
                        if base=="rho_total", dz=dz-(wvt.rho0/wvt.g)*reshape(wvt.N2,1,1,[]); end
                    else
                        dz=wvt.diffZF(value);
                    end
                    testCase.near(report.fields.(field).x,surfaceValue(dx,surface),field+" dx")
                    testCase.near(report.fields.(field).y,surfaceValue(dy,surface),field+" dy")
                    testCase.near(report.fields.(field).z,surfaceValue(dz,surface),field+" dz")
                end
                evolved={ap.*wvt.phase,am.*wvt.conjPhase,a0};
                for j=1:3
                    testCase.near(complexOutput(report.roundtrip.(names(j))),initial{j},names(j)+" round trip")
                    testCase.near(complexOutput(report.project.(names(j))),projected{j},names(j)+" projection")
                    testCase.near(complexOutput(report.flux.(names(j))),flux{j},names(j)+" flux")
                    testCase.near(complexOutput(report.evolved.(names(j))),evolved{j},names(j)+" phase")
                    testCase.near(complexOutput(report.trajectory.(names(j))),trajectory{j},names(j)+" trajectory")
                end
                testCase.near(report.enstrophy,wvt.totalEnstrophy(),"enstrophy")
                components=["all","wave","inertial","geostrophic","meanDensityAnomaly"];
                masks={wave|inertial,geo|mda; wave,false(size(geo)); inertial,false(size(geo)); false(size(wave)),geo; false(size(wave)),mda};
                for component=1:5
                    wvt.Ap=ap.*masks{component,1}; wvt.Am=am.*masks{component,1}; wvt.A0=a0.*masks{component,2};
                    part=report.components.(components(component));
                    testCase.near(part.energy,wvt.totalEnergy,"component energy "+components(component))
                    testCase.near(part.spatialEnergy,wvt.totalEnergySpatiallyIntegrated,"component spatial energy "+components(component))
                    for field=fields(~ismember(fields,["rho_total","ssh","ssu","ssv"]))
                        testCase.near(part.(field),wvt.variableWithName(char(field)),components(component)+" "+field)
                    end
                end
                wvt.Ap=ap; wvt.Am=am; wvt.A0=a0;
                q=reshape(input.q,wvt.spatialMatrixSize);
                for family=["F","G"]
                    part=report.raw.(family);
                    raw=ap; raw(inertial)=real(raw(inertial));
                    inverse=wvt.("transformToSpatialDomainWith"+family)(Apm=raw);
                    forward=wvt.("transformFromSpatialDomainWith"+family+"g")(wvt.transformFromSpatialDomainWithFourier(q));
                    testCase.near(part.inverse,inverse,family+" inverse")
                    testCase.near(complexOutput(part.forward),forward,family+" forward")
                    for order=1:4
                        testCase.near(part.("d"+order),wvt.("diffZ"+family)(q,n=order),family+" derivative "+order)
                    end
                    testCase.near(part.integral,wvt.("intZ"+family)(q),family+" integral")
                end
                cp=(.1+.2i)*(wave|inertial); cm=(.1+.2i)*wave+(.1-.2i)*inertial; c0=(.1+.2i)*geo+.1*mda;
                testCase.near(complexOutput(report.constrained.Ap),cp,"constrained Ap")
                testCase.near(complexOutput(report.constrained.Am),cm,"constrained Am")
                testCase.near(complexOutput(report.constrained.A0),c0,"constrained A0")
            end
        end
        function result=trajectory(~,wvt,initial)
            result=initial; dt=.5;
            for step=0:3
                slopes=cell(4,3);
                for stage=1:4
                    fractions=[0 .5 .5 1]; f=fractions(stage);
                    values=result;
                    if stage>1
                        for j=1:3, values{j}=values{j}+dt*f*slopes{stage-1,j}; end
                    end
                    wvt.Ap=values{1}; wvt.Am=values{2}; wvt.A0=values{3}; wvt.t=83+dt*(step+f);
                    [slopes{stage,1},slopes{stage,2},slopes{stage,3}]=wvt.nonlinearFlux();
                end
                for j=1:3
                    result{j}=result{j}+dt/6*(slopes{1,j}+2*slopes{2,j}+2*slopes{3,j}+slopes{4,j});
                end
            end
        end
        function near(testCase,actual,expected,label)
            testCase.assertEqual(numel(actual),numel(expected),label)
            testCase.verifyLessThanOrEqual(max(abs(actual(:)-expected(:))),2e-10*max(abs(expected),[],"all")+1e-18,label)
        end
    end
end
function value=surfaceValue(value,surface)
if surface, value=value(:,:,end); end
end
function value=complexOutput(record)
value=record.real+1i*record.imag;
end
function value=shellQuote(value)
value="'"+replace(string(value),"'","'""'""'")+"'";
end
function [status,output]=systemWithoutMatlabRuntime(command)
if isunix && ~ismac
    command="env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH "+command;
end
[status,output]=system(command);
end

function schedule = expectedHorizontalSchedule(provider)
if endsWith(provider,["-pruned","-compact"]) && startsWith(provider,"native")
    schedule="fftw-streaming-pruned-tile16";
else
    schedule="full-fft-gather";
end
end
