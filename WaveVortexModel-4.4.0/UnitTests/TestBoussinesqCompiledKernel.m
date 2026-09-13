classdef TestBoussinesqCompiledKernel < matlab.unittest.TestCase
    properties (SetAccess = private)
        folder (1,1) string
        hydrostaticExecutable (1,1) string
        executable (1,1) string
        columnOrder (1,:) double = []
        providers (1,:) string
    end
    methods (TestClassSetup)
        function buildKernel(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.folder = string(fixture.Folder);
            testCase.executable = string(getenv("WV_BOUSS_KERNEL_DUMP"));
            testCase.providers = ["reference","reference-pruned","reference-compact"];
            if getenv("WV_BOUSS_TEST_NATIVE") == "1"
                testCase.providers = [testCase.providers,"native","native-accelerate","native-pruned","native-accelerate-pruned","native-compact","native-accelerate-compact"];
            end
            if testCase.executable == ""
                build = fullfile(testCase.folder,"build");
                [status,text] = systemWithoutMatlabRuntime("cmake -S " + shellQuote(fullfile(root,"PortableRuntime")) + " -B " + shellQuote(build) + " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,text)
                [status,text] = systemWithoutMatlabRuntime("cmake --build " + shellQuote(build) + " --parallel 4 --target WVBoussinesqKernelDump WVHydrostaticKernelDump");
                testCase.assertEqual(status,0,text)
                testCase.executable = fullfile(build,"WVBoussinesqKernelDump");
            end
            testCase.hydrostaticExecutable=fullfile(fileparts(testCase.executable),"WVHydrostaticKernelDump");
            testCase.assertTrue(isfile(testCase.executable))
            testCase.assertTrue(isfile(testCase.hydrostaticExecutable))
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
                        wvt = WVTransformBoussinesq([19000 11000 1000],sz(1:3),Nj=sz(4),N2Function=profiles{iProfile},shouldAntialias=antialias);
                        testCase.compareKernel(wvt)
                    end
                end
            end
        end
        function weakSurfaceStratificationMatchesMatlab(testCase)
            wvt=WVTransformBoussinesq([19000 11000 1000],[8 6 9],Nj=4,N2Function=@(z) 1e-10+1e-4*(z/1000).^2,shouldAntialias=false);
            testCase.verifyLessThan(wvt.N2(end),wvt.f^2)
            testCase.compareKernel(wvt)
        end
        function discontiguousGroupsAndNonuniformGrid(testCase)
            wvt=WVTransformBoussinesq([12000 12000 1000],[8 8 11],Nj=5,z=-1000*linspace(1,0,11)'.^1.3,N2Function=@(z) 5e-5*(1+.4*tanh((z+400)/180)),shouldAntialias=false);
            % MATLAB sorts groups contiguously. Permute the temporary source
            % and coefficient columns to exercise exact discontiguous batches.
            testCase.columnOrder=[1,2:2:wvt.Nkl,3:2:wvt.Nkl];
            membership=wvt.iK2unique(testCase.columnOrder);
            separated=false;
            for group=1:numel(wvt.K2unique), separated=separated || any(diff(find(membership==group))>1); end
            testCase.verifyTrue(separated,"Fixture must exercise discontiguous persisted membership.")
            testCase.compareKernel(wvt)
        end
        function approachesHydrostaticLimit(testCase)
            % Small f/N and decreasing horizontal wavenumber isolate the
            % nonhydrostatic correction without assuming finite-f bases agree.
            lengths=[12000 1200000]; errors=zeros(1,2);
            for index=1:2
                L=lengths(index); profile=@(z) 1e-4*exp(z/700);
                b=WVTransformBoussinesq([L L 1000],[6 6 9],Nj=4,N2Function=profile,rotationRate=1e-7,shouldAntialias=false);
                h=WVTransformHydrostatic([L L 1000],[6 6 9],Nj=4,N2Function=profile,rotationRate=1e-7,shouldAntialias=false,z=b.z);
                ap=zeros(b.spectralMatrixSize); mode=find(b.Kh(1,:)>0,1); ap(2,mode)=.001;
                input=struct(t=0,t0=0);
                for name=["Ap","Am","A0"]
                    a=zeros(size(ap)); if name=="Ap", a=ap; end
                    input.(name)=struct('real',real(a(:)),'imag',imag(a(:)));
                end
                for name=["u","v","w","eta","q"], input.(name)=zeros(prod(b.spatialMatrixSize),1); end
                br=testCase.runKernel(b,input,"reference");
                hr=testCase.runKernel(h,input,"reference",testCase.hydrostaticExecutable);
                difference=0;
                for field=["u","v","w","eta","p"]
                    bv=br.fields.(field).value; hv=hr.fields.(field).value;
                    difference=max(difference,max(abs(bv-hv))/max(abs(hv)));
                end
                sample=sub2ind(b.spectralMatrixSize,2,mode);
                difference=max(difference,abs(br.factors.Omega(sample)/hr.factors.Omega(sample)-1));
                difference=max(difference,abs(br.components.wave.energy/hr.components.wave.energy-1));
                errors(index)=difference;
            end
            testCase.verifyLessThan(errors(2),1e-5)
            testCase.verifyLessThan(errors(2),errors(1)/1000)
        end
        function retainedSubsetAndSouthernLatitude(testCase)
            base = WVTransformBoussinesq([19000 11000 1000],[9 7 10],Nj=6,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=false,latitude=-33);
            keep = [1 3 5];
            wvt = WVTransformBoussinesq([base.Lx base.Ly base.Lz],[base.Nx base.Ny base.Nz],N2Function=base.N2Function,shouldAntialias=false,latitude=-33,z=base.z,j=base.j(keep),Nj=numel(keep),dLnN2=base.dLnN2,PF0inv=base.PF0inv(:,keep),QG0inv=base.QG0inv(:,keep),PF0=base.PF0(keep,:),QG0=base.QG0(keep,:),P0=base.P0(keep),Q0=base.Q0(keep),h_0=base.h_0(keep),z_int=base.z_int,PFpmInv=base.PFpmInv(:,keep,:),QGpmInv=base.QGpmInv(:,keep,:),PFpm=base.PFpm(keep,:,:),QGpm=base.QGpm(keep,:,:),Ppm=base.Ppm(keep,:),Qpm=base.Qpm(keep,:),QGwg=base.QGwg(keep,keep,:),h_pm=base.h_pm(keep,:),K2unique=base.K2unique,iK2unique=base.iK2unique);
            testCase.verifyEqual(wvt.j,[0;2;4])
            testCase.compareKernel(wvt)
        end
    end
    methods (Access = private)
        function report = runKernel(testCase,wvt,input,provider,probe)
            if nargin<5, probe=testCase.executable; end
            path = fullfile(testCase.folder,"boussinesq.nc");
            file = wvt.writeToFile(char(path),shouldOverwriteExisting=true);
            file.close();
            order=testCase.columnOrder;
            if ~isempty(order)
                for name=["k","l","iK2unique"]
                    values=ncread(path,name); ncwrite(path,name,values(order));
                end
                values=ncread(path,"h_pm"); ncwrite(path,"h_pm",values(:,order));
                for name=["Ap","Am","A0"], input.(name)=reorderModalRecord(input.(name),wvt.Nj,order); end
            end
            inputPath = fullfile(testCase.folder,"input.json");
            outputPath = fullfile(testCase.folder,"output.json");
            fid = fopen(inputPath,'w');
            cleanup = onCleanup(@()fclose(fid));
            fwrite(fid,jsonencode(input),'char');
            clear cleanup
            [status,text] = systemWithoutMatlabRuntime(shellQuote(probe) + " " + shellQuote(path) + " " + shellQuote(inputPath) + " " + shellQuote(outputPath) + " " + provider);
            testCase.assertEqual(status,0,text)
            report = jsondecode(fileread(outputPath));
            if ~isempty(order)
                inverse(order)=1:numel(order);
                report.waveGroup=report.waveGroup(inverse);
                for name=reshape(string(fieldnames(report.factors)),1,[])
                    report.factors.(name)=reorderModalRecord(report.factors.(name),wvt.Nj,inverse);
                end
                for result=["roundtrip","roundtrip4","project","project4","evolved","flux","trajectory","constrained"]
                    for name=["Ap","Am","A0"], report.(result).(name)=reorderModalRecord(report.(result).(name),wvt.Nj,inverse); end
                end
                for family=["F","G","Fw","Gw"], report.raw.(family).forward=reorderModalRecord(report.raw.(family).forward,wvt.Nj,inverse); end
            end
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
            input.w=reshape(.03*sin(.23*n),[],1);
            input.eta=reshape(.3*sin(.09*n),[],1);
            % Irregular full-grid values retain horizontal frequencies outside
            % the wave-vortex spectrum in the raw vertical calculus checks.
            input.q=reshape(sin(.07*n)+.2*cos(.19*n),[],1);
            fields=["u","v","w","eta","pi","p","psi","qgpv","rho_e","rho_total","zeta_x","zeta_y","zeta_z","ssh","ssu","ssv"];
            [fp,fm,f0]=wvt.nonlinearFlux(); flux={fp,fm,f0};
            [pp,pm,p0]=wvt.transformUVEtaToWaveVortex(reshape(input.u,wvt.spatialMatrixSize),reshape(input.v,wvt.spatialMatrixSize),reshape(input.eta,wvt.spatialMatrixSize)); projected={pp,pm,p0};
            [pp,pm,p0]=wvt.transformUVWEtaToWaveVortex(reshape(input.u,wvt.spatialMatrixSize),reshape(input.v,wvt.spatialMatrixSize),reshape(input.w,wvt.spatialMatrixSize),reshape(input.eta,wvt.spatialMatrixSize)); projected4={pp,pm,p0};
            [pp,pm,p0]=wvt.transformUVEtaToWaveVortex(wvt.u,wvt.v,wvt.eta); roundtrip={pp,pm,p0};
            [pp,pm,p0]=wvt.transformUVWEtaToWaveVortex(wvt.u,wvt.v,wvt.w,wvt.eta); roundtrip4={pp,pm,p0};
            trajectory=testCase.trajectory(wvt,initial);
            wvt.Ap=ap; wvt.Am=am; wvt.A0=a0; wvt.t=83;
            for provider=testCase.providers
                report=testCase.runKernel(wvt,input,provider);
                testCase.verifyEqual(string(report.contract),"wave-vortex-boussinesq-kernel-v1")
                schedule=string(report.horizontalSchedule);
                testCase.verifyEqual(schedule,expectedHorizontalSchedule(provider))
                testCase.verifyEqual(report.streamedNonlinear,endsWith(provider,["-pruned","-compact"]))
                testCase.verifyEqual(report.compactSplitViews,endsWith(provider,"-compact"))
                testCase.verifyEqual(report.pointwiseWorkers,1+double(endsWith(provider,"-compact")))
                testCase.verifyEqual(report.preparedAllocations,0)
                testCase.verifyTrue(report.inputPreserved && report.storageStable && report.ownerReleased)
                expectedScratchVolumes=11;
                if endsWith(provider,["-pruned","-compact"]), expectedScratchVolumes=6; end
                testCase.verifyEqual(report.realScratchBytes,expectedScratchVolumes*wvt.Nx*wvt.Ny*wvt.Nz*8)
                testCase.verifyEqual(report.waveGroup(:),wvt.iK2unique(:)-1)
                testCase.verifyEqual(report.K2unique(:),wvt.K2unique(:))
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
                    testCase.near(complexOutput(report.roundtrip.(names(j))),roundtrip{j},names(j)+" round trip")
                    testCase.near(complexOutput(report.roundtrip4.(names(j))),roundtrip4{j},names(j)+" four-field round trip")
                    testCase.near(complexOutput(report.project4.(names(j))),projected4{j},names(j)+" four-field projection")
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
                for family=["F","G","Fw","Gw"]
                    part=report.raw.(family);
                    raw=ap; raw(inertial)=real(raw(inertial));
                    if family=="F" || family=="G"
                        inverse=wvt.("transformToSpatialDomainWith"+family)(A0=raw);
                        forward=wvt.("transformFromSpatialDomainWith"+family+"g")(wvt.transformFromSpatialDomainWithFourier(q));
                    else
                        inverse=wvt.transformToSpatialDomainWithFourier(wvt.("transformToSpatialDomainWith"+family)(raw));
                        fourier=wvt.transformFromSpatialDomainWithFourier(q); forward=zeros(wvt.spectralMatrixSize);
                        for group=1:numel(wvt.K2unique)
                            columns=wvt.iK2unique==group;
                            if family=="Fw", forward(:,columns)=(wvt.PFpm(:,:,group)*fourier(:,columns))./wvt.Ppm(:,group);
                            else, forward(:,columns)=(wvt.QGpm(:,:,group)*fourier(:,columns))./wvt.Qpm(:,group); end
                        end
                    end
                    testCase.near(part.inverse,inverse,family+" inverse")
                    testCase.near(complexOutput(part.forward),forward,family+" forward")
                    if family=="Fw" || family=="Gw", continue; end
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

function value=reorderModalRecord(value,Nj,order)
if isstruct(value)
    value.real=reorderModalRecord(value.real,Nj,order); value.imag=reorderModalRecord(value.imag,Nj,order);
else
    value=reshape(value,Nj,[]); value=reshape(value(:,order),[],1);
end
end

function schedule = expectedHorizontalSchedule(provider)
if endsWith(provider,["-pruned","-compact"]) && startsWith(provider,"native")
    schedule="fftw-streaming-pruned-tile16";
else
    schedule="full-fft-gather";
end
end
