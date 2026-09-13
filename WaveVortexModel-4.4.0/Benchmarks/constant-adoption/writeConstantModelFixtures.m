function manifest = writeConstantModelFixtures(folder,options)
% Write deterministic constant-stratification adoption fixtures and references.
% profileSet="large" writes six four-step cases at 256/512 horizontal sizes;
% composite cases save full u at endpoints and once off-step, plus particles/dye.
% Existing default small fixtures and integration/output schedules are preserved.
arguments
    folder (1,1) string = string(tempname)
    options.profileSet (1,1) string {mustBeMember(options.profileSet,["small","large"])} = "small"
end
entries = dir(folder); entries = entries(~ismember({entries.name},{'.','..'}));
if isfolder(folder) && ~isempty(entries)
    error("WaveVortexModel:FixtureFolderExists","Refusing to overwrite nonempty fixture folder: %s",folder)
end
if ~isfolder(folder), mkdir(folder); end

profiles = ["constant-hydrostatic-coefficient-only" "constant-hydrostatic-composite" ...
    "constant-nonhydrostatic-coefficient-only" "constant-nonhydrostatic-composite"];
large = options.profileSet == "large";
Nxyz = [32 24 33]; finalTime = 32; outputInterval = 2; denseInterval = .5; fields = {'u','v','w','eta'};
if large
    profiles = ["constant-hydrostatic-coefficient-only-256" "constant-nonhydrostatic-coefficient-only-256" "constant-nonhydrostatic-composite-256" ...
        "constant-hydrostatic-coefficient-only-512" "constant-nonhydrostatic-coefficient-only-512" "constant-nonhydrostatic-composite-512"];
    finalTime = 2; outputInterval = 2; denseInterval = .75; fields = {'u'};
end
generatorPath = string(which("writeConstantModelFixtures"));
sourceHead = string(strtrim(runCommand("git -C "+quote(fileparts(fileparts(fileparts(generatorPath))))+" rev-parse HEAD")));
records = repmat(struct("id","","sourcePath","","referencePath","","requestPath","","sourceSHA256","","referenceSHA256","","requestSHA256","","isHydrostatic",false,"composite",false,"outputCounts",struct()),numel(profiles),1);
for index = 1:numel(profiles)
    id = profiles(index); composite = contains(id,"composite"); hydrostatic = startsWith(id,"constant-hydrostatic-");
    sourcePath = fullfile(folder,id+"-source.nc"); referencePath = fullfile(folder,id+"-reference.nc"); requestPath = fullfile(folder,id+"-request.json");
    if large
        if endsWith(id,"256"), Nxyz = [256 256 129]; else, Nxyz = [512 512 257]; end
    end
    if large
        availableBytes = double(java.io.File(char(folder)).getUsableSpace());
        fprintf("MODEL_FIXTURE_START %s freeGiB=%.3f\n",id,availableBytes/2^30);
        assert(availableBytes>=20*2^30,"WaveVortexModel:FixtureDiskSpace","Large model authoring requires 20 GiB free before each profile; completed profiles are preserved.");
    end
    model = makeModel(hydrostatic,composite,Nxyz,fields); writeInitial(model,sourcePath,composite,finalTime,outputInterval,denseInterval,fields,large); model.closeNetCDFFile();
    clear model
    copyfile(sourcePath,referencePath); reference = WVModel.modelFromFile(referencePath); reference.outputTimesForIntegrationPeriod(0,finalTime); reference.setupIntegrator(integratorType="fixed",deltaT=.5); reference.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]); reference.closeNetCDFFile();
    if large
        assert(isequal(ncread(referencePath,"/wave-vortex/t"),[0;2]),"WaveVortexModel:FixtureSchedule","Unexpected coefficient output times.");
        if composite
            assert(isequal(ncread(referencePath,"/dense/t"),.75),"WaveVortexModel:FixtureSchedule","Expected exactly one off-step dense record.");
            assert(isequal(ncread(referencePath,"/particles/t"),[0;2]) && isequal(ncread(referencePath,"/tracers/t"),[0;2]),"WaveVortexModel:FixtureSchedule","Unexpected dynamic-state output times.");
        end
    end
    WVModel.writePortableRunRequest(requestPath,sourcePath,method="fixed-rk4",finalTime=finalTime,initialStep=.5,fftProvider="native-fftw",threads=16,reportPath=fullfile(folder,id+"-report.json"));
    records(index)=struct("id",id,"sourcePath",sourcePath,"referencePath",referencePath,"requestPath",requestPath,"sourceSHA256",sha256(sourcePath),"referenceSHA256",sha256(referencePath),"requestSHA256",sha256(requestPath),"isHydrostatic",hydrostatic,"composite",composite,"outputCounts",counts(composite,finalTime,outputInterval,denseInterval,large));
    clear model reference
    if large
        saveJSON(fullfile(folder,id+"-fixture.json"),records(index));
        fprintf("MODEL_FIXTURE_COMPLETE %s freeGiB=%.3f\n",id,double(java.io.File(char(folder)).getUsableSpace())/2^30);
    end
end
manifest=struct("schema","wvm-constant-model-fixtures-v1","generatorPath",generatorPath,"generatorSHA256",sha256(generatorPath),"sourceHEAD",sourceHead,"matlabVersion",string(version),"parameters",struct("Nxyz",[32 24 33],"Lxyz",[150e3 100e3 1000],"N0",5.2e-3,"shouldAntialias",true,"finalTime",32,"initialStep",.5,"outputInterval",2,"denseOutputInterval",.5,"threads",16,"forcing","WVNonlinearAdvection","fields",{{"u","v","w","eta"}},"particleX",[10e3 70e3],"particleY",[9e3 65e3],"particleZ",[-250 -850],"tracer","dye"),"profiles",records);
if large
    manifest.schema = "wvm-constant-model-fixtures-large-v1";
    manifest.parameters.Nxyz = [256 256 129;512 512 257];
    manifest.parameters.finalTime = finalTime;
    manifest.parameters.outputInterval = outputInterval;
    manifest.parameters.denseOutputInterval = denseInterval;
    manifest.parameters.denseOutputTimes = .75;
    manifest.parameters.fields = fields;
    manifest.parameters.fixedStepCount = 4;
    manifest.parameters.rhsEvaluationCount = 16;
    manifest.parameters.compositeScope = "Nonhydrostatic at each size: full u at endpoints, one dense u at .75, two moving particles and active dye tracer; no hydrostatic composite claim.";
end
saveJSON(fullfile(folder,"manifest.json"),manifest);
end

function model=makeModel(hydrostatic,composite,Nxyz,fields)
wvt=WVTransformConstantStratification([150e3 100e3 1000],Nxyz,N0=5.2e-3,isHydrostatic=hydrostatic,shouldAntialias=true);
n=reshape(1:numel(wvt.A0),size(wvt.A0)); wvt.A0=1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J); wvt.A0(:,wvt.k==0 & wvt.l==0)=0;
wave=wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0; inertial=(wvt.K.^2+wvt.L.^2)==0; wvt.Ap=.001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial); wvt.Am=.002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave; wvt.Am(inertial)=conj(wvt.Ap(inertial)); wvt.A0(inertial & wvt.J>0)=.003*sin(n(inertial & wvt.J>0)); wvt.t=0; wvt.t0=0; wvt.setForcing(WVNonlinearAdvection(wvt));
model=WVModel(wvt);
if composite, model.eulerianObservingSystem.addNetCDFOutputVariables(fields{:}); model.setFloatPositions([10e3 70e3],[9e3 65e3],[-250 -850],'u'); model.addTracer(sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly),"dye"); end
end

function writeInitial(model,path,composite,finalTime,outputInterval,denseInterval,fields,large)
file=model.createNetCDFFileForModelOutput(path,outputInterval=outputInterval,shouldOverwriteExisting=true);
if composite
    group=file.outputGroupWithName(model.defaultOutputGroupName()); particles=model.fluxedObservingSystemWithName("float"); tracer=model.fluxedObservingSystemWithName("dye"); group.removeObservingSystem([particles tracer]);
    denseInitial = 0; denseFinal = finalTime;
    if large, denseInitial = .75; denseFinal = .75; end
    dense=file.addNewEvenlySpacedOutputGroup("dense",outputInterval=denseInterval,initialTime=denseInitial,finalTime=denseFinal); dense.addObservingSystem(WVEulerianFields(model,fieldNames=fields));
    particleGroup=file.addNewEvenlySpacedOutputGroup("particles",outputInterval=outputInterval,initialTime=0,finalTime=finalTime); particleGroup.addObservingSystem(particles);
    tracerGroup=file.addNewEvenlySpacedOutputGroup("tracers",outputInterval=outputInterval,initialTime=0,finalTime=finalTime); tracerGroup.addObservingSystem(tracer);
end
file.outputTimesForIntegrationPeriod(0,finalTime); file.writeTimeStepToOutputFile(0);
end

function value=counts(composite,finalTime,outputInterval,denseInterval,large)
ordinary = 1+floor(finalTime/outputInterval); dense = 1+floor(finalTime/denseInterval);
if large, dense = 1; end
value=struct("waveVortex",ordinary,"dense",double(composite)*dense,"particles",double(composite)*ordinary,"tracers",double(composite)*ordinary);
end

function value=sha256(path)
value=extractBefore(string(strtrim(runCommand("/usr/bin/shasum -a 256 "+quote(path))))," ");
end

function value=runCommand(command)
[status,value]=system(command); if status~=0, error("WaveVortexModel:FixtureCommand","Command failed: %s",value); end
end

function value=quote(path)
value="'"+replace(string(path),"'","'""'""'")+"'";
end

function saveJSON(path,value)
file=fopen(path,"w"); if file<0, error("WaveVortexModel:FixtureWrite","Unable to write %s",path); end
cleanup=onCleanup(@()fclose(file)); fwrite(file,jsonencode(value)); fwrite(file,newline);
end
