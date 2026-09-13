function records = hydroCompositePrepare(outputDirectory,options)
% Recreate the two scientifically matched hydrostatic audit benchmark fixtures.
arguments
    outputDirectory (1,1) string
    options.dependencyPaths (1,:) string = strings(1,0)
end

repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@()restoreState(originalPath,originalRng));
addRepositoryPaths(repositoryRoot,options.dependencyPaths);
if ~isfolder(outputDirectory), mkdir(outputDirectory); end

Lxyz = [150e3 150e3 1300];
Nxyz = [256 256 129];
wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2Function=@(z)2e-5*exp(2*z/Lxyz(3)),latitude=45,shouldAntialias=true);
[state,evidence] = initializeThreeInterfaceIntegratorState(wvt,4001,gmEnergyLevel=0.5,initialConditionId="gm0p5-red-geostrophic-j1-v1");
wvt.throwErrorIfDensityViolation(A0=wvt.A0,Ap=wvt.Ap,Am=wvt.Am);
delete(wvt);

workloads = ["coefficient-endpoint" "composite-dense-output"];
records = repmat(struct("path","","workload","","sha256","","bytes",0,"initialCondition",struct()),2,1);
for iWorkload = 1:2
    pathname = fullfile(outputDirectory,"matched-hydrostatic-exponential-"+workloads(iWorkload)+"-model.nc");
    if isfile(pathname), error("WaveVortexBenchmark:HydroAuditFixtureExists","Refusing to overwrite fixture: %s",pathname); end
    createFixture(pathname,workloads(iWorkload),Lxyz,Nxyz,state);
    information = dir(pathname);
    records(iWorkload) = struct("path",pathname,"workload",workloads(iWorkload),"sha256",sha256File(pathname),"bytes",information.bytes,"initialCondition",evidence);
end
clear stateCleanup
end

function createFixture(pathname,workload,Lxyz,Nxyz,state)
wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2Function=@(z)2e-5*exp(2*z/Lxyz(3)),latitude=45,shouldAntialias=true);
advanceWaveVortexBenchmarkState(wvt,state,0);
model = WVModel(wvt);
modelCleanup = onCleanup(@()closeModel(model));
isComposite = workload == "composite-dense-output";
if isComposite
    model.setFloatPositions([10e3 70e3],[9e3 65e3],[-250 -850],'u');
    model.addTracer(sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly),"dye");
end
outputFile = model.createNetCDFFileForModelOutput(pathname,outputInterval=7168,shouldOverwriteExisting=true);
if isComposite
    group = outputFile.outputGroupWithName(model.defaultOutputGroupName());
    particles = model.fluxedObservingSystemWithName("float");
    tracer = model.fluxedObservingSystemWithName("dye");
    group.removeObservingSystem([particles tracer]);
    denseGroup = outputFile.addNewEvenlySpacedOutputGroup("dense",outputInterval=32,initialTime=0,finalTime=128);
    denseGroup.addObservingSystem(WVEulerianFields(model,fieldNames={'u'}));
    denseGroup.addObservingSystem(WVMooring(model,name="mooring",x=[0 50e3],y=[0 40e3],trackedFieldNames={'u'}));
    particleGroup = outputFile.addNewEvenlySpacedOutputGroup("particles",outputInterval=32,initialTime=0,finalTime=128);
    particleGroup.addObservingSystem(particles);
    tracerGroup = outputFile.addNewEvenlySpacedOutputGroup("tracers",outputInterval=32,initialTime=0,finalTime=128);
    tracerGroup.addObservingSystem(tracer);
end
model.outputTimesForIntegrationPeriod(model.t,model.t);
model.writeTimeStepToNetCDFFile(model.t);
model.recordNetCDFFileHistory();
model.closeNetCDFFile();
clear modelCleanup
end

function addRepositoryPaths(repositoryRoot,dependencyPaths)
for dependencyPath = reshape(dependencyPaths,1,[])
    if ~isfolder(dependencyPath), error("WaveVortexBenchmark:HydroAuditDependency","Dependency folder does not exist: %s",dependencyPath); end
    addpath(dependencyPath);
end
metadata = jsondecode(fileread(fullfile(repositoryRoot,"resources","mpackage.json")));
for item = reshape(metadata.folders,1,[])
    folder = fullfile(repositoryRoot,item.path);
    if isfolder(folder), addpath(folder); end
end
addpath(repositoryRoot);
addpath(fullfile(repositoryRoot,"Benchmarks"));
end

function value = sha256File(pathname)
[status,output] = system(sprintf('/usr/bin/shasum -a 256 %s',shellQuote(pathname)));
if status ~= 0, error("WaveVortexBenchmark:HydroAuditHashFailed","%s",output); end
value = string(extractBefore(strtrim(output),65));
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function closeModel(model)
if isempty(model) || ~isvalid(model), return, end
try
    model.closeNetCDFFile();
catch
end
delete(model);
end

function restoreState(originalPath,originalRng)
path(originalPath);
rng(originalRng);
end
