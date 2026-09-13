function manifest = writeVariableModelAdoptionFixtures(folder,options)
% Write immutable SQG and Hydrostatic full-model qualification fixtures.
arguments (Input)
    folder (1,1) string
    options.profileSet (1,1) string {mustBeMember(options.profileSet,["small","large"])} = "small"
end
arguments (Output)
    manifest (1,1) struct
end

assert(~isfolder(folder),"Refusing to overwrite a fixture directory.");
folder = string(java.io.File(char(folder)).getCanonicalPath());
mkdir(folder);
generatorPath = string(mfilename("fullpath"))+".m";
root = string(fileparts(fileparts(fileparts(fileparts(generatorPath)))));
originalFolder = string(pwd);
cd(root);
folderCleanup = onCleanup(@()cd(originalFolder));
[status,sourceHead] = system("git -C "+shellQuote(root)+" rev-parse HEAD");
assert(status==0,"Unable to identify the WVM source commit.");

if options.profileSet=="small"
    grid = [32 24 33];
    finalTime = 2;
    outputInterval = 1;
    denseInterval = .5;
    denseTime = [.75 1.25];
    minimumFreeBytes = 0;
else
    grid = [256 256 129];
    finalTime = .5;
    outputInterval = .5;
    denseInterval = .25;
    denseTime = .25;
    % This is an authoring guard, not a model-memory prediction. Exact file
    % sizes and free space are recorded after each generated profile.
    minimumFreeBytes = 8*2^30;
end
step = .5;
families = ["stratified-qg","hydrostatic"];
profileCells = cell(1,numel(families));

for index = 1:numel(families)
    family = families(index);
    probe = makeTransform(family,[8 6 9]);
    if family=="stratified-qg"
        probeFlux = probe.nonlinearFlux();
    else
        [probeFp,probeFm,probeF0] = probe.nonlinearFlux();
        probeFlux = probeFp;
    end
    assert(all(isfinite(probeFlux),"all"), ...
        "WaveVortexModel:FixturePreflight","Small nonlinear preflight failed.");
    clear probe probeFlux probeFp probeFm probeF0
    usableBefore = double(java.io.File(char(folder)).getUsableSpace());
    if usableBefore < minimumFreeBytes
        error("WaveVortexModel:FixtureCapacity", ...
            "Large fixture authoring requires at least %.3f GiB free; %.3f GiB is available.", ...
            minimumFreeBytes/2^30,usableBefore/2^30);
    end
    identity = "variable-"+family+"-composite-"+options.profileSet;
    sourcePath = fullfile(folder,identity+"-source.nc");
    referencePath = fullfile(folder,identity+"-matlab-reference.nc");
    model = makeModel(family,grid);
    Nj = model.wvt.Nj;
    Nkl = model.wvt.Nkl;
    writeInitial(model,family,sourcePath,finalTime,outputInterval,denseInterval,denseTime);
    model.closeNetCDFFile();
    clear model

    copyfile(sourcePath,referencePath);
    reference = WVModel.modelFromFile(referencePath);
    reference.outputTimesForIntegrationPeriod(0,finalTime);
    reference.setupIntegrator(integratorType="fixed",deltaT=step);
    reference.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
    reference.closeNetCDFFile();
    clear reference

    assert(isequal(ncread(referencePath,"/wave-vortex/t"), ...
        (0:outputInterval:finalTime)'),"WaveVortexModel:FixtureSchedule", ...
        "Unexpected coefficient output times.");
    assert(isequal(ncread(referencePath,"/dense/t"),denseTime(:)), ...
        "WaveVortexModel:FixtureSchedule","Unexpected dense output times.");
    assert(isequal(ncread(referencePath,"/tracers/t"), ...
        (0:outputInterval:finalTime)'),"WaveVortexModel:FixtureSchedule", ...
        "Unexpected tracer output times.");
    hasParticles = family=="hydrostatic";
    if hasParticles
        assert(isequal(ncread(referencePath,"/particles/t"), ...
            (0:outputInterval:finalTime)'),"WaveVortexModel:FixtureSchedule", ...
            "Unexpected particle output times.");
    end

    sourceInfo = dir(sourcePath);
    referenceInfo = dir(referencePath);
    realVolumeBytes = prod(grid)*8;
    modalComplexBytes = Nj*Nkl*16;
    coefficientFamilyCount = 1+2*double(family=="hydrostatic");
    coefficientRecordBytes = coefficientFamilyCount*modalComplexBytes;
    ordinaryRecords = numel(0:outputInterval:finalTime);
    estimatedOutputPayloadBytes = ordinaryRecords*(coefficientRecordBytes+realVolumeBytes) + ...
        numel(denseTime)*realVolumeBytes + ordinaryRecords*realVolumeBytes;
    profileCells{index} = struct( ...
        id=identity,family=family,grid=grid,Nj=Nj,Nkl=Nkl,sourcePath=sourcePath, ...
        referencePath=referencePath,sourceSHA256=portableCompatibilitySHA256(sourcePath), ...
        referenceSHA256=portableCompatibilitySHA256(referencePath), ...
        sourceBytes=sourceInfo.bytes,referenceBytes=referenceInfo.bytes, ...
        usableBytesBefore=usableBefore, ...
        usableBytesAfter=double(java.io.File(char(folder)).getUsableSpace()), ...
        capabilities=struct(fullXYZTracer=true,tracerAdvection=conditional( ...
            family=="stratified-qg","horizontal-xy","three-dimensional-xyz"), ...
            movingParticles=family=="hydrostatic",denseEulerianU=true), ...
        capacityEstimate=struct(realVolumeBytes=realVolumeBytes, ...
            modalComplexBytes=modalComplexBytes, ...
            coefficientFamilyCount=coefficientFamilyCount, ...
            coefficientRecordBytes=coefficientRecordBytes, ...
            outputPayloadLowerBoundBytes=estimatedOutputPayloadBytes), ...
        expectedStepCount=round(finalTime/step), ...
        expectedRHSEvaluationCount=4*round(finalTime/step), ...
        expectedCoefficientRecords=numel(0:outputInterval:finalTime), ...
        expectedDenseRecords=numel(denseTime), ...
        expectedParticleRecords=double(hasParticles)*numel(0:outputInterval:finalTime), ...
        expectedTracerRecords=numel(0:outputInterval:finalTime));
    saveJSON(fullfile(folder,identity+"-fixture.json"),profileCells{index});
    fprintf("VARIABLE_MODEL_FIXTURE_COMPLETE %s sourceGiB=%.3f referenceGiB=%.3f freeGiB=%.3f\n", ...
        identity,sourceInfo.bytes/2^30,referenceInfo.bytes/2^30,profileCells{index}.usableBytesAfter/2^30);
end
profiles = [profileCells{:}];

realVolumeBytes = prod(grid)*8;
halfComplexVolumeUpperBoundBytes = (floor(grid(1)/2)+1)*grid(2)*grid(3)*16;
manifest = struct( ...
    schema="wvm-variable-model-fixtures-"+options.profileSet+"-v1", ...
    generatorPath=generatorPath,generatorSHA256=portableCompatibilitySHA256(generatorPath), ...
    sourceHEAD=strtrim(string(sourceHead)),matlabVersion=string(version), ...
    parameters=struct(grid=grid,Lxyz=[19000 11000 1000], ...
        N2Function="1e-4*exp(z/700)",shouldAntialias=true, ...
        forcing="WVNonlinearAdvection",finalTime=finalTime,initialStep=step, ...
        outputInterval=outputInterval,denseOutputInterval=denseInterval, ...
        denseOutputTimes=denseTime, ...
        fields={{"u"}},particleX=[1900 11400],particleY=[1100 7700], ...
        particleZ=[-250 -750],tracer="dye"), ...
    capacityPreflight=struct( ...
        policy="refuse rather than substitute a smaller grid", ...
        minimumAuthoringFreeBytes=minimumFreeBytes, ...
        singleRealVolumeBytes=realVolumeBytes, ...
        halfComplexVolumeUpperBoundBytes=halfComplexVolumeUpperBoundBytes, ...
        scopedFamilies={{"stratified-qg","hydrostatic"}}, ...
        scopeExclusion=struct(family="boussinesq",grid=[256 256 129], ...
            reason="Boussinesq 256 is feasible by capacity but requires a separate full-model qualification fixture."), ...
        capacityExclusion=struct(grid=[512 512 257], ...
            reason="The measured Boussinesq grouped scientific arrays alone exceed available disk; four retained model outputs also approach or exceed the remaining capacity.", ...
            boussinesqGroupedScientificArrayBytes=20.286523699760437*2^30, ...
            availableBytesAtEstimate=15.881553649902344*2^30, ...
            sqgFourOutputPayloadLowerBoundBytes=4*2838951360, ...
            hydroOrBoussFourOutputPayloadLowerBoundBytes=4*3127173440), ...
        boussinesq256NarrowQualification=struct( ...
            groupedScientificArrayBytes=1.290598213672638*2^30, ...
            modelOutputPayloadLowerBoundBytes=392274720, ...
            firstBlockFourOutputPayloadLowerBoundBytes=4*392274720, ...
            proposedWork="one fixed RK4 step with XYZ tracer, particles, endpoint coefficients/u, and one dense u event")), ...
    profiles=profiles);
saveJSON(fullfile(folder,"manifest.json"),manifest);
end

function model = makeModel(family,grid)
model = WVModel(makeTransform(family,grid));
wvt = model.wvt;
model.eulerianObservingSystem.addNetCDFOutputVariables("u");
phi = sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly);
if family=="hydrostatic"
    model.setFloatPositions([1900 11400],[1100 7700],[-250 -750],'u');
    model.addTracer(phi,"dye");
else
    tracer = WVTracer(model,name="dye",phi=phi,isXYOnly=true);
    model.addFluxedObservingSystem(tracer);
end
end

function wvt = makeTransform(family,grid)
constructors = dictionary( ...
    "stratified-qg",{@WVTransformStratifiedQG}, ...
    "hydrostatic",{@WVTransformHydrostatic});
constructor = constructors{family};
wvt = constructor([19000 11000 1000],grid, ...
    N2Function=@(z)1e-4*exp(z/700),shouldAntialias=true);
n = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
inertial = wvt.K.^2+wvt.L.^2==0;
wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*~inertial;
if family=="hydrostatic"
    wave = wvt.J>0 & ~inertial;
    wvt.Ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
    wvt.Am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
    wvt.Am(inertial) = conj(wvt.Ap(inertial));
    mda = inertial & wvt.J>0;
    wvt.A0(mda) = .003*sin(n(mda));
end
wvt.t = 0;
wvt.t0 = 0;
wvt.setForcing(WVNonlinearAdvection(wvt));
end

function writeInitial(model,family,path,finalTime,outputInterval,denseInterval,denseTime)
file = model.createNetCDFFileForModelOutput(path, ...
    outputInterval=outputInterval,shouldOverwriteExisting=false);
defaultGroup = file.outputGroupWithName(model.defaultOutputGroupName());
tracer = model.fluxedObservingSystemWithName("dye");
if family=="hydrostatic"
    particles = model.fluxedObservingSystemWithName("float");
    defaultGroup.removeObservingSystem([particles tracer]);
else
    defaultGroup.removeObservingSystem(tracer);
end
dense = file.addNewEvenlySpacedOutputGroup("dense", ...
    outputInterval=denseInterval, ...
    initialTime=denseTime(1),finalTime=denseTime(end));
dense.addObservingSystem(WVEulerianFields(model,fieldNames="u"));
if family=="hydrostatic"
    particleGroup = file.addNewEvenlySpacedOutputGroup("particles", ...
        outputInterval=outputInterval,initialTime=0,finalTime=finalTime);
    particleGroup.addObservingSystem(particles);
end
tracerGroup = file.addNewEvenlySpacedOutputGroup("tracers", ...
    outputInterval=outputInterval,initialTime=0,finalTime=finalTime);
tracerGroup.addObservingSystem(tracer);
file.outputTimesForIntegrationPeriod(0,finalTime);
file.writeTimeStepToOutputFile(0);
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function saveJSON(path,value)
file = fopen(path,"w");
if file<0, error("WaveVortexModel:FixtureWrite","Unable to write %s",path); end
cleanup = onCleanup(@()fclose(file));
fwrite(file,jsonencode(value,PrettyPrint=true));
fwrite(file,newline);
end

function value = conditional(condition,ifTrue,ifFalse)
if condition, value=ifTrue; else, value=ifFalse; end
end
