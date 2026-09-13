function manifest = writeVariableAdoptionFixtures(folder,options)
% Write deterministic MATLAB states and independent variable-transform fluxes.
arguments (Input)
    folder (1,1) string
    options.grid (1,3) double {mustBeInteger,mustBePositive} = [256 256 129]
    options.families (1,:) string {mustBeMember(options.families,["stratified-qg","hydrostatic","boussinesq"])} = ["stratified-qg","hydrostatic"]
end
arguments (Output)
    manifest (1,1) struct
end
assert(~isfolder(folder),"Refusing to overwrite a fixture directory.");
mkdir(folder);
root = string(fileparts(fileparts(fileparts(mfilename("fullpath")))));
[status,commit] = system("git -C "+shellQuote(root)+" rev-parse HEAD"); assert(status==0);
manifest = struct(schema="wvm-variable-screening-fixtures-v1",sourceCommit=strtrim(string(commit)),matlabRelease=string(version),generatorSHA256=portableCompatibilitySHA256(mfilename("fullpath")+".m"),cases={{}});
for family = options.families
    grid = options.grid;
    identity = family+"-"+join(string(grid),"x");
    sourcePath = fullfile(folder,identity+".nc");
    fluxPath = fullfile(folder,identity+"-matlab-flux.bin");
    constructors = dictionary("stratified-qg",{@WVTransformStratifiedQG},"hydrostatic",{@WVTransformHydrostatic},"boussinesq",{@WVTransformBoussinesq});
    constructor = constructors{family};
    wvt = constructor([19000 11000 1000],grid,N2Function=@(z) 1e-4*exp(z/700),shouldAntialias=true);
    n = reshape(1:prod(wvt.spectralMatrixSize),wvt.spectralMatrixSize);
    inertial = wvt.K.^2+wvt.L.^2==0;
    wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*~inertial;
    if family~="stratified-qg"
        wave = wvt.J>0 & ~inertial;
        ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
        am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
        am(inertial) = conj(ap(inertial));
        wvt.Ap = ap; wvt.Am = am;
        mda = inertial & wvt.J>0; wvt.A0(mda) = .003*sin(n(mda));
    end
    clear n inertial wave ap am mda
    wvt.t0=17; wvt.t=83;
    wvt.setForcing(WVNonlinearAdvection(wvt));
    % Persist expensive modal matrices before evaluating the oracle, so a
    % later diagnostic failure leaves a recoverable scientific source.
    model=WVModel(wvt);
    output=model.createNetCDFFileForModelOutput(sourcePath,outputInterval=1,shouldOverwriteExisting=false);
    output.outputTimesForIntegrationPeriod(wvt.t,wvt.t); output.writeTimeStepToOutputFile(wvt.t); model.closeNetCDFFile();
    if family=="stratified-qg"
        values={wvt.nonlinearFlux()};
    else
        [Fp,Fm,F0] = wvt.nonlinearFlux();
        values={Fp,Fm,F0};
    end
    file = fopen(fluxPath,"w","ieee-le"); assert(file>=0); fileCleanup=onCleanup(@()fclose(file));
    for value=values
        array=value{1}; assert(all(isfinite(array),"all"));
        assert(fwrite(file,[real(array(:))';imag(array(:))'],"double")==2*numel(array));
    end
    clear fileCleanup Fp Fm F0 values value array
    record=struct(id=identity,family=family,grid=grid,Nj=wvt.Nj,Nkl=wvt.Nkl,t=wvt.t,t0=wvt.t0,sourcePath=sourcePath,sourceSHA256=portableCompatibilitySHA256(sourcePath),matlabFluxPath=fluxPath,matlabFluxSHA256=portableCompatibilitySHA256(fluxPath),maximumScaleNormalizedTolerance=1e-10,relativeL2Tolerance=1e-10);
    manifest.cases{end+1}=record;
    fprintf("Wrote %s with Nj=%d Nkl=%d\n",identity,wvt.Nj,wvt.Nkl);
    clear model wvt output
end
file=fopen(fullfile(folder,"manifest.json"),"w"); assert(file>=0); fileCleanup=onCleanup(@()fclose(file));
fwrite(file,jsonencode(manifest,PrettyPrint=true));
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
