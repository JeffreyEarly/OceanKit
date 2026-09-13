function manifest = writeConstantAdoptionFixtures(folder,options)
% Write current MATLAB states and independent complete-flux references for #358.
arguments (Input)
    folder (1,1) string
    options.grids (:,3) double {mustBeInteger,mustBePositive} = [256 256 129;512 512 257]
    options.hydrostatic (1,:) logical = [false true]
end
arguments (Output)
    manifest (1,1) struct
end
if ~isfolder(folder), mkdir(folder); end
root = string(fileparts(fileparts(fileparts(mfilename("fullpath")))));
[status,commit] = system("git -C "+shellQuote(root)+" rev-parse HEAD"); assert(status==0);
manifest = struct(schema="wvm-constant-adoption-fixtures-v1",sourceCommit=strtrim(string(commit)), ...
    matlabRelease=string(version),generatorSHA256=portableCompatibilitySHA256(mfilename("fullpath")+".m"),cases={{}});
originalRng = rng; cleanup = onCleanup(@()rng(originalRng));
for gridIndex = 1:size(options.grids,1)
    grid = options.grids(gridIndex,:);
    for isHydrostatic = options.hydrostatic
        families = ["nonhydrostatic","hydrostatic"];
        identity = "constant-"+families(double(isHydrostatic)+1);
        identity = identity+"-"+join(string(grid),"x");
        sourcePath = fullfile(folder,identity+".nc"); fluxPath = fullfile(folder,identity+"-matlab-flux.bin");
        assert(~isfile(sourcePath) && ~isfile(fluxPath),"Refusing to overwrite an existing frozen fixture.");
        rng(20020,"twister");
        wvt = WVTransformConstantStratification([15000 15000 1300],grid,N0=5.2e-3,latitude=45,rotationRate=7.2921e-5,g=9.81,shouldAntialias=true,isHydrostatic=isHydrostatic,computationalBackend="matlab");
        shape = wvt.spectralMatrixSize;
        wvt.Ap = complex((2*rand(shape)-1)/8,(2*rand(shape)-1)/8);
        wvt.Am = complex((2*rand(shape)-1)/8,(2*rand(shape)-1)/8);
        wvt.A0 = complex((2*rand(shape)-1)/8,(2*rand(shape)-1)/8);
        dc = find(wvt.kMode_wv==0 & wvt.lMode_wv==0,1);
        wvt.Am(:,dc) = conj(wvt.Ap(:,dc)); wvt.A0(:,dc) = real(wvt.A0(:,dc));
        wvt.t = wvt.t0+123.5; wvt.setForcing(WVNonlinearAdvection(wvt));
        [Fp,Fm,F0] = wvt.nonlinearFlux();
        assert(all(isfinite(Fp),"all") && all(isfinite(Fm),"all") && all(isfinite(F0),"all"));
        file = fopen(fluxPath,"w","ieee-le"); assert(file>=0); fileCleanup = onCleanup(@()fclose(file));
        for value = {Fp,Fm,F0}
            array = value{1}; interleaved = [real(array(:))';imag(array(:))'];
            assert(fwrite(file,interleaved,"double")==2*numel(array));
        end
        clear fileCleanup Fp Fm F0 array interleaved value
        model = WVModel(wvt);
        output = model.createNetCDFFileForModelOutput(sourcePath,outputInterval=1,shouldOverwriteExisting=false);
        output.outputTimesForIntegrationPeriod(wvt.t,wvt.t); output.writeTimeStepToOutputFile(wvt.t); model.closeNetCDFFile();
        record = struct(id=identity,grid=grid,isHydrostatic=isHydrostatic,antialias=true,Nj=wvt.Nj,Nkl=wvt.Nkl, ...
            t=wvt.t,t0=wvt.t0,seed=20020,sourcePath=sourcePath,sourceSHA256=fileHash(sourcePath), ...
            matlabFluxPath=fluxPath,matlabFluxSHA256=fileHash(fluxPath),fluxLayout="Fp,Fm,F0; each canonical column-major complex interleaved float64 little-endian", ...
            maximumScaleNormalizedTolerance=1e-10,relativeL2Tolerance=1e-10);
        manifest.cases{end+1} = record;
        fprintf("Wrote %s with Nj=%d Nkl=%d\n",identity,wvt.Nj,wvt.Nkl);
        clear model wvt output
    end
end
file = fopen(fullfile(folder,"manifest.json"),"w"); assert(file>=0); fileCleanup = onCleanup(@()fclose(file));
fwrite(file,jsonencode(manifest,PrettyPrint=true));
end
function hash = fileHash(path)
[status,value] = system("shasum -a 256 "+shellQuote(path)); assert(status==0); hash=extractBefore(string(value)," ");
end
function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
