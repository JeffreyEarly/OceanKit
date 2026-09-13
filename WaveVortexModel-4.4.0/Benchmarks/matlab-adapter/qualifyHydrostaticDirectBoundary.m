function artifact = qualifyHydrostaticDirectBoundary(options)
% Compare one small Hydrostatic state across MATLAB, MEX, and direct native calls.
arguments
    options.executable (1,1) string = string(getenv("WV_HYDRO_KERNEL_DUMP"))
    options.moduleDirectory (1,1) string = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))),".compiled-backend-cache","stage")
    options.workDirectory (1,1) string = fullfile(tempdir,"wvm503","direct-boundary")
    options.outputPath (1,1) string = fullfile(tempdir,"wvm503","hydrostatic-direct-boundary.json")
end
if ~isfile(options.executable)
    error("WaveVortexModel:HydrostaticDirectBoundaryExecutable","The Hydrostatic kernel dump executable is missing: %s",options.executable)
end
modulePath = fullfile(options.moduleDirectory,"wv_compiled_backend_mex."+mexext);
if ~isfile(modulePath)
    error("WaveVortexModel:HydrostaticDirectBoundaryModule","The staged MEX module is missing: %s",modulePath)
end
if ~isfolder(options.workDirectory)
    mkdir(options.workDirectory)
end

rng(503,"twister");
N2Function = @(z)1e-4*exp(z/700);
wvt = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=N2Function,shouldAntialias=false);
wvt.initWithRandomFlow(uvMax=0.01);
wvt.t0 = 2;
wvt.t = 7;
wvtCleanup = onCleanup(@()delete(wvt));
names = ["u" "v" "w" "eta"];
matlabFields = cell(1,numel(names));
fieldNames = cellstr(names);
[matlabFields{:}] = wvt.variableWithName(fieldNames{:});
[matlabFp,matlabFm,matlabF0] = wvt.nonlinearFlux();
matlabFlux = {matlabFp,matlabFm,matlabF0};

recordPath = fullfile(options.workDirectory,"hydrostatic.nc");
record = wvt.writeToFile(char(recordPath),shouldOverwriteExisting=true);
record.close();
inputPath = fullfile(options.workDirectory,"input.json");
directPath = fullfile(options.workDirectory,"direct.json");
input = directInput(wvt,matlabFields);
writeText(inputPath,jsonencode(input));
command = shellQuote(options.executable)+" "+shellQuote(recordPath)+" "+shellQuote(inputPath)+" "+shellQuote(directPath)+" native-accelerate-compact";
[status,output] = system(command);
if status ~= 0
    error("WaveVortexModel:HydrostaticDirectBoundaryExecution","The direct Hydrostatic kernel failed: %s",output)
end
direct = jsondecode(fileread(directPath));

addpath(options.moduleDirectory,"-begin")
pathCleanup = onCleanup(@()rmpath(options.moduleDirectory));
module = "wv_compiled_backend_mex";
clearModule(module)
handle = feval(module,'transformCreate',modalConfiguration(wvt)); %#ok<FVAL>
handleCleanup = onCleanup(@()deleteTransform(module,handle));
feval(module,'transformPrepare',handle,cellstr(names(:))); %#ok<FVAL>
mexResult = feval(module,'transformEvaluate',handle,complex(wvt.Ap),complex(wvt.Am),complex(wvt.A0),wvt.t,wvt.t0,cellstr(names(:)),true); %#ok<FVAL>

fieldComparisons = repmat(struct("name","","matlabDirectRelativeError",NaN,"matlabMexRelativeError",NaN,"directMexRelativeError",NaN,"absoluteTolerancePassed",false),numel(names),1);
for iName = 1:numel(names)
    directValue = reshape(direct.fields.(names(iName)).value,size(matlabFields{iName}));
    fieldComparisons(iName) = comparison(names(iName),matlabFields{iName},directValue,mexResult.values{iName});
end
familyNames = ["Ap" "Am" "A0"];
fluxComparisons = repmat(struct("name","","matlabDirectRelativeError",NaN,"matlabMexRelativeError",NaN,"directMexRelativeError",NaN,"absoluteTolerancePassed",false),3,1);
for iFamily = 1:3
    directValue = complexValue(direct.flux.(familyNames(iFamily)),size(matlabFlux{iFamily}));
    fluxComparisons(iFamily) = comparison(familyNames(iFamily),matlabFlux{iFamily},directValue,mexResult.flux{iFamily});
end
directLifecyclePassed = direct.inputPreserved && direct.storageStable && direct.ownerReleased && direct.preparedAllocations == 0;
allPassed = all([fieldComparisons.absoluteTolerancePassed]) && all([fluxComparisons.absoluteTolerancePassed]) && directLifecyclePassed;
artifact = struct( ...
    "schemaVersion","wvm503-hydrostatic-direct-boundary-v1", ...
    "classification","small direct correctness comparison; not performance qualification", ...
    "configuration",struct("Lxyz",[4000 3000 1000],"Nxyz",[6 6 5],"Nj",3,"N2","1e-4*exp(z/700)","shouldAntialias",false,"seed",503,"uvMax",0.01,"t0",2,"t",7), ...
    "boundaries",struct("matlab","WVTransformHydrostatic","mex","in-memory WVStratifiedModalArrays copied into WVOwnedStratifiedModalSource","direct","NetCDF fixture read by WVHydrostaticKernelDump"), ...
    "directPolicy",struct("engine",string(direct.engine),"matrixBackend",string(direct.backend),"horizontalSchedule",string(direct.horizontalSchedule),"compactSplitViews",direct.compactSplitViews,"pointwiseWorkers",direct.pointwiseWorkers), ...
    "directLifecycle",struct("inputPreserved",direct.inputPreserved,"storageStable",direct.storageStable,"recordOwnerReleased",direct.ownerReleased,"preparedAllocations",direct.preparedAllocations,"passed",directLifecyclePassed), ...
    "mexPolicy",struct("matrixBackend",string(mexResult.metrics.matrixBackend),"policy",string(mexResult.metrics.policy),"horizontalWorkers",mexResult.metrics.horizontalWorkers,"pointwiseWorkers",mexResult.metrics.pointwiseWorkers), ...
    "policyDifference","Both use FFTW, Accelerate, compact split views, and the pruned tiled schedule; the direct test driver fixes horizontal/pointwise workers at 2 while the MEX selects qualified host-topology worker counts.", ...
    "fields",fieldComparisons,"flux",fluxComparisons,"absoluteTolerance","2e-10*max(abs(expected))+1e-18 per existing TestHydrostaticCompiledKernel","passed",allPassed);
writeText(options.outputPath,jsonencode(artifact,PrettyPrint=true));
if ~allPassed
    error("WaveVortexModel:HydrostaticDirectBoundaryMismatch","The small direct Hydrostatic boundary comparison failed; evidence was written to %s.",options.outputPath)
end
clear handleCleanup pathCleanup wvtCleanup
end

function input = directInput(wvt,fields)
input = struct("t",wvt.t,"t0",wvt.t0);
names = ["Ap" "Am" "A0"];
for iName = 1:numel(names)
    value = wvt.(names(iName));
    input.(names(iName)) = struct("real",real(value(:)),"imag",imag(value(:)));
end
input.u = fields{1}(:);
input.v = fields{2}(:);
input.eta = fields{4}(:);
input.q = reshape(sin(0.07*(0:prod(wvt.spatialMatrixSize)-1))+0.2*cos(0.19*(0:prod(wvt.spatialMatrixSize)-1)),[],1);
end

function configuration = modalConfiguration(wvt)
configuration = struct("schemaVersion",'wv-matlab-stratified-modal-source-v1',"transformClass",char(class(wvt)),"modelVersion",char(wvt.version),"Nx",wvt.Nx,"Ny",wvt.Ny,"Nz",wvt.Nz,"Nj",wvt.Nj,"Nkl",wvt.Nkl,"Lx",wvt.Lx,"Ly",wvt.Ly,"Lz",wvt.Lz,"g",wvt.g,"rho0",wvt.rho0,"latitude",wvt.latitude,"rotationRate",wvt.rotationRate,"planetaryRadius",wvt.planetaryRadius,"shouldAntialias",wvt.shouldAntialias,"x",wvt.x,"y",wvt.y,"z",wvt.z,"j",wvt.j,"k",wvt.k,"l",wvt.l,"kMode",wvt.kMode_wv,"lMode",wvt.lMode_wv,"N2",wvt.N2,"rho_nm0",wvt.rho_nm0,"dLnN2",wvt.dLnN2,"z_int",wvt.z_int,"P0",wvt.P0,"Q0",wvt.Q0,"h_0",wvt.h_0,"PF0inv",wvt.PF0inv,"QG0inv",wvt.QG0inv,"PF0",wvt.PF0,"QG0",wvt.QG0);
end

function result = comparison(name,matlabValue,directValue,mexValue)
scale = max(abs(matlabValue),[],"all");
tolerance = 2e-10*scale+1e-18;
result = struct("name",name,"matlabDirectRelativeError",relativeError(matlabValue,directValue),"matlabMexRelativeError",relativeError(matlabValue,mexValue),"directMexRelativeError",relativeError(directValue,mexValue),"absoluteTolerancePassed",max(abs(matlabValue(:)-directValue(:)))<=tolerance && max(abs(matlabValue(:)-mexValue(:)))<=tolerance && max(abs(directValue(:)-mexValue(:)))<=tolerance);
end

function value = complexValue(record,shape)
value = reshape(record.real+1i*record.imag,shape);
end

function value = relativeError(first,second)
value = max(abs(first(:)-second(:)))/max(max(abs(first(:))),realmin);
end

function deleteTransform(module,handle)
feval(module,'transformDelete',handle);
clearModule(module)
end

function clearModule(module)
eval("clear "+module)
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function writeText(pathname,value)
fid = fopen(pathname,"w");
if fid < 0
    error("WaveVortexModel:HydrostaticDirectBoundaryWrite","Unable to write %s.",pathname)
end
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,"%s",value);
clear cleanup
end
