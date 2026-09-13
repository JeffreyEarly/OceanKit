function constantAdoptionMexWorker(configurationPath)
% Execute a frozen state in one explicitly named isolated MEX, then destroy it.
config = jsondecode(fileread(configurationPath));
addpath(config.moduleFolder,"-begin");
module = string(config.moduleName);
assert(string(which(module))==string(config.modulePath),"Unexpected MEX resolution.");
assert(portableCompatibilitySHA256(config.modulePath)==string(config.moduleSHA256),"Module changed after campaign freeze.");
moduleInfo = feval(module,'moduleInfo');
assert(contains(string(moduleInfo.baseLibrary),string(config.providerRoot)) && contains(string(moduleInfo.threadLibrary),string(config.providerRoot)),"Unexpected FFTW provider.");
assert(string(moduleInfo.openMPRuntimeLibrary)=="","Unexpected OpenMP runtime.");
timer = tic;
[wvt,reader] = WVTransform.waveVortexTransformFromFile(config.sourcePath,iTime=Inf,shouldReadOnly=true);
reader.close(); clear reader
readSeconds = toc(timer);
configuration = struct("Nx",wvt.Nx,"Ny",wvt.Ny,"Nz",wvt.Nz,"Nj",wvt.Nj,"Lx",wvt.Lx,"Ly",wvt.Ly,"Lz",wvt.Lz,"N0",wvt.N0,"rho0",wvt.rho0,"g",wvt.g,"planetaryRadius",wvt.planetaryRadius,"rotationRate",wvt.rotationRate,"latitude",wvt.latitude,"isHydrostatic",wvt.isHydrostatic,"shouldAntialias",wvt.shouldAntialias);
Ap=wvt.Ap; Am=wvt.Am; A0=wvt.A0; t=wvt.t; t0=wvt.t0;
before = feval(module,'moduleMetrics');
timer = tic; handle = feval(module,'create',configuration,config.threads); constructionSeconds=toc(timer);
cleanup = onCleanup(@()feval(module,'delete',handle));
for i=1:config.warmups
    [Fp,Fm,F0] = feval(module,'nonlinearFlux',handle,Ap,Am,A0,t,t0);
end
raw = zeros(1,config.samples);
for i=1:config.samples
    timer = tic;
    [Fp,Fm,F0] = feval(module,'nonlinearFlux',handle,Ap,Am,A0,t,t0);
    raw(i)=toc(timer);
end
file = fopen(config.fluxPath,"w","ieee-le"); assert(file>=0); fileCleanup=onCleanup(@()fclose(file));
for value={Fp,Fm,F0}
    array=value{1}; interleaved=[real(array(:))';imag(array(:))'];
    assert(fwrite(file,interleaved,"double")==2*numel(array));
end
clear fileCleanup value array interleaved
% The saved payload is the last authoritative uninstrumented result.
feval(module,'setStageInstrumentation',handle,true);
[~,~,~] = feval(module,'nonlinearFlux',handle,Ap,Am,A0,t,t0);
metrics = feval(module,'metrics',handle);
clear Fp Fm F0 cleanup
after=feval(module,'moduleMetrics');
assert(after.kernelCount==before.kernelCount && after.activePlans==before.activePlans && after.outstandingPlanningBytes==0,"Kernel/provider ownership did not return to baseline.");
result=struct(status="complete",interface="matlab-loaded",timing=struct(medianSeconds=median(raw),samplesSeconds=raw),readSeconds=readSeconds,constructionSeconds=constructionSeconds,metrics=metrics,moduleInfo=moduleInfo,lifetimeBefore=before,lifetimeAfter=after);
file=fopen(config.resultPath,"w"); assert(file>=0); fileCleanup=onCleanup(@()fclose(file));
fwrite(file,jsonencode(result,PrettyPrint=true));
end
