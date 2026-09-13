function folder = writePortableDiagnosticBaselineFixtures(folder)
% Write deterministic, complete-integration inputs for diagnostic nonregression.
arguments
    folder (1,1) string = string(tempname)
end
if ~isfolder(folder), mkdir(folder); end
for family = ["constant","hydrostatic","boussinesq"]
    switch family
        case "constant"
            wvt = WVTransformConstantStratification([150e3 100e3 1000],[32 24 33],N0=5.2e-3,isHydrostatic=false,shouldAntialias=true);
        case "hydrostatic"
            wvt = WVTransformHydrostatic([150e3 100e3 1000],[24 18 25],Nj=12,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=true);
        case "boussinesq"
            wvt = WVTransformBoussinesq([150e3 100e3 1000],[16 12 17],Nj=8,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=true);
    end
    n = reshape(1:numel(wvt.A0),size(wvt.A0));
    wvt.A0 = 1e-6*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J); wvt.A0(:,wvt.k==0 & wvt.l==0)=0;
    wave = wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0; inertial = (wvt.K.^2+wvt.L.^2)==0;
    wvt.Ap = .001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
    wvt.Am = .002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave; wvt.Am(inertial)=conj(wvt.Ap(inertial));
    wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]);
    model = WVModel(wvt);
    model.eulerianObservingSystem.addNetCDFOutputVariables('u','v','w','eta','p','qgpv','ssu','ssv','ssh');
    source = fullfile(folder,family+"-source.nc");
    file = model.createNetCDFFileForModelOutput(source,outputInterval=2,shouldOverwriteExisting=true);
    dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.5,initialTime=0,finalTime=32);
    dense.addObservingSystem(WVEulerianFields(model,fieldNames={'u','eta'}));
    file.outputTimesForIntegrationPeriod(0,32); file.writeTimeStepToOutputFile(0); model.closeNetCDFFile();
    request = fullfile(folder,family+"-request.json");
    WVModel.writePortableRunRequest(request,source,method="fixed-rk4",finalTime=32,initialStep=.5,fftProvider="native-fftw",reportPath=family+"-report.json");
end
end
