%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
Lz = 4000;
Lxy = 500e3;
Nxy = 128;
Nz = WVStratification.verticalResolutionForHorizontalResolution(Lxy,Lz,Nxy,N2=N2,latitude=27);
wvt = WVTransformHydrostatic([Lxy, Lxy, Lz],[Nxy, Nxy, Nz], N2=N2,latitude=27);
wvt.addForcing(WVAdaptiveDamping(wvt));

wvt.initWithGMSpectrum();
model = WVModel(wvt);

model.createNetCDFFileForModelOutput("igw-simulation-day0-4.nc",outputInterval=wvt.inertialPeriod/4,shouldOverwriteExisting=true);

model.integrateToTime(wvt.t + 4*wvt.inertialPeriod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Integrate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


model = WVModel.modelFromFile("igw-simulation-day0-6.nc");
model.integrateToTime(model.wvt.t + 2*model.wvt.inertialPeriod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Diagnostics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wvd = WVDiagnostics("igw-simulation-day0-4.nc");
wvd.createDiagnosticsFile();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Diagnostics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wvd = WVDiagnostics("igw-simulation-day0-6.nc");
wvd.createDiagnosticsFile();