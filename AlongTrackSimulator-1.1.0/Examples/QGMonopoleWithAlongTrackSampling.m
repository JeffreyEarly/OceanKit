%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 2000e3;
Ly = 1000e3;

Nx = 2*256;
Ny = 2*128;

latitude = 24;

wvt = WVTransformBarotropicQG([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2),shouldRemoveMeanPressure=1 );

figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
wvt.addForcing(WVAdaptiveDamping(wvt));
wvt.addForcing(WVBetaPlanePVAdvection(wvt));
model = WVModel(wvt);

outputFile = model.createNetCDFFileForModelOutput('QGMonopoleWithAlongTrack.nc',outputInterval=86400,shouldOverwriteExisting=1);
model.addNetCDFOutputVariables("ssh","zeta_z")
ats = AlongTrackSimulator();
ats.addMissionsToOutputFile(outputFile);


%%
model.integrateToTime(365*86400);
ncfile = model.ncfile;

%%
[obs.x,obs.y,obs.t,obs.ssh] = ncfile.readVariables("s6a/track_x","s6a/track_y","s6a/t","s6a/ssh");
figure, scatter3(obs.x/1e3,obs.y/1e3,obs.ssh*1e2,[],obs.ssh*1e2,'filled'), colorbar('eastoutside')

[obs.x,obs.y,obs.t,obs.ssh] = ncfile.readVariables("j3n/track_x","j3n/track_y","j3n/t","j3n/ssh");
figure, scatter3(obs.x/1e3,obs.y/1e3,obs.ssh*1e2,[],obs.ssh*1e2,'filled'), colorbar('eastoutside')

[obs.x,obs.y,obs.t,obs.ssh] = ncfile.readVariables("alg/track_x","alg/track_y","alg/t","alg/ssh");
figure, scatter3(obs.x/1e3,obs.y/1e3,obs.ssh*1e2,[],obs.ssh*1e2,'filled'), colorbar('eastoutside')
