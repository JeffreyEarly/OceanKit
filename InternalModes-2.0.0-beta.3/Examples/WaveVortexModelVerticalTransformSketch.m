%% Sketch vertical transforms for wave-vortex-model integration
% This sketch is intentionally not a public API adapter. It shows how the
% wave-vortex model can assemble its own vertical matrices from generic
% InternalModesBasis and InternalModesTransform objects while keeping the
% full wave-vortex coefficient projectors in wave-vortex-model.

%% Build a hydrostatic grid and transform
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
Nz = 64;
Nj = floor(2*(Nz - 1)/3);
nEVP = max(256,ceil(2.1*Nz));
zDomain = [-Lz 0];

zReference = linspace(zDomain(1),zDomain(2),10*Nz).';
verticalModes = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,nModes=Nz-1,g=g);
verticalModes.normalization = Normalization.geostrophic;
verticalModes.upperBoundary = UpperBoundary.rigidLid;

hydroBasis = InternalModesBasis.fromSolverAtFrequency(verticalModes,0,nModes=Nz-1,useModeAdaptedGrid=true,g=g);
hydroTransform = hydroBasis.modelTransform(nModes=Nj,component="both",nonlinearAliasingPolicy="quadratic",projectionTolerance=1e-2);

PF0 = hydroTransform.forward(component="F");
QG0 = hydroTransform.forward(component="G");
PF0inv = hydroTransform.inverse(component="F");
QG0inv = hydroTransform.inverse(component="G");

fprintf('Hydrostatic vertical transform\n');
fprintf('  PF0 size:    [%d %d]\n', size(PF0,1), size(PF0,2));
fprintf('  QG0 size:    [%d %d]\n', size(QG0,1), size(QG0,2));
fprintf('  PF0inv size: [%d %d]\n', size(PF0inv,1), size(PF0inv,2));
fprintf('  QG0inv size: [%d %d]\n', size(QG0inv,1), size(QG0inv,2));
fprintf('  selection reason: %s\n', hydroTransform.selectionReason);

%% Build a nonzero-kappa wave transform on the hydrostatic grid
kappa = 2*pi/110;
waveModes = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=hydroTransform.z,latitude=latitude,nEVP=nEVP,nModes=Nz-1,g=g);
waveModes.normalization = Normalization.kConstant;
waveModes.upperBoundary = UpperBoundary.rigidLid;

waveBasis = InternalModesBasis.fromSolverAtWavenumber(waveModes,kappa,nModes=Nz-1,g=g);
waveTransform = waveBasis.fixedGridTransform(hydroTransform,component="both",projectionTolerance=1e-2,preserveSize=true);

% G is canonical for nonzero-kappa IGW modes. F is available for numerical
% reconstruction/model mechanics, but it is not a canonical wave-F modal
% projection and must be requested explicitly.
PFpm = waveTransform.forward(component="F",allowNoncanonical=true);
QGpm = waveTransform.forward(component="G");
PFpminv = waveTransform.inverse(component="F");
QGpminv = waveTransform.inverse(component="G");
QGwg = waveTransform.crossTransformTo(hydroTransform,component="G");

fprintf('\nFixed-grid wave vertical transform at kappa %.3e rad/m\n', kappa);
fprintf('  PFpm size:    [%d %d], status %s\n', size(PFpm,1), size(PFpm,2), waveTransform.transformStatusF);
fprintf('  QGpm size:    [%d %d], status %s\n', size(QGpm,1), size(QGpm,2), waveTransform.transformStatusG);
fprintf('  PFpminv size: [%d %d]\n', size(PFpminv,1), size(PFpminv,2));
fprintf('  QGpminv size: [%d %d]\n', size(QGpminv,1), size(QGpminv,2));
fprintf('  QGwg size:    [%d %d]\n', size(QGwg,1), size(QGwg,2));
fprintf('  retained wave G modes: %d\n', length(waveTransform.retainedModesG));
