%% Free-surface gravity wave modes at fixed wavenumber

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

D = 4000;
N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g = 9.81;
zDomain = [-D 0];

k = 2*pi/1e3;
nModes = 4;
nEVP = 128;
z = linspace(zDomain(1), zDomain(2), 512).';

background = IMBackground.exponentialStratification(N0=N0,b=b,zDomain=zDomain,rotationRate=abs(f0)/(2*sind(33)),latitude=33*sign(f0),g=g,surfaceCondition="freeSurface");
evp = IMInternalModes.waveModesAtWavenumber(background,k=k);

solver = IMSolverSpectral(nEVP=nEVP, coordinateKind="wkb");
basisSet = solver.solveEVP(evp, nModes=nModes);

%%
basisSet.normalization = "surfaceF";


G = basisSet.G(z);
F = basisSet.F(z);

h = basisSet.h(:);
phaseSpeed = sqrt(g*h);
omega = sqrt(f0*f0 + g*h*k*k);
periodHours = 2*pi./omega/3600;
surfaceResidual = G(end,:).' - (h/D).*F(end,:).';
bottomResidual = G(1,:).';
diagnostics = table(basisSet.modeNumber(:), h, phaseSpeed, omega, periodHours, ...
    surfaceResidual, bottomResidual, ...
    VariableNames=["modeNumber" "h" "phaseSpeed" "omega" "periodHours" ...
    "surfaceResidual" "bottomResidual"]);

fprintf("Free-surface gravity wave modes at fixed wavenumber\n");
fprintf("k = %.6g m^-1, wavelength = %.3g km\n\n", k, 2*pi/k/1000);
disp(diagnostics)

modeLabels = "mode " + string(basisSet.modeNumber) + ", T = " + compose("%.2g", periodHours.') + " hr";

figure(Name="V2 free-surface gravity wave modes: exponential stratification", Color="w");
tiledlayout(1, 2, TileSpacing="compact", Padding="compact");

nexttile
plot(G, z, LineWidth=1.2)
grid on
xlabel("G")
ylabel("z (m)")
title("G modes")
legend(modeLabels, Location="best")

nexttile
plot(F, z, LineWidth=1.2)
grid on
xlabel("F")
title("F diagnostics")
legend(modeLabels, Location="best")
