%% Hydrostatic modes in standard exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

N0 = 5.2e-3;
L_gm = 1.3e3;
zDomain = [-5000 0];
N2 = @(z) N0*N0*exp(2*z/L_gm);

nModes = 4;
nEVP = 96;
z = linspace(zDomain(1), zDomain(2), 256).';

evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, surfaceBoundary=IMBoundaryCondition.dirichlet(), ...
    bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=nEVP);
basisSet = solver.solveEVP(evp, nModes=nModes);

G = basisSet.G(z);
F = basisSet.F(z);
modeLabels = "mode " + string(1:nModes);

figure(Name="V2 hydrostatic modes: exponential stratification", Color="w");
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
