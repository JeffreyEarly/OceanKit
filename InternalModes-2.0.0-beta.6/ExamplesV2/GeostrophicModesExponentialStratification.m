%% Geostrophic APV modes with exponential stratification
% Generalized-energy APV modes can have positive, zero, or negative
% eigendepths. The public factory defines their F/G relationship, endpoint
% conditions, metadata, and depth normalization together.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Define the water column and APV endpoint parameters
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
g0 = 0.03;
gd = 0.01;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

evp = IMInternalModes.geostrophicAPVModes( ...
    N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd, ...
    surfaceBoundary="freeSurface");

%% Solve and evaluate the aligned F/G family
solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=6);
z = linspace(zDomain(1),zDomain(2),512).';
F = basisSet.F(z);
G = basisSet.G(z);

modeLabels = "j = " + string(basisSet.modeNumber) + ...
    ", h = " + compose("%.3g m",basisSet.h);

figure(Name="V2 geostrophic APV modes",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(F,z,LineWidth=1.2)
grid on
xlabel("F")
ylabel("z (m)")
title("APV F modes")
legend(modeLabels,Location="best")

nexttile
plot(G,z,LineWidth=1.2)
grid on
xlabel("G")
title("Aligned displacement modes")
legend(modeLabels,Location="best")

fprintf("Normalization: %s\n",string(basisSet.normalization));
fprintf("Mode numbers: %s\n",join(string(basisSet.modeNumber),", "));
fprintf("Equivalent depths (m): %s\n",join(compose("%.6g",basisSet.h),", "));
