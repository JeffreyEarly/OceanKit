%% Surface gravity wave mode penetration depth in exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g = 9.81;
zDomain = [-4000 0];
N2 = @(z) N0*N0*exp(2*z/b);

wavelength = logspace(log10(100e3), log10(10), 80);
wavelengthKm = wavelength/1000;
k = 2*pi./wavelength;
nModes = 1;
nEVP = 128;
z = linspace(zDomain(1), zDomain(2), 4000).';
threshold = 0.01;

surfaceBoundary = IMBoundaryCondition(a=0, b=1, c=1, d=0);
bottomBoundary = IMBoundaryCondition.dirichlet();
solver = IMSolverSpectral(nEVP=nEVP, coordinateKind="wkb");

h = zeros(size(k));
periodSeconds = zeros(size(k));
penetrationDepth = zeros(size(k));
fullDepth = zDomain(2) - zDomain(1);

for iMode = 1:length(k)
    evp = IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=k(iMode), f0=f0, g=g, ...
        surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary);
    basisSet = solver.solveEVP(evp, nModes=nModes);
    basisSet.normalization = "surfacePressure";

    h(iMode) = basisSet.h(1);
    omega = sqrt(f0*f0 + g*h(iMode)*k(iMode)*k(iMode));
    periodSeconds(iMode) = 2*pi/omega;

    F = basisSet.F(z);
    relativeAmplitude = abs(F(:,1))/abs(F(end,1));

    iBelow = find(relativeAmplitude <= threshold, 1, "last");
    if isempty(iBelow)
        penetrationDepth(iMode) = fullDepth;
        continue
    end

    if iBelow == length(z)
        penetrationDepth(iMode) = 0;
        continue
    end

    zBelow = z(iBelow);
    zAbove = z(iBelow+1);
    amplitudeBelow = relativeAmplitude(iBelow);
    amplitudeAbove = relativeAmplitude(iBelow+1);
    if amplitudeAbove == amplitudeBelow
        zCrossing = zBelow;
    else
        fraction = (threshold - amplitudeBelow)/(amplitudeAbove - amplitudeBelow);
        zCrossing = zBelow + fraction*(zAbove - zBelow);
    end
    penetrationDepth(iMode) = zDomain(2) - zCrossing;
end

figure(Name="Surface gravity wave mode penetration depth", Color="w");
tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

topAxes = nexttile;
loglog(wavelengthKm, periodSeconds, LineWidth=1.5)
grid on
ylabel("period (s)")
topAxes.XTickLabel = [];
% title("Surface wave-mode period and penetration depth")

bottomAxes = nexttile;
semilogx(wavelengthKm, -penetrationDepth, LineWidth=1.5)
grid on
xlabel("wavelength (km)")
ylabel("penetration depth (m)")
ylim(zDomain)
set([topAxes bottomAxes], XDir="reverse")
linkaxes([topAxes bottomAxes], "x")
xlim(bottomAxes, [min(wavelengthKm) max(wavelengthKm)])
