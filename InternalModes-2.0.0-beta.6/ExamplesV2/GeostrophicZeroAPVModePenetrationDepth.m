%% Geostrophic zero-APV mode penetration depth in exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g0 = -N0*N0*b;
zDomain = [-4000 0];
N2 = @(z) N0*N0*exp(2*z/b);

wavelength = logspace(log10(1000e3), log10(100), 80);
wavelengthKm = wavelength/1000;
k = 2*pi./wavelength;
nEVP = 512;
z = linspace(zDomain(1), zDomain(2), 4000).';
threshold = 0.01;

problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,f0=f0,k=k,endpoints="surface");
solver = IMSolverSpectral(nEVP=nEVP);
boundaryModes = solver.solveGeostrophicZeroAPVModes(problem);
basisSet = boundaryModes.rotateBoundaryDepth(g0=g0);

h = basisSet.h0;
iSignChange = find(h(1:end-1).*h(2:end) < 0);
kSignChange = zeros(size(iSignChange));
for iCross = 1:length(iSignChange)
    iLeft = iSignChange(iCross);
    fraction = -h(iLeft)/(h(iLeft+1) - h(iLeft));
    logKCrossing = log(k(iLeft)) + fraction*(log(k(iLeft+1)) - log(k(iLeft)));
    kSignChange(iCross) = exp(logKCrossing);
end
wavelengthSignChangeKm = 2*pi./kSignChange/1000;

F = squeeze(basisSet.F(z));
surfaceAmplitude = abs(F(end,:));
relativeAmplitude = abs(F)./surfaceAmplitude;
penetrationDepth = zeros(size(k));
fullDepth = zDomain(2) - zDomain(1);

for iMode = 1:length(k)
    iBelow = find(relativeAmplitude(:,iMode) <= threshold, 1, "last");
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
    amplitudeBelow = relativeAmplitude(iBelow,iMode);
    amplitudeAbove = relativeAmplitude(iBelow+1,iMode);
    if amplitudeAbove == amplitudeBelow
        zCrossing = zBelow;
    else
        fraction = (threshold - amplitudeBelow)/(amplitudeAbove - amplitudeBelow);
        zCrossing = zBelow + fraction*(zAbove - zBelow);
    end
    penetrationDepth(iMode) = zDomain(2) - zCrossing;
end

figure(Name="Geostrophic zero-APV mode penetration depth",Color="w");
tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

topAxes = nexttile;
hPositive = abs(h);
hPositive(h <= 0) = NaN;
hNegative = abs(h);
hNegative(h >= 0) = NaN;
semilogx(wavelengthKm, hPositive, LineWidth=1.5)
hold on
semilogx(wavelengthKm, hNegative, "--", LineWidth=1.5)
hold off
set(topAxes, YScale="log")
grid on
for iCross = 1:length(kSignChange)
    xline(wavelengthSignChangeKm(iCross), "-", LineWidth=0.75, Color=[0.35 0.35 0.35])
end
ylabel("|h| (m)")
topAxes.XTickLabel = [];
% title("Surface-mode eigendepth and penetration depth")
legend(["h > 0" "h < 0"], Location="best")

bottomAxes = nexttile;
semilogx(wavelengthKm, -penetrationDepth, LineWidth=1.5)
grid on
for iCross = 1:length(kSignChange)
    xline(wavelengthSignChangeKm(iCross), "-", LineWidth=0.75, Color=[0.35 0.35 0.35])
end
xlabel("wavelength (km)")
ylabel("penetration depth (m)")
ylim(zDomain)
set([topAxes bottomAxes], XDir="reverse")
linkaxes([topAxes bottomAxes], "x")
xlim(bottomAxes, [min(wavelengthKm) max(wavelengthKm)])
