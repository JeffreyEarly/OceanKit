%% Mean-density-anomaly modes and their sampled transform
% Mean-density-anomaly (MDA) modes represent the horizontally averaged
% interior displacement. Their G structures are the prognostic expansion
% functions. The aligned F structures diagnose mean pressure by integrating
% N^2 G downward from the surface.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Define an exponentially stratified water column
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
g0 = 0.03;
gd = 0;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

% Both endpoint parameters are finite, so both endpoints are active. The
% bottom value gd=0 is the Neumann limit. Positive infinity would instead
% impose G=0 and remove that endpoint from the generalized energy.
evp = IMInternalModes.meanDensityAnomalyModes( ...
    N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd);
nAvailableModes = 20;
solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);

fprintf("Active endpoints: %s\n",join(basisSet.activeEndpoints,", "));
fprintf("Mode numbers: %s\n",join(string(basisSet.modeNumber),", "));
fprintf("Generalized-energy signatures: %s\n",join(string(basisSet.signatures),", "));

%% Inspect the aligned G and F structures
% The solved G modes satisfy the generalized-energy endpoint conditions.
% The diagnostic pressure modes are
%
%   F_j(z) = (1/g) integral_z^0 N^2(z') G_j(z') dz'.
%
% Consequently every F mode is exactly zero at the surface, including the
% pressure structure paired with the constant zero-eigenvalue G mode.
zPlot = linspace(zDomain(1),zDomain(2),801).';
nPlotModes = 5;
GPlot = basisSet.G(zPlot);
FPlot = basisSet.F(zPlot);

figure(Name="Mean-density-anomaly modes",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(GPlot(:,1:nPlotModes),zPlot,LineWidth=1.4)
grid on
xlabel("G_j(z)")
ylabel("z (m)")
title("Displacement modes")
legend("j="+string(basisSet.modeNumber(1:nPlotModes)),Location="best")

nexttile
plot(FPlot(:,1:nPlotModes),zPlot,LineWidth=1.4)
grid on
xlabel("F_j(z)")
ylabel("z (m)")
title("Surface-referenced mean pressure")
legend("j="+string(basisSet.modeNumber(1:nPlotModes)),Location="best")

%% Build the sampled MDA transform
% Omitting variables selects the directly projectable channel. For MDA
% modes this is G. F remains available for synthesis from the same aligned
% coefficient vector, but it does not define a separate projection metric.
z = linspace(zDomain(1),zDomain(2),97).';
nTransformModes = 18;
[transform,assessment] = basisSet.discreteTransform(z=z,nModes=nTransformModes);

fprintf("The transform uses %d points and retains %d of %d candidate modes.\n", ...
    length(transform.z),assessment.retainedModeCount,assessment.candidateModeCount);
fprintf("Direct projection channel: %s\n",join(transform.availableVariables,", "));
fprintf("F forward transform available: %d\n",transform.hasForwardTransform(variable="F"));

%% Project G and synthesize the aligned pressure field
% A sampled displacement is projected through G. The recovered coefficients
% can reconstruct both displacement and mean pressure on the same grid.
j = transform.modeNumber(:);
coefficients = exp(-abs(j)/4).*cos(pi*j/7);
G = transform.inverseMatrix(variable="G");
F = transform.inverseMatrix(variable="F");

displacement = G*coefficients;
recoveredCoefficients = transform.transformForward(displacement,variable="G");
meanPressure = transform.transformBack(recoveredCoefficients,variable="F");
coefficientError = norm(recoveredCoefficients-coefficients)/norm(coefficients);

fprintf("G coefficient round-trip error: %.3e\n",coefficientError);
fprintf("Maximum reconstructed F value at the surface: %.3e\n",max(abs(meanPressure(end,:))));

figure(Name="MDA sampled reconstruction",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(displacement,transform.z,"k-",LineWidth=1.6)
grid on
xlabel("mean displacement")
ylabel("z (m)")
title("G-channel synthesis")

nexttile
plot(meanPressure,transform.z,"Color",[0.85 0.33 0.1],LineWidth=1.6)
grid on
xlabel("mean pressure structure")
ylabel("z (m)")
title("Aligned F-channel synthesis")
