%% Carry aligned F and G mode families on one discrete transform
% Hydrostatic internal modes come in paired F and G variables. A discrete
% transform uses one set of sample points and quadrature weights for both
% variables while preserving their common physical mode labels and
% equivalent depths. Constant stratification makes this alignment visible:
% F is a cosine family, G is the corresponding sine family, and the
% barotropic member has nonzero F but identically zero G.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve a paired constant-stratification family
D = 4000;
N0 = 5.2e-3;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*ones(size(z));

nAvailableModes = 16;
nModes = 8;
evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain,g=g);
solver = IMSolverSpectral(nEVP=160);
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Put both variables on one endpoint-inclusive rule
% The uniform points and trapezoidal weights are the familiar DCT-I/DST-I
% rule for this exact cosine/sine family. Supplying the weights explicitly
% bypasses fitting and isolates the aligned-transform structure.
M = 12;
z = linspace(zDomain(1),zDomain(2),M+1).';
weights = [0.5;ones(M-1,1);0.5]*(D/M);
variables = ["F","G"];
[transform, assessment] = basisSet.discreteTransform(z=z,weights=weights,nModes=nModes,variables=variables,gramTolerance=1e-10);

FActive = transform.activeModeMask(variable="F").';
GActive = transform.activeModeMask(variable="G").';
family = table(transform.modeNumber.',transform.h.',FActive,GActive,VariableNames=["modeNumber" "h" "FActive" "GActive"]);

fprintf("\nAvailable direct channels: %s\n",join(transform.availableVariables,", "));
fprintf("Both channels use %d points and one weight vector summing to %.1f m.\n\n",length(transform.z),sum(transform.weights));
disp(family)

%% See the barotropic active projector
% Every F column is active. The barotropic G column is identically zero, so
% its direct forward row is also exactly zero. The family is not shortened:
% G keeps the inactive column so all physical labels remain aligned.
FInverse = transform.inverseMatrix(variable="F");
GInverse = transform.inverseMatrix(variable="G");
FProjector = diag(double(transform.activeModeMask(variable="F")));
GProjector = diag(double(transform.activeModeMask(variable="G")));
FProjectorError = norm(transform.forwardMatrix(variable="F")*FInverse-FProjector,2);
GProjectorError = norm(transform.forwardMatrix(variable="G")*GInverse-GProjector,2);

fprintf("F active-projector error: %.3e\n",FProjectorError);
fprintf("G active-projector error: %.3e\n",GProjectorError);
fprintf("Maximum magnitude of the barotropic G column: %.3e\n",max(abs(GInverse(:,1))));

%% Use one coefficient vector for both physical variables
% A coefficient vector addresses the shared modal family. F synthesis sees
% its barotropic coefficient; G synthesis does not because G_0 vanishes.
% Accordingly F recovers the complete vector, while G recovers its active
% projection with the first row set to zero.
coefficients = [1;-0.7;0.45;-0.25;0.16;-0.10;0.06;-0.03];
FValues = transform.transformBack(coefficients,variable="F");
GValues = transform.transformBack(coefficients,variable="G");
FRecovered = transform.transformForward(FValues,variable="F");
GRecovered = transform.transformForward(GValues,variable="G");

expectedG = GProjector*coefficients;
fprintf("Relative F coefficient error: %.3e\n",norm(FRecovered-coefficients)/norm(coefficients));
fprintf("Relative active-G coefficient error: %.3e\n",norm(GRecovered-expectedG)/norm(expectedG));

figure(Name="V2 aligned F-G synthesis",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(FValues,z,"o-",LineWidth=1.3,MarkerSize=4)
grid on
xlabel("F")
ylabel("z (m)")
title("F sees the barotropic coefficient")

nexttile
plot(GValues,z,"o-",LineWidth=1.3,MarkerSize=4)
grid on
xlabel("G")
ylabel("z (m)")
title("G contains only active columns")

%% Inspect endpoint traces and the detached physical snapshot
% Endpoint traces are always ordered surface then bottom. The symbolic
% endpoint identities remain valid for any physical domain; their physical
% coordinates follow from `zDomain`. The transform also carries the
% physical state needed to interpret its matrices without retaining the
% source basis object.
FEndpoints = transform.endpointValues(variable="F");
GEndpoints = transform.endpointValues(variable="G");
endpointZ = [transform.zDomain(2);transform.zDomain(1)];
endpointSummary = table(transform.endpointLocations,endpointZ,FEndpoints(:,1),GEndpoints(:,1), ...
    VariableNames=["endpoint" "z" "barotropicF" "barotropicG"]);
disp(endpointSummary)

fprintf("Snapshot: depth %.1f m, g %.2f m s^-2, family %s, normalization %s.\n", ...
    transform.depth,transform.g,transform.modeFamily,transform.normalization);
fprintf("The assessment retained all %d explicitly requested aligned modes: %d.\n",nModes,assessment.retainedModeCount == nModes);

figure(Name="V2 aligned F-G mode columns",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(FInverse,z,LineWidth=1.0)
grid on
xlabel("F")
ylabel("z (m)")
title("Aligned F columns")

nexttile
plot(GInverse,z,LineWidth=1.0)
grid on
xlabel("G")
ylabel("z (m)")
title("Aligned G columns, including G_0=0")
