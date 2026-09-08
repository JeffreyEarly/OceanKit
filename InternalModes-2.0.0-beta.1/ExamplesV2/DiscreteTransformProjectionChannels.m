%% Distinguish modal pairings, coefficients, and diagnostic channels
% A sampled inner product first produces modal pairings. Extracting modal
% coefficients requires solving the sampled Gram system. This distinction
% matters whenever the quadrature rule is not exactly orthogonal. A second
% example then shows how an internal-mode variable can remain available for
% synthesis and endpoint diagnostics even when it has no direct forward
% projection.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Build an intentionally imperfect exponential-stratification rule
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nModes = 6;
evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain,g=g);
solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=10);
basisSet.normalization = "geostrophic";

nPoints = 10;
sigma = linspace(0,1,nPoints).';
z = zDomain(1)+D*(1-(1-sigma).^2);
edges = [zDomain(1);0.5*(z(1:end-1)+z(2:end));zDomain(2)];
geometricWeights = diff(edges);
transform = basisSet.discreteTransform(z=z,weights=geometricWeights,nModes=nModes,variables=["F","G"],gramTolerance=1);

%% Compare the projection functional with coefficient extraction
% For sampled values X in variable V, `modeProjectionFunctional` returns
%
%   P_V(X) = (A_i^V)' W_V X.
%
% If X=A_i^V a, these pairings equal Gamma_V a rather than a itself. The
% forward transform solves the active Gram system and recovers a.
coefficients = [0;0.8;-0.55;0.30;-0.16;0.07];
GValues = transform.transformBack(coefficients,variable="G");
GPairings = transform.modeProjectionFunctional(GValues,variable="G");
GRecovered = transform.transformForward(GValues,variable="G");
GGramPrediction = transform.gramMatrix(variable="G")*coefficients;

comparison = table(transform.modeNumber.',coefficients,GPairings,GGramPrediction,GRecovered, ...
    VariableNames=["modeNumber" "coefficient" "pairing" "gramTimesCoefficient" "recoveredCoefficient"]);
disp(comparison)

fprintf("Pairing identity error: %.3e\n",norm(GPairings-GGramPrediction));
fprintf("Coefficient recovery error: %.3e\n",norm(GRecovered-coefficients));
fprintf("Distance between raw pairings and coefficients: %.3e\n",norm(GPairings-coefficients));

figure(Name="V2 projection functionals and coefficients",Color="w");
stem(transform.modeNumber,coefficients,"filled",LineWidth=1.1)
hold on
stem(transform.modeNumber,GPairings,"x",LineWidth=1.1)
stem(transform.modeNumber,GRecovered,"o",LineWidth=1.1)
hold off
grid on
xlabel("physical mode number")
ylabel("modal value")
title("Pairings require a Gram solve to become coefficients")
legend(["true coefficients" "raw pairings" "forward-transform coefficients"],Location="best")

%% Keep a companion wave variable for synthesis without projecting it
% For fixed-frequency wave modes, this EVP catalog supplies a direct G
% inner product but no continuous F inner product. Requesting G therefore
% creates a one-channel forward transform. The sampled F modes and their
% endpoint traces remain part of the aligned physical snapshot.
waveN2 = @(z) 1e-4*ones(size(z));
waveEVP = IMInternalModes.waveModesAtFrequency(N2=waveN2,zDomain=zDomain,omega=2e-3);
waveSolver = IMSolverSpectral(nEVP=128);
waveBasis = waveSolver.solveEVP(waveEVP,nModes=4);
waveTransform = waveBasis.discreteTransform(z=z,weights=geometricWeights,nModes=3,variables="G",gramTolerance=1);

channel = ["F";"G"];
hasForwardTransform = [waveTransform.hasForwardTransform(variable="F");waveTransform.hasForwardTransform(variable="G")];
reason = [waveTransform.forwardTransformReason(variable="F");waveTransform.forwardTransformReason(variable="G")];
channelAvailability = table(channel,hasForwardTransform,reason);
disp(channelAvailability)

waveCoefficients = [1;-0.4;0.2];
diagnosticF = waveTransform.transformBack(waveCoefficients,variable="F");
projectedG = waveTransform.transformBack(waveCoefficients,variable="G");
FEndpoints = waveTransform.endpointValues(variable="F");

fprintf("Diagnostic F synthesis matrix size: %d x %d.\n",size(waveTransform.inverseMatrix(variable="F")));
fprintf("F endpoint trace matrix size: %d x %d, ordered surface then bottom.\n",size(FEndpoints));

figure(Name="V2 direct and diagnostic transform channels",Color="w");
plot(diagnosticF,z,"-",LineWidth=1.4)
hold on
plot(projectedG,z,"--",LineWidth=1.4)
hold off
grid on
xlabel("normalized amplitude")
ylabel("z (m)")
title("Both variables synthesize; only G projects directly")
legend(["diagnostic F synthesis" "direct G channel"],Location="best")
