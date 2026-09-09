%% Build aligned internal-mode transforms from an exact point budget
% Suppose an instrument or numerical model fixes the number of available
% vertical samples. This walkthrough asks how many aligned hydrostatic F/G
% modes an exact physical point budget can support and lets the transform
% assessment choose one common production family band.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve rigid-boundary exponential-stratification G modes
% This is the scalar G formulation of the hydrostatic internal-mode EVP.
% Both endpoints use Dirichlet conditions, so every mode satisfies
%
%   G(-D) = G(0) = 0.
%
% The modes are solved spectrally in the WKB coordinate and then placed in
% geostrophic normalization. A profile reconstructed from this basis should
% therefore obey the same zero-value boundary conditions.
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 24;
nEVP = 192;
surfaceBoundary = IMBoundaryCondition.dirichlet();
bottomBoundary = IMBoundaryCondition.dirichlet();
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=surfaceBoundary,bottomBoundary=bottomBoundary);
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";
fprintf("Mode basis: hydrostatic G modes with G(-D) = G(0) = 0; normalization: %s.\n",basisSet.normalization);

%% Ask for an exact mode-root point count
% `nPoints` is exact, not an upper bound. The builder searches the available
% mode-root grids, finds the largest candidate band with exactly this many
% physical points, fits one weight vector for that full band, and assesses
% every leading modal prefix without refitting the weights.
nPoints = 18;
variables = ["F","G"];
[transform, assessment] = basisSet.discreteTransform(nPoints=nPoints,variables=variables);

requestedPointCount = assessment.requestedPointCount;
actualPointCount = assessment.actualPointCount;
candidateModeCount = assessment.candidateModeCount;
retainedModeCount = assessment.retainedModeCount;
gramTolerance = assessment.gramPolicy.tolerance;
limitingPolicy = assessment.limitingPolicy;
transformSummary = table(requestedPointCount,actualPointCount,candidateModeCount,retainedModeCount,gramTolerance,limitingPolicy);
disp(transformSummary)
fprintf("Directly projected channels: %s. Both channels remain aligned to the same family labels.\n",join(transform.availableVariables,", "));

% `weightFit.transform` and `candidateTransform` contain the complete band
% used to fit the rule. `transform` and `assessment.transform` contain only
% the production prefix accepted by the active policies. They use exactly
% the same physical points and weights.
candidateTransform = assessment.weightFit.transform;
fprintf("Candidate physical mode labels: %s\n",mat2str(candidateTransform.modeNumber));
fprintf("Retained physical mode labels:  %s\n",mat2str(transform.modeNumber));
fprintf("Candidate and production rules share points and weights: %d\n", ...
    isequal(candidateTransform.z,transform.z) && isequal(candidateTransform.weights,transform.weights));

prefixDiagnostics = assessment.prefixDiagnostics(:,["modeCount" "lastModeNumber" "gramError" ...
    "gramLimitingVariable" "gramAccepted" "combinedAccepted"]);
disp(prefixDiagnostics)
disp(assessment.variablePrefixDiagnostics(variable="G"))

%% Recover known modal coefficients
% A field assembled from the retained modes should transform back to sample
% space and forward to coefficients with only numerical round-off error.
mode = (1:retainedModeCount).';
coefficientsTrue = [1./mode (-1).^(mode - 1)./mode.^2];
valuesFromModes = transform.transformBack(coefficientsTrue,variable="G");
coefficientsRecovered = transform.transformForward(valuesFromModes,variable="G");
coefficientRoundTripError = norm(coefficientsRecovered - coefficientsTrue,2)/norm(coefficientsTrue,2);
fprintf("Relative coefficient round-trip error: %.3e\n",coefficientRoundTripError);

%% Project a smooth profile onto the retained band
% The test profile includes exponential surface intensification but vanishes
% at both boundaries, matching the G-mode endpoint conditions. With 18
% points, the two zero endpoints leave 16 interior samples for 16 retained
% modes. The reconstruction is consequently almost interpolatory on the
% transform grid, so accuracy is assessed on an independent fine grid.
z = transform.z;
x = (z - zDomain(1))/D;
profile = sin(pi*x).*exp(z/1800).*(1 + 0.15*cos(2*pi*x));
profileCoefficients = transform.transformForward(profile,variable="G");
profileReconstruction = transform.transformBack(profileCoefficients,variable="G");
profileResidual = profile - profileReconstruction;
GMetric = transform.metricMatrix(variable="G");
profileNorm = sqrt(profile.'*GMetric*profile);
relativeSampledResidual = sqrt(profileResidual.'*GMetric*profileResidual)/profileNorm;

zFine = linspace(zDomain(1),zDomain(2),801).';
xFine = (zFine - zDomain(1))/D;
profileFine = sin(pi*xFine).*exp(zFine/1800).*(1 + 0.15*cos(2*pi*xFine));
GFine = basisSet.G(zFine);
profileReconstructionFine = GFine(:,1:retainedModeCount)*profileCoefficients;
profileResidualFine = profileFine - profileReconstructionFine;
continuousWeight = N2(zFine)/g;
relativeFineGridResidual = sqrt(trapz(zFine,continuousWeight.*profileResidualFine.^2) ...
    /trapz(zFine,continuousWeight.*profileFine.^2));
fprintf("Relative residual on the transform points: %.3e\n",relativeSampledResidual);
fprintf("Relative N2/g-weighted residual on a fine grid: %.3e\n",relativeFineGridResidual);

%% Inspect the point-limited rule
figure(Name="V2 point-limited scalar discrete transform",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(transform.inverseMatrix(variable="G"),z,LineWidth=1.1)
grid on
xlabel("G")
ylabel("z (m)")
title(sprintf("%d Dirichlet G modes on %d points",retainedModeCount,nPoints))

nexttile
plot(transform.weights,z,"o-",LineWidth=1.1,MarkerSize=4)
grid on
xlabel("weight (m)")
ylabel("z (m)")
title("Fitted mode-root quadrature")

figure(Name="V2 point-limited profile projection",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(profileFine,zFine,"k-",LineWidth=1.5)
hold on
plot(profileReconstructionFine,zFine,"--",LineWidth=1.5)
plot(profile,z,"o",MarkerSize=4)
hold off
grid on
xlabel("amplitude")
ylabel("z (m)")
title("Boundary-compatible G-mode reconstruction")
legend(["target profile" "fine-grid reconstruction" "transform points"],Location="best")

nexttile
plot(profileResidualFine,zFine,LineWidth=1.3)
grid on
xlabel("target - reconstruction")
ylabel("z (m)")
title("Error between transform points")

%% Compare several exact point budgets
% Not every integer must occur as a mode-root point count. An unattainable
% request reports nearby attainable counts. This sweep deliberately uses
% attainable budgets so the tutorial remains runnable from start to finish.
pointBudgets = [6 10 14 18 22].';
candidateCounts = zeros(size(pointBudgets));
retainedCounts = zeros(size(pointBudgets));
maximumGramError = zeros(size(pointBudgets));
for iBudget = 1:length(pointBudgets)
    [~, budgetAssessment] = basisSet.discreteTransform(nPoints=pointBudgets(iBudget),variables=variables);
    candidateCounts(iBudget) = budgetAssessment.candidateModeCount;
    retainedCounts(iBudget) = budgetAssessment.retainedModeCount;
    maximumGramError(iBudget) = max(budgetAssessment.prefixDiagnostics.gramError);
end
budgetSummary = table(pointBudgets,candidateCounts,retainedCounts,maximumGramError);
disp(budgetSummary)

figure(Name="V2 exact point budgets",Color="w");
plot(pointBudgets,candidateCounts,"o-",LineWidth=1.3,MarkerSize=6)
hold on
plot(pointBudgets,retainedCounts,"s--",LineWidth=1.3,MarkerSize=6)
hold off
grid on
xlabel("exact physical point count")
ylabel("mode count")
title("Candidate and retained bands grow with the point budget")
legend(["candidate" "retained"],Location="northwest")
