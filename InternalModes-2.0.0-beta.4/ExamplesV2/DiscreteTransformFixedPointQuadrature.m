%% Fit aligned F/G transforms on caller-chosen physical points
% Suppose measurements already exist on a fixed vertical grid. This
% walkthrough fits one quadrature rule on those points, builds aligned F/G
% Galerkin transforms, and separates the diagnostics that describe shared
% quadrature accuracy from those that describe each variable's linear
% algebra.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve paired exponential-stratification modes
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 12;
nModes = 8;
nEVP = 160;
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Choose a surface-refined grid and fit its quadrature weights
% The points include both boundaries and become more closely spaced toward
% the surface. Supplying `z` and `nModes` makes the requested eight-mode
% band strict: every requested mode must pass the enabled policies.
nPoints = 24;
sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

variables = ["F","G"];
[transform, assessment] = basisSet.discreteTransform(z=z,nModes=nModes,variables=variables);
fit = assessment.weightFit;

% Supplying weights bypasses fitting. Use the geometric control-volume
% weights stored by the fit to build a comparison on exactly the same grid.
geometricTransform = fit.geometricTransform;

% Every normalized-Gram objective row retains its variable and mode-pair
% provenance. F has one active row for every upper-triangle pair. G omits
% any identically zero family column while preserving its aligned label.
objectiveVariable = variables.';
objectiveRows = [sum(fit.objectiveRowVariables == "F");sum(fit.objectiveRowVariables == "G")];
fittedResidual = [fit.variableResidualNorm(variable="F");fit.variableResidualNorm(variable="G")];
geometricResidual = [fit.variableGeometricResidualNorm(variable="F");fit.variableGeometricResidualNorm(variable="G")];
stackedObjectiveSummary = table(objectiveVariable,objectiveRows,fittedResidual,geometricResidual);
disp(stackedObjectiveSummary)

%% Separate fit, transform, and conditioning diagnostics
% The fitted objective is the Frobenius norm of the normalized Gram error;
% it aggregates errors over all mode pairs. `relativeGramOperatorError` is
% the worst error over normalized modal combinations. `roundTripError`
% measures algebraic coefficient recovery and may remain tiny even when the
% quadrature rule has visible Gram error.
rule = ["fitted";"fitted";"geometric";"geometric"];
variable = ["F";"G";"F";"G"];
frobeniusResidual = [fit.variableResidualNorm(variable="F");fit.variableResidualNorm(variable="G"); ...
    fit.variableGeometricResidualNorm(variable="F");fit.variableGeometricResidualNorm(variable="G")];
transforms = {transform;transform;geometricTransform;geometricTransform};
relativeGramOperatorError = zeros(4,1);
roundTripError = zeros(4,1);
forwardMatrixConditionNumber = zeros(4,1);
inverseMatrixConditionNumber = zeros(4,1);
gramConditionNumber = zeros(4,1);
for iRow = 1:4
    current = transforms{iRow};
    currentVariable = variable(iRow);
    relativeGramOperatorError(iRow) = current.relativeGramOperatorError(variable=currentVariable);
    roundTripError(iRow) = current.roundTripError(variable=currentVariable);
    forwardMatrixConditionNumber(iRow) = cond(current.forwardMatrix(variable=currentVariable));
    inverseMatrixConditionNumber(iRow) = current.inverseMatrixConditionNumber(variable=currentVariable);
    gramConditionNumber(iRow) = current.gramConditionNumber(variable=currentVariable);
end
depthSum = [sum(transform.weights);sum(transform.weights);sum(geometricTransform.weights);sum(geometricTransform.weights)];
minimumWeight = [min(transform.weights);min(transform.weights);min(geometricTransform.weights);min(geometricTransform.weights)];
maximumWeight = [max(transform.weights);max(transform.weights);max(geometricTransform.weights);max(geometricTransform.weights)];
diagnostics = table(rule,variable,frobeniusResidual,relativeGramOperatorError,roundTripError, ...
    forwardMatrixConditionNumber,inverseMatrixConditionNumber,gramConditionNumber, ...
    depthSum,minimumWeight,maximumWeight);

fprintf("\nFixed-point transform for aligned exponential hydrostatic F/G modes\n");
fprintf("Requested modes: %d; fixed points: %d; normalization: %s\n\n",nModes,nPoints,transform.normalization);
disp(diagnostics)

% Prefix diagnostics show how the same fitted rule behaves as modes are
% added. The Gram policy is the acceptance test; the other columns diagnose
% algebraic consistency and sensitivity.
prefixDiagnostics = assessment.prefixDiagnostics(:,["modeCount" "lastModeNumber" "gramError" "gramLimitingVariable" "gramAccepted"]);
disp(prefixDiagnostics)
fprintf("\nPer-variable prefix diagnostics for G:\n");
disp(assessment.variablePrefixDiagnostics(variable="G"))

%% Recover one coefficient matrix through both physical variables
% The same coefficient matrix addresses the aligned F/G family. A field
% made entirely from retained modes should round-trip through either direct
% channel to numerical precision.
mode = (1:nModes).';
coefficientsTrue = [1./mode (-1).^(mode - 1)./mode.^2];
FValuesFromModes = transform.transformBack(coefficientsTrue,variable="F");
GValuesFromModes = transform.transformBack(coefficientsTrue,variable="G");
FCoefficientsRecovered = transform.transformForward(FValuesFromModes,variable="F");
GCoefficientsRecovered = transform.transformForward(GValuesFromModes,variable="G");
variable = ["F";"G"];
coefficientRoundTripError = [norm(FCoefficientsRecovered-coefficientsTrue,2);norm(GCoefficientsRecovered-coefficientsTrue,2)]/norm(coefficientsTrue,2);
disp(table(variable,coefficientRoundTripError))

%% Project a smooth sampled profile
% Galerkin projection returns the retained modal field whose residual is
% orthogonal to the retained basis in the sampled metric W. An arbitrary
% profile is therefore projected, not expected to round-trip exactly.
profile = sin(pi*(z - zDomain(1))/D).*exp(z/900).*(1 + 0.2*cos(2*pi*(z - zDomain(1))/D));
profileCoefficients = transform.transformForward(profile,variable="G");
profileReconstruction = transform.transformBack(profileCoefficients,variable="G");
profileResidual = profile - profileReconstruction;
GMetric = transform.metricMatrix(variable="G");
profileNorm = sqrt(profile.'*GMetric*profile);
relativeProfileResidual = sqrt(profileResidual.'*GMetric*profileResidual)/profileNorm;
fprintf("Relative sampled-metric residual of the smooth profile: %.3e\n",relativeProfileResidual);

%% Inspect the sampled modes, weights, and Gram errors
targetGramG = transform.targetGramMatrix(variable="G");
sampledGramG = transform.gramMatrix(variable="G");
geometricTargetGramG = geometricTransform.targetGramMatrix(variable="G");
geometricSampledGramG = geometricTransform.gramMatrix(variable="G");
targetNorms = diag(targetGramG);
gramScale = 1./sqrt(abs(targetNorms));
fittedGramError = gramScale.*(sampledGramG-targetGramG).*gramScale.';
geometricGramError = gramScale.*(geometricSampledGramG-geometricTargetGramG).*gramScale.';
gramColorLimit = max(abs([fittedGramError geometricGramError]),[],"all");

figure(Name="V2 fixed-point aligned F-G discrete transform",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");

nexttile
plot(transform.inverseMatrix(variable="G"),z,LineWidth=1.1)
grid on
xlabel("G")
ylabel("z (m)")
title("Retained sampled modes")

nexttile
plot(fit.weights,z,"o-",LineWidth=1.1,MarkerSize=4)
hold on
plot(fit.geometricWeights,z,".-",LineWidth=1.1,MarkerSize=10)
hold off
grid on
xlabel("weight (m)")
ylabel("z (m)")
title("Quadrature weights")
legend(["fitted" "geometric"],Location="best")

nexttile
imagesc(transform.modeNumber,transform.modeNumber,geometricGramError)
axis xy image
clim([-gramColorLimit gramColorLimit])
colorbar
xlabel("mode number")
ylabel("mode number")
title("Geometric normalized Gram error")

nexttile
imagesc(transform.modeNumber,transform.modeNumber,fittedGramError)
axis xy image
clim([-gramColorLimit gramColorLimit])
colorbar
xlabel("mode number")
ylabel("mode number")
title("Fitted normalized Gram error")

figure(Name="V2 fixed-point profile projection",Color="w");
plot(profile,z,"k-",LineWidth=1.5)
hold on
plot(profileReconstruction,z,"--",LineWidth=1.5)
plot(profileResidual,z,":",LineWidth=1.2)
hold off
grid on
xlabel("amplitude")
ylabel("z (m)")
title("Smooth profile and retained-mode projection")
legend(["sampled profile" "reconstruction" "residual"],Location="best")
