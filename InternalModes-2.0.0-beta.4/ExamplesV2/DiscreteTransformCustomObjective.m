%% Investigation: custom quadrature objectives on fixed physical points
% The default fixed-point fit minimizes a normalized Gram residual. This
% walkthrough appends a small penalty on relative movement away from the
% geometric weights, illustrating how custom objectives can trade a small
% amount of Parseval accuracy for a less extreme quadrature rule.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Build the aligned hydrostatic F/G mode basis
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 12;
nModes = 8;
nPoints = 24;
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

%% Compute the default normalized Gram Frobenius fit
variables = ["F","G"];
[~, defaultFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables=variables);

%% Add relative regularization toward geometric weights
% A custom callback receives the default normalized Gram system in
% `context.normalizedGramA` and `context.normalizedGramB`. This is the
% stacked F/G system; `context.normalizedGramVariables(q)` and
% `context.normalizedGramModePairs(q,:)` identify the channel and retained
% pair for row `q`. Scalar-style channel contexts are also available as
% `context.variableContexts.F` and `.G`. Returning A and b replaces the
% objective, so callers may reweight rows or append residuals. Here the
% additional rows penalize
%
%   sqrt(lambda) * (w_i/wGeometric_i - 1).
lambda = 1e-6;
objective = @(context) regularizedGramObjective(context,lambda);
[~, regularizedFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables=variables,objective=objective);

% Both fits retain the physical defaults `nonnegative=true` and
% `constrainDepth=true`. Algebraic investigations may disable them explicitly:
% algebraicWeights = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,nonnegative=false);
% freeDepthWeights = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,constrainDepth=false);

%% Compare every rule under the same Gram metric
% The custom objective contains extra regularization rows, so compare all
% candidates using the default normalized-Gram system rather than comparing
% their different full objective norms.
rule = ["geometric"; "normalized Gram"; "regularized Gram"];
weights = [defaultFit.geometricWeights defaultFit.weights regularizedFit.weights];
gramResiduals = defaultFit.objectiveMatrix*weights - defaultFit.objectiveTarget;
normalizedGramResidual = vecnorm(gramResiduals,2,1).';
relativeGramOperatorError = [worstGramError(defaultFit.geometricTransform,variables); ...
    worstGramError(defaultFit.transform,variables);worstGramError(regularizedFit.transform,variables)];
relativeWeightDisplacement = vecnorm((weights - defaultFit.geometricWeights)./defaultFit.geometricWeights,2,1).'/sqrt(nPoints);
roundTripError = [worstRoundTripError(defaultFit.geometricTransform,variables); ...
    worstRoundTripError(defaultFit.transform,variables);worstRoundTripError(regularizedFit.transform,variables)];
depthSum = sum(weights,1).';
minimumWeight = min(weights,[],1).';
maximumWeight = max(weights,[],1).';
diagnostics = table(rule,normalizedGramResidual,relativeGramOperatorError,relativeWeightDisplacement, ...
    roundTripError,depthSum,minimumWeight,maximumWeight);

fprintf("\nCustom fixed-point quadrature objective\n");
fprintf("Relative geometric-weight regularization: lambda = %.1e\n\n",lambda);
disp(diagnostics)

%% Inspect the Parseval-versus-regularity tradeoff
plotLabels = ["geom" "Gram" "Gram + reg."];
figure(Name="V2 custom fixed-point quadrature objective",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
plot(weights(:,1),z,".-",LineWidth=1.1,MarkerSize=10)
hold on
plot(weights(:,2),z,"o-",LineWidth=1.1,MarkerSize=4)
plot(weights(:,3),z,"s-",LineWidth=1.1,MarkerSize=4)
hold off
grid on
xlabel("weight (m)")
ylabel("z (m)")
title("Quadrature weights")
legend(rule,Location="best")

nexttile
bar(relativeGramOperatorError)
set(gca,YScale="log",XTick=1:3,XTickLabel=plotLabels,XTickLabelRotation=20)
grid on
ylabel("relative Gram operator error")
title("Parseval accuracy")

nexttile
bar(max(relativeWeightDisplacement,eps))
set(gca,YScale="log",XTick=1:3,XTickLabel=plotLabels,XTickLabelRotation=20)
grid on
ylabel("RMS relative displacement")
title("Distance from geometric weights")

%% Local objective builder
function specification = regularizedGramObjective(context,lambda)
relativeScale = 1./context.geometricWeights;
regularizationMatrix = sqrt(lambda)*diag(relativeScale);
regularizationTarget = sqrt(lambda)*ones(length(context.z),1);
specification = struct("A",[context.normalizedGramA; regularizationMatrix], ...
    "b",[context.normalizedGramB; regularizationTarget], ...
    "name","normalizedGramWithGeometricRegularization");
end

function value = worstGramError(transform,variables)
value = max(arrayfun(@(variable) transform.relativeGramOperatorError(variable=variable),variables));
end

function value = worstRoundTripError(transform,variables)
value = max(arrayfun(@(variable) transform.roundTripError(variable=variable),variables));
end
