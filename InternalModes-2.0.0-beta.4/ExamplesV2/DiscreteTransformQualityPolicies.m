%% Select an aligned F/G transform band with quality policies
% A fixed point count can represent a candidate modal band without making
% every prefix equally useful for an application. This walkthrough applies
% Gram, rejected-mode leakage, and coupled quadratic-aliasing policies to
% one exponential-stratification mode-root rule.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve a basis wider than the transform candidate band
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 24;
nEVP = 192;
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Apply three illustrative quality budgets
% These tolerances make all three limits visible for this example. They are
% pedagogical application budgets, not universal defaults or recommendations.
% The weights are fitted once for the complete candidate band; every prefix
% below is assessed on those same points and weights without refitting.
nPoints = 18;
nCheckModes = 20;
gramTolerance = 1e-2;
leakageTolerance = 4e-2;
quadraticAliasingTolerance = 0.1;
variables = ["F","G"];
[transform, assessment] = basisSet.discreteTransform(nPoints=nPoints,gramTolerance=gramTolerance, ...
    leakageTolerance=leakageTolerance,nCheckModes=nCheckModes, ...
    quadraticAliasingTolerance=quadraticAliasingTolerance,variables=variables);

fprintf("\nCandidate modes: %d; retained modes: %d; limiting policy: %s\n", ...
    assessment.candidateModeCount,assessment.retainedModeCount,assessment.limitingPolicy);
fprintf("Candidate physical mode labels: %s\n",mat2str(assessment.candidateModeNumber));
fprintf("Retained physical mode labels:  %s\n\n",mat2str(assessment.retainedModeNumber));

%% Read each policy result
% The Gram error is the worst normalized Gram distortion over F and G.
% Leakage projects rejected F modes into retained F and rejected G modes
% into retained G, then reports the worse channel through `nCheckModes`.
% Coupled quadratic aliasing compares independently integrated continuous
% projections in three physical product channels:
%
%   F_i F_j -> F,   F_i G_j -> G,   G_i G_j -> F.
%
% `DiscreteTransformCoupledProducts.m` shows these channels and their first
% exact DCT-I/DST-I alias explicitly. Here the same definition is applied
% to a nonuniform exponential-stratification family.
%
% This scalar metric over modal products is more specific than projecting
% one arbitrary smooth profile through the discrete transform.
policyName = ["worst-channel Gram";"worst-channel rejected-mode leakage";"coupled quadratic aliasing"];
tolerance = [assessment.gramPolicy.tolerance; assessment.leakagePolicy.tolerance; ...
    assessment.quadraticAliasingPolicy.tolerance];
maximumAcceptedModeCount = [assessment.gramPolicy.maximumAcceptedModeCount; ...
    assessment.leakagePolicy.maximumAcceptedModeCount; ...
    assessment.quadraticAliasingPolicy.maximumAcceptedModeCount];
limitingValue = [assessment.gramPolicy.limitingValue; assessment.leakagePolicy.limitingValue; ...
    assessment.quadraticAliasingPolicy.limitingValue];
reason = [assessment.gramPolicy.reason; assessment.leakagePolicy.reason; assessment.quadraticAliasingPolicy.reason];
policySummary = table(policyName,tolerance,maximumAcceptedModeCount,limitingValue,reason);
disp(policySummary)

policyDiagnostics = assessment.prefixDiagnostics(:,["modeCount" "lastModeNumber" "gramError" ...
    "gramLimitingVariable" "leakageError" "leakageLimitingVariable" "quadraticAliasingError" "quadraticLimitingChannel" "gramAccepted" "leakageAccepted" ...
    "quadraticAccepted" "combinedAccepted"]);
disp(policyDiagnostics)

firstQuadraticRejection = find(~assessment.prefixDiagnostics.quadraticAccepted,1);
if ~isempty(firstQuadraticRejection)
    rejected = assessment.prefixDiagnostics(firstQuadraticRejection,:);
    fprintf("First rejected quadratic product: %s from physical modes %g and %g.\n", ...
        rejected.quadraticLimitingChannel,rejected.quadraticLimitingModeNumberI,rejected.quadraticLimitingModeNumberJ);
end

% The combined rule retains the intersection of all cumulative decisions.
% The common family band is the minimum accepted prefix over every enabled
% policy and variable. The printed table makes the calibrated limits and
% the channel responsible for each one explicit for this solver setup.
candidateTransform = assessment.weightFit.transform;
sameFixedRule = isequal(candidateTransform.z,transform.z) && isequal(candidateTransform.weights,transform.weights);
fprintf("Candidate and production prefixes use the same fixed rule: %d\n",sameFixedRule);
fprintf("Leakage checks rejected source modes through physical mode %g.\n", ...
    basisSet.modeNumber(assessment.leakagePolicy.nCheckModes));

%% Inspect where each policy begins to limit the band
diagnostics = assessment.prefixDiagnostics;
modeCount = diagnostics.modeCount;

figure(Name="V2 aligned F-G transform quality policies",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
semilogy(modeCount,max(diagnostics.gramError,eps),"o-",LineWidth=1.2,MarkerSize=4)
hold on
yline(gramTolerance,"--",LineWidth=1.0)
xline(assessment.gramPolicy.maximumAcceptedModeCount,":",LineWidth=1.2)
hold off
grid on
xlabel("candidate prefix count")
ylabel("normalized Gram error")
title(sprintf("Gram accepts %d",assessment.gramPolicy.maximumAcceptedModeCount))

nexttile
semilogy(modeCount,max(diagnostics.leakageError,eps),"o-",LineWidth=1.2,MarkerSize=4)
hold on
yline(leakageTolerance,"--",LineWidth=1.0)
xline(assessment.leakagePolicy.maximumAcceptedModeCount,":",LineWidth=1.2)
hold off
grid on
xlabel("candidate prefix count")
ylabel("rejected-mode leakage")
title(sprintf("Leakage accepts %d",assessment.leakagePolicy.maximumAcceptedModeCount))

nexttile
semilogy(modeCount,max(diagnostics.quadraticAliasingError,eps),"o-",LineWidth=1.2,MarkerSize=4)
hold on
yline(quadraticAliasingTolerance,"--",LineWidth=1.0)
xline(assessment.quadraticAliasingPolicy.maximumAcceptedModeCount,":",LineWidth=1.2)
hold off
grid on
xlabel("candidate prefix count")
ylabel("quadratic-aliasing error")
title(sprintf("Quadratic policy accepts %d",assessment.quadraticAliasingPolicy.maximumAcceptedModeCount))

%% Understand automatic and strict selection
% Omitting `nModes`, as above, returns the largest prefix accepted by every
% enabled policy. With explicit points, an explicit `nModes` is strict and
% throws `IMBasisSet:StrictDiscreteModeCountRejected` instead of silently
% reducing the request. For example, the following commented call asks for
% the first failing quadratic prefix on this already fitted rule:
%
% basisSet.discreteTransform(z=candidateTransform.z,weights=candidateTransform.weights,nModes=11, ...
%     gramTolerance=gramTolerance,leakageTolerance=leakageTolerance,nCheckModes=nCheckModes, ...
%     quadraticAliasingTolerance=quadraticAliasingTolerance,variables=variables);
%
% Leakage and quadratic aliasing use norm ratios and therefore require a
% positive-definite target Gram matrix. Signed or indefinite target metrics
% remain supported by the Gram policy alone.
