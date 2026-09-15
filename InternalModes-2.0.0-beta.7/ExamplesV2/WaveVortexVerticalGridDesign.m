%% Design a vertical grid for free-surface Wave-Vortex modes
% Given a stratification profile and a fixed number of vertical points,
% this example builds APV and mean-density-anomaly (MDA) transforms on one
% WKB-stretched Chebyshev--Lobatto quadrature rule. The two families share
% physical points and weights but choose their usable mode counts
% independently.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Specify the stratification, endpoints, and point count
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
Nz = 65;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);
g0 = -integral(N2,zDomain(1),zDomain(2));
gd = Inf;
gramTolerance = 1e-2;
quadraticAliasingTolerance = 0.1;

%% Solve the continuous APV and MDA mode families
nAvailableModes = Nz+4;
nEVP = max(96,3*nAvailableModes);
modeSolver = IMSolverSpectral(nEVP=nEVP);
apvEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd,surfaceBoundary="freeSurface");
mdaEVP = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd);
apvBasis = modeSolver.solveEVP(apvEVP,nModes=nAvailableModes);
mdaBasis = modeSolver.solveEVP(mdaEVP,nModes=nAvailableModes);

%% Construct one WKB-stretched quadrature and differentiation rule
% `nativeDifferentiationRule` returns the points, weights, and physical
% first-derivative matrix together. For this
% solver, the points are Chebyshev--Lobatto points in
% $$x(z)=\int_{-D}^{z}N(s)\,ds$$. The weights already act on values sampled
% in physical $$z$$, and `Dz*values` differentiates values in the same
% increasing-$$z$$ order. A final scalar adjustment makes the
% constant-function integral equal to the exact depth.
gridSolver = IMSolverSpectral(nEVP=Nz,coordinateKind="wkb").configuredForEVP(apvEVP);
[z,weights,Dz] = gridSolver.nativeDifferentiationRule(zDomain);
weights = weights*(D/sum(weights));

%% Select APV and MDA modes on the fixed rule
% Explicit weights are never replaced. Each call assesses all leading mode
% sets on this one rule and returns the largest set satisfying its errors.
% APV uses F/G Gram errors and the coupled quadratic-product error. MDA uses
% its G Gram error and may choose a different number of modes.
[apvTransform,apvAssessment] = apvBasis.discreteTransform(z=z,weights=weights,variables=["F","G"], ...
    gramTolerance=gramTolerance,quadraticAliasingTolerance=quadraticAliasingTolerance);
[mdaTransform,mdaAssessment] = mdaBasis.discreteTransform(z=z,weights=weights,variables="G",gramTolerance=gramTolerance);

apvRow = length(apvTransform.modeNumber);
mdaRow = length(mdaTransform.modeNumber);
family = ["APV";"MDA"];
modeCount = [apvRow;mdaRow];
lastModeNumber = [apvTransform.modeNumber(end);mdaTransform.modeNumber(end)];
gramError = [apvAssessment.prefixDiagnostics.gramError(apvRow);mdaAssessment.prefixDiagnostics.gramError(mdaRow)];
quadraticError = [apvAssessment.prefixDiagnostics.quadraticAliasingError(apvRow);NaN];
summary = table(family,modeCount,lastModeNumber,gramError,quadraticError);
disp(summary)
fprintf("Shared rule: %d points, minimum weight %.6g m, depth error %.3e.\n",Nz,min(weights),abs(sum(weights)-D)/D);
fprintf("Maximum linear-derivative error: %.3e.\n",max(abs(Dz*z-1)));

%% Inspect the quadrature and lowest modes
% Physical mode labels control line colors in both family panels. Negative
% labels use red shades, the zero mode is gray, and equal positive labels
% use the same color in APV and MDA.
zPlot = linspace(zDomain(1),zDomain(2),801).';
nModeShapes = 4;
apvShapeCount = min(nModeShapes,length(apvTransform.modeNumber));
mdaShapeCount = min(nModeShapes,length(mdaTransform.modeNumber));
apvGReference = apvBasis.G(zPlot);
apvGReference = apvGReference(:,1:apvShapeCount);
mdaGReference = mdaBasis.G(zPlot);
mdaGReference = mdaGReference(:,1:mdaShapeCount);
apvGPlot = normalizeColumns(apvGReference);
apvGPoints = apvBasis.G(z);
apvGPoints = normalizeWithReference(apvGPoints(:,1:apvShapeCount),apvGReference);
mdaGPlot = normalizeColumns(mdaGReference);
mdaGPoints = mdaBasis.G(z);
mdaGPoints = normalizeWithReference(mdaGPoints(:,1:mdaShapeCount),mdaGReference);

figure(Name="Free-surface Wave-Vortex vertical quadrature",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");

nexttile
plot(sqrt(N2(zPlot)),zPlot,"k-",LineWidth=1.4)
hold on
plot(sqrt(N2(z)),z,"o",MarkerSize=4)
hold off
grid on
xlabel("N (s^{-1})")
ylabel("z (m)")
title("WKB-stretched points")

nexttile
plot(weights,z,"o-",LineWidth=1.1,MarkerSize=4)
grid on
xlabel("physical weight (m)")
ylabel("z (m)")
title("Shared APV/MDA weights")

nexttile
plotModeFamily(apvGPlot,apvGPoints,zPlot,z,apvTransform.modeNumber(1:apvShapeCount),"APV")

nexttile
plotModeFamily(mdaGPlot,mdaGPoints,zPlot,z,mdaTransform.modeNumber(1:mdaShapeCount),"MDA")

function values = normalizeColumns(values)
scale = max(abs(values),[],1);
values = values./scale;
end

function values = normalizeWithReference(values,reference)
scale = max(abs(reference),[],1);
values = values./scale;
end

function plotModeFamily(modeValues,pointValues,zPlot,z,modeNumbers,family)
hold on
labels = strings(1,length(modeNumbers));
for iMode = 1:length(modeNumbers)
    color = colorForModeNumber(modeNumbers(iMode));
    plot(modeValues(:,iMode),zPlot,Color=color,LineWidth=1.3)
    plot(pointValues(:,iMode),z,"o",Color=color,MarkerSize=3)
    labels(iMode) = sprintf("%s %g",family,modeNumbers(iMode));
end
hold off
grid on
xlabel("normalized G")
ylabel("z (m)")
title(family+" modes on the shared rule")
legend(labels,Location="best")
end

function color = colorForModeNumber(modeNumber)
positiveColors = lines(7);
negativeColors = [0.64 0.08 0.18;0.85 0.33 0.10];
if modeNumber < 0
    color = negativeColors(min(abs(modeNumber),size(negativeColors,1)),:);
elseif modeNumber == 0
    color = [0.35 0.35 0.35];
else
    color = positiveColors(mod(modeNumber-1,size(positiveColors,1))+1,:);
end
end
