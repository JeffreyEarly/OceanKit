%% Project partial-depth observations onto resolvable vertical modes
% Observations are an observation-specific inverse problem. The modal basis
% and Parseval normalization are fixed before looking at the data. The
% observation operator decides which coefficients can be recovered stably.
%
% This example uses hydrostatic G modes, configurable partial-depth sampling
% scenarios, Gaussian modal-coefficient realizations, and an equal-weight
% observation metric. QR with column pivoting is used only to identify a
% well-conditioned subset of modes; the final projection is built from the
% physically normalized modal columns.
%
% For these G modes,
%
%     int N2(z) eta(z)^2 dz = g sum_j eta_j^2,
%
% so g eta_j^2 is the potential-energy spectrum. The observation operator
% changes what part of that spectrum is recoverable.

%% Choose sampling and realization scenarios
% Edit these values, or define them before running the script.
if ~exist('samplingScenario','var')
    samplingScenario = "irregularPartialDepth";
end
if ~exist('realizationScenario','var')
    realizationScenario = "gaussianSeedA";
end
if ~exist('jStar','var')
    jStar = 3;
end
if ~exist('nObs','var')
    nObs = 26;
end

%% Build a candidate hydrostatic G-mode basis
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
omega = 0;
nCandidateModes = 48;
nEVP = max(256,ceil(2.1*(nCandidateModes + 1)));
zDomain = [-Lz 0];
zPlot = linspace(zDomain(1),zDomain(2),512).';

imPlot = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zPlot,latitude=latitude,nEVP=nEVP,nModes=nCandidateModes,g=g);
imPlot.normalization = Normalization.geostrophic;
imPlot.upperBoundary = UpperBoundary.rigidLid;
[~,GPlot] = imPlot.modesAtFrequency(omega);
PhiPlot = GPlot(:,1:nCandidateModes);
expectedPotentialEnergySpectrum = targetModeSpectrum(nCandidateModes, jStar);

%% Define a partial-depth observing grid
[zObs,samplingDescription] = observationGridForScenario(samplingScenario, nObs);

imObs = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zObs,latitude=latitude,nEVP=nEVP,nModes=nCandidateModes,g=g);
imObs.normalization = Normalization.geostrophic;
imObs.upperBoundary = UpperBoundary.rigidLid;
[~,GObs] = imObs.modesAtFrequency(omega);
B = GObs(:,1:nCandidateModes);
% B is nObs-by-nCandidateModes. Row r is the modal prediction at observation
% r, and column j is the observational fingerprint of mode j.

%% Select modes that are identifiable by the observation operator
% The observation matrix B = H Phi maps modal coefficients to samples. The
% observation weights define the data-space metric, not the physical
% Parseval quadrature. Here equal weights represent equal observation quality.
weightsObs = ones(size(zObs));
% weightsObs is the diagonal of W_obs. Use inverse variances here when
% observations have different error levels.
rankTolerance = 5e-2;
maxConditionNumber = 10;
[retainedModes,selectionDiagnostics] = selectResolvableModes(B, weightsObs, rankTolerance, maxConditionNumber);
prefixModes = 1:length(retainedModes);

retainedDiagnostics = projectionDiagnostics(B, weightsObs, retainedModes);
prefixDiagnostics = projectionDiagnostics(B, weightsObs, prefixModes);
observedSingularValues = svd(sqrt(weightsObs).*B,'econ');
normalizedObservedSingularValues = observedSingularValues/observedSingularValues(1);

fprintf('Sampling scenario: %s (%s)\n', char(string(samplingScenario)), char(samplingDescription));
fprintf('Realization scenario: %s\n', char(string(realizationScenario)));
fprintf('Target ensemble spectrum H(j)=A/(j_*^2+j^2), j_*=%.3g, sum(H)=%.12f.\n', jStar, sum(expectedPotentialEnergySpectrum));
fprintf('Observation grid has %d samples over %.0f m to %.0f m.\n', length(zObs), min(zObs), max(zObs));
fprintf('Candidate hydrostatic G modes: %d.\n', nCandidateModes);
fprintf('Observed modal matrix B has size [%d observations x %d candidate modes].\n', size(B,1), size(B,2));
fprintf('  Row r: modal prediction at observation r; column j: fingerprint of mode j.\n');
fprintf('weightsObs defines diagonal W_obs in observation space; it is not physical quadrature.\n');
fprintf('The retained set I contains mode fingerprints that remain independent on this grid.\n');
fprintf('Rank-revealing retained modes: %s\n', mat2str(retainedModes));
fprintf('Contiguous prefix comparison:  %s\n', mat2str(prefixModes));
fprintf('QR relative pivot at final retained mode: %.3e\n', selectionDiagnostics.finalRelativePivot);
fprintf('Condition number after QR column normalization: %.3e\n', selectionDiagnostics.conditionNumber);
fprintf('Final projection is built from unscaled physical columns B(:,retainedModes).\n');

fprintf('\nProjection diagnostics\n');
printProjectionSummary("rank-revealing",retainedDiagnostics);
printProjectionSummary("same-count prefix",prefixDiagnostics);

%% Recover coefficients from synthetic observations
% The true coefficients are a Gaussian draw from an ensemble with
% E[g*eta_j^2]=H(j). Multiplying one realization by g gives its realized
% potential-energy spectrum.
trueCoefficients = gaussianSpectrumCoefficients(expectedPotentialEnergySpectrum, g, realizationScenario);

trueProfile = PhiPlot*trueCoefficients;
observations = B*trueCoefficients;

recoveredRetainedCoefficients = retainedDiagnostics.forwardProjection*observations;
recoveredCoefficients = zeros(nCandidateModes,1);
recoveredCoefficients(retainedModes) = recoveredRetainedCoefficients;
projectedObservations = retainedDiagnostics.observationProjectionMatrix*observations;
observationResidual = observations - projectedObservations;
observationVariance = weightedObservationVariance(observations, weightsObs);
observationResidualFraction = weightedObservationVariance(observationResidual, weightsObs)/observationVariance;

predictedRecoveredCoefficients = retainedDiagnostics.resolutionMatrix*trueCoefficients;
unresolvedBias = retainedDiagnostics.resolutionMatrix(:,retainedDiagnostics.rejectedModes)*trueCoefficients(retainedDiagnostics.rejectedModes);
selectedRelativeError = norm(recoveredRetainedCoefficients - trueCoefficients(retainedModes))/norm(trueCoefficients(retainedModes));

truePotentialEnergySpectrum = g*trueCoefficients.^2;
recoveredPotentialEnergySpectrum = g*recoveredCoefficients.^2;
expectedRecoveredSpectrum = zeros(nCandidateModes,1);
expectedRecoveredSpectrum(retainedModes) = retainedDiagnostics.spectralWindow*expectedPotentialEnergySpectrum;

fprintf('\nSynthetic recovery\n');
fprintf('Relative error in retained coefficients: %.3e\n', selectedRelativeError);
fprintf('Relative unresolved-mode bias: %.3e\n', norm(unresolvedBias)/norm(trueCoefficients(retainedModes)));
fprintf('Deterministic resolution check: %.3e\n', norm(predictedRecoveredCoefficients - retainedDiagnostics.forwardProjection*observations));
fprintf('Synthetic observation residual fraction: %.3e\n', observationResidualFraction);
fprintf('Realized true potential energy: %.3e\n', sum(truePotentialEnergySpectrum));
fprintf('Recovered potential energy: %.3e\n', sum(recoveredPotentialEnergySpectrum));
fprintf('Expected recovered ensemble potential energy: %.3e\n', sum(expectedRecoveredSpectrum));
fprintf('Expected spectrum assumes uncorrelated coefficients; deterministic recovery can include coherent aliasing.\n');

%% Plot true and recovered coefficients with resolution diagnostics
figure('Name',sprintf('Partial-depth observational mode projection: %s, %s', char(string(samplingScenario)), char(string(realizationScenario))))
tiledlayout(4,2,'TileSpacing','compact','Padding','compact')

nexttile
plot(trueProfile,zPlot,'k','LineWidth',1.5), hold on
plot(observations,zObs,'ko','MarkerFaceColor','w')
xlabel('profile value')
ylabel('z (m)')
title('Synthetic partial-depth observations')
legend('truth','samples','Location','southwest')
grid on

nexttile
stem(1:nCandidateModes,trueCoefficients,'k','DisplayName','true'), hold on
stem(retainedModes,recoveredCoefficients(retainedModes),'r','filled','DisplayName','recovered')
xlabel('mode number')
ylabel('coefficient')
title('Hydrostatic G-mode coefficients')
legend('Location','northeast')
grid on

nexttile
stem(1:nCandidateModes,expectedPotentialEnergySpectrum,'Color',[0.45 0.45 0.45],'DisplayName','target H'), hold on
stem(1:nCandidateModes,truePotentialEnergySpectrum,'k','DisplayName','realization')
stem(retainedModes,expectedRecoveredSpectrum(retainedModes),'b','filled','DisplayName','expected recovered')
stem(retainedModes,recoveredPotentialEnergySpectrum(retainedModes),'r','DisplayName','deterministic recovered')
xlabel('mode number')
ylabel('potential energy')
title('Resolvable potential-energy spectrum')
legend('Location','northeast')
grid on

nexttile
imagesc(1:nCandidateModes,retainedModes,abs(retainedDiagnostics.resolutionMatrix))
set(gca,'YDir','normal')
xlabel('true mode number')
ylabel('recovered mode number')
title('|A B| resolution/aliasing matrix')
colorbar

nexttile
stem(retainedDiagnostics.rejectedModes,retainedDiagnostics.aliasingColumnNorm,'filled')
xlabel('rejected mode number')
ylabel('aliasing norm')
title('Rejected-mode leakage into retained coefficients')
grid on

nexttile
imagesc(1:nCandidateModes,retainedModes,retainedDiagnostics.spectralWindow)
set(gca,'YDir','normal')
xlabel('true mode number')
ylabel('recovered mode number')
title('Potential-energy spectral window')
colorbar

nexttile
semilogy(1:length(normalizedObservedSingularValues),normalizedObservedSingularValues,'ko-'), hold on
xline(length(retainedModes),'r','DisplayName','retained count')
yline(rankTolerance,'b','DisplayName','rank tolerance')
xlabel('singular value index')
ylabel('normalized value')
title('Observed modal singular values')
legend('singular values','retained count','rank tolerance','Location','southwest')
grid on

nexttile
stem(1:nCandidateModes,retainedDiagnostics.observationResidualColumnFraction,'filled')
xlabel('candidate mode number')
ylabel('residual fraction')
title('Unexplained observed fingerprint variance')
grid on

%% Local helpers
function [zObs,description] = observationGridForScenario(samplingScenario, nObs)
switch string(samplingScenario)
    case "irregularPartialDepth"
        zObsUniform = linspace(-1600,-80,nObs).';
        zObs = zObsUniform + 35*sin((1:length(zObsUniform)).');
        zObs = sort(zObs);
        description = string(sprintf('irregular samples from %.0f m to %.0f m', min(zObs), max(zObs)));
    case "uniformUpper300m"
        zObs = linspace(-380,-80,nObs).';
        description = "uniform 300 m window from -380 m to -80 m";
    otherwise
        error("ObservationalProjectionPartialDepth:UnknownSamplingScenario", ...
            "Unknown samplingScenario ""%s"".", char(string(samplingScenario)));
end
end

function spectrum = targetModeSpectrum(nCandidateModes, jStar)
if jStar <= 0
    error("ObservationalProjectionPartialDepth:InvalidJStar", "jStar must be positive.");
end

modeNumbers = (1:nCandidateModes).';
shape = 1./(jStar*jStar + modeNumbers.^2);
spectrum = shape/sum(shape);
end

function coefficients = gaussianSpectrumCoefficients(expectedPotentialEnergySpectrum, g, realizationScenario)
switch string(realizationScenario)
    case "gaussianSeedA"
        rng(42,"twister");
    case "gaussianSeedB"
        rng(84,"twister");
    case "gaussianSeedC"
        rng(126,"twister");
    otherwise
        error("ObservationalProjectionPartialDepth:UnknownRealizationScenario", ...
            "Unknown realizationScenario ""%s"".", char(string(realizationScenario)));
end

coefficients = sqrt(expectedPotentialEnergySpectrum/g).*randn(size(expectedPotentialEnergySpectrum));
end

function [retainedModes,diagnostics] = selectResolvableModes(B, weights, rankTolerance, maxConditionNumber)
weightedB = sqrt(weights) .* B;
columnNorms = vecnorm(weightedB,2,1);
normalizedB = weightedB ./ columnNorms;
[~,R,pivotOrder] = qr(normalizedB,0);
relativePivots = abs(diag(R))/abs(R(1,1));
nRank = find(relativePivots >= rankTolerance,1,'last');
if isempty(nRank)
    nRank = 0;
end

selectedInPivotOrder = zeros(1,nRank);
nSelected = 0;
for iPivot = 1:nRank
    candidateModes = sort([selectedInPivotOrder(1:nSelected) pivotOrder(iPivot)]);
    candidateConditionNumber = cond(normalizedB(:,candidateModes).'*normalizedB(:,candidateModes));
    if candidateConditionNumber <= maxConditionNumber
        nSelected = nSelected + 1;
        selectedInPivotOrder(nSelected) = pivotOrder(iPivot);
    end
end

selectedInPivotOrder = selectedInPivotOrder(1:nSelected);
retainedModes = sort(selectedInPivotOrder);
diagnostics.relativePivots = relativePivots;
if isempty(retainedModes)
    diagnostics.finalRelativePivot = NaN;
    diagnostics.conditionNumber = NaN;
else
    diagnostics.finalRelativePivot = relativePivots(nSelected);
    diagnostics.conditionNumber = cond(normalizedB(:,retainedModes).'*normalizedB(:,retainedModes));
end
end

function diagnostics = projectionDiagnostics(B, weights, retainedModes)
nCandidateModes = size(B,2);
diagnostics.retainedModes = retainedModes;
diagnostics.rejectedModes = setdiff(1:nCandidateModes, retainedModes);
[diagnostics.forwardProjection,diagnostics.resolutionMatrix,diagnostics.gramMatrix] = weightedProjection(B, weights, retainedModes);
diagnostics.observationProjectionMatrix = B(:,retainedModes)*diagnostics.forwardProjection;
observationResidualMatrix = B - diagnostics.observationProjectionMatrix*B;
diagnostics.observationVariance = weightedObservationVariance(B, weights);
diagnostics.observationResidualVariance = weightedObservationVariance(observationResidualMatrix, weights);
diagnostics.observationExplainedFraction = 1 - diagnostics.observationResidualVariance/diagnostics.observationVariance;
diagnostics.observationResidualFraction = diagnostics.observationResidualVariance/diagnostics.observationVariance;
observationResidualNumerator = sum(weights .* observationResidualMatrix.^2,1);
observationResidualDenominator = sum(weights .* B.^2,1);
diagnostics.observationResidualColumnFraction = observationResidualNumerator ./ observationResidualDenominator;
diagnostics.retainedError = norm(diagnostics.resolutionMatrix(:,retainedModes) - eye(length(retainedModes)),'fro')/sqrt(length(retainedModes));
diagnostics.aliasingMatrix = diagnostics.resolutionMatrix(:,diagnostics.rejectedModes);
diagnostics.aliasingColumnNorm = vecnorm(diagnostics.aliasingMatrix,2,1);
diagnostics.maxAliasing = max([0 diagnostics.aliasingColumnNorm]);
diagnostics.conditionNumber = cond(diagnostics.gramMatrix);
diagnostics.spectralWindow = diagnostics.resolutionMatrix.^2;
rejectedWindow = diagnostics.spectralWindow(:,diagnostics.rejectedModes);
diagnostics.windowLeakage = sum(rejectedWindow(:))/size(diagnostics.spectralWindow,1);
end

function [A,resolutionMatrix,gramMatrix] = weightedProjection(B, weights, retainedModes)
BS = B(:,retainedModes);
weightedBS = weights .* BS;
gramMatrix = BS.' * weightedBS;
A = gramMatrix \ weightedBS.';
resolutionMatrix = A*B;
end

function variance = weightedObservationVariance(values, weights)
variance = sum(weights .* sum(values.^2,2));
end

function printProjectionSummary(label,diagnostics)
fprintf('%-18s retained=%2d cond=%9.3e retained err=%9.3e max alias=%9.3e window leak=%9.3e\n', ...
    label,length(diagnostics.retainedModes),diagnostics.conditionNumber, ...
    diagnostics.retainedError,diagnostics.maxAliasing,diagnostics.windowLeakage);
fprintf('  observed candidate variance explained=%9.3e residual=%9.3e\n', ...
    diagnostics.observationExplainedFraction,diagnostics.observationResidualFraction);
end
