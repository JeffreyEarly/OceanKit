%% Compare hydrostatic partial-window formulations
% This example mirrors the main formulations in Section 2.4 of
% "vertical modes and how to use them" and adds one diagnostic variant:
%
% 1. F window/HKE formulation using int u^2 dz.
% 2. Standard G window/PE formulation using int N2 eta^2 dz.
% 3. Signal-dependent signed G interior formulation with the Section 2.4
%    boundary term evaluated from the synthetic displacement eta.
% 4. Signal-dependent norm-like G interior variant using absolute eta at the
%    window boundaries.
%
% The main diagnostic is Parseval consistency: the spatial or quadratic-form
% energy should match the sum of the plotted window spectrum.

%% User settings
if ~exist('stratificationCase','var')
    stratificationCase = "exponential";
end
if ~exist('windowBounds','var')
    if string(stratificationCase) == "constant"
        windowBounds = [-650 0];
    else
        windowBounds = [-1500 -500];
    end
end
if ~exist('nBaroclinicModes','var')
    nBaroclinicModes = 100;
end
if ~exist('nShow','var')
    nShow = 6;
end
if ~exist('nPositiveSignedModes','var')
    nPositiveSignedModes = 5;
end
if ~exist('fModePlotLimits','var')
    fModePlotLimits = [];
end
if ~exist('windowModeNormalization','var')
    if exist('fWindowEigenmodes','var')
        switch string(fWindowEigenmodes)
            case "normalized"
                windowModeNormalization = "unit";
            case "raw"
                windowModeNormalization = "equivalentDepth";
            otherwise
                error("HydrostaticWindowModeComparison:UnknownFWindowEigenmodes", ...
                    "fWindowEigenmodes must be ""raw"" or ""normalized"".");
        end
    else
        windowModeNormalization = "equivalentDepth";
    end
end
if ~exist('nSpectrumPlot','var')
    nSpectrumPlot = 40;
end
if ~exist('nReconstructionModes','var')
    nReconstructionModes = 20;
end
if ~exist('nPlot','var')
    nPlot = 1024;
end
if ~exist('nQuad','var')
    nQuad = 4097;
end
if ~exist('jStar','var')
    jStar = 3;
end
if ~exist('etaPeakMeters','var')
    etaPeakMeters = 300;
end
if ~exist('testVelocityEquivalentDepth','var')
    testVelocityEquivalentDepth = 0.80;
end
if ~exist('etaSignConvention','var')
    etaSignConvention = "positiveLowerWindow";
end
if ~exist('realizationSeed','var')
    realizationSeed = 41;
end

latitude = 31;
g = 9.81;

%% Build high-quality hydrostatic modes on plotting and window grids
modeData = buildHydrostaticModeData(stratificationCase,windowBounds,nBaroclinicModes,nPlot,nQuad,latitude,g);
nF = nBaroclinicModes + 1;
nG = nBaroclinicModes;
nShowF = min(nShow,nF);
nShowG = min(nShow,nG);
nSpectrumPlotF = min(nSpectrumPlot,nF);
nSpectrumPlotG = min(nSpectrumPlot,nG);
nReconstructionF = min(nReconstructionModes,nF);
nReconstructionG = min(nReconstructionModes,nG);

gammaF0 = [modeData.depth; modeData.h(:)];
invSqrtGammaF0 = 1./sqrt(gammaF0);
switch string(windowModeNormalization)
    case "unit"
        fScale = invSqrtGammaF0;
        gScale = ones(nG,1);
        fFullDepthMetric = ones(nF,1);
        gFullDepthMetric = ones(nG,1);
        normalizationLabel = "unit";
    case "equivalentDepth"
        fScale = ones(nF,1);
        gScale = sqrt(modeData.h);
        fFullDepthMetric = gammaF0;
        gFullDepthMetric = modeData.h;
        normalizationLabel = "equivalent-depth";
    otherwise
        error("HydrostaticWindowModeComparison:UnknownWindowModeNormalization", ...
            "windowModeNormalization must be ""unit"" or ""equivalentDepth"".");
end

weightsFWindow = trapzWeights(modeData.zWindow);
weightsGWindow = weightsFWindow .* modeData.N2Window/g;
weightsFPlot = trapzWeights(modeData.zPlot);
weightsPEPlot = weightsFPlot .* modeData.N2Plot;

PhiFPlot = modeData.PhiFPlot;
PhiFWindow = modeData.PhiFWindow;
PhiGPlot = modeData.PhiGPlot;
PhiGWindow = modeData.PhiGWindow;
PhiFAnalysisPlot = PhiFPlot.*fScale.';
PhiFAnalysisWindow = PhiFWindow.*fScale.';
PhiGAnalysisPlot = PhiGPlot.*gScale.';
PhiGAnalysisWindow = PhiGWindow.*gScale.';

GammaFWin = PhiFAnalysisWindow.'*(weightsFWindow.*PhiFAnalysisWindow);
GammaFWin = (GammaFWin + GammaFWin.')/2;

GammaGWin = PhiGAnalysisWindow.'*(weightsGWindow.*PhiGAnalysisWindow);
GammaGWin = (GammaGWin + GammaGWin.')/2;

%% Synthetic fields
rng(realizationSeed,"twister");
etaShape = targetModeSpectrum((1:nBaroclinicModes).',jStar);
velocityCoefficientScale = sqrt(g/testVelocityEquivalentDepth);

etaHat = sqrt(etaShape/g).*randn(nG,1);
etaPlotUnscaled = PhiGPlot*etaHat;
etaScale = etaPeakMeters/max(abs(etaPlotUnscaled));
etaHat = etaScale*etaHat;
etaPlot = PhiGPlot*etaHat;
etaWindow = PhiGWindow*etaHat;
[etaHat,etaPlot,etaWindow] = applyEtaSignConvention(etaHat,etaPlot,etaWindow,etaSignConvention);
uHat = [0; velocityCoefficientScale*etaHat];
uPlot = PhiFPlot*uHat;
uWindow = PhiFWindow*uHat;
uAnalysisHat = uHat./fScale;
etaAnalysisHat = etaHat./gScale;

%% Signal-dependent G interior forms
GAtLower = PhiGAnalysisWindow(1,:).';
GAtUpper = PhiGAnalysisWindow(end,:).';
% The signed bracket convention is [f]_{z_a}^{z_b} = f(z_b) - f(z_a).
boundaryTermSigned = (modeData.N2Window(end)*etaWindow(end)*(GAtUpper*GAtUpper.') - modeData.N2Window(1)*etaWindow(1)*(GAtLower*GAtLower.'));
% boundaryTermSigned = -modeData.N2Window(end)*abs(etaWindow(end))*(GAtUpper*GAtUpper.') - modeData.N2Window(1)*abs(etaWindow(1))*(GAtLower*GAtLower.');
boundaryTermAbsEta = modeData.N2Window(end)*abs(etaWindow(end))*(GAtUpper*GAtUpper.') + modeData.N2Window(1)*abs(etaWindow(1))*(GAtLower*GAtLower.');
GammaGIntSigned = GammaGWin + boundaryTermSigned/g;
GammaGIntAbsEta = GammaGWin + boundaryTermAbsEta/g;
GammaGIntSigned = (GammaGIntSigned + GammaGIntSigned.')/2;
GammaGIntAbsEta = (GammaGIntAbsEta + GammaGIntAbsEta.')/2;

%% Diagonalize the quadratic forms
[QF,lambdaF] = sortedEigenvectors(GammaFWin,"descend");
[QG,lambdaG] = sortedEigenvectors(GammaGWin,"descend");
[QGIntSigned,lambdaGIntSigned] = sortedEigenvectors(GammaGIntSigned,"descend");
[QGIntAbsEta,lambdaGIntAbsEta] = sortedEigenvectors(GammaGIntAbsEta,"descend");
QFfromG = blkdiag(1,QG);
QFfromGIntSigned = blkdiag(1,QGIntSigned);
GammaFfromGWin = QFfromG.'*GammaFWin*QFfromG;
GammaFfromGWin = (GammaFfromGWin + GammaFfromGWin.')/2;
GammaFfromGIntSignedWin = QFfromGIntSigned.'*GammaFWin*QFfromGIntSigned;
GammaFfromGIntSignedWin = (GammaFfromGIntSignedWin + GammaFfromGIntSignedWin.')/2;
offDiagonalFfromG = offDiagonalNorm(GammaFfromGWin)/max(norm(GammaFfromGWin,'fro'),eps);
offDiagonalFfromGIntSigned = offDiagonalNorm(GammaFfromGIntSignedWin)/max(norm(GammaFfromGIntSignedWin,'fro'),eps);

fVisibility = modeVisibility(QF,lambdaF,fFullDepthMetric);
gVisibility = modeVisibility(QG,lambdaG,gFullDepthMetric);
gIntSignedVisibility = modeVisibility(QGIntSigned,lambdaGIntSigned,gFullDepthMetric);
gIntAbsEtaVisibility = modeVisibility(QGIntAbsEta,lambdaGIntAbsEta,gFullDepthMetric);
fLabel = "F " + normalizationLabel + " window HKE";
fFromGLabel = "F from standard G window modes";
fFromGIntSignedLabel = "F from signed G window modes";
gLabel = "G " + normalizationLabel + " window PE";

plotIndexF = (1:min(nShowF,length(lambdaF))).';
plotIndexFfromG = (1:min(nShowF,size(QFfromG,2))).';
plotIndexG = (1:min(nShowG,length(lambdaG))).';
plotIndexGIntSigned = signedModePlotIndices(lambdaGIntSigned,nPositiveSignedModes);
plotIndexFfromGIntSigned = [1; plotIndexGIntSigned + 1];
plotIndexGIntAbsEta = (1:min(nShowG,length(lambdaGIntAbsEta))).';
plotLabelGIntSigned = signedModeLegendLabels(lambdaGIntSigned,plotIndexGIntSigned);
plotLabelFfromGIntSigned = ["F0"; plotLabelGIntSigned];

PsiFPlot = PhiFAnalysisPlot*QF(:,plotIndexF);
PsiFfromGPlot = PhiFAnalysisPlot*QFfromG(:,plotIndexFfromG);
PsiFfromGIntSignedPlot = PhiFAnalysisPlot*QFfromGIntSigned(:,plotIndexFfromGIntSigned);
PsiGPlot = PhiGAnalysisPlot*QG(:,plotIndexG);
PsiGIntSignedPlot = PhiGAnalysisPlot*QGIntSigned(:,plotIndexGIntSigned);
PsiGIntAbsEtaPlot = PhiGAnalysisPlot*QGIntAbsEta(:,plotIndexGIntAbsEta);
if isempty(fModePlotLimits)
    fModePlotLimits = symmetricProfileLimits([PsiFPlot(:,2:end), ...
        PsiFfromGPlot(:,2:end),PsiFfromGIntSignedPlot(:,2:end)],1.15);
end

cF = QF.'*uAnalysisHat;
cG = QG.'*etaAnalysisHat;
cGIntSigned = QGIntSigned.'*etaAnalysisHat;
cGIntAbsEta = QGIntAbsEta.'*etaAnalysisHat;

spectrumF = lambdaF.*cF.^2;
spectrumG = g*lambdaG.*cG.^2;
spectrumGIntSigned = g*lambdaGIntSigned.*cGIntSigned.^2;
spectrumGIntAbsEta = g*lambdaGIntAbsEta.*cGIntAbsEta.^2;

reconstructionIndexF = dominantSpectrumIndices(spectrumF,nReconstructionF);
reconstructionIndexFfromG = (1:nReconstructionF).';
reconstructionIndexFfromGIntSigned = (1:nReconstructionF).';
reconstructionIndexG = dominantSpectrumIndices(spectrumG,nReconstructionG);
reconstructionIndexGIntSigned = dominantSignedContributionIndices(spectrumGIntSigned,nReconstructionG);
reconstructionIndexGIntAbsEta = dominantSignedContributionIndices(spectrumGIntAbsEta,nReconstructionG);

uReconstructionHat = fScale.*(QF(:,reconstructionIndexF)*cF(reconstructionIndexF));
uFromGReconstructionHat = fWindowGalerkinReconstruction(QFfromG,reconstructionIndexFfromG, ...
    fScale,uWindow,weightsFWindow,modeData);
uFromGIntSignedReconstructionHat = fWindowGalerkinReconstruction(QFfromGIntSigned, ...
    reconstructionIndexFfromGIntSigned,fScale,uWindow,weightsFWindow,modeData);
etaReconstructionHat = gScale.*(QG(:,reconstructionIndexG)*cG(reconstructionIndexG));
etaInteriorSignedReconstructionHat = gScale.*(QGIntSigned(:,reconstructionIndexGIntSigned)*cGIntSigned(reconstructionIndexGIntSigned));
etaInteriorAbsEtaReconstructionHat = gScale.*(QGIntAbsEta(:,reconstructionIndexGIntAbsEta)*cGIntAbsEta(reconstructionIndexGIntAbsEta));

uReconstructionPlot = PhiFPlot*uReconstructionHat;
uFromGReconstructionPlot = PhiFPlot*uFromGReconstructionHat;
uFromGIntSignedReconstructionPlot = PhiFPlot*uFromGIntSignedReconstructionHat;
etaReconstructionPlot = PhiGPlot*etaReconstructionHat;
etaInteriorSignedReconstructionPlot = PhiGPlot*etaInteriorSignedReconstructionHat;
etaInteriorAbsEtaReconstructionPlot = PhiGPlot*etaInteriorAbsEtaReconstructionHat;

uResidualPlot = uPlot - uReconstructionPlot;
uFromGResidualPlot = uPlot - uFromGReconstructionPlot;
uFromGIntSignedResidualPlot = uPlot - uFromGIntSignedReconstructionPlot;
etaResidualPlot = etaPlot - etaReconstructionPlot;
etaInteriorSignedResidualPlot = etaPlot - etaInteriorSignedReconstructionPlot;
etaInteriorAbsEtaResidualPlot = etaPlot - etaInteriorAbsEtaReconstructionPlot;
uResidualWindow = uWindow - PhiFWindow*uReconstructionHat;
uFromGResidualWindow = uWindow - PhiFWindow*uFromGReconstructionHat;
uFromGIntSignedResidualWindow = uWindow - PhiFWindow*uFromGIntSignedReconstructionHat;
etaResidualWindow = etaWindow - PhiGWindow*etaReconstructionHat;
etaInteriorSignedResidualWindow = etaWindow - PhiGWindow*etaInteriorSignedReconstructionHat;
etaInteriorAbsEtaResidualWindow = etaWindow - PhiGWindow*etaInteriorAbsEtaReconstructionHat;
uResidualHat = uHat - uReconstructionHat;
uFromGResidualHat = uHat - uFromGReconstructionHat;
uFromGIntSignedResidualHat = uHat - uFromGIntSignedReconstructionHat;
etaResidualHat = etaHat - etaReconstructionHat;
etaInteriorSignedResidualHat = etaHat - etaInteriorSignedReconstructionHat;
etaInteriorAbsEtaResidualHat = etaHat - etaInteriorAbsEtaReconstructionHat;

duDz = fCoefficientDerivative(uHat,modeData,g);
duResidualDz = fCoefficientDerivative(uResidualHat,modeData,g);
duFromGResidualDz = fCoefficientDerivative(uFromGResidualHat,modeData,g);
duFromGIntSignedResidualDz = fCoefficientDerivative(uFromGIntSignedResidualHat,modeData,g);
duDzWindow = fCoefficientDerivativeOnGrid(uHat,modeData.PhiGWindow,modeData.N2Window,g);
duResidualDzWindow = fCoefficientDerivativeOnGrid(uResidualHat,modeData.PhiGWindow,modeData.N2Window,g);
duFromGResidualDzWindow = fCoefficientDerivativeOnGrid(uFromGResidualHat, ...
    modeData.PhiGWindow,modeData.N2Window,g);
duFromGIntSignedResidualDzWindow = fCoefficientDerivativeOnGrid(uFromGIntSignedResidualHat, ...
    modeData.PhiGWindow,modeData.N2Window,g);
detaDz = gCoefficientDerivative(etaHat,modeData);
detaResidualDz = gCoefficientDerivative(etaResidualHat,modeData);
detaDzWindow = gCoefficientDerivativeOnGrid(etaHat,modeData.PhiFWindow,modeData.h);
detaResidualDzWindow = gCoefficientDerivativeOnGrid(etaResidualHat,modeData.PhiFWindow,modeData.h);
detaInteriorSignedResidualDzWindow = gCoefficientDerivativeOnGrid(etaInteriorSignedResidualHat,modeData.PhiFWindow,modeData.h);
detaInteriorAbsEtaResidualDzWindow = gCoefficientDerivativeOnGrid(etaInteriorAbsEtaResidualHat,modeData.PhiFWindow,modeData.h);
detaInteriorSignedResidualDz = gCoefficientDerivative(etaInteriorSignedResidualHat,modeData);
detaInteriorAbsEtaResidualDz = gCoefficientDerivative(etaInteriorAbsEtaResidualHat,modeData);

spatialEnergyF = weightedEnergy(uWindow,weightsFWindow);
spectralEnergyF = sum(spectrumF);
quadraticEnergyF = uAnalysisHat.'*GammaFWin*uAnalysisHat;

spatialEnergyG = weightedEnergy(etaWindow,modeData.N2Window.*weightsFWindow);
spectralEnergyG = sum(spectrumG);
quadraticEnergyG = g*(etaAnalysisHat.'*GammaGWin*etaAnalysisHat);

quadraticEnergyGIntSigned = g*(etaAnalysisHat.'*GammaGIntSigned*etaAnalysisHat);
spectralEnergyGIntSigned = sum(spectrumGIntSigned);
quadraticEnergyGIntAbsEta = g*(etaAnalysisHat.'*GammaGIntAbsEta*etaAnalysisHat);
spectralEnergyGIntAbsEta = sum(spectrumGIntAbsEta);

fullEnergyF = weightedEnergy(uPlot,weightsFPlot);
residualEnergyF = weightedEnergy(uResidualPlot,weightsFPlot);
residualEnergyFfromG = weightedEnergy(uFromGResidualPlot,weightsFPlot);
residualEnergyFfromGIntSigned = weightedEnergy(uFromGIntSignedResidualPlot,weightsFPlot);
windowResidualEnergyF = weightedEnergy(uResidualWindow,weightsFWindow);
windowResidualEnergyFfromG = weightedEnergy(uFromGResidualWindow,weightsFWindow);
windowResidualEnergyFfromGIntSigned = weightedEnergy(uFromGIntSignedResidualWindow,weightsFWindow);
fullEnergyG = weightedEnergy(etaPlot,weightsPEPlot);
residualEnergyG = weightedEnergy(etaResidualPlot,weightsPEPlot);
residualEnergyGIntSigned = weightedEnergy(etaInteriorSignedResidualPlot,weightsPEPlot);
residualEnergyGIntAbsEta = weightedEnergy(etaInteriorAbsEtaResidualPlot,weightsPEPlot);
windowResidualEnergyG = weightedEnergy(etaResidualWindow,modeData.N2Window.*weightsFWindow);
windowResidualEnergyGIntSigned = weightedEnergy(etaInteriorSignedResidualWindow,modeData.N2Window.*weightsFWindow);
windowResidualEnergyGIntAbsEta = weightedEnergy(etaInteriorAbsEtaResidualWindow,modeData.N2Window.*weightsFWindow);
windowCapturedF = capturedEnergy(windowResidualEnergyF,spatialEnergyF);
windowCapturedFfromG = capturedEnergy(windowResidualEnergyFfromG,spatialEnergyF);
windowCapturedFfromGIntSigned = capturedEnergy(windowResidualEnergyFfromGIntSigned,spatialEnergyF);
windowCapturedG = capturedEnergy(windowResidualEnergyG,spatialEnergyG);
windowCapturedGIntSigned = capturedEnergy(windowResidualEnergyGIntSigned,spatialEnergyG);
windowCapturedGIntAbsEta = capturedEnergy(windowResidualEnergyGIntAbsEta,spatialEnergyG);
fullCapturedF = capturedEnergy(residualEnergyF,fullEnergyF);
fullCapturedFfromG = capturedEnergy(residualEnergyFfromG,fullEnergyF);
fullCapturedFfromGIntSigned = capturedEnergy(residualEnergyFfromGIntSigned,fullEnergyF);
fullCapturedG = capturedEnergy(residualEnergyG,fullEnergyG);
fullCapturedGIntSigned = capturedEnergy(residualEnergyGIntSigned,fullEnergyG);
fullCapturedGIntAbsEta = capturedEnergy(residualEnergyGIntAbsEta,fullEnergyG);
fullDerivativeEnergyF = inverseN2DerivativeEnergy(duDz,weightsFPlot,modeData.N2Plot,g);
residualDerivativeEnergyF = inverseN2DerivativeEnergy(duResidualDz,weightsFPlot,modeData.N2Plot,g);
residualDerivativeEnergyFfromG = inverseN2DerivativeEnergy(duFromGResidualDz,weightsFPlot,modeData.N2Plot,g);
residualDerivativeEnergyFfromGIntSigned = inverseN2DerivativeEnergy(duFromGIntSignedResidualDz, ...
    weightsFPlot,modeData.N2Plot,g);
windowDerivativeEnergyF = inverseN2DerivativeEnergy(duDzWindow,weightsFWindow,modeData.N2Window,g);
windowResidualDerivativeEnergyF = inverseN2DerivativeEnergy(duResidualDzWindow,weightsFWindow,modeData.N2Window,g);
windowResidualDerivativeEnergyFfromG = inverseN2DerivativeEnergy(duFromGResidualDzWindow, ...
    weightsFWindow,modeData.N2Window,g);
windowResidualDerivativeEnergyFfromGIntSigned = inverseN2DerivativeEnergy(duFromGIntSignedResidualDzWindow, ...
    weightsFWindow,modeData.N2Window,g);
fullDerivativeEnergyG = weightedEnergy(detaDz,weightsFPlot);
residualDerivativeEnergyG = weightedEnergy(detaResidualDz,weightsFPlot);
windowDerivativeEnergyG = weightedEnergy(detaDzWindow,weightsFWindow);
residualDerivativeEnergyGIntSigned = weightedEnergy(detaInteriorSignedResidualDz,weightsFPlot);
residualDerivativeEnergyGIntAbsEta = weightedEnergy(detaInteriorAbsEtaResidualDz,weightsFPlot);
windowResidualDerivativeEnergyG = weightedEnergy(detaResidualDzWindow,weightsFWindow);
windowResidualDerivativeEnergyGIntSigned = weightedEnergy(detaInteriorSignedResidualDzWindow,weightsFWindow);
windowResidualDerivativeEnergyGIntAbsEta = weightedEnergy(detaInteriorAbsEtaResidualDzWindow,weightsFWindow);
convergenceF = fResidualConvergence(QF,cF,spectrumF,fScale,uHat,uPlot,weightsFPlot, ...
    uWindow,weightsFWindow,modeData,g,nF,fullEnergyF,spatialEnergyF, ...
    fullDerivativeEnergyF,windowDerivativeEnergyF);
convergenceFfromG = fConstructionOrderConvergence(QFfromG,fScale,uHat, ...
    uPlot,weightsFPlot,uWindow,weightsFWindow,modeData,g,nF, ...
    fullEnergyF,spatialEnergyF,fullDerivativeEnergyF,windowDerivativeEnergyF);
convergenceFfromGIntSigned = fConstructionOrderConvergence(QFfromGIntSigned, ...
    fScale,uHat,uPlot,weightsFPlot,uWindow,weightsFWindow,modeData,g,nF, ...
    fullEnergyF,spatialEnergyF,fullDerivativeEnergyF,windowDerivativeEnergyF);
convergenceG = gResidualConvergence(QG,cG,spectrumG,gScale,etaHat,etaPlot,weightsPEPlot, ...
    etaWindow,modeData.N2Window.*weightsFWindow,weightsFPlot,weightsFWindow, ...
    modeData,nG,fullEnergyG,spatialEnergyG,fullDerivativeEnergyG,windowDerivativeEnergyG);
convergenceGIntSigned = gResidualConvergence(QGIntSigned,cGIntSigned,spectrumGIntSigned, ...
    gScale,etaHat,etaPlot,weightsPEPlot,etaWindow,modeData.N2Window.*weightsFWindow, ...
    weightsFPlot,weightsFWindow,modeData,nG,fullEnergyG,spatialEnergyG, ...
    fullDerivativeEnergyG,windowDerivativeEnergyG);
convergenceGIntAbsEta = gResidualConvergence(QGIntAbsEta,cGIntAbsEta,spectrumGIntAbsEta, ...
    gScale,etaHat,etaPlot,weightsPEPlot,etaWindow,modeData.N2Window.*weightsFWindow, ...
    weightsFPlot,weightsFWindow,modeData,nG,fullEnergyG,spatialEnergyG, ...
    fullDerivativeEnergyG,windowDerivativeEnergyG);

%% Print diagnostics
fprintf('Hydrostatic window formulation comparison\n');
fprintf('  stratification: %s\n', string(stratificationCase));
fprintf('  depth: %.0f m, window: [%.0f %.0f] m\n', modeData.depth, windowBounds(1), windowBounds(2));
fprintf('  F modes: %d = barotropic + %d baroclinic; G modes: %d baroclinic\n', nF, nBaroclinicModes, nG);
fprintf('  window mode normalization: %s\n', normalizationLabel);
fprintf('  eta sign convention: %s\n', string(etaSignConvention));
fprintf('  target spectrum A/(j_*^2+j^2), j_* = %.3g\n', jStar);
fprintf('  velocity coefficients: u_j = sqrt(g/%.3g) eta_j for j >= 1, u_0 = 0\n', testVelocityEquivalentDepth);
fprintf('  eta range: [%.3f %.3f] m, max(abs(eta)) = %.3f m\n', min(etaPlot), max(etaPlot), max(abs(etaPlot)));
fprintf('  eta at window boundaries: eta(z_a)=%.3f m, eta(z_b)=%.3f m\n', etaWindow(1), etaWindow(end));

fprintf('\nEnergy checks\n');
printEnergyCheck(fLabel, spatialEnergyF, quadraticEnergyF, spectralEnergyF);
printEnergyCheck(gLabel, spatialEnergyG, quadraticEnergyG, spectralEnergyG);
printEnergyCheck("G interior signed form", quadraticEnergyGIntSigned, quadraticEnergyGIntSigned, spectralEnergyGIntSigned);
printEnergyCheck("G interior abs-eta form", quadraticEnergyGIntAbsEta, quadraticEnergyGIntAbsEta, spectralEnergyGIntAbsEta);

fprintf('\nEigenvalue diagnostics\n');
fprintf('  F window eigenvalue range:    [%.4e %.4e]\n', min(lambdaF), max(lambdaF));
fprintf('  F visibility fraction range:  [%.4e %.4e]\n', min(fVisibility), max(fVisibility));
fprintf('  F-from-G window Gram offdiag: %.4e of Frobenius norm\n', offDiagonalFfromG);
fprintf('  F-from-signed-G Gram offdiag: %.4e of Frobenius norm\n', offDiagonalFfromGIntSigned);
fprintf('  G window lambda range:        [%.4e %.4e]\n', min(lambdaG), max(lambdaG));
fprintf('  G visibility fraction range:  [%.4e %.4e]\n', min(gVisibility), max(gVisibility));
printInteriorEigenvalueDiagnostics("G interior signed",lambdaGIntSigned);
printInteriorEigenvalueDiagnostics("G interior abs-eta",lambdaGIntAbsEta);

fprintf('\nSigned-form mode selection diagnostics\n');
printSignedSelectionDiagnostics("G interior signed",lambdaGIntSigned,spectrumGIntSigned,reconstructionIndexGIntSigned);
printSignedSelectionDiagnostics("G interior abs-eta",lambdaGIntAbsEta, ...
    spectrumGIntAbsEta,reconstructionIndexGIntAbsEta);

fprintf('\nFull-depth partial reconstructions from selected window modes\n');
printReconstructionCheck(fLabel, reconstructionIndexF, windowCapturedF, fullCapturedF);
printReconstructionCheck(fFromGLabel, reconstructionIndexFfromG, windowCapturedFfromG, fullCapturedFfromG);
printReconstructionCheck(fFromGIntSignedLabel, reconstructionIndexFfromGIntSigned, ...
    windowCapturedFfromGIntSigned, fullCapturedFfromGIntSigned);
printReconstructionCheck(gLabel, reconstructionIndexG, windowCapturedG, fullCapturedG);
printReconstructionCheck("G interior signed", reconstructionIndexGIntSigned, windowCapturedGIntSigned, fullCapturedGIntSigned);
printReconstructionCheck("G interior abs-eta", reconstructionIndexGIntAbsEta, windowCapturedGIntAbsEta, fullCapturedGIntAbsEta);

fprintf('\nFull-depth residual energies using straight definitions\n');
printResidualEnergy("F residual", residualEnergyF, fullEnergyF, "int r_u^2 dz");
printResidualEnergy("F-from-G residual", residualEnergyFfromG, fullEnergyF, "int r_u^2 dz");
printResidualEnergy("F-from-signed-G residual", residualEnergyFfromGIntSigned, fullEnergyF, "int r_u^2 dz");
printResidualEnergy("G residual", residualEnergyG, fullEnergyG, "int N2 r_eta^2 dz");
printResidualEnergy("G interior signed residual", residualEnergyGIntSigned, fullEnergyG, "int N2 r_eta^2 dz");
printResidualEnergy("G interior abs-eta residual", residualEnergyGIntAbsEta, fullEnergyG, "int N2 r_eta^2 dz");

fprintf('\nFull-depth derivative residual energies computed spectrally\n');
printResidualEnergy("F derivative residual", residualDerivativeEnergyF, fullDerivativeEnergyF, "g int r_u_z^2/N2 dz");
printResidualEnergy("F-from-G derivative residual", residualDerivativeEnergyFfromG, fullDerivativeEnergyF, "g int r_u_z^2/N2 dz");
printResidualEnergy("F-from-signed-G derivative", residualDerivativeEnergyFfromGIntSigned, ...
    fullDerivativeEnergyF, "g int r_u_z^2/N2 dz");
printResidualEnergy("G derivative residual", residualDerivativeEnergyG, fullDerivativeEnergyG, "int r_eta_z^2 dz");
printResidualEnergy("G signed derivative residual", residualDerivativeEnergyGIntSigned, fullDerivativeEnergyG, "int r_eta_z^2 dz");
printResidualEnergy("G abs-eta derivative residual", residualDerivativeEnergyGIntAbsEta, fullDerivativeEnergyG, "int r_eta_z^2 dz");

fprintf('\nFigure 1 displays ordinary panels in descending partial-window eigenvalue order.\n');
fprintf('  Signed-G panels show the first %d positive signed directions plus all negative signed directions.\n', ...
    nPositiveSignedModes);
fprintf('  Visibility fractions use the selected full-depth normalization as the denominator.\n');
printPlotOrder(fLabel, plotIndexF, fullDepthMetricScore(QF,fFullDepthMetric));
printPlotOrder(fFromGLabel, plotIndexFfromG, fullDepthMetricScore(QFfromG,fFullDepthMetric));
printPlotOrder(fFromGIntSignedLabel, plotIndexFfromGIntSigned, ...
    fullDepthMetricScore(QFfromGIntSigned,fFullDepthMetric));
printPlotOrder(gLabel, plotIndexG, fullDepthMetricScore(QG,gFullDepthMetric));
printPlotOrder("G interior signed", plotIndexGIntSigned, fullDepthMetricScore(QGIntSigned,gFullDepthMetric));
printPlotOrder("G interior abs-eta", plotIndexGIntAbsEta, fullDepthMetricScore(QGIntAbsEta,gFullDepthMetric));

%% Figure 1: leading modes
figure('Name',sprintf('Hydrostatic window modes: %s', string(stratificationCase)))
tiledlayout(2,3,TileSpacing="compact",Padding="compact")

nexttile
plotProfilesWithWindow(PsiFPlot,modeData.zPlot,windowBounds, ...
    "F window modes, eigenvalue order",compose('%d: %.2g',plotIndexF(:),lambdaF(plotIndexF)),fModePlotLimits)

nexttile
plotProfilesWithWindow(PsiFfromGPlot,modeData.zPlot,windowBounds, ...
    "F from standard G window modes",compose('%d',plotIndexFfromG(:)),fModePlotLimits)

nexttile
plotProfilesWithWindow(PsiFfromGIntSignedPlot,modeData.zPlot,windowBounds, ...
    "F from signed G window modes",plotLabelFfromGIntSigned,fModePlotLimits)

nexttile
plotProfilesWithWindow(PsiGPlot,modeData.zPlot,windowBounds, ...
    "G window modes, eigenvalue order",compose('%d: %.2g',plotIndexG(:),lambdaG(plotIndexG)))

nexttile
plotProfilesWithWindow(PsiGIntSignedPlot,modeData.zPlot,windowBounds, ...
    "G signed interior modes",plotLabelGIntSigned)

nexttile
plotProfilesWithWindow(PsiGIntAbsEtaPlot,modeData.zPlot,windowBounds, ...
    "G abs-eta interior modes, eigenvalue order",compose('%d: %.2g',plotIndexGIntAbsEta(:),lambdaGIntAbsEta(plotIndexGIntAbsEta)))

%% Figure 2: window spectra
figure('Name',sprintf('Hydrostatic window spectra: %s', string(stratificationCase)))
tiledlayout(2,2,TileSpacing="compact",Padding="compact")

nexttile
plotWindowSpectrum(spectrumF,nSpectrumPlotF,fLabel,spatialEnergyF,spectralEnergyF,reconstructionIndexF)

nexttile
plotWindowSpectrum(spectrumG,nSpectrumPlotG,gLabel,spatialEnergyG,spectralEnergyG,reconstructionIndexG)

nexttile
plotWindowSpectrum(spectrumGIntSigned,nSpectrumPlotG,"G signed interior form", ...
    quadraticEnergyGIntSigned,spectralEnergyGIntSigned,reconstructionIndexGIntSigned)

nexttile
plotWindowSpectrum(spectrumGIntAbsEta,nSpectrumPlotG,"G abs-eta interior form", ...
    quadraticEnergyGIntAbsEta,spectralEnergyGIntAbsEta,reconstructionIndexGIntAbsEta)

%% Figure 3: test profiles
figure('Name',sprintf('Hydrostatic window test profiles: %s', string(stratificationCase)))
tiledlayout(1,2,TileSpacing="compact",Padding="compact")

nexttile
plotProfileWithWindow(uPlot,modeData.zPlot,windowBounds,"test velocity u(z)")

nexttile
plotProfileWithWindow(etaPlot,modeData.zPlot,windowBounds,"test displacement \eta(z)")

%% Figure 4: full-depth partial reconstructions
figure('Name',sprintf('Hydrostatic full-depth partial reconstructions: %s', string(stratificationCase)))
tiledlayout(2,3,TileSpacing="compact",Padding="compact")

nexttile
plotFullDepthReconstruction(modeData.zPlot,uPlot,uReconstructionPlot,windowBounds, ...
    "F reconstruction",windowCapturedF,fullCapturedF)

nexttile
plotFullDepthReconstruction(modeData.zPlot,uPlot,uFromGReconstructionPlot,windowBounds, ...
    "F-from-G reconstruction",windowCapturedFfromG,fullCapturedFfromG)

nexttile
plotFullDepthReconstruction(modeData.zPlot,uPlot,uFromGIntSignedReconstructionPlot,windowBounds, ...
    "F-from-signed-G reconstruction",windowCapturedFfromGIntSigned,fullCapturedFfromGIntSigned)

nexttile
plotFullDepthReconstruction(modeData.zPlot,etaPlot,etaReconstructionPlot,windowBounds, ...
    "G reconstruction",windowCapturedG,fullCapturedG)

nexttile
plotFullDepthReconstruction(modeData.zPlot,etaPlot,etaInteriorSignedReconstructionPlot,windowBounds, ...
    "G signed interior reconstruction",windowCapturedGIntSigned,fullCapturedGIntSigned)

nexttile
plotFullDepthReconstruction(modeData.zPlot,etaPlot,etaInteriorAbsEtaReconstructionPlot,windowBounds, ...
    "G abs-eta interior reconstruction",windowCapturedGIntAbsEta,fullCapturedGIntAbsEta)

%% Figure 5: residuals
figure('Name',sprintf('Hydrostatic residuals: %s', string(stratificationCase)))
tiledlayout(2,3,TileSpacing="compact",Padding="compact")

nexttile
plotResidualProfile(modeData.zPlot,uResidualPlot,windowBounds, ...
    sprintf('F residual\nwindow %.3e (full %.3e)', windowResidualEnergyF, residualEnergyF))

nexttile
plotResidualProfile(modeData.zPlot,uFromGResidualPlot,windowBounds, ...
    sprintf('F-from-G residual\nwindow %.3e (full %.3e)', windowResidualEnergyFfromG, residualEnergyFfromG))

nexttile
plotResidualProfile(modeData.zPlot,uFromGIntSignedResidualPlot,windowBounds, ...
    sprintf('F-from-signed-G residual\nwindow %.3e (full %.3e)', ...
    windowResidualEnergyFfromGIntSigned, residualEnergyFfromGIntSigned))

nexttile
plotResidualProfile(modeData.zPlot,etaResidualPlot,windowBounds, ...
    sprintf('G residual\nwindow %.3e (full %.3e)', windowResidualEnergyG, residualEnergyG))

nexttile
plotResidualProfile(modeData.zPlot,etaInteriorSignedResidualPlot,windowBounds, ...
    sprintf('G signed interior residual\nwindow %.3e (full %.3e)', windowResidualEnergyGIntSigned, residualEnergyGIntSigned))

nexttile
plotResidualProfile(modeData.zPlot,etaInteriorAbsEtaResidualPlot,windowBounds, ...
    sprintf('G abs-eta interior residual\nwindow %.3e (full %.3e)', windowResidualEnergyGIntAbsEta, residualEnergyGIntAbsEta))

%% Figure 6: derivative residuals
figure('Name',sprintf('Hydrostatic derivative residuals: %s', string(stratificationCase)))
tiledlayout(2,3,TileSpacing="compact",Padding="compact")

nexttile
plotResidualProfile(modeData.zPlot,duResidualDz,windowBounds, ...
    sprintf('F derivative residual\nwindow %.3e (full %.3e)', windowResidualDerivativeEnergyF, residualDerivativeEnergyF))

nexttile
plotResidualProfile(modeData.zPlot,duFromGResidualDz,windowBounds, ...
    sprintf('F-from-G derivative residual\nwindow %.3e (full %.3e)', ...
    windowResidualDerivativeEnergyFfromG, residualDerivativeEnergyFfromG))

nexttile
plotResidualProfile(modeData.zPlot,duFromGIntSignedResidualDz,windowBounds, ...
    sprintf('F-from-signed-G derivative\nwindow %.3e (full %.3e)', ...
    windowResidualDerivativeEnergyFfromGIntSigned, residualDerivativeEnergyFfromGIntSigned))

nexttile
plotResidualProfile(modeData.zPlot,detaResidualDz,windowBounds, ...
    sprintf('G derivative residual\nwindow %.3e (full %.3e)', windowResidualDerivativeEnergyG, residualDerivativeEnergyG))

nexttile
plotResidualProfile(modeData.zPlot,detaInteriorSignedResidualDz,windowBounds, ...
    sprintf('G signed derivative residual\nwindow %.3e (full %.3e)', ...
    windowResidualDerivativeEnergyGIntSigned, residualDerivativeEnergyGIntSigned))

nexttile
plotResidualProfile(modeData.zPlot,detaInteriorAbsEtaResidualDz,windowBounds, ...
    sprintf('G abs-eta derivative residual\nwindow %.3e (full %.3e)', ...
    windowResidualDerivativeEnergyGIntAbsEta, residualDerivativeEnergyGIntAbsEta))

%% Figure 7: residual convergence
figure('Name',sprintf('Hydrostatic residual convergence: %s', string(stratificationCase)))
tiledlayout(2,2,TileSpacing="compact",Padding="compact")

nexttile
plotConvergenceCurves({convergenceF.modeCount, convergenceF.modeCount, ...
    convergenceFfromG.modeCount, convergenceFfromG.modeCount, ...
    convergenceFfromGIntSigned.modeCount, convergenceFfromGIntSigned.modeCount}, ...
    {convergenceF.profileResidualWindow, convergenceF.profileResidualFull, ...
    convergenceFfromG.profileResidualWindow, convergenceFfromG.profileResidualFull, ...
    convergenceFfromGIntSigned.profileResidualWindow, convergenceFfromGIntSigned.profileResidualFull}, ...
    ["F window", "F window", "F from G", "F from G", "F from signed G", "F from signed G"], ...
    ["-", "--", "-", "--", "-", "--"], ...
    "F profile residual", "relative int r_u^2 dz")

nexttile
plotConvergenceCurves( ...
    {convergenceG.modeCount, convergenceG.modeCount, convergenceGIntSigned.modeCount, ...
    convergenceGIntSigned.modeCount, convergenceGIntAbsEta.modeCount, convergenceGIntAbsEta.modeCount}, ...
    {convergenceG.profileResidualWindow, convergenceG.profileResidualFull, ...
    convergenceGIntSigned.profileResidualWindow, convergenceGIntSigned.profileResidualFull, ...
    convergenceGIntAbsEta.profileResidualWindow, convergenceGIntAbsEta.profileResidualFull}, ...
    ["standard", "standard", "signed", "signed", "abs-eta", "abs-eta"], ...
    ["-", "--", "-", "--", "-", "--"], ...
    "G profile residual", "relative int N^2 r_\eta^2 dz")

nexttile
plotConvergenceCurves({convergenceF.modeCount, convergenceF.modeCount, ...
    convergenceFfromG.modeCount, convergenceFfromG.modeCount, ...
    convergenceFfromGIntSigned.modeCount, convergenceFfromGIntSigned.modeCount}, ...
    {convergenceF.derivativeResidualWindow, convergenceF.derivativeResidualFull, ...
    convergenceFfromG.derivativeResidualWindow, convergenceFfromG.derivativeResidualFull, ...
    convergenceFfromGIntSigned.derivativeResidualWindow, convergenceFfromGIntSigned.derivativeResidualFull}, ...
    ["F window", "F window", "F from G", "F from G", "F from signed G", "F from signed G"], ...
    ["-", "--", "-", "--", "-", "--"], ...
    "F derivative residual", "relative g int r_{u,z}^2/N^2 dz")

nexttile
plotConvergenceCurves( ...
    {convergenceG.modeCount, convergenceG.modeCount, convergenceGIntSigned.modeCount, ...
    convergenceGIntSigned.modeCount, convergenceGIntAbsEta.modeCount, convergenceGIntAbsEta.modeCount}, ...
    {convergenceG.derivativeResidualWindow, convergenceG.derivativeResidualFull, ...
    convergenceGIntSigned.derivativeResidualWindow, convergenceGIntSigned.derivativeResidualFull, ...
    convergenceGIntAbsEta.derivativeResidualWindow, convergenceGIntAbsEta.derivativeResidualFull}, ...
    ["standard", "standard", "signed", "signed", "abs-eta", "abs-eta"], ...
    ["-", "--", "-", "--", "-", "--"], ...
    "G derivative residual", "relative int r_{\eta,z}^2 dz")

drawnow

%% Local helpers
function data = buildHydrostaticModeData(stratificationCase,windowBounds,nBaroclinicModes,nPlot,nQuad,latitude,g)
caseName = string(stratificationCase);
switch caseName
    case "exponential"
        depth = 4000;
        N0 = 3*2*pi/3600;
        L_gm = 1300;
        N2Function = @(z) N0*N0*exp(2*z/L_gm);
        zDomain = [-depth 0];
        nEVP = max(256,ceil(2.5*(nBaroclinicModes + 1)));
        zPlot = linspace(zDomain(1),zDomain(2),nPlot).';
        zWindow = linspace(windowBounds(1),windowBounds(2),nQuad).';

        imPlot = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=zPlot, ...
            latitude=latitude,nEVP=nEVP,nModes=nBaroclinicModes,g=g);
        imPlot.normalization = Normalization.geostrophic;
        imPlot.upperBoundary = UpperBoundary.rigidLid;

        imWindow = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=zWindow, ...
            latitude=latitude,nEVP=nEVP,nModes=nBaroclinicModes,g=g);
        imWindow.normalization = Normalization.geostrophic;
        imWindow.upperBoundary = UpperBoundary.rigidLid;
    case "constant"
        depth = 1300;
        N0 = 5.2e-3;
        N2Function = @(z) N0*N0 + 0*z;
        zDomain = [-depth 0];
        zPlot = linspace(zDomain(1),zDomain(2),nPlot).';
        zWindow = linspace(windowBounds(1),windowBounds(2),nQuad).';

        imPlot = InternalModesConstantStratification(N0=N0,zIn=zDomain,zOut=zPlot, ...
            latitude=latitude,nModes=nBaroclinicModes,g=g);
        imPlot.upperBoundary = UpperBoundary.rigidLid;
        imPlot.normalization = Normalization.kConstant;

        imWindow = InternalModesConstantStratification(N0=N0,zIn=zDomain,zOut=zWindow, ...
            latitude=latitude,nModes=nBaroclinicModes,g=g);
        imWindow.upperBoundary = UpperBoundary.rigidLid;
        imWindow.normalization = Normalization.kConstant;
    otherwise
        error("HydrostaticWindowModeComparison:UnknownStratification", "Unknown stratificationCase ""%s"".", caseName);
end

if windowBounds(1) < zDomain(1) || windowBounds(2) > zDomain(2) || windowBounds(1) >= windowBounds(2)
    error("HydrostaticWindowModeComparison:InvalidWindow", ...
        "windowBounds must lie inside [%.3g %.3g].", zDomain(1), zDomain(2));
end

[FPlot,GPlot,h] = imPlot.modesAtFrequency(0);
[FWindow,GWindow] = imWindow.modesAtFrequency(0);

if caseName == "constant"
    geostrophicScale = sqrt((N0*N0 - imPlot.f0*imPlot.f0)/(N0*N0));
    FPlot = geostrophicScale*FPlot;
    GPlot = geostrophicScale*GPlot;
    FWindow = geostrophicScale*FWindow;
    GWindow = geostrophicScale*GWindow;
end

data.depth = depth;
data.g = g;
data.zPlot = zPlot;
data.zWindow = zWindow;
data.N2Plot = N2Function(zPlot);
data.N2Window = N2Function(zWindow);
data.h = h(1:nBaroclinicModes);
data.h = data.h(:);
data.PhiFPlot = cat(2,ones(size(zPlot)),FPlot(:,1:nBaroclinicModes));
data.PhiFWindow = cat(2,ones(size(zWindow)),FWindow(:,1:nBaroclinicModes));
data.PhiGPlot = GPlot(:,1:nBaroclinicModes);
data.PhiGWindow = GWindow(:,1:nBaroclinicModes);
end

function weights = trapzWeights(z)
weights = zeros(size(z));
weights(1) = (z(2) - z(1))/2;
weights(end) = (z(end) - z(end-1))/2;
weights(2:end-1) = (z(3:end) - z(1:end-2))/2;
end

function spectrum = targetModeSpectrum(modeNumbers,jStar)
shape = 1./(jStar*jStar + modeNumbers.^2);
spectrum = shape/sum(shape);
end

function [etaHat,etaPlot,etaWindow] = applyEtaSignConvention(etaHat,etaPlot,etaWindow,etaSignConvention)
switch string(etaSignConvention)
    case "positiveLowerWindow"
        if etaWindow(1) < 0
            [etaHat,etaPlot,etaWindow] = flipEtaSign(etaHat,etaPlot,etaWindow);
        end
    case "negativeLowerWindow"
        if etaWindow(1) > 0
            [etaHat,etaPlot,etaWindow] = flipEtaSign(etaHat,etaPlot,etaWindow);
        end
    case "flip"
        [etaHat,etaPlot,etaWindow] = flipEtaSign(etaHat,etaPlot,etaWindow);
    case "none"
    otherwise
        error("HydrostaticWindowModeComparison:UnknownEtaSignConvention", ...
            "etaSignConvention must be ""positiveLowerWindow"", ""negativeLowerWindow"", ""flip"", or ""none"".");
end
end

function [etaHat,etaPlot,etaWindow] = flipEtaSign(etaHat,etaPlot,etaWindow)
etaHat = -etaHat;
etaPlot = -etaPlot;
etaWindow = -etaWindow;
end

function index = dominantSpectrumIndices(spectrum,nModes)
[~,index] = sort(abs(spectrum),"descend");
index = index(1:min(nModes,length(index)));
end

function index = dominantSignedContributionIndices(signedContribution,nModes)
% Indefinite quadratic forms have signed terms. Retain the largest-magnitude
% terms, and report the selected sign structure separately.
[~,index] = sort(abs(signedContribution),"descend");
index = index(1:min(nModes,length(index)));
end

function index = signedModePlotIndices(lambda,nPositive)
tolerance = 100*eps(max(abs(lambda)));
positiveIndex = find(lambda > tolerance,nPositive,"first");
negativeIndex = find(lambda < -tolerance);
index = [positiveIndex(:); negativeIndex(:)];
end

function labels = signedModeLegendLabels(lambda,index)
tolerance = 100*eps(max(abs(lambda)));
labels = strings(size(index));
for iIndex = 1:length(index)
    iMode = index(iIndex);
    if lambda(iMode) < -tolerance
        labels(iIndex) = sprintf('neg %d: %.2g', iMode, lambda(iMode));
    else
        labels(iIndex) = sprintf('pos %d: %.2g', iMode, lambda(iMode));
    end
end
end

function visibility = modeVisibility(Q,lambda,gamma)
fullDepthNorm = (gamma(:).'*Q.^2).';
visibility = lambda(:)./fullDepthNorm;
end

function score = fullDepthMetricScore(Q,gamma)
score = (gamma(:).'*Q.^2).';
end

function value = offDiagonalNorm(A)
value = norm(A - diag(diag(A)),'fro');
end

function [Q,lambda] = sortedEigenvectors(A,sortDirection)
A = (A + A.')/2;
[Q,Lambda] = eig(A);
lambda = real(diag(Lambda));
[lambda,index] = sort(lambda,sortDirection);
Q = orientColumns(real(Q(:,index)));
end

function Q = orientColumns(Q)
for iMode = 1:size(Q,2)
    [~,iMax] = max(abs(Q(:,iMode)));
    Q(:,iMode) = sign(Q(iMax,iMode))*Q(:,iMode);
end
end

function energy = weightedEnergy(field,weights)
energy = sum(weights.*field.^2);
end

function derivative = fCoefficientDerivative(coefficients,modeData,g)
derivative = fCoefficientDerivativeOnGrid(coefficients,modeData.PhiGPlot,modeData.N2Plot,g);
end

function derivative = fCoefficientDerivativeOnGrid(coefficients,PhiG,N2,g)
derivative = -(N2/g).*(PhiG*coefficients(2:end));
end

function derivative = gCoefficientDerivative(coefficients,modeData)
derivative = gCoefficientDerivativeOnGrid(coefficients,modeData.PhiFPlot,modeData.h);
end

function derivative = gCoefficientDerivativeOnGrid(coefficients,PhiF,h)
derivative = PhiF(:,2:end)*(coefficients./h);
end

function energy = inverseN2DerivativeEnergy(derivative,weights,N2,g)
energy = g*sum(weights.*derivative.^2./N2);
end

function convergence = fResidualConvergence(Q,c,spectrum,fScale,uHat,uPlot,weightsFull, ...
    uWindow,weightsWindow,modeData,g,nModes,fullEnergy,windowEnergy, ...
    fullDerivativeEnergy,windowDerivativeEnergy)
modeOrder = dominantSpectrumIndices(spectrum,nModes);
convergence.modeCount = (1:length(modeOrder)).';
convergence.profileResidualFull = zeros(length(modeOrder),1);
convergence.profileResidualWindow = zeros(length(modeOrder),1);
convergence.derivativeResidualFull = zeros(length(modeOrder),1);
convergence.derivativeResidualWindow = zeros(length(modeOrder),1);
for iMode = 1:length(modeOrder)
    index = modeOrder(1:iMode);
    reconstructionHat = fScale.*(Q(:,index)*c(index));
    residualHat = uHat - reconstructionHat;
    residualFull = uPlot - modeData.PhiFPlot*reconstructionHat;
    residualWindow = uWindow - modeData.PhiFWindow*reconstructionHat;
    derivativeResidualFull = fCoefficientDerivative(residualHat,modeData,g);
    derivativeResidualWindow = fCoefficientDerivativeOnGrid(residualHat, ...
        modeData.PhiGWindow,modeData.N2Window,g);
    convergence.profileResidualFull(iMode) = weightedEnergy(residualFull,weightsFull)/max(fullEnergy,eps);
    convergence.profileResidualWindow(iMode) = weightedEnergy(residualWindow,weightsWindow)/max(windowEnergy,eps);
    convergence.derivativeResidualFull(iMode) = inverseN2DerivativeEnergy(derivativeResidualFull, ...
        weightsFull,modeData.N2Plot,g)/max(fullDerivativeEnergy,eps);
    convergence.derivativeResidualWindow(iMode) = inverseN2DerivativeEnergy(derivativeResidualWindow, ...
        weightsWindow,modeData.N2Window,g)/max(windowDerivativeEnergy,eps);
end
end

function convergence = fConstructionOrderConvergence(Q,fScale,uHat,uPlot,weightsFull, ...
    uWindow,weightsWindow,modeData,g,nModes,fullEnergy,windowEnergy, ...
    fullDerivativeEnergy,windowDerivativeEnergy)
nModes = min(nModes,size(Q,2));
convergence.modeCount = (1:nModes).';
convergence.profileResidualFull = zeros(nModes,1);
convergence.profileResidualWindow = zeros(nModes,1);
convergence.derivativeResidualFull = zeros(nModes,1);
convergence.derivativeResidualWindow = zeros(nModes,1);
for iMode = 1:nModes
    index = (1:iMode).';
    reconstructionHat = fWindowGalerkinReconstruction(Q,index,fScale,uWindow,weightsWindow,modeData);
    residualHat = uHat - reconstructionHat;
    residualFull = uPlot - modeData.PhiFPlot*reconstructionHat;
    residualWindow = uWindow - modeData.PhiFWindow*reconstructionHat;
    derivativeResidualFull = fCoefficientDerivative(residualHat,modeData,g);
    derivativeResidualWindow = fCoefficientDerivativeOnGrid(residualHat, ...
        modeData.PhiGWindow,modeData.N2Window,g);
    convergence.profileResidualFull(iMode) = weightedEnergy(residualFull,weightsFull)/max(fullEnergy,eps);
    convergence.profileResidualWindow(iMode) = weightedEnergy(residualWindow,weightsWindow)/max(windowEnergy,eps);
    convergence.derivativeResidualFull(iMode) = inverseN2DerivativeEnergy(derivativeResidualFull, ...
        weightsFull,modeData.N2Plot,g)/max(fullDerivativeEnergy,eps);
    convergence.derivativeResidualWindow(iMode) = inverseN2DerivativeEnergy(derivativeResidualWindow, ...
        weightsWindow,modeData.N2Window,g)/max(windowDerivativeEnergy,eps);
end
end

function reconstructionHat = fWindowGalerkinReconstruction(Q,index,fScale,uWindow,weightsWindow,modeData)
QRetained = Q(:,index);
basisWindow = modeData.PhiFWindow*(fScale.*QRetained);
gram = basisWindow.'*(weightsWindow.*basisWindow);
rhs = basisWindow.'*(weightsWindow.*uWindow);
tolerance = max(size(gram))*eps(norm(gram,2));
coefficients = pinv(gram,tolerance)*rhs;
reconstructionHat = fScale.*(QRetained*coefficients);
end

function convergence = gResidualConvergence(Q,c,spectrum,gScale,etaHat,etaPlot,weightsPEFull, ...
    etaWindow,weightsPEWindow,weightsFFull,weightsFWindow,modeData,nModes,fullEnergy, ...
    windowEnergy,fullDerivativeEnergy,windowDerivativeEnergy)
modeOrder = dominantSpectrumIndices(spectrum,nModes);
convergence.modeCount = (1:length(modeOrder)).';
convergence.profileResidualFull = zeros(length(modeOrder),1);
convergence.profileResidualWindow = zeros(length(modeOrder),1);
convergence.derivativeResidualFull = zeros(length(modeOrder),1);
convergence.derivativeResidualWindow = zeros(length(modeOrder),1);
for iMode = 1:length(modeOrder)
    index = modeOrder(1:iMode);
    reconstructionHat = gScale.*(Q(:,index)*c(index));
    residualHat = etaHat - reconstructionHat;
    residualFull = etaPlot - modeData.PhiGPlot*reconstructionHat;
    residualWindow = etaWindow - modeData.PhiGWindow*reconstructionHat;
    derivativeResidualFull = gCoefficientDerivative(residualHat,modeData);
    derivativeResidualWindow = gCoefficientDerivativeOnGrid(residualHat,modeData.PhiFWindow,modeData.h);
    convergence.profileResidualFull(iMode) = weightedEnergy(residualFull,weightsPEFull)/max(fullEnergy,eps);
    convergence.profileResidualWindow(iMode) = weightedEnergy(residualWindow,weightsPEWindow)/max(windowEnergy,eps);
    convergence.derivativeResidualFull(iMode) = weightedEnergy(derivativeResidualFull,weightsFFull)/max(fullDerivativeEnergy,eps);
    convergence.derivativeResidualWindow(iMode) = weightedEnergy(derivativeResidualWindow,weightsFWindow)/max(windowDerivativeEnergy,eps);
end
end

function fraction = capturedEnergy(residualEnergy,totalEnergy)
fraction = 1 - residualEnergy/max(totalEnergy,eps);
end

function printEnergyCheck(label,spatialEnergy,quadraticEnergy,spectralEnergy)
fprintf('  %-32s spatial/quadratic %.6e, modal quadratic %.6e, spectral %.6e, rel err %.3e\n', ...
    label, spatialEnergy, quadraticEnergy, spectralEnergy, relativeDifference(spatialEnergy,spectralEnergy));
end

function printReconstructionCheck(label,index,windowCaptured,fullCaptured)
fprintf('  %-18s %2d modes, window captured %.3f, full-depth captured %.3f, indices %s\n', ...
    label, length(index), windowCaptured, fullCaptured, mat2str(index(:).'));
end

function printResidualEnergy(label,residualEnergy,totalEnergy,quantity)
fprintf('  %-30s %-20s %.6e, rel %.3e\n', ...
    label, quantity, residualEnergy, residualEnergy/max(totalEnergy,eps));
end

function printInteriorEigenvalueDiagnostics(label,lambda)
negativityTolerance = -100*eps(max(abs(lambda)));
nNegative = nnz(lambda < negativityTolerance);
fprintf('  %-22s lambda range [%.4e %.4e], negative count %d\n', ...
    label, min(lambda), max(lambda), nNegative);
if nNegative > 0
    fprintf('    %s has signed boundary-correction directions; negative eigenvalues are expected for this form.\n', label);
end
end

function printSignedSelectionDiagnostics(label,lambda,signedContribution,index)
negativityTolerance = -100*eps(max(abs(lambda)));
isNegativeDirection = lambda < negativityTolerance;
nNegativeTotal = nnz(isNegativeDirection);
nNegativeSelected = nnz(isNegativeDirection(index));
absoluteTotal = sum(abs(signedContribution));
absoluteSelected = sum(abs(signedContribution(index)));
signedSelected = sum(signedContribution(index));
fprintf(['  %-22s selected %d negative signed directions out of %d total; ' ...
    '|selected|/|total| %.3f; signed selected %.6e\n'], ...
    label,nNegativeSelected,nNegativeTotal,absoluteSelected/max(absoluteTotal,eps),signedSelected);
end

function printPlotOrder(label,index,score)
fprintf('  %-18s indices %s, full-depth score %s\n', ...
    label, mat2str(index(:).'), mat2str(score(index).',4));
end

function value = relativeDifference(a,b)
value = abs(a - b)/max([abs(a),abs(b),eps]);
end

function limits = symmetricProfileLimits(profiles,padding)
if isempty(profiles)
    limits = [-1 1];
    return
end
maxAmplitude = max(abs(profiles(:)));
if maxAmplitude == 0
    limits = [-1 1];
else
    limits = padding*maxAmplitude*[-1 1];
end
end

function plotProfilesWithWindow(profiles,z,windowBounds,plotTitle,legendLabels,xLimits)
if nargin < 6 || isempty(xLimits)
    shouldUseProfileLimits = true;
else
    shouldUseProfileLimits = false;
end
if isempty(profiles)
    xLimits = [-1 1];
elseif shouldUseProfileLimits
    xLimits = [min(profiles(:)) max(profiles(:))];
end
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none",HandleVisibility="off")
hold on
if isempty(profiles)
    text(0,mean(windowBounds),"no modes retained",HorizontalAlignment="center")
else
    lineHandles = plot(profiles,z,LineWidth=1.1);
    legend(lineHandles,legendLabels,Location="best")
end
xlabel('mode amplitude')
ylabel('z (m)')
title(plotTitle)
xlim(xLimits)
grid on
end

function plotProfileWithWindow(profile,z,windowBounds,plotTitle)
xLimits = [min(profile(:)) max(profile(:))];
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none")
hold on
plot(profile,z,'k',LineWidth=1.4)
xlabel('profile value')
ylabel('z (m)')
title(plotTitle)
grid on
end

function plotWindowSpectrum(spectrum,nPlot,plotTitle,totalEnergy,spectralEnergy,selectedIndex)
iPlot = 1:min(nPlot,length(spectrum));
stem(iPlot,spectrum(iPlot),'filled')
hold on
selectedIndex = selectedIndex(selectedIndex <= max(iPlot));
if ~isempty(selectedIndex)
    stem(selectedIndex,spectrum(selectedIndex),'r','filled',DisplayName='reconstruction modes')
end
xlabel('window mode index')
ylabel('spectrum')
title(sprintf('%s\ntotal %.3e, spectral %.3e, rel %.2e', ...
    plotTitle, totalEnergy, spectralEnergy, relativeDifference(totalEnergy,spectralEnergy)))
grid on
end

function plotFullDepthReconstruction(z,trueProfile,reconstruction,windowBounds,plotTitle,windowCaptured,fullCaptured)
xLimits = [min([trueProfile(:); reconstruction(:)]) max([trueProfile(:); reconstruction(:)])];
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none",DisplayName='window')
hold on
plot(trueProfile,z,'k',LineWidth=1.4,DisplayName='truth')
plot(reconstruction,z,'r--',LineWidth=1.2,DisplayName='partial reconstruction')
xlabel('profile value')
ylabel('z (m)')
title(sprintf('%s\nwindow captured %.2f, full captured %.2f', plotTitle, windowCaptured, fullCaptured))
legend(Location="best")
grid on
end

function plotResidualProfile(z,residual,windowBounds,plotTitle)
xLimits = [min(residual(:)) max(residual(:))];
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none",DisplayName='window')
hold on
plot(residual,z,'k',LineWidth=1.2,DisplayName='residual')
xline(0,'Color',[0.45 0.45 0.45],LineStyle=':')
xlabel('residual')
ylabel('z (m)')
title(plotTitle)
legend(Location="best")
grid on
end

function plotConvergenceCurves(modeCounts,residuals,labels,lineStyles,plotTitle,yAxisLabel)
hold on
uniqueLabels = unique(labels,"stable");
colors = lines(numel(uniqueLabels));
for iCurve = 1:length(residuals)
    colorIndex = find(uniqueLabels == labels(iCurve),1);
    color = colors(colorIndex,:);
    legendLabel = labels(iCurve) + " " + lineStyleLabel(lineStyles(iCurve));
    semilogy(modeCounts{iCurve},max(residuals{iCurve},eps),LineWidth=1.2, ...
        Color=color,LineStyle=char(lineStyles(iCurve)),DisplayName=legendLabel)
end
hold off
xlabel('number of modes included')
ylabel(yAxisLabel)
title(plotTitle)
legend(Location="best")
grid on
end

function label = lineStyleLabel(lineStyle)
if lineStyle == "-"
    label = "window";
else
    label = "full depth";
end
end
