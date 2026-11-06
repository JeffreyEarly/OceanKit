%% Antialias k-dependent IGW modes on a fixed hydrostatic quadrature grid
% This example mimics the vertical-mode setup used by the Boussinesq
% wave-vortex model without calling wave-vortex-model code. The hydrostatic
% omega=0 modes define a single vertical quadrature grid. Nonzero-kappa IGW
% modes are then evaluated on that fixed grid, even though each kappa would
% prefer its own mode-adapted quadrature points.
%
% The important diagnostic is not whether the retained columns can be
% inverted on the fixed grid. A weighted pseudoinverse makes that true by
% construction. The useful diagnostic is whether excluded higher modes leak
% into the retained coefficients. When the leakage is too large, those modes
% should be discarded from the forward transform.

%% Use the same stratification as ProjectionOperatorComparison
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
Nz = 64;
Nj = floor(2*(Nz - 1)/3);
nInternalTarget = Nj - 1;
nEVP = max(256,ceil(2.1*Nz));
zDomain = [-Lz 0];

%% Build the hydrostatic quadrature grid and projections
zReference = linspace(zDomain(1),zDomain(2),10*Nz).';
imReference = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,g=g);
imReference.normalization = Normalization.geostrophic;
imReference.upperBoundary = UpperBoundary.rigidLid;
[~,~,~,zHydro] = imReference.modesAtQuadraturePoints(nPoints=Nz,omega=0);

imHydro = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zHydro,latitude=latitude,nEVP=nEVP,nModes=Nz-1,g=g);
imHydro.normalization = Normalization.geostrophic;
imHydro.upperBoundary = UpperBoundary.rigidLid;

[FinvHydro,GinvHydro,hHydro] = imHydro.modesAtFrequency(0);
[~,~,~,PF0,QGinv0,QG0,~,wHydro] = rigidLidProjectionOperators(FinvHydro,GinvHydro,hHydro,Nj,Lz);

fprintf('Hydrostatic grid: Nz=%d, WVM-style Nj=%d (%d internal G modes)\n', Nz, Nj, nInternalTarget);
fprintf('  size(PF0)=[%d %d], size(QG0)=[%d %d]\n', size(PF0,1), size(PF0,2), size(QG0,1), size(QG0,2));

%% Scan nonzero kappas on the fixed hydrostatic grid
imWave = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zHydro,latitude=latitude,nEVP=nEVP,nModes=Nz-1,g=g);
imWave.normalization = Normalization.kConstant;
imWave.upperBoundary = UpperBoundary.rigidLid;

gWeights = wHydro(2:end-1) .* (N2(zHydro(2:end-1)) - imWave.f0*imWave.f0)/g;
aliasTolerance = 1e-2;
maxConditionNumber = 1e6;
nTailCheck = 8;

wavelength = (160:-5:80).';
kappa = 2*pi./wavelength;
nResolved = zeros(size(kappa));
maxTailLeakage = zeros(size(kappa));
gramConditionNumber = zeros(size(kappa));

for iK = 1:length(kappa)
    [~,GinvWave] = imWave.modesAtWavenumber(kappa(iK));
    Phi = normalizedInteriorGModes(GinvWave);
    [nResolved(iK),maxTailLeakage(iK),gramConditionNumber(iK)] = resolvedModeCount(Phi, gWeights, nInternalTarget, nTailCheck, aliasTolerance, maxConditionNumber);
end

targetResolved = round(nInternalTarget/2);
[~,iSelected] = min(abs(nResolved - targetResolved));
kappaSelected = kappa(iSelected);
wavelengthSelected = wavelength(iSelected);

scanTable = table(wavelength,kappa,nResolved,maxTailLeakage,gramConditionNumber);
disp(scanTable);

fprintf('Selected kappa %.3e rad/m (wavelength %.1f m): keeps %d of %d internal modes.\n', ...
    kappaSelected, wavelengthSelected, nResolved(iSelected), nInternalTarget);

%% Build the full legacy-style and antialiased wave projections
[FinvWave,GinvWave,hWave] = imWave.modesAtWavenumber(kappaSelected);
[~,QGpmInv,QGpm,~] = rigidLidGProjectionOperators(GinvWave,hWave,Nj);

Phi = normalizedInteriorGModes(GinvWave);
[nInternalResolved,acceptedLeakage,acceptedCondition,nextLeakage] = resolvedModeCount(Phi, gWeights, nInternalTarget, nTailCheck, aliasTolerance, maxConditionNumber);
[~,PFpmInvAntialiased,PFpmAntialiased] = antialiasedFProjectionOperators(FinvWave, wHydro, nInternalResolved, Nj);
[QGpmAntialiased,QGpmInvAntialiased] = antialiasedGProjectionOperators(GinvWave, gWeights, nInternalResolved, Nj);

legacyRoundTripError = relativeFrobeniusError(QGpm(2:end,:)*QGpmInv(:,2:end), eye(Nj-1));
retainedFRoundTripError = relativeFrobeniusError( ...
    PFpmAntialiased(1:nInternalResolved+1,:)*PFpmInvAntialiased(:,1:nInternalResolved+1), ...
    eye(nInternalResolved+1));
retainedRoundTripError = relativeFrobeniusError( ...
    QGpmAntialiased(2:nInternalResolved+1,:)*QGpmInvAntialiased(:,2:nInternalResolved+1), ...
    eye(nInternalResolved));

discardedPFRowNorm = max(vecnorm(PFpmAntialiased(nInternalResolved+2:end,:),2,2));
discardedQGRowNorm = max(vecnorm(QGpmAntialiased(nInternalResolved+2:end,:),2,2));
QGwgAntialiased = QGpmAntialiased*QGinv0;
discardedCrossRowNorm = max(vecnorm(QGwgAntialiased(nInternalResolved+2:end,:),2,2));

fprintf('\nWave projection at selected kappa\n');
fprintf('  full direct inverse, retained-block round-trip: %.3e\n', legacyRoundTripError);
fprintf('  antialiased F retained-mode round-trip:        %.3e\n', retainedFRoundTripError);
fprintf('  antialiased G retained-mode round-trip:        %.3e\n', retainedRoundTripError);
fprintf('  accepted tail leakage:                         %.3e\n', acceptedLeakage);
fprintf('  next-mode tail leakage:                        %.3e\n', nextLeakage);
fprintf('  accepted Gram condition number:                %.3e\n', acceptedCondition);
fprintf('  discarded PF-row norm:                         %.3e\n', discardedPFRowNorm);
fprintf('  discarded QG-row norm:                         %.3e\n', discardedQGRowNorm);
fprintf('  discarded QGwg-row norm:                       %.3e\n', discardedCrossRowNorm);

%% Local helpers
function [P,Q,PFinv,PF,QGinv,QG,h,w] = rigidLidProjectionOperators(Finv,Ginv,h,Nj,Lz)
Nz = size(Finv,1);
nModes = size(Finv,2);

FWithBarotropic = cat(2,ones(Nz,1),Finv);
GInterior = Ginv(2:end-1,1:end-1);

PAll = max(abs(FWithBarotropic),[],1);
QAll = max(abs(GInterior),[],1);

PFinvFull = FWithBarotropic./PAll;
QGinvFull = GInterior./QAll;
PFFull = inv(PFinvFull);
QGFull = inv(QGinvFull);

b = zeros(Nz,1);
b(1) = Lz;
w = (PFinvFull.')\b;

PFinv = PFinvFull(:,1:end-1);
PF = PFFull(1:end-1,:);
QGinv = cat(2,zeros(Nz,1),cat(1,zeros(1,nModes-1),QGinvFull,zeros(1,nModes-1)));
QG = cat(2,zeros(nModes,1),cat(1,zeros(1,Nz-2),QGFull),zeros(nModes,1));
h = cat(1,1,reshape(h(1:end-1),[],1));

P = reshape(PAll(1:end-1),[],1);
Q = reshape(cat(2,1,QAll),[],1);

PFinv = PFinv(:,1:Nj);
PF = PF(1:Nj,:);
P = P(1:Nj);
QGinv = QGinv(:,1:Nj);
QG = QG(1:Nj,:);
Q = Q(1:Nj);
h = h(1:Nj);
end

function Phi = normalizedInteriorGModes(Ginv)
GInterior = Ginv(2:end-1,1:end-1);
Q = max(abs(GInterior),[],1);
Phi = GInterior./Q;
end

function [Q,QGinv,QG,h] = rigidLidGProjectionOperators(Ginv,h,Nj)
Nz = size(Ginv,1);
nModes = size(Ginv,2);

GInterior = Ginv(2:end-1,1:end-1);
QAll = max(abs(GInterior),[],1);
QGinvFull = GInterior./QAll;
QGFull = inv(QGinvFull);

QGinv = cat(2,zeros(Nz,1),cat(1,zeros(1,nModes-1),QGinvFull,zeros(1,nModes-1)));
QG = cat(2,zeros(nModes,1),cat(1,zeros(1,Nz-2),QGFull),zeros(nModes,1));
h = cat(1,1,reshape(h(1:end-1),[],1));
Q = reshape(cat(2,1,QAll),[],1);

QGinv = QGinv(:,1:Nj);
QG = QG(1:Nj,:);
Q = Q(1:Nj);
h = h(1:Nj);
end

function [nResolved,maxLeakage,gramConditionNumber,nextLeakage] = resolvedModeCount(Phi, weights, nTarget, nTailCheck, aliasTolerance, maxConditionNumber)
nAvailable = min(nTarget, size(Phi,2));
iTailEnd = min(size(Phi,2), nAvailable + nTailCheck);
nResolved = 0;
maxLeakage = NaN;
gramConditionNumber = NaN;
nextLeakage = NaN;

for n = 1:nAvailable
    [candidateLeakage,candidateCondition] = projectionQuality(Phi, weights, n, iTailEnd);
    if candidateCondition > maxConditionNumber || candidateLeakage > aliasTolerance
        nextLeakage = candidateLeakage;
        break;
    end

    nResolved = n;
    maxLeakage = candidateLeakage;
    gramConditionNumber = candidateCondition;
end

if nResolved == nAvailable
    nextLeakage = 0;
elseif isnan(nextLeakage)
    [nextLeakage,~] = projectionQuality(Phi, weights, nResolved + 1, iTailEnd);
end
end

function [maxLeakage,gramConditionNumber] = projectionQuality(Phi, weights, nModes, iTailEnd)
PhiN = Phi(:,1:nModes);
weightedPhiN = weights .* PhiN;
gram = PhiN.' * weightedPhiN;
A = gram \ weightedPhiN.';
gramConditionNumber = cond(gram);

if nModes < iTailEnd
    leakage = A * Phi(:,nModes+1:iTailEnd);
    maxLeakage = max(vecnorm(leakage,2,1));
else
    maxLeakage = 0;
end
end

function [QG,QGinv] = antialiasedGProjectionOperators(Ginv, weights, nInternalResolved, Nj)
Nz = size(Ginv,1);
Phi = normalizedInteriorGModes(Ginv);
PhiResolved = Phi(:,1:nInternalResolved);
weightedPhi = weights .* PhiResolved;
GForward = (PhiResolved.' * weightedPhi) \ weightedPhi.';

QG = zeros(Nj,Nz);
QG(2:nInternalResolved+1,2:end-1) = GForward;

QGinv = zeros(Nz,Nj);
QGinv(2:end-1,2:nInternalResolved+1) = PhiResolved;
end

function [P,PFinv,PF] = antialiasedFProjectionOperators(Finv, weights, nInternalResolved, Nj)
Nz = size(Finv,1);
nRetained = nInternalResolved + 1;
FWithBarotropic = cat(2,ones(Nz,1),Finv);
PAll = max(abs(FWithBarotropic),[],1);
Phi = FWithBarotropic./PAll;
PhiRetained = Phi(:,1:nRetained);
weightedPhi = weights .* PhiRetained;
FForward = (PhiRetained.' * weightedPhi) \ weightedPhi.';

P = ones(Nj,1);
P(1:nRetained) = PAll(1:nRetained).';
PFinv = zeros(Nz,Nj);
PFinv(:,1:nRetained) = PhiRetained;
PF = zeros(Nj,Nz);
PF(1:nRetained,:) = FForward;
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,'fro')/norm(B,'fro');
end
