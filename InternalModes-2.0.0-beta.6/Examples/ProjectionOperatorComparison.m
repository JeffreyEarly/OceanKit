%% Compare direct and weighted projection operators for vertical modes
% This example follows the projection-operator experiment in
% wave-vortex-model/UnitTests/MiscTests/QuadratureWeightsReRevisited.m, but
% only uses the InternalModes package.
%
% The main point is that the weighted-adjoint operator Phi' W is only a
% forward transform when the sampled modes are exactly orthonormal under
% the discrete weights. The weighted pseudoinverse (Phi' W Phi)\Phi' W and
% the direct inverse agree on a square, full-rank native grid, while the
% reduced weighted pseudoinverse gives the best weighted fit to a smaller
% set of retained modes.

%% Compute modes on mode-adapted quadrature points
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
k = 0;
nPoints = 64;
nModes = nPoints - 1;
nEVP = max(256,ceil(2.1*nPoints));
zDomain = [-Lz 0];

zReference = linspace(zDomain(1),zDomain(2),1024).';
imReference = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,g=g);
imReference.normalization = Normalization.kConstant;
imReference.upperBoundary = UpperBoundary.rigidLid;

z = imReference.GaussQuadraturePointsForModesAtWavenumber(nPoints,k);
im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=latitude,nEVP=nEVP,nModes=nModes,g=g);
im.normalization = Normalization.kConstant;
im.upperBoundary = UpperBoundary.rigidLid;

[Finv,Ginv] = im.modesAtWavenumber(k);

%% Build square F and G inverse-transform matrices
Nz = size(Finv,1);

% Add the barotropic F mode so F is square. For rigid-lid G modes, remove
% the zero boundary rows and drop the highest G mode so G is square.
Finv = cat(2,ones(Nz,1),Finv);
Ginv = Ginv(2:end-1,1:end-1);

P = max(abs(Finv),[],1);
Q = max(abs(Ginv),[],1);

PFinv = Finv./P;
QGinv = Ginv./Q;

%% Compute quadrature weights from int F dz = 0
% With the barotropic mode first, the weights w satisfy
%
%     PFinv.' * w = [Lz; 0; ...; 0].
%
% This uses the fact that the barotropic F mode integrates to Lz and the
% baroclinic F modes integrate to zero.
b = zeros(Nz,1);
b(1) = Lz;
w = (PFinv.')\b;

zInterior = z(2:end-1);
weights = w(2:end-1) .* (N2(zInterior) - im.f0*im.f0)/g;

%% Compare three forward operators for the full square G basis
T = eye(size(QGinv,2));

QGDirect = inv(QGinv);
QGWeightedAdjoint = ((weights .* Ginv) .* Q).';
QGWeightedInverse = weightedPseudoinverse(QGinv, weights);

directError = relativeFrobeniusError(QGDirect*QGinv,T);
weightedAdjointError = relativeFrobeniusError(QGWeightedAdjoint*QGinv,T);
weightedInverseError = relativeFrobeniusError(QGWeightedInverse*QGinv,T);

fprintf('Full square G basis with %d interior points\n', size(QGinv,1));
fprintf('  direct inverse round-trip error:       %.3e\n', directError);
fprintf('  weighted-adjoint round-trip error:    %.3e\n', weightedAdjointError);
fprintf('  weighted-pseudoinverse round-trip:    %.3e\n', weightedInverseError);

discreteGram = Ginv.' * (weights .* Ginv);
fprintf('  unscaled discrete Gram error:         %.3e\n', relativeFrobeniusError(discreteGram,T));
fprintf('  condition number of QGinv:            %.3e\n', cond(QGinv));

%% Build a reduced-mode weighted projection operator
nReducedModes = 10;
PhiReduced = QGinv(:,1:nReducedModes);
TReduced = eye(nReducedModes);

QGReduced = weightedPseudoinverse(PhiReduced, weights);
QGDirectTruncated = QGDirect(1:nReducedModes,:);

reducedError = relativeFrobeniusError(QGReduced*PhiReduced,TReduced);
truncatedError = relativeFrobeniusError(QGDirectTruncated*PhiReduced,TReduced);

fprintf('\nReduced %d-mode operator on the same grid\n', nReducedModes);
fprintf('  reduced weighted-pseudoinverse error: %.3e\n', reducedError);
fprintf('  first rows of full inverse error:     %.3e\n', truncatedError);

smoothProfile = exp(zInterior/L_gm);
coefficientsFromReducedFit = QGReduced*smoothProfile;
coefficientsFromFullInverse = QGDirectTruncated*smoothProfile;
coefficientDifference = norm(coefficientsFromReducedFit - coefficientsFromFullInverse)/norm(coefficientsFromReducedFit);

fprintf('  coefficient difference for exp(z/L):  %.3e\n', coefficientDifference);

%% Local helpers
function A = weightedPseudoinverse(Phi, weights)
weightedPhi = weights .* Phi;
A = (Phi.' * weightedPhi) \ weightedPhi.';
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,'fro')/norm(B,'fro');
end
