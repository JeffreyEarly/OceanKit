%% Compare shared dz and effective N2 quadrature for hydrostatic modes
% This example compares several ways of building the discrete inner products
% for hydrostatic F and G modes on a mode-adapted rigid-lid grid.
%
% The F modes use a shared vertical increment dz chosen so that
%
%     PhiF.' * dz = [D; 0; ...; 0].
%
% The G modes depend on the product qG = N2_i dz_i/g. If dz is fixed by the
% F constraints, qG may still be improved by interpreting N2_i as an effective
% quadrature value rather than the pointwise physical stratification.

%% Single-resolution comparison
baseNPoints = 64;
base = computeDiagnostics(baseNPoints);

fprintf('\nEffective N2 quadrature comparison at nPoints = %d\n', baseNPoints);
fprintf('%-30s %12s %12s %12s %12s %12s %12s %12s\n', ...
    'weights','F Gram','G Gram','F int','sum dz','min ratio','max ratio','rel N2eff');
printDiagnosticRow('geometric', base.geometric);
printDiagnosticRow('F-compatible', base.fCompatible);
printDiagnosticRow('Christoffel dz', base.christoffelDz);
printDiagnosticRow('F-compatible + effective N2', base.effectiveN2);

%% Resolution dependence of effective N2
% If N2Effective/N2(z_i) approaches one with resolution, the effective
% quadrature is mostly correcting finite-resolution sampling error. If it
% does not, the modal quadrature rule remains distinct from pointwise
% stratification sampling.
nPointsList = [24 32 48 64 96 128];

fprintf('\nResolution dependence of effective N2\n');
fprintf('%8s %10s %18s %18s %12s %12s %12s %12s\n', ...
    'nPoints','nInterior','G Gram physical','G Gram effective','rel N2eff','min ratio','max ratio','mean ratio');
for iPoint = 1:length(nPointsList)
    diagnostics = computeDiagnostics(nPointsList(iPoint));
    row = diagnostics.effectiveN2;
    fprintf('%8d %10d %18.6e %18.6e %12.6e %12.6e %12.6e %12.6e\n', ...
        diagnostics.nPoints,diagnostics.nInterior,diagnostics.fCompatible.gGramError, ...
        row.gGramError,row.relativeN2Difference,row.minN2Ratio,row.maxN2Ratio,row.meanN2Ratio);
end

%% Local helpers
function diagnostics = computeDiagnostics(nPoints)
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
omega = 0;
nModes = nPoints - 1;
nEVP = max(256,ceil(2.1*nPoints));
zDomain = [-Lz 0];

zReference = linspace(zDomain(1),zDomain(2),1024).';
imReference = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,g=g);
imReference.normalization = Normalization.kConstant;
imReference.upperBoundary = UpperBoundary.rigidLid;

z = imReference.GaussQuadraturePointsForModesAtFrequency(nPoints,omega);
im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=latitude,nEVP=nEVP,nModes=nModes,g=g);
im.normalization = Normalization.kConstant;
im.upperBoundary = UpperBoundary.rigidLid;

[F,G,h] = im.modesAtFrequency(omega);

Nz = length(z);
PhiF = cat(2,ones(Nz,1),F);
GammaF = diag([Lz; h(:)]);
bF = zeros(Nz,1);
bF(1) = Lz;

dzF = (PhiF.')\bF;
dzGeometric = geometricIncrements(z);

zInterior = z(2:end-1);
N2Physical = N2(zInterior);
Gactive = G(2:end-1,1:end-1);
nInterior = size(Gactive,1);

qGChristoffel = 1./sum(Gactive.^2,2);
dzChristoffel = zeros(Nz,1);
dzChristoffel(2:end-1) = g*qGChristoffel./N2Physical;
missingDepth = Lz - sum(dzChristoffel);
dzChristoffel(1) = missingDepth/2;
dzChristoffel(end) = missingDepth/2;

N2Effective = g*qGChristoffel./dzF(2:end-1);
ratio = N2Effective./N2Physical;

diagnostics.nPoints = nPoints;
diagnostics.nInterior = nInterior;
diagnostics.geometric = rowDiagnostics(PhiF,Gactive,GammaF,bF,dzGeometric,N2Physical,N2Physical,g,[]);
diagnostics.fCompatible = rowDiagnostics(PhiF,Gactive,GammaF,bF,dzF,N2Physical,N2Physical,g,[]);
diagnostics.christoffelDz = rowDiagnostics(PhiF,Gactive,GammaF,bF,dzChristoffel,N2Physical,N2Physical,g,[]);
diagnostics.effectiveN2 = rowDiagnostics(PhiF,Gactive,GammaF,bF,dzF,N2Effective,N2Physical,g,ratio);
end

function row = rowDiagnostics(PhiF,Gactive,GammaF,bF,dz,N2Quadrature,N2Physical,g,n2Ratio)
qG = dz(2:end-1).*N2Quadrature/g;
targetG = eye(size(Gactive,2));

row.fGramError = relativeFrobeniusError(PhiF.'*(dz.*PhiF),GammaF);
row.gGramError = relativeFrobeniusError(Gactive.'*(qG.*Gactive),targetG);
row.fIntegralError = norm(PhiF.'*dz - bF)/norm(bF);
row.sumDz = sum(dz);

if isempty(n2Ratio)
    row.minN2Ratio = NaN;
    row.maxN2Ratio = NaN;
    row.meanN2Ratio = NaN;
    row.relativeN2Difference = NaN;
else
    row.minN2Ratio = min(n2Ratio);
    row.maxN2Ratio = max(n2Ratio);
    row.meanN2Ratio = mean(n2Ratio);
    row.relativeN2Difference = norm(N2Quadrature - N2Physical)/norm(N2Physical);
end
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2)-z(1))/2,(z(3:end)-z(1:end-2))/2,(z(end)-z(end-1))/2);
end

function printDiagnosticRow(name,row)
fprintf('%-30s %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n', ...
    name,row.fGramError,row.gGramError,row.fIntegralError,row.sumDz, ...
    row.minN2Ratio,row.maxN2Ratio,row.relativeN2Difference);
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,'fro')/norm(B,'fro');
end
