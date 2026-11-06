%% Hydrostatic transform numerics on a mode-adapted grid
% This example tests how well hydrostatic F and G modes preserve Parseval on
% a fixed G^{N+1} quadrature grid. The primary metric is the worst-case
% relative Parseval error over retained modal coefficients,
%
%   E_spec = || Gamma0^{-1/2} (Gamma-Gamma0) Gamma0^{-1/2} ||_2.
%
% The example also reports coefficient round-trip errors. Those can be close
% to machine precision even when the Parseval metric is not, so they are
% useful diagnostics but not the primary quadrature-quality target.

nPoints = 64;
nPointsList = [24 32 48 64 96 128];
latitude = 31;
g = 9.81;

cases = [
    hydrostaticCase("exponential",4000,3*2*pi/3600)
    hydrostaticCase("constant",1300,5.2e-3)
    ];

for iCase = 1:length(cases)
    data = modeData(cases(iCase),nPoints,latitude,g);
    rows = diagnosticRows(data);

    fprintf('\n%s stratification, nPoints = %d\n', data.name, nPoints);
    printComparisonHeader();
    for iRow = 1:length(rows)
        printComparisonRow(rows(iRow));
    end

    if data.name == "constant"
        printConstantStratificationCheck(data,rows(1).dz);
    end

    practicalRow = practicalWeightRow(rows);
    fprintf('\nBand-limited Parseval diagnostics using %s weights\n', practicalRow.dzRule);
    printBandLimitHeader();
    bandRows = bandLimitedRows(data,practicalRow.dz);
    for iRow = 1:length(bandRows)
        printBandLimitRow(bandRows(iRow));
    end
end

fprintf('\nResolution sweep on fixed G^{N+1} nodes\n');
fprintf('%-13s %8s %-18s %10s %10s %10s %10s %10s\n', ...
    'case','nPoints','dz rule','F spec','G spec','max spec','cond F','cond G');
for iCase = 1:length(cases)
    for iPoint = 1:length(nPointsList)
        data = modeData(cases(iCase),nPointsList(iPoint),latitude,g);
        rows = diagnosticRows(data);
        row = practicalWeightRow(rows);
        fprintf('%-13s %8d %-18s %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
            data.name,nPointsList(iPoint),row.dzRule,row.f.specError,row.g.specError, ...
            max(row.f.specError,row.g.specError),row.f.conditionNumber,row.g.conditionNumber);
    end
end

%% Local helpers
function definition = hydrostaticCase(name,depth,N0)
definition.name = string(name);
definition.depth = depth;
definition.N0 = N0;
end

function data = modeData(definition,nPoints,latitude,g)
omega = 0;
data.name = definition.name;
data.depth = definition.depth;
data.g = g;
data.nPoints = nPoints;

if definition.name == "exponential"
    L_gm = 1300;
    N2Function = @(z) definition.N0*definition.N0*exp(2*z/L_gm);
    nEVP = max(256,ceil(2.1*nPoints));
    zDomain = [-definition.depth 0];
    zReference = linspace(zDomain(1),zDomain(2),1024).';
    imReference = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,g=g);
    imReference.normalization = Normalization.geostrophic;
    imReference.upperBoundary = UpperBoundary.rigidLid;
    z = imReference.GaussQuadraturePointsForModesAtFrequency(nPoints,omega);

    im = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=z,latitude=latitude,nEVP=nEVP,nModes=nPoints-1,g=g);
    im.normalization = Normalization.geostrophic;
    im.upperBoundary = UpperBoundary.rigidLid;
elseif definition.name == "constant"
    z = linspace(-definition.depth,0,nPoints).';
    im = InternalModesConstantStratification(N0=definition.N0,zIn=[-definition.depth 0],zOut=z,latitude=latitude,nModes=nPoints-1,g=g);
    im.upperBoundary = UpperBoundary.rigidLid;
    im.normalization = Normalization.kConstant;
else
    error("HydrostaticTransformNumerics:UnknownCase", "Unknown case %s.", definition.name);
end

[F,G,h] = im.modesAtFrequency(omega);
if definition.name == "constant"
    % The analytical constant-stratification class does not expose
    % Normalization.geostrophic, but for omega=0 it differs from kConstant by
    % this scalar factor. After rescaling, the physical hydrostatic targets
    % int F^2 dz = h and (1/g) int N2 G^2 dz = 1 are recovered.
    geostrophicScale = sqrt((definition.N0*definition.N0 - im.f0*im.f0)/(definition.N0*definition.N0));
    F = geostrophicScale*F;
    G = geostrophicScale*G;
end

nInterior = nPoints - 2;
data.z = z;
data.N2 = im.N2(:);
data.N2Interior = data.N2(2:end-1);
data.PhiFParseval = cat(2,ones(nPoints,1),F(:,1:nInterior));
data.PhiFInverse = cat(2,ones(nPoints,1),F(:,1:(nInterior+1)));
data.PhiG = G(2:end-1,1:nInterior);
data.GammaF0 = diag([definition.depth; h(1:nInterior).']);
data.GammaF0Inverse = diag([definition.depth; h(1:(nInterior+1)).']);
data.GammaG0 = eye(nInterior);
data.bF = zeros(size(data.PhiFParseval,2),1);
data.bF(1) = definition.depth;
end

function rows = diagnosticRows(data)
rows = diagnosticRow(data,"geometric",geometricIncrements(data.z));
rows(end+1) = diagnosticRow(data,"F-compatible",leastSquaresSolution(data.PhiFParseval.',data.bF));

dz = jointLeastSquaresIncrements(data);
if ~isempty(dz)
    rows(end+1) = diagnosticRow(data,"joint LS",dz);
end
end

function row = diagnosticRow(data,dzRule,dz)
row.dzRule = string(dzRule);
row.dz = dz;
weightsF = dz;
weightsG = dz(2:end-1).*data.N2Interior/data.g;

row.f = componentDiagnostics(data.PhiFParseval,weightsF,data.GammaF0);
row.g = componentDiagnostics(data.PhiG,weightsG,data.GammaG0);
row.fInverseError = directInverseError(data.PhiFInverse,data.GammaF0Inverse);
row.gInverseError = directInverseError(data.PhiG,data.GammaG0);
row.sumDz = sum(dz);
row.minDz = min(dz);
end

function diagnostics = componentDiagnostics(Phi,weights,Gamma0)
Gamma = Phi.'*(weights.*Phi);
modeCount = size(Phi,2);
normalizedDifference = normalizedMatrix(Gamma - Gamma0,Gamma0);
normalizedGram = normalizedMatrix(Gamma,Gamma0);
offdiag = normalizedGram - diag(diag(normalizedGram));

diagnostics.specError = norm(normalizedDifference,2);
diagnostics.frobeniusError = norm(Gamma - Gamma0,'fro')/norm(Gamma0,'fro');
diagnostics.maxDiagError = max(abs(diag(normalizedGram) - 1));
diagnostics.offdiagError = norm(offdiag,2);
diagnostics.adjointError = diagnostics.specError;
diagnostics.galerkinError = normalizedOperatorError(Gamma \ Gamma,Gamma0);
diagnostics.conditionNumber = cond(Phi);
diagnostics.modeCount = modeCount;
end

function value = directInverseError(Phi,Gamma0)
if size(Phi,1) ~= size(Phi,2)
    value = NaN;
    return
end

value = normalizedOperatorError(Phi\Phi,Gamma0);
end

function value = normalizedOperatorError(operatorProduct,Gamma0)
scale = diag(sqrt(diag(Gamma0)));
inverseScale = diag(1./sqrt(diag(Gamma0)));
value = norm(scale*(operatorProduct - eye(size(operatorProduct)))*inverseScale,2);
end

function value = normalizedMatrix(matrix,Gamma0)
scale = diag(1./sqrt(diag(Gamma0)));
value = scale*matrix*scale;
end

function row = practicalWeightRow(rows)
score = zeros(size(rows));
for iRow = 1:length(rows)
    score(iRow) = max(rows(iRow).f.specError,rows(iRow).g.specError);
end
[~,index] = min(score);
row = rows(index);
end

function dz = jointLeastSquaresIncrements(data)
[A,rhs] = normalizedJointSystem(data);
Aeq = ones(1,length(data.z));
beq = data.depth;
dz = equalityConstrainedLeastSquares(A,rhs,Aeq,beq);
if all(dz > -1e-10)
    dz(dz < 0) = 0;
    return
end

try
    options = optimoptions("lsqlin",Display="off");
    dz = lsqlin(A,rhs,[],[],Aeq,beq,zeros(length(data.z),1),[],[],options);
catch
    dz = [];
end
end

function x = equalityConstrainedLeastSquares(A,b,Aeq,beq)
x0 = Aeq.'*((Aeq*Aeq.')\beq);
Z = null(Aeq);
x = x0 + Z*leastSquaresSolution(A*Z,b - A*x0);
end

function [A,rhs] = normalizedJointSystem(data)
nGrid = length(data.z);
nF = size(data.PhiFParseval,2);
nG = size(data.PhiG,2);
nRows = nF*nF + nG*nG;
A = zeros(nRows,nGrid);
rhs = zeros(nRows,1);
row = 0;

fScale = sqrt(diag(data.GammaF0));
for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        A(row,:) = (data.PhiFParseval(:,iMode).*data.PhiFParseval(:,jMode)./(fScale(iMode)*fScale(jMode))).';
        rhs(row) = double(iMode == jMode);
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(nGrid,1);
        weights(2:end-1) = (data.N2Interior/data.g).*data.PhiG(:,iMode).*data.PhiG(:,jMode);
        A(row,:) = weights.';
        rhs(row) = double(iMode == jMode);
    end
end
end

function rows = bandLimitedRows(data,dz)
maxModes = min(size(data.PhiFParseval,2)-1,size(data.PhiG,2));
modeCounts = unique([8 16 24 32 48 maxModes]);
modeCounts = modeCounts(modeCounts <= maxModes);
rows = repmat(struct("nModes",0,"fSpecError",0,"gSpecError",0,"maxSpecError",0,"fConditionNumber",0,"gConditionNumber",0),size(modeCounts));

for iCount = 1:length(modeCounts)
    nModes = modeCounts(iCount);
    phiF = data.PhiFParseval(:,1:(nModes+1));
    gammaF = data.GammaF0(1:(nModes+1),1:(nModes+1));
    phiG = data.PhiG(:,1:nModes);
    gammaG = data.GammaG0(1:nModes,1:nModes);

    fDiagnostics = componentDiagnostics(phiF,dz,gammaF);
    gDiagnostics = componentDiagnostics(phiG,dz(2:end-1).*data.N2Interior/data.g,gammaG);
    rows(iCount).nModes = nModes;
    rows(iCount).fSpecError = fDiagnostics.specError;
    rows(iCount).gSpecError = gDiagnostics.specError;
    rows(iCount).maxSpecError = max(fDiagnostics.specError,gDiagnostics.specError);
    rows(iCount).fConditionNumber = cond(phiF);
    rows(iCount).gConditionNumber = cond(phiG);
end
end

function printConstantStratificationCheck(data,dz)
nInterior = size(data.PhiG,2);
x = (data.z + data.depth)/data.depth;
cosineBasis = cos(x*(0:nInterior)*pi);
sineBasis = sin(x(2:end-1)*(1:nInterior)*pi);

cosineBasis = scaleColumnsToMatch(cosineBasis,data.PhiFParseval);
sineBasis = scaleColumnsToMatch(sineBasis,data.PhiG);
weightsF = dz;
weightsG = dz(2:end-1).*data.N2Interior/data.g;

fShapeError = relativeFrobeniusError(cosineBasis,data.PhiFParseval);
gShapeError = relativeFrobeniusError(sineBasis,data.PhiG);
fAdjointAgreement = operatorAgreement(cosineBasis,data.PhiFParseval,weightsF,data.GammaF0,"adjoint");
gAdjointAgreement = operatorAgreement(sineBasis,data.PhiG,weightsG,data.GammaG0,"adjoint");
fGalerkinAgreement = operatorAgreement(cosineBasis,data.PhiFParseval,weightsF,data.GammaF0,"galerkin");
gGalerkinAgreement = operatorAgreement(sineBasis,data.PhiG,weightsG,data.GammaG0,"galerkin");

fprintf('\nConstant stratification DCT/DST-I sanity check using geometric weights\n');
fprintf('  F cosine-basis shape error:       %.3e\n', fShapeError);
fprintf('  G sine-basis shape error:         %.3e\n', gShapeError);
fprintf('  F adjoint operator agreement:     %.3e\n', fAdjointAgreement);
fprintf('  G adjoint operator agreement:     %.3e\n', gAdjointAgreement);
fprintf('  F Galerkin operator agreement:    %.3e\n', fGalerkinAgreement);
fprintf('  G Galerkin operator agreement:    %.3e\n', gGalerkinAgreement);
end

function value = operatorAgreement(referencePhi,modePhi,weights,Gamma0,method)
referenceOperator = projectionOperator(referencePhi,weights,Gamma0,method);
modeOperator = projectionOperator(modePhi,weights,Gamma0,method);
value = relativeFrobeniusError(referenceOperator,modeOperator);
end

function operator = projectionOperator(Phi,weights,Gamma0,method)
if method == "adjoint"
    operator = Gamma0 \ (Phi.' .* weights.');
elseif method == "galerkin"
    gram = Phi.'*(weights.*Phi);
    operator = gram \ (Phi.' .* weights.');
else
    error("HydrostaticTransformNumerics:UnknownProjection", "Unknown projection method %s.", method);
end
end

function scaledReference = scaleColumnsToMatch(reference,target)
scaledReference = zeros(size(reference));
for iColumn = 1:size(reference,2)
    scale = reference(:,iColumn)\target(:,iColumn);
    scaledReference(:,iColumn) = scale*reference(:,iColumn);
end
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2)-z(1))/2,(z(3:end)-z(1:end-2))/2,(z(end)-z(end-1))/2);
end

function printComparisonHeader()
fprintf('%-18s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'dz rule','F spec','G spec','max spec','F adj','G adj','F gal','G gal','F inv','G inv','cond F','cond G','sum dz','min dz');
end

function printComparisonRow(row)
fprintf('%-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.dzRule,row.f.specError,row.g.specError,max(row.f.specError,row.g.specError), ...
    row.f.adjointError,row.g.adjointError,row.f.galerkinError,row.g.galerkinError, ...
    row.fInverseError,row.gInverseError,row.f.conditionNumber,row.g.conditionNumber,row.sumDz,row.minDz);
end

function printBandLimitHeader()
fprintf('%8s %10s %10s %10s %10s %10s\n', 'nModes','F spec','G spec','max spec','cond F','cond G');
end

function printBandLimitRow(row)
fprintf('%8d %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.nModes,row.fSpecError,row.gSpecError,row.maxSpecError,row.fConditionNumber,row.gConditionNumber);
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,'fro')/norm(B,'fro');
end

function x = leastSquaresSolution(A,b)
if exist("lsqminnorm","file") == 2
    x = lsqminnorm(A,b);
else
    warningState = warning("off","MATLAB:rankDeficientMatrix");
    cleanup = onCleanup(@() warning(warningState));
    x = A\b;
end
end
