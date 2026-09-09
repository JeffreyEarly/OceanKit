%% Investigate hydrostatic quadrature weights for Parseval accuracy
% This example keeps the hydrostatic mode-adapted nodes fixed and asks which
% positive increment weights dz_i best preserve Parseval for both F and G
% modes when N2(z_i) is sampled physically.
%
% The primary metric is the worst-case relative Parseval error
%
%   E_spec = || Gamma0^{-1/2} (Gamma-Gamma0) Gamma0^{-1/2} ||_2.
%
% Galerkin and direct-inverse round trips are reported as coefficient-recovery
% checks, but they are not quadrature-quality metrics.

baseNPoints = 64;
nPointsList = [24 32 48 64 96 128];
regularizationParameters = [1e-6 1e-4 1e-2 1 1e2];
latitude = 31;
g = 9.81;

cases = [
    hydrostaticCase("exponential",4000,3*2*pi/3600)
    hydrostaticCase("constant",1300,5.2e-3)
    ];

for iCase = 1:length(cases)
    data = modeData(cases(iCase),baseNPoints,latitude,g);
    rows = candidateRows(data,regularizationParameters);

    fprintf('\n%s stratification, nPoints = %d\n', data.name, baseNPoints);
    fprintf('Fixed nodes, physical N2, shared dz weights\n');
    printParsevalHeader();
    for iRow = 1:length(rows)
        printParsevalRow(rows(iRow));
    end

    fprintf('\nCoefficient recovery and conditioning\n');
    printTransformHeader();
    for iRow = 1:length(rows)
        printTransformRow(rows(iRow));
    end

    fprintf('\nWorst normalized Gram directions\n');
    printWorstHeader();
    for iRow = 1:length(rows)
        printWorstRow(rows(iRow));
    end

    if data.name == "constant"
        printConstantStratificationCheck(data,geometricIncrements(data.z));
    end

    fprintf('\nBand-limited full-space weights vs prefix-refit weights\n');
    printBandLimitHeader();
    bestSimple = bestSimpleRow(rows);
    printBandLimitedComparison(data,bestSimple);

    printRecommendationSummary(data,rows);
end

fprintf('\nResolution sweep on fixed nodes\n');
printSweepHeader();
for iCase = 1:length(cases)
    for iPoint = 1:length(nPointsList)
        data = modeData(cases(iCase),nPointsList(iPoint),latitude,g);
        rows = candidateRows(data,regularizationParameters);
        sweepRows = selectedSweepRows(rows);
        for iRow = 1:length(sweepRows)
            printSweepRow(data.nPoints,sweepRows(iRow));
        end
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
    error("HydrostaticQuadratureWeightInvestigation:UnknownCase", "Unknown case %s.", definition.name);
end

[F,G,h] = im.modesAtFrequency(omega);
if definition.name == "constant"
    geostrophicScale = sqrt((definition.N0*definition.N0 - im.f0*im.f0)/(definition.N0*definition.N0));
    F = geostrophicScale*F;
    G = geostrophicScale*G;
end

nInterior = nPoints - 2;
h = h(:);
data.z = z;
data.N2 = im.N2(:);
data.N2Interior = data.N2(2:end-1);
data.PhiFParseval = cat(2,ones(nPoints,1),F(:,1:nInterior));
data.PhiFInverse = cat(2,ones(nPoints,1),F(:,1:(nInterior+1)));
data.PhiG = G(2:end-1,1:nInterior);
data.GammaF0 = diag([definition.depth; h(1:nInterior)]);
data.GammaF0Inverse = diag([definition.depth; h(1:(nInterior+1))]);
data.GammaG0 = eye(nInterior);
data.bF = zeros(size(data.PhiFParseval,2),1);
data.bF(1) = definition.depth;
end

function rows = candidateRows(data,regularizationParameters)
dzGeometric = geometricIncrements(data.z);
rows = diagnosticRow(data,"geometric",dzGeometric,dzGeometric,NaN);
rows(end+1) = diagnosticRow(data,"F-compatible",leastSquaresSolution(data.PhiFParseval.',data.bF),dzGeometric,NaN);

[dzJoint,exitflag] = jointLeastSquaresIncrements(data,maximumModeCount(data));
rows(end+1) = diagnosticRow(data,"joint LS",dzJoint,dzGeometric,exitflag);

[dzPositive,exitflag] = positiveJointLeastSquaresIncrements(data,maximumModeCount(data));
if ~isempty(dzPositive)
    rows(end+1) = diagnosticRow(data,"positive joint LS",dzPositive,dzGeometric,exitflag);
end

for iParameter = 1:length(regularizationParameters)
    lambda = regularizationParameters(iParameter);
    [dzRegularized,exitflag] = regularizedLeastSquaresIncrements(data,maximumModeCount(data),dzGeometric,lambda);
    if ~isempty(dzRegularized)
        rows(end+1) = diagnosticRow(data,sprintf("reg LS %.0e",lambda),dzRegularized,dzGeometric,exitflag);
    end
end

[dzReweighted,exitflag] = reweightedLeastSquaresIncrements(data,1);
if ~isempty(dzReweighted)
    rows(end+1) = diagnosticRow(data,"reweighted LS 1",dzReweighted,dzGeometric,exitflag);
end

[dzReweighted,exitflag] = reweightedLeastSquaresIncrements(data,3);
if ~isempty(dzReweighted)
    rows(end+1) = diagnosticRow(data,"reweighted LS 3",dzReweighted,dzGeometric,exitflag);
end
end

function row = diagnosticRow(data,dzRule,dz,dzReference,exitflag)
row.dzRule = string(dzRule);
row.dz = dz;
row.exitflag = exitflag;
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.negativeDzCount = sum(dz < -1e-10);
row.relativeDzDistance = norm(dz - dzReference)/norm(dzReference);

weightsF = dz;
weightsG = dz(2:end-1).*data.N2Interior/data.g;
row.f = componentDiagnostics(data.PhiFParseval,weightsF,data.GammaF0);
row.g = componentDiagnostics(data.PhiG,weightsG,data.GammaG0);
row.fInverseError = directInverseError(data.PhiFInverse,data.GammaF0Inverse);
row.gInverseError = directInverseError(data.PhiG,data.GammaG0);
row.maxSpecError = max(row.f.specError,row.g.specError);
row.nodeCase = data.name;
row.nPoints = data.nPoints;
end

function diagnostics = componentDiagnostics(Phi,weights,Gamma0)
Gamma = Phi.'*(weights.*Phi);
normalizedDifference = normalizedMatrix(Gamma - Gamma0,Gamma0);
normalizedGram = normalizedMatrix(Gamma,Gamma0);
offdiag = normalizedDifference - diag(diag(normalizedDifference));

diagnostics.specError = norm(symmetrize(normalizedDifference),2);
diagnostics.frobeniusError = norm(normalizedDifference,"fro");
diagnostics.maxDiagError = max(abs(diag(normalizedGram) - 1));
diagnostics.offdiagSpecError = norm(symmetrize(offdiag),2);
diagnostics.galerkinError = normalizedOperatorError(Gamma\Gamma,Gamma0);
diagnostics.conditionNumber = cond(Phi);
diagnostics.worst = worstDirectionDiagnostics(normalizedDifference);
end

function diagnostics = worstDirectionDiagnostics(normalizedDifference)
[V,D] = eig(symmetrize(normalizedDifference));
[value,index] = max(abs(diag(D)));
vector = V(:,index);
weights = abs(vector).^2/sum(abs(vector).^2);
nModes = length(vector);
topCount = max(1,ceil(0.1*nModes));

diagnostics.value = value;
[diagnostics.dominantFraction,diagnostics.dominantMode] = max(weights);
diagnostics.highModeFraction = sum(weights((nModes-topCount+1):nModes));
diagnostics.isHighModeDominated = diagnostics.highModeFraction > 0.5 || diagnostics.dominantMode >= nModes - 1;
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

function matrix = normalizedMatrix(matrix,Gamma0)
scale = diag(1./sqrt(diag(Gamma0)));
matrix = scale*matrix*scale;
end

function matrix = symmetrize(matrix)
matrix = 0.5*(matrix + matrix.');
end

function nModes = maximumModeCount(data)
nModes = min(size(data.PhiFParseval,2)-1,size(data.PhiG,2));
end

function [dz,exitflag] = jointLeastSquaresIncrements(data,nModes)
[A,rhs] = normalizedJointSystem(data,nModes);
[Aeq,beq] = depthConstraint(data);
dz = equalityConstrainedLeastSquares(A,rhs,Aeq,beq);
exitflag = NaN;
end

function [dz,exitflag] = positiveJointLeastSquaresIncrements(data,nModes)
[A,rhs] = normalizedJointSystem(data,nModes);
[Aeq,beq] = depthConstraint(data);
[dz,exitflag] = positiveConstrainedLeastSquares(A,rhs,Aeq,beq);
end

function [dz,exitflag] = regularizedLeastSquaresIncrements(data,nModes,dzReference,lambda)
[A,rhs] = normalizedJointSystem(data,nModes);
relativeScale = 1./max(abs(dzReference),eps);
A = cat(1,A,sqrt(lambda)*diag(relativeScale));
rhs = cat(1,rhs,sqrt(lambda)*(relativeScale.*dzReference));
[Aeq,beq] = depthConstraint(data);
[dz,exitflag] = positiveConstrainedLeastSquares(A,rhs,Aeq,beq);
end

function [dz,exitflag] = reweightedLeastSquaresIncrements(data,nIterations)
nModes = maximumModeCount(data);
[baseA,baseRhs] = normalizedJointSystem(data,nModes);
[Aeq,beq] = depthConstraint(data);
[dz,exitflag] = positiveConstrainedLeastSquares(baseA,baseRhs,Aeq,beq);
if isempty(dz)
    return
end

A = baseA;
rhs = baseRhs;
reweightStrength = sqrt(0.01*size(baseA,1));
for iIteration = 1:nIterations
    [extraA,extraRhs] = worstDirectionRows(data,dz,nModes);
    A = cat(1,A,reweightStrength*extraA);
    rhs = cat(1,rhs,reweightStrength*extraRhs);
    [dz,exitflag] = positiveConstrainedLeastSquares(A,rhs,Aeq,beq);
    if isempty(dz)
        return
    end
end
end

function [A,rhs] = normalizedJointSystem(data,nModes)
nGrid = length(data.z);
PhiF = data.PhiFParseval(:,1:(nModes+1));
GammaF0 = data.GammaF0(1:(nModes+1),1:(nModes+1));
PhiG = data.PhiG(:,1:nModes);
GammaG0 = data.GammaG0(1:nModes,1:nModes);
nF = size(PhiF,2);
nG = size(PhiG,2);

A = zeros(nF*nF + nG*nG,nGrid);
rhs = zeros(size(A,1),1);
row = 0;
gammaF = diag(GammaF0);
gammaG = diag(GammaG0);

for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        A(row,:) = (PhiF(:,iMode).*PhiF(:,jMode)./sqrt(gammaF(iMode)*gammaF(jMode))).';
        rhs(row) = GammaF0(iMode,jMode)/sqrt(gammaF(iMode)*gammaF(jMode));
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(nGrid,1);
        weights(2:end-1) = (data.N2Interior/data.g).*PhiG(:,iMode).*PhiG(:,jMode)./sqrt(gammaG(iMode)*gammaG(jMode));
        A(row,:) = weights.';
        rhs(row) = GammaG0(iMode,jMode)/sqrt(gammaG(iMode)*gammaG(jMode));
    end
end
end

function [A,rhs] = worstDirectionRows(data,dz,nModes)
PhiF = data.PhiFParseval(:,1:(nModes+1));
GammaF0 = data.GammaF0(1:(nModes+1),1:(nModes+1));
weightsF = dz;
gramF = PhiF.'*(weightsF.*PhiF);
vectorF = worstDirection(normalizedMatrix(gramF - GammaF0,GammaF0));
scaleF = diag(1./sqrt(diag(GammaF0)));
modeShapeF = PhiF*scaleF*vectorF;

PhiG = data.PhiG(:,1:nModes);
GammaG0 = data.GammaG0(1:nModes,1:nModes);
weightsG = dz(2:end-1).*data.N2Interior/data.g;
gramG = PhiG.'*(weightsG.*PhiG);
vectorG = worstDirection(normalizedMatrix(gramG - GammaG0,GammaG0));
scaleG = diag(1./sqrt(diag(GammaG0)));
modeShapeG = PhiG*scaleG*vectorG;

A = zeros(2,length(data.z));
A(1,:) = (modeShapeF.^2).';
A(2,2:end-1) = ((data.N2Interior/data.g).*modeShapeG.^2).';
rhs = [1; 1];
end

function vector = worstDirection(normalizedDifference)
[V,D] = eig(symmetrize(normalizedDifference));
[~,index] = max(abs(diag(D)));
vector = V(:,index);
end

function [Aeq,beq] = depthConstraint(data)
Aeq = ones(1,length(data.z));
beq = data.depth;
end

function [x,exitflag] = positiveConstrainedLeastSquares(A,b,Aeq,beq)
x = equalityConstrainedLeastSquares(A,b,Aeq,beq);
if all(x >= -1e-10)
    x(x < 0) = 0;
    exitflag = 0;
    return
end

if exist("lsqlin","file") == 2
    try
        options = optimoptions("lsqlin",Display="off");
        [x,~,~,exitflag] = lsqlin(A,b,[],[],Aeq,beq,zeros(size(A,2),1),[],[],options);
        if ~isempty(x) && exitflag > 0
            return
        end
    catch
    end
end

x = [];
exitflag = NaN;
end

function x = equalityConstrainedLeastSquares(A,b,Aeq,beq)
x0 = Aeq.'*((Aeq*Aeq.')\beq);
Z = null(Aeq);
x = x0 + Z*leastSquaresSolution(A*Z,b - A*x0);
end

function rows = selectedSweepRows(rows)
geometric = find([rows.dzRule] == "geometric",1);
fCompatible = find([rows.dzRule] == "F-compatible",1);
joint = find([rows.dzRule] == "positive joint LS",1);
if isempty(joint)
    joint = find([rows.dzRule] == "joint LS",1);
end
regularized = bestRegularizedIndex(rows);
reweighted = find([rows.dzRule] == "reweighted LS 3",1);
if isempty(reweighted)
    reweighted = find([rows.dzRule] == "reweighted LS 1",1);
end
indices = unique([geometric fCompatible joint regularized reweighted],"stable");
rows = rows(indices);
end

function index = bestRegularizedIndex(rows)
isRegularized = startsWith([rows.dzRule],"reg LS");
indices = find(isRegularized);
if isempty(indices)
    index = [];
    return
end
[~,localIndex] = min([rows(indices).maxSpecError]);
index = indices(localIndex);
end

function row = bestSimpleRow(rows)
simpleRules = ["geometric","F-compatible","positive joint LS","joint LS"];
isSimple = false(size(rows));
for iRule = 1:length(simpleRules)
    isSimple = isSimple | [rows.dzRule] == simpleRules(iRule);
end
indices = find(isSimple);
[~,localIndex] = min([rows(indices).maxSpecError]);
row = rows(indices(localIndex));
end

function row = bestOverallRow(rows)
valid = [rows.negativeDzCount] == 0;
indices = find(valid);
[~,localIndex] = min([rows(indices).maxSpecError]);
row = rows(indices(localIndex));
end

function printBandLimitedComparison(data,fullRow)
modeCounts = bandLimitedModeCounts(data);
for iCount = 1:length(modeCounts)
    nModes = modeCounts(iCount);
    full = prefixDiagnostics(data,fullRow.dz,nModes);
    [dzFit,exitflag] = positiveJointLeastSquaresIncrements(data,nModes);
    if isempty(dzFit)
        [dzFit,exitflag] = jointLeastSquaresIncrements(data,nModes);
    end
    fit = prefixDiagnostics(data,dzFit,nModes);
    fprintf('%8d %-16s %10.3e %10.3e %-16s %10.3e %10.3e %10.3e %10.3e %6.0f\n', ...
        nModes,fullRow.dzRule,full.fSpecError,full.gSpecError,"prefix fit", ...
        fit.fSpecError,fit.gSpecError,min(dzFit),norm(dzFit - fullRow.dz)/norm(fullRow.dz),exitflag);
end
end

function diagnostics = prefixDiagnostics(data,dz,nModes)
PhiF = data.PhiFParseval(:,1:(nModes+1));
GammaF0 = data.GammaF0(1:(nModes+1),1:(nModes+1));
PhiG = data.PhiG(:,1:nModes);
GammaG0 = data.GammaG0(1:nModes,1:nModes);
weightsG = dz(2:end-1).*data.N2Interior/data.g;
fDiagnostics = componentDiagnostics(PhiF,dz,GammaF0);
gDiagnostics = componentDiagnostics(PhiG,weightsG,GammaG0);
diagnostics.fSpecError = fDiagnostics.specError;
diagnostics.gSpecError = gDiagnostics.specError;
end

function modeCounts = bandLimitedModeCounts(data)
maxModes = maximumModeCount(data);
modeCounts = unique([8 16 24 32 48 maxModes]);
modeCounts = modeCounts(modeCounts <= maxModes);
end

function printRecommendationSummary(data,rows)
bestSimple = bestSimpleRow(rows);
bestOverall = bestOverallRow(rows);
improvement = (bestSimple.maxSpecError - bestOverall.maxSpecError)/bestSimple.maxSpecError;
fprintf('\nRecommendation summary for %s stratification\n', data.name);
fprintf('  best simple rule:   %-18s max E_spec %.3e\n', bestSimple.dzRule,bestSimple.maxSpecError);
fprintf('  best overall rule:  %-18s max E_spec %.3e\n', bestOverall.dzRule,bestOverall.maxSpecError);
fprintf('  relative improvement from advanced rules: %.3e\n', improvement);
if data.name == "constant" && bestOverall.maxSpecError > 1e-10
    fprintf('  caution: constant stratification no longer looks DCT/DST-like.\n');
elseif improvement < 0.1
    fprintf('  advanced rules do not materially improve the simple baseline.\n');
elseif bestOverall.relativeDzDistance > 0.25
    fprintf('  advanced rules improve E_spec, but noticeably distort dz.\n');
else
    fprintf('  advanced rules materially improve E_spec without obvious dz distortion.\n');
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
    error("HydrostaticQuadratureWeightInvestigation:UnknownProjection", "Unknown projection method %s.", method);
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

function printParsevalHeader()
fprintf('%-18s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %9s %10s %6s\n', ...
    'dz rule','F spec','G spec','max spec','F fro','G fro','F diag','G diag','F off','G off','min dz','max dz','dz dist','neg');
end

function printParsevalRow(row)
fprintf('%-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %9.3e %10.3e %6d\n', ...
    row.dzRule,row.f.specError,row.g.specError,row.maxSpecError,row.f.frobeniusError,row.g.frobeniusError, ...
    row.f.maxDiagError,row.g.maxDiagError,row.f.offdiagSpecError,row.g.offdiagSpecError, ...
    row.minDz,row.maxDz,row.relativeDzDistance,row.negativeDzCount);
end

function printTransformHeader()
fprintf('%-18s %10s %10s %10s %10s %10s %10s %10s %10s %10s %6s\n', ...
    'dz rule','F gal','G gal','F inv','G inv','cond F','cond G','sum dz','min dz','max dz','exit');
end

function printTransformRow(row)
fprintf('%-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.0f\n', ...
    row.dzRule,row.f.galerkinError,row.g.galerkinError,row.fInverseError,row.gInverseError, ...
    row.f.conditionNumber,row.g.conditionNumber,row.sumDz,row.minDz,row.maxDz,row.exitflag);
end

function printWorstHeader()
fprintf('%-18s %10s %8s %8s %8s %10s %8s %8s %8s\n', ...
    'dz rule','F eig','F mode','F frac','F high','G eig','G mode','G frac','G high');
end

function printWorstRow(row)
fprintf('%-18s %10.3e %8d %8.3f %8.3f %10.3e %8d %8.3f %8.3f\n', ...
    row.dzRule,row.f.worst.value,row.f.worst.dominantMode,row.f.worst.dominantFraction,row.f.worst.highModeFraction, ...
    row.g.worst.value,row.g.worst.dominantMode,row.g.worst.dominantFraction,row.g.worst.highModeFraction);
end

function printBandLimitHeader()
fprintf('%8s %-16s %10s %10s %-16s %10s %10s %10s %10s %6s\n', ...
    'nModes','full dz','full F','full G','refit dz','refit F','refit G','min dz','dz dist','exit');
end

function printSweepHeader()
fprintf('%-13s %8s %-18s %10s %10s %10s %10s %10s %10s\n', ...
    'case','nPoints','dz rule','F spec','G spec','max spec','min dz','dz dist','cond F');
end

function printSweepRow(nPoints,row)
fprintf('%-13s %8d %-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.nodeCase,nPoints,row.dzRule,row.f.specError,row.g.specError,row.maxSpecError, ...
    row.minDz,row.relativeDzDistance,row.f.conditionNumber);
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,"fro")/norm(B,"fro");
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
