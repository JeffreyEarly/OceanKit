%% Optimize physical-sampling quadrature nodes from Chebyshev mode vectors
% This example solves the hydrostatic EVP once, keeps the Chebyshev
% eigenvectors, and evaluates the modes spectrally at trial quadrature
% locations. The optimization keeps physical N2(z_i) and asks whether moving
% z_i improves one shared positive increment vector dz_i for F and G Parseval.
%
% The objective is now the worst-case relative Parseval error
%
%     J = E_{F,spec}^2 + E_{G,spec}^2,
%
% where E_{spec}=||Gamma^{-1/2}(Phi' W Phi-Gamma)Gamma^{-1/2}||_2.
% The top endpoint/Nyquist F mode is kept for square mechanics, but excluded
% from the Parseval diagnostics because its discrete square aliases.

%% Single-resolution node optimization
baseNPoints = 64;
[baseRows, validationError] = compareOptimizedNodes(baseNPoints,true);

fprintf('\nChebyshev evaluator validation error at nPoints = %d: %.3e\n', baseNPoints, validationError);
if validationError > 1e-11
    error("PhysicalSamplingChebyshevNodeOptimization:ValidationFailed", ...
        "The local Chebyshev evaluator does not match modesAtFrequency(0).");
end

fprintf('\nPhysical N2 node optimization at nPoints = %d\n', baseNPoints);
printComparisonHeader();
for iRow = 1:length(baseRows)
    printComparisonRow(baseRows(iRow));
end

fprintf('\nRecommended Parseval diagnostics at nPoints = %d\n', baseNPoints);
printMetricHeader();
for iRow = 1:length(baseRows)
    printMetricRow(baseRows(iRow));
end

fprintf('\nBand-limited worst-case Parseval errors at nPoints = %d\n', baseNPoints);
printBandLimitHeader();
for iRow = 1:length(baseRows)
    if any(baseRows(iRow).dzRule == ["F-compatible","joint positive LS","total energy LS"])
        printBandLimitRows(baseRows(iRow));
    end
end

fprintf('\nTransform-level diagnostics at nPoints = %d\n', baseNPoints);
printTransformHeader();
for iRow = 1:length(baseRows)
    if any(baseRows(iRow).dzRule == ["joint positive LS","total energy LS"])
        printTransformRow(baseRows(iRow));
    end
end

printOptimizationFooter(baseRows);

%% Resolution dependence
nPointsList = [24 32 48 64];

fprintf('\nResolution dependence with fixed and optimized G^{N+1} nodes\n');
printSweepHeader();
for iPoint = 1:length(nPointsList)
    rows = compareOptimizedNodes(nPointsList(iPoint),false);
    for iRow = 1:length(rows)
        if any(rows(iRow).dzRule == ["joint positive LS","total energy LS"])
            printSweepRow(nPointsList(iPoint),rows(iRow));
        end
    end
end

%% Local helpers
function [rows, validationError] = compareOptimizedNodes(nPoints,shouldValidate)
reference = chebyshevReference(nPoints);
zInitial = reference.zInitial;
validationError = NaN;
if shouldValidate
    validationError = validateChebyshevEvaluator(reference,zInitial);
end

fixedData = modeData(reference,zInitial,"G^{N+1}");
rows = rowDiagnostics(fixedData,"geometric",geometricIncrements(zInitial),NaN,0);
rows(end+1) = rowDiagnostics(fixedData,"F-compatible",leastSquaresSolution(fixedData.PhiF.',fixedData.bF),NaN,0);
[dzJoint,dzRule,exitflag] = bestJointIncrements(fixedData);
rows(end+1) = rowDiagnostics(fixedData,dzRule,dzJoint,exitflag,0);
[dzEnergy,dzRule,exitflag] = bestTotalEnergyIncrements(fixedData);
rows(end+1) = rowDiagnostics(fixedData,dzRule,dzEnergy,exitflag,0);

[zFull,fullExitflag] = optimizeFullNodes(reference,zInitial);
fullData = modeData(reference,zFull,"optimized full");
[dzFull,dzRule,incrementExitflag] = bestJointIncrements(fullData);
rows(end+1) = rowDiagnostics(fullData,dzRule,dzFull,fullExitflag,max(abs(zFull - zInitial)));
rows(end).incrementExitflag = incrementExitflag;
[dzFullEnergy,dzRule,incrementExitflag] = bestTotalEnergyIncrements(fullData);
rows(end+1) = rowDiagnostics(fullData,dzRule,dzFullEnergy,fullExitflag,max(abs(zFull - zInitial)));
rows(end).incrementExitflag = incrementExitflag;

[zSmooth,smoothExitflag] = optimizeSmoothNodes(reference,zInitial);
smoothData = modeData(reference,zSmooth,"optimized smooth");
[dzSmooth,dzRule,incrementExitflag] = bestJointIncrements(smoothData);
rows(end+1) = rowDiagnostics(smoothData,dzRule,dzSmooth,smoothExitflag,max(abs(zSmooth - zInitial)));
rows(end).incrementExitflag = incrementExitflag;
[dzSmoothEnergy,dzRule,incrementExitflag] = bestTotalEnergyIncrements(smoothData);
rows(end+1) = rowDiagnostics(smoothData,dzRule,dzSmoothEnergy,smoothExitflag,max(abs(zSmooth - zInitial)));
rows(end).incrementExitflag = incrementExitflag;
end

function reference = chebyshevReference(nPoints)
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

im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,nModes=nModes,g=g);
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;
im.gridFrequency = omega;

[A,B] = im.EigenmatricesForFrequency(omega);
[VRaw,hRaw] = solveGEPForExample(im,A,B);
zInitial = quadraturePointsForModesForExample(im,nPoints,VRaw,hRaw);
[V,h] = normalizedChebyshevVectors(im,VRaw,hRaw,nModes);

reference.im = im;
reference.N2 = N2;
reference.Lz = Lz;
reference.g = g;
reference.nPoints = nPoints;
reference.nModes = nModes;
reference.VRaw = VRaw;
reference.hRaw = hRaw;
reference.V = V;
reference.dV = differentiatedChebyshevVectors(im,V);
reference.h = h;
reference.zInitial = zInitial;
end

function [VCheb,h] = solveGEPForExample(im,A,B)
if any(isnan(A(:))) || any(isnan(B(:)))
    error("PhysicalSamplingChebyshevNodeOptimization:NaNInMatrix", ...
        "EVP setup failed. Found at least one NaN in the generalized EVP matrices.");
end

[V,D] = eig(A,B);
[h,permutation] = sort(real(im.hFromLambda(diag(D))),'descend');
VCheb = V(:,permutation);
end

function z = quadraturePointsForModesForExample(im,nPoints,VCheb,h)
resolvedModes = ceil(find(h > 0,1,"last")/2);
if isempty(resolvedModes) || resolvedModes < nPoints
    resolvedModeCount = 0;
    if ~isempty(resolvedModes)
        resolvedModeCount = resolvedModes;
    end
    error("PhysicalSamplingChebyshevNodeOptimization:NeedMorePoints", ...
        "Returned %d valid modes (%d quadrature points requested) using nEVP=%d.", ...
        resolvedModeCount,nPoints,im.nEVP);
end

rootMode = nPoints - 1;
rootsVar = InternalModesSpectral.FindRootsFromChebyshevVector(VCheb(:,rootMode),im.xDomain);
rootsVar(rootsVar < im.xMin) = im.xMin;
rootsVar(rootsVar > im.xMax) = im.xMax;
rootsVar = cat(1,im.xMin,rootsVar,im.xMax);
rootsVar = unique(rootsVar,"stable");

while length(rootsVar) > nPoints
    rootsVar = sort(rootsVar);
    F = InternalModesSpectral.IntegrateChebyshevVector(VCheb(:,rootMode));
    value = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(rootsVar,im.xDomain,F);
    dv = diff(value);
    [~,minIndex] = min(abs(dv));
    rootsVar(minIndex+1) = [];
end

if length(rootsVar) < nPoints
    error("PhysicalSamplingChebyshevNodeOptimization:NeedMorePoints", ...
        "Returned %d unique roots (requested %d). Maybe need more EVP points.", length(rootsVar),nPoints);
end

z = InternalModesSpectral.fInverseBisection(im.x_function,reshape(rootsVar,[],1),min(im.zDomain),max(im.zDomain),1e-12);
end

function [V,h] = normalizedChebyshevVectors(im,VRaw,hRaw,nModes)
V = VRaw(:,1:nModes);
h = reshape(hRaw(1:nModes),1,[]);
maxIndexZ = find(im.N2_xLobatto - im.gridFrequency*im.gridFrequency > 0, 1, "first");
if maxIndexZ > 1
    maxIndexZ = maxIndexZ - 1;
elseif isempty(maxIndexZ)
    maxIndexZ = 1;
end

for iMode = 1:nModes
    Fj = im.FFromVCheb(V(:,iMode),h(iMode));
    Gj = im.GFromVCheb(V(:,iMode),h(iMode));
    switch im.normalization
        case Normalization.uMax
            A = max(abs(Fj));
        case Normalization.wMax
            A = max(abs(Gj));
        case Normalization.kConstant
            A = sqrt(im.GNorm(Gj));
        case Normalization.geostrophic
            A = sqrt(im.GeostrophicNorm(Gj));
    end
    if Fj(maxIndexZ) < 0
        A = -A;
    end
    V(:,iMode) = V(:,iMode)/A;
end
end

function dV = differentiatedChebyshevVectors(im,V)
dV = zeros(size(V));
for iMode = 1:size(V,2)
    dV(:,iMode) = im.Diff1_xCheb(V(:,iMode));
end
end

function validationError = validateChebyshevEvaluator(reference,z)
reference.im.z = z;
[FReference,GReference] = reference.im.modesAtFrequency(0);
[~,~,F,G] = evaluatedModeMatrices(reference,z);
validationError = max(relativeFrobeniusError(F,FReference),relativeFrobeniusError(G,GReference));
end

function data = modeData(reference,z,nodeRule)
[PhiF,PhiG] = evaluatedModeMatrices(reference,z);
data.nodeRule = string(nodeRule);
data.Lz = reference.Lz;
data.g = reference.g;
data.z = z;
data.N2 = reference.N2(z);
data.N2Interior = reference.N2(z(2:end-1));
data.PhiF = PhiF;
data.PhiFParseval = PhiF(:,1:end-1);
data.PhiG = PhiG;
data.GammaF = diag([reference.Lz; reference.h(:)]);
data.GammaFParseval = data.GammaF(1:end-1,1:end-1);
data.GammaG = eye(size(PhiG,2));
data.bF = zeros(size(PhiF,2),1);
data.bF(1) = reference.Lz;
data.bFParseval = data.bF(1:end-1);
end

function [PhiF,PhiG,F,G] = evaluatedModeMatrices(reference,z)
x = reference.im.x_function(z);
x(x > max(reference.im.xLobatto)) = max(reference.im.xLobatto);
x(x < min(reference.im.xLobatto)) = min(reference.im.xLobatto);
T = InternalModesSpectral.ChebyshevTransformForGrid(reference.im.xLobatto,x);

Nz = length(z);
F = zeros(Nz,reference.nModes);
G = zeros(Nz,reference.nModes);
N = sqrt(reference.N2(z));
for iMode = 1:reference.nModes
    G(:,iMode) = T(reference.V(:,iMode));
    F(:,iMode) = reference.h(iMode)*N.*T(reference.dV(:,iMode));
end

PhiF = cat(2,ones(Nz,1),F);
PhiG = G(2:end-1,1:end-1);
end

function [zOptimized,exitflag] = optimizeFullNodes(reference,zInitial)
if exist("fmincon","file") ~= 2
    zOptimized = zInitial;
    exitflag = NaN;
    return;
end

yInitial = normalizedXInterior(reference,zInitial);
starts = fullNodeStarts(yInitial);
[yOptimized,exitflag] = bestFullNodeRun(reference,yInitial,starts);
zOptimized = zFromNormalizedX(reference,yOptimized);
end

function [zOptimized,exitflag] = optimizeSmoothNodes(reference,zInitial)
if exist("fmincon","file") ~= 2
    zOptimized = zInitial;
    exitflag = NaN;
    return;
end

yInitial = normalizedXInterior(reference,zInitial);
basis = smoothPerturbationBasis(length(yInitial),min(6,max(2,floor(length(yInitial)/4))));
[coefficients,exitflag] = bestSmoothRun(reference,yInitial,basis);
yOptimized = yInitial + basis*coefficients;
zOptimized = zFromNormalizedX(reference,yOptimized);
end

function [bestY,bestExitflag] = bestFullNodeRun(reference,yInitial,starts)
nInterior = length(yInitial);
minSpacing = 1e-6;
Aineq = zeros(nInterior - 1,nInterior);
for iRow = 1:(nInterior - 1)
    Aineq(iRow,iRow) = 1;
    Aineq(iRow,iRow+1) = -1;
end
bineq = -minSpacing*ones(nInterior - 1,1);
lb = minSpacing*ones(nInterior,1);
ub = (1 - minSpacing)*ones(nInterior,1);
options = optimoptions("fmincon", ...
    Algorithm="sqp", ...
    Display="off", ...
    MaxIterations=12, ...
    MaxFunctionEvaluations=max(140,3*nInterior), ...
    OptimalityTolerance=1e-8, ...
    StepTolerance=1e-8);

bestY = yInitial;
bestObjective = objectiveForNormalizedX(reference,yInitial);
bestExitflag = NaN;
for iStart = 1:size(starts,2)
    objective = @(yInterior) objectiveForNormalizedX(reference,yInterior);
    [candidateY,candidateObjective,exitflag] = fmincon(objective,starts(:,iStart),Aineq,bineq,[],[],lb,ub,[],options);
    if isfinite(candidateObjective) && candidateObjective < bestObjective
        bestY = candidateY;
        bestObjective = candidateObjective;
        bestExitflag = exitflag;
    end
end
end

function [bestCoefficients,bestExitflag] = bestSmoothRun(reference,yInitial,basis)
nCoefficients = size(basis,2);
amplitude = 0.15*min(diff(cat(1,0,yInitial,1)));
starts = zeros(nCoefficients,1);
for iCoefficient = 1:nCoefficients
    start = zeros(nCoefficients,1);
    start(iCoefficient) = amplitude*(-1)^iCoefficient;
    starts(:,end+1) = start;
end

options = optimoptions("fmincon", ...
    Algorithm="sqp", ...
    Display="off", ...
    MaxIterations=24, ...
    MaxFunctionEvaluations=220, ...
    OptimalityTolerance=1e-8, ...
    StepTolerance=1e-8);
lb = -4*amplitude*ones(nCoefficients,1);
ub = 4*amplitude*ones(nCoefficients,1);
constraint = @(coefficients) smoothCoordinateConstraint(yInitial,basis,coefficients);

bestCoefficients = zeros(nCoefficients,1);
bestObjective = objectiveForSmoothCoefficients(reference,yInitial,basis,bestCoefficients);
bestExitflag = NaN;
for iStart = 1:size(starts,2)
    objective = @(coefficients) objectiveForSmoothCoefficients(reference,yInitial,basis,coefficients);
    [candidateCoefficients,candidateObjective,exitflag] = fmincon( ...
        objective,starts(:,iStart),[],[],[],[],lb,ub,constraint,options);
    if isfinite(candidateObjective) && candidateObjective < bestObjective
        bestCoefficients = candidateCoefficients;
        bestObjective = candidateObjective;
        bestExitflag = exitflag;
    end
end
end

function starts = fullNodeStarts(yInitial)
spacing = min(diff(cat(1,0,yInitial,1)));
amplitude = 0.10*spacing;
t = linspace(0,1,length(yInitial)).';
starts = yInitial;
for iMode = 1:3
    starts(:,end+1) = yInitial + amplitude*sin(iMode*pi*t);
end
end

function basis = smoothPerturbationBasis(nInterior,nBasis)
t = linspace(0,1,nInterior).';
basis = zeros(nInterior,nBasis);
for iBasis = 1:nBasis
    basis(:,iBasis) = sin(iBasis*pi*t);
end
end

function [cineq,ceq] = smoothCoordinateConstraint(yInitial,basis,coefficients)
minSpacing = 1e-6;
y = yInitial + basis*coefficients(:);
allY = cat(1,0,y,1);
cineq = cat(1,minSpacing - diff(allY),minSpacing - y,y - (1 - minSpacing));
ceq = [];
end

function objective = objectiveForSmoothCoefficients(reference,yInitial,basis,coefficients)
y = yInitial + basis*coefficients(:);
objective = objectiveForNormalizedX(reference,y);
end

function objective = objectiveForNormalizedX(reference,yInterior)
if any(diff(cat(1,0,yInterior(:),1)) <= 0)
    objective = realmax;
    return;
end

z = zFromNormalizedX(reference,yInterior);
data = modeData(reference,z,"optimized");
[dz,~,~] = bestJointIncrements(data);
objective = gramObjective(data,dz);
if ~isfinite(objective)
    objective = realmax;
end
end

function y = normalizedXInterior(reference,z)
x = reference.im.x_function(z);
y = (x(2:end-1) - reference.im.xMin)/reference.im.Lx;
end

function z = zFromNormalizedX(reference,yInterior)
x = reference.im.xMin + reference.im.Lx*cat(1,0,yInterior(:),1);
z = InternalModesSpectral.fInverseBisection( ...
    reference.im.x_function,x,min(reference.im.zDomain),max(reference.im.zDomain),1e-12);
end

function [dz,dzRule,exitflag] = bestJointIncrements(data)
[A,rhs] = jointLeastSquaresSystem(data);
Aeq = ones(1,length(data.z));
beq = data.Lz;
if exist("lsqlin","file") == 2
    try
        options = optimoptions("lsqlin",Display="off");
        [dz,~,~,exitflag] = lsqlin(A,rhs,[],[],Aeq,beq,zeros(length(data.z),1),[],[],options);
        if ~isempty(dz) && exitflag > 0
            dzRule = "joint positive LS";
            return;
        end
    catch
    end
end

dz = equalityConstrainedLeastSquares(A,rhs,Aeq,beq);
dzRule = "joint LS fallback";
exitflag = NaN;
end

function [dz,dzRule,exitflag] = bestTotalEnergyIncrements(data)
[A,rhs] = totalEnergyLeastSquaresSystem(data);
Aeq = ones(1,length(data.z));
beq = data.Lz;
if exist("lsqlin","file") == 2
    try
        options = optimoptions("lsqlin",Display="off");
        [dz,~,~,exitflag] = lsqlin(A,rhs,[],[],Aeq,beq,zeros(length(data.z),1),[],[],options);
        if ~isempty(dz) && exitflag > 0
            dzRule = "total energy LS";
            return;
        end
    catch
    end
end

dz = equalityConstrainedLeastSquares(A,rhs,Aeq,beq);
dzRule = "energy LS fallback";
exitflag = NaN;
end

function [A,rhs] = jointLeastSquaresSystem(data)
Nz = length(data.z);
nF = size(data.PhiFParseval,2);
nG = size(data.PhiG,2);

nRows = nF*nF + nG*nG;
A = zeros(nRows,Nz);
rhs = zeros(nRows,1);
row = 0;
gammaF = diag(data.GammaFParseval);

for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        scale = sqrt(gammaF(iMode)*gammaF(jMode));
        A(row,:) = (data.PhiFParseval(:,iMode).*data.PhiFParseval(:,jMode)).'/scale;
        rhs(row) = data.GammaFParseval(iMode,jMode)/scale;
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(Nz,1);
        weights(2:end-1) = (data.N2Interior/data.g).*data.PhiG(:,iMode).*data.PhiG(:,jMode);
        A(row,:) = weights.';
        rhs(row) = data.GammaG(iMode,jMode);
    end
end
end

function [A,rhs] = totalEnergyLeastSquaresSystem(data)
Nz = length(data.z);
nModes = min(size(data.PhiFParseval,2) - 1,size(data.PhiG,2));
gammaF = diag(data.GammaFParseval);

A = zeros(nModes*nModes,Nz);
rhs = zeros(size(A,1),1);
row = 0;
for iMode = 1:nModes
    iF = iMode + 1;
    for jMode = 1:nModes
        jF = jMode + 1;
        row = row + 1;
        fWeights = data.PhiFParseval(:,iF).*data.PhiFParseval(:,jF)/sqrt(gammaF(iF)*gammaF(jF));
        gWeights = zeros(Nz,1);
        gWeights(2:end-1) = (data.N2Interior/data.g).*data.PhiG(:,iMode).*data.PhiG(:,jMode);
        A(row,:) = (fWeights + gWeights).';
        rhs(row) = 2*(iMode == jMode);
    end
end
end

function row = rowDiagnostics(data,dzRule,dz,exitflag,maxDeltaZ)
gramF = data.PhiFParseval.'*(dz.*data.PhiFParseval);
weightsG = dz(2:end-1).*data.N2Interior/data.g;
gramG = data.PhiG.'*(weightsG.*data.PhiG);
[fGramError,fDiagError,fOffError,fSpecError] = gramErrors(gramF,data.GammaFParseval);
[gGramError,gDiagError,gOffError,gSpecError] = gramErrors(gramG,data.GammaG);
[~,~,~,~,fMaxDiagError,fOffSpecError] = gramErrors(gramF,data.GammaFParseval);
[~,~,~,~,gMaxDiagError,gOffSpecError] = gramErrors(gramG,data.GammaG);
totalEnergyError = totalEnergyGramError(data,dz);

row.nodeRule = data.nodeRule;
row.dzRule = string(dzRule);
row.fGramError = fGramError;
row.fDiagError = fDiagError;
row.fOffError = fOffError;
row.fSpecError = fSpecError;
row.fMaxDiagError = fMaxDiagError;
row.fOffSpecError = fOffSpecError;
row.gGramError = gGramError;
row.gDiagError = gDiagError;
row.gOffError = gOffError;
row.gSpecError = gSpecError;
row.gMaxDiagError = gMaxDiagError;
row.gOffSpecError = gOffSpecError;
row.totalEnergyError = totalEnergyError;
row.fIntegralError = norm(data.PhiFParseval.'*dz - data.bFParseval)/norm(data.bFParseval);
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.negativeDzCount = sum(dz < 0);
row.objective = row.fSpecError^2 + row.gSpecError^2;
row.maxDeltaZ = maxDeltaZ;
row.exitflag = exitflag;
row.incrementExitflag = exitflag;
row.transform = transformDiagnostics(data,dz);
row.bandLimits = bandLimitDiagnostics(data,dz);
end

function objective = gramObjective(data,dz)
gramF = data.PhiFParseval.'*(dz.*data.PhiFParseval);
weightsG = dz(2:end-1).*data.N2Interior/data.g;
gramG = data.PhiG.'*(weightsG.*data.PhiG);
fGramError = normalizedSpectralGramError(gramF,data.GammaFParseval);
gGramError = normalizedSpectralGramError(gramG,data.GammaG);
objective = fGramError^2 + gGramError^2;
end

function diagnostics = transformDiagnostics(data,dz)
weightsG = dz(2:end-1).*data.N2Interior/data.g;
AFCan = data.GammaFParseval \ (data.PhiFParseval.'.*dz.');
AGCan = data.PhiG.'.*weightsG.';
AFLS = weightedPseudoinverse(data.PhiFParseval,dz);
AGLS = weightedPseudoinverse(data.PhiG,weightsG);

identityF = eye(size(data.PhiFParseval,2));
identityG = eye(size(data.PhiG,2));
diagnostics.fCanRoundTrip = relativeFrobeniusError(AFCan*data.PhiFParseval,identityF);
diagnostics.gCanRoundTrip = relativeFrobeniusError(AGCan*data.PhiG,identityG);
diagnostics.fLSRoundTrip = relativeFrobeniusError(AFLS*data.PhiFParseval,identityF);
diagnostics.gLSRoundTrip = relativeFrobeniusError(AGLS*data.PhiG,identityG);
diagnostics.fParsevalError = parsevalErrorF(data,dz);
diagnostics.gParsevalError = parsevalErrorG(data,weightsG);
diagnostics.fConditionNumber = cond(data.PhiFParseval);
diagnostics.gConditionNumber = cond(data.PhiG);
end

function error = parsevalErrorF(data,dz)
nF = size(data.PhiFParseval,2);
a = deterministicCoefficients(nF);
u = data.PhiFParseval*a;
left = u.'*(dz.*u);
right = a.'*data.GammaFParseval*a;
error = abs(left - right)/abs(right);
end

function error = parsevalErrorG(data,weightsG)
nG = size(data.PhiG,2);
a = deterministicCoefficients(nG);
eta = data.PhiG*a;
left = eta.'*(weightsG.*eta);
right = a.'*a;
error = abs(left - right)/abs(right);
end

function a = deterministicCoefficients(n)
j = (1:n).';
a = cos(j) + sin(sqrt(2)*j);
end

function A = weightedPseudoinverse(Phi,weights)
weightedPhi = weights.*Phi;
A = (Phi.'*weightedPhi)\weightedPhi.';
end

function [totalError,diagError,offError,specError,maxDiagError,offSpecError] = gramErrors(gram,target)
totalError = relativeFrobeniusError(gram,target);
diagError = norm(diag(gram) - diag(target))/norm(diag(target));
offDifference = (gram - target) - diag(diag(gram - target));
offError = norm(offDifference,"fro")/norm(target,"fro");
specError = normalizedSpectralGramError(gram,target);
maxDiagError = max(abs(diag(gram)./diag(target) - 1));
offSpecError = normalizedSpectralOffdiagError(gram,target);
end

function value = normalizedSpectralGramError(gram,target)
scale = diag(1./sqrt(diag(target)));
value = norm(scale*(gram - target)*scale,2);
end

function value = normalizedSpectralOffdiagError(gram,target)
difference = gram - target;
offDifference = difference - diag(diag(difference));
scale = diag(1./sqrt(diag(target)));
value = norm(scale*offDifference*scale,2);
end

function diagnostics = bandLimitDiagnostics(data,dz)
maxModes = min(size(data.PhiFParseval,2),size(data.PhiG,2));
modeCounts = unique([8 16 24 32 48 maxModes]);
modeCounts = modeCounts(modeCounts <= maxModes);
weightsG = dz(2:end-1).*data.N2Interior/data.g;

diagnostics.modeCounts = modeCounts;
diagnostics.fSpecError = zeros(size(modeCounts));
diagnostics.gSpecError = zeros(size(modeCounts));
diagnostics.totalEnergyError = zeros(size(modeCounts));
diagnostics.objective = zeros(size(modeCounts));
for iCount = 1:length(modeCounts)
    nModes = modeCounts(iCount);
    phiF = data.PhiFParseval(:,1:nModes);
    gammaF = data.GammaFParseval(1:nModes,1:nModes);
    gramF = phiF.'*(dz.*phiF);

    phiG = data.PhiG(:,1:nModes);
    gammaG = data.GammaG(1:nModes,1:nModes);
    gramG = phiG.'*(weightsG.*phiG);

    diagnostics.fSpecError(iCount) = normalizedSpectralGramError(gramF,gammaF);
    diagnostics.gSpecError(iCount) = normalizedSpectralGramError(gramG,gammaG);
    diagnostics.totalEnergyError(iCount) = totalEnergyGramError(data,dz,nModes);
    diagnostics.objective(iCount) = diagnostics.fSpecError(iCount)^2 + diagnostics.gSpecError(iCount)^2;
end
end

function value = totalEnergyGramError(data,dz,nModes)
if nargin < 3
    nModes = min(size(data.PhiFParseval,2) - 1,size(data.PhiG,2));
end

phiF = data.PhiFParseval(:,2:(nModes + 1));
gammaF = data.GammaFParseval(2:(nModes + 1),2:(nModes + 1));
gramF = phiF.'*(dz.*phiF);

weightsG = dz(2:end-1).*data.N2Interior/data.g;
phiG = data.PhiG(:,1:nModes);
gramG = phiG.'*(weightsG.*phiG);

combinedEnergy = normalizeGram(gramF,gammaF) + gramG;
value = 0.5*norm(combinedEnergy - 2*eye(nModes),2);
end

function normalizedGram = normalizeGram(gram,target)
scale = diag(1./sqrt(diag(target)));
normalizedGram = scale*gram*scale;
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2)-z(1))/2,(z(3:end)-z(1:end-2))/2,(z(end)-z(end-1))/2);
end

function printComparisonHeader()
fprintf('%-16s %-18s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %6s\n', ...
    'node rule','dz rule','F fro','F spec','F off','G fro','G spec','G off', ...
    'F int','min dz','J spec','max |dz|','exit');
end

function printComparisonRow(row)
fprintf('%-16s %-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %6.0f\n', ...
    row.nodeRule,row.dzRule,row.fGramError,row.fSpecError,row.fOffError,row.gGramError, ...
    row.gSpecError,row.gOffError,row.fIntegralError,row.minDz,row.objective,row.maxDeltaZ,row.exitflag);
end

function printTransformHeader()
fprintf('%-16s %-16s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'node rule','dz rule','F can','G can','F LS','G LS','F par','G par','cond F','cond G');
end

function printMetricHeader()
fprintf('%-16s %-18s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'node rule','dz rule','F spec','F maxdiag','F offsp','G spec','G maxdiag','G offsp','E total','cond F','cond G');
end

function printMetricRow(row)
transform = row.transform;
fprintf('%-16s %-18s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.nodeRule,row.dzRule,row.fSpecError,row.fMaxDiagError,row.fOffSpecError, ...
    row.gSpecError,row.gMaxDiagError,row.gOffSpecError,row.totalEnergyError, ...
    transform.fConditionNumber,transform.gConditionNumber);
end

function printBandLimitHeader()
fprintf('%-16s %-18s %8s %10s %10s %10s %10s\n', ...
    'node rule','dz rule','nModes','F spec','G spec','E total','J spec');
end

function printBandLimitRows(row)
for iCount = 1:length(row.bandLimits.modeCounts)
    fprintf('%-16s %-18s %8d %10.3e %10.3e %10.3e %10.3e\n', ...
        row.nodeRule,row.dzRule,row.bandLimits.modeCounts(iCount), ...
        row.bandLimits.fSpecError(iCount),row.bandLimits.gSpecError(iCount), ...
        row.bandLimits.totalEnergyError(iCount),row.bandLimits.objective(iCount));
end
end

function printTransformRow(row)
transform = row.transform;
fprintf('%-16s %-16s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.nodeRule,row.dzRule,transform.fCanRoundTrip,transform.gCanRoundTrip,transform.fLSRoundTrip, ...
    transform.gLSRoundTrip,transform.fParsevalError,transform.gParsevalError, ...
    transform.fConditionNumber,transform.gConditionNumber);
end

function printSweepHeader()
fprintf('%8s %-16s %-16s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'nPoints','node rule','dz rule','F spec','G spec','E total','F off','G off','F par','G par','J spec','max |dz|');
end

function printSweepRow(nPoints,row)
transform = row.transform;
fprintf('%8d %-16s %-16s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    nPoints,row.nodeRule,row.dzRule,row.fSpecError,row.gSpecError,row.totalEnergyError,row.fOffError,row.gOffError, ...
    transform.fParsevalError,transform.gParsevalError,row.objective,row.maxDeltaZ);
end

function printOptimizationFooter(rows)
fixedIndex = find([rows.nodeRule] == "G^{N+1}" & [rows.dzRule] == "joint positive LS",1);
candidateIndex = find([rows.dzRule] == "joint positive LS" & [rows.nodeRule] ~= "G^{N+1}");
if isempty(fixedIndex) || isempty(candidateIndex)
    return;
end

[bestObjective,bestLocalIndex] = min([rows(candidateIndex).objective]);
bestIndex = candidateIndex(bestLocalIndex);
fixedObjective = rows(fixedIndex).objective;
improvement = (fixedObjective - bestObjective)/fixedObjective;
fprintf('\nBest optimized row: %s, J improvement relative to fixed G^{N+1}: %.3e\n', ...
    rows(bestIndex).nodeRule,improvement);
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

function x = equalityConstrainedLeastSquares(A,b,Aeq,beq)
normalMatrix = A.'*A;
rhs = A.'*b;
kkt = [normalMatrix Aeq.'; Aeq zeros(size(Aeq,1))];
solution = kkt\[rhs; beq];
x = solution(1:size(A,2));
end
