%% Generalized Gaussian quadrature for retained hydrostatic bands
% This example asks whether a retained modal band benefits from generalized
% Gaussian quadrature: moving the quadrature nodes z_i and fitting one shared
% positive increment vector dz_i. The score is the worst-case relative
% Parseval error on the retained hydrostatic F/G band.
%
% The EVP is solved once per resolution for variable stratification. Trial
% nodes are evaluated from the Chebyshev vectors, so the outer optimization
% explores quadrature nodes without repeatedly solving the Sturm-Liouville
% problem.

baseNPoints = 64;
nPointsList = [24 32 48 64 96 128];
retainedFractions = [1/2 2/3];
latitude = 31;
g = 9.81;
rng(42,"twister");

definitions(1) = quadratureCase("exponential",4000,3*2*pi/3600);
definitions(2) = quadratureCase("constant",1300,5.2e-3);

for iCase = 1:length(definitions)
    reference = quadratureReference(definitions(iCase),baseNPoints,latitude,g);
    fprintf('\n%s stratification, nPoints = %d\n', reference.name, baseNPoints);
    if reference.name == "exponential"
        fprintf('Chebyshev evaluator validation error: %.3e\n', validateEvaluator(reference,reference.zInitial));
    end

    for iFraction = 1:length(retainedFractions)
        retainedFraction = retainedFractions(iFraction);
        r = retainedModeCount(reference,retainedFraction);
        fprintf('\nRetained band: r = %d modes, fraction = %.3f\n', r, retainedFraction);
        printComparisonHeader();
        rows = retainedBandRows(reference,r,retainedFraction,optimizationBudget("base"));
        for iRow = 1:length(rows)
            printComparisonRow(rows(iRow));
        end
        printRecommendation(reference,rows);
    end

    if reference.name == "constant"
        printConstantStratificationCheck(reference,retainedModeCount(reference,2/3));
    end

    fprintf('\nBand-growth comparison at nPoints = %d\n', baseNPoints);
    printBandGrowthHeader();
    printBandGrowthRows(reference,optimizationBudget("band"));
end

fprintf('\nResolution sweep\n');
printSweepHeader();
for iCase = 1:length(definitions)
    for iPoint = 1:length(nPointsList)
        reference = quadratureReference(definitions(iCase),nPointsList(iPoint),latitude,g);
        for iFraction = 1:length(retainedFractions)
            retainedFraction = retainedFractions(iFraction);
            r = retainedModeCount(reference,retainedFraction);
            rows = retainedBandRows(reference,r,retainedFraction,optimizationBudget("sweep"));
            fixedIndex = find([rows.nodeRule] == "G^{N+1}" & [rows.dzRule] ~= "geometric",1);
            smoothIndex = find([rows.nodeRule] == "optimized smooth",1);
            printSweepRow(reference,retainedFraction,rows(fixedIndex),rows(smoothIndex));
        end
    end
end

%% Local helpers
function definition = quadratureCase(name,depth,N0)
definition.name = string(name);
definition.depth = depth;
definition.N0 = N0;
end

function budget = optimizationBudget(kind)
budget.kind = string(kind);
switch budget.kind
    case "base"
        budget.maxIterations = 14;
        budget.maxFunctionEvaluations = 140;
        budget.nSmoothModes = 6;
        budget.nRandomStarts = 3;
        budget.shouldRunFullNodes = false;
    case "band"
        budget.maxIterations = 8;
        budget.maxFunctionEvaluations = 90;
        budget.nSmoothModes = 5;
        budget.nRandomStarts = 1;
        budget.shouldRunFullNodes = false;
    case "sweep"
        budget.maxIterations = 6;
        budget.maxFunctionEvaluations = 70;
        budget.nSmoothModes = 4;
        budget.nRandomStarts = 0;
        budget.shouldRunFullNodes = false;
    otherwise
        error("HydrostaticRetainedBandGeneralizedQuadrature:UnknownBudget", ...
            "Unknown optimization budget %s.", kind);
end
end

function reference = quadratureReference(definition,nPoints,latitude,g)
omega = 0;
reference.name = definition.name;
reference.depth = definition.depth;
reference.N0 = definition.N0;
reference.g = g;
reference.nPoints = nPoints;
reference.nInterior = nPoints - 2;
reference.nModes = nPoints - 1;
reference.latitude = latitude;

if definition.name == "exponential"
    L_gm = 1300;
    N2 = @(z) definition.N0*definition.N0*exp(2*z/L_gm);
    nEVP = max(256,ceil(2.1*nPoints));
    zDomain = [-definition.depth 0];
    zReference = linspace(zDomain(1),zDomain(2),1024).';
    im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,nModes=reference.nModes,g=g);
    im.normalization = Normalization.geostrophic;
    im.upperBoundary = UpperBoundary.rigidLid;
    im.gridFrequency = omega;

    [A,B] = im.EigenmatricesForFrequency(omega);
    [VRaw,hRaw] = solveGEPForExample(im,A,B);
    zInitial = quadraturePointsForModesForExample(im,nPoints,VRaw,hRaw);
    [V,h] = normalizedChebyshevVectors(im,VRaw,hRaw,reference.nModes);

    reference.kind = "chebyshev";
    reference.im = im;
    reference.N2 = N2;
    reference.h = h(:);
    reference.V = V;
    reference.dV = differentiatedChebyshevVectors(im,V);
    reference.zInitial = zInitial;
elseif definition.name == "constant"
    zInitial = linspace(-definition.depth,0,nPoints).';
    im = InternalModesConstantStratification(N0=definition.N0,zIn=[-definition.depth 0],zOut=zInitial,latitude=latitude,nModes=reference.nModes,g=g);
    im.upperBoundary = UpperBoundary.rigidLid;
    im.normalization = Normalization.kConstant;
    [~,~,h] = im.modesAtFrequency(omega);

    reference.kind = "constant";
    reference.N2 = @(z) definition.N0*definition.N0*ones(size(z));
    reference.h = h(:);
    reference.zInitial = zInitial;
    reference.f0 = im.f0;
else
    error("HydrostaticRetainedBandGeneralizedQuadrature:UnknownCase", ...
        "Unknown stratification case %s.", definition.name);
end
end

function [VCheb,h] = solveGEPForExample(im,A,B)
if any(isnan(A(:))) || any(isnan(B(:)))
    error("HydrostaticRetainedBandGeneralizedQuadrature:NaNInMatrix", ...
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
    error("HydrostaticRetainedBandGeneralizedQuadrature:NeedMorePoints", ...
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
    error("HydrostaticRetainedBandGeneralizedQuadrature:NeedMorePoints", ...
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

function validationError = validateEvaluator(reference,z)
if reference.kind ~= "chebyshev"
    validationError = 0;
    return
end

reference.im.z = z;
[FReference,GReference] = reference.im.modesAtFrequency(0);
[~,~,F,G] = evaluatedModeMatrices(reference,z);
validationError = max(relativeFrobeniusError(F,FReference),relativeFrobeniusError(G,GReference));
end

function rows = retainedBandRows(reference,r,retainedFraction,budget)
fixedData = modeData(reference,reference.zInitial,"G^{N+1}");
dzGeometric = geometricIncrements(fixedData.z);
rows = diagnosticRow(fixedData,"geometric",dzGeometric,dzGeometric,r,retainedFraction,0,NaN,NaN);

[dzFixed,dzRule,incrementExitflag] = bestRetainedIncrements(fixedData,r);
rows(end+1) = diagnosticRow(fixedData,dzRule,dzFixed,dzGeometric,r,retainedFraction,0,NaN,incrementExitflag);

[zSmooth,nodeExitflag] = optimizeSmoothNodes(reference,r,budget);
smoothData = modeData(reference,zSmooth,"optimized smooth");
[dzSmooth,dzRule,incrementExitflag] = bestRetainedIncrements(smoothData,r);
rows(end+1) = diagnosticRow(smoothData,dzRule,dzSmooth,dzGeometric,r,retainedFraction,max(abs(zSmooth - reference.zInitial)),nodeExitflag,incrementExitflag);

if budget.shouldRunFullNodes && reference.nPoints <= 32
    [zFull,nodeExitflag] = optimizeFullNodes(reference,r,budget);
    fullData = modeData(reference,zFull,"optimized full");
    [dzFull,dzRule,incrementExitflag] = bestRetainedIncrements(fullData,r);
    rows(end+1) = diagnosticRow(fullData,dzRule,dzFull,dzGeometric,r,retainedFraction,max(abs(zFull - reference.zInitial)),nodeExitflag,incrementExitflag);
end
end

function data = modeData(reference,z,nodeRule)
[PhiF,PhiG,F,G] = evaluatedModeMatrices(reference,z);
data.name = reference.name;
data.kind = reference.kind;
data.nodeRule = string(nodeRule);
data.depth = reference.depth;
data.g = reference.g;
data.N2 = reference.N2(z);
data.N2Interior = data.N2(2:end-1);
data.z = z;
data.h = reference.h;
data.PhiF = PhiF(:,1:(reference.nInterior + 1));
data.PhiG = PhiG(:,1:reference.nInterior);
data.F = F;
data.G = G;
end

function [PhiF,PhiG,F,G] = evaluatedModeMatrices(reference,z)
switch reference.kind
    case "chebyshev"
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
    case "constant"
        im = InternalModesConstantStratification(N0=reference.N0,zIn=[-reference.depth 0],zOut=z,latitude=reference.latitude,nModes=reference.nModes,g=reference.g);
        im.upperBoundary = UpperBoundary.rigidLid;
        im.normalization = Normalization.kConstant;
        [F,G] = im.modesAtFrequency(0);
        geostrophicScale = sqrt((reference.N0*reference.N0 - reference.f0*reference.f0)/(reference.N0*reference.N0));
        F = geostrophicScale*F;
        G = geostrophicScale*G;
    otherwise
        error("HydrostaticRetainedBandGeneralizedQuadrature:UnknownKind", ...
            "Unknown reference kind %s.", reference.kind);
end

PhiF = cat(2,ones(length(z),1),F(:,1:reference.nInterior));
PhiG = G(2:end-1,1:reference.nInterior);
end

function r = retainedModeCount(reference,retainedFraction)
r = max(1,floor(retainedFraction*reference.nInterior));
end

function row = diagnosticRow(data,dzRule,dz,dzReference,r,retainedFraction,maxNodeShift,nodeExitflag,incrementExitflag)
row.nodeCase = data.name;
row.nodeRule = data.nodeRule;
row.dzRule = string(dzRule);
row.retainedModes = r;
row.retainedFraction = retainedFraction;
row.dz = dz;
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.relativeDzDistance = norm(dz - dzReference)/norm(dzReference);
row.maxNodeShift = maxNodeShift;
row.nodeExitflag = nodeExitflag;
row.incrementExitflag = incrementExitflag;

f = componentDiagnostics(data.PhiF(:,1:(r+1)),dz,gammaF0(data,r));
g = componentDiagnostics(data.PhiG(:,1:r),dz(2:end-1).*data.N2Interior/data.g,eye(r));
row.fSpecError = f.specError;
row.gSpecError = g.specError;
row.maxSpecError = max(f.specError,g.specError);
row.fFrobeniusError = f.frobeniusError;
row.gFrobeniusError = g.frobeniusError;
row.fGalerkinError = f.galerkinError;
row.gGalerkinError = g.galerkinError;
row.conditionF = f.conditionNumber;
row.conditionG = g.conditionNumber;
end

function diagnostics = componentDiagnostics(Phi,weights,Gamma0)
Gamma = Phi.'*(weights.*Phi);
normalizedDifference = normalizedMatrix(Gamma - Gamma0,Gamma0);
diagnostics.specError = norm(symmetrize(normalizedDifference),2);
diagnostics.frobeniusError = norm(normalizedDifference,"fro");
diagnostics.galerkinError = normalizedOperatorError(Gamma\Gamma,Gamma0);
diagnostics.conditionNumber = cond(Phi);
end

function GammaF0 = gammaF0(data,r)
GammaF0 = diag([data.depth; data.h(1:r)]);
end

function [dz,dzRule,exitflag] = bestRetainedIncrements(data,r)
dzGeometric = geometricIncrements(data.z);
[dzFit,exitflag] = positiveRetainedIncrements(data,r);
if isempty(dzFit)
    dz = dzGeometric;
    dzRule = "geometric fallback";
    exitflag = NaN;
    return
end

geometricError = retainedMaxSpecError(data,dzGeometric,r);
fitError = retainedMaxSpecError(data,dzFit,r);
if fitError <= geometricError*(1 + 1e-8)
    dz = dzFit;
    dzRule = "retained positive LS";
else
    dz = dzGeometric;
    dzRule = "geometric fallback";
    exitflag = -1;
end
end

function value = retainedMaxSpecError(data,dz,r)
f = componentDiagnostics(data.PhiF(:,1:(r+1)),dz,gammaF0(data,r));
g = componentDiagnostics(data.PhiG(:,1:r),dz(2:end-1).*data.N2Interior/data.g,eye(r));
value = max(f.specError,g.specError);
end

function [dz,exitflag] = positiveRetainedIncrements(data,r)
[A,rhs] = normalizedRetainedSystem(data,r);
Aeq = ones(1,length(data.z));
beq = data.depth;
[dz,exitflag] = positiveConstrainedLeastSquares(A,rhs,Aeq,beq);
end

function [A,rhs] = normalizedRetainedSystem(data,r)
nGrid = length(data.z);
PhiF = data.PhiF(:,1:(r+1));
GammaF0 = gammaF0(data,r);
PhiG = data.PhiG(:,1:r);
GammaG0 = eye(r);
nF = size(PhiF,2);
nG = size(PhiG,2);

A = zeros(nF*nF + nG*nG,nGrid);
rhs = zeros(size(A,1),1);
gammaF = diag(GammaF0);
gammaG = diag(GammaG0);
row = 0;

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

function [zOptimized,exitflag] = optimizeSmoothNodes(reference,r,budget)
if exist("fmincon","file") ~= 2
    zOptimized = reference.zInitial;
    exitflag = NaN;
    return
end

yInitial = normalizedCoordinate(reference,reference.zInitial);
nBasis = min(budget.nSmoothModes,max(2,floor(length(yInitial)/4)));
basis = smoothPerturbationBasis(length(yInitial),nBasis);
[coefficients,exitflag] = bestSmoothRun(reference,r,yInitial,basis,budget);
yOptimized = yInitial + basis*coefficients;
zOptimized = zFromNormalizedCoordinate(reference,yOptimized);
end

function [zOptimized,exitflag] = optimizeFullNodes(reference,r,budget)
if exist("fmincon","file") ~= 2
    zOptimized = reference.zInitial;
    exitflag = NaN;
    return
end

yInitial = normalizedCoordinate(reference,reference.zInitial);
[yOptimized,exitflag] = bestFullRun(reference,r,yInitial,budget);
zOptimized = zFromNormalizedCoordinate(reference,yOptimized);
end

function [bestCoefficients,bestExitflag] = bestSmoothRun(reference,r,yInitial,basis,budget)
nCoefficients = size(basis,2);
spacing = min(diff(cat(1,0,yInitial,1)));
amplitude = 0.12*spacing;
starts = zeros(nCoefficients,1 + nCoefficients + budget.nRandomStarts);
for iCoefficient = 1:nCoefficients
    start = zeros(nCoefficients,1);
    start(iCoefficient) = amplitude*(-1)^iCoefficient;
    starts(:,iCoefficient+1) = start;
end
for iStart = 1:budget.nRandomStarts
    starts(:,1+nCoefficients+iStart) = amplitude*randn(nCoefficients,1);
end

options = optimoptions("fmincon", ...
    Algorithm="sqp", ...
    Display="off", ...
    MaxIterations=budget.maxIterations, ...
    MaxFunctionEvaluations=budget.maxFunctionEvaluations, ...
    OptimalityTolerance=1e-7, ...
    StepTolerance=1e-8);
lb = -5*amplitude*ones(nCoefficients,1);
ub = 5*amplitude*ones(nCoefficients,1);
constraint = @(coefficients) smoothCoordinateConstraint(yInitial,basis,coefficients);

bestCoefficients = zeros(nCoefficients,1);
bestObjective = objectiveForSmoothCoefficients(reference,r,yInitial,basis,bestCoefficients);
bestExitflag = NaN;
for iStart = 1:size(starts,2)
    objective = @(coefficients) objectiveForSmoothCoefficients(reference,r,yInitial,basis,coefficients);
    [candidateCoefficients,candidateObjective,exitflag] = fmincon(objective,starts(:,iStart),[],[],[],[],lb,ub,constraint,options);
    if isfinite(candidateObjective) && candidateObjective < bestObjective
        bestCoefficients = candidateCoefficients;
        bestObjective = candidateObjective;
        bestExitflag = exitflag;
    end
end
end

function [bestY,bestExitflag] = bestFullRun(reference,r,yInitial,budget)
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
starts = fullNodeStarts(yInitial);

options = optimoptions("fmincon", ...
    Algorithm="sqp", ...
    Display="off", ...
    MaxIterations=min(8,budget.maxIterations), ...
    MaxFunctionEvaluations=min(90,budget.maxFunctionEvaluations), ...
    OptimalityTolerance=1e-7, ...
    StepTolerance=1e-8);

bestY = yInitial;
bestObjective = objectiveForNormalizedCoordinate(reference,r,yInitial);
bestExitflag = NaN;
for iStart = 1:size(starts,2)
    objective = @(yInterior) objectiveForNormalizedCoordinate(reference,r,yInterior);
    [candidateY,candidateObjective,exitflag] = fmincon(objective,starts(:,iStart),Aineq,bineq,[],[],lb,ub,[],options);
    if isfinite(candidateObjective) && candidateObjective < bestObjective
        bestY = candidateY;
        bestObjective = candidateObjective;
        bestExitflag = exitflag;
    end
end
end

function objective = objectiveForSmoothCoefficients(reference,r,yInitial,basis,coefficients)
y = yInitial + basis*coefficients(:);
objective = objectiveForNormalizedCoordinate(reference,r,y);
end

function objective = objectiveForNormalizedCoordinate(reference,r,yInterior)
if any(diff(cat(1,0,yInterior(:),1)) <= 0) || any(yInterior <= 0) || any(yInterior >= 1)
    objective = realmax;
    return
end

z = zFromNormalizedCoordinate(reference,yInterior);
data = modeData(reference,z,"trial");
[dz,~,~] = bestRetainedIncrements(data,r);
if isempty(dz)
    objective = realmax;
    return
end

f = componentDiagnostics(data.PhiF(:,1:(r+1)),dz,gammaF0(data,r));
g = componentDiagnostics(data.PhiG(:,1:r),dz(2:end-1).*data.N2Interior/data.g,eye(r));
objective = f.specError*f.specError + g.specError*g.specError;
if ~isfinite(objective)
    objective = realmax;
end
end

function basis = smoothPerturbationBasis(nInterior,nBasis)
t = linspace(0,1,nInterior).';
basis = zeros(nInterior,nBasis);
for iBasis = 1:nBasis
    basis(:,iBasis) = sin(iBasis*pi*t);
end
end

function starts = fullNodeStarts(yInitial)
spacing = min(diff(cat(1,0,yInitial,1)));
amplitude = 0.08*spacing;
t = linspace(0,1,length(yInitial)).';
starts = zeros(length(yInitial),5);
starts(:,1) = yInitial;
for iMode = 1:2
    starts(:,2*iMode) = yInitial + amplitude*sin(iMode*pi*t);
    starts(:,2*iMode+1) = yInitial - amplitude*sin(iMode*pi*t);
end
end

function [cineq,ceq] = smoothCoordinateConstraint(yInitial,basis,coefficients)
minSpacing = 1e-6;
y = yInitial + basis*coefficients(:);
allY = cat(1,0,y,1);
cineq = cat(1,minSpacing - diff(allY),minSpacing - y,y - (1 - minSpacing));
ceq = [];
end

function y = normalizedCoordinate(reference,z)
switch reference.kind
    case "chebyshev"
        x = reference.im.x_function(z);
        y = (x(2:end-1) - reference.im.xMin)/reference.im.Lx;
    case "constant"
        y = (z(2:end-1) + reference.depth)/reference.depth;
    otherwise
        error("HydrostaticRetainedBandGeneralizedQuadrature:UnknownKind", ...
            "Unknown reference kind %s.", reference.kind);
end
end

function z = zFromNormalizedCoordinate(reference,yInterior)
switch reference.kind
    case "chebyshev"
        x = reference.im.xMin + reference.im.Lx*cat(1,0,yInterior(:),1);
        z = InternalModesSpectral.fInverseBisection( ...
            reference.im.x_function,x,min(reference.im.zDomain),max(reference.im.zDomain),1e-12);
    case "constant"
        z = -reference.depth + reference.depth*cat(1,0,yInterior(:),1);
    otherwise
        error("HydrostaticRetainedBandGeneralizedQuadrature:UnknownKind", ...
            "Unknown reference kind %s.", reference.kind);
end
end

function matrix = normalizedMatrix(matrix,Gamma0)
scale = diag(1./sqrt(diag(Gamma0)));
matrix = scale*matrix*scale;
end

function value = normalizedOperatorError(operatorProduct,Gamma0)
scale = diag(sqrt(diag(Gamma0)));
inverseScale = diag(1./sqrt(diag(Gamma0)));
value = norm(scale*(operatorProduct - eye(size(operatorProduct)))*inverseScale,2);
end

function matrix = symmetrize(matrix)
matrix = 0.5*(matrix + matrix.');
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

function x = leastSquaresSolution(A,b)
[Q,R,E] = qr(A,0);
tol = max(size(A))*eps(norm(R,inf));
rankA = sum(abs(diag(R)) > tol);
if rankA == 0
    x = zeros(size(A,2),1);
else
    x = E(:,1:rankA)*(R(1:rankA,1:rankA)\(Q(:,1:rankA)'*b));
end
end

function dz = geometricIncrements(z)
dz = zeros(size(z));
dz(1) = 0.5*(z(2) - z(1));
dz(end) = 0.5*(z(end) - z(end-1));
dz(2:end-1) = 0.5*(z(3:end) - z(1:end-2));
end

function printBandGrowthRows(reference,budget)
maxModes = reference.nInterior;
modeCounts = unique([8 16 24 32 floor(maxModes/2) floor(2*maxModes/3) maxModes]);
modeCounts = modeCounts(modeCounts >= 1 & modeCounts <= maxModes);
fixedData = modeData(reference,reference.zInitial,"G^{N+1}");
dzGeometric = geometricIncrements(reference.zInitial);
for iMode = 1:length(modeCounts)
    r = modeCounts(iMode);
    [dzFixed,dzRule,~] = bestRetainedIncrements(fixedData,r);
    fixed = diagnosticRow(fixedData,dzRule,dzFixed,dzGeometric,r,r/maxModes,0,NaN,NaN);
    [zSmooth,nodeExitflag] = optimizeSmoothNodes(reference,r,budget);
    smoothData = modeData(reference,zSmooth,"optimized smooth");
    [dzSmooth,dzRule,incrementExitflag] = bestRetainedIncrements(smoothData,r);
    smooth = diagnosticRow(smoothData,dzRule,dzSmooth,dzGeometric,r,r/maxModes,max(abs(zSmooth-reference.zInitial)),nodeExitflag,incrementExitflag);
    improvement = fixed.maxSpecError/smooth.maxSpecError;
    fprintf('%5d %8.3f %10.3e %10.3e %10.3e %10.3e %10.3f %10.3e\n', ...
        r,r/maxModes,fixed.fSpecError,fixed.gSpecError,smooth.fSpecError,smooth.gSpecError,improvement,smooth.maxNodeShift);
end
end

function printConstantStratificationCheck(reference,r)
data = modeData(reference,reference.zInitial,"G^{N+1}");
dz = geometricIncrements(data.z);
x = (data.z + data.depth)/data.depth;
cosineBasis = cos(x*(0:r)*pi);
sineBasis = sin(x(2:end-1)*(1:r)*pi);
cosineBasis = scaleColumnsToMatch(cosineBasis,data.PhiF(:,1:(r+1)));
sineBasis = scaleColumnsToMatch(sineBasis,data.PhiG(:,1:r));
weightsG = dz(2:end-1).*data.N2Interior/data.g;

fprintf('\nConstant stratification DCT/DST-I check, r = %d\n', r);
fprintf('  F cosine shape error:      %.3e\n', relativeFrobeniusError(cosineBasis,data.PhiF(:,1:(r+1))));
fprintf('  G sine shape error:        %.3e\n', relativeFrobeniusError(sineBasis,data.PhiG(:,1:r)));
fprintf('  F adjoint agreement:       %.3e\n', operatorAgreement(cosineBasis,data.PhiF(:,1:(r+1)),dz,gammaF0(data,r)));
fprintf('  G adjoint agreement:       %.3e\n', operatorAgreement(sineBasis,data.PhiG(:,1:r),weightsG,eye(r)));
end

function value = operatorAgreement(referencePhi,modePhi,weights,Gamma0)
referenceOperator = Gamma0 \ (referencePhi.' .* weights.');
modeOperator = Gamma0 \ (modePhi.' .* weights.');
value = relativeFrobeniusError(referenceOperator,modeOperator);
end

function scaledReference = scaleColumnsToMatch(reference,target)
scaledReference = reference;
for iColumn = 1:size(reference,2)
    scale = (reference(:,iColumn)'*target(:,iColumn))/(reference(:,iColumn)'*reference(:,iColumn));
    scaledReference(:,iColumn) = scale*reference(:,iColumn);
end
end

function printComparisonHeader()
fprintf('%-17s %-20s %5s %6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %9s %9s %9s %6s\n', ...
    'node rule','dz rule','r','frac','F spec','G spec','max spec','F fro','G fro','min dz','max dz','sum dz','rel dz','cond F','cond G','node m','exit');
end

function printComparisonRow(row)
fprintf('%-17s %-20s %5d %6.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %9.2e %9.2e %9.2e %6.0f\n', ...
    char(row.nodeRule),char(row.dzRule),row.retainedModes,row.retainedFraction,row.fSpecError,row.gSpecError, ...
    row.maxSpecError,row.fFrobeniusError,row.gFrobeniusError,row.minDz,row.maxDz,row.sumDz,row.relativeDzDistance, ...
    row.conditionF,row.conditionG,row.maxNodeShift,row.nodeExitflag);
end

function printBandGrowthHeader()
fprintf('%5s %8s %10s %10s %10s %10s %10s %10s\n', ...
    'r','r/N','fixed F','fixed G','smooth F','smooth G','improve','node m');
end

function printSweepHeader()
fprintf('%-12s %7s %6s %5s %10s %10s %10s %10s %10s %10s\n', ...
    'case','nPoints','frac','r','fixed max','smooth max','improve','fixed G','smooth G','node m');
end

function printSweepRow(reference,retainedFraction,fixed,smooth)
improvement = fixed.maxSpecError/smooth.maxSpecError;
fprintf('%-12s %7d %6.3f %5d %10.3e %10.3e %10.3f %10.3e %10.3e %10.3e\n', ...
    char(reference.name),reference.nPoints,retainedFraction,fixed.retainedModes, ...
    fixed.maxSpecError,smooth.maxSpecError,improvement,fixed.gSpecError,smooth.gSpecError,smooth.maxNodeShift);
end

function printRecommendation(reference,rows)
fixedIndex = find([rows.nodeRule] == "G^{N+1}" & [rows.dzRule] ~= "geometric",1);
smoothIndex = find([rows.nodeRule] == "optimized smooth",1);
fixed = rows(fixedIndex);
smooth = rows(smoothIndex);
improvement = (fixed.maxSpecError - smooth.maxSpecError)/fixed.maxSpecError;
fprintf('Recommendation: fixed-node retained fit max E_spec %.3e, smooth generalized max E_spec %.3e.\n', ...
    fixed.maxSpecError,smooth.maxSpecError);
if reference.name == "constant" && smooth.maxSpecError > 1e-10
    fprintf('  Smooth nodes damage the constant-stratification DCT/DST limit; prefer fixed nodes here.\n');
elseif improvement < 0.1
    fprintf('  Smooth node optimization is not materially better than fixed G^{N+1} nodes.\n');
elseif smooth.relativeDzDistance > 0.25
    fprintf('  Smooth nodes improve E_spec but require noticeably different increments.\n');
else
    fprintf('  Smooth generalized nodes materially improve the retained-band Parseval metric.\n');
end
end

function value = relativeFrobeniusError(A,B)
value = norm(A - B,"fro")/norm(B,"fro");
end
