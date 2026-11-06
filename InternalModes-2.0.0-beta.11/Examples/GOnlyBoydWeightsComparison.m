%% G-only Boyd-style quadrature weights
% This example asks what happens if we ignore the hydrostatic F modes and
% try to build the best possible G-mode quadrature on fixed G^{N+1} roots.
%
% For active interior samples
%
%     B_ij = G^j(z_i),
%
% the G-only Parseval target is
%
%     B' diag(q) B = I,      q_i = N2(z_i) dz_i/g.
%
% Boyd's classical Gaussian-quadrature story suggests the Christoffel weight
%
%     q_i = 1/sum_j G^j(z_i)^2.
%
% For variable stratification this is not guaranteed to be exact, so we also
% solve a positive least-squares fit for q_i. Finally we convert q_i back to a
% physical shared increment dz_i = g q_i/N2(z_i). The rigid-lid G modes do
% not constrain endpoint increments, so the remaining endpoint depth gives
% one free split between the bottom and top endpoint. The table reports both
% a symmetric split and the split that minimizes the Parseval-safe F-mode
% error. The highest endpoint F mode belongs to square interpolation
% mechanics and is not included in the reported F spectrum.

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
    printHeader();
    for iRow = 1:length(rows)
        printRow(rows(iRow));
    end
end

fprintf('\nResolution sweep: fixed G^{N+1} nodes\n');
fprintf('%-13s %8s %-20s %10s %10s %10s %10s %10s %10s\n', ...
    'case','nPoints','weight rule','G spec','F spec','F int','max spec','sum dz','min dz');
for iCase = 1:length(cases)
    for iPoint = 1:length(nPointsList)
        data = modeData(cases(iCase),nPointsList(iPoint),latitude,g);
        rows = diagnosticRows(data);
        for iRow = 1:length(rows)
            fprintf('%-13s %8d %-20s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
                data.name,nPointsList(iPoint),rows(iRow).rule,rows(iRow).gSpecError, ...
                rows(iRow).fSpecError,rows(iRow).fIntegralError, ...
                max(rows(iRow).gSpecError,rows(iRow).fSpecError),rows(iRow).sumDz,rows(iRow).minDz);
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
    error("GOnlyBoydWeightsComparison:UnknownCase", "Unknown case %s.", definition.name);
end

[F,G,h] = im.modesAtFrequency(omega);
if definition.name == "constant"
    geostrophicScale = sqrt((definition.N0*definition.N0 - im.f0*im.f0)/(definition.N0*definition.N0));
    F = geostrophicScale*F;
    G = geostrophicScale*G;
end

nInterior = nPoints - 2;
data.z = z;
data.N2 = im.N2(:);
data.N2Interior = data.N2(2:end-1);
data.PhiF = cat(2,ones(nPoints,1),F(:,1:nInterior));
data.PhiFIntegral = cat(2,ones(nPoints,1),F(:,1:(nInterior+1)));
data.PhiG = G(2:end-1,1:nInterior);
data.GammaF0 = diag([definition.depth; h(1:nInterior).']);
data.GammaFNative = diag([definition.depth; h(1:(nInterior+1)).']);
data.GammaFNative(end,end) = 2*data.GammaFNative(end,end);
data.GammaG0 = eye(nInterior);
data.bFIntegral = zeros(size(data.PhiFIntegral,2),1);
data.bFIntegral(1) = definition.depth;
end

function rows = diagnosticRows(data)
qChristoffel = christoffelWeights(data.PhiG);
qPositiveLS = positiveGOnlyWeights(data.PhiG);
dzFChristoffel = fChristoffelIncrements(data);

rows = [
    rowFromDz(data,"geometric dz",geometricIncrements(data.z))
    rowFromDz(data,"F-compatible dz",leastSquaresSolution(data.PhiFIntegral.',data.bFIntegral))
    rowFromDz(data,"F Christoffel dz",dzFChristoffel)
    rowFromDz(data,"joint F/G LS dz",jointLeastSquaresIncrements(data))
    rowFromQ(data,"G Christoffel sym",qChristoffel,"symmetric")
    rowFromQ(data,"G Christoffel Fopt",qChristoffel,"fOptimal")
    rowFromQ(data,"G positive LS sym",qPositiveLS,"symmetric")
    rowFromQ(data,"G positive LS Fopt",qPositiveLS,"fOptimal")
    ];
end

function row = rowFromDz(data,rule,dz)
q = data.N2Interior.*dz(2:end-1)/data.g;
row = diagnosticsFromDzAndQ(data,rule,dz,q);
end

function row = rowFromQ(data,rule,q,endpointRule)
dzInterior = data.g*q./data.N2Interior;
endpointTotal = data.depth - sum(dzInterior);
if endpointRule == "symmetric"
    bottomDz = endpointTotal/2;
elseif endpointRule == "fOptimal"
    bottomDz = bestBottomEndpointDz(data,dzInterior,endpointTotal);
else
    error("GOnlyBoydWeightsComparison:UnknownEndpointRule", "Unknown endpoint rule %s.", endpointRule);
end
topDz = endpointTotal - bottomDz;
dz = [bottomDz; dzInterior; topDz];
row = diagnosticsFromDzAndQ(data,rule,dz,q);
end

function row = diagnosticsFromDzAndQ(data,rule,dz,q)
row.rule = string(rule);
row.gSpecError = specError(data.PhiG.'*(q.*data.PhiG),data.GammaG0);
row.gFroError = norm(data.PhiG.'*(q.*data.PhiG) - data.GammaG0,'fro')/norm(data.GammaG0,'fro');
row.fSpecError = specError(data.PhiF.'*(dz.*data.PhiF),data.GammaF0);
row.fIntegralError = norm(data.PhiFIntegral.'*dz - data.bFIntegral)/norm(data.bFIntegral);
row.sumDz = sum(dz);
row.minDz = min(dz);
row.bottomDz = dz(1);
row.topDz = dz(end);
row.minQ = min(q);
row.maxQ = max(q);
row.relDzToGeometric = norm(dz - geometricIncrements(data.z))/norm(geometricIncrements(data.z));
row.condG = cond(data.PhiG);
row.condF = cond(data.PhiF);
end

function bottomDz = bestBottomEndpointDz(data,dzInterior,endpointTotal)
if endpointTotal <= 0
    bottomDz = endpointTotal/2;
    return
end

objective = @(candidateBottomDz) fSpecForEndpointSplit(data,dzInterior,candidateBottomDz,endpointTotal);
bottomDz = fminbnd(objective,0,endpointTotal);
end

function value = fSpecForEndpointSplit(data,dzInterior,bottomDz,endpointTotal)
dz = [bottomDz; dzInterior; endpointTotal - bottomDz];
value = specError(data.PhiF.'*(dz.*data.PhiF),data.GammaF0);
end

function q = christoffelWeights(PhiG)
q = 1./sum(PhiG.^2,2);
end

function dz = fChristoffelIncrements(data)
scale = diag(1./sqrt(diag(data.GammaFNative)));
normalizedF = data.PhiFIntegral*scale;
dz = 1./sum(normalizedF.^2,2);
dz = data.depth*dz/sum(dz);
end

function q = positiveGOnlyWeights(PhiG)
[A,b] = gOnlyLeastSquaresSystem(PhiG);
try
    options = optimoptions("lsqlin",Display="off");
    q = lsqlin(A,b,[],[],[],[],zeros(size(PhiG,1),1),[],[],options);
catch
    q = leastSquaresSolution(A,b);
end
end

function dz = jointLeastSquaresIncrements(data)
[A,b] = jointLeastSquaresSystem(data);
Aeq = ones(1,length(data.z));
beq = data.depth;
try
    options = optimoptions("lsqlin",Display="off");
    dz = lsqlin(A,b,[],[],Aeq,beq,zeros(length(data.z),1),[],[],options);
catch
    dz = equalityConstrainedLeastSquares(A,b,Aeq,beq);
end
end

function [A,b] = gOnlyLeastSquaresSystem(PhiG)
nGrid = size(PhiG,1);
nModes = size(PhiG,2);
A = zeros(nModes*nModes,nGrid);
b = zeros(nModes*nModes,1);
row = 0;
for iMode = 1:nModes
    for jMode = 1:nModes
        row = row + 1;
        A(row,:) = (PhiG(:,iMode).*PhiG(:,jMode)).';
        b(row) = double(iMode == jMode);
    end
end
end

function [A,b] = jointLeastSquaresSystem(data)
nGrid = length(data.z);
nF = size(data.PhiF,2);
nG = size(data.PhiG,2);
A = zeros(nF*nF + nG*nG,nGrid);
b = zeros(size(A,1),1);
row = 0;

fScale = sqrt(diag(data.GammaF0));
for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        A(row,:) = (data.PhiF(:,iMode).*data.PhiF(:,jMode)./(fScale(iMode)*fScale(jMode))).';
        b(row) = double(iMode == jMode);
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(nGrid,1);
        weights(2:end-1) = data.N2Interior.*data.PhiG(:,iMode).*data.PhiG(:,jMode)/data.g;
        A(row,:) = weights.';
        b(row) = double(iMode == jMode);
    end
end
end

function value = specError(Gamma,Gamma0)
scale = diag(1./sqrt(diag(Gamma0)));
value = norm(scale*(Gamma - Gamma0)*scale,2);
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2)-z(1))/2,(z(3:end)-z(1:end-2))/2,(z(end)-z(end-1))/2);
end

function x = equalityConstrainedLeastSquares(A,b,Aeq,beq)
x0 = Aeq.'*((Aeq*Aeq.')\beq);
Z = null(Aeq);
x = x0 + Z*leastSquaresSolution(A*Z,b - A*x0);
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

function printHeader()
fprintf('%-20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'weight rule','G spec','G fro','F spec','F int','sum dz','min dz','bot dz','top dz','rel dz','min q','max q');
end

function printRow(row)
fprintf('%-20s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    row.rule,row.gSpecError,row.gFroError,row.fSpecError,row.fIntegralError,row.sumDz, ...
    row.minDz,row.bottomDz,row.topDz,row.relDzToGeometric,row.minQ,row.maxQ);
end
