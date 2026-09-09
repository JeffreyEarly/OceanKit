%% Evaluate physical-sampling quadrature for hydrostatic modes
% This example keeps the physical stratification value N2(z_i) fixed and asks
% how well one shared grid z_i and one shared positive increment vector dz_i
% can support both F and G Parseval identities.
%
% The diagnostics are
%
%     E_F   = ||PhiF' WF PhiF - GammaF|| / ||GammaF||,
%     E_G   = ||PhiG' WG PhiG - I|| / ||I||,
%     E_int = ||PhiF' dz - [D;0;...]|| / D,
%
% with WF = diag(dz) and WG = diag(N2(z_i) dz_i/g).

%% Single-resolution comparison
baseNPoints = 64;
baseRows = compareQuadratureRules(baseNPoints,true);

fprintf('\nPhysical N2 quadrature comparison at nPoints = %d\n', baseNPoints);
printComparisonHeader();
for iRow = 1:length(baseRows)
    printComparisonRow(baseRows(iRow));
end

%% Resolution dependence on the mode-adapted G grid
nPointsList = [24 32 48 64 96 128];

fprintf('\nResolution dependence with physical N2 on G^{N+1} nodes\n');
fprintf('%8s %-15s %12s %12s %12s %12s %12s %8s %12s\n', ...
    'nPoints','dz rule','F Gram','G Gram','F int','sum dz','min dz','neg dz','J');
for iPoint = 1:length(nPointsList)
    rows = compareQuadratureRules(nPointsList(iPoint),false);
    for iRow = 1:length(rows)
        if rows(iRow).nodeRule == "G^{N+1}" && any(rows(iRow).dzRule == ["geometric","F-compatible","joint LS"])
            fprintf('%8d %-15s %12.6e %12.6e %12.6e %12.6e %12.6e %8d %12.6e\n', ...
                nPointsList(iPoint),rows(iRow).dzRule,rows(iRow).fGramError,rows(iRow).gGramError, ...
                rows(iRow).fIntegralError,rows(iRow).sumDz,rows(iRow).minDz,rows(iRow).negativeDzCount,rows(iRow).objective);
        end
    end
end

%% Local helpers
function rows = compareQuadratureRules(nPoints,shouldIncludeAlternates)
data = modeData(nPoints,"G^{N+1}");
rows = diagnosticRowsForData(data);

if shouldIncludeAlternates
    chebyshevData = modeData(nPoints,"Chebyshev-Lobatto");
    rows = [rows diagnosticRowsForData(chebyshevData)];
end
end

function rows = diagnosticRowsForData(data)
dzGeometric = geometricIncrements(data.z);
dzF = leastSquaresSolution(data.PhiF.',data.bF);
dzJoint = jointLeastSquaresIncrements(data);

rows = rowDiagnostics(data,"geometric",dzGeometric);
rows(end+1) = rowDiagnostics(data,"F-compatible",dzF);
rows(end+1) = rowDiagnostics(data,"joint LS",dzJoint);

if exist("lsqlin","file") == 2
    dzPositive = positiveLeastSquaresIncrements(data);
    if ~isempty(dzPositive)
        rows(end+1) = rowDiagnostics(data,"positive joint LS",dzPositive);
    end
end
end

function data = modeData(nPoints,nodeRule)
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

if nodeRule == "G^{N+1}"
    zReference = linspace(zDomain(1),zDomain(2),1024).';
    imReference = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zReference,latitude=latitude,nEVP=nEVP,g=g);
    imReference.normalization = Normalization.geostrophic;
    imReference.upperBoundary = UpperBoundary.rigidLid;
    z = imReference.GaussQuadraturePointsForModesAtFrequency(nPoints,omega);
elseif nodeRule == "Chebyshev-Lobatto"
    theta = pi*(0:(nPoints-1)).'/(nPoints-1);
    z = -Lz*(1 + cos(theta))/2;
else
    error("PhysicalSamplingQuadratureOptimization:InvalidNodeRule", ...
        "nodeRule must be G^{N+1} or Chebyshev-Lobatto.");
end

im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=latitude,nEVP=nEVP,nModes=nModes,g=g);
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;

[F,G,h] = im.modesAtFrequency(omega);
Nz = length(z);

data.nodeRule = string(nodeRule);
data.Lz = Lz;
data.g = g;
data.z = z;
data.N2 = N2(z);
data.N2Interior = N2(z(2:end-1));
data.PhiF = cat(2,ones(Nz,1),F);
data.PhiG = G(2:end-1,1:end-1);
data.GammaF = diag([Lz; h(:)]);
data.GammaG = eye(size(data.PhiG,2));
data.bF = zeros(Nz,1);
data.bF(1) = Lz;
end

function dz = jointLeastSquaresIncrements(data)
[A,rhs] = jointLeastSquaresSystem(data);
dz = leastSquaresSolution(A,rhs);
end

function dz = positiveLeastSquaresIncrements(data)
[A,rhs] = jointLeastSquaresSystem(data);
try
    options = optimoptions("lsqlin",Display="off");
    dz = lsqlin(A,rhs,[],[],[],[],zeros(length(data.z),1),[],[],options);
catch
    dz = [];
end
end

function [A,rhs] = jointLeastSquaresSystem(data)
Nz = length(data.z);
nF = size(data.PhiF,2);
nG = size(data.PhiG,2);
fScale = norm(data.GammaF,'fro');
gScale = norm(data.GammaG,'fro');
intScale = data.Lz;

nRows = nF*nF + nG*nG + nF;
A = zeros(nRows,Nz);
rhs = zeros(nRows,1);
row = 0;

for iMode = 1:nF
    for jMode = 1:nF
        row = row + 1;
        A(row,:) = (data.PhiF(:,iMode).*data.PhiF(:,jMode)).'/fScale;
        rhs(row) = data.GammaF(iMode,jMode)/fScale;
    end
end

for iMode = 1:nG
    for jMode = 1:nG
        row = row + 1;
        weights = zeros(Nz,1);
        weights(2:end-1) = (data.N2Interior/data.g).*data.PhiG(:,iMode).*data.PhiG(:,jMode);
        A(row,:) = weights.'/gScale;
        rhs(row) = data.GammaG(iMode,jMode)/gScale;
    end
end

for iMode = 1:nF
    row = row + 1;
    A(row,:) = data.PhiF(:,iMode).'/intScale;
    rhs(row) = data.bF(iMode)/intScale;
end
end

function row = rowDiagnostics(data,dzRule,dz)
row.nodeRule = data.nodeRule;
row.dzRule = string(dzRule);
row.fGramError = relativeFrobeniusError(data.PhiF.'*(dz.*data.PhiF),data.GammaF);
row.gGramError = relativeFrobeniusError(data.PhiG.'*((dz(2:end-1).*data.N2Interior/data.g).*data.PhiG),data.GammaG);
row.fIntegralError = norm(data.PhiF.'*dz - data.bF)/norm(data.bF);
row.sumDz = sum(dz);
row.minDz = min(dz);
row.maxDz = max(dz);
row.negativeDzCount = sum(dz < 0);
row.objective = row.fGramError^2 + row.gGramError^2 + row.fIntegralError^2;
end

function dz = geometricIncrements(z)
dz = cat(1,(z(2)-z(1))/2,(z(3:end)-z(1:end-2))/2,(z(end)-z(end-1))/2);
end

function printComparisonHeader()
fprintf('%-18s %-18s %12s %12s %12s %12s %12s %12s %8s %12s\n', ...
    'node rule','dz rule','F Gram','G Gram','F int','sum dz','min dz','max dz','neg dz','J');
end

function printComparisonRow(row)
fprintf('%-18s %-18s %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %8d %12.6e\n', ...
    row.nodeRule,row.dzRule,row.fGramError,row.gGramError,row.fIntegralError, ...
    row.sumDz,row.minDz,row.maxDz,row.negativeDzCount,row.objective);
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
