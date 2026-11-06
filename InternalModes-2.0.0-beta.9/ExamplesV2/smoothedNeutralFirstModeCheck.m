function result = smoothedNeutralFirstModeCheck(options)
% Check the neutral first mode against independent piecewise shooting.
%
% This supplements bandwidthCoordinateComparison(includeSmoothed=true).
% The linear cap and constant-N deep solution are exact; ode113 integrates
% only the quintic ramp. Two ODE tolerances check the shooting reference.
% This first-mode diagnostic does not replace the primary spectral-reference
% gate or establish a complete resolved-mode count.
arguments (Input)
    options.outputDirectory (1,1) string = ""
    options.smoothingLength (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 147.77523446738553
end
arguments (Output)
    result (1,1) struct
end
outputDirectory = options.outputDirectory;
if outputDirectory == "", outputDirectory = string(tempname)+"-neutral-first-mode"; end
outputNames = ["comparison","measurements","validation","provenance"];
for name = outputNames
    extension = ".csv";
    if name == "provenance", extension = ".json"; end
    if isfile(fullfile(outputDirectory,"neutral-first-mode-shooting-"+name+extension))
        error("IMExample:OutputExists","Choose an outputDirectory without previous first-mode shooting results.");
    end
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
g = 9.81;
evp = IMInternalModes.hydrostaticGModes(N2=@neutralN2,zDomain=[-1000 0],rotationRate=7.2921e-5,latitude=45,g=g,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
rule = IMSolverSpectral(nEVP=1025,coordinateKind="z").configuredForEVP(evp);
[z,w] = rule.nativeQuadratureRule(evp.zDomain);
low = shootingMode(z,1e-11,g);
high = shootingMode(z,1e-13,g);
identity = struct(family=string(evp.name),columnLabels="1",normalization="rmsF",zDomain=evp.zDomain);
lo = prepared(low,identity,"piecewise-shooting-1e-11");
hi = prepared(high,identity,"piecewise-shooting-1e-13");
validation = assessModeConvergence(lo,hi,z,w);
rows = cell(0,1);
allMeasurements = cell(0,1);
% Eight bounded diagnostic rows, independent of the full mode inventory.
for coordinate = ["z","sqrtB","smoothedB"]
    if coordinate == "z", resolutions = [64 128 512 1024]; else, resolutions = [64 128]; end
    for n = resolutions
        if coordinate == "smoothedB"
            solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=n,smoothingLength=options.smoothingLength);
        elseif coordinate == "sqrtB"
            solver = IMSolverSpectral(nEVP=n,coordinateKind="wkb");
        else
            solver = IMSolverSpectral(nEVP=n,coordinateKind="z");
        end
        b = solver.solveEVP(evp,nModes=1);
        scale = b.normalizationFactors();
        coefficients = b.nativeModes/scale(1);
        data = struct(identity=identity,values=struct(F=b.F(z),G=b.G(z)), ...
            derivatives=struct(F=b.h(1)*b.solver.evaluatePhysicalDerivative(coefficients,z,2),G=b.solver.evaluatePhysicalDerivative(coefficients,z,1)), ...
            equivalentDepths=b.h,eigenvalues=b.eigenvalues,provenance=struct(coordinate=coordinate,nEVP=n));
        report = assessModeConvergence(data,hi,z,w);
        measurements = report.measurements;
        measurements = addvars(measurements,repmat(coordinate,height(measurements),1),repmat(n,height(measurements),1),Before=1,NewVariableNames=["coordinate","nEVP"]);
        allMeasurements{end+1,1} = measurements; %#ok<AGROW>
        m = report.measurements;
        rows{end+1,1} = struct(coordinate=coordinate,nEVP=n,hRelativeError=m.value(m.quantity=="equivalentDepth"),FH1=m.value(m.quantity=="h1" & m.variable=="F"),GH1=m.value(m.quantity=="h1" & m.variable=="G")); %#ok<AGROW>
    end
end
summary = struct2table(vertcat(rows{:}));
writetable(summary,fullfile(outputDirectory,"neutral-first-mode-shooting-comparison.csv"));
writetable(vertcat(allMeasurements{:}),fullfile(outputDirectory,"neutral-first-mode-shooting-measurements.csv"));
writetable(validation.measurements,fullfile(outputDirectory,"neutral-first-mode-shooting-validation.csv"));
provenance = struct(smoothingLength=options.smoothingLength,method="Exact linear neutral cap and constant-N deep solution; ode113 on quintic ramp; bracketed first-mode shooting in c=sqrt(g*h)", ...
    cBracket=[1.4 1.6],relativeTolerances=[1e-11 1e-13],shootingH=high.h, ...
    derivativeRule="Gz from integrated state; F=h*Gz; Fz=-N2*G/g exactly", ...
    normalization="rmsF: integral F^2/D=1; hydrostatic identity converts the independent N2*G^2/g integral from augmented ODE and analytic deep layer", ...
    referenceIndependent="No spectral eigenvalues or eigenvectors used to compute shooting mode; bracket chosen for the first positive mode, checked for no interior zeros", ...
    noInteriorZeros=all(high.G(2:end-1)>0),bottomResidual=high.bottomResidual, ...
    quadrature="Same 1025-point physical-z rule as primary study; quadrature convergence unverified", ...
    scope="Supplemental first-mode diagnostic only; does not replace primary 384/512 then 512/1024 reference gate or certify a complete mode count");
fid = fopen(fullfile(outputDirectory,"neutral-first-mode-shooting-provenance.json"),"w");
if fid == -1, error("IMExample:WriteFailed","Cannot write first-mode shooting provenance in %s.",outputDirectory); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",jsonencode(provenance,PrettyPrint=true));
result = struct(summary=summary,validation=validation,provenance=provenance,z=z,weights=w,shooting=high);
disp(summary)
disp(validation.measurements)
end

function data = prepared(mode,identity,method)
data = struct(identity=identity,values=struct(F=mode.F,G=mode.G), ...
    derivatives=struct(F=mode.Fz,G=mode.Gz),equivalentDepths=mode.h,eigenvalues=1/mode.h, ...
    provenance=struct(method=method));
end

function mode = shootingMode(z,tolerance,g)
options = odeset(RelTol=tolerance,AbsTol=[tolerance tolerance/200 tolerance]);
c = fzero(@(c) bottomValue(c,g,options),[1.4 1.6],optimset(TolX=5e-14));
h = c^2/g;
solution = ode113(@(d,y) rampEquation(d,y,c,g),[0 200],[200;1;0],options);
endValues = deval(solution,200);
A = endValues(1);
m = .005/c;
C = endValues(2)/m;
deepLength = 600;
deepIntegral = A^2*(deepLength/2+sin(2*m*deepLength)/(4*m)) ...
    + C^2*(deepLength/2-sin(2*m*deepLength)/(4*m)) ...
    + A*C*sin(m*deepLength)^2/m;
normalization = sqrt((h/1000)*(endValues(3)+.005^2/g*deepIntegral));
G = zeros(size(z));
Gd = zeros(size(z));
cap = z >= -200;
ramp = z > -400 & z < -200;
deep = z <= -400;
G(cap) = -z(cap);
Gd(cap) = 1;
y = deval(solution,-z(ramp)-200);
G(ramp) = y(1,:).';
Gd(ramp) = y(2,:).';
s = -z(deep)-400;
G(deep) = A*cos(m*s)+C*sin(m*s);
Gd(deep) = -m*A*sin(m*s)+m*C*cos(m*s);
mode = struct(h=h,G=G/normalization,Gz=-Gd/normalization, ...
    F=-h*Gd/normalization,Fz=-neutralN2(z).*G/(g*normalization), ...
    bottomResidual=(A*cos(m*deepLength)+C*sin(m*deepLength))/normalization);
end

function value = bottomValue(c,g,options)
solution = ode113(@(d,y) rampEquation(d,y,c,g),[0 200],[200;1;0],options);
y = deval(solution,200);
m = .005/c;
value = y(1)*cos(m*600)+y(2)*sin(m*600)/m;
end

function dy = rampEquation(d,y,c,g)
t = d/200;
N2 = .005^2*(6*t^5-15*t^4+10*t^3);
dy = [y(2);-N2*y(1)/c^2;N2*y(1)^2/g];
end

function values = neutralN2(z)
t = min(max((-z-200)/200,0),1);
values = .005^2*(6*t.^5-15*t.^4+10*t.^3);
end
