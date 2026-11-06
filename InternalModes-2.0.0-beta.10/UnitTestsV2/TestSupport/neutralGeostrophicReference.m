function reference = neutralGeostrophicReference(problem,options)
% Build independent shooting references for the observable balanced families.
% Positive roots use regular first-order shooting. Negative roots use local
% transfer matrices and a multiple-shooting matching system, so exponentially
% small opposite-end traces do not determine the reconstructed mode.
arguments
    problem
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 40
    options.joins (1,:) double {mustBeReal,mustBeFinite} = [-1000 -400 -200 0]
    options.tolerance (1,1) double {mustBePositive} = 1e-11
    options.rootBrackets (:,2) double = zeros(0,2)
    options.includeSignedSectors (1,1) logical = false
end
configuration = describeProblem(problem,options);
[zq,wq] = neutralGeostrophicQuadrature(configuration.joins,513);
configuration.integratedN2 = sum(wq.*configuration.N2(zq));
configuration.maximumN2 = max(configuration.N2(zq));
configuration.phaseLength = sum(wq.*sqrt(configuration.N2(zq)/configuration.g));
if configuration.maximumN2 == 0
    error("NeutralGeostrophicReference:FullyNeutralReference","Use the exact finite-family tests for identically zero stratification.");
end
if configuration.family == "zeroAPV"
    reference = zeroAPVReference(configuration,zq,wq);
    return
end
if isempty(options.rootBrackets)
    [negativeCount,hasZero] = rootInertia(configuration);
    nRoots = options.nModes;
    if options.includeSignedSectors && configuration.family ~= "gpe"
        nRoots = nRoots+negativeCount+double(hasZero);
    end
    [roots,brackets,scan] = findRoots(configuration,nRoots);
else
    nRoots = options.nModes;
    brackets = options.rootBrackets;
    roots = zeros(size(brackets,1),1);
    for i = 1:size(brackets,1)
        roots(i) = refineRoot(configuration,brackets(i,:));
    end
    scan = struct("reusedIndependentBrackets",true);
end
[roots,order] = sort(roots);
brackets = brackets(order,:);
roots = roots(1:min(nRoots,numel(roots)));
brackets = brackets(1:numel(roots),:);
profiles = cell(1,numel(roots));
for i = 1:numel(roots)
    profiles{i} = makeProfile(configuration,roots(i));
end
raw = evaluateProfiles(configuration,profiles,zq);
factors = normalizationFactors(configuration,raw,wq);
reference.family = configuration.family;
reference.inverseDepth = roots.';
reference.h = 1./reference.inverseDepth;
reference.modeNumber = rootLabels(configuration,roots);
reference.rootBrackets = brackets;
reference.scan = scan;
reference.tolerance = options.tolerance;
reference.joins = configuration.joins;
reference.evaluate = @(z) normalizedProfiles(configuration,profiles,factors,z);
reference.eigenvalues = physicalEigenvalues(configuration,raw,wq,roots);
reference.configuration = configuration;
end

function configuration = describeProblem(problem,options)
configuration.N2 = problem.N2;
configuration.zDomain = problem.zDomain;
configuration.g = problem.g;
configuration.tolerance = options.tolerance;
configuration.joins = unique([problem.zDomain options.joins]);
configuration.joins = configuration.joins(configuration.joins >= problem.zDomain(1) & configuration.joins <= problem.zDomain(2));
configuration.sigma = 1;
configuration.k = 0;
configuration.f0 = problem.f0;
if isa(problem,"IMGeostrophicZeroAPVModes")
    if numel(problem.k) ~= 1
        error("NeutralGeostrophicReference:ScalarWavenumber","Build one zero-APV reference per distinct wavenumber.");
    end
    configuration.family = "zeroAPV";
    configuration.k = problem.k;
    configuration.sigma = double(problem.surfaceBoundary == "freeSurface");
elseif string(problem.name) == "geostrophicAPVModes"
    configuration.family = "apv";
    configuration.g0 = problem.parameters.g0;
    configuration.gd = problem.parameters.gd;
    configuration.sigma = double(string(problem.parameters.surfaceBoundary) == "freeSurface");
elseif string(problem.name) == "geostrophicGeneralizedPotentialEnstrophyModes"
    configuration.family = "gpe";
    configuration.k = problem.parameters.k;
    configuration.alpha0 = problem.parameters.alpha0;
    configuration.alphaD = problem.parameters.alphaD;
elseif isa(problem,"IMMeanDensityAnomalyModes")
    configuration.family = "mda";
    configuration.g0 = problem.g0;
    configuration.gd = problem.gd;
else
    error("NeutralGeostrophicReference:UnsupportedFamily","No independent balanced reference is defined for %s.",class(problem));
end
end

function [roots,brackets,scan] = findRoots(c,nModes)
[expectedNegative,expectedZero] = rootInertia(c);
roots = zeros(0,1);
brackets = zeros(0,2);
% Endpoint inertia certifies the number of negative branches. Refine the
% logarithmic search if two endpoint roots initially occupy one scan cell.
negativeScale = 1/max(c.phaseLength^2,realmin);
negativeMaximum = 100*negativeScale;
if c.family == "gpe"
    negativeMaximum = c.g*c.k^2/c.f0^2;
else
    accelerations = abs([c.g0 c.gd]);
    accelerations = accelerations(isfinite(accelerations) & accelerations > 0);
    if ~isempty(accelerations)
        negativeMaximum = max(negativeMaximum,100*c.g*c.maximumN2/min(accelerations)^2);
    end
end
negativeRoots = zeros(0,1);
negativeBrackets = zeros(0,2);
nScan = 0;
if expectedNegative > 0
    for nScan = [161 321 641 1281]
        negativeGrid = -exp(linspace(log(negativeMaximum),log(1e-12*negativeScale),nScan));
        [negativeRoots,negativeBrackets] = scanRoots(c,negativeGrid);
        if numel(negativeRoots) == expectedNegative
            break
        end
    end
end
if numel(negativeRoots) ~= expectedNegative
    error("NeutralGeostrophicReference:UnresolvedNegativeBranches", ...
        "Independent shooting found %d negative roots; endpoint inertia requires %d.",numel(negativeRoots),expectedNegative);
end
roots = [roots;negativeRoots];
brackets = [brackets;negativeBrackets];
if expectedZero
    roots(end+1,1) = 0;
    brackets(end+1,:) = [0 0];
end
neededPositive = max(0,nModes-numel(roots));
tStep = pi/(8*sqrt(c.maximumN2/c.g)*diff(c.zDomain));
tMaximum = (neededPositive+5)*pi/c.phaseLength;
positiveGrid = (linspace(1e-6*tStep,tMaximum,ceil(tMaximum/tStep)+1)).^2;
[positiveRoots,positiveBrackets] = scanRoots(c,positiveGrid);
if numel(positiveRoots) < neededPositive
    error("NeutralGeostrophicReference:InsufficientPositiveRoots", ...
        "Independent shooting found %d positive roots; %d are required.",numel(positiveRoots),neededPositive);
end
roots = [roots;positiveRoots(1:neededPositive)];
brackets = [brackets;positiveBrackets(1:neededPositive,:)];
scan = struct("negativeExpected",expectedNegative,"negativeFound",numel(negativeRoots), ...
    "zeroExpected",expectedZero,"positiveFound",numel(positiveRoots),"negativeSamples",nScan, ...
    "positiveSamples",numel(positiveGrid),"reusedIndependentBrackets",false);
end

function [negativeCount,hasZero] = rootInertia(c)
I = c.integratedN2;
if c.family == "mda"
    negativeCount = double(c.g0 < 0)+double(c.gd < 0);
    hasZero = isfinite(c.g0) && isfinite(c.gd);
    if hasZero
        constantNorm = I+c.g0+c.gd;
        if abs(constantNorm) <= 1e3*eps*(I+abs(c.g0)+abs(c.gd))
            error("NeutralGeostrophicReference:NullConstantMode","The independent MDA reference has a null constant mode.");
        end
        negativeCount = negativeCount-double(constantNorm < 0);
    end
    return
end
if c.family == "apv"
    beta0 = c.sigma/c.g+1/c.g0;
    betaD = 1/c.gd;
    K = [1+I*beta0 -1;-1 1+I*betaD];
    retained = [c.g0 ~= 0 c.gd ~= 0];
    K = K(retained,retained);
else
    K = [1+(I/c.g)*(1-c.g*c.k^2/c.alpha0) -1;-1 1-I*c.k^2/c.alphaD];
end
if isempty(K)
    negativeCount = 0;
    hasZero = false;
else
    eigenvalues = eig(K);
    threshold = 1e3*eps*norm(K,2);
    negativeCount = nnz(eigenvalues < -threshold);
    hasZero = any(abs(eigenvalues) <= threshold);
end
end

function [roots,brackets] = scanRoots(c,grid)
values = zeros(size(grid));
for i = 1:numel(grid)
    values(i) = endpointDeterminant(c,grid(i));
end
roots = zeros(0,1);
brackets = zeros(0,2);
for i = 1:numel(grid)-1
    if sign(values(i))*sign(values(i+1)) < 0
        bracket = grid(i:i+1);
        root = refineRoot(c,bracket);
        if isempty(roots) || min(abs(root-roots)) > 1e-8*max(abs(root),realmin)
            roots(end+1,1) = root; %#ok<AGROW>
            brackets(end+1,:) = bracket; %#ok<AGROW>
        end
    end
end
end

function root = refineRoot(c,bracket)
if bracket(1) == bracket(2)
    root = bracket(1);
else
    % Solve in a dimensionless local parameter so fzero's absolute stopping
    % rule never sets the accuracy of a small inverse depth.
    scale = max(abs(bracket));
    root = scale*fzero(@(t) endpointDeterminant(c,scale*t),bracket/scale,optimset("TolX",2e-14));
end
end

function determinant = endpointDeterminant(c,r)
if r < 0
    [diagonal,lower,upper] = matchingSystem(c,r);
    previous = 1;
    current = diagonal(1);
    for i = 2:numel(diagonal)
        next = diagonal(i)*current-lower(i-1)*upper(i-1)*previous;
        scale = max([abs(next) abs(current) realmin]);
        previous = current/scale;
        current = next/scale;
    end
    determinant = current;
else
    initial = bottomState(c,r);
    state = initial/norm(initial);
    for i = 1:numel(c.joins)-1
        state = propagatePair(c,r,c.joins(i:i+1),state);
        state = state/max(norm(state),realmin);
    end
    row = surfaceRow(c,r);
    determinant = (row*state)/max(norm(row),realmin);
end
end

function initial = bottomState(c,r)
if c.family == "mda"
    if isinf(c.gd)
        initial = [0;1];
    else
        initial = [1;-c.gd*r/c.g];
    end
elseif c.family == "apv"
    if isinf(c.gd)
        initial = [1;0];
    else
        initial = [c.gd;-c.g];
    end
else
    if isinf(c.alphaD)
        initial = [1;0];
    else
        initial = [c.alphaD;c.g*c.k^2+c.f0^2*r];
    end
end
end

function row = surfaceRow(c,r)
if c.family == "mda"
    if isinf(c.g0)
        row = [1 0];
    else
        row = [-c.g0*r/c.g 1];
    end
elseif c.family == "apv"
    if isinf(c.g0)
        row = [-c.sigma 1];
    else
        row = [-(c.g+c.sigma*c.g0) c.g0];
    end
else
    if isinf(c.alpha0)
        row = [-1 1];
    else
        row = [c.g*c.k^2+c.f0^2*r-c.alpha0 c.alpha0];
    end
end
end

function final = propagatePair(c,r,interval,initial)
a = c.N2([interval(1);mean(interval);interval(2)])/c.g;
if all(a == a(1))
    if c.family == "mda"
        A = [0 1;-r*a(1) 0];
    else
        A = [0 -a(1);r 0];
    end
    final = constantTransfer(A,diff(interval))*initial;
else
    odeOptions = odeset("RelTol",c.tolerance,"AbsTol",c.tolerance*1e-2);
    [~,values] = ode113(@(z,y) pairDerivative(c,r,z,y),interval,initial,odeOptions);
    final = values(end,:).';
end
end

function derivative = pairDerivative(c,r,z,y)
if c.family == "mda"
    derivative = [y(2);-r*c.N2(z)*y(1)/c.g];
else
    derivative = [-c.N2(z)*y(2)/c.g;r*y(1)];
end
end

function T = constantTransfer(A,L)
q = -A(1,2)*A(2,1);
if q == 0
    T = eye(2)+L*A;
elseif q > 0
    frequency = sqrt(q);
    T = cos(frequency*L)*eye(2)+(sin(frequency*L)/frequency)*A;
else
    frequency = sqrt(-q);
    T = cosh(frequency*L)*eye(2)+(sinh(frequency*L)/frequency)*A;
end
end

function [diagonal,lower,upper,nodes,transfers] = matchingSystem(c,r)
% Exact local Dirichlet-to-Neumann matching for G''=-r*N2*G/g.
maximumLength = min(25,4/max(sqrt(abs(r)*c.maximumN2/c.g),realmin));
nodes = c.joins(1);
for i = 1:numel(c.joins)-1
    interval = c.joins(i:i+1);
    zi = linspace(interval(1),interval(2),ceil(diff(interval)/maximumLength)+1);
    nodes = [nodes zi(2:end)]; %#ok<AGROW>
end
n = numel(nodes);
diagonal = zeros(n,1);
lower = zeros(n-1,1);
upper = zeros(n-1,1);
transfers = cell(1,n-1);
for i = 1:n-1
    interval = nodes(i:i+1);
    T = gradientTransfer(c,r,interval);
    transfers{i} = T;
    diagonal(i) = diagonal(i)+T(1,1)/T(1,2);
    diagonal(i+1) = diagonal(i+1)+T(2,2)/T(1,2);
    lower(i) = -1/T(1,2);
    upper(i) = lower(i);
end
if c.family == "mda"
    if isfinite(c.gd)
        diagonal(1) = diagonal(1)-r*c.gd/c.g;
    else
        diagonal(1) = 1;
        upper(1) = 0;
    end
    if isfinite(c.g0)
        diagonal(end) = diagonal(end)-r*c.g0/c.g;
    else
        diagonal(end) = 1;
        lower(end) = 0;
    end
elseif c.family == "apv"
    if isinf(c.gd)
        bottomCoefficient = 0;
    else
        bottomCoefficient = c.g/c.gd;
    end
    if c.gd ~= 0
        if bottomCoefficient == 0
            diagonal(1) = 1;
            upper(1) = 0;
        else
            diagonal(1) = diagonal(1)-r/bottomCoefficient;
        end
    end
    if c.g0 == 0
        surfaceCoefficient = Inf;
    else
        surfaceCoefficient = c.sigma+c.g/c.g0;
    end
    if surfaceCoefficient == 0
        diagonal(end) = 1;
        lower(end) = 0;
    else
        diagonal(end) = diagonal(end)-r/surfaceCoefficient;
    end
else
    LambdaG = c.g*c.k^2+c.f0^2*r;
    bottomCoefficient = -LambdaG/c.alphaD;
    surfaceCoefficient = 1-LambdaG/c.alpha0;
    if isinf(c.alphaD)
        diagonal(1) = 1;
        upper(1) = 0;
    else
        diagonal(1) = bottomCoefficient*diagonal(1)-r;
        upper(1) = bottomCoefficient*upper(1);
    end
    diagonal(end) = surfaceCoefficient*diagonal(end)-r;
    lower(end) = surfaceCoefficient*lower(end);
end
rowScale = max([abs(diagonal) [abs(upper);0] [0;abs(lower)]],[],2);
rowScale = max(rowScale,realmin);
diagonal = diagonal./rowScale;
upper = upper./rowScale(1:end-1);
lower = lower./rowScale(2:end);
end

function T = gradientTransfer(c,r,interval)
a = c.N2([interval(1);mean(interval);interval(2)])/c.g;
if all(a == a(1))
    T = constantTransfer([0 1;-r*a(1) 0],diff(interval));
else
    odeOptions = odeset("RelTol",c.tolerance,"AbsTol",c.tolerance*1e-2);
    initial = reshape(eye(2),4,1);
    [~,values] = ode113(@(z,y) reshape([0 1;-r*c.N2(z)/c.g 0]*reshape(y,2,2),4,1),interval,initial,odeOptions);
    T = reshape(values(end,:),2,2);
end
end

function profile = makeProfile(c,r)
profile.r = r;
if r < 0
    [d,l,u,nodes,transfers] = matchingSystem(c,r);
    matrix = diag(d)+diag(l,-1)+diag(u,1);
    [~,~,V] = svd(matrix,"econ");
    G = V(:,end);
    profile.nodes = nodes;
    profile.segments = cell(1,numel(nodes)-1);
    for i = 1:numel(nodes)-1
        T = transfers{i};
        P = (G(i+1)-T(1,1)*G(i))/T(1,2);
        profile.segments{i} = denseGradientSegment(c,r,nodes(i:i+1),[G(i);P]);
    end
    profile.kind = "multipleShooting";
else
    initial = bottomState(c,r);
    initial = initial/max(norm(initial),realmin);
    if c.family ~= "mda"
        if r == 0
            gradientInitial = [initial(2);0];
        else
            gradientInitial = [initial(2);r*initial(1)];
        end
    else
        gradientInitial = initial;
    end
    profile.nodes = c.joins;
    profile.segments = cell(1,numel(c.joins)-1);
    for i = 1:numel(c.joins)-1
        segment = denseGradientSegment(c,r,c.joins(i:i+1),gradientInitial);
        profile.segments{i} = segment;
        gradientInitial = evaluateGradientSegment(segment,segment.interval(2));
    end
    profile.kind = "shooting";
    profile.zeroF = initial(1);
end
if c.family == "mda"
    profile.surfaceGradient = evaluateGradientSegment(profile.segments{end},c.zDomain(2));
end
end

function segment = denseGradientSegment(c,r,interval,initial)
segment.interval = interval;
segment.initial = initial;
a = c.N2([interval(1);mean(interval);interval(2)])/c.g;
if all(a == a(1))
    segment.A = [0 1;-r*a(1) 0];
    segment.solution = [];
else
    segment.A = [];
    odeOptions = odeset("RelTol",c.tolerance,"AbsTol",c.tolerance*1e-2);
    segment.solution = ode113(@(z,y) [y(2);-r*c.N2(z)*y(1)/c.g],interval,initial,odeOptions);
end
end

function values = evaluateGradientSegment(segment,z)
if isempty(z)
    values = zeros(2,0);
    return
end
if isempty(segment.solution)
    values = zeros(2,numel(z));
    for i = 1:numel(z)
        values(:,i) = constantTransfer(segment.A,z(i)-segment.interval(1))*segment.initial;
    end
else
    values = deval(segment.solution,z);
end
end

function values = evaluateProfiles(c,profiles,z)
z = z(:);
nModes = numel(profiles);
values = struct("F",zeros(numel(z),nModes),"G",zeros(numel(z),nModes), ...
    "Fz",zeros(numel(z),nModes),"Gz",zeros(numel(z),nModes));
for j = 1:nModes
    profile = profiles{j};
    gradient = zeros(2,numel(z));
    for i = 1:numel(profile.segments)
        segment = profile.segments{i};
        selected = z >= segment.interval(1) & z <= segment.interval(2);
        gradient(:,selected) = evaluateGradientSegment(segment,z(selected));
    end
    values.G(:,j) = gradient(1,:).';
    values.Gz(:,j) = gradient(2,:).';
    values.Fz(:,j) = -(c.N2(z)/c.g).*values.G(:,j);
    if c.family == "mda"
        if profile.r == 0
            values.F(:,j) = neutralIntegral(c,values.G(:,j),z,profile);
        else
            values.F(:,j) = (values.Gz(:,j)-profile.surfaceGradient(2))/profile.r;
        end
    elseif profile.r == 0
        values.F(:,j) = profile.zeroF;
        for i = 1:numel(z)
            selectedJoins = c.joins(c.joins < z(i));
            intervals = [selectedJoins z(i)];
            for k = 1:numel(intervals)-1
                values.F(i,j) = values.F(i,j)-values.G(i,j)*integral(@(x) c.N2(x)/c.g, ...
                    intervals(k),intervals(k+1),"RelTol",c.tolerance,"AbsTol",c.tolerance*1e-3);
            end
        end
    else
        values.F(:,j) = values.Gz(:,j)/profile.r;
    end
end
end

function F = neutralIntegral(c,~,z,profile)
% The MDA zero mode is constant for finite endpoint accelerations.
amplitude = profile.surfaceGradient(1);
F = zeros(size(z));
for i = 1:numel(z)
    selectedJoins = c.joins(c.joins > z(i));
    intervals = [z(i) selectedJoins];
    for j = 1:numel(intervals)-1
        F(i) = F(i)+amplitude*integral(@(x) c.N2(x)/c.g,intervals(j),intervals(j+1), ...
            "RelTol",c.tolerance,"AbsTol",c.tolerance*1e-3);
    end
end
end

function factors = normalizationFactors(c,values,weights)
if c.family == "mda"
    normSquared = sum(weights.*(c.N2(weightedZ(c,weights))/c.g).*values.G.^2,1);
    if isfinite(c.g0)
        normSquared = normSquared+(c.g0/c.g)*values.G(end,:).^2;
    end
    if isfinite(c.gd)
        normSquared = normSquared+(c.gd/c.g)*values.G(1,:).^2;
    end
    factors = sqrt(abs(normSquared));
else
    normSquared = sum(weights.*values.F.^2,1);
    if c.family == "gpe"
        normSquared = normSquared+(c.f0^2/c.alpha0)*values.F(end,:).^2+(c.f0^2/c.alphaD)*values.F(1,:).^2;
    end
    factors = sqrt(normSquared/diff(c.zDomain));
end
end

function z = weightedZ(c,~)
[z,~] = neutralGeostrophicQuadrature(c.joins,513);
end

function values = normalizedProfiles(c,profiles,factors,z)
values = evaluateProfiles(c,profiles,z);
names = ["F","G","Fz","Gz"];
for name = names
    values.(name) = values.(name)./factors;
end
end

function eigenvalues = physicalEigenvalues(c,values,weights,roots)
eigenvalues = roots.';
if c.family == "gpe"
    z = weightedZ(c,weights);
    denominator = sum(weights.*values.F.^2,1)+(c.f0^2/c.alpha0)*values.F(end,:).^2+(c.f0^2/c.alphaD)*values.F(1,:).^2;
    numerator = c.k^2*sum(weights.*values.F.^2,1)+(c.f0^2/c.g^2)*sum(weights.*c.N2(z).*values.G.^2,1)+(c.f0^2/c.g)*values.F(end,:).^2;
    eigenvalues = numerator./denominator;
    eigenvalues(roots == 0) = c.k^2;
end
end

function labels = rootLabels(c,roots)
if c.family == "gpe"
    labels = 1:numel(roots);
else
    nNegative = nnz(roots < 0);
    nZero = nnz(roots == 0);
    labels = [-1:-1:-nNegative zeros(1,nZero) 1:(numel(roots)-nNegative-nZero)];
end
end

function reference = zeroAPVReference(c,zq,~)
r = -c.g*c.k^2/c.f0^2;
% Prescribe the two canonical endpoint responses directly in the local
% matching system. Its columns are independent of any production solve.
c.family = "mda";
c.g0 = Inf;
c.gd = Inf;
[d,l,u,nodes,transfers] = matchingSystem(c,r);
matrix = diag(d)+diag(l,-1)+diag(u,1);
% For zero APV, B_s=G_s-sigma*F_s with F_s=G'_s/r.
T = transfers{end};
surfaceGradientRow = zeros(1,numel(nodes));
surfaceGradientRow(end-1:end) = [-1/T(1,2) T(2,2)/T(1,2)];
matrix(end,:) = 0;
matrix(end,end) = 1;
matrix(end,:) = matrix(end,:)-(c.sigma/r)*surfaceGradientRow;
rhs = zeros(numel(nodes),2);
rhs(end,1) = 1;
rhs(1,2) = 1;
G = matrix\rhs;
profiles = cell(1,2);
for j = 1:2
    profile.r = r;
    profile.nodes = nodes;
    profile.segments = cell(1,numel(nodes)-1);
    for i = 1:numel(nodes)-1
        T = transfers{i};
        P = (G(i+1,j)-T(1,1)*G(i,j))/T(1,2);
        profile.segments{i} = denseGradientSegment(c,r,nodes(i:i+1),[G(i,j);P]);
    end
    profiles{j} = profile;
end
c.family = "zeroAPV";
reference.family = "zeroAPV";
reference.inverseDepth = [r r];
reference.h = [1/r 1/r];
reference.eigenvalues = [NaN NaN];
reference.modeNumber = [1 2];
reference.rootBrackets = zeros(0,2);
reference.tolerance = c.tolerance;
reference.joins = c.joins;
reference.evaluate = @(z) evaluateProfiles(c,profiles,z);
reference.configuration = c;
reference.scan = struct("endpointResponses",["surface","bottom"],"localIntervals",numel(nodes)-1);
% Force evaluation now so a malformed reference fails before readiness work.
reference.evaluate(zq);
end
