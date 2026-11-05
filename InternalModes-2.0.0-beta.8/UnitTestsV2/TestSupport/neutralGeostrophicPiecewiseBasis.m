function basis = neutralGeostrophicPiecewiseBasis(problem,nEVP,options)
% Solve the neutral-cap trial independently on joined Chebyshev intervals.
% This study-only solver preserves the total point budget and imposes
% continuity of the regular pair (or G/G' for the scalar MDA problem).
arguments
    problem
    nEVP (1,1) double {mustBeInteger,mustBePositive}
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 40
    options.joins (1,:) double {mustBeReal,mustBeFinite} = [-1000 -400 -200 0]
end
c = describeProblem(problem,options.joins);
segments = makeSegments(c,nEVP);
N = sum([segments.n]);
[zq,wq] = neutralGeostrophicQuadrature(c.joins,513);
c.integratedN2 = sum(wq.*c.N2(zq));
if c.integratedN2<=0
    error("NeutralGeostrophicPiecewise:NeutralDomain","This trial requires a nonzero stratified interval.");
end
rScale = c.g/(c.integratedN2*diff(c.zDomain));
if c.family=="mda"
    [A,B] = assembleMDA(c,segments,N);
    fields = ones(N,1);
else
    [A,B] = assemblePair(c,segments,N);
    fields = [c.integratedN2/c.g*ones(N,1);ones(N,1)];
end
if c.family=="zeroAPV"
    rhs = zeros(2*N,2);
    rhs(segments(end).index(1),1) = 1;
    rhs(N+segments(1).index(end),2) = 1;
    r = -c.g*c.k^2/c.f0^2;
    matrix = A-r*B;
    columns = fields;
    rows = max(abs(matrix.*columns.'),[],2);
    rows(rows==0) = 1;
    factorization = decomposition((matrix.*columns.')./rows);
    V = columns.*(factorization\(rhs./rows));
    inverseDepth = [r r];
    eigenvalues = inverseDepth;
    labels = [1 2];
    qualification = true(1,2);
else
    Bs = rScale*B;
    rows = max(max(abs(A.*fields.'),[],2),max(abs(Bs.*fields.'),[],2));
    rows(rows==0) = 1;
    Ar = A./rows;
    Br = Bs./rows;
    degreeScale = zeros(size(fields));
    for i = 1:numel(segments)
        index = segments(i).index;
        degreeScale(index) = max(1,0:segments(i).n-1).'.^(-2);
        if c.family~="mda", degreeScale(N+index) = degreeScale(index); end
    end
    columns = fields.*degreeScale;
    [V,L] = eig(Ar.*columns.',Br.*columns.');
    V = columns.*V;
    lambda = diag(L);
    finite = isfinite(real(lambda)) & isfinite(imag(lambda));
    V = V(:,finite);
    lambda = lambda(finite);
    realPair = abs(imag(lambda))<=1e-10*max(abs(real(lambda)),1);
    for j = 1:numel(lambda)
        [~,k] = max(abs(V(:,j)));
        V(:,j) = V(:,j)*exp(-1i*angle(V(k,j)));
    end
    realPair = realPair & vecnorm(imag(V),2,1).'<=1e-8*vecnorm(real(V),2,1).';
    V = real(V(:,realPair));
    lambda = real(lambda(realPair));
    balanced = V./fields;
    action = vecnorm(Br*V,2,1)./max(norm(Br.*fields.',1)*vecnorm(balanced,2,1),realmin);
    observable = action>1e3*eps;
    if c.family~="mda"
        observable = observable & vecnorm(V(1:N,:),2,1)/(c.integratedN2/c.g)>1e3*eps*vecnorm(balanced,2,1);
    end
    V = V(:,observable);
    lambda = lambda(observable);
    [negativeCount,hasZero] = modeInertia(c,problem);
    if hasZero && ~isempty(lambda)
        [~,j] = min(abs(lambda));
        zeroResidual = norm(Ar*V(:,j))/max(norm(Ar.*fields.',1)*norm(V(:,j)./fields),realmin);
        if abs(lambda(j))<1e-7 && zeroResidual<1e-9
            lambda(j) = 0;
            if c.family=="mda"
                V(:,j) = 0;
                for i = 1:numel(segments), V(segments(i).index(1),j) = 1; end
            elseif c.family=="gpe" && isinf(c.alphaD) ...
                    && (c.alpha0==c.g*c.k^2 || c.alpha0==(c.g*c.k)*c.k)
                V(:,j) = 0;
                for i = 1:numel(segments), V(segments(i).index(1),j) = 1; end
            end
        end
    end
    inverseDepth = (rScale*lambda).';
    eigenvalues = inverseDepth;
    if c.family=="gpe"
        rawSegments = storeModes(c,segments,V,N);
        values = evaluateModes(c,rawSegments,zq);
        volumeF = sum(wq.*values.F.^2,1);
        numerator = c.k^2*volumeF+c.f0^2/c.g^2*sum(wq.*c.N2(zq).*values.G.^2,1)+c.f0^2/c.g*values.F(end,:).^2;
        denominator = volumeF+c.f0^2/c.alpha0*values.F(end,:).^2+c.f0^2/c.alphaD*values.F(1,:).^2;
        eigenvalues = numerator./denominator;
        eigenvalues(inverseDepth==0) = c.k^2;
    end
    residual = vecnorm(Ar*V-(Br*V).*lambda.',2,1) ...
        ./max((norm(Ar.*fields.',1)+abs(lambda.')*norm(Br.*fields.',1)).*vecnorm(V./fields,2,1),realmin);
    qualification = residual<=1e-8;
    if c.family=="gpe"
        affine = c.k^2+c.f0^2/c.g*inverseDepth;
        scale = c.k^2+c.f0^2/c.g*abs(inverseDepth);
        qualification = qualification & eigenvalues>0 & isfinite(eigenvalues) & abs(eigenvalues-affine)<=1e-7*scale;
        eigenvalues(~qualification) = max(affine(~qualification),realmin);
    end
    [index,labels] = selectModes(c,eigenvalues,negativeCount,hasZero,options.nModes);
    if any(~qualification(index))
        error("NeutralGeostrophicPiecewise:UnresolvedMode","A requested piecewise eigenpair failed pencil qualification.");
    end
    V = V(:,index);
    inverseDepth = inverseDepth(index);
    eigenvalues = eigenvalues(index);
    qualification = qualification(index);
end
cleanup = struct("algorithm","none","relativeTolerance",eps,"cutoffDegrees",[], ...
    "relativeCoefficientPerturbation",zeros(1,size(V,2)),"postCleanupPencilResidual",[]);
if c.family=="mda"
    original = V;
    cutoffs = zeros(numel(segments),size(V,2));
    for i = 1:numel(segments)
        [V(segments(i).index,:),cutoffs(i,:)] = chopCoefficients(V(segments(i).index,:));
    end
    perturbation = vecnorm(V-original,2,1)./max(vecnorm(original,2,1),realmin);
    normalizedR = inverseDepth/rScale;
    postResidual = vecnorm(Ar*V-(Br*V).*normalizedR,2,1) ...
        ./max((norm(Ar,1)+abs(normalizedR)*norm(Br,1)).*vecnorm(V,2,1),realmin);
    if any(postResidual>1e-8)
        error("NeutralGeostrophicPiecewise:CleanupResidual","Machine-precision MDA coefficient cleanup failed the original pencil residual check.");
    end
    cleanup = struct("algorithm","Chebfun standardChop","relativeTolerance",eps, ...
        "cutoffDegrees",cutoffs,"relativeCoefficientPerturbation",perturbation,"postCleanupPencilResidual",postResidual);
end
stored = storeModes(c,segments,V,N);
raw = evaluateModes(c,stored,zq);
if c.family=="mda"
    signedNorm = sum(wq.*(c.N2(zq)/c.g).*raw.G.^2,1);
    if isfinite(c.g0), signedNorm = signedNorm+c.g0/c.g*raw.G(end,:).^2; end
    if isfinite(c.gd), signedNorm = signedNorm+c.gd/c.g*raw.G(1,:).^2; end
    majorant = sum(wq.*(c.N2(zq)/c.g).*raw.G.^2,1);
    if isfinite(c.g0), majorant = majorant+abs(c.g0)/c.g*raw.G(end,:).^2; end
    if isfinite(c.gd), majorant = majorant+abs(c.gd)/c.g*raw.G(1,:).^2; end
    if any(majorant<=0 | abs(signedNorm)./majorant<=1e3*eps)
        error("NeutralGeostrophicPiecewise:NullNorm","A retained MDA mode has zero signed norm relative to its majorant.");
    end
    factors = sqrt(abs(signedNorm));
    signatures = sign(signedNorm);
elseif c.family=="zeroAPV"
    factors = ones(1,2);
    signatures = [];
else
    normSquared = sum(wq.*raw.F.^2,1);
    if c.family=="gpe"
        normSquared = normSquared+c.f0^2/c.alpha0*raw.F(end,:).^2+c.f0^2/c.alphaD*raw.F(1,:).^2;
    end
    factors = sqrt(normSquared/diff(c.zDomain));
    signatures = [];
end
if any(~isfinite(factors) | factors<=0)
    error("NeutralGeostrophicPiecewise:NullNorm","The retained piecewise basis contains a null normalization.");
end
% Signs do not affect the pair or its interface equations. Comparison with
% shooting matches signs by the physical F/G overlap independently.
if c.family~="zeroAPV"
    [~,where] = max(abs(raw.G),[],1);
    signs = ones(size(factors));
    for j = 1:numel(signs)
        signs(j) = sign(raw.G(where(j),j));
        if signs(j)==0
            [~,k] = max(abs(raw.F(:,j)));
            signs(j) = sign(raw.F(k,j));
        end
    end
    factors = factors.*signs;
end
stored = scaleModes(stored,factors);
interface = interfaceChecks(c,stored);
basis = struct("family",c.family,"configuration",c,"h",1./inverseDepth, ...
    "inverseDepth",inverseDepth,"eigenvalues",eigenvalues,"modeNumber",labels, ...
    "signatures",signatures,"normalization","rmsF", ...
    "evaluate",@(z) evaluateModes(c,stored,z), ...
    "metadata",struct("discretization","piecewiseChebyshev","nEVP",nEVP, ...
    "nPerInterval",[segments.n],"pointsPerInterval",[segments.n],"joins",c.joins,"uniquePhysicalPointCount",N-numel(segments)+1, ...
    "interfaceResiduals",interface.maximumResidual,"interfaceJumps",{interface.jumps}, ...
    "qualifiedPencilModes",qualification,"referenceRootsUsed",false,"coefficientCleanup",cleanup));
if c.family=="mda", basis.normalization = "meanDensityAnomalyEnergy"; end
end

function c = describeProblem(problem,joins)
c = struct("N2",problem.N2,"zDomain",problem.zDomain,"g",problem.g,"f0",problem.f0, ...
    "sigma",1,"k",0,"joins",unique([problem.zDomain joins]));
c.joins = c.joins(c.joins>=c.zDomain(1) & c.joins<=c.zDomain(2));
if isa(problem,"IMGeostrophicZeroAPVModes")
    if numel(problem.k)~=1, error("NeutralGeostrophicPiecewise:ScalarWavenumber","Supply one zero-APV wavenumber."); end
    c.family = "zeroAPV";
    c.k = problem.k;
    c.sigma = double(problem.surfaceBoundary=="freeSurface");
elseif isa(problem,"IMMeanDensityAnomalyModes")
    c.family = "mda";
    c.g0 = problem.g0;
    c.gd = problem.gd;
elseif problem.name=="geostrophicAPVModes"
    c.family = "apv";
    c.g0 = problem.parameters.g0;
    c.gd = problem.parameters.gd;
    c.sigma = double(problem.parameters.surfaceBoundary=="freeSurface");
elseif problem.name=="geostrophicGeneralizedPotentialEnstrophyModes"
    c.family = "gpe";
    c.k = problem.parameters.k;
    c.alpha0 = problem.parameters.alpha0;
    c.alphaD = problem.parameters.alphaD;
else
    error("NeutralGeostrophicPiecewise:UnsupportedFamily","Unsupported balanced family.");
end
end

function segments = makeSegments(c,budget)
lengths = diff(c.joins);
m = numel(lengths);
if budget<16*m, error("NeutralGeostrophicPiecewise:InsufficientBudget","Each interval requires at least 16 points."); end
ideal = budget*lengths/sum(lengths);
counts = max(16,floor(ideal));
while sum(counts)<budget
    [~,j] = max(ideal-counts);
    counts(j) = counts(j)+1;
end
while sum(counts)>budget
    excess = counts-ideal;
    excess(counts<=16) = -Inf;
    [~,j] = max(excess);
    counts(j) = counts(j)-1;
end
segments = struct("solver",{},"interval",{},"n",{},"index",{},"T",{},"D",{},"D2",{});
offset = 0;
for i = 1:m
    interval = c.joins(i:i+1);
    probe = IMEigenvalueProblem(zDomain=interval,p=1,q=0,r=1);
    solver = IMSolverSpectral(nEVP=counts(i),coordinateKind="z").configuredForEVP(probe);
    segments(i) = struct("solver",solver,"interval",interval,"n",counts(i), ...
        "index",offset+(1:counts(i)),"T",solver.physicalDerivativeMatrix(0), ...
        "D",solver.physicalDerivativeMatrix(1),"D2",solver.physicalDerivativeMatrix(2));
    offset = offset+counts(i);
end
end

function [A,B] = assemblePair(c,segments,N)
A = zeros(2*N);
B = zeros(2*N);
for i = 1:numel(segments)
    s = segments(i);
    index = s.index;
    N2 = c.N2(s.solver.zNative);
    if ~isreal(N2) || any(~isfinite(N2)) || any(N2<0), error("NeutralGeostrophicPiecewise:InvalidN2","Invalid sampled N2."); end
    A(index,index) = s.D;
    A(index,N+index) = diag(N2/c.g)*s.T;
    A(N+index,N+index) = s.D;
    B(N+index,index) = s.T;
end
for i = 1:numel(segments)-1
    lower = segments(i);
    upper = segments(i+1);
    rowF = lower.index(1);
    rowG = N+upper.index(end);
    A([rowF rowG],:) = 0;
    B([rowF rowG],:) = 0;
    A(rowF,lower.index) = lower.T(1,:);
    A(rowF,upper.index) = -upper.T(end,:);
    A(rowG,N+lower.index) = lower.T(1,:);
    A(rowG,N+upper.index) = -upper.T(end,:);
end
surface = segments(end);
bottom = segments(1);
rs = surface.index(1);
rb = N+bottom.index(end);
A([rs rb],:) = 0;
B([rs rb],:) = 0;
if c.family=="apv"
    if isinf(c.g0)
        A(rs,surface.index) = -c.sigma*surface.T(1,:);
        A(rs,N+surface.index) = surface.T(1,:);
    else
        A(rs,surface.index) = -(c.g+c.sigma*c.g0)*surface.T(1,:);
        A(rs,N+surface.index) = c.g0*surface.T(1,:);
    end
    if isinf(c.gd)
        A(rb,N+bottom.index) = bottom.T(end,:);
    else
        A(rb,bottom.index) = c.g*bottom.T(end,:);
        A(rb,N+bottom.index) = c.gd*bottom.T(end,:);
    end
elseif c.family=="gpe"
    A(rs,surface.index) = (c.g*c.k^2/c.alpha0-1)*surface.T(1,:);
    A(rs,N+surface.index) = surface.T(1,:);
    B(rs,surface.index) = -c.f0^2/c.alpha0*surface.T(1,:);
    A(rb,bottom.index) = -c.g*c.k^2/c.alphaD*bottom.T(end,:);
    A(rb,N+bottom.index) = bottom.T(end,:);
    B(rb,bottom.index) = c.f0^2/c.alphaD*bottom.T(end,:);
else
    A(rs,surface.index) = -c.sigma*surface.T(1,:);
    A(rs,N+surface.index) = surface.T(1,:);
    A(rb,N+bottom.index) = bottom.T(end,:);
end
end

function [A,B] = assembleMDA(c,segments,N)
A = zeros(N);
B = zeros(N);
for i = 1:numel(segments)
    s = segments(i);
    N2 = c.N2(s.solver.zNative);
    A(s.index,s.index) = -s.D2;
    B(s.index,s.index) = diag(N2/c.g)*s.T;
end
for i = 1:numel(segments)-1
    lower = segments(i);
    upper = segments(i+1);
    rows = [lower.index(1),upper.index(end)];
    A(rows,:) = 0;
    B(rows,:) = 0;
    A(rows(1),lower.index) = lower.T(1,:);
    A(rows(1),upper.index) = -upper.T(end,:);
    A(rows(2),lower.index) = lower.D(1,:);
    A(rows(2),upper.index) = -upper.D(end,:);
end
surface = segments(end);
bottom = segments(1);
rs = surface.index(1);
rb = bottom.index(end);
A([rs rb],:) = 0;
B([rs rb],:) = 0;
if isinf(c.g0)
    A(rs,surface.index) = surface.T(1,:);
else
    A(rs,surface.index) = surface.D(1,:);
    B(rs,surface.index) = (c.g0/c.g)*surface.T(1,:);
end
if isinf(c.gd)
    A(rb,bottom.index) = bottom.T(end,:);
else
    A(rb,bottom.index) = bottom.D(end,:);
    B(rb,bottom.index) = -(c.gd/c.g)*bottom.T(end,:);
end
end

function [negative,zero] = modeInertia(c,problem)
if c.family=="mda"
    negative = double(c.g0<0)+double(c.gd<0);
    zero = isfinite(c.g0) && isfinite(c.gd);
    if zero
        total = c.integratedN2+c.g0+c.gd;
        if abs(total)<=1e3*eps*(c.integratedN2+abs(c.g0)+abs(c.gd))
            error("NeutralGeostrophicPiecewise:NullConstant","The MDA constant has zero signed norm.");
        end
        negative = negative-double(total<0);
    end
else
    probe = IMSolverSpectral(nEVP=16).configuredForEVP(problem);
    diagnostics = problem.modeSelectionDiagnostics(probe,[]);
    negative = diagnostics.minNegativeEigenvalueCount;
    zero = diagnostics.zeroInverseDepthCount>0;
end
end

function [index,labels] = selectModes(c,lambda,negativeCount,hasZero,nModes)
if c.family=="gpe"
    [~,index] = sort(lambda);
    labels = 1:min(nModes,numel(index));
    index = index(1:numel(labels));
    return
end
negative = find(lambda<0);
[~,order] = sort(lambda(negative));
if numel(negative)<min(nModes,negativeCount)
    error("NeutralGeostrophicPiecewise:MissingSector","A required negative sector was not found.");
end
negative = negative(order(1:min(negativeCount,numel(order))));
zero = find(lambda==0,1);
if ~hasZero, zero = []; elseif isempty(zero) && nModes>negativeCount
    error("NeutralGeostrophicPiecewise:MissingZero","The required zero mode was not found.");
end
positive = find(lambda>0);
[~,order] = sort(lambda(positive));
positive = positive(order);
index = [negative(:);zero(:);positive(:)];
labels = [-1:-1:-numel(negative),zeros(1,numel(zero)),1:numel(positive)];
count = min(nModes,numel(index));
index = index(1:count);
labels = labels(1:count);
end

function stored = storeModes(c,segments,V,N)
stored = segments;
for i = 1:numel(segments)
    index = segments(i).index;
    if c.family=="mda"
        stored(i).G = V(index,:);
        nSample = 2*segments(i).n;
        theta = (0:nSample-1)'*pi/(nSample-1);
        z = mean(segments(i).interval)+diff(segments(i).interval)*cos(theta)/2;
        G = segments(i).solver.evaluateNativeModes(stored(i).G,z);
        integrand = (c.N2(z)/c.g).*G;
        coefficients = chopCoefficients(lobattoCoefficients(integrand));
        primitive = zeros(size(coefficients,1)+1,size(coefficients,2));
        for j = 1:size(coefficients,2)
            polynomial = (diff(segments(i).interval)/2)*InternalModesSpectral.IntegrateChebyshevVector(coefficients(:,j));
            polynomial = -polynomial;
            polynomial(1) = polynomial(1)-sum(polynomial);
            primitive(:,j) = polynomial;
        end
        stored(i).F = primitive;
    else
        stored(i).F = V(index,:);
        stored(i).G = V(N+index,:);
    end
end
if c.family=="mda"
    accumulated = zeros(1,size(V,2));
    for i = numel(stored):-1:1
        stored(i).F(1,:) = stored(i).F(1,:)+accumulated;
        accumulated = stored(i).solver.evaluateNativeModes(stored(i).F,stored(i).interval(1));
    end
end
end

function [coefficients,degrees] = chopCoefficients(coefficients)
% The same machine-precision representation cleanup used by Chebfun.
% Padding retains the allocated pencil size; no N2 value is modified.
degrees = zeros(1,size(coefficients,2));
for j = 1:size(coefficients,2)
    cutoff = standardChop(coefficients(:,j),eps);
    coefficients(cutoff+1:end,j) = 0;
    degrees(j) = cutoff-1;
end
end

function coefficients = lobattoCoefficients(values)
n = size(values,1);
extended = [values;values(end-1:-1:2,:)];
coefficients = real(fft(extended,[],1))/(n-1);
coefficients = coefficients(1:n,:);
coefficients([1 end],:) = coefficients([1 end],:)/2;
end

function values = evaluateModes(c,stored,z)
z = z(:);
nModes = size(stored(1).G,2);
values = struct("F",zeros(numel(z),nModes),"G",zeros(numel(z),nModes), ...
    "Fz",zeros(numel(z),nModes),"Gz",zeros(numel(z),nModes),"Gzz",zeros(numel(z),nModes));
for i = 1:numel(stored)
    if i==numel(stored)
        selected = z>=stored(i).interval(1) & z<=stored(i).interval(2);
    else
        selected = z>=stored(i).interval(1) & z<stored(i).interval(2);
    end
    solver = stored(i).solver;
    points = z(selected);
    values.F(selected,:) = solver.evaluateNativeModes(stored(i).F,points);
    values.G(selected,:) = solver.evaluateNativeModes(stored(i).G,points);
    values.Fz(selected,:) = solver.evaluatePhysicalDerivative(stored(i).F,points,1);
    values.Gz(selected,:) = solver.evaluatePhysicalDerivative(stored(i).G,points,1);
    values.Gzz(selected,:) = solver.evaluatePhysicalDerivative(stored(i).G,points,2);
end
if any(z<c.zDomain(1) | z>c.zDomain(2))
    error("NeutralGeostrophicPiecewise:OutsideDomain","Evaluation depths must lie in the physical domain.");
end
end

function stored = scaleModes(stored,factors)
for i = 1:numel(stored)
    stored(i).F = stored(i).F./factors;
    stored(i).G = stored(i).G./factors;
end
end

function result = interfaceChecks(c,stored)
jumps = cell(1,numel(stored)-1);
maximum = 0;
FScale = zeros(1,size(stored(1).F,2));
GScale = FScale;
GzScale = FScale;
for i = 1:numel(stored)
    solver = stored(i).solver;
    FScale = max(FScale,max(abs(solver.evaluateNativeModes(stored(i).F,solver.zNative)),[],1));
    GScale = max(GScale,max(abs(solver.evaluateNativeModes(stored(i).G,solver.zNative)),[],1));
    GzScale = max(GzScale,diff(c.zDomain)*max(abs(solver.evaluatePhysicalDerivative(stored(i).G,solver.zNative,1)),[],1));
end
pairScale = max(FScale,GScale);
FScale(FScale<=1e3*eps*pairScale) = pairScale(FScale<=1e3*eps*pairScale);
GScale(GScale<=1e3*eps*pairScale) = pairScale(GScale<=1e3*eps*pairScale);
for i = 1:numel(stored)-1
    z = stored(i).interval(2);
    lowerF = stored(i).solver.evaluateNativeModes(stored(i).F,z);
    upperF = stored(i+1).solver.evaluateNativeModes(stored(i+1).F,z);
    lowerG = stored(i).solver.evaluateNativeModes(stored(i).G,z);
    upperG = stored(i+1).solver.evaluateNativeModes(stored(i+1).G,z);
    if c.family=="mda"
        lower = [lowerF;lowerG;diff(c.zDomain)*stored(i).solver.evaluatePhysicalDerivative(stored(i).G,z,1)];
        upper = [upperF;upperG;diff(c.zDomain)*stored(i+1).solver.evaluatePhysicalDerivative(stored(i+1).G,z,1)];
        scale = [FScale;GScale;max(GScale,GzScale)];
    else
        lower = [lowerF;lowerG];
        upper = [upperF;upperG];
        scale = [FScale;GScale];
    end
    jumps{i} = lower-upper;
    maximum = max(maximum,max(abs(jumps{i})./max(scale,realmin),[],"all"));
end
result = struct("jumps",{jumps},"maximumResidual",maximum);
end
