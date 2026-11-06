function results = HydrostaticRootTransformComparison(options)
% Compare the hydrostatic root transform with the v4 algorithm on identical points.
%
% For each profile and point budget, `IMInternalModesDiscreteTransform.hydrostaticRoot` builds
% the root grid, fitted weights and qualified band, with its default
% truncated-collocation operator and, separately, the Galerkin operator. The
% v4 algorithm is reproduced from the same converged modes on the same
% points: square collocation with F = [1, F_1..F_{Nz-1}] and the moment
% weights w4 = F^{-T}[D;0;...;0]. On identical points and modes the
% truncated-collocation rows equal v4's. Reported quantities:
%
% - usable band: the largest r with F and G energy errors <= tolerance under
%   the v4 moment weights and under the fitted weights (fitted per band);
% - omitted-mode gain: retained energy amplitude per unit energy amplitude
%   of internal modes r+1..r+8 under truncated collocation and Galerkin;
% - nonmodal input: truncated-collocation, Galerkin and prescribed-dual F
%   coefficients of u(z) = exp(z/300), relative to the continuous projection;
% - timing: transform construction, and repeated batched forward/back
%   transforms of an Nz-by-nBatch array (dense multiplies in both cases).
%
% ```matlab
% results = HydrostaticRootTransformComparison(nPoints=[17 33 65]);
% ```
%
% - Topic: Compare discrete transforms
% - Declaration: results = HydrostaticRootTransformComparison(options)
% - Parameter options.nPoints: point budgets
% - Parameter options.gramTolerance: energy-error tolerance
% - Parameter options.nBatch: columns in the timing batch
% - Parameter options.outputFile: optional CSV path
% - Parameter options.verbose: print progress for each case
% - Returns results: one table row per profile and point budget
arguments
    options.nPoints (1,:) double {mustBeInteger} = [17 33 65]
    options.gramTolerance (1,1) double = 1e-2
    options.nBatch (1,1) double {mustBeInteger,mustBePositive} = 128^2
    options.outputFile (1,1) string = ""
    options.verbose (1,1) logical = false
end
profiles = {
    "constant",       @(z) 2.5e-5*ones(size(z)),            1000
    "gmExponential",  @(z) (5.2e-3)^2*exp(2*z/1300),         4000
    "exponential400", @(z) 1e-4*exp(2*z/400),                1000
    "pycnocline",     @(z) 1e-6+4e-5*exp(-((z+300)/80).^2),  2000 };
rows = {};
for p = 1:size(profiles,1)
    [name,N2,D] = profiles{p,:};
    background = IMBackground(N2=N2,zDomain=[-D 0],latitude=33,surfaceCondition="rigidLid");
    for nPoints = options.nPoints
        if options.verbose, fprintf("%s, %d points: building transforms\n",name,nPoints); end
        t0 = tic;
        [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints,gramTolerance=options.gramTolerance);
        setupTime = toc(t0);
        galerkinTransform = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints,gramTolerance=options.gramTolerance,forwardOperator="galerkin");
        % Continuous reference modes at the converged resolution, beyond the retained band.
        solver = IMSolverSpectral(nEVP=2*assessment.convergence.candidateNEVP,coordinateKind=assessment.convergence.coordinateKind);
        nModes = nPoints+8;
        evp = IMInternalModes.hydrostaticGModes(background);
        basis = IMHydrostaticModesBasis(solver.solveEVP(evp,nModes=nModes),normalization="rmsF");
        z = grid.z; g = background.g; h = basis.h(2:end);
        F = basis.F(z); G = basis.G(z);
        HF = D*ones(1,nModes+1); HG = D./h; GMetric = N2(z)/g;
        % v4 algorithm on the same points.
        Fsquare = F(:,1:nPoints);
        w4 = Fsquare.'\[D;zeros(nPoints-1,1)];
        collocationForward = inv(Fsquare); % v4 forms the explicit inverse
        usableV4 = 0;
        for r = nPoints-2:-1:1
            [eF,eG] = energyErrors(F(:,1:r+1),G(:,2:r+1),w4,HF(1:r+1),HG(1:r),GMetric);
            if eF <= options.gramTolerance && eG <= options.gramTolerance, usableV4 = r; break, end
        end
        r = assessment.retainedInternalModeCount;
        [eF4,eG4] = energyErrors(F(:,1:r+1),G(:,2:r+1),w4,HF(1:r+1),HG(1:r),GMetric);
        omitted = r+2:min(r+9,nModes+1);
        collocationF = transform.forwardMatrix(variable="F");
        galerkinF = galerkinTransform.forwardMatrix(variable="F");
        gain = @(A) norm(sqrt(HF(1:r+1)).'.*(A*F(:,omitted))./sqrt(HF(omitted)),2);
        v4Rows = collocationForward(1:r+1,:);
        v4Difference = norm(collocationF-v4Rows)/norm(v4Rows);
        % Nonmodal input u = exp(z/300).
        fine = IMVerticalGrid.physicalDepth([-D 0],2049);
        Ffine = basis.F(fine.z);
        u = @(zz) exp(zz/300);
        continuous = (Ffine(:,1:r+1).'*(fine.weights.*u(fine.z)))./HF(1:r+1).';
        galerkin = galerkinF*u(z);
        dual = galerkinTransform.modeProjectionFunctional(u(z),variable="F")./HF(1:r+1).';
        collocation = collocationF*u(z);
        relative = @(c) norm(sqrt(HF(1:r+1)).'.*(c-continuous))/norm(sqrt(HF(1:r+1)).'.*continuous);
        % Timing of repeated batched transforms.
        X = randn(nPoints,options.nBatch);
        GForward = transform.forwardMatrix(variable="G"); FBack = transform.inverseMatrix(variable="F");
        forwardTime = timeit(@() collocationF*X);
        backTime = timeit(@() FBack*(GForward*X));
        v4Time = timeit(@() v4Rows*X);
        if options.verbose, fprintf("%s, %d points: done (r = %d)\n",name,nPoints,assessment.retainedInternalModeCount); end
        rows(end+1,:) = {name,nPoints,r,usableV4,assessment.final.energyErrorF,assessment.final.energyErrorG,eF4,eG4, ...
            gain(collocationF),gain(galerkinF),v4Difference,relative(collocation),relative(galerkin),relative(dual),min(grid.weights),min(w4), ...
            assessment.convergence.candidateNEVP,setupTime,forwardTime,backTime,v4Time}; %#ok<AGROW>
    end
end
results = cell2table(rows,VariableNames=["profile","nPoints","retainedNew","usableV4Weights","energyErrorFNew","energyErrorGNew", ...
    "energyErrorFV4WeightsSameBand","energyErrorGV4WeightsSameBand","omittedGainFCollocation","omittedGainFGalerkin", ...
    "collocationRelativeDifferenceFromV4","nonmodalErrorCollocation","nonmodalErrorGalerkin","nonmodalErrorPrescribedDual","minimumWeightNew","minimumWeightV4", ...
    "nEVP","setupSeconds","forwardSeconds","forwardBackSeconds","v4ForwardSeconds"]);
if options.outputFile ~= ""
    writetable(results,options.outputFile);
end
end

function [epsF,epsG] = energyErrors(F,G,w,HF,HG,GMetric)
gram = F.'*(w.*F);
epsF = norm((gram-diag(HF))./sqrt(HF(:)*HF(:).'),2);
gram = G.'*((GMetric.*w).*G);
epsG = norm((gram-diag(HG))./sqrt(HG(:)*HG(:).'),2);
end
