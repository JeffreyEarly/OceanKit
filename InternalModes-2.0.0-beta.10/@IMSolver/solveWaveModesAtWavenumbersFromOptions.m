function [collection,diagnostics] = solveWaveModesAtWavenumbersFromOptions(self,kappa,options)
% Solve continuous wave bases for an ordered collection of wavenumbers.
%
% Each distinct requested wavenumber is solved exactly once. Coordinates,
% stratification samples, derivative matrices and the shared pencil terms
% are prepared once for this call. The wave pencil is
% $$A(\kappa)=A_0+\kappa^2 A_2,\qquad B(\kappa)=B_0,$$
% with zero endpoint rows in $$A_2$$ because boundary equations replace the
% interior operator there. No approximate sharing or global cache is used.
% Shared preparation supports the exact built-in IMSolverSpectral and
% IMSolverFiniteDifference classes. Custom solvers must use independent
% solveEVP calls until they expose a validated preparation-reuse capability.
%
% Counts align with every requested page. Shared solves retain the largest
% count at each kappa; distinct requested prefixes preserve their identities.
% - Developer: true
arguments (Input)
    self (1,1) IMSolver
    kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.N2 (1,1) function_handle
    options.zDomain (1,2) double {mustBeReal,mustBeFinite}
    options.rotationRate (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 7.2921e-5
    options.latitude (1,1) double {mustBeReal,mustBeFinite,mustBeGreaterThanOrEqual(options.latitude,-90),mustBeLessThanOrEqual(options.latitude,90)} = 33
    options.g (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 9.81
    options.surfaceBoundary (1,1) IMBoundaryCondition
    options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
    options.nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    options.normalization = []
    options.background = []
end
arguments (Output)
    collection (1,1) IMBasisCollection
    diagnostics (1,1) struct
end
if ~ismember(string(class(self)),["IMSolverSpectral","IMSolverFiniteDifference"])
    error("IMSolver:UnsupportedBulkSolver","Shared wave preparation supports built-in IMSolverSpectral and IMSolverFiniteDifference only. Use independent solveEVP calls for custom solver classes.");
end
plan = IMModeSelection.wavePagePlan(kappa,options.nModes);
if any(~isfinite(kappa.^2))
    error("IMSolver:InvalidSquaredWavenumber","Every squared wavenumber must remain finite.");
end
totalTimer = tic;
distinctKappa = plan.solvedKappa;
basisIndex = plan.basisIndex;
factoryOptions = rmfield(options,["nModes","normalization","background"]);
factoryArguments = namedargs2cell(factoryOptions);
if ~isempty(options.background), factoryArguments = [{options.background} factoryArguments]; end
setupTimer = tic;
baseEVP = IMInternalModes.waveModesAtWavenumber(factoryArguments{:},k=0);
solver = self.configuredForEVP(baseEVP);
setupSeconds = toc(setupTimer);
assemblyTimer = tic;
[A0,B,samples] = baseEVP.assembleConfigured(solver);
A2 = solver.physicalDerivativeMatrix(0);
A2(solver.boundaryIndex("surface"),:) = 0;
A2(solver.boundaryIndex("bottom"),:) = 0;
sharedAssemblySeconds = toc(assemblyTimer);
assemblyStorage = whos("A0","A2","B","samples");
bases = cell(1,numel(distinctKappa));
requestedCount = plan.solveCounts;
pageAssemblySeconds = zeros(size(distinctKappa));
eigensolveSeconds = zeros(size(distinctKappa));
finalizationSeconds = zeros(size(distinctKappa));
for iK = 1:numel(distinctKappa)
    assemblyTimer = tic;
    evp = IMInternalModes.waveModesAtWavenumber(factoryArguments{:},k=distinctKappa(iK));
    A = A0+distinctKappa(iK)^2*A2;
    pageSamples = samples;
    pageSamples.q = distinctKappa(iK)^2*ones(size(samples.q));
    pageAssemblySeconds(iK) = toc(assemblyTimer);
    finalizationTimer = tic;
    selectionDiagnostics = evp.preparedModeSelectionDiagnostics(pageSamples,A);
    diagnosticSeconds = toc(finalizationTimer);
    request = struct(counted=true,nModes=requestedCount(iK),branch="");
    [basis,cost] = solver.solveAssembledEVP(evp,A,B,max(1,requestedCount(iK)),selectionDiagnostics,request);
    if ~isempty(options.normalization), basis.normalization = options.normalization; end
    if numel(basis.modeNumber) ~= requestedCount(iK)
        error("IMSolver:InsufficientWaveModes","At kappa=%g the solve returned %d columns, but %d were explicitly requested. Increase solver resolution or revise the requested count.",distinctKappa(iK),numel(basis.modeNumber),requestedCount(iK));
    end
    bases{iK} = basis.withConstructionRequest(struct(nModes=requestedCount(iK),branch=""),self);
    eigensolveSeconds(iK) = cost.eigensolveSeconds;
    finalizationSeconds(iK) = diagnosticSeconds+cost.finalizationSeconds;
end
collectionTimer = tic;
collection = IMModeSelection.wavePrefixCollection(bases,plan,self);
collectionSeconds = toc(collectionTimer);
diagnostics = struct("requestedKappa",kappa,"solvedKappa",distinctKappa,"basisIndex",basisIndex,"requestedPageCount",numel(kappa), ...
    "solveCount",numel(distinctKappa),"reusedPageCount",numel(kappa)-numel(distinctKappa),"preparationCount",1,"reuseKind","exactEqualKappa", ...
    "solverClass",string(class(solver)),"nModes",options.nModes,"requestedCounts",plan.pageCounts,"solveAccuracy","unverified", ...
    "setupSeconds",setupSeconds,"sharedAssemblySeconds",sharedAssemblySeconds,"assemblySeconds",sharedAssemblySeconds+sum(pageAssemblySeconds), ...
    "eigensolveSeconds",sum(eigensolveSeconds),"finalizationSeconds",sum(finalizationSeconds),"collectionSeconds",collectionSeconds,"totalSeconds",0, ...
    "sharedAssemblyStorageBytes",sum([assemblyStorage.bytes]));
diagnostics.solves = table(distinctKappa(:),requestedCount(:),pageAssemblySeconds(:),eigensolveSeconds(:),finalizationSeconds(:), ...
    VariableNames=["kappa","requestedColumnCount","assemblySeconds","eigensolveSeconds","finalizationSeconds"]);
diagnostics.totalSeconds = toc(totalTimer);
end
