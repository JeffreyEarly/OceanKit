function [collection,diagnostics] = solveWaveModesAtWavenumbers(self,kappa,options)
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
% Positive wavenumbers request `nModes` columns. A zero-wavenumber request
% requires an explicit independent `nInertialModes` count. Both counts are
% exact requests: an insufficient returned family raises an error. Stored
% bases retain the same scientific labels, normalization and boundary
% conditions as independent `solveEVP` calls. Frequency signs and physical
% wave polarizations belong to the caller, not this vertical-basis solve.
%
% ```matlab
% [bases,costs] = solver.solveWaveModesAtWavenumbers([2e-4 0 1e-4 2e-4],N2=N2,zDomain=[-1000 0],nModes=8,nInertialModes=5);
% G = bases.evaluate(z,variable="G",pages=[1 3]);
% ```
%
% - Topic: Solve EVPs
% - Declaration: [collection,diagnostics] = solveWaveModesAtWavenumbers(solver,kappa,options)
% - Parameter kappa: nonempty row of requested nonnegative wavenumbers, in radians per meter
% - Parameter options.N2: buoyancy frequency squared evaluator
% - Parameter options.zDomain: increasing physical vertical bounds
% - Parameter options.f0: Coriolis parameter
% - Parameter options.g: gravitational acceleration
% - Parameter options.surfaceBoundary: surface condition of the G EVP
% - Parameter options.bottomBoundary: bottom condition of the G EVP
% - Parameter options.nModes: uniform count for positive wavenumbers
% - Parameter options.nInertialModes: independent count, required when zero is requested
% - Returns collection: continuous bases and exact requested-page mapping
% - Returns diagnostics: construction costs and solve provenance, not accuracy qualification
arguments (Input)
    self (1,1) IMSolver
    kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.N2 (1,1) function_handle
    options.zDomain (1,2) double {mustBeReal,mustBeFinite}
    options.f0 (1,1) double {mustBeReal,mustBeFinite} = 0
    options.g (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 9.81
    options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
    options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 100
    options.nInertialModes double {mustBeInteger,mustBePositive} = []
end
arguments (Output)
    collection (1,1) IMBasisCollection
    diagnostics (1,1) struct
end
if ~ismember(string(class(self)),["IMSolverSpectral","IMSolverFiniteDifference"])
    error("IMSolver:UnsupportedBulkSolver","Shared wave preparation supports built-in IMSolverSpectral and IMSolverFiniteDifference only. Use independent solveEVP calls for custom solver classes.");
end
if isempty(kappa)
    error("IMSolver:EmptyWavenumberCollection","Supply at least one requested wavenumber.");
end
if ~isempty(options.nInertialModes) && ~isscalar(options.nInertialModes)
    error("IMSolver:InvalidInertialModeCount","nInertialModes must be a scalar count.");
end
if any(kappa == 0) && isempty(options.nInertialModes)
    error("IMSolver:MissingInertialModeCount","A zero-wavenumber page requires an explicit nInertialModes count independent of nModes.");
end
if any(~isfinite(kappa.^2))
    error("IMSolver:InvalidSquaredWavenumber","Every squared wavenumber must remain finite.");
end
totalTimer = tic;
[distinctKappa,~,basisIndex] = unique(kappa,"stable");
basisIndex = reshape(basisIndex,1,[]);
factoryOptions = rmfield(options,["nModes","nInertialModes"]);
factoryArguments = namedargs2cell(factoryOptions);
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
requestedCount = repmat(options.nModes,size(distinctKappa));
if any(distinctKappa == 0)
    requestedCount(distinctKappa == 0) = options.nInertialModes;
end
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
    [basis,cost] = solver.solveAssembledEVP(evp,A,B,requestedCount(iK),selectionDiagnostics);
    if numel(basis.modeNumber) ~= requestedCount(iK)
        error("IMSolver:InsufficientWaveModes","At kappa=%g the solve returned %d columns, but %d were explicitly requested. Increase solver resolution or revise the requested count.",distinctKappa(iK),numel(basis.modeNumber),requestedCount(iK));
    end
    bases{iK} = basis;
    eigensolveSeconds(iK) = cost.eigensolveSeconds;
    finalizationSeconds(iK) = diagnosticSeconds+cost.finalizationSeconds;
end
collectionTimer = tic;
collection = IMBasisCollection(bases,kappa=kappa,basisIndex=basisIndex);
collectionSeconds = toc(collectionTimer);
diagnostics = struct("requestedKappa",kappa,"solvedKappa",distinctKappa,"basisIndex",basisIndex,"requestedPageCount",numel(kappa), ...
    "solveCount",numel(distinctKappa),"reusedPageCount",numel(kappa)-numel(distinctKappa),"preparationCount",1,"reuseKind","exactEqualKappa", ...
    "solverClass",string(class(solver)),"nModes",options.nModes,"nInertialModes",options.nInertialModes,"solveAccuracy","unverified", ...
    "setupSeconds",setupSeconds,"sharedAssemblySeconds",sharedAssemblySeconds,"assemblySeconds",sharedAssemblySeconds+sum(pageAssemblySeconds), ...
    "eigensolveSeconds",sum(eigensolveSeconds),"finalizationSeconds",sum(finalizationSeconds),"collectionSeconds",collectionSeconds,"totalSeconds",0, ...
    "sharedAssemblyStorageBytes",sum([assemblyStorage.bytes]));
diagnostics.solves = table(distinctKappa(:),requestedCount(:),pageAssemblySeconds(:),eigensolveSeconds(:),finalizationSeconds(:), ...
    VariableNames=["kappa","requestedColumnCount","assemblySeconds","eigensolveSeconds","finalizationSeconds"]);
diagnostics.totalSeconds = toc(totalTimer);
end
