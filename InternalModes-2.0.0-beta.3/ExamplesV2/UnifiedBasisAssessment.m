% Evaluate and assess mode families through one explicit fixed-grid contract.
% See Documentation/AssessmentAPI.md for capabilities and prepared recipes.
% These examples measure sampled transforms; they do not qualify EVP accuracy.

zDomain = [-1000 0];
N0 = 5e-3;
N2 = @(z) N0^2*ones(size(z));
solver = IMSolverSpectral(nEVP=64);
z = linspace(zDomain(1),zDomain(2),65).';
weights = [0.5;ones(63,1);0.5]*diff(zDomain)/64;

% Scalar basis: explicit count, measurements, then an independent decision.
scalarProblem = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.neumann(),bottomBoundary=IMBoundaryCondition.neumann());
scalar = solver.solveEVP(scalarProblem,nModes=4);
scalarCollection = IMBasisCollection({scalar});
scalarAssessment = scalarCollection.assess(z,weights,prefixColumnCounts=1:4);
scalarDecision = checkBasisAssessment(scalarAssessment,gramTolerance=1e-2);
assert(scalarDecision.requestedColumnCount == 4);

% Aligned F/G, APV, mean-density-anomaly, waves, and inertial kappa=0.
% Each uses its own metric and independently requested scientific count.
problems = {IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain), ...
    IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g0=0,gd=0), ...
    IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=zDomain,g0=0,gd=0), ...
    IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=1e-4), ...
    IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=0)};
counts = [4 4 4 4 3]; % Inertial count is independent of the wave count.
familyAssessments = cell(size(problems));
for iFamily = 1:numel(problems)
    basis = solver.solveEVP(problems{iFamily},nModes=counts(iFamily));
    collection = IMBasisCollection({basis});
    variable = string(basis.evp.formulation);
    values = collection.evaluate(z,variable=variable);
    derivatives = collection.evaluate(z,variable=variable,derivativeOrder=1);
    familyAssessments{iFamily} = collection.assess(z,weights,variable=variable);
    assert(size(values,2) == counts(iFamily) && isequal(size(values),size(derivatives)));
end

% Closed-form bases use the same assessment entry without a spectral solver.
solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain);
analytical = solution.internalModes(problems{1},nModes=4);
analyticalCollection = IMBasisCollection({analytical});
analyticalAssessment = analyticalCollection.assess(z,weights,variable="G");
% Its target Gram integration provenance is explicit, not an exactness claim.
analyticalProvenance = analyticalAssessment.identity.projectionProvenance;

% Repeated wave requests share the same continuous stored basis exactly.
wave1 = solver.solveEVP(problems{4},nModes=4);
wave2 = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=2e-4),nModes=4);
waves = IMBasisCollection({wave1,wave2},kappa=[2e-4 1e-4 2e-4],basisIndex=[2 1 2]);
waveG = waves.evaluate(z,variable="G",pages=[3 1]);
assert(isequal(waveG(:,:,1),waveG(:,:,2)));

% Numerical and analytical zero-APV bases preserve endpoint/page identity.
boundaryProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,f0=1e-4,k=[1e-4 2e-4]);
numericalBoundary = solver.solveGeostrophicZeroAPVModes(boundaryProblem);
analyticalBoundary = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=1e-4).geostrophicZeroAPVModesAtWavenumber([1e-4 2e-4]);
for boundary = {numericalBoundary,analyticalBoundary}
    collection = IMBasisCollection(boundary,basisIndex=[1 1],sourcePage=[2 1]);
    boundaryF = collection.evaluate(z,variable="F");
    boundaryFz = collection.evaluate(z,variable="F",derivativeOrder=1);
    % Scalar projection is unavailable for boundary coordinates; their physical
    % source dual must be supplied explicitly when an assessment needs it.
    assert(collection.metadata.columnKind == "endpoint");
end
