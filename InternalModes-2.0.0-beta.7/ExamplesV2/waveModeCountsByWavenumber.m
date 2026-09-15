function result = waveModeCountsByWavenumber(options)
% Plot linear wave-mode resolution on a fixed physical sampling grid.
%
% Run with InternalModes and its dependencies on the path:
%   addpath ExamplesV2
%   result = waveModeCountsByWavenumber();
%   result = waveModeCountsByWavenumber(Nz=33,outputDirectory="my-wave-counts");
%
% Every distinct positive Fourier radius of the stated periodic domain is
% solved at two explicit EVP resolutions. Counts include the external mode.
% All counts are contiguous prefixes of the requested scientific modes;
% neither solving nor assessment changes those modes. The grid criterion is
% the G Gram discrepancy in its physical wave metric, including the surface
% endpoint. It measures linear sampling, not quadratic-product aliasing.
arguments (Input)
    options.Lxy (1,2) double {mustBePositive,mustBeFinite} = [1000 1000]
    options.Nxy (1,2) double {mustBeInteger,mustBePositive} = [16 16]
    options.D (1,1) double {mustBePositive,mustBeFinite} = 1000
    options.N0 (1,1) double {mustBePositive,mustBeFinite} = 1e-2
    options.scaleDepth (1,1) double {mustBePositive,mustBeFinite} = 400
    options.f0 (1,1) double {mustBeReal,mustBeFinite} = 1e-4
    options.g (1,1) double {mustBePositive,mustBeFinite} = 9.81
    options.Nz (1,1) double {mustBeInteger,mustBePositive} = 25
    options.candidateCount (1,1) double {mustBeInteger,mustBePositive} = 24
    options.inertialCount (1,1) double {mustBeInteger,mustBePositive} = 12
    options.nEVP (1,2) double {mustBeInteger,mustBePositive} = [96 144]
    options.convergenceTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-5
    options.gramTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-2
    options.outputDirectory (1,1) string = ""
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on","off"])} = "on"
end
if any(mod(options.Nxy,2)) || any(options.Nxy < 2) || options.Nz < 3
    error("IMExample:InvalidGrid","Use even horizontal counts of at least two and at least three vertical points.");
end
if options.nEVP(2) <= options.nEVP(1) || options.nEVP(1) <= max(options.candidateCount,options.inertialCount)+2
    error("IMExample:InvalidRefinement","Use increasing EVP resolutions exceeding both requested mode counts by at least three.");
end
outputDirectory = options.outputDirectory;
if outputDirectory == "", outputDirectory = string(tempname)+"-wave-counts"; end
if isfolder(outputDirectory) || isfile(outputDirectory)
    error("IMExample:OutputExists","Choose a new outputDirectory; existing results are never overwritten.");
end
N2 = @(z) options.N0^2*exp(2*z/options.scaleDepth);
zDomain = [-options.D 0];
if N2(-options.D) <= options.f0^2
    error("IMExample:UnsupportedStratification","Require N2 > f0^2 throughout this example's domain.");
end
[nx,ny] = ndgrid(-options.Nxy(1)/2:options.Nxy(1)/2-1,-options.Nxy(2)/2:options.Nxy(2)/2-1);
kappa = unique(2*pi*sqrt((nx(:)/options.Lxy(1)).^2+(ny(:)/options.Lxy(2)).^2)).';
surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,f0=options.f0,g=options.g,k=0,surfaceBoundary=surface);
gridSolver = IMSolverSpectral(nEVP=options.Nz,coordinateKind="wkb").configuredForEVP(evp);
[z,weights] = gridSolver.nativeQuadratureRule(zDomain);
weights = weights*options.D/sum(weights);
referenceGrid = IMSolverSpectral(nEVP=2*options.nEVP(2),coordinateKind="wkb").configuredForEVP(evp);
[zReference,referenceWeights] = referenceGrid.nativeQuadratureRule(zDomain);
collections = cell(1,2);
costs = cell(1,2);
for i = 1:2
    solver = IMSolverSpectral(nEVP=options.nEVP(i),coordinateKind="wkb");
    [collections{i},costs{i}] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=zDomain,f0=options.f0,g=options.g, ...
        surfaceBoundary=surface,nModes=options.candidateCount,nInertialModes=options.inertialCount);
end
counts = table(kappa(:),zeros(numel(kappa),1),zeros(numel(kappa),1),zeros(numel(kappa),1),false(numel(kappa),1),false(numel(kappa),1), ...
    VariableNames=["kappa","convergedCount","gridSupportedCount","combinedCount","atCandidateCeiling","inconclusive"]);
measurements = cell(numel(kappa),1);
gramMeasurements = cell(numel(kappa),1);
for page = 1:numel(kappa)
    candidate = collections{1}.bases{collections{1}.basisIndex(page)};
    reference = collections{2}.bases{collections{2}.basisIndex(page)};
    n = numel(candidate.modeNumber);
    coarse = sampledModes(candidate,zReference);
    fine = sampledModes(reference,zReference);
    report = assessModeConvergence(coarse,fine,zReference,referenceWeights);
    measurements{page} = addvars(report.measurements,repmat(kappa(page),height(report.measurements),1),Before=1,NewVariableNames="kappa");
    converged = false(n,1);
    for j = 1:n
        rows = report.measurements.columnLabel == string(candidate.modeNumber(j)) & ismember(report.measurements.quantity,["equivalentDepth","h1"]);
        converged(j) = nnz(rows) == 3 && all(report.measurements.status(rows) == "measured") && all(report.measurements.value(rows) <= options.convergenceTolerance);
    end
    % kappa=0 is an independent inertial inventory with its F metric.
    if kappa(page) == 0
        % The wave catalog has no F diagnostic recipe at general kappa.
        % At zero kappa its physical inertial identity is integral(F_i F_j)=h_i delta_ij.
        projection = IMProjection(candidate.F(z),diag(weights),diag(candidate.h),columnLabels=string(candidate.modeNumber));
        gridReport = projection.assess(prefixColumnCounts=1:n,identity=coarse.identity);
    else
        gridReport = collections{1}.assess(z,weights,page=page,variable="G",prefixColumnCounts=1:n);
    end
    rows = gridReport.measurements.quantity == "gram";
    gram = gridReport.measurements(rows,:);
    gramMeasurements{page} = addvars(gram,repmat(kappa(page),height(gram),1),Before=1,NewVariableNames="kappa");
    gridAccepted = gram.status == "measured" & gram.value <= options.gramTolerance;
    counts.convergedCount(page) = leadingCount(converged);
    counts.gridSupportedCount(page) = leadingCount(gridAccepted);
    counts.combinedCount(page) = min(counts.convergedCount(page),counts.gridSupportedCount(page));
    counts.atCandidateCeiling(page) = counts.combinedCount(page) == n;
    required = ismember(report.measurements.quantity,["equivalentDepth","h1"]);
    counts.inconclusive(page) = any(report.matches.status ~= "matched") || any(report.measurements.status(required) ~= "measured") || any(gram.status ~= "measured");
end
inertial = counts(counts.kappa == 0,:);
waves = counts(counts.kappa > 0,:);
fig = figure(Name="Wave mode counts by wavenumber",Color="w",Visible=options.figureVisible,Position=[100 100 1000 600]);
ax = axes(fig);
plot(ax,waves.kappa,waves.convergedCount,"-",LineWidth=1.8,DisplayName="EVP converged");
hold(ax,"on")
plot(ax,waves.kappa,waves.gridSupportedCount,"--",LineWidth=1.8,DisplayName="Grid supported (G metric)");
plot(ax,waves.kappa,waves.combinedCount,"o-",LineWidth=1.2,MarkerSize=4,DisplayName="Combined usable prefix");
ceiling = waves.convergedCount == options.candidateCount;
plot(ax,waves.kappa(ceiling),waves.convergedCount(ceiling),"^",Color=[0.2 0.2 0.2],MarkerSize=5,DisplayName="EVP candidate ceiling (lower bound)");
if any(waves.atCandidateCeiling)
    plot(ax,waves.kappa(waves.atCandidateCeiling),waves.combinedCount(waves.atCandidateCeiling),"s",Color="k",MarkerSize=7,DisplayName="Combined candidate ceiling");
end
if any(waves.inconclusive)
    plot(ax,waves.kappa(waves.inconclusive),waves.combinedCount(waves.inconclusive),"x",Color="r",MarkerSize=9,DisplayName="Inconclusive evidence");
end
grid(ax,"on")
xlabel(ax,"\kappa (rad m^{-1})")
ylabel(ax,"Leading wave modes, including external mode")
ylim(ax,[0 options.candidateCount+2])
title(ax,sprintf("%.3g × %.3g km domain; %d × %d horizontal points; %d WKB–Chebyshev vertical points",options.Lxy/1000,options.Nxy,options.Nz));
subtitle(ax,sprintf("EVP %d → %d, tolerance %.1g; Gram tolerance %.1g. Inertial: %d/%d combined (separate).",options.nEVP,options.convergenceTolerance,options.gramTolerance,inertial.combinedCount,options.inertialCount));
legend(ax,Location="southwest")
provenance = struct(configuration=options,profile="N2(z)=N0^2 exp(2z/scaleDepth)",surfaceBoundary="G(0)=h Gz(0)", ...
    bottomBoundary="G(-D)=0",sampling="fixed WKB-stretched Chebyshev-Lobatto physical quadrature", ...
    interpretation="Linear mode convergence and sampled Gram evidence only; no quadratic aliasing certification", ...
    referenceQuadraturePointCount=numel(zReference),matlabVersion=string(version),computer=string(computer), ...
    createdAtUTC=string(datetime("now",TimeZone="UTC")),source=sourceProvenance());
mkdir(outputDirectory);
writetable(waves,fullfile(outputDirectory,"wave-counts.csv"));
writetable(inertial,fullfile(outputDirectory,"inertial-counts.csv"));
writetable(vertcat(measurements{:}),fullfile(outputDirectory,"convergence-measurements.csv"));
writetable(vertcat(gramMeasurements{:}),fullfile(outputDirectory,"grid-measurements.csv"));
writetable(table(z,weights),fullfile(outputDirectory,"physical-grid.csv"));
fid = fopen(fullfile(outputDirectory,"provenance.json"),"w");
if fid == -1, error("IMExample:WriteFailed","Cannot write provenance.json in %s.",outputDirectory); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",jsonencode(provenance,PrettyPrint=true));
exportgraphics(fig,fullfile(outputDirectory,"wave-counts.png"),Resolution=180);
exportgraphics(fig,fullfile(outputDirectory,"wave-counts.pdf"),ContentType="vector");
result = struct(waves=waves,inertial=inertial,z=z,weights=weights,provenance=provenance,costs={costs},figure=fig,outputDirectory=outputDirectory);
disp(waves)
disp(inertial)
fprintf("Saved figure, counts, measurements and provenance to %s\n",outputDirectory);
end

function data = sampledModes(basis,z)
factors = basis.normalizationFactors();
dG = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,1)./factors;
dF = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./factors.*basis.h(:).';
identity = struct(family=string(basis.evp.name),columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain,kappa=basis.evp.parameters.k);
provenance = struct(solverClass=string(class(basis.solver)),nEVP=basis.solver.nEVP,coordinateKind="wkb");
data = struct(identity=identity,values=struct(F=basis.F(z),G=basis.G(z)),derivatives=struct(F=dF,G=dG), ...
    equivalentDepths=basis.h,eigenvalues=basis.eigenvalues,provenance=provenance);
end

function count = leadingCount(accepted)
firstFailure = find(~accepted,1);
if isempty(firstFailure), count = numel(accepted); else, count = firstFailure-1; end
end

function source = sourceProvenance()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
manifest = jsondecode(fileread(fullfile(repoRoot,"resources","mpackage.json")));
source = struct(packageVersion=string(manifest.version),revision="unavailable",worktreeStatus="unavailable",solverPath=string(which("IMSolverSpectral")));
if isfolder(fullfile(repoRoot,".git")) || isfile(fullfile(repoRoot,".git"))
    previousDirectory = pwd;
    cleanup = onCleanup(@() cd(previousDirectory));
    cd(repoRoot);
    [status,revision] = system("git rev-parse HEAD");
    if status == 0, source.revision = strtrim(string(revision)); end
    [status,changes] = system("git status --porcelain");
    if status == 0, source.worktreeStatus = strtrim(string(changes)); end
end
end
