function result = neutralPhysicalZComparison(options)
% Add physical-z candidates to the saved neutral-coordinate comparison.
%
% Add ExamplesV2 to the path after InternalModes and its dependencies:
%   result = neutralPhysicalZComparison();
%   result = neutralPhysicalZComparison(studyDirectory="saved-study");
%
% Solve physical z at 64/128 points against a fresh z1024 reference, using
% the saved 1025-point quadrature and unchanged 512/1024 acceptance mask.
% Candidate-only prefixes and plots omit that mask and remain diagnostics.
% The original study files are never overwritten.
arguments (Input)
    options.studyDirectory (1,1) string = ""
end
arguments (Output)
    result (1,1) struct
end
repoRoot = fileparts(fileparts(mfilename("fullpath")));
studyDirectory = options.studyDirectory;
if studyDirectory == ""
    studyDirectory = string(fullfile(repoRoot,"Documentation","SmoothedBandwidthCoordinateComparison"));
end
inputNames = ["convergence-measurements.csv","reference-validation.csv","physical-quadrature.csv","provenance.json"];
for name = inputNames
    if ~isfile(fullfile(studyDirectory,name))
        error("IMExample:MissingStudyFile","Required original study file is missing: %s.",fullfile(studyDirectory,name));
    end
end
outputNames = ["neutral-physical-z-measurements.csv","neutral-physical-z-matches.csv", ...
    "neutral-physical-z-summary.csv","neutral-physical-z-provenance.json", ...
    "neutral-coordinate-comparison-with-z.png","neutral-mode-errors-with-z.png"];
for name = outputNames
    if isfile(fullfile(studyDirectory,name)) || isfolder(fullfile(studyDirectory,name))
        error("IMExample:OutputExists","Supplement output already exists: %s. Preserve it or choose another studyDirectory.",fullfile(studyDirectory,name));
    end
end
priorProvenance = jsondecode(fileread(fullfile(studyDirectory,"provenance.json")));
validatePriorProvenance(priorProvenance);
priorMeasurements = readtable(fullfile(studyDirectory,"convergence-measurements.csv"),TextType="string");
referenceValidation = readtable(fullfile(studyDirectory,"reference-validation.csv"),TextType="string");
quadrature = readtable(fullfile(studyDirectory,"physical-quadrature.csv"));
priorMeasurements.columnLabel = string(priorMeasurements.columnLabel);
referenceValidation.columnLabel = string(referenceValidation.columnLabel);
z = quadrature.z;
weights = quadrature.weights;
if numel(z) ~= 1025 || numel(weights) ~= 1025 || any(~isfinite(z) | ~isfinite(weights)) ...
        || any(diff(z) <= 0) || any(weights <= 0) || abs(z(1)+1000) > 1e-10 ...
        || abs(z(end)) > 1e-10 || abs(sum(weights)-1000) > 1e-8
    error("IMExample:InvalidStudyQuadrature","The saved rule must contain 1025 increasing physical-z points on [-1000,0] with positive weights summing to 1000 m.");
end
reference = referenceValidation(referenceValidation.profile == "neutralSurfaceLayer",:);
if height(reference) ~= 126 || any(reference.family ~= "hydrostaticGModes") ...
        || any(~isnan(reference.wavelength) | ~isnan(reference.kappa)) ...
        || any(reference.coordinate ~= "reference" | reference.nEVP ~= 512 | reference.referenceNEVP ~= 1024) ...
        || ~isequal(reference.columnLabel,string((1:126).')) ...
        || any(~ismember(reference.accepted,[0 1]) | ~ismember(reference.measured,[0 1]))
    error("IMExample:InvalidReferenceMask","The original neutral reference mask must contain modes 1:126 from the unchanged 512/1024 comparison.");
end
originalCurves = neutralErrorCurves(priorMeasurements);
evp = IMInternalModes.hydrostaticGModes(N2=@neutralN2,zDomain=[-1000 0],rotationRate=7.2921e-5,latitude=45,g=9.81, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
fprintf("Solving the neutral physical-z reference (nEVP=1024, nModes=126).\n");
high = IMSolverSpectral(nEVP=1024,coordinateKind="z").solveEVP(evp,nModes=126);
highSamples = sampledModes(high,z);
measurementTables = cell(2,1);
matchTables = cell(2,1);
summaryRows = cell(2,1);
physicalCurves = NaN(30,2);
for iResolution = 1:2
    resolution = [64 128];
    resolution = resolution(iResolution);
    requestedCount = resolution-2;
    fprintf("Comparing the neutral physical-z candidate (nEVP=%d).\n",resolution);
    basis = IMSolverSpectral(nEVP=resolution,coordinateKind="z").solveEVP(evp,nModes=requestedCount);
    report = assessModeConvergence(sampledModes(basis,z),highSamples,z,weights);
    labels = string(basis.modeNumber(:));
    [accepted,measured] = acceptedModes(report,labels,1e-5);
    [found,index] = ismember(labels,reference.columnLabel);
    referenceValid = false(size(labels));
    referenceValid(found) = logical(reference.accepted(index(found)));
    resolvedPrefix = leadingCount(accepted & referenceValid);
    returnedCount = numel(labels);
    incomplete = returnedCount < requestedCount;
    inconclusive = any(~measured | ~referenceValid) || incomplete;
    tags = ["original reference uncertainty","incomplete candidate inventory","unmeasured candidate modes"];
    flags = [any(~referenceValid),incomplete,any(~measured)];
    qualification = "";
    if any(flags), qualification = join(tags(flags),"; "); end
    atCandidateCeiling = returnedCount == requestedCount && resolvedPrefix == requestedCount;
    summaryRows{iResolution} = struct(profile="neutralSurfaceLayer",family="hydrostaticGModes",coordinate="z", ...
        nEVP=resolution,referenceNEVP=1024,requestedCount=requestedCount,returnedCount=returnedCount, ...
        resolvedPrefix=resolvedPrefix,candidateOnlyDiagnosticPrefix=leadingCount(accepted), ...
        referenceValidatedPrefix=leadingCount(referenceValid),atCandidateCeiling=atCandidateCeiling, ...
        inconclusive=inconclusive,qualification=qualification);
    measurementTables{iResolution} = annotate(report.measurements,resolution);
    matchTables{iResolution} = annotate(report.matches,resolution);
    physicalCurves(:,iResolution) = fH1Curve(report.measurements);
end
measurements = vertcat(measurementTables{:});
matches = vertcat(matchTables{:});
summary = struct2table(vertcat(summaryRows{:}));
raw = IMSolverSpectral(nEVP=128,coordinateKind="wkb").configuredForEVP(evp);
smoothed = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=128,smoothingLength=priorProvenance.smoothing.length).configuredForEVP(evp);
combinedFigure = fullfile(studyDirectory,outputNames(5));
zoomFigure = fullfile(studyDirectory,outputNames(6));
plotComparison(originalCurves,physicalCurves,raw,smoothed,combinedFigure,zoomFigure);
provenance = struct(profile="neutralSurfaceLayer",family="hydrostaticGModes",rotationRate=7.2921e-5,latitude=45,g=9.81, ...
    zDomain=[-1000 0],profileFormula=string(priorProvenance.profileFormulas(end)),boundary="G=0 at both physical endpoints", ...
    candidateCoordinate="z",candidateNEVP=[64 128],candidateCount="nEVP-2",freshReferenceNEVP=1024,freshReferenceCount=126, ...
    convergenceTolerance=1e-5,referenceTolerance=1e-7,referenceMask="Original neutral 512/1024 mask reused without refinement or threshold changes", ...
    countRule="Contiguous prefix; resolvedPrefix applies original mask; candidateOnlyDiagnosticPrefix omits it", ...
    quadrature="Original saved 1025-point physical-z rule reused exactly",quadratureConvergence="unverified", ...
    derivativeRule="Same normalized sampledModes as primary study: physical Gz and Fz=h*Gzz", ...
    plotRule="First 30 F H1 candidate errors against z1024; reference mask omitted; physical z added only to error panels", ...
    originalInputFiles=inputNames,originalCreatedAtUTC=string(priorProvenance.createdAtUTC),smoothing=priorProvenance.smoothing,matlabVersion=string(version), ...
    createdAtUTC=string(datetime("now",TimeZone="UTC")),solverPath=string(which("IMSolverSpectral")), ...
    scope="Neutral-case supplement only; no new reference certification, full matrix rerun, or change to original results");
writetable(measurements,fullfile(studyDirectory,outputNames(1)));
writetable(matches,fullfile(studyDirectory,outputNames(2)));
writetable(summary,fullfile(studyDirectory,outputNames(3)));
fid = fopen(fullfile(studyDirectory,outputNames(4)),"w");
if fid == -1, error("IMExample:WriteFailed","Cannot write supplement provenance in %s.",studyDirectory); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",jsonencode(provenance,PrettyPrint=true));
result = struct(summary=summary,measurements=measurements,matches=matches,provenance=provenance, ...
    referenceValidation=reference,z=z,weights=weights,studyDirectory=studyDirectory, ...
    combinedFigure=string(combinedFigure),zoomFigure=string(zoomFigure));
disp(summary)
end

function validatePriorProvenance(p)
neutral = find(string(p.profileNames) == "neutralSurfaceLayer");
formula = "N2=(5e-3)^2 s(t), t=clip((-z-200)/200,0,1), s=6t^5-15t^4+10t^3";
if p.rotationRate ~= 7.2921e-5 || p.latitude ~= 45 || p.g ~= 9.81 ...
        || ~isequal(p.zDomain(:),[-1000;0]) || ~isequal(p.nEVP(:),[64;128]) ...
        || numel(neutral) ~= 1 || string(p.profileFormulas(neutral)) ~= formula ...
        || string(p.neutralBoundary) ~= "hydrostatic G=0 at both endpoints" ...
        || string(p.candidateCount) ~= "nEVP-2" || p.convergenceTolerance ~= 1e-5 ...
        || p.referenceTolerance ~= 1e-7 || p.referenceCount ~= 126 ...
        || ~isequal(p.referenceRefinements,[384 512;512 1024]) ...
        || ~isscalar(p.smoothing.length) || ~isreal(p.smoothing.length) || ~isfinite(p.smoothing.length) || p.smoothing.length <= 0 ...
        || string(p.smoothing.equation) ~= "Q-Ls^2 Qzz=B; q=sqrt(Q)" ...
        || string(p.smoothing.boundary) ~= "Qz=0 at both physical endpoints"
    error("IMExample:IncompatibleStudy","Use the original Earth-rate, 45-degree, 1000-m neutral hydrostatic study with unchanged counts, gates and smoothing definition.");
end
end

function curves = neutralErrorCurves(measurements)
curves = NaN(30,4);
coordinates = ["sqrtB","smoothedB"];
resolutions = [64 128];
for iCoordinate = 1:2
    for iResolution = 1:2
        rows = measurements.profile == "neutralSurfaceLayer" & measurements.coordinate == coordinates(iCoordinate) ...
            & measurements.nEVP == resolutions(iResolution);
        values = measurements(rows,:);
        if isempty(values) || any(values.family ~= "hydrostaticGModes" | values.referenceNEVP ~= 1024) ...
                || any(~isnan(values.wavelength) | ~isnan(values.kappa))
            error("IMExample:InvalidOriginalMeasurements","Saved neutral raw/smoothed measurements must use the original z1024 reference and hydrostatic case.");
        end
        curves(:,2*(iCoordinate-1)+iResolution) = fH1Curve(values);
    end
end
end

function curve = fH1Curve(measurements)
values = measurements(measurements.quantity == "h1" & measurements.variable == "F",:);
[found,index] = ismember(string((1:30).'),string(values.columnLabel));
if any(~found) || numel(unique(values.columnLabel)) ~= height(values) ...
        || any(values.status(index) ~= "measured" | ~isfinite(values.value(index)))
    error("IMExample:MissingPlotMeasurement","A unique measured F H1 error is required for each of the first 30 modes.");
end
curve = values.value(index);
end

function values = neutralN2(z)
t = min(max((-z-200)/200,0),1);
values = (5e-3)^2*(6*t.^5-15*t.^4+10*t.^3);
end

function data = sampledModes(basis,z)
factors = basis.normalizationFactors();
dG = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,1)./factors;
dF = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./factors.*basis.h(:).';
identity = struct(family=string(basis.evp.name),columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain);
data = struct(identity=identity,values=struct(F=basis.F(z),G=basis.G(z)),derivatives=struct(F=dF,G=dG), ...
    equivalentDepths=basis.h,eigenvalues=basis.eigenvalues, ...
    provenance=struct(solverClass=string(class(basis.solver)),nEVP=basis.solver.nEVP,coordinate="z"));
end

function [accepted,measured] = acceptedModes(report,labels,tolerance)
accepted = false(numel(labels),1);
measured = false(size(accepted));
for j = 1:numel(labels)
    rows = report.measurements.columnLabel == labels(j) & ismember(report.measurements.quantity,["equivalentDepth","h1"]);
    matches = report.matches.columnLabel == labels(j);
    measured(j) = nnz(rows) == 3 && nnz(matches) == 1 && report.matches.status(matches) == "matched" ...
        && all(report.measurements.status(rows) == "measured");
    accepted(j) = measured(j) && all(report.measurements.value(rows) <= tolerance);
end
end

function count = leadingCount(accepted)
firstFailure = find(~accepted,1);
if isempty(firstFailure), count = numel(accepted); else, count = firstFailure-1; end
end

function values = annotate(values,resolution)
metadata = struct(profile="neutralSurfaceLayer",family="hydrostaticGModes",wavelength=NaN,kappa=NaN, ...
    coordinate="z",nEVP=resolution,referenceNEVP=1024);
fields = string(fieldnames(metadata)).';
for field = fliplr(fields)
    values = addvars(values,repmat(metadata.(field),height(values),1),Before=1,NewVariableNames=field);
end
end

function plotComparison(original,physical,raw,smoothed,combinedPath,zoomPath)
colors = [0 .447 .741;.85 .325 .098;.2 .6 .2];
curves = [original physical];
names = ["Raw, n = 64","Raw, n = 128","Smoothed, n = 64","Smoothed, n = 128", ...
    "Physical z, n = 64","Physical z, n = 128"];
figureHandle = figure(Visible="off",Color="w",Position=[100 100 1300 620]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,1,2,TileSpacing="compact",Padding="compact");
monitorAxes = nexttile(layout);
hold(monitorAxes,"on");
patch(monitorAxes,[1e-4 6e-3 6e-3 1e-4],[-200 -200 0 0],[.93 .93 .93],EdgeColor="none",HandleVisibility="off");
monitorLines = gobjects(2,1);
monitorLines(1) = semilogx(monitorAxes,raw.qReference,raw.zReference,Color=colors(1,:),LineWidth=2);
monitorLines(2) = semilogx(monitorAxes,smoothed.qReference,smoothed.zReference,Color=colors(2,:),LineWidth=2);
set(monitorAxes,XScale="log",XLim=[1e-4 6e-3],YLim=[-1000 0],FontSize=12);
grid(monitorAxes,"on");
text(monitorAxes,1.2e-4,-90,"Neutral cap",Color=[.3 .3 .3]);
xlabel(monitorAxes,"Monitor q = dx/dz (s^{-1})");
ylabel(monitorAxes,"Depth z (m)");
title(monitorAxes,"Original raw and smoothed monitors");
legend(monitorAxes,monitorLines,["Raw bandwidth","Smoothed, L_s = "+round(smoothed.smoothingLength)+" m"],Location="southwest");
errorAxes = nexttile(layout);
plotErrorCurves(errorAxes,curves,colors,names);
ylim(errorAxes,[1e-8 1e3]);
title(errorAxes,"Candidate errors: raw, smoothed and physical z");
title(layout,"Neutral surface layer: reference gate omitted from error curves",FontSize=16);
exportgraphics(figureHandle,combinedPath,Resolution=180);
clear cleanup
zoomHandle = figure(Visible="off",Color="w",Position=[100 100 900 600]);
cleanup = onCleanup(@() close(zoomHandle));
zoomAxes = axes(zoomHandle);
zoomCurves = curves(:,3:6);
plotErrorCurves(zoomAxes,zoomCurves,colors(2:3,:),names(3:6));
positive = zoomCurves(isfinite(zoomCurves) & zoomCurves > 0);
limits = [min(1e-6,10^floor(log10(min(positive)))),max(1e-1,10^ceil(log10(max(positive))))];
ylim(zoomAxes,limits);
title(zoomAxes,"Smoothed and physical-z errors; raw curves omitted");
exportgraphics(zoomHandle,zoomPath,Resolution=180);
end

function plotErrorCurves(ax,curves,colors,names)
hold(ax,"on");
lines = gobjects(size(curves,2),1);
styles = ["--","-"];
for j = 1:size(curves,2)
    lines(j) = semilogy(ax,1:30,curves(:,j),Color=colors(ceil(j/2),:), ...
        LineStyle=styles(mod(j-1,2)+1),LineWidth=2);
end
set(ax,YScale="log",XLim=[1 30],FontSize=12);
grid(ax,"on");
grid(ax,"minor");
yline(ax,1e-5,":","10^{-5} target",Color=[.25 .25 .25],HandleVisibility="off");
xlabel(ax,"Hydrostatic mode number");
ylabel(ax,"F H^1 error against z1024 reference");
legend(ax,lines,names,Location="best",FontSize=10);
end
