function result = bandwidthCoordinateComparison(options)
% Compare contiguous resolved mode prefixes for N and sqrt(B) stretches.
%
% Add ExamplesV2 to the path after InternalModes and its dependencies:
%   result = bandwidthCoordinateComparison();
%   result = bandwidthCoordinateComparison(outputDirectory="coordinate-study");
%   result = bandwidthCoordinateComparison(includeSmoothed=true);
%
% Three fixed-wavenumber wave profiles use a free surface and rigid bottom.
% An exactly neutral surface layer uses hydrostatic G modes with rigid
% endpoints, avoiding a separate signed wave-branch ordering question.
% The legacy neutral coordinate is recorded as unsupported, not as zero
% resolved modes. An unavailable exact requested inventory is recorded
% explicitly without a smaller retry. Counts stop at the first unaccepted
% scientific mode.
% This study measures EVP convergence; it does not assess a simulation's
% vertical sampling, nonlinear aliasing, or the nontraditional wave EVP.
% The opt-in smoothed comparison uses one fixed physical smoothing length
% for all profiles and adds coordinate layout and calculus diagnostics.
arguments (Input)
    options.outputDirectory (1,1) string = ""
    options.includeSmoothed (1,1) logical = false
    options.smoothingLength (1,1) double {mustBeReal,mustBePositive,mustBeFinite} = 147.77523446738553
end
arguments (Output)
    result (1,1) struct
end
outputDirectory = options.outputDirectory;
if outputDirectory == "", outputDirectory = string(tempname)+"-bandwidth-coordinate"; end
if isfolder(outputDirectory) || isfile(outputDirectory)
    error("IMExample:OutputExists","Choose a new outputDirectory; existing results are never overwritten.");
end
D = 1000;
zDomain = [-D 0];
rotationRate = 7.2921e-5;
latitude = 45;
g = 9.81;
nEVP = [64 128];
wavelengths = [1e4 1e3];
convergenceTolerance = 1e-5;
referenceTolerance = 1e-7;
referenceCount = max(nEVP)-2;
coordinates = ["N","sqrtB"];
if options.includeSmoothed, coordinates = ["N","sqrtB","smoothedB"]; end
profiles = struct(name={"constant","standardExponential","weakDeepExponential","neutralSurfaceLayer"}, ...
    N2={@(z) 1e-4*ones(size(z)),@(z) 1e-4*exp(2*z/400),@(z) 9e-6*exp(2*z/300),@neutralSurfaceN2}, ...
    formula={"N2=1e-4","N2=1e-4 exp(2z/400)","N2=9e-6 exp(2z/300)", ...
    "N2=(5e-3)^2 s(t), t=clip((-z-200)/200,0,1), s=6t^5-15t^4+10t^3"});
quadratureEVP = IMInternalModes.hydrostaticGModes(N2=profiles(1).N2,zDomain=zDomain,rotationRate=rotationRate,latitude=latitude,g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
rule = IMSolverSpectral(nEVP=1025,coordinateKind="z").configuredForEVP(quadratureEVP);
[z,weights] = rule.nativeQuadratureRule(zDomain);
% Aggregations are bounded by 42 candidate rows, seven cases and one retry.
summaryRows = cell(0,1);
measurementTables = cell(0,1);
matchTables = cell(0,1);
referenceChecks = cell(0,1);
coordinateControlRows = cell(0,1);
for iProfile = 1:numel(profiles)
    profile = profiles(iProfile);
    isNeutral = iProfile == numel(profiles);
    if isNeutral, profileWavelengths = NaN; else, profileWavelengths = wavelengths; end
    for wavelength = profileWavelengths
        if isNeutral
            kappa = NaN;
            evp = IMInternalModes.hydrostaticGModes(N2=profile.N2,zDomain=zDomain,rotationRate=rotationRate,latitude=latitude,g=g, ...
                surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
        else
            kappa = 2*pi/wavelength;
            surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            evp = IMInternalModes.waveModesAtWavenumber(N2=profile.N2,zDomain=zDomain,k=kappa, ...
                rotationRate=rotationRate,latitude=latitude,g=g,surfaceBoundary=surface,bottomBoundary=IMBoundaryCondition.dirichlet());
        end
        fprintf("Comparing %s, wavelength %.3g m\n",profile.name,wavelength);
        referenceResolutions = [384 512];
        low = IMSolverSpectral(nEVP=referenceResolutions(1),coordinateKind="z").solveEVP(evp,nModes=referenceCount);
        high = IMSolverSpectral(nEVP=referenceResolutions(2),coordinateKind="z").solveEVP(evp,nModes=referenceCount);
        lowSamples = sampledModes(low,z,"z");
        highSamples = sampledModes(high,z,"z");
        referenceReport = assessModeConvergence(lowSamples,highSamples,z,weights);
        [referenceAccepted,referenceMeasured] = acceptedModes(referenceReport,string(low.modeNumber),referenceTolerance);
        metadata = comparisonMetadata(profile,evp,wavelength,kappa,"reference",referenceResolutions(1),referenceResolutions(2));
        measurementTables{end+1,1} = annotate(referenceReport.measurements,metadata); %#ok<AGROW>
        matchTables{end+1,1} = annotate(referenceReport.matches,metadata); %#ok<AGROW>
        referenceRefined = any(~referenceAccepted) || numel(low.modeNumber) < referenceCount;
        if referenceRefined
            low = high;
            lowSamples = highSamples;
            referenceResolutions = [512 1024];
            high = IMSolverSpectral(nEVP=referenceResolutions(2),coordinateKind="z").solveEVP(evp,nModes=referenceCount);
            highSamples = sampledModes(high,z,"z");
            referenceReport = assessModeConvergence(lowSamples,highSamples,z,weights);
            [referenceAccepted,referenceMeasured] = acceptedModes(referenceReport,string(low.modeNumber),referenceTolerance);
            metadata = comparisonMetadata(profile,evp,wavelength,kappa,"reference",referenceResolutions(1),referenceResolutions(2));
            measurementTables{end+1,1} = annotate(referenceReport.measurements,metadata); %#ok<AGROW>
            matchTables{end+1,1} = annotate(referenceReport.matches,metadata); %#ok<AGROW>
        end
        referenceLabels = string(low.modeNumber(:));
        check = table(referenceLabels,referenceAccepted,referenceMeasured, ...
            VariableNames=["columnLabel","accepted","measured"]);
        referenceChecks{end+1,1} = annotate(check,metadata); %#ok<AGROW>
        for resolution = nEVP
            requestedCount = resolution-2;
            for coordinate = coordinates
                if coordinate == "N"
                    solver = IMBandwidthComparison.LegacyWKBSolver(nEVP=resolution,coordinateKind="wkb");
                elseif coordinate == "sqrtB"
                    solver = IMSolverSpectral(nEVP=resolution,coordinateKind="bandwidth");
                else
                    solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=resolution,smoothingLength=options.smoothingLength);
                end
                row = struct(profile=string(profile.name),family=string(evp.name),wavelength=wavelength,kappa=kappa, ...
                    nEVP=resolution,coordinate=coordinate,requestedCount=requestedCount,returnedCount=NaN, ...
                    resolvedPrefix=NaN,atCandidateCeiling=false,inconclusive=false,status="unsupported", ...
                    referenceLow=referenceResolutions(1),referenceHigh=referenceResolutions(2), ...
                    referenceRefined=referenceRefined,referenceValidatedPrefix=NaN,errorIdentifier="",errorMessage="");
                if options.includeSmoothed
                    row.candidateOnlyDiagnosticPrefix = NaN;
                    control = coordinateControlRow(profile,evp,wavelength,kappa,coordinate,resolution);
                end
                try
                    basis = solver.solveEVP(evp,nModes=requestedCount);
                catch exception
                    if string(exception.identifier) == "IMSolver:InsufficientModes"
                        row.status = "insufficientInventory";
                        row.inconclusive = true;
                    elseif ~(isNeutral && coordinate == "N" && string(exception.identifier) == "IMSolverSpectral:InvalidCoordinate")
                        rethrow(exception)
                    end
                    row.errorIdentifier = string(exception.identifier);
                    row.errorMessage = string(exception.message);
                    summaryRows{end+1,1} = row; %#ok<AGROW>
                    if options.includeSmoothed
                        control.status = row.status;
                        control.errorIdentifier = row.errorIdentifier;
                        coordinateControlRows{end+1,1} = control; %#ok<AGROW>
                    end
                    continue
                end
                candidate = sampledModes(basis,z,coordinate);
                report = assessModeConvergence(candidate,highSamples,z,weights);
                labels = string(basis.modeNumber(:));
                [accepted,measured] = acceptedModes(report,labels,convergenceTolerance);
                [found,index] = ismember(labels,referenceLabels);
                referenceValid = false(size(labels));
                referenceValid(found) = referenceAccepted(index(found));
                row.returnedCount = numel(labels);
                row.resolvedPrefix = leadingCount(accepted & referenceValid);
                row.referenceValidatedPrefix = leadingCount(referenceValid);
                row.atCandidateCeiling = row.returnedCount == requestedCount && row.resolvedPrefix == requestedCount;
                row.inconclusive = any(~measured | ~referenceValid) || row.returnedCount < requestedCount;
                row.status = "assessed";
                if row.inconclusive, row.status = "inconclusive"; end
                if options.includeSmoothed
                    row.candidateOnlyDiagnosticPrefix = leadingCount(accepted);
                    control = assessCoordinateControl(control,basis.solver,z);
                    coordinateControlRows{end+1,1} = control; %#ok<AGROW>
                end
                summaryRows{end+1,1} = row; %#ok<AGROW>
                metadata = comparisonMetadata(profile,evp,wavelength,kappa,coordinate,resolution,referenceResolutions(2));
                measurementTables{end+1,1} = annotate(report.measurements,metadata); %#ok<AGROW>
                matchTables{end+1,1} = annotate(report.matches,metadata); %#ok<AGROW>
            end
        end
    end
end
zeroRotationEVP = IMInternalModes.hydrostaticGModes(N2=profiles(2).N2,zDomain=zDomain,rotationRate=0,latitude=latitude,g=g,surfaceBoundary=IMBoundaryCondition.dirichlet());
legacy = IMBandwidthComparison.LegacyWKBSolver(nEVP=nEVP(1),coordinateKind="wkb").configuredForEVP(zeroRotationEVP);
bandwidth = IMSolverSpectral(nEVP=nEVP(1),coordinateKind="bandwidth").configuredForEVP(zeroRotationEVP);
zeroRotationRecovery = struct(nEVP=nEVP(1),profile=string(profiles(2).name), ...
    identicalCoordinate=isequal(legacy.zNative,bandwidth.zNative) && isequal(legacy.qReference,bandwidth.qReference) && isequal(legacy.qzReference,bandwidth.qzReference));
summary = struct2table(vertcat(summaryRows{:}));
paired = pairedResults(summary);
if options.includeSmoothed
    smoothed = smoothedResults(summary,paired);
    coordinateControls = struct2table(vertcat(coordinateControlRows{:}));
end
measurements = vertcat(measurementTables{:});
matches = vertcat(matchTables{:});
referenceValidation = vertcat(referenceChecks{:});
provenance = struct(rotationRate=rotationRate,latitude=latitude,f0=2*rotationRate*sind(latitude),g=g,zDomain=zDomain, ...
    nEVP=nEVP,wavelengths=wavelengths,profileNames=string({profiles.name}),profileFormulas=string({profiles.formula}), ...
    oldCoordinate="dx/dz=N; only coordinate rotationRate forced to zero", ...
    newCoordinate="dx/dz=sqrt(B), B=sqrt((N2+fHorizontal^2-f0^2)^2+4*f0^2*fHorizontal^2), fHorizontal=2*rotationRate*cosd(latitude)", ...
    waveBoundary="surface Gz=G/h, bottom G=0",neutralBoundary="hydrostatic G=0 at both endpoints", ...
    candidateCount="nEVP-2",countRule="contiguous prefix before first failed or unverified mode", ...
    unavailableInventory="insufficientInventory status with the exact solver error; returned count and resolved prefix remain NaN; no smaller retry", ...
    convergenceTolerance=convergenceTolerance,referenceTolerance=referenceTolerance, ...
    referenceCount=referenceCount,referenceRefinements=[384 512;512 1024], ...
    quadrature="1025-point physical-z Chebyshev-Lobatto rule shared by both coordinates", ...
    quadratureConvergence="unverified",referenceAccuracy="checked by one refinement; no absolute accuracy guarantee", ...
    interpretation="EVP convergence only; no simulation-grid, nonlinear-aliasing, or nontraditional-EVP qualification", ...
    zeroRotationRecovery=zeroRotationRecovery,matlabVersion=string(version),computer=string(computer), ...
    createdAtUTC=string(datetime("now",TimeZone="UTC")),source=sourceProvenance());
if options.includeSmoothed
    seedH = 0.224486019760322;
    seedOnsetCoefficient = 10*(5e-3)^2/200^3;
    provenance.smoothing = struct(length=options.smoothingLength,nominalLength=147.77523446738553, ...
        equation="Q-Ls^2 Qzz=B; q=sqrt(Q)",boundary="Qz=0 at both physical endpoints", ...
        seedEquivalentDepth=seedH,seedGravity=g,seedOnsetCoefficient=seedOnsetCoefficient, ...
        nominalLengthRule="(g*seedEquivalentDepth/seedOnsetCoefficient)^(1/5)", ...
        selection="Previously measured first neutral mode; one requested length fixed across all profiles and resolutions; no accuracy tuning");
    provenance.coordinateControls = struct(qRange="sampled on each configured solver's zReference", ...
        depthBands="cap z>=-200; ramp -400<z<-200; deep z<=-400; geometric bands for nonneutral profiles", ...
        calculus="interpolate physical z on the native grid, then evaluate first and second physical derivatives on common quadrature points excluding native nodes", ...
        weights="absolute and relative error in native physical quadrature weights summing to the domain depth", ...
        diagnosticPrefixes="candidate-only prefixes omit the reference-accuracy gate and are not certified counts");
end
mkdir(outputDirectory);
writetable(summary,fullfile(outputDirectory,"coordinate-comparison.csv"));
writetable(paired,fullfile(outputDirectory,"paired-coordinate-comparison.csv"));
if options.includeSmoothed
    writetable(smoothed,fullfile(outputDirectory,"smoothed-coordinate-comparison.csv"));
    writetable(coordinateControls,fullfile(outputDirectory,"coordinate-controls.csv"));
end
writetable(measurements,fullfile(outputDirectory,"convergence-measurements.csv"));
writetable(matches,fullfile(outputDirectory,"mode-matches.csv"));
writetable(referenceValidation,fullfile(outputDirectory,"reference-validation.csv"));
writetable(table(z,weights),fullfile(outputDirectory,"physical-quadrature.csv"));
fid = fopen(fullfile(outputDirectory,"provenance.json"),"w");
if fid == -1, error("IMExample:WriteFailed","Cannot write provenance.json in %s.",outputDirectory); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",jsonencode(provenance,PrettyPrint=true));
result = struct(paired=paired,summary=summary,measurements=measurements,matches=matches,referenceValidation=referenceValidation, ...
    z=z,weights=weights,provenance=provenance,outputDirectory=outputDirectory);
if options.includeSmoothed
    result.smoothed = smoothed;
    result.coordinateControls = coordinateControls;
end
disp(paired)
if options.includeSmoothed, disp(smoothed); end
fprintf("Saved comparison, measurements and provenance to %s\n",outputDirectory);
end

function values = neutralSurfaceN2(z)
t = min(max((-z-200)/200,0),1);
values = (5e-3)^2*(6*t.^5-15*t.^4+10*t.^3);
end

function data = sampledModes(basis,z,coordinate)
dG = basis.Gz(z);
dF = basis.evaluate(z,variable="F",derivativeOrder=1);
identity = struct(family=string(basis.evp.name),columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain);
if isfield(basis.evp.parameters,"k"), identity.kappa = basis.evp.parameters.k; end
provenance = struct(solverClass=string(class(basis.solver)),nEVP=basis.solver.nEVP,coordinate=coordinate);
data = struct(identity=identity,values=struct(F=basis.F(z),G=basis.G(z)),derivatives=struct(F=dF,G=dG), ...
    equivalentDepths=basis.h,eigenvalues=basis.eigenvalues,provenance=provenance);
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

function paired = pairedResults(summary)
old = summary(summary.coordinate == "N",:);
new = summary(summary.coordinate == "sqrtB",:);
keys = ["profile","family","wavelength","kappa","nEVP"];
if ~isequaln(old(:,keys),new(:,keys))
    error("IMExample:UnpairedCoordinates","Old and new coordinate rows must describe the same ordered cases.");
end
qualification = strings(height(old),1);
for j = 1:height(old)
    % Each paired case has at most eight qualification tags.
    tags = strings(0,1);
    if old.status(j) == "unsupported", tags(end+1,1) = "legacy unsupported"; end %#ok<AGROW>
    if new.status(j) == "unsupported", tags(end+1,1) = "new coordinate unsupported"; end %#ok<AGROW>
    if old.status(j) == "insufficientInventory" || new.status(j) == "insufficientInventory"
        tags(end+1,1) = "requested candidate inventory unavailable"; %#ok<AGROW>
    end
    if old.inconclusive(j) || new.inconclusive(j), tags(end+1,1) = "inconclusive evidence"; end %#ok<AGROW>
    if old.referenceValidatedPrefix(j) < old.returnedCount(j) || new.referenceValidatedPrefix(j) < new.returnedCount(j)
        tags(end+1,1) = "reference uncertainty"; %#ok<AGROW>
    end
    if old.returnedCount(j) < old.requestedCount(j) || new.returnedCount(j) < new.requestedCount(j)
        tags(end+1,1) = "incomplete candidate inventory"; %#ok<AGROW>
    end
    if old.atCandidateCeiling(j), tags(end+1,1) = "old at candidate ceiling"; end %#ok<AGROW>
    if new.atCandidateCeiling(j), tags(end+1,1) = "new at candidate ceiling"; end %#ok<AGROW>
    if isempty(tags), qualification(j) = "assessed"; else, qualification(j) = join(tags,"; "); end
end
paired = table(old.profile,old.wavelength,old.kappa,old.nEVP,old.resolvedPrefix,new.resolvedPrefix, ...
    new.resolvedPrefix-old.resolvedPrefix,qualification, ...
    VariableNames=["profile","wavelength","kappa","nEVP","oldResolvedCount","newResolvedCount","difference","qualification"]);
end

function smoothed = smoothedResults(summary,paired)
old = summary(summary.coordinate == "N",:);
raw = summary(summary.coordinate == "sqrtB",:);
smooth = summary(summary.coordinate == "smoothedB",:);
keys = ["profile","family","wavelength","kappa","nEVP"];
if ~isequaln(old(:,keys),smooth(:,keys)) || ~isequaln(raw(:,keys),smooth(:,keys))
    error("IMExample:UnpairedCoordinates","All three coordinate rows must describe the same ordered cases.");
end
qualification = strings(height(smooth),1);
for j = 1:height(smooth)
    % At most six smoothed/diagnostic tags extend the old/raw qualification.
    tags = paired.qualification(j);
    if smooth.inconclusive(j), tags(end+1,1) = "smoothed evidence inconclusive"; end %#ok<AGROW>
    if smooth.status(j) == "insufficientInventory", tags(end+1,1) = "smoothed requested inventory unavailable"; end %#ok<AGROW>
    if smooth.referenceValidatedPrefix(j) < smooth.returnedCount(j)
        tags(end+1,1) = "smoothed reference uncertainty"; %#ok<AGROW>
    end
    if smooth.returnedCount(j) < smooth.requestedCount(j)
        tags(end+1,1) = "smoothed incomplete candidate inventory"; %#ok<AGROW>
    end
    if smooth.atCandidateCeiling(j), tags(end+1,1) = "smoothed at candidate ceiling"; end %#ok<AGROW>
    tags(end+1,1) = "candidate-only prefixes are diagnostic; reference gate omitted"; %#ok<AGROW>
    qualification(j) = join(tags,"; ");
end
smoothed = table(old.profile,old.wavelength,old.kappa,old.nEVP,old.resolvedPrefix,raw.resolvedPrefix,smooth.resolvedPrefix, ...
    raw.resolvedPrefix-old.resolvedPrefix,smooth.resolvedPrefix-raw.resolvedPrefix,smooth.resolvedPrefix-old.resolvedPrefix, ...
    old.candidateOnlyDiagnosticPrefix,raw.candidateOnlyDiagnosticPrefix,smooth.candidateOnlyDiagnosticPrefix,qualification, ...
    VariableNames=["profile","wavelength","kappa","nEVP","oldResolvedCount","rawResolvedCount","smoothedResolvedCount", ...
    "rawMinusOld","smoothedMinusRaw","smoothedMinusOld","oldCandidateOnlyDiagnosticPrefix", ...
    "rawCandidateOnlyDiagnosticPrefix","smoothedCandidateOnlyDiagnosticPrefix","qualification"]);
end

function control = coordinateControlRow(profile,evp,wavelength,kappa,coordinate,nEVP)
control = struct(profile=string(profile.name),family=string(evp.name),wavelength=wavelength,kappa=kappa, ...
    coordinate=coordinate,nEVP=nEVP,status="unsupported",minQ=NaN,maxQ=NaN, ...
    neutralBandNodeCount=NaN,rampBandNodeCount=NaN,deepBandNodeCount=NaN,capCoordinateFraction=NaN, ...
    offGridPointCount=NaN,linearFirstDerivativeMaxError=NaN,linearSecondDerivativeMaxError=NaN, ...
    weightsSumAbsoluteError=NaN,weightsSumRelativeError=NaN,errorIdentifier="");
end

function control = assessCoordinateControl(control,solver,z)
control.status = "measured";
control.minQ = min(solver.qReference);
control.maxQ = max(solver.qReference);
control.neutralBandNodeCount = nnz(solver.zNative >= -200);
control.rampBandNodeCount = nnz(solver.zNative > -400 & solver.zNative < -200);
control.deepBandNodeCount = nnz(solver.zNative <= -400);
xBottom = solver.xOfZ(solver.zDomain(1));
xSurface = solver.xOfZ(solver.zDomain(2));
control.capCoordinateFraction = (xSurface-solver.xOfZ(-200))/(xSurface-xBottom);
gridTolerance = 100*eps(max(1,max(abs(z))));
offGrid = all(abs(z(:)-solver.zNative(:).') > gridTolerance,2);
zControl = z(offGrid);
control.offGridPointCount = numel(zControl);
linearCoefficients = solver.T\solver.zNative;
d1 = solver.evaluatePhysicalDerivative(linearCoefficients,zControl,1);
d2 = solver.evaluatePhysicalDerivative(linearCoefficients,zControl,2);
control.linearFirstDerivativeMaxError = max(abs(d1-1));
control.linearSecondDerivativeMaxError = max(abs(d2));
[~,nativeWeights] = solver.nativeQuadratureRule(solver.zDomain);
D = diff(solver.zDomain);
control.weightsSumAbsoluteError = abs(sum(nativeWeights)-D);
control.weightsSumRelativeError = control.weightsSumAbsoluteError/D;
end

function metadata = comparisonMetadata(profile,evp,wavelength,kappa,coordinate,nEVP,referenceNEVP)
metadata = struct(profile=string(profile.name),family=string(evp.name),wavelength=wavelength,kappa=kappa, ...
    coordinate=coordinate,nEVP=nEVP,referenceNEVP=referenceNEVP);
end

function values = annotate(values,metadata)
fields = string(fieldnames(metadata)).';
for field = fliplr(fields)
    values = addvars(values,repmat(metadata.(field),height(values),1),Before=1,NewVariableNames=field);
end
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
