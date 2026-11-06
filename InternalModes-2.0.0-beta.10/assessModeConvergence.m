function assessment = assessModeConvergence(candidate,reference,z,weights)
% Compare explicitly prepared resolved modes without solving or selecting modes.
%
% Each input is a scalar struct with identity, values, and provenance. The
% identity contains family, columnLabels, normalization, and zDomain; include
% kappa when the family has a horizontal-wavenumber identity. All identities
% except columnLabels must agree exactly. Labels match exactly in candidate
% order, without nearest-eigenvalue matching. Missing or duplicate labels
% produce inconclusive rows. Callers are responsible for preparing the same
% physical problem, units, normalization, and physical sample points.
%
% values is a struct of real nZ-by-nMode matrices, for example F and G.
% Optional derivatives contains the corresponding physical z-derivatives.
% Optional eigenvalues and equivalentDepths are vectors in column order.
% provenance is a caller-supplied struct recording the provider revision,
% solver, coordinate kind, eigenproblem resolution, and reference rule when
% applicable. The function preserves provenance without certifying it.
%
% One common sign aligns each matched physical mode across every variable
% and derivative. Shape and derivative errors are weighted relative L2 norms.
% The joint positive H1 error is
% $$\frac{\sqrt{\|v-v_r\|_w^2+D^2\|v_z-v_{r,z}\|_w^2}}{\sqrt{\|v_r\|_w^2+D^2\|v_{r,z}\|_w^2}},$$
% where D is the physical domain depth. H1 remains meaningful for a nearly
% constant mode. A zero reference norm yields zero only for zero discrepancy;
% otherwise the error is Inf. Derivative-only relative errors can therefore
% be large for roundoff-sized derivatives; apply H1 tolerances when judging
% joint field/derivative convergence. No cross-mode smallness cutoff is used.
%
% Matching signed infinities in scalar quantities yield zero error. Finite
% versus infinite values, NaNs, and opposite infinities are inconclusive.
% Equal zero scalar quantities yield zero, and nonzero versus zero yields Inf.
% Two-resolution agreement is convergence evidence, not an absolute accuracy
% bound. This operation applies no tolerances, alters no modes, and does not
% qualify the supplied quadrature or reference solve.
%
% ```matlab
% identity = struct(family="waves",columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain,kappa=kappa);
% candidate = struct(identity=identity,values=struct(G=basis.G(z)),equivalentDepths=basis.h,provenance=struct(nEVP=64));
% reference = struct(identity=identity,values=struct(G=refined.G(z)),equivalentDepths=refined.h,provenance=struct(nEVP=128));
% assessment = assessModeConvergence(candidate,reference,z,weights);
% ```
%
% - Topic: Measure mode convergence
% - Declaration: assessment = assessModeConvergence(candidate,reference,z,weights)
% - Parameter candidate: prepared resolved modes, identity, and provenance
% - Parameter reference: independently prepared comparison with explicit provenance
% - Parameter z: increasing physical sample coordinates, nZ by 1
% - Parameter weights: positive physical integration weights, nZ by 1
% - Returns assessment: per-mode measurements, matching, identity, provenance, coverage, and cost
arguments (Input)
    candidate (1,1) struct
    reference (1,1) struct
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite,mustBePositive}
end
arguments (Output)
    assessment (1,1) struct
end
started = tic;
if numel(z) < 2 || numel(z) ~= numel(weights) || any(diff(z) <= 0)
    error("assessModeConvergence:InvalidQuadrature","Supply at least two increasing physical points and one positive weight per point.");
end
candidate = validatePrepared(candidate,numel(z));
reference = validatePrepared(reference,numel(z));
identityFields = union(string(fieldnames(candidate.identity)),string(fieldnames(reference.identity)));
for field = reshape(identityFields,1,[])
    if field ~= "columnLabels" && (~isfield(candidate.identity,field) || ~isfield(reference.identity,field) || ~isequaln(candidate.identity.(field),reference.identity.(field)))
        error("assessModeConvergence:IdentityMismatch","Candidate and reference must have identical %s identity.",field);
    end
end
zDomain = candidate.identity.zDomain;
if any(z < zDomain(1) | z > zDomain(2))
    error("assessModeConvergence:InvalidQuadrature","Physical sample points must lie within identity.zDomain.");
end
variables = string(fieldnames(candidate.values)).';
if ~isequal(sort(variables),sort(string(fieldnames(reference.values)).'))
    error("assessModeConvergence:VariableMismatch","Candidate and reference values must contain the same physical variables.");
end
labels = candidate.identity.columnLabels(:);
nModes = numel(labels);
referenceColumn = nan(nModes,1);
matchStatus = repmat("inconclusive",nModes,1);
orientation = nan(nModes,1);
quantities = ["eigenvalue";"equivalentDepth";repmat(["shape";"derivative";"h1"],numel(variables),1)];
rowVariables = ["";"";reshape(repmat(variables,3,1),[],1)];
rowsPerMode = numel(quantities);
columnLabel = reshape(repelem(labels,rowsPerMode),[],1);
quantity = repmat(quantities,nModes,1);
variable = repmat(rowVariables,nModes,1);
value = nan(nModes*rowsPerMode,1);
status = repmat("inconclusive",size(value));
% Resolve unique labels once, preserving candidate order and ambiguity rules.
[~,~,candidateGroup] = unique(labels);
[~,~,referenceGroup] = unique(reference.identity.columnLabels(:));
candidateMultiplicity = accumarray(candidateGroup,1);
referenceMultiplicity = accumarray(referenceGroup,1);
[found,match] = ismember(labels,reference.identity.columnLabels);
matched = find(found);
matched = matched(candidateMultiplicity(candidateGroup(matched)) == 1 & referenceMultiplicity(referenceGroup(match(matched))) == 1);
match = match(matched);
referenceColumn(matched) = match;
matchStatus(matched) = "matched";
orientation(matched) = commonOrientation(candidate,reference,matched,match,variables,weights);
rows = (matched(:).'-1)*rowsPerMode;
for iScalar = 1:2
    fields = ["eigenvalues","equivalentDepths"];
    field = fields(iScalar);
    if isfield(candidate,field) && isfield(reference,field)
        [value(rows+iScalar),status(rows+iScalar)] = scalarError(candidate.(field)(matched),reference.(field)(match));
    else
        status(rows+iScalar) = "notRequested";
    end
end
for iVariable = 1:numel(variables)
    field = variables(iVariable);
    a = candidate.values.(field)(:,matched);
    b = reference.values.(field)(:,match).*reshape(orientation(matched),1,[]);
    numerator = sum(weights.*abs(a-b).^2,1);
    denominator = sum(weights.*abs(b).^2,1);
    row = rows+3*(iVariable-1)+3;
    value(row) = normRatio(numerator,denominator);
    status(row) = "measured";
    if isfield(candidate.derivatives,field) && isfield(reference.derivatives,field)
        da = candidate.derivatives.(field)(:,matched);
        db = reference.derivatives.(field)(:,match).*reshape(orientation(matched),1,[]);
        derivativeNumerator = sum(weights.*abs(da-db).^2,1);
        derivativeDenominator = sum(weights.*abs(db).^2,1);
        value(row+1) = normRatio(derivativeNumerator,derivativeDenominator);
        value(row+2) = normRatio(numerator+diff(zDomain)^2*derivativeNumerator,denominator+diff(zDomain)^2*derivativeDenominator);
        status(row+1) = "measured";
        status(row+2) = "measured";
    else
        status(row+1) = "notRequested";
        status(row+2) = "notRequested";
    end
end
assessment = struct();
assessment.identity = candidate.identity;
assessment.matches = table(labels,referenceColumn,matchStatus,orientation,VariableNames=["columnLabel","referenceColumn","status","orientation"]);
assessment.measurements = table(columnLabel,quantity,variable,value,status);
assessment.provenance = struct(candidate=candidate.provenance,reference=reference.provenance);
assessment.coverage = struct(z=z,weights=weights,variables=variables,requestedColumnCount=nModes,matchedColumnCount=nnz(matchStatus=="matched"),referenceAccuracy="unverified",quadratureConvergence="unverified",absoluteAccuracyGuarantee=false);
assessment.costs = struct(assessmentSeconds=toc(started));
end

function prepared = validatePrepared(prepared,nZ)
if ~all(isfield(prepared,["identity","values","provenance"])) || ~isstruct(prepared.identity) || ~isscalar(prepared.identity) || ~isstruct(prepared.values) || ~isscalar(prepared.values) || ~isstruct(prepared.provenance) || ~isscalar(prepared.provenance)
    error("assessModeConvergence:InvalidPreparedModes","Each input requires scalar identity, values, and provenance structs.");
end
identity = prepared.identity;
if ~all(isfield(identity,["family","columnLabels","normalization","zDomain"]))
    error("assessModeConvergence:InvalidIdentity","Identity requires family, columnLabels, normalization, and zDomain.");
end
for field = ["family","normalization"]
    if ~(isstring(identity.(field)) || ischar(identity.(field))) || ~isscalar(string(identity.(field))) || ismissing(string(identity.(field))) || strlength(string(identity.(field))) == 0
        error("assessModeConvergence:InvalidIdentity","Identity %s must be a nonempty text scalar.",field);
    end
    identity.(field) = string(identity.(field));
end
if ~isstring(identity.columnLabels) || ~isvector(identity.columnLabels) || isempty(identity.columnLabels) || any(ismissing(identity.columnLabels) | strlength(identity.columnLabels)==0)
    error("assessModeConvergence:InvalidIdentity","columnLabels must be a nonempty string vector of explicit scientific labels.");
end
identity.columnLabels = reshape(identity.columnLabels,1,[]);
if ~isa(identity.zDomain,"double") || ~isequal(size(identity.zDomain),[1 2]) || ~isreal(identity.zDomain) || any(~isfinite(identity.zDomain)) || diff(identity.zDomain) <= 0
    error("assessModeConvergence:InvalidIdentity","zDomain must contain two finite increasing physical endpoints.");
end
if isfield(identity,"kappa") && (~isa(identity.kappa,"double") || ~isscalar(identity.kappa) || ~isreal(identity.kappa) || ~isfinite(identity.kappa) || identity.kappa < 0)
    error("assessModeConvergence:InvalidIdentity","kappa must be an exact finite nonnegative scalar wavenumber.");
end
prepared.identity = identity;
nModes = numel(identity.columnLabels);
if isempty(fieldnames(prepared.values))
    error("assessModeConvergence:InvalidPreparedModes","Supply at least one physical variable in values.");
end
if ~isfield(prepared,"derivatives")
    prepared.derivatives = struct();
elseif ~isstruct(prepared.derivatives) || ~isscalar(prepared.derivatives) || ~all(isfield(prepared.values,fieldnames(prepared.derivatives)))
    error("assessModeConvergence:InvalidPreparedModes","derivatives must be a scalar struct with matching values fields.");
end
for group = ["values","derivatives"]
    for field = string(fieldnames(prepared.(group))).'
        samples = prepared.(group).(field);
        if ~isa(samples,"double") || ~isreal(samples) || ~isequal(size(samples),[nZ nModes]) || any(~isfinite(samples),"all")
            error("assessModeConvergence:InvalidSamples","%s.%s must contain finite real nZ-by-nMode samples.",group,field);
        end
    end
end
for field = ["eigenvalues","equivalentDepths"]
    if isfield(prepared,field)
        values = prepared.(field);
        if ~isa(values,"double") || ~isreal(values) || ~isvector(values) || numel(values) ~= nModes
            error("assessModeConvergence:InvalidScalarQuantity","%s must contain one real value per mode; exceptional values are reported explicitly.",field);
        end
        prepared.(field) = reshape(values,1,[]);
    end
end
end

function orientation = commonOrientation(candidate,reference,iMode,iReference,variables,weights)
% Choose the strongest normalized overlap, using derivatives only if every
% shape overlap is zero. One sign then applies to the entire physical mode.
orientation = ones(1,numel(iMode));
strongest = zeros(1,numel(iMode));
for group = ["values","derivatives"]
    eligible = strongest == 0;
    if ~any(eligible), break; end
    for field = variables
        if ~isfield(candidate.(group),field) || ~isfield(reference.(group),field)
            continue
        end
        a = candidate.(group).(field)(:,iMode);
        b = reference.(group).(field)(:,iReference);
        scale = sqrt(sum(weights.*a.^2,1).*sum(weights.*b.^2,1));
        overlap = sum(weights.*a.*b,1)./scale;
        update = eligible & scale > 0 & abs(overlap) > strongest;
        strongest(update) = abs(overlap(update));
        orientation(update) = sign(overlap(update));
    end
end
end

function [value,status] = scalarError(a,b)
status = repmat("measured",size(a));
value = abs(a-b)./abs(b);
value(a==b) = 0;
inconclusive = isnan(a) | isnan(b) | ((~isfinite(a) | ~isfinite(b)) & a~=b);
value(inconclusive) = NaN;
status(inconclusive) = "inconclusive";
end

function value = normRatio(numerator,denominator)
value = sqrt(numerator./denominator);
value(~(denominator>0) & numerator==0) = 0;
value(~(denominator>0) & numerator~=0) = Inf;
end
