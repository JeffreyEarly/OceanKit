function assessment = assess(projection,options)
% Measure fixed projection columns without changing them or applying tolerances.
%
% Return a struct containing projection, identity, columnKind,
% prefixColumnCounts, measurements, coverage, and costs. Measurement rows
% record Gram, leakage, and quadraticAliasing values and reference status.
% Products supply signed referencePairings; every prefix solves its own target
% system. Error magnitudes use the separately supplied positive majorant.
% - Topic: Measure projection quality
% - Parameter options.prefixColumnCounts: increasing column counts ending at the full band
% - Parameter options.identity: scientific family and page identity
% - Parameter options.columnKind: mode or endpoint; endpoints are never prefixed
% - Parameter options.leakage: sampledValues, normSquared, columnLabels
% - Parameter options.products: sampledValues, referencePairings, normSquared, inputColumns (2 by n), inputLabels (2 by n), referenceStatus; zero inputColumns denote another family
% - Parameter options.coverage: caller supplied sampled coverage
% - Parameter options.constructionSeconds: separate construction cost
% - Returns assessment: measurements and evidence; checkBasisAssessment applies tolerances
arguments (Input)
    projection (1,1) IMProjection
    options.identity (1,1) struct = struct()
    options.columnKind (1,1) string {mustBeMember(options.columnKind,["mode","endpoint"])} = "mode"
    options.prefixColumnCounts (1,:) double {mustBeInteger,mustBePositive} = []
    options.leakage struct = struct.empty
    options.products struct = struct.empty
    options.coverage (1,1) struct = struct()
    options.constructionSeconds (1,1) double = NaN
end
arguments (Output)
    assessment (1,1) struct
end
assessment = struct();
started = tic;
nColumns = projection.columnCount;
if nColumns == 0
    error("IMProjection:EmptyProjection","Supply a projection with at least one scientific column.");
end
counts = options.prefixColumnCounts;
if isempty(counts)
    counts = nColumns;
end
if any(counts > nColumns) || numel(unique(counts)) ~= numel(counts) || any(diff(counts) <= 0) || counts(end) ~= nColumns
    error("IMProjection:InvalidPrefixes","prefixColumnCounts must increase without duplicates, end at the requested column count, and stay within that count.");
end
if options.columnKind == "endpoint" && ~isequal(counts,nColumns)
    error("IMProjection:EndpointPrefix","Endpoint coordinates are assessed as the complete requested set, never as modal prefixes.");
end
leakage = options.leakage;
products = options.products;
if ~isempty(leakage)
    requireFields(leakage,["sampledValues","normSquared","columnLabels"],"leakage");
    if size(leakage.sampledValues,1) ~= projection.sampleCount || numel(leakage.normSquared) ~= size(leakage.sampledValues,2) || numel(leakage.columnLabels) ~= size(leakage.sampledValues,2)
        error("IMProjection:InvalidLeakage","Leakage samples, norms, and labels must describe the same check columns.");
    end
    leakage.columnLabels = reshape(string(leakage.columnLabels),1,[]);
    leakage.normSquared = reshape(leakage.normSquared,1,[]);
    if any(~isfinite(leakage.normSquared) | leakage.normSquared <= 0) || numel(unique(leakage.columnLabels)) ~= numel(leakage.columnLabels)
        error("IMProjection:InvalidLeakage","Leakage check columns require distinct labels and finite positive norm squares.");
    end
end
if ~isempty(products)
    requireFields(products,["sampledValues","referencePairings","normSquared","inputColumns","inputLabels","referenceStatus"],"products");
    nProducts = size(products.sampledValues,2);
    if size(products.sampledValues,1) ~= projection.sampleCount || ~isequal(size(products.referencePairings),[nColumns nProducts]) || numel(products.normSquared) ~= nProducts || ~isequal(size(products.inputColumns),[2 nProducts]) || ~isequal(size(products.inputLabels),[2 nProducts])
        error("IMProjection:InvalidProducts","Product samples, reference pairings, norms, and input identities must describe the same products and requested output columns.");
    end
    if any(~isfinite(products.inputColumns) | products.inputColumns < 0 | products.inputColumns > nColumns | fix(products.inputColumns) ~= products.inputColumns,"all")
        error("IMProjection:InvalidProducts","Product inputColumns must contain ordinals within the output family; zero denotes an input outside this output family.");
    end
    if any(~isfinite(products.sampledValues),"all") || any(~isfinite(products.referencePairings),"all") || ~isreal(products.normSquared) || any(~isfinite(products.normSquared) | products.normSquared < 0,"all")
        error("IMProjection:InvalidProducts","Products and signed reference pairings must be finite; norm squares must be real, finite, and nonnegative.");
    end
    products.normSquared = reshape(products.normSquared,1,[]);
    products.inputLabels = string(products.inputLabels);
    products.referenceStatus = string(products.referenceStatus);
    if ~isscalar(products.referenceStatus) || ~ismember(products.referenceStatus,["qualified","inconclusive","unverified"])
        error("IMProjection:InvalidReferenceStatus","Product referenceStatus must be qualified, inconclusive, or unverified.");
    end
    if ~isfield(products,"referenceProvenance")
        products.referenceProvenance = "unspecified";
    end
end
nRows = 3*numel(counts);
columnCount = reshape(repelem(counts(:),3),[],1);
quantity = repmat(["gram";"leakage";"quadraticAliasing"],numel(counts),1);
value = nan(nRows,1);
status = repmat("notRequested",nRows,1);
limitingInputI = strings(nRows,1);
limitingInputJ = strings(nRows,1);
referenceStatus = repmat("notRequired",nRows,1);
referenceProvenance = strings(nRows,1);
examinedCount = zeros(nRows,1);
for iPrefix = 1:numel(counts)
    count = counts(iPrefix);
    p = projection.prefix(count);
    row = 3*(iPrefix-1)+1;
    value(row) = p.gramError;
    status(row) = "measured";
    if ~p.supportsGramAssessment
        status(row) = "unsupported";
    end
    examinedCount(row) = count;
    if ~isempty(leakage)
        rejected = ~ismember(leakage.columnLabels,p.columnLabels);
        examinedCount(row+1) = nnz(rejected);
        if any(rejected)
            errors = p.leakage(leakage.sampledValues(:,rejected),leakage.normSquared(rejected));
            [value(row+1),limiting] = max(errors);
            labels = leakage.columnLabels(rejected);
            limitingInputI(row+1) = labels(limiting);
            status(row+1) = "measured";
        else
            status(row+1) = "inconclusive";
        end
    end
    if ~isempty(products)
        selected = all(products.inputColumns <= count,1);
        examinedCount(row+2) = nnz(selected);
        referenceStatus(row+2) = products.referenceStatus;
        referenceProvenance(row+2) = string(products.referenceProvenance);
        if any(selected)
            active = p.activeColumnMask;
            coefficients = zeros(count,nnz(selected));
            pairings = products.referencePairings(1:count,selected);
            if any(pairings(~active,:) ~= 0,"all")
                error("IMProjection:InactiveReferencePairing","Inactive output columns require zero reference pairings; supplied reference data are inconsistent with the projection.");
            end
            target = p.targetGramMatrix(active,active);
            if rank(target) < nnz(active)
                error("IMProjection:SingularReferenceTarget","Continuous reference projection requires a nonsingular target on active columns.");
            end
            coefficients(active,:) = target \ pairings(active,:);
            errors = p.productError(products.sampledValues(:,selected),coefficients,products.normSquared(selected));
            [value(row+2),limiting] = max(errors);
            labels = products.inputLabels(:,selected);
            limitingInputI(row+2) = labels(1,limiting);
            limitingInputJ(row+2) = labels(2,limiting);
            status(row+2) = "measured";
        else
            status(row+2) = "inconclusive";
        end
        if products.referenceStatus ~= "qualified"
            status(row+2) = "inconclusive";
        end
    end
end
assessment.projection = projection;
assessment.identity = options.identity;
assessment.columnKind = options.columnKind;
assessment.prefixColumnCounts = counts;
assessment.measurements = table(columnCount,quantity,value,status,limitingInputI,limitingInputJ,referenceStatus,referenceProvenance,examinedCount);
assessment.coverage = options.coverage;
assessment.coverage.exhaustiveGuarantee = false;
assessment.coverage.superpositionGuarantee = false;
assessment.costs = struct(constructionSeconds=options.constructionSeconds,assessmentSeconds=toc(started));
end

function requireFields(value,fields,name)
if ~isscalar(value) || ~all(isfield(value,fields))
    error("IMProjection:InvalidPreparedData","%s must be a scalar struct containing %s.",name,join(fields,", "));
end
end
