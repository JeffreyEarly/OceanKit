function decision = checkBasisAssessment(assessment,options)
% Apply explicit tolerances without changing the assessed basis.
%
% Disabled quantities do not enter the decision. Missing or
% unqualified enabled measurements are inconclusive. A measured
% failure rejects even when another enabled quantity is inconclusive.
% Prefix recommendations are diagnostic only; requestedColumnCount
% and the stored projection are never reduced.
%
% - Topic: Apply acceptance policies
% - Declaration: decision = checkBasisAssessment(assessment,options)
arguments (Input)
    assessment (1,1) struct
    options.gramTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
    options.leakageTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
    options.quadraticAliasingTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
end
arguments (Output)
    decision (1,1) struct
end
if ~all(isfield(assessment,["projection","prefixColumnCounts","measurements"])) || ~isa(assessment.projection,"IMProjection") || ~istable(assessment.measurements)
    error("checkBasisAssessment:InvalidAssessment","Supply the result of projection.assess or collection.assess.");
end
counts = assessment.prefixColumnCounts;
measurements = assessment.measurements;
if ~isnumeric(counts) || ~isrow(counts) || isempty(counts) || ~isreal(counts) || any(~isfinite(counts) | counts < 1 | fix(counts) ~= counts) || any(diff(counts) <= 0) || counts(end) ~= assessment.projection.columnCount || ~all(ismember(["columnCount","quantity","value","status"],string(measurements.Properties.VariableNames)))
    error("checkBasisAssessment:InvalidAssessment","Preserve the assessment's increasing prefixColumnCounts and measurement columns.");
end
names = ["gram","leakage","quadraticAliasing"];
for count = counts
    for quantity = names
        rows = measurements.columnCount == count & measurements.quantity == quantity;
        if nnz(rows) ~= 1
            error("checkBasisAssessment:InvalidAssessment","Supply exactly one measurement row per examined column count and quantity.");
        end
    end
end
tolerances = {options.gramTolerance,options.leakageTolerance,options.quadraticAliasingTolerance};
if any(cellfun(@numel,tolerances) > 1) || all(cellfun(@isempty,tolerances))
    error("checkBasisAssessment:InvalidPolicy","Supply at least one scalar tolerance; omit quantities that are not acceptance requirements.");
end
outcomes = repmat("accepted",numel(assessment.prefixColumnCounts),1);
for iPrefix = 1:numel(assessment.prefixColumnCounts)
    for iQuantity = 1:numel(names)
        if isempty(tolerances{iQuantity})
            continue
        end
        row = assessment.measurements(assessment.measurements.columnCount == assessment.prefixColumnCounts(iPrefix) & assessment.measurements.quantity == names(iQuantity),:);
        if row.status ~= "measured" || isnan(row.value)
            if outcomes(iPrefix) ~= "rejected"
                outcomes(iPrefix) = "inconclusive";
            end
        elseif row.value > tolerances{iQuantity}
            outcomes(iPrefix) = "rejected";
        end
    end
end
consecutive = cumprod(double(outcomes == "accepted")) ~= 0;
acceptedCounts = assessment.prefixColumnCounts(consecutive);
if isempty(acceptedCounts)
    largestExaminedAcceptedPrefix = 0;
else
    largestExaminedAcceptedPrefix = acceptedCounts(end);
end
overallStatus = outcomes(end);
if any(outcomes == "rejected")
    overallStatus = "rejected";
elseif any(outcomes == "inconclusive")
    overallStatus = "inconclusive";
end
decision = struct(status=overallStatus,requestedBandStatus=outcomes(end),requestedColumnCount=assessment.projection.columnCount,largestExaminedAcceptedPrefix=largestExaminedAcceptedPrefix,allPrefixesExamined=isequal(assessment.prefixColumnCounts,1:assessment.projection.columnCount),prefixDecisions=table(assessment.prefixColumnCounts(:),outcomes,consecutive,VariableNames=["columnCount","status","consecutiveAccepted"]),tolerances=options);
end
