function kdeComparison = buildSite1BootstrapStrainKdeComparison(comparison, timeDays, options)
arguments
    comparison struct
    timeDays {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonempty, mustBeVector}
    options.pctTarget (:,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.pctTarget, 0), mustBeLessThanOrEqual(options.pctTarget, 1)} = (0.8:-0.1:0.1).'
    options.summaryMass (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.summaryMass, 0), mustBeLessThanOrEqual(options.summaryMass, 1)} = 0.68
    options.boundPaddingFactor (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1.1
    options.minimumHalfWidth (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1e-6
end

timeDays = reshape(timeDays, [], 1);
site = requiredField(comparison, "site", "comparison");
queryDays = reshape(requiredField(site, "queryTimes", "comparison.site") / 86400, [], 1);
if isempty(queryDays)
    error("GriddedStreamfunction:MissingComparisonQueryTimes", ...
        "comparison.site.queryTimes must be nonempty.");
end

if any(timeDays < queryDays(1) | timeDays > queryDays(end))
    error("GriddedStreamfunction:ComparisonTimeOutOfRange", ...
        "Requested timeDays must lie inside [%.6g, %.6g] days.", queryDays(1), queryDays(end));
end

timeIndices = zeros(size(timeDays));
for iTime = 1:numel(timeDays)
    [~, timeIndices(iTime)] = min(abs(queryDays - timeDays(iTime)));
end

if numel(unique(timeIndices)) ~= numel(timeIndices)
    error("GriddedStreamfunction:DuplicateComparisonTimeIndex", ...
        "Requested timeDays must map to distinct query samples.");
end

matchedTimeDays = queryDays(timeIndices);
methodKeys = ["old"; "new"];
methodLabels = ["Oscroft"; "Gridded Streamfunction"];

allCoordinates = zeros(0, 1);
for iMethod = 1:numel(methodKeys)
    methodKey = char(methodKeys(iMethod));
    methodComparison = requiredField(comparison, methodKey, "comparison");
    sharedComparison = requiredField(methodComparison, "shared", sprintf("comparison.%s", methodKey));
    commonSummary = requiredField(sharedComparison, "commonSummary", sprintf("comparison.%s.shared", methodKey));
    fullCommonSummary = requiredField(sharedComparison, "fullCommonSummary", sprintf("comparison.%s.shared", methodKey));

    sigma_n = requiredField(commonSummary, "sigma_n", sprintf("comparison.%s.shared.commonSummary", methodKey));
    sigma_s = requiredField(commonSummary, "sigma_s", sprintf("comparison.%s.shared.commonSummary", methodKey));
    fullSigmaN = requiredField(fullCommonSummary, "sigma_n", sprintf("comparison.%s.shared.fullCommonSummary", methodKey));
    fullSigmaS = requiredField(fullCommonSummary, "sigma_s", sprintf("comparison.%s.shared.fullCommonSummary", methodKey));

    sigmaNSelection = sigma_n(timeIndices, :);
    sigmaSSelection = sigma_s(timeIndices, :);
    allCoordinates = [allCoordinates; sigmaNSelection(:); sigmaSSelection(:); fullSigmaN(timeIndices); fullSigmaS(timeIndices)]; %#ok<AGROW>
end

halfWidth = options.boundPaddingFactor * max(abs(allCoordinates), [], "all");
if ~isfinite(halfWidth) || halfWidth <= 0
    halfWidth = options.minimumHalfWidth;
end
halfWidth = max(halfWidth, options.minimumHalfWidth);

minimum = [-halfWidth, -halfWidth];
maximum = [halfWidth, halfWidth];
panels = repmat(emptyStrainKdePanel(), numel(timeDays), numel(methodKeys));

for iTime = 1:numel(timeDays)
    for iMethod = 1:numel(methodKeys)
        methodKey = char(methodKeys(iMethod));
        sharedComparison = comparison.(methodKey).shared;
        commonSummary = sharedComparison.commonSummary;
        fullCommonSummary = sharedComparison.fullCommonSummary;

        samples = [ ...
            reshape(commonSummary.sigma_n(timeIndices(iTime), :), [], 1), ...
            reshape(commonSummary.sigma_s(timeIndices(iTime), :), [], 1)];
        fullPoint = [ ...
            reshape(fullCommonSummary.sigma_n(timeIndices(iTime)), 1, 1), ...
            reshape(fullCommonSummary.sigma_s(timeIndices(iTime)), 1, 1)];

        statistics = KernelDensityEstimate.planarStatisticsFromData( ...
            samples, ...
            referencePoint=fullPoint, ...
            minimum=minimum, ...
            maximum=maximum, ...
            summaryMass=options.summaryMass, ...
            boundPaddingFactor=options.boundPaddingFactor);

        panels(iTime, iMethod) = struct( ...
            "methodKey", methodKeys(iMethod), ...
            "methodLabel", methodLabels(iMethod), ...
            "timeIndex", timeIndices(iTime), ...
            "requestedTimeDays", timeDays(iTime), ...
            "matchedTimeDays", matchedTimeDays(iTime), ...
            "samples", samples, ...
            "statistics", statistics);
    end
end

kdeComparison = struct( ...
    "requestedTimeDays", timeDays, ...
    "matchedTimeDays", matchedTimeDays, ...
    "timeIndices", timeIndices, ...
    "methodKeys", methodKeys, ...
    "methodLabels", methodLabels, ...
    "bounds", struct("minimum", minimum, "maximum", maximum, "halfWidth", halfWidth), ...
    "pctTarget", options.pctTarget, ...
    "summaryMass", options.summaryMass, ...
    "panels", panels);
end

function value = requiredField(source, fieldName, sourceName)
if ~isfield(source, fieldName)
    error("GriddedStreamfunction:MissingComparisonField", ...
        "Missing required field %s.%s.", sourceName, fieldName);
end

value = source.(fieldName);
end

function panel = emptyStrainKdePanel()
panel = struct( ...
    "methodKey", "", ...
    "methodLabel", "", ...
    "timeIndex", NaN, ...
    "requestedTimeDays", NaN, ...
    "matchedTimeDays", NaN, ...
    "samples", zeros(0, 2), ...
    "statistics", struct());
end
