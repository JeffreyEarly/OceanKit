function comparison = compareProfileHotspots(baseline, candidate, options)
% Compare two hotspot analyses or profiler outputs.
%
% Use `compareProfileHotspots` to compare same-workload before-and-after
% profiling runs with absolute and percent deltas for functions and
% executed lines.
%
% ```matlab
% diff = compareProfileHotspots(baselineAnalysis, candidateAnalysis);
% ```

arguments
    baseline (1,1) struct
    candidate (1,1) struct
    options.projectRoots = string(pwd)
    options.maxFunctions (1,1) double {mustBeInteger, mustBePositive} = 25
    options.maxLines (1,1) double {mustBeInteger, mustBePositive} = 25
    options.shouldIncludeToolboxInPriorityTargets (1,1) logical = false
end

baselineAnalysis = coerceAnalysis(baseline, options.projectRoots, ...
    options.maxFunctions, options.maxLines, options.shouldIncludeToolboxInPriorityTargets, "baseline");
candidateAnalysis = coerceAnalysis(candidate, options.projectRoots, ...
    options.maxFunctions, options.maxLines, options.shouldIncludeToolboxInPriorityTargets, "candidate");

functionDiffs = compareFunctionTables(baselineAnalysis.functionMetrics, candidateAnalysis.functionMetrics);
lineDiffs = compareLineTables(baselineAnalysis.lineMetrics, candidateAnalysis.lineMetrics);

comparison = struct( ...
    "baselineLabel", chooseLabel(baselineAnalysis.label, "baseline"), ...
    "candidateLabel", chooseLabel(candidateAnalysis.label, "candidate"), ...
    "functionDiffs", functionDiffs, ...
    "lineDiffs", lineDiffs, ...
    "regressions", buildComparisonTargets(functionDiffs, lineDiffs, true, ...
        options.maxFunctions, options.maxLines, options.shouldIncludeToolboxInPriorityTargets), ...
    "improvements", buildComparisonTargets(functionDiffs, lineDiffs, false, ...
        options.maxFunctions, options.maxLines, options.shouldIncludeToolboxInPriorityTargets), ...
    "summary", strings(0, 1));

comparison.summary = buildComparisonSummary(comparison);

if nargout == 0
    disp("Profile comparison summary:")
    for iSummary = 1:numel(comparison.summary)
        disp("  - " + comparison.summary(iSummary))
    end
    if ~isempty(comparison.regressions)
        disp("Regressions:")
        disp(comparison.regressions)
    end
    if ~isempty(comparison.improvements)
        disp("Improvements:")
        disp(comparison.improvements)
    end
end
end

function analysis = coerceAnalysis(inputStruct, projectRoots, maxFunctions, maxLines, ...
        shouldIncludeToolboxInPriorityTargets, fallbackLabel)
if isfield(inputStruct, "functionMetrics") && isfield(inputStruct, "lineMetrics")
    analysis = inputStruct;
    if ~isfield(analysis, "label") || strlength(strtrim(string(analysis.label))) == 0
        analysis.label = fallbackLabel;
    end
    return
end

analysis = analyzeProfileInfo(inputStruct, ...
    projectRoots=projectRoots, ...
    maxFunctions=maxFunctions, ...
    maxLines=maxLines, ...
    shouldIncludeToolboxInPriorityTargets=shouldIncludeToolboxInPriorityTargets, ...
    label=fallbackLabel);
end

function functionDiffs = compareFunctionTables(baselineTable, candidateTable)
if isempty(baselineTable) && isempty(candidateTable)
    functionDiffs = emptyFunctionDiffTable();
    return
end

baselineKeys = string(baselineTable.completeName);
candidateKeys = string(candidateTable.completeName);
keys = unique([baselineKeys; candidateKeys], "stable");

[isBaseline, baselineIndex] = ismember(keys, baselineKeys);
[isCandidate, candidateIndex] = ismember(keys, candidateKeys);

baselineTotalTime = selectNumericRows(baselineTable.totalTime, baselineIndex, isBaseline);
candidateTotalTime = selectNumericRows(candidateTable.totalTime, candidateIndex, isCandidate);
baselineSelfTime = selectNumericRows(baselineTable.selfTime, baselineIndex, isBaseline);
candidateSelfTime = selectNumericRows(candidateTable.selfTime, candidateIndex, isCandidate);
baselineNumCalls = selectNumericRows(baselineTable.numCalls, baselineIndex, isBaseline);
candidateNumCalls = selectNumericRows(candidateTable.numCalls, candidateIndex, isCandidate);
baselineTotalTimePerCall = selectNumericRows(baselineTable.totalTimePerCall, baselineIndex, isBaseline);
candidateTotalTimePerCall = selectNumericRows(candidateTable.totalTimePerCall, candidateIndex, isCandidate);
baselineSelfTimePerCall = selectNumericRows(baselineTable.selfTimePerCall, baselineIndex, isBaseline);
candidateSelfTimePerCall = selectNumericRows(candidateTable.selfTimePerCall, candidateIndex, isCandidate);

functionName = chooseTextRows(candidateTable.functionName, candidateIndex, isCandidate, ...
    baselineTable.functionName, baselineIndex, isBaseline);
fileName = chooseTextRows(candidateTable.fileName, candidateIndex, isCandidate, ...
    baselineTable.fileName, baselineIndex, isBaseline);
isProjectCode = chooseLogicalRows(candidateTable.isProjectCode, candidateIndex, isCandidate, ...
    baselineTable.isProjectCode, baselineIndex, isBaseline);
changeType = buildChangeType(isBaseline, isCandidate);

functionDiffs = table(keys, functionName, fileName, isProjectCode, changeType, ...
    baselineNumCalls, candidateNumCalls, ...
    deltaValue(candidateNumCalls, baselineNumCalls), percentDelta(candidateNumCalls, baselineNumCalls), ...
    baselineTotalTime, candidateTotalTime, ...
    deltaValue(candidateTotalTime, baselineTotalTime), percentDelta(candidateTotalTime, baselineTotalTime), ...
    baselineSelfTime, candidateSelfTime, ...
    deltaValue(candidateSelfTime, baselineSelfTime), percentDelta(candidateSelfTime, baselineSelfTime), ...
    baselineTotalTimePerCall, candidateTotalTimePerCall, ...
    deltaValue(candidateTotalTimePerCall, baselineTotalTimePerCall), ...
    percentDelta(candidateTotalTimePerCall, baselineTotalTimePerCall), ...
    baselineSelfTimePerCall, candidateSelfTimePerCall, ...
    deltaValue(candidateSelfTimePerCall, baselineSelfTimePerCall), ...
    percentDelta(candidateSelfTimePerCall, baselineSelfTimePerCall), ...
    'VariableNames', ["completeName", "functionName", "fileName", "isProjectCode", ...
    "changeType", "baselineNumCalls", "candidateNumCalls", "deltaNumCalls", ...
    "percentDeltaNumCalls", "baselineTotalTime", "candidateTotalTime", ...
    "deltaTotalTime", "percentDeltaTotalTime", "baselineSelfTime", ...
    "candidateSelfTime", "deltaSelfTime", "percentDeltaSelfTime", ...
    "baselineTotalTimePerCall", "candidateTotalTimePerCall", ...
    "deltaTotalTimePerCall", "percentDeltaTotalTimePerCall", ...
    "baselineSelfTimePerCall", "candidateSelfTimePerCall", ...
    "deltaSelfTimePerCall", "percentDeltaSelfTimePerCall"]);
end

function lineDiffs = compareLineTables(baselineTable, candidateTable)
if isempty(baselineTable) && isempty(candidateTable)
    lineDiffs = emptyLineDiffTable();
    return
end

baselineKeys = composeLineKeys(baselineTable);
candidateKeys = composeLineKeys(candidateTable);
keys = unique([baselineKeys; candidateKeys], "stable");

[isBaseline, baselineIndex] = ismember(keys, baselineKeys);
[isCandidate, candidateIndex] = ismember(keys, candidateKeys);

baselineLineTime = selectNumericRows(baselineTable.lineTime, baselineIndex, isBaseline);
candidateLineTime = selectNumericRows(candidateTable.lineTime, candidateIndex, isCandidate);
baselineLineFraction = selectNumericRows(baselineTable.lineFractionOfFunction, baselineIndex, isBaseline);
candidateLineFraction = selectNumericRows(candidateTable.lineFractionOfFunction, candidateIndex, isCandidate);

fileName = chooseTextRows(candidateTable.fileName, candidateIndex, isCandidate, ...
    baselineTable.fileName, baselineIndex, isBaseline);
functionName = chooseTextRows(candidateTable.functionName, candidateIndex, isCandidate, ...
    baselineTable.functionName, baselineIndex, isBaseline);
lineNumber = chooseNumericRows(candidateTable.lineNumber, candidateIndex, isCandidate, ...
    baselineTable.lineNumber, baselineIndex, isBaseline);
isProjectCode = chooseLogicalRows(candidateTable.isProjectCode, candidateIndex, isCandidate, ...
    baselineTable.isProjectCode, baselineIndex, isBaseline);
changeType = buildChangeType(isBaseline, isCandidate);

lineDiffs = table(fileName, functionName, lineNumber, isProjectCode, changeType, ...
    baselineLineTime, candidateLineTime, ...
    deltaValue(candidateLineTime, baselineLineTime), percentDelta(candidateLineTime, baselineLineTime), ...
    baselineLineFraction, candidateLineFraction, ...
    deltaValue(candidateLineFraction, baselineLineFraction), ...
    percentDelta(candidateLineFraction, baselineLineFraction), ...
    'VariableNames', ["fileName", "functionName", "lineNumber", "isProjectCode", ...
    "changeType", "baselineLineTime", "candidateLineTime", "deltaLineTime", ...
    "percentDeltaLineTime", "baselineLineFractionOfFunction", ...
    "candidateLineFractionOfFunction", "deltaLineFractionOfFunction", ...
    "percentDeltaLineFractionOfFunction"]);
end

function targets = buildComparisonTargets(functionDiffs, lineDiffs, isRegression, ...
        maxFunctions, maxLines, shouldIncludeToolbox)
targets = emptyComparisonTargetsTable();

functionRows = functionDiffs;
lineRows = lineDiffs;
if ~shouldIncludeToolbox
    functionRows = functionRows(functionRows.isProjectCode, :);
    lineRows = lineRows(lineRows.isProjectCode, :);
end

direction = ternary(isRegression, 1, -1);
functionRows = functionRows(direction*functionRows.deltaSelfTime > 0, :);
lineRows = lineRows(direction*lineRows.deltaLineTime > 0, :);
functionRows = sortComparisonRows(functionRows, "deltaSelfTime", isRegression);
lineRows = sortComparisonRows(lineRows, "deltaLineTime", isRegression);

targetTables = {};
numTables = 0;

if ~isempty(functionRows)
    rowCount = min(maxFunctions, height(functionRows));
    selectedRows = functionRows(1:rowCount, :);
    numTables = numTables + 1;
    targetTables{numTables} = table( ...
        repmat("function", rowCount, 1), ...
        selectedRows.completeName, ...
        selectedRows.completeName, ...
        selectedRows.fileName, ...
        selectedRows.functionName, ...
        NaN(rowCount, 1), ...
        selectedRows.changeType, ...
        repmat("selfTime", rowCount, 1), ...
        selectedRows.baselineSelfTime, ...
        selectedRows.candidateSelfTime, ...
        selectedRows.deltaSelfTime, ...
        selectedRows.percentDeltaSelfTime, ...
        buildComparisonReasons(selectedRows.changeType, "function", isRegression), ...
        'VariableNames', comparisonTargetVariableNames());
end

if ~isempty(lineRows)
    rowCount = min(maxLines, height(lineRows));
    selectedRows = lineRows(1:rowCount, :);
    numTables = numTables + 1;
    targetTables{numTables} = table( ...
        repmat("line", rowCount, 1), ...
        selectedRows.fileName + ":" + string(selectedRows.lineNumber), ...
        repmat("", rowCount, 1), ...
        selectedRows.fileName, ...
        selectedRows.functionName, ...
        selectedRows.lineNumber, ...
        selectedRows.changeType, ...
        repmat("lineTime", rowCount, 1), ...
        selectedRows.baselineLineTime, ...
        selectedRows.candidateLineTime, ...
        selectedRows.deltaLineTime, ...
        selectedRows.percentDeltaLineTime, ...
        buildComparisonReasons(selectedRows.changeType, "line", isRegression), ...
        'VariableNames', comparisonTargetVariableNames());
end

if numTables == 0
    return
end

targets = vertcat(targetTables{1:numTables});
targets = sortComparisonRows(targets, "deltaValue", isRegression);
end

function summary = buildComparisonSummary(comparison)
summary = strings(0, 1);
summary(end + 1) = comparison.baselineLabel + " -> " + comparison.candidateLabel + ".";
summary(end + 1) = compose("Compared %d functions and %d executed lines.", ...
    height(comparison.functionDiffs), height(comparison.lineDiffs));

if ~isempty(comparison.regressions)
    summary(end + 1) = formatComparisonSummaryItem(comparison.regressions(1, :));
end
if ~isempty(comparison.improvements)
    summary(end + 1) = formatComparisonSummaryItem(comparison.improvements(1, :));
end
end

function text = formatComparisonSummaryItem(targetRow)
switch targetRow.changeType
    case "new"
        text = targetRow.targetName + " is newly executed in the candidate profile.";
    case "removed"
        text = targetRow.targetName + " is no longer executed in the candidate profile.";
    otherwise
        directionText = ternary(targetRow.deltaValue >= 0, "regressed", "improved");
        text = targetRow.targetName + " " + directionText + " in " + targetRow.metricName + ...
            " by " + formatSeconds(abs(targetRow.deltaValue)) + ".";
end
end

function keys = composeLineKeys(lineTable)
if isempty(lineTable)
    keys = strings(0, 1);
    return
end

keys = lineTable.fileName + "|" + lineTable.functionName + "|" + string(lineTable.lineNumber);
end

function text = chooseLabel(label, fallback)
text = string(label);
if strlength(strtrim(text)) == 0
    text = string(fallback);
end
end

function values = selectNumericRows(sourceValues, sourceIndex, isPresent)
values = zeros(numel(sourceIndex), 1);
values(isPresent) = sourceValues(sourceIndex(isPresent));
end

function values = chooseNumericRows(candidateValues, candidateIndex, isCandidate, ...
        baselineValues, baselineIndex, isBaseline)
values = NaN(numel(candidateIndex), 1);
values(isCandidate) = candidateValues(candidateIndex(isCandidate));
onlyBaseline = ~isCandidate & isBaseline;
values(onlyBaseline) = baselineValues(baselineIndex(onlyBaseline));
end

function values = chooseLogicalRows(candidateValues, candidateIndex, isCandidate, ...
        baselineValues, baselineIndex, isBaseline)
values = false(numel(candidateIndex), 1);
values(isCandidate) = candidateValues(candidateIndex(isCandidate));
onlyBaseline = ~isCandidate & isBaseline;
values(onlyBaseline) = baselineValues(baselineIndex(onlyBaseline));
end

function values = chooseTextRows(candidateValues, candidateIndex, isCandidate, ...
        baselineValues, baselineIndex, isBaseline)
values = strings(numel(candidateIndex), 1);
values(isCandidate) = candidateValues(candidateIndex(isCandidate));
onlyBaseline = ~isCandidate & isBaseline;
values(onlyBaseline) = baselineValues(baselineIndex(onlyBaseline));
end

function types = buildChangeType(isBaseline, isCandidate)
types = repmat("changed", numel(isBaseline), 1);
types(isCandidate & ~isBaseline) = "new";
types(isBaseline & ~isCandidate) = "removed";
types(isBaseline & isCandidate) = "changed";
end

function delta = deltaValue(candidateValues, baselineValues)
delta = candidateValues - baselineValues;
isMissing = isnan(candidateValues) | isnan(baselineValues);
delta(isMissing) = NaN;
end

function delta = percentDelta(candidateValues, baselineValues)
delta = NaN(size(candidateValues));
isValid = isfinite(candidateValues) & isfinite(baselineValues) & baselineValues ~= 0;
delta(isValid) = 100*(candidateValues(isValid) - baselineValues(isValid)) ./ baselineValues(isValid);
sameZero = candidateValues == 0 & baselineValues == 0;
delta(sameZero) = 0;
end

function reasons = buildComparisonReasons(changeTypes, targetType, isRegression)
reasons = strings(numel(changeTypes), 1);
for iReason = 1:numel(changeTypes)
    changeType = changeTypes(iReason);
    switch changeType
        case "new"
            reasons(iReason) = "This " + targetType + " is newly executed in the candidate profile.";
        case "removed"
            reasons(iReason) = "This " + targetType + " is no longer executed in the candidate profile.";
        otherwise
            if targetType == "function"
                reasons(iReason) = ternary(isRegression, ...
                    "Self time increased in this function.", ...
                    "Self time decreased in this function.");
            else
                reasons(iReason) = ternary(isRegression, ...
                    "This executed line got slower.", ...
                    "This executed line got faster.");
            end
    end
end
end

function priority = changeTypePriority(changeTypes)
priority = zeros(numel(changeTypes), 1);
priority(changeTypes == "changed") = 0;
priority(changeTypes == "new") = 1;
priority(changeTypes == "removed") = 2;
end

function rows = sortComparisonRows(rows, metricName, isRegression)
if isempty(rows)
    return
end

rows.changePriority = changeTypePriority(rows.changeType);
metricValues = rows.(metricName);
if isRegression
    rows.metricPriority = -metricValues;
else
    rows.metricPriority = metricValues;
end

rows = sortrows(rows, ["changePriority", "metricPriority"], "ascend", MissingPlacement="last");
rows.changePriority = [];
rows.metricPriority = [];
end

function names = comparisonTargetVariableNames()
names = ["targetType", "targetName", "completeName", "fileName", ...
    "functionName", "lineNumber", "changeType", "metricName", ...
    "baselineValue", "candidateValue", "deltaValue", "percentDelta", "reason"];
end

function value = ternary(condition, whenTrue, whenFalse)
if condition
    value = whenTrue;
else
    value = whenFalse;
end
end

function text = formatSeconds(value)
text = string(compose("%.3g s", value));
end

function functionDiffs = emptyFunctionDiffTable()
functionDiffs = table('Size', [0 25], ...
    'VariableTypes', ["string", "string", "string", "logical", "string", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double"], ...
    'VariableNames', ["completeName", "functionName", "fileName", ...
    "isProjectCode", "changeType", "baselineNumCalls", "candidateNumCalls", ...
    "deltaNumCalls", "percentDeltaNumCalls", "baselineTotalTime", ...
    "candidateTotalTime", "deltaTotalTime", "percentDeltaTotalTime", ...
    "baselineSelfTime", "candidateSelfTime", "deltaSelfTime", ...
    "percentDeltaSelfTime", "baselineTotalTimePerCall", ...
    "candidateTotalTimePerCall", "deltaTotalTimePerCall", ...
    "percentDeltaTotalTimePerCall", "baselineSelfTimePerCall", ...
    "candidateSelfTimePerCall", "deltaSelfTimePerCall", ...
    "percentDeltaSelfTimePerCall"]);
end

function lineDiffs = emptyLineDiffTable()
lineDiffs = table('Size', [0 13], ...
    'VariableTypes', ["string", "string", "double", "logical", "string", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double"], ...
    'VariableNames', ["fileName", "functionName", "lineNumber", ...
    "isProjectCode", "changeType", "baselineLineTime", "candidateLineTime", ...
    "deltaLineTime", "percentDeltaLineTime", ...
    "baselineLineFractionOfFunction", "candidateLineFractionOfFunction", ...
    "deltaLineFractionOfFunction", "percentDeltaLineFractionOfFunction"]);
end

function targets = emptyComparisonTargetsTable()
targets = table('Size', [0 13], ...
    'VariableTypes', ["string", "string", "string", "string", "string", ...
    "double", "string", "string", "double", "double", "double", ...
    "double", "string"], ...
    'VariableNames', comparisonTargetVariableNames());
end
