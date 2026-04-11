function analysis = analyzeProfileInfo(profileInfo, options)
% Analyze MATLAB profiler output and surface hotspot tables.
%
% Use `analyzeProfileInfo` to convert `profile("info")` output into ranked
% function and line tables plus action-oriented hotspot guidance.
%
% ```matlab
% profileInfo = profile("info");
% analysis = analyzeProfileInfo(profileInfo, projectRoots=pwd);
% ```

arguments
    profileInfo (1,1) struct
    options.projectRoots = string(pwd)
    options.maxFunctions (1,1) double {mustBeInteger, mustBePositive} = 25
    options.maxLines (1,1) double {mustBeInteger, mustBePositive} = 25
    options.label = ""
    options.elapsedTime (1,1) double = NaN
    options.shouldIncludeToolboxInPriorityTargets (1,1) logical = false
end

if isfinite(options.elapsedTime) && options.elapsedTime < 0
    error("analyzeProfileInfo:InvalidElapsedTime", ...
        "elapsedTime must be nonnegative when it is finite.");
end

projectRoots = normalizeProjectRoots(options.projectRoots);
[functionMetrics, functionTable] = buildFunctionMetrics(profileInfo, projectRoots);
lineMetrics = buildLineMetrics(functionTable, functionMetrics);

projectFunctionMask = functionMetrics.isProjectCode;
externalFunctionMask = functionMetrics.isToolboxCode | functionMetrics.isBuiltIn;
projectActionableLineMask = lineMetrics.isProjectCode & lineMetrics.isActionableLine;
projectCallsiteLineMask = lineMetrics.isProjectCode & lineMetrics.isCallBoundary;

analysis = struct( ...
    "label", string(options.label), ...
    "elapsedTime", options.elapsedTime, ...
    "profileInfo", profileInfo, ...
    "functionMetrics", functionMetrics, ...
    "lineMetrics", lineMetrics, ...
    "topByTotalTime", topFunctionRows(functionMetrics, "totalTime", options.maxFunctions), ...
    "topBySelfTime", topFunctionRows(functionMetrics, "selfTime", options.maxFunctions), ...
    "topByNumCalls", topFunctionRows(functionMetrics, "numCalls", options.maxFunctions), ...
    "topByTotalTimePerCall", topFunctionRows(functionMetrics, "totalTimePerCall", options.maxFunctions), ...
    "topBySelfTimePerCall", topFunctionRows(functionMetrics, "selfTimePerCall", options.maxFunctions), ...
    "topLines", topLineRows(lineMetrics, "lineTime", options.maxLines), ...
    "topProjectByTotalTime", topFunctionRows(functionMetrics(projectFunctionMask, :), "totalTime", options.maxFunctions), ...
    "topProjectBySelfTime", topFunctionRows(functionMetrics(projectFunctionMask, :), "selfTime", options.maxFunctions), ...
    "topProjectByNumCalls", topFunctionRows(functionMetrics(projectFunctionMask, :), "numCalls", options.maxFunctions), ...
    "topProjectByTotalTimePerCall", topFunctionRows(functionMetrics(projectFunctionMask, :), "totalTimePerCall", options.maxFunctions), ...
    "topProjectBySelfTimePerCall", topFunctionRows(functionMetrics(projectFunctionMask, :), "selfTimePerCall", options.maxFunctions), ...
    "topExternalBySelfTime", topFunctionRows(functionMetrics(externalFunctionMask, :), "selfTime", options.maxFunctions), ...
    "topActionableLines", topLineRows(lineMetrics(projectActionableLineMask, :), "lineTime", options.maxLines), ...
    "topCallsiteLines", topLineRows(lineMetrics(projectCallsiteLineMask, :), "lineTime", options.maxLines), ...
    "priorityTargets", buildPriorityTargets(functionMetrics, lineMetrics, ...
        options.maxFunctions, options.maxLines, options.shouldIncludeToolboxInPriorityTargets), ...
    "summary", strings(0, 1));

analysis.summary = buildAnalysisSummary(analysis, options.shouldIncludeToolboxInPriorityTargets);
end

function [functionMetrics, functionTable] = buildFunctionMetrics(profileInfo, projectRoots)
functionMetrics = emptyFunctionMetricsTable();
functionTable = struct([]);

if ~isfield(profileInfo, "FunctionTable") || isempty(profileInfo.FunctionTable)
    return
end

functionTable = profileInfo.FunctionTable;
numFunctions = numel(functionTable);

completeName = strings(numFunctions, 1);
functionName = strings(numFunctions, 1);
fileName = strings(numFunctions, 1);
type = strings(numFunctions, 1);
numCalls = zeros(numFunctions, 1);
totalTime = zeros(numFunctions, 1);
childTime = NaN(numFunctions, 1);
selfTime = NaN(numFunctions, 1);
totalTimePerCall = NaN(numFunctions, 1);
selfTimePerCall = NaN(numFunctions, 1);
selfFraction = NaN(numFunctions, 1);
isLeaf = false(numFunctions, 1);
isRecursive = false(numFunctions, 1);
numChildren = zeros(numFunctions, 1);
numParents = zeros(numFunctions, 1);
numExecutedLines = zeros(numFunctions, 1);
isBuiltIn = false(numFunctions, 1);
isToolboxCode = false(numFunctions, 1);
isProjectCode = false(numFunctions, 1);

for iFunction = 1:numFunctions
    entry = functionTable(iFunction);

    completeName(iFunction) = getFunctionText(entry, "CompleteName", ...
        getFunctionText(entry, "FunctionName", ""));
    functionName(iFunction) = getFunctionText(entry, "FunctionName", completeName(iFunction));
    fileName(iFunction) = getFunctionText(entry, "FileName", "");
    type(iFunction) = getFunctionText(entry, "Type", "");
    numCalls(iFunction) = getFunctionDouble(entry, "NumCalls", 0);
    totalTime(iFunction) = getFunctionDouble(entry, "TotalTime", 0);

    children = getFunctionField(entry, "Children", []);
    parents = getFunctionField(entry, "Parents", []);
    executedLines = getFunctionField(entry, "ExecutedLines", []);

    isLeaf(iFunction) = isempty(children);
    numChildren(iFunction) = countConnectionEntries(children);
    numParents(iFunction) = countConnectionEntries(parents);
    numExecutedLines(iFunction) = countExecutedLines(executedLines);

    childTime(iFunction) = computeChildTime(children);
    if ~isnan(childTime(iFunction))
        selfTime(iFunction) = max(totalTime(iFunction) - childTime(iFunction), 0);
    end

    totalTimePerCall(iFunction) = safeDivide(totalTime(iFunction), numCalls(iFunction));
    selfTimePerCall(iFunction) = safeDivide(selfTime(iFunction), numCalls(iFunction));
    selfFraction(iFunction) = safeDivide(selfTime(iFunction), totalTime(iFunction));
    isRecursive(iFunction) = logical(getFunctionDouble(entry, "IsRecursive", 0));

    [isProjectCode(iFunction), isToolboxCode(iFunction), isBuiltIn(iFunction)] = ...
        classifyFunctionLocation(fileName(iFunction), type(iFunction), projectRoots);
end

functionMetrics = table(completeName, functionName, fileName, type, ...
    numCalls, totalTime, childTime, selfTime, totalTimePerCall, ...
    selfTimePerCall, selfFraction, isLeaf, isRecursive, numChildren, ...
    numParents, numExecutedLines, isBuiltIn, isToolboxCode, isProjectCode);
end

function lineMetrics = buildLineMetrics(functionTable, functionMetrics)
lineMetrics = emptyLineMetricsTable();

if isempty(functionTable)
    return
end

lineTables = cell(height(functionMetrics), 1);
numTables = 0;
for iFunction = 1:numel(functionTable)
    entry = functionTable(iFunction);
    lineTable = parseExecutedLines(entry, functionMetrics(iFunction, :));
    if isempty(lineTable)
        continue
    end
    numTables = numTables + 1;
    lineTables{numTables} = lineTable;
end

if numTables == 0
    return
end

lineMetrics = vertcat(lineTables{1:numTables});
end

function lineTable = parseExecutedLines(entry, functionRow)
lineTable = emptyLineMetricsTable();

executedLines = getFunctionField(entry, "ExecutedLines", []);
[lineNumber, lineTime] = extractExecutedLineColumns(executedLines);
if isempty(lineNumber)
    return
end

numLines = numel(lineNumber);
fileName = repmat(functionRow.fileName, numLines, 1);
functionName = repmat(functionRow.functionName, numLines, 1);
lineFractionOfFunction = NaN(numLines, 1);
if functionRow.totalTime > 0
    lineFractionOfFunction = lineTime ./ functionRow.totalTime;
end

functionTotalTime = repmat(functionRow.totalTime, numLines, 1);
functionSelfTime = repmat(functionRow.selfTime, numLines, 1);
functionSelfFraction = repmat(functionRow.selfFraction, numLines, 1);
lineText = strings(numLines, 1);
lineKind = repmat("unknown", numLines, 1);
isCallBoundary = false(numLines, 1);
isActionableLine = false(numLines, 1);
calledTargetName = strings(numLines, 1);
isProjectCode = repmat(functionRow.isProjectCode, numLines, 1);

sourceLines = readSourceLines(functionRow.fileName);
for iLine = 1:numLines
    sourceLineText = sourceLineAt(sourceLines, lineNumber(iLine));
    lineText(iLine) = logicalStatementAt(sourceLines, lineNumber(iLine));
    if strlength(lineText(iLine)) == 0
        lineText(iLine) = sourceLineText;
    end
    [lineKind(iLine), calledTargetName(iLine)] = classifyExecutedLine( ...
        sourceLineText, lineText(iLine), functionSelfFraction(iLine), lineFractionOfFunction(iLine));
    isCallBoundary(iLine) = lineKind(iLine) == "callBoundary";
    isActionableLine(iLine) = lineKind(iLine) == "compute";
end

lineTable = table(fileName, functionName, lineNumber, lineTime, ...
    lineFractionOfFunction, functionTotalTime, functionSelfTime, ...
    functionSelfFraction, lineText, lineKind, isCallBoundary, ...
    isActionableLine, calledTargetName, isProjectCode);
end

function rankedFunctions = topFunctionRows(functionMetrics, metricName, maxFunctions)
rankedFunctions = emptyFunctionMetricsTable();

if isempty(functionMetrics)
    return
end

rankedFunctions = sortrows(functionMetrics, metricName, "descend", MissingPlacement="last");
rankedFunctions = rankedFunctions(1:min(maxFunctions, height(rankedFunctions)), :);
end

function rankedLines = topLineRows(lineMetrics, metricName, maxLines)
rankedLines = emptyLineMetricsTable();

if isempty(lineMetrics)
    return
end

rankedLines = sortrows(lineMetrics, metricName, "descend", MissingPlacement="last");
rankedLines = rankedLines(1:min(maxLines, height(rankedLines)), :);
end

function priorityTargets = buildPriorityTargets(functionMetrics, lineMetrics, ...
        maxFunctions, maxLines, shouldIncludeToolbox)
projectFunctionRows = functionMetrics(functionMetrics.isProjectCode, :);
projectLineRows = lineMetrics(lineMetrics.isProjectCode, :);
projectTargets = buildScopedPriorityTargets(projectFunctionRows, projectLineRows, maxFunctions, maxLines);

priorityTargets = projectTargets;
if shouldIncludeToolbox
    externalFunctionRows = functionMetrics(functionMetrics.isToolboxCode, :);
    externalLineRows = lineMetrics(~lineMetrics.isProjectCode, :);
    externalTargets = buildScopedPriorityTargets(externalFunctionRows, externalLineRows, maxFunctions, maxLines);
    if ~isempty(externalTargets)
        priorityTargets = [priorityTargets; externalTargets];
    end
end

if isempty(priorityTargets)
    return
end

priorityTargets = unique(priorityTargets, "rows", "stable");
priorityTargets = priorityTargets(1:min(max(maxFunctions, maxLines), height(priorityTargets)), :);
end

function targets = buildScopedPriorityTargets(functionRows, lineRows, maxFunctions, maxLines)
targets = emptyPriorityTargetsTable();
numFunctionTargets = min(3, maxFunctions);
numLineTargets = min(3, maxLines);

targetTables = {};
numTables = 0;

selfRows = sortrows(functionRows, "selfTime", "descend", MissingPlacement="last");
selfRows = selfRows(isfinite(selfRows.selfTime) & selfRows.selfTime > 0, :);
if ~isempty(selfRows)
    numTables = numTables + 1;
    targetTables{numTables} = buildFunctionTargets( ...
        selfRows(1:min(numFunctionTargets, height(selfRows)), :), ...
        "selfFunction", ...
        "High self time; optimize the work inside this function first.", ...
        "selfTime");
end

computeLineRows = sortrows(lineRows(lineRows.isActionableLine, :), "lineTime", "descend", MissingPlacement="last");
computeLineRows = computeLineRows(isfinite(computeLineRows.lineTime) & computeLineRows.lineTime > 0, :);
if ~isempty(computeLineRows)
    numTables = numTables + 1;
    targetTables{numTables} = buildLineTargets( ...
        computeLineRows(1:min(numLineTargets, height(computeLineRows)), :), ...
        "computeLine", ...
        "High line time on a compute line; optimize the work performed here.", ...
        "lineTime");
end

callRows = sortrows(functionRows, ["numCalls", "totalTime"], ["descend", "descend"]);
callRows = callRows(callRows.numCalls > 0 & callRows.totalTime > 0, :);
if ~isempty(callRows)
    numTables = numTables + 1;
    targetTables{numTables} = buildFunctionTargets( ...
        callRows(1:min(numFunctionTargets, height(callRows)), :), ...
        "callCountFunction", ...
        "High call count; reduce repeated calls or batch work if possible.", ...
        "numCalls");
end

fanOutRows = functionRows;
fanOutRows.childFraction = safeDivide(fanOutRows.childTime, fanOutRows.totalTime);
fanOutRows = sortrows(fanOutRows, "childTime", "descend", MissingPlacement="last");
fanOutRows = fanOutRows(isfinite(fanOutRows.childFraction) & fanOutRows.childFraction >= 0.5 & ...
    fanOutRows.childTime > 0, :);
if ~isempty(fanOutRows)
    numTables = numTables + 1;
    targetTables{numTables} = buildFunctionTargets( ...
        fanOutRows(1:min(numFunctionTargets, height(fanOutRows)), :), ...
        "fanOutFunction", ...
        "Most time is spent in child calls; inspect callees or reduce fan-out from this parent.", ...
        "childTime");
end

callsiteLineRows = sortrows(lineRows(lineRows.isCallBoundary, :), "lineTime", "descend", MissingPlacement="last");
callsiteLineRows = callsiteLineRows(isfinite(callsiteLineRows.lineTime) & callsiteLineRows.lineTime > 0, :);
if ~isempty(callsiteLineRows)
    numTables = numTables + 1;
    targetTables{numTables} = buildLineTargets( ...
        callsiteLineRows(1:min(numLineTargets, height(callsiteLineRows)), :), ...
        "callsiteLine", ...
        "", ...
        "lineTime");
end

if numTables == 0
    return
end

targets = vertcat(targetTables{1:numTables});
end

function summary = buildAnalysisSummary(analysis, shouldIncludeToolbox)
maxSummaryItems = min(4, height(analysis.priorityTargets));
summary = strings(3 + maxSummaryItems + 2, 1);
summaryCount = 0;

if strlength(strtrim(analysis.label)) > 0
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "Profile label: " + analysis.label;
end
if isfinite(analysis.elapsedTime)
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "Wall time: " + formatSeconds(analysis.elapsedTime) + ".";
end
summaryCount = summaryCount + 1;
summary(summaryCount) = compose("Profiled %d functions and %d executed lines.", ...
    height(analysis.functionMetrics), height(analysis.lineMetrics));

if ~isempty(analysis.priorityTargets)
    for iTarget = 1:maxSummaryItems
        summaryCount = summaryCount + 1;
        summary(summaryCount) = formatPriorityTarget(analysis.priorityTargets(iTarget, :));
    end
elseif isempty(analysis.functionMetrics)
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "Profiler output did not contain any function entries.";
elseif shouldIncludeToolbox
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "No actionable project or toolbox hotspots were identified.";
else
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "No project-code hotspots were identified in the profiler output.";
    summaryCount = summaryCount + 1;
    summary(summaryCount) = "Set shouldIncludeToolboxInPriorityTargets=true to surface toolbox or external-code hotspots.";
end

summary = summary(1:summaryCount);
end

function targetTable = buildFunctionTargets(functionRows, targetSubtype, reasonText, metricName)
targetType = repmat("function", height(functionRows), 1);
targetSubtype = repmat(string(targetSubtype), height(functionRows), 1);
targetName = functionRows.completeName;
reason = repmat(string(reasonText), height(functionRows), 1);
metricName = repmat(string(metricName), height(functionRows), 1);
metricValue = functionRows.(char(metricName(1)));
targetTable = table(targetType, targetSubtype, targetName, reason, metricName, metricValue);
end

function targetTable = buildLineTargets(lineRows, targetSubtype, reasonText, metricName)
targetType = repmat("line", height(lineRows), 1);
targetSubtype = repmat(string(targetSubtype), height(lineRows), 1);
targetName = lineRows.fileName + ":" + string(lineRows.lineNumber);
metricName = repmat(string(metricName), height(lineRows), 1);
metricValue = lineRows.(char(metricName(1)));
reason = repmat(string(reasonText), height(lineRows), 1);

if targetSubtype(1) == "callsiteLine"
    for iLine = 1:height(lineRows)
        calleeName = lineRows.calledTargetName(iLine);
        if strlength(calleeName) > 0
            reason(iLine) = "Callee " + calleeName + " dominates this call site; inspect that callee path first.";
        else
            reason(iLine) = "Callee dominates this call site; inspect the callee path first.";
        end
    end
end

targetTable = table(targetType, targetSubtype, targetName, reason, metricName, metricValue);
end

function text = formatPriorityTarget(targetRow)
if targetRow.targetType == "line"
    metricText = formatSeconds(targetRow.metricValue);
else
    metricText = formatMetricValue(targetRow.metricValue);
end

text = targetRow.targetName + " is a priority because " + lower(targetRow.reason) + ...
    " (" + targetRow.metricName + " = " + metricText + ").";
end

function [lineNumber, lineTime] = extractExecutedLineColumns(executedLines)
lineNumber = zeros(0, 1);
lineTime = zeros(0, 1);

if isempty(executedLines)
    return
end

if isnumeric(executedLines) && size(executedLines, 2) >= 3
    lineNumber = double(executedLines(:, 1));
    lineTime = double(executedLines(:, 3));
elseif isstruct(executedLines)
    lineNumber = extractStructField(executedLines, ["LineNumber", "Line", "Number"]);
    lineTime = extractStructField(executedLines, ["Time", "LineTime", "TotalTime"]);
else
    return
end

isValid = isfinite(lineNumber) & isfinite(lineTime);
lineNumber = lineNumber(isValid);
lineTime = lineTime(isValid);
end

function values = extractStructField(structArray, candidateNames)
values = zeros(0, 1);
for iName = 1:numel(candidateNames)
    if isfield(structArray, char(candidateNames(iName)))
        values = double(transpose([structArray.(char(candidateNames(iName)))]));
        return
    end
end
end

function count = countExecutedLines(executedLines)
[lineNumber, ~] = extractExecutedLineColumns(executedLines);
count = numel(lineNumber);
end

function sourceLines = readSourceLines(fileName)
sourceLines = strings(0, 1);
fileName = normalizePathString(fileName);
if strlength(fileName) == 0 || ~isfile(fileName)
    return
end

try
    sourceLines = splitlines(string(fileread(fileName)));
catch
    sourceLines = strings(0, 1);
end
end

function text = sourceLineAt(sourceLines, lineNumber)
text = "";
if isempty(sourceLines)
    return
end

lineIndex = round(lineNumber);
if lineIndex < 1 || lineIndex > numel(sourceLines)
    return
end

text = strtrim(stripLineComment(sourceLines(lineIndex)));
end

function text = logicalStatementAt(sourceLines, lineNumber)
text = "";
if isempty(sourceLines)
    return
end

lineIndex = round(lineNumber);
if lineIndex < 1 || lineIndex > numel(sourceLines)
    return
end

startIndex = lineIndex;
while startIndex > 1 && endsWith(strtrim(stripLineComment(sourceLines(startIndex - 1))), "...")
    startIndex = startIndex - 1;
end

endIndex = lineIndex;
currentLine = strtrim(stripLineComment(sourceLines(endIndex)));
while endsWith(currentLine, "...") && endIndex < numel(sourceLines)
    endIndex = endIndex + 1;
    currentLine = strtrim(stripLineComment(sourceLines(endIndex)));
end

parts = strings(0, 1);
for iPart = startIndex:endIndex
    part = strtrim(stripLineComment(sourceLines(iPart)));
    if strlength(part) == 0
        continue
    end
    part = regexprep(char(part), "\.\.\.\s*$", "");
    parts(end + 1, 1) = strtrim(string(part)); %#ok<AGROW>
end

if isempty(parts)
    return
end

text = strjoin(parts, " ");
end

function text = stripLineComment(text)
text = string(text);
commentIndex = strfind(text, "%");
if isempty(commentIndex)
    return
end

text = extractBefore(text, commentIndex(1));
end

function [lineKind, calledTargetName] = classifyExecutedLine(sourceLineText, statementText, ...
        functionSelfFraction, lineFractionOfFunction)
lineKind = "unknown";
calledTargetName = "";
sourceLineText = strtrim(string(sourceLineText));
statementText = strtrim(string(statementText));

if strlength(statementText) == 0
    statementText = sourceLineText;
end

if strlength(statementText) == 0
    return
end

calledTargetName = extractCallBoundaryTarget(sourceLineText, statementText);
isCallLike = isCallLikeStatement(sourceLineText) || isCallLikeStatement(statementText);
hasExplicitComputeOperator = hasExplicitComputeOperatorText(sourceLineText) || ...
    hasExplicitComputeOperatorText(statementText);
isCallBoundary = strlength(calledTargetName) > 0 && ...
    isCallLike && ~hasExplicitComputeOperator && ...
    isfinite(functionSelfFraction) && functionSelfFraction <= 0.40 && ...
    isfinite(lineFractionOfFunction) && lineFractionOfFunction >= 0.03;

if isCallBoundary
    lineKind = "callBoundary";
    return
end

lineKind = "compute";
calledTargetName = "";
end

function calledTargetName = extractCallBoundaryTarget(sourceLineText, statementText)
candidateTarget = lastCallTarget(sourceLineText);
if strlength(candidateTarget) > 0
    calledTargetName = candidateTarget;
    return
end

calledTargetName = lastCallTarget(statementText);
end

function targetName = lastCallTarget(text)
targetName = "";
text = regexprep(char(strtrim(string(text))), "\.\.\.\s*$", "");
if isempty(text)
    return
end

tokens = regexp(text, "([A-Za-z]\w*(?:\.[A-Za-z]\w*)*)\s*\(", "tokens");
if isempty(tokens)
    return
end

targetName = string(tokens{end}{1});
end

function tf = isCallLikeStatement(text)
text = strtrim(string(text));
if strlength(text) == 0
    tf = false;
    return
end

pattern = "^\s*(?:(?:\[[^\]]*\]|[A-Za-z]\w*(?:\.[A-Za-z]\w*)*)\s*=\s*)?[+\-~]?\s*[A-Za-z]\w*(?:\.[A-Za-z]\w*)*\s*\(.*\)\s*;?\s*$";
tf = ~isempty(regexp(char(text), pattern, "once"));
end

function tf = hasExplicitComputeOperatorText(text)
text = regexprep(char(string(text)), "'[^']*'", "");
text = regexprep(text, '"[^"]*"', "");
tf = contains(text, "\") || contains(text, "*") || contains(text, "/") || contains(text, "^");
end

function [isProjectCode, isToolboxCode, isBuiltIn] = classifyFunctionLocation(fileName, type, projectRoots)
fileName = normalizePathString(fileName);
type = string(type);

isBuiltIn = strlength(fileName) == 0 || contains(lower(type), "builtin");
isProjectCode = false;
if ~isBuiltIn
    for iRoot = 1:numel(projectRoots)
        root = normalizePathString(projectRoots(iRoot));
        if fileName == root || startsWith(fileName, root + "/")
            isProjectCode = true;
            break
        end
    end
end
isToolboxCode = ~isBuiltIn && ~isProjectCode && strlength(fileName) > 0;
end

function roots = normalizeProjectRoots(projectRoots)
roots = string(projectRoots);
roots = roots(:);
roots = roots(strlength(strtrim(roots)) > 0);
if isempty(roots)
    roots = string(pwd);
end

normalizedRoots = strings(size(roots));
for iRoot = 1:numel(roots)
    normalizedRoots(iRoot) = normalizePathString(roots(iRoot));
end

roots = unique(normalizedRoots, "stable");
end

function normalizedPath = normalizePathString(pathText)
normalizedPath = string(pathText);
normalizedPath = strtrim(normalizedPath);
normalizedPath = replace(normalizedPath, "\", "/");
while endsWith(normalizedPath, "/")
    normalizedPath = extractBefore(normalizedPath, strlength(normalizedPath));
end
end

function value = getFunctionText(entry, fieldName, defaultValue)
value = string(defaultValue);
if isfield(entry, fieldName)
    value = string(entry.(fieldName));
end
end

function value = getFunctionDouble(entry, fieldName, defaultValue)
value = defaultValue;
if isfield(entry, fieldName)
    value = double(entry.(fieldName));
end
end

function value = getFunctionField(entry, fieldName, defaultValue)
value = defaultValue;
if isfield(entry, fieldName)
    value = entry.(fieldName);
end
end

function count = countConnectionEntries(entries)
if isempty(entries)
    count = 0;
else
    count = numel(entries);
end
end

function childTime = computeChildTime(children)
if isempty(children)
    childTime = 0;
    return
end

childTime = NaN;
if ~isstruct(children)
    return
end

if isfield(children, "TotalTime")
    childTime = sum([children.TotalTime]);
elseif isfield(children, "Time")
    childTime = sum([children.Time]);
end
end

function value = safeDivide(numerator, denominator)
value = NaN(size(numerator));
isValid = denominator ~= 0 & isfinite(denominator);
value(isValid) = numerator(isValid) ./ denominator(isValid);
end

function text = formatSeconds(value)
text = string(compose("%.3g s", value));
end

function text = formatMetricValue(value)
if isnan(value)
    text = "NaN";
elseif abs(value - round(value)) < 10*eps(max(1, abs(value)))
    text = string(compose("%.0f", value));
else
    text = string(compose("%.3g", value));
end
end

function functionMetrics = emptyFunctionMetricsTable()
functionMetrics = table('Size', [0 19], ...
    'VariableTypes', ["string", "string", "string", "string", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "logical", "logical", "double", "double", "double", ...
    "logical", "logical", "logical"], ...
    'VariableNames', ["completeName", "functionName", "fileName", "type", ...
    "numCalls", "totalTime", "childTime", "selfTime", ...
    "totalTimePerCall", "selfTimePerCall", "selfFraction", "isLeaf", ...
    "isRecursive", "numChildren", "numParents", "numExecutedLines", ...
    "isBuiltIn", "isToolboxCode", "isProjectCode"]);
end

function lineMetrics = emptyLineMetricsTable()
lineMetrics = table('Size', [0 14], ...
    'VariableTypes', ["string", "string", "double", "double", "double", ...
    "double", "double", "double", "string", "string", "logical", ...
    "logical", "string", "logical"], ...
    'VariableNames', ["fileName", "functionName", "lineNumber", ...
    "lineTime", "lineFractionOfFunction", "functionTotalTime", ...
    "functionSelfTime", "functionSelfFraction", "lineText", "lineKind", ...
    "isCallBoundary", "isActionableLine", "calledTargetName", ...
    "isProjectCode"]);
end

function priorityTargets = emptyPriorityTargetsTable()
priorityTargets = table('Size', [0 6], ...
    'VariableTypes', ["string", "string", "string", "string", "string", "double"], ...
    'VariableNames', ["targetType", "targetSubtype", "targetName", ...
    "reason", "metricName", "metricValue"]);
end
